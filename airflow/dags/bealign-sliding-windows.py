import yaml
import datetime
from dateutil.relativedelta import *
from dateutil.rrule import *
from dateutil.parser import *

# The DAG object; we'll need this to instantiate a DAG
from airflow import DAG

# Operators; we need this to operate!
from airflow.utils.task_group import TaskGroup
from airflow.operators.bash import BashOperator
from airflow.models.baseoperator import cross_downstream
from airflow.operators.python import PythonOperator
from airflow.utils.dates import days_ago
from airflow.hooks.base import BaseHook
from airflow.models import Variable

# These args will get passed on to each operator
# You can override them on a per-task basis during operator initialization

from libs.callbacks import task_fail_slack_alert, task_success_slack_alert, dag_fail_slack_alert, dag_success_slack_alert

import os
import sys
import pathlib

p = os.path.abspath(str(pathlib.Path(__file__).parent.absolute()) + '/../../python/')
if p not in sys.path:
    sys.path.append(p)

from export_sequences import export_sequences, export_premsa_sequences, export_bealign_sequences
from export_duplicates import export_duplicates
from export_meta import export_meta
from remove_seq import remove_reference_seq, reserve_only_original_input
from merge_duplicates import merge_duplicates
from fix_duplicates import fix_duplicates
from update_fasta_duplicates import update_fasta_duplicates
from get_raw_duplicates import write_nuc_raw_duplicates

WORKING_DIR = Variable.get("WORKING_DIR")

# ["2020-08", "2020-09", "2020-10"]  # this is used for the baseline analysis forecasting Wave 3
# ["2021-02", "2021-03"]  # this is used for forecasting the next wave

def create_dag(dag_id, schedule, window, default_args):
    with DAG(
        dag_id,
        default_args=default_args,
        description='creates sliding windows based on months',
        schedule_interval=schedule,
        start_date=datetime.datetime(2021, 4, 30),
        # on_failure_callback=dag_fail_slack_alert,
        # on_success_callback=dag_success_slack_alert,
        tags=['selection','sliding'],
        ) as dag:

        OUTPUT_DIR = WORKING_DIR + "/data/sliding-windows-bealign/" + '_'.join(window)
        default_args["params"]["output-dir"] = OUTPUT_DIR
        default_args["params"]["meta-output"] = OUTPUT_DIR + '/master-no-sequences.json'
        default_args["params"]["sequence-output"] = OUTPUT_DIR + '/sequences'

        with open(dag.params["region_cfg"], 'r') as stream:
            regions = yaml.safe_load(stream)

        mk_dir_task = BashOperator(
            task_id='make_directory',
            bash_command='mkdir -p {{params.output}}',
            params={'output': default_args['params']['output-dir']},
            dag=dag,
        )

        export_meta_task = PythonOperator(
                task_id='export_meta',
                python_callable=export_meta,
                op_kwargs={ "config" : default_args['params'] },
                pool='mongo',
                dag=dag,
            )

        export_meta_task.set_upstream(mk_dir_task)
        export_sequences_task = PythonOperator(
                task_id='export_sequences',
                python_callable=export_sequences,
                op_kwargs={ "config" : default_args['params'] },
                pool='mongo',
                dag=dag,
            )

        export_sequences_task.set_upstream(mk_dir_task)

        # For each region
        export_by_gene = []

        for gene in regions.keys():

            reference_filepath = WORKING_DIR + 'reference_genes/reference.' + gene + '_protein.fas'
            filepath_prefix = OUTPUT_DIR + '/sequences.' + gene

            nuc_sequence_output = filepath_prefix + '_nuc.fas'
            uniques_fn = filepath_prefix + '_nuc.uniques.fas'
            duplicate_output = filepath_prefix + '.duplicates.json'
            map_output = filepath_prefix + '.map.json'

            variants_csv_output = filepath_prefix + '.variants.csv'
            variants_json_output = filepath_prefix + '.variants.json'
            filtered_fasta_output = filepath_prefix + '.compressed.filtered.fas'
            filtered_json_output = filepath_prefix + '.filtered.json'
            output_edits_fn = filepath_prefix + '.filtered.edits.json'

            compressed_output_filepath =  filepath_prefix + '.compressed.fas'
            compressor_duplicate_out = filepath_prefix + '.duplicates.variants.json'

            tree_output = filepath_prefix + '.compressed.filtered.fas.rapidnj.bestTree'
            sto_output = filepath_prefix + '.compressed.filtered.sto';

            tmp_output_fn = filepath_prefix + '.tmp.msa'
            output_fn = filepath_prefix + '.msa'

            slac_output_fn = filepath_prefix + '.SLAC.json'
            fel_output_fn = filepath_prefix + '.FEL.json'
            meme_output_fn = filepath_prefix + '.MEME.json'

            summary_output_fn = filepath_prefix + '.json'

            default_args["params"]["nuc-sequence-output"] = nuc_sequence_output
            default_args["params"]["duplicate-output"] = duplicate_output

            with TaskGroup(f"alignment_{gene}") as alignment:

                export_bealign_task = PythonOperator(
                    task_id=f'export_bealign',
                    python_callable=export_bealign_sequences,
                    op_kwargs={ "config" : default_args['params'], 'nuc_output_fn':  nuc_sequence_output, 'gene' : gene },
                    dag=dag,
                )

                # Occasional errors when cleaning up tmp files, so or'ing true
                cleanup_task = BashOperator(
                    task_id=f'cleanup',
                    bash_command="sed -i '/^>/! s/[^ACTG-]/N/g' $NUC_OUTPUT_FN || true",
                    env={'NUC_OUTPUT_FN': nuc_sequence_output, **os.environ},
                    dag=dag
                )

                export_bealign_task >> cleanup_task

            with TaskGroup(f"duplicates_{gene}") as duplicates_group:

                compute_duplicates_task = PythonOperator(
                    task_id=f'write_raw_duplicates',
                    python_callable=write_nuc_raw_duplicates,
                    op_kwargs={ "input" : nuc_sequence_output, "duplicate_output" : duplicate_output, 'uniques_output': uniques_fn },
                    dag=dag,
                )

                compute_duplicates_task

            # $HYPHY LIBPATH=$HYPHYLIBPATH $COMPRESSOR --msa ${FILE}.${GENE}.compressed.fas --regexp "epi_isl_([0-9]+)" --duplicates ${FILE}.${GENE}.duplicates.json --output ${FILE}.${GENE}.variants.csv --json ${FILE}.${GENE}.variants.json --duplicate-out ${FILE}.${GENE}.duplicates.variants.json

            with TaskGroup(f"filter_{gene}") as filter:
                COMPRESSOR = """
                {{ params.hyphy }} LIBPATH={{params.hyphy_lib_path}} {{ params.compressor }} --msa $FASTA_FN --regexp "epi_isl_([0-9]+)" --duplicates $DUPLICATE_FN --output $VARIANTS_CSV_FN  --json $VARIANTS_JSON_FN --duplicate-out $COMPRESSOR_DUPLICATE_OUT
                """
                compressor_task = BashOperator(
                    task_id=f'compressor',
                    bash_command=COMPRESSOR,
                    env={'FASTA_FN': uniques_fn, 'DUPLICATE_FN': duplicate_output, 'VARIANTS_CSV_FN': variants_csv_output, 'VARIANTS_JSON_FN': variants_json_output, 'COMPRESSOR_DUPLICATE_OUT': compressor_duplicate_out, **os.environ},
                    dag=dag
                )

                # --output-edits ${FILE}.${GENE}.filtered.edits.json
                COMPRESSOR2 = """
                {{ params.hyphy }} LIBPATH={{params.hyphy_lib_path}} {{ params.compressor2 }} --msa $FASTA_FN --duplicates $DUPLICATE_FN --csv $VARIANTS_CSV_FN  --byseq $VARIANTS_JSON_FN --p 0.95 --output $FILTERED_FASTA_FN --json $FILTERED_JSON_FN --output-edits ${OUTPUT_EDITS}
                """
                compressor_two_task = BashOperator(
                    task_id=f'compressor_two',
                    bash_command=COMPRESSOR2,
                    env={'FASTA_FN': uniques_fn, 'DUPLICATE_FN': compressor_duplicate_out, 'VARIANTS_CSV_FN': variants_csv_output, 'VARIANTS_JSON_FN': variants_json_output, 'FILTERED_FASTA_FN': filtered_fasta_output, 'FILTERED_JSON_FN': filtered_json_output, 'OUTPUT_EDITS': output_edits_fn, **os.environ},
                    dag=dag
                )

                compressor_task >> compressor_two_task

            INFER_TREE = """
            seqmagick convert $FILTERED_FASTA_FN $STO_OUTPUT;
            rapidnj $STO_OUTPUT -i sth > $TREE_OUTPUT
            sed -i "s/'//g" $TREE_OUTPUT;
            """

            infer_tree_task = BashOperator(
                task_id=f'infer_tree_{gene}',
                bash_command=INFER_TREE,
                env={'FILTERED_FASTA_FN': filtered_fasta_output, 'STO_OUTPUT': sto_output, 'TREE_OUTPUT': tree_output, **os.environ},
                dag=dag
            )

            slac_task = BashOperator(
                task_id=f'slac_{gene}',
                bash_command="{{ params.hyphy }} LIBPATH={{params.hyphy_lib_path}} slac --kill-zero-lengths Constrain ENV='_DO_TREE_REBALANCE_=1' --alignment $FILTERED_FASTA_FN --tree $TREE_OUTPUT --branches All --samples 0 --output $SLAC_OUTPUT",
                env={'FILTERED_FASTA_FN': filtered_fasta_output, 'TREE_OUTPUT': tree_output, 'SLAC_OUTPUT': slac_output_fn, **os.environ},
                dag=dag,
            )

            big_data_flags='--full-model No'

            fel_task = BashOperator(
                task_id=f'fel_{gene}',
                bash_command="{{ params.hyphy }} LIBPATH={{params.hyphy_lib_path}} fel --kill-zero-lengths Constrain ENV='_DO_TREE_REBALANCE_=1' $BIG_DATA_FLAGS --alignment $FILTERED_FASTA_FN --tree $TREE_OUTPUT --branches Internal --output $FEL_OUTPUT",
                env={'FILTERED_FASTA_FN': filtered_fasta_output, 'TREE_OUTPUT': tree_output, 'FEL_OUTPUT': fel_output_fn, 'BIG_DATA_FLAGS': big_data_flags, **os.environ},
                dag=dag,
            )

            meme_task = BashOperator(
                task_id=f'meme_{gene}',
                bash_command="{{ params.hyphy }} LIBPATH={{params.hyphy_lib_path}} meme --kill-zero-lengths Constrain ENV='_DO_TREE_REBALANCE_=1' $BIG_DATA_FLAGS --alignment $FILTERED_FASTA_FN --tree $TREE_OUTPUT --branches Internal --output $MEME_OUTPUT",
                env={'FILTERED_FASTA_FN': filtered_fasta_output, 'TREE_OUTPUT': tree_output, 'MEME_OUTPUT': meme_output_fn, 'BIG_DATA_FLAGS': big_data_flags, **os.environ},
                dag=dag,
            )

            # fubar_task = BashOperator(
            #     task_id='fubar_{gene}',
            #     bash_command='mkdir -p {{params.working_dir}}/data/fasta/{{params.date}}',
            #     dag=dag,
            # )

            # prime_task = BashOperator(
            #     task_id='prime_{gene}',
            #     bash_command='mkdir -p {{params.working_dir}}/data/fasta/{{params.date}}',
            #     dag=dag,
            # )

            annotation_file = filepath_prefix + '.annotation.json'
            copy_annotation_task = BashOperator(
                task_id=f'copy_annotation_{gene}',
                bash_command='cp {{params.working_dir}}/data/comparative-annotation.json {{params.annotation_file}}',
                params={'annotation_file': annotation_file, 'working_dir': WORKING_DIR},
                dag=dag
            )

            summarize_gene_task = BashOperator(
                task_id=f'summarize_gene_{gene}',
                bash_command='{{ params.python }} {{params.working_dir}}/python/summarize_gene.py -T {{params.working_dir}}/data/ctl/epitopes.json -B {{params.working_dir}}/data/single_mut_effects.csv -D $MASTERNOFASTA -d $DUPLICATES -s $SLAC_OUTPUT -f $FEL_OUTPUT -m $MEME_OUTPUT -P 0.1 --output  $SUMMARY_OUTPUT -c $COMPRESSED_OUTPUT_FN -E {{params.working_dir}}/data/evo_annotation.json -A {{params.working_dir}}/data/mafs.csv -V {{params.working_dir}}/data/evo_freqs.csv -F $FRAGMENT --frame_shift $ADDSHIFT --fragment_shift $SHIFT -S $OFFSET -O $ANNOTATION',
                params={'python': default_args['params']['python'], 'working_dir': WORKING_DIR},
                env={
                    'MASTERNOFASTA': default_args["params"]["meta-output"],
                    'DUPLICATES': duplicate_output,
                    'SLAC_OUTPUT': slac_output_fn,
                    'FEL_OUTPUT': fel_output_fn,
                    'MEME_OUTPUT': meme_output_fn,
                    'SUMMARY_OUTPUT': summary_output_fn,
                    'COMPRESSED_OUTPUT_FN': filtered_fasta_output,
                    'FRAGMENT': str(regions[gene]['fragment']),
                    'ADDSHIFT': str(regions[gene]['add_one']),
                    'SHIFT': str(regions[gene]['shift']),
                    'OFFSET': str(regions[gene]['offset']),
                    'ANNOTATION': annotation_file,
                    **os.environ},
                dag=dag,
            )

            summarize_gene_task.set_upstream(export_meta_task)
            alignment.set_upstream(export_sequences_task)
            export_by_gene.append(alignment >> duplicates_group >> filter >> infer_tree_task >> [slac_task, fel_task, meme_task] >> copy_annotation_task >> summarize_gene_task)

        dag.doc_md = __doc__

        # Add export meta and export sequence tasks to be executed in parallel
        cross_downstream([export_meta_task, export_sequences_task], export_by_gene)

        return dag

sliding_windows = [
                    ("2020-08-01", "2020-10-31"),
                    ("2021-02-01", "2021-03-31"),
                    ("2019-12-01", "2020-02-29"),
                    ("2020-01-01", "2020-03-31"),
                    ("2020-02-01", "2020-04-30"),
                    ("2020-03-01", "2020-05-31"),
                    ("2020-04-01", "2020-06-30"),
                    ("2020-05-01", "2020-07-31"),
                    ("2020-06-01", "2020-08-31"),
                    ("2020-07-01", "2020-09-30")
                  ]

# Supplement with 3 month sliding windows since beginning of pandemic
TODAY = datetime.date.today()
LASTMONTH = TODAY-relativedelta(months=+1, day=31)
THREEMONTHSAGO = TODAY-relativedelta(months=+3, day=1)

starts = [dt.strftime('%Y-%m-%d') for dt in rrule(MONTHLY, interval=1,bymonthday=(1),dtstart=parse("20191201T000000"), until=THREEMONTHSAGO)]
ends = [dt.strftime('%Y-%m-%d') for dt in rrule(MONTHLY, interval=1,bymonthday=(-1),dtstart=parse("20200228T000000"), until=LASTMONTH)]
sliding_windows = set(list(zip(starts,ends)) + sliding_windows)

for window in sliding_windows:
    dag_id = 'bealign_sliding_windows_{}'.format(str('_'.join(window)))
    default_args = {
        'owner': 'sweaver',
        'depends_on_past': False,
        'email': ['sweaver@temple.edu'],
        'email_on_failure': False,
        'email_on_retry': False,
        'params' : {
            'working_dir' : WORKING_DIR,
            'num_procs': 64,
            'region_cfg' : "/data/shares/veg/SARS-CoV-2/SARS-CoV-2-devel/airflow/libs/regions.yaml",
            'python': "/data/shares/veg/SARS-CoV-2/SARS-CoV-2-devel/env/bin/python3",
            'hyphy': "/data/shares/veg/SARS-CoV-2/hyphy/hyphy",
            'hyphy_mpi': "/data/shares/veg/SARS-CoV-2/hyphy/HYPHYMPI",
            'hyphy_lib_path': "/data/shares/veg/SARS-CoV-2/hyphy/res",
            'mafft': "/usr/local/bin/mafft",
            'pre_msa' : "/data/shares/veg/SARS-CoV-2/hyphy-analyses-devel/codon-msa/pre-msa.bf",
            'compressor' : "/data/shares/veg/SARS-CoV-2/SARS-CoV-2-devel/scripts/compressor.bf",
            'compressor2' : "/data/shares/veg/SARS-CoV-2/SARS-CoV-2-devel/scripts/compressor-2.bf",
            'collection-date-range': window,
            'date' : datetime.date.today().strftime('%Y-%m-%d')
        },
        'retries': 1,
        'retry_delay': datetime.timedelta(minutes=5),
        'task_concurrency' : 5,
        # 'on_failure_callback': task_fail_slack_alert,
        # 'on_success_callback': task_success_slack_alert
        # 'queue': 'bash_queue',
        # 'pool': 'backfill',
        # 'priority_weight': 10,
        # 'end_date': datetime(2016, 1, 1),
        # 'wait_for_downstream': False,
        # 'dag': dag,
        # 'sla': timedelta(hours=2),
        # 'execution_timeout': timedelta(seconds=300),
        # 'on_retry_callback': another_function,
        # 'sla_miss_callback': yet_another_function,
        # 'trigger_rule': 'all_success'
    }
    schedule = None
    globals()[dag_id] = create_dag(dag_id,
                                  schedule,
                                  window,
                                  default_args)

