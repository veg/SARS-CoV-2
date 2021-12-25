import os
import sys
import pathlib

import yaml
import json
import datetime
from dateutil.relativedelta import *
from dateutil.rrule import *
from dateutil.parser import *

# The DAG object; we'll need this to instantiate a DAG
from airflow import DAG

# Operators; we need this to operate!
from airflow.utils.task_group import TaskGroup
from airflow.operators.bash import BashOperator
from airflow.operators.python import PythonOperator
from airflow.models.baseoperator import cross_downstream, chain
from airflow.utils.dates import days_ago
from airflow.hooks.base import BaseHook
from airflow.models import Variable
from airflow.utils.trigger_rule import TriggerRule

# These args will get passed on to each operator
# You can override them on a per-task basis during operator initialization
from libs.callbacks import dag_fail_slack_alert, dag_success_slack_alert


p = os.path.abspath(str(pathlib.Path(__file__).parent.absolute()) + "/../../python/")

if p not in sys.path:
    sys.path.append(p)

from export_sequences import export_sequences, export_bealign_sequences
from summarize_rascl import rascl_summary_report
from export_meta import export_meta
from get_raw_duplicates import write_nuc_raw_duplicates

WORKING_DIR = Variable.get("WORKING_DIR")

def export_with_context(func, config, **context):
    # Get last 90 days
    execution_date = datetime.datetime.strptime(context['ds'], '%Y-%m-%d')
    ninety_days_from_execution_date = execution_date-relativedelta(days=+90)
    window = (ninety_days_from_execution_date.strftime('%Y-%m-%d'), execution_date.strftime('%Y-%m-%d'))
    config['collection-date-range'] = window
    func(config)

def export_bealign_with_context(func, config, nuc_output_fn, gene, **context):
    # Get last 90 days
    execution_date = datetime.datetime.strptime(context['ds'], '%Y-%m-%d')
    ninety_days_from_execution_date = execution_date-relativedelta(days=+90)
    window = (ninety_days_from_execution_date.strftime('%Y-%m-%d'), execution_date.strftime('%Y-%m-%d'))
    config['collection-date-range'] = window
    func(config, nuc_output_fn, gene)

def create_dag(dag_id, schedule, clade, default_args):
    with DAG(
        dag_id,
        default_args=default_args,
        description='RASCL',
        schedule_interval=schedule,
        start_date=datetime.datetime(2021, 8, 1),
        on_failure_callback=dag_fail_slack_alert,
        on_success_callback=dag_success_slack_alert,
        tags=['selection','clade'],
        ) as dag:

        last_exec_date = dag.get_latest_execution_date()

        if last_exec_date is None:
            last_exec_date = datetime.datetime(year=1970, month=1, day=1)

        unique_id = str(round(last_exec_date.timestamp()))
        priority = 100

        OUTPUT_DIR = WORKING_DIR + '/data/rascl/' + clade.strip() + '/' + unique_id
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
                python_callable=lambda config, **context: export_with_context(export_meta,config, **context),
                op_kwargs={ "config" : default_args['params'] },
                provide_context=True,
                pool='mongo',
                dag=dag,
            )

        export_meta_task.set_upstream(mk_dir_task)

        export_sequences_task = PythonOperator(
                task_id='export_sequences',
                python_callable=lambda config, **context: export_with_context(export_sequences,config, **context),
                op_kwargs={ "config" : default_args['params'] },
                provide_context=True,
                pool='mongo',
                dag=dag,
            )

        export_sequences_task.set_upstream(mk_dir_task)

        # For each region
        export_by_gene = []

        for gene in regions.keys():

            filepath_prefix = OUTPUT_DIR + '/sequences.' + gene

            nuc_sequence_output = filepath_prefix + '_nuc.fas'
            fasta_no_ambigs_fn = filepath_prefix + '.no_ambigs.fas'
            uniques_fn = filepath_prefix + '_nuc.uniques.fas'
            duplicate_output = filepath_prefix + '.duplicates.json'


            fade_output_fn = filepath_prefix + '.FADE.json'
            fel_output_fn = filepath_prefix + '.FEL.json'
            meme_output_fn = filepath_prefix + '.MEME.json'
            meme_full_output_fn = filepath_prefix + '.FULL.MEME.json'
            prime_output_fn = filepath_prefix + '.PRIME.json'
            absrel_output_fn = filepath_prefix + '.ABSREL.json'
            busted_output_fn = filepath_prefix + '.BUSTED.json'
            relax_output_fn = filepath_prefix + '.RELAX.json'
            cfel_output_fn = filepath_prefix + '.CFEL.json'
            slac_output_fn = filepath_prefix + '.SLAC.json'
            bgm_output_fn = filepath_prefix + '.BGM.json'
            summary_output_fn = filepath_prefix + '.json'

            tn93_cluster_fn = filepath_prefix + '.tn93.cluster.json'
            centroid_fn = filepath_prefix + '.centroids.fas'

            ref_combined_fn = filepath_prefix + '.combined.fas'
            ref_combined_protein_fn = filepath_prefix + '.combined.prot.fas'
            tree_output = ref_combined_fn + '.raxml.bestTree'

            out_int_tree_path = os.path.join(OUTPUT_DIR, "sequences." + gene + ".int.nwk")
            out_clade_tree_path = os.path.join(OUTPUT_DIR, "sequences." + gene + ".clade.nwk")
            out_full_tree_path = os.path.join(OUTPUT_DIR, "sequences." + gene + ".full.nwk")

            input_ref_fn =  WORKING_DIR + "/rascl_reference/processed/" + gene + ".reference.compressed.fas"
            in_gene_ref_seq =  default_args["params"]["rascl_path"] + "/data/ReferenceSeq/" + gene + ".fas"


            with TaskGroup(f"alignment_{gene}") as alignment:

                export_bealign_task = PythonOperator(
                    task_id=f'export_bealign',
                    python_callable=lambda config, nuc_output_fn, gene, **context: export_bealign_with_context(export_bealign_sequences, config, nuc_output_fn, gene, **context),
                    op_kwargs={ "config" : default_args['params'], 'nuc_output_fn':  nuc_sequence_output, 'gene' : gene },
                    provide_context=True,
                    pool='mongo',
                    dag=dag,
                )

                # Occasional errors when cleaning up tmp files, so or'ing true
                cleanup_task = BashOperator(
                    task_id=f'cleanup',
                    bash_command="sed -i '/^>/! s/[^ACTG-]/N/g' $NUC_OUTPUT_FN || true",
                    env={'NUC_OUTPUT_FN': nuc_sequence_output, **os.environ},
                    dag=dag
                )

                strike_ambigs_task = BashOperator(
                    task_id=f'strike_ambigs',
                    bash_command="{{ params.hyphy }} LIBPATH={{params.hyphy_lib_path}} {{ params.strike_ambigs_path }} --alignment {{ params.input_fn }} --output {{ params.output_fn }}",
                    params={'input_fn': nuc_sequence_output, 'output_fn': fasta_no_ambigs_fn },
                    execution_timeout=datetime.timedelta(minutes=5),
                    dag=dag
                )

                export_bealign_task >> cleanup_task >> strike_ambigs_task

            # Use tn93-cluster
            with TaskGroup(f"sampling_{gene}") as sampling:

                TN93_CLUSTER = """
                {{ params.tn93_cluster }} -t {{ params.threshold_query }} -o {{ params.output_fn }} {{ params.input_fn }}
                """

                tn93_cluster_task = BashOperator(
                    task_id=f'tn93_cluster',
                    bash_command=TN93_CLUSTER,
                    params={'input_fn': fasta_no_ambigs_fn, 'output_fn': tn93_cluster_fn},
                    priority_weight=priority,
                    dag=dag
                )

                WRITE_CENTROIDS= """
                cat {{ params.input_fn }} | /usr/local/bin/jq -r '.[].centroid' > {{ params.output_fn }}
                """

                write_centroids_task = BashOperator(
                    task_id=f'write_centroids',
                    bash_command=WRITE_CENTROIDS,
                    params={ 'input_fn': tn93_cluster_fn, 'output_fn': centroid_fn },
                    priority_weight=priority,
                    dag=dag
                )

                write_centroids_task.set_upstream(tn93_cluster_task)

            combine_task = BashOperator(
                task_id=f'combine_{gene}',
                bash_command="python3 {{ params.combine_path }} --input {{ params.input_fn }} -o {{ params.output_fn}} --threshold {{ params.threshold_query }} --msa {{ params.input_ref_fn }} --reference_seq {{ params.in_gene_ref_seq }}",
                params={'input_fn': centroid_fn, 'output_fn': ref_combined_fn, 'input_ref_fn' : input_ref_fn, 'in_gene_ref_seq' : in_gene_ref_seq },
                env={**os.environ},
                dag=dag
            )

            protein_conv_task = BashOperator(
                task_id=f'protein_conv_{gene}',
                bash_command="{{ params.hyphy }} LIBPATH={{params.hyphy_lib_path}} conv Universal 'Keep Deletions' {{ params.input_fn }} {{ params.output_fn }}",
                params={'input_fn': ref_combined_fn, 'output_fn': ref_combined_protein_fn },
                dag=dag
            )

            INFER_TREE = """
            raxml-ng --model GTR --msa {{ params.combined_fas }} --threads {{ params.THREADS }} --tree pars{3} --force
            """

            infer_tree_task = BashOperator(
                task_id=f'infer_tree_{gene}',
                bash_command=INFER_TREE,
                params={'combined_fas': ref_combined_fn, 'tree_output': tree_output},
                execution_timeout=datetime.timedelta(hours=24),
                dag=dag
            )

            annotation_task = BashOperator(
                task_id=f'annotation_{gene}',
                bash_command="bash {{ params.annotate_path }} {{ params.input_fn }} 'REFERENCE' {{ params.in_compressed_fas }} {{ params.LABEL }} {{ params.output_dir }}",
                params={'input_fn': tree_output, 'in_compressed_fas': centroid_fn, 'output_dir': OUTPUT_DIR },
                execution_timeout=datetime.timedelta(hours=24),
                dag=dag
            )

            with TaskGroup(f"selection_analyses_{gene}") as selection_analyses:

                big_data_flags='--full-model No'

                fade_task = BashOperator(
                    task_id=f'fade',
                    bash_command="{{ params.hyphy }} LIBPATH={{params.hyphy_lib_path}} FADE --alignment $FILTERED_FASTA_FN  --tree $TREE_OUTPUT --output $FADE_OUTPUT --branches {{ params.LABEL }}",
                    env={'FILTERED_FASTA_FN': ref_combined_protein_fn, 'TREE_OUTPUT': out_clade_tree_path, 'FADE_OUTPUT': fade_output_fn, 'BIG_DATA_FLAGS': big_data_flags, **os.environ},
                    pool='hyphy',
                    dag=dag,
                )

                fel_task = BashOperator(
                    task_id=f'fel',
                    bash_command="{{ params.hyphy }} LIBPATH={{params.hyphy_lib_path}} fel --kill-zero-lengths Constrain ENV='_DO_TREE_REBALANCE_=1' $BIG_DATA_FLAGS --alignment $FILTERED_FASTA_FN --tree $TREE_OUTPUT --branches Internal --output $FEL_OUTPUT",
                    env={'FILTERED_FASTA_FN': ref_combined_fn, 'TREE_OUTPUT': out_int_tree_path, 'FEL_OUTPUT': fel_output_fn, 'BIG_DATA_FLAGS': big_data_flags, **os.environ},
                    pool='hyphy',
                    dag=dag,
                )

                meme_task = BashOperator(
                    task_id=f'meme',
                    bash_command="{{ params.hyphy }} LIBPATH={{params.hyphy_lib_path}} meme --kill-zero-lengths Constrain ENV='_DO_TREE_REBALANCE_=1' $BIG_DATA_FLAGS --alignment $FILTERED_FASTA_FN --tree $TREE_OUTPUT --branches Internal --output $MEME_OUTPUT",
                    env={'FILTERED_FASTA_FN': ref_combined_fn, 'TREE_OUTPUT': out_int_tree_path, 'MEME_OUTPUT': meme_output_fn, 'BIG_DATA_FLAGS': big_data_flags, **os.environ},
                    pool='hyphy',
                    dag=dag,
                )

                meme_full_task = BashOperator(
                    task_id=f'meme_full',
                    bash_command="{{ params.hyphy }} LIBPATH={{params.hyphy_lib_path}} MEME --alignment $FILTERED_FASTA_FN --tree $TREE_OUTPUT --output $MEME_OUTPUT --branches {{ params.LABEL }}",
                    env={'FILTERED_FASTA_FN': ref_combined_fn, 'TREE_OUTPUT': out_full_tree_path, 'MEME_OUTPUT': meme_full_output_fn, 'BIG_DATA_FLAGS': big_data_flags, **os.environ},
                    pool='hyphy',
                    dag=dag,
                )

                prime_task = BashOperator(
                    task_id=f'prime',
                    bash_command="{{ params.hyphy }} LIBPATH={{params.hyphy_lib_path}} PRIME --alignment $FILTERED_FASTA_FN --tree $TREE_OUTPUT  --output $PRIME_OUTPUT --branches {{ params.LABEL }}",
                    env={'FILTERED_FASTA_FN': ref_combined_fn, 'TREE_OUTPUT': out_int_tree_path, 'PRIME_OUTPUT': prime_output_fn, 'BIG_DATA_FLAGS': big_data_flags, **os.environ},
                    pool='hyphy',
                    dag=dag,
                )

                contrast_fel_task = BashOperator(
                    task_id=f'cfel',
                    bash_command="{{ params.hyphy }} LIBPATH={{params.hyphy_lib_path}} contrast-fel --alignment $FILTERED_FASTA_FN --tree $TREE_OUTPUT  --output $CFEL_OUTPUT --branch-set {{ params.LABEL }} --branch-set Reference",
                    env={'FILTERED_FASTA_FN': ref_combined_fn, 'TREE_OUTPUT': out_clade_tree_path, 'CFEL_OUTPUT': cfel_output_fn, 'BIG_DATA_FLAGS': big_data_flags, **os.environ},
                    pool='hyphy',
                    dag=dag,
                )

                relax_task = BashOperator(
                    task_id=f'relax',
                    bash_command="{{ params.hyphy }} LIBPATH={{params.hyphy_lib_path}} RELAX --alignment $FILTERED_FASTA_FN --models Minimal --tree $TREE_OUTPUT --output $RELAX_OUTPUT --test {{ params.LABEL }} --reference Reference --starting-points 10 --srv Yes",
                    env={'FILTERED_FASTA_FN': ref_combined_fn, 'TREE_OUTPUT': out_clade_tree_path, 'RELAX_OUTPUT': relax_output_fn, 'BIG_DATA_FLAGS': big_data_flags, **os.environ},
                    pool='hyphy',
                    dag=dag,
                )

                absrel_task = BashOperator(
                    task_id=f'absrel',
                    bash_command="{{ params.hyphy }} LIBPATH={{params.hyphy_lib_path}} ABSREL --alignment $FILTERED_FASTA_FN --tree $TREE_OUTPUT --output $ABSREL_OUTPUT --branches {{ params.LABEL }}",
                    env={'FILTERED_FASTA_FN': ref_combined_fn, 'TREE_OUTPUT': out_int_tree_path, 'ABSREL_OUTPUT': absrel_output_fn, 'BIG_DATA_FLAGS': big_data_flags, **os.environ},
                    pool='hyphy',
                    dag=dag,
                )

                busted_task = BashOperator(
                    task_id=f'busted',
                    bash_command="{{ params.hyphy }} LIBPATH={{params.hyphy_lib_path}} BUSTED --alignment $FILTERED_FASTA_FN --tree $TREE_OUTPUT --output $BUSTED_OUTPUT --branches {{ params.LABEL }} --starting-points 10",
                    env={'FILTERED_FASTA_FN': ref_combined_fn, 'TREE_OUTPUT': out_clade_tree_path, 'BUSTED_OUTPUT': busted_output_fn, 'BIG_DATA_FLAGS': big_data_flags, **os.environ},
                    pool='hyphy',
                    dag=dag,
                )

                bgm_task = BashOperator(
                    task_id=f'bgm',
                    bash_command="{{ params.hyphy }} LIBPATH={{params.hyphy_lib_path}} BGM --alignment $FILTERED_FASTA_FN --tree $TREE_OUTPUT --output $BGM_OUTPUT --branches {{ params.LABEL }}",
                    env={'FILTERED_FASTA_FN': ref_combined_fn, 'TREE_OUTPUT': out_int_tree_path, 'BGM_OUTPUT': bgm_output_fn, 'BIG_DATA_FLAGS': big_data_flags, **os.environ},
                    pool='hyphy',
                    dag=dag,
                )

                slac_task = BashOperator(
                    task_id=f'slac',
                    bash_command="{{ params.hyphy }} LIBPATH={{params.hyphy_lib_path}} slac --kill-zero-lengths Constrain ENV='_DO_TREE_REBALANCE_=1' --alignment $FILTERED_FASTA_FN --tree $TREE_OUTPUT --branches All --samples 0 --output $SLAC_OUTPUT",
                    env={'FILTERED_FASTA_FN': ref_combined_fn, 'TREE_OUTPUT': out_int_tree_path, 'SLAC_OUTPUT': slac_output_fn, **os.environ},
                    pool='hyphy',
                    dag=dag,
                )

                [fade_task, fel_task, meme_task, meme_full_task, prime_task, contrast_fel_task, relax_task, absrel_task, busted_task, bgm_task, slac_task]

            alignment.set_upstream(export_sequences_task)
            export_by_gene.append(alignment >> sampling >> combine_task >> protein_conv_task >> infer_tree_task >> annotation_task >> selection_analyses)

        with TaskGroup(f"summary_reports") as summary_reports:

           simple_summary_report_task = PythonOperator(
                task_id=f'simple_summary_report',
                trigger_rule=TriggerRule.ALL_DONE,
                python_callable=rascl_summary_report,
                op_kwargs={ "input_dir" : OUTPUT_DIR, 'output_fn':  OUTPUT_DIR + '/report.csv' },
                provide_context=True,
                dag=dag,
            )

           simple_summary_report_task

        with TaskGroup(f"release") as release:

            mk_release_dir_task = BashOperator(
                task_id='make_directory',
                trigger_rule=TriggerRule.ALL_DONE,
                bash_command='mkdir -p /data/shares/web/web/covid-19/selection-analyses/rascl/{{ params.clade }}/{{ run_id }}/',
                params={'clade': clade },
                dag=dag,
            )

            copy_results_task = BashOperator(
                task_id=f'copy_results',
                trigger_rule=TriggerRule.ALL_DONE,
                bash_command='cp {{params.output}}/sequences.*.json {{params.output}}/report.csv /data/shares/web/web/covid-19/selection-analyses/rascl/{{ params.clade }}/{{ run_id }}/',
                params={'output': default_args['params']['output-dir'], 'clade': clade.strip() },
                dag=dag,
            )

            mk_release_dir_task >> copy_results_task

        export_by_gene >> summary_reports >> release

        dag.doc_md = __doc__

        # Add export meta and export sequence tasks to be executed in parallel
        # cross_downstream([export_meta_task, export_sequences_task], export_by_gene)

        return dag

clades = [
    "B.1.2",
    "B.1.596",
    "B.1",
    "B.1.1.519",
    "B.1.243",
    "B.1.234",
    "B.1.526.1",
    "B.1.1",
    "B.1.1.529",
    "B.1.526.2",
    "B.1.575",
    "C.37",
    "R.1",
    "B.1.1.7",
    "B.1.429",
    "B.1.427",
    "B.1.351",
    "B.1.351.2",
    "B.1.351.3",
    "P.1",
    "P.1.1",
    "P.1.2",
    "B.1.526",
    "P.2",
    "B.1.525",
    "B.1.617",
    "B.1.617.1",
    "B.1.617.2",
    "B.1.617.3",
	"AY.1",
	"AY.2",
	"AY.3",
	"AY.4",
	"AY.5",
	"AY.6",
	"AY.7",
	"AY.8",
	"AY.9",
	"AY.10",
	"AY.11",
	"AY.12"
    ]

# Add VOCs from WHO config
with open(WORKING_DIR + "/airflow/libs/voc.json", 'r') as stream:
    vocs = json.load(stream)

clades = list(set(vocs['clades']).union(clades))

for clade in clades:
    dag_id = 'rascl_{}'.format(clade)
    default_args = {
        'owner': 'sweaver',
        'depends_on_past': False,
        'email': ['sweaver@temple.edu'],
        'email_on_failure': False,
        'email_on_retry': False,
        'params' : {
            'working_dir' : WORKING_DIR,
            'num_procs': 64,
            'region_cfg' : WORKING_DIR + "/airflow/libs/regions.yaml",
            'python': WORKING_DIR + "/env/bin/python3",
            'hyphy': "/data/shares/veg/SARS-CoV-2/hyphy/hyphy",
            'hyphy_mpi': "/data/shares/veg/SARS-CoV-2/hyphy/HYPHYMPI",
            'hyphy_lib_path': "/data/shares/veg/SARS-CoV-2/hyphy/res",
            'tn93_cluster': "/usr/local/bin/tn93-cluster",
            'rascl_path': "/data/shares/veg/SARS-CoV-2/SARS-CoV-2_Clades",
            'annotate_path': "/data/shares/veg/SARS-CoV-2/SARS-CoV-2_Clades/scripts/annotate.sh",
            'combine_path': "/data/shares/veg/SARS-CoV-2/SARS-CoV-2_Clades/scripts/combine.py",
            'strike_ambigs_path': "/data/shares/veg/SARS-CoV-2/SARS-CoV-2_Clades/scripts/strike-ambigs.bf",
            "LABEL":"TEST",
            "GISAID_WG":"TEST.fasta",
            "max_ref":"200",
            "max_query":"500",
            "threshold_query":"0.0005",
            "threshold_ref":"0.001",
            'THREADS' : 8,
            'clades' : [clade],
            'date' : datetime.date.today().strftime('%Y-%m-%d')
        },
        'retries': 1,
        'retry_delay': datetime.timedelta(minutes=5),
        'task_concurrency' : 5,
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
    schedule = '@monthly'
    globals()[dag_id] = create_dag(dag_id,
                                  schedule,
                                  clade,
                                  default_args)

