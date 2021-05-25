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
from airflow.models import Variable

# These args will get passed on to each operator
# You can override them on a per-task basis during operator initialization

from libs.callbacks import dag_fail_slack_alert, dag_success_slack_alert

import os
import sys
import pathlib

p = os.path.abspath(str(pathlib.Path(__file__).parent.absolute()) + '/../../python/')
if p not in sys.path:
    sys.path.append(p)

from export_sequences import export_sequences, export_bealign_sequences
from export_meta import export_meta
from get_raw_duplicates import write_nuc_raw_duplicates

WORKING_DIR = Variable.get("WORKING_DIR")

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
        'get-latest-by-collection-date': 100000,
        'date' : datetime.date.today().strftime('%Y-%m-%d'),
        'python': WORKING_DIR + "/env/bin/python3",
        'hyphy': "/data/shares/veg/SARS-CoV-2/hyphy/hyphy",
        'hyphy_mpi': "/data/shares/veg/SARS-CoV-2/hyphy/HYPHYMPI",
        'hyphy_lib_path': "/data/shares/veg/SARS-CoV-2/hyphy/res",
        'mafft': "/usr/local/bin/mafft",
        'pre_msa' : "/data/shares/veg/SARS-CoV-2/hyphy-analyses-devel/codon-msa/pre-msa.bf",
        'post_msa' : "/data/shares/veg/SARS-CoV-2/hyphy-analyses-devel/codon-msa/post-msa.bf",
        'compressor' : "/data/shares/veg/SARS-CoV-2/SARS-CoV-2-devel/scripts/compressor.bf",
        'compressor2' : "/data/shares/veg/SARS-CoV-2/SARS-CoV-2-devel/scripts/compressor-2.bf",
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

with DAG(
    'initiate_run_last_100k_by_collection_date',
    default_args=default_args,
    description='initiates run with sequences from latest 100k sequences',
    schedule_interval='@weekly',
    start_date=datetime.datetime(2021, 3, 25),
    on_failure_callback=dag_fail_slack_alert,
    on_success_callback=dag_success_slack_alert,
    tags=['selection'],
    ) as dag:

    OUTPUT_DIR = WORKING_DIR + "/data/fasta/last100k-by-collection-date"
    default_args["params"]["output-dir"] = OUTPUT_DIR
    default_args["params"]["meta-output"] = OUTPUT_DIR + '/master-no-sequences.json'
    default_args["params"]["sequence-output"] = OUTPUT_DIR + '/sequences'

    with open(dag.params["region_cfg"], 'r') as stream:
        regions = yaml.safe_load(stream)

    def my_task():
        context = get_current_context()
        ds = context["ds"]
        OUTPUT_DIR = WORKING_DIR + "/data/fasta/last100k-by-collection-date"
        default_args["params"]["output-dir"] = OUTPUT_DIR
        default_args["params"]["meta-output"] = OUTPUT_DIR + '/master-no-sequences.json'
        default_args["params"]["sequence-output"] = OUTPUT_DIR + '/sequences'

    test_task = PythonOperator(
            task_id=f'test_task',
            python_callable=my_task,
            dag=dag,
        )

    mk_dir_task = BashOperator(
        task_id='make_directory',
        bash_command='mkdir -p {{params.output}}',
        params={'output': default_args['params']['output-dir']},
        dag=dag,
    )

    mk_dir_task.set_upstream(test_task)

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

        filepath_prefix = OUTPUT_DIR + '/sequences.' + gene

        nuc_sequence_output = filepath_prefix + '_nuc.fas'
        uniques_fn = filepath_prefix + '_nuc.uniques.fas'
        duplicate_output = filepath_prefix + '.duplicates.json'

        variants_csv_output = filepath_prefix + '.variants.csv'
        variants_json_output = filepath_prefix + '.variants.json'
        filtered_fasta_output = filepath_prefix + '.compressed.filtered.fas'
        filtered_json_output = filepath_prefix + '.filtered.json'
        output_edits_fn = filepath_prefix + '.filtered.edits.json'

        tn93_cluster_fn = filepath_prefix + '.tn93.cluster.json'
        centroid_fn = filepath_prefix + '.centroids.fas'

        compressor_duplicate_out = filepath_prefix + '.duplicates.variants.json'

        tree_output = filepath_prefix + '.compressed.filtered.fas.rapidnj.bestTree'
        sto_output = filepath_prefix + '.compressed.filtered.sto';

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

        export_bealign_task >> cleanup_task

    with TaskGroup(f"duplicates_{gene}") as duplicates_group:

        compute_duplicates_task = PythonOperator(
            task_id=f'write_raw_duplicates',
            python_callable=write_nuc_raw_duplicates,
            op_kwargs={ "input" : nuc_sequence_output, "duplicate_output" : duplicate_output, 'uniques_output': uniques_fn },
            dag=dag,
        )

        compute_duplicates_task

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

        # bpsh 0 tn93-cluster -t 0.001 sequences.S.compressed.filtered.fas
        TN93_CLUSTER = """
        {{ params.tn93_cluster }} -t {{ params.threshold }} -o {{ params.output_fn }} {{ params.input_fn }}
        """
        tn93_cluster_task = BashOperator(
            task_id=f'tn93_cluster',
            bash_command=TN93_CLUSTER,
            params={'input_fn': filtered_fasta_output, 'threshold': str(regions[gene]['cluster_threshold']), 'output_fn': tn93_cluster_fn},
            dag=dag
        )

        WRITE_CENTROIDS= """
        cat {{ params.input_fn }} | /usr/local/bin/jq -r '.[].centroid' > {{ params.output_fn }}
        """
        write_centroids_task = BashOperator(
            task_id=f'write_centroids',
            bash_command=WRITE_CENTROIDS,
            params={ 'input_fn': tn93_cluster_fn, 'output_fn': centroid_fn },
            dag=dag
        )


        compressor_task >> compressor_two_task >> tn93_cluster_task >> write_centroids_task

    INFER_TREE = """
    seqmagick convert $FILTERED_FASTA_FN $STO_OUTPUT;
    rapidnj $STO_OUTPUT -i sth > $TREE_OUTPUT
    sed -i "s/'//g" $TREE_OUTPUT;
    """

    infer_tree_task = BashOperator(
        task_id=f'infer_tree_{gene}',
        bash_command=INFER_TREE,
        env={'FILTERED_FASTA_FN': centroid_fn, 'STO_OUTPUT': sto_output, 'TREE_OUTPUT': tree_output, **os.environ},
        dag=dag
    )

    slac_task = BashOperator(
        task_id=f'slac_{gene}',
        bash_command="{{ params.hyphy }} LIBPATH={{params.hyphy_lib_path}} slac --kill-zero-lengths Constrain ENV='_DO_TREE_REBALANCE_=1' --alignment $FILTERED_FASTA_FN --tree $TREE_OUTPUT --branches All --samples 0 --output $SLAC_OUTPUT",
        env={'FILTERED_FASTA_FN': centroid_fn, 'TREE_OUTPUT': tree_output, 'SLAC_OUTPUT': slac_output_fn, **os.environ},
        pool='hyphy',
        dag=dag,
    )

    big_data_flags='--full-model No'

    fel_task = BashOperator(
        task_id=f'fel_{gene}',
        bash_command="{{ params.hyphy }} LIBPATH={{params.hyphy_lib_path}} fel --kill-zero-lengths Constrain ENV='_DO_TREE_REBALANCE_=1' $BIG_DATA_FLAGS --alignment $FILTERED_FASTA_FN --tree $TREE_OUTPUT --branches Internal --output $FEL_OUTPUT",
        env={'FILTERED_FASTA_FN': centroid_fn, 'TREE_OUTPUT': tree_output, 'FEL_OUTPUT': fel_output_fn, 'BIG_DATA_FLAGS': big_data_flags, **os.environ},
        pool='hyphy',
        dag=dag,
    )

    meme_task = BashOperator(
        task_id=f'meme_{gene}',
        bash_command="{{ params.hyphy }} LIBPATH={{params.hyphy_lib_path}} meme --kill-zero-lengths Constrain ENV='_DO_TREE_REBALANCE_=1' $BIG_DATA_FLAGS --alignment $FILTERED_FASTA_FN --tree $TREE_OUTPUT --branches Internal --output $MEME_OUTPUT",
        env={'FILTERED_FASTA_FN': centroid_fn, 'TREE_OUTPUT': tree_output, 'MEME_OUTPUT': meme_output_fn, 'BIG_DATA_FLAGS': big_data_flags, **os.environ},
        pool='hyphy',
        dag=dag,
    )

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
    # cross_downstream([export_meta_task, export_sequences_task], export_by_gene)
    cross_downstream([export_meta_task, export_sequences_task], export_by_gene)
