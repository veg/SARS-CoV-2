import yaml
import datetime

# The DAG object; we'll need this to instantiate a DAG
from airflow import DAG

# Operators; we need this to operate!
from airflow.models.baseoperator import cross_downstream, chain
from airflow.operators.bash import BashOperator
from airflow.operators.python import PythonOperator
from airflow.utils.dates import days_ago
from airflow.hooks.base import BaseHook
from airflow.models import Variable

# These args will get passed on to each operator
# You can override them on a per-task basis during operator initialization

from libs.callbacks import task_fail_slack_alert, task_success_slack_alert

import os
import sys
import pathlib

p = os.path.abspath(str(pathlib.Path(__file__).parent.absolute()) + '/../../python/')
if p not in sys.path:
    sys.path.append(p)

from export_sequences import export_sequences, export_premsa_sequences
from export_duplicates import export_duplicates
from export_meta import export_meta
from remove_seq import remove_reference_seq, reserve_only_original_input
from merge_duplicates import merge_duplicates
from fix_duplicates import fix_duplicates
from update_fasta_duplicates import update_fasta_duplicates

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
        'meta-output' : WORKING_DIR + '/master-no-fasta.json',
        'mafft': "/usr/local/bin/mafft",
        'get-latest-by-collection-date': 100000,
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

default_args["params"]["meta-output"] = WORKING_DIR + "/data/fasta/" + default_args["params"]["date"]  + '/master-no-sequences.json'
default_args["params"]["sequence-output"] = WORKING_DIR + "/data/fasta/" + default_args["params"]["date"] + '/sequences'

dag = DAG(
    'initiate_run',
    default_args=default_args,
    description='initiates run',
    schedule_interval='0 11 * * *',
    start_date=datetime.datetime(2021, 2, 10),
    tags=['selection'],
)

with open(dag.params["region_cfg"], 'r') as stream:
    regions = yaml.safe_load(stream)

MAFFT = """
{{ params.mafft }} --auto --thread -1 --mapout --addfragments $INPUT_FN $REFERENCE_FILEPATH >| $TMP_OUTPUT_FN
"""

POSTMSA = """
{{ params.hyphy }} LIBPATH={{params.hyphy_lib_path}} {{ params.post_msa }} --protein-msa $INPUT_FN --nucleotide-sequences $NUC_INPUT_FN --output $COMPRESSED_OUTPUT_FN --duplicates $DUPLICATE_OUTPUT_FN
"""

COMPRESSOR = """
{{ params.hyphy }} LIBPATH={{params.hyphy_lib_path}} {{ params.compressor }} --msa $COMPRESSED_FN --duplicates $DUPLICATE_FN --output $VARIANTS_CSV_FN  --json $VARIANTS_JSON_FN
"""

COMPRESSOR2 = """
{{ params.hyphy }} LIBPATH={{params.hyphy_lib_path}} {{ params.compressor2 }} --msa $COMPRESSED_FN --duplicates $DUPLICATE_FN --csv $VARIANTS_CSV_FN  --byseq $VARIANTS_JSON_FN --p 0.9 --output $FILTERED_FASTA_FN --json FILTERED_JSON_FN
"""


mk_dir = BashOperator(
    task_id='make_directory',
    bash_command='mkdir -p {{params.working_dir}}/data/fasta/{{params.date}}',
    dag=dag,
)

export_meta = PythonOperator(
        task_id='export_meta',
        python_callable=export_meta,
        op_kwargs={ "config" : default_args['params'] },
        dag=dag,
    )

export_sequences = PythonOperator(
        task_id='export_sequences',
        python_callable=export_sequences,
        op_kwargs={ "config" : default_args['params'] },
        dag=dag,
    )

# For each region
export_by_gene = []

for gene in regions.keys():

    reference_filepath = WORKING_DIR + 'reference_genes/reference.' + gene + '_protein.fas'
    filepath_prefix = WORKING_DIR + "/data/fasta/" + default_args["params"]["date"] + '/sequences.' + gene

    nuc_sequence_output = filepath_prefix + '_nuc.fas'
    prot_sequence_output = filepath_prefix + '_protein.fas'

    initial_duplicate_output = filepath_prefix + '.initial.duplicates.json'
    protein_duplicate_output = filepath_prefix + '.protein.duplicates.json'
    duplicate_output = filepath_prefix + '.duplicates.json'
    map_output = filepath_prefix + '.map.json'

    variants_csv_output = filepath_prefix + '.variants.csv'
    variants_json_output = filepath_prefix + '.variants.json'
    filtered_fasta_output = filepath_prefix + '.compressed.filtered.fas'
    filtered_json_output = filepath_prefix + '.filtered.json'

    compressed_output_filepath =  filepath_prefix + '.compressed.fas'

    tree_output = filepath_prefix + '.compressed.filtered.fas.rapidnj.bestTree'
    sto_output = filepath_prefix + '.compressed.filtered.sto';

    tmp_output_fn = filepath_prefix + '.tmp.msa'
    output_fn = filepath_prefix + '.msa'

    slac_output_fn = filepath_prefix + '.SLAC.json'
    fel_output_fn = filepath_prefix + '.FEL.json'
    meme_output_fn = filepath_prefix + '.MEME.json'


    default_args["params"]["nuc-sequence-output"] = nuc_sequence_output
    default_args["params"]["prot-sequence-output"] = prot_sequence_output
    default_args["params"]["duplicate-output"] = duplicate_output
    default_args["params"]["protein-duplicate-output"] = duplicate_output
    default_args["params"]["inital-duplicate-output"] = initial_duplicate_output

    export_premsa_sequence_task = PythonOperator(
            task_id=f'export_premsa_sequences_{gene}',
            python_callable=export_premsa_sequences,
            op_kwargs={ "config" : default_args['params'], 'nuc_output_fn':  nuc_sequence_output, 'prot_output_fn' : prot_sequence_output, 'gene' : gene },
            dag=dag,
        )

    export_duplicates_task = PythonOperator(
        task_id=f'export_duplicates_{gene}',
        python_callable=export_duplicates,
        op_kwargs={ 'output_fn' : initial_duplicate_output, 'gene': gene },
        dag=dag,
    )

    mafft_task = BashOperator(
        task_id=f'mafft_{gene}',
        bash_command=MAFFT,
        params={'mafft': default_args['params']['mafft']},
        env={'INPUT_FN': prot_sequence_output, 'TMP_OUTPUT_FN': tmp_output_fn, 'REFERENCE_FILEPATH': reference_filepath },
        dag=dag
    )

    # input_fn, reference_fn, output_fn
    remove_ref_task = PythonOperator(
        task_id=f'remove_ref_{gene}',
        python_callable=reserve_only_original_input,
        op_kwargs={ "input_fn" : tmp_output_fn, "original_fn" : prot_sequence_output, "output_fn": output_fn },
        dag=dag,
    )

    # Run POST-MSA on cancatenated dataset to translate back to nucleotides
    reverse_translate_task = BashOperator(
        task_id=f'post_msa_{gene}',
        bash_command=POSTMSA,
        env={'INPUT_FN': output_fn, 'NUC_INPUT_FN': nuc_sequence_output , 'COMPRESSED_OUTPUT_FN': compressed_output_filepath, 'DUPLICATE_OUTPUT_FN': protein_duplicate_output, **os.environ},
        dag=dag
    )

    cleanup_task = BashOperator(
        task_id=f'cleanup_{gene}',
        bash_command="sed -i '/^>/! s/[^ACTG-]/N/g' $COMPRESSED_OUTPUT_FN",
        env={'COMPRESSED_OUTPUT_FN': compressed_output_filepath, **os.environ},
        dag=dag
    )

    merge_duplicate_task = PythonOperator(
        task_id=f'merge_duplicates_{gene}',
        python_callable=merge_duplicates,
        op_kwargs={ 'protein_duplicates' : protein_duplicate_output, 'nuc_duplicates': initial_duplicate_output, 'output':  duplicate_output},
        dag=dag,
    )

    # Fix duplicates
    fix_duplicate_task = PythonOperator(
        task_id=f'fix_duplicates_{gene}',
        python_callable=fix_duplicates,
        op_kwargs={ 'duplicates' : duplicate_output, 'map': map_output, 'overwrite': True },
        dag=dag,
    )

    # # Fix header files
    # echo "$PYTHON python/update_fasta_duplicates.py -f ${FILE}.${GENE}.compressed.fas -m ${FILE}.${GENE}.map.json"
    # $PYTHON python/update_fasta_duplicates.py -f ${FILE}.${GENE}.compressed.fas -m ${FILE}.${GENE}.map.json

    update_fasta_duplicates_task = PythonOperator(
        task_id=f'update_fasta_duplicates_{gene}',
        python_callable=update_fasta_duplicates,
        op_kwargs={ 'fasta_file' : compressed_output_filepath, 'map_file': gene },
        dag=dag,
    )

    compressor_task = BashOperator(
        task_id=f'compressor_{gene}',
        bash_command=COMPRESSOR,
        env={'COMPRESSED_FN': compressed_output_filepath, 'DUPLICATE_FN': duplicate_output, 'VARIANTS_CSV_FN': variants_csv_output, 'VARIANTS_JSON_FN': variants_json_output, **os.environ},
        dag=dag
    )

    compressor_two_task = BashOperator(
        task_id=f'compressor_two_{gene}',
        bash_command=COMPRESSOR2,
        env={'COMPRESSED_FN': compressed_output_filepath, 'DUPLICATE_FN': duplicate_output, 'VARIANTS_CSV_FN': variants_csv_output, 'VARIANTS_JSON_FN': variants_json_output, 'FILTERED_FASTA_FN': filtered_fasta_output, 'FILTERED_JSON_FN': filtered_json_output, **os.environ},
        dag=dag
    )

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

    big_data_flags=''
    if(gene == 'nsp3'):
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

    # summarize_gene_task = BashOperator(
    #     task_id='summarize_gene_{gene}',
    #     bash_command='mkdir -p {{params.working_dir}}/data/fasta/{{params.date}}',
    #     dag=dag,
    # )

    selection_flow = export_premsa_sequence_task >> export_duplicates_task >> mafft_task >> remove_ref_task >> reverse_translate_task >> cleanup_task >> merge_duplicate_task >> fix_duplicate_task >> update_fasta_duplicates_task >> compressor_task >> compressor_two_task >> infer_tree_task
    selection_flow.set_downstream(slac_task)
    selection_flow.set_downstream(fel_task)
    selection_flow.set_downstream(meme_task)
    export_by_gene.append(selection_flow)

dag.doc_md = __doc__

export_sequences.doc_md = """\
#### Task Documentation
Creates a directory and exports selected sequences
"""

# Add export meta and export sequence tasks to be executed in parallel
export_by_gene.extend([export_meta, export_sequences])
mk_dir >> export_by_gene


