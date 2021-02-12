import yaml
import datetime
from datetime import timedelta

# The DAG object; we'll need this to instantiate a DAG
from airflow import DAG

# Operators; we need this to operate!
from airflow.operators.bash import BashOperator
from airflow.operators.python import PythonOperator
from airflow.utils.dates import days_ago
# These args will get passed on to each operator
# You can override them on a per-task basis during operator initialization

import os
import sys
import pathlib

p = os.path.abspath(str(pathlib.Path(__file__).parent.absolute()) + '/../../python/')
if p not in sys.path:
    sys.path.append(p)


from export_sequences_without_premsa import export_sequences
from store_premsa import store_premsa_file
from premsa_log_parse import mark_troubled

WORKING_DIR = "/data/shares/veg/SARS-CoV-2/SARS-CoV-2-devel/"

default_args = {
    'owner': 'sweaver',
    'depends_on_past': False,
    'email': ['sweaver@temple.edu'],
    'email_on_failure': False,
    'email_on_retry': False,
    'params' : {
        'working_dir' : WORKING_DIR,
        'num_procs': 16,
        'python': "/data/shares/veg/SARS-CoV-2/SARS-CoV-2-devel/env/bin/python3",
        'hyphy': "/data/shares/veg/SARS-CoV-2/hyphy/hyphy",
        'hyphy_mpi': "/data/shares/veg/SARS-CoV-2/hyphy/HYPHYMPI",
        'hyphy_lib_path': "/data/shares/veg/SARS-CoV-2/hyphy/res",
        'pre_msa' : "/data/shares/veg/SARS-CoV-2/hyphy-analyses/codon-msa/pre-msa.bf",
        'compressor' : "/data/shares/veg/SARS-CoV-2/SARS-CoV-2-devel/scripts/compressor.bf",
        'compressor2' : "/data/shares/veg/SARS-CoV-2/SARS-CoV-2-devel/scripts/compressor-2.bf",
        'region_cfg' : "/data/shares/veg/SARS-CoV-2/SARS-CoV-2-devel/airflow/libs/regions.yaml",
        'zero_length_flags' : '--kill-zero-lengths Constrain ENV="_DO_TREE_REBALANCE_=1"'
    },
    'retries': 3,
    'retry_delay': timedelta(minutes=5),
    'concurrency': 10,
    # 'execution_timeout': timedelta(minutes=30),
    # 'queue': 'bash_queue',
    # 'pool': 'backfill',
    # 'priority_weight': 10,
    # 'end_date': datetime(2016, 1, 1),
    # 'wait_for_downstream': False,
    # 'dag': dag,
    # 'sla': timedelta(hours=2),
    # 'on_failure_callback': some_function,
    # 'on_success_callback': some_other_function,
    # 'on_retry_callback': another_function,
    # 'sla_miss_callback': yet_another_function,
    # 'trigger_rule': 'all_success'
}

dag = DAG(
    'populate_pre_msa',
    default_args=default_args,
    description='performs selection analysis',
    schedule_interval=timedelta(days=1),
    start_date=days_ago(2),
    tags=['selection'],
)

with open(dag.params["region_cfg"], 'r') as stream:
    regions = yaml.safe_load(stream)

PREMSA = """
bpsh {{ params.node }} mpirun -np {{ params.num_procs }} {{ params.hyphy_mpi }} LIBPATH={{ params.hyphy_lib_path}} {{ params.pre_msa }} --input {{ params.filepath }} --reference {{ params.working_dir }}/{{ params.regions[params["gene"]]["reference"] }} --trim-from {{ params.regions[params.gene]["trim_from"] }} --trim-to {{ params.regions[params.gene]["trim_to"] }} --E 0.01 --N-fraction {{ params.regions[params["gene"]]["fraction"] }} --remove-stop-codons Yes > {{ params.stdout }}
"""

DATE_STRING = datetime.date.today().strftime('%Y-%m-%d')

pre_msa_tasks = []
i = 0

for gene in regions.keys():

    filepath = WORKING_DIR + 'data/premsa-processor/' + gene + '/sequences.' + DATE_STRING + '.fasta'
    stdout = WORKING_DIR + 'data/premsa-processor/' + gene + '/sequences.' + DATE_STRING + '.stdout.log'

    export_missing = PythonOperator(
        task_id=f'export_missing_premsa_{gene}',
        python_callable=export_sequences,
        op_kwargs={ "gene" : gene, "output_fn" : filepath },
        dag=dag,
    )

    pre_msa = BashOperator(
        task_id=f'pre_msa_{gene}',
        bash_command=PREMSA,
        params={'regions': regions, 'filepath': filepath, 'gene': gene, 'node' : i % 8, 'stdout' : stdout },
        dag=dag,
    )

    # Store nuc_input, prot_input, type
    import_premsa_seqs = PythonOperator(
        task_id=f'store_premsa_{gene}',
        python_callable=store_premsa_file,
        op_kwargs={ "nuc_input" : filepath, "prot_input" : filepath, "gene": gene },
        dag=dag,
    )

    mark_troubled_task = PythonOperator(
        task_id=f'mark_troubled_{gene}',
        python_callable=mark_troubled,
        op_kwargs={ "log_file" : stdout, "gene": gene },
        dag=dag,
    )

    i += 1
    pre_msa_tasks.append(export_missing >> pre_msa >> [import_premsa_seqs, mark_troubled_task])

dag.doc_md = __doc__
pre_msa_tasks
