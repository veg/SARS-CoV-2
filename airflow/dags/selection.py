import yaml
from datetime import timedelta

# The DAG object; we'll need this to instantiate a DAG
from airflow import DAG

# Operators; we need this to operate!
from airflow.operators.bash import BashOperator
from airflow.utils.dates import days_ago
# These args will get passed on to each operator
# You can override them on a per-task basis during operator initialization

WORKING_DIR = "/data/shares/veg/SARS-CoV-2/SARS-CoV-2-devel/"

default_args = {
    'owner': 'sweaver',
    'depends_on_past': False,
    'email': ['sweaver@temple.edu'],
    'email_on_failure': False,
    'email_on_retry': False,
    'params' : {
        'working_dir' : WORKING_DIR,
        'num_procs': 64,
        'python': "/data/shares/veg/SARS-CoV-2/SARS-CoV-2-devel/env/bin/python3",
        'hyphy': "/data/shares/veg/SARS-CoV-2/hyphy/hyphy",
        'hyphy_mpi': "/data/shares/veg/SARS-CoV-2/hyphy/HYPHYMPI",
        'hyphy_lib_path': "/data/shares/veg/SARS-CoV-2/hyphy/res",
        'pre_msa' : "/data/shares/veg/SARS-CoV-2/hyphy-analyses/codon-msa/pre-msa.bf",
        'compressor' : "/data/shares/veg/SARS-CoV-2/SARS-CoV-2-devel/scripts/compressor.bf",
        'compressor2' : "/data/shares/veg/SARS-CoV-2/SARS-CoV-2-devel/scripts/compressor-2.bf",
        'region_cfg' : "/data/shares/veg/SARS-CoV-2/SARS-CoV-2-devel/airflow/libs/regions.yaml",
        'gene': 'leader',
        'trim_from' : 1,
        'trim_to' : 1000,
        'reference' :"/data/shares/veg/SARS-CoV-2/SARS-CoV-2-devel/reference_genes/leader.fas",
        'fraction' : 0.001,
        'input': '/data/shares/veg/SARS-CoV-2/SARS-CoV-2-devel/data/fasta/2020-04-01/sequences',
        'zero_length_flags' : '--kill-zero-lengths Constrain ENV="_DO_TREE_REBALANCE_=1"'
    },
    'retries': 1,
    'retry_delay': timedelta(minutes=5),
    # 'queue': 'bash_queue',
    # 'pool': 'backfill',
    # 'priority_weight': 10,
    # 'end_date': datetime(2016, 1, 1),
    # 'wait_for_downstream': False,
    # 'dag': dag,
    # 'sla': timedelta(hours=2),
    # 'execution_timeout': timedelta(seconds=300),
    # 'on_failure_callback': some_function,
    # 'on_success_callback': some_other_function,
    # 'on_retry_callback': another_function,
    # 'sla_miss_callback': yet_another_function,
    # 'trigger_rule': 'all_success'
}

dag = DAG(
    'selection',
    default_args=default_args,
    description='performs selection analysis',
    schedule_interval=timedelta(days=1),
    start_date=days_ago(2),
    tags=['selection'],
)

with open(dag.params["region_cfg"], 'r') as stream:
    regions = yaml.safe_load(stream)

PREMSA = """
mpirun -np {{ params.num_procs }} {{ params.hyphy_mpi }} LIBPATH={{ params.hyphy_lib_path}} {{ params.pre_msa }} --input {{ params.input}}.{{ params.gene }}.tmp --reference {{ params.working_dir }}/{{ params.regions[params["gene"]]["reference"] }} --trim-from {{ params.regions[params.gene]["trim_from"] }} --trim-to {{ params.regions[params.gene]["trim_to"] }} --E 0.01 --N-fraction {{ params.regions[params["gene"]]["fraction"] }} --remove-stop-codons Yes
"""

# t1, t2 and t3 are examples of tasks created by instantiating operators
mk_tmp = BashOperator(
    task_id='mk_tmp',
    bash_command='cp {{ params.input }} {{ params.input }}.{{ params.gene }}.tmp',
    dag=dag,
)

pre_msa = BashOperator(
    task_id='pre_msa',
    bash_command=PREMSA,
    params={'regions': regions},
    dag=dag,
)

dag.doc_md = __doc__

pre_msa.doc_md = """\
#### Task Documentation
PREMSA
"""

mk_tmp >> pre_msa
