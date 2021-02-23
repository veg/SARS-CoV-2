import yaml
import datetime
from datetime import timedelta

# The DAG object; we'll need this to instantiate a DAG
from airflow import DAG

# Operators; we need this to operate!
from airflow.operators.bash import BashOperator
from airflow.operators.python import PythonOperator, ShortCircuitOperator
from airflow.utils.dates import days_ago
# These args will get passed on to each operator
# You can override them on a per-task basis during operator initialization

from libs.callbacks import task_fail_slack_alert, task_success_slack_alert

import os
import sys
import pathlib
from pathlib import Path

p = os.path.abspath(str(pathlib.Path(__file__).parent.absolute()) + '/../../python/')
if p not in sys.path:
    sys.path.append(p)


from export_sequences_without_premsa import export_sequences
from store_premsa import store_premsa_file
from premsa_log_parse import mark_troubled
from mark_premsa_dupes import mark_premsa_dupes

WORKING_DIR = "/data/shares/veg/SARS-CoV-2/SARS-CoV-2-devel/"
DATE_STRING = datetime.date.today().strftime('%Y-%m-%d')

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
        'zero_length_flags' : '--kill-zero-lengths Constrain ENV="_DO_TREE_REBALANCE_=1"',
        'date_string': DATE_STRING
    },
    'retries': 3,
    'retry_delay': timedelta(minutes=5),
    'on_failure_callback': task_fail_slack_alert,
    'concurrency': 50,
    'dag_concurrency' : 50,
	'max_active_runs': 1,
    # 'execution_timeout': timedelta(minutes=30),
    # 'queue': 'bash_queue',
    # 'pool': 'backfill',
    # 'priority_weight': 10,
    # 'end_date': datetime(2016, 1, 1),
    # 'wait_for_downstream': False,
    # 'dag': dag,
    # 'sla': timedelta(hours=2),
    # 'on_retry_callback': another_function,
    # 'sla_miss_callback': yet_another_function,
    # 'trigger_rule': 'all_success'
}

dag = DAG(
    'populate_pre_msa',
    default_args=default_args,
    description='performs selection analysis',
    schedule_interval='0 8 * * *',
    start_date=days_ago(2),
    tags=['selection'],
)

with open(dag.params["region_cfg"], 'r') as stream:
    regions = yaml.safe_load(stream)

PREMSA = """
bpsh {{ params.node }} mpirun -np {{ params.num_procs }} {{ params.hyphy_mpi }} LIBPATH={{ params.hyphy_lib_path}} {{ params.pre_msa }} --input {{ params.filepath }} --reference {{ params.working_dir }}/{{ params.regions[params["gene"]]["reference"] }} --trim-from {{ params.regions[params.gene]["trim_from"] }} --trim-to {{ params.regions[params.gene]["trim_to"] }} --E 0.01 --N-fraction {{ params.regions[params["gene"]]["fraction"] }} --remove-stop-codons Yes > {{ params.stdout }}
"""

def is_export_populated(filepath):
	return Path(filepath).stat() > 0

pre_msa_tasks = []
i = 0

for gene in regions.keys():

    filepath = WORKING_DIR + 'data/premsa-processor/' + gene + '/sequences.' + default_args['params']['date_string'] + '.fasta'
    stdout = WORKING_DIR + 'data/premsa-processor/' + gene + '/sequences.' + default_args['params']['date_string'] + '.stdout.log'
    nuc_input_filepath = filepath + '_nuc.fas'
    prot_input_filepath = filepath + '_protein.fas'
    dupe_input_filepath = filepath + '_copies.json'

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

    populated_check_task = ShortCircuitOperator(
        task_id=f'check_if_populated_{gene}',
        python_callable=is_export_populated,
        params={ 'filepath': filepath },
        dag=dag
    )

    # Store nuc_input, prot_input, type
    import_premsa_seqs = PythonOperator(
        task_id=f'store_premsa_{gene}',
        python_callable=store_premsa_file,
        op_kwargs={ "nuc_input" : nuc_input_filepath, "prot_input" : prot_input_filepath, "gene": gene },
        dag=dag,
    )

    mark_troubled_task = PythonOperator(
        task_id=f'mark_troubled_{gene}',
        python_callable=mark_troubled,
        op_kwargs={ "log_file" : stdout, "gene": gene },
        dag=dag,
    )

    mark_premsa_dupes_task = PythonOperator(
        task_id=f'mark_premsa_duplicates_{gene}',
        python_callable=mark_premsa_dupes,
        op_kwargs={ "dupe_input" : dupe_input_filepath, "gene": gene },
        dag=dag,
    )

    i += 1
    pre_msa_tasks.append(export_missing >> populated_check_task >> pre_msa >> [import_premsa_seqs, mark_troubled_task] >> mark_premsa_dupes_task)

dag.doc_md = __doc__
pre_msa_tasks
