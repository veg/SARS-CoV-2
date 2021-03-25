import yaml
import datetime

# The DAG object; we'll need this to instantiate a DAG
from airflow import DAG

# Operators; we need this to operate!
from airflow.operators.bash import BashOperator
from airflow.operators.python import PythonOperator
from airflow.utils.dates import days_ago
from airflow.hooks.base import BaseHook
from airflow.models import Variable
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
        'gene' : 'S',
        'date' : datetime.date.today().strftime('%Y-%m-%d')
    },
    'retries': 1,
    'retry_delay': datetime.timedelta(minutes=5),
    'concurrency' : 50,
    'on_failure_callback': task_fail_slack_alert,
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

gene = default_args["params"]["gene"]
directory_output = WORKING_DIR + "/data/" + gene + "/" + default_args["params"]["date"]

default_args["params"]["directory"] = directory_output
default_args["params"]["sequence-output"] = directory_output + '/sequences'
default_args["params"]["meta-output"] = directory_output  + '/master-no-sequences.json'
default_args["params"]["duplicate-output"] = directory_output  + '/duplicates.json'

nuc_sequence_output = default_args["params"]["directory"] + '/sequences_nuc.fas'
prot_sequence_output = default_args["params"]["directory"] + '/sequences_protein.fas'

default_args["params"]["nuc-sequence-output"] = nuc_sequence_output
default_args["params"]["prot-sequence-output"] = prot_sequence_output

dag = DAG(
    'export_by_gene',
    default_args=default_args,
    description='export_by_gene',
    schedule_interval='@weekly',
    start_date=datetime.datetime(2021, 2, 25),
    tags=['selection'],
)

with open(dag.params["region_cfg"], 'r') as stream:
    regions = yaml.safe_load(stream)

mk_dir_task = BashOperator(
    task_id='make_directory',
    bash_command='mkdir -p {{params.directory}}',
    dag=dag,
)

export_meta_task = PythonOperator(
        task_id='export_meta',
        python_callable=export_meta,
        op_kwargs={ "config" : default_args['params'] },
        dag=dag,
    )

export_duplicates_task = PythonOperator(
        task_id=f'export_duplicates_{default_args["params"]["gene"]}',
        python_callable=export_duplicates,
        op_kwargs={ 'output_fn' : default_args['params']['duplicate-output'], 'gene': default_args['params']['gene'] },
        dag=dag,
    )

# For each region
export_by_gene_task = PythonOperator(
        task_id=f'export_premsa_sequences_{default_args["params"]["gene"]}',
        python_callable=export_premsa_sequences,
        op_kwargs={ "config" : default_args['params'], 'nuc_output_fn':  nuc_sequence_output, 'prot_output_fn' : prot_sequence_output, 'gene' : gene },
        dag=dag,
    )

dag.doc_md = __doc__

export_by_gene_task.doc_md = """\
# Task Documentation
Exports by specific gene
"""

# Add export meta and export sequence tasks to be executed in parallel
mk_dir_task >> [export_meta_task, export_duplicates_task, export_by_gene_task]

