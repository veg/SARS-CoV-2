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
        'date' : datetime.date.today().strftime('%Y-%m-%d')
    },
    'retries': 1,
    'retry_delay': datetime.timedelta(minutes=5),
    'task_concurrency' : 5,
    'on_failure_callback': task_fail_slack_alert,
    'on_success_callback': task_success_slack_alert
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

OUTPUT_DIR = WORKING_DIR + "/data/fasta/" + default_args["params"]["date"]  + '-last100k-by-collection-date'
default_args["params"]["output-dir"] = OUTPUT_DIR
default_args["params"]["meta-output"] = OUTPUT_DIR + '/master-no-sequences.json'
default_args["params"]["sequence-output"] = OUTPUT_DIR + '/sequences'

dag = DAG(
    'initiate_run_last_100k_by_collection_date',
    default_args=default_args,
    description='initiates run with sequences from latest 100k sequences',
    schedule_interval='@weekly',
    start_date=datetime.datetime(2021, 3, 25),
    tags=['selection'],
)

with open(dag.params["region_cfg"], 'r') as stream:
    regions = yaml.safe_load(stream)

mk_dir = BashOperator(
    task_id='make_directory',
    bash_command='mkdir -p {{params.output}}',
    params={'output': default_args['params']['output-dir']},
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


dag.doc_md = __doc__

# Add export meta and export sequence tasks to be executed in parallel
mk_dir >> [export_meta, export_sequences]

