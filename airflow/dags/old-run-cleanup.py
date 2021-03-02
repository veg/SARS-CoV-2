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
from airflow.models import Variable

p = os.path.abspath(str(pathlib.Path(__file__).parent.absolute()) + '/../../python/')
if p not in sys.path:
    sys.path.append(p)

from export_sequences_without_premsa import export_sequences
from store_premsa import store_premsa_file
from premsa_log_parse import mark_troubled
from mark_premsa_dupes import mark_premsa_dupes

WORKING_DIR = Variable.get("WORKING_DIR")
DATE_STRING = datetime.date.today().strftime('%Y-%m-%d')

default_args = {
    'owner': 'sweaver',
    'depends_on_past': False,
    'email': ['sweaver@temple.edu'],
    'email_on_failure': False,
    'email_on_retry': False,
    'params' : {
        'working_dir' : WORKING_DIR,
        'older_than' : 14,
        'region_cfg' : "/data/shares/veg/SARS-CoV-2/SARS-CoV-2-devel/airflow/libs/regions.yaml",
    },
    'retries': 3,
    'retry_delay': timedelta(minutes=5),
    'on_failure_callback': task_fail_slack_alert,
    'on_success_callback': task_success_slack_alert,
    'concurrency': 50,
    'dag_concurrency' : 50,
	'max_active_runs': 1
}

dag = DAG(
    'old_run_cleanup',
    default_args=default_args,
    description='performs selection analysis',
    schedule_interval='0 8 0 * *',
    start_date=days_ago(2),
    tags=['selection'],
)

with open(dag.params["region_cfg"], 'r') as stream:
    regions = yaml.safe_load(stream)

ZIP_OLD_RUNS = """
find {{ params.working_dir }}/data/fasta -maxdepth 1 -type d -mtime +{{ params.older_than }} -exec tar -czvf {}.tar.gz {} \;
"""

REMOVE_OLD_RUNS = """
find {{ params.working_dir }}/data/fasta -maxdepth 1 -type d -mtime +{{ params.older_than }} -exec rm -R {} \;
"""


zip_old_runs_task = BashOperator(
    task_id='zip_old_runs',
    bash_command=ZIP_OLD_RUNS,
    dag=dag,
)

remove_tarred_dirs_task = BashOperator(
    task_id='remove_tarred_dirs',
    bash_command=REMOVE_OLD_RUNS,
    dag=dag,
)

zip_old_runs_task >> remove_tarred_dirs_task
