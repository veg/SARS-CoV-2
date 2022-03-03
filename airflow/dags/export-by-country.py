import yaml
import json
import datetime

# The DAG object; we'll need this to instantiate a DAG
from airflow import DAG

# Operators; we need this to operate!
from airflow.operators.bash import BashOperator
from airflow.operators.python import PythonOperator
from airflow.utils.dates import days_ago
from airflow.hooks.base import BaseHook
from airflow.models import Variable
from libs.callbacks import task_fail_slack_alert, task_success_slack_alert, dag_fail_slack_alert, dag_success_slack_alert

# These args will get passed on to each operator
# You can override them on a per-task basis during operator initialization

from airflow.contrib.operators.slack_webhook_operator import SlackWebhookOperator

import os
import sys
import pathlib

p = os.path.abspath(str(pathlib.Path(__file__).parent.absolute()) + '/../../python/')
if p not in sys.path:
    sys.path.append(p)

from export_sequences import export_sequences
from export_meta import export_meta

WORKING_DIR = Variable.get("WORKING_DIR")

SLACK_CONN_ID = 'slack'

default_args = {
    'owner': 'sweaver',
    'depends_on_past': False,
    'email': ['sweaver@temple.edu'],
    'email_on_failure': False,
    'email_on_retry': False,
    'params' : {
        'working_dir' : WORKING_DIR,
        'region_cfg' : "/data/shares/veg/SARS-CoV-2/SARS-CoV-2-devel/airflow/libs/regions.yaml",
        'date' : datetime.date.today().strftime('%Y-%m-%d')
    },
    'retries': 1,
    'retry_delay': datetime.timedelta(minutes=5)
}

with DAG(
    'export_by_country',
    default_args=default_args,
    description='exports clades',
    schedule_interval='@weekly',
    start_date=datetime.datetime(2021, 6, 2),
    on_failure_callback=dag_fail_slack_alert,
    on_success_callback=dag_success_slack_alert,
    tags=['export'],
    ) as dag:

    with open(dag.params["region_cfg"], 'r') as stream:
        regions = yaml.safe_load(stream)

    last_exec_date = dag.get_latest_execution_date()

    if last_exec_date is None:
        last_exec_date = datetime.datetime(year=1970, month=1, day=1)

    unique_id = str(round(last_exec_date.timestamp()))
    directory_output = WORKING_DIR + "/data/exports/whole-genome-clades/" + unique_id + "/"


    mk_dir_task = BashOperator(
        task_id='make_directory',
        bash_command='mkdir -p {{params.directory_output}}',
        params={"directory_output": directory_output},
        dag=dag,
    )


    countries = [
        "Argentina",
        "Australia",
        "Austria",
        "Belgium",
        "Brazil",
        "Bulgaria",
        "Cambodia",
        "Canada",
        "Chile",
        "Colombia",
        "Croatia",
        "Czech Republic",
        "Denmark",
        "Estonia",
        "Finland",
        "France",
        "Germany",
        "Greece",
        "Iceland",
        "India",
        "Indonesia",
        "Ireland",
        "Israel",
        "Italy",
        "Japan",
        "Latvia",
        "Lithuania",
        "Luxembourg",
        "Malaysia",
        "Mexico",
        "Netherlands",
        "New Zealand",
        "Norway",
        "Peru",
        "Philippines",
        "Poland",
        "Portugal",
        "Romania",
        "Russia",
        "Singapore",
        "Slovakia",
        "South Africa",
        "South Korea",
        "Spain",
        "Sweden",
        "Switzerland",
        "Thailand",
        "Turkey",
        "USA",
        "United Kingdom"
        ]

    for country in countries:

        params = {}

        sanitized_country = country.lower().replace(' ', '_')
        params['meta-output'] = directory_output + '/' + sanitized_country + '-no-fasta.json'
        params["sequence-output"] = directory_output + '/' + sanitized_country + '.fas'
        params['only-uniques'] = False
        params["countries"] = [country]

        PRIORITY_RANK = 9999

        export_meta_task = PythonOperator(
                task_id=f'export_meta_{sanitized_country}',
                python_callable=export_meta,
                op_kwargs={ "config" : params },
                priority_weight=PRIORITY_RANK,
                dag=dag,
            )

        export_meta_task.set_upstream(mk_dir_task)

        export_sequences_task = PythonOperator(
                task_id=f'export_sequences_{sanitized_country}',
                python_callable=export_sequences,
                op_kwargs={ "config" : params },
                priority_weight=PRIORITY_RANK,
                dag=dag,
            )

        export_sequences_task.set_upstream(mk_dir_task)
