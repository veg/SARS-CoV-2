import yaml
from datetime import timedelta

# The DAG object; we'll need this to instantiate a DAG
from airflow import DAG

# Operators; we need this to operate!
from airflow.operators.bash import BashOperator
from airflow.operators.python import PythonOperator
from airflow.utils.dates import days_ago
from airflow.hooks.base import BaseHook
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

WORKING_DIR = "/data/shares/veg/SARS-CoV-2/SARS-CoV-2-devel/"

SLACK_CONN_ID = 'slack'

def task_fail_slack_alert(context):
    slack_webhook_token = BaseHook.get_connection(SLACK_CONN_ID).password
    print(slack_webhook_token)
    slack_msg = """
            :red_circle: Task Failed.
            *Task*: {task}
            *Dag*: {dag}
            *Execution Time*: {exec_date}
            *Log Url*: {log_url}
            """.format(
            task=context.get('task_instance').task_id,
            dag=context.get('task_instance').dag_id,
            ti=context.get('task_instance'),
            exec_date=context.get('execution_date'),
            log_url=context.get('task_instance').log_url,
        )
    failed_alert = SlackWebhookOperator(
        task_id='slack_test',
        http_conn_id='slack',
        webhook_token=slack_webhook_token,
        message=slack_msg,
        username='airflow')
    return failed_alert.execute(context=context)

def task_success_slack_alert(context):
    slack_webhook_token = BaseHook.get_connection(SLACK_CONN_ID).password
    slack_msg = """
            :large_green_circle: Task Succeeded.
            *Task*: {task}
            *Dag*: {dag}
            *Execution Time*: {exec_date}
            *Log Url*: {log_url}
            """.format(
            task=context.get('task_instance').task_id,
            dag=context.get('task_instance').dag_id,
            ti=context.get('task_instance'),
            exec_date=context.get('execution_date'),
            log_url=context.get('task_instance').log_url,
        )
    failed_alert = SlackWebhookOperator(
        task_id='slack_test',
        http_conn_id='slack',
        webhook_token=slack_webhook_token,
        message=slack_msg,
        username='airflow')
    return failed_alert.execute(context=context)


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
        'input': '/data/shares/veg/SARS-CoV-2/SARS-CoV-2-devel/data/to-import/metadata.tsv',
        'fasta_input': '/data/shares/veg/SARS-CoV-2/SARS-CoV-2-devel/data/to-import/sequences.fasta',
        'imported_dir': '/data/shares/veg/SARS-CoV-2/SARS-CoV-2-devel/data/imported/',
        'zero_length_flags' : '--kill-zero-lengths Constrain ENV="_DO_TREE_REBALANCE_=1"',
        'meta-output' : WORKING_DIR + '/master-no-sequences.json',
        'sequence-output' : WORKING_DIR + '/sequences.fasta',
        'clades' : ["B.1.351", "P.1"],
        'clade-type' : 'pangolinLineage'
    },
    'retries': 1,
    'retry_delay': timedelta(minutes=5),
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
    # 'on_failure_callback': some_function,
    # 'on_success_callback': some_other_function,
    # 'on_retry_callback': another_function,
    # 'sla_miss_callback': yet_another_function,
    # 'trigger_rule': 'all_success'
}

dag = DAG(
    'export_clades',
    default_args=default_args,
    description='exports clades',
    schedule_interval=timedelta(days=1),
    start_date=days_ago(2),
    tags=['selection'],
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

export_sequences.doc_md = """\
#### Task Documentation
IMPORT TSV FROM GISAID
"""

[export_meta, export_sequences]
