import yaml
import datetime

# The DAG object; we'll need this to instantiate a DAG
from airflow import DAG

# Operators; we need this to operate!
from airflow.operators.bash import BashOperator
from airflow.utils.dates import days_ago
from airflow.hooks.base import BaseHook
# These args will get passed on to each operator
# You can override them on a per-task basis during operator initialization

from airflow.contrib.operators.slack_webhook_operator import SlackWebhookOperator

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
        'import_dir' : '/data/shares/veg/SARS-CoV-2/SARS-CoV-2-devel/data/to-import/',
        'imported_dir': '/data/shares/veg/SARS-CoV-2/SARS-CoV-2-devel/data/imported/',
    },
    'retries': 5,
    'retry_delay': datetime.timedelta(minutes=5),
    'execution_timeout': datetime.timedelta(minutes=30),
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
    # 'on_failure_callback': some_function,
    # 'on_success_callback': some_other_function,
    # 'on_retry_callback': another_function,
    # 'sla_miss_callback': yet_another_function,
    # 'trigger_rule': 'all_success'
}

dag = DAG(
    'import',
    default_args=default_args,
    description='performs selection analysis',
    schedule_interval='0 2 * * *',
    start_date=datetime.datetime(2021, 2, 10),
    tags=['selection'],
)

retrieve_meta_from_gisaid = BashOperator(
    task_id='retrieve_meta_from_gisaid',
    bash_command='node {{ params.working_dir }}/js/get_metadata.js',
    dag=dag,
)

retrieve_fasta_from_gisaid = BashOperator(
    task_id='retrieve_fasta_from_gisaid',
    bash_command='node {{ params.working_dir }}/js/get_seqs.js',
    dag=dag,
)

gunzip_files = BashOperator(
    task_id='gunzip_files',
    bash_command='gunzip {{ params.import_dir}}/*.gz',
    dag=dag,
)

import_tsv = BashOperator(
    task_id='import_tsv',
    bash_command='for x in $(ls {{ params.import_dir }}/*.tsv); do node {{ params.working_dir }}/js/submit-tsv-to-mongo.js $x; done;',
    dag=dag,
)

update_mongo_with_sequences = BashOperator(
    task_id='update_with_sequences',
    bash_command='for x in $(ls {{ params.import_dir }}/*.fasta); do python3 {{ params.working_dir }}/python/update-with-sequence-name.py -i $x; done;',
    dag=dag,
)

mv_files = BashOperator(
    task_id='move_files',
    bash_command='mv {{ params.import_dir }}/*.tsv {{ params.import_dir }}/*.fasta {{ params.imported_dir }}',
    dag=dag,
)

dag.doc_md = __doc__

import_tsv.doc_md = """\
#### Task Documentation
IMPORT TSV FROM GISAID
"""

[retrieve_meta_from_gisaid, retrieve_fasta_from_gisaid] >> gunzip_files >> import_tsv >> update_mongo_with_sequences >> mv_files
