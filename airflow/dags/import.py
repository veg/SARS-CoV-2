import yaml
import datetime

# The DAG object; we'll need this to instantiate a DAG
from airflow import DAG

# Operators; we need this to operate!
from airflow.operators.bash import BashOperator
from airflow.operators.python import PythonOperator, ShortCircuitOperator
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

from filter_gisaid_exports import filter_gisaid_exports_by_dir

WORKING_DIR = Variable.get("WORKING_DIR")

default_args = {
    'owner': 'sweaver',
    'depends_on_past': False,
    'email': ['sweaver@temple.edu'],
    'email_on_failure': False,
    'email_on_retry': False,
    'params' : {
        'working_dir' : WORKING_DIR,
        'import_dir' : WORKING_DIR + '/data/to-import/',
        'imported_dir': WORKING_DIR + '/data/imported/',
    },
    'retries': 5,
    'retry_delay': datetime.timedelta(minutes=5),
    'execution_timeout': datetime.timedelta(minutes=480),
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

# Airflow uses UTC time
dag = DAG(
    'import',
    default_args=default_args,
    description='performs selection analysis',
    schedule_interval='0 4 * * *',
    start_date=datetime.datetime(2021, 2, 10),
    tags=['selection'],
)

retrieve_meta_from_gisaid = BashOperator(
    task_id='retrieve_meta_from_gisaid',
    bash_command='/home/airflow/.nvm/versions/node/v13.14.0/bin/node {{ params.working_dir }}/js/get_metadata.js',
    env={**os.environ},
    dag=dag,
)

retrieve_fasta_from_gisaid = BashOperator(
    task_id='retrieve_fasta_from_gisaid',
    bash_command='/home/airflow/.nvm/versions/node/v13.14.0/bin/node {{ params.working_dir }}/js/get_seqs.js',
    env={**os.environ},
    dag=dag,
)

gunzip_files = BashOperator(
    task_id='gunzip_files',
    bash_command='gunzip {{ params.import_dir}}/*.gz',
    dag=dag,
)

new_meta = default_args['params']['import_dir'] + 'new.tsv'
new_fasta = default_args['params']['import_dir'] + 'new.fasta'

# Split out items from
split_out_new_task = PythonOperator(
    task_id='split_out_new',
    python_callable=filter_gisaid_exports_by_dir,
    op_kwargs={ "dir": default_args['params']['import_dir'], "fasta_output" : new_fasta, "meta_output" : new_meta },
    pool='mongo',
    dag=dag,
)

import_tsv = BashOperator(
    task_id='import_tsv',
    bash_command='node {{ params.working_dir }}/js/submit-tsv-to-mongo.js {{ params.meta_tsv }}',
    params={'meta_tsv': new_meta},
    dag=dag,
)

update_mongo_with_sequences = BashOperator(
    task_id='update_with_sequences',
    bash_command='python3 {{ params.working_dir }}/python/update_with_sequence_name.py -i {{ params.fasta }}',
    params={'fasta': new_fasta},
    dag=dag,
)

mv_files = BashOperator(
    task_id='move_files',
    bash_command='mv {{ params.import_dir }}/*.tsv {{ params.import_dir }}/*.fasta {{ params.imported_dir }}',
    dag=dag,
)

dag.doc_md = __doc__

# import_tsv.doc_md = """\
# #### Task Documentation
# IMPORT TSV FROM GISAID
# """

[retrieve_meta_from_gisaid, retrieve_fasta_from_gisaid] >> gunzip_files >> split_out_new_task >> import_tsv >> update_mongo_with_sequences >> mv_files
