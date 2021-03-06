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

default_args["params"]["meta-output"] = WORKING_DIR + "/data/fasta/" + default_args["params"]["date"]  + '/master-no-sequences.json'
default_args["params"]["sequence-output"] = WORKING_DIR + "/data/fasta/" + default_args["params"]["date"] + '/sequences'


dag = DAG(
    'initiate_run',
    default_args=default_args,
    description='initiates run',
    schedule_interval='0 11 * * *',
    start_date=datetime.datetime(2021, 2, 10),
    tags=['selection'],
)

with open(dag.params["region_cfg"], 'r') as stream:
    regions = yaml.safe_load(stream)

mk_dir = BashOperator(
    task_id='make_directory',
    bash_command='mkdir -p {{params.working_dir}}/data/fasta/{{params.date}}',
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

# For each region
export_by_gene = []

for gene in regions.keys():

    nuc_sequence_output = WORKING_DIR + "/data/fasta/" + default_args["params"]["date"] + '/sequences.' + gene + '_nuc.fas'
    prot_sequence_output = WORKING_DIR + "/data/fasta/" + default_args["params"]["date"] + '/sequences.' + gene + '_protein.fas'
    default_args["params"]["nuc-sequence-output"] = nuc_sequence_output
    default_args["params"]["prot-sequence-output"] = prot_sequence_output

    export_premsa_sequence = PythonOperator(
            task_id=f'export_premsa_sequences_{gene}',
            python_callable=export_premsa_sequences,
            op_kwargs={ "config" : default_args['params'], 'nuc_output_fn':  nuc_sequence_output, 'prot_output_fn' : prot_sequence_output, 'gene' : gene },
            dag=dag,
        )

    export_by_gene.append(export_premsa_sequence)


dag.doc_md = __doc__

export_sequences.doc_md = """\
#### Task Documentation
Creates a directory and exports selected sequences
"""

# Add export meta and export sequence tasks to be executed in parallel
export_by_gene.extend([export_meta, export_sequences])

mk_dir >> export_by_gene

