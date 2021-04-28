import yaml
import datetime

# The DAG object; we'll need this to instantiate a DAG
from airflow import DAG

# Operators; we need this to operate!
from airflow.models.baseoperator import cross_downstream, chain
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
from remove_seq import remove_reference_seq, reserve_only_original_input
from merge_duplicates import merge_duplicates
from fix_duplicates import fix_duplicates
from update_fasta_duplicates import update_fasta_duplicates

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
        'mafft': "/usr/local/bin/mafft",
        'date' : datetime.date.today().strftime('%Y-%m-%d')
    },
    'retries': 1,
    'retry_delay': datetime.timedelta(minutes=5),
    'task_concurrency' : 5,
    'on_failure_callback': task_fail_slack_alert,
    'on_success_callback': task_success_slack_alert,
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

export_directory = WORKING_DIR + "/data/exports/" + default_args["params"]["date"]

default_args["params"]["meta-output"] = export_directory + '/master-no-sequences.json'
default_args["params"]["sequence-output"] = export_directory + '/sequences'

dag = DAG(
    'export_directory',
    default_args=default_args,
    description='exports mongodb to FASTA to and JSON metadata',
    schedule_interval=None,
    start_date=datetime.datetime(2021, 4, 14),
    tags=['selection'],
)

with open(dag.params["region_cfg"], 'r') as stream:
    regions = yaml.safe_load(stream)

mk_dir = BashOperator(
    task_id='make_directory',
    bash_command='mkdir -p {{params.export_directory}}',
    params={ 'export_directory': export_directory },
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

    reference_filepath = WORKING_DIR + 'reference_genes/reference.' + gene + '_protein.fas'
    filepath_prefix = export_directory + '/sequences.' + gene

    nuc_sequence_output = filepath_prefix + '_nuc.fas'
    prot_sequence_output = filepath_prefix + '_protein.fas'

    initial_duplicate_output = filepath_prefix + '.initial.duplicates.json'
    protein_duplicate_output = filepath_prefix + '.protein.duplicates.json'
    duplicate_output = filepath_prefix + '.duplicates.json'

    default_args["params"]["nuc-sequence-output"] = nuc_sequence_output
    default_args["params"]["prot-sequence-output"] = prot_sequence_output
    default_args["params"]["duplicate-output"] = duplicate_output
    default_args["params"]["protein-duplicate-output"] = duplicate_output
    default_args["params"]["inital-duplicate-output"] = initial_duplicate_output

    export_premsa_sequence_task = PythonOperator(
            task_id=f'export_premsa_sequences_{gene}',
            python_callable=export_premsa_sequences,
            op_kwargs={ "config" : default_args['params'], 'nuc_output_fn':  nuc_sequence_output, 'prot_output_fn' : prot_sequence_output, 'gene' : gene },
            dag=dag,
        )

    export_duplicates_task = PythonOperator(
        task_id=f'export_duplicates_{gene}',
        python_callable=export_duplicates,
        op_kwargs={ 'output_fn' : initial_duplicate_output, 'gene': gene },
        dag=dag,
    )

    export_by_gene.append(export_premsa_sequence_task)
    export_by_gene.append(export_duplicates_task)

dag.doc_md = __doc__

export_sequences.doc_md = """\
#### Task Documentation
Creates a directory and exports selected sequences
"""

# Add export meta and export sequence tasks to be executed in parallel
export_by_gene.extend([export_meta, export_sequences])
mk_dir >> export_by_gene


