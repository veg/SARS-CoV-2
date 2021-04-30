import yaml
import datetime
import csv

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
from get_stats import get_all_unique_haplos, get_clade_counts
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
        'python': "/data/shares/veg/SARS-CoV-2/SARS-CoV-2-devel/env/bin/python3",
        'hyphy': "/data/shares/veg/SARS-CoV-2/hyphy/hyphy",
        'hyphy_mpi': "/data/shares/veg/SARS-CoV-2/hyphy/HYPHYMPI",
        'hyphy_lib_path': "/data/shares/veg/SARS-CoV-2/hyphy/res",
        'mafft': "/usr/local/bin/mafft",
        'post_msa' : "/data/shares/veg/SARS-CoV-2/hyphy-analyses-devel/codon-msa/post-msa.bf",
        'compressor' : "/data/shares/veg/SARS-CoV-2/SARS-CoV-2-devel/scripts/compressor.bf",
        'compressor2' : "/data/shares/veg/SARS-CoV-2/SARS-CoV-2-devel/scripts/compressor-2.bf",
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

dag = DAG(
    'get_stats',
    default_args=default_args,
    description='get stats',
    schedule_interval='@daily',
    start_date=datetime.datetime(2021, 4, 30),
    tags=['selection'],
)

def write_unique_haplo_counts(genes, output_fn):
    counts = get_all_unique_haplos(genes)
    try:
        with open(output_fn, "wt") as csvfile:
            spamwriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            for x in counts:
                print(x)
                spamwriter.writerow(x)
    except Exception as e:
        log.error(e)
        raise AirflowException(e)

def write_clade_counts(output_fn):
    counts = get_clade_counts()
    try:
        with open(output_fn, "wt") as csvfile:
            spamwriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            for x in counts:
                print(x)
                spamwriter.writerow(x)
    except Exception as e:
        log.error(e)
        raise AirflowException(e)


with open(dag.params["region_cfg"], 'r') as stream:
    regions = yaml.safe_load(stream)

# For each region
genes = regions.keys()
output_fn = WORKING_DIR + "/data/exports/all_haplo_counts.csv"

get_unique_haplos_task = PythonOperator(
        task_id=f'get_unique_haplos',
        python_callable=write_unique_haplo_counts,
        op_kwargs={ "genes" : genes, "output_fn" : output_fn },
        dag=dag,
    )

clade_output_fn = WORKING_DIR + "/data/exports/clade_counts.csv"
get_clade_counts_task = PythonOperator(
        task_id=f'get_clade_counts',
        python_callable=write_clade_counts,
        op_kwargs={ "output_fn" : clade_output_fn},
        dag=dag,
    )


dag.doc_md = __doc__

get_all_unique_haplos.doc_md = """\
#### Task Documentation
Creates a directory and exports selected sequences
"""

# Add export meta and export sequence tasks to be executed in parallel
[get_unique_haplos_task, get_clade_counts_task]


