import yaml
import datetime

# The DAG object; we'll need this to instantiate a DAG
from airflow import DAG

# Operators; we need this to operate!
from airflow.operators.bash import BashOperator
from airflow.models.baseoperator import cross_downstream
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

from export_sequences import export_sequences, export_premsa_sequences, export_bealign_sequences
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
    'concurrency' : 50
}


with DAG(
    'export_by_gene',
    default_args=default_args,
    description='export_by_gene',
    schedule_interval='@weekly',
    start_date=datetime.datetime(2021, 6, 2),
    tags=['selection'],
    ) as dag:

    with open(dag.params["region_cfg"], 'r') as stream:
        regions = yaml.safe_load(stream)

    last_exec_date = dag.get_latest_execution_date()

    if last_exec_date is None:
        last_exec_date = datetime.datetime(year=1970, month=1, day=1)

    unique_id = str(round(last_exec_date.timestamp()))

    directory_output = WORKING_DIR + "/data/exports/" + unique_id + "/"
    default_args['meta-output'] = directory_output + '/master-no-fasta.json'

    mk_dir_task = BashOperator(
        task_id='make_directory',
        bash_command='mkdir -p {{params.directory_output}}',
        params={"directory_output": directory_output},
        dag=dag,
    )

    export_meta_task = PythonOperator(
            task_id='export_meta',
            python_callable=export_meta,
            op_kwargs={ "config" : default_args['params'] },
            pool='mongo',
            dag=dag,
        )

    export_meta_task.set_upstream(mk_dir_task)


    export_by_gene = []

    for gene in regions.keys():

        gene_directory_output = directory_output + "/" + gene

        default_args["params"]["directory"] = gene_directory_output
        default_args["params"]["sequence-output"] = gene_directory_output + '/sequences'
        default_args["params"]["meta-output"] = gene_directory_output  + '/master-no-sequences.json'
        default_args["params"]["duplicate-output"] = gene_directory_output  + '/duplicates.json'

        nuc_sequence_output = default_args["params"]["directory"] + '/sequences_nuc.fas'
        bealign_nuc_sequence_output = default_args["params"]["directory"] + '/sequences_nuc.bealign.fas'
        prot_sequence_output = default_args["params"]["directory"] + '/sequences_protein.fas'

        default_args["params"]["nuc-sequence-output"] = nuc_sequence_output
        default_args["params"]["prot-sequence-output"] = prot_sequence_output

        gene_mk_dir_task = BashOperator(
            task_id=f'make_directory_{gene}',
            bash_command='mkdir -p {{params.directory_output}}',
            params={"directory_output": gene_directory_output},
            dag=dag,
        )


        export_duplicates_task = PythonOperator(
                task_id=f'export_duplicates_{gene}',
                python_callable=export_duplicates,
                op_kwargs={ 'output_fn' : default_args['params']['duplicate-output'], 'gene': default_args['params']['gene'] },
                pool='mongo',
                dag=dag,
            )

        export_duplicates_task.set_upstream(gene_mk_dir_task)

        # For each region
        export_by_gene_task = PythonOperator(
                task_id=f'export_premsa_sequences_{gene}',
                python_callable=export_premsa_sequences,
                op_kwargs={ "config" : default_args['params'], 'nuc_output_fn':  nuc_sequence_output, 'prot_output_fn' : prot_sequence_output, 'gene' : gene },
                pool='mongo',
                dag=dag,
            )

        export_by_gene_task.set_upstream(gene_mk_dir_task)

        export_bealign_task = PythonOperator(
            task_id=f'export_bealign_{gene}',
            python_callable=export_bealign_sequences,
            op_kwargs={ "config" : default_args['params'], 'nuc_output_fn':  bealign_nuc_sequence_output, 'gene' : gene },
            provide_context=True,
            pool='mongo',
            dag=dag,
        )

        export_bealign_task.set_upstream(gene_mk_dir_task)

        export_by_gene.extend([export_duplicates_task, export_by_gene_task, export_bealign_sequences])

    dag.doc_md = __doc__

    export_by_gene_task.doc_md = """\
        # Task Documentation
        Exports by specific gene
    """
