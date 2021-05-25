import yaml
import datetime
from dateutil.relativedelta import *
from dateutil.rrule import *
from dateutil.parser import *

# The DAG object; we'll need this to instantiate a DAG
from airflow import DAG

# Operators; we need this to operate!
from airflow.utils.task_group import TaskGroup
from airflow.operators.bash import BashOperator
from airflow.operators.python import PythonOperator
from airflow.utils.dates import days_ago
from airflow.hooks.base import BaseHook
from airflow.models import Variable
from libs.callbacks import dag_fail_slack_alert, dag_success_slack_alert

import os
import sys
import pathlib

p = os.path.abspath(str(pathlib.Path(__file__).parent.absolute()) + '/../../python/')
if p not in sys.path:
    sys.path.append(p)

from get_raw_duplicates import write_nuc_raw_duplicates
from random_sample import export_random

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
        'python': "/data/shares/veg/SARS-CoV-2/SARS-CoV-2-devel/env/bin/python3",
        'hyphy': "/data/shares/veg/SARS-CoV-2/hyphy/hyphy",
        'hyphy_mpi': "/data/shares/veg/SARS-CoV-2/hyphy/HYPHYMPI",
        'hyphy_lib_path': "/data/shares/veg/SARS-CoV-2/hyphy/res",
        'compressor' : "/data/shares/veg/SARS-CoV-2/SARS-CoV-2-devel/scripts/compressor.bf",
        'compressor2' : "/data/shares/veg/SARS-CoV-2/SARS-CoV-2-devel/scripts/compressor-2.bf",
        'tn93_cluster' : "/usr/local/bin/tn93-cluster",
        'meta-output' : WORKING_DIR + '/master-no-fasta.json',
        'gene' : 'S',
        'date' : datetime.date.today().strftime('%Y-%m-%d')
    },
    'retries': 1,
    'retry_delay': datetime.timedelta(minutes=5),
    'concurrency' : 50,
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

directory_output = WORKING_DIR + "/data/random-sample/"
default_args["params"]["output-dir"] = directory_output

def export_random_by_windows(windows, dir, size):

    for window in windows:
        # Output directory
        output_dir = os.path.join(dir, '_'.join(window))
        # Make directory
        try:
            os.mkdir(output_dir)
        except:
            print('Directory already exists')
        # Create config with directory and window
        config = {"output-dir" : output_dir }
        config["collection-date-range"] = window
        config["size"] = 1000

        export_random(config)

with DAG(
    'random_sample',
    default_args=default_args,
    description='random_sample',
    schedule_interval=None,
    start_date=datetime.datetime(2021, 2, 25),
    on_failure_callback=dag_fail_slack_alert,
    on_success_callback=dag_success_slack_alert,
    tags=['selection'],
    ) as dag:

    with open(dag.params["region_cfg"], 'r') as stream:
        regions = yaml.safe_load(stream)

    # Write random samples for each month
    TODAY = datetime.date.today()
    LASTMONTH = TODAY-relativedelta(months=+1, day=31)
    ONEMONTHAGO = TODAY-relativedelta(months=+1, day=1)

    starts = [dt.strftime('%Y-%m-%d') for dt in rrule(MONTHLY, interval=1,bymonthday=(1),dtstart=parse("20191201T000000"), until=ONEMONTHAGO)]
    ends = [dt.strftime('%Y-%m-%d') for dt in rrule(MONTHLY, interval=1,bymonthday=(-1),dtstart=parse("20191231T000000"), until=LASTMONTH)]
    months = set(list(zip(starts,ends)))

    mk_dir_task = BashOperator(
        task_id='make_directory',
        bash_command='mkdir -p {{params["output-dir"]}}',
        dag=dag,
    )

    # export_random_by_windows(windows, dir, size)
    export_random_samples_task = PythonOperator(
            task_id='export_random_samples',
            python_callable=export_random_by_windows,
            op_kwargs={ "windows" : months, "dir" : dag.params["output-dir"], "size": 1000 },
            pool='mongo',
            dag=dag,
        )

    # For each region
    export_by_gene = []

    for month in months:

        month_str = '_'.join(month)
        OUTPUT_DIR = os.path.join(dag.params["output-dir"], month_str)

        for gene in regions.keys():

            filepath_prefix = OUTPUT_DIR + '/sequences.' + gene

            nuc_sequence_output = filepath_prefix + '.fas'
            uniques_fn = filepath_prefix + '.uniques.fas'
            duplicate_output = filepath_prefix + '.duplicates.json'

            variants_csv_output = filepath_prefix + '.variants.csv'
            variants_json_output = filepath_prefix + '.variants.json'
            filtered_fasta_output = filepath_prefix + '.compressed.filtered.fas'
            filtered_json_output = filepath_prefix + '.filtered.json'
            output_edits_fn = filepath_prefix + '.filtered.edits.json'
            tn93_output = filepath_prefix + '.tn93.csv'

            compressor_duplicate_out = filepath_prefix + '.duplicates.variants.json'

            tree_output = filepath_prefix + '.compressed.filtered.fas.rapidnj.bestTree'
            sto_output = filepath_prefix + '.compressed.filtered.sto';


            with TaskGroup(f"alignment_{month_str}_{gene}") as alignment:

                # Occasional errors when cleaning up tmp files, so or'ing true
                cleanup_task = BashOperator(
                    task_id=f'cleanup',
                    bash_command="sed -i '/^>/! s/[^ACTG-]/N/g' $NUC_OUTPUT_FN || true",
                    env={'NUC_OUTPUT_FN': nuc_sequence_output, **os.environ},
                    dag=dag
                )

                cleanup_task

            with TaskGroup(f"duplicates_{month_str}_{gene}") as duplicates_group:

                compute_duplicates_task = PythonOperator(
                    task_id=f'write_raw_duplicates',
                    python_callable=write_nuc_raw_duplicates,
                    op_kwargs={ "input" : nuc_sequence_output, "duplicate_output" : duplicate_output, 'uniques_output': uniques_fn },
                    dag=dag,
                )

                compute_duplicates_task

            # $HYPHY LIBPATH=$HYPHYLIBPATH $COMPRESSOR --msa ${FILE}.${GENE}.compressed.fas --regexp "epi_isl_([0-9]+)" --duplicates ${FILE}.${GENE}.duplicates.json --output ${FILE}.${GENE}.variants.csv --json ${FILE}.${GENE}.variants.json --duplicate-out ${FILE}.${GENE}.duplicates.variants.json

            with TaskGroup(f"filter_{month_str}_{gene}") as filter:
                COMPRESSOR = """
                {{ params.hyphy }} LIBPATH={{params.hyphy_lib_path}} {{ params.compressor }} --msa $FASTA_FN --regexp "epi_isl_([0-9]+)" --duplicates $DUPLICATE_FN --output $VARIANTS_CSV_FN  --json $VARIANTS_JSON_FN --duplicate-out $COMPRESSOR_DUPLICATE_OUT
                """
                compressor_task = BashOperator(
                    task_id=f'compressor',
                    bash_command=COMPRESSOR,
                    env={'FASTA_FN': uniques_fn, 'DUPLICATE_FN': duplicate_output, 'VARIANTS_CSV_FN': variants_csv_output, 'VARIANTS_JSON_FN': variants_json_output, 'COMPRESSOR_DUPLICATE_OUT': compressor_duplicate_out, **os.environ},
                    dag=dag
                )

                # --output-edits ${FILE}.${GENE}.filtered.edits.json
                COMPRESSOR2 = """
                {{ params.hyphy }} LIBPATH={{params.hyphy_lib_path}} {{ params.compressor2 }} --msa $FASTA_FN --duplicates $DUPLICATE_FN --csv $VARIANTS_CSV_FN  --byseq $VARIANTS_JSON_FN --p 0.95 --output $FILTERED_FASTA_FN --json $FILTERED_JSON_FN --output-edits ${OUTPUT_EDITS}
                """
                compressor_two_task = BashOperator(
                    task_id=f'compressor_two',
                    bash_command=COMPRESSOR2,
                    env={'FASTA_FN': uniques_fn, 'DUPLICATE_FN': compressor_duplicate_out, 'VARIANTS_CSV_FN': variants_csv_output, 'VARIANTS_JSON_FN': variants_json_output, 'FILTERED_FASTA_FN': filtered_fasta_output, 'FILTERED_JSON_FN': filtered_json_output, 'OUTPUT_EDITS': output_edits_fn, **os.environ},
                    dag=dag
                )

                compressor_task >> compressor_two_task

            tn93_task = BashOperator(
                task_id=f'tn93_{month_str}_{gene}',
                bash_command='tn93 -o {{params.tn93_output}} {{params.filtered_fasta_fn}}',
                params={'filtered_fasta_fn': filtered_fasta_output, 'tn93_output': tn93_output, **os.environ},
                dag=dag
            )

            export_by_gene.append(alignment >> duplicates_group >> filter >> tn93_task)

    dag.doc_md = __doc__

    # Add export meta and export sequence tasks to be executed in parallel
    mk_dir_task >> export_random_samples_task

