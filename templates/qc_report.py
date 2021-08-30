#!/usr/bin/env python3

import datetime
import re
import pandas as pd
import skbio
import json
import os
import csv

qc_report_path = "$report_dir"

def load_new_report():
    f = open("$report",)
    data = json.load(f)

    report = [int(data['summary']['before_filtering']['total_reads']), 
            int(data['summary']['after_filtering']['total_reads']), 
            int(data['filtering_result']['low_quality_reads']), 
            int(data['filtering_result']['too_short_reads'])
            #int(data['filtering_result']['too_long_reads'])
            ]
    f.close()

    return report

#qc_report_path = "$projectDir/viz_webapp/data/$barcode/qc_report.csv"
qc_report_path = "$report_dir"

if(os.path.exists(qc_report_path)):
    #Open the last report
    last_report = pd.read_csv(qc_report_path).values.tolist()[0]
    #Open the new qc data
    new_report = load_new_report()
    #Merge data
    data = [int(x) + int(y) for x, y in zip(last_report, new_report)]
    data_dict = {'Reads before QC': [data[0]], 'Reads after QC': [data[1]], 'Low quality reads': [data[2]], 'Filtered short reads': [data[3]]}
    qc_df = pd.DataFrame.from_dict(data_dict)
    qc_df.to_csv("qc_report.csv", index=False)

else:
    new_report = load_new_report()

    data_dict = {'Reads before QC': [new_report[0]], 'Reads after QC': [new_report[1]], 'Low quality reads': [new_report[2]], 'Filtered short reads': [new_report[3]]}
    qc_df = pd.DataFrame.from_dict(data_dict)
    qc_df.to_csv("qc_report.csv", index=False)

    



