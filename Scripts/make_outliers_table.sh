#!/bin/bash
updown="up"
location_of_py_file="outlier_analysis/make_outliers_table.py"
location_of_data_file="phospho/WJ_Clean_phospho_iP.tsv"
iqrs_over_median=0.5 #Note 1.5 IQRs is suggested, this is just for test data.
gene_column_name="Gene"
output_prefix="phospho/outliers_output_${updown}"
sample_names_file="phospho/Samplelist.txt"
aggregate=True
write_frac_table="True"

python3 ${location_of_py_file} \
--input_df ${location_of_data_file} \
--iqrs_over_median ${iqrs_over_median} \
--gene_column_name ${gene_column_name} \
--output_prefix ${output_prefix} \
--sample_names_file ${sample_names_file} \
--aggregate ${aggregate} \
--up_or_down ${updown} \
--write_frac_table ${write_frac_table}
