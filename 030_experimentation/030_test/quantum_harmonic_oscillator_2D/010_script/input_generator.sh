#!/bin/bash

computer="notebook"
if [ ${computer} == notebook ]; then
    kernel_dir="mendez"
elif [ ${computer} == pcfamaf ]; then
    kernel_dir="martin"
elif [ ${computer} == ccad ]; then
    kernel_dir="martinmendez"
fi

path_repository="/home/${kernel_dir}/github_repositories/my_repositories/FEMTISE_TUTORIAL/"
path_input_file="${path_repository}030_experimentation/030_test/quantum_harmonic_oscillator_2D/020_input/"
path_output_file="${path_repository}030_experimentation/030_test/quantum_harmonic_oscillator_2D/040_output/"

old_input_file_name=${path_input_file}"model_input"
new_input_file_name=${path_input_file}"input"

old_string="value"
new_string=${path_output_file}
row=1
sed -e ""${row}" s|"${old_string}"|"${new_string}"|g" ${old_input_file_name}.dat > 001_aux.dat

old_string="value"
new_string="/home/${kernel_dir}/github_repositories/my_repositories/FEMTISE.jl/src/functions/default_running_miscellaneous_functions"
row=6
sed -e ""${row}" s|"${old_string}"|"${new_string}"|g" 001_aux.dat > ${new_input_file_name}.dat

rm 001_aux.dat