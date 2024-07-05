#!/bin/bash

computer=$1
if [ ${computer} == notebook ]; then
    kernel_dir="mendez"
elif [ ${computer} == pcfamaf ]; then
    kernel_dir="martin"
elif [ ${computer} == ccad ]; then
    kernel_dir="martinmendez"
fi

path_repository="/home/${kernel_dir}/github_repositories/my_repositories/FEMTISE_TUTORIAL/"
path_input_file="${path_repository}030_experimentation/030_test/coulomb_potential_2D_parameter_variation/020_input/"
path_src_file="${path_repository}030_experimentation/030_test/coulomb_potential_2D_parameter_variation/030_src/"
input_file_name="${path_input_file}input"
model_run_file_name="${path_input_file}run"
run_code_file_name="${path_src_file}run"

develop_package="false"
sed -e "2 s|value|"${develop_package}"|g" -e "3 s|value|"${computer}"|g" ${model_run_file_name}.dat > aux.dat

sed "31 s|value|"${input_file_name}"|g" aux.dat > ${run_code_file_name}.jl

rm -f aux.dat