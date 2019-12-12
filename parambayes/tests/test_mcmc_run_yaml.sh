#!/bin/bash
test -e ssshtest ||  wget -q http://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest
. ssshtest

run test_mcmc_run_yaml_help run_mcmc_sim -h
assert_exit_code 0
assert_in_stdout usage:

gen_runfile --n_steps 1000 --date False

run test_gen_run_normal run_mcmc_sim -f runfiles/N2_All_prod_2019-12-11.yml
assert_exit_code 0
assert_in_stdout 'Simulation Done!'