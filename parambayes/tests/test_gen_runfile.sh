#!/bin/bash
test -e ssshtest ||  wget -q http://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest
. ssshtest

run test_gen_runfile_help gen_runfile -h
assert_exit_code 0
assert_in_stdout usage:

run test_gen_runfile_normal gen_runfile
assert_exit_code 0
assert_in_stdout specifications

run test_invalid_compound gen_runfile --compound H2H2
assert_exit_code 2
assert_in_stderr H2H2

run test_trange_too_large gen_runfile --trange -300 1000
assert_exit_code 1
assert_in_stderr unphysical

run test_nsteps_negative gen_runfile --n_steps -300
assert_exit_code 1
assert_in_stderr integer

run test_n_data_points_negative gen_runfile --n_data_points -300
assert_exit_code 1
assert_in_stderr integer

run test_bad_json_file gen_runfile --priors_JSON NOTAJSON
assert_exit_code 1
assert_in_stderr json

run test_bad_json_dist gen_runfile --priors_JSON '{"epsilon":["gamma",[134.665,0,0.2683]],"sigma":["gamma",[2646.618,0,0.0001246]],"L":["gamma",[53.1725,0,0.002059]],"Q":["gaussian",[0,1]]}'
assert_exit_code 1
assert_in_stderr gaussian

run test_bad_json_dist_params  gen_runfile --priors_JSON '{"epsilon":["gamma",[134.665,0,0.2683]],"sigma":["gamma",[2646.618,0,0.0001246]],"L":["gamma",[53.1725,0,0.002059]],"Q":["exponential",[1,2,3,4,0,1]]}'
assert_exit_code 1
assert_in_stderr parameters