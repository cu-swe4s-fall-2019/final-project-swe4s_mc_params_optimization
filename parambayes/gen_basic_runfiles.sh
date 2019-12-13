#Generate runfiles for 'rhol+Psat' cases


gen_runfile \
--compound O2\
 --prop rhol+Psat\
 --trange 0.55 0.95\
 --n_steps 2000000\
 --n_data_points 10\
 --priors_JSON '{"epsilon":["exponential",[0,100]],"sigma":["exponential",[0,5]],"L":["exponential",[0,3]],"Q":["exponential",[0,1]]}'\
 --save_traj 'True'\
 --date 'False'

gen_runfile \
--compound N2\
 --prop rhol+Psat\
 --trange 0.55 0.95\
 --n_steps 2000000\
 --n_data_points 10\
 --priors_JSON '{"epsilon":["exponential",[0,100]],"sigma":["exponential",[0,5]],"L":["exponential",[0,3]],"Q":["exponential",[0,1]]}'\
 --save_traj 'True'\
 --date 'False'

gen_runfile \
--compound Br2\
 --prop rhol+Psat\
 --trange 0.55 0.95\
 --n_steps 2000000\
 --n_data_points 10\
 --priors_JSON '{"epsilon":["exponential",[0,100]],"sigma":["exponential",[0,5]],"L":["exponential",[0,3]],"Q":["exponential",[0,1]]}'\
 --save_traj 'True'\
 --date 'False'

gen_runfile \
--compound F2\
 --prop rhol+Psat\
 --trange 0.50 0.6\
 --n_steps 2000000\
 --n_data_points 10\
 --priors_JSON '{"epsilon":["exponential",[0,100]],"sigma":["exponential",[0,5]],"L":["exponential",[0,3]],"Q":["exponential",[0,1]]}'\
 --save_traj 'True'\
 --date 'False'

gen_runfile \
--compound C2H2\
 --prop rhol+Psat\
 --trange 0.55 0.95\
 --n_steps 2000000\
 --n_data_points 10\
 --priors_JSON '{"epsilon":["exponential",[0,100]],"sigma":["exponential",[0,5]],"L":["exponential",[0,3]],"Q":["exponential",[0,1]]}'\
 --save_traj 'True'\
 --date 'False'

gen_runfile \
--compound C2H4\
 --prop rhol+Psat\
 --trange 0.55 0.95\
 --n_steps 2000000\
 --n_data_points 10\
 --priors_JSON '{"epsilon":["exponential",[0,100]],"sigma":["exponential",[0,5]],"L":["exponential",[0,3]],"Q":["exponential",[0,1]]}'\
 --save_traj 'True'\
 --date 'False'

gen_runfile \
--compound C2H6\
 --prop rhol+Psat\
 --trange 0.55 0.95\
 --n_steps 2000000\
 --n_data_points 10\
 --priors_JSON '{"epsilon":["exponential",[0,100]],"sigma":["exponential",[0,5]],"L":["exponential",[0,3]],"Q":["exponential",[0,1]]}'\
 --save_traj 'True'\
 --date 'False'

gen_runfile \
--compound C2F4\
 --prop rhol+Psat\
 --trange 0.55 0.95\
 --n_steps 2000000\
 --n_data_points 10\
 --priors_JSON '{"epsilon":["exponential",[0,100]],"sigma":["exponential",[0,5]],"L":["exponential",[0,3]],"Q":["exponential",[0,1]]}'\
 --save_traj 'True'\
 --date 'False'

#Generate 'ALL' runfiles

gen_runfile \
--compound O2\
 --prop All\
 --trange 0.55 0.95\
 --n_steps 2000000\
 --n_data_points 10\
 --priors_JSON '{"epsilon":["exponential",[0,100]],"sigma":["exponential",[0,5]],"L":["exponential",[0,3]],"Q":["exponential",[0,1]]}'\
 --save_traj 'True'\
 --date 'False'

gen_runfile \
--compound N2\
 --prop All\
 --trange 0.55 0.95\
 --n_steps 2000000\
 --n_data_points 10\
 --priors_JSON '{"epsilon":["exponential",[0,100]],"sigma":["exponential",[0,5]],"L":["exponential",[0,3]],"Q":["exponential",[0,1]]}'\
 --save_traj 'True'\
 --date 'False'

gen_runfile \
--compound Br2\
 --prop All\
 --trange 0.47 0.55\
 --n_steps 2000000\
 --n_data_points 10\
 --priors_JSON '{"epsilon":["exponential",[0,100]],"sigma":["exponential",[0,5]],"L":["exponential",[0,3]],"Q":["exponential",[0,1]]}'\
 --save_traj 'True'\
 --date 'False'

gen_runfile \
--compound F2\
 --prop All\
 --trange 0.5 0.55\
 --n_steps 2000000\
 --n_data_points 10\
 --priors_JSON '{"epsilon":["exponential",[0,100]],"sigma":["exponential",[0,5]],"L":["exponential",[0,3]],"Q":["exponential",[0,1]]}'\
 --save_traj 'True'\
 --date 'False'

gen_runfile \
--compound C2H2\
 --prop All\
 --trange 0.62 0.7\
 --n_steps 2000000\
 --n_data_points 10\
 --priors_JSON '{"epsilon":["exponential",[0,100]],"sigma":["exponential",[0,5]],"L":["exponential",[0,3]],"Q":["exponential",[0,1]]}'\
 --save_traj 'True'\
 --date 'False'

gen_runfile \
--compound C2H4\
 --prop All\
 --trange 0.41 0.65\
 --n_steps 2000000\
 --n_data_points 10\
 --priors_JSON '{"epsilon":["exponential",[0,100]],"sigma":["exponential",[0,5]],"L":["exponential",[0,3]],"Q":["exponential",[0,1]]}'\
 --save_traj 'True'\
 --date 'False'

gen_runfile \
--compound C2H6\
 --prop All\
 --trange 0.55 0.95\
 --n_steps 2000000\
 --n_data_points 10\
 --priors_JSON '{"epsilon":["exponential",[0,100]],"sigma":["exponential",[0,5]],"L":["exponential",[0,3]],"Q":["exponential",[0,1]]}'\
 --save_traj 'True'\
 --date 'False'
 
#run diatomics rhol+Psat

run_mcmc_sim -f runfiles/O2_rhol+Psat_prod.yml &
run_mcmc_sim -f runfiles/N2_rhol+Psat_prod.yml &
run_mcmc_sim -f runfiles/Br2_rhol+Psat_prod.yml &
run_mcmc_sim -f runfiles/F2_rhol+Psat_prod.yml


#run diatomics All
run_mcmc_sim -f runfiles/O2_All_prod.yml &
run_mcmc_sim -f runfiles/N2_All_prod.yml &
run_mcmc_sim -f runfiles/Br2_All_prod.yml &
run_mcmc_sim -f runfiles/F2_All_prod.yml &


#run C2XX rhol+Psat

run_mcmc_sim -f runfiles/C2H2_rhol+Psat_prod.yml &
run_mcmc_sim -f runfiles/C2H4_rhol+Psat_prod.yml &
run_mcmc_sim -f runfiles/C2H6_rhol+Psat_prod.yml &
run_mcmc_sim -f runfiles/C2F4_rhol+Psat_prod.yml


#run C2XX All
run_mcmc_sim -f runfiles/C2H2_All_prod.yml &
run_mcmc_sim -f runfiles/C2H4_All_prod.yml &
run_mcmc_sim -f runfiles/C2H6_All_prod.yml 

