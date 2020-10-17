load('check_alg_tan_stiff.mat')


mp = [Gshear,v,gamma_i,A,B];

update_variables(sigma_old,ep_eff_old,delta_eps,Dstar,mp)
