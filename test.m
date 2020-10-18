load('check_update_variables.mat')


mp = [Gshear,v,gamma_i,A,B];

K  = 2*(1+v)/(3*(1-2*v))*Gshear;


sigtrial       = sigma_old + Dstar*delta_eps;
J2trial        = stress_invariant_J2(sigtrial);

update_variables(sigma_old,ep_eff_old,delta_eps,Dstar,mp);
