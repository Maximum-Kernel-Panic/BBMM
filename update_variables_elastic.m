function [sigma, dlambda, ep_eff] = update_variables_elastic(sigma_old,ep_eff_old,delta_eps,Dstar,mp) 
G = mp(1);
v= mp(2);
K  = 2*(1+v)/(3*(1-2*v))*G;

dsigma = Dstar*delta_eps;

sigma = sigma_old + dsigma;
dlambda = 0;
ep_eff = ep_eff_old;

