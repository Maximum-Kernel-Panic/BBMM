clc 
clear all

load('check_alg_tan_stiff.mat');

G      = 3e8; %Young modulus, GPa
nu     = 0.3; %Poisson ratio
gamma  = 30*pi/180;  %How to intrepet? Which angle? Need to change to radians
A      = 0.0067;
B      = 48.2;

mp     = [G,nu,gamma,A,B];
Dats_2 = alg_tan_stiff(sigma, dlambda, ep_eff,Dstar,mp)
% [sigma_2, dlambda_2, ep_eff_2] = update_variables(sigma_old, ep_eff_old, delta_eps, Dstar, mp);