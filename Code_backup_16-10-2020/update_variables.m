function [sigma,dlambda,ep_eff] = update_variables(sigma_old,ep_eff_old,delta_eps,Dstar,mp)

G = mp(1);
nu= mp(2);
gamma = mp(3);
A  = mp(4);
B = mp(5);


f = yield(sigma,ep_eff_old,mp);
plastic = true;
G  = mp(1);

sigkko= sigma_old(1)+sigma_old(2)+sigma_old(3);
epskko= [sigkko,sigkko,sigkko,0]';
sold  = sigma_old-(1/3)*epskko;

delta_eps_new = [delta_eps(1),delta_eps(2),0,delta_eps(3)];

I1trial = stress_invariant_I1(sigma_old + Dstar*delta_eps_new);
J2trial = stress_invariant_J2(sigma_old + Dstar*delta_eps_new);


if f < 0
    plastic = false;
end

if plastic
    strial = sold + Dstar*delta_eps-1/3*Dstar*(delta_eps(1)+delta_eps(2))*[1,1,1,0];
    sigma = strial-3G*dlambda*strial/sqrt(3*J2trial)+ (1/3)*(I1trial-9K*dlambda*alpha)*sigma_old    
else
    sigma = sigma_old + Dstar*delta_eps_new;
end

