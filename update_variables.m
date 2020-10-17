function [sigma,dlambda,ep_eff] = update_variables(sigma_old,ep_eff_old,delta_eps,Dstar,mp)

G = mp(1);
v= mp(2);

K  = 2*(1+v)/(3*(1-2*v))*G;


sigtrial       = sigma_old + Dstar*delta_eps;
sigkktrial     = sigtrial(1)+sigtrial(2)+sigtrial(3);
sigkktrialvec  = [sigkktrial,sigkktrial,sigkktrial,0]';
I1trial        = stress_invariant_I1(sigtrial);
J2trial        = stress_invariant_J2(sigtrial);

strial         = sigtrial - 1/3* sigkktrialvec;



f = yield(sigtrial,ep_eff_old,mp);

% I1 = stress_invariant_I1(sigma_old);
% J2 = stress_invariant_J2(sigma_old);




if f < 0
    sigma = sigma_old + Dstar*delta_eps;
else
    %f_tr = @(dl) sqrt(3*J2) + 1/3*(A*(B*(ep_eff_old+dl)^2 + ep_eff_old + dl)/(1e-4+(ep_eff_old+dl)^2) + tan(gamma*pi/180))*I1;
    %dlambda = fzero(f_tr,1e-6)
    %ep_eff = ep_eff_old + dlambda
    
    dlambda = 1.0067e-4;
    ep_eff  = ep_eff_old + dlambda;
    
    alpha   = alpha_fun(ep_eff,mp);
    
    term = (1/3)*(I1trial-9*K*dlambda*alpha);
    termv= [term,term,term,0]';
    sigma  = strial-3*G*dlambda*strial/sqrt(3*J2trial)+ termv;
end

