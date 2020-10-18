function [sigma,dlambda,ep_eff] = update_variables(sigma_old,ep_eff_old,delta_eps,Dstar,mp)
e1 = [1 0 0 0];
e2 = [0 1 0 0];
e3 = [0 0 1 0];
e4 = [0 0 0 1];


G = mp(1);
v= mp(2);

K  = 2*(1+v)/(3*(1-2*v))*G;


sigtrial       = sigma_old + Dstar*delta_eps
sigkktrial     = sigtrial(1)+sigtrial(2)+sigtrial(3);
sigkktrialvec  = [sigkktrial,sigkktrial,sigkktrial,0]';
I1trial        = stress_invariant_I1(sigtrial);
J2trial        = stress_invariant_J2(sigtrial);
strial  = sigtrial - 1/3*sigkktrialvec;
s_2 = @(dl) strial - 3.*dl*G*strial/sqrt(3*J2trial);
J2_2 = @(dl) 0.5*((e1*s_2(dl))^2 + (e2*s_2(dl))^2 + (e3*s_2(dl))^2 + 2*(e4*s_2(dl))^2);
I1_2 = @(dl) I1trial - 9*K.*dl * alpha_fun(ep_eff_old+dl,mp);

f = yield(sigtrial,ep_eff_old,mp);

if f < 0
    sigma = sigma_old + Dstar*delta_eps;
else

    f_tr = @(dl) sqrt(3*J2_2(dl)) + alpha_fun(ep_eff_old+dl,mp)*I1_2(dl); 
    dlambda  = fzero(f_tr, 1e-5);
    ep_eff  = ep_eff_old + dlambda;
    
    alpha   = alpha_fun(ep_eff,mp);
    
    term = (1/3)*(I1trial-9*K*dlambda*alpha);
    termv= [term,term,term,0]';
    sigma  = strial - 3*G*dlambda*strial/sqrt(3*J2trial) + termv;

end

