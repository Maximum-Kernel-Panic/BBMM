function [sigma,dlambda,ep_eff] = update_variables(sigma_old,ep_eff_old,delta_eps,Dstar,mp)
e1 = [1 0 0 0];
e2 = [0 1 0 0];
e3 = [0 0 1 0];
e4 = [0 0 0 1];


G = mp(1);

v= mp(2);
gamma = mp(3)*pi/180;
A  = mp(4);
B = mp(5);
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

I1 = stress_invariant_I1(sigma_old);
J2 = stress_invariant_J2(sigma_old);

% delta_eps_new = [delta_eps(1),delta_eps(2),0,delta_eps(3)]';

% f = yield(sigma_old+Dstar*delta_eps,ep_eff_old,mp);
% I1 = stress_invariant_I1(sigma_old);
% J2 = stress_invariant_J2(sigma_old);
% I1trial = stress_invariant_I1(sigma_old + Dstar*delta_eps)
% J2trial = stress_invariant_J2(sigma_old + Dstar*delta_eps)




if f < 0
    sigma = sigma_old + Dstar*delta_eps;
else

    f_tr = @(dl) sqrt(3*J2_2(dl)) + alpha_fun(ep_eff_old+dl,mp)*I1_2(dl); 
    dlambda  = fzero(f_tr, 1e-4)
    xstart = -1e-3;
    xend = 1e-2;
    N = 1000;
    x = linspace(xstart,xend, N);
    f_tr_val = zeros(N,1);
    for i = 1:N
        f_tr_val(i) = f_tr(x(i));
    end

    plot(x,f_tr_val);
    dlambda = 1.0067e-4;
    ep_eff  = ep_eff_old + dlambda;
    
    alpha   = alpha_fun(ep_eff,mp);
    
    term = (1/3)*(I1trial-9*K*dlambda*alpha);
    termv= [term,term,term,0]';
    sigma  = strial-3*G*dlambda*strial/sqrt(3*J2trial)+ termv;

end

