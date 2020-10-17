function [sigma,dlambda,ep_eff] = update_variables(sigma_old,ep_eff_old,delta_eps,Dstar,mp)

G = mp(1);
v= mp(2);
gamma = mp(3)*pi/180;
A  = mp(4);
B = mp(5);

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

% delta_eps_new = [delta_eps(1),delta_eps(2),0,delta_eps(3)]';




if f < 0
    sigma = sigma_old + Dstar*delta_eps;
else
    alph = @(dl) (1/3)*(A*(B*(ep_eff_old+dl).^2 + ep_eff_old + dl)./(10^(-4)+(ep_eff_old+dl).^2)+ tan(gamma));
    f_tr = @(dl) sqrt(3*J2) + alph(dl)*I1;
    dlambda = fzero(f_tr,1e-4)
%    x = linspace(-0.5, 0.5, 1000);
%    f_tr_val = f_tr(x);
    ep_eff  = ep_eff_old + dlambda;
    
    alpha   = alpha_fun(ep_eff,mp);
    
    term = (1/3)*(I1trial-9*K*dlambda*alpha);
    termv= [term,term,term,0]';
    sigma  = strial-3*G*dlambda*strial/sqrt(3*J2trial)+ termv;

end

