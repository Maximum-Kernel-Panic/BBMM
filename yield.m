function f = yield(sigma, ep_eff, mp)


I1      = stress_invariant_I1(sigma);
J2      = stress_invariant_J2(sigma);
alpha   = alpha_fun(ep_eff,mp);

f       = sqrt(3*J2) + alpha*I1;

end