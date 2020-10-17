function f = yield(sigma, ep_eff, mp)

G     = mp(1);
nu    = mp(2);
gamma = mp(3);
A     = mp(4);
B     = mp(5);


I1      = stress_invariant_I1(sigma);
J2      = stress_invariant_J2(sigma);
alpha   = alpha_fun(ep_eff,mp);

f       = sqrt(3*J2) + alpha*I1;

end