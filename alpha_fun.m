function alpha = alpha_fun(ep_eff,mp)

gamma = mp(3);
A  = mp(4);
B = mp(5);


alpha = (1/3)*tan(gamma)+ (A/3)*((B*ep_eff^2+ep_eff)/(10^(-4)+ep_eff^2));
end
