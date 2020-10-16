function J2 = stress_invariant_J2(sigma)

sigkk = sigma(1)+sigma(2)+sigma(3);
epskk = [sigkk,sigkk,sigkk,0]';
s     = sigma-(1/3)*epskk;

J2    = 1/2*(s(1)^2+s(2)^2+s(3)^2+2*s(4)^2);

end