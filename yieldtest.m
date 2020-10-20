x = linspace(-500,500,1000);
f = zeros(1000,1);

G      = 3e8; %Young modulus, GPa
nu     = 0.3; %Poisson ratio
gamma  = 30*pi/180;  %How to intrepet? Which angle? Need to change to radians
A      = 0.0067;
B      = 48.2;

mp     = [G,nu,gamma,A,B];
hold on;
for j = 0:0
    ep_eff = j*0.001;
    for i = 1:1000
        f(i) = yield([-1000,-1000,-1000,x(i)]',ep_eff,mp);
    
    end
    plot(x,f);
end
