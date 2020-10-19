x = linspace(0,1e9,200);
y = linspace(-1e9,0,200);
[X, Y] = meshgrid(x,y);

G      = 3e8; %Young modulus, GPa
nu     = 0.3; %Poisson ratio
gamma  = 30*pi/180;  %How to intrepet? Which angle? Need to change to radians
A      = 0.0067;
B      = 48.2;

mp     = [G,nu,gamma,A,B];
ep_eff = 0.01;
f = @(J,I) sqrt(3*J) + alpha_fun(ep_eff,mp)*I;
Z = f(X,Y);

surf(X,Y,Z,'LineStyle', 'none');