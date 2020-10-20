load('check_update_variables.mat')


G      = 3e8; %Young modulus, GPa
nu     = 0.3; %Poisson ratio
gamma  = 30*pi/180;  %How to intrepet? Which angle? Need to change to radians
A      = 0.0067;
B      = 48.2;

mp     = [G,nu,gamma,A,B];
K  = 2*(1+v)/(3*(1-2*v))*Gshear;


% sigtrial       = sigma_old + Dstar*delta_eps;
% J2trial        = stress_invariant_J2(sigtrial);
% update_variables(sigma_old,ep_eff_old,delta_eps,Dstar,mp);

nH = linspace(0, 1e8, 1000);
nS = linspace(0, 1e8, 1000);
f = zeros(1000);
[pH, +pS] = meshgrid(nH,nS);

for i = 1:1000
    for j = 1:1000
        sig = sigGen(i,j);
        f(i,j) = yield(sig, 0.01, mp);
    end
end

hold on
surf(f)
shading interp



% testnorm = zeros(1000,1);
% J2 = zeros(1000,1);
% I1 = zeros(1000,1);
% 
% for i = 1:1000
%    testsig(i,:) = sigbase*i*1e6; 
%    testnorm(i) = norm(testsig(i,:));
%    J2(i) = stress_invariant_J2(testsig(i,:)');
%    I1(i) = stress_invariant_I1(testsig(i,:)');
%    f(i) = yield(testsig(i,:)',1,mp);
% end
% hold on;
% % plot(testnorm,sqrt(3.*J2),'.');
% % plot(testnorm,I1,'o');
% plot(testnorm,f);
% Dats = alg_tan_stiff(sigmatest,0,0,Dstar,mp);
% Dstar-Dats;

% x = linspace(0, 1e-2, 1000);
% dx = 1e-7;
% K = @(k) A/3 * (B*k.^2+k)/(1e-4+k.^2);
% dKdk_alg_fun= @(dlambda) A/(3*(10^(-4)+dlambda^2)^2)*(-dlambda^2+2*10^(-4)*B*dlambda+10^(-4));
% dKdk_num = zeros(1000,1);
% dKdk_alg = zeros(1000,1);
% K_alg = zeros(1000);
% for i = 1:1000
%     dKdk_num(i) = (K(x(i)+dx) - K(x(i)))/dx;
%     K_alg = K(x(i));
%     dKdk_alg(i) = dKdk_alg_fun(x(i));
% end
% 
% hold on;
% plot(x,dKdk_num,'o');
% plot(x,dKdk_alg,'.');

    