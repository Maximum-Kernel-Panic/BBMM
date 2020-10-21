clear all
clc
load('check_alg_tan_stiff.mat')



G      = 3e8; %Young modulus, GPa
nu     = 0.3; %Poisson ratio
gamma  = 30*pi/180;  %How to intrepet? Which angle? Need to change to radians
A      = 0.0067;
B      = 48.2;

mp     = [G,nu,gamma,A,B];
K  = 2*(1+v)/(3*(1-2*v))*Gshear;

mp = [Gshear,v,gamma_i,A,B];
%K  = 2*(1+v)/(3*(1-2*v))*Gshear;
>>>>>>> 41b2e5cb0913a7d8ad0aea81b1a194906c376a49

dx = 10e-5
%Correct
I1    = stress_invariant_I1(sigma);
%Correct
J2    = stress_invariant_J2(sigma);

<<<<<<< HEAD
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

    
=======
%Correct
sigkk = sigma(1)+sigma(2)+sigma(3);
epskk = [sigkk,sigkk,sigkk,0]';
s     = sigma-(1/3)*epskk;
sigmadx = [sigma(1)+dx,sigma(2),sigma(3),sigma(4)]';
sigmadx1 = [sigma(1)-dx/2,sigma(2),sigma(3),sigma(4)]';
%dfds is now correct. Write in report why 2*
dfdsii  = (3*s/(2*sqrt(3*J2)) + alpha_fun(ep_eff,mp)*[1,1,1,0]');
dfds    = [dfdsii(1),dfdsii(2),dfdsii(3),2*dfdsii(4)]';


df2ds1test = (3*(sigma(1)+dx-(1/3)*(sigma(1)+sigma(2)+sigma(3)+dx))/(2*sqrt(3*stress_invariant_J2(sigmadx))) +alpha_fun(ep_eff,mp)-dfds(1))/dx
(1/(2*J2))*(2/3*sqrt(3*J2)-3*s(1)*s(1)/(2*sqrt(3*J2)))
(3/2)/sqrt(3*J2)*(-3/2/(3*J2)*s(1)*s(1)+2/3)

for ij = 1:3
    for lk = 1:4
        if ij == 1 || ij ==2 || ij ==3
            i= ij;
            j= ij;
        end
        if ij == 4
            i = 1;
            j = 2;
        end
        if lk == 1 || lk ==2 || lk ==3
            l = lk;
            k = lk;
        end
        if lk == 4
            l = 1;
            k = 2;
        end
        %df2ds2(ij,lk) = (2-delta(l,k))/(4*J2)*((delta(i,l)*delta(j,k)+delta(i,k)*delta(j,l)- ... 
        %    2/3*delta(i,j)*delta(l,k))*sqrt(3*J2)+3*(2-delta(i,j))*s(lk)*s(ij)/sqrt(3*J2));
        %df2ds2(ij,lk)  = 3/(2*sqrt(J2))*(-1/(2*J2)*s(ij)*s(lk)+(delta(i,k)*delta(j,l)-1/3*delta(i,j)*delta(k,l)));    
        df2ds2(ij,lk)  = (3/2)/sqrt(3*J2)*(-(3/2)/(3*J2)*s(ij)*s(lk)+(delta(i,k)*delta(j,l)-1/3*delta(i,j)*delta(k,l)));
    end
