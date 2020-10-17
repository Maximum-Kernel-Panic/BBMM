load('check_alg_tan_stiff.mat')


mp = [Gshear,v,gamma_i,A,B];

%Correct
I1    = stress_invariant_I1(sigma);
%Correct
J2    = stress_invariant_J2(sigma);

%Correct
sigkk = sigma(1)+sigma(2)+sigma(3);
epskk = [sigkk,sigkk,sigkk,0]';
s     = sigma-(1/3)*epskk;

%Correct
dfdK  = I1;
%dfds is now correct. Write in report why 2*
dfdsii  = (3*s/(2*sqrt(3*J2)) + alpha_fun(ep_eff,mp)*[1,1,1,0]');
dfds    = [dfdsii(1),dfdsii(2),dfdsii(3),2*dfdsii(4)];


df2dsdK = [1,1,1,0];

dKdk  = A/(3*(10^(-4)+dlambda^2)^2)*(-dlambda^2+2*10^(-4)*B*dlambda+10^(-4));
f  = yield(sigma,ep_eff,mp);

dx = 1E-2;

sigma_dis = [sigma(1)+dx,sigma(2),sigma(3),sigma(4)]';
sigma_dis_1 = [sigma(1)-dx,sigma(2),sigma(3),sigma(4)]';
sigma = [sigma(1),sigma(2),sigma(3),sigma(4)]';

