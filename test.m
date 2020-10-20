clear all
clc
load('check_alg_tan_stiff.mat')


mp = [Gshear,v,gamma_i,A,B];
%K  = 2*(1+v)/(3*(1-2*v))*Gshear;

dx = 10e-5
%Correct
I1    = stress_invariant_I1(sigma);
%Correct
J2    = stress_invariant_J2(sigma);

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
end