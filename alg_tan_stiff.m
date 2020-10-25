function Dats = alg_tan_stiff(sigma,dlambda,ep_eff,Dstar,mp)

A  = mp(4);
B  = mp(5);

%Correct
I1    = stress_invariant_I1(sigma);
%Correct
J2    = stress_invariant_J2(sigma);

%Correct
sigkk = sigma(1)+sigma(2)+sigma(3);
epskk = [sigkk,sigkk,sigkk,0]';
s     = sigma-(1/3)*epskk;

%Correct
dfdK    = I1;
%dfds is now correct. Write in report why 2*
dfdsii  = (3*s/(2*sqrt(3*J2)) + alpha_fun(ep_eff,mp)*[1,1,1,0]');
dfds    = [dfdsii(1),dfdsii(2),dfdsii(3),2*dfdsii(4)]';
df2dsdK = [1,1,1,0]';

%dKdk  = A/(3*(10^(-4)+dlambda^2)^2)*(-dlambda^2+2*10^(-4)*B*dlambda+10^(-4));

dKdk  = A/(3*(10^(-4)+(dlambda+ep_eff)^2)^2)*(-(dlambda+ep_eff)^2+2*10^(-4)*B*(dlambda+ep_eff)+10^(-4));


for ij = 1:4
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
        if ij == 4 && lk == 4
            factor = 4; 
            factor1 =2;
        elseif ij == 4 || lk == 4
            factor = 2;
            factor1=1;
        else
            factor = 1;
            factor1 = 1;
        end
        df2ds2(ij,lk)  = (3/2)/sqrt(3*J2)*(-(3/2)/(3*J2)*s(ij)*s(lk)*factor + factor1*(delta(i,l)*delta(j,k)-1/3*delta(i,j)*delta(k,l)));
    end
end


Da      = dKdk;


Davec   = ((Dstar\eye(4) + dlambda*df2ds2))\eye(4);


dgdsiga = dfds  + dlambda*df2dsdK*Da;
H       = -dfdK*Da;
Aa      = dfds'*Davec*dfds+H;

if dlambda == 0
    Dats       = Dstar;
else
    Dats       = (Davec - 1/Aa*Davec*dgdsiga*dfds'*Davec);
end

end