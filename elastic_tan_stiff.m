function Dstar = elastic_tan_stiff(mp)

G = mp(1);
v = mp(2);


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
        Dstar(ij,lk) = 2*G*(0.5*(delta(i,k)*delta(j,l)+delta(i,l)*delta(j,k))+v/(1-2*v)*delta(i,j)*delta(k,l));
    end
end
end