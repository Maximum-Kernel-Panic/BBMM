function [ev, mult] = heltalsev(A,tol)


if nargin < 2
    tol = 5e-2;
end

evf = eig(A);
multm = [];
evm = [];

for i = 1:length(evf)
    if abs(round(evf(i))-evf(i)) < tol
        evm(i) = round(evf(i)); 
        multm(i) = length(find(abs((evf-evm(i)))<tol));
    end
end

ev = [];
mult = [];

if length(evm) > 1
   for i = 1:length(evm) 
       if length(find(abs(ev-evm(i))<tol)) == 0
           ev(i) = evm(i);
           mult(i) = multm(i);
       end
   end
end

