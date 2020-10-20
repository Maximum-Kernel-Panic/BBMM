function sig = sigGen(sigH, sigS)

sig = [-1 -1 -1 0]'*sigH + [0 0 0 1]'*sigS;

end