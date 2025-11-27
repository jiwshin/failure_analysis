global N pr q q2 k1

N = 6;
pr = 0.31;

q = 25.7; %pA
q2 = q*q; 
k1 = 4; %1/s

% Note that Min of pv = pr
pvw = pr:0.01:1;  [covi, covir] = covar_(pvw);
figure(1); plot(pvw, covi, pvw, covir);

% cov = N pr N pr (1 - pv) q2 / N
%        = N pr2 (1 - pv) q2 
% For const N, pr and q, cov is a linear function of pv
% Since N = 6, q = 25.7, pr = 0.312;
% N*pr*pr*q*q = 386 pA2
