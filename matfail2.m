p0_data = 0.106;
p00_data = 0.062;
p01_data = 0.044;

N=6;
pf1 = p0_data;
pr = 1 - pf1^(1/N);  %Since pf1 = (1-pr)^N

%%
f1w = 1:0.5:8;
pvw = 0.8:0.01:1;
nf = length(f1w);
np = length(pvw);

p00 = zeros(nf, np);
p0 = zeros(nf, np);
rsd =  zeros(nf, np);

for x = 1:nf
    for y = 1:np
        [p0(x,y), p00(x,y)] = calcfail2(f1w(x), pvw(y), pr);
        rsd(x,y) = abs(p00_data - p00(x,y));
    end
end

figure(1); clf; surf(pvw, f1w, rsd);view([0 90]); ylabel('k1'); xlabel('pv'); title('residual'); colorbar;
figure(3); clf; surf(pvw, f1w, p00);view([0 90]); ylabel('k1'); xlabel('pv'); title('P(F1,F2)'); colorbar;
