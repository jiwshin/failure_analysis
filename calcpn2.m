function [pn2] = calcpn2(n1, n2, k1, b1, N, dt)
%n1, SV num at 1st pulse
%n2, SV num at 2nd pulse
%k1,b1: forward and backward rate const
%N, docking sites
%dt, inter-spike interval
%Calc  P(n2|n1)

    ne = N-n1; % empty sites
    nmv = n2 - n1; %nmv SVs are transferred during dt

    if nmv>0
        pn2 = nchoosek(ne, nmv)*(1-exp(-k1*dt));
    elseif nmv<0
        pn2 = nchoosek(n1, -nmv)*(1-exp(-b1*dt));
    else % no move
        pn2 = exp(-ne*k1*dt)*exp(-n1*b1*dt);
    end
end