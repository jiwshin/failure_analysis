function [p0, p00] = calcfail1(f1, pv, pr)
% f1 = forward rate const (1/s)
% pv = vesicle fusion prob
% pr = pv pocc
% assuming that at most one SV transfer during isi

    N = 6;
    pocc = pr/pv; 
    %since pocc = f1 / (f1+b1)
    b1 = f1*(1-pocc)/pocc; %backward rate const (1/s)
    isi = 0.02;
    
    % F1 = fail at 1st pulse
    % F2 = fail at 2nd pulse
    % k = RRP at 1st pulse

    %%
    if pv==1
        pocc_ = pr;
        f1_ = f1;  %forward rate const (1/s)
        fail1_ = (1-pocc_)^N;
        fail2_ = exp(-N*f1_*isi); % prob for no refill
        p00 = fail1_*fail2_;
        p0 = fail1_;
    else
        
        fail1w = zeros(N+1,1); %P(F1| n1=k)*P(n1=k)
        fail2w = zeros(N+1,1);  %P(F2| n1=k)
        
        pn1w = zeros(N+1,1); %prob for RRP = n at 1st pulse
        pin = zeros(N+1,1); % prob for n inc by 1
        pout = zeros(N+1,1); % prob for n dec by 1
        pnull = zeros(N+1,1); % prob for no change in n 
        
        %Calc fail1w = P(F1|n1=k) P(n1=k)
        for idx = 1:N+1
            k = idx-1;
            pn1w(idx) = nchoosek(N,k)*pocc^k*(1-pocc)^(N-k);
            fail1w(idx) =  pn1w(idx)*(1-pv)^k;
        end
    
        efdt = exp(-f1*isi);
        ebdt = exp(-b1*isi);
        
        % Calc fail2w = P(F2|n1=k)
        for idx = 1:N+1
            k = idx-1; % n2 = n1 = k = RRP size at 1st and 2nd pulse d/t no rls at 1st pulse
            pfnull = efdt^(N-k); % no forward transition
            pbnull = ebdt^k;   % no backward transition 
            pin(idx) = 1-pfnull; % n inc 
            pout(idx) = 1-pbnull; % n dec
            pnull(idx) = pfnull*pbnull; % n stays
            fail2w(idx) =  pin(idx)*(1-pv)^(k+1) + pout(idx)*(1-pv)^(k-1) + pnull(idx)*(1-pv)^k;
        end
       
               
        p00 = sum(fail1w.*fail2w);
        p0 = sum(fail1w);
    end

end