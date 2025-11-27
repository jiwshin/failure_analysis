function [p0, p00] = calcfail2(f1, pv, pr)
% f1 = forward rate const (1/s)
% pv = vesicle fusion prob
% pr = pv pocc

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
        pn2 = zeros(N+1, N+1); % p(n1, n2) matrix = p(n2|n1)

        % Calc fail2w = P(F2|n1=k)
        for idx = 1:N+1
            k = idx-1; % n1 = k = RRP size at 1st  pulse 
            n2 = k;
            while n2<=N
                 pn2(k+1, n2+1) = calcpn2(k, n2, f1, b1, N, isi); %=p(n2|n1=k)
                 n2 = n2+1;
            end
            n2 = k-1;
            while n2>=0
                 pn2(k+1, n2+1) = calcpn2(k, n2, f1, b1, N, isi);
                 n2 = n2-1;
            end
        end

        for idx=1:N+1  %n1
            for jdx = 1:N+1  %n2
                n2 = jdx - 1;
                fail2w(idx) = fail2w(idx) + pn2(idx,jdx)*(1-pv)^n2;
            end
        end
             
        p00 = sum(fail1w.*fail2w);
        p0 = sum(fail1w);        
  
    end

end