function[optiontree,sharetree,u,d] = price_s0_fixed_Condor(N,T,S0,K,r,sigma)

delta_t = T/N;

%% Tree calibration
u = exp(delta_t*(r-(sigma^2)/2)+sigma*sqrt(delta_t));
d = exp(delta_t*(r-(sigma^2)/2)-sigma*sqrt(delta_t));
p = 0.5;

%% Cox-ross calibration
% u = exp(sigma*sqrt(delta_t));
% d = exp(-sigma*sqrt(delta_t));
% p = (exp(r*delta_t)-d)/(u-d);

%% Final condition
for i = 1:N+1
    S(N+1,i) = S0*(u^(i-1))*(d^(N+1-i)); % Share price law
    
    if S(N+1,i)>=2*K && S(N+1,i)<3*K
       V(N+1,i) = K;
    elseif S(N+1,i)>=3*K && S(N+1,i)<4*K
       V(N+1,i) = max(4*K-S(N+1,i),0);
    elseif S(N+1,i)>=4*K || S(N+1,i)<K
       V(N+1,i) = 0;
    else
       V(N+1,i) = max(S(N+1,i)-K,0); % Pay-off function call
       % V(N+1,i) = max(K-S(N+1,i),0); % pay-off function put
    end;
    
end;

%% Theoretical solution
for n = N:-1:1
    for i = 1:1:n
        S(n,i) = S0*(u^(i-1))*(d^(n-i));
        V(n,i) = exp(-r*delta_t)*(p*V(n+1,i+1)+(1-p)*V(n+1,i));
    end;
end;

%% Return
optiontree = V;
sharetree = S;

return