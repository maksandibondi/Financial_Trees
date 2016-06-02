function[optiontree,sharetree,u,d] = price_s0_fixed_EU(N,T,S0,K,r,sigma0)

delta_t = T/N;
sigma_hat = 0.3;

%% Calibration
d = 1+r*delta_t-sigma_hat*sqrt(delta_t);
u = ((1+r*delta_t)^2)/d;
for n = N+1:-1:1
    for i = 1:1:2*n+1
        S(n,i) = S0*(u^((i-1)/2))*(d^(n-1-(i-1)/2));
        sigma(n,i) = min(sigma_hat, sigma0/sqrt(S(n,i)));
        p(n,i) = (sigma(n,i)^2)*delta_t/((u-d)*(1-r*delta_t-1));
        q(n,i) = (sigma(n,i)^2)*delta_t/((u-d)*(1+r*delta_t-d));
    end;
end;

%% Final condition
for i = 1:2*(N+1)+1
    S(N+1,i) = S0*(u^((i-1)/2))*(d^((N-1)-(i-1)/2)); % Share price law
    % V(N+1,i) = max(S(N+1,i)-K,0); % Pay-off function call
    V(N+1,i) = max(K-S(N+1,i),0); % pay-off function put
end;

%% Tree filling
for n = N:-1:1
    for i = 1:1:2*n+1
        V(n,i) = exp(-r*delta_t)*(p(n,i)*V(n+1,i+2)+(1-p(n,i)-q(n,i))*V(n+1,i+1)+q(n,i)*V(n+1,i));
    end;
end;

%% Return
optiontree = V;
sharetree = S;

return