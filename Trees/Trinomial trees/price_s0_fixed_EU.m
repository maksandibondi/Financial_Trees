function[optiontree,sharetree,u,d,p,q] = price_s0_fixed_EU(N,T,S0,K,r,sigma0)

delta_t = T/N;
sigma_hat = 0.3;

%% Calibration
d = 1+r*delta_t-sigma_hat*sqrt(delta_t);
u = ((1+r*delta_t)^2)/d;
S(1,1) = S0;
for n = 2:N+1
        for i = 1:1:2*n-1
            S(n,i) = S0*(u^((i-1)/2))*(d^(n-(i+1)/2));
        	sigma(n,i) = min(sigma_hat, sigma0/sqrt(S(n,i)));
            p(n,i) = (sigma(n,i)^2)*delta_t/((u-d)*(u-r*delta_t-1));
            q(n,i) = (sigma(n,i)^2)*delta_t/((u-d)*(1+r*delta_t-d));
        end;    
end;
sigma(1,1) = min(sigma_hat, sigma0/sqrt(S(1,1)));
p0 = (sigma(1,1)^2)*delta_t/((u-d)*(u-r*delta_t-1));
q0 = (sigma(1,1)^2)*delta_t/((u-d)*(1+r*delta_t-d));
temp = [p0 zeros(1,2*N)];
p = [temp; p(2:end,:)];
temp = [q0 zeros(1,2*N)];
q = [temp; q(2:end,:)];

%% Final condition
for i = 1:2*(N+1)-1
    V(N+1,i) = max(S(N+1,i)-K,0); % Pay-off function call
    % V(N+1,i) = max(K-S(N+1,i),0); % pay-off function put
end;

%% Tree filling
for n = N:-1:1
    for i = 1:1:2*n-1
        V(n,i) = exp(-r*delta_t)*(p(n,i)*V(n+1,i+2)+(1-p(n,i)-q(n,i))*V(n+1,i+1)+q(n,i)*V(n+1,i));
    end;
end;

%% Return
optiontree = V;
sharetree = S;

return