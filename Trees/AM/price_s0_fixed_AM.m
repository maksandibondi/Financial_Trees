function[optiontree,sharetree,u,d] = price_s0_fixed_AM(N,T,S0,K,r,sigma)

delta_t = T/N;

%% Jarrow-Rudd calibration
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
    V(N+1,i) = max(S(N+1,i)-K,0); % Pay-off function call
    % V(N+1,i) = max(K-S(N+1,i),0); % Pay-off function put
end;

%% Tree filling
for n = N:-1:1
    for i = 1:1:n
       S(n,i) = S0*(u^(i-1))*(d^(n-i));
       V(n,i) = max(max(S(n,i)-K,0),exp(-r*delta_t)*(p*V(n+1,i+1)+(1-p)*V(n+1,i))); % Call
       % V(n,i) = max(max(K-S(n,i),0),exp(-r*delta_t)*(p*V(n+1,i+1)+(1-p)*V(n+1,i))); % Put
    end;
end;

%% Return
optiontree = V;
sharetree = S;

return