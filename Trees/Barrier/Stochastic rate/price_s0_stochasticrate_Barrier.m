function[optiontree,sharetree,u,d] = price_s0_stochasticrate_Barrier(N,T,S0,K,R,sigma,B1,B2)

delta_t = T/N;

%% Jarrow-Rudd calibration
u = exp(delta_t*(R(1,1)-(sigma^2)/2)+sigma*sqrt(delta_t));
d = exp(delta_t*(R(1,1)-(sigma^2)/2)-sigma*sqrt(delta_t));
p = 0.5;

%% Cox-ross calibration
% u = exp(sigma*sqrt(delta_t));
% d = exp(-sigma*sqrt(delta_t));
% p = (exp(r*delta_t)-d)/(u-d);

%% Final condition
for i = 1:N+1
    S(N+1,i) = S0*(u^(i-1))*(d^(N+1-i)); % Share price law
    if S(N+1,i) >= B1
        V(N+1,i) = 0;
    elseif S(N+1,i) <= B2
        V(N+1,i) = 0;
    else    
    % V(N+1,i) = max(S(N+1,i)-K,0); % Pay-off function call
    V(N+1,i) = max(K-S(N+1,i),0); % pay-off function put
    end;
end;

%% Tree filling
for n = N:-1:1
    for i = 1:1:n
        S(n,i) = S0*(u^(i-1))*(d^(n-i));
        if S(n,i) > B1
            V(n,i) = 0;
        elseif S(n,i) < B2
            V(n,i) = 0;
        else
        V(n,i) = exp(-R(n,i)*delta_t)*(p*V(n+1,i+1)+(1-p)*V(n+1,i));
        end
       
    end;
end;

%% Return
optiontree = V;
sharetree = S;

return