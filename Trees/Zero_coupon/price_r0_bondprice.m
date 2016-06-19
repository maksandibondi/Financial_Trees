function[bondpricetree] = price_r0_bondprice(N,T,r0)

delta_t = T/N;

p = 1/2;
R = Interest_rate_tree(r0,N);

%% Final condition
for i = 1:N+1
    P(N+1,i) = 1;
end;

%% Tree filling
for n = N:-1:1
    for i = 1:1:n
        P(n,i) = exp(-R(n,i)*delta_t)*(p*P(n+1,i+1)+(1-p)*P(n+1,i));
    end;
end;

%% Return
bondpricetree = P;

return