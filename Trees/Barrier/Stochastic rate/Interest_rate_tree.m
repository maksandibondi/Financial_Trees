function [r] = Interest_rate_tree(r0,N)

%% Jarrow-Rudd calibration
p = 1/2; u = 1.1; d = 0.9;

%% Tree filling
for n = N+1:-1:1
    for i = 1:1:n
    r(n,i) = r0*(u^(n-1))*(d^(n-i));
    end;
end;

return;
    