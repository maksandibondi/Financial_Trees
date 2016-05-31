clear; clc;
N = 15; T = 0.5;
delta_t = T/N;
vector_t = 0:delta_t:T;

K = 10; S0 = 10; r = 0.3; sigma0 = 0.4;

%% Main program
[V,S,u,d] = price_s0_fixed_EU(N,T,S0,K,r,sigma0);
Initial_price = V(1,1);
display(Initial_price);

%% Greeks
h = 0.5;
V1 = price_s0_fixed_EU(N,T,S0+h,K,r,sigma0);
V2 = price_s0_fixed_EU(N,T,S0-h,K,r,sigma0);
V0 = price_s0_fixed_EU(N,T,S0,K,r,sigma0);
DELTA = (V1(1,1)-V2(1,1))/(2*h);
GAMMA = (V1(1,1)+V2(1,1)-2*V0(1,1))/(h^2);
display(DELTA);
display(GAMMA);

%% Checking final condition at t = T and dependence S of V at t = 0
S01 = 0:1:20;
for i = 1:size(S01,2);
    [V1,S1,u1,d1] = price_s0_fixed_EU(N,T,S01(i),K,r,sigma0);
    Initial_price1(i) = V1(1,1);
end;
plot(S01,Initial_price1); hold on;
plot(S(N+1,:),V(N+1,:)); xlabel('variable S'); ylabel('variable V');






%% Plotting

%% Share price tree
figure; 
for n = 1:N
    for i = 1:2*n+1
        scatter(vector_t(n),S(n,i)); xlabel('variable t'); ylabel('variable S');
        hold on;
    end;
end;

%% Random Path creation and plotting over share price tree
num_of_iter  = 10;
for q = 1:num_of_iter
random = rand(N,1);    
path(q,1) = S0;

for n = 2:N
    
    if random(n) >= 2/3
        path(q,n) = u*path(q,n-1);
    elseif random(n) < 2/3 && random(n) >= 1/3
        path(q,n) = path(q,n-1);
    elseif random(n) < 1/3
        path(q,n) = d*path(q,n-1);
    end;
end;
    
end;

hold on
for m = 1:num_of_iter
plot(vector_t(1:end-1),path(m,:));
end;

%% Surface
figure;
surf(S(N+1,:),vector_t,V); xlabel('variable S'); ylabel('variable t'); zlabel('variable V');
