clear; clc;
N = 20; T = 0.5; 
delta_t = T/N;
vector_t = 0:delta_t:T;

K = 10; S0 = 10; r = 0.3; sigma = 0.4;
B1 = 17; B2 = 9;

%% Main program
[V,S,u,d] = price_s0_fixed_Barrier(N,T,S0,K,r,sigma,B1,B2);
Initial_price = V(1,1);
display(Initial_price);

%% Greeks
i = 1; step = 0.1;
for newS0 = S0/2:step:1.5*S0
h = 0.5;
V1 = price_s0_fixed_Barrier(N,T,newS0+h,K,r,sigma,B1,B2);
V2 = price_s0_fixed_Barrier(N,T,newS0-h,K,r,sigma,B1,B2);
V0 = price_s0_fixed_Barrier(N,T,newS0,K,r,sigma,B1,B2);
DELTA(i) = (V1(1,1)-V2(1,1))/(2*h);
GAMMA(i) = (V1(1,1)+V2(1,1)-2*V0(1,1))/(h^2);
i = i+1;
end;
display(DELTA(S0/(2*step)+1));
display(GAMMA(S0/(2*step)+1));
figure; plot(S0/2:step:1.5*S0,DELTA); xlabel('variable S0'); ylabel('variable DELTA'); title('DELTA');
figure; plot(S0/2:step:1.5*S0,GAMMA); xlabel('variable S0'); ylabel('variable GAMMA'); title('GAMMA');

%% Checking final condition at t = T and dependence S of V at t = 0
S01 = 0:1:20;
for i = 1:size(S01,2);
    [V1,S1,u1,d1] = price_s0_fixed_Barrier(N,T,S01(i),K,r,sigma,B1,B2);
    Initial_price1(i) = V1(1,1);
end;
figure; plot(S01,Initial_price1); hold on;
plot(S(N+1,:),V(N+1,:)); xlabel('variable S'); ylabel('variable V');



%% Plotting

%% Share price tree
figure; 
for n = 1:N
    for i = 1:n
        scatter(vector_t(n),S(n,i)); xlabel('variable t'); ylabel('variable S');
        hold on;
    end;
end;

%% Random Path creation and plotting over share price tree
num_of_iter  = 50;
for q = 1:num_of_iter
random = rand(N,1);    
path(q,1) = S0;

for n = 2:N
    
    if random(n)>=0.5
        path(q,n) = u*path(q,n-1);
    else
        path(q,n) = d*path(q,n-1);
    end;
    
end;

end;
hold on
for m = 1:50
plot(vector_t(1:end-1),path(m,:));
end;
hold on
plot(vector_t(1:end-1),B1*ones(size(vector_t,2)-1),'black');
plot(vector_t(1:end-1),B2*ones(size(vector_t,2)-1),'black');

%% Surface
figure;
surf(S(N+1,:),vector_t,V); xlabel('variable S'); ylabel('variable t'); zlabel('variable V');
