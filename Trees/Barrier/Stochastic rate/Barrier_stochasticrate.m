clear; clc;
N = 20; T = 0.5; 
delta_t = T/N;
vector_t = 0:delta_t:T;

K = 10; S0 = 10; r0 = 0.3; sigma = 0.4;
B1 = 13; B2 = 7;

%% Main program
R = Interest_rate_tree(r0,N);
[V,S,u,d] = price_s0_stochasticrate_Barrier(N,T,S0,K,R,sigma,B1,B2);
Initial_price = V(1,1);
display(Initial_price);


%% Greeks
i = 1; step = 0.1;
for newS0 = S0/2:step:1.5*S0
h = 2.5;
V1 = price_s0_stochasticrate_Barrier(N,T,newS0+h,K,R,sigma,B1,B2);
V2 = price_s0_stochasticrate_Barrier(N,T,newS0-h,K,R,sigma,B1,B2);
V0 = price_s0_stochasticrate_Barrier(N,T,newS0,K,R,sigma,B1,B2);
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
    [V1,S1,u1,d1] = price_s0_stochasticrate_Barrier(N,T,S01(i),K,R,sigma,B1,B2);
    Initial_price1(i) = V1(1,1);
end;
figure; plot(S01,Initial_price1); hold on;
plot(S(N+1,:),V(N+1,:)); xlabel('variable S'); ylabel('variable V'); legend('initial price', 'final condition');



%% Plotting

%% Share price tree
figure; 
for n = 1:N
    for i = 1:n
        scatter(vector_t(n),S(n,i)); xlabel('variable t'); ylabel('variable S'); title('share price tree');
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
for m = 1:num_of_iter
plot(vector_t(1:end-1),path(m,:));
end;
hold on
plot(vector_t(1:end-1),B1*ones(size(vector_t,2)-1),'black');
plot(vector_t(1:end-1),B2*ones(size(vector_t,2)-1),'black');

%% Interest rate price tree plotting
figure; 
for n = 1:N
    for i = 1:n
        scatter(vector_t(n),R(n,i)); xlabel('variable t'); ylabel('variable R'); title ('interest rate tree');
        hold on;
    end;
end;

%% Surface

figure;
S_new(1) = 0; 
for k = 2:N+1
    S_new(k) = S_new(k-1)+S0*3/N;
end;

for i = 1:N+1
    for k = 1:N+1
        [V,S,u,d] = price_s0_stochasticrate_Barrier(N,T-vector_t(i),S_new(k),K,R,sigma,B1,B2);
        V_new(i,k) = V(1,1);
        k = k+1;
    end;
i = i+1;
end;

surf(vector_t,S_new,V_new); title('price surface'); xlabel('t'); ylabel('S'); zlabel('V');
