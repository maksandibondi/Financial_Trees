clear; clc;
N = 20; T = 0.5;
delta_t = T/N;
vector_t = 0:delta_t:T;

r0 = 0.5;

%% Main program
P = price_r0_bondprice(N,T,r0);
Initial_price = P(1,1);
display(Initial_price);

%% Checking final condition at t = T and dependence S of V at t = 0
T1 = 0:0.05:6;
for i = 1:size(T1,2);
    temp = price_r0_bondprice(N,T1(i),r0);
    bond_price1(i) = temp(1,1);
    bond_fixe_price(i) = exp(-r0*T1(i));
end;
figure; plot(T1,bond_price1,'black'); hold on;
plot(T1,bond_fixe_price,'red'); title('bond price curve'); xlabel('T'); ylabel('price');

%% Plotting

%% Bond price tree
figure; 
for n = 1:N+1
    for i = 1:n
        scatter(vector_t(n),P(n,i)); xlabel('variable t'); ylabel('variable P'); title('bond price tree');
        hold on;
    end;
end;

