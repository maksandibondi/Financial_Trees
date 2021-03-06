clear; clc;
N = 15; T = 0.5;
delta_t = T/N;
vector_t = 0:delta_t:T;

K = 10; S0 = 10; r = 0.3; sigma0 = 0.4; sigmahat = 0.3;

%% Main program
[V,S,u,d,p,q] = price_s0_fixed_EU(N,T,S0,K,r,sigma0);
Initial_price = V(1,1);
display(Initial_price);

%% Greeks
i = 1; step = 0.1;
for newS0 = S0/2:step:1.5*S0
h = 2.5;
V1 = price_s0_fixed_EU(N,T,newS0+h,K,r,sigma0);
V2 = price_s0_fixed_EU(N,T,newS0-h,K,r,sigma0);
V0 = price_s0_fixed_EU(N,T,newS0,K,r,sigma0);
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
    [V1,S1,u1,d1] = price_s0_fixed_EU(N,T,S01(i),K,r,sigma0);
    Initial_price1(i) = V1(1,1);
end;
figure; plot(S01,Initial_price1); hold on;
plot(S(N+1,:),V(N+1,:)); xlabel('variable S'); ylabel('variable V'); legend('initial price', 'final condition');



%% Plotting

%% Share price tree
figure; 
for n = 1:N
    for i = 1:2*n+1
        if S(n,i)==0
           S(n,i)=S(n,n);
        end;
        scatter(vector_t(n),S(n,i)); xlabel('variable t'); ylabel('variable S');
        hold on;
    end;
end;

%% Random Path creation and plotting over share price tree
num_of_iter  = 1000;
for k = 1:num_of_iter
random = rand(N,1);    
path(k,1) = S0;

for n = 2:N
     
   sigma1(n-1) = min(sigmahat,sigma0/sqrt(path(k,n-1)));
   p1(n-1) = sigma1(n-1)^2*delta_t/((u-d)*(u-1-r*delta_t));
   q1(n-1) = sigma1(n-1)^2*delta_t/((u-d)*(1+r*delta_t-d));
   
       if random(n) < p1(n-1) 
           path(k,n)=path(k,n-1)*u;
       elseif random(n) < q1(n-1) 
           path(k,n)=path(k,n-1)*u^((1/2))*d^((1/2));
       else
           path(k,n)=path(k,n-1)*d;
       end;
       
    
%     if random(n) >= 2/3
%         path(k,n) = u*path(k,n-1);
%     elseif random(n) < 2/3 && random(n) >= 1/3
%         path(k,n) = path(k,n-1);
%         % path(k,n) = path(k,n-1)*(u^(1/2))*(d^(1/2));
%     elseif random(n) < 1/3
%         path(k,n) = d*path(k,n-1);
%     end;
    

end;
    
end;


hold on
for m = 1:num_of_iter
plot(vector_t(1:end-1),path(m,:));
end;

%% VaR calculation and plotting
[VaR, pdf, ddf] = VaR(path,num_of_iter,0.95);
display(VaR); step = 0.1;
plot(vector_t(1:end-1),(S0-VaR)*ones(size(vector_t,2)-1),'black'); % Line of share price dropdown maximum
plot(pdf/10+1.1*T,S0-path(1,1):step:S0+2*path(1,1),'LineWidth',0.4); % plotting pdf values on x-axis against share price on y-axis, normalized
plot(ddf/10+T,S0-path(1,1):step:S0+2*path(1,1),'Marker','o','LineWidth',0.2); % plotting ddf values on x-axis against share price on y-axis, normalized
title ('VaR and Random Path');

figure; plot(-path(1,1):step:2*path(1,1),pdf,'black'); hold on;
plot(-path(1,1):step:2*path(1,1),ddf,'LineWidth',0.2); title('VaR'); xlabel('P/L'); ylabel('pdf/ddf');
plot(-VaR*ones(size(pdf,2)-1),pdf(1:end-1),'y'); legend('pdf', 'ddf', '-VaR');

%% Surface

figure;
S_new(1) = 0; 
for k = 2:N+1
    S_new(k) = S_new(k-1)+S0*3/N;
end;

for i = 1:N+1
    for k = 1:N+1
        [V,S,u,d] = price_s0_fixed_EU(N,T-vector_t(i),S_new(k),K,r,sigma0);
        V_new(i,k) = V(1,1);
        k = k+1;
    end;
i = i+1;
end;

surf(vector_t,S_new,V_new); title('price surface'); xlabel('t'); ylabel('S'); zlabel('V');
    