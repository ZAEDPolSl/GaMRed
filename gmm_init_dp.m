function [alpha,mu,sigma] = gmm_init_dp(x,K)
% GMM_INIT_DP(data,K)
% Compute initial conditions for GMM by using dynamic programming 
% Input:
% data - sample to partition [1xn]
% K - number of partitions [1x1]
% Output:
% alpha - weights
% mu - means
% sigma - standard deviations

%parameters
par1 = 0.1; %for robustness (fine for data in range 0-20)
par2 = 5;   %minimum no. of points in signal fragment

% initialize
x = sort(x);  %input vector must be sorted
K = K - 1;
N = length(x);
p_opt_idx = zeros(1,N);
p_aux = zeros(1,N);
opt_pals = zeros(K,N);
for a=1:N
    invec = x(a:N);
    if length(invec) < par2
        p_opt_idx(a) = inf;
    else
        p_opt_idx(a) = (par1+std(invec))/(max(invec)-min(invec));
    end
end

% aux_mx
aux_mx = zeros(N,N);
for a=1:N-1
   for b=a+1:N
       invec = x(a:b-1);
       if length(invec) < par2
            aux_mx(a,b) = inf;
        else
            aux_mx(a,b) = (par1+std(invec))/(max(invec)-min(invec));
        end
   end
end

%iterate
for kster = 1:K
   % kster
   for a=1:N-kster
       for b=a+1:N-kster+1
           p_aux(b) =  aux_mx(a,b) + p_opt_idx(b);
       end
       [mm,ix] = min(p_aux(a+1:N-kster+1));
       p_opt_idx(a) = mm;
       opt_pals(kster,a) = a + ix(1);
   end
end

%restore optimal decisions
opt_part = zeros(1,K);
opt_part(1) = opt_pals(K,1);
for kster = K-1:-1:1
   opt_part(K-kster+1) = opt_pals(kster,opt_part(K-kster));
end

%find initial conditions
opt_part=[1 opt_part N+1];             
alpha = zeros(K,1); mu = alpha; sigma = alpha;
for a=1:K+1
    alpha(a) = (opt_part(a+1)-opt_part(a))/N;
    mu(a) = mean(x(opt_part(a):opt_part(a+1)-1));   
    sigma(a) = std(x(opt_part(a):opt_part(a+1)-1));  
end