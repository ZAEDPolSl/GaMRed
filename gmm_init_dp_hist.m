function [alpha,mu,sigma] = gmm_init_dp_hist(data,K)
% GMM_INIT_DP(data,K)
% Compute initial conditions for GMM by using dynamic programming for
% approximate signal (by operation of binning).
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
n_bins = round(max(10,min(100,sqrt(length(data)))));
[y,x] = hist(data,n_bins);      %create binned data (approximate signal)
s_corr = ((x(2) - x(1))^2)/12;  %sheppards correction for binned data
K = K - 1; 
N = length(x);
p_opt_idx = zeros(1,N);
p_aux = zeros(1,N);
opt_pals = zeros(K,N);
for a=1:N
    invec = x(a:N);
    yinvec = y(a:N);
    if sum(yinvec) <= par2
        p_opt_idx(a) = inf;
    else
        wwec=yinvec/(sum(yinvec));
        var_bin = sum(((invec-sum(invec.*wwec)).^2).*wwec);
        if var_bin > s_corr
            p_opt_idx(a)=(par1+sqrt(var_bin-s_corr))/(max(invec)-min(invec));
        else
            p_opt_idx(a) = inf;
        end
   end
end

% aux_mx
aux_mx = zeros(N,N);
for a=1:N-1
    for b=a+1:N
        invec = x(a:b-1);
        yinvec = y(a:b-1);
        if sum(yinvec) <= par2
            aux_mx(a,b) = inf;
        else
            wwec = yinvec/(sum(yinvec));
            var_bin = sum(((invec-sum(invec.*wwec)).^2).*wwec);
            if var_bin > s_corr
                aux_mx(a,b) = (par1+sqrt(var_bin-s_corr))/(max(invec)-min(invec));
            else
                aux_mx(a,b) = inf;
            end
            
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
alpha = zeros(1,K); mu = alpha; sigma = alpha;
for a=1:K+1
    invec = x(opt_part(a):opt_part(a+1)-1);
    yinvec = y(opt_part(a):opt_part(a+1)-1);
    wwec = yinvec/(sum(yinvec));
    alpha(a) = sum(yinvec)/sum(y);
    mu(a) = sum(invec.*wwec);
    sigma(a)= sqrt(sum(((invec-sum(invec.*wwec)).^2).*wwec)-s_corr);
end