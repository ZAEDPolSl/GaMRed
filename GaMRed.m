function [thr,bic,stats] = GaMRed(data,K,K_noise,draw,SW)
%GaMRed  Estimating noise components using Gaussian mixture model
% Input:
% data - vector of data (1xn)
% K - number of Gaussian components
% K_noise - number of noise components from left
% draw - if draw results
% Output:
% thr - noise threshold
% bic - Bayesian Information Criterion for estimated model
% stats - additional statistics
% 
% Author: Michal Marczyk
% Michal.Marczyk@polsl.pl

if nargin < 3
    error('Insufficient number of arguments.')
end
if nargin < 4
    draw = 0;
end

data = sort(data(:));
it = 1;                 %no. of EM iterations
N = length(data);       %nb of measurements
epsilon = 1e-6;         %EM stop criterion threshold
max_iter = 10000;       %max. no. of iterations
bic = Inf;
% SW = 0.001;

while bic == Inf || bic == 0
    
    %initial conditions
    if K==1
        alpha = 1;
        mi = mean(data);
        sigma = std(data);
    else
        [alpha,mi,sigma] = gmm_init_dp_hist(data,K);
    end
    
    f = zeros(N,K); p = f; L_old = 1;
    for k=1:K; f(:,k) = alpha(k) * normpdf(data,mi(k),sigma(k)); end
    f(isnan(f)) = 1e-45; f(f==0) = 1e-45;   %for numerical stability
        
    if draw
        disp('Starting values')
        disp(num2str([alpha; mi; sigma]))
    end

    %main loop
    while it < max_iter
        px = sum(f,2);
	px(isnan(px) | px==0) = 5e-324;
        L_new = sum(log(px),1);
        
        %stop criterion
        if 100*abs((L_new - L_old)/L_old) < epsilon || isnan(L_new) || isinf(L_new)
            break
        end
        for k=1:K
            p(:,k) = f(:,k)./px;
            sum_pk = sum(p(:,k),1);
            alpha(k) = sum_pk/N;
            mi(k) = sum(data.*p(:,k),1)/sum_pk;
            sigma(k) = SW + sqrt(sum(((data-mi(k)).^2) .* p(:,k),1)/sum_pk);
            if sigma(k) <= 0; sigma = 1e-5; disp('Sigma too low. Too many components.');end
            f(:,k) = alpha(k) * normpdf(data,mi(k),sigma(k)); 
        end
        f(isnan(f)) = 1e-45; f(f==0) = 1e-45;
        L_old = L_new;
        it = it + 1;       
    end
    
    %calculating BIC
    bic = -2*L_old + (3*K - 1)*log(N);
    if bic == Inf || bic == 0
        disp('EM crash. Repeat calculations.')
        it = 1;
    end
end

[mi,ind] = sort(mi); 
alpha = alpha(ind); 
sigma = sigma(ind);

if K_noise >= K
    thr = max(data) + 1e-10;
elseif K_noise >= 0 && K_noise < K
    if K == 1
        K_noise = 0;
    elseif K == 2
        K_noise = 1;
        thr = find_thr(data,alpha,mi,sigma,[0;1],draw);
    elseif K_noise == 0
        temp3 = [alpha;mi;sigma]';
        idx = kmeans(temp3,2,'emptyaction','singleton','replicates',20);
        K_noise = sum(idx == idx(1));
        thr = find_thr(data,alpha,mi,sigma,idx-1,draw);
    elseif K > 2 && K_noise >0
        thr = find_thr(data,alpha,mi,sigma,[0,ones(1,K-1)],draw);
    end
elseif K_noise < 0
    thr = min(data) - 1e-10;     
end

if ~exist('thr','var')
    thr = NaN;
end

stats.thr = thr;
stats.alpha = alpha;
stats.mu = mi;
stats.K_noise = K_noise;
stats.K = K;
stats.sigma = sigma;
stats.logL = L_old;

%drawing results
if draw
    %histogram of input data
    [n,x] = hist(data,min(30,round(sqrt(N))));
    disp('After EM:')
    disp(num2str([alpha; mi; sigma]))
    disp(['Iterations: ' num2str(it) ' Stop crit: ' num2str(100*abs((L_new - L_old)/L_old))])
    
    f_temp = zeros(1e5,length(mi)); 
    x_temp = linspace(min(data),max(data),1e5)';
    step_x_temp = [min(diff(x_temp)); diff(x_temp)];
    for k=1:K; f_temp(:,k) = alpha(k)*normpdf(x_temp,mi(k),sigma(k)).*step_x_temp; end 
    
    figure; hold on; box on
    bar(x,n,[min(data) max(data)],'hist');
    plot(x_temp',mean(diff(x))*N*sum(f_temp./repmat(step_x_temp,1,K),2),'b.');
    colors(1:K_noise) = 'r'; colors(K_noise+1:K) = 'g'; 
    for a=1:K
        plot(x_temp',mean(diff(x))*N*f_temp(:,a)./step_x_temp,colors(a));
    end
    if ~isnan(thr); plot([thr,thr],ylim,'r--'); end
    title('After EM')
    lines = findobj(gca,'Type','Line');
    set(lines,'LineWidth',2)
    set(get(gca,'Ylabel'),'FontSize',14)
    xlabel('Variable'); ylabel('Counts');
end

function y = normpdf(x,mu,sigma)
y = exp(-0.5 * ((x - mu)./sigma).^2) ./ (sqrt(2*pi) .* sigma);

function thr = find_thr(data,alpha,mi,sigma,idx,draw)
% alpha,mi,sigma - components parameters
% idx index for informative/non-informative components 

idx = logical(idx);

%generate data with better precision
K = length(mi);
f_temp = zeros(1e7,K); 
x_temp = linspace(min(data),max(data),1e7)';
for k=1:K; f_temp(:,k) = alpha(k)*normpdf(x_temp,mi(k),sigma(k)); end 

%find GMM for informative and non-informative components
f1 = sum(f_temp(:,~idx),2);
f2 = sum(f_temp(:,idx),2);

%calculate difference of f1 and f2 and find its global minimum
f_diff = abs(f1-f2);
% [~,ind1] = max(f1);
% [~,ind2] = max(f2);
[~,ind] = sort(f_diff);

a = 1;
thr_ind = ind(a);
while thr_ind == 1 || thr_ind == length(ind)
    thr_ind = ind(a+1);
end
thr = x_temp(thr_ind);

if draw
    figure; subplot(2,1,1);
    plot(x_temp,f1,'g',x_temp,f2,'r')
    lines = findobj(gca,'Type','Line');
    set(lines,'LineWidth',2)
    set(get(gca,'Ylabel'),'FontSize',14)
    xlabel('Variable'); ylabel('Model')
    
    subplot(2,1,2);
    plot(x_temp,f_diff,'r')
    lines = findobj(gca,'Type','Line');
    set(lines,'LineWidth',2)
    set(get(gca,'Ylabel'),'FontSize',14)
    xlabel('Variable'); ylabel('Models difference')
end