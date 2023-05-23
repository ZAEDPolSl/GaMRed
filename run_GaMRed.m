clc; clearvars; close all

%% Load data (single column .txt file with average or variance of expression;
%$ use log-transform for skewed distributions)
name = 'Stockholm_RNAseq';
data_table = readtable(['data/',name,'.txt'],'Delimiter','\t','HeaderLines', 0, 'ReadVariableNames', true);
data = data_table.Value;   % second column named 'Value' is an average gene expression
data = data(data<16.5);

%% Parameters
K = 10;     % maximum no. of Gaussian components
SW = 1e-2;  % minimum standard deviation of components (to avoid spikes, when you have multiple duplicated values)
ifshow = true;  % if plot results

%% run GaMRed K times
BIC = nan(K,1); % BIC for model fitting accuracy
thr = BIC;  % filtering threshold
stats = cell(K,1);  %GMM statistics
parfor a=1:K
    if ifshow;disp(['K=' num2str(a) '/' num2str(K)]);end
   [thr(a),BIC(a),stats{a}] = GaMRed(data,a,0,false,SW);
end

%% Minimize BIC to get the optimal no. of GMM components
[~,n_opt] = min(abs(BIC));
thr_opt = thr(n_opt);   % Final threshold for filtering genes
thr_opt = 8.03;

if ifshow
    disp(['Optimal threshold for ',num2str(n_opt),' components GMM model is ',num2str(thr_opt)])
    
    figure; hold on; box on;
    plot(1:K,BIC,'b*');
    plot([n_opt n_opt],get(gca,'Ylim'),'Color','red')
    xlabel('Number of Gaussian components');
    ylabel('BIC');
    title(['Dataset: ',name],'Interpreter','none')
    
    figure; draw_hist_pdf(sort(data),stats{n_opt}.mu,stats{n_opt}.sigma,stats{n_opt}.alpha)
    title([num2str(n_opt) ' components model'])
    plot([thr_opt,thr_opt],get(gca,'Ylim'),'r');
end

%% Filter data and save
del = data < thr_opt;     %remove features with value lower than the threshold
data_filt = data_table(~del,:);
disp([num2str(sum(del)),' (',num2str(round(100*sum(del)/length(del))), '%) features filtered. ', ...
    num2str(sum(~del)), ' features remained.'])
writetable(data_filt,[name,'_filt.txt'],'Delimiter','\t')