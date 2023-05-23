function [res,res2,res3] = filter_type(data,grp,draw)
% FILTER_TYPE(DATA,GRP)
% Choose filter type based on LR slope coefficient.
% Input:
% data - data matrix [mxn]
% grp - group labels [1xn]
% Output:
% res,res2 - chosen filter

%calculate values for different filter type
S = mean(data,2);
V = log(var(pow2(data),0,2));
LV = log(var(data,0,2));

%perform statistical testing
grp_un = unique(grp);
grp_nb = length(grp_un);
fea_nb = length(S);
samp_nb = length(grp);
stat = zeros(samp_nb,1);
if grp_nb==2    %test T
    g1 = grp==grp_un(1); 
    g2 = grp==grp_un(2);
    for a=1:fea_nb
        stat(a) = abs(mean(data(a,g1))-mean(data(a,g2)))/...
            sqrt(var(data(a,g1))/sum(g1) + var(data(a,g2))/sum(g2));
    end
%     for a=1:fea_nb
%         [~,~,stats] = ranksum(data(a,g1),data(a,g2));
%         stat(a) = abs(stats.zval);
%     end
    
elseif grp_nb>2 %test ANOVA
    for a=1:fea_nb
        tot_av = mean(data(a,:));
        SS_among = 0; SS_within = 0;
        for b=1:grp_nb
            SS_among = SS_among + sum(grp==grp_un(b))*(mean(data(a,grp==grp_un(b)))-tot_av)^2;
            SS_within = SS_within + sum((data(a,grp==grp_un(b))-mean(data(a,grp==grp_un(b)))).^2);
        end
        stat(a) = (SS_among/(grp_nb-1))/(SS_within/(samp_nb-grp_nb));
    end
else
    error('Wrong number of groups.')
end
stat(isinf(stat)) = 0;
stat(isnan(stat)) = 0;

%find sorting indices for each filter type
[~,ind_S] = sort(S);
[~,ind_V] = sort(V);
[~,ind_LV] = sort(LV);

% calc linear regression coefficient b1
X = (1:fea_nb)';   %vector with ranks
% s = std(stat)/std(X);
b_S = corr(X,stat(ind_S));%*s;
b_V = corr(X,stat(ind_V));%*s;
b_LV = corr(X,stat(ind_LV));%*s;
% w kodzie mo¿esz pomin¹æ zmienn¹ s, bo dla ka¿dego typu filtra jest to ta
% sama wartoœæ, ale zostawi³em ¿eby pokazaæ ideê

[~,res] = max([b_S,b_V,b_LV]);
switch res
    case 1
        res2 = 'S';
    case 2
        res2 = 'V';
    case 3
        res2 = 'LV';
end
res3 = max([b_S,b_V,b_LV],0);
res3 = 100*res3/sum(res3); 

if draw
    figure; subplot(3,1,1); hold on; box on;
    plot(X,stat(ind_S),'.');h = lsline; set(h(1),'color','r'); xlim([1,fea_nb])
    xlabel('Rank of S'); ylabel('Statistic')
    subplot(3,1,2); hold on; box on;
    plot(X,stat(ind_V),'.');h = lsline; set(h(1),'color','r'); xlim([1,fea_nb])
    xlabel('Rank of V'); ylabel('Statistic')
    subplot(3,1,3); hold on; box on;
    plot(X,stat(ind_LV),'.');h = lsline; set(h(1),'color','r'); xlim([1,fea_nb])
    xlabel('Rank of LV'); ylabel('Statistic')
    
    disp(['The best filter type is ' res2])
end