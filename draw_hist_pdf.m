function ok=draw_hist_pdf(data,mu_gmm,sig_gmm,pp_gmm)
% draw histogram versus pdf

hist_nb = max(10,min(100,sqrt(length(data))));
[hn,xx]=hist(data,round(hist_nb));
N=length(data);
shn=hn/(N*(xx(2)-xx(1)));
hold off
bar(xx,shn,1,'w');
hold on

KS=length(pp_gmm);
fit_pdf=0*data;
for kks=1:KS
    fit_pdf=fit_pdf+pp_gmm(kks)*normpdf(data,mu_gmm(kks),sig_gmm(kks));
    plot(data,pp_gmm(kks)*normpdf(data,mu_gmm(kks),sig_gmm(kks)),'b');
end

plot(data,fit_pdf,'r');
grid on; 
ylabel('PDF');
xlabel('X'); 