clear all
load('/home/danielmessenger/Dropbox/Boulder/research/data/WSINDy_OL/wave2Du3c9/wave2DVCspeedu3_c0.9.mat_sweep_16-Aug-2022.mat')

T = length(Wscell{1});
lap_inds = find(ismember(tags_pde_G,{'u^{1}_{yy}','u^{1}_{xx}'}));
func = @(x) mean(x);
dat = cell2mat(cellfun(@(x) cellfun(@(y) func(y(lap_inds)) ,x)', Wscell, 'uni',0));
true_weights = cell2mat(cellfun(@(x) func(x(lap_inds)), axi((Kmem+1)/2:(Kmem+1)/2+tOL-1),'uni',0));

meandat=mean(dat,2);
maxdat=max(dat,[],2);
mindat=min(dat,[],2);

figure(1);clf

t=1:T;
a=area(t,[mindat maxdat-mindat]);
a(1).FaceColor = [1 1 1];
a(2).FaceColor = [0.5137    0.8627    0.9882];
hold on;
h3=plot(t,dat,'k-','linewidth',0.25);
h2=plot(t,true_weights,'k.','linewidth',2.5);
h1=plot(t,meandat,'r--','linewidth',2.5);
h4=plot(t,t*0,'k--')
% set(h1,'color',[ 0.8706   0 0.8118])
ylim([-1 3])
xlim([1 length(t)])

set(gca,'ticklabelinterpreter','latex','fontsize',14)
xlabel('iterations','interpreter','latex','fontsize',12)
ylabel('$c$','interpreter','latex','fontsize',12)
legend([h1;a(2);h2],{'mean $\hat{c}$','range $\hat{c}$','c'},'interpreter','latex','location','se','fontsize',12)
saveas(gcf,['~/Desktop/learnedc.png'])