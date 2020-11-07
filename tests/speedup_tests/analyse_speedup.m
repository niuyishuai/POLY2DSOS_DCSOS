%% draw speedup
prefix='dcsos_ip';
resultfile='_speedup_1-60.mat';
load([prefix,resultfile]);
figure(1)
hold all;
plot(nbcpus,timelst(:,1)./timelst,'-o','LineWidth',1.5);
xlabel('number of processors','FontWeight','bold','FontSize',12);
ylabel('speedup','FontWeight','bold','FontSize',12);
setupfig(1);
grid on;
labels={'n=6,d=3','n=8,d=4','n=10,d=6','n=12,d=6','n=16,d=6'};
legend(labels,'FontSize',12,'Location','northwest');
%legend(labels,'FontSize',12,'Location','best');
savefig([prefix,'_speedup'])

%%
figure(2)
hold all;
plot(nbcpus,(timelst(:,1)./timelst)./nbcpus,'-o','LineWidth',1.5);
xlabel('number of processors','FontWeight','bold','FontSize',12);
ylabel('efficiency','FontWeight','bold','FontSize',12);
setupfig(1);
grid on;
legend(labels,'FontSize',12);
savefig([prefix,'_efficiency'])
