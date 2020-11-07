% load data from result files
%outputdir='result_ppdsos';
setN=2:3:20;
setd=2:6;
setdensity=[0.2,0.4,0.6,0.8,1];
setidx=1:10;
T_pp=zeros(length(setN),length(setd),length(setdensity),length(setidx));
T_ip=T_pp;
T_dbs=T_pp;
T_mbs=T_pp;
MDEG_pp=uint32(T_pp);
MDEG_ip=MDEG_pp;
MDEG_dbs=MDEG_pp;
MDEG_mbs=MDEG_pp;
%ii=1;
for l=1:length(setidx)
    for k=1:length(setdensity)
        for j=1:length(setd)
            for i=1:length(setN)
                %    if ii==876
                %        disp([i j k l]);
                %    end
                %    ii=ii+1;
                [T_pp(i,j,k,l),MDEG_pp(i,j,k,l)]=loadresults('result_ppdsos',setN(i),setd(j),setdensity(k),setidx(l));
                [T_ip(i,j,k,l),MDEG_ip(i,j,k,l)]=loadresults('result_ipdsos',setN(i),setd(j),setdensity(k),setidx(l));
                [T_dbs(i,j,k,l),MDEG_dbs(i,j,k,l)]=loadresults('result_dbsdsos',setN(i),setd(j),setdensity(k),setidx(l));
                [T_mbs(i,j,k,l),MDEG_mbs(i,j,k,l)]=loadresults('result_mbsdsos',setN(i),setd(j),setdensity(k),setidx(l));
                disp([i,j,k,l]);
            end
        end
    end
end

%%
% plot n v.s. time
X=setN;
Y=zeros(4,length(setN));
for i=1:length(setN)
    Y(1,i)=log(sum(T_pp(i,:,:,:),'all'));
    Y(2,i)=log(sum(T_ip(i,:,:,:),'all'));
    Y(3,i)=log(sum(T_dbs(i,:,:,:),'all'));
    Y(4,i)=log(sum(T_mbs(i,:,:,:),'all'));
end
figure(1);
setupfig;
plot(X,Y,'-o','LineWidth',1.5);
xticks(setN);
legend({'PP-DSOS','IP-DSOS','DBS-DSOS','MBS-DSOS'},'Location','northwest');
xlabel('n');
ylabel('log time');
savefig('nvstime');

%%
% plot d v.s. time
X=setd;
Y=zeros(4,length(setd));
for i=1:length(setd)
    Y(1,i)=log(sum(T_pp(:,i,:,:),'all'));
    Y(2,i)=log(sum(T_ip(:,i,:,:),'all'));
    Y(3,i)=log(sum(T_dbs(:,i,:,:),'all'));
    Y(4,i)=log(sum(T_mbs(:,i,:,:),'all'));
end
figure(2)
setupfig;
plot(X,Y,'-o','LineWidth',1.5)
xticks(setd);
legend({'PP-DSOS','IP-DSOS','DBS-DSOS','MBS-DSOS'},'Location','northwest');
xlabel('d');
ylabel('log time');
savefig('dvstime');

%%
% plot density v.s. time
X=setdensity;
Y=zeros(4,length(setdensity));
for i=1:length(setdensity)
    Y(1,i)=log(sum(T_pp(:,:,i,:),'all'));
    Y(2,i)=log(sum(T_ip(:,:,i,:),'all'));
    Y(3,i)=log(sum(T_dbs(:,:,i,:),'all'));
    Y(4,i)=log(sum(T_mbs(:,:,i,:),'all'));
end
figure(3)
setupfig;
plot(X,Y,'-o','LineWidth',1.5)
xticks(setdensity);
legend({'PP-DSOS','IP-DSOS','DBS-DSOS','MBS-DSOS'},'Location','northwest');
xlabel('density');
ylabel('log time');
savefig('densityvstime');

%%
% plot algo v.s. max degree
Y=zeros(numel(MDEG_pp),4);
Y(:,1)=MDEG_pp(:);
Y(:,2)=MDEG_ip(:);
Y(:,3)=MDEG_dbs(:);
Y(:,4)=MDEG_mbs(:);
figure(4);
setupfig;
boxplot(Y,{'PP-DSOS','IP-DSOS','DBS-DSOS','MBS-DSOS'});
xlabel('D-SOS decomposition algorithms');
ylabel('max degree of D-SOS components');
savefig('maxdeg');
