% load data from result files
setN=2:3:20;
setd=2:6;
setdensity=[0.2,0.4,0.6,0.8,1];
setidx=1:10;
T_ip=zeros(length(setN),length(setd),length(setdensity),length(setidx));
T_md=T_ip;
T_fmd=T_ip;
MDEG_ip=uint32(T_ip);
MDEG_md=MDEG_ip;
MDEG_fmd=MDEG_ip;
%ii=1;
for l=1:length(setidx)
    for k=1:length(setdensity)
        for j=1:length(setd)
            for i=1:length(setN)
                [T_ip(i,j,k,l),MDEG_ip(i,j,k,l)]=loadresults('result_ipdcsos',setN(i),setd(j),setdensity(k),setidx(l));
                [T_md(i,j,k,l),MDEG_md(i,j,k,l)]=loadresults('result_mddcsos',setN(i),setd(j),setdensity(k),setidx(l));
                [T_fmd(i,j,k,l),MDEG_fmd(i,j,k,l)]=loadresults('result_fmddcsos',setN(i),setd(j),setdensity(k),setidx(l));
                if MDEG_ip(i,j,k,l)==0
                    MDEG_ip(i,j,k,l)=2;
                    MDEG_md(i,j,k,l)=2;
                    MDEG_fmd(i,j,k,l)=2;
                end
                disp([i,j,k,l]);
            end
        end
    end
end

%%
% plot n v.s. time
X=setN;
Y=zeros(3,length(setN));
for i=1:length(setN)
    Y(1,i)=log(sum(T_ip(i,:,:,:),'all'));
    Y(2,i)=log(sum(T_md(i,:,:,:),'all'));
    Y(3,i)=log(sum(T_fmd(i,:,:,:),'all'));
end
figure(1);
plot(X,Y,'-o','LineWidth',1.5);
xticks(setN);
legend({'IP-DCSOS','MD-DCSOS','FMD-DCSOS'},'Location','northwest','FontSize',12);
xlabel('n','FontWeight','bold','FontSize',15);
ylabel('log time','FontWeight','bold','FontSize',15);
setupfig;
savefig('dcsos_nvstime');

%%
% plot d v.s. time
X=setd;
Y=zeros(3,length(setd));
for i=1:length(setd)
    Y(1,i)=log(sum(T_ip(:,i,:,:),'all'));
    Y(2,i)=log(sum(T_md(:,i,:,:),'all'));
    Y(3,i)=log(sum(T_fmd(:,i,:,:),'all'));
end
figure(2)
plot(X,Y,'-o','LineWidth',1.5)
xticks(setd);
legend({'IP-DCSOS','MD-DCSOS','FMD-DCSOS'},'Location','northwest','FontSize',12);
xlabel('d','FontWeight','bold','FontSize',15);
ylabel('log time','FontWeight','bold','FontSize',15);
setupfig;
savefig('dcsos_dvstime');

%%
% plot density v.s. time
X=setdensity;
Y=zeros(3,length(setdensity));
for i=1:length(setdensity)
    Y(1,i)=log(sum(T_ip(:,:,i,:),'all'));
    Y(2,i)=log(sum(T_md(:,:,i,:),'all'));
    Y(3,i)=log(sum(T_fmd(:,:,i,:),'all'));
end
figure(3)
plot(X,Y,'-o','LineWidth',1.5)
xticks(setdensity);
legend({'IP-DCSOS','MD-DCSOS','FMD-DCSOS'},'Location','northwest','FontSize',12);
xlabel('density','FontWeight','bold','FontSize',15);
ylabel('log time','FontWeight','bold','FontSize',15);
setupfig;
savefig('dcsos_densityvstime');

%%
% plot algo v.s. max degree
Y=zeros(numel(MDEG_ip),3);
Y(:,1)=MDEG_ip(:);
Y(:,2)=MDEG_md(:);
Y(:,3)=MDEG_fmd(:);
figure(4);
boxplot(Y,{'IP-DCSOS','MD-DCSOS','FMD-DCSOS'});
xlabel('DC-SOS decomposition algorithms','FontWeight','bold','FontSize',15);
ylabel('max degree of DC-SOS components','FontWeight','bold','FontSize',15);
setupfig;
box off;
savefig('dcsos_maxdeg');
