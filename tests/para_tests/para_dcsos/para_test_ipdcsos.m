% Parallel tests on DSOS algorithms
inputdir='dataset';
outputdir='result_ipdcsos';
mkdir(outputdir);
count=0;
setN=2:3:20;
setd=2:6;
setdensity=[0.2,0.4,0.6,0.8,1];
nbproblems=10*length(setN)*length(setd)*length(setdensity);
% parameter settings
method = 'IP'; % 'IP'|'MD'|'FMD'
isparal = true;
for N=setN
    for density=setdensity
        for d=setd
            %% Call IP-DCSOS
            for i=1:10
                % loading random polynomial
                filename=['P_',num2str(N),'_',num2str(d),'_',num2str(density),'_',num2str(i),'.mat'];
                data=load([inputdir,'/',filename]);
                p=data.p;
                fprintf('* load polynomial of %d variables and of degree %d density %.2f from %s ... %d / %d\n',N,d,density,filename,count+i,nbproblems);
                
                tic
                [PSOS,~] = poly2dcsos(p,method,isparal);
                t=toc;
                
                fprintf('IP-DCSOS within %.3f seconds -> %s.\n',t,datetime);
                parsave([outputdir,'/',filename],p,N,d,density,PSOS,[],t);
            end
            count=count+10;
        end
    end
end
