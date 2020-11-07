% Parallel tests on DSOS algorithms
inputdir='dataset';
outputdir='result_ppdsos';
mkdir(outputdir);
count=0;
setN=2;%2:3:20;
setd=2;%2:6;
setdensity=0.2;%[0.2,0.4,0.6,0.8,1];
nbproblems=10*length(setN)*length(setd)*length(setdensity);
% parameter settings
method = 'PP'; % 'PP'|'IP'|'DBS'|'MBS'
ismatrixform=false;
for N=setN
    for density=setdensity
        for d=setd
            %% Call PP-DSOS
            parfor i=1:10
                % loading random polynomial
                filename=['P_',num2str(N),'_',num2str(d),'_',num2str(density),'_',num2str(i),'.mat'];
                data=load([inputdir,'/',filename]);
                p=data.p;
                fprintf('* load polynomial of %d variables and of degree %d density %.2f from %s ... %d / %d\n',N,d,density,filename,count+i,nbproblems);
                
                tic
                [PSOS,~] = poly2dsos(p,method,ismatrixform);
                t=toc;
                
                fprintf('PP-DSOS within %.3f seconds.\n',t);
                parsave([outputdir,'/',filename],p,N,d,density,PSOS,[],t);
            end
            count=count+10;
        end
    end
end