function [t,maxdeg]=loadresults(dir,N,d,density,idx)
    filename=['P_',num2str(N),'_',num2str(d),'_',num2str(density),'_',num2str(idx),'.mat'];
    if exist([dir,'/',filename],'file')~=0
        data=load([dir,'/',filename]);
    else
        data.time=0;
        data.maxdeg=0;
    end
    t=data.time;
    maxdeg=data.maxdeg;
end