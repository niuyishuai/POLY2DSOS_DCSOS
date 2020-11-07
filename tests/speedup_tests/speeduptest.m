function speeduptest(fh,prefix)
N=[6,8,10,12,16];
d=[3,4,6,6,6];
polytype=0; %0: polylab, 1: multipoly, 2: yalmip, 3: syms, 4: sostools
density=1;
nbcpus=0:5:60; % testing up to 60 cpus, make sure you have enough cpu, otherwise, change it to fit your case
nbcpus(1)=1;
timelst=zeros(numel(N),numel(nbcpus)); % each row are excution time for all processors

p=MPOLY.zeros(1,1,numel(N));
% generate polynomial
for idx=1:length(N)
    p(idx)=genpoly(N(idx),d(idx),polytype,density);
end

% for a setting of processors
for i=length(nbcpus):-1:1
    if nbcpus(i)==1
        isparal=false;
        kmax=1;
    else
        isparal=true;
        kmax=3; % number of repeatations, for average solution time
        parpool(nbcpus(i));
    end
    % test for all polynomials
    for idx=1:length(N)
        tic
        for k=1:kmax
            fh(p(idx),isparal); % speed only
        end
        timelst(idx,i)=toc;
        timelst(idx,i)=timelst(idx,i)/kmax; % compute average time
        fprintf('%s decomposition of a polynomial with %d variables and %d degree using %d CPUs in %f (s.) -> %s\n',...
            prefix,N(idx),d(idx),nbcpus(i),timelst(idx,i),datetime);
    end
    delete(gcp);
end
save([prefix,'_speedup_1-60.mat'],'p','timelst','nbcpus');
end
