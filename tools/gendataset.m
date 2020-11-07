datadir = 'dataset';
mkdir(datadir);
polytype=0; %0: polylab, 1: multipoly, 2: yalmip, 3: syms, 4: sostools
setN=2:3:20;
setd=2:6;
setdensity=[0.2,0.4,0.6,0.8,1];
fprintf('Start generating random polynomial dataset\n');
fprintf('* N = [%d, %d], d = [%d, %d], density = [%.1f, %.1f]\n', min(setN), max(setN), min(setd), max(setd), min(setdensity),max(setdensity));
nbproblems=10*length(setN)*length(setd)*length(setdensity);
iter=1;
for N=setN
    for density=setdensity
        for d=setd
            for i=1:10
                p=genpoly(N,d,polytype,density);
                filename=['P_',num2str(N),'_',num2str(d),'_',num2str(density),'_',num2str(i),'.mat'];
                save([datadir,'\\',filename],'p','N','d','density');
                fprintf('* %s is generated... %d / %d \n',filename, iter, nbproblems);
                iter = iter+1;
            end
        end
    end
end
