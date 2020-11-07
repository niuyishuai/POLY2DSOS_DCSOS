% script for testing all DCSOS algorithms
%% generate polynomial
if 1
    clc;
    clear;
    N=6; % number of variables
    d=6; % degree
    polytype=0; %polynomial type: 0: polylab, 1: multipoly, 2: yalmip, 3: syms, 4: sostools
    density=1; %[0,1]
    p=genpoly(N,d,polytype,density);
end
%% Test DCSOS decomposition
isparal=true;
tic
%method='IP'; % 'IP|MD|FMD'
%[pos,neg]=poly2dcsos(p,method,isparal);
[pos,neg]=dcsos_ip(p,isparal);
%[pos,neg]=dcsos_md(p,isparal);
%[pos,neg]=dcsos_fmd(p,isparal);
t=toc;
if isparal
    pool=gcp('nocreate');
    fprintf('DCSOS decomposition using %d workers in %f (s.)\n',pool.NumWorkers,t);
else    
    fprintf('DCSOS decomposition in %f (s.)\n',t);
end

if 1
    fprintf('verify results:\n');
    tic
    pos=pos.simplify;
    toc
    tic
    neg=neg.simplify;
    toc
    degree_info=[pos.degree, neg.degree];
    sdisp(simplify(sum(pos-neg)-p));
end
