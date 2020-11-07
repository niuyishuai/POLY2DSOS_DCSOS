% script for testing all DSOS algorithms
%% generate polynomial
if 1
    clc;
    clear;
    N=5; % number of variables
    d=5; % degree
    polytype=0; %polynomial type 0: polylab, 1: multipoly, 2: yalmip, 3: syms, 4: sostools
    density=1; %[0,1]
    p=genpoly(N,d,polytype,density);
end
%% Test DSOS_MBS decompositions
method='PP';
isparal=true;
tic
[pos,neg]=poly2dsos(p,method,[],isparal);
%[pos,neg]=dsos_pp(p,isparal);
t=toc;

if isparal
    pool=gcp('nocreate');
    fprintf('PP_DSOS decomposition using %d workers in %f (s.)\n',pool.NumWorkers,t);
else    
    fprintf('PP_DSOS decomposition in %f (s.)\n',t);
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

%% Test DSOS_MBS decompositions
method='IP';
isparal=true;
tic
[pos,neg]=poly2dsos(p,method,[],isparal);
%[pos,neg]=dsos_ip(p,isparal);
t=toc;

if isparal
    pool=gcp('nocreate');
    fprintf('IP_DSOS decomposition using %d workers in %f (s.)\n',pool.NumWorkers,t);
else    
    fprintf('IP_DSOS decomposition in %f (s.)\n',t);
end
if 1
    tic
    pos=pos.simplify;
    toc
    tic
    neg=neg.simplify;
    toc
    degree_info=[pos.degree, neg.degree];
    sdisp(simplify(sum(pos-neg)-p));
end

%% Test DSOS_DBS decompositions
method='DBS';
ismatrixform=true; % only for MBS|DBS
tic
[pos,neg]=poly2dsos(p,method,ismatrixform);
%[pos,neg]=dsos_dbs(p,ismatrixform);
t=toc;

fprintf('DBS_DSOS decomposition in %f (s.)\n',t);
if 1
    fprintf('verify results:\n');
    if ~ismatrixform
        tic
        pos=pos.simplify;
        toc
        tic
        neg=neg.simplify;
        toc
        degree_info=[pos.degree, neg.degree];
        sdisp(simplify(sum(pos-neg)-p));
    else
        tic
        [pos,neg]=dbsmatrixform2poly(pos,neg);
        toc
        degree_info=[pos.degree, neg.degree];
        sdisp(simplify(pos-neg-p));
    end
end

%% Test DSOS_MBS decompositions
method='MBS';
ismatrixform=true; % only for MBS|DBS
tic
[pos,neg,dsos]=poly2dsos(p,method,ismatrixform);
%[pos,neg,dsos]=dsos_mbs(p,ismatrixform);
t=toc;

fprintf('MBS_DSOS decomposition in %f (s.)\n',t);
if 1
    fprintf('verify results:\n');
    if ~ismatrixform
        tic
        pos=pos.simplify;
        toc
        tic
        neg=neg.simplify;
        toc
        degree_info=[pos.degree, neg.degree];
        sdisp(simplify(sum(pos-neg)-p));
    else
        tic
        [pos,neg]=mbsmatrixform2poly(dsos);
        toc
        degree_info=[pos.degree, neg.degree];
        sdisp(simplify(pos-neg-p));
    end
end
