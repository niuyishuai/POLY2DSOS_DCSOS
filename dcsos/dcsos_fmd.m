function [pc,nc]=dcsos_fmd(p,isparal)
    % Formulation based minimal degree DCSOS decomposition
    % Syntax:
    %   [g,h]=dcsos_fmd(p,isparal)
    % Inputs:
    % p: mpoly polynomial
    % isparal: proceed monomials in parallel
    % Outputs:
    % g: csos polynomial
    % h: csos polynomial
    
    if nargin<2
        isparal=false;
    end
    pc=MPOLY.zeros(p.n,p.k,1);
    nc=pc;
    [c,m]=coefficients(p);    
    if isparal
        parfor i=1:length(c)
            % decomposition for monomial in parallel
            [pc(i),nc(i)]=dcsos_mono(c(i),m(i));
        end
    else
        for i=1:length(c)
            % decomposition for monomial
            [pc(i),nc(i)]=dcsos_mono(c(i),m(i));
        end
    end
%    g=pc; 
%    h=nc;
%    g=sum(pc);
%    h=sum(nc);
end

%%
function [pcsos,ncsos]=dcsos_mono(coeff,mono)
    nn=mono.n;
    pow=mono.pow;
    idx=find(pow);
    x=MPOLY.mpolyvars(nn);
    if isempty(idx)
        pctmp=MPOLY.ones(nn,1,1);
        nctmp=MPOLY.zeros(nn,1,1);
        if coeff>0
            pcsos=coeff*pctmp;
            ncsos=coeff*nctmp;
        else
            pcsos=-coeff*nctmp;
            ncsos=-coeff*pctmp;
        end
    else
        %fomulate monomial as x1*x1*x1....xn*xn*xn
        V=[];
        for k=1:length(idx)
            V=[V;MPOLY.ones(nn,pow(idx(k)),1)*x(idx(k))];
        end
        %add 1 for odd length
        if mod(length(V),2)~=0
            V=[V;MPOLY.ones(nn,1,1)];
        end
        
        m=length(V); % length of V
        Aset=ff2n(m);
        a=sum(Aset,2);
        
        Vecx=powsum(Aset,V,m);
        %Vecx=(Aset*V).^m;
        %Vecx=simplify(Aset*V).^m;
        pidx=mod(a,2)==0;
        cc=coeff/factorial(m);
        if coeff>0
            pcsos=cc*sum(Vecx(pidx));%produce pcsos without coef
            ncsos=cc*sum(Vecx(~pidx));%produce ncsos without coef
        else
            pcsos=(-cc)*sum(Vecx(~pidx));%produce pcsos without coef
            ncsos=(-cc)*sum(Vecx(pidx));%produce ncsos without coef
        end
        pcsos=pcsos.simplify; % to avoid out of memory, but slow down a little
        ncsos=ncsos.simplify;
    end
%    sdisp(simplify(pcsos-ncsos-mono*coeff))
end

function R=powsum(C,V,m)
    %R=(C*V).^m;
    % loops
    polyzero = MPOLY.zeros(V(1).n,1,1);
    R = MPOLY.zeros(V(1).n,size(C,1),1);
    
    for i=1:size(C,1)
        idx=logical(C(i,:));
        if ~any(idx)
            R(i)=polyzero;
        else
            R(i)=simplify(sum(V(idx)))^m; % simplify for faster power
            %s1=sum(V(idx));
            %s2=simplify(s1);
            %R(i)=s2^m;
            % compute s2^m
            %[~,~,B]=MPOLY.transmatconvex(s2.k,m);
            %bincoef=floor(exp(gammaln(m+1) - sum(gammaln(B+1),2)) + .5);
            %R(i)=MPOLY(V(1).n, bincoef.*prod(s2.coef'.^B,2), B*s2.pow);
            
        end
    end
end

