function [pc,nc]=dcsos_md(p,isparal)
    % minimal degree DCSOS decomposition (MPOLY)
    % Syntax:
    %   [g,h]=dcsos_md(f)
    %
    % Inputs:
    % f: MPOLY polynomial
    %
    % Outputs:
    % g: csos polynomial(MPOLY)
    % h: csos polynomial(MPOLY)
    %global polyzero;
    %global polyone;
    if nargin<2
        isparal=false;
    end
    polyzero = MPOLY.zeros(p.n,1,1);
    pc=repmat(polyzero,p.k,1);
    nc=pc;
    % STEP 1: extract monomials
    [c,m]=coefficients(p);
    x=MPOLY.mpolyvars(p.n);
    % for each monomial, do DCSOS decomposition
    if isparal
        parfor i=1:length(c)
            [pc(i),nc(i)]=dcsos_mono(c(i),m(i),x);
        end
    else
        for i=1:length(c)
            [pc(i),nc(i)]=dcsos_mono(c(i),m(i),x);
        end
    end
    %g=simplify(sum(pcmat));
    %h=simplify(sum(ncmat));
    % sum the results for DCSOS of polynomial
 %   g=sum(pc);
 %   h=sum(nc);
end

function [pcsos,ncsos]=dcsos_mono(coeff,mono,x)
    % STEP 2: Procedure F
    [listo,liste]=ProcedureF(mono,x);
    % STEP 3: Procedure M
    if isempty(listo) && isempty(liste)
        % only for constant
        pcsos=mono;
        ncsos=mono;
        if coeff>0
            pcsos.coef=coeff;
            ncsos.coef=0;
        else
            pcsos.coef=0;
            ncsos.coef=-coeff;
        end
        return;
    elseif ~isempty(listo)
        % O non-empty
        len=length(listo)/2;
        A=ff2n(len);
        B=[A,~A];
        idx=mod(sum(~A,2),2)==0;
        Bp=B(idx,:);% extract the positive product list
        Bn=B(~idx,:);% extract the negative product list
        rr=size(Bp,1); %size(B,1)/2;
        oddB=repmat(listo,1,rr);
        oddp=reshape(oddB(Bp'>0),len,rr);% extract the nonzero positive sos list of o(x)
        oddn=reshape(oddB(Bn'>0),len,rr);% extract the nonzero negative sos list of o(x)
        if ~isempty(liste)
            even=repmat(liste,1,rr);
            oddp=[oddp;even];
            oddn=[oddn;even];
        end
        m=len+length(liste);
        [pm,nm]=getsubsets(m);
        % loops       
        [pc,nc]=powsum2(oddp,oddn,pm,nm,m);
        %pc=sum(sum((pm*oddp).^m+(nm*oddn).^m))/factorial(m);
        %nc=sum(sum((pm*oddn).^m+(nm*oddp).^m))/factorial(m);
    else
        % O empty
        m=length(liste);
        [pm,nm]=getsubsets(m);
        [pc,nc]=powsum1(liste,pm,nm,m);
        %pc=sum((pm*liste).^m)/factorial(m);
        %nc=sum((nm*liste).^m)/factorial(m);
    end
    if coeff>0
        pcsos=coeff*pc;
        ncsos=coeff*nc;
    else
        pcsos=(-coeff)*nc;
        ncsos=(-coeff)*pc;
    end
end

%%
function [listo,liste]=ProcedureF(mono,x)
    %global polyone
    polyone = MPOLY.ones(mono.n,1,1);

    listo=[];
    liste=[];
    pow=mono.pow';%extract degree vector of all variables
   
    %extract o(x) and e(x)^2
    % get odd degree variables and even degree variables
    oddpow=mod(pow,2);
    idxo=find(oddpow);% find index of odd degree variables
    evenpow=pow-oddpow; % power for e^2
    idxe=find(evenpow);% find index of even degree variables
    % get couple list of odd degree variables as listo 
    if ~isempty(idxo)
        oddx=x(idxo);
        if mod(length(idxo),2)~=0
            oddx=[oddx;polyone]; % with additional 1 if the length is odd
        end
        listo=getlisto(oddx);
    end
    % get couple list of even degree variables as liste
    if ~isempty(idxe)
        evenx=x(idxe);% extract the even variables
        evendeg=evenpow(idxe); %extract the even degree;
        liste=getliste(evenx,evendeg);
    end
end

%%
function Listo=getlisto(oddx)
    %global polyzero
    polyzero = MPOLY.zeros(oddx(1).n,1,1);
    
    len=length(oddx);
    id1=1:len;
    id2=randperm(len,len/2);
    id1(id2)=[];
    SOS1=repmat(polyzero,len/2,1);
    SOS2=SOS1;
    for i=1:len/2
        SOS1(i).pow=[2*oddx(id1(i)).pow;2*oddx(id2(i)).pow;...
            oddx(id1(i)).pow+oddx(id2(i)).pow];
        SOS1(i).coef=[1;1;2]/4;
        SOS1(i).k=3;
        SOS2(i).pow=SOS1(i).pow;
        SOS2(i).coef=[1;1;-2]/4;
        SOS2(i).k=3;
    end
    Listo=[SOS1;SOS2];
end

%%
function Liste=getliste(evenx,evendeg)
    %global polyzero
    polyzero = MPOLY.zeros(evenx(1).n,1,1);

    k=1;
    Liste=repmat(polyzero,sum(evendeg)/2,1);   % only g part in e^2
    for i=1:length(evenx)
        Liste(k:k+(evendeg(i)/2)-1)=evenx(i)^2;
        k=k+evendeg(i)/2;
    end
end

%%
function [pm,nm]=getsubsets(m)
    A=ff2n(m);%generate all subsets
    a=sum(A,2);%divide subsets
    if mod(m,2)==0
        pm=A(mod(a,2)==0,:);%get the subsets to compute pcsos
        nm=A(mod(a,2)~=0,:);%get the subsets to compute ncsos
    else
        pm=A(mod(a,2)~=0,:);
        nm=A(mod(a,2)==0,:);
    end
end

%%
function [pc,nc]=powsum2(oddp,oddn,pm,nm,m)
    %        pc=sum(sum((pm*oddp).^m+(nm*oddn).^m))/factorial(m);
    %        nc=sum(sum((pm*oddn).^m+(nm*oddp).^m))/factorial(m);
    % loops
    %global polyzero
    polyzero = MPOLY.zeros(oddp(1).n,1,1);
    
    pc=polyzero;
    nc=polyzero;
    for i=1:size(pm,1) 
        idx=logical(pm(i,:));
        if ~any(idx)
            a=polyzero;
            a1=polyzero;
        elseif sum(idx)==1
            a=sum(oddp(idx,:).^m);
            a1=sum(oddn(idx,:).^m);
        else
            a=sum(sum(oddp(idx,:)).^m);
            a1=sum(sum(oddn(idx,:)).^m);            
        end
        idx=logical(nm(i,:));
        if ~any(idx)
            b=polyzero;
            b1=polyzero;
        elseif sum(idx)==1
            b=sum(oddn(idx,:).^m);
            b1=sum(oddp(idx,:).^m);
        else    
            b=sum(sum(oddn(idx,:)).^m);
            b1=sum(sum(oddp(idx,:)).^m);           
        end
        pc=pc+a+b;
        nc=nc+a1+b1;
    end
    facm=factorial(m);
    pc=pc/facm;
    nc=nc/facm;
end

%%
function [pc,nc]=powsum1(list,pm,nm,m)
    %pc=sum((pm*liste).^m)/factorial(m);
    %nc=sum((nm*liste).^m)/factorial(m);
    % loops
    %global polyzero
    polyzero = MPOLY.zeros(list(1).n,1,1);
    
    pc=polyzero;
    nc=polyzero;
    for i=1:size(pm,1)
        idx=logical(pm(i,:));
        if ~any(idx)
            a=polyzero;
        else
            a=sum(list(idx))^m;
        end
        idx=logical(nm(i,:));
        if ~any(idx)
            b=polyzero;
        else
            b=sum(list(idx))^m;
        end
        pc=pc+a;
        nc=nc+b;
    end
    facm=factorial(m);
    pc=pc/facm;
    nc=nc/facm;
end

