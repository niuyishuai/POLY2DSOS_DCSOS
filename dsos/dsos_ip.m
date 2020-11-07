function [PSOS,NSOS]=dsos_ip(p,isparal)
    % Parity DSOS decomposition for polynomial (MPOLY version)
    % Syntax:
    %   [PSOS,NSOS]=DSOS_IP(p,isparal)
    %
    % Inputs:
    % p: MPOLY polynomial
    % isparal: true for parallel mode (false by default)
    %
    % Outputs:
    % PSOS: sos polynomial (MPOLY)
    % NSOS: sos polynomial (MPOLY)
    %
    if nargin<2
        isparal=false;
    end
    P=MPOLY.zeros(p.n,p.k,1);
    N=MPOLY.zeros(p.n,p.k,1);
    % Extract coefs and monomials of f
    [c,m]=coefficients(p);
    x=MPOLY.mpolyvars(p.n);
    if isparal
        parfor i=1:length(c)
            % for each monomial, compute its DSOS components
            [P(i),N(i)]=parity_mono(c(i),m(i),x);
        end
    else
        for i=1:length(c)
            % for each monomial, compute its DSOS components
            [P(i),N(i)]=parity_mono(c(i),m(i),x);
        end
    end
    PSOS=sum(P);
    NSOS=sum(N);
end

function [PSOS,NSOS]=parity_mono(coeff,mono,x)
    % Parity DSOS decomposition for a single monomial
    
    % STEP 2: extract o(x) and e(x)^2 for mono
    % get odd and even degree variables
    deg=mono.pow;
    modd=mod(deg,2);
    idxo=find(modd);% find index of odd degree variables, i.e., O(m)
    edeg=deg-modd; % degrees of e^2
    % produce e2 = c*e^2, where c=coeff if coeff>0 and c=-coeff if coeff<0
    e2=mono;
    e2.pow=edeg;
    if coeff>0
        e2.coef=coeff;
    else
        e2.coef=-coeff;
    end
    
    % STEP 3: Procedure D for DSOS decomposition of o(x)
    [Psos,Nsos]=procedureD(idxo,mono.n,x);
    
    % STEP 4: Get DSOS decomposition for mono
    if  coeff>0
        PSOS=e2*Psos;
        NSOS=e2*Nsos;
    else
        PSOS=e2*Nsos;
        NSOS=e2*Psos;
    end
end

function [Psos,Nsos]=procedureD(idxo,n,x)
    % procedure D for DSOS decomposition of o(x)
    if isempty(idxo) % if index set O(m) is empty
        Psos=1;
        Nsos=0;
        return;
    else
        oddx=x(idxo);% extract odd degree variables
        if mod(length(idxo),2)>0
            oddx=[oddx;MPOLY.ones(n)]; % add 1 if O(m) is odd
        end
        SOS=PROC_D(oddx,n);
        len=length(oddx)/2;
        A=ff2n(len);
        B=[A,~A];
        Bp=B(mod(sum(~A,2),2)==0,:);
        Bn=B(mod(sum(~A,2),2)~=0,:);
        rr=size(B,1)/2;
        SOSB=repmat(SOS,1,rr);
        SOS1=reshape(SOSB(Bp'>0),len,rr);
        SOS2=reshape(SOSB(Bn'>0),len,rr);
        Psos=sum(SOS1.prod);
        Nsos=sum(SOS2.prod);
    end
end

function SOS=PROC_D(oddx,nn)
    len=length(oddx);
    id1=1:len;
    id2=randperm(len,len/2);
    id1(id2)=[];
    method=1; % change method for speeding polynomial computations
    switch method
        case 1
            % fast computation
            SOS1=MPOLY.zeros(nn,len/2,1);
            SOS2=SOS1;
            for i=1:len/2
                % SOS1 = (xi+xj)^2/2
                SOS1(i).pow=[2*oddx(id1(i)).pow;2*oddx(id2(i)).pow;...
                    oddx(id1(i)).pow+oddx(id2(i)).pow];
                SOS1(i).coef=[1;1;2]/4;
                SOS1(i).k=3;
                % SOS2 = (xi-xj)^2/2
                SOS2(i).pow=SOS1(i).pow;
                SOS2(i).coef=[1;1;-2]/4;
                SOS2(i).k=3;
            end
            SOS=[SOS1;SOS2];
        case 2
            % using polynomial operations
            SOS=[oddx(id1)+oddx(id2);oddx(id1)-oddx(id2)].^2/4;
    end
end
