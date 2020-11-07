function [pcsos,ncsos]=dcsos_ip(p,isparal)
    if nargin<2
        isparal=false;
    end
    pcsos=MPOLY.zeros(p.n,p.k,1);
    ncsos=pcsos;
    %extract the coefficients and monomials
    [c,m]=coefficients(p);
    x=MPOLY.mpolyvars(p.n);
    if isparal
        parfor i=1:length(c)
            % do decomposition for monomial
            [pcsos(i),ncsos(i)]=dcsos_mono(c(i),m(i),x);
        end
    else
        for i=1:length(c)
            % do decomposition for monomial
            [pcsos(i),ncsos(i)]=dcsos_mono(c(i),m(i),x);
        end
    end
    % g=simplify(sum(pcsos));
    % h=simplify(sum(ncsos));
    %g=sum(pcsos);
    %h=sum(ncsos);
end
%%
function [pcsos,ncsos]=dcsos_mono(coeff,mono,x)
    [List,deg]=ProcedureF(mono,x);
    List=ProcedureM(List,deg,1);
    if coeff>0
        pcsos=List(1,1)*coeff;
        ncsos=List(1,2)*coeff;
    else
        pcsos=List(1,2)*(-coeff);
        ncsos=List(1,1)*(-coeff);
    end
end
%%
function [List,deg]=ProcedureF(mono,x)
    nn=mono.n;
    pow=mono.pow';%extract degree vector of all variables
    
    %extract o(x) and e(x)^2
    %get odd degree variables and even degree variables
    oddpow=mod(pow,2);
    idxo=find(oddpow);% find index of odd degree variables
    evenpow=pow-oddpow;
    idxe=find(evenpow);% find index of even degree variables
    
    if isempty(idxo) && isempty(idxe)
        List=[MPOLY.ones(nn),MPOLY.zeros(nn)];
        deg=[0,0];
        return;
    else
        List=[];
        deg=[];
        if ~isempty(idxo)
            oddx=x(idxo);
            if mod(length(idxo),2)~=0 %produce even length of Vector
                oddx=[oddx;MPOLY.ones(nn)];
            end
            List=getlisto(oddx,nn,1);
            deg=ones(length(oddx)/2,1)*2;
            deg=[deg,deg];
        end
        if ~isempty(idxe)
            evenx=x(idxe);% extract the even variables
            evendeg=evenpow(idxe); %extract the even degree;
            Liste=getliste(evenx,evendeg,nn,2);
            len=size(Liste,1);
            edeg=[2*ones(len,1),zeros(len,1)];
            List=[List;Liste];
            deg=[deg;edeg];
        end
    end
end

%%
function Listo=getlisto(oddx,nn,method)
    len=length(oddx);
    id1=1:len;
    id2=randperm(len,len/2);
    id1(id2)=[];
    switch method
        case 1
            SOS1=MPOLY.zeros(nn,len/2,1);
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
            Listo=[SOS1,SOS2];
        case 2
            Listo=[oddx(id1)+oddx(id2);oddx(id1)-oddx(id2)].^2/4;
        otherwise
            error('input error!');
    end
end

%%
function Liste=getliste(evenx,evendeg,nn,method)
    switch method
        case 1
            L=[];
            for i=1:length(evenx)
                L=[L;MPOLY.ones(nn,evendeg(i)/2,1)*evenx(i)^2];
            end
            Liste=[L,MPOLY.zeros(nn,length(L),1)];
        case 2
            k=1;
            Liste=MPOLY.zeros(nn,sum(evendeg)/2,1);
            zerol=Liste;
            for i=1:length(evenx)
                for j=1:evendeg(i)/2
                    Liste(k)=evenx(i);
                    Liste(k).pow=2*Liste(k).pow;
                    k=k+1;
                end
            end
            Liste=[Liste,zerol];
    end
end

%%
function List=ProcedureM(List,deg,method)
    switch method
        case 1
            while size(List,1)>1
                [M,I]=mink(deg(:,1),2);%find first two minimal degree terms
                if deg(I(1),2)==0 && deg(I(2),2)==0
                    temp=[(sum(List(I,1)))^2,sum((List(I,1).^2))]/2;
                elseif deg(I(1),2)~=0 && deg(I(2),2)==0
                    temp=[(sum(List(I,1)))^2+List(I(1),2)^2,...
                        List(I(1),1)^2+(List(I(1),2)+List(I(2),1))^2]/2;
                elseif deg(I(1),2)==0 && deg(I(2),2)~=0
                    temp=[(sum(List(I,1)))^2+List(I(2),2)^2,...
                        List(I(2),1)^2+(List(I(1),1)+List(I(2),2))^2]/2;
                else
                    temp=[(sum(List(I,1)))^2+(sum(List(I,2)))^2,...
                        (List(I(1),1)+List(I(2),2))^2+(List(I(1),2)+List(I(2),1))^2]/2;
                end
                %List(I(1),:)=temp;
                List(I(1),:)=simplify(temp); % this simplification makes polynomials shorter for speedup
                maxdeg=max(M)^2;
                deg(I(1),:)=[maxdeg,maxdeg];
                List(I(2),:)=[];
                deg(I(2),:)=[];
            end
        case 2
            while size(List,1)>1
                [M,I]=mink(deg(:,1),2);%find first two minimal degree terms
                
                temp=[sum(List(I(:),1))^2+sum(List(I(:),2))^2,...
                    (List(I(1),1)+List(I(2),2))^2+(List(I(1),2)+List(I(2),1))^2];
                List(I(1),:)=temp;
                maxdeg=max(M)^2;
                deg(I(1),:)=[maxdeg,maxdeg];
                List(I(2),:)=[];
                deg(I(2),:)=[];
            end
    end
end
