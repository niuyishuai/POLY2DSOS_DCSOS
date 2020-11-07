function [PSOS,NSOS]=dsos_dbs(p,ismatrixform)
    % Direct basis spectral DSOS decomposition
    % Syntax:
    %   [PSOS,NSOS]=dsos_dbs(p)
    % Inputs:
    % p: MPOLY polynomial
    % ismatrixform: results in matrix form or in polynomial form
    % Outputs:
    % PSOS: a SOS polynomial
    % NSOS: a SOS polynomial
    [c,m]=coefficients(p);
    if sum(m(1).pow)~=0 %~isempty(find(m(1).pow,1))
        c=[0;c];
        m=[MPOLY.ones(m(1).n);m];
    end
    delta=norm(c);
    lambda1=(c(1)+delta)/2;
    lambda2=(c(1)-delta)/2;
    delta1=delta^2-c(1)^2;
    v1=full([delta1;c(2:end)*(delta-c(1))]);
    v2=full([delta1;-c(2:end)*(delta+c(1))]);
    %Q1=tall((lambda1/(v1'*v1))*(v1*v1'));
    %Q2=tall((-lambda2/(v2'*v2))*(v2*v2'));    
    if ismatrixform
        % results in matrix format (fastest)
        %PSOS.Q=tall((lambda1/(v1'*v1))*(v1*v1'));
        %NSOS.Q=tall((-lambda2/(v2'*v2))*(v2*v2'));
        PSOS.v=v1;
        PSOS.lambda=lambda1;
        PSOS.b=m;
        NSOS.v=v2;
        NSOS.lambda=-lambda2;
        NSOS.b=m;
    else
        % quadprod computation (faster)
        PSOS=MPOLY.quadprod((lambda1/(v1'*v1))*(v1*v1'),m);
        NSOS=MPOLY.quadprod((-lambda2/(v2'*v2))*(v2*v2'),m);
        % direct computation (slowest)
        %PSOS=lambda1*(m'*v1/norm(v1))^2;
        %NSOS=-lambda2*(m'*v2/norm(v2))^2;
    end
end
