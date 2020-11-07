function [PSOS,NSOS]=dsos_pp(p,isparal)
    % Parameterized Parity DSOS decomposition for polynomial (MPOLY version)
    % Syntax:
    %   [PSOS,NSOS]=DSOS_PP(p)
    %
    % Inputs:
    % p: MPOLY polynomial
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
    % Extract coefs and monomials of p
    [c,m]=coefficients(p);
    if isparal
        parfor i=1:length(c)
            % for each monomial, compute its DSOS components
            [P(i),N(i)]=parity_mono(c(i),m(i));
        end
    else
        for i=1:length(c)
            % for each monomial, compute its DSOS components
            [P(i),N(i)]=parity_mono(c(i),m(i));
        end
    end
    PSOS=sum(P);
    NSOS=sum(N);
end

function [PSOS,NSOS]=parity_mono(coeff,mono)
    % Parity DSOS decomposition for a single monomial
    method = 1;
    if method == 1
        % set parameter s
        s = 1; % s can be any positive number
        % STEP 2: Procedure O to extract o(x) and e(x)^2 for mono
        deg=mono.pow; % degrees of mono
        odeg=mod(deg,2); % degrees of o
        e2deg=deg-odeg; % degrees of e^2
        n=mono.n;
        
        % STEP 4: DSOS decomposition
        if coeff>0
            [PSOS,NSOS] = dsos_mono(n,odeg,e2deg,deg,coeff,s);
        else
            [NSOS,PSOS] = dsos_mono(n,odeg,e2deg,deg,-coeff,s);
        end
    else
        % set parameter s
        s = 1; % s can be any positive number
        % STEP 2: Procedure O to extract o(x) and e(x)^2 for mono
        deg=mono.pow; % degrees of mono
        odeg=mod(deg,2); % degrees of o
        edeg=(deg-odeg)/2; % degrees of e
        e=mono;
        e.pow=edeg;
        o=mono;
        o.pow=odeg;
        
        % STEP 4: DSOS decomposition
        if coeff>0
            PSOS = (e*o+s*e)^2*coeff/(2*s);
            NSOS = (e^2*o^2+s^2*e^2)*coeff/(2*s);
        else
            PSOS = (e^2*o^2+s^2*e^2)*(-coeff)/(2*s);
            NSOS = (e*o+s*e)^2*(-coeff)/(2*s);
        end
    end
end

function [L,R]=dsos_mono(n,odeg,e2deg,mdeg,c,s)
    L=MPOLY.zeros(n);
    R=L;
    L.k=3;
    L.coef=[c/2*s;c*s/2;c];
    L.pow=[mdeg+odeg;e2deg;mdeg];
    R.k=2;
    R.coef=L.coef(1:2);
    R.pow=L.pow(1:2,:);
end
