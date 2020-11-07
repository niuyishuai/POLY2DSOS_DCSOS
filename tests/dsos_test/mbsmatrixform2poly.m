function [pos,neg] = mbsmatrixform2poly(dsos)
    % convert dsos_mbs matrix form output to polylab
    % dsos: mbs matrix format polynomial
    lambda1=dsos.L;
    lambda2=dsos.L;
    lambda1(lambda1<0)=0;
    lambda2(lambda2>0)=0;
    pos=MPOLY.quadprod(dsos.P*lambda1*dsos.P',dsos.b);
    neg=MPOLY.quadprod(dsos.P*(-lambda2)*dsos.P',dsos.b);
end