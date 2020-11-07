function [pos,neg] = dbsmatrixform2poly(pos,neg)
    % convert dsos_dbs matrix form output to polylab
    % pos,neg: dbs matrix format polynomial
    pos=simplify(pos.lambda*(pos.b'*pos.v/norm(pos.v))^2);
    neg=simplify(neg.lambda*(neg.b'*neg.v/norm(neg.v))^2);
 end