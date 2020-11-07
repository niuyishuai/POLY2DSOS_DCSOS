function [pos,neg,dsos]=poly2dsos(p,method,ismatrixform,isparal)
%POLY2DSOS 
%% DESCRIPTION:
% DSOS decomposition of multivariate polynomials (POLYLAB version)
%% SYNTAX:
%   [pos,neg] = poly2dsos(p,method,ismatrixform,isparal)
%   return two polynomials such that p=pos-neg
%
%% INPUTS:
%   p: polynomial (MPOLY object)
%   method: DSOS decomposition methods 'PP'|'IP'|'DBS'|'MBS'
%   ismatrixform: output in matrix form (only for spectral decompositions)
%   isparal: using parallel computing
%   0: spectral decomposition
%   1: solving sdp
%   2: parity decomposition
%% OUTPUTS:
%   pos: SOS polynomial
%   neg: SOS polynomial
%   two output formats:
%   i)  polynomial format if pos and neg are POLYLAB polynomials with p=pos-neg
%   ii) matrix format (only for spectral decomposition) if pos and neg are structure 
%       with entries Q (double matrix) and b (monomial vector) such that 
%       the polynomial (pos or neg) equals to b'*Q*b.
%
%% COPYRIGHT:
% Copyright since 2017, Yi-Shuai NIU. All Rights Reserved.
% 2017/12/1     Initial code 
% 2020/10/15    Improved using Polylab 
% 2020/10/25    Introduce parallel mode

if nargin<4
    isparal=false;
end
%% Construct DSOS decompositions
dsos=[];
switch method
    case 'PP' % parametrized parity DSOS decomposition
        [pos,neg]=dsos_pp(p,isparal);
    case 'IP' % improved parity DSOS decomposition
        [pos,neg]=dsos_ip(p,isparal);
    case 'DBS'% direct basis spectral DSOS decomposition
        [pos,neg]=dsos_dbs(p,ismatrixform);
    case 'MBS'% minimal basis spectral DSOS decomposition
        [pos,neg,dsos]=dsos_mbs(p,ismatrixform);
end
end
