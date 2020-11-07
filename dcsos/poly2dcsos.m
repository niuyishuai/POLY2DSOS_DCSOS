function [pos,neg]=poly2dcsos(p,method,isparal)
%POLY2DCSOS 
%% DESCRIPTION:
% DCSOS decomposition of multivariate polynomials (POLYLAB version)
%% SYNTAX:
%   [pos,neg] = poly2dcsos(p,method,isparal)
%   return two polynomials such that p=pos-neg
%
%% INPUTS:
%   p: polynomial (MPOLY object)
%   method: DCSOS decomposition methods 'IP'|'MD'|'FMD'
%   isparal: using parallel computing
%% OUTPUTS:
%   pos: CSOS polynomial
%   neg: CSOS polynomial
%
%% COPYRIGHT:
% Copyright since 2017, Yi-Shuai NIU. All Rights Reserved.
% 2018/10/1     Initial code 
% 2020/10/20    Improved using Polylab 
% 2020/10/25    Introduce parallel mode

if nargin<3
    isparal=false;
end
%% Construct DSOS decompositions
switch method
    case 'IP' % improved parity DCSOS decomposition
        [pos,neg]=dcsos_ip(p,isparal);
    case 'MD'% minimal degree DCSOS decomposition
        [pos,neg]=dcsos_md(p,isparal);
    case 'FMD'% formulation based minimal degree DCSOS decomposition
        [pos,neg]=dcsos_fmd(p,isparal);
end
end
