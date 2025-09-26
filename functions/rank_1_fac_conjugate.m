%% extracting conjugate factor from a rank-1 matrix
function [fac,res] = rank_1_fac_conjugate(X)
% clc; clear; close all;
% fac_gt = randn(15,1) + 1i * randn(15,1);
% X = fac_gt * fac_gt';

    [eigvec,eigval] = eig(X);
    [max_val, max_idx] = max(diag(abs(eigval)));
    

    fac = sqrt(max_val) * eigvec(:,max_idx);
    
    Xest = fac * fac';
    
    res = frob(X-Xest)^2;

end