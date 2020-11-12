function [M, H, snr] = hyperMnf(M, h, w)
% HYPERMNF Performs the minimum noise fraction (MNF) transform
%  hyperMnf performs the minimum noise fraction (MNF) transform on the data
% and uses spatial (row) offsets of the data to estimate the covariance
% matrix of the data.
%
% Usage
%   M = hyperMnf(M, h, w)
% Inputs
%   M - 2D matrix (p x N)
%   h - height of image in pixels
%   w - width of image in pixels
% Outputs
%   M - 2D transformed data
%   H - 2D transformation matrix
%
% References
%   C-I Change and Q Du, "Interference and Noise-Adjusted Pricipal 
% Components Analysis," IEEE TGRS, Vol 36, No 5, September 1999.

[p, N] = size(M);

sigmaZ = hyperCov(M);
M = hyperConvert3d(M, h, w, p);

% 1. Estimate the covariance of the noise.
dX = zeros(h-1, w, p);
for i=1:(h-1)
    dX(i, :, :) = M(i, :, :) - M(i+1, :, :);
end
dX = hyperConvert2d(dX);
sigmaN = hyperCov(dX);
[U,S,E] = svd(sigmaN);
F = E*inv(sqrt(S));

sigmaAdj = F'*sigmaZ*F;
[U,S,G] = svd(sigmaAdj);
H = G*F;

% Compute SNR
% At this point, all noise powers are one.  So, the variance in each band is 
% from signal + noise where noise is at unity power.
snr = diag(S)-1;

% Perform transform
M = H*hyperConvert2d(M);