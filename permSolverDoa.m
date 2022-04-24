function [est, perm] = permSolverDoa(demixMat, mix, micPos, sampFreq)
% permSolverDoa solves frequency-wise permutation problem using direction
% of arrivals (DOA)
%
% [Syntax]
%  [est, perm] = permSolverDoa(mix, micPos, sampFreq)
%
% [Input]
%    demixMat: demixing matrix estimated by former FDICA (nSrc x nCh, nFreq)
%         mix: complex-valued spectrograms with permutation problem (nFreq x nTime x nSrc)
%      micPos: position of each microphone [m] (1 x nSrc)
%    sampFreq: sampling frequency of observed signal [Hz] (scalar)
%
% [Output]
%         est: permutation-aligned complex-valued estimated spectrograms (nFreq x nTime x nSrc)
%        perm: estimated permutation (nFreq x nSrc)
%
% [Note]
%    This function requires Statistics and Machine Learning Toolbox for
%    kmeans function
%
% Reference
%   H. Saruwatari, T. Kawamura, T. Nishikawa, A. Lee and K. Shikano, "Blind
%   source separation based on a fastconvergence algorithm combining ICA 
%   and beamforming," IEEE Trans. ASLP, vol. 14, no. 2, pp.666–678, 2006.
%

% Arguments check and set default values
arguments
    demixMat (:,:,:) {mustBeNumeric}
    mix (:,:,:) double
    micPos (1,:) {mustBeNonnegative}
    sampFreq (1,1) {mustBePositive}
end
[nFreq, nTime, nSrc] = size(mix, [1, 2, 3]);
nCh = size(demixMat, 2);
if isreal(mix); error("'mix' must be complex-valued spectrograms.\n"); end
if nCh ~= nSrc; error("nCh must be equal to nSrc.\n"); end
if numel(micPos) ~= nSrc; error("numel(micPos) must be equal to size(mix, 3).\n"); end
if nSrc >= 3; error("permSolverDoa for nSrc >= 3 has not been implemented yet.\n"); end

% Calculate DOA
soundSpd = 340; % sound speed [m/s]
freqAx = linspace(0, sampFreq/2, nFreq).'; % nFreq x 1
sinDoa = zeros(nFreq, nSrc);
doa = zeros(nFreq, nSrc);
mixMat = zeros(nCh, nSrc, nFreq);
for iFreq = 1:nFreq
    mixMat(:, :, iFreq) = inv(demixMat(:,:,iFreq)); % mixMat(:, n, iFreq) is a steering vector of n-th source
end
isValid = zeros(nFreq, nSrc);
for iSrc = 1:nSrc
    zeroCheck(1, 1) = 1;
    zeroCheck(2:nFreq, 1) = (mixMat(iSrc, 2, 2:nFreq)==0).*isinf(mixMat(iSrc, 2, 2:nFreq));
    nonZero = find(~zeroCheck);
    srcAngle = permute(angle(mixMat(1, iSrc, nonZero)./mixMat(2, iSrc, nonZero)), [3, 1, 2]); % nFreq x nCh x nSrc
    sinDoa(nonZero, iSrc)  = srcAngle ./ (2*pi*freqAx(nonZero) * abs(micPos(1)-micPos(2))) * soundSpd;
    doa(nonZero, iSrc)     = real(asin(sinDoa(nonZero, iSrc))) / pi*180;
    isValid(nonZero, iSrc) = abs(sinDoa(nonZero, iSrc)) < 1;
end

% Apply k-means clustering
isValid = (isValid == 1); % convert to logical values
validDOA = doa(isValid);
[~, centroids] = kmeans(validDOA, nSrc);
boundary = mean(centroids); % boundary DOA of clusters

% Align estimated spectrogram based on estimated permutation
est = zeros(nFreq, nTime, nSrc);
perm = zeros(nFreq, nSrc);
for iFreq = 1:nFreq
    if doa(iFreq, 1) <= boundary && doa(iFreq, 2) >= boundary
        perm(iFreq, :) = [1, 2];
    else
        perm(iFreq, :) = [2, 1];
    end
    est(iFreq, :, :) = mix(iFreq, :, perm(iFreq, :));
end

fprintf("Permutation solver (DOA) done.\n");
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%