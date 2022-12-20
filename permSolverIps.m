function [est, perm] = permSolverIps(mix, src)
% permSolverIps solves frequency-wise permutation problem using oracle
% source signals (ideal permutation solver)
%
% [Syntax]
%  [est, perm] = permSolverIps(mix, src)
%
% [Input]
%         mix: complex-valued spectrograms with permutation problem (nFreq x nTime x nSrc)
%         src: oracle source spectrogram (nFreq x nTime x nSrc)
%
% [Output]
%         est: permutation-aligned complex-valued estimated spectrograms (nFreq x nTime x nSrc)
%        perm: estimated permutation (nFreq x nSrc)
%

% Arguments check and set default values
arguments
    mix (:,:,:) double
    src (:,:,:) double
end
[nFreq, nTime, nSrc] = size(mix, [1, 2, 3]);
if isreal(mix); error("'mix' must be complex-valued spectrograms.\n"); end
if isreal(src); error("'src' must be complex-valued spectrograms.\n"); end
if ~isequal(size(mix), size(src)); error("Sizes of 'mix' and 'src' must be equal.\n"); end

% Align estimated spectrogram using oracle source spectrogram
est = zeros(nFreq, nTime, nSrc);
perm = zeros(nFreq, nSrc);
permAll = perms(1:nSrc);
nPerm = size(permAll, 1); % = factorial(nSrc)
for iFreq = 1:nFreq
    for iPerm = 1:nPerm
        err(iPerm) = 0;
        for iSrc = 1:nSrc
            mixInd = permAll(iPerm, iSrc);
            err(iPerm) = err(iPerm) + sum(abs(mix(iFreq, :, mixInd) - src(iFreq, :, iSrc)).^2, "all");
        end
    end
    [~, ind] = min(err);
    perm(iFreq, :) = permAll(ind, :);
    est(iFreq, :, :) = mix(iFreq, :, perm(iFreq,:));
end

fprintf("Permutation solver (IPS) done.\n");
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%