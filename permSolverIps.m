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
if nSrc >= 4; error("permSolverIps for nSrc >= 4 has not been implemented yet.\n"); end

% Align estimated spectrogram using oracle source spectrogram
est = zeros(nFreq, nTime, nSrc);
perm = zeros(nFreq, nSrc);
if nSrc == 2
    for iFreq = 1:nFreq
        err(1) = sum(abs(mix(iFreq, :, 1) - src(iFreq, :, 1)).^2 + abs(mix(iFreq, :, 2) - src(iFreq, :, 2)).^2, "all");
        err(2) = sum(abs(mix(iFreq, :, 2) - src(iFreq, :, 1)).^2 + abs(mix(iFreq, :, 1) - src(iFreq, :, 2)).^2, "all");
        if err(1) <= err(2)
            perm(iFreq, :) = [1, 2];
        else
            perm(iFreq, :) = [2, 1];
        end
        est(iFreq, :, :) = mix(iFreq, :, perm(iFreq,:));
    end
else
    for iFreq = 1:nFreq
        err(1) = sum(abs(mix(iFreq, :, 1) - src(iFreq, :, 1)).^2 + abs(mix(iFreq, :, 2) - src(iFreq, :, 2)).^2 + abs(mix(iFreq, :, 3) - src(iFreq, :, 3)).^2, "all");
        err(2) = sum(abs(mix(iFreq, :, 1) - src(iFreq, :, 1)).^2 + abs(mix(iFreq, :, 3) - src(iFreq, :, 2)).^2 + abs(mix(iFreq, :, 2) - src(iFreq, :, 3)).^2, "all");
        err(3) = sum(abs(mix(iFreq, :, 2) - src(iFreq, :, 1)).^2 + abs(mix(iFreq, :, 1) - src(iFreq, :, 2)).^2 + abs(mix(iFreq, :, 3) - src(iFreq, :, 3)).^2, "all");
        err(4) = sum(abs(mix(iFreq, :, 2) - src(iFreq, :, 1)).^2 + abs(mix(iFreq, :, 3) - src(iFreq, :, 2)).^2 + abs(mix(iFreq, :, 1) - src(iFreq, :, 3)).^2, "all");
        err(5) = sum(abs(mix(iFreq, :, 3) - src(iFreq, :, 1)).^2 + abs(mix(iFreq, :, 1) - src(iFreq, :, 2)).^2 + abs(mix(iFreq, :, 2) - src(iFreq, :, 3)).^2, "all");
        err(6) = sum(abs(mix(iFreq, :, 3) - src(iFreq, :, 1)).^2 + abs(mix(iFreq, :, 2) - src(iFreq, :, 2)).^2 + abs(mix(iFreq, :, 1) - src(iFreq, :, 3)).^2, "all");
        [~, ind] = min(err);
        switch(ind)
            case 1
                perm(iFreq, :) = [1, 2, 3];
            case 2
                perm(iFreq, :) = [1, 3, 2];
            case 3
                perm(iFreq, :) = [2, 1, 3];
            case 4
                perm(iFreq, :) = [2, 1, 3];
            case 5
                perm(iFreq, :) = [3, 1, 2];
            case 6
                perm(iFreq, :) = [3, 2, 1];
        end
        est(iFreq, :, :) = mix(iFreq, :, perm(iFreq,:));
    end
end

fprintf("Permutation solver (IPS) done.\n");
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%