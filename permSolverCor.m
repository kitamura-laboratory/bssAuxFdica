function [est, perm] = permSolverCor(mix, isPowRatio, type, deltaFreq, ratioFreq)
% permSolverCor solves frequency-wise permutation problem using power
% ratios and clustering algorithm
%
% [Syntax]
%  [est, perm] = permSolverCor(mix)
%  [est, perm] = permSolverCor(mix, isPowRatio)
%  [est, perm] = permSolverCor(mix, isPowRatio, type)
%  [est, perm] = permSolverCor(mix, isPowRatio, type, deltaFreq)
%  [est, perm] = permSolverCor(mix, isPowRatio, type, deltaFreq, ratioFreq)
%
% [Input]
%         mix: complex-valued spectrograms with permutation problem (nFreq x nTime x nSrc)
%  isPowRatio: use power ratio feature for clustering (true/false, false uses amplitude spectrogram for clustering, default: true)
%        type: type of cost function ("Gl", "Lo", or "Gl+Lo", default: "Gl+Lo")
%   deltaFreq: adjacent frequencies for "Lo" cost (default: 3, 0 means adjacent cost is not used)
%   ratioFreq: harmoinc frequencies for "Lo" cost (if ratioFreq=3, round(iFreq/2), round(iFreq/3), 2*iFreq, and 3*iFreq and their adjacent frequencies (e.g., 2*iFreq-1 and 2*iFreq+1) are considered, default: 2, 0 means harmonic cost is not used)
%
% [Output]
%         est: permutation-aligned complex-valued estimated spectrograms (nFreq x nTime x nSrc)
%        perm: estimated permutation (nFreq x nSrc)
%
% Reference
%   H. Sawada, S. Araki, and S. Makino, "Measuring dependence of bin-wise
%   separated signals for permutation alignment in frequency-domain BSS,"
%   Proc. ISCAS, pp. 3247-3250, 2007.
%

% Arguments check and set default values
arguments
    mix (:,:,:) double
    isPowRatio (1,1) logical = true
    type {mustBeMember(type,{'Gl', 'Lo', 'Gl+Lo'})} = "Gl+Lo"
    deltaFreq (1,1) {mustBeNonnegative} = 3
    ratioFreq (1,1) {mustBeNonnegative} = 2
end
[nFreq, nTime, nSrc] = size(mix, [1, 2, 3]);
if isreal(mix); error("'mix' must be complex-valued spectrograms.\n"); end

% Calculate feature vector for clustering
if isPowRatio
    v = abs(mix).^2 ./ sum(abs(mix).^2, 3); % power ratio Eq. (14), nFreq x nTime x nSrc
else
    v = abs(mix); % amplitude spectrogram Eq. (12), nFreq x nTime x nSrc
end

% Iterative k-means clustering
nPerm = factorial(nSrc); % number of patterns in permutation
allPerm = perms((1:nSrc)); % all permutation patterns in nSrc sources case, nPerm x nSrc
perm = repmat(allPerm(end,:), [nFreq, 1]); % initial permutation so that vPerm = v, nFreq x nSrc
vPerm = v; % initial permutation-fixed feature vector
sumRho = zeros(nPerm, 1); % variable for storing cost in Eq. (18)
fprintf("Iteration:    "); iIter = 1;
while(true)
    fprintf("\b\b\b\b%4d", iIter);

    % Store current permutation
    permOld = perm; % permutation of previous iteration

    if type == "Gl"
        % Calculate centroid vector
        c = squeeze(mean(vPerm, 1)); % Eq. (17), nTime x nSrc
        
        for iFreq = 1:nFreq
            % Calculate feature vector
            vf = squeeze(v(iFreq, :, :)); % nTime x nSrc

            % Calculate correlations between iFreq and centroid
            rhoGl = corr(vf, c); % correlation between feature and centroid vectors

            % Calucule cost function value
            for iPerm = 1:nPerm % calc Eq. (18) for all permutation patterns
                sumRho(iPerm, 1) = sum(diag(rhoGl(:, allPerm(iPerm, :)))); % diagonal elements of "rho(:, allPerm(iPerm, :))" are permuted combination
            end

            % Update permutation by maximizing cost
            [perm, vPerm] = local_updatePerm(sumRho, allPerm, perm, v, vPerm, iFreq);
        end

    elseif type == "Lo"
        for iFreq = 1:nFreq
            % Calculate feature vector
            vf = squeeze(v(iFreq, :, :)); % nTime x nSrc

            % Define set of local frequency for vg, i.e., R(f)
            localFreqSet = local_produceLocalFreqSet(iFreq, nFreq, deltaFreq, ratioFreq);
            
            % Calculate correlations between iFreq and local frequency set components
            rhoLoFreqwise = local_calcCorrInLocalFreqSet(localFreqSet, nSrc, vPerm, vf);

            % Calucule cost function value
            for iPerm = 1:nPerm % calc Eq. (19) for all permutation patterns
                rhoLo = squeeze(mean(rhoLoFreqwise(:, :, allPerm(iPerm, :)), 1)); % diagonal elements of "rho(:, :, allPerm(iPerm, :))" are permuted combination
                sumRho(iPerm, 1) = sum(diag(rhoLo));
            end

            % Update permutation by maximizing cost
            [perm, vPerm] = local_updatePerm(sumRho, allPerm, perm, v, vPerm, iFreq);
        end
        
    else % type == "Gl+Lo"
        % Calculate centroid vector
        c = squeeze(mean(vPerm, 1)); % Eq. (17), nTime x nSrc
        
        for iFreq = 1:nFreq
            % Calculate feature vector
            vf = squeeze(v(iFreq, :, :)); % nTime x nSrc

            % Calculate correlations between iFreq and centroid
            rhoGl = corr(vf, c); % correlation between feature and centroid vectors

            % Define set of local frequency for vg, i.e., R(f)
            localFreqSet = local_produceLocalFreqSet(iFreq, nFreq, deltaFreq, ratioFreq);
            
            % Calculate correlations between iFreq and local frequency set components
            rhoLoFreqwise = local_calcCorrInLocalFreqSet(localFreqSet, nSrc, vPerm, vf);

            % Calucule cost function value
            for iPerm = 1:nPerm % calc Eq. (19) for all permutation patterns
                rhoLo = squeeze(mean(rhoLoFreqwise(:, :, allPerm(iPerm, :)), 1)); % diagonal elements of "rho(:, :, allPerm(iPerm, :))" are permuted combination
                sumRho(iPerm, 1) = sum(diag(rhoGl(:, allPerm(iPerm, :)))) + sum(diag(rhoLo)); % Sum of global and local costs
            end

            % Update permutation by maximizing cost
            [perm, vPerm] = local_updatePerm(sumRho, allPerm, perm, v, vPerm, iFreq);
        end
    end

    % Check convergence
    if all(permOld==perm, 'all'); break; end
    iIter = iIter + 1;
end

% Align signal based on estimated permutation
est = zeros(nFreq, nTime, nSrc);
for iFreq = 1:nFreq
    est(iFreq, :, :) = mix(iFreq, :, perm(iFreq,:));
end

fprintf(" Permutation solver (COR) done.\n");
end

%% Local functions
function localSet = local_produceLocalFreqSet(f, F, delta, ratio)
adjSet = [f-delta:f-1, f+1:f+delta]; % set of adjacent local frequency, i.e., A(f)
harSet = [];
for iRatio = 2:ratio % set of harmonic local frequency, i.e., H(f)
    harSet = [harSet, round(f/iRatio)-1:round(f/iRatio)+1, f*iRatio-1:f*iRatio+1];
end
localSet = unique([adjSet, harSet]); % Union of A(f) and H(f) (and sorting)
localSet = localSet(localSet>=1 & localSet<=F); % frequency index must be in the range [1:nFreq]
end

function rhoFreqwise = local_calcCorrInLocalFreqSet(localSet, N, vPerm, vf)
rhoFreqwise = zeros(numel(localSet), N, N);
for f = localSet
    vg = squeeze(vPerm(f, :, :)); % feature vector, nTime x nSrc
    rhoFreqwise(f, :, :) = corr(vf, vg); % correlation between feature and vg vectors
end
end

function [perm, vPerm] = local_updatePerm(cost, allPerm, perm, v, vPerm, f)
[~, idx] = max(cost); % find index of maximum value
perm(f, :) = allPerm(idx, :); % permutation that maximizes Eq. (18)
vPerm(f, :, :) = v(f, :, perm(f, :)); % update permutation-fixed v for calculating Eq. (17)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%