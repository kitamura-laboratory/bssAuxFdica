function [estSig, cost] = bssAuxFdica(obsSig, nSrc, args)
% bssAuxFdica: blind source separation based on frequency-domain 
%              independent component analysis with auxiliary-function
%              technique
% [Syntax]
%   Using traditional name-value pair expression:
%   [estSig, cost] = bssAuxFdica(obsSig, nSrc, "fftSize", 1024, 
%                    "shiftSize", 512, "nIter", 50, "isWhiten", true, 
%                    "srcModel", "LAP", "refMic", 1 "permSolver", "COR", 
%                    "isDraw", false, "sampFreq", 16000)
%   Using pythonic expression (later R2021a):
%   [estSig, cost] = bssAuxFdica(obsSig, nSrc, fftSize=1024, 
%                    shiftSize=512, nIter=50, isWhiten=true, 
%                    srcModel="LAP", refMic=1, permSolver="COR", 
%                    isDraw=false, sampFreq=16000)
%
% [Input]
%      obsSig: time-domain observed signal (sample x channel)
%        nSrc: number of sources (scalar)
%     fftSize: window length in STFT (scalar, default: 1024)
%   shiftSize: window shift length in STFT (scalar, default: fftSize/2)
%       nIter: number of iterations for FDICA optimization 
%              (scalar, default: 50)
%    isWhiten: apply whitening before BSS (true/false, default: true)
%    srcModel: generative model of each source 
%              ("LAP" or "TVG", default: "LAP")
%              "LAP": isotropic complex Laplace distribution
%              "TVG": isotropic time-varying complex Gaussian distribution
%      refMic: reference microphone onto which estimated spectrogram is
%              projected by projection back technique 
%              (scalar or row vector, default: 1) 
%  permSolver: type of permutation solver 
%              ("none", "COR", "DOA", or "IPS", default: "COR")
%              "none": do not apply permutation solver after FDICA
%              "COR": correlation-based permutation solver
%              "DOA": direction-of-arrivals-based permutation solver
%              "IPS": ideal permutation solver using oracle source signals 
%                     (for checking upper bound performance of FDICA)
%      isDraw: draw spectrograms and cost function behavior 
%              (true/false, default: false)
%   sampFreq*: sampling frequency of observed signal [Hz]
%              (scalar, default: 16000, used for plotting spectrogram and 
%              DOA-based permutation solver)
% isPowRatio*: use power ratio feature for clustering (true/false, 
%              false uses raw amplitude spectrogram, default: true, used 
%              for COR-based permutation solver)
%    typeCor*: type of cost function 
%              ("Gl", "Lo", or "Gl+Lo", default: "Gl+Lo", used for
%              COR-based permutation solver)
%              "Gl": use global correlation
%              "Lo": use local correlation
%              "Gl+Lo": use both global and local correlations
%  deltaFreq*: adjacent frequencies for "Lo" correlation 
%              (scalar, default: 3, 0 means adjacent cost is not used, used
%              for COR-based permutation solver)
%  ratioFreq*: harmoinc frequencies for "Lo" correlation (scalar, 
%              if ratioFreq=3, round(iFreq/2), round(iFreq/3), 2*iFreq, 
%              and 3*iFreq and their adjacent frequencies 
%              (e.g., 2*iFreq-1 and 2*iFreq+1) are considered, default: 2,
%              0 means harmonic cost is not used, used for COR-based 
%              permutation solver)
%     micPos*: position of each microphone [m] (row vector, used for 
%              DOR-based permutation solver)
%     srcSig*: oracle source image signals used for "IPS" permutation
%              solver (sample x channel x source)
%     isFilt*: apply time-domain demixing filter to avoid circular
%              convoluation
% Arguments with * are not necessary
%
% [Output]
%      estSig: estimated signal obtained by FDICA (sample x source)
%        cost: cost function values of FDICA in each iteration (nIter x 1)
%    demixMat: estimated demixing matrix (nSrc x channel x fftSize/2+1)
%
% [Note]
%    This function requires the following functions:
%        DGTtool.m (https://github.com/KoheiYatabe/DGTtool)
%        permSolverCor.m
%        permSolverDoa.m
%        permSolverIps.m
%

% Check arguments and set default values
arguments
    obsSig (:,:) double
    nSrc (1,1) double {mustBeInteger, mustBePositive}
    args.fftSize (1,1) double {mustBeInteger, mustBePositive} = 1024
    args.shiftSize (1,1) double {mustBeInteger, mustBePositive} = 512
    args.nIter (1,1) double {mustBeInteger, mustBePositive} = 50
    args.isWhiten (1,1) logical = true
    args.srcModel (1,1) string {mustBeMember(args.srcModel, ["LAP", "TVG"])} = "LAP"
    args.refMic (1,:) double {mustBeInteger, mustBePositive} = 1
    args.permSolver (1,1) string {mustBeMember(args.permSolver, ["none", "COR", "DOA", "IPS"])} = "COR"
    args.isDraw (1,1) logical = false
    args.sampFreq (1,1) double {mustBePositive} = 16000
    args.isPowRatio (1,1) logical = true
    args.typeCor (1,1) string {mustBeMember(args.typeCor, ["Gl", "Lo", "Gl+Lo"])} = "Gl+Lo"
    args.deltaFreq (1,1) double {mustBeInteger, mustBeNonnegative} = 3
    args.ratioFreq (1,1) double {mustBeInteger, mustBeNonnegative} = 2
    args.micPos (1,:) double {mustBeNonnegative}
    args.srcSig (:,:,:) double {mustBeNumeric}
    args.isFilt (1,1) logical = false
end
fftSize = args.fftSize;
shiftSize = args.shiftSize;
nIter = args.nIter;
isWhiten = args.isWhiten;
srcModel = args.srcModel;
refMic = args.refMic;
permSolver = args.permSolver;
isDraw = args.isDraw;
isFilt = args.isFilt;

% Check argument errors
[sigLen, nCh] = size(obsSig, [1, 2]);
if nSrc > nCh; error("'nSrc' must be equal or less than size(obsSig, 2).\n"); end
if fftSize < shiftSize; error("'shiftSize' must be equal or less than fftSize.\n"); end
if numel(refMic) > nCh; error("numel(refMic) must be equal or less than size(obsSig, 2).\n"); end

% Caluculate STFT
F = DGTtool("windowName", "b", "windowLength", fftSize, "windowShift", shiftSize); % create DGTtool instance
obsSpec = F.DGT(obsSig); % STFT
nFreq = size(obsSpec, 1);

% Apply whitening (decorrelation and normalization of observed signals)
if isWhiten
    obsSpecInput = local_whitening(obsSpec, nSrc);
else
    obsSpecInput = obsSpec(:, :, 1:nSrc); % discard unnecessary channels
end

% Apply FDICA
[estSpecFdica, demixMat, cost] = local_auxFdica(obsSpecInput, nIter, srcModel, isDraw);

% Apply projection back technique
[estSpecFdicaFix, demixMatFix] = local_projectionBack(estSpecFdica, obsSpec(:,:,refMic), demixMat);

% Apply permutation solver
if permSolver == "none"
    estSpec = estSpecFdicaFix;
    estPerm = repmat(1:nSrc, [nFreq, 1]);
elseif permSolver == "COR"
    [estSpec, estPerm] = permSolverCor(estSpecFdicaFix, args.isPowRatio, args.typeCor, args.deltaFreq, args.ratioFreq);
elseif permSolver == "DOA"
    [estSpec, estPerm] = permSolverDoa(demixMatFix, estSpecFdicaFix, args.micPos, args.sampFreq);
else % IPS
    srcSpect = F.DGT(squeeze(args.srcSig(:, args.refMic, :)));
    [estSpec, estPerm] = permSolverIps(estSpecFdicaFix, srcSpect);
end
for iFreq = 1:nFreq
    demixMatFix(:, :, iFreq) = demixMatFix(estPerm(iFreq, :), :, iFreq);
end

% Calculate estimated time-domain signal
if isFilt
    % Apply demixing filter in time domain to avoid circular convolution
    obsSigInput = F.pinv(obsSpecInput); % observed signal
    W = cat(3, demixMatFix, flip(conj(demixMatFix(:, :, 2:end-1)), 3)); % produce beyond Nyquist components
    demixFilt = real(ifft(W, fftSize, 3)); % fftSize x nSrc x nMic
    demixFilt = circshift(demixFilt, fftSize/2+1, 3); % move peak to center by circular shifting
    for iSrc = 1:nSrc
        for iCh = 1:nCh
            f = squeeze(demixFilt(iSrc, iCh, :));
            tmp(:, iCh) = conv(obsSigInput(:, iCh), f); % linear convolution
        end
        estSig(:, iSrc) = sum(tmp, 2);
    end
    estSig(1:fftSize/2+1,:) = []; % cut initial components caused by group delay (circular shifting)
else
    % Calculate inverse STFT
    estSig = F.pinv(estSpec);
end
estSig = estSig(1:sigLen, :);

% Plot spectrograms and cost function behavior
if isDraw
    F.plot(obsSig, args.sampFreq); % observed signal
    F.plot(F.pinv(obsSpecInput), args.sampFreq); % observed signal input to FDICA
    F.plot(F.pinv(estSpecFdica), args.sampFreq); % estimated spectrogram before projection-back technique
    if permSolver ~= "none"
        F.plot(F.pinv(estSpecFdicaFix), args.sampFreq); % estimated spectrogram before permutation solver
    end
    F.plot(estSig, args.sampFreq); % estimated signal
    local_plotCost(cost, nIter); % cost function behavior
end
end

%% Local functions
%--------------------------------------------------------------------------
function [Y, dP] = local_whitening(X, N)
% Whitening based on frequency-wise principal component analysis
%
% [inputs]
%    X: input spectrogram (I x J x M, nFreq x nTime x nCh)
%    N: number of sources (dimensions to which X(i,:,:) is projected)
%
% [outputs]
%    Y: output matrix (I x J x N)
%

% Initialize
[I, J, ~] = size(X, [1, 2, 3]); % nFreq x nTime x nCh
Y = zeros(I, J, N);

% Apply frequency-wise whitening
Xp = permute(X, [3, 2, 1]); % I x J x M -> M x J x I
for i = 1:I
    Xi = Xp(:, :, i); % M x J
    V = Xi*(Xi')/J; % covariance matrix of data matrix X (K x K)
    [P, D] = eig(V); % eigenvalue decomposition (V = P*D*inv(P), P includes eigenvectors and D is a diagonal matrix with eigenvalues)
    [~, idx] = sort(diag(D), "descend"); % sort eigenvalues in descending order
    D = D(idx, idx); % sorted D
    P = P(:, idx); % sorted P
    dP = P(:, 1:N); % top-d eigenvectors
    Yi = sqrt(D)\(dP')*Xi; % whitened vector (N x J)
    Y(i, :, :) = Yi.'; % J x N
end
end

%--------------------------------------------------------------------------
function [Y, W, cost] = local_auxFdica(X, nIter, srcModel, isDraw)
% BSS using FDICA
%
% [inputs]
%        X: observed spectrogram (I x J x M, nFreq x nTime x nCh, nCh=nSrc)
%    nIter: number of iterations
% srcModel: generative model of each source ("LAP" or "TVG")
%   isDraw: draw cost function behavior or not
%
% [outputs]
%        Y: estimated spectrogram (I x J x N, nFreq x nTime x nSrc)
%     cost: cost function values of FDICA in each iteration (nIter x 1)
%

% Initialize
[I, J, M] = size(X, [1,2,3]); % nFreq x nTime x nCh
N = M; % number of sources
E = repmat(eye(M), [1, 1, I]);
W = E; % initial demixing matrix (N x M x I)
Y = X; % initial estimated spectrogram
Xp = permute(X, [3, 2, 1]); % M x J x I
Xph = pagectranspose(Xp); % J x M x I, pagewise Hermitian transpose (Xp')
Yp = permute(Y, [3, 2, 1]); % N x J x I
cost = zeros(nIter, 1);
if isDraw
    cost(1,1) = local_calcFdicaCost(Yp, W, srcModel, I, J);
end

% Optimize
fprintf("Iteration:    ");
for iIter = 1:nIter
    fprintf("\b\b\b\b%4d", iIter);
    if srcModel == "LAP"
        Rp = max(abs(Yp), 10000*eps);
    elseif srcModel == "TVG"
        Rp = max(abs(Yp).^2, 10000*eps);
    end
    
    invRp = 1./Rp; % N x J x I
    for n = 1:N
        D = repmat(invRp(n, :, :), [M, 1, 1]); % M x J x I
        Vk = pagemtimes(D.*Xp, Xph)/J; % M x M x I, pagewise matrix multiplication ((D(:,:,i).*Xp(:,:,i))*Xp(:,:,i)'/J)
        wn = pagemldivide(pagemtimes(W, Vk), E(:, n, :)); % M x 1 x I, pagewise operation ((W(:,:,i)*Vk(:,:,i)) \ E(:, n, :))
        wn = wn ./ sqrt( pagemtimes(pagemtimes(wn, "ctranspose", Vk, "none"), wn) ); % M x 1 x I, pagewise operation (wn(:,:,i)/sqrt(wn(:,:,i)'*Vk(:,:,i)*wn(:,:,i)))
        wnh = pagectranspose(wn); % 1 x M x I, pagewise Hermitian transpose (wn(:,:,i)')
        Yp(n, :, :) = pagemtimes(wnh, Xp); % 1 x J x I, pagewise matrix multiplication (wnh(:,:,i)*Xp(:,:,i))
        W(n, :, :) = wnh;
    end
% Readable implimentation
%     for i = 1:I
%         for n = 1:N
%             rn = Rp(n, :, i); % 1 x J
%             dg = ones(M, 1)*(1./rn); % M x J
%             Vk = (dg.*Xp(:, :, i))*Xp(:, :, i)'/J; % M x M
%             wn = (W(:, :, i)*Vk) \ E(:, n, i);
%             wn = wn/sqrt((wn')*Vk*wn);
%             Yp(n, :, i) = (wn')*Xp(:, :, i);
%             W(n, :, i) = wn';
%         end
%     end
    if isDraw
        cost(iIter+1, 1) = local_calcFdicaCost(Yp, W, srcModel, I, J);
    end
end
Y = permute(Yp, [3, 2, 1]);
fprintf(" FDICA done.\n");
end

%--------------------------------------------------------------------------
function costVal = local_calcFdicaCost(Yp, W, srcModel, I, J)
detW = zeros(I, 1);
for i = 1:I
    detW(i,1) = det(W(:,:,i));
end
if srcModel == "LAP"
    costVal = sum(abs(Yp), "all") - 2*J*sum(log(abs(detW)));
elseif srcModel == "TVG"
    costVal = sum(log(max(abs(Yp).^2, eps)), "all") - 2*J*sum(log(abs(detW)));
end
end

%--------------------------------------------------------------------------
function [fixY, fixW] = local_projectionBack(Y, S, W)
% Projection back technique to fix frequency-wise scales of estimated
% spectrogram obtained by FDICA
%
% [inputs]
%      Y: estimated spectrograms (I x J x N, nFreq x nTime x nSrc)
%      S: reference channel of observed spectrogram (I x J x 1)
%         or observed multichannel spectrogram (I x J x M, nFreq x nTime x nMic)
%      W: estimated emixing matrix (N x N x I, nSrc x nCh x nFreq)
%
% [outputs]
%   fixY: scale-fixed estimated spectrograms (I x J x N)
%         or scale-fitted estimated source images (I x J x N x M)
%   fixW: scale-fixed demixing matrix (N x N x I)
%         or scale-fitted demixing matrix for source images (N x N x I x M)
%

% Projection back
Yp = permute(Y, [3, 2, 1]); % N x J x I
Sp = permute(S, [3, 2, 1]); % 1 x J x 1 or M x J x I
Wp = permute(W, [4, 1, 2, 3]); % 1 x N x N x I
Yph = pagectranspose(Yp); % J x N x I, pagewise Hermitian transpose (Yp')
YpYph = pagemtimes(Yp, Yph); % N x N x I, pagewise matrix multiplication (Yp*Yp')
YphOnYpYph = pagemrdivide(Yph, YpYph); % J x N x I, pagewise matrix right-division (Yp'/(Yp*Yp'))
A = pagemtimes(Sp, YphOnYpYph); % 1 x N x I or M x N x I, pagewise matrix multiplication (Sp * Yp'/(Yp*Yp'))
Ap = permute(A, [1, 2, 4, 3]); % M x N x 1 x I
Ypp = permute(Yp, [4, 1, 2, 3]); % 1 x N x J x I
fixY = Ap .* Ypp; % M x N x J x I, using implicit expansion
fixY = permute(fixY, [4, 3, 2, 1]); % I x J x N x M
fixW = Ap .* Wp; % M x N x N x I, using implicit expansion
fixW = permute(fixW, [2, 3, 4, 1]); % N x N x I x M

% Readable implementation
% [I, J, N] = size(Y, [1, 2, 3]); % nFreq x nTime x nSrc
% M = size(S, 3); % nCh
% A = zeros(M, N, I); % frequency-wise projection matrix
% for i = 1:I
%     for m = 1:M
%         Yi = permute(Y(i, :, :), [3, 2, 1]); % I x J x N -> N x J x 1
%         A(m, :, i) = S(i, :, m)*Yi'/(Yi*Yi');
%     end
% end
% fixY = zeros(I, J, N, M); % scale-fixed estimated spectrograms
% fixW = zeros(N, N, I, M); % scale-fixed demixing matrix
% for n = 1:N
%     for m = 1:M
%         for i = 1:I
%             fixY(i, :, n, m) = A(m, n, i)*Y(i, :, n);
%             fixW(n, :, i, m) = A(m, n, i)*W(n, :, i);
%         end
%     end
% end
end

%--------------------------------------------------------------------------
function local_plotCost(cost, nIter)
figure; plot(0:nIter, cost);
set(gca, "FontSize", 12);
xlabel("Number of iterations"); ylabel("Value of cost function");
grid on;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%