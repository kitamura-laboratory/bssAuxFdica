% Main script for frequency-domain independent component analysis (FDICA)
% Corded by D. Kitamura (d-kitamura@ieee.org) on April 23rd, 2022

clear; close all; clc;
addpath("./bss_eval/");

% Set parameters
seed = 1; % pseudorandom seed
fftSize = 4096; % window length in STFT [points]
shiftSize = fftSize/2; % window shift length in STFT [points]
nSrc = 2; % number of sources in observed signal
nIter = 50; % number of iterations of FDICA
isWhiten = true; % apply whitening before FDICA or not (true/false)
srcModel = "LAP"; % generative model of each source ("LAP" or "TGV")
refMic = 1; % index of reference microphone for projection back technique
permSolver = "COR"; % type of permutation solver ("none", "COR", "DOA", or "PPS")
isDraw = true; % plot spectrograms and cost function behavior for debug (true/false)
micPos(1) = 0; % position of the first microphone [m]
micPos(2) = 0.0566; % position of the second microphone [m]
dataNo = 1; % file number of input data (see getInputFileNames) (1-8)

%% Preprocessing
% Set pseudorandom seed
rng(seed);

% Get input file names
[dirPath, fileName] = getInputFileNames(dataNo);

% Read input source image files
for iSrc = 1:nSrc
    filePath = dirPath + fileName(iSrc);
    [srcSig(:,:,iSrc), fs] = audioread(filePath); % srcSig: sample x mic x source
end

% Mix source images
obsSig = sum(srcSig, 3); % obsSig: sample x mic

% Check wave clipping
peakVal = max(abs(obsSig), [], "all");
if  peakVal > 1 % clipped
    obsSig = 0.99 * obsSig / peakValue; % maximum value is set to 0.99
    refSig = 0.99 * squeeze(srcSig(:, refMic, :)) / peakValue; % refSig: sample x source
    fprintf('Observed signal is normalized during mixture.\n');
else
    refSig = squeeze(srcSig(:, refMic, :)); % refSig: sample x source
end

%% BSS based on FDICA and permutation solver
% Sample for permSolver="COR"
estSig = bssAuxFdica(obsSig, nSrc, ...
    "fftSize", fftSize, "shiftSize", shiftSize, "nIter", nIter, ...
    "isWhiten", isWhiten, "srcModel", srcModel, "refMic", refMic, ...
    "permSolver", permSolver, "isDraw", isDraw, "sampFreq", fs);

% Sample for permSolver="DOA"
% estSig = bssAuxFdica(obsSig, nSrc, ...
%     "fftSize", fftSize, "shiftSize", shiftSize, "nIter", nIter, ...
%     "isWhiten", isWhiten, "srcModel", srcModel, "refMic", refMic, ...
%     "permSolver", permSolver, "isDraw", isDraw, "sampFreq", fs, "micPos", micPos);

% Sample for permSolver="IPS"
% estSig = bssAuxFdica(obsSig, nSrc, ...
%     "fftSize", fftSize, "shiftSize", shiftSize, "nIter", nIter, ...
%     "isWhiten", isWhiten, "srcModel", srcModel, "refMic", refMic, ...
%     "permSolver", permSolver, "isDraw", isDraw, "sampFreq", fs, "srcSig", srcSig);

%% Evaluation of BSS performance

% Calculate input SDR and SIR
[inSdr, inSir, inSar] = bss_eval_sources(repmat(obsSig(:, refMic), [1, nSrc]).', refSig.');

% Calculate output SDR, SIR, and SAR
[outSdr, outSir, outSar] = bss_eval_sources(estSig.', refSig.');

% Display improvements of SDR and SIR and raw SAR
for iSrc = 1:nSrc
    impSdr(iSrc, 1) = outSdr(iSrc, 1) - inSdr(iSrc, 1);
    impSir(iSrc, 1) = outSir(iSrc, 1) - inSir(iSrc, 1);
    rawSar(iSrc, 1) = outSar(iSrc, 1);
    fprintf('Source %d\n  SDRi: %.2f[dB], SIRi: %.2f[dB], SAR: %.2f[dB]\n', iSrc, impSdr(iSrc, 1), impSir(iSrc, 1), rawSar(iSrc, 1));
end

%% Output estimated wave files
outDir = "./output/";
if ~isfolder(outDir); mkdir(outDir); end
audiowrite(outDir+sprintf("data%d", dataNo)+"_obs.wav", obsSig, fs); % observed signal
audiowrite(outDir+sprintf("data%d", dataNo)+"_src1.wav", refSig(:, 1), fs); % source signal
audiowrite(outDir+sprintf("data%d", dataNo)+"_src2.wav", refSig(:, 2), fs); % source signal
audiowrite(outDir+sprintf("data%d", dataNo)+"_est1.wav", estSig(:, 1), fs); % estimated signal
audiowrite(outDir+sprintf("data%d", dataNo)+"_est2.wav", estSig(:, 2), fs); % estimated signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%