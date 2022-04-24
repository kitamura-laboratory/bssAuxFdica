classdef DGTtool < handle
    %DGTtool: A simple and user-friendly tool for computing STFT/DGT.
    %   (MATLAB 2020b or later and Signal Processing Toolbox are required.)
    %
    %
    %   --- Quick start ---
    %
    %   DGTtool provides an easy way to compute STFT/DGT. 
    %    - 1st step: Construct a DGTtool object.
    %         F = DGTtool
    %    - 2nd step: Compute spectrogram X from a time-domain signal x.
    %         X = F(x);
    %   That's all for computing STFT/DGT.
    %
    %   The following command reconstructs the time-domain signal.
    %         x = F.pinv(X);
    %   Note: The reconstructed signal may be longer than the original one.
    %
    %
    %   --- Parameter settings ---
    %
    %   Parameters of STFT/DGT can be set as follows:
    %      F = DGTtool('windowShift',100,'windowLength',1000,'FFTnum',500)
    %
    %   The acceptable parameters are
    %      'windowShift'  (positive integer)
    %      'windowLength' (positive integer)
    %      'FFTnum'       (positive integer)
    %      'windowName'   (char or string)
    %      'windowVector' (column vector)
    %
    %   List of window names is given by typing DGTtool.windowList at the command line.
    %   Partial matching of leading characters is supported (case-insensitive).
    %
    %   The order of parameters can be altered as follows:
    %      F = DGTtool('windowName','b','windowLength',1000,'windowShift',100)
    %
    %   Note: MATLAB 2021a or later allows the Name=Value syntax:
    %      F = DGTtool(windowName='b', windowLength=1000, windowShift=100)
    %
    %
    %   --- Methods ---
    %
    %   The following methods can be called as F.methodName(inputVars).
    %
    %   DGTtool Methods:
    %   Forward Transforms
    %      subsref  - Implemented for shortcut notation of DGT (F(x))
    %      DGT      - Compute STFT/DGT (F(x) is shortcut of F.DGT(x))
    %      reassign - Compute sparse (reassigned) spectrogram
    %
    %   Inverse transforms
    %      H    - Inverse STFT/DGT (complex conjugate transpose of F)
    %      pinv - Inverse STFT/DGT with perfect reconstruction (pseudo-inverse of F)
    %
    %   Plot functions
    %      plot         - Draw spectrogram
    %      plotPhase    - Visualize phase spectrogram
    %      plotReassign - Draw sparse (reassigned) spectrogram
    %
    %   Window utilities
    %      setWindow       - Change window
    %      makeWindowTight - Compute canonical tight window
    %      plotWin         - Draw windows
    %
    %   Phase manipulation
    %      makeZeroPhase     - Remove linear phase component of window
    %      undoMakeZeroPhase - Cancel the effect of makeZeroPhase
    %      changeDGTdef      - Convert DGT definition
    %      undoChangeDGTdef  - Cancel the effect of changeDGTdef
    %
    %
    %   --- Static methods ---
    %
    %   Static methods can be used without constructing a DGTtool object.
    %   To use them, type DGTtool.methodName with input/output arguments.
    %
    %   DGTtool Methods:
    %   Window utilities
    %      windowList              - Returns acceptable window names
    %      getWindow               - Compute window
    %      computeCanonicalDual    - Compute canonical dual window
    %      computeCanonicalTight   - Compute canonical tight window
    %      computeNumericalDiffWin - Compute numerical differential
    %      isdual                  - Check whether a pair of windows is dual
    %
    %   Zero-padding
    %      zeroPad               - Add zero at the end of signals
    %      extendSignalByZeroPad - Add zero at the end of signals
    %      zeroPadForFactorDGT   - Add zero at the end of signals
    %      sigLenForFactorDGT    - Compute required signal length
    %
    %
    %   See the associated demo file (demo.m) for explanation and usage.
    %   See also DGTtool, DGT, pinv, plot
    
    %   Author: Kohei Yatabe (2021)
    
    % [Memo]
    % 2020b is required because of pagemtimes.
    % 2019b might be enough for performing STFT/DGT (depending on settings).
    % 2019b is required because of arguments.
    % Signal Processing Toolbox is required because of buffer.
    % Signal Processing Toolbox is required for Slepian and Chebyshev windows.
    
    properties (Dependent)
        redundancy % Ratio of number of time-frequency bins to signal length.
    end
    
    properties
        shift   % Amount of shift of window
        FFTnum  % Number of frequency bins in time-frequency domain
        win     % Column vector of window
        dualWin % Column vector of dual window of win
    end
    
    properties (SetAccess = private)
        diffWin             % Differential of win
        sigLen              % Length of signal that is inputted last time
        isDual = false      % True when dualWin is dual of win
        isCanonical = false % True when dualWin is canonical dual of win
    end
    
    properties (Hidden = true)
        OLAindex
        factorIdx
        defConverter
        zeroPhaseConverter
        isNotCompDual = true
        isNotWinCalcInDGT = true
    end
    
    methods % main
        
        function obj = DGTtool(options)
            %Constructor: Create a DGTtool object with given parameters.
            %   Parameters of STFT/DGT can be set through Name-Value pairs.
            %
            %   Usage:
            %      F = DGTtool
            %      F = DGTtool(Name,Value)
            %      F = DGTtool(Name1,Value1,Name2,Value2,...)
            %
            %   Name of options:
            %      'windowShift'  (positive integer, default = 256)
            %      'windowLength' (positive integer, default = 2048)
            %      'FFTnum'       (positive integer, default = length(win))
            %      'windowName'   (char or string,   default = '4termC5Nuttall')
            %      'windowVector' (column vector)
            %
            %   Acceptable window name:
            %      'Hann'
            %      'Blackman'
            %      '3termC1Nuttall'
            %      '3termC3Nuttall'
            %      '4termC1Nuttall'
            %      '4termC3Nuttall'
            %      '4termC5Nuttall'
            %      'Gauss'
            %      'Slepian'
            %      'Chebyshev'
            %
            %   Note: If windowVector is given, windowLength and windowName are ignored.
            %
            %   See also windowList, setWindow, plotWin
            arguments
                options.windowShift  (1,1) {mustBePositive,mustBeInteger} = 256
                options.windowLength (1,1) {mustBePositive,mustBeInteger} = 2048
                options.FFTnum       (1,1) {mustBePositive,mustBeInteger}
                options.windowVector (:,1) double
                options.windowName   char
            end
            
            if isfield(options,'windowVector')
                obj.win = options.windowVector;
            elseif isfield(options,'windowName')
                [obj.win,obj.diffWin] = DGTtool.getWindow(options.windowLength,options.windowName);
            else
                [obj.win,obj.diffWin] = DGTtool.getWindow(options.windowLength,'4termC5Nuttall');
            end
            
            if isfield(options,'FFTnum')
                obj.FFTnum = options.FFTnum;
            else
                obj.FFTnum = length(obj.win);
            end
            
            obj.shift  = options.windowShift;
            obj.factorIdx = struct('c',[],'d',[],'p',[],'q',[],'wIdx',[],'xIdx',[],'cIdx',[]);
        end
        
        function varargout = subsref(obj,data)
            %SUBSREF: Shortcut for DGT.
            %   Parenthesis notation (operator form) calls DGT.
            %
            %   Usage:
            %      X = F(x)
            %
            %   See also DGT
            
            % Note: F(x) is shortcut of subsref(F,struct('type','()','subs',{{x}}))
            
            switch data(1).type
                case '()' % Call DGT
                    % Many of validations are skipped for speed.
                    if length(data) ~= 1
                        error 'DGT must be performed as F(x)'
                    end
                    if numel(data.subs) ~= 1
                        warning 'Only one input is allowed. Others will be ignored.'
                    end
                    
                    x = data.subs{1};
                    
                    [varargout{1:nargout}] = DGT(obj,x);
                    
                case '.' % implemented for access to properties and methods
                    switch length(data)
                        case 1
                            [varargout{1:nargout}] = obj.(data.subs);
                        case 2
                            [varargout{1:nargout}] = obj.(data(1).subs)(data(2).subs{:});
                        otherwise % not implemented for now
                            error 'Unexpected usage!'
                    end
                    
                otherwise % not implemented for now
                    error 'Unexpected usage!'
            end
        end
        
        function [X,f,t] = DGT(obj,x)
            %DGT: Compute spectrogram by DGT.
            %   Parenthesis notation can be used for shortcut.
            %   Normalized frequency f and sample index t can be returned.
            %
            %   Usage:
            %      X = F(x)
            %      [X,f,t] = F(x)
            %
            %   Note: The following notation gives the same result.
            %      X = F.DGT(x)
            %      X = DGT(F,x)
            %
            %   See also DGTtool, DGTtool/DGTtool, reassign
            
            if obj.FFTnum < length(obj.win)
                x = zeroPadForFactorAlg(obj,x,length(x));
            end
            
            if size(x,1) < length(obj.win)
                setFlag_winCalcInDGT(obj)
                obj.win = zeroPadForFactorAlg(obj,obj.win,length(obj.win));
                x = zeroPadForFactorAlg(obj,x,length(obj.win));
            end
            
            useFactorizationAlgorithm = size(x,1) == length(obj.win);
            
            if ~isequal(obj.sigLen,size(x,1))
                if useFactorizationAlgorithm
                    c = lcm(obj.shift,obj.FFTnum);
                else
                    c = obj.shift;
                end
                obj.sigLen = ceil(size(x,1)/c)*c;
            end
            
            if useFactorizationAlgorithm
                if ~isequal(obj.sigLen,size(x,1))
                    x = zeroPadForFactorAlg(obj,x,obj.sigLen);
                end
                if ~isequal(obj.sigLen,length(obj.win))
                    obj.win = zeroPadForFactorAlg(obj,obj.win,obj.sigLen);
                end
                if factorIdxMismatch(obj.shift,obj.FFTnum,size(x,1)/obj.shift,obj.factorIdx)
                    computeIndexForFactorAlg(obj,size(x,1)/obj.shift,length(x),size(x,2));
                end
                
                if nargout > 1
                    [X,f,t] = DGT_factorAlg(x,obj.win,obj.shift,obj.FFTnum,obj.factorIdx);
                else
                    X = DGT_factorAlg(x,obj.win,obj.shift,obj.FFTnum,obj.factorIdx);
                end
            else
                if nargout > 1
                    [X,f,t] = DGT_usualAlg(x,obj.win,obj.shift,obj.FFTnum,obj.sigLen);
                else
                    X = DGT_usualAlg(x,obj.win,obj.shift,obj.FFTnum,obj.sigLen);
                end
            end
        end
        
        function x = H(obj,X,sWin)
            %H: Inverse DGT.
            %   Synthesis window can be specified. This function is complex
            %   conjugate transpose of DGT if a window is not specified.
            %
            %   Usage:
            %      x = F.H(X)
            %      x = F.H(X,synthesisWindow)
            %
            %   See also pinv, makeWindowTight
            
            % Validation of input arguments is skipped for speed.
            % By default, this function uses analysis window.
            if ~exist('sWin','var')
                sWin = obj.win;
            end
            
            signalLength = size(X,2)*obj.shift;
            useFactorizationAlgorithm = signalLength <= length(sWin);
            
            if useFactorizationAlgorithm
                if length(sWin) ~= length(obj.win)
                    X = changeDGTdef(obj,X);
                end
                if length(sWin) ~= signalLength
                    sWin = zeroPadForFactorAlg(obj,sWin,length(sWin));
                    X = [X, X(:,1:(length(sWin)/obj.shift)-size(X,2),:)];
                end
                if isempty(obj.factorIdx) || sizeMismatch(obj.factorIdx,X) || factorIdxMismatch(obj.shift,obj.FFTnum,size(X,2),obj.factorIdx)
                    computeIndexForFactorAlg(obj,size(X,2),length(sWin),size(X,3));
                end
                
                x = invDGT_factorAlg(X,sWin,obj.shift,obj.FFTnum,obj.factorIdx);
            else
                if ~isequal(size(obj.OLAindex,1:3),[length(sWin) size(X,[2 3])])
                    computeIndexForOLA(obj,X,sWin);
                end
                
                x = invDGT_usualAlg(X,sWin,obj.FFTnum,obj.OLAindex);
            end
            
            obj.sigLen = size(x,1);
        end
        
        function x = pinv(obj,X)
            %PINV: Inverse DGT with perfect reconstruction.
            %   Signal is reconstructed using canonical dual window.
            %
            %   Usage:
            %      x = F.pinv(X)
            %
            %   See also DGTtool, DGT, H, computeCanonicalDual
            
            signalLength = size(X,2)*obj.shift;
            
            if isempty(obj.dualWin) || ~obj.isDual || signalLength < length(obj.dualWin)
                newLen = max(signalLength,length(obj.dualWin));
                newLen = DGTtool.sigLenForFactorDGT(newLen,obj.shift,obj.FFTnum);
                obj.dualWin = getCanonicalDualWin(obj,newLen);
            end
            if ~obj.isCanonical
                warning 'This is not pseudoinverse because the given dual window is not canonical. Consider using F.H(x,F.dualWin) for non-canonical dual.'
            end
            
            useFactorizationAlgorithm = signalLength <= length(obj.dualWin);
            
            x = H(obj,X,obj.dualWin);
            
            if useFactorizationAlgorithm && length(obj.win) ~= length(obj.dualWin)
                x = circshift(x,-length(obj.win)+obj.shift);
            end
        end
        
        function [reassignedS,f,t,X,IF,GD,dXdt,dXdf] = reassign(obj,x,epsilon)
            %REASSIGN: Compute reassigned spectrogram.
            %   Reassigned spectrogram is a sparse time-frequency representation.
            %
            %   Usage:
            %      rS = F.reassign(x)
            %      rS = F.reassign(x,epsilon)
            %      [rS,f,t,X,IF,GD,dXdt,dXdf] = F.reassign(x)
            %
            %   Input:
            %      x       - Time-domain signal (column vectors)
            %      epsilon - Small number used for avoiding zero-division
            %
            %   Output:
            %      rS   - Reassigned spectrogram (non-negative valued)
            %      f    - Normalized frequency (column vector)
            %      t    - Time indices (row vector)
            %      X    - Spectrogram (complex-valued)
            %      IF   - Instantaneous frequency (time-derivative of phase)
            %      GD   - Group delay (frequency-derivative of phase)
            %      dXdt - Time-derivative of spectrogram
            %      dXdf - Frequency-derivative of spectrogram
            %
            %   See also plotReassign
            arguments
                obj
                x       (:,:) double
                epsilon (1,1) double {mustBePositive} = 1e-12;
            end
            
            [~,maxIdx] = max(obj.win);
            if (maxIdx - length(obj.win)/2) > 1
                warning 'Window seems improper. Consider using F.setWindow.'
            end
            
            nWin = obj.win / sum(obj.win);
            
            cg = centerOfGravity(nWin);
            tRamp = (-length(nWin)/2:length(nWin)/2-1)' + mod(length(nWin),2)/2; % mod for odd/even cases
            tRamp = circshift(tRamp,floor(cg)-floor(length(nWin)/2)-1);
            tWin = tRamp .* nWin;
            
            if isempty(obj.diffWin)
                obj.sigLen = DGTtool.sigLenForFactorDGT(size(x,1),obj.shift,obj.FFTnum);
                dWin = DGTtool.computeNumericalDiffWin(nWin,obj.sigLen) * obj.FFTnum / length(obj.win);
                
                x    = DGTtool.extendSignalByZeroPad(x,obj.sigLen);
                tWin = DGTtool.extendSignalByZeroPad(tWin,obj.sigLen);
                nWin = DGTtool.extendSignalByZeroPad(nWin,obj.sigLen);
            else
                obj.sigLen = ceil(size(x,1)/obj.shift) * obj.shift;
                dWin = obj.diffWin * obj.FFTnum / length(obj.win) / sum(obj.win);
            end
            
            if size(x,1) < length(nWin) || size(x,1) < length(dWin)
                error 'Signal is shorter than window.'
            end
            useFactorizationAlgorithm = obj.sigLen == length(dWin);
            
            if useFactorizationAlgorithm
                if factorIdxMismatch(obj.shift,obj.FFTnum,size(x,1)/obj.shift,obj.factorIdx)
                    computeIndexForFactorAlg(obj,size(x,1)/obj.shift,length(x),size(x,2));
                end
                [X,f,t] = DGT_factorAlg(x,nWin,obj.shift,obj.FFTnum,obj.factorIdx);
                dXdt = DGT_factorAlg(x,dWin,obj.shift,obj.FFTnum,obj.factorIdx);
                dXdf = DGT_factorAlg(x,tWin,obj.shift,obj.FFTnum,obj.factorIdx);
            else
                [X,f,t] = DGT_usualAlg(x,nWin,obj.shift,obj.FFTnum,obj.sigLen);
                dXdt = DGT_usualAlg(x,dWin,obj.shift,obj.FFTnum,obj.sigLen);
                dXdf = DGT_usualAlg(x,tWin,obj.shift,obj.FFTnum,obj.sigLen);
            end
            
            S = abs(X).^2;
            Splus = S + epsilon;
            
            IF = -imag(dXdt.*conj(X)./Splus);
            GD =  real(dXdf.*conj(X)./Splus);
            
            fNew = f*obj.FFTnum + IF;
            tNew = t + GD;
            
            fNewIdx = round(fNew(:)) + 1;
            tNewIdx = round(tNew(:)/obj.shift);
            tNewIdx = mod(tNewIdx,size(X,2)) + 1;
            
            idx = (fNewIdx >= 1) & (fNewIdx <= size(X,1));
            chIdx = ones(size(X)) .* reshape(1:size(X,3),1,1,[]);
            reassignedS = accumarray([fNewIdx(idx) tNewIdx(idx) chIdx(idx)], S(idx), size(X,1:3));
            
            if nargin > 1
                t = mod(t,size(X,2)*obj.shift);
                [~,tRot] = min(t);
                X = circshift(X,-tRot+1,2);
                IF = circshift(IF,-tRot+1,2);
                GD = circshift(GD,-tRot+1,2);
                dXdt = circshift(dXdt,-tRot+1,2);
                dXdf = circshift(dXdf,-tRot+1,2);
                t = circshift(t,-tRot+1);
            end
        end
        
        function setWindow(obj,windowLength,windowName,varargin)
            %setWindow: Set a new window (current window is deleted).
            %   All input variables are directly passed to DGTtool.getWindow.
            %
            %   Usage:
            %      F.setWindow(windowLength,windowName)
            %      F.setWindow(___,Name,Value)
            %
            %   See also getWindow, plotWin
            
            [obj.win,obj.diffWin] = DGTtool.getWindow(windowLength,windowName,varargin{:});
        end
        
        function makeWindowTight(obj,signalLength)
            %makeWindowTight: Replace window by its canonical tight window.
            %   Canonical tight window makes F.H(X) and F.pinv(X) equal.
            %   Signal length must be specified whenever FFTnum < winLen.
            %
            %   Usage:
            %      F.makeWindowTight
            %      F.makeWindowTight(signalLength)
            %
            %   See also computeCanonicalTight, plotWin
            
            if ~exist('signalLength','var')
                signalLength = max(obj.sigLen,length(obj.win));
            else
                signalLength = max(signalLength,length(obj.win));
                obj.sigLen = ceil(signalLength/obj.shift)*obj.shift;
            end
            if isempty(signalLength)
                if obj.FFTnum < length(obj.win)
                    error 'Must specify signal length, e.g., F.makeWindowTight(length(x)).'
                else
                    signalLength = length(obj.win);
                end
            end
            obj.win = DGTtool.computeCanonicalTight(obj.win,obj.shift,obj.FFTnum,signalLength);
            obj.dualWin = obj.win;
        end
        
        function X = changeDGTdef(obj,X)
            %changeDGTdef: Phase is modified to change definition of DGT.
            %
            %   Usage:
            %      X = F.changeDGTdef(X)
            %
            %   See also undoChangeDGTdef
            
            if ~isequal(size(obj.defConverter),size(X,1:2))
                obj.defConverter = calculateDefConverter(X,obj.shift,obj.FFTnum);
            end
            X = X .* obj.defConverter;
        end
        
        function X = undoChangeDGTdef(obj,X)
            %undoChangeDGTdef: Reset phase modified by changeDGTdef.
            %
            %   Usage:
            %      X = F.undoChangeDGTdef(X)
            %
            %   See also changeDGTdef
            
            if ~isequal(size(obj.defConverter),size(X,1:2))
                error 'Parameters seem different.'
            end
            X = X .* conj(obj.defConverter);
        end
        
        function X = makeZeroPhase(obj,X)
            %makeZeroPhase: Phase is modified to remove linear phase component of window.
            %
            %   Usage:
            %      X = F.makeZeroPhase(X)
            %
            %   See also undoMakeZeroPhase
            
            if ~isequal(size(obj.zeroPhaseConverter,1),size(X,1))
                obj.zeroPhaseConverter = calculateZeroPhaseConverter(X,obj.FFTnum,obj.win);
            end
            X = X .* obj.zeroPhaseConverter;
        end
        
        function X = undoMakeZeroPhase(obj,X)
            %undoMakeZeroPhase: Reset phase modified by makeZeroPhase.
            %
            %   Usage:
            %      X = F.undoMakeZeroPhase(X)
            %
            %   See also makeZeroPhase
            
            if ~isequal(size(obj.zeroPhaseConverter,1),size(X,1))
                error 'Parameters seem different.'
            end
            X = X .* conj(obj.zeroPhaseConverter);
        end
    end
    
    methods % plot
        
        function plot(obj,x,fs,options)
            %PLOT: Compute spectrogram and display it.
            %   Signal and spectrum are also displayed.
            %
            %   Usage:
            %      F.plot(x)
            %      F.plot(x,fs)
            %      F.plot(___,Name,Value)
            %
            %   Options:
            %      'range'     (number, default = 80 [dB])
            %      'trunc'     (number, default = 0  [dB])
            %      'normalize' (true/false, default = true)
            %
            %   See also plotPhase, plotReassign, DGT
            arguments
                obj
                x   (:,:) {mustBeNumeric,mustBeSkinny}
                fs  (1,1) {mustBePositive} = 1
                options.range = 80
                options.trunc = 0
                options.normalize = true
            end
            
            if options.normalize
                normConst = sum(obj.win) / 2;
            else
                normConst = 1;
            end
            
            [X,f,t] = DGT(obj,x);
            X = X / normConst;
            
            t = mod(t,size(X,2)*obj.shift);
            [~,tRot] = min(t);
            X = circshift(X,-tRot+1,2);
            t = circshift(t,-tRot+1);
            
            if fs == 1
                fsPlot = fs;
                fUnit = '[periods/sample]';
                tUnit = '[samples]';
            else
                fsPlot = fs/1000;
                f = f*fsPlot;
                t = t/fs;
                fUnit = '[kHz]';
                tUnit = '[s]';
            end
            
            for n = 1:size(X,3)
                s = 20*log10(abs(fft(x(:,n)/normConst)));
                xLonger = buffer(x(:,n),size(X,2)*obj.shift);
                
                figure
                tiledlayout(10,14,'TileSpacing','none','Padding','compact')
                
                ax_s = nexttile(1,[8 2]);
                plot(s(1:floor(length(s)/2)+1),(0:floor(length(s)/2))/length(s)*fsPlot)
                xlim(max(s)-[options.range 0]-options.trunc)
                ylim([0 floor(length(s)/2)/length(s)*fsPlot])
                set(gca,'FontSize',10,'xDir','reverse')
                ylabel(['Frequency ' fUnit],'FontSize',12)
                
                ax_X = nexttile(3,[8 12]);
                imagesc(t,f,20*log10(abs(X(:,:,n))))
                axis xy off
                caxis(max(caxis)-[options.range 0]-options.trunc)
                xlim([0 size(x,1)-1]/fs)
                ylim([0 floor(length(s)/2)/length(s)*fsPlot])
                
                ax_x = nexttile(115,[2 12]);
                plot((0:length(xLonger)-1)/fs,xLonger)
                axis tight
                xlim([0 size(x,1)-1]/fs)
                set(gca,'FontSize',10)
                yticks(0)
                yline(0)
                xlabel(['Time ' tUnit],'FontSize',12)
                
                nexttile(127,[1 1])
                image([-1 1],ax_X.CLim,reshape(colormap,[],1,3))
                axis xy
                set(gca,'FontSize',10)
                xticks([])
                ylabel('Power [dB]','FontSize',11)
                
                linkaxes([ax_X ax_x],'x')
                linkaxes([ax_X ax_s],'y')
            end
        end
        
        function plotPhase(obj,x,fs,options)
            %plotPhase: Visualize phase.
            %   Color and brightness represent phase and magnitude, respectively.
            %   F.makeZeroPhase is applied for better visibility.
            %
            %   Usage:
            %      F.plotPhase(x)
            %      F.plotPhase(x,fs)
            %      F.plotPhase(___,Name,Value)
            %
            %   Options:
            %      'range'     (number, default = 40 [dB])
            %      'trunc'     (number, default = 15 [dB])
            %      'normalize' (true/false, default = true)
            %
            %   See also plot, plotReassign, makeZeroPhase
            arguments
                obj
                x   (:,:) {mustBeNumeric,mustBeSkinny}
                fs  (1,1) {mustBePositive} = 1
                options.range = 40
                options.trunc = 15
                options.normalize = true
            end
            
            if options.normalize
                normConst = sum(obj.win) / 2;
            else
                normConst = 1;
            end
            
            [X,f,t] = DGT(obj,x);
            X = X / normConst;
            X = makeZeroPhase(obj,X);
            
            t = mod(t,size(X,2)*obj.shift);
            [~,tRot] = min(t);
            X = circshift(X,-tRot+1,2);
            t = circshift(t,-tRot+1);
            
            if fs == 1
                fsPlot = fs;
                fUnit = '[periods/sample]';
                tUnit = '[samples]';
            else
                fsPlot = fs/1000;
                f = f*fsPlot;
                t = t/fs;
                fUnit = '[kHz]';
                tUnit = '[s]';
            end
            
            for n = 1:size(X,3)
                A = 20*log10(abs(X(:,:,n)));
                maxC = max(A(:)) - options.trunc;
                minC = maxC - options.range;
                
                A = rescale(A,'InputMin',minC,'InputMax',maxC);
                P = rescale(angle(X(:,:,n)),'InputMin',-pi,'InputMax',pi);
                C = hsv2rgb(cat(3,P,ones(size(X(:,:,n))),A));
                
                s = 20*log10(abs(fft(x(:,n)/normConst)));
                xLonger = buffer(x(:,n),size(X,2)*obj.shift);
                
                figure
                tiledlayout(10,14,'TileSpacing','none','Padding','compact')
                
                ax_s = nexttile(1,[8 2]);
                plot(s(1:floor(length(s)/2)+1),(0:floor(length(s)/2))/length(s)*fsPlot)
                xlim(max(s)-[options.range 0]-options.trunc)
                ylim([0 floor(length(s)/2)/length(s)*fsPlot])
                set(gca,'FontSize',10,'xDir','reverse')
                ylabel(['Frequency ' fUnit],'FontSize',12)
                
                ax_X = nexttile(3,[8 12]);
                image(t,f,C)
                axis xy off
                xlim([0 size(x,1)-1]/fs)
                ylim([0 floor(length(s)/2)/length(s)*fsPlot])
                
                ax_x = nexttile(115,[2 12]);
                plot((0:length(xLonger)-1)/fs,xLonger)
                axis tight
                xlim([0 size(x,1)-1]/fs)
                set(gca,'FontSize',10)
                yticks(0)
                yline(0)
                xlabel(['Time ' tUnit],'FontSize',12)
                
                nexttile(127,[1 1])
                image([-1 1],[minC maxC],hsv2rgb(cat(3,repmat(linspace(0,1,64),128,1),ones(128,64),repmat(linspace(0,1,128)',1,64))))
                axis xy
                set(gca,'FontSize',10)
                xticks(-1:1)
                xticklabels({'-\pi',0,'\pi'})
                xtickangle(0)
                xlabel('Phase [rad]','FontSize',11)
                ylabel('Power [dB]','FontSize',11)
                
                linkaxes([ax_X ax_x],'x')
                linkaxes([ax_X ax_s],'y')
            end
        end
        
        function plotReassign(obj,x,fs,options)
            %plotReassign: Compute reassigned spectrogram and display it.
            %   Resolution of reassigned spectrogram depends on shift and FFTnum.
            %
            %   Usage:
            %      F.plotReassign(x)
            %      F.plotReassign(x,fs)
            %      F.plotReassign(___,Name,Value)
            %
            %   Options:
            %      'range'     (number, default = 100 [dB])
            %      'trunc'     (number, default = 0   [dB])
            %      'normalize' (true/false, default = true)
            %
            %   See also plot, plotPhase, reassign
            arguments
                obj
                x   (:,:) {mustBeNumeric,mustBeSkinny}
                fs  (1,1) {mustBePositive} = 1
                options.range = 100
                options.trunc = 0
                options.normalize = true
            end
            
            [X,f,t] = reassign(obj,x);
            
            if options.normalize
                normConst = sum(obj.win);
            else
                normConst = 1;
                X = X * sum(obj.win);
            end
            
            t = mod(t,size(X,2)*obj.shift);
            [~,tRot] = min(t);
            X = circshift(X,-tRot+1,2);
            t = circshift(t,-tRot+1);
            
            if fs == 1
                fsPlot = fs;
                fUnit = '[periods/sample]';
                tUnit = '[samples]';
            else
                fsPlot = fs/1000;
                f = f*fsPlot;
                t = t/fs;
                fUnit = '[kHz]';
                tUnit = '[s]';
            end
            
            for n = 1:size(X,3)
                s = 20*log10(abs(fft(x(:,n)/normConst)));
                xLonger = buffer(x(:,n),size(X,2)*obj.shift);
                
                figure
                tiledlayout(10,14,'TileSpacing','none','Padding','compact')
                
                ax_s = nexttile(1,[8 2]);
                plot(s(1:floor(length(s)/2)+1),(0:floor(length(s)/2))/length(s)*fsPlot)
                xlim(max(s)-[options.range 0]-options.trunc)
                ylim([0 floor(length(s)/2)/length(s)*fsPlot])
                set(gca,'FontSize',10,'xDir','reverse')
                ylabel(['Frequency ' fUnit],'FontSize',12)
                
                ax_X = nexttile(3,[8 12]);
                imagesc(t,f,20*log10(abs(X(:,:,n))))
                axis xy off
                caxis(max(caxis)-[options.range 0]-options.trunc)
                xlim([0 size(x,1)-1]/fs)
                ylim([0 floor(length(s)/2)/length(s)*fsPlot])
                
                ax_x = nexttile(115,[2 12]);
                plot((0:length(xLonger)-1)/fs,xLonger)
                axis tight
                xlim([0 size(x,1)-1]/fs)
                set(gca,'FontSize',10)
                yticks(0)
                yline(0)
                xlabel(['Time ' tUnit],'FontSize',12)
                
                nexttile(127,[1 1])
                image([-1 1],ax_X.CLim,reshape(colormap,[],1,3))
                axis xy
                set(gca,'FontSize',10)
                xticks([])
                ylabel('Power [dB]','FontSize',11)
                
                linkaxes([ax_X ax_x],'x')
                linkaxes([ax_X ax_s],'y')
            end
        end
        
        function plotWin(obj)
            %plotWin: Display currently available windows.
            %   Some windows may not appear if they are not calculated yet.
            %
            %   Usage:
            %      F.plotWin
            %
            %   See also setWindow, getWindow
            
            h = figure;
            h.Position(4) = h.Position(4)/2;
            tiledlayout(1,3,'TileSpacing','compact','Padding','compact')
            
            nexttile
            windowStylePlot(obj.win)
            title('window','FontSize',12)
            
            nexttile
            windowStylePlot(obj.dualWin)
            title('dual window','FontSize',12)
            
            nexttile
            windowStylePlot(obj.diffWin)
            title('differential window','FontSize',12)
        end
        
    end
    
    methods (Static)
        
        function winList = windowList
            %winList: List of window names acceptable in DGTtool.
            %
            %   Usage:
            %      c = DGTtool.windowList
            %
            %   See also getWindow, DGTtool, DGTtool/DGTtool
            
            winList = {
                'ForKitamuraHamming'
                'Hann'
                'Blackman'
                '3termC1Nuttall'
                '3termC3Nuttall'
                '4termC1Nuttall'
                '4termC3Nuttall'
                '4termC5Nuttall'
                'Gauss'
                'Slepian'
                'Chebyshev'
                };
        end
        
        function [win,diffWin] = getWindow(windowLength,windowName,options)
            %getWindow: Compute specified window.
            %   Windows are returned as column vectors.
            %
            %   Usage:
            %      win = DGTtool.getWindow(windowLength,windowName)
            %      [win,diffWin] = DGTtool.getWindow(___)
            %
            %   Acceptable window name (shortcut):
            %      'Hann'           ('h')
            %      'Blackman'       ('b')
            %      '3termC1Nuttall' ('3termC1')
            %      '3termC3Nuttall' ('3termC3')
            %      '4termC1Nuttall' ('4termC1')
            %      '4termC3Nuttall' ('4termC3')
            %      '4termC5Nuttall' ('4termC5')
            %      'Gauss'          ('g')
            %      'Slepian'        ('s')
            %      'Chebyshev'      ('c')
            %
            %   For 'Gauss', 'Slepian' and 'Chebyshev', options can be set.
            %      win = DGTtool.getWindow(___,Name,Value)
            %
            %   Options (width parameters):
            %      'Gauss'     (number, default = 0.02)
            %      'Slepian'   (number, default = 12)
            %      'Chebyshev' (number, default = 340 [dB])
            %
            %   see also DGTtool, DGTtool/DGTtool, setWindow
            arguments
                windowLength (1,1) {mustBePositive,mustBeInteger}
                windowName   char
                options.Gauss     (1,1) {mustBePositive} = 0.02
                options.Slepian   (1,1) {mustBePositive} = 12
                options.Chebyshev (1,1) {mustBePositive} = 340
            end
            
            winList = DGTtool.windowList;
            winName = validatestring(windowName,winList);
            
            isodd = mod(windowLength,2); % 1 (if odd) or 0 (if even)
            K = windowLength + isodd;    % always even
            
            if ismember(winName,winList(1:8)) % cosine window case
                
                switch winName
                    case 'ForKitamuraHamming'
                        c = [25; 21]./46;
                    case 'Hann'
                        c = [0.5; 0.5];
                    case 'Blackman'
                        c = [0.42; 0.5; 0.08];
                    case '3termC1Nuttall'
                        c = [0.40897; 0.5; 0.09103];
                    case '3termC3Nuttall'
                        c = [3; 4; 1]/8;
                    case '4termC1Nuttall'
                        c = [0.355768; 0.487396; 0.144232; 0.012604];
                    case '4termC3Nuttall'
                        c = [0.338946; 0.481973; 0.161054; 0.018027];
                    case '4termC5Nuttall'
                        c = [10; 15; 6; 1]/32;
                end
                
                t = (-(K/2-1):(K/2-1))' / K; % always odd (-0.5 < t < 0.5)
                win = cos(2*pi*t.*(0:length(c)-1)) * c;
                
                if nargout == 2
                    diffC = (0:length(c)-1)' .* c;
                    diffWin = -sin(2*pi*t.*(0:length(c)-1)) * diffC;
                end
                
            else % other case
                
                switch winName
                    case 'Gauss'
                        t = (-(K/2-1):(K/2-1))' / K;
                        win = exp(-pi*t.^2/options.Gauss);
                    case 'Slepian'
                        win = dpss(K-1,options.Slepian,1);    % Signal Processing Toolbox
                        win = win / max(win);
                    case 'Chebyshev'
                        win = chebwin(K-1,options.Chebyshev); % Signal Processing Toolbox
                end
                
                if nargout == 2
                    if isequal(winName,'Gauss')
                        diffWin = -exp(-pi*t.^2/options.Gauss).*t/options.Gauss;
                    elseif win(1) < eps * length(win)
                        diffWin = DGTtool.computeNumericalDiffWin(win,2^nextpow2(2*length(win)));
                        diffWin = diffWin(1:length(win));
                    else
                        diffWin = [];
                    end
                end
                
            end
            
            if ~isodd % even
                win = [0; win];
                if exist('diffWin','var') && ~isempty(diffWin)
                    diffWin = [0; diffWin];
                end
            end
        end
        
        function dualWin = computeCanonicalDual(win,shift,FFTnum,sigLen)
            %computeCanonicalDual: Compute canonical dual of given window.
            %   Canonical dual ensures perfect reconstruction of signal.
            %   Signal length must be specified whenever FFTnum < winLen.
            %
            %   Usage:
            %      dualWin = DGTtool.computeCanonicalDual(win,shift,FFTnum)
            %      dualWin = DGTtool.computeCanonicalDual(win,shift,FFTnum,sigLen)
            %
            %   See also pinv
            arguments
                win    (:,1) {mustBeNumeric}
                shift  (1,1) {mustBePositive,mustBeInteger}
                FFTnum (1,1) {mustBePositive,mustBeInteger}
                sigLen (1,1) {mustBeInteger} = -1
            end
            
            if FFTnum >= length(win)
                dualWin = buffer(win,shift);
                dualWin = dualWin ./ sum(abs(dualWin).^2,2);
                dualWin = reshape(dualWin(1:length(win)),[],1);
            else
                if sigLen < length(win)
                    error 'sigLen must be specified whenever FFTnum < winLen.'
                end
                sigLen = DGTtool.sigLenForFactorDGT(sigLen,shift,FFTnum);
                [c,~,p,q,d,k,l,r,s] = getConstantsForFacAlg(shift,FFTnum,sigLen/shift);
                idx = getWinIdxForFacAlg(k,l,r,s,c,p,q,d);
                
                win = DGTtool.extendSignalByZeroPad(win,sigLen);
                phi = fft(win(idx),[],4);
                
                S = pagemtimes(phi,'none',phi,'ctranspose');
                for m = 1:size(phi,4)
                    for n = 1:size(phi,3)
                        phi(:,:,n,m) = S(:,:,n,m) \ phi(:,:,n,m);
                    end
                end
                
                phi = ifft(phi,[],4,'symmetric');
                dualWin = zeros(size(win));
                dualWin(idx) = phi;
            end
        end
        
        function tightWin = computeCanonicalTight(win,shift,FFTnum,sigLen)
            %computeCanonicalTight: Compute canonical tight window of given window.
            %   Using canonical tight window for both analysis and synthesis
            %   results in perfect reconstruction of signal.
            %   Signal length must be specified whenever FFTnum < winLen.
            %
            %   Usage:
            %      tightWin = DGTtool.computeCanonicalTight(win,shift,FFTnum)
            %      tightWin = DGTtool.computeCanonicalTight(win,shift,FFTnum,sigLen)
            %
            %   See also makeWindowTight
            arguments
                win    (:,1) {mustBeNumeric}
                shift  (1,1) {mustBePositive,mustBeInteger}
                FFTnum (1,1) {mustBePositive,mustBeInteger}
                sigLen (1,1) {mustBeInteger} = -1
            end
            
            if FFTnum >= length(win)
                tightWin = buffer(win,shift);
                tightWin = tightWin ./ sqrt(sum(abs(tightWin).^2,2));
                tightWin = reshape(tightWin(1:length(win)),[],1);
            else
                if sigLen < length(win)
                    error 'sigLen (>= winLen) must be specified whenever FFTnum < winLen.'
                end
                sigLen = DGTtool.sigLenForFactorDGT(sigLen,shift,FFTnum);
                [c,~,p,q,d,k,l,r,s] = getConstantsForFacAlg(shift,FFTnum,sigLen/shift);
                idx = getWinIdxForFacAlg(k,l,r,s,c,p,q,d);
                
                win = DGTtool.extendSignalByZeroPad(win,sigLen);
                phi = fft(win(idx),[],4);
                
                if p == 1
                    phi = phi ./ sqrt(sum(abs(phi).^2,[1 2]));
                else
                    for m = 1:size(phi,4)
                        for n = 1:size(phi,3)
                            [U,~,V] = svd(phi(:,:,n,m),'econ');
                            phi(:,:,n,m) = U*V';
                        end
                    end
                end
                
                phi = ifft(phi,[],4,'symmetric');
                tightWin = zeros(size(win));
                tightWin(idx) = phi;
            end
        end
        
        function diffWin = computeNumericalDiffWin(win,sigLen)
            %computeNumericalDiffWin: Compute numerical differential.
            %   Differential of window is necessary for reassignment.
            %   If sigLen is specified, output is extended to sigLen.
            %
            %   Usage:
            %      diffWin = DGTtool.computeNumericalDiffWin(win)
            %      diffWin = DGTtool.computeNumericalDiffWin(win,sigLen)
            %
            %   See also reassign, plotReassign
            arguments
                win    (:,1) {mustBeNumeric}
                sigLen (1,1) {mustBePositive,mustBeInteger} = length(win)
            end
            
            origWinLen = length(win);
            win = DGTtool.extendSignalByZeroPad(win,sigLen);
            
            M = floor((length(win)-1)/2);
            fftIdx = ifftshift([zeros(mod(length(win)-1,2)),-M:M]); % mod for odd/even cases
            fftIdx = fftIdx(:)*origWinLen/length(win);
            
            diffWin = ifft(1i*fftIdx.*fft(win),'symmetric'); % spectral method
        end
        
        function [tf,reconstErrorBound] = isdual(win,dualWin,shift,FFTnum)
            %ISDUAL: Check whether two windows are dual of each other.
            %   Dual window pair can perfectly reconstruct signal.
            %   This function can return upper bound of (relative) reconstruction error.
            %
            %   Usage:
            %      trueFalse = DGTtool.isdual(win1,win2,shift,FFTnum)
            %      [trueFalse,relReconstErrorBound] = DGTtool.isdual(___)
            %
            %   See also computeCanonicalDual, computeCanonicalTight
            arguments
                win     (:,1) {mustBeNumeric}
                dualWin (:,1) {mustBeNumeric}
                shift   (1,1) {mustBePositive,mustBeInteger}
                FFTnum  (1,1) {mustBePositive,mustBeInteger}
            end
            
            maxLen  = max(length(win),length(dualWin));
            win     = DGTtool.zeroPadForFactorDGT(win,    shift,FFTnum,maxLen);
            dualWin = DGTtool.zeroPadForFactorDGT(dualWin,shift,FFTnum,maxLen);
            idx = getIndexForFactorAlg(FFTnum,shift,length(win)/FFTnum,length(win),1);
            
            WexlerRaz = DGT_factorAlg(dualWin,win,FFTnum,shift,idx);
            WexlerRaz(1,1) = WexlerRaz(1,1) - shift;
            absWR = abs(WexlerRaz)/FFTnum;
            reconstErrorBound = sum(absWR,'all') + sum(absWR(2:end-1+mod(shift,2),:),'all');
            
            tf = reconstErrorBound < eps(max(abs(win))) * length(win);
        end
        
        function y = zeroPad(x,zeroNum)
            %zeroPad: Zero-padding by specifying number of added zero.
            %   This function allows multi-channel signal.
            %
            %   Usage:
            %      x = DGTtool.zeroPad(x,zeroNum)
            %
            %   See also extendSignalByZeroPad, zeroPadForFactorDGT
            arguments
                x       (:,:) {mustBeNumeric,mustBeSkinny}
                zeroNum (1,1) {mustBeNonnegative,mustBeInteger}
            end
            
            newLength = size(x,1) + zeroNum;
            y = zeros(newLength,size(x,2));
            for n = 1:size(x,2)
                y(:,n) = buffer(x(:,n),newLength);
            end
        end
        
        function y = extendSignalByZeroPad(x,outputLength)
            %extendSignalByZeroPad: Zero-padding by specifying length of output.
            %   This function allows multi-channel signal.
            %
            %   Usage:
            %      x = DGTtool.extendSignalByZeroPad(x,outputLength)
            %
            %   See also zeroPad, zeroPadForFactorDGT
            
            if outputLength < size(x,1)
                error 'Must satisfy outputLength >= size(x,1)'
            end
            y = DGTtool.zeroPad(x,outputLength-size(x,1));
        end
        
        function y = zeroPadForFactorDGT(x,shift,FFTnum,minLen)
            %zeroPadForFactorDGT: Zero-padding for factorization algorithm.
            %   Number of zeros is determined to satisfy L = aN = bM.
            %   Signal can be further extended by specifying lower bound of length.
            %
            %   Usage:
            %      x = DGTtool.zeroPadForFactorDGT(x,shift,FFTnum)
            %      x = DGTtool.zeroPadForFactorDGT(x,shift,FFTnum,minLen)
            %
            %   See also sigLenForFactorDGT, zeroPad, extendSignalByZeroPad
            
            % This zero-padding is necessary for using DGT_factorAlg.
            % The result may be unnecessarily long for DGT_usualAlg.
            arguments
                x      (:,:) {mustBeNumeric,mustBeSkinny}
                shift  (1,1) {mustBePositive,mustBeInteger}
                FFTnum (1,1) {mustBePositive,mustBeInteger}
                minLen (1,1) {mustBePositive,mustBeInteger} = size(x,1)
            end
            
            if size(x,1) > minLen
                error 'Must satisfy minLen >= size(x,1)'
            end
            
            newLength = DGTtool.sigLenForFactorDGT(minLen,shift,FFTnum);
            y = DGTtool.extendSignalByZeroPad(x,newLength);
        end
        
        function newLength = sigLenForFactorDGT(xLen,shift,FFTnum)
            %sigLenForFactorDGT: Compute signal length necessary for factorization algorithm.
            %   Obtained length L satisfies L = aN = bM.
            %
            %   Usage:
            %      newLen = DGTtool.sigLenForFactorDGT(sigLen,shift,FFTnum)
            %
            %   See also zeroPadForFactorDGT, extendSignalByZeroPad
            arguments
                xLen   (1,1) {mustBePositive,mustBeInteger}
                shift  (1,1) {mustBePositive,mustBeInteger}
                FFTnum (1,1) {mustBePositive,mustBeInteger}
            end
            
            c = lcm(shift,FFTnum);
            newLength = ceil(xLen/c) * c;
        end
    end
    
    methods % set and get (property access methods)
        
        function set.shift(obj,shift)
            arguments
                obj
                shift (1,1) {mustBePositive,mustBeInteger}
            end
            obj.shift = shift;
            checkRedundancy(obj)
            checkWindowAndShift(obj)
            deleteInternalParameters(obj)
        end
        
        function set.FFTnum(obj,FFTnum)
            arguments
                obj
                FFTnum (1,1) {mustBePositive,mustBeInteger}
            end
            obj.FFTnum = FFTnum;
            checkRedundancy(obj)
            deleteInternalParameters(obj)
        end
        
        function set.redundancy(obj,redundancy)
            obj.FFTnum = ceil(obj.shift * redundancy);
            deleteInternalParameters(obj)
        end
        
        function redundancy = get.redundancy(obj)
            redundancy = obj.FFTnum / obj.shift;
        end
        
        function set.win(obj,win)
            arguments
                obj
                win {mustBeSkinny}
            end
            obj.win = win;
            deleteDualWin(obj)
            deleteDiffWin(obj)
        end
        
        function set.dualWin(obj,win)
            arguments
                obj
                win {mustBeSkinny}
            end
            obj.dualWin = win;
            if ~isempty(win)
                checkDualWin(obj)
            end
        end
    end
    
    methods (Hidden) % internal functions
        
        function computeIndexForFactorAlg(obj,segNum,sigLen,chNum)
            obj.factorIdx = getIndexForFactorAlg(obj.shift,obj.FFTnum,segNum,sigLen,chNum);
            disp 'index computed (factor DGT)'
        end
        
        function computeIndexForOLA(obj,X,win)
            obj.OLAindex = getIndexForOLA(X,win,obj.shift);
            disp 'index computed (OLA)'
        end
        
        function x = zeroPadForFactorAlg(obj,x,newLength)
            x = DGTtool.zeroPadForFactorDGT(x,obj.shift,obj.FFTnum,newLength);
        end
        
        function checkRedundancy(obj)
            if ~isempty(obj.shift) && ~isempty(obj.FFTnum)
                if obj.redundancy <= 1
                    error(['FFTnum > shift is required for reconstruction! ' ...
                        '(FFTnum = ' num2str(obj.FFTnum) ', shift = ' num2str(obj.shift) ')'])
                end
            end
        end
        
        function checkWindowAndShift(obj)
            if ~isempty(obj.shift) && ~isempty(obj.win)
                if obj.shift >= nnz(obj.win)
                    error(['winLen > shift is required for reconstruction! ' ...
                        '(winLen = ' num2str(nnz(obj.win)) ', shift = ' num2str(obj.shift) ')'])
                end
            end
        end
        
        function deleteInternalParameters(obj)
            if ~isempty(obj.dualWin)
                deleteDualWin(obj)
            end
            if ~isempty(obj.sigLen)
                obj.sigLen = [];
            end
            if ~isempty(obj.OLAindex)
                obj.OLAindex = [];
            end
            if ~isempty(obj.factorIdx)
                obj.factorIdx = [];
            end
            if ~isempty(obj.defConverter)
                obj.defConverter = [];
            end
            if ~isempty(obj.zeroPhaseConverter)
                obj.zeroPhaseConverter = [];
            end
        end
        
        function deleteDualWin(obj)
            if ~isempty(obj.dualWin) && obj.isNotWinCalcInDGT
                obj.dualWin = [];
                obj.isDual = false;
                obj.isCanonical = false;
                disp 'dualWin deleted'
            elseif ~obj.isNotWinCalcInDGT
                obj.isNotWinCalcInDGT = true;
            end
        end
        
        function deleteDiffWin(obj)
            if ~isempty(obj.diffWin)
                obj.diffWin = [];
                disp 'diffWin deleted'
            end
        end
        
        function canonicalDual = getCanonicalDualWin(obj,sigLen)
            canonicalDual = DGTtool.computeCanonicalDual(obj.win,obj.shift,obj.FFTnum,sigLen);
            obj.isNotCompDual = false;
            obj.isDual = true;
            obj.isCanonical = true;
        end
        
        function checkDualWin(obj)
            if obj.isNotCompDual
                obj.isDual = DGTtool.isdual(obj.win,obj.dualWin,obj.shift,obj.FFTnum);
                setCanonicalFlag(obj)
                disp 'dual window checked'
            else
                obj.isNotCompDual = true;
            end
        end
        
        function setCanonicalFlag(obj)
            if obj.isDual
                canonicalDual = DGTtool.computeCanonicalDual(obj.win,obj.shift,obj.FFTnum,length(obj.dualWin));
                winError = norm(obj.dualWin - canonicalDual);
                if winError < eps(max(abs(canonicalDual))) * length(canonicalDual)
                    obj.isCanonical = true;
                else
                    obj.isCanonical = false;
                end
            else
                obj.isCanonical = false;
            end
        end
        
        function setFlag_winCalcInDGT(obj)
            obj.isNotWinCalcInDGT = false;
        end
    end
end



% ------------------------------------------------------------
% DGT/IDGT
% ------------------------------------------------------------

function [X,f,t] = DGT_usualAlg(x,win,shift,FFTnum,paddedSiglen)
% paddedSiglen = ceil(size(x,1)/shift)*shift;
% This computation is performed outside for speed.
% paddedSiglen must be an integer multiple of shift (otherwise, error occurs at zeros).

winLen = length(win);
segNum = paddedSiglen / shift; % must be integer
overlap = winLen - shift;
rotNum = overlap - (paddedSiglen - size(x,1));

wx = zeros(winLen,segNum,size(x,2));
for n = 1:size(x,2)
    wx(:,:,n) = buffer(x(:,n),winLen,overlap,buffer(x(end-rotNum+1:end,n),overlap));
end
wx = win .* wx;

if winLen <= FFTnum
    Xfull = fft(wx,FFTnum);
else
    Xfull = zeros(FFTnum,segNum,size(x,2));
    for n = 1:size(wx,3)
        for m = 1:size(wx,2)
            Xfull(:,m,n) = fft(wrapData(wx(:,m,n),FFTnum));
        end
    end
end
X = Xfull(1:floor(FFTnum/2)+1,:,:);

if nargout > 1
    f = (0:size(X,1)-1)'/FFTnum;
    cg = centerOfGravity(win) - 1;
    t = (0:size(X,2)-1)*shift + cg - overlap;
end
end

function x = invDGT_usualAlg(X,win,FFTnum,OLAidx)
% Over-lap add (OLA) is performed by accumarray and OLAidx.

wx = ifft([X; zeros(FFTnum-size(X,1),size(X,2),size(X,3))],'symmetric');

if size(wx,1) == length(win)
    wx = win .* wx;
elseif size(wx,1) > length(win)
    wx = win .* wx(1:length(win),:,:);
else
    xrep = repmat(wx,ceil(length(win)/size(wx,1)),1,1);
    wx = win .* xrep(1:length(win),:,:);
end

x = accumarray(OLAidx(:),wx(:)); % OLA
x = reshape(x,[],size(wx,3));
end

function OLAidx = getIndexForOLA(X,win,shift)
winLen = length(win);
overlap = winLen - shift;

signalLength = size(X,2) * shift;
idx = uint32(1:signalLength);
OLAindex = buffer(idx,winLen,overlap,idx(end-overlap+1:end));
OLAidx = OLAindex + reshape(uint32(signalLength*(0:size(X,3)-1)),1,1,[]);
end



% ------------------------------------------------------------
% DGT/IDGT (Factorization Algorithm)
% ------------------------------------------------------------

function [X,f,t] = DGT_factorAlg(x,win,shift,FFTnum,in)
% Peter L. Sondergaard
% Efficient algorithms for the discrete Gabor transform with a long FIR window
% Journal of Fourier Analysis and Applications, vol.18, pp.456-470, 2012.

if isNotInteger(size(x,1)/lcm(shift,FFTnum))
    error 'Length of input signal L must satisfy L = aN = bM.'
end

phi = fft(win(in.wIdx),[],4);

psi = fft(x(in.xIdx),[],4);
C = pagemtimes(phi,'ctranspose',psi,'none');
C = ifft(C,[],4,'symmetric');
wx = C(in.cIdx);

if size(in.cIdx,2) == 1 && size(in.cIdx,3) == 1
    Xfull = fft(wx(:));
    X = Xfull(1:floor(FFTnum/2)+1);
else
    Xfull = fft(wx);
    X = Xfull(1:floor(FFTnum/2)+1,:,:);
end

if nargout > 1
    f = (0:size(X,1)-1)'/FFTnum;
    cg = centerOfGravity(win) - 1;
    t = mod((0:size(X,2)-1)*shift + cg, size(x,1));
end
end

function x = invDGT_factorAlg(X,win,shift,FFTnum,in)
% Peter L. Sondergaard
% Efficient algorithms for the discrete Gabor transform with a long FIR window
% Journal of Fourier Analysis and Applications, vol.18, pp.456-470, 2012.

chNum = size(X,3);
sigLen = size(X,2) * shift;

x = zeros(sigLen,size(X,3));
C = zeros([in.q in.q*chNum in.c in.d]);

phi = fft(win(in.wIdx),[],4);

wx = ifft([X; zeros(FFTnum-size(X,1),size(X,2),size(X,3))],'symmetric');
C(in.cIdx) = wx;
C = fft(C,[],4);

psi = pagemtimes(phi,C);
psi = ifft(psi,[],4,'symmetric');

x(in.xIdx) = psi;
end

function factorIdx = getIndexForFactorAlg(shift,FFTnum,segNum,sigLen,chNum)
[c,ha,p,q,d,k,l,r,s] = getConstantsForFacAlg(shift,FFTnum,segNum);
[wIdx,xIdx,cIdx] = getIndicesForFacAlg(k,l,r,s,c,p,q,d,ha,sigLen,shift,FFTnum,segNum,chNum);
factorIdx = struct('c',c,'d',d,'p',p,'q',q,'wIdx',wIdx,'xIdx',xIdx,'cIdx',cIdx);
end

function [c,ha,p,q,d,k,l,r,s] = getConstantsForFacAlg(shift,FFTnum,segNum)
[c,ha,~] = gcd(shift,FFTnum);
c = int32(c);

p = shift/c;
q = FFTnum/c;
d = segNum/q;

k = reshape(0:p-1, [],1,1,1);
l = reshape(0:q-1, 1,[],1,1);
r = reshape(1:c  , 1,1,[],1);
s = reshape(0:d-1, 1,1,1,[]);
end

function [wIdx,xIdx,cIdx] = getIndicesForFacAlg(k,l,r,s,c,p,q,d,ha,L,shift,FFTnum,segNum,chNum)
wIdx = getWinIdxForFacAlg(k,l,r,s,c,p,q,d);
xIdx = getSigIdxForFacAlg(k,l,r,s,p,q,ha,L,shift,FFTnum,chNum);
cIdx = getSpecIdxForFacAlg(c,q,d,ha,FFTnum,segNum,chNum);
end

function idx = getWinIdxForFacAlg(k,l,r,s,c,p,q,d)
idx = r + c*mod(k*q - l*p + s*(p*q), d*p*q);
end

function idx = getSigIdxForFacAlg(k,l,r,s,p,q,ha,L,shift,FFTnum,chNum)
idx = r + mod(k*FFTnum + s*(p*FFTnum) + l*(ha*shift), L);
if chNum > 1
    idx = repmat(idx,1,chNum) + int32(repelem(L*(0:chNum-1),q));
end
end

function idx = getSpecIdxForFacAlg(c,q,d,ha,FFTnum,segNum,chNum)
sizeC = [q q*chNum c d];
idx = reshape(permute(reshape(int32(1:prod(sizeC)),sizeC),[3 1 4 2]),c,q*d,q,chNum);
for n = 1:q
    idx(:,:,n,:) = circshift(idx(:,:,n,:),(n-1)*ha,2);
end
idx = reshape(permute(idx,[1 3 2 4]),FFTnum,segNum,chNum);
end



% ------------------------------------------------------------
% Phase manipulation
% ------------------------------------------------------------

function defConverter = calculateDefConverter(X,shift,FFTnum)
f = (0:size(X,1)-1)';
t = (0:size(X,2)-1) * shift;
idx = mod(f*t,FFTnum);
defConverter = exp(-2i*pi*idx/FFTnum);
end

function zeroPhaseConverter = calculateZeroPhaseConverter(X,FFTnum,win)
rotNum = centerOfGravity(win) - 1;
f = (0:size(X,1)-1)';
zeroPhaseConverter = exp(2i*pi*f*rotNum/FFTnum);
end



% ------------------------------------------------------------
% Helper functions
% ------------------------------------------------------------

function y = wrapData(x,M)
%wrapData: Simplified version of datawrap in Signal Processing Toolbox.
%   Validation of input arguments is skipped for speed.
%   x must be a vector (:,1) or (1,:)
%   M must be a positive integer (1,1)
y = sum(buffer(x,M),2);
end

function cg = centerOfGravity(x)
arguments
    x (:,1) double
end
cg = sum((1:length(x))'.*x) / sum(x);
end

function windowStylePlot(win)
plot(win,'linewidth',2)
xticks([])
yticks([])
yline(0)
box on
if ~isempty(win)
    xlim([1 length(win)])
    ylim(max(abs(win))*1.4*[-1 1])
end
end

function tf = isNotInteger(x)
tf = x ~= round(x);
end

function tf = sizeMismatch(in,X)
cIdxSize = size(in.cIdx);
cIdxSize(1) = floor(cIdxSize(1)/2) + 1;
tf = ~isequal(cIdxSize,size(X));
end

function tf = factorIdxMismatch(shift,FFTnum,segNum,in)
[c,~,~] = gcd(shift,FFTnum);
c = int32(c);
p = shift/c;
q = FFTnum/c;
d = segNum/q;

checkTF = zeros(5,1);
checkTF(1) = ~isequal(in.c,c);
checkTF(2) = ~isequal(in.d,d);
checkTF(3) = ~isequal(in.p,p);
checkTF(4) = ~isequal(in.q,q);
tf = any(checkTF);
end

function mustBeSkinny(x)
% The first dimension of input matrix x must be greater than the second dimension.
if size(x,1) ~= 0
    if size(x,1) < size(x,2)
        eidType = 'mustBeSkinny:notSkinny';
        msgType = 'Input matrix must be skinny (tall), i.e., size(x,1) >= size(x,2)';
        throwAsCaller(MException(eidType,msgType))
    end
end
end