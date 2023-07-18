function [simTS,simMarkov,simGauss] = ...
    sim_bold_markov(targSz,targTR,targFq,hrfLen,targMarkov,targGausMStd)

if nargin < 4
    hrfLen = 5 ;
end

if nargin<5
    targMarkov = [0.95 0.05 ; 0.2 0.8] ;
end

if nargin<6
    targGausMStd = [0 0] ; % no effect
end

% need econometric toolbox
march = dtmc(targMarkov) ;

% generate random data
hrf = getcanonicalhrf(hrfLen,targTR) ;
% simulation length pre-convolve
simLen = targSz(1)+length(hrf)-1 ; 

% simulate markov signals
rs_markov = nan(simLen,targSz(2)) ;
for idx = 1:targSz(2)
    rs_markov(:,idx) = simulate(march,simLen-1) - 1 ;
end

% randomize gaussian noise
rs_gauss = targGausMStd(1) + ...
    (targGausMStd(2) .* randn(simLen,targSz(2))) ;

% add the two together
rs_ = rs_gauss + rs_markov ; 

mm = mean(rs_) ;

% conv w/ hrf
rs_hrf = cell2mat(...
    arrayfun(@(x_) ...
    conv(rs_(:,x_),hrf,'valid'),1:targSz(2), ...
    'UniformOutput',false)) ;

% bandpass the randsignal
rs_hrf_bp = cell2mat(...
    arrayfun(@(x_) ...
    bpfilter(rs_hrf(:,x_),targFq,1/targTR),1:targSz(2), ...
    'UniformOutput',false)')' ;

% put the orignal mean back in
simTS = rs_hrf_bp + mm ;

if nargout > 1
    simMarkov = rs_markov(length(hrf):end) ;
    simGauss = rs_gauss(length(hrf):end) ;
end