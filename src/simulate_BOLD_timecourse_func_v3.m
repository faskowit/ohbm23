function [ts_sim_cov,ts_sim_nocov]= simulate_BOLD_timecourse_func_v3(...
    timelength,TRin,TRout,cov_target,P_target)
% Inputs:
% timelength - desired time length of output simulated data
% TRin - TR of real time series from which cov_target and P_target are generated
% TRout - desired TR of output time series
% cov_target - target covariance matrix 
% P_target - target power spectrum
% 
% Outputs:
% ts_sim_cov - covariance and spectrally matched simulated time series
% ts_sim_nocov - before the rotation step
%
% TOL, 06/2016
% joshedit

%TRout must be greater than or equal to TRin
if TRin>TRout
    error('TRout must be greater than or equal to TRin')
end

%%Generate random normal deviates
wn = randn(timelength,size(cov_target,1));

%%Normalize to unit variance
wn_std = std(wn);
wn_mean = mean(wn);
wn_unit = (wn-repmat(wn_mean,timelength,1))./repmat(wn_std,timelength,1);

%Identify eigenvalues and eigenvectors of target covariance
[V,D] = eig(cov_target);

%Check that matrix is positive definite
if sum(diag(D)<0)>0
    error('Covariance matrix is not positive definite')
end

%%Power spectral density resampling
if mod(timelength,2)
    midpoint = ceil(timelength./2);
else
    midpoint = timelength./2+1;
end

% Find frequency indices to interpolate
%if mod(length(P_target),2)
    P_target_length = ceil(length(P_target)/2)+1;
%else
%    P_target_length = length(P_target)/2+1;
%end

% count = 1;
% for k = 0:timelength/2
%     Vq(count) = (k*TRin*length(P_target))/(TRout*(timelength))+1;
%     count = count + 1;
% end

% create resampling inds
ii = 0:timelength/2 ;
ll = length(P_target) ;
Vq = nan(1,length(ii)) ;
for k = ii
    Vq(k+1) = (k*TRin*ll)/(TRout*(timelength))+1;
end

% Interpolate 
resamp_P_target = interp1(P_target(1:P_target_length),Vq);

% Flip over nyquist
resamp_P_target = [resamp_P_target fliplr(resamp_P_target(2:(midpoint-1)))]';
resamp_P_target(1) = 0;

%%Enforce power spectra
fft_timecourse_sim = fft(wn_unit);
resamp_P_target_rep = repmat(resamp_P_target,[1 size(fft_timecourse_sim,2)]);
F_sim = fft_timecourse_sim.*sqrt(resamp_P_target_rep);

ts_sim_nocov = ifft(F_sim);
ts_sim_nocov = (ts_sim_nocov-repmat(mean(ts_sim_nocov),size(ts_sim_nocov,1),1))./repmat(std(ts_sim_nocov),size(ts_sim_nocov,1),1);

%%Enforce covariance
%Multiply eigenvectors from real data into random normal deviates
X = ts_sim_nocov;
V_use = V*sqrt(D);
Y = V_use*X';
ts_sim_cov = Y';
