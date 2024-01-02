function [itc,cmphase] = itc_map3(lfp,cycle)
lfp = lfp';
% wavelet parameters
num_frex = 60;
min_freq =  1;
max_freq = 60;



% set range for variable number of wavelet cycles
range_cycles = [cycle cycle];

% other wavelet parameters
% frex = logspace(log10(min_freq),log10(max_freq),num_frex);
frex = linspace(min_freq,max_freq,num_frex);
% wavecycles = logspace(log10(range_cycles(1)),log10(range_cycles(end)),num_frex);
wavecycles = linspace(range_cycles(1),range_cycles(end),num_frex);
time = -0.25:1/1000:0.25;
half_wave_size = (length(time)-1)/2;

% FFT parameters
nWave = length(time);
nData = size(lfp,1)*size(lfp,2);
nConv = nWave+nData-1;
A1 = reshape(lfp,1,nData);

% FFT of data (doesn't change on frequency iteration)
dataX = fft( A1 ,nConv);

% initialize output time-frequency data
% tf = zeros(num_frex,size(LFP_fix_in,1));

% loop over frequencies
for fi=1:num_frex
    
    % create wavelet and get its FFT
    s = wavecycles(fi)/(2*pi*frex(fi));
    wavelet  = exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*s^2));
%     [psi,x] = cmorwavf(-0.2,0.2,401,cycle/1000,fi);
% wavelet = psi;
    waveletX = fft(wavelet,nConv);
    
    % run convolution
    as = ifft(waveletX.*dataX,nConv);
    as = as(half_wave_size+1:end-half_wave_size);
    as = reshape(as,size(lfp,1),size(lfp,2));
    cmphase{1,fi} = exp(1i*angle(as'));
    %  cmphase{1,fi} = as';
    % compute ITPC
    
    tf(fi,:) = abs(mean(cmphase{1,fi}));
    %  tf(fi,:) = abs(mean(exp(1i*angle(cmphase{1,fi}))));
end
itc = tf;
end

