% calculating the ground-truth power spectra of an input signal 
function ps = powerspectra_gt(seq,num_fft)
    
    ps = zeros(1,num_fft);
    auto_corr = (1/(2*num_fft-1)) * xcorr(seq.',seq.',num_fft-1).';
    
    for ff = 1:num_fft
        
        ps(ff) = (exp(-1i*2*pi*((ff-1)/num_fft)*[-(num_fft-1):(num_fft-1)])*auto_corr.')/num_fft;
    end


end