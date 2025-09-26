function xcorr_spectra = xcross_spectra(Seq1,Seq2,num_fft)
%% returning the cross-correlation spectra of two received sequences

    M1 = size(Seq1,1);
    M2 = size(Seq2,1);
    
    if size(Seq1,2) ~= size(Seq2,2)
        error("Unequivalent sequence length");
    end
    
    xcorr_spectra = zeros(M1,M2,num_fft);
    for m = 1:M1
        Seq1_m = Seq1(m,:);
        for mm = 1:M2
            Seq2_mm = Seq2(mm,:);
            cross_corr = (1/(2*num_fft-1)) * xcorr(Seq1_m.',Seq2_mm.',num_fft-1).';
            for ff = 1:num_fft
                
                xcorr_spectra(m,mm,ff) = (exp(-1i*2*pi*((ff-1)/num_fft)*[-(num_fft-1):(num_fft-1)])*cross_corr.')/num_fft;
            end
        end
    
    end


end