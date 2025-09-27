function cross_spectra = denoiser_nd(y_multicoset, scenario_setup)
    
  % y_multicoset: cell
    
    Ts = scenario_setup.Ts;
    lag_channels = scenario_setup.lag_channels;
    num_fft = scenario_setup.num_fft;
    R = scenario_setup.R;
    freq_subNyquist = scenario_setup.freq_subNyquist;


    Mantenna = size(y_multicoset{1},1);    
    Ncomb = nchoosek(1:numel(lag_channels),2);
    num_comb = size(Ncomb,1);

    
    coupled_tensor_G0 = xcross_spectra_nd(y_multicoset{1}, y_multicoset{1}, num_fft);
    num_ccpd = length(coupled_tensor_G0);
    
    %% ccpd 

    options_ft = 1;
    Factor_ft = alg_ccpd_mu_fib(coupled_tensor_G0, R, options_ft);
    
    xauto_spectra = Factor_ft{end};
    
    %% estimating H1H2
    krHHconj = [];
    for cc = 1:num_ccpd 
        
        hfac1 = Factor_ft{2*(cc-1)+1};
        hfac2 = Factor_ft{2*cc};
        hconj12 = kr(hfac2, hfac1);
        krHHconj = [krHHconj; hconj12];
   
    end

    %% flipping negative PSDs, in case there's possibly a scaling
    for rr = 1:R
        [absmax, pos] = max(abs(real(xauto_spectra(:,rr))));
        if absmax/real(xauto_spectra(pos,rr)) < 0
            xauto_spectra(:,rr) = -xauto_spectra(:,rr);
            krHHconj(:,rr) = -krHHconj(:,rr);
        end
    end

    cross_spectra{1} = fftshift(xauto_spectra.',2);

    for cc = 1:num_comb
    
        receive_1 = y_multicoset{Ncomb(cc,1)};
        receive_2 = y_multicoset{Ncomb(cc,2)};
        
        coupled_tensor_Gq = xcross_spectra_nd(receive_1, receive_2,num_fft ...
            );
        
        m3_Gq = [];
    
        for ccc = 1:num_ccpd
                
            m3_Gq = [m3_Gq; tens2mat(coupled_tensor_Gq{ccc},[],3)];
        
        end
    
        %% calibrating phase rotation 
        delta_t = lag_channels(Ncomb(cc,2)) - lag_channels(Ncomb(cc,1)); 
        time_shift = delta_t * Ts; 
        phase_correction = exp(1i * 2 * pi * freq_subNyquist * time_shift);
    
        cross_spectra_cc = fftshift(pinv(krHHconj) * m3_Gq,2);
        cross_spectra{cc+1} = cross_spectra_cc.*phase_correction;
    
    end

end