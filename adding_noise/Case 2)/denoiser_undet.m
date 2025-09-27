function cross_spectra = denoiser_undet(y_multicoset, scenario_setup)
    
    % y_multicoset: cell

    Ts = scenario_setup.Ts;
    lag_channels = scenario_setup.lag_channels;
    num_fft = scenario_setup.num_fft;
    R = scenario_setup.R;
    freq_subNyquist = scenario_setup.freq_subNyquist;


    Mantenna = size(y_multicoset{1},1);
    Ncomb = nchoosek(1:numel(lag_channels),2);
    num_comb = size(Ncomb,1);


    tensor_G0 = xcross_spectra(y_multicoset{1}, y_multicoset{1}, num_fft); %%cross-correlation tensor of direct path

    for mm = 1:Mantenna

        tensor_G0_fib = squeeze(abs(tensor_G0(mm,mm,:)));
        noise_power_est(mm) = min(tensor_G0_fib);
        tensor_G0(mm,mm,:) = tensor_G0(mm,mm,:) - noise_power_est(mm);

    end
    
    %% ccpd 
    options_ft = 1;
    Factor_ft = alg_ccpd_m3_fib(tensor_G0, R, options_ft);

    xauto_spectra = Factor_ft{end};

    %% estimating H (not necessary)
    
    H_est = zeros(Mantenna,R);
    m3_G0 = tens2mat(tensor_G0,[],3);

    %% flipping negative PSDs, in case there's possibly a scaling
    for rr = 1:R
        [absmax, pos] = max(abs(real(xauto_spectra(:,rr))));
        if absmax/real(xauto_spectra(pos,rr)) < 0
            xauto_spectra(:,rr) = -xauto_spectra(:,rr);
        end
    end

    krHHconj = m3_G0/(xauto_spectra.');
    
    for rr = 1:R
        rank1r = krHHconj(:,rr);
        [H_est(:,rr),res] = rank_1_fac_conjugate(reshape(rank1r,[Mantenna,Mantenna]));
        % [Hestn(:,rr),res] = rank_1_phase_correction(reshape(rank1r,[Mantenna,Mantenna]));
    end
   

    cross_spectra{1} = fftshift(xauto_spectra.',2);

    %% cross-correlation tensor Gq

    for cc = 1:num_comb
    
        receive_1 = y_multicoset{Ncomb(cc,1)};
        receive_2 = y_multicoset{Ncomb(cc,2)};
        
        tensor_Gq = xcross_spectra(receive_1, receive_2,num_fft ...
            );
        
        %% calibrating phase rotation 
        delta_t = lag_channels(Ncomb(cc,2)) - lag_channels(Ncomb(cc,1)); 
        time_shift = delta_t * Ts; 
        phase_correction = exp(1i * 2 * pi * freq_subNyquist * time_shift);
    
        %% denoising
        for mm = 1:Mantenna                
            tensor_Gq_fib = squeeze(tensor_Gq(mm,mm,:)).';
            noise_power_Gq = min(abs(tensor_Gq_fib));
            tensor_Gq(mm,mm,:) = tensor_Gq_fib - fftshift(noise_power_Gq * conj(phase_correction));
        end
    
        m3_Gq = tens2mat(tensor_Gq,[],3);
    
        cross_spectra_cc = fftshift(pinv(krHHconj) * m3_Gq,2);
        cross_spectra{cc+1} = cross_spectra_cc.*phase_correction;

    end



end