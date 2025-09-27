function [cross_spectra, Ccol] = denoiser_mr_undet(y_multicoset, scenario_setup)
    
    Ts = scenario_setup.Ts;
    lag_channels = scenario_setup.lag_channels;
    num_fft = scenario_setup.num_fft;
    R = scenario_setup.R;
    freq_subNyquist = scenario_setup.freq_subNyquist;
    Dreceivers = scenario_setup.Dreceivers;
    multicoset_N = scenario_setup.multicoset_N; 

    Mantenna = size(y_multicoset{1}{1}, 1);
    noise_power_est = zeros(Dreceivers,Mantenna);

    tensor_G0 = cell(1,Dreceivers);
    
    for dd = 1:Dreceivers
    
        tensor_G0_d = xcross_spectra(y_multicoset{dd}{1}, y_multicoset{dd}{1}, num_fft); %%cross-correlation tensor direct path
    
        for mm = 1:Mantenna
            tensor_G0_d_fib = squeeze(abs(tensor_G0_d(mm,mm,:)));
            noise_power_est(dd, mm) = min(tensor_G0_d_fib);
            tensor_G0_d(mm,mm,:) = tensor_G0_d(mm,mm,:) - noise_power_est(dd, mm);
        end
    
        tensor_G0{dd} = tensor_G0_d;
    
    end


    %% ccpd 
    options_ft = 1;
    Factor_ft = alg_ccpd_mu_fib(tensor_G0, R, options_ft);
    
    xauto_spectra = Factor_ft{end};

    
    %% flipping negative PSDs, in case there's possibly a scaling
    for rr = 1:R
        [absmax, pos] = max(abs(real(xauto_spectra(:,rr))));
        if absmax/real(xauto_spectra(pos,rr)) < 0
            xauto_spectra(:,rr) = -xauto_spectra(:,rr);
        end
    end



    cross_spectra{1} = fftshift(xauto_spectra.',2);

    num_Q = 2;
    %% cross-correlation tensor Gq && establishing matrix C
    Ccol(1) = 1;

    for dd = 1:Dreceivers
        
        m3_G0d = tens2mat(tensor_G0{dd},[],3);
    
        krHHconj_d = m3_G0d*pinv(xauto_spectra)';
        
        lag_comb = nchoosek(1:numel(lag_channels{dd}), 2);
        
    
        for cc = 1:size(lag_comb,1)
    
            idx1 = lag_comb(1);
            idx2 = lag_comb(2);
            
            Seq1 = y_multicoset{dd}{idx1};
            Seq2 = y_multicoset{dd}{idx2};
    
            tensor_Gqd = xcross_spectra(Seq1, Seq2, num_fft);
    
            %% calibrating phase rotation 
            delta_t = lag_channels{dd}(idx2) - lag_channels{dd}(idx1); 
            time_shift = delta_t * Ts; 
            phase_correction = exp(1i * 2 * pi * freq_subNyquist * time_shift);
    
            Ccol(num_Q) = exp(1i*2*pi*(-delta_t)/multicoset_N);
            
            for mm = 1:Mantenna
                tensor_Gqd_fib = squeeze(tensor_Gqd(mm,mm,:)).';
                noise_power_Gqd = min(abs(tensor_Gqd_fib));
                tensor_Gqd(mm,mm,:) = tensor_Gqd_fib - fftshift(noise_power_Gqd * conj(phase_correction));
            end
            
            m3_Gqd = tens2mat(tensor_Gqd,[],3);
            cross_spectra{num_Q} = fftshift(pinv(krHHconj_d) * m3_Gqd,2) .*phase_correction ;
            num_Q = num_Q + 1;
        end
    
    end
    
end