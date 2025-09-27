function [cross_spectra, Ccol] = denoiser_mr_nd(y_multicoset, scenario_setup)
    
    % y_multicoset: double cell

    Ts = scenario_setup.Ts;
    lag_channels = scenario_setup.lag_channels;
    num_fft = scenario_setup.num_fft;
    R = scenario_setup.R;
    freq_subNyquist = scenario_setup.freq_subNyquist;
    Dreceivers = scenario_setup.Dreceivers;
    multicoset_N = scenario_setup.multicoset_N; 

    Mantenna = size(y_multicoset{1}{1}, 1);
    
    coupled_tensor_G_mr = [];

    for dd = 1:Dreceivers
        
        coupled_tensor_G0 = xcross_spectra_nd(y_multicoset{dd}{1}, y_multicoset{dd}{1}, num_fft);
    
        coupled_tensor_G_mr = [coupled_tensor_G_mr, coupled_tensor_G0];
    
    end

    %% ccpd 
    options_ft = 1;
    Factor_ft = alg_ccpd_mu_fib(coupled_tensor_G_mr, R, options_ft);
    xauto_spectra = Factor_ft{end};
    
    num_ccpd = (length(Factor_ft) - 1)/(2*Dreceivers);

    %% estimating kr(H, H*)
    for dd = 1:Dreceivers
        
        idx_start = 1 + 2*num_ccpd*(dd-1);
        krHHconj_d = [];
        for cc = 1:num_ccpd
    
            hfac1 = Factor_ft{2*(cc-1) + idx_start};
            hfac2 = Factor_ft{2*cc + idx_start - 1};
            hconj12 = kr(hfac2, hfac1);
            krHHconj_d = [krHHconj_d; hconj12];
    
        end
    
        krHHconj{dd} = krHHconj_d;
    
    end


    %% flipping negative PSDs, in case there's possibly a scaling
    for rr = 1:R
        [absmax, pos] = max(abs(real(xauto_spectra(:,rr))));
        if absmax/real(xauto_spectra(pos,rr)) < 0
            xauto_spectra(:,rr) = -xauto_spectra(:,rr);
            for dd = 1:Dreceivers
                krHHconj{dd}(:,rr) = -krHHconj{dd}(:,rr);
            end
    
        end
    end


    cross_spectra{1} = fftshift(xauto_spectra.',2);

    num_Q = 2;
    %% cross-correlation tensor Gq && establishing matrix C
    Ccol(1) = 1;
    
    for dd = 1:Dreceivers
        
        Ncomb = nchoosek(1:numel(lag_channels{dd}),2);
        num_comb = size(Ncomb,1);
        krHHconj_d = krHHconj{dd};
        
        for cc = 1:num_comb
    
            idx1 = Ncomb(cc,1);
            idx2 = Ncomb(cc, 2);
            
            
            
            Seq1 = y_multicoset{dd}{idx1};
            Seq2 = y_multicoset{dd}{idx2};
    
            coupled_tensor_Gqd = xcross_spectra_nd(Seq1, Seq2, num_fft);
            
            m3_Gqd = [];
            for ccc = 1:num_ccpd
                
                m3_Gqd = [m3_Gqd; tens2mat(coupled_tensor_Gqd{ccc}, [] ,3)];
    
            end
    
            %% calibrating phase rotation 
            delta_t = lag_channels{dd}(idx2) - lag_channels{dd}(idx1); 
            time_shift = delta_t * Ts; 
            phase_correction = exp(1i * 2 * pi * freq_subNyquist * time_shift);
    
            Ccol(num_Q) = exp(1i*2*pi*(-delta_t)/multicoset_N);
        
            cross_spectra{num_Q} = fftshift(pinv(krHHconj_d) * m3_Gqd,2) .*phase_correction ;
            num_Q = num_Q + 1;
        end

    end

end