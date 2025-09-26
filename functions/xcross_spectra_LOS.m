function tensor_G_LOS = xcross_spectra_LOS(Seq1, Seq2, R, num_fft, denoise)

    M1 = size(Seq1,1);

    M2 = size(Seq2,1);

    if M1 ~= M2
        error("Seq1 should has same size with Seq2");
    end

    M = M1;

    if 2*M < R+1
        error("Insufficient amount of lags");
    end

    space_lag = -(M-1):(M-1);

    mat_G_LOS = zeros(2*M-1,num_fft);

    for ll = 1:2*M-1
        
        if space_lag(ll) < 0
            pick1 = 1;
            pick2 = 1-space_lag(ll);
        else
            pick1 = space_lag(ll) + 1;
            pick2 = 1;
        end
        
        Seq1_l = Seq1(pick1,:);
        Seq2_l = Seq2(pick2,:); 

        cross_corr = (1/(2*num_fft-1)) * xcorr(Seq1_l.',Seq2_l.',num_fft-1).';

        for ff = 1:num_fft

            mat_G_LOS(ll,ff) = (exp(-1i*2*pi*((ff-1)/num_fft)*[-(num_fft-1):(num_fft-1)])*cross_corr.')/num_fft;

        end
        
        

    end

    %% denoising

    if denoise
        
        noise_power = min(abs(mat_G_LOS(M,:)));

        mat_G_LOS(M,:) = mat_G_LOS(M,:) - noise_power;
    
    end

    %% extending mat to tens
    block_len = 2*M - R;

    for rr = 1:R

        sel_row = rr:rr+block_len - 1;
        tensor_G_LOS(:,:,rr) = mat_G_LOS(sel_row,:);

    end


end