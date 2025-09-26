function coupled_tensor = xcross_spectra_nd(Seq1, Seq2, num_fft)
    [M1, data_len_1] = size(Seq1);
    [M2, data_len_2] = size(Seq2);

    if data_len_1 ~= data_len_2
        
        error('Unequivalent data length');

    end

    row_comb = nchoosek(1:M1,2);
    num_ccpd = 0;


    for cc = 1:size(row_comb,1)

        row_comb_idx = row_comb(cc,:);
        col_idx = 1:M2;
        col_idx(row_comb_idx) = [];

        if numel(col_idx) <= 1
        
            error('Unable to create sufficient non-diagonal combiniations');
            
        end
        col_comb_cc = nchoosek(col_idx, 2);

        for ccc = 1:size(col_comb_cc,1)

            col_comb_idx = col_comb_cc(ccc,:);
            seq_row = Seq1(row_comb_idx, :);
            seq_col = Seq2(col_comb_idx, :);

            num_ccpd = num_ccpd + 1;

            coupled_tensor{num_ccpd} = xcross_spectra(seq_row, seq_col, num_fft);
        end

    end



end