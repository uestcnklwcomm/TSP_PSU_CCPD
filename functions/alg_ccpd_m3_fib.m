function Fac = alg_ccpd_m3_fib(tensor, rank, options_opt)

    [M1,M2,M3] = size(tensor);

    if M3 < rank
        error('Mode-3 factor is rank deficient');
    end

    if nchoosek(M1,2) * nchoosek(M2,2) < nchoosek(rank,2)

        error('Amont of coupled fibers is insufficient');

    end

    %% constructing coupled tensors

    Num_ccpd = 1;
    Comb_row = nchoosek([1:M1],2);
    Comb_col = nchoosek([1:M2],2);

    for c_row = 1:nchoosek(M1,2)
        row_idx = Comb_row(c_row,:);
        for c_col = 1:nchoosek(M2,2)
    
            col_idx = Comb_col(c_col,:);
            ccpd{Num_ccpd} = tensor(row_idx,col_idx,:);
    
            Num_ccpd = Num_ccpd + 1;
        end
    
    end
    
    Num_ccpd = Num_ccpd - 1;


    model_ccpd = struct;

    for nn = 1:Num_ccpd
        model_ccpd.factorizations{nn}.data = ccpd{nn};
        model_ccpd.factorizations{nn}.cpd = [1 + 2*(nn-1), 2*nn ,2*Num_ccpd+ 1];
    end

    %% initialization via sd
    Factor_init = ccpd_sd(ccpd, rank);
    Fac = Factor_init;
    % 
    if options_opt 
        model_ccpd.variables = Factor_init;
        model_ccpd.factors = 1:2*Num_ccpd + 1;
        
    
        options_nls.TolFun = 1e-8;
        options_nls.MaxIter = 1000;
        options_nls.Weights = ones(1,Num_ccpd);
        
        [Factor_ft,~] = ccpd_nls(model_ccpd,options_nls);
        
        Fac = Factor_ft;
    end


end