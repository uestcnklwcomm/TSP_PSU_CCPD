%% ccpd for multi-receiver scheme
function Fac = alg_ccpd_mu_fib(coupled_tensor, rank, options_opt)

    num_coupled = length(coupled_tensor);

    num_compound_row = 0;
    for nn = 1:num_coupled
        
        [M1,M2,~] = size(coupled_tensor{nn});

        num_compound_row = num_compound_row + nchoosek(M1,2) * nchoosek(M2,2);
        

    end

    if num_compound_row < nchoosek(rank,2)

        error('Amount of coupled fibers is insufficient');

    end


    Factor_init = ccpd_sd(coupled_tensor, rank);
    Fac = Factor_init;
    %% fine-tunning

    model_ccpd = struct;
    for cc = 1:num_coupled
        model_ccpd.factorizations{cc}.data = coupled_tensor{cc};
        model_ccpd.factorizations{cc}.cpd = [1 + 2*(cc-1), 2*cc, 2*num_coupled + 1];
    end

    if options_opt

        model_ccpd.variables = Factor_init;
        model_ccpd.factors = 1:2*num_coupled + 1;
        options_nls.TolFun = 1e-8;
        options_nls.MaxIter = 1000;
        options_nls.Weights = ones(1,num_coupled);
        [Factor_ft,~] = ccpd_nls(model_ccpd,options_nls);
        Fac = Factor_ft;
    
    end



end
