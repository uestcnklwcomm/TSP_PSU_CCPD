function support = bruteforce(y,Q,k)

    [~,n] = size(Q);
    
    support_list = nchoosek(1:n, k);
    num_supports = size(support_list, 1);
    min_residual = inf;
    
    %% Exhausting
    for i = 1:num_supports
        S = support_list(i,:);
        Qs = Q(:, S);
        b_s = Qs \ y;
        residual = norm(y - Qs * b_s);
    
        if residual < min_residual
            min_residual = residual;
            best_support = S;
            best_b_s = b_s;
        end
    end

    support = best_support;

end