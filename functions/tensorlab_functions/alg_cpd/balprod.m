function [ind,perm] = balprod(a)
%BALPROD greedy balanced product
%   [ind,perm] = BALPROD(a) splits the entries of the vector a in two sets
%   such that the products of the entries of each set are approximately
%   equal using a greedy approach.

%   Authors:    Martijn Boussé          (Martijn.Bousse@kuleuven.be)        
%
%   Version History:
%       - 2020/09/18    MB      Initial version.
    
    if ~isvector(a)
        error('balprod:a','a should be a vector.')
    end
    N = length(a);

    [val,perm] = sort(a);
    n = 1; 
    check = true;
    while check
        ind = 1:n;        
        % Determine which values to combine in a greedy way.
        if prod(val(ind)) <= prod(val(setdiff(1:N,ind)))
            n = n+1;
        else
            r1 = prod(val(ind))/prod(val(setdiff(1:N,ind)));
            r2 = prod(val(setdiff(1:N,ind(1:end-1))))/prod(val(ind(1:end-1)));
            if r1 >= r2
                ind = ind(1:end-1);
            end
            check = false;
        end           
    end
end