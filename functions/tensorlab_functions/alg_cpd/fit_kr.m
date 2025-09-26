function [U,V,output] = fit_kr(A, size_tens, varargin)
%FIT_KR Khatri-Rao structure
%   [U,V,output] = FIT_KR(A, size_tens) computes the factor matrices U{1},
%   U{2}, ...,U{N} of dimensions size_tens(n)-by-size(A,2) belonging to a 
%   Khatri-Rao decomposition of the matrix A by computing a rank-1 
%   approximation of each tensorized column of A of dimensions size_tens. 
%   V contains kr(U).
%
%   FIT_KR(A, size_tens, options) allow the user to set the following
%   options:
%
%       - Flip = [{false},true] Flip the order of the dimensions.
%       - Refinement =          By default the fit_kr uses mlsvd to compute
%         [{false},true]        the rank-1 approximation. If true, cpd is
%                               used, which entails additional refinement.
%       - AlgorithmOptions      Struct to pass options to underlying
%                               algorithm.
%
%   The structure output returns additional information:
%
%       output.Refinement   - True if refinement was used.
%       output.Algorithm    - Information on the algorithm that was used 
%                             for the rank-1 approximation.
%       output.Flip         - True if order was flipped.
%
%   Example:
%       size_tens = [10 9 8]; R = 3;
%       U = cpd_rnd(size_tens,R);
%       A = kr(U);
%       [Uest,Aest] = fit_kr(A,size_tens);
%       cpderr(U,Uest)
%       frob(A-Aest)/frob(A)
%
%   See also fit_vander, fit_vander_mix, fit_rational, fit_poly, fit_orth,
%   fit_eqfac, fit_toeplitz

%   Authors: Martijn Bouss          (Martijn.Bousse@esat.kuleuven.be)
%            Wouter Deleersnyder    (Wouter.Deleersnyder@kuleuven.be)
%            Lieven De Lathauwer    (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   Versions: 
%       - 31/10/2019    MB+WD   Original version.
%       - 11/12/2019    MB      Improved implementation.
%       - 16/03/2020    MB      Changed output structure.
%       - 25/02/2021    MB      Added optional refinement step and flip.
    
    % Check inputs.
    if size(A,1) ~= prod(size_tens)
        error('fit_kr:A','size(A,1) should be equal to prod(size_tens).')
    end
    
    % Set some parameters.
    N = length(size_tens);
    R = size(A,2);

    % Check options structure.
    p = inputParser;
    p.addOptional('Flip', false);
    p.addOptional('Refinement', false);
    p.addOptional('AlgorithmOptions', struct);    
    p.parse(varargin{:});   
    
    fn = fieldnames(p.Results);
    data = struct2cell(p.Results);
    options = cell2struct(data, fn);
    
    % Flip the order of size_tens if requested.
    if options.Flip
        size_tens = fliplr(size_tens);
    end
    
    % Compute rank-1 approximation ...
    U = arrayfun(@(i) zeros(i,R), size_tens, 'UniformOutput', false);
    if options.Refinement % ... using CPD                
        for r = 1:R            
            % Compute approximation.
            [Ur,out] = cpd(reshape(A(:,r),fliplr(size_tens)), 1, ...
                        options.AlgorithmOptions);
                
            % Set factor matrices.            
            for n = 1:N
                U{n}(:,r) = Ur{N-n+1};
            end                    
        end
        
        % Set output structure.
        output.Refinement = true;
        output.Algorithm.Name = 'cpd';
        output.Algorithm.Output = out;
        
    else % ... using MLSVD                
        for r = 1:R   
            % Compute approximation.
            [Ur,Sr,out] = mlsvd(reshape(A(:,r),fliplr(size_tens)), ...
                                    ones(1,N), options.AlgorithmOptions);

            % Set factor matrices.
            U{1}(:,r) = Ur{N}*Sr;
            for n = 2:N
                U{n}(:,r) = Ur{N-n+1};
            end
        end
        
        % Set output structure.
        output.Refinement = false;
        output.Algorithm.Name = 'mlsvd';
        output.Algorithm.Output = out;
        
    end
    
    % Optionally, return kr(U).
    if nargout > 1
        V = kr(U);
    end
        
    % Flip back if needed.
    if options.Flip
        U = U(end:-1:1);
        output.Flip = true;
    else
        output.Flip = false;
    end

end