function [U,cidx,output] = ccpd_sd(datasets,R,varargin)
%CCPD_SD Coupled CPD using simultaneous diagonalization.
%   [U,cidx,output] = CCPD_SD(datasets,R) computes the factor matrices 
%   U{1}, U{2}, ..., U{N} belonging to a coupled canonical polyadic 
%   decomposition of the Nth-order tensors given in datasets with coupling 
%   in the last mode and rank R. The problem is reformulated as 
%   simultaneous matrix diagonalization, which is solved by computing the 
%   CPD of an auxiliary tensor. cidx stores the coupling indices.
%
%   CCPD_SD(datasets,R,options) allows the user to set the following 
%   options:
%
%       Set options.Compression = [{True}|False] to false to skip the
%       dimensionality reduction step that is executed before constructing 
%       the SD problem. Only do so when the data is already compressed,
%       i.e. the size of the shared Nth mode in each tensor in datasets has
%       been compressed to size R.
% 
%       options.Initialization = [{@cpd_gevd}|@cpd3_extQZ|@cpd_rnd]
%       determines what algorithm is used to initialize the CPD of the
%       auxiliary tensor.
%       
%       options.InitializationOptions can be used to specify options for
%       the initialization algorithm mentioned above.
%
%       options.Algorithm = [{@cpd_nls}|@cpd_als|@cpd_minf|...|[]] 
%       determines what algorithm is used to improve the initial CPD of the 
%       auxiliary tensor.
%
%       options.AlgorithmOptions can be used to specify options for the 
%       algorithm mentioned above.
%
%       If options.IndirectP = [{True}|False] is true, the algorithm uses 
%       an indirect approach for working with the rank-1 detection device, 
%       which is more efficient for large tensors.
%   
%       Set options.OnlyReturnF = [{False}|True] to true if only the matrix 
%       F is needed, instead of the complete CCPD.
%
%   The structure output returns additional information on the computation
%   of the factor matrices:
%
%      output.fval       - The value of the objective function in every 
%                          iteration.
%      output.info       - The circumstances under which the procedure
%                          terminated:
%                             1: Objective function tolerance reached.
%                             2: Step size tolerance reached.
%                             3: Maximum number of iterations reached.
%                             4: Absolute objective function tolerance
%                                reached.
%      output.iterations - The number of iterations.
%
%   Example:
%       size_tens = {[5 4 3 10],[10 3 3 10]}; R = 8;
%       [U,cidx] = ccpd_rnd(size_tens,R);
%       datasets = ccpdgen(U,cidx);
%       Uest = ccpd_sd(datasets,R);
%       cpderr(U,Uest)
%
%   See also ccpd_nls, ccpd_minf, ccpd_als.

%   Authors:    Mikael Sorensen         
%               Martijn Boussé          (Martijn.Bousse@kuleuven.be)
%               Stijn Hendrikx          (Stijn.Hendrikx@kuleuven.be)
%               Nico Vervliet           (Nico.Vervliet@kuleuven.be)
%               Lieven De Lathauwer     (Lieven.DeLathauwer@kuleuven.be) 
%
%   Version History:
%       - 2021/01/18    SH      Extra options.
%       - 2020/08/19    MB      Higher-order extension.
%       - 2020/08/04    MB      Improved implementation and documentation.
%       - 2020/06/27    MS      Initial version.
%
%   References:
%   [1] M. Sorensen, L. De Lathauwer, "Coupled canonical polyadic 
%       decompositions and (coupled) decompositions in multilinear 
%       rank-(L_{r,n}, L_{r,n}, 1) terms --- Part II: Algorithms'', SIAM 
%       Journal on Matrix Analysis and Applications, vol. 36, no. 3, July 
%       2015, pp. 1015-1045.

    % Check inputs.
    if ~iscell(datasets) || ~isvector(datasets)
        error('ccpd_sd:datasets','Datasets should be a 1D cell array.')
    end
    % If datasets is a column cell array, reshape it into a row cell array.
    if size(datasets, 1) ~= 1
        datasets = datasets.';
    end
    P = length(datasets);
    if P == 1
        error('ccpd_sd:datasets','Use regular cpd algorithms.')
    end
    N = getorder(datasets{1});
    if any (cellfun('ndims', datasets) ~= N)
        error('ccpd_sd:datasets','getorder(datasets{p}) should be the same for all p.')
    end
    if any(cellfun('size', datasets, N) ~= size(datasets{1},N))
        error('ccpd_sd:datasets','size(datasets{p},N) should be the same for all p.')
    end
    size_tens = cellfun(@(t) size(t), datasets, 'UniformOutput', false);
    if R > size(datasets{1},N)
        error('ccpd_sd:R','R should be smaller than or equal to size(datasets{p},N).')
    end
    
    % Check the options structure
    xsfunc = @(f)isa(f,'function_handle')&&exist(func2str(f),'file');
    p = inputParser;
    p.addOptional('IndirectP', true);
    p.addOptional('Initialization', @cpd_gevd, xsfunc);
    p.addOptional('InitializationOptions', struct(), @isstruct);
    p.addOptional('Algorithm', @cpd_nls);
    p.addOptional('AlgorithmOptions', struct(), @isstruct);    
    p.addOptional('OnlyReturnF', false);
    p.addOptional('Compression', true);
    p.KeepUnmatched = true;
    p.parse(varargin{:});

    fn = [fieldnames(p.Results); fieldnames(p.Unmatched)];
    data = [struct2cell(p.Results); struct2cell(p.Unmatched)];
    options = cell2struct(data, fn);        
    
    if ~isempty(options.Algorithm)
        if ~xsfunc(options.Algorithm)
            error('cpd3_sd_fs:Algorithm', ['options.Algorithm should either '...
                'be a function or empty.']);
        end
    end
    
    % Combine modes using greedy approach, if getorder(T) > 3.     
    if N > 3
        comb = cell(1,P);
        perm = cell(1,P);
        for p = 1:P
            % Determine mode combination
            [ind,perm{p}] = balprod(size_tens{p}(1:N-1));

            % Permute and combine modes
            datasets{p} = permute(datasets{p},[perm{p} N]);
            datasets{p} = reshape(datasets{p},prod(size_tens{p}(perm{p}(ind))),[],size_tens{p}(end));
            comb{p} = ind;    
        end
    end
    
    % Dimensionality reduction    
    Y = cellfun(@(y) tens2mat(y,[],3), datasets, 'UniformOutput', false);
    Y = cell2mat(Y.');
    if options.Compression
        [Uy,Sy,Vy] = svd(Y,0);    
        Ey = Uy(:,1:R)*Sy(1:R,1:R);  
    else
        if ~all(cellfun(@(d) size(d,N), datasets) == R)
            error('ccpd_sd:options',['Only disable compression if ' ...
                'the data has already been compressed: size(datasets{p}' ...
                ',N) = R) for all p.']);
        end
        Ey = Y;
        Vy = eye(size(Y,2));
    end
    
    if R > 1
        
        % Check if method can be applied (after combining modes).    
        sz = cellfun(@(t) [size(t,1) size(t,2)], datasets, 'UniformOutput', false);
        maxsz = sum(cellfun(@(s) nchoosek(s(1),2)*nchoosek(s(2),2), sz));    
        if nchoosek(R,2) > maxsz
            error('ccpd_sd:R','R is too large.')
        end

        % Construct Ptot
        size_tensc = cellfun(@(t) size(t), datasets, 'UniformOutput', false);
        I = cellfun('size', datasets, 1); 
        J = cellfun('size', datasets, 2);
        if options.IndirectP % indirect approach (P'*P(w kron w)=0, useful for large P)
            Ptot = 0;
            for p = 1:P
                ind1 = I(1:p-1)*J(1:p-1).'+1;
                ind2 = I(1:p)*J(1:p).';
                Eyp = Ey(ind1:ind2,:);
                G = rank1_mapping(Eyp,size_tensc{p}(1:2),true);
                Ptot = Ptot + G / p;
            end
        else % direct approach (P*(w kron w)=0)
            Ptot = zeros(sum((I.*J).^2),R*(R+1)/2);
            offset = cumsum([1 (I.*J).^2]);
            for p = 1:P
                ind1 = I(1:p-1)*J(1:p-1).'+1;
                ind2 = I(1:p)*J(1:p).';            
                Eyp = Ey(ind1:ind2,:);             
                indp = offset(p):offset(p+1)-1;
                Ptot(indp,:) = rank1_mapping(Eyp,size_tensc{p}(1:2));
            end
        end
    
        % Construct auxiliary tensor.
        M = symmetric_null(Ptot,R);    

        % Compute CPD of auxiliary tensor.
        % Initial solution.
        if abs(nargin(options.Initialization)) == 2
            % Initialization does not need rank as an argument, like cpd3_extQZ.
            [Um, output.Initialization] = options.Initialization(M, ...
                options.InitializationOptions);
        elseif strcmpi(func2str(options.Initialization), 'cpd_rnd')
            % Initialization is random.
            Um = options.Initialization([R R R], R, options.InitializationOptions);
        else
            % Initialization needs rank as an argument, like cpd_gevd.
            [Um, output.Initialization] = options.Initialization(M, R, ...
                options.InitializationOptions);
        end
        % Optimization algorithm.
        if ~isempty(options.Algorithm)
            [Um,output.Algorithm] = options.Algorithm(M, Um, options.AlgorithmOptions);
        end
        F = Um{1};
    else
        output = struct;
        F = 1;        
    end
    
    % Return only F or compute all factor matrices.
    if options.OnlyReturnF
        U = F;
        cidx = [];
    else    
        % Determine coupled factor matrix
        V = cell(1,2*P+1);
        V{end} = conj(Vy(:,1:R))/F.';

        % Determine other factor matrices (via best rank-1 approximation)    
        tmp = [reshape(1:length(V)-1,2,P); repmat(length(V),1,P)];
        cidx = mat2cell(tmp(:).',1,repmat(3,1,P));
        for p = 1:P        
            X = tens2mat(datasets{p},[],3)/V{end}.';
            u = fit_kr(X, size(datasets{p},[2 1]));       
            V(cidx{p}([1 2])) = u([2 1]);
        end    
        
        % Compute rank-1 approximations of combined modes, if needed.    
        if N > 3    
            U = cell(1,(N-1)*P+1);
            for p = 1:P
                sz = sort(size_tens{p}(1:N-1));
                Utmp = cell(1,N-1);
                idxp = (p-1)*2+1; 
                if length(comb{p}) == 1
                    % Vest{idxp} does not contain combined modes.
                    Utmp{1} = V{idxp};

                    % Extract modes from Vest{idxp+1}
                    ind = fliplr(setdiff(1:N-1,comb{p}));
                    Utmp(ind) =  fit_kr(V{idxp+1},sz(ind));
                else                    
                    % Extract modes from Vest{idxp}
                    ind = fliplr(comb{p});
                    Utmp(ind) = fit_kr(V{idxp},sz(ind));
                    if length(setdiff(1:N-1,comb{1})) == 1
                        % Vest{idxp+1} does not contain combines modes.
                        Utmp{end} = V{idxp+1};
                    else
                        % Extract modes from Vest{idxp+1}
                        indc = fliplr(setdiff(1:N-1,comb{p}));
                        Utmp(indc) = fit_kr(V{idxp+1},sz(indc));
                    end
                end

                % Permute back.
                idx = (p-1)*(N-1)+1:p*(N-1);
                U(idx(perm{p})) = Utmp;                        
            end

            % Assign coupled mode.
            U{end} = V{end};
        else
            U = V;
        end        

        % Determine cidx for output.
        tmp = [reshape(1:length(U)-1,N-1,P); repmat(length(U),1,P)];
        cidx = mat2cell(tmp(:).',1,repmat(N,1,P));
    end     
end