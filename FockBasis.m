% Filippos Tzimkas-Dakis   Virginia Tech   January 2024
%
% Any feedback and suggestions are much appreciated! 
%    
%     ----->  dakisfilippos@vt.edu  <-------
%
% Version V 1.1
% 
% The present script defines a class whose objects behave as Fock/Number states.
% Every object of this class has all the basic features of a quantum Fock state .
% For instance one can define a state, normalize it, compute overlaps between states,  
% add quantum states, apply displacement operators, compute the Wigner function, calculate 
% the photon distribution, average photo number.
%
% Please have a look at the examples accompanying this script 
%
% I might add more features in the future.


classdef FockBasis
    properties (Hidden)
        Lnk     % Matrix containing the coefficients of the associated Laguerre Polynomials  L_{n}^{k} ,up to  m-th order
                % For more info check the constructor  " function obj = Nbasis(n,N_space) "
                % We compute and store the coefficients to make our code faster !!  .
        Parity  % Parity matrix diag(1,-1,1,-1 ....)
    end
    properties
        Coeff {mustBeNumeric}        % column vector that contains the coefficient of every Fock/Number state.
        n     {mustBeNonnegative}    % column vector depicting the populated number states, ie n = [1 4 5] 
                                     % means that the |1> |4> and |5> have nonzero coefficients.  .
    end 
    methods
        
        function obj = FockBasis(n,N_space) 
            % constructor of the Number Basis class
            % inputs:   n : [vector] column vector containing the coefficients of every state number state, ie n = [0 1 1i] = |1> + 1i|2> .
            %           N_space: [number] truncates the infinite Hilbert space to N_space levels. 
            if nargin==1
                N_space = 30;       % if N_space is not defined, is set to  N_space = 30  which is a good approximation           
            end                    
            if length(n)>=N_space   % N_space must always be larger than length(n)
                N_space = length(n) + 15;
            end
            obj.Coeff = [n(:);zeros(N_space -length(n),1)]; % assign the input values to class properties
            obj.n     = find(obj.Coeff)-1;                  % assign the populated number states
            
            % ------------------------------------------------------------
            % Here we calculate coefficients of the associated Laguerre Polynomials,  L_{n}^{k} ,up to  m-th order.
            % for more info check   https://mathworld.wolfram.com/AssociatedLaguerrePolynomial.html
            % We will use then in the Displacement operator. To that end, we calculate them once and then use their values through the matrix    .
            % We do so, because they contain factorials, which is a time consuming calculation, and we aim to use factorials a less as possible .
            obj.Lnk = zeros(N_space,N_space,N_space);       
            for k = 0:N_space-1
                for nn= 0:N_space-1
                    for m=0:nn
                        obj.Lnk(nn+1,m+1,k+1) = 1/factorial(m) *nchoosek(k+nn,nn-m)*(-1)^m;   
                        % L(nn+1, m+1, k+1)   = m-th coefficient of L_{n}^{k}    .
                        % SOS --- EVERY CELL CONTAINS ONLY THE m-th ORDER COEFFICIENT OF THE SERIES, the value of  (x) is NOT included .
                    end
                end
            end
            % -------------------------------------------------
        obj.Parity = sparse(diag(cos((0:(N_space-1))*pi)));   % Parity matrix
        % obj.Parity = diag(cos((0:(N_space-1))*pi));
        end

        function obj = normalize(obj)                % This function normalizes the state (the obj.Coeff vector).
            obj.Coeff = obj.Coeff/norm(obj.Coeff);
        end

        function q = braket(obj1,obj2)        % This function computes the dot product, namely   q = <ψ_1|ψ_2>  . 
            if nargin == 1                    % if nargin = 1 q = <ψ_1|ψ_1>  
                obj2  = obj1;
            else
                l = max(length(obj1.Coeff),length(obj2.Coeff));
                obj1.Coeff = [obj1.Coeff; zeros(l-length(obj1.Coeff),1)];    % make sure both vectors have the same Hilbert space
                obj2.Coeff = [obj2.Coeff; zeros(l-length(obj2.Coeff),1)];    % we keep the largest Hilbert space.

             end
            q = obj1.Coeff'*obj2.Coeff;
        end

        function obj = plus(obj1,obj2)       % This function adds up two objects,  |ψ_3> = |ψ_1> + |ψ_2>   (it DOES NOT normalize the final vector/object) .
            l1 = length(obj1.Coeff);         % If a Number state is included in both |ψ_1>,|ψ_2> we merge its coefficients 
            l2 = length(obj2.Coeff);
            if l1 == l2
                s = obj1.Coeff + obj2.Coeff;
            elseif l1 > l2
                s = obj1.Coeff + [obj2.Coeff(:); zeros(l1-l2,1)];
            else
                s = obj2.Coeff + [obj1.Coeff(:); zeros(l2-l1,1)];
            end
            obj = FockBasis(s,max(l1,l2));    % the final object/vector has dimensions  max(l1,l2)
        end

        function [obj, D] = D_(obj,a)
            % Displacement operator function 
            % We calculate D(a) in the number basis. We truncate the operator in the N_space^2 Hilbert space .
            % For large Hilbert space D(a)|ψ> approaches the exact result but the code becomes slower - KEEP THAT IN MIND.
            % We calculate    <m|D(a)|n> = D_{m,n}(a) using Glauber's formula, or this link  https://physics.stackexchange.com/questions/553225/representation-of-the-displacement-operator-in-number-basis.
            D = zeros(length(obj.Coeff));      % initialize D
            for nn = 0:(length(obj.Coeff)-1)

                for m = nn:(length(obj.Coeff)-1)
                    % This part calculates the      m >= n    elements
                    % LaguerreLnk = L_{n}^{k}(x)  where x = |a|^2 in our case.  
                    % We compute  L_{n}^{k}(x)  using the matrix we defined in the constructor    " function obj = Nbasis(n,N_space) "  .
                    LaguerreLnk = sum( times(obj.Lnk(nn+1,1:(nn+1),m-nn+1), abs(a).^(2*(0:nn))) );

                    D(m+1,nn+1) = sqrt(factorial(nn)/factorial(m))* a^(m-nn) *exp(-1/2 *abs(a)^2) *...
                                    LaguerreLnk;
                end

                for m = 0:(nn-1)
                    % This part calculates the      m < n    elements
                    % LaguerreLnk = L_{n}^{k}(x)  where x = |a|^2 in our case.  
                    % We compute  L_{n}^{k}(x)  using the matrix we defined in the constructor    " function obj = Nbasis(n,N_space) "  .
                    LaguerreLnk = sum( times(obj.Lnk(m+1,1:(m+1),nn-m+1), abs(a).^(2*(0:m))) );

                    D(m+1,nn+1) = sqrt(factorial(m)/factorial(nn))* (-a')^(nn-m) *exp(-1/2 *abs(a)^2) *...
                                    LaguerreLnk;
                end
            end
            obj.Coeff = D*obj.Coeff;     % returns the displaced number state.
            obj.n     = find(obj.Coeff);

        end
        function [W] = WignerFunction(obj,x_max,N)    % computes the Wigner distribution

            obj   = normalize(obj);                   % normalizes |ψ>
            x     = linspace(-x_max,x_max,N);         % x-component of the "indipendent variable"
            y     = x;                                % y-component of the "indipendent variable"
            W     = zeros(N,N);                       % initializes the matrix
            for i = 1:N
                for j = 1:N
                    b      = x(j) + 1i*y(i);          % indipendent variable of Wigner function
                    W(i,j) =  (obj.D_(-b).Coeff)' * obj.Parity * (obj.D_(-b).Coeff);     % W(b) = Trace(\rho * D(-b) * Parity * D(-b))
                    % 
                    %   <ψ|D†(-b) = (obj.D_(-b).Coeff)' 
                    %    D(-b)|ψ> = (obj.D_(-b).Coeff)
                end  
            end
            W = W * 2/pi;   % Wigner Distribution (ready for plot)

        end
        
        function [average_num, P] = PhotonNumber(obj)   % This function calculates the Photon Number in the state.
            obj = normalize(obj); 
            P   = times(conj(obj.Coeff),obj.Coeff)';       % Photon distribution P(n)
            average_num = P * (0:(length(obj.Coeff)-1))';  % <n> average photon number

        end
    end
end
