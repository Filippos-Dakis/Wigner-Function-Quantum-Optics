% Filippos Tzimkas-Dakis   Virginia Tech   MARCH 2024
%
% Any feedback and suggestions are much appreciated! 
%    
%     ----->  dakisfilippos@vt.edu  <-------
%
% Version V 1.2.2
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
        sqrt_n_over_m
        Parity  % Parity matrix diag(1,-1,1,-1 ....)
    end
    properties
        Coeff {mustBeNumeric}        % column vector that contains the coefficient of every Fock/Number state.
        Kets  {mustBeNumeric}        % column vector that contains the contributing Fock states.
        n     {mustBeNonnegative}    % column vector depicting the populated number states, ie n = [1 4 5] 
                                     % means that the |1> |4> and |5> have nonzero coefficients.  .
        N_Hilbert {mustBeNonnegative} % a number depicting the size of the Hilbert space
    end 
    methods
        
        function obj = FockBasis(nn,N_space) 
            % constructor of the Number Basis class
            % inputs:   n : [vector] column vector containing the coefficients of every state number state, ie n = [0 1 1i] = |1> + 1i|2> .
            %           N_space: [number] truncates the infinite Hilbert space to N_space levels. 
            if nargin==1
                N_space = 30;       % if N_space is not defined, is set to  N_space = 30  which is a good approximation           
            end  
            if ( find(nn, 1, 'last') + 5) >=N_space   % N_space must always be larger than length(n)
                N_space = length(nn) + 10;
            else
                nn = nn(1:find(nn, 1, 'last' ));
            end
      
            obj.Coeff = [nn(:);zeros(N_space -length(nn),1)]; % assign the input values to class properties
            obj.Kets  = (0:(N_space-1))';                % because the basis starts from |0>
            obj.n     = find(obj.Coeff)-1;                  % assign the populated number states
            obj.N_Hilbert = N_space;                        % 
            
            % ------------------------------------------------------------
            % Here we calculate coefficients of the associated Laguerre Polynomials,  L_{n}^{k} ,up to  m-th order.
            % for more info check   https://mathworld.wolfram.com/AssociatedLaguerrePolynomial.html
            % We will use then in the Displacement operator. To that end, we calculate them once and then use their values through the matrix   .
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
            % ------------------------------------------------------------
            % The following matrix "obj.sqrt_n_over_m" is a lower-triangular matrix that contains the values of sqrt(n!/m!) needed for
            % the Dispalcement operator written in the the Fock basis, see https://physics.stackexchange.com/questions/553225/representation-of-the-displacement-operator-in-number-basis.
			obj.sqrt_n_over_m = tril(zeros(N_space,N_space));
            for m = 0:(N_space-1)
                for nn = 0:m
                    obj.sqrt_n_over_m(m+1,nn+1) = sqrt(factorial(nn)/factorial(m));

                end
            end
            % -------------------------------------------------
        obj.Parity = sparse(diag(cos((0:(N_space-1))*pi)));   % Parity matrix
        % obj.Parity = diag(cos((0:(N_space-1))*pi));
        end



        function obj = normalize(obj)                % This function normalizes the state (the obj.Coeff vector).
            obj.Coeff = obj.Coeff/norm(obj.Coeff);
        end


        function [s1,s2] = print_state(obj,s)
            if nargin == 1
                s = '';
            else
                s = ['_',num2str(s)];
            end
            s1 = ['|n',s,'> = '];
            s2 = '';
            for ii = obj.n'+1
                coeff_str = compact_complex(obj.Coeff(ii));
                ket_str   = num2str(obj.Kets(ii));
    
                if real(obj.Coeff(ii))>1e-15 && imag(obj.Coeff(ii))>1e-15
                    s2 = [s2,'(',coeff_str,') |',ket_str,'>'];
                else
                    s2 = [s2, coeff_str,' |',ket_str,'>'];
                end

                if ii < max(obj.n)+1
                    s2 = [s2,'  +  '];
                end
            end
            s1 = [s1, s2];
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


        function obj = A(obj)           % Annhilation operator  A|n> = sqrt(n)|n-1>
            new_coeffs = sqrt(obj.Kets(2:end)).* obj.Coeff(2:end);
            obj        = FockBasis([new_coeffs; 0],obj.N_Hilbert);
        end


        function obj = A_dagger(obj)    % Creation operator  A_dagger|n> = sqrt(n+1)|n+1>
            new_coeffs = sqrt(obj.Kets(2:end)).* obj.Coeff(1:end-1);
            obj        = FockBasis([0;new_coeffs],obj.N_Hilbert);
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


        function obj = mtimes(scalar,obj)     % multiplies the coefficients in fornt of the kets with a scalar .
            obj.Coeff = scalar*obj.Coeff;
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

                    if m==nn % diagonal elements
                        LaguerreLnk = sum( times(obj.Lnk(nn+1,1:(nn+1),m-nn+1), abs(a).^(2*(0:nn))) );

                        D(m+1,nn+1) = obj.sqrt_n_over_m(m+1,nn+1) * a^(m-nn) *exp(-1/2 *abs(a)^2) *...
                                      LaguerreLnk;
                    else    
                        LaguerreLnk = sum( times(obj.Lnk(nn+1,1:(nn+1),m-nn+1), abs(a).^(2*(0:nn))) );
                        
                        % lower triangular matrix elements
                        D(m+1,nn+1) = obj.sqrt_n_over_m(m+1,nn+1) * a^(m-nn) *exp(-1/2 *abs(a)^2) *...
                                      LaguerreLnk;
                        % upper trangular matrix elements
                        D(nn+1,m+1) = obj.sqrt_n_over_m(m+1,nn+1)* (-a')^(m-nn) *exp(-1/2 *abs(a)^2) *...
                                      LaguerreLnk;

                        % coefficient = obj.sqrt_n_over_m(m+1,nn+1) * exp(-1/2 *abs(a)^2) * LaguerreLnk;
                        % D(m+1,nn+1) = coefficient * a^(m-nn);
                        % D(nn+1,m+1) = coefficient * (-a')^(m-nn);
                    end
                end

            end
            obj.Coeff = D*obj.Coeff;     % returns the displaced number state.
            obj.n     = find(obj.Coeff);
        end
        

        function obj = S_(obj,z) 
            % Squeeze operator function
            a     = sparse(diag(sqrt(obj.Kets(2:end)),+1));
            a_dag = a';
            obj.Coeff = expm(1/2 * (z'*a^2 - z*a_dag^2) ) * obj.Coeff;
            % size(obj.Coeff)
            % size(expm(1/2 * (z'*a^2 - z*a_dag^2) ) )
            % q = expm(1/2 * (z'*a^2 - z*a_dag^2) );
            % q(3,:)
            obj = FockBasis(obj.Coeff,obj.N_Hilbert);
        end



        function [W] = WignerFunction(obj,x_max,N)    % computes the Wigner distribution
            % Wigner Distribution
            % inputs : obj   = the object/state to calculate the Wigner function .
            %          x_max = maximum x (and y -- square grid y_max = x_max). W-function will be computed between [(-x_max,x_max),(-y_max,y_max)] .
            %          N     = number of points in each direction. The final matrix W will have dimensions of NxN.
            obj   = normalize(obj);                   % normalizes |ψ>
            x     = linspace(-x_max,x_max,N);         % x-component of the "indipendent variable"
            y     = x;                                % y-component of the "indipendent variable"
            W     = zeros(N,N);                       % initializes the matrix
            for i = 1:N
                for j = 1:N
                    b       = x(j) + 1i*y(i);           % indipendent variable of Wigner function
                    D_psi   = obj.D_(- b ).Coeff;         %  D(-b)|ψ> = (obj.D_(-b).Coeff)
                    % [D_psi, D_psi_coj] = obj.D__(- (x(j) + 1i*y(i)) );
                    % D_psi' = (obj.D_(-b).Coeff)' ;    % <ψ|D†(-b) = (obj.D_(-b).Coeff)'
                    
                    %   W(b) = Trace(\rho * D(-b) * Parity * D(-b)) = <ψ|D†(-b) * Π * D(-b)|ψ>   . 
                    W(i,j)  =  D_psi' * obj.Parity * D_psi; 
                    % W(i,j)  =  D_psi_coj * obj.Parity * D_psi; 
                    %
                    % W(i,j) =  (obj.D_(-b).Coeff)' * obj.Parity * (obj.D_(-b).Coeff);     % W(b) = Trace(\rho * D(-b) * Parity * D(-b))
                    % 
                    %  <ψ|D†(-b) = (obj.D_(-b).Coeff)' 
                    %   D(-b)|ψ> = (obj.D_(-b).Coeff)
                end  
            end
            W = W * 2/pi;   % Wigner Distribution (ready for plot)

        end
        


        function [Q] = Q_function(obj,x_max,N)
            % Husimi-Q function
            % Inputs : obj   = the object/state to calculate the Q-function .
            %          x_max = maximum x (and y -- square grid y_max = x_max). Q-function will be computed between [(-x_max,x_max),(-y_max,y_max)] .
            %          N     = number of points in each direction. The final matrix Q will have dimensions of NxN.
            % Outputs: Q = Husimi distribution, NxN matrix

            obj   = normalize(obj);                   % normalizes |ψ>
            k     = 0:(length(obj.Coeff)-1);          % size of Hilbert space
            k     = k(:);
            x     = linspace(-x_max,x_max,N);         % creates the meshgrid
            [X,Y] = meshgrid(x);                      % assings the grid to matrices
            B     = X + 1i*Y;                         % indipendent variable  W = W(b)
            Q     = zeros(N,N);                       % initializes the matrix
            
            for i = 1:N
                for j = 1:N

                    Q(i,j) = sum( ((B(i,j)').^(k))./sqrt(factorial(k)) .*obj.Coeff ) ;

                end
            end
            Q = 1/pi * abs( exp(-1/2 *abs(B).^2).*Q ).^2; % Q function (ready for plot)
            
        end


        
        function [average_num, P] = PhotonNumber(obj)      % This function calculates the Photon Number in the state.
            obj = normalize(obj); 
            P   = times(conj(obj.Coeff),obj.Coeff)';       % Photon distribution P(n)
            average_num = P * (0:(length(obj.Coeff)-1))';  % <n> average photon number

        end

        function [D_psi, D_psi_conj] = D__(obj,a)
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

                    if m==nn % diagonal elements
                        LaguerreLnk = sum( times(obj.Lnk(nn+1,1:(nn+1),m-nn+1), abs(a).^(2*(0:nn))) );

                        D(m+1,nn+1) = obj.sqrt_n_over_m(m+1,nn+1) * a^(m-nn) *exp(-1/2 *abs(a)^2) *...
                                      LaguerreLnk;
                    else    
                        LaguerreLnk = sum( times(obj.Lnk(nn+1,1:(nn+1),m-nn+1), abs(a).^(2*(0:nn))) );
                        
                        % lower triangular matrix elements
                        D(m+1,nn+1) = obj.sqrt_n_over_m(m+1,nn+1) * a^(m-nn) *exp(-1/2 *abs(a)^2) *...
                                      LaguerreLnk;
                        % upper trangular matrix elements
                        D(nn+1,m+1) = obj.sqrt_n_over_m(m+1,nn+1)* (-a')^(m-nn) *exp(-1/2 *abs(a)^2) *...
                                      LaguerreLnk;

                        % coefficient = obj.sqrt_n_over_m(m+1,nn+1) * exp(-1/2 *abs(a)^2) * LaguerreLnk;
                        % D(m+1,nn+1) = coefficient * a^(m-nn);
                        % D(nn+1,m+1) = coefficient * (-a')^(m-nn);
                    end
                end

            end
            D_psi = D*obj.Coeff;     % returns the displaced number state.
            D_psi_conj = D_psi';
        end

        
    end
end
