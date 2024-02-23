% Filippos Tzimkas-Dakis   Virginia Tech   February 2024
%
% Version V 1.2
% 
% The present script defines a class whose objects behave as coherent states.
% Every object of this class has all the basic features of a quantum coherent state .
% Namely, one can define a state, normalize it, compute overlaps between states, apply, 
% add quantum states,  displacement operators, compute the Wigner function, calculate 
% the photon distribution, average photo number.
%
% Please have a look at the examples accompanying this script 
%
% I might add more features in the future.

classdef CoherentBasis
    properties
        Coeff {mustBeNumeric}    % column vector containing the coefficient of every Coherent state.
        Kets  {mustBeNumeric}    % column vector containing the coherent state arguments, namely
                                 % [0; 2.4; 1+1i]  ------->  |0> + |2.4> +  |1+1i>
    end
    methods
        function obj = CoherentBasis(c,s)
            % constructor of the CoherentBasis class
            % inputs:   c : [vector] column vector containing the coefficients of every coherent state, ie c = [1 1i] = |s(1)> + 1i|s(s)> .
            %           s : [vector] column vector containing the arguments of input coherent states, ie s= [0; 2.4; 1+1i]  ------->  |0> + |2.4> +  |1+1i>   s
            %     
            if nargin==2
                obj.Coeff = c(:);
                obj.Kets  = s(:);
            end
        end

        function obj = normalize(obj)                % normalizes the input state
            obj.Coeff = obj.Coeff/sqrt(braket(obj));
        end
        function q = braket(obj,objj)   % computes <\psi_1|\psi_2>
            if (nargin==1)              % if nargin=1 computes the 1/N normalization factor
                objj = obj;
            end
            q         = 0;
            objj.Coeff = conj(objj.Coeff);
            % The following loop computes all the dot products between <\psi_1|psi_2> = (k_1*<a_1| + k_2*<a_2| + ...+ k_n*<a_n|) (c_1*|b_1> + c_2*|b_2> + .... + c_m*|b_m>) .    p
            for m = 1:length(obj.Coeff(:,1))
                for n = 1:length(objj.Coeff(:,1))
                    % we use the fact that  <a|b>  = exp(-1/2 *(|a|^2 + |b|^2 - 2*conj(a)*b)) 
                    q = q + obj.Coeff(m)*objj.Coeff(n)* exp(-1/2 *( abs(obj.Kets(m))^2 + abs(objj.Kets(n))^2 - 2*obj.Kets(m)'*objj.Kets(n) ));
                end
            end
        end

        function obj = D_(obj,a)                       % displacement operator D(b)|a> = exp(i*Im(a^*·b))|a+b>
            obj.Coeff = times(obj.Coeff,exp(1i*imag(a * obj.Kets) ));  % computes the phase needed
            obj.Kets = obj.Kets + a;                                   % displaces the argument
        end

        function [W] = WignerFunction(obj,x_max,N)    % computes the Wigner distribution
            % Wigner Distribution
            % inputs : obj   = the object/state to calculate the Wigner function .
            %          x_max = maximum x (and y -- square grid y_max = x_max). W-function will be computed between [(-x_max,x_max),(-y_max,y_max)] .
            %          N     = number of points in each direction. The final matrix W will have dimensions of NxN.
            obj   = normalize(obj);                   % normalizes |ψ>

            x     = linspace(-x_max,x_max,N);         % creates the meshgrid
            [X,Y] = meshgrid(x);                      % assings the grid to matrices
            B     = X + 1i*Y;                         % indipendent variable  W = W(b)
            W     = zeros(N,N);                       % initializes the matrix

            for m = 1:length(obj.Coeff)
                for n = 1:length(obj.Coeff)
                    % we compute all the possible combinations of     coeff(a_2^*)*coeff(a_1) <a_2|D†(-b) Π D(-b)|a_1> .
                    % In the following formula we have already did some algebra, in detail we used Parity's properties to obtain   .
                    % <a_2|D†(-b) Π D(-b)|a_1> = <a_2|D†(-b) D(b) Π |a_1> = <a_2|D†(-b) D(b)|- a_1> = 
                    % <a_2|D†(-b) = exp(1i * Im(conj(a_2) * b)) <a_2 -b|
                    % D(b)|- a_1> = exp(1i * Im(conj(-a_1) * b)) |b-a_1>
                    % <a_2 -b|b-a_1> = exp(-1/2 *(|a_2 - b|^2 + |b - a_1|^2  - 2*coj(a_2 - b)*(a_1 - b)  )  )
                    W = W + obj.Coeff(m)'*obj.Coeff(n)*...
                        exp(1i *( imag(obj.Kets(m)'*B) + imag(-obj.Kets(n)'*B) ) ) .*...
                        exp(-1/2 *( abs(obj.Kets(m) - B).^2 + abs(obj.Kets(n) - B).^2 - 2*conj(obj.Kets(m) - B).*(B-obj.Kets(n)) ) );
                end
            end
            W = W * 2/pi;  % Wigner function (ready for plot)
        end

        function [Q] = Q_function(obj,x_max,N)
            % Husimi-Q function
            % Inputs : obj   = the object/state to calculate the Q-function .
            %          x_max = maximum x (and y -- square grid y_max = x_max). Q-function will be computed between [(-x_max,x_max),(-y_max,y_max)] .
            %          N     = number of points in each direction. The final matrix Q will have dimensions of NxN.
            % Outputs  Q     = Husimi distribution, NxN matrix
            obj   = normalize(obj);                   % normalizes |ψ>

            x     = linspace(-x_max,x_max,N);         % creates the meshgrid
            [X,Y] = meshgrid(x);                      % assings the grid to matrices
            B     = X + 1i*Y;                         % indipendent variable  W = W(b)
            Q     = zeros(N,N);                       % initializes the matrix

            for m = 1:length(obj.Coeff)

                Q = Q + obj.Coeff(m) * exp(-1/2 *( abs(B).^2 + abs(obj.Kets(m)).^2 - 2*conj(B).*(obj.Kets(m)) ) );
            
            end
            
            Q = 1/pi * abs(Q).^2;  % Q function (ready for plot)

        end

        function obj = plus(obj1,obj2)          % adds two vectors and combine the repeated states.
            coeff = [obj1.Coeff; obj2.Coeff];
            kets  = [obj1.Kets; obj2.Kets];
            
            [kets,indx] = sort(kets);   % we sort "kets" to make for coding reasons only
            coeff      = coeff(indx);   % sort the "coeff"s accordingly
            
            j        = 1;
            kets2(j) = kets(1);
            coeff2(j) = coeff(1);
            for i = 1:(length(kets)-1)                  % This loop locates the "reapeated" coherent states and add their coefficients
                if kets(i) == kets(i+1)                 % For instance, if we have 2*|a_1> + (3+1i)|a_1> = (5 + 1i)|a_1>  in the final state
                    coeff2(j) = coeff2(j) + coeff(i+1); % Taking care of the above cases speeds-up our code.
                else
                    j = j+1;
                    kets2(j) = kets(i+1);
                    coeff2(j) = coeff(i+1);
                end
            end
            obj = CoherentBasis(coeff2(:),kets2(:));   % final state/object
        end

        function [average_num, Photon_Distribution] = PhotonNumber(obj) 
            % This function calculates the average photon in the state and the photon distribution.
            % Inputs:  obj = object of the class  .
            % Outputs: average_num = average photon number
            %          Photon_distribution = number of photons in every Fock State |n>      .
            %          
            obj       = normalize(obj);        % normalizes |ψ>
            N_hilbert = 30;                    % truncates the infinte Hilbert to first 30 Fock States 
            Fock_state = zeros(N_hilbert,1);
            for i = 1:N_hilbert
                % projection on Fock states   <n|a> = exp(-1/2 |a|^2) (a^n)/sqrt(n!)        .
                Fock_state(i) = sum( times(obj.Coeff, times(exp(-1/2 * abs(obj.Kets).^2), obj.Kets.^(i-1))) )/sqrt(factorial(i-1));
            end
            Photon_Distribution = abs(Fock_state).^2;                   % P(n) = |<n|ψ>|^2     .
            average_num = sum(Photon_Distribution'.*(0:(N_hilbert-1))); % <n> = Sum_n P(n)*n   .

        end

    end
end
