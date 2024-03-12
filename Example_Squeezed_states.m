% Filippos Tzimkas-Dakis   Virginia Tech  February 2024
%
% Any feedback and suggestions are much appreciated! 
%    
%     ----->  dakisfilippos@vt.edu  <-------
% 
% This script was developed on MATLAB 2023a
%
% Example on Squeezed States using the FockBasis class
% You should better execute its section one after the other while also reading the comments.
% In this way, you will understand the Squeezed state properties and fetuares.
%
% This script produces four frigures.
% Figure 1: population distribution for the the Squeezed vacuum state S(z)|0>
% Figure 2: Q function for vacuum |0> and Squeezed vacuum S(z)|0> states to compare and contrast 
% Figure 3: population distribution for the the Squeezed coherent state D(a)S(z)|0>
% Figure 4: Q function for coherent |a> and Squeezed coherent D(a)S(z)|0> states to compare and contrast 
%
% The runtime of this script is ~1 seconds on a gaming laptop.
%
% Version V 1.2.2

%% ----------------------- Section 1 --------------------------------------
% Define your quantum states using the Fock basis
warning('off','MATLAB:nchoosek:LargeCoefficient');     % 

close all
clear all
clc
tic
% First, we define the ground state in the Fock/Number basis
fprintf('\n\n --------------- Section 1 ----------------------\n')

N_hilbert = 25;                 % truncates the Hilbert space up to first 25 states

c0  = 2+2i;                     % coefficient in front of the \ket
n_0 = FockBasis(c0,N_hilbert);  % |n_1> = (2+2i)|0>
fprintf(['\n Before normalization   |n_1> = (',num2str(n_0.Coeff(1)),')|',num2str(n_0.Kets(1)),'>\n'])
n_0 = n_0.normalize;            % normalize our state
n_0.Coeff;                      % display normalized coefficients 
%                                 n_1.Coeff column vector of length = 20
fprintf(['\n After normalization    |n_1> = (',num2str(n_0.Coeff(1)),')|',num2str(n_0.Kets(1)),'>\n'])

%% ----------------------- Section 2 --------------------------------------
fprintf('\n\n --------------- Section 2 ----------------------\n')
% Here we create the Squeezed Vaccum state 
% 
% Recall the SQUEEZE OPERATOR       S(z) = exp(1/2 * (z'*a^2 - z*(a†)^2) )
%                                  

r     = 1;                      % Squeeze magnitude
theta = 0;                      % Squeeze angle
z     = r * exp(1i*theta);      % total Squeeze parameter

Squeezed_Vac = n_0.S_(z);       % Squeezed Vacuum state
Squeezed_Vac.Coeff;             % Print the coefficients that multiply each Fock state
Squeezed_Vac.n;                 % Print the kets that have NON-zero Coefficients

% Photon number 
[N_sq,P_sq] = Squeezed_Vac.PhotonNumber;     % N_ = average photon number
%                                              P_n = photon distribution as a funtion of n


% ---- PLOTS ---------------------------------
f1 = figure(1);
bar(0:length(P_sq)-1,P_sq,'EdgeColor','none')
text(14,max(P_sq)*0.9,sprintf('$\\langle n \\rangle $= %.2f',N_sq), ...
    'Interpreter', 'latex', 'fontsize', 10);
xlabel('n')
ylabel('P(n)')
title('$|\psi\rangle = S(z)|0\rangle$','Interpreter','latex')
f1.Units = 'normalized';
f1.OuterPosition = [0.1273    0.6131    0.3086    0.2850];
% ---- END OF PLOTS ---------------------------




x_max = 5;                                    % needed for the square grid
N     = 151;                                    % N points across each direction  (affects the run time)
[X,Y] = meshgrid(linspace(-x_max,x_max,N));     % grid
dxdy  = (X(1,2) - X(1,1)) * (Y(2,1) - Y(1,1));  % surface differential


% -----  Vacuum state ----------------------
Q_0     = n_0.Q_function(x_max,N);              % Q-function for even cat
check_Q = sum(sum(Q_0)) * dxdy;                 % the integral of Q function must be 1
fprintf(['\n ∫∫Q_0*dada^* = ',num2str(check_Q,3),'\n'])
% 
% ----- Squeezed Vacumm state ---------------
Q_sq    = Squeezed_Vac.Q_function(x_max,N);              % Q-function for even cat
check_Q = sum(sum(Q_sq)) * dxdy;                         % the integral of Q function must be 1
fprintf(['\n ∫∫Q_sq*dada^* = ',num2str(check_Q,3),'\n\n'])



% ----------------- PLOTS -------------------------------------------------
cmap_1 = 'turbo';                               % plot color map
cmap_2 = 'turbo';                               % plot color map
face_alpha = 1;                                 % affects the black lines on the plot

M = max(max(abs(Q_0)));                         % maximum value of Q function, used for display purposes


f2 = figure(2);

subplot(1,2,1)                                  % plot Q Function
[~, hc]     =   contourf(X,Y,Q_0,100,'EdgeAlpha',0,'FaceAlpha',face_alpha);
clim([-1 1]*M);
axis square
a1        =  gca;
a2        =  axes('Parent', f2, 'Position', a1.Position);
hs        =  surf(X, Y, Q_0, 'Parent', a2, 'EdgeColor','none');
a1.Color  =  'none';
a2.Color  =  'none';
a1.ZLim   =  [0 1];
a2.ZLim   =  [-9 9];
a1.XTick  =  [];
a1.YTick  =  [];
a1.ZTick  =  [];
a1.Box    =  'off';
a2.Box    =  'off';
a2.View   = [-47.8105, 25.4282];
% Call after setting desired view on a2 (surf plot) 
a1.View   =   a2.View;
zlim([-1 1]*M*1.1)
xlabel('Re(\beta)')
ylabel('Im(\beta)')
zlabel('Q_{0}(\beta)')
colormap(a1,cmap_1)
colormap(a2,cmap_2)
clim([-1 1]*M);
axis square
%  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
x=X(1,:);
y=Y(:,1);
%               Divide the lengths by the number of lines needed
xnumlines = 14;                %               10 lines
ynumlines = 14;                %               10 partitions
xspacing = round(length(x)/xnumlines);
yspacing = round(length(y)/ynumlines);

% Plot the mesh lines 
% Plotting lines in the X-Z plane
hold on

for i = 1:yspacing:length(y)
    Y1 = y(i)*ones(size(x));    % a constant vector
    Z1 = Q_0(i,:);
    plot3(x,Y1,Z1,'Color',[0,0,0,0.3]);    % [0,0,0,Transparency]  0,0,0 = black
end

% Plotting lines in the Y-Z plane

for i = 1:xspacing:length(x)
    X2 = x(i)*ones(size(y)); % a constant vector
    Z2 = Q_0(:,i);
    plot3(X2,y,Z2,'Color',[0,0,0,0.3]);
end

hold off
title('Q Function, Vacuum state')
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 


M = max(max(Q_sq));                     % maximum value of Q function, used for display purposes

subplot(1,2,2)                               % plots Q fucntion of Squeezed Vacuum state

[~, hc]     =   contourf(X,Y,Q_sq,100,'EdgeAlpha',0,'FaceAlpha',face_alpha);
clim([-1 1]*M);
axis square
a1        =  gca;
a2        =  axes('Parent', f2, 'Position', a1.Position);
hs        =  surf(X, Y, Q_sq, 'Parent', a2, 'EdgeColor','none');
a1.Color  =  'none';
a2.Color  =  'none';
a1.ZLim   =  [0 1];
a2.ZLim   =  [-9 9];
a1.XTick  =  [];
a1.YTick  =  [];
a1.ZTick  =  [];
a1.Box    =  'off';
a2.Box    =  'off';
a2.View   = [-47.8105, 25.4282];
% Call after setting desired view on a2 (surf plot) 
a1.View   =   a2.View;
zlim([-1 1]*M*1.1)
xlabel('Re(\beta)')
ylabel('Im(\beta)')
zlabel('Q_{Sq}(\beta)')
colormap(a1,cmap_1)
colormap(a2,cmap_2)
clim([-1 1]*M);
axis square
%  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
x=X(1,:);
y=Y(:,1);
%               Divide the lengths by the number of lines needed
xnumlines = 14;                %               10 lines
ynumlines = 14;                %               10 partitions
xspacing = round(length(x)/xnumlines);
yspacing = round(length(y)/ynumlines);

% Plot the mesh lines 
% Plotting lines in the X-Z plane
hold on

for i = 1:yspacing:length(y)
    Y1 = y(i)*ones(size(x));    % a constant vector
    Z1 = Q_sq(i,:);
    plot3(x,Y1,Z1,'Color',[0,0,0,0.3]);    % [0,0,0,Transparency]  0,0,0 = black
end

% Plotting lines in the Y-Z plane

for i = 1:xspacing:length(x)
    X2 = x(i)*ones(size(y)); % a constant vector
    Z2 = Q_sq(:,i);
    plot3(X2,y,Z2,'Color',[0,0,0,0.3]);
end

hold off
title('Q Fucntion, Squeezed Vacuum state')
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
f2.Units = 'normalized';
f2.OuterPosition = [0.1148    0.1044    0.3367    0.4756];
% ------ END OF PLOTS -----------------------------------------------------



%% ----------------------- Section 3 --------------------------------------
fprintf('\n\n --------------- Section 3 ----------------------\n')
% Here we create a Squeezed Coherent state
%                                  
% Recal that  |alpha,z> = D(alpha)*S(z)|0>
% 

r     = 1;                      % Squeeze magnitude
theta = 0;                      % Squeeze angle
z     = r * exp(1i*theta);      % total Squeeze parameter

alpha = -2 - 2i;                 % Displacement operator coefficient

Squeezed_Coh = Squeezed_Vac.D_(alpha);       % Squeezed Vacuum state Coherent state
Squeezed_Coh.Coeff;                          % Print the coefficients that multiply each Fock state
Squeezed_Coh.n;                              % Print the kets that have NON-zero Coefficients

Coh_state = n_0.D_(alpha);

% Photon number 
[N_sq_coh,P_sq_coh] = Squeezed_Coh.PhotonNumber;     % N_ = average photon number
%                                                      P_n = photon distribution as a funtion of n


% ---- PLOTS ---------------------------------
f3 = figure(3);
bar(0:length(P_sq_coh)-1,P_sq_coh,'EdgeColor','none')
text(14,max(P_sq_coh)*0.9,sprintf('$\\langle n \\rangle $= %.2f',N_sq_coh), ...
    'Interpreter', 'latex', 'fontsize', 10);
xlabel('n')
ylabel('P(n)')
title('$|\alpha,z\rangle = D(\alpha)S(z)|0\rangle$','Interpreter','latex')
f3.Units = 'normalized';
f3.OuterPosition = [ 0.5625    0.6144    0.3086    0.2850];
% ---- END OF PLOTS ---------------------------


% -----  Vacuum state ----------------------
Q_a     = Coh_state.Q_function(x_max,N);              % Q-function for even cat
check_Q = sum(sum(Q_a)) * dxdy;                 % the integral of Q function must be 1
fprintf(['\n ∫∫Q_0*dada^* = ',num2str(check_Q,3),'\n'])
% 
% ----- Squeezed Vacumm state ---------------
Q_sq_coh    = Squeezed_Coh.Q_function(x_max,N);              % Q-function for even cat
check_Q = sum(sum(Q_sq_coh)) * dxdy;                         % the integral of Q function must be 1
fprintf(['\n ∫∫Q_sq*dada^* = ',num2str(check_Q,3),'\n\n'])



% ----------------- PLOTS -------------------------------------------------
cmap_1 = 'turbo';                               % plot color map
cmap_2 = 'turbo';                               % plot color map
face_alpha = 1;                                 % affects the black lines on the plot

M = max(max(abs(Q_a)));                         % maximum value of Q function, used for display purposes


f4 = figure(4);

subplot(1,2,1)                                  % plot Q Function
[~, hc]     =   contourf(X,Y,Q_a,100,'EdgeAlpha',0,'FaceAlpha',face_alpha);
clim([-1 1]*M);
axis square
a1        =  gca;
a2        =  axes('Parent', f4, 'Position', a1.Position);
hs        =  surf(X, Y, Q_a, 'Parent', a2, 'EdgeColor','none');
a1.Color  =  'none';
a2.Color  =  'none';
a1.ZLim   =  [0 1];
a2.ZLim   =  [-9 9];
a1.XTick  =  [];
a1.YTick  =  [];
a1.ZTick  =  [];
a1.Box    =  'off';
a2.Box    =  'off';
a2.View   = [-47.8105, 25.4282];
% Call after setting desired view on a2 (surf plot) 
a1.View   =   a2.View;
zlim([-1 1]*M*1.1)
xlabel('Re(\beta)')
ylabel('Im(\beta)')
zlabel('Q_{0}(\beta)')
colormap(a1,cmap_1)
colormap(a2,cmap_2)
clim([-1 1]*M);
axis square
%  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
x=X(1,:);
y=Y(:,1);
%               Divide the lengths by the number of lines needed
xnumlines = 14;                %               10 lines
ynumlines = 14;                %               10 partitions
xspacing = round(length(x)/xnumlines);
yspacing = round(length(y)/ynumlines);

% Plot the mesh lines 
% Plotting lines in the X-Z plane
hold on

for i = 1:yspacing:length(y)
    Y1 = y(i)*ones(size(x));    % a constant vector
    Z1 = Q_a(i,:);
    plot3(x,Y1,Z1,'Color',[0,0,0,0.3]);    % [0,0,0,Transparency]  0,0,0 = black
end

% Plotting lines in the Y-Z plane

for i = 1:xspacing:length(x)
    X2 = x(i)*ones(size(y)); % a constant vector
    Z2 = Q_a(:,i);
    plot3(X2,y,Z2,'Color',[0,0,0,0.3]);
end

hold off
title('Q Function, Vacuum state')
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 


M = max(max(Q_sq_coh));                     % maximum value of Q function, used for display purposes

subplot(1,2,2)                               % plots Q fucntion of Squeezed Vacuum state

[~, hc]     =   contourf(X,Y,Q_sq_coh,100,'EdgeAlpha',0,'FaceAlpha',face_alpha);
clim([-1 1]*M);
axis square
a1        =  gca;
a2        =  axes('Parent', f4, 'Position', a1.Position);
hs        =  surf(X, Y, Q_sq_coh, 'Parent', a2, 'EdgeColor','none');
a1.Color  =  'none';
a2.Color  =  'none';
a1.ZLim   =  [0 1];
a2.ZLim   =  [-9 9];
a1.XTick  =  [];
a1.YTick  =  [];
a1.ZTick  =  [];
a1.Box    =  'off';
a2.Box    =  'off';
a2.View   = [-47.8105, 25.4282];
% Call after setting desired view on a2 (surf plot) 
a1.View   =   a2.View;
zlim([-1 1]*M*1.1)
xlabel('Re(\beta)')
ylabel('Im(\beta)')
zlabel('Q_{SqCoh}(\beta)')
colormap(a1,cmap_1)
colormap(a2,cmap_2)
clim([-1 1]*M);
axis square
%  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
x=X(1,:);
y=Y(:,1);
%               Divide the lengths by the number of lines needed
xnumlines = 14;                %               10 lines
ynumlines = 14;                %               10 partitions
xspacing = round(length(x)/xnumlines);
yspacing = round(length(y)/ynumlines);

% Plot the mesh lines 
% Plotting lines in the X-Z plane
hold on

for i = 1:yspacing:length(y)
    Y1 = y(i)*ones(size(x));    % a constant vector
    Z1 = Q_sq_coh(i,:);
    plot3(x,Y1,Z1,'Color',[0,0,0,0.3]);    % [0,0,0,Transparency]  0,0,0 = black
end

% Plotting lines in the Y-Z plane

for i = 1:xspacing:length(x)
    X2 = x(i)*ones(size(y)); % a constant vector
    Z2 = Q_sq_coh(:,i);
    plot3(X2,y,Z2,'Color',[0,0,0,0.3]);
end

hold off
title('Q Fucntion, Squeezed Coherent state')
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
f4.Units = 'normalized';
f4.OuterPosition = [0.5473    0.1081    0.3367    0.4756];
% ------ END OF PLOTS -----------------------------------------------------


toc










