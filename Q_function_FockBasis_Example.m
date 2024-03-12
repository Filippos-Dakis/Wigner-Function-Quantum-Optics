% Filippos Tzimkas-Dakis   Virginia Tech  MARCH 2024
%
% Any feedback and suggestions are much appreciated! 
%    
%     ----->  dakisfilippos@vt.edu  <-------
% 
% This script was developed on MATLAB 2023a
%
% Example on Q-Husimi Function
%
% You should better execute this script upon running FockBasis_Example_1.m 
% and CoherentBasis_Examlpe_1.m . 
%
% The main feature of this script is to calculate and plot the Q-fuction of
% a quantum state written in the Fock basis. For sake of comparison,
% next to Q-fuction we also plot the Winger function.
% This script produces one figure.
% Figure 1: Wigner and Q fucntion of the quantum state |ψ>  = |0> + |1> 
%
% The runtime of this script is ~4 seconds on a gaming laptop.
% N  and N_hilbert affect the run time. 
%
% Version V 1.2.2
%%
warning('off','MATLAB:nchoosek:LargeCoefficient');     % 
close all
clear all
clc
tic

x_max = 3.5;                                    % needed for the square grid
N     = 100;                                    % N points across each direction  (affects the run time)
[X,Y] = meshgrid(linspace(-x_max,x_max,N));     % grid
dxdy  = (X(1,2) - X(1,1)) * (Y(2,1) - Y(1,1));  % surface differential

N_hilbert  = 15;                                % truncates the Hilbert space up to first 15 states (affects the run time)
n_2        = FockBasis([1;1],N_hilbert);        % |n_2> = |0> + |1>
n_2        = n_2.normalize;                     % normalize the state
  
Q_2     = n_2.Q_function(x_max,N);              % Q-function for even cat
check_Q = sum(sum(Q_2)) * dxdy;                 % the integral of Q function must be 1
fprintf(['\n ∫∫Q*dada^* = ',num2str(check_Q,3),'\n'])

W_2 = n_2.WignerFunction(x_max,N);              % Wigner function for even cat
W_2 = real(W_2);
check_W    = sum(sum(W_2)) * dxdy;       % the integral of W function must be 1
fprintf(['\n ∫∫W*dada^* = ',num2str(check_W,3)])
if ( check_W<=0.97 || check_W>= 1.05)
    fprintf(' ≠ 1 \n\n')
    fprintf(' Increase the Hilbert space ( N_hilbert ) \n to get a better approximation in Wigner Distribution !\n\n\n')
else
    fprintf('\n\n\n')
end


% ---- plots ---------------
cmap_1 = 'turbo';                               % plot color map
cmap_2 = 'turbo';                               % plot color map
face_alpha = 1;                                 % affects the black lines on the plot


M_2 = max(max(abs(Q_2)));                       % maximum value of Q function, used for display purposes


f1 = figure(1);

subplot(1,2,1)                                  % plot Q Function
[~, hc]     =   contourf(X,Y,Q_2,100,'EdgeAlpha',0,'FaceAlpha',face_alpha);
clim([-1 1]*M_2);
axis square
a1        =  gca;
a2        =  axes('Parent', f1, 'Position', a1.Position);
hs        =  surf(X, Y, Q_2, 'Parent', a2, 'EdgeColor','none');
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
zlim([-1 1]*M_2*1.1)
xlabel('Re(\beta)')
ylabel('Im(\beta)')
zlabel('Q_{even}(\beta)')
colormap(a1,cmap_1)
colormap(a2,cmap_2)
clim([-1 1]*M_2);
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
    Z1 = Q_2(i,:);
    plot3(x,Y1,Z1,'Color',[0,0,0,0.3]);    % [0,0,0,Transparency]  0,0,0 = black
end

% Plotting lines in the Y-Z plane

for i = 1:xspacing:length(x)
    X2 = x(i)*ones(size(y)); % a constant vector
    Z2 = Q_2(:,i);
    plot3(X2,y,Z2,'Color',[0,0,0,0.3]);
end

hold off
title('Q Function')
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 


M_2 = max(max(abs(W_2)));                    % maximum value of Wigner function, used for display purposes

subplot(1,2,2)                               % plots Wigner Function

[~, hc]     =   contourf(X,Y,W_2,100,'EdgeAlpha',0,'FaceAlpha',face_alpha);
clim([-1 1]*M_2);
axis square
a1        =  gca;
a2        =  axes('Parent', f1, 'Position', a1.Position);
hs        =  surf(X, Y, W_2, 'Parent', a2, 'EdgeColor','none');
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
zlim([-1 1]*M_2*1.1)
xlabel('Re(\beta)')
ylabel('Im(\beta)')
zlabel('W_{even}(\beta)')
colormap(a1,cmap_1)
colormap(a2,cmap_2)
clim([-1 1]*M_2);
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
    Z1 = W_2(i,:);
    plot3(x,Y1,Z1,'Color',[0,0,0,0.3]);    % [0,0,0,Transparency]  0,0,0 = black
end

% Plotting lines in the Y-Z plane

for i = 1:xspacing:length(x)
    X2 = x(i)*ones(size(y)); % a constant vector
    Z2 = W_2(:,i);
    plot3(X2,y,Z2,'Color',[0,0,0,0.3]);
end

hold off
title('Wigner Function')
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 

toc
