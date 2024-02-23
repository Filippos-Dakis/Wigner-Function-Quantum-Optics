% Filippos Tzimkas-Dakis   Virginia Tech  February 2024
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
% a quantum state written in the Coherent basis. For sake of comparison,
% next to Q-fuction we also plot the Winger function.
% This script produces one figure.
% Figure 1: Wigner and Q fucntion of the EVEN Cat state. 
%
% The runtime of this script is ~1 seconds on a gaming laptop.
%
% Version V 1.2
%%
close all
clear all
clc
tic

x_max = 6;                                      % needed for the square grid
N     = 400;                                    % N points across each direction
[X,Y] = meshgrid(linspace(-x_max,x_max,N));     % grid
dxdy  = (X(1,2) - X(1,1)) * (Y(2,1) - Y(1,1));  % surface differential


even_cat   = CoherentBasis([1;1],[3;-3]);       % create an even cat 

Q_even_cat = even_cat.Q_function(x_max,N);      % Q-function for even cat
check_Q    = sum(sum(Q_even_cat)) * dxdy;       % the integral of Q function must be 1


W_even_cat = even_cat.WignerFunction(x_max,N);  % Wigner function for even cat
W_even_cat = real(W_even_cat);
 



% ---- plots ---------------
cmap_1 = 'turbo';                               % plot color map
cmap_2 = 'turbo';                               % plot color map
face_alpha = 1;                                 % affects the black lines on the plot


M_2 = max(max(abs(Q_even_cat)));                % maximum value of Q function, used for display purposes


f1 = figure(1);

subplot(1,2,1)                                  % plot Q Function
[~, hc]     =   contourf(X,Y,Q_even_cat,100,'EdgeAlpha',0,'FaceAlpha',face_alpha);
clim([-1 1]*M_2);
axis square
a1        =  gca;
a2        =  axes('Parent', f1, 'Position', a1.Position);
hs        =  surf(X, Y, Q_even_cat, 'Parent', a2, 'EdgeColor','none');
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
    Z1 = Q_even_cat(i,:);
    plot3(x,Y1,Z1,'Color',[0,0,0,0.3]);    % [0,0,0,Transparency]  0,0,0 = black
end

% Plotting lines in the Y-Z plane

for i = 1:xspacing:length(x)
    X2 = x(i)*ones(size(y)); % a constant vector
    Z2 = Q_even_cat(:,i);
    plot3(X2,y,Z2,'Color',[0,0,0,0.3]);
end

hold off
title('Q Function')
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 


M_2 = max(max(abs(W_even_cat)));             % maximum value of Wigner function, used for display purposes

subplot(1,2,2)                               % plots Wigner Function

[~, hc]     =   contourf(X,Y,W_even_cat,100,'EdgeAlpha',0,'FaceAlpha',face_alpha);
clim([-1 1]*M_2);
axis square
a1        =  gca;
a2        =  axes('Parent', f1, 'Position', a1.Position);
hs        =  surf(X, Y, W_even_cat, 'Parent', a2, 'EdgeColor','none');
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
    Z1 = W_even_cat(i,:);
    plot3(x,Y1,Z1,'Color',[0,0,0,0.3]);    % [0,0,0,Transparency]  0,0,0 = black
end

% Plotting lines in the Y-Z plane

for i = 1:xspacing:length(x)
    X2 = x(i)*ones(size(y)); % a constant vector
    Z2 = W_even_cat(:,i);
    plot3(X2,y,Z2,'Color',[0,0,0,0.3]);
end

hold off
title('Wigner Function')
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 

toc