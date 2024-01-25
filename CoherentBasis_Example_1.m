% Filippos Tzimkas-Dakis   Virginia Tech  January 2024
%
% Any feedback and suggestions are much appreciated! 
%    
%     ----->  dakisfilippos@vt.edu  <-------
% 
% This script was developed on MATLAB 2023a
%
% Example on CoherentBasis class
% You should better execute its section one after the other while also reading the comments.
% In this way, you will understand the properties, the methods and the capabilities of the 
% Coherent Basis class.
%
% This script produces four frigures.
% Figure 1: population distribution in Fock basis for four different states written in the coherent basis
% Figure 2: population distribution in Fock basis for an even and an odd cat state 
% Figure 3: Wigner functions for (i) vacumm |0> , (ii) coherent |\alpha>,  (iii) superposition of two coherent states |\alpha> + |\beta> 
% Figure 4: Wigner functions for (i) even cat state, (ii) odd cat state , (iii) compass cat state 
%
% The runtime of this script is ~5 seconds on a gaming laptop.
%
% Version V 1.1

%% Define your quantum states |psi_1>, |psi_2> using coherent states
close all
clear all
clc
tic
% First we define a simple coherent state of argument \alpha = 1
c1 = 2+2i;                      % coefficient in front of the \ket
s1 = 1;                         % argument inside the \ket 
psi_1 = CoherentBasis(c1,s1);   % |psi_1> = (2+2i)|1>,  here |1> = |α=1> is a coherent state, NOT a Fock state!
psi_1 = psi_1.normalize;        % normalize our state
psi_1.Coeff                     % display normalized coefficients

% Now let's define a more general state
c2 = [1;  2];                   % coefficient in front of the \ket
s2 = [-1; 3i];                  % argument inside the \ket 
psi_2 = CoherentBasis(c2,s2);   % |psi_2> = |-1> + 2|3i>, here |-1>, |3i> are coherent states
psi_2 = psi_2.normalize;        % normalize the state
psi_2.Coeff                     % display normalized coefficients

%% Dot product, quantum state addition, Displacement operator

% Dot probucts 
a = braket(psi_1);              % a and b are the same
b = braket(psi_1,psi_1);        % 

c1 = braket(psi_1,psi_2);       % c1 = <psi_1|psi_2>
c2 = braket(psi_2,psi_1);       % c2 = <psi_2|psi_1> = conj(c1)

% Add two quantum states or two objects
psi_3 = psi_1 + psi_2;          % |psi_3> = |psi_1> + |psi_2>    NOT NORMALIZED
psi_3.Coeff                     % display coefficients
psi_3 = psi_3.normalize;        % normalize the new state
psi_3.Coeff                     % display normalized coefficients 
d = braket(psi_3)               % d = <psi_3|psi_3>

psi_4 = psi_1 + psi_1;          % |psi_4> = 2*psi_1 = 2*psi_1.Coeff * |psi_1.Kets> 

% Displacement operator
psi_0 = CoherentBasis(1,0);     % |psi_0> = |0>
psi_0.Kets                      % display the \kets of this state
z     = 1 + 2i;                 % argument of Displacement operator
psi_z = psi_0.D_(z);            % |psi_z> = D(z)|psi_0> = |z>                 
psi_z.Kets                      % display the \kets of this state  

%% Photon number 
[N_0,P_n0] = psi_0.PhotonNumber;     % N_ = average photon number
[N_1,P_n1] = psi_1.PhotonNumber;     % P_n = photon distribution as a funtion of n
[N_2,P_n2] = psi_2.PhotonNumber;
[N_3,P_n3] = psi_3.PhotonNumber;

n = 0:length(P_n0)-1;                % for plot reasons
% ---- plots ---------------------------------
if 1 
figure(1)
subplot(2,2,1)
bar(n,P_n0,'EdgeColor','none')
text(5,0.75,sprintf('$\\langle n \\rangle $= %d',N_0), ...
    'Interpreter', 'latex', 'fontsize', 10);
xlabel('n')
ylabel('P(n)')
title('$|\psi_0\rangle = |0\rangle$','Interpreter','latex')
subplot(2,2,2)
bar(n,P_n1,'EdgeColor','none')
text(5,0.3,sprintf('$\\langle n \\rangle $= %d',N_1), ...
    'Interpreter', 'latex', 'fontsize', 10);
xlabel('n')
ylabel('P(n)')
title('$|\psi_1\rangle = |1\rangle$','Interpreter','latex')
subplot(2,2,3)
bar(n,P_n2,'EdgeColor','none')
text(20,0.08,sprintf('$\\langle n \\rangle $= %.2f',N_2), ...
    'Interpreter', 'latex', 'fontsize', 10);
xlabel('n')
ylabel('P(n)')
title('$|\psi_2\rangle \propto |-1\rangle + 2|3i\rangle$','Interpreter','latex')
subplot(2,2,4)
bar(n,P_n3,'EdgeColor','none')
text(20,0.3,sprintf('$\\langle n \\rangle $= %.2f',N_3), ...
    'Interpreter', 'latex', 'fontsize', 10);
xlabel('n')
ylabel('P(n)')
title('$|\psi_3\rangle \propto |\psi_1\rangle + |\psi_2\rangle$','Interpreter','latex')
end
% ---------------------------------------------


even_cat = CoherentBasis([1;1],[3;-3]);   % create an even cat 
odd_cat  = CoherentBasis([1;-1],[3;-3]);  % create an odd cat

[N_even,P_even] = even_cat.PhotonNumber;  % N_ average photon number
[N_odd,P_odd]   = odd_cat.PhotonNumber;   % P_ photon number distribution 

% ------ plots ---------------------------------
if 1
figure(2)
subplot(1,2,1)
bar(n,P_even,'EdgeColor','none')
text(20,0.25,sprintf('$\\langle n \\rangle $= %d',round(N_even)), ...
    'Interpreter', 'latex', 'fontsize', 10);
xlabel('n')
ylabel('P(n)')
title('$|\psi_{\rm even}\rangle \propto |3\rangle + |-3\rangle$','Interpreter','latex')
subplot(1,2,2)
bar(n,P_odd,'EdgeColor','none')
text(20,0.25,sprintf('$\\langle n \\rangle $= %d',round(N_odd)), ...
    'Interpreter', 'latex', 'fontsize', 10);
xlabel('n')
ylabel('P(n)')
title('$|\psi_{\rm odd}\rangle \propto |3\rangle - |-3\rangle$','Interpreter','latex')
end
% ----------------------------------------------

%%  Wigner function 

% Wigner function of Simple Coherent sates
x_max = 5;                             % needed for the square grid
N     = 600;                           % N points across each direction

W_0 = psi_0.WignerFunction(x_max,N);   % Wigner function W(\alpha, \alpha*) NxN maxtrix 
W_0 = real(W_0);                       % Taking the real part because there will always be a small remenant imaginary component
W_1 = psi_1.WignerFunction(x_max,N);   % coherent state
W_1 = real(W_1);                       % Taking the real part because there will always be a small remenant imaginary component
W_2 = psi_2.WignerFunction(x_max,N);   % superposition of coherent states
W_2 = real(W_2);

W_even_cat = even_cat.WignerFunction(x_max,N);  % Wigner function for even cat
W_even_cat = real(W_even_cat);
W_odd_cat  = odd_cat.WignerFunction(x_max,N);   % Wigner function for odd cat
W_odd_cat = real(W_odd_cat);

compass   = CoherentBasis([1;1;1;1],[3.5;-3.5;3.5i;-3.5i]) ;
W_compass = compass.WignerFunction(x_max,N);   % Wigner function for even cat
W_compass = real(W_compass);

% ---- plots ---------------
cmap_1 = 'turbo';
cmap_2 = 'turbo';
face_alpha = 1;

[X,Y] = meshgrid(linspace(-x_max,x_max,N));


f3 = figure(3);
% --- Plot W_0 ----
if 1 
M_0 = max(max(abs(W_0)));

subplot(1,3,1) 
[~, hc]     =   contourf(X,Y,W_0,100,'EdgeAlpha',0,'FaceAlpha',face_alpha);
clim([-1 1]*M_0);
axis square
a1        =  gca;
a2        =  axes('Parent', f3, 'Position', a1.Position);
hs        =  surf(X, Y, W_0, 'Parent', a2, 'EdgeColor','none');
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
zlim([-1 1]*M_0*1.1)
xlabel('Re(\beta)')
ylabel('Im(\beta)')
zlabel('W_0(\beta)')
colormap(a1,cmap_1)
colormap(a2,cmap_2)
clim([-1 1]*M_0);
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
    Z1 = W_0(i,:);
    plot3(x,Y1,Z1,'Color',[0,0,0,0.3]);    % [0,0,0,Transparency]  0,0,0 = black
end

% Plotting lines in the Y-Z plane

for i = 1:xspacing:length(x)
    X2 = x(i)*ones(size(y)); % a constant vector
    Z2 = W_0(:,i);
    plot3(X2,y,Z2,'Color',[0,0,0,0.3]);
end

hold off

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
end

% --- Plot W_1 ----
if 1
M_1 = max(max(abs(W_1)));

subplot(1,3,2) 
[~, hc]     =   contourf(X,Y,W_1,100,'EdgeAlpha',0,'FaceAlpha',face_alpha);
clim([-1 1]*M_1);
axis square
a1        =  gca;
a2        =  axes('Parent', f3, 'Position', a1.Position);
hs        =  surf(X, Y, W_1, 'Parent', a2, 'EdgeColor','none');
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
zlim([-1 1]*M_1*1.1)
xlabel('Re(\beta)')
ylabel('Im(\beta)')
zlabel('W_1(\beta)')
colormap(a1,cmap_1)
colormap(a2,cmap_2)
clim([-1 1]*M_1);
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
    Z1 = W_1(i,:);
    plot3(x,Y1,Z1,'Color',[0,0,0,0.3]);    % [0,0,0,Transparency]  0,0,0 = black
end

% Plotting lines in the Y-Z plane

for i = 1:xspacing:length(x)
    X2 = x(i)*ones(size(y)); % a constant vector
    Z2 = W_1(:,i);
    plot3(X2,y,Z2,'Color',[0,0,0,0.3]);
end

hold off
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
end

% --- Plot W_2 ----
if 1
M_2 = max(max(abs(W_2)));

subplot(1,3,3) 
[~, hc]     =   contourf(X,Y,W_2,100,'EdgeAlpha',0,'FaceAlpha',face_alpha);
clim([-1 1]*M_2);
axis square
a1        =  gca;
a2        =  axes('Parent', f3, 'Position', a1.Position);
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
zlabel('W_2(\beta)')
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
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 

end

hold off

f4 = figure(4);
% --- Plot W_even_cat ----
if 1
M_2 = max(max(abs(W_even_cat)));

subplot(1,3,1) 
[~, hc]     =   contourf(X,Y,W_even_cat,100,'EdgeAlpha',0,'FaceAlpha',face_alpha);
clim([-1 1]*M_2);
axis square
a1        =  gca;
a2        =  axes('Parent', f4, 'Position', a1.Position);
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
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
end

% --- Plot W_odd_cat ----
if 1
M_2 = max(max(abs(W_odd_cat)));

subplot(1,3,2) 
[~, hc]     =   contourf(X,Y,W_odd_cat,100,'EdgeAlpha',0,'FaceAlpha',face_alpha);
clim([-1 1]*M_2);
axis square
a1        =  gca;
a2        =  axes('Parent', f4, 'Position', a1.Position);
hs        =  surf(X, Y, W_odd_cat, 'Parent', a2, 'EdgeColor','none');
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
zlabel('W_{odd}(\beta)')
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
    Z1 = W_odd_cat(i,:);
    plot3(x,Y1,Z1,'Color',[0,0,0,0.3]);    % [0,0,0,Transparency]  0,0,0 = black
end

% Plotting lines in the Y-Z plane

for i = 1:xspacing:length(x)
    X2 = x(i)*ones(size(y)); % a constant vector
    Z2 = W_odd_cat(:,i);
    plot3(X2,y,Z2,'Color',[0,0,0,0.3]);
end

hold off
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 

end
hold off

% --- Plot W_odd_cat ----
if 1
M_2 = max(max(abs(W_compass)));

subplot(1,3,3) 
[~, hc]     =   contourf(X,Y,W_compass,100,'EdgeAlpha',0,'FaceAlpha',face_alpha);
clim([-1 1]*M_2);
axis square
a1        =  gca;
a2        =  axes('Parent', f4, 'Position', a1.Position);
hs        =  surf(X, Y, W_compass, 'Parent', a2, 'EdgeColor','none');
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
zlabel('W_{odd}(\beta)')
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
    Z1 = W_compass(i,:);
    plot3(x,Y1,Z1,'Color',[0,0,0,0.3]);    % [0,0,0,Transparency]  0,0,0 = black
end

% Plotting lines in the Y-Z plane

for i = 1:xspacing:length(x)
    X2 = x(i)*ones(size(y)); % a constant vector
    Z2 = W_compass(:,i);
    plot3(X2,y,Z2,'Color',[0,0,0,0.3]);
end

hold off
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 

end
hold off

toc
