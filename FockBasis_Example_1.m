% Filippos Tzimkas-Dakis   Virginia Tech  January 2024
%
% Any feedback and suggestions are much appreciated! 
%    
%     ----->  dakisfilippos@vt.edu  <-------
% 
% This script was developed on MATLAB 2023a
%
% Example on FockBasis class
% You should better execute its section one after the other while also reading the comments.
% In this way, you will understand the properties, the methods and the capabilities of the 
% Fock Basis class.
%
% This script produces two frigures. The first figure reveals the population distribution 
% in the Fock basis for each defined state. The second figure plots the Wigner function for 
% three different states defined at the beginning of the last section.
%
% The runtime of this script is ~70 seconds on a gaming laptop. 
% The default Hilbert space used here is  N_hilbert = 30. If you want to reduce the runtime further 
% you must define the states the way I define n_1 in the section below. However, in this way  
% the accuracy/resolution of the Wigner fuction will be reduced. 
% Also, you can change the variable "N_hilbert" at the beginning of the last section, 
% named "Wigner function". I have set it to   N_hilbert = 15   but you change this value and see how
% this affects the runtime and the resolution/accuracy of the Wigner distribution
% 
%
% Version V 1.1

%% Define your quantum states |n_1>, |n_2> usinge number basis
warning('off','MATLAB:nchoosek:LargeCoefficient');     % 

close all
clear all
clc
% First we define a simple coherent state of argument \alpha = 1
c1  = 2+2i;                     % coefficient in front of the \ket
N_hilbert = 20;                 % truncates the Hilbert space up to first 20 states
n_1 = FockBasis(c1,N_hilbert);  % |n_1> = (2+2i)|0>
n_1 = n_1.normalize;            % normalize our state
n_1.Coeff                       % display normalized coefficients 
%                                 n_1.Coeff column vector of length = 20

% Now let's define a more general state
c2 = [1; 2; 0; 2i];             % coefficient in front of the \ket
n_2 = FockBasis(c2);            % |n_2> = |0> + 2|1> + 2i|4>
%                                 if you dont put N_hilbert as input, the
%                                 default value is N_hilbert = 30;
n_2 = n_2.normalize;            % normalize the state
n_2.Coeff                       % display normalized coefficients
%                                 n_2.Coeff column vector of length = 30

%% Dot product, quantum state addition, Displacement operator

% Dot probucts 
a = braket(n_1);               % a and b are the same
b = braket(n_1,n_1);           % 

d1 = braket(n_1,n_2);          % c1 = <n_1|n_2>
d2 = braket(n_2,n_1);          % c2 = <n_2|n_1> = conj(c1)

% Add two quantum states or two objects
n_3 = n_1 + n_2;               % |n_3> = |n_1> + |n_2>    NOT NORMALIZED
n_3.Coeff                      % display coefficients
n_3 = n_3.normalize;           % normalize the new state
n_3.Coeff;                     % display normalized coefficients 
d = braket(n_3);               % d = <n_3|n_3>

n_4 = n_1 + n_1;               % |n_4> = 2*n_1 = 2*n_1.Coeff * |n_1>   NOT NORMALIZED

% Displacement operator
n_0 = FockBasis(1);            % |n_0> = |0>
n_0.Coeff;                     % display the \kets of this state
z  = 1 + 2i;                   % argument of Displacement operator
n_z = n_0.D_(z);               % |n_z> = D(z)|n_0> = D(z)|0> = |z>                 
n_z.Coeff;                     % display the \kets of this state 
%                                n_z.Coeff is an approximate result because
%                                we have truncated the infinite Hilbert space to finite size.

%% Photon number 
[N_0,P_n0] = n_0.PhotonNumber;     % N_ = average photon number
[N_1,P_n1] = n_1.PhotonNumber;     % P_n = photon distribution as a funtion of n
[N_2,P_n2] = n_2.PhotonNumber;
[N_3,P_n3] = n_3.PhotonNumber;


% ---- plots ---------------------------------
if 1 
figure(1)
subplot(2,2,1)
bar(0:length(P_n0)-1,P_n0,'EdgeColor','none')
text(5,0.75,sprintf('$\\langle n \\rangle $= %d',N_0), ...
    'Interpreter', 'latex', 'fontsize', 10);
xlabel('n')
ylabel('P(n)')
title('$|\psi_0\rangle = |0\rangle$','Interpreter','latex')
subplot(2,2,2)
bar(0:length(P_n1)-1,P_n1,'EdgeColor','none')
text(5,0.3,sprintf('$\\langle n \\rangle $= %d',N_1), ...
    'Interpreter', 'latex', 'fontsize', 10);
xlabel('n')
ylabel('P(n)')
title('$|\psi_1\rangle = |0\rangle$','Interpreter','latex')
subplot(2,2,3)
bar(0:length(P_n2)-1,P_n2,'EdgeColor','none')
text(20,0.08,sprintf('$\\langle n \\rangle $= %.2f',N_2), ...
    'Interpreter', 'latex', 'fontsize', 10);
xlabel('n')
ylabel('P(n)')
title('$|\psi_2\rangle \propto |0\rangle + 2|1\rangle + 2i|4\rangle$','Interpreter','latex')
subplot(2,2,4)
bar(0:length(P_n3)-1,P_n3,'EdgeColor','none')
text(20,0.3,sprintf('$\\langle n \\rangle $= %.2f',N_3), ...
    'Interpreter', 'latex', 'fontsize', 10);
xlabel('n')
ylabel('P(n)')
title('$|\psi_3\rangle \propto |\psi_1\rangle + |\psi_2\rangle$','Interpreter','latex')
end
% ---------------------------------------------


%%  Wigner function 

% Wigner function of Simple Coherent sates
x_max = 2.5;                          % needed for the square grid
N     = 100;                          % N points across each direction
%                                       do not increase N too much  because the code slowes down.
N_hilbert = 15;                       % increase N_hilbert to get better approximation, however it gets slower !!
% try N_hilbert = 10 15 20 25 30

n_0 = FockBasis(1,N_hilbert);         % |n_0> = |0> 
n_1 = FockBasis([0;1],N_hilbert);     % |n_1> = |1> 
n_2 = FockBasis([1;1],N_hilbert);     % |n_2> = |0> + |1>

tic
W_0 = n_0.WignerFunction(x_max,N);    % Wigner function W(\alpha, \alpha*) NxN maxtrix 
W_0 = real(W_0);                      % Taking the real part because there will always be a small imaginary component
W_1 = n_1.WignerFunction(x_max,N);    % coherent state
W_1 = real(W_1);
W_2 = n_2.WignerFunction(x_max,N);    % 
W_2 = real(W_2);
toc

% ---- plots ---------------
cmap_1 = 'turbo';
cmap_2 = 'turbo';
face_alpha = 1;

[X,Y] = meshgrid(linspace(-x_max,x_max,N));


f2 = figure(2);
% --- Plot W_0 ----
if 1 
M_0 = max(max(abs(W_0)));

subplot(1,3,1) 
[~, hc]     =   contourf(X,Y,W_0,100,'EdgeAlpha',0,'FaceAlpha',face_alpha);
clim([-1 1]*M_0);
axis square
a1        =  gca;
a2        =  axes('Parent', f2, 'Position', a1.Position);
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
a2        =  axes('Parent', f2, 'Position', a1.Position);
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
a2        =  axes('Parent', f2, 'Position', a1.Position);
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



