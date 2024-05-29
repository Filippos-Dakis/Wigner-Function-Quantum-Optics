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
% a quantum state written in the Coherent basis. For sake of comparison,
% next to Q-fuction we also plot the Winger function.
% This script produces one figure.
% Figure 1: Wigner and Q fucntion of the EVEN Cat state. 
%
% The runtime of this script is ~1 seconds on a gaming laptop.
%
% Version V 1.2.2

Functions_path = fullfile(fileparts(mfilename('fullpath')), '..','Functions');
addpath(Functions_path);
Classes_path = fullfile(fileparts(mfilename('fullpath')), '..','Classes');
addpath(Classes_path);
warning('off','MATLAB:nchoosek:LargeCoefficient'); 
%%
%% ----------------------- Section 1 --------------------------------------
% Define your quantum states |n_1>, |n_2> usinge number basis

close all
clear all
clc

% First, we define the ground state in the Fock/Number basis
fprintf('\n\n --------------- Section 1 ----------------------\n')

N_hilbert = 20;                 % truncates the Hilbert space up to first 20 states

c1  = 2+2i;                     % coefficient in front of the \ket
n_1 = FockBasis(c1,N_hilbert);  % |n_1> = (2+2i)|0>
fprintf(['\n Before normalization   |n_1> = (',num2str(n_1.Coeff(1)),')|',num2str(n_1.Kets(1)),'>\n'])

% Here we use a function defined in the class to print the state 
fprintf(['\n Before normalization    ',n_1.print_state(1),'\n\n'])


n_1 = n_1.normalize;            % normalize our state
n_1.Coeff;                      % display normalized coefficients 
%                                 n_1.Coeff column vector of length = 20
fprintf(['\n After normalization    |n_1> = (',num2str(n_1.Coeff(1)),')|',num2str(n_1.Kets(1)),'>\n'])

% Here we use a function defined in the class to print the state 
[~,s2] = n_1.print_state(1);
fprintf(['\n After normalization    |n_1> = ',s2,'\n'])

%% ----------------------- Section 2 --------------------------------------
% Now let's define a more general state
fprintf('\n\n --------------- Section 2 ----------------------\n')

c2  = [1; 2; 0; 2i];            % coefficient in front of the \ket
n_2 = FockBasis(c2);            % |n_2> = |0> + 2|1> + 2i|4>
%                                 if you dont put N_hilbert as input, the
%                                 default value is N_hilbert = 30;
fprintf(['\n Before Normalization   |n_2> = ',num2str(n_2.Coeff(1),3),'|',num2str(n_2.n(1)),'> + ',...
                                              num2str(n_2.Coeff(2),3),'|',num2str(n_2.n(2)),'> + ',...
                                              '(',num2str(n_2.Coeff(4),3),')|',num2str(n_2.n(3)),'>\n'])

% Here we use a function defined in the class to print the state 
fprintf(['\n Before normalization    ',n_2.print_state(2),'\n\n'])


n_2 = n_2.normalize;            % normalize the state
n_2.Coeff ;                     % display normalized coefficients
%                                 n_2.Coeff column vector of length = 30
fprintf(['\n After Normalization    |n_2> = ',num2str(n_2.Coeff(1),3),'|',num2str(n_2.Kets(1)),'> + ',...
                                             num2str(n_2.Coeff(2),3),'|',num2str(n_2.Kets(2)),'> + ',...
                                         '(',num2str(n_2.Coeff(4),3),')|',num2str(n_2.Kets(3)),'>\n'])
fprintf(['\n After Normalization    ',n_2.print_state(2),'\n\n'])
%% ----------------------- Section 3 -------------------------------------- 
% Dot product
fprintf('\n\n --------------- Section 3 ----------------------\n')

% Dot probucts 
a = braket(n_1);                % a, b, and c are the same
b = braket(n_1,n_1);            % 
c = n_1.braket();
fprintf(['\n n_1.braket()  =  braket(n_1)  =  braket(n_1,n_1) =  ',num2str(a),'\n'])

d1 = braket(n_1,n_2);           % c1 = <n_1|n_2>
d2 = braket(n_2,n_1);           % c2 = <n_2|n_1> = conj(c1)
fprintf(['\n braket(n_1,n_2) = ',num2str(d1,3),'\n'])
fprintf(['\n braket(n_2,n_1) = ',num2str(d2,3),'\n\n'])

%% ----------------------- Section 4 -------------------------------------- 
fprintf('\n\n --------------- Section 4 ----------------------\n')
% Add two quantum states or two objects
n_3 = n_1 + n_2;                % |n_3> = |n_1> + |n_2>    NOT NORMALIZED
n_3.Coeff;                      % display coefficients
fprintf('\n Before normalization   |n_3> = |n_1> + |n_2>\n')
fprintf(['\n                              = (',num2str(n_3.Coeff(1),3),')|',num2str(n_3.Kets(1)),'> + ',...
                                             num2str(n_3.Coeff(2),3),'|',num2str(n_3.Kets(2)),'> + ',...
                                          '(',num2str(n_3.Coeff(4),3),')|',num2str(n_3.Kets(3)),'>\n'])
fprintf(['\n <n_3|n_3>= ',num2str(braket(n_3),3),'\n'])
n_3 = n_3.normalize;            % normalize the new state
n_3.Coeff;                      % display normalized coefficients 
fprintf(['\n After normalization   |n_3> = (',num2str(n_3.Coeff(1),3),')|',num2str(n_3.Kets(1)),'> + ',...
                                             num2str(n_3.Coeff(2),3),'|',num2str(n_3.Kets(2)),'> + ',...
                                          '(',num2str(n_3.Coeff(4),3),')|',num2str(n_3.Kets(3)),'>\n'])
fprintf(['\n After normalization   ',n_3.print_state(3),'\n'])
d = braket(n_3);                % d = <n_3|n_3>
fprintf(['\n <n_3|n_3>= ',num2str(d,3),'\n\n'])

n_4 = n_1 + n_1;                % |n_4> = 2*n_1 = 2*n_1.Coeff * |n_1>   NOT NORMALIZED
[~,s2] = n_4.print_state(4); 
fprintf(['\n    |n_4> = |n_1> + n_1> = ',s2,'\n'])
% scalar multiplication         k*|ψ> = k(c_1|n_1> + c_2|n_2> + ..... c_j|n_j>)
%           ↓
n_5 = n_1 + 2*n_3 + exp(-1i*pi/4)*n_4;   % |n_5> = |n_1> + 2*|n_3> + exp(-1i*pi/4)*|n_4>     NOT NORMALIZED
fprintf(['\n    ',n_5.print_state(5),' \n'])

%% ----------------------- Section 5 --------------------------------------
fprintf('\n\n --------------- Section 5 ----------------------\n')
% Displacement operator
n_0 = FockBasis(1);             % |n_0> = |0>
n_0.Coeff;                      % display the \kets of this state
z   = 1 + 2i;                   % argument of Displacement operator
n_z = n_0.D_(z);                % |n_z> = D(z)|n_0> = D(z)|0> = |z> 
fprintf('\n The coefficients of the displaced number state are given \n')
n_z.Coeff                       % display the \kets of this state 
%                                n_z.Coeff is an approximate result because
%                                we have truncated the infinite Hilbert space to finite size.

%% ----------------------- Section 6 --------------------------------------
% Annihilation (a) and Creation (a^†) operators
fprintf('\n\n --------------- Section 6 ----------------------\n')
fprintf('\n Annihilation(a) and Creation (a^†) operators, have a look in the script.\n\n')

n_ = n_1.A;                     % n_1.A == a|n_1> = a c_0|0> = c_0*Sqrt(0) = 0
[~,s_] = n_.print_state();
fprintf(['\n    ',n_1.print_state(1),'\n\n'])
fprintf(['\n   a|n_1> = ',s_,' (no state)\n\n\n'])
n_.Coeff;                       % print the coefficients

n_ = n_2.A;                     % n_2.A == a|n_2> = a (c_0|0> + c_1|1> + c_2|2> + c_3|3>)
%                                                 = c_1*Sqrt(1)|0> + c_2*sqrt(2)|1> + c_3*sqrt(3)|2> .                                              
[~,s_] = n_.print_state();
fprintf(['\n    ',n_2.print_state(1),'\n\n'])
fprintf(['\n   a|n_2> = ',s_,'\n\n'])
n_.Coeff;                       % print the coefficients
n_ = n_.normalize();
fprintf(['\n    ',n_.print_state('2`'),'   (upon normalization)\n\n\n'])



n_ = n_1.A_dagger;              % n_1.A_dagger == a^†|n_1> = a^† c_0|0> = c_0*Sqrt(1)|1>
n_.Coeff;                       % print the coefficients
fprintf(['\n        ',n_1.print_state(1),'\n\n'])
[~,s_] = n_.print_state();
fprintf(['\n   (a^†)|n_1> = ',s_,' \n\n'])

n_ = n_2.A_dagger;              % n_2.A_dagger == a^†|n_2> = a^† (c_0|0> + c_1|1> + c_2|2> + c_3|3>)
%                                                 = c_0*Sqrt(1)|1> + c_1*sqrt(2)|2> + c_2*sqrt(3)|3> + c_3*sqrt(4)|4> .                                               
n_.Coeff;                       % print the coefficients
fprintf(['\n        ',n_2.print_state(1),'\n\n'])
[~,s_] = n_.print_state();
fprintf(['\n   (a^†)|n_2> = ',s_,' \n\n'])

%% ----------------------- Section 7 --------------------------------------
fprintf('\n\n --------------- Section 7 ----------------------\n')
fprintf('\n Photon number plots\n')
% Photon number 
[N_0,P_n0] = n_0.PhotonNumber;     % N_ = average photon number
[N_1,P_n1] = n_1.PhotonNumber;     % P_n = photon distribution as a funtion of n
[N_2,P_n2] = n_2.PhotonNumber;
[N_3,P_n3] = n_3.PhotonNumber;
[N_z,P_nz] = n_z.PhotonNumber;


% ---- plots ---------------------------------
f1 = figure(1);
subplot(3,2,1)
bar(0:length(P_n0)-1,P_n0,'EdgeColor','none')
text(5,0.75,sprintf('$\\langle n \\rangle $= %d',N_0), ...
    'Interpreter', 'latex', 'fontsize', 10);
xlabel('n')
ylabel('P(n)')
title('$|\psi_0\rangle = |0\rangle$','Interpreter','latex')
subplot(3,2,2)
bar(0:length(P_n1)-1,P_n1,'EdgeColor','none')
text(5,0.3,sprintf('$\\langle n \\rangle $= %d',N_1), ...
    'Interpreter', 'latex', 'fontsize', 10);
xlabel('n')
ylabel('P(n)')
title('$|\psi_1\rangle = |0\rangle$','Interpreter','latex')
subplot(3,2,3)
bar(0:length(P_n2)-1,P_n2,'EdgeColor','none')
text(20,0.08,sprintf('$\\langle n \\rangle $= %.2f',N_2), ...
    'Interpreter', 'latex', 'fontsize', 10);
xlabel('n')
ylabel('P(n)')
title('$|\psi_2\rangle \propto |0\rangle + 2|1\rangle + 2i|4\rangle$','Interpreter','latex')
subplot(3,2,4)
bar(0:length(P_n3)-1,P_n3,'EdgeColor','none')
text(20,0.3,sprintf('$\\langle n \\rangle $= %.2f',N_3), ...
    'Interpreter', 'latex', 'fontsize', 10);
xlabel('n')
ylabel('P(n)')
title('$|\psi_3\rangle \propto |\psi_1\rangle + |\psi_2\rangle$','Interpreter','latex')
subplot(3,2,5:6)
bar(0:length(P_nz)-1,P_nz,'EdgeColor','none')
text(20,0.1,sprintf('$\\langle n \\rangle $= %.2f',N_z), ...
    'Interpreter', 'latex', 'fontsize', 10);
xlabel('n')
ylabel('P(n)')
title(['$|\psi_z\rangle = \hat{D}(z)|\psi_0\rangle = |',num2str(z),'\rangle$'],'Interpreter','latex')
f1.Units = 'normalized';
f1.OuterPosition = [0 0.4188 0.48 0.58];
% ---------------------------------------------
%sgtitle('Photon Number Distribution')

%% ----------------------- Section 8 --------------------------------------   
%  Wigner function 
fprintf('\n\n --------------- Section 8 ----------------------\n')
fprintf('\n Wigner Function plots\n\n')

x_max = 2.5;                          % needed for the square grid
N     = 100;                          % N points across each direction
%                                       do not increase N too much  because the code slows down.
N_hilbert = 35;                       % increase N_hilbert to get better approximation, however it gets slower !!
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
cmap_1 = 'turbo';     % for the redblue colormap download https://www.mathworks.com/matlabcentral/fileexchange/25536-red-blue-colormap 
cmap_2 = 'turbo';     % and use 'redblue'
face_alpha = 1;

[X,Y] = meshgrid(linspace(-x_max,x_max,N));


f2 = figure(2);
% --- Plot W_0 ----

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
title('W_0')
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 



% --- Plot W_1 -------------------------------------
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
title('W_1')
hold off
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 


% --- Plot W_2 --------------------------------------
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
    X2 = x(i)*ones(size(y));   % a constant vector
    Z2 = W_2(:,i);
    plot3(X2,y,Z2,'Color',[0,0,0,0.3]);
end

hold off
title('W_2')
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
f2.Units = 'normalized';
f2.OuterPosition = [0.4738 0.041 0.5262 0.53];
%sgtitle('Wigner Functions')
toc
