% Filippos Tzimkas-Dakis   Virginia Tech  MARCH 2024
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
% CoherentBasis class.
%
% This script produces four frigures.
% Figure 1: population distribution in Fock basis for four different states written in the coherent basis
% Figure 2: population distribution in Fock basis for an even and an odd cat state 
% Figure 3: Wigner functions for (i) vacumm |0> , (ii) coherent |\alpha>,  (iii) superposition of two coherent states |\alpha> + |\beta> 
% Figure 4: Wigner functions for (i) even cat state, (ii) odd cat state , (iii) compass cat state 
%
% The runtime of this script is ~3 seconds on a gaming laptop.
%
% Version V 1.2.2

Functions_path = fullfile(fileparts(mfilename('fullpath')), '..','Functions');
addpath(Functions_path);
Classes_path = fullfile(fileparts(mfilename('fullpath')), '..','Classes');
addpath(Classes_path);
warning('off','MATLAB:nchoosek:LargeCoefficient'); 



%% ----------------------- Section 1 --------------------------------------
% Define your quantum states |psi_1>, |psi_2> using coherent states
close all
clear all
clc
tic

% First we define a simple coherent state of argument \alpha = 1
fprintf('\n\n --------------- Section 1 ----------------------\n')

c1    = 2+2i;                   % coefficient in front of the \ket
s1    = 1;                      % argument inside the \ket 
psi_1 = CoherentBasis(c1,s1);   % |psi_1> = (2+2i)|1>,  here |1> = |α=1> is a coherent state, NOT a Fock state!
fprintf(['\n Before normalization   |ψ_1> = (',num2str(psi_1.Coeff),')|',num2str(psi_1.Kets),'>\n'])
psi_1 = psi_1.normalize;        % normalize our state
psi_1.Coeff;                    % display normalized coefficients
fprintf(['\n After normalization    |ψ_1> = (',num2str(psi_1.Coeff,3),')|',num2str(psi_1.Kets),'>\n'])

% Here we use a function defined in the class to print the state 
fprintf(['\n After normalization    ',psi_1.print_state(1),'\n\n'])


% Now let's define a more general state
c2 = [1;  2];                   % coefficient in front of the \ket
s2 = [-1; 3i];                  % argument inside the \ket 
psi_2 = CoherentBasis(c2,s2);   % |psi_2> = |-1> + 2|3i>, here |-1>, |3i> are coherent states
fprintf(['\n Before normalization   |ψ_2> = ',num2str(psi_2.Coeff(1)),'|',num2str(psi_2.Kets(1)),'> + '...
                                           '',num2str(psi_2.Coeff(2)),'|',num2str(psi_2.Kets(2)),'>\n'])
psi_2 = psi_2.normalize;        % normalize the state
psi_2.Coeff;                    % display coefficients upon normalization
fprintf(['\n After normalization    |ψ_2> = ',num2str(psi_2.Coeff(1),3),'|',num2str(psi_2.Kets(1)),'> + '...
                                           '',num2str(psi_2.Coeff(2),3),'|',num2str(psi_2.Kets(2)),'>\n'])

% Here we use a function defined in the class to print the state 
[~,s2] = psi_2.print_state(2);
fprintf(['\n After normalization    |ψ_2> = ',s2,'\n'])


gg = psi_2.print_state(2);
%% ----------------------- Section 2 --------------------------------------
% Dot product, quantum state addition, Displacement operator

fprintf('\n\n --------------- Section 2 ----------------------\n')
% Dot probucts 
a = braket(psi_1);              % a, b,and c are the same
b = braket(psi_1,psi_1);        % 
c = psi_1.braket();             %
fprintf(['\n braket(psi_1)  =  braket(psi_1,psi_1) =  ',num2str(a),'\n'])

c1 = braket(psi_1,psi_2);       % c1 = <psi_1|psi_2>
c2 = braket(psi_2,psi_1);       % c2 = <psi_2|psi_1> = conj(c1)
fprintf(['\n braket(psi_1,psi_2) = ',num2str(c1,3),'\n'])
fprintf(['\n braket(psi_2,psi_1) = ',num2str(c2,3),'\n\n'])

%% ----------------------- Section 3 --------------------------------------
% Add two quantum states or two objects

fprintf('\n\n --------------- Section 3 ----------------------\n')
psi_3 = psi_1 + psi_2;          % |psi_3> = |psi_1> + |psi_2>    NOT NORMALIZED
psi_3.Coeff;                    % display coefficients
fprintf('\n Before normalization   |ψ_3> = |ψ_1> + |ψ_2>\n')
fprintf(['\n                              = (',num2str(psi_1.Coeff,3),')|',num2str(psi_1.Kets),'> + ',...
                                              num2str(psi_2.Coeff(1),3),'|',num2str(psi_2.Kets(1)),'> + ',...
                                              num2str(psi_2.Coeff(2),3),'|',num2str(psi_2.Kets(2)),'>\n'])
[~,s2] = psi_3.print_state();
fprintf('\n Before normalization   |ψ_3> = |ψ_1> + |ψ_2>\n')
fprintf(['\n                              = ',s2,'\n'])

psi_3 = psi_3.normalize;        % normalize the new state
psi_3.Coeff;                    % display normalized coefficients 
psi_3.Kets;                     % display kets 
fprintf(['\n After normalization    ',psi_3.print_state(3),'\n\n'])
d = braket(psi_3);              % d = <psi_3|psi_3>
fprintf(['\n <ψ_3|ψ_3> = ',num2str(d,3),' = ',num2str(abs(d)),' \n'])
fprintf(['\n <ψ_3|ψ_3> = ',compact_complex(d),' \n\n'])


% Addition of the same state
psi_4 = psi_1 + psi_1;          % |psi_4> = 2*psi_1 = 2*psi_1.Coeff * |psi_1.Kets>       NOT NORMALIZED
fprintf(['\n Before normalization  ',psi_4.print_state(4),'\n\n']);


% scalar multiplication         k*|ψ> = k(c_1|a_1> + c_2|a_2> + ..... c_j|a_j>)
%           ↓
psi_5 = psi_1 + 2*psi_3 + exp(-1i*pi/4)*psi_4;   % |psi_5> = |psi_1> + 2*|psi_3> + exp(-1i*pi/4)*|psi_4>     NOT NORMALIZED
fprintf(['\n Before normalization  ',psi_5.print_state(5),'\n\n']);

%% ----------------------- Section 4 --------------------------------------
% Displacement operator

fprintf('\n\n --------------- Section 4 ----------------------\n')
psi_0 = CoherentBasis(1,0);     % |psi_0> = |0>
psi_0.Kets;                     % display the \kets of this state
[s1,s0] = psi_0.print_state(0);
fprintf(['\n ',s1,'\n'])


z     = 1 + 2i;                 % argument of Displacement operator
psi_z = psi_0.D_(z);            % |psi_z> = D(z)|psi_0> = |z>    
[s1,s2] = psi_z.print_state('z');
fprintf([' |ψ_z> = D(z)|',num2str(psi_0.Kets),'>  = D(',num2str(z),') ',s0,'\n'])
fprintf([' ',psi_z.print_state('z'),'\n'])
psi_z.Kets;                     % display the \kets of this state 

%% ----------------------- Section 5 --------------------------------------
% Annihilation (a) and Creation (a^†) operators

fprintf('\n\n --------------- Section 5 ----------------------\n')
fprintf('\n Annihilation(a) and Creation (a^†) operators, have a look in the script.\n\n')

fprintf(['\n      ',psi_2.print_state(1),'\n\n'])
psi_ = psi_2.A;                  % psi_1.A == a|psi_2> = (c_1 + a_1)|a_1> + (c_2 + a_2)|a_2>  
psi_.Coeff;                      % print the coefficients
fprintf(['\n     a|ψ_1> = ', psi_.print_state('1`'),'\n\n\n'])

psi_0  = CoherentBasis(1+1i,0);   % (1+1i) |0>  (vacuum state in Coherent basis)
n_11 = psi_0.A_dagger; 
fprintf(['\n      ',n_11.print_state('11'),'   This is in Fock basis !!!\n\n'])

psi_ = psi_1.A_dagger;           % psi_1.A_dagger == a^† |psi_2>  see ------>  Phys. Rev. A 43, 492 . 
% ↑  
% psi_   is a state written in the FOCK/NUMBER basis, also known as "Agrawal State" !!!!!!
psi_.Coeff;                      % print the coefficients of Fock states
psi_.n;                          % print kets (Fock basis) with NON-zero coefficients

%% ----------------------- Section 6 -------------------------------------- 
fprintf('\n\n --------------- Section 6 ----------------------\n')
fprintf('\n Photon number plots\n')
% Photon number ------
[N_0,P_n0] = psi_0.PhotonNumber;     % N_ = average photon number
[N_1,P_n1] = psi_1.PhotonNumber;     % P_n = photon distribution as a funtion of n
[N_2,P_n2] = psi_2.PhotonNumber;
[N_3,P_n3] = psi_3.PhotonNumber;

n = 0:length(P_n0)-1;                % for plot reasons
% ---- plots ---------------------------------
f1 = figure(1);
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

f1.Units = 'normalized';
f1.OuterPosition = [0 0.52 0.3367 0.4756]; 
sgtitle('Photon Number Distribution')
% ---------------------------------------------


even_cat = CoherentBasis([1;1],[3;-3]);   % create an even cat 
odd_cat  = CoherentBasis([1;-1],[3;-3]);  % create an odd cat

[N_even,P_even] = even_cat.PhotonNumber;  % N_ average photon number
[N_odd,P_odd]   = odd_cat.PhotonNumber;   % P_ photon number distribution 

% ------ plots ---------------------------------
f2 = figure(2);
subplot(1,2,1)
bar(n,P_even,'EdgeColor','none')
text(20,0.25,sprintf('$\\langle n \\rangle $= %.2f',N_even), ...
    'Interpreter', 'latex', 'fontsize', 10);
xlabel('n')
ylabel('P(n)')
title('$|\psi_{\rm even}\rangle \propto |3\rangle + |-3\rangle$','Interpreter','latex')
subplot(1,2,2)
bar(n,P_odd,'EdgeColor','none')
text(20,0.25,sprintf('$\\langle n \\rangle $= %.2f',N_odd), ...
    'Interpreter', 'latex', 'fontsize', 10);
xlabel('n')
ylabel('P(n)')
title('$|\psi_{\rm odd}\rangle \propto |3\rangle - |-3\rangle$','Interpreter','latex')

f2.Units = 'normalized';
f2.OuterPosition = [0.6621 0.52 0.3367 0.4756]; 
sgtitle('Photon Number Distribution')
% ----------------------------------------------

%% ----------------------- Section 7 --------------------------------------   
% Wigner function 
fprintf('\n\n --------------- Section 7 ----------------------\n')
fprintf('\n Wigner Function plots\n\n')
% Wigner function of Simple Coherent sates
x_max = 5;                             % needed for the square grid
N     = 300;                           % N points across each direction

W_0 = psi_0.WignerFunction(x_max,N);   % Wigner function W(\alpha, \alpha*) NxN maxtrix 
W_0 = real(W_0);                       % Taking the real part because there will always be a small remenant imaginary component
W_1 = psi_1.WignerFunction(x_max,N);   % coherent state
W_1 = real(W_1);                       % Taking the real part because there will always be a small remenant imaginary component
W_2 = psi_2.WignerFunction(x_max,N);   % superposition of coherent states
W_2 = real(W_2);

W_even_cat = even_cat.WignerFunction(x_max,N);  % Wigner function for even cat
W_even_cat = real(W_even_cat);
W_odd_cat  = odd_cat.WignerFunction(x_max,N);   % Wigner function for odd cat
W_odd_cat  = real(W_odd_cat);

compass   = CoherentBasis([1;1;1;1],[3.5;-3.5;3.5i;-3.5i]) ;
W_compass = compass.WignerFunction(x_max,N);   % Wigner function for even cat
W_compass = real(W_compass);

% ---- plots ---------------
cmap_1 = 'turbo';     % for the redblue colormap download https://www.mathworks.com/matlabcentral/fileexchange/25536-red-blue-colormap 
cmap_2 = 'turbo';     % and use 'redblue'
face_alpha = 1;

[X,Y] = meshgrid(linspace(-x_max,x_max,N));


f3 = figure(3);
% --- Plot W_0 ---- 
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
title('W_0')
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 

% --- Plot W_1 ----------
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
title('W_1')
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 


% --- Plot W_2 ----
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
x = X(1,:);
y = Y(:,1);
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
title('W_2')
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 

f3.Units = 'normalized';
f3.OuterPosition = [0    0.05    0.45    0.47]; 

hold off

f4 = figure(4);
% --- Plot W_even_cat ----
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
box off
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
title('Even Cat')
hold off
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 


% --- Plot W_odd_cat ----
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
title('Odd Cat')
hold off

% --- Plot Compass ----------
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
zlabel('W_{compass}(\beta)')
colormap(a1,cmap_1)
colormap(a2,cmap_2)
clim([-1 1]*M_2);
axis square
%  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
x = X(1,:);
y = Y(:,1);
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
    X2 = x(i)*ones(size(y));    % a constant vector
    Z2 = W_compass(:,i);
    plot3(X2,y,Z2,'Color',[0,0,0,0.3]);
end

hold off
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
title('Compass')
hold off
f4.Units = 'normalized';
f4.OuterPosition = [f3.OuterPosition(3)  f3.OuterPosition(2)    1-f3.OuterPosition(3)    0.47]; 
toc
