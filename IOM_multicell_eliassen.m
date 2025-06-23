% last updated 6/23/25
% partially vectorized inbr function to save time
% added colorbar to graph to highlight hub cells

clc; clear; close all;

num_cells = 42; 
% works best when num_cells is a multiple of 3

% randomize initial conditions using Gaussian distribution w/ 10 percent SD
% noise set to 0 for all cells

yinit = zeros(1,num_cells*8);
for i = 1:num_cells
    yinit(8*i - 7) = normrnd(-64.3581, 6.43581); %V
    yinit(8*i - 6) = normrnd(0.0001, 0.00001); %n
    yinit(8*i - 5) = normrnd(460.4832, 46.04832); %cer
    yinit(8*i - 4) = normrnd(0.0959, 0.00959); %c
    yinit(8*i - 3) = normrnd(839.3863, 83.93863); % adp
    yinit(8*i - 2) = normrnd(60.1374,6.01374); %f6p
    yinit(8*i - 1) = normrnd(110.4577,11.04577); %fbp
end 

% alternative idea (TODO?): run 1 simulation of 1 cell, generate random timesteps and
% take initial conditions to be the values as these timesteps

tspan = [0 900000]; % t in msec
dt = 1;

% Ornstein-Uhlenbeck noise parameters
lambda = 0.00001;
sigma = 1.0;

layout = WattsStrogatz(num_cells,3,0.9);
simplify(layout); % removes self-loops and multiple edges between two vertices
disp("Graph generated")

% creates unweighted graph with num_cells vertices, where each cell is, on average,
% connected to 2*(num_cells/3) other vertices. The last input represents
% the probability (beta) that an existing edge is swapped, according to the Watts
% Strogatz algorithm
% beta = 0 -> ring lattice
% beta = 1 -> random graph
% 2/3 to represent that around 65% of cell pairs are coupled 
% (Smolen, Rinzel, & Sherman 1993)

nbrs = adjacency(layout); 


% adjacency matrix for a graph:
% nbrs(i,j) = 1 iff there exists an edge between 
% vertices i and j in the graph, 0 otherwise
% all 0's along main diagonal

num_nbrs = zeros(1,num_cells); % number of neighbors per cell, used in inbr

for i = 1:num_cells
    num_nbrs(i) = length(neighbors(layout,i));
end

toss, hub = max(nbrs); % identify which cell is a "hub", 
% i.e. is connected to the most cells


gc = 175; %gc = 0, gc = 25, gc = 400
taun = 30; % default: 20 
sdpct = 10;

% heterogeneous g, gca, gna
G0 = 15;

sd = sdpct/100*G0;
G_array = G0 + sd*randn(num_cells, 1);
G_array = max(G_array, 1);

gca0 = 900;
sd = sdpct/100*gca0;
gca_array = gca0 + sd*randn(num_cells, 1);

gna0 = 50; %75;
sd = sdpct/100*gna0;
gna_array = gna0 + sd*randn(num_cells, 1);

% TODO: add options for ode15s, ode45
rhs = @(t,y) iom17(t,y,nbrs, gc, taun, G_array, gca_array, gna_array, ...
    lambda, sigma, dt, num_nbrs, num_cells);

tic
disp("Solving ODE with method heuns")
[t, y] = heuns(rhs, tspan, yinit, dt);
disp("Solved!")
toc

addpath('/Users/eethb/Documents/MATLAB');
% plot of the graph, vertex = cell, edge = coupling connection

figure(1);
colormap copper
deg = degree(layout);
nSizes = 2*sqrt(deg-min(deg)+0.2);
nColors = deg;
plot(layout,'MarkerSize',nSizes,'NodeCData',nColors,'EdgeAlpha',0.1)
title('Watts-Strogatz Beta Cell with $N = 120$ nodes, $K = 42$, and $\beta = 0.75$', ...
    'Interpreter','latex')
colorbar

% select 4 cells at random, plot Vm and C graphs
to_graph = randperm(num_cells);

figure(2); % voltage graphs
fig_v = tiledlayout(3,3);
for i = 1:9
    nexttile
    plot(t/1000, y(:,8*to_graph(i)-7))
    title(append('Vm - Cell ',string(to_graph(i))))
end

figure(3); % calcium graphs
fig_c = tiledlayout(3,3);
for j = 1:9
    nexttile
    plot(t/1000, y(:,8*to_graph(j)-4))
    title(append('C - Cell ',string(to_graph(j))))
end

figure(4); % noise current graphs
fig_noise = tiledlayout(3,3);
for k = 1:9
    nexttile
    plot(t/1000, y(:,8*to_graph(k)))
    title(append('Noise - Cell ',string(to_graph(k))))
end

%===================================================================
function [t, yout] = heuns(rhs, tspan, yinit, dt)

num_vars = length(yinit); 
numeq = num_vars;

ybar = NaN(numeq, 1);

y = yinit';

numstep = tspan(2)/dt;
yout = NaN(numstep+1, numeq);  %too big: 838 GB
t = linspace(tspan(1), tspan(2), numstep + 1);
t = t';

yout(1, :) = yinit;

for i = 2:numstep
    dydt1 = rhs(t(i-1), y);
    ybar = y + dt.*dydt1;
    dydt2 = rhs(t(i), ybar);
    y = y + 0.5*dt.*(dydt1 + dydt2);
    yout(i+1, :) = y';
    if mod(i, 100000) == 1
        disp("step " + i + " t = " + num2str(t(i)))
    end
end

end
%===================================================================
function dydt = iom17(t,y,nbrs, gc, taun, G, gca, gna, lambda, sigma, dt, num_nbrs, num_cells)

% Some constants

num_vars = 8;
% y, like yinit, is a row vector with cell1, cell2, ...
% make y a 2D array
y = reshape(y, [num_vars num_cells]); 
y = y';
% make dydt a 2D array
dydt = zeros(num_cells, num_vars);

atot = 3000;
vgk=0.005;
kgk=8.0;
ngk=1.7;
Jgk=vgk.*G.^ngk./(G.^ngk + kgk^ngk);

% taun passed as input
taua=300000;
fca=0.01;
sigmaer=30;
Cm=5300;
kgo=1;
kadp=1;

v = y(:,1);
n = y(:,2);
cer = y(:,3);
c = y(:,4);
adp = y(:,5);
f6p = y(:,6);
fbp = y(:,7);
noise = y(:,8);


atp = 0.5*(atot-adp + sqrt(adp*atot).*sqrt(-2+atot./adp-3*adp/atot));

%###########################
%       Vector field	   #
%###########################

dydt(:,1) = (-(ica(v, gca) + ina(v, n, gna) + ik(v,n) + ikca(c,v) + ikatp(adp,atp,v) ...
         + noise) - inbr(nbrs,num_nbrs, gc,v, num_cells))/Cm; % v

dydt(:,2) = (ninf(v) - n)./taun; % n


dydt(:,3) = fca*sigmaer.*jer(c,cer); % cer


dydt(:,4) = fca*(jmem(c,v,gca) - jer(c,cer)); % c


dydt(:,5) = kadp*(vhydtot(c).*atp-jprod(fbp,c,adp))/taua; % adp
% NOTE: jhyd := vhydtot * atp

dydt(:,6) = kgo*0.3*(Jgk-jpfk(atp,f6p,fbp,adp)); % f6p

dydt(:,7) = kgo*(jpfk(atp,f6p,fbp,adp)-0.5*jpdh(fbp,c)); % fbp

dydt(:,8) = -lambda*noise + sqrt(dt)*sigma*randn(num_cells,1); % noise

% make dydt a column vector
dydt = reshape(dydt', [num_cells*num_vars 1]);

end

% function x = nbrs(num_cells)
% % Creates 5x5x5 cube. Can be modified to create other configurations
% 
% edge = nthroot(num_cells,3);
% face = edge^2;
% x = zeros(num_cells,num_cells);
% 
% for j = 1:num_cells
%     % right
%     if mod(j,edge) ~= 0
%         x(j,j+1) = 1;
%     end
%     % left
%     if mod(j-1,edge) ~= 0
%         x(j,j-1) = 1;
%     end
%     % up
%     if j <= (num_cells - face)
%         x(j,j+face) = 1;
%     end
%     %down
%     if j > face
%         x(j,j-face) = 1;
%     end
%     % forward
%     if mod(j,face) > edge || mod(j,face) == 0
%         x(j,j - edge) = 1;
%     end
%     % backward
%     if mod(j,face) <= face - edge && mod(j,face) ~= 0
%         x(j, j + edge) = 1;
%     end
% 
%     x(j,j) = 0;
% end
% end

% Isoc
% vsoc=-20;
% gsoc=0; %10;
% ksoc=300;
% osoc = ksoc^4/(ksoc^4 + cer^4);
% Isoc = gsoc * osoc * (v-vsoc);

%===================================================================
function x = ica(v, gca)
    vca=25;
    %gca=1000; passed as an array now
    vm=-20; 
    sm=12;
    minf = 1./(1+exp((vm-v)/sm));
    x = gca.*minf.*(v-vca);
end

%===================================================================
function x = ina(v, n, gna)
    vna=25;

    vm=-40;

    sm=9;

    minf = 1./(1+exp((vm-v)/sm));

    minf3 = minf.^3;

    h=0.1-0.5*(n-0.8);

    x = gna.*minf3.*h.*(v-vna);

end
%===================================================================
function x = jpmca(c)
    vpmca=0.042;
    kpmca=0.1;
    x = vpmca*c.^2./(c.^2 + kpmca^2);
end

%===================================================================
function x = jmem(c,v,gca)
    alpha=4.5e-6;
    jin = -alpha*(ica(v,gca));
    
    x = jin - jpmca(c);
end

%===================================================================
function x = jserca(c)
    vserca = 0.125;
    kserca = 0.2;
    x = vserca*c.^2./(c.^2 + kserca^2);
end

%===================================================================
function x = jer(c,cer)
    pleak = 0.000075; %pleak = 0.00015;
    jleak = pleak*(cer-c);

    x = jserca(c) - jleak;
end

%===================================================================
function x = jpdh(fbp,c)
    vpdh=0.005;
    kCaPDH=300;
    kcam=400;

    Jgapdh = sqrt(fbp);
    cam = kcam*c;
    sinfty = cam./(cam + kCaPDH);

    x = vpdh*sinfty.*Jgapdh;
end

%===================================================================
function x = jprod(fbp,c,adp)    
    % basal and glucose-stimulated ADP phosphorylation rate
    %par vg=2
    %par kg=0.05
    %y = vg*(Jpdh/(kg+Jpdh))
    %alt: linearly boost JPDH to r's scale.
    ky=20;
    yprime = ky*jpdh(fbp,c);
    
    % Calcium inhibition of ATP production (Keizer-Magnus effect)
    kkm=0;
    cainh = 1-kkm*c;
    
    r=0.7;
    freezeprod = 0; 
    Jprod0 = 0;
    vphos = exp((r + yprime) .* cainh);
    x = (1 - freezeprod)*vphos.*adp + freezeprod*Jprod0;
end

%===================================================================
function x = vhydtot(c)
    vhyd=5;
    vhydbas=0.8;
    x = vhydbas + vhyd*(jpmca(c) + jserca(c)/2);
end

%===================================================================
function x = ninf(v)
    vn=-16; 
    sn=5; 
    x = 1./(1+exp((vn-v)/sn));
end

%===================================================================
function x = ik(v,n)
    vk = -75;
    gk=2700;
    x = gk*n.*(v-vk);
end

%===================================================================
function x = ikca(c,v)
    vk = -75;
    gkca = 50; 
    %gkca=1500; % for fast
    kd=0.3;
    nkca=4;

    qinf = c.^nkca./(kd^nkca+c.^nkca);
    x = gkca*qinf.*(v-vk);
end

%===================================================================
function x = ikatp(adp,atp,v)
    gkatpbar=27000;
    ktt=1;
    kdd=17; 
    ktd=26;
    vk = -75;

    mgadp = 0.165*adp;
    adp3m = 0.135*adp;
    atp4m = 0.05*atp;
    topo = 0.08*(1+2.*mgadp/kdd) + 0.89.*(mgadp/kdd).^2;
    bottomo = (1+mgadp/kdd).^2 .* (1+adp3m./ktd+atp4m/ktt);
    katpo = topo./bottomo;
    x = gkatpbar*katpo.*(v-vk);
end

%===================================================================
function output = jpfk(atp,f6p,fbp,adp)

% Sent by Patrick 05/31/24
% Gives slightly different output

kpfk=0.05;
vpfk=0.01;
k1=30.000;
k2=1.000;
k3=50000;
k4=1000.000;
famp=0.0200;
fatp=20.00;
ffbp=0.200;
fbt=20.00;
fmt=20.00;

amp = adp.^2./atp;


ampk1 = amp/k1;
fbpk2 = fbp/k2;
f6p2k3 = f6p.^2/k3;
atp2k4 = atp.^2/k4;
wmax = ampk1.*fbpk2.*f6p2k3./(famp*ffbp);
wf6p = f6p2k3.*(1 + atp2k4/fatp + fbpk2/ffbp + fbpk2./fbt.*atp2k4./(ffbp*fatp) + ampk1/famp + ampk1/fmt.*atp2k4/(famp*fatp) + ampk1./fmt.*fbpk2./fbt.*atp2k4./(famp*ffbp*fatp));
wnof6p = 1 + atp2k4 + fbpk2 + fbpk2./fbt.*atp2k4 + ampk1 + ampk1./fmt.*atp2k4 + ampk1.*fbpk2 + ampk1./fmt.*fbpk2./fbt.*atp2k4;
% jpfk:
% output = vpfk.*(wmax + kpfkbas.*wf6p)./( wmax + wf6p + wnof6p);
output = vpfk.*(wmax + kpfk.*wf6p)./( wmax + wf6p + wnof6p);

end

%===================================================================

function inbr = inbr(nbrs,num_nbrs, gc,v, num_cells)
    % Parameter to represent intercellular electrical coupling
    % gc = 0 => completely independent cell activity ("uncoupled")
    % gc -> infinity => unit acts as one 
    % 25 is (for now) considered "low" coupling, 175 is "high" coupling

    % TODO: further vectorize (if possible?)
    inbr = zeros(num_cells,1);
    for i = 1:num_cells
        inbr(i) = (num_nbrs(i) * v(i) - dot(nbrs(i,:), v)); 
    end

    inbr = gc * inbr;
end