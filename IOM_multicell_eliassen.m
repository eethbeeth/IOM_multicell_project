% by: Ethan Eliassen
% last updated 8/11/25
% un-looped inbr, runs much faster (approx. 5 mins for 125 cells)
% TODO: detect communities within graph?

clc; clear; close all;

num_cells = 125; % num_cells = 25

% randomize initial conditions using Gaussian distribution w/ 10 percent SD
% noise IC set to 0 for all cells

yinit = zeros(1,num_cells*8);
% for i = 1:num_cells
%     yinit(8*i - 7) = normrnd(-64.3581, 6.43581); %V
%     yinit(8*i - 6) = normrnd(0.0001, 0.00001); %n
%     yinit(8*i - 5) = normrnd(460.4832, 46.04832); %cer
%     yinit(8*i - 4) = normrnd(0.0959, 0.00959); %c
%     yinit(8*i - 3) = 839.3863; % adp
%     yinit(8*i - 2) = normrnd(60.1374,6.01374); %f6p
%     yinit(8*i - 1) = 110.4577; %fbp
% end 

for i = 1:num_cells
    yinit(8*i - 7) = -64.3581; %V
    yinit(8*i - 6) = 0.0001; %n
    yinit(8*i - 5) = 460.4832; %cer
    yinit(8*i - 4) = 0.0959; %c
    yinit(8*i - 3) = 839.3863; % adp
    yinit(8*i - 2) = 60.1374; %f6p
    yinit(8*i - 1) = 110.4577; %fbp
end 

% alternative idea (TODO?): run 1 simulation of 1 cell, generate random timesteps and
% take initial conditions to be the values as these timesteps

tspan = [0 900000]; % t in msec
dt = 1;

% Ornstein-Uhlenbeck noise parameters
lambda = 0.00001;
sigma = 1.0;

% creates unweighted graph with num_cells vertices, where each cell is, on average,
% connected to 6 other vertices. The last input represents
% the probability (beta) that an existing edge is swapped, according to the Watts
% Strogatz algorithm
% beta = 0 -> ring lattice
% beta = 1 -> random graph

% layout = WattsStrogatz(num_cells,3,0.9);

% creates unweighted graph according to the Barabási-Albert model
% starts (for simplicity) with a base graph of 3 interconnected nodes
% for every node added, it is connected to 3 existing nodes
% the probability that a new node connects to existing node i is 
% (number of edges incident to i)/(total number of edges in graph)
% over time more well-connectedd nodes tend to accumulate more connections
% ("rich-get-richer" phenomenon)
% stops adding once num_cells nodes have been added
% for simplicity a new node cannot connect to itself or to one cell
% multiple times

layout = BarabasiAlbert(3,num_cells);
simplify(layout); % removes self-loops and multiple edges between two vertices
disp("graph generated")

nbrs = adjacency(layout); 

% adjacency matrix for a graph:
% nbrs(i,j) = 1 iff there exists an edge between 
% vertices/nodes i and j in the graph, 0 otherwise
% all 0's along main diagonal (since no self-loops)

num_nbrs = zeros(1,num_cells); % number of neighbors per cell, used in inbr

for i = 1:num_cells
    num_nbrs(i) = length(neighbors(layout,i));
end

% identify which cells are the main "hubs", 
% i.e. are coupled to the most other cells
% after running many simulations (W-S), 2-4 cells usually have >=10 connections

[~, hub] = maxk(num_nbrs, 4); 
disp(hub + ": hub cell")


gc = 25.*nbrs; %gc = 0, gc = 25, gc = 175, gc = 400
% islet very hard to excite when coupling strength is high...
taun = 30; % default: 20 
sdpct = 10;

gc(hub,:) = 0;
gc(:,hub) = 0;
% SILENCES (hub) cells

% heterogeneous g, gca, gna
G0 = 20;

sd = sdpct/100*G0;
G_array = G0 + sd*randn(num_cells, 1);
G_array = max(G_array, 1);

gca0 = 900;
sd = sdpct/100*gca0;
gca_array = gca0 + sd*randn(num_cells, 1);

gna0 = 50; %75;
sd = sdpct/100*gna0;
gna_array = gna0 + sd*randn(num_cells, 1);

% would a more precise solver be helpful/efficient?
rhs = @(t,y) iom17(t,y, gc, taun, G_array, gca_array, gna_array, ...
    lambda, sigma, dt, num_cells, hub);

tic
[t, y] = heuns(rhs, tspan, yinit, dt);
disp("Solved!")
toc

addpath('/Users/eethb/Documents/MATLAB');

% plot of the graph, vertex = cell, edge = coupling connection
figure;
colormap copper
deg = degree(layout);
nSizes = 2*sqrt(deg-min(deg)+0.2);
nColors = deg;
plot(layout,'MarkerSize',nSizes,'NodeCData',nColors,'EdgeAlpha',0.1, 'Nodelabel', 1:num_cells)
title('Barabasi-Albert Islet with $N = 125$ nodes, $K = 3$, and $\beta = 0.9$', ...
    'Interpreter','latex')
colorbar

% select 6 cells at random, plot Vm and C graphs
to_graph = randperm(num_cells);

figure; % voltage graphs
fig_v = tiledlayout(3,3);
for i = 1:6
    nexttile
    plot(t/1000, y(:,8*to_graph(i)-7))
    title(append('Vm - Cell ',string(to_graph(i))))
end

for j = 1:3 % plot hub cell voltage
    nexttile
    plot(t/1000,y(:,8*hub(j)-7))
    title(append('Vm - Cell ',string(hub(j))))
end

figure; % calcium graphs
fig_c = tiledlayout(3,3);
for i = 1:6
    nexttile
    plot(t/1000, y(:,8*to_graph(i)-4))
    title(append('C - Cell ',string(to_graph(i))))
end

for j = 1:3 % plot hub cell calcium
    nexttile
    plot(t/1000,y(:,8*hub(j)-4))
    title(append('C - Cell ',string(hub(j))))
end

figure; % noise current graphs
fig_noise = tiledlayout(3,3);
for k = 1:6
    nexttile
    plot(t/1000, y(:,8*to_graph(k)))
    title(append('Noise - Cell ',string(to_graph(k))))
end

for j = 1:3 % plot hub cell noise
    nexttile
    plot(t/1000,y(:,8*hub(j)))
    title(append('Noise - Cell ',string(hub(j))))
end

threshold = -45; % threshold that determines when a cell is idle/bursting, could be lower?
% computes times at which the cells burst and the time since the first cell bursts
[crossings, diff] = v_threshold(y,num_cells,threshold,gc,hub); 
vbar = compute_vbar(y,num_cells); % computed average voltage of entire islet at each timestep
figure; plot(t/1000,vbar); % plot vbar

% create movie of bursting progression (black = not bursting, red =
% bursting)
M = plot_bursting(layout, y, num_cells, threshold, crossings, hub); 

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
function dydt = iom17(t,y, gc, taun, G, gca, gna, lambda, sigma, dt, num_cells, hub)

% Some constants

num_vars = 8;
% y, like yinit, is a row vector with cell1, cell2, ...
% make y a 2D array
y = reshape(y, [num_vars num_cells]); 
y = y';
% make dydt a 2D array
dydt = zeros(num_cells, num_vars);

atot = 3000;
% hub cells observed increase glucokinase rxn rate
% TODO: figure out + implement different (better?) 
% ways to hyperpolarize/silence cells
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

% heterogeneous gc

% if t > 500000
%     gc = 50;
% end

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

dydt(:,1) = -(ica(v, gca) + ina(v, n, gna) + ik(v,n) + ikca(c,v) + ikatp(adp,atp,v) ...
         + noise - inbr(gc,v))/Cm; % v

% dydt(hub,1) = 0; % HYPERPOLARIZE (hub) cells

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
function inbr = inbr(gc,v)
    % Parameter to represent intercellular electrical coupling
    % gc = 0 => completely independent cell activity ("uncoupled")
    % gc -> infinity => unit acts as one 
    % 25 is (for now) considered "low" coupling, 175 is "high" coupling
   
    inbr = gc*v  - v .* sum(gc,2);
end
%===================================================================
% computes the average value of v across all cells for each timestep
function vbar = compute_vbar(y, num_cells)
    v_only = var_extract(y,num_cells, "v");
    vbar = NaN(length(y),1);

    for j = 1:length(y)
        vbar(j) = mean(v_only(j,:));
    end
end

%===================================================================
% creates an array of the values of one particular variable across time
function one_var = var_extract(y, num_cells, var)
    one_var = NaN(length(y), num_cells);
    switch(var)
        case "v"
            for i = 1:num_cells
                one_var(:,i) = y(:,8*i - 7);
            end
        case "n"
            for i = 1:num_cells
                one_var(:,i) = y(:,8*i - 6);
            end
        case "cer"
            for i = 1:num_cells
                one_var(:,i) = y(:,8*i - 5);
            end
        case "c"
            for i = 1:num_cells
                one_var(:,i) = y(:,8*i - 4);
            end
        case "adp"
            for i = 1:num_cells
                one_var(:,i) = y(:,8*i - 3);
            end
        case "f6p"
            for i = 1:num_cells
                one_var(:,i) = y(:,8*i - 2);
            end
        case "fbp"
            for i = 1:num_cells
                one_var(:,i) = y(:,8*i - 1);
            end
        case "noise"
            for i = 1:num_cells
                one_var(:,i) = y(:,8*i);
            end
        otherwise
            disp("invalid variable")
    end 
end

%===================================================================
function [crossings,diff] = v_threshold(y,num_cells, dta,gc,hub)
    v = var_extract(y,num_cells,"v");
    vbar = compute_vbar(y,num_cells);
    vbar_crosstime = find(vbar > dta, 1, "first");
    crossings = NaN(1, num_cells);
    for i = 1:num_cells
        if v(vbar_crosstime, i) >= dta % cell is bursting at this time
            % search backwards from this time until the cell begins its
            % burst
            cross_time = find(v(1:vbar_crosstime,i) > dta, 1, "last");

            if isempty(cross_time)
                cross_time = NaN;
            else 
                cross_time = cross_time + 1;
            end
        
        else % cell is not bursting, look forward until it bursts
            cross_time = find(v(vbar_crosstime:end,i) > dta, 1, "first");
            if isempty(cross_time)
                cross_time = NaN;
            else
                cross_time = cross_time + vbar_crosstime;
            end
        end
       
        crossings(i) = cross_time;
    end
    
    % check if any cells in the islet burst
    if all(isnan(crossings))
        disp("This islet did not burst")
        diff = crossings;
        return;
    end
    
    % identify cells that did not burst
    if anynan(crossings)
        for j = 1:length(crossings)
            if isnan(crossings(j))
                disp("Cell " + j + " did not burst")
            end
        end
    end

    if gc(hub,:) == 0
        crossings(hub) = NaN; % disregard hub cells that have been shut off
        disp("Hub cells have been SILENCED")
    end
    
    % identify first responder
    [min_time, first_responder] = min(crossings);
    disp("the first responder is cell " + first_responder + ", bursts at " + min_time/1000 + " seconds")

    % comute average time for a given cell to burst AFTER the first
    % responder

    diff = (crossings - min_time)/1000;
    avg_wait = mean(diff, 'omitnan');
    sd_wait = std(diff, 'omitnan');
    disp("the average time for a cell to respond is " + avg_wait + " seconds")
    disp("the standard deviation is " + sd_wait);
    
    figure; 
    histogram(diff(~isnan(diff)),'BinWidth', 0.02);
    
    % identify last responder
    [max_time, last_responder] = max(crossings);
    disp("the last responder is cell " + last_responder);
    
    % calculate time it takes for the cell to go from completely idle to bursting 
    duration = (max_time - min_time)/1000;
    disp("the time it takes for the islet to go from idle to bursting is " + duration + " seconds")
 
end
%===================================================================
% Copyright 2015 The MathWorks, Inc.
function h = WattsStrogatz(N,K,beta)

% H = WattsStrogatz(N,K,beta) returns a Watts-Strogatz model graph with N
% nodes, N*K edges, mean node degree 2*K, and rewiring probability beta.
%
% beta = 0 is a ring lattice, and beta = 1 is a random graph.

% Connect each node to its K next and previous neighbors. This constructs
% indices for a ring lattice.
s = repelem((1:N)',1,K);
t = s + repmat(1:K,N,1);
t = mod(t-1,N)+1;

% Rewire the target node of each edge with probability beta
for source=1:N    
    switchEdge = rand(K, 1) < beta;
    
    newTargets = rand(N, 1);
    newTargets(source) = 0;
    newTargets(s(t==source)) = 0;
    newTargets(t(source, ~switchEdge)) = 0;
    
    [~, ind] = sort(newTargets, 'descend');
    t(source, switchEdge) = ind(1:nnz(switchEdge));
end

h = graph(s,t);
end
%===================================================================
function M = plot_bursting(G, y, num_cells, threshold, crossings, hub)
% TODO: ensure a movie can be made even if not all cells burst
    if all(isnan(crossings))
        disp("No movie made")
        return
    end

    [min_time, ~] = min(crossings);
    [max_time, ~] = max(crossings);
    fr = 30;
    elapsed_time = max_time - min_time;
   
    if elapsed_time <= 1000
        fr = 5;
    elseif elapsed_time > 1000 && elapsed_time <= 100000
        fr = 10;
    end

    v = var_extract(y, num_cells, "v");
    nColors = zeros(num_cells,3); % RGB data
    nColors(hub,2) = 1; % highlight hubs
    
    h = figure;  
    deg = degree(G);
    nSizes = 2*sqrt(deg-min(deg)+0.2);
    
    plot(G,'MarkerSize',nSizes,'NodeColor',nColors,'EdgeAlpha',0.1, 'Nodelabel', 1:num_cells);
    axis tight manual
    ax = gca;
    ax.NextPlot = 'replaceChildren';

    loops = length((min_time-100):10:(max_time+100));
    M(loops) = struct('cdata',[],'colormap', []);
    h.Visible = 'off';
    
    cd  '/Users/eethb/Documents/MATLAB'
    vw = VideoWriter('burst_progression.avi');
    vw.FrameRate = fr;
    open(vw);
    set(gcf, 'renderer', 'zbuffer');

    for j = (min_time-100):10:(max_time+100)
        idx = 1;
        for cell = 1:num_cells
            if v(j,cell) > threshold
                nColors(cell, 1) = 1; % red
                nColors(hub,3) = 1; % ensures "quiet" hubs are a different color
            end
        end
    
    plot(G,'MarkerSize',nSizes,'NodeColor',nColors,'EdgeAlpha',0.1, 'Nodelabel', 1:num_cells);
    drawnow
    M(idx) = getframe(gcf);
    writeVideo(vw,M(idx));
    idx = idx + 1;
    end

    h.Visible = 'on';
    movie(h,M,1,vw.FrameRate);
    close(vw);
end
%===================================================================
function G = BarabasiAlbert(m, num_cells)
    % start with 3 interconnected nodes ("triangle graph")
    G = graph([1 2 3], [2 3 1]); 
    connections = zeros(num_cells,1); % number of edges per node
    connections(1:3) = 2; % 2 edges per node to start
    for i = 4:num_cells
        total_conns = sum(connections);
        new_conns = randsample(num_cells,m,true, connections./total_conns);
        % add node i to G with connections to nodes in new_conns
        G = addedge(G, i*ones(m,1), new_conns); 
        connections(i) = m;
        connections(new_conns) = connections(new_conns) + 1;
    end
end