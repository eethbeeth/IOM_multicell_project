num_cells = 125; 

% randomize initial conditions using Gaussian distribution w/ 10 percent SD
yinit = zeros(num_cells*7, 1);
for i = 1:125
    yinit(7*i - 6) = normrnd(-64.3581, 6.43581);
    yinit(7*i - 5) = normrnd(0.0001, 0.00001);
    yinit(7*i - 4) = normrnd(460.4832, 46.04832);
    yinit(7*i - 3) = normrnd(0.0959, 0.00959);
    yinit(7*i - 2) = normrnd(839.3863, 83.93863);
    yinit(7*i - 1) = normrnd(60.1374,6.01374);
    yinit(7*i) = normrnd(110.4577,11.04577);
end 

% alternative idea: run 1 simulation of 1 cell, generate random numbers and
% take initial conditions based on the timesteps they correspond to 

tspan = [0 900000]; % t in msec

[t, y] = ode45(@(t,y) iom17(t,y), tspan, yinit);

figure;
plot(t/1000, y(:,4))
xlabel('t (sec)')
ylabel('Ca (\muM)')

figure;
plot(t/1000, y(:,221))
xlabel('t (sec)')
ylabel('Ca (\muM)')

figure;
plot(t/1000, y(:,438))
xlabel('t (sec)')
ylabel('Ca (\muM)')

%===================================================================
function dydt = iom17(t,y)

% Some constants

num_cells = 125;
dydt = zeros(7*num_cells,1);

G = 15; %G = 10
atot = 3000;
vgk=0.005;
kgk=8.0;
ngk=1.7;
Jgk=vgk*G^ngk/(G^ngk+kgk^ngk);

taun=20;
taua=300000;
fca=0.01;
sigmaer=30;
Cm=5300;
kgo=1;
kadp=1;

v = zeros(num_cells,1);
n = v;
cer = v;
c = v;
adp = v;
f6p = v;
fbp = v;
atp = v;
neighbors = nbrs(num_cells);

for i = 1:num_cells
v(i) = y(1);
n(i) = y(2);
cer(i) = y(3);
c(i) = y(4);
adp(i) = y(5);
f6p(i) = y(6);
fbp(i) = y(7);
end

for i = 1:num_cells
atp(i) = 0.5*(atot-adp(i) + sqrt(adp(i)*atot)*sqrt(-2+atot/adp(i)-3*adp(i)/atot));
% ratio(i) = atp(i)/adp(i); commented out to speed up computation time

%###########################
%       Vector field	   #
%###########################

dydt(7*i-6) = -(ica(v(i)) + ik(v(i), n(i)) + ikca(c(i),v(i)) + ikatp(adp(i),atp(i),v(i)) + inbr(neighbors, i, v))/Cm; % V
dydt(7*i-5) = (ninf(v(i))-n(i))/taun; % n
dydt(7*i-4) = fca*sigmaer*jer(c(i),cer(i)); % cer
dydt(7*i-3) = fca*(jmem(c(i),v(i)) - jer(c(i),cer(i))); % c
dydt(7*i-2) = kadp*((vhydtot(c(i))*atp(i))-jprod(fbp(i),c(i),adp(i)))/taua; % adp
% NOTE THAT vhydtot*atp = jhyd
dydt(7*i-1) = kgo*0.3*(Jgk-jpfk(atp(i),f6p(i),fbp(i),adp(i))); % f6p
dydt(7*i) = kgo*(jpfk(atp(i),f6p(i),fbp(i),adp(i))-0.5*jpdh(fbp(i),c(i))); % fbp

end
end

function x = nbrs(num_cells)
% Creates 5x5x5 cube. Can be modified to create other configurations

edge = nthroot(num_cells,3);
face = edge^2;
x = zeros(num_cells,num_cells);

for j = 1:num_cells
    % right
    if mod(j,edge) ~= 0
        x(j,j+1) = 1;
    end
    % left
    if mod(j-1,edge) ~= 0
        x(j,j-1) = 1;
    end
    % up
    if j <= (num_cells - face)
        x(j,j+face) = 1;
    end
    %down
    if j > face
        x(j,j-face) = 1;
    end
    % forward
    if mod(j,face) > edge || mod(j,face) == 0
        x(j,j - edge) = 1;
    end
    % backward
    if mod(j,face) <= face - edge && mod(j,face) ~= 0
        x(j, j + edge) = 1;
    end

    x(j,j) = 0;
end
end

% Isoc
% vsoc=-20;
% gsoc=0; %10;
% ksoc=300;
% osoc = ksoc^4/(ksoc^4 + cer^4);
% Isoc = gsoc * osoc * (v-vsoc);

function x = ica(v)
    vca=25;
    gca=1000;
    vm=-20; 
    sm=12;
    minf = 1/(1+exp((vm-v)/sm));
    x = gca*minf*(v-vca);
end

function x = jpmca(c)
    vpmca=0.042;
    kpmca=0.1;
    x = vpmca*c^2/(c^2+kpmca^2);
end

function x = jmem(c,v)
    alpha=4.5e-6;
    jin = -alpha*(ica(v));
    
    x = jin - jpmca(c);
end

function x = jserca(c)
    vserca = 0.125;
    kserca = 0.2;
    x = vserca*c^2/(c^2+kserca^2);
end

function x = jer(c,cer)
    pleak = 0.000075; %pleak = 0.00015;
    jleak = pleak*(cer-c);

    x = jserca(c) - jleak;
end

function x = jpdh(fbp,c)
    vpdh=0.005;
    kCaPDH=300;
    kcam=400;

    Jgapdh = sqrt(fbp);
    cam = kcam*c;
    sinfty = cam/(cam + kCaPDH);

    x = vpdh*sinfty*Jgapdh;
end

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
    freezeprod = 0; Jprod0 = 0;
    vphos = exp((r + yprime) * cainh);
    x = (1 - freezeprod)*vphos*adp + freezeprod*Jprod0;
end

function x = vhydtot(c)
    vhyd=5;
    vhydbas=0.8;
    x = vhydbas + vhyd*(jpmca(c) + jserca(c)/2);
end

function x = ninf(v)
    vn=-16; sn=5; 
    x = 1/(1+exp((vn-v)/sn));
end

function x = ik(v,n)
    vk = -75;
    gk=2700;
    x = gk*n*(v-vk);
end

function x = ikca(c,v)
    vk = -75;
    gkca = 50; 
    %gkca=1500; % for fast
    kd=0.3;
    nkca=4;

    qinf = c^nkca/(kd^nkca+c^nkca);
    x = gkca*qinf*(v-vk);
end

function x = ikatp(adp,atp,v)
    gkatpbar=27000;
    ktt=1;
    kdd=17; ktd=26;
    vk = -75;

    mgadp = 0.165*adp;
    adp3m = 0.135*adp;
    atp4m = 0.05*atp;
    topo = 0.08*(1+2*mgadp/kdd) + 0.89*(mgadp/kdd)^2;
    bottomo = (1+mgadp/kdd)^2 * (1+adp3m/ktd+atp4m/ktt);
    katpo = topo/bottomo;
    x = gkatpbar*katpo*(v-vk);
end

function x = inbr(nbrs,i,v)
    % Parameter to represent intercellular electrical coupling
    % gc = 0 => completely independent cell activity ("uncoupled")
    % gc -> infinity => unit acts as one 
    % 25 is (for now) considered "low" coupling, 175 is "high" coupling

    gc = 0;

    inbr = 0;
    for j = 1:length(nbrs)
        if nbrs(i,j) ~= 0
            inbr = inbr + v(i) - v(j);
        end
    end

    x = gc * inbr;
end
    