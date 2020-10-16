%% Periodic flow
%
% This code reproduces figures 3 and 4 from the paper [1]. The code
% computes the complex potential for flow past a periodic array of slits.
% The position and lengths of the slits can be changed by changing the 
% centers and radii of the excised circles. The angle of the slits can be
% changed by changing "Chi" and the angle of attack can be changed
% by changing "aoa". You should change the number of grid points so that
% the figure looks nice.
%
% In order to run the code you need to have the SKPrime toolbox [2] in the
% MATLAB search path. You can download the toolbox at 
%
% https://github.com/ACCA-Imperial/SKPrime
%
% References:
%
% [1] "A calculus for flows in periodic domains", P. J. Baddoo & L. J.
% Ayton, Theoretical and Computational Fluid Dynamics.
%
% [2] "The Schottky–Klein prime function: A theoretical and computational 
% tool for applications", ﻿D. G. Crowdy, E. H. Kropf, C. C. Green & M. M. 
% S. Nasser, IMA Journal of Applied Mathematics.

%% Setup the geometry

% Define the centers of the excised discs
DV=[.5;.5;.5].*exp(1i*(.4+2*pi/3*(0:2)'));
% Define the radii of the excised discs
QV=[.3;.3;.3]/2;

% Define the pre-image of infinity
aInfP = .6;
aInfM = -aInfP;

%% Points on the boundary.

% Define grids. Use more grid points to make the figure look nice.
sinGrid = @(np) (1+sin(linspace(-1,1,np)*pi/2))/2;
sinGrid2 = @(np) (sin(linspace(0,1,np)*pi/2));

% Points on the boundary.
thet1 = pi*sinGrid(30).';
thet2 = pi*(1+sinGrid(30)).';
thet = [thet1;thet2];

zetax = [-1+(1-aInfP)*(sinGrid2(50)),-aInfP+2*aInfP*sinGrid(50),aInfP+(1-aInfP)*(1-flip(sinGrid2(50)))];
zetay = [sinGrid(100)-1,sinGrid(100)];

zb1 = exp(1i*thet);

for k = 1:2
%
ZETA = zetax+1i*zetay.';

for M = 0:3

    dv = DV(1:M);
    qv = QV(1:M);
    D = skpDomain(dv,qv);    
    zb = bsxfun(@plus, [0; dv].', bsxfun(@times, [1;qv].', zb1));


for j = 1:numel(dv)
    ZETA(abs(ZETA - dv(j)) < qv(j)+eps(2)) = nan;
end
%
ZETA(abs(ZETA)>1-eps(2))=nan;
ZETA(abs(ZETA)<=aInfP & abs(mod(angle(ZETA),pi))==0) = nan;
zb(abs(mod(angle(zb),pi))==0) = nan;

%% Define prime functions

wa= skprime(aInfP,D); % Prime function with parameter a_{\infty^+}
dwadz= diff(wa); % Derivative of prime function 
wia = invParam(wa); % Prime function with inverse parameter
dwiadz= diff(wia); % Derivative of prime function with inverse parameter
wb = skprime(aInfM,wa); % Prime function with parameter a_{\infty^-}
dwbdz= diff(wb); % Derivative of prime function
wib = invParam(wb); % Prime function with inverse parameter
dwibdz= diff(wib); % Derivative of prime function with inverse parameter

%% Complex Potential

Chi=-pi/4;
period = 1;

% We now define the complex potential for uniform flow.
% For technical reasons relating to the skprime toolbox, we use different
% representations for the simply connected (M=0) and multiply connected
% (M>0) cases. Also, care is taken to define sensible branch cuts so that
% the figures look nice.

if M==0

  potential = @(z,ang) period/(2i*pi)*(...
     exp(-1i*ang)*(log((z-aInfP)./(z-aInfM))) ...
    -exp( 1i*ang)*(log(-(z-1./conj(aInfP))./(z-1./conj(aInfM)))));
    % This differentiates the potential so that we can apply the Kutta
    % condition later.
    uniDeriv = @(z,ang) period/(2i*pi)*(exp(-1i*ang).*(1./(z-aInfP)-1./(z-aInfM))...
                               - exp( 1i*ang).*(1./(z-1/conj(aInfP))-1./(z-1/conj(aInfM))));  
else
potential = @(z,ang) period/(2i*pi)*(...
     exp(-1i*ang)*(log(hat(wa,z)./hat(wb,z))+log((z-aInfP)./(z-aInfM))) ...
    -exp( 1i*ang)*(log(hat(wia,z)./hat(wib,z))+log(-(z-1./conj(aInfP))./(z-1./conj(aInfM)))));

    % This differentiates the potential so that we can apply the Kutta
    % condition later.
    uniDeriv = @(z,ang) period/(2i*pi)*(exp(-1i*ang).*(dwadz(z)./wa(z)-dwbdz(z)./wb(z))...
                               - exp( 1i*ang).*(dwiadz(z)./wia(z)-dwibdz(z)./wib(z)));
end                           

% We can express the conformal map in terms of the uniform flow potential.
f = @(z) exp(1i*Chi).*potential(z,Chi);

% We now define the complex potential due to vortices at the preimages of
% a_{\infty^-}.
cg0 = greensC0(aInfM,D);
v1 = diff(cg0);
if M>0
    cg1 = greensCj(aInfM,1,D);
    v2 = diff(cg1);
else
    cg1 = @(t) 0*t; v2 = cg1;
end
if M>1
    cg2 = greensCj(aInfM,2,D);
    v3 = diff(cg2);
else
    cg2 = @(t) 0*t; v3 = cg2;
end
if M>2
    cg3 = greensCj(aInfM,3,D);
    v4 = diff(cg3);
else
    cg3= @(t) 0*t; v4 = cg3;
end

% Define the position of the other vortices. If you change them then make
% sure that they are inside the flow domain.
vortPosa =.8+.2i;
vortPosb =-.3;

% Now define the complex potential and complex velocities
vortFun1a = greensC0(vortPosa,D);
vortFunDeriv1a = diff(vortFun1a);
vortFun2a = greensC0(aInfM,D);
vortFunDeriv2a = diff(vortFun2a);
vortFun1b = greensC0(vortPosb,D);
vortFunDeriv1b = diff(vortFun1b);
vortFun2b = greensC0(aInfM,D);
vortFunDeriv2b = diff(vortFun2b);

% Define the circulation of the vortices
kappaa = -1; kappab = -1;
vortFun = @(zVar) kappaa*vortFun1a(zVar)...
                + kappab*vortFun1b(zVar);
vortFunDeriv = @(zVar) kappaa*vortFunDeriv1a(zVar)...
                     + kappab*vortFunDeriv1b(zVar);
                 
% Find edges of slits to apply the Kutta condition                 
cornerAngs = zeros(1,M+1);
cornerLocs = zeros(1,M+1);
for j = 1:M+1
    if j == 1
        cornerAngs(j) = fminbnd(@(tVar) imag(f((1 - 1e-5)*exp(1i*tVar))),0,2*pi);
        cornerLocs(j) = (1 - 1e-5)*exp(1i*cornerAngs(1));
    else
        cornerAngs(j) = fminbnd(@(tVar) imag(f(dv(j-1)+(1+0e-15)*qv(j-1)*exp(1i*tVar))),0,2*pi);
        cornerLocs(j) = dv(j-1)+qv(j-1)*exp(1i*cornerAngs(j));
    end
end

%% Define the full complex potential 
aoa = -pi/8; % Angle of attack
uniPotential = @(z) potential(z,aoa) ;
% Set up the matrix and vector for the Kutta condition
A = ([v1(cornerLocs);v2(cornerLocs);v3(cornerLocs);v4(cornerLocs)]);
bVec = (uniDeriv(cornerLocs,aoa) + heaviside(k-1.5)*vortFunDeriv(cornerLocs));
% Find the circulations around each plate such that the Kutta condition is
% satisfied.
circs = -(imag(A.')\imag(bVec.'));

% Calculate the full complex potential
comPot = uniPotential(ZETA) + circs(1)*cg0(ZETA) + circs(2)*cg1(ZETA) + circs(3)*cg2(ZETA) + circs(4)*cg3(ZETA) +heaviside(k-1.5)*vortFun(ZETA);

% Calculate the image of the circular domain
Z = f(ZETA);

%% Plots
ns = 3; % Number of periods to plot.

figure
clf
hold on
vals = linspace(-3,3,400);

for n = -ns:ns
contour(real(Z) + n*period,imag(Z),imag(comPot-period*n*exp(1i*(aoa))),vals, 'color', 'k','LineWidth',1)
plot(f(zb(:,1))+n*period,'b','LineWidth',4)
    if size(zb,2)>1
        plot(f(zb(:,2))+n*period,'r','LineWidth',4);
    end
    if size(zb,2)>2
        plot(f(zb(:,3))+n*period,'g','LineWidth',4);
    end
    if size(zb,2)>3
        plot(f(zb(:,4))+n*period,'m','LineWidth',4);
    end
    if k>1
        plot(f(vortPosa)+n*period,'ko','MarkerFaceColor','k')
        plot(f(vortPosb)+n*period,'ko','MarkerFaceColor','k')
    end
end
axis off
hold off

ylimits= [-.55;.6];
xlimits=[-2*period,3*period];

xlim(xlimits)
ylim(ylimits)

rat = (ylimits(1)-ylimits(2))./(xlimits(1)-xlimits(2));

axis off
box off
axis equal
drawnow
end

end