%% Particles with red noise and sample psi with various scales
clear all, close all
PLOT_NOW=0;
gpuOn = 1
%% physical parameters
% domain size
L=10;
% number of particles
numparts=10000
% time step (defined early since it is needed for the noise memory)
dt=1e-3;
%% The parameters for the red noise as realized by Bartosch's algorithm 
memt=0.2; %sets the memory of the noise.  Need at least ten time steps in memt
rho=exp(-dt/memt); rhoc=sqrt(1-rho*rho); % parameters for Bartosch's method
mysig=0.06; % set the standard deviation of the noise
prevrandx=mysig*randn(numparts,1); % initial random number for x pert.
prevrandy=mysig*randn(numparts,1); % initial random number for y pert.

% Linear shear, U(y)=u0*(y/L), U'=u0/L, psi=u0*y^2/(2*L)
%u0=0.1;
% Parabolic flow
% the domain is [0,L]x[0,L] so want the parabolic flow to be zero at the
% boundaries U(y)=-u0*y*(y-L) umax=U(L/2)=u0*L^2/4;
% so a better way to write things is U(y)=-4*umax/(L^2) * y *(y-L);
umax=0.15;
u0=umax*4/(L*L);
% to get the streamfunction integrate u=psi_y
% psi=-u0*y^3/3+u0*L*y^2/2

% parameters for perturbations
lengthst=(2*pi/L)*4;
ampst=5e-18;
lengthunst=(2*pi/L)*10;
ampunst=0.5;
perunst=10;
%steady part, Linear profile 
%psistbac=@(x,y) u0*y.^2/(2*L);
%steady part, parabolic profile 
psistbac=@(x,y) -u0*y.^3/3+u0*L*y.^2/2;
%steady part, perturbations
psistpert=@(x,y) ampst*cos(lengthst*x-pi/4).*cos(lengthst*y+pi/3);
psist=@(x,y) psistbac(x,y)+psistpert(x,y);
% unsteady part
psiunst=@(x,y,t) ampunst*sin(2*pi*t/perunst).*sin(lengthunst*x).*sin(lengthunst*y);
psi=@(x,y,t) psist(x,y)+psiunst(x,y,t);

dxnum=1e-8*L;twodxnum=2*dxnum;
dynum=dxnum;twodynum=2*dynum;
% time stepping parameters
 numsteps=200; numouts=800;
% ICs
rng('default');

if gpuOn == 0
    x=0.05*L+0.3*L*rand(numparts,1);
    y=0.05*L+0.9*L*rand(numparts,1);
    t=0;
    u=zeros(size(x));v=u;
else
    x=gpuArray(0.05*L+0.3*L*rand(numparts,1));
    y=gpuArray(0.05*L+0.9*L*rand(numparts,1));
    t=0;
    u=gpuArray(zeros(size(x))); v=u;
end

if PLOT_NOW==1
    figure(1)
    clf
    colormap darkjet
    figure(2)
    clf
    colormap hot
end
tic
%% main time step loop
for ii=1:numouts
    ts(ii)=t;
    partsxsave(ii,:)=x;
    partsysave(ii,:)=y;
    partsusave(ii,:)=u;
    partsvsave(ii,:)=v;
    % inner loop
    for jj=1:numsteps   
        % random perturbations
        % generate a new bunch of random numbers
	if gpuOn == 0
		nowrand=randn(numparts,1)*mysig;
	else
        	nowrand=randn(numparts,1, 'gpuArray')*mysig;
	end 
        % use Bartosch's memory to get the right "memory"
        pertnowx=(nowrand*rhoc+rho*prevrandx); prevrandx=pertnowx;
        % repeat for y
	if gpuOn == 0
		nowrand=randn(numparts,1)*mysig;
	else
        	nowrand=randn(numparts,1, 'gpuArray')*mysig;
	end 
        % use Bartosch's memory to get the right "memory"
        pertnowy=(nowrand*rhoc+rho*prevrandy); prevrandy=pertnowy;
        % velocities
        % numerical derivatives for velocity from sreamfunction
        unow=(psi(x,y+dynum,t)-psi(x,y-dynum,t))/twodynum;
        x=x+dt*unow+dt*pertnowx;
        t=t+dt;
        % because x is updated this is symplectic Euler
        vnow=-(psi(x+dxnum,y,t)-psi(x-dxnum,y,t))/twodxnum;
        y=y+dt*vnow+dt*pertnowy;      
    end
    if PLOT_NOW==1
        pcolor(xx,yy,psi(xx,yy,t)),shading flat,hold on, plot(x,y,'ko'),title(t),drawnow
    end
end
% one last save
ts(numouts+1)=t;
partsxsave(numouts+1,:)=x;
partsysave(numouts+1,:)=y;
partsusave(numouts+1,:)=unow;
partsvsave(numouts+1,:)=vnow;
toc
