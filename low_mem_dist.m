%% A test for a binned version of close distance test
clear all,close all
% physical parameters

rng('default')
Lx=10;Ly=10;
% Nx=512;Ny=512;
% x1d=(0.5:Nx-0.5)*Lx/Nx;
% y1d=(0.5:Ny-0.5)*Ly/Ny;
% [xx,yy]=meshgrid(x1d,y1d);

%{
numside = 1000
x1d = linspace(0, 1000, numside + 1);
[xx, yy] = meshgrid(x1d, x1d);
x1 = xx(:);
y1 = yy(:);
%}


% initial particles
numparts=5e5; % Number of particles
usegpu = 1;

if usegpu == 0
    contacts=zeros(numparts,1); % Number of nearby particles, per particle (init 0)
    x1=0.375*Lx+rand(numparts,1)*Lx*0.25; % Random x location
    y1=0.375*Ly+rand(numparts,1)*Ly*0.25; % Random y location
else
    contacts = gpuArray(zeros(numparts,1)); % Number of nearby particles, per particle (init 0)
    x1 = gpuArray(0.375*Lx+rand(numparts,1)*Lx*0.25); % Random x location
    y1 = gpuArray(0.375*Ly+rand(numparts,1)*Ly*0.25); % Random y location
end


% info for the bins and the cutoff
tic
% Bin particle locations
minx=min(x1); maxx=max(x1); delx=maxx-minx; % Overall size, x
miny=min(y1); maxy=max(y1); dely=maxy-miny; % Overall size, y
mycutoff=0.1; % Threshold for 'close' particle distance
numbins=10; % Number of bins to use
binxi=1+floor((x1-minx)*(numbins-1)/delx); % Bin index per particle, x
binyi=1+floor((y1-miny)*(numbins-1)/dely); % Bin index per particle, y
% loop over bin numbers
for xi=1:numbins % Loop over x-bin
    % needed later for the frame around current bin
    xil=xi-1+numbins*(xi==1); % The number of the x-left bin
    %xir=mod(xi+1,numbins);    % Number of the x-right bin
    % Corrected formula.  Bins are 1-indexed, so if we have 10 bins then
    % the bin right of 9 is 10, and the bin right of 10 is 1.
    xir=1+mod(xi,numbins);    % Number of the x-right bin
    % are you in the right x bin number
    dummyxbin=find(binxi==xi); % Select all particles in this x-bin
    for yi=1:numbins % Loop over y-bin
        % are you in the right ybin number
        dummyybin=find(binyi==yi); % Select all particles in this y-bin
        % make sure you are in both
        dummybin=intersect(dummyxbin,dummyybin); % Select all particles in both bins
        % initialize xnow and ynow
        xnow=[];
        ynow=[];
        % build the frame around the current bin
        yil=yi-1+numbins*(yi==1); % Number of the y-down bin
        %yir=mod(yi+1,numbins); % Number of the y-up bin
        % Corrected formula, see computation of xir
        yir = 1+ mod(yi,numbins); % Number of the y-up bin
        xis=[xil xi xir]; 
        yis=[yil yi yir];
        for xii=1:3 % Loop through left, here, right bins
            dummyx=find(binxi==xis(xii));
            for yii=1:3 % Loop through down, here, up bins
                dummyy=find(binyi==yis(yii));
                dummy=intersect(dummyx,dummyy);
                xbin=x1(dummybin);ybin=y1(dummybin);
                xnow=[xnow ; x1(dummy)];
                ynow=[ynow ; y1(dummy)];
            end
        end
        % calculate the distance from particles in current bin
        % and the frame
        mydists=pdist2([xbin ybin],[xnow ynow]);
        % update the contacts
        myclosebin=sum((mydists<mycutoff),2);
        contacts(dummybin)=myclosebin;
    end
end
bintime=toc

%{
%% now do it again with pdist
tic
mydist0=squareform(pdist([x1 y1]));
time0=toc
contacts0=sum(mydist0<mycutoff,2);
%% Test for differences
fprintf("%d differences between bin and pdist\n",sum(abs(contacts-contacts0)))

%}
