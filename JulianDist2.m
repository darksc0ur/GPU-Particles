
% Define Parameters
L = 5; % Physical Length
cutoff = 0.5; % Cutoff length for a 'close' particle, also the bin size
N = L / cutoff; % Number of Bins
%lngth = 100000000; %100 mil
numparts = 20; % Number of Particles
rng('default')

gpuOn = 1

% Define initial conditions
if gpuOn == 0
    %x1 = [1.25, 3.25, 2.25, 4.75, 4.25] + 0.01;
    %y1 = [1.25, 3.25, 2.25, 4.75, 4.25] + 0.01;

    x1 = unifrnd(0, L, numparts, 1);
    y1 = unifrnd(0, L, numparts, 1);
    
    int_mat = zeros(numparts, numparts);
    int_mat_sh = zeros(numparts, numparts);
else
    %x1 = gpuArray([1.25, 3.25, 2.25, 4.75, 4.25] + 0.01);
    %y1 = gpuArray([1.25, 3.25, 2.25, 4.75, 4.25] + 0.01);

    x1 = gpuArray(unifrnd(0, L, numparts, 1));
    y1 = gpuArray(unifrnd(0, L, numparts, 1));
    
    int_mat = gpuArray(zeros(numparts, numparts));
    int_mat_sh = gpuArray(zeros(numparts, numparts));
end


tic
% Create the bins, arrays return the the bin of each particle (each
% particle corresponds to the index)
closex = ceil(x1*N/L); 
closey = ceil(y1*N/L);

closex_sh = ceil(N*mod(x1+(L/(2*N)), L) / L); %Same as above but for the shift
closey_sh = ceil(N*mod(y1+(L/(2*N)), L) / L);

%{
if gpuOn ~= 0
    closex = gather(closex);
    closey = gather(closey);
    closex_sh = gather(closex_sh);
    closey_sh = gather(closey_sh);
end
%}

Nbins_x = max(closex);
Nbins_y = max(closey);

% Gather particles into cell arrays with cell{i} corresponding to all the
% particles in bin i
for i = 1:Nbins_x
    if gpuOn == 0
        xparts{i} = find(closex == i);
    else
        xparts{i} = gpuArray(find(closex == i));
    end
end

%this is needed to get around empty bins
x_empty = cellfun('isempty', xparts);
xparts(x_empty) = {0};

for i = 1:Nbins_y
    if gpuOn == 0
        yparts{i} = find(closey == i);
    else
        yparts{i} = gpuArray(find(closey == i));
    end
end

y_empty = cellfun('isempty', yparts);
yparts(y_empty) = {0};

store_arr = zeros(2, Nbins_x*Nbins_y);

for i = 1:Nbins_x
    for j = 1:Nbins_y
        if length(intersect(xparts{i}, yparts{j})) > 1
            idx = intersect(xparts{i}, yparts{j});
            int_mat(idx, idx) = 1;
            store_arr(1, i*j) = i;
            store_arr(2, i*j) = j;
        end
    end
end
store_arr(:, any(store_arr == 0)) = [];


% Repeat for shifted
for i = 1:Nbins_x
    if gpuOn == 0
        xparts_sh{i} = find(closex_sh == i);
    else
        xparts_sh{i} = gpuArray(find(closex_sh == i));
    end
end

%this is needed to get around empty bins
x_empty_sh = cellfun('isempty', xparts_sh);
xparts_sh(x_empty_sh) = {0};

for i = 1:Nbins_y
    if gpuOn == 0
        yparts_sh{i} = find(closey_sh == i);
    else
        yparts_sh{i} = gpuArray(find(closey_sh == i));
    end
end

y_empty_sh = cellfun('isempty', yparts_sh);
yparts_sh(y_empty_sh) = {0};


store_arr_sh = zeros(2, Nbins_x*Nbins_y);

for i = 1:Nbins_x
    for j = 1:Nbins_y
        if length(intersect(xparts_sh{i}, yparts_sh{j})) > 1
            idx = intersect(xparts_sh{i}, yparts_sh{j});
            int_mat_sh(idx, idx) = 1;
            store_arr_sh(1, i*j) = i;
            store_arr_sh(2, i*j) = j;
        end
    end
end
store_arr_sh(:, any(store_arr_sh == 0)) = [];
toc

for i = 1:numparts
    int_mat(i, i) = 0;
    it_mat_sh(i, i) = 0;
end



%{
figure(1)
for i = 0:cutoff:L
    for j = 0:cutoff:L
        rectangle('Position', [i, j, cutoff, cutoff]) 
    end
end

for k = 1:length(store_arr(1,:))
    %disp(k)
    rectangle('Position', [store_arr(1,k)*cutoff - cutoff, store_arr(2, k)*cutoff - cutoff, cutoff, cutoff], 'FaceColor', [0, 0, 0])
end

figure(2)
for i = (0:cutoff:L) + (cutoff / 2)
    for j = (0:cutoff:L) + (cutoff / 2)
        rectangle('Position', [i, j, cutoff, cutoff])
    end
end

for k = 1:length(store_arr_sh(1, :))
    rectangle('Position', [store_arr_sh(1, k)*cutoff - (cutoff/2), store_arr_sh(2, k)*cutoff - (cutoff/2), cutoff, cutoff], 'FaceColor', [0, 0, 0])
end
%}