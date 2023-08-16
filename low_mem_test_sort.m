%% Close distance test, array sorting method

% Theory of operation: the L2 norm (distance) is bounded from above by the
% L∞ (max) norm.  For each particle (i), we can cull most of the distance
% tests by selecting candidates (j) with 
%   xi - thresh ≤ xj ≤ xi + thresh and
%   yi - thresh ≤ yj ≤ yi + thresh

dev = gpuDevice;


%{
numside = 1300
x1d = linspace(0, 1000, numside + 1);
%x1d = gpuArray(linspace(0, 1000, numside + 1));
[xx, yy] = meshgrid(x1d, x1d);
x1 = xx(:);
y1 = yy(:);
%}

gpuOn = 1;

Lx = 10; Ly = 10;
rng('default');
numparts = 5e6;

if gpuOn == 0
    disp('CPU')
    x1 = 0.375*Lx+rand(numparts,1)*Lx*0.25; % Random x location
    y1 = 0.375*Ly+rand(numparts,1)*Ly*0.25; % Random y location
else
    disp('GPU')
    x1 = gpuArray(0.375*Lx+rand(numparts,1)*Lx*0.25); % Random x location
    y1 = gpuArray(0.375*Ly+rand(numparts,1)*Ly*0.25); % Random y location
end

mycutoff = 0.1;

% Selecting these candidates is made much easier if we have a pre-sorted
% array.  Finding the bounds inside a sorted array is fast with a binary
% search.

% Sort the raw coordinates, keeping the indexes around for
% set intersections
tic

[x1_sorted, idx_x1_sorted] = sort(x1);
[y1_sorted, idx_y1_sorted] = sort(y1);

contacts_sorted = zeros(length(x1),1);

mycutoff_sq = mycutoff^2; % Squared distance threshold

% In theory, we can use the sorted arrays to find candidates,
% point-by-point.  However, this involves an outer loop over all points and
% very little vectorization, which results in a slow calculation.  A
% compromise between this and the binning approach is to select points to
% _calculate_ based on bins, but to select _candidates_ based on
% neighbours. 

% Number of bins, per dimension.  Estimate we want something like sqrt(N)
% bins in total to compromise between batch size and a minimal neighbour
% set.  Points sharing a bin will always be each other's neighbours, so if
% points are equally distributed we'd expect O(sqrt(N)) points per bin;
% pdist is O(Npt^2) work, so the overall running time is O(N^3/2).

NBins_sort = max(1,floor(length(x1)^0.25)); %Automatically select the number of bins

% Bin bounds
xbin_sort_bounds = linspace(x1_sorted(1),x1_sorted(end),NBins_sort+1); %set the bounds of each bin
ybin_sort_bounds = linspace(y1_sorted(1),y1_sorted(end),NBins_sort+1);

% Cell array containing particle indices, by 1D bin
xbin_indices = cell(NBins_sort,1);
ybin_indices = cell(NBins_sort,1);

wait(dev);

% Loop over each implied bin, segmenting the points by 1D bin index
% This section we fill each bin cell with the indices that should be in
% each bin.


for ibin = 1:NBins_sort
    % Use the sorted arrays to select points for computation
    xlb = xbin_sort_bounds(ibin);
    xub = xbin_sort_bounds(ibin+1);

    % Search over the x-sorted values to find all indices
    % within x ∈ [xlb, xub] 
    xsel_lb = binary_array_search(x1_sorted,xlb);
    % Debugging note: binary_array_search returns the smallest index
    % greater than _or equal to_ the provided value.  We know that xub
    % is inside x1_sorted because we're choosing real particles, so when
    % xub = max(x1) binary_array_search will return the last index.
    xsel_ub = binary_array_search(x1_sorted,xub);
    xbin_indices{ibin} = sort(idx_x1_sorted(xsel_lb:xsel_ub));
end

wait(dev);

% Repeat for y
for jbin = 1:NBins_sort
    ylb = ybin_sort_bounds(jbin);
    yub = ybin_sort_bounds(jbin+1);
    ysel_lb = binary_array_search(y1_sorted,ylb);
    ysel_ub = binary_array_search(y1_sorted,yub);
    ybin_indices{jbin} = sort(idx_y1_sorted(ysel_lb:ysel_ub));
end

% Now, loop over all bins to write to the contacts array
for ibin = 1:NBins_sort
    for jbin = 1:NBins_sort
        % Find all particles that are in the correct x-bin and y-bin
        selected_indices = intersect(xbin_indices{ibin},ybin_indices{jbin},'sorted');
        if (isempty(selected_indices))
            % No particles, skip
            continue
        end

        % Find the min/max of the selected points to bound the search for
        % neighbours
        xlb = min(x1(selected_indices)) - mycutoff;
        xub = max(x1(selected_indices)) + mycutoff;
        ylb = min(y1(selected_indices)) - mycutoff;
        yub = max(y1(selected_indices)) + mycutoff;

        % Search the sorted x/y coordinates for the indices of all
        % candidate points inside _either_ the x and y limits

        % Debugging note: on the upper bound search, we're only interested
        % in points that are strictly less than max(x1(this)) + mycutoff,
        % so even if we get a real index that satisfies strict equality, we
        % don't need to consider it.
        xsel_lb = binary_array_search(x1_sorted,xlb);
        xsel_ub = binary_array_search(x1_sorted,xub)-1;
        ysel_lb = binary_array_search(y1_sorted,ylb);
        ysel_ub = binary_array_search(y1_sorted,yub)-1;
        

        % Intersect the indices to find the points inside _both_ the x and
        % y limits
        candidate_indices = intersect(idx_x1_sorted(xsel_lb:xsel_ub),...
                                      idx_y1_sorted(ysel_lb:ysel_ub));
        % Rather than use pdist2 to calculate distances, we can use array
        %disp(jbin)
        % broadcasting to calculate squared distance a bit more quickly:
        distances_sq = (x1(selected_indices) - x1(candidate_indices)').^2 + (y1(selected_indices) - y1(candidate_indices)').^2;
        % And count how many are lower than the (squared) threshold,
        % assigning to contacts_sorted
        contacts_sorted(selected_indices) = sum(distances_sq < mycutoff_sq,2);
        %{
        clear distances_sq;
        clear xlb;
        clear xub;
        clear ylb;
        clear yub;
        clear candidate_indices;
        %}
        
    end
end
%{
[max_contacts, max_contacts_idx] = max(contacts_sorted);
x_most=x1(max_contacts_idx);
y_most=y1(max_contacts_idx);
dummy=sqrt((x1-x_most).^2+(y1-y_most).^2);
myneighbours=find(dummy<mycutoff);
%}
wait(dev);

timesort = toc

% fprintf("%d differences between sort and pdist\n",sum(abs(contacts_sorted-contacts0)))


function lbound = binary_array_search(array, val)
    % binary_array_search(array,val)
    % searches ascending-sorted array for value, returning the smallest
    % index that is greater or equal to val. 

    lbound = 1;
    ubound = length(array)+1;

    while (ubound > lbound)
        test = floor(0.5*(lbound + ubound));
        if (array(test) >= val) % Decrease ubound
            ubound = test;
        else % Increase lbound
            lbound = test+1;
        end
    end
end
