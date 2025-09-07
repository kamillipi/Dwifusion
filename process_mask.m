function [sorted_data, values_loc, orig_size, cat_dim] = process_mask(signal_data, mask, mask_action)
%PROCESS_MASK Apply averaging or calculation based on a mask to the input signal data.
%
%   [sorted_data, values_location, D, ndimensions] = process_mask(signal_data, mask, mask_action)
%
%   INPUTS:
%       signal_data : 2D / 3D / 4D matrix
%                     Last dimension usually corresponds to b-values or time points.
%       mask        : Binary or numerical mask
%                     (0 = excluded, >0 = included).
%       mask_action : String, either "average" or "calculate".
%
%   OUTPUTS:
%       sorted_data     : Reshaped data in 2D [rows x columns].
%       values_loc      : Indices of voxels selected by mask.
%       orig_size       : Spatial dimensions of original image.
%       cat_dim         : Dimension number for concatenating data.
%
%   The function performs one of two actions:
%       - "average": averages data within the mask across spatial dimensions.
%       - "calculate": applies thresholded or provided mask to extract voxel locations.
%
%   Supported input dimensions: 2D, 3D, 4D.

signal_data=double(signal_data);
cat_dim = ndims(signal_data);
if or(ischar(mask),isstring(mask))
    mask=niftiread(mask);
end
%% -------------------------------
%  Case 1: Averaging within mask
%  -------------------------------
if mask_action == "average"
    orig_size=[1 1];
    if isequal(mask, 0)
        % Automatic mask: threshold at 0.1% of max
        mask = max(signal_data, [], cat_dim) > 0.001 * max(signal_data, [], "all");
    elseif isequal(numel(mask), 1) && mask > 1
        mask = max(signal_data, [], cat_dim) > mask;
    end
    avmask = double(mask);
    mask(mask==0)=NaN;

    switch cat_dim
        case 2
            signal_data = squeeze(mean(signal_data .* avmask, [1], "omitnan"));
        case 3
            signal_data = squeeze(mean(signal_data .* avmask, [1 2], "omitnan"));
        case 4
            signal_data = squeeze(mean(signal_data .* avmask, [1 2 3], "omitnan"));
        otherwise
            disp("Please provide proper data: 2D, 3D, or 4D where the last dimension corresponds with bvals")
            sorted_data = -1;
            values_loc = [];
            orig_size = [];
            cat_dim = -1;
            return
    end
    
    % Ensure row vector
    if ~isrow(signal_data)
        signal_data = signal_data';
    end

    dimensions = size(signal_data);
    cat_dim = numel(dimensions);
    values_loc = 1;
else
    dimensions = size(signal_data);
end

%% -------------------------------
%  Reshape data for table format
%  -------------------------------
switch cat_dim
    case 2
        table_cols = dimensions(2);
        table_rows = dimensions(1);
       
    case 3
        table_cols = dimensions(3);
        table_rows = dimensions(1) * dimensions(2);
       
    case 4
        table_cols = dimensions(4);
        table_rows = dimensions(1) * dimensions(2) * dimensions(3);
        
    otherwise
        disp("Please provide proper data: 2D, 3D, or 4D where the last dimension corresponds with bvals")
        sorted_data = -1;
        values_loc = [];
        orig_size = [];
        cat_dim = -1;
        return
end

% Reshape signal_data into [rows x columns] for further calculations
sorted_data = reshape(signal_data, table_rows, table_cols);

%% -------------------------------
%  Case 2: Calculating based on mask
%  -------------------------------
if mask_action == "calculate"
    orig_size=size(signal_data);
    orig_size=orig_size(1:end-1);
    if isequal(mask, 0)
        % Automatic mask: threshold at 0.1% of max
        auto_mask = max(signal_data, [], cat_dim) > 0.001 * max(signal_data, [], "all");
        sorted_mask = reshape(auto_mask, table_rows, 1);
        [values_loc, ~] = find(sorted_mask > 0);
        disp("Using thresholded mask, calculating " + sum(auto_mask,'all')/numel(auto_mask)*100 + "% of data, which is " + sum(auto_mask,'all') + " voxels");

    elseif isequal(numel(mask), 1) && mask > 1
        % Threshold defined by user
        thresh_mask = max(signal_data, [], cat_dim) > mask;
        sorted_mask = reshape(thresh_mask, table_rows, 1);
        [values_loc, ~] = find(sorted_mask > 0);
        disp("Using thresholded mask, calculating " + sum(thresh_mask,'all')/numel(thresh_mask)*100 + "% of data, which is " + sum(thresh_mask,'all') + " voxels");

    else
        % User-provided mask
        sorted_mask = reshape(mask, numel(mask), 1);
        if (table_rows == 1 || table_cols == 1)
            sorted_data = sorted_data';
        else
            sorted_mask = repmat(sorted_mask, table_rows/numel(mask), 1);
        end
        [values_loc, ~] = find(sorted_mask > 0);
        disp("Using provided mask, calculating " + sum(mask,'all')/numel(mask)*100 + "% of data, which is " + sum(mask,'all') + " voxels");
    end
end
end
