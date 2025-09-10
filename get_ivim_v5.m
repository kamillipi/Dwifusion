function [results] = get_ivim_v5(bvals_filename, nifti_filename, varargin)

%GET_IVIM_V5  Fit IVIM model parameters from diffusion data.
%
%   RESULTS = GET_IVIM_V5(BVALS_FILENAME, NIFTI_FILENAME, ...)
%
%   This function estimates IntraVoxel Incoherent Motion (IVIM) parameters
%   from diffusion-weighted MRI data. Supply a text file with b-values and a
%   NIfTI file containing the diffusion volumes. Optional name-value pairs
%   configure the fitting method and related settings.
%
%   Required inputs
%   ---------------
%   BVALS_FILENAME   : char/string
%       Path to a text file containing one b-value per volume (same order as
%       in the NIfTI). Units are s/mm^2.
%
%   NIFTI_FILENAME   : char/string
%       Path to a 3D/4D NIfTI file with diffusion data. Supported shapes:
%         - 2D: [voxels, b]   (rare; vectorized)
%         - 3D: [x, y, b]
%         - 4D: [x, y, z, b]
%
%   Name–Value pairs (optional)
%   ---------------------------
%   'method'         : char/string  (default: 'grid')
%       Fitting algorithm for IVIM parameters [D, D*, f, S0].
%         '1step'     - Nonlinear fit of all parameters simultaneously.
%         'segmented'       - Segmented fit: estimate D and S0 using high-b data,
%                       then fit D* and f.
%         'grid'      - Grid search for D and S0, then refine D* and f.
%         'IVCM_grid' - (Variant) Grid approach for IVCM data.
%         'IVCM_seg'  - (Variant) Segmented approach for IVCM data.
%
%   'bsplit'         : positive scalar (default: 250)
%       Threshold b-value (s/mm^2) separating pseudo-diffusion (perfusion)
%       from true diffusion when performing segmented fit.
%
%
%   'large_delta'    : positive scalar (default: 0.043192) representing
%       gradient separation
%   'small_delta'    : positive scalar (default: 0.023848) representing
%       gradient duration
%  Diffusion timing parameters Δ (large_delta) and δ (small_delta) in
%       seconds, used when fitting IVCM (requires gradient timing)
%
%   'mask'           : numeric/logic (default: 0)
%       Either 0 (no mask provided) or a logical/uint8 array matching the
%       spatial dimensions of the data. Non-zero/true entries indicate voxels
%       to process.
%
%   'mask_action'    : "calculate" | "use" (default: "calculate")
%       How to handle masking:
%         "calculate" - Calculate voxel by voxel
%         "average"   - Average signal across all mask
%
%   'save_results'   : "y" | "n" (default: "n")
%       Whether to write output maps to disk alongside the input NIfTI.
%
%   'prefix'         : char/string (default: "IVIM_")
%       Filename prefix for any saved result maps.
%
%   Output
%   ------
%   RESULTS : 4D matrix
%       Container for fitted parameter maps. 
%       Include: 
%        vol 1 S0, 
%        vol 2 f
%        vol 3 Dstar
%        vol 4 D
%
%   Notes
%   -----
%   - Ensure the number of b-values equals the number of diffusion volumes.
%   - Units: b in s/mm^2; diffusion coefficients in mm^2/s.
%
%   Example
%   -------
%     results = get_ivim_v5('sub01.bval', 'sub01_dwi.nii', ...
%                 'method','seg', 'bsplit',200, 'save_results',"y");
%
%

% ---- Input parsing -------------------------------------------------------

p = inputParser;

% Required: paths to b-values file and NIfTI data
addRequired(p, 'bvals_filename');
addRequired(p, 'nifti_filename');

% Allowed fitting methods
allowedmethods = {'1step','segmented','grid','IVCM_grid','IVCM_seg'};

% Validator for strictly positive scalars
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);

% Fitting method (default: 'grid'); validate against allowed set
addParameter(p, 'method', 'grid', @(x) any(validatestring(x, allowedmethods)));

% Segmented fit threshold b (s/mm^2)
addParameter(p, 'bsplit', 250, validScalarPosNum);

% Diffusion timing parameters (seconds)
addParameter(p, 'large_delta', 0.043192, validScalarPosNum);
addParameter(p, 'small_delta', 0.023848, validScalarPosNum);

% Masking options
addParameter(p, 'mask', 0);                 % logical/uint8 array or 0
addParameter(p, 'mask_action', "calculate");% "calculate" | "average"

% Output control
addParameter(p, 'save_results', "n");       % "y" to write maps to disk
addParameter(p, 'prefix', "IVIM_");         % filename prefix for saved maps

% Parse all inputs
parse(p, bvals_filename, nifti_filename, varargin{:});

% Initialize parameters
number_of_points1=int16(600);
number_of_points2=int16(600);
Dstar_min=0.003; D_min=0.000; f_min=0.001; v_min=0;
Dstar_max=0.400; D_max=0.003; f_max=1.000; v_max=3;

%% Process
if or(isstring(bvals_filename),ischar(bvals_filename))
    bvals=load(bvals_filename);                        % Load b-values from file
elseif isvector(bvals_filename)
    bvals=bvals_filename;
else
    error("Invalid bvals")
end

if or(isstring(nifti_filename),ischar(nifti_filename))
    signal_data = niftiread(nifti_filename);                         % Load nifti from file
elseif ismatrix(bvals_filename)
    signal_data = nifti_filename;  
else
    error("Invalid data")
end

if ~iscolumn(bvals)                                % Ensure bvals is a column vector
    bvals=bvals';
end

% Apply/prepare mask, get masked data sorted to calculation
[sorted_data, values_location, orig_size, ndimensions] = process_mask(signal_data, p.Results.mask, p.Results.mask_action);

to_calculation=double(sorted_data(values_location,:));  % Extract voxel signals to calculate

if iscolumn(to_calculation)                           % Ensure proper orientation of the data
else
    to_calculation=to_calculation';
end

calculated_values=zeros(numel(values_location),4);    % Initialize result array [D, D*, f, S0]

top_signal=double(max(signal_data,[],'all'));         % Max signal intensity (for constraints)

disp("Started at " + string(datetime('now')));        % Display start time

switch p.Results.method                               % Select fitting method

    case "1step"
        % Fit all IVIM parameters simultaneously
        calculated_values = fit_IVIM_1step(bvals, to_calculation, ...
            D_min, D_max, Dstar_min, Dstar_max, ...
            f_min, f_max, top_signal);

    case "segmented"
        % Segmented fitting: fit D, S0 first, then D* and f
        calculated_values = fit_IVIM_segmented(bvals, to_calculation, ...
            D_min, D_max, Dstar_min, Dstar_max, ...
            top_signal, p.Results.bsplit);

    case "IVCM_seg"
        % Segmented fitting variant for IVCM model
        calculated_values = fit_IVCM_segmented(bvals, to_calculation, ...
            D_min, D_max, top_signal, ...
            v_min, v_max, p.Results.large_delta, p.Results.small_delta, p.Results.bsplit);

    case "grid"
        % Grid search fitting with noise estimation
        calculated_values = fit_IVIM_grid(bvals, to_calculation, ...
            D_min, D_max, Dstar_min, Dstar_max, ...
            number_of_points1, number_of_points2, p.Results.bsplit);

    case "IVCM_grid"
        % Grid search fitting variant for IVCM model
        calculated_values = fit_IVCM_grid(bvals, to_calculation, ...
            D_min, D_max, v_min, v_max, ...
            number_of_points1, number_of_points2, ...
            p.Results.large_delta, p.Results.small_delta, p.Results.bsplit);

end

imags=zeros(size(sorted_data,1),4);
imags(values_location,:)=calculated_values(1:numel(values_location),:);

S0=reshape(imags(:,1),orig_size);
f=reshape(imags(:,2),orig_size);
Dstar=reshape(imags(:,3),orig_size);
Dstar(f==0)=0;
D=reshape(imags(:,4),orig_size);
results=cat(ndimensions,S0,f,Dstar,D);

if p.Results.save_results=="y"
    save_IVIM(results,nifti_filename,prefix=p.Results.prefix);
end

disp("Ended at " + string(datetime('now')));
end