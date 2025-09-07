function calculated_values = fit_IVCM_segmented(bvals, to_calculation, ...
    D_min, D_max, top_signal, ...
    v_min, v_max, large_delta, small_delta, bsplit)

%%FIT_IVCM_SEGMENTED  Segmented fitting of IVCM parameters [S0, f, V, D].
%
%   CALCULATED_VALUES = FIT_IVCM_SEGMENTED(BVALS, TO_CALCULATION, ...
%                       D_MIN, D_MAX, TOP_SIGNAL, ...
%                       V_MIN, V_MAX, LARGE_DELTA, SMALL_DELTA, BSPLIT)
%
%   This function estimates intravoxel coherent motion (IVCM) parameters
%   using a segmented two-step fit: (1) diffusion coefficient D and tissue
%   baseline S0 from high-b data, followed by (2) vascular velocity V and
%   perfusion fraction f from low-b residuals, with LARGE_DELTA (Δ) and
%   SMALL_DELTA (δ) as problem parameters in the vascular model.
%
%   Required inputs
%   ---------------
%   BVALS           : [Nb x 1] vector
%       b-values (s/mm^2).
%
%   TO_CALCULATION  : [Nb x Nsel] matrix
%       Diffusion-weighted signal per voxel (one column per voxel).
%
%   D_MIN, D_MAX    : scalar
%       Bounds for diffusion coefficient D (mm^2/s).
%
%   TOP_SIGNAL      : scalar
%       Maximum expected baseline signal intensity (used as S0 upper bound).
%
%   V_MIN, V_MAX    : scalar
%       Bounds for vascular velocity parameter V (mm/s).
%
%   LARGE_DELTA (Δ) : scalar
%       Gradient separation time (s).
%
%   SMALL_DELTA (δ) : scalar
%       Gradient duration (s).
%
%   BSPLIT          : scalar
%       Threshold b-value (s/mm^2) separating high-b (diffusion) and
%       low-b (perfusion/vascular) regimes.
%
%   Output
%   ------
%   CALCULATED_VALUES : [Nsel x 4] matrix
%       Fitted IVCM parameters per voxel, columns correspond to:
%         1. S0      (baseline signal)
%         2. f       (perfusion fraction)
%         3. V       (vascular velocity parameter, mm/s)
%         4. D       (diffusion coefficient, mm^2/s)
%
%   Notes
%   -----
%   - Step 1: Tissue fit of D and S0 from high-b data.
%   - Step 2: Vascular fit of V and f from residual low-b signal,
%             using model exp(-√(b/(LARGE_DELTA - SMALL_DELTA/3)) * V).
%   - Parallelized across voxels with `parfor`.
%



    % Split b-values
    ivimb    = bvals(bvals < bsplit);  
    nonivimb = bvals(bvals > bsplit); 

    % Step 1: fit tissue (D, S0_tissue)
    fo1 = fitoptions('Method','NonlinearLeastSquares', ...
        'StartPoint',[1e-3 0.6*top_signal], ...
        'Lower', [D_min 0], ...
        'Upper', [D_max top_signal]);
    ft1 = fittype('S0*exp(-x*D)', 'options', fo1);

    % Step 2: fit vascular term (V, S0_perf) with problem params (Delta, delta)
    fo2 = fitoptions('Method','NonlinearLeastSquares', ...
        'StartPoint',[0.005*top_signal (v_min+v_max)/2], ...
        'Lower', [0 v_min], ...
        'Upper', [0.25*top_signal v_max]);
    ft2 = fittype('S0*exp(-sqrt(x/(large_delta - small_delta/3))*V)', ...
        'independent','x', 'dependent','y', ...
        'problem', {'large_delta','small_delta'}, ...
        'coefficients', {'S0','V'}, ...
        'options', fo2);

    % Progress + output
    %progbar = progressBar(size(to_calculation,2), 'pname','Calculating segmented');
    calculated_values = zeros(length(to_calculation), 4);

    parfor i = 1:length(to_calculation)
        %progbar.progress

        nonivimy       = to_calculation(bvals > bsplit, i);
        tissueandivimy = to_calculation(bvals < bsplit, i);

        % Step 1: tissue
        [fit1, ~, ~] = fit(nonivimb, nonivimy, ft1);

        % Residual for vascular component
        ivimy = tissueandivimy - fit1.S0 * exp(-fit1.D * ivimb);

        % Step 2: vascular (with problem params)
        [fit2, ~, ~] = fit(ivimb, ivimy, ft2, 'problem', {large_delta, small_delta});

        % Output: [S0_total, f, V, D]
        calculated_values(i,:) = [ ...
            (fit1.S0 + fit2.S0), ...
            (fit2.S0 / (fit1.S0 + fit2.S0)), ...
            fit2.V, ...
            fit1.D ];
    end
end
