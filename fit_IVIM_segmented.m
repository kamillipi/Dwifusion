function calculated_values = fit_IVIM_segmented( ...
    bvals, to_calculation, ...
    D_min, D_max, Dstar_min, Dstar_max, top_signal, bsplit)
%FIT_IVIM_SEGMENTED  Segmented fitting of IVIM parameters [S0, f, D*, D].
%
%   CALCULATED_VALUES = FIT_IVIM_SEGMENTED(BVALS, TO_CALCULATION, ...
%                       D_MIN, D_MAX, DSTAR_MIN, DSTAR_MAX, TOP_SIGNAL, BSPLIT)
%
%   This function estimates IVIM parameters using a segmented two-step fit:
%   (1) diffusion coefficient D and baseline signal S0 from high-b data,
%   followed by (2) pseudo-diffusion coefficient D* and perfusion fraction f
%   from low-b residuals.
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
%   DSTAR_MIN, DSTAR_MAX : scalar
%       Bounds for pseudo-diffusion coefficient D* (mm^2/s).
%
%   TOP_SIGNAL      : scalar
%       Maximum expected baseline signal intensity (used as S0 upper bound).
%
%   BSPLIT          : scalar
%       Threshold b-value (s/mm^2) separating high-b (diffusion) and
%       low-b (perfusion) regimes.
%
%   Output
%   ------
%   CALCULATED_VALUES : [Nsel x 4] matrix
%       Fitted IVIM parameters per voxel, columns correspond to:
%         1. S0      (baseline signal)
%         2. f       (perfusion fraction)
%         3. D*      (pseudo-diffusion coefficient, mm^2/s)
%         4. D       (diffusion coefficient, mm^2/s)
%
%   Notes
%   -----
%   - Step 1: Fit D and S0 from high-b data.
%   - Step 2: Fit D* and f from residual low-b signal.
%   - Parallelized across voxels with `parfor`.
%


    % Separate IVIM components
    ivimb    = bvals(bvals < bsplit);  ivimb    = ivimb(:);
    nonivimb = bvals(bvals > bsplit);  nonivimb = nonivimb(:);

    % Step 1 fit: tissue component (D, S0_tissue)
    fo1 = fitoptions('Method','NonlinearLeastSquares', ...
        'StartPoint',[1e-3 0.6*top_signal], ...
        'Lower', [D_min 0], ...
        'Upper', [D_max top_signal]);
    ft1 = fittype('S0*exp(-x*D)', 'options', fo1);

    % Step 2 fit: perfusion component (D*, S0_perf)
    fo2 = fitoptions('Method','NonlinearLeastSquares', ...
        'StartPoint',[1e-2 0.2*top_signal], ...
        'Lower', [Dstar_min 0], ...
        'Upper', [Dstar_max 0.5*top_signal]);
    ft2 = fittype('S0*exp(-x*Dstar)', 'options', fo2);

    % Progress bar and output allocation
    %progbar = progressBar(size(to_calculation,2), 'pname','Calculating segmented');
    calculated_values = zeros(length(to_calculation), 4);

    % Parallel segmented fitting
    parfor i = 1:length(to_calculation)
        %progbar.progress

        nonivimy      = to_calculation(bvals > bsplit, i);
        tissueandivimy= to_calculation(bvals < bsplit, i);

        % Step 1: fit tissue contribution
        [fit1, ~, ~] = fit(nonivimb, nonivimy, ft1);

        % Residual signal for perfusion
        ivimy = tissueandivimy - fit1.S0 * exp(-fit1.D * ivimb);

        % Step 2: fit perfusion contribution
        [fit2, ~, ~] = fit(ivimb, ivimy, ft2);

        % Store [S0_total, f, D*, D]
        calculated_values(i,:) = [ ...
            (fit1.S0 + fit2.S0), ...
            (fit2.S0 / (fit1.S0 + fit2.S0)), ...
            fit2.Dstar, ...
            fit1.D];
    end
end
