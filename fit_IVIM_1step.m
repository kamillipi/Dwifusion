function calculated_values = fit_IVIM_1step(bvals, to_calculation, ...
     D_min, D_max, Dstar_min, Dstar_max, f_min, f_max, ...
    top_signal)

%FIT_IVIM_1STEP  One-step nonlinear fit of IVIM parameters [S0, f, D*, D].
%
%   CALCULATED_VALUES = FIT_IVIM_1STEP(BVALS, TO_CALCULATION, VALUES_LOCATION, ...
%                       D_MIN, D_MAX, DSTAR_MIN, DSTAR_MAX, F_MIN, F_MAX, TOP_SIGNAL)
%
%   This function performs voxel-wise nonlinear least-squares fitting of the
%   IVIM model, estimating signal baseline (S0), perfusion fraction (f),
%   pseudo-diffusion coefficient (D*), and diffusion coefficient (D).
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
%   F_MIN, F_MAX    : scalar
%       Bounds for perfusion fraction f (unitless, 0–1).
%
%   TOP_SIGNAL      : scalar
%       Upper bound for S0 (max signal intensity).
%
%   Output
%   ------
%   CALCULATED_VALUES : [Nsel x 4] matrix
%       Fitted IVIM parameters per voxel, columns correspond to:
%         1. S0
%         2. f
%         3. D*
%         4. D
%
%   Notes
%   -----
%   - Uses MATLAB's nonlinear least-squares `fit` with exponential model:
%         S(b) = S0 * [ f*exp(-b*D*) + (1-f)*exp(-b*D) ]
%   - Parallelized across voxels with `parfor`.
%


    % Fit options and model (D, Dstar, S0, f order in options)
    fo3 = fitoptions('Method','NonlinearLeastSquares', ...
        'StartPoint',[1e-3 1e-4 top_signal 0.1], ...
        'Lower',     [D_min Dstar_min 0        f_min], ...
        'Upper',     [D_max Dstar_max top_signal f_max]);

    ft3 = fittype('S0*f*exp(-x*Dstar)+(1-f)*S0*exp(-x*D)', 'options', fo3);

    % Progress bar (if available) and output allocation
    %progbar = progressBar(size(to_calculation,2), 'pname','Calculating 1step');
    calculated_values = zeros(length(to_calculation), 4);

    % Parallel fit over selected voxels (columns)
    parfor i = 1:length(to_calculation)
        [fit3, ~, ~] = fit(bvals, squeeze(to_calculation(:,i)), ft3);
        calculated_values(i,:) = [fit3.S0, fit3.f, fit3.Dstar, fit3.D];
        %progbar.progress
    end
end