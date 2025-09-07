function calculated_values = fit_IVIM_grid(bvals, to_calculation, ...
    D_min, D_max, Dstar_min, Dstar_max, ...
    number_of_points1, number_of_points2, bsplit)

%FIT_IVIM_GRID  Grid search fitting of IVIM parameters [S0, f, D*, D].
%
%   CALCULATED_VALUES = FIT_IVIM_GRID(BVALS, TO_CALCULATION, ...
%                       D_MIN, D_MAX, DSTAR_MIN, DSTAR_MAX, ...
%                       NUMBER_OF_POINTS1, NUMBER_OF_POINTS2, BSPLIT)
%
%   This function estimates IVIM model parameters using a two-step grid
%   search: (1) diffusion coefficient D and baseline signal S0 from
%   high-b data, and (2) pseudo-diffusion coefficient D* and perfusion
%   fraction f from low-b residuals.
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
%   NUMBER_OF_POINTS1 : integer
%       Grid resolution for D and diffusion component S0.
%
%   NUMBER_OF_POINTS2 : integer
%       Grid resolution for D* and perfusion component S0.
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
%   - Step 1: Fit D and S0 by grid search on high-b subset.
%   - Step 2: Fit D* and f from low-b residuals, conditional on signal.
%   - Uses precomputed exponent tables for efficiency.
%   - Parallelized across voxels with `parfor`.
%

    % Split b-values (precompute masks and sizes)
    mask_hi = (bvals > bsplit);
    mask_lo = (bvals < bsplit);
    nonivimb = bvals(mask_hi);
    ivimb    = bvals(mask_lo);
    %n1 = numel(nonivimb);
    %n2 = numel(ivimb);

    % Grids for D and D*
    w2D     = linspace(D_min,    D_max,     number_of_points1);  % D
    w2Dstar = linspace(Dstar_min, Dstar_max, number_of_points2);  % D*

% Precompute exponent tables once (outer products!)
E_hi = exp( - (nonivimb(:)) * (w2D(:).') );        % [n1 x nD]
E_lo = exp( - (ivimb(:))    * (w2Dstar(:).') );    % [n2 x nD*]

    % Also precompute components that don't depend on the signal vector:
    % For each column j in E: s1_j = sum(E(:,j).^2)
    s1_hi = sum(E_hi.^2, 1).';   % [nD x 1]
    s1_lo = sum(E_lo.^2, 1).';   % [nD* x 1]

    % Progress + output allocation
    %progbar = progressBar(size(to_calculation,2),'pname','Calculating grid search');
    calculated_values = zeros(length(to_calculation),4);

    parfor i = 1:length(to_calculation)
        %progbar.progress

        nonivimy = to_calculation(mask_hi, i);
        ivimy    = to_calculation(mask_lo, i);

        % ---------- Step 1: Fit D, S0 from high b-values ----------
        % Build S0 candidate range from data as before
        %max_attenuation = max(nonivimy)/max(ivimy);
        top_signal      = max(ivimy);%/max_attenuation;
        w1 = linspace(0.6*top_signal, top_signal, number_of_points1); % S0 candidates [1 x nS0]

        % Closed-form SSE surface over (D_j, S0_k):
        % SSE(j,k) = sum(y.^2) - 2*S0_k*sum(y.*E(:,j)) + S0_k^2*sum(E(:,j).^2)
        y    = nonivimy(:);                       % [n1 x 1]
        y2   = sum(y.^2);                         % scalar
        s2   = E_hi.' * y;                        % [nD x 1], sum(y.*E(:,j))

        % Build SSE matrix without creating n1-by-(nD*nS0) arrays
        % Use outer products: (s1_hi * w1.^2) and (s2 * w1)
        SSE_hi = y2 - 2*(s2 * w1) + (s1_hi * (w1.^2));  % [nD x nS0]

        % Find the same argmax as original (max likelihood == min SSE)
        SSE_min = min(SSE_hi(:));
        ind_all = find(SSE_hi == SSE_min);
        ind     = round(median(ind_all));    % tie-break to match original
        [xindex,yindex] = ind2sub(size(SSE_hi), ind);

        S0 = w1(yindex);
        Dp = w2D(xindex);

        % ---------- Step 2: Fit D*, f from residual low b-values ----------
        onlyivimy = ivimy - S0*exp(-Dp*ivimb);
        max_ivim_signal = max(onlyivimy);

        if max_ivim_signal > 0.02*S0
            w1b = linspace(0, 1.5*max_ivim_signal, number_of_points2); % S0_ivim candidates

            yb   = onlyivimy(:);             % [n2 x 1]
            yb2  = sum(yb.^2);               % scalar
            s2lo = E_lo.' * yb;              % [nD* x 1]

            SSE_lo = yb2 - 2*(s2lo * w1b) + (s1_lo * (w1b.^2));  % [nD* x nS0b]

            SSE_min2 = min(SSE_lo(:));
            ind_all2 = find(SSE_lo == SSE_min2);
            ind2     = round(median(ind_all2));
            [xindex2,yindex2] = ind2sub(size(SSE_lo), ind2);

            S0iv  = w1b(yindex2);
            Dstar = w2Dstar(xindex2);
            f     = S0iv / (S0 + S0iv);

            calculated_values(i,:) = [ (S0 + S0iv), f, Dstar, Dp ];
        else
            calculated_values(i,:) = [ S0, 0, 0, Dp ];
        end
    end
end
