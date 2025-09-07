function calculated_values = fit_IVCM_grid(bvals, to_calculation, ...
    D_min, D_max, v_min, v_max, ...
    number_of_points1, number_of_points2, ...
    large_delta, small_delta, bsplit)

%FIT_IVCM_GRID  Grid search fitting of IVCM parameters [S0, f, V, D].
%
%   CALCULATED_VALUES = FIT_IVCM_GRID(BVALS, TO_CALCULATION, ...
%                       D_MIN, D_MAX, V_MIN, V_MAX, ...
%                       NUMBER_OF_POINTS1, NUMBER_OF_POINTS2, ...
%                       LARGE_DELTA, SMALL_DELTA, BSPLIT)
%
%   This function extends the IVIM grid search to an intravoxel coherent
%   motion (IVCM) model, estimating diffusion coefficient (D), vascular
%   flow velocity (V), perfusion fraction (f), and baseline signal (S0).
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
%   V_MIN, V_MAX    : scalar
%       Bounds for vascular velocity parameter V (mm/s).
%
%   NUMBER_OF_POINTS1 : integer
%       Grid resolution for D and S0.
%
%   NUMBER_OF_POINTS2 : integer
%       Grid resolution for V and perfusion component.
%
%   LARGE_DELTA     : scalar
%       Gradient separation time Δ (s).
%
%   SMALL_DELTA     : scalar
%       Gradient duration δ (s).
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
%   - Step 1: Fit D and S0 by grid search on high-b subset.
%   - Step 2: Fit V and f from low-b residuals using Δ and δ.
%   - Vascular term uses model: exp(-Δ * sqrt(b/(Δ - δ/3)) * V).
%   - Parallelized across voxels with `parfor`.
%


    % ---- Split b-values (masks & vectors) ----
    mask_hi = (bvals > bsplit);
    mask_lo = (bvals < bsplit);
    nonivimb = bvals(mask_hi);
    ivimb    = bvals(mask_lo);
    %n1 = numel(nonivimb);
    %n2 = numel(ivimb);

    % ---- Grids ----
    w2D = linspace(D_min, D_max, number_of_points1);   % D
    w2V = linspace(v_min, v_max, number_of_points2);   % V

    % ---- Precompute exponent tables (outer products, correct shapes) ----
    % High-b (diffusion-only): exp(-b * D)
    E_hi = exp( - (nonivimb(:)) * (w2D(:).') );        % [n1 x nD]
    s1_hi = sum(E_hi.^2, 1).';                         % [nD x 1]

    % Low-b (vascular term): exp(-Δ * sqrt(b/(Δ-δ/3)) * V)
    c  = large_delta .* sqrt( ivimb(:) ./ (large_delta - small_delta/3) ); % [n2 x 1]
    E_lo  = exp( - c * (w2V(:).') );                   % [n2 x nV]
    s1_lo = sum(E_lo.^2, 1).';                         % [nV x 1]

    % ---- Progress + output allocation ----
    %progbar = progressBar(size(to_calculation,2),'pname','Calculating grid search');
    calculated_values = zeros(length(to_calculation),4);

    parfor i = 1:length(to_calculation)
        %progbar.progress
        nonivimy = to_calculation(mask_hi, i);
        ivimy    = to_calculation(mask_lo, i);

        % ---------- Step 1: Fit D, S0 from high b-values ----------
        max_attenuation = max(nonivimy)/max(ivimy);
        top_signal      = max(nonivimy)/max_attenuation;

        w1 = linspace(0.75*top_signal, top_signal, number_of_points1); % S0 candidates [1 x nS0]

        y   = nonivimy(:);                     % [n1 x 1]
        y2  = sum(y.^2);                       % scalar
        s2  = E_hi.' * y;                      % [nD x 1], sum(y.*E(:,j))

        SSE_hi = y2 - 2*(s2 * w1) + (s1_hi * (w1.^2));   % [nD x nS0]

        SSE_min = min(SSE_hi(:));
        ind_all = find(SSE_hi == SSE_min);
        ind     = round(median(ind_all));                 % tie-break to match original
        [xindex,yindex] = ind2sub(size(SSE_hi), ind);

        S0 = w1(yindex);
        Dp = w2D(xindex);

        % ---------- Step 2: Fit vascular V, f from residual low b-values ----------
        onlyivimy = ivimy - S0 * exp( - Dp * ivimb(:) );  % [n2 x 1]
        max_ivim_signal = max(onlyivimy);

        if max_ivim_signal > 0.02*S0
            w1b = linspace(0, 1.5*max_ivim_signal, number_of_points2); % S0_ivim candidates

            yb   = onlyivimy(:);               % [n2 x 1]
            yb2  = sum(yb.^2);                 % scalar
            s2lo = E_lo.' * yb;                % [nV x 1]

            SSE_lo = yb2 - 2*(s2lo * w1b) + (s1_lo * (w1b.^2));  % [nV x nS0b]

            SSE_min2 = min(SSE_lo(:));
            ind_all2 = find(SSE_lo == SSE_min2);
            ind2     = round(median(ind_all2));
            [xindex2,yindex2] = ind2sub(size(SSE_lo), ind2);

            S0iv = w1b(yindex2);
            V    = w2V(xindex2);
            f    = S0iv / (S0 + S0iv);

            calculated_values(i,:) = [ (S0 + S0iv), f, V, Dp ];
        else
            calculated_values(i,:) = [ S0, 0, 0, Dp ];
        end
    end
end
