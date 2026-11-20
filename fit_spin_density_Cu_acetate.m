%% FIT_SPIN_DENSITY_CU_ACETATE
%  Fit a real-space spin-density model to the 3D inelastic neutron
%  scattering (INS) intensity of the triplet excitation in copper acetate.
%
% PURPOSE
%   Copper acetate, Cu2(CH3COO)4 (Cu2Ac4), is a system of isolated Cu(II)
%   dimers with a singlet ground state and a triplet excitation. This
%   script fits the intensity I(Q) of the triplet excitation, measured in
%   a neutron time-of-flight experiment on a single crystal and stored on a
%   3D grid in reciprocal space, with the Fourier transform of a real-space
%   spin-density model made of anisotropic Gaussians centred on the atoms.
%   It is a first attempt to reconstruct the real-space spin density
%   directly from INS data.
%
%   The output of this script is shown in Fig. 3 of
%   https://arxiv.org/abs/2604.13751
%
% MODEL
%   * The unit cell contains two differently oriented dimers:
%     Cu1-Cu2 and Cu3-Cu4. Each Cu carries a cluster of 5 atoms: the Cu
%     itself and its 4 nearest oxygens.
%   * Each atom carries a normalised, axially symmetric 3D Gaussian spin
%     density, so its Fourier transform is
%         f(Q) = exp(-1/2 * (s_perp^2*Q_perp^2 + s_par^2*Q_par^2)),
%     where Q_par is the component of Q along the atom's symmetry axis:
%         Cu : the Cu-Cu direction of its dimer
%         O  : the O -> Cu bond direction
%   * Spin is shared between the atoms of a cluster: a fraction rho_Cu on
%     the Cu and (1 - rho_Cu)/4 on each of the 4 oxygens.
%   * Cluster structure factor: S_c(Q) = sum_atoms amp * f(Q) * exp(i Q.r)
%   * Triplet-excitation intensity of a dimer (Cu_a, Cu_b) is proportional
%     to |S_a(Q) - S_b(Q)|^2 (singlet-triplet transition); the two dimers
%     in the cell add incoherently:
%         I(Q) = scale * ( |S1-S2|^2 + |S3-S4|^2 ) + background
%     No polarisation factor or separate magnetic form factor is applied.
%   * The measured intensity is an average over a voxel of the 3D grid. The
%     model is therefore evaluated on nsub^3 points per voxel and averaged.
%
% FREE PARAMETERS (vector p, lengths in Angstrom)
%   p(1) rho_Cu      fraction of the spin on each Cu 
%   p(2) sigma_Cu_perp   Cu Gaussian width perpendicular to the Cu-Cu axis
%                        (i.e. within the CuO4 plane)
%   p(3) sigma_Cu_par    Cu Gaussian width along the Cu-Cu axis
%   p(4) sigma_O_perp    O  Gaussian width perpendicular to the O-Cu bond
%   p(5) sigma_O_par     O  Gaussian width along the O-Cu bond
%   p(6) scale       overall intensity scale factor
%   The flat background is fixed in this implementation (cfg.background, default 0).
%
% COORDINATES
%   The data axes (H, K, L) are NOT Miller indices. A point of the grid
%   corresponds to the reciprocal-lattice vector
%         hkl = H*u1 + L*u2 + K*[0 1 0]        (r.l.u.)
%   so H and L run along the in-plane directions u1 and u2 (cfg.u1, cfg.u2)
%   and K along b*. The Cartesian momentum transfer is Q = B*hkl' with
%   B = 2*pi*inv(A), A being the matrix of lattice vectors (rows).
%
% INPUT DATA (folder cfg.dataDir, one .mat file per variable)
%   I_background_subtracted_dQ_<dQ>.mat      variable 'sub'      3D intensity (Qx,Qy,Qz)
%   I_err_background_subtracted_dQ_<dQ>.mat  variable 'sub_err'  its uncertainty
%   Qx_background_subtracted_dQ_<dQ>.mat     variable 'Qx'       grid axis (-> H)
%   Qy_background_subtracted_dQ_<dQ>.mat     variable 'Qy'       grid axis (-> L)
%   Qz_background_subtracted_dQ_<dQ>.mat     variable 'Qz'       grid axis (-> K)
%   The data are background-subtracted. Empty voxels have zero intensity
%   or zero error and are excluded from the fit.
%
% OUTPUT
%   Fitted parameters with uncertainties (command window, and optionally
%   cfg.resultsFile), and figures comparing data and best-fit model in the
%   K = 0 plane and along three 1D cuts.
%
% REQUIREMENTS
%   MATLAB R2022a or newer (clim; R2017b is enough without the plots),
%   Optimization Toolbox (lsqnonlin). Parallel Computing Toolbox is
%   optional: parfor runs serially without it.

clear; close all; clc

%% 1. Configuration -------------------------------------------------------

cfg.dataDir = fullfile(pwd,'data');   % folder containing the .mat files
cfg.dQ      = 0.15;     % voxel size of the data grid, in units of the
                        % H, K, L axes (also part of the file names)

% In-plane directions of the data grid in r.l.u. (see COORDINATES above)
cfg.u1 = [0.9091, 0,  1];
cfg.u2 = [0.9786, 0, -1];

% Region of the grid that is loaded: [min max] on each axis
cfg.crop.H = [-7.5,  7.5];
cfg.crop.K = [-1.4,  2.8];
cfg.crop.L = [-10.5, 10.5];

% Regions excluded from the fit (defined in the H-L plane, all K)
cfg.mask.beam     = [1.9, 3.1];     % ellipse semi-axes [H, L] around the
                                    % direct beam (inner exclusion)
cfg.mask.outer    = [7, 9.5];       % ellipse semi-axes [H, L]; everything
                                    % outside is excluded
cfg.mask.spurions = [4.65, -4.75, 0.4, 0.75];
                                    % rows: [H0 L0 radiusH radiusL]. The
                                    % point-symmetric ellipse at (-H0,-L0)
                                    % is excluded as well.
cfg.mask.corner   = [2.1, 5];       % [|H|, |L|]: opposite corners, to mask other spurions
                                    % (H > 2.1, L < -5) and (H < -2.1, L > 5)


% Model evaluation
cfg.nsub       = 3;     % voxel sub-sampling nsub x nsub x nsub
cfg.background = 0;     % fixed flat background

% Fit: parameter order as in the header
cfg.fit.names = {'rho_Cu', 'sigma_Cu_perp', 'sigma_Cu_par', ...
                 'sigma_O_perp', 'sigma_O_par', 'scale'};
cfg.fit.p0 = [0.7, 0.2, 0.15, 0.25, 0.40, 6.5e-5];   % start values
cfg.fit.lb = [0,    0.001, 0.001, 0.001, 0.001, 0];    % lower bounds
cfg.fit.ub = [1,    2,     2,     2,     2,     1];   % upper bounds
cfg.fit.maxIter = 30;

% Output
cfg.resultsFile = 'fit_results.mat';   % set to '' to skip saving
cfg.plot.climMax = 3e-4;               % colour scale of the intensity maps
cfg.plot.nsub    = cfg.nsub;           % sub-sampling of the plotted model;
                                       % 1 = model at voxel centres only

% 1D cuts shown in the figures (the nearest grid point is used)
cfg.cut.H0 = 2;         % cut along L at this H
cfg.cut.L0 = 3.3;       % cut along H at this L
cfg.cut.K  = [-2, 2];   % cut along K at this [H, L]

%% 2. Load data and define the voxels used in the fit ---------------------

[H, K, L, Z, Zerr] = load_dataset(cfg);          % Z, Zerr: nH x nK x nL

mask = build_mask(H, L, Z, Zerr, cfg.mask);      % true = voxel is fitted
fprintf('Voxels used in the fit: %d of %d\n', nnz(mask), numel(mask));

geom = build_geometry(cfg);

% Experimental data as vectors, restricted to the fitted voxels
[iH, iK, iL] = ind2sub(size(mask), find(mask));
hFit = H(iH);  kFit = K(iK);  lFit = L(iL);
yExp = Z(mask);
w    = 1 ./ Zerr(mask);                          % weights 1/sigma

%% 3. Fit -----------------------------------------------------------------

model       = @(p) voxel_average(p, hFit, kFit, lFit, geom, cfg.dQ, cfg.nsub);
residualFun = @(p) w .* (yExp - model(p));

opts = optimoptions('lsqnonlin', 'Display', 'iter', ...
    'MaxIterations', cfg.fit.maxIter, 'FunctionTolerance', 1e-10);
[pFit, chi2, resid, exitflag, ~, ~, J] = ...
    lsqnonlin(residualFun, cfg.fit.p0, cfg.fit.lb, cfg.fit.ub, opts);

if exitflag <= 0
    warning('lsqnonlin did not converge (exitflag = %d). Increase maxIter?', exitflag);
end

% Parameter uncertainties from the Jacobian, with the covariance rescaled
% by the reduced chi^2. (Not meaningful for parameters sitting on a bound.)
nData  = numel(resid);
nPar   = numel(pFit);
redChi2 = chi2 / (nData - nPar);
J      = full(J);
covP   = redChi2 * ((J' * J) \ eye(nPar));
pErr   = sqrt(diag(covP));

fprintf('\nReduced chi^2 = %.3f  (%d points, %d parameters)\n', redChi2, nData, nPar);
disp(table(cfg.fit.names(:), pFit(:), pErr(:), ...
    'VariableNames', {'Parameter', 'FittedValue', 'Uncertainty'}));

if ~isempty(cfg.resultsFile)
    results = struct('names', {cfg.fit.names}, 'p', pFit, 'pErr', pErr, ...
                     'redChi2', redChi2, 'cfg', cfg);
    save(cfg.resultsFile, 'results');
end

%% 4. Best-fit model on the full grid -------------------------------------

[Hg, Kg, Lg] = ndgrid(H, K, L);
bestFit = reshape(voxel_average(pFit, Hg(:), Kg(:), Lg(:), geom, ...
                                cfg.dQ, cfg.plot.nsub), size(Z));
clear Hg Kg Lg

%% 5. Figures -------------------------------------------------------------

% Axis lengths in Angstrom^-1, used to label the plot axes
lenU1 = norm(geom.B * geom.u1');
lenU2 = norm(geom.B * geom.u2');
lenB  = norm(geom.B(:, 2));
labH  = [cut_label([0 0 0], geom.u1, 'h') sprintf(' in %.4f Angstrom^{-1}', lenU1)];
labL  = [cut_label([0 0 0], geom.u2, 'h') sprintf(' in %.4f Angstrom^{-1}', lenU2)];

[~, iK0] = min(abs(K));            % equatorial plane (K = 0)
[~, iH0] = min(abs(H));
[~, iL0] = min(abs(L));
clim0 = [0, cfg.plot.climMax];

% 5a. Overview of the data: three orthogonal slices through the origin
figure('Name', 'Data overview');
subplot(1, 3, 1); imagesc(H, L, squeeze(Z(:, iK0, :))'); axis xy; clim(clim0)
xlabel(labH); ylabel(labL); title('K = 0')
subplot(1, 3, 2); imagesc(H, K, squeeze(Z(:, :, iL0))'); axis xy; clim(clim0)
xlabel(labH); ylabel('K'); title('L = 0')
subplot(1, 3, 3); imagesc(L, K, squeeze(Z(iH0, :, :))); axis xy; clim(clim0)
xlabel(labL); ylabel('K'); title('H = 0')

% 5b. Data, best-fit model and residuals in the K = 0 plane
planeExp = squeeze(Z(:, iK0, :));
planeFit = squeeze(bestFit(:, iK0, :));
figure('Name', 'Fit, K = 0 plane');
subplot(1, 3, 1); imagesc(H, L, planeExp');  axis xy; clim(clim0); title('Experiment')
xlabel(labH); ylabel(labL)
subplot(1, 3, 2); imagesc(H, L, planeFit');  axis xy; clim(clim0); title('Fit')
xlabel(labH); ylabel(labL)
subplot(1, 3, 3); imagesc(H, L, (planeExp - planeFit)'); axis xy
clim(cfg.plot.climMax * [-1, 1]); title('Residuals (experiment - fit)')
xlabel(labH); ylabel(labL)

% 5c. Cut along L (along u2), at fixed H
[~, ih] = min(abs(H - cfg.cut.H0));
lab = [cut_label(H(ih) * geom.u1, geom.u2, 't') ...
       sprintf(', t in units of %.4f Angstrom^{-1}', lenU2)];
plot_cut(L, Z(ih, iK0, :), Zerr(ih, iK0, :), mask(ih, iK0, :), bestFit(ih, iK0, :), lab);

% 5d. Cut along H (along u1), at fixed L
[~, il] = min(abs(L - cfg.cut.L0));
lab = [cut_label(L(il) * geom.u2, geom.u1, 't') ...
       sprintf(', t in units of %.4f Angstrom^{-1}', lenU1)];
plot_cut(H, Z(:, iK0, il), Zerr(:, iK0, il), mask(:, iK0, il), bestFit(:, iK0, il), lab);

% 5e. Cut along K (out of plane), at fixed H and L
[~, ih] = min(abs(H - cfg.cut.K(1)));
[~, il] = min(abs(L - cfg.cut.K(2)));
lab = [cut_label(H(ih) * geom.u1 + L(il) * geom.u2, [0 1 0], 't') ...
       sprintf(', t in units of %.4f Angstrom^{-1}', lenB)];
plot_cut(K, Z(ih, :, il), Zerr(ih, :, il), mask(ih, :, il), bestFit(ih, :, il), lab);


%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

function [H, K, L, Z, Zerr] = load_dataset(cfg)
% LOAD_DATASET  Read the 3D intensity map and crop it to cfg.crop.
%   Returns the axes H, K, L (column vectors) and the intensity Z and its
%   error Zerr, both ordered (H, K, L). The files store the data ordered
%   (Qx, Qy, Qz), which corresponds to (H, L, K).
    tag = sprintf('dQ_%g.mat', cfg.dQ);
    f   = @(prefix) fullfile(cfg.dataDir, [prefix tag]);

    S1 = load(f('I_background_subtracted_'),     'sub');
    S2 = load(f('I_err_background_subtracted_'), 'sub_err');
    Sx = load(f('Qx_background_subtracted_'),    'Qx');
    Sy = load(f('Qy_background_subtracted_'),    'Qy');
    Sz = load(f('Qz_background_subtracted_'),    'Qz');

    Qx = Sx.Qx(:);  Qy = Sy.Qy(:);  Qz = Sz.Qz(:);
    ix = crop_indices(Qx, cfg.crop.H);
    iy = crop_indices(Qy, cfg.crop.L);
    iz = crop_indices(Qz, cfg.crop.K);

    H = Qx(ix);  L = Qy(iy);  K = Qz(iz);
    Z    = permute(S1.sub(ix, iy, iz),     [1 3 2]);   % (Qx,Qy,Qz) -> (H,K,L)
    Zerr = permute(S2.sub_err(ix, iy, iz), [1 3 2]);
end

function idx = crop_indices(q, range)
% CROP_INDICES  Contiguous indices of the (sorted) axis q within
%   lower < q <= upper.
    inside = find(q > range(1) & q <= range(2));
    idx = inside(1):inside(end);
end

function mask = build_mask(H, L, Z, Zerr, m)
% BUILD_MASK  Logical 3D array, true for the voxels used in the fit.
%   A voxel is used if it contains data (non-zero intensity and error) and
%   lies outside the regions excluded by hand (see cfg.mask). The excluded
%   regions are defined in the H-L plane and apply to all K.
    [Hg, Lg] = ndgrid(H, L);
    ellipse  = @(h0, l0, rh, rl) ((Hg - h0) / rh).^2 + ((Lg - l0) / rl).^2;

    bad = (Hg >  m.corner(1) & Lg < -m.corner(2)) | ...
          (Hg < -m.corner(1) & Lg >  m.corner(2));       % opposite corners
    bad = bad | ellipse(0, 0, m.beam(1),  m.beam(2))  < 1;   % direct beam
    bad = bad | ellipse(0, 0, m.outer(1), m.outer(2)) > 1;   % outer region
    for s = 1:size(m.spurions, 1)                            % spurions
        bad = bad | ellipse( m.spurions(s,1),  m.spurions(s,2), ...
                             m.spurions(s,3),  m.spurions(s,4)) < 1;
        bad = bad | ellipse(-m.spurions(s,1), -m.spurions(s,2), ...
                             m.spurions(s,3),  m.spurions(s,4)) < 1;
    end

    hasData = Z ~= 0 & Zerr ~= 0 & isfinite(Z) & isfinite(Zerr);
    mask = hasData & ~reshape(bad, size(bad, 1), 1, size(bad, 2));
end

function geom = build_geometry(cfg)
% BUILD_GEOMETRY  Precompute everything in the model that does not depend
%   on the fit parameters.
%   geom.B        3x3 matrix, Q = B*hkl' (columns are a*, b*, c*)
%   geom.u1, u2   in-plane directions of the grid (row vectors)
%   geom.cluster  1x4 struct array, one per Cu (Cu1..Cu4), with fields
%       .pos   5x3  positions (Angstrom) of the Cu and its 4 nearest O
%       .axis  5x3  unit symmetry axis of each atom's Gaussian
%                   (Cu: along the Cu-Cu direction; O: along the O->Cu bond)
%       .isCu  5x1  true for the Cu, false for the oxygens
%   Clusters 1,2 form the first dimer and clusters 3,4 the second.
    [pos, elem, A] = copper_acetate_structure();
    unit = @(v) v / norm(v);

    geom.B  = 2 * pi * inv(A);
    geom.u1 = cfg.u1(:).';
    geom.u2 = cfg.u2(:).';
    geom.background = cfg.background;

    cuIdx = find(strcmp(elem, 'Cu'));
    oIdx  = find(strcmp(elem, 'O'));
    dimers = [1 2; 3 4];                  % pairs of Cu forming a dimer

    for d = 1:size(dimers, 1)
        pair = cuIdx(dimers(d, :));
        dimerAxis = unit(pos(pair(2), :) - pos(pair(1), :));

        for m = 1:2
            c  = 2 * (d - 1) + m;
            cu = pair(m);

            % four nearest oxygens of this Cu
            dist = vecnorm(pos(oIdx, :) - pos(cu, :), 2, 2);
            [dist, order] = sort(dist);
            assert(all(dist(1:4) < 2.5) && dist(5) > 2.5, ...
                   'Unexpected Cu-O coordination for Cu %d', c);
            oNear = oIdx(order(1:4));

            axes = zeros(5, 3);
            axes(1, :) = dimerAxis;
            for s = 1:4
                axes(s + 1, :) = unit(pos(cu, :) - pos(oNear(s), :));
            end
            geom.cluster(c).pos  = [pos(cu, :); pos(oNear, :)];
            geom.cluster(c).axis = axes;
            geom.cluster(c).isCu = [true; false(4, 1)];
        end
    end
end

function [pos, elem, A] = copper_acetate_structure()
% COPPER_ACETATE_STRUCTURE  Unit cell and atomic positions of Cu2Ac4.
%   A    lattice vectors as rows (Angstrom), a along x, b along y, c in x-z
%        (monoclinic cell).
%   pos  Cartesian positions (Angstrom) of the 4 Cu and 16 O of the cell
%        that enter the model, from single-crystal x-ray diffraction.
%   elem cell array with the element of each row of pos.
    a = 13.0866;  b = 8.4743;  c = 13.6453;
    beta = 116.8909 * pi / 180;
    A = [a,            0, 0;
         0,            b, 0;
         c*cos(beta),  0, c*sin(beta)];

    pos = [
        2.51732   3.51760   6.63122;   % Cu1
        4.39761   4.95670   5.53861;   % Cu2
        7.48344   3.51760  -0.54630;   % Cu3
        5.60316   4.95670   0.54630;   % Cu4
        1.72546   5.26373   7.00872;   % four O near Cu1
        3.55022   1.97680   6.06690;
        3.68044   3.71183   8.22936;
        1.65288   3.47379   4.85260;
        3.36471   6.49750   6.10292;   % four O near Cu2
        5.18947   3.21057   5.16110;
        5.26205   5.00052   7.31723;
        3.23448   4.76247   3.94047;
        8.27530   5.26373  -0.92381;   % four O near Cu3
        6.45054   1.97680   0.01801;
        6.32032   3.71183  -2.14444;
        8.34788   3.47379   1.23232;
        6.63606   6.49750  -0.01801;   % four O near Cu4
        4.81130   3.21057   0.92381;
        4.73872   5.00052  -1.23232;
        6.76628   4.76247   2.14444];

    elem = [repmat({'Cu'}, 4, 1); repmat({'O'}, 16, 1)];
end

function par = unpack_params(p)
% UNPACK_PARAMS  Name the entries of the fit-parameter vector (see header).
    par.rhoCu   = p(1);            % spin fraction on each Cu
    par.rhoO    = (1 - p(1)) / 4;  % spin fraction on each of the 4 O
    par.sCuPerp = p(2);            % Cu width perpendicular to Cu-Cu axis
    par.sCuPar  = p(3);            % Cu width along Cu-Cu axis
    par.sOPerp  = p(4);            % O width perpendicular to O-Cu bond
    par.sOPar   = p(5);            % O width along O-Cu bond
    par.scale   = p(6);
end

function I = model_intensity(p, h, k, l, geom)
% MODEL_INTENSITY  Model intensity at the points (h, k, l) (column vectors
%   of equal length, in the units of the data axes):
%       I = scale * (|S1-S2|^2 + |S3-S4|^2) + background
    par = unpack_params(p);
    hkl = h(:) * geom.u1 + l(:) * geom.u2 + k(:) * [0 1 0];   % N x 3, r.l.u.
    Q   = hkl * geom.B.';                                     % N x 3, 1/Angstrom

    S = zeros(numel(h), numel(geom.cluster));
    for c = 1:numel(geom.cluster)
        S(:, c) = cluster_structure_factor(Q, geom.cluster(c), par);
    end
    I = par.scale * (abs(S(:,1) - S(:,2)).^2 + abs(S(:,3) - S(:,4)).^2) ...
        + geom.background;
end

function S = cluster_structure_factor(Q, cl, par)
% CLUSTER_STRUCTURE_FACTOR  Fourier transform of the spin density of one
%   Cu + 4 O cluster, S(Q) = sum_atoms amp * f(Q) * exp(i Q.r).
%   Q is N x 3 (1/Angstrom); returns an N x 1 complex vector.
%   Each atom is an axially symmetric Gaussian, so the form factor only
%   depends on the components of Q along (Q_par) and perpendicular to
%   (Q_perp) its symmetry axis:
%       f = exp(-1/2 * (s_perp^2 * Q_perp^2 + s_par^2 * Q_par^2))
    S  = zeros(size(Q, 1), 1);
    Q2 = sum(Q.^2, 2);
    for s = 1:numel(cl.isCu)
        if cl.isCu(s)
            amp = par.rhoCu;  sPerp = par.sCuPerp;  sPar = par.sCuPar;
        else
            amp = par.rhoO;   sPerp = par.sOPerp;   sPar = par.sOPar;
        end
        qPar = Q * cl.axis(s, :).';
        f    = exp(-0.5 * (sPerp^2 * (Q2 - qPar.^2) + sPar^2 * qPar.^2));
        S    = S + amp * f .* exp(1i * (Q * cl.pos(s, :).'));
    end
end

function I = voxel_average(p, h, k, l, geom, dQ, nsub)
% VOXEL_AVERAGE  Model intensity averaged over a voxel of size dQ, using
%   nsub^3 equally weighted sample points per voxel (nsub = 1: voxel centre).
%   h, k, l are the voxel centres (column vectors).
%   Note: the sample points span the full voxel including its faces
%   (offsets -dQ/2 ... +dQ/2)
    if nsub == 1
        offsets = 0;
    else
        offsets = linspace(-0.5, 0.5, nsub) * dQ;
    end
    n    = numel(offsets);
    nTot = n^3;
    I    = zeros(size(h));
    parfor idx = 1:nTot
        [ia, ib, ic] = ind2sub([n n n], idx);
        I = I + model_intensity(p, h + offsets(ia), k + offsets(ib), ...
                                   l + offsets(ic), geom);
    end
    I = I / nTot;
end

function s = cut_label(base, dir, var)
% CUT_LABEL  Text label of the line  hkl = base + dir * var  in r.l.u.,
%   e.g. '[1.818+0.9786t, 0, 2-1t]'.
    c = cell(1, 3);
    for i = 1:3
        if dir(i) == 0
            c{i} = sprintf('%.4g', base(i));
        elseif base(i) == 0
            c{i} = sprintf('%.4g%s', dir(i), var);
        else
            c{i} = sprintf('%.4g%+.4g%s', base(i), dir(i), var);
        end
    end
    s = ['[' strjoin(c, ', ') ']'];
end

function plot_cut(x, y, yerr, used, ymodel, xlab)
% PLOT_CUT  1D cut: data used in the fit (blue circles), masked data (red
%   crosses) and best-fit model (dashed line).
    x = x(:);  y = y(:);  yerr = yerr(:);  ymodel = ymodel(:);
    used = logical(used(:));
    figure; hold on; box on
    errorbar(x(used),  y(used),  yerr(used),  'bo');
    errorbar(x(~used), y(~used), yerr(~used), 'rx');
    plot(x, ymodel, 'k--', 'LineWidth', 2);
    xlabel(xlab); ylabel('Intensity');
    legend({'data (fitted)', 'data (masked)', 'model'}, 'Location', 'best');
end
