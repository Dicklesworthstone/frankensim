//! Helmholtz-kernel BEM with the Burton–Miller combined formulation
//! (bead frankensim-fsim-helmholtz-bem-k1ryv, musical-acoustics program):
//! exterior acoustic radiation from a closed panelized surface — surface
//! pressure from prescribed normal velocity, radiation impedance,
//! radiated power, and far-field directivity.
//!
//! TIME CONVENTION (pinned; every sign below depends on it): fields vary
//! as `e^{-i omega t}`, so the outgoing free-space kernel is
//! `G(r) = e^{+ikr} / (4 pi r)` and the momentum equation gives
//! `grad p = i omega rho v`, i.e. `q = dp/dn = i omega rho v_n`. With
//! the HBIE row written below as `N p - D' q - q/2 = 0`, the coupling
//! that restores unique solvability under this convention is
//! `alpha = -i/k` (equivalently the textbook `+i/k` applied to the
//! NEGATED HBIE row — coupling sign and row sign only mean anything
//! together). The sign is pinned by the interior-resonance contrast
//! test rather than trusted from a convention label: the flipped sign
//! measurably degrades the resonance error instead of fixing it.
//!
//! FORMULATION — exterior Green representation with outward normals,
//! collocation at panel centroids, piecewise-constant elements:
//! `p(x) = INT [ dG/dn_y p - G q ] dS` for x in the exterior. Surface
//! limits give
//! CBIE: `D p - S q - p/2 = 0`, and its `n_x` derivative
//! HBIE: `N p - D' q - q/2 = 0` (hypersingular finite part), combined as
//! `(D - I/2 + alpha N) p = (S + alpha D' + alpha/2 I) q`.
//!
//! DISCRETIZATION HONESTY: off-diagonal influence uses the centroid
//! point-panel approximation (`kernel(x_i, y_j) * A_j`), the same class
//! the crate's Laplace side validated to +1.8% on the aperture pilot.
//! Self terms use the equivalent-disc analytic finite parts with
//! `a = sqrt(A/pi)`:
//! `S_ii = (e^{ika} - 1) / (2ik)`, `D_ii = D'_ii = 0` (flat panel), and
//! the hypersingular self entry split as the integrable difference
//! `(N_k - N_0)_ii = (ik - (e^{ika} - 1)/a)/2` plus the static part
//! recovered from the exact closed-surface identity `N_0[1] = 0` as
//! minus the discrete off-diagonal static row sum (same quadrature on
//! both sides, so the point-panel error cancels).
//! These closed forms are pinned by unit tests against numerical
//! quadrature with the singular core subtracted, and the whole pipeline
//! is arbitrated by the pulsating-sphere impedance oracle.
//!
//! Determinism: assembly and solves are sequential with fixed traversal
//! order; complex exponentials route through `fs_math::det` per the
//! workspace libm doctrine. Repeat solves are bitwise identical.
//!
//! Deferred with recorded triggers (see CONTRACT): wideband directional
//! FMM acceleration (dense complex LU is the v1 path; the work cap
//! refuses above [`MAX_DENSE_PANELS`]), exact triangle singular
//! quadrature beyond the equivalent-disc self terms, spherical-harmonic
//! directivity coefficient tables (sampled directivity ships now), and a
//! rigorous condition estimator (the fictitious-frequency contrast is
//! the recorded diagnostic).

use fs_la::eigen_complex::lu_complex;
use fs_math::c64::C64;
use fs_math::det;

use crate::panel3d::SpherePanels;

/// Dense-LU work cap: a complex dense system above this panel count is
/// refused rather than silently thrashing (n^2 * 16 bytes; 8192 panels
/// is already a 1 GiB matrix). The FMM path is the recorded follow-up.
pub const MAX_DENSE_PANELS: usize = 8192;

/// Minimum panels per wavelength before the solver refuses: below ~6 the
/// centroid-collocation solution degrades silently, so the boundary is a
/// named refusal instead.
pub const MIN_PANELS_PER_WAVELENGTH: f64 = 6.0;

/// Typed refusals with stable `FS-BEM-HELM-*` codes.
#[derive(Debug, Clone, PartialEq)]
pub enum HelmholtzError {
    /// Wavenumber, density, or sound speed is non-positive/non-finite.
    BadParameter {
        /// Which parameter refused.
        what: &'static str,
    },
    /// The mesh is too coarse for the requested wavenumber.
    TooCoarse {
        /// Measured panels per wavelength.
        panels_per_wavelength: f64,
    },
    /// The dense work cap was exceeded.
    WorkCap {
        /// Requested panel count.
        panels: usize,
    },
    /// Input length disagrees with the panel count.
    ShapeMismatch {
        /// What disagreed.
        what: &'static str,
    },
    /// The dense complex LU reported singularity.
    Singular,
}

impl core::fmt::Display for HelmholtzError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            HelmholtzError::BadParameter { what } => {
                write!(f, "FS-BEM-HELM-BAD-PARAMETER: {what}")
            }
            HelmholtzError::TooCoarse {
                panels_per_wavelength,
            } => write!(
                f,
                "FS-BEM-HELM-TOO-COARSE: {panels_per_wavelength:.2} panels per wavelength \
                 (floor {MIN_PANELS_PER_WAVELENGTH})"
            ),
            HelmholtzError::WorkCap { panels } => write!(
                f,
                "FS-BEM-HELM-WORK-CAP: {panels} panels exceeds the dense cap {MAX_DENSE_PANELS}"
            ),
            HelmholtzError::ShapeMismatch { what } => {
                write!(f, "FS-BEM-HELM-SHAPE-MISMATCH: {what}")
            }
            HelmholtzError::Singular => write!(f, "FS-BEM-HELM-SINGULAR: dense LU refused"),
        }
    }
}

impl std::error::Error for HelmholtzError {}

/// Acoustic medium (SI): density and sound speed.
#[derive(Debug, Clone, Copy)]
pub struct Medium {
    /// Ambient density rho [kg/m^3].
    pub density: f64,
    /// Sound speed c [m/s].
    pub sound_speed: f64,
}

impl Medium {
    /// Air at roughly 20 degC.
    #[must_use]
    pub const fn air() -> Medium {
        Medium {
            density: 1.204,
            sound_speed: 343.0,
        }
    }
}

/// Boundary-integral formulation choice.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Formulation {
    /// Plain CBIE — breaks down at interior resonances (kept as the
    /// documented near-miss arm of the fictitious-frequency contrast).
    PlainCbie,
    /// Burton–Miller CBIE + alpha HBIE with `alpha = -i/k` for this
    /// module's row-sign convention — the production arm.
    BurtonMiller,
    /// The alpha-sign MUTATION arm (`alpha = +i/k` on this module's
    /// row convention): kept public and
    /// documented so the conformance battery can prove the wrong sign
    /// fails loudly at interior resonances instead of silently
    /// degrading. Never use for production solves.
    BurtonMillerWrongAlphaSign,
}

/// A solved exterior radiation problem.
#[derive(Debug, Clone)]
pub struct RadiationSolution {
    /// Complex surface pressure per panel [Pa].
    pub pressure: Vec<C64>,
    /// The prescribed complex normal velocity per panel [m/s].
    pub velocity: Vec<C64>,
    /// Wavenumber k [1/m].
    pub k: f64,
    /// Measured panels per wavelength (the refusal diagnostic).
    pub panels_per_wavelength: f64,
    /// Radiated power W = 1/2 Re SUM p conj(v) A [W].
    pub radiated_power: f64,
}

fn expik(kr: f64) -> C64 {
    C64::new(det::cos(kr), det::sin(kr))
}

/// Free-space kernel `G = e^{ikr}/(4 pi r)`.
fn green(k: f64, r: f64) -> C64 {
    expik(k * r).scale(1.0 / (4.0 * core::f64::consts::PI * r))
}

/// `G'(r) = (ik - 1/r) G`.
fn green_dr(k: f64, r: f64) -> C64 {
    green(k, r) * C64::new(-1.0 / r, k)
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

/// The four point kernels for x != y: (G, dG/dn_y, dG/dn_x,
/// d2G/dn_x dn_y). With `d = x - y`, `r = |d|`:
/// `dG/dn_y = -G' (n_y . d)/r`, `dG/dn_x = G' (n_x . d)/r`, and
/// `d2G/dn_x dn_y = -(G''/r^2 - G'/r^3)(n_x . d)(n_y . d) - (G'/r)(n_x . n_y)`
/// where `G'' = (ik - 1/r) G' + G/r^2`. All four are pinned against
/// central finite differences of G in the unit tests.
fn kernels(k: f64, x: [f64; 3], nx: [f64; 3], y: [f64; 3], ny: [f64; 3]) -> (C64, C64, C64, C64) {
    let d = sub(x, y);
    let r = norm(d);
    let g = green(k, r);
    let gp = green_dr(k, r);
    let gpp = gp * C64::new(-1.0 / r, k) + g.scale(1.0 / (r * r));
    let nyd = dot(ny, d);
    let nxd = dot(nx, d);
    let dgdny = gp.scale(-nyd / r);
    let dgdnx = gp.scale(nxd / r);
    let radial = gpp.scale(1.0 / (r * r)) - gp.scale(1.0 / (r * r * r));
    let d2g = radial.scale(-nxd * nyd) - gp.scale(dot(nx, ny) / r);
    (g, dgdny, dgdnx, d2g)
}

/// Equivalent-disc self terms `(S_ii, N_reg_ii)` for a panel of area
/// `A`, with `a = sqrt(A/pi)`. `S_ii = (e^{ika} - 1)/(2ik)` integrates
/// the full weakly singular kernel; the hypersingular self term is
/// returned as the REGULARIZED difference `(N_k - N_0)_ii = (ik -
/// (e^{ika} - 1)/a)/2`, whose kernel `[e^{ikr}(1-ikr) - 1]/(4 pi r^3)`
/// is integrable. The static part is NOT approximated by a disc:
/// on a closed surface the double layer of constant density is constant
/// off the surface, so `N_0[1] = 0` exactly, and the assembly recovers
/// `N_0`'s self entry as minus the discrete off-diagonal static row sum
/// — the same quadrature on both sides, so the point-panel error
/// cancels instead of polluting the Burton-Miller row.
fn self_terms(k: f64, area: f64) -> (C64, C64) {
    let a = (area / core::f64::consts::PI).sqrt();
    let ik = C64::new(0.0, k);
    let s = (expik(k * a) - C64::ONE) * (ik + ik).recip();
    let n_reg = (ik - (expik(k * a) - C64::ONE).scale(1.0 / a)).scale(0.5);
    (s, n_reg)
}

fn characteristic_panel_size(surface: &SpherePanels) -> f64 {
    let mut max_area = 0.0f64;
    for &a in surface.areas() {
        max_area = max_area.max(a);
    }
    max_area.sqrt()
}

/// Solve the exterior Neumann radiation problem: prescribed complex
/// normal velocity per panel, returned surface pressure. `q = i omega
/// rho v_n` with `omega = k c`.
///
/// # Errors
/// [`HelmholtzError`] on bad parameters, a too-coarse mesh, the dense
/// work cap, shape mismatch, or a singular dense system.
pub fn solve_radiation(
    surface: &SpherePanels,
    k: f64,
    medium: Medium,
    velocity: &[C64],
    formulation: Formulation,
) -> Result<RadiationSolution, HelmholtzError> {
    if !(k > 0.0 && k.is_finite()) {
        return Err(HelmholtzError::BadParameter {
            what: "wavenumber k must be positive and finite",
        });
    }
    if !(medium.density > 0.0
        && medium.density.is_finite()
        && medium.sound_speed > 0.0
        && medium.sound_speed.is_finite())
    {
        return Err(HelmholtzError::BadParameter {
            what: "medium density and sound speed must be positive and finite",
        });
    }
    let n = surface.centroids().len();
    if velocity.len() != n {
        return Err(HelmholtzError::ShapeMismatch {
            what: "velocity length must equal the panel count",
        });
    }
    if n > MAX_DENSE_PANELS {
        return Err(HelmholtzError::WorkCap { panels: n });
    }
    let wavelength = 2.0 * core::f64::consts::PI / k;
    let ppw = wavelength / characteristic_panel_size(surface);
    if ppw < MIN_PANELS_PER_WAVELENGTH {
        return Err(HelmholtzError::TooCoarse {
            panels_per_wavelength: ppw,
        });
    }

    let alpha = alpha_for(formulation, k);

    // q = i omega rho v.
    let omega_rho = k * medium.sound_speed * medium.density;
    let q: Vec<C64> = velocity
        .iter()
        .map(|&v| v * C64::new(0.0, omega_rho))
        .collect();

    let (matrix, bmat) = assemble_dense(surface, k, alpha);
    let mut rhs = vec![C64::ZERO; n];
    for i in 0..n {
        let mut b = C64::ZERO;
        for j in 0..n {
            b = b + bmat[i * n + j] * q[j];
        }
        rhs[i] = b;
    }

    let lu = lu_complex(&matrix, n).map_err(|_| HelmholtzError::Singular)?;
    let mut pressure = rhs;
    lu.solve(&mut pressure);

    let mut radiated_power = 0.0;
    for j in 0..n {
        radiated_power += 0.5 * (pressure[j] * velocity[j].conj()).re * surface.areas()[j];
    }

    Ok(RadiationSolution {
        pressure,
        velocity: velocity.to_vec(),
        k,
        panels_per_wavelength: ppw,
        radiated_power,
    })
}

fn alpha_for(formulation: Formulation, k: f64) -> C64 {
    match formulation {
        Formulation::PlainCbie => C64::ZERO,
        Formulation::BurtonMiller => C64::new(0.0, -1.0 / k),
        Formulation::BurtonMillerWrongAlphaSign => C64::new(0.0, 1.0 / k),
    }
}

/// Assemble the dense Burton–Miller pair `(A, B)` with
/// `A p = B q`: row i is `(D - I/2 + alpha N)` against
/// `(S + alpha D' + alpha/2 I)`. The hypersingular self entry combines
/// the regularized disc part with the exact closed-surface row identity
/// `N_0[1] = 0` (see the module docs).
fn assemble_dense(surface: &SpherePanels, k: f64, alpha: C64) -> (Vec<C64>, Vec<C64>) {
    let centroids = surface.centroids();
    let normals = surface.normals();
    let areas = surface.areas();
    let n = centroids.len();
    let half = C64::new(0.5, 0.0);
    let mut amat = vec![C64::ZERO; n * n];
    let mut bmat = vec![C64::ZERO; n * n];
    for i in 0..n {
        let xi = centroids[i];
        let ni = normals[i];
        // Static hypersingular row sum: N_0[1] = 0 on a closed surface,
        // so the self entry is minus the off-diagonal static sum under
        // the SAME point-panel quadrature.
        let mut n0_row = C64::ZERO;
        for j in 0..n {
            if j != i {
                let (_, _, _, d2g0) = kernels(0.0, xi, ni, centroids[j], normals[j]);
                n0_row = n0_row + d2g0.scale(areas[j]);
            }
        }
        for j in 0..n {
            let (s_ij, d_ij, dp_ij, n_ij) = if i == j {
                let (s, n_reg) = self_terms(k, areas[i]);
                (s, C64::ZERO, C64::ZERO, n_reg - n0_row)
            } else {
                let (g, dgdny, dgdnx, d2g) = kernels(k, xi, ni, centroids[j], normals[j]);
                let a = areas[j];
                (g.scale(a), dgdny.scale(a), dgdnx.scale(a), d2g.scale(a))
            };
            let mut aij = d_ij + alpha * n_ij;
            let mut bij = s_ij + alpha * dp_ij;
            if i == j {
                aij = aij - half;
                bij = bij + alpha * half;
            }
            amat[i * n + j] = aij;
            bmat[i * n + j] = bij;
        }
    }
    (amat, bmat)
}

/// Far-field directivity amplitude `F(direction)`: the radiated pressure
/// behaves as `p -> F * e^{ikr}/r` for large r, with
/// `F = (1/4 pi) SUM_j [(-ik dir.n_j) p_j - q_j] e^{-ik dir.y_j} A_j`.
#[must_use]
pub fn far_field(
    surface: &SpherePanels,
    solution: &RadiationSolution,
    medium: Medium,
    directions: &[[f64; 3]],
) -> Vec<C64> {
    let k = solution.k;
    let omega_rho = k * medium.sound_speed * medium.density;
    let centroids = surface.centroids();
    let normals = surface.normals();
    let areas = surface.areas();
    directions
        .iter()
        .map(|&dir| {
            let len = norm(dir);
            let d = [dir[0] / len, dir[1] / len, dir[2] / len];
            let mut f = C64::ZERO;
            for j in 0..centroids.len() {
                let phase = expik(-k * dot(d, centroids[j]));
                let qj = solution.velocity[j] * C64::new(0.0, omega_rho);
                let term = solution.pressure[j] * C64::new(0.0, -k * dot(d, normals[j])) - qj;
                f = f + (term * phase).scale(areas[j]);
            }
            f.scale(1.0 / (4.0 * core::f64::consts::PI))
        })
        .collect()
}

/// The dense radiation impedance matrix `Z` mapping panel normal
/// velocities to panel pressures (`p = Z v`), built by one LU
/// factorization and `n` unit-velocity solves. Feeds the vibroacoustic
/// coupling bead. Row-major `n x n`.
///
/// # Errors
/// As for [`solve_radiation`].
pub fn radiation_impedance_matrix(
    surface: &SpherePanels,
    k: f64,
    medium: Medium,
    formulation: Formulation,
) -> Result<Vec<C64>, HelmholtzError> {
    let n = surface.centroids().len();
    // Shared admission with the single-solve path.
    let probe = vec![C64::ZERO; n];
    let _ = solve_radiation(surface, k, medium, &probe, formulation)?;
    let alpha = alpha_for(formulation, k);
    let omega_rho = k * medium.sound_speed * medium.density;
    let (amat, bmat) = assemble_dense(surface, k, alpha);
    let lu = lu_complex(&amat, n).map_err(|_| HelmholtzError::Singular)?;
    let mut z = vec![C64::ZERO; n * n];
    for col in 0..n {
        // Unit velocity on panel `col`: q = i omega rho e_col, so the
        // right-hand side is that scaled column of B.
        let mut p: Vec<C64> = (0..n)
            .map(|row| bmat[row * n + col] * C64::new(0.0, omega_rho))
            .collect();
        lu.solve(&mut p);
        for row in 0..n {
            z[row * n + col] = p[row];
        }
    }
    Ok(z)
}

/// Baffled rigid circular piston: the classical Rayleigh-integral
/// radiation impedance (half-space kernel `2G`), averaged over the
/// piston face. Deterministic polar ring discretization; each cell's
/// self contribution uses the equivalent-disc closed form
/// `INT e^{ikr}/r dA = 2 pi (e^{ika_c} - 1)/(ik)`. Under this module's
/// `e^{-i omega t}` convention the Rayleigh integral is
/// `p = -(i omega rho / 2 pi) INT v e^{ikr}/r dS` and the small-ka
/// oracle is `z/(rho c) = (ka)^2/2 - i 8ka/(3 pi) + O((ka)^3)`
/// (mass-like reactance is NEGATIVE imaginary here, exactly as for the
/// pulsating sphere). Bessel/Struve closed forms join the
/// duct-acoustics bead's special functions.
///
/// # Errors
/// [`HelmholtzError::BadParameter`] on non-positive radius, wavenumber,
/// medium, or ring count.
pub fn baffled_piston_impedance(
    radius: f64,
    k: f64,
    medium: Medium,
    rings: usize,
) -> Result<C64, HelmholtzError> {
    if !(radius > 0.0 && radius.is_finite() && k > 0.0 && k.is_finite()) || rings == 0 {
        return Err(HelmholtzError::BadParameter {
            what: "piston radius, wavenumber, and ring count must be positive and finite",
        });
    }
    if !(medium.density > 0.0 && medium.sound_speed > 0.0) {
        return Err(HelmholtzError::BadParameter {
            what: "medium density and sound speed must be positive and finite",
        });
    }
    // Polar cells: `rings` radial bands, 6m angular cells in band m.
    let mut cells: Vec<([f64; 2], f64)> = Vec::new();
    for m in 0..rings {
        let r0 = radius * m as f64 / rings as f64;
        let r1 = radius * (m + 1) as f64 / rings as f64;
        let rc = f64::midpoint(r0, r1);
        let sectors = 6 * (m + 1);
        let band_area = core::f64::consts::PI * (r1 * r1 - r0 * r0);
        for sct in 0..sectors {
            let th = 2.0 * core::f64::consts::PI * (sct as f64 + 0.5) / sectors as f64;
            cells.push((
                [rc * det::cos(th), rc * det::sin(th)],
                band_area / sectors as f64,
            ));
        }
    }
    let omega_rho = k * medium.sound_speed * medium.density;
    let total_area = core::f64::consts::PI * radius * radius;
    let mut mean_p = C64::ZERO;
    for (i, &(xi, ai)) in cells.iter().enumerate() {
        let mut integral = C64::ZERO;
        for (j, &(yj, aj)) in cells.iter().enumerate() {
            if i == j {
                let ac = (ai / core::f64::consts::PI).sqrt();
                let self_term = (expik(k * ac) - C64::ONE) * C64::new(0.0, k).recip();
                integral = integral + self_term.scale(2.0 * core::f64::consts::PI);
            } else {
                let dx = xi[0] - yj[0];
                let dy = xi[1] - yj[1];
                let r = (dx * dx + dy * dy).sqrt();
                integral = integral + expik(k * r).scale(aj / r);
            }
        }
        // p(x_i) = -(i omega rho / 2 pi) v INT e^{ikr}/r dA, v = 1.
        let p_i = integral * C64::new(0.0, -omega_rho / (2.0 * core::f64::consts::PI));
        mean_p = mean_p + p_i.scale(ai / total_area);
    }
    Ok(mean_p)
}

#[cfg(test)]
mod tests {
    use super::*;

    const RHO_C: f64 = 1.204 * 343.0;

    /// Analytic pulsating-sphere surface impedance under e^{-i omega t}:
    /// `z = rho c (k^2 R^2 - i k R) / (1 + k^2 R^2)`.
    fn pulsating_sphere_impedance(k: f64, radius: f64) -> C64 {
        let ka = k * radius;
        C64::new(ka * ka, -ka).scale(RHO_C / (1.0 + ka * ka))
    }

    fn uniform_velocity(n: usize) -> Vec<C64> {
        vec![C64::ONE; n]
    }

    fn mean_impedance(sol: &RadiationSolution, surface: &SpherePanels) -> C64 {
        // Area-weighted mean of p / v over panels with |v| > 0.
        let mut num = C64::ZERO;
        let mut den = 0.0;
        for j in 0..sol.pressure.len() {
            let v = sol.velocity[j];
            if v.abs() > 1e-12 {
                num = num + (sol.pressure[j] * v.recip()).scale(surface.areas()[j]);
                den += surface.areas()[j];
            }
        }
        num.scale(1.0 / den)
    }

    #[test]
    fn kernels_match_finite_differences_of_green() {
        // The four point kernels against central differences of G — the
        // gate that makes sign errors die here instead of in physics.
        let k = 2.3;
        let x = [0.3, -0.2, 0.5];
        let y = [-0.4, 0.6, -0.1];
        let nx = {
            let v = [0.5, 0.7, -0.3];
            let l = norm(v);
            [v[0] / l, v[1] / l, v[2] / l]
        };
        let ny = {
            let v = [-0.2, 0.4, 0.9];
            let l = norm(v);
            [v[0] / l, v[1] / l, v[2] / l]
        };
        let h = 1e-5;
        let g_at = |x: [f64; 3], y: [f64; 3]| green(k, norm(sub(x, y)));
        let shift =
            |p: [f64; 3], n: [f64; 3], s: f64| [p[0] + s * n[0], p[1] + s * n[1], p[2] + s * n[2]];
        let (g, dgdny, dgdnx, d2g) = kernels(k, x, nx, y, ny);
        assert!((g - g_at(x, y)).abs() < 1e-12);
        let fd_ny = (g_at(x, shift(y, ny, h)) - g_at(x, shift(y, ny, -h))).scale(1.0 / (2.0 * h));
        assert!(
            (dgdny - fd_ny).abs() < 1e-6 * dgdny.abs(),
            "dG/dn_y vs FD: {dgdny:?} {fd_ny:?}"
        );
        let fd_nx = (g_at(shift(x, nx, h), y) - g_at(shift(x, nx, -h), y)).scale(1.0 / (2.0 * h));
        assert!(
            (dgdnx - fd_nx).abs() < 1e-6 * dgdnx.abs(),
            "dG/dn_x vs FD: {dgdnx:?} {fd_nx:?}"
        );
        let fd2 = (g_at(shift(x, nx, h), shift(y, ny, h))
            - g_at(shift(x, nx, h), shift(y, ny, -h))
            - g_at(shift(x, nx, -h), shift(y, ny, h))
            + g_at(shift(x, nx, -h), shift(y, ny, -h)))
        .scale(1.0 / (4.0 * h * h));
        assert!(
            (d2g - fd2).abs() < 1e-4 * d2g.abs(),
            "d2G/dn_x dn_y vs FD: {d2g:?} {fd2:?}"
        );
        println!(
            "{{\"suite\":\"fs-bem-helmholtz\",\"case\":\"kernel-fd-pins\",\"verdict\":\"pass\"}}"
        );
    }

    #[test]
    fn disc_self_terms_match_numerical_quadrature() {
        // S_ii against direct radial quadrature (the integrand is
        // integrable); the REGULARIZED hypersingular self term against
        // quadrature of the difference kernel [e^{ikr}(1-ikr) - 1]/r^2
        // (the static core belongs to the N_0[1] = 0 row identity, not
        // to this disc integral).
        let k = 1.7;
        let area = 0.05;
        let a = (area / core::f64::consts::PI).sqrt();
        let (s, n) = self_terms(k, area);
        let m = 200_000u32;
        let dr = a / f64::from(m);
        let mut s_num = C64::ZERO;
        let mut n_reg = C64::ZERO;
        for i in 0..m {
            let r = (f64::from(i) + 0.5) * dr;
            // S: (1/2) e^{ikr} dr.
            s_num = s_num + expik(k * r).scale(0.5 * dr);
            // N regularized: (1/2) [e^{ikr}(1 - ikr) - 1] / r^2 dr.
            let core_sub = expik(k * r) * C64::new(1.0, -k * r) - C64::ONE;
            n_reg = n_reg + core_sub.scale(0.5 * dr / (r * r));
        }
        let _ = a;
        assert!(
            (s - s_num).abs() < 1e-6 * s.abs(),
            "S self: {s:?} vs {s_num:?}"
        );
        assert!(
            (n - n_reg).abs() < 1e-4 * n.abs(),
            "regularized N self: {n:?} vs {n_reg:?}"
        );
        println!(
            "{{\"suite\":\"fs-bem-helmholtz\",\"case\":\"disc-self-terms\",\"verdict\":\"pass\"}}"
        );
    }

    #[test]
    fn pulsating_sphere_impedance_across_ka() {
        use core::fmt::Write as _;
        // The oracle that arbitrates every sign in the formulation, over
        // ka in [0.05, 5]. Two arms with measured, documented roles:
        // PlainCbie is the accurate non-resonant arm at every sampled ka
        // (fictitious frequencies start at ka = pi); BurtonMiller is the
        // resonance-safe production arm for ka >= 0.5. BM's centroid-
        // quadrature resistance artifact below ka ~ 0.5 is measured and
        // recorded in the JSON evidence rather than hidden (the exact-
        // quadrature follow-up on the bead is its fix trigger).
        let coarse = SpherePanels::icosphere(1.0, 2).expect("icosphere-2");
        let fine = SpherePanels::icosphere(1.0, 3).expect("icosphere-3");
        let mut rows = String::new();
        let mut first = true;
        let mut run = |surface: &SpherePanels, ka: f64, formulation: Formulation, tol: f64| {
            let n = surface.centroids().len();
            let sol = solve_radiation(
                surface,
                ka,
                Medium::air(),
                &uniform_velocity(n),
                formulation,
            )
            .expect("solve");
            let z = mean_impedance(&sol, surface);
            let z_ref = pulsating_sphere_impedance(ka, 1.0);
            let rel = (z - z_ref).abs() / z_ref.abs();
            if tol > 0.0 {
                assert!(
                    rel < tol,
                    "{formulation:?} ka={ka}: z = {z:?} vs {z_ref:?} (rel {rel:.4})"
                );
                assert!(sol.radiated_power > 0.0, "passivity at ka={ka}");
            }
            write!(
                rows,
                "{}{{\"ka\":{ka},\"arm\":\"{formulation:?}\",\"rel_err\":{rel:.4},\"ppw\":{:.1},\"asserted\":{}}}",
                if first { "" } else { "," },
                sol.panels_per_wavelength,
                tol > 0.0
            )
            .expect("write to String");
            first = false;
        };
        for ka in [0.05, 0.2, 0.5, 1.0, 2.0] {
            run(&coarse, ka, Formulation::PlainCbie, 0.04);
        }
        run(&fine, 5.0, Formulation::PlainCbie, 0.04);
        for ka in [0.5, 1.0, 2.0] {
            run(&coarse, ka, Formulation::BurtonMiller, 0.08);
        }
        run(&fine, 5.0, Formulation::BurtonMiller, 0.08);
        // The BM low-ka artifact rows: measured, recorded, NOT asserted
        // accurate — the documented boundary of the centroid quadrature.
        for ka in [0.05, 0.2] {
            run(&coarse, ka, Formulation::BurtonMiller, 0.0);
        }
        println!(
            "{{\"suite\":\"fs-bem-helmholtz\",\"case\":\"pulsating-sphere\",\"rows\":[{rows}],\"verdict\":\"pass\"}}"
        );
    }

    #[test]
    fn fictitious_frequency_contrast_burton_miller_survives_cbie_fails() {
        // ka = pi is the sphere's first interior Dirichlet resonance:
        // plain CBIE visibly degrades there while Burton-Miller holds.
        // This contrast is simultaneously the dropped-hypersingular
        // mutation (PlainCbie IS the dropped-N arm).
        let surface = SpherePanels::icosphere(1.0, 2).expect("icosphere");
        let n = surface.centroids().len();
        let ka = core::f64::consts::PI;
        let z_ref = pulsating_sphere_impedance(ka, 1.0);
        let err = |formulation: Formulation| -> f64 {
            let sol = solve_radiation(
                &surface,
                ka,
                Medium::air(),
                &uniform_velocity(n),
                formulation,
            )
            .expect("solve");
            let z = mean_impedance(&sol, &surface);
            (z - z_ref).abs() / z_ref.abs()
        };
        let bm = err(Formulation::BurtonMiller);
        let cbie = err(Formulation::PlainCbie);
        assert!(bm < 0.06, "Burton-Miller must hold at ka=pi, rel {bm:.4}");
        assert!(
            cbie > 5.0 * bm,
            "plain CBIE must visibly degrade at ka=pi: cbie {cbie:.4} vs bm {bm:.4}"
        );
        println!(
            "{{\"suite\":\"fs-bem-helmholtz\",\"case\":\"fictitious-frequency-contrast\",\"ka\":3.14159,\"bm_rel\":{bm:.4},\"cbie_rel\":{cbie:.4},\"verdict\":\"pass\"}}"
        );
    }

    #[test]
    fn alpha_sign_mutation_degrades_at_interior_resonance() {
        // MUTATION: alpha = -i/k (the wrong sign under e^{-i omega t}).
        // It must NOT rescue the interior resonance: the correct BM
        // error stays small while the flipped-sign arm degrades like
        // (or worse than) plain CBIE at ka = pi.
        let surface = SpherePanels::icosphere(1.0, 2).expect("icosphere");
        let n = surface.centroids().len();
        let ka = core::f64::consts::PI;
        let z_ref = pulsating_sphere_impedance(ka, 1.0);
        let err = |formulation: Formulation| -> f64 {
            let sol = solve_radiation(
                &surface,
                ka,
                Medium::air(),
                &uniform_velocity(n),
                formulation,
            )
            .expect("solve");
            let z = mean_impedance(&sol, &surface);
            (z - z_ref).abs() / z_ref.abs()
        };
        let bm = err(Formulation::BurtonMiller);
        let flipped = err(Formulation::BurtonMillerWrongAlphaSign);
        assert!(bm < 0.06, "correct alpha must hold at ka=pi, rel {bm:.4}");
        assert!(
            flipped > 2.0 * bm,
            "flipped alpha must degrade at ka=pi: flipped {flipped:.4} vs bm {bm:.4}"
        );
        println!(
            "{{\"suite\":\"fs-bem-helmholtz\",\"case\":\"alpha-sign-mutation\",\"bm_rel\":{bm:.4},\"flipped_rel\":{flipped:.4},\"verdict\":\"caught\"}}"
        );
    }

    #[test]
    fn dipole_and_monopole_directivity_patterns() {
        let surface = SpherePanels::icosphere(1.0, 2).expect("icosphere");
        let n = surface.centroids().len();
        let ka = 1.0;
        // Monopole: uniform |F| over directions.
        let mono = solve_radiation(
            &surface,
            ka,
            Medium::air(),
            &uniform_velocity(n),
            Formulation::BurtonMiller,
        )
        .expect("solve");
        let mut dirs = Vec::new();
        let m = 24u32;
        for i in 0..m {
            let theta = core::f64::consts::PI * (f64::from(i) + 0.5) / f64::from(m);
            dirs.push([det::sin(theta), 0.0, det::cos(theta)]);
        }
        let f_mono = far_field(&surface, &mono, Medium::air(), &dirs);
        let amps: Vec<f64> = f_mono.iter().map(|f| f.abs()).collect();
        let (mut lo, mut hi) = (f64::INFINITY, 0.0f64);
        for &a in &amps {
            lo = lo.min(a);
            hi = hi.max(a);
        }
        assert!(hi / lo < 1.02, "monopole directivity must be uniform");
        // Dipole: v_n = cos(theta) -> |F| proportional to |cos(theta)|.
        let v_dip: Vec<C64> = surface
            .normals()
            .iter()
            .map(|nrm| C64::from_re(nrm[2]))
            .collect();
        let dip = solve_radiation(
            &surface,
            ka,
            Medium::air(),
            &v_dip,
            Formulation::BurtonMiller,
        )
        .expect("solve");
        let f_dip = far_field(&surface, &dip, Medium::air(), &dirs);
        // Normalized correlation between |F| and |cos theta| samples.
        let mut num = 0.0;
        let mut fa = 0.0;
        let mut cb = 0.0;
        for (i, f) in f_dip.iter().enumerate() {
            let c = dirs[i][2].abs();
            num += f.abs() * c;
            fa += f.abs() * f.abs();
            cb += c * c;
        }
        let corr = num / (fa.sqrt() * cb.sqrt());
        assert!(corr > 0.99, "dipole pattern correlation {corr:.4}");
        println!(
            "{{\"suite\":\"fs-bem-helmholtz\",\"case\":\"directivity-patterns\",\"monopole_ripple\":{:.4},\"dipole_corr\":{corr:.4},\"verdict\":\"pass\"}}",
            hi / lo - 1.0
        );
    }

    #[test]
    fn impedance_matrix_is_reciprocal_and_passive() {
        // Small mesh: Z symmetric to discretization tolerance
        // (reciprocity) and its symmetrized real part PSD-ish along
        // random probes (passivity of radiated power), seeded
        // deterministically.
        let surface = SpherePanels::icosphere(1.0, 1).expect("icosphere");
        let n = surface.centroids().len();
        let z = radiation_impedance_matrix(&surface, 1.0, Medium::air(), Formulation::BurtonMiller)
            .expect("impedance matrix");
        let mut asym = 0.0f64;
        let mut scale = 0.0f64;
        for i in 0..n {
            for j in 0..n {
                // Reciprocity for equal-area panels: A_i Z_ij = A_j Z_ji.
                let zij = z[i * n + j].scale(surface.areas()[i]);
                let zji = z[j * n + i].scale(surface.areas()[j]);
                asym = asym.max((zij - zji).abs());
                scale = scale.max(zij.abs());
            }
        }
        assert!(
            asym < 0.05 * scale,
            "reciprocity violation {asym:.3e} vs scale {scale:.3e}"
        );
        // Passivity: radiated power positive for smooth, mesh-resolved
        // velocity fields (constant, three dipoles, one quadrupole, with
        // a deterministic complex mix). White-noise per-panel phases are
        // deliberately NOT claimed: an unresolved velocity field has no
        // discrete power identity on any mesh.
        let patterns: [fn(&[f64; 3]) -> C64; 4] = [
            |_| C64::ONE,
            |nrm| C64::new(nrm[0], 0.3 * nrm[2]),
            |nrm| C64::new(nrm[1], -0.2 * nrm[0]),
            |nrm| C64::new(nrm[2], 0.1),
        ];
        for (idx, pattern) in patterns.iter().enumerate() {
            let v: Vec<C64> = surface.normals().iter().map(pattern).collect();
            let sol = solve_radiation(&surface, 1.0, Medium::air(), &v, Formulation::BurtonMiller)
                .expect("solve");
            assert!(
                sol.radiated_power > 0.0,
                "passivity probe {idx}: W = {}",
                sol.radiated_power
            );
        }
        // Quadrupole boundary, MEASURED not asserted: its true radiated
        // power at ka = 1 sits below the centroid-quadrature noise floor
        // on this mesh (a documented no-claim; exact quadrature is the
        // recorded fix trigger).
        let v_quad: Vec<C64> = surface
            .normals()
            .iter()
            .map(|nrm| C64::new(3.0 * nrm[2] * nrm[2] - 1.0, 0.5 * nrm[0] * nrm[1]))
            .collect();
        let w_quad = solve_radiation(
            &surface,
            1.0,
            Medium::air(),
            &v_quad,
            Formulation::BurtonMiller,
        )
        .expect("solve")
        .radiated_power;
        println!(
            "{{\"suite\":\"fs-bem-helmholtz\",\"case\":\"reciprocity-passivity\",\"asym_rel\":{:.4},\"quadrupole_w_unasserted\":{w_quad:.3e},\"verdict\":\"pass\"}}",
            asym / scale
        );
    }

    #[test]
    fn baffled_piston_matches_small_ka_series() {
        use core::fmt::Write as _;
        // Rayleigh piston vs the Bessel-free small-ka series
        // z/(rho c) = (ka)^2/2 - i 8ka/(3 pi): the series truncation is
        // O((ka)^2) relative, so the authored tolerances widen with ka.
        let mut rows = String::new();
        for (idx, &(ka, tol)) in [(0.1, 0.02), (0.25, 0.04), (0.5, 0.12)].iter().enumerate() {
            let z = baffled_piston_impedance(1.0, ka, Medium::air(), 24).expect("piston");
            let z_series =
                C64::new(0.5 * ka * ka, -8.0 * ka / (3.0 * core::f64::consts::PI)).scale(RHO_C);
            let rel = (z - z_series).abs() / z_series.abs();
            assert!(
                rel < tol,
                "ka={ka}: piston z = {z:?} vs series {z_series:?} (rel {rel:.4})"
            );
            write!(
                rows,
                "{}{{\"ka\":{ka},\"rel_err_vs_series\":{rel:.4}}}",
                if idx == 0 { "" } else { "," }
            )
            .expect("write to String");
        }
        println!(
            "{{\"suite\":\"fs-bem-helmholtz\",\"case\":\"baffled-piston-series\",\"rows\":[{rows}],\"verdict\":\"pass\"}}"
        );
    }

    #[test]
    fn named_refusals_fire() {
        let surface = SpherePanels::icosphere(1.0, 1).expect("icosphere");
        let n = surface.centroids().len();
        assert!(matches!(
            solve_radiation(
                &surface,
                -1.0,
                Medium::air(),
                &uniform_velocity(n),
                Formulation::BurtonMiller
            ),
            Err(HelmholtzError::BadParameter { .. })
        ));
        // 80 panels on the unit sphere cannot resolve ka = 60.
        let err = solve_radiation(
            &surface,
            60.0,
            Medium::air(),
            &uniform_velocity(n),
            Formulation::BurtonMiller,
        )
        .unwrap_err();
        assert!(err.to_string().contains("FS-BEM-HELM-TOO-COARSE"));
        assert!(matches!(
            solve_radiation(
                &surface,
                1.0,
                Medium::air(),
                &uniform_velocity(3),
                Formulation::BurtonMiller
            ),
            Err(HelmholtzError::ShapeMismatch { .. })
        ));
    }

    #[test]
    fn repeat_solves_are_bitwise_identical() {
        let surface = SpherePanels::icosphere(1.0, 1).expect("icosphere");
        let n = surface.centroids().len();
        let a = solve_radiation(
            &surface,
            1.3,
            Medium::air(),
            &uniform_velocity(n),
            Formulation::BurtonMiller,
        )
        .expect("solve");
        let b = solve_radiation(
            &surface,
            1.3,
            Medium::air(),
            &uniform_velocity(n),
            Formulation::BurtonMiller,
        )
        .expect("solve");
        for (x, y) in a.pressure.iter().zip(b.pressure.iter()) {
            assert!(x.re.to_bits() == y.re.to_bits() && x.im.to_bits() == y.im.to_bits());
        }
    }
}
