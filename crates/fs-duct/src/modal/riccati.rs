//! Continuous multimodal (m = 0) impedance propagation on smooth radius
//! profiles (IRP M2.1-W3b), after Pagneux, Amir & Kergomard (1996), JASA
//! 100, 2034-2048, Part I, and Amir, Pagneux & Kergomard (1997), JASA 101,
//! 2504-2517. Unlike the staircase in the parent module, no radius step is
//! ever introduced; the staircase remains the O(dx) comparison rung.
//!
//! **Modes and amplitudes** are the staircase's: `psi_n(r; R) =
//! J0(gamma_n r/R)/|J0(gamma_n)|` (orthonormal under `(1/S) int dS`),
//! `P_n = (1/S) int p psi_n dS`, `U_n = int v_x psi_n dS`, so `U_0` is the
//! volume flow and `P = Z U`.
//!
//! **Coupled equations.** Under `e^{-i omega t}` (`grad p = i omega rho v`,
//! `div v = i omega p/(rho c^2)`, rigid sloped wall `v_r = R' v_x` at
//! `r = R`), projecting onto `psi_m` gives, with `s_n = psi_n(R) =
//! sgn J0(gamma_n)` (`s_0 = 1`, `gamma_0 = 0`):
//!
//! ```text
//! P' = (R'/R) A P + Z' U,         Z'  = i k_0 Zc_0           (= i omega rho/S lossless)
//! U' = Y' P - (R'/R) A^T U,       Y'_n = i k_n^2 / (k_0 Zc_0),  k_n^2 = k_0^2 - (gamma_n/R)^2
//! A_mn = 2 gamma_n^2 s_m s_n / (gamma_n^2 - gamma_m^2)  (m != n),   A_nn = 0
//! ```
//!
//! Derivation of `A`, in `eta = r/R`: differentiating `(1/S) int p psi_m dS`
//! leaves the wall term `2 s_m s_n` and the moving-basis term
//! `-2 int_0^1 eta^2 psi_n psi_m' d eta`. The latter integral equals
//! `gamma_m^2 s_m s_n / (gamma_m^2 - gamma_n^2)` for `m != n` and 0 for
//! `m = n`: integrate `d/d eta [eta^2 psi_m' psi_n']` with the Bessel ODE and
//! `psi'(1) = 0`, and `int_0^1 (eta^2 psi_m psi_n)' d eta = s_m s_n` with
//! orthogonality. In the flow equation the boundary terms cancel; the
//! transverse divergence gives `-(gamma_m/R)^2` (Green's identity, exact
//! even though `dp/dr != 0` at a sloped wall) and the moving basis gives
//! the same integral with its indices exchanged, `-A^T`.
//!
//! Consequences used as tests: `N = 1` is Webster's equation exactly; with
//! `R' = 0` the modes decouple into the staircase's uniform waves; and the
//! generator `[[C, D], [E, -C^T]]` with symmetric `D`, `E` is Hamiltonian,
//! so the continuous input impedance is reciprocal by construction.
//!
//! **Losses.** The plane mode takes `segment_wave` verbatim. Higher modes
//! follow the staircase's disclosed rule (`k_n` from the lossy `k_0`,
//! `Zc_n = Zc_0 k_0/k_n`); a uniform section of this system reproduces both
//! exactly. Per-mode boundary-layer losses belong to M2.1-W6.
//!
//! **Integration.** From the mouth to the inlet in fixed steps. Each step is
//! the fourth-order two-node Gauss-Magnus exponential of the `2N x 2N`
//! generator (scaled by `Z_ref = rho c / S_inlet`), applied as a Moebius map
//! to the scaled impedance `zeta = Z/Z_ref`:
//! `zeta <- (Phi11 zeta + Phi12)(Phi21 zeta + Phi22)^-1`. The exact local
//! exponential carries evanescent modes without explicit Riccati marching:
//! integrated toward the source the impedance relaxes onto the decaying
//! branch (local rate `-2 kappa_n`), and only ratios enter the map.
//!
//! **Mouth.** Closed and ideal pressure-release boundaries act on every
//! retained mode. Scalar radiation and tabulated loads act on the plane mode,
//! with the disclosed matched higher modes. [`MouthLoad::Matched`] is the
//! infinite-pipe continuation of every mode. The explicit `mm_profile_*_flanged`
//! APIs impose the complete original Zorumski radiation matrix, retaining
//! signed mutual loads and the exterior gas.
//!
//! **Field.** [`mm_profile_field`] keeps the scaled impedance at every step
//! station and each step's inverse flow map, then carries a prescribed
//! inlet modal flow back toward the mouth (`u_j = Y_j^-1 u_{j+1}`, which
//! contracts evanescent components: they decay away from the source) and
//! sets `P = Z U` station by station. The pressure on the axis is
//! `sum_n P_n psi_n(0)`, `psi_n(0) = 1/|J0(gamma_n)|`.
//!
//! **Bounds.** Mode counts are admitted up to [`certified_mode_cap`]: every
//! `gamma_n` used lies below [`BESSEL_CERTIFIED_X`], the range over which the
//! tests certify `bessel_j0_certified` and `J1` against an independent
//! quadrature. (The staircase's legacy `bessel_j0` is not used here.)

use super::{ModalResponse, bessel_j0_certified, j1_roots_unchecked, mat_inv, modal_wave};
use crate::profile::{RadiusProfile, RiccatiSteps};
use crate::{
    DuctError, LossModel, TabulatedLoad, Termination, segment_wave, termination_impedance,
};
use fs_material::gas::GasState;
use fs_math::c64::C64;
use fs_math::det;

/// Upper end of the Bessel range certified for continuous propagation.
pub const BESSEL_CERTIFIED_X: f64 = 64.0;
/// Largest admitted evanescent decay per step, `max_n |Im k_n| |h|`. The step
/// exponential spans a dynamic range of about `e^(kappa |h|)`, and scaling
/// and squaring errors are relative to its norm, so once modes couple the
/// O(1) propagating entries keep roughly `eps e^(2 kappa |h|)` absolute
/// accuracy. A larger step refuses ("refine Riccati steps") instead.
pub const MAX_EVANESCENT_STEP: f64 = 8.0;
/// Squarings allowed in one step exponential before the step must be refined.
const MAX_SQUARINGS: u32 = 40;
/// Taylor terms allowed for the scaled exponential (converges by ~20).
const MAX_TAYLOR_TERMS: u32 = 40;

fn bad(what: &'static str) -> DuctError {
    DuctError::BadParameter { what }
}
fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}

/// Largest admitted mode count: the plane mode plus every `J1` root below
/// [`BESSEL_CERTIFIED_X`]. Computed from the roots, never a hand constant.
#[must_use]
pub fn certified_mode_cap() -> usize {
    // gamma_n ~ (n + 1/4) pi; extend until the last root clears the range.
    let mut count = (BESSEL_CERTIFIED_X / core::f64::consts::PI) as usize + 1;
    loop {
        let roots = j1_roots_unchecked(count);
        if roots.last().is_some_and(|&g| g >= BESSEL_CERTIFIED_X) {
            return 1 + roots.iter().filter(|&&g| g < BESSEL_CERTIFIED_X).count();
        }
        count += 1;
    }
}

/// Retained mode count, including the plane mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ModeCount(pub usize);

impl ModeCount {
    /// Admit `1..=certified_mode_cap()`; anything else refuses by name.
    ///
    /// # Errors
    /// [`DuctError::BadParameter`].
    pub fn validate(self) -> Result<(), DuctError> {
        if self.0 == 0 || self.0 > certified_mode_cap() {
            return Err(bad(
                "mode count outside 1..=certified_mode_cap() (certified Bessel range)",
            ));
        }
        Ok(())
    }
}

/// Mouth load: exact closed/pressure-release boundaries act on all modes;
/// scalar radiation and tabulated loads retain matched higher modes.
#[derive(Debug, Clone, Copy)]
pub enum MouthLoad<'a> {
    /// Exact ideal boundaries on all modes, or scalar radiation with its refusals.
    Analytic(Termination),
    /// A tabulated bake; out-of-table frequencies refuse.
    Tabulated(&'a TabulatedLoad),
    /// Every mode into its own characteristic impedance: the semi-infinite
    /// continuation of the mouth cross-section (Pagneux's infinite-pipe end).
    Matched,
}

#[derive(Clone, Copy)]
enum ProfileLoad<'a> {
    Scalar(MouthLoad<'a>),
    Flanged(&'a GasState),
}

impl<'a> From<MouthLoad<'a>> for ProfileLoad<'a> {
    fn from(load: MouthLoad<'a>) -> Self {
        Self::Scalar(load)
    }
}

/// Row-major complex `n x n` product with fixed summation order.
fn mul(a: &[C64], b: &[C64], n: usize) -> Vec<C64> {
    let mut out = vec![C64::ZERO; n * n];
    for i in 0..n {
        for l in 0..n {
            let ail = a[i * n + l];
            if ail.re == 0.0 && ail.im == 0.0 {
                continue;
            }
            for j in 0..n {
                out[i * n + j] = out[i * n + j] + ail * b[l * n + j];
            }
        }
    }
    out
}

fn norm1(a: &[C64], n: usize) -> f64 {
    (0..n)
        .map(|j| (0..n).map(|i| a[i * n + j].abs()).sum::<f64>())
        .fold(0.0, f64::max)
}

/// `exp(a)` by scaling and squaring with a Taylor core at `||a/2^s||_1 <= 1/2`.
fn expm(a: &[C64], n: usize) -> Result<Vec<C64>, DuctError> {
    let norm = norm1(a, n);
    if !norm.is_finite() {
        return Err(DuctError::Singular);
    }
    let mut squarings = 0_u32;
    let mut scale = 1.0;
    while norm * scale > 0.5 {
        scale *= 0.5;
        squarings += 1;
        if squarings > MAX_SQUARINGS {
            return Err(bad("step exponential norm too large; refine Riccati steps"));
        }
    }
    let b: Vec<C64> = a.iter().map(|v| v.scale(scale)).collect();
    let mut sum = vec![C64::ZERO; n * n];
    for i in 0..n {
        sum[i * n + i] = C64::ONE;
    }
    let mut term = sum.clone();
    let mut converged = false;
    for k in 1..=MAX_TAYLOR_TERMS {
        term = mul(&term, &b, n);
        let inv = 1.0 / f64::from(k);
        for v in &mut term {
            *v = v.scale(inv);
        }
        for (s, t) in sum.iter_mut().zip(&term) {
            *s = *s + *t;
        }
        if norm1(&term, n) <= f64::EPSILON * 0.25 * norm1(&sum, n) {
            converged = true;
            break;
        }
    }
    if !converged {
        return Err(bad("step exponential Taylor series did not converge"));
    }
    for _ in 0..squarings {
        sum = mul(&sum, &sum, n);
    }
    if sum.iter().any(|&v| !finite(v)) {
        return Err(DuctError::Singular);
    }
    Ok(sum)
}

/// Mode data shared by every generator evaluation of one propagation.
pub(super) struct Basis {
    n: usize,
    pub(super) gammas: Vec<f64>,
    /// Real coupling matrix `A`, row-major `n x n`.
    pub(super) coupling: Vec<f64>,
}

impl Basis {
    pub(super) fn new(n: usize) -> Self {
        let mut gammas = vec![0.0];
        gammas.extend(j1_roots_unchecked(n - 1));
        let signs: Vec<f64> = gammas
            .iter()
            .enumerate()
            .map(|(i, &g)| {
                if i == 0 || bessel_j0_certified(g) > 0.0 {
                    1.0
                } else {
                    -1.0
                }
            })
            .collect();
        let mut coupling = vec![0.0; n * n];
        for m in 0..n {
            for k in 0..n {
                if m != k {
                    let (gm2, gk2) = (gammas[m] * gammas[m], gammas[k] * gammas[k]);
                    coupling[m * n + k] = 2.0 * gk2 * signs[m] * signs[k] / (gk2 - gm2);
                }
            }
        }
        Self {
            n,
            gammas,
            coupling,
        }
    }
}

/// Local generator data at one quadrature node.
struct Local {
    /// Scaled generator `M(x)` on `[P; Z_ref U]`, row-major `2n x 2n`.
    m: Vec<C64>,
    shear: f64,
    radius: f64,
    /// Largest modal decay rate `max_n Im k_n` [1/m].
    kappa: f64,
}

fn generator(
    profile: &RadiusProfile,
    state: &GasState,
    omega: f64,
    loss: LossModel,
    basis: &Basis,
    zref: f64,
    x: f64,
) -> Result<Local, DuctError> {
    let n = basis.n;
    let w = 2 * n;
    let radius = profile.radius(x)?;
    let flare = profile.slope(x)? / radius;
    let wave = segment_wave(state, radius, omega, loss)?;
    let (k0, zc0) = (wave.wavenumber, wave.characteristic_impedance);
    let k0zc0 = k0 * zc0;
    if k0zc0.abs() == 0.0 || !finite(k0zc0) {
        return Err(DuctError::Singular);
    }
    let i = C64::new(0.0, 1.0);
    let series = (i * k0zc0).scale(1.0 / zref);
    let shunt_scale = (i * k0zc0.recip()).scale(zref);
    let k0sq = k0 * k0;
    let mut m = vec![C64::ZERO; w * w];
    let mut kappa = 0.0_f64;
    for r in 0..n {
        for c in 0..n {
            let a = flare * basis.coupling[r * n + c];
            m[r * w + c] = C64::from_re(a);
            m[(n + c) * w + (n + r)] = C64::from_re(-a);
        }
        m[r * w + n + r] = series;
        let kt = basis.gammas[r] / radius;
        let kn2 = k0sq - C64::from_re(kt * kt);
        m[(n + r) * w + r] = shunt_scale * kn2;
        kappa = kappa.max(kn2.sqrt().im.abs());
    }
    Ok(Local {
        m,
        shear: wave.shear_number,
        radius,
        kappa,
    })
}

/// Per-station record of a tracked propagation, mouth (index 0) to inlet:
/// the scaled impedance at every station and each step's inverse flow map,
/// `u_j = Y_j^-1 u_{j+1}`.
struct Stations {
    zeta: Vec<Vec<C64>>,
    flow_back: Vec<Vec<C64>>,
}

/// Inlet impedance (scaled), diagnostics, and optionally the flow transfer
/// `u_in = G u_out` with the station record (only for loads with a finite
/// plane-mode impedance).
struct Propagated {
    zeta: Vec<C64>,
    zref: f64,
    min_shear: f64,
    throat: f64,
    mouth_ka: f64,
    flow: Option<Vec<C64>>,
    stations: Option<Stations>,
}

#[allow(clippy::too_many_arguments)]
fn propagate<'a>(
    profile: &RadiusProfile,
    state: &GasState,
    omega: f64,
    loss: LossModel,
    load: impl Into<ProfileLoad<'a>>,
    n: usize,
    steps: RiccatiSteps,
    track_flow: bool,
) -> Result<Propagated, DuctError> {
    let load = load.into();
    steps.validate()?;
    if !(omega.is_finite() && omega > 0.0) {
        return Err(bad("omega must be finite and positive"));
    }
    if !(state.density > 0.0 && state.sound_speed > 0.0) {
        return Err(bad(
            "modal profile gas must have positive density and sound speed",
        ));
    }
    let basis = Basis::new(n);
    let w = 2 * n;
    let r_in = profile.inlet_radius();
    let zref = state.density * state.sound_speed / (core::f64::consts::PI * r_in * r_in);
    if !(zref.is_finite() && zref > 0.0) {
        return Err(bad("modal profile reference impedance is not finite"));
    }
    let length = profile.length();
    let r_mouth = profile.radius(length)?;
    // Columns of [X; Y] span admissible mouth states: P/Zref = X c, U = Y c.
    let mut x_blk = vec![C64::ZERO; n * n];
    let mut y_blk = vec![C64::ZERO; n * n];
    let mouth_ka = match load {
        ProfileLoad::Scalar(load) => {
            let mouth = modal_wave(state, r_mouth, omega, loss, &basis.gammas[1..])?;
            let (plane, ka) = match load {
                MouthLoad::Analytic(t) => termination_impedance(t, state, r_mouth, omega)?,
                MouthLoad::Tabulated(table) => (
                    Some(table.z_at(omega)?),
                    omega * r_mouth / state.sound_speed,
                ),
                MouthLoad::Matched => (Some(mouth.zc[0]), omega * r_mouth / state.sound_speed),
            };
            if track_flow && plane.is_none() {
                return Err(bad("flow tracking needs a finite plane-mode mouth impedance"));
            }
            for k in 0..n {
                let (xk, yk) = match (k, plane) {
                    // A rigid cap annihilates every modal normal velocity.
                    (_, None) => (C64::ONE, C64::ZERO),
                    // Zero pressure is likewise a whole-aperture condition.
                    (_, Some(z))
                        if k == 0
                            || matches!(load, MouthLoad::Analytic(Termination::IdealOpen)) =>
                    {
                        (z.scale(1.0 / zref), C64::ONE)
                    }
                    _ => (mouth.zc[k].scale(1.0 / zref), C64::ONE),
                };
                if !finite(xk) {
                    return Err(DuctError::Singular);
                }
                x_blk[k * n + k] = xk;
                y_blk[k * n + k] = yk;
            }
            ka
        }
        ProfileLoad::Flanged(exterior) => {
            crate::exterior::admit(state, exterior)?;
            let matrix = crate::flanged::flanged_modal_impedance(exterior, r_mouth, omega, n)?;
            for (x, z) in x_blk.iter_mut().zip(matrix) {
                *x = z.scale(1.0 / zref);
                if !finite(*x) {
                    return Err(DuctError::Singular);
                }
            }
            for k in 0..n {
                y_blk[k * n + k] = C64::ONE;
            }
            // No Zc_n is formed here: the finite aperture matrix and the
            // original generator stay regular at a bare modal cutoff.
            omega * r_mouth / exterior.sound_speed
        }
    };
    let mut flow = track_flow.then(|| {
        let mut g = vec![C64::ZERO; n * n];
        for k in 0..n {
            g[k * n + k] = C64::ONE;
        }
        g
    });
    let mut stations = track_flow.then(|| Stations {
        zeta: vec![x_blk.clone()],
        flow_back: Vec::with_capacity(steps.0),
    });
    let h = -length / steps.0 as f64;
    let root3 = det::sqrt(3.0);
    let (c1, c2) = (0.5 - root3 / 6.0, 0.5 + root3 / 6.0);
    let mut min_shear = f64::INFINITY;
    let mut throat = r_mouth.min(r_in);
    for j in 0..steps.0 {
        let x0 = length + j as f64 * h;
        let l1 = generator(profile, state, omega, loss, &basis, zref, x0 + c1 * h)?;
        let l2 = generator(profile, state, omega, loss, &basis, zref, x0 + c2 * h)?;
        if l1.kappa.max(l2.kappa) * h.abs() > MAX_EVANESCENT_STEP {
            return Err(bad(
                "evanescent decay per step exceeds MAX_EVANESCENT_STEP; refine Riccati steps",
            ));
        }
        min_shear = min_shear.min(l1.shear).min(l2.shear);
        throat = throat.min(l1.radius).min(l2.radius);
        let (m1, m2) = (l1.m, l2.m);
        // Omega = h/2 (M1 + M2) - sqrt(3)/12 h^2 [M1, M2].
        let p12 = mul(&m1, &m2, w);
        let p21 = mul(&m2, &m1, w);
        let (half, comm) = (0.5 * h, -root3 * h * h / 12.0);
        let omega_m: Vec<C64> = (0..w * w)
            .map(|idx| (m1[idx] + m2[idx]).scale(half) + (p12[idx] - p21[idx]).scale(comm))
            .collect();
        let phi = expm(&omega_m, w)?;
        let mut x_new = vec![C64::ZERO; n * n];
        let mut y_new = vec![C64::ZERO; n * n];
        for r in 0..n {
            for c in 0..n {
                let (mut xs, mut ys) = (C64::ZERO, C64::ZERO);
                for l in 0..n {
                    xs = xs
                        + phi[r * w + l] * x_blk[l * n + c]
                        + phi[r * w + n + l] * y_blk[l * n + c];
                    ys = ys
                        + phi[(n + r) * w + l] * x_blk[l * n + c]
                        + phi[(n + r) * w + n + l] * y_blk[l * n + c];
                }
                x_new[r * n + c] = xs;
                y_new[r * n + c] = ys;
            }
        }
        if let Some(g) = &mut flow {
            *g = mul(&y_new, g, n);
        }
        // Moebius normalization: zeta = X' Y'^-1, then [X; Y] = [zeta; I].
        let y_inv = mat_inv(&y_new, n)?;
        x_blk = mul(&x_new, &y_inv, n);
        if x_blk.iter().any(|&v| !finite(v)) {
            return Err(DuctError::Singular);
        }
        if let Some(record) = &mut stations {
            record.zeta.push(x_blk.clone());
            record.flow_back.push(y_inv);
        }
        y_blk.fill(C64::ZERO);
        for k in 0..n {
            y_blk[k * n + k] = C64::ONE;
        }
    }
    Ok(Propagated {
        zeta: x_blk,
        zref,
        min_shear,
        throat,
        mouth_ka,
        flow,
        stations,
    })
}

/// Continuous multimodal input impedance of a smooth hole-free profile.
///
/// `modes` and `steps` are separate refinable axes, both reported. The
/// response's `tail_bound_estimate` is `|Z00(N) - Z00(N/2)|` from a second
/// propagation at half the mode count (none for `N = 1`).
///
/// # Errors
/// [`DuctError`]: a mode count outside `1..=certified_mode_cap()`, steps
/// outside `1..=MAX_RICCATI_STEPS`, a step exponential that needs more
/// resolution, the scalar termination's or table's refusals, and singular
/// Moebius algebra.
pub fn mm_profile_impedance(
    profile: &RadiusProfile,
    state: &GasState,
    omega: f64,
    loss: LossModel,
    load: MouthLoad<'_>,
    modes: ModeCount,
    steps: RiccatiSteps,
) -> Result<ModalResponse, DuctError> {
    profile_impedance_with_load(profile, state, omega, loss, ProfileLoad::Scalar(load), modes, steps)
}

/// Continuous multimodal response with complete flanged-mouth radiation.
///
/// The exterior impedance matrix is rebuilt for each retained-mode rung,
/// including the half-order comparison used by `tail_bound_estimate`. Every
/// original signed cross-mode entry loads the continuous propagation. The bore
/// retains `state`; radiation uses `exterior`, with a common static pressure.
/// The two mode/step axes and the existing evanescent-step bound are unchanged.
///
/// # Errors
/// As [`mm_profile_impedance`] and
/// [`crate::flanged::flanged_modal_impedance`], or incompatible static pressure.
#[allow(clippy::too_many_arguments)]
pub fn mm_profile_impedance_flanged(
    profile: &RadiusProfile,
    state: &GasState,
    exterior: &GasState,
    omega: f64,
    loss: LossModel,
    modes: ModeCount,
    steps: RiccatiSteps,
) -> Result<ModalResponse, DuctError> {
    profile_impedance_with_load(profile, state, omega, loss, ProfileLoad::Flanged(exterior), modes, steps)
}

fn profile_impedance_with_load(
    profile: &RadiusProfile,
    state: &GasState,
    omega: f64,
    loss: LossModel,
    load: ProfileLoad<'_>,
    modes: ModeCount,
    steps: RiccatiSteps,
) -> Result<ModalResponse, DuctError> {
    modes.validate()?;
    let n = modes.0;
    let full = propagate(profile, state, omega, loss, load, n, steps, false)?;
    let z: Vec<C64> = full.zeta.iter().map(|v| v.scale(full.zref)).collect();
    let tail_bound_estimate = if n >= 2 {
        let half = propagate(profile, state, omega, loss, load, n / 2, steps, false)?;
        Some((z[0] - half.zeta[0].scale(half.zref)).abs())
    } else {
        None
    };
    let cutoffs = Basis::new(n).gammas[1..]
        .iter()
        .map(|g| state.sound_speed * g / (core::f64::consts::TAU * full.throat))
        .collect();
    Ok(ModalResponse {
        omega,
        n_modes: n,
        plane_impedance: z[0],
        impedance_matrix: z,
        min_shear_number: full.min_shear,
        mouth_ka: full.mouth_ka,
        cutoffs_at_throat_hz: cutoffs,
        steps_used: steps.0,
        tail_bound_estimate,
    })
}

/// Modal pressure and flow at every step station of a continuous
/// propagation driven by a prescribed inlet modal flow.
#[derive(Debug, Clone)]
pub struct ModalField {
    /// Retained mode count, including the plane mode.
    pub n_modes: usize,
    /// Station distances from the inlet [m], ascending: the step nodes.
    pub x: Vec<f64>,
    /// Modal pressures `P_n` [Pa], one row of `n_modes` per station.
    pub pressure: Vec<Vec<C64>>,
    /// Modal flows `U_n` [m^3/s], one row of `n_modes` per station.
    pub flow: Vec<Vec<C64>>,
    /// Pressure on the axis `r = 0` [Pa] per station,
    /// `sum_n P_n psi_n(0)` with `psi_n(0) = 1/|J0(gamma_n)|`.
    pub axis_pressure: Vec<C64>,
}

fn mat_vec(a: &[C64], v: &[C64], n: usize) -> Vec<C64> {
    (0..n)
        .map(|i| (0..n).fold(C64::ZERO, |acc, l| acc + a[i * n + l] * v[l]))
        .collect()
}

/// Continuous multimodal field of a smooth hole-free profile driven at the
/// inlet by the modal flow `inlet_flow` (`U_n` [m^3/s]). A plane piston of
/// volume flow `q` is `[q, 0, ..]`: every higher mode integrates to zero over
/// the section. Memory is `O(steps N^2)`.
///
/// # Errors
/// As [`mm_profile_impedance`], plus an inlet flow that is not one finite
/// value per retained mode, and a mouth load with no finite plane-mode
/// impedance (a closed end).
#[allow(clippy::too_many_arguments)]
pub fn mm_profile_field(
    profile: &RadiusProfile,
    state: &GasState,
    omega: f64,
    loss: LossModel,
    load: MouthLoad<'_>,
    modes: ModeCount,
    steps: RiccatiSteps,
    inlet_flow: &[C64],
) -> Result<ModalField, DuctError> {
    profile_field_with_load(profile, state, omega, loss, ProfileLoad::Scalar(load), modes, steps, inlet_flow)
}

/// Continuous modal pressure/flow field with complete flanged-mouth radiation.
///
/// `inlet_flow` and all returned fields retain the original signed modal basis.
/// At the mouth `P = Z_rad U`; radiated power is `Re(U^H P)/2`, and every
/// higher-mode flow remains available to physical observers.
///
/// # Errors
/// As [`mm_profile_field`] and [`mm_profile_impedance_flanged`].
#[allow(clippy::too_many_arguments)]
pub fn mm_profile_field_flanged(
    profile: &RadiusProfile,
    state: &GasState,
    exterior: &GasState,
    omega: f64,
    loss: LossModel,
    modes: ModeCount,
    steps: RiccatiSteps,
    inlet_flow: &[C64],
) -> Result<ModalField, DuctError> {
    profile_field_with_load(profile, state, omega, loss, ProfileLoad::Flanged(exterior), modes, steps, inlet_flow)
}

#[allow(clippy::too_many_arguments)]
fn profile_field_with_load(
    profile: &RadiusProfile,
    state: &GasState,
    omega: f64,
    loss: LossModel,
    load: ProfileLoad<'_>,
    modes: ModeCount,
    steps: RiccatiSteps,
    inlet_flow: &[C64],
) -> Result<ModalField, DuctError> {
    modes.validate()?;
    let n = modes.0;
    if inlet_flow.len() != n || !inlet_flow.iter().all(|&v| finite(v)) {
        return Err(bad(
            "inlet flow must hold one finite modal flow per retained mode",
        ));
    }
    let run = propagate(profile, state, omega, loss, load, n, steps, true)?;
    let Some(record) = run.stations else {
        return Err(bad(
            "flow tracking needs a finite plane-mode mouth impedance",
        ));
    };
    let axis_weight: Vec<f64> = Basis::new(n)
        .gammas
        .iter()
        .map(|&g| {
            if g == 0.0 {
                1.0
            } else {
                1.0 / bessel_j0_certified(g).abs()
            }
        })
        .collect();
    let s = record.flow_back.len();
    let (length, h) = (profile.length(), -profile.length() / s as f64);
    let mut field = ModalField {
        n_modes: n,
        x: Vec::with_capacity(s + 1),
        pressure: Vec::with_capacity(s + 1),
        flow: Vec::with_capacity(s + 1),
        axis_pressure: Vec::with_capacity(s + 1),
    };
    // Station j counts steps from the mouth; walk back from the inlet (j = s).
    let mut u = inlet_flow.to_vec();
    for j in (0..=s).rev() {
        if j < s {
            u = mat_vec(&record.flow_back[j], &u, n);
        }
        let p: Vec<C64> = mat_vec(&record.zeta[j], &u, n)
            .iter()
            .map(|v| v.scale(run.zref))
            .collect();
        let axis = p
            .iter()
            .zip(&axis_weight)
            .fold(C64::ZERO, |acc, (pn, w)| acc + pn.scale(*w));
        if !(finite(axis) && u.iter().all(|&v| finite(v))) {
            return Err(DuctError::Singular);
        }
        field
            .x
            .push(if j == s { 0.0 } else { length + j as f64 * h });
        field.pressure.push(p);
        field.flow.push(u.clone());
        field.axis_pressure.push(axis);
    }
    Ok(field)
}

#[cfg(test)]
mod tests;
