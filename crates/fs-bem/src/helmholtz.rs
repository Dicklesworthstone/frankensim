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
//! the HBIE row written below as `N p - D' q - q/2 = 0`, the production
//! coupling is `alpha = +i/k`. Any `alpha` with a nonzero imaginary part
//! restores unique solvability; the sign sets the conditioning and how the
//! hypersingular discretization error enters. It is pinned by measurement,
//! not by a convention label: on the pulsating sphere with the exact static
//! rows (`N_0[1] = 0`, `(D_0 + I/2)1 = 0`), `+i/k` converges at the plain
//! CBIE rate away from interior resonances and is 3-5x more accurate than
//! `-i/k` for `ka >= 2`, while both rescue the resonances CBIE loses.
//!
//! FORMULATION — exterior Green representation with outward normals,
//! collocation at panel centroids, piecewise-constant elements:
//! `p(x) = INT [ dG/dn_y p - G q ] dS` for x in the exterior. Surface
//! limits give
//! CBIE: `D p - S q - p/2 = 0`, and its `n_x` derivative
//! HBIE: `N p - D' q - q/2 = 0` (hypersingular finite part), combined as
//! `(D - I/2 + alpha N) p = (S + alpha D' + alpha/2 I) q`.
//!
//! DISCRETIZATION HONESTY: Plain CBIE surfaces retaining oriented triangles
//! integrate static single/double layers analytically for every triangle, then
//! use tensor Gauss integration for their nonsingular Helmholtz remainders:
//! 8x8 near a panel and 4x4 for well-separated interactions. This resolves
//! opposite faces of thin shells without sampling a narrow static kernel spike,
//! preserves the closed-surface static double-layer identity, and retains the
//! far-panel spatial moments needed by weakly radiating bending motions.
//! Burton-Miller and legacy centroid-only surfaces retain equivalent-disc self
//! terms and point-panel off-diagonal approximation. Triangle CBIE self terms use
//! a centroid fan plus a Duffy transform, which cancels the weak `1/r`
//! singularity and the regularized hypersingular `1/r` remainder. The legacy
//! equivalent-disc finite parts use
//! `a = sqrt(A/pi)`:
//! `S_ii = (e^{ika} - 1) / (2ik)`, `D_ii = D'_ii = 0` (flat panel), and
//! the hypersingular self entry split as the integrable difference
//! `(N_k - N_0)_ii = (ik - (e^{ika} - 1)/a)/2` plus the static part
//! recovered from the exact closed-surface identity `N_0[1] = 0` as
//! minus the discrete off-diagonal static row sum (same quadrature on
//! both sides, so the point-panel error cancels). Burton–Miller also closes
//! `(D_0+I/2)1=0` with its static double-layer self coefficient and uses the
//! same diagonal for `D'_0`, its area-weighted adjoint. Without this paired
//! consistency, the `1/k` coupling amplifies a spurious static net flux.
//! These paths are pinned by analytic/independent quadrature tests and the
//! whole pipeline is arbitrated by the pulsating-sphere impedance oracle.
//!
//! Determinism: assembly and solves are sequential with fixed traversal
//! order; complex exponentials route through `fs_math::det` per the
//! workspace libm doctrine. Repeat solves are bitwise identical.
//!
//! Deferred with recorded triggers (see CONTRACT): wideband directional
//! FMM acceleration (dense complex LU is the v1 path; the work cap
//! refuses above [`MAX_DENSE_PANELS`]) and a Galerkin/higher-order boundary
//! basis beyond the constant-panel triangle integration. Spherical-harmonic
//! directivity tables ([`directivity_sh_table`]), radiation efficiency
//! ([`radiation_efficiency`]), and per-solve condition/work-cap
//! diagnostics (on [`RadiationSolution`]) ship in this module; the
//! condition number is an honest probe-based LOWER BOUND, not a
//! rigorous estimator (the fictitious-frequency contrast remains the
//! physics-level diagnostic).

use fs_la::eigen_complex::{LuComplex, lu_complex};
use fs_material::gas::{GasSpec, GasState};
use fs_math::c64::C64;
use fs_math::det;

use crate::panel3d::SpherePanels;

#[path = "helmholtz_laplace.rs"]
pub mod laplace;

/// One source-owned molecular gas for radiation reaction and finite receivers.
#[path = "helmholtz_molecular.rs"]
pub mod molecular;

#[path = "cbie_geometry.rs"]
mod cbie_geometry;
pub use cbie_geometry::PreparedCbieGeometry;

/// FMM Burton-Miller operator and restart-resumable, right-preconditioned GMRES.
#[path = "helmholtz_iterative.rs"]
pub mod iterative;
pub use iterative::{
    IterativeError, IterativeOptions, IterativeRadiation, IterativeState, Preconditioner,
    ResourceReceipt, solve_radiation_iterative,
};

/// Locally/nonlocally reacting boundaries on the same Burton-Miller owner.
#[path = "helmholtz_bc.rs"]
pub mod bc;

/// Declared triangle near rule for the closed Burton-Miller rows.
#[path = "helmholtz_triangle_rule.rs"]
mod triangle_rule;
pub use triangle_rule::{PanelRule, solve_radiation_batch_with_rule, solve_radiation_with_rule};

/// Reciprocal Rayleigh loading of complete planar modal velocity fields.
#[path = "rayleigh_modal.rs"]
pub mod rayleigh_modal;

/// Self and mutual radiation impedance of declared apertures on one body.
#[path = "aperture_impedance.rs"]
pub mod aperture;

/// Dense-LU work cap: a complex dense system above this panel count is
/// refused rather than silently thrashing (n^2 * 16 bytes; 8192 panels
/// is already a 1 GiB matrix). The FMM path is the recorded follow-up.
pub const MAX_DENSE_PANELS: usize = 8192;

std::thread_local! {
    static DENSE_WORKERS: core::cell::Cell<usize> = const { core::cell::Cell::new(1) };
}

/// Run `f` with `workers` (at least one) threads declared for this thread's
/// dense operator assembly and its LU factorization
/// ([`fs_la::eigen_complex::with_lu_workers`]), and for the frequencies of a
/// reciprocal modal-loading sweep
/// ([`crate::modal_loading::modal_radiation_loading_reciprocal`]), then restore
/// the previous declaration, also on unwind. This is a scheduling budget only:
/// every matrix entry and factor is bit-identical for every worker count.
pub fn with_dense_workers<R>(workers: usize, f: impl FnOnce() -> R) -> R {
    struct Restore(usize);
    impl Drop for Restore {
        fn drop(&mut self) {
            DENSE_WORKERS.with(|w| w.set(self.0));
        }
    }
    let _restore = Restore(DENSE_WORKERS.with(|w| w.replace(workers.max(1))));
    fs_la::eigen_complex::with_lu_workers(workers, f)
}

/// This thread's declared [`with_dense_workers`] count; 1 when undeclared.
pub(crate) fn dense_workers() -> usize {
    DENSE_WORKERS.with(core::cell::Cell::get)
}

/// Maximum number of boundary-velocity fields solved against one shared
/// Helmholtz factorization. This bounds the aggregate output and triangular
/// solve work while still covering modal/directivity batches without forcing
/// callers back to one factorization per field.
pub const MAX_RADIATION_FIELDS_PER_BATCH: usize = 256;

/// Minimum panels per wavelength before the solver refuses: below ~6 the
/// centroid-collocation solution degrades silently, so the boundary is a
/// named refusal instead.
pub const MIN_PANELS_PER_WAVELENGTH: f64 = 6.0;

const ROOM_TEMPERATURE_K: f64 = 293.15;
const STANDARD_ATMOSPHERE_PA: f64 = 101_325.0;

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
    /// The optional reusable static-geometry table exceeds its byte budget.
    CacheBudget {
        /// Logical bytes needed for the complete table.
        required_bytes: usize,
        /// Caller-admitted maximum logical bytes.
        maximum_bytes: usize,
    },
    /// Bounded geometry-cache storage could not be reserved.
    AllocationFailed {
        /// Which allocation refused.
        what: &'static str,
    },
    /// Geometry preparation was cancelled between target rows.
    Cancelled,
    /// The dense complex LU reported singularity.
    Singular,
    /// A declared [`PanelRule::TriangleNear`] admits fewer corrected panel
    /// pairs than the surface's near zone holds.
    NearRuleBudget {
        /// Near (target, source) pairs the surface requires.
        required_pairs: usize,
        /// Declared maximum.
        maximum_pairs: usize,
    },
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
            HelmholtzError::CacheBudget {
                required_bytes,
                maximum_bytes,
            } => write!(
                f,
                "FS-BEM-HELM-CACHE-BUDGET: {required_bytes} static-table bytes exceeds the admitted {maximum_bytes} bytes"
            ),
            HelmholtzError::AllocationFailed { what } => {
                write!(f, "FS-BEM-HELM-ALLOCATION: could not reserve {what}")
            }
            HelmholtzError::Cancelled => {
                write!(f, "FS-BEM-HELM-CANCELLED: geometry preparation cancelled")
            }
            HelmholtzError::Singular => write!(f, "FS-BEM-HELM-SINGULAR: dense LU refused"),
            HelmholtzError::NearRuleBudget {
                required_pairs,
                maximum_pairs,
            } => write!(
                f,
                "FS-BEM-HELM-NEAR-RULE-BUDGET: {required_pairs} near panel pairs exceeds the declared {maximum_pairs}"
            ),
        }
    }
}

impl std::error::Error for HelmholtzError {}

/// Acoustic medium (SI): density and sound speed.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Medium {
    /// Ambient density rho [kg/m^3].
    pub density: f64,
    /// Sound speed c [m/s].
    pub sound_speed: f64,
}

impl From<&GasState> for Medium {
    /// Copy the resolved medium's SI values exactly. Source authority remains
    /// with its producer; this view evaluates no gas law or default atmosphere.
    fn from(gas: &GasState) -> Self {
        Self {
            density: gas.density,
            sound_speed: gas.sound_speed,
        }
    }
}

impl Medium {
    /// Copy a checked gas or liquid acoustic state into the radiation medium.
    /// The source retains its original EOS, complete query and usage receipts;
    /// this numerical view introduces no gas law or atmosphere default.
    ///
    /// # Errors
    /// Refuses invalid or incoherent material acoustic states.
    pub fn from_material<M: fs_material::fluid_acoustics::FluidAcoustics + ?Sized>(
        material: &M,
    ) -> Result<Self, fs_material::MaterialError> {
        let state = material.acoustic_state()?;
        Ok(Self {
            density: state.density_kg_m3(),
            sound_speed: state.sound_speed_m_s(),
        })
    }

    /// Dry USSA-1976 air at the declared room state of 293.15 K and 101325 Pa.
    /// Parameterized studies should construct their own [`GasState`] for their
    /// actual temperature, pressure, and gas identity instead of using this
    /// convenience.
    #[must_use]
    pub fn air() -> Medium {
        let gas = GasState::try_new(
            &GasSpec::dry_air_ussa1976(),
            ROOM_TEMPERATURE_K,
            STANDARD_ATMOSPHERE_PA,
        )
        .expect("declared dry-air room state is within GasState's validity window");
        Medium {
            density: gas.density,
            sound_speed: gas.sound_speed,
        }
    }
}

/// Boundary-integral formulation choice.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Formulation {
    /// Plain CBIE — breaks down at interior resonances (kept as the
    /// documented near-miss arm of the fictitious-frequency contrast).
    PlainCbie,
    /// Burton–Miller CBIE + alpha HBIE with `alpha = +i/k` for this
    /// module's row-sign convention — the production arm.
    BurtonMiller,
    /// The alpha-sign MUTATION arm (`alpha = -i/k` on this module's
    /// row convention): kept public and documented so the conformance
    /// battery can prove the wrong sign measurably degrades accuracy and
    /// convergence instead of passing silently. Never use for production
    /// solves.
    BurtonMillerWrongAlphaSign,
}

/// Scale diagnostics for the discrete surface-power quadrature.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RadiationPowerDiagnostics {
    /// Raw complex surface power `1/2 SUM p conj(v) A` [W].
    pub surface_power: C64,
    /// Sum of non-negative real panel-power terms [W].
    pub positive_real_power: f64,
    /// Sum of negative real panel-power terms [W].
    pub negative_real_power: f64,
    /// Sum of panel-power magnitudes `1/2 SUM |p| |v| A` [W].
    pub apparent_power: f64,
    /// Plane-wave reference `1/2 rho c SUM |v|^2 A` [W].
    pub plane_wave_reference_power: f64,
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
    /// Acoustic medium used to turn velocity into Neumann data.
    pub medium: Medium,
    /// Source-panel fingerprint: a cross-wiring guard, not a certificate.
    pub surface_fingerprint: u64,
    /// Measured panels per wavelength (the refusal diagnostic).
    pub panels_per_wavelength: f64,
    /// Radiated power W = 1/2 Re SUM p conj(v) A [W].
    pub radiated_power: f64,
    /// Outward-rounded enclosure of the exact panel-power dot product [W].
    ///
    /// This covers only floating-point evaluation of already-computed
    /// pressure, velocity, and area values. It does not bound BEM solve,
    /// quadrature, discretization, or model error. If this interval straddles
    /// zero, a negative rounded sum is published as neutral zero rather than
    /// as evidence of non-passivity.
    pub radiated_power_roundoff_interval: (f64, f64),
    /// Signed, reactive, cancellation, and plane-wave reference scales.
    pub power_diagnostics: RadiationPowerDiagnostics,
    /// LOWER BOUND on the 1-norm condition number of the assembled
    /// system, from `||A||_1` times the largest `||A^{-1} b||_1` over a
    /// fixed set of deterministic unit-1-norm probe solves. A large
    /// value is a reliable warning (plain CBIE near a fictitious
    /// frequency measurably inflates it); a small value is NOT a
    /// conditioning certificate — the probes can miss the worst
    /// direction. The iterative owner reports the universal lower bound 1
    /// without dense probes. Rigorous estimation (Hager/Higham with adjoint
    /// solves) is the recorded follow-up.
    pub condition_lower_bound: f64,
    /// Work-cap utilization: panel count over [`MAX_DENSE_PANELS`], in
    /// (0, 1] for dense solves, potentially >1 for iterative solves.
    /// The headroom diagnostic the pilot intel asked to surface
    /// instead of refusing opaquely at the cap.
    pub dense_cap_utilization: f64,
}

fn expik(kr: f64) -> C64 {
    C64::new(det::cos(kr), det::sin(kr))
}

/// Free-space kernel `G = e^{ikr}/(4 pi r)`.
pub(crate) fn green<K: WaveNumber>(k: K, r: f64) -> C64 {
    k.outgoing(r).scale(1.0 / (4.0 * core::f64::consts::PI * r))
}

/// `G'(r) = (ik - 1/r) G`.
fn green_dr<K: WaveNumber>(k: K, r: f64) -> C64 {
    green(k, r) * k.radial_factor(r)
}

// Real-frequency solves retain their original operation graph. The damped
// extension supplies only the analytic scalar kernels; matrix geometry, static
// closures, normal derivatives and LU assembly remain this owner's code.
pub(crate) trait WaveNumber: Copy + Send + Sync {
    fn value(self) -> C64;
    fn outgoing(self, r: f64) -> C64;
    fn radial_factor(self, r: f64) -> C64;
    fn disc_self(self, area: f64) -> (C64, C64);
    fn cbie_real(self) -> f64;
}
impl WaveNumber for f64 {
    fn value(self) -> C64 {
        C64::from_re(self)
    }
    fn outgoing(self, r: f64) -> C64 {
        expik(self * r)
    }
    fn radial_factor(self, r: f64) -> C64 {
        C64::new(-1.0 / r, self)
    }
    fn disc_self(self, area: f64) -> (C64, C64) {
        real_self_terms(self, area)
    }
    fn cbie_real(self) -> f64 {
        self
    }
}

/// Analytic Laplace-subtracted single layer on the SAME WaveNumber seam.
/// The finite coincident limit is ik/(4pi); the series avoids subtracting
/// nearly equal real values on graded edge panels or damped contours.
pub(crate) fn green_remainder<K: WaveNumber>(wave: K, r: f64) -> C64 {
    let ik = C64::new(0.0, 1.0) * wave.value();
    if r == 0.0 {
        return ik.scale(1.0 / (4.0 * core::f64::consts::PI));
    }
    let z = ik.scale(r);
    if z.abs() < 0.25 {
        let mut term = ik;
        let mut sum = term;
        for n in 2..=24 {
            term = (term * z).scale(1.0 / f64::from(n));
            sum = sum + term;
        }
        sum.scale(1.0 / (4.0 * core::f64::consts::PI))
    } else {
        green(wave, r) - C64::from_re(1.0 / (4.0 * core::f64::consts::PI * r))
    }
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

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn affine3(origin: [f64; 3], a: [f64; 3], wa: f64, b: [f64; 3], wb: f64) -> [f64; 3] {
    [
        origin[0] + wa * a[0] + wb * b[0],
        origin[1] + wa * a[1] + wb * b[1],
        origin[2] + wa * a[2] + wb * b[2],
    ]
}

/// The four point kernels for x != y: (G, dG/dn_y, dG/dn_x,
/// d2G/dn_x dn_y). With `d = x - y`, `r = |d|`:
/// `dG/dn_y = -G' (n_y . d)/r`, `dG/dn_x = G' (n_x . d)/r`, and
/// `d2G/dn_x dn_y = -(G''/r^2 - G'/r^3)(n_x . d)(n_y . d) - (G'/r)(n_x . n_y)`
/// where `G'' = (ik - 1/r) G' + G/r^2`. All four are pinned against
/// central finite differences of G in the unit tests.
fn kernels<K: WaveNumber>(
    k: K,
    x: [f64; 3],
    nx: [f64; 3],
    y: [f64; 3],
    ny: [f64; 3],
) -> (C64, C64, C64, C64) {
    let d = sub(x, y);
    let r = norm(d);
    let g = green(k, r);
    let gp = green_dr(k, r);
    let gpp = gp * k.radial_factor(r) + g.scale(1.0 / (r * r));
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
fn self_terms<K: WaveNumber>(k: K, area: f64) -> (C64, C64) {
    k.disc_self(area)
}

#[cfg(test)]
mod real_kernel_bits {
    use super::*;
    #[test]
    fn g0_real_k_point_kernel_bits_match_preextension_arithmetic() {
        // Verbatim pre-extension scalar expressions. The generic f64 path
        // must retain these bits: dense and iterative assembly share them.
        fn original(k: f64, x: [f64; 3], nx: [f64; 3], y: [f64; 3], ny: [f64; 3]) -> [C64; 4] {
            let d = sub(x, y);
            let r = norm(d);
            let g = expik(k * r).scale(1.0 / (4.0 * core::f64::consts::PI * r));
            let gp = g * C64::new(-1.0 / r, k);
            let gpp = gp * C64::new(-1.0 / r, k) + g.scale(1.0 / (r * r));
            let nyd = dot(ny, d);
            let nxd = dot(nx, d);
            let radial = gpp.scale(1.0 / (r * r)) - gp.scale(1.0 / (r * r * r));
            [
                g,
                gp.scale(-nyd / r),
                gp.scale(nxd / r),
                radial.scale(-nxd * nyd) - gp.scale(dot(nx, ny) / r),
            ]
        }
        for k in [0.0, 0.0001, 0.5, core::f64::consts::PI, 137.0] {
            for (x, y) in [
                ([0.2, 0.7, -0.4], [0.0, 0.1, 0.3]),
                ([2.0, -0.3, 0.7], [-0.2, 0.0, 0.1]),
            ] {
                let nx = [0.0, 0.6, 0.8];
                let ny = [1.0, 0.0, 0.0];
                let old = original(k, x, nx, y, ny);
                let new = kernels(k, x, nx, y, ny);
                for (old, new) in old.into_iter().zip([new.0, new.1, new.2, new.3]) {
                    assert_eq!(old.re.to_bits(), new.re.to_bits());
                    assert_eq!(old.im.to_bits(), new.im.to_bits());
                }
            }
            for area in [1e-8, 0.007, 1.0] {
                let a = self_terms(k, area);
                let b = real_self_terms(k, area);
                assert_eq!(
                    (
                        a.0.re.to_bits(),
                        a.0.im.to_bits(),
                        a.1.re.to_bits(),
                        a.1.im.to_bits()
                    ),
                    (
                        b.0.re.to_bits(),
                        b.0.im.to_bits(),
                        b.1.re.to_bits(),
                        b.1.im.to_bits()
                    )
                );
            }
        }
    }
}
fn real_self_terms(k: f64, area: f64) -> (C64, C64) {
    let a = (area / core::f64::consts::PI).sqrt();
    if k == 0.0 {
        return (C64::new(0.5 * a, 0.0), C64::ZERO);
    }
    let ik = C64::new(0.0, k);
    let s = (expik(k * a) - C64::ONE) * (ik + ik).recip();
    let n_reg = (ik - (expik(k * a) - C64::ONE).scale(1.0 / a)).scale(0.5);
    (s, n_reg)
}

// Eight-point Gauss-Legendre rule mapped to [0, 1]. The fixed table keeps the
// integration order and floating-point operation graph deterministic.
const GL8_X: [f64; 8] = [
    0.019_855_071_751_231_884,
    0.101_666_761_293_186_64,
    0.237_233_795_041_835_5,
    0.408_282_678_752_175_1,
    0.591_717_321_247_824_9,
    0.762_766_204_958_164_5,
    0.898_333_238_706_813_4,
    0.980_144_928_248_768_1,
];
const GL8_W: [f64; 8] = [
    0.050_614_268_145_188_13,
    0.111_190_517_226_687_24,
    0.156_853_322_938_943_66,
    0.181_341_891_689_181,
    0.181_341_891_689_181,
    0.156_853_322_938_943_66,
    0.111_190_517_226_687_24,
    0.050_614_268_145_188_13,
];

// Four-point Gauss-Legendre rule on [0, 1]. Well-separated regular remainders
// do not need the near-panel 8x8 rule, but a centroid loses their spatial
// moments. Those errors can overwhelm a bending mode's small radiation power.
const GL4_X: [f64; 4] = [
    0.069_431_844_202_973_71,
    0.330_009_478_207_571_87,
    0.669_990_521_792_428_1,
    0.930_568_155_797_026_2,
];
const GL4_W: [f64; 4] = [
    0.173_927_422_568_726_92,
    0.326_072_577_431_273_07,
    0.326_072_577_431_273_07,
    0.173_927_422_568_726_92,
];

fn regularized_coplanar_hypersingular(k: f64, r: f64) -> C64 {
    if k == 0.0 {
        return C64::ZERO;
    }
    let z = k * r;
    let numerator = if z.abs() < 1.0e-3 {
        let z2 = z * z;
        let z3 = z2 * z;
        let z4 = z2 * z2;
        let z5 = z4 * z;
        let z6 = z3 * z3;
        C64::new(z2 / 2.0 - z4 / 8.0 + z6 / 144.0, z3 / 3.0 - z5 / 30.0)
    } else {
        expik(z) * C64::new(1.0, -z) - C64::ONE
    };
    numerator.scale(1.0 / (4.0 * core::f64::consts::PI * r * r * r))
}

fn regularized_single_layer(k: f64, r: f64) -> C64 {
    let z = k * r;
    // cos(z)-1 = -2*sin(z/2)^2 remains accurate for nearly touching skins.
    let half_sine = det::sin(0.5 * z);
    C64::new(-2.0 * half_sine * half_sine, det::sin(z))
        .scale(1.0 / (4.0 * core::f64::consts::PI * r))
}

/// Static single and double layers, with the same outward source normal as the
/// Helmholtz kernels. The solid angle is the van Oosterom-Strackee expression
/// (doi:10.1109/TBME.1983.325207); the single layer follows the polygon edge
/// integral of Wilton et al. (doi:10.1109/TAP.1984.1143304).
///
/// In the triangle plane, div((y-x_parallel)/r) = 1/r + h^2/r^3. Its boundary
/// integral therefore gives S = (sum edge_distance * integral_edge(1/r)
/// - h * solid_angle)/(4*pi), while D = solid_angle/(4*pi). This separation
/// retains the double-layer jump as h tends to zero instead of asking a fixed
/// area quadrature to find an increasingly narrow h/r^3 peak.
fn triangle_static_influence(x: [f64; 3], triangle: [[f64; 3]; 3], ny: [f64; 3]) -> (f64, f64) {
    let relative = triangle.map(|p| sub(p, x));
    let radii = relative.map(norm);
    let height = -dot(relative[0], ny);
    let double_area = norm(cross(
        sub(triangle[1], triangle[0]),
        sub(triangle[2], triangle[0]),
    ));
    let denominator = radii[0] * radii[1] * radii[2]
        + dot(relative[0], relative[1]) * radii[2]
        + dot(relative[1], relative[2]) * radii[0]
        + dot(relative[2], relative[0]) * radii[1];
    // Coplanar off-panel targets have zero ordinary double layer. The self
    // panel has its separately specified trace/jump in assemble_dense.
    let solid_angle = if height == 0.0 {
        0.0
    } else {
        2.0 * det::atan2(height * double_area, denominator)
    };
    let mut edge_sum = 0.0;
    for edge in 0..3 {
        let next = (edge + 1) % 3;
        let direction = sub(triangle[next], triangle[edge]);
        let length = norm(direction);
        let tangent = direction.map(|v| v / length);
        let outward = cross(tangent, ny);
        let distance = dot(relative[edge], outward);
        // The product tends to zero even when the projection lies on this
        // edge and its logarithmic line integral diverges.
        if distance == 0.0 {
            continue;
        }
        let start = dot(relative[edge], tangent);
        let end = start + length;
        // Avoid subtracting almost equal endpoint radii for a target close to
        // an edge. These are algebraic forms of asinh(end/d)-asinh(start/d).
        let integral = if start >= 0.0 {
            det::ln((radii[next] + end) / (radii[edge] + start))
        } else if end <= 0.0 {
            det::ln((radii[edge] - start) / (radii[next] - end))
        } else {
            let perpendicular = norm(cross(tangent, relative[edge]));
            det::ln((radii[edge] - start) / perpendicular)
                + det::ln((radii[next] + end) / perpendicular)
        };
        edge_sum += distance * integral;
    }
    let four_pi = 4.0 * core::f64::consts::PI;
    (
        (edge_sum - height * solid_angle) / four_pi,
        solid_angle / four_pi,
    )
}

/// Integrate the weakly singular single-layer and ordinary double-layer
/// kernels over a flat source triangle after subtracting their exact static
/// singularities. The target must not lie on the source triangle; the self
/// panel uses [`triangle_self_terms`]. No geometry refinement is implied.
fn triangle_weak_influence(
    k: f64,
    x: [f64; 3],
    triangle: [[f64; 3]; 3],
    ny: [f64; 3],
) -> (C64, C64) {
    triangle_weak_influence_with_rule(k, x, triangle, ny, &GL8_X, &GL8_W)
}

fn triangle_weak_influence_with_rule<const ORDER: usize>(
    k: f64,
    x: [f64; 3],
    triangle: [[f64; 3]; 3],
    ny: [f64; 3],
    nodes: &[f64; ORDER],
    weights: &[f64; ORDER],
) -> (C64, C64) {
    triangle_weak_influence_with_static(
        k,
        x,
        triangle,
        ny,
        (nodes, weights),
        triangle_static_influence(x, triangle, ny),
    )
}

fn regularized_double_layer_remainder(k: f64, r: f64) -> C64 {
    if k == 0.0 {
        return C64::ZERO;
    }
    let z = k * r;
    if z.abs() < 1.0e-2 {
        let z2 = z * z;
        let z3 = z2 * z;
        return C64::new(
            -z / 8.0 + z3 / 144.0,
            1.0 / 3.0 - z2 / 30.0 + z2 * z2 / 840.0,
        )
        .scale(k * k * k / (4.0 * core::f64::consts::PI));
    }
    let numerator = expik(z) * C64::new(1.0, -z) - C64::new(1.0 + 0.5 * z * z, 0.0);
    numerator.scale(1.0 / (4.0 * core::f64::consts::PI * r * r * r))
}

fn triangle_weak_influence_with_static<const ORDER: usize>(
    k: f64,
    x: [f64; 3],
    triangle: [[f64; 3]; 3],
    ny: [f64; 3],
    rule: (&[f64; ORDER], &[f64; ORDER]),
    static_layers: (f64, f64),
) -> (C64, C64) {
    let (nodes, weights) = rule;
    let [a, b, c] = triangle;
    let ab = sub(b, a);
    let ac = sub(c, a);
    let double_area = norm(cross(ab, ac));
    let (static_single, static_double) = static_layers;
    let height = dot(sub(a, x), ny);
    let normal_offset = -height;
    let mut single = C64::from_re(static_single);
    let mut double = C64::from_re(static_double + 0.5 * k * k * normal_offset * static_single);
    let x_proj = [
        x[0] + height * ny[0],
        x[1] + height * ny[1],
        x[2] + height * ny[2],
    ];
    let sub_areas: [f64; 3] = core::array::from_fn(|edge| {
        let e0 = sub(triangle[edge], x_proj);
        let e1 = sub(triangle[(edge + 1) % 3], x_proj);
        dot(cross(e0, e1), ny)
    });
    if sub_areas.iter().all(|&sa| sa > 0.0) {
        for (edge, &sub_double_area) in sub_areas.iter().enumerate() {
            let e0 = sub(triangle[edge], x_proj);
            let e1 = sub(triangle[(edge + 1) % 3], x_proj);
            for (&radial, &wr) in nodes.iter().zip(weights) {
                for (&tangent, &wt) in nodes.iter().zip(weights) {
                    let ray = [
                        (1.0 - tangent).mul_add(e0[0], tangent * e1[0]),
                        (1.0 - tangent).mul_add(e0[1], tangent * e1[1]),
                        (1.0 - tangent).mul_add(e0[2], tangent * e1[2]),
                    ];
                    let y = [
                        radial.mul_add(ray[0], x_proj[0]),
                        radial.mul_add(ray[1], x_proj[1]),
                        radial.mul_add(ray[2], x_proj[2]),
                    ];
                    let d = sub(x, y);
                    let r = norm(d);
                    let weight = wr * wt * radial * sub_double_area;
                    let g = regularized_single_layer(k, r);
                    let dgdny = regularized_double_layer_remainder(k, r).scale(normal_offset);
                    single = single + g.scale(weight);
                    double = double + dgdny.scale(weight);
                }
            }
        }
        return (single, double);
    }
    for (&u, &wu) in nodes.iter().zip(weights) {
        for (&v, &wv) in nodes.iter().zip(weights) {
            let y = affine3(a, ab, u, ac, (1.0 - u) * v);
            let weight = wu * wv * (1.0 - u) * double_area;
            let d = sub(x, y);
            let r = norm(d);
            let g = regularized_single_layer(k, r);
            let dgdny = regularized_double_layer_remainder(k, r).scale(normal_offset);
            single = single + g.scale(weight);
            double = double + dgdny.scale(weight);
        }
    }
    (single, double)
}

/// `(S_ii, (N_k-N_0)_ii)` for a flat triangle collocated at its centroid. The
/// static S term uses the same exact integral as a nearby opposite face; a
/// quadrature mismatch between those terms would dominate their thin-gap limit.
/// The dynamic remainders use a centroid fan and Duffy transform, cancelling
/// the hypersingular remainder's integrable `1/r` behavior.
fn triangle_self_terms(k: f64, x: [f64; 3], triangle: [[f64; 3]; 3]) -> (C64, C64) {
    triangle_self_terms_with_static(k, x, triangle, triangle_self_static(x, triangle))
}

fn triangle_self_static(x: [f64; 3], triangle: [[f64; 3]; 3]) -> f64 {
    let area_normal = cross(sub(triangle[1], triangle[0]), sub(triangle[2], triangle[0]));
    let area_scale = norm(area_normal);
    let normal = area_normal.map(|v| v / area_scale);
    triangle_static_influence(x, triangle, normal).0
}

fn triangle_self_terms_with_static(
    k: f64,
    x: [f64; 3],
    triangle: [[f64; 3]; 3],
    static_single: f64,
) -> (C64, C64) {
    let mut single_layer = C64::from_re(static_single);
    let mut hypersingular_regular = C64::ZERO;
    for edge in 0..3 {
        let e0 = sub(triangle[edge], x);
        let e1 = sub(triangle[(edge + 1) % 3], x);
        let double_area = norm(cross(e0, e1));
        for (&radial, &wr) in GL8_X.iter().zip(&GL8_W) {
            for (&tangent, &wt) in GL8_X.iter().zip(&GL8_W) {
                let ray = [
                    (1.0 - tangent).mul_add(e0[0], tangent * e1[0]),
                    (1.0 - tangent).mul_add(e0[1], tangent * e1[1]),
                    (1.0 - tangent).mul_add(e0[2], tangent * e1[2]),
                ];
                let r = radial * norm(ray);
                let weight = wr * wt * radial * double_area;
                single_layer = single_layer + regularized_single_layer(k, r).scale(weight);
                hypersingular_regular =
                    hypersingular_regular + regularized_coplanar_hypersingular(k, r).scale(weight);
            }
        }
    }
    (single_layer, hypersingular_regular)
}

fn triangle_is_near(x: [f64; 3], centroid: [f64; 3], triangle: [[f64; 3]; 3]) -> bool {
    let radius = triangle
        .iter()
        .map(|&vertex| norm(sub(vertex, centroid)))
        .fold(0.0_f64, f64::max);
    norm(sub(x, centroid)) <= 4.0 * radius
}

/// Green rows at a well-separated receiver. Integrate the complete smooth
/// kernels here: subtracting an exact static potential at a very distant
/// receiver would needlessly introduce cancellation in its edge formula.
fn triangle_far_pressure_influence(
    k: f64,
    x: [f64; 3],
    triangle: [[f64; 3]; 3],
    normal: [f64; 3],
) -> (C64, C64) {
    let [a, b, c] = triangle;
    let ab = sub(b, a);
    let ac = sub(c, a);
    let double_area = norm(cross(ab, ac));
    let (mut single, mut double) = (C64::ZERO, C64::ZERO);
    for (&u, &wu) in GL4_X.iter().zip(&GL4_W) {
        for (&v, &wv) in GL4_X.iter().zip(&GL4_W) {
            let y = affine3(a, ab, u, ac, (1.0 - u) * v);
            let weight = wu * wv * (1.0 - u) * double_area;
            let delta = sub(x, y);
            let r = norm(delta);
            let g = green(k, r);
            let derivative = (g * C64::new(-1.0 / r, k)).scale(-dot(normal, delta) / r);
            single = single + g.scale(weight);
            double = double + derivative.scale(weight);
        }
    }
    (single, double)
}

/// Integral of exp(-i k direction.y) over one triangle. Pressure and normal
/// velocity are constant on a panel, but this directional phase is not.
fn triangle_far_phase(k: f64, direction: [f64; 3], triangle: [[f64; 3]; 3]) -> C64 {
    let [a, b, c] = triangle;
    let ab = sub(b, a);
    let ac = sub(c, a);
    let double_area = norm(cross(ab, ac));
    let mut integral = C64::ZERO;
    for (&u, &wu) in GL4_X.iter().zip(&GL4_W) {
        for (&v, &wv) in GL4_X.iter().zip(&GL4_W) {
            let y = affine3(a, ab, u, ac, (1.0 - u) * v);
            let weight = wu * wv * (1.0 - u) * double_area;
            integral = integral + expik(-k * dot(direction, y)).scale(weight);
        }
    }
    integral
}

fn characteristic_panel_size(surface: &SpherePanels) -> f64 {
    let mut max_area = 0.0f64;
    for &a in surface.areas() {
        max_area = max_area.max(a);
    }
    max_area.sqrt()
}

// Also used before a reusable geometry table allocates or integrates anything.
pub(crate) fn admit_resolution(surface: &SpherePanels, k: f64) -> Result<f64, HelmholtzError> {
    if !(k > 0.0 && k.is_finite()) {
        return Err(HelmholtzError::BadParameter {
            what: "wavenumber k must be positive and finite",
        });
    }
    let n = surface.centroids().len();
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
    Ok(ppw)
}

// Dense radiation-operator admission: frequency, medium, batch shape and
// resolution, in this order, before any assembly.
fn admit_operator(
    surface: &SpherePanels,
    k: f64,
    medium: Medium,
    velocity_fields: &[&[C64]],
) -> Result<f64, HelmholtzError> {
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
    if velocity_fields.is_empty() || velocity_fields.len() > MAX_RADIATION_FIELDS_PER_BATCH {
        return Err(HelmholtzError::BadParameter {
            what: "velocity-field batch must be nonempty and within MAX_RADIATION_FIELDS_PER_BATCH",
        });
    }
    if velocity_fields.iter().any(|velocity| velocity.len() != n) {
        return Err(HelmholtzError::ShapeMismatch {
            what: "velocity length must equal the panel count",
        });
    }
    admit_resolution(surface, k)
}

fn surface_fingerprint(surface: &SpherePanels) -> u64 {
    let mut hash = 0xcbf2_9ce4_8422_2325u64;
    hash ^= surface.centroids().len() as u64;
    for j in 0..surface.centroids().len() {
        for value in surface.centroids()[j]
            .into_iter()
            .chain(surface.normals()[j])
            .chain([surface.areas()[j]])
        {
            hash = (hash ^ value.to_bits()).wrapping_mul(0x0000_0100_0000_01b3);
        }
        if let Some(triangles) = surface.triangles() {
            for value in triangles[j].into_iter().flatten() {
                hash = (hash ^ value.to_bits()).wrapping_mul(0x0000_0100_0000_01b3);
            }
        }
    }
    hash
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
    let fields = [velocity];
    let operator = RadiationOperator::prepare(surface, k, medium, &fields, formulation)?;
    Ok(operator.solve(velocity))
}

/// Solve a nonempty batch of exterior Neumann radiation problems at one
/// frequency. All fields share the exact same validated Helmholtz matrix, LU
/// factorization, condition diagnostic, and panels-per-wavelength diagnostic;
/// returned solutions remain in input order and otherwise have exactly the
/// semantics of [`solve_radiation`].
///
/// # Errors
/// [`HelmholtzError`] on the same admissions as [`solve_radiation`], if any
/// velocity field has the wrong panel count, or if the batch is empty or
/// exceeds [`MAX_RADIATION_FIELDS_PER_BATCH`].
pub fn solve_radiation_batch(
    surface: &SpherePanels,
    k: f64,
    medium: Medium,
    velocity_fields: &[&[C64]],
    formulation: Formulation,
) -> Result<Vec<RadiationSolution>, HelmholtzError> {
    let operator = RadiationOperator::prepare(surface, k, medium, velocity_fields, formulation)?;
    Ok(operator.solve_all(velocity_fields))
}

/// One admitted and factorized exterior-radiation operator. Keeping this
/// private makes it impossible to mix a factorization with a different
/// surface, medium, frequency, or formulation.
struct RadiationOperator<'a> {
    surface: &'a SpherePanels,
    k: f64,
    medium: Medium,
    panels_per_wavelength: f64,
    bmat: Vec<C64>,
    lu: LuComplex,
    condition_lower_bound: f64,
}

impl<'a> RadiationOperator<'a> {
    fn prepare(
        surface: &'a SpherePanels,
        k: f64,
        medium: Medium,
        velocity_fields: &[&[C64]],
        formulation: Formulation,
    ) -> Result<Self, HelmholtzError> {
        Self::prepare_with_static(surface, k, medium, velocity_fields, formulation, None)
    }

    fn prepare_with_static(
        surface: &'a SpherePanels,
        k: f64,
        medium: Medium,
        velocity_fields: &[&[C64]],
        formulation: Formulation,
        static_layers: Option<&[(f64, f64)]>,
    ) -> Result<Self, HelmholtzError> {
        let ppw = admit_operator(surface, k, medium, velocity_fields)?;
        let n = surface.centroids().len();

        let alpha = alpha_for(formulation, k);
        let (matrix, bmat) = assemble_dense_with_static(surface, k, alpha, static_layers);
        let lu = lu_complex(&matrix, n).map_err(|_| HelmholtzError::Singular)?;
        let condition_lower_bound = condition_lower_bound(&matrix, n, &lu);
        Ok(Self {
            surface,
            k,
            medium,
            panels_per_wavelength: ppw,
            bmat,
            lu,
            condition_lower_bound,
        })
    }

    fn solve(&self, velocity: &[C64]) -> RadiationSolution {
        let n = self.surface.centroids().len();
        let omega_rho = self.k * self.medium.sound_speed * self.medium.density;
        // q = i omega rho v.
        let q: Vec<C64> = velocity
            .iter()
            .map(|&v| v * C64::new(0.0, omega_rho))
            .collect();
        let mut pressure = vec![C64::ZERO; n];
        for i in 0..n {
            let mut b = C64::ZERO;
            for (j, &qj) in q.iter().enumerate() {
                b = b + self.bmat[i * n + j] * qj;
            }
            pressure[i] = b;
        }
        self.lu.solve(&mut pressure);

        let (radiated_power, radiated_power_roundoff_interval, power_diagnostics) =
            evaluate_radiated_power(&pressure, velocity, self.surface.areas(), self.medium);

        RadiationSolution {
            pressure,
            velocity: velocity.to_vec(),
            k: self.k,
            medium: self.medium,
            surface_fingerprint: surface_fingerprint(self.surface),
            panels_per_wavelength: self.panels_per_wavelength,
            radiated_power,
            radiated_power_roundoff_interval,
            power_diagnostics,
            condition_lower_bound: self.condition_lower_bound,
            dense_cap_utilization: n as f64 / MAX_DENSE_PANELS as f64,
        }
    }

    /// Every field against the one factorization, in input order, on this
    /// thread's declared dense workers ([`with_dense_workers`]). Each field
    /// runs exactly [`Self::solve`], so the worker count never changes a bit.
    fn solve_all(&self, fields: &[&[C64]]) -> Vec<RadiationSolution> {
        let workers = DENSE_WORKERS.with(core::cell::Cell::get).min(fields.len());
        if workers <= 1 {
            return fields
                .iter()
                .map(|&velocity| self.solve(velocity))
                .collect();
        }
        let mut solved: Vec<Option<RadiationSolution>> = fields.iter().map(|_| None).collect();
        std::thread::scope(|scope| {
            let mut rest = solved.as_mut_slice();
            for p in 0..workers {
                let (start, end) = (p * fields.len() / workers, (p + 1) * fields.len() / workers);
                let (head, tail) = core::mem::take(&mut rest).split_at_mut(end - start);
                rest = tail;
                scope.spawn(move || {
                    for (slot, &velocity) in head.iter_mut().zip(&fields[start..end]) {
                        *slot = Some(self.solve(velocity));
                    }
                });
            }
        });
        solved
            .into_iter()
            .map(|s| s.expect("every scoped worker fills its own fields"))
            .collect()
    }
}

fn outward_product(a: f64, b: f64) -> (f64, f64) {
    let product = a * b;
    (product.next_down(), product.next_up())
}

fn outward_add(left: (f64, f64), right: (f64, f64)) -> (f64, f64) {
    ((left.0 + right.0).next_down(), (left.1 + right.1).next_up())
}

fn outward_scale_positive(value: (f64, f64), scale: f64) -> (f64, f64) {
    ((value.0 * scale).next_down(), (value.1 * scale).next_up())
}

fn evaluate_radiated_power(
    pressure: &[C64],
    velocity: &[C64],
    areas: &[f64],
    medium: Medium,
) -> (f64, (f64, f64), RadiationPowerDiagnostics) {
    let mut rounded = 0.0;
    let mut reactive = 0.0;
    let mut positive = 0.0;
    let mut negative = 0.0;
    let mut apparent = 0.0;
    let mut mean_square_velocity = 0.0;
    let mut enclosure = (0.0, 0.0);
    for ((&pressure, &velocity), &area) in pressure.iter().zip(velocity).zip(areas) {
        let dot = pressure * velocity.conj();
        let real_power = 0.5 * dot.re * area;
        rounded += real_power;
        reactive += 0.5 * dot.im * area;
        if real_power >= 0.0 {
            positive += real_power;
        } else {
            negative += real_power;
        }
        apparent += 0.5 * pressure.abs() * velocity.abs() * area;
        mean_square_velocity += velocity.norm_sq() * area;
        let real_dot = outward_add(
            outward_product(pressure.re, velocity.re),
            outward_product(pressure.im, velocity.im),
        );
        let panel_power = outward_scale_positive(outward_scale_positive(real_dot, 0.5), area);
        enclosure = outward_add(enclosure, panel_power);
    }
    let raw_real = rounded;
    let finite_enclosure = enclosure.0.is_finite() && enclosure.1.is_finite();
    if rounded < 0.0 && finite_enclosure && enclosure.0 <= 0.0 && enclosure.1 >= 0.0 {
        rounded = 0.0;
    }
    let diagnostics = RadiationPowerDiagnostics {
        surface_power: C64::new(raw_real, reactive),
        positive_real_power: positive,
        negative_real_power: negative,
        apparent_power: apparent,
        plane_wave_reference_power: 0.5
            * medium.density
            * medium.sound_speed
            * mean_square_velocity,
    };
    (rounded, enclosure, diagnostics)
}

/// Evaluate complex acoustic pressure [Pa] at finite exterior points [m].
/// `p(x) = SUM_j [dG(x,y_j)/dn_y p_j - G(x,y_j) q_j] A_j`, with
/// `q_j = i omega rho v_j`; output order matches `points` order. This entry
/// point retains the original centroid rule even when triangles are present.
/// See [`exterior_pressure_at_points_integrated`] for triangle integration.
/// Its discrete solid-angle guard is not a topology certificate.
///
/// # Errors
/// [`HelmholtzError`] for mismatched/malformed solution data or a non-finite,
/// centroid-coincident, or non-exterior point.
pub fn exterior_pressure_at_points(
    surface: &SpherePanels,
    solution: &RadiationSolution,
    medium: Medium,
    points: &[[f64; 3]],
) -> Result<Vec<C64>, HelmholtzError> {
    exterior_pressure_at_points_impl(surface, solution, medium, points, false)
}

/// Evaluate the same finite Green representation using retained triangles:
/// `p(x) = SUM_j INT_j [dG(x,y)/dn_y p_j - G(x,y) q_j] dS_y`.
/// Well-separated triangles use 4x4 Gauss integration; nearby triangles use
/// the static-subtracted 8x8 rule. Centroid-only surfaces keep their original
/// point-panel rule. Propagation and spreading are already included once.
///
/// This fixed quadrature is not an error certificate. [`crate::near_field`]
/// supplies the separately controlled close observer. Surface, solution and
/// receiver admission are identical to [`exterior_pressure_at_points`].
///
/// # Errors
/// The same malformed/cross-wired solution and point refusals as
/// [`exterior_pressure_at_points`].
pub fn exterior_pressure_at_points_integrated(
    surface: &SpherePanels,
    solution: &RadiationSolution,
    medium: Medium,
    points: &[[f64; 3]],
) -> Result<Vec<C64>, HelmholtzError> {
    exterior_pressure_at_points_impl(surface, solution, medium, points, true)
}

fn exterior_pressure_at_points_impl(
    surface: &SpherePanels,
    solution: &RadiationSolution,
    medium: Medium,
    points: &[[f64; 3]],
    integrate_triangles: bool,
) -> Result<Vec<C64>, HelmholtzError> {
    let n = surface.centroids().len();
    if solution.pressure.len() != n || solution.velocity.len() != n {
        return Err(HelmholtzError::ShapeMismatch {
            what: "solution fields must match the panel count",
        });
    }
    if solution.surface_fingerprint != surface_fingerprint(surface) {
        return Err(HelmholtzError::ShapeMismatch {
            what: "solution surface fingerprint must match the evaluation surface",
        });
    }
    if solution.medium.density.to_bits() != medium.density.to_bits()
        || solution.medium.sound_speed.to_bits() != medium.sound_speed.to_bits()
    {
        return Err(HelmholtzError::BadParameter {
            what: "evaluation medium must exactly match the solved medium",
        });
    }
    if !(medium.density > 0.0
        && medium.density.is_finite()
        && medium.sound_speed > 0.0
        && medium.sound_speed.is_finite()
        && solution.k > 0.0
        && solution.k.is_finite())
        || solution
            .pressure
            .iter()
            .chain(&solution.velocity)
            .any(|z| !z.re.is_finite() || !z.im.is_finite())
    {
        return Err(HelmholtzError::BadParameter {
            what: "radiation solution must contain finite fields and positive finite k",
        });
    }

    let k = solution.k;
    let omega_rho = k * medium.sound_speed * medium.density;
    let four_pi = 4.0 * core::f64::consts::PI;
    let triangles = if integrate_triangles {
        surface.triangles()
    } else {
        None
    };
    let mut pressure = Vec::with_capacity(points.len());
    for &x in points {
        if x.iter().any(|coordinate| !coordinate.is_finite()) {
            return Err(HelmholtzError::BadParameter {
                what: "evaluation points must be finite",
            });
        }
        let mut value = C64::ZERO;
        let mut solid_angle = 0.0;
        for j in 0..n {
            let d = sub(x, surface.centroids()[j]);
            let r = norm(d);
            if !(r > 0.0 && r.is_finite()) {
                return Err(HelmholtzError::BadParameter {
                    what: "evaluation points must not coincide with panel centroids",
                });
            }
            let normal_projection = dot(surface.normals()[j], d);
            solid_angle += normal_projection * surface.areas()[j] / (four_pi * r * r * r);
            let q = solution.velocity[j] * C64::new(0.0, omega_rho);
            if let Some(triangles) = triangles {
                let (single, double) = if triangle_is_near(x, surface.centroids()[j], triangles[j])
                {
                    triangle_weak_influence(k, x, triangles[j], surface.normals()[j])
                } else {
                    triangle_far_pressure_influence(k, x, triangles[j], surface.normals()[j])
                };
                value = value + (double * solution.pressure[j] - single * q);
            } else {
                let dgdny = green_dr(k, r).scale(-normal_projection / r);
                value = value
                    + (dgdny * solution.pressure[j] - green(k, r) * q).scale(surface.areas()[j]);
            }
        }
        if !solid_angle.is_finite()
            || solid_angle.abs() >= 0.5
            || !value.re.is_finite()
            || !value.im.is_finite()
        {
            return Err(HelmholtzError::BadParameter {
                what: "evaluation points must lie in the discrete exterior",
            });
        }
        pressure.push(value);
    }
    Ok(pressure)
}

/// Probe-based lower bound on `cond_1(A) = ||A||_1 ||A^{-1}||_1`:
/// `||A||_1` exactly from the assembled matrix, `||A^{-1}||_1` bounded
/// below by `||A^{-1} b||_1` over five deterministic canonical columns
/// (spread across the index range) plus the uniform vector, each with
/// `||b||_1 = 1`. Cost is O(n^2) per probe against the O(n^3)
/// factorization already paid.
fn condition_lower_bound(matrix: &[C64], n: usize, lu: &fs_la::eigen_complex::LuComplex) -> f64 {
    let mut a_norm1 = 0.0f64;
    for j in 0..n {
        let mut col = 0.0;
        for i in 0..n {
            col += matrix[i * n + j].abs();
        }
        a_norm1 = a_norm1.max(col);
    }
    let mut inv_norm1_lb = 0.0f64;
    let mut probe = |b: Vec<C64>| {
        let mut x = b;
        lu.solve(&mut x);
        let norm1: f64 = x.iter().map(|c| c.abs()).sum();
        inv_norm1_lb = inv_norm1_lb.max(norm1);
    };
    for t in 0..5usize {
        let j = if n == 1 { 0 } else { t * (n - 1) / 4 };
        let mut e = vec![C64::ZERO; n];
        e[j] = C64::ONE;
        probe(e);
    }
    probe(vec![C64::from_re(1.0 / n as f64); n]);
    a_norm1 * inv_norm1_lb
}

fn alpha_for(formulation: Formulation, k: f64) -> C64 {
    match formulation {
        Formulation::PlainCbie => C64::ZERO,
        Formulation::BurtonMiller => C64::new(0.0, 1.0 / k),
        Formulation::BurtonMillerWrongAlphaSign => C64::new(0.0, -1.0 / k),
    }
}

/// Assemble the dense Burton–Miller pair `(A, B)` with
/// `A p = B q`: row i is `(D - I/2 + alpha N)` against
/// `(S + alpha D' + alpha/2 I)`. The hypersingular self entry combines
/// the regularized disc part with the exact closed-surface row identity
/// `N_0[1] = 0` (see the module docs).
fn assemble_dense(surface: &SpherePanels, k: f64, alpha: C64) -> (Vec<C64>, Vec<C64>) {
    assemble_dense_with_static(surface, k, alpha, None)
}

/// Fill a row-major n x n pair row by row on this thread's declared dense
/// workers ([`with_dense_workers`]). `fill_row(i, a_row, b_row)` must write
/// row i from shared data only, so the worker count never changes a bit.
fn fill_rows(
    amat: &mut [C64],
    bmat: &mut [C64],
    n: usize,
    fill_row: &(impl Fn(usize, &mut [C64], &mut [C64]) + Sync),
) {
    let workers = DENSE_WORKERS.with(core::cell::Cell::get).min(n);
    if workers <= 1 {
        for (i, (arow, brow)) in amat
            .chunks_exact_mut(n.max(1))
            .zip(bmat.chunks_exact_mut(n.max(1)))
            .enumerate()
        {
            fill_row(i, arow, brow);
        }
    } else {
        let (mut a_rest, mut b_rest) = (amat, bmat);
        std::thread::scope(|scope| {
            for p in 0..workers {
                let (start, end) = (p * n / workers, (p + 1) * n / workers);
                let entries = (end - start) * n;
                let (a_head, a_tail) = core::mem::take(&mut a_rest).split_at_mut(entries);
                let (b_head, b_tail) = core::mem::take(&mut b_rest).split_at_mut(entries);
                (a_rest, b_rest) = (a_tail, b_tail);
                scope.spawn(move || {
                    for (r, (arow, brow)) in a_head
                        .chunks_exact_mut(n)
                        .zip(b_head.chunks_exact_mut(n))
                        .enumerate()
                    {
                        fill_row(start + r, arow, brow);
                    }
                });
            }
        });
    }
}

fn assemble_dense_with_static<K: WaveNumber>(
    surface: &SpherePanels,
    k: K,
    alpha: C64,
    static_layers: Option<&[(f64, f64)]>,
) -> (Vec<C64>, Vec<C64>) {
    let centroids = surface.centroids();
    let normals = surface.normals();
    let areas = surface.areas();
    let triangles = surface.triangles();
    // The triangle rule below is sound for the weakly singular single-layer
    // and ordinary double-layer operators used by Plain CBIE.  Burton-Miller
    // also contains a finite-part hypersingular operator; applying an ordinary
    // area rule to its edge-near panels changes that operator and regresses the
    // analytic pulsating-sphere impedance.  Preserve the established
    // equivalent-disc/centroid discretization for that arm until a genuine
    // finite-part Galerkin rule is implemented.
    let weak_triangles = if alpha == C64::ZERO { triangles } else { None };
    let n = centroids.len();
    debug_assert!(static_layers.is_none_or(|layers| layers.len() == n * n));
    let half = C64::new(0.5, 0.0);
    let mut amat = vec![C64::ZERO; n * n];
    let mut bmat = vec![C64::ZERO; n * n];
    // Each row reads only shared geometry and writes only its own entries, so
    // rows split across the declared workers give the same bits.
    let fill_row = |i: usize, arow: &mut [C64], brow: &mut [C64]| {
        let xi = centroids[i];
        let ni = normals[i];
        // Static hypersingular row sum: N_0[1] = 0 on a closed surface,
        // so the self entry is minus the off-diagonal static sum under
        // the SAME point-panel quadrature.
        let mut n0_row = C64::ZERO;
        let mut d0_row = C64::ZERO;
        if alpha != C64::ZERO {
            for j in 0..n {
                if j != i {
                    let (_, d0, _, n0) = kernels(0.0, xi, ni, centroids[j], normals[j]);
                    n0_row = n0_row + n0.scale(areas[j]);
                    d0_row = d0_row + d0.scale(areas[j]);
                }
            }
        }
        for j in 0..n {
            let (s_ij, d_ij, dp_ij, n_ij) = if i == j {
                let (s, n_reg) = if let Some(triangles) = weak_triangles {
                    let single = if let Some(layers) = static_layers {
                        triangle_self_terms_with_static(
                            k.cbie_real(),
                            xi,
                            triangles[i],
                            layers[i * n + i].0,
                        )
                        .0
                    } else {
                        triangle_self_terms(k.cbie_real(), xi, triangles[i]).0
                    };
                    (single, C64::ZERO)
                } else {
                    self_terms(k, areas[i])
                };
                // The SAME static quadrature must also satisfy (D_0+I/2)1=0.
                // Otherwise its weighted adjoint D'_0 has a spurious net flux:
                // A^T(D'_0+I/2) != 0, while A^T N_0 = 0. Multiplication by
                // alpha (|alpha| = 1/k) then amplifies this incompatible RHS at low k,
                // creating an O(panel-size) real pressure at k -> 0.
                // This geometric Gauss correction closes the static identity,
                // just as the N_0 row correction does. It uses no computed
                // resistance or pressure and preserves the outgoing kernel.
                let d_self = if alpha == C64::ZERO {
                    C64::ZERO
                } else {
                    -half - d0_row
                };
                (s, d_self, d_self, n_reg - n0_row)
            } else if let Some(triangles) = weak_triangles {
                // Static S/D must use the same triangle geometry everywhere.
                // Integrate the regular far remainder too: centroid moments
                // can swamp the radiation of nearly cancelling skin velocities.
                let (s, d) = if triangle_is_near(xi, centroids[j], triangles[j]) {
                    if let Some(layers) = static_layers {
                        triangle_weak_influence_with_static(
                            k.cbie_real(),
                            xi,
                            triangles[j],
                            normals[j],
                            (&GL8_X, &GL8_W),
                            layers[i * n + j],
                        )
                    } else {
                        triangle_weak_influence(k.cbie_real(), xi, triangles[j], normals[j])
                    }
                } else if let Some(layers) = static_layers {
                    triangle_weak_influence_with_static(
                        k.cbie_real(),
                        xi,
                        triangles[j],
                        normals[j],
                        (&GL4_X, &GL4_W),
                        layers[i * n + j],
                    )
                } else {
                    triangle_weak_influence_with_rule(
                        k.cbie_real(),
                        xi,
                        triangles[j],
                        normals[j],
                        &GL4_X,
                        &GL4_W,
                    )
                };
                (s, d, C64::ZERO, C64::ZERO)
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
            arow[j] = aij;
            brow[j] = bij;
        }
    };
    fill_rows(&mut amat, &mut bmat, n, &fill_row);
    (amat, bmat)
}

/// Far-field directivity amplitude `F(direction)`: the radiated pressure
/// behaves as `p -> F * e^{ikr}/r` for large r, with
/// `F = (1/4 pi) SUM_j [(-ik dir.n_j) p_j - q_j] e^{-ik dir.y_j} A_j`.
/// This entry point retains the original centroid phase even when triangles
/// are present. See [`far_field_integrated`] for triangle integration.
#[must_use]
pub fn far_field(
    surface: &SpherePanels,
    solution: &RadiationSolution,
    medium: Medium,
    directions: &[[f64; 3]],
) -> Vec<C64> {
    far_field_impl(surface, solution, medium, directions, false)
}

/// Far-field amplitude with the same traces and convention as [`far_field`],
/// integrating `exp(-i k direction.y)` over every retained triangle with a
/// 4x4 Gauss rule. Centroid-only surfaces keep the original point-panel rule.
/// The result still requires `exp(i k r)/r` exactly once to obtain pressure;
/// this does not replace a far-field receiver by a finite or near receiver.
#[must_use]
pub fn far_field_integrated(
    surface: &SpherePanels,
    solution: &RadiationSolution,
    medium: Medium,
    directions: &[[f64; 3]],
) -> Vec<C64> {
    far_field_impl(surface, solution, medium, directions, true)
}

fn far_field_impl(
    surface: &SpherePanels,
    solution: &RadiationSolution,
    medium: Medium,
    directions: &[[f64; 3]],
    integrate_triangles: bool,
) -> Vec<C64> {
    let k = solution.k;
    let omega_rho = k * medium.sound_speed * medium.density;
    let centroids = surface.centroids();
    let normals = surface.normals();
    let areas = surface.areas();
    let triangles = if integrate_triangles {
        surface.triangles()
    } else {
        None
    };
    directions
        .iter()
        .map(|&dir| {
            let len = norm(dir);
            let d = [dir[0] / len, dir[1] / len, dir[2] / len];
            let mut f = C64::ZERO;
            for j in 0..centroids.len() {
                let qj = solution.velocity[j] * C64::new(0.0, omega_rho);
                let term = solution.pressure[j] * C64::new(0.0, -k * dot(d, normals[j])) - qj;
                let contribution = if let Some(triangles) = triangles {
                    term * triangle_far_phase(k, d, triangles[j])
                } else {
                    (term * expik(-k * dot(d, centroids[j]))).scale(areas[j])
                };
                f = f + contribution;
            }
            f.scale(1.0 / (4.0 * core::f64::consts::PI))
        })
        .collect()
}

/// Solver-neutral far-field table: `p → F(ω, dir) e^{ikr}/r`.
///
/// Couple (and any other L3 observer) consumes these samples without
/// taking a production dependency on this crate.
#[derive(Debug, Clone, PartialEq)]
pub struct FarFieldTable {
    /// Producer identity retained by the consumer.
    pub source_id: String,
    /// Angular frequencies [rad/s], strictly increasing.
    pub omegas: Vec<f64>,
    /// Unit listener directions.
    pub directions: Vec<[f64; 3]>,
    /// Row-major amplitudes `F[i_ω, i_dir]`.
    pub amplitudes: Vec<C64>,
    /// Sound speed used to form `k = ω/c` [m/s].
    pub sound_speed: f64,
    /// Ambient density [kg/m³].
    pub density: f64,
}

impl FarFieldTable {
    /// Amplitude `F` at frequency `i` and direction `j`.
    #[must_use]
    pub fn amplitude(&self, i_omega: usize, i_dir: usize) -> Option<C64> {
        if i_omega >= self.omegas.len() || i_dir >= self.directions.len() {
            return None;
        }
        self.amplitudes
            .get(i_omega * self.directions.len() + i_dir)
            .copied()
    }
}

/// Tabulate [`far_field`] on an explicit `(ω, direction)` grid.
///
/// # Errors
/// [`HelmholtzError`] from [`solve_radiation`], or a bad frequency /
/// direction / velocity admission.
pub fn tabulate_far_field(
    surface: &SpherePanels,
    omegas: &[f64],
    directions: &[[f64; 3]],
    medium: Medium,
    velocity: &[C64],
    formulation: Formulation,
) -> Result<FarFieldTable, HelmholtzError> {
    if omegas.is_empty() || directions.is_empty() {
        return Err(HelmholtzError::BadParameter {
            what: "far-field table needs at least one frequency and direction",
        });
    }
    let mut last = 0.0;
    for (i, &w) in omegas.iter().enumerate() {
        if !(w > 0.0 && w.is_finite()) {
            return Err(HelmholtzError::BadParameter {
                what: "far-field frequencies must be finite and positive",
            });
        }
        if i > 0 && w <= last {
            return Err(HelmholtzError::BadParameter {
                what: "far-field frequencies must be strictly increasing",
            });
        }
        last = w;
    }
    for dir in directions {
        let n2 = dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2];
        if !(n2 > 0.0 && n2.is_finite()) {
            return Err(HelmholtzError::BadParameter {
                what: "far-field directions must be finite and nonzero",
            });
        }
    }
    let mut amplitudes = Vec::with_capacity(omegas.len() * directions.len());
    for &omega in omegas {
        let k = omega / medium.sound_speed;
        let sol = solve_radiation(surface, k, medium, velocity, formulation)?;
        amplitudes.extend(far_field(surface, &sol, medium, directions));
    }
    Ok(FarFieldTable {
        source_id: "fs-bem::helmholtz::tabulate_far_field".to_string(),
        omegas: omegas.to_vec(),
        directions: directions.to_vec(),
        amplitudes,
        sound_speed: medium.sound_speed,
        density: medium.density,
    })
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
    // Share admission without assembling/factoring a zero-source problem;
    // the one factorization below serves every coherent source column.
    let probe = vec![C64::ZERO; n];
    admit_operator(surface, k, medium, &[&probe])?;
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

/// Mean piston-face impedance `p/v` [Pa s/m] in an infinite rigid baffle.
/// Converts the shared Bessel–Struve acoustic `p/U` owner by the face area.
///
/// # Errors
/// [`HelmholtzError::BadParameter`] on non-positive radius, wavenumber,
/// medium, or unrepresentable derived impedance.
pub fn baffled_piston_face_impedance(
    radius: f64,
    k: f64,
    medium: Medium,
) -> Result<C64, HelmholtzError> {
    let (r, x) = fs_phs::baffled_piston_impedance(
        medium.density,
        medium.sound_speed,
        radius,
        k * medium.sound_speed,
    )
    .map_err(|_| HelmholtzError::BadParameter {
        what: "baffled piston parameters or derived scale",
    })?;
    let result = C64::new(r, x).scale(core::f64::consts::PI * radius * radius);
    if !result.re.is_finite() || !result.im.is_finite() {
        return Err(HelmholtzError::BadParameter {
            what: "baffled piston face impedance overflow",
        });
    }
    Ok(result)
}

/// One historical polar-cell rule, retained as a G3 discretization comparison.
/// The former bem, phs and thin_plate copies had the same cells and kernel;
/// their equivalent-disc self terms and p/U conversion were algebraically
/// equivalent. They are one algorithm, not three independent oracles.
#[cfg(test)]
fn historical_ring_piston(
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

/// Cap on the spherical-harmonic table degree: the projection grid is
/// `(l_max + 1) x (2 l_max + 1)` far-field evaluations and the
/// normalized-recurrence Legendre path is validated to this degree.
pub const MAX_SH_DEGREE: usize = 64;

/// Far-field directivity as spherical-harmonic coefficients: the
/// rendering-facing table promised to the musical-acoustics program
/// (runtime stereo/space rendering evaluates `F(direction)` from
/// `(l_max + 1)^2` complex numbers instead of resampling the BEM
/// far field).
///
/// Convention: orthonormal complex harmonics with Condon–Shortley
/// phase, `Y_lm(theta, phi) = Pbar_l^m(cos theta) e^{i m phi}` for
/// `m >= 0` and `Y_{l,-m} = (-1)^m conj(Y_lm)`, where `Pbar` carries
/// the full `sqrt((2l+1)/(4 pi) (l-m)!/(l+m)!)` normalization inside
/// the recurrence (no factorial overflow). `F(dir) = SUM a_lm Y_lm`.
#[derive(Debug, Clone)]
pub struct DirectivityTable {
    /// Wavenumber the table was built at [1/m].
    pub k: f64,
    /// Maximum spherical-harmonic degree l.
    pub l_max: usize,
    /// Coefficients `a_lm`, indexed `l (l + 1) + m` for `m` in
    /// `-l ..= l` — length `(l_max + 1)^2`.
    pub coefficients: Vec<C64>,
    /// `SUM |a_lm|^2` over the quadrature estimate of
    /// `INT |F|^2 dOmega`: the fraction of far-field power the
    /// truncated table captures. Near 1 for mesh-resolved fields; a low
    /// value means the truncation (not the BEM solve) is the limiting
    /// approximation.
    pub captured_fraction: f64,
}

impl DirectivityTable {
    /// The coefficient `a_lm`.
    ///
    /// # Panics
    /// If `l > l_max` or `|m| > l` (programmer error, not a refusal).
    #[must_use]
    pub fn coefficient(&self, l: usize, m: i64) -> C64 {
        let m_abs = usize::try_from(m.unsigned_abs()).expect("|m| fits in usize");
        assert!(l <= self.l_max && m_abs <= l);
        // Index l (l + 1) + m: for the asserted |m| <= l this is
        // l^2 + (l + m) with both terms non-negative.
        let offset = usize::try_from(i64::try_from(l).expect("l bounded by MAX_SH_DEGREE") + m)
            .expect("l + m is non-negative for |m| <= l");
        self.coefficients[l * l + offset]
    }

    /// Per-degree power `SUM_m |a_lm|^2` — the multipole spectrum.
    #[must_use]
    pub fn power_by_degree(&self) -> Vec<f64> {
        let mut out = vec![0.0; self.l_max + 1];
        for (l, slot) in out.iter_mut().enumerate() {
            let lm = i64::try_from(l).expect("l bounded by MAX_SH_DEGREE");
            for m in -lm..=lm {
                *slot += self.coefficient(l, m).norm_sq();
            }
        }
        out
    }

    /// Evaluate the table at a direction (need not be normalized):
    /// the far-field amplitude `F` with `p -> F e^{ikr}/r`.
    ///
    /// # Panics
    /// If the direction is the zero vector or non-finite (programmer
    /// error, not a refusal).
    #[must_use]
    pub fn evaluate(&self, direction: [f64; 3]) -> C64 {
        let len = norm(direction);
        assert!(
            len > 0.0 && len.is_finite(),
            "direction must be a nonzero finite vector"
        );
        let x = direction[2] / len;
        let phi = det::atan2(direction[1], direction[0]);
        let pbar = norm_assoc_legendre(self.l_max, x);
        let idx = |l: usize, m: usize| l * (l + 1) / 2 + m;
        let mut f = C64::ZERO;
        for l in 0..=self.l_max {
            // m = 0 term.
            f = f + self.coefficient(l, 0).scale(pbar[idx(l, 0)]);
            for m in 1..=l {
                let e = expik(m as f64 * phi);
                let y_pos = e.scale(pbar[idx(l, m)]);
                // Y_{l,-m} = (-1)^m conj(Y_lm).
                let sign = if m % 2 == 0 { 1.0 } else { -1.0 };
                let y_neg = y_pos.conj().scale(sign);
                let ms = i64::try_from(m).expect("m bounded by MAX_SH_DEGREE");
                f = f + self.coefficient(l, ms) * y_pos + self.coefficient(l, -ms) * y_neg;
            }
        }
        f
    }
}

/// Project the BEM far field onto spherical harmonics up to `l_max`.
/// Quadrature is exact for band-limited integrands: Gauss–Legendre in
/// `cos theta` with `l_max + 1` nodes (degree `2 l_max + 1`
/// polynomials) times a `2 l_max + 1`-point uniform trapezoid in `phi`
/// (exact for `e^{i m phi}`, `|m| <= 2 l_max`); the residual truncation
/// error is reported as `1 - captured_fraction`, not hidden.
///
/// # Errors
/// [`HelmholtzError::BadParameter`] when `l_max` exceeds
/// [`MAX_SH_DEGREE`].
pub fn directivity_sh_table(
    surface: &SpherePanels,
    solution: &RadiationSolution,
    medium: Medium,
    l_max: usize,
) -> Result<DirectivityTable, HelmholtzError> {
    directivity_sh_table_from_evaluator(solution.k, l_max, |directions| {
        Ok(far_field(surface, solution, medium, directions))
    })
}

/// Project any caller-supplied far-field evaluator onto the canonical
/// spherical-harmonic grid used by [`directivity_sh_table`].
///
/// This is the representation-neutral radiation seam: boundary-integral,
/// Rayleigh-integral, analytic, and future fast-multipole producers can all
/// publish the same checked [`DirectivityTable`] without pretending to be one
/// another. The evaluator must return one finite complex amplitude for every
/// requested unit direction, using the convention `p -> F(dir) exp(i k r)/r`.
///
/// # Errors
/// [`HelmholtzError::BadParameter`] if `k` is not finite and positive, the
/// degree exceeds [`MAX_SH_DEGREE`], or the evaluator returns the wrong number
/// of samples or a non-finite amplitude. Evaluator refusals propagate.
pub fn directivity_sh_table_from_evaluator(
    k: f64,
    l_max: usize,
    evaluator: impl FnOnce(&[[f64; 3]]) -> Result<Vec<C64>, HelmholtzError>,
) -> Result<DirectivityTable, HelmholtzError> {
    if !(k.is_finite() && k > 0.0) {
        return Err(HelmholtzError::BadParameter {
            what: "directivity wavenumber must be finite and positive",
        });
    }
    if l_max > MAX_SH_DEGREE {
        return Err(HelmholtzError::BadParameter {
            what: "spherical-harmonic degree exceeds MAX_SH_DEGREE",
        });
    }
    let n_theta = l_max + 1;
    let n_phi = 2 * l_max + 1;
    let (nodes, weights) = gauss_legendre(n_theta);
    let mut dirs = Vec::with_capacity(n_theta * n_phi);
    for &x in &nodes {
        let s = (1.0 - x * x).max(0.0).sqrt();
        for p in 0..n_phi {
            let phi = 2.0 * core::f64::consts::PI * p as f64 / n_phi as f64;
            dirs.push([s * det::cos(phi), s * det::sin(phi), x]);
        }
    }
    let f = evaluator(&dirs)?;
    if f.len() != dirs.len()
        || f.iter()
            .any(|sample| !sample.re.is_finite() || !sample.im.is_finite())
    {
        return Err(HelmholtzError::BadParameter {
            what: "far-field evaluator returned malformed directivity samples",
        });
    }
    let w_phi = 2.0 * core::f64::consts::PI / n_phi as f64;
    let mut coefficients = vec![C64::ZERO; (l_max + 1) * (l_max + 1)];
    let mut quadrature_power = 0.0f64;
    let idx = |l: usize, m: usize| l * (l + 1) / 2 + m;
    for (ti, &x) in nodes.iter().enumerate() {
        let pbar = norm_assoc_legendre(l_max, x);
        for p in 0..n_phi {
            let w = weights[ti] * w_phi;
            let fval = f[ti * n_phi + p];
            quadrature_power += fval.norm_sq() * w;
            let phi = 2.0 * core::f64::consts::PI * p as f64 / n_phi as f64;
            for l in 0..=l_max {
                let base = l * (l + 1);
                coefficients[base] = coefficients[base] + fval.scale(pbar[idx(l, 0)] * w);
            }
            for m in 1..=l_max {
                let e = expik(m as f64 * phi);
                let sign = if m % 2 == 0 { 1.0 } else { -1.0 };
                for l in m..=l_max {
                    let y_pos = e.scale(pbar[idx(l, m)]);
                    let y_neg = y_pos.conj().scale(sign);
                    let base = l * (l + 1);
                    coefficients[base + m] =
                        coefficients[base + m] + (fval * y_pos.conj()).scale(w);
                    coefficients[base - m] =
                        coefficients[base - m] + (fval * y_neg.conj()).scale(w);
                }
            }
        }
    }
    let table_power: f64 = coefficients.iter().map(|c| c.norm_sq()).sum();
    let captured_fraction = if quadrature_power > 0.0 {
        table_power / quadrature_power
    } else {
        0.0
    };
    Ok(DirectivityTable {
        k,
        l_max,
        coefficients,
        captured_fraction,
    })
}

/// Radiation efficiency of the solved velocity pattern:
/// `sigma = W / (1/2 rho c INT |v_n|^2 dS)` — the radiated-power ratio
/// per structural mode promised to the vibroacoustic-coupling bead.
/// Oracles: pulsating sphere `sigma = (ka)^2 / (1 + (ka)^2)`,
/// oscillating (dipole) sphere `sigma = (ka)^4 / (4 + (ka)^4)`.
/// The value inherits the accuracy boundary of the solve arm: with
/// `BurtonMiller` below ka ~ 0.5 the radiated power (and hence sigma)
/// can even go negative — the documented centroid-quadrature
/// resistance artifact; use `PlainCbie` below the first interior
/// resonance.
///
/// # Errors
/// [`HelmholtzError::BadParameter`] when the stored velocity field has
/// zero or non-finite mean-square (no efficiency is defined).
pub fn radiation_efficiency(
    surface: &SpherePanels,
    solution: &RadiationSolution,
    medium: Medium,
) -> Result<f64, HelmholtzError> {
    let areas = surface.areas();
    if solution.velocity.len() != areas.len() {
        return Err(HelmholtzError::ShapeMismatch {
            what: "solution velocity length must equal the panel count",
        });
    }
    let mut msv = 0.0f64;
    for (v, &a) in solution.velocity.iter().zip(areas.iter()) {
        msv += v.norm_sq() * a;
    }
    if !(msv > 0.0 && msv.is_finite()) {
        return Err(HelmholtzError::BadParameter {
            what: "radiation efficiency needs a nonzero finite velocity field",
        });
    }
    Ok(solution.radiated_power / (0.5 * medium.density * medium.sound_speed * msv))
}

/// Legendre `P_n(x)` and its derivative by the three-term recurrence.
fn legendre_pair(n: usize, x: f64) -> (f64, f64) {
    if n == 0 {
        return (1.0, 0.0);
    }
    let mut p0 = 1.0f64;
    let mut p1 = x;
    for l in 2..=n {
        let lf = l as f64;
        let p2 = ((2.0 * lf - 1.0) * x * p1 - (lf - 1.0) * p0) / lf;
        p0 = p1;
        p1 = p2;
    }
    let dp = n as f64 * (x * p1 - p0) / (x * x - 1.0);
    (p1, dp)
}

/// Gauss–Legendre nodes and weights on [-1, 1]: cosine initial guesses
/// refined by a fixed number of Newton steps (deterministic; converges
/// to machine precision well inside the budget for every degree the
/// [`MAX_SH_DEGREE`] cap admits).
fn gauss_legendre(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut nodes = vec![0.0; n];
    let mut weights = vec![0.0; n];
    for i in 0..n {
        let mut x = det::cos(core::f64::consts::PI * (i as f64 + 0.75) / (n as f64 + 0.5));
        let mut dp = 1.0;
        for _ in 0..16 {
            let (p, d) = legendre_pair(n, x);
            dp = d;
            let step = p / d;
            x -= step;
            if step.abs() < 1e-15 {
                break;
            }
        }
        nodes[i] = x;
        weights[i] = 2.0 / ((1.0 - x * x) * dp * dp);
    }
    (nodes, weights)
}

/// Orthonormal (Condon–Shortley) associated Legendre values
/// `Pbar_l^m(x)` for all `0 <= m <= l <= l_max`, packed at
/// `l (l + 1) / 2 + m`. The `4 pi` normalization lives inside the
/// recurrence, so `INT |Pbar_l^m(cos theta) e^{i m phi}|^2 dOmega = 1`
/// and no factorial ratio is ever formed.
fn norm_assoc_legendre(l_max: usize, x: f64) -> Vec<f64> {
    let s = (1.0 - x * x).max(0.0).sqrt();
    let count = (l_max + 1) * (l_max + 2) / 2;
    let mut p = vec![0.0f64; count];
    let idx = |l: usize, m: usize| l * (l + 1) / 2 + m;
    p[0] = 0.5 / core::f64::consts::PI.sqrt();
    for m in 1..=l_max {
        let mf = m as f64;
        p[idx(m, m)] = -((2.0 * mf + 1.0) / (2.0 * mf)).sqrt() * s * p[idx(m - 1, m - 1)];
    }
    for m in 0..l_max {
        let mf = m as f64;
        p[idx(m + 1, m)] = (2.0 * mf + 3.0).sqrt() * x * p[idx(m, m)];
    }
    for m in 0..=l_max {
        for l in (m + 2)..=l_max {
            let lf = l as f64;
            let mf = m as f64;
            let a = ((4.0 * lf * lf - 1.0) / (lf * lf - mf * mf)).sqrt();
            let b = (((lf - 1.0) * (lf - 1.0) - mf * mf) / (4.0 * (lf - 1.0) * (lf - 1.0) - 1.0))
                .sqrt();
            p[idx(l, m)] = a * (x * p[idx(l - 1, m)] - b * p[idx(l - 2, m)]);
        }
    }
    p
}

#[cfg(test)]
mod tests {
    #[test]
    fn g3_resolved_medium_preserves_density_and_sound_speed_bits() {
        let mut gas =
            super::GasState::try_new(&super::GasSpec::dry_air_ussa1976(), 308.15, 90000.0).unwrap();
        gas.density = f64::from_bits(gas.density.to_bits() + 1);
        gas.sound_speed = f64::from_bits(gas.sound_speed.to_bits() + 1);
        let medium = super::Medium::from(&gas);
        assert_eq!(medium.density.to_bits(), gas.density.to_bits());
        assert_eq!(medium.sound_speed.to_bits(), gas.sound_speed.to_bits());
    }

    use super::*;

    #[test]
    fn g3_checked_material_acoustics_preserve_gas_and_reject_stale_state() {
        let mut gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100_000.0).unwrap();
        let medium = Medium::from_material(&gas).unwrap();
        assert_eq!(medium.density.to_bits(), gas.density.to_bits());
        assert_eq!(medium.sound_speed.to_bits(), gas.sound_speed.to_bits());
        gas.sound_speed *= 1.01;
        assert!(Medium::from_material(&gas).is_err());
    }

    #[test]
    fn g1_g3_material_acoustics_drive_liquid_radiation_without_a_gas_law() {
        use fs_material::fluid_acoustics::FluidAcousticState;
        // Manufactured liquid coefficients, not measured material data.
        let liquid = FluidAcousticState::from_adiabatic_bulk_modulus(1000.0, 2.25e9).unwrap();
        let medium = Medium::from_material(&liquid).unwrap();
        assert_eq!(medium.sound_speed, 1500.0);
        let surface = SpherePanels::icosphere(1.0, 2).unwrap();
        let velocity = uniform_velocity(surface.centroids().len());
        let k = 0.2;
        let solution =
            solve_radiation(&surface, k, medium, &velocity, Formulation::BurtonMiller).unwrap();
        let expected =
            C64::new(k * k, -k).scale(medium.density * medium.sound_speed / (1.0 + k * k));
        assert!((mean_impedance(&solution, &surface) - expected).abs() / expected.abs() < 0.04);
        let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100_000.0).unwrap();
        let gas_medium = Medium::from_material(&gas).unwrap();
        let reference = solve_radiation(
            &surface,
            k,
            gas_medium,
            &velocity,
            Formulation::BurtonMiller,
        )
        .unwrap();
        let ratio =
            medium.density * medium.sound_speed / (gas_medium.density * gas_medium.sound_speed);
        assert!((solution.radiated_power / reference.radiated_power / ratio - 1.0).abs() < 1e-10);
    }

    fn rho_c() -> f64 {
        let medium = Medium::air();
        medium.density * medium.sound_speed
    }

    /// Analytic pulsating-sphere surface impedance under e^{-i omega t}:
    /// `z = rho c (k^2 R^2 - i k R) / (1 + k^2 R^2)`.
    fn pulsating_sphere_impedance(k: f64, radius: f64) -> C64 {
        let ka = k * radius;
        C64::new(ka * ka, -ka).scale(rho_c() / (1.0 + ka * ka))
    }

    #[test]
    fn medium_air_is_bitwise_the_declared_room_temperature_gas_state() {
        let gas = GasState::try_new(
            &GasSpec::dry_air_ussa1976(),
            ROOM_TEMPERATURE_K,
            STANDARD_ATMOSPHERE_PA,
        )
        .expect("declared dry-air room state");
        let medium = Medium::air();

        assert_eq!(gas.temperature.to_bits(), ROOM_TEMPERATURE_K.to_bits());
        assert_eq!(gas.pressure.to_bits(), STANDARD_ATMOSPHERE_PA.to_bits());
        assert_eq!(medium.density.to_bits(), gas.density.to_bits());
        assert_eq!(medium.sound_speed.to_bits(), gas.sound_speed.to_bits());
    }

    fn uniform_velocity(n: usize) -> Vec<C64> {
        vec![C64::ONE; n]
    }

    fn thin_box_surface() -> SpherePanels {
        // A closed plate-sized boundary with opposing skins 1 mm apart.
        // Six divisions across each skin leave both near and far interactions;
        // the coarse wavelength guard still admits the 1640 Hz endpoint.
        let extents = [0.30, 0.24, 0.001];
        let divisions = [6, 6, 1];
        let mut triangles = Vec::new();
        for axis in 0..3 {
            let u = (axis + 1) % 3;
            let v = (axis + 2) % 3;
            for side in [-1.0, 1.0] {
                let point = |i: usize, j: usize| {
                    let mut p = [0.0; 3];
                    p[axis] = side * extents[axis] / 2.0;
                    p[u] = (i as f64 / divisions[u] as f64 - 0.5) * extents[u];
                    p[v] = (j as f64 / divisions[v] as f64 - 0.5) * extents[v];
                    p
                };
                for i in 0..divisions[u] {
                    for j in 0..divisions[v] {
                        let a = point(i, j);
                        let b = point(i + 1, j);
                        let c = point(i + 1, j + 1);
                        let d = point(i, j + 1);
                        for mut triangle in [[a, b, c], [a, c, d]] {
                            if side < 0.0 {
                                triangle.swap(1, 2);
                            }
                            triangles.push(triangle);
                        }
                    }
                }
            }
        }
        SpherePanels::from_triangles(triangles).expect("closed thin box")
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

    fn assert_solution_bitwise_equal(actual: &RadiationSolution, expected: &RadiationSolution) {
        let assert_complex_bits = |actual: C64, expected: C64| {
            assert_eq!(actual.re.to_bits(), expected.re.to_bits());
            assert_eq!(actual.im.to_bits(), expected.im.to_bits());
        };
        assert_eq!(actual.pressure.len(), expected.pressure.len());
        for (&actual, &expected) in actual.pressure.iter().zip(expected.pressure.iter()) {
            assert_complex_bits(actual, expected);
        }
        assert_eq!(actual.velocity.len(), expected.velocity.len());
        for (&actual, &expected) in actual.velocity.iter().zip(expected.velocity.iter()) {
            assert_complex_bits(actual, expected);
        }
        assert_eq!(actual.k.to_bits(), expected.k.to_bits());
        assert_eq!(actual.medium, expected.medium);
        assert_eq!(actual.surface_fingerprint, expected.surface_fingerprint);
        assert_eq!(
            actual.panels_per_wavelength.to_bits(),
            expected.panels_per_wavelength.to_bits()
        );
        assert_eq!(
            actual.radiated_power.to_bits(),
            expected.radiated_power.to_bits()
        );
        assert_eq!(
            actual.radiated_power_roundoff_interval.0.to_bits(),
            expected.radiated_power_roundoff_interval.0.to_bits()
        );
        assert_eq!(
            actual.radiated_power_roundoff_interval.1.to_bits(),
            expected.radiated_power_roundoff_interval.1.to_bits()
        );
        assert_eq!(actual.power_diagnostics, expected.power_diagnostics);
        assert_eq!(
            actual.condition_lower_bound.to_bits(),
            expected.condition_lower_bound.to_bits()
        );
        assert_eq!(
            actual.dense_cap_utilization.to_bits(),
            expected.dense_cap_utilization.to_bits()
        );
    }

    #[test]
    fn radiated_power_roundoff_clamp_is_scale_aware_and_fail_closed() {
        let cancelling_pressure = [
            C64::from_re(-2.0e16),
            C64::from_re(2.0),
            C64::from_re(2.0e16),
            C64::from_re(-2.0),
        ];
        let velocity = [C64::ONE; 4];
        let areas = [1.0; 4];
        let medium = Medium::air();
        let (neutral, neutral_interval, diagnostics) =
            evaluate_radiated_power(&cancelling_pressure, &velocity, &areas, medium);
        assert_eq!(neutral.to_bits(), 0.0_f64.to_bits());
        assert!(neutral_interval.0 <= 0.0 && neutral_interval.1 >= 0.0);
        assert_eq!(diagnostics.surface_power, C64::from_re(-1.0));
        assert_eq!(diagnostics.positive_real_power, 1.0e16 + 1.0);
        assert_eq!(diagnostics.negative_real_power, -1.0e16 - 1.0);
        assert_eq!(diagnostics.apparent_power, 2.0e16 + 2.0);
        assert_eq!(
            diagnostics.plane_wave_reference_power,
            2.0 * medium.density * medium.sound_speed
        );

        let (negative, negative_interval, _) =
            evaluate_radiated_power(&[C64::from_re(-2.0)], &[C64::ONE], &[1.0], medium);
        assert_eq!(negative.to_bits(), (-1.0_f64).to_bits());
        assert!(negative_interval.1 < 0.0);
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
    fn triangle_duffy_self_terms_obey_helmholtz_scale_covariance() {
        let triangle = [[-0.8, -0.4, 0.0], [1.1, -0.3, 0.0], [0.2, 1.3, 0.0]];
        let centroid = [
            (triangle[0][0] + triangle[1][0] + triangle[2][0]) / 3.0,
            (triangle[0][1] + triangle[1][1] + triangle[2][1]) / 3.0,
            0.0,
        ];
        let k = 1.7;
        let (single, regularized_hypersingular) = triangle_self_terms(k, centroid, triangle);
        assert!(single.re > 0.0 && single.im > 0.0);
        assert!(regularized_hypersingular.re > 0.0);

        let scale = 3.25;
        let scaled_triangle = triangle.map(|point| point.map(|value| scale * value));
        let scaled_centroid = centroid.map(|value| scale * value);
        let (scaled_single, scaled_regularized_hypersingular) =
            triangle_self_terms(k / scale, scaled_centroid, scaled_triangle);
        assert!((scaled_single - single.scale(scale)).abs() < 2.0e-13 * scaled_single.abs());
        assert!(
            (scaled_regularized_hypersingular - regularized_hypersingular.scale(1.0 / scale)).abs()
                < 2.0e-13 * scaled_regularized_hypersingular.abs()
        );
    }

    #[test]
    fn thin_opposing_panels_retain_the_exact_static_solid_angle_and_potential() {
        // An independent rectangular-panel integral checks both triangle
        // orientation and the near-face limit. A plain 8x8 area rule misses
        // most of the normal-derivative peak at the cymbal-sized gaps below.
        let a = 0.010_f64;
        let b = 0.007_f64;
        let triangles = [
            [[-a, -b, 0.0], [a, -b, 0.0], [a, b, 0.0]],
            [[-a, -b, 0.0], [a, b, 0.0], [-a, b, 0.0]],
        ];
        for height in [0.0016, 0.0005, 0.000001, -0.0005] {
            let h = f64::abs(height);
            let radius = (a * a + b * b + h * h).sqrt();
            let angle = det::atan2(a * b, h * radius);
            let single_ref = (a * det::ln((b + radius) / (a * a + h * h).sqrt())
                + b * det::ln((a + radius) / (b * b + h * h).sqrt())
                - h * angle)
                / core::f64::consts::PI;
            let double_ref = angle.copysign(height) / core::f64::consts::PI;
            let (mut single, mut double) = (C64::ZERO, C64::ZERO);
            for triangle in triangles {
                let (s, d) =
                    triangle_weak_influence(0.0, [0.0, 0.0, height], triangle, [0.0, 0.0, 1.0]);
                single = single + s;
                double = double + d;
            }
            assert!(
                (single.re - single_ref).abs() < 2.0e-12 * single_ref,
                "gap {height}: single={single:?}, exact={single_ref}"
            );
            assert!(
                (double.re - double_ref).abs() < 2.0e-12 * double_ref.abs(),
                "gap {height}: double={double:?}, exact={double_ref}"
            );
            assert_eq!(single.im, 0.0);
            assert_eq!(double.im, 0.0);
        }
    }

    #[test]
    fn closed_thin_surface_preserves_the_static_double_layer_identity() {
        // The exterior trace of the static double layer of constant density
        // vanishes: (D_0 + I/2) 1 = 0, so the CBIE matrix has row sum -1.
        // Mixing exact near panels with centroid far panels broke this identity.
        let surface = thin_box_surface();
        let n = surface.centroids().len();
        let triangles = surface.triangles().expect("triangle geometry");
        assert!((0..n).any(|i| (0..n).any(|j| {
            !triangle_is_near(surface.centroids()[i], surface.centroids()[j], triangles[j])
        })));
        let (matrix, _) = assemble_dense(&surface, 0.0, C64::ZERO);
        for row in matrix.chunks_exact(n) {
            let sum = row.iter().fold(C64::ZERO, |sum, &value| sum + value);
            assert!((sum.re + 1.0).abs() < 2.0e-12, "static row sum: {sum:?}");
            assert_eq!(sum.im, 0.0);
        }
    }

    #[test]
    fn g0_burton_miller_static_rows_and_weighted_flux_are_compatible() {
        // At k=0 all kernels are real. alpha=i exposes N_0 in Im(A),
        // D_0-I/2 in Re(A), and D'_0+I/2 in Im(B), independently of S_0.
        // These are the actual assembled operators, not a repaired solution.
        let sphere = SpherePanels::icosphere(1.0, 1).unwrap();
        let triangles = sphere
            .triangles()
            .unwrap()
            .iter()
            .map(|triangle| triangle.map(|p| [1.7 * p[0], 0.8 * p[1], 0.6 * p[2]]))
            .collect();
        let ellipsoid = SpherePanels::from_triangles(triangles).unwrap();
        for surface in [&sphere, &ellipsoid] {
            let n = surface.centroids().len();
            let (a, b) = assemble_dense(surface, 0.0, C64::new(0.0, 1.0));
            let mut old_defect = 0.0_f64;
            for i in 0..n {
                let row = &a[i * n..(i + 1) * n];
                let d_sum: f64 = row.iter().map(|z| z.re).sum();
                let n_sum: f64 = row.iter().map(|z| z.im).sum();
                let n_norm: f64 = row.iter().map(|z| z.im.abs()).sum();
                assert!((d_sum + 1.0).abs() < 2e-13);
                assert!(n_sum.abs() < 2e-13 * n_norm);
                // The omitted old diagonal really violates the static law.
                old_defect = old_defect.max((row[i].re + 0.5).abs());
                let mut flux = 0.0;
                let mut flux_norm = 0.0;
                for j in 0..n {
                    let term = surface.areas()[j] * b[j * n + i].im;
                    flux += term;
                    flux_norm += term.abs();
                    let dji = a[i * n + j].re + if i == j { 1.0 } else { 0.0 };
                    let adjoint = surface.areas()[i] * dji;
                    assert!((term - adjoint).abs() < 2e-13 * surface.areas()[i]);
                }
                assert!(flux.abs() < 2e-13 * flux_norm);
            }
            assert!(old_defect > 1e-3, "static-flux falsifier must be live");
        }
    }

    #[test]
    fn thin_plate_bending_fields_radiate_positive_power_across_the_band() {
        // G1 physical falsifier: the two skins have opposite normal velocities
        // for the same transverse bending motion. At low frequency their small
        // radiated field is especially sensitive to lost far-panel moments.
        // The previous centroid far rule produced six negative powers at 40 Hz
        // and one at 1640 Hz in this 16-field family, despite an admitted mesh.
        let surface = thin_box_surface();
        let medium = Medium {
            density: 1.2,
            sound_speed: 343.0,
        };
        let mut fields = Vec::new();
        for mx in 1..=4 {
            for my in 1..=4 {
                fields.push(
                    surface
                        .centroids()
                        .iter()
                        .zip(surface.normals())
                        .map(|(point, normal)| {
                            let x_phase =
                                mx as f64 * core::f64::consts::PI * (point[0] / 0.30 + 0.5);
                            let y_phase =
                                my as f64 * core::f64::consts::PI * (point[1] / 0.24 + 0.5);
                            C64::from_re(normal[2] * det::sin(x_phase) * det::sin(y_phase))
                        })
                        .collect::<Vec<_>>(),
                );
            }
        }
        let field_refs: Vec<&[C64]> = fields.iter().map(Vec::as_slice).collect();
        for frequency in [40.0, 1640.0] {
            let k = core::f64::consts::TAU * frequency / medium.sound_speed;
            let solutions =
                solve_radiation_batch(&surface, k, medium, &field_refs, Formulation::PlainCbie)
                    .expect("mesh-resolved thin-plate radiation");
            for (index, solution) in solutions.iter().enumerate() {
                assert!(
                    solution.radiated_power_roundoff_interval.0 > 0.0,
                    "{frequency} Hz, bending field ({}, {}): {:?}",
                    index / 4 + 1,
                    index % 4 + 1,
                    solution.radiated_power_roundoff_interval
                );
            }
        }
    }

    #[test]
    fn far_panel_phase_matches_the_analytic_rectangular_aperture() {
        // Independent Fourier integral over a rectangle: area times two sinc
        // factors, with the physical translation phase. Sampling only the two
        // triangle centroids loses more than 2% at the admitted-band scale.
        let a = 0.025;
        let b = 0.020;
        let z = 0.011;
        let triangles = [
            [[-a, -b, z], [a, -b, z], [a, b, z]],
            [[-a, -b, z], [a, b, z], [-a, b, z]],
        ];
        let raw_direction = [1.5, 0.7, 1.5];
        let direction = raw_direction.map(|v| v / norm(raw_direction));
        let sinc = |x: f64| if x == 0.0 { 1.0 } else { det::sin(x) / x };
        for k in [0.0, 1.0, 30.0] {
            let expected = expik(-k * direction[2] * z)
                .scale(4.0 * a * b * sinc(k * direction[0] * a) * sinc(k * direction[1] * b));
            let actual = triangles.iter().fold(C64::ZERO, |sum, &triangle| {
                sum + triangle_far_phase(k, direction, triangle)
            });
            assert!(
                (actual - expected).abs() < 1.0e-7 * expected.abs(),
                "k={k}: integrated={actual:?}, rectangle={expected:?}"
            );
            if k == 30.0 {
                let centroid = triangles.iter().fold(C64::ZERO, |sum, triangle| {
                    let point = core::array::from_fn(|axis| {
                        (triangle[0][axis] + triangle[1][axis] + triangle[2][axis]) / 3.0
                    });
                    sum + expik(-k * dot(direction, point)).scale(2.0 * a * b)
                });
                assert!((centroid - expected).abs() > 0.02 * expected.abs());
            }
        }
    }

    #[test]
    fn default_finite_microphone_preserves_a_thin_plate_bending_transfer() {
        // The default percussion microphone compared with the independent,
        // adaptive Green-row owner. Its source solution and position are held
        // fixed: this checks observation quadrature, not another boundary solve.
        let surface = thin_box_surface();
        let medium = Medium {
            density: 1.2,
            sound_speed: 343.0,
        };
        let k = core::f64::consts::TAU * 1640.0 / medium.sound_speed;
        let velocity: Vec<_> = surface
            .centroids()
            .iter()
            .zip(surface.normals())
            .map(|(p, normal)| {
                C64::from_re(
                    normal[2]
                        * det::sin(3.0 * core::f64::consts::PI * (p[0] / 0.30 + 0.5))
                        * det::sin(2.0 * core::f64::consts::PI * (p[1] / 0.24 + 0.5)),
                )
            })
            .collect();
        let solution = solve_radiation(&surface, k, medium, &velocity, Formulation::PlainCbie)
            .expect("thin-plate boundary solution");
        let points = [[0.08, 0.05, 0.35]];
        let geometry = crate::near_field::Geometry::new(&surface, &points, 1.0e-4)
            .expect("default microphone is outside the closed plate");
        let reference = geometry
            .prepare(
                k,
                medium,
                crate::near_field::Options {
                    relative_tolerance: 1.0e-9,
                    ..crate::near_field::Options::default()
                },
            )
            .expect("independent adaptive observer")
            .evaluate(&solution)
            .expect("same boundary trace")
            .pressure[0];
        let actual = exterior_pressure_at_points_integrated(&surface, &solution, medium, &points)
            .expect("default finite-point observer")[0];
        assert!(
            (actual - reference).abs() < 1.0e-6 * reference.abs(),
            "finite microphone: {actual:?}, adaptive reference: {reference:?}"
        );
        // Explicit centroid consumers retain their original numerical arm.
        // The discrepancy is material here, so accidentally routing both APIs
        // to the same rule cannot pass this comparison.
        let centroid = exterior_pressure_at_points(&surface, &solution, medium, &points)
            .expect("legacy centroid observer")[0];
        assert!((centroid - reference).abs() > 0.02 * reference.abs());
        let directions = [[1.5, 0.7, 1.5]];
        let integrated_far = far_field_integrated(&surface, &solution, medium, &directions)[0];
        let centroid_far = far_field(&surface, &solution, medium, &directions)[0];
        assert!(integrated_far.re.is_finite() && integrated_far.im.is_finite());
        assert!((integrated_far - centroid_far).abs() > 0.02 * integrated_far.abs());
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
        // PlainCbie is the non-resonant arm (fictitious frequencies start
        // at ka = pi). BurtonMiller must also remain accurate and passive
        // below ka=0.5, where inconsistent static D/D' self coefficients
        // previously produced a negative resistance amplified by |alpha| = 1/k.
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
        for ka in [0.05, 0.2, 0.5, 1.0, 2.0] {
            run(&coarse, ka, Formulation::BurtonMiller, 0.08);
        }
        run(&fine, 5.0, Formulation::BurtonMiller, 0.08);
        println!(
            "{{\"suite\":\"fs-bem-helmholtz\",\"case\":\"pulsating-sphere\",\"rows\":[{rows}],\"verdict\":\"pass\"}}"
        );
    }

    #[test]
    fn finite_exterior_pressure_matches_pulsating_sphere_and_far_field() {
        let surface = SpherePanels::icosphere(1.0, 2).expect("icosphere");
        let medium = Medium::air();
        let k = 1.0;
        let solution = solve_radiation(
            &surface,
            k,
            medium,
            &uniform_velocity(surface.centroids().len()),
            Formulation::PlainCbie,
        )
        .expect("solve");
        let radii = [2.0, 8.0, 64.0];
        let points = radii.map(|r| [0.0, 0.0, r]);
        let finite = exterior_pressure_at_points(&surface, &solution, medium, &points)
            .expect("finite exterior pressure");
        let surface_pressure = pulsating_sphere_impedance(k, 1.0);
        for (&r, &actual) in radii.iter().zip(&finite) {
            let expected = surface_pressure.scale(1.0 / r) * expik(k * (r - 1.0));
            let relative_error = (actual - expected).abs() / expected.abs();
            assert!(
                relative_error < 0.08,
                "r={r}: finite pressure {actual:?} vs analytic {expected:?}, rel {relative_error:.4}"
            );
        }
        let far = far_field(&surface, &solution, medium, &[[0.0, 0.0, 1.0]])[0];
        let scaled = |index: usize| finite[index].scale(radii[index]) * expik(-k * radii[index]);
        let near_error = (scaled(1) - far).abs();
        let distant_error = (scaled(2) - far).abs();
        assert!(distant_error < near_error);
        assert!(distant_error < 0.02 * far.abs());
        assert!(matches!(
            exterior_pressure_at_points(&surface, &solution, medium, &[[0.0, 0.0, 0.0]]),
            Err(HelmholtzError::BadParameter { .. })
        ));
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
        // MUTATION: alpha = -i/k (the wrong sign under e^{-i omega t}
        // for this module's HBIE row sign). Any nonzero Im(alpha) keeps
        // the system uniquely solvable, so the mutation still beats plain
        // CBIE here; what it loses is accuracy: its ka = pi error must
        // exceed twice the production arm's.
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
    fn burton_miller_sign_tracks_cbie_accuracy_off_resonance() {
        // G2 + G3: away from interior resonances plain CBIE is exact in the
        // continuum, so a consistent coupling must not add discretization
        // error. Measured on the pulsating sphere (icosphere-2, 320 panels):
        // production/CBIE 0.58, 0.60, 0.78 and mutation/production 2.8, 3.2,
        // 2.7 at ka = 2, 2.5, 4. Under refinement to icosphere-3 at ka = 2.5
        // the production error falls 6.9x (CBIE 3.8x, mutation 2.7x).
        let err = |surface: &SpherePanels, ka: f64, formulation: Formulation| -> f64 {
            let n = surface.centroids().len();
            let sol = solve_radiation(
                surface,
                ka,
                Medium::air(),
                &uniform_velocity(n),
                formulation,
            )
            .expect("solve");
            let z_ref = pulsating_sphere_impedance(ka, 1.0);
            (mean_impedance(&sol, surface) - z_ref).abs() / z_ref.abs()
        };
        let coarse = SpherePanels::icosphere(1.0, 2).expect("icosphere-2");
        for ka in [2.0, 2.5, 4.0] {
            let bm = err(&coarse, ka, Formulation::BurtonMiller);
            let cbie = err(&coarse, ka, Formulation::PlainCbie);
            let mutated = err(&coarse, ka, Formulation::BurtonMillerWrongAlphaSign);
            assert!(
                bm < cbie,
                "ka={ka}: production {bm:.3e} must beat plain CBIE {cbie:.3e}"
            );
            assert!(
                mutated > 2.0 * bm,
                "ka={ka}: mutation {mutated:.3e} must exceed twice production {bm:.3e}"
            );
            println!(
                "{{\"suite\":\"fs-bem-helmholtz\",\"case\":\"alpha-sign-off-resonance\",\"ka\":{ka},\"bm_rel\":{bm:.4e},\"cbie_rel\":{cbie:.4e},\"mutation_rel\":{mutated:.4e}}}"
            );
        }
        let fine = SpherePanels::icosphere(1.0, 3).expect("icosphere-3");
        let ratio = err(&coarse, 2.5, Formulation::BurtonMiller)
            / err(&fine, 2.5, Formulation::BurtonMiller);
        assert!(
            ratio > 3.0,
            "production must converge at least at the CBIE rate: {ratio:.2}x"
        );
        println!(
            "{{\"suite\":\"fs-bem-helmholtz\",\"case\":\"alpha-sign-refinement\",\"ka\":2.5,\"error_ratio_ico2_over_ico3\":{ratio:.3}}}"
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
            let pressure_scale = sol
                .pressure
                .iter()
                .map(|p| p.abs())
                .fold(0.0_f64, f64::max);
            for (row, &pressure) in z.chunks_exact(n).zip(&sol.pressure) {
                let coherent = row
                    .iter()
                    .zip(&v)
                    .fold(C64::ZERO, |sum, (&zij, &vj)| sum + zij * vj);
                assert!(
                    (coherent - pressure).abs() < 1e-11 * pressure_scale,
                    "matrix columns lost coherent source phase in probe {idx}"
                );
            }
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
            let z = baffled_piston_face_impedance(1.0, ka, Medium::air()).expect("piston");
            let z_series =
                C64::new(0.5 * ka * ka, -8.0 * ka / (3.0 * core::f64::consts::PI)).scale(rho_c());
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

    /// G3: compare the one historical polar rule and an independent smooth
    /// disk-overlap integral with the closed form. Equivalent-disc replacement
    /// of the singular near field is not a smooth midpoint rule: the leading
    /// O(h) contribution can cancel other errors on coarse meshes. In particular
    /// 8 -> 16 rings reverses before 16 -> 32 -> 64 improves. This is not an
    /// observed error floor, nor evidence for three independent implementations.
    #[test]
    fn historical_ring_and_overlap_integral_match_piston() {
        let start = std::time::Instant::now();
        let air = Medium::air();
        let closed = baffled_piston_face_impedance(1.0, 1.0, air).unwrap();
        let relative_error = |candidate: C64| (candidate - closed).abs() / closed.abs();
        let mut errors = Vec::new();
        for rings in [8, 16, 32, 64] {
            let value = historical_ring_piston(1.0, 1.0, air, rings).unwrap();
            let err = relative_error(value);
            errors.push(err);
            println!(
                "{{\"case\":\"historical-ring-piston\",\"inputs\":\"unit disk, ka=1; one historical algorithm\",\"axis\":\"rings\",\"rung\":{rings},\"observed\":{err:e},\"oracle\":\"Bessel-Struve\",\"tol\":\"record coarse reversal; final rung < 0.002\",\"verdict\":\"measured\",\"elapsed_ms\":{}}}",
                start.elapsed().as_millis()
            );
            if rings == 64 {
                assert!(err < 0.002, "final historical error: {errors:?}");
                // The same acceptance band must reject a planted 1% load error.
                assert!(relative_error(value.scale(1.01)) > 0.002);
            }
        }
        assert!(errors[1] > errors[0], "coarse cancellation: {errors:?}");
        assert!(errors[3] < errors[2] && errors[2] < errors[1]);
        let ladder = fs_mms::RefinementLadder::new(
            vec![1.0 / 16.0, 1.0 / 32.0, 1.0 / 64.0],
            errors[1..].to_vec(),
        )
        .unwrap();
        let fit = fs_mms::fit_order(&ladder);
        println!(
            "{{\"case\":\"historical-ring-order\",\"inputs\":\"unit disk, ka=1; 16/32/64 rings\",\"axis\":\"rings\",\"observed\":{},\"rms_residual\":{},\"tol\":\"record-only\",\"verdict\":\"measured\",\"elapsed_ms\":{}}}",
            fit.observed,
            fit.rms_residual,
            start.elapsed().as_millis()
        );

        // Independently reduce the double Rayleigh integral by convolution:
        // A(r) = area of overlap of two unit disks separated by r,
        // Z_face = -i k rho c / pi * INT_0^2 A(r) exp(i k r) dr.
        // Set r=2 sin(theta): A=2(pi/2-theta-sin(theta)cos(theta)).
        // This removes the square-root endpoint and the 1/r singularity;
        // composite Simpson now has a smooth, fourth-order integrand. It calls
        // neither J1/H1 nor the historical equivalent-disc cell kernel.
        let overlap_integral = |intervals: usize| {
            let step = core::f64::consts::FRAC_PI_2 / intervals as f64;
            let mut sum = C64::ZERO;
            for index in 0..=intervals {
                let theta = index as f64 * step;
                let (s, c) = (det::sin(theta), det::cos(theta));
                let overlap = 2.0 * (core::f64::consts::FRAC_PI_2 - theta - s * c);
                let weight = if index == 0 || index == intervals {
                    1.0
                } else if index % 2 == 1 {
                    4.0
                } else {
                    2.0
                };
                sum = sum + expik(2.0 * s).scale(overlap * 2.0 * c * weight);
            }
            sum * C64::new(0.0, -rho_c() * step / (3.0 * core::f64::consts::PI))
        };
        let mut overlap_errors = Vec::new();
        for intervals in [64, 128, 256, 512] {
            let value = overlap_integral(intervals);
            let err = relative_error(value);
            if let Some(previous) = overlap_errors.last() {
                assert!(err < *previous);
            }
            overlap_errors.push(err);
            println!(
                "{{\"case\":\"disk-overlap-piston\",\"inputs\":\"unit disk, ka=1; independent Rayleigh convolution\",\"axis\":\"Simpson intervals\",\"rung\":{intervals},\"observed\":{err:e},\"oracle\":\"Bessel-Struve\",\"tol\":\"decreasing; final rung < 1e-10\",\"verdict\":\"measured\",\"elapsed_ms\":{}}}",
                start.elapsed().as_millis()
            );
            if intervals == 512 {
                assert!(err < 1e-10);
                assert!(relative_error(value.scale(1.01)) > 1e-10);
                assert!((value - closed.scale(1.01)).abs() / value.abs() > 1e-10);
            }
        }
        let ladder = fs_mms::RefinementLadder::new(
            vec![1.0 / 64.0, 1.0 / 128.0, 1.0 / 256.0, 1.0 / 512.0],
            overlap_errors,
        )
        .unwrap();
        let fit = fs_mms::fit_order(&ladder);
        assert!((3.8..4.2).contains(&fit.observed), "{fit:?}");
        println!(
            "{{\"case\":\"disk-overlap-order\",\"inputs\":\"unit disk, ka=1\",\"axis\":\"Simpson intervals\",\"observed\":{},\"rms_residual\":{},\"tol\":\"3.8..4.2\",\"verdict\":\"pass\",\"elapsed_ms\":{}}}",
            fit.observed,
            fit.rms_residual,
            start.elapsed().as_millis()
        );
    }

    #[test]
    fn far_field_table_matches_pulsating_sphere() {
        let surface = SpherePanels::icosphere(1.0, 1).expect("icosphere");
        let n = surface.centroids().len();
        let medium = Medium::air();
        let ka = 0.4;
        let omega = ka * medium.sound_speed;
        let table = tabulate_far_field(
            &surface,
            &[omega],
            &[[0.0, 0.0, 1.0]],
            medium,
            &uniform_velocity(n),
            Formulation::PlainCbie,
        )
        .expect("table");
        let f = table.amplitude(0, 0).expect("F");
        let z = pulsating_sphere_impedance(ka, 1.0);
        let f_ref = z.scale(1.0) * C64::new(det::cos(-ka), det::sin(-ka));
        let rel = (f - f_ref).abs() / f_ref.abs().max(1.0e-18);
        assert!(
            rel < 0.08,
            "tabulated F {f:?} vs analytic {f_ref:?} (rel {rel})"
        );
        assert!(
            tabulate_far_field(
                &surface,
                &[],
                &[[1.0, 0.0, 0.0]],
                medium,
                &uniform_velocity(n),
                Formulation::PlainCbie
            )
            .is_err()
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
        assert_eq!(
            a.condition_lower_bound.to_bits(),
            b.condition_lower_bound.to_bits()
        );
    }

    #[test]
    fn g0_batched_radiation_matches_independent_solves_bitwise() {
        let surface = SpherePanels::icosphere(1.0, 1).expect("icosphere");
        let n = surface.centroids().len();
        let uniform = uniform_velocity(n);
        let dipole: Vec<C64> = surface
            .normals()
            .iter()
            .map(|normal| C64::from_re(normal[0]))
            .collect();
        let mixed: Vec<C64> = surface
            .normals()
            .iter()
            .map(|normal| C64::new(normal[1] + 0.25 * normal[2], -0.3 * normal[0]))
            .collect();
        let fields = [uniform.as_slice(), dipole.as_slice(), mixed.as_slice()];
        let batched = solve_radiation_batch(
            &surface,
            1.3,
            Medium::air(),
            &fields,
            Formulation::BurtonMiller,
        )
        .expect("batched solve");
        let replay = solve_radiation_batch(
            &surface,
            1.3,
            Medium::air(),
            &fields,
            Formulation::BurtonMiller,
        )
        .expect("deterministic batched replay");
        assert_eq!(batched.len(), fields.len());
        for (index, &velocity) in fields.iter().enumerate() {
            let independent = solve_radiation(
                &surface,
                1.3,
                Medium::air(),
                velocity,
                Formulation::BurtonMiller,
            )
            .expect("independent solve");
            assert_solution_bitwise_equal(&batched[index], &independent);
            assert_solution_bitwise_equal(&replay[index], &batched[index]);
        }
        assert!(
            batched
                .windows(2)
                .all(|pair| pair[0].pressure != pair[1].pressure)
        );

        let empty: [&[C64]; 0] = [];
        assert!(matches!(
            solve_radiation_batch(
                &surface,
                1.3,
                Medium::air(),
                &empty,
                Formulation::BurtonMiller
            ),
            Err(HelmholtzError::BadParameter { .. })
        ));
        let short = vec![C64::ZERO; n - 1];
        assert!(matches!(
            solve_radiation_batch(
                &surface,
                1.3,
                Medium::air(),
                &[uniform.as_slice(), short.as_slice()],
                Formulation::BurtonMiller
            ),
            Err(HelmholtzError::ShapeMismatch { .. })
        ));
        let oversized = vec![uniform.as_slice(); MAX_RADIATION_FIELDS_PER_BATCH + 1];
        assert!(matches!(
            solve_radiation_batch(
                &surface,
                1.3,
                Medium::air(),
                &oversized,
                Formulation::BurtonMiller
            ),
            Err(HelmholtzError::BadParameter { .. })
        ));
    }

    #[test]
    fn prepared_cbie_geometry_reuses_static_terms_without_changing_frequency_solutions() {
        // A flattened, triangulated source exercises closely spaced skins,
        // near/far rows, and non-axis-aligned self normals on a small mesh.
        let sphere = SpherePanels::icosphere(0.12, 1).unwrap();
        let triangles = sphere
            .triangles()
            .unwrap()
            .iter()
            .map(|triangle| triangle.map(|[x, y, z]| [x, y, z / 240.0]))
            .collect();
        let surface = SpherePanels::from_triangles(triangles).unwrap();
        let n = surface.centroids().len();
        let prepared = PreparedCbieGeometry::new(&surface, 10.0, 16 * n * n).unwrap();
        assert_eq!(prepared.cache_bytes(), 16 * n * n);
        let uniform = vec![C64::ONE; n];
        let bending: Vec<_> = surface
            .centroids()
            .iter()
            .zip(surface.normals())
            .map(|(p, normal)| C64::new(normal[2] * (1.0 - p[0] * p[0] / 0.0144), 0.0))
            .collect();
        let complex: Vec<_> = surface
            .normals()
            .iter()
            .enumerate()
            .map(|(i, normal)| C64::new(i as f64 * 0.01, -normal[0]))
            .collect();
        let fields = [uniform.as_slice(), bending.as_slice(), complex.as_slice()];
        for (k, medium) in [
            (
                0.7,
                Medium {
                    density: 1.2,
                    sound_speed: 343.0,
                },
            ),
            (
                10.0,
                Medium {
                    density: 1.4,
                    sound_speed: 320.0,
                },
            ),
        ] {
            let actual = prepared.solve_batch(k, medium, &fields).unwrap();
            let expected =
                solve_radiation_batch(&surface, k, medium, &fields, Formulation::PlainCbie)
                    .unwrap();
            for (actual, expected) in actual.iter().zip(&expected) {
                assert_solution_bitwise_equal(actual, expected);
            }
        }
        assert!(matches!(
            prepared.solve_batch(1.0, Medium::air(), &[&uniform[..n - 1]]),
            Err(HelmholtzError::ShapeMismatch { .. })
        ));
    }

    #[test]
    fn prepared_cbie_geometry_admits_resolution_memory_and_cancellation_before_publication() {
        let surface = SpherePanels::icosphere(0.1, 0).unwrap();
        let bytes = 16 * surface.centroids().len().pow(2);
        // Wavelength refusal wins before the zero memory budget or any rows.
        assert!(matches!(
            PreparedCbieGeometry::new(&surface, 1e6, 0),
            Err(HelmholtzError::TooCoarse { .. })
        ));
        assert!(
            matches!(PreparedCbieGeometry::new(&surface, 1.0, bytes - 1),
            Err(HelmholtzError::CacheBudget { required_bytes, maximum_bytes })
                if required_bytes == bytes && maximum_bytes == bytes - 1)
        );
        let mut polls = 0;
        let cancelled = PreparedCbieGeometry::new_with_cancel(&surface, 1.0, bytes, || {
            polls += 1;
            polls == 3
        });
        assert!(matches!(cancelled, Err(HelmholtzError::Cancelled)));
        assert_eq!(polls, 3);
        let prepared = PreparedCbieGeometry::new(&surface, 1.0, bytes).unwrap();
        assert!(matches!(
            prepared.solve_batch(
                1e6,
                Medium::air(),
                &[&vec![C64::ONE; surface.centroids().len()]]
            ),
            Err(HelmholtzError::TooCoarse { .. })
        ));
    }

    #[test]
    fn gauss_legendre_and_normalized_legendre_pin() {
        // Quadrature identities: SUM w = 2, INT x^2 = 2/3, INT x^10 =
        // 2/11 (all exact for the node counts used); then the
        // orthonormality of Pbar under that same rule:
        // 2 pi INT Pbar_l^m Pbar_l'^m dx = delta_ll'.
        let (nodes, weights) = gauss_legendre(9);
        let wsum: f64 = weights.iter().sum();
        assert!((wsum - 2.0).abs() < 1e-14, "sum of GL weights {wsum}");
        let x2: f64 = nodes
            .iter()
            .zip(weights.iter())
            .map(|(&x, &w)| w * x * x)
            .sum();
        assert!((x2 - 2.0 / 3.0).abs() < 1e-14);
        let x10: f64 = nodes
            .iter()
            .zip(weights.iter())
            .map(|(&x, &w)| w * det::powi(x, 10))
            .sum();
        assert!((x10 - 2.0 / 11.0).abs() < 1e-13);
        let l_max = 8usize;
        let (nodes, weights) = gauss_legendre(l_max + 1);
        let idx = |l: usize, m: usize| l * (l + 1) / 2 + m;
        for m in 0..=3usize {
            for la in m..=l_max {
                for lb in m..=l_max {
                    let mut acc = 0.0;
                    for (i, &x) in nodes.iter().enumerate() {
                        let p = norm_assoc_legendre(l_max, x);
                        acc += weights[i] * p[idx(la, m)] * p[idx(lb, m)];
                    }
                    acc *= 2.0 * core::f64::consts::PI;
                    let expected = if la == lb { 1.0 } else { 0.0 };
                    assert!(
                        (acc - expected).abs() < 1e-12,
                        "orthonormality (l={la},l'={lb},m={m}): {acc}"
                    );
                }
            }
        }
        println!(
            "{{\"suite\":\"fs-bem-helmholtz\",\"case\":\"sh-basis-pins\",\"verdict\":\"pass\"}}"
        );
    }

    #[test]
    fn generic_far_field_projection_preserves_an_analytic_monopole() {
        let amplitude = C64::new(2.5, -0.75);
        let table = directivity_sh_table_from_evaluator(3.0, 6, |directions| {
            Ok(vec![amplitude; directions.len()])
        })
        .expect("analytic far field");
        assert!(table.captured_fraction > 1.0 - 1.0e-12);
        for direction in [[1.0, 0.0, 0.0], [0.3, -0.4, 0.5], [0.0, 0.0, -1.0]] {
            let reconstructed = table.evaluate(direction);
            assert!((reconstructed - amplitude).abs() < 1.0e-11);
        }
        assert!(matches!(
            directivity_sh_table_from_evaluator(3.0, 2, |_| Ok(Vec::new())),
            Err(HelmholtzError::BadParameter { .. })
        ));
    }

    #[test]
    fn sh_table_concentrates_and_reconstructs() {
        // Monopole power must land in a_00, dipole power in (l=1, m=0),
        // the truncation must capture essentially all quadrature power,
        // and evaluate() must reproduce the sampled far field it was
        // built from. Repeat construction is bitwise identical.
        let surface = SpherePanels::icosphere(1.0, 2).expect("icosphere");
        let n = surface.centroids().len();
        let ka = 1.0;
        let l_max = 8;
        let mono = solve_radiation(
            &surface,
            ka,
            Medium::air(),
            &uniform_velocity(n),
            Formulation::BurtonMiller,
        )
        .expect("solve");
        let t_mono =
            directivity_sh_table(&surface, &mono, Medium::air(), l_max).expect("monopole table");
        let spec = t_mono.power_by_degree();
        let total: f64 = spec.iter().sum();
        assert!(
            spec[0] / total > 0.999,
            "monopole must concentrate at l=0: {:.6}",
            spec[0] / total
        );
        assert!(
            t_mono.captured_fraction > 0.999,
            "captured fraction {:.6}",
            t_mono.captured_fraction
        );
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
        let t_dip =
            directivity_sh_table(&surface, &dip, Medium::air(), l_max).expect("dipole table");
        let spec = t_dip.power_by_degree();
        let total: f64 = spec.iter().sum();
        let a10 = t_dip.coefficient(1, 0).norm_sq();
        assert!(
            a10 / total > 0.99,
            "dipole must concentrate at (l=1, m=0): {:.6}",
            a10 / total
        );
        // Reconstruction against directly sampled far field at
        // directions the projection grid never used.
        let mut dirs = Vec::new();
        for i in 0..17u32 {
            let theta = core::f64::consts::PI * (f64::from(i) + 0.37) / 17.0;
            let phi = 2.61 * f64::from(i);
            dirs.push([
                det::sin(theta) * det::cos(phi),
                det::sin(theta) * det::sin(phi),
                det::cos(theta),
            ]);
        }
        let direct = far_field(&surface, &dip, Medium::air(), &dirs);
        let mut num = 0.0f64;
        let mut den = 0.0f64;
        for (i, d) in dirs.iter().enumerate() {
            let rec = t_dip.evaluate(*d);
            num += (rec - direct[i]).norm_sq();
            den += direct[i].norm_sq();
        }
        let rel = (num / den).sqrt();
        assert!(rel < 1e-3, "SH reconstruction relative L2 error {rel:.2e}");
        assert!(
            t_dip.captured_fraction > 0.999,
            "dipole captured fraction {:.6}",
            t_dip.captured_fraction
        );
        let t_repeat =
            directivity_sh_table(&surface, &dip, Medium::air(), l_max).expect("repeat table");
        for (x, y) in t_dip.coefficients.iter().zip(t_repeat.coefficients.iter()) {
            assert!(x.re.to_bits() == y.re.to_bits() && x.im.to_bits() == y.im.to_bits());
        }
        println!(
            "{{\"suite\":\"fs-bem-helmholtz\",\"case\":\"sh-directivity-table\",\"mono_l0_fraction\":{:.6},\"dip_l1m0_fraction\":{:.6},\"captured_fraction\":{:.6},\"reconstruction_rel\":{rel:.2e},\"verdict\":\"pass\"}}",
            t_mono.power_by_degree()[0] / t_mono.power_by_degree().iter().sum::<f64>(),
            a10 / total,
            t_dip.captured_fraction
        );
    }

    #[test]
    fn off_axis_dipole_pins_nonaxisymmetric_harmonics() {
        // FALSIFIER for one-sided sign/phase bugs the axisymmetric
        // fixtures cannot see: an x-axis dipole has
        // F proportional to sin(theta) cos(phi)
        //   = sqrt(2 pi / 3) (Y_{1,-1} - Y_{1,1}),
        // so the table must put its power into (1, +/-1) with the exact
        // relation a_{1,-1} = -a_{1,1}. A missing Condon-Shortley phase
        // or a conjugation error in EITHER projection or evaluation
        // breaks that sign while leaving every m = 0 test green.
        let surface = SpherePanels::icosphere(1.0, 2).expect("icosphere");
        let v_x: Vec<C64> = surface
            .normals()
            .iter()
            .map(|nrm| C64::from_re(nrm[0]))
            .collect();
        let sol = solve_radiation(
            &surface,
            1.0,
            Medium::air(),
            &v_x,
            Formulation::BurtonMiller,
        )
        .expect("solve");
        let table = directivity_sh_table(&surface, &sol, Medium::air(), 8).expect("table");
        let spec = table.power_by_degree();
        let total: f64 = spec.iter().sum();
        assert!(
            spec[1] / total > 0.99,
            "x-dipole power must land at l = 1: {:.6}",
            spec[1] / total
        );
        let a_pos = table.coefficient(1, 1);
        let a_neg = table.coefficient(1, -1);
        let a_zero = table.coefficient(1, 0);
        assert!(
            a_pos.abs() > 100.0 * a_zero.abs(),
            "x-dipole must be m = +/-1, not m = 0: |a11| {:.3e} vs |a10| {:.3e}",
            a_pos.abs(),
            a_zero.abs()
        );
        let sign_defect = (a_neg + a_pos).abs() / a_pos.abs();
        assert!(
            sign_defect < 0.01,
            "a_(1,-1) = -a_(1,1) must hold: relative defect {sign_defect:.2e}"
        );
        // Evaluation must reproduce the m != 0 field at held-out
        // directions (this closes the loop through BOTH conventions).
        let mut dirs = Vec::new();
        for i in 0..17u32 {
            let theta = core::f64::consts::PI * (f64::from(i) + 0.61) / 17.0;
            let phi = 1.13 * f64::from(i) + 0.4;
            dirs.push([
                det::sin(theta) * det::cos(phi),
                det::sin(theta) * det::sin(phi),
                det::cos(theta),
            ]);
        }
        let direct = far_field(&surface, &sol, Medium::air(), &dirs);
        let mut num = 0.0f64;
        let mut den = 0.0f64;
        for (i, d) in dirs.iter().enumerate() {
            num += (table.evaluate(*d) - direct[i]).norm_sq();
            den += direct[i].norm_sq();
        }
        let rel = (num / den).sqrt();
        assert!(rel < 1e-3, "m != 0 reconstruction error {rel:.2e}");
        assert!(table.captured_fraction > 0.999);
        println!(
            "{{\"suite\":\"fs-bem-helmholtz\",\"case\":\"off-axis-dipole-sh\",\"l1_fraction\":{:.6},\"sign_defect\":{sign_defect:.2e},\"reconstruction_rel\":{rel:.2e},\"verdict\":\"pass\"}}",
            spec[1] / total
        );
    }

    #[test]
    fn radiation_efficiency_matches_sphere_oracles() {
        use core::fmt::Write as _;
        // sigma = W / (1/2 rho c INT |v|^2 dS) against the analytic
        // monopole (ka^2/(1+ka^2)) and dipole (ka^4/(4+ka^4)) sphere
        // efficiencies. PlainCbie is the accurate arm at these
        // non-resonant ka (the BM low-ka resistance artifact is the
        // recorded no-claim).
        let surface = SpherePanels::icosphere(1.0, 2).expect("icosphere");
        let n = surface.centroids().len();
        let mut rows = String::new();
        let mut first = true;
        for &ka in &[0.5, 1.0, 2.0] {
            let sol = solve_radiation(
                &surface,
                ka,
                Medium::air(),
                &uniform_velocity(n),
                Formulation::PlainCbie,
            )
            .expect("solve");
            let sigma = radiation_efficiency(&surface, &sol, Medium::air()).expect("efficiency");
            let oracle = ka * ka / (1.0 + ka * ka);
            let rel = (sigma - oracle).abs() / oracle;
            assert!(
                rel < 0.05,
                "monopole efficiency ka={ka}: {sigma:.4} vs {oracle:.4}"
            );
            write!(
                rows,
                "{}{{\"mode\":\"monopole\",\"ka\":{ka},\"sigma\":{sigma:.4},\"oracle\":{oracle:.4},\"rel\":{rel:.4}}}",
                if first { "" } else { "," }
            )
            .expect("write to String");
            first = false;
        }
        let ka = 1.0;
        let v_dip: Vec<C64> = surface
            .normals()
            .iter()
            .map(|nrm| C64::from_re(nrm[2]))
            .collect();
        let sol = solve_radiation(&surface, ka, Medium::air(), &v_dip, Formulation::PlainCbie)
            .expect("solve");
        let sigma = radiation_efficiency(&surface, &sol, Medium::air()).expect("efficiency");
        let ka4 = det::powi(ka, 4);
        let oracle = ka4 / (4.0 + ka4);
        let rel = (sigma - oracle).abs() / oracle;
        assert!(
            rel < 0.10,
            "dipole efficiency ka={ka}: {sigma:.4} vs {oracle:.4} (rel {rel:.4})"
        );
        write!(
            rows,
            ",{{\"mode\":\"dipole\",\"ka\":{ka},\"sigma\":{sigma:.4},\"oracle\":{oracle:.4},\"rel\":{rel:.4}}}"
        )
        .expect("write to String");
        // Refusal: a zero velocity field has no efficiency.
        let zero = RadiationSolution {
            pressure: vec![C64::ZERO; n],
            velocity: vec![C64::ZERO; n],
            k: ka,
            medium: Medium::air(),
            surface_fingerprint: surface_fingerprint(&surface),
            panels_per_wavelength: 10.0,
            radiated_power: 0.0,
            radiated_power_roundoff_interval: (0.0, 0.0),
            power_diagnostics: RadiationPowerDiagnostics {
                surface_power: C64::ZERO,
                positive_real_power: 0.0,
                negative_real_power: 0.0,
                apparent_power: 0.0,
                plane_wave_reference_power: 0.0,
            },
            condition_lower_bound: 1.0,
            dense_cap_utilization: 0.01,
        };
        assert!(matches!(
            radiation_efficiency(&surface, &zero, Medium::air()),
            Err(HelmholtzError::BadParameter { .. })
        ));
        println!(
            "{{\"suite\":\"fs-bem-helmholtz\",\"case\":\"radiation-efficiency\",\"rows\":[{rows}],\"verdict\":\"pass\"}}"
        );
    }

    #[test]
    fn condition_diagnostic_inflates_at_fictitious_frequency() {
        // The probe-based lower bound must SEE the interior resonance.
        // The flat-panel discretization SHIFTS the sphere's first
        // interior Dirichlet resonance away from the continuum ka = pi
        // (centroids sit at radius < 1), so the claim is scanned over
        // the resonance neighborhood: somewhere in ka in [3.0, 3.35]
        // plain CBIE's bound must spike visibly above Burton-Miller's
        // at the same ka and above its own off-resonance value.
        // Measured on this mesh (2026-08-08): the CBIE bound peaks near
        // ka ~ 3.20 at ~21x its own off-resonance value, while at
        // exactly pi the cross-arm ratio is only ~1.6x — asserting at a
        // single hardcoded ka would test the continuum resonance, not
        // the discrete operator actually solved.
        let surface = SpherePanels::icosphere(1.0, 2).expect("icosphere");
        let n = surface.centroids().len();
        let cond = |ka: f64, formulation: Formulation| -> f64 {
            solve_radiation(
                &surface,
                ka,
                Medium::air(),
                &uniform_velocity(n),
                formulation,
            )
            .expect("solve")
            .condition_lower_bound
        };
        let cbie_off = cond(1.0, Formulation::PlainCbie);
        // Stage 1: coarse scan locates the discrete resonance.
        let mut coarse_ka = 0.0f64;
        let mut coarse_peak = 0.0f64;
        for step in 0..=14u32 {
            let ka = 3.0 + 0.025 * f64::from(step);
            let cbie = cond(ka, Formulation::PlainCbie);
            assert!(cbie.is_finite());
            if cbie > coarse_peak {
                coarse_peak = cbie;
                coarse_ka = ka;
            }
        }
        // Stage 2: refine — the near-singular spike is much narrower
        // than the coarse step.
        let mut best_ka = coarse_ka;
        let mut peak_cbie = coarse_peak;
        for step in 0..=20u32 {
            let ka = coarse_ka - 0.025 + 0.0025 * f64::from(step);
            let cbie = cond(ka, Formulation::PlainCbie);
            if cbie > peak_cbie {
                peak_cbie = cbie;
                best_ka = ka;
            }
        }
        let bm_at_peak = cond(best_ka, Formulation::BurtonMiller);
        let bm_off = cond(1.0, Formulation::BurtonMiller);
        let ratio_vs_bm = peak_cbie / bm_at_peak;
        assert!(cbie_off >= 1.0 && bm_at_peak >= 1.0 && bm_off >= 1.0);
        // Measured on this mesh (2026-08-08): peak_cbie ~ 97 at
        // ka ~ 3.20, bm_at_peak ~ 33 — the 320-panel discretization
        // regularizes the resonance, so CBIE's bound rises ABOVE BM's
        // (~2.9x) but not by an order of magnitude. The discriminating
        // property is RELATIVE: CBIE inflates across the band while
        // BM stays flat — that contrast is what makes the diagnostic
        // usable as a fictitious-frequency warning.
        assert!(
            ratio_vs_bm > 1.5,
            "CBIE must exceed BM at the discrete resonance: \
             ratio {ratio_vs_bm:.1} at ka={best_ka:.4} (cbie {peak_cbie:.1}, bm {bm_at_peak:.1})"
        );
        assert!(
            peak_cbie > 3.0 * cbie_off,
            "CBIE resonance peak must dominate its own off-resonance value: \
             {peak_cbie:.1} vs {cbie_off:.1}"
        );
        let cbie_spike = peak_cbie / cbie_off;
        let bm_spike = bm_at_peak / bm_off;
        assert!(
            cbie_spike > 2.0 * bm_spike,
            "the diagnostic must inflate for the resonant arm, not both: \
             CBIE spike {cbie_spike:.1}x vs BM spike {bm_spike:.1}x"
        );
        // The cap-utilization diagnostic is the panel count over the cap.
        let sol = solve_radiation(
            &surface,
            1.0,
            Medium::air(),
            &uniform_velocity(n),
            Formulation::BurtonMiller,
        )
        .expect("solve");
        assert!((sol.dense_cap_utilization - n as f64 / MAX_DENSE_PANELS as f64).abs() < 1e-15);
        println!(
            "{{\"suite\":\"fs-bem-helmholtz\",\"case\":\"condition-diagnostic\",\"peak_cbie\":{peak_cbie:.1},\"ratio_vs_bm\":{ratio_vs_bm:.1},\"best_ka\":{best_ka:.4},\"cbie_off\":{cbie_off:.1},\"cbie_spike\":{cbie_spike:.1},\"bm_spike\":{bm_spike:.1},\"verdict\":\"pass\"}}"
        );
    }

    #[test]
    fn sh_degree_cap_refuses() {
        let surface = SpherePanels::icosphere(1.0, 1).expect("icosphere");
        let n = surface.centroids().len();
        let sol = solve_radiation(
            &surface,
            1.0,
            Medium::air(),
            &uniform_velocity(n),
            Formulation::BurtonMiller,
        )
        .expect("solve");
        let err = directivity_sh_table(&surface, &sol, Medium::air(), MAX_SH_DEGREE + 1)
            .expect_err("degree cap");
        assert!(err.to_string().contains("FS-BEM-HELM-BAD-PARAMETER"));
    }

    #[test]
    fn declared_dense_workers_leave_every_entry_and_solution_bitwise() {
        fn bits(values: &[C64]) -> Vec<(u64, u64)> {
            values
                .iter()
                .map(|v| (v.re.to_bits(), v.im.to_bits()))
                .collect()
        }
        let surface = SpherePanels::icosphere(1.0, 2).expect("icosphere");
        let n = surface.centroids().len();
        let velocity: Vec<C64> = surface
            .normals()
            .iter()
            .map(|v| C64::new(v[2], 0.25))
            .collect();
        for (k, formulation) in [
            (0.5, Formulation::PlainCbie),
            (2.0, Formulation::BurtonMiller),
        ] {
            let alpha = alpha_for(formulation, k);
            let (a, b) = assemble_dense(&surface, k, alpha);
            let solved = solve_radiation(&surface, k, Medium::air(), &velocity, formulation)
                .expect("solve")
                .pressure;
            for workers in [3, 7, n + 5] {
                let (wa, wb, ws) = with_dense_workers(workers, || {
                    let (wa, wb) = assemble_dense(&surface, k, alpha);
                    let ws = solve_radiation(&surface, k, Medium::air(), &velocity, formulation)
                        .expect("solve")
                        .pressure;
                    (wa, wb, ws)
                });
                assert_eq!(bits(&wa), bits(&a), "{formulation:?} workers {workers}");
                assert_eq!(bits(&wb), bits(&b), "{formulation:?} workers {workers}");
                assert_eq!(
                    bits(&ws),
                    bits(&solved),
                    "{formulation:?} workers {workers}"
                );
            }
            // A batch splits its fields over the workers; each field is still its own one-field solve, in order.
            let fields: Vec<Vec<C64>> = (0..5)
                .map(|f| {
                    let s = f as f64;
                    velocity
                        .iter()
                        .zip(surface.normals())
                        .map(|(&v, nrm)| v.scale(1.0 + s) + C64::new(s * nrm[0], -0.5 * s * nrm[1]))
                        .collect()
                })
                .collect();
            let refs: Vec<&[C64]> = fields.iter().map(Vec::as_slice).collect();
            let one: Vec<_> = refs
                .iter()
                .map(|&f| {
                    solve_radiation(&surface, k, Medium::air(), f, formulation).expect("solve")
                })
                .collect();
            let rule = PanelRule::Triangle { max_pairs: n * n };
            let ruled = solve_radiation_batch_with_rule(
                &surface,
                k,
                Medium::air(),
                &refs,
                formulation,
                rule,
            )
            .expect("ruled batch");
            for workers in [2, 3, 7, n + 5] {
                let (batch, ruled_batch) = with_dense_workers(workers, || {
                    (
                        solve_radiation_batch(&surface, k, Medium::air(), &refs, formulation)
                            .expect("batch"),
                        solve_radiation_batch_with_rule(
                            &surface,
                            k,
                            Medium::air(),
                            &refs,
                            formulation,
                            rule,
                        )
                        .expect("ruled batch"),
                    )
                });
                for (expected, actual) in
                    one.iter().zip(&batch).chain(ruled.iter().zip(&ruled_batch))
                {
                    assert_eq!(
                        bits(&actual.pressure),
                        bits(&expected.pressure),
                        "{formulation:?} workers {workers}"
                    );
                    assert_eq!(
                        actual.radiated_power.to_bits(),
                        expected.radiated_power.to_bits()
                    );
                }
                assert_eq!((batch.len(), ruled_batch.len()), (5, 5));
            }
        }
        // The declaration is scoped and restored, also after a panic.
        assert!(std::panic::catch_unwind(|| with_dense_workers(4, || panic!("inside"))).is_err());
        assert_eq!(DENSE_WORKERS.with(core::cell::Cell::get), 1);
        assert_eq!(
            with_dense_workers(4, || DENSE_WORKERS.with(core::cell::Cell::get)),
            4
        );
    }
}
