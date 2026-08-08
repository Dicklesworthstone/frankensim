//! Modal vibroacoustic coupling: structure <-> interior cavity <->
//! exterior radiation (bead frankensim-fsim-vibroacoustic-wgkq7,
//! musical-acoustics program) — fs-couple's first numerical engine.
//!
//! FORMULATION (frequency domain, `e^{-i omega t}`, matching
//! `fs_bem::helmholtz`): the coupled problem is projected onto
//! mass-normalized structural modes `phi_r` (unit modal mass,
//! `phi^T M phi = 1`, as produced by `fs_modal`) and rigid-wall cavity
//! pressure modes `psi_q` with norms `Lambda_q = INT psi_q^2 dV`. With
//! structural modal displacement `b_r` and cavity modal pressure `a_q`:
//!
//! structure row r:
//!   `[omega_r^2 (1 - i eta_s) - omega^2] b_r - SUM_q C_rq a_q
//!      - i omega SUM_s Zm_rs b_s = F_r`
//! cavity row q:
//!   `[omega_q^2 (1 - i eta_a) - omega^2] a_q
//!      - (rho0 c0^2 / Lambda_q) omega^2 SUM_r C_rq b_r = 0`
//!
//! (Under `e^{-i omega t}` the DISSIPATIVE hysteretic stiffness is
//! `k (1 - i eta)`; the opposite sign pumps energy, and the power
//! audit below is the executable proof.)
//!
//! where `C_rq = INT_S phi_r psi_q dA` over the interface with
//! positive structural deflection pointing AWAY from the cavity
//! (toward the exterior — the convention that composes directly with
//! the BEM's outward panel velocities; the sealed-cavity limit
//! `omega_q = 0` then gives `a = -(rho0 c0^2 / Lambda) C^T b`:
//! deflection away from the cavity rarefies it). One convention must
//! be shared by BOTH rows; an inconsistent flip between them is the
//! "sign-flipped interface normal" mutation, caught by the
//! interface-equals-cavity audit residual, not the (tautological)
//! input balance. `Zm` is the exterior radiation impedance matrix
//! projected onto the structural modes
//! ([`project_radiation_impedance`]; produced upstream by
//! `fs_bem::helmholtz::radiation_impedance_matrix`). Hysteretic loss
//! factors `eta_s`, `eta_a` model structural and cavity dissipation.
//!
//! POWER ACCOUNTING (the audit surface this module exists to make
//! real): at every solved frequency the response carries the complete
//! per-channel steady-state power breakdown —
//! `input = structural + interface + radiated` and
//! `interface = cavity` — with both residuals reported, never hidden.
//! The casebook integrates these over one period into an fs-couple
//! [`crate::WindowAuditReport`].
//!
//! TRUNCATION HONESTY: modal truncation error is OBSERVED, not
//! assumed. [`VibroacousticModel::frf_with_convergence`] re-solves on
//! the half-retained bases and refuses to report a converged value
//! when the relative delta exceeds the caller's tolerance.
//!
//! Determinism: dense assembly and solves in fixed traversal order via
//! `fs_la::eigen_complex`; repeat solves are bitwise identical.
//!
//! Deferred with recorded triggers (see CONTRACT): time-domain
//! realization (vector-fitting bead pvv40), non-rectangular cavities
//! beyond the lumped Helmholtz mode (numeric Laplacian cavity bases
//! plug into the same [`CavityModes`] carrier), and modal-density /
//! statistical-energy regimes.

use fs_la::eigen_complex::{eig, lu_complex};
use fs_math::c64::C64;

/// Typed refusals with stable `FS-COUPLE-VIBRO-*` codes.
#[derive(Debug, Clone, PartialEq)]
pub enum VibroError {
    /// A parameter is non-finite, non-positive, or out of range.
    BadParameter {
        /// Which parameter refused.
        what: &'static str,
    },
    /// Input lengths disagree.
    ShapeMismatch {
        /// What disagreed.
        what: &'static str,
    },
    /// A basis is empty.
    EmptyBasis {
        /// Which basis.
        what: &'static str,
    },
    /// The coupled dense system was singular at the requested frequency.
    Singular,
    /// The eigenvalue path failed to converge.
    EigenFailure,
    /// The truncation convergence check exceeded the caller tolerance.
    TruncationNotConverged {
        /// Measured relative delta between full and halved bases.
        relative_delta: f64,
        /// The caller's tolerance.
        tolerance: f64,
    },
}

impl core::fmt::Display for VibroError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            VibroError::BadParameter { what } => {
                write!(f, "FS-COUPLE-VIBRO-BAD-PARAMETER: {what}")
            }
            VibroError::ShapeMismatch { what } => {
                write!(f, "FS-COUPLE-VIBRO-SHAPE-MISMATCH: {what}")
            }
            VibroError::EmptyBasis { what } => {
                write!(f, "FS-COUPLE-VIBRO-EMPTY-BASIS: {what}")
            }
            VibroError::Singular => {
                write!(f, "FS-COUPLE-VIBRO-SINGULAR: coupled dense solve refused")
            }
            VibroError::EigenFailure => {
                write!(f, "FS-COUPLE-VIBRO-EIGEN-FAILURE: eigenvalue path refused")
            }
            VibroError::TruncationNotConverged {
                relative_delta,
                tolerance,
            } => write!(
                f,
                "FS-COUPLE-VIBRO-TRUNCATION-NOT-CONVERGED: relative delta \
                 {relative_delta:.3e} exceeds tolerance {tolerance:.3e}"
            ),
        }
    }
}

impl std::error::Error for VibroError {}

/// Mass-normalized structural modes restricted to the coupling
/// interface: `shapes[r][i]` is the transverse component of mode `r`
/// at interface sample point `i` (unit modal mass, `fs_modal`
/// convention `phi^T M phi = 1`).
#[derive(Debug, Clone)]
pub struct StructuralModes {
    /// Natural angular frequencies [rad/s], ascending.
    pub omegas: Vec<f64>,
    /// Interface-sampled mode shapes, one inner vector per mode.
    pub shapes: Vec<Vec<f64>>,
    /// Hysteretic structural loss factor eta_s (>= 0).
    pub loss_factor: f64,
}

/// Rigid-wall cavity pressure modes restricted to the same interface
/// sample points: `interface[q][i]` is `psi_q` at point `i`.
#[derive(Debug, Clone)]
pub struct CavityModes {
    /// Cavity natural angular frequencies [rad/s] (the bulk-compliance
    /// `omega = 0` mode of a sealed cavity is admitted).
    pub omegas: Vec<f64>,
    /// Mode norms `Lambda_q = INT psi_q^2 dV` [m^3].
    pub lambdas: Vec<f64>,
    /// Interface-sampled pressure mode shapes.
    pub interface: Vec<Vec<f64>>,
    /// Hysteretic acoustic loss factor eta_a (>= 0).
    pub loss_factor: f64,
    /// Ambient density rho0 [kg/m^3].
    pub rho0: f64,
    /// Sound speed c0 [m/s].
    pub c0: f64,
}

/// Acoustic medium (SI): ambient density and sound speed.
#[derive(Debug, Clone, Copy)]
pub struct AcousticMedium {
    /// Ambient density rho0 [kg/m^3].
    pub rho0: f64,
    /// Sound speed c0 [m/s].
    pub c0: f64,
}

impl AcousticMedium {
    /// Air at roughly 20 degC (matches `fs_bem::helmholtz::Medium::air`).
    /// These constants are `fs_material::gas::GasState` evaluated at
    /// (293.15 K, 101325 Pa) — parameterized studies should DERIVE the
    /// medium from that first-principles primitive (any temperature,
    /// pressure, or gas) instead of using this fixed convenience; the
    /// casebook asserts the equivalence.
    #[must_use]
    pub const fn air() -> AcousticMedium {
        AcousticMedium {
            rho0: 1.204,
            c0: 343.0,
        }
    }
}

/// Analytic rigid-wall modes of a rectangular cavity `lx x ly x lz`,
/// sampled on the interface plane `z = 0` at the given `(x, y)` points:
/// `psi_lmn = cos(l pi x / lx) cos(m pi y / ly) cos(n pi z / lz)`,
/// `omega = c0 pi sqrt((l/lx)^2 + (m/ly)^2 + (n/lz)^2)`,
/// `Lambda = V / (e_l e_m e_n)` with `e_0 = 1` and `e_k = 2` for
/// `k >= 1`. Modes are returned ascending in `omega` (ties broken by
/// `(l, m, n)` lexicographic order), truncated to `count`.
///
/// # Errors
/// [`VibroError`] on non-finite/non-positive geometry, zero count, or
/// empty interface.
pub fn rectangular_cavity_modes(
    lx: f64,
    ly: f64,
    lz: f64,
    medium: AcousticMedium,
    loss_factor: f64,
    count: usize,
    interface_points: &[[f64; 2]],
) -> Result<CavityModes, VibroError> {
    let AcousticMedium { rho0, c0 } = medium;
    for (value, what) in [
        (lx, "cavity dimension lx must be positive and finite"),
        (ly, "cavity dimension ly must be positive and finite"),
        (lz, "cavity dimension lz must be positive and finite"),
        (c0, "sound speed must be positive and finite"),
        (rho0, "density must be positive and finite"),
    ] {
        if !(value > 0.0 && value.is_finite()) {
            return Err(VibroError::BadParameter { what });
        }
    }
    if !(loss_factor >= 0.0 && loss_factor.is_finite()) {
        return Err(VibroError::BadParameter {
            what: "cavity loss factor must be finite and non-negative",
        });
    }
    if count == 0 {
        return Err(VibroError::BadParameter {
            what: "cavity mode count must be positive",
        });
    }
    if interface_points.is_empty() {
        return Err(VibroError::EmptyBasis {
            what: "cavity interface sample points",
        });
    }
    let pi = core::f64::consts::PI;
    // Enumerate a generous index box, sort by omega, truncate.
    // Enumeration box with an explicit SUFFICIENCY proof instead of an
    // isotropic-count heuristic (which silently dropped the lowest
    // modes of elongated cavities — executed counterexample in the
    // tests): any mode OUTSIDE the `0..=max_index` box has, along the
    // violating axis alone, `omega >= c0 pi (max_index + 1) / L_axis`.
    // The box is therefore complete for the lowest `count` modes once
    // that single-axis exclusion bound exceeds the count-th smallest
    // frequency found INSIDE the box, for every axis (the binding axis
    // is the longest one); double the box and retry otherwise.
    let mut max_index = 2usize;
    let mut all: Vec<(f64, usize, usize, usize)>;
    loop {
        all = Vec::with_capacity((max_index + 1).pow(3));
        for l in 0..=max_index {
            for m in 0..=max_index {
                for n in 0..=max_index {
                    let k2 =
                        (l as f64 / lx).powi(2) + (m as f64 / ly).powi(2) + (n as f64 / lz).powi(2);
                    all.push((c0 * pi * k2.sqrt(), l, m, n));
                }
            }
        }
        all.sort_by(|p, q| {
            p.0.partial_cmp(&q.0)
                .expect("finite frequencies")
                .then(p.1.cmp(&q.1))
                .then(p.2.cmp(&q.2))
                .then(p.3.cmp(&q.3))
        });
        if all.len() >= count {
            let omega_cut = all[count - 1].0;
            let excluded_floor = c0 * pi * (max_index + 1) as f64 / lx.max(ly).max(lz);
            if excluded_floor > omega_cut {
                break;
            }
        }
        max_index *= 2;
    }
    all.truncate(count);
    let volume = lx * ly * lz;
    let eps = |k: usize| if k == 0 { 1.0 } else { 2.0 };
    let mut omegas = Vec::with_capacity(count);
    let mut lambdas = Vec::with_capacity(count);
    let mut interface = Vec::with_capacity(count);
    for &(omega, l, m, n) in &all {
        omegas.push(omega);
        lambdas.push(volume / (eps(l) * eps(m) * eps(n)));
        let mut row = Vec::with_capacity(interface_points.len());
        for &[x, y] in interface_points {
            // psi on the z = 0 face: cos(n pi 0 / lz) = 1.
            row.push(det_cos(l as f64 * pi * x / lx) * det_cos(m as f64 * pi * y / ly));
        }
        interface.push(row);
    }
    Ok(CavityModes {
        omegas,
        lambdas,
        interface,
        loss_factor,
        rho0,
        c0,
    })
}

fn det_cos(x: f64) -> f64 {
    fs_math::det::cos(x)
}

/// Lumped Helmholtz-resonator mode for a cavity of volume `V` vented
/// through a cylindrical neck (radius `a`, physical length `l_neck`):
/// `omega_h = c0 sqrt(S / (V l_eff))` with the unflanged-aperture end
/// correction `l_eff = l_neck + 2 (8 / 3 pi) a` — the Laplace-BEM
/// closed-body recipe the epic's aperture pilot validated to +1.8%.
/// Returned as a one-mode [`CavityModes`] with uniform interior
/// pressure (`psi = 1`, `Lambda = V`): the resonator IS the cavity's
/// lowest mode in the modal-coupling picture.
///
/// # Errors
/// [`VibroError::BadParameter`] on non-positive geometry.
pub fn helmholtz_resonator_mode(
    volume: f64,
    neck_radius: f64,
    neck_length: f64,
    medium: AcousticMedium,
    loss_factor: f64,
    interface_count: usize,
) -> Result<CavityModes, VibroError> {
    let AcousticMedium { rho0, c0 } = medium;
    for (value, what) in [
        (volume, "resonator volume must be positive and finite"),
        (neck_radius, "neck radius must be positive and finite"),
        (c0, "sound speed must be positive and finite"),
        (rho0, "density must be positive and finite"),
    ] {
        if !(value > 0.0 && value.is_finite()) {
            return Err(VibroError::BadParameter { what });
        }
    }
    if !(neck_length >= 0.0 && neck_length.is_finite()) {
        return Err(VibroError::BadParameter {
            what: "neck length must be finite and non-negative",
        });
    }
    if !(loss_factor >= 0.0 && loss_factor.is_finite()) {
        return Err(VibroError::BadParameter {
            what: "resonator loss factor must be finite and non-negative",
        });
    }
    if interface_count == 0 {
        return Err(VibroError::EmptyBasis {
            what: "resonator interface sample points",
        });
    }
    let pi = core::f64::consts::PI;
    let area = pi * neck_radius * neck_radius;
    let l_eff = neck_length + 2.0 * (8.0 / (3.0 * pi)) * neck_radius;
    let omega = c0 * (area / (volume * l_eff)).sqrt();
    Ok(CavityModes {
        omegas: vec![omega],
        lambdas: vec![volume],
        interface: vec![vec![1.0; interface_count]],
        loss_factor,
        rho0,
        c0,
    })
}

/// Assemble the modal coupling matrix `C_rq = INT_S phi_r psi_q dA` by
/// the interface quadrature `SUM_i phi_r(i) psi_q(i) A_i`, row-major
/// `n_s x n_a`. The deflection sign convention (positive AWAY from the
/// cavity, matching the module formulation) is carried by the SIGN of
/// the structural shapes; use one convention for every projection
/// touching the same interface.
///
/// # Errors
/// [`VibroError::ShapeMismatch`] when the sample counts disagree.
pub fn assemble_coupling(
    structure: &StructuralModes,
    cavity: &CavityModes,
    node_areas: &[f64],
) -> Result<Vec<f64>, VibroError> {
    let n_i = node_areas.len();
    if structure.shapes.iter().any(|s| s.len() != n_i)
        || cavity.interface.iter().any(|s| s.len() != n_i)
    {
        return Err(VibroError::ShapeMismatch {
            what: "interface sample counts must match node_areas",
        });
    }
    let n_s = structure.shapes.len();
    let n_a = cavity.interface.len();
    let mut c = vec![0.0f64; n_s * n_a];
    for (r, phi) in structure.shapes.iter().enumerate() {
        for (q, psi) in cavity.interface.iter().enumerate() {
            let mut acc = 0.0;
            for i in 0..n_i {
                acc += phi[i] * psi[i] * node_areas[i];
            }
            c[r * n_a + q] = acc;
        }
    }
    Ok(c)
}

/// Project a panel-space radiation impedance matrix (`p = Z v`,
/// row-major over wet panels, from
/// `fs_bem::helmholtz::radiation_impedance_matrix`) onto the
/// structural modes: `Zm_rs = SUM_ij phi_r(i) A_i Z_ij phi_s(j)`, so
/// the modal radiation force is `f_r = -SUM_s Zm_rs v_s` with
/// `v_s = -i omega b_s`.
///
/// # Errors
/// [`VibroError::ShapeMismatch`] when panel counts disagree.
pub fn project_radiation_impedance(
    z_panel: &[C64],
    panel_areas: &[f64],
    shapes_at_panels: &[Vec<f64>],
) -> Result<Vec<C64>, VibroError> {
    let n_p = panel_areas.len();
    if z_panel.len() != n_p * n_p || shapes_at_panels.iter().any(|s| s.len() != n_p) {
        return Err(VibroError::ShapeMismatch {
            what: "radiation panel counts must agree between Z, areas, and shapes",
        });
    }
    let n_s = shapes_at_panels.len();
    let mut zm = vec![C64::ZERO; n_s * n_s];
    // t[j][s] = SUM over columns applied to shapes first: (Z * phi_s)(i).
    for s in 0..n_s {
        for i in 0..n_p {
            let mut acc = C64::ZERO;
            for j in 0..n_p {
                acc = acc + z_panel[i * n_p + j].scale(shapes_at_panels[s][j]);
            }
            let weight = panel_areas[i];
            for (r, shape_r) in shapes_at_panels.iter().enumerate() {
                zm[r * n_s + s] = zm[r * n_s + s] + acc.scale(shape_r[i] * weight);
            }
        }
    }
    Ok(zm)
}

/// Per-frequency steady-state power breakdown [W], with the two audit
/// residuals reported rather than hidden.
#[derive(Debug, Clone, Copy)]
pub struct PowerBreakdown {
    /// Input power `1/2 Re SUM F_r conj(v_r)`.
    pub input: f64,
    /// Structural hysteretic dissipation.
    pub structural: f64,
    /// Power flowing through the interface into the cavity.
    pub interface: f64,
    /// Cavity hysteretic dissipation.
    pub cavity: f64,
    /// Exterior radiated power `1/2 Re (v^H Zm v)` (zero without a
    /// radiation matrix).
    pub radiated: f64,
    /// `input - structural - interface - radiated`, relative to input.
    pub balance_residual: f64,
    /// `interface - cavity`, relative to max(interface, cavity).
    pub interface_residual: f64,
}

/// A solved coupled response at one frequency.
#[derive(Debug, Clone)]
pub struct VibroResponse {
    /// Angular frequency [rad/s].
    pub omega: f64,
    /// Structural modal displacements.
    pub b: Vec<C64>,
    /// Cavity modal pressures.
    pub a: Vec<C64>,
    /// Steady-state power accounting.
    pub power: PowerBreakdown,
}

/// A response whose modal-truncation error was observed by re-solving
/// on the half-retained bases.
#[derive(Debug, Clone)]
pub struct ConvergedResponse {
    /// The full-basis response.
    pub response: VibroResponse,
    /// Relative input-power delta between full and halved bases.
    pub truncation_delta: f64,
    /// Mode counts used: (structure full, cavity full, structure half,
    /// cavity half).
    pub mode_counts: (usize, usize, usize, usize),
}

/// The assembled modal vibroacoustic model.
#[derive(Debug, Clone)]
pub struct VibroacousticModel {
    omegas_s: Vec<f64>,
    eta_s: f64,
    omegas_a: Vec<f64>,
    lambdas: Vec<f64>,
    eta_a: f64,
    rho_c2: f64,
    coupling: Vec<f64>,
    z_modal: Option<Vec<C64>>,
}

impl VibroacousticModel {
    /// Assemble from structural/cavity bases and a coupling matrix
    /// (row-major `n_s x n_a`; see [`assemble_coupling`]). An optional
    /// modal radiation impedance (row-major `n_s x n_s`; see
    /// [`project_radiation_impedance`]) attaches exterior radiation.
    ///
    /// # Errors
    /// [`VibroError`] on empty bases, shape mismatches, or non-finite
    /// parameters.
    pub fn try_new(
        structure: &StructuralModes,
        cavity: &CavityModes,
        coupling: Vec<f64>,
        z_modal: Option<Vec<C64>>,
    ) -> Result<Self, VibroError> {
        let n_s = structure.omegas.len();
        let n_a = cavity.omegas.len();
        if n_s == 0 || structure.shapes.len() != n_s {
            return Err(if n_s == 0 {
                VibroError::EmptyBasis {
                    what: "structural modes",
                }
            } else {
                VibroError::ShapeMismatch {
                    what: "structural omegas and shapes must agree",
                }
            });
        }
        if n_a == 0 || cavity.interface.len() != n_a || cavity.lambdas.len() != n_a {
            return Err(if n_a == 0 {
                VibroError::EmptyBasis {
                    what: "cavity modes",
                }
            } else {
                VibroError::ShapeMismatch {
                    what: "cavity omegas, lambdas, and interface shapes must agree",
                }
            });
        }
        if coupling.len() != n_s * n_a {
            return Err(VibroError::ShapeMismatch {
                what: "coupling matrix must be n_s x n_a",
            });
        }
        if let Some(z) = &z_modal
            && z.len() != n_s * n_s
        {
            return Err(VibroError::ShapeMismatch {
                what: "modal radiation impedance must be n_s x n_s",
            });
        }
        for &w in structure.omegas.iter().chain(cavity.omegas.iter()) {
            if !(w >= 0.0 && w.is_finite()) {
                return Err(VibroError::BadParameter {
                    what: "natural frequencies must be finite and non-negative",
                });
            }
        }
        for &l in &cavity.lambdas {
            if !(l > 0.0 && l.is_finite()) {
                return Err(VibroError::BadParameter {
                    what: "cavity mode norms Lambda must be positive and finite",
                });
            }
        }
        for &(eta, what) in &[
            (
                structure.loss_factor,
                "structural loss factor must be finite and non-negative",
            ),
            (
                cavity.loss_factor,
                "cavity loss factor must be finite and non-negative",
            ),
        ] {
            if !(eta >= 0.0 && eta.is_finite()) {
                return Err(VibroError::BadParameter { what });
            }
        }
        if coupling.iter().any(|c| !c.is_finite()) {
            return Err(VibroError::BadParameter {
                what: "coupling entries must be finite",
            });
        }
        Ok(Self {
            omegas_s: structure.omegas.clone(),
            eta_s: structure.loss_factor,
            omegas_a: cavity.omegas.clone(),
            lambdas: cavity.lambdas.clone(),
            eta_a: cavity.loss_factor,
            rho_c2: cavity.rho0 * cavity.c0 * cavity.c0,
            coupling,
            z_modal,
        })
    }

    /// Structural mode count.
    #[must_use]
    pub fn structure_count(&self) -> usize {
        self.omegas_s.len()
    }

    /// Cavity mode count.
    #[must_use]
    pub fn cavity_count(&self) -> usize {
        self.omegas_a.len()
    }

    /// Undamped coupled natural frequencies [rad/s], ascending: the
    /// roots of the coupled polynomial problem, obtained by the exact
    /// linearization `A z = x B z` in `x = omega^2` with
    /// `A = [[Omega_s^2, -C], [0, Omega_a^2]]`,
    /// `B = [[I, 0], [-(rho0 c0^2 / Lambda) C^T, I]]` (`B` is unit
    /// lower-triangular, so `B^{-1} A` is formed exactly) and a dense
    /// eigensolve. Radiation and loss factors are ignored — this is
    /// the conservative pencil the two-oscillator closed form pins.
    ///
    /// # Errors
    /// [`VibroError::EigenFailure`] if the dense eigensolve refuses.
    pub fn undamped_natural_frequencies(&self) -> Result<Vec<f64>, VibroError> {
        let n_s = self.omegas_s.len();
        let n_a = self.omegas_a.len();
        let n = n_s + n_a;
        // From the undamped rows: omega_r^2 b - C a = x b and
        // omega_q^2 a = x (a + s C^T b) with s = rho0 c0^2 / Lambda_q,
        // i.e. A z = x B z with A = [[Omega_s^2, -C], [0, Omega_a^2]],
        // B = [[I, 0], [s C^T, I]]. B is unit lower-triangular, so
        // B^{-1} = [[I, 0], [-s C^T, I]] exactly and
        // M = B^{-1} A has bottom rows -s C^T [Omega_s^2, -C]
        // + [0, Omega_a^2]. Two-mode check: trace(M) =
        // omega_s^2 + omega_a^2 + s C^2, det(M) = omega_s^2 omega_a^2
        // — exactly the closed-form split polynomial.
        let mut m = vec![C64::ZERO; n * n];
        for r in 0..n_s {
            m[r * n + r] = C64::from_re(self.omegas_s[r] * self.omegas_s[r]);
            for q in 0..n_a {
                m[r * n + (n_s + q)] = C64::from_re(-self.coupling[r * n_a + q]);
            }
        }
        for q in 0..n_a {
            let scale = self.rho_c2 / self.lambdas[q];
            for col_s in 0..n_s {
                // (-s C^T Omega_s^2)_{q, col_s} = -s C_{col_s, q} omega_s^2.
                m[(n_s + q) * n + col_s] = C64::from_re(
                    -scale
                        * self.coupling[col_s * n_a + q]
                        * self.omegas_s[col_s]
                        * self.omegas_s[col_s],
                );
            }
            for col_a in 0..n_a {
                // +s (C^T C)_{q, col_a} from (-s C^T)(-C).
                let mut acc = 0.0;
                for r in 0..n_s {
                    acc += self.coupling[r * n_a + q] * self.coupling[r * n_a + col_a];
                }
                let mut entry = scale * acc;
                if col_a == q {
                    entry += self.omegas_a[q] * self.omegas_a[q];
                }
                m[(n_s + q) * n + (n_s + col_a)] = C64::from_re(entry);
            }
        }
        let values = eig(&m, n).map_err(|_| VibroError::EigenFailure)?;
        let mut roots: Vec<f64> = values
            .iter()
            .map(|x| {
                // Physical coupled squares are real and non-negative;
                // tiny imaginary parts are eigensolver roundoff.
                x.re.max(0.0).sqrt()
            })
            .collect();
        roots.sort_by(|p, q| p.partial_cmp(q).expect("finite frequencies"));
        Ok(roots)
    }

    /// Solve the coupled frequency response at `omega` for the given
    /// structural modal forces, with the complete power breakdown.
    ///
    /// # Errors
    /// [`VibroError`] on bad frequency, force-shape mismatch, or a
    /// singular coupled system.
    pub fn frf(&self, omega: f64, modal_force: &[C64]) -> Result<VibroResponse, VibroError> {
        self.frf_truncated(omega, modal_force, self.omegas_s.len(), self.omegas_a.len())
    }

    /// [`Self::frf`] restricted to the first `keep_s` structural and
    /// `keep_a` cavity modes — the observability hook the convergence
    /// check uses.
    ///
    /// # Errors
    /// As for [`Self::frf`], plus a refusal on zero retained modes.
    pub fn frf_truncated(
        &self,
        omega: f64,
        modal_force: &[C64],
        keep_s: usize,
        keep_a: usize,
    ) -> Result<VibroResponse, VibroError> {
        if !(omega > 0.0 && omega.is_finite()) {
            return Err(VibroError::BadParameter {
                what: "frequency must be positive and finite",
            });
        }
        let n_s_full = self.omegas_s.len();
        let n_a_full = self.omegas_a.len();
        if modal_force.len() != n_s_full {
            return Err(VibroError::ShapeMismatch {
                what: "modal force length must equal the structural mode count",
            });
        }
        if keep_s == 0 || keep_s > n_s_full || keep_a == 0 || keep_a > n_a_full {
            return Err(VibroError::BadParameter {
                what: "retained mode counts must be in 1..=basis size",
            });
        }
        let n = keep_s + keep_a;
        let w2 = omega * omega;
        let mut mat = vec![C64::ZERO; n * n];
        for r in 0..keep_s {
            let wr2 = self.omegas_s[r] * self.omegas_s[r];
            // Under e^{-i omega t} the DISSIPATIVE hysteretic stiffness
            // is k (1 - i eta); the +i eta sign feeds energy in (the
            // power audit below is the executable proof).
            mat[r * n + r] = C64::new(wr2 - w2, -wr2 * self.eta_s);
            for q in 0..keep_a {
                mat[r * n + (keep_s + q)] = C64::from_re(-self.coupling[r * n_a_full + q]);
            }
            if let Some(z) = &self.z_modal {
                // Radiation force -Zm v = i omega Zm b moves to the
                // left as -i omega Zm b.
                for s in 0..keep_s {
                    let zc = z[r * n_s_full + s];
                    mat[r * n + s] = mat[r * n + s] - C64::new(0.0, omega) * zc;
                }
            }
        }
        for q in 0..keep_a {
            let wq2 = self.omegas_a[q] * self.omegas_a[q];
            let row = keep_s + q;
            mat[row * n + row] = C64::new(wq2 - w2, -wq2 * self.eta_a);
            let scale = self.rho_c2 / self.lambdas[q] * w2;
            for r in 0..keep_s {
                mat[row * n + r] = C64::from_re(-scale * self.coupling[r * n_a_full + q]);
            }
        }
        let mut rhs = vec![C64::ZERO; n];
        rhs[..keep_s].copy_from_slice(&modal_force[..keep_s]);
        let lu = lu_complex(&mat, n).map_err(|_| VibroError::Singular)?;
        lu.solve(&mut rhs);
        let b: Vec<C64> = rhs[..keep_s].to_vec();
        let a: Vec<C64> = rhs[keep_s..].to_vec();

        // Power accounting (e^{-i omega t}: v_r = -i omega b_r).
        let vel = |b_r: C64| C64::new(0.0, -omega) * b_r;
        let mut input = 0.0;
        let mut structural = 0.0;
        let mut interface = 0.0;
        for r in 0..keep_s {
            let v = vel(b[r]);
            input += 0.5 * (modal_force[r] * v.conj()).re;
            let wr2 = self.omegas_s[r] * self.omegas_s[r];
            structural += 0.5 * omega * self.eta_s * wr2 * b[r].norm_sq();
            // Pressure force on structure = +SUM_q C_rq a_q; power INTO
            // the cavity is minus the power that force delivers to the
            // structure.
            let mut f_p = C64::ZERO;
            for (q, aq) in a.iter().enumerate() {
                f_p = f_p + aq.scale(self.coupling[r * n_a_full + q]);
            }
            interface -= 0.5 * (f_p * v.conj()).re;
        }
        let mut cavity = 0.0;
        for (q, aq) in a.iter().enumerate() {
            // Cavity oscillator mass Lambda_q / (rho0 c0^2 omega^2)
            // in the pressure coordinate: hysteretic dissipation
            // 1/2 omega eta_a (Lambda/(rho c^2)) (omega_q/omega)^2 |a_q|^2.
            let wq2 = self.omegas_a[q] * self.omegas_a[q];
            cavity += 0.5 * self.eta_a * self.lambdas[q] / self.rho_c2 * wq2 / omega * aq.norm_sq();
        }
        let mut radiated = 0.0;
        if let Some(z) = &self.z_modal {
            for r in 0..keep_s {
                let vr = vel(b[r]);
                for s in 0..keep_s {
                    let vs = vel(b[s]);
                    radiated += 0.5 * (vr.conj() * z[r * n_s_full + s] * vs).re;
                }
            }
        }
        let balance = input - structural - interface - radiated;
        let denom = input.abs().max(1e-300);
        let iface_denom = interface.abs().max(cavity.abs()).max(1e-300);
        Ok(VibroResponse {
            omega,
            b,
            a,
            power: PowerBreakdown {
                input,
                structural,
                interface,
                cavity,
                radiated,
                balance_residual: balance / denom,
                interface_residual: (interface - cavity) / iface_denom,
            },
        })
    }

    /// Solve with OBSERVED truncation error: the response is
    /// re-computed on the half-retained bases and the relative
    /// input-power delta must sit inside `tolerance`, else the solve
    /// refuses (`FS-COUPLE-VIBRO-TRUNCATION-NOT-CONVERGED`) rather
    /// than reporting an unconverged value as converged.
    ///
    /// # Errors
    /// As for [`Self::frf`], plus the truncation refusal.
    pub fn frf_with_convergence(
        &self,
        omega: f64,
        modal_force: &[C64],
        tolerance: f64,
    ) -> Result<ConvergedResponse, VibroError> {
        if !(tolerance > 0.0 && tolerance.is_finite()) {
            return Err(VibroError::BadParameter {
                what: "truncation tolerance must be positive and finite",
            });
        }
        let n_s = self.omegas_s.len();
        let n_a = self.omegas_a.len();
        let full = self.frf(omega, modal_force)?;
        let half_s = (n_s / 2).max(1);
        let half_a = (n_a / 2).max(1);
        let half = self.frf_truncated(omega, modal_force, half_s, half_a)?;
        let denom = full.power.input.abs().max(1e-300);
        let truncation_delta = (full.power.input - half.power.input).abs() / denom;
        if truncation_delta > tolerance {
            return Err(VibroError::TruncationNotConverged {
                relative_delta: truncation_delta,
                tolerance,
            });
        }
        Ok(ConvergedResponse {
            response: full,
            truncation_delta,
            mode_counts: (n_s, n_a, half_s, half_a),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Deterministic small model used across tests: n_s structure
    /// modes, n_a cavity modes, fixed synthetic data.
    fn synthetic_model(
        n_s: usize,
        n_a: usize,
        coupling_scale: f64,
        eta_s: f64,
        eta_a: f64,
    ) -> VibroacousticModel {
        let structure = StructuralModes {
            omegas: (0..n_s).map(|r| 100.0 + 60.0 * r as f64).collect(),
            shapes: (0..n_s).map(|r| vec![1.0 + 0.1 * r as f64]).collect(),
            loss_factor: eta_s,
        };
        let cavity = CavityModes {
            omegas: (0..n_a).map(|q| 140.0 + 90.0 * q as f64).collect(),
            lambdas: (0..n_a).map(|q| 0.02 + 0.005 * q as f64).collect(),
            interface: (0..n_a).map(|q| vec![1.0 - 0.05 * q as f64]).collect(),
            loss_factor: eta_a,
            rho0: 1.204,
            c0: 343.0,
        };
        let coupling: Vec<f64> = (0..n_s * n_a)
            .map(|i| coupling_scale * (1.0 + 0.3 * (i % 5) as f64))
            .collect();
        VibroacousticModel::try_new(&structure, &cavity, coupling, None).expect("model")
    }

    #[test]
    fn coupling_overlap_matches_closed_form() {
        // Plate mode sin(p pi x/a) sin(q pi y/b) against cavity mode
        // cos(l pi x/a) cos(m pi y/b): the overlap factorizes into 1D
        // integrals with the closed form
        // INT_0^L sin(p pi x/L) cos(l pi x/L) dx
        //   = (L p / pi) (1 - (-1)^(p+l)) / (p^2 - l^2)   for p != l,
        //   = 0 for p = l.
        let (a, b) = (0.42, 0.31);
        let one_d = |p: usize, l: usize, len: f64| -> f64 {
            if p == l {
                0.0
            } else {
                let pf = p as f64;
                let lf = l as f64;
                let parity = if (p + l).is_multiple_of(2) { 0.0 } else { 2.0 };
                len * pf * parity / (core::f64::consts::PI * (pf * pf - lf * lf))
            }
        };
        let pi = core::f64::consts::PI;
        let (nx, ny) = (200usize, 200usize);
        let (dx, dy) = (a / nx as f64, b / ny as f64);
        let mut points = Vec::new();
        for i in 0..nx {
            for j in 0..ny {
                points.push([(i as f64 + 0.5) * dx, (j as f64 + 0.5) * dy]);
            }
        }
        let areas = vec![dx * dy; points.len()];
        for &(p, q, l, m) in &[(1usize, 1usize, 0usize, 0usize), (2, 1, 1, 0), (3, 2, 2, 1)] {
            let phi: Vec<f64> = points
                .iter()
                .map(|&[x, y]| {
                    fs_math::det::sin(p as f64 * pi * x / a)
                        * fs_math::det::sin(q as f64 * pi * y / b)
                })
                .collect();
            let psi: Vec<f64> = points
                .iter()
                .map(|&[x, y]| det_cos(l as f64 * pi * x / a) * det_cos(m as f64 * pi * y / b))
                .collect();
            let structure = StructuralModes {
                omegas: vec![1.0],
                shapes: vec![phi],
                loss_factor: 0.0,
            };
            let cavity = CavityModes {
                omegas: vec![1.0],
                lambdas: vec![1.0],
                interface: vec![psi],
                loss_factor: 0.0,
                rho0: 1.204,
                c0: 343.0,
            };
            let c = assemble_coupling(&structure, &cavity, &areas).expect("coupling");
            let exact = one_d(p, l, a) * one_d(q, m, b);
            let err = (c[0] - exact).abs();
            let scale = exact.abs().max(1e-6);
            assert!(
                err < 1e-4 * scale.max(1e-3),
                "overlap ({p},{q})x({l},{m}): {} vs {exact} (err {err:.2e})",
                c[0]
            );
        }
        println!(
            "{{\"suite\":\"fs-couple-vibro\",\"case\":\"coupling-overlap\",\"verdict\":\"pass\"}}"
        );
    }

    #[test]
    #[allow(clippy::too_many_lines)] // one coherent split gate + mutation arm
    fn two_oscillator_split_matches_closed_form_and_dropped_coupling_is_caught() {
        use core::fmt::Write as _;
        // One plate mode + one cavity mode: the coupled squares solve
        // x^2 - x (ws^2 + wc^2 + K) + ws^2 wc^2 = 0 with
        // K = rho0 c0^2 C^2 / Lambda. The engine's undamped pencil must
        // reproduce both roots to 1e-8 relative across a coupling sweep;
        // the DROPPED-COUPLING mutation (C = 0) collapses the split to
        // the uncoupled frequencies and is caught by the same gate.
        let ws = 2.0 * core::f64::consts::PI * 180.0;
        let wc = 2.0 * core::f64::consts::PI * 210.0;
        let rho0 = 1.204;
        let c0 = 343.0;
        let lambda = 0.012;
        let mut rows = String::new();
        let mut first = true;
        for &c_val in &[1e-5, 1e-4, 1e-3, 3e-3, 1e-2] {
            let structure = StructuralModes {
                omegas: vec![ws],
                shapes: vec![vec![1.0]],
                loss_factor: 0.0,
            };
            let cavity = CavityModes {
                omegas: vec![wc],
                lambdas: vec![lambda],
                interface: vec![vec![1.0]],
                loss_factor: 0.0,
                rho0,
                c0,
            };
            let model =
                VibroacousticModel::try_new(&structure, &cavity, vec![c_val], None).expect("model");
            let freqs = model.undamped_natural_frequencies().expect("pencil");
            assert_eq!(freqs.len(), 2);
            let k = rho0 * c0 * c0 * c_val * c_val / lambda;
            let s = ws * ws + wc * wc + k;
            let disc = (s * s - 4.0 * ws * ws * wc * wc).sqrt();
            let x_lo = 0.5 * (s - disc);
            let x_hi = f64::midpoint(s, disc);
            let rel_lo = (freqs[0] - x_lo.sqrt()).abs() / x_lo.sqrt();
            let rel_hi = (freqs[1] - x_hi.sqrt()).abs() / x_hi.sqrt();
            assert!(
                rel_lo < 1e-8 && rel_hi < 1e-8,
                "split at C={c_val}: engine ({}, {}) vs closed form ({}, {}) rel ({rel_lo:.2e}, {rel_hi:.2e})",
                freqs[0],
                freqs[1],
                x_lo.sqrt(),
                x_hi.sqrt()
            );
            write!(
                rows,
                "{}{{\"c\":{c_val},\"lo_rel\":{rel_lo:.2e},\"hi_rel\":{rel_hi:.2e}}}",
                if first { "" } else { "," }
            )
            .expect("write");
            first = false;
        }
        // Strong coupling widens the split beyond the uncoupled gap.
        let strong = {
            let structure = StructuralModes {
                omegas: vec![ws],
                shapes: vec![vec![1.0]],
                loss_factor: 0.0,
            };
            let cavity = CavityModes {
                omegas: vec![wc],
                lambdas: vec![lambda],
                interface: vec![vec![1.0]],
                loss_factor: 0.0,
                rho0,
                c0,
            };
            VibroacousticModel::try_new(&structure, &cavity, vec![0.2], None)
                .expect("model")
                .undamped_natural_frequencies()
                .expect("pencil")
        };
        let dropped = {
            let structure = StructuralModes {
                omegas: vec![ws],
                shapes: vec![vec![1.0]],
                loss_factor: 0.0,
            };
            let cavity = CavityModes {
                omegas: vec![wc],
                lambdas: vec![lambda],
                interface: vec![vec![1.0]],
                loss_factor: 0.0,
                rho0,
                c0,
            };
            VibroacousticModel::try_new(&structure, &cavity, vec![0.0], None)
                .expect("model")
                .undamped_natural_frequencies()
                .expect("pencil")
        };
        let split_strong = strong[1] - strong[0];
        let split_dropped = dropped[1] - dropped[0];
        assert!(
            split_strong > 1.2 * (wc - ws),
            "strong coupling must widen the split: {split_strong} vs uncoupled {}",
            wc - ws
        );
        assert!(
            (split_dropped - (wc - ws)).abs() < 1e-8 * wc,
            "dropped coupling must collapse to the uncoupled gap"
        );
        println!(
            "{{\"suite\":\"fs-couple-vibro\",\"case\":\"two-oscillator-split\",\"rows\":[{rows}],\"verdict\":\"pass\"}}"
        );
    }

    #[test]
    fn light_fluid_added_mass_downshift_is_first_order() {
        // Perturbation limit: the structural root of the split
        // polynomial sits at x = ws^2 + K ws^2 / (ws^2 - wc^2) + O(K^2)
        // (a DOWNSHIFT for wc > ws). First-order falsifier: halving K
        // must reduce the defect against the linear formula by ~4x.
        let ws = 2.0 * core::f64::consts::PI * 120.0;
        let wc = 2.0 * core::f64::consts::PI * 400.0;
        let rho0 = 1.204;
        let c0 = 343.0;
        let lambda = 0.02;
        let defect = |c_val: f64| -> (f64, f64) {
            let structure = StructuralModes {
                omegas: vec![ws],
                shapes: vec![vec![1.0]],
                loss_factor: 0.0,
            };
            let cavity = CavityModes {
                omegas: vec![wc],
                lambdas: vec![lambda],
                interface: vec![vec![1.0]],
                loss_factor: 0.0,
                rho0,
                c0,
            };
            let freqs = VibroacousticModel::try_new(&structure, &cavity, vec![c_val], None)
                .expect("model")
                .undamped_natural_frequencies()
                .expect("pencil");
            let k = rho0 * c0 * c0 * c_val * c_val / lambda;
            let x_linear = ws * ws + k * ws * ws / (ws * ws - wc * wc);
            let x_engine = freqs[0] * freqs[0];
            ((x_engine - x_linear).abs(), k)
        };
        let (d1, k1) = defect(2e-3);
        // K scales as C^2, so C / sqrt(2) halves K and the O(K^2)
        // defect against the LINEAR formula must drop ~4x.
        let (d2, _) = defect(2e-3 / core::f64::consts::SQRT_2);
        // Downshift direction: engine root below the uncoupled one.
        let (base, _) = defect(0.0);
        assert!(base < 1e-9 * ws * ws, "zero coupling must be exact");
        assert!(d1 > 0.0 && d2 > 0.0, "finite-K defects must be nonzero");
        let ratio = d1 / d2;
        assert!(
            (3.0..5.0).contains(&ratio),
            "halving K must reduce the O(K^2) defect ~4x, got {ratio:.2} (d1 {d1:.3e}, d2 {d2:.3e}, K {k1:.3e})"
        );
        println!(
            "{{\"suite\":\"fs-couple-vibro\",\"case\":\"added-mass-first-order\",\"quartering\":{ratio:.2},\"verdict\":\"pass\"}}"
        );
    }

    #[test]
    fn energy_balance_is_exact_and_one_sided_normal_flip_is_caught() {
        use core::fmt::Write as _;
        // The per-frequency identities input = structural + interface
        // and interface = cavity are EXACT for the discrete system (the
        // derivation lives in the power-accounting code), so the
        // residuals must sit at solver roundoff across frequencies and
        // mode subsets. The MUTATION arm rebuilds the 2x2 system with
        // the interface normal flipped in the STRUCTURE row only (the
        // engine's shared coupling matrix makes that inconsistency
        // unrepresentable through the public API, so the arm is
        // hand-assembled here) — its balance residual must blow up.
        let model = synthetic_model(3, 4, 2e-3, 0.02, 0.008);
        let force: Vec<C64> = vec![C64::ONE, C64::new(0.3, 0.2), C64::new(-0.4, 0.1)];
        let mut rows = String::new();
        let mut first = true;
        for &omega in &[60.0, 101.0, 139.5, 220.0, 400.0] {
            let sol = model.frf(omega, &force).expect("frf");
            assert!(
                sol.power.balance_residual.abs() < 1e-10,
                "balance residual {:.3e} at omega {omega}",
                sol.power.balance_residual
            );
            assert!(
                sol.power.interface_residual.abs() < 1e-10,
                "interface residual {:.3e} at omega {omega}",
                sol.power.interface_residual
            );
            assert!(sol.power.input > 0.0 && sol.power.structural >= 0.0);
            write!(
                rows,
                "{}{{\"omega\":{omega},\"balance\":{:.2e},\"interface\":{:.2e}}}",
                if first { "" } else { "," },
                sol.power.balance_residual,
                sol.power.interface_residual
            )
            .expect("write");
            first = false;
        }
        // Hand-assembled 2x2 with a one-sided normal flip. NOTE (a
        // finding from executing this arm): the input-balance identity
        // closes TAUTOLOGICALLY for whatever system was assembled,
        // because the audit's interface power uses the same assembled
        // row — a globally flipped normal is merely the other (equally
        // consistent) convention. What a one-SIDED flip breaks is the
        // cross-row identity interface == cavity dissipation, so THAT
        // residual is the discriminating alarm.
        let (ws, wc, lambda, rho_c2, c_val, eta_s, eta_a, omega) = (
            100.0f64,
            140.0f64,
            0.02f64,
            1.204 * 343.0 * 343.0,
            2e-3f64,
            0.02f64,
            0.01f64,
            120.0f64,
        );
        let w2 = omega * omega;
        // Returns (balance_residual, interface_vs_cavity_residual).
        let solve2 = |c_structure_row: f64| -> (f64, f64) {
            let mat = [
                C64::new(ws * ws - w2, -ws * ws * eta_s),
                C64::from_re(-c_structure_row),
                C64::from_re(-rho_c2 / lambda * w2 * c_val),
                C64::new(wc * wc - w2, -wc * wc * eta_a),
            ];
            let lu = lu_complex(&mat, 2).expect("2x2");
            let mut rhs = [C64::ONE, C64::ZERO];
            lu.solve(&mut rhs);
            let (b, aq) = (rhs[0], rhs[1]);
            let v = C64::new(0.0, -omega) * b;
            let input = 0.5 * (C64::ONE * v.conj()).re;
            let structural = 0.5 * omega * eta_s * ws * ws * b.norm_sq();
            let interface = -0.5 * (aq.scale(c_structure_row) * v.conj()).re;
            let cavity = 0.5 * eta_a * lambda / rho_c2 * wc * wc / omega * aq.norm_sq();
            (
                (input - structural - interface) / input.abs().max(1e-300),
                (interface - cavity) / interface.abs().max(cavity.abs()).max(1e-300),
            )
        };
        let (bal_ok, iface_ok) = solve2(c_val);
        let (bal_flip, iface_flip) = solve2(-c_val);
        assert!(
            bal_ok.abs() < 1e-10 && iface_ok.abs() < 1e-10,
            "consistent normals must pass both residuals: {bal_ok:.3e}, {iface_ok:.3e}"
        );
        // The tautological balance stays green even for the mutation —
        // asserted so the finding stays executable...
        assert!(
            bal_flip.abs() < 1e-10,
            "the balance identity alone CANNOT see the one-sided flip"
        );
        // ...and the cross-row residual is the alarm that fires.
        assert!(
            iface_flip.abs() > 0.5,
            "one-sided normal flip must break interface == cavity: {iface_flip:.3e}"
        );
        let flipped = iface_flip.abs();
        println!(
            "{{\"suite\":\"fs-couple-vibro\",\"case\":\"energy-balance\",\"rows\":[{rows}],\"flipped_residual\":{flipped:.2e},\"verdict\":\"pass\"}}"
        );
    }

    #[test]
    fn structural_frf_is_reciprocal() {
        // The cavity-condensed structural dynamic stiffness is complex
        // symmetric, so H_rs = H_sr for unit modal forces.
        let model = synthetic_model(3, 4, 2e-3, 0.015, 0.005);
        let omega = 150.0;
        let unit = |r: usize| -> Vec<C64> {
            let mut f = vec![C64::ZERO; 3];
            f[r] = C64::ONE;
            f
        };
        for r in 0..3 {
            for s in (r + 1)..3 {
                let b_sr = model.frf(omega, &unit(r)).expect("frf").b[s];
                let b_rs = model.frf(omega, &unit(s)).expect("frf").b[r];
                let defect = (b_sr - b_rs).abs() / b_sr.abs().max(1e-300);
                assert!(
                    defect < 1e-10,
                    "reciprocity ({r},{s}): {b_sr:?} vs {b_rs:?} (defect {defect:.2e})"
                );
            }
        }
        println!("{{\"suite\":\"fs-couple-vibro\",\"case\":\"reciprocity\",\"verdict\":\"pass\"}}");
    }

    #[test]
    fn truncation_convergence_is_observed_not_assumed() {
        // Drive near the HIGHEST retained cavity mode: the halved basis
        // cannot represent it, so the observed delta is large and a
        // tight tolerance must REFUSE. The plain frf() path makes no
        // convergence claim — that asymmetry is the API contract, so
        // the "truncation-check disabled" mutation is unrepresentable
        // through frf_with_convergence.
        let model = synthetic_model(2, 6, 3e-3, 0.02, 0.01);
        let force = vec![C64::ONE, C64::ZERO];
        // Highest cavity mode of the synthetic model: 140 + 90*5 = 590.
        let omega_near_top = 589.0;
        let err = model
            .frf_with_convergence(omega_near_top, &force, 1e-9)
            .expect_err("tight tolerance near an unconverged frequency must refuse");
        let VibroError::TruncationNotConverged { relative_delta, .. } = err else {
            panic!("wrong refusal: {err}");
        };
        assert!(relative_delta > 1e-9);
        // Far below every truncated mode the delta is small and the
        // report carries it.
        let ok = model
            .frf_with_convergence(30.0, &force, 0.05)
            .expect("low-frequency solve converges");
        assert!(ok.truncation_delta < 0.05);
        assert_eq!(ok.mode_counts, (2, 6, 1, 3));
        println!(
            "{{\"suite\":\"fs-couple-vibro\",\"case\":\"truncation-honesty\",\"refused_delta\":{relative_delta:.2e},\"low_freq_delta\":{:.2e},\"verdict\":\"pass\"}}",
            ok.truncation_delta
        );
    }

    #[test]
    fn rectangular_cavity_and_resonator_pins() {
        // Cavity: first mode is the omega = 0 bulk mode with Lambda = V;
        // the next is c pi / max_dimension; psi samples match the cosine
        // closed form.
        let (lx, ly, lz) = (0.5, 0.4, 0.3);
        let points = [[0.1, 0.1], [0.25, 0.2], [0.4, 0.35]];
        let cav = rectangular_cavity_modes(lx, ly, lz, AcousticMedium::air(), 0.01, 6, &points)
            .expect("cavity modes");
        assert_eq!(cav.omegas.len(), 6);
        assert_eq!(cav.omegas[0].to_bits(), 0);
        assert!((cav.lambdas[0] - lx * ly * lz).abs() < 1e-15);
        let expected_second = 343.0 * core::f64::consts::PI / lx;
        assert!(
            (cav.omegas[1] - expected_second).abs() < 1e-9 * expected_second,
            "second mode {} vs {expected_second}",
            cav.omegas[1]
        );
        // Second mode is (1,0,0): psi = cos(pi x / lx), Lambda = V/2.
        for (i, &[x, _y]) in points.iter().enumerate() {
            let expected = det_cos(core::f64::consts::PI * x / lx);
            assert!((cav.interface[1][i] - expected).abs() < 1e-12);
        }
        assert!((cav.lambdas[1] - lx * ly * lz / 2.0).abs() < 1e-15);
        // Resonator: 1 liter, 5 mm neck radius, 10 mm neck length ->
        // f = omega / 2 pi ~ 112.5 Hz with the 2 * 8/(3 pi) a end
        // correction.
        let res = helmholtz_resonator_mode(1e-3, 0.005, 0.010, AcousticMedium::air(), 0.01, 4)
            .expect("resonator");
        let a = 0.005f64;
        let s_neck = core::f64::consts::PI * a * a;
        let l_eff = 0.010 + 2.0 * (8.0 / (3.0 * core::f64::consts::PI)) * a;
        let omega_expected = 343.0 * (s_neck / (1e-3 * l_eff)).sqrt();
        assert!(
            (res.omegas[0] - omega_expected).abs() < 1e-12 * omega_expected,
            "resonator omega {} vs {omega_expected}",
            res.omegas[0]
        );
        assert!((res.omegas[0] / (2.0 * core::f64::consts::PI) - 112.5).abs() < 1.0);
        println!(
            "{{\"suite\":\"fs-couple-vibro\",\"case\":\"cavity-resonator-pins\",\"resonator_hz\":{:.2},\"verdict\":\"pass\"}}",
            res.omegas[0] / (2.0 * core::f64::consts::PI)
        );
    }

    #[test]
    fn named_refusals_fire() {
        let model = synthetic_model(2, 2, 1e-3, 0.01, 0.01);
        let force = vec![C64::ONE, C64::ZERO];
        assert!(
            model
                .frf(-1.0, &force)
                .unwrap_err()
                .to_string()
                .contains("FS-COUPLE-VIBRO-BAD-PARAMETER")
        );
        assert!(
            model
                .frf(100.0, &[C64::ONE])
                .unwrap_err()
                .to_string()
                .contains("FS-COUPLE-VIBRO-SHAPE-MISMATCH")
        );
        assert!(matches!(
            rectangular_cavity_modes(0.0, 1.0, 1.0, AcousticMedium::air(), 0.0, 3, &[[0.1, 0.1]]),
            Err(VibroError::BadParameter { .. })
        ));
        assert!(matches!(
            helmholtz_resonator_mode(1e-3, -0.01, 0.01, AcousticMedium::air(), 0.0, 1),
            Err(VibroError::BadParameter { .. })
        ));
        let structure = StructuralModes {
            omegas: vec![],
            shapes: vec![],
            loss_factor: 0.0,
        };
        let cavity = CavityModes {
            omegas: vec![1.0],
            lambdas: vec![1.0],
            interface: vec![vec![1.0]],
            loss_factor: 0.0,
            rho0: 1.2,
            c0: 343.0,
        };
        assert!(matches!(
            VibroacousticModel::try_new(&structure, &cavity, vec![], None),
            Err(VibroError::EmptyBasis { .. })
        ));
    }

    #[test]
    fn repeat_solves_are_bitwise_identical() {
        let model = synthetic_model(3, 4, 2e-3, 0.02, 0.008);
        let force: Vec<C64> = vec![C64::ONE, C64::new(0.3, 0.2), C64::new(-0.4, 0.1)];
        let x = model.frf(151.7, &force).expect("frf");
        let y = model.frf(151.7, &force).expect("frf");
        for (p, q) in x.b.iter().zip(y.b.iter()).chain(x.a.iter().zip(y.a.iter())) {
            assert!(p.re.to_bits() == q.re.to_bits() && p.im.to_bits() == q.im.to_bits());
        }
        assert_eq!(
            x.power.balance_residual.to_bits(),
            y.power.balance_residual.to_bits()
        );
    }
}

#[cfg(test)]
mod review_regressions {
    use super::*;
    use fs_la::eigen_complex::det_complex;

    #[test]
    fn elongated_cavity_keeps_its_lowest_axial_modes() {
        // REGRESSION (fresh-eyes review, executed counterexample): the
        // old isotropic index-box heuristic returned (0,0,1)/(0,1,0)
        // modes at ~10776 rad/s as the 5th/6th modes of a 1.0 x 0.1 x
        // 0.1 duct, silently dropping the true (4,0,0) and (5,0,0)
        // axial modes. The sufficiency-checked box must return the
        // pure axial ladder omega_l = l pi c0 / lx.
        let cav = rectangular_cavity_modes(
            1.0,
            0.1,
            0.1,
            AcousticMedium::air(),
            0.0,
            6,
            &[[0.25, 0.05]],
        )
        .expect("elongated cavity");
        let pi = core::f64::consts::PI;
        for (l, &omega) in cav.omegas.iter().enumerate() {
            let expected = l as f64 * pi * 343.0 / 1.0;
            assert!(
                (omega - expected).abs() <= 1e-9 * expected.max(1.0),
                "mode {l}: {omega} vs axial ladder {expected}"
            );
        }
        println!(
            "{{\"suite\":\"fs-couple-vibro\",\"case\":\"elongated-cavity-regression\",\"verdict\":\"pass\"}}"
        );
    }

    #[test]
    fn multi_mode_pencil_roots_zero_the_independent_determinant() {
        // COVERAGE (fresh-eyes review): every prior pencil test was
        // 1x1, where the C / C^T index arithmetic collapses. Here a
        // 2x2 model with four DISTINCT couplings is checked against an
        // independently assembled determinant of the ORIGINAL
        // (non-linearized) coupled block matrix
        //   D(x) = [[Omega_s^2 - x I, -C], [-s_q x C^T, Omega_a^2 - x I]]
        // — a transpose-index typo in the engine's linearization gives
        // roots that fail to zero this determinant.
        let structure = StructuralModes {
            omegas: vec![100.0, 170.0],
            shapes: vec![vec![1.0], vec![1.0]],
            loss_factor: 0.0,
        };
        let cavity = CavityModes {
            omegas: vec![140.0, 260.0],
            lambdas: vec![0.02, 0.011],
            interface: vec![vec![1.0], vec![1.0]],
            loss_factor: 0.0,
            rho0: 1.204,
            c0: 343.0,
        };
        // Row-major n_s x n_a with four distinct entries.
        let coupling = vec![2.0e-3, -1.1e-3, 0.7e-3, 1.9e-3];
        let model = VibroacousticModel::try_new(&structure, &cavity, coupling.clone(), None)
            .expect("model");
        let freqs = model.undamped_natural_frequencies().expect("pencil");
        assert_eq!(freqs.len(), 4);
        for pair in freqs.windows(2) {
            assert!(pair[0] <= pair[1], "roots must be ascending");
        }
        let rho_c2 = 1.204 * 343.0 * 343.0;
        let det_at = |x: f64| -> f64 {
            let n = 4usize;
            let mut d = vec![C64::ZERO; n * n];
            for r in 0..2 {
                d[r * n + r] = C64::from_re(structure.omegas[r] * structure.omegas[r] - x);
                for q in 0..2 {
                    d[r * n + (2 + q)] = C64::from_re(-coupling[r * 2 + q]);
                }
            }
            for q in 0..2 {
                let s_q = rho_c2 / cavity.lambdas[q];
                d[(2 + q) * n + (2 + q)] = C64::from_re(cavity.omegas[q] * cavity.omegas[q] - x);
                for r in 0..2 {
                    d[(2 + q) * n + r] = C64::from_re(-s_q * x * coupling[r * 2 + q]);
                }
            }
            det_complex(&d, n).abs()
        };
        // Scale: determinant magnitude away from any root.
        let scale =
            det_at(f64::midpoint(freqs[0] * freqs[0], freqs[1] * freqs[1])).max(det_at(1.0));
        for &f in &freqs {
            let residual = det_at(f * f) / scale;
            assert!(
                residual < 1e-6,
                "pencil root {f} must zero the independent determinant: {residual:.3e}"
            );
        }
        println!(
            "{{\"suite\":\"fs-couple-vibro\",\"case\":\"multi-mode-pencil\",\"roots\":{freqs:?},\"verdict\":\"pass\"}}"
        );
    }
}
