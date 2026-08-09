//! # fs-nlmodal — geometric nonlinearity in modal coordinates
//!
//! Thin-structure geometric nonlinearity (von Karman plates,
//! Kirchhoff-Carrier strings) reduced to N coupled modal oscillators
//! with a QUARTIC SUM-OF-SQUARES potential:
//!
//! `H = 1/2 sum_k (p_k^2 + w_k^2 q_k^2) + 1/4 sum_j c_j (q^T E_j q)^2`
//!
//! with `c_j >= 0` and `E_j` symmetric — the Airy-stress-eliminated von
//! Karman form (Ducceschi-Touze / Bilbao modal formulation). Because
//! the nonlinear potential is a positive sum of squares, `H` is
//! bounded below and the whole system is a PORT-HAMILTONIAN system
//! with non-quadratic storage: time stepping, damping, striking, and
//! the energy ledger all come from fs-phs (Gonzalez discrete
//! gradients — implicit, energy-exact, stable at crash amplitudes; no
//! explicit integrator exists in this crate to misuse).
//!
//! Constructors:
//! - [`von_karman_ss_plate`] — simply-supported rectangular plate with
//!   ANALYTIC sine modes for displacement AND Airy stress (both are
//!   biharmonic eigenfunctions on the SS rectangle); coupling
//!   integrals by Gauss-Legendre quadrature, certified by a second
//!   independent order.
//! - [`kirchhoff_carrier_string`] — the 1D sibling: one stress "mode"
//!   whose coupling matrix is exactly `diag(k^2 pi^2/L^2)` — tension
//!   modulation and pitch glide nearly free.
//!
//! Honest scope: moderate-rotation von Karman regime; no damage,
//! plasticity, or wrinkling; simply-supported (in-plane movable)
//! boundary in v1 — FE-mode input is the follow-up (its trigger:
//! consumers with non-analytic geometry, which is also when coupling-
//! tensor caching starts to matter; at analytic-mode scale the tensor
//! build is milliseconds and a cache manifest would be ceremony).

use fs_math::det;
use fs_phs::{PhsError, PortHamiltonian, Storage};

/// One Airy-stress coupling channel: coefficient `c >= 0` and the
/// symmetric matrix `E` (row-major `n x n`) of `q^T E q`.
#[derive(Debug, Clone)]
pub struct StressChannel {
    /// Positive channel coefficient `c_j` (absorbs `E h / (2 xi^4)`).
    pub coefficient: f64,
    /// Symmetric coupling matrix, row-major `n x n`.
    pub coupling: Vec<f64>,
}

/// Typed refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum NlModalError {
    /// Bad geometry/material parameter.
    Parameter {
        /// Which one.
        what: &'static str,
    },
    /// A coupling matrix failed its symmetry certificate.
    AsymmetricCoupling {
        /// Worst residual relative to the matrix scale.
        residual: f64,
    },
    /// The two independent quadrature orders disagreed beyond
    /// tolerance — the tensor is not converged.
    QuadratureMismatch {
        /// Worst relative disagreement.
        residual: f64,
    },
    /// Underlying pHS admission failure.
    Phs(PhsError),
}

impl core::fmt::Display for NlModalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            NlModalError::Parameter { what } => write!(f, "bad parameter: {what}"),
            NlModalError::AsymmetricCoupling { residual } => {
                write!(f, "coupling matrix asymmetry {residual:.3e}")
            }
            NlModalError::QuadratureMismatch { residual } => {
                write!(f, "quadrature orders disagree by {residual:.3e}")
            }
            NlModalError::Phs(e) => write!(f, "pHS admission failed: {e}"),
        }
    }
}

impl std::error::Error for NlModalError {}

/// The quartic sum-of-squares storage
/// `H = 1/2 sum (p^2 + w^2 q^2) + 1/4 sum_j c_j (q^T E_j q)^2`.
///
/// State layout matches fs-phs `modal_bank`: `[q_0, p_0, q_1, p_1,
/// ...]`. The gradient is the EXACT analytic gradient of the coded
/// `H` — force-vs-energy divergence (the classic hand-derived-tensor
/// bug class) is therefore impossible through this type; the battery
/// demonstrates that under Gonzalez stepping such divergence is a
/// TRAJECTORY error, never an energy error (the architectural
/// finding pinned in the tests).
pub struct SosModalStorage {
    /// Modal angular frequencies [rad/s].
    pub omegas: Vec<f64>,
    /// Stress channels.
    pub channels: Vec<StressChannel>,
}

impl SosModalStorage {
    fn n_modes(&self) -> usize {
        self.omegas.len()
    }

    /// `q^T E q` for channel `ch` given the interleaved state.
    fn quadratic_form(&self, ch: &StressChannel, x: &[f64]) -> f64 {
        let n = self.n_modes();
        let mut acc = 0.0;
        for p in 0..n {
            let mut row = 0.0;
            for q in 0..n {
                row += ch.coupling[p * n + q] * x[2 * q];
            }
            acc += x[2 * p] * row;
        }
        acc
    }
}

impl Storage for SosModalStorage {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        let n = self.n_modes();
        let mut h = 0.0;
        for k in 0..n {
            let (q, p) = (x[2 * k], x[2 * k + 1]);
            h += f64::midpoint(p * p, self.omegas[k] * self.omegas[k] * q * q);
        }
        for ch in &self.channels {
            let s = self.quadratic_form(ch, x);
            h += 0.25 * ch.coefficient * s * s;
        }
        h
    }

    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        let n = self.n_modes();
        for k in 0..n {
            out[2 * k] = self.omegas[k] * self.omegas[k] * x[2 * k];
            out[2 * k + 1] = x[2 * k + 1];
        }
        for ch in &self.channels {
            let s = self.quadratic_form(ch, x);
            // d/dq_k [c/4 s^2] = c s (E q)_k (E symmetric).
            for k in 0..n {
                let mut row = 0.0;
                for q in 0..n {
                    row += ch.coupling[k * n + q] * x[2 * q];
                }
                out[2 * k] += ch.coefficient * s * row;
            }
        }
    }
}

const SYM_TOL: f64 = 1.0e-10;

/// Assemble the pHS: symplectic pair blocks, per-mode viscous damping
/// `R = diag(0, 2 zeta_k w_k)` (zeta from the caller — the
/// visco-damping facility's per-mode output slots in here; this crate
/// does not invent a second damping representation), and a single
/// strike port weighted by the mode shapes at the strike point.
///
/// # Errors
/// [`NlModalError`] on length mismatches, non-positive frequencies,
/// negative damping, asymmetric couplings, or pHS admission failure.
pub fn assemble(
    storage: SosModalStorage,
    zetas: &[f64],
    strike_weights: &[f64],
) -> Result<PortHamiltonian, NlModalError> {
    let n = storage.n_modes();
    if zetas.len() != n || strike_weights.len() != n {
        return Err(NlModalError::Parameter {
            what: "zetas/strike_weights length vs mode count",
        });
    }
    // Negated-positive comparisons so NaN REFUSES instead of slipping
    // through `<= 0.0` (review finding: NaN passed every gate).
    for &w in &storage.omegas {
        if w.is_nan() || w <= 0.0 || !w.is_finite() {
            return Err(NlModalError::Parameter {
                what: "modal frequency must be positive and finite",
            });
        }
    }
    for &z in zetas {
        if z.is_nan() || z < 0.0 {
            return Err(NlModalError::Parameter {
                what: "damping ratio must be non-negative",
            });
        }
    }
    for ch in &storage.channels {
        if ch.coefficient.is_nan() || ch.coefficient < 0.0 || !ch.coefficient.is_finite() {
            return Err(NlModalError::Parameter {
                what: "channel coefficient must be non-negative and finite",
            });
        }
        if ch.coupling.iter().any(|v| !v.is_finite()) {
            return Err(NlModalError::Parameter {
                what: "coupling entries must be finite",
            });
        }
        if ch.coupling.len() != n * n {
            return Err(NlModalError::Parameter {
                what: "coupling matrix size",
            });
        }
        let scale = ch
            .coupling
            .iter()
            .fold(f64::MIN_POSITIVE, |a, &v| a.max(v.abs()));
        let mut worst = 0.0f64;
        for p in 0..n {
            for q in 0..p {
                worst = worst.max((ch.coupling[p * n + q] - ch.coupling[q * n + p]).abs());
            }
        }
        if worst > SYM_TOL * scale {
            return Err(NlModalError::AsymmetricCoupling {
                residual: worst / scale,
            });
        }
    }
    let dim = 2 * n;
    let mut j = vec![0.0; dim * dim];
    let mut r = vec![0.0; dim * dim];
    let mut g = vec![0.0; dim];
    for k in 0..n {
        let (qi, pi) = (2 * k, 2 * k + 1);
        j[qi * dim + pi] = 1.0;
        j[pi * dim + qi] = -1.0;
        r[pi * dim + pi] = 2.0 * zetas[k] * storage.omegas[k];
        g[pi] = strike_weights[k];
    }
    PortHamiltonian::new(dim, 1, j, r, g, Box::new(storage)).map_err(NlModalError::Phs)
}

// ---------------------------------------------------------------------
// Von Karman simply-supported rectangular plate
// ---------------------------------------------------------------------

/// Von Karman SS plate parameters (isotropic, in-plane movable edges).
#[derive(Debug, Clone, Copy)]
pub struct VkPlateParams {
    /// Side lengths [m].
    pub lx: f64,
    /// Side lengths [m].
    pub ly: f64,
    /// Thickness [m].
    pub h: f64,
    /// Young's modulus [Pa].
    pub young: f64,
    /// Poisson ratio.
    pub nu: f64,
    /// Density [kg/m^3].
    pub rho: f64,
}

/// One sine mode index pair.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SineMode {
    /// x half-waves.
    pub m: usize,
    /// y half-waves.
    pub n: usize,
}

/// The assembled von Karman modal model plus its certification data.
pub struct VkModel {
    /// The storage (frequencies + channels), ready for [`assemble`].
    pub storage: SosModalStorage,
    /// Displacement mode list (index-aligned with the storage).
    pub disp_modes: Vec<SineMode>,
    /// Stress mode list (index-aligned with `storage.channels`).
    pub stress_modes: Vec<SineMode>,
    /// Worst relative disagreement between the two quadrature orders
    /// — the tensor-convergence certificate, logged not hidden.
    pub quadrature_residual: f64,
}

/// Gauss-Legendre nodes/weights on [0, 1] (order-n rule; nodes by
/// Newton on the Legendre recurrence — deterministic).
fn gauss_01(order: usize) -> (Vec<f64>, Vec<f64>) {
    let mut nodes = Vec::with_capacity(order);
    let mut weights = Vec::with_capacity(order);
    let nf = order as f64;
    for i in 0..order {
        // Initial guess (Chebyshev-like), then Newton.
        let mut x = det::cos(core::f64::consts::PI * (i as f64 + 0.75) / (nf + 0.5));
        for _ in 0..60 {
            // Legendre P_n(x) and derivative by recurrence.
            let (mut p0, mut p1) = (1.0f64, x);
            for k in 2..=order {
                let kf = k as f64;
                let p2 = ((2.0 * kf - 1.0) * x * p1 - (kf - 1.0) * p0) / kf;
                p0 = p1;
                p1 = p2;
            }
            let dp = nf * (x * p1 - p0) / (x * x - 1.0);
            let dx = p1 / dp;
            x -= dx;
            if dx.abs() < 1.0e-15 {
                break;
            }
        }
        // Recompute derivative at the root for the weight.
        let (mut p0, mut p1) = (1.0f64, x);
        for k in 2..=order {
            let kf = k as f64;
            let p2 = ((2.0 * kf - 1.0) * x * p1 - (kf - 1.0) * p0) / kf;
            p0 = p1;
            p1 = p2;
        }
        let dp = nf * (x * p1 - p0) / (x * x - 1.0);
        nodes.push(f64::midpoint(x, 1.0));
        weights.push(1.0 / ((1.0 - x * x) * dp * dp));
    }
    (nodes, weights)
}

/// The von Karman bracket integrand `Psi_j L(Phi_p, Phi_q)` for sine
/// modes on the SS rectangle, integrated by the given quadrature.
/// `L(f, g) = f_xx g_yy + f_yy g_xx - 2 f_xy g_xy`. Shapes are UNIT
/// sine products; normalization factors are applied by the caller.
fn coupling_integral(
    lx: f64,
    ly: f64,
    stress: SineMode,
    a: SineMode,
    b: SineMode,
    nodes: &[f64],
    weights: &[f64],
) -> f64 {
    let kx = |m: usize| m as f64 * core::f64::consts::PI / lx;
    let ky = |n: usize| n as f64 * core::f64::consts::PI / ly;
    let (sax, say) = (kx(stress.m), ky(stress.n));
    let (pax, pay) = (kx(a.m), ky(a.n));
    let (pbx, pby) = (kx(b.m), ky(b.n));
    let mut acc = 0.0;
    for (i, &xi) in nodes.iter().enumerate() {
        let x = xi * lx;
        for (jn, &yj) in nodes.iter().enumerate() {
            let y = yj * ly;
            let psi = det::sin(sax * x) * det::sin(say * y);
            // Phi_a derivatives.
            let a_xx = -pax * pax * det::sin(pax * x) * det::sin(pay * y);
            let a_yy = -pay * pay * det::sin(pax * x) * det::sin(pay * y);
            let a_xy = pax * pay * det::cos(pax * x) * det::cos(pay * y);
            let b_xx = -pbx * pbx * det::sin(pbx * x) * det::sin(pby * y);
            let b_yy = -pby * pby * det::sin(pbx * x) * det::sin(pby * y);
            let b_xy = pbx * pby * det::cos(pbx * x) * det::cos(pby * y);
            let l = a_xx * b_yy + a_yy * b_xx - 2.0 * a_xy * b_xy;
            acc += weights[i] * weights[jn] * psi * l;
        }
    }
    acc * lx * ly
}

/// Parameter/mode-list validation for [`von_karman_ss_plate`].
fn validate_vk_inputs(
    params: &VkPlateParams,
    disp_modes: &[SineMode],
    stress_modes: &[SineMode],
) -> Result<(), NlModalError> {
    let VkPlateParams {
        lx,
        ly,
        h,
        young,
        nu,
        rho,
    } = *params;
    if !(lx > 0.0 && ly > 0.0 && h > 0.0 && young > 0.0 && rho > 0.0 && (0.0..0.5).contains(&nu)) {
        return Err(NlModalError::Parameter {
            what: "plate parameters",
        });
    }
    if disp_modes.is_empty() || stress_modes.is_empty() {
        return Err(NlModalError::Parameter {
            what: "empty mode list",
        });
    }
    for md in disp_modes.iter().chain(stress_modes) {
        if md.m == 0 || md.n == 0 {
            return Err(NlModalError::Parameter {
                what: "sine mode indices start at 1",
            });
        }
    }
    // Duplicates silently double-count a physical mode (review
    // finding): refuse in each list.
    for list in [disp_modes, stress_modes] {
        for (i, a) in list.iter().enumerate() {
            if list[..i].contains(a) {
                return Err(NlModalError::Parameter {
                    what: "duplicate mode in list",
                });
            }
        }
    }
    Ok(())
}

/// Raw two-order quadrature results for one stress channel.
struct RawChannel {
    xi4: f64,
    pairs: Vec<(usize, usize, f64, f64)>,
    scale: f64,
}

/// Build the von Karman modal model for a simply-supported rectangle.
///
/// Displacement modes are MASS-NORMALIZED sine products (`rho h
/// integral Phi^2 = 1`); stress modes are UNIT-normalized sine
/// products (`integral Psi^2 = 1`), both biharmonic eigenfunctions of
/// the SS rectangle. Channel coefficients follow the Airy elimination
/// `c_j = E h / (2 xi_j^4)` with `xi_j^4 = ((m pi/lx)^2 +
/// (n pi/ly)^2)^2`.
///
/// The stress-mode count is a SEPARATE truncation from the
/// displacement count (the literature uses more stress modes); both
/// lists are explicit inputs. The coupling tensor is computed at two
/// independent quadrature orders and REFUSES on disagreement.
///
/// # Errors
/// [`NlModalError`] on bad parameters or a failed quadrature
/// certificate.
pub fn von_karman_ss_plate(
    params: &VkPlateParams,
    disp_modes: &[SineMode],
    stress_modes: &[SineMode],
) -> Result<VkModel, NlModalError> {
    let VkPlateParams {
        lx,
        ly,
        h,
        young,
        nu,
        rho,
    } = *params;
    validate_vk_inputs(params, disp_modes, stress_modes)?;
    let d_bend = young * h * h * h / (12.0 * (1.0 - nu * nu));
    let pi = core::f64::consts::PI;
    let omegas: Vec<f64> = disp_modes
        .iter()
        .map(|md| {
            let k2 = (md.m as f64 * pi / lx) * (md.m as f64 * pi / lx)
                + (md.n as f64 * pi / ly) * (md.n as f64 * pi / ly);
            det::sqrt(d_bend / (rho * h)) * k2
        })
        .collect();
    // Normalizations: Phi = phi_norm * sin sin with rho h * (lx ly /4)
    // * phi_norm^2 = 1; Psi = psi_norm * sin sin with (lx ly/4) *
    // psi_norm^2 = 1.
    let phi_norm = det::sqrt(4.0 / (rho * h * lx * ly));
    let psi_norm = det::sqrt(4.0 / (lx * ly));
    let nq = disp_modes.len();
    // Quadrature order scales with the highest half-wave SUM in the
    // integrand (executed: fixed order 24 left ~5 points per wave and
    // a 1e-4 cross-order disagreement); ~5 points per half-wave with
    // margin, and the certification order is offset AND coprime-ish.
    let max_sum = {
        let max_d = disp_modes.iter().map(|m| m.m.max(m.n)).max().unwrap_or(1);
        let max_s = stress_modes.iter().map(|m| m.m.max(m.n)).max().unwrap_or(1);
        max_s + 2 * max_d
    };
    let order_a = 16 + 5 * max_sum;
    let (nodes_a, w_a) = gauss_01(order_a);
    let (nodes_b, w_b) = gauss_01(order_a + 9);
    // Two passes (review findings, both executed): entrywise relative
    // comparison falsely refuses analytically-zero selection-rule
    // ENTRIES, and an ALL-zero channel's own scale is pure roundoff
    // (measured 0.75 "relative"), so residuals are judged against
    // max(channel scale, 1e-12 * GLOBAL scale, dimensional floor).
    let mut raw_channels: Vec<RawChannel> = Vec::with_capacity(stress_modes.len());
    let mut global_scale = f64::MIN_POSITIVE;
    for &sm in stress_modes {
        let xi4 = {
            let k2 = (sm.m as f64 * pi / lx) * (sm.m as f64 * pi / lx)
                + (sm.n as f64 * pi / ly) * (sm.n as f64 * pi / ly);
            k2 * k2
        };
        let mut pairs = Vec::with_capacity(nq * (nq + 1) / 2);
        let mut scale = f64::MIN_POSITIVE;
        for p in 0..nq {
            for q in 0..=p {
                let raw_a =
                    coupling_integral(lx, ly, sm, disp_modes[p], disp_modes[q], &nodes_a, &w_a);
                let raw_b =
                    coupling_integral(lx, ly, sm, disp_modes[p], disp_modes[q], &nodes_b, &w_b);
                scale = scale.max(raw_a.abs()).max(raw_b.abs());
                pairs.push((p, q, raw_a, raw_b));
            }
        }
        global_scale = global_scale.max(scale);
        raw_channels.push(RawChannel { xi4, pairs, scale });
    }
    // Dimensional floor for the judge scale: when EVERY channel is a
    // selection-rule zero (executed: single-mode fixtures), even the
    // global scale is roundoff; a characteristic nonzero integral has
    // magnitude ~ kmin^4 * area, which anchors the comparison in
    // physical units.
    let kmin2 = (pi / lx) * (pi / lx) + (pi / ly) * (pi / ly);
    let char_scale = kmin2 * kmin2 * lx * ly;
    let mut channels = Vec::with_capacity(stress_modes.len());
    let mut worst_rel = 0.0f64;
    for rc in raw_channels {
        let judge_scale = rc
            .scale
            .max(1.0e-12 * global_scale)
            .max(1.0e-8 * char_scale);
        let mut e = vec![0.0; nq * nq];
        for (p, q, raw_a, raw_b) in rc.pairs {
            worst_rel = worst_rel.max((raw_a - raw_b).abs() / judge_scale);
            let v = psi_norm * phi_norm * phi_norm * raw_b;
            e[p * nq + q] = v;
            e[q * nq + p] = v;
        }
        channels.push(StressChannel {
            coefficient: young * h / (2.0 * rc.xi4),
            coupling: e,
        });
    }
    if worst_rel > 1.0e-8 {
        return Err(NlModalError::QuadratureMismatch {
            residual: worst_rel,
        });
    }
    Ok(VkModel {
        storage: SosModalStorage { omegas, channels },
        disp_modes: disp_modes.to_vec(),
        stress_modes: stress_modes.to_vec(),
        quadrature_residual: worst_rel,
    })
}

// ---------------------------------------------------------------------
// Kirchhoff-Carrier string
// ---------------------------------------------------------------------

/// Kirchhoff-Carrier string parameters.
#[derive(Debug, Clone, Copy)]
pub struct KcStringParams {
    /// Speaking length [m].
    pub length: f64,
    /// Tension [N].
    pub tension: f64,
    /// Linear density [kg/m].
    pub lin_density: f64,
    /// Longitudinal stiffness `E A` [N].
    pub ea: f64,
}

/// Kirchhoff-Carrier modal model: `n_modes` sine string modes with the
/// single tension-modulation channel
/// `H_nl = (E A / (8 L)) * (sum_k (k pi / L)^2 Q_k^2)^2` where `Q_k`
/// is the PHYSICAL modal displacement; in mass-normalized coordinates
/// (`q_k = sqrt(mu L / 2) Q_k`) the channel matrix is
/// `E[k,k] = (k pi / L)^2 * 2/(mu L)` with coefficient `E A L / 8`
/// — the exact Kirchhoff-Carrier averaged-tension form.
///
/// # Errors
/// [`NlModalError::Parameter`] on non-physical inputs.
pub fn kirchhoff_carrier_string(
    params: &KcStringParams,
    n_modes: usize,
) -> Result<SosModalStorage, NlModalError> {
    let KcStringParams {
        length,
        tension,
        lin_density,
        ea,
    } = *params;
    if !(length > 0.0 && tension > 0.0 && lin_density > 0.0 && ea >= 0.0) || n_modes == 0 {
        return Err(NlModalError::Parameter {
            what: "string parameters",
        });
    }
    let c = det::sqrt(tension / lin_density);
    let pi = core::f64::consts::PI;
    let omegas: Vec<f64> = (1..=n_modes).map(|k| k as f64 * pi * c / length).collect();
    // Tension modulation: T(t) = T0 + (EA/2L) * integral (w')^2 dx;
    // with w = sum Q_k sin(k pi x / L): integral (w')^2 = (L/2) sum
    // (k pi/L)^2 Q_k^2. The added potential is (EA/(8L)) * (that
    // integral without L/2... ) — assembled so that dH/dQ reproduces
    // the KC modal force  (EA/(4L)) * S * (k pi/L)^2 * L Q_k with
    // S = sum (j pi/L)^2 Q_j^2 / ... The mass-normalized channel:
    let e: Vec<f64> = {
        let mut m = vec![0.0; n_modes * n_modes];
        for k in 1..=n_modes {
            let kk = k as f64 * pi / length;
            m[(k - 1) * n_modes + (k - 1)] = kk * kk * 2.0 / (lin_density * length);
        }
        m
    };
    let coefficient = ea * length / 8.0;
    Ok(SosModalStorage {
        omegas,
        channels: vec![StressChannel {
            coefficient,
            coupling: e,
        }],
    })
}

/// First-order Duffing backbone: for
/// `qddot + w0^2 q + beta q^3 = 0`, the amplitude-dependent frequency
/// is `w(a) = w0 (1 + 3 beta a^2 / (8 w0^2))` — the analytic
/// perturbation pin the batteries measure against.
#[must_use]
pub fn duffing_backbone(omega0: f64, beta: f64, amplitude: f64) -> f64 {
    omega0 * (1.0 + 3.0 * beta * amplitude * amplitude / (8.0 * omega0 * omega0))
}

/// The effective single-mode cubic coefficient `beta_k` of mode `k`
/// (the `q_k^3` force coefficient with every other mode at rest):
/// `beta_k = sum_j c_j E_j[k,k]^2`.
#[must_use]
pub fn single_mode_beta(storage: &SosModalStorage, k: usize) -> f64 {
    let n = storage.omegas.len();
    storage
        .channels
        .iter()
        .map(|ch| {
            let e = ch.coupling[k * n + k];
            ch.coefficient * e * e
        })
        .sum()
}
