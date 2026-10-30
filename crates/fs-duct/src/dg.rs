//! Nodal discontinuous-Galerkin (DG) time-domain discretization of the lossless
//! variable-area plane-wave system (IRP M2.2-W3a, bead `frankensim-irp-v5-instrument-realism-pexfp.83`).
//!
//! # Governing equations
//!
//! On a 1-D axial domain $x \in [0, L]$ with positive cross-sectional area $A(x) = \pi r(x)^2$,
//! prescribed equilibrium gas density $\rho(x)$, and adiabatic sound speed $c(x)$, acoustic pressure $p(x, t)$ [Pa]
//! and volume velocity $q(x, t)$ [$\text{m}^3/\text{s}$] satisfy the lossless Webster conservation
//! system:
//!
//! $$\frac{A(x)}{\rho c^2} \partial_t p + \partial_x q = 0, \qquad \frac{\rho}{A(x)} \partial_t q + \partial_x p = 0.$$
//!
//! The stored acoustic energy is:
//!
//! $$E(t) = \frac{1}{2} \int_0^L \left( \frac{A(x)}{\rho c^2} p(x, t)^2 + \frac{\rho}{A(x)} q(x, t)^2 \right) dx.$$
//!
//! # Spatial discretization (LGL nodal DG with alternating energy-conserving flux)
//!
//! Each element $e \in \{0, \dots, K-1\}$ of length $h_e = x_R^{(e)} - x_L^{(e)}$ carries
//! Legendre–Gauss–Lobatto (LGL) nodal polynomials of degree $p \in [1, \text{MAX\_DG\_ORDER}]$.
//! Element mass matrices $M_p^{(e)}$ and $M_q^{(e)}$ integrate the exact segment area
//! (`Segment::Cylinder`, `Segment::Cone`, `Segment::Profile`) by Gauss–Legendre rules doubled
//! from $2p + 4$ points until they repeat to roundoff, so the rational weight $\rho / A(x)$ of a
//! taper is integrated to roundoff, not to a fixed rule.
//!
//! At every interior element interface $e + 1/2$, the numerical flux uses the energy-conserving
//! alternating trace (Xing, Chou & Shu 2013; Chou, Shu & Xing 2014):
//!
//! $$\hat{q}_{e+1/2} = q_R^{(e)}, \qquad \hat{p}_{e+1/2} = p_L^{(e+1)},$$
//!
//! which satisfies $K_q = -K_p^\top$ identically in the discrete energy inner product while
//! delivering optimal $L^2$ convergence order $p + 1$ for both even and odd polynomial degrees $p$.
//!
//! # Exact exponential time propagation (Plan §3.3)
//!
//! In energy-weighted coordinates $\tilde{p} = M_p^{1/2} p$, $\tilde{q} = M_q^{1/2} q$, the
//! semi-discrete operator with lossless boundaries is skew,
//!
//! $$\tilde{L} = \begin{bmatrix} 0 & -G \\ G^\top & 0 \end{bmatrix}, \qquad G = M_p^{-1/2} K_p M_q^{-1/2},$$
//!
//! so $E = \tfrac12 \|\tilde x\|^2$ is conserved exactly. [`DgDuct::propagator`] diagonalizes it
//! into $2 \times 2$ rotations through the SVD of $G$, taken from the symmetric Jordan–Wielandt
//! matrix, and admits it only under a certificate (singular residual and basis orthogonality) at
//! most `tol`. A `Matched` boundary adds a dissipative block; that operator is exponentiated by a
//! scaled Taylor series with squaring, certified by a rigorous truncation bound and an
//! a-posteriori halving estimate. No Runge–Kutta, trapezoid, implicit-midpoint or
//! generalized-alpha step is used. There is no Krylov path: propagators above
//! [`MAX_DG_PROPAGATOR_DOFS`] refuse.
//!
//! # Energy ledger and ports
//!
//! A port holds its effort (inflow volume velocity) over a step and returns the exact impulse
//! $\Delta X = \int p\,dt$ of its conjugate pressure, so the exchanged work is exactly
//! $e\,\Delta X$. Every ledger row reports the stored energy of the reconstructed physical
//! state. For lossless boundaries the balance $E(t) - E(0) - W$ is a conservation check; for a
//! `Matched` boundary the absorbed energy is inferred from the deficit, so that balance is not a
//! check (`balance_checked` is false).
//!
//! [`DgDuct::with_side_ports`] exposes existing interior faces as point flow
//! junctions, using the same alternating pressure trace and exact work integral.
//! Complete side loads interact through [`DgPropagator::multiport_response`];
//! finite chimney geometry and aperture corrections belong to those loads.
//!
//! [`DgDuct::with_viscothermal`] optionally appends the exact-Bessel-derived
//! positive viscous and thermal memory supplied by `fs_phs::ZkDiffusive`.
//! Identical circular cylinders in a uniform material gas are supported. These
//! states share the same dense exponential and energy coordinates; no midpoint
//! relaxation or duplicated inviscid storage is introduced. Total energy then
//! includes wall memory, and the dense ledger still reports inferred loss.
//!
//! [`DgDuct::with_compliant_wall`] adds the displacement and velocity of a
//! locally reacting cylindrical wall to that SAME operator. Its pressure work
//! is reciprocal and its elastic/kinetic storage is part of the native state.

use crate::profile::gas_field::GasFieldDuctError;
use crate::{Duct, DuctError, Segment};
use fs_blake3::{ContentHash, DomainHasher};
use fs_la::eigen::{admit_jacobi_eigh, householder_eigh, jacobi_eigh};
use fs_material::gas::GasState;
use fs_material::gas_field::AxialGasField;
use fs_math::det;

pub mod continuous;
pub use continuous::{DgContinuous, DgHarmonicPort, DgHarmonicResponse};

mod viscothermal;
pub use viscothermal::{DgWallAdmission, DgWallResolution};
mod compliant_wall;
pub use compliant_wall::{DgCompliantWallState, DgWallMaterialUpdate, DgWallRelaxationTransfer};
mod port_response;
pub use port_response::{DgMultiportResponse, DgPortResponse};
mod side_ports;
pub use side_ports::{DgSidePort, DgSidePortSource, MAX_DG_SIDE_PORTS};
mod gas_update;
pub use gas_update::DgGasMaterialUpdate;
mod molecular;
pub use molecular::heat::DgMolecularHeat;

/// Maximum admitted nodal polynomial order $p$ per element.
pub const MAX_DG_ORDER: usize = 8;

/// Maximum admitted element count $K$ in a single 1-D DG duct.
pub const MAX_DG_ELEMENTS: usize = 512;

/// Largest total DOF count $N = 2 K (p + 1)$ a dense propagator admits. There is no Krylov
/// exp-action path, so larger discretizations refuse by name.
pub const MAX_DG_PROPAGATOR_DOFS: usize = 2048;

/// Jordan–Wielandt eigenvalues within this fraction of the spectral radius are null
/// singular values of `G` (the closed–closed duct's uniform-pressure mode).
const NULL_SINGULAR_REL: f64 = 1e-9;

/// Mass-matrix quadrature stops doubling once both matrices repeat to this many ulps of their
/// largest entry.
const MASS_QUADRATURE_ULPS: f64 = 64.0;

/// Gauss–Legendre points beyond which the mass quadrature refuses as unconverged.
const MAX_MASS_QUADRATURE_POINTS: usize = 512;

/// Stable identifier for a boundary port coupled to an external valve, reed, lip, or radiation load.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct PortId(pub u32);

/// Refinable spatial discretization axes $(h, p)$ for a 1-D DG duct.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DgMesh {
    /// Total number of axial elements $K \in [1, \text{MAX\_DG\_ELEMENTS}]$.
    pub elements: usize,
    /// Polynomial order $p \in [1, \text{MAX\_DG\_ORDER}]$ on each element ($p + 1$ LGL nodes).
    pub order: usize,
}

/// Boundary condition at the inlet ($x = 0$) or outlet ($x = L$) of a [`DgDuct`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DgBoundary {
    /// Rigid wall: zero volume flow ($\hat{q} = 0$, $\hat{p} = p_{\text{int}}$).
    Closed,
    /// Ideal pressure-release termination ($\hat{p} = 0$, $\hat{q} = q_{\text{int}}$).
    IdealOpen,
    /// Characteristic non-reflecting termination with local impedance $Z_c = \rho c / A$.
    Matched,
    /// Typed external port carrying held volume-flow effort $q_{\text{in}}$ [$\text{m}^3/\text{s}$]
    /// (into the domain) and returning conjugate boundary pressure $p$ [Pa] and integrated impulse
    /// $\Delta X = \int_0^{\Delta t} p(\tau)\,d\tau$ [$\text{Pa}\cdot\text{s}$] so port work is
    /// exact ($W = q_{\text{in}} \cdot \Delta X$).
    Port(PortId),
}

/// Held effort supplied to a [`DgBoundary::Port`] over one time step $\Delta t$.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PortInput {
    /// Target port identifier.
    pub id: PortId,
    /// Held inflow volume velocity $q_{\text{in}}$ [$\text{m}^3/\text{s}$] over $[t_n, t_n + \Delta t]$.
    pub effort: f64,
}

/// Conjugate response emitted by a [`DgBoundary::Port`] over one time step $\Delta t$.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PortOutput {
    /// Port identifier.
    pub id: PortId,
    /// Step-averaged conjugate pressure $\bar{p} = \Delta x / \Delta t$ [Pa].
    pub conjugate_rate: f64,
    /// ZOH-integrated conjugate impulse $\Delta x = \int_0^{\Delta t} p(\tau)\,d\tau$ [$\text{Pa}\cdot\text{s}$].
    /// Exchanged port work over the step is exactly `effort * delta_x`.
    pub delta_x: f64,
    /// Instantaneous end-of-step boundary pressure $p(t_{n+1})$ [Pa].
    pub end_conjugate: f64,
}

/// Per-step energy ledger entry emitted by [`DgPropagator::step`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EnergyLedgerRow {
    /// Completed step index (1-indexed after the first step).
    pub step: u64,
    /// Stored acoustic and any viscothermal memory energy $E(t_n)$ [J].
    pub stored_energy_j: f64,
    /// Exact port work [J], `sum(e * delta_x)` for held flow; a composed
    /// continuous trajectory instead supplies its actual `sum(integral p Q dt)`.
    pub port_work_step_j: f64,
    /// Cumulative work [J] injected by ports since initialization.
    pub cumulative_port_work_j: f64,
    /// Cumulative signed work [J] supplied by accepted material-parameter
    /// changes. This is separate from port work and irreversible heat.
    pub cumulative_parameter_work_j: f64,
    /// Cumulative energy [J] absorbed by `Matched` boundaries and any wall memory,
    /// The standalone path infers this from the energy deficit; a composed exact
    /// trajectory supplies independently integrated loss (zero when lossless).
    pub cumulative_dissipated_j: f64,
    /// Signed residual $E(t_n) - E(0) - W_{\text{port}} - W_{\text{parameter}} + D$ [J].
    pub energy_balance_residual_j: f64,
    /// Absolute residual divided by initial energy plus the magnitudes of
    /// accumulated port and parameter work (with a 1e-300 J floor).
    pub relative_energy_drift: f64,
    /// True for the lossless image or an independently integrated composed
    /// trajectory. False when boundary/wall loss is inferred from the deficit.
    pub balance_checked: bool,
}

/// Return value of a single [`DgPropagator::step`].
#[derive(Debug, Clone, PartialEq)]
pub struct DgStepOutput {
    /// Conjugate outputs for each active port on the duct (in inlet-then-outlet order).
    pub port_outputs: Vec<PortOutput>,
    /// Energy balance ledger row for the completed step.
    pub ledger: EnergyLedgerRow,
}

/// Legendre–Gauss–Lobatto quadrature rule and nodal differentiation matrix on $[-1, 1]$.
#[derive(Debug, Clone, PartialEq)]
pub struct LglRule {
    /// Polynomial degree $p$ ($p + 1$ nodes).
    pub order: usize,
    /// Ascending nodes $\xi_0 = -1 < \xi_1 < \dots < \xi_p = +1$.
    pub nodes: Vec<f64>,
    /// Positive quadrature weights $w_0, \dots, w_p$ summing to $2$.
    pub weights: Vec<f64>,
    /// Row-major $(p + 1) \times (p + 1)$ nodal differentiation matrix $D_{ij} = \ell_j'(\xi_i)$.
    pub diff_matrix: Vec<f64>,
    /// Barycentric interpolation weights $\beta_j = 1 / P_p(\xi_j)$.
    bary_weights: Vec<f64>,
}

/// Evaluate $P_n(\xi)$, $P'_n(\xi)$, and $P''_n(\xi)$ at $\xi \in [-1, 1]$ by deterministic
/// three-term recurrence.
fn legendre_p_and_derivs(n: usize, xi: f64) -> (f64, f64, f64) {
    if n == 0 {
        return (1.0, 0.0, 0.0);
    }
    if n == 1 {
        return (xi, 1.0, 0.0);
    }
    let mut p_prev = 1.0;
    let mut p_curr = xi;
    let mut dp_prev = 0.0;
    let mut dp_curr = 1.0;
    let mut d2p_prev = 0.0;
    let mut d2p_curr = 0.0;
    for k in 1..n {
        let kf = k as f64;
        let two_k_plus_1 = 2.0 * kf + 1.0;
        let k_plus_1 = kf + 1.0;
        let p_next = (two_k_plus_1 * xi * p_curr - kf * p_prev) / k_plus_1;
        let dp_next = (two_k_plus_1 * (p_curr + xi * dp_curr) - kf * dp_prev) / k_plus_1;
        let d2p_next = (two_k_plus_1 * (2.0 * dp_curr + xi * d2p_curr) - kf * d2p_prev) / k_plus_1;
        p_prev = p_curr;
        p_curr = p_next;
        dp_prev = dp_curr;
        dp_curr = dp_next;
        d2p_prev = d2p_curr;
        d2p_curr = d2p_next;
    }
    (p_curr, dp_curr, d2p_curr)
}

/// Construct the order-$p$ Legendre–Gauss–Lobatto rule on $[-1, 1]$.
///
/// # Errors
/// Returns [`DuctError::BadParameter`] if `order == 0` or `order > MAX_DG_ORDER`.
pub fn lgl_rule(order: usize) -> Result<LglRule, DuctError> {
    if order == 0 || order > MAX_DG_ORDER {
        return Err(DuctError::BadParameter {
            what: "DG polynomial order must be in 1..=MAX_DG_ORDER",
        });
    }
    let n_nodes = order + 1;
    let mut nodes = vec![0.0; n_nodes];
    nodes[0] = -1.0;
    nodes[order] = 1.0;

    let pf = order as f64;
    let half = order / 2;
    for j in 1..=half {
        if 2 * j == order {
            nodes[j] = 0.0;
            continue;
        }
        // Asymptotic interior root seed between Chebyshev-Lobatto knots.
        let theta = core::f64::consts::PI * (j as f64) / pf;
        let mut xi = -det::cos(theta);
        for _iter in 0..64 {
            let (_p, dp, d2p) = legendre_p_and_derivs(order, xi);
            let step = dp / d2p;
            xi -= step;
            if step.abs() <= 1e-16 * (1.0 + xi.abs()) {
                break;
            }
        }
        nodes[j] = xi;
        nodes[order - j] = -xi;
    }

    let mut p_vals = vec![0.0; n_nodes];
    let mut weights = vec![0.0; n_nodes];
    let mut bary_weights = vec![0.0; n_nodes];
    let denom_factor = pf * (pf + 1.0);
    for j in 0..n_nodes {
        let (pj, _, _) = legendre_p_and_derivs(order, nodes[j]);
        p_vals[j] = pj;
        weights[j] = 2.0 / (denom_factor * pj * pj);
        bary_weights[j] = 1.0 / pj;
    }

    let mut diff_matrix = vec![0.0; n_nodes * n_nodes];
    for i in 0..n_nodes {
        for j in 0..n_nodes {
            if i != j {
                diff_matrix[i * n_nodes + j] = p_vals[i] / (p_vals[j] * (nodes[i] - nodes[j]));
            }
        }
    }
    diff_matrix[0] = -0.25 * denom_factor;
    diff_matrix[order * n_nodes + order] = 0.25 * denom_factor;

    Ok(LglRule {
        order,
        nodes,
        weights,
        diff_matrix,
        bary_weights,
    })
}

impl LglRule {
    /// Evaluate the $p + 1$ Lagrange cardinal basis polynomials $\ell_0(\xi), \dots, \ell_p(\xi)$
    /// at $\xi \in [-1, 1]$.
    #[must_use]
    pub fn eval_basis(&self, xi: f64) -> Vec<f64> {
        let n = self.order + 1;
        for j in 0..n {
            if (xi - self.nodes[j]).abs() <= 1e-15 {
                let mut out = vec![0.0; n];
                out[j] = 1.0;
                return out;
            }
        }
        let mut numer = vec![0.0; n];
        let mut denom = 0.0;
        for j in 0..n {
            let t = self.bary_weights[j] / (xi - self.nodes[j]);
            numer[j] = t;
            denom += t;
        }
        for v in &mut numer {
            *v /= denom;
        }
        numer
    }
}

/// Compute an $n$-point Gauss–Legendre rule on $[-1, 1]$ (for exact variable-area mass matrices
/// and continuous $L^2$ error norms).
fn gauss_legendre_rule(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut nodes = vec![0.0; n];
    let mut weights = vec![0.0; n];
    let half = n.div_ceil(2);
    let nf = n as f64;
    for i in 0..half {
        let theta = core::f64::consts::PI * (i as f64 + 0.75) / (nf + 0.5);
        let mut xi = -det::cos(theta);
        for _ in 0..64 {
            let (p, dp, _) = legendre_p_and_derivs(n, xi);
            let step = p / dp;
            xi -= step;
            if step.abs() <= 1e-16 {
                break;
            }
        }
        let (_, dp, _) = legendre_p_and_derivs(n, xi);
        let w = 2.0 / ((1.0 - xi * xi) * dp * dp);
        nodes[i] = if 2 * i + 1 == n { 0.0 } else { xi };
        nodes[n - 1 - i] = -nodes[i];
        weights[i] = w;
        weights[n - 1 - i] = w;
    }
    (nodes, weights)
}

/// Compute $A^{1/2}$ and $A^{-1/2}$ for a symmetric positive-definite $m \times m$ matrix $A$.
fn spd_sqrt_and_inv_sqrt(a: &[f64], m: usize) -> Result<(Vec<f64>, Vec<f64>), DuctError> {
    let (evals, evecs) = jacobi_eigh(a, m);
    let mut a_sqrt = vec![0.0; m * m];
    let mut a_inv_sqrt = vec![0.0; m * m];
    for k in 0..m {
        let lam = evals[k];
        if !(lam.is_finite() && lam > 0.0) {
            return Err(DuctError::Singular);
        }
        let s = det::sqrt(lam);
        let inv_s = 1.0 / s;
        for i in 0..m {
            let vi = evecs[i * m + k];
            for j in 0..m {
                let vj = evecs[j * m + k];
                a_sqrt[i * m + j] += s * vi * vj;
                a_inv_sqrt[i * m + j] += inv_s * vi * vj;
            }
        }
    }
    Ok((a_sqrt, a_inv_sqrt))
}

/// Element mass matrices `M_p = ∫ A/(ρc²) ℓ_i ℓ_j dx` and `M_q = ∫ ρ/A ℓ_i ℓ_j dx` over
/// `[s_left, s_left + h]` of `span`, by Gauss–Legendre rules doubled from `2p + 4` points until
/// both repeat to `MASS_QUADRATURE_ULPS` ulps of their largest entry.
fn element_masses<E: From<DuctError>>(
    span: &AxialSpan,
    s_left: f64,
    h: f64,
    lgl: &LglRule,
    gas_at: &mut impl FnMut(f64) -> Result<GasState, E>,
) -> Result<(Vec<f64>, Vec<f64>), E> {
    let n_local = lgl.order + 1;
    let jac = 0.5 * h;
    let mut rule = |points: usize| -> Result<(Vec<f64>, Vec<f64>), E> {
        let (xi, w) = gauss_legendre_rule(points);
        let mut mp = vec![0.0; n_local * n_local];
        let mut mq = vec![0.0; n_local * n_local];
        for (&x, &wq) in xi.iter().zip(&w) {
            let s = s_left + jac * (x + 1.0);
            let r = span.radius_at_local(s)?;
            let area = core::f64::consts::PI * r * r;
            if !(area.is_finite() && area > 0.0) {
                return Err(DuctError::BadParameter {
                    what: "cross-sectional area A(x) must be positive and finite",
                }
                .into());
            }
            let gas = gas_at(span.x_start + s)?;
            let rho = gas.density;
            let bulk = rho * gas.sound_speed * gas.sound_speed;
            if !(rho.is_finite() && rho > 0.0 && bulk.is_finite() && bulk > 0.0) {
                return Err(DuctError::BadParameter {
                    what: "DG gas density and bulk modulus must be finite and positive",
                }
                .into());
            }
            let phi = lgl.eval_basis(x);
            let (cp, lq) = (jac * wq * (area / bulk), jac * wq * (rho / area));
            for i in 0..n_local {
                for j in 0..n_local {
                    mp[i * n_local + j] += cp * phi[i] * phi[j];
                    mq[i * n_local + j] += lq * phi[i] * phi[j];
                }
            }
        }
        Ok((mp, mq))
    };
    let settled = |a: &[f64], b: &[f64]| {
        let scale = b.iter().fold(0.0_f64, |acc, v| acc.max(v.abs()));
        let diff = a
            .iter()
            .zip(b)
            .fold(0.0_f64, |acc, (x, y)| acc.max((x - y).abs()));
        diff <= MASS_QUADRATURE_ULPS * f64::EPSILON * scale
    };
    let mut points = 2 * lgl.order + 4;
    let mut current = rule(points)?;
    loop {
        points *= 2;
        if points > MAX_MASS_QUADRATURE_POINTS {
            return Err(DuctError::BadParameter {
                what: "DG mass quadrature of A(x) did not converge",
            }
            .into());
        }
        let finer = rule(points)?;
        if settled(&current.0, &finer.0) && settled(&current.1, &finer.1) {
            return Ok(finer);
        }
        current = finer;
    }
}

#[derive(Debug, Clone)]
struct AxialSpan {
    x_start: f64,
    length: f64,
    /// Offset into the original geometry after a material station splits it.
    geometry_start: f64,
    segment: Segment,
}

impl AxialSpan {
    fn radius_at_local(&self, s: f64) -> Result<f64, DuctError> {
        let s_clamped = (self.geometry_start + s.clamp(0.0, self.length))
            .clamp(0.0, self.segment.axial_length());
        match &self.segment {
            Segment::Cylinder { radius, .. } => Ok(*radius),
            Segment::Cone {
                inlet_radius,
                outlet_radius,
                length,
            } => {
                let t = s_clamped / *length;
                Ok(*inlet_radius + (*outlet_radius - *inlet_radius) * t)
            }
            Segment::Profile { profile, .. } => profile.radius(s_clamped),
            Segment::ToneHole { .. } => Err(DuctError::BadParameter {
                what: "tone-hole side branches require M2.2-W3c coupling",
            }),
        }
    }
}

#[derive(Debug, Clone)]
struct DgElement {
    x_left: f64,
    length: f64,
    nodes_x: Vec<f64>,
    areas: Vec<f64>,
    mp_sqrt: Vec<f64>,
    mp_inv_sqrt: Vec<f64>,
    mq_sqrt: Vec<f64>,
    mq_inv_sqrt: Vec<f64>,
}

/// Assembled 1-D variable-area discontinuous-Galerkin duct, optionally carrying
/// admitted circular-tube viscothermal memory through the same exponential.
#[derive(Debug, Clone)]
pub struct DgDuct {
    operator_identity: ContentHash,
    mesh: DgMesh,
    gas: GasState,
    outlet_gas: GasState,
    gas_field: Option<AxialGasField>,
    /// Exact source geometry eligibility, not inferred from sampled node radii.
    /// Original geometry retained for material-point wall integration.
    geometry: Duct,
    uniform_radius_m: Option<f64>,
    viscothermal: Option<viscothermal::ViscothermalMemory>,
    molecular: Option<molecular::MolecularMemory>,
    compliant_wall: Option<compliant_wall::CompliantWall>,
    memory_dofs: usize,
    inlet: DgBoundary,
    outlet: DgBoundary,
    total_length: f64,
    elements: Vec<DgElement>,
    lgl: LglRule,
    /// Number of scalar DOFs per field ($M = K (p + 1)$).
    field_dofs: usize,
    /// Energy-weighted operator, including every admitted wall-memory coordinate.
    l_tilde: Vec<f64>,
    /// Energy-weighted off-diagonal block $G = M_p^{-1/2} K_p M_q^{-1/2}$ ($M \times M$, row-major)
    /// when boundaries are purely skew-symmetric (`Closed`, `IdealOpen`, `Port`).
    g_block: Option<Vec<f64>>,
    /// Active port descriptors `(PortId, vector_of_len_total_dofs)`.
    ports: Vec<(PortId, Vec<f64>)>,
    /// Point junctions on original element faces, in attachment order.
    side_ports: Vec<DgSidePort>,
}

/// Bind the physical coordinate maps as well as the energy-normalized operator.
/// In particular, changing density can preserve the normalized wave operator
/// while changing the meaning of every pressure/flow energy coordinate.
fn hash_operator(
    mesh: DgMesh,
    inlet: DgBoundary,
    outlet: DgBoundary,
    elements: &[DgElement],
    l_tilde: &[f64],
    ports: &[(PortId, Vec<f64>)],
) -> ContentHash {
    fn floats(hash: &mut DomainHasher, values: &[f64]) {
        hash.update(&(values.len() as u64).to_le_bytes());
        for value in values {
            hash.update(&value.to_bits().to_le_bytes());
        }
    }
    let mut hash = DomainHasher::new("org.frankensim.fs-duct.dg-operator.v1");
    for count in [mesh.elements, mesh.order, ports.len()] {
        hash.update(&(count as u64).to_le_bytes());
    }
    for boundary in [inlet, outlet] {
        let (tag, port) = match boundary {
            DgBoundary::Closed => (0, 0),
            DgBoundary::IdealOpen => (1, 0),
            DgBoundary::Matched => (2, 0),
            DgBoundary::Port(id) => (3, id.0),
        };
        hash.update(&[tag]);
        hash.update(&port.to_le_bytes());
    }
    for element in elements {
        floats(&mut hash, &[element.x_left, element.length]);
        for values in [
            &element.nodes_x,
            &element.areas,
            &element.mp_sqrt,
            &element.mp_inv_sqrt,
            &element.mq_sqrt,
            &element.mq_inv_sqrt,
        ] {
            floats(&mut hash, values);
        }
    }
    floats(&mut hash, l_tilde);
    for (id, vector) in ports {
        hash.update(&id.0.to_le_bytes());
        floats(&mut hash, vector);
    }
    hash.finalize()
}

impl DgDuct {
    /// Construct a nodal DG discretization of `duct` in `gas` with `mesh` and boundary conditions.
    ///
    /// # Errors
    /// - [`DuctError::EmptyDuct`] if `duct.segments` is empty.
    /// - [`DuctError::BadParameter`] if `mesh.order` is outside `1..=MAX_DG_ORDER`,
    ///   `mesh.elements` is outside `duct.segments.len()..=MAX_DG_ELEMENTS`, the DOF count
    ///   `2 K (p + 1)` exceeds [`MAX_DG_PROPAGATOR_DOFS`], the mass quadrature of `A(x)` does not
    ///   converge, any segment has
    ///   non-positive/non-finite radius or length, `Segment::ToneHole` is present, or two ports
    ///   share the same [`PortId`].
    pub fn new(
        duct: &Duct,
        gas: &GasState,
        mesh: DgMesh,
        inlet: DgBoundary,
        outlet: DgBoundary,
    ) -> Result<Self, DuctError> {
        Self::assemble(duct, gas, mesh, inlet, outlet, None, &mut |_| Ok(*gas))
    }

    /// Bind the shared primitive thermochemical field to the time-domain DG
    /// operator. Temperature and composition are reevaluated by the material
    /// owner at every mass quadrature point; derived coefficients are never
    /// interpolated. Matched boundaries use their own endpoint gas state.
    ///
    /// Geometry boundaries and field stations are exact element boundaries.
    /// `mesh.elements` is the total requested count, not a count per span; it
    /// must cover all resulting spans. A uniform field delegates [`Self::new`]
    /// bit-for-bit and retains the original material field for inspection.
    ///
    /// This is a frozen, lossless, no-mean-flow field. It does not evolve
    /// temperature or composition, add viscothermal losses, or model flow.
    ///
    /// # Errors
    /// Refuses unequal field/duct lengths, insufficient mesh resolution,
    /// material evaluations outside their source/phase domain, and the same
    /// geometry, quadrature and DOF limits as [`Self::new`].
    pub fn new_with_gas_field(
        duct: &Duct,
        field: &AxialGasField,
        mesh: DgMesh,
        inlet: DgBoundary,
        outlet: DgBoundary,
    ) -> Result<Self, GasFieldDuctError> {
        let length: f64 = duct.segments.iter().map(Segment::axial_length).sum();
        if !length.is_finite() || length <= 0.0 || length != field.length_m() {
            return Err(DuctError::BadParameter {
                what: "gas field must cover the exact complete DG duct length",
            }
            .into());
        }
        if let Some(state) = field.uniform_state() {
            let mut duct = Self::new(duct, state.state(), mesh, inlet, outlet)?;
            duct.gas_field = Some(field.clone());
            return Ok(duct);
        }
        let gas = *field.state_at(0.0)?.state();
        Self::assemble(duct, &gas, mesh, inlet, outlet, Some(field), &mut |x| {
            Ok(*field.state_at(x)?.state())
        })
    }

    fn assemble<E: From<DuctError>>(
        duct: &Duct,
        gas: &GasState,
        mesh: DgMesh,
        inlet: DgBoundary,
        outlet: DgBoundary,
        field: Option<&AxialGasField>,
        gas_at: &mut impl FnMut(f64) -> Result<GasState, E>,
    ) -> Result<Self, E> {
        if duct.segments.is_empty() {
            return Err(DuctError::EmptyDuct.into());
        }
        if mesh.order == 0 || mesh.order > MAX_DG_ORDER {
            return Err(DuctError::BadParameter {
                what: "DG order must be in 1..=MAX_DG_ORDER",
            }
            .into());
        }
        if mesh.elements == 0 || mesh.elements > MAX_DG_ELEMENTS {
            return Err(DuctError::BadParameter {
                what: "DG element count must be in 1..=MAX_DG_ELEMENTS",
            }
            .into());
        }
        if 2 * mesh.elements * (mesh.order + 1) > MAX_DG_PROPAGATOR_DOFS {
            return Err(DuctError::BadParameter {
                what: "DG DOFs 2 K (p + 1) exceed MAX_DG_PROPAGATOR_DOFS (no Krylov path)",
            }
            .into());
        }
        crate::validate_gas(gas)?;
        if let (DgBoundary::Port(id_in), DgBoundary::Port(id_out)) = (inlet, outlet)
            && id_in == id_out
        {
            return Err(DuctError::BadParameter {
                what: "inlet and outlet ports must have distinct PortIds",
            }
            .into());
        }

        let mut spans = Vec::with_capacity(duct.segments.len());
        let mut x_cursor = 0.0;
        for seg in &duct.segments {
            let len = match seg {
                Segment::Cylinder { radius, length } => {
                    if !(*radius > 0.0 && radius.is_finite() && *length > 0.0 && length.is_finite())
                    {
                        return Err(DuctError::BadParameter {
                            what: "cylinder radius and length must be positive and finite",
                        }
                        .into());
                    }
                    *length
                }
                Segment::Cone {
                    inlet_radius,
                    outlet_radius,
                    length,
                } => {
                    if !(*inlet_radius > 0.0
                        && inlet_radius.is_finite()
                        && *outlet_radius > 0.0
                        && outlet_radius.is_finite()
                        && *length > 0.0
                        && length.is_finite())
                    {
                        return Err(DuctError::BadParameter {
                            what: "cone radii and length must be positive and finite",
                        }
                        .into());
                    }
                    *length
                }
                Segment::Profile { profile, .. } => {
                    let len = profile.length();
                    if !(len > 0.0 && len.is_finite()) {
                        return Err(DuctError::BadParameter {
                            what: "profile length must be positive and finite",
                        }
                        .into());
                    }
                    len
                }
                Segment::ToneHole { .. } => {
                    return Err(DuctError::BadParameter {
                        what: "tone-hole side branches require M2.2-W3c coupling",
                    }
                    .into());
                }
            };
            spans.push(AxialSpan {
                x_start: x_cursor,
                length: len,
                geometry_start: 0.0,
                segment: seg.clone(),
            });
            x_cursor += len;
        }
        let total_length = x_cursor;
        if !total_length.is_finite() || total_length <= 0.0 {
            return Err(DuctError::BadParameter {
                what: "total DG duct length must be finite and positive",
            }
            .into());
        }
        if let Some(field) = field {
            let mut split_spans = Vec::new();
            for span in spans {
                let end = span.x_start + span.length;
                let mut left = span.x_start;
                for right in field
                    .stations_m()
                    .iter()
                    .copied()
                    .filter(|x| *x > span.x_start && *x < end)
                    .chain(core::iter::once(end))
                {
                    split_spans.push(AxialSpan {
                        x_start: left,
                        length: right - left,
                        geometry_start: left - span.x_start,
                        segment: span.segment.clone(),
                    });
                    left = right;
                }
            }
            spans = split_spans;
        }
        if mesh.elements < spans.len() {
            return Err(DuctError::BadParameter {
                what: "DG mesh must have at least one element per geometry/material span",
            }
            .into());
        }

        // Partition mesh.elements proportionally across spans, guaranteeing >= 1 element per span.
        let n_spans = spans.len();
        let mut elem_counts = vec![1usize; n_spans];
        let remaining = mesh.elements - n_spans;
        if remaining > 0 {
            let mut assigned = 0usize;
            let mut remainders = Vec::with_capacity(n_spans);
            for (idx, sp) in spans.iter().enumerate() {
                let exact = (remaining as f64) * (sp.length / total_length);
                let whole = exact.floor() as usize;
                elem_counts[idx] += whole;
                assigned += whole;
                remainders.push((exact - (whole as f64), idx));
            }
            remainders.sort_by(|a, b| b.0.total_cmp(&a.0).then(a.1.cmp(&b.1)));
            for &(_, idx) in remainders.iter().take(remaining - assigned) {
                elem_counts[idx] += 1;
            }
        }

        let lgl = lgl_rule(mesh.order)?;
        let n_local = mesh.order + 1;

        let outlet_gas = gas_at(total_length)?;

        let mut elements = Vec::with_capacity(mesh.elements);
        for (span_index, sp) in spans.iter().enumerate() {
            let k_span = elem_counts[span_index];
            let h = sp.length / (k_span as f64);
            let jac = 0.5 * h;
            for sub in 0..k_span {
                let s_left = (sub as f64) * h;
                let x_left = sp.x_start + s_left;
                let mut nodes_x = Vec::with_capacity(n_local);
                let mut areas = Vec::with_capacity(n_local);
                for &xi in &lgl.nodes {
                    let s_loc = s_left + jac * (xi + 1.0);
                    let r = sp.radius_at_local(s_loc)?;
                    let area = core::f64::consts::PI * r * r;
                    if !(area.is_finite() && area > 0.0) {
                        return Err(DuctError::BadParameter {
                            what: "cross-sectional area A(x) must be positive and finite",
                        }
                        .into());
                    }
                    nodes_x.push(sp.x_start + s_loc);
                    areas.push(area);
                }

                let (mp, mq) = element_masses(sp, s_left, h, &lgl, gas_at)?;
                let (mp_sqrt, mp_inv_sqrt) = spd_sqrt_and_inv_sqrt(&mp, n_local)?;
                let (mq_sqrt, mq_inv_sqrt) = spd_sqrt_and_inv_sqrt(&mq, n_local)?;

                elements.push(DgElement {
                    x_left,
                    length: h,
                    nodes_x,
                    areas,
                    mp_sqrt,
                    mp_inv_sqrt,
                    mq_sqrt,
                    mq_inv_sqrt,
                });
            }
        }

        let field_dofs = mesh.elements * n_local;
        let total_dofs = 2 * field_dofs;

        // Reference stiffness matrix S_{ij} = w_i D_{ij} on [-1, 1].
        let mut s_ref = vec![0.0; n_local * n_local];
        for i in 0..n_local {
            for j in 0..n_local {
                s_ref[i * n_local + j] = lgl.weights[i] * lgl.diff_matrix[i * n_local + j];
            }
        }

        // Assemble unscaled weak operators:
        // M_p dp/dt = -K_p q - D_pp p + B_p u
        // M_q dq/dt = -K_q p - D_qq q + B_q u
        // Using the alternating interior flux:
        //   q_hat_{e+1/2} = q_R^{(e)}  (from left)
        //   p_hat_{e+1/2} = p_L^{(e+1)} (from right)
        let mut kp = vec![0.0; field_dofs * field_dofs];
        let mut kq = vec![0.0; field_dofs * field_dofs];
        let mut d_pp = vec![0.0; field_dofs * field_dofs];
        let mut d_qq = vec![0.0; field_dofs * field_dofs];

        for e in 0..mesh.elements {
            let row0 = e * n_local;
            for i in 0..n_local {
                for j in 0..n_local {
                    let s_ij = s_ref[i * n_local + j];
                    kp[(row0 + i) * field_dofs + (row0 + j)] += s_ij;
                    kq[(row0 + i) * field_dofs + (row0 + j)] += s_ij;
                }
            }
            // Right face of element e (node n_local - 1):
            let r_node = row0 + (n_local - 1);
            if e + 1 < mesh.elements {
                let next_l_node = (e + 1) * n_local;
                // q_hat_R = q_R^{(e)} => (q_hat_R - q_R^{(e)}) = 0
                // p_hat_R = p_L^{(e+1)} => +e_p (p_L^{(e+1)} - p_R^{(e)})
                kq[r_node * field_dofs + next_l_node] += 1.0;
                kq[r_node * field_dofs + r_node] -= 1.0;
            }
            // Left face of element e (node 0):
            let l_node = row0;
            if e > 0 {
                let prev_r_node = (e - 1) * n_local + (n_local - 1);
                // q_hat_L = q_R^{(e-1)} => -e_0 (q_R^{(e-1)} - q_L^{(e)})
                kp[l_node * field_dofs + prev_r_node] -= 1.0;
                kp[l_node * field_dofs + l_node] += 1.0;
                // p_hat_L = p_L^{(e)} => (p_hat_L - p_L^{(e)}) = 0
            }
        }

        let mut ports = Vec::new();

        // Inlet boundary at x = 0 (element 0, node 0, outward normal -1):
        // Weak boundary terms are -e_0 (q_hat_L - q_0) in K_p and -e_0 (p_hat_L - p_0) in K_q.
        let in_node = 0usize;
        match inlet {
            DgBoundary::Closed => {
                // q_hat_L = 0, p_hat_L = p_0 => -e_0 (0 - q_0) = +e_0 q_0 in K_p
                kp[in_node * field_dofs + in_node] += 1.0;
            }
            DgBoundary::IdealOpen => {
                // p_hat_L = 0, q_hat_L = q_0 => -e_0 (0 - p_0) = +e_0 p_0 in K_q
                kq[in_node * field_dofs + in_node] += 1.0;
            }
            DgBoundary::Matched => {
                // Upwind characteristic flux exiting at x = 0 (p_hat_L = 0.5*(p_0 - Zc*q_0), q_hat_L = -0.5*(p_0/Zc - q_0)):
                let zc0 = gas.density * gas.sound_speed / elements[0].areas[0];
                kp[in_node * field_dofs + in_node] += 0.5;
                d_pp[in_node * field_dofs + in_node] += 0.5 / zc0;
                kq[in_node * field_dofs + in_node] += 0.5;
                d_qq[in_node * field_dofs + in_node] += 0.5 * zc0;
            }
            DgBoundary::Port(port_id) => {
                // Homogeneous part is Closed (q_hat_L = u_in, p_hat_L = p_0):
                // M_p dp/dt has -(-e_0 q_hat_L) = +e_0 u_in, so tilde_b_p = M_p^{-1/2} e_0.
                kp[in_node * field_dofs + in_node] += 1.0;
                let mut tilde_b = vec![0.0; total_dofs];
                for i in 0..n_local {
                    tilde_b[i] = elements[0].mp_inv_sqrt[i * n_local];
                }
                ports.push((port_id, tilde_b));
            }
        }

        // Outlet boundary at x = L (element K - 1, node p, outward normal +1):
        // Weak boundary terms are +e_p (q_hat_R - q_R) in K_p and +e_p (p_hat_R - p_R) in K_q.
        let last_elem = mesh.elements - 1;
        let out_node = last_elem * n_local + (n_local - 1);
        match outlet {
            DgBoundary::Closed => {
                // q_hat_R = 0, p_hat_R = p_R => +e_p (0 - q_R) = -e_p q_R in K_p
                kp[out_node * field_dofs + out_node] -= 1.0;
            }
            DgBoundary::IdealOpen => {
                // p_hat_R = 0, q_hat_R = q_R => +e_p (0 - p_R) = -e_p p_R in K_q
                kq[out_node * field_dofs + out_node] -= 1.0;
            }
            DgBoundary::Matched => {
                let zc_l = outlet_gas.density * outlet_gas.sound_speed
                    / elements[last_elem].areas[n_local - 1];
                kp[out_node * field_dofs + out_node] -= 0.5;
                d_pp[out_node * field_dofs + out_node] += 0.5 / zc_l;
                kq[out_node * field_dofs + out_node] -= 0.5;
                d_qq[out_node * field_dofs + out_node] += 0.5 * zc_l;
            }
            DgBoundary::Port(port_id) => {
                // Inflow into domain at x = L is u_out = -q_hat_R, p_hat_R = p_R:
                // M_p dp/dt has -e_p q_hat_R = +e_p u_out, so tilde_b_p = M_p^{-1/2} e_p.
                kp[out_node * field_dofs + out_node] -= 1.0;
                let mut tilde_b = vec![0.0; total_dofs];
                let base = last_elem * n_local;
                for i in 0..n_local {
                    tilde_b[base + i] =
                        elements[last_elem].mp_inv_sqrt[i * n_local + (n_local - 1)];
                }
                ports.push((port_id, tilde_b));
            }
        }

        // Transform to energy-normalized coordinates:
        // G = M_p^{-1/2} K_p M_q^{-1/2},  H = M_q^{-1/2} K_q M_p^{-1/2}
        let g_mat = Self::scale_operator(&elements, n_local, field_dofs, &kp, true);
        let h_mat = Self::scale_operator(&elements, n_local, field_dofs, &kq, false);
        let d_pp_tilde = Self::scale_diag_operator(&elements, n_local, field_dofs, &d_pp, true);
        let d_qq_tilde = Self::scale_diag_operator(&elements, n_local, field_dofs, &d_qq, false);

        let is_lossless_bcs =
            !matches!(inlet, DgBoundary::Matched) && !matches!(outlet, DgBoundary::Matched);
        if is_lossless_bcs {
            Self::admit_lossless_transpose_pair(&g_mat, &h_mat, field_dofs, n_local)?;
        }
        let mut l_tilde = vec![0.0; total_dofs * total_dofs];
        for i in 0..field_dofs {
            for j in 0..field_dofs {
                l_tilde[i * total_dofs + j] = -d_pp_tilde[i * field_dofs + j];
                l_tilde[i * total_dofs + (field_dofs + j)] = -g_mat[i * field_dofs + j];
                // The existing SkewModal time owner propagates [0,-G;G^T,0].
                // Store that SAME generator, after checking the independently
                // assembled SBP pair, so roundoff is not exported as indefinite
                // physical damping to a continuous or composed-time consumer.
                l_tilde[(field_dofs + i) * total_dofs + j] = if is_lossless_bcs {
                    g_mat[j * field_dofs + i]
                } else {
                    -h_mat[i * field_dofs + j]
                };
                l_tilde[(field_dofs + i) * total_dofs + (field_dofs + j)] =
                    -d_qq_tilde[i * field_dofs + j];
            }
        }

        let g_block = is_lossless_bcs.then_some(g_mat);
        let operator_identity = hash_operator(mesh, inlet, outlet, &elements, &l_tilde, &ports);

        Ok(Self {
            operator_identity,
            mesh,
            gas: *gas,
            outlet_gas,
            gas_field: field.cloned(),
            geometry: duct.clone(),
            uniform_radius_m: viscothermal::uniform_radius(duct),
            viscothermal: None,
            molecular: None,
            compliant_wall: None,
            memory_dofs: 0,
            inlet,
            outlet,
            total_length,
            elements,
            lgl,
            field_dofs,
            l_tilde,
            g_block,
            ports,
            side_ports: Vec::new(),
        })
    }

    /// Check the SBP identity H=-G^T before selecting the original lossless
    /// time owner's exact transpose. Two local mass-normalization products and
    /// the LGL/face arithmetic get a conservative 64*n_local^2*epsilon normwise
    /// roundoff allowance (n_local <= 9); larger defects are source errors, not
    /// damping to project away. Scaling before subtraction also avoids overflow.
    fn admit_lossless_transpose_pair(
        g: &[f64],
        h: &[f64],
        field_dofs: usize,
        n_local: usize,
    ) -> Result<(), DuctError> {
        let mut scale = 0.0_f64;
        for &value in g.iter().chain(h) {
            if !value.is_finite() {
                return Err(DuctError::BadParameter {
                    what: "lossless DG weak operators must be finite",
                });
            }
            scale = scale.max(value.abs());
        }
        if scale == 0.0 {
            return Ok(());
        }
        let allowance = 64.0 * (n_local * n_local) as f64 * f64::EPSILON;
        for i in 0..field_dofs {
            for j in 0..field_dofs {
                let defect = h[i * field_dofs + j] / scale + g[j * field_dofs + i] / scale;
                if defect.abs() > allowance {
                    return Err(DuctError::BadParameter {
                        what: "lossless DG weak operators exceed the roundoff SBP transpose allowance",
                    });
                }
            }
        }
        Ok(())
    }

    fn scale_operator(
        elements: &[DgElement],
        n_local: usize,
        field_dofs: usize,
        op: &[f64],
        left_is_p: bool,
    ) -> Vec<f64> {
        let mut temp = vec![0.0; field_dofs * field_dofs];
        // Multiply on the left by block-diagonal M_left^{-1/2}:
        for (e_row, elem) in elements.iter().enumerate() {
            let inv_sqrt = if left_is_p {
                &elem.mp_inv_sqrt
            } else {
                &elem.mq_inv_sqrt
            };
            let r0 = e_row * n_local;
            for i in 0..n_local {
                for col in 0..field_dofs {
                    let mut acc = 0.0;
                    for k in 0..n_local {
                        acc += inv_sqrt[i * n_local + k] * op[(r0 + k) * field_dofs + col];
                    }
                    temp[(r0 + i) * field_dofs + col] = acc;
                }
            }
        }
        let mut out = vec![0.0; field_dofs * field_dofs];
        // Multiply on the right by block-diagonal M_right^{-1/2}:
        for (e_col, elem) in elements.iter().enumerate() {
            let inv_sqrt = if left_is_p {
                &elem.mq_inv_sqrt
            } else {
                &elem.mp_inv_sqrt
            };
            let c0 = e_col * n_local;
            for row in 0..field_dofs {
                for j in 0..n_local {
                    let mut acc = 0.0;
                    for k in 0..n_local {
                        acc += temp[row * field_dofs + (c0 + k)] * inv_sqrt[k * n_local + j];
                    }
                    out[row * field_dofs + (c0 + j)] = acc;
                }
            }
        }
        out
    }

    fn scale_diag_operator(
        elements: &[DgElement],
        n_local: usize,
        field_dofs: usize,
        op: &[f64],
        is_p: bool,
    ) -> Vec<f64> {
        let mut temp = vec![0.0; field_dofs * field_dofs];
        for (e_row, elem) in elements.iter().enumerate() {
            let inv_sqrt = if is_p {
                &elem.mp_inv_sqrt
            } else {
                &elem.mq_inv_sqrt
            };
            let r0 = e_row * n_local;
            for i in 0..n_local {
                for col in 0..field_dofs {
                    let mut acc = 0.0;
                    for k in 0..n_local {
                        acc += inv_sqrt[i * n_local + k] * op[(r0 + k) * field_dofs + col];
                    }
                    temp[(r0 + i) * field_dofs + col] = acc;
                }
            }
        }
        let mut out = vec![0.0; field_dofs * field_dofs];
        for (e_col, elem) in elements.iter().enumerate() {
            let inv_sqrt = if is_p {
                &elem.mp_inv_sqrt
            } else {
                &elem.mq_inv_sqrt
            };
            let c0 = e_col * n_local;
            for row in 0..field_dofs {
                for j in 0..n_local {
                    let mut acc = 0.0;
                    for k in 0..n_local {
                        acc += temp[row * field_dofs + (c0 + k)] * inv_sqrt[k * n_local + j];
                    }
                    out[row * field_dofs + (c0 + j)] = acc;
                }
            }
        }
        out
    }

    /// Mesh parameters `(elements, order)`.
    #[must_use]
    pub const fn mesh(&self) -> DgMesh {
        self.mesh
    }

    /// Total scalar DOFs, including all retained viscothermal memory.
    #[must_use]
    pub const fn total_dofs(&self) -> usize {
        2 * self.field_dofs + self.memory_dofs
    }

    /// Scalar DOFs per field $M = K (p + 1)$.
    #[must_use]
    pub const fn field_dofs(&self) -> usize {
        self.field_dofs
    }

    /// Total axial length $L$ [m].
    #[must_use]
    pub const fn total_length(&self) -> f64 {
        self.total_length
    }

    /// Inlet reference state. This is the uniform medium for [`Self::new`];
    /// a field-backed duct instead uses its local material at each point.
    #[must_use]
    pub const fn gas(&self) -> &GasState {
        &self.gas
    }

    /// Exact inlet state used by the inlet characteristic boundary.
    #[must_use]
    pub const fn inlet_gas(&self) -> &GasState {
        &self.gas
    }

    /// Exact outlet state used by the outlet characteristic boundary.
    #[must_use]
    pub const fn outlet_gas(&self) -> &GasState {
        &self.outlet_gas
    }

    /// Physical inlet area from the original boundary geometry [m²].
    #[must_use]
    pub fn inlet_area_m2(&self) -> f64 {
        self.elements[0].areas[0]
    }

    /// Physical outlet area from the original boundary geometry [m²].
    #[must_use]
    pub fn outlet_area_m2(&self) -> f64 {
        let element = &self.elements[self.elements.len() - 1];
        element.areas[element.areas.len() - 1]
    }

    /// Original primitive field and material/source receipts, when supplied.
    #[must_use]
    pub const fn gas_field(&self) -> Option<&AxialGasField> {
        self.gas_field.as_ref()
    }

    /// Instantaneous conjugate pressures of every actual declared port [Pa].
    /// This uses the current native pressure loading columns and physical state,
    /// including after a material update; it neither advances time nor infers a
    /// new integrated port work. A foreign state or nonfinite observation refuses.
    pub fn port_pressures(&self, state: &DgState) -> Result<Vec<(PortId, f64)>, DuctError> {
        if state.operator_identity != self.operator_identity()
            || state.tilde.len() != self.total_dofs()
        {
            return Err(DuctError::BadParameter { what: "DG port observation requires the current operator state" });
        }
        self.ports.iter().map(|(id, column)| {
            let pressure: f64 = column.iter().zip(&state.tilde).map(|(b, x)| b * x).sum();
            if !pressure.is_finite() {
                return Err(DuctError::BadParameter { what: "DG port pressure is unrepresentable" });
            }
            Ok((*id, pressure))
        }).collect()
    }

    /// Inlet boundary condition at $x = 0$.
    #[must_use]
    pub const fn inlet(&self) -> DgBoundary {
        self.inlet
    }

    /// Outlet boundary condition at $x = L$.
    #[must_use]
    pub const fn outlet(&self) -> DgBoundary {
        self.outlet
    }

    /// Physical axial coordinates $x_i^{(e)}$ [m] of all $M = K (p + 1)$ LGL nodes in element order.
    #[must_use]
    pub fn node_positions(&self) -> Vec<f64> {
        let mut out = Vec::with_capacity(self.field_dofs);
        for elem in &self.elements {
            out.extend_from_slice(&elem.nodes_x);
        }
        out
    }

    /// Map physical nodal fields $(p, q)$ to energy-weighted coordinates
    /// $\tilde{x} = [M_p^{1/2} p; M_q^{1/2} q; 0_{memory}]$.
    /// Wall memories start at zero; this is initialization, not a state transfer.
    #[must_use]
    pub fn to_energy_coords(&self, p: &[f64], q: &[f64]) -> Vec<f64> {
        assert_eq!(p.len(), self.field_dofs);
        assert_eq!(q.len(), self.field_dofs);
        let n_local = self.mesh.order + 1;
        let mut tilde = vec![0.0; self.total_dofs()];
        for (e, elem) in self.elements.iter().enumerate() {
            let base = e * n_local;
            for i in 0..n_local {
                let mut sp = 0.0;
                let mut sq = 0.0;
                for j in 0..n_local {
                    sp += elem.mp_sqrt[i * n_local + j] * p[base + j];
                    sq += elem.mq_sqrt[i * n_local + j] * q[base + j];
                }
                tilde[base + i] = sp;
                tilde[self.field_dofs + base + i] = sq;
            }
        }
        tilde
    }

    /// Map energy-weighted coordinates $\tilde{x}$ back to physical nodal fields $(p, q)$.
    #[must_use]
    pub fn from_energy_coords(&self, tilde: &[f64]) -> (Vec<f64>, Vec<f64>) {
        assert_eq!(tilde.len(), self.total_dofs());
        let n_local = self.mesh.order + 1;
        let mut p = vec![0.0; self.field_dofs];
        let mut q = vec![0.0; self.field_dofs];
        for (e, elem) in self.elements.iter().enumerate() {
            let base = e * n_local;
            for i in 0..n_local {
                let mut sp = 0.0;
                let mut sq = 0.0;
                for j in 0..n_local {
                    sp += elem.mp_inv_sqrt[i * n_local + j] * tilde[base + j];
                    sq += elem.mq_inv_sqrt[i * n_local + j] * tilde[self.field_dofs + base + j];
                }
                p[base + i] = sp;
                q[base + i] = sq;
            }
        }
        (p, q)
    }

    /// Initialize a [`DgState`] by element-wise Gauss–Radau projection $(\Pi^- p_0, \Pi^+ q_0)$
    /// compatible with the alternating numerical flux ($\hat{p}_{e-1/2} = p_L^{(e)}$,
    /// $\hat{q}_{e+1/2} = q_R^{(e)}$), using $(2p + 4)$-node Gauss–Legendre quadrature.
    #[must_use]
    pub fn project_initial_state<Fp, Fq>(&self, mut f_p: Fp, mut f_q: Fq) -> DgState
    where
        Fp: FnMut(f64) -> f64,
        Fq: FnMut(f64) -> f64,
    {
        let p_ord = self.mesh.order;
        let n_local = p_ord + 1;
        let (gl_xi, gl_w) = gauss_legendre_rule(2 * p_ord + 4);
        // Precompute Legendre polynomials P_k(xi_q) at Gauss-Legendre quadrature points for k = 0..p_ord,
        // and P_k(xi_i) at LGL nodes for k = 0..=p_ord.
        let gl_leg: Vec<Vec<f64>> = gl_xi
            .iter()
            .map(|&xi| (0..p_ord).map(|k| legendre_p_and_derivs(k, xi).0).collect())
            .collect();
        let lgl_leg: Vec<Vec<f64>> = self
            .lgl
            .nodes
            .iter()
            .map(|&xi| {
                (0..=p_ord)
                    .map(|k| legendre_p_and_derivs(k, xi).0)
                    .collect()
            })
            .collect();

        let mut p = vec![0.0; self.field_dofs];
        let mut q = vec![0.0; self.field_dofs];
        for (e, elem) in self.elements.iter().enumerate() {
            let jac = 0.5 * elem.length;
            let mut c_p = vec![0.0; n_local];
            let mut c_q = vec![0.0; n_local];
            for (q_idx, (&xi, &wq)) in gl_xi.iter().zip(gl_w.iter()).enumerate() {
                let x = elem.x_left + jac * (xi + 1.0);
                let vp = f_p(x);
                let vq = f_q(x);
                let leg = &gl_leg[q_idx];
                for k in 0..p_ord {
                    let scale = 0.5 * (2 * k + 1) as f64 * wq;
                    c_p[k] += scale * leg[k] * vp;
                    c_q[k] += scale * leg[k] * vq;
                }
            }
            // Pin Pi^- p at xi = -1 (x_left) and Pi^+ q at xi = +1 (x_right):
            let p_left = f_p(elem.x_left);
            let q_right = f_q(elem.x_left + elem.length);
            let mut sum_p_left = 0.0;
            let mut sum_q_right = 0.0;
            for k in 0..p_ord {
                let sign = if k % 2 == 0 { 1.0 } else { -1.0 };
                sum_p_left += c_p[k] * sign;
                sum_q_right += c_q[k];
            }
            let sign_p = if p_ord % 2 == 0 { 1.0 } else { -1.0 };
            c_p[p_ord] = sign_p * (p_left - sum_p_left);
            c_q[p_ord] = q_right - sum_q_right;

            let base = e * n_local;
            for (i, leg_i) in lgl_leg.iter().enumerate() {
                let mut acc_p = 0.0;
                let mut acc_q = 0.0;
                for k in 0..=p_ord {
                    acc_p += c_p[k] * leg_i[k];
                    acc_q += c_q[k] * leg_i[k];
                }
                p[base + i] = acc_p;
                q[base + i] = acc_q;
            }
        }
        let tilde = self.to_energy_coords(&p, &q);
        let e0 = 0.5 * tilde.iter().map(|v| v * v).sum::<f64>();
        DgState {
            operator_identity: self.operator_identity,
            p,
            q,
            tilde,
            modal_anchor_p: Vec::new(),
            modal_anchor_q: Vec::new(),
            modal_anchor_dt_bits: None,
            unforced_steps: 0,
            step_count: 0,
            initial_energy_j: e0,
            cumulative_port_work_j: 0.0,
            cumulative_parameter_work_j: 0.0,
            cumulative_dissipated_j: 0.0,
            nodal_dirty: false,
        }
    }

    /// Identity of the assembled spatial operator, physical mass transforms,
    /// geometry, and typed boundaries. Independent identical assemblies share it.
    #[must_use]
    pub const fn operator_identity(&self) -> ContentHash {
        self.operator_identity
    }

    /// Initialize a zero [`DgState`].
    #[must_use]
    pub fn zero_state(&self) -> DgState {
        DgState {
            operator_identity: self.operator_identity,
            p: vec![0.0; self.field_dofs],
            q: vec![0.0; self.field_dofs],
            tilde: vec![0.0; self.total_dofs()],
            modal_anchor_p: Vec::new(),
            modal_anchor_q: Vec::new(),
            modal_anchor_dt_bits: None,
            unforced_steps: 0,
            step_count: 0,
            initial_energy_j: 0.0,
            cumulative_port_work_j: 0.0,
            cumulative_parameter_work_j: 0.0,
            cumulative_dissipated_j: 0.0,
            nodal_dirty: false,
        }
    }

    /// Compute the continuous $L^2(0, L)$ pressure error
    /// $\left(\int_0^L (p_h(x) - p_{\text{exact}}(x))^2\,dx\right)^{1/2}$ using
    /// $(2p + 4)$-node Gauss–Legendre quadrature over each element.
    ///
    /// # Errors
    /// Refuses a state bound to a different spatial operator or physical mass map.
    pub fn l2_pressure_error<F>(&self, state: &DgState, mut p_exact: F) -> Result<f64, DuctError>
    where
        F: FnMut(f64) -> f64,
    {
        if state.operator_identity != self.operator_identity {
            return Err(DuctError::BadParameter {
                what: "DG state does not belong to this duct",
            });
        }
        let p_nodal = if state.nodal_dirty {
            let (p, _) = self.from_energy_coords(&state.tilde);
            p
        } else {
            state.p.clone()
        };
        let n_local = self.mesh.order + 1;
        let (gl_xi, gl_w) = gauss_legendre_rule(2 * self.mesh.order + 4);
        let gl_basis: Vec<Vec<f64>> = gl_xi.iter().map(|&xi| self.lgl.eval_basis(xi)).collect();
        let mut err_sq = 0.0;
        for (e, elem) in self.elements.iter().enumerate() {
            let jac = 0.5 * elem.length;
            let base = e * n_local;
            for (q_idx, (&xi, &wq)) in gl_xi.iter().zip(gl_w.iter()).enumerate() {
                let x = elem.x_left + jac * (xi + 1.0);
                let phi = &gl_basis[q_idx];
                let mut ph = 0.0;
                for i in 0..n_local {
                    ph += p_nodal[base + i] * phi[i];
                }
                let diff = ph - p_exact(x);
                err_sq += jac * wq * diff * diff;
            }
        }
        Ok(det::sqrt(err_sq))
    }

    /// Build the certified exact-exponential propagator for step `dt` [s] with certificate
    /// ceiling `tol`.
    ///
    /// Lossless boundaries (`Closed`, `IdealOpen`, `Port`) make the energy-coordinate operator
    /// skew. It is diagonalized through the SVD of its coupling block `G`, taken from the
    /// symmetric Jordan–Wielandt matrix `[[0, G], [G^T, 0]]` (eigenpairs `±σ`, `(u; ±v)/√2`),
    /// so no condition number is squared. The certificate is the larger of the relative
    /// singular residual `sqrt(‖GV − UΣ‖_F² + ‖GᵀU − VΣ‖_F²) / (√2 ‖G‖_F)` and the
    /// orthogonality defects `‖UᵀU − I‖_F`, `‖VᵀV − I‖_F`.
    ///
    /// A `Matched` boundary or attached viscothermal memory makes the operator dissipative
    /// and non-normal. Its exponential is
    /// a scaled Taylor series with squaring. The certificate is the larger of the rigorous
    /// truncation bound carried through the squarings in the energy norm (where the
    /// semigroup is contractive) and an a-posteriori estimate, the change under one more
    /// halving. That estimate is not an enclosure.
    ///
    /// # Errors
    /// [`DuctError::BadParameter`]: `dt` or `tol` not positive and finite, more than
    /// [`MAX_DG_PROPAGATOR_DOFS`] DOFs (there is no Krylov path), or the dense eigensolver's
    /// work cap. [`DuctError::PropagatorCertificate`]: the certificate exceeds `tol`.
    pub fn propagator(&self, dt: f64, tol: f64) -> Result<DgPropagator, DuctError> {
        if !(dt.is_finite() && dt > 0.0) {
            return Err(DuctError::BadParameter {
                what: "DG propagator step dt must be positive and finite",
            });
        }
        if !(tol.is_finite() && tol > 0.0) {
            return Err(DuctError::BadParameter {
                what: "DG propagator certificate tolerance must be positive and finite",
            });
        }
        if self.total_dofs() > MAX_DG_PROPAGATOR_DOFS {
            return Err(DuctError::BadParameter {
                what: "DG propagator DOFs exceed MAX_DG_PROPAGATOR_DOFS (no Krylov path)",
            });
        }
        match &self.g_block {
            Some(g) => self.build_skew_propagator(g, dt, tol),
            None => self.build_series_propagator(dt, tol),
        }
    }

    fn build_skew_propagator(
        &self,
        g: &[f64],
        dt: f64,
        tol: f64,
    ) -> Result<DgPropagator, DuctError> {
        let m = self.field_dofs;
        let n = 2 * m;
        admit_jacobi_eigh(n).map_err(|_| DuctError::BadParameter {
            what: "DG propagator exceeds the dense symmetric eigensolver's work cap",
        })?;
        let mut jw = vec![0.0; n * n];
        for i in 0..m {
            for k in 0..m {
                let v = g[i * m + k];
                jw[i * n + m + k] = v;
                jw[(m + k) * n + i] = v;
            }
        }
        let (lambda, w) = householder_eigh(&jw, n);
        let spectral = lambda.iter().fold(0.0_f64, |acc, l| acc.max(l.abs()));
        let threshold = NULL_SINGULAR_REL * spectral.max(f64::MIN_POSITIVE);
        let mut sigmas = Vec::with_capacity(m);
        let mut u_list: Vec<Vec<f64>> = Vec::with_capacity(m);
        let mut v_list: Vec<Vec<f64>> = Vec::with_capacity(m);
        let (mut null_u, mut null_v) = (Vec::new(), Vec::new());
        for k in 0..n {
            let u: Vec<f64> = (0..m).map(|i| w[i * n + k]).collect();
            let v: Vec<f64> = (0..m).map(|i| w[(m + i) * n + k]).collect();
            if lambda[k] > threshold {
                sigmas.push(lambda[k]);
                u_list.push(u.iter().map(|x| x * core::f64::consts::SQRT_2).collect());
                v_list.push(v.iter().map(|x| x * core::f64::consts::SQRT_2).collect());
            } else if lambda[k] >= -threshold {
                null_u.push(u);
                null_v.push(v);
            }
        }
        if sigmas.len() > m {
            return Err(DuctError::Singular);
        }
        let null_rank = m - sigmas.len();
        for (u, v) in orthonormal_span(&null_u, null_rank, m)?
            .into_iter()
            .zip(orthonormal_span(&null_v, null_rank, m)?)
        {
            sigmas.push(0.0);
            u_list.push(u);
            v_list.push(v);
        }
        let mut u_cols = vec![0.0; m * m];
        let mut v_cols = vec![0.0; m * m];
        for k in 0..m {
            for i in 0..m {
                u_cols[i * m + k] = u_list[k][i];
                v_cols[i * m + k] = v_list[k][i];
            }
        }
        orthonormalize_columns(&mut u_cols, m);
        orthonormalize_columns(&mut v_cols, m);
        // Rayleigh singular values on the final bases.
        for k in 0..m {
            let mut dot = 0.0;
            for i in 0..m {
                let gv: f64 = (0..m).map(|j| g[i * m + j] * v_cols[j * m + k]).sum();
                dot += u_cols[i * m + k] * gv;
            }
            sigmas[k] = dot;
        }

        let g_norm = det::sqrt(g.iter().map(|v| v * v).sum::<f64>()).max(f64::MIN_POSITIVE);
        let mut res_sq = 0.0;
        for i in 0..m {
            for k in 0..m {
                let mut gv = 0.0;
                let mut gtu = 0.0;
                for j in 0..m {
                    gv += g[i * m + j] * v_cols[j * m + k];
                    gtu += g[j * m + i] * u_cols[j * m + k];
                }
                let r1 = gv - u_cols[i * m + k] * sigmas[k];
                let r2 = gtu - v_cols[i * m + k] * sigmas[k];
                res_sq += r1 * r1 + r2 * r2;
            }
        }
        let residual = det::sqrt(res_sq) / (core::f64::consts::SQRT_2 * g_norm);
        let certificate = residual
            .max(orthogonality_defect(&u_cols, m))
            .max(orthogonality_defect(&v_cols, m));
        if !(certificate <= tol) {
            return Err(DuctError::PropagatorCertificate {
                residual: certificate,
                tol,
            });
        }

        let mut modal_ports = Vec::with_capacity(self.ports.len());
        for (port_id, tilde_b) in &self.ports {
            let mut bp = vec![0.0; m];
            let mut bq = vec![0.0; m];
            for k in 0..m {
                for i in 0..m {
                    bp[k] += u_cols[i * m + k] * tilde_b[i];
                    bq[k] += v_cols[i * m + k] * tilde_b[m + i];
                }
            }
            modal_ports.push((*port_id, bp, bq));
        }

        Ok(DgPropagator {
            operator_identity: self.operator_identity,
            dt,
            tol,
            certificate_residual: certificate,
            field_dofs: m,
            total_dofs: n,
            elements: self.elements.clone(),
            order: self.mesh.order,
            kernel: PropagatorKernel::SkewModal {
                sigmas,
                u_cols,
                v_cols,
                modal_ports,
            },
        })
    }

    fn build_series_propagator(&self, dt: f64, tol: f64) -> Result<DgPropagator, DuctError> {
        let m = self.field_dofs;
        let n = self.total_dofs();
        let shift = series_shift(&self.l_tilde, n, dt);
        let (phi0, phi1, phi2, truncation) = series_phis(&self.l_tilde, n, dt, shift)?;
        let (phi0_finer, _, _, _) = series_phis(&self.l_tilde, n, dt, shift + 1)?;
        let halving: f64 = det::sqrt(
            phi0.iter()
                .zip(&phi0_finer)
                .map(|(a, b)| (a - b) * (a - b))
                .sum::<f64>(),
        );
        // Rigorous: the scaled Taylor remainder (1-norm), carried to the 2-norm (sqrt n) and
        // through 2^s squarings of a contraction. Estimated: one more halving's change.
        let carried = f64::from(1u32 << shift) * det::sqrt(n as f64) * truncation;
        let certificate = carried.max(halving);
        if !(certificate <= tol) {
            return Err(DuctError::PropagatorCertificate {
                residual: certificate,
                tol,
            });
        }

        let mut port_cols_phi1 = Vec::with_capacity(self.ports.len());
        let mut port_cols_phi2 = Vec::with_capacity(self.ports.len());
        for (_, tilde_b) in &self.ports {
            let mut c1 = vec![0.0; n];
            let mut c2 = vec![0.0; n];
            for i in 0..n {
                for j in 0..n {
                    c1[i] += dt * phi1[i * n + j] * tilde_b[j];
                    c2[i] += dt * dt * phi2[i * n + j] * tilde_b[j];
                }
            }
            port_cols_phi1.push(c1);
            port_cols_phi2.push(c2);
        }

        Ok(DgPropagator {
            operator_identity: self.operator_identity,
            dt,
            tol,
            certificate_residual: certificate,
            field_dofs: m,
            total_dofs: n,
            elements: self.elements.clone(),
            order: self.mesh.order,
            kernel: PropagatorKernel::DenseSeries {
                phi0,
                dt_phi1: scale_matrix(&phi1, dt),
                port_cols_phi1,
                port_cols_phi2,
                ports: self.ports.clone(),
            },
        })
    }
}

/// Halvings `s` that bring `‖L dt / 2^s‖_1` to at most 1/2. The Taylor remainder there is
/// still far below roundoff (`series_phis` bounds it), and each halving saved is one squaring
/// less to amplify rounding.
fn series_shift(l: &[f64], n: usize, dt: f64) -> u32 {
    let mut norm1 = 0.0f64;
    for j in 0..n {
        let col: f64 = (0..n).map(|i| l[i * n + j].abs() * dt).sum();
        norm1 = norm1.max(col);
    }
    let mut shift = 0u32;
    while norm1 > 0.5 && shift < 30 {
        norm1 *= 0.5;
        shift += 1;
    }
    shift
}

/// `phi_0(L dt)`, `phi_1(L dt)`, `phi_2(L dt)` by the Taylor series of the matrix scaled by
/// `2^-shift`, then the exact doubling recurrences. Returns the series' rigorous 1-norm
/// remainder bound for the scaled exponential.
#[allow(clippy::type_complexity)]
fn series_phis(
    l: &[f64],
    n: usize,
    dt: f64,
    shift: u32,
) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>, f64), DuctError> {
    if shift > 30 {
        return Err(DuctError::BadParameter {
            what: "DG series propagator needs more than 30 halvings; refine dt",
        });
    }
    let dt_sub = dt / f64::from(1u32 << shift);
    let a_sub: Vec<f64> = l.iter().map(|v| v * dt_sub).collect();
    let mut a_norm = 0.0f64;
    for j in 0..n {
        a_norm = a_norm.max((0..n).map(|i| a_sub[i * n + j].abs()).sum());
    }
    let mut phi0 = eye_matrix(n);
    let mut phi1 = eye_matrix(n);
    let mut phi2 = scale_matrix(&eye_matrix(n), 0.5);
    let mut term = eye_matrix(n);
    let mut remainder = f64::INFINITY;
    let mut power = 1.0;
    let mut factorial = 1.0;
    for k in 1..=64usize {
        term = mat_mul(&term, &a_sub, n);
        let inv_k = 1.0 / (k as f64);
        let inv_k1 = 1.0 / ((k + 1) as f64);
        let inv_k2 = 1.0 / (((k + 1) * (k + 2)) as f64);
        for i in 0..(n * n) {
            // Multiplication above produced A^k/(k-1)!. All three phi
            // series require A^k/k!, including their held-input integrals.
            term[i] *= inv_k;
            phi0[i] += term[i];
            phi1[i] += term[i] * inv_k1;
            phi2[i] += term[i] * inv_k2;
        }
        // ‖Σ_{j>k} A^j/j!‖ ≤ ‖A‖^{k+1}/(k+1)! · 1/(1 − ‖A‖/(k+2)).
        power *= a_norm;
        factorial *= k as f64;
        let next = power * a_norm / (factorial * (k + 1) as f64);
        remainder = next / (1.0 - a_norm / (k + 2) as f64);
        if k >= 12 && remainder <= 1e-18 {
            break;
        }
    }
    for _ in 0..shift {
        let next_phi0 = mat_mul(&phi0, &phi0, n);
        let mut i_plus_phi0 = phi0.clone();
        for d in 0..n {
            i_plus_phi0[d * n + d] += 1.0;
        }
        let next_phi1 = scale_matrix(&mat_mul(&i_plus_phi0, &phi1, n), 0.5);
        let phi1_sq = mat_mul(&phi1, &phi1, n);
        let next_phi2: Vec<f64> = phi1_sq
            .iter()
            .zip(&phi2)
            .map(|(a, b)| 0.25 * (a + 2.0 * b))
            .collect();
        phi0 = next_phi0;
        phi1 = next_phi1;
        phi2 = next_phi2;
    }
    Ok((phi0, phi1, phi2, remainder))
}

/// An orthonormal basis of `rank` vectors for the span of `vectors` (length `m`), by
/// two-pass modified Gram–Schmidt in order; candidates whose remainder falls below 1e-3 of
/// their norm are skipped.
fn orthonormal_span(
    vectors: &[Vec<f64>],
    rank: usize,
    m: usize,
) -> Result<Vec<Vec<f64>>, DuctError> {
    let mut basis: Vec<Vec<f64>> = Vec::with_capacity(rank);
    for candidate in vectors {
        if basis.len() == rank {
            break;
        }
        let norm0 = det::sqrt(candidate.iter().map(|x| x * x).sum::<f64>());
        if norm0 == 0.0 {
            continue;
        }
        let mut v = candidate.clone();
        for _pass in 0..2 {
            for b in &basis {
                let dot: f64 = (0..m).map(|i| b[i] * v[i]).sum();
                for i in 0..m {
                    v[i] -= dot * b[i];
                }
            }
        }
        let norm = det::sqrt(v.iter().map(|x| x * x).sum::<f64>());
        if norm > 1e-3 * norm0 {
            basis.push(v.iter().map(|x| x / norm).collect());
        }
    }
    if basis.len() == rank {
        Ok(basis)
    } else {
        Err(DuctError::Singular)
    }
}

/// `‖CᵀC − I‖_F` for the `m` columns of the row-major `m × m` matrix `C`.
fn orthogonality_defect(cols: &[f64], m: usize) -> f64 {
    let mut sum = 0.0;
    for a in 0..m {
        for b in 0..m {
            let dot: f64 = (0..m).map(|i| cols[i * m + a] * cols[i * m + b]).sum();
            let target = if a == b { 1.0 } else { 0.0 };
            sum += (dot - target) * (dot - target);
        }
    }
    det::sqrt(sum)
}

fn orthonormalize_columns(cols: &mut [f64], m: usize) {
    for _pass in 0..2 {
        for k in 0..m {
            for j in 0..k {
                let mut dot = 0.0;
                for i in 0..m {
                    dot += cols[i * m + j] * cols[i * m + k];
                }
                for i in 0..m {
                    cols[i * m + k] -= dot * cols[i * m + j];
                }
            }
            let nrm = det::sqrt((0..m).map(|i| cols[i * m + k].powi(2)).sum::<f64>());
            if nrm > 0.0 {
                for i in 0..m {
                    cols[i * m + k] /= nrm;
                }
            }
        }
    }
}

fn eye_matrix(n: usize) -> Vec<f64> {
    let mut out = vec![0.0; n * n];
    for i in 0..n {
        out[i * n + i] = 1.0;
    }
    out
}

fn scale_matrix(a: &[f64], s: f64) -> Vec<f64> {
    a.iter().map(|&v| v * s).collect()
}

fn mat_mul(a: &[f64], b: &[f64], n: usize) -> Vec<f64> {
    let mut c = vec![0.0; n * n];
    for i in 0..n {
        for k in 0..n {
            let aik = a[i * n + k];
            if aik == 0.0 {
                continue;
            }
            for j in 0..n {
                c[i * n + j] += aik * b[k * n + j];
            }
        }
    }
    c
}

/// Evaluate $(\operatorname{sinc}(\theta), \operatorname{cosc}(\theta), \varphi_{2,\text{diag}}(\theta), \varphi_{2,\text{off}}(\theta))$
/// where:
/// - $\operatorname{sinc}(\theta) = \sin(\theta)/\theta$
/// - $\operatorname{cosc}(\theta) = (1 - \cos(\theta))/\theta$
/// - $\varphi_{2,\text{diag}}(\theta) = (1 - \cos(\theta))/\theta^2$
/// - $\varphi_{2,\text{off}}(\theta) = (\theta - \sin(\theta))/\theta^2$
fn modal_phi_coeffs(theta: f64) -> (f64, f64, f64, f64) {
    if theta.abs() < 1e-4 {
        let t2 = theta * theta;
        let sinc = 1.0 - t2 / 6.0 + (t2 * t2) / 120.0;
        let cosc = 0.5 * theta * (1.0 - t2 / 12.0 + (t2 * t2) / 360.0);
        let p2_diag = 0.5 - t2 / 24.0 + (t2 * t2) / 720.0;
        let p2_off = theta * (1.0 / 6.0 - t2 / 120.0 + (t2 * t2) / 5040.0);
        (sinc, cosc, p2_diag, p2_off)
    } else {
        let s = det::sin(theta);
        let one_minus_c = 2.0 * det::powi(det::sin(0.5 * theta), 2);
        let sinc = s / theta;
        let cosc = one_minus_c / theta;
        let p2_diag = one_minus_c / (theta * theta);
        let p2_off = (theta - s) / (theta * theta);
        (sinc, cosc, p2_diag, p2_off)
    }
}

#[derive(Debug, Clone)]
enum PropagatorKernel {
    SkewModal {
        sigmas: Vec<f64>,
        u_cols: Vec<f64>,
        v_cols: Vec<f64>,
        modal_ports: Vec<(PortId, Vec<f64>, Vec<f64>)>,
    },
    DenseSeries {
        phi0: Vec<f64>,
        dt_phi1: Vec<f64>,
        port_cols_phi1: Vec<Vec<f64>>,
        port_cols_phi2: Vec<Vec<f64>>,
        ports: Vec<(PortId, Vec<f64>)>,
    },
}

/// Time-domain state of a [`DgDuct`]. The energy coordinates `tilde` are current after
/// every step; the nodal `(p, q)` arrays follow on [`DgPropagator::sync_nodal`].
#[derive(Debug, Clone, PartialEq)]
pub struct DgState {
    operator_identity: ContentHash,
    p: Vec<f64>,
    q: Vec<f64>,
    tilde: Vec<f64>,
    modal_anchor_p: Vec<f64>,
    modal_anchor_q: Vec<f64>,
    modal_anchor_dt_bits: Option<u64>,
    unforced_steps: u64,
    step_count: u64,
    initial_energy_j: f64,
    cumulative_port_work_j: f64,
    cumulative_parameter_work_j: f64,
    cumulative_dissipated_j: f64,
    nodal_dirty: bool,
}

impl DgState {
    /// Physical nodal pressures $p_i^{(e)}$ [Pa] across all elements (current after
    /// [`DgPropagator::sync_nodal`]).
    #[must_use]
    pub fn pressure(&self) -> &[f64] {
        &self.p
    }

    /// Physical nodal volume velocities $q_i^{(e)}$ [$\text{m}^3/\text{s}$] across all
    /// elements (current after [`DgPropagator::sync_nodal`]).
    #[must_use]
    pub fn volume_flow(&self) -> &[f64] {
        &self.q
    }

    /// Stored acoustic and wall-memory energy $E(t) = \tfrac12 \|\tilde x\|^2$ [J].
    #[must_use]
    pub fn stored_energy_j(&self) -> f64 {
        0.5 * self.tilde.iter().map(|v| v * v).sum::<f64>()
    }

    /// Initial stored acoustic energy $E(0)$ [J].
    #[must_use]
    pub const fn initial_energy_j(&self) -> f64 {
        self.initial_energy_j
    }

    /// Signed external parameter work accepted without advancing the clock [J].
    /// Material updates do not reset the initial energy or reclassify work as heat.
    #[must_use]
    pub const fn cumulative_parameter_work_j(&self) -> f64 {
        self.cumulative_parameter_work_j
    }

    /// Original accumulated boundary/wall loss [J]. For the dense kernel this
    /// is inferred from the energy deficit, as declared by its step ledger.
    #[must_use]
    pub const fn cumulative_dissipated_j(&self) -> f64 {
        self.cumulative_dissipated_j
    }

    /// Total steps executed on this state.
    #[must_use]
    pub const fn step_count(&self) -> u64 {
        self.step_count
    }
}

/// Certified exact matrix-exponential propagator for a [`DgDuct`].
#[derive(Debug, Clone)]
pub struct DgPropagator {
    operator_identity: ContentHash,
    dt: f64,
    tol: f64,
    certificate_residual: f64,
    field_dofs: usize,
    total_dofs: usize,
    elements: Vec<DgElement>,
    order: usize,
    kernel: PropagatorKernel,
}

/// Held efforts per declared port, in declaration order.
fn resolve_efforts(port_inputs: &[PortInput], ids: &[PortId]) -> Result<Vec<f64>, DuctError> {
    let mut efforts = vec![0.0; ids.len()];
    for pin in port_inputs {
        if !pin.effort.is_finite() {
            return Err(DuctError::BadParameter {
                what: "port input effort must be finite",
            });
        }
        let Some(slot) = ids.iter().position(|id| *id == pin.id) else {
            return Err(DuctError::BadParameter {
                what: "port input references an undeclared PortId",
            });
        };
        efforts[slot] += pin.effort;
    }
    Ok(efforts)
}

impl DgPropagator {
    /// Time step $\Delta t$ [s].
    #[must_use]
    pub const fn dt(&self) -> f64 {
        self.dt
    }

    /// Admitted certificate tolerance ceiling.
    #[must_use]
    pub const fn tolerance(&self) -> f64 {
        self.tol
    }

    /// The certificate this propagator passed: the larger of the relative singular residual
    /// and the basis orthogonality defects (lossless), or of the carried truncation bound and
    /// the halving estimate (matched).
    #[must_use]
    pub const fn certificate_residual(&self) -> f64 {
        self.certificate_residual
    }

    /// Synchronize the nodal arrays `(state.pressure(), state.volume_flow())` from the
    /// current energy coordinates.
    ///
    /// # Errors
    /// Refuses a state from a different spatial operator without changing it.
    pub fn sync_nodal(&self, state: &mut DgState) -> Result<(), DuctError> {
        self.validate_state(state)?;
        if !state.nodal_dirty {
            return Ok(());
        }
        let m = self.field_dofs;
        let n_local = self.order + 1;
        for (e, elem) in self.elements.iter().enumerate() {
            let base = e * n_local;
            for i in 0..n_local {
                let mut sp = 0.0;
                let mut sq = 0.0;
                for j in 0..n_local {
                    sp += elem.mp_inv_sqrt[i * n_local + j] * state.tilde[base + j];
                    sq += elem.mq_inv_sqrt[i * n_local + j] * state.tilde[m + base + j];
                }
                state.p[base + i] = sp;
                state.q[base + i] = sq;
            }
        }
        state.nodal_dirty = false;
        Ok(())
    }

    fn validate_state(&self, state: &DgState) -> Result<(), DuctError> {
        if state.operator_identity != self.operator_identity
            || state.tilde.len() != self.total_dofs
        {
            return Err(DuctError::BadParameter {
                what: "DG state does not belong to this propagator",
            });
        }
        Ok(())
    }

    /// Advance `state` by one step $\Delta t$ under held `port_inputs`, returning conjugate port
    /// outputs and the step's [`EnergyLedgerRow`].
    ///
    /// The lossless kernel keeps modal amplitudes anchored at the last forced step and
    /// evaluates the exact rotation at the elapsed time, so unforced steps accumulate no
    /// rounding; the ledger's energy is always that of the reconstructed physical state.
    /// Switching to another time step for the same spatial operator reanchors the
    /// modal cache at the current state, preserving the elapsed physical evolution.
    ///
    /// # Errors
    /// [`DuctError::BadParameter`] if any `PortInput` has a non-finite effort, references an
    /// unknown [`PortId`], or the state does not belong to this propagator.
    pub fn step(
        &self,
        state: &mut DgState,
        port_inputs: &[PortInput],
    ) -> Result<DgStepOutput, DuctError> {
        self.validate_state(state)?;
        match &self.kernel {
            PropagatorKernel::SkewModal {
                sigmas,
                u_cols,
                v_cols,
                modal_ports,
            } => self.step_skew(state, port_inputs, sigmas, u_cols, v_cols, modal_ports),
            PropagatorKernel::DenseSeries {
                phi0,
                dt_phi1,
                port_cols_phi1,
                port_cols_phi2,
                ports,
            } => self.step_series(
                state,
                port_inputs,
                phi0,
                dt_phi1,
                port_cols_phi1,
                port_cols_phi2,
                ports,
            ),
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn step_skew(
        &self,
        state: &mut DgState,
        port_inputs: &[PortInput],
        sigmas: &[f64],
        u_cols: &[f64],
        v_cols: &[f64],
        modal_ports: &[(PortId, Vec<f64>, Vec<f64>)],
    ) -> Result<DgStepOutput, DuctError> {
        let m = self.field_dofs;
        let ids: Vec<PortId> = modal_ports.iter().map(|p| p.0).collect();
        let efforts = resolve_efforts(port_inputs, &ids)?;
        if state.modal_anchor_p.is_empty() || state.modal_anchor_dt_bits != Some(self.dt.to_bits())
        {
            let mut ap = vec![0.0; m];
            let mut aq = vec![0.0; m];
            for k in 0..m {
                for i in 0..m {
                    ap[k] += u_cols[i * m + k] * state.tilde[i];
                    aq[k] += v_cols[i * m + k] * state.tilde[m + i];
                }
            }
            state.modal_anchor_p = ap;
            state.modal_anchor_q = aq;
            state.modal_anchor_dt_bits = Some(self.dt.to_bits());
            state.unforced_steps = 0;
        }
        let mut port_outputs = Vec::with_capacity(modal_ports.len());
        let mut step_work = 0.0;
        if modal_ports.is_empty() {
            state.unforced_steps += 1;
        } else {
            // Materialize the anchor at the start of the step, then one exact forced step.
            let t_prev = (state.unforced_steps as f64) * self.dt;
            for k in 0..m {
                let (s, c) = det::sin_cos(sigmas[k] * t_prev);
                let (ap, aq) = (state.modal_anchor_p[k], state.modal_anchor_q[k]);
                state.modal_anchor_p[k] = c * ap - s * aq;
                state.modal_anchor_q[k] = s * ap + c * aq;
            }
            state.unforced_steps = 0;
            let mut int_p = vec![0.0; m];
            let mut int_q = vec![0.0; m];
            for k in 0..m {
                let mut bk_p = 0.0;
                let mut bk_q = 0.0;
                for (slot, (_, bp, bq)) in modal_ports.iter().enumerate() {
                    bk_p += bp[k] * efforts[slot];
                    bk_q += bq[k] * efforts[slot];
                }
                let theta = sigmas[k] * self.dt;
                let (s, c) = det::sin_cos(theta);
                let (sinc, cosc, p2_diag, p2_off) = modal_phi_coeffs(theta);
                let p0 = state.modal_anchor_p[k];
                let q0 = state.modal_anchor_q[k];
                // ∫_0^dt of the modal state: dt phi_1(theta) z_0 + dt^2 phi_2(theta) b_k.
                int_p[k] = self.dt * (sinc * p0 - cosc * q0)
                    + self.dt * self.dt * (p2_diag * bk_p - p2_off * bk_q);
                int_q[k] = self.dt * (cosc * p0 + sinc * q0)
                    + self.dt * self.dt * (p2_off * bk_p + p2_diag * bk_q);
                // End of step: phi_0(theta) z_0 + dt phi_1(theta) b_k.
                state.modal_anchor_p[k] = (c * p0 - s * q0) + self.dt * (sinc * bk_p - cosc * bk_q);
                state.modal_anchor_q[k] = (s * p0 + c * q0) + self.dt * (cosc * bk_p + sinc * bk_q);
            }
            for (slot, (port_id, bp, bq)) in modal_ports.iter().enumerate() {
                let mut delta_x = 0.0;
                let mut end_conj = 0.0;
                for k in 0..m {
                    delta_x += bp[k] * int_p[k] + bq[k] * int_q[k];
                    end_conj += bp[k] * state.modal_anchor_p[k] + bq[k] * state.modal_anchor_q[k];
                }
                step_work += efforts[slot] * delta_x;
                port_outputs.push(PortOutput {
                    id: *port_id,
                    conjugate_rate: delta_x / self.dt,
                    delta_x,
                    end_conjugate: end_conj,
                });
            }
        }
        // The physical state at the end of the step, from the anchor rotated to now.
        let t_now = (state.unforced_steps as f64) * self.dt;
        let mut tilde = vec![0.0; 2 * m];
        for k in 0..m {
            let (s, c) = det::sin_cos(sigmas[k] * t_now);
            let (ap, aq) = (state.modal_anchor_p[k], state.modal_anchor_q[k]);
            let (a, b) = (c * ap - s * aq, s * ap + c * aq);
            for i in 0..m {
                tilde[i] += u_cols[i * m + k] * a;
                tilde[m + i] += v_cols[i * m + k] * b;
            }
        }
        state.tilde = tilde;
        state.nodal_dirty = true;
        state.step_count += 1;
        state.cumulative_port_work_j += step_work;
        let energy_now = state.stored_energy_j();
        let residual = energy_now - state.initial_energy_j - state.cumulative_port_work_j
            - state.cumulative_parameter_work_j;
        let scale = (state.initial_energy_j
            + state.cumulative_port_work_j.abs()
            + state.cumulative_parameter_work_j.abs())
        .max(1e-300);
        Ok(DgStepOutput {
            port_outputs,
            ledger: EnergyLedgerRow {
                step: state.step_count,
                stored_energy_j: energy_now,
                port_work_step_j: step_work,
                cumulative_port_work_j: state.cumulative_port_work_j,
                cumulative_parameter_work_j: state.cumulative_parameter_work_j,
                cumulative_dissipated_j: 0.0,
                energy_balance_residual_j: residual,
                relative_energy_drift: residual.abs() / scale,
                balance_checked: true,
            },
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn step_series(
        &self,
        state: &mut DgState,
        port_inputs: &[PortInput],
        phi0: &[f64],
        dt_phi1: &[f64],
        port_cols_phi1: &[Vec<f64>],
        port_cols_phi2: &[Vec<f64>],
        ports: &[(PortId, Vec<f64>)],
    ) -> Result<DgStepOutput, DuctError> {
        let n = self.total_dofs;
        let ids: Vec<PortId> = ports.iter().map(|p| p.0).collect();
        let efforts = resolve_efforts(port_inputs, &ids)?;
        let e_prev = state.stored_energy_j();
        let mut next_tilde = vec![0.0; n];
        let mut int_tilde = vec![0.0; n];
        for i in 0..n {
            let mut acc_next = 0.0;
            let mut acc_int = 0.0;
            for j in 0..n {
                acc_next += phi0[i * n + j] * state.tilde[j];
                acc_int += dt_phi1[i * n + j] * state.tilde[j];
            }
            for (slot, &u) in efforts.iter().enumerate() {
                acc_next += port_cols_phi1[slot][i] * u;
                acc_int += port_cols_phi2[slot][i] * u;
            }
            next_tilde[i] = acc_next;
            int_tilde[i] = acc_int;
        }
        let mut port_outputs = Vec::with_capacity(ports.len());
        let mut step_work = 0.0;
        for (slot, (port_id, tilde_b)) in ports.iter().enumerate() {
            let mut delta_x = 0.0;
            let mut end_conj = 0.0;
            for i in 0..n {
                delta_x += tilde_b[i] * int_tilde[i];
                end_conj += tilde_b[i] * next_tilde[i];
            }
            step_work += efforts[slot] * delta_x;
            port_outputs.push(PortOutput {
                id: *port_id,
                conjugate_rate: delta_x / self.dt,
                delta_x,
                end_conjugate: end_conj,
            });
        }
        state.tilde = next_tilde;
        state.nodal_dirty = true;
        state.step_count += 1;
        state.cumulative_port_work_j += step_work;
        let e_now = state.stored_energy_j();
        // Absorbed energy in boundaries and wall memory is inferred as the deficit, unclamped: a
        // negative value would expose a non-passive step.
        state.cumulative_dissipated_j += e_prev + step_work - e_now;
        let residual = e_now - state.initial_energy_j - state.cumulative_port_work_j
            - state.cumulative_parameter_work_j
            + state.cumulative_dissipated_j;
        let scale = (state.initial_energy_j
            + state.cumulative_port_work_j.abs()
            + state.cumulative_parameter_work_j.abs())
        .max(1e-300);
        Ok(DgStepOutput {
            port_outputs,
            ledger: EnergyLedgerRow {
                step: state.step_count,
                stored_energy_j: e_now,
                port_work_step_j: step_work,
                cumulative_port_work_j: state.cumulative_port_work_j,
                cumulative_parameter_work_j: state.cumulative_parameter_work_j,
                cumulative_dissipated_j: state.cumulative_dissipated_j,
                energy_balance_residual_j: residual,
                relative_energy_drift: residual.abs() / scale,
                balance_checked: false,
            },
        })
    }
}

#[cfg(test)]
mod tests;
