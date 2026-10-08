//! Planar nonlinear frame analysis with DISTRIBUTED PLASTICITY (plan §15.2
//! step 3: "fiber-section beam-columns … nonlinear time history").
//!
//! ## Element: force-based beam-column (Spacone–Filippou–Taucer 1996)
//!
//! In the corotational BASIC system an element carries basic forces
//! `q = [N, Mᵢ, Mⱼ]` and basic deformations `v = [Δ, θᵢ, θⱼ]`. Equilibrium is
//! interpolated EXACTLY: the section forces at `ξ = x/L` are
//! `D(ξ) = b(ξ) q = [N, (ξ − 1)Mᵢ + ξMⱼ]`, so there is no displacement
//! shape-function error and one element per member resolves spreading
//! plasticity. Sections (fiber or elastic) sit at Gauss–Lobatto points (a
//! point at each end, where hinges form). State determination is the
//! iterative scheme of Neuenhofer & Filippou (J. Struct. Eng. 1997): from
//! the committed state, alternate the element compatibility correction
//! `Δq = F⁻¹(v − ∫bᵀd)` with section Newton corrections
//! `Δd = f_s (b q − D_r(d))` until compatibility AND section equilibrium
//! hold; the element tangent is `K_b = F⁻¹`,
//! `F = ∫ bᵀ f_s b dx`. For an elastic prismatic member this reproduces the
//! exact Euler–Bernoulli flexibility in one pass.
//!
//! ## Geometry
//!
//! - [`Geometry::Linear`]: small displacements (`v = B₀ u`).
//! - [`Geometry::Corotational`]: exact large rigid-body rotations of the
//!   chord (Crisfield), with the consistent geometric stiffness
//!   `N zzᵀ/Lₙ + (Mᵢ + Mⱼ)(rzᵀ + zrᵀ)/Lₙ²`; chord angles are unwrapped against
//!   the committed angle so members may rotate past ±π.
//!
//! ## Analyses
//!
//! Load-controlled and displacement-controlled (pushover) static Newton,
//! modal periods (massless DOFs statically condensed), and implicit
//! Newmark-β dynamics under uniform base acceleration with Rayleigh damping
//! and an energy ledger (input, kinetic, damping, internal work).
//!
//! Determinism: fixed traversal order, dense LU with partial pivoting, no
//! threads — bit-identical for the same input on the same build/ISA.
//!
//! No-claim: planar Euler–Bernoulli members (no shear deformation, no
//! warping/torsion, no out-of-plane buckling); lumped translational mass;
//! no element-level shear failure, bond-slip, or joint panel flexibility;
//! the corotational strains are small (large rotation, small strain).

// Dense structural kernels index `[row][col]` by design.
#![allow(clippy::needless_range_loop)]

use crate::SolidError;
use crate::fiber::{Section, SectionState};

/// A section constitutive model for frame elements.
#[derive(Debug, Clone)]
pub enum SectionModel {
    /// Linear elastic `N = EA·ε₀`, `M = EI·κ`.
    Elastic {
        /// Axial rigidity.
        ea: f64,
        /// Flexural rigidity.
        ei: f64,
    },
    /// A fiber section (fs-material uniaxial laws per fiber).
    Fiber(Section),
}

impl SectionModel {
    /// Response at `(ε₀, κ)` from the committed state (pure).
    #[must_use]
    pub fn respond(&self, eps0: f64, kappa: f64) -> SectionState {
        match self {
            SectionModel::Elastic { ea, ei } => SectionState {
                n: ea * eps0,
                m: ei * kappa,
                tangent: [[*ea, 0.0], [0.0, *ei]],
            },
            SectionModel::Fiber(s) => s.respond(eps0, kappa),
        }
    }

    /// Commit the state at `(ε₀, κ)`.
    pub fn commit(&mut self, eps0: f64, kappa: f64) {
        if let SectionModel::Fiber(s) = self {
            s.commit(eps0, kappa);
        }
    }
}

/// Geometric transformation of the element chord.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Geometry {
    /// Small-displacement theory.
    Linear,
    /// Exact chord rotation (large displacement, small strain).
    Corotational,
}

/// Gauss–Lobatto abscissae on `[0, 1]` and weights (sum 1).
fn lobatto(np: usize) -> Option<(Vec<f64>, Vec<f64>)> {
    let (t, w): (Vec<f64>, Vec<f64>) = match np {
        3 => (vec![-1.0, 0.0, 1.0], vec![1.0 / 3.0, 4.0 / 3.0, 1.0 / 3.0]),
        4 => {
            let a = 1.0 / 5.0f64.sqrt();
            (
                vec![-1.0, -a, a, 1.0],
                vec![1.0 / 6.0, 5.0 / 6.0, 5.0 / 6.0, 1.0 / 6.0],
            )
        }
        5 => {
            let a = (3.0f64 / 7.0).sqrt();
            (
                vec![-1.0, -a, 0.0, a, 1.0],
                vec![0.1, 49.0 / 90.0, 32.0 / 45.0, 49.0 / 90.0, 0.1],
            )
        }
        6 => {
            let (a, b) = (0.765_055_323_929_464_7, 0.285_231_516_480_645_1);
            let (wa, wb) = (0.378_474_956_297_847, 0.554_858_377_035_486_3);
            (
                vec![-1.0, -a, -b, b, a, 1.0],
                vec![1.0 / 15.0, wa, wb, wb, wa, 1.0 / 15.0],
            )
        }
        7 => {
            let (a, b) = (0.830_223_896_278_567, 0.468_848_793_470_714_2);
            let (wa, wb) = (0.276_826_047_361_565_9, 0.431_745_381_209_862_7);
            (
                vec![-1.0, -a, -b, 0.0, b, a, 1.0],
                vec![1.0 / 21.0, wa, wb, 256.0 / 525.0, wb, wa, 1.0 / 21.0],
            )
        }
        _ => return None,
    };
    Some((
        t.iter().map(|x| 0.5 * (x + 1.0)).collect(),
        w.iter().map(|x| 0.5 * x).collect(),
    ))
}

/// Converged element state: basic forces, basic tangent, section
/// deformations.
type BasicState = ([f64; 3], [[f64; 3]; 3], Vec<(f64, f64)>);

/// Global internal force, dense tangent, and per-element trials.
type Assembled = (Vec<f64>, Vec<f64>, Vec<ElementTrial>);

/// Committed element state.
#[derive(Debug, Clone)]
struct ElementCommit {
    q: [f64; 3],
    d: Vec<(f64, f64)>,
    alpha: f64,
}

/// A force-based beam-column between two nodes.
#[derive(Debug, Clone)]
pub struct FrameElement {
    /// Start node.
    pub i: usize,
    /// End node.
    pub j: usize,
    /// Mass per unit length (lumped half to each end).
    pub mass_per_length: f64,
    sections: Vec<SectionModel>,
    xi: Vec<f64>,
    w: Vec<f64>,
    l0: f64,
    cos0: f64,
    sin0: f64,
    commit: ElementCommit,
}

/// One element's trial state at a global configuration.
#[derive(Debug, Clone)]
struct ElementTrial {
    q: [f64; 3],
    kb: [[f64; 3]; 3],
    d: Vec<(f64, f64)>,
    alpha: f64,
    /// Global internal force (6) and tangent (6×6).
    p: [f64; 6],
    k: [[f64; 6]; 6],
}

fn inv2(t: [[f64; 2]; 2]) -> [[f64; 2]; 2] {
    let scale = t[0][0].abs() + t[1][1].abs();
    let mut det = t[0][0] * t[1][1] - t[0][1] * t[1][0];
    let floor = 1e-14 * scale * scale;
    if det.abs() < floor.max(f64::MIN_POSITIVE) {
        det = floor
            .max(f64::MIN_POSITIVE)
            .copysign(if det == 0.0 { 1.0 } else { det });
    }
    [
        [t[1][1] / det, -t[0][1] / det],
        [-t[1][0] / det, t[0][0] / det],
    ]
}

fn inv3(m: [[f64; 3]; 3]) -> Option<[[f64; 3]; 3]> {
    let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    if !det.is_finite() || det == 0.0 {
        return None;
    }
    let c = |a: usize, b: usize, c2: usize, d: usize| m[a][b] * m[c2][d] - m[a][d] * m[c2][b];
    Some([
        [
            c(1, 1, 2, 2) / det,
            -c(0, 1, 2, 2) / det,
            c(0, 1, 1, 2) / det,
        ],
        [
            -c(1, 0, 2, 2) / det,
            c(0, 0, 2, 2) / det,
            -c(0, 0, 1, 2) / det,
        ],
        [
            c(1, 0, 2, 1) / det,
            -c(0, 0, 2, 1) / det,
            c(0, 0, 1, 1) / det,
        ],
    ])
}

impl FrameElement {
    /// Basic flexibility interpolation `b(ξ)`: rows (N, M), columns (N, Mᵢ, Mⱼ).
    fn b(xi: f64) -> [[f64; 3]; 2] {
        [[1.0, 0.0, 0.0], [0.0, xi - 1.0, xi]]
    }

    /// Element state determination for trial basic deformations `v`: the
    /// direct iteration from the committed state, falling back to
    /// sub-incrementing `v_c → v` in 2, 4, … 64 substeps (each substep
    /// starting from the previous substep's converged trial state) when the
    /// direct iteration fails — the standard remedy at fiber load reversals.
    fn basic_state(&self, v: [f64; 3]) -> Result<BasicState, SolidError> {
        let first = match self.iterate(self.commit.q, self.commit.d.clone(), v) {
            Ok(s) => return Ok(s),
            Err(e) => e,
        };
        let vc = self.committed_deformation();
        let mut parts = 2usize;
        while parts <= 64 {
            let mut state: Option<BasicState> = None;
            let mut ok = true;
            for k in 1..=parts {
                let t = k as f64 / parts as f64;
                let vk = [
                    vc[0] + t * (v[0] - vc[0]),
                    vc[1] + t * (v[1] - vc[1]),
                    vc[2] + t * (v[2] - vc[2]),
                ];
                let (q0, d0) = match &state {
                    Some((q, _, d)) => (*q, d.clone()),
                    None => (self.commit.q, self.commit.d.clone()),
                };
                if let Ok(s) = self.iterate(q0, d0, vk) {
                    state = Some(s);
                } else {
                    ok = false;
                    break;
                }
            }
            if let (true, Some(s)) = (ok, state) {
                return Ok(s);
            }
            parts *= 2;
        }
        Err(first)
    }

    /// Basic deformations implied by the committed section deformations.
    fn committed_deformation(&self) -> [f64; 3] {
        let mut v = [0.0; 3];
        for (p, &(e0, k0)) in self.commit.d.iter().enumerate() {
            let b = Self::b(self.xi[p]);
            let wl = self.w[p] * self.l0;
            for a in 0..3 {
                v[a] += wl * (b[0][a] * e0 + b[1][a] * k0);
            }
        }
        v
    }

    /// Neuenhofer–Filippou iteration from a given element state.
    fn iterate(
        &self,
        q0: [f64; 3],
        d0: Vec<(f64, f64)>,
        v: [f64; 3],
    ) -> Result<BasicState, SolidError> {
        let l = self.l0;
        let np = self.xi.len();
        let mut q = q0;
        let mut d = d0;
        let vnorm = v.iter().fold(0.0f64, |m, x| m.max(x.abs()));
        let mut history = Vec::new();
        for _ in 0..100 {
            // Section responses, flexibilities, element flexibility and v_r.
            let mut f = [[0.0; 3]; 3];
            let mut vr = [0.0; 3];
            let mut resp = Vec::with_capacity(np);
            for p in 0..np {
                let s = self.sections[p].respond(d[p].0, d[p].1);
                let fs = inv2(s.tangent);
                let b = Self::b(self.xi[p]);
                let wl = self.w[p] * l;
                for a in 0..3 {
                    vr[a] += wl * (b[0][a] * d[p].0 + b[1][a] * d[p].1);
                    for c in 0..3 {
                        let mut acc = 0.0;
                        for r in 0..2 {
                            for s2 in 0..2 {
                                acc += b[r][a] * fs[r][s2] * b[s2][c];
                            }
                        }
                        f[a][c] += wl * acc;
                    }
                }
                resp.push((s, fs));
            }
            let kb = inv3(f).ok_or_else(|| SolidError::InternalInvariant {
                what: "singular element flexibility".into(),
            })?;
            // Section unbalance at the current q.
            let mut unbal = 0.0f64;
            let mut dscale = 1.0f64;
            for p in 0..np {
                let b = Self::b(self.xi[p]);
                let dn = b[0][0] * q[0];
                let dm = b[1][1] * q[1] + b[1][2] * q[2];
                let (s, _) = resp[p];
                unbal = unbal.max((dn - s.n).abs()).max((dm - s.m).abs());
                dscale = dscale.max(dn.abs()).max(dm.abs());
            }
            let dv: Vec<f64> = (0..3).map(|a| v[a] - vr[a]).collect();
            let dvn = dv.iter().fold(0.0f64, |m, x| m.max(x.abs()));
            history.push(dvn);
            if dvn <= 1e-15 + 1e-13 * vnorm && unbal <= 1e-12 * dscale {
                return Ok((q, kb, d));
            }
            // Compatibility correction, then section corrections.
            for a in 0..3 {
                q[a] += (0..3).map(|c| kb[a][c] * dv[c]).sum::<f64>();
            }
            for p in 0..np {
                let b = Self::b(self.xi[p]);
                let dn = b[0][0] * q[0];
                let dm = b[1][1] * q[1] + b[1][2] * q[2];
                let (s, fs) = resp[p];
                let (rn, rm) = (dn - s.n, dm - s.m);
                d[p].0 += fs[0][0] * rn + fs[0][1] * rm;
                d[p].1 += fs[1][0] * rn + fs[1][1] * rm;
            }
        }
        Err(SolidError::NewtonStalled { history })
    }
}

/// A planar frame model.
#[derive(Debug, Clone)]
pub struct Frame2d {
    nodes: Vec<[f64; 2]>,
    elements: Vec<FrameElement>,
    fixed: Vec<bool>,
    nodal_mass: Vec<f64>,
    geometry: Geometry,
    u: Vec<f64>,
    vel: Vec<f64>,
    acc: Vec<f64>,
}

/// One static step record.
#[derive(Debug, Clone, PartialEq)]
pub struct StaticStep {
    /// Load factor reached.
    pub lambda: f64,
    /// Displacements (3 per node: ux, uy, θ).
    pub u: Vec<f64>,
    /// Newton iterations used.
    pub iterations: usize,
}

/// Rayleigh damping `C = a₀M + a₁K₀` (initial stiffness).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Rayleigh {
    /// Mass-proportional coefficient.
    pub a0: f64,
    /// Initial-stiffness-proportional coefficient.
    pub a1: f64,
}

impl Rayleigh {
    /// Coefficients giving damping ratio `zeta` at circular frequencies
    /// `w1` and `w2`.
    #[must_use]
    pub fn from_modes(zeta: f64, w1: f64, w2: f64) -> Rayleigh {
        Rayleigh {
            a0: zeta * 2.0 * w1 * w2 / (w1 + w2),
            a1: zeta * 2.0 / (w1 + w2),
        }
    }
}

/// A dynamic analysis record.
#[derive(Debug, Clone, PartialEq)]
pub struct DynamicHistory {
    /// Times.
    pub t: Vec<f64>,
    /// Displacements per step (3 per node).
    pub u: Vec<Vec<f64>>,
    /// Base shear (sum of horizontal reactions) per step.
    pub base_shear: Vec<f64>,
    /// Cumulative input energy (relative formulation).
    pub input_energy: Vec<f64>,
    /// Kinetic energy.
    pub kinetic_energy: Vec<f64>,
    /// Cumulative damping dissipation.
    pub damping_energy: Vec<f64>,
    /// Cumulative internal work (recoverable strain + hysteretic).
    pub internal_work: Vec<f64>,
    /// Newton iterations per step.
    pub iterations: Vec<usize>,
}

impl DynamicHistory {
    /// Largest relative energy-balance error `|E_in − E_k − E_d − E_int|`
    /// over the record, normalized by the peak input energy.
    #[must_use]
    pub fn energy_balance_error(&self) -> f64 {
        let peak = self
            .input_energy
            .iter()
            .fold(0.0f64, |m, e| m.max(e.abs()))
            .max(f64::MIN_POSITIVE);
        (0..self.t.len())
            .map(|k| {
                (self.input_energy[k]
                    - self.kinetic_energy[k]
                    - self.damping_energy[k]
                    - self.internal_work[k])
                    .abs()
            })
            .fold(0.0, f64::max)
            / peak
    }
}

/// Newton acceptance: the force residual meets `tol`, or the iteration has
/// reached the round-off floor — a correction negligible against the
/// displacement scale while the residual is within `floor_tol`.
fn newton_done(rn: f64, tol: f64, floor_tol: f64, du_norm: f64, u_norm: f64) -> bool {
    rn <= tol || (rn <= floor_tol && du_norm <= 1e-13 * u_norm.max(1e-30))
}

fn inf_norm(v: &[f64]) -> f64 {
    v.iter().fold(0.0f64, |m, x| m.max(x.abs()))
}

/// Dense LU solve (partial pivoting). `None` if singular.
fn lu_solve(a: &[f64], n: usize, b: &[f64]) -> Option<Vec<f64>> {
    let mut m = a.to_vec();
    let mut x = b.to_vec();
    let scale = m.iter().fold(0.0f64, |s, v| s.max(v.abs()));
    if !(scale > 0.0) {
        return None;
    }
    for col in 0..n {
        let mut piv = col;
        for r in (col + 1)..n {
            if m[r * n + col].abs() > m[piv * n + col].abs() {
                piv = r;
            }
        }
        if m[piv * n + col].abs() <= 1e-15 * scale {
            return None;
        }
        if piv != col {
            for k in 0..n {
                m.swap(col * n + k, piv * n + k);
            }
            x.swap(col, piv);
        }
        let dg = m[col * n + col];
        for r in (col + 1)..n {
            let f = m[r * n + col] / dg;
            if f == 0.0 {
                continue;
            }
            for k in col..n {
                m[r * n + k] -= f * m[col * n + k];
            }
            x[r] -= f * x[col];
        }
    }
    for r in (0..n).rev() {
        let mut s = x[r];
        for k in (r + 1)..n {
            s -= m[r * n + k] * x[k];
        }
        x[r] = s / m[r * n + r];
    }
    x.iter().all(|v| v.is_finite()).then_some(x)
}

impl Frame2d {
    /// An empty frame with the given geometric transformation.
    #[must_use]
    pub fn new(geometry: Geometry) -> Frame2d {
        Frame2d {
            nodes: Vec::new(),
            elements: Vec::new(),
            fixed: Vec::new(),
            nodal_mass: Vec::new(),
            geometry,
            u: Vec::new(),
            vel: Vec::new(),
            acc: Vec::new(),
        }
    }

    /// Add a node; returns its index.
    pub fn add_node(&mut self, x: f64, y: f64) -> usize {
        self.nodes.push([x, y]);
        self.fixed.extend([false; 3]);
        self.nodal_mass.extend([0.0; 3]);
        self.u.extend([0.0; 3]);
        self.vel.extend([0.0; 3]);
        self.acc.extend([0.0; 3]);
        self.nodes.len() - 1
    }

    /// Fix DOFs of a node (`[ux, uy, θ]`).
    pub fn fix(&mut self, node: usize, dofs: [bool; 3]) {
        for k in 0..3 {
            self.fixed[3 * node + k] |= dofs[k];
        }
    }

    /// Add a lumped mass to a node's translational DOFs.
    pub fn add_mass(&mut self, node: usize, mass: f64) {
        self.nodal_mass[3 * node] += mass;
        self.nodal_mass[3 * node + 1] += mass;
    }

    /// Add a force-based element with `np` Gauss–Lobatto sections built by
    /// `make_section`; returns its index.
    ///
    /// # Errors
    /// [`SolidError::InvalidInput`] for unknown nodes, coincident nodes, an
    /// unsupported `np` (3..=7), or non-finite mass.
    pub fn add_element(
        &mut self,
        i: usize,
        j: usize,
        np: usize,
        mass_per_length: f64,
        make_section: &dyn Fn() -> SectionModel,
    ) -> Result<usize, SolidError> {
        if i >= self.nodes.len() || j >= self.nodes.len() || i == j {
            return Err(SolidError::InvalidInput {
                what: format!("element nodes ({i}, {j}) invalid"),
            });
        }
        let (xi, w) = lobatto(np).ok_or_else(|| SolidError::InvalidInput {
            what: format!("{np} Lobatto points unsupported (3..=7)"),
        })?;
        if !(mass_per_length.is_finite() && mass_per_length >= 0.0) {
            return Err(SolidError::InvalidInput {
                what: "mass per length must be finite and non-negative".into(),
            });
        }
        let (a, b) = (self.nodes[i], self.nodes[j]);
        let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
        let l0 = dx.hypot(dy);
        if !(l0 > 0.0) {
            return Err(SolidError::InvalidInput {
                what: "zero-length element".into(),
            });
        }
        let half = 0.5 * mass_per_length * l0;
        self.add_mass(i, half);
        self.add_mass(j, half);
        self.elements.push(FrameElement {
            i,
            j,
            mass_per_length,
            sections: (0..np).map(|_| make_section()).collect(),
            commit: ElementCommit {
                q: [0.0; 3],
                d: vec![(0.0, 0.0); np],
                alpha: 0.0,
            },
            xi,
            w,
            l0,
            cos0: dx / l0,
            sin0: dy / l0,
        });
        Ok(self.elements.len() - 1)
    }

    /// Number of DOFs (3 per node).
    #[must_use]
    pub fn ndof(&self) -> usize {
        3 * self.nodes.len()
    }

    /// Committed displacements.
    #[must_use]
    pub fn displacements(&self) -> &[f64] {
        &self.u
    }

    /// Committed basic forces `[N, Mᵢ, Mⱼ]` of an element.
    #[must_use]
    pub fn element_forces(&self, e: usize) -> [f64; 3] {
        self.elements[e].commit.q
    }

    /// Committed section deformations `(ε₀, κ)` of an element.
    #[must_use]
    pub fn section_deformations(&self, e: usize) -> &[(f64, f64)] {
        &self.elements[e].commit.d
    }

    fn element_trial(&self, e: &FrameElement, u: &[f64]) -> Result<ElementTrial, SolidError> {
        let dofs = [
            3 * e.i,
            3 * e.i + 1,
            3 * e.i + 2,
            3 * e.j,
            3 * e.j + 1,
            3 * e.j + 2,
        ];
        let ue: Vec<f64> = dofs.iter().map(|&d| u[d]).collect();
        let (c, s, ln, alpha) = match self.geometry {
            Geometry::Linear => (e.cos0, e.sin0, e.l0, 0.0),
            Geometry::Corotational => {
                let dx = e.l0 * e.cos0 + ue[3] - ue[0];
                let dy = e.l0 * e.sin0 + ue[4] - ue[1];
                let ln = dx.hypot(dy);
                let (c, s) = (dx / ln, dy / ln);
                // Rigid chord rotation relative to the initial chord,
                // unwrapped to the branch nearest the committed angle.
                let raw = (e.cos0 * s - e.sin0 * c).atan2(e.cos0 * c + e.sin0 * s);
                let two_pi = 2.0 * std::f64::consts::PI;
                let k = ((e.commit.alpha - raw) / two_pi).round();
                (c, s, ln, raw + k * two_pi)
            }
        };
        let r = [-c, -s, 0.0, c, s, 0.0];
        let z = [s, -c, 0.0, -s, c, 0.0];
        let mut bmat = [[0.0; 6]; 3];
        bmat[0] = r;
        for k in 0..6 {
            bmat[1][k] = -z[k] / ln;
            bmat[2][k] = -z[k] / ln;
        }
        bmat[1][2] += 1.0;
        bmat[2][5] += 1.0;
        let v = match self.geometry {
            Geometry::Linear => {
                let mut v = [0.0; 3];
                for a in 0..3 {
                    v[a] = (0..6).map(|k| bmat[a][k] * ue[k]).sum();
                }
                v
            }
            Geometry::Corotational => [ln - e.l0, ue[2] - alpha, ue[5] - alpha],
        };
        let (q, kb, d) = e.basic_state(v)?;
        let mut p = [0.0; 6];
        let mut k = [[0.0; 6]; 6];
        for a in 0..6 {
            p[a] = (0..3).map(|m| bmat[m][a] * q[m]).sum();
            for b in 0..6 {
                let mut acc = 0.0;
                for m in 0..3 {
                    for n in 0..3 {
                        acc += bmat[m][a] * kb[m][n] * bmat[n][b];
                    }
                }
                k[a][b] = acc;
            }
        }
        if self.geometry == Geometry::Corotational {
            let mm = (q[1] + q[2]) / (ln * ln);
            for a in 0..6 {
                for b in 0..6 {
                    k[a][b] += q[0] / ln * z[a] * z[b] + mm * (r[a] * z[b] + z[a] * r[b]);
                }
            }
        }
        Ok(ElementTrial {
            q,
            kb,
            d,
            alpha,
            p,
            k,
        })
    }

    /// Internal force vector, tangent (dense row-major) and element trials.
    fn assemble(&self, u: &[f64]) -> Result<Assembled, SolidError> {
        let n = self.ndof();
        let mut f = vec![0.0; n];
        let mut kg = vec![0.0; n * n];
        let mut trials = Vec::with_capacity(self.elements.len());
        for e in &self.elements {
            let t = self.element_trial(e, u)?;
            let dofs = [
                3 * e.i,
                3 * e.i + 1,
                3 * e.i + 2,
                3 * e.j,
                3 * e.j + 1,
                3 * e.j + 2,
            ];
            for a in 0..6 {
                f[dofs[a]] += t.p[a];
                for b in 0..6 {
                    kg[dofs[a] * n + dofs[b]] += t.k[a][b];
                }
            }
            trials.push(t);
        }
        Ok((f, kg, trials))
    }

    fn commit_trials(&mut self, trials: Vec<ElementTrial>) {
        for (e, t) in self.elements.iter_mut().zip(trials) {
            for (p, &(e0, k0)) in t.d.iter().enumerate() {
                e.sections[p].commit(e0, k0);
            }
            e.commit = ElementCommit {
                q: t.q,
                d: t.d,
                alpha: t.alpha,
            };
            let _ = t.kb;
        }
    }

    /// Solve `K Δu = r` on the free DOFs (fixed DOFs get zero).
    fn solve_free(&self, k: &[f64], r: &[f64]) -> Option<Vec<f64>> {
        let n = self.ndof();
        let free: Vec<usize> = (0..n).filter(|&d| !self.fixed[d]).collect();
        let nf = free.len();
        let mut kf = vec![0.0; nf * nf];
        for (a, &da) in free.iter().enumerate() {
            for (b, &db) in free.iter().enumerate() {
                kf[a * nf + b] = k[da * n + db];
            }
        }
        let rf: Vec<f64> = free.iter().map(|&d| r[d]).collect();
        let x = lu_solve(&kf, nf, &rf)?;
        let mut out = vec![0.0; n];
        for (a, &d) in free.iter().enumerate() {
            out[d] = x[a];
        }
        Some(out)
    }

    fn free_norm(&self, r: &[f64]) -> f64 {
        r.iter()
            .enumerate()
            .filter(|(d, _)| !self.fixed[*d])
            .fold(0.0f64, |m, (_, v)| m.max(v.abs()))
    }

    /// Backtracking line search along `du` on the residual norm returned by
    /// `residual`: the first step length in 1, ½, … that decreases the norm
    /// below `rn`, or the smallest tried (2⁻⁹) so Newton can still leave a
    /// kink. Fiber laws are non-smooth at reversals and at the concrete
    /// tension cut-off (zero tangent at zero strain), where full Newton steps
    /// overshoot and cycle.
    fn line_search(
        &self,
        u: &[f64],
        du: &[f64],
        rn: f64,
        residual: impl Fn(&[f64], &Assembled) -> f64,
    ) -> Result<(Vec<f64>, Assembled), SolidError> {
        let mut alpha = 1.0f64;
        let mut last_err = None;
        for _ in 0..10 {
            let ut: Vec<f64> = u.iter().zip(du).map(|(a, b)| a + alpha * b).collect();
            match self.assemble(&ut) {
                Ok(asm) => {
                    if residual(&ut, &asm) < rn || alpha <= 1.0 / 512.0 {
                        return Ok((ut, asm));
                    }
                }
                Err(e) => last_err = Some(e),
            }
            alpha *= 0.5;
        }
        Err(last_err.unwrap_or(SolidError::NewtonStalled { history: vec![rn] }))
    }

    /// Load-controlled static analysis: apply `load` (one entry per DOF) in
    /// `steps` equal increments of the load factor, Newton at each.
    ///
    /// # Errors
    /// [`SolidError::NewtonStalled`] if an increment fails to converge
    /// (refine `steps`); [`SolidError::InvalidInput`] on a wrong-length load.
    pub fn static_load(
        &mut self,
        load: &[f64],
        steps: usize,
    ) -> Result<Vec<StaticStep>, SolidError> {
        let n = self.ndof();
        if load.len() != n || steps == 0 {
            return Err(SolidError::InvalidInput {
                what: "load vector length must equal ndof and steps ≥ 1".into(),
            });
        }
        let pref = load.iter().fold(0.0f64, |m, v| m.max(v.abs())).max(1.0);
        let mut out = Vec::with_capacity(steps);
        for step in 1..=steps {
            let lambda = step as f64 / steps as f64;
            let p: Vec<f64> = load.iter().map(|v| lambda * v).collect();
            let mut u = self.u.clone();
            let mut history = Vec::new();
            let mut done = None;
            let (mut f, mut k, mut trials) = self.assemble(&u)?;
            for it in 0..60 {
                let r: Vec<f64> = (0..n).map(|d| p[d] - f[d]).collect();
                let rn = self.free_norm(&r);
                history.push(rn);
                if rn <= 1e-9 * pref {
                    done = Some((it, trials));
                    break;
                }
                let du = self
                    .solve_free(&k, &r)
                    .ok_or_else(|| SolidError::NewtonStalled {
                        history: history.clone(),
                    })?;
                if newton_done(rn, 1e-9 * pref, 1e-6 * pref, inf_norm(&du), inf_norm(&u)) {
                    done = Some((it, trials));
                    break;
                }
                let residual = |_: &[f64], a: &Assembled| -> f64 {
                    let rr: Vec<f64> = (0..n).map(|d| p[d] - a.0[d]).collect();
                    self.free_norm(&rr)
                };
                let (ut, asm) = self.line_search(&u, &du, rn, residual)?;
                u = ut;
                (f, k, trials) = asm;
            }
            let Some((iterations, trials)) = done else {
                return Err(SolidError::NewtonStalled { history });
            };
            self.u = u;
            self.commit_trials(trials);
            out.push(StaticStep {
                lambda,
                u: self.u.clone(),
                iterations,
            });
        }
        Ok(out)
    }

    /// Displacement-controlled static analysis (pushover): scale the
    /// reference load pattern `pattern` so DOF `control` reaches each of
    /// `targets` in turn (bordered Newton; captures softening branches a
    /// load-controlled run cannot).
    ///
    /// # Errors
    /// [`SolidError`] on bad input or a non-converging increment.
    pub fn displacement_control(
        &mut self,
        pattern: &[f64],
        control: usize,
        targets: &[f64],
    ) -> Result<Vec<StaticStep>, SolidError> {
        let n = self.ndof();
        if pattern.len() != n || control >= n || self.fixed[control] {
            return Err(SolidError::InvalidInput {
                what: "pattern length, or control DOF free and in range".into(),
            });
        }
        let pref = pattern.iter().fold(0.0f64, |m, v| m.max(v.abs())).max(1.0);
        let mut lambda = 0.0f64;
        let mut out = Vec::with_capacity(targets.len());
        for &target in targets {
            let mut u = self.u.clone();
            let mut history = Vec::new();
            let mut done = None;
            for it in 0..60 {
                let (f, k, trials) = self.assemble(&u)?;
                let r: Vec<f64> = (0..n).map(|d| lambda * pattern[d] - f[d]).collect();
                let rn = self.free_norm(&r);
                let gap = target - u[control];
                history.push(rn + gap.abs());
                if rn <= 1e-9 * pref * lambda.abs().max(1.0)
                    && gap.abs() <= 1e-12 * (1.0 + target.abs())
                {
                    done = Some((it, trials));
                    break;
                }
                let ut = self
                    .solve_free(&k, pattern)
                    .ok_or_else(|| SolidError::NewtonStalled {
                        history: history.clone(),
                    })?;
                let ur = self
                    .solve_free(&k, &r)
                    .ok_or_else(|| SolidError::NewtonStalled {
                        history: history.clone(),
                    })?;
                if ut[control].abs() < 1e-300 {
                    return Err(SolidError::NewtonStalled { history });
                }
                let dl = (gap - ur[control]) / ut[control];
                let step: Vec<f64> = (0..n).map(|d| ur[d] + dl * ut[d]).collect();
                let scale = pref * lambda.abs().max(1.0);
                if gap.abs() <= 1e-12 * (1.0 + target.abs())
                    && newton_done(
                        rn,
                        1e-9 * scale,
                        1e-6 * scale,
                        inf_norm(&step),
                        inf_norm(&u),
                    )
                {
                    done = Some((it, trials));
                    break;
                }
                lambda += dl;
                for d in 0..n {
                    u[d] += step[d];
                }
            }
            let Some((iterations, trials)) = done else {
                return Err(SolidError::NewtonStalled { history });
            };
            self.u = u;
            self.commit_trials(trials);
            out.push(StaticStep {
                lambda,
                u: self.u.clone(),
                iterations,
            });
        }
        Ok(out)
    }

    /// Free DOFs carrying mass (the dynamic DOFs) and the massless free ones.
    fn mass_partition(&self) -> (Vec<usize>, Vec<usize>) {
        let n = self.ndof();
        let free: Vec<usize> = (0..n).filter(|&d| !self.fixed[d]).collect();
        let massive = free
            .iter()
            .copied()
            .filter(|&d| self.nodal_mass[d] > 0.0)
            .collect();
        let massless = free
            .iter()
            .copied()
            .filter(|&d| self.nodal_mass[d] <= 0.0)
            .collect();
        (massive, massless)
    }

    /// Circular natural frequencies `ω` (ascending) of the current tangent,
    /// with massless DOFs statically condensed out.
    ///
    /// # Errors
    /// [`SolidError`] if the tangent cannot be formed or condensed.
    pub fn natural_frequencies(&self) -> Result<Vec<f64>, SolidError> {
        let n = self.ndof();
        let (_, k, _) = self.assemble(&self.u)?;
        let (t, r) = self.mass_partition();
        let (nt, nr) = (t.len(), r.len());
        // K_cc = K_tt − K_tr K_rr⁻¹ K_rt.
        let mut krr = vec![0.0; nr * nr];
        for (a, &da) in r.iter().enumerate() {
            for (b, &db) in r.iter().enumerate() {
                krr[a * nr + b] = k[da * n + db];
            }
        }
        let mut kcc = vec![0.0; nt * nt];
        for (a, &da) in t.iter().enumerate() {
            for (b, &db) in t.iter().enumerate() {
                kcc[a * nt + b] = k[da * n + db];
            }
        }
        if nr > 0 {
            for (b, &db) in t.iter().enumerate() {
                let col: Vec<f64> = r.iter().map(|&dr| k[dr * n + db]).collect();
                let x = lu_solve(&krr, nr, &col).ok_or_else(|| SolidError::InternalInvariant {
                    what: "singular massless block in modal condensation".into(),
                })?;
                for (a, &da) in t.iter().enumerate() {
                    let s: f64 = r
                        .iter()
                        .enumerate()
                        .map(|(c, &dr)| k[da * n + dr] * x[c])
                        .sum();
                    kcc[a * nt + b] -= s;
                }
            }
        }
        // M^{-1/2} K_cc M^{-1/2} (M diagonal), symmetric eigenvalues.
        let mut a = vec![vec![0.0; nt]; nt];
        for i in 0..nt {
            for j in 0..nt {
                let mi = self.nodal_mass[t[i]].sqrt();
                let mj = self.nodal_mass[t[j]].sqrt();
                a[i][j] = 0.5 * (kcc[i * nt + j] + kcc[j * nt + i]) / (mi * mj);
            }
        }
        let mut ev = jacobi_eigenvalues(a);
        ev.sort_by(f64::total_cmp);
        Ok(ev.iter().map(|l| l.max(0.0).sqrt()).collect())
    }

    /// Implicit Newmark-β (average acceleration, γ = ½, β = ¼) response to a
    /// uniform base acceleration record `ground` (sampled every `dt`) acting
    /// along the global direction `(dir_x, dir_y)`, with Rayleigh damping on
    /// the initial stiffness. Displacements are RELATIVE to the base.
    ///
    /// # Errors
    /// [`SolidError::NewtonStalled`] if a step fails (reduce `dt`).
    #[allow(clippy::too_many_lines)]
    pub fn newmark(
        &mut self,
        ground: &[f64],
        dt: f64,
        dir: (f64, f64),
        damping: Rayleigh,
    ) -> Result<DynamicHistory, SolidError> {
        let n = self.ndof();
        if !(dt.is_finite() && dt > 0.0) {
            return Err(SolidError::InvalidInput {
                what: "time step must be finite and positive".into(),
            });
        }
        let (gamma, beta) = (0.5, 0.25);
        let m = self.nodal_mass.clone();
        let infl: Vec<f64> = (0..n)
            .map(|d| match d % 3 {
                0 => dir.0,
                1 => dir.1,
                _ => 0.0,
            })
            .collect();
        let (_, k0, _) = self.assemble(&self.u)?;
        let cmat: Vec<f64> = (0..n * n)
            .map(|idx| {
                let (i, j) = (idx / n, idx % n);
                damping.a1 * k0[idx] + if i == j { damping.a0 * m[i] } else { 0.0 }
            })
            .collect();
        let cv = |v: &[f64]| -> Vec<f64> {
            (0..n)
                .map(|i| (0..n).map(|j| cmat[i * n + j] * v[j]).sum())
                .collect()
        };
        let pext = |ag: f64| -> Vec<f64> { (0..n).map(|d| -m[d] * infl[d] * ag).collect() };
        let fref = (0..n)
            .map(|d| (m[d] * infl[d]).abs())
            .fold(0.0f64, f64::max)
            * ground.iter().fold(0.0f64, |a, g| a.max(g.abs()))
            + 1.0;
        let mut hist = DynamicHistory {
            t: vec![0.0],
            u: vec![self.u.clone()],
            base_shear: vec![0.0],
            input_energy: vec![0.0],
            kinetic_energy: vec![
                0.5 * (0..n)
                    .map(|d| m[d] * self.vel[d] * self.vel[d])
                    .sum::<f64>(),
            ],
            damping_energy: vec![0.0],
            internal_work: vec![0.0],
            iterations: vec![0],
        };
        // Consistent initial acceleration on massive DOFs.
        let (f0, _, _) = self.assemble(&self.u)?;
        let p0 = pext(ground.first().copied().unwrap_or(0.0));
        let c0 = cv(&self.vel);
        for d in 0..n {
            self.acc[d] = if m[d] > 0.0 && !self.fixed[d] {
                (p0[d] - c0[d] - f0[d]) / m[d]
            } else {
                0.0
            };
        }
        let mut f_prev = f0;
        let mut p_prev = p0;
        let mut c_prev = c0;
        for step in 1..ground.len() {
            let p1 = pext(ground[step]);
            let (u0, v0, a0) = (self.u.clone(), self.vel.clone(), self.acc.clone());
            let mut u = u0.clone();
            let mut history = Vec::new();
            let mut done = None;
            let kin = |u: &[f64]| -> (Vec<f64>, Vec<f64>) {
                let acc: Vec<f64> = (0..n)
                    .map(|d| {
                        (u[d] - u0[d]) / (beta * dt * dt)
                            - v0[d] / (beta * dt)
                            - (0.5 / beta - 1.0) * a0[d]
                    })
                    .collect();
                let vel: Vec<f64> = (0..n)
                    .map(|d| v0[d] + dt * ((1.0 - gamma) * a0[d] + gamma * acc[d]))
                    .collect();
                (acc, vel)
            };
            let dyn_residual = |u: &[f64], f: &[f64]| -> Vec<f64> {
                let (acc, vel) = kin(u);
                let cvel = cv(&vel);
                (0..n)
                    .map(|d| p1[d] - m[d] * acc[d] - cvel[d] - f[d])
                    .collect()
            };
            let (mut f, mut kt, mut trials) = self.assemble(&u)?;
            for it in 0..60 {
                let r = dyn_residual(&u, &f);
                let rn = self.free_norm(&r);
                history.push(rn);
                let (acc, vel) = kin(&u);
                if rn <= 1e-9 * fref {
                    let cvel = cv(&vel);
                    done = Some((it, trials, f, acc, vel, cvel));
                    break;
                }
                let mut keff = kt;
                let (cm, cc) = (1.0 / (beta * dt * dt), gamma / (beta * dt));
                for i in 0..n {
                    keff[i * n + i] += cm * m[i];
                    for j in 0..n {
                        keff[i * n + j] += cc * cmat[i * n + j];
                    }
                }
                let du = self
                    .solve_free(&keff, &r)
                    .ok_or_else(|| SolidError::NewtonStalled {
                        history: history.clone(),
                    })?;
                if newton_done(rn, 1e-9 * fref, 1e-6 * fref, inf_norm(&du), inf_norm(&u)) {
                    let cvel = cv(&vel);
                    done = Some((it, trials, f, acc, vel, cvel));
                    break;
                }
                let (ut, asm) = self.line_search(&u, &du, rn, |ut: &[f64], a: &Assembled| {
                    self.free_norm(&dyn_residual(ut, &a.0))
                })?;
                u = ut;
                (f, kt, trials) = asm;
            }
            let Some((iters, trials, f, acc, vel, cvel)) = done else {
                return Err(SolidError::NewtonStalled { history });
            };
            // Energy increments (trapezoidal over the step).
            let du: Vec<f64> = (0..n).map(|d| u[d] - u0[d]).collect();
            let e_in = 0.5 * (0..n).map(|d| du[d] * (p_prev[d] + p1[d])).sum::<f64>();
            let e_int = 0.5 * (0..n).map(|d| du[d] * (f_prev[d] + f[d])).sum::<f64>();
            let e_d = 0.5 * (0..n).map(|d| du[d] * (c_prev[d] + cvel[d])).sum::<f64>();
            let kin = 0.5 * (0..n).map(|d| m[d] * vel[d] * vel[d]).sum::<f64>();
            let base: f64 = (0..n)
                .filter(|&d| self.fixed[d] && d % 3 == 0)
                .map(|d| -f[d])
                .sum();
            self.u = u;
            self.vel = vel;
            self.acc = acc;
            self.commit_trials(trials);
            let last = hist.t.len() - 1;
            hist.t.push(step as f64 * dt);
            hist.u.push(self.u.clone());
            hist.base_shear.push(base);
            hist.input_energy.push(hist.input_energy[last] + e_in);
            hist.kinetic_energy.push(kin);
            hist.damping_energy.push(hist.damping_energy[last] + e_d);
            hist.internal_work.push(hist.internal_work[last] + e_int);
            hist.iterations.push(iters);
            f_prev = f;
            p_prev = p1;
            c_prev = cvel;
        }
        Ok(hist)
    }
}

/// Eigenvalues of a small symmetric matrix by cyclic Jacobi.
fn jacobi_eigenvalues(mut a: Vec<Vec<f64>>) -> Vec<f64> {
    let n = a.len();
    let scale = a
        .iter()
        .flat_map(|r| r.iter())
        .fold(0.0f64, |m, v| m.max(v.abs()))
        .max(f64::MIN_POSITIVE);
    for _ in 0..80 {
        let mut off = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                off += a[i][j] * a[i][j];
            }
        }
        if off.sqrt() <= 1e-15 * scale {
            break;
        }
        for p in 0..n {
            for q in (p + 1)..n {
                if a[p][q].abs() <= 1e-300 {
                    continue;
                }
                let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                let t = if theta == 0.0 {
                    1.0
                } else {
                    theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt())
                };
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..n {
                    let (akp, akq) = (a[k][p], a[k][q]);
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for k in 0..n {
                    let (apk, aqk) = (a[p][k], a[q][k]);
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
            }
        }
    }
    (0..n).map(|i| a[i][i]).collect()
}
