//! SDF-native barrier contact (bead tfz.16, plan §8.2 [F]): the
//! Incremental Potential Contact family — a log barrier on the
//! CERTIFIED distance to the obstacle that diverges before
//! penetration, so trajectories are INTERSECTION-FREE BY CONSTRUCTION
//! (no penetrate-then-push-out cheats). Distances are CutSdf one-liner
//! queries at deformed node positions — where generic IPC fights
//! mesh-distance computation, the chart hands us φ and ∇φ.
//!
//! Newton runs on E(u) = E_elastic(u) + Σ b(φ(x+u)) with a FILTERED
//! line search: the step is first capped by the SDF normal-closing
//! component, then every trial is audited with an actual φ > 0 query
//! before Armijo on the total energy. Rigorous interval CCD for curved
//! SDFs is recorded successor scope. Every accepted iterate keeps the
//! minimum gap ever seen as part of the solution record. Friction is
//! the lagged smoothed Coulomb of IPC: normal forces from the previous
//! outer round, a Huber-smoothed tangential potential — stick/slip is
//! gated by the battery, not asserted here. The barrier Hessian keeps
//! the b''·∇φ∇φᵀ term and drops the b'·∇²φ curvature term (standard
//! IPC PSD practice; exact for planar obstacles).

use crate::SolidError;
use crate::hyper2d::HyperProblem;
use fs_cutfem::CutSdf;

/// IPC log-barrier: `b(d) = −κ(d−d̂)²·ln(d/d̂)` for `0 < d < d̂`,
/// zero for `d ≥ d̂`, divergent as `d → 0⁺`.
#[derive(Debug, Clone, Copy)]
pub struct Barrier {
    /// Stiffness κ.
    pub kappa: f64,
    /// Activation distance d̂.
    pub dhat: f64,
}

impl Barrier {
    /// Validate barrier parameters before entering a solver.
    ///
    /// # Errors
    /// [`SolidError::InvalidInput`] when κ or d̂ is non-finite or
    /// non-positive.
    pub fn validate(&self) -> Result<(), SolidError> {
        if !(self.kappa.is_finite() && self.kappa > 0.0) {
            return Err(SolidError::InvalidInput {
                what: format!(
                    "contact barrier kappa={} must be finite and positive",
                    self.kappa
                ),
            });
        }
        if !(self.dhat.is_finite() && self.dhat > 0.0) {
            return Err(SolidError::InvalidInput {
                what: format!(
                    "contact barrier dhat={} must be finite and positive",
                    self.dhat
                ),
            });
        }
        Ok(())
    }

    /// Barrier value.
    #[must_use]
    pub fn value(&self, d: f64) -> f64 {
        if d >= self.dhat {
            return 0.0;
        }
        if !(d > 0.0) {
            return f64::INFINITY;
        }
        let g = d - self.dhat;
        -self.kappa * g * g * (d / self.dhat).ln()
    }

    /// First derivative b′(d).
    #[must_use]
    pub fn d1(&self, d: f64) -> f64 {
        if d >= self.dhat {
            return 0.0;
        }
        if !(d > 0.0) {
            return f64::NAN;
        }
        let g = d - self.dhat;
        let l = (d / self.dhat).ln();
        -self.kappa * (2.0 * g).mul_add(l, g * g / d)
    }

    /// Second derivative b″(d).
    #[must_use]
    pub fn d2(&self, d: f64) -> f64 {
        if d >= self.dhat {
            return 0.0;
        }
        if !(d > 0.0) {
            return f64::NAN;
        }
        let g = d - self.dhat;
        let l = (d / self.dhat).ln();
        -self.kappa * ((2.0 * l + 4.0 * g / d) - g * g / (d * d))
    }
}

/// Lagged smoothed-Coulomb friction settings.
#[derive(Debug, Clone, Copy)]
pub struct Friction {
    /// Coulomb coefficient μ.
    pub mu: f64,
    /// Huber smoothing half-width (slip below this is quadratic).
    pub eps_v: f64,
    /// Outer lag rounds.
    pub rounds: u32,
}

impl Friction {
    /// Validate lagged friction parameters.
    ///
    /// # Errors
    /// [`SolidError::InvalidInput`] when μ, eps_v, or the round count
    /// is outside the documented contact surface.
    pub fn validate(&self) -> Result<(), SolidError> {
        if !(self.mu.is_finite() && self.mu >= 0.0) {
            return Err(SolidError::InvalidInput {
                what: format!("friction mu={} must be finite and non-negative", self.mu),
            });
        }
        if !(self.eps_v.is_finite() && self.eps_v > 0.0) {
            return Err(SolidError::InvalidInput {
                what: format!("friction eps_v={} must be finite and positive", self.eps_v),
            });
        }
        if self.rounds == 0 {
            return Err(SolidError::InvalidInput {
                what: "friction rounds must be at least one".to_string(),
            });
        }
        Ok(())
    }
}

/// The contact-augmented static problem.
pub struct ContactProblem<'a> {
    /// The elastic core (materials, Dirichlet, tractions).
    pub hyper: &'a HyperProblem<'a>,
    /// The obstacle chart (φ > 0 = free space).
    pub sdf: &'a dyn CutSdf,
    /// Barrier parameters.
    pub barrier: Barrier,
    /// Optional lagged friction.
    pub friction: Option<Friction>,
    /// Pinned DOFs `(dof index, target displacement at load 1)` —
    /// the elastic probe is unconstrained, so rigid modes the contact
    /// does not ground (e.g. tangential float on a frictionless
    /// plane) are pinned here.
    pub pins: Vec<(usize, f64)>,
    /// Newton iteration cap.
    pub max_newton: u32,
    /// Residual tolerance (max-norm).
    pub tol: f64,
}

/// The converged record: displacements plus the auditable evidence.
pub struct ContactSolution {
    /// Nodal displacements (flat, 2 per node).
    pub u: Vec<f64>,
    /// Newton iterations used (all rounds).
    pub iterations: u32,
    /// Minimum gap at the SOLUTION.
    pub min_gap: f64,
    /// Minimum gap over EVERY accepted iterate — the intersection-free
    /// audit trail (must be > 0 always).
    pub min_gap_ever: f64,
    /// Barrier energy at the solution (the ledger row).
    pub barrier_energy: f64,
    /// Per-node contact reactions −∇b (zero outside the active set).
    pub reactions: Vec<[f64; 2]>,
    /// Final residual max-norm.
    pub residual: f64,
}

impl ContactProblem<'_> {
    fn node_count(&self) -> usize {
        self.hyper.mesh.nodes.len()
    }

    fn ndof(&self) -> usize {
        2 * self.node_count()
    }

    fn validate_inputs(&self, load: f64) -> Result<(), SolidError> {
        if self.node_count() == 0 {
            return Err(SolidError::InvalidInput {
                what: "contact problem needs at least one node".to_string(),
            });
        }
        if !load.is_finite() {
            return Err(SolidError::InvalidInput {
                what: format!("load={load} must be finite"),
            });
        }
        if !(self.tol.is_finite() && self.tol > 0.0) {
            return Err(SolidError::InvalidInput {
                what: format!("contact tolerance={} must be finite and positive", self.tol),
            });
        }
        if self.max_newton == 0 {
            return Err(SolidError::InvalidInput {
                what: "contact max_newton must be at least one".to_string(),
            });
        }
        self.barrier.validate()?;
        if let Some(friction) = self.friction {
            friction.validate()?;
        }
        let ndof = self.ndof();
        for &(dof, target) in &self.pins {
            if dof >= ndof {
                return Err(SolidError::InvalidInput {
                    what: format!("contact pin dof {dof} is outside 0..{ndof}"),
                });
            }
            if !target.is_finite() {
                return Err(SolidError::InvalidInput {
                    what: format!("contact pin target for dof {dof} is not finite"),
                });
            }
        }
        Ok(())
    }

    fn validate_state_len(&self, u: &[f64]) -> Result<(), SolidError> {
        let ndof = self.ndof();
        if u.len() != ndof {
            return Err(SolidError::InvalidInput {
                what: format!("contact state has {} entries, expected {ndof}", u.len()),
            });
        }
        if let Some((i, v)) = u.iter().enumerate().find(|(_, v)| !v.is_finite()) {
            return Err(SolidError::InvalidInput {
                what: format!("contact state entry {i} is not finite ({v})"),
            });
        }
        Ok(())
    }

    /// Deformed position of node i.
    fn deformed(&self, u: &[f64], i: usize) -> [f64; 2] {
        let p = self.hyper.mesh.nodes[i];
        [p[0] + u[2 * i], p[1] + u[2 * i + 1]]
    }

    fn sdf_value(&self, u: &[f64], i: usize) -> Result<f64, SolidError> {
        let x = self.deformed(u, i);
        let d = self.sdf.value(x);
        if d.is_finite() {
            Ok(d)
        } else {
            Err(SolidError::InvalidInput {
                what: format!("contact SDF returned non-finite value {d} at node {i}"),
            })
        }
    }

    fn sdf_gradient(&self, x: [f64; 2]) -> Result<[f64; 2], SolidError> {
        let g = self.sdf.gradient(x);
        if g[0].is_finite() && g[1].is_finite() {
            Ok(g)
        } else {
            Err(SolidError::InvalidInput {
                what: format!(
                    "contact SDF returned non-finite gradient [{}, {}]",
                    g[0], g[1]
                ),
            })
        }
    }

    /// Minimum gap over all nodes at state u.
    fn min_gap(&self, u: &[f64]) -> Result<f64, SolidError> {
        self.validate_state_len(u)?;
        let mut min_gap = f64::INFINITY;
        for i in 0..self.node_count() {
            min_gap = min_gap.min(self.sdf_value(u, i)?);
        }
        Ok(min_gap)
    }

    /// Barrier energy at state u (plus friction if anchors given).
    fn barrier_energy(&self, u: &[f64]) -> Result<f64, SolidError> {
        self.validate_state_len(u)?;
        let mut total = 0.0;
        for i in 0..self.node_count() {
            total += self.barrier.value(self.sdf_value(u, i)?);
        }
        Ok(total)
    }

    /// Total energy for the line search (elastic may refuse a state).
    fn total_energy(
        &self,
        u: &[f64],
        load: f64,
        lag: Option<&FrictionLag>,
    ) -> Result<f64, SolidError> {
        let e = self
            .hyper
            .potential_energy(u, load)
            .ok_or(SolidError::SolveFailed {
                iters: 0,
                rel_residual: f64::INFINITY,
            })?;
        let mut total = e + self.barrier_energy(u)?;
        if let Some(lag) = lag {
            total += lag.energy(self, u)?;
        }
        if total.is_finite() {
            Ok(total)
        } else {
            Err(SolidError::SolveFailed {
                iters: 0,
                rel_residual: f64::INFINITY,
            })
        }
    }

    /// Solve the barrier-augmented equilibrium at `load`.
    ///
    /// # Errors
    /// [`SolidError::NewtonStalled`] when the filtered Newton cannot
    /// reach tolerance; [`SolidError::SolveFailed`] on a singular
    /// system; [`SolidError::InvalidInput`] for infeasible or
    /// non-finite contact data; material refusals propagate.
    pub fn solve(&self, load: f64) -> Result<ContactSolution, SolidError> {
        self.validate_inputs(load)?;
        let n = self.node_count();
        let ndof = self.ndof();
        let mut u = vec![0.0f64; ndof];
        let initial_gap = self.min_gap(&u)?;
        if initial_gap <= 0.0 {
            return Err(SolidError::InvalidInput {
                what: format!(
                    "IPC requires an intersection-free initial state; min gap is {initial_gap:.3e}"
                ),
            });
        }
        let mut min_gap_ever = initial_gap;
        let mut iterations = 0u32;
        let rounds = self.friction.map_or(1, |f| f.rounds);
        // The lag starts from the INITIAL state: the barrier is
        // already active there in any sensible fixture, so friction
        // grounds tangential rigid modes from the first round (a
        // frictionless first round with no pins was measured singular).
        let mut lag: Option<FrictionLag> = match self.friction {
            Some(fr) => Some(FrictionLag::capture(self, &u, fr)?),
            None => None,
        };
        let mut residual = f64::INFINITY;
        for _round in 0..rounds {
            let mut stalled = true;
            let mut history = Vec::new();
            for _ in 0..self.max_newton {
                iterations += 1;
                let (mut r, k) = self.hyper.residual_and_tangent(&u, load)?;
                let mut h = k.to_dense();
                // Barrier terms (+ friction if lagged).
                for i in 0..n {
                    let x = self.deformed(&u, i);
                    let d = self.sdf_value(&u, i)?;
                    if d < self.barrier.dhat {
                        let gphi = self.sdf_gradient(x)?;
                        let b1 = self.barrier.d1(d);
                        let b2 = self.barrier.d2(d);
                        for a in 0..2 {
                            r[2 * i + a] += b1 * gphi[a];
                            for b in 0..2 {
                                h[(2 * i + a) * ndof + (2 * i + b)] += b2 * gphi[a] * gphi[b];
                            }
                        }
                    }
                }
                if let Some(l) = &lag {
                    l.add_terms(self, &u, &mut r, &mut h)?;
                }
                // Pins: identity rows driving u[d] to load·target.
                for &(d, target) in &self.pins {
                    r[d] = u[d] - load * target;
                    for c in 0..ndof {
                        h[d * ndof + c] = if c == d { 1.0 } else { 0.0 };
                    }
                }
                residual = r.iter().fold(0.0f64, |m, &v| m.max(v.abs()));
                history.push(residual);
                if residual < self.tol {
                    stalled = false;
                    break;
                }
                // Newton direction.
                let f = fs_la::factor::lu(&h, ndof).map_err(|_| SolidError::SolveFailed {
                    iters: 0,
                    rel_residual: f64::INFINITY,
                })?;
                let mut du: Vec<f64> = r.iter().map(|v| -v).collect();
                f.solve(&mut du);
                // Normal-closing filter: tangential slip should not be
                // throttled by the normal gap. The subsequent φ > 0
                // audit rejects any curved-SDF trial this first-order
                // cap does not already exclude.
                let mut alpha = 1.0f64;
                for i in 0..n {
                    let x = self.deformed(&u, i);
                    let d = self.sdf_value(&u, i)?;
                    let gphi = self.sdf_gradient(x)?;
                    let closing = -(gphi[0] * du[2 * i] + gphi[1] * du[2 * i + 1]);
                    if closing > 0.0 {
                        alpha = alpha.min(0.9 * d / closing);
                    }
                }
                alpha = alpha.min(1.0);
                // Armijo on the total energy with the φ > 0 audit.
                let e0 = self.total_energy(&u, load, lag.as_ref())?;
                let slope: f64 = r.iter().zip(&du).map(|(a, b)| a * b).sum();
                let mut accepted = false;
                for _ in 0..40 {
                    let trial: Vec<f64> = u
                        .iter()
                        .zip(&du)
                        .map(|(a, b)| alpha.mul_add(*b, *a))
                        .collect();
                    let gap = self.min_gap(&trial)?;
                    if gap > 0.0 {
                        if let Ok(e1) = self.total_energy(&trial, load, lag.as_ref()) {
                            if e1 <= (1e-4 * alpha).mul_add(slope, e0) {
                                u = trial;
                                min_gap_ever = min_gap_ever.min(gap);
                                accepted = true;
                                break;
                            }
                        }
                    }
                    alpha *= 0.5;
                }
                if !accepted {
                    return Err(SolidError::NewtonStalled { history });
                }
            }
            if stalled {
                return Err(SolidError::NewtonStalled { history });
            }
            // Refresh the friction lag from the converged normal set.
            if let Some(fr) = self.friction {
                lag = Some(FrictionLag::capture(self, &u, fr)?);
            } else {
                break;
            }
        }
        let mut reactions = vec![[0.0f64; 2]; n];
        for (i, reaction) in reactions.iter_mut().enumerate() {
            let x = self.deformed(&u, i);
            let d = self.sdf_value(&u, i)?;
            if d < self.barrier.dhat {
                let gphi = self.sdf_gradient(x)?;
                let b1 = self.barrier.d1(d);
                *reaction = [-b1 * gphi[0], -b1 * gphi[1]];
            }
        }
        Ok(ContactSolution {
            min_gap: self.min_gap(&u)?,
            min_gap_ever,
            barrier_energy: self.barrier_energy(&u)?,
            reactions,
            iterations,
            residual,
            u,
        })
    }

    /// Equilibrium-constrained gradient of `J = j·u` with respect to a
    /// rigid TRANSLATION of the obstacle along unit direction `e`
    /// (differentiable contact): solve `Hᵀλ = j` at the solution and
    /// return `−λᵀ ∂r/∂h`, with `∂r/∂h = −b″(∇φ·e)∇φ` per active node
    /// (exact for planar obstacles; first-order in curvature
    /// otherwise — documented, gated against FD by the battery).
    ///
    /// # Errors
    /// [`SolidError::SolveFailed`] on a singular adjoint system;
    /// [`SolidError::InvalidInput`] for non-finite or unsupported
    /// contact-adjoint inputs.
    pub fn translation_gradient(
        &self,
        sol: &ContactSolution,
        load: f64,
        j: &[f64],
        e: [f64; 2],
    ) -> Result<f64, SolidError> {
        self.validate_inputs(load)?;
        self.validate_state_len(&sol.u)?;
        let n = self.node_count();
        let ndof = self.ndof();
        if j.len() != ndof {
            return Err(SolidError::InvalidInput {
                what: format!(
                    "contact adjoint seed has {} entries, expected {ndof}",
                    j.len()
                ),
            });
        }
        if !j.iter().all(|v| v.is_finite()) {
            return Err(SolidError::InvalidInput {
                what: "contact adjoint seed contains a non-finite entry".to_string(),
            });
        }
        if !(e[0].is_finite() && e[1].is_finite()) {
            return Err(SolidError::InvalidInput {
                what: format!(
                    "contact translation direction [{}, {}] is non-finite",
                    e[0], e[1]
                ),
            });
        }
        if self.friction.is_some() {
            return Err(SolidError::InvalidInput {
                what: "contact translation_gradient is currently frictionless-only".to_string(),
            });
        }
        let (_, k) = self.hyper.residual_and_tangent(&sol.u, load)?;
        let mut h = k.to_dense();
        for i in 0..n {
            let x = self.deformed(&sol.u, i);
            let d = self.sdf_value(&sol.u, i)?;
            if d < self.barrier.dhat {
                let gphi = self.sdf_gradient(x)?;
                let b2 = self.barrier.d2(d);
                for a in 0..2 {
                    for b in 0..2 {
                        h[(2 * i + a) * ndof + (2 * i + b)] += b2 * gphi[a] * gphi[b];
                    }
                }
            }
        }
        for &(d, _) in &self.pins {
            for c in 0..ndof {
                h[d * ndof + c] = if c == d { 1.0 } else { 0.0 };
            }
        }
        // Transpose for the adjoint (H is symmetric here, but keep the
        // transpose explicit against future nonsymmetric terms).
        let mut ht = vec![0.0f64; ndof * ndof];
        for r in 0..ndof {
            for c in 0..ndof {
                ht[c * ndof + r] = h[r * ndof + c];
            }
        }
        let f = fs_la::factor::lu(&ht, ndof).map_err(|_| SolidError::SolveFailed {
            iters: 0,
            rel_residual: f64::INFINITY,
        })?;
        let mut lam = j.to_vec();
        f.solve(&mut lam);
        let mut dj = 0.0f64;
        for i in 0..n {
            let x = self.deformed(&sol.u, i);
            let d = self.sdf_value(&sol.u, i)?;
            if d < self.barrier.dhat {
                let gphi = self.sdf_gradient(x)?;
                let b2 = self.barrier.d2(d);
                let de_dh = -(gphi[0] * e[0] + gphi[1] * e[1]);
                for a in 0..2 {
                    let dr_dh = b2 * de_dh * gphi[a];
                    dj -= lam[2 * i + a] * dr_dh;
                }
            }
        }
        Ok(dj)
    }
}

/// The lagged friction state: anchors and normal forces captured at
/// the last converged round.
struct FrictionLag {
    settings: Friction,
    /// Per node: (anchor position, lagged normal force magnitude).
    anchors: Vec<([f64; 2], f64)>,
}

impl FrictionLag {
    fn capture(
        cp: &ContactProblem<'_>,
        u: &[f64],
        settings: Friction,
    ) -> Result<FrictionLag, SolidError> {
        let n = cp.node_count();
        let mut anchors = Vec::with_capacity(n);
        for i in 0..n {
            let x = cp.deformed(u, i);
            let d = cp.sdf_value(u, i)?;
            let lambda = if d < cp.barrier.dhat {
                -cp.barrier.d1(d)
            } else {
                0.0
            };
            anchors.push((x, lambda.max(0.0)));
        }
        Ok(FrictionLag { settings, anchors })
    }

    /// Tangential slip of node i relative to its anchor.
    fn slip(
        &self,
        cp: &ContactProblem<'_>,
        u: &[f64],
        i: usize,
    ) -> Result<([f64; 2], f64), SolidError> {
        let x = cp.deformed(u, i);
        let (anchor, _) = self.anchors[i];
        let dx = [x[0] - anchor[0], x[1] - anchor[1]];
        let gphi = cp.sdf_gradient(x)?;
        let gn = gphi[0].hypot(gphi[1]).max(1e-30);
        let nrm = [gphi[0] / gn, gphi[1] / gn];
        let dn = dx[0].mul_add(nrm[0], dx[1] * nrm[1]);
        let t = [dx[0] - dn * nrm[0], dx[1] - dn * nrm[1]];
        Ok((t, t[0].hypot(t[1])))
    }

    /// Huber value: quadratic below eps_v, linear above (C¹, convex).
    fn huber(&self, s: f64) -> f64 {
        let e = self.settings.eps_v;
        if s <= e {
            s * s / (2.0 * e)
        } else {
            s - e / 2.0
        }
    }

    fn huber_d1(&self, s: f64) -> f64 {
        let e = self.settings.eps_v;
        if s <= e { s / e } else { 1.0 }
    }

    fn energy(&self, cp: &ContactProblem<'_>, u: &[f64]) -> Result<f64, SolidError> {
        let mut total = 0.0;
        for i in 0..cp.node_count() {
            let (_, s) = self.slip(cp, u, i)?;
            total += self.settings.mu * self.anchors[i].1 * self.huber(s);
        }
        Ok(total)
    }

    /// Add friction gradient and a convex (Gauss–Newton) Hessian.
    fn add_terms(
        &self,
        cp: &ContactProblem<'_>,
        u: &[f64],
        r: &mut [f64],
        h: &mut [f64],
    ) -> Result<(), SolidError> {
        let ndof = r.len();
        for i in 0..cp.node_count() {
            let lam = self.anchors[i].1;
            if lam <= 0.0 {
                continue;
            }
            let (t, s) = self.slip(cp, u, i)?;
            let coef = self.settings.mu * lam;
            let k_iso = coef / self.settings.eps_v.max(1e-12);
            if s < 1e-30 {
                // Zero slip: the Huber quadratic branch's Hessian is
                // the ISOTROPIC tangential projector k·(I − nnᵀ) —
                // well-defined with no slip direction, and exactly what
                // grounds tangential rigid modes at the first iterate
                // (skipping it was measured singular).
                let x = cp.deformed(u, i);
                let gphi = cp.sdf_gradient(x)?;
                let gn = gphi[0].hypot(gphi[1]).max(1e-30);
                let nrm = [gphi[0] / gn, gphi[1] / gn];
                for a in 0..2 {
                    for b in 0..2 {
                        let proj = f64::from(u8::from(a == b)) - nrm[a] * nrm[b];
                        h[(2 * i + a) * ndof + (2 * i + b)] += k_iso * proj;
                    }
                }
                continue;
            }
            let dir = [t[0] / s, t[1] / s];
            let g1 = coef * self.huber_d1(s);
            for a in 0..2 {
                r[2 * i + a] += g1 * dir[a];
            }
            // Convex surrogate Hessian: (coef/eps_v)·P_t (the quadratic
            // branch's exact Hessian; kept PSD in the linear branch).
            let k = coef / self.settings.eps_v.max(1e-12);
            for a in 0..2 {
                for b in 0..2 {
                    h[(2 * i + a) * ndof + (2 * i + b)] += k * dir[a] * dir[b];
                }
            }
        }
        Ok(())
    }
}
