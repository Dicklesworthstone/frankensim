//! SDF-native barrier contact (bead tfz.16, plan §8.2 [F]): the
//! Incremental Potential Contact family — a log barrier on the obstacle
//! level-set value that diverges as a sampled node approaches zero gap.
//! This fixed-obstacle lane first requires a finite,
//! strictly positive nodal gap, then accepts only trial iterates whose
//! sampled SDF queries remain positive (no penetrate-then-push-out cheats).
//! That is a conditional discrete-iterate guarantee, not an unconditional
//! continuous-trajectory nonintersection theorem. Distances are CutSdf
//! one-liner queries at deformed node positions — where generic IPC fights
//! mesh-distance computation, the chart hands us φ and ∇φ.
//!
//! Newton runs on E(u) = E_elastic(u) + Σ b(φ(x+u)) with a FILTERED
//! line search: the step is first capped by the SDF normal-closing
//! component, then every trial is audited with an actual φ > 0 query
//! before Armijo on the total energy. Rigorous interval CCD for curved
//! SDFs, continuous paths between iterates, and reusable contact-candidate
//! certification are recorded successor scope. Every accepted iterate keeps
//! the minimum sampled nodal gap seen as part of the solution record. Friction is
//! the lagged smoothed Coulomb of IPC: normal forces and tangent frames from
//! the previous outer round, a fixed physical displacement reference, and a
//! Huber-smoothed tangential potential — stick/slip is
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
        if !d.is_finite() || d <= 0.0 {
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
        if !d.is_finite() || d <= 0.0 {
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
        if !d.is_finite() || d <= 0.0 {
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
    /// Huber smoothing half-width in displacement units. Slip below this
    /// is quadratic; a velocity threshold must be multiplied by the step time.
    pub eps_v: f64,
    /// Maximum outer lag rounds; the refreshed friction residual must
    /// satisfy the problem tolerance before a solution is accepted.
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
    /// Minimum sampled nodal gap over every accepted iterate — the
    /// positive-gap audit trail (must be > 0 always).
    pub min_gap_ever: f64,
    /// Barrier energy at the solution (the ledger row).
    pub barrier_energy: f64,
    /// Per-node normal barrier reactions −∇b (zero outside the active set).
    pub reactions: Vec<[f64; 2]>,
    /// Per-node tangential reactions opposing the incremental slip.
    pub friction_reactions: Vec<[f64; 2]>,
    /// Nonnegative endpoint friction work, −Σ f_t·(u−reference_u).
    /// This incremental diagnostic is not trajectory-integrated heat.
    pub friction_work: f64,
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
        if let Some(&(dof, _)) = self.pins.iter().find(|(dof, _)| *dof >= ndof) {
            return Err(SolidError::InvalidInput {
                what: format!("contact pin dof {dof} is outside 0..{ndof}"),
            });
        }
        if let Some(&(dof, _)) = self.pins.iter().find(|(_, target)| !target.is_finite()) {
            return Err(SolidError::InvalidInput {
                what: format!("contact pin target for dof {dof} is not finite"),
            });
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

    /// Barrier energy at state u, excluding friction.
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
            total += lag.energy(u);
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

    fn residual_norm_at(
        &self,
        u: &[f64],
        load: f64,
        lag: Option<&FrictionLag>,
    ) -> Result<f64, SolidError> {
        self.validate_state_len(u)?;
        let n = self.node_count();
        let ndof = self.ndof();
        let mut r = self.hyper.unconstrained_residual(u, load)?;
        for i in 0..n {
            let x = self.deformed(u, i);
            let d = self.sdf_value(u, i)?;
            if d < self.barrier.dhat {
                let gphi = self.sdf_gradient(x)?;
                let b1 = self.barrier.d1(d);
                for a in 0..2 {
                    r[2 * i + a] += b1 * gphi[a];
                }
            }
        }
        if let Some(lag) = lag {
            lag.add_residual(u, &mut r);
        }
        for &(d, target) in &self.pins {
            r[d] = u[d] - load * target;
        }
        debug_assert_eq!(r.len(), ndof);
        Self::finite_max_norm(&r, "contact trial residual")
    }

    fn finite_max_norm(values: &[f64], label: &str) -> Result<f64, SolidError> {
        if let Some(i) = values.iter().position(|value| !value.is_finite()) {
            let value = values[i];
            return Err(SolidError::InvalidInput {
                what: format!("{label} entry {i} is not finite ({value})"),
            });
        }
        let mut norm = 0.0f64;
        for &value in values {
            norm = norm.max(value.abs());
        }
        Ok(norm)
    }

    /// Keep the exact Newton step when it descends. Sliding can remove
    /// tangential stiffness and expose an indefinite elastic tangent;
    /// then regularize the search matrix, never the physical residual.
    fn descent_direction(
        h: &[f64],
        r: &[f64],
        pins: &[(usize, f64)],
    ) -> Result<Vec<f64>, SolidError> {
        let ndof = r.len();
        Self::finite_max_norm(h, "contact tangent")?;
        let solve = |matrix: &[f64]| -> Option<Vec<f64>> {
            let factor = fs_la::factor::lu(matrix, ndof).ok()?;
            let mut direction: Vec<f64> = r.iter().map(|v| -v).collect();
            factor.solve(&mut direction);
            let slope: f64 = r.iter().zip(&direction).map(|(a, b)| a * b).sum();
            (direction.iter().all(|v| v.is_finite()) && slope.is_finite() && slope < 0.0)
                .then_some(direction)
        };
        if let Some(direction) = solve(h) {
            return Ok(direction);
        }

        let mut fixed = vec![false; ndof];
        for &(d, _) in pins {
            fixed[d] = true;
        }
        let scale = (0..ndof)
            .filter(|&i| !fixed[i])
            .map(|i| h[i * ndof + i].abs())
            .fold(0.0_f64, f64::max);
        let diagonal: Vec<f64> = (0..ndof)
            .map(|i| h[i * ndof + i].abs().max(scale * 1e-12))
            .collect();
        // A dimensionless shift follows local stiffness units. Pins keep
        // their identity equations; no artificial support enters equilibrium.
        let mut shift = 1e-6;
        for _ in 0..16 {
            let mut regularized = h.to_vec();
            for i in 0..ndof {
                if !fixed[i] {
                    regularized[i * ndof + i] += shift * diagonal[i];
                }
            }
            if let Some(direction) = solve(&regularized) {
                return Ok(direction);
            }
            shift *= 10.0;
        }
        Err(SolidError::SolveFailed {
            iters: 0,
            rel_residual: f64::INFINITY,
        })
    }

    /// Solve the barrier-augmented equilibrium at `load`.
    ///
    /// # Errors
    /// [`SolidError::NewtonStalled`] when the filtered Newton cannot
    /// reach tolerance; [`SolidError::SolveFailed`] on a singular
    /// system; [`SolidError::InvalidInput`] for infeasible or
    /// non-finite contact data; material refusals propagate.
    pub fn solve(&self, load: f64) -> Result<ContactSolution, SolidError> {
        self.solve_from(load, &vec![0.0; self.ndof()])
    }

    /// Solve an incremental equilibrium from a previous displacement.
    /// The previous state is both the initial iterate and the physical
    /// friction reference, held fixed throughout every Newton/lag round.
    ///
    /// # Errors
    /// The same refusals as [`Self::solve`], including invalid, non-finite,
    /// or intersecting reference states and unconverged friction lags.
    #[allow(clippy::too_many_lines)] // contact Newton + audited line search is one proof path
    pub fn solve_from(
        &self,
        load: f64,
        reference_u: &[f64],
    ) -> Result<ContactSolution, SolidError> {
        self.validate_inputs(load)?;
        self.validate_state_len(reference_u)?;
        let n = self.node_count();
        let ndof = self.ndof();
        let mut u = reference_u.to_vec();
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
            Some(fr) => Some(FrictionLag::capture(self, reference_u, &u, fr)?),
            None => None,
        };
        let mut residual = f64::INFINITY;
        let mut lag_history = Vec::new();
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
                    l.add_terms(&u, &mut r, &mut h);
                }
                // Pins: identity rows driving u[d] to load·target.
                for &(d, target) in &self.pins {
                    r[d] = u[d] - load * target;
                    for c in 0..ndof {
                        h[d * ndof + c] = if c == d { 1.0 } else { 0.0 };
                    }
                }
                residual = Self::finite_max_norm(&r, "contact residual")?;
                history.push(residual);
                if residual < self.tol {
                    stalled = false;
                    break;
                }
                // A non-descent Newton direction must not consume forty
                // expensive trial residuals before the gradient fallback.
                let du = Self::descent_direction(&h, &r, &self.pins)?;
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
                        if let Ok(e1) = self.total_energy(&trial, load, lag.as_ref())
                            && e1 <= (1e-4 * alpha).mul_add(slope, e0)
                        {
                            u = trial;
                            min_gap_ever = min_gap_ever.min(gap);
                            accepted = true;
                            break;
                        }
                        if self
                            .residual_norm_at(&trial, load, lag.as_ref())
                            .is_ok_and(|trial_residual| trial_residual <= 0.9 * residual)
                        {
                            u = trial;
                            min_gap_ever = min_gap_ever.min(gap);
                            accepted = true;
                            break;
                        }
                    }
                    alpha *= 0.5;
                }
                if !accepted {
                    // The contact Hessian is intentionally a PSD
                    // approximation (the SDF curvature term is
                    // dropped). If that
                    // Newton direction is not an energy descent, fall
                    // back to the true residual-gradient direction.
                    let du_sd: Vec<f64> = r.iter().map(|v| -v).collect();
                    let slope_sd: f64 = -r.iter().map(|v| v * v).sum::<f64>();
                    let mut alpha_sd = 1.0f64;
                    for i in 0..n {
                        let x = self.deformed(&u, i);
                        let d = self.sdf_value(&u, i)?;
                        let gphi = self.sdf_gradient(x)?;
                        let closing = -(gphi[0] * du_sd[2 * i] + gphi[1] * du_sd[2 * i + 1]);
                        if closing > 0.0 {
                            alpha_sd = alpha_sd.min(0.9 * d / closing);
                        }
                    }
                    alpha_sd = alpha_sd.min(1.0);
                    for _ in 0..40 {
                        let trial: Vec<f64> = u
                            .iter()
                            .zip(&du_sd)
                            .map(|(a, b)| alpha_sd.mul_add(*b, *a))
                            .collect();
                        let gap = self.min_gap(&trial)?;
                        if gap > 0.0 {
                            if let Ok(e1) = self.total_energy(&trial, load, lag.as_ref())
                                && e1 <= (1e-4 * alpha_sd).mul_add(slope_sd, e0)
                            {
                                u = trial;
                                min_gap_ever = min_gap_ever.min(gap);
                                accepted = true;
                                break;
                            }
                            if self
                                .residual_norm_at(&trial, load, lag.as_ref())
                                .is_ok_and(|trial_residual| trial_residual <= 0.9 * residual)
                            {
                                u = trial;
                                min_gap_ever = min_gap_ever.min(gap);
                                accepted = true;
                                break;
                            }
                        }
                        alpha_sd *= 0.5;
                    }
                }
                if !accepted {
                    return Err(SolidError::NewtonStalled { history });
                }
            }
            if stalled {
                return Err(SolidError::NewtonStalled { history });
            }
            // Refresh only the normal load and tangent frame. Numerical
            // iterations must never advance the physical slip reference.
            if let Some(fr) = self.friction {
                lag = Some(FrictionLag::capture(self, reference_u, &u, fr)?);
                residual = self.residual_norm_at(&u, load, lag.as_ref())?;
                lag_history.push(residual);
                if residual < self.tol {
                    break;
                }
            } else {
                break;
            }
        }
        if residual >= self.tol {
            return Err(SolidError::NewtonStalled {
                history: lag_history,
            });
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
        let mut friction_residual = vec![0.0; ndof];
        let friction_work = if let Some(lag) = &lag {
            lag.add_residual(&u, &mut friction_residual);
            lag.work(&u)
        } else {
            0.0
        };
        if !friction_work.is_finite() {
            return Err(SolidError::InvalidInput {
                what: "contact friction work is not finite".to_string(),
            });
        }
        let friction_reactions = friction_residual
            .chunks_exact(2)
            .map(|r| [-r[0], -r[1]])
            .collect();
        Ok(ContactSolution {
            min_gap: self.min_gap(&u)?,
            min_gap_ever,
            barrier_energy: self.barrier_energy(&u)?,
            reactions,
            friction_reactions,
            friction_work,
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

/// Normal force and tangent frame lagged at the previous outer iterate.
/// The physical displacement reference is fixed for the entire solve.
struct FrictionLag<'a> {
    settings: Friction,
    reference_u: &'a [f64],
    /// Per node: (unit tangent, lagged normal force magnitude).
    frames: Vec<([f64; 2], f64)>,
}

impl<'a> FrictionLag<'a> {
    fn capture(
        cp: &ContactProblem<'_>,
        reference_u: &'a [f64],
        u: &[f64],
        settings: Friction,
    ) -> Result<Self, SolidError> {
        let mut frames = Vec::with_capacity(cp.node_count());
        for i in 0..cp.node_count() {
            let d = cp.sdf_value(u, i)?;
            if d >= cp.barrier.dhat {
                frames.push(([0.0; 2], 0.0));
                continue;
            }
            let g = cp.sdf_gradient(cp.deformed(u, i))?;
            let gn = g[0].hypot(g[1]);
            if !gn.is_finite() || gn <= 0.0 {
                return Err(SolidError::InvalidInput {
                    what: format!("active contact normal at node {i} has invalid norm {gn}"),
                });
            }
            // CutSdf permits level sets with non-unit gradients: the
            // physical barrier force magnitude includes the gradient norm.
            let lambda = -cp.barrier.d1(d) * gn;
            if !lambda.is_finite() || lambda < 0.0 {
                return Err(SolidError::InvalidInput {
                    what: format!("contact normal force at node {i} is invalid ({lambda})"),
                });
            }
            frames.push(([-g[1] / gn, g[0] / gn], lambda));
        }
        Ok(Self {
            settings,
            reference_u,
            frames,
        })
    }

    /// Signed scalar slip in the frozen 2D tangent frame.
    fn slip(&self, u: &[f64], i: usize) -> f64 {
        let t = self.frames[i].0;
        t[0].mul_add(
            u[2 * i] - self.reference_u[2 * i],
            t[1] * (u[2 * i + 1] - self.reference_u[2 * i + 1]),
        )
    }

    /// Huber value: quadratic below eps_v, linear above (C¹, convex).
    fn huber(&self, s: f64) -> f64 {
        let e = self.settings.eps_v;
        let s = s.abs();
        if s <= e {
            s * s / (2.0 * e)
        } else {
            s - e / 2.0
        }
    }

    fn huber_d1(&self, s: f64) -> f64 {
        (s / self.settings.eps_v).clamp(-1.0, 1.0)
    }

    fn energy(&self, u: &[f64]) -> f64 {
        self.frames
            .iter()
            .enumerate()
            .map(|(i, (_, lambda))| self.settings.mu * lambda * self.huber(self.slip(u, i)))
            .sum()
    }

    fn add_residual(&self, u: &[f64], r: &mut [f64]) {
        for (i, &(t, lambda)) in self.frames.iter().enumerate() {
            let g = self.settings.mu * lambda * self.huber_d1(self.slip(u, i));
            for a in 0..2 {
                r[2 * i + a] += g * t[a];
            }
        }
    }

    /// Exact gradient and Hessian of the frozen-frame Huber potential.
    /// In 2D its sliding branch is linear and has zero curvature.
    fn add_terms(&self, u: &[f64], r: &mut [f64], h: &mut [f64]) {
        let ndof = r.len();
        for (i, &(t, lambda)) in self.frames.iter().enumerate() {
            let slip = self.slip(u, i);
            let coef = self.settings.mu * lambda;
            let g = coef * self.huber_d1(slip);
            let k = if slip.abs() < self.settings.eps_v {
                coef / self.settings.eps_v
            } else {
                0.0
            };
            for a in 0..2 {
                r[2 * i + a] += g * t[a];
                for b in 0..2 {
                    h[(2 * i + a) * ndof + 2 * i + b] += k * t[a] * t[b];
                }
            }
        }
    }

    fn work(&self, u: &[f64]) -> f64 {
        self.frames
            .iter()
            .enumerate()
            .map(|(i, (_, lambda))| {
                let slip = self.slip(u, i);
                self.settings.mu * lambda * self.huber_d1(slip) * slip
            })
            .sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hyper2d::NewtonSettings;
    use crate::mesh2::Mesh2;
    use fs_cutfem::sdf::Circle;
    use fs_material::hyper::{Hyperelastic, HyperelasticModel};

    #[test]
    fn contact_direction_preserves_newton_and_regularizes_ascent() {
        let positive = [4.0, 1.0, 1.0, 3.0];
        let residual = [1.0, -2.0];
        let direction = ContactProblem::descent_direction(&positive, &residual, &[])
            .expect("exact Newton direction");
        for i in 0..2 {
            let defect =
                positive[2 * i] * direction[0] + positive[2 * i + 1] * direction[1] + residual[i];
            assert!(defect.abs() < 1e-14, "well-posed Newton must stay exact");
        }

        // The unregularized solve would ascend: [1, -0.1]·[1, 0.2] > 0.
        let indefinite = [-1.0, 0.0, 0.0, 2.0];
        let residual = [1.0, 0.2];
        let direction = ContactProblem::descent_direction(&indefinite, &residual, &[])
            .expect("regularized descent");
        assert!(direction[0] < 0.0 && direction[1] < 0.0);
        assert!(
            residual
                .iter()
                .zip(&direction)
                .map(|(r, d)| r * d)
                .sum::<f64>()
                < 0.0
        );
        assert_eq!(
            indefinite,
            [-1.0, 0.0, 0.0, 2.0],
            "physical tangent unchanged"
        );

        let pinned = [1.0, 0.0, 0.0, 3.0, -1.0, 0.0, 0.0, 0.0, 2.0];
        let residual = [0.25, 1.0, 0.2];
        let direction = ContactProblem::descent_direction(&pinned, &residual, &[(0, 0.0)])
            .expect("pin-preserving descent");
        assert!(
            (direction[0] + 0.25).abs() < 1e-14,
            "pin equation unchanged"
        );
        assert!(
            residual
                .iter()
                .zip(&direction)
                .map(|(r, d)| r * d)
                .sum::<f64>()
                < 0.0
        );
    }

    #[test]
    fn friction_curved_frame_derivatives() {
        let mesh = Mesh2::quads(1.0, 1.0, 1, 1);
        let material = Hyperelastic::new(
            HyperelasticModel::NeoHookean {
                mu: 1.0,
                lambda: 3.0,
            },
            3.0,
        )
        .expect("material");
        let hyper = HyperProblem {
            mesh: &mesh,
            material: &material,
            dirichlet: vec![],
            traction: vec![],
            settings: NewtonSettings::default(),
        };
        let circle = Circle {
            center: [-1.0, -1.0],
            radius: 0.5,
        };
        let cp = ContactProblem {
            hyper: &hyper,
            sdf: &circle,
            barrier: Barrier {
                kappa: 1.0,
                dhat: 3.0,
            },
            friction: None,
            pins: vec![],
            max_newton: 30,
            tol: 1e-9,
        };
        let reference = vec![0.0; cp.ndof()];
        let lag = FrictionLag::capture(
            &cp,
            &reference,
            &reference,
            Friction {
                mu: 0.5,
                eps_v: 0.01,
                rounds: 1,
            },
        )
        .expect("lag");
        for slip in [0.0, 0.001, 0.04, -0.04] {
            let mut u = reference.clone();
            for (i, &(t, _)) in lag.frames.iter().enumerate() {
                u[2 * i] = slip * t[0];
                u[2 * i + 1] = slip * t[1];
            }
            let n = u.len();
            let mut r = vec![0.0; n];
            let mut h = vec![0.0; n * n];
            lag.add_terms(&u, &mut r, &mut h);
            let delta = 1e-7;
            for j in 0..n {
                let mut plus = u.clone();
                let mut minus = u.clone();
                plus[j] += delta;
                minus[j] -= delta;
                let gradient = (lag.energy(&plus) - lag.energy(&minus)) / (2.0 * delta);
                assert!(
                    (gradient - r[j]).abs() < 1e-7 * r[j].abs().max(1.0),
                    "slip={slip}, dof={j}: energy derivative {gradient} != {}",
                    r[j]
                );
                let mut rp = vec![0.0; n];
                let mut rm = vec![0.0; n];
                lag.add_residual(&plus, &mut rp);
                lag.add_residual(&minus, &mut rm);
                for i in 0..n {
                    let derivative = (rp[i] - rm[i]) / (2.0 * delta);
                    let hij = h[i * n + j];
                    assert!(
                        (derivative - hij).abs() < 1e-7 * hij.abs().max(1.0),
                        "slip={slip}, ({i},{j}): residual derivative {derivative} != {hij}"
                    );
                }
            }
            assert!(lag.work(&u) >= 0.0);
        }
    }
}
