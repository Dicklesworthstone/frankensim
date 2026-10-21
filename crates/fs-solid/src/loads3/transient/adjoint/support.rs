//! Full-state pullbacks for prescribed-support Newmark dynamics.
//!
//! With b=h²/4, d=h/2, the reduced multiplier solves the same transposed
//! M_ff+d*C_ff+b*K_ff equation as the force-only adjoint. Its contributions to
//! support acceleration, velocity and displacement are -M_fc^T*lambda,
//! -C_fc^T*lambda and -K_fc^T*lambda. Reverse the actual support recurrence
//! a_c1=(s1-s0-h*v_c0-b*a_c0)/b and v_c1=v_c0+d*(a_c0+a_c1).
//! This includes direct objective dependence on constrained motion. No primal
//! tape is needed for this linear map at fixed steps, geometry and material.
//!
//! These are derivatives of the admitted discrete equations, not derivatives
//! of PCG, acceptance gates or rounding. Stress and work/energy report
//! objectives, frame motion and event times require their own additional
//! pullbacks. `reaction` supplies support-force observations and
//! `TetTransientRun3::material_scale_pullback` adds full nodal material and
//! damping contractions. Displacement impulses are not supported.

use super::*;

mod material;
pub mod reaction;
pub mod trajectory;

/// Full node-major objective cotangents, including constrained components.
/// Units are objective/m, objective/(m/s), and objective/(m/s²).
#[derive(Debug, Clone, PartialEq)]
pub struct NodalMotionCotangent3 {
    pub q: Vec<[f64; 3]>,
    pub v: Vec<[f64; 3]>,
    pub a: Vec<[f64; 3]>,
}

impl NodalMotionCotangent3 {
    fn validate(&self, nodes: usize, cx: &Cx<'_>) -> Result<(), TransientAdjointError3> {
        check(cx)?;
        for field in [&self.q, &self.v, &self.a] {
            if field.len() != nodes {
                return Err(TransientAdjointError3::Invalid(
                    "full motion cotangent needs one xyz vector per node",
                ));
            }
            for chunk in field.chunks(256) {
                check(cx)?;
                for &value in chunk.iter().flatten() {
                    if !value.is_finite() {
                        return Err(TransientAdjointError3::Invalid("nonfinite nodal cotangent"));
                    }
                }
            }
        }
        Ok(())
    }
}

/// Bounds for full-node scratch/output and each canonical element sweep.
/// The reduced solve and sparse transpose retain `TransientAdjointConfig3`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SupportAdjointBudget3 {
    pub maximum_full_dofs: usize,
    pub maximum_elements: usize,
}

impl Default for SupportAdjointBudget3 {
    fn default() -> Self {
        Self {
            maximum_full_dofs: 150_000,
            maximum_elements: 100_000,
        }
    }
}

/// Reverse of one fixed-duration prescribed-support interval.
#[derive(Debug, Clone, PartialEq)]
pub struct PrescribedStepPullback3 {
    /// Full previous q/v/a, treated as independent at this boundary. Chain
    /// the preceding interval or `pullback_initialization` to continue reverse.
    pub previous: NodalMotionCotangent3,
    /// Endpoint external/equivalent thermal force, in `free_dofs()` order.
    /// Fixed-node external force has no effect on the admitted motion map.
    pub force: Vec<f64>,
    /// Derivative with respect to absolute endpoint support displacement [m],
    /// in sorted `fixed_dofs()` order. Includes direct constrained-state seeds.
    /// An omitted primal row prescribed zero; vary it by explicitly adding it.
    pub prescribed_displacement: Vec<f64>,
    pub iterations: usize,
    pub true_relative_residual: f64,
}

/// Reverse of the initial coupled mass solve with explicit support kinematics.
#[derive(Debug, Clone, PartialEq)]
pub struct PrescribedInitialPullback3 {
    /// Derivatives of full initial q/v and prescribed initial acceleration.
    /// Free acceleration was solved, so its input cotangent is exactly zero.
    pub initial: NodalMotionCotangent3,
    /// Initial external/equivalent thermal force, in `free_dofs()` order.
    pub force: Vec<f64>,
    pub iterations: usize,
    pub true_relative_residual: f64,
}

/// Frozen transposed operators and original element owner for one exact run.
/// Clones share sparse matrices and fixed-DOF indexing; calls publish no partial
/// gradient and do not modify the primal trajectory or this companion.
#[derive(Clone)]
pub struct PrescribedTransientAdjoint3<'g, 'm, 'p, 's> {
    reduced: TransientAdjoint3<'g, 'm>,
    problem: &'p TetLinearElasticProblem<'s>,
    fixed: Arc<[usize]>,
    full_dofs: usize,
}

impl std::fmt::Debug for PrescribedTransientAdjoint3<'_, '_, '_, '_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PrescribedTransientAdjoint3")
            .field("free_dofs", &self.free_dofs())
            .field("fixed_dofs", &self.fixed_dofs())
            .finish_non_exhaustive()
    }
}

impl<'g, 'm, 'p, 's> TetTransientRun3<'g, 'm, 'p, 's> {
    /// Prepare the explicit full-state/support-control transpose companion.
    /// Uses the original immutable reference mesh, fixed world component set,
    /// materials, damping and frozen thermal strain. Extra full-node and element
    /// allowances are checked before allocating the companion or its scratch.
    pub fn prescribed_adjoint(
        &self,
        config: TransientAdjointConfig3,
        budget: SupportAdjointBudget3,
        cx: &Cx<'_>,
    ) -> Result<PrescribedTransientAdjoint3<'g, 'm, 'p, 's>, TransientAdjointError3> {
        check(cx)?;
        if self.prescribed.is_none() {
            return Err(TransientAdjointError3::Invalid(
                "prescribed-support adjoint requires a support-enabled trajectory",
            ));
        }
        let full_dofs = self
            .problem
            .nodes_m
            .len()
            .checked_mul(3)
            .ok_or(TransientAdjointError3::Budget)?;
        if full_dofs > budget.maximum_full_dofs
            || self.problem.tetrahedra.len() > budget.maximum_elements
        {
            return Err(TransientAdjointError3::Budget);
        }
        let mut fixed = Vec::with_capacity(self.problem.fixed_dofs.len());
        let mut free = self.free_dofs().iter().copied().peekable();
        for d in 0..full_dofs {
            if d % 256 == 0 {
                check(cx)?;
            }
            if free.peek() == Some(&d) {
                free.next();
            } else {
                fixed.push(d);
            }
        }
        let reduced = self.prepare_adjoint(config, cx)?;
        check(cx)?;
        Ok(PrescribedTransientAdjoint3 {
            reduced,
            problem: self.problem,
            fixed: Arc::from(fixed),
            full_dofs,
        })
    }
}

impl PrescribedTransientAdjoint3<'_, '_, '_, '_> {
    pub fn free_dofs(&self) -> &[usize] {
        self.reduced.free_dofs()
    }
    pub fn fixed_dofs(&self) -> &[usize] {
        &self.fixed
    }

    fn reduce(
        &self,
        seed: &NodalMotionCotangent3,
        cx: &Cx<'_>,
    ) -> Result<MotionCotangent3, TransientAdjointError3> {
        seed.validate(self.problem.nodes_m.len(), cx)?;
        let mut reduced = MotionCotangent3 {
            q: Vec::with_capacity(self.free_dofs().len()),
            v: Vec::with_capacity(self.free_dofs().len()),
            a: Vec::with_capacity(self.free_dofs().len()),
        };
        for &d in self.free_dofs() {
            check(cx)?;
            reduced.q.push(seed.q[d / 3][d % 3]);
            reduced.v.push(seed.v[d / 3][d % 3]);
            reduced.a.push(seed.a[d / 3][d % 3]);
        }
        Ok(reduced)
    }

    fn scatter(
        &self,
        reduced: &MotionCotangent3,
        cx: &Cx<'_>,
    ) -> Result<NodalMotionCotangent3, TransientAdjointError3> {
        let nodes = self.problem.nodes_m.len();
        let mut out = NodalMotionCotangent3 {
            q: vec![[0.0; 3]; nodes],
            v: vec![[0.0; 3]; nodes],
            a: vec![[0.0; 3]; nodes],
        };
        for (i, &d) in self.free_dofs().iter().enumerate() {
            if i % 256 == 0 {
                check(cx)?;
            }
            out.q[d / 3][d % 3] = reduced.q[i];
            out.v[d / 3][d % 3] = reduced.v[i];
            out.a[d / 3][d % 3] = reduced.a[i];
        }
        Ok(out)
    }

    // Embed only free equation multipliers, then transpose the original tiles.
    // The resulting fixed entries are exactly the cross-block contractions.
    fn cross_transpose(
        &self,
        lambda: &[f64],
        cx: &Cx<'_>,
    ) -> Result<(Vec<f64>, Vec<f64>), TransientAdjointError3> {
        validate_values(lambda, self.free_dofs().len(), cx)?;
        let mut full = vec![0.0; self.full_dofs];
        for (i, &d) in self.free_dofs().iter().enumerate() {
            if i % 256 == 0 {
                check(cx)?;
            }
            full[d] = lambda[i];
        }
        self.problem
            .transpose_full_kinematic_actions(&full, cx)
            .map_err(TransientError3::from)
            .map_err(TransientAdjointError3::from)
    }

    /// Reverse a fixed positive-duration step, including every constrained
    /// motion component and the full inertial coupling. Initial/previous states
    /// are independent at this boundary; chain preceding pullbacks to recover
    /// their true dependencies. No primal state or solver replay is required.
    pub fn pullback_step(
        &self,
        h_s: f64,
        seed: &NodalMotionCotangent3,
        cx: &Cx<'_>,
    ) -> Result<PrescribedStepPullback3, TransientAdjointError3> {
        let reduced_seed = self.reduce(seed, cx)?;
        let reduced = self.reduced.pullback_step(h_s, &reduced_seed, cx)?;
        let (mt, kt) = self.cross_transpose(&reduced.force, cx)?;
        let mut previous = self.scatter(&reduced.previous, cx)?;
        let mut prescribed_displacement = Vec::with_capacity(self.fixed.len());
        let b = number(0.25 * h_s * h_s)?;
        let d = 0.5 * h_s;
        let eta = self.reduced.primal_config.mass_damping_per_s;
        let zeta = self.reduced.primal_config.stiffness_damping_s;
        for &fixed in self.fixed_dofs() {
            check(cx)?;
            let node = fixed / 3;
            let component = fixed % 3;
            let q = number(seed.q[node][component] - kt[fixed])?;
            let v = number(seed.v[node][component] - number(eta * mt[fixed] + zeta * kt[fixed])?)?;
            let a = number(seed.a[node][component] - mt[fixed])?;
            let acceleration_seed = number(a + d * v)?;
            let displacement_seed = number(acceleration_seed / b)?;
            prescribed_displacement.push(number(q + displacement_seed)?);
            previous.q[node][component] = -displacement_seed;
            previous.v[node][component] = number(v - h_s * displacement_seed)?;
            previous.a[node][component] = number(d * v - acceleration_seed)?;
        }
        check(cx)?;
        Ok(PrescribedStepPullback3 {
            previous,
            force: reduced.force,
            prescribed_displacement,
            iterations: reduced.iterations,
            true_relative_residual: reduced.true_relative_residual,
        })
    }

    /// Reverse the initial mass solve with independent full initial q/v and
    /// constrained a. The primal initializer requires free a input to be zero;
    /// that unused input receives exactly zero cotangent. This is not a load
    /// jump API and cannot manufacture an impulse at the initial endpoint.
    pub fn pullback_initialization(
        &self,
        seed: &NodalMotionCotangent3,
        cx: &Cx<'_>,
    ) -> Result<PrescribedInitialPullback3, TransientAdjointError3> {
        let reduced_seed = self.reduce(seed, cx)?;
        let reduced = self.reduced.pullback_acceleration(&reduced_seed, cx)?;
        let (mt, kt) = self.cross_transpose(&reduced.force, cx)?;
        let mut initial = self.scatter(&reduced.previous, cx)?;
        let eta = self.reduced.primal_config.mass_damping_per_s;
        let zeta = self.reduced.primal_config.stiffness_damping_s;
        for &fixed in self.fixed_dofs() {
            check(cx)?;
            let node = fixed / 3;
            let component = fixed % 3;
            initial.q[node][component] = number(seed.q[node][component] - kt[fixed])?;
            initial.v[node][component] =
                number(seed.v[node][component] - number(eta * mt[fixed] + zeta * kt[fixed])?)?;
            initial.a[node][component] = number(seed.a[node][component] - mt[fixed])?;
        }
        check(cx)?;
        Ok(PrescribedInitialPullback3 {
            initial,
            force: reduced.force,
            iterations: reduced.iterations,
            true_relative_residual: reduced.true_relative_residual,
        })
    }

    /// Chain free-force cotangents through the existing load integration owner.
    pub fn contract_load_direction(
        &self,
        force: &[f64],
        direction: &TetLoadCase3<'_, '_>,
        cx: &Cx<'_>,
    ) -> Result<f64, TransientAdjointError3> {
        self.reduced.contract_load_direction(force, direction, cx)
    }
}
