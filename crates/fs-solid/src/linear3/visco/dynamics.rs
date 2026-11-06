//! Midpoint inertia with the original material's exact work-mean force.
//!
//! With w=(v1+v0)/2 and u1-u0=h*w, the original affine strain ramp gives
//! f_mean=f_hold_mean+h*K_mean*w. The sparse solve is
//! (M+h^2*K_mean/2)w=M*v0+h*(f_external_mean-f_hold_mean)/2.
//! External nodal work plus prescribed eigenstrain work equals kinetic/storage
//! change plus the original physical Maxwell/Kelvin heat, up to the independently
//! checked solve error. Frozen coefficients admit thermal forcing through the
//! SAME mean tangent; thermal storage/capacity is not part of this work port.
//! Endpoint stress is retained for observation, never substituted for f_mean.
//! Explicit support motion uses the original full consistent mass, including
//! free/support cross-couplings. Its independently recovered mean reactions
//! contribute their signed work; neither support work nor a solve defect is heat.

use super::quasistatic::{assemble_tangents, check_residual, solve_linear};
use super::*;
use crate::linear::Jacobi;
use crate::linear3::{TetStaticSolveConfig, euclidean_norm};
use fs_math::dd::Dd;
use fs_math::quadratic::{QuadraticError, increment};
use fs_solver::op::CsrOp;

mod source;
pub use source::ViscoDynamicMaterialChange;

mod impulse;
pub use impulse::{ViscoDynamicImpulse, ViscoDynamicImpulseProposal};

mod ports;
pub use ports::{
    ViscoDynamicPortBudget, ViscoDynamicPortOperator, ViscoDynamicPortResponse,
    ViscoDynamicPortStep,
};

/// Explicit trajectory, sparse solve and energy-balance limits.
#[derive(Clone, Copy, Debug)]
pub struct ViscoDynamicConfig {
    /// Maximum accepted intervals over this trajectory.
    pub maximum_steps: usize,
    /// Maximum admitted free displacement/velocity coordinates.
    pub maximum_free_dofs: usize,
    /// Shared sparse CG iteration and true-residual controls.
    pub solve: TetStaticSolveConfig,
    /// Numerical work-balance allowance; never converted into material heat.
    pub absolute_energy_tolerance_j: f64,
    /// Relative allowance against interval work, storage increments and heat.
    /// An unchanged kinetic or material inventory cannot enlarge this allowance.
    pub relative_energy_tolerance: f64,
}
impl Default for ViscoDynamicConfig {
    fn default() -> Self {
        Self {
            maximum_steps: 10_000,
            maximum_free_dofs: 50_000,
            solve: TetStaticSolveConfig::default(),
            absolute_energy_tolerance_j: 1e-12,
            relative_energy_tolerance: 1e-10,
        }
    }
}

/// Reusable fixed-duration sparse operator bound to the complete original solid.
/// Free bodies are admitted: consistent mass removes the acceleration nullspace.
pub struct ViscoDynamicOperator<'a> {
    system: &'a ViscoElements,
    dt_s: f64,
    config: ViscoDynamicConfig,
    reduced: Vec<Option<usize>>,
    mean_tangents: Vec<Mandel6>,
    operator: CsrOp,
    preconditioner: Jacobi,
}

/// Complete accepted motion and original constitutive memory.
#[derive(Clone)]
pub struct ViscoDynamicState<'p, 'a> {
    prepared: &'p ViscoDynamicOperator<'a>,
    material: ViscoElementState<'a>,
    velocity_mps: Vec<f64>,
    kinetic_j: f64,
    external_work_j: f64,
    support_work_j: f64,
    eigenstrain_work_j: f64,
    parameter_work_j: f64,
    time_s: f64,
    steps: usize,
}

/// One accepted interval. Work and heat are independently integrated.
pub struct ViscoDynamicStep {
    /// Accepted time since starting this inertial trajectory [s].
    pub time_s: f64,
    /// Sparse CG iterations used by this interval.
    pub iterations: usize,
    /// Independently recomputed mean momentum-equation residual [N].
    pub force_residual_n: f64,
    /// Endpoint kinetic energy from the original consistent mass [J].
    pub kinetic_energy_j: f64,
    /// Signed kinetic storage increment from the actual endpoint velocities
    /// and original consistent mass [J], including moving-support cross terms.
    pub kinetic_energy_change_j: f64,
    /// Signed material storage increment, reduced from the original per-element
    /// constitutive increments [J]. Unchanged stored inventories cancel locally.
    pub material_stored_energy_change_j: f64,
    /// Actual interval-mean applied nodal force times actual displacement [J].
    pub external_work_j: f64,
    /// Signed mean support-reaction work over the actual displacement [J].
    pub support_work_j: f64,
    /// Full mean inertia plus internal force minus applied force [N]. Fixed
    /// rows are forces exerted by supports on the solid; free rows retain the
    /// physical residual. None when preparation omitted supported mass rows.
    pub mean_reaction_force_n: Option<Vec<f64>>,
    /// Applied-load, support and imposed-eigenstrain work minus kinetic/storage
    /// changes and original physical heat [J].
    pub balance_residual_j: f64,
    /// Endpoint and exact mean internal forces; original storage/work/heat.
    pub constitutive: ViscoForceStep,
}

fn invalid(what: &'static str) -> ViscoElementError {
    ViscoElementError::InvalidHistory { what }
}
fn finite(value: f64) -> Result<f64, ViscoElementError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ViscoElementError::NonFiniteResult)
    }
}
fn mass_action(system: &ViscoElements, velocity: &[f64]) -> Vec<f64> {
    let mut momentum = vec![0.0; velocity.len()];
    mass_action_into(system, velocity, &mut momentum);
    momentum
}
fn mass_action_into(system: &ViscoElements, velocity: &[f64], momentum: &mut [f64]) {
    if let Some(mass) = system.full_reference_mass() {
        mass.spmv(velocity, momentum);
    } else {
        // Legacy stationary supports need only the existing reduced matrix.
        // These zero support rows are never published as full reactions.
        let free = &system.equilibrium.free_dofs;
        let reduced: Vec<_> = free.iter().map(|&d| velocity[d]).collect();
        let mut result = vec![0.0; free.len()];
        system.equilibrium.mass.spmv(&reduced, &mut result);
        momentum.fill(0.0);
        for (&d, value) in free.iter().zip(result) {
            momentum[d] = value;
        }
    }
}
fn kinetic(system: &ViscoElements, velocity: &[f64]) -> Result<f64, ViscoElementError> {
    let momentum = mass_action(system, velocity);
    let value = finite(
        0.5 * velocity
            .iter()
            .zip(momentum)
            .map(|(v, p)| v * p)
            .sum::<f64>(),
    )?;
    if value < 0.0 {
        return Err(invalid("negative consistent-mass kinetic energy"));
    }
    Ok(value)
}

impl ViscoElements {
    /// Prepare energy-consistent dynamics on the existing spatial operators.
    /// All material bands must include DC. The caller still admits the actual
    /// transient loading spectrum and refines time to resolve wave/relaxation
    /// dynamics; a discrete energy identity is not an accuracy certificate.
    pub fn dynamic_operator(
        &self,
        dt_s: f64,
        config: ViscoDynamicConfig,
        cx: &Cx<'_>,
    ) -> Result<ViscoDynamicOperator<'_>, ViscoElementError> {
        checkpoint(cx)?;
        if !dt_s.is_finite()
            || dt_s <= 0.0
            || config.maximum_steps == 0
            || config.maximum_free_dofs == 0
            || self.equilibrium.free_dofs.len() > config.maximum_free_dofs
            || config.solve.maximum_iterations == 0
            || !config.solve.relative_residual_tolerance.is_finite()
            || !(0.0..1.0).contains(&config.solve.relative_residual_tolerance)
            || config.solve.relative_residual_tolerance == 0.0
            || !config.absolute_energy_tolerance_j.is_finite()
            || config.absolute_energy_tolerance_j < 0.0
            || !config.relative_energy_tolerance.is_finite()
            || config.relative_energy_tolerance <= 0.0
            || config.relative_energy_tolerance >= 1.0
        {
            return Err(invalid(
                "finite positive duration and bounded dynamic solve/energy controls required",
            ));
        }
        let coefficient = finite((0.5 * dt_s) * dt_s)?;
        if coefficient == 0.0 {
            return Err(invalid("dynamic duration squared is unrepresentable"));
        }
        let mut tangents = Vec::with_capacity(self.elements.len());
        for element in &self.elements {
            checkpoint(cx)?;
            let _ = element.law.complex_modulus(0.0)?;
            tangents.push(
                element
                    .law
                    .mean_tangent_mandel(dt_s)
                    .map_err(ViscoElementError::Projection)?,
            );
        }
        let (reduced, stiffness) = assemble_tangents(self, &tangents, cx)?;
        let mass = &self.equilibrium.mass;
        let n = mass.nrows();
        let mut effective = Coo::new(n, n);
        for i in 0..n {
            checkpoint(cx)?;
            for (matrix, scale) in [(mass, 1.0), (&stiffness, coefficient)] {
                let (columns, values) = matrix.row(i);
                for (&j, &value) in columns.iter().zip(values) {
                    effective.push(i, j, finite(scale * value)?);
                }
            }
        }
        let effective = effective.assemble();
        if (0..n).any(|i| {
            effective.get(i, i) <= 0.0 || effective.row(i).1.iter().any(|v| !v.is_finite())
        }) {
            return Err(invalid("unrepresentable positive dynamic operator"));
        }
        let preconditioner = Jacobi::new(&effective);
        checkpoint(cx)?;
        Ok(ViscoDynamicOperator {
            system: self,
            dt_s,
            config,
            reduced,
            mean_tangents: tangents,
            operator: CsrOp::symmetric(effective),
            preconditioner,
        })
    }
}

impl<'a> ViscoDynamicOperator<'a> {
    /// Immutable original solid whose tensors, quadrature and inertia are used.
    pub const fn system(&self) -> &'a ViscoElements {
        self.system
    }

    /// Fixed positive duration represented by this prepared operator [s].
    pub const fn duration_s(&self) -> f64 {
        self.dt_s
    }
    /// Consistent mass plus h²/2 times the exact work-mean material tangent.
    pub fn effective_mass(&self) -> &Csr {
        self.operator.matrix()
    }

    /// Continue any accepted original material history, including a quasistatic
    /// preload. Full node-major velocity must be exactly zero on fixed DOFs.
    /// Ordinary `advance` holds their displacement; explicit support-motion
    /// calls may subsequently prescribe motion using retained full mass.
    pub fn start<'p>(
        &'p self,
        material: ViscoElementState<'a>,
        velocity_mps: &[f64],
        cx: &Cx<'_>,
    ) -> Result<ViscoDynamicState<'p, 'a>, ViscoElementError> {
        self.start_impl(material, velocity_mps, false, cx)
    }

    /// Start with explicit physical velocity on supports as well as free nodes.
    /// Requires `ViscoElements::new_with_full_reference_mass` when constrained.
    /// Subsequent support endpoint velocities are supplied to
    /// `advance_with_support_velocity`; no inertial cross term is discarded.
    pub fn start_with_support_motion<'p>(
        &'p self,
        material: ViscoElementState<'a>,
        velocity_mps: &[f64],
        cx: &Cx<'_>,
    ) -> Result<ViscoDynamicState<'p, 'a>, ViscoElementError> {
        self.start_impl(material, velocity_mps, true, cx)
    }

    fn start_impl<'p>(
        &'p self,
        material: ViscoElementState<'a>,
        velocity_mps: &[f64],
        moving_supports: bool,
        cx: &Cx<'_>,
    ) -> Result<ViscoDynamicState<'p, 'a>, ViscoElementError> {
        checkpoint(cx)?;
        if !core::ptr::eq(self.system, material.system)
            || velocity_mps.len() != self.system.full_dofs
            || (moving_supports && self.system.full_reference_mass().is_none())
            || velocity_mps.iter().enumerate().any(|(d, v)| {
                !v.is_finite() || (!moving_supports && self.reduced[d].is_none() && *v != 0.0)
            })
        {
            return Err(invalid(
                "same-system history, finite velocity and full reference mass for moving supports required",
            ));
        }
        let kinetic_j = kinetic(self.system, velocity_mps)?;
        checkpoint(cx)?;
        Ok(ViscoDynamicState {
            prepared: self,
            material,
            velocity_mps: velocity_mps.to_vec(),
            kinetic_j,
            external_work_j: 0.0,
            support_work_j: 0.0,
            eigenstrain_work_j: 0.0,
            parameter_work_j: 0.0,
            time_s: 0.0,
            steps: 0,
        })
    }
}

impl<'p, 'a> ViscoDynamicState<'p, 'a> {
    /// Original fixed-duration sparse operator; no material memory is rebuilt.
    pub const fn operator(&self) -> &'p ViscoDynamicOperator<'a> {
        self.prepared
    }

    /// Original accepted quadrature history and full nodal displacement.
    pub fn material_history(&self) -> &ViscoElementState<'_> {
        &self.material
    }
    /// Full node-major physical velocity, including prescribed supports [m/s].
    pub fn velocity_mps(&self) -> &[f64] {
        &self.velocity_mps
    }
    /// Time since starting this inertial trajectory [s].
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    /// Accepted intervals charged to the original trajectory budget.
    pub const fn steps(&self) -> usize {
        self.steps
    }
    /// Current kinetic energy from the original consistent mass [J].
    pub const fn kinetic_energy_j(&self) -> f64 {
        self.kinetic_j
    }
    /// Cumulative applied-load work since starting this trajectory [J].
    pub const fn external_work_j(&self) -> f64 {
        self.external_work_j
    }

    /// Cumulative signed support work, separate from applied loads and heat [J].
    pub const fn support_work_j(&self) -> f64 {
        self.support_work_j
    }

    /// Cumulative signed imposed-eigenstrain work since starting this trajectory
    /// [J]. This work is separate from applied-load work and physical heat.
    pub const fn eigenstrain_work_j(&self) -> f64 {
        self.eigenstrain_work_j
    }

    /// Signed coefficient storage work since this trajectory started [J].
    /// It excludes material changes applied before `start` and is separate
    /// from applied-load work, imposed-eigenstrain work and physical heat.
    pub const fn parameter_work_j(&self) -> f64 {
        self.parameter_work_j
    }

    /// Advance one interval using the actual interval-mean external nodal load.
    /// For an affine force history this is its endpoint average. No material
    /// endpoint-force trapezoid is used. All fixed supports stay stationary,
    /// so their reaction work is zero. Nonzero accepted support velocity needs
    /// explicit support-motion input instead. Any accepted free strain is held.
    /// Failure preserves the entire state.
    pub fn advance(
        &mut self,
        mean_force_n: &[f64],
        cx: &Cx<'_>,
    ) -> Result<ViscoDynamicStep, ViscoElementError> {
        self.advance_impl(mean_force_n, None, None, None, cx)
            .map(|r| r.0)
    }

    /// Advance the original inertial trajectory while imposing source-integrated
    /// thermal strain at elements or original quadrature points. The mean
    /// tangent, consistent reference mass,
    /// source-domain checks and original quadrature memories are unchanged.
    /// Free strain interpolates affinely over the interval and its exact signed
    /// work enters the existing energy gate separately from nodal loading.
    /// This does not transport constitutive coefficients or solve caloric heat.
    pub fn advance_with_thermal_strain<'t>(
        &mut self,
        mean_force_n: &[f64],
        thermal: impl Into<super::super::QuadraticThermalStrainField<'t>>,
        cx: &Cx<'_>,
    ) -> Result<ViscoDynamicStep, ViscoElementError> {
        self.advance_impl(mean_force_n, None, Some(thermal.into()), None, cx)
            .map(|r| r.0)
    }

    /// Advance with prescribed support endpoint velocities [m/s]. The full
    /// node-major input must be zero on free DOFs. Support displacement follows
    /// `u1-u0 = h*(v0+v1)/2`, so the same midpoint kinematics closes momentum
    /// and reaction work. Unrepresentable support increments refuse. The
    /// original full reference mass must have been retained at preparation.
    pub fn advance_with_support_velocity(
        &mut self,
        mean_force_n: &[f64],
        support_velocity_mps: &[f64],
        cx: &Cx<'_>,
    ) -> Result<ViscoDynamicStep, ViscoElementError> {
        self.advance_impl(mean_force_n, Some(support_velocity_mps), None, None, cx)
            .map(|r| r.0)
    }

    /// The same imposed support motion combined with original thermal strain.
    /// Nodal, support and eigenstrain work remain separate signed inputs.
    pub fn advance_with_support_velocity_and_thermal_strain<'t>(
        &mut self,
        mean_force_n: &[f64],
        support_velocity_mps: &[f64],
        thermal: impl Into<super::super::QuadraticThermalStrainField<'t>>,
        cx: &Cx<'_>,
    ) -> Result<ViscoDynamicStep, ViscoElementError> {
        self.advance_impl(
            mean_force_n,
            Some(support_velocity_mps),
            Some(thermal.into()),
            None,
            cx,
        )
        .map(|r| r.0)
    }

    fn advance_impl(
        &mut self,
        mean_force_n: &[f64],
        support_velocity_mps: Option<&[f64]>,
        thermal: Option<super::super::QuadraticThermalStrainField<'_>>,
        port: Option<(&ViscoDynamicPortOperator<'_, '_>, &[f64])>,
        cx: &Cx<'_>,
    ) -> Result<(ViscoDynamicStep, Option<ViscoDynamicPortResponse>), ViscoElementError> {
        checkpoint(cx)?;
        let prepared = self.prepared;
        let system = prepared.system;
        let h = prepared.dt_s;
        if let Some((port, history)) = port {
            port.admit(prepared, history, cx)?;
        }
        let time = self.time_s + h;
        if self.steps >= prepared.config.maximum_steps
            || !time.is_finite()
            || time <= self.time_s
            || mean_force_n.len() != system.full_dofs
            || mean_force_n.iter().any(|v| !v.is_finite())
        {
            return Err(invalid(
                "dynamic step budget, representable time and complete finite mean load required",
            ));
        }
        if let Some(support) = support_velocity_mps {
            if system.full_reference_mass().is_none()
                || support.len() != system.full_dofs
                || support
                    .iter()
                    .enumerate()
                    .any(|(d, v)| !v.is_finite() || (prepared.reduced[d].is_some() && *v != 0.0))
            {
                return Err(invalid(
                    "moving supports require full reference mass and finite fixed-DOF-only endpoint velocity",
                ));
            }
        } else if self
            .velocity_mps
            .iter()
            .enumerate()
            .any(|(d, v)| prepared.reduced[d].is_none() && *v != 0.0)
        {
            return Err(invalid(
                "nonzero support velocity requires explicit support motion",
            ));
        }
        let thermal = thermal
            .map(|field| self.material.prepare_thermal_field(field, cx))
            .transpose()?;
        let free = &system.equilibrium.free_dofs;
        let mut displacement = self.material.displacement_m.clone();
        let mut next_velocity = self.velocity_mps.clone();
        let mut momentum_velocity = self.velocity_mps.clone();
        let mut support_increment = support_velocity_mps.map(|_| vec![0.0; system.full_dofs]);
        if let (Some(support), Some(increments)) =
            (support_velocity_mps, support_increment.as_mut())
        {
            for d in 0..system.full_dofs {
                if prepared.reduced[d].is_none() {
                    let mean = 0.5 * self.velocity_mps[d] + 0.5 * support[d];
                    let requested_increment = finite(h * mean)?;
                    displacement[d] = finite(displacement[d] + requested_increment)?;
                    let increment = finite(displacement[d] - self.material.displacement_m[d])?;
                    if mean != 0.0 && (requested_increment == 0.0 || increment == 0.0) {
                        return Err(invalid("support midpoint displacement is unrepresentable"));
                    }
                    increments[d] = increment;
                    next_velocity[d] = support[d];
                    // M_ff v_f0 - M_fc (v_c1-v_c0)/2. Applying the native
                    // full mass also retains every off-diagonal inertial term.
                    momentum_velocity[d] = 0.5 * self.velocity_mps[d] - 0.5 * support[d];
                }
            }
        }
        let mut held = {
            let mut trial = self.material.clone();
            trial.step_internal_force(self.material.displacement_m(), h, cx)?
        };
        self.material.add_thermal_mean_force(
            thermal.as_deref(),
            &prepared.mean_tangents,
            &prepared.reduced,
            &mut held.reduced_mean_force_n,
            cx,
        )?;
        if let Some(increment) = &support_increment {
            // Add K_mean,fc delta_u_c with the original B/B^T pair. This
            // affine offset is not an accepted support-only intermediate
            // strain: the final complete motion must pass source admission.
            for (element, tangent) in system.elements.iter().zip(&prepared.mean_tangents) {
                for point in &element.points {
                    checkpoint(cx)?;
                    let strain = point.strain(&element.dofs, increment);
                    let stress: [f64; 6] =
                        core::array::from_fn(|i| (0..6).map(|j| tangent[i][j] * strain[j]).sum());
                    point.for_each_force(stress, |local, contribution| {
                        if let Some(row) = prepared.reduced[element.dofs[local]] {
                            held.reduced_mean_force_n[row] += contribution;
                        }
                    });
                }
            }
        }
        let momentum = mass_action(system, &momentum_velocity);
        // A prescribed support impulse can cancel the previous free momentum.
        // Retain its uncancelled scale when checking the physical equation.
        let support_momentum_scale = match (support_velocity_mps, system.full_reference_mass()) {
            (Some(_), Some(mass)) => euclidean_norm(free.iter().map(|&d| {
                let (columns, values) = mass.row(d);
                columns
                    .iter()
                    .zip(values)
                    .filter(|(column, _)| prepared.reduced[**column].is_none())
                    .map(|(&column, &value)| value * momentum_velocity[column])
                    .sum::<f64>()
            })),
            _ => 0.0,
        };
        let rhs: Vec<_> = free
            .iter()
            .enumerate()
            .map(|(i, &d)| momentum[d] + 0.5 * h * (mean_force_n[d] - held.reduced_mean_force_n[i]))
            .collect();
        let scale = euclidean_norm(rhs.iter().copied())
            .max(euclidean_norm(free.iter().map(|&d| momentum[d])))
            .max(support_momentum_scale)
            .max(0.5 * h * euclidean_norm(free.iter().map(|&d| mean_force_n[d])))
            .max(0.5 * h * euclidean_norm(held.reduced_mean_force_n.iter().copied()))
            .max(f64::MIN_POSITIVE);
        if !scale.is_finite() || rhs.iter().any(|v| !v.is_finite()) {
            return Err(ViscoElementError::NonFiniteResult);
        }
        let (mut mean_velocity, iterations) = solve_linear(
            &prepared.operator,
            &prepared.preconditioner,
            &rhs,
            scale,
            prepared.config.solve,
            cx,
        )?;
        if let Some((port, history)) = port {
            port.correct(&mut mean_velocity, history, cx)?;
        }
        for (i, &d) in free.iter().enumerate() {
            displacement[d] = finite(displacement[d] + h * mean_velocity[i])?;
            // Recover the actual representable increment, ensuring constitutive
            // work and midpoint inertia use the same kinematics.
            let actual_mean = (displacement[d] - self.material.displacement_m[d]) / h;
            next_velocity[d] = finite(2.0 * actual_mean - self.velocity_mps[d])?;
        }
        // The source law sees the actual representable displacement increment.
        // Reconstruct its force before the original momentum, reaction and
        // material-work checks; a predicted port effort never bypasses them.
        let port_result = port
            .map(|(port, history)| {
                port.response(
                    mean_force_n,
                    self.material.displacement_m(),
                    &displacement,
                    history,
                    cx,
                )
            })
            .transpose()?;
        let (mean_force_n, port_response, scale) = match port_result.as_ref() {
            Some((net_force, response)) => {
                let boundary_scale =
                    euclidean_norm(free.iter().map(|&d| mean_force_n[d] - net_force[d]));
                let net_scale = euclidean_norm(free.iter().map(|&d| net_force[d]));
                (
                    net_force.as_slice(),
                    Some(response),
                    finite(scale.max(0.5 * h * boundary_scale).max(0.5 * h * net_scale))?,
                )
            }
            None => (mean_force_n, None, scale),
        };
        let mut candidate = self.material.clone();
        let constitutive =
            candidate.step_internal_force_impl(&displacement, h, thermal.as_ref(), cx)?;
        let delta_velocity: Vec<_> = next_velocity
            .iter()
            .zip(&self.velocity_mps)
            .map(|(new, old)| new - old)
            .collect();
        let impulse = mass_action(system, &delta_velocity);
        let residual = euclidean_norm(free.iter().enumerate().map(|(i, &d)| {
            impulse[d] + h * (constitutive.reduced_mean_force_n[i] - mean_force_n[d])
        }));
        check_residual(residual / (2.0 * scale), iterations, prepared.config.solve)?;
        let force_residual_n = finite(residual / h)?;
        let mean_reaction_force_n = if system.full_reference_mass().is_some() {
            Some(
                (0..system.full_dofs)
                    .map(|d| {
                        finite(impulse[d] / h + constitutive.full_mean_force_n[d] - mean_force_n[d])
                    })
                    .collect::<Result<Vec<_>, _>>()?,
            )
        } else {
            None
        };
        let kinetic_change = increment(
            &self.velocity_mps,
            &next_velocity,
            |velocity, momentum, _| {
                mass_action_into(system, velocity, momentum);
                Ok::<(), ViscoElementError>(())
            },
            &mut || cx.checkpoint().is_err(),
        )
        .map_err(|error| match error {
            QuadraticError::Cancelled => TetElasticError::Cancelled.into(),
            QuadraticError::Operator(error) => error,
            QuadraticError::Dimension | QuadraticError::NonFinite => {
                ViscoElementError::NonFiniteResult
            }
        })?
        .change;
        let kinetic_energy_change_j = finite(kinetic_change.to_f64())?;
        let mut load_work = Dd::ZERO;
        let mut support_work = Dd::ZERO;
        for d in 0..system.full_dofs {
            if d % 256 == 0 {
                checkpoint(cx)?;
            }
            let displacement_change = Dd::from_f64(displacement[d])
                - Dd::from_f64(self.material.displacement_m[d]);
            load_work = load_work + Dd::from_f64(mean_force_n[d]) * displacement_change;
            if prepared.reduced[d].is_none() {
                if let Some(reaction) = &mean_reaction_force_n {
                    support_work = support_work
                        + Dd::from_f64(reaction[d]) * displacement_change;
                }
            }
        }
        let support_work_j = finite(support_work.to_f64())?;
        let kinetic_j = kinetic(system, &next_velocity)?;
        let work = finite(load_work.to_f64())?;
        let mut material_change = Dd::ZERO;
        for element in &constitutive.element_steps {
            checkpoint(cx)?;
            material_change = material_change
                + Dd::from_f64(element.stored_energy_change_j);
        }
        let material_stored_energy_change_j = finite(material_change.to_f64())?;
        let defect = finite(
            (load_work + support_work + Dd::from_f64(constitutive.step_eigenstrain_work_j)
                - kinetic_change - material_change
                - Dd::from_f64(constitutive.step_dissipated_j))
                .to_f64(),
        )?;
        let allowance = finite(
            prepared.config.absolute_energy_tolerance_j
                + prepared.config.relative_energy_tolerance
                    * (work.abs()
                        + support_work_j.abs()
                        + constitutive.step_eigenstrain_work_j.abs()
                        + kinetic_energy_change_j.abs()
                        + material_stored_energy_change_j.abs()
                        + constitutive.step_dissipated_j.abs()),
        )?;
        if defect.abs() > allowance {
            return Err(ViscoElementError::DynamicEnergy {
                residual_j: defect,
                allowance_j: allowance,
            });
        }
        if let Some(port) = port_response {
            let residual_j = finite(port.applied_work_j - port.boundary_work_j - work)?;
            let allowance_j = finite(
                allowance
                    + prepared.config.relative_energy_tolerance
                        * (port.applied_work_j.abs() + port.boundary_work_j.abs()),
            )?;
            if residual_j.abs() > allowance_j {
                return Err(ViscoElementError::DynamicEnergy {
                    residual_j,
                    allowance_j,
                });
            }
        }
        let external_work = finite(self.external_work_j + work)?;
        let support_work = finite(self.support_work_j + support_work_j)?;
        let eigenstrain_work =
            finite(self.eigenstrain_work_j + constitutive.step_eigenstrain_work_j)?;
        checkpoint(cx)?;
        self.material = candidate;
        self.velocity_mps = next_velocity;
        self.kinetic_j = kinetic_j;
        self.external_work_j = external_work;
        self.support_work_j = support_work;
        self.eigenstrain_work_j = eigenstrain_work;
        self.time_s = time;
        self.steps += 1;
        Ok((
            ViscoDynamicStep {
                time_s: time,
                iterations,
                force_residual_n,
                kinetic_energy_j: kinetic_j,
                kinetic_energy_change_j,
                material_stored_energy_change_j,
                external_work_j: work,
                support_work_j,
                mean_reaction_force_n,
                balance_residual_j: defect,
                constitutive,
            },
            port_result.map(|(_, response)| response),
        ))
    }
}
