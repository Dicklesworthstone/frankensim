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
//! Supports remain stationary; the original reduced consistent mass is exact
//! on that subspace. Moving supports need full mass coupling and are absent.

use super::quasistatic::{assemble_tangents, check_residual, solve_linear};
use super::*;
use crate::linear::Jacobi;
use crate::linear3::{TetStaticSolveConfig, euclidean_norm};
use fs_solver::op::CsrOp;

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
    /// Relative allowance against actual work, stored energy and heat.
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
    eigenstrain_work_j: f64,
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
    /// Actual interval-mean applied nodal force times actual displacement [J].
    pub external_work_j: f64,
    /// Exactly zero for this stationary-support owner [J].
    pub support_work_j: f64,
    /// External nodal plus imposed-eigenstrain work minus kinetic/storage
    /// changes and physical heat [J].
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
fn kinetic(mass: &Csr, velocity: &[f64]) -> Result<f64, ViscoElementError> {
    let mut momentum = vec![0.0; velocity.len()];
    mass.spmv(velocity, &mut momentum);
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
    /// Their accepted displacement remains constant for this entire trajectory.
    pub fn start<'p>(
        &'p self,
        material: ViscoElementState<'a>,
        velocity_mps: &[f64],
        cx: &Cx<'_>,
    ) -> Result<ViscoDynamicState<'p, 'a>, ViscoElementError> {
        checkpoint(cx)?;
        if !core::ptr::eq(self.system, material.system)
            || velocity_mps.len() != self.system.full_dofs
            || velocity_mps
                .iter()
                .enumerate()
                .any(|(d, v)| !v.is_finite() || (self.reduced[d].is_none() && *v != 0.0))
        {
            return Err(invalid(
                "same-system material history and stationary-support finite velocity required",
            ));
        }
        let velocity: Vec<_> = self
            .system
            .equilibrium
            .free_dofs
            .iter()
            .map(|&d| velocity_mps[d])
            .collect();
        let kinetic_j = kinetic(&self.system.equilibrium.mass, &velocity)?;
        checkpoint(cx)?;
        Ok(ViscoDynamicState {
            prepared: self,
            material,
            velocity_mps: velocity_mps.to_vec(),
            kinetic_j,
            external_work_j: 0.0,
            eigenstrain_work_j: 0.0,
            time_s: 0.0,
            steps: 0,
        })
    }
}

impl ViscoDynamicState<'_, '_> {
    /// Original accepted quadrature history and full nodal displacement.
    pub fn material_history(&self) -> &ViscoElementState<'_> {
        &self.material
    }
    /// Full node-major physical velocity; fixed support entries are zero [m/s].
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

    /// Cumulative signed imposed-eigenstrain work since starting this trajectory
    /// [J]. This work is separate from applied-load work and physical heat.
    pub const fn eigenstrain_work_j(&self) -> f64 {
        self.eigenstrain_work_j
    }

    /// Advance one interval using the actual interval-mean external nodal load.
    /// For an affine force history this is its endpoint average. No material
    /// endpoint-force trapezoid is used. All fixed supports stay stationary,
    /// so their reaction work is zero; full dynamic support reactions are not
    /// reconstructed from the reduced mass. Any accepted free strain is held.
    /// Failure preserves the entire state.
    pub fn advance(
        &mut self,
        mean_force_n: &[f64],
        cx: &Cx<'_>,
    ) -> Result<ViscoDynamicStep, ViscoElementError> {
        self.advance_impl(mean_force_n, None, cx)
    }

    /// Advance the original inertial trajectory while imposing source-integrated
    /// elementwise thermal strain. The mean tangent, consistent reference mass,
    /// source-domain checks and original quadrature memories are unchanged.
    /// Free strain interpolates affinely over the interval and its exact signed
    /// work enters the existing energy gate separately from nodal loading.
    /// This does not transport constitutive coefficients or solve caloric heat.
    pub fn advance_with_thermal_strain(
        &mut self,
        mean_force_n: &[f64],
        thermal: super::super::TetThermalStrainField<'_>,
        cx: &Cx<'_>,
    ) -> Result<ViscoDynamicStep, ViscoElementError> {
        self.advance_impl(mean_force_n, Some(thermal), cx)
    }

    fn advance_impl(
        &mut self,
        mean_force_n: &[f64],
        thermal: Option<super::super::TetThermalStrainField<'_>>,
        cx: &Cx<'_>,
    ) -> Result<ViscoDynamicStep, ViscoElementError> {
        checkpoint(cx)?;
        let prepared = self.prepared;
        let system = prepared.system;
        let h = prepared.dt_s;
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
        self.material.validate_thermal_field(thermal, cx)?;
        let free = &system.equilibrium.free_dofs;
        let mut held = {
            let mut trial = self.material.clone();
            trial.step_internal_force(self.material.displacement_m(), h, cx)?
        };
        self.material.add_thermal_mean_force(
            thermal, &prepared.mean_tangents, &prepared.reduced,
            &mut held.reduced_mean_force_n, cx,
        )?;
        let velocity: Vec<_> = free.iter().map(|&d| self.velocity_mps[d]).collect();
        let mut momentum = vec![0.0; free.len()];
        system.equilibrium.mass.spmv(&velocity, &mut momentum);
        let rhs: Vec<_> = free
            .iter()
            .enumerate()
            .map(|(i, &d)| momentum[i] + 0.5 * h * (mean_force_n[d] - held.reduced_mean_force_n[i]))
            .collect();
        let scale = euclidean_norm(rhs.iter().copied())
            .max(euclidean_norm(momentum.iter().copied()))
            .max(0.5 * h * euclidean_norm(free.iter().map(|&d| mean_force_n[d])))
            .max(0.5 * h * euclidean_norm(held.reduced_mean_force_n.iter().copied()))
            .max(f64::MIN_POSITIVE);
        if !scale.is_finite() || rhs.iter().any(|v| !v.is_finite()) {
            return Err(ViscoElementError::NonFiniteResult);
        }
        let (mean_velocity, iterations) = solve_linear(
            &prepared.operator,
            &prepared.preconditioner,
            &rhs,
            scale,
            prepared.config.solve,
            cx,
        )?;
        let mut displacement = self.material.displacement_m.clone();
        let mut next_velocity = self.velocity_mps.clone();
        for (i, &d) in free.iter().enumerate() {
            displacement[d] = finite(displacement[d] + h * mean_velocity[i])?;
            // Recover the actual representable increment, ensuring constitutive
            // work and midpoint inertia use the same kinematics.
            let actual_mean = (displacement[d] - self.material.displacement_m[d]) / h;
            next_velocity[d] = finite(2.0 * actual_mean - self.velocity_mps[d])?;
        }
        let mut candidate = self.material.clone();
        let constitutive = candidate.step_internal_force_impl(&displacement, h, thermal, cx)?;
        let delta_velocity: Vec<_> = free
            .iter()
            .map(|&d| next_velocity[d] - self.velocity_mps[d])
            .collect();
        let mut impulse = vec![0.0; free.len()];
        system.equilibrium.mass.spmv(&delta_velocity, &mut impulse);
        let residual = euclidean_norm(free.iter().enumerate().map(|(i, &d)| {
            impulse[i] + h * (constitutive.reduced_mean_force_n[i] - mean_force_n[d])
        }));
        check_residual(residual / (2.0 * scale), iterations, prepared.config.solve)?;
        let force_residual_n = finite(residual / h)?;
        let new_free_velocity: Vec<_> = free.iter().map(|&d| next_velocity[d]).collect();
        let kinetic_j = kinetic(&system.equilibrium.mass, &new_free_velocity)?;
        let work = finite(
            mean_force_n
                .iter()
                .zip(displacement.iter().zip(&self.material.displacement_m))
                .map(|(f, (new, old))| f * (new - old))
                .sum(),
        )?;
        let defect = finite(
            work + constitutive.step_eigenstrain_work_j - (kinetic_j - self.kinetic_j)
                - (constitutive.energy.stored_j - self.material.energy.stored_j)
                - constitutive.step_dissipated_j,
        )?;
        let allowance = finite(
            prepared.config.absolute_energy_tolerance_j
                + prepared.config.relative_energy_tolerance
                    * (work.abs()
                        + constitutive.step_eigenstrain_work_j.abs()
                        + kinetic_j
                        + self.kinetic_j
                        + constitutive.energy.stored_j.abs()
                        + self.material.energy.stored_j.abs()
                        + constitutive.step_dissipated_j.abs()),
        )?;
        if defect.abs() > allowance {
            return Err(ViscoElementError::DynamicEnergy {
                residual_j: defect,
                allowance_j: allowance,
            });
        }
        let external_work = finite(self.external_work_j + work)?;
        let eigenstrain_work = finite(
            self.eigenstrain_work_j + constitutive.step_eigenstrain_work_j,
        )?;
        checkpoint(cx)?;
        self.material = candidate;
        self.velocity_mps = next_velocity;
        self.kinetic_j = kinetic_j;
        self.external_work_j = external_work;
        self.eigenstrain_work_j = eigenstrain_work;
        self.time_s = time;
        self.steps += 1;
        Ok(ViscoDynamicStep {
            time_s: time,
            iterations,
            force_residual_n,
            kinetic_energy_j: kinetic_j,
            external_work_j: work,
            support_work_j: 0.0,
            balance_residual_j: defect,
            constitutive,
        })
    }
}
