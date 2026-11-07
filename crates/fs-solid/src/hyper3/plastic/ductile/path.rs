//! Adaptive physical intervals on the original finite ductile endpoint owners.
//!
//! One controller serves quasistatic and moving specimens. Each proposal uses
//! the same immutable last accepted configuration, velocity and plastic/damage
//! history. Supports interpolate from their actual accepted coordinates; the
//! caller supplies both endpoints of dead nodal and referential body loads.
//! The accepted fraction supplies the physical duration to BOTH initiation
//! and native momentum. Temperatures stay held for the complete interval.
//!
//! Damage, summed absolute material work defects and kinetic numerical loss
//! control subdivision with separate whole-interval budgets. These gates are
//! not a trajectory-error estimator, arc-length continuation through a force
//! limit point, a fracture-topology model or an experimental failure claim.
//! Loading reversals require separate intervals. Only the last full endpoint
//! and bounded scalar records survive; a failed request publishes no prefix.

use super::dynamics::{
    FiniteDuctileMotionControls, FiniteDuctileMotionStep, FiniteDuctileTetMotion,
};
use super::*;

/// Bounded subdivision and whole-interval physical-work allowances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteDuctilePathControls {
    pub initial_fraction: f64,
    /// Smallest attempted fraction, except one final shorter remainder.
    pub minimum_fraction: f64,
    pub maximum_fraction: f64,
    /// Bounds the accepted records and native endpoint advances.
    pub maximum_accepted_steps: usize,
    /// Bounds all native endpoint solves, including refused attempts.
    pub maximum_trials: usize,
    /// Maximum damage increase of any element in an accepted subinterval.
    /// The original material's own, possibly stricter, limit remains active.
    pub maximum_damage_increment: f64,
    /// Quasistatic initial free-force infinity norm [N]. Moving specimens can
    /// start with physical acceleration and do not impose this equilibrium.
    pub initial_equilibrium_tolerance_n: f64,
    /// The material defect and equilibrium-work tolerances apply to the WHOLE
    /// interval. Each trial receives its actual fraction of these allowances;
    /// the accepted yield tolerance is unchanged.
    pub work: FinitePlasticStepControls,
}

/// Independent whole-interval kinetic-loss allowance for moving specimens.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteDuctileMotionPathControls {
    pub path: FiniteDuctilePathControls,
    pub maximum_numerical_kinetic_loss_j: f64,
}

/// Original material, optional native momentum and complete physical balance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteDuctilePathEnergy {
    pub material: FiniteDuctileTetEnergy,
    /// Absent for quasistatic paths. Kinetic change is the compensated sum of
    /// native work-conjugate increments, avoiding subtraction of large K.
    pub inertia: Option<InertialEnergy>,
    pub balance_error_j: f64,
}

/// One accepted physical endpoint, without copying its full fields.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteDuctilePathRecord {
    pub fraction: f64,
    pub duration_s: f64,
    pub maximum_damage: f64,
    pub maximum_damage_increment: f64,
    pub energy: FiniteDuctilePathEnergy,
}

/// Complete interval. The original native endpoint and its element receipts
/// describe the LAST substep; `energy` and `records` cover the whole interval.
#[derive(Debug)]
pub struct FiniteDuctilePathStep<S> {
    pub endpoint: S,
    pub records: Vec<FiniteDuctilePathRecord>,
    pub trials: usize,
    pub rejected_trials: usize,
    pub duration_s: f64,
    pub energy: FiniteDuctilePathEnergy,
}

#[derive(Debug)]
pub enum FiniteDuctileRefinementCause {
    Endpoint(FiniteDuctileTetError),
    DamageIncrement { observed: f64, maximum: f64 },
}

/// A refusal never returns an internally accepted partial interval.
#[derive(Debug)]
pub enum FiniteDuctilePathError {
    Invalid(&'static str),
    Allocation,
    Endpoint(FiniteDuctileTetError),
    InitialEquilibrium {
        maximum_residual_n: f64,
        tolerance_n: f64,
    },
    WorkLimit {
        accepted_steps: usize,
        trials: usize,
    },
    MinimumInterval {
        accepted_fraction: f64,
        attempted_fraction: f64,
        cause: FiniteDuctileRefinementCause,
    },
}
impl core::fmt::Display for FiniteDuctilePathError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite ductile physical path refused: {self:?}")
    }
}
impl std::error::Error for FiniteDuctilePathError {}
impl From<FiniteDuctileTetError> for FiniteDuctilePathError {
    fn from(error: FiniteDuctileTetError) -> Self {
        Self::Endpoint(error)
    }
}
impl From<SmallStrainTetError> for FiniteDuctilePathError {
    fn from(error: SmallStrainTetError) -> Self {
        Self::Endpoint(error.into())
    }
}

impl<'g, 'm, L> FiniteDuctileTetSolid<'g, 'm, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    /// Resolve one affine physical interval through the original `advance`.
    /// `prescribed_end_m` fixes the support set and final absolute coordinates;
    /// starts are read from this accepted specimen. Dead nodal loads are
    /// `(dof, [initial, final])`; body loads are `[initial_vector, final_vector]`
    /// per original reference volume. Initial free equations must already be
    /// in equilibrium, including the declared initial body and nodal loads.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_path(
        &self,
        prescribed_end_m: &[(usize, f64)],
        force_intervals_n: &[(usize, [f64; 2])],
        body_force_intervals_n_m3: [[f64; 3]; 2],
        conditions: DuctileStepConditions<'_>,
        settings: IncrementSettings,
        controls: FiniteDuctilePathControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteDuctilePathStep<FiniteDuctileTetStep<'g, 'm, L>>, FiniteDuctilePathError>
    {
        advance_path(
            self,
            prescribed_end_m,
            force_intervals_n,
            body_force_intervals_n_m3,
            conditions,
            settings,
            controls,
            None,
            cx,
        )
    }
}

impl<'g, 'm, L> FiniteDuctileTetMotion<'g, 'm, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    /// Resolve an affine load/support interval using the original consistent
    /// reference mass and backward-Euler momentum at every accepted duration.
    /// Nonzero incoming velocity is retained. A prescribed velocity jump can
    /// exhaust its kinetic-loss budget; subdivision does not supply an impulse.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_path(
        &self,
        prescribed_end_m: &[(usize, f64)],
        force_intervals_n: &[(usize, [f64; 2])],
        body_force_intervals_n_m3: [[f64; 3]; 2],
        conditions: DuctileStepConditions<'_>,
        settings: IncrementSettings,
        controls: FiniteDuctileMotionPathControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteDuctilePathStep<FiniteDuctileMotionStep<'g, 'm, L>>, FiniteDuctilePathError>
    {
        advance_path(
            self,
            prescribed_end_m,
            force_intervals_n,
            body_force_intervals_n_m3,
            conditions,
            settings,
            controls.path,
            Some(controls.maximum_numerical_kinetic_loss_j),
            cx,
        )
    }
}

// Private adapters select the existing endpoint owner; they own no mechanics,
// return mapping, inertia, sparse matrix, Newton iteration or heat partition.
trait PathBody<'g, 'm, L>: Sized
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    type Step;
    fn solid(&self) -> &FiniteDuctileTetSolid<'g, 'm, L>;
    fn body(step: &Self::Step) -> &Self;
    fn energy(step: &Self::Step) -> FiniteDuctilePathEnergy;
    #[allow(clippy::too_many_arguments)]
    fn endpoint(
        &self,
        supports: &[(usize, f64)],
        forces: &[(usize, f64)],
        body: [f64; 3],
        conditions: DuctileStepConditions<'_>,
        settings: IncrementSettings,
        work: FinitePlasticStepControls,
        kinetic_limit: Option<f64>,
        cx: &Cx<'_>,
    ) -> Result<Self::Step, FiniteDuctileTetError>;
}
impl<'g, 'm, L> PathBody<'g, 'm, L> for FiniteDuctileTetSolid<'g, 'm, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    type Step = FiniteDuctileTetStep<'g, 'm, L>;
    fn solid(&self) -> &Self {
        self
    }
    fn body(step: &Self::Step) -> &Self {
        &step.solid
    }
    fn energy(step: &Self::Step) -> FiniteDuctilePathEnergy {
        FiniteDuctilePathEnergy {
            material: step.energy,
            inertia: None,
            balance_error_j: step.energy.balance_error_j,
        }
    }
    fn endpoint(
        &self,
        supports: &[(usize, f64)],
        forces: &[(usize, f64)],
        body: [f64; 3],
        conditions: DuctileStepConditions<'_>,
        settings: IncrementSettings,
        work: FinitePlasticStepControls,
        _: Option<f64>,
        cx: &Cx<'_>,
    ) -> Result<Self::Step, FiniteDuctileTetError> {
        self.advance(supports, forces, body, conditions, settings, work, cx)
    }
}
impl<'g, 'm, L> PathBody<'g, 'm, L> for FiniteDuctileTetMotion<'g, 'm, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    type Step = FiniteDuctileMotionStep<'g, 'm, L>;
    fn solid(&self) -> &FiniteDuctileTetSolid<'g, 'm, L> {
        self.solid()
    }
    fn body(step: &Self::Step) -> &Self {
        &step.motion
    }
    fn energy(step: &Self::Step) -> FiniteDuctilePathEnergy {
        FiniteDuctilePathEnergy {
            material: step.material_energy,
            inertia: Some(step.inertia),
            balance_error_j: step.balance_error_j,
        }
    }
    fn endpoint(
        &self,
        supports: &[(usize, f64)],
        forces: &[(usize, f64)],
        body: [f64; 3],
        conditions: DuctileStepConditions<'_>,
        settings: IncrementSettings,
        work: FinitePlasticStepControls,
        kinetic_limit: Option<f64>,
        cx: &Cx<'_>,
    ) -> Result<Self::Step, FiniteDuctileTetError> {
        self.advance(
            supports,
            forces,
            body,
            conditions,
            settings,
            FiniteDuctileMotionControls {
                material: work,
                maximum_numerical_kinetic_loss_j: kinetic_limit
                    .expect("moving path has kinetic allowance"),
            },
            cx,
        )
    }
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn advance_path<'g, 'm, L, B>(
    initial: &B,
    prescribed_end_m: &[(usize, f64)],
    force_intervals_n: &[(usize, [f64; 2])],
    body_intervals: [[f64; 3]; 2],
    conditions: DuctileStepConditions<'_>,
    settings: IncrementSettings,
    controls: FiniteDuctilePathControls,
    kinetic_limit: Option<f64>,
    cx: &Cx<'_>,
) -> Result<FiniteDuctilePathStep<B::Step>, FiniteDuctilePathError>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq + 'm,
    B: PathBody<'g, 'm, L>,
{
    checkpoint(cx)?;
    validate_controls(controls, kinetic_limit)?;
    let specimen = initial.solid();
    let old_u = specimen.displacement_m().as_flattened();
    let covered = match conditions.temperatures_k {
        DuctileTemperatures::Uniform(_) => true,
        DuctileTemperatures::PerPoint(values) => values.len() == specimen.materials.len(),
        DuctileTemperatures::PerHexElement(_) => false,
    };
    if !covered
        || !conditions.duration_s.is_finite()
        || conditions.duration_s <= 0.0
        || body_intervals.as_flattened().iter().any(|v| !v.is_finite())
    {
        return Err(FiniteDuctilePathError::Invalid(
            "complete held tetrahedral temperatures, positive duration and finite body-load endpoints required",
        ));
    }
    let mut supports = buffer(prescribed_end_m.len())?;
    for &(dof, end) in prescribed_end_m {
        checkpoint(cx)?;
        let Some(&start) = old_u.get(dof) else {
            return Err(FiniteDuctilePathError::Invalid(
                "prescribed DOF is outside original mesh",
            ));
        };
        if !end.is_finite() {
            return Err(FiniteDuctilePathError::Invalid(
                "finite prescribed endpoint required",
            ));
        }
        supports.push((dof, start));
    }
    let mut forces = buffer(force_intervals_n.len())?;
    for &(dof, endpoints) in force_intervals_n {
        checkpoint(cx)?;
        if dof >= old_u.len() || endpoints.iter().any(|v| !v.is_finite()) {
            return Err(FiniteDuctilePathError::Invalid(
                "finite force endpoints on original DOFs required",
            ));
        }
        forces.push((dof, endpoints[0]));
    }
    if kinetic_limit.is_none() {
        initial_equilibrium(
            specimen,
            &supports,
            &forces,
            body_intervals[0],
            conditions,
            settings,
            controls.initial_equilibrium_tolerance_n,
            cx,
        )?;
    }
    let mut current: Option<B::Step> = None;
    let mut records = buffer(controls.maximum_accepted_steps.min(64))?;
    // Original material losses/work plus native kinetic work/change/loss.
    let mut sums = [Sum::default(); 10];
    let mut progress = 0.0_f64;
    let mut fraction = controls.initial_fraction;
    let mut trials = 0;
    let mut rejected_trials = 0;
    while progress < 1.0 {
        checkpoint(cx)?;
        if trials >= controls.maximum_trials || records.len() >= controls.maximum_accepted_steps {
            return Err(FiniteDuctilePathError::WorkLimit {
                accepted_steps: records.len(),
                trials,
            });
        }
        // Close a roundoff-sized final remainder in this actual interval;
        // otherwise a nearly exact division can invent a vanishing last dt.
        // The resulting full span still supplies every physical/loss budget.
        let target = progress + fraction;
        let end = if 1.0 - target <= 16.0 * f64::EPSILON * target.abs() {
            1.0
        } else {
            target
        };
        let span = end - progress;
        let duration_s = conditions.duration_s * span;
        if end <= progress || !duration_s.is_finite() || duration_s <= 0.0 {
            return Err(FiniteDuctilePathError::Invalid(
                "physical subinterval cannot be represented",
            ));
        }
        for (index, (entry, &(dof, target))) in
            supports.iter_mut().zip(prescribed_end_m).enumerate()
        {
            if index % 256 == 0 {
                checkpoint(cx)?;
            }
            entry.1 = interpolate(old_u[dof], target, end);
        }
        for (index, (entry, &(_, [start, target]))) in
            forces.iter_mut().zip(force_intervals_n).enumerate()
        {
            if index % 256 == 0 {
                checkpoint(cx)?;
            }
            entry.1 = interpolate(start, target, end);
        }
        let body =
            core::array::from_fn(|i| interpolate(body_intervals[0][i], body_intervals[1][i], end));
        let allowance = controls.work.maximum_endpoint_work_defect_j * span;
        let kinetic_allowance = kinetic_limit.map(|limit| limit * span);
        let work = FinitePlasticStepControls {
            equilibrium_work_tolerance_j: controls.work.equilibrium_work_tolerance_j * span,
            maximum_endpoint_work_defect_j: allowance,
            ..controls.work
        };
        let accepted = current.as_ref().map_or(initial, B::body);
        trials += 1;
        let candidate = match accepted.endpoint(
            &supports,
            &forces,
            body,
            DuctileStepConditions {
                duration_s,
                ..conditions
            },
            settings,
            work,
            kinetic_allowance,
            cx,
        ) {
            Ok(candidate) => candidate,
            Err(error) if retryable(&error) => {
                fraction = refine(
                    progress,
                    span,
                    controls.minimum_fraction,
                    FiniteDuctileRefinementCause::Endpoint(error),
                )?;
                rejected_trials += 1;
                continue;
            }
            Err(error) => return Err(error.into()),
        };
        let mut maximum_damage = 0.0_f64;
        let mut maximum_increment = 0.0_f64;
        let mut damage_ratio = 0.0_f64;
        for (element, ((old, new), material)) in accepted
            .solid()
            .element_states()
            .iter()
            .zip(B::body(&candidate).solid().element_states())
            .zip(specimen.materials)
            .enumerate()
        {
            if element % 128 == 0 {
                checkpoint(cx)?;
            }
            let increment = new.damage() - old.damage();
            if !increment.is_finite() || increment < 0.0 {
                return Err(FiniteDuctilePathError::Invalid(
                    "irreversible finite damage increment required",
                ));
            }
            maximum_damage = maximum_damage.max(new.damage());
            maximum_increment = maximum_increment.max(increment);
            let limit = controls
                .maximum_damage_increment
                .min(material.settings().maximum_damage_increment);
            damage_ratio = damage_ratio.max(increment / limit);
        }
        if maximum_increment > controls.maximum_damage_increment {
            fraction = refine(
                progress,
                span,
                controls.minimum_fraction,
                FiniteDuctileRefinementCause::DamageIncrement {
                    observed: maximum_increment,
                    maximum: controls.maximum_damage_increment,
                },
            )?;
            rejected_trials += 1;
            continue;
        }
        let energy = B::energy(&candidate);
        let material = energy.material;
        let inertia = energy.inertia;
        for (sum, value) in sums.iter_mut().zip([
            material.plastic_dissipation_j,
            material.damage_release_j,
            material.dead_load_endpoint_work_j,
            material.support_endpoint_work_j,
            material.internal_endpoint_work_j,
            material.endpoint_work_defect_j,
            material.absolute_endpoint_work_defect_j,
            inertia.map_or(0.0, |e| e.inertial_work_j),
            inertia.map_or(0.0, |e| e.kinetic_change_j),
            inertia.map_or(0.0, |e| e.numerical_loss_j),
        ]) {
            sum.add(value);
            if !sum.value().is_finite() {
                return Err(FiniteDuctilePathError::Invalid(
                    "whole-interval work cannot be represented",
                ));
            }
        }
        if records.len() == records.capacity() {
            let growth = records
                .capacity()
                .max(1)
                .min(controls.maximum_accepted_steps - records.len());
            records
                .try_reserve_exact(growth)
                .map_err(|_| FiniteDuctilePathError::Allocation)?;
        }
        records.push(FiniteDuctilePathRecord {
            fraction: end,
            duration_s,
            maximum_damage,
            maximum_damage_increment: maximum_increment,
            energy,
        });
        current = Some(candidate);
        progress = end;
        let damage_growth = if damage_ratio > 0.0 {
            (0.8 / damage_ratio).clamp(0.5, 2.0)
        } else {
            2.0
        };
        let material_growth = growth(allowance, material.absolute_endpoint_work_defect_j);
        let kinetic_growth = match (kinetic_allowance, inertia) {
            (Some(limit), Some(account)) => growth(limit, account.numerical_loss_j),
            _ => 2.0,
        };
        fraction = (span * damage_growth.min(material_growth).min(kinetic_growth))
            .clamp(controls.minimum_fraction, controls.maximum_fraction);
    }
    let endpoint = current.ok_or(FiniteDuctilePathError::Invalid(
        "no physical endpoint accepted",
    ))?;
    let values = sums.map(Sum::value);
    let stored_before = specimen.stored_energy_j();
    let stored_after = B::body(&endpoint).solid().stored_energy_j();
    let material = FiniteDuctileTetEnergy {
        stored_before_j: stored_before,
        stored_after_j: stored_after,
        plastic_dissipation_j: values[0],
        damage_release_j: values[1],
        dead_load_endpoint_work_j: values[2],
        support_endpoint_work_j: values[3],
        internal_endpoint_work_j: values[4],
        endpoint_work_defect_j: values[5],
        absolute_endpoint_work_defect_j: values[6],
        equilibrium_work_residual_j: (values[2] + values[3]) - values[4] - values[7],
        balance_error_j: ((values[2] + values[3]) - (stored_after - stored_before))
            - values[0]
            - values[1]
            - values[5]
            - values[7],
    };
    let inertia = match (records[0].energy.inertia, B::energy(&endpoint).inertia) {
        (Some(first), Some(last)) => Some(InertialEnergy {
            reference_mass_kg: first.reference_mass_kg,
            kinetic_before_j: first.kinetic_before_j,
            kinetic_after_j: last.kinetic_after_j,
            kinetic_change_j: values[8],
            numerical_loss_j: values[9],
            inertial_work_j: values[7],
            balance_error_j: values[7] - values[8] - values[9],
        }),
        (None, None) => None,
        _ => {
            return Err(FiniteDuctilePathError::Invalid(
                "original momentum ownership changed",
            ));
        }
    };
    let balance_error_j = material.balance_error_j + inertia.map_or(0.0, |e| e.balance_error_j);
    for (name, value) in [
        (
            "path constitutive work",
            values[4] - (stored_after - stored_before) - values[0] - values[1] - values[5],
        ),
        (
            "path support/dead-load equilibrium work",
            material.equilibrium_work_residual_j,
        ),
        (
            "path material physical work balance",
            material.balance_error_j,
        ),
        (
            "path kinetic work identity",
            inertia.map_or(0.0, |e| e.balance_error_j),
        ),
        ("path complete physical work balance", balance_error_j),
    ] {
        gate(name, value, controls.work.equilibrium_work_tolerance_j)?;
    }
    gate(
        "path sum of absolute material integration defects",
        values[6],
        controls.work.maximum_endpoint_work_defect_j,
    )?;
    if let Some(limit) = kinetic_limit {
        gate("path numerical kinetic loss", values[9], limit)?;
    }
    checkpoint(cx)?;
    Ok(FiniteDuctilePathStep {
        endpoint,
        records,
        trials,
        rejected_trials,
        duration_s: conditions.duration_s,
        energy: FiniteDuctilePathEnergy {
            material,
            inertia,
            balance_error_j,
        },
    })
}

/// Read original retained Piola response through its native force assembler;
/// no return map, strain operator or body-load integrator is duplicated.
#[allow(clippy::too_many_arguments)]
fn initial_equilibrium<L>(
    solid: &FiniteDuctileTetSolid<'_, '_, L>,
    supports: &[(usize, f64)],
    forces: &[(usize, f64)],
    body: [f64; 3],
    conditions: DuctileStepConditions<'_>,
    settings: IncrementSettings,
    tolerance: f64,
    cx: &Cx<'_>,
) -> Result<(), FiniteDuctilePathError>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    let mut point_conditions = buffer(solid.materials.len())?;
    let mut frozen = buffer(solid.materials.len())?;
    for (element, (model, state)) in solid
        .materials
        .iter()
        .zip(&solid.state.element_states)
        .enumerate()
    {
        if element % 128 == 0 {
            checkpoint(cx)?;
        }
        let temperature_k = match conditions.temperatures_k {
            DuctileTemperatures::Uniform(value) => value,
            DuctileTemperatures::PerPoint(values) => values[element],
            DuctileTemperatures::PerHexElement(_) => unreachable!("coverage admitted above"),
        };
        point_conditions.push(DuctileJ2Conditions {
            duration_s: conditions.duration_s,
            temperature_k,
        });
        frozen.push(FrozenPoint {
            law: model.native(),
            accepted: state.native_state(),
        });
    }
    let problem = HyperTetProblem {
        nodes_m: solid.state.nodes_m,
        tetrahedra: solid.state.tetrahedra,
        material: HyperelasticLawField::PerElement(&frozen),
        prescribed_m: supports,
        nodal_forces_n: forces,
        body_force_n_m3: body,
        budget: solid.state.budget,
        settings: HyperTetSettings {
            load_steps: 1,
            max_newton_iterations: settings.max_newton_iterations,
            relative_tolerance: settings.relative_tolerance,
            max_backtracks: settings.max_backtracks,
            linear_tolerance: settings.linear_tolerance,
            max_linear_iterations: settings.max_linear_iterations,
        },
    };
    problem
        .validate_inputs()
        .map_err(FiniteDuctileTetError::from)?;
    let external = problem.external(&solid.state.elements, 1.0);
    let operator = Operator {
        elements: &solid.state.elements,
        materials: solid.materials,
        accepted: &solid.state.element_states,
        conditions: &point_conditions,
        inertia: None,
        dofs: solid.state.displacement_m.as_flattened().len(),
        refusal: RefCell::new(None),
    };
    let (mut residual, _) = operator
        .internal(
            &solid.state.element_states,
            solid.state.displacement_m.as_flattened(),
            None,
            cx,
        )
        .map_err(|error| operator.endpoint_refusal(error))?;
    for (dof, (residual, force)) in residual.iter_mut().zip(external).enumerate() {
        if dof % 256 == 0 {
            checkpoint(cx)?;
        }
        *residual -= force;
        if !residual.is_finite() {
            return Err(FiniteDuctilePathError::Invalid(
                "initial free-force residual cannot be represented",
            ));
        }
    }
    for (index, &(dof, _)) in supports.iter().enumerate() {
        if index % 256 == 0 {
            checkpoint(cx)?;
        }
        residual[dof] = 0.0;
    }
    let mut maximum = 0.0_f64;
    for (dof, value) in residual.into_iter().enumerate() {
        if dof % 256 == 0 {
            checkpoint(cx)?;
        }
        maximum = maximum.max(value.abs());
    }
    if maximum > tolerance {
        return Err(FiniteDuctilePathError::InitialEquilibrium {
            maximum_residual_n: maximum,
            tolerance_n: tolerance,
        });
    }
    checkpoint(cx)?;
    Ok(())
}

fn validate_controls(
    c: FiniteDuctilePathControls,
    kinetic: Option<f64>,
) -> Result<(), FiniteDuctilePathError> {
    if ![
        c.minimum_fraction,
        c.initial_fraction,
        c.maximum_fraction,
        c.maximum_damage_increment,
    ]
    .iter()
    .all(|v| v.is_finite() && *v > 0.0)
        || c.minimum_fraction > c.initial_fraction
        || c.initial_fraction > c.maximum_fraction
        || c.maximum_fraction > 1.0
        || c.maximum_damage_increment > 1.0
        || c.maximum_accepted_steps == 0
        || c.maximum_trials == 0
        || ![
            c.initial_equilibrium_tolerance_n,
            c.work.equilibrium_work_tolerance_j,
            c.work.maximum_endpoint_work_defect_j,
            c.work.yield_tolerance_pa,
        ]
        .iter()
        .all(|v| v.is_finite() && *v >= 0.0)
        || kinetic.is_some_and(|v| !v.is_finite() || v < 0.0)
    {
        return Err(FiniteDuctilePathError::Invalid(
            "ordered positive fractions, positive work counts and finite nonnegative tolerances required",
        ));
    }
    Ok(())
}
fn buffer<T>(capacity: usize) -> Result<Vec<T>, FiniteDuctilePathError> {
    let mut result = Vec::new();
    result
        .try_reserve_exact(capacity)
        .map_err(|_| FiniteDuctilePathError::Allocation)?;
    Ok(result)
}
fn interpolate(start: f64, end: f64, fraction: f64) -> f64 {
    if fraction == 1.0 {
        end
    } else {
        fraction.mul_add(end, (1.0 - fraction) * start)
    }
}
fn growth(allowance: f64, observed: f64) -> f64 {
    if observed > 0.0 {
        (0.8 * allowance / observed).clamp(0.5, 2.0)
    } else {
        2.0
    }
}
fn retryable(error: &FiniteDuctileTetError) -> bool {
    matches!(
        error,
        FiniteDuctileTetError::Solve(
            SmallStrainTetError::NewtonStalled { .. }
                | SmallStrainTetError::LinearSolveFailed { .. }
        )
    ) || matches!(error, FiniteDuctileTetError::Material {
        source: FiniteDuctileError::DamageIncrement { proposed, .. }, ..
    } if proposed.is_finite())
        || matches!(error, FiniteDuctileTetError::Energy {
        quantity: "sum of absolute element integration defects" | "backward-Euler numerical kinetic loss", value_j, ..
    } if value_j.is_finite())
}
fn refine(
    progress: f64,
    span: f64,
    minimum: f64,
    cause: FiniteDuctileRefinementCause,
) -> Result<f64, FiniteDuctilePathError> {
    let smaller = 0.5 * span;
    if smaller < minimum || progress + smaller <= progress {
        Err(FiniteDuctilePathError::MinimumInterval {
            accepted_fraction: progress,
            attempted_fraction: span,
            cause,
        })
    } else {
        Ok(smaller)
    }
}
