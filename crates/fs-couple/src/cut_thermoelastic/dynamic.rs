//! Actual original-mass inertia coupled to the same finite thermal timeline.
//!
//! The existing generalized-alpha owner at rho_inf=1 supplies average-
//! acceleration kinematics. Native old/new endpoint resistance is averaged,
//! including temperature-dependent stiffness and integrated free strain.
//! No static equilibrium solve replaces the accepted dynamic displacement.
//!
//! Reversible heat uses only physical bulk source exchange. Represented kinetic
//! change, native resistance work and external work have an independent gate;
//! a timestep/roundoff defect cannot be deposited as heat. Ghost work remains a
//! separately bounded numerical account. Damping and geometric G are not added
//! to this small-strain evolution. Optional G is recovered at the ACTUAL motion
//! endpoint for downstream frozen observations, with no evolving-G energy claim.

mod heat;

use fs_cutfem::elastic3::mass::MaterialMass3;
use fs_cutfem::elastic3::thermal::ThermoelasticState3;
use fs_solid::cutfront::dynamics3::{
    CutSecondOrderProblem3, CutThermalMidpointSettings3,
};
use fs_solver::FlexiblePreconditioner;
use fs_time::galpha::{
    ImplicitStepTelemetry, SecondOrderOperatorWeights, SecondOrderProblem, SecondOrderState,
};
use fs_time::galpha::initialization::{InitialSolveConfig, second_order_acceleration};

use super::*;

/// Native inertia controls, independent of source/heat/energy allowances.
/// The native step bounds retained history and Krylov workspace before cloning.
#[derive(Clone, Copy, Debug)]
pub struct CutThermoelasticDynamicsSettings3 {
    pub initial: InitialSolveConfig,
    pub midpoint: CutThermalMidpointSettings3,
}

struct DynamicHistory {
    temperatures_k: Vec<f64>,
    specific_h_j_kg: Vec<f64>,
    external_n: Vec<f64>,
    motion: SecondOrderState,
    endpoint: ThermoelasticState3,
    native_velocity_m_s: Vec<f64>,
    native_acceleration_m_s2: Vec<f64>,
    kinetic_energy_j: f64,
}

/// Accepted dynamic displacement, velocity and acceleration on the original
/// caloric/material specimen. Geometry and constraint maps are never copied.
pub struct DynamicCutThermoelasticSolid3<'m> {
    bindings: Bindings<'m>,
    history: DynamicHistory,
    scales: Vec<f64>,
    fixed: Vec<bool>,
    occupancy: Vec<f64>,
    settings: ThermoelasticEquilibriumSettings3,
    dynamics: CutThermoelasticDynamicsSettings3,
}

/// Closed inertial/caloric interval. Absolute kinetic values are diagnostics;
/// the balance uses its independently polarized retained increment.
#[derive(Debug)]
pub struct DynamicCutThermoelasticReport3 {
    pub coupled: CutThermoelasticReport3,
    pub native: ImplicitStepTelemetry,
    pub kinetic_energy_j: [f64; 2],
    pub kinetic_energy_change_j: f64,
    pub initial_force_residual_norm_n: f64,
    pub endpoint_force_residual_norm_n: f64,
    retained_mechanical_change: Dd,
}
impl DynamicCutThermoelasticReport3 {
    /// Kinetic plus physical elastic internal change; no caloric energy or
    /// numerical work is hidden in this value.
    pub const fn retained_mechanical_energy_change_j(&self) -> Dd {
        self.retained_mechanical_change
    }
}

struct Accepted<'t, O> {
    thermal: DepositionSimulation<'t>,
    observer: O,
    history: DynamicHistory,
    report: DynamicCutThermoelasticReport3,
}
enum Trial<'t, O> {
    Accepted(Accepted<'t, O>),
    Retry { heat_j: Vec<f64>, heat_defect_j: f64, combined_residual_j: f64 },
}

struct MassDiagonal<'a, 'b> {
    model: &'a CutSecondOrderProblem3<'b>,
    displacement: &'a [f64],
}
impl FlexiblePreconditioner for MassDiagonal<'_, '_> {
    fn apply(&self, iteration: usize, residual: &[f64], output: &mut [f64]) {
        self.model.preconditioner_apply(self.displacement,
            SecondOrderOperatorWeights { mass: 1.0, damping: 0.0, tangent: 0.0 },
            0, iteration, residual, output);
    }
}

fn dynamics_controls<E>(settings: CutThermoelasticDynamicsSettings3)
    -> Result<(), CutThermoelasticError3<E>>
{
    if settings.midpoint.maximum_steps == 0
        || settings.midpoint.maximum_workspace_components == 0
        || !settings.midpoint.maximum_force_residual_n.is_finite()
        || settings.midpoint.maximum_force_residual_n < 0.0
    {
        return Err(CutThermoelasticError3::Invalid("bounded native inertia history/workspace and finite force allowance required"));
    }
    Ok(())
}

impl<'m> CutThermoelasticSolid3<'m> {
    /// Launch physical inertia from the accepted loaded displacement and an
    /// explicit native velocity and external force. The new force may release
    /// the previous equilibrium load at this same initial instant: acceleration
    /// is solved consistently with that supplied force by the original mass
    /// owner, rather than retaining the equilibrium's zero acceleration.
    ///
    /// This consumes the static history but changes neither the native material
    /// nor thermal timeline. Fixed-node velocity must be exactly zero. Initial
    /// velocity is an authored physical initial condition, not a heat source.
    #[allow(clippy::too_many_arguments)]
    pub fn begin_transient(
        self, native: CutElasticityRef3<'_>, native_velocity_m_s: &[f64],
        initial_external_n: &[f64], dynamics: CutThermoelasticDynamicsSettings3, cx: &Cx<'_>,
    ) -> Result<DynamicCutThermoelasticSolid3<'m>, CutThermoelasticError3> {
        poll(cx)?; dynamics_controls(dynamics)?;
        forces(native, initial_external_n, self.settings.maximum_native_dofs)?;
        native::validate(native, &self.history.equilibrium)?;
        self.bindings.admit_native(native, cx)?;
        if native::scales(native) != self.scales || native::fixed(native) != self.fixed
            || native_velocity_m_s.len() != native::n(native)
            || native_velocity_m_s.iter().enumerate().any(|(i, v)|
                !v.is_finite() || (self.fixed[i / 3] && *v != 0.0))
        {
            return Err(CutThermoelasticError3::Invalid("unchanged native supports/scales and finite compatible initial velocity required"));
        }
        let equilibrium = self.history.equilibrium;
        let endpoint = native::recover_state(native, equilibrium.solution.coefficients(),
            equilibrium.eigenstrain, equilibrium.mass, self.settings, cx)?;
        let occupancy = vec![1.0; self.bindings.cells.len()];
        let (motion, acceleration) = {
            let model = CutSecondOrderProblem3::new_with_eigenstrain(
                native, &endpoint.mass, &occupancy, &endpoint.eigenstrain,
                self.settings.maximum_kinematic_strain, self.settings.maximum_points,
                self.settings.maximum_native_dofs, cx.cancel_gate(),
            )?;
            let q = model.restrict_force(&endpoint.displacement)?;
            let v = model.restrict_force(native_velocity_m_s)?;
            let load = model.restrict_force(initial_external_n)?;
            let preconditioner = MassDiagonal { model: &model, displacement: &q };
            let initial = second_order_acceleration(&model, &q, &v, &load,
                &preconditioner, dynamics.initial,
                dynamics.midpoint.maximum_workspace_components,
                &mut || cx.is_cancel_requested())?;
            let mut inertia = vec![0.0; q.len()];
            model.mass_apply(&initial.value, &mut inertia);
            poll(cx)?;
            let resistance = model.restrict_force(&endpoint.resistance)?;
            let mut residual = 0.0_f64;
            for i in 0..q.len() {
                if i % 192 == 0 { poll(cx)?; }
                let value = finite((Dd::from_f64(inertia[i])
                    + Dd::from_f64(resistance[i]) - Dd::from_f64(load[i])).to_f64())?;
                residual = residual.hypot(value);
            }
            if finite(residual)? > dynamics.midpoint.maximum_force_residual_n {
                return Err(CutThermoelasticError3::Invalid("consistent initial inertia failed the absolute force allowance"));
            }
            let acceleration = model.expand_field(&initial.value)?;
            (SecondOrderState::new(self.history.time_s, &q, &v, &initial.value), acceleration)
        };
        let kinetic_energy_j = kinetic_energy(&endpoint.mass, native_velocity_m_s, cx)?;
        let history = DynamicHistory {
            temperatures_k: self.history.temperatures_k,
            specific_h_j_kg: self.history.specific_h_j_kg,
            external_n: initial_external_n.to_vec(), motion, endpoint,
            native_velocity_m_s: native_velocity_m_s.to_vec(),
            native_acceleration_m_s2: acceleration, kinetic_energy_j,
        };
        poll(cx)?;
        Ok(DynamicCutThermoelasticSolid3 {
            bindings: self.bindings, history, scales: self.scales, fixed: self.fixed,
            occupancy, settings: self.settings, dynamics,
        })
    }
}

impl DynamicCutThermoelasticSolid3<'_> {
    pub fn time_s(&self) -> f64 { self.history.motion.t }
    pub fn step(&self) -> usize { self.history.motion.steps }
    /// Actual independent free-coordinate checkpoint of the original integrator.
    pub fn motion(&self) -> &SecondOrderState { &self.history.motion }
    /// Physical source response/preparations at actual dynamic displacement.
    /// Its resistance is instantaneous elastic resistance; it is NOT an
    /// inertia-balanced support reaction or a solved static equilibrium.
    pub fn state(&self) -> &ThermoelasticState3 { &self.history.endpoint }
    pub fn native_velocity_m_s(&self) -> &[f64] { &self.history.native_velocity_m_s }
    pub fn native_acceleration_m_s2(&self) -> &[f64] { &self.history.native_acceleration_m_s2 }
    pub fn temperatures_k(&self) -> &[f64] { &self.history.temperatures_k }
    pub fn nodal_forces_n(&self) -> &[f64] { &self.history.external_n }
    pub fn kinetic_energy_j(&self) -> f64 { self.history.kinetic_energy_j }

    /// Advance real inertia and finite heating as one accepted state.
    /// External forces use endpoint-average work; changing endpoint forces
    /// describe that time quadrature, not an unreported impulse. Every thermal
    /// Picard trial restarts the same accepted q/v/a and observer history.
    ///
    /// The native step conservatively admits its old, stage and new mechanical
    /// states against the trial source; a cross-temperature domain exit refuses.
    /// The independent source exchange additionally admits both physical and
    /// constitutive endpoint strains at their actual endpoint temperatures.
    #[allow(clippy::too_many_arguments)]
    pub fn advance<'t, O: DepositionObserver>(
        &mut self, mut native: CutElasticityMut3<'_>,
        thermal: &mut DepositionSimulation<'t>, observer: &mut O,
        target_s: f64, external_n: &[f64], thermal_config: StepConfig,
        config: CutThermoelasticConfig3, cx: &Cx<'_>,
    ) -> Result<DynamicCutThermoelasticReport3, CutThermoelasticError3<O::Error>> {
        let accepted = self.attempt(&mut native, thermal, observer, target_s,
            external_n, thermal_config, config, |_| 0.0, cx)?;
        // The material callback has passed its final cancellation poll. Retain
        // displaced values until every accepted owner moves, including an O
        // with an arbitrary destructor. No new callback/validation follows.
        let old_history = core::mem::replace(&mut self.history, accepted.history);
        let old_thermal = core::mem::replace(thermal, accepted.thermal);
        let old_observer = core::mem::replace(observer, accepted.observer);
        drop((old_history, old_thermal, old_observer));
        Ok(accepted.report)
    }

    #[allow(clippy::too_many_arguments)]
    fn attempt<'t, O: DepositionObserver>(
        &self, native: &mut CutElasticityMut3<'_>, thermal: &DepositionSimulation<'t>,
        observer: &O, target_s: f64, external_n: &[f64],
        thermal_config: StepConfig, config: CutThermoelasticConfig3,
        radiation_net: impl Fn(&O) -> f64, cx: &Cx<'_>,
    ) -> Result<Accepted<'t, O>, CutThermoelasticError3<O::Error>> {
        poll(cx)?; controls(config)?; dynamics_controls(self.dynamics)?;
        let view = native.as_ref();
        forces(view, external_n, self.settings.maximum_native_dofs)?;
        if self.bindings.cells.len() > config.maximum_cells
            || native::scales(view) != self.scales || native::fixed(view) != self.fixed
            || self.history.motion.steps >= self.dynamics.midpoint.maximum_steps
        {
            return Err(CutThermoelasticError3::Invalid("fixed bounded original cells/supports/scales and remaining inertial history budget required"));
        }
        native::validate_state(view, &self.history.endpoint)?;
        self.bindings.admit_native(view, cx)?;
        let (old_t, old_h) = self.bindings.admit(thermal, cx)?;
        if thermal.time_s().to_bits() != self.time_s().to_bits()
            || old_t.iter().zip(&self.history.temperatures_k).any(|(a,b)| a.to_bits()!=b.to_bits())
            || old_h.iter().zip(&self.history.specific_h_j_kg).any(|(a,b)| a.to_bits()!=b.to_bits())
        {
            return Err(CutThermoelasticError3::Invalid("accepted inertial and caloric histories must remain exactly paired"));
        }
        let dt = target_s - self.time_s();
        if !target_s.is_finite() || !dt.is_finite() || dt <= 0.0
            || !thermal_config.energy_tolerance_j.is_finite()
            || thermal_config.energy_tolerance_j < 0.0
            || thermal_config.energy_tolerance_j > config.energy_tolerance_j
        {
            return Err(CutThermoelasticError3::Invalid("positive finite inertia interval and compatible caloric energy tolerance required"));
        }
        let mut applied = vec![0.0; self.bindings.cells.len()];
        let mut last = [f64::INFINITY; 2];
        for iteration in 1..=config.maximum_iterations {
            poll(cx)?;
            let coupled = ReversibleHeat::new(&self.bindings, observer.clone(), &applied, dt, cx)?;
            let (thermal_next, coupled_next, audit) = thermal.advance_observed(
                &coupled, target_s, thermal_config,
                || cx.checkpoint().map_err(|_| DepositionError::Cancelled),
            )?;
            drop(coupled);
            let observer_next = coupled_next.inner;
            let applied_actual = coupled_next.applied_j;
            let radiation = finite(radiation_net(&observer_next))?;
            let (temperatures, enthalpies) = self.bindings.admit(&thermal_next, cx)?;
            let materials = self.bindings.materials(&temperatures, cx)?;
            let field = |bounds| self.bindings.material_at(bounds, &materials);
            let decision = native.refresh(ElasticMaterialField3::ThermoelasticByCell(&field),
                self.settings.refresh, |candidate| -> Result<_, CutThermoelasticError3<O::Error>> {
                    let (eigenstrain, mass) = native::thermal_parts(candidate, self.settings, cx)?;
                    let (step, displacement, velocity, acceleration) = {
                        let model = CutSecondOrderProblem3::new_with_eigenstrain(
                            candidate, &mass, &self.occupancy, &eigenstrain,
                            self.settings.maximum_kinematic_strain, self.settings.maximum_points,
                            self.settings.maximum_native_dofs, cx.cancel_gate(),
                        )?;
                        let step = model.thermal_midpoint_candidate(
                            &self.history.motion, &self.history.endpoint.resistance,
                            &self.history.external_n, external_n, dt, self.dynamics.midpoint,
                        )?;
                        let displacement = model.expand_field(&step.state.q)?;
                        let velocity = model.expand_field(&step.state.v)?;
                        let acceleration = model.expand_field(&step.state.a)?;
                        (step, displacement, velocity, acceleration)
                    };
                    if step.state.t.to_bits() != target_s.to_bits() {
                        return Err(CutThermoelasticError3::Invalid("native inertia and caloric endpoint clocks differ"));
                    }
                    let endpoint = native::recover_state(candidate, &displacement,
                        eigenstrain, mass, self.settings, cx)?;
                    let exchanges = native::exchange(candidate, &self.bindings.laws,
                        &self.history.temperatures_k, &self.history.endpoint.displacement,
                        &endpoint.displacement, self.settings.maximum_kinematic_strain,
                        config.exchange, cx)?;
                    let kinetic_change = kinetic_increment(&endpoint.mass,
                        &self.history.native_velocity_m_s, &velocity, cx)?;
                    let kinetic_after = kinetic_energy(&endpoint.mass, &velocity, cx)?;
                    let external = work(&self.history.external_n, external_n,
                        &self.history.endpoint.displacement, &endpoint.displacement, cx)?;
                    let native_work = work(&self.history.endpoint.resistance, &endpoint.resistance,
                        &self.history.endpoint.displacement, &endpoint.displacement, cx)?;
                    let (coupled, proposed, internal_change) = account_exchange(
                        &self.bindings, &self.history.temperatures_k, &temperatures,
                        exchanges, &applied_actual, external, native_work, kinetic_change,
                        audit, radiation, iteration, config, cx)?;
                    let converged = coupled.heat_defect_j <= config.heat_tolerance_j
                        && coupled.thermal_residual_j.abs() <= config.energy_tolerance_j
                        && coupled.combined_residual_j.abs() <= config.energy_tolerance_j;
                    if converged {
                        if coupled.numerical_work_j.abs() > config.maximum_numerical_work_j {
                            return Err(CutThermoelasticError3::NumericalWork {
                                work_j: coupled.numerical_work_j,
                                allowance_j: config.maximum_numerical_work_j,
                            });
                        }
                        let report = DynamicCutThermoelasticReport3 {
                            coupled, native: step.telemetry,
                            kinetic_energy_j: [self.history.kinetic_energy_j, kinetic_after],
                            kinetic_energy_change_j: finite(kinetic_change.to_f64())?,
                            initial_force_residual_norm_n: step.initial_residual_norm_n,
                            endpoint_force_residual_norm_n: step.endpoint_residual_norm_n,
                            retained_mechanical_change: kinetic_change + internal_change,
                        };
                        let history = DynamicHistory {
                            temperatures_k: temperatures, specific_h_j_kg: enthalpies,
                            external_n: external_n.to_vec(), motion: step.state, endpoint,
                            native_velocity_m_s: velocity, native_acceleration_m_s2: acceleration,
                            kinetic_energy_j: kinetic_after,
                        };
                        Ok(MaterialRefreshDecision3::Accept(Trial::Accepted(Accepted {
                            thermal: thermal_next, observer: observer_next, history, report,
                        })))
                    } else {
                        drop((step, endpoint, observer_next, thermal_next));
                        Ok(MaterialRefreshDecision3::Reject(Trial::Retry {
                            heat_j: proposed, heat_defect_j: coupled.heat_defect_j,
                            combined_residual_j: coupled.combined_residual_j,
                        }))
                    }
                }, || checkpoint(cx))?;
            match decision {
                MaterialRefreshDecision3::Accept(Trial::Accepted(accepted)) => return Ok(accepted),
                MaterialRefreshDecision3::Reject(Trial::Retry {
                    heat_j, heat_defect_j, combined_residual_j,
                }) => {
                    last = [heat_defect_j, combined_residual_j];
                    for (i, (used, target)) in applied.iter_mut().zip(heat_j).enumerate() {
                        if i % 64 == 0 { poll(cx)?; }
                        *used = finite(config.heat_relaxation.mul_add(target - *used, *used))?;
                    }
                }
                _ => unreachable!("native material decision preserves its supplied payload"),
            }
        }
        Err(CutThermoelasticError3::NotConverged {
            iterations: config.maximum_iterations, heat_defect_j: last[0],
            combined_residual_j: last[1],
        })
    }
}

// Pair represented endpoints before summing; no rounded absolute-work or
// kinetic-energy subtraction. If time-update rounding loses a physical
// increment, the independent mechanical balance must refuse it, never heat it.
fn work<E>(f0: &[f64], f1: &[f64], q0: &[f64], q1: &[f64], cx: &Cx<'_>)
    -> Result<Dd, CutThermoelasticError3<E>>
{
    if f0.len()!=f1.len() || f0.len()!=q0.len() || q0.len()!=q1.len() {
        return Err(CutThermoelasticError3::Invalid("paired original force/displacement dimensions required"));
    }
    let mut result = Dd::ZERO;
    for i in 0..q0.len() {
        if i % 192 == 0 { poll(cx)?; }
        result = result + Dd::from_f64(0.5) * (Dd::from_f64(f0[i]) + Dd::from_f64(f1[i]))
            * (Dd::from_f64(q1[i]) - Dd::from_f64(q0[i]));
    }
    finite(result.to_f64())?;
    Ok(result)
}

fn kinetic_increment<E>(mass: &MaterialMass3, before: &[f64], after: &[f64], cx: &Cx<'_>)
    -> Result<Dd, CutThermoelasticError3<E>>
{
    if before.len()!=mass.n() || after.len()!=mass.n() {
        return Err(CutThermoelasticError3::Invalid("kinetic field must use original native mass coordinates"));
    }
    let mut result = Dd::ZERO;
    for (index, (_, [a,b], value)) in mass.reference_scalar_entries().enumerate() {
        if index % 64 == 0 { poll(cx)?; }
        for component in 0..3 {
            let i = 3*a+component;
            let j = 3*b+component;
            let average_i = Dd::from_f64(0.5) * (Dd::from_f64(after[i]) + Dd::from_f64(before[i]));
            let change_j = Dd::from_f64(after[j]) - Dd::from_f64(before[j]);
            result = result + Dd::from_f64(value) * average_i * change_j;
            if a!=b {
                let average_j = Dd::from_f64(0.5) * (Dd::from_f64(after[j]) + Dd::from_f64(before[j]));
                let change_i = Dd::from_f64(after[i]) - Dd::from_f64(before[i]);
                result = result + Dd::from_f64(value) * average_j * change_i;
            }
        }
    }
    finite(result.to_f64())?;
    Ok(result)
}

fn kinetic_energy<E>(mass: &MaterialMass3, velocity: &[f64], cx: &Cx<'_>)
    -> Result<f64, CutThermoelasticError3<E>>
{
    if velocity.len()!=mass.n() {
        return Err(CutThermoelasticError3::Invalid("kinetic field must use original native mass coordinates"));
    }
    let mut result = Dd::ZERO;
    for (index, (_, [a,b], value)) in mass.reference_scalar_entries().enumerate() {
        if index % 64 == 0 { poll(cx)?; }
        let factor = Dd::from_f64(if a==b {0.5}else{1.0}) * Dd::from_f64(value);
        for component in 0..3 {
            result = result + factor * Dd::from_f64(velocity[3*a+component])
                * Dd::from_f64(velocity[3*b+component]);
        }
    }
    let value = finite(result.to_f64())?;
    if value < 0.0 {
        return Err(CutThermoelasticError3::Invalid("original consistent mass produced negative kinetic energy"));
    }
    Ok(value)
}
