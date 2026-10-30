//! Native mechanics of a finite-thickness, quasi-static shear layer.
//!
//! This owner declares a first-order director interpolation through the actual
//! layer thickness. If `d = x_b - x_a`, its mean shear rate is the tangential
//! part of `v_b - v_a - (omega_a + omega_b) cross d / 2`. It follows by
//! integrating the shear strain with linearly interpolated face directors.
//! The existing retained elastic-viscoplastic law supplies its mean shear
//! effort. Virtual work then gives force `+F` on face A, `-F` on face B, and
//! face couples `d cross F / 2` on BOTH faces. These couples are part of this
//! declared interpolation; they are not arbitrary repairs to a point contact.
//! Together the actual face wrenches conserve total orbital plus spin momentum
//! and do no work under a common rigid translation or rotation.
//!
//! This is only the shear branch of a massless mechanical layer. It supplies
//! no normal, bending, torsional or layer-inertia law. The geometry is frozen
//! during the momentum phase; poses are not advanced. Finite rotations of the
//! retained material frame, geometry changes, opening, wear and transport need
//! their own state transfer. The two real face points are never moved together.

use super::{
    ThirdBodyHeatInput, ThirdBodyHeatReport, ThirdBodyHeatRoute, ThirdBodyHeatStep,
    ThirdBodyHeatTolerances, ThirdBodyThermoelasticInput, half_distances,
};
use fs_mbd::{DynamicsError, MassProperties, RigidBodyState, Vec3};
use fs_material::state_point::ShearLayerThermalSlopes;
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress,
};
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};
use fs_tribo::layer_contact::{ThirdBodyLayerContact, ThirdBodyLayerState};
use fs_tribo::partial_slip::TangentFrame;
use fs_tribo::{
    TangentialLaw, TangentialSolveConfig, TangentialStep, TangentialStepInput, TriboError,
};
use std::cell::{Cell, RefCell};

/// Source-pinned thermoelastic layer feedback on the original rigid/heat solve.
pub mod thermal;

/// Full source tensor shear with original rigid momenta and finite caloric cells.
pub mod tensor;

/// The actual opposing face centres and unchanged material tangent frame.
/// The normal points from face A to face B. Arms are COM-relative in each
/// original body's principal frame. Area and thickness are physical dimensions.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidLayerGeometry {
    pub arms_body_m: [Vec3; 2],
    pub frame: TangentFrame,
    pub thickness_m: f64,
    pub area_m2: f64,
}

/// Original retained layer deformation, bound to its geometry and clock.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidLayerHistory {
    time_s: f64,
    geometry: RigidLayerGeometry,
    law: ThirdBodyLayerContact,
    deformation: ThirdBodyLayerState,
}
impl RigidLayerHistory {
    pub fn new(
        law: &ThirdBodyLayerContact,
        geometry: RigidLayerGeometry,
        time_s: f64,
    ) -> Result<Self, RigidLayerError> {
        geometry.validate()?;
        if !time_s.is_finite() || time_s < 0.0 {
            return Err(RigidLayerError::Invalid(
                "finite nonnegative layer clock required",
            ));
        }
        Ok(Self {
            time_s,
            geometry,
            law: law.clone(),
            deformation: law.initial_state(),
        })
    }
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub const fn geometry(&self) -> RigidLayerGeometry {
        self.geometry
    }
    pub const fn law(&self) -> &ThirdBodyLayerContact {
        &self.law
    }
    pub const fn deformation(&self) -> &ThirdBodyLayerState {
        &self.deformation
    }
}

/// Independent nonlinear, geometry, force and native conservation gates.
#[derive(Clone, Copy, Debug)]
pub struct RigidLayerSolve {
    pub maximum_trials: usize,
    pub maximum_newton_iterations: usize,
    pub difference_step_n: f64,
    pub geometry_tolerance_m: f64,
    pub force_tolerance_n: f64,
    pub linear_momentum_tolerance_ns: f64,
    pub angular_momentum_tolerance_nms: f64,
    pub energy_tolerance_j: f64,
    pub newton: NewtonKrylovConfig,
}

#[derive(Debug)]
pub enum RigidLayerError {
    Invalid(&'static str),
    Dynamics(DynamicsError),
    Law(TriboError),
    Newton(NewtonError),
    NotConverged(NewtonReport),
    Interrupted(SolverRunProgress),
    TrialBudget,
    Cancelled,
}
impl From<DynamicsError> for RigidLayerError {
    fn from(error: DynamicsError) -> Self {
        Self::Dynamics(error)
    }
}
impl From<TriboError> for RigidLayerError {
    fn from(error: TriboError) -> Self {
        Self::Law(error)
    }
}
impl From<NewtonError> for RigidLayerError {
    fn from(error: NewtonError) -> Self {
        Self::Newton(error)
    }
}
impl core::fmt::Display for RigidLayerError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native finite-layer shear: {self:?}")
    }
}
impl std::error::Error for RigidLayerError {}

/// Native work and momentum accounts on the two separated material faces.
#[derive(Debug)]
pub struct RigidLayerReport {
    pub constitutive: TangentialStep<ThirdBodyLayerState>,
    pub applied_force_on_a_n: Vec3,
    pub face_couple_nm: [Vec3; 2],
    /// Actual mean shear velocity, including the interpolated directors.
    pub midpoint_shear_velocity_mps: [f64; 2],
    /// Original pair mobility in the material tangent basis, kg^-1.
    pub inverse_shear_mass: [[f64; 2]; 2],
    pub force_residual_n: [f64; 2],
    pub kinetic_energy_change_j: f64,
    pub mechanical_work_j: f64,
    pub work_energy_residual_j: f64,
    pub total_energy_residual_j: f64,
    pub linear_momentum_change_ns: Vec3,
    pub angular_momentum_change_nms: Vec3,
}

/// Both original native bodies and retained layer history accept together.
pub struct RigidLayerCandidate {
    bodies: [RigidBodyState; 2],
    history: RigidLayerHistory,
    report: RigidLayerReport,
    nonlinear: NewtonReport,
    trials: usize,
}
impl RigidLayerCandidate {
    pub const fn bodies(&self) -> &[RigidBodyState; 2] {
        &self.bodies
    }
    pub const fn history(&self) -> &RigidLayerHistory {
        &self.history
    }
    pub const fn report(&self) -> &RigidLayerReport {
        &self.report
    }
    pub const fn nonlinear(&self) -> &NewtonReport {
        &self.nonlinear
    }
    pub const fn trials(&self) -> usize {
        self.trials
    }
    pub fn into_parts(self) -> ([RigidBodyState; 2], RigidLayerHistory, RigidLayerReport) {
        (self.bodies, self.history, self.report)
    }
}

impl RigidLayerGeometry {
    fn validate(self) -> Result<(), RigidLayerError> {
        if self.arms_body_m.iter().any(|v| !v.is_finite())
            || [self.thickness_m, self.area_m2]
                .iter()
                .any(|x| !x.is_finite() || *x <= 0.0)
        {
            return Err(RigidLayerError::Invalid(
                "finite arms, positive layer thickness and area required",
            ));
        }
        Ok(())
    }
}
impl RigidLayerSolve {
    fn validate(self) -> Result<(), RigidLayerError> {
        self.newton.validate()?;
        if self.maximum_trials == 0
            || self.maximum_newton_iterations == 0
            || [
                self.difference_step_n,
                self.geometry_tolerance_m,
                self.force_tolerance_n,
                self.linear_momentum_tolerance_ns,
                self.angular_momentum_tolerance_nms,
                self.energy_tolerance_j,
            ]
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.0)
        {
            return Err(RigidLayerError::Invalid(
                "positive bounded native layer solve required",
            ));
        }
        Ok(())
    }
}

/// Resolve held layer effort from the original two-body momentum response.
///
/// The original layer law uses midpoint elastic strain and physical Bingham
/// plastic work. Here its two shear coordinates are jointly solved against the
/// full original inertial response, including off-diagonal rotation coupling.
/// No isotropic effective mass replaces that response. The returned work closes
/// `Delta(K + U_layer) + D_plastic = 0`; no numerical deficit becomes heat.
///
/// The supplied normal load stays positive throughout this frozen phase by the
/// surrounding normal owner's admission. It is not solved here. The law's
/// fixed temperature is explicit; thermal evolution uses the heated consumer.
#[allow(clippy::too_many_arguments)]
pub fn step_rigid_layer_candidate(
    bodies: &[RigidBodyState; 2],
    properties: [MassProperties; 2],
    history: &RigidLayerHistory,
    law: &ThirdBodyLayerContact,
    target_s: f64,
    normal_n: f64,
    law_solve: TangentialSolveConfig,
    controls: RigidLayerSolve,
    force_guess_n: [f64; 2],
    cancelled: impl FnMut() -> bool,
) -> Result<RigidLayerCandidate, RigidLayerError> {
    let cancel = RefCell::new(cancelled);
    let poll = || {
        if (*cancel.borrow_mut())() {
            Err(RigidLayerError::Cancelled)
        } else {
            Ok(())
        }
    };
    poll()?;
    controls.validate()?;
    let dt = duration(history, target_s)?;
    if force_guess_n.iter().any(|x| !x.is_finite()) {
        return Err(RigidLayerError::Invalid(
            "finite shear force guess required",
        ));
    }
    // This consumer supplies no coefficient-work port. A changing potential
    // refuses even if a loose energy tolerance could conceal its work.
    let continued = law.continue_state(&history.deformation)?;
    let prepared = Prepared::new(
        bodies,
        properties,
        history.geometry,
        controls.geometry_tolerance_m,
    )?;
    let trials = Cell::new(0usize);
    let evaluate = |x: &[f64]| {
        poll()?;
        spend(&trials, controls.maximum_trials)?;
        if x.len() != 2 {
            return Err(RigidLayerError::Invalid(
                "two shear force coordinates required",
            ));
        }
        prepared.evaluate(law, &continued, dt, normal_n, law_solve, [x[0], x[1]])
    };
    let problem = FallibleFiniteDifference::new(2, controls.difference_step_n, |x: &[f64]| {
        evaluate(x).map(|stage| stage.report.force_residual_n.to_vec())
    })
    .map_err(numerical)?;
    let creation = NewtonKrylovState::new(&problem, force_guess_n.to_vec(), controls.newton);
    if let Some(error) = problem.take_failure() {
        return Err(numerical(error));
    }
    let mut solver = creation?;
    let result = solver.run_controlled(&problem, controls.maximum_newton_iterations, &mut || {
        (*cancel.borrow_mut())() || problem.has_failed()
    });
    if let Some(error) = problem.take_failure() {
        return Err(numerical(error));
    }
    poll()?;
    if result.progress != SolverRunProgress::Complete {
        return Err(RigidLayerError::Interrupted(result.progress));
    }
    if !result.report.converged {
        return Err(RigidLayerError::NotConverged(result.report));
    }
    let stage = evaluate(&solver.x)?;
    stage.report.admit(controls)?;
    poll()?;
    Ok(RigidLayerCandidate {
        bodies: stage.bodies,
        history: RigidLayerHistory {
            time_s: target_s,
            geometry: history.geometry,
            law: law.clone(),
            deformation: stage.report.constitutive.state.clone(),
        },
        report: stage.report,
        nonlinear: result.report,
        trials: trials.get(),
    })
}

/// Original temperature-resolved layer law and thermal transport. Only its
/// yield force and plastic resistance may vary: this owner has no elastic
/// coefficient-work supply. The callback retains material/source authority.
pub struct RigidLayerHeatSource<P> {
    pub law: ThirdBodyLayerContact,
    /// Original body A, body B, layer conductivities [W/(m K)].
    pub conductivity_w_m_k: [f64; 3],
    /// Original two area-specific interface resistances [m² K/W].
    pub area_resistance_m2_k_w: [f64; 2],
    /// OTHER sources and links; physical layer plastic work is supplied once
    /// by this owner. No authored split of that work is accepted.
    pub transport: StepCoefficients<P>,
}

#[derive(Clone, Copy, Debug)]
pub struct HeatedRigidLayerSolve {
    pub mechanical: RigidLayerSolve,
    pub maximum_cells: usize,
    /// Numerical energy-coordinate scale [J/N], never a thermal capacity.
    pub energy_coordinate_j_per_n: f64,
    pub heat: ThirdBodyHeatTolerances,
    pub combined_energy_tolerance_j: f64,
    /// Audit original caloric entropy integrals at the final trial only.
    /// Every active thermal law must supply its original entropy integral.
    pub audit_entropy: bool,
}

#[derive(Debug)]
pub enum HeatedRigidLayerError<E> {
    Mechanical(RigidLayerError),
    Thermal(DepositionError),
    Source(E),
}
impl<E> From<RigidLayerError> for HeatedRigidLayerError<E> {
    fn from(error: RigidLayerError) -> Self {
        Self::Mechanical(error)
    }
}
impl<E> From<DepositionError> for HeatedRigidLayerError<E> {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl<E> From<NewtonError> for HeatedRigidLayerError<E> {
    fn from(error: NewtonError) -> Self {
        Self::Mechanical(RigidLayerError::Newton(error))
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for HeatedRigidLayerError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "heated native finite-layer shear: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for HeatedRigidLayerError<E> {}

pub struct HeatedRigidLayerReport<P> {
    pub mechanical: RigidLayerReport,
    pub heat: ThirdBodyHeatReport,
    pub nonlinear: NewtonReport,
    pub trials: usize,
    pub source: P,
    /// Present only for the explicitly selected source-derived affine potential.
    pub thermoelastic: Option<RigidLayerThermoelasticReport>,
    /// Delta(K_bodies + U_layer_shear + H_original_network) - Q_external.
    pub combined_energy_residual_j: f64,
}

/// Actual two-temperature mechanical storage; none of these terms is plastic
/// or Maxwell dissipation. The native constitutive report retains its original
/// frozen-coefficient mechanical work account alongside this correction.
#[derive(Clone, Copy, Debug)]
pub struct RigidLayerThermoelasticReport {
    pub coefficient_free_energy_change_j: f64,
    pub free_energy_change_j: f64,
    pub mechanical_entropy_change_j_k: f64,
    pub entropy_energy_change_j: f64,
    pub internal_energy_change_j: f64,
    pub reversible_heat_j: f64,
}

/// Original body, retained deformation and finite enthalpy endpoints. Publish
/// these parts together; all input owners remain unchanged on every refusal.
pub struct HeatedRigidLayerCandidate<'t, P> {
    bodies: [RigidBodyState; 2],
    history: RigidLayerHistory,
    thermal: DepositionSimulation<'t>,
    report: HeatedRigidLayerReport<P>,
}
impl<'t, P> HeatedRigidLayerCandidate<'t, P> {
    pub const fn bodies(&self) -> &[RigidBodyState; 2] {
        &self.bodies
    }
    pub const fn history(&self) -> &RigidLayerHistory {
        &self.history
    }
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub const fn report(&self) -> &HeatedRigidLayerReport<P> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        [RigidBodyState; 2],
        RigidLayerHistory,
        DepositionSimulation<'t>,
        HeatedRigidLayerReport<P>,
    ) {
        (self.bodies, self.history, self.thermal, self.report)
    }
}
struct HeatStage<P> {
    mechanical: Stage,
    law: ThirdBodyLayerContact,
    source: P,
    thermoelastic: Option<RigidLayerThermoelasticReport>,
}

/// Joint finite-layer force, retained slip and ORIGINAL finite-cell heat solve.
///
/// Two shear forces and the original cell energy increments form ONE nonlinear
/// system. Every trial starts from accepted body momenta and strain; there is
/// no nested mechanical solve or sequential heat deposit. Dissipated plastic
/// work enters the existing layer cell, whose real caloric law/mass retains it
/// and whose original material/geometry conductances exchange heat with both
/// bodies. Changing dissipative coefficients sees the same trial layer Kelvin.
///
/// Native face points must equal the route's two physical face centres; area
/// and thickness must also match. The layer's mechanical inertia is neglected
/// by the explicitly quasi-static shear interpolation, even though its actual
/// caloric storage is retained. This is not a complete dynamic solid layer.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn step_heated_rigid_layer_candidate<'t, P, E>(
    bodies: &[RigidBodyState; 2],
    properties: [MassProperties; 2],
    history: &RigidLayerHistory,
    thermal: &DepositionSimulation<'t>,
    route: &ThirdBodyHeatRoute,
    target_s: f64,
    normal_n: f64,
    law_solve: TangentialSolveConfig,
    controls: HeatedRigidLayerSolve,
    force_guess_n: [f64; 2],
    energy_guess_j: &[f64],
    resolve: impl Fn(&DepositionSimulation<'t>, f64) -> Result<RigidLayerHeatSource<P>, E>,
    cancelled: impl FnMut() -> bool,
) -> Result<HeatedRigidLayerCandidate<'t, P>, HeatedRigidLayerError<E>> {
    step_heated_impl(bodies, properties, history, thermal, route, target_s,
        normal_n, law_solve, controls, force_guess_n, energy_guess_j, None, resolve, cancelled)
}

// The optional derivative is private to the source-bound thermal ingress. The
// original public callback cannot supply an unadmitted elastic temperature law.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn step_heated_impl<'t, P, E>(
    bodies: &[RigidBodyState; 2], properties: [MassProperties; 2],
    history: &RigidLayerHistory, thermal: &DepositionSimulation<'t>,
    route: &ThirdBodyHeatRoute, target_s: f64, normal_n: f64,
    law_solve: TangentialSolveConfig, controls: HeatedRigidLayerSolve,
    force_guess_n: [f64; 2], energy_guess_j: &[f64], thermal_slopes: Option<&ShearLayerThermalSlopes>,
    resolve: impl Fn(&DepositionSimulation<'t>, f64) -> Result<RigidLayerHeatSource<P>, E>,
    cancelled: impl FnMut() -> bool,
) -> Result<HeatedRigidLayerCandidate<'t, P>, HeatedRigidLayerError<E>> {
    let cancel = RefCell::new(cancelled);
    let poll = || {
        if (*cancel.borrow_mut())() {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    };
    poll()?;
    controls.mechanical.validate()?;
    let dt = duration(history, target_s)?;
    if history.time_s != thermal.time_s()
        || history.law.temperature_kelvin().to_bits()
            != thermal.temperature_k(route.layer_cell())?.to_bits()
        || controls.maximum_cells == 0
        || thermal.cells().len() > controls.maximum_cells
        || energy_guess_j.len() != thermal.cells().len()
        || energy_guess_j
            .iter()
            .chain(force_guess_n.iter())
            .any(|x| !x.is_finite())
        || [
            controls.energy_coordinate_j_per_n,
            controls.combined_energy_tolerance_j,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
    {
        return Err(RigidLayerError::Invalid("original synchronized layer temperature, clocks and bounded thermal coordinates required").into());
    }
    let prepared = Prepared::new(
        bodies,
        properties,
        history.geometry,
        controls.mechanical.geometry_tolerance_m,
    )?;
    prepared.admit_thermal_route(route, controls.mechanical.geometry_tolerance_m)?;
    let heat_step = ThirdBodyHeatStep::new(thermal, target_s, route, &poll)?;
    let old_thermodynamics = thermal_slopes.map(|slopes|
        history.law.thermodynamics_with_relaxation(&history.deformation,
            slopes.equilibrium_n_m_k, &slopes.relaxation_n_m_k)
            .map_err(RigidLayerError::from)).transpose()?;
    let count = heat_step.cell_ids().len();
    let scale = controls.energy_coordinate_j_per_n;
    let trials = Cell::new(0usize);
    let evaluate = |x: &[f64], entropy: bool| {
        poll()?;
        spend(&trials, controls.mechanical.maximum_trials)?;
        if x.len() != count + 2 {
            return Err(HeatedRigidLayerError::Mechanical(RigidLayerError::Invalid(
                "incomplete joint force/enthalpy iterate",
            )));
        }
        let energy: Vec<_> = x[2..].iter().map(|value| scale * value).collect();
        let coefficients =
            |endpoint: &DepositionSimulation<'t>,
             temperature: f64|
             -> Result<ThirdBodyThermoelasticInput<HeatStage<P>>, HeatedRigidLayerError<E>> {
                let source =
                    resolve(endpoint, temperature).map_err(HeatedRigidLayerError::Source)?;
                if source.law.temperature_kelvin().to_bits() != temperature.to_bits() {
                    return Err(RigidLayerError::Invalid(
                        "layer source must resolve at its actual trial Kelvin",
                    )
                    .into());
                }
                let continuation = if thermal_slopes.is_some() {
                    Some(source.law.continue_state_with_work(&history.deformation)
                        .map_err(RigidLayerError::from)?)
                } else { None };
                let continued = match &continuation {
                    Some(value) => value.state.clone(),
                    None => source.law.continue_state(&history.deformation)
                        .map_err(RigidLayerError::from)?,
                };
                let mechanical = prepared.evaluate(
                    &source.law,
                    &continued,
                    dt,
                    normal_n,
                    law_solve,
                    [x[0], x[1]],
                )?;
                let thermoelastic = if let (Some(slopes), Some(old), Some(continuation)) =
                    (thermal_slopes, old_thermodynamics, continuation) {
                    let new = source.law.thermodynamics_with_relaxation(
                        &mechanical.report.constitutive.state, slopes.equilibrium_n_m_k,
                        &slopes.relaxation_n_m_k)
                        .map_err(RigidLayerError::from)?;
                    let coefficient = continuation.coefficient_free_energy_change_j;
                    let entropy_energy_change_j = temperature * new.entropy_j_k
                        - history.law.temperature_kelvin() * old.entropy_j_k;
                    Some(RigidLayerThermoelasticReport {
                        coefficient_free_energy_change_j: coefficient,
                        free_energy_change_j: mechanical.report.constitutive.stored_change_j + coefficient,
                        mechanical_entropy_change_j_k: new.entropy_j_k - old.entropy_j_k,
                        entropy_energy_change_j,
                        internal_energy_change_j: new.internal_energy_j - old.internal_energy_j,
                        reversible_heat_j: -coefficient - entropy_energy_change_j,
                    })
                } else { None };
                Ok(ThirdBodyThermoelasticInput {
                    reversible_heat_j: thermoelastic.map_or(0.0, |value| value.reversible_heat_j),
                    mechanical_entropy_change_j_k: thermoelastic.map_or(0.0,
                        |value| value.mechanical_entropy_change_j_k),
                    heat: ThirdBodyHeatInput {
                    dissipated_j: mechanical.report.constitutive.dissipated_j,
                    conductivity_w_m_k: source.conductivity_w_m_k,
                    area_resistance_m2_k_w: source.area_resistance_m2_k_w,
                    transport: StepCoefficients {
                        sources_w: source.transport.sources_w,
                        conductances_w_k: source.transport.conductances_w_k,
                        payload: HeatStage {
                            mechanical,
                            law: source.law,
                            source: source.transport.payload,
                            thermoelastic,
                        },
                    },
                }})
            };
        heat_step.evaluate_thermoelastic(&energy, coefficients, entropy, &poll)
    };
    let problem = FallibleFiniteDifference::new(
        count + 2,
        controls.mechanical.difference_step_n,
        |x: &[f64]| {
            evaluate(x, false).map(|trial| {
                let mut residual = trial.payload().mechanical.report.force_residual_n.to_vec();
                residual.extend(trial.residual_j().iter().map(|value| value / scale));
                residual
            })
        },
    )
    .map_err(heat_numerical)?;
    let mut initial = force_guess_n.to_vec();
    initial.extend(energy_guess_j.iter().map(|value| value / scale));
    let construction = NewtonKrylovState::new(&problem, initial, controls.mechanical.newton);
    if let Some(error) = problem.take_failure() {
        return Err(heat_numerical(error));
    }
    let mut solver = construction?;
    let result = solver.run_controlled(
        &problem,
        controls.mechanical.maximum_newton_iterations,
        &mut || (*cancel.borrow_mut())() || problem.has_failed(),
    );
    if let Some(error) = problem.take_failure() {
        return Err(heat_numerical(error));
    }
    poll()?;
    if result.progress != SolverRunProgress::Complete {
        return Err(RigidLayerError::Interrupted(result.progress).into());
    }
    if !result.report.converged {
        return Err(RigidLayerError::NotConverged(result.report).into());
    }
    let trial = evaluate(&solver.x, controls.audit_entropy)?;
    let mechanical = &trial.payload().mechanical.report;
    mechanical.admit(controls.mechanical)?;
    let mechanical_storage = trial.payload().thermoelastic.map_or(
        mechanical.constitutive.stored_change_j, |value| value.internal_energy_change_j);
    let combined = mechanical.kinetic_energy_change_j
        + mechanical_storage
        + trial.report().thermal.stored_energy_change_j
        - trial.report().external_heat_j;
    if !combined.is_finite() || combined.abs() > controls.combined_energy_tolerance_j {
        return Err(RigidLayerError::Invalid(
            "original native layer and caloric storage do not balance",
        )
        .into());
    }
    let (thermal, staged, heat) = trial.accept(controls.heat, &poll)?;
    poll()?;
    Ok(HeatedRigidLayerCandidate {
        bodies: staged.mechanical.bodies,
        history: RigidLayerHistory {
            time_s: target_s,
            geometry: history.geometry,
            law: staged.law,
            deformation: staged.mechanical.report.constitutive.state.clone(),
        },
        thermal,
        report: HeatedRigidLayerReport {
            mechanical: staged.mechanical.report,
            heat,
            nonlinear: result.report,
            trials: trials.get(),
            source: staged.source,
            thermoelastic: staged.thermoelastic,
            combined_energy_residual_j: combined,
        },
    })
}

fn heat_numerical<E>(
    error: FiniteDifferenceError<HeatedRigidLayerError<E>>,
) -> HeatedRigidLayerError<E> {
    match error {
        FiniteDifferenceError::Source(error) => error,
        FiniteDifferenceError::Invalid(message) => RigidLayerError::Invalid(message).into(),
    }
}

struct Prepared {
    before: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    geometry: RigidLayerGeometry,
    axes: [Vec3; 2],
    separation: Vec3,
    /// Actual point force torque plus the interpolation's face couple, in each
    /// original principal frame. Cached once; no geometry work inside Newton.
    angular_loads: [[Vec3; 2]; 2],
    mobility: [[f64; 2]; 2],
    shear_before: [f64; 2],
    kinetic_before: f64,
}
struct Stage {
    bodies: [RigidBodyState; 2],
    report: RigidLayerReport,
}
impl Prepared {
    fn admit_thermal_route(
        &self,
        route: &ThirdBodyHeatRoute,
        tolerance: f64,
    ) -> Result<(), DepositionError> {
        for side in 0..2 {
            let area = route.patches[side].link.area_m2;
            if (area - self.geometry.area_m2).abs()
                > 128.0 * f64::EPSILON * area.max(self.geometry.area_m2)
            {
                return Err(DepositionError::InvalidInput(
                    "mechanical and thermal layer must have the same original area",
                ));
            }
            let (_, face) = half_distances(area, route.geometry[side], route.geometry[2])?;
            let point = self.before[side]
                .pose()
                .point_world_from_body(self.geometry.arms_body_m[side])
                .map_err(|_| DepositionError::Arithmetic)?;
            if norm(point.sub(vector(face))) > tolerance {
                return Err(DepositionError::InvalidInput(
                    "native layer face must equal the original thermal face centre",
                ));
            }
        }
        Ok(())
    }

    fn new(
        bodies: &[RigidBodyState; 2],
        properties: [MassProperties; 2],
        geometry: RigidLayerGeometry,
        tolerance: f64,
    ) -> Result<Self, RigidLayerError> {
        geometry.validate()?;
        let points = [
            bodies[0].point_kinematics(properties[0], geometry.arms_body_m[0])?,
            bodies[1].point_kinematics(properties[1], geometry.arms_body_m[1])?,
        ];
        let separation = points[1].point_world.sub(points[0].point_world);
        let normal = vector(geometry.frame.normal());
        if norm(separation.sub(normal.scale(geometry.thickness_m))) > tolerance {
            return Err(RigidLayerError::Invalid(
                "actual face separation differs from the original layer thickness and normal",
            ));
        }
        let axes = [
            vector(geometry.frame.longitudinal()),
            vector(geometry.frame.lateral()),
        ];
        let angular_loads =
            std::array::from_fn(|side| {
                std::array::from_fn(|axis| {
                    let point_torque = points[side]
                        .arm_world
                        .cross(axes[axis].scale(if side == 0 { 1.0 } else { -1.0 }));
                    bodies[side].pose().orientation().rotate_world_to_body(
                        point_torque.add(separation.cross(axes[axis]).scale(0.5)),
                    )
                })
            });
        let mobility = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                let mut value = (properties[0].mass().recip() + properties[1].mass().recip())
                    * axes[i].dot(axes[j]);
                for side in 0..2 {
                    value += angular_loads[side][i]
                        .dot(properties[side].angular_velocity_body(angular_loads[side][j]));
                }
                value
            })
        });
        if mobility.iter().flatten().any(|x| !x.is_finite())
            || angular_loads.iter().flatten().any(|x| !x.is_finite())
        {
            return Err(RigidLayerError::Invalid("non-finite native layer mobility"));
        }
        let relative = points[1]
            .point_velocity_world
            .sub(points[0].point_velocity_world)
            .sub(
                points[0]
                    .angular_velocity_world
                    .add(points[1].angular_velocity_world)
                    .scale(0.5)
                    .cross(separation),
            );
        let shear_before = axes.map(|axis| relative.dot(axis));
        let kinetic_before =
            bodies[0].kinetic_energy(properties[0])? + bodies[1].kinetic_energy(properties[1])?;
        Ok(Self {
            before: *bodies,
            properties,
            geometry,
            axes,
            separation,
            angular_loads,
            mobility,
            shear_before,
            kinetic_before,
        })
    }

    fn evaluate(
        &self,
        law: &ThirdBodyLayerContact,
        state: &ThirdBodyLayerState,
        dt: f64,
        normal: f64,
        solve: TangentialSolveConfig,
        force: [f64; 2],
    ) -> Result<Stage, RigidLayerError> {
        if force.iter().any(|x| !x.is_finite()) {
            return Err(RigidLayerError::Invalid(
                "finite candidate shear force required",
            ));
        }
        let world_force = self.axes[0]
            .scale(force[0])
            .add(self.axes[1].scale(force[1]));
        let face_couple = self.separation.cross(world_force).scale(0.5);
        let mut bodies = self.before;
        for side in 0..2 {
            let sign = if side == 0 { 1.0 } else { -1.0 };
            let angular = self.angular_loads[side][0]
                .scale(force[0])
                .add(self.angular_loads[side][1].scale(force[1]));
            bodies[side] = RigidBodyState::new(
                self.before[side].pose(),
                self.before[side]
                    .linear_momentum_world()
                    .add(world_force.scale(sign * dt)),
                self.before[side]
                    .angular_momentum_body()
                    .add(angular.scale(dt)),
            )?;
        }
        let points = [
            bodies[0].point_kinematics(self.properties[0], self.geometry.arms_body_m[0])?,
            bodies[1].point_kinematics(self.properties[1], self.geometry.arms_body_m[1])?,
        ];
        let relative = points[1]
            .point_velocity_world
            .sub(points[0].point_velocity_world)
            .sub(
                points[0]
                    .angular_velocity_world
                    .add(points[1].angular_velocity_world)
                    .scale(0.5)
                    .cross(self.separation),
            );
        let midpoint =
            std::array::from_fn(|i| 0.5 * (self.shear_before[i] + relative.dot(self.axes[i])));
        let slip = midpoint.map(|v| dt * v);
        let constitutive = law.advance(
            state,
            &TangentialStepInput {
                normal_n: normal,
                slip_increment_m: slip,
                dt_s: dt,
                contact_temperature_k: Some(law.temperature_kelvin()),
                nominal_pressure_pa: Some(normal / self.geometry.area_m2),
                compliance_m_per_n: 0.0,
                solve,
            },
        )?;
        let force_residual_n = std::array::from_fn(|i| force[i] - constitutive.traction_n[i]);
        let kinetic_after = bodies[0].kinetic_energy(self.properties[0])?
            + bodies[1].kinetic_energy(self.properties[1])?;
        let kinetic_energy_change_j = kinetic_after - self.kinetic_before;
        let mechanical_work_j = -(force[0] * slip[0] + force[1] * slip[1]);
        let total_energy_residual_j =
            kinetic_energy_change_j + constitutive.stored_change_j + constitutive.dissipated_j;
        let linear_momentum_change_ns = total_linear(&bodies).sub(total_linear(&self.before));
        let angular_momentum_change_nms = total_angular(&bodies).sub(total_angular(&self.before));
        Ok(Stage {
            bodies,
            report: RigidLayerReport {
                constitutive,
                applied_force_on_a_n: world_force,
                face_couple_nm: [face_couple; 2],
                midpoint_shear_velocity_mps: midpoint,
                inverse_shear_mass: self.mobility,
                force_residual_n,
                kinetic_energy_change_j,
                mechanical_work_j,
                work_energy_residual_j: kinetic_energy_change_j - mechanical_work_j,
                total_energy_residual_j,
                linear_momentum_change_ns,
                angular_momentum_change_nms,
            },
        })
    }
}
impl RigidLayerReport {
    fn admit(&self, controls: RigidLayerSolve) -> Result<(), RigidLayerError> {
        let residuals = [
            self.work_energy_residual_j,
            self.total_energy_residual_j,
            self.constitutive.closure_residual_j,
        ];
        if residuals
            .iter()
            .any(|x| !x.is_finite() || x.abs() > controls.energy_tolerance_j)
            || self
                .force_residual_n
                .iter()
                .any(|x| !x.is_finite() || x.abs() > controls.force_tolerance_n)
            || norm(self.linear_momentum_change_ns) > controls.linear_momentum_tolerance_ns
            || norm(self.angular_momentum_change_nms) > controls.angular_momentum_tolerance_nms
            || !self.linear_momentum_change_ns.is_finite()
            || !self.angular_momentum_change_nms.is_finite()
        {
            return Err(RigidLayerError::Invalid(
                "original finite-layer force, work or momentum balance failed",
            ));
        }
        Ok(())
    }
}
fn duration(history: &RigidLayerHistory, target: f64) -> Result<f64, RigidLayerError> {
    let dt = target - history.time_s;
    if !target.is_finite() || !dt.is_finite() || dt <= 0.0 {
        Err(RigidLayerError::Invalid(
            "strictly advancing finite layer clock required",
        ))
    } else {
        Ok(dt)
    }
}
fn spend(trials: &Cell<usize>, maximum: usize) -> Result<(), RigidLayerError> {
    if trials.get() >= maximum {
        return Err(RigidLayerError::TrialBudget);
    }
    trials.set(trials.get() + 1);
    Ok(())
}
fn numerical(error: FiniteDifferenceError<RigidLayerError>) -> RigidLayerError {
    match error {
        FiniteDifferenceError::Source(error) => error,
        FiniteDifferenceError::Invalid(message) => RigidLayerError::Invalid(message),
    }
}
fn vector(v: [f64; 3]) -> Vec3 {
    Vec3::new(v[0], v[1], v[2])
}
fn norm(v: Vec3) -> f64 {
    v.x.hypot(v.y).hypot(v.z)
}
fn total_linear(bodies: &[RigidBodyState; 2]) -> Vec3 {
    bodies[0]
        .linear_momentum_world()
        .add(bodies[1].linear_momentum_world())
}
fn total_angular(bodies: &[RigidBodyState; 2]) -> Vec3 {
    bodies.iter().fold(Vec3::ZERO, |sum, body| {
        sum.add(
            body.pose()
                .position_world()
                .cross(body.linear_momentum_world()),
        )
        .add(
            body.pose()
                .orientation()
                .rotate_body_to_world(body.angular_momentum_body()),
        )
    })
}
