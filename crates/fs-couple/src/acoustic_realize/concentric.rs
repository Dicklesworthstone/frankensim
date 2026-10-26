//! Continuing layered-string motion and per-layer heat on the native owners.
//!
//! Perfectly bonded circular annuli keep their own source, temperature, caloric
//! cell, free strain and bending loss. The ordinary string modal propagator,
//! gas resistance and pressure observer consume their section integrals. Heat
//! uses the existing finite thermal network and its backward-Euler residual.
//! There is no equivalent material card or independent vibration/heat solver.
//!
//! A step advances motion once at the accepted material state, then solves
//! layer temperatures at frozen new displacement/momentum. Axial and modal
//! Helmholtz energy and T*s supply the reversible exchange into EACH layer.
//! This is a first-order partitioned method; refine physical time and retained
//! modes independently of nonlinear residual tolerances. Fixed circular
//! geometry, constant density, scalar alpha, affine E(T), solid caloric charts
//! and Kelvin-Voigt bending are explicit limits. Helical slip, evolving geometry,
//! melting, nonlinear stretching and a coupled instrument body are not supplied.

/// Finite source gas, reciprocal surface heat and retained acoustic history.
pub mod finite_gas;

use std::cell::Cell;
use std::collections::{BTreeMap, BTreeSet};
use std::convert::Infallible;

use super::{AcousticRealizeError, LinearMember, MaterialStringFrame, ModalAcousticMode,
    ModalAcousticState, Pluck, PrestressedString, StringDissipation, StringPointForce,
    StringPointMotion, linear_accelerations, linear_string_member, next_string_epoch,
    partition_string_dissipation, validate_external_gas, validate_incremental_string,
    validate_string, validate_string_station, string_point_participation};
use crate::material_assembly::concentric::with_concentric_kelvin_voigt_section;
use crate::string_specimen::{KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY, StringPrestress};
use fs_exec::Cx;
use fs_material::conductor::{ConductorError, OhmicDrive};
use fs_material::gas::GasState;
use fs_material::phase::{MAX_HEAT_CAPACITY_ENTHALPY_KNOTS, PhaseStateError, ResolvedHeatCapacityCurve};
use fs_material::thermochemical_gas::energy::GasEnergyError;
use fs_material::section::ConcentricKelvinVoigtSection;
use fs_material::section::transport::ElectrothermalSectionResponse;
use fs_material::section_temperature::thermoelastic::{ConcentricThermoelasticLaw,
    ConcentricThermoelasticState};
use fs_material::state_point::{DENSITY_PROPERTY, MaterialStatePointError,
    SPECIFIC_HEAT_CAPACITY_PROPERTY};
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport,
    SolverRunProgress};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation, EnergyAudit,
    backward_euler::{BackwardEulerStep, StepCoefficients}};

/// One original source-backed caloric cell for one annulus, in spatial order.
#[derive(Clone, Copy)]
pub struct ConcentricStringThermalCell<'a> {
    pub cell: CellId,
    pub caloric: &'a ResolvedHeatCapacityCurve,
}

/// Explicit finite thermal work and acceptance limits; no default physics.
#[derive(Clone, Copy, Debug)]
pub struct ConcentricStringThermalSolve {
    pub maximum_cells: usize,
    pub maximum_thermal_trials: usize,
    pub maximum_newton_iterations: usize,
    pub difference_step_j: f64,
    pub maximum_temperature_increment_k: f64,
    pub cell_energy_tolerance_j: f64,
    pub thermal_energy_tolerance_j: f64,
    pub total_energy_tolerance_j: f64,
    pub newton: NewtonKrylovConfig,
}

/// Failure is transactional across material, motion, observer and thermal state.
#[derive(Debug)]
pub enum ConcentricStringError<E = Infallible> {
    Invalid(&'static str),
    Material(MaterialStatePointError),
    Electrical(ConductorError),
    Acoustic(AcousticRealizeError),
    Thermal(DepositionError),
    Caloric(PhaseStateError),
    Gas(GasEnergyError),
    ExchangeOutsideBracket { lower_residual_j: f64, upper_residual_j: f64 },
    ExchangeNotConverged { iterations: usize, residual_j: f64 },
    Source(E),
    Newton(NewtonError),
    NotConverged(NewtonReport),
    Interrupted(SolverRunProgress),
    TrialBudget,
}
impl<E> From<MaterialStatePointError> for ConcentricStringError<E> {
    fn from(e: MaterialStatePointError) -> Self { Self::Material(e) }
}
impl<E> From<ConductorError> for ConcentricStringError<E> {
    fn from(e: ConductorError) -> Self { Self::Electrical(e) }
}
impl<E> From<AcousticRealizeError> for ConcentricStringError<E> {
    fn from(e: AcousticRealizeError) -> Self { Self::Acoustic(e) }
}
impl<E> From<DepositionError> for ConcentricStringError<E> {
    fn from(e: DepositionError) -> Self { Self::Thermal(e) }
}
impl<E> From<PhaseStateError> for ConcentricStringError<E> {
    fn from(e: PhaseStateError) -> Self { Self::Caloric(e) }
}
impl<E> From<GasEnergyError> for ConcentricStringError<E> {
    fn from(e: GasEnergyError) -> Self { Self::Gas(e) }
}
impl<E> From<NewtonError> for ConcentricStringError<E> {
    fn from(e: NewtonError) -> Self { Self::Newton(e) }
}
impl<E: core::fmt::Display> core::fmt::Display for ConcentricStringError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(s) => write!(f, "concentric thermal string: {s}"),
            Self::Material(e) => write!(f, "concentric thermal material: {e}"),
            Self::Electrical(e) => write!(f, "concentric thermal electrical transport: {e}"),
            Self::Acoustic(e) => write!(f, "concentric thermal mechanics: {e}"),
            Self::Thermal(e) => write!(f, "concentric thermal network: {e}"),
            Self::Caloric(e) => write!(f, "concentric caloric source: {e}"),
            Self::Gas(e) => write!(f, "concentric finite gas: {e}"),
            Self::ExchangeOutsideBracket { lower_residual_j, upper_residual_j } =>
                write!(f, "concentric gas exchange is not bracketed: {lower_residual_j}, {upper_residual_j} J"),
            Self::ExchangeNotConverged { iterations, residual_j } =>
                write!(f, "concentric gas exchange did not converge in {iterations} iterations: {residual_j} J"),
            Self::Source(e) => write!(f, "concentric thermal source: {e}"),
            Self::Newton(e) => write!(f, "concentric thermal solve: {e:?}"),
            Self::NotConverged(r) => write!(f, "concentric thermal solve did not converge: {r:?}"),
            Self::Interrupted(p) => write!(f, "concentric thermal solve interrupted: {p:?}"),
            Self::TrialBudget => f.write_str("concentric thermal trial budget exhausted"),
        }
    }
}
impl<E: core::fmt::Debug + core::fmt::Display> std::error::Error for ConcentricStringError<E> {}

/// Actual per-layer thermodynamic transfers for the accepted sample [J].
#[derive(Clone, Debug, PartialEq)]
pub struct ConcentricStringLayerExchange {
    pub cell: CellId,
    pub material_heat_j: f64,
    pub axial_free_energy_change_j: f64,
    pub modal_parameter_work_j: f64,
    pub mechanical_entropy_energy_change_j: f64,
    pub reversible_heat_j: f64,
}

/// One accepted physical clock and full stored-energy account.
#[derive(Debug)]
pub struct ConcentricStringThermalFrame<S> {
    pub vibration: MaterialStringFrame,
    pub thermal: EnergyAudit,
    pub source: S,
    pub nonlinear: NewtonReport,
    pub thermal_trials: usize,
    pub maximum_cell_residual_j: f64,
    pub layers: Vec<ConcentricStringLayerExchange>,
    pub external_heat_j: f64,
    pub mechanical_internal_energy_before_j: f64,
    pub mechanical_internal_energy_after_j: f64,
    pub total_energy_residual_j: f64,
    pub roundoff_tolerance_j: f64,
}

/// Actual endpoint electrical response and the supply delivered during a step.
/// Joule heat is evaluated implicitly at the accepted layer temperatures.
#[derive(Debug)]
pub struct ConcentricStringOhmicSource<S> {
    pub drive: OhmicDrive,
    pub electrical: ElectrothermalSectionResponse,
    pub joule_heat_j: f64,
    pub layer_joule_heat_j: Vec<f64>,
    /// Authored boundary and independent callback supply, excluding Joule heat.
    pub other_external_heat_j: f64,
    pub transport: S,
}

/// The native string mode/observer history with original layered constitutive data.
#[derive(Clone)]
pub struct ThermalConcentricStringRuntime<'a> {
    law: ConcentricThermoelasticLaw<'a>,
    bindings: Vec<ConcentricStringThermalCell<'a>>,
    template: PrestressedString,
    reference_length_m: f64,
    strain_limit: f64,
    source: ConcentricThermoelasticState,
    section: ConcentricKelvinVoigtSection,
    string: PrestressedString,
    member: LinearMember,
    gas: GasState,
    listener_m: f64,
    sample_rate_hz: u32,
    initial_time_s: f64,
    accepted_samples: u64,
}

#[derive(Clone, Copy)]
struct LayerEnergy {
    axial_j: f64,
    modal_j: f64,
    entropy_j_k: f64,
}

impl<'a> ThermalConcentricStringRuntime<'a> {
    /// Initialize a fixed-support layered string from the SAME thermal cells.
    ///
    /// The source law retains independent layer cards and expansion references.
    /// `reference_length_m` is their common stress-free support separation.
    /// Every cell must carry exactly that annulus's mass, volume and original
    /// source caloric law. Density and source pins remain invariant. A node can
    /// represent a whole annulus; radial gradients inside it are not resolved.
    /// The ordinary source-bound Kelvin-Voigt selector is used for every layer.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        cx: &Cx<'_>, template: PrestressedString,
        law: ConcentricThermoelasticLaw<'a>,
        bindings: Vec<ConcentricStringThermalCell<'a>>,
        thermal: &DepositionSimulation<'_>, reference_length_m: f64,
        strain_limit: f64, pluck: Option<Pluck>, gas: GasState,
        listener_m: f64, sample_rate_hz: u32,
    ) -> Result<Self, ConcentricStringError> {
        checkpoint(cx)?;
        validate_external_gas(&gas)?;
        if !positive(reference_length_m) || !positive(listener_m)
            || bindings.len() != law.reference().layers().len()
            || !thermal.time_s().is_finite()
        { return Err(invalid("invalid support geometry, listener or layer/cell count")); }
        let temperatures = temperatures(&bindings, thermal)?;
        let strain = finite((template.length_m-reference_length_m)/reference_length_m)?;
        let source = law.at_temperatures_k(&temperatures, strain, strain_limit)?;
        let (section, string) = bind_string(&template, &source)?;
        validate_incremental_string(&string)?;
        validate_string(&string, pluck, None)?;
        let member = linear_string_member(&string, pluck, None, &gas,
            listener_m, sample_rate_hz, 1.0)?;
        let result = Self { law, bindings, template, reference_length_m, strain_limit,
            source, section, string, member, gas, listener_m, sample_rate_hz,
            initial_time_s: thermal.time_s(), accepted_samples: 0 };
        result.admit_cells(cx, thermal)?;
        result.admit_capacity(cx, &result.source, result.states())?;
        checkpoint(cx)?;
        Ok(result)
    }

    /// Accepted source bundle and separate layer thermal/elastic responses.
    pub const fn source(&self) -> &ConcentricThermoelasticState { &self.source }
    pub const fn section(&self) -> &ConcentricKelvinVoigtSection { &self.section }
    pub const fn string(&self) -> &PrestressedString { &self.string }
    pub const fn accepted_samples(&self) -> u64 { self.accepted_samples }
    pub fn states(&self) -> &[ModalAcousticState] { self.member.model.states() }
    pub fn modes(&self) -> &[ModalAcousticMode] { self.member.model.modes() }
    pub const fn gas(&self) -> &GasState { &self.gas }
    pub fn thermal_bindings(&self) -> &[ConcentricStringThermalCell<'a>] { &self.bindings }

    /// Settle physical point loads before the simulation window, then hold or
    /// release them through the ordinary force projection. Only an unexcited
    /// runtime without accepted samples or observer history can be initialized.
    /// The native modal equilibrium owner supplies the displacement; existing
    /// motion is never replaced and no loading work is invented inside the run.
    pub fn initialize_point_load_equilibrium(&mut self, cx: &Cx<'_>, loads: &[StringPointForce])
        -> Result<(), ConcentricStringError>
    {
        checkpoint(cx)?;
        if self.accepted_samples != 0 || self.member.observer.previous_acceleration.is_some()
            || self.states().iter().any(|state| state.displacement_m_sqrt_kg != 0.0
                || state.velocity_m_sqrt_kg_per_s != 0.0)
        { return Err(invalid("settled loads require an unexcited layered string before stepping")); }
        let forces = self.project_point_forces(cx, loads)?;
        let mut member = self.member.clone();
        member.model.initialize_static_equilibrium(&forces).map_err(AcousticRealizeError::Modal)?;
        self.admit_capacity(cx, &self.source, member.model.states())?;
        member.observer.previous_acceleration = Some(vec![0.0; forces.len()]);
        checkpoint(cx)?;
        self.member = member;
        Ok(())
    }

    /// Project held physical point forces using the accepted mass-normalized basis.
    pub fn project_point_forces(&self, cx: &Cx<'_>, loads: &[StringPointForce])
        -> Result<Vec<f64>, ConcentricStringError>
    {
        if loads.len() > MAX_POINT_FORCES { return Err(invalid("too many point forces")); }
        checkpoint(cx)?;
        let normalization = fs_math::det::sqrt(0.5*self.string.lin_density_kg_m*self.string.length_m);
        let mut forces = vec![0.0; self.string.n_modes];
        for load in loads {
            checkpoint(cx)?;
            validate_string_station(load.station_frac)?;
            if !load.force_n.is_finite() { return Err(invalid("invalid physical point force")); }
            for (mode, force) in forces.iter_mut().enumerate() {
                if mode % 64 == 0 { checkpoint(cx)?; }
                let phi = string_point_participation(mode, load.station_frac, normalization);
                *force = finite(*force+phi*load.force_n)?;
            }
        }
        Ok(forces)
    }

    /// Physical motion from the same retained sine coordinates used by forces.
    pub fn point_motion(&self, cx: &Cx<'_>, station_frac: f64) -> Result<StringPointMotion, ConcentricStringError> {
        checkpoint(cx)?;
        validate_string_station(station_frac)?;
        let normalization = fs_math::det::sqrt(0.5*self.string.lin_density_kg_m*self.string.length_m);
        let (mut displacement_m, mut velocity_m_s) = (0.0, 0.0);
        for (mode, state) in self.states().iter().enumerate() {
            if mode % 64 == 0 { checkpoint(cx)?; }
            let phi = string_point_participation(mode, station_frac, normalization);
            displacement_m = finite(displacement_m+phi*state.displacement_m_sqrt_kg)?;
            velocity_m_s = finite(velocity_m_s+phi*state.velocity_m_sqrt_kg_per_s)?;
        }
        Ok(StringPointMotion { epoch: self.accepted_samples, displacement_m, velocity_m_s })
    }

    /// Kinetic + axial/modal Helmholtz + sum(T_i*s_i), excluding caloric storage.
    pub fn mechanical_internal_energy_j(&self, cx: &Cx<'_>) -> Result<f64, ConcentricStringError> {
        internal_energy(cx, &self.source, self.states(), self.reference_length_m, self.string.length_m)
    }

    /// Continue the same motion and finite thermal network atomically.
    ///
    /// `resolve` supplies independent sources and original-link conductances at
    /// the trial endpoint. Electrical adapters can evaluate their original
    /// conductor/section there. Do not add the solid bending loss or reversible
    /// exchange: this owner routes them once to their actual material layers.
    /// Other thermal cells and boundary reservoirs stay in the original network.
    /// The callback must not mutate external physical history during a trial.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_thermal_network<'t, S, E>(
        &mut self, cx: &Cx<'_>, thermal: &mut DepositionSimulation<'t>,
        controls: ConcentricStringThermalSolve, energy_guess_j: &[f64], forces: &[f64],
        resolve: impl Fn(&DepositionSimulation<'t>) -> Result<StepCoefficients<S>, E>,
    ) -> Result<ConcentricStringThermalFrame<S>, ConcentricStringError<E>> {
        self.step_with_resolved_thermal_network(cx, thermal, controls, energy_guess_j, forces,
            |endpoint, _, _, _| resolve(endpoint).map_err(ConcentricStringError::Source))
    }

    /// Drive current or voltage through the original layered conductor and
    /// continue its heat, thermal prestress, vibration and sound in one solve.
    ///
    /// The existing parallel DC owner supplies each annulus's current and
    /// Joule power from its OWN resistivity/conductivity at the SAME trial
    /// temperature used by mechanics. Total current is shared according to
    /// actual conductance; it is never applied independently to every layer.
    /// The original material bundles must carry all required electrical and
    /// thermal transport properties. No material name or missing coefficient
    /// supplies a default. Radial currents, terminal contact losses and skin
    /// effects are outside this common-equipotential-terminal reduction.
    ///
    /// `resolve` supplies only OTHER heat and original-link conductances. Joule
    /// heat, solid bending loss and reversible exchange are each routed once.
    /// `frame.external_heat_j` includes the electrical supply; the payload
    /// separates it from other sources. End heat flow is zero in this
    /// axially uniform layer model; it is not added as a net heat source.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_ohmic_heating<'t, S, E>(
        &mut self, cx: &Cx<'_>, thermal: &mut DepositionSimulation<'t>,
        controls: ConcentricStringThermalSolve, energy_guess_j: &[f64], forces: &[f64],
        drive: OhmicDrive,
        resolve: impl Fn(&DepositionSimulation<'t>) -> Result<StepCoefficients<S>, E>,
    ) -> Result<ConcentricStringThermalFrame<ConcentricStringOhmicSource<S>>, ConcentricStringError<E>> {
        let length_m = self.string.length_m;
        let cells: Vec<_> = self.bindings.iter().map(|binding| binding.cell).collect();
        self.step_with_resolved_thermal_network(cx, thermal, controls, energy_guess_j, forces,
            |endpoint, source, dt, _| add_ohmic_heat(cx, endpoint, source, length_m,
                dt, &cells, drive, controls.total_energy_tolerance_j,
                resolve(endpoint).map_err(ConcentricStringError::Source)?))
    }

    // The electrical adapter receives exactly the source already resolved for
    // this thermal/mechanical candidate; it must not sample another material epoch.
    #[allow(clippy::too_many_arguments)]
    fn step_with_resolved_thermal_network<'t, S, E>(
        &mut self, cx: &Cx<'_>, thermal: &mut DepositionSimulation<'t>,
        controls: ConcentricStringThermalSolve, energy_guess_j: &[f64], forces: &[f64],
        resolve: impl Fn(&DepositionSimulation<'t>, &ConcentricThermoelasticState, f64, &StringDissipation)
            -> Result<StepCoefficients<S>, ConcentricStringError<E>>,
    ) -> Result<ConcentricStringThermalFrame<S>, ConcentricStringError<E>> {
        checkpoint(cx)?;
        controls.newton.validate()?;
        let next_sample = next_string_epoch(self.accepted_samples)?;
        let current_s = finite(self.initial_time_s+self.accepted_samples as f64/f64::from(self.sample_rate_hz))?;
        let target_s = finite(self.initial_time_s+next_sample as f64/f64::from(self.sample_rate_hz))?;
        if thermal.time_s() != current_s || controls.maximum_cells == 0
            || thermal.cells().len() > controls.maximum_cells
            || thermal.cells().len() != energy_guess_j.len()
            || controls.maximum_thermal_trials == 0 || controls.maximum_newton_iterations == 0
            || !positive(controls.maximum_temperature_increment_k)
            || [controls.cell_energy_tolerance_j, controls.thermal_energy_tolerance_j,
                controls.total_energy_tolerance_j].iter().any(|v| !v.is_finite() || *v < 0.0)
        { return Err(invalid("thermal clocks, controls or energy vector are incompatible")); }
        self.admit_cells(cx, thermal)?;
        let thermal_checkpoint = || cx.checkpoint().map_err(|_| DepositionError::Cancelled);
        let before_energy = internal_energy(cx, &self.source, self.states(), self.reference_length_m, self.string.length_m)?;
        let before_layers = layer_energies(cx, &self.source, self.states(), self.reference_length_m, self.string.length_m)?;
        // Actual motion advances ONCE, at accepted material temperature. All
        // nonlinear thermal trials rebind this same candidate without stepping it.
        let mut advanced = self.member.clone();
        let mut acoustic = advanced.model.step(forces).map_err(AcousticRealizeError::Modal)?;
        let dissipation = partition_string_dissipation(&acoustic, &advanced.loss_parts)?;
        let fixed_layers = layer_energies(cx, &self.source, advanced.model.states(), self.reference_length_m, self.string.length_m)?;
        let initial_acceleration = self.member.observer.previous_acceleration.is_none()
            .then(|| linear_accelerations(&self.member.model, forces));
        let material_loss = layer_losses(&self.section, &dissipation)?;
        let step = BackwardEulerStep::new(thermal, target_s, thermal_checkpoint)?;
        let dt = step.duration_s();
        // The native mechanical clock is fixed; an unrepresentable thermal
        // endpoint must not silently use a different physical interval.
        let clock_tolerance = 16.0*f64::EPSILON*(target_s.abs()+current_s.abs()+dt);
        if (dt-advanced.model.sample_period_s()).abs() > clock_tolerance {
            return Err(invalid("thermal and mechanical sample durations differ"));
        }
        let trials = Cell::new(0usize);
        let evaluate = |x: &[f64]| {
            step.evaluate(x, |endpoint| {
                checkpoint(cx)?;
                if trials.get() >= controls.maximum_thermal_trials { return Err(ConcentricStringError::TrialBudget); }
                trials.set(trials.get()+1);
                let temperatures = temperatures(&self.bindings, endpoint)?;
                if temperatures.iter().zip(&self.source.layers).any(|(t, layer)|
                    (*t-layer.temperature_k).abs() > controls.maximum_temperature_increment_k) {
                    return Err(invalid("layer temperature increment exceeds its physical step budget"));
                }
                let strain = (self.string.length_m-self.reference_length_m)/self.reference_length_m;
                let source = self.law.at_temperatures_k(&temperatures, strain, self.strain_limit)?;
                let (section, string) = bind_string(&self.template, &source)?;
                if string.lin_density_kg_m.to_bits() != self.string.lin_density_kg_m.to_bits() {
                    return Err(invalid("thermal update changed the fixed-mass modal basis"));
                }
                self.admit_capacity(cx, &source, advanced.model.states())?;
                let mut member = linear_string_member(&string, None, None, &self.gas,
                    self.listener_m, self.sample_rate_hz, 1.0)?;
                member.model.restore_states(advanced.model.states()).map_err(AcousticRealizeError::Modal)?;
                member.observer.previous_acceleration = self.member.observer.previous_acceleration.clone();
                let after_layers = layer_energies(cx, &source, member.model.states(), self.reference_length_m, string.length_m)?;
                let mut exchanges = Vec::with_capacity(self.bindings.len());
                let StepCoefficients { mut sources_w, conductances_w_k, payload } =
                    resolve(endpoint, &source, dt, &dissipation)?;
                let external_heat_j = finite(dt*endpoint.external_step_power_w(&sources_w,
                    &conductances_w_k, thermal_checkpoint)?)?;
                for (i, binding) in self.bindings.iter().enumerate() {
                    checkpoint(cx)?;
                    let axial = finite(after_layers[i].axial_j-before_layers[i].axial_j)?;
                    let modal = finite(after_layers[i].modal_j-fixed_layers[i].modal_j)?;
                    let entropy = finite(source.layers[i].temperature_k*after_layers[i].entropy_j_k
                        -self.source.layers[i].temperature_k*before_layers[i].entropy_j_k)?;
                    let reversible = finite(-axial-modal-entropy)?;
                    deposit(&mut sources_w, binding.cell, finite(material_loss[i]+reversible)?, dt,
                        controls.total_energy_tolerance_j)?;
                    exchanges.push(ConcentricStringLayerExchange { cell: binding.cell,
                        material_heat_j: material_loss[i], axial_free_energy_change_j: axial,
                        modal_parameter_work_j: modal, mechanical_entropy_energy_change_j: entropy,
                        reversible_heat_j: reversible });
                }
                Ok::<_, ConcentricStringError<E>>(StepCoefficients { sources_w, conductances_w_k,
                    payload: (source, section, string, member, exchanges, external_heat_j, payload) })
            }, thermal_checkpoint)
        };
        let problem = FallibleFiniteDifference::new(step.cell_ids().len(), controls.difference_step_j,
            |x: &[f64]| evaluate(x).map(|t| t.residual_j().to_vec())).map_err(numerical)?;
        let construction = NewtonKrylovState::new(&problem, energy_guess_j.to_vec(), controls.newton);
        if let Some(error) = problem.take_failure() { return Err(numerical(error)); }
        let mut solver = construction?;
        let result = solver.run_controlled(&problem, controls.maximum_newton_iterations,
            &mut || cx.checkpoint().is_err() || problem.has_failed());
        if let Some(error) = problem.take_failure() { return Err(numerical(error)); }
        checkpoint(cx)?;
        if result.progress != SolverRunProgress::Complete { return Err(ConcentricStringError::Interrupted(result.progress)); }
        if !result.report.converged { return Err(ConcentricStringError::NotConverged(result.report)); }
        let trial = evaluate(&solver.x)?;
        let maximum_cell_residual_j = trial.maximum_cell_residual_j();
        let (next, (source, section, string, mut member, layers, external_heat_j, payload), audit) =
            trial.accept(controls.cell_energy_tolerance_j, controls.thermal_energy_tolerance_j, thermal_checkpoint)?;
        acoustic.modal_energy_j = member.model.modes().iter().zip(member.model.states())
            .map(|(mode, state)| 0.5*(state.velocity_m_sqrt_kg_per_s.powi(2)
                +(mode.angular_frequency_rad_s*state.displacement_m_sqrt_kg).powi(2))).collect();
        acoustic.total_modal_energy_j = acoustic.modal_energy_j.iter().sum();
        acoustic.observer_pressure_pa = member.observer.observe(
            linear_accelerations(&member.model, forces), initial_acceleration, member.model.sample_period_s())?;
        let after_energy = internal_energy(cx, &source, member.model.states(), self.reference_length_m, string.length_m)?;
        let residual = finite(acoustic.input_work_j+external_heat_j
            -(after_energy-before_energy+audit.stored_energy_change_j+dissipation.air_loss_j
                +dissipation.authored_loss_j+dissipation.roundoff_residual_j))?;
        let scale = before_energy.abs()+after_energy.abs()+acoustic.input_work_j.abs()
            +external_heat_j.abs()+audit.stored_energy_change_j.abs()
            +layers.iter().map(|l| l.material_heat_j.abs()+l.axial_free_energy_change_j.abs()
                +l.modal_parameter_work_j.abs()+l.mechanical_entropy_energy_change_j.abs()).sum::<f64>();
        let roundoff = finite(dissipation.roundoff_tolerance_j+256.0*f64::EPSILON*scale)?;
        if residual.abs() > controls.total_energy_tolerance_j+roundoff {
            return Err(invalid("layered modal, thermoelastic and caloric energy do not close"));
        }
        checkpoint(cx)?;
        self.source = source;
        self.section = section;
        self.string = string;
        self.member = member;
        self.accepted_samples = next_sample;
        *thermal = next;
        Ok(ConcentricStringThermalFrame { vibration: MaterialStringFrame {
                epoch: next_sample, acoustic, dissipation }, thermal: audit, source: payload,
            nonlinear: result.report, thermal_trials: trials.get(), maximum_cell_residual_j,
            layers, external_heat_j, mechanical_internal_energy_before_j: before_energy,
            mechanical_internal_energy_after_j: after_energy, total_energy_residual_j: residual,
            roundoff_tolerance_j: roundoff })
    }

    fn admit_cells<E>(&self, cx: &Cx<'_>, thermal: &DepositionSimulation<'_>) -> Result<(), ConcentricStringError<E>> {
        let mut ids = BTreeSet::new();
        for ((binding, layer), geometry) in self.bindings.iter()
            .zip(self.source.section.section().layers()).zip(self.source.section.section().layer_integrals()) {
            checkpoint(cx)?;
            if !ids.insert(binding.cell) { return Err(invalid("layer thermal cells must be distinct")); }
            let cell = thermal.cells().get(&binding.cell).ok_or(invalid("missing layer thermal cell"))?;
            let mass = geometry.line_density_kg_m()*self.string.length_m;
            let volume = geometry.area_m2()*self.string.length_m;
            if !same_inventory(cell.mass_kg, mass) || !same_inventory(cell.geometry.volume_m3, volume)
                || !thermal.uses_material_law(cell.material_id, binding.caloric.curve())
                || binding.caloric.curve().material_card_identity() != layer.material.card_identity()
            { return Err(invalid("thermal cell differs from its original material, annular mass, volume or caloric law")); }
            let temperature = thermal.temperature_k(binding.cell)?;
            let actual = layer.material.query_point().iter().find(|(key, _)| key == "T")
                .ok_or(invalid("layer source lacks canonical absolute temperature"))?.1;
            if temperature.to_bits() != actual.to_bits() { return Err(invalid("thermal and mechanical layer temperatures differ")); }
            for pair in binding.caloric.receipts().chunks_exact(2) {
                for (receipt, name) in pair.iter().zip([SPECIFIC_HEAT_CAPACITY_PROPERTY, DENSITY_PROPERTY]) {
                    let property = layer.material.property(name).ok_or(invalid("layer lacks its original caloric property"))?;
                    if receipt.selected != property.answer().receipt.selected
                        || receipt.axis_quantities != *layer.material.axis_quantities()
                        || !receipt.query_point.iter().filter(|(key, _)| key != "T")
                            .eq(layer.material.query_point().iter().filter(|(key, _)| key != "T")) {
                        return Err(invalid("caloric and mechanical layer source claims or fixed conditions differ"));
                    }
                }
            }
        }
        Ok(())
    }

    fn admit_capacity<E>(&self, cx: &Cx<'_>, source: &ConcentricThermoelasticState, states: &[ModalAcousticState])
        -> Result<(), ConcentricStringError<E>>
    {
        let gradient = modal_gradient(cx, states, self.string.length_m)?;
        let mu = self.string.lin_density_kg_m;
        for (binding, layer) in self.bindings.iter().zip(&source.layers) {
            checkpoint(cx)?;
            let knots = binding.caloric.curve().knots();
            if knots.len() < 2 || knots.len() > MAX_HEAT_CAPACITY_ENTHALPY_KNOTS
                || knots.iter().any(|k| k.liquid_mass_fraction != 0.0 || k.bulk_density_kg_m3 != layer.density_kg_m3) {
                return Err(invalid("layer heat requires a bounded single-solid constant-density chart"));
            }
            let a = layer.alpha_per_k;
            let b = layer.young_modulus_derivative_pa_k;
            let curvature = finite(self.reference_length_m*layer.free_energy_temperature_second_derivative_j_m_k2
                -layer.area_m2/mu*b*a*gradient)?;
            let slope = finite(3.0*self.reference_length_m*layer.area_m2*b*a*a)?;
            let penalty_at = |t: f64| finite::<E>(t*(curvature+slope*(t-layer.temperature_k)));
            for pair in knots.windows(2) {
                checkpoint(cx)?;
                let (lo, hi) = (pair[0].temperature_k, pair[1].temperature_k);
                let capacity = finite(layer.density_kg_m3*layer.area_m2*self.string.length_m
                    *(pair[1].specific_enthalpy_j_kg-pair[0].specific_enthalpy_j_kg)/(hi-lo))?;
                let mut penalty = penalty_at(lo)?.max(penalty_at(hi)?);
                if slope != 0.0 {
                    let stationary = 0.5*(layer.temperature_k-curvature/slope);
                    if stationary.is_finite() && stationary > lo && stationary < hi {
                        penalty = penalty.max(penalty_at(stationary)?);
                    }
                }
                if !positive(hi-lo) || !positive(capacity) || !positive(capacity-penalty) {
                    return Err(invalid("layer caloric capacity does not exceed its axial and modal thermoelastic penalty"));
                }
            }
        }
        Ok(())
    }
}

const MAX_POINT_FORCES: usize = 4096;

fn bind_string<E>(template: &PrestressedString, source: &ConcentricThermoelasticState)
    -> Result<(ConcentricKelvinVoigtSection, PrestressedString), ConcentricStringError<E>>
{
    // Each layer books its own modal heat from its bending viscosity. A tension
    // loss factor would dissipate energy that no layer owns, so it refuses
    // until that heat has an owner. Zero is the absent factor.
    if template.tension_loss_factor.is_some_and(|eta| eta != 0.0) {
        return Err(invalid("a layered thermoelastic string has no owner for tension-loss heat"));
    }
    let section = source.section.section().clone()
        .with_kelvin_voigt_bending(KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY)?;
    let force: f64 = source.layers.iter().map(|l| l.young_modulus_pa*l.area_m2*l.elastic_strain).sum();
    let string = with_concentric_kelvin_voigt_section(template.clone(), &section,
        StringPrestress::FixedTension(force))?;
    Ok((section, string))
}

fn temperatures<E>(bindings: &[ConcentricStringThermalCell<'_>], thermal: &DepositionSimulation<'_>)
    -> Result<Vec<f64>, ConcentricStringError<E>>
{ bindings.iter().map(|b| thermal.temperature_k(b.cell).map_err(ConcentricStringError::Thermal)).collect() }

fn layer_energies<E>(cx: &Cx<'_>, source: &ConcentricThermoelasticState, states: &[ModalAcousticState],
    reference_length_m: f64, length_m: f64) -> Result<Vec<LayerEnergy>, ConcentricStringError<E>>
{
    let mu = source.section.section().line_density_kg_m();
    let mut result = Vec::with_capacity(source.layers.len());
    for layer in &source.layers {
        checkpoint(cx)?;
        let mut modal_j = 0.0;
        let mut modal_entropy = 0.0;
        let e = layer.elastic_strain;
        let modulus = layer.young_modulus_pa;
        let b = layer.young_modulus_derivative_pa_k;
        for (mode, state) in states.iter().enumerate() {
            if mode % 64 == 0 { checkpoint(cx)?; }
            let k = (mode+1) as f64*core::f64::consts::PI/length_m;
            let q2 = state.displacement_m_sqrt_kg.powi(2);
            modal_j = finite(modal_j+0.5/mu*(modulus*layer.second_moment_m4*k.powi(4)
                +modulus*layer.area_m2*e*k*k)*q2)?;
            modal_entropy = finite(modal_entropy-0.5/mu*(b*layer.second_moment_m4*k.powi(4)
                +layer.area_m2*(b*e-modulus*layer.alpha_per_k)*k*k)*q2)?;
        }
        result.push(LayerEnergy {
            axial_j: finite(reference_length_m*layer.axial_free_energy_per_reference_length_j_m)?,
            modal_j,
            entropy_j_k: finite(reference_length_m*layer.entropy_per_reference_length_j_m_k+modal_entropy)?,
        });
    }
    Ok(result)
}

fn internal_energy<E>(cx: &Cx<'_>, source: &ConcentricThermoelasticState, states: &[ModalAcousticState],
    reference_length_m: f64, length_m: f64) -> Result<f64, ConcentricStringError<E>>
{
    let mut result = finite(states.iter().map(|s| 0.5*s.velocity_m_sqrt_kg_per_s.powi(2)).sum())?;
    for (layer, energy) in source.layers.iter().zip(layer_energies::<E>(cx, source, states, reference_length_m, length_m)?) {
        result = finite(result+energy.axial_j+energy.modal_j+layer.temperature_k*energy.entropy_j_k)?;
    }
    Ok(result)
}

fn modal_gradient<E>(cx: &Cx<'_>, states: &[ModalAcousticState], length: f64) -> Result<f64, ConcentricStringError<E>> {
    let mut result = 0.0;
    for (mode, state) in states.iter().enumerate() {
        if mode % 64 == 0 { checkpoint(cx)?; }
        let kq = finite((mode+1) as f64*core::f64::consts::PI/length*state.displacement_m_sqrt_kg)?;
        result = finite(result+kq*kq)?;
    }
    Ok(result)
}

fn layer_losses<E>(section: &ConcentricKelvinVoigtSection, loss: &StringDissipation)
    -> Result<Vec<f64>, ConcentricStringError<E>>
{
    let total = section.viscous_stiffness_n_m2_s();
    let values: Vec<_> = section.layer_viscous_stiffness_n_m2_s().iter()
        .map(|&eta_i| if total == 0.0 { 0.0 } else { loss.material_heat_j*(eta_i/total) }).collect();
    let sum: f64 = values.iter().sum();
    if !sum.is_finite() || (sum-loss.material_heat_j).abs()
        > 32.0*f64::EPSILON*loss.material_heat_j {
        return Err(invalid("source-derived layer bending losses do not sum to the native solid heat"));
    }
    Ok(values)
}

// Shared by the prescribed-medium and finite-gas paths. Electrical and
// mechanical consumers receive one already resolved physical material state.
#[allow(clippy::too_many_arguments)]
fn add_ohmic_heat<'t, S, E>(
    cx: &Cx<'_>, endpoint: &DepositionSimulation<'t>, source: &ConcentricThermoelasticState,
    length_m: f64, dt: f64, cells: &[CellId], drive: OhmicDrive, tolerance_j: f64,
    coefficients: StepCoefficients<S>,
) -> Result<StepCoefficients<ConcentricStringOhmicSource<S>>, ConcentricStringError<E>> {
    checkpoint(cx)?;
    let electrical = source.section.section().thermal_transport()?
        .with_ohmic_transport(length_m)?.evaluate(drive, 0.0)?;
    let StepCoefficients { mut sources_w, conductances_w_k, payload } = coefficients;
    let check = || cx.checkpoint().map_err(|_| DepositionError::Cancelled);
    let other_external_heat_j = finite(dt*endpoint.external_step_power_w(
        &sources_w, &conductances_w_k, check)?)?;
    let mut layer_joule_heat_j = Vec::with_capacity(cells.len());
    let joule_heat_j = finite(dt*electrical.joule_power_w)?;
    for (&cell, layer) in cells.iter().zip(&electrical.layers) {
        checkpoint(cx)?;
        let heat = finite(dt*layer.joule_power_w)?;
        deposit(&mut sources_w, cell, heat, dt, tolerance_j)?;
        layer_joule_heat_j.push(heat);
    }
    let actual_heat_j = finite(dt*endpoint.external_step_power_w(
        &sources_w, &conductances_w_k, check)?)?;
    let split_residual = finite(actual_heat_j-other_external_heat_j-joule_heat_j)?;
    let roundoff = 32.0*f64::EPSILON*(actual_heat_j.abs()+other_external_heat_j.abs()+joule_heat_j);
    if split_residual.abs() > tolerance_j+roundoff {
        return Err(invalid("layer Joule heat does not match actual electrical supply into the thermal network"));
    }
    Ok(StepCoefficients { sources_w, conductances_w_k,
        payload: ConcentricStringOhmicSource { drive, electrical, joule_heat_j,
            layer_joule_heat_j, other_external_heat_j, transport: payload } })
}

fn deposit<E>(sources: &mut BTreeMap<CellId, f64>, cell: CellId, heat_j: f64, dt: f64,
    tolerance_j: f64) -> Result<(), ConcentricStringError<E>>
{
    let power = finite(heat_j/dt)?;
    let old = sources.get(&cell).copied().unwrap_or(0.0);
    let next = finite(old+power)?;
    let actual = finite((next-old)*dt)?;
    if heat_j != 0.0 && (actual == 0.0 || actual.signum() != heat_j.signum())
        || (actual-heat_j).abs() > tolerance_j {
        return Err(invalid("layer heat cannot be represented in its actual thermal source"));
    }
    sources.insert(cell, next);
    Ok(())
}

fn numerical<E>(e: FiniteDifferenceError<ConcentricStringError<E>>) -> ConcentricStringError<E> {
    match e {
        FiniteDifferenceError::Source(error) => error,
        FiniteDifferenceError::Invalid(message) => invalid(message),
    }
}
fn checkpoint<E>(cx: &Cx<'_>) -> Result<(), ConcentricStringError<E>> {
    cx.checkpoint().map_err(|_| ConcentricStringError::Acoustic(AcousticRealizeError::Cancelled))
}
fn invalid<E>(message: &'static str) -> ConcentricStringError<E> { ConcentricStringError::Invalid(message) }
fn finite<E>(value: f64) -> Result<f64, ConcentricStringError<E>> {
    if value.is_finite() { Ok(value) } else { Err(invalid("unrepresentable layered thermomechanical arithmetic")) }
}
fn positive(value: f64) -> bool { value.is_finite() && value > 0.0 }
fn same_inventory(a: f64, b: f64) -> bool {
    positive(a) && positive(b) && (a-b).abs() <= 64.0*f64::EPSILON*a.max(b)
}

#[cfg(test)]
mod tests;
