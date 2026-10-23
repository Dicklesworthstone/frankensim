//! Closed retained-film heating on the original modal/contact transaction.
//!
//! Only the layer's plastic loss becomes a source to its SAME fs-thermal cell.
//! Conductive links, other sources and boundary losses keep the existing finite
//! enthalpy owner. At each nonlinear endpoint, the original material claims are
//! evaluated at that cell's temperature, and the original simultaneous contact
//! solve starts from accepted elastic/plastic history. No trial advances time.
//!
//! This bounded continuation requires temperature-independent elastic modulus
//! and reference density, unchanged geometry, retained compression and solid
//! caloric state. Only yield and plastic viscosity can change with temperature;
//! changing elastic storage needs a thermoelastic law and is refused. The grid
//! is stationary and callers must justify its uniform cells and conduction.
//! Moving footprints, layer advection, expansion, melting and flash-temperature
//! approximations are not supplied. Endpoint thermal feedback is first order;
//! solver convergence is not temporal/spatial convergence or a real-rosin claim.

use super::*;
use fs_matdb::MaterialCard;
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, PhaseStateError};
use fs_material::state_point::MaterialStatePointError;
use fs_qty::{Density, QuantitySpec};
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress,
};
use fs_thermal::enthalpy::backward_euler::{BackwardEulerStep, StepCoefficients};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation, EnergyAudit};
use std::cell::{Cell, RefCell};
use std::collections::BTreeSet;
use std::convert::Infallible;

/// An existing finite layer cell, with its original material/caloric owners.
/// Every thermally continued layer has one distinct cell; no receiving volume
/// or mass is created here. The layer's resolved bundle must include density.
#[derive(Clone, Copy, Debug)]
pub struct LayerThermalCell<'a> {
    pub port: usize,
    pub cell: CellId,
    pub material: &'a MaterialCard,
    pub caloric: &'a EquilibriumEnthalpyPhaseCurve,
    pub temperature_axis: &'a str,
}

/// Independent nonlinear work, cell balance and combined energy admission.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LayerThermalSolve {
    pub maximum_cells: usize,
    pub maximum_mechanical_trials: usize,
    pub maximum_newton_iterations: usize,
    /// Finite-difference perturbation in cell ENERGY [J], not temperature.
    pub difference_step_j: f64,
    pub cell_energy_tolerance_j: f64,
    pub thermal_energy_tolerance_j: f64,
    pub total_energy_tolerance_j: f64,
    pub newton: NewtonKrylovConfig,
}

/// The accepted original pressure/contact observation and finite heat account.
pub struct LayerThermalFrame<S> {
    pub mechanics: MultiFrictionFrame,
    /// Includes the internal transfer of layer plastic heat.
    pub thermal: EnergyAudit,
    pub source: S,
    pub nonlinear: NewtonReport,
    pub mechanical_trials: usize,
    pub maximum_cell_residual_j: f64,
    /// Actual increment of the original cell sources after floating addition.
    pub deposited_heat_j: f64,
    /// Original physical plastic loss removed from the combined loss account.
    pub routed_plastic_heat_j: f64,
    /// Independently evaluated other heat, excluding the internal layer source.
    pub external_heat_j: f64,
    /// Mechanical + caloric storage + unrouted physical loss - external supplies.
    pub total_energy_residual_j: f64,
}

#[derive(Debug)]
pub enum LayerThermalError<E = Infallible> {
    Invalid(&'static str),
    Material(MaterialStatePointError),
    Phase(PhaseStateError),
    Mechanics(ModalFrictionError),
    Thermal(DepositionError),
    Source(E),
    Newton(NewtonError),
    NotConverged(NewtonReport),
    Interrupted(SolverRunProgress),
    MechanicalBudget,
    Cancelled,
}
impl<E> From<MaterialStatePointError> for LayerThermalError<E> {
    fn from(error: MaterialStatePointError) -> Self {
        Self::Material(error)
    }
}
impl<E> From<PhaseStateError> for LayerThermalError<E> {
    fn from(error: PhaseStateError) -> Self {
        Self::Phase(error)
    }
}
impl<E> From<ModalFrictionError> for LayerThermalError<E> {
    fn from(error: ModalFrictionError) -> Self {
        Self::Mechanics(error)
    }
}
impl<E> From<DepositionError> for LayerThermalError<E> {
    fn from(error: DepositionError) -> Self {
        if error == DepositionError::Cancelled {
            Self::Cancelled
        } else {
            Self::Thermal(error)
        }
    }
}
impl<E> From<NewtonError> for LayerThermalError<E> {
    fn from(error: NewtonError) -> Self {
        Self::Newton(error)
    }
}
impl<E: core::fmt::Display> core::fmt::Display for LayerThermalError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(s) => write!(f, "layer thermal admission: {s}"),
            Self::Material(e) => write!(f, "layer thermal material: {e}"),
            Self::Phase(e) => write!(f, "layer caloric state: {e}"),
            Self::Mechanics(e) => write!(f, "layer thermal mechanics: {e}"),
            Self::Thermal(e) => write!(f, "layer thermal network: {e}"),
            Self::Source(e) => write!(f, "layer thermal transport source: {e}"),
            Self::Newton(e) => write!(f, "layer thermal numerical solve: {e:?}"),
            Self::NotConverged(r) => write!(f, "layer thermal solve did not converge: {r:?}"),
            Self::Interrupted(p) => write!(f, "layer thermal solve interrupted: {p:?}"),
            Self::MechanicalBudget => write!(f, "layer thermal mechanical trial budget exhausted"),
            Self::Cancelled => write!(f, "layer thermal step cancelled"),
        }
    }
}
impl<E: std::error::Error + 'static> std::error::Error for LayerThermalError<E> {}

fn number<E>(value: f64) -> Result<f64, LayerThermalError<E>> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(LayerThermalError::Invalid(
            "unrepresentable layer thermal coordinate, power or energy",
        ))
    }
}
fn numerical<E>(error: FiniteDifferenceError<LayerThermalError<E>>) -> LayerThermalError<E> {
    match error {
        FiniteDifferenceError::Source(error) => error,
        FiniteDifferenceError::Invalid(s) => LayerThermalError::Invalid(s),
    }
}
fn same_size(actual: f64, expected: f64) -> bool {
    actual.is_finite()
        && expected.is_finite()
        && expected > 0.0
        && (actual - expected).abs() <= 32.0 * f64::EPSILON * expected
}

impl LayerThermalCell<'_> {
    fn at<E>(
        &self,
        original: &MaterialShearLayer,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<MaterialShearLayer, LayerThermalError<E>> {
        let invalid = LayerThermalError::Invalid;
        let cell = thermal
            .cells()
            .get(&self.cell)
            .ok_or(invalid("layer thermal cell is inactive"))?;
        if !thermal.uses_material_law(cell.material_id, self.caloric) {
            return Err(invalid(
                "layer must retain the original finite-cell caloric object",
            ));
        }
        let phase = self
            .caloric
            .state_at_specific_enthalpy(cell.specific_h_j_kg)?;
        if !phase.is_fully_solid() {
            return Err(invalid(
                "retained elastic layer does not admit melting or a liquid fraction",
            ));
        }
        let source = &original.source;
        if source.resolved().card_identity() != self.material.content_hash()
            || phase.material_card_identity() != self.material.content_hash()
        {
            return Err(invalid(
                "layer, caloric chart and thermal cell must share one material card",
            ));
        }
        for name in ["shear_modulus", "density"] {
            let property = source.resolved().property(name).ok_or(invalid(
                "thermal layer needs its original modulus and density claims",
            ))?;
            let claim = self
                .material
                .claims()
                .claim(property.answer().receipt.selected)
                .ok_or(invalid(
                    "layer thermal source claim is absent from its original card",
                ))?;
            if claim.value.axes().next().is_some() {
                return Err(invalid(
                    "thermal layer continuation requires constant elastic modulus and reference density",
                ));
            }
        }
        let density = source
            .resolved()
            .property("density")
            .expect("required above");
        if density.requirement().quantity() != QuantitySpec::dimensional(Density::DIMS)
            || !density.value_si().is_finite()
            || density.value_si() <= 0.0
        {
            return Err(invalid(
                "layer reference density must be a positive static mass density",
            ));
        }
        if !same_size(phase.bulk_density_kg_m3(), density.value_si())
            || self
                .caloric
                .knots()
                .iter()
                .any(|knot| !same_size(knot.bulk_density_kg_m3, density.value_si()))
        {
            return Err(invalid(
                "fixed-volume layer requires its caloric chart to retain the same source density",
            ));
        }
        let volume = number(source.area_m2() * source.thickness_m())?;
        let mass = number(volume * density.value_si())?;
        if !same_size(cell.geometry.volume_m3, volume) || !same_size(cell.mass_kg, mass) {
            return Err(invalid(
                "layer thermal cell must contain exactly its source density times area times thickness",
            ));
        }
        let query = source.resolved().query();
        if !query.axes().contains_key(self.temperature_axis) {
            return Err(invalid(
                "layer temperature must update an existing source coordinate",
            ));
        }
        let query = query
            .with(self.temperature_axis, phase.temperature_k())
            .map_err(|_| invalid("layer temperature coordinate refused its original schema"))?;
        let resolved = source.resolved().at_query_point(self.material, &query)?;
        let refreshed = ResolvedShearLayer::new(
            resolved,
            source.area_m2(),
            source.thickness_m(),
            source.maximum_elastic_shear_strain(),
        )?;
        if refreshed.temperature_kelvin().to_bits() != phase.temperature_k().to_bits()
            || refreshed.stiffness_n_m().to_bits() != source.stiffness_n_m().to_bits()
            || refreshed
                .resolved()
                .property("density")
                .expect("refreshed same bundle")
                .value_si()
                .to_bits()
                != density.value_si().to_bits()
        {
            return Err(invalid(
                "thermal layer may update only temperature and dissipative response",
            ));
        }
        let mut next = MaterialShearLayer::new(refreshed)?;
        next.state = next
            .law
            .continue_state(&original.state)
            .map_err(ModalFrictionError::from)?;
        Ok(next)
    }
}

impl MultiFrictionModalSystem {
    /// Advance the SAME modal/contact state and finite layer network together.
    ///
    /// Every layer must have a distinct explicit receiving cell. Its complete
    /// source bundle must additionally contain positive static `density`;
    /// temperature-dependent modulus/density, material replacement and phase
    /// change refuse. The selected claims and geometry remain pinned. All
    /// thermal cells, links, source models and heat capacity retain their owners.
    ///
    /// `resolve` supplies only OTHER heat and existing-link conductances at the
    /// proposed endpoint; it must not include layer plastic loss or advance
    /// external history. Its payload is returned from the accepted evaluation.
    /// Energy guesses use the original thermal owner's sorted CellId order.
    /// Source/topology events inside the sample refuse; endpoint events remain
    /// pending for explicit processing after this transaction.
    ///
    /// Actual solved normal reactions, series hair storage/loss, material
    /// memory, body/cavity reaction and the pressure observer use the original
    /// simultaneous solve. The exposed ordinary dry/isothermal step is unchanged.
    /// Dynamic footprint motion and separate nonlinear geometric storage are
    /// outside this fixed-layer adapter. A refusal leaves accepted mechanics,
    /// all source/history/heat, last pressure and both clocks unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_layer_thermal_network<'t, S, E>(
        &mut self,
        thermal: &mut DepositionSimulation<'t>,
        bindings: &[LayerThermalCell<'_>],
        controls: LayerThermalSolve,
        energy_guess_j: &[f64],
        external: &[f64],
        drives: &[FrictionDrive],
        gate: &CancelGate,
        resolve: impl Fn(&DepositionSimulation<'t>) -> Result<StepCoefficients<S>, E>,
    ) -> Result<LayerThermalFrame<S>, LayerThermalError<E>> {
        let checkpoint = || {
            if gate.is_requested() {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        };
        checkpoint()?;
        controls.newton.validate()?;
        let current_s = number(self.samples_rendered() as f64 * self.sample_period_s())?;
        let next_sample =
            self.samples_rendered()
                .checked_add(1)
                .ok_or(LayerThermalError::Invalid(
                    "layer thermal sample clock overflow",
                ))?;
        let target_s = number(next_sample as f64 * self.sample_period_s())?;
        if thermal.time_s() != current_s
            || controls.maximum_cells == 0
            || thermal.cells().len() > controls.maximum_cells
            || thermal.cells().len() != energy_guess_j.len()
            || controls.maximum_mechanical_trials == 0
            || controls.maximum_newton_iterations == 0
            || [
                controls.cell_energy_tolerance_j,
                controls.thermal_energy_tolerance_j,
                controls.total_energy_tolerance_j,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v < 0.0)
            || self.core.has_nonlinear_potentials()
        {
            return Err(LayerThermalError::Invalid(
                "incompatible layer clocks, thermal unknowns, work controls or nonlinear storage",
            ));
        }
        let original_layers = self.shear_layers.clone();
        if bindings.is_empty() || bindings.len() != original_layers.iter().flatten().count() {
            return Err(LayerThermalError::Invalid(
                "bind every retained layer to one original finite thermal cell",
            ));
        }
        let (mut ports, mut cells) = (BTreeSet::new(), BTreeSet::new());
        for binding in bindings {
            checkpoint()?;
            if !ports.insert(binding.port) || !cells.insert(binding.cell) {
                return Err(LayerThermalError::Invalid(
                    "layer thermal ports and receiving cells must be distinct",
                ));
            }
            let original = original_layers
                .get(binding.port)
                .and_then(Option::as_ref)
                .ok_or(LayerThermalError::Invalid(
                    "thermal binding requires an existing material layer",
                ))?;
            let admitted = binding.at(original, thermal)?;
            if admitted.source != original.source {
                return Err(LayerThermalError::Invalid(
                    "initial layer source must match its accepted thermal cell",
                ));
            }
        }
        self.prepare_step(drives, None, None, Some(gate))?;
        let step = BackwardEulerStep::new(thermal, target_s, checkpoint)?;
        let dt = step.duration_s();
        let trials = Cell::new(0usize);
        let mechanics = RefCell::new(self);
        let evaluate = |x: &[f64]| {
            step.evaluate(
                x,
                |endpoint| {
                    checkpoint()?;
                    if trials.get() >= controls.maximum_mechanical_trials {
                        return Err(LayerThermalError::MechanicalBudget);
                    }
                    trials.set(trials.get() + 1);
                    // Every thermal/Jacobian/contact trial restarts ACCEPTED history.
                    let mut layers = original_layers.clone();
                    for binding in bindings {
                        layers[binding.port] = Some(binding.at(
                            original_layers[binding.port].as_ref().expect("admitted"),
                            endpoint,
                        )?);
                    }
                    let mut system = mechanics.borrow_mut();
                    let original = std::mem::replace(&mut system.shear_layers, layers);
                    let staged = system.stage_multiple(external, drives, None, Some(gate));
                    let candidate_layers = std::mem::replace(&mut system.shear_layers, original);
                    // Accepted source/strain was restored even if mechanics refused.
                    let frame = staged?;
                    let physical_heat = bindings
                        .iter()
                        .map(|binding| {
                            system.staged_ports[binding.port]
                                .shear_layer
                                .expect("staged material layer")
                                .dissipated_j
                        })
                        .collect::<Vec<_>>();
                    drop(system);
                    let StepCoefficients {
                        mut sources_w,
                        conductances_w_k,
                        payload,
                    } = resolve(endpoint).map_err(LayerThermalError::Source)?;
                    let external_heat_j = number(
                        dt * endpoint.external_step_power_w(
                            &sources_w,
                            &conductances_w_k,
                            checkpoint,
                        )?,
                    )?;
                    let (mut deposited_heat_j, mut plastic_heat_j) = (0.0, 0.0);
                    for (binding, heat_j) in bindings.iter().zip(physical_heat) {
                        checkpoint()?;
                        let power = number(heat_j / dt)?;
                        if heat_j < 0.0 || heat_j != 0.0 && power == 0.0 {
                            return Err(LayerThermalError::Invalid(
                                "layer physical heat cannot vanish in power conversion",
                            ));
                        }
                        plastic_heat_j = number(plastic_heat_j + heat_j)?;
                        let other = sources_w.entry(binding.cell).or_insert(0.0);
                        let merged = number(*other + power)?;
                        let increment = number(merged - *other)?;
                        if heat_j != 0.0 && increment <= 0.0 {
                            return Err(LayerThermalError::Invalid(
                                "layer physical heat vanished when added to another cell source",
                            ));
                        }
                        let delivered = number(dt * increment)?;
                        if number(delivered - heat_j)?.abs() > controls.total_energy_tolerance_j {
                            return Err(LayerThermalError::Invalid(
                                "actual layer source increment does not preserve its physical heat",
                            ));
                        }
                        deposited_heat_j = number(deposited_heat_j + delivered)?;
                        *other = merged;
                    }
                    if number(deposited_heat_j - plastic_heat_j)?.abs()
                        > controls.total_energy_tolerance_j
                    {
                        return Err(LayerThermalError::Invalid(
                            "layer plastic loss and finite-cell heat do not balance",
                        ));
                    }
                    Ok::<_, LayerThermalError<E>>(StepCoefficients {
                        sources_w,
                        conductances_w_k,
                        payload: (
                            frame,
                            candidate_layers,
                            payload,
                            deposited_heat_j,
                            plastic_heat_j,
                            external_heat_j,
                        ),
                    })
                },
                checkpoint,
            )
        };
        let problem = FallibleFiniteDifference::new(
            step.cell_ids().len(),
            controls.difference_step_j,
            |x: &[f64]| evaluate(x).map(|t| t.residual_j().to_vec()),
        )
        .map_err(numerical)?;
        let construction =
            NewtonKrylovState::new(&problem, energy_guess_j.to_vec(), controls.newton);
        if let Some(error) = problem.take_failure() {
            return Err(numerical(error));
        }
        let mut solver = construction?;
        let result =
            solver.run_controlled(&problem, controls.maximum_newton_iterations, &mut || {
                gate.is_requested() || problem.has_failed()
            });
        if let Some(error) = problem.take_failure() {
            return Err(numerical(error));
        }
        checkpoint()?;
        if result.progress != SolverRunProgress::Complete {
            return Err(LayerThermalError::Interrupted(result.progress));
        }
        if !result.report.converged {
            return Err(LayerThermalError::NotConverged(result.report));
        }
        // Last evaluation leaves the original mechanical candidate at precisely
        // this thermal endpoint. Nothing evaluates another trial after it.
        let trial = evaluate(&solver.x)?;
        let maximum_cell_residual_j = trial.maximum_cell_residual_j();
        let (
            next,
            (frame, layers, source, deposited_heat_j, plastic_heat_j, external_heat_j),
            thermal_audit,
        ) = trial.accept(
            controls.cell_energy_tolerance_j,
            controls.thermal_energy_tolerance_j,
            checkpoint,
        )?;
        let total = number(
            frame.energy_residual_j + thermal_audit.stored_energy_change_j
                - plastic_heat_j
                - external_heat_j,
        )?;
        if total.abs() > controls.total_energy_tolerance_j {
            return Err(LayerThermalError::Invalid(
                "modal, layer and caloric energy do not close",
            ));
        }
        checkpoint()?;
        let system = mechanics.into_inner();
        // All admission is complete. Publish source+strain, original mechanical
        // state/observer, and the existing thermal candidate without fallible work.
        system.shear_layers = layers;
        let mechanics = system.publish_step(frame).clone();
        *thermal = next;
        Ok(LayerThermalFrame {
            mechanics,
            thermal: thermal_audit,
            source,
            nonlinear: result.report,
            mechanical_trials: trials.get(),
            maximum_cell_residual_j,
            deposited_heat_j,
            routed_plastic_heat_j: plastic_heat_j,
            external_heat_j,
            total_energy_residual_j: total,
        })
    }
}
