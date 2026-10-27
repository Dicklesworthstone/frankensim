//! One finite, fixed-volume source gas exchanges heat with spatial solid cells.
//!
//! The existing thermal backward-Euler owner supplies every local storage/rate
//! residual. The existing gas owner supplies NASA9 internal energy and EOS; the
//! existing material-layer owner supplies each source-bound conductance. One
//! Newton--Krylov solve closes all cell and gas joule balances at the same
//! endpoint. No thermal body is collapsed into a mean temperature.
//!
//! The gas is well mixed and its volume/composition stay fixed. Thermal geometry
//! stays fixed. Internal links may retain authored conductances or resolve the
//! original prepared material network at every proposed endpoint.
//! The returned solid field is a sequential quasistatic thermal load: reciprocal
//! thermoelastic heat, gas pressure loads, gas flow and moving surfaces require
//! their own closures. The caloric account does not include elastic energy.

use super::*;
pub use crate::cavity_phs::piston::thermal::boundary::network::material::GasConductionGeometry;
use crate::cavity_phs::piston::thermal::layers::MaterialConductionLayers;
use fs_matdb::PropertyKey;
use fs_material::thermochemical_gas::energy::{
    GasEnergyError, GasEnergySolve, GasHeatStep, IsochoricGasParcel,
};
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress,
};
use fs_thermal::enthalpy::backward_euler::{BackwardEulerStep, StepCoefficients};
use std::cell::Cell;

mod transport;
/// Declared chemistry, finite gas heat and the same spatial solid state.
pub mod reactive;
use crate::circuit_enthalpy::material::{
    MaterialNetworkError, MaterialNetworkSample, PreparedMaterialNetwork,
};

/// A distinct exposed patch on an existing structural cell. Several disjoint
/// patches may name one cell; their combined area must fit its exposed area.
#[derive(Clone, Debug)]
pub struct SolidGasContact {
    pub id: u64,
    pub cell: CellId,
    pub geometry: GasConductionGeometry,
    pub conductivity: PropertyKey,
    /// Select at the original stress-free reference; keep that claim in later intervals.
    pub selection: MaterialPropertySelection,
}

/// Numerical work and independent physical energy admission for one interval.
#[derive(Clone, Copy, Debug)]
pub struct SolidGasSolve {
    pub maximum_cells: usize,
    pub maximum_contacts: usize,
    /// Complete cell/gas/source evaluations, including Jacobian probes and final admission.
    pub maximum_trials: usize,
    pub maximum_newton_iterations: usize,
    /// Maximum-component numerical Jacobian perturbation in joules.
    pub difference_step_j: f64,
    pub cell_energy_tolerance_j: f64,
    pub thermal_energy_tolerance_j: f64,
    pub gas_energy_tolerance_j: f64,
    pub total_energy_tolerance_j: f64,
    pub maximum_gas_temperature_increment_k: f64,
    pub gas: GasEnergySolve,
    pub newton: NewtonKrylovConfig,
}

#[derive(Debug)]
pub enum SolidGasError {
    Invalid(&'static str),
    Solid(SolidEnthalpyError),
    Thermal(DepositionError),
    Gas(GasEnergyError),
    Reaction(crate::cavity_phs::piston::reactive::GasReactionError),
    Source(ThermalMaterialError),
    Transport(MaterialNetworkError),
    Newton(NewtonError),
    NotConverged(NewtonReport),
    Interrupted(SolverRunProgress),
    TrialBudget,
    Cancelled,
}
impl std::fmt::Display for SolidGasError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "spatial solid/finite gas coupling refused: {self:?}")
    }
}
impl std::error::Error for SolidGasError {}
impl From<DepositionError> for SolidGasError {
    fn from(error: DepositionError) -> Self {
        if error == DepositionError::Cancelled {
            Self::Cancelled
        } else {
            Self::Thermal(error)
        }
    }
}
impl From<SolidEnthalpyError> for SolidGasError {
    fn from(error: SolidEnthalpyError) -> Self {
        match error {
            SolidEnthalpyError::Cancelled => Self::Cancelled,
            other => Self::Solid(other),
        }
    }
}
impl From<GasEnergyError> for SolidGasError {
    fn from(error: GasEnergyError) -> Self {
        match error {
            GasEnergyError::Cancelled => Self::Cancelled,
            other => Self::Gas(other),
        }
    }
}
impl From<crate::cavity_phs::piston::reactive::GasReactionError> for SolidGasError {
    fn from(error: crate::cavity_phs::piston::reactive::GasReactionError) -> Self {
        match error {
            crate::cavity_phs::piston::reactive::GasReactionError::Cancelled => Self::Cancelled,
            other => Self::Reaction(other),
        }
    }
}
impl From<ThermalMaterialError> for SolidGasError {
    fn from(error: ThermalMaterialError) -> Self {
        Self::Source(error)
    }
}
impl From<NewtonError> for SolidGasError {
    fn from(error: NewtonError) -> Self {
        Self::Newton(error)
    }
}
impl From<MaterialNetworkError> for SolidGasError {
    fn from(error: MaterialNetworkError) -> Self {
        match error {
            MaterialNetworkError::Thermal(DepositionError::Cancelled) => Self::Cancelled,
            other => Self::Transport(other),
        }
    }
}

/// Exact endpoint material snapshot and reciprocal transfer of one patch.
#[derive(Clone, Debug)]
pub struct SolidGasContactSample {
    pub id: u64,
    pub cell: CellId,
    pub solid_temperature_k: f64,
    pub conductance_w_k: f64,
    /// Positive from this cell into the shared gas; negative heats this cell.
    pub solid_to_gas_heat_j: f64,
    pub source: ResolvedMaterialStatePoint,
}

/// Both original caloric accounts, with internal interface heats excluded from supply.
#[derive(Clone, Debug)]
pub struct SolidGasReport {
    pub thermal: EnergyAudit,
    pub gas: GasHeatStep,
    pub contacts: Vec<SolidGasContactSample>,
    /// Exact accepted cell and internal-link transport sources, when supplied.
    pub transport: Option<MaterialNetworkSample>,
    pub external_thermal_heat_j: f64,
    pub external_gas_heat_j: f64,
    pub solid_to_gas_heat_j: f64,
    pub gas_energy_residual_j: f64,
    pub combined_energy_residual_j: f64,
    pub maximum_cell_residual_j: f64,
    pub nonlinear: NewtonReport,
    pub trials: usize,
}

/// Immutable proposal. Reduce/solve the returned solid field BEFORE publishing
/// any endpoint when a downstream aperture, stability or support gate is needed.
pub struct SolidGasCandidate<'m, 't> {
    thermal: DepositionSimulation<'t>,
    gas: IsochoricGasParcel,
    state: SolidEnthalpyState<'m>,
    report: SolidGasReport,
}
impl<'m, 't> SolidGasCandidate<'m, 't> {
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub const fn gas(&self) -> &IsochoricGasParcel {
        &self.gas
    }
    pub const fn state(&self) -> &SolidEnthalpyState<'m> {
        &self.state
    }
    pub const fn report(&self) -> &SolidGasReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        IsochoricGasParcel,
        SolidEnthalpyState<'m>,
        SolidGasReport,
    ) {
        (self.thermal, self.gas, self.state, self.report)
    }
}

fn finite(value: f64) -> Result<f64, SolidGasError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(SolidGasError::Invalid("unrepresentable coupled energy"))
    }
}
fn numerical(error: FiniteDifferenceError<SolidGasError>) -> SolidGasError {
    match error {
        FiniteDifferenceError::Source(error) => error,
        FiniteDifferenceError::Invalid(message) => SolidGasError::Invalid(message),
    }
}

struct PreparedContact<'a> {
    contact: &'a SolidGasContact,
    element: usize,
    selection: MaterialPropertySelection,
}
impl PreparedContact<'_> {
    fn layers<'a>(&'a self, source: &'a TetCellMaterial<'_>) -> MaterialConductionLayers<'a> {
        MaterialConductionLayers {
            material: source.point.material,
            context: source.point.context,
            temperature_axis: source.point.temperature_axis,
            liquid_fraction_axis: source.point.liquid_fraction_axis,
            conductivity: &self.contact.conductivity,
            selection: self.selection.clone(),
            area_m2: self.contact.geometry.area_m2,
            gas_distance_m: self.contact.geometry.gas_distance_m,
            body_distance_m: self.contact.geometry.cell_distance_m,
        }
    }

    fn resolve(
        &self,
        source: &TetCellMaterial<'_>,
        gas: &fs_material::thermochemical_gas::ResolvedThermochemicalGasState,
        phase: EquilibriumPhaseState,
        transport: Option<&MaterialNetworkSample>,
    ) -> Result<
        crate::cavity_phs::piston::thermal::ConductanceSample<ResolvedMaterialStatePoint>,
        SolidGasError,
    > {
        let layers = self.layers(source);
        Ok(
            match transport::cell_material(transport, self.contact.cell)? {
                Some(material) => layers.resolve_from_material(gas, phase, material.clone())?,
                None => layers.resolve(gas, phase)?,
            },
        )
    }
}

impl<'m> SolidEnthalpy<'m> {
    /// Propose one simultaneous cell/gas thermal step and its source-resolved
    /// local solid field. Unknowns are sorted CellId energy increments followed
    /// by the TOTAL gas energy increment (including external_gas_heat_j).
    ///
    /// Q_i = dt G_i(T_cell,end - T_gas,end) is removed from its own cell and
    /// added to the ONE gas. Even zero net gas heat retains all signed Q_i.
    /// Each trial starts from the accepted inventories; the original thermal
    /// owner evaluates all existing conduction and external boundaries. Every
    /// trial also admits all structural phase/stiffness/expansion sources.
    ///
    /// Source selection is pinned at the original stress-free reference, including
    /// repeated calls on this specimen. Missing properties,
    /// supported-domain exits, duplicate/overlapping surface budgets, numerical
    /// failures and cancellation return no candidate. The geometric non-overlap
    /// of patches within the area budget remains the caller's assertion.
    ///
    /// Process start events first and split at interior events. Target events
    /// remain pending. Backward Euler is first order; convergence and energy
    /// gates do not establish time accuracy or uniqueness. All numerical probes
    /// must lie inside the declared source domains, including Jacobian stencils.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_isochoric_gas<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        gas: &IsochoricGasParcel,
        target_s: f64,
        external_gas_heat_j: f64,
        contacts: &[SolidGasContact],
        controls: SolidGasSolve,
        initial_energy_increments_j: &[f64],
        cx: &Cx<'_>,
    ) -> Result<SolidGasCandidate<'m, 't>, SolidGasError> {
        self.step_gas_with_transport(
            thermal,
            gas,
            target_s,
            external_gas_heat_j,
            contacts,
            None,
            controls,
            initial_energy_increments_j,
            cx,
        )
    }

    /// The same joint gas/solid solve with temperature-dependent internal heat
    /// transport from the original, resistor-free prepared material network.
    ///
    /// Prepare once and retain that source selection across intervals. Every
    /// trial refreshes cell k(T, phase), lowers existing half-cell resistances
    /// through the original owner, and resolves the gas patches at that same
    /// endpoint. Structural cells must share the original material and physical
    /// coordinates. A gas-facing property must use the same claim as the cell's
    /// transport bundle when one is present. Resistors refuse instead of losing
    /// their electrical work; use the circuit coupling owner for that case.
    ///
    /// All fixed-geometry, source-domain, cancellation and energy admission of
    /// `step_with_isochoric_gas` remains active. The exact accepted transport
    /// snapshots are retained in `report().transport`.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_isochoric_gas_transport<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        gas: &IsochoricGasParcel,
        target_s: f64,
        external_gas_heat_j: f64,
        contacts: &[SolidGasContact],
        transport: &PreparedMaterialNetwork<'_, '_>,
        controls: SolidGasSolve,
        initial_energy_increments_j: &[f64],
        cx: &Cx<'_>,
    ) -> Result<SolidGasCandidate<'m, 't>, SolidGasError> {
        self.step_gas_with_transport(
            thermal,
            gas,
            target_s,
            external_gas_heat_j,
            contacts,
            Some(transport),
            controls,
            initial_energy_increments_j,
            cx,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_gas_with_transport<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        gas: &IsochoricGasParcel,
        target_s: f64,
        external_gas_heat_j: f64,
        contacts: &[SolidGasContact],
        transport: Option<&PreparedMaterialNetwork<'_, '_>>,
        controls: SolidGasSolve,
        initial_energy_increments_j: &[f64],
        cx: &Cx<'_>,
    ) -> Result<SolidGasCandidate<'m, 't>, SolidGasError> {
        self.step_gas_with_endpoint(
            thermal, gas, target_s, external_gas_heat_j, contacts, transport,
            controls, initial_energy_increments_j, cx,
            |heat| {
                let mut next = gas.clone();
                let report = next.add_heat(heat, controls.gas, || cx.checkpoint().is_err())?;
                Ok((next, report, ()))
            },
        ).map(|(candidate, ())| candidate)
    }

    // One original spatial heat residual and publication boundary. The gas
    // owner may resolve composition together with its total-energy increment;
    // its source/kinetic payload travels with the EXACT accepted thermal trial.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    fn step_gas_with_endpoint<'t, P>(
        &self,
        thermal: &DepositionSimulation<'t>,
        gas: &IsochoricGasParcel,
        target_s: f64,
        external_gas_heat_j: f64,
        contacts: &[SolidGasContact],
        transport: Option<&PreparedMaterialNetwork<'_, '_>>,
        controls: SolidGasSolve,
        initial_energy_increments_j: &[f64],
        cx: &Cx<'_>,
        resolve_gas_endpoint: impl Fn(f64) -> Result<(IsochoricGasParcel, GasHeatStep, P), SolidGasError>,
    ) -> Result<(SolidGasCandidate<'m, 't>, P), SolidGasError> {
        let checkpoint = || cx.checkpoint().map_err(|_| DepositionError::Cancelled);
        checkpoint()?;
        controls.newton.validate()?;
        controls.gas.validate()?;
        if controls.maximum_cells == 0
            || thermal.cells().len() > controls.maximum_cells
            || contacts.is_empty()
            || contacts.len() > controls.maximum_contacts
            || initial_energy_increments_j.len() != thermal.cells().len().saturating_add(1)
            || controls.maximum_trials == 0
            || controls.maximum_newton_iterations == 0
            || !external_gas_heat_j.is_finite()
            || !controls.maximum_gas_temperature_increment_k.is_finite()
            || controls.maximum_gas_temperature_increment_k <= 0.0
            || [
                controls.cell_energy_tolerance_j,
                controls.thermal_energy_tolerance_j,
                controls.gas_energy_tolerance_j,
                controls.total_energy_tolerance_j,
            ]
            .iter()
            .any(|x| !x.is_finite() || *x < 0.0)
        {
            return Err(SolidGasError::Invalid(
                "bounded nonempty contacts, complete energy guess and valid physical controls required",
            ));
        }
        let initial = self.resolve(thermal, cx)?;
        let (_, initial_transport) = self.resolve_gas_transport(transport, thermal, cx)?;
        let initial_gas_temperature_k = gas.state().state().temperature;
        let mut ordered: Vec<_> = contacts.iter().collect();
        ordered.sort_by_key(|contact| contact.id);
        let mut areas = BTreeMap::<CellId, f64>::new();
        let mut prepared: Vec<PreparedContact<'_>> = Vec::with_capacity(contacts.len());
        for (index, contact) in ordered.into_iter().enumerate() {
            checkpoint()?;
            if index > 0 && prepared[index - 1].contact.id == contact.id {
                return Err(SolidGasError::Invalid("unique gas patch IDs required"));
            }
            let element = self
                .cells
                .iter()
                .position(|bound| bound.source.cell == contact.cell)
                .ok_or(SolidGasError::Invalid(
                    "gas patch must name a source-bound structural cell",
                ))?;
            let cell = &thermal.cells()[&contact.cell];
            let area = areas.entry(contact.cell).or_insert(0.0);
            *area = finite(*area + contact.geometry.area_m2)?;
            if !contact.geometry.area_m2.is_finite()
                || contact.geometry.area_m2 < 0.0
                || *area > thermal.exposed_area_m2(contact.cell)?
                || cell.boundary.htc_w_m2_k != 0.0
            {
                return Err(SolidGasError::Invalid(
                    "gas patches exceed exposed area or duplicate ambient convection",
                ));
            }
            let mut entry = PreparedContact {
                contact,
                element,
                selection: contact.selection.clone(),
            };
            let bound = &self.cells[element];
            let reference_phase = bound
                .source
                .curve
                .state_at_specific_enthalpy(bound.reference_cell.specific_h_j_kg)
                .map_err(SolidEnthalpyError::from)?;
            let sample = entry
                .layers(&bound.source)
                .resolve(gas.state(), reference_phase)?;
            let claim = sample
                .source
                .property_by_key(&contact.conductivity)
                .ok_or(SolidGasError::Invalid(
                    "layer conductivity snapshot lacks its selected property",
                ))?
                .answer()
                .receipt
                .selected;
            entry.selection = MaterialPropertySelection::PinnedByProperty(vec![(
                contact.conductivity.name().to_owned(),
                claim,
            )]);
            // The accepted start must still support the same source, even if
            // a supplied initial guess happens to return to an earlier domain.
            entry.resolve(
                &bound.source,
                gas.state(),
                initial.samples()[element].phase,
                initial_transport.as_ref(),
            )?;
            prepared.push(entry);
        }
        let step = BackwardEulerStep::new(thermal, target_s, checkpoint)?;
        let dt = step.duration_s();
        let n = step.cell_ids().len();
        let trials = Cell::new(0usize);
        let evaluate = |x: &[f64]| {
            checkpoint()?;
            if x.len() != n + 1 {
                return Err(SolidGasError::Invalid(
                    "complete cell/gas energy vector required",
                ));
            }
            if trials.get() >= controls.maximum_trials {
                return Err(SolidGasError::TrialBudget);
            }
            trials.set(trials.get() + 1);
            let (next_gas, gas_step, gas_payload) = resolve_gas_endpoint(x[n])?;
            step.evaluate(
                &x[..n],
                |endpoint| {
                    let state = self.resolve(endpoint, cx)?;
                    let (conductances, transport_sample) =
                        self.resolve_gas_transport(transport, endpoint, cx)?;
                    let empty_sources = BTreeMap::new();
                    let external = finite(
                        dt * endpoint.external_step_power_w(
                            &empty_sources,
                            &conductances,
                            checkpoint,
                        )?,
                    )?;
                    let mut heats = BTreeMap::<CellId, f64>::new();
                    let mut samples = Vec::with_capacity(prepared.len());
                    let mut net = 0.0;
                    for entry in &prepared {
                        checkpoint()?;
                        let phase = state.samples()[entry.element].phase;
                        let sample = entry.resolve(
                            &self.cells[entry.element].source,
                            next_gas.state(),
                            phase,
                            transport_sample.as_ref(),
                        )?;
                        let difference =
                            phase.temperature_k() - next_gas.state().state().temperature;
                        let power = finite(sample.conductance_w_k * difference)?;
                        let heat = finite(dt * power)?;
                        if power != 0.0 && heat == 0.0 {
                            return Err(SolidGasError::Invalid("gas patch heat underflows"));
                        }
                        net = finite(net + heat)?;
                        let cell_heat = heats.entry(entry.contact.cell).or_insert(0.0);
                        *cell_heat = finite(*cell_heat + heat)?;
                        samples.push(SolidGasContactSample {
                            id: entry.contact.id,
                            cell: entry.contact.cell,
                            solid_temperature_k: phase.temperature_k(),
                            conductance_w_k: sample.conductance_w_k,
                            solid_to_gas_heat_j: heat,
                            source: sample.source,
                        });
                    }
                    let mut sources = BTreeMap::new();
                    let mut transferred = 0.0;
                    for (cell, heat) in heats {
                        let power = finite(-heat / dt)?;
                        let actual = finite(-dt * power)?;
                        let authored = endpoint.cells()[&cell].boundary.source_w;
                        let combined = finite(authored + power)?;
                        if (heat != 0.0 && power == 0.0)
                            || (power != 0.0 && combined == authored)
                            || (authored != 0.0 && combined == power)
                            || (actual - heat).abs() > controls.total_energy_tolerance_j
                        {
                            return Err(SolidGasError::Invalid(
                                "gas heat cannot be represented at the thermal source",
                            ));
                        }
                        transferred = finite(transferred + actual)?;
                        sources.insert(cell, power);
                    }
                    if (transferred - net).abs() > controls.total_energy_tolerance_j {
                        return Err(SolidGasError::Invalid(
                            "surface heat differs from applied cell transfer",
                        ));
                    }
                    let gas_residual =
                        finite(gas_step.stored_energy_change_j - external_gas_heat_j - net)?;
                    Ok::<_, SolidGasError>(StepCoefficients {
                        sources_w: sources,
                        conductances_w_k: conductances,
                        payload: (
                            next_gas,
                            gas_step,
                            state,
                            samples,
                            net,
                            external,
                            gas_residual,
                            transport_sample,
                            gas_payload,
                        ),
                    })
                },
                checkpoint,
            )
        };
        let problem =
            FallibleFiniteDifference::new(n + 1, controls.difference_step_j, |x: &[f64]| {
                let trial = evaluate(x)?;
                let mut residual = trial.residual_j().to_vec();
                residual.push(trial.payload().6);
                Ok::<_, SolidGasError>(residual)
            })
            .map_err(numerical)?;
        let construction = NewtonKrylovState::new(
            &problem,
            initial_energy_increments_j.to_vec(),
            controls.newton,
        );
        if let Some(error) = problem.take_failure() {
            return Err(numerical(error));
        }
        let mut solver = construction?;
        let result =
            solver.run_controlled(&problem, controls.maximum_newton_iterations, &mut || {
                cx.checkpoint().is_err() || problem.has_failed()
            });
        if let Some(error) = problem.take_failure() {
            return Err(numerical(error));
        }
        checkpoint()?;
        if result.progress != SolverRunProgress::Complete {
            return Err(SolidGasError::Interrupted(result.progress));
        }
        if !result.report.converged {
            return Err(SolidGasError::NotConverged(result.report));
        }
        let trial = evaluate(&solver.x)?;
        let maximum_cell_residual_j = trial.maximum_cell_residual_j();
        let (
            thermal,
            (gas, gas_step, state, contacts, net, external, gas_energy_residual_j, transport, gas_payload),
            audit,
        ) = trial.accept(
            controls.cell_energy_tolerance_j,
            controls.thermal_energy_tolerance_j,
            checkpoint,
        )?;
        let total = finite(
            audit.stored_energy_change_j + gas_step.stored_energy_change_j
                - external
                - external_gas_heat_j,
        )?;
        if gas_energy_residual_j.abs() > controls.gas_energy_tolerance_j
            || total.abs() > controls.total_energy_tolerance_j
            || (audit.external_heat_j + net - external).abs() > controls.total_energy_tolerance_j
            || (gas.state().state().temperature - initial_gas_temperature_k).abs()
                > controls.maximum_gas_temperature_increment_k
        {
            return Err(SolidGasError::Invalid(
                "local gas, thermal transfer or combined caloric energy exceeds its allowance",
            ));
        }
        checkpoint()?;
        Ok((SolidGasCandidate {
            thermal,
            gas,
            state,
            report: SolidGasReport {
                thermal: audit,
                gas: gas_step,
                contacts,
                transport,
                external_thermal_heat_j: external,
                external_gas_heat_j,
                solid_to_gas_heat_j: net,
                gas_energy_residual_j,
                combined_energy_residual_j: total,
                maximum_cell_residual_j,
                nonlinear: result.report,
                trials: trials.get(),
            },
        }, gas_payload))
    }
}
