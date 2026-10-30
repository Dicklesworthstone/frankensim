//! Closed retained-film heating on the original modal/contact transaction.
//!
//! Plastic loss and reversible thermoelastic exchange reach the SAME thermal cell.
//! Explicit series-drive routes can deposit physical dashpot loss into other
//! original cells; their existing conductive links close the contact feedback.
//! Conductive links, other sources and boundary losses keep the existing finite
//! enthalpy owner. At each nonlinear endpoint, the original material claims are
//! evaluated at that cell's temperature, and the original simultaneous contact
//! solve starts from accepted elastic/plastic history. No trial advances time.
//!
//! This bounded continuation admits affine G(T), fixed reference density,
//! unchanged film geometry, retained compression and solid caloric state. Mechanical
//! free energy is psi=k(T)e^2/2, entropy s=-k'e^2/2 and internal energy psi+T*s.
//! Source knots are checked across the complete caloric interval; a kink in G
//! refuses because its entropy jump is absent from this thermal model. The grid
//! retains reference mass/volume: the stress-free caloric enthalpy changes are
//! internal-energy changes here, with no evolving pressure-volume work.
//! Affine G contributes zero fixed-strain heat-capacity correction. Endpoint
//! thermal stepping can produce numerical entropy; it is not plastic heating.
//! Cells retain material coordinates; callers must justify their uniform state.
//! A trial-aware source may conduct into the actual sliding counterface footprint.
//! Independent normal supports may move through the original contact-work owner.
//! Moving mechanical film maps, layer advection, expansion, melting and flash-temperature
//! approximations are not supplied. Endpoint thermal feedback is first order;
//! solver convergence is not temporal/spatial convergence or a real-rosin claim.

/// Material-coordinate finite-cell dry friction heat and contact conduction.
pub mod dry;
/// Moving counterface heat around an original material-attached shear film.
pub mod film_surface;
mod source;

use super::*;
use dry::{DryContactHeat, DryContactHeatFrame};
use source::DryTemperatureUnknown;
use fs_material::thermal::source::ThermalMaterialError;
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
use std::collections::{BTreeMap, BTreeSet};
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

/// Deposit ONE existing series-drive dashpot's physical loss into an existing
/// finite cell. No partition, receiving mass, heat capacity or conduction is
/// inferred. Several drives may feed one declared uniform receiving body, but
/// each port may be routed only once. Unrouted losses remain in the mechanical
/// report's external loss account.
///
/// The declared spring/dashpot retains its fixed coefficients. This route does
/// not identify a hair material or infer temperature dependence from k and c;
/// the caller must justify that constitutive approximation over the receiving
/// cell's temperature range. Source-backed films still refresh their own laws.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SeriesDriveHeatDestination {
    pub port: usize,
    pub cell: CellId,
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
    /// Original held-coefficient mechanical work/pressure report. Its shear
    /// storage is free energy; the terms below complete the thermal account.
    pub mechanics: MultiFrictionFrame,
    /// Includes routed plastic/drive heat and signed reversible layer exchange.
    pub thermal: EnergyAudit,
    pub source: S,
    pub nonlinear: NewtonReport,
    pub mechanical_trials: usize,
    pub maximum_cell_residual_j: f64,
    /// Actual increment of the original cell sources after floating addition.
    pub deposited_heat_j: f64,
    /// Original physical plastic loss removed from the combined loss account.
    pub routed_plastic_heat_j: f64,
    /// Maxwell-arm relaxation heat, transferred to the SAME original film
    /// cells. It may be positive below plastic yield or at held total slip.
    pub routed_relaxation_heat_j: f64,
    /// Physical series-drive dashpot loss deposited into its named cells.
    /// Already included in the mechanical report, so removed once from the
    /// combined external loss account. Never includes stored spring energy.
    pub routed_series_drive_heat_j: f64,
    /// Original dry friction loss delivered to the same finite thermal cells.
    pub routed_dry_heat_j: f64,
    /// Physical contact temperatures and signed body fluxes, including zero-slip conduction.
    pub dry_contacts: Vec<DryContactHeatFrame>,
    /// Signed reversible caloric exchange, never classified as dissipation.
    pub routed_reversible_heat_j: f64,
    /// psi(T_new,e_old)-psi(T_old,e_old), in addition to the mechanical report.
    pub coefficient_free_energy_change_j: f64,
    /// Change of the mechanical entropy contribution T*s to internal energy.
    pub mechanical_entropy_energy_change_j: f64,
    /// Actual change of layer internal energy psi+T*s across both endpoints.
    pub layer_internal_energy_change_j: f64,
    /// Independently evaluated other heat, excluding all internal routed sources.
    pub external_heat_j: f64,
    /// Mechanical + caloric storage + unrouted physical loss - external supplies.
    pub total_energy_residual_j: f64,
}

#[derive(Debug)]
pub enum LayerThermalError<E = Infallible> {
    Invalid(&'static str),
    Material(MaterialStatePointError),
    ThermalMaterial(ThermalMaterialError),
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
            Self::ThermalMaterial(e) => write!(f, "dry contact thermal material: {e}"),
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
    ) -> Result<(MaterialShearLayer, f64), LayerThermalError<E>> {
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
        for name in ["density"] {
            let property = source.resolved().property(name).ok_or(invalid(
                "thermal layer needs its original density claim",
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
                    "thermal layer continuation requires constant reference density",
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
        let refreshed = source.at_query_point(self.material, &query)?;
        if refreshed.temperature_kelvin().to_bits() != phase.temperature_k().to_bits()
            || refreshed
                .resolved()
                .property("density")
                .expect("refreshed same bundle")
                .value_si()
                .to_bits()
                != density.value_si().to_bits()
        {
            return Err(invalid(
                "thermal layer must retain its reference density and actual temperature",
            ));
        }
        let mut next = MaterialShearLayer::new(refreshed)?;
        let continuation = next
            .law
            .continue_state_with_work(&original.state)
            .map_err(ModalFrictionError::from)?;
        next.state = continuation.state;
        Ok((next, continuation.coefficient_free_energy_change_j))
    }
}

/// Original source and retained strain at one proposed finite-cell endpoint.
/// This is a proposal carrier; no source or mechanical history is published.
pub(in crate::render::schedule::force::coupled::contact::friction) struct ContinuedMaterialLayer {
    pub source: ResolvedShearLayer,
    pub law: ThirdBodyLayerContact,
    pub state: ThirdBodyLayerState,
    pub coefficient_free_energy_change_j: f64,
}

/// Reuse the exact material/card/caloric/geometry admission for vector consumers.
pub(in crate::render::schedule::force::coupled::contact::friction) fn continue_material_layer<E>(
    source: &ResolvedShearLayer,
    state: &ThirdBodyLayerState,
    binding: LayerThermalCell<'_>,
    thermal: &DepositionSimulation<'_>,
) -> Result<ContinuedMaterialLayer, LayerThermalError<E>> {
    let mut original = MaterialShearLayer::new(source.clone())?;
    original.state = state.clone();
    if !original.energy_j().is_finite() || original.energy_j() < 0.0 {
        return Err(LayerThermalError::Invalid(
            "material layer and accepted constitutive history must have the same owner",
        ));
    }
    let (next, coefficient_free_energy_change_j) = binding.at(&original, thermal)?;
    Ok(ContinuedMaterialLayer {
        source: next.source,
        law: next.law,
        state: next.state,
        coefficient_free_energy_change_j,
    })
}

impl MultiFrictionModalSystem {
    /// Advance the SAME modal/contact state and finite layer network together.
    ///
    /// Every layer must have a distinct explicit receiving cell. Its complete
    /// source bundle must additionally contain positive static `density`;
    /// density changes, non-affine modulus, material replacement and phase
    /// change refuse. Affine G(T) retains mechanical entropy and coefficient
    /// work. The selected claims and geometry remain pinned. All
    /// thermal cells, links, source models and heat capacity retain their owners.
    ///
    /// `resolve` supplies only OTHER heat and existing-link conductances at the
    /// proposed endpoint; it must not include layer plastic/reversible heat or advance
    /// external history. Its payload is returned from the accepted evaluation.
    /// Energy guesses use the original thermal owner's sorted CellId order.
    /// Source/topology events inside the sample refuse; endpoint events remain
    /// pending for explicit processing after this transaction.
    ///
    /// Actual solved normal reactions, series hair storage/loss, material
    /// memory, body/cavity reaction and the pressure observer use the original
    /// simultaneous solve. Bound nonlinear mechanical storage uses its original
    /// discrete-gradient iteration at each thermal endpoint; its material and
    /// strain law stay frozen while the declared film/receiving cells evolve.
    /// Dynamic footprint motion remains outside this fixed-layer adapter.
    /// A refusal leaves accepted mechanics,
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
        self.step_with_contact_thermal_network(thermal, bindings, &[], controls,
            energy_guess_j, external, drives, gate, resolve)
    }

    /// Advance the original contact/pressure/thermal transaction, including
    /// explicit heat destinations for existing series-compliant drives.
    ///
    /// Each routed source is the actual accepted c*(delta eta)^2/dt, not total
    /// contact work or spring-energy release. Original thermal links conduct
    /// it to other cells, including temperature-dependent material films, in
    /// the SAME nonlinear endpoint solve. Each trial restarts accepted stretch,
    /// strain, mechanics and caloric state; rejection publishes none of them.
    ///
    /// Every material film must still have its original binding. Empty film
    /// bindings are admitted for a drive-only thermal route. Unrouted modal,
    /// dry-contact, obstacle and drive losses keep their existing accounts.
    /// `resolve` supplies only OTHER heat: it must not add either routed loss.
    /// Receiving cells already own their physical volume, mass and caloric law.
    /// This transfers fixed-coefficient drive loss, without inventing material
    /// temperature laws for those drives. See `SeriesDriveHeatDestination`.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_contact_thermal_network<'t, S, E>(
        &mut self,
        thermal: &mut DepositionSimulation<'t>,
        bindings: &[LayerThermalCell<'_>],
        series_heat: &[SeriesDriveHeatDestination],
        controls: LayerThermalSolve,
        energy_guess_j: &[f64],
        external: &[f64],
        drives: &[FrictionDrive],
        gate: &CancelGate,
        resolve: impl Fn(&DepositionSimulation<'t>) -> Result<StepCoefficients<S>, E>,
    ) -> Result<LayerThermalFrame<S>, LayerThermalError<E>> {
        self.step_with_contact_thermal_network_and_obstacles(thermal, bindings, series_heat,
            controls, energy_guess_j, external, drives, None, None, gate, resolve)
    }

    /// Continue the same contact/enthalpy transaction with prescribed normal
    /// support geometry and translation. Omitted controls hold their accepted
    /// values. The original contact owner supplies the exact product-rule
    /// geometry work and normal loss; neither is deposited as film plastic heat.
    ///
    /// Geometry is prepared once per sample. Every thermal/Jacobian trial starts
    /// from the same accepted contact maps, offsets and constitutive histories.
    /// New maps, offsets, mechanics, pressure and heat publish together only
    /// after the complete energy balance accepts. A refusal publishes none.
    ///
    /// Independent finger/fret supports can migrate. A normal bound to a
    /// tangential contact must retain its map: this entry point cannot transport
    /// the matching tangent or film history. Film area/thickness, thermal cells
    /// and solid-state restrictions remain those of the fixed-film owner.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_contact_thermal_network_and_obstacles<'t, S, E>(
        &mut self,
        thermal: &mut DepositionSimulation<'t>,
        bindings: &[LayerThermalCell<'_>],
        series_heat: &[SeriesDriveHeatDestination],
        controls: LayerThermalSolve,
        energy_guess_j: &[f64],
        external: &[f64],
        drives: &[FrictionDrive],
        obstacle_geometry: Option<&[obstacles::geometry::ObstacleGeometry]>,
        obstacle_offsets: Option<&[f64]>,
        gate: &CancelGate,
        resolve: impl Fn(&DepositionSimulation<'t>) -> Result<StepCoefficients<S>, E>,
    ) -> Result<LayerThermalFrame<S>, LayerThermalError<E>> {
        self.step_with_contact_thermal_feedback(thermal, bindings, series_heat, controls,
            energy_guess_j, external, drives, obstacle_geometry, obstacle_offsets, gate,
            |endpoint, _ports| resolve(endpoint))
    }

    /// Resolve other heat from the SAME candidate thermal state and mechanical
    /// port displacements. This closes moving-counterface heat feedback without
    /// advancing a separate stroke cursor or using the previous sample's slip.
    /// Retained films keep their original patch, mass and shear history. The
    /// callback must not redeposit plastic/drive loss already owned here, and
    /// internal transfers must balance over the original receiving cells.
    /// Failed source/coverage checks publish neither motion nor thermal history.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_contact_thermal_feedback<'t, S, E>(
        &mut self,
        thermal: &mut DepositionSimulation<'t>,
        bindings: &[LayerThermalCell<'_>],
        series_heat: &[SeriesDriveHeatDestination],
        controls: LayerThermalSolve,
        energy_guess_j: &[f64],
        external: &[f64],
        drives: &[FrictionDrive],
        obstacle_geometry: Option<&[obstacles::geometry::ObstacleGeometry]>,
        obstacle_offsets: Option<&[f64]>,
        gate: &CancelGate,
        resolve: impl Fn(&DepositionSimulation<'t>, &[FrictionPortFrame]) -> Result<StepCoefficients<S>, E>,
    ) -> Result<LayerThermalFrame<S>, LayerThermalError<E>> {
        self.step_with_contact_thermal_inner(thermal, bindings, series_heat, &[], controls,
            energy_guess_j, external, drives, obstacle_geometry, obstacle_offsets, false, gate, resolve)
    }

    #[allow(clippy::too_many_arguments)]
    fn step_with_contact_thermal_inner<'t, S, E>(
        &mut self,
        thermal: &mut DepositionSimulation<'t>,
        bindings: &[LayerThermalCell<'_>],
        series_heat: &[SeriesDriveHeatDestination],
        dry: &[DryContactHeat<'_>],
        controls: LayerThermalSolve,
        energy_guess_j: &[f64],
        external: &[f64],
        drives: &[FrictionDrive],
        obstacle_geometry: Option<&[obstacles::geometry::ObstacleGeometry]>,
        obstacle_offsets: Option<&[f64]>,
        contact_patch_motion: bool,
        gate: &CancelGate,
        resolve: impl Fn(&DepositionSimulation<'t>, &[FrictionPortFrame]) -> Result<StepCoefficients<S>, E>,
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
        {
            return Err(LayerThermalError::Invalid(
                "incompatible layer clocks, thermal unknowns or work controls",
            ));
        }
        let original_layers = self.shear_layers.clone();
        if bindings.is_empty() && series_heat.is_empty() && dry.is_empty()
            || bindings.len() != original_layers.iter().flatten().count() {
            return Err(LayerThermalError::Invalid(
                "bind every retained layer to one original finite thermal cell",
            ));
        }
        let mut heated_dry = BTreeSet::new();
        for route in dry {
            checkpoint()?;
            if !heated_dry.insert(route.port) {
                return Err(LayerThermalError::Invalid("one dry heat route per original contact is required"));
            }
            route.admit(self, thermal, controls.maximum_cells)?;
        }
        if (0..self.port_count()).any(|port|
            self.interface(port).is_some_and(|source| source.temperature_source().is_some())
                && !heated_dry.contains(&port)) {
            return Err(LayerThermalError::Invalid(
                "every temperature-dependent friction source needs its original finite dry heat route",
            ));
        }
        let mut heated_drives = BTreeSet::new();
        for route in series_heat {
            checkpoint()?;
            if !heated_drives.insert(route.port)
                || self.series_drives.get(route.port).and_then(Option::as_ref).is_none()
                || !thermal.cells().contains_key(&route.cell) {
                return Err(LayerThermalError::Invalid(
                    "each heat route needs one distinct existing series drive and an existing receiving cell",
                ));
            }
        }
        let (mut ports, mut cells) = (BTreeSet::new(), BTreeSet::new());
        let mut thermodynamics = Vec::with_capacity(bindings.len());
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
            let (admitted, _) = binding.at(original, thermal)?;
            // Refresh pins the already selected claims. Compare two pinned
            // evaluations so selection-policy receipts cannot masquerade as a
            // changed physical state; a stale temperature still fails exactly.
            let original_pinned = original.source.at_query_point(binding.material,
                &original.source.resolved().query())?;
            if admitted.source != original_pinned {
                return Err(LayerThermalError::Invalid(
                    "initial layer source must match its accepted thermal cell",
                ));
            }
            let knots = binding.caloric.knots();
            let span = (knots.first().expect("admitted caloric chart").temperature_k,
                knots.last().expect("admitted caloric chart").temperature_k);
            let slope = original.source.thermal_stiffness_slope_n_m_k(
                binding.material, binding.temperature_axis, span)?;
            let initial = original.law.thermodynamics(&original.state, slope)
                .map_err(ModalFrictionError::from)?;
            thermodynamics.push((slope, initial));
        }
        self.prepare_step(drives, obstacle_offsets, obstacle_geometry, contact_patch_motion, Some(gate))?;
        let step = BackwardEulerStep::new(thermal, target_s, checkpoint)?;
        let dt = step.duration_s();
        let temperatures = DryTemperatureUnknown::prepare(self, dry, thermal, dt)?;
        let cell_count = step.cell_ids().len();
        let unknown_count = cell_count + temperatures.len();
        let trials = Cell::new(0usize);
        let mechanics = RefCell::new(self);
        let evaluate = |x: &[f64]| {
            if x.len() != unknown_count {
                return Err(LayerThermalError::Invalid("incomplete cell and interface-temperature iterate"));
            }
            step.evaluate(
                &x[..cell_count],
                |endpoint| {
                    checkpoint()?;
                    if trials.get() >= controls.maximum_mechanical_trials {
                        return Err(LayerThermalError::MechanicalBudget);
                    }
                    trials.set(trials.get() + 1);
                    // Every thermal/Jacobian/contact trial restarts ACCEPTED history.
                    let mut layers = original_layers.clone();
                    let mut coefficient_changes = Vec::with_capacity(bindings.len());
                    for binding in bindings {
                        let (layer, coefficient_work) = binding.at(
                            original_layers[binding.port].as_ref().expect("admitted"),
                            endpoint,
                        )?;
                        layers[binding.port] = Some(layer);
                        coefficient_changes.push(coefficient_work);
                    }
                    let proposed_sources = temperatures.iter().zip(&x[cell_count..])
                        .map(|(source, &coordinate)| source.binding(coordinate))
                        .collect::<Result<Vec<_>, LayerThermalError<E>>>()?;
                    let mut system = mechanics.borrow_mut();
                    let original = std::mem::replace(&mut system.shear_layers, layers);
                    let staged = system.stage_temperature_sources(
                        &proposed_sources, external, drives, obstacle_offsets, gate,
                    );
                    let candidate_layers = std::mem::replace(&mut system.shear_layers, original);
                    // Accepted source/strain was restored even if mechanics refused.
                    let frame = staged?;
                    let mut physical_heat = Vec::with_capacity(bindings.len());
                    let (mut coefficient_work_j, mut entropy_energy_j,
                        mut internal_change_j, mut reversible_heat_j) = (0.0, 0.0, 0.0, 0.0);
                    for ((binding, &(slope, initial)), coefficient_work) in bindings.iter()
                        .zip(&thermodynamics).zip(coefficient_changes) {
                        checkpoint()?;
                        let layer = candidate_layers[binding.port].as_ref().expect("admitted");
                        let staged = layer.staged.as_ref().expect("staged material layer");
                        let current = layer.law.thermodynamics(&staged.state, slope)
                            .map_err(ModalFrictionError::from)?;
                        let old_temperature = original_layers[binding.port].as_ref()
                            .expect("admitted").source.temperature_kelvin();
                        let entropy_energy_change = if slope == 0.0 { 0.0 } else {
                            number(layer.source.temperature_kelvin() * current.entropy_j_k
                                - old_temperature * initial.entropy_j_k)?
                        };
                        // Native mechanics solves with endpoint k. Its work is
                        // delta(psi at that k)+D. The missing coefficient work
                        // and delta(T*s) are reversible exchange with caloric H.
                        let reversible = number(-coefficient_work - entropy_energy_change)?;
                        coefficient_work_j = number(coefficient_work_j + coefficient_work)?;
                        entropy_energy_j = number(entropy_energy_j + entropy_energy_change)?;
                        internal_change_j = number(internal_change_j
                            + current.internal_energy_j - initial.internal_energy_j)?;
                        reversible_heat_j = number(reversible_heat_j + reversible)?;
                        let relaxation = layer.law.step_dissipation(&layer.state,
                            &staged.state, dt).map_err(ModalFrictionError::from)?.relaxation_j;
                        physical_heat.push((number(staged.dissipated_j-relaxation)?, relaxation,
                            number(staged.dissipated_j + reversible)?));
                    }
                    let mut series_loss = Vec::with_capacity(series_heat.len());
                    let mut series_heat_j = 0.0;
                    for route in series_heat {
                        checkpoint()?;
                        let heat_j = number(system.staged_ports[route.port].series_drive_dissipation_j)?;
                        if heat_j < 0.0 {
                            return Err(LayerThermalError::Invalid("series drive physical heat must be nonnegative"));
                        }
                        series_heat_j = number(series_heat_j + heat_j)?;
                        series_loss.push((route.cell, heat_j));
                    }
                    let mut dry_frames = Vec::with_capacity(dry.len());
                    let mut dry_sources = Vec::with_capacity(dry.len());
                    let mut dry_heat_j = 0.0;
                    let mut interface_residual_j = 0.0;
                    for route in dry {
                        checkpoint()?;
                        let port = &system.staged_ports[route.port];
                        let report = if let Some((index, source)) = temperatures.iter().enumerate()
                            .find(|(_, source)| source.port == route.port) {
                            let report = route.evaluate_at_temperature(endpoint, port.normal_force_n,
                                port.dissipation_j, dt, source.temperature(x[cell_count + index])?)?;
                            interface_residual_j = number(interface_residual_j + report.interface_energy_residual_j)?;
                            report
                        } else {
                            route.evaluate(endpoint, port.normal_force_n,
                                port.dissipation_j, dt, controls.total_energy_tolerance_j)?
                        };
                        if let Some(temperature) = report.interface_temperature_k {
                            let source = system.validate_dry_temperature(route.port, temperature, drives[route.port])?;
                            dry_sources.push((route.port, source));
                        }
                        dry_heat_j = number(dry_heat_j + report.friction_heat_j)?;
                        dry_frames.push(report);
                    }
                    let StepCoefficients {
                        mut sources_w,
                        conductances_w_k,
                        payload,
                    } = resolve(endpoint, &system.staged_ports).map_err(LayerThermalError::Source)?;
                    drop(system);
                    let external_heat_j = number(
                        dt * endpoint.external_step_power_w(
                            &sources_w,
                            &conductances_w_k,
                            checkpoint,
                        )?,
                    )?;
                    let (mut deposited_heat_j, mut plastic_heat_j, mut relaxation_heat_j) = (0.0, 0.0, 0.0);
                    for (binding, (plastic_j, relaxation_j, heat_j)) in bindings.iter().zip(physical_heat) {
                        checkpoint()?;
                        if plastic_j < 0.0 {
                            return Err(LayerThermalError::Invalid(
                                "layer plastic loss must be nonnegative",
                            ));
                        }
                        plastic_heat_j = number(plastic_heat_j + plastic_j)?;
                        relaxation_heat_j = number(relaxation_heat_j + relaxation_j)?;
                        let delivered = deposit_heat(&mut sources_w, binding.cell,
                            heat_j, dt, controls.total_energy_tolerance_j)?;
                        deposited_heat_j = number(deposited_heat_j + delivered)?;
                    }
                    for (cell, heat_j) in series_loss {
                        checkpoint()?;
                        let delivered = deposit_heat(&mut sources_w, cell,
                            heat_j, dt, controls.total_energy_tolerance_j)?;
                        deposited_heat_j = number(deposited_heat_j + delivered)?;
                    }
                    for report in &dry_frames {
                        for &(cell, heat_j) in &report.cell_heat_j {
                            checkpoint()?;
                            let delivered = deposit_heat(&mut sources_w, cell,
                                heat_j, dt, controls.total_energy_tolerance_j)?;
                            deposited_heat_j = number(deposited_heat_j + delivered)?;
                        }
                    }
                    if number(deposited_heat_j - plastic_heat_j - relaxation_heat_j - reversible_heat_j - series_heat_j
                        - dry_heat_j - interface_residual_j)?.abs()
                        > controls.total_energy_tolerance_j
                    {
                        return Err(LayerThermalError::Invalid(
                            "contact plastic loss, reversible exchange, series loss and finite-cell heat do not balance",
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
                            relaxation_heat_j,
                            external_heat_j,
                            reversible_heat_j,
                            coefficient_work_j,
                            entropy_energy_j,
                            internal_change_j,
                            series_heat_j,
                            dry_heat_j, dry_frames, dry_sources,
                        ),
                    })
                },
                checkpoint,
            )
        };
        let problem = FallibleFiniteDifference::new(
            unknown_count,
            controls.difference_step_j,
            |x: &[f64]| evaluate(x).map(|t| {
                let mut residual = t.residual_j().to_vec();
                for (index, source) in temperatures.iter().enumerate() {
                    let frame = t.payload().13.iter().find(|frame| frame.port == source.port)
                        .expect("one report per admitted dry source");
                    residual.push(if frame.interface_temperature_k.is_some() {
                        frame.interface_energy_residual_j
                    } else {
                        // No physical interface exists while released. Anchor
                        // its unused numerical coordinate, without heat/storage.
                        x[cell_count + index]
                    });
                }
                residual
            }),
        )
        .map_err(numerical)?;
        let mut initial = energy_guess_j.to_vec();
        initial.resize(unknown_count, 0.0);
        let construction =
            NewtonKrylovState::new(&problem, initial, controls.newton);
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
        if trial.payload().13.iter().any(|contact|
            contact.interface_energy_residual_j.abs() > controls.total_energy_tolerance_j) {
            return Err(LayerThermalError::Invalid(
                "dry interface fluxes do not close the same temperature-dependent friction work",
            ));
        }
        let (
            next,
            (frame, layers, source, deposited_heat_j, plastic_heat_j, relaxation_heat_j, external_heat_j,
                reversible_heat_j, coefficient_work_j, entropy_energy_j, internal_change_j, series_heat_j,
                dry_heat_j, dry_contacts, dry_sources),
            thermal_audit,
        ) = trial.accept(
            controls.cell_energy_tolerance_j,
            controls.thermal_energy_tolerance_j,
            checkpoint,
        )?;
        let total = number(
            frame.energy_residual_j + thermal_audit.stored_energy_change_j
                - plastic_heat_j - relaxation_heat_j
                - series_heat_j - dry_heat_j - external_heat_j + coefficient_work_j + entropy_energy_j,
        )?;
        if total.abs() > controls.total_energy_tolerance_j {
            return Err(LayerThermalError::Invalid(
                "modal, layer and caloric energy do not close",
            ));
        }
        checkpoint()?;
        let system = mechanics.into_inner();
        system.validate_obstacle_sources(Some(gate))?;
        // All admission is complete. Publish source+strain, original mechanical
        // state/observer, and the existing thermal candidate without fallible work.
        system.shear_layers = layers;
        for (port, source) in dry_sources {
            if source.temperature_source().is_some() { system.install_dry_source(port, source); }
            else if port == 0 { system.core.interface = Some(source); }
            else { system.interfaces[port - 1] = Some(source); }
        }
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
            routed_relaxation_heat_j: relaxation_heat_j,
            routed_series_drive_heat_j: series_heat_j,
            routed_dry_heat_j: dry_heat_j,
            dry_contacts,
            routed_reversible_heat_j: reversible_heat_j,
            coefficient_free_energy_change_j: coefficient_work_j,
            mechanical_entropy_energy_change_j: entropy_energy_j,
            layer_internal_energy_change_j: internal_change_j,
            external_heat_j,
            total_energy_residual_j: total,
        })
    }
}

// One source insertion rule for signed reversible exchange and positive
// physical losses. It checks the ACTUAL representable increment in each cell,
// including multiple mechanical sources sharing one finite receiving body.
pub(in crate::render::schedule::force::coupled::contact::friction) fn deposit_heat<E>(sources_w: &mut BTreeMap<CellId, f64>, cell: CellId,
    heat_j: f64, dt: f64, tolerance_j: f64) -> Result<f64, LayerThermalError<E>> {
    let power = number(heat_j / dt)?;
    if heat_j != 0.0 && power == 0.0 {
        return Err(LayerThermalError::Invalid("contact heat vanished in power conversion"));
    }
    let other = sources_w.entry(cell).or_insert(0.0);
    let merged = number(*other + power)?;
    let increment = number(merged - *other)?;
    if heat_j != 0.0 && (increment == 0.0 || (increment > 0.0) != (heat_j > 0.0)) {
        return Err(LayerThermalError::Invalid("contact heat vanished in the receiving-cell source sum"));
    }
    let delivered = number(dt * increment)?;
    if number(delivered - heat_j)?.abs() > tolerance_j {
        return Err(LayerThermalError::Invalid("actual contact source increment does not preserve its physical heat"));
    }
    *other = merged;
    Ok(delivered)
}
