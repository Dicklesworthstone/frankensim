//! Native rod motion and complete layer thermodynamics on original finite cells.
//!
//! First advance the accepted-temperature native rod. At its unchanged
//! endpoint pose solve the original finite enthalpy residual, with each layer's
//! Q_rev = -Delta_T Psi - Delta(T S). This equals frozen-temperature strain work
//! minus the layer's complete internal-energy change, without thermalizing the
//! native rod's signed integration defect. Source E/G/alpha then refresh the
//! SAME rod, preserving pose, momentum, reference inertia, loads and clock.
//!
//! The split is first order and requires timestep refinement. Constant source
//! density/alpha, affine E/G, complete concentric layers and wholly solid caloric
//! charts are admitted. Fixed reference geometry, original small fibre strains,
//! stationary supports and separately bound unheated links retain their native
//! limits. Optional source-bound torsional viscosity deposits its native physical
//! heat into these same original layer cells. No phase transition, cross-section
//! expansion or alternate rod/thermal integrator is introduced.

/// Intrinsic torsional shear heat on the same finite-cell timeline.
pub mod torsional;
pub use torsional::{RodLayerTorsionalHeat, RodTorsionalHeating, RodTorsionalViscosity};

use super::{
    PreparedRodLoads, RodMotion, RodMotionControls, RodMotionError, RodMotionInertia,
    RodMotionStep, RodNodalWrench, StageLoads,
};
use crate::circuit_enthalpy::material::{
    MaterialNetwork, MaterialNetworkError, MaterialNetworkSample, PreparedMaterialNetwork,
};
use crate::circuit_enthalpy::{
    CircuitEnthalpyError, CircuitEnthalpySolve, DrivenCoefficients, DriverAccount, ThermalDriver,
    cell_driver_candidate_entropy_option,
};
use crate::material_circuit::MaterialCircuitSolve;
use fs_exec::CancelGate;
use fs_material::section_temperature::thermoelastic::rod::{
    ConcentricRodThermoelasticLaw, ConcentricRodThermoelasticState,
};
use fs_material::state_point::MaterialStatePointError;
use fs_solid::rod::material_field::motion::thermal::RodSegmentStrainState;
use fs_solid::rod::material_field::thermal::source::RodMaterialSource;
use fs_solid::rod::material_field::{SpatialRodError, motion::RodNodalInertia};
use fs_solver::{NewtonKrylovConfig, NewtonReport};
use fs_thermal::enthalpy::{
    CellId, DepositedCell, DepositionError, DepositionSimulation, EnergyAudit,
};
use std::collections::{BTreeMap, BTreeSet};

/// Spatial correspondence to one whole original caloric inventory.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RodLayerCell {
    pub segment: usize,
    pub layer: usize,
    pub cell: CellId,
}

/// Existing joule-coordinate nonlinear solve with independent physical gates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RodThermoelasticSolve {
    pub maximum_cells: usize,
    pub maximum_material_trials: usize,
    pub maximum_newton_iterations: usize,
    pub difference_step_j: f64,
    pub cell_energy_tolerance_j: f64,
    pub thermal_energy_tolerance_j: f64,
    /// Constitutive/global coupling defect, separate from native motion error.
    pub coupling_energy_tolerance_j: f64,
    pub newton: NewtonKrylovConfig,
}
impl RodThermoelasticSolve {
    fn driver(self) -> CircuitEnthalpySolve {
        CircuitEnthalpySolve {
            maximum_cells: self.maximum_cells,
            maximum_circuit_trials: self.maximum_material_trials,
            maximum_newton_iterations: self.maximum_newton_iterations,
            difference_step_j: self.difference_step_j,
            cell_energy_tolerance_j: self.cell_energy_tolerance_j,
            thermal_energy_tolerance_j: self.thermal_energy_tolerance_j,
            total_energy_tolerance_j: self.coupling_energy_tolerance_j,
            circuit: MaterialCircuitSolve {
                energy_tolerance_j: 0.0,
                power_tolerance_w: 0.0,
            },
            newton: self.newton,
        }
    }
}

#[derive(Debug)]
pub enum RodThermoelasticSourceError {
    Motion(RodMotionError),
    Material(MaterialStatePointError),
    Network(MaterialNetworkError),
}
impl core::fmt::Display for RodThermoelasticSourceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "rod thermoelastic source: {self:?}")
    }
}
impl std::error::Error for RodThermoelasticSourceError {}
pub type RodThermoelasticError = CircuitEnthalpyError<RodThermoelasticSourceError>;
type Result<T> = core::result::Result<T, RodThermoelasticError>;
fn invalid(s: &'static str) -> RodThermoelasticError {
    CircuitEnthalpyError::Invalid(s)
}
fn motion(e: RodMotionError) -> RodThermoelasticError {
    match e {
        RodMotionError::Solid(SpatialRodError::Cancelled) => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(RodThermoelasticSourceError::Motion(other)),
    }
}
fn solid(e: SpatialRodError) -> RodThermoelasticError {
    motion(e.into())
}
fn material(e: MaterialStatePointError) -> RodThermoelasticError {
    CircuitEnthalpyError::Source(RodThermoelasticSourceError::Material(e))
}
fn network(e: MaterialNetworkError) -> RodThermoelasticError {
    match e {
        MaterialNetworkError::Thermal(DepositionError::Cancelled) => {
            CircuitEnthalpyError::Cancelled
        }
        other => CircuitEnthalpyError::Source(RodThermoelasticSourceError::Network(other)),
    }
}
fn check(gate: &CancelGate) -> Result<()> {
    if gate.is_requested() {
        Err(CircuitEnthalpyError::Cancelled)
    } else {
        Ok(())
    }
}
fn finite(x: f64) -> Result<f64> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(invalid("unrepresentable rod thermoelastic arithmetic"))
    }
}
fn positive(x: f64) -> Result<f64> {
    if x.is_finite() && x > 0.0 {
        Ok(x)
    } else {
        Err(invalid("positive finite rod caloric capacity required"))
    }
}
fn sum(x: impl IntoIterator<Item = f64>) -> Result<f64> {
    x.into_iter().try_fold(0.0, |s, x| finite(s + x))
}
fn same_measure(a: f64, b: f64) -> bool {
    a.is_finite() && b.is_finite() && a > 0.0 && b > 0.0 && (a - b).abs() <= 1e-10 * b
}

struct Receiver {
    binding: RodLayerCell,
    original: DepositedCell,
    minimum_caloric_capacity_j_k: f64,
    temperature_range_k: [f64; 2],
}
struct Segment<'a> {
    law: ConcentricRodThermoelasticLaw<'a>,
    receivers: Vec<Receiver>,
    length_m: f64,
}

/// Whole-layer potentials at one original material/temperature/pose state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RodThermoelasticLayerState {
    pub destination: RodLayerCell,
    pub temperature_k: f64,
    pub free_energy_j: f64,
    pub entropy_j_k: f64,
    pub internal_energy_j: f64,
    pub mechanical_heat_capacity_j_k: f64,
}

/// Signed reversible exchange, never classified as physical dissipation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RodThermoelasticExchange {
    pub states: [RodThermoelasticLayerState; 2],
    /// Native conservative potential difference during frozen-T motion [J].
    /// The native integrator's separate global defect is not included here.
    pub frozen_strain_work_j: f64,
    pub coefficient_free_energy_change_j: f64,
    pub mechanical_entropy_energy_change_j: f64,
    pub reversible_heat_j: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RodThermoelasticEntropy {
    pub caloric_storage_change_j_k: f64,
    pub mechanical_storage_change_j_k: f64,
    pub reservoir_entropy_input_j_k: f64,
    pub conduction_generation_j_k: f64,
    /// Nonnegative endpoint sum Q_torsion/T_end [J/K], a production quadrature
    /// under the first-order split, not the exact continuous integral of P/T.
    pub torsional_generation_j_k: f64,
    pub external_source_at_cell_j_k: f64,
    pub has_unclassified_external_sources: bool,
    pub reversible_transfer_entropy_j_k: f64,
    /// Signed Delta S_mechanical + sum Q_rev/T_end; not entropy production.
    pub thermoelastic_temporal_defect_j_k: f64,
    pub thermal_temporal_defect_j_k: f64,
    pub algebraic_residual_j_k: f64,
    pub accounting_residual_j_k: f64,
}

/// One jointly accepted physical interval, including original numerical errors.
pub struct RodThermoelasticReport {
    /// Original frozen-temperature motion interval, before endpoint source refresh.
    pub mechanical: RodMotionStep,
    pub thermal: EnergyAudit,
    pub source: MaterialNetworkSample,
    pub sections: Vec<ConcentricRodThermoelasticState>,
    pub layers: Vec<RodThermoelasticExchange>,
    /// Frozen viscosity sources and actual original-layer heat deposits, if enabled.
    pub torsional: Option<RodTorsionalHeating>,
    pub entropy: RodThermoelasticEntropy,
    pub reversible_heat_j: f64,
    pub external_thermal_heat_j: f64,
    pub mechanical_internal_energy_before_j: f64,
    pub mechanical_internal_energy_after_j: f64,
    /// Raw complete energy defect; includes the native signed motion error.
    pub total_energy_residual_j: f64,
    /// Raw total minus the independently reported native integration defect.
    pub coupling_energy_residual_j: f64,
    pub nonlinear: NewtonReport,
    pub material_trials: usize,
    pub maximum_cell_residual_j: f64,
}

/// One rod's source binding and accepted mechanics; the original thermal owner
/// is supplied to each atomic step and can include additional passive cells.
pub struct RodThermoelasticSystem<'r, 'n, 'a, I = RodNodalInertia> {
    source: &'r RodMaterialSource<'a>,
    mechanical: RodMotion<I>,
    prepared: PreparedMaterialNetwork<'n, 'a>,
    segments: Vec<Segment<'a>>,
    states: Vec<RodThermoelasticLayerState>,
    sections: Vec<ConcentricRodThermoelasticState>,
    source_energy_tolerance_j: f64,
}

fn temperatures(segments: &[Segment<'_>], thermal: &DepositionSimulation<'_>) -> Result<Vec<f64>> {
    segments
        .iter()
        .flat_map(|s| &s.receivers)
        .map(|r| {
            let cell = thermal
                .cells()
                .get(&r.binding.cell)
                .ok_or_else(|| invalid("rod caloric cell disappeared"))?;
            let old = &r.original;
            if cell.id != old.id
                || cell.material_id != old.material_id
                || cell.body_id != old.body_id
                || cell.segment_id != old.segment_id
                || cell.mass_kg.to_bits() != old.mass_kg.to_bits()
                || cell.geometry != old.geometry
            {
                return Err(invalid(
                    "rod caloric inventory or reference geometry changed",
                ));
            }
            let t = thermal.temperature_k(cell.id)?;
            if t < r.temperature_range_k[0] || t > r.temperature_range_k[1] {
                return Err(invalid(
                    "rod temperature left its original solid caloric interval",
                ));
            }
            Ok(t)
        })
        .collect()
}

fn evaluate(
    segments: &[Segment<'_>],
    temperatures: &[f64],
    strains: &[RodSegmentStrainState],
    limit: f64,
    gate: &CancelGate,
) -> Result<(
    Vec<ConcentricRodThermoelasticState>,
    Vec<RodThermoelasticLayerState>,
)> {
    if strains.len() != segments.len()
        || temperatures.len() != segments.iter().map(|s| s.receivers.len()).sum::<usize>()
    {
        return Err(invalid(
            "complete original rod strain/temperature layout required",
        ));
    }
    let mut sections = Vec::with_capacity(segments.len());
    let mut states = Vec::with_capacity(temperatures.len());
    let mut offset = 0;
    for (segment, strain) in segments.iter().zip(strains) {
        check(gate)?;
        let end = offset + segment.receivers.len();
        let source = segment
            .law
            .at_strain(
                &temperatures[offset..end],
                strain.gamma,
                strain.curvature,
                limit,
            )
            .map_err(material)?;
        for (layer, receiver) in source.layers.iter().zip(&segment.receivers) {
            let l = segment.length_m;
            let mechanical_heat_capacity_j_k = finite(l * layer.heat_capacity_j_m_k)?;
            positive(receiver.minimum_caloric_capacity_j_k + mechanical_heat_capacity_j_k)?;
            states.push(RodThermoelasticLayerState {
                destination: receiver.binding,
                temperature_k: layer.temperature_k,
                free_energy_j: finite(l * layer.free_energy_j_m)?,
                entropy_j_k: finite(l * layer.entropy_j_m_k)?,
                internal_energy_j: finite(l * layer.internal_energy_j_m)?,
                mechanical_heat_capacity_j_k,
            });
        }
        sections.push(source);
        offset = end;
    }
    Ok((sections, states))
}

fn admit_native(
    segments: &[Segment<'_>],
    states: &[RodThermoelasticLayerState],
    native: &[RodSegmentStrainState],
    tolerance: f64,
) -> Result<()> {
    let mut offset = 0;
    for (segment, original) in segments.iter().zip(native) {
        let end = offset + segment.receivers.len();
        let energy = sum(states[offset..end].iter().map(|s| s.free_energy_j))?;
        if finite(energy - original.strain_energy_j)?.abs() > tolerance {
            return Err(invalid(
                "rod thermal potential disagrees with its original section energy",
            ));
        }
        offset = end;
    }
    Ok(())
}

impl<'r, 'n, 'a, I: RodMotionInertia> RodThermoelasticSystem<'r, 'n, 'a, I> {
    /// Bind every original segment/layer to its distinct whole caloric cell.
    /// Reference source contexts and cards, fixed density, actual mass/volume,
    /// a fully solid chart and positive combined capacity must agree. Additional
    /// passive cells are allowed. Circuits require an explicit electrical owner.
    pub fn try_new(
        mechanical: RodMotion<I>,
        source: &'r RodMaterialSource<'a>,
        thermal: &DepositionSimulation<'_>,
        network_source: &'n MaterialNetwork<'a>,
        destinations: &[RodLayerCell],
        source_energy_tolerance_j: f64,
        gate: &CancelGate,
    ) -> Result<Self> {
        check(gate)?;
        if !source_energy_tolerance_j.is_finite()
            || source_energy_tolerance_j < 0.0
            || !network_source.resistors.is_empty()
            || mechanical.time_s().to_bits() != thermal.time_s().to_bits()
            || destinations.len() != source.reference_temperatures_k().len()
        {
            return Err(invalid(
                "rod/heat requires one common clock, complete layers, thermal-only network and finite source tolerance",
            ));
        }
        let prepared = network_source
            .prepare(thermal, || gate.is_requested())
            .map_err(network)?;
        let routes: BTreeMap<_, _> = destinations
            .iter()
            .map(|d| ((d.segment, d.layer), d.cell))
            .collect();
        if routes.len() != destinations.len() {
            return Err(invalid("duplicate rod layer destination"));
        }
        let laws = source.thermoelastic_laws(gate).map_err(solid)?;
        let length = positive(source.reference_rod().reference_length_m() / laws.len() as f64)?;
        let mut segments = Vec::with_capacity(laws.len());
        let mut used = BTreeSet::new();
        for (i, law) in laws.into_iter().enumerate() {
            let mut receivers = Vec::with_capacity(law.axial_law().reference().layers().len());
            for (j, geometry) in law
                .axial_law()
                .reference()
                .layer_integrals()
                .iter()
                .enumerate()
            {
                check(gate)?;
                let cell = *routes
                    .get(&(i, j))
                    .ok_or_else(|| invalid("foreign or missing rod layer destination"))?;
                if !used.insert(cell) {
                    return Err(invalid(
                        "rod layers need distinct whole-cell thermal inventories",
                    ));
                }
                let bound = network_source
                    .cells
                    .iter()
                    .find(|b| b.cell == cell)
                    .ok_or_else(|| invalid("rod layer has no original caloric source"))?;
                let original = thermal
                    .cells()
                    .get(&cell)
                    .ok_or_else(|| invalid("inactive rod caloric cell"))?;
                let layer_source = &source.segments()[i].layers[j];
                let rho = positive(geometry.line_density_kg_m() / geometry.area_m2())?;
                if bound.point.material.content_hash() != layer_source.card.content_hash()
                    || bound.point.context != &layer_source.reference_point
                    || bound.point.temperature_axis != "T"
                    || !same_measure(original.mass_kg, length * geometry.line_density_kg_m())
                    || !same_measure(original.geometry.volume_m3, length * geometry.area_m2())
                    || bound.curve.knots().iter().any(|k| {
                        k.liquid_mass_fraction != 0.0 || !same_measure(k.bulk_density_kg_m3, rho)
                    })
                {
                    return Err(invalid(
                        "rod layer and caloric cell must share source context, whole fixed mass/volume and solid density",
                    ));
                }
                let knots = bound.curve.knots();
                let mut capacity = f64::INFINITY;
                for pair in knots.windows(2) {
                    capacity = capacity.min(positive(
                        original.mass_kg
                            * (pair[1].specific_enthalpy_j_kg - pair[0].specific_enthalpy_j_kg)
                            / positive(pair[1].temperature_k - pair[0].temperature_k)?,
                    )?);
                }
                receivers.push(Receiver {
                    binding: RodLayerCell {
                        segment: i,
                        layer: j,
                        cell,
                    },
                    original: original.clone(),
                    minimum_caloric_capacity_j_k: positive(capacity)?,
                    temperature_range_k: [
                        knots[0].temperature_k,
                        knots[knots.len() - 1].temperature_k,
                    ],
                });
            }
            segments.push(Segment {
                law,
                receivers,
                length_m: length,
            });
        }
        let temperatures = temperatures(&segments, thermal)?;
        if mechanical.model().section_temperatures_k().map_err(solid)? != temperatures {
            return Err(invalid(
                "rod and heat must begin at the same actual local temperatures",
            ));
        }
        // Replay exact original source selection before admitting a continuing
        // trajectory. This also gives unheated reference rods their zero binding.
        let (mechanical, _) = mechanical
            .rebind_material_temperatures_candidate(source, &temperatures, gate)
            .map_err(motion)?;
        let native = mechanical
            .model()
            .segment_strain_states(mechanical.pose(), gate)
            .map_err(solid)?;
        let (sections, states) = evaluate(
            &segments,
            &temperatures,
            &native,
            source.linear_strain_limit(),
            gate,
        )?;
        admit_native(&segments, &states, &native, source_energy_tolerance_j)?;
        // Bound the potentially negative axial capacity over EVERY layer's
        // whole caloric interval and admitted elastic strain, not just the start.
        for (segment, state) in segments.iter().zip(&sections) {
            for (receiver, layer) in segment.receivers.iter().zip(&state.axial.layers) {
                let [lo, hi] = receiver.temperature_range_k;
                let at = |t: f64| {
                    finite(
                        layer.young_modulus_pa
                            + layer.young_modulus_derivative_pa_k * (t - layer.temperature_k),
                    )
                };
                let emax = at(lo)?.max(at(hi)?).max(0.0);
                let curvature = finite(
                    (1.0 + 2048.0 * f64::EPSILON)
                        * receiver.original.geometry.volume_m3
                        * (emax * layer.alpha_per_k.powi(2)
                            + 2.0
                                * (layer.young_modulus_derivative_pa_k * layer.alpha_per_k).abs()
                                * source.linear_strain_limit()),
                )?;
                positive(receiver.minimum_caloric_capacity_j_k - hi * curvature)?;
            }
        }
        check(gate)?;
        Ok(Self {
            source,
            mechanical,
            prepared,
            segments,
            states,
            sections,
            source_energy_tolerance_j,
        })
    }

    pub const fn mechanical(&self) -> &RodMotion<I> {
        &self.mechanical
    }
    pub fn layers(&self) -> &[RodThermoelasticLayerState] {
        &self.states
    }
    pub fn sections(&self) -> &[ConcentricRodThermoelasticState] {
        &self.sections
    }
    /// Native kinetic, link/gravity/hydrostatic potential and source Psi+TS.
    /// The distinct original caloric cell energy is not included here.
    pub fn mechanical_internal_energy_j(&self) -> Result<f64> {
        let native = self.mechanical.observe([0.0; 3]).map_err(motion)?;
        finite(
            native.kinetic_energy_j
                + native.potential_energy_j
                + sum(self.states.iter().map(|s| s.temperature_k * s.entropy_j_k))?,
        )
    }

    /// One original representable time interval. Each thermal trial starts from
    /// the accepted heat history; mechanics is stepped once at accepted T. Both
    /// owners publish only after the full native/source/caloric/energy audit.
    /// A late failure changes neither owner and permits the same exact retry.
    pub fn step(
        &mut self,
        thermal: &mut DepositionSimulation<'_>,
        target_s: f64,
        loads: &[RodNodalWrench],
        mechanical_controls: RodMotionControls,
        controls: RodThermoelasticSolve,
        energy_guess_j: &[f64],
        gate: &CancelGate,
    ) -> Result<RodThermoelasticReport> {
        self.step_with_loads(
            thermal,
            target_s,
            StageLoads::World(loads),
            mechanical_controls,
            controls,
            energy_guess_j,
            gate,
        )
    }

    /// Advance heated rod motion with prepared material-point forces and
    /// body/world-held free couples. Both mechanical kicks evaluate the actual
    /// attachment frames through the original spatial-load owner.
    ///
    /// This shares the complete native/source/caloric solve and atomic commit
    /// with [`Self::step`]. Applied work and the signed mechanical integration
    /// defect retain their original accounts; neither becomes fictitious heat.
    /// Source refresh preserves the accepted endpoint pose and momenta. The
    /// thermoelastic split remains first order and requires refinement.
    pub fn step_spatial(
        &mut self,
        thermal: &mut DepositionSimulation<'_>,
        target_s: f64,
        loads: &PreparedRodLoads,
        mechanical_controls: RodMotionControls,
        controls: RodThermoelasticSolve,
        energy_guess_j: &[f64],
        gate: &CancelGate,
    ) -> Result<RodThermoelasticReport> {
        self.step_with_loads(
            thermal,
            target_s,
            StageLoads::Spatial(loads),
            mechanical_controls,
            controls,
            energy_guess_j,
            gate,
        )
    }

    fn step_with_loads(
        &mut self,
        thermal: &mut DepositionSimulation<'_>,
        target_s: f64,
        loads: StageLoads<'_>,
        mechanical_controls: RodMotionControls,
        controls: RodThermoelasticSolve,
        energy_guess_j: &[f64],
        gate: &CancelGate,
    ) -> Result<RodThermoelasticReport> {
        self.step_with_torsion(
            thermal,
            target_s,
            loads,
            None,
            mechanical_controls,
            controls,
            energy_guess_j,
            gate,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_with_torsion(
        &mut self,
        thermal: &mut DepositionSimulation<'_>,
        target_s: f64,
        loads: StageLoads<'_>,
        viscosity: Option<RodTorsionalViscosity<'_>>,
        mechanical_controls: RodMotionControls,
        controls: RodThermoelasticSolve,
        energy_guess_j: &[f64],
        gate: &CancelGate,
    ) -> Result<RodThermoelasticReport> {
        check(gate)?;
        let before = self.mechanical.time_s();
        if thermal.time_s().to_bits() != before.to_bits()
            || !target_s.is_finite()
            || target_s <= before
        {
            return Err(invalid(
                "rod and caloric history need the same advancing physical interval",
            ));
        }
        let old_temperatures = temperatures(&self.segments, thermal)?;
        if old_temperatures
            .iter()
            .zip(&self.states)
            .any(|(t, s)| t.to_bits() != s.temperature_k.to_bits())
        {
            return Err(invalid(
                "rod and caloric temperature histories were changed independently",
            ));
        }
        let dt = target_s - before;
        let internal_before = self.mechanical_internal_energy_j()?;
        let prepared_loss = viscosity
            .map(|selection| {
                torsional::prepare(
                    self.mechanical.model(),
                    self.source,
                    &self.segments,
                    selection,
                    gate,
                )
            })
            .transpose()?;
        let (frozen, mechanical) = self
            .mechanical
            .step_with_torsion_loads(
                dt,
                loads,
                prepared_loss.as_ref().map(|p| &p.loss),
                mechanical_controls,
                gate,
            )
            .map_err(motion)?;
        let torsional = prepared_loss
            .as_ref()
            .map(|prepared| {
                prepared.partition(
                    mechanical
                        .torsional_heat
                        .as_ref()
                        .ok_or_else(|| invalid("native torsional step lost its heat receipt"))?,
                    controls.coupling_energy_tolerance_j,
                    gate,
                )
            })
            .transpose()?;
        if frozen.time_s().to_bits() != target_s.to_bits() {
            return Err(invalid(
                "rod and thermal endpoint clocks do not represent the same interval",
            ));
        }
        let strains = frozen
            .model()
            .segment_strain_states(frozen.pose(), gate)
            .map_err(solid)?;
        let (_, frozen_states) = evaluate(
            &self.segments,
            &old_temperatures,
            &strains,
            self.source.linear_strain_limit(),
            gate,
        )?;
        admit_native(
            &self.segments,
            &frozen_states,
            &strains,
            self.source_energy_tolerance_j,
        )?;
        let (next_thermal, (), audit, entropy) = cell_driver_candidate_entropy_option(
            thermal,
            target_s,
            controls.driver(),
            energy_guess_j,
            |trial, _remaining| {
                check(gate)?;
                let mut transport = self
                    .prepared
                    .resolve(trial, || gate.is_requested())
                    .map_err(network)?;
                let mut external_source_entropy = 0.0;
                let mut has_external_source = false;
                for (&id, cell) in trial.cells() {
                    for power in [
                        cell.boundary.source_w,
                        transport.sources_w.get(&id).copied().unwrap_or(0.0),
                    ] {
                        external_source_entropy = finite(
                            external_source_entropy + dt * power / trial.temperature_k(id)?,
                        )?;
                        has_external_source |= power != 0.0;
                    }
                }
                let (sections, states) = evaluate(
                    &self.segments,
                    &temperatures(&self.segments, trial)?,
                    &strains,
                    self.source.linear_strain_limit(),
                    gate,
                )?;
                let mut exchanges = Vec::with_capacity(states.len());
                for (index, ((old, frozen), new)) in self
                    .states
                    .iter()
                    .zip(&frozen_states)
                    .zip(&states)
                    .enumerate()
                {
                    let parameter = finite(new.free_energy_j - frozen.free_energy_j)?;
                    let entropy_energy = finite(
                        new.temperature_k * new.entropy_j_k - old.temperature_k * old.entropy_j_k,
                    )?;
                    let reversible = finite(-parameter - entropy_energy)?;
                    let work = finite(frozen.free_energy_j - old.free_energy_j)?;
                    if finite(work - (new.internal_energy_j - old.internal_energy_j) - reversible)?
                        .abs()
                        > controls.coupling_energy_tolerance_j
                    {
                        return Err(invalid(
                            "rod reversible source does not close its original thermodynamic potential",
                        ));
                    }
                    let power = finite(reversible / dt)?;
                    if reversible != 0.0 && power == 0.0 {
                        return Err(invalid("rod reversible heat underflow"));
                    }
                    let p = transport
                        .sources_w
                        .entry(new.destination.cell)
                        .or_insert(0.0);
                    *p = finite(*p + power)?;
                    if let Some(heating) = &torsional {
                        let layer = &heating.layers[index];
                        if layer.destination != new.destination {
                            return Err(invalid(
                                "torsional heat cannot move to another caloric layer",
                            ));
                        }
                        if layer.heat_j != 0.0 {
                            let loss_power = finite(layer.heat_j / dt)?;
                            if loss_power <= 0.0 {
                                return Err(invalid(
                                    "positive torsional heat rate is unrepresentable",
                                ));
                            }
                            *p = finite(*p + loss_power)?;
                        }
                    }
                    exchanges.push(RodThermoelasticExchange {
                        states: [*old, *new],
                        frozen_strain_work_j: work,
                        coefficient_free_energy_change_j: parameter,
                        mechanical_entropy_energy_change_j: entropy_energy,
                        reversible_heat_j: reversible,
                    });
                }
                Ok(DrivenCoefficients {
                    circuit: CaloricDriver,
                    sources_w: transport.sources_w,
                    conductances_w_k: transport.conductances_w_k,
                    source: (
                        transport.source,
                        sections,
                        states,
                        exchanges,
                        external_source_entropy,
                        has_external_source,
                    ),
                    circuit_trials: 1,
                })
            },
            |_, _, checkpoint| {
                checkpoint()?;
                Ok(BTreeMap::new())
            },
            true,
            || gate.is_requested(),
        )?;
        let (source, sections, states, layers, external_source_entropy, has_external_source) =
            audit.source;
        let temperatures = states.iter().map(|s| s.temperature_k).collect::<Vec<_>>();
        let (candidate, _) = frozen
            .rebind_material_temperatures_candidate(self.source, &temperatures, gate)
            .map_err(motion)?;
        // The interval's loss/partition remains frozen. Check its SAME selected
        // viscosity claims at the actual thermal endpoint before publication.
        if let (Some(selection), Some(heating)) = (viscosity, torsional.as_ref()) {
            torsional::admit_endpoint(
                candidate.model(),
                self.source,
                &self.segments,
                selection,
                heating,
                gate,
            )?;
        }
        let native = candidate
            .model()
            .segment_strain_states(candidate.pose(), gate)
            .map_err(solid)?;
        admit_native(
            &self.segments,
            &states,
            &native,
            self.source_energy_tolerance_j,
        )?;
        let endpoint = candidate
            .observe(mechanical_controls.balance_origin_m)
            .map_err(motion)?;
        let internal_after = finite(
            endpoint.kinetic_energy_j
                + endpoint.potential_energy_j
                + sum(states.iter().map(|s| s.temperature_k * s.entropy_j_k))?,
        )?;
        let reversible_heat_j = sum(layers.iter().map(|l| l.reversible_heat_j))?;
        let mut external_thermal_heat_j =
            finite(audit.thermal.external_heat_j - reversible_heat_j)?;
        if let Some(heating) = &torsional {
            external_thermal_heat_j = finite(external_thermal_heat_j - heating.deposited_heat_j)?;
        }
        let total = finite(
            internal_after - internal_before + audit.thermal.stored_energy_change_j
                - external_thermal_heat_j
                - mechanical.source_work_quadrature_j,
        )?;
        let coupling = finite(total - mechanical.integration_defect_j)?;
        if coupling.abs() > controls.coupling_energy_tolerance_j
            || total.abs()
                > mechanical_controls.energy_defect_limit_j + controls.coupling_energy_tolerance_j
        {
            return Err(invalid(
                "rod/caloric storage does not close independently of native motion error",
            ));
        }
        let entropy = entropy
            .ok_or_else(|| invalid("rod thermoelastic step requires native caloric entropy"))?;
        let mechanical_entropy = sum(layers
            .iter()
            .map(|l| l.states[1].entropy_j_k - l.states[0].entropy_j_k))?;
        let reversible_entropy = sum(layers
            .iter()
            .map(|l| l.reversible_heat_j / l.states[1].temperature_k))?;
        let thermoelastic_defect = finite(mechanical_entropy + reversible_entropy)?;
        let torsional_generation = if let Some(heating) = &torsional {
            sum(heating
                .layers
                .iter()
                .zip(&states)
                .map(|(heat, state)| heat.heat_j / state.temperature_k))?
        } else {
            0.0
        };
        if torsional_generation < 0.0 {
            return Err(invalid(
                "negative intrinsic torsional entropy production quadrature",
            ));
        }
        let mut accounting = finite(
            entropy.thermal_storage_change_j_k + mechanical_entropy
                - entropy.reservoir_entropy_input_j_k
                - entropy.conduction_generation_j_k
                - external_source_entropy
                - entropy.temporal_defect_j_k
                - thermoelastic_defect
                - entropy.algebraic_residual_j_k,
        )?;
        if torsional.is_some() {
            accounting = finite(accounting - torsional_generation)?;
        }
        check(gate)?;
        let report = RodThermoelasticReport {
            mechanical,
            thermal: audit.thermal,
            source,
            sections: sections.clone(),
            layers,
            torsional,
            entropy: RodThermoelasticEntropy {
                caloric_storage_change_j_k: entropy.thermal_storage_change_j_k,
                mechanical_storage_change_j_k: mechanical_entropy,
                reservoir_entropy_input_j_k: entropy.reservoir_entropy_input_j_k,
                conduction_generation_j_k: entropy.conduction_generation_j_k,
                torsional_generation_j_k: torsional_generation,
                external_source_at_cell_j_k: external_source_entropy,
                has_unclassified_external_sources: has_external_source,
                reversible_transfer_entropy_j_k: reversible_entropy,
                thermoelastic_temporal_defect_j_k: thermoelastic_defect,
                thermal_temporal_defect_j_k: entropy.temporal_defect_j_k,
                algebraic_residual_j_k: entropy.algebraic_residual_j_k,
                accounting_residual_j_k: accounting,
            },
            reversible_heat_j,
            external_thermal_heat_j,
            mechanical_internal_energy_before_j: internal_before,
            mechanical_internal_energy_after_j: internal_after,
            total_energy_residual_j: total,
            coupling_energy_residual_j: coupling,
            nonlinear: audit.nonlinear,
            material_trials: audit.circuit_trials,
            maximum_cell_residual_j: audit.maximum_cell_residual_j,
        };
        // All fallible work and cancellation precede the three infallible moves.
        self.mechanical = candidate;
        self.states = states;
        self.sections = sections;
        *thermal = next_thermal;
        Ok(report)
    }
}

// This driver shares only the existing finite-cell residual. Native rod motion
// has already run once and is audited independently above; zero electrical
// storage/loss introduces no fictitious circuit or converted numerical heat.
struct CaloricDriver;
impl ThermalDriver for CaloricDriver {
    type Model = ();
    type Step = ();
    fn account(&self) -> DriverAccount<'_> {
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: 0.0,
            numerical_dissipation_j: 0.0,
            stored_change_j: 0.0,
            supplied_j: 0.0,
        }
    }
    fn into_parts(self) -> ((), ()) {
        ((), ())
    }
}
