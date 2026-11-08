//! Loaded native AC electrode circuits on the original finite enthalpy solve.
//!
//! Conductivity and external resistor material laws follow the same proposed
//! temperature. The affine field/circuit response supplies both spatial Joule
//! heat and the actual loaded resistor heat. Only original graph sources and
//! impressed field excitation supply external work; electrode work is internal.

use super::super::{
    SpatialMagneticMaterialSample, SpatialMagneticRegion, SpatialMagneticThermalNetwork, induction,
    native_error,
    terminals::circuit::{
        SpatialTerminalCircuitThermalNetwork, admit_initial_graph_sources, pinned_graph_sources,
    },
};
use super::{
    SpatialHarmonicHeatingStep, SpatialHarmonicThermalSample,
    field::{HarmonicFieldMaterials, recover_field_heat},
};
use crate::{
    circuit_enthalpy::{
        CircuitEnthalpyError, CircuitEnthalpyReport, CircuitEnthalpySolve, DrivenCoefficients,
        DriverAccount, ThermalDriver, admit_driver_unknowns, cell_driver_candidate,
        contacts::heat::{HeatAssembly, HeatNetworkError, closed_balance},
        finite,
        magnetic::{MagneticThermalError, MagneticThermalNetwork, check, invalid},
    },
    material_circuit::{
        MaterialCircuitError, admit_graph_with_field_inductors,
        harmonic::{MaterialHarmonicResistorHeat, account_harmonic_resistor_heat},
        lower_resistors_with_field_inductors,
    },
};
use fs_em::{
    planar::eddy::harmonic::circuit::HarmonicCircuitResponse,
    spatial::{
        SpatialMagneticMaterial,
        vector_potential::induction::{
            harmonic::terminals::{
                HarmonicCircuitControls, SpatialHarmonicTerminalControls,
                SpatialHarmonicTerminalError,
            },
            terminals::current::SpatialCurrentInduction,
        },
    },
};
use fs_exec::Cx;
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_material::conductor::spatial::SpatialConductorMaterial;
use fs_math::{c64::C64, dd::Dd};
use fs_phs::circuit::{CircuitGraph, MutualInductance};
use fs_thermal::enthalpy::{CellId, DepositionSimulation};
use std::{cell::RefCell, collections::BTreeMap};

/// Original field/material, loaded AC circuit, or complete heat-provider refusal.
#[derive(Debug)]
pub enum SpatialHarmonicCircuitThermalError<E = HeatNetworkError> {
    Native(MagneticThermalError),
    Circuit(SpatialHarmonicTerminalError),
    Heat(E),
}
impl<E> From<MagneticThermalError> for SpatialHarmonicCircuitThermalError<E> {
    fn from(error: MagneticThermalError) -> Self {
        Self::Native(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for SpatialHarmonicCircuitThermalError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "settled electrode circuit/finite heat: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for SpatialHarmonicCircuitThermalError<E> {}
type Error<E> = CircuitEnthalpyError<SpatialHarmonicCircuitThermalError<E>>;

fn electrical<E>(error: SpatialHarmonicTerminalError) -> Error<E> {
    match error {
        SpatialHarmonicTerminalError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(SpatialHarmonicCircuitThermalError::Circuit(other)),
    }
}
fn heat_error<E: From<HeatNetworkError>>(error: HeatNetworkError) -> Error<E> {
    match error {
        HeatNetworkError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => {
            CircuitEnthalpyError::Source(SpatialHarmonicCircuitThermalError::Heat(other.into()))
        }
    }
}

/// One actual loaded response at the accepted thermal endpoint.
/// Native field and circuit mean storage remain diagnostics, not changes in
/// transient energy. `circuit.field` retains the original conductor heat and
/// separate terminal/impressed account; total external work is reported here.
pub struct SpatialHarmonicCircuitHeatingStep {
    pub circuit: HarmonicCircuitResponse<SpatialHarmonicHeatingStep>,
    pub resistors: Vec<MaterialHarmonicResistorHeat>,
    /// Conductor plus original material-resistor heat over the interval [J].
    pub heat_j: f64,
    /// Original circuit-source plus impressed-field mean supply times dt [J].
    pub supplied_work_j: f64,
    /// Native port probes, optional independent excitation, and actual recovery.
    pub field_trials: usize,
}

/// Unpublished thermal endpoint and the complete settled field/circuit account.
pub struct SpatialHarmonicCircuitThermalCandidate<'t, S> {
    thermal: DepositionSimulation<'t>,
    report:
        CircuitEnthalpyReport<SpatialHarmonicThermalSample<S>, SpatialHarmonicCircuitHeatingStep>,
}
impl<'t, S> SpatialHarmonicCircuitThermalCandidate<'t, S> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn report(
        &self,
    ) -> &CircuitEnthalpyReport<SpatialHarmonicThermalSample<S>, SpatialHarmonicCircuitHeatingStep>
    {
        &self.report
    }
    /// Every reported conductor and resistor loss is already deposited once.
    /// The borrowed magnetic and graph templates have no new transient state.
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        CircuitEnthalpyReport<SpatialHarmonicThermalSample<S>, SpatialHarmonicCircuitHeatingStep>,
    ) {
        (self.thermal, self.report)
    }
}

struct CircuitHeatDriver(SpatialHarmonicCircuitHeatingStep);
impl ThermalDriver for CircuitHeatDriver {
    type Model = ();
    type Step = SpatialHarmonicCircuitHeatingStep;
    fn account(&self) -> DriverAccount<'_> {
        DriverAccount {
            // Peak complex receipts retain their own type. The shared cell
            // deposition callback below routes their physical mean powers.
            resistors: &[],
            contacts: &[],
            heat_j: self.0.heat_j,
            numerical_dissipation_j: 0.0,
            stored_change_j: 0.0,
            supplied_j: self.0.supplied_work_j,
        }
    }
    fn into_parts(self) -> ((), Self::Step) {
        ((), self.0)
    }
}

/// Heat a loaded physical-electrode AC circuit with complete material feedback.
///
/// Every original graph resistor must have one source-bound `CellResistor` in
/// the provider's material network. Its own whole-volume caloric body is the
/// unique heat receiver and must be separate from mapped magnetic bodies.
/// Numeric graph resistances are placeholders: the accepted thermal state and
/// every proposal re-resolve the SAME original material laws and geometry.
/// This permits consecutive intervals using the same immutable graph template.
///
/// `terminal_branches` identifies the native field-owned inductor placeholders
/// in original non-reference electrode order. Their numeric L is ignored;
/// actual ordinary L/C, mutual links, transformers and source indices retain
/// the original native circuit meaning. All phasors are peaks, exp(-i omega t).
/// Independent impressed excitation and graph sources are solved coherently.
///
/// The shared thermal residual receives both native cell Joule heat and the
/// original resistor's independently evaluated AC heat exactly once. Complete
/// contacts, radiation, layers and coolant use the SAME proposed temperatures.
/// Electrical interface contacts and unowned field-dependent heat boundaries
/// refuse. Ordinary resistor sources must be static scalars at held external
/// conditions; instantaneous E/J dependence cannot be frozen into this AC law.
///
/// `maximum_circuit_trials` counts every native field probe, optional impressed
/// trial and final actual recovery, including discarded Newton stencils. One
/// thermal proposal reserves that whole group BEFORE material/field preparation;
/// native descriptor work has its own separate budget. Material graph limits
/// remain 32 nodes, 64 original branches and 16 transformers.
///
/// This is stationary, fixed-permeability, settled-carrier/slow-thermal feedback.
/// Mean reactive storage is retained but its envelope change is not integrated.
/// Startup, fast modulation, magnetic parameter work, nonlinear waveforms and
/// moving thermal occupancy require their existing transient/physical owners.
/// No failure publishes a thermal endpoint or changes either borrowed template.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn heat_assembly_candidate<'t, A: HeatAssembly + ?Sized, C: SpatialConductorMaterial>(
    provider: &A,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    template: &SpatialCurrentInduction<'_, C>,
    graph: &CircuitGraph,
    terminal_branches: &[usize],
    external_mutual: &[MutualInductance],
    source_phasors: &[C64],
    impressed_current_a_m2: &[[C64; 3]],
    omega_rad_s: f64,
    regions: &[SpatialMagneticRegion],
    volume_tolerance_m3: f64,
    controls: CircuitEnthalpySolve,
    terminal_controls: SpatialHarmonicTerminalControls,
    circuit_controls: HarmonicCircuitControls,
    assembly: WeightedAssemblyLimits,
    initial_energy_increments_j: &[f64],
) -> Result<SpatialHarmonicCircuitThermalCandidate<'t, A::Sample>, Error<A::Error>> {
    let cancel = RefCell::new(|| cx.checkpoint().is_err());
    check(&cancel).map_err(native_error)?;
    admit_driver_unknowns(thermal, controls, initial_energy_increments_j)?;
    admit_graph_with_field_inductors(graph, terminal_branches)?;
    let dt = finite(target_s - thermal.time_s())?;
    let problem = template.problem();
    if dt <= 0.0
        || !omega_rad_s.is_finite()
        || omega_rad_s <= 0.0
        || impressed_current_a_m2.len() != problem.mesh().tets.len()
    {
        return Err(native_error(invalid(
            "positive AC thermal interval/frequency and complete impressed source required",
        )));
    }
    let mut excited = false;
    for cell in impressed_current_a_m2 {
        check(&cancel).map_err(native_error)?;
        for value in cell {
            if !value.re.is_finite() || !value.im.is_finite() {
                return Err(native_error(invalid("finite impressed AC source required")));
            }
            excited |= *value != C64::ZERO;
        }
    }
    let ports = template
        .electrode_references()
        .iter()
        .enumerate()
        .filter(|&(terminal, reference)| terminal != *reference)
        .count();
    if ports == 0 || ports != terminal_branches.len() || ports > terminal_controls.maximum_ports {
        return Err(native_error(invalid(
            "complete bounded original electrode port mapping required",
        )));
    }
    let field_trials = ports
        .checked_add(usize::from(excited))
        .and_then(|n| n.checked_add(1))
        .ok_or_else(|| native_error(invalid("AC field trial count overflow")))?;
    let prepared = provider.prepared_heat_network();
    let bulk = prepared.field_material_network().map_err(heat_error)?;
    prepared
        .admit_field_independent_boundaries(cx, &[])
        .map_err(heat_error)?;
    if bulk.cells.len() != thermal.cells().len() {
        return Err(native_error(invalid(
            "complete loaded AC thermal cells required",
        )));
    }
    let destinations = bulk
        .resistors
        .iter()
        .map(|r| r.destination)
        .collect::<Vec<_>>();
    let layout = SpatialTerminalCircuitThermalNetwork {
        magnetic: SpatialMagneticThermalNetwork {
            cells: bulk.cells,
            links: bulk.links,
            magnetic_regions: regions,
            volume_tolerance_m3,
        },
        graph_resistors: &destinations,
    };
    let graph_routes = layout
        .routes(thermal, thermal.time_s(), graph, &cancel)
        .map_err(native_error)?;
    layout
        .admit_graph_material_bindings(&[], Some(bulk.resistors), &graph_routes, &cancel)
        .map_err(native_error)?;
    let regions = layout.magnetic.regions();
    let network = MagneticThermalNetwork {
        cells: bulk.cells,
        links: bulk.links,
        magnetic_regions: &regions,
        volume_tolerance_m3,
    };
    let field_materials = HarmonicFieldMaterials::new(problem, template.conductors(), &cancel)
        .map_err(native_error)?;
    let field_routes = field_materials
        .routes(&network, thermal, &cancel)
        .map_err(native_error)?;
    let resolve = |state: &DepositionSimulation<'_>| {
        check(&cancel).map_err(native_error)?;
        let coefficients = provider
            .resolve_heat(cx, state)
            .map_err(|error| error.map_source(SpatialHarmonicCircuitThermalError::Heat))?;
        if !coefficients.contacts.is_empty() {
            return Err(native_error(invalid(
                "loaded AC heat cannot discard electrical contact work",
            )));
        }
        // Reject a frozen nonlinear/AC-fitted resistor before any field work.
        // The independent receipt audit repeats this source admission after
        // recovery; neither path substitutes a constant for the original law.
        for resistor in &coefficients.bulk.resistors {
            check(&cancel).map_err(native_error)?;
            if !resistor
                .conductor
                .source_property()
                .is_static_at_fixed_external_conditions()
            {
                return Err(MaterialCircuitError::Invalid(
                    "loaded AC heat requires instantaneous scalar resistor sources",
                )
                .into());
            }
        }
        let (sources, magnetic) = field_materials
            .resolve(
                &network,
                state,
                &coefficients.bulk.source.as_ref().material,
                &field_routes,
                &cancel,
            )
            .map_err(native_error)?;
        Ok((coefficients.bulk, sources, magnetic))
    };
    // The template is an AC topology, not a retained resistor-temperature state.
    // Lower from the accepted caloric sources before admitting their identity.
    let (accepted, _, _) = resolve(thermal)?;
    let (accepted_graph, accepted_sources) =
        lower_resistors_with_field_inductors(graph, accepted.resistors, terminal_branches)?;
    admit_initial_graph_sources(
        &accepted_graph,
        &[],
        graph_routes.len(),
        &accepted_sources,
        &cancel,
    )
    .map_err(native_error)?;
    let velocity = vec![[0.0; 3]; problem.mesh().vertex_count];
    let (thermal, (), mut report) = cell_driver_candidate(
        thermal,
        target_s,
        controls,
        initial_energy_increments_j,
        |trial, remaining| {
            if remaining < field_trials {
                return Err(CircuitEnthalpyError::CircuitBudget);
            }
            let (coefficients, sources, magnetic) = resolve(trial)?;
            pinned_graph_sources(&accepted_sources, &coefficients.resistors, &cancel)
                .map_err(native_error)?;
            let (loaded_graph, resistors) = lower_resistors_with_field_inductors(
                graph,
                coefficients.resistors,
                terminal_branches,
            )?;
            let mut refreshed = vec![None; field_materials.cell_count()];
            for (cell, source) in sources {
                refreshed[cell] = Some(source);
            }
            let field = template
                .prepare_harmonic_materials_trial(
                    magnetic
                        .iter()
                        .map(|s| SpatialMagneticMaterial::Linear(s.material.clone()))
                        .collect(),
                    refreshed,
                    assembly,
                    &mut || (*cancel.borrow_mut())(),
                )
                .map_err(induction)
                .map_err(native_error)?;
            let frequency = field
                .prepare_harmonic_frequency(
                    impressed_current_a_m2,
                    omega_rad_s,
                    &velocity,
                    terminal_controls,
                    &mut || (*cancel.borrow_mut())(),
                )
                .map_err(electrical)?;
            if frequency.preparation_field_trials() + 1 != field_trials {
                return Err(native_error(invalid(
                    "AC field preparation changed its admitted work count",
                )));
            }
            let response = frequency
                .circuit_response(
                    &loaded_graph,
                    terminal_branches,
                    external_mutual,
                    source_phasors,
                    circuit_controls,
                    &mut || (*cancel.borrow_mut())(),
                )
                .map_err(electrical)?;
            let field_step = recover_field_heat(&response.field.field, dt, &field_routes, &cancel)
                .map_err(native_error)?;
            let (resistors, graph_heat_j) = account_harmonic_resistor_heat(
                &loaded_graph,
                &resistors,
                &response.branches,
                response.mean_resistor_power_w,
                dt,
                controls.circuit,
                circuit_controls.current_tolerance_a,
                &mut || (*cancel.borrow_mut())(),
            )?;
            let heat_j = finite(field_step.heat_j + graph_heat_j)?;
            let supplied_work_j = finite(dt * response.supplied_power_va.re)?;
            let step = SpatialHarmonicCircuitHeatingStep {
                circuit: response.map_field(|_| field_step),
                resistors,
                heat_j,
                supplied_work_j,
                field_trials,
            };
            Ok(DrivenCoefficients {
                circuit: CircuitHeatDriver(step),
                sources_w: coefficients.sources_w,
                conductances_w_k: coefficients.conductances_w_k,
                circuit_trials: field_trials,
                source: SpatialHarmonicThermalSample {
                    network: coefficients.source,
                    magnetic: magnetic
                        .into_iter()
                        .map(|s| SpatialMagneticMaterialSample {
                            tetrahedron: s.triangle,
                            cell: s.cell,
                            material: s.material,
                        })
                        .collect(),
                },
            })
        },
        |driver: &CircuitHeatDriver, _, checkpoint| {
            let mut powers = BTreeMap::<CellId, Dd>::new();
            for cell in &driver.0.circuit.field.cells {
                checkpoint()?;
                let power = powers.entry(cell.cell).or_insert(Dd::ZERO);
                *power = *power + Dd::from_f64(cell.mean_ohmic_power_w);
            }
            for loss in &driver.0.resistors {
                checkpoint()?;
                let cell = *graph_routes.get(&loss.material.branch).ok_or_else(|| {
                    native_error(invalid(
                        "AC resistor heat lost its original finite receiver",
                    ))
                })?;
                let power = powers.entry(cell).or_insert(Dd::ZERO);
                *power = *power + Dd::from_f64(loss.mean_power_w);
            }
            powers
                .into_iter()
                .map(|(cell, power)| Ok((cell, finite(power.to_f64())?)))
                .collect()
        },
        || (*cancel.borrow_mut())(),
    )?;
    let sample = report.source.network.as_ref();
    closed_balance(
        dt,
        report.total_energy_residual_j,
        sample,
        controls.total_energy_tolerance_j,
    )
    .map_err(heat_error)?;
    let internal_net = finite(dt * sample.radiation_net_power_w())?;
    report.non_electrical_heat_j = finite(report.non_electrical_heat_j - internal_net)?;
    report.total_energy_residual_j = finite(report.total_energy_residual_j + internal_net)?;
    check(&cancel).map_err(native_error)?;
    Ok(SpatialHarmonicCircuitThermalCandidate { thermal, report })
}
