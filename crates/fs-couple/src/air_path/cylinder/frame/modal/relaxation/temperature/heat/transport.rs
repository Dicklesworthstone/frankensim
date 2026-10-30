//! Complete material heat transport on the existing relaxation transaction.
//! Source laws are resolved once per proposed endpoint. The original enthalpy,
//! circuit and mechanical owners still advance and publish every physical state.
use super::*;
use crate::circuit_enthalpy::contacts::{ContactThermalCoefficients, heat::{
    HeatNetworkError, HeatNetworkSample,
}};
use crate::circuit_enthalpy::CircuitThermalCoefficients;

/// The same borrowed prepared heat assembly and execution authority used by
/// source-defined thermoelastic frames. Retain its preparation to retain pins.
pub use crate::frame_vibration::thermoelastic::{
    FrameHeatInput as RelaxationHeatInput, FrameHeatTransfer,
};

fn source_error(error: HeatNetworkError, gate: &CancelGate) -> RelaxationHeatError {
    match error {
        HeatNetworkError::Cancelled => {
            gate.request();
            RelaxationHeatError::Cancelled
        }
        other => RelaxationHeatError::Source(RelaxationHeatSourceError::Heat(other)),
    }
}

/// A complete evaluation by an existing physical heat owner, plus its own
/// accepted source payload. The first item retains the original heat network.
pub type RelaxationHeatCoefficients<S> = ContactThermalCoefficients<(HeatNetworkSample, S)>;

/// Keep the caller's typed transport failure distinct from native failures.
#[derive(Debug)]
pub enum RelaxationTransportError<E> {
    Source(E),
    Relaxation(RelaxationHeatError),
}
impl<E: core::fmt::Debug> core::fmt::Display for RelaxationTransportError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "relaxation transport: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for RelaxationTransportError<E> {}

// This private adapter is not a solver or new material interface. It allows
// higher-layer owners to supply the original residual without a dependency cycle.
pub(super) trait Source<S> {
    fn checkpoint(&self, gate: &CancelGate) -> Result<()>;
    fn resolve(&self, thermal: &DepositionSimulation<'_>, gate: &CancelGate)
        -> Result<RelaxationHeatCoefficients<S>>;
}

pub(super) fn checkpoint<S>(input: Option<&dyn Source<S>>, gate: &CancelGate) -> Result<()> {
    match input { Some(input) => input.checkpoint(gate), None => check(gate) }
}

const SOURCE_REFUSAL: &str = "original relaxation heat source callback refused";
struct Callback<'p, 'n, 'a, 'scope, R, E> {
    heat: RelaxationHeatInput<'p, 'n, 'a, 'scope>,
    resolve: R,
    failure: std::cell::RefCell<Option<E>>,
}
impl<S, R, E> Source<S> for Callback<'_, '_, '_, '_, R, E>
where R: Fn(&DepositionSimulation<'_>) -> core::result::Result<RelaxationHeatCoefficients<S>, E> {
    fn checkpoint(&self, gate: &CancelGate) -> Result<()> {
        if self.heat.cx.checkpoint().is_err() { gate.request(); }
        check(gate)
    }
    fn resolve(&self, thermal: &DepositionSimulation<'_>, gate: &CancelGate)
        -> Result<RelaxationHeatCoefficients<S>> {
        self.checkpoint(gate)?;
        let _ = self.failure.borrow_mut().take();
        let result = (self.resolve)(thermal);
        self.checkpoint(gate)?;
        result.map_err(|error| {
            *self.failure.borrow_mut() = Some(error);
            invalid(SOURCE_REFUSAL)
        })
    }
}

impl<F> FrameRelaxationThermalSystem<'_, '_, '_, '_, F> {
    /// Continue the SAME material-memory/fluid/circuit/caloric histories with
    /// source-defined thermal contacts, nonlinear conducting layers and closed
    /// gray-diffuse radiation. The prepared assembly replaces bulk-only source
    /// evaluation; it is not added as a second heat solve or second conductance.
    ///
    /// Original cell/link slices, accepted electrical sources, physical queries,
    /// geometry and two-sided Joule partitions must agree before motion starts.
    /// An electrically bound owner requires `Some(input)`, even at zero drive.
    /// An unpowered owner takes `None`; electrical work cannot be dropped.
    ///
    /// All transports and R(T) see the same proposed cell temperatures. Maxwell
    /// loss is independently integrated at accepted start rates; endpoint tau(T)
    /// acts next interval. Enclosure roundoff is not external energy, and its
    /// material Q/T transfer is separate from unclassified authored sources.
    /// This is not a photon-entropy model or a new time-integration scheme.
    ///
    /// Geometry, view factors and mechanical storage moduli/density stay fixed.
    /// Source C(T)/L(T) requires the original field/thermal binding; otherwise
    /// authored electrical storage remains fixed rather than being inferred.
    /// Thermal paths do not infer deformation, opening contacts or melting.
    /// Source/solver/budget failures and either cancellation authority return
    /// before the existing joint publication. Explicitly choosing a different
    /// prepared heat assembly is a caller intervention, not a source fallback.
    #[allow(clippy::too_many_arguments)]
    pub fn step_heat_network(
        &mut self,
        thermal: &mut DepositionSimulation<'_>,
        external: &[f64],
        heat: RelaxationHeatInput<'_, '_, '_, '_>,
        electrical_input: Option<RelaxationElectricalInput<'_>>,
        solve: FrameThermalSolve,
        energy_guess_j: &[f64],
        gate: &CancelGate,
    ) -> Result<RelaxationHeatReport> {
        match self.step_heat_network_with(thermal, external, heat, electrical_input,
            solve, energy_guess_j, gate, |state| {
                let ContactThermalCoefficients { bulk, contacts } = heat.network.resolve(heat.cx, state)?;
                Ok::<_, HeatNetworkError>(ContactThermalCoefficients { contacts,
                    bulk: CircuitThermalCoefficients { resistors: bulk.resistors,
                        sources_w: bulk.sources_w, conductances_w_k: bulk.conductances_w_k,
                        source: (bulk.source, ()) } })
            }) {
            Ok(report) => Ok(report),
            Err(RelaxationTransportError::Source(error)) => Err(source_error(error, gate)),
            Err(RelaxationTransportError::Relaxation(error)) => Err(error),
        }
    }

    /// Compose another existing heat-source owner with this SAME mechanical,
    /// field, electrical and finite-enthalpy transaction. The source adapter
    /// must resolve `heat.network` at the supplied trial, add each physical
    /// exchange once, and retain its own payload beside the original sample.
    /// For example, fs-airflow uses this for material-resolved reservoir cooling.
    ///
    /// The callback is repeatable, has no accepted mutable physical state, and
    /// must preserve original source identity, context, geometry and support.
    /// It runs once at admission and once per original thermal source evaluation,
    /// NOT after publication. Its final accepted payload is returned in
    /// `transport_source`; rejected trials are never presented as accepted data.
    /// Extra sources are signed watts; extra finite transfers belong in the
    /// original reciprocal link map, never in two independent source entries.
    /// This seam does not authorize a second solver or skipped thermal storage.
    ///
    /// Original circuit and source correspondence, source-evaluation budgets,
    /// field reversible heat, radiation and joint rollback remain unchanged.
    /// Authored/callback external heat remains explicitly unclassified in the
    /// generic entropy report; a physical source payload may describe its actual
    /// reservoirs. Q/T_cell is not automatically reservoir entropy production.
    #[allow(clippy::too_many_arguments)]
    pub fn step_heat_network_with<S, E>(
        &mut self,
        thermal: &mut DepositionSimulation<'_>,
        external: &[f64],
        heat: RelaxationHeatInput<'_, '_, '_, '_>,
        electrical_input: Option<RelaxationElectricalInput<'_>>,
        solve: FrameThermalSolve,
        energy_guess_j: &[f64],
        gate: &CancelGate,
        resolve: impl Fn(&DepositionSimulation<'_>)
            -> core::result::Result<RelaxationHeatCoefficients<S>, E>,
    ) -> core::result::Result<RelaxationHeatReport<S>, RelaxationTransportError<E>> {
        let callback = Callback { heat, resolve, failure: std::cell::RefCell::new(None) };
        let result = (|| -> Result<RelaxationHeatReport<S>> {
            Source::<S>::checkpoint(&callback, gate)?;
            crate::circuit_enthalpy::admit_driver_unknowns(thermal, controls(solve), energy_guess_j)?;
            if self.electrical.is_some() != electrical_input.is_some() {
                return Err(invalid("electrically bound relaxation requires explicit held electrical inputs"));
            }
            let original = heat.network.field_storage_material_network();
            if !core::ptr::eq(self.source.cells, original.cells)
                || !core::ptr::eq(self.source.links, original.links) {
                return Err(invalid("relaxation and complete heat require original cell/link slices"));
            }
            let mut initial = Source::<S>::resolve(&callback, thermal, gate)?;
            let destinations = heat.network.resistor_destinations();
            let routes = heat.network.heat_destinations(thermal)
                .map_err(|error| source_error(error, gate))?;
            electrical::admit_heat_sources(self.electrical.as_ref(), &mut initial,
                &destinations, &routes, thermal, gate)?;
            drop(initial);
            Source::<S>::checkpoint(&callback, gate)?;
            self.step_inner(thermal, external, electrical_input, Some(&callback),
                solve, energy_guess_j, gate)
        })();
        match result {
            Err(CircuitEnthalpyError::Invalid(SOURCE_REFUSAL)) => match callback.failure.into_inner() {
                Some(error) => Err(RelaxationTransportError::Source(error)),
                None => Err(RelaxationTransportError::Relaxation(
                    invalid("heat-source refusal lost its original typed cause"))),
            },
            other => other.map_err(RelaxationTransportError::Relaxation),
        }
    }
}

pub(super) struct Sample {
    transfer: FrameHeatTransfer,
    independent_entropy_j_k: f64,
    has_independent_sources: bool,
}

pub(super) fn resolve<S>(
    input: Option<&dyn Source<S>>,
    electrical: Option<&electrical::Binding<'_, '_>>,
    prepared: &PreparedMaterialNetwork<'_, '_>,
    thermal: &DepositionSimulation<'_>,
    dt: f64,
    gate: &CancelGate,
) -> Result<(ContactThermalCoefficients<MaterialNetworkSample>, Option<Sample>, Option<S>)> {
    checkpoint(input, gate)?;
    let Some(input) = input else {
        return Ok((electrical::resolve(electrical, prepared, thermal, gate)?, None, None));
    };
    let resolved = input.resolve(thermal, gate)?;
    let ContactThermalCoefficients { bulk, contacts } = resolved;
    let (HeatNetworkSample { material, thermal_contacts, integrated_conduction, radiation }, payload) = bulk.source;
    // Capture authored heat separately BEFORE the native driver combines it
    // with Maxwell/Joule losses. Internal radiative transfer is already in the
    // source map, but never establishes an independent external heat source.
    let mut independent_entropy_j_k = 0.0;
    let mut has_independent_sources = false;
    for (&cell, state) in thermal.cells() {
        checkpoint(Some(input), gate)?;
        let radiative = radiation.as_ref()
            .and_then(|r| r.sources_w().get(&cell)).copied().unwrap_or(0.0);
        let independent = finite(bulk.sources_w.get(&cell).copied().unwrap_or(0.0) - radiative)?;
        let temperature = thermal.temperature_k(cell)?;
        for power in [state.boundary.source_w, independent] {
            if power != 0.0 { has_independent_sources = true; }
            independent_entropy_j_k = finite(independent_entropy_j_k
                + source_entropy(power, dt, temperature)?)?;
        }
    }
    checkpoint(Some(input), gate)?;
    Ok((ContactThermalCoefficients { contacts,
        bulk: CircuitThermalCoefficients { resistors: bulk.resistors,
            sources_w: bulk.sources_w, conductances_w_k: bulk.conductances_w_k,
            source: material },
    }, Some(Sample { transfer: FrameHeatTransfer {
        thermal_contacts, integrated_conduction, radiation,
    }, independent_entropy_j_k, has_independent_sources }), Some(payload)))
}

pub(super) fn source_entropy(power: f64, dt: f64, temperature: f64) -> Result<f64> {
    let heat = finite(dt * power)?;
    let entropy = finite(heat / temperature)?;
    if power != 0.0 && (heat == 0.0 || entropy == 0.0) {
        return Err(invalid("material heat-transfer entropy underflow"));
    }
    Ok(entropy)
}

pub(super) fn finish(sample: Option<Sample>, thermal: &DepositionSimulation<'_>,
    dt: f64, entropy: &mut FrameHeatEntropy, gate: &CancelGate,
) -> Result<(Option<FrameHeatTransfer>, f64, f64)> {
    let Some(sample) = sample else { return Ok((None, 0.0, 0.0)); };
    let mut radiation_entropy = 0.0;
    let radiation_net_j = match &sample.transfer.radiation {
        Some(radiation) => {
            for (&cell, &power) in radiation.sources_w() {
                check(gate)?;
                radiation_entropy = finite(radiation_entropy
                    + source_entropy(power, dt, thermal.temperature_k(cell)?)?)?;
            }
            finite(dt * radiation.net_heat_w())?
        }
        None => 0.0,
    };
    entropy.unclassified_source_at_cell_j_k = sample.independent_entropy_j_k;
    entropy.has_unclassified_sources = sample.has_independent_sources;
    // Direct original terms, not a residual relabeled as physical entropy.
    // The separate report radiation term completes this material account.
    entropy.accounting_residual_j_k = [entropy.thermal_storage_change_j_k,
        -entropy.reservoir_entropy_input_j_k, -entropy.conduction_generation_j_k,
        -entropy.material_loss_generation_j_k, -sample.independent_entropy_j_k,
        -radiation_entropy, -entropy.temporal_defect_j_k, -entropy.algebraic_residual_j_k,
    ].into_iter().try_fold(0.0, |sum, value| finite(sum + value))?;
    check(gate)?;
    Ok((Some(sample.transfer), radiation_net_j, radiation_entropy))
}
