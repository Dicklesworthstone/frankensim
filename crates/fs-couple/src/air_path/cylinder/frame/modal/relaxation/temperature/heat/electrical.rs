//! Original circuit storage, Joule heat and source relaxation on ONE thermal solve.
//!
//! The circuit supplies an unpublished candidate for each actual caloric trial.
//! Its accepted charge/flux never advances in a Newton stencil. The existing
//! caloric driver adds its original work/storage/heat to the mechanical account;
//! no RC integrator, material evaluator or second temperature is introduced.

mod contact;
/// Original material C(T)/L(T), reversible field heat and conserved coordinates.
pub mod field;
use crate::circuit_enthalpy::field_storage::FieldThermalBinding;
use super::*;
use crate::circuit_enthalpy::contacts::ContactThermalCoefficients;
use crate::material_circuit::contacts::SeriesElectricalContact;
use crate::material_circuit::{
    MaterialCircuit, MaterialCircuitCandidate, MaterialCircuitStep, MaterialResistor,
};
use fs_material::conductor::ElectricalPropertyBasis;
use fs_thermal::enthalpy::{CellId, DepositionError};

/// Reuse the existing held-port and native electrical-tolerance input.
/// Empty ports mean the circuit is powered solely by its finite C/L stores.
pub use crate::frame_vibration::thermoelastic::FrameElectricalInput as RelaxationElectricalInput;

pub(super) struct Binding<'n, 'a> {
    pub(super) circuit: MaterialCircuit,
    destinations: BTreeMap<usize, CellId>,
    contacts: Option<contact::Binding<'n, 'a>>,
    fields: Vec<FieldThermalBinding>,
}

impl<'m, 'f, 'n, 'a, F> FrameRelaxationThermalSystem<'m, 'f, 'n, 'a, F> {
    /// Cold-bind native electrical storage to the SAME original caloric cells.
    /// `source` adds resistors to the original network, retaining the exact
    /// borrowed cell/link slices. Resistors may heat selected relaxation layers
    /// or other existing cells; the original heat network conducts between them.
    ///
    /// Every bulk resistor must match the initial resolved source and geometry.
    /// Prepare the circuit using `source.prepare_with_bases(...).resolve(...)`.
    /// The bases are explicit; conductivity-only material cards need no invented
    /// resistivity claim. Any refusal leaves this owner and all histories intact.
    ///
    /// Explicitly authored C/L and electrical geometry are held fixed. Material-
    /// bound field storage refuses without its reciprocal thermal-work owner.
    /// This binding does not infer
    /// temperature-dependent field storage, motional EMF or electromagnetic
    /// forces. Electrical contacts need their own two-sided source/heat owner
    /// and use `bind_material_contact_circuit`; this bulk-only entry point refuses
    /// them rather than silently discarding interface heat.
    /// Rate-only solid storage and the original first-order thermal split apply.
    pub fn bind_material_circuit(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        source: &'n MaterialNetwork<'a>,
        circuit: MaterialCircuit,
        electrical_bases: &[(usize, ElectricalPropertyBasis)],
        gate: &CancelGate,
    ) -> Result<()> {
        self.bind_bulk_circuit(thermal, source, circuit, electrical_bases, Vec::new(), gate)
    }

    #[allow(clippy::too_many_arguments)]
    fn bind_bulk_circuit(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        source: &'n MaterialNetwork<'a>,
        circuit: MaterialCircuit,
        electrical_bases: &[(usize, ElectricalPropertyBasis)],
        fields: Vec<FieldThermalBinding>,
        gate: &CancelGate,
    ) -> Result<()> {
        self.admit_electrical_binding(thermal, source, &circuit, &fields, gate)?;
        if !circuit.contacts().is_empty() {
            return Err(invalid("electrical contact heat requires explicit original two-sided receivers"));
        }
        let prepared = source.prepare_with_bases(thermal, electrical_bases,
            || gate.is_requested()).map_err(network)?;
        let mut resolved = prepared.resolve(thermal, || gate.is_requested())
            .map_err(network)?.resistors;
        resolved.sort_unstable_by_key(|resistor| resistor.branch);
        if resolved.as_slice() != circuit.resistors() {
            return Err(invalid("electrical relaxation must retain every original initial resistor source and geometry"));
        }
        let destinations = source.resistors.iter()
            .map(|binding| (binding.destination.branch, binding.destination.cell)).collect();
        check(gate)?;
        // No fallible work follows these moves; material/source/clock identity
        // is shared, not rebuilt by a second network or exposed mutable circuit.
        self.prepared = prepared;
        self.source = source;
        self.electrical = Some(Binding { circuit, destinations, contacts: None, fields });
        Ok(())
    }

    /// Read-only accepted electrical state. It cannot be independently stepped.
    pub fn electrical(&self) -> Option<&MaterialCircuit> {
        self.electrical.as_ref().map(|binding| &binding.circuit)
    }

    /// One original mechanical sample and ONE caloric residual containing
    /// independently integrated Maxwell heat and native circuit Joule heat.
    /// R(T) and conduction use the actual proposed endpoint; source tau(T) is
    /// installed only after complete acceptance and acts on the NEXT sample.
    ///
    /// Every electrical trial begins at the same accepted charge/flux. Failure
    /// publishes no circuit, finite-cell, material-memory, fluid or clock state.
    /// Ordinary `step` refuses an electrically bound owner even at zero drive:
    /// a zero port value does not imply zero stored electrical energy.
    #[allow(clippy::too_many_arguments)]
    pub fn step_electrical(
        &mut self,
        thermal: &mut DepositionSimulation<'_>,
        external: &[f64],
        input: RelaxationElectricalInput<'_>,
        solve: FrameThermalSolve,
        energy_guess_j: &[f64],
        gate: &CancelGate,
    ) -> Result<RelaxationHeatReport> {
        self.step_inner(thermal, external, Some(input), None, solve, energy_guess_j, gate)
    }
}

// A complete heat evaluator may change thermal transport, never electrical
// source identity, accepted physical state, geometry or heat destinations.
pub(super) fn admit_heat_sources<S>(
    binding: Option<&Binding<'_, '_>>,
    coefficients: &mut ContactThermalCoefficients<S>,
    destinations: &[crate::circuit_enthalpy::ResistorCell],
    routes: &[crate::circuit_enthalpy::contacts::ContactHeatCells],
    thermal: &DepositionSimulation<'_>,
    gate: &CancelGate,
) -> Result<()> {
    use crate::circuit_enthalpy::contacts;
    check(gate)?;
    match binding {
        Some(binding) => {
            coefficients.bulk.resistors.sort_unstable_by_key(|r| r.branch);
            if coefficients.bulk.resistors.as_slice() != binding.circuit.resistors() {
                return Err(invalid("complete heat changed accepted relaxation resistor sources, state or geometry"));
            }
            contacts::admit_sources(&binding.circuit, &coefficients.contacts, || gate.is_requested())?;
            coefficients.contacts.sort_unstable_by_key(|c| c.id());
            for (actual, expected) in binding.circuit.contacts().iter().zip(&coefficients.contacts) {
                check(gate)?;
                // Pinning policy receipts may differ; the original selected
                // source (above), physical query and reduction may not.
                if actual.id() != expected.id() || actual.branch() != expected.branch()
                    || actual.area_m2().to_bits() != expected.area_m2().to_bits()
                    || actual.resistance_ohm().to_bits() != expected.resistance_ohm().to_bits()
                    || actual.source().resolved().query() != expected.source().resolved().query()
                {
                    return Err(invalid("complete heat changed accepted relaxation electrical interface state or area"));
                }
            }
            let map:BTreeMap<_,_> = destinations.iter().map(|r| (r.branch, r.cell)).collect();
            if map.len() != destinations.len() || map != binding.destinations {
                return Err(invalid("complete heat changed original relaxation Joule receivers"));
            }
            let routes = contacts::admit_routes(thermal, &binding.circuit, routes, || {
                if gate.is_requested() { Err(DepositionError::Cancelled) } else { Ok(()) }
            })?;
            let same = match &binding.contacts {
                Some(original) => routes == original.routes,
                None => routes.is_empty(),
            };
            if !same {
                return Err(invalid("complete heat changed original relaxation contact heat partitions"));
            }
        }
        None => {
            if !coefficients.bulk.resistors.is_empty() || !coefficients.contacts.is_empty()
                || !destinations.is_empty() || !routes.is_empty()
            {
                return Err(invalid("unpowered relaxation cannot discard heat-network electrical work"));
            }
        }
    }
    check(gate)
}

// Keep original source selection and contact context evaluation in their
// published owner. Both paths feed the same native circuit/caloric transaction.
pub(super) fn resolve(
    binding: Option<&Binding<'_, '_>>,
    prepared: &PreparedMaterialNetwork<'_, '_>,
    thermal: &DepositionSimulation<'_>,
    gate: &CancelGate,
) -> Result<ContactThermalCoefficients<MaterialNetworkSample>> {
    check(gate)?;
    match binding.and_then(|binding| binding.contacts.as_ref()) {
        Some(contact) => contact.prepared.resolve(thermal, || gate.is_requested())
            .map_err(contact::source_error),
        None => Ok(ContactThermalCoefficients {
            bulk: prepared.resolve(thermal, || gate.is_requested()).map_err(network)?,
            contacts: Vec::new(),
        }),
    }
}

pub(super) fn advance(
    binding: Option<&Binding<'_, '_>>,
    input: Option<RelaxationElectricalInput<'_>>,
    resistors: Vec<MaterialResistor>,
    contacts: Vec<SeriesElectricalContact>,
    duration: f64,
    gate: &CancelGate,
) -> Result<Option<MaterialCircuitCandidate>> {
    check(gate)?;
    match (binding, input) {
        (Some(binding), Some(input)) => Ok(Some(
            crate::circuit_enthalpy::contacts::advance_circuit(
                &binding.circuit, resistors, contacts, input.ports, duration,
                input.solve, || gate.is_requested(),
            )?,
        )),
        (None, None) if resistors.is_empty() && contacts.is_empty() => Ok(None),
        _ => Err(invalid("electrical source and accepted relaxation owner disagree")),
    }
}

// Add ORIGINAL driver accounts; the shared caloric residual independently
// checks all deposited powers against the sum and the global work balance.
pub(super) struct ElectricalDriver {
    mechanical: Driver,
    electrical: Option<MaterialCircuitCandidate>,
}
impl ElectricalDriver {
    pub(super) fn new(mechanical: Driver, electrical: Option<MaterialCircuitCandidate>) -> Result<Self> {
        let driver = Self { mechanical, electrical };
        let a = driver.account();
        finite(a.heat_j)?;
        finite(a.stored_change_j)?;
        finite(a.supplied_j)?;
        finite(a.numerical_dissipation_j)?;
        Ok(driver)
    }
}
impl ThermalDriver for ElectricalDriver {
    type Model = Option<MaterialCircuit>;
    type Step = (CoupledModalFrame, Option<MaterialCircuitStep>);
    fn account(&self) -> DriverAccount<'_> {
        let mut account = self.mechanical.account();
        if let Some(candidate) = &self.electrical {
            let electrical = candidate.account();
            account.resistors = electrical.resistors;
            account.contacts = electrical.contacts;
            account.heat_j += electrical.heat_j;
            account.stored_change_j += electrical.stored_change_j;
            account.supplied_j += electrical.supplied_j;
            account.numerical_dissipation_j += electrical.numerical_dissipation_j;
        }
        account
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        let (_, mechanical) = self.mechanical.into_parts();
        match self.electrical {
            Some(candidate) => {
                let (circuit, report) = candidate.into_parts();
                (Some(circuit), (mechanical, Some(report)))
            }
            None => (None, (mechanical, None)),
        }
    }
}

pub(super) fn deposit(
    binding: Option<&Binding<'_, '_>>,
    driver: &ElectricalDriver,
    duration: f64,
    powers: &mut BTreeMap<CellId, f64>,
    checkpoint: &mut dyn FnMut() -> core::result::Result<(), DepositionError>,
) -> Result<()> {
    match (binding, &driver.electrical) {
        (Some(binding), Some(candidate)) => {
            if candidate.step().resistors.len() != binding.destinations.len()
            {
                return Err(invalid("electrical heat coverage changed in the shared relaxation trial"));
            }
            for branch in &candidate.step().resistors {
                checkpoint()?;
                let cell = binding.destinations.get(&branch.material.branch)
                    .ok_or_else(|| invalid("electrical heat lost its original caloric receiver"))?;
                let power = finite(branch.heat_j / duration)?;
                let previous = powers.get(cell).copied().unwrap_or(0.0);
                let next = finite(previous + power)?;
                if branch.heat_j > 0.0 && (power == 0.0 || next == previous) {
                    return Err(invalid("electrical heat disappears in the shared receiving power"));
                }
                powers.insert(*cell, next);
            }
            match &binding.contacts {
                Some(contact) => crate::circuit_enthalpy::contacts::deposit(
                    &candidate.step().contacts, &contact.routes, duration, powers, checkpoint,
                ),
                None if candidate.step().contacts.is_empty() => Ok(()),
                None => Err(invalid("electrical interface losses lack their original two-sided receivers")),
            }
        }
        (None, None) => Ok(()),
        _ => Err(invalid("electrical candidate and heat receiver owner disagree")),
    }
}
