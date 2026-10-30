//! Original ordered interface resistance and two-sided Joule heat.
use super::*;
use crate::circuit_enthalpy::contacts::{
    self, ContactHeatCells,
    material::{ElectricalContactNetwork, ElectricalContactNetworkError,
        PreparedElectricalContactNetwork},
};

pub(super) struct Binding<'n, 'a> {
    pub(super) prepared: PreparedElectricalContactNetwork<'n, 'a>,
    pub(super) routes: BTreeMap<u64, ContactHeatCells>,
}

pub(super) fn source_error(error: ElectricalContactNetworkError) -> RelaxationHeatError {
    match error {
        ElectricalContactNetworkError::Bulk(error) => network(error),
        ElectricalContactNetworkError::Thermal(error) => error.into(),
        ElectricalContactNetworkError::Circuit(error) => error.into(),
        other => RelaxationHeatError::Source(RelaxationHeatSourceError::Contact(other)),
    }
}

impl<'m, 'f, 'n, 'a, F> FrameRelaxationThermalSystem<'m, 'f, 'n, 'a, F> {
    // One cold-source/clock admission for bulk and interface-bearing circuits.
    pub(super) fn admit_electrical_binding(
        &self,
        thermal: &DepositionSimulation<'_>,
        source: &MaterialNetwork<'_>,
        circuit: &MaterialCircuit,
        fields: &[FieldThermalBinding],
        gate: &CancelGate,
    ) -> Result<()> {
        check(gate)?;
        if self.electrical.is_some() || self.mechanical.samples_rendered() != 0
            || thermal.time_s().to_bits() != self.origin_s.to_bits()
            || !core::ptr::eq(source.cells, self.source.cells)
            || !core::ptr::eq(source.links, self.source.links)
        {
            return Err(invalid("electrical relaxation binding requires its original cold clock and same cell/link owners"));
        }
        if fields.is_empty() && !circuit.storage().is_empty() {
            return Err(invalid("source-bound electrical storage requires its original field/thermal work owner"));
        }
        if !fields.is_empty() { field::admit(source, thermal, circuit, fields, gate)?; }
        self.admit_temperatures(thermal, gate)?;
        propose_shared(&self.materials, &self.mechanical, self.source, &self.receivers,
            &self.layer_ranges, thermal, gate)?;
        Ok(())
    }

    /// Bind source-defined interface resistance and BOTH original heat receivers.
    /// The source retains original ordered material sides, interface claims,
    /// thermal-link identity/area and an explicit declared heat partition.
    /// Nothing infers that partition from bulk resistance or conductivity.
    ///
    /// Every thermal trial queries the existing contact network at both proposed
    /// side temperatures. Bulk and interface Joule losses remain separate in
    /// the native circuit report. Their actual heat is deposited once alongside
    /// each layer's independent Maxwell loss. All histories publish together.
    ///
    /// This is fixed electrical geometry and prescribed partition, not changing
    /// contact pressure/area, opening, arcing, contact capacitance, finite-surface
    /// thermomechanics or field-force feedback. Fixed storage and the original
    /// first-order rate/thermal split remain binding. Failed binding is atomic.
    pub fn bind_material_contact_circuit(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        source: &'n ElectricalContactNetwork<'a>,
        circuit: MaterialCircuit,
        gate: &CancelGate,
    ) -> Result<()> {
        self.bind_contact_circuit(thermal, source, circuit, Vec::new(), gate)
    }

    pub(super) fn bind_contact_circuit(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        source: &'n ElectricalContactNetwork<'a>,
        circuit: MaterialCircuit,
        fields: Vec<FieldThermalBinding>,
        gate: &CancelGate,
    ) -> Result<()> {
        self.admit_electrical_binding(thermal, source.bulk, &circuit, &fields, gate)?;
        if source.contacts.is_empty() {
            return Err(invalid("contact relaxation binding needs explicit original interfaces"));
        }
        let prepared = source.prepare(thermal, || gate.is_requested()).map_err(source_error)?;
        let mut initial = prepared.resolve(thermal, || gate.is_requested()).map_err(source_error)?;
        initial.bulk.resistors.sort_unstable_by_key(|resistor| resistor.branch);
        if initial.bulk.resistors.as_slice() != circuit.resistors() {
            return Err(invalid("contact relaxation circuit changed original initial bulk sources or geometry"));
        }
        contacts::admit_sources(&circuit, &initial.contacts, || gate.is_requested())?;
        initial.contacts.sort_unstable_by_key(|contact| contact.id());
        for (actual, expected) in circuit.contacts().iter().zip(&initial.contacts) {
            check(gate)?;
            // Pinning changes the policy receipt, not the selected source.
            // Source identity is admitted above; state/geometry must still agree.
            if actual.id() != expected.id() || actual.branch() != expected.branch()
                || actual.area_m2().to_bits() != expected.area_m2().to_bits()
                || actual.resistance_ohm().to_bits() != expected.resistance_ohm().to_bits()
                || actual.source().resolved().query() != expected.source().resolved().query()
            {
                return Err(invalid("initial electrical interface state/area differs from its original thermal link"));
            }
        }
        let routes = prepared.heat_destinations(thermal).map_err(source_error)?;
        let routes = contacts::admit_routes(thermal, &circuit, &routes, || {
            if gate.is_requested() { Err(DepositionError::Cancelled) } else { Ok(()) }
        })?;
        let destinations = prepared.resistor_destinations().into_iter()
            .map(|route| (route.branch, route.cell)).collect();
        check(gate)?;
        // Original thermal-only source/preparation stays intact. The contact
        // network owns the added bulk/interface source queries for every trial.
        self.electrical = Some(super::Binding {
            circuit, destinations, contacts: Some(Binding { prepared, routes }), fields,
        });
        Ok(())
    }

    /// Original fixed ordered A/B heat receivers, in stable interface-ID order.
    /// Actual interval work is in `report.electrical.contacts`, not these routes.
    pub fn electrical_contact_heat_destinations(&self) -> Vec<ContactHeatCells> {
        self.electrical.as_ref().and_then(|binding| binding.contacts.as_ref())
            .map_or_else(Vec::new, |contact| contact.routes.values().copied().collect())
    }
}
