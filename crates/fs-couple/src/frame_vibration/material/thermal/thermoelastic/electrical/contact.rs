//! Ordered interface resistance and two-sided heat on the shared frame timeline.
use super::*;
use crate::circuit_enthalpy::contacts::{
    self, ContactHeatCells,
    material::{
        ElectricalContactNetwork, ElectricalContactNetworkError,
        PreparedElectricalContactNetwork,
    },
};

pub(super) struct Binding<'n, 'a> {
    pub(super) prepared: PreparedElectricalContactNetwork<'n, 'a>,
    pub(super) routes: BTreeMap<u64, ContactHeatCells>,
}

// Preserve existing typed source failures instead of hiding them in an opaque
// callback error. Cancellation remains distinguishable at the atomic boundary.
pub(super) fn source_error(error: ElectricalContactNetworkError) -> FrameThermalError {
    match error {
        ElectricalContactNetworkError::Invalid(message) => invalid(message),
        ElectricalContactNetworkError::Bulk(error) => network(error),
        ElectricalContactNetworkError::State(error) => {
            network(MaterialNetworkError::Material(error))
        }
        ElectricalContactNetworkError::Electrical(error) => material(error),
        ElectricalContactNetworkError::Circuit(error) => error.into(),
        ElectricalContactNetworkError::Thermal(error) => error.into(),
    }
}

impl<'n, 'a> FrameThermoelasticSystem<'n, 'a> {
    // Bulk-only and contact-bearing circuits share the same original cell,
    // geometry, clock and initial frame-temperature admission.
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
            || !core::ptr::eq(self.source.cells, source.cells)
            || !core::ptr::eq(self.source.links, source.links)
            || source.resistors.is_empty()
            || (fields.is_empty() && !circuit.storage().is_empty())
        {
            return Err(invalid(
                "frame circuit requires one cold binding, original cell/link slices, bulk resistors and explicit field thermal ownership",
            ));
        }
        for (member, state) in self.members.iter().zip(&self.states) {
            check(gate)?;
            if temperatures(&member.region, thermal)?.iter().zip(&state.layers)
                .any(|(&t, layer)| t.to_bits() != layer.temperature_k.to_bits())
            {
                return Err(invalid("frame circuit starts from the jointly accepted material temperatures"));
            }
        }
        if !fields.is_empty() {
            super::field::admit(self, source, thermal, circuit, fields, gate)?;
        }
        Ok(())
    }

    /// Connect the original source-bound electrical interfaces to the running
    /// frame/thermal composition, before its first mechanical sample.
    ///
    /// Each contact retains its ordered interface card, selected claim, original
    /// thermal-link area and A/B receivers. Bulk resistors and every interface
    /// are evaluated together at each proposed thermal endpoint. The native
    /// circuit's separately integrated contact heat enters both finite cells
    /// through the existing contact-deposition owner, exactly once. A receiving
    /// frame then refreshes its own E(T), expansion and viscosity normally.
    ///
    /// The existing prepared contact network owns source continuity and fixed
    /// interface geometry. Caller-specified heat fractions remain explicit model
    /// inputs, not an inferred effusivity or contact-mechanics solution. Opening,
    /// arcing, contact capacitance, moving area and changing C/L are not modeled.
    /// Mechanical motion and material temperature retain the documented first-
    /// order split. This adds no thermal integrator or contact-resistance law.
    pub fn with_material_contact_circuit(
        self,
        thermal: &DepositionSimulation<'_>,
        source: &'n ElectricalContactNetwork<'a>,
        circuit: MaterialCircuit,
        gate: &CancelGate,
    ) -> Result<Self> {
        self.bind_material_contact_circuit(thermal, source, circuit, Vec::new(), gate)
    }

    pub(super) fn bind_material_contact_circuit(
        mut self,
        thermal: &DepositionSimulation<'_>,
        source: &'n ElectricalContactNetwork<'a>,
        circuit: MaterialCircuit,
        fields: Vec<FieldThermalBinding>,
        gate: &CancelGate,
    ) -> Result<Self> {
        self.admit_electrical_binding(thermal, source.bulk, &circuit, &fields, gate)?;
        if source.contacts.is_empty() {
            return Err(invalid("contact frame binding requires explicit original interfaces"));
        }
        let prepared = source.prepare(thermal, || gate.is_requested()).map_err(source_error)?;
        let mut initial = prepared.resolve(thermal, || gate.is_requested()).map_err(source_error)?;
        initial.bulk.resistors.sort_unstable_by_key(|resistor| resistor.branch);
        if initial.bulk.resistors.as_slice() != circuit.resistors() {
            return Err(invalid("initial contact circuit changed its original bulk sources or geometry"));
        }
        contacts::admit_sources(&circuit, &initial.contacts, || gate.is_requested())?;
        initial.contacts.sort_unstable_by_key(|contact| contact.id());
        for (actual, expected) in circuit.contacts().iter().zip(&initial.contacts) {
            check(gate)?;
            // Preparation replaces the policy receipt with an explicit pin.
            // The physical query, actual selected source, reduction and state
            // must match; policy receipt equality would reject a valid binding.
            if actual.id() != expected.id() || actual.branch() != expected.branch()
                || actual.area_m2().to_bits() != expected.area_m2().to_bits()
                || actual.resistance_ohm().to_bits() != expected.resistance_ohm().to_bits()
                || actual.source().resolved().query() != expected.source().resolved().query()
            {
                return Err(invalid("initial interface geometry or physical state differs from its thermal link"));
            }
        }
        let routes = prepared.heat_destinations(thermal).map_err(source_error)?;
        let routes = contacts::admit_routes(thermal, &circuit, &routes, || {
            if gate.is_requested() { Err(DepositionError::Cancelled) } else { Ok(()) }
        })?;
        let destinations = prepared.resistor_destinations().into_iter()
            .map(|route| (route.branch, route.cell)).collect();
        check(gate)?;
        self.electrical = Some(super::Binding {
            circuit, destinations, contacts: Some(Binding { prepared, routes }), fields,
        });
        Ok(self)
    }

    /// Original fixed A/B heat partitions, in stable interface-ID order.
    /// Exact source contact work remains in the electrical step's contact rows.
    pub fn electrical_contact_heat_destinations(&self) -> Vec<ContactHeatCells> {
        self.electrical.as_ref().and_then(|binding| binding.contacts.as_ref())
            .map_or_else(Vec::new, |binding| binding.routes.values().copied().collect())
    }
}
