//! Kirchhoff storage and frame materials on the original shared heat residual.
use super::*;
use crate::material_circuit::{
    MaterialCircuit, MaterialCircuitCandidate, MaterialCircuitSolve, MaterialCircuitStep,
};
use fs_material::conductor::ElectricalPropertyBasis;

/// Held electrical source ports and independent native circuit tolerances.
/// Use an empty port slice for a circuit powered only by its finite C/L stores.
#[derive(Clone, Copy)]
pub struct FrameElectricalInput<'p> {
    pub ports: &'p [f64],
    pub solve: MaterialCircuitSolve,
}

pub(super) struct Binding {
    pub(super) circuit: MaterialCircuit,
    destinations: BTreeMap<usize, CellId>,
}

impl<'n, 'a> FrameThermoelasticSystem<'n, 'a> {
    /// Bind a circuit to the SAME material cells and conduction links as the
    /// complete thermoelastic assembly. `source` adds whole-cell resistors to
    /// the original network; its `cells` and `links` must borrow the original
    /// slices, not independent copies or a second set of material laws.
    ///
    /// Every circuit resistor must equal its freshly resolved initial source
    /// and geometry. Conductivity/resistivity bases are explicit per branch.
    /// The one prepared network replaces the thermal-only preparation, so each
    /// residual trial resolves mechanical temperature, R(T), and conduction at
    /// one endpoint. No additional material evaluator or thermal solve is used.
    ///
    /// The original charges and fluxes are retained. This is cold binding only:
    /// the frame sample clock must still be zero, and thermal cells must match
    /// the accepted frame state. Electrical contacts and material-derived C/L
    /// currently require their own reciprocal thermal laws and refuse here;
    /// explicitly authored fixed C/L retain their native storage accounting.
    /// Fixed geometry means no motional EMF, magnetic force, or changing wire
    /// length is implied. The mechanical/thermal split remains first order.
    #[allow(clippy::too_many_arguments)]
    pub fn with_material_circuit(
        mut self,
        thermal: &DepositionSimulation<'_>,
        source: &'n MaterialNetwork<'a>,
        circuit: MaterialCircuit,
        electrical_bases: &[(usize, ElectricalPropertyBasis)],
        gate: &CancelGate,
    ) -> Result<Self> {
        check(gate)?;
        if self.electrical.is_some() || self.mechanical.samples_rendered() != 0
            || thermal.time_s().to_bits() != self.origin_s.to_bits()
            || !core::ptr::eq(self.source.cells, source.cells)
            || !core::ptr::eq(self.source.links, source.links)
            || source.resistors.is_empty()
            || !circuit.contacts().is_empty() || !circuit.storage().is_empty()
        {
            return Err(invalid(
                "frame circuit requires one cold binding, original cell/link slices, bulk resistors and fixed authored C/L",
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
        let prepared = source.prepare_with_bases(thermal, electrical_bases,
            || gate.is_requested()).map_err(network)?;
        let mut resolved = prepared.resolve(thermal, || gate.is_requested()).map_err(network)?.resistors;
        resolved.sort_unstable_by_key(|resistor| resistor.branch);
        if resolved.as_slice() != circuit.resistors() {
            return Err(invalid("initial circuit must retain every original material resistor and its geometry"));
        }
        let destinations = source.resistors.iter().map(|r|
            (r.destination.branch, r.destination.cell)).collect();
        check(gate)?;
        self.prepared = prepared;
        self.electrical = Some(Binding { circuit, destinations });
        Ok(self)
    }

    /// Accepted native charge, flux, topology, and source material snapshots.
    /// No mutable circuit reference can advance it independently of frame heat.
    pub fn electrical(&self) -> Option<&MaterialCircuit> {
        self.electrical.as_ref().map(|binding| &binding.circuit)
    }

    /// Advance mechanics, circuit storage, and all finite material cells as one
    /// transaction. Every electrical trial starts from the same accepted charge
    /// and flux. Resistance and conduction use endpoint temperature; mechanical
    /// coefficients use the accepted start temperature over this interval.
    /// New source stiffness/expansion/viscosity act in the next interval.
    /// Any electrical, thermal, material, mechanical, or cancellation failure
    /// publishes none of the three histories. Sources are held over the sample.
    #[allow(clippy::too_many_arguments)]
    pub fn step_electrical(
        &mut self,
        thermal: &mut DepositionSimulation<'_>,
        external_modal_forces: &[f64],
        input: FrameElectricalInput<'_>,
        controls: FrameThermalSolve,
        energy_guess_j: &[f64],
        gate: &CancelGate,
    ) -> Result<FrameThermoelasticReport> {
        self.step_inner(thermal, external_modal_forces, Some(input), controls, energy_guess_j, gate)
    }
}

impl Binding {
    pub(super) fn deposit(
        &self,
        candidate: &MaterialCircuitCandidate,
        dt: f64,
        powers: &mut BTreeMap<CellId, f64>,
        checkpoint: &mut dyn FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<()> {
        for branch in &candidate.step().resistors {
            checkpoint()?;
            let cell = self.destinations.get(&branch.material.branch)
                .ok_or_else(|| invalid("electrical loss has no original frame-network receiver"))?;
            let power = finite(branch.heat_j / dt)?;
            if branch.heat_j > 0.0 && power == 0.0 {
                return Err(invalid("frame electrical heat power underflow"));
            }
            let total = powers.entry(*cell).or_insert(0.0);
            *total = finite(*total + power)?;
        }
        Ok(())
    }
}

// Add admitted native driver accounts, not synthetic heat rates. The shared
// caloric owner still independently compares integrated source powers with
// total physical dissipation and closes all supplied work and storage.
pub(super) struct Driver {
    mechanical: MechanicalCandidate,
    pub(super) electrical: Option<MaterialCircuitCandidate>,
}
impl Driver {
    pub(super) fn new(
        mechanical: MechanicalCandidate,
        electrical: Option<MaterialCircuitCandidate>,
    ) -> Result<Self> {
        let driver = Self { mechanical, electrical };
        let a = driver.account();
        finite(a.heat_j)?;
        finite(a.stored_change_j)?;
        finite(a.supplied_j)?;
        Ok(driver)
    }
}
impl ThermalDriver for Driver {
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
        }
        account
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        let (_, mechanical) = self.mechanical.into_parts();
        match self.electrical {
            Some(candidate) => {
                let (circuit, electrical) = candidate.into_parts();
                (Some(circuit), (mechanical, Some(electrical)))
            }
            None => (None, (mechanical, None)),
        }
    }
}
