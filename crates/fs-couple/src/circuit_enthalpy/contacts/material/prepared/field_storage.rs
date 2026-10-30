//! Source-dependent electrical storage and two-sided contact heat, one residual.
use super::*;
use crate::circuit_enthalpy::field_storage::{
    FieldThermalBinding, FieldThermalCandidate, FieldThermalError, PreparedFieldThermal,
    close_report,
};

type Error = CircuitEnthalpyError<FieldThermalError>;

// The original field error already retains the shared material, circuit and
// thermal failures. Preserve each contact cause in that same typed vocabulary.
fn contact_error(error: ElectricalContactNetworkError) -> Error {
    use ElectricalContactNetworkError::*;
    match error {
        Thermal(error) => error.into(),
        Circuit(error) => error.into(),
        Bulk(MaterialNetworkError::Thermal(DepositionError::Cancelled)) =>
            CircuitEnthalpyError::Cancelled,
        Bulk(error) => CircuitEnthalpyError::Source(FieldThermalError::Network(error)),
        State(error) => CircuitEnthalpyError::Source(FieldThermalError::Network(
            MaterialNetworkError::Material(error),
        )),
        Electrical(error) => CircuitEnthalpyError::Source(FieldThermalError::Network(
            MaterialNetworkError::Material(ThermalMaterialError::Material(error)),
        )),
        Invalid(message) => CircuitEnthalpyError::Source(FieldThermalError::Invalid(message)),
    }
}

impl PreparedElectricalContactNetwork<'_, '_> {
    /// Advance the original material C(T)/L(T), bulk resistors, ordered electrical
    /// interfaces and finite caloric bodies on ONE accepted thermal endpoint.
    ///
    /// This retains this prepared network's source pins across calls. Every
    /// nonlinear trial starts from the same accepted charge/flux, rebinds the
    /// original field coefficients at fixed coordinates, then runs the existing
    /// circuit step. The original field thermodynamics supplies signed reversible
    /// heat separately from bulk and two-sided contact Joule heat. All receiving
    /// cells, source conductances and existing boundaries share one Newton solve.
    ///
    /// All C/L branches need whole-volume thermal bindings, including constant
    /// laws. The original positive combined-capacity, coordinate-envelope, source
    /// context and native voltage-constraint checks apply. Interface area, heat
    /// partitions, topology and geometry are fixed; no arc, contact storage,
    /// hysteresis or mechanical force is inferred. Coupling is first order.
    ///
    /// No input is mutated. Publish both returned histories together. The report
    /// retains native per-contact work and the original field entropy/internal
    /// energy account, rather than treating changes in field free energy as heat.
    #[allow(clippy::too_many_arguments)]
    pub fn field_storage_candidate<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        circuit: &MaterialCircuit,
        ports: &[f64],
        bindings: &[FieldThermalBinding],
        controls: CircuitEnthalpySolve,
        guess_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<FieldThermalCandidate<'m>, Error> {
        let cancel = RefCell::new(cancelled);
        let checkpoint = || if (*cancel.borrow_mut())() {
            Err(DepositionError::Cancelled)
        } else { Ok(()) };
        checkpoint()?;
        crate::circuit_enthalpy::admit_driver_unknowns(thermal, controls, guess_j)?;
        let mut initial = self.resolve(thermal, || (*cancel.borrow_mut())())
            .map_err(contact_error)?;
        initial.bulk.resistors.sort_unstable_by_key(|r| r.branch);
        if initial.bulk.resistors.as_slice() != circuit.resistors() {
            return Err(CircuitEnthalpyError::Invalid(
                "field/contact circuit must retain its original bulk sources and geometry",
            ));
        }
        crate::circuit_enthalpy::contacts::admit_sources(
            circuit, &initial.contacts, || (*cancel.borrow_mut())(),
        )?;
        initial.contacts.sort_unstable_by_key(|c| c.id());
        for (actual, expected) in circuit.contacts().iter().zip(&initial.contacts) {
            checkpoint()?;
            // Pinning changes the policy receipt, not the physical source/query.
            if actual.id() != expected.id() || actual.branch() != expected.branch()
                || actual.area_m2().to_bits() != expected.area_m2().to_bits()
                || actual.resistance_ohm().to_bits() != expected.resistance_ohm().to_bits()
                || actual.source().resolved().query() != expected.source().resolved().query()
            {
                return Err(CircuitEnthalpyError::Invalid(
                    "field/contact circuit and original interface state or geometry disagree",
                ));
            }
        }
        let fields = PreparedFieldThermal::new(
            self.network.bulk, thermal, circuit.into(), bindings, &cancel,
        )?;
        let routes = self.heat_destinations(thermal).map_err(contact_error)?;
        let destinations = self.resistor_destinations();
        let dt = target_s - thermal.time_s();
        let (next, electrical, coupled) = crate::circuit_enthalpy::coupled_driver_candidate(
            thermal, target_s,
            CircuitEnthalpyInput { circuit, ports, destinations: &destinations },
            &routes, controls, guess_j,
            |trial, _remaining| {
                checkpoint()?;
                let coefficients = self.resolve(trial, || (*cancel.borrow_mut())())
                    .map_err(contact_error)?;
                let mut bulk = coefficients.bulk;
                let (storage, samples) = fields.refresh(trial, &cancel)?;
                let rebound = circuit.rebind_storage_candidate(
                    storage, ports, || (*cancel.borrow_mut())(),
                )?;
                let candidate = crate::circuit_enthalpy::contacts::advance_circuit(
                    rebound.circuit(), bulk.resistors, coefficients.contacts, ports, dt,
                    controls.circuit, || (*cancel.borrow_mut())(),
                )?;
                let source = fields.account(
                    samples, candidate.circuit().into(), rebound.parameter_work_j(), dt,
                    bulk.source, &mut bulk.sources_w, controls.total_energy_tolerance_j, &cancel,
                )?;
                Ok(DrivenCoefficients {
                    circuit: candidate, sources_w: bulk.sources_w,
                    conductances_w_k: bulk.conductances_w_k, source, circuit_trials: 1,
                })
            },
            || (*cancel.borrow_mut())(),
        )?;
        let work = coupled.electrical.native.supplied;
        let report = close_report(coupled, work, controls.total_energy_tolerance_j)?;
        checkpoint()?;
        Ok(FieldThermalCandidate::new(next, electrical, report))
    }
}
