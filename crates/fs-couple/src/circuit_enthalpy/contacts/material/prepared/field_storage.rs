//! Source-dependent electrical storage and two-sided contact heat, one residual.
use super::*;
use crate::circuit_enthalpy::field_storage::{
    FieldThermalBinding, FieldThermalCandidate, FieldThermalError, driven_candidate,
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
        let routes = self.heat_destinations(thermal).map_err(contact_error)?;
        let destinations = self.resistor_destinations();
        driven_candidate(
            self.network.bulk, thermal, target_s, circuit, ports, bindings,
            &destinations, &routes, controls, guess_j,
            |trial| self.resolve(trial, || (*cancel.borrow_mut())()).map_err(contact_error),
            |_| 0.0,
            || (*cancel.borrow_mut())(),
        )
    }
}
