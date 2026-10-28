//! Live resistance changes on the original complete modal propagator.
//!
//! A dashpot stores no energy. This compatibility API retains the original
//! spring/rest parameters and delegates to the SAME coefficient transaction.
//! Source applicability and the finite thermal update belong to the caller.

use super::coefficients::{CoefficientEndpoint, PreparedConnectionCoefficients};
use super::*;

/// One replacement resistance in the original physical connection order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConnectionDamping {
    pub connection: usize,
    /// New nonnegative resistance [N s/m] on the unchanged work coordinate.
    pub damping_n_s_m: f64,
}

pub(crate) type PreparedConnectionDamping = PreparedConnectionCoefficients;

impl CoupledModalSystem {
    /// Current resistance of an original physical connection, if it exists.
    pub fn connection_damping_n_s_m(&self, connection: usize) -> Option<f64> {
        self.connections
            .get(connection)
            .map(|law| law.damping_n_s_m)
    }

    /// Replace selected original resistances without resetting any history.
    ///
    /// Indices are strictly increasing and unique; coefficients are finite and
    /// nonnegative, with positive resistance on existing constant-loss anchors.
    /// All spring/rest coefficients, q/v, clock, Maxwell/acoustic states, ports
    /// and previous accepted reports remain unchanged. The original complete
    /// ordinary/prepared propagator and its budgets are retained. Unchanged
    /// batches allocate and rebuild nothing; refusal is atomic.
    ///
    /// An ideal dashpot has no parameter storage/work. Its subsequent physical
    /// dissipation must still be assigned by the owning material/thermal law.
    pub fn set_connection_damping(
        &mut self,
        updates: &[ConnectionDamping],
        gate: &CancelGate,
    ) -> Result<(), ModalCouplingError> {
        let prepared = self.prepare_damping(updates, CoefficientEndpoint::Accepted, gate)?;
        poll(Some(gate))?;
        self.install_connection_coefficients(prepared);
        Ok(())
    }

    pub(crate) fn prepare_staged_connection_damping(
        &self,
        updates: &[ConnectionDamping],
        gate: &CancelGate,
    ) -> Result<PreparedConnectionDamping, ModalCouplingError> {
        self.prepare_damping(updates, CoefficientEndpoint::Staged, gate)
    }

    fn prepare_damping(
        &self,
        updates: &[ConnectionDamping],
        endpoint: CoefficientEndpoint,
        gate: &CancelGate,
    ) -> Result<PreparedConnectionDamping, ModalCouplingError> {
        self.prepare_coefficients_iter(
            updates.iter().map(|update| {
                let mut law = self
                    .connection_coefficients(update.connection)
                    .ok_or_else(|| invalid("resistance update names an absent connection"))?;
                law.damping_n_s_m = update.damping_n_s_m;
                Ok(law)
            }),
            endpoint,
            gate,
        )
        .map_err(|error| match error {
            ModalCouplingError::Budget {
                what: "connection reaction at coefficient update",
                value,
                limit,
            } => ModalCouplingError::Budget {
                what: "connection reaction at resistance update",
                value,
                limit,
            },
            other => other,
        })
    }

    pub(crate) fn install_connection_damping(&mut self, prepared: PreparedConnectionDamping) {
        self.install_connection_coefficients(prepared);
    }
}
