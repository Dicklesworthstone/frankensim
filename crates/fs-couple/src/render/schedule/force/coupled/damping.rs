//! Live resistance changes on the original complete modal propagator.
//!
//! A dashpot stores no energy. This compatibility API retains the original
//! spring/rest parameters and delegates to the SAME coefficient transaction.
//! Source applicability and the finite thermal update belong to the caller.

use super::coefficients::rates::mobility::ModalDampingRate;
use super::coefficients::{CoefficientEndpoint, PreparedConnectionCoefficients};
use super::*;

/// One replacement resistance in the original physical connection order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConnectionDamping {
    pub connection: usize,
    /// New nonnegative resistance [N s/m] on the unchanged work coordinate.
    pub damping_n_s_m: f64,
}

/// Complete replacement of one original component's intrinsic viscous rate.
/// Connection, Maxwell and radiation losses remain separate and unchanged.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ComponentDampingRate {
    /// Original component index, in construction order.
    pub component: usize,
    /// Coordinate within that component's unchanged mass-normalized basis.
    pub mode: usize,
    /// Nonnegative coefficient `c` [1/s] in `q'' + c*q' + omega^2*q = F`.
    /// For dimensional resistance `b` on mass `m`, this is `b/m`.
    pub damping_rate_per_s: f64,
}

pub(crate) type PreparedConnectionDamping = PreparedConnectionCoefficients;

impl CoupledModalSystem {
    /// Replace selected intrinsic rates at the accepted mechanical boundary.
    ///
    /// Entries must be unique and strictly ordered by `(component, mode)`;
    /// rates are finite and nonnegative, including on free coordinates. Each
    /// value replaces the complete intrinsic rate of that coordinate, so its
    /// caller must retain every intended material/ambient contribution. Use
    /// the source-bound thermal owner for selective viscosity feedback.
    ///
    /// The existing material-mobility transaction rebuilds the SAME complete
    /// ordinary/prepared/harmonic-only operator. All component states, Maxwell
    /// and acoustic histories, springs, connections, sample clock and previous
    /// interval reports survive. An ideal viscous rate has no stored energy:
    /// this instantaneous change supplies no impulse, parameter work or heat.
    /// Subsequent intervals use the new rate for both motion and physical heat.
    ///
    /// Empty/unchanged batches allocate and rebuild nothing. Invalid input,
    /// cancellation or failed operator preparation preserves the original
    /// assembly. The caller owns the rate law, its physical domain and temporal
    /// resolution; split prescribed changes at accepted sample boundaries.
    pub fn set_component_damping_rates(
        &mut self,
        updates: &[ComponentDampingRate],
        gate: &CancelGate,
    ) -> Result<(), ModalCouplingError> {
        poll(Some(gate))?;
        if updates.len() > self.mode_count() {
            return Err(invalid(
                "intrinsic damping updates exceed original mode count",
            ));
        }
        let mut previous = None;
        let mut changed = false;
        for update in updates {
            poll(Some(gate))?;
            let original = self
                .models
                .get(update.component)
                .and_then(|model| model.damping_rates_per_s().get(update.mode))
                .ok_or_else(|| {
                    invalid("intrinsic damping update names no original component coordinate")
                })?;
            let key = (update.component, update.mode);
            if previous.is_some_and(|p| key <= p)
                || !update.damping_rate_per_s.is_finite()
                || update.damping_rate_per_s < 0.0
            {
                return Err(invalid(
                    "intrinsic damping needs ordered unique coordinates and nonnegative finite rates",
                ));
            }
            previous = Some(key);
            changed |= *original != update.damping_rate_per_s;
        }
        if !changed {
            poll(Some(gate))?;
            return Ok(());
        }
        let damping: Vec<_> = updates
            .iter()
            .map(|update| ModalDampingRate {
                mode: self.offsets[update.component] + update.mode,
                previous: self.models[update.component].damping_rates_per_s()[update.mode],
                next: update.damping_rate_per_s,
            })
            .collect();
        let prepared =
            self.prepare_material_mobility(&[], &damping, CoefficientEndpoint::Accepted, gate)?;
        poll(Some(gate))?;
        self.install_material_relaxation_times(prepared);
        Ok(())
    }

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
