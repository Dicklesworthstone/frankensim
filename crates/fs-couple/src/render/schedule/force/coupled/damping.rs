//! Live resistance changes on the original complete modal propagator.
//!
//! A dashpot stores no energy. Changing its nonnegative coefficient at fixed
//! coordinates changes the subsequent force and loss, with no parameter-work
//! impulse. Spring, Maxwell and acoustic storage remain with their original
//! owners. Source applicability and the temperature update belong to the caller.

use super::*;

/// One replacement resistance in the original physical connection order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConnectionDamping {
    pub connection: usize,
    /// New nonnegative resistance [N s/m] on the unchanged work coordinate.
    pub damping_n_s_m: f64,
}

/// Prepared only against the owner whose unchanged topology will install it.
/// The private value permits a surrounding thermal transaction to refuse before
/// publishing either its mechanical step or its new resistance coefficients.
pub(crate) struct PreparedConnectionDamping {
    replacement: Option<(Vec<ConnectionLaw>, Option<exact::ExactNetwork>)>,
}

impl CoupledModalSystem {
    /// Current resistance of an original physical connection, if it exists.
    pub fn connection_damping_n_s_m(&self, connection: usize) -> Option<f64> {
        self.connections
            .get(connection)
            .map(|law| law.damping_n_s_m)
    }

    /// Replace selected physical resistances without resetting running history.
    ///
    /// Updates must have distinct, strictly increasing original indices. All
    /// modal coordinates, Maxwell internal strains, acoustic states, sample
    /// clock, springs, attachments, source ports and previous accepted reports
    /// remain unchanged. The SAME ordinary/prepared propagator is rebuilt with
    /// every retained extension. Connected islands may merge or split, but no
    /// physical coordinate is removed. Existing preparation budgets still apply.
    /// Equal coefficients (including an empty list) need no new allocation or
    /// matrix preparation. A failed rebuild or cancellation changes no state.
    ///
    /// No energy or entropy is assigned to changing an ideal dashpot itself.
    /// Temperature-dependent material applications must provide the original
    /// source law and account for the subsequent dissipated work as heat.
    /// This does not change elastic/material-memory parameters or radiation.
    /// The last frame/transfer ledger still describes the preceding interval;
    /// the new resistance acts from this boundary forward. The ordinary next
    /// step retains its full-interval force and energy gates.
    pub fn set_connection_damping(
        &mut self,
        updates: &[ConnectionDamping],
        gate: &CancelGate,
    ) -> Result<(), ModalCouplingError> {
        let prepared = self.prepare_connection_damping(updates, &self.models, gate)?;
        poll(Some(gate))?;
        self.install_connection_damping(prepared);
        Ok(())
    }

    /// For a native stage/admit/publish callback: the refreshed instantaneous
    /// reaction gate uses the actual staged endpoint, not the previous state.
    pub(crate) fn prepare_staged_connection_damping(
        &self,
        updates: &[ConnectionDamping],
        gate: &CancelGate,
    ) -> Result<PreparedConnectionDamping, ModalCouplingError> {
        self.prepare_connection_damping(updates, &self.candidates, gate)
    }

    fn prepare_connection_damping(
        &self,
        updates: &[ConnectionDamping],
        endpoint: &[ModalAcousticTimeModel],
        gate: &CancelGate,
    ) -> Result<PreparedConnectionDamping, ModalCouplingError> {
        poll(Some(gate))?;
        if updates.len() > self.connections.len() {
            return Err(invalid(
                "resistance updates exceed original connection count",
            ));
        }
        let mut previous = None;
        let mut changed = false;
        for update in updates {
            poll(Some(gate))?;
            let law = self
                .connections
                .get(update.connection)
                .ok_or_else(|| invalid("resistance update names an absent connection"))?;
            if previous.is_some_and(|index| update.connection <= index)
                || !update.damping_n_s_m.is_finite()
                || update.damping_n_s_m < 0.0
                || (law.constant_loss_reference_omega_rad_s.is_some()
                    && update.damping_n_s_m == 0.0)
            {
                return Err(invalid(
                    "ordered unique finite nonnegative resistances required; constant-loss anchors stay positive",
                ));
            }
            previous = Some(update.connection);
            changed |= update.damping_n_s_m != law.damping_n_s_m;
        }
        if !changed {
            return Ok(PreparedConnectionDamping { replacement: None });
        }
        let mut connections: Vec<_> = self
            .connections
            .iter()
            .map(|law| ConnectionLaw {
                stiffness_n_m: law.stiffness_n_m,
                damping_n_s_m: law.damping_n_s_m,
                rest_extension_m: law.rest_extension_m,
                constant_loss_reference_omega_rad_s: law.constant_loss_reference_omega_rad_s,
            })
            .collect();
        for update in updates {
            connections[update.connection].damping_n_s_m = update.damping_n_s_m;
        }
        let exact = if let Some(budget) = self.prepared_budget {
            Some(exact::ExactNetwork::new_prepared(
                &self.models,
                &self.columns,
                &connections,
                self.mode_count(),
                self.dt,
                budget,
                self.radiation.as_ref(),
                gate,
            )?)
        } else if connections
            .iter()
            .any(|law| law.stiffness_n_m != 0.0 || law.damping_n_s_m != 0.0)
            || !self.relaxation_branches.is_empty()
            || self.radiation.is_some()
        {
            Some(exact::ExactNetwork::new_with_extensions(
                &self.models,
                &self.columns,
                &connections,
                self.mode_count(),
                &self.relaxation_branches,
                self.radiation.as_ref(),
                gate,
            )?)
        } else {
            None
        };
        for update in updates {
            poll(Some(gate))?;
            let index = update.connection;
            let law = &connections[index];
            let column = &self.columns[index];
            let dashpot = exact
                .as_ref()
                .and_then(|operator| operator.damping_column(index))
                .unwrap_or(column);
            let velocity = endpoint
                .iter()
                .flat_map(|model| model.states())
                .zip(dashpot)
                .try_fold(0.0, |sum, (state, weight)| {
                    finite(sum + weight * state.velocity_m_sqrt_kg_per_s)
                })?;
            let reaction = finite(
                -law.stiffness_n_m * extension(endpoint, column, law.rest_extension_m)?
                    - law.damping_n_s_m * velocity,
            )?;
            limit(
                "connection reaction at resistance update",
                reaction.abs(),
                self.config.maximum_abs_connection_force_n,
            )?;
        }
        poll(Some(gate))?;
        Ok(PreparedConnectionDamping {
            replacement: Some((connections, exact)),
        })
    }

    /// Infallible publication after the enclosing complete physical transaction.
    pub(crate) fn install_connection_damping(&mut self, prepared: PreparedConnectionDamping) {
        if let Some((connections, exact)) = prepared.replacement {
            self.connections = connections;
            self.exact = exact;
        }
    }
}
