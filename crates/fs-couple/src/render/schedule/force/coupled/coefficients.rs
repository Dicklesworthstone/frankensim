//! Material-dependent spring, rest-position and resistance changes on the SAME
//! coupled propagator. Fixed-coordinate spring-energy change is explicit
//! parameter work; a surrounding material/thermal owner supplies its physical
//! internal-energy and entropy interpretation. No impulse or heat is invented.

use super::*;
use fs_math::det;

/// Complete replacement of one original physical connection law.
/// Its work-conjugate column and source port identity remain unchanged.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConnectionCoefficients {
    pub connection: usize,
    pub stiffness_n_m: f64,
    pub damping_n_s_m: f64,
    pub rest_extension_m: f64,
}

/// Spring energy at the SAME mechanical coordinates before/after a law change.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConnectionParameterWork {
    pub connection: usize,
    pub free_energy_before_j: f64,
    pub free_energy_after_j: f64,
    /// Positive means energy supplied to mechanical spring storage [J].
    pub parameter_work_j: f64,
}

/// Instantaneous coefficient work, separate from every accepted time interval.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ConnectionCoefficientReport {
    pub parameter_work_j: f64,
    /// One row per requested index, in original order. An entirely unchanged
    /// batch is an allocation-free no-op and returns no rows.
    pub changes: Vec<ConnectionParameterWork>,
}

#[derive(Clone, Copy)]
pub(super) enum CoefficientEndpoint {
    Accepted,
    Staged,
}

struct Replacement {
    connections: Vec<ConnectionLaw>,
    exact: Option<exact::ExactNetwork>,
    mechanical_bound: f64,
    frequency_bound: f64,
}

/// Prepared only for the unchanged owner and mechanical endpoint that admitted
/// it. The surrounding native transaction installs it by infallible moves.
pub(crate) struct PreparedConnectionCoefficients {
    replacement: Option<Replacement>,
    report: ConnectionCoefficientReport,
}
impl PreparedConnectionCoefficients {
    pub(crate) fn parameter_work_j(&self) -> f64 {
        self.report.parameter_work_j
    }
    pub(crate) fn parameter_changes(&self) -> &[ConnectionParameterWork] {
        &self.report.changes
    }
}

impl CoupledModalSystem {
    pub fn connection_coefficients(&self, connection: usize) -> Option<ConnectionCoefficients> {
        self.connections
            .get(connection)
            .map(|law| ConnectionCoefficients {
                connection,
                stiffness_n_m: law.stiffness_n_m,
                damping_n_s_m: law.damping_n_s_m,
                rest_extension_m: law.rest_extension_m,
            })
    }

    /// Atomically change original spring/rest/dashpot laws at fixed q and v.
    ///
    /// Indices must be distinct and strictly increasing. Stiffness/resistance
    /// are finite and nonnegative; rest positions are finite. All original
    /// coordinates, clock, material memory, acoustic state, physical columns,
    /// ports and previous frame/transfer ledger survive. The same complete
    /// ordinary/prepared operator is rebuilt, retaining every extension and its
    /// original budgets. Actual new spring storage, complete instantaneous
    /// stiffness/Nyquist admission and ALL endpoint reactions are checked.
    ///
    /// The returned fixed-coordinate spring-energy change belongs to the
    /// parameter source. Changing a rest position may supply or extract energy
    /// even without motion. A thermal consumer must close physical U and S
    /// using its original constitutive potential; this API does not call that
    /// energy heat, dissipation or time-integrated mechanical work.
    /// The last frame continues to describe the preceding interval.
    pub fn set_connection_coefficients(
        &mut self,
        updates: &[ConnectionCoefficients],
        gate: &CancelGate,
    ) -> Result<ConnectionCoefficientReport, ModalCouplingError> {
        let prepared = self.prepare_coefficients_iter(
            updates.iter().copied().map(Ok),
            CoefficientEndpoint::Accepted,
            gate,
        )?;
        let report = prepared.report.clone();
        poll(Some(gate))?;
        self.install_connection_coefficients(prepared);
        Ok(report)
    }

    pub(crate) fn prepare_staged_connection_coefficients(
        &self,
        updates: &[ConnectionCoefficients],
        gate: &CancelGate,
    ) -> Result<PreparedConnectionCoefficients, ModalCouplingError> {
        self.prepare_coefficients_iter(
            updates.iter().copied().map(Ok),
            CoefficientEndpoint::Staged,
            gate,
        )
    }

    // A cloneable iterator keeps the original damping-only API's no-op path
    // allocation-free while sharing every preparation/publication rule.
    pub(super) fn prepare_coefficients_iter<I>(
        &self,
        updates: I,
        endpoint: CoefficientEndpoint,
        gate: &CancelGate,
    ) -> Result<PreparedConnectionCoefficients, ModalCouplingError>
    where
        I: ExactSizeIterator<Item = Result<ConnectionCoefficients, ModalCouplingError>> + Clone,
    {
        poll(Some(gate))?;
        if updates.len() > self.connections.len() {
            return Err(invalid(
                "coefficient updates exceed original connection count",
            ));
        }
        let mut previous = None;
        let mut changed = false;
        let mut stiffness_changed = false;
        for update in updates.clone() {
            poll(Some(gate))?;
            let update = update?;
            let old = self
                .connections
                .get(update.connection)
                .ok_or_else(|| invalid("coefficient update names an absent connection"))?;
            if previous.is_some_and(|index| update.connection <= index)
                || !update.stiffness_n_m.is_finite()
                || update.stiffness_n_m < 0.0
                || !update.damping_n_s_m.is_finite()
                || update.damping_n_s_m < 0.0
                || !update.rest_extension_m.is_finite()
                || (old.constant_loss_reference_omega_rad_s.is_some()
                    && update.damping_n_s_m == 0.0)
            {
                return Err(invalid(
                    "ordered unique finite nonnegative K/C and finite rest required; constant-loss resistance stays positive",
                ));
            }
            previous = Some(update.connection);
            stiffness_changed |= update.stiffness_n_m != old.stiffness_n_m;
            changed |= update.stiffness_n_m != old.stiffness_n_m
                || update.damping_n_s_m != old.damping_n_s_m
                || update.rest_extension_m != old.rest_extension_m;
        }
        if !changed {
            return Ok(PreparedConnectionCoefficients {
                replacement: None,
                report: ConnectionCoefficientReport::default(),
            });
        }
        let models = match endpoint {
            CoefficientEndpoint::Accepted => &self.models,
            CoefficientEndpoint::Staged => &self.candidates,
        };
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
        let mut report = ConnectionCoefficientReport {
            parameter_work_j: 0.0,
            changes: Vec::with_capacity(updates.len()),
        };
        for update in updates {
            poll(Some(gate))?;
            let update = update?;
            let old = &mut connections[update.connection];
            let column = &self.columns[update.connection];
            let before = extension(models, column, old.rest_extension_m)?;
            let after = extension(models, column, update.rest_extension_m)?;
            let free_energy_before_j = finite(0.5 * old.stiffness_n_m * before * before)?;
            let free_energy_after_j = finite(0.5 * update.stiffness_n_m * after * after)?;
            let parameter_work_j = finite(free_energy_after_j - free_energy_before_j)?;
            report.parameter_work_j = finite(report.parameter_work_j + parameter_work_j)?;
            report.changes.push(ConnectionParameterWork {
                connection: update.connection,
                free_energy_before_j,
                free_energy_after_j,
                parameter_work_j,
            });
            old.stiffness_n_m = update.stiffness_n_m;
            old.damping_n_s_m = update.damping_n_s_m;
            old.rest_extension_m = update.rest_extension_m;
        }
        let old_mechanical_bound = self
            .radiation
            .as_ref()
            .map_or(self.linear_frequency_upper_bound_rad_s, |r| {
                r.mechanical_bound
            });
        let mechanical_bound = if stiffness_changed {
            self.refreshed_mechanical_bound(&connections, gate)?
        } else {
            old_mechanical_bound
        };
        let frequency_bound = self.radiation.as_ref().map_or(Ok(mechanical_bound), |r| {
            r.frequency_bound(mechanical_bound)
        })?;
        if frequency_bound > core::f64::consts::PI / self.dt * self.config.nyquist_guard_fraction {
            return Err(invalid(
                "refreshed complete stiffness exceeds the Nyquist guard",
            ));
        }
        let (material_energy, radiation_energy) = match endpoint {
            CoefficientEndpoint::Accepted => {
                (self.material_energy_j()?, self.radiation_energy_j()?)
            }
            CoefficientEndpoint::Staged => (
                self.staged_frame.material_energy_j,
                self.staged_frame.radiation_energy_j,
            ),
        };
        let mut energy = component_energy(models)?;
        for (column, law) in self.columns.iter().zip(&connections) {
            let x = extension(models, column, law.rest_extension_m)?;
            energy = finite(energy + 0.5 * law.stiffness_n_m * x * x)?;
        }
        energy = finite(energy + material_energy + radiation_energy)?;
        limit(
            "total storage at coefficient update",
            energy,
            self.config.maximum_total_energy_j,
        )?;
        let mut refreshed_radiation = stiffness_changed.then(|| self.radiation.clone()).flatten();
        if let Some(r) = &mut refreshed_radiation {
            r.mechanical_bound = mechanical_bound;
        }
        let radiation = refreshed_radiation.as_ref().or(self.radiation.as_ref());
        let exact = if let Some(budget) = self.prepared_budget {
            Some(exact::ExactNetwork::new_prepared_with_relaxation(
                &self.models,
                &self.columns,
                &connections,
                self.mode_count(),
                self.dt,
                budget,
                &self.relaxation_branches,
                radiation,
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
                radiation,
                gate,
            )?)
        } else {
            None
        };
        // A changed stiffness can also change a different constant-loss port's
        // projected damping column. Admit ALL actual new reactions.
        for (index, law) in connections.iter().enumerate() {
            poll(Some(gate))?;
            let column = &self.columns[index];
            let damping_column = exact
                .as_ref()
                .and_then(|e| e.damping_column(index))
                .unwrap_or(column);
            let velocity = models
                .iter()
                .flat_map(|model| model.states())
                .zip(damping_column)
                .try_fold(0.0, |sum, (state, b)| {
                    finite(sum + b * state.velocity_m_sqrt_kg_per_s)
                })?;
            let force = finite(
                -law.stiffness_n_m * extension(models, column, law.rest_extension_m)?
                    - law.damping_n_s_m * velocity,
            )?;
            limit(
                "connection reaction at coefficient update",
                force.abs(),
                self.config.maximum_abs_connection_force_n,
            )?;
        }
        poll(Some(gate))?;
        Ok(PreparedConnectionCoefficients {
            replacement: Some(Replacement {
                connections,
                exact,
                mechanical_bound,
                frequency_bound,
            }),
            report,
        })
    }

    fn refreshed_mechanical_bound(
        &self,
        connections: &[ConnectionLaw],
        gate: &CancelGate,
    ) -> Result<f64, ModalCouplingError> {
        let mut squared = 0.0_f64;
        for model in &self.models {
            poll(Some(gate))?;
            for mode in model.modes() {
                squared = squared.max(finite(mode.angular_frequency_rad_s.powi(2))?);
            }
        }
        for (column, law) in self.columns.iter().zip(connections) {
            poll(Some(gate))?;
            let norm = column.iter().try_fold(0.0, |sum, b| finite(sum + b * b))?;
            squared = finite(squared + law.stiffness_n_m * norm)?;
        }
        // Exactly the original instantaneous Maxwell-stiffness contribution;
        // relaxation memory itself is preserved, never relaxed or reseeded.
        for branch in &self.relaxation_branches {
            poll(Some(gate))?;
            let root = det::sqrt(branch.stiffness);
            for b in branch.projection.iter().step_by(2) {
                let value = finite(root * b)?;
                squared = finite(squared + value * value)?;
            }
        }
        finite(det::sqrt(squared))
    }

    pub(crate) fn install_connection_coefficients(
        &mut self,
        prepared: PreparedConnectionCoefficients,
    ) {
        if let Some(replacement) = prepared.replacement {
            self.connections = replacement.connections;
            self.exact = replacement.exact;
            self.linear_frequency_upper_bound_rad_s = replacement.frequency_bound;
            if let Some(r) = &mut self.radiation {
                r.mechanical_bound = replacement.mechanical_bound;
            }
        }
    }
}
