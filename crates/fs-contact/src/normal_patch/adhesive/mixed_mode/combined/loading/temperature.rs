//! Original thermal sources on the held proportional I/II/III loading family.
//!
//! The accepted approach, vector shear, twist and load factor remain fixed.
//! Attached continuation recomputes the first boundary and the actual relaxed
//! common radius; a new limit cannot silently clamp or detach accepted history.
//! Detached continuation needs no hypothetical attached preload and retains
//! every joule already released. All four signed storage changes are explicit
//! coefficient work, not caloric heat, dissipation or a mechanical impulse.
//!
//! This endpoint operation makes no stability claim along an unspecified
//! temperature path and supplies no thermal expansion, snap or rebonding law.

use super::{
    JkrCombinedError, JkrCombinedLoading, JkrCombinedLoadingError, JkrCombinedLoadingLimit,
    JkrCombinedLoadingState, JkrCombinedResponse, SourceJkrMixedMode, Work, dd, finite, poll,
};
use crate::normal_patch::adhesive::{
    JkrBranch,
    temperature::{AdhesiveThermalControl, SourceJkrThermalRefresh},
};
use fs_blake3::ContentHash;
use fs_exec::CancelGate;
use fs_math::dd::Dd;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrCombinedLoadingTemperatureReport {
    pub original_model_identity: ContentHash,
    pub refreshed_model_identity: ContentHash,
    pub before: Option<JkrCombinedResponse>,
    pub after: Option<JkrCombinedResponse>,
    pub previous_limit: Option<JkrCombinedLoadingLimit>,
    pub refreshed_limit: Option<JkrCombinedLoadingLimit>,
    pub normal_elastic_change_j: f64,
    pub tangential_elastic_change_j: f64,
    pub torsional_elastic_change_j: f64,
    pub surface_energy_change_j: f64,
    /// Signed external supply from changing the original source coefficients.
    pub parameter_work_j: f64,
    pub storage_change_j: f64,
    pub energy_residual_j: f64,
    pub energy_tolerance_j: f64,
    /// Original normal witness plus every successor preparation/endpoint probe.
    pub evaluations: usize,
}

/// Adopt the successor source, loading owner, state and work together. The
/// original caller-owned history is unchanged by failed or cancelled proposals.
#[derive(Clone, Copy, Debug)]
pub struct JkrCombinedLoadingTemperatureContinuation<'a> {
    previous: JkrCombinedLoadingState,
    loading: JkrCombinedLoading<'a>,
    state: JkrCombinedLoadingState,
    report: JkrCombinedLoadingTemperatureReport,
}
impl<'a> JkrCombinedLoadingTemperatureContinuation<'a> {
    pub const fn previous_state(&self) -> JkrCombinedLoadingState {
        self.previous
    }
    pub const fn loading(&self) -> &JkrCombinedLoading<'a> {
        &self.loading
    }
    pub const fn state(&self) -> JkrCombinedLoadingState {
        self.state
    }
    pub const fn report(&self) -> &JkrCombinedLoadingTemperatureReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        JkrCombinedLoading<'a>,
        JkrCombinedLoadingState,
        JkrCombinedLoadingTemperatureReport,
    ) {
        (self.loading, self.state, self.report)
    }
}

impl JkrCombinedLoading<'_> {
    /// Continue the original loading history at exactly held coordinates.
    ///
    /// `refresh` must be the original normal owner's held-indentation refresh,
    /// at this approach and on the same attached/detached branch as `previous`.
    /// Construct `successor` from its source and the unchanged mixed-mode
    /// envelope. The selected bulk/interface claims, independent temperatures,
    /// original geometry and all applicability limits remain authoritative.
    ///
    /// The old normal witness, new zero-load equilibrium, scalar first-boundary
    /// solve and actual loaded endpoint share one original evaluation budget.
    /// No-op continuation preserves the exact owner and accepted state after
    /// authenticating that witness. A changed-source detached successor skips
    /// all attached preparation, retains its accumulated release, and cannot
    /// initialize an attached state. Neither advances the loading step count.
    ///
    /// An unavailable attached endpoint refuses; it is never replaced by the
    /// pure-normal response, a clamped load, a fresh preload or an unowned snap.
    pub fn continue_at_temperature<'new>(
        &self,
        previous: JkrCombinedLoadingState,
        refresh: &SourceJkrThermalRefresh,
        successor: &'new SourceJkrMixedMode,
        gate: &CancelGate,
    ) -> Result<JkrCombinedLoadingTemperatureContinuation<'new>, JkrCombinedLoadingError> {
        self.admit(previous, gate)?;
        let normal = self.source.normal_source();
        if refresh.original_model_identity() != normal.identity()
            || refresh.before().normal_state_identity() != normal.normal_state().identity()
        {
            return Err(JkrCombinedLoadingError::StateMismatch);
        }
        let branch = if previous.contact.is_some() {
            JkrBranch::Attached
        } else {
            JkrBranch::Detached
        };
        if refresh.control() != AdhesiveThermalControl::HeldIndentation
            || refresh.mechanical_work_j() != 0.0
            || refresh.before().branch != branch
            || refresh.response().branch != branch
            || refresh.before().indentation_m.to_bits() != self.indentation_m.to_bits()
            || refresh.response().indentation_m.to_bits() != self.indentation_m.to_bits()
            || successor.normal_source() != refresh.source()
            || successor.envelope() != self.source.envelope()
            || refresh.source().effective_radius_m() != normal.effective_radius_m()
            || refresh.source().interaction_range_m() != normal.interaction_range_m()
            || refresh.source().envelope() != normal.envelope()
        {
            return Err(JkrCombinedLoadingError::Invalid(
                "original held-coordinate thermal successor and unchanged geometry required",
            ));
        }
        let mut work = Work::new(self.controls.maximum_evaluations, gate);
        work.charge()?;
        let canonical = match branch {
            JkrBranch::Attached => normal
                .evaluate_attached(self.indentation_m)
                .map_err(JkrCombinedError::from)?,
            JkrBranch::Detached => {
                normal
                    .detached_state(self.indentation_m)
                    .map_err(JkrCombinedError::from)?;
                normal.detached_response(self.indentation_m)
            }
        };
        if refresh.before() != &canonical {
            return Err(JkrCombinedLoadingError::StateMismatch);
        }
        poll(gate)?;
        let (loading, contact) = if successor == self.source {
            // Preserve cached boundary approximations, budgets and history bits.
            (
                JkrCombinedLoading {
                    source: successor,
                    indentation_m: self.indentation_m,
                    direction: self.direction,
                    controls: self.controls,
                    identity: self.identity,
                    zero: self.zero,
                    limit: self.limit,
                    preparation_evaluations: self.preparation_evaluations,
                },
                previous.contact,
            )
        } else {
            let mut loading = JkrCombinedLoading::unprepared(
                successor,
                self.indentation_m,
                self.direction,
                self.controls,
            )?;
            let contact = if previous.contact.is_some() {
                loading.prepare(&mut work)?;
                Some(loading.evaluate(previous.load_factor, &mut work)?)
            } else {
                None
            };
            (loading, contact)
        };
        let state = JkrCombinedLoadingState {
            identity: loading.identity,
            source_identity: successor.identity(),
            contact,
            ..previous
        };
        let changes = match (previous.contact, contact) {
            (Some(before), Some(after)) => [
                (dd(after.normal_elastic_energy_j) - dd(before.normal_elastic_energy_j)).to_f64(),
                (dd(after.tangential_elastic_energy_j) - dd(before.tangential_elastic_energy_j))
                    .to_f64(),
                (dd(after.torsional_elastic_energy_j) - dd(before.torsional_elastic_energy_j))
                    .to_f64(),
                (dd(after.surface_energy_j) - dd(before.surface_energy_j)).to_f64(),
            ],
            (None, None) => [0.0; 4],
            _ => {
                return Err(JkrCombinedLoadingError::Invalid(
                    "thermal continuation cannot change the accepted attachment branch",
                ));
            }
        };
        for value in changes {
            finite(value, "finite thermal storage change")?;
        }
        let parameter_work_j = finite(
            changes
                .into_iter()
                .fold(Dd::ZERO, |total, value| total + dd(value))
                .to_f64(),
            "finite thermal parameter work",
        )?;
        let storage_change_j = finite(
            (dd(state.stored_energy_j()) - dd(previous.stored_energy_j())).to_f64(),
            "finite thermal storage change",
        )?;
        let energy_residual_j = finite(
            (dd(storage_change_j) - dd(parameter_work_j)).to_f64(),
            "finite thermal storage balance",
        )?;
        let mut scale = parameter_work_j.abs();
        for point in [previous.contact, contact].into_iter().flatten() {
            scale += point.normal_elastic_energy_j
                + point.tangential_elastic_energy_j
                + point.torsional_elastic_energy_j
                - point.surface_energy_j;
        }
        let energy_tolerance_j = finite(
            self.controls.energy_absolute_tolerance_j
                + self.controls.energy_relative_tolerance * scale,
            "finite thermal balance tolerance",
        )?;
        if energy_residual_j.abs() > energy_tolerance_j {
            return Err(JkrCombinedLoadingError::EnergyBalance {
                residual_j: energy_residual_j,
                tolerance_j: energy_tolerance_j,
            });
        }
        poll(gate)?;
        Ok(JkrCombinedLoadingTemperatureContinuation {
            previous,
            loading,
            state,
            report: JkrCombinedLoadingTemperatureReport {
                original_model_identity: self.source.identity(),
                refreshed_model_identity: successor.identity(),
                before: previous.contact,
                after: contact,
                previous_limit: self.limit,
                refreshed_limit: loading.limit,
                normal_elastic_change_j: changes[0],
                tangential_elastic_change_j: changes[1],
                torsional_elastic_change_j: changes[2],
                surface_energy_change_j: changes[3],
                parameter_work_j,
                storage_change_j,
                energy_residual_j,
                energy_tolerance_j,
                evaluations: work.evaluations,
            },
        })
    }
}
