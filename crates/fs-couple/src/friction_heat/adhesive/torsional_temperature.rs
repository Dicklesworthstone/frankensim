//! Original loaded torsional adhesion follows finite-cell temperatures.
//!
//! The original native source owns every mechanical storage change. Thermal
//! transport and the solved massless-interface temperature use the same shared
//! finite-cell/area backend as the other admitted adhesive modes. Signed
//! coefficient work remains an explicitly named external supply.

use super::temperature_coupling::{
    self, NativeContinuation, NativeEndpoint, NativeTemperature, NativeView, ParameterWork,
};
use fs_contact::normal_patch::adhesive::mixed_mode::{
    JkrMixedModeError, SourceJkrMixedMode,
    torsion::motion::temperature::JkrTorsionalTemperatureReport,
};
use fs_contact::normal_patch::adhesive::temperature::{
    AdhesiveTemperatureError, AdhesiveTemperatures, SourceJkrThermalRefresh,
};
use fs_exec::CancelGate;
use fs_mbd::tangential_friction::adhesive::jkr_torsional::{
    PreparedRigidJkrTorsionalPair, RigidJkrTorsionalError, RigidJkrTorsionalState,
};
use fs_solver::{NewtonError, NewtonReport, SolverRunProgress};
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};

pub use super::temperature_coupling::{
    JkrAdhesiveTemperatureHeatBinding as JkrTorsionalTemperatureHeatBinding,
    JkrAdhesiveTemperatureHeatControls as JkrTorsionalTemperatureHeatControls,
    JkrAdhesiveTemperatureHeatSource as JkrTorsionalTemperatureHeatSource,
    JkrAdhesiveThermalControls as JkrTorsionalThermalControls,
};
#[derive(Debug)]
pub enum JkrTorsionalTemperatureHeatError<E> {
    Invalid(&'static str),
    Native(RigidJkrTorsionalError),
    Temperature(AdhesiveTemperatureError),
    MixedMode(JkrMixedModeError),
    Thermal(DepositionError),
    Source(E),
    Newton(NewtonError),
    NotConverged(NewtonReport),
    Interrupted(SolverRunProgress),
    TrialBudget,
    EnergyBalance {
        residual_j: f64,
        absolute_stage_residual_j: f64,
    },
}
impl<E> From<RigidJkrTorsionalError> for JkrTorsionalTemperatureHeatError<E> {
    fn from(e: RigidJkrTorsionalError) -> Self {
        Self::Native(e)
    }
}
impl<E> From<AdhesiveTemperatureError> for JkrTorsionalTemperatureHeatError<E> {
    fn from(e: AdhesiveTemperatureError) -> Self {
        Self::Temperature(e)
    }
}
impl<E> From<JkrMixedModeError> for JkrTorsionalTemperatureHeatError<E> {
    fn from(e: JkrMixedModeError) -> Self {
        Self::MixedMode(e)
    }
}
impl<E> From<DepositionError> for JkrTorsionalTemperatureHeatError<E> {
    fn from(e: DepositionError) -> Self {
        Self::Thermal(e)
    }
}
impl<E> From<NewtonError> for JkrTorsionalTemperatureHeatError<E> {
    fn from(e: NewtonError) -> Self {
        Self::Newton(e)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for JkrTorsionalTemperatureHeatError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "JKR torsional temperature heat refused: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for JkrTorsionalTemperatureHeatError<E> {}

pub type JkrTorsionalTemperatureHeatReport<P> =
    temperature_coupling::JkrAdhesiveTemperatureHeatReport<JkrTorsionalTemperatureReport, P>;
pub type JkrTorsionalTemperatureHeatCandidate<'t, P> =
    temperature_coupling::JkrAdhesiveTemperatureHeatCandidate<
        't,
        RigidJkrTorsionalState,
        JkrTorsionalTemperatureReport,
        P,
    >;

// Preserve the original public enum (including variant imports) while sharing
// the numerical driver. This mapping changes neither the error nor its payload.
impl<E> From<temperature_coupling::JkrAdhesiveTemperatureHeatError<RigidJkrTorsionalError, E>>
    for JkrTorsionalTemperatureHeatError<E>
{
    fn from(
        error: temperature_coupling::JkrAdhesiveTemperatureHeatError<RigidJkrTorsionalError, E>,
    ) -> Self {
        use temperature_coupling::JkrAdhesiveTemperatureHeatError as Shared;
        match error {
            Shared::Invalid(message) => Self::Invalid(message),
            Shared::Native(error) => Self::Native(error),
            Shared::Temperature(error) => Self::Temperature(error),
            Shared::MixedMode(error) => Self::MixedMode(error),
            Shared::Thermal(error) => Self::Thermal(error),
            Shared::Source(error) => Self::Source(error),
            Shared::Newton(error) => Self::Newton(error),
            Shared::NotConverged(report) => Self::NotConverged(report),
            Shared::Interrupted(progress) => Self::Interrupted(progress),
            Shared::TrialBudget => Self::TrialBudget,
            Shared::EnergyBalance {
                residual_j,
                absolute_stage_residual_j,
            } => Self::EnergyBalance {
                residual_j,
                absolute_stage_residual_j,
            },
        }
    }
}

impl NativeTemperature for PreparedRigidJkrTorsionalPair<'_> {
    type State = RigidJkrTorsionalState;
    type Report = JkrTorsionalTemperatureReport;
    type Error = RigidJkrTorsionalError;

    fn view(&self) -> NativeView<'_> {
        NativeView {
            source: self.contact().source,
            normal: self.contact().frame.normal(),
            reaction_arms: self.reaction_arms_body(),
        }
    }
    fn endpoint(previous: Self::State) -> NativeEndpoint {
        NativeEndpoint {
            time_s: previous.time_s(),
            approach_m: previous.coordinates()[0],
            contact_area_m2: previous.motion_state().contact().contact_area_m2,
        }
    }
    fn refresh(
        &self,
        previous: Self::State,
        refresh: &SourceJkrThermalRefresh,
        successor: &SourceJkrMixedMode,
        gate: &CancelGate,
    ) -> Result<NativeContinuation<Self::State, Self::Report>, Self::Error> {
        let continuation = self.continue_at_temperature(previous, refresh, successor, gate)?;
        let report = *continuation.report();
        Ok(NativeContinuation {
            bodies: *continuation.bodies(),
            state: continuation.state(),
            area_m2: report.after.contact_area_m2,
            work: ParameterWork {
                storage_change_j: report.storage_change_j,
                parameter_work_j: report.parameter_work_j,
                energy_residual_j: report.energy_residual_j,
            },
            report,
        })
    }
}

/// Advance finite calorics to the already staged native endpoint, solving the
/// original two bulk temperatures, massless interface temperature and common
/// relaxed contact area together. Native geometry must already admit the state.
/// Incoming bulk source temperatures must agree with their original receivers.
///
/// The source callback stages original conductivity and other thermal supplies
/// at each actual refreshed endpoint. It must not publish side effects, reinsert
/// parameter work as heat, or replace this route's conduction link. The returned
/// candidate owns source, native history, thermal state and route together.
/// No impulse, geometry correction, entropy law or thermal snap is inferred.
#[allow(clippy::too_many_arguments)]
pub fn propose_jkr_torsional_temperature_heat<'t, P, E>(
    native: &PreparedRigidJkrTorsionalPair<'_>,
    previous: RigidJkrTorsionalState,
    thermal: &DepositionSimulation<'t>,
    binding: JkrTorsionalTemperatureHeatBinding<'_>,
    controls: JkrTorsionalTemperatureHeatControls,
    energy_guess_j: &[f64],
    temperature_guess_k: f64,
    resolve: impl Fn(
        &DepositionSimulation<'t>,
        AdhesiveTemperatures,
        &SourceJkrMixedMode,
    ) -> Result<JkrTorsionalTemperatureHeatSource<P>, E>,
    gate: &CancelGate,
) -> Result<JkrTorsionalTemperatureHeatCandidate<'t, P>, JkrTorsionalTemperatureHeatError<E>> {
    temperature_coupling::propose(
        native,
        previous,
        thermal,
        binding,
        controls,
        energy_guess_j,
        temperature_guess_k,
        resolve,
        gate,
    )
    .map_err(Into::into)
}
