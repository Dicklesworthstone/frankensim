//! Original loaded combined adhesion follows finite-cell temperatures.
//!
//! The original native source owns every mechanical storage change. Thermal
//! transport and the solved massless-interface temperature use the same shared
//! finite-cell/area backend as the other admitted adhesive modes. Signed
//! coefficient work remains an explicitly named external supply.

use super::temperature_coupling::{
    self, NativeContinuation, NativeEndpoint, NativeTemperature, NativeView, ParameterWork,
};
use fs_contact::normal_patch::adhesive::mixed_mode::{
    SourceJkrMixedMode, combined::motion::temperature::JkrCombinedTemperatureReport,
};
use fs_contact::normal_patch::adhesive::temperature::{
    AdhesiveTemperatures, SourceJkrThermalRefresh,
};
use fs_exec::CancelGate;
use fs_mbd::tangential_friction::adhesive::jkr_combined::{
    PreparedRigidJkrCombinedPair, RigidJkrCombinedError, RigidJkrCombinedState,
};
use fs_thermal::enthalpy::DepositionSimulation;

pub use super::temperature_coupling::{
    JkrAdhesiveTemperatureHeatBinding as JkrCombinedTemperatureHeatBinding,
    JkrAdhesiveTemperatureHeatControls as JkrCombinedTemperatureHeatControls,
    JkrAdhesiveTemperatureHeatSource as JkrCombinedTemperatureHeatSource,
    JkrAdhesiveThermalControls as JkrCombinedThermalControls,
};
pub type JkrCombinedTemperatureHeatError<E> =
    temperature_coupling::JkrAdhesiveTemperatureHeatError<RigidJkrCombinedError, E>;
pub type JkrCombinedTemperatureHeatReport<P> =
    temperature_coupling::JkrAdhesiveTemperatureHeatReport<JkrCombinedTemperatureReport, P>;
pub type JkrCombinedTemperatureHeatCandidate<'t, P> =
    temperature_coupling::JkrAdhesiveTemperatureHeatCandidate<
        't,
        RigidJkrCombinedState,
        JkrCombinedTemperatureReport,
        P,
    >;

impl<E> From<RigidJkrCombinedError> for JkrCombinedTemperatureHeatError<E> {
    fn from(error: RigidJkrCombinedError) -> Self {
        Self::Native(error)
    }
}

impl NativeTemperature for PreparedRigidJkrCombinedPair<'_> {
    type State = RigidJkrCombinedState;
    type Report = JkrCombinedTemperatureReport;
    type Error = RigidJkrCombinedError;

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
pub fn propose_jkr_combined_temperature_heat<'t, P, E>(
    native: &PreparedRigidJkrCombinedPair<'_>,
    previous: RigidJkrCombinedState,
    thermal: &DepositionSimulation<'t>,
    binding: JkrCombinedTemperatureHeatBinding<'_>,
    controls: JkrCombinedTemperatureHeatControls,
    energy_guess_j: &[f64],
    temperature_guess_k: f64,
    resolve: impl Fn(
        &DepositionSimulation<'t>,
        AdhesiveTemperatures,
        &SourceJkrMixedMode,
    ) -> Result<JkrCombinedTemperatureHeatSource<P>, E>,
    gate: &CancelGate,
) -> Result<JkrCombinedTemperatureHeatCandidate<'t, P>, JkrCombinedTemperatureHeatError<E>> {
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
}
