//! Adaptive finite caloric feedback with original temperature/depth hardness.
//!
//! The ordinary prescribed-velocity controller retains the two native half
//! steps. In addition to its mechanical, thermal-energy and work defects, this
//! owner checks each body's temperature, each original exposed H(T,z), and the
//! native effective pressure cap. A stronger body's evolving hardness therefore
//! cannot disappear behind the softer body's pressure limit. These are local
//! split-error estimates, not global error or constitutive-validity bounds.

use super::super::adaptive::{self as shared, MeasuredDefect, MechanicalState, NativeStep};
use super::*;
use shared::{PlasticSlidingMarchControls, PlasticSlidingMarchDefect, PlasticSlidingMarchError};

/// Independent physical scales supplementing the original native march gates.
#[derive(Clone, Copy, Debug)]
pub struct DepthSlidingMarchControls {
    pub march: PlasticSlidingMarchControls,
    /// Absolute [K] scale applied independently to each original receiving body.
    pub body_temperature_k: f64,
    /// Absolute [Pa] scale applied to each body's exposed original hardness and
    /// each cell's effective cap. The shared relative tolerance uses changes
    /// from the accepted interval origin, never the background hardness.
    pub hardness_pa: f64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct DepthSlidingMarchDefect {
    pub common: PlasticSlidingMarchDefect,
    pub body_temperature_k: f64,
    /// Maximum original body-hardness or effective-cap difference [Pa].
    pub hardness_pa: f64,
    pub normalized: f64,
}

pub type DepthSlidingMarchError<E = Infallible> =
    PlasticSlidingMarchError<E, DepthSlidingMarchDefect>;
pub type DepthSlidingMarchReport<P> =
    shared::PlasticSlidingMarchReport<P, DepthSlidingHeatMechanics>;

/// Complete accepted successors and original fine-path source/work receipts.
pub struct DepthSlidingMarchCandidate<'a, 't, P> {
    contact: DepthSlidingHeatContact<'a>,
    thermal: DepositionSimulation<'t>,
    report: DepthSlidingMarchReport<P>,
}

impl<'a, 't, P> DepthSlidingMarchCandidate<'a, 't, P> {
    pub const fn contact(&self) -> &DepthSlidingHeatContact<'a> {
        &self.contact
    }
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub const fn report(&self) -> &DepthSlidingMarchReport<P> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepthSlidingHeatContact<'a>,
        DepositionSimulation<'t>,
        DepthSlidingMarchReport<P>,
    ) {
        (self.contact, self.thermal, self.report)
    }
}

impl<'a> DepthSlidingHeatContact<'a> {
    /// Propose a whole constant-actual-velocity span at one original normal
    /// boundary. Each native interval has zero tangential body compliance.
    ///
    /// Only original fine candidates contribute accepted work, wear, plastic
    /// depth, slip, caloric state and source payloads. Coarse/rejected solves
    /// still consume budgets. A refusal after private accepted steps returns no
    /// partial span. Original source, contact-resistance, footprint, removal,
    /// work and nonlinear gates remain strict; they are never retried as an
    /// apparent local-error rejection. Thermal events and worn-mass export
    /// remain explicit endpoint operations.
    #[allow(clippy::too_many_arguments)]
    pub fn march_prescribed_velocity_candidate<'t, P, E>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        velocity_m_s: &[[f64; 2]],
        boundary: PlasticSlidingHeatBoundary,
        step_controls: PlasticSlidingHeatControls,
        controls: DepthSlidingMarchControls,
        interface_temperature_guess_k: f64,
        resolve: impl Fn(&DepositionSimulation<'t>, f64) -> Result<PlasticSlidingHeatSource<P>, E>,
        cx: &Cx<'_>,
    ) -> Result<DepthSlidingMarchCandidate<'a, 't, P>, DepthSlidingMarchError<E>> {
        checkpoint(cx)?;
        self.admit_thermal(thermal)?;
        if [controls.body_temperature_k, controls.hardness_pa]
            .iter()
            .any(|value| !value.is_finite() || *value < 0.0)
        {
            return Err(PlasticSlidingMarchError::Invalid(
                "finite nonnegative body-temperature and hardness scales required",
            ));
        }
        let result = shared::march(
            self,
            thermal,
            target_s,
            velocity_m_s,
            self.source.tangential_states().len(),
            step_controls.maximum_cells,
            controls.march,
            interface_temperature_guess_k,
            |contact, state, end, motion, energy_guess, guess| {
                contact
                    .step_candidate(
                        state,
                        end,
                        motion,
                        boundary,
                        step_controls,
                        energy_guess,
                        guess,
                        &resolve,
                        cx,
                    )
                    .map(|candidate| {
                        let (contact, thermal, report) = candidate.into_parts();
                        NativeStep {
                            contact,
                            thermal,
                            report,
                        }
                    })
            },
            |old, old_thermal, coarse, first, fine, _, cx| {
                compare(old, old_thermal, coarse, first, fine, controls, cx)
            },
            accounts,
            cx,
        )?;
        Ok(DepthSlidingMarchCandidate {
            contact: result.contact,
            thermal: result.thermal,
            report: result.report,
        })
    }
}

fn accounts<P>(report: &DepthSlidingHeatReport<P>) -> [f64; 13] {
    // The depth owner has THREE native normal accounts: thermal refresh,
    // recession, and exposure of a new hardness cap. None may be omitted.
    let elastic = report.mechanics.wear_normals().iter().fold(
        report
            .mechanics
            .temperature_normal()
            .elastic_energy_change_j,
        |sum, continuation| sum + continuation.normal.elastic_energy_change_j,
    );
    shared::account_values(report, elastic)
}

fn view<'s>(contact: &'s DepthSlidingHeatContact<'_>) -> MechanicalState<'s> {
    MechanicalState {
        normal: contact.source.wear_state().contact_state(),
        recession_m: contact.source.wear_state().recession_m(),
        tangential: contact.source.tangential_states(),
    }
}

#[allow(clippy::too_many_arguments)]
fn compare<P, E>(
    old: &DepthSlidingHeatContact<'_>,
    old_thermal: &DepositionSimulation<'_>,
    coarse: &NativeStep<'_, DepthSlidingHeatContact<'_>, P, DepthSlidingHeatMechanics>,
    first: &NativeStep<'_, DepthSlidingHeatContact<'_>, P, DepthSlidingHeatMechanics>,
    fine: &NativeStep<'_, DepthSlidingHeatContact<'_>, P, DepthSlidingHeatMechanics>,
    controls: DepthSlidingMarchControls,
    cx: &Cx<'_>,
) -> Result<MeasuredDefect<DepthSlidingMarchDefect>, DepthSlidingMarchError<E>> {
    let common = shared::compare_common(
        view(old),
        view(&coarse.contact),
        view(&fine.contact),
        old_thermal,
        coarse,
        first,
        fine,
        controls.march.tolerance,
        accounts,
        cx,
    )?;
    let mut out = DepthSlidingMarchDefect {
        common,
        normalized: common.normalized,
        ..DepthSlidingMarchDefect::default()
    };
    let relative = controls.march.tolerance.relative;
    for id in old.route.cells() {
        let (error, normalized) = shared::difference(
            coarse.thermal.temperature_k(id)?,
            fine.thermal.temperature_k(id)?,
            old_thermal.temperature_k(id)?,
            controls.body_temperature_k,
            relative,
        )?;
        out.body_temperature_k = out.body_temperature_k.max(error);
        out.normalized = out.normalized.max(normalized);
    }
    let o = old.source.wear_state();
    let a = coarse.contact.source.wear_state();
    let b = fine.contact.source.wear_state();
    for i in 0..o.contact_state().pressure_pa().len() {
        if i % 64 == 0 {
            checkpoint(cx)?;
        }
        for side in 0..2 {
            let at = |state: &fs_contact::normal_patch::plastic_source::depth_wear::SourceDepthHardnessWear| {
                state.hardness_profiles()[side]
                    .hardness_pa_at(state.recession_m()[side][i])
                    .map_err(|error| PlasticSlidingHeatError::DepthContact(
                        HardnessSlidingError::Depth(error.into()),
                    ))
            };
            let (error, normalized) =
                shared::difference(at(a)?, at(b)?, at(o)?, controls.hardness_pa, relative)?;
            out.hardness_pa = out.hardness_pa.max(error);
            out.normalized = out.normalized.max(normalized);
        }
        let (error, normalized) = shared::difference(
            a.model().config().hardness_pa[i],
            b.model().config().hardness_pa[i],
            o.model().config().hardness_pa[i],
            controls.hardness_pa,
            relative,
        )?;
        out.hardness_pa = out.hardness_pa.max(error);
        out.normalized = out.normalized.max(normalized);
    }
    checkpoint(cx)?;
    Ok(MeasuredDefect {
        normalized: out.normalized,
        value: out,
    })
}
