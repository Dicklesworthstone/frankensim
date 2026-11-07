//! Whole-duration heterogeneous diffusion with one original finite solution.

use super::radial::{observe_radial, radial_water_distance};
use super::*;
use crate::moisture::{
    radial::layered::MoistureLayeredRadial,
    solution::{
        FiniteMoistureSolution, SolutionMoistureError, SolutionMoistureSolve, SolutionMoistureStep,
    },
};

mod sdirk;

impl<R: LayeredTransport> LayeredTransport for SolutionMoistureStep<R> {
    fn transport(&self) -> &MoistureLayeredReport {
        self.transport.transport()
    }
}

#[allow(clippy::too_many_arguments, clippy::type_complexity)]
impl FiniteMoistureSolution {
    /// Complete finite-solution exchange with an original heterogeneous slab.
    /// Full and half backward-Euler probes each solve the native solvent-volume
    /// root and all original nonlinear material/interface/surface laws. Only
    /// accepted fine transfers enter the per-layer and exterior water accounts.
    ///
    /// Accuracy resolves each local moisture and source activity, solution
    /// activity and actual solvent mass, plus the L1 sum of local bound-water
    /// discrepancies weighted by each original layer's dry cell masses. The
    /// latter uses `reservoir_water_absolute_tolerance_kg`, so opposing internal
    /// errors cannot hide behind zero net water or constant pure-water activity.
    /// Relative scales refer to changes from the last accepted state.
    ///
    /// The shared bounded controller checks original-to-final bound/solvent
    /// water once against the fixed whole-duration tolerance. Source and phase
    /// admission stay native. Kinetic dissolution retains its first-order
    /// endpoint-volume law; `time.with_kinetic_phase_accuracy(...)` additionally
    /// controls every actual dissolved/crystal stock with an explicit molar scale.
    /// Geometry, temperature and exterior laws stay fixed, and the original
    /// equal-pressure vapor bridge does not imply liquid immersion. Any failure
    /// discards both endpoints, including successful private prefixes. These
    /// local asymptotic comparisons do not assert a global error guarantee.
    pub fn advance_layered_adaptive_to(
        &self,
        moisture: &MoistureLayeredSlab,
        target_s: f64,
        exposed: [bool; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        physical: SolutionMoistureSolve,
        time: impl Into<SolutionMoistureTimeControls>,
        cx: &Cx<'_>,
    ) -> Result<
        (
            Self,
            MoistureLayeredSlab,
            LayeredMoistureAdvanceReport<SolutionMoistureStep<MoistureLayeredReport>>,
        ),
        MoistureAdvanceError<SolutionMoistureError>,
    > {
        checkpoint(cx)?;
        let accuracy = time.into();
        accuracy.admit(self, cx)?;
        let time = accuracy.transport;
        let initial = (self.clone(), moisture.clone());
        let ((reservoir, field), report) = evolve_fields_with_method(
            &initial,
            target_s,
            time,
            TimeMethod::FirstOrder,
            |state, dt| {
                let (bath, field, step) = state
                    .0
                    .advance_layered(&state.1, dt, exposed, transport, nonlinear, physical, cx)?
                    .into_parts();
                Ok(((bath, field), LayeredStep(step)))
            },
            |state| {
                let mut snapshot = observe(&state.1);
                snapshot.reservoir_water_kg = Some(state.0.solvent_mass_kg());
                snapshot.activity = Some(state.0.solution().water_activity()?);
                Ok(snapshot)
            },
            |old, coarse, fine| {
                let water = water_distance(
                    old.1
                        .layers()
                        .iter()
                        .zip(coarse.1.layers())
                        .zip(fine.1.layers())
                        .map(|((a, b), c)| (a.slab().field(), b.slab().field(), c.slab().field())),
                    time,
                    cx,
                )?;
                Ok(water.max(accuracy.phase_error(&old.0, &coarse.0, &fine.0, cx)?))
            },
            cx,
        )?;
        Ok((
            reservoir,
            field,
            finish(report, time.water_balance_tolerance_kg, cx)?,
        ))
    }

    /// Complete finite-solution exchange with a coated cylinder or sphere.
    /// Each probe retains exact curved paths, the actual unequal shell masses
    /// and the original sorption/diffusivity of every layer. Only real exposed
    /// inner/outer surfaces contact the solution; material joins keep their own
    /// partition and resistance laws. A solid center cannot be exposed.
    ///
    /// The water/activity comparisons, bounded retries, kinetic-phase accuracy
    /// limits and atomic whole-duration publication are those documented by
    /// [`Self::advance_layered_adaptive_to`]. Every shell is weighted by its own
    /// dry mass in both the temporal discrepancy and final inventory. Each
    /// internal interface debit has an equal adjacent-layer credit; the whole
    /// body's summary exposes only the two exterior water transfers.
    pub fn advance_layered_radial_adaptive_to(
        &self,
        moisture: &MoistureLayeredRadial,
        target_s: f64,
        exposed: [bool; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        physical: SolutionMoistureSolve,
        time: impl Into<SolutionMoistureTimeControls>,
        cx: &Cx<'_>,
    ) -> Result<
        (
            Self,
            MoistureLayeredRadial,
            LayeredMoistureAdvanceReport<SolutionMoistureStep<MoistureLayeredReport>>,
        ),
        MoistureAdvanceError<SolutionMoistureError>,
    > {
        checkpoint(cx)?;
        let accuracy = time.into();
        accuracy.admit(self, cx)?;
        let time = accuracy.transport;
        let initial = (self.clone(), moisture.clone());
        let ((reservoir, field), report) = evolve_fields_with_method(
            &initial,
            target_s,
            time,
            TimeMethod::FirstOrder,
            |state, dt| {
                let (bath, field, step) = state
                    .0
                    .advance_layered_radial(
                        &state.1, dt, exposed, transport, nonlinear, physical, cx,
                    )?
                    .into_parts();
                Ok(((bath, field), LayeredStep(step)))
            },
            |state| {
                let mut snapshot = observe_radial(&state.1);
                snapshot.reservoir_water_kg = Some(state.0.solvent_mass_kg());
                snapshot.activity = Some(state.0.solution().water_activity()?);
                Ok(snapshot)
            },
            |old, coarse, fine| {
                Ok(radial_water_distance(&old.1, &coarse.1, &fine.1, time, cx)?
                    .max(accuracy.phase_error(&old.0, &coarse.0, &fine.0, cx)?))
            },
            cx,
        )?;
        Ok((
            reservoir,
            field,
            finish(report, time.water_balance_tolerance_kg, cx)?,
        ))
    }
}
