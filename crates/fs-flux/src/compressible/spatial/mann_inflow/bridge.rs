//! One replayable OU path across native inviscid and molecular stage times.

use super::*;
use crate::compressible::spatial::{StepControl, StepReport};
use fs_qty::Time;
use std::cell::RefCell;

mod molecular;

/// Native accepted physical balances and actual cumulative atmosphere work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MannReservoirBridgeReport3<R> {
    pub flow: R,
    /// All faces observed across stages and discarded complete trials.
    pub evaluated_faces: usize,
    pub mode_samples: usize,
    /// Conditional normal variates, including repeated paths after retries.
    /// Preparing the bridge's original right endpoint is budgeted separately.
    pub normal_draws: usize,
    /// One amplitude-frame update per newly requested stage time. All faces
    /// at the same consecutive time share that update without random draws.
    pub frame_updates: usize,
    /// Requested duration clipped to this bridge's right endpoint and the
    /// advective spatial-phase bound; the original fluid may shorten it more.
    pub maximum_duration: Time,
}

#[derive(Clone, Copy)]
struct NativeGeometry {
    grid: CartesianGrid,
    boundaries: Boundaries3d,
    selected: [bool; 6],
    time_s: f64,
}
impl NativeGeometry {
    fn frozen(flow: &Euler3d) -> Self {
        Self {
            grid: *flow.grid(),
            boundaries: flow.boundaries(),
            selected: flow.prepared_reservoir_profile_sides(),
            time_s: flow.time().value(),
        }
    }
    fn mixture(flow: &MixtureEuler3d) -> Self {
        Self {
            grid: *flow.grid(),
            boundaries: flow.boundaries(),
            selected: flow.prepared_reservoir_profile_sides(),
            time_s: flow.time().value(),
        }
    }
}

struct BridgeSampler<'a> {
    inlet: &'a MannReservoirInflow3,
    bridge: &'a MannOuBridge,
    frame: &'a mut MannOuBridgeFrame,
    maximum_duration: Time,
    maximum_face_evaluations: usize,
    maximum_normal_draws: usize,
    sampled_time_s: Option<f64>,
    evaluated_faces: usize,
    mode_count: usize,
    normal_draws: usize,
    frame_updates: usize,
}

impl<'a> BridgeSampler<'a> {
    fn new(
        inlet: &'a MannReservoirInflow3,
        bridge: &'a MannOuBridge,
        frame: &'a mut MannOuBridgeFrame,
        native: NativeGeometry,
        maximum_duration: Time,
        maximum_normal_draws: usize,
    ) -> Result<Self, MannReservoirError3> {
        inlet.validate_geometry(&native.grid, native.boundaries)?;
        let requested = maximum_duration.value();
        let time = native.time_s;
        if native.selected != inlet.selected
            || !requested.is_finite()
            || requested <= 0.0
            || !time.is_finite()
            || time < bridge.start_time_s()
            || time >= bridge.end_time_s()
        {
            return Err(MannReservoirError3::Input(
                "prepared Mann sides, positive duration and fluid time within this bridge required",
            ));
        }
        let mode_count = bridge.model().mode_count();
        if mode_count != inlet.mode_samples / inlet.faces {
            return Err(MannReservoirError3::Input(
                "Mann bridge mode count differs from prepared inlet",
            ));
        }
        let phase_rate = inlet.k_max_per_m * speed(bridge.advection_mps());
        if !phase_rate.is_finite() {
            return Err(MannReservoirError3::AdvectiveNyquist {
                phase_bound_rad: phase_rate,
            });
        }
        let mut duration = requested.min(bridge.end_time_s() - time);
        if phase_rate > 0.0 {
            duration =
                duration.min(core::f64::consts::PI * (1.0 - 8.0 * f64::EPSILON) / phase_rate);
        }
        // Roundoff in endpoint subtraction must never observe beyond the
        // retained path. Two downward ulps suffice for one subtraction/addition.
        for _ in 0..2 {
            if time + duration > bridge.end_time_s() {
                duration = duration.next_down();
            }
        }
        if !duration.is_finite() || time + duration <= time || time + duration > bridge.end_time_s()
        {
            return Err(MannReservoirError3::Input(
                "Mann bridge duration is not representable",
            ));
        }
        let required = inlet
            .mode_samples
            .checked_mul(2)
            .ok_or(MannReservoirError3::Input("Mann bridge face work overflow"))?;
        if required > inlet.config.maximum_mode_samples {
            return Err(MannReservoirError3::Work {
                required,
                maximum: inlet.config.maximum_mode_samples,
            });
        }
        Ok(Self {
            inlet,
            bridge,
            frame,
            maximum_duration: Time::new(duration),
            maximum_face_evaluations: inlet.config.maximum_mode_samples / mode_count,
            maximum_normal_draws,
            sampled_time_s: None,
            evaluated_faces: 0,
            mode_count,
            normal_draws: 0,
            frame_updates: 0,
        })
    }

    fn fluctuation(
        &mut self,
        time: f64,
        side: BoundarySide,
        face: usize,
        cancelled: impl FnMut() -> bool,
    ) -> Result<[f64; 3], MannReservoirError3> {
        self.evaluated_faces += 1; // Native owner admits/charges every callback.
        if self.sampled_time_s != Some(time) {
            let draws = self.bridge.required_normal_draws(time)?;
            let required =
                self.normal_draws
                    .checked_add(draws)
                    .ok_or(MannReservoirError3::Input(
                        "Mann bridge draw counter overflow",
                    ))?;
            if required > self.maximum_normal_draws {
                return Err(MannReservoirError3::BridgeWork {
                    required,
                    maximum: self.maximum_normal_draws,
                });
            }
            let report = self.bridge.sample_frame(
                time,
                self.frame,
                self.maximum_normal_draws - self.normal_draws,
                cancelled,
            )?;
            debug_assert_eq!(report.normal_draws, draws);
            self.normal_draws = required;
            self.frame_updates += 1;
            self.sampled_time_s = Some(time);
        }
        self.inlet
            .fluctuation_at(side, face, [0.0; 3], |point, filter| {
                Ok(self.frame.sample_filtered(point, filter)?.u)
            })
    }

    fn report<R>(&self, flow: R) -> MannReservoirBridgeReport3<R> {
        MannReservoirBridgeReport3 {
            flow,
            evaluated_faces: self.evaluated_faces,
            mode_samples: self.evaluated_faces * self.mode_count,
            normal_draws: self.normal_draws,
            frame_updates: self.frame_updates,
            maximum_duration: self.maximum_duration,
        }
    }
}

type FrozenProvider<'a> =
    dyn FnMut(f64, BoundarySide, usize, [f64; 3]) -> Result<ReservoirFaceState, Error> + 'a;
type MixtureProvider<'a> = dyn FnMut(
        f64,
        BoundarySide,
        usize,
        [f64; 3],
        &mut [f64],
    ) -> Result<MixtureReservoirFaceState, MixtureError>
    + 'a;

impl MannReservoirInflow3 {
    /// Advance original frozen-composition Euler flow along one conditional
    /// OU realization, observing the original full Mann tensor at every native
    /// SSP stage and retry. The bridge retains both original fixed-clock
    /// endpoints; neither this call nor a failed attempt advances that owner.
    ///
    /// Prepare exactly these sides with `apply_frozen(BridgeFrame, Prepare)`
    /// before stepping. `baseline` must be pure and replayable and supplies
    /// density, pressure and mean velocity. One face/mode budget and one draw
    /// budget span the complete native call. A reusable frame caches only the
    /// most recent observation; it is scratch and may change on flow failure.
    /// Retrying still recovers the same path. No stochastic convergence order
    /// or temporal bandlimit for OU noise is asserted.
    #[allow(clippy::too_many_arguments)]
    pub fn step_frozen_bridge(
        &self,
        flow: &mut Euler3d,
        bridge: &MannOuBridge,
        frame: &mut MannOuBridgeFrame,
        maximum_duration: Time,
        controls: StepControl,
        maximum_normal_draws: usize,
        baseline: impl FnMut(f64, BoundarySide, usize, [f64; 3]) -> Result<ReservoirFaceState, Error>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MannReservoirBridgeReport3<StepReport>, MannReservoirError3> {
        self.drive_frozen_bridge(
            bridge,
            frame,
            NativeGeometry::frozen(flow),
            maximum_duration,
            maximum_normal_draws,
            baseline,
            cancelled,
            |duration, budget, sample, cancel| {
                flow.step_with_reservoir_profiles(duration, controls, budget, sample, cancel)
                    .map_err(Into::into)
            },
        )
    }

    /// Same conditional path for the original evolving-composition Euler
    /// owner. The baseline fills every canonical species density. This adapter
    /// adds velocity fluctuations only; calorics and fluxes stay source-owned.
    #[allow(clippy::too_many_arguments)]
    pub fn step_mixture_bridge(
        &self,
        flow: &mut MixtureEuler3d,
        bridge: &MannOuBridge,
        frame: &mut MannOuBridgeFrame,
        maximum_duration: Time,
        controls: StepControl,
        maximum_normal_draws: usize,
        baseline: impl FnMut(
            f64,
            BoundarySide,
            usize,
            [f64; 3],
            &mut [f64],
        ) -> Result<MixtureReservoirFaceState, MixtureError>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MannReservoirBridgeReport3<super::super::mixture::StepReport>, MannReservoirError3>
    {
        self.drive_mixture_bridge(
            bridge,
            frame,
            NativeGeometry::mixture(flow),
            maximum_duration,
            maximum_normal_draws,
            baseline,
            cancelled,
            |duration, budget, sample, cancel| {
                flow.step_with_reservoir_profiles(duration, controls, budget, sample, cancel)
                    .map_err(Into::into)
            },
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn drive_frozen_bridge<R>(
        &self,
        bridge: &MannOuBridge,
        frame: &mut MannOuBridgeFrame,
        native: NativeGeometry,
        maximum_duration: Time,
        maximum_normal_draws: usize,
        mut baseline: impl FnMut(
            f64,
            BoundarySide,
            usize,
            [f64; 3],
        ) -> Result<ReservoirFaceState, Error>,
        mut cancelled: impl FnMut() -> bool,
        step: impl FnOnce(
            Time,
            usize,
            &mut FrozenProvider<'_>,
            &mut dyn FnMut() -> bool,
        ) -> Result<R, MannReservoirError3>,
    ) -> Result<MannReservoirBridgeReport3<R>, MannReservoirError3> {
        poll(&mut cancelled)?;
        let mut sampler = BridgeSampler::new(
            self,
            bridge,
            frame,
            native,
            maximum_duration,
            maximum_normal_draws,
        )?;
        let maximum_duration = sampler.maximum_duration;
        let budget = sampler.maximum_face_evaluations;
        let cancelled = RefCell::new(cancelled);
        let mut sample_error = None;
        let result = step(
            maximum_duration,
            budget,
            &mut |time, side, face, point| {
                let u = sampler
                    .fluctuation(time, side, face, || (cancelled.borrow_mut())())
                    .map_err(|error| {
                        sample_error = Some(error);
                        Error::Invalid("Mann bridge stage observation refused")
                    })?;
                let mut mean = baseline(time, side, face, point)?;
                for (value, fluctuation) in mean.velocity.iter_mut().zip(u) {
                    *value += fluctuation;
                }
                Ok(mean)
            },
            &mut || (cancelled.borrow_mut())(),
        );
        match result {
            Ok(report) => Ok(sampler.report(report)),
            Err(error) => Err(sample_error.unwrap_or(error)),
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn drive_mixture_bridge<R>(
        &self,
        bridge: &MannOuBridge,
        frame: &mut MannOuBridgeFrame,
        native: NativeGeometry,
        maximum_duration: Time,
        maximum_normal_draws: usize,
        mut baseline: impl FnMut(
            f64,
            BoundarySide,
            usize,
            [f64; 3],
            &mut [f64],
        ) -> Result<MixtureReservoirFaceState, MixtureError>,
        mut cancelled: impl FnMut() -> bool,
        step: impl FnOnce(
            Time,
            usize,
            &mut MixtureProvider<'_>,
            &mut dyn FnMut() -> bool,
        ) -> Result<R, MannReservoirError3>,
    ) -> Result<MannReservoirBridgeReport3<R>, MannReservoirError3> {
        poll(&mut cancelled)?;
        let mut sampler = BridgeSampler::new(
            self,
            bridge,
            frame,
            native,
            maximum_duration,
            maximum_normal_draws,
        )?;
        let maximum_duration = sampler.maximum_duration;
        let budget = sampler.maximum_face_evaluations;
        let cancelled = RefCell::new(cancelled);
        let mut sample_error = None;
        let result = step(
            maximum_duration,
            budget,
            &mut |time, side, face, point, densities| {
                let u = sampler
                    .fluctuation(time, side, face, || (cancelled.borrow_mut())())
                    .map_err(|error| {
                        sample_error = Some(error);
                        MixtureError::Invalid("Mann bridge stage observation refused")
                    })?;
                let mut mean = baseline(time, side, face, point, densities)?;
                for (value, fluctuation) in mean.velocity.iter_mut().zip(u) {
                    *value += fluctuation;
                }
                Ok(mean)
            },
            &mut || (cancelled.borrow_mut())(),
        );
        match result {
            Ok(report) => Ok(sampler.report(report)),
            Err(error) => Err(sample_error.unwrap_or(error)),
        }
    }
}
