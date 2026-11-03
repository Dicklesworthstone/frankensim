//! Molecular scattering survival along an already connected optical branch.
//!
//! The original deterministic solver owns geometry, Fresnel power, true Beer
//! absorption and the receiving-area Jacobian. One additional native bundle
//! walk replays its solved launch while independently measuring ONLY survival
//! against scattering. Uniform rates integrate their actual segment arc exactly;
//! bounded heterogeneous fields use Poisson majorant ratio tracking. No density
//! smoothness is inferred and no scattering rate is disguised as absorption.
//!
//! For a fixed admitted path, E[product(1-beta_s/M)] = exp(-integral beta_s ds).
//! See Novak, Selle and Jarosz (2014), *Residual Ratio Tracking for Estimating
//! Attenuation in Participating Media*, <https://doi.org/10.1145/2661229.2661292>.
//! This conditional estimator does not certify the numerical trajectory, global
//! optical-image completeness, or convergence through a refused work limit.
//! Replay retains native step/query budgets, contact collars, interfaces and
//! source identities; its final receiver and ordered branch are revalidated.

use super::*;
use crate::dielectric::DielectricEvent;
use crate::polarization::dielectric::scene::curved::bundle::CurvedBundleSceneError;
use crate::polarization::dielectric::scene::curved::bundle::connection::{
    BundleBranchPolicy, BundleConnection, BundleConnectionError, BundleConnectionLaunch,
    BundleConnectionOptions, BundleConnectionTarget,
};
use core::cell::RefCell;

/// Distinct optical images with independently tracked molecular survival.
pub mod images;

/// Cumulative work limits for one source-leg replay. Zero budgets permit
/// uniform/clear paths; refusal occurs only before an actual draw or candidate.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CurvedTransmittanceOptions {
    pub max_free_flight_draws: u32,
    pub max_tracking_candidates: u32,
}

impl Default for CurvedTransmittanceOptions {
    fn default() -> Self {
        Self {
            max_free_flight_draws: MAX_POLARIZED_PATH_EVENTS,
            max_tracking_candidates: MAX_POLARIZED_PATH_EVENTS,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum CurvedTransmittanceError {
    Connection(BundleConnectionError),
    Scene(CurvedBundleSceneError),
    Volume(VolumeTransportError),
    InvalidOptions,
    DrawLimit,
    CandidateLimit,
    ReplayBranchChanged,
    NumericalRange,
    Cancelled,
}

impl From<BundleConnectionError> for CurvedTransmittanceError {
    fn from(error: BundleConnectionError) -> Self {
        Self::Connection(error)
    }
}
impl From<CurvedBundleSceneError> for CurvedTransmittanceError {
    fn from(error: CurvedBundleSceneError) -> Self {
        Self::Scene(error)
    }
}
impl From<VolumeTransportError> for CurvedTransmittanceError {
    fn from(error: VolumeTransportError) -> Self {
        Self::Volume(error)
    }
}
impl From<fs_exec::Cancelled> for CurvedTransmittanceError {
    fn from(_: fs_exec::Cancelled) -> Self {
        Self::Cancelled
    }
}
impl core::fmt::Display for CurvedTransmittanceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "curved source scattering survival: {self:?}")
    }
}
impl core::error::Error for CurvedTransmittanceError {}

/// Original solved connection and separate source-leg scattering survival.
/// The log preserves very small positive weights. Negative infinity occurs
/// ONLY for an exact ratio-tracking zero; it is not an underflow convention.
/// Work includes every deterministic shooting trial and the tracked replay.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CurvedTransmittanceConnection {
    connection: BundleConnection,
    log_scattering_survival: f64,
    ratio_tracking_zero: bool,
    free_flight_draws: u32,
    tracking_candidates: u32,
    trace_evaluations: u32,
    field_evaluations: u64,
    query_work: u64,
}

impl CurvedTransmittanceConnection {
    pub const fn connection(&self) -> &BundleConnection {
        &self.connection
    }
    pub const fn log_scattering_survival(&self) -> f64 {
        self.log_scattering_survival
    }
    pub const fn ratio_tracking_zero(&self) -> bool {
        self.ratio_tracking_zero
    }
    pub const fn free_flight_draws(&self) -> u32 {
        self.free_flight_draws
    }
    pub const fn tracking_candidates(&self) -> u32 {
        self.tracking_candidates
    }
    pub const fn trace_evaluations(&self) -> u32 {
        self.trace_evaluations
    }
    /// Original optical/jet evaluations, including the replay. Local scattering
    /// queries are counted separately by `tracking_candidates`.
    pub const fn field_evaluations(&self) -> u64 {
        self.field_evaluations
    }
    pub const fn query_work(&self) -> u64 {
        self.query_work
    }
}

impl CurvedRayleighScene<'_, '_, '_> {
    /// Independent counter-based source flights. Camera surface and volume
    /// domains are untouched; callers choose a source-connection path identity.
    #[allow(clippy::too_many_arguments)]
    pub fn connect_bundle_with_transmittance(
        &self,
        cx: &Cx<'_>,
        launch: BundleConnectionLaunch,
        target: BundleConnectionTarget,
        connection: BundleConnectionOptions,
        policy: BundleBranchPolicy<'_>,
        tracking: CurvedTransmittanceOptions,
        seed: u64,
        path_id: u64,
    ) -> Result<CurvedTransmittanceConnection, CurvedTransmittanceError> {
        self.connect_bundle_with_transmittance_samples(
            cx,
            launch,
            target,
            connection,
            policy,
            tracking,
            |draw| {
                let block = philox4x32_10(
                    [path_id as u32, (path_id >> 32) as u32, draw, 0x706f_6c35],
                    [seed as u32, (seed >> 32) as u32],
                );
                (f64::from(block[0]) + 0.5) / 4_294_967_296.0
            },
        )
    }

    /// Solve first without stochastic subdivisions, then replay the converged
    /// launch once through the original bundle event owner. Each heterogeneous
    /// flight draws one independent u in [0,1); interrupted flights still count.
    /// `connection.max_trace_evaluations` includes the replay, so one trace is
    /// reserved before Newton begins. Entirely absent/zero-density media retain
    /// the original solve with no replay, new draws, or changed result bits.
    #[allow(clippy::too_many_arguments)]
    pub fn connect_bundle_with_transmittance_samples(
        &self,
        cx: &Cx<'_>,
        launch: BundleConnectionLaunch,
        target: BundleConnectionTarget,
        options: BundleConnectionOptions,
        policy: BundleBranchPolicy<'_>,
        tracking: CurvedTransmittanceOptions,
        samples: impl FnMut(u32) -> f64,
    ) -> Result<CurvedTransmittanceConnection, CurvedTransmittanceError> {
        cx.checkpoint()?;
        options.validate()?;
        let clear = self.scattering_known_clear(cx, launch.wavelength_nm)?;
        if !clear && options.max_trace_evaluations < 2 {
            return Err(BundleConnectionError::TraceBudgetExceeded.into());
        }
        let solve_options = BundleConnectionOptions {
            max_trace_evaluations: options.max_trace_evaluations - u32::from(!clear),
            ..options
        };
        let connection =
            self.scene
                .connect_bundle_with_policy(cx, launch, target, solve_options, policy)?;
        self.replay_connected_bundle_with_transmittance(
            cx, launch, target, options, policy, tracking, connection, clear, samples,
        )
    }

    // Only native solves owned by this exact scene call the replay below.
    // A public replay accepting an arbitrary connection would lose that source
    // and scene identity guarantee.
    fn scattering_known_clear(
        &self,
        cx: &Cx<'_>,
        wavelength_nm: f64,
    ) -> Result<bool, CurvedTransmittanceError> {
        for index in 0..self.scene.media.len() {
            if index % 1024 == 0 {
                cx.checkpoint()?;
            }
            let known_clear = match self.media.get(index) {
                None => true,
                Some(CurvedRayleighBinding::Density(law)) => {
                    matches!(law.density, CurvedDensityProfile::Uniform(density) if density == 0.0)
                        && law.scattering.at(wavelength_nm).is_ok()
                }
                Some(CurvedRayleighBinding::Field(_)) => false,
            };
            if !known_clear {
                // An unsupported spectrum is not vacuum. Defer its refusal to
                // the replay so an unrelated, unvisited medium is not queried
                // as though the connected branch had traversed it.
                return Ok(false);
            }
        }
        Ok(true)
    }

    #[allow(clippy::too_many_arguments)]
    fn replay_connected_bundle_with_transmittance(
        &self,
        cx: &Cx<'_>,
        launch: BundleConnectionLaunch,
        target: BundleConnectionTarget,
        options: BundleConnectionOptions,
        policy: BundleBranchPolicy<'_>,
        tracking: CurvedTransmittanceOptions,
        connection: BundleConnection,
        clear: bool,
        samples: impl FnMut(u32) -> f64,
    ) -> Result<CurvedTransmittanceConnection, CurvedTransmittanceError> {
        let mut output = CurvedTransmittanceConnection {
            connection,
            log_scattering_survival: 0.0,
            ratio_tracking_zero: false,
            free_flight_draws: 0,
            tracking_candidates: 0,
            trace_evaluations: connection.trace_evaluations(),
            field_evaluations: connection.field_evaluations(),
            query_work: connection.query_work(),
        };
        if clear {
            cx.checkpoint()?;
            return Ok(output);
        }
        let ray = connection.path().bundle().propagation().initial_ray();
        let bundle = self.scene.launch_bundle(
            ray,
            launch.transverse_axis,
            launch.wavelength_nm,
            launch.normalized_source_stokes,
            launch.initial_medium,
        )?;
        let state = RefCell::new(Tracking::new(self.media, tracking, samples));
        let mut consumed = 0_usize;
        let replay = self.scene.trace_bundle_with_receiver_flights(
            cx,
            bundle,
            launch.initial_medium,
            options.scene,
            true,
            target.native_receiver(),
            |_, context| match policy {
                BundleBranchPolicy::TransmitOrReflectTir => {
                    Ok(if context.total_internal_reflection {
                        DielectricEvent::Reflection
                    } else {
                        DielectricEvent::Transmission
                    })
                }
                BundleBranchPolicy::Sequence(events) => {
                    let event = events.get(consumed).copied().ok_or(
                        CurvedBundleSceneError::BranchPlanExhausted {
                            dielectric_event: consumed,
                        },
                    )?;
                    consumed += 1;
                    Ok(event)
                }
            },
            |medium, wavelength| {
                let mut state = state.borrow_mut();
                match state.flight(medium, wavelength) {
                    Ok(value) => Ok(value),
                    Err(error) => {
                        state.error = Some(error);
                        Err(CurvedBundleSceneError::NumericalRange)
                    }
                }
            },
            |medium, wavelength, distance, position, candidate| {
                let mut state = state.borrow_mut();
                match state.observe(medium, wavelength, distance, position, candidate) {
                    Ok(()) => Ok(()),
                    Err(error) => {
                        state.error = Some(error);
                        Err(CurvedBundleSceneError::NumericalRange)
                    }
                }
            },
        );
        let state = state.into_inner();
        if let Some(error) = state.error {
            return Err(error);
        }
        let replay = replay?;
        if let BundleBranchPolicy::Sequence(events) = policy
            && (replay.termination() == super::super::PathTermination::Detected
                || replay.is_diffuse_receiver()
                || replay.is_volume_receiver())
            && consumed != events.len()
        {
            return Err(BundleConnectionError::BranchPlanNotConsumed {
                planned: events.len(),
                consumed,
            }
            .into());
        }
        connection
            .validate_replay(replay, target, options.projection)
            .map_err(|error| match error {
                BundleConnectionError::BranchChanged => {
                    CurvedTransmittanceError::ReplayBranchChanged
                }
                other => CurvedTransmittanceError::Connection(other),
            })?;
        output.log_scattering_survival = state.log_survival;
        output.ratio_tracking_zero = state.ratio_zero;
        output.free_flight_draws = state.draws;
        output.tracking_candidates = state.candidates;
        output.trace_evaluations = output
            .trace_evaluations
            .checked_add(1)
            .ok_or(CurvedTransmittanceError::NumericalRange)?;
        output.field_evaluations = output
            .field_evaluations
            .checked_add(replay.field_evaluations())
            .ok_or(CurvedTransmittanceError::NumericalRange)?;
        output.query_work = output
            .query_work
            .checked_add(replay.query_work())
            .ok_or(CurvedTransmittanceError::NumericalRange)?;
        cx.checkpoint()?;
        Ok(output)
    }
}

#[derive(Clone, Copy)]
enum TrackingRate<'a> {
    Clear,
    Uniform(f64),
    Density {
        majorant: f64,
        field: &'a dyn CurvedDensityField,
        density_bound: f64,
    },
    Field {
        majorant: f64,
        field: &'a dyn CurvedRayleighField,
    },
}

struct Tracking<'a, 'optics, S> {
    media: CurvedRayleighBindings<'a, 'optics>,
    options: CurvedTransmittanceOptions,
    samples: S,
    cached: Option<(usize, f64, TrackingRate<'optics>)>,
    log_survival: f64,
    ratio_zero: bool,
    draws: u32,
    candidates: u32,
    error: Option<CurvedTransmittanceError>,
}

impl<'a, 'optics, S: FnMut(u32) -> f64> Tracking<'a, 'optics, S> {
    fn new(
        media: CurvedRayleighBindings<'a, 'optics>,
        options: CurvedTransmittanceOptions,
        samples: S,
    ) -> Self {
        Self {
            media,
            options,
            samples,
            cached: None,
            log_survival: 0.0,
            ratio_zero: false,
            draws: 0,
            candidates: 0,
            error: None,
        }
    }

    fn rate(
        &mut self,
        medium: usize,
        wavelength: f64,
    ) -> Result<TrackingRate<'optics>, CurvedTransmittanceError> {
        if let Some((index, lambda, rate)) = self.cached
            && index == medium
            && lambda == wavelength
        {
            return Ok(rate);
        }
        let rate = match self.media.get(medium) {
            None => TrackingRate::Clear,
            Some(CurvedRayleighBinding::Density(law)) => {
                let mass = law.scattering.at(wavelength)?;
                let density = law.density_kg_m3();
                let majorant = mass * density;
                if !majorant.is_finite() || (mass > 0.0 && density > 0.0 && majorant == 0.0) {
                    return Err(VolumeTransportError::Unrepresentable.into());
                }
                if majorant == 0.0 {
                    TrackingRate::Clear
                } else {
                    match law.density {
                        CurvedDensityProfile::Uniform(_) => TrackingRate::Uniform(majorant),
                        CurvedDensityProfile::Field {
                            field,
                            majorant_kg_m3,
                            ..
                        } => TrackingRate::Density {
                            majorant,
                            field,
                            density_bound: majorant_kg_m3,
                        },
                    }
                }
            }
            Some(CurvedRayleighBinding::Field(field)) => {
                let majorant = field.scattering_majorant_per_m(wavelength)?;
                if !majorant.is_finite() || majorant < 0.0 {
                    return Err(VolumeTransportError::ScatteringSpectrum.into());
                }
                if majorant == 0.0 {
                    TrackingRate::Clear
                } else {
                    TrackingRate::Field { majorant, field }
                }
            }
        };
        self.cached = Some((medium, wavelength, rate));
        Ok(rate)
    }

    fn flight(
        &mut self,
        medium: usize,
        wavelength: f64,
    ) -> Result<Option<f64>, CurvedTransmittanceError> {
        if self.ratio_zero {
            return Ok(None);
        }
        let majorant = match self.rate(medium, wavelength)? {
            TrackingRate::Clear | TrackingRate::Uniform(_) => return Ok(None),
            TrackingRate::Density { majorant, .. } | TrackingRate::Field { majorant, .. } => {
                majorant
            }
        };
        if self.draws >= self.options.max_free_flight_draws {
            return Err(CurvedTransmittanceError::DrawLimit);
        }
        let sample = (self.samples)(self.draws);
        self.draws += 1;
        if !sample.is_finite() || !(0.0..1.0).contains(&sample) {
            return Err(VolumeTransportError::Sample.into());
        }
        let complement = 1.0 - sample;
        let depth = if complement == 1.0 {
            sample
        } else {
            -det::ln(complement) * (sample / (1.0 - complement))
        };
        let distance = depth / majorant;
        if distance.is_nan() || distance < 0.0 || (depth > 0.0 && distance == 0.0) {
            return Err(VolumeTransportError::Unrepresentable.into());
        }
        Ok(Some(distance))
    }

    fn observe(
        &mut self,
        medium: usize,
        wavelength: f64,
        distance: f64,
        position: Point3,
        candidate: bool,
    ) -> Result<(), CurvedTransmittanceError> {
        if self.ratio_zero {
            return Ok(());
        }
        let rate = self.rate(medium, wavelength)?;
        if let TrackingRate::Uniform(rate) = rate {
            let depth = rate * distance;
            if !depth.is_finite() || depth < 0.0 || (rate > 0.0 && distance > 0.0 && depth == 0.0) {
                return Err(VolumeTransportError::Unrepresentable.into());
            }
            return self.add_log(-depth);
        }
        if !candidate {
            return Ok(());
        }
        if self.candidates >= self.options.max_tracking_candidates {
            return Err(CurvedTransmittanceError::CandidateLimit);
        }
        self.candidates += 1;
        let ratio = match rate {
            TrackingRate::Density {
                field,
                density_bound,
                ..
            } => {
                let density = field.density_kg_m3(position)?;
                if !density.is_finite() || density < 0.0 || density > density_bound {
                    return Err(VolumeTransportError::Density { index: 0 }.into());
                }
                let ratio = density / density_bound;
                if density > 0.0 && ratio == 0.0 {
                    return Err(VolumeTransportError::Unrepresentable.into());
                }
                ratio
            }
            TrackingRate::Field { majorant, field } => {
                let local = field.sample_scattering(position, wavelength)?;
                if !local.scattering_per_m.is_finite()
                    || local.scattering_per_m < 0.0
                    || local.scattering_per_m > majorant
                {
                    return Err(VolumeTransportError::ScatteringSpectrum.into());
                }
                PhaseFunction::rayleigh_depolarized(local.depolarization_ratio)?;
                let ratio = local.scattering_per_m / majorant;
                if local.scattering_per_m > 0.0 && ratio == 0.0 {
                    return Err(VolumeTransportError::Unrepresentable.into());
                }
                ratio
            }
            TrackingRate::Clear | TrackingRate::Uniform(_) => {
                return Err(CurvedTransmittanceError::NumericalRange);
            }
        };
        if ratio == 1.0 {
            self.ratio_zero = true;
            self.log_survival = f64::NEG_INFINITY;
            return Ok(());
        }
        // Compensate 1-ratio near one, retaining a small positive extinction
        // instead of rounding its log survival to zero.
        let complement = 1.0 - ratio;
        let logarithm = if complement == 1.0 {
            -ratio
        } else {
            det::ln(complement) * (ratio / (1.0 - complement))
        };
        self.add_log(logarithm)
    }

    fn add_log(&mut self, value: f64) -> Result<(), CurvedTransmittanceError> {
        let next = self.log_survival + value;
        if !value.is_finite() || value > 0.0 || !next.is_finite() || next > 0.0 {
            return Err(CurvedTransmittanceError::NumericalRange);
        }
        self.log_survival = next;
        Ok(())
    }
}

#[cfg(test)]
mod tests;
