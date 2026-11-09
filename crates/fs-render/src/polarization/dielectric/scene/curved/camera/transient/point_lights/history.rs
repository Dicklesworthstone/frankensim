//! Source-owned point histories and the original shared retarded-clock law.

use super::super::exposure::{checked_exposure_close, piecewise_linear::PiecewiseLinearHistory};
use super::super::{ConnectionClockError, checked_retarded_time, connected_group_delay};
use super::*;
use crate::photometry::source::PhotometricSource;
use fs_blake3::ContentHash;
use fs_geom::{Point3, Vec3};

/// One exact calibrated source on its actual emitted connection direction.
/// Identity belongs to the authored immutable history; the source reference
/// retains the original calibration and placement without a fabricated object ID.
#[derive(Clone, Copy, Debug)]
pub struct PointEmissionQuery<'a> {
    pub source: &'a PhotometricSource,
    pub source_identity: ContentHash,
    pub position_m: Point3,
    pub outgoing_direction: Vec3,
    pub wavelength_nm: f64,
    pub observation_time_s: Dd,
    pub emission_time_s: Dd,
    pub group_delay_s: Dd,
}

#[derive(Clone, Copy, Debug)]
pub struct PointExposureQuery<'a> {
    pub opening: PointEmissionQuery<'a>,
    pub observation_close_s: Dd,
    pub emission_close_s: Dd,
    pub duration_s: Dd,
}

/// Immutable intensity multiplier on the original calibrated fixture. This
/// changes neither its spectral/polarization state nor any proposal or survival
/// probability. A finite shutter requires a source-owned interval mean.
pub trait PointEmissionHistory: Sync {
    fn sample(
        &self,
        cx: &Cx<'_>,
        query: PointEmissionQuery<'_>,
    ) -> Result<f64, SurfaceHistoryError>;
    fn average(
        &self,
        cx: &Cx<'_>,
        _query: PointExposureQuery<'_>,
    ) -> Result<f64, SurfaceHistoryError> {
        cx.checkpoint()?;
        Err(SurfaceHistoryError::ExposureUnsupported)
    }
}

/// Bind the exact original fixture, not a clone or equal-valued calibration.
#[derive(Clone, Copy)]
pub struct PointEmissionHistoryBinding<'a> {
    pub source: &'a PhotometricSource,
    pub source_identity: ContentHash,
    pub history: &'a dyn PointEmissionHistory,
}

impl PointEmissionHistory for PiecewiseLinearHistory<'_> {
    fn sample(
        &self,
        cx: &Cx<'_>,
        query: PointEmissionQuery<'_>,
    ) -> Result<f64, SurfaceHistoryError> {
        self.value_at(cx, query.emission_time_s)
    }
    fn average(
        &self,
        cx: &Cx<'_>,
        query: PointExposureQuery<'_>,
    ) -> Result<f64, SurfaceHistoryError> {
        self.mean_between(cx, query.opening.emission_time_s, query.emission_close_s)
    }
}

pub(super) struct PointHistoryObserver<'a, 'h, 's, 'm, 'l> {
    pub(super) camera: &'a TransientCurvedCameraScene<'h, 's, 'm>,
    histories: &'a [PointEmissionHistoryBinding<'l>],
}

impl<'a, 'h, 's, 'm, 'l> PointHistoryObserver<'a, 'h, 's, 'm, 'l> {
    pub(super) fn try_new(
        cx: &Cx<'_>,
        camera: &'a TransientCurvedCameraScene<'h, 's, 'm>,
        histories: &'a [PointEmissionHistoryBinding<'l>],
        mut contains_source: impl FnMut(&PhotometricSource) -> bool,
    ) -> Result<Self, TransientPointLightCameraError> {
        cx.checkpoint().map_err(TransientCameraError::from)?;
        if histories.len() > super::super::super::point_lights::MAX_POINT_LIGHT_BRANCHES {
            return Err(TransientPointLightCameraError::InvalidHistoryBinding {
                index: super::super::super::point_lights::MAX_POINT_LIGHT_BRANCHES,
            });
        }
        for (index, binding) in histories.iter().enumerate() {
            cx.checkpoint().map_err(TransientCameraError::from)?;
            if binding.source_identity == ContentHash([0; 32]) {
                return Err(TransientPointLightCameraError::MissingHistoryIdentity { index });
            }
            if !contains_source(binding.source)
                || histories[..index]
                    .iter()
                    .any(|earlier| core::ptr::eq(earlier.source, binding.source))
            {
                return Err(TransientPointLightCameraError::InvalidHistoryBinding { index });
            }
        }
        Ok(Self { camera, histories })
    }

    pub(super) fn scale(
        &self,
        cx: &Cx<'_>,
        observation: Observation,
        launch_delay_s: Dd,
        path: PointSourcePath<'_>,
    ) -> Result<f64, TransientPointLightCameraError> {
        let (source, camera, connection) = match path {
            PointSourcePath::Surface(path) => {
                return Ok(match observation {
                    Observation::Instant(time) => {
                        self.camera
                            .source_scale_with_delay(cx, time, launch_delay_s, path)?
                    }
                    Observation::Exposure(window) => {
                        self.camera
                            .exposure_scale_with_delay(cx, window, launch_delay_s, path)?
                    }
                });
            }
            PointSourcePath::Connected {
                source,
                camera,
                connection,
            } => (source, camera, connection),
        };
        let Some(binding) = self
            .histories
            .iter()
            .find(|binding| core::ptr::eq(binding.source, source))
        else {
            return Ok(1.0);
        };
        cx.checkpoint().map_err(TransientCameraError::from)?;
        let source_identity = binding.source_identity;
        let bundle = connection.path().bundle();
        let group_delay_s = connected_group_delay(
            camera.group_optical_path_length_m(),
            bundle
                .propagation()
                .propagation()
                .group_optical_path_length_m,
            launch_delay_s,
        )
        .map_err(|error| match error {
            ConnectionClockError::MissingGroupDelay => {
                TransientPointLightCameraError::MissingGroupDelay { source_identity }
            }
            ConnectionClockError::UnresolvedClock => {
                TransientPointLightCameraError::UnresolvedClock { source_identity }
            }
        })?;
        let observation_time_s = observation.open_s();
        let emission_time_s = checked_retarded_time(observation_time_s, group_delay_s)
            .ok_or(TransientPointLightCameraError::UnresolvedClock { source_identity })?;
        let query = PointEmissionQuery {
            source,
            source_identity,
            position_m: source.position_m(),
            outgoing_direction: bundle.propagation().initial_ray().dir,
            wavelength_nm: camera.beam().wavelength_nm(),
            observation_time_s,
            emission_time_s,
            group_delay_s,
        };
        let scale = match observation {
            Observation::Instant(_) => binding.history.sample(cx, query),
            Observation::Exposure(window) => {
                let emission_close_s =
                    checked_exposure_close(emission_time_s, group_delay_s, window).ok_or(
                        TransientPointLightCameraError::UnresolvedClock { source_identity },
                    )?;
                binding.history.average(
                    cx,
                    PointExposureQuery {
                        opening: query,
                        observation_close_s: window.close_s(),
                        emission_close_s,
                        duration_s: window.duration_s(),
                    },
                )
            }
        }
        .map_err(|error| TransientPointLightCameraError::History {
            source_identity,
            error,
        })?;
        cx.checkpoint().map_err(TransientCameraError::from)?;
        if !scale.is_finite() || scale < 0.0 {
            return Err(TransientPointLightCameraError::InvalidHistoryScale { source_identity });
        }
        Ok(scale)
    }
}
