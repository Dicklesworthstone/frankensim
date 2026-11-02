//! Point-fixture next-event terms at original accepted molecular vertices.
//!
//! A transverse source bundle supplies W/(m² nm). The reciprocal Rayleigh
//! contraction supplies 1/sr, including both optical subpaths' polarization.
//! The sampled camera collision already supplies scattering density and its
//! survival probability; neither beta nor a phase-sampling PDF belongs here.

use super::super::{
    CameraRadianceError, CurvedPathResult, CurvedPointLightCameraScene, CurvedSceneError,
    FixtureConnectionRequest, FixtureConnector, FixtureContributions, MAX_POINT_LIGHT_BRANCHES,
    PointLightCameraError, nonnegative, phase_index,
};
use crate::charts::Ray;
use crate::photometry::source::PhotometricSource;
use crate::photometry::source::bundle::{BundlePhotometryError, PlanarIrradianceReceiver};
use crate::polarization::Stokes;
use crate::polarization::dielectric::scene::PathTermination;
use crate::polarization::dielectric::scene::curved::SceneMedium;
use crate::polarization::dielectric::scene::curved::bundle::connection::{
    BundleBranchPolicy, BundleConnectionLaunch, BundleConnectionTarget,
};
use crate::polarization::dielectric::transport::PolarizationFrame;
use crate::volumes::transport::scattering::PhaseFunction;
use fs_exec::Cx;
use fs_geom::Vec3;
use fs_math::det;

/// One original fixture and one prescribed source-to-volume optical branch.
/// Every accepted molecular collision in `receiver_medium` attempts this local
/// seed. The native solver must resolve the actual sampled point and medium;
/// the seed is not a completeness claim for other optical images. Each source
/// pointer/receiving-medium pair may appear only once.
#[derive(Clone, Copy, Debug)]
pub struct MolecularPointLightBranch<'a> {
    pub source: &'a PhotometricSource,
    pub receiver_medium: usize,
    pub initial_medium: usize,
    pub initial_direction: Vec3,
    pub source_transverse_axis: Vec3,
    pub policy: BundleBranchPolicy<'a>,
    pub position_tolerance_m: f64,
}

pub(super) fn validate_branches(
    cx: &Cx<'_>,
    owner: &CurvedPointLightCameraScene<'_, '_, '_, '_>,
    branches: &[MolecularPointLightBranch<'_>],
) -> Result<(), PointLightCameraError> {
    if branches.len() > MAX_POINT_LIGHT_BRANCHES.saturating_sub(owner.branches.len()) {
        return Err(PointLightCameraError::TooManyBindings);
    }
    for (slot, branch) in branches.iter().enumerate() {
        cx.checkpoint()?;
        let index = owner.branches.len() + slot;
        let invalid = || PointLightCameraError::InvalidBinding { index };
        for medium in [branch.initial_medium, branch.receiver_medium] {
            if !matches!(
                owner.camera.scene.media.get(medium),
                Some(SceneMedium::Homogeneous(_) | SceneMedium::Differentiable { .. })
            ) {
                return Err(invalid());
            }
        }
        if !branch.position_tolerance_m.is_finite() || branch.position_tolerance_m <= 0.0 {
            return Err(invalid());
        }
        PolarizationFrame::try_new(branch.initial_direction, branch.source_transverse_axis)
            .map_err(|_| invalid())?;
        for (first, earlier) in branches[..slot].iter().enumerate() {
            if earlier.receiver_medium == branch.receiver_medium
                && core::ptr::eq(earlier.source, branch.source)
            {
                return Err(PointLightCameraError::DuplicateBinding {
                    first: owner.branches.len() + first,
                    second: index,
                });
            }
        }
    }
    Ok(())
}

impl FixtureContributions<'_, '_, '_, '_, '_> {
    pub(super) fn observe_volume(
        &mut self,
        cx: &Cx<'_>,
        path: &CurvedPathResult,
        phase: &PhaseFunction,
        branches: &[MolecularPointLightBranch<'_>],
        connector: &mut impl FixtureConnector,
    ) -> Result<(), CurvedSceneError> {
        self.add_volume(cx, path, phase, branches, connector)
            .map_err(|error| {
                self.error = Some(error);
                CurvedSceneError::NumericalRange
            })
    }

    fn add_volume(
        &mut self,
        cx: &Cx<'_>,
        path: &CurvedPathResult,
        phase: &PhaseFunction,
        branches: &[MolecularPointLightBranch<'_>],
        connector: &mut impl FixtureConnector,
    ) -> Result<(), PointLightCameraError> {
        let beam = *path.beam();
        for (slot, branch) in branches.iter().enumerate() {
            if branch.receiver_medium != path.medium_index() || branch.source.flux_lm() == 0.0 {
                continue;
            }
            cx.checkpoint()?;
            let index = self.owner.branches.len() + slot;
            if self.work.connections >= self.settings.max_connections {
                return Err(PointLightCameraError::ConnectionLimit);
            }
            let remaining = self
                .settings
                .max_trace_evaluations
                .checked_sub(self.work.traces)
                .ok_or(PointLightCameraError::TraceEvaluationLimit)?;
            if remaining == 0 {
                return Err(PointLightCameraError::TraceEvaluationLimit);
            }
            let mut options = self.settings.connection;
            options.max_trace_evaluations = options.max_trace_evaluations.min(remaining);
            // This section only localizes the endpoint for Newton's two
            // transverse residuals. It does not define the physical receiving
            // area; that is perpendicular to the actual solved arrival below.
            let target = BundleConnectionTarget::try_volume_receiver(
                branch.receiver_medium,
                beam.ray().origin,
                branch.source.position_m().delta_from(beam.ray().origin),
                branch.position_tolerance_m,
            )
            .map_err(|error| PointLightCameraError::Connection {
                binding: index,
                error,
            })?;
            let request = FixtureConnectionRequest {
                binding: index,
                ordinal: self.work.connections,
                launch: BundleConnectionLaunch {
                    ray: Ray {
                        origin: branch.source.position_m(),
                        dir: branch.initial_direction,
                    },
                    transverse_axis: branch.source_transverse_axis,
                    wavelength_nm: beam.wavelength_nm(),
                    normalized_source_stokes: Stokes::unpolarized(1.0),
                    initial_medium: branch.initial_medium,
                },
                target,
                options,
                policy: branch.policy,
            };
            self.work.connections += 1;
            let connection = connector.connect(cx, request, self.work)?;
            self.work.add(connection.work)?;
            if connection.ratio_tracking_zero {
                if connection.log_scattering_survival != f64::NEG_INFINITY {
                    return Err(PointLightCameraError::NumericalRange);
                }
                continue;
            }
            if !connection.log_scattering_survival.is_finite()
                || connection.log_scattering_survival > 0.0
            {
                return Err(PointLightCameraError::NumericalRange);
            }
            let bundle = connection.original.path().bundle();
            let endpoint = bundle.propagation().propagation();
            let source_frame =
                PolarizationFrame::try_new(endpoint.ray.dir, endpoint.transverse_axis).map_err(
                    |error| PointLightCameraError::MolecularPhase {
                        binding: index,
                        error: error.into(),
                    },
                )?;
            let angular = phase
                .rayleigh_connection_factor(beam, source_frame, bundle.interface_state())
                .map_err(|error| PointLightCameraError::MolecularPhase {
                    binding: index,
                    error,
                })?;
            // A physical polarized phase can be exactly zero. The complete
            // source geometry has already been validated, so this is a valid
            // dark contribution, not a failed shooting attempt or retry.
            if angular == 0.0 {
                continue;
            }
            let receiver = PlanarIrradianceReceiver::try_new(
                target.position_m(),
                endpoint.ray.dir.scale(-1.0),
                target.position_tolerance_m(),
            )
            .map_err(|error| PointLightCameraError::Photometry {
                binding: index,
                error,
            })?;
            let irradiance = branch
                .source
                .spectral_irradiance_from_refractive_bundle(bundle, receiver, options.projection)
                .map_err(|error| PointLightCameraError::Photometry {
                    binding: index,
                    error,
                })?;
            if irradiance == 0.0 {
                if bundle.interface_log_power() != f64::NEG_INFINITY {
                    let intensity = branch
                        .source
                        .spectral_intensity_w_sr_nm(
                            bundle.propagation().initial_ray().dir,
                            beam.wavelength_nm(),
                        )
                        .map_err(|error| PointLightCameraError::Photometry {
                            binding: index,
                            error: BundlePhotometryError::from(error),
                        })?;
                    if intensity > 0.0 {
                        return Err(PointLightCameraError::NumericalRange);
                    }
                }
                continue;
            }
            if beam.power().i == 0.0 {
                return Err(CameraRadianceError::Unresolved(PathTermination::Extinguished).into());
            }
            let local_index = phase_index(beam.medium(), beam.wavelength_nm())?;
            let contribution = volume_contribution(
                beam.power().i,
                self.initial_index,
                local_index,
                angular,
                irradiance,
                connection.log_scattering_survival,
            )?;
            self.radiance = nonnegative(self.radiance + contribution)?;
        }
        Ok(())
    }
}

fn volume_contribution(
    prefix: f64,
    initial_index: f64,
    local_index: f64,
    angular_per_sr: f64,
    irradiance: f64,
    log_survival: f64,
) -> Result<f64, PointLightCameraError> {
    let index_ratio = initial_index / local_index;
    let value = prefix
        * (index_ratio * index_ratio)
        * (angular_per_sr * irradiance)
        * det::exp(log_survival);
    if value.is_finite() && value > 0.0 {
        return Ok(value);
    }
    // A tiny surviving source term may still have representable radiance once
    // all calibrated factors are combined. Exact phase/tracking zeros were
    // handled above and never enter this positive logarithmic product.
    let logarithm = det::ln(prefix)
        + 2.0 * (det::ln(initial_index) - det::ln(local_index))
        + det::ln(angular_per_sr)
        + det::ln(irradiance)
        + log_survival;
    let value = det::exp(logarithm);
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(PointLightCameraError::NumericalRange)
    }
}
