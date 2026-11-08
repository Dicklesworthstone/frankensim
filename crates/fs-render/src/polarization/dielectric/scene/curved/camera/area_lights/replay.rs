//! Evaluate directional thermal power on an already discovered native image.
//!
//! Image support and MIS stay geometric. A bundle retains one transported
//! Stokes state, so each nonzero thermal s/p source gets one fixed-ray replay
//! of the same branch. No Newton correction, second image search, angular
//! source fit, or reconstruction of a weak channel from rounded I/Q occurs.

use super::{
    AreaEmission, AreaEmitterError, AreaEmitterSample, CameraRadianceError, CurvedAreaEmitter,
};
use crate::gradient_index::event::BoundarySide;
use crate::polarization::Stokes;
use crate::polarization::dielectric::scene::curved::bundle::connection::{
    BundleBranchPolicy, BundleConnection, BundleConnectionError, BundleConnectionLaunch,
    BundleConnectionOptions, BundleConnectionTarget, BundleSurfaceSource,
};
use fs_exec::Cx;

#[derive(Default)]
pub(super) struct ImageEmission {
    pub(super) channels: [Option<(BundleConnection, f64)>; 2],
    pub(super) traces: u32,
    pub(super) fields: u64,
    pub(super) geometry: u64,
}

impl CurvedAreaEmitter<'_, '_, '_> {
    #[allow(clippy::too_many_arguments)]
    pub(super) fn image_emission(
        &self,
        cx: &Cx<'_>,
        sample: AreaEmitterSample<'_>,
        side: BoundarySide,
        launch: BundleConnectionLaunch,
        original: &BundleConnection,
        target: BundleConnectionTarget,
        options: BundleConnectionOptions,
        policy: BundleBranchPolicy<'_>,
        uniform_radiance: Option<f64>,
        remaining_traces: u32,
    ) -> Result<ImageEmission, AreaEmitterError> {
        if let Some(radiance) = uniform_radiance {
            return Ok(ImageEmission {
                channels: [Some((*original, radiance)), None],
                ..ImageEmission::default()
            });
        }
        let path = original.path();
        let contact = path
            .source_contact()
            .ok_or(AreaEmitterError::SourceContactMismatch)?;
        let initial = path.bundle().propagation().initial_ray();
        if !core::ptr::eq(sample.emitter, &self.emitter)
            || !core::ptr::eq(sample.instance, self.instance)
            || !core::ptr::eq(sample.support, &self.support)
            || !self.emitter.supports_side(side)
            || contact.object_id != self.emitter.object_id()
            || contact.geometry_identity != self.instance.geometry_identity()
            || contact.frame_identity != self.instance.frame_identity()
            || initial.origin != sample.position_m()
            || contact.witness.normal.dot(sample.normal()) < 1.0 - 1e-8
        {
            return Err(AreaEmitterError::SourceContactMismatch);
        }
        self.emitter
            .admit_phase_index(
                path.bundle()
                    .propagation()
                    .propagation()
                    .initial_sample
                    .phase_index(),
            )
            .map_err(AreaEmitterError::Camera)?;
        if let AreaEmission::Thermal(emitter) = &self.emitter {
            // EON emission is unpolarized but directional. Evaluate it on
            // the converged original source ray without another image search
            // or polarization replay. Proposal geometry and MIS are unchanged.
            let radiance = emitter
                .radiance_along(launch.wavelength_nm, contact.witness.normal, initial.dir)
                .map_err(CameraRadianceError::Thermal)
                .map_err(AreaEmitterError::Camera)?;
            return Ok(ImageEmission {
                channels: [Some((*original, radiance)), None],
                ..ImageEmission::default()
            });
        }
        let AreaEmission::Coating(emitter) = &self.emitter else {
            return Err(AreaEmitterError::InvalidEmitter {
                object_id: self.emitter.object_id(),
            });
        };
        // The coherent Planck/work solve runs once at the converged source
        // direction, never at the seed guess. Its frame is physical s/p.
        let (radiance, frame) = emitter
            .sample(
                launch.wavelength_nm,
                initial.dir,
                contact.witness.normal,
                launch.transverse_axis,
            )
            .map_err(AreaEmitterError::Camera)?;
        if radiance == [0.0; 2] {
            return Ok(ImageEmission::default());
        }
        if radiance[0] == radiance[1] {
            return Ok(ImageEmission {
                channels: [Some((*original, radiance[0] + radiance[1])), None],
                ..ImageEmission::default()
            });
        }
        let required = radiance.iter().filter(|&&value| value > 0.0).count() as u32;
        if required > remaining_traces {
            return Err(AreaEmitterError::ReplayTraceLimit {
                required,
                remaining: remaining_traces,
            });
        }
        let source = BundleSurfaceSource::try_new(emitter.object_id(), side)
            .map_err(AreaEmitterError::Source)?;
        let mut replay_options = options;
        // The solver must accept its FIRST native trace of this converged
        // ray. Its existing cap refuses any Newton/backtracking retrace.
        // This original connection route sets continue_extinguished=true:
        // exact-zero interface power reaches the same receiver with log power
        // -infinity. Generic trace errors and positive underflow still refuse.
        replay_options.max_trace_evaluations = 1;
        let mut result = ImageEmission::default();
        for (channel, radiance) in radiance.into_iter().enumerate() {
            if radiance == 0.0 {
                continue;
            }
            cx.checkpoint()?;
            let replay = self
                .scene
                .connect_bundle_from_surface_with_policy(
                    cx,
                    BundleConnectionLaunch {
                        ray: initial,
                        transverse_axis: frame.horizontal(),
                        normalized_source_stokes: Stokes {
                            i: 1.0,
                            q: if channel == 0 { 1.0 } else { -1.0 },
                            u: 0.0,
                            v: 0.0,
                        },
                        ..launch
                    },
                    source,
                    target,
                    replay_options,
                    policy,
                )
                .map_err(|error| AreaEmitterError::PolarizationReplay { channel, error })?;
            if replay.iterations() != 0 || replay.trace_evaluations() != 1 {
                return Err(AreaEmitterError::PolarizationReplay {
                    channel,
                    error: BundleConnectionError::BranchChanged,
                });
            }
            original
                .validate_replay(*replay.path(), target, options.projection)
                .map_err(|error| AreaEmitterError::PolarizationReplay { channel, error })?;
            result.traces = result
                .traces
                .checked_add(replay.trace_evaluations())
                .ok_or(AreaEmitterError::NumericalRange)?;
            result.fields = result
                .fields
                .checked_add(replay.field_evaluations())
                .ok_or(AreaEmitterError::NumericalRange)?;
            result.geometry = result
                .geometry
                .checked_add(replay.query_work())
                .ok_or(AreaEmitterError::NumericalRange)?;
            result.channels[channel] = Some((replay, radiance));
        }
        Ok(result)
    }
}
