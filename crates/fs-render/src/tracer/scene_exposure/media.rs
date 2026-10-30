//! Accepted interior and ambient field owners selected with their scene frame.

use super::*;

impl<'a> SceneExposure<'a> {
    /// Bind one already-admitted physical scene owner to every supplied frame.
    ///
    /// The owner must borrow that frame's exact scene. This prevents primitive
    /// IDs from accidentally binding density/temperature fields to another
    /// geometry state. Each sample holds the selected scene and its field
    /// owners for the complete optical path, with the same shutter time used
    /// by the camera. No state interpolation or field copying occurs.
    ///
    /// These owners include their own optional ambient field, so a separately
    /// attached [`MediumTimeline`] is refused in either attachment order.
    /// [`SceneMedia`] currently admits static world-space geometry/fields;
    /// successive accepted frames may change either, but continuous motion of
    /// an interior grid is not supplied by this method.
    pub fn with_scene_media(
        mut self,
        cx: &Cx<'_>,
        media: &'a [&'a SceneMedia<'a, 'a>],
    ) -> Result<Self, TracerError> {
        cx.checkpoint()?;
        if self.medium.is_some() {
            return Err(SceneExposureError::ConflictingMediumOwners.into());
        }
        if media.len() != self.frames.len() {
            return Err(SceneExposureError::MediumFrameCount {
                expected: self.frames.len(),
                actual: media.len(),
            }
            .into());
        }
        let mut hash = DomainHasher::new("org.frankensim.render.scene-exposure-media.v1");
        hash.update(&(media.len() as u64).to_le_bytes());
        for (index, (frame, owner)) in self.frames.iter().zip(media).enumerate() {
            cx.checkpoint()?;
            if !core::ptr::eq(frame.scene, owner.scene()) {
                return Err(SceneExposureError::MediumSceneMismatch { index }.into());
            }
            hash.update(owner.content_identity().as_bytes());
        }
        // SceneMedia already prepared this exact scene's source proposals.
        // Release duplicate exposure metadata once; selection later borrows
        // the owner's lighting without allocation or hierarchy rebuilding.
        for lighting in &mut self.lighting {
            cx.checkpoint()?;
            *lighting = None;
        }
        self.scene_media_identity = Some(hash.finalize());
        self.scene_media = Some(media);
        Ok(self)
    }

    pub(in crate::tracer) fn media_at_time(
        &self,
        absolute_time_s: f64,
    ) -> Result<Option<&'a SceneMedia<'a, 'a>>, TracerError> {
        if !absolute_time_s.is_finite()
            || absolute_time_s < self.shutter.open_s()
            || absolute_time_s > self.shutter.close_s()
        {
            return Err(SceneExposureError::TimeOutsideExposure.into());
        }
        Ok(self
            .scene_media
            .map(|media| media[self.frame_index(absolute_time_s)]))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::motion::{ShutterConvention, ShutterDistribution};
    use crate::volumes::VolumeGrid;
    use crate::volumes::transport::timeline::{MediumFrame, MediumTimelineBudget, PhysicalMedium};
    use crate::volumes::transport::{AbsorbingGrid, AbsorbingGridBudget, MassAbsorptionSpectrum};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

    fn with_cx<R>(operation: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
        let gate = CancelGate::new_clock_free();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            operation(
                &Cx::new(
                    &gate,
                    arena,
                    StreamKey {
                        seed: 53,
                        kernel_id: 1,
                        tile: 0,
                        iteration: 0,
                    },
                    Budget::INFINITE,
                    ExecMode::Deterministic,
                ),
                &gate,
            )
        })
    }

    fn scene() -> Scene {
        Scene {
            primitives: Vec::new(),
            lights: Vec::new(),
            environment: None,
            camera: Camera {
                eye: Point3::new(0.0, 0.0, 1.0),
                forward: Vec3::new(0.0, 0.0, -1.0),
                up: Vec3::new(0.0, 1.0, 0.0),
                half_tan: 0.0,
            },
        }
    }

    fn shutter() -> ShutterInterval {
        ShutterInterval::try_from_canonical_parts(
            0.0,
            2.0,
            ShutterConvention::FrontLoaded,
            ShutterDistribution::UniformCounterV1,
        )
        .unwrap()
    }

    fn exposure_budget() -> SceneExposureBudget {
        SceneExposureBudget {
            max_frames: 2,
            max_total_primitives: 0,
            max_total_light_triangles: 0,
        }
    }

    fn media_budget() -> SceneMediaBudget {
        SceneMediaBudget {
            max_interior_bindings: 0,
            max_primitives: 0,
            max_light_triangles: 0,
        }
    }

    fn field<'a>(cx: &Cx<'_>, density: &'a [f64]) -> AbsorbingGrid<'a> {
        AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1; 3], density, [-1.0; 3], [2.0; 3]),
            MassAbsorptionSpectrum::constant(0.5).unwrap(),
            None,
            AbsorbingGridBudget {
                max_grid_cells: 1,
                max_steps_per_segment: 16,
            },
        )
        .unwrap()
    }

    #[test]
    fn g0_frame_count_and_scene_ownership_must_match() {
        with_cx(|cx, _| {
            let scene = scene();
            let other_scene = self::scene();
            let frames = [
                SceneFrame {
                    start_s: 0.0,
                    end_s: 1.0,
                    scene: &scene,
                },
                SceneFrame {
                    start_s: 1.0,
                    end_s: 2.0,
                    scene: &scene,
                },
            ];
            let correct = SceneMedia::try_new(cx, &scene, None, &[], media_budget()).unwrap();
            let incorrect =
                SceneMedia::try_new(cx, &other_scene, None, &[], media_budget()).unwrap();
            let exposure =
                || SceneExposure::try_new(cx, &frames, shutter(), exposure_budget()).unwrap();
            assert_eq!(
                exposure().with_scene_media(cx, &[&correct]).err(),
                Some(
                    SceneExposureError::MediumFrameCount {
                        expected: 2,
                        actual: 1
                    }
                    .into()
                )
            );
            assert_eq!(
                exposure()
                    .with_scene_media(cx, &[&correct, &incorrect])
                    .err(),
                Some(SceneExposureError::MediumSceneMismatch { index: 1 }.into())
            );
        });
    }

    #[test]
    fn g0_scene_fields_and_prepared_lights_share_exact_accepted_state_selection() {
        with_cx(|cx, _| {
            let scene = scene();
            let density = [1.0];
            let field = field(cx, &density);
            let owners = [
                SceneMedia::try_new(cx, &scene, None, &[], media_budget()).unwrap(),
                SceneMedia::try_new(
                    cx,
                    &scene,
                    Some(PhysicalMedium::Absorbing(&field)),
                    &[],
                    media_budget(),
                )
                .unwrap(),
            ];
            let bindings = [&owners[0], &owners[1]];
            let frames = [
                SceneFrame {
                    start_s: 0.0,
                    end_s: 1.0,
                    scene: &scene,
                },
                SceneFrame {
                    start_s: 1.0,
                    end_s: 2.0,
                    scene: &scene,
                },
            ];
            let exposure = SceneExposure::try_new(cx, &frames, shutter(), exposure_budget())
                .unwrap()
                .with_scene_media(cx, &bindings)
                .unwrap();
            assert!(exposure.lighting.iter().all(Option::is_none));
            for (time, index) in [(0.0, 0), (0.75, 0), (1.0, 1), (2.0, 1)] {
                assert!(core::ptr::eq(
                    exposure.media_at_time(time).unwrap().unwrap(),
                    &owners[index]
                ));
                assert!(core::ptr::eq(
                    exposure.at_time(time).unwrap().1,
                    &owners[index].lighting
                ));
            }
            for time in [-0.1, 2.1, f64::NAN] {
                assert_eq!(
                    exposure.media_at_time(time).err(),
                    Some(SceneExposureError::TimeOutsideExposure.into())
                );
            }
        });
    }

    #[test]
    fn g0_ambient_history_and_frame_owners_refuse_ambiguous_ownership_in_either_order() {
        with_cx(|cx, _| {
            let scene = scene();
            let density = [1.0];
            let field = field(cx, &density);
            let physical_frames = [MediumFrame {
                start_s: 0.0,
                end_s: 2.0,
                medium: PhysicalMedium::Absorbing(&field),
            }];
            let history = MediumTimeline::try_new(
                cx,
                &physical_frames,
                MediumTimelineBudget { max_frames: 1 },
            )
            .unwrap();
            let frames = [SceneFrame {
                start_s: 0.0,
                end_s: 2.0,
                scene: &scene,
            }];
            let owner = SceneMedia::try_new(cx, &scene, None, &[], media_budget()).unwrap();
            let bindings = [&owner];
            let exposure =
                || SceneExposure::try_new(cx, &frames, shutter(), exposure_budget()).unwrap();
            assert_eq!(
                exposure()
                    .with_medium_history(cx, &history)
                    .unwrap()
                    .with_scene_media(cx, &bindings)
                    .err(),
                Some(SceneExposureError::ConflictingMediumOwners.into())
            );
            assert_eq!(
                exposure()
                    .with_scene_media(cx, &bindings)
                    .unwrap()
                    .with_medium_history(cx, &history)
                    .err(),
                Some(SceneExposureError::ConflictingMediumOwners.into())
            );
        });
    }

    #[test]
    fn g3_progressive_binding_uses_numeric_field_state_instead_of_borrowed_addresses() {
        with_cx(|cx, _| {
            let identity = |density: f64| {
                let scene = scene();
                let densities = [density];
                let field = field(cx, &densities);
                let owner = SceneMedia::try_new(
                    cx,
                    &scene,
                    Some(PhysicalMedium::Absorbing(&field)),
                    &[],
                    media_budget(),
                )
                .unwrap();
                let bindings = [&owner];
                let frames = [SceneFrame {
                    start_s: 0.0,
                    end_s: 2.0,
                    scene: &scene,
                }];
                SceneExposure::try_new(cx, &frames, shutter(), exposure_budget())
                    .unwrap()
                    .with_scene_media(cx, &bindings)
                    .unwrap()
                    .binding_identity()
            };
            assert_eq!(identity(1.0), identity(1.0));
            assert_ne!(identity(1.0), identity(1.1));
        });
    }

    #[test]
    fn g4_cancelled_attachment_publishes_no_accepted_owner_history() {
        with_cx(|cx, gate| {
            let scene = scene();
            let owner = SceneMedia::try_new(cx, &scene, None, &[], media_budget()).unwrap();
            let bindings = [&owner];
            let frames = [SceneFrame {
                start_s: 0.0,
                end_s: 2.0,
                scene: &scene,
            }];
            let exposure =
                SceneExposure::try_new(cx, &frames, shutter(), exposure_budget()).unwrap();
            gate.request();
            assert_eq!(
                exposure.with_scene_media(cx, &bindings).err(),
                Some(TracerError::Cancelled)
            );
        });
    }
}
