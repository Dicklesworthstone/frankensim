//! Native chart departures for straight camera, light and shadow paths.
//!
//! The shaded point stays physical. Only a spawned ray is moved, and only a
//! used hemisphere observes a geometric refusal. Optional normal-neighborhood
//! certificates strengthen the historical chart path without pretending that
//! an unsupported custom chart supplied a self-separation theorem.

use super::*;
use crate::charts::proximity::{ProximityBudget, ProximityError, box_clearance};
use crate::instances::proximity::PreparedNativeProximity;
use fs_ivl::Interval;

const BUILD_WORK: u64 = 64 * 1024 * 1024;
const DEPARTURE_WORK: u64 = 16 * 1024;

#[derive(Clone, Copy)]
pub(super) struct ChartOrigin {
    outward: Vec3,
    positive: Result<Point3, TracerError>,
    negative: Result<Point3, TracerError>,
}

impl ChartOrigin {
    pub(super) fn spawn(self, direction: Vec3) -> Result<Point3, TracerError> {
        // Opaque callers may pass a face-forwarded shading normal. Stored
        // hemispheres always refer to the native outward geometric normal.
        if self.outward.dot(direction) >= 0.0 {
            self.positive
        } else {
            self.negative
        }
    }
}

pub(super) struct Prepared<'scene> {
    scene: &'scene Scene,
    geometry: Vec<PreparedNativeProximity<'scene>>,
}

pub(super) fn has_chart(shape: &Shape) -> bool {
    match shape {
        Shape::Chart(_) => true,
        Shape::Instance(instance) => matches!(instance.geometry(), SharedGeometry::Chart(_)),
        Shape::AnimatedInstance(instance) => matches!(instance.geometry(), SharedGeometry::Chart(_)),
        Shape::Mesh(_) => false,
    }
}

impl<'scene> Prepared<'scene> {
    pub(super) fn for_scene(
        scene: &'scene Scene,
        cx: &Cx<'_>,
    ) -> Result<Option<Self>, TracerError> {
        if !scene.primitives.iter().any(|primitive| has_chart(&primitive.shape)) {
            return Ok(None);
        }
        let mut geometry = Vec::new();
        geometry.try_reserve_exact(scene.primitives.len()).map_err(|_| TracerError::InvalidInput)?;
        let mut budget = ProximityBudget::new(BUILD_WORK);
        let mut meshes = std::collections::BTreeMap::new();
        for primitive in &scene.primitives {
            budget.charge(cx).map_err(error)?;
            let mesh = match &primitive.shape {
                Shape::Mesh(mesh) => Some(mesh),
                Shape::Instance(instance) => match instance.geometry() {
                    SharedGeometry::Mesh(mesh) => Some(mesh.as_ref()),
                    SharedGeometry::Chart(_) => None,
                },
                Shape::AnimatedInstance(instance) => match instance.geometry() {
                    SharedGeometry::Mesh(mesh) => Some(mesh.as_ref()),
                    SharedGeometry::Chart(_) => None,
                },
                Shape::Chart(_) => None,
            };
            let view = if let Some(mesh) = mesh {
                let key = std::ptr::from_ref(mesh) as usize;
                if let Some(view) = meshes.get(&key) {
                    *view
                } else {
                    let view = PreparedNativeProximity::mesh(mesh, cx, &mut budget).map_err(error)?;
                    meshes.insert(key, view);
                    view
                }
            } else {
                match &primitive.shape {
                    Shape::Chart(chart) => PreparedNativeProximity::chart(chart.as_ref()),
                    Shape::Instance(instance) => PreparedNativeProximity::shared(
                        instance.geometry(), cx, &mut budget,
                    ).map_err(error)?,
                    Shape::AnimatedInstance(instance) => PreparedNativeProximity::shared(
                        instance.geometry(), cx, &mut budget,
                    ).map_err(error)?,
                    Shape::Mesh(_) => return Err(TracerError::InvalidInput),
                }
            };
            let view = if let Shape::Instance(instance) = &primitive.shape {
                view.at_transform(instance.transform()).map_err(error)?
            } else {
                view
            };
            geometry.push(view);
        }
        Ok(Some(Self { scene, geometry }))
    }

    fn placed(
        &self,
        index: usize,
        cx: &Cx<'_>,
        time: Option<&PathTime>,
    ) -> Result<PreparedNativeProximity<'scene>, TracerError> {
        let view = *self.geometry.get(index).ok_or(TracerError::InvalidInput)?;
        let Shape::AnimatedInstance(animated) = &self.scene.primitives[index].shape else {
            return Ok(view);
        };
        let time = time.ok_or(TracerError::MissingRayTime)?;
        let transform = if let Some(cached) = time.cached_animated.iter().flatten()
            .find(|cached| cached.primitive_index == index)
        {
            cached.instance.transform()
        } else {
            let ray = TimedRay::at_normalized(
                Ray { origin: Point3::new(0.0, 0.0, 0.0), dir: Vec3::new(0.0, 0.0, 1.0) },
                time.interval, time.normalized,
            );
            animated.instance_at(cx, ray.absolute_time_s())?.transform()
        };
        view.at_transform(transform).map_err(error)
    }

    fn clear_other_surfaces(
        &self,
        index: usize,
        cx: &Cx<'_>,
        time: Option<&PathTime>,
        origin: Point3,
        radius: f64,
        budget: &mut ProximityBudget,
    ) -> Result<(), TracerError> {
        // The known-size box pass is separately bounded and cancellable. A
        // large scene of distant instances must not consume the native-query
        // allowance simply because those cached boxes have to be inspected.
        let mut boxes = ProximityBudget::new(self.geometry.len() as u64);
        for other in 0..self.geometry.len() {
            if other == index {
                continue;
            }
            boxes.charge(cx).map_err(error)?;
            let geometry = self.placed(other, cx, time)?;
            if geometry.bounds().is_some_and(|bounds|
                box_clearance(origin, bounds).is_ok_and(|clearance| clearance > radius))
            {
                continue;
            }
            if geometry.query(cx, origin, budget).map_err(error)?.clearance_m <= radius {
                return Err(error(ProximityError::AmbiguousContact));
            }
        }
        Ok(())
    }

    fn side(
        &self,
        index: usize,
        cx: &Cx<'_>,
        time: Option<&PathTime>,
        point: Point3,
        normal: Vec3,
        side: f64,
        native: PreparedNativeProximity<'_>,
        signed: Option<Interval>,
        certified_neighborhood: bool,
    ) -> Result<Point3, TracerError> {
        let mut budget = ProximityBudget::new(DEPARTURE_WORK);
        let tolerance = native.trace_tolerance(TRACE_EPS);
        let uncertainty = signed.map_or(tolerance, |value| value.abs().hi());
        let direction = normal.scale(side);
        let origin = if certified_neighborhood {
            let value = signed.ok_or(error(ProximityError::UnsupportedDeparture))?;
            let deficit = if side > 0.0 { tolerance - value.lo() } else { tolerance + value.hi() };
            let mut offset = (Interval::point(deficit.max(tolerance * 0.25))
                * Interval::point(1.25)).hi();
            let mut accepted = None;
            for _ in 0..4 {
                let candidate = MeshBoundaryOrigin { point, normal_error: offset }.spawn(normal, direction);
                let value = native.signed_distance(cx, candidate, &mut budget).map_err(error)?;
                let native_value = native.native_trace_signed_distance(cx, candidate, &mut budget)
                    .map_err(error)?;
                if clear_side(value, tolerance, side)
                    && native_value.is_none_or(|value| clear_side(value, TRACE_EPS, side))
                {
                    if !monotone_departure(native, cx, point, candidate, side, &mut budget)? {
                        return Err(error(ProximityError::UnsupportedDeparture));
                    }
                    accepted = Some(candidate);
                    break;
                }
                offset = (Interval::point(offset) * Interval::point(2.0)).hi();
            }
            accepted.ok_or(error(ProximityError::UnsupportedDeparture))?
        } else {
            // Compatibility path for a chart that supplied no neighborhood
            // certificate at this hit. Keep its original arithmetic and
            // resolution; the other-surface exclusion below is still required.
            // This path does not certify uniqueness of its own local sheet.
            let candidate = dielectric_spawn_origin(point, normal, direction);
            if signed.is_some()
                && (!clear_side(native.signed_distance(cx, candidate, &mut budget).map_err(error)?,
                    tolerance, side)
                    || native.native_trace_signed_distance(cx, candidate, &mut budget).map_err(error)?
                        .is_some_and(|value| !clear_side(value, TRACE_EPS, side)))
            {
                return Err(error(ProximityError::UnsupportedDeparture));
            }
            candidate
        };
        let (center, radius) = segment_ball(point, origin, Interval::new(0.0, 1.0));
        let radius = (Interval::point(radius) + Interval::point(uncertainty)).hi();
        if !radius.is_finite() {
            return Err(error(ProximityError::NumericalRange));
        }
        // This one-sided ball encloses the actual skipped segment and the
        // original hit uncertainty. A nearby surface behind the hit need not
        // invalidate a reflection leaving on the opposite side.
        self.clear_other_surfaces(index, cx, time, center, radius, &mut budget)?;
        Ok(origin)
    }

    fn origin(
        &self,
        intersection: SceneIntersection,
        cx: &Cx<'_>,
        time: Option<&PathTime>,
    ) -> Result<ChartOrigin, TracerError> {
        let index = intersection.primitive_index;
        let point = intersection.hit.point;
        let normal = intersection.hit.normal.ok_or(TracerError::MissingNormal)?;
        let native = self.placed(index, cx, time)?;
        if !native.is_chart() {
            return Err(TracerError::InvalidInput);
        }
        let mut budget = ProximityBudget::new(DEPARTURE_WORK);
        let signed = native.is_exact_distance()
            .then(|| native.signed_distance(cx, point, &mut budget).map_err(error)).transpose()?;
        let certified_neighborhood = native.normal_cone(
            cx, point, native.trace_tolerance(TRACE_EPS) / 16.0, &mut budget,
        ).map_err(error)?.is_some();
        let positive = self.side(index, cx, time, point, normal, 1.0, native, signed, certified_neighborhood);
        if matches!(positive, Err(TracerError::Cancelled)) {
            return Err(TracerError::Cancelled);
        }
        let negative = self.side(index, cx, time, point, normal, -1.0, native, signed, certified_neighborhood);
        if matches!(negative, Err(TracerError::Cancelled)) {
            return Err(TracerError::Cancelled);
        }
        Ok(ChartOrigin { outward: normal, positive, negative })
    }
}

pub(super) fn prepare(
    mut intersection: SceneIntersection,
    cache: Option<&Prepared<'_>>,
    cx: &Cx<'_>,
    time: Option<&PathTime>,
) -> Result<SceneIntersection, TracerError> {
    if !matches!(intersection.boundary_origin, Some(BoundaryOrigin::UnpreparedChart)) {
        return Ok(intersection);
    }
    let result = cache.ok_or(TracerError::InvalidInput)
        .and_then(|cache| cache.origin(intersection, cx, time));
    intersection.boundary_origin = Some(BoundaryOrigin::Chart(match result {
        Ok(origin) => origin,
        Err(TracerError::Cancelled) => return Err(TracerError::Cancelled),
        Err(failure) => ChartOrigin {
            outward: intersection.hit.normal.unwrap_or(Vec3::new(0.0, 0.0, 0.0)),
            positive: Err(failure),
            negative: Err(failure),
        },
    }));
    Ok(intersection)
}

fn error(error: ProximityError) -> TracerError {
    if error == ProximityError::Cancelled {
        TracerError::Cancelled
    } else {
        TracerError::BoundaryDeparture(error)
    }
}

fn clear_side(value: Interval, tolerance: f64, side: f64) -> bool {
    if side > 0.0 { value.lo() > tolerance } else { value.hi() < -tolerance }
}

fn displacement(a: Point3, b: Point3) -> [Interval; 3] {
    let a = [a.x, a.y, a.z];
    let b = [b.x, b.y, b.z];
    std::array::from_fn(|i| if a[i] == b[i] {
        Interval::point(0.0)
    } else {
        Interval::point(b[i]) - Interval::point(a[i])
    })
}

fn length(vector: [Interval; 3]) -> Interval {
    let v = vector.map(Interval::abs);
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

fn segment_ball(start: Point3, end: Point3, parameter: Interval) -> (Point3, f64) {
    let delta = displacement(start, end);
    let start = [start.x, start.y, start.z];
    let segment: [Interval; 3] = std::array::from_fn(|axis| {
        if delta[axis].lo() == 0.0 && delta[axis].hi() == 0.0 {
            Interval::point(start[axis])
        } else {
            Interval::point(start[axis]) + delta[axis] * parameter
        }
    });
    let center = segment.map(Interval::midpoint);
    let radius = length(std::array::from_fn(|i| {
        if segment[i].lo() == center[i] && segment[i].hi() == center[i] {
            Interval::point(0.0)
        } else {
            segment[i] - Interval::point(center[i])
        }
    })).hi();
    (Point3::new(center[0], center[1], center[2]), radius)
}

fn monotone_departure(
    native: PreparedNativeProximity<'_>,
    cx: &Cx<'_>,
    start: Point3,
    end: Point3,
    side: f64,
    budget: &mut ProximityBudget,
) -> Result<bool, TracerError> {
    let delta = displacement(start, end);
    let chord_length = length(delta);
    // Dyadic parameter cells cover the original exact chord, not a polyline
    // through rounded midpoints. Ordinary surfaces need only the first cone.
    for level in 0..=5 {
        let cells = 1_u32 << level;
        let mut proved = true;
        for cell in 0..cells {
            let parameter = Interval::new(f64::from(cell) / f64::from(cells),
                f64::from(cell + 1) / f64::from(cells));
            let (center, radius) = segment_ball(start, end, parameter);
            let cone = native.normal_cone(cx, center, radius, budget).map_err(error)?;
            let Some(cone) = cone else {
                proved = false;
                continue;
            };
            let axis = [cone.axis.x, cone.axis.y, cone.axis.z];
            let slope = (0..3).fold(Interval::point(0.0), |sum, i|
                sum + Interval::point(side * axis[i]) * delta[i]);
            if slope.lo() <= (Interval::point(cone.max_deviation) * chord_length).hi() {
                proved = false;
            }
        }
        if proved {
            return Ok(true);
        }
    }
    Ok(false)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::animated_instances::{RigidTransformTrajectory, TransformKeyframe};
    use crate::motion::{ShotTimeBounds, ShutterConvention, ShutterDistribution};
    use crate::spectral::lift_rgb;
    use fs_blake3::hash_domain;
    use fs_exec::{CancelGate, StreamKey};
    use fs_geom::fixtures::SphereChart;

    fn with_cx(operation: impl FnOnce(&Cx<'_>)) {
        let gate = CancelGate::new();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(&gate, arena, StreamKey {
                seed: 73, kernel_id: 1, tile: 0, iteration: 0,
            }, Budget::INFINITE, ExecMode::Deterministic);
            operation(&cx);
        });
    }

    fn quad(x: f64, z: f64) -> TriMesh {
        TriMesh::new(vec![
            [x - 1.0, -1.0, z], [x + 1.0, -1.0, z],
            [x + 1.0, 1.0, z], [x - 1.0, 1.0, z],
        ], vec![[0, 1, 2], [0, 2, 3]])
    }

    fn fixture(placement: u32, x: f64, gap: f64, side: f64) -> (Scene, Option<PathTime>) {
        let sphere = SphereChart { center: Point3::new(0.0, 0.0, 0.0), radius: 1.0 };
        let transform = |x| RigidTransform::try_new(
            [0.0, 0.0, 0.0, 1.0], [x, 0.0, -1.0],
        ).unwrap();
        let identity = hash_domain("chart-departure-test", b"sphere");
        let shape = match placement {
            0 => Shape::Chart(Box::new(SphereChart {
                center: Point3::new(x, 0.0, -1.0), ..sphere
            })),
            1 => Shape::Instance(GeometryInstance::try_new(
                1, identity, SharedGeometry::chart(sphere), transform(x),
            ).unwrap()),
            _ => Shape::AnimatedInstance(AnimatedGeometryInstance::try_new(
                1, identity, SharedGeometry::chart(sphere),
                RigidTransformTrajectory::try_new(vec![
                    TransformKeyframe::try_new(0.0, transform(x - 1.0), [2.0, 0.0, 0.0]).unwrap(),
                    TransformKeyframe::try_new(1.0, transform(x + 1.0), [2.0, 0.0, 0.0]).unwrap(),
                ]).unwrap(),
            ).unwrap()),
        };
        let time = (placement == 2).then(|| PathTime {
            interval: ShutterInterval::resolve(0.5, 1.0, ShutterConvention::Centered,
                ShutterDistribution::UniformCounterV1, ShotTimeBounds::try_new(0.0, 1.0).unwrap(),
            ).unwrap(),
            normalized: NormalizedShutterTime::try_new(0.5).unwrap(),
            cached_animated: std::array::from_fn(|_| None),
        });
        let white = lift_rgb([1.0; 3]);
        (Scene {
            primitives: vec![
                Primitive {
                    shape,
                    material: Material::Conductor {
                        optics: ConductorOptics::representative_tungsten(),
                        surface: ConductorSurface::SMOOTH,
                    },
                    emission: None,
                },
                Primitive {
                    shape: Shape::Mesh(quad(x, side * gap)),
                    material: Material::SpectralDiffuse {
                        reflectance: crate::spectral::SpectralReflectance::constant(0.0).unwrap(),
                    },
                    emission: Some((white, 3.0)),
                },
            ],
            lights: vec![RectLight {
                corner: Point3::new(x - 1.0, -1.0, side * gap),
                edge_u: Vec3::new(2.0, 0.0, 0.0), edge_v: Vec3::new(0.0, 2.0, 0.0),
                prim: 1, emission: (white, 3.0),
            }], environment: None,
            camera: Camera {
                eye: Point3::new(x, 0.0, side * gap * 0.9),
                forward: Vec3::new(0.0, 0.0, -side),
                up: Vec3::new(0.0, 1.0, 0.0), half_tan: 0.0,
            },
        }, time)
    }

    #[test]
    fn chart_reflections_retain_submicron_sources_and_native_sides() {
        with_cx(|cx| {
            let settings = Settings { width: 1, height: 1, spp: 1, max_depth: 2,
                sampler: Sampler::Iid, strategy: DirectStrategy::BsdfOnly, seed: 73 };
            let wavelengths = hero_wavelength_packet::<PACKET>(
                0.5 * (LAMBDA_MIN + LAMBDA_MAX), LAMBDA_MIN, LAMBDA_MAX,
            );
            let optics = ConductorOptics::representative_tungsten();
            let expected = packet_radiance_to_xyz(wavelengths.map(|lambda|
                3.0 * lift_rgb([1.0; 3]).eval(lambda) * optics.fresnel(lambda, 1.0, 1.0).unwrap()),
                &wavelengths, 1.0 / y_integral());
            for placement in 0..3 {
                for side in [-1.0, 1.0] {
                    // The first gap also rejects the unnecessarily large ball
                    // centered at the departure endpoint instead of its midpoint.
                    for gap in [1.5e-7, 5.0e-7, 0.1] {
                        let (scene, mut time) = fixture(placement, 3.0, gap, side);
                        if let Some(time) = &mut time {
                            let Shape::AnimatedInstance(animated) = &scene.primitives[0].shape else {
                                unreachable!();
                            };
                            if gap == 5.0e-7 {
                                time.cached_animated[0] = Some(CachedAnimatedInstance {
                                    primitive_index: 0, instance: animated.instance_at(cx, 0.5).unwrap(),
                                });
                            }
                        }
                        let lighting = validate_scene(&scene, cx, time.as_ref().map(|t| t.interval)).unwrap();
                        let ray = Ray { origin: scene.camera.eye, dir: scene.camera.forward };
                        let hit = intersect(&scene, cx, &ray, time.as_ref()).unwrap().unwrap();
                        assert_eq!(hit.primitive_index, 0);
                        let physical = hit.hit.point;
                        let hit = prepare(hit, lighting.chart_departure.as_ref(), cx, time.as_ref()).unwrap();
                        assert_eq!(hit.hit.point, physical);
                        let outward = hit.hit.normal.unwrap();
                        let outgoing = outward.scale(side);
                        // Inside reflection and outward transmission choose the
                        // native hemisphere even when the caller flips its normal.
                        let origin = hit.spawn_boundary_ray(outward.scale(-1.0), outgoing).unwrap();
                        assert_eq!(origin, hit.spawn_boundary_ray(outward, outgoing).unwrap());
                        let continuation = Ray { origin, dir: outgoing };
                        assert_eq!(intersect(&scene, cx, &continuation, time.as_ref()).unwrap()
                            .unwrap().primitive_index, 1);
                        let actual = trace_path(&scene, &lighting, cx, &settings, 1.0 / y_integral(),
                            0, 0, 0.5, 0.5, 0.5, time.as_ref(), CameraPath::Legacy,
                            false, None, None).unwrap().xyz;
                        for channel in 0..3 {
                            assert!((actual[channel] - expected[channel]).abs()
                                <= 3.0e-12 * expected[channel].abs().max(1.0),
                                "placement={placement}, side={side}, gap={gap}, actual={actual:?}");
                        }
                    }
                }
            }
            // Large translation parallel to the departure must not invent a
            // normal thickness or make a resolvable nearby reflector disappear.
            let (scene, _) = fixture(0, 1.0e8, 5.0e-7, 1.0);
            let lighting = validate_scene(&scene, cx, None).unwrap();
            let actual = trace_path(&scene, &lighting, cx, &settings, 1.0 / y_integral(),
                0, 0, 0.5, 0.5, 0.5, None, CameraPath::Legacy, false, None, None).unwrap().xyz;
            assert!(actual.iter().all(|value| *value > 0.0));
        });
    }

    #[test]
    fn chart_departure_refusals_are_lazy_and_shadow_origins_do_not_skip_blockers() {
        with_cx(|cx| {
            let (mut scene, _) = fixture(0, 0.0, 5.0e-7, 1.0);
            scene.primitives[1].emission = None;
            scene.lights.clear();
            scene.environment = Some(EnvironmentMap::try_from_linear_srgb(
                2, 2, vec![[1.0; 3]; 4], 0.0,
            ).unwrap());
            scene.primitives.push(Primitive {
                shape: Shape::Mesh(quad(0.0, -0.5 * TRACE_EPS)),
                material: scene.primitives[1].material.clone(), emission: None,
            });
            let lighting = validate_scene(&scene, cx, None).unwrap();
            let ray = Ray { origin: scene.camera.eye, dir: scene.camera.forward };
            let hit = intersect(&scene, cx, &ray, None).unwrap().unwrap();
            let hit = prepare(hit, lighting.chart_departure.as_ref(), cx, None).unwrap();
            let normal = hit.hit.normal.unwrap();
            let origin = hit.spawn_boundary_ray(normal, normal).unwrap();
            let shadow = Ray { origin, dir: normal };
            assert_eq!(intersect(&scene, cx, &shadow, None).unwrap().unwrap().primitive_index, 1);
            assert!(matches!(hit.spawn_boundary_ray(normal, normal.scale(-1.0)),
                Err(TracerError::BoundaryDeparture(ProximityError::AmbiguousContact))));
            // An unused unresolved side does not prevent terminal source emission.
            drop(lighting);
            scene.primitives[0].emission = Some((lift_rgb([1.0; 3]), 2.0));
            let settings = Settings { width: 1, height: 1, spp: 1, max_depth: 1,
                sampler: Sampler::Iid, strategy: DirectStrategy::BsdfOnly, seed: 73 };
            assert!(render(&scene, cx, &settings).unwrap().xyz[0].iter().all(|v| *v > 0.0));
        });
    }

    struct WithoutNeighborhood(SphereChart);
    impl Chart for WithoutNeighborhood {
        fn name(&self) -> &'static str { "sphere-without-neighborhood" }
        fn eval(&self, point: Point3, cx: &Cx<'_>) -> fs_geom::ChartSample {
            self.0.eval(point, cx)
        }
        fn support(&self) -> fs_geom::Aabb { self.0.support() }
        fn trace_step_claim(&self) -> fs_geom::TraceStepClaim {
            fs_geom::TraceStepClaim::ExactDistance
        }
    }

    #[test]
    fn unsupported_chart_neighborhood_preserves_legacy_resolution_with_clearance_guard() {
        with_cx(|cx| {
            for gap in [5.0e-7, 0.1] {
                let (mut scene, _) = fixture(0, 0.0, gap, 1.0);
                scene.primitives[0].shape = Shape::Chart(Box::new(WithoutNeighborhood(
                    SphereChart { center: Point3::new(0.0, 0.0, -1.0), radius: 1.0 },
                )));
                let lighting = validate_scene(&scene, cx, None).unwrap();
                let ray = Ray { origin: scene.camera.eye, dir: scene.camera.forward };
                let hit = intersect(&scene, cx, &ray, None).unwrap().unwrap();
                let hit = prepare(hit, lighting.chart_departure.as_ref(), cx, None).unwrap();
                let normal = hit.hit.normal.unwrap();
                if gap < RAY_EPS {
                    assert!(matches!(hit.spawn_boundary_ray(normal, normal),
                        Err(TracerError::BoundaryDeparture(ProximityError::AmbiguousContact))));
                } else {
                    assert_eq!(hit.spawn_boundary_ray(normal, normal).unwrap(),
                        dielectric_spawn_origin(hit.hit.point, normal, normal));
                }
            }
        });
    }

    #[test]
    fn rotated_chart_departures_clear_the_actual_native_inverse_arithmetic() {
        with_cx(|cx| {
            let transform = RigidTransform::try_new(
                [0.0, 0.3826834323650898, 0.0, 0.9238795325112867],
                [1.0e8, -1.0e8, 1.0e8],
            ).unwrap();
            let sphere = SphereChart { center: Point3::new(0.0, 0.0, 0.0), radius: 1.0 };
            let instance = GeometryInstance::try_new(1, hash_domain("chart-departure-test", b"large"),
                SharedGeometry::chart(sphere), transform).unwrap();
            let direction = transform.transform_vector(Vec3::new(0.0, 0.0, 1.0));
            let scene = Scene {
                primitives: vec![Primitive { shape: Shape::Instance(instance),
                    material: Material::Conductor {
                        optics: ConductorOptics::representative_tungsten(), surface: ConductorSurface::SMOOTH,
                    }, emission: None }],
                lights: Vec::new(), environment: None,
                camera: Camera { eye: transform.transform_point(Point3::new(0.0, 0.0, 2.0)),
                    forward: direction.scale(-1.0), up: Vec3::new(0.0, 1.0, 0.0), half_tan: 0.0 },
            };
            let cache = Prepared::for_scene(&scene, cx).unwrap();
            let ray = Ray { origin: scene.camera.eye, dir: scene.camera.forward };
            let hit = intersect(&scene, cx, &ray, None).unwrap().unwrap();
            let hit = prepare(hit, cache.as_ref(), cx, None).unwrap();
            let normal = hit.hit.normal.unwrap();
            for side in [-1.0, 1.0] {
                let origin = hit.spawn_boundary_ray(normal, normal.scale(side)).unwrap();
                let local = transform.inverse().unwrap().transform_point(origin);
                let value = sphere.eval(local, cx);
                let enclosure = sphere.trace_value_enclosure(local, &value, cx);
                assert!(if side > 0.0 { enclosure.lo > TRACE_EPS } else { enclosure.hi < -TRACE_EPS });
                let continuation = Ray { origin, dir: normal.scale(side) };
                let next = intersect(&scene, cx, &continuation, None).unwrap();
                if side > 0.0 {
                    assert!(next.is_none());
                } else {
                    assert!(next.unwrap().hit.t > 1.0);
                }
            }
        });
    }
}
