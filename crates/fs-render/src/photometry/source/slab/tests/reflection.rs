use super::*;

fn precise() -> SlabHistoryOptions {
    SlabHistoryOptions {
        absolute_tolerance_w_m2_nm: 1e-16,
        relative_tolerance: 1e-11,
        ..SlabHistoryOptions::default()
    }
}

fn receiver(x: f64) -> IlluminanceProbe {
    IlluminanceProbe::try_new(Point3::new(x, 0.0, -0.6), Vec3::new(0.0, 0.0, 1.0)).unwrap()
}

fn face() -> Point3 {
    Point3::new(0.0, 0.0, 0.0)
}

#[test]
fn g1_reflected_axial_series_includes_front_mirror_internal_returns_and_retroreception() {
    scope(|cx, _| {
        let source = source([0.0, 0.0, -1.0], AngularDistribution::isotropic());
        // The selected reflected field is finite even when the endpoints
        // coincide. Direct point-source light is deliberately not included.
        let receiver =
            IlluminanceProbe::try_new(source.position_m(), Vec3::new(0.3, 0.4, 1.0)).unwrap();
        for n in [1.0_f64, 1.5, 2.4] {
            let scene = scene(CauchyIor::try_constant(n).unwrap(), 0.7);
            let result = source
                .reflected_slab_history_irradiance(
                    cx,
                    &scene,
                    0,
                    face(),
                    receiver,
                    TransportOptions::default(),
                    SlabHistoryOptions {
                        // A fixed per-bin spectral allowance can exceed the
                        // relative lux accuracy below for this dim reflection.
                        // Request the corresponding relative-only solve.
                        absolute_tolerance_w_m2_nm: 0.0,
                        ..precise()
                    },
                )
                .unwrap();
            assert_eq!(result.family(), SlabHistoryFamily::Reflected);
            let r = ((n - 1.0) / (n + 1.0)).powi(2);
            let front = r / 4.0;
            let reference = receiver.normal().z
                * (front
                    + (1..128)
                        .map(|k| {
                            let depth = 0.4 * f64::from(k);
                            (1.0 - r).powi(2) * r.powi(2 * k - 1) * (-0.7 * depth).exp()
                                / (2.0 + depth / n).powi(2)
                        })
                        .sum::<f64>());
            close(result.reading().illuminance_lux, reference, 2e-10);
            let omitted = reference - result.reading().illuminance_lux;
            assert!(omitted >= -2e-13 * reference);
            assert!(omitted <= result.omitted_reading_bound().illuminance_lux + 2e-13 * reference);
            let work: u64 = result
                .images_per_wavelength()
                .iter()
                .map(|n| u64::from(*n).pow(2))
                .sum();
            assert_eq!(result.reserved_interface_interactions(), work);
            for bin in 0..SPECTRAL_BINS {
                assert!(
                    result.omitted_spectral_irradiance_bound_w_m2_nm()[bin]
                        <= precise().absolute_tolerance_w_m2_nm
                            + precise().relative_tolerance
                                * result.spectral_irradiance_w_m2_nm()[bin]
                );
            }
            if n == 1.0 {
                assert_eq!(result.reading().illuminance_lux, 0.0);
                assert!(result.images_per_wavelength().iter().all(|n| *n == 0));
                assert_eq!(result.omitted_reading_bound().illuminance_lux, 0.0);
            } else {
                assert!(result.reading().illuminance_lux > receiver.normal().z * front);
            }
        }
    });
}

// Independent forward construction. No mirrored receiver, inverse Snell
// connector or production differential is used: reflect the actual z direction
// at each actual face, then intersect the actual receiving plane.
fn endpoint(direction: Vec3, n: f64, receiver: IlluminanceProbe, order: u32) -> Point3 {
    let mut point = Point3::new(0.0, 0.0, -1.0).offset(direction.scale(1.0 / direction.z));
    if order > 0 {
        let mut inside = Vec3::new(
            direction.x / n,
            direction.y / n,
            (1.0 - (direction.x * direction.x + direction.y * direction.y) / (n * n)).sqrt(),
        );
        for crossing in 0..2 * order {
            point = point.offset(inside.scale(0.2 / inside.z.abs()));
            if crossing + 1 != 2 * order {
                inside.z = -inside.z;
            }
        }
    }
    let outgoing = Vec3::new(direction.x, direction.y, -direction.z);
    let parameter = receiver
        .position_m()
        .delta_from(point)
        .dot(receiver.normal())
        / outgoing.dot(receiver.normal());
    point.offset(outgoing.scale(parameter))
}

fn area(direction: Vec3, n: f64, receiver: IlluminanceProbe, order: u32) -> f64 {
    let epsilon = 2e-5;
    let derivative = |axis: Vec3| {
        let plus = Vec3::new(
            direction.x + epsilon * axis.x,
            direction.y + epsilon * axis.y,
            direction.z + epsilon * axis.z,
        );
        let minus = Vec3::new(
            direction.x - epsilon * axis.x,
            direction.y - epsilon * axis.y,
            direction.z - epsilon * axis.z,
        );
        endpoint(plus.scale(1.0 / plus.norm()), n, receiver, order)
            .delta_from(endpoint(
                minus.scale(1.0 / minus.norm()),
                n,
                receiver,
                order,
            ))
            .scale(0.5 / epsilon)
    };
    let a = derivative(Vec3::new(direction.z, 0.0, -direction.x));
    let b = derivative(Vec3::new(0.0, 1.0, 0.0));
    Vec3::new(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    )
    .norm()
}

fn bare_power(direction: Vec3, n: f64, sigma: f64, order: u32) -> f64 {
    let cg = (1.0 - direction.x * direction.x / (n * n)).sqrt();
    let c = direction.z;
    let rs = ((c - n * cg) / (c + n * cg)).powi(2);
    let rp = ((n * c - cg) / (n * c + cg)).powi(2);
    if order == 0 {
        0.5 * (rs + rp)
    } else {
        let ts = 4.0 * n * c * cg / (c + n * cg).powi(2);
        let tp = 4.0 * n * c * cg / (n * c + cg).powi(2);
        0.5 * (ts * ts * rs.powi(2 * order as i32 - 1) + tp * tp * rp.powi(2 * order as i32 - 1))
            * (-sigma * 0.4 * f64::from(order) / cg).exp()
    }
}

#[test]
fn g1_reflected_images_match_folded_snell_differentials_and_polarized_power() {
    scope(|cx, _| {
        let law = CauchyIor::try_new(1.42, 0.018, 0.0).unwrap();
        let scene = scene(law, 0.3);
        let receiver =
            IlluminanceProbe::try_new(receiver(3.6).position_m(), Vec3::new(0.1, 0.2, 1.0))
                .unwrap();
        let connection = SlabConnection::try_reflected(
            cx,
            &scene,
            0,
            Point3::new(0.0, 0.0, -1.0),
            face(),
            receiver,
            TransportOptions::default(),
        )
        .unwrap();
        for bin in [0, 40, 79] {
            let lambda = wavelength_nm(bin);
            let n = law.eval(lambda).unwrap();
            let mut previous_cosine = 0.0;
            for order in [0, 1, 3] {
                let image = connection.image(cx, lambda, order).unwrap();
                let ray = connection.sample_image(cx, &image).unwrap();
                let d = ray.launch_direction;
                assert!(ray.unoccluded);
                assert!(d.z > previous_cosine);
                previous_cosine = d.z;
                assert!(
                    endpoint(d, n, receiver, order)
                        .delta_from(receiver.position_m())
                        .norm()
                        < 2e-11
                );
                close(ray.area_per_solid_angle, area(d, n, receiver, order), 8e-8);
                close(ray.power_transmission, bare_power(d, n, 0.3, order), 3e-12);
            }
        }
    });
}

#[test]
fn g3_hidden_front_mirror_and_blocked_internal_return_do_not_erase_later_light() {
    scope(|cx, _| {
        let mut scene = scene(CauchyIor::try_constant(1.5).unwrap(), 0.1);
        let receiver = receiver(3.6);
        let connection = SlabConnection::try_reflected(
            cx,
            &scene,
            0,
            Point3::new(0.0, 0.0, -1.0),
            face(),
            receiver,
            TransportOptions::default(),
        )
        .unwrap();
        let first = connection
            .sample_image(cx, &connection.image(cx, 550.0, 0).unwrap())
            .unwrap();
        let internal = connection
            .sample_image(cx, &connection.image(cx, 550.0, 1).unwrap())
            .unwrap();
        let cutoff = (0.5 * (first.launch_direction.z + internal.launch_direction.z)).acos();
        let source = source(
            [0.0, 0.0, -1.0],
            AngularDistribution::uniform_cone(cutoff).unwrap(),
        );
        assert_eq!(source.intensity_cd(first.launch_direction).unwrap(), 0.0);
        let clear = source
            .reflected_slab_history_irradiance(
                cx,
                &scene,
                0,
                face(),
                receiver,
                TransportOptions::default(),
                precise(),
            )
            .unwrap();
        assert!(clear.reading().illuminance_lux > 0.0);
        // Neither the face discovery segment nor the opposite-face normal
        // probe is an optical segment. Opaque objects on them do not shade it.
        scene.primitives.push(panel(0.0, -0.5, 0.001));
        scene.primitives.push(panel(0.0, 0.1, 0.001));
        let same = source
            .reflected_slab_history_irradiance(
                cx,
                &scene,
                0,
                face(),
                receiver,
                TransportOptions::default(),
                precise(),
            )
            .unwrap();
        assert_eq!(clear, same);
        let second_lux = source.intensity_cd(internal.launch_direction).unwrap()
            / internal.area_per_solid_angle
            * internal.power_transmission;
        let blocker_x = receiver.position_m().x
            - 0.3 * internal.launch_direction.x / internal.launch_direction.z;
        scene.primitives.push(panel(blocker_x, -0.3, 0.001));
        let blocked = source
            .reflected_slab_history_irradiance(
                cx,
                &scene,
                0,
                face(),
                receiver,
                TransportOptions::default(),
                precise(),
            )
            .unwrap();
        assert!(blocked.reading().illuminance_lux > 0.0);
        assert!(
            (clear.reading().illuminance_lux - blocked.reading().illuminance_lux - second_lux)
                .abs()
                <= clear.omitted_reading_bound().illuminance_lux
                    + blocked.omitted_reading_bound().illuminance_lux
                    + 2e-12 * second_lux
        );
    });
}

#[test]
fn g1_reflected_receiver_projection_can_select_only_later_images() {
    scope(|cx, _| {
        let source = source([0.0, 0.0, -1.0], AngularDistribution::isotropic());
        let scene = scene(CauchyIor::try_constant(1.5).unwrap(), 0.0);
        let original = receiver(3.6);
        let connection = SlabConnection::try_reflected(
            cx,
            &scene,
            0,
            source.position_m(),
            face(),
            original,
            TransportOptions::default(),
        )
        .unwrap();
        let d0 = connection
            .sample_image(cx, &connection.image(cx, 550.0, 0).unwrap())
            .unwrap()
            .launch_direction;
        let d1 = connection
            .sample_image(cx, &connection.image(cx, 550.0, 1).unwrap())
            .unwrap()
            .launch_direction;
        let slope = 0.5 * (d0.x / d0.z + d1.x / d1.z);
        let tilted =
            IlluminanceProbe::try_new(original.position_m(), Vec3::new(1.0, 0.0, slope)).unwrap();
        assert!(-tilted.normal().x * d0.x + tilted.normal().z * d0.z < 0.0);
        assert!(-tilted.normal().x * d1.x + tilted.normal().z * d1.z > 0.0);
        let lit = source
            .reflected_slab_history_irradiance(
                cx,
                &scene,
                0,
                face(),
                tilted,
                TransportOptions::default(),
                precise(),
            )
            .unwrap();
        assert!(lit.reading().illuminance_lux > 0.0);
        let grazing =
            IlluminanceProbe::try_new(original.position_m(), Vec3::new(0.0, 1.0, 0.0)).unwrap();
        let dark = source
            .reflected_slab_history_irradiance(
                cx,
                &scene,
                0,
                face(),
                grazing,
                TransportOptions::default(),
                precise(),
            )
            .unwrap();
        assert_eq!(dark.reading().illuminance_lux, 0.0);
        assert!(dark.images_per_wavelength().iter().all(|n| *n > 1));
    });
}

// Direct single-film amplitude summation, independent of the native coating
// transfer owner. Positive imaginary index uses exp(+i k z) attenuation.
fn film_power(incident: f64, transmitted: f64, wavelength_nm: f64) -> (f64, f64) {
    type C = (f64, f64);
    fn add(a: C, b: C) -> C {
        (a.0 + b.0, a.1 + b.1)
    }
    fn mul(a: C, b: C) -> C {
        (a.0 * b.0 - a.1 * b.1, a.0 * b.1 + a.1 * b.0)
    }
    fn div(a: C, b: C) -> C {
        let norm = b.0 * b.0 + b.1 * b.1;
        (
            (a.0 * b.0 + a.1 * b.1) / norm,
            (a.1 * b.0 - a.0 * b.1) / norm,
        )
    }
    fn squared(a: C) -> f64 {
        a.0 * a.0 + a.1 * a.1
    }
    let index = (2.0, 0.15);
    let r01 = div(
        (incident - index.0, -index.1),
        (incident + index.0, index.1),
    );
    let r12 = div(
        (index.0 - transmitted, index.1),
        (index.0 + transmitted, index.1),
    );
    let phase = 2.0 * core::f64::consts::PI * 80.0 / wavelength_nm;
    let e = (
        (phase * index.0).cos() * (-phase * index.1).exp(),
        (phase * index.0).sin() * (-phase * index.1).exp(),
    );
    let e2 = mul(e, e);
    let denominator = add((1.0, 0.0), mul(mul(r01, r12), e2));
    let reflection = div(add(r01, mul(r12, e2)), denominator);
    let t01 = div((2.0 * incident, 0.0), (incident + index.0, index.1));
    let t12 = div(
        (2.0 * index.0, 2.0 * index.1),
        (index.0 + transmitted, index.1),
    );
    let transmission = div(mul(mul(t01, t12), e), denominator);
    (
        squared(reflection),
        transmitted / incident * squared(transmission),
    )
}

#[test]
fn g1_absorbing_coatings_use_outer_reflection_and_reversed_internal_interfaces() {
    use crate::conductor::{
        CauchyCoatingLayer, CoherentPeriodicCoating, ConductorDataStatus, ConductorSource,
    };
    use crate::dielectric::DielectricCoating;
    scope(|cx, _| {
        let mut scene = scene(CauchyIor::try_constant(1.5).unwrap(), 12.0);
        let provenance = ConductorSource::try_new(
            fs_blake3::hash_domain("fs-render:slab-reflection-test", b"absorbing-film"),
            ConductorDataStatus::Representative,
        )
        .unwrap();
        let layer = CauchyCoatingLayer::try_new(CauchyIor::try_constant(2.0).unwrap(), 0.15, 80e-9)
            .unwrap();
        if let Material::Dielectric { surface, .. } = &mut scene.primitives[0].material {
            *surface = DielectricSurface::SMOOTH.with_coating(
                DielectricCoating::try_new(
                    CoherentPeriodicCoating::try_new(&[layer], 1, provenance).unwrap(),
                )
                .unwrap(),
            );
        }
        let source = source([0.0, 0.0, -1.0], AngularDistribution::isotropic());
        let result = source
            .reflected_slab_history_irradiance(
                cx,
                &scene,
                0,
                face(),
                receiver(0.0),
                TransportOptions::default(),
                precise(),
            )
            .unwrap();
        for bin in [0, 40, 79] {
            let lambda = wavelength_nm(bin);
            let (r1, t1) = film_power(1.0, 1.5, lambda);
            let (r2, t2) = film_power(1.5, 1.0, lambda);
            assert!((r1 - r2).abs() > 0.01);
            let reference = r1 / 1.6_f64.powi(2)
                + (1..64)
                    .map(|k| {
                        let depth = 0.4 * f64::from(k);
                        t1 * t2 * r2.powi(2 * k - 1) * (-12.0 * depth).exp()
                            / (1.6 + depth / 1.5).powi(2)
                    })
                    .sum::<f64>();
            let expected = reference
                * source
                    .spectral_intensity_w_sr_nm(Vec3::new(0.0, 0.0, 1.0), lambda)
                    .unwrap();
            let actual = result.spectral_irradiance_w_m2_nm()[bin];
            close(actual, expected, 2e-9);
            assert!(
                expected - actual
                    <= result.omitted_spectral_irradiance_bound_w_m2_nm()[bin] + 2e-12 * expected
            );
        }
    });
}

#[test]
fn g1_bright_reflected_internal_image_survives_normalized_beer_underflow() {
    scope(|cx, _| {
        let scene = scene(CauchyIor::try_constant(1.5).unwrap(), 2200.0);
        let receiver = receiver(3.6);
        let connection = SlabConnection::try_reflected(
            cx,
            &scene,
            0,
            Point3::new(0.0, 0.0, -1.0),
            face(),
            receiver,
            TransportOptions::default(),
        )
        .unwrap();
        let front = connection
            .sample_image(cx, &connection.image(cx, 550.0, 0).unwrap())
            .unwrap();
        let internal = connection
            .sample_image(cx, &connection.image(cx, 550.0, 1).unwrap())
            .unwrap();
        assert_eq!(internal.power_transmission, 0.0);
        let cutoff = (0.5 * (front.launch_direction.z + internal.launch_direction.z)).acos();
        let source = PhotometricSource::try_new(
            RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [0.0, 0.0, -1.0]).unwrap(),
            4.0 * core::f64::consts::PI * 1e295,
            1.0,
            AngularDistribution::uniform_cone(cutoff).unwrap(),
            SourceSpectrum::equal_energy(),
        )
        .unwrap();
        let result = source
            .reflected_slab_history_irradiance(
                cx,
                &scene,
                0,
                face(),
                receiver,
                TransportOptions::default(),
                SlabHistoryOptions {
                    absolute_tolerance_w_m2_nm: 1e-210,
                    relative_tolerance: 0.0,
                    ..precise()
                },
            )
            .unwrap();
        let d = internal.launch_direction;
        let geometric_power = bare_power(d, 1.5, 0.0, 1);
        let cg = (1.0 - d.x * d.x / 2.25).sqrt();
        for bin in 0..SPECTRAL_BINS {
            let intensity = source
                .spectral_intensity_w_sr_nm(d, wavelength_nm(bin))
                .unwrap();
            let expected = (intensity.ln() + geometric_power.ln()
                - 2200.0 * 0.4 / cg
                - internal.area_per_solid_angle.ln())
            .exp();
            assert!(expected > 1e-210);
            close(
                result.spectral_irradiance_w_m2_nm()[bin] / expected,
                1.0,
                3e-12,
            );
            assert_eq!(result.images_per_wavelength()[bin], 2);
        }
    });
}

#[test]
fn g0_grazing_front_bound_does_not_hide_an_unrepresentable_receiving_jacobian() {
    scope(|cx, _| {
        let scene = scene(CauchyIor::try_constant(1.5).unwrap(), 10000.0);
        let receiver =
            IlluminanceProbe::try_new(Point3::new(0.0, 0.0, -1.0), Vec3::new(1.0, 0.0, 1e-320))
                .unwrap();
        let connection = SlabConnection::try_reflected(
            cx,
            &scene,
            0,
            receiver.position_m(),
            face(),
            receiver,
            TransportOptions::default(),
        )
        .unwrap();
        let image = connection.image(cx, 550.0, 0).unwrap();
        let bound = connection.omitted_bound(cx, &image, 1e295).unwrap();
        let expected = (1e295 * receiver.normal().z) * 0.04 / 4.0;
        close(bound / expected, 1.0, 3e-12);
        assert!(bound > 1e-29);
        // The geometric sample retains its explicit finite-J admission; the
        // bound cannot erase a representable scaled front signal before that
        // refusal is reached. Internal histories here are Beer-negligible.
        assert!(matches!(
            connection.sample_image(cx, &image),
            Err(PhotometryError::NonFiniteResult)
        ));
    });
}

#[test]
fn g0_reflected_history_refuses_invalid_face_media_unfinished_work_and_cancellation() {
    scope(|cx, gate| {
        let source = source([0.0, 0.0, -1.0], AngularDistribution::isotropic());
        let scene = scene(CauchyIor::try_constant(1.5).unwrap(), 0.0);
        let measure = |witness, receiver, transport, history| {
            source.reflected_slab_history_irradiance(
                cx, &scene, 0, witness, receiver, transport, history,
            )
        };
        for invalid in [
            Point3::new(0.0, 0.0, 0.1),
            Point3::new(0.0, 0.0, 0.2),
            Point3::new(7.0, 0.0, 0.0),
            Point3::new(f64::NAN, 0.0, 0.0),
        ] {
            assert!(
                measure(
                    invalid,
                    receiver(0.0),
                    TransportOptions::default(),
                    precise()
                )
                .is_err()
            );
        }
        assert!(measure(face(), probe(0.0), TransportOptions::default(), precise()).is_err());
        let inside =
            IlluminanceProbe::try_new(Point3::new(0.0, 0.0, 0.1), Vec3::new(0.0, 0.0, -1.0))
                .unwrap();
        assert!(measure(face(), inside, TransportOptions::default(), precise()).is_err());
        assert!(matches!(
            measure(
                face(),
                receiver(0.0),
                TransportOptions::default(),
                SlabHistoryOptions {
                    max_images_per_wavelength: 1,
                    ..precise()
                }
            ),
            Err(SlabHistoryError::ImageLimitExceeded {
                bin: 0,
                images: 1,
                ..
            })
        ));
        assert!(matches!(
            measure(
                face(),
                receiver(0.0),
                TransportOptions::default(),
                SlabHistoryOptions {
                    max_interface_interactions: 2,
                    ..precise()
                }
            ),
            Err(SlabHistoryError::InterfaceBudgetExceeded {
                bin: 0,
                round_trips: 1,
                required: 3,
                remaining: 1
            })
        ));
        assert!(matches!(
            measure(
                face(),
                receiver(0.0),
                TransportOptions {
                    max_distance_m: 1.8,
                    ..TransportOptions::default()
                },
                precise()
            ),
            Err(SlabHistoryError::Transport(
                PhotometryError::DistanceLimitExceeded
            ))
        ));
        assert!(matches!(
            measure(
                face(),
                receiver(0.0),
                TransportOptions::default(),
                SlabHistoryOptions {
                    relative_tolerance: f64::NAN,
                    ..precise()
                }
            ),
            Err(SlabHistoryError::InvalidOptions)
        ));
        gate.request();
        assert!(matches!(
            measure(
                face(),
                receiver(0.0),
                TransportOptions::default(),
                precise()
            ),
            Err(SlabHistoryError::Transport(PhotometryError::Cancelled))
        ));
    });
}
