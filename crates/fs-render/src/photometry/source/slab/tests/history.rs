use super::*;

fn precise() -> SlabHistoryOptions {
    SlabHistoryOptions {
        absolute_tolerance_w_m2_nm: 1e-16,
        relative_tolerance: 1e-11,
        ..SlabHistoryOptions::default()
    }
}

#[test]
fn g1_history_matches_independent_axial_series_and_bounds_its_omitted_power() {
    scope(|cx, _| {
        let source = source([0.0, 0.0, -1.0], AngularDistribution::isotropic());
        let receiver =
            IlluminanceProbe::try_new(probe(0.0).position_m(), Vec3::new(0.3, 0.4, -1.0)).unwrap();
        for index in [1.0_f64, 1.5, 2.4] {
            let scene = scene(CauchyIor::try_constant(index).unwrap(), 0.7);
            let result = source
                .transmitted_slab_history_irradiance(
                    cx,
                    &scene,
                    0,
                    receiver,
                    TransportOptions::default(),
                    precise(),
                )
                .unwrap();
            let r = ((index - 1.0) / (index + 1.0)).powi(2);
            let term = |m: i32| {
                let depth = f64::from(2 * m + 1) * 0.2;
                -receiver.normal().z * (1.0 - r).powi(2) * r.powi(2 * m) * (-0.7 * depth).exp()
                    / (2.0 + depth / index).powi(2)
            };
            let reference: f64 = (0..256).map(term).sum();
            let omitted = reference - result.reading().illuminance_lux;
            assert!(omitted >= -2e-13 * reference);
            assert!(omitted <= result.omitted_reading_bound().illuminance_lux + 2e-13 * reference);
            close(result.reading().illuminance_lux, reference, 2e-10);
            for bin in 0..SPECTRAL_BINS {
                assert!(
                    result.omitted_spectral_irradiance_bound_w_m2_nm()[bin]
                        <= precise().absolute_tolerance_w_m2_nm
                            + precise().relative_tolerance
                                * result.spectral_irradiance_w_m2_nm()[bin]
                );
            }
            let reserved: u64 = result
                .images_per_wavelength()
                .iter()
                .map(|n| u64::from(*n) * (u64::from(*n) + 1))
                .sum();
            assert_eq!(result.reserved_interface_interactions(), reserved);
            if index == 1.0 {
                assert!(result.images_per_wavelength().iter().all(|n| *n == 1));
                assert_eq!(result.omitted_reading_bound().illuminance_lux, 0.0);
            } else {
                assert!(result.reading().illuminance_lux > term(0));
            }
        }
    });
}

// Independent folded forward ray: each return actually reverses z, then the
// final outgoing ray intersects the supplied tilted receiver plane.
fn folded_endpoint(direction: Vec3, index: f64, receiver: IlluminanceProbe, order: u32) -> Point3 {
    let mut inside = Vec3::new(
        direction.x / index,
        direction.y / index,
        (1.0 - (direction.x * direction.x + direction.y * direction.y) / (index * index)).sqrt(),
    );
    let mut point = Point3::new(0.0, 0.0, -1.0).offset(direction.scale(1.0 / direction.z));
    for crossing in 0..=2 * order {
        point = point.offset(inside.scale(0.2 / inside.z.abs()));
        if crossing != 2 * order {
            inside.z = -inside.z;
        }
    }
    let parameter = receiver
        .position_m()
        .delta_from(point)
        .dot(receiver.normal())
        / direction.dot(receiver.normal());
    point.offset(direction.scale(parameter))
}

fn folded_area(direction: Vec3, index: f64, receiver: IlluminanceProbe, order: u32) -> f64 {
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
        folded_endpoint(plus.scale(1.0 / plus.norm()), index, receiver, order)
            .delta_from(folded_endpoint(
                minus.scale(1.0 / minus.norm()),
                index,
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

#[test]
fn g1_return_images_have_independent_snell_focusing_polarization_and_beer_history() {
    scope(|cx, _| {
        let law = CauchyIor::try_new(1.42, 0.018, 0.0).unwrap();
        let scene = scene(law, 0.3);
        let receiver =
            IlluminanceProbe::try_new(probe(3.6).position_m(), Vec3::new(0.1, 0.2, -1.0)).unwrap();
        let connection = SlabConnection::try_new(
            cx,
            &scene,
            0,
            Point3::new(0.0, 0.0, -1.0),
            receiver,
            TransportOptions::default(),
        )
        .unwrap();
        for bin in [0, 40, 79] {
            let wavelength = wavelength_nm(bin);
            let n = law.eval(wavelength).unwrap();
            let mut previous_cosine = 0.0;
            for order in [0, 1, 3] {
                let image = connection.image(cx, wavelength, order).unwrap();
                let ray = connection.sample_image(cx, &image).unwrap();
                let d = ray.launch_direction;
                assert!(d.z > previous_cosine);
                previous_cosine = d.z;
                assert!(
                    folded_endpoint(d, n, receiver, order)
                        .delta_from(receiver.position_m())
                        .norm()
                        < 2e-11
                );
                close(
                    ray.area_per_solid_angle,
                    folded_area(d, n, receiver, order),
                    8e-8,
                );
                let cg = (1.0 - d.x * d.x / (n * n)).sqrt();
                let rs = ((d.z - n * cg) / (d.z + n * cg)).powi(2);
                let rp = ((n * d.z - cg) / (n * d.z + cg)).powi(2);
                let ts = 4.0 * n * d.z * cg / (d.z + n * cg).powi(2);
                let tp = 4.0 * n * d.z * cg / (n * d.z + cg).powi(2);
                let power = 0.5
                    * (ts * ts * rs.powi(2 * order as i32) + tp * tp * rp.powi(2 * order as i32))
                    * (-0.3 * f64::from(2 * order + 1) * 0.2 / cg).exp();
                close(ray.power_transmission, power, 2e-12);
                assert!(ray.unoccluded);
            }
        }
    });
}

#[test]
fn g3_dark_or_occluded_early_image_does_not_discard_later_transmitted_light() {
    scope(|cx, _| {
        let mut scene = scene(CauchyIor::try_constant(1.5).unwrap(), 0.1);
        let receiver = probe(3.6);
        let connection = SlabConnection::try_new(
            cx,
            &scene,
            0,
            Point3::new(0.0, 0.0, -1.0),
            receiver,
            TransportOptions::default(),
        )
        .unwrap();
        let first = connection.sample(cx, 550.0).unwrap();
        let second = connection
            .sample_image(cx, &connection.image(cx, 550.0, 1).unwrap())
            .unwrap();
        let cutoff = (0.5 * (first.launch_direction.z + second.launch_direction.z)).acos();
        let source = source(
            [0.0, 0.0, -1.0],
            AngularDistribution::uniform_cone(cutoff).unwrap(),
        );
        let direct = source
            .transmitted_slab_irradiance(cx, &scene, 0, receiver, TransportOptions::default())
            .unwrap();
        assert_eq!(direct.reading().illuminance_lux, 0.0);
        let clear = source
            .transmitted_slab_history_irradiance(
                cx,
                &scene,
                0,
                receiver,
                TransportOptions::default(),
                precise(),
            )
            .unwrap();
        assert!(clear.reading().illuminance_lux > 0.0);
        let second_lux = source.intensity_cd(second.launch_direction).unwrap()
            / second.area_per_solid_angle
            * second.power_transmission;
        let x = 3.6 - 0.5 * second.launch_direction.x / second.launch_direction.z;
        scene.primitives.push(panel(x, 0.7, 0.001));
        let blocked = source
            .transmitted_slab_history_irradiance(
                cx,
                &scene,
                0,
                receiver,
                TransportOptions::default(),
                precise(),
            )
            .unwrap();
        assert!(blocked.reading().illuminance_lux > 0.0);
        let lost = clear.reading().illuminance_lux - blocked.reading().illuminance_lux;
        assert!(
            (lost - second_lux).abs()
                <= clear.omitted_reading_bound().illuminance_lux
                    + blocked.omitted_reading_bound().illuminance_lux
                    + 2e-12 * second_lux
        );
    });
}

#[test]
fn g1_receiver_can_face_only_later_images_or_be_grazing_to_all_images() {
    scope(|cx, _| {
        let source = source([0.0, 0.0, -1.0], AngularDistribution::isotropic());
        let scene = scene(CauchyIor::try_constant(1.5).unwrap(), 0.0);
        let connection = SlabConnection::try_new(
            cx,
            &scene,
            0,
            source.position_m(),
            probe(3.6),
            TransportOptions::default(),
        )
        .unwrap();
        let d0 = connection.sample(cx, 550.0).unwrap().launch_direction;
        let d1 = connection
            .sample_image(cx, &connection.image(cx, 550.0, 1).unwrap())
            .unwrap()
            .launch_direction;
        let slope = 0.5 * (d0.x / d0.z + d1.x / d1.z);
        let receiver =
            IlluminanceProbe::try_new(probe(3.6).position_m(), Vec3::new(1.0, 0.0, -slope))
                .unwrap();
        assert!(receiver.normal().dot(d0) > 0.0);
        assert!(receiver.normal().dot(d1) < 0.0);
        assert!(
            source
                .transmitted_slab_irradiance(cx, &scene, 0, receiver, TransportOptions::default())
                .is_err()
        );
        let result = source
            .transmitted_slab_history_irradiance(
                cx,
                &scene,
                0,
                receiver,
                TransportOptions::default(),
                precise(),
            )
            .unwrap();
        assert!(result.reading().illuminance_lux > 0.0);
        let grazing =
            IlluminanceProbe::try_new(receiver.position_m(), Vec3::new(0.0, 1.0, 0.0)).unwrap();
        let dark = source
            .transmitted_slab_history_irradiance(
                cx,
                &scene,
                0,
                grazing,
                TransportOptions::default(),
                precise(),
            )
            .unwrap();
        assert_eq!(dark.reading().illuminance_lux, 0.0);
        assert!(dark.images_per_wavelength().iter().all(|n| *n > 1));
    });
}

#[test]
fn g1_tail_keeps_rescaled_power_when_unscaled_reflection_or_beer_underflows() {
    scope(|cx, _| {
        for (extinction, order) in [(0.0, 128), (5000.0, 0)] {
            let scene = scene(CauchyIor::try_constant(1.5).unwrap(), extinction);
            let connection = SlabConnection::try_new(
                cx,
                &scene,
                0,
                Point3::new(0.0, 0.0, -1.0),
                probe(0.0),
                TransportOptions::default(),
            )
            .unwrap();
            let image = connection.image(cx, 550.0, order).unwrap();
            let bound = connection.omitted_bound(cx, &image, 1e295).unwrap();
            let depth = f64::from(2 * order + 1) * 0.2;
            let width = 2.0 + depth / 1.5;
            let log_term = 295.0 * 10.0_f64.ln() - extinction * depth
                + 2.0 * f64::from(order) * 0.04_f64.ln()
                + 2.0 * 0.96_f64.ln()
                - 2.0 * width.ln();
            let first_omitted = log_term.exp();
            assert!(first_omitted > 1e-200);
            assert!(bound >= first_omitted, "{bound:e} < {first_omitted:e}");
            assert!(bound.is_finite());
        }
    });
}

#[test]
fn g1_passive_coated_history_preserves_reversed_film_power_and_bounds_its_tail() {
    use crate::conductor::{
        CauchyCoatingLayer, CoherentPeriodicCoating, ConductorDataStatus, ConductorSource,
    };
    use crate::dielectric::DielectricCoating;
    use crate::polarization::Stokes;
    scope(|cx, _| {
        let mut scene = scene(CauchyIor::try_constant(1.5).unwrap(), 12.0);
        let provenance = ConductorSource::try_new(
            fs_blake3::hash_domain("fs-render:slab-history-test", b"absorbing-surface-film"),
            ConductorDataStatus::Representative,
        )
        .unwrap();
        let layer = CauchyCoatingLayer::try_new(CauchyIor::try_constant(2.0).unwrap(), 0.15, 80e-9)
            .unwrap();
        let surface = DielectricSurface::SMOOTH.with_coating(
            DielectricCoating::try_new(
                CoherentPeriodicCoating::try_new(&[layer], 1, provenance).unwrap(),
            )
            .unwrap(),
        );
        if let Material::Dielectric {
            surface: authored, ..
        } = &mut scene.primitives[0].material
        {
            *authored = surface;
        }
        let source = source([0.0, 0.0, -1.0], AngularDistribution::isotropic());
        let result = source
            .transmitted_slab_history_irradiance(
                cx,
                &scene,
                0,
                probe(0.0),
                TransportOptions::default(),
                precise(),
            )
            .unwrap();
        for bin in [0, 40, 79] {
            let wavelength = wavelength_nm(bin);
            let entry = surface
                .interface(wavelength, false)
                .unwrap()
                .fresnel_polarized(1.0, 1.0, 1.5)
                .unwrap();
            let inside = surface
                .interface(wavelength, true)
                .unwrap()
                .fresnel_polarized(1.0, 1.5, 1.0)
                .unwrap();
            let mut reference = 0.0;
            for q in [-1.0, 1.0] {
                let polarization = Stokes {
                    i: 1.0,
                    q,
                    u: 0.0,
                    v: 0.0,
                };
                let transmission = entry.transmission.apply(polarization).i
                    * inside.transmission.apply(polarization).i;
                let reflection = inside.reflection.apply(polarization).i;
                for order in 0..64 {
                    let depth = f64::from(2 * order + 1) * 0.2;
                    reference +=
                        0.5 * transmission * reflection.powi(2 * order) * (-12.0 * depth).exp()
                            / (2.0 + depth / 1.5).powi(2);
                }
            }
            reference *= source
                .spectral_intensity_w_sr_nm(Vec3::new(0.0, 0.0, 1.0), wavelength)
                .unwrap();
            let omitted = reference - result.spectral_irradiance_w_m2_nm()[bin];
            assert!(omitted >= -2e-12 * reference);
            assert!(
                omitted
                    <= result.omitted_spectral_irradiance_bound_w_m2_nm()[bin] + 2e-12 * reference
            );
            close(result.spectral_irradiance_w_m2_nm()[bin], reference, 2e-9);
        }
    });
}

#[test]
fn g1_bright_source_recovers_representable_irradiance_after_beer_underflow() {
    scope(|cx, _| {
        let source = PhotometricSource::try_new(
            RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [0.0, 0.0, -1.0]).unwrap(),
            4.0 * core::f64::consts::PI * 1e295,
            1.0,
            AngularDistribution::isotropic(),
            SourceSpectrum::equal_energy(),
        )
        .unwrap();
        let scene = scene(CauchyIor::try_constant(1.5).unwrap(), 5000.0);
        let direct = source
            .transmitted_slab_irradiance(cx, &scene, 0, probe(0.0), TransportOptions::default())
            .unwrap();
        let history = source
            .transmitted_slab_history_irradiance(
                cx,
                &scene,
                0,
                probe(0.0),
                TransportOptions::default(),
                SlabHistoryOptions {
                    absolute_tolerance_w_m2_nm: 1e-160,
                    relative_tolerance: 0.0,
                    ..precise()
                },
            )
            .unwrap();
        for bin in 0..SPECTRAL_BINS {
            let intensity = source
                .spectral_intensity_w_sr_nm(Vec3::new(0.0, 0.0, 1.0), wavelength_nm(bin))
                .unwrap();
            let expected =
                (intensity.ln() + 2.0 * 0.96_f64.ln() - 1000.0 - 2.0 * (2.0_f64 + 0.2 / 1.5).ln())
                    .exp();
            assert!(expected > 0.0);
            close(
                direct.spectral_irradiance_w_m2_nm()[bin] / expected,
                1.0,
                2e-12,
            );
            close(
                history.spectral_irradiance_w_m2_nm()[bin] / expected,
                1.0,
                2e-12,
            );
            assert_eq!(history.images_per_wavelength()[bin], 1);
        }
        assert!(direct.reading().illuminance_lux > 0.0);
        assert!(history.reading().illuminance_lux > 0.0);
    });
}

#[test]
fn g0_history_refuses_incomplete_image_or_global_work_budgets_distance_and_cancellation() {
    scope(|cx, gate| {
        let source = source([0.0, 0.0, -1.0], AngularDistribution::isotropic());
        let scene = scene(CauchyIor::try_constant(1.5).unwrap(), 0.0);
        assert!(matches!(
            source.transmitted_slab_history_irradiance(
                cx,
                &scene,
                0,
                probe(0.0),
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
            source.transmitted_slab_history_irradiance(
                cx,
                &scene,
                0,
                probe(0.0),
                TransportOptions::default(),
                SlabHistoryOptions {
                    max_interface_interactions: 3,
                    ..precise()
                }
            ),
            Err(SlabHistoryError::InterfaceBudgetExceeded {
                bin: 0,
                round_trips: 1,
                required: 4,
                remaining: 1
            })
        ));
        assert!(matches!(
            source.transmitted_slab_history_irradiance(
                cx,
                &scene,
                0,
                probe(0.0),
                TransportOptions {
                    max_distance_m: 2.4,
                    ..TransportOptions::default()
                },
                precise()
            ),
            Err(SlabHistoryError::Transport(
                PhotometryError::DistanceLimitExceeded
            ))
        ));
        let matched = super::scene(CauchyIor::try_constant(1.0).unwrap(), 0.0);
        assert!(matches!(
            source.transmitted_slab_history_irradiance(
                cx,
                &matched,
                0,
                probe(0.0),
                TransportOptions::default(),
                SlabHistoryOptions {
                    max_interface_interactions: 2,
                    ..precise()
                }
            ),
            Err(SlabHistoryError::InterfaceBudgetExceeded {
                bin: 1,
                round_trips: 0,
                ..
            })
        ));
        gate.request();
        assert!(matches!(
            source.transmitted_slab_history_irradiance(
                cx,
                &scene,
                0,
                probe(0.0),
                TransportOptions::default(),
                precise()
            ),
            Err(SlabHistoryError::Transport(PhotometryError::Cancelled))
        ));
    });
}
