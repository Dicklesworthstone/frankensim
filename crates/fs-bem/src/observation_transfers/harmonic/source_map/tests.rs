use super::*;

const LIMITS: HarmonicSourceMapLimits = HarmonicSourceMapLimits {
    maximum_sources: 8,
    maximum_products: 10_000,
    maximum_storage_bytes: 100_000,
};
const CONTROL: HarmonicControlLimits = HarmonicControlLimits {
    maximum_quadrature_points: 8,
    maximum_products: 10_000,
    maximum_workspace_bytes: 100_000,
};

fn fixture() -> PreparedHarmonicObservation {
    let mut samples = Vec::new();
    for source in 0..3 {
        for receiver in 0..2 {
            for component in 0..4 {
                let x = (1 + 3 * source + 2 * receiver + component) as f64;
                samples.push(C64::new(0.2 * x, (x * x - 3.0) * 0.07));
            }
        }
    }
    PreparedHarmonicObservation::from_samples(
        7.0,
        Medium {
            density: 1.2,
            sound_speed: 340.0,
        },
        &(0..3)
            .map(|s| SourceId::ApertureVolumeFlow(format!("opening-{s}")))
            .collect::<Vec<_>>(),
        &[
            ReceiverSpec {
                id: "first".into(),
                position_m: [1.0, 0.2, -0.3],
            },
            ReceiverSpec {
                id: "second".into(),
                position_m: [-0.4, 0.6, 1.5],
            },
        ],
        ContentHash([7; 32]),
        samples,
    )
    .unwrap()
}
fn sources() -> [SourceId; 2] {
    [
        SourceId::OpenCircuitPressure("pressure".into()),
        SourceId::PrescribedVolumeFlow("flow".into()),
    ]
}
fn matrix() -> [C64; 6] {
    [
        C64::new(0.3, -0.6),
        C64::new(0.1, 0.4),
        C64::new(-0.2, 0.7),
        C64::new(-0.4, 0.2),
        C64::new(0.8, -0.1),
        C64::new(0.5, 0.3),
    ]
}
fn close(a: C64, b: C64) {
    assert!(
        (a - b).abs() <= 2e-12 * b.abs().max(1e-10),
        "{a:?} != {b:?}"
    );
}

#[test]
fn original_fields_and_hermitian_pullback_compose_with_complete_complex_source_law() {
    let field = fixture();
    let map = matrix();
    let mapped = field
        .compose_sources(&sources(), &map, LIMITS, || false)
        .unwrap();
    assert_eq!(mapped.sources(), sources());
    assert_eq!(mapped.receivers(), field.receivers());
    assert_eq!(
        mapped.angular_frequency_rad_s(),
        field.angular_frequency_rad_s()
    );
    assert_eq!(mapped.medium.density, field.medium.density);
    assert_eq!(mapped.medium.sound_speed, field.medium.sound_speed);
    assert_ne!(mapped.producer_hash(), field.producer_hash());
    let a = [C64::new(0.8, -0.3), C64::new(-0.4, 0.2)];
    let old: Vec<_> = (0..3).map(|s| map[s] * a[0] + map[3 + s] * a[1]).collect();
    let (mut expected, mut actual) = ([HarmonicField::default(); 2], [HarmonicField::default(); 2]);
    field
        .evaluate_coherent(&old, &mut expected, || false)
        .unwrap();
    mapped.evaluate_coherent(&a, &mut actual, || false).unwrap();
    for (a, b) in actual.iter().zip(expected) {
        close(a.pressure_pa, b.pressure_pa);
        for i in 0..3 {
            close(a.particle_velocity_m_s[i], b.particle_velocity_m_s[i]);
        }
    }
    let covector = [HarmonicFieldCovector {
        pressure: C64::new(0.3, 0.6),
        particle_velocity: [
            C64::new(-0.2, 0.4),
            C64::new(0.5, -0.3),
            C64::new(-0.1, -0.2),
        ],
    }; 2];
    let (mut old_gradient, mut gradient) = ([C64::ZERO; 3], [C64::ZERO; 2]);
    let budget = HarmonicWorkBudget {
        maximum_products: 10_000,
    };
    field
        .source_pullback(&covector, budget, &mut old_gradient, || false)
        .unwrap();
    mapped
        .source_pullback(&covector, budget, &mut gradient, || false)
        .unwrap();
    for s in 0..2 {
        let chain = (0..3).fold(C64::ZERO, |sum, i| {
            sum + map[3 * s + i].conj() * old_gradient[i]
        });
        close(gradient[s], chain);
    }
}

#[test]
fn product_cancellation_is_retained_before_rounding_composed_columns() {
    // Ordinary f64 multiplication rounds (1+d)(1-d) to one. The physical
    // residual is -d²; forming a rounded intermediate column would erase it.
    let d = 2.0_f64.powi(-27);
    let mut field = fixture();
    field.sources.truncate(2);
    field.receivers.truncate(1);
    field.transfers = vec![
        [C64::from_re(1.0 + d), C64::ZERO, C64::ZERO, C64::ZERO],
        [C64::ONE, C64::ZERO, C64::ZERO, C64::ZERO],
    ];
    let composed = field
        .compose_sources(
            &[sources()[0].clone()],
            &[C64::from_re(1.0 - d), -C64::ONE],
            LIMITS,
            || false,
        )
        .unwrap();
    let mut output = [HarmonicField::default()];
    composed
        .evaluate_coherent(&[C64::ONE], &mut output, || false)
        .unwrap();
    assert_eq!(output[0].pressure_pa, C64::from_re(-d * d));
}

#[test]
fn transformation_identity_binds_map_labels_units_and_original_physics() {
    let field = fixture();
    let sources = sources();
    let map = matrix();
    let hash = field
        .compose_sources(&sources, &map, LIMITS, || false)
        .unwrap()
        .producer_hash();
    assert_eq!(
        hash,
        field
            .compose_sources(&sources, &map, LIMITS, || false)
            .unwrap()
            .producer_hash()
    );
    for change in 0..5 {
        let mut changed = fixture();
        let mut ids = sources.clone();
        let mut matrix = map;
        match change {
            0 => matrix[0].im += 0.125,
            1 => ids[0] = SourceId::PrescribedVolumeFlow("pressure".into()),
            2 => ids[0] = SourceId::OpenCircuitPressure("other-pressure".into()),
            3 => changed.producer_hash = ContentHash([8; 32]),
            _ => changed.receivers[0].position_m[1] += 0.1,
        }
        assert_ne!(
            hash,
            changed
                .compose_sources(&ids, &matrix, LIMITS, || false)
                .unwrap()
                .producer_hash()
        );
    }
}

#[test]
fn silent_maps_still_admit_identity_extent_work_finiteness_and_cancellation() {
    let field = fixture();
    let sources = sources();
    let zero = [C64::ZERO; 6];
    assert!(
        field
            .compose_sources(&sources, &zero[..5], LIMITS, || false)
            .is_err()
    );
    let duplicate = [
        sources[0].clone(),
        SourceId::PrescribedVolumeFlow("pressure".into()),
    ];
    assert!(
        field
            .compose_sources(&duplicate, &zero, LIMITS, || false)
            .is_err()
    );
    assert!(
        field
            .compose_sources(
                &[
                    SourceId::OpenCircuitPressure(String::new()),
                    sources[1].clone()
                ],
                &zero,
                LIMITS,
                || false
            )
            .is_err()
    );
    let mut invalid = zero;
    invalid[5].im = f64::NAN;
    assert!(
        field
            .compose_sources(&sources, &invalid, LIMITS, || false)
            .is_err()
    );
    for limits in [
        HarmonicSourceMapLimits {
            maximum_sources: 1,
            ..LIMITS
        },
        HarmonicSourceMapLimits {
            maximum_products: 47,
            ..LIMITS
        },
        HarmonicSourceMapLimits {
            maximum_storage_bytes: 1,
            ..LIMITS
        },
    ] {
        assert!(
            field
                .compose_sources(&sources, &zero, limits, || false)
                .is_err()
        );
    }
    let mut polls = 0;
    assert!(matches!(
        field.compose_sources(&sources, &zero, LIMITS, || {
            polls += 1;
            polls >= 15
        }),
        Err(TransferError::Cancelled)
    ));
    let silent = field
        .compose_sources(&sources, &zero, LIMITS, || false)
        .unwrap();
    let mut output = [HarmonicField::default(); 2];
    silent
        .evaluate_coherent(&[C64::ONE; 2], &mut output, || false)
        .unwrap();
    assert_eq!(output, [HarmonicField::default(); 2]);
}

#[test]
fn composition_cannot_invent_nondispersive_storage_for_a_relaxing_exterior() {
    let field = fixture().without_radiation_momentum();
    let composed = field
        .compose_sources(&sources(), &matrix(), LIMITS, || false)
        .unwrap();
    let quadrature = [HarmonicVolumeQuadrature {
        receiver: 0,
        volume_m3: 0.2,
    }];
    assert!(
        PreparedHarmonicEnergyControl::new(&composed, &quadrature, None, CONTROL, || false)
            .is_err()
    );
    // A Hermitian map pullback remains valid without asserting a storage law.
    let mut output = [C64::ONE; 2];
    composed
        .source_pullback(
            &[HarmonicFieldCovector::default(); 2],
            HarmonicWorkBudget {
                maximum_products: 100,
            },
            &mut output,
            || false,
        )
        .unwrap();
    assert_eq!(output, [C64::ZERO; 2]);
}
