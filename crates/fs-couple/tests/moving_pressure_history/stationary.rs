//! G1/G3/G4: actual accepted source amplitudes through vector fields and PCM.
use super::*;
use fs_aeroac::curle3d::{CompactForceObserver, IdealAcousticPlane, IdealPlaneBoundary};
use fs_aeroac::curle3d::moving::trajectory::IdealReflectionObserver;
use fs_bem::near_field::FirstOrder;
use fs_couple::pcm_wav::observation::moving::stationary::{
    StationaryCompactRenderer, StationaryCompactSample,
};

fn impulse(t: f64, scale: f64) -> [f64; 3] {
    [
        scale * (0.7 + t + t * t - t.powi(3) + t.powi(4)),
        scale * (-0.2 - 2.0 * t + 1.5 * t * t + t.powi(3) / 3.0 - 0.5 * t.powi(4)),
        scale * (0.4 + 0.5 * t - 0.5 * t * t + 2.0 * t.powi(3) / 3.0 + 0.75 * t.powi(4)),
    ]
}

fn build_fields<'a>(
    sources: &[CompactSourceHistory<'a>],
    initial: &[[f64; 3]],
    microphone: FirstOrder,
    config: MovingPressureConfig,
) -> StationaryCompactRenderer<'a> {
    StationaryCompactRenderer::new(
        observer([12.0, 0.0, 0.0]),
        sources,
        initial,
        microphone,
        config,
        budget(),
        &mut || false,
    )
    .unwrap()
}

#[test]
fn original_cubic_history_drives_force_impulse_monopole_velocity_and_microphone() {
    let accepted = history(&[-0.25, -0.1, 0.05, 0.2, 0.5], 0.0, 1.0);
    let sources = [source(&accepted, [0.0; 3], [0.0; 3])];
    let mic = FirstOrder::new(0.5, [-1.0, 0.0, 0.0]).unwrap();
    let mut renderer = build_fields(&sources, &[impulse(-0.25, 1.0)], mic, config());
    let mut fields = vec![StationaryCompactSample::default(); 81];
    renderer.field_block(&mut fields).unwrap();
    let mut legacy = build(&sources, config());
    let mut pressure_only = [0.0; 81];
    legacy.block(&mut pressure_only).unwrap();
    let force_observer = CompactForceObserver::new(
        observer([12.0, 0.0, 0.0]).medium(),
        [12.0, 0.0, 0.0],
        [0.0; 3],
    )
    .unwrap();
    for (i, field) in fields.iter().enumerate() {
        let time = renderer.receiver_time_at(i as u64).unwrap();
        let tau = time - 12.0 / 300.0;
        let source = strength(tau, 1.0);
        let force = force_observer
            .field_at_retarded(
                impulse(tau, 1.0),
                source.force_on_fluid_n,
                source.force_rate_on_fluid_n_s,
            )
            .unwrap();
        let pressure = force.pressure_pa + 1.2 * source.volume_flow_rate_m3_s2 / (4.0 * PI * 12.0);
        let mut velocity = force.particle_velocity_m_s;
        velocity[0] += (source.volume_flow_m3_s / (12.0 * 12.0)
            + source.volume_flow_rate_m3_s2 / (300.0 * 12.0))
            / (4.0 * PI);
        assert_eq!(field.receiver_time_s, time);
        close(field.pressure_pa, pressure);
        close(field.pressure_pa, pressure_only[i]);
        for axis in 0..3 {
            close(field.particle_velocity_m_s[axis], velocity[axis]);
            close(
                field.instantaneous_intensity_w_m2[axis],
                pressure * velocity[axis],
            );
        }
        close(
            field.microphone_pressure_pa,
            0.5 * (pressure + 1.2 * 300.0 * velocity[0]),
        );
    }
}

#[test]
fn explicit_prehistory_changes_particle_motion_without_inventing_pressure() {
    let accepted: Vec<_> = [0.0, 0.25, 0.5]
        .into_iter()
        .map(|source_time_s| TimedCompactSource {
            source_time_s,
            strength: CompactSourceStrength::default(),
        })
        .collect();
    let sources = [source(&accepted, [0.0; 3], [0.0; 3])];
    let dipole = FirstOrder::new(0.0, [0.0, 1.0, 0.0]).unwrap();
    let mut renderer = build_fields(&sources, &[[0.0, 2.0, 0.0]], dipole, config());
    let mut fields = [StationaryCompactSample::default(); 81];
    renderer.field_block(&mut fields).unwrap();
    let transverse = -2.0 / (4.0 * PI * 1.2 * 12.0_f64.powi(3));
    for field in fields {
        assert_eq!(field.pressure_pa, 0.0);
        close(field.particle_velocity_m_s[1], transverse);
        assert_eq!(field.instantaneous_intensity_w_m2, [0.0; 3]);
        close(field.microphone_pressure_pa, -1.2 * 300.0 * transverse);
    }
    let mut zero_history = build_fields(&sources, &[[0.0; 3]], dipole, config());
    let mut silence = [1.0; 81];
    zero_history.block(&mut silence).unwrap();
    assert_eq!(silence, [0.0; 81]);
}

#[test]
fn coherent_amplitudes_are_summed_before_quadratic_flux_and_never_normalized() {
    let positive = history(&[-0.25, 0.0, 0.25, 0.5], 0.0, 1.0);
    let negative = history(&[-0.25, 0.0, 0.25, 0.5], 0.0, -1.0);
    let first = source(&positive, [0.0; 3], [0.0; 3]);
    let second = source(&negative, [0.0; 3], [0.0; 3]);
    let mut one = build_fields(
        &[first],
        &[impulse(-0.25, 1.0)],
        FirstOrder::default(),
        config(),
    );
    let mut two = build_fields(
        &[first, first],
        &[impulse(-0.25, 1.0); 2],
        FirstOrder::default(),
        config(),
    );
    let mut cancelling = build_fields(
        &[first, second],
        &[impulse(-0.25, 1.0), impulse(-0.25, -1.0)],
        FirstOrder::default(),
        config(),
    );
    let mut baseline = [StationaryCompactSample::default(); 81];
    let mut doubled = baseline;
    let mut null = baseline;
    one.field_block(&mut baseline).unwrap();
    two.field_block(&mut doubled).unwrap();
    cancelling.field_block(&mut null).unwrap();
    for ((base, double), zero) in baseline.iter().zip(doubled).zip(null) {
        close(double.pressure_pa, 2.0 * base.pressure_pa);
        close(
            double.microphone_pressure_pa,
            2.0 * base.microphone_pressure_pa,
        );
        for axis in 0..3 {
            close(
                double.particle_velocity_m_s[axis],
                2.0 * base.particle_velocity_m_s[axis],
            );
            close(
                double.instantaneous_intensity_w_m2[axis],
                4.0 * base.instantaneous_intensity_w_m2[axis],
            );
        }
        assert_eq!(zero.pressure_pa, 0.0);
        assert_eq!(zero.particle_velocity_m_s, [0.0; 3]);
        assert_eq!(zero.microphone_pressure_pa, 0.0);
        assert_eq!(zero.instantaneous_intensity_w_m2, [0.0; 3]);
    }
}

#[test]
fn field_and_pcm_callbacks_share_history_and_are_partition_invariant() {
    let accepted = history(&[-0.25, -0.1, 0.05, 0.2, 0.5], 0.0, 1.0);
    let sources = [source(&accepted, [0.0; 3], [0.0; 3])];
    let mic = FirstOrder::new(0.5, [-1.0, 0.0, 0.0]).unwrap();
    let make = || build_fields(&sources, &[impulse(-0.25, 1.0)], mic, config());
    let mut whole = make();
    let mut expected = [StationaryCompactSample::default(); 81];
    whole.field_block(&mut expected).unwrap();
    let mut split = make();
    let mut first = [StationaryCompactSample::default(); 7];
    split.field_block(&mut first).unwrap();
    assert_eq!(&first, &expected[..7]);
    let mut rest = [0.0; 74];
    for chunk in rest.chunks_mut(11) {
        PressureRenderer::block(&mut split, chunk).unwrap();
    }
    let expected_pressure: Vec<_> = expected.iter().map(|f| f.microphone_pressure_pa).collect();
    assert_eq!(&rest, &expected_pressure[7..]);
    assert_eq!(split.samples_rendered(), 81);
    assert!(split.block(&mut [0.0]).is_err());
    let mut direct = make();
    let mut pressure = [0.0; 81];
    direct.block(&mut pressure).unwrap();
    assert_eq!(
        encode_pcm16_wav(&pressure, 1000, 1.0).unwrap(),
        encode_pcm16_wav(&expected_pressure, 1000, 1.0).unwrap()
    );
}

#[test]
fn cancellation_is_atomic_and_additional_history_storage_is_admitted() {
    let accepted = history(&[-0.25, 0.0, 0.25, 0.5], 0.0, 1.0);
    let sources = [source(&accepted, [0.0; 3], [0.0; 3])];
    let initial = [impulse(-0.25, 1.0)];
    let mut renderer = build_fields(&sources, &initial, FirstOrder::default(), config());
    let sentinel = StationaryCompactSample {
        pressure_pa: 123.0,
        ..Default::default()
    };
    let mut output = [sentinel; 81];
    let mut polls = 0;
    assert!(
        renderer
            .field_block_controlled(&mut output, &mut || {
                polls += 1;
                polls == 8
            })
            .is_err()
    );
    assert_eq!(output, [sentinel; 81]);
    assert_eq!(renderer.samples_rendered(), 0);
    renderer.field_block(&mut output).unwrap();
    let mut fresh = build_fields(&sources, &initial, FirstOrder::default(), config());
    let mut expected = [StationaryCompactSample::default(); 81];
    fresh.field_block(&mut expected).unwrap();
    assert_eq!(output, expected);
    let admit = |sources: &[CompactSourceHistory<'_>], impulse: &[[f64; 3]], budget| {
        StationaryCompactRenderer::new(
            observer([12.0, 0.0, 0.0]),
            sources,
            impulse,
            FirstOrder::default(),
            config(),
            budget,
            &mut || false,
        )
        .is_ok()
    };
    assert!(!admit(&sources, &[], budget()));
    assert!(!admit(&sources, &[[f64::NAN; 3]], budget()));
    let moving = [source(&accepted, [0.0; 3], [1.0, 0.0, 0.0])];
    assert!(!admit(&moving, &initial, budget()));
    let mut too_small = budget();
    too_small.maximum_payload_bytes = 1;
    assert!(!admit(&sources, &initial, too_small));
}

fn reflected_strength(time: f64, scale: f64) -> CompactSourceStrength {
    let mut value = strength(time, scale);
    value.volume_flow_m3_s = scale * (0.3 + 0.7 * time - 0.2 * time * time + 0.1 * time.powi(3));
    value.volume_flow_rate_m3_s2 = scale * (0.7 - 0.4 * time + 0.3 * time * time);
    value
}

fn reflected_history(times: &[f64], scale: f64) -> Vec<TimedCompactSource> {
    times
        .iter()
        .map(|&source_time_s| TimedCompactSource {
            source_time_s,
            strength: reflected_strength(source_time_s, scale),
        })
        .collect()
}

fn reflected_observer(position_m: [f64; 3]) -> MovingCompactObserver {
    MovingCompactObserver::new(
        FluidAcousticState::try_new(2.0, 5.0).unwrap(),
        position_m,
        MovingCompactLimits {
            maximum_mach: 0.9,
            minimum_distance_m: 0.01,
        },
    )
    .unwrap()
}

fn reflected_config() -> MovingPressureConfig {
    MovingPressureConfig {
        first_receiver_time_s: 1.4,
        samples: 9,
        max_block: 9,
        ..config()
    }
}

// Independent SOURCE images, with explicit F/Fdot/impulse reflection. The
// implementation instead observes unmodified histories at virtual receivers.
// The full Green derivatives retain the 1/r, 1/r² and 1/r³ velocity terms.
fn reflected_reference(
    receiver: [f64; 3],
    source: [f64; 3],
    time: f64,
    scale: f64,
    planes: &[IdealAcousticPlane],
) -> (f64, [f64; 3]) {
    let dot = |a: [f64; 3], b: [f64; 3]| (0..3).map(|i| a[i] * b[i]).sum::<f64>();
    let mut pressure = 0.0;
    let mut velocity = [0.0; 3];
    for mask in 0..1_usize << planes.len() {
        let mut image = source;
        let mut sign = 1.0;
        for (index, plane) in planes.iter().enumerate() {
            if mask & (1 << index) != 0 {
                let normal = plane.inward_normal();
                let distance = dot(
                    core::array::from_fn(|i| image[i] - plane.point_m()[i]),
                    normal,
                );
                image = core::array::from_fn(|i| image[i] - 2.0 * distance * normal[i]);
                if plane.boundary() == IdealPlaneBoundary::PressureRelease {
                    sign = -sign;
                }
            }
        }
        let delta: [f64; 3] = core::array::from_fn(|i| receiver[i] - image[i]);
        let radius = dot(delta, delta).sqrt();
        let direction = delta.map(|x| x / radius);
        let tau = time - radius / 5.0;
        let source = reflected_strength(tau, scale);
        let reflect = |mut vector: [f64; 3]| {
            for (index, plane) in planes.iter().enumerate() {
                if mask & (1 << index) != 0 {
                    let normal = plane.inward_normal();
                    let projection = dot(vector, normal);
                    vector = core::array::from_fn(|i| vector[i] - 2.0 * projection * normal[i]);
                }
            }
            vector
        };
        let force = reflect(source.force_on_fluid_n);
        let rate = reflect(source.force_rate_on_fluid_n_s);
        let impulse = reflect(impulse(tau, scale));
        let radial_force = dot(direction, force);
        let radial_rate = dot(direction, rate);
        let radial_impulse = dot(direction, impulse);
        let r2 = radius * radius;
        pressure += sign
            * (2.0 * source.volume_flow_rate_m3_s2 / radius
                + radial_force / r2
                + radial_rate / (5.0 * radius))
            / (4.0 * PI);
        for axis in 0..3 {
            velocity[axis] += sign
                * (direction[axis]
                    * (source.volume_flow_m3_s / r2
                        + source.volume_flow_rate_m3_s2 / (5.0 * radius))
                    + (3.0 * direction[axis] * radial_impulse - impulse[axis])
                        / (2.0 * r2 * radius)
                    + (3.0 * direction[axis] * radial_force - force[axis]) / (10.0 * r2)
                    + direction[axis] * radial_rate / (50.0 * radius))
                / (4.0 * PI);
        }
    }
    (pressure, velocity)
}

#[test]
fn reflected_cubic_fields_match_independent_source_images_in_rotated_corners() {
    let first = reflected_history(&[-0.25, 0.05, 0.4, 0.8, 1.2, 1.6, 2.0], 1.0);
    let second = reflected_history(&[-0.2, 0.15, 0.5, 0.9, 1.3, 1.65, 2.0], -0.4);
    let positions = [[0.28, 2.34, 0.8], [0.7, 1.8, 0.6]];
    let sources = [
        source(&first, positions[0], [0.0; 3]),
        source(&second, positions[1], [0.0; 3]),
    ];
    let initial = [impulse(-0.25, 1.0), impulse(-0.2, -0.4)];
    let receiver = [0.39, 1.32, 1.3];
    let point = [0.2, -0.1, 0.4];
    let planes = [
        IdealAcousticPlane::new(point, [0.6, 0.8, 0.0], IdealPlaneBoundary::Rigid).unwrap(),
        IdealAcousticPlane::new(point, [-0.8, 0.6, 0.0], IdealPlaneBoundary::PressureRelease)
            .unwrap(),
        IdealAcousticPlane::new(point, [0.0, 0.0, 1.0], IdealPlaneBoundary::Rigid).unwrap(),
    ];
    let microphone = FirstOrder::new(0.25, [0.0, 0.0, 1.0]).unwrap();
    for count in 1..=3 {
        let images = IdealReflectionObserver::new(
            reflected_observer(receiver),
            &planes[..count],
            1 << count,
        )
        .unwrap();
        let mut limits = budget();
        limits.maximum_sources = 2;
        limits.maximum_source_samples = first.len() + second.len();
        limits.maximum_source_receiver_pairs = (2 * (1_u64 << count)) * 9;
        let mut renderer = StationaryCompactRenderer::with_ideal_reflections(
            &images,
            &sources,
            &initial,
            microphone,
            reflected_config(),
            limits,
            &mut || false,
        )
        .unwrap();
        assert_eq!(renderer.source_path_count(), 2 * (1 << count));
        let mut fields = [StationaryCompactSample::default(); 9];
        renderer.field_block(&mut fields).unwrap();
        for field in fields {
            let (p1, u1) = reflected_reference(
                receiver,
                positions[0],
                field.receiver_time_s,
                1.0,
                &planes[..count],
            );
            let (p2, u2) = reflected_reference(
                receiver,
                positions[1],
                field.receiver_time_s,
                -0.4,
                &planes[..count],
            );
            let pressure = p1 + p2;
            close(field.pressure_pa, pressure);
            for axis in 0..3 {
                let velocity = u1[axis] + u2[axis];
                close(field.particle_velocity_m_s[axis], velocity);
                close(field.instantaneous_intensity_w_m2[axis], pressure * velocity);
            }
            close(field.microphone_pressure_pa, 0.25 * pressure - 7.5 * (u1[2] + u2[2]));
        }
    }
}

#[test]
fn rotated_wall_boundary_conditions_retain_directional_signal_at_a_pressure_node() {
    let accepted = reflected_history(&[-0.25, 0.05, 0.4, 0.8, 1.2, 1.6, 2.0], 1.0);
    let position = [0.96, 1.78, 0.4];
    let sources = [source(&accepted, position, [0.0; 3])];
    let normal = [0.6, 0.8, 0.0];
    // Commuted products give exactly zero represented signed plane distance.
    let receiver = [-0.8, 0.6, 0.9];
    for boundary in [
        IdealPlaneBoundary::Rigid,
        IdealPlaneBoundary::PressureRelease,
    ] {
        let plane = IdealAcousticPlane::new([0.0; 3], normal, boundary).unwrap();
        let images =
            IdealReflectionObserver::new(reflected_observer(receiver), &[plane], 2).unwrap();
        let mut renderer = StationaryCompactRenderer::with_ideal_reflections(
            &images,
            &sources,
            &[impulse(-0.25, 1.0)],
            FirstOrder::new(0.0, normal).unwrap(),
            reflected_config(),
            budget(),
            &mut || false,
        )
        .unwrap();
        let mut fields = [StationaryCompactSample::default(); 9];
        renderer.field_block(&mut fields).unwrap();
        for field in fields {
            let normal_velocity = normal[0] * field.particle_velocity_m_s[0]
                + normal[1] * field.particle_velocity_m_s[1];
            if boundary == IdealPlaneBoundary::Rigid {
                let (direct, _) =
                    reflected_reference(receiver, position, field.receiver_time_s, 1.0, &[]);
                close(field.pressure_pa, 2.0 * direct);
                assert!(normal_velocity.abs() < 1e-13);
            } else {
                assert!(field.pressure_pa.abs() < 1e-13);
                let tangent_velocity = -0.8 * field.particle_velocity_m_s[0]
                    + 0.6 * field.particle_velocity_m_s[1];
                assert!(tangent_velocity.abs() < 1e-13);
                assert!(field.particle_velocity_m_s[2].abs() < 1e-13);
                assert!(normal_velocity.abs() > 1e-5);
                close(field.microphone_pressure_pa, -10.0 * normal_velocity);
                assert!(
                    field
                        .instantaneous_intensity_w_m2
                        .iter()
                        .all(|x| x.abs() < 1e-13)
                );
            }
        }
    }
}

#[test]
fn reflected_history_work_and_callbacks_use_all_paths_without_duplicating_sources() {
    let accepted = reflected_history(&[-0.25, 0.05, 0.4, 0.8, 1.2, 1.6, 2.0], 1.0);
    let sources = [source(&accepted, [1.0, 0.0, 0.0], [0.0; 3])];
    let initial = [impulse(-0.25, 1.0)];
    let plane =
        IdealAcousticPlane::new([0.0; 3], [1.0, 0.0, 0.0], IdealPlaneBoundary::Rigid).unwrap();
    let images =
        IdealReflectionObserver::new(reflected_observer([3.0, 0.0, 0.0]), &[plane], 2).unwrap();
    let microphone = FirstOrder::new(0.5, [-1.0, 0.0, 0.0]).unwrap();
    let mut limits = budget();
    limits.maximum_sources = 1;
    limits.maximum_source_samples = accepted.len();
    limits.maximum_source_receiver_pairs = 18;
    let make = |limits| {
        StationaryCompactRenderer::with_ideal_reflections(
            &images,
            &sources,
            &initial,
            microphone,
            reflected_config(),
            limits,
            &mut || false,
        )
    };
    let mut renderer = make(limits).unwrap();
    let sentinel = StationaryCompactSample {
        pressure_pa: 123.0,
        ..Default::default()
    };
    let mut output = [sentinel; 9];
    let mut polls = 0;
    assert!(
        renderer
            .field_block_controlled(&mut output, &mut || {
                polls += 1;
                polls == 7
            })
            .is_err()
    );
    assert_eq!(output, [sentinel; 9]);
    assert_eq!(renderer.samples_rendered(), 0);
    renderer.field_block(&mut output).unwrap();
    let mut split = make(limits).unwrap();
    let mut first = [StationaryCompactSample::default(); 3];
    split.field_block(&mut first).unwrap();
    assert_eq!(&first, &output[..3]);
    let mut rest = [0.0; 6];
    PressureRenderer::block(&mut split, &mut rest).unwrap();
    assert_eq!(
        rest,
        core::array::from_fn(|i| output[i + 3].microphone_pressure_pa)
    );
    assert_eq!(split.samples_rendered(), 9);
    limits.maximum_source_receiver_pairs = 17;
    assert!(make(limits).is_err());

    // Direct sound is covered, but the longer image path needs earlier data.
    let short = reflected_history(&[0.8, 1.05, 1.3], 1.0);
    let short_sources = [source(&short, [1.0, 0.0, 0.0], [0.0; 3])];
    let short_initial = [impulse(0.8, 1.0)];
    assert!(
        StationaryCompactRenderer::new(
            *images.observer(),
            &short_sources,
            &short_initial,
            microphone,
            reflected_config(),
            budget(),
            &mut || false,
        )
        .is_ok()
    );
    assert!(
        StationaryCompactRenderer::with_ideal_reflections(
            &images,
            &short_sources,
            &short_initial,
            microphone,
            reflected_config(),
            budget(),
            &mut || false,
        )
        .is_err()
    );
    for x in [-1.0, 0.0] {
        let outside = [source(&accepted, [x, 0.0, 0.0], [0.0; 3])];
        assert!(
            StationaryCompactRenderer::with_ideal_reflections(
                &images,
                &outside,
                &initial,
                microphone,
                reflected_config(),
                budget(),
                &mut || false,
            )
            .is_err()
        );
    }
}
