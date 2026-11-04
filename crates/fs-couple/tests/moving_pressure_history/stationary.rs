//! G1/G3/G4: actual accepted source amplitudes through vector fields and PCM.
use super::*;
use fs_aeroac::curle3d::CompactForceObserver;
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
