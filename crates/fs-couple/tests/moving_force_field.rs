//! Full moving force histories: physical near fields, accepted geometry and audio.

use fs_aeroac::curle3d::{CompactForceObserver, IdealAcousticPlane, IdealPlaneBoundary};
use fs_aeroac::curle3d::moving::{
    CompactSourceStrength, MovingCompactLimits, MovingCompactObserver, SourceKinematics,
    TimedCompactSource, UniformSourceMotion,
};
use fs_aeroac::curle3d::moving::trajectory::{
    IdealReflectionObserver, PreparedReflectedSourceTrajectory, PreparedSourceTrajectory,
    ReflectedTrajectoryBudget, TimedSourceKinematics, TrajectoryBudget,
};
use fs_bem::near_field::FirstOrder;
use fs_couple::pcm_wav::observation::moving::field::CompactFieldSample;
use fs_couple::pcm_wav::observation::moving::force::{
    ForceHistoryQuadrature, MovingForceRenderer,
};
use fs_couple::pcm_wav::observation::moving::{
    CompactSourceHistory, MovingCompactRenderer, MovingPressureBudget, MovingPressureConfig,
};
use fs_material::fluid_acoustics::FluidAcousticState;

fn medium() -> FluidAcousticState {
    FluidAcousticState::try_new(1.2, 300.0).unwrap()
}
fn observer(position: [f64; 3]) -> MovingCompactObserver {
    MovingCompactObserver::new(
        medium(), position,
        MovingCompactLimits { maximum_mach: 0.9, minimum_distance_m: 0.01 },
    ).unwrap()
}
fn config() -> MovingPressureConfig {
    MovingPressureConfig {
        first_receiver_time_s: 0.08, sample_rate_hz: 1000, samples: 33, max_block: 33,
        receiver_bandwidth_hz: 20.0, maximum_time_residual_s: 1e-12,
    }
}
fn budget() -> MovingPressureBudget {
    MovingPressureBudget {
        maximum_sources: 8, maximum_source_samples: 2048, maximum_output_samples: 4096,
        maximum_payload_bytes: 1_000_000, maximum_source_receiver_pairs: 65_536,
    }
}
fn quadrature() -> ForceHistoryQuadrature {
    ForceHistoryQuadrature {
        absolute_tolerance_n_s_m3: 1e-13,
        maximum_subdivision_depth: 20,
        maximum_evaluations: 1_000_000,
    }
}
fn trajectory_budget() -> TrajectoryBudget {
    TrajectoryBudget {
        maximum_segments: 16, maximum_payload_bytes: 131_072, maximum_bound_nodes: 8192,
        maximum_subdivision_depth: 16, maximum_root_iterations: 64,
    }
}
fn force(t: f64, direction: [f64; 3]) -> CompactSourceStrength {
    CompactSourceStrength {
        force_on_fluid_n: direction.map(|d| d * (t * t + 0.3 * t * t * t)),
        force_rate_on_fluid_n_s: direction.map(|d| d * (2.0 * t + 0.9 * t * t)),
        ..Default::default()
    }
}
fn history(direction: [f64; 3]) -> [TimedCompactSource; 4] {
    [0.0, 0.1, 0.2, 0.3].map(|source_time_s| TimedCompactSource {
        source_time_s, strength: force(source_time_s, direction),
    })
}
fn close(actual: f64, expected: f64, absolute: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * expected.abs(),
        "actual={actual:e}, expected={expected:e}",
    );
}

#[test]
fn stationary_transverse_pressure_null_retains_impulse_velocity_and_directional_audio() {
    let samples = history([0.0, 1.0, 0.0]);
    let source = CompactSourceHistory {
        motion: UniformSourceMotion::new(0.0, [0.0; 3], [0.0; 3]).unwrap(),
        samples: &samples, source_bandwidth_hz: 1.0,
    };
    let mut renderer = MovingForceRenderer::new_with_zero_force_prehistory(
        observer([12.0, 0.0, 0.0]), &[source],
        FirstOrder::new(0.0, [0.0, 1.0, 0.0]).unwrap(),
        config(), budget(), quadrature(), &mut || false,
    ).unwrap();
    let oracle = CompactForceObserver::new(medium(), [12.0, 0.0, 0.0], [0.0; 3]).unwrap();
    let mut fields = [CompactFieldSample::default(); 33];
    renderer.field_block(&mut fields).unwrap();
    for field in fields {
        let t = field.receiver_time_s - oracle.delay_s();
        let strength = force(t, [0.0, 1.0, 0.0]);
        let impulse = [0.0, t.powi(3) / 3.0 + 0.075 * t.powi(4), 0.0];
        let expected = oracle.field_at_retarded(
            impulse, strength.force_on_fluid_n, strength.force_rate_on_fluid_n_s,
        ).unwrap();
        assert_eq!(field.pressure_pa, 0.0);
        assert_eq!(field.instantaneous_intensity_w_m2, [0.0; 3]);
        assert!(field.particle_velocity_m_s[1] < 0.0);
        assert!(field.microphone_pressure_pa > 0.0);
        for axis in 0..3 {
            close(field.particle_velocity_m_s[axis], expected.particle_velocity_m_s[axis],
                1e-20, 2e-11);
        }
        close(field.microphone_pressure_pa, -360.0 * expected.particle_velocity_m_s[1],
            1e-18, 2e-11);
    }
}

fn curved(t: f64) -> SourceKinematics {
    SourceKinematics {
        position_m: [20.0 * t + 60.0 * t * t, 1.0 + 4.0 * t - 10.0 * t * t, 0.2],
        velocity_m_s: [20.0 + 120.0 * t, 4.0 - 20.0 * t, 0.0],
        acceleration_m_s2: [120.0, -20.0, 0.0],
    }
}
fn norm(v: [f64; 3]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}
// Independent original analytic motion/force and a dense fixed Simpson rule.
// It neither reads compiled prefixes nor uses the production adaptive scheme.
fn reference_history(position: [f64; 3], end: f64) -> [f64; 3] {
    let steps = 4096;
    let h = end / steps as f64;
    let mut sum = [0.0; 3];
    for i in 0..=steps {
        let time = i as f64 * h;
        let source = curved(time);
        let r: [f64; 3] = core::array::from_fn(|j| position[j] - source.position_m[j]);
        let distance = norm(r);
        let n = r.map(|x| x / distance);
        let force = force(time, [1.0, -0.4, 0.7]).force_on_fluid_n;
        let radial: f64 = (0..3).map(|j| n[j] * force[j]).sum();
        let weight = if i == 0 || i == steps { 1.0 } else if i % 2 == 0 { 2.0 } else { 4.0 };
        for j in 0..3 {
            sum[j] += weight * (3.0 * n[j] * radial - force[j]) / distance.powi(3);
        }
    }
    sum.map(|x| x * h / 3.0)
}
fn reference_emission(position: [f64; 3], receiver_time: f64) -> f64 {
    let mut lo: f64 = 0.0;
    let mut hi: f64 = 0.3;
    for _ in 0..80 {
        let t = lo.midpoint(hi);
        let source = curved(t);
        let distance = norm(core::array::from_fn(|j| position[j] - source.position_m[j]));
        if t + distance / 300.0 < receiver_time { lo = t; } else { hi = t; }
    }
    lo.midpoint(hi)
}

#[test]
fn curved_force_prefixes_split_geometry_knots_inside_one_original_force_interval() {
    let position = [6.0, 2.0, 1.0];
    let microphone = observer(position);
    let knots = [0.0, 0.03, 0.07, 0.11, 0.18, 0.3].map(|source_time_s| {
        TimedSourceKinematics { source_time_s, kinematics: curved(source_time_s) }
    });
    let trajectory = PreparedSourceTrajectory::new(
        microphone, &knots, trajectory_budget(), &mut || false,
    ).unwrap();
    let samples = [0.0, 0.3].map(|source_time_s| TimedCompactSource {
        source_time_s, strength: force(source_time_s, [1.0, -0.4, 0.7]),
    });
    let source = CompactSourceHistory {
        motion: &trajectory, samples: &samples, source_bandwidth_hz: 1.0,
    };
    let clock = MovingPressureConfig { first_receiver_time_s: 0.12, ..config() };
    let mut pressure_owner = MovingCompactRenderer::new(
        microphone, &[source], clock, budget(), &mut || false,
    ).unwrap();
    let mut original_pressure = [0.0; 33];
    pressure_owner.block(&mut original_pressure).unwrap();
    let mut errors = [0.0_f64; 3];
    for (rung, tolerance) in [1e-6, 1e-10, 1e-14].into_iter().enumerate() {
        let mut renderer = MovingForceRenderer::new_with_zero_force_prehistory(
            microphone, &[source], FirstOrder::default(), clock, budget(),
            ForceHistoryQuadrature { absolute_tolerance_n_s_m3: tolerance, ..quadrature() },
            &mut || false,
        ).unwrap();
        let mut fields = [CompactFieldSample::default(); 33];
        renderer.field_block(&mut fields).unwrap();
        for (i, field) in fields.iter().enumerate() {
            close(field.pressure_pa, original_pressure[i], 1e-18, 2e-12);
            if i % 8 == 0 {
                let time = reference_emission(position, field.receiver_time_s);
                let expected = microphone.prepare_force_emission(time, curved(time)).unwrap()
                    .evaluate(force(time, [1.0, -0.4, 0.7]), reference_history(position, time))
                    .unwrap();
                for axis in 0..3 {
                    errors[rung] = errors[rung].max(
                        (field.particle_velocity_m_s[axis] - expected.particle_velocity_m_s[axis]).abs(),
                    );
                    close(field.instantaneous_intensity_w_m2[axis],
                        field.pressure_pa * field.particle_velocity_m_s[axis], 0.0, 0.0);
                }
            }
        }
    }
    assert!(errors[2] <= errors[0] + 1e-14, "{errors:?}");
    assert!(errors[2] < 2e-12, "{errors:?}");
}

#[test]
fn fixed_images_transform_full_moving_force_velocity_and_pressure_once() {
    let physical = observer([0.0, 3.0, 1.0]);
    let knots = [0.0, 0.1, 0.3].map(|source_time_s| TimedSourceKinematics {
        source_time_s,
        kinematics: SourceKinematics {
            position_m: [2.0, 1.0 + 2.0 * source_time_s, 0.5],
            velocity_m_s: [0.0, 2.0, 0.0],
            acceleration_m_s2: [0.0; 3],
        },
    });
    let trajectory = PreparedSourceTrajectory::new(
        physical, &knots, trajectory_budget(), &mut || false,
    ).unwrap();
    let samples = history([1.0, 0.4, -0.3]);
    for boundary in [IdealPlaneBoundary::Rigid, IdealPlaneBoundary::PressureRelease] {
        let walls = IdealReflectionObserver::new(
            physical,
            &[IdealAcousticPlane::new([0.0; 3], [1.0, 0.0, 0.0], boundary).unwrap()],
            2,
        ).unwrap();
        let reflected = PreparedReflectedSourceTrajectory::new(
            &trajectory, &walls,
            ReflectedTrajectoryBudget {
                maximum_payload_bytes: 1_000_000,
                maximum_bound_nodes: 4096, maximum_subdivision_depth: 16,
            },
            &mut || false,
        ).unwrap();
        let source = CompactSourceHistory {
            motion: &reflected, samples: &samples, source_bandwidth_hz: 1.0,
        };
        let mut renderer = MovingForceRenderer::new_with_zero_force_prehistory(
            &walls, &[source], FirstOrder::new(0.0, [1.0, 0.0, 0.0]).unwrap(),
            config(), budget(), quadrature(), &mut || false,
        ).unwrap();
        assert_eq!(renderer.source_count(), 1);
        assert_eq!(renderer.source_path_count(), 2);
        let mut fields = [CompactFieldSample::default(); 33];
        renderer.field_block(&mut fields).unwrap();
        for field in fields {
            match boundary {
                IdealPlaneBoundary::Rigid => {
                    assert!(field.pressure_pa.abs() > 1e-10);
                    assert_eq!(field.particle_velocity_m_s[0], 0.0);
                    assert_eq!(field.microphone_pressure_pa, 0.0);
                }
                IdealPlaneBoundary::PressureRelease => {
                    assert_eq!(field.pressure_pa, 0.0);
                    assert_eq!(&field.particle_velocity_m_s[1..], &[0.0; 2]);
                    assert!(field.particle_velocity_m_s[0].abs() > 1e-10);
                    close(field.microphone_pressure_pa, -360.0 * field.particle_velocity_m_s[0],
                        1e-18, 2e-12);
                }
            }
            assert_eq!(field.instantaneous_intensity_w_m2[0], 0.0);
        }
    }
}

#[test]
fn fixed_prefixes_preserve_block_partition_and_cancelled_work_is_not_refunded() {
    let samples = history([1.0, 0.4, -0.3]);
    let source = CompactSourceHistory {
        motion: UniformSourceMotion::new(0.0, [0.0; 3], [30.0, 2.0, 0.0]).unwrap(),
        samples: &samples, source_bandwidth_hz: 1.0,
    };
    let make = || MovingForceRenderer::new_with_zero_force_prehistory(
        observer([12.0, 2.0, 1.0]), &[source], FirstOrder::new(0.5, [1.0, 0.0, 0.0]).unwrap(),
        config(), budget(), quadrature(), &mut || false,
    ).unwrap();
    let mut whole = make();
    let mut expected = [CompactFieldSample::default(); 33];
    whole.field_block(&mut expected).unwrap();
    let mut split = make();
    let mut actual = [CompactFieldSample::default(); 33];
    for chunk in actual.chunks_mut(7) { split.field_block(chunk).unwrap(); }
    assert_eq!(actual, expected);
    assert_eq!(split.quadrature_evaluations(), whole.quadrature_evaluations());

    let mut retry = make();
    let before = retry.quadrature_evaluations();
    let sentinel = CompactFieldSample { pressure_pa: 123.0, ..Default::default() };
    let mut cancelled_output = [sentinel; 33];
    let mut polls = 0;
    assert!(retry.field_block_controlled(&mut cancelled_output, &mut || {
        polls += 1; polls == 25
    }).is_err());
    assert_eq!(cancelled_output, [sentinel; 33]);
    assert_eq!(retry.samples_rendered(), 0);
    assert!(retry.quadrature_evaluations() > before);
    retry.field_block(&mut cancelled_output).unwrap();
    assert_eq!(cancelled_output, expected);
    assert!(retry.quadrature_evaluations() > whole.quadrature_evaluations());
}

#[test]
fn causal_onset_complete_past_clearance_and_lifetime_work_are_real_admission_gates() {
    let samples = [0.0, 1.0].map(|source_time_s| TimedCompactSource {
        source_time_s, strength: force(source_time_s, [1.0, 0.0, 0.0]),
    });
    let source = CompactSourceHistory {
        motion: UniformSourceMotion::new(0.0, [0.0; 3], [20.0, 0.0, 0.0]).unwrap(),
        samples: &samples, source_bandwidth_hz: 0.1,
    };
    let late = MovingPressureConfig { first_receiver_time_s: 0.7, ..config() };
    // Requested retarded events are clear, but earlier history crosses x=1.
    assert!(MovingCompactRenderer::new(
        observer([1.0, 0.0, 0.0]), &[source], late, budget(), &mut || false,
    ).is_ok());
    assert!(MovingForceRenderer::new_with_zero_force_prehistory(
        observer([1.0, 0.0, 0.0]), &[source], FirstOrder::default(),
        late, budget(), quadrature(), &mut || false,
    ).is_err());

    let clear = CompactSourceHistory {
        motion: UniformSourceMotion::new(0.0, [0.0; 3], [0.0; 3]).unwrap(), ..source
    };
    for nonzero_rate in [false, true] {
        let mut invalid = samples;
        if nonzero_rate { invalid[0].strength.force_rate_on_fluid_n_s[0] = 1.0; }
        else { invalid[0].strength.force_on_fluid_n[0] = 1.0; }
        assert!(MovingForceRenderer::new_with_zero_force_prehistory(
            observer([12.0, 0.0, 0.0]), &[CompactSourceHistory { samples: &invalid, ..clear }],
            FirstOrder::default(), config(), budget(), quadrature(), &mut || false,
        ).is_err());
    }

    // No original knot precedes the render window, so preparation stores zero
    // without spending integrand work; four evaluations cannot finish Simpson.
    let mut limited = MovingForceRenderer::new_with_zero_force_prehistory(
        observer([12.0, 0.0, 0.0]), &[clear], FirstOrder::default(),
        config(), budget(),
        ForceHistoryQuadrature { maximum_evaluations: 4, ..quadrature() },
        &mut || false,
    ).unwrap();
    let mut output = [123.0; 1];
    assert!(limited.block(&mut output).is_err());
    assert_eq!(limited.quadrature_evaluations(), 4);
    assert_eq!(limited.samples_rendered(), 0);
    assert_eq!(output, [123.0]);
    assert!(limited.block(&mut output).is_err());
    assert_eq!(limited.quadrature_evaluations(), 4);
}

#[test]
fn geometric_refinement_recovers_a_narrow_transverse_close_pass_missed_by_five_nodes() {
    let speed: f64 = 280.0;
    let closest_time: f64 = 0.317;
    let clearance: f64 = 0.01;
    let amplitude: f64 = 1e-6;
    let microphone = MovingCompactObserver::new(
        medium(), [0.0; 3],
        MovingCompactLimits { maximum_mach: 0.99, minimum_distance_m: 0.001 },
    ).unwrap();
    let motion = UniformSourceMotion::new(
        0.0, [-speed * closest_time, clearance, 0.0], [speed, 0.0, 0.0],
    ).unwrap();
    let samples = [0.0, 1.0].map(|source_time_s| TimedCompactSource {
        source_time_s,
        strength: CompactSourceStrength {
            force_on_fluid_n: [0.0, 0.0, amplitude * source_time_s * source_time_s],
            force_rate_on_fluid_n_s: [0.0, 0.0, 2.0 * amplitude * source_time_s],
            ..Default::default()
        },
    });
    let last = microphone.prepare_force_emission(1.0, motion.kinematics_at(1.0).unwrap()).unwrap();
    let clock = MovingPressureConfig {
        first_receiver_time_s: last.pressure_sample().receiver_time_s(),
        samples: 1, max_block: 1, ..config()
    };
    let mut renderer = MovingForceRenderer::new_with_zero_force_prehistory(
        microphone,
        &[CompactSourceHistory { motion, samples: &samples, source_bandwidth_hz: 0.1 }],
        FirstOrder::new(0.0, [0.0, 0.0, 1.0]).unwrap(), clock, budget(),
        ForceHistoryQuadrature {
            absolute_tolerance_n_s_m3: 1e-9, maximum_subdivision_depth: 24, ..quadrature()
        },
        &mut || false,
    ).unwrap();
    // Analytic integral of s²/[v²(s-a)²+b²]^(3/2), without quadrature.
    let primitive = |s: f64| {
        let u = speed * (s - closest_time);
        let r = u.hypot(clearance);
        (closest_time * closest_time * u / (clearance * clearance * r)
            - (2.0 * closest_time / speed) / r
            + ((u / clearance).asinh() - u / r) / (speed * speed)) / speed
    };
    let j = -amplitude * (primitive(1.0) - primitive(0.0));
    assert!(j.abs() > 7e-6);
    let expected = last.evaluate(samples[1].strength, [0.0, 0.0, j]).unwrap();
    let mut field = [CompactFieldSample::default(); 1];
    renderer.field_block(&mut field).unwrap();
    assert_eq!(field[0].pressure_pa, 0.0);
    assert!(field[0].particle_velocity_m_s[2].abs() > 4e-7);
    close(field[0].particle_velocity_m_s[2], expected.particle_velocity_m_s[2],
        5e-11, 1e-5);
    assert!(renderer.quadrature_evaluations() > 5);
}

#[test]
fn force_history_refuses_geometry_missing_before_the_requested_retarded_window() {
    let microphone = observer([12.0, 0.0, 0.0]);
    let knots = [0.02, 0.3].map(|source_time_s| TimedSourceKinematics {
        source_time_s,
        kinematics: SourceKinematics {
            position_m: [0.0; 3], velocity_m_s: [0.0; 3], acceleration_m_s2: [0.0; 3],
        },
    });
    let trajectory = PreparedSourceTrajectory::new(
        microphone, &knots, trajectory_budget(), &mut || false,
    ).unwrap();
    let samples = history([0.0, 1.0, 0.0]);
    let source = CompactSourceHistory {
        motion: &trajectory, samples: &samples, source_bandwidth_hz: 1.0,
    };
    assert!(MovingCompactRenderer::new(
        microphone, &[source], config(), budget(), &mut || false,
    ).is_ok());
    assert!(MovingForceRenderer::new_with_zero_force_prehistory(
        microphone, &[source], FirstOrder::default(), config(), budget(), quadrature(),
        &mut || false,
    ).is_err());
}
