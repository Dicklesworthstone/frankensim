//! Volume sources excite the native field, finite ducts and shared pressure observer.
//! Coarse one-cell fixtures use adiabatic walls: constraining every vertex to
//! zero temperature would remove the P1 heating coordinates being exercised.
use super::super::super::dg::{
    ThermoviscousDgFrame, ThermoviscousVolumeDg,
    performance::{
        ThermoviscousDgControl, ThermoviscousDgPerformance, ThermoviscousDgPerformanceConfig,
        ThermoviscousDgVolumeControl, ThermoviscousVolumeDgPerformance,
    },
};
use super::*;
use crate::bernoulli_aperture::performance::ApertureObservation;
use crate::pcm_wav::observation::PressureRenderer;
use fs_feec::acoustics::thermoviscous::{
    body_force::ThermoviscousBodyForceSource, heat_source::ThermoviscousHeatSource,
};
use std::sync::Arc;

fn profiles<'a, 'c>(
    interior: &'a ThermoviscousInterior<'c>,
) -> (
    ThermoviscousBodyForceSource<'a, 'c>,
    ThermoviscousHeatSource<'a, 'c>,
) {
    let volume = LENGTH * HEIGHT * WIDTH;
    let force = interior
        .volume_force_source(
            4,
            &mut |_, x| {
                Ok([
                    (1.0 + x[0] / LENGTH) / volume,
                    -0.2 * x[1] / (HEIGHT * volume),
                    0.0,
                ])
            },
            &mut || false,
        )
        .unwrap();
    let heat = interior
        .volume_heat_source(4, &mut |_, _| Ok(1.0 / volume), &mut || false)
        .unwrap();
    (force, heat)
}

fn system<'a, 'c>(
    interior: &'a ThermoviscousInterior<'c>,
    force: &ThermoviscousBodyForceSource<'_, '_>,
    heat: &ThermoviscousHeatSource<'_, '_>,
    duration: f64,
    limits: ThermoviscousDgLimits,
    gate: &CancelGate,
) -> ThermoviscousVolumeDg<'a, 'c> {
    ThermoviscousDg::new(
        interior,
        interior
            .prepare_time_with_volume_sources(
                duration,
                &[force],
                &[heat],
                ThermoviscousTimeLimits::default(),
                &mut || false,
            )
            .unwrap(),
        time_dg_load(interior.gas(), duration, DgBoundary::Closed),
        1,
        limits,
        gate,
    )
    .unwrap()
}

fn near(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 2e-10 * scale.max(1e-30),
        "{actual:e} != {expected:e}, scale={scale:e}"
    );
}

fn source_balance(frame: &ThermoviscousDgFrame) {
    let external = frame.field.work.body_force_work_j + frame.field.work.heat_source_work_j;
    near(
        frame.authored_external_work_j,
        external,
        frame.field.work.body_force_work_j.abs() + frame.field.work.heat_source_work_j.abs(),
    );
    assert!(frame.balance_residual_j.abs() <= frame.balance_tolerance_j);
    assert!(frame.interface_work_residual_j.abs() <= frame.balance_tolerance_j);
    assert!(
        (frame.storage_change_j + frame.material_dissipation_j + frame.terminal_work_j - external)
            .abs()
            <= frame.balance_tolerance_j
    );
}

#[test]
fn original_force_and_heat_drive_actual_dg_fields_and_retained_source_observations() {
    let (mesh, cube) = kuhn_cube(1);
    let xyz = points(&cube);
    let gas = gas();
    let interior = ThermoviscousInterior::from_relaxation(
        build(&mesh, &xyz, &gas, false),
        source(gas, true, 1.0),
        BAND,
    )
    .unwrap();
    let (force, heat) = profiles(&interior);
    let dt = 1.0 / 96_000.0;
    let gate = CancelGate::new_clock_free();
    // Each physical mechanism must produce duct sound from rest without a
    // traction, structural force, initial pulse or additive output source.
    for [f, q] in [[1e-6, 0.0], [0.0, 1e-3]] {
        let mut coupled = system(
            &interior,
            &force,
            &heat,
            dt,
            ThermoviscousDgLimits::default(),
            &gate,
        );
        assert!(coupled.last_volume_observations().is_none());
        let frame = coupled
            .step_with_volume_sources(&[0.0; 2], &[], &[f], &[q], &gate)
            .unwrap()
            .clone();
        assert!(frame.mean_pressure_pa != 0.0 && frame.mean_flow_m3_s != 0.0);
        assert!(frame.duct.ledger.stored_energy_j > 0.0);
        assert!(coupled.load().state().pressure().iter().any(|p| *p != 0.0));
        source_balance(&frame);
    }
    let mut coupled = system(
        &interior,
        &force,
        &heat,
        dt,
        ThermoviscousDgLimits::default(),
        &gate,
    );
    for [f, q] in [[1e-6, 0.0], [0.0, 1e-3], [-5e-7, -3e-4], [0.0, 0.0]] {
        let initial = coupled.state().clone();
        let mut original_duct = coupled.load().clone();
        let frame = coupled
            .step_with_volume_sources(&[0.0; 2], &[], &[f], &[q], &gate)
            .unwrap()
            .clone();
        let replay = coupled
            .prepared()
            .step(
                &initial,
                &[0.0, frame.mean_pressure_pa],
                &[],
                &[f],
                &[q],
                &mut || false,
            )
            .unwrap();
        let duct = original_duct
            .preview_step(frame.mean_flow_m3_s)
            .unwrap()
            .accept();
        assert_eq!(frame.field.state, replay.field.state);
        assert_eq!(
            frame.field.mean_inward_flow_m3_s,
            replay.field.mean_inward_flow_m3_s
        );
        assert_eq!(frame.duct, duct);
        assert_eq!(coupled.load().state(), original_duct.state());
        let observed = coupled.last_volume_observations().unwrap();
        assert_eq!(
            observed.mean_body_force_velocity_m_s,
            replay.mean_body_force_velocity_m_s
        );
        assert_eq!(
            observed.body_force_work_by_source_j,
            replay.body_force_work_by_source_j
        );
        assert_eq!(
            observed.mean_heat_source_temperature_ratio,
            replay.mean_heat_source_temperature_ratio
        );
        assert_eq!(
            observed.heat_source_work_by_source_j,
            replay.heat_source_work_by_source_j
        );
        source_balance(&frame);
    }
    assert_eq!(coupled.accepted_steps(), 4);
    assert!(
        coupled
            .state()
            .molecular_temperature_k
            .iter()
            .flatten()
            .any(|t| *t != 0.0)
    );
}

#[test]
fn volume_sources_keep_the_full_two_duct_pressure_closure_and_every_native_history() {
    let (mesh, cube) = kuhn_cube(1);
    let xyz = points(&cube);
    let gas = gas();
    let interior =
        ThermoviscousInterior::from_gas(build(&mesh, &xyz, &gas, false), gas, BAND).unwrap();
    let (force, heat) = profiles(&interior);
    let dt = 1.0 / 96_000.0;
    let gate = CancelGate::new_clock_free();
    let mut coupled = system(
        &interior,
        &force,
        &heat,
        dt,
        ThermoviscousDgLimits::default(),
        &gate,
    )
    .with_additional_duct(0, time_dg_load(&gas, dt, DgBoundary::Matched), &gate)
    .unwrap();
    for [f, q] in [[1e-6, 1e-3], [-3e-7, 4e-4], [0.0, 0.0]] {
        let initial = coupled.state().clone();
        let mut ducts = [
            coupled.load().clone(),
            coupled.load_at_port(0).unwrap().clone(),
        ];
        let (free, mobility) = coupled
            .prepared()
            .average_relations(&initial, &[0.0; 2], &[], &[f], &[q], &[1, 0], &mut || false)
            .unwrap();
        assert!(mobility[1] != 0.0 && mobility[2] != 0.0);
        let relations = [
            ducts[0].port_response().unwrap(),
            ducts[1].port_response().unwrap(),
        ];
        let z = [
            relations[0].impedance_pa_s_m3,
            relations[1].impedance_pa_s_m3,
        ];
        // Independent explicit 2x2 elimination retains the off-diagonal
        // response; a pair of separate scalar closures does not satisfy it.
        let m = [
            1.0 + z[0] * mobility[0],
            z[0] * mobility[1],
            z[1] * mobility[2],
            1.0 + z[1] * mobility[3],
        ];
        let rhs = [
            relations[0].history_pressure_pa - z[0] * free[0],
            relations[1].history_pressure_pa - z[1] * free[1],
        ];
        let det = m[0] * m[3] - m[1] * m[2];
        let expected = [
            (rhs[0] * m[3] - m[1] * rhs[1]) / det,
            (m[0] * rhs[1] - rhs[0] * m[2]) / det,
        ];
        let frame = coupled
            .step_with_volume_sources(&[0.0; 2], &[], &[f], &[q], &gate)
            .unwrap()
            .clone();
        let pressures = [
            frame.mean_pressure_pa,
            frame.additional_ducts[0].mean_pressure_pa,
        ];
        let flows = [
            frame.mean_flow_m3_s,
            frame.additional_ducts[0].mean_flow_m3_s,
        ];
        let scale = pressures
            .iter()
            .chain(&expected)
            .map(|p| p.abs())
            .fold(1e-30, f64::max);
        for i in 0..2 {
            near(pressures[i], expected[i], scale);
        }
        let replay = coupled
            .prepared()
            .step(
                &initial,
                &[pressures[1], pressures[0]],
                &[],
                &[f],
                &[q],
                &mut || false,
            )
            .unwrap();
        assert_eq!(frame.field.state, replay.field.state);
        for (i, port) in [1, 0].into_iter().enumerate() {
            assert_eq!(
                flows[i].to_bits(),
                (-replay.field.mean_inward_flow_m3_s[port]).to_bits()
            );
            let original = ducts[i].preview_step(flows[i]).unwrap().accept();
            assert_eq!(
                coupled.load_at_port(port).unwrap().state(),
                ducts[i].state()
            );
            assert_eq!(
                if i == 0 {
                    &frame.duct
                } else {
                    &frame.additional_ducts[0].duct
                },
                &original
            );
        }
        source_balance(&frame);
    }
    assert_eq!(coupled.coupled_ports(), &[1, 0]);
    assert_eq!(coupled.load_at_port(0).unwrap().state().step_count(), 3);
}

#[test]
fn volume_pressure_playback_continues_both_ducts_and_complete_fields_across_blocks() {
    let (mesh, cube) = kuhn_cube(1);
    let xyz = points(&cube);
    let gas = gas();
    let interior = ThermoviscousInterior::from_relaxation(
        build(&mesh, &xyz, &gas, false),
        source(gas, true, 1.0),
        BAND,
    )
    .unwrap();
    let (force, heat) = profiles(&interior);
    let dt = 1.0 / 96_000.0;
    let gate = Arc::new(CancelGate::new_clock_free());
    let controls: Vec<_> = (0..48)
        .map(|i| ThermoviscousDgVolumeControl {
            traction_pa: vec![0.0; 2],
            wall_forces_n_per_sqrt_kg: vec![],
            force_drive_n: vec![if i < 8 {
                1e-6
            } else if i < 16 {
                -5e-7
            } else {
                0.0
            }],
            heat_drive_w: vec![if i < 5 {
                1e-3
            } else if (12..16).contains(&i) {
                -2e-4
            } else {
                0.0
            }],
        })
        .collect();
    let make = |maximum_block| {
        let coupled = system(
            &interior,
            &force,
            &heat,
            dt,
            ThermoviscousDgLimits::default(),
            &gate,
        )
        .with_additional_duct(0, time_dg_load(&gas, dt, DgBoundary::Matched), &gate)
        .unwrap();
        ThermoviscousVolumeDgPerformance::new_at_port(
            coupled,
            0,
            controls.clone(),
            ApertureObservation::Inlet,
            ThermoviscousDgPerformanceConfig {
                sample_rate_hz: 96_000,
                maximum_block,
                maximum_control_values: 4 * controls.len(),
            },
            Arc::clone(&gate),
        )
        .unwrap()
    };
    let mut whole = make(48);
    let mut divided = make(5);
    let mut a = [0.0; 48];
    let mut b = [0.0; 48];
    whole.block(&mut a).unwrap();
    for block in b.chunks_mut(5) {
        divided.block(block).unwrap();
    }
    assert_eq!(a, b);
    assert!(a.iter().any(|p| *p != 0.0));
    assert_eq!(whole.system().state(), divided.system().state());
    for port in [0, 1] {
        assert_eq!(
            whole.system().load_at_port(port).unwrap().state(),
            divided.system().load_at_port(port).unwrap().state()
        );
    }
    assert_eq!(
        whole.last_volume_observations(),
        divided.last_volume_observations()
    );
    assert_eq!(whole.controls(), controls.as_slice());
    assert_eq!(whole.samples_rendered(), 48);
    assert_eq!(whole.remaining_samples(), 0);
    assert_eq!(
        a[47].to_bits(),
        whole.last_frame().unwrap().additional_ducts[0]
            .mean_pressure_pa
            .to_bits()
    );
    assert!(whole.validate_sample_count(1).is_err());
}

#[test]
fn volume_source_playback_reaches_the_original_baffled_microphone_from_terminal_flows() {
    use crate::bernoulli_aperture::radiation::environment::multiport::source::ComputedBaffledInput;
    use crate::pcm_wav::baffled::CircularOutletReceiver;
    use fs_material::gas_species::input::read_gas_record;
    let exterior =
        read_gas_record(&mut "293.15 101325 0".split_ascii_whitespace(), None, None).unwrap();
    let gas = *exterior.state();
    let (mesh, cube) = kuhn_cube(1);
    let xyz = points(&cube);
    let interior =
        ThermoviscousInterior::from_gas(build(&mesh, &xyz, &gas, false), gas, BAND).unwrap();
    let (force, heat) = profiles(&interior);
    let gate = Arc::new(CancelGate::new_clock_free());
    let rate = 8000;
    let dt = 1.0 / f64::from(rate);
    let duct = DgDuct::new(
        &fs_duct::Duct {
            segments: vec![Segment::Cylinder {
                radius: (HEIGHT * WIDTH / std::f64::consts::PI).sqrt(),
                length: 0.05,
            }],
        },
        &gas,
        DgMesh {
            elements: 2,
            order: 2,
        },
        DgBoundary::Port(PortId(7)),
        DgBoundary::Port(PortId(1)),
    )
    .unwrap();
    let load = DgPortLoad::new_with_computed_baffled(
        duct,
        rate,
        1e-11,
        MAX_DG_PROPAGATOR_DOFS,
        &exterior,
        ComputedBaffledInput::parse("20 65 4 8 12 20 0.01 67108864 2 100 1e-8 64 1048576 257")
            .unwrap(),
        &gate,
    )
    .unwrap();
    let receiver = CircularOutletReceiver {
        position_m: [0.05, 0.02, 0.1],
        radial_rings: 4,
        angular_points: 16,
        maximum_frequency_hz: 1000.0,
    };
    let mut original_receiver = load.baffled_receiver(rate, receiver).unwrap();
    let coupled = ThermoviscousDg::new(
        &interior,
        interior
            .prepare_time_with_volume_sources(
                dt,
                &[&force],
                &[&heat],
                ThermoviscousTimeLimits::default(),
                &mut || false,
            )
            .unwrap(),
        load,
        1,
        ThermoviscousDgLimits::default(),
        &gate,
    )
    .unwrap();
    let controls: Vec<_> = (0..256)
        .map(|i| ThermoviscousDgVolumeControl {
            traction_pa: vec![0.0; 2],
            wall_forces_n_per_sqrt_kg: vec![],
            force_drive_n: vec![if i < 4 { 1e-6 } else { 0.0 }],
            heat_drive_w: vec![if i < 4 { 1e-3 } else { 0.0 }],
        })
        .collect();
    let mut playback = ThermoviscousVolumeDgPerformance::new(
        coupled,
        controls,
        ApertureObservation::TubeBaffled(receiver),
        ThermoviscousDgPerformanceConfig {
            sample_rate_hz: rate,
            maximum_block: 16,
            maximum_control_values: 4 * 256,
        },
        Arc::clone(&gate),
    )
    .unwrap();
    let mut audible_energy = 0.0;
    for _ in 0..256 {
        let mut sample = [0.0];
        playback.block(&mut sample).unwrap();
        let terminal = playback.system().load().last_terminal_frames();
        assert_eq!(terminal.len(), 1);
        let expected = original_receiver
            .step_trace(&[terminal[0].outward_flow_m3_s])
            .unwrap();
        assert_eq!(sample[0].to_bits(), expected.to_bits());
        audible_energy += sample[0] * sample[0];
        source_balance(playback.last_frame().unwrap());
    }
    assert!(audible_energy > 0.0);
    assert!(playback.baffled_receiver().is_some());
    assert_eq!(playback.samples_rendered(), 256);
}

#[test]
fn volume_controls_and_refused_joint_steps_preserve_accepted_histories_and_observations() {
    let (mesh, cube) = kuhn_cube(1);
    let xyz = points(&cube);
    let gas = gas();
    let interior =
        ThermoviscousInterior::from_gas(build(&mesh, &xyz, &gas, false), gas, BAND).unwrap();
    let (force, heat) = profiles(&interior);
    let dt = 1.0 / 96_000.0;
    let gate = Arc::new(CancelGate::new_clock_free());
    let make = |limits| {
        system(&interior, &force, &heat, dt, limits, &gate)
            .with_additional_duct(0, time_dg_load(&gas, dt, DgBoundary::Matched), &gate)
            .unwrap()
    };
    let mut coupled = make(ThermoviscousDgLimits {
        maximum_total_energy_j: 1e-12,
        ..ThermoviscousDgLimits::default()
    });
    coupled
        .step_with_volume_sources(&[0.0; 2], &[], &[1e-8], &[1e-6], &gate)
        .unwrap();
    let field = coupled.state().clone();
    let ducts = [
        coupled.load().state().clone(),
        coupled.load_at_port(0).unwrap().state().clone(),
    ];
    let observed = coupled.last_volume_observations().cloned();
    for (f, q) in [
        (&[][..], &[0.0][..]),
        (&[f64::NAN][..], &[0.0][..]),
        (&[0.1][..], &[1000.0][..]),
    ] {
        assert!(
            coupled
                .step_with_volume_sources(&[0.0; 2], &[], f, q, &gate)
                .is_err()
        );
        assert_eq!(coupled.state(), &field);
        assert_eq!(coupled.load().state(), &ducts[0]);
        assert_eq!(coupled.load_at_port(0).unwrap().state(), &ducts[1]);
        assert_eq!(coupled.last_volume_observations(), observed.as_ref());
        assert_eq!(coupled.accepted_steps(), 1);
    }
    let valid = ThermoviscousDgVolumeControl {
        traction_pa: vec![0.0; 2],
        wall_forces_n_per_sqrt_kg: vec![],
        force_drive_n: vec![1e-6],
        heat_drive_w: vec![1e-3],
    };
    let config = ThermoviscousDgPerformanceConfig {
        sample_rate_hz: 96_000,
        maximum_block: 8,
        maximum_control_values: 16,
    };
    for bad in [
        ThermoviscousDgVolumeControl {
            force_drive_n: vec![],
            ..valid.clone()
        },
        ThermoviscousDgVolumeControl {
            heat_drive_w: vec![f64::INFINITY],
            ..valid.clone()
        },
        ThermoviscousDgVolumeControl {
            traction_pa: vec![1.0, 0.0],
            ..valid.clone()
        },
    ] {
        assert!(
            ThermoviscousVolumeDgPerformance::new(
                make(ThermoviscousDgLimits::default()),
                vec![bad],
                ApertureObservation::Inlet,
                config,
                Arc::clone(&gate),
            )
            .is_err()
        );
    }
    assert!(
        ThermoviscousVolumeDgPerformance::new(
            make(ThermoviscousDgLimits::default()),
            vec![valid.clone()],
            ApertureObservation::Inlet,
            ThermoviscousDgPerformanceConfig {
                maximum_control_values: 3,
                ..config
            },
            Arc::clone(&gate),
        )
        .is_err()
    );
    // Legacy controls cannot silently omit the prepared force and heat columns.
    assert!(
        ThermoviscousDgPerformance::new(
            make(ThermoviscousDgLimits::default()),
            vec![ThermoviscousDgControl {
                traction_pa: vec![0.0; 2],
                wall_forces_n_per_sqrt_kg: vec![]
            }],
            ApertureObservation::Inlet,
            config,
            Arc::clone(&gate),
        )
        .is_err()
    );
    let mut stream = ThermoviscousVolumeDgPerformance::new(
        make(ThermoviscousDgLimits::default()),
        vec![valid; 4],
        ApertureObservation::Inlet,
        config,
        Arc::clone(&gate),
    )
    .unwrap();
    stream.block(&mut [0.0; 2]).unwrap();
    let state = stream.system().state().clone();
    let observation = stream.last_volume_observations().cloned();
    gate.request();
    assert!(stream.block(&mut [0.0]).is_err());
    assert_eq!(stream.samples_rendered(), 2);
    assert_eq!(stream.system().accepted_steps(), 2);
    assert_eq!(stream.system().state(), &state);
    assert_eq!(stream.last_volume_observations(), observation.as_ref());
    assert!(stream.validate_sample_count(0).is_err());
}
