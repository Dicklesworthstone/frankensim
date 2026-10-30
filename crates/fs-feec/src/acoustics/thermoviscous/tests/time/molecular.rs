use super::*;
use crate::acoustics::thermoviscous::time::ThermoviscousCaloricReservoir;

fn reservoirs() -> [ThermoviscousCaloricReservoir; 3] {
    [
        ThermoviscousCaloricReservoir {
            heat_capacity_j_kg_k: 110.0,
            relaxation_time_s: 0.001,
        },
        // Keep the original source ordering, without a singular ghost state.
        ThermoviscousCaloricReservoir {
            heat_capacity_j_kg_k: 0.0,
            relaxation_time_s: 0.0004,
        },
        ThermoviscousCaloricReservoir {
            heat_capacity_j_kg_k: 70.0,
            relaxation_time_s: 0.0002,
        },
    ]
}

#[test]
fn g1_isothermal_translational_wall_does_not_erase_internal_molecular_state() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let pos = [
        [0.0; 3],
        [0.002, 0.0, 0.0],
        [0.0, 0.001, 0.0],
        [0.0, 0.0, 0.001],
    ];
    let boundary = ThermoviscousBoundary {
        no_slip_faces: mesh.faces.clone(),
        isothermal_faces: mesh.faces.clone(),
        ports: vec![],
    };
    let model =
        ThermoviscousFe::new(&mesh, &pos, air(), &boundary, limits(), &mut || false).unwrap();
    assert!(model.temperature_map.iter().all(Option::is_none));
    let reservoir = reservoirs()[0];
    let dt = 0.4 * reservoir.relaxation_time_s;
    let prepared = model
        .prepare_time_with_caloric_reservoirs(
            dt,
            air().cp_j_kg_k - reservoir.heat_capacity_j_kg_k,
            &[reservoir],
            ThermoviscousTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let mut initial = prepared.zero_state();
    assert_eq!(initial.molecular_temperature_k[0].len(), model.space.ndof);
    initial.molecular_temperature_k[0][0] = 1.0;
    initial.molecular_temperature_k[0][model.space.ndof - 1] = -0.5;
    let original = initial.clone();
    let energy = prepared.stored_energy_j(&initial).unwrap();
    assert!(energy > 0.0);
    let frame = prepared.step(&initial, &[], &[], &mut || false).unwrap();
    assert_eq!(initial, original);
    let decay = det::exp(-dt / reservoir.relaxation_time_s);
    for (before, after) in initial.molecular_temperature_k[0]
        .iter()
        .zip(&frame.state.molecular_temperature_k[0])
    {
        assert!((after - decay * before).abs() < 2e-12);
    }
    assert!(frame.state.temperature_k.iter().all(|v| *v == 0.0));
    assert!(frame.state.pressure_pa.iter().all(|v| v.abs() < 1e-12));
    assert!(relative(frame.work.final_energy_j, energy * decay * decay) < 2e-11);
    assert!(
        relative(
            frame.work.molecular_relaxation_heat_j,
            energy * (1.0 - decay * decay)
        ) < 2e-11
    );
    assert!(frame.work.balance_residual_j().abs() < 1e-11 * energy);
    assert_eq!(frame.work.fluid_conductive_heat_j, 0.0);
    assert_eq!(frame.work.fluid_viscous_heat_j, 0.0);
}

#[test]
fn g3_original_complex_cp_fields_and_loss_reappear_as_the_force_clock_refines() {
    let size = [0.004, 0.0008, 0.001];
    let (mesh, pos) = box_mesh([1, 1, 1], size);
    let spectrum = reservoirs();
    let frozen = air().cp_j_kg_k - 180.0;
    let frequency = 400.0;
    let omega = std::f64::consts::TAU * frequency;
    let mut cp = C64::from_re(frozen);
    for r in &spectrum {
        cp = cp
            + C64::new(1.0, omega * r.relaxation_time_s)
                .recip()
                .scale(r.heat_capacity_j_kg_k);
    }
    let drive = [C64::new(1.0, 0.2), C64::new(0.5, -0.3)];
    for isothermal in [false, true] {
        let boundary = boundaries(&mesh, &pos, size, isothermal);
        let model =
            ThermoviscousFe::new(&mesh, &pos, air(), &boundary, limits(), &mut || false).unwrap();
        let harmonic = model
            .prepare_with_heat_capacity(frequency, cp, &mut || false)
            .unwrap()
            .response(&drive, &mut || false)
            .unwrap();
        let mut initial = phasor_state(&harmonic, omega);
        initial.molecular_temperature_k = spectrum
            .iter()
            .map(|r| {
                if r.heat_capacity_j_kg_k == 0.0 {
                    return vec![];
                }
                let lag = C64::new(1.0, omega * r.relaxation_time_s).recip();
                harmonic
                    .temperature_k
                    .iter()
                    .map(|t| (*t * lag).re)
                    .collect()
            })
            .collect();
        let pressure_scale = harmonic
            .pressure_pa
            .iter()
            .map(|p| p.abs())
            .fold(0.0, f64::max);
        let mut errors = Vec::new();
        for samples in [80, 320] {
            let dt = 1.0 / (frequency * samples as f64);
            let prepared = model
                .prepare_time_with_caloric_reservoirs(
                    dt,
                    frozen,
                    &spectrum,
                    ThermoviscousTimeLimits::default(),
                    &mut || false,
                )
                .unwrap();
            assert_eq!(prepared.caloric_reservoirs(), &spectrum);
            assert!(prepared.zero_state().molecular_temperature_k[1].is_empty());
            let mut state = initial.clone();
            let mut molecular_heat = 0.0;
            for k in 0..samples {
                let phase = omega * (k as f64 + 0.5) * dt;
                let rotation = C64::new(phase.cos(), phase.sin());
                let held = drive.map(|p| (p * rotation).re);
                let frame = prepared.step(&state, &held, &[], &mut || false).unwrap();
                molecular_heat += frame.work.molecular_relaxation_heat_j;
                state = frame.state;
            }
            let error = state
                .pressure_pa
                .iter()
                .zip(&initial.pressure_pa)
                .map(|(a, b)| (a - b).abs() / pressure_scale)
                .fold(0.0, f64::max);
            errors.push(error);
            if samples == 320 {
                assert!(relative(molecular_heat, harmonic.power.caloric_w / frequency) < 2e-3);
                for (actual, expected) in state
                    .molecular_temperature_k
                    .iter()
                    .zip(&initial.molecular_temperature_k)
                {
                    let scale = expected
                        .iter()
                        .map(|v| v.abs())
                        .fold(0.0, f64::max)
                        .max(1e-30);
                    assert!(
                        actual
                            .iter()
                            .zip(expected)
                            .all(|(a, b)| (a - b).abs() < 2e-3 * scale)
                    );
                }
            }
        }
        assert!(errors[1] < 0.2 * errors[0], "{errors:?}");
        assert!(errors[1] < 2e-4, "{errors:?}");
    }
}

#[test]
fn g0_original_spectrum_state_capacity_and_cancellation_are_transactional() {
    let size = [0.004, 0.0008, 0.001];
    let (mesh, pos) = box_mesh([1, 1, 1], size);
    let boundary = boundaries(&mesh, &pos, size, true);
    let model =
        ThermoviscousFe::new(&mesh, &pos, air(), &boundary, limits(), &mut || false).unwrap();
    let spectrum = reservoirs();
    let frozen = air().cp_j_kg_k - 180.0;
    let limits0 = ThermoviscousTimeLimits::default();
    assert!(
        model
            .prepare_time_with_caloric_reservoirs(
                1e-5,
                frozen + 1.0,
                &spectrum,
                limits0,
                &mut || false
            )
            .is_err()
    );
    let mut bad = spectrum;
    bad[0].relaxation_time_s = -1.0;
    assert!(
        model
            .prepare_time_with_caloric_reservoirs(1e-5, frozen, &bad, limits0, &mut || false)
            .is_err()
    );
    let small = ThermoviscousTimeLimits {
        max_states: model.n,
        ..limits0
    };
    assert!(
        model
            .prepare_time_with_caloric_reservoirs(1e-5, frozen, &spectrum, small, &mut || false)
            .is_err()
    );
    let prepared = model
        .prepare_time_with_caloric_reservoirs(1e-5, frozen, &spectrum, limits0, &mut || false)
        .unwrap();
    let mut initial = prepared.zero_state();
    initial.molecular_temperature_k[0][0] = 1.0;
    let original = initial.clone();
    let mut calls = 0;
    assert!(matches!(
        prepared.step(&initial, &[0.0, 0.0], &[], &mut || {
            calls += 1;
            calls == 2
        }),
        Err(ThermoviscousTimeError::Fluid(ThermoviscousError::Cancelled))
    ));
    assert_eq!(initial, original);
    initial.molecular_temperature_k[1].push(1.0);
    assert!(prepared.stored_energy_j(&initial).is_err());
}

#[test]
fn g1_zero_capacity_spectrum_preserves_the_original_full_field_trajectory() {
    let size = [0.004, 0.0008, 0.001];
    let (mesh, pos) = box_mesh([1, 1, 1], size);
    let boundary = boundaries(&mesh, &pos, size, true);
    let model =
        ThermoviscousFe::new(&mesh, &pos, air(), &boundary, limits(), &mut || false).unwrap();
    let plain = model
        .prepare_time(1e-5, ThermoviscousTimeLimits::default(), &mut || false)
        .unwrap();
    let empty = model
        .prepare_time_with_caloric_reservoirs(
            1e-5,
            air().cp_j_kg_k,
            &[reservoirs()[1]],
            ThermoviscousTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(plain.system().state_dim(), empty.system().state_dim());
    let a = plain
        .step(&plain.zero_state(), &[1.0, 0.0], &[], &mut || false)
        .unwrap();
    let b = empty
        .step(&empty.zero_state(), &[1.0, 0.0], &[], &mut || false)
        .unwrap();
    assert_eq!(a.state.pressure_pa, b.state.pressure_pa);
    assert_eq!(a.state.temperature_k, b.state.temperature_k);
    assert_eq!(a.state.velocity_m_s, b.state.velocity_m_s);
    assert_eq!(a.work.pressure_work_j, b.work.pressure_work_j);
    assert_eq!(b.work.molecular_relaxation_heat_j, 0.0);
    assert_eq!(b.state.molecular_temperature_k, vec![Vec::<f64>::new()]);
}

#[test]
fn g1_original_wall_inertia_and_molecular_heat_share_one_complete_state_and_work_window() {
    let size = [0.002, 0.001, 0.001];
    let volume = size.iter().product::<f64>();
    let density = 1000.0;
    let amplitude = 1.0 / det::sqrt(density * volume);
    // The existing zero-expansion liquid fixture has an exact translating
    // added-mass solution. Its independent molecular relaxation must neither
    // replace that inertia nor alias the appended wall displacement coordinate.
    let (mesh, pos) = box_mesh([1, 1, 1], size);
    let model = translating(&mesh, &pos, density, amplitude, 0.0, 0.0);
    let spectrum = reservoirs();
    let dt = 1e-4;
    let prepared = model
        .prepare_time_with_caloric_reservoirs(
            dt,
            model.medium.cp_j_kg_k - 180.0,
            &spectrum,
            ThermoviscousTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let mut initial = prepared.zero_state();
    let acceleration = 1.0 / (1.0 + density * volume * amplitude * amplitude);
    for (p, xyz) in initial.pressure_pa.iter_mut().zip(&pos) {
        *p = -density * amplitude * acceleration * (xyz[0] - 0.5 * size[0]);
    }
    initial.wall_displacement_m_sqrt_kg[0] = 1e-6;
    for t in &mut initial.molecular_temperature_k[0][..mesh.vertex_count] {
        *t = 1.0;
    }
    let original = initial.clone();
    let frame = prepared.step(&initial, &[], &[1.0], &mut || false).unwrap();
    assert_eq!(initial, original);
    assert!(
        relative(
            frame.state.wall_velocity_m_sqrt_kg_per_s[0],
            acceleration * dt
        ) < 2e-8
    );
    assert!(
        relative(
            frame.state.wall_displacement_m_sqrt_kg[0] - initial.wall_displacement_m_sqrt_kg[0],
            0.5 * acceleration * dt * dt,
        ) < 2e-8
    );
    assert!(
        relative(
            frame.mean_wall_velocity_m_sqrt_kg_per_s[0],
            0.5 * acceleration * dt
        ) < 2e-8
    );
    assert!(relative(frame.work.structural_work_j, 0.5 * acceleration * dt * dt) < 2e-8);
    assert!(frame.work.molecular_relaxation_heat_j > 0.0);
    assert!(frame.state.temperature_k.iter().any(|t| *t > 0.0));
    assert!(
        relative(
            prepared.stored_energy_j(&frame.state).unwrap(),
            frame.work.final_energy_j
        ) < 1e-13
    );
    assert!(frame.work.balance_residual_j().abs() < 1e-10 * frame.work.initial_energy_j);
}
