//! Original volume-force trajectories, work, complete state and source ownership.
use super::*;
use crate::acoustics::thermoviscous::body_force::{
    MAX_BODY_FORCE_SOURCES, ThermoviscousBodyForceSource,
};
use crate::acoustics::thermoviscous::heat_source::ThermoviscousHeatSource;
use crate::acoustics::thermoviscous::time::{
    PreparedThermoviscousVolumeTime, ThermoviscousCaloricReservoir,
};
use std::cell::Cell;

fn prepare<'a, 'c>(
    model: &'a ThermoviscousFe<'c>,
    duration: f64,
    forces: &[&ThermoviscousBodyForceSource<'_, '_>],
    heat: &[&ThermoviscousHeatSource<'_, '_>],
) -> PreparedThermoviscousVolumeTime<'a, 'c> {
    model
        .prepare_time_with_volume_sources(
            duration,
            model.medium().cp_j_kg_k,
            &[],
            None,
            None,
            forces,
            heat,
            ThermoviscousTimeLimits::default(),
            &mut || false,
        )
        .unwrap()
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-8 * actual.abs().max(expected.abs()).max(1e-30),
        "{actual:e} != {expected:e}"
    );
}

fn same(a: &[f64], b: &[f64]) {
    assert_eq!(a.len(), b.len());
    let scale = a.iter().chain(b).map(|v| v.abs()).fold(1e-12_f64, f64::max);
    assert!(a.iter().zip(b).all(|(x, y)| (x - y).abs() < 2e-8 * scale));
}

fn same_state(a: &ThermoviscousTimeState, b: &ThermoviscousTimeState) {
    same(&a.pressure_pa, &b.pressure_pa);
    same(
        &a.velocity_m_s.iter().flatten().copied().collect::<Vec<_>>(),
        &b.velocity_m_s.iter().flatten().copied().collect::<Vec<_>>(),
    );
    same(&a.temperature_k, &b.temperature_k);
    for (a, b) in [
        (&a.molecular_temperature_k, &b.molecular_temperature_k),
        (
            &a.conductive_wall_temperature_k,
            &b.conductive_wall_temperature_k,
        ),
    ] {
        assert_eq!(a.len(), b.len());
        for (a, b) in a.iter().zip(b) {
            same(a, b);
        }
    }
    same(
        &a.wall_displacement_m_sqrt_kg,
        &b.wall_displacement_m_sqrt_kg,
    );
    same(
        &a.wall_velocity_m_sqrt_kg_per_s,
        &b.wall_velocity_m_sqrt_kg_per_s,
    );
}

#[test]
fn g1_uniform_vector_force_uses_actual_fluid_mass_and_reverses_signed_work() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let positions = [
        [0.0, 0.0, 0.0],
        [0.01, 0.0, 0.0],
        [0.0, 0.008, 0.0],
        [0.0, 0.0, 0.006],
    ];
    let volume = 0.01 * 0.008 * 0.006 / 6.0;
    let gas = air();
    let model = ThermoviscousFe::new(
        &mesh,
        &positions,
        gas,
        &ThermoviscousBoundary {
            no_slip_faces: vec![],
            isothermal_faces: vec![],
            ports: vec![BoundaryMode::uniform(&mesh.faces)],
        },
        limits(),
        &mut || false,
    )
    .unwrap();
    let direction = [1.0, -2.0, 0.5];
    let calls = Cell::new(0);
    let source = model
        .volume_force_source(
            4,
            &mut |_, _| {
                calls.set(calls.get() + 1);
                Ok::<_, ThermoviscousError>(direction.map(|v| v / volume))
            },
            &mut || false,
        )
        .unwrap();
    let sampled = calls.get();
    let duration = 2e-5;
    let drive = 2e-6;
    let mass = gas.density_kg_m3 * volume;
    let prepared = prepare(&model, duration, &[&source], &[]);
    assert_eq!(prepared.body_force_source_count(), 1);
    assert_eq!(prepared.heat_source_count(), 0);
    assert_eq!(prepared.system().port_dim(), 2);
    let zero = prepared.zero_state();
    let original = zero.clone();
    let whole = prepared
        .step(&zero, &[0.0], &[], &[drive], &[], &mut || false)
        .unwrap();
    assert_eq!(zero, original);
    // Independent continuum solution: rho*v_dot=f, div(v)=strain(v)=0,
    // p=T'=0. Hierarchical P2 represents the uniform velocity at vertices.
    let speed = drive * duration / mass;
    for velocity in &whole.field.state.velocity_m_s[..mesh.vertex_count] {
        for axis in 0..3 {
            close(velocity[axis], direction[axis] * speed);
        }
    }
    for velocity in &whole.field.state.velocity_m_s[mesh.vertex_count..] {
        assert!(velocity.iter().all(|v| v.abs() < 2e-8 * speed));
    }
    let sound_speed = ThermoviscousThermodynamics::ideal_gas(gas)
        .unwrap()
        .sound_speed_m_s(gas)
        .unwrap();
    let pressure_scale = gas.density_kg_m3 * sound_speed * speed;
    assert!(
        whole
            .field
            .state
            .pressure_pa
            .iter()
            .all(|p| p.abs() < 2e-8 * pressure_scale)
    );
    assert!(
        whole
            .field
            .state
            .temperature_k
            .iter()
            .all(|t| { t.abs() < 2e-8 * gas.temperature_k / gas.pressure_pa * pressure_scale })
    );
    let direction2: f64 = direction.iter().map(|v| v * v).sum();
    let expected = 0.5 * mass * direction2 * speed * speed;
    close(whole.field.work.body_force_work_j, expected);
    close(whole.body_force_work_by_source_j[0], expected);
    close(
        whole.mean_body_force_velocity_m_s[0],
        0.5 * direction2 * speed,
    );
    close(whole.field.work.final_energy_j, expected);
    close(
        prepared.stored_energy_j(&whole.field.state).unwrap(),
        expected,
    );
    assert!(whole.field.work.dissipated_j().abs() < 2e-8 * expected);
    assert!(whole.field.work.heat_source_work_j == 0.0);

    let reversed = prepared
        .step(&whole.field.state, &[0.0], &[], &[-drive], &[], &mut || {
            false
        })
        .unwrap();
    close(reversed.field.work.body_force_work_j, -expected);
    assert!(
        reversed
            .field
            .state
            .velocity_m_s
            .iter()
            .flatten()
            .all(|v| v.abs() < 2e-8 * speed)
    );
    assert!(reversed.field.work.final_energy_j.abs() < 2e-8 * expected);

    // The original harmonic force owner has v_hat=F/(i*omega*rho*V).
    // This checks that time and harmonic source units/maps agree.
    let frequency = 123.0;
    let omega = std::f64::consts::TAU * frequency;
    let harmonic = model
        .prepare(frequency, &mut || false)
        .unwrap()
        .response_with_volume_sources(
            &[C64::ZERO],
            &[],
            &[(&source, C64::from_re(drive))],
            &[],
            &mut || false,
        )
        .unwrap();
    for (steady, causal) in harmonic.velocity_m_s[..mesh.vertex_count]
        .iter()
        .zip(&whole.field.state.velocity_m_s)
    {
        for axis in 0..3 {
            close(-omega * duration * steady[axis].im, causal[axis]);
            assert!(steady[axis].re.abs() < 2e-8 * steady[axis].abs());
        }
    }
    assert_eq!(
        calls.get(),
        sampled,
        "preparation/steps never resample a source"
    );
}

#[test]
fn g1_volume_force_moves_the_original_structure_and_added_fluid_mass() {
    let size = [0.002, 0.001, 0.001];
    let volume = size.iter().product::<f64>();
    let (mesh, positions) = box_mesh([1, 1, 1], size);
    let amplitude = 1.0 / det::sqrt(1000.0 * volume);
    let duration = 1e-4;
    let drive = 0.001;
    for density in [1000.0, 7000.0] {
        let model = translating(&mesh, &positions, density, amplitude, 0.0, 0.0);
        let source = model
            .volume_force_source(
                4,
                &mut |_, _| Ok::<_, ThermoviscousError>([1.0 / volume, 0.0, 0.0]),
                &mut || false,
            )
            .unwrap();
        let prepared = prepare(&model, duration, &[&source], &[]);
        let half = prepare(&model, 0.5 * duration, &[&source], &[]);
        let acceleration = amplitude * drive / (1.0 + density * volume * amplitude * amplitude);
        let gradient = drive / volume - density * amplitude * acceleration;
        let mut initial = prepared.zero_state();
        // Exact fixed linear pressure with uniform acceleration. No structural
        // force is authored: the volume force reaches it through the same fluid.
        for (p, point) in initial.pressure_pa.iter_mut().zip(&positions) {
            *p = gradient * (point[0] - 0.5 * size[0]);
        }
        let original = initial.clone();
        let whole = prepared
            .step(&initial, &[], &[0.0], &[drive], &[], &mut || false)
            .unwrap();
        let first = half
            .step(&initial, &[], &[0.0], &[drive], &[], &mut || false)
            .unwrap();
        let second = half
            .step(&first.field.state, &[], &[0.0], &[drive], &[], &mut || {
                false
            })
            .unwrap();
        assert_eq!(initial, original);
        same_state(&whole.field.state, &second.field.state);
        close(
            whole.field.state.wall_velocity_m_sqrt_kg_per_s[0],
            acceleration * duration,
        );
        close(
            whole.field.state.wall_displacement_m_sqrt_kg[0],
            0.5 * acceleration * duration * duration,
        );
        close(
            whole.mean_body_force_velocity_m_s[0],
            amplitude * 0.5 * acceleration * duration,
        );
        for velocity in &whole.field.state.velocity_m_s[..mesh.vertex_count] {
            close(velocity[0], amplitude * acceleration * duration);
        }
        let expected = 0.5 * drive * amplitude * acceleration * duration * duration;
        close(whole.field.work.body_force_work_j, expected);
        close(whole.field.work.delta_energy_j, expected);
        close(
            whole.field.work.body_force_work_j,
            first.field.work.body_force_work_j + second.field.work.body_force_work_j,
        );
        assert_eq!(whole.field.work.structural_work_j, 0.0);
        assert!(whole.field.work.dissipated_j().abs() < 2e-8 * expected);
    }
}

#[test]
fn g3_force_heat_and_molecular_fields_share_the_same_complete_semigroup() {
    let size = [0.004, 0.002, 0.002];
    let volume = size.iter().product::<f64>();
    let (mesh, positions) = box_mesh([1, 1, 1], size);
    let gas = air();
    let model = ThermoviscousFe::new(
        &mesh,
        &positions,
        gas,
        &boundaries(&mesh, &positions, size, false),
        limits(),
        &mut || false,
    )
    .unwrap();
    let force = model
        .volume_force_source(
            4,
            &mut |_, x| {
                Ok::<_, ThermoviscousError>([
                    (1.0 + x[0] / size[0]) / volume,
                    -x[1] / (size[1] * volume),
                    0.0,
                ])
            },
            &mut || false,
        )
        .unwrap();
    let heat = model
        .volume_heat_source(4, &mut |_, _| Ok(1.0 / volume), &mut || false)
        .unwrap();
    let reservoirs = [ThermoviscousCaloricReservoir {
        heat_capacity_j_kg_k: 50.0,
        relaxation_time_s: 2e-5,
    }];
    let frozen_cp = gas.cp_j_kg_k - 50.0;
    let duration = 2e-5;
    let build = |duration| {
        model
            .prepare_time_with_volume_sources(
                duration,
                frozen_cp,
                &reservoirs,
                None,
                None,
                &[&force],
                &[&heat],
                ThermoviscousTimeLimits::default(),
                &mut || false,
            )
            .unwrap()
    };
    let full = build(duration);
    let half = build(duration * 0.5);
    let pressure = [0.2, -0.1];
    let force_drive = [3e-6];
    let heat_drive = [3e-4];
    let zero = full.zero_state();
    let whole = full
        .step(
            &zero,
            &pressure,
            &[],
            &force_drive,
            &heat_drive,
            &mut || false,
        )
        .unwrap();
    let first = half
        .step(
            &zero,
            &pressure,
            &[],
            &force_drive,
            &heat_drive,
            &mut || false,
        )
        .unwrap();
    let second = half
        .step(
            &first.field.state,
            &pressure,
            &[],
            &force_drive,
            &heat_drive,
            &mut || false,
        )
        .unwrap();
    same_state(&whole.field.state, &second.field.state);
    assert_eq!(full.caloric_reservoirs().len(), 1);
    assert_eq!(whole.field.state.molecular_temperature_k.len(), 1);
    assert!(
        whole.field.state.molecular_temperature_k[0]
            .iter()
            .any(|v| *v != 0.0)
    );
    for (actual, a, b) in [
        (
            whole.field.work.body_force_work_j,
            first.field.work.body_force_work_j,
            second.field.work.body_force_work_j,
        ),
        (
            whole.field.work.heat_source_work_j,
            first.field.work.heat_source_work_j,
            second.field.work.heat_source_work_j,
        ),
        (
            whole.field.work.fluid_viscous_heat_j,
            first.field.work.fluid_viscous_heat_j,
            second.field.work.fluid_viscous_heat_j,
        ),
        (
            whole.field.work.fluid_conductive_heat_j,
            first.field.work.fluid_conductive_heat_j,
            second.field.work.fluid_conductive_heat_j,
        ),
        (
            whole.field.work.molecular_relaxation_heat_j,
            first.field.work.molecular_relaxation_heat_j,
            second.field.work.molecular_relaxation_heat_j,
        ),
    ] {
        close(actual, a + b);
    }
    close(
        whole.mean_body_force_velocity_m_s[0],
        0.5 * (first.mean_body_force_velocity_m_s[0] + second.mean_body_force_velocity_m_s[0]),
    );
    close(
        whole.mean_heat_source_temperature_ratio[0],
        0.5 * (first.mean_heat_source_temperature_ratio[0]
            + second.mean_heat_source_temperature_ratio[0]),
    );
    close(
        whole.body_force_work_by_source_j[0],
        whole.field.work.body_force_work_j,
    );
    close(
        whole.heat_source_work_by_source_j[0],
        whole.field.work.heat_source_work_j,
    );

    // With no force columns the original heat path is bitwise unchanged.
    let legacy = model
        .prepare_time_with_heat_sources(
            duration,
            frozen_cp,
            &reservoirs,
            None,
            None,
            &[&heat],
            ThermoviscousTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let empty_force = model
        .prepare_time_with_volume_sources(
            duration,
            frozen_cp,
            &reservoirs,
            None,
            None,
            &[],
            &[&heat],
            ThermoviscousTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let a = legacy
        .step(
            &legacy.zero_state(),
            &pressure,
            &[],
            &heat_drive,
            &mut || false,
        )
        .unwrap();
    let b = empty_force
        .step(
            &empty_force.zero_state(),
            &pressure,
            &[],
            &[],
            &heat_drive,
            &mut || false,
        )
        .unwrap();
    assert_eq!(a.field.state, b.field.state);
    assert_eq!(a.field.mean_inward_flow_m3_s, b.field.mean_inward_flow_m3_s);
    assert_eq!(
        a.mean_heat_source_temperature_ratio,
        b.mean_heat_source_temperature_ratio
    );
    assert_eq!(
        a.heat_source_work_by_source_j,
        b.heat_source_work_by_source_j
    );
    assert_eq!(a.field.work.supplied_j(), b.field.work.supplied_j());
    assert_eq!(b.field.work.body_force_work_j, 0.0);
    assert!(b.body_force_work_by_source_j.is_empty());
}

#[test]
fn g0_volume_time_preserves_original_identity_complete_counts_and_transactional_refusals() {
    let size = [0.004, 0.002, 0.002];
    let (mesh, positions) = box_mesh([1, 1, 1], size);
    let boundary = boundaries(&mesh, &positions, size, false);
    let make = || {
        ThermoviscousFe::new(&mesh, &positions, air(), &boundary, limits(), &mut || false).unwrap()
    };
    let model = make();
    let other = make();
    let mut law = |_, _| Ok::<_, ThermoviscousError>([1.0, -2.0, 3.0]);
    let source = model
        .volume_force_source(4, &mut law, &mut || false)
        .unwrap();
    let foreign = other
        .volume_force_source(4, &mut law, &mut || false)
        .unwrap();
    for sources in [&[][..], &[&foreign][..]] {
        assert!(
            model
                .prepare_time_with_volume_sources(
                    1e-5,
                    air().cp_j_kg_k,
                    &[],
                    None,
                    None,
                    sources,
                    &[],
                    ThermoviscousTimeLimits::default(),
                    &mut || false,
                )
                .is_err()
        );
    }
    for count in [MAX_BODY_FORCE_SOURCES, MAX_BODY_FORCE_SOURCES + 1] {
        // The first is within the body-source count cap but exceeds the original
        // complete exponential port/block budget; the second exceeds both.
        let sources = vec![&source; count];
        assert!(
            model
                .prepare_time_with_volume_sources(
                    1e-5,
                    air().cp_j_kg_k,
                    &[],
                    None,
                    None,
                    &sources,
                    &[],
                    ThermoviscousTimeLimits::default(),
                    &mut || false,
                )
                .is_err()
        );
    }
    assert!(matches!(
        model.prepare_time_with_volume_sources(
            1e-5,
            air().cp_j_kg_k,
            &[],
            None,
            None,
            &[&source],
            &[],
            ThermoviscousTimeLimits::default(),
            &mut || true,
        ),
        Err(ThermoviscousTimeError::Fluid(ThermoviscousError::Cancelled))
    ));
    let prepared = prepare(&model, 1e-5, &[&source], &[]);
    let state = prepared.zero_state();
    let original = state.clone();
    for drive in [&[][..], &[f64::NAN][..], &[1.0, 2.0][..]] {
        assert!(
            prepared
                .step(&state, &[0.0; 2], &[], drive, &[], &mut || false)
                .is_err()
        );
    }
    assert!(
        prepared
            .step(&state, &[0.0; 2], &[], &[1.0], &[0.0], &mut || false)
            .is_err()
    );
    let mut invalid = state.clone();
    let fixed = model.velocity_map.iter().position(Option::is_none).unwrap();
    invalid.velocity_m_s[fixed][0] = 1e-10;
    assert!(
        prepared
            .step(&invalid, &[0.0; 2], &[], &[1.0], &[], &mut || false)
            .is_err()
    );
    let mut polls = 0;
    assert!(matches!(
        prepared.step(&state, &[0.0; 2], &[], &[1.0], &[], &mut || {
            polls += 1;
            polls == 2
        }),
        Err(ThermoviscousTimeError::Fluid(ThermoviscousError::Cancelled))
    ));
    assert_eq!(state, original);
}
