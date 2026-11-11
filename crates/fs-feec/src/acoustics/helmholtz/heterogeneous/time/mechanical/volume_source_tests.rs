use super::*;
use crate::acoustics::helmholtz::HelmholtzFe;
use crate::acoustics::helmholtz::heterogeneous::MAX_VOLUME_SOURCES;
use crate::acoustics::helmholtz::heterogeneous::time::tests::{MEDIUM, brick, p2, plane, uniform};
use fs_math::det;
use std::cell::Cell;

fn trace(prepared: &PreparedMechanicalTime<'_, '_, '_>, field: &[f64]) -> Vec<f64> {
    prepared
        .wall_dofs(0)
        .unwrap()
        .iter()
        .map(|&i| field[i])
        .collect()
}

fn close(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() < 5e-10 * scale,
        "actual={actual:e}, expected={expected:e}, scale={scale:e}"
    );
}

#[test]
fn distributed_injection_and_extraction_drive_complete_viscoelastic_wall_motion() {
    let pin = WallPin {
        surface_density: 0.4,
        stiffness_per_area: 1.7,
        resistance: 0.3,
    };
    let arms = vec![
        AcousticWallRelaxation {
            stiffness_per_area: 2.5,
            relaxation_time_s: 0.04,
        },
        AcousticWallRelaxation {
            stiffness_per_area: 0.8,
            relaxation_time_s: 0.2,
        },
    ];
    let drag = pin.resistance
        + arms
            .iter()
            .map(|arm| arm.stiffness_per_area * arm.relaxation_time_s)
            .sum::<f64>();
    let beta = 1.0 / (MEDIUM.density_kg_m3 * MEDIUM.sound_speed_m_s.powi(2));
    let h = 0.0125;
    for height in [0.1, 0.3] {
        let lengths = [1.0, height, 0.2];
        let area = lengths[1] * lengths[2];
        let volume = area * lengths[0];
        let (complex, positions) = brick(2, lengths);
        let geometry = HelmholtzFe::new(&complex, &positions, 2, MEDIUM).unwrap();
        let model = uniform(&geometry);
        let calls = Cell::new(0);
        let source = model
            .volume_velocity_source(
                4,
                &mut |_, _| {
                    calls.set(calls.get() + 1);
                    Ok::<_, HeterogeneousError>(1.0 / volume)
                },
                &mut || false,
            )
            .unwrap();
        let prepared_calls = calls.get();
        close(source.integrated_shape(), 1.0, 1.0);
        let owner = model
            .prepare_time_with_walls(
                h,
                &[],
                &[MechanicalBoundary {
                    faces: plane(&geometry, 1.0),
                    wall: pin,
                    relaxation: arms.clone(),
                }],
                MechanicalTimeLimits::default(),
                &mut || false,
            )
            .unwrap();
        for amplitude in [1.0, 1e-8] {
            for fraction in [-2.5, 1.5] {
                // Exact continuum polynomial, with all original wall memories:
                // w=w0+v*t, z_j=sqrt(k_j)*(w-tau_j*v), p=drag*v+K*w.
                // beta*p_dot+div(v_fluid)=q/V, v_fluid=-(a*x+b)/rho.
                // The inlet and volume source divide the physical input work.
                let v = 0.02 * amplitude;
                let w0 = 0.03 * amplitude;
                let q = fraction * beta * pin.stiffness_per_area * v * volume;
                let a = MEDIUM.density_kg_m3 * (beta * pin.stiffness_per_area * v - q / volume);
                let b = -MEDIUM.density_kg_m3 * v - a;
                let p0 = drag * v + pin.stiffness_per_area * w0;
                let potential = p2(&geometry, |x| 0.5 * a * x[0].powi(2) + b * x[0]);
                let pressure = p2(&geometry, |_| p0);
                let displacement = trace(&owner, &p2(&geometry, |_| w0));
                let velocity = trace(&owner, &p2(&geometry, |_| v));
                let memory = arms
                    .iter()
                    .map(|arm| {
                        trace(
                            &owner,
                            &p2(&geometry, |_| {
                                det::sqrt(arm.stiffness_per_area) * (w0 - arm.relaxation_time_s * v)
                            }),
                        )
                    })
                    .collect::<Vec<_>>();
                let mut state = owner
                    .state_with_wall_memory(
                        &potential,
                        &pressure,
                        &[displacement],
                        &[velocity],
                        &[memory],
                    )
                    .unwrap();
                let initial = owner.energy(&state).unwrap().total_j();
                let inlet_velocity = -b / MEDIUM.density_kg_m3;
                let inlet = plane(&geometry, 0.0)
                    .into_iter()
                    .map(|face| (face, [inlet_velocity; 3]))
                    .collect::<Vec<_>>();
                let mut total_input = 0.0;
                let mut total_heat = 0.0;
                for index in 1..=6 {
                    let next = owner
                        .step_with_volume_sources(&state, &inlet, &[(&source, q)], &mut || false)
                        .unwrap();
                    let time = f64::from(index) * h;
                    let w = w0 + v * time;
                    let pressure = drag * v + pin.stiffness_per_area * w;
                    let midpoint_pressure = pressure - 0.5 * h * pin.stiffness_per_area * v;
                    let expected_p = p2(&geometry, |_| pressure);
                    let expected_psi = p2(&geometry, |x| {
                        0.5 * a * x[0].powi(2)
                            + b * x[0]
                            + p0 * time
                            + 0.5 * pin.stiffness_per_area * v * time.powi(2)
                    });
                    for ((&p, &expected_p), (&psi, &expected_psi)) in
                        next.state.pressure_pa().iter().zip(&expected_p).zip(
                            next.state
                                .pressure_potential_pa_s()
                                .iter()
                                .zip(&expected_psi),
                        )
                    {
                        close(p, expected_p, pressure);
                        close(psi, expected_psi, amplitude);
                    }
                    let energy = owner.energy(&next.state).unwrap();
                    for (actual, expected) in [
                        (
                            energy.acoustic.kinetic_j,
                            area * (a * a / 3.0 + a * b + b * b) / (2.0 * MEDIUM.density_kg_m3),
                        ),
                        (
                            energy.acoustic.compression_j,
                            0.5 * beta * volume * pressure.powi(2),
                        ),
                        (
                            energy.wall_kinetic_j,
                            0.5 * pin.surface_density * area * v * v,
                        ),
                        (
                            energy.wall_elastic_j,
                            0.5 * pin.stiffness_per_area * area * w * w,
                        ),
                        (
                            energy.wall_relaxation_j,
                            0.5 * area
                                * v
                                * v
                                * arms
                                    .iter()
                                    .map(|arm| {
                                        arm.stiffness_per_area * arm.relaxation_time_s.powi(2)
                                    })
                                    .sum::<f64>(),
                        ),
                    ] {
                        close(actual, expected, expected);
                    }
                    let volume_work = h * midpoint_pressure * q;
                    close(
                        next.work.volume_source_work_j,
                        volume_work,
                        volume_work.abs(),
                    );
                    assert_eq!(next.work.volume_source_work_j.is_sign_negative(), q < 0.0);
                    close(
                        next.work.input_work_j,
                        h * midpoint_pressure * (area * inlet_velocity + q),
                        initial,
                    );
                    close(
                        next.work.wall_work_j,
                        h * midpoint_pressure * area * v,
                        initial,
                    );
                    close(
                        next.work.wall_dissipation_j,
                        h * drag * area * v * v,
                        initial,
                    );
                    close(
                        next.work.wall_relaxation_dissipation_j,
                        h * (drag - pin.resistance) * area * v * v,
                        initial,
                    );
                    assert_eq!(next.work.boundary_dissipation_j, 0.0);
                    total_input += next.work.input_work_j;
                    total_heat += next.work.wall_dissipation_j;
                    close(
                        energy.total_j() + total_heat,
                        initial + total_input,
                        initial,
                    );
                    state = next.state;
                }
            }
        }
        assert_eq!(
            calls.get(),
            prepared_calls,
            "warm steps must reuse original source quadrature"
        );
    }
}

#[test]
fn zero_integral_volume_source_moves_release_adjacent_walls_and_retains_source_work() {
    let lengths = [1.0, 0.2, 0.3];
    let volume = lengths.iter().product::<f64>();
    let (complex, positions) = brick(2, lengths);
    let geometry = HelmholtzFe::new(&complex, &positions, 2, MEDIUM).unwrap();
    let model = uniform(&geometry);
    let source = model
        .volume_velocity_source(
            4,
            &mut |_, x| Ok::<_, HeterogeneousError>((2.0 * x[0] - 1.0) / volume),
            &mut || false,
        )
        .unwrap();
    assert!(source.integrated_shape().abs() < 1e-13);
    let walls = geometry
        .boundary_faces()
        .into_iter()
        .filter(|face| face.iter().all(|&i| positions[i as usize][1] == 0.0))
        .collect();
    let owner = model
        .prepare_time_with_walls_and_pressure_release(
            0.01,
            &[],
            &[MechanicalBoundary {
                faces: walls,
                wall: WallPin {
                    surface_density: 0.6,
                    stiffness_per_area: 3.0,
                    resistance: 0.2,
                },
                relaxation: vec![AcousticWallRelaxation {
                    stiffness_per_area: 2.0,
                    relaxation_time_s: 0.04,
                }],
            }],
            &plane(&geometry, 1.0),
            MechanicalTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let mut state = owner.zero_state();
    let mut input = 0.0_f64;
    let mut heat = 0.0;
    for index in 0..12 {
        let q = if index < 4 { 0.002 } else { 0.0 };
        let next = owner
            .step_with_volume_sources(&state, &[], &[(&source, q)], &mut || false)
            .unwrap();
        close(
            next.work.input_work_j,
            next.work.volume_source_work_j,
            input.abs().max(next.work.input_work_j.abs()).max(1e-15),
        );
        input += next.work.input_work_j;
        heat += next.work.wall_dissipation_j;
        close(
            owner.energy(&next.state).unwrap().total_j() + heat,
            input,
            input,
        );
        for (dof, &fixed) in owner
            .condensed
            .pressure_release
            .constrained
            .iter()
            .enumerate()
        {
            if fixed {
                assert_eq!(next.state.pressure_pa()[dof], 0.0);
                assert_eq!(next.state.pressure_potential_pa_s()[dof], 0.0);
            }
        }
        state = next.state;
    }
    let energy = owner.energy(&state).unwrap();
    assert!(input > 0.0 && heat > 0.0);
    assert!(
        energy.wall_kinetic_j > 0.0
            && energy.wall_elastic_j > 0.0
            && energy.wall_relaxation_j > 0.0
    );
    assert!(
        state
            .wall_relaxation_memory(0, 0)
            .unwrap()
            .iter()
            .any(|&z| z != 0.0)
    );
}

#[test]
fn mechanical_volume_sources_are_bound_transactional_and_empty_path_is_identical() {
    let (complex, positions) = brick(1, [1.0, 0.2, 0.3]);
    let geometry = HelmholtzFe::new(&complex, &positions, 2, MEDIUM).unwrap();
    let model = uniform(&geometry);
    let foreign = uniform(&geometry);
    let source = model
        .volume_velocity_source(4, &mut |_, _| Ok::<_, HeterogeneousError>(1.0), &mut || {
            false
        })
        .unwrap();
    let foreign_source = foreign
        .volume_velocity_source(4, &mut |_, _| Ok::<_, HeterogeneousError>(1.0), &mut || {
            false
        })
        .unwrap();
    let owner = model
        .prepare_time_with_walls(
            0.01,
            &[],
            &[MechanicalBoundary {
                faces: plane(&geometry, 1.0),
                wall: WallPin {
                    surface_density: 0.6,
                    stiffness_per_area: 3.0,
                    resistance: 0.2,
                },
                relaxation: Vec::new(),
            }],
            MechanicalTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let initial = owner
        .state(&vec![0.0; model.ndof()], &p2(&geometry, |_| 0.1))
        .unwrap();
    for invalid in [&[(&foreign_source, 0.1)][..], &[(&source, f64::NAN)][..]] {
        assert!(
            owner
                .step_with_volume_sources(&initial, &[], invalid, &mut || false)
                .is_err()
        );
    }
    assert!(matches!(
        owner.step_with_volume_sources(
            &initial,
            &[],
            &vec![(&source, 0.0); MAX_VOLUME_SOURCES + 1],
            &mut || false
        ),
        Err(HeterogeneousError::Limit { .. })
    ));
    assert!(matches!(
        owner.step_with_volume_sources(&initial, &[], &[(&source, 0.1)], &mut || true),
        Err(HeterogeneousError::Cancelled)
    ));
    let old = owner.step(&initial, &[], &mut || false).unwrap();
    let empty = owner
        .step_with_volume_sources(&initial, &[], &[], &mut || false)
        .unwrap();
    let cancelled = owner
        .step_with_volume_sources(
            &initial,
            &[],
            &[(&source, 0.125), (&source, -0.125)],
            &mut || false,
        )
        .unwrap();
    for other in [&empty, &cancelled] {
        for (a, b) in old
            .state
            .pressure_pa()
            .iter()
            .zip(other.state.pressure_pa())
        {
            assert_eq!(a.to_bits(), b.to_bits());
        }
        for (a, b) in old
            .state
            .pressure_potential_pa_s()
            .iter()
            .zip(other.state.pressure_potential_pa_s())
        {
            assert_eq!(a.to_bits(), b.to_bits());
        }
        assert_eq!(
            old.state.wall_displacement_m(0),
            other.state.wall_displacement_m(0)
        );
        assert_eq!(
            old.state.wall_normal_velocity_m_s(0),
            other.state.wall_normal_velocity_m_s(0)
        );
        assert_eq!(
            old.work.input_work_j.to_bits(),
            other.work.input_work_j.to_bits()
        );
        assert_eq!(other.work.volume_source_work_j, 0.0);
    }
    assert_eq!(initial.pressure_pa(), p2(&geometry, |_| 0.1));
}
