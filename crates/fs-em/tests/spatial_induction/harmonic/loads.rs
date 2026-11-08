//! AC mechanical forcing from the original solved field, including magnetic hum.
use super::*;
use fs_em::spatial::vector_potential::induction::{
    harmonic::SpatialHarmonicInductionResponse,
    mechanical::{SpatialLorentzLoadError, SpatialLorentzLoadLimits},
};
use fs_ga::{Se3, Twist, Vec3};
use fs_math::det;

const LOAD_LIMITS: SpatialLorentzLoadLimits = SpatialLorentzLoadLimits {
    max_vertices: 8,
    max_cells: 6,
};

fn virtual_velocity(p: [f64; 3]) -> [f64; 3] {
    [
        0.3 + 2.0 * p[0] - 0.2 * p[2],
        -0.4 + 0.7 * p[1] + 0.1 * p[2],
        0.8 - 0.5 * p[0] + 1.3 * p[2],
    ]
}

fn instantaneous_reference(
    response: &SpatialHarmonicInductionResponse<'_, '_>,
    selected: &[usize],
    phase: f64,
) -> (Vec<[f64; 3]>, [f64; 6], f64) {
    instantaneous_superposition_reference(&[(response, phase)], selected)
}

fn instantaneous_superposition_reference(
    responses: &[(&SpatialHarmonicInductionResponse<'_, '_>, f64)],
    selected: &[usize],
) -> (Vec<[f64; 3]>, [f64; 6], f64) {
    let problem = responses[0].0.problem();
    let mut nodal = vec![[0.0; 3]; problem.mesh().vertex_count];
    let mut wrench = [0.0; 6];
    let mut power = 0.0;
    // Independent four-point degree-two tetrahedral quadrature. Reconstruct
    // real J(t), B(t) BEFORE the cross product; do not use Fourier force formulas
    // or the production load integrator's mass weights.
    let high = (5.0 + 3.0 * 5.0_f64.sqrt()) / 20.0;
    let low = (5.0 - 5.0_f64.sqrt()) / 20.0;
    for &cell in selected {
        let record = responses[0].0.cells().iter().find(|c| c.cell == cell).unwrap();
        let tet = problem.mesh().tets[cell].map(|v| v as usize);
        let mut b = [0.0; 3];
        let mut current = [[0.0; 3]; 4];
        for &(response, phase) in responses {
            let (sine, cosine) = det::sin_cos(phase);
            let at_phase = |z: C64| z.re * cosine + z.im * sine;
            let record = response.cells().iter().find(|c| c.cell == cell).unwrap();
            for axis in 0..3 {
                b[axis] += at_phase(response.flux_density_t()[cell][axis]);
                for vertex in 0..4 {
                    current[vertex][axis] += at_phase(record.current_density_a_m2[vertex][axis]);
                }
            }
        }
        for point in 0..4 {
            let shape: [f64; 4] = std::array::from_fn(|i| if i == point { high } else { low });
            let j = std::array::from_fn(|axis| {
                (0..4)
                    .map(|i| shape[i] * current[i][axis])
                    .sum()
            });
            let position = std::array::from_fn(|axis| {
                (0..4)
                    .map(|i| shape[i] * problem.positions_m()[tet[i]][axis])
                    .sum()
            });
            let force = cross(j, b);
            let moment = cross(position, force);
            let velocity = virtual_velocity(position);
            let weight = record.volume_m3 / 4.0;
            for axis in 0..3 {
                wrench[axis] += weight * moment[axis];
                wrench[axis + 3] += weight * force[axis];
                power += weight * force[axis] * velocity[axis];
                for i in 0..4 {
                    nodal[tet[i]][axis] += weight * shape[i] * force[axis];
                }
            }
        }
    }
    (nodal, wrench, power)
}

#[test]
fn mean_and_twice_frequency_loads_reconstruct_instantaneous_force_torque_and_work() {
    let (problem, sources) = cube(1, None);
    let initial = problem
        .solve(
            &[[0.0; 3]; 6],
            &trace(&problem, [0.0; 3], false),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    let state =
        SpatialInduction::new(&initial, sources, preparation(), 2e-9, &mut || false).unwrap();
    let virtual_velocity: Vec<_> = problem
        .positions_m()
        .iter()
        .copied()
        .map(virtual_velocity)
        .collect();
    for alpha in [0.0, 0.35] {
        let velocity: Vec<_> = problem
            .positions_m()
            .iter()
            .map(|p| p.map(|x| alpha * x))
            .collect();
        let response = state
            .solve_harmonic(
                &[[C64::new(3.2, -2.4); 3]; 6],
                &velocity,
                3.0,
                harmonic_controls(),
                &mut || false,
            )
            .unwrap();
        for selected in [&[0, 1, 2, 3, 4, 5][..], &[5, 0, 3][..]] {
            let loads = response
                .lorentz_loads(selected, LOAD_LIMITS, &mut || false)
                .unwrap();
            near(loads.second_harmonic_omega_rad_s(), 6.0, 1e-15);
            let expected_mean_power = response
                .cells()
                .iter()
                .filter(|c| selected.contains(&c.cell))
                .map(|c| c.mean_mechanical_power_w)
                .sum();
            near(
                loads.mean().power_w(&velocity, &mut || false).unwrap(),
                expected_mean_power,
                2e-10,
            );
            let coefficients = [
                loads.mean(),
                loads.second_harmonic_real(),
                loads.second_harmonic_imaginary(),
            ];
            let wrenches =
                coefficients.map(|c| c.wrench(&Se3::identity(), [0.0; 3]).unwrap().to_array());
            let powers = coefficients.map(|c| c.power_w(&virtual_velocity, &mut || false).unwrap());
            assert!(
                coefficients[1..]
                    .iter()
                    .flat_map(|c| c.nodal_forces_n())
                    .flatten()
                    .any(|v| v.abs() > 1e-5)
            );
            if alpha == 0.0 {
                // A standing-wave mode has zero mean force yet nonzero 2omega
                // forcing. A mean-only export would entirely lose this load.
                assert!(
                    loads
                        .mean()
                        .nodal_forces_n()
                        .iter()
                        .flatten()
                        .all(|v| v.abs() < 2e-9)
                );
            }
            for phase in [0.0, 0.31, 1.2, 2.7] {
                let (reference, wrench, power) =
                    instantaneous_reference(&response, selected, phase);
                let (sine, cosine) = det::sin_cos(2.0 * phase);
                for (node, expected) in reference.iter().enumerate() {
                    for axis in 0..3 {
                        near(
                            coefficients[0].nodal_forces_n()[node][axis]
                                + cosine * coefficients[1].nodal_forces_n()[node][axis]
                                + sine * coefficients[2].nodal_forces_n()[node][axis],
                            expected[axis],
                            2e-10,
                        );
                    }
                }
                for axis in 0..6 {
                    near(
                        wrenches[0][axis] + cosine * wrenches[1][axis] + sine * wrenches[2][axis],
                        wrench[axis],
                        2e-10,
                    );
                }
                near(
                    powers[0] + cosine * powers[1] + sine * powers[2],
                    power,
                    2e-10,
                );
            }
        }
    }
}

#[test]
fn harmonic_loads_preserve_squared_source_phase_body_frames_and_atomic_selection() {
    let (problem, sources) = cube(1, None);
    let initial = problem
        .solve(
            &[[0.0; 3]; 6],
            &trace(&problem, [0.0; 3], false),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    let state =
        SpatialInduction::new(&initial, sources, preparation(), 2e-9, &mut || false).unwrap();
    let velocity: Vec<_> = problem
        .positions_m()
        .iter()
        .map(|p| p.map(|x| 0.35 * x))
        .collect();
    let response = state
        .solve_harmonic(
            &[[C64::new(4.0, 0.0); 3]; 6],
            &velocity,
            3.0,
            harmonic_controls(),
            &mut || false,
        )
        .unwrap();
    let loads = response
        .lorentz_loads(&[5, 0, 3], LOAD_LIMITS, &mut || false)
        .unwrap();
    let phase = C64::new(0.6, 0.8).scale(1.3);
    let changed = state
        .solve_harmonic(
            &[[phase.scale(4.0); 3]; 6],
            &velocity,
            3.0,
            harmonic_controls(),
            &mut || false,
        )
        .unwrap();
    let shifted = changed
        .lorentz_loads(&[3, 5, 0], LOAD_LIMITS, &mut || false)
        .unwrap();
    for node in 0..problem.mesh().vertex_count {
        for axis in 0..3 {
            near(
                shifted.mean().nodal_forces_n()[node][axis],
                phase.norm_sq() * loads.mean().nodal_forces_n()[node][axis],
                2e-9,
            );
            let expected = C64::new(
                loads.second_harmonic_real().nodal_forces_n()[node][axis],
                loads.second_harmonic_imaginary().nodal_forces_n()[node][axis],
            ) * phase
                * phase;
            complex_near(
                C64::new(
                    shifted.second_harmonic_real().nodal_forces_n()[node][axis],
                    shifted.second_harmonic_imaginary().nodal_forces_n()[node][axis],
                ),
                expected,
                2e-9,
            );
        }
    }
    let reordered = response
        .lorentz_loads(&[0, 3, 5], LOAD_LIMITS, &mut || false)
        .unwrap();
    for (a, b) in [
        loads.mean(),
        loads.second_harmonic_real(),
        loads.second_harmonic_imaginary(),
    ]
    .into_iter()
    .zip([
        reordered.mean(),
        reordered.second_harmonic_real(),
        reordered.second_harmonic_imaginary(),
    ]) {
        assert_eq!(a.nodal_forces_n(), b.nodal_forces_n());
    }

    let placement = Se3::exp(Twist::new(
        Vec3::new(0.2, -0.4, 0.3),
        Vec3::new(2.0, -1.0, 3.0),
    ))
    .unwrap();
    let origin = Vec3::new(4.0, -2.0, 1.0);
    let twist = Twist::new(Vec3::new(0.7, -0.2, 0.4), Vec3::new(-0.5, 0.3, 1.1));
    let inverse = placement.inverse().unwrap();
    let original_velocity: Vec<_> = problem
        .positions_m()
        .iter()
        .map(|p| {
            let arm = placement
                .transform_point(Vec3::new(p[0], p[1], p[2]))
                .unwrap()
                - origin;
            let v = inverse
                .transform_vector(twist.linear + twist.angular.cross(arm))
                .unwrap();
            [v.x, v.y, v.z]
        })
        .collect();
    for coefficient in [
        loads.mean(),
        loads.second_harmonic_real(),
        loads.second_harmonic_imaginary(),
    ] {
        let wrench = coefficient
            .wrench(&placement, [origin.x, origin.y, origin.z])
            .unwrap();
        near(
            wrench.pairing(twist),
            coefficient
                .power_w(&original_velocity, &mut || false)
                .unwrap(),
            2e-10,
        );
    }
    for selection in [&[][..], &[0, 0][..], &[usize::MAX][..]] {
        assert!(matches!(
            response.lorentz_loads(selection, LOAD_LIMITS, &mut || false),
            Err(SpatialLorentzLoadError::Invalid(_))
        ));
    }
    let all = [0, 1, 2, 3, 4, 5];
    for limits in [
        SpatialLorentzLoadLimits {
            max_vertices: 7,
            ..LOAD_LIMITS
        },
        SpatialLorentzLoadLimits {
            max_cells: 5,
            ..LOAD_LIMITS
        },
    ] {
        assert!(response.lorentz_loads(&all, limits, &mut || false).is_err());
    }
    let before = response.potential_wb().to_vec();
    for stop in [0, 12, 27, 42] {
        let mut polls = 0;
        assert!(matches!(
            response.lorentz_loads(&all, LOAD_LIMITS, &mut || {
                polls += 1;
                polls > stop
            }),
            Err(SpatialLorentzLoadError::Cancelled)
        ));
    }
    assert_eq!(response.potential_wb(), before);
    assert_eq!(state.time_s(), 0.0);
}

#[test]
fn cross_carrier_sum_and_beat_forces_reconstruct_original_fields() {
    let (problem, sources) = cube(1, None);
    let initial = problem
        .solve(
            &[[0.0; 3]; 6],
            &trace(&problem, [0.0; 3], false),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    let state =
        SpatialInduction::new(&initial, sources, preparation(), 2e-9, &mut || false).unwrap();
    let velocity: Vec<_> = problem
        .positions_m()
        .iter()
        .map(|p| p.map(|x| 0.35 * x))
        .collect();
    let virtual_velocity: Vec<_> = problem
        .positions_m()
        .iter()
        .copied()
        .map(virtual_velocity)
        .collect();
    let first = state
        .solve_harmonic(
            &[[C64::new(3.2, -2.4); 3]; 6],
            &velocity,
            3.0,
            harmonic_controls(),
            &mut || false,
        )
        .unwrap();
    // A true third harmonic, an inharmonic beat, and a coherent equal-frequency
    // source all use the same original spatial operator and retain their phase.
    for omega in [9.0, 3.125, 3.0] {
        let second = state
            .solve_harmonic(
                &[[C64::new(-1.8, 0.7); 3]; 6],
                &velocity,
                omega,
                harmonic_controls(),
                &mut || false,
            )
            .unwrap();
        for selected in [&[0, 1, 2, 3, 4, 5][..], &[5, 0, 3][..]] {
            let a = first
                .lorentz_loads(selected, LOAD_LIMITS, &mut || false)
                .unwrap();
            let b = second
                .lorentz_loads(selected, LOAD_LIMITS, &mut || false)
                .unwrap();
            let mixed = first
                .cross_lorentz_loads(&second, selected, LOAD_LIMITS, &mut || false)
                .unwrap();
            near(mixed.difference_omega_rad_s(), (omega - 3.0).abs(), 1e-15);
            near(mixed.sum_omega_rad_s(), omega + 3.0, 1e-15);
            let swapped = second
                .cross_lorentz_loads(&first, selected, LOAD_LIMITS, &mut || false)
                .unwrap();
            let cross_coefficients = [
                mixed.difference_real(),
                mixed.difference_imaginary(),
                mixed.sum_real(),
                mixed.sum_imaginary(),
            ];
            for (left, right) in cross_coefficients.into_iter().zip([
                swapped.difference_real(),
                swapped.difference_imaginary(),
                swapped.sum_real(),
                swapped.sum_imaginary(),
            ]) {
                assert_eq!(left.nodal_forces_n(), right.nodal_forces_n());
            }
            assert!(
                cross_coefficients
                    .iter()
                    .flat_map(|c| c.nodal_forces_n())
                    .flatten()
                    .any(|value| value.abs() > 1e-5)
            );
            if omega == 3.0 {
                assert!(
                    mixed
                        .difference_imaginary()
                        .nodal_forces_n()
                        .iter()
                        .flatten()
                        .all(|value| *value == 0.0)
                );
            }
            let coefficients = [
                a.mean(),
                a.second_harmonic_real(),
                a.second_harmonic_imaginary(),
                b.mean(),
                b.second_harmonic_real(),
                b.second_harmonic_imaginary(),
                mixed.difference_real(),
                mixed.difference_imaginary(),
                mixed.sum_real(),
                mixed.sum_imaginary(),
            ];
            let wrenches =
                coefficients.map(|c| c.wrench(&Se3::identity(), [0.0; 3]).unwrap().to_array());
            let powers = coefficients.map(|c| c.power_w(&virtual_velocity, &mut || false).unwrap());
            for time in [0.0, 0.137, 0.91, 3.7, 13.0] {
                let phase_a = 3.0 * time;
                let phase_b = omega * time;
                let (nodal, wrench, power) = instantaneous_superposition_reference(
                    &[(&first, phase_a), (&second, phase_b)],
                    selected,
                );
                let weights = [
                    1.0,
                    det::cos(2.0 * phase_a),
                    det::sin(2.0 * phase_a),
                    1.0,
                    det::cos(2.0 * phase_b),
                    det::sin(2.0 * phase_b),
                    det::cos(mixed.difference_omega_rad_s() * time),
                    det::sin(mixed.difference_omega_rad_s() * time),
                    det::cos(mixed.sum_omega_rad_s() * time),
                    det::sin(mixed.sum_omega_rad_s() * time),
                ];
                for (node, expected) in nodal.iter().enumerate() {
                    for axis in 0..3 {
                        let actual = coefficients
                            .iter()
                            .zip(weights)
                            .map(|(c, w)| w * c.nodal_forces_n()[node][axis])
                            .sum();
                        near(actual, expected[axis], 3e-9);
                    }
                }
                for axis in 0..6 {
                    near(
                        wrenches.iter().zip(weights).map(|(c, w)| w * c[axis]).sum(),
                        wrench[axis],
                        3e-9,
                    );
                }
                near(
                    powers.iter().zip(weights).map(|(p, w)| p * w).sum(),
                    power,
                    3e-9,
                );
            }
        }
    }
    assert_eq!(state.time_s(), 0.0);
}

#[test]
fn cross_carrier_loads_refuse_incompatible_sources_motion_and_partial_work() {
    let (problem, sources) = cube(1, None);
    let initial = problem
        .solve(
            &[[0.0; 3]; 6],
            &trace(&problem, [0.0; 3], false),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    let other_state = SpatialInduction::new(
        &initial,
        sources.clone(),
        preparation(),
        2e-9,
        &mut || false,
    )
    .unwrap();
    let state =
        SpatialInduction::new(&initial, sources, preparation(), 2e-9, &mut || false).unwrap();
    let stationary = [[0.0; 3]; 8];
    let moving: Vec<_> = problem
        .positions_m()
        .iter()
        .map(|p| p.map(|x| 0.35 * x))
        .collect();
    let source = [[C64::new(4.0, 0.0); 3]; 6];
    let first = state
        .solve_harmonic(
            &source,
            &stationary,
            3.0,
            harmonic_controls(),
            &mut || false,
        )
        .unwrap();
    let second = state
        .solve_harmonic(
            &source,
            &stationary,
            9.0,
            harmonic_controls(),
            &mut || false,
        )
        .unwrap();
    let changed_motion = state
        .solve_harmonic(&source, &moving, 9.0, harmonic_controls(), &mut || false)
        .unwrap();
    let changed_source = other_state
        .solve_harmonic(
            &source,
            &stationary,
            9.0,
            harmonic_controls(),
            &mut || false,
        )
        .unwrap();
    for incompatible in [&changed_motion, &changed_source] {
        assert!(matches!(
            first.cross_lorentz_loads(incompatible, &[0], LOAD_LIMITS, &mut || false),
            Err(SpatialLorentzLoadError::Invalid(_))
        ));
    }
    for selection in [&[][..], &[0, 0][..], &[usize::MAX][..]] {
        assert!(matches!(
            first.cross_lorentz_loads(&second, selection, LOAD_LIMITS, &mut || false),
            Err(SpatialLorentzLoadError::Invalid(_))
        ));
    }
    let selected = [5, 0, 3];
    for limits in [
        SpatialLorentzLoadLimits {
            max_vertices: 7,
            ..LOAD_LIMITS
        },
        SpatialLorentzLoadLimits {
            max_cells: 2,
            ..LOAD_LIMITS
        },
    ] {
        assert!(
            first
                .cross_lorentz_loads(&second, &selected, limits, &mut || false)
                .is_err()
        );
    }
    let before = first.potential_wb().to_vec();
    for stop in [0, 5, 14, 30, 50] {
        let mut polls = 0;
        assert!(matches!(
            first.cross_lorentz_loads(&second, &selected, LOAD_LIMITS, &mut || {
                polls += 1;
                polls > stop
            }),
            Err(SpatialLorentzLoadError::Cancelled)
        ));
    }
    assert_eq!(first.potential_wb(), before);
    assert_eq!(state.time_s(), 0.0);
}
