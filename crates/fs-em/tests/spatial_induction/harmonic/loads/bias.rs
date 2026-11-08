//! Static/AC superposition on original solved tetrahedral induction fields.
use super::*;

#[test]
fn static_bias_adds_fundamental_force_without_changing_ac_fields_or_state() {
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
    let response = state
        .solve_harmonic(
            &[[C64::new(3.2, -2.4); 3]; 6],
            &[[0.0; 3]; 8],
            3.0,
            harmonic_controls(),
            &mut || false,
        )
        .unwrap();
    let before = response.potential_wb().to_vec();
    let velocity: Vec<_> = problem
        .positions_m()
        .iter()
        .copied()
        .map(virtual_velocity)
        .collect();
    let mut previous: Option<Vec<[f64; 3]>> = None;
    for gauge in [false, true] {
        let bias = problem
            .solve(
                &[[0.0; 3]; 6],
                &trace(&problem, [0.4, -0.7, 0.9], gauge),
                static_controls(),
                &mut || false,
            )
            .unwrap();
        for selected in [&[0, 1, 2, 3, 4, 5][..], &[5, 0, 3][..]] {
            let ac = response
                .lorentz_loads(selected, LOAD_LIMITS, &mut || false)
                .unwrap();
            let dc = response
                .static_bias_lorentz_loads(&bias, selected, LOAD_LIMITS, &mut || false)
                .unwrap();
            assert_eq!(dc.omega_rad_s(), response.omega_rad_s());
            assert!(
                dc.real()
                    .nodal_forces_n()
                    .iter()
                    .flatten()
                    .any(|v| v.abs() > 1e-5)
            );
            let coefficients = [
                ac.mean(),
                ac.second_harmonic_real(),
                ac.second_harmonic_imaginary(),
                dc.real(),
                dc.imaginary(),
            ];
            let wrenches =
                coefficients.map(|c| c.wrench(&Se3::identity(), [0.0; 3]).unwrap().to_array());
            let powers = coefficients.map(|c| c.power_w(&velocity, &mut || false).unwrap());
            for phase in [0.0, 0.31, 1.2, 2.7] {
                let (nodal, wrench, power) = instantaneous_superposition_with_bias(
                    &[(&response, phase)],
                    selected,
                    Some(bias.flux_density_t()),
                );
                let weights = [
                    1.0,
                    det::cos(2.0 * phase),
                    det::sin(2.0 * phase),
                    det::cos(phase),
                    det::sin(phase),
                ];
                for (node, expected) in nodal.iter().enumerate() {
                    for axis in 0..3 {
                        near(
                            coefficients
                                .iter()
                                .zip(weights)
                                .map(|(c, w)| w * c.nodal_forces_n()[node][axis])
                                .sum(),
                            expected[axis],
                            3e-9,
                        );
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
            if selected.len() == 6 {
                if let Some(previous) = &previous {
                    for (a, b) in previous
                        .iter()
                        .flatten()
                        .zip(dc.real().nodal_forces_n().iter().flatten())
                    {
                        near(*a, *b, 3e-10);
                    }
                } else {
                    previous = Some(dc.real().nodal_forces_n().to_vec());
                }
            }
            let mut canonical = selected.to_vec();
            canonical.sort_unstable();
            let reordered = response
                .static_bias_lorentz_loads(&bias, &canonical, LOAD_LIMITS, &mut || false)
                .unwrap();
            assert_eq!(
                dc.real().nodal_forces_n(),
                reordered.real().nodal_forces_n()
            );
            assert_eq!(
                dc.imaginary().nodal_forces_n(),
                reordered.imaginary().nodal_forces_n()
            );
        }
    }
    let zero = response
        .static_bias_lorentz_loads(&initial, &[0, 3, 5], LOAD_LIMITS, &mut || false)
        .unwrap();
    assert!(
        zero.real()
            .nodal_forces_n()
            .iter()
            .flatten()
            .all(|v| *v == 0.0)
    );
    assert!(
        zero.imaginary()
            .nodal_forces_n()
            .iter()
            .flatten()
            .all(|v| *v == 0.0)
    );
    assert_eq!(response.potential_wb(), before);
    assert_eq!(state.time_s(), 0.0);
    assert_eq!(state.potential_wb(), initial.potential_wb());
}

#[test]
fn static_bias_refuses_missing_dc_physics_and_incompatible_source_or_trace() {
    let (problem, sources) = cube(1, None);
    let initial = problem
        .solve(
            &[[0.0; 3]; 6],
            &trace(&problem, [0.2, 0.1, -0.3], false),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    let state =
        SpatialInduction::new(&initial, sources, preparation(), 2e-9, &mut || false).unwrap();
    let response = state
        .solve_harmonic(
            &[[C64::new(4.0, 0.0); 3]; 6],
            &[[0.0; 3]; 8],
            3.0,
            harmonic_controls(),
            &mut || false,
        )
        .unwrap();
    let moving: Vec<_> = problem
        .positions_m()
        .iter()
        .map(|p| p.map(|x| 0.35 * x))
        .collect();
    let moving = state
        .solve_harmonic(
            &[[C64::new(4.0, 0.0); 3]; 6],
            &moving,
            3.0,
            harmonic_controls(),
            &mut || false,
        )
        .unwrap();
    assert!(matches!(
        moving.static_bias_lorentz_loads(&initial, &[0], LOAD_LIMITS, &mut || false),
        Err(SpatialLorentzLoadError::Invalid(_))
    ));
    let dc_current = problem
        .solve(
            &[[4.0; 3]; 6],
            &trace(&problem, [0.0; 3], false),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    let different_trace = problem
        .solve(
            &[[0.0; 3]; 6],
            &VectorMagneticBoundary::default(),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    let (other, _) = cube(1, None);
    let other_field = other
        .solve(
            &[[0.0; 3]; 6],
            &trace(&other, [0.2, 0.1, -0.3], false),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    for rejected in [&dc_current, &different_trace, &other_field] {
        assert!(matches!(
            response.static_bias_lorentz_loads(rejected, &[0], LOAD_LIMITS, &mut || false),
            Err(SpatialLorentzLoadError::Invalid(_))
        ));
    }
    for selected in [&[][..], &[0, 0][..], &[usize::MAX][..]] {
        assert!(
            response
                .static_bias_lorentz_loads(&initial, selected, LOAD_LIMITS, &mut || false)
                .is_err()
        );
    }
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
            response
                .static_bias_lorentz_loads(&initial, &[0, 3, 5], limits, &mut || false)
                .is_err()
        );
    }
    let before = response.potential_wb().to_vec();
    for stop in [0, 5, 12, 25, 35] {
        let mut polls = 0;
        assert!(matches!(
            response.static_bias_lorentz_loads(&initial, &[0, 3, 5], LOAD_LIMITS, &mut || {
                polls += 1;
                polls > stop
            }),
            Err(SpatialLorentzLoadError::Cancelled)
        ));
    }
    assert_eq!(response.potential_wb(), before);
    assert_eq!(state.time_s(), 0.0);
}
