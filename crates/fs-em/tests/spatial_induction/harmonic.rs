//! G1/G3: independently integrated cube mode and time-domain AC agreement.
use super::*;
use fs_em::spatial::vector_potential::induction::harmonic::SpatialHarmonicInductionControls;
use fs_math::c64::C64;

fn harmonic_controls() -> SpatialHarmonicInductionControls {
    SpatialHarmonicInductionControls {
        relative_tolerance: 2e-11,
        absolute_residual_a: 1e-11,
        restart: 64,
        maximum_cycles: 100,
        maximum_krylov_scalars: 2_000_000,
        maximum_field_scalars: 2_000_000,
        current_tolerance_a: 2e-8,
        flux_tolerance_wb: 2e-9,
        power_tolerance_w: 2e-8,
        maximum_cell_magnetic_peclet: 100.0,
        maximum_cell_skin_ratio: 1000.0,
    }
}

fn complex_near(a: C64, b: C64, tolerance: f64) {
    near(a.re, b.re, tolerance);
    near(a.im, b.im, tolerance);
}

#[test]
fn harmonic_cube_matches_exact_diffusion_phase_charge_heat_and_reactive_storage() {
    // Independently integrated constants: K=4/mu, M=3*sigma/40 after
    // eliminating nodal charge. This fixture has one magnetic unknown.
    let (problem, sources) = cube(1, None);
    let initial = problem
        .solve(
            &vec![[0.0; 3]; 6],
            &trace(&problem, [0.0; 3], false),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    let state =
        SpatialInduction::new(&initial, sources, preparation(), 2e-9, &mut || false).unwrap();
    let (m, k) = (3.0 * 5.0 / 40.0, 4.0 / 2.0);
    let phase = C64::new(0.8, -0.6);
    for omega in [0.03, 3.0, 40.0, 400.0] {
        let response = state
            .solve_harmonic(
                &vec![[phase.scale(4.0); 3]; 6],
                &[[0.0; 3]; 8],
                omega,
                harmonic_controls(),
                &mut || false,
            )
            .unwrap();
        let a = phase.scale(k) / C64::new(k, -omega * m);
        complex_near(response.potential_wb()[diagonal(&problem)], a, 2e-9);
        let square = a.re * a.re + a.im * a.im;
        let balance = response.balance();
        near(balance.mean_magnetic_energy_j, 0.25 * k * square, 2e-9);
        near(
            balance.mean_ohmic_power_w,
            0.5 * m * omega * omega * square,
            2e-9,
        );
        near(
            balance.mean_supplied_power_w,
            balance.mean_ohmic_power_w,
            2e-9,
        );
        near(balance.mean_mechanical_power_w, 0.0, 1e-15);
        near(
            balance.reactive_supplied_power_var,
            -0.5 * omega * k * square,
            2e-9,
        );
        near(balance.power_residual_w, 0.0, 2e-9);
        near(balance.reactive_power_residual_var, 0.0, 2e-9);
        near(
            response.maximum_cell_skin_ratio(),
            (3.0_f64 * omega * 2.0 * 5.0 / 2.0).sqrt(),
            1e-12,
        );
        assert!(response.report().converged_euclidean());
        assert!(response.maximum_ampere_residual_a() < 2e-8);
        assert!(response.maximum_charge_residual_a() < 2e-8);
        let q = |v| {
            if v == 0 {
                0.0
            } else if v == 7 {
                0.5
            } else {
                0.25
            }
        };
        for cell in response.cells() {
            let tet = problem.mesh().tets[cell.cell];
            for vertex in 0..4 {
                complex_near(
                    cell.potential_v[vertex] - cell.potential_v[0],
                    (C64::new(0.0, omega) * a).scale(q(tet[vertex]) - q(tet[0])),
                    2e-9,
                );
                for axis in 0..3 {
                    complex_near(
                        cell.current_density_a_m2[vertex][axis],
                        cell.conductor_electric_field_v_m[vertex][axis].scale(5.0),
                        2e-9,
                    );
                }
            }
            assert!(cell.mean_ohmic_power_w >= 0.0);
        }
        near(
            response.cells().iter().map(|c| c.mean_ohmic_power_w).sum(),
            balance.mean_ohmic_power_w,
            2e-9,
        );
        assert_eq!(response.conductors(), state.conductors());
        assert_eq!(state.time_s(), 0.0);
        assert_eq!(state.potential_wb(), initial.potential_wb());
    }
}

#[test]
fn moving_harmonic_mode_keeps_mechanical_power_out_of_ohmic_heat() {
    let (problem, sources) = cube(1, None);
    let initial = problem
        .solve(
            &vec![[0.0; 3]; 6],
            &trace(&problem, [0.0; 3], false),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    let state =
        SpatialInduction::new(&initial, sources, preparation(), 2e-9, &mut || false).unwrap();
    let (m, k, omega) = (3.0 * 5.0 / 40.0, 4.0 / 2.0, 3.0);
    for alpha in [0.6, -0.3] {
        let velocity: Vec<_> = problem
            .positions_m()
            .iter()
            .map(|p| p.map(|x| alpha * x))
            .collect();
        let response = state
            .solve_harmonic(
                &[[C64::new(4.0, 0.0); 3]; 6],
                &velocity,
                omega,
                harmonic_controls(),
                &mut || false,
            )
            .unwrap();
        // Projected v cross curl(W) = -2 alpha (W - grad q).
        let a = C64::new(k, 0.0) / C64::new(k + 2.0 * alpha * m, -omega * m);
        let square = a.re * a.re + a.im * a.im;
        complex_near(response.potential_wb()[diagonal(&problem)], a, 2e-9);
        let b = response.balance();
        near(
            b.mean_ohmic_power_w,
            0.5 * m * (omega * omega + 4.0 * alpha * alpha) * square,
            2e-9,
        );
        near(
            b.mean_mechanical_power_w,
            -2.0 * alpha * alpha * m * square,
            2e-9,
        );
        near(
            b.reactive_mechanical_power_var,
            -alpha * omega * m * square,
            2e-9,
        );
        near(
            b.mean_supplied_power_w,
            0.5 * m * omega * omega * square,
            2e-9,
        );
        near(b.power_residual_w, 0.0, 2e-9);
        near(b.reactive_power_residual_var, 0.0, 2e-9);
        assert!(b.mean_ohmic_power_w > b.mean_supplied_power_w);
        near(
            response
                .cells()
                .iter()
                .map(|c| c.mean_mechanical_power_w)
                .sum(),
            b.mean_mechanical_power_w,
            2e-9,
        );
    }
}

#[test]
fn harmonic_response_agrees_with_refined_transient_induction_after_settling() {
    let (problem, sources) = cube(1, None);
    let initial = problem
        .solve(
            &[[4.0; 3]; 6],
            &trace(&problem, [0.0; 3], false),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    let original =
        SpatialInduction::new(&initial, sources, preparation(), 2e-9, &mut || false).unwrap();
    let omega = 3.0;
    let ac = original
        .solve_harmonic(
            &[[C64::new(4.0, 0.0); 3]; 6],
            &[[0.0; 3]; 8],
            omega,
            harmonic_controls(),
            &mut || false,
        )
        .unwrap();
    let expected = ac.potential_wb()[diagonal(&problem)].re;
    let period = std::f64::consts::TAU / omega;
    let mut errors = vec![];
    for samples in [24, 48, 96] {
        let mut state = original.clone();
        let dt = period / samples as f64;
        for step in 1..=3 * samples {
            let j = 4.0 * (omega * dt * step as f64).cos();
            state = state
                .step_candidate(&[[j; 3]; 6], &[[0.0; 3]; 8], dt, controls(), &mut || false)
                .unwrap()
                .into_parts()
                .0;
        }
        errors.push((state.potential_wb()[diagonal(&problem)] - expected).abs());
    }
    for pair in errors.windows(2) {
        let order = (pair[0] / pair[1]).log2();
        assert!(
            (1.8..2.2).contains(&order),
            "errors {errors:?}, order {order}"
        );
    }
    assert!(errors[2] < 5e-4, "{errors:?}");
}

#[test]
fn harmonic_phase_and_dc_gauge_do_not_change_cycle_mean_physics() {
    let (problem, sources) = cube(2, None);
    let zero = vec![[0.0; 3]; problem.mesh().tets.len()];
    let current = vec![[C64::new(0.4, -0.2); 3]; problem.mesh().tets.len()];
    let velocity = vec![[0.0; 3]; problem.mesh().vertex_count];
    let mut outcomes = vec![];
    for gauge in [false, true] {
        let initial = problem
            .solve(
                &zero,
                &trace(&problem, [0.1, 0.2, 0.3], gauge),
                static_controls(),
                &mut || false,
            )
            .unwrap();
        let state =
            SpatialInduction::new(&initial, sources.clone(), preparation(), 2e-9, &mut || {
                false
            })
            .unwrap();
        let phase = if gauge { C64::new(0.0, 1.0) } else { C64::ONE };
        let rotated: Vec<_> = current.iter().map(|j| j.map(|v| v * phase)).collect();
        let response = state
            .solve_harmonic(&rotated, &velocity, 2.0, harmonic_controls(), &mut || false)
            .unwrap();
        outcomes.push((
            response
                .flux_density_t()
                .iter()
                .map(|b| b.map(|v| v / phase))
                .collect::<Vec<_>>(),
            *response.balance(),
        ));
    }
    for (a, b) in outcomes[0].0.iter().zip(&outcomes[1].0) {
        for j in 0..3 {
            complex_near(a[j], b[j], 2e-8);
        }
    }
    near(
        outcomes[0].1.mean_ohmic_power_w,
        outcomes[1].1.mean_ohmic_power_w,
        2e-9,
    );
    near(
        outcomes[0].1.mean_magnetic_energy_j,
        outcomes[1].1.mean_magnetic_energy_j,
        2e-9,
    );
}

#[test]
fn harmonic_resolution_budgets_cancellation_and_unsupported_materials_refuse_atomically() {
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
    let current = [[C64::new(4.0, 0.0); 3]; 6];
    let velocity = [[0.0; 3]; 8];
    for omega in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            state
                .solve_harmonic(&current, &velocity, omega, harmonic_controls(), &mut || {
                    false
                })
                .is_err()
        );
    }
    for settings in [
        SpatialHarmonicInductionControls {
            maximum_cell_skin_ratio: 0.1,
            ..harmonic_controls()
        },
        SpatialHarmonicInductionControls {
            maximum_krylov_scalars: 1,
            ..harmonic_controls()
        },
        SpatialHarmonicInductionControls {
            maximum_field_scalars: 1,
            ..harmonic_controls()
        },
        SpatialHarmonicInductionControls {
            restart: 0,
            ..harmonic_controls()
        },
    ] {
        assert!(
            state
                .solve_harmonic(&current, &velocity, 3.0, settings, &mut || false)
                .is_err()
        );
    }
    let mut polls = 0;
    assert!(matches!(
        state.solve_harmonic(&current, &velocity, 3.0, harmonic_controls(), &mut || {
            polls += 1;
            polls >= 25
        }),
        Err(SpatialMagneticError::Cancelled)
    ));
    assert_eq!(state.time_s(), 0.0);
    assert_eq!(state.potential_wb(), initial.potential_wb());
    let (recoil, conductors) = cube(1, Some(0.2));
    let initial = recoil
        .solve(
            &[[0.0; 3]; 6],
            &trace(&recoil, [0.0; 3], false),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    let state =
        SpatialInduction::new(&initial, conductors, preparation(), 2e-9, &mut || false).unwrap();
    assert!(
        state
            .solve_harmonic(&current, &velocity, 3.0, harmonic_controls(), &mut || false)
            .is_err()
    );
}
