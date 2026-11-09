//! G0/G1 physical dissipation-to-heat checks using the existing independent
//! dense thermal oracle, not the production heat assembly or budget helpers.

use super::*;
use crate::rans_solver::{RansAdaptiveControls, RansAdaptiveWorkspace, RansFieldTolerances};

fn load() -> RansHeatLoad<'static> {
    RansHeatLoad {
        coupling: RansHeatCoupling::Dissipative,
        ..RansHeatLoad::default()
    }
}

pub(super) fn physical_heating(
    solver: &RansSolver,
    grid: &RansChannelGrid,
    old: &RansFieldState,
    new: &RansFieldState,
) -> (Vec<f64>, f64, f64) {
    let n = old.u.len();
    let mut volumes = vec![0.0; n];
    let mut power = vec![0.0; n];
    let mut molecular = 0.0;
    for face in 0..n - 1 {
        let dy = grid.y_nodes[face + 1] - grid.y_nodes[face];
        let gradient = (new.u[face + 1] - new.u[face]) / dy;
        let loss = solver.rho * solver.nu_laminar * gradient.powi(2) * dy;
        molecular += loss;
        for i in [face, face + 1] {
            volumes[i] += dy / 2.0;
            power[i] += loss / 2.0;
        }
    }
    let mut turbulent = 0.0;
    for i in 1..n {
        if old.k[i] > 0.0 {
            // Differentiate the interpolating quadratic in sqrt(k). The
            // centerline derivative is zero by the authored symmetry trace.
            let derivative = if i == n - 1 {
                0.0
            } else {
                let a = grid.y_nodes[i] - grid.y_nodes[i - 1];
                let b = grid.y_nodes[i + 1] - grid.y_nodes[i];
                -b * old.k[i - 1].sqrt() / (a * (a + b))
                    + (b - a) * old.k[i].sqrt() / (a * b)
                    + a * old.k[i + 1].sqrt() / (b * (a + b))
            };
            let physical_epsilon = old.eps[i] + 2.0 * solver.nu_laminar * derivative.powi(2);
            let loss = solver.rho * physical_epsilon * new.k[i] / old.k[i] * volumes[i];
            turbulent += loss;
            power[i] += loss;
        }
    }
    for (q, volume) in power.iter_mut().zip(volumes) {
        *q /= volume;
    }
    (power, molecular, turbulent)
}

fn inventory_changes(
    solver: &RansSolver,
    grid: &RansChannelGrid,
    old: &RansFieldState,
    new: &RansFieldState,
) -> (f64, f64) {
    let mechanical: Vec<_> = (0..old.u.len())
        .map(|i| 0.5 * (new.u[i] + old.u[i]) * (new.u[i] - old.u[i]) + new.k[i] - old.k[i])
        .collect();
    let thermal: Vec<_> = old
        .temp_k
        .iter()
        .zip(&new.temp_k)
        .map(|(a, b)| b - a)
        .collect();
    (
        solver.rho * integral(grid, &mechanical),
        solver.rho * solver.cp * integral(grid, &thermal),
    )
}

#[test]
fn local_physical_heat_matches_dense_temperature_and_total_energy() {
    let grid = grid();
    let solver = solver();
    for turbulent in [false, true] {
        for dt in [0.02, 0.5, 2.0] {
            let mut initial = seeded(&grid);
            initial.temp_k.fill(5.0);
            if !turbulent {
                initial.k.fill(0.0);
                initial.eps.fill(0.0);
            }
            let pressure = -0.8;
            let authored = [0.2, -0.1, 0.3, 0.0, -0.05, 0.1, 0.25];
            let wall = 0.07;
            let mut native = initial.clone();
            let mut work = RansTransportWorkspace::try_new(&grid, 7).unwrap();
            let mechanics = solver
                .advance_launder_sharma(&grid, &mut native, pressure, dt, &mut work, || false)
                .unwrap();
            let (mut heating, molecular, destruction) =
                physical_heating(&solver, &grid, &initial, &native);
            for (q, external) in heating.iter_mut().zip(authored) {
                *q += external;
            }
            let expected = dense_heat(
                &solver,
                &grid,
                &initial.temp_k,
                &native.nu_t,
                wall,
                Some(&heating),
                dt,
            );
            let mut actual = initial.clone();
            let report = solver
                .advance_launder_sharma_with_heat(
                    &grid,
                    &mut actual,
                    pressure,
                    dt,
                    RansHeatLoad {
                        wall_heat_flux_w_m2: wall,
                        volumetric_heating_w_m3: Some(&authored),
                        ..load()
                    },
                    &mut work,
                    || false,
                )
                .unwrap();
            assert_eq!(&bits(&actual)[..4], &bits(&native)[..4]);
            assert_eq!(report.mechanical, mechanics);
            for (&value, expected) in actual.temp_k.iter().zip(expected) {
                close(value, expected, 5e-13);
            }
            close(report.molecular_heat_input_w_m2, molecular, 5e-14);
            close(report.turbulent_heat_input_w_m2, destruction, 5e-14);
            let (mechanical_change, thermal_change) =
                inventory_changes(&solver, &grid, &initial, &actual);
            let imposed = wall + integral(&grid, &authored);
            close(
                thermal_change / dt,
                imposed + molecular + destruction,
                5e-11,
            );

            // Independent full inventory check: pressure work and authored
            // heating enter; temporal loss and k diffusion through the wall
            // leave. Neither modeled production nor dissipation is external.
            let old_nut = old_viscosities(&solver, &initial);
            let wall_k_loss = solver.rho
                * (solver.nu_laminar + old_nut[1] / (2.0 * solver.coeffs.sigma_k))
                * actual.k[1]
                / (grid.y_nodes[1] - grid.y_nodes[0]);
            let du_squared: Vec<_> = initial
                .u
                .iter()
                .zip(&actual.u)
                .map(|(a, b)| (b - a).powi(2))
                .collect();
            let temporal = solver.rho * integral(&grid, &du_squared) / (2.0 * dt);
            let net = -pressure * integral(&grid, &actual.u) + imposed - wall_k_loss - temporal;
            assert!(((mechanical_change + thermal_change) / dt - net).abs() < 2e-10);
            assert!(report.thermal_balance_defect_w_m2.abs() < 2e-10);
        }
    }
}

#[test]
fn eddy_transfer_and_numerical_or_wall_losses_are_not_deposited_again() {
    let solver = solver();
    let grid = grid();
    let mut initial = seeded(&grid);
    initial.temp_k.fill(5.0);
    let mut passive = initial.clone();
    let mut active = initial.clone();
    let mut work = RansTransportWorkspace::try_new(&grid, 7).unwrap();
    let p = solver
        .advance_launder_sharma_with_heat(
            &grid,
            &mut passive,
            -0.8,
            0.5,
            RansHeatLoad::default(),
            &mut work,
            || false,
        )
        .unwrap();
    let a = solver
        .advance_launder_sharma_with_heat(&grid, &mut active, -0.8, 0.5, load(), &mut work, || {
            false
        })
        .unwrap();
    assert_eq!(p.mechanical, a.mechanical);
    assert_eq!(&bits(&passive)[..4], &bits(&active)[..4]);
    assert_eq!(passive.temp_k, initial.temp_k);
    assert_eq!(p.molecular_heat_input_w_m2, 0.0);
    assert_eq!(p.turbulent_heat_input_w_m2, 0.0);
    let (_, heat) = inventory_changes(&solver, &grid, &initial, &active);
    let physical =
        solver.rho * (a.mechanical.molecular_dissipation + a.mechanical.turbulent_destruction);
    close(heat / 0.5, physical, 5e-11);
    for excluded in [
        a.mechanical.turbulence_production,
        a.mechanical.temporal_dissipation,
        a.mechanical.turbulent_wall_loss,
    ] {
        assert!(excluded > 1e-6);
        assert!((heat / 0.5 - (physical + solver.rho * excluded)).abs() > 1e-6);
    }
    // With no mechanical motion or turbulence, both coupling choices coincide.
    let mut quiet = laminar(&grid, 5.0);
    let unchanged = bits(&quiet);
    let q = solver
        .advance_launder_sharma_with_heat(&grid, &mut quiet, 0.0, 0.5, load(), &mut work, || false)
        .unwrap();
    assert_eq!(bits(&quiet), unchanged);
    assert_eq!(q.molecular_heat_input_w_m2, 0.0);
    assert_eq!(q.turbulent_heat_input_w_m2, 0.0);
}

#[test]
fn adaptive_rejected_trials_do_not_deposit_extra_dissipation() {
    let solver = solver();
    let grid = grid();
    let mut initial = seeded(&grid);
    initial.temp_k.fill(5.0);
    let mut state = initial.clone();
    let duration = 0.25;
    let source = [0.1; 7];
    let heat = RansHeatLoad {
        wall_heat_flux_w_m2: 0.05,
        volumetric_heating_w_m3: Some(&source),
        ..load()
    };
    let controls = RansAdaptiveControls {
        initial_interval_s: duration,
        minimum_interval_s: 1e-10,
        maximum_interval_s: duration,
        maximum_attempts: 20_000,
        relative_tolerance: 1e-3,
        absolute: RansFieldTolerances {
            velocity_m_s: 1e-5,
            kinetic_energy_m2_s2: 1e-7,
            dissipation_m2_s3: 1e-7,
            eddy_viscosity_m2_s: 1e-7,
            temperature_k: 1e-5,
        },
        maximum_fractional_change: 0.3,
    };
    let mut work = RansAdaptiveWorkspace::try_new(&grid, 7).unwrap();
    let report = solver
        .advance_launder_sharma_adaptive(
            &grid,
            &mut state,
            -0.8,
            duration,
            Some(heat),
            controls,
            &mut work,
            || false,
        )
        .unwrap();
    assert!(report.rejected_intervals > 0);
    assert!(report.accepted_intervals > 1);
    assert!(report.maximum_accepted_error_ratio <= 1.0);
    let h = report.heat.unwrap();
    let m = report.mechanical;
    close(
        h.molecular_heat_j_m2,
        solver.rho * duration * m.molecular_dissipation,
        5e-14,
    );
    close(
        h.turbulent_heat_j_m2,
        solver.rho * duration * m.turbulent_destruction,
        5e-14,
    );
    let (mechanical_change, thermal_change) = inventory_changes(&solver, &grid, &initial, &state);
    let generated = h.molecular_heat_j_m2 + h.turbulent_heat_j_m2;
    let imposed = duration * (heat.wall_heat_flux_w_m2 + integral(&grid, &source));
    assert!((thermal_change - imposed - generated).abs() < 2e-10);
    let net = imposed
        + solver.rho
            * duration
            * (m.pressure_power - m.temporal_dissipation - m.turbulent_wall_loss);
    assert!((mechanical_change + thermal_change - net).abs() < 2e-10);
    assert!(h.balance_defect_j_m2.abs() < 2e-10);
}

#[test]
fn every_dissipative_cancellation_rolls_back_and_reuses_all_nine_buffers() {
    let solver = solver();
    let grid = grid();
    let initial = seeded(&grid);
    let mut work = RansTransportWorkspace::try_new(&grid, 7).unwrap();
    let addresses = |w: &RansTransportWorkspace| {
        [
            &w.velocity,
            &w.kinetic,
            &w.dissipation,
            &w.viscosity,
            &w.production,
            &w.lower,
            &w.diagonal,
            &w.upper,
            &w.rhs,
        ]
        .map(|v| (v.as_ptr() as usize, v.len(), v.capacity()))
    };
    let before = addresses(&work);
    let mut expected = initial.clone();
    let mut polls = 0;
    let report = solver
        .advance_launder_sharma_with_heat(
            &grid,
            &mut expected,
            -0.8,
            0.5,
            load(),
            &mut work,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    for fail_at in 0..polls {
        let mut state = initial.clone();
        let mut seen = 0;
        let error = solver.advance_launder_sharma_with_heat(
            &grid,
            &mut state,
            -0.8,
            0.5,
            load(),
            &mut work,
            || {
                let stop = seen == fail_at;
                seen += 1;
                stop
            },
        );
        assert_eq!(error, Err(RansTransportError::Cancelled));
        assert_eq!(bits(&state), bits(&initial));
        assert_eq!(addresses(&work), before);
        let retried = solver
            .advance_launder_sharma_with_heat(
                &grid,
                &mut state,
                -0.8,
                0.5,
                load(),
                &mut work,
                || false,
            )
            .unwrap();
        assert_eq!(retried, report);
        assert_eq!(bits(&state), bits(&expected));
        assert_eq!(addresses(&work), before);
    }
}

#[test]
fn unrepresentable_generated_heat_refuses_without_publishing_mechanics() {
    let mut solver = solver();
    let grid = grid();
    solver.rho = 1e300;
    solver.cp = 1e-300; // finite ordinary heat capacity per volume
    let mut initial = laminar(&grid, 5.0);
    for (u, y) in initial.u.iter_mut().zip(&grid.y_nodes) {
        *u = 1e8 * y;
    }
    let mut work = RansTransportWorkspace::try_new(&grid, 7).unwrap();
    let mut state = initial.clone();
    assert!(matches!(
        solver.advance_launder_sharma_with_heat(
            &grid,
            &mut state,
            0.0,
            0.5,
            load(),
            &mut work,
            || false
        ),
        Err(RansTransportError::Arithmetic { .. })
    ));
    assert_eq!(bits(&state), bits(&initial));
    let passive = solver
        .advance_launder_sharma_with_heat(
            &grid,
            &mut state,
            0.0,
            0.5,
            RansHeatLoad::default(),
            &mut work,
            || false,
        )
        .unwrap();
    assert_eq!(state.temp_k, initial.temp_k);
    assert!(passive.mechanical.molecular_dissipation.is_finite());
    assert!(passive.mechanical.molecular_dissipation > 0.0);
}
