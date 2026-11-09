use super::*;
use crate::rans_solver::{
    RansAdaptiveControls, RansConjugateAdaptiveWorkspace, RansFieldTolerances,
};

fn force(gravity: f64) -> RansBoussinesqForce {
    RansBoussinesqForce::try_new(0.002, 300.0, gravity, 0.1).unwrap()
}

// Independent integrated physical momentum equations, including actual nodal
// buoyancy and lagged LS viscosity. Dense pivoting shares no native row solver.
fn momentum(
    s: &RansSolver,
    g: &RansChannelGrid,
    old: &RansFieldState,
    gravity: f64,
    pressure: f64,
    dt: f64,
) -> Vec<f64> {
    let n = g.y_nodes.len();
    let nut = old_viscosities(s, old);
    let mut volume = vec![0.0; n];
    let mut a = vec![vec![0.0; n - 1]; n - 1];
    for i in 0..n - 1 {
        let dy = g.y_nodes[i + 1] - g.y_nodes[i];
        volume[i] += 0.5 * dy;
        volume[i + 1] += 0.5 * dy;
        let conductance = (s.nu_laminar + 0.5 * (nut[i] + nut[i + 1])) / dy;
        a[i][i] += conductance;
        if i > 0 {
            a[i - 1][i - 1] += conductance;
            a[i - 1][i] -= conductance;
            a[i][i - 1] -= conductance;
        }
    }
    let mut rhs: Vec<_> = (1..n)
        .map(|i| {
            a[i - 1][i - 1] += volume[i] / dt;
            volume[i]
                * (old.u[i] / dt - pressure / s.rho - gravity * 0.002 * (old.temp_k[i] - 300.0))
        })
        .collect();
    for col in 0..n - 1 {
        let pivot = (col..n - 1)
            .max_by(|&i, &j| a[i][col].abs().total_cmp(&a[j][col].abs()))
            .unwrap();
        a.swap(col, pivot);
        rhs.swap(col, pivot);
        let diagonal = a[col][col];
        for j in col..n - 1 {
            a[col][j] /= diagonal;
        }
        rhs[col] /= diagonal;
        for row in 0..n - 1 {
            if row != col {
                let factor = a[row][col];
                for j in col..n - 1 {
                    a[row][j] -= factor * a[col][j];
                }
                rhs[row] -= factor * rhs[col];
            }
        }
    }
    let mut u = vec![0.0];
    u.extend(rhs);
    u
}

#[test]
fn g0_buoyancy_drives_native_momentum_and_closes_external_work_without_heat_duplication() {
    let s = solver();
    let g = grid();
    for turbulent in [false, true] {
        let mut initial = seeded(&g);
        if !turbulent {
            initial.k.fill(0.0);
            initial.eps.fill(0.0);
        }
        for gravity in [-9.81, 9.81] {
            let mut state = initial.clone();
            let expected = momentum(&s, &g, &initial, gravity, -0.8, 0.15);
            let report = s
                .advance_launder_sharma_with_heat(
                    &g,
                    &mut state,
                    -0.8,
                    0.15,
                    RansHeatLoad {
                        coupling: RansHeatCoupling::Dissipative,
                        buoyancy: Some(force(gravity)),
                        ..Default::default()
                    },
                    &mut RansTransportWorkspace::try_new(&g, 7).unwrap(),
                    || false,
                )
                .unwrap();
            for (&actual, expected) in state.u.iter().zip(expected) {
                close(actual, expected, 3e-12);
            }
            let power: Vec<_> = initial
                .temp_k
                .iter()
                .zip(&state.u)
                .map(|(&t, &u)| -gravity * 0.002 * (t - 300.0) * u)
                .collect();
            close(
                report.mechanical.buoyancy_power,
                integral(&g, &power),
                2e-13,
            );
            assert!(report.mechanical.mean_balance_defect.abs() < 3e-13);
            assert!(report.mechanical.turbulent_balance_defect.abs() < 3e-13);
            let generated = s.rho
                * (report.mechanical.molecular_dissipation
                    + report.mechanical.turbulent_destruction);
            close(report.thermal_energy_change_j_m2, 0.15 * generated, 5e-10);
            assert!(report.thermal_balance_defect_w_m2.abs() < 3e-11);
            let mut passive = initial.clone();
            s.advance_launder_sharma_with_heat(
                &g,
                &mut passive,
                -0.8,
                0.15,
                RansHeatLoad {
                    coupling: RansHeatCoupling::Dissipative,
                    ..Default::default()
                },
                &mut RansTransportWorkspace::try_new(&g, 7).unwrap(),
                || false,
            )
            .unwrap();
            assert!((state.u[6] - passive.u[6]).abs() > 0.02);
            if turbulent {
                assert_ne!(state.k, passive.k);
                assert_ne!(state.nu_t, passive.nu_t);
                assert_ne!(state.temp_k, passive.temp_k);
            }
        }
    }
}

#[test]
fn g3_zero_force_is_exact_and_reference_pressure_shift_preserves_flow() {
    let s = solver();
    let g = grid();
    let initial = seeded(&g);
    let mut outputs = Vec::new();
    for buoyancy in [
        None,
        Some(RansBoussinesqForce::try_new(0.0, 300.0, -9.81, 0.1).unwrap()),
    ] {
        let mut state = initial.clone();
        let report = s
            .advance_launder_sharma_with_heat(
                &g,
                &mut state,
                -0.8,
                0.15,
                RansHeatLoad {
                    buoyancy,
                    ..Default::default()
                },
                &mut RansTransportWorkspace::try_new(&g, 7).unwrap(),
                || false,
            )
            .unwrap();
        outputs.push((bits(&state), report));
    }
    assert_eq!(outputs[0], outputs[1]);
    let mut states = Vec::new();
    for reference in [300.0, 304.0] {
        let mut state = initial.clone();
        s.advance_launder_sharma_with_heat(
            &g,
            &mut state,
            -0.8 + s.rho * -9.81 * 0.002 * (reference - 300.0),
            0.15,
            RansHeatLoad {
                buoyancy: Some(RansBoussinesqForce::try_new(0.002, reference, -9.81, 0.1).unwrap()),
                ..Default::default()
            },
            &mut RansTransportWorkspace::try_new(&g, 7).unwrap(),
            || false,
        )
        .unwrap();
        states.push(state);
    }
    for (a, b) in [
        &states[0].u,
        &states[0].k,
        &states[0].eps,
        &states[0].nu_t,
        &states[0].temp_k,
    ]
    .into_iter()
    .zip([
        &states[1].u,
        &states[1].k,
        &states[1].eps,
        &states[1].nu_t,
        &states[1].temp_k,
    ]) {
        for (&a, &b) in a.iter().zip(b) {
            close(a, b, 2e-12);
        }
    }
}

#[test]
fn g1_vertical_thermal_channel_converges_to_quartic_velocity() {
    let s = solver();
    let mut errors = Vec::new();
    for cells in [16, 32, 64] {
        let g = RansChannelGrid {
            num_cells: cells,
            half_height_m: 1.0,
            y_nodes: (0..=cells).map(|i| i as f64 / cells as f64).collect(),
        };
        let mut state = laminar(&g, 300.0);
        let exact: Vec<_> = g
            .y_nodes
            .iter()
            .map(|&y| (0.02 / s.nu_laminar) * (2.0 * y / 3.0 - y.powi(3) / 3.0 + y.powi(4) / 12.0))
            .collect();
        state.u.clone_from(&exact);
        for (&y, t) in g.y_nodes.iter().zip(&mut state.temp_k) {
            *t = 300.0 + 10.0 * y * (2.0 - y);
        }
        let heating = vec![20.0 * s.k_thermal; cells + 1];
        s.advance_launder_sharma_with_heat(
            &g,
            &mut state,
            0.0,
            1e5,
            RansHeatLoad {
                wall_heat_flux_w_m2: -20.0 * s.k_thermal,
                volumetric_heating_w_m3: Some(&heating),
                buoyancy: Some(force(-1.0)),
                ..Default::default()
            },
            &mut RansTransportWorkspace::try_new(&g, cells + 1).unwrap(),
            || false,
        )
        .unwrap();
        errors.push(
            state
                .u
                .iter()
                .zip(exact)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0_f64, f64::max),
        );
        for (&y, &t) in g.y_nodes.iter().zip(&state.temp_k) {
            assert!((t - (300.0 + 10.0 * y * (2.0 - y))).abs() < 1e-8);
        }
        assert!(state.k.iter().chain(&state.eps).all(|&x| x == 0.0));
    }
    for pair in errors.windows(2) {
        let order = (pair[0] / pair[1]).log2();
        assert!(
            order > 1.95 && order < 2.05,
            "errors {errors:?}, order {order}"
        );
    }
}

fn controls(duration: f64) -> RansAdaptiveControls {
    RansAdaptiveControls {
        initial_interval_s: duration,
        minimum_interval_s: 1e-12,
        maximum_interval_s: duration,
        maximum_attempts: 10000,
        relative_tolerance: 0.0,
        absolute: RansFieldTolerances {
            velocity_m_s: 1.0,
            kinetic_energy_m2_s2: 1.0,
            dissipation_m2_s3: 1.0,
            eddy_viscosity_m2_s: 1.0,
            temperature_k: 1.0,
        },
        maximum_fractional_change: 0.5,
    }
}

#[test]
fn g0_conjugate_adaptive_buoyancy_replays_thermal_feedback_and_integrates_only_accepted_work() {
    let s = solver();
    let g = grid();
    let initial = seeded(&g);
    let wall = || {
        RansConjugateWall::try_new(
            &[RansWallLayer {
                thickness_m: 0.1,
                volumetric_heat_capacity_j_m3_k: 10.0,
                conductivity_w_m_k: 0.2,
            }],
            &[330.0],
            1,
        )
        .unwrap()
    };
    let load = RansConjugateHeatLoad {
        exterior_heat_flux_w_m2: 0.5,
        coupling: RansHeatCoupling::Dissipative,
        buoyancy: Some(force(-9.81)),
        ..Default::default()
    };
    let duration = 0.01;
    let mut expected = initial.clone();
    let mut expected_wall = wall();
    let mut kernel = RansConjugateWorkspace::try_new(&g, &expected_wall, 7, 1).unwrap();
    let first = s
        .advance_launder_sharma_with_conjugate_wall(
            &g,
            &mut expected,
            &mut expected_wall,
            -0.8,
            duration / 2.0,
            load,
            &mut kernel,
            || false,
        )
        .unwrap();
    let second = s
        .advance_launder_sharma_with_conjugate_wall(
            &g,
            &mut expected,
            &mut expected_wall,
            -0.8,
            duration / 2.0,
            load,
            &mut kernel,
            || false,
        )
        .unwrap();
    assert_ne!(
        first.fluid.mechanical.buoyancy_power,
        second.fluid.mechanical.buoyancy_power
    );
    let mut state = initial;
    let mut solid = wall();
    let mut work = RansConjugateAdaptiveWorkspace::try_new(&g, &solid, 7, 1).unwrap();
    let report = s
        .advance_launder_sharma_conjugate_adaptive(
            &g,
            &mut state,
            &mut solid,
            -0.8,
            duration,
            load,
            controls(duration),
            &mut work,
            || false,
        )
        .unwrap();
    assert_eq!(report.fluid.accepted_intervals, 1);
    assert_eq!(bits(&state), bits(&expected));
    assert_eq!(solid.temperatures_k(), expected_wall.temperatures_k());
    close(
        report.fluid.mechanical.buoyancy_power,
        0.5 * (first.fluid.mechanical.buoyancy_power + second.fluid.mechanical.buoyancy_power),
        2e-14,
    );
    assert!(report.fluid.mechanical.mean_balance_defect.abs() < 1e-12);
}

#[test]
fn g4_buoyancy_admission_endpoint_limit_and_every_cancellation_leave_state_untouched() {
    for (beta, reference, gravity, maximum) in [
        (-0.01, 300.0, 1.0, 0.1),
        (0.01, 0.0, 1.0, 0.1),
        (0.01, 300.0, f64::INFINITY, 0.1),
        (0.01, 300.0, 1.0, 0.11),
        (f64::NAN, 300.0, 1.0, 0.1),
        (0.01, 300.0, 1.0, 0.0),
    ] {
        assert!(RansBoussinesqForce::try_new(beta, reference, gravity, maximum).is_err());
    }
    let s = solver();
    let g = grid();
    let initial = seeded(&g);
    let load = RansHeatLoad {
        buoyancy: Some(force(-9.81)),
        ..Default::default()
    };
    let mut work = RansTransportWorkspace::try_new(&g, 7).unwrap();
    let mut expected = initial.clone();
    let mut polls = 0;
    let expected_report = s
        .advance_launder_sharma_with_heat(&g, &mut expected, -0.8, 0.15, load, &mut work, || {
            polls += 1;
            false
        })
        .unwrap();
    for fail in 0..polls {
        let mut state = initial.clone();
        let mut seen = 0;
        assert_eq!(
            s.advance_launder_sharma_with_heat(&g, &mut state, -0.8, 0.15, load, &mut work, || {
                let stop = seen == fail;
                seen += 1;
                stop
            },),
            Err(RansTransportError::Cancelled)
        );
        assert_eq!(bits(&state), bits(&initial));
        let report = s
            .advance_launder_sharma_with_heat(&g, &mut state, -0.8, 0.15, load, &mut work, || false)
            .unwrap();
        assert_eq!(bits(&state), bits(&expected));
        assert_eq!(report, expected_report);
    }
    for temperature in [300.0, 351.0] {
        let mut state = laminar(&g, temperature);
        let before = bits(&state);
        let source = vec![6000.0; 7];
        assert!(matches!(
            s.advance_launder_sharma_with_heat(
                &g,
                &mut state,
                0.0,
                0.1,
                RansHeatLoad {
                    volumetric_heating_w_m3: Some(&source),
                    ..load
                },
                &mut work,
                || false,
            ),
            Err(RansTransportError::Arithmetic {
                equation: "Boussinesq density-departure limit",
                ..
            })
        ));
        assert_eq!(bits(&state), before);
    }
}
