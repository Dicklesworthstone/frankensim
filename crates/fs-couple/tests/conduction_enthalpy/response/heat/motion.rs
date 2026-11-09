//! Independent two-body work/heat closure with nonlinear moving material laws.
use super::*;
use fs_couple::circuit_enthalpy::conduction::response::motion::{
    MovingConductionDrive, step_moving_conduction_response_candidate,
};
use fs_em::conduction::motion::ConductionMotion;

// At prescribed positive -Jx, rho_i=r_i+alpha_i*(T_i-300)+beta_i*|E'_i|
// gives |E'_i|=J*(r_i+alpha_i*(T_i-300))/(1-beta_i*J).
// Eliminate the second temperature through TOTAL heat and solve the first
// body's balance, including its original interface and gray-radiation losses.
fn endpoint(j: f64, old: [f64; 2], dt: f64) -> ([f64; 2], [f64; 2], f64) {
    let d = BETA.map(|beta| 1.0 - beta * j);
    let a = [0, 1].map(|i| CP[i] - dt * j * j * ALPHA[i] / d[i]);
    let supply = dt * j * j * (RHO[0] / d[0] + RHO[1] / d[1])
        + CP[0] * (old[0] - 300.0)
        + CP[1] * (old[1] - 300.0);
    let equation = |first: f64| {
        let second = (supply - a[0] * first) / a[1];
        let t = [300.0 + first, 300.0 + second];
        let e = [0, 1].map(|i| j * (RHO[i] + ALPHA[i] * (t[i] - 300.0)) / d[i]);
        let g = 1.0 / (1.0 + 0.2 * (1.0 + 0.001 * first));
        let emissivity = [0.2, 0.6];
        let emissivity = [0, 1].map(|i| emissivity[i] * (1.0 + 0.0001 * (t[i] - 300.0)));
        let transfer = g * (first - second)
            + 5.670_374_419e-8 * (t[0].powi(4) - t[1].powi(4))
                / (1.0 / emissivity[0] + 1.0 / emissivity[1] - 1.0);
        let residual = CP[0] * (t[0] - old[0]) - dt * (j * e[0] - transfer);
        (t, e, transfer, residual)
    };
    let (mut lo, mut hi) = (-20.0, 40.0);
    assert!(equation(lo).3 <= 0.0 && equation(hi).3 >= 0.0);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if equation(mid).3 < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let (t, e, transfer, _) = equation(0.5 * (lo + hi));
    (t, e, transfer)
}

fn reference_moving(
    voltage: f64,
    velocity_y: f64,
    old: [f64; 2],
    dt: f64,
) -> (f64, [f64; 2], [f64; 2], f64) {
    // The two unit-length bodies have Bz=1 T. Their material-frame voltage
    // drop is V_right - V_left - integral(u_y Bz dx), including its sign.
    let driving = voltage - 2.0 * velocity_y;
    assert!(driving > 0.0);
    let (mut lo, mut hi) = (0.0, 3.0);
    assert!(endpoint(hi, old, dt).1.iter().sum::<f64>() > driving);
    for _ in 0..70 {
        let mid = 0.5 * (lo + hi);
        if endpoint(mid, old, dt).1.iter().sum::<f64>() < driving {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let j = 0.5 * (lo + hi);
    let (t, e, transfer) = endpoint(j, old, dt);
    (j, t, e, transfer)
}

#[test]
fn nonlinear_motoring_generation_and_braking_heat_the_original_two_bodies() {
    for (voltage, velocity_y) in [(8.0, 0.5), (-4.0, -5.0), (0.0, -3.0)] {
        with_heat(
            0.2,
            0.2,
            BETA,
            Support::Ordinary,
            |cx, _, assembly, mut thermal, mut field| {
                let selected = assembly.prepare(cx, &thermal).unwrap();
                let velocities = [[0.0, velocity_y]; 4];
                let flux = [1.0; 4];
                let boundary = ends(voltage, 0.0);
                let drive = MovingConductionDrive {
                    boundary_voltage_v: &boundary,
                    current_traces: &[],
                    motion: ConductionMotion {
                        velocity_m_s: &velocities,
                        normal_flux_density_t: &flux,
                    },
                };
                for interval in 1..=2 {
                    let old =
                        [CellId(30), CellId(10)].map(|cell| thermal.temperature_k(cell).unwrap());
                    let old_energy = thermal.energy_j();
                    let (j, expected, e, transfer) =
                        reference_moving(voltage, velocity_y, old, 0.5);
                    let next = step_moving_conduction_response_candidate(
                        &selected,
                        cx,
                        &thermal,
                        0.5 * interval as f64,
                        &field,
                        drive,
                        &[],
                        &REGIONS,
                        1e-12,
                        limits(),
                        field_controls(),
                        controls(),
                        config(),
                        &[0.0; 2],
                    )
                    .unwrap();
                    let report = next.report();
                    let step = &report.electrical;
                    for (i, cell) in [CellId(30), CellId(10)].into_iter().enumerate() {
                        near(next.thermal().temperature_k(cell).unwrap(), expected[i]);
                        // Each body's balance is independent, so swapped heat
                        // routes cannot hide behind a correct total account.
                        let loss = if i == 0 { transfer } else { -transfer };
                        near(CP[i] * (expected[i] - old[i]), 0.5 * (j * e[i] - loss));
                    }
                    for triangle in 0..4 {
                        let body = triangle / 2;
                        near(step.field.field.current_density_a_m2[triangle][0], -j);
                        near(step.field.motion.conductor_field_v_m[triangle][0], -e[body]);
                        near(step.field.motion.lorentz_force_density_n_m3[triangle][1], j);
                        near(
                            step.field.field.evaluated[triangle]
                                .conductor
                                .material()
                                .query()
                                .axes()["operating"],
                            e[body],
                        );
                    }
                    near(step.terminal_work_j, 0.5 * voltage * j);
                    near(step.mechanical_work_j, velocity_y * j);
                    near(step.heat_j, 0.5 * j * (e[0] + e[1]));
                    near(step.net_supplied_work_j, step.heat_j);
                    near(report.deposited_joule_j, step.heat_j);
                    near(next.thermal().energy_j() - old_energy, step.heat_j);
                    assert!(step.heat_j > 0.0);
                    if voltage < 0.0 {
                        assert!(step.terminal_work_j < 0.0);
                        assert!(step.mechanical_work_j < 0.0);
                    } else if voltage == 0.0 {
                        near(step.terminal_work_j, 0.0);
                        near(step.heat_j, -step.mechanical_work_j);
                    } else {
                        assert!(step.mechanical_work_j > 0.0);
                        assert!(step.terminal_work_j > step.heat_j);
                    }
                    assert!(report.circuit_trials > 1);
                    assert!(
                        (report.total_energy_residual_j
                            + 0.5 * report.source.radiation_net_power_w())
                        .abs()
                            < config().total_energy_tolerance_j
                    );
                    assert!(report.source.radiation.is_some());
                    let (accepted_thermal, accepted_field, _) = next.into_parts();
                    thermal = accepted_thermal;
                    field = accepted_field;
                }
                near(thermal.time_s(), 1.0);
            },
        );
    }
}

#[test]
fn moving_heat_refusals_and_budget_do_not_publish_either_history() {
    for support in [
        Support::Ordinary,
        Support::ShortEmissivity,
        Support::EmissivityField,
        Support::ContactField,
    ] {
        with_heat(
            0.2,
            0.2,
            BETA,
            support,
            |cx, gate, assembly, thermal, field| {
                let selected = assembly.prepare(cx, &thermal).unwrap();
                let original_energy = thermal.energy_j();
                let velocity = [[0.0, -3.0]; 4];
                let flux = [1.0; 4];
                let boundary = ends(0.0, 0.0);
                let drive = MovingConductionDrive {
                    boundary_voltage_v: &boundary,
                    current_traces: &[],
                    motion: ConductionMotion {
                        velocity_m_s: &velocity,
                        normal_flux_density_t: &flux,
                    },
                };
                let call = |solve| {
                    step_moving_conduction_response_candidate(
                        &selected,
                        cx,
                        &thermal,
                        0.5,
                        &field,
                        drive,
                        &[],
                        &REGIONS,
                        1e-12,
                        limits(),
                        field_controls(),
                        controls(),
                        solve,
                        &[0.0; 2],
                    )
                };
                if support == Support::Ordinary {
                    let mut exhausted = config();
                    exhausted.maximum_circuit_trials = 1;
                    assert!(matches!(
                        call(exhausted),
                        Err(CircuitEnthalpyError::CircuitBudget)
                    ));
                    // A failed numerical attempt must not poison a valid retry.
                    assert!(call(config()).is_ok());
                } else {
                    assert!(call(config()).is_err());
                }
                gate.request();
                assert!(matches!(
                    call(config()),
                    Err(CircuitEnthalpyError::Cancelled)
                ));
                assert_eq!(thermal.time_s(), 0.0);
                assert_eq!(thermal.energy_j(), original_energy);
                assert!(
                    field.conductors().iter().all(|conductor| {
                        conductor.material().query().axes()["operating"] == 0.0
                    })
                );
            },
        );
    }
}
