//! Independent moving nonlinear source, work and finite-body energy oracles.
use super::*;
use fs_em::spatial::conduction::motion::SpatialConductionMotion;

// Uniform unit-area cells in series. The imposed laboratory drop is V,
// while the conductor-frame drop is V-2uB. Original source laws and two-cell
// backward Euler balance determine J; no production field/heat result is used.
fn moving_endpoint(
    voltage: f64,
    speed: f64,
    old: [f64; 2],
    r: [f64; 2],
    coordinate: ConductorMagnitude,
) -> (f64, [f64; 2], [f64; 2]) {
    let effective = voltage - 2.0 * speed;
    let field_law = matches!(coordinate, ConductorMagnitude::ElectricField);
    let (mut low, mut high) = (0.0, 4.0);
    for _ in 0..100 {
        let j = (low + high) / 2.0;
        let (_, e) = endpoint(j, 0.1, old, r, [0.05, 0.08], field_law);
        if e[0] + e[1] < effective.abs() {
            low = j;
        } else {
            high = j;
        }
    }
    let magnitude = (low + high) / 2.0;
    let (theta, e) = endpoint(magnitude, 0.1, old, r, [0.05, 0.08], field_law);
    let sign = -effective.signum();
    (sign * magnitude, theta, e.map(|value| sign * value))
}

#[test]
fn moving_nonlinear_sources_change_current_force_and_local_heat_across_two_intervals() {
    let mut observations = Vec::new();
    for (coordinate, rho, speed) in [
        (ConductorMagnitude::ElectricField, 1.0, 0.0),
        (ConductorMagnitude::ElectricField, 1.0, 1.0),
        (ConductorMagnitude::ElectricField, 1.0, 5.0),
        (ConductorMagnitude::ElectricField, 2.0, 1.0),
        (ConductorMagnitude::CurrentDensity, 1.0, 1.0),
    ] {
        let r = [rho, 2.0];
        let (a, ca) = nonlinear_material(r[0], CP[0], ALPHA[0], 0.05, coordinate);
        let (b, cb) = nonlinear_material(r[1], CP[1], ALPHA[1], 0.08, coordinate);
        let (cards, curves) = ([a, b], [ca, cb]);
        let context = query(300.0, coordinate);
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let cells = bindings(&cards, &curves, &context, &key);
        let regions = regions();
        let layout = SpatialConductionThermalNetwork {
            cells: &cells,
            links: &LINKS,
            regions: &regions,
            volume_tolerance_m3: 1e-12,
        };
        let mut thermal = initial(&curves, 999.0, None);
        let mut problem = problem(&cards, coordinate);
        let velocities = [[0.0, speed, 0.0]; 12];
        let magnetic = [[0.0, 0.0, 1.0]; 12];
        let motion = SpatialConductionMotion {
            velocity_m_s: &velocities,
            flux_density_t: &magnetic,
        };
        let mut solve = config();
        solve.maximum_circuit_trials = 4096;
        let mut old = [0.0; 2];
        let mut prior_current: Option<f64> = None;
        for interval in 1..=2 {
            let drive = boundary_response(&problem, 0.0, Some(6.0));
            let before = thermal.energy_j();
            let original = problem.conductors().to_vec();
            let target = interval as f64 * 0.1;
            let (j, theta, e) = moving_endpoint(6.0, speed, old, r, coordinate);
            let mut polls = 0;
            let next = layout
                .step_moving_response_candidate(
                    &thermal,
                    target,
                    &problem,
                    &drive,
                    motion,
                    limits(),
                    field_controls(),
                    controls(),
                    solve,
                    &[0.0; 2],
                    || {
                        polls += 1;
                        false
                    },
                )
                .unwrap();
            let report = next.report();
            let step = &report.electrical;
            let field = &step.field.field;
            let heat = 0.1 * j * (e[0] + e[1]);
            let mechanical = -0.2 * speed * j;
            let supplied = -0.6 * j;
            for (i, cell) in [CellId(30), CellId(10)].into_iter().enumerate() {
                near(
                    next.thermal().temperature_k(cell).unwrap(),
                    300.0 + theta[i],
                );
                // Independent signed per-body balance includes finite
                // conduction to the other body, not only aggregate energy.
                near(
                    CP[i] * (theta[i] - old[i]),
                    0.1 * (j * e[i] - (theta[i] - theta[1 - i])),
                );
                for tet in 6 * i..6 * i + 6 {
                    near(field.current_density_a_m2()[tet][0], j);
                    near(field.conductor_field_v_m()[tet][0], e[i]);
                    near(field.electric_field_v_m()[tet][0], e[i] - speed);
                    near(field.cell_joule_power_w()[tet], j * e[i] / 6.0);
                    near(field.lorentz_force_density_n_m3()[tet][1], -j);
                    near(field.cell_mechanical_power_w()[tet], -speed * j / 6.0);
                    let source = &step.field.evaluated[tet].conductor;
                    near(source.material().query().axes()["T"], 300.0 + theta[i]);
                    let magnitude = match coordinate {
                        ConductorMagnitude::ElectricField => e[i].abs(),
                        ConductorMagnitude::CurrentDensity => j.abs(),
                    };
                    near(source.material().query().axes()["operating"], magnitude);
                    assert_eq!(source.source_key().name(), "original_nonlinear_transport");
                    assert!(source.has_same_electrical_source(&original[tet]));
                }
            }
            near(step.heat_j, heat);
            near(step.supplied_work_j, supplied);
            near(step.mechanical_work_j, mechanical);
            near(step.net_supplied_work_j, supplied - mechanical);
            near(step.net_supplied_work_j, heat);
            near(field.terminal_power_w(), supplied / 0.1);
            near(field.power_balance_residual_w(), 0.0);
            near(next.thermal().energy_j() - before, heat);
            near(report.deposited_joule_j, heat);
            near(report.total_energy_residual_j, 0.0);
            assert!(step.field.field_solves > 1);
            assert!(report.circuit_trials > step.field.field_solves);
            assert_eq!(next.problem().relaxation(), 0.5);
            assert_eq!(problem.conductors(), original);
            near(thermal.time_s(), target - 0.1);
            if speed == 5.0 {
                assert!(step.supplied_work_j < 0.0);
                assert!(step.mechanical_work_j < 0.0);
                assert!(step.heat_j > 0.0);
            }
            if let Some(previous) = prior_current {
                assert!((j - previous).abs() > 1e-6);
            }
            prior_current = Some(j);
            if interval == 1 {
                observations.push([j, heat, -j]);
            }
            // Repeat the same real candidate with a bounded work refusal and
            // cancellation at its last checkpoint: neither history may publish.
            if interval == 1 && speed == 1.0 && rho == 1.0 {
                let mut bounded = solve;
                bounded.maximum_circuit_trials = 1;
                assert!(matches!(
                    layout.step_moving_response_candidate(
                        &thermal,
                        target,
                        &problem,
                        &drive,
                        motion,
                        limits(),
                        field_controls(),
                        controls(),
                        bounded,
                        &[0.0; 2],
                        || false,
                    ),
                    Err(CircuitEnthalpyError::CircuitBudget)
                ));
                let mut seen = 0;
                assert!(matches!(
                    layout.step_moving_response_candidate(
                        &thermal,
                        target,
                        &problem,
                        &drive,
                        motion,
                        limits(),
                        field_controls(),
                        controls(),
                        solve,
                        &[0.0; 2],
                        || {
                            seen += 1;
                            seen == polls
                        },
                    ),
                    Err(CircuitEnthalpyError::Cancelled)
                ));
                assert_eq!(thermal.energy_j(), before);
                assert_eq!(problem.conductors(), original);
            }
            (thermal, problem, _) = next.into_parts();
            old = theta;
        }
        assert_eq!(context.axes()["operating"], 0.0);
    }
    // Speed, the original resistivity claim and the source coordinate each
    // change actual current, native force and the heat reaching finite bodies.
    for alternative in [0, 3, 4] {
        for quantity in 0..3 {
            assert!((observations[1][quantity] - observations[alternative][quantity]).abs() > 1e-6);
        }
    }
}
