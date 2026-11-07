use super::*;
use fs_em::conduction::ConductionPowerSource;

// Independent scalar elimination of the TWO finite-capacity heat equations
// and each original material law, followed by the physical relation P=I V.
fn power_reference(
    power: f64,
    dt: f64,
    old: [f64; 2],
    field_law: bool,
) -> (f64, [f64; 2], [f64; 2]) {
    let (mut low, mut high) = (0.0, 4.0);
    for _ in 0..100 {
        let j = (low + high) / 2.0;
        let (_, e) = endpoint_at_current(j, dt, old, [1.0, 2.0], [0.05, 0.08], field_law);
        if j * (e[0] + e[1]) < power {
            low = j;
        } else {
            high = j;
        }
    }
    let j = (low + high) / 2.0;
    let (theta, e) = endpoint_at_current(j, dt, old, [1.0, 2.0], [0.05, 0.08], field_law);
    (j, theta, e)
}

#[test]
fn nonlinear_regulated_power_heats_the_actual_material_and_continues_its_history() {
    for field_law in [true, false] {
        let coordinate = if field_law {
            ConductorMagnitude::ElectricField
        } else {
            ConductorMagnitude::CurrentDensity
        };
        let (a, ca) = material(1.0, CP[0], ALPHA[0], 0.05, coordinate, false);
        let (b, cb) = material(2.0, CP[1], ALPHA[1], 0.08, coordinate, false);
        let (cards, curves) = ([a, b], [ca, cb]);
        let context = point(300.0, coordinate);
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let bindings = cells(&cards, &curves, &context, &key);
        let layout = network(&bindings);
        let mut thermal = initial(&curves, 999.0, None);
        let mut field = problem(&cards, coordinate).with_relaxation(0.5).unwrap();
        let source = ConductionPowerSource {
            power_w: 12.0,
            maximum_voltage_v: 10.0,
            maximum_current_a: 4.0,
        };
        let drive = ConductionDrive::Power {
            positive: &[2, 5],
            negative: &[0, 3],
            source,
        };
        let mut solve = config();
        solve.maximum_circuit_trials = 4096;
        let mut old = [0.0; 2];
        let mut previous_current = f64::INFINITY;
        for step in 1..=2 {
            let (j, theta, e) = power_reference(source.power_w, 0.1, old, field_law);
            let before = thermal.energy_j();
            let candidate = layout
                .step_response_candidate(
                    &thermal,
                    step as f64 * 0.1,
                    &field,
                    drive,
                    limits(),
                    field_controls(),
                    controls(),
                    solve,
                    &[0.0; 2],
                    || false,
                )
                .unwrap();
            let report = candidate.report();
            near(report.electrical.heat_j, 1.2);
            near(report.electrical.supplied_work_j, 1.2);
            near(candidate.thermal().energy_j() - before, 1.2);
            near(report.total_energy_residual_j, 0.0);
            for (i, cell) in [CellId(30), CellId(10)].into_iter().enumerate() {
                near(
                    candidate.thermal().temperature_k(cell).unwrap(),
                    300.0 + theta[i],
                );
                for triangle in 2 * i..2 * i + 2 {
                    near(
                        report.electrical.field.electric_field_v_m[triangle][0],
                        -e[i],
                    );
                    near(
                        report.electrical.field.current_density_a_m2[triangle][0],
                        -j,
                    );
                    near(
                        report.electrical.field.cell_power_w[triangle],
                        0.5 * e[i] * j,
                    );
                }
            }
            assert!(report.electrical.field.field_solves > 1);
            assert!(report.circuit_trials > report.electrical.field.field_solves);
            assert!(j < previous_current);
            previous_current = j;
            old = theta;
            (thermal, field, _) = candidate.into_parts();
        }
    }
}

#[test]
fn power_compliance_at_the_hot_nonlinear_endpoint_refuses_without_publishing() {
    let coordinate = ConductorMagnitude::ElectricField;
    let (a, ca) = material(1.0, CP[0], ALPHA[0], 0.05, coordinate, false);
    let (b, cb) = material(2.0, CP[1], ALPHA[1], 0.08, coordinate, false);
    let (cards, curves) = ([a, b], [ca, cb]);
    let context = point(300.0, coordinate);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let bindings = cells(&cards, &curves, &context, &key);
    let thermal = initial(&curves, 999.0, None);
    let field = problem(&cards, coordinate).with_relaxation(0.5).unwrap();
    let (_, _, cold) = power_reference(12.0, 0.0, [0.0; 2], true);
    let (_, _, hot) = power_reference(12.0, 1.0, [0.0; 2], true);
    let source = ConductionPowerSource {
        power_w: 12.0,
        maximum_voltage_v: (cold.iter().sum::<f64>() + hot.iter().sum::<f64>()) / 2.0,
        maximum_current_a: 4.0,
    };
    // The cold nonlinear electrical problem is physically feasible.
    field
        .solve_power_candidate(
            &[2, 5],
            &[0, 3],
            source,
            limits(),
            field_controls(),
            controls(),
            &mut || false,
        )
        .unwrap();
    let before = thermal.energy_j();
    let mut solve = config();
    solve.maximum_circuit_trials = 4096;
    assert!(
        network(&bindings)
            .step_response_candidate(
                &thermal,
                1.0,
                &field,
                ConductionDrive::Power {
                    positive: &[2, 5],
                    negative: &[0, 3],
                    source
                },
                limits(),
                field_controls(),
                controls(),
                solve,
                &[0.0; 2],
                || false,
            )
            .is_err()
    );
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.energy_j().to_bits(), before.to_bits());
    assert_eq!(field.conductors()[0].material().query().axes()["T"], 300.0);
    assert_eq!(
        field.conductors()[0].material().query().axes()["operating"],
        0.0
    );
}
