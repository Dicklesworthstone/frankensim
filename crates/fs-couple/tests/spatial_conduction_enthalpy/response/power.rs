use super::*;
use fs_em::conduction::ConductionPowerSource;

fn power_reference(
    power: f64,
    dt: f64,
    old: [f64; 2],
    field_law: bool,
) -> (f64, [f64; 2], [f64; 2]) {
    let (mut low, mut high) = (0.0, 4.0);
    for _ in 0..100 {
        let j = (low + high) / 2.0;
        let (_, e) = endpoint(j, dt, old, [1.0, 2.0], [0.05, 0.08], field_law);
        if j * (e[0] + e[1]) < power {
            low = j;
        } else {
            high = j;
        }
    }
    let j = (low + high) / 2.0;
    let (theta, e) = endpoint(j, dt, old, [1.0, 2.0], [0.05, 0.08], field_law);
    (j, theta, e)
}

fn power_faces(problem: &SpatialConductionResponse) -> [Vec<usize>; 2] {
    [0.0, 2.0].map(|x| {
        problem
            .exterior_faces()
            .filter(|&face| {
                problem.mesh().faces[face].map(|v| problem.positions_m()[v as usize][0]) == [x; 3]
            })
            .collect()
    })
}

#[test]
fn spatial_nonlinear_power_preserves_local_heat_source_work_and_continuation() {
    for field_law in [true, false] {
        let coordinate = if field_law {
            ConductorMagnitude::ElectricField
        } else {
            ConductorMagnitude::CurrentDensity
        };
        let (a, ca) = nonlinear_material(1.0, CP[0], ALPHA[0], 0.05, coordinate);
        let (b, cb) = nonlinear_material(2.0, CP[1], ALPHA[1], 0.08, coordinate);
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
        let mut field = problem(&cards, coordinate);
        let faces = power_faces(&field);
        let source = ConductionPowerSource {
            power_w: 12.0,
            maximum_voltage_v: 10.0,
            maximum_current_a: 4.0,
        };
        let drive = SpatialConductionResponseDrive::Power {
            positive_faces: &faces[1],
            negative_faces: &faces[0],
            source,
            reference_voltage_v: 0.5,
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
                for tet in 6 * i..6 * i + 6 {
                    near(
                        report.electrical.field.field.electric_field_v_m[tet][0],
                        -e[i],
                    );
                    near(
                        report.electrical.field.field.current_density_a_m2[tet][0],
                        -j,
                    );
                    near(
                        report.electrical.field.field.cell_power_w[tet],
                        e[i] * j / 6.0,
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
fn hot_spatial_power_compliance_rejects_both_candidate_histories() {
    let coordinate = ConductorMagnitude::ElectricField;
    let (a, ca) = nonlinear_material(1.0, CP[0], ALPHA[0], 0.05, coordinate);
    let (b, cb) = nonlinear_material(2.0, CP[1], ALPHA[1], 0.08, coordinate);
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
    let thermal = initial(&curves, 999.0, None);
    let field = problem(&cards, coordinate);
    let faces = power_faces(&field);
    let (_, _, cold) = power_reference(12.0, 0.0, [0.0; 2], true);
    let (_, _, hot) = power_reference(12.0, 1.0, [0.0; 2], true);
    let source = ConductionPowerSource {
        power_w: 12.0,
        maximum_voltage_v: (cold.iter().sum::<f64>() + hot.iter().sum::<f64>()) / 2.0,
        maximum_current_a: 4.0,
    };
    let drive = SpatialConductionResponseDrive::Power {
        positive_faces: &faces[1],
        negative_faces: &faces[0],
        source,
        reference_voltage_v: 0.5,
    };
    field
        .solve_candidate(drive, limits(), field_controls(), controls(), &mut || false)
        .unwrap();
    let before = thermal.energy_j();
    let mut solve = config();
    solve.maximum_circuit_trials = 4096;
    assert!(
        layout
            .step_response_candidate(
                &thermal,
                1.0,
                &field,
                drive,
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
