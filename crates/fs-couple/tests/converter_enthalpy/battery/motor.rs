//! Complete battery-to-shaft power and its ten physical heat receivers.
use super::*;
use fs_couple::circuit_enthalpy::converter::battery::motor::*;
use fs_phs::circuit::{
    converter::battery_fed::motor::BatteryMotorDrive,
    motor::{DcMotor, DcMotorParameters},
};

const DRIVE_IDS: [CellId; 10] = [
    ALL_IDS[0],
    ALL_IDS[1],
    ALL_IDS[2],
    ALL_IDS[3],
    ALL_IDS[4],
    ALL_IDS[5],
    ALL_IDS[6],
    ALL_IDS[7],
    CellId(139),
    CellId(140),
];
const DRIVE_R: [f64; 9] = [
    RESISTANCES[0],
    RESISTANCES[1],
    RESISTANCES[2],
    RESISTANCES[3],
    RESISTANCES[4],
    RESISTANCES[5],
    RESISTANCES[6],
    RESISTANCES[7],
    0.4,
];

struct DriveFixture {
    cards: [MaterialCard; 10],
    curves: [ResolvedHeatCapacityCurve; 10],
    cp: [f64; 10],
}
impl DriveFixture {
    fn new(cp: [f64; 10], alpha: f64) -> Self {
        let pairs: [_; 10] = core::array::from_fn(|i| {
            // The final body supplies only caloric data; bearing drag is the
            // separate fixed mechanical coefficient, not this unused R claim.
            material(if i < 9 { DRIVE_R[i] } else { 1.0 }, cp[i], alpha)
        });
        Self {
            cards: pairs.each_ref().map(|(card, _)| card.clone()),
            curves: pairs.map(|(_, curve)| curve),
            cp,
        }
    }
    fn body(&self) -> DepositionSimulation<'_> {
        let events = (0..10)
            .map(|i| TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id: DRIVE_IDS[i],
                    segment_id: i as u64,
                    body_id: i as u64,
                    material_id: MaterialId(i as u64),
                    mass_kg: 1.0,
                    specific_h_j_kg: 100.0 * self.cp[i],
                    geometry: CellGeometry {
                        start_m: [i as f64, 0.0, 0.0],
                        end_m: [i as f64 + 1.0, 0.0, 0.0],
                        volume_m3: 1.0,
                        surface_area_m2: 6.0,
                    },
                    boundary: Boundary {
                        ambient_k: 300.0,
                        htc_w_m2_k: 0.0,
                        reservoir: None,
                        source_w: 0.0,
                    },
                }),
            })
            .collect();
        let laws = self
            .curves
            .iter()
            .enumerate()
            .map(|(i, curve)| (MaterialId(i as u64), curve.curve() as &dyn ThermalLaw))
            .collect::<Vec<_>>();
        DepositionSimulation::new(
            events,
            &laws,
            DepositionBudget {
                max_cells: 10,
                max_links: 10,
                max_events: 20,
                max_substeps: 1000,
            },
        )
        .unwrap()
        .advance_to(
            0.0,
            StepConfig {
                max_dt_s: 0.1,
                energy_tolerance_j: 1e-8,
            },
            || Ok(()),
        )
        .unwrap()
        .0
    }
    fn cells<'a>(&'a self, context: &'a QueryPoint) -> [MaterialCell<'a>; 10] {
        core::array::from_fn(|i| MaterialCell {
            cell: DRIVE_IDS[i],
            material_id: MaterialId(i as u64),
            curve: self.curves[i].curve(),
            point: ThermalMaterialPoint {
                material: &self.cards[i],
                context,
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            conductivity: None,
        })
    }
}
fn drive_paths(key: &PropertyKey) -> [ConverterThermalChannel<'_>; 9] {
    core::array::from_fn(|i| ConverterThermalChannel {
        cell: DRIVE_IDS[i],
        resistance: ConverterResistance::Resistance {
            key,
            selection: MaterialPropertySelection::SingleClaimOnly,
        },
    })
}
fn drive_controls() -> CircuitEnthalpySolve {
    let mut solve = controls();
    solve.maximum_cells = 10;
    solve.newton.linear_restart = 10;
    solve
}
fn motor_parameters() -> DcMotorParameters {
    DcMotorParameters {
        inductance_h: 0.05,
        winding_resistance_ohm: DRIVE_R[8],
        rotor_inertia_kg_m2: 0.4,
        stator_inertia_kg_m2: 0.8,
        coupling_nm_per_a: 0.2,
        bearing_drag_nm_s_per_rad: 0.0003,
    }
}

// Independent steady power path: held physical torques keep both shafts at
// finite absolute speed; only the finite battery's chemical coordinate drifts.
fn steady_drive(
    converter_current: f64,
    auxiliary_current: f64,
) -> (
    BatteryMotorDrive,
    [f64; 9],
    [f64; 2],
    Option<RegulatedPowerLoad>,
    [f64; 10],
) {
    let (electrical, _, _, _) = steady();
    let (a, b) = (0.7, 0.6);
    let battery_current = a * converter_current;
    let polarization = battery_current * DRIVE_R[6];
    let battery_bus = 12.0 - polarization - battery_current * DRIVE_R[5];
    let input_bus = battery_bus - battery_current * DRIVE_R[7];
    let factors = [1.0, a, 1.0 - a, b, 1.0 - b];
    let resistance: f64 = (0..5).map(|i| factors[i] * DRIVE_R[i]).sum();
    let output_bus = (a * input_bus - resistance * converter_current) / b;
    let source_state = electrical
        .state(
            0.8,
            polarization,
            battery_bus,
            converter_current,
            [input_bus, output_bus],
        )
        .unwrap();
    let p = motor_parameters();
    let current = b * converter_current - auxiliary_current;
    let relative_speed = (output_bus - p.winding_resistance_ohm * current) / p.coupling_nm_per_a;
    let speeds = [relative_speed + 3.0, 3.0];
    let torque = p.bearing_drag_nm_s_per_rad * relative_speed - p.coupling_nm_per_a * current;
    let motor = DcMotor::new(p).unwrap();
    let motor_state = motor.state(current, speeds).unwrap();
    let native = electrical.with_dc_motor(motor).unwrap();
    let state = native.state(&source_state, &motor_state).unwrap();
    let load = (auxiliary_current > 0.0).then(|| {
        RegulatedPowerLoad::try_new(output_bus * auxiliary_current, 10.0, 1.0, 50.0).unwrap()
    });
    let heat = core::array::from_fn(|i| match i {
        0..=4 => factors[i] * DRIVE_R[i] * converter_current * converter_current,
        5..=7 => DRIVE_R[i] * battery_current * battery_current,
        8 => p.winding_resistance_ohm * current * current,
        _ => p.bearing_drag_nm_s_per_rad * relative_speed * relative_speed,
    });
    (native, state, [torque, -torque], load, heat)
}

#[test]
fn ten_channel_calorimetry_excludes_shaft_work_and_auxiliary_load_export() {
    let fixture = DriveFixture::new([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0], 0.0);
    let context = point(300.0);
    let cells = fixture.cells(&context);
    let resistance = key();
    let channels = drive_paths(&resistance);
    let network = BatteryMotorThermalNetwork {
        cells: &cells,
        links: &[],
        channels: &channels,
        bearing_cell: DRIVE_IDS[9],
    };
    let body = fixture.body();
    let (native, before, torques, load, heat_w) = steady_drive(1.0, 0.1);
    let drive = network.bind(&body, &native, before, || false).unwrap();
    let dt = 0.02;
    let result = network
        .step_candidate(
            &body,
            dt,
            &drive,
            native.electrical().modulation(),
            &torques,
            load,
            electrical_controls(),
            drive_controls(),
            &[0.0; 10],
            || false,
        )
        .unwrap();
    let step = &result.report().electrical;
    for i in 0..10 {
        let expected = dt * heat_w[i];
        let caloric =
            fixture.cp[i] * (result.thermal().temperature_k(DRIVE_IDS[i]).unwrap() - 300.0);
        assert!(
            (caloric - expected).abs() < 2e-9,
            "channel {i}: {caloric} != {expected}"
        );
        near(step.heat_joules()[i], expected);
    }
    near(step.electrical.battery.chemical_energy_change_j, -dt * 8.4);
    near(step.electrical.load_work_j, dt * 1.3385);
    let p = motor_parameters();
    near(
        step.shaft_work_j[0],
        dt * torques[0] * before[7] / p.rotor_inertia_kg_m2,
    );
    near(
        step.shaft_work_j[1],
        dt * torques[1] * before[8] / p.stator_inertia_kg_m2,
    );
    assert!(step.shaft_work_j.iter().sum::<f64>() < 0.0);
    assert!(step.shaft_work_j[0] < 0.0 && step.shaft_work_j[1] > 0.0);
    near(step.motor_electrical_work_j, dt * 13.385 * 0.5);
    near(
        step.electrical.external_work_j,
        -step.motor_electrical_work_j,
    );
    near(
        step.native.supplied,
        step.shaft_work_j.iter().sum::<f64>() - step.electrical.load_work_j,
    );
    near(
        result.report().deposited_joule_j,
        dt * heat_w.iter().sum::<f64>(),
    );
    near(
        result.thermal().energy_j() - body.energy_j(),
        result.report().deposited_joule_j,
    );
    near(result.report().non_electrical_heat_j, 0.0);
    near(result.report().total_energy_residual_j, 0.0);
    for i in 1..9 {
        near(result.drive().state()[i], before[i]);
    }
    assert_eq!(result.report().source.bearing_cell, DRIVE_IDS[9]);
    near(
        result.report().source.bearing_drag_nm_s_per_rad,
        p.bearing_drag_nm_s_per_rad,
    );
}

#[test]
fn finite_rotor_regeneration_charges_battery_and_shared_receiver_sums_both_motor_losses() {
    let fixture = DriveFixture::new([1.0; 10], 0.0);
    let context = point(300.0);
    let cells = fixture.cells(&context);
    let resistance = key();
    let channels = drive_paths(&resistance);
    let network = BatteryMotorThermalNetwork {
        cells: &cells,
        links: &[],
        channels: &channels,
        bearing_cell: DRIVE_IDS[8], // Explicit winding/bearing homogenization.
    };
    let body = fixture.body();
    let (native, before, _, _, _) = steady_drive(-1.0, 0.0);
    let drive = network.bind(&body, &native, before, || false).unwrap();
    let result = network
        .step_candidate(
            &body,
            0.005,
            &drive,
            native.electrical().modulation(),
            &[0.0; 2],
            None,
            electrical_controls(),
            drive_controls(),
            &[0.0; 10],
            || false,
        )
        .unwrap();
    let step = &result.report().electrical;
    let after = result.drive().state();
    assert!(after[0] > before[0], "finite chemical charge must increase");
    assert!(step.electrical.battery.chemical_energy_change_j > 0.0);
    assert!(step.motor_electrical_work_j < 0.0);
    assert!(step.kinetic_energy_change_j.iter().sum::<f64>() < 0.0);
    assert!(after[7] < before[7] && after[8] > before[8]);
    near(after[7] + after[8], before[7] + before[8]);
    near(step.angular_momentum_residual_nms, 0.0);
    near(step.native.supplied, 0.0);
    near(step.electrical.load_work_j, 0.0);
    near(step.shaft_work_j.iter().sum(), 0.0);
    near(
        result.thermal().temperature_k(DRIVE_IDS[8]).unwrap() - 300.0,
        step.copper_heat_j + step.bearing_heat_j,
    );
    near(result.thermal().temperature_k(DRIVE_IDS[9]).unwrap(), 300.0);
    assert!(step.heat_joules().iter().all(|q| *q >= 0.0));
    near(
        result.report().deposited_joule_j,
        step.heat_joules().iter().sum(),
    );
    near(result.report().total_energy_residual_j, 0.0);
}

#[test]
fn nine_resistance_sources_and_fixed_bearing_heat_use_the_same_endpoint_solution() {
    let (_, _, _, _, heat_w) = steady_drive(1.0, 0.1);
    let cp = heat_w.map(|power| 0.02 * power);
    let mut endpoints = Vec::new();
    for alpha in [0.0, 0.0039] {
        let fixture = DriveFixture::new(cp, alpha);
        let context = point(300.0);
        let cells = fixture.cells(&context);
        let resistance = key();
        let channels = drive_paths(&resistance);
        let network = BatteryMotorThermalNetwork {
            cells: &cells,
            links: &[],
            channels: &channels,
            bearing_cell: DRIVE_IDS[9],
        };
        let body = fixture.body();
        let (native, before, torques, load, _) = steady_drive(1.0, 0.1);
        let drive = network.bind(&body, &native, before, || false).unwrap();
        let dt = 0.05;
        let result = network
            .step_candidate(
                &body,
                dt,
                &drive,
                native.electrical().modulation(),
                &torques,
                load,
                electrical_controls(),
                drive_controls(),
                &[0.0; 10],
                || false,
            )
            .unwrap();
        let next = *result.drive().state();
        let r = core::array::from_fn::<_, 9, _>(|i| result.drive().channels()[i].resistance_ohm());
        let p = motor_parameters();
        let battery = native.electrical().battery_parameters();
        let converter_current = f64::midpoint(before[3], next[3]) / parameters().inductance_h;
        let cell_current = -(next[0] - before[0]) / dt;
        let vp = f64::midpoint(before[1], next[1]) / battery.polarization_capacitance_f;
        let vb = f64::midpoint(before[2], next[2]) / battery.bus_capacitance_f;
        let vi = f64::midpoint(before[4], next[4]) / parameters().bus_capacitance_f[0];
        let lead_current = (vb - vi) / r[7];
        let motor_current = f64::midpoint(before[6], next[6]) / p.inductance_h;
        let relative_speed = f64::midpoint(before[7], next[7]) / p.rotor_inertia_kg_m2
            - f64::midpoint(before[8], next[8]) / p.stator_inertia_kg_m2;
        let factors = [1.0, 0.7, 0.3, 0.6, 0.4];
        for i in 0..10 {
            let temperature = result.thermal().temperature_k(DRIVE_IDS[i]).unwrap();
            let power = match i {
                0..=4 => factors[i] * r[i] * converter_current * converter_current,
                5 => r[i] * cell_current * cell_current,
                6 => vp * vp / r[i],
                7 => r[i] * lead_current * lead_current,
                8 => r[i] * motor_current * motor_current,
                _ => p.bearing_drag_nm_s_per_rad * relative_speed * relative_speed,
            };
            assert!(
                (cp[i] * (temperature - 300.0) - dt * power).abs() < 2e-9,
                "alpha={alpha}, channel={i}: independent caloric/source balance"
            );
            if i < 9 {
                near(r[i], DRIVE_R[i] * (1.0 + alpha * (temperature - 300.0)));
                assert_eq!(
                    drive.channels()[i].material().unwrap().card_identity(),
                    result.drive().channels()[i]
                        .material()
                        .unwrap()
                        .card_identity()
                );
            }
        }
        near(
            result.drive().motor_parameters().bearing_drag_nm_s_per_rad,
            p.bearing_drag_nm_s_per_rad,
        );
        near(result.report().total_energy_residual_j, 0.0);
        endpoints.push(next);
    }
    assert!(
        (endpoints[1][6] - endpoints[0][6]).abs() > 1e-7,
        "actual endpoint resistance must alter solved motor flux/current"
    );
}

#[test]
fn refused_receiver_load_and_cancelled_trials_preserve_all_physical_histories() {
    let fixture = DriveFixture::new([1.0; 10], 0.0);
    let context = point(300.0);
    let cells = fixture.cells(&context);
    let resistance = key();
    let channels = drive_paths(&resistance);
    let network = BatteryMotorThermalNetwork {
        cells: &cells,
        links: &[],
        channels: &channels,
        bearing_cell: DRIVE_IDS[9],
    };
    let body = fixture.body();
    let (native, before, torques, load, _) = steady_drive(1.0, 0.1);
    let drive = network.bind(&body, &native, before, || false).unwrap();
    let invalid_bearing = BatteryMotorThermalNetwork {
        cells: &cells,
        links: &[],
        channels: &channels,
        bearing_cell: CellId(999),
    };
    assert!(
        invalid_bearing
            .bind(&body, &native, before, || false)
            .is_err()
    );
    let changed_bearing = BatteryMotorThermalNetwork {
        cells: &cells,
        links: &[],
        channels: &channels,
        bearing_cell: DRIVE_IDS[8],
    };
    assert!(
        changed_bearing
            .step_candidate(
                &body,
                0.02,
                &drive,
                native.electrical().modulation(),
                &torques,
                load,
                electrical_controls(),
                drive_controls(),
                &[0.0; 10],
                || false,
            )
            .is_err()
    );
    let unsupported = RegulatedPowerLoad::try_new(1.3385, 10.0, 1.0, 5.0).unwrap();
    assert!(
        network
            .step_candidate(
                &body,
                0.02,
                &drive,
                native.electrical().modulation(),
                &torques,
                Some(unsupported),
                electrical_controls(),
                drive_controls(),
                &[0.0; 10],
                || false,
            )
            .is_err()
    );
    let polls = Cell::new(0);
    assert!(matches!(
        network.step_candidate(
            &body,
            0.02,
            &drive,
            native.electrical().modulation(),
            &torques,
            load,
            electrical_controls(),
            drive_controls(),
            &[0.0; 10],
            || {
                polls.set(polls.get() + 1);
                polls.get() > 20
            },
        ),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    assert_eq!(drive.state(), &before);
    assert_eq!(drive.time_s(), 0.0);
    assert_eq!(body.time_s(), 0.0);
    for id in DRIVE_IDS {
        near(body.temperature_k(id).unwrap(), 300.0);
    }
    let retry = network
        .step_candidate(
            &body,
            0.02,
            &drive,
            native.electrical().modulation(),
            &torques,
            load,
            electrical_controls(),
            drive_controls(),
            &[0.0; 10],
            || false,
        )
        .unwrap();
    let (advanced_body, advanced_drive, _) = retry.into_parts();
    assert_eq!(advanced_body.time_s(), 0.02);
    assert_eq!(advanced_drive.time_s(), 0.02);
    assert!(advanced_drive.state()[0] < before[0]);
    assert!(
        network
            .step_candidate(
                &advanced_body,
                0.04,
                &drive,
                native.electrical().modulation(),
                &torques,
                load,
                electrical_controls(),
                drive_controls(),
                &[0.0; 10],
                || false,
            )
            .is_err()
    );
}
