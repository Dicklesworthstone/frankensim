//! Full Q(V) source continuity, integral storage and finite material heat.
use super::*;
use fs_couple::circuit_enthalpy::converter::commutation::battery_motor::charge::*;
use fs_phs::circuit::capacitance::CapacitorChargeCurve;
use fs_phs::circuit::converter::commutation::capacitance::{
    battery_motor::PreparedChargeBatteryMotorCycle, charge::ChargeCommutationParameters,
};
use fs_qty::ElectricCharge;

const SCALE: [f64; 4] = [0.025, 0.04, 0.035, 0.05];
const VOLTAGE: [f64; 8] = [-100.0, -5.0, -0.5, 0.0, 0.5, 2.0, 20.0, 100.0];
const CHARGE: [f64; 8] = [-50.0, -5.0, -1.0, 0.0, 1.0, 2.5, 7.0, 15.0];
const CHORD_BUDGET: usize = 4096;

fn context(t: f64) -> QueryPoint {
    common::point(t)
        .with_quantity("Vds", QuantitySpec::dimensional(Voltage::DIMS), 0.0)
        .unwrap()
}
fn source_points(device: usize, tail: f64) -> Vec<(f64, f64)> {
    VOLTAGE
        .into_iter()
        .zip(CHARGE)
        .enumerate()
        .map(|(i, (v, q))| (v, SCALE[device] * q * if i == 7 { tail } else { 1.0 }))
        .collect()
}
fn curves(tail: f64) -> [CapacitorChargeCurve; 4] {
    core::array::from_fn(|i| {
        CapacitorChargeCurve::try_new(
            source_points(i, tail)
                .into_iter()
                .map(|(v, q)| (q, v))
                .collect(),
        )
        .unwrap()
    })
}
struct ChargeFixture {
    base: Fixture,
    keys: [PropertyKey; 4],
}
fn charge_fixture(
    cp: f64,
    threshold_slope: f64,
    upper_k: f64,
    temperature_law: bool,
    tail: f64,
) -> ChargeFixture {
    charge_fixture_reference(cp, threshold_slope, upper_k, temperature_law, tail, true)
}
fn charge_fixture_reference(
    cp: f64,
    threshold_slope: f64,
    upper_k: f64,
    temperature_law: bool,
    tail: f64,
    explicit_zero: bool,
) -> ChargeFixture {
    let f = fixture(cp, threshold_slope, 350.0, 100.0);
    let keys = core::array::from_fn(|i| {
        PropertyKey::new(format!("device_output_charge_{i}"), ElectricCharge::DIMS)
    });
    let mut claims = ClaimSet::new();
    for (_, claim) in f.card.claims().claims_ordered() {
        claims.insert_claim(claim.clone()).unwrap();
    }
    for (device, key) in keys.iter().enumerate() {
        let mut knots = source_points(device, tail);
        if !explicit_zero {
            knots.retain(|(v, _)| *v != 0.0);
        }
        let value = if temperature_law {
            let mut samples = Vec::new();
            for t in [299.0, upper_k] {
                samples.extend(
                    knots
                        .iter()
                        .map(|(v, q)| ([*v, t], q * (1.0 + 0.01 * (t - 300.0)))),
                );
            }
            let mut triangles = Vec::new();
            let count = knots.len() as u32;
            for i in 0..count - 1 {
                triangles.extend([[i, i + 1, i + count], [i + 1, i + count + 1, i + count]]);
            }
            PropertyValue::Surface2D {
                axes: ["Vds".into(), "T".into()],
                axis_dims: [Voltage::DIMS, Temperature::DIMS],
                samples,
                triangles,
                dims: ElectricCharge::DIMS,
            }
        } else {
            PropertyValue::Curve {
                abscissa: "Vds".into(),
                abscissa_dims: Voltage::DIMS,
                knots,
                dims: ElectricCharge::DIMS,
            }
        };
        claims
            .insert_claim(PropertyClaim {
                key: key.clone(),
                value,
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "T",
                        QuantitySpec::dimensional(Temperature::DIMS),
                        299.0,
                        upper_k,
                    )
                    .with_quantity(
                        "Vds",
                        QuantitySpec::dimensional(Voltage::DIMS),
                        -100.0,
                        100.0,
                    ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic monotone asymmetric Q(V) fixture, not qualified device data"
                        .into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic-device".into(),
            phase: "solid".into(),
            process: "authored fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &common::point(299.0),
            upper: &common::point(350.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-6,
    )
    .unwrap();
    ChargeFixture {
        base: Fixture {
            card,
            caloric,
            keys: f.keys,
            cp,
        },
        keys,
    }
}
impl ChargeFixture {
    fn bindings(&self) -> [DeviceChargeBinding<'_>; 4] {
        core::array::from_fn(|i| DeviceChargeBinding {
            cell: IDS[i + 1],
            key: &self.keys[i],
            voltage_axis: "Vds",
            voltage_range_v: [-100.0, 100.0],
            selection: MaterialPropertySelection::SingleClaimOnly,
        })
    }
}
fn coupled<'a>(
    cells: &'a [MaterialCell<'a>],
    resistors: &'a [ConverterThermalChannel<'a>; 5],
    devices: &'a [DeviceConductionBinding<'a>; 4],
    charges: &'a [DeviceChargeBinding<'a>; 4],
) -> ChargeBatteryMotorThermalNetwork<'a> {
    ChargeBatteryMotorThermalNetwork {
        drive: CommutatedBatteryMotorThermalNetwork {
            cells,
            links: &[],
            resistors,
            devices,
            bearing_cell: IDS[0],
        },
        charges,
        maximum_cells_per_chord: CHORD_BUDGET,
    }
}
fn initial(steady: &SteadyDrive, node_offset_v: [f64; 2]) -> [f64; 11] {
    let p = parameters();
    let current = steady.state[3] / p.inductance_h;
    let input = steady.state[4] / p.bus_capacitance_f[0];
    let output = steady.state[5] / p.bus_capacitance_f[1];
    let left = if steady.gates.on[0] {
        (10.0 * input - current) / 10.01
    } else {
        (5.01 * input + 2.5 - current) / 5.02
    };
    let right = left - 0.1 * current;
    let c = ChargeCommutationParameters {
        converter: p,
        output_charge_curves: curves(1.0),
        maximum_cells_per_chord: CHORD_BUDGET,
    }
    .state(
        current,
        [input, output],
        [left + node_offset_v[0], right + node_offset_v[1]],
    )
    .unwrap();
    let mut state = [0.0; 11];
    state[..3].copy_from_slice(&steady.state[..3]);
    state[3..8].copy_from_slice(&c);
    state[8..].copy_from_slice(&steady.state[6..]);
    state
}
// Independent source-table evaluation and exact per-segment integral V dQ.
// Neither QV/2 nor a pointwise differential capacitance defines this energy.
fn q(device: usize, v: f64) -> f64 {
    let points = source_points(device, 1.0);
    let pair = points
        .windows(2)
        .find(|p| p[0].0 <= v && v <= p[1].0)
        .unwrap();
    pair[0].1 + (v - pair[0].0) / (pair[1].0 - pair[0].0) * (pair[1].1 - pair[0].1)
}
fn energy(device: usize, v: f64) -> f64 {
    let (lo, hi) = (v.min(0.0), v.max(0.0));
    let mut integral = 0.0;
    for pair in source_points(device, 1.0).windows(2) {
        let a = lo.max(pair[0].0);
        let b = hi.min(pair[1].0);
        if a < b {
            integral += 0.5 * (a + b) * (q(device, b) - q(device, a));
        }
    }
    if v < 0.0 { -integral } else { integral }
}
fn voltage(state: &[f64; 11], leg: usize) -> [f64; 2] {
    // Eliminate the bare bus capacitor, then independently bisect monotone
    // Ql(Vn)-Qh((qb+qn-Ql(Vn))/Cb-Vn)=qn on this fixture's known support.
    let cb = parameters().bus_capacitance_f[leg];
    let (qb, qn) = (state[4 + leg], state[6 + leg]);
    let high = 2 * leg;
    let low = high + 1;
    let residual = |node| {
        let bus = (qb + qn - q(low, node)) / cb;
        q(low, node) - q(high, bus - node) - qn
    };
    let (mut lo, mut hi) = (-20.0, 20.0);
    assert!(residual(lo) < 0.0 && residual(hi) > 0.0);
    for _ in 0..80 {
        let mid = f64::midpoint(lo, hi);
        if residual(mid) > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    let node = f64::midpoint(lo, hi);
    [(qb + qn - q(low, node)) / cb, node]
}
fn device_voltages(state: &[f64; 11]) -> [f64; 4] {
    core::array::from_fn(|i| {
        let [bus, node] = voltage(state, i / 2);
        if i % 2 == 0 { bus - node } else { node }
    })
}

#[test]
fn nonlinear_charge_sources_preserve_independent_steady_motoring_and_regenerative_heat() {
    for regeneration in [false, true] {
        let f = charge_fixture(1.0, 0.0, 350.0, false, 1.0);
        let thermal = f.base.thermal();
        let context = context(300.0);
        let cells = f.base.cells(&context);
        let r = drive_resistors(&f.base);
        let d = f.base.devices();
        let charges = f.bindings();
        let network = coupled(&cells, &r, &d, &charges);
        let steady = steady_drive(regeneration);
        let state = initial(&steady, [0.0; 2]);
        let authored = curves(1.0);
        let drive = network
            .bind(&thermal, &steady.native, &authored, state, || false)
            .unwrap();
        let dt = 1.0 / 64.0;
        let next = network
            .step_candidate(
                &thermal,
                dt,
                &drive,
                &[CommutationInterval {
                    gates: steady.gates,
                    duration_s: dt,
                }],
                &[steady.command],
                4,
                Some(IDS[3]),
                power_controls(),
                cycle_controls(),
                drive_controls(7),
                &[0.0; 7],
                || false,
            )
            .unwrap();
        let step = &next.report().electrical;
        let heat = steady.heat_w.map(|p| dt * p);
        let load = steady.command.auxiliary_load.unwrap().requested_power_w() * dt;
        for (a, b) in step.heat_joules().into_iter().zip(heat) {
            near(a, b);
        }
        for value in step.converter.device_capacitor_energy_change_j {
            near(value, 0.0);
        }
        for i in 1..11 {
            near(next.drive().state()[i], state[i]);
        }
        near(
            next.drive().state()[0],
            state[0] - dt * steady.cell_current_a,
        );
        near(step.storage_change_j[0], -12.0 * dt * steady.cell_current_a);
        near(step.load_work_j, load);
        for (i, expected) in per_body(heat, load).into_iter().enumerate() {
            near(
                next.thermal().temperature_k(IDS[i]).unwrap() - 300.0,
                expected,
            );
        }
        near(
            next.report().deposited_joule_j,
            heat.iter().sum::<f64>() + load,
        );
        near(next.report().total_energy_residual_j, 0.0);
        for i in 0..4 {
            let accepted = &next.drive().charges()[i];
            assert_eq!(accepted.curve(), &authored[i]);
            assert_eq!(
                drive.charges()[i].source().card_identity(),
                accepted.source().card_identity()
            );
            near(
                accepted
                    .source()
                    .query_point()
                    .iter()
                    .find(|(axis, _)| axis == "T")
                    .unwrap()
                    .1,
                next.thermal().temperature_k(IDS[i + 1]).unwrap(),
            );
            assert_eq!(accepted.voltage_range_v(), [-100.0, 100.0]);
            assert_eq!(accepted.source_path().knots.len(), VOLTAGE.len());
            assert_eq!(&next.report().source.charges[i], accepted);
        }
        assert_eq!(next.report().source.maximum_cells_per_chord, CHORD_BUDGET);
    }
}

#[test]
fn integral_device_energy_remains_storage_while_hot_diode_heat_and_all_electrode_flows_close() {
    let f = charge_fixture(0.01, 0.02, 350.0, false, 1.0);
    let thermal = f.base.thermal();
    let context = context(300.0);
    let cells = f.base.cells(&context);
    let r = drive_resistors(&f.base);
    let d = f.base.devices();
    let charges = f.bindings();
    let network = coupled(&cells, &r, &d, &charges);
    let steady = steady_drive(false);
    let state = initial(&steady, [-0.45, 0.4]);
    let authored = curves(1.0);
    let drive = network
        .bind(&thermal, &steady.native, &authored, state, || false)
        .unwrap();
    let dt = 1.0 / 1024.0;
    let next = network
        .step_candidate(
            &thermal,
            dt,
            &drive,
            &[CommutationInterval {
                gates: steady.gates,
                duration_s: dt,
            }],
            &[steady.command],
            1,
            None,
            power_controls(),
            cycle_controls(),
            drive_controls(7),
            &[0.0; 7],
            || false,
        )
        .unwrap();
    let step = &next.report().electrical;
    let end = next.drive().state();
    let mid: [f64; 11] = core::array::from_fn(|i| f64::midpoint(state[i], end[i]));
    let p = next.drive().parameters();
    let b = next.drive().battery_parameters();
    let m = next.drive().motor_parameters();
    // The nonlinear native owner reports its actual accepted discrete effort.
    // Single-step KCL below independently ties it to every physical charge update.
    let [input, output] = step.converter.mean_bus_voltage_v;
    let [left, right] = step.converter.mean_node_voltage_v;
    let vds = [input - left, left, output - right, right];
    let current = mid[3] / p.inductance_h;
    let pol = mid[1] / b.polarization_capacitance_f;
    let battery_bus = mid[2] / b.bus_capacitance_f;
    let cell_current = (12.0 - pol - battery_bus) / b.series_resistance_ohm;
    let lead_current = (battery_bus - input) / next.drive().lead_resistance_ohm();
    let motor_current = mid[8] / m.inductance_h;
    let omega = [
        mid[9] / m.rotor_inertia_kg_m2,
        mid[10] / m.stator_inertia_kg_m2,
    ];
    let relative = omega[0] - omega[1];
    let mut heat = [0.0; 14];
    let mut branch_current = [0.0; 4];
    heat[0] = dt * p.winding_resistance_ohm * current * current;
    for i in 0..4 {
        let device = p.devices[i];
        let temperature = next.thermal().temperature_k(IDS[i + 1]).unwrap();
        near(
            device.body_diode.forward_voltage_v,
            0.5 + 0.02 * (temperature - 300.0),
        );
        let diode = ((-vds[i] - device.body_diode.forward_voltage_v)
            / device.body_diode.differential_resistance_ohm)
            .max(0.0);
        let r = if steady.gates.on[i] {
            device.on_resistance_ohm
        } else {
            device.off_resistance_ohm
        };
        branch_current[i] = vds[i] / r - diode;
        heat[1 + i] = dt * vds[i] * vds[i] / r;
        heat[10 + i] = dt * -vds[i] * diode;
        near(step.converter.last_devices[i].diode_current_a, diode);
    }
    heat[5] = dt * b.series_resistance_ohm * cell_current * cell_current;
    heat[6] = dt * pol * pol / b.polarization_resistance_ohm;
    heat[7] = dt * next.drive().lead_resistance_ohm() * lead_current * lead_current;
    heat[8] = dt * m.winding_resistance_ohm * motor_current * motor_current;
    heat[9] = dt * m.bearing_drag_nm_s_per_rad * relative * relative;
    for (a, b) in step.heat_joules().into_iter().zip(heat) {
        near(a, b);
    }
    let aux = steady.command.auxiliary_load.unwrap().requested_power_w() / output;
    let derivatives = [
        -cell_current,
        cell_current - pol / b.polarization_resistance_ohm,
        cell_current - lead_current,
        left - right - p.winding_resistance_ohm * current,
        lead_current - branch_current[0],
        -motor_current - aux - branch_current[2],
        branch_current[0] - branch_current[1] - current,
        branch_current[2] - branch_current[3] + current,
        output - m.winding_resistance_ohm * motor_current - m.coupling_nm_per_a * relative,
        m.coupling_nm_per_a * motor_current - m.bearing_drag_nm_s_per_rad * relative
            + steady.command.shaft_torque_nm[0],
        -m.coupling_nm_per_a * motor_current
            + m.bearing_drag_nm_s_per_rad * relative
            + steady.command.shaft_torque_nm[1],
    ];
    for i in 0..11 {
        near(end[i] - state[i], dt * derivatives[i]);
    }

    let before_voltage = device_voltages(&state);
    let after_voltage = device_voltages(end);
    assert!(
        before_voltage[0] > 0.5 && after_voltage[0] < 0.5,
        "the input device must traverse an original charge-law knot"
    );
    let changes: [f64; 4] =
        core::array::from_fn(|i| energy(i, after_voltage[i]) - energy(i, before_voltage[i]));
    let incorrect_half_qv: [f64; 4] = core::array::from_fn(|i| {
        0.5 * (q(i, after_voltage[i]) * after_voltage[i]
            - q(i, before_voltage[i]) * before_voltage[i])
    });
    for (a, b) in step
        .converter
        .device_capacitor_energy_change_j
        .into_iter()
        .zip(changes)
    {
        near(a, b);
    }
    assert!(changes.iter().map(|v| v.abs()).sum::<f64>() > 1e-4);
    assert!(
        changes
            .iter()
            .zip(incorrect_half_qv)
            .map(|(a, b)| (a - b).abs())
            .sum::<f64>()
            > 1e-5,
        "the fixture must distinguish the original integral from QV/2"
    );
    for (i, expected) in per_body(heat, 0.0).into_iter().enumerate() {
        near(
            f.base.cp * (next.thermal().temperature_k(IDS[i]).unwrap() - 300.0),
            expected,
        );
        if i < 4 {
            assert_eq!(next.drive().charges()[i].curve(), &authored[i]);
        }
    }
    near(next.report().deposited_joule_j, heat.iter().sum());
    near(next.report().total_energy_residual_j, 0.0);
    assert!(
        next.drive().devices()[2]
            .parameters()
            .body_diode
            .forward_voltage_v
            > 0.500001
    );

    let other = DeviceGateState {
        on: [false, false, true, false],
    };
    let continued = network
        .step_candidate(
            next.thermal(),
            2.0 * dt,
            next.drive(),
            &[CommutationInterval {
                gates: other,
                duration_s: dt,
            }],
            &[steady.command],
            1,
            None,
            power_controls(),
            cycle_controls(),
            drive_controls(7),
            &[0.0; 7],
            || false,
        )
        .unwrap();
    let accepted = continued.drive();
    let original = replay_original(&steady.native, accepted);
    let mut replay = PreparedChargeBatteryMotorCycle::new_charge_controlled(
        original,
        accepted.parameters().devices,
        authored.clone(),
        CHORD_BUDGET,
        &[CommutationInterval {
            gates: other,
            duration_s: dt,
        }],
        1,
        power_controls(),
        cycle_controls(),
        || false,
    )
    .unwrap();
    let mut replay_end = [0.0; 11];
    let expected = replay
        .step_controlled(
            next.drive().state(),
            steady.command,
            &mut replay_end,
            || false,
        )
        .unwrap();
    for (a, b) in accepted.state().iter().zip(replay_end) {
        near(*a, b);
    }
    for (a, b) in continued
        .report()
        .electrical
        .heat_joules()
        .into_iter()
        .zip(expected.heat_joules())
    {
        near(a, b);
    }
    assert!((accepted.state()[6] - next.drive().state()[6]).abs() > 1e-8);
}
fn replay_original(
    original: &BatteryMotorDrive,
    accepted: &MaterialChargeBatteryMotorDrive,
) -> BatteryMotorDrive {
    let p = accepted.parameters();
    AveragedBuckBoost::new(
        BuckBoostParameters {
            inductance_h: p.inductance_h,
            bus_capacitance_f: p.bus_capacitance_f,
            winding_resistance_ohm: p.winding_resistance_ohm,
            input_on_resistance_ohm: [
                p.devices[0].on_resistance_ohm,
                p.devices[1].on_resistance_ohm,
            ],
            output_on_resistance_ohm: [
                p.devices[2].on_resistance_ohm,
                p.devices[3].on_resistance_ohm,
            ],
        },
        BridgeModulation {
            input_high_fraction: 0.0,
            output_high_fraction: 0.0,
        },
    )
    .unwrap()
    .with_battery_source(
        BatterySource::new(
            original.electrical().battery_curve().clone(),
            accepted.battery_parameters(),
        )
        .unwrap(),
        accepted.lead_resistance_ohm(),
    )
    .unwrap()
    .with_dc_motor(DcMotor::new(accepted.motor_parameters()).unwrap())
    .unwrap()
}

#[test]
fn full_original_law_support_and_temperature_authority_are_required_for_charge_history() {
    let steady = steady_drive(false);
    let state = initial(&steady, [0.0; 2]);
    let authored = curves(1.0);
    for temperature_law in [false, true] {
        let f = charge_fixture(1.0, 0.0, 350.0, temperature_law, 1.0);
        let thermal = f.base.thermal();
        let context = context(300.0);
        let cells = f.base.cells(&context);
        let r = drive_resistors(&f.base);
        let d = f.base.devices();
        let bindings = f.bindings();
        let network = coupled(&cells, &r, &d, &bindings);
        let bound = network.bind(&thermal, &steady.native, &authored, state, || false);
        if temperature_law {
            assert!(bound.is_err());
            continue;
        }
        let drive = bound.unwrap();

        // The modified tail has identical Q(0), the same local slope and the
        // same charge throughout the actual initial voltage range.
        let altered = curves(1.1);
        for i in 0..4 {
            near(
                authored[i].charge_c(0.0).unwrap(),
                altered[i].charge_c(0.0).unwrap(),
            );
            near(
                authored[i].charge_c(12.0).unwrap(),
                altered[i].charge_c(12.0).unwrap(),
            );
        }
        assert!(
            network
                .bind(&thermal, &steady.native, &altered, state, || false)
                .is_err()
        );
        let missing = PropertyKey::new("undeclared_output_charge", ElectricCharge::DIMS);
        for key in [&missing, &f.base.keys[0]] {
            let mut bad = f.bindings();
            bad[0].key = key;
            assert!(
                coupled(&cells, &r, &d, &bad)
                    .bind(&thermal, &steady.native, &authored, state, || false)
                    .is_err()
            );
        }
        let mut wrong_package = f.bindings();
        wrong_package[0].cell = IDS[0];
        assert!(
            coupled(&cells, &r, &d, &wrong_package)
                .bind(&thermal, &steady.native, &authored, state, || false)
                .is_err()
        );
        let untyped = common::point(300.0).with("Vds", 0.0).unwrap();
        let untyped_cells = f.base.cells(&untyped);
        assert!(
            coupled(&untyped_cells, &r, &d, &bindings)
                .bind(&thermal, &steady.native, &authored, state, || false)
                .is_err()
        );
        let mut zero_budget = coupled(&cells, &r, &d, &bindings);
        zero_budget.maximum_cells_per_chord = 0;
        assert!(
            zero_budget
                .bind(&thermal, &steady.native, &authored, state, || false)
                .is_err()
        );

        let dt = 1.0 / 64.0;
        let mut restricted = f.bindings();
        restricted[0].voltage_range_v[0] = -99.0;
        assert!(
            coupled(&cells, &r, &d, &restricted)
                .step_candidate(
                    &thermal,
                    dt,
                    &drive,
                    &[CommutationInterval {
                        gates: steady.gates,
                        duration_s: dt
                    }],
                    &[steady.command],
                    4,
                    None,
                    power_controls(),
                    cycle_controls(),
                    drive_controls(7),
                    &[0.0; 7],
                    || false,
                )
                .is_err()
        );
        let changed = charge_fixture(1.0, 0.0, 350.0, false, 1.1);
        let changed_cells = changed.base.cells(&context);
        let changed_r = drive_resistors(&changed.base);
        let changed_d = changed.base.devices();
        let changed_bindings = changed.bindings();
        assert!(
            coupled(&changed_cells, &changed_r, &changed_d, &changed_bindings)
                .step_candidate(
                    &thermal,
                    dt,
                    &drive,
                    &[CommutationInterval {
                        gates: steady.gates,
                        duration_s: dt
                    }],
                    &[steady.command],
                    4,
                    None,
                    power_controls(),
                    cycle_controls(),
                    drive_controls(7),
                    &[0.0; 7],
                    || false,
                )
                .is_err()
        );
        assert_eq!(drive.state(), &state);
        assert_eq!(drive.time_s(), 0.0);
    }

    // A source segment already passing through (0,0) may acquire that exact
    // canonical knot without changing the original law or its integral.
    let f = charge_fixture_reference(1.0, 0.0, 350.0, false, 1.0, false);
    let thermal = f.base.thermal();
    let reference = context(300.0);
    let cells = f.base.cells(&reference);
    let r = drive_resistors(&f.base);
    let d = f.base.devices();
    let bindings = f.bindings();
    let drive = coupled(&cells, &r, &d, &bindings)
        .bind(&thermal, &steady.native, &authored, state, || false)
        .unwrap();
    for i in 0..4 {
        assert_eq!(drive.charges()[i].curve(), &authored[i]);
        assert_eq!(drive.charges()[i].source_path().knots.len(), VOLTAGE.len());
    }

    let f = charge_fixture(1.0, 0.0, 300.000001, false, 1.0);
    let thermal = f.base.thermal();
    let reference = context(300.0);
    let cells = f.base.cells(&reference);
    let r = drive_resistors(&f.base);
    let d = f.base.devices();
    let bindings = f.bindings();
    let network = coupled(&cells, &r, &d, &bindings);
    let drive = network
        .bind(&thermal, &steady.native, &authored, state, || false)
        .unwrap();
    let before = thermal.cells().clone();
    let dt = 1.0 / 64.0;
    assert!(
        network
            .step_candidate(
                &thermal,
                dt,
                &drive,
                &[CommutationInterval {
                    gates: steady.gates,
                    duration_s: dt
                }],
                &[steady.command],
                4,
                None,
                power_controls(),
                cycle_controls(),
                drive_controls(7),
                &[0.0; 7],
                || false,
            )
            .is_err()
    );
    assert_eq!(thermal.cells(), &before);
    assert_eq!(drive.state(), &state);
}

#[test]
fn nonlinear_charge_drive_uses_the_original_two_inventory_heat_assembly() {
    let f = charge_fixture(1.0, 0.0, 350.0, false, 1.0);
    let base = f.base.thermal();
    let mut events = base
        .cells()
        .values()
        .map(|cell| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell.clone()),
        })
        .collect::<Vec<_>>();
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 77,
            a: IDS[0],
            b: IDS[6],
            area_m2: 1.0,
            conductance_w_k: 0.0,
        }),
    });
    let thermal = DepositionSimulation::new(
        events,
        &[(MaterialId(0), f.base.caloric.curve() as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 7,
            max_links: 1,
            max_events: 8,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.1,
            energy_tolerance_j: 1e-9,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    let reference = context(300.0);
    let cells = f.base.cells(&reference);
    let bulk = MaterialNetwork {
        cells: &cells,
        links: &[],
        resistors: &[],
    };
    let electrical = ElectricalContactNetwork {
        bulk: &bulk,
        contacts: &[],
        electrical_bases: &[],
    };
    let heat = MaterialHeatNetwork {
        electrical: &electrical,
        thermal_contacts: &[],
        radiation: None,
    };
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 41,
                kernel_id: 93,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        let prepared = heat
            .prepare_with_link_conductances(&cx, &thermal, &BTreeMap::from([(77, 0.0)]))
            .unwrap();
        let transport = FixedCooling {
            prepared: &prepared,
            conductance_w_k: 0.25,
        };
        let r = drive_resistors(&f.base);
        let d = f.base.devices();
        let bindings = f.bindings();
        let network = ChargeBatteryMotorHeatNetwork {
            drive: CommutatedBatteryMotorHeatNetwork {
                assembly: &transport,
                resistors: &r,
                devices: &d,
                bearing_cell: IDS[0],
            },
            charges: &bindings,
            maximum_cells_per_chord: CHORD_BUDGET,
        };
        let steady = steady_drive(false);
        let state = initial(&steady, [0.0; 2]);
        let authored = curves(1.0);
        let drive = network
            .bind(&cx, &thermal, &steady.native, &authored, state)
            .unwrap();
        let dt = 1.0 / 64.0;
        let next = network
            .step_candidate(
                &cx,
                &thermal,
                dt,
                &drive,
                &[CommutationInterval {
                    gates: steady.gates,
                    duration_s: dt,
                }],
                &[steady.command],
                4,
                None,
                power_controls(),
                cycle_controls(),
                drive_controls(7),
                &[0.0; 7],
            )
            .unwrap();
        let expected = per_body(steady.heat_w.map(|p| dt * p), 0.0);
        let exchange = dt * 0.25 * (expected[0] - expected[6]) / (1.0 + 2.0 * dt * 0.25);
        near(
            next.thermal().temperature_k(IDS[0]).unwrap(),
            300.0 + expected[0] - exchange,
        );
        near(
            next.thermal().temperature_k(IDS[6]).unwrap(),
            300.0 + expected[6] + exchange,
        );
        for i in 1..6 {
            near(
                next.thermal().temperature_k(IDS[i]).unwrap(),
                300.0 + expected[i],
            );
        }
        near(
            next.thermal().energy_j() - thermal.energy_j(),
            expected.iter().sum(),
        );
        near(next.report().deposited_joule_j, expected.iter().sum());
        near(next.report().total_energy_residual_j, 0.0);
        assert_eq!(next.report().source.drive.network.material.cells.len(), 7);
        for i in 0..4 {
            assert_eq!(next.report().source.charges[i].curve(), &authored[i]);
        }
        for value in next
            .report()
            .electrical
            .converter
            .device_capacitor_energy_change_j
        {
            near(value, 0.0);
        }
    });
}

// Charge laws bound under SingleClaimOnly are refreshed under a pin of the same claim by the first
// step, which must be accepted; different charge laws are different sources. Both source forms
// step: an authored 0 V knot, and a law whose 0 V knot the library inserts into the path.
#[test]
fn charge_laws_bound_under_single_claim_only_take_their_first_step_and_differ_from_each_other() {
    for explicit_zero in [true, false] {
        let f = charge_fixture_reference(1.0, 0.0, 350.0, false, 1.0, explicit_zero);
        let thermal = f.base.thermal();
        let ctx = context(300.0);
        let cells = f.base.cells(&ctx);
        let r = drive_resistors(&f.base);
        let d = f.base.devices();
        let charges = f.bindings();
        let network = coupled(&cells, &r, &d, &charges);
        let steady = steady_drive(false);
        let state = initial(&steady, [0.0; 2]);
        let authored = curves(1.0);
        let drive = network
            .bind(&thermal, &steady.native, &authored, state, || false)
            .unwrap();
        for charge in drive.charges() {
            assert_eq!(charge.source_path().knots.len(), VOLTAGE.len());
        }
        let dt = 1.0 / 64.0;
        match network.step_candidate(
            &thermal,
            dt,
            &drive,
            &[CommutationInterval {
                gates: steady.gates,
                duration_s: dt,
            }],
            &[steady.command],
            4,
            Some(IDS[3]),
            power_controls(),
            cycle_controls(),
            drive_controls(7),
            &[0.0; 7],
            || false,
        ) {
            Ok(_) => {}
            Err(error) => panic!(
                "a step from the bound state is refused (explicit 0 V knot: {explicit_zero}): {error:?}"
            ),
        }
        // (Sources refreshed by the step sit at the stepped state, so their query point differs
        // from the bound one by design; the accepted step above is the positive.)
        // Four charge laws of four packages: another key, cell and curve is another source.
        assert_ne!(drive.charges()[0], drive.charges()[1]);
    }
}

// Planted negative on the same check (see the converter's): the drive stays bound at 300 K while the
// accepted heat state is another one at the same time, and the diode threshold depends on temperature.
#[test]
fn charge_drive_bound_at_one_state_refuses_a_step_from_another_accepted_state() {
    let f = charge_fixture(1.0, 0.02, 350.0, false, 1.0);
    let thermal = f.base.thermal();
    let hot = f.base.thermal_at(320.0);
    assert!((hot.temperature_k(IDS[1]).unwrap() - 320.0).abs() < 1e-5);
    let ctx = context(300.0);
    let cells = f.base.cells(&ctx);
    let r = drive_resistors(&f.base);
    let d = f.base.devices();
    let charges = f.bindings();
    let network = coupled(&cells, &r, &d, &charges);
    let steady = steady_drive(false);
    let state = initial(&steady, [0.0; 2]);
    let authored = curves(1.0);
    let drive = network
        .bind(&thermal, &steady.native, &authored, state, || false)
        .unwrap();
    let dt = 1.0 / 64.0;
    let attempt = |at: &DepositionSimulation<'_>| {
        network
            .step_candidate(
                at,
                dt,
                &drive,
                &[CommutationInterval {
                    gates: steady.gates,
                    duration_s: dt,
                }],
                &[steady.command],
                4,
                Some(IDS[3]),
                power_controls(),
                cycle_controls(),
                drive_controls(7),
                &[0.0; 7],
                || false,
            )
            .map(|_| ())
            .map_err(|error| format!("{error:?}"))
    };
    assert_eq!(attempt(&thermal), Ok(()));
    let Err(text) = attempt(&hot) else {
        panic!("a step from another accepted state was taken");
    };
    assert!(
        text.contains(
            "commutated drive must retain every accepted resistance, conduction and capacitance source"
        ),
        "{text}"
    );
}
