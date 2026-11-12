//! Original node-charge storage, semiconductor calorimetry and source continuity.
use super::*;
use fs_couple::circuit_enthalpy::converter::commutation::battery_motor::capacitance::*;
use fs_phs::circuit::converter::commutation::capacitance::{
    CapacitiveCommutationParameters,
    battery_motor::PreparedCapacitiveBatteryMotorCycle,
};
use fs_qty::Capacitance;

const COSS: [f64; 4] = [0.025, 0.04, 0.035, 0.05];

struct CapacitiveFixture { base: Fixture, keys: [PropertyKey; 4] }
fn capacitive_fixture(cp: f64, threshold_slope: f64, cap_upper_k: f64, variable: bool, scale: f64) -> CapacitiveFixture {
    let f = fixture(cp, threshold_slope, 350.0, 100.0);
    let keys = core::array::from_fn(|i| PropertyKey::new(format!("device_output_capacitance_{i}"), Capacitance::DIMS));
    let mut claims = ClaimSet::new();
    for (_, claim) in f.card.claims().claims_ordered() { claims.insert_claim(claim.clone()).unwrap(); }
    for (i, key) in keys.iter().enumerate() {
        let value = COSS[i] * scale;
        claims.insert_claim(PropertyClaim {
            key: key.clone(),
            value: if variable {
                // Even an accidentally flat T curve is not an original constant charge law.
                PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                    knots: vec![(299.0, value), (cap_upper_k, value)], dims: Capacitance::DIMS }
            } else { PropertyValue::Scalar { value, dims: Capacitance::DIMS } },
            validity: ValidityDomain::unconstrained()
                .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 299.0, cap_upper_k),
            uncertainty: UncertaintyModel::Unstated,
            // The Scalar is a constant law: LinearInside has no knot span on a scalar.
            interpolation: if variable { InterpolationPolicy::LinearInside } else { InterpolationPolicy::ConstantWithinValidity },
            observations: vec![],
            provenance: Provenance {
                source: "synthetic constant device capacitor fixture, not qualified device data".into(),
                license: "CC0-1.0".into(), artifact: None,
            },
        }).unwrap();
    }
    let card = MaterialCard::assemble(MaterialStateId {
        chemistry: "synthetic-device".into(), phase: "solid".into(),
        process: "authored fixture".into(), revision: 0,
    }, claims, vec![]).unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(
        &card, HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS), temperature_axis: "T",
            lower: &common::point(299.0), upper: &common::point(350.0),
            selection: SelectionPolicy::SingleClaimOnly,
        }, SolidLiquidPhase::Solid, 0.0, 1e-6,
    ).unwrap();
    CapacitiveFixture { base: Fixture { card, caloric, keys: f.keys, cp }, keys }
}
impl CapacitiveFixture {
    fn bindings(&self) -> [DeviceCapacitanceBinding<'_>; 4] {
        core::array::from_fn(|i| DeviceCapacitanceBinding {
            cell: IDS[i + 1], key: &self.keys[i], selection: MaterialPropertySelection::SingleClaimOnly,
        })
    }
}
fn coupled<'a>(
    cells: &'a [MaterialCell<'a>], resistors: &'a [ConverterThermalChannel<'a>; 5],
    devices: &'a [DeviceConductionBinding<'a>; 4], capacitances: &'a [DeviceCapacitanceBinding<'a>; 4],
) -> CapacitiveBatteryMotorThermalNetwork<'a> {
    CapacitiveBatteryMotorThermalNetwork {
        drive: CommutatedBatteryMotorThermalNetwork { cells, links: &[], resistors, devices, bearing_cell: IDS[0] },
        capacitances,
    }
}
fn initial(steady: &SteadyDrive, node_offset_v: [f64; 2]) -> [f64; 11] {
    let p = parameters();
    let current = steady.state[3] / p.inductance_h;
    let input = steady.state[4] / p.bus_capacitance_f[0];
    let output = steady.state[5] / p.bus_capacitance_f[1];
    // Independent DC branch KCL from the steady fixture; only this explicitly
    // authored initialization is steady. Subsequent nodes are dynamic charges.
    let left = if steady.gates.on[0] { (10.0 * input - current) / 10.01 }
        else { (5.01 * input + 2.5 - current) / 5.02 };
    let right = left - 0.1 * current;
    let c = CapacitiveCommutationParameters { converter: p, output_capacitance_f: COSS }
        .state(current, [input, output], [left + node_offset_v[0], right + node_offset_v[1]]).unwrap();
    let mut state = [0.0; 11];
    state[..3].copy_from_slice(&steady.state[..3]);
    state[3..8].copy_from_slice(&c);
    state[8..].copy_from_slice(&steady.state[6..]);
    state
}
// Independent two-electrode charge equations for each bus/node pair.
// q_b=(Cb+Ch)*Vb-Ch*Vn; q_n=-Ch*Vb+(Ch+Cl)*Vn.
fn voltage(state: &[f64; 11], leg: usize) -> [f64; 2] {
    let cb = parameters().bus_capacitance_f[leg];
    let (ch, cl) = (COSS[2 * leg], COSS[2 * leg + 1]);
    let det = cb * (ch + cl) + ch * cl;
    let (qb, qn) = (state[4 + leg], state[6 + leg]);
    [((ch + cl) * qb + ch * qn) / det, (ch * qb + (cb + ch) * qn) / det]
}
fn capacitor_energy(state: &[f64; 11]) -> [f64; 4] {
    core::array::from_fn(|i| {
        let [bus, node] = voltage(state, i / 2);
        let v = if i % 2 == 0 { bus - node } else { node };
        0.5 * COSS[i] * v * v
    })
}

#[test]
fn constant_capacitors_preserve_independent_steady_motoring_and_regenerative_calorimetry() {
    for regeneration in [false, true] {
        let f = capacitive_fixture(1.0, 0.0, 350.0, false, 1.0);
        let thermal = f.base.thermal();
        let context = common::point(300.0);
        let cells = f.base.cells(&context);
        let resistors = drive_resistors(&f.base);
        let devices = f.base.devices();
        let caps = f.bindings();
        let network = coupled(&cells, &resistors, &devices, &caps);
        let steady = steady_drive(regeneration);
        let state = initial(&steady, [0.0; 2]);
        let drive = network.bind(&thermal, &steady.native, COSS, state, || false).unwrap();
        let dt = 1.0 / 64.0;
        let next = network.step_candidate(
            &thermal, dt, &drive,
            &[CommutationInterval { gates: steady.gates, duration_s: dt }], &[steady.command],
            4, Some(IDS[3]), power_controls(), cycle_controls(), drive_controls(7), &[0.0; 7], || false,
        ).unwrap();
        let step = &next.report().electrical;
        let heat = steady.heat_w.map(|p| dt * p);
        let load = steady.command.auxiliary_load.unwrap().requested_power_w() * dt;
        for (a, b) in step.heat_joules().into_iter().zip(heat) { near(a, b); }
        for q in step.converter.device_capacitor_energy_change_j { near(q, 0.0); }
        for i in 1..11 { near(next.drive().state()[i], state[i]); }
        near(next.drive().state()[0], state[0] - dt * steady.cell_current_a);
        near(step.storage_change_j[0], -12.0 * dt * steady.cell_current_a);
        near(step.load_work_j, load);
        for (i, q) in per_body(heat, load).into_iter().enumerate() {
            near(next.thermal().temperature_k(IDS[i]).unwrap() - 300.0, q);
        }
        near(next.report().deposited_joule_j, heat.iter().sum::<f64>() + load);
        near(next.report().total_energy_residual_j, 0.0);
        for i in 0..4 {
            let old = &drive.capacitances()[i];
            let new = &next.drive().capacitances()[i];
            assert_eq!(new.capacitance_f().to_bits(), COSS[i].to_bits());
            assert_eq!(old.source().card_identity(), new.source().card_identity());
            near(new.source().query_point().iter().find(|(axis, _)| axis == "T").unwrap().1, next.thermal().temperature_k(IDS[i + 1]).unwrap());
            assert_eq!(next.report().source.capacitances[i], *new);
        }
    }
}

#[test]
fn dynamic_capacitor_storage_is_separate_from_hot_device_heat_and_nodes_continue_across_gates() {
    let f = capacitive_fixture(0.01, 0.02, 350.0, false, 1.0);
    let thermal = f.base.thermal();
    let context = common::point(300.0);
    let cells = f.base.cells(&context);
    let resistors = drive_resistors(&f.base);
    let devices = f.base.devices();
    let caps = f.bindings();
    let network = coupled(&cells, &resistors, &devices, &caps);
    let steady = steady_drive(false);
    let state = initial(&steady, [-0.3, 0.4]);
    let drive = network.bind(&thermal, &steady.native, COSS, state, || false).unwrap();
    let dt = 1.0 / 1024.0;
    let gates = steady.gates;
    let next = network.step_candidate(
        &thermal, dt, &drive, &[CommutationInterval { gates, duration_s: dt }], &[steady.command],
        1, None, power_controls(), cycle_controls(), drive_controls(7), &[0.0; 7], || false,
    ).unwrap();
    let step = &next.report().electrical;
    let end = next.drive().state();
    let mid = core::array::from_fn(|i| f64::midpoint(state[i], end[i]));
    let [input, left] = voltage(&mid, 0);
    let [output, right] = voltage(&mid, 1);
    let p = next.drive().parameters();
    let b = next.drive().battery_parameters();
    let m = next.drive().motor_parameters();
    let vds = [input - left, left, output - right, right];
    let current = mid[3] / p.inductance_h;
    let pol = mid[1] / b.polarization_capacitance_f;
    let battery_bus = mid[2] / b.bus_capacitance_f;
    let cell_current = (12.0 - pol - battery_bus) / b.series_resistance_ohm;
    let lead_current = (battery_bus - input) / next.drive().lead_resistance_ohm();
    let motor_current = mid[8] / m.inductance_h;
    let relative = mid[9] / m.rotor_inertia_kg_m2 - mid[10] / m.stator_inertia_kg_m2;
    let mut heat = [0.0; 14];
    heat[0] = dt * p.winding_resistance_ohm * current * current;
    for i in 0..4 {
        let device = p.devices[i];
        let temperature = next.thermal().temperature_k(IDS[i + 1]).unwrap();
        near(device.body_diode.forward_voltage_v, 0.5 + 0.02 * (temperature - 300.0));
        let rd = device.body_diode.differential_resistance_ohm;
        let diode = ((-vds[i] - device.body_diode.forward_voltage_v) / rd).max(0.0);
        let r = if gates.on[i] { device.on_resistance_ohm } else { device.off_resistance_ohm };
        heat[1 + i] = dt * vds[i] * vds[i] / r;
        heat[10 + i] = dt * -vds[i] * diode;
        near(step.converter.last_devices[i].diode_current_a, diode);
    }
    heat[5] = dt * b.series_resistance_ohm * cell_current * cell_current;
    heat[6] = dt * pol * pol / b.polarization_resistance_ohm;
    heat[7] = dt * next.drive().lead_resistance_ohm() * lead_current * lead_current;
    heat[8] = dt * m.winding_resistance_ohm * motor_current * motor_current;
    heat[9] = dt * m.bearing_drag_nm_s_per_rad * relative * relative;
    for (a, b) in step.heat_joules().into_iter().zip(heat) { near(a, b); }
    let before_cap = capacitor_energy(&state);
    let after_cap = capacitor_energy(end);
    let cap_change = core::array::from_fn::<_, 4, _>(|i| after_cap[i] - before_cap[i]);
    for (actual, expected) in step.converter.device_capacitor_energy_change_j.into_iter().zip(cap_change) {
        near(actual, expected);
    }
    assert!(cap_change.iter().map(|q| q.abs()).sum::<f64>() > 1e-4,
        "the fixture must exchange reversible device-capacitor energy");
    assert!(next.drive().devices()[2].parameters().body_diode.forward_voltage_v > 0.500001);
    for (i, q) in per_body(heat, 0.0).into_iter().enumerate() {
        near(f.base.cp * (next.thermal().temperature_k(IDS[i]).unwrap() - 300.0), q);
    }
    near(next.report().deposited_joule_j, heat.iter().sum());
    near(next.report().total_energy_residual_j, 0.0);

    // Continue from the actual eleven accepted charges with a different gate.
    // The native replay freezes the accepted endpoint's dissipative coefficients
    // and starts from precisely that state, so a hidden algebraic node reset fails.
    let other_gates = DeviceGateState { on: [false, false, true, false] };
    let continued = network.step_candidate(
        next.thermal(), 2.0 * dt, next.drive(),
        &[CommutationInterval { gates: other_gates, duration_s: dt }], &[steady.command],
        1, None, power_controls(), cycle_controls(), drive_controls(7), &[0.0; 7], || false,
    ).unwrap();
    let accepted = continued.drive();
    let original = replay_original(&steady.native, accepted);
    let mut replay = PreparedCapacitiveBatteryMotorCycle::new_controlled(
        original, accepted.parameters().devices, COSS,
        &[CommutationInterval { gates: other_gates, duration_s: dt }],
        1, power_controls(), cycle_controls(), || false,
    ).unwrap();
    let mut replay_end = [0.0; 11];
    let expected = replay.step_controlled(next.drive().state(), steady.command, &mut replay_end, || false).unwrap();
    for (a, b) in accepted.state().iter().zip(replay_end) { near(*a, b); }
    for (a, b) in continued.report().electrical.heat_joules().into_iter().zip(expected.heat_joules()) { near(a, b); }
    assert!((accepted.state()[6] - next.drive().state()[6]).abs() > 1e-8);
}
fn replay_original(original: &BatteryMotorDrive, accepted: &MaterialCapacitiveBatteryMotorDrive) -> BatteryMotorDrive {
    let p = accepted.parameters();
    AveragedBuckBoost::new(BuckBoostParameters {
        inductance_h: p.inductance_h, bus_capacitance_f: p.bus_capacitance_f,
        winding_resistance_ohm: p.winding_resistance_ohm,
        input_on_resistance_ohm: [p.devices[0].on_resistance_ohm, p.devices[1].on_resistance_ohm],
        output_on_resistance_ohm: [p.devices[2].on_resistance_ohm, p.devices[3].on_resistance_ohm],
    }, BridgeModulation { input_high_fraction: 0.0, output_high_fraction: 0.0 }).unwrap()
        .with_battery_source(BatterySource::new(original.electrical().battery_curve().clone(),
            accepted.battery_parameters()).unwrap(), accepted.lead_resistance_ohm()).unwrap()
        .with_dc_motor(DcMotor::new(accepted.motor_parameters()).unwrap()).unwrap()
}

#[test]
fn only_original_constant_capacitance_with_matching_package_charge_and_validity_is_admitted() {
    let steady = steady_drive(false);
    let state = initial(&steady, [0.0; 2]);
    for variable in [false, true] {
        let f = capacitive_fixture(1.0, 0.0, 350.0, variable, 1.0);
        let thermal = f.base.thermal();
        let context = common::point(300.0);
        let cells = f.base.cells(&context);
        let r = drive_resistors(&f.base);
        let d = f.base.devices();
        let caps = f.bindings();
        let network = coupled(&cells, &r, &d, &caps);
        let bound = network.bind(&thermal, &steady.native, COSS, state, || false);
        if variable { assert!(bound.is_err()); continue; }
        let drive = bound.unwrap();
        let mut altered = COSS; altered[0] *= 1.01;
        assert!(network.bind(&thermal, &steady.native, altered, state, || false).is_err());
        let missing = PropertyKey::new("undeclared_output_capacitance", Capacitance::DIMS);
        for key in [&missing, &f.base.keys[0]] {
            let mut bad = f.bindings(); bad[0].key = key;
            assert!(coupled(&cells, &r, &d, &bad).bind(&thermal, &steady.native, COSS, state, || false).is_err());
        }
        let mut wrong_package = f.bindings(); wrong_package[0].cell = IDS[0];
        assert!(coupled(&cells, &r, &d, &wrong_package).bind(&thermal, &steady.native, COSS, state, || false).is_err());
        let changed = capacitive_fixture(1.0, 0.0, 350.0, false, 1.1);
        let changed_cells = changed.base.cells(&context);
        let changed_r = drive_resistors(&changed.base);
        let changed_d = changed.base.devices();
        let changed_caps = changed.bindings();
        let changed_network = coupled(&changed_cells, &changed_r, &changed_d, &changed_caps);
        let dt = 1.0 / 64.0;
        assert!(changed_network.step_candidate(
            &thermal, dt, &drive, &[CommutationInterval { gates: steady.gates, duration_s: dt }],
            &[steady.command], 4, None, power_controls(), cycle_controls(),
            drive_controls(7), &[0.0; 7], || false,
        ).is_err());
        assert_eq!(drive.state(), &state);
        assert_eq!(drive.time_s(), 0.0);
    }

    // Constant C is still qualified only within its ORIGINAL thermal domain.
    let f = capacitive_fixture(1.0, 0.0, 300.000001, false, 1.0);
    let thermal = f.base.thermal();
    let context = common::point(300.0);
    let cells = f.base.cells(&context);
    let r = drive_resistors(&f.base); let d = f.base.devices(); let caps = f.bindings();
    let network = coupled(&cells, &r, &d, &caps);
    let drive = network.bind(&thermal, &steady.native, COSS, state, || false).unwrap();
    let before = thermal.cells().clone();
    let dt = 1.0 / 64.0;
    assert!(network.step_candidate(
        &thermal, dt, &drive, &[CommutationInterval { gates: steady.gates, duration_s: dt }],
        &[steady.command], 4, None, power_controls(), cycle_controls(),
        drive_controls(7), &[0.0; 7], || false,
    ).is_err());
    assert_eq!(thermal.cells(), &before);
    assert_eq!(drive.state(), &state);
}

#[test]
fn capacitive_drive_heat_assembly_retains_both_finite_thermal_inventories() {
    let f = capacitive_fixture(1.0, 0.0, 350.0, false, 1.0);
    let base = f.base.thermal();
    let mut events = base.cells().values().map(|cell| TimedEvent {
        time_s: 0.0, event: Event::Birth(cell.clone()),
    }).collect::<Vec<_>>();
    events.push(TimedEvent { time_s: 0.0, event: Event::Contact(ThermalLink {
        id: 77, a: IDS[0], b: IDS[6], area_m2: 1.0, conductance_w_k: 0.0,
    }) });
    let thermal = DepositionSimulation::new(
        events, &[(MaterialId(0), f.base.caloric.curve() as &dyn ThermalLaw)],
        DepositionBudget { max_cells: 7, max_links: 1, max_events: 8, max_substeps: 1000 },
    ).unwrap().advance_to(0.0, StepConfig { max_dt_s: 0.1, energy_tolerance_j: 1e-9 }, || Ok(())).unwrap().0;
    let context = common::point(300.0);
    let cells = f.base.cells(&context);
    let bulk = MaterialNetwork { cells: &cells, links: &[], resistors: &[] };
    let electrical = ElectricalContactNetwork { bulk: &bulk, contacts: &[], electrical_bases: &[] };
    let heat = MaterialHeatNetwork { electrical: &electrical, thermal_contacts: &[], radiation: None };
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(&gate, arena, fs_exec::StreamKey {
            seed: 31, kernel_id: 83, tile: 0, iteration: 0,
        }, fs_exec::Budget::INFINITE, fs_exec::ExecMode::Deterministic);
        let prepared = heat.prepare_with_link_conductances(&cx, &thermal, &BTreeMap::from([(77, 0.0)])).unwrap();
        let transport = FixedCooling { prepared: &prepared, conductance_w_k: 0.25 };
        let resistors = drive_resistors(&f.base); let devices = f.base.devices(); let caps = f.bindings();
        let network = CapacitiveBatteryMotorHeatNetwork {
            drive: CommutatedBatteryMotorHeatNetwork {
                assembly: &transport, resistors: &resistors, devices: &devices, bearing_cell: IDS[0],
            }, capacitances: &caps,
        };
        let steady = steady_drive(false);
        let state = initial(&steady, [0.0; 2]);
        let drive = network.bind(&cx, &thermal, &steady.native, COSS, state).unwrap();
        let dt = 1.0 / 64.0;
        let next = network.step_candidate(
            &cx, &thermal, dt, &drive,
            &[CommutationInterval { gates: steady.gates, duration_s: dt }], &[steady.command],
            4, None, power_controls(), cycle_controls(), drive_controls(7), &[0.0; 7],
        ).unwrap();
        let q = per_body(steady.heat_w.map(|p| dt * p), 0.0);
        // Independent finite-inventory elimination: both C=1J/K, T0=300K.
        // X=h*G*(Q0-Q6)/(1+2*h*G); transport neither creates nor discards heat.
        let exchange = dt * 0.25 * (q[0] - q[6]) / (1.0 + 2.0 * dt * 0.25);
        near(next.thermal().temperature_k(IDS[0]).unwrap(), 300.0 + q[0] - exchange);
        near(next.thermal().temperature_k(IDS[6]).unwrap(), 300.0 + q[6] + exchange);
        for i in 1..6 { near(next.thermal().temperature_k(IDS[i]).unwrap(), 300.0 + q[i]); }
        near(next.thermal().energy_j() - thermal.energy_j(), q.iter().sum());
        near(next.report().deposited_joule_j, q.iter().sum());
        near(next.report().total_energy_residual_j, 0.0);
        assert_eq!(next.report().source.drive.network.material.cells.len(), 7);
        assert_eq!(next.report().source.capacitances.len(), 4);
        for value in next.report().electrical.converter.device_capacitor_energy_change_j { near(value, 0.0); }
    });
}

// Constant capacitors bound under SingleClaimOnly are refreshed under a pin of the same claim by the
// first step, which must be accepted; different capacitors are different sources.
#[test]
fn capacitors_bound_under_single_claim_only_take_their_first_step_and_differ_from_each_other() {
    let f = capacitive_fixture(1.0, 0.0, 350.0, false, 1.0);
    let thermal = f.base.thermal();
    let context = common::point(300.0);
    let cells = f.base.cells(&context);
    let resistors = drive_resistors(&f.base);
    let devices = f.base.devices();
    let caps = f.bindings();
    let network = coupled(&cells, &resistors, &devices, &caps);
    let steady = steady_drive(false);
    let state = initial(&steady, [0.0; 2]);
    let drive = network.bind(&thermal, &steady.native, COSS, state, || false).unwrap();
    let dt = 1.0 / 64.0;
    match network.step_candidate(
        &thermal, dt, &drive,
        &[CommutationInterval { gates: steady.gates, duration_s: dt }], &[steady.command],
        4, Some(IDS[3]), power_controls(), cycle_controls(), drive_controls(7), &[0.0; 7], || false,
    ) {
        Ok(_) => {}
        Err(error) => panic!("a step from the bound state is refused: {error:?}"),
    }
    // (Sources refreshed by the step sit at the stepped state, so their query point differs from the
    // bound one by design; the accepted step above is the positive.)
    // Four capacitors of four packages: another key, cell and value is another source.
    assert_ne!(drive.capacitances()[0], drive.capacitances()[1]);
}

// Planted negative on the same check (see the converter's): the drive stays bound at 300 K while the
// accepted heat state is another one at the same time, and the diode threshold depends on temperature.
#[test]
fn capacitive_drive_bound_at_one_state_refuses_a_step_from_another_accepted_state() {
    let f = capacitive_fixture(1.0, 0.02, 350.0, false, 1.0);
    let thermal = f.base.thermal();
    let hot = f.base.thermal_at(320.0);
    assert!((hot.temperature_k(IDS[1]).unwrap() - 320.0).abs() < 1e-5);
    let context = common::point(300.0);
    let cells = f.base.cells(&context);
    let resistors = drive_resistors(&f.base);
    let devices = f.base.devices();
    let caps = f.bindings();
    let network = coupled(&cells, &resistors, &devices, &caps);
    let steady = steady_drive(false);
    let state = initial(&steady, [0.0; 2]);
    let drive = network.bind(&thermal, &steady.native, COSS, state, || false).unwrap();
    let dt = 1.0 / 64.0;
    let attempt = |at: &DepositionSimulation<'_>| {
        network.step_candidate(
            at, dt, &drive,
            &[CommutationInterval { gates: steady.gates, duration_s: dt }], &[steady.command],
            4, Some(IDS[3]), power_controls(), cycle_controls(), drive_controls(7), &[0.0; 7], || false,
        ).map(|_| ()).map_err(|error| format!("{error:?}"))
    };
    assert_eq!(attempt(&thermal), Ok(()));
    let Err(text) = attempt(&hot) else { panic!("a step from another accepted state was taken"); };
    assert!(
        text.contains("commutated drive must retain every accepted resistance, conduction and capacitance source"),
        "{text}"
    );
}
