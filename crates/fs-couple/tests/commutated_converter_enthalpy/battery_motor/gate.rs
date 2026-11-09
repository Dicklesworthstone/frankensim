//! Finite gate-source work, original G(Vgs,T), Miller storage and distinct caloric receivers.
use super::*;

#[path = "gate/current.rs"]
mod current;
use fs_couple::circuit_enthalpy::converter::commutation::battery_motor::gate::*;
use fs_phs::circuit::capacitance::CapacitorChargeCurve;
use fs_phs::circuit::converter::commutation::{
    capacitance::{battery_motor::PreparedGateBatteryMotorCycle, gate::GateCommutationParameters},
    gate::{GateConductanceCurve, GateDriver},
};
use fs_qty::{Capacitance, Conductance, ElectricCharge};

const DRIVER: [CellId; 4] = [CellId(80), CellId(90), CellId(100), CellId(110)];
const CDS: [f64; 4] = [0.025, 0.04, 0.035, 0.05];
const G0: [f64; 4] = [10.0, 0.01, 0.01, 0.01];

fn gate_context() -> QueryPoint {
    common::point(300.0)
        .with_quantity("Vds", QuantitySpec::dimensional(Voltage::DIMS), 0.0)
        .unwrap()
        .with_quantity("Vgs", QuantitySpec::dimensional(Voltage::DIMS), 0.0)
        .unwrap()
        .with("lot", 7.0)
        .unwrap()
}
fn g(device: usize, gate_v: f64, temperature_k: f64, feedback: bool) -> f64 {
    G0[device]
        * if feedback {
            1.0 + 0.1 * gate_v + 0.05 * (temperature_k - 300.0)
        } else {
            1.0
        }
}
fn drivers() -> [GateDriver; 4] {
    core::array::from_fn(|i| GateDriver {
        on_voltage_v: if i == 0 { 2.0 } else { 0.0 },
        off_voltage_v: if i == 0 { 1.0 } else { 0.0 },
        on_resistance_ohm: 0.25,
        off_resistance_ohm: 2.0,
        maximum_current_a: 20.0,
    })
}
fn authored(miller_f: f64, feedback: bool) -> GateCommutationParameters {
    GateCommutationParameters {
        converter: parameters(),
        intrinsic_drain_charge_curves: core::array::from_fn(|i| {
            CapacitorChargeCurve::try_new(vec![
                (-100.0 * CDS[i], -100.0),
                (0.0, 0.0),
                (100.0 * CDS[i], 100.0),
            ])
            .unwrap()
        }),
        gate_source_capacitance_f: [1.0; 4],
        gate_drain_capacitance_f: [miller_f; 4],
        channel_conductance: core::array::from_fn(|i| {
            GateConductanceCurve::try_new(
                [-4.0, 0.0, 4.0]
                    .map(|v| (v, g(i, v, 300.0, feedback)))
                    .to_vec(),
            )
            .unwrap()
        }),
        gate_drivers: drivers(),
        maximum_cells_per_chord: 4096,
    }
}
struct GateFixture {
    base: Fixture,
    charge: [PropertyKey; 4],
    capacitance: [PropertyKey; 2],
    conductance: [PropertyKey; 4],
    driver_resistance: [PropertyKey; 2],
}
fn gate_fixture(
    miller_f: f64,
    feedback: bool,
    resistance_slope: f64,
    variable_storage: bool,
    bias_condition: bool,
) -> GateFixture {
    let base = fixture(1.0, 0.0, 350.0, 100.0);
    let charge = core::array::from_fn(|i| {
        PropertyKey::new(format!("intrinsic_drain_charge_{i}"), ElectricCharge::DIMS)
    });
    let capacitance = [
        "intrinsic_gate_source_capacitance",
        "intrinsic_gate_drain_capacitance",
    ]
    .map(|key| PropertyKey::new(key, Capacitance::DIMS));
    let conductance = core::array::from_fn(|i| {
        PropertyKey::new(format!("gate_channel_conductance_{i}"), Conductance::DIMS)
    });
    let driver_resistance = ["gate_driver_on_resistance", "gate_driver_off_resistance"]
        .map(|key| PropertyKey::new(key, Resistance::DIMS));
    let thermal_domain = || {
        ValidityDomain::unconstrained().with_quantity(
            "T",
            QuantitySpec::dimensional(Temperature::DIMS),
            299.0,
            350.0,
        )
    };
    let mut claims = ClaimSet::new();
    for (_, claim) in base.card.claims().claims_ordered() {
        claims.insert_claim(claim.clone()).unwrap();
    }
    let mut add = |key: &PropertyKey, value, validity| {
        claims
            .insert_claim(PropertyClaim {
                key: key.clone(),
                value,
                validity,
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source:
                        "synthetic intrinsic gate/storage fixture, not qualified transistor data"
                            .into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    };
    for i in 0..4 {
        add(
            &charge[i],
            PropertyValue::Curve {
                abscissa: "Vds".into(),
                abscissa_dims: Voltage::DIMS,
                knots: vec![
                    (-100.0, -100.0 * CDS[i]),
                    (0.0, 0.0),
                    (100.0, 100.0 * CDS[i]),
                ],
                dims: ElectricCharge::DIMS,
            },
            thermal_domain().with_quantity(
                "Vds",
                QuantitySpec::dimensional(Voltage::DIMS),
                -100.0,
                100.0,
            ),
        );
        let value = if feedback {
            let mut samples = Vec::new();
            for t in [299.0, 350.0] {
                samples.extend([-4.0, 0.0, 4.0].map(|v| ([v, t], g(i, v, t, true))));
            }
            PropertyValue::Surface2D {
                axes: ["Vgs".into(), "T".into()],
                axis_dims: [Voltage::DIMS, Temperature::DIMS],
                samples,
                triangles: vec![[0, 1, 3], [1, 4, 3], [1, 2, 4], [2, 5, 4]],
                dims: Conductance::DIMS,
            }
        } else {
            PropertyValue::Curve {
                abscissa: "Vgs".into(),
                abscissa_dims: Voltage::DIMS,
                knots: vec![(-4.0, G0[i]), (0.0, G0[i]), (4.0, G0[i])],
                dims: Conductance::DIMS,
            }
        };
        let mut domain = thermal_domain().with_quantity(
            "Vgs",
            QuantitySpec::dimensional(Voltage::DIMS),
            -4.0,
            4.0,
        );
        if bias_condition {
            domain =
                domain.with_quantity("Vds", QuantitySpec::dimensional(Voltage::DIMS), 0.0, 0.0);
        }
        add(&conductance[i], value, domain);
    }
    for (key, value) in capacitance.iter().zip([1.0, miller_f]) {
        let law = if variable_storage {
            PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(299.0, value), (350.0, value)],
                dims: Capacitance::DIMS,
            }
        } else {
            PropertyValue::Scalar {
                value,
                dims: Capacitance::DIMS,
            }
        };
        add(key, law, thermal_domain());
    }
    for (i, key) in driver_resistance.iter().enumerate() {
        let (value, slope) = if i == 0 {
            (0.25, resistance_slope)
        } else {
            (2.0, 0.0)
        };
        add(
            key,
            PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(299.0, value - slope), (350.0, value + 50.0 * slope)],
                dims: Resistance::DIMS,
            },
            thermal_domain(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic-device".into(),
            phase: "solid".into(),
            process: "authored intrinsic gate decomposition".into(),
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
    GateFixture {
        base: Fixture {
            card,
            caloric,
            keys: base.keys,
            cp: 1.0,
        },
        charge,
        capacitance,
        conductance,
        driver_resistance,
    }
}
impl GateFixture {
    fn thermal(&self) -> DepositionSimulation<'_> {
        let mut events: Vec<_> = self
            .base
            .thermal()
            .cells()
            .values()
            .map(|cell| TimedEvent {
                time_s: 0.0,
                event: Event::Birth(cell.clone()),
            })
            .collect();
        for (i, id) in DRIVER.into_iter().enumerate() {
            events.push(TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id,
                    segment_id: 7 + i as u64,
                    body_id: 7 + i as u64,
                    material_id: MaterialId(0),
                    mass_kg: 1.0,
                    specific_h_j_kg: 1.0,
                    geometry: CellGeometry {
                        start_m: [20.0 + i as f64, 0.0, 0.0],
                        end_m: [21.0 + i as f64, 0.0, 0.0],
                        volume_m3: 1.0,
                        surface_area_m2: 6.0,
                    },
                    boundary: common::boundary(),
                }),
            });
        }
        DepositionSimulation::new(
            events,
            &[(MaterialId(0), self.base.caloric.curve() as &dyn ThermalLaw)],
            DepositionBudget {
                max_cells: 11,
                max_links: 0,
                max_events: 11,
                max_substeps: 1000,
            },
        )
        .unwrap()
        .advance_to(
            0.0,
            StepConfig {
                max_dt_s: 1.0,
                energy_tolerance_j: 1e-9,
            },
            || Ok(()),
        )
        .unwrap()
        .0
    }
    fn cells<'a>(&'a self, context: &'a QueryPoint) -> Vec<MaterialCell<'a>> {
        let mut cells = self.base.cells(context).into_iter().collect::<Vec<_>>();
        cells.extend(DRIVER.map(|cell| MaterialCell {
            cell,
            material_id: MaterialId(0),
            curve: self.base.caloric.curve(),
            point: ThermalMaterialPoint {
                material: &self.base.card,
                context,
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            conductivity: None,
        }));
        cells
    }
    fn charges(&self) -> [DeviceChargeBinding<'_>; 4] {
        core::array::from_fn(|i| DeviceChargeBinding {
            cell: IDS[i + 1],
            key: &self.charge[i],
            voltage_axis: "Vds",
            voltage_range_v: [-100.0, 100.0],
            selection: MaterialPropertySelection::SingleClaimOnly,
        })
    }
    fn caps(&self, miller: bool) -> [DeviceCapacitanceBinding<'_>; 4] {
        core::array::from_fn(|i| DeviceCapacitanceBinding {
            cell: IDS[i + 1],
            key: &self.capacitance[usize::from(miller)],
            selection: MaterialPropertySelection::SingleClaimOnly,
        })
    }
    fn conductances(&self) -> [GateConductanceBinding<'_>; 4] {
        core::array::from_fn(|i| GateConductanceBinding {
            cell: IDS[i + 1],
            key: &self.conductance[i],
            gate_voltage_axis: "Vgs",
            gate_voltage_range_v: [-4.0, 4.0],
            selection: MaterialPropertySelection::SingleClaimOnly,
        })
    }
    fn driver_resistors(&self) -> [ConverterThermalChannel<'_>; 8] {
        core::array::from_fn(|i| ConverterThermalChannel {
            cell: DRIVER[i / 2],
            resistance: ConverterResistance::Resistance {
                key: &self.driver_resistance[i % 2],
                selection: MaterialPropertySelection::SingleClaimOnly,
            },
        })
    }
}
fn coupled<'a>(
    cells: &'a [MaterialCell<'a>],
    resistors: &'a [ConverterThermalChannel<'a>; 5],
    devices: &'a [DeviceConductionBinding<'a>; 4],
    gate: GateMaterialBindings<'a>,
) -> GateBatteryMotorThermalNetwork<'a> {
    GateBatteryMotorThermalNetwork {
        drive: CommutatedBatteryMotorThermalNetwork {
            cells,
            links: &[],
            resistors,
            devices,
            bearing_cell: IDS[0],
        },
        gate,
        maximum_cells_per_chord: 4096,
    }
}
fn initial(steady: &SteadyDrive, p: &GateCommutationParameters) -> [f64; 15] {
    let current = steady.state[3] / p.converter.inductance_h;
    let input = steady.state[4] / p.converter.bus_capacitance_f[0];
    let output = steady.state[5] / p.converter.bus_capacitance_f[1];
    let left = (10.0 * input - current) / 10.01;
    let right = left - 0.1 * current;
    let converter = p
        .state(current, [input, output], [left, right], [0.0; 4])
        .unwrap();
    let mut state = [0.0; 15];
    state[..3].copy_from_slice(&steady.state[..3]);
    state[3..12].copy_from_slice(&converter);
    state[12..].copy_from_slice(&steady.state[6..]);
    state
}
fn gate_bindings<'a>(
    charge: &'a [DeviceChargeBinding<'a>; 4],
    cgs: &'a [DeviceCapacitanceBinding<'a>; 4],
    cgd: &'a [DeviceCapacitanceBinding<'a>; 4],
    conductance: &'a [GateConductanceBinding<'a>; 4],
    resistors: &'a [ConverterThermalChannel<'a>; 8],
) -> GateMaterialBindings<'a> {
    GateMaterialBindings {
        intrinsic_drain_charge: charge,
        gate_source_capacitance: cgs,
        gate_drain_capacitance: cgd,
        channel_conductance: conductance,
        driver_resistors: resistors,
        drivers: drivers(),
    }
}

#[test]
fn finite_driver_endpoint_feedback_and_recovered_gate_work_have_distinct_receivers() {
    let f = gate_fixture(0.0, false, 0.125, false, false);
    let thermal = f.thermal();
    let context = gate_context();
    let cells = f.cells(&context);
    let resistors = drive_resistors(&f.base);
    let devices = f.base.devices();
    let q = f.charges();
    let cgs = f.caps(false);
    let cgd = f.caps(true);
    let conductance = f.conductances();
    let driver_r = f.driver_resistors();
    let network = coupled(
        &cells,
        &resistors,
        &devices,
        gate_bindings(&q, &cgs, &cgd, &conductance, &driver_r),
    );
    let steady = steady_drive(false);
    let p = authored(0.0, false);
    let state = initial(&steady, &p);
    let drive = network
        .bind(&thermal, &steady.native, &p, state, || false)
        .unwrap();
    let next = network
        .step_candidate(
            &thermal,
            1.0,
            &drive,
            &[CommutationInterval {
                gates: steady.gates,
                duration_s: 1.0,
            }],
            &[steady.command],
            1,
            None,
            power_controls(),
            cycle_controls(),
            drive_controls(11),
            &[0.0; 11],
            || false,
        )
        .unwrap();
    let step = &next.report().electrical;
    let gate = step.converter.gate.as_ref().unwrap();
    // Independent implicit-midpoint/caloric solution: C=1, U=2, h=1,
    // R(T)=.25+.125*(T-300), caloric capacity1. At T=302: R=.5,
    // v1=2, Ibar=2, Wsource=4, Eg=2 and physical driver heat=2.
    near(next.drive().state()[8], 2.0);
    near(
        next.drive().gate_material().drivers()[0].on_resistance_ohm,
        0.5,
    );
    near(next.thermal().temperature_k(DRIVER[0]).unwrap(), 302.0);
    near(gate.source_work_j[0], 4.0);
    near(gate.resistor_heat_j[0], 2.0);
    near(gate.gate_source_energy_change_j[0], 2.0);
    for i in 1..4 {
        near(gate.source_work_j[i], 0.0);
        near(gate.resistor_heat_j[i], 0.0);
        near(next.thermal().temperature_k(DRIVER[i]).unwrap(), 300.0);
    }
    for (actual, expected) in step.heat_joules().into_iter().zip(steady.heat_w) {
        near(actual, expected);
    }
    for (i, expected) in per_body(steady.heat_w, 0.0).into_iter().enumerate() {
        near(
            next.thermal().temperature_k(IDS[i]).unwrap() - 300.0,
            expected,
        );
    }
    near(
        next.report().deposited_joule_j,
        steady.heat_w.iter().sum::<f64>() + 2.0,
    );
    near(step.storage_change_j[0], -12.0 * steady.cell_current_a);
    near(next.report().total_energy_residual_j, 0.0);
    // OFF is a finite U=1,R=2 source branch, not a gate reset or forced Ron/Roff.
    // Starting at2V, midpoint yields1.6V, Wsource=-.4J and resistor heat.32J.
    let off = network
        .step_candidate(
            next.thermal(),
            2.0,
            next.drive(),
            &[CommutationInterval {
                gates: DeviceGateState { on: [false; 4] },
                duration_s: 1.0,
            }],
            &[steady.command],
            1,
            None,
            power_controls(),
            cycle_controls(),
            drive_controls(11),
            &[0.0; 11],
            || false,
        )
        .unwrap();
    let og = off.report().electrical.converter.gate.as_ref().unwrap();
    near(off.drive().state()[8], 1.6);
    near(og.source_work_j[0], -0.4);
    near(og.resistor_heat_j[0], 0.32);
    near(og.gate_source_energy_change_j[0], -0.72);
    near(off.thermal().temperature_k(DRIVER[0]).unwrap(), 302.32);
    for (actual, expected) in off
        .report()
        .electrical
        .heat_joules()
        .into_iter()
        .zip(steady.heat_w)
    {
        near(actual, expected);
    }
    assert_eq!(drive.state(), &state);
    assert_eq!(thermal.time_s(), 0.0);
    for (a, b) in next
        .drive()
        .gate_material()
        .intrinsic_drain_charge()
        .iter()
        .zip(off.drive().gate_material().intrinsic_drain_charge())
    {
        assert_eq!(a.curve(), b.curve());
        assert_eq!(a.source().card_identity(), b.source().card_identity());
    }
    near(off.report().total_energy_residual_j, 0.0);
}

// Independent linear intrinsic branch electrostatics, not native energy queries.
// qg=(Cgs+Cgd)*Vgs-Cgd*Vds. Eliminating the two gates preserves the
// reciprocal bus/node charge matrix with Cds+Cgs*Cgd/(Cgs+Cgd).
// The native grouped node contains the high-side gate electrode as well.
fn voltages(state: &[f64; 15], miller: f64) -> ([f64; 4], [f64; 4]) {
    let mut drain = [0.0; 4];
    let mut gate = [0.0; 4];
    let sum = 1.0 + miller;
    let alpha = miller / sum;
    for leg in 0..2 {
        let hi = 2 * leg;
        let lo = hi + 1;
        let ch = CDS[hi] + miller / sum;
        let cl = CDS[lo] + miller / sum;
        let cb = parameters().bus_capacitance_f[leg];
        let qb = state[4 + leg] + alpha * state[8 + hi];
        let qn = state[6 + leg] - alpha * state[8 + hi] + alpha * state[8 + lo];
        let det = cb * (ch + cl) + ch * cl;
        let bus = ((ch + cl) * qb + ch * qn) / det;
        let node = (ch * qb + (cb + ch) * qn) / det;
        drain[hi] = bus - node;
        drain[lo] = node;
        gate[hi] = (state[8 + hi] + miller * drain[hi]) / sum;
        gate[lo] = (state[8 + lo] + miller * drain[lo]) / sum;
    }
    (drain, gate)
}
fn replay_model(drive: &MaterialGateBatteryMotorDrive) -> BatteryMotorDrive {
    let p = drive.parameters();
    let battery =
        BatterySource::new(drive.battery_curve().clone(), drive.battery_parameters()).unwrap();
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
    .with_battery_source(battery, drive.lead_resistance_ohm())
    .unwrap()
    .with_dc_motor(DcMotor::new(drive.motor_parameters()).unwrap())
    .unwrap()
}

#[test]
fn thermal_gate_surface_and_reciprocal_miller_charge_reach_the_same_native_motion() {
    let f = gate_fixture(0.05, true, 0.125, false, false);
    let thermal = f.thermal();
    let context = gate_context();
    let cells = f.cells(&context);
    let resistors = drive_resistors(&f.base);
    let devices = f.base.devices();
    let q = f.charges();
    let cgs = f.caps(false);
    let cgd = f.caps(true);
    let conductance = f.conductances();
    let driver_r = f.driver_resistors();
    let network = coupled(
        &cells,
        &resistors,
        &devices,
        gate_bindings(&q, &cgs, &cgd, &conductance, &driver_r),
    );
    let steady = steady_drive(false);
    let p = authored(0.05, true);
    let state = initial(&steady, &p);
    let drive = network
        .bind(&thermal, &steady.native, &p, state, || false)
        .unwrap();
    let h = 0.125;
    let intervals = [CommutationInterval {
        gates: steady.gates,
        duration_s: h,
    }];
    let next = network
        .step_candidate(
            &thermal,
            h,
            &drive,
            &intervals,
            &[steady.command],
            1,
            None,
            power_controls(),
            cycle_controls(),
            drive_controls(11),
            &[0.0; 11],
            || false,
        )
        .unwrap();
    let step = &next.report().electrical;
    let gate = step.converter.gate.as_ref().unwrap();
    let (v0, g0) = voltages(&state, 0.05);
    let (v1, g1) = voltages(next.drive().state(), 0.05);
    let mut miller_storage = 0.0;
    for i in 0..4 {
        let temp = next.thermal().temperature_k(IDS[i + 1]).unwrap();
        let source = &next.drive().gate_material().channel_conductance()[i];
        near(
            source
                .curve()
                .evaluate(0.5 * (g0[i] + g1[i]))
                .unwrap()
                .conductance_s,
            g(i, 0.5 * (g0[i] + g1[i]), temp, true),
        );
        assert_eq!(
            source.source().card_identity(),
            drive.gate_material().channel_conductance()[i]
                .source()
                .card_identity()
        );
        near(
            source
                .source()
                .query_point()
                .iter()
                .find(|(axis, _)| axis == "T")
                .unwrap()
                .1,
            temp,
        );
        let vd = 0.5 * (v0[i] + v1[i]);
        let vg = 0.5 * (g0[i] + g1[i]);
        near(
            step.heat_joules()[1 + i],
            h * g(i, vg, temp, true) * vd * vd,
        );
        let driver = next.drive().gate_material().drivers()[i];
        let (u, r) = if steady.gates.on[i] {
            (driver.on_voltage_v, driver.on_resistance_ohm)
        } else {
            (driver.off_voltage_v, driver.off_resistance_ohm)
        };
        let ig = (u - vg) / r;
        near(gate.source_work_j[i], h * u * ig);
        near(gate.resistor_heat_j[i], h * r * ig * ig);
        near(
            next.thermal().temperature_k(DRIVER[i]).unwrap() - 300.0,
            gate.resistor_heat_j[i],
        );
        near(
            step.converter.device_capacitor_energy_change_j[i],
            0.5 * CDS[i] * (v1[i] - v0[i]) * (v1[i] + v0[i]),
        );
        near(
            gate.gate_source_energy_change_j[i],
            0.5 * (g1[i] - g0[i]) * (g1[i] + g0[i]),
        );
        let before = g0[i] - v0[i];
        let after = g1[i] - v1[i];
        let de = 0.025 * (after - before) * (after + before);
        near(gate.gate_drain_energy_change_j[i], de);
        miller_storage += de.abs();
    }
    assert!(miller_storage > 1e-6);
    assert!((v1[1] - v0[1]).abs() > 1e-6);
    let actual = next.drive().gate_parameters(4096);
    let mut replay = PreparedGateBatteryMotorCycle::new_gate_controlled(
        replay_model(next.drive()),
        actual,
        &intervals,
        1,
        power_controls(),
        cycle_controls(),
        || false,
    )
    .unwrap();
    let mut end = [0.0; 15];
    let raw = replay
        .step_interval_commands_controlled(&state, &[steady.command], &mut end, || false)
        .unwrap();
    for (a, b) in next.drive().state().iter().zip(end) {
        near(*a, b);
    }
    for (a, b) in raw
        .all_heat_joules()
        .into_iter()
        .zip(step.all_heat_joules())
    {
        near(a, b);
    }
    near(
        next.report().deposited_joule_j,
        step.all_heat_joules().iter().sum(),
    );
    near(next.report().total_energy_residual_j, 0.0);
}

#[test]
fn storage_authority_gate_bias_and_driver_receivers_are_retained_without_state_reinterpretation() {
    for (variable_storage, bias_condition) in [(true, false), (false, true)] {
        let f = gate_fixture(0.0, false, 0.0, variable_storage, bias_condition);
        let thermal = f.thermal();
        let context = gate_context();
        let cells = f.cells(&context);
        let resistors = drive_resistors(&f.base);
        let devices = f.base.devices();
        let q = f.charges();
        let cgs = f.caps(false);
        let cgd = f.caps(true);
        let conductance = f.conductances();
        let driver_r = f.driver_resistors();
        let network = coupled(
            &cells,
            &resistors,
            &devices,
            gate_bindings(&q, &cgs, &cgd, &conductance, &driver_r),
        );
        let steady = steady_drive(false);
        let p = authored(0.0, false);
        let state = initial(&steady, &p);
        // A flat sampled C(T) is still not an original constant storage law.
        // A G(Vgs) measurement restricted to Vds=0 cannot own a varying drain bias.
        assert!(
            network
                .bind(&thermal, &steady.native, &p, state, || false)
                .is_err()
        );
        assert_eq!(thermal.time_s(), 0.0);
    }
    let f = gate_fixture(0.0, false, 0.0, false, false);
    let thermal = f.thermal();
    let context = gate_context();
    let cells = f.cells(&context);
    let resistors = drive_resistors(&f.base);
    let devices = f.base.devices();
    let q = f.charges();
    let cgs = f.caps(false);
    let cgd = f.caps(true);
    let conductance = f.conductances();
    let driver_r = f.driver_resistors();
    let bindings = gate_bindings(&q, &cgs, &cgd, &conductance, &driver_r);
    let network = coupled(&cells, &resistors, &devices, bindings);
    let steady = steady_drive(false);
    let p = authored(0.0, false);
    let state = initial(&steady, &p);
    let drive = network
        .bind(&thermal, &steady.native, &p, state, || false)
        .unwrap();
    for gate_source in [false, true] {
        let mut changed = p.clone();
        if gate_source {
            changed.gate_source_capacitance_f[0] = 2.0;
        } else {
            changed.gate_drain_capacitance_f[0] = 0.01;
        }
        assert!(
            network
                .bind(&thermal, &steady.native, &changed, state, || false)
                .is_err()
        );
    }
    let mut changed = p.clone();
    changed.intrinsic_drain_charge_curves[0] = CapacitorChargeCurve::try_new(vec![
        (-100.0 * CDS[0], -100.0),
        (0.0, 0.0),
        (100.1 * CDS[0], 100.0),
    ])
    .unwrap();
    assert!(
        network
            .bind(&thermal, &steady.native, &changed, state, || false)
            .is_err()
    );

    let mut wrong_g = f.conductances();
    wrong_g[0].key = &f.base.keys[0];
    let wrong_quantity = coupled(
        &cells,
        &resistors,
        &devices,
        GateMaterialBindings {
            channel_conductance: &wrong_g,
            ..bindings
        },
    );
    assert!(
        wrong_quantity
            .bind(&thermal, &steady.native, &p, state, || false)
            .is_err()
    );
    let mut wrong_package = f.conductances();
    wrong_package[0].cell = IDS[2];
    let permuted = coupled(
        &cells,
        &resistors,
        &devices,
        GateMaterialBindings {
            channel_conductance: &wrong_package,
            ..bindings
        },
    );
    assert!(
        permuted
            .bind(&thermal, &steady.native, &p, state, || false)
            .is_err()
    );
    let mut split = f.driver_resistors();
    split[1].cell = DRIVER[1];
    let ambiguous_heat = coupled(
        &cells,
        &resistors,
        &devices,
        GateMaterialBindings {
            driver_resistors: &split,
            ..bindings
        },
    );
    assert!(
        ambiguous_heat
            .bind(&thermal, &steady.native, &p, state, || false)
            .is_err()
    );

    let changed_context = context.clone().with("lot", 8.0).unwrap();
    let changed_cells = f.cells(&changed_context);
    let changed_network = coupled(&changed_cells, &resistors, &devices, bindings);
    assert!(
        changed_network
            .step_candidate(
                &thermal,
                0.125,
                &drive,
                &[CommutationInterval {
                    gates: steady.gates,
                    duration_s: 0.125
                }],
                &[steady.command],
                1,
                None,
                power_controls(),
                cycle_controls(),
                drive_controls(11),
                &[0.0; 11],
                || false
            )
            .is_err()
    );
    assert_eq!(drive.state(), &state);
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn source_bound_driver_heat_and_finite_liquid_coolant_use_the_same_heat_assembly() {
    const COOLANT: CellId = CellId(120);
    let f = gate_fixture(0.0, false, 0.175, false, false);
    let mut fluid_claims = ClaimSet::new();
    for (_, claim) in f.base.card.claims().claims_ordered() {
        fluid_claims.insert_claim(claim.clone()).unwrap();
    }
    let fluid = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic-coolant".into(),
            phase: "liquid".into(),
            process: "authored finite driver-cooling inventory".into(),
            revision: 0,
        },
        fluid_claims,
        vec![],
    )
    .unwrap();
    let fluid_curve = ResolvedHeatCapacityCurve::try_new(
        &fluid,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &common::point(299.0),
            upper: &common::point(350.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Liquid,
        0.0,
        1e-6,
    )
    .unwrap();
    let base = f.thermal();
    let mut events: Vec<_> = base
        .cells()
        .values()
        .map(|cell| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell.clone()),
        })
        .collect();
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Birth(DepositedCell {
            id: COOLANT,
            segment_id: 11,
            body_id: 11,
            material_id: MaterialId(1),
            mass_kg: 2.0,
            specific_h_j_kg: 1.0,
            geometry: CellGeometry {
                start_m: [30.0, 0.0, 0.0],
                end_m: [32.0, 0.0, 0.0],
                volume_m3: 2.0,
                surface_area_m2: 10.0,
            },
            boundary: common::boundary(),
        }),
    });
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 77,
            a: DRIVER[0],
            b: COOLANT,
            area_m2: 1.0,
            conductance_w_k: 0.0,
        }),
    });
    let thermal = DepositionSimulation::new(
        events,
        &[
            (MaterialId(0), f.base.caloric.curve() as &dyn ThermalLaw),
            (MaterialId(1), fluid_curve.curve() as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 12,
            max_links: 1,
            max_events: 13,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-9,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    let context = gate_context();
    let fluid_context = common::point(300.0);
    let mut cells = f.cells(&context);
    cells.push(MaterialCell {
        cell: COOLANT,
        material_id: MaterialId(1),
        curve: fluid_curve.curve(),
        point: ThermalMaterialPoint {
            material: &fluid,
            context: &fluid_context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: None,
    });
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
    let cancel = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &cancel,
            arena,
            fs_exec::StreamKey {
                seed: 47,
                kernel_id: 97,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        let prepared = heat
            .prepare_with_link_conductances(&cx, &thermal, &BTreeMap::from([(77, 0.0)]))
            .unwrap();
        let cooling = FixedCooling {
            prepared: &prepared,
            conductance_w_k: 0.5,
        };
        let r = drive_resistors(&f.base);
        let d = f.base.devices();
        let q = f.charges();
        let cgs = f.caps(false);
        let cgd = f.caps(true);
        let conductance = f.conductances();
        let driver_r = f.driver_resistors();
        let network = GateBatteryMotorHeatNetwork {
            drive: CommutatedBatteryMotorHeatNetwork {
                assembly: &cooling,
                resistors: &r,
                devices: &d,
                bearing_cell: IDS[0],
            },
            gate: gate_bindings(&q, &cgs, &cgd, &conductance, &driver_r),
            maximum_cells_per_chord: 4096,
        };
        let steady = steady_drive(false);
        let p = authored(0.0, false);
        let state = initial(&steady, &p);
        let drive = network
            .bind(&cx, &thermal, &steady.native, &p, state)
            .unwrap();
        let next = network
            .step_candidate(
                &cx,
                &thermal,
                1.0,
                &drive,
                &[CommutationInterval {
                    gates: steady.gates,
                    duration_s: 1.0,
                }],
                &[steady.command],
                1,
                None,
                power_controls(),
                cycle_controls(),
                drive_controls(12),
                &[0.0; 12],
            )
            .unwrap();
        // Independent backward-Euler elimination: capacities1 and2, g=.5,h=1.
        // Qdriver=2 gives dTdriver=10/7, dTcoolant=2/7; R=.25+.175*(10/7)=.5.
        // This closes the same midpoint gate solution U2,C1 -> v1=2,W4,Eg2,Q2.
        near(
            next.thermal().temperature_k(DRIVER[0]).unwrap(),
            300.0 + 10.0 / 7.0,
        );
        near(
            next.thermal().temperature_k(COOLANT).unwrap(),
            300.0 + 2.0 / 7.0,
        );
        near(
            next.drive().gate_material().drivers()[0].on_resistance_ohm,
            0.5,
        );
        let gate = next.report().electrical.converter.gate.as_ref().unwrap();
        near(gate.resistor_heat_j[0], 2.0);
        near(gate.source_work_j[0], 4.0);
        near(next.drive().state()[8], 2.0);
        let whole_heat = steady.heat_w.iter().sum::<f64>() + 2.0;
        near(next.thermal().energy_j() - thermal.energy_j(), whole_heat);
        near(next.report().deposited_joule_j, whole_heat);
        near(next.report().non_electrical_heat_j, 0.0);
        near(next.report().total_energy_residual_j, 0.0);
        assert_eq!(next.report().source.drive.network.material.cells.len(), 12);
        for (i, expected) in per_body(steady.heat_w, 0.0).into_iter().enumerate() {
            near(
                next.thermal().temperature_k(IDS[i]).unwrap() - 300.0,
                expected,
            );
        }
    });
}
