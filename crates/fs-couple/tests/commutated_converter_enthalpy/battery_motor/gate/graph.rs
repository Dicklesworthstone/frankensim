//! The gate fixtures as declared graphs (`battery_motor/graph.rs`), the
//! test module of the gate lane that `battery_motor::gate` coupled by hand
//! until it retired.
//!
//! The drive runs in the `GateStorage` lane. Its storage is the package
//! card's, read once ([`card_storage`]), since the module held it constant.
//! Each switch's channel follows the card's G(Vgs, T) line at its package
//! ([`CardChannel`]); each gate's on and off resistances follow the card at
//! its driver cell. The four driver cells join the timeline after the seven,
//! each receiving its gate's driver heat, and eight `Force` sources hold the
//! rails at the module's driver voltages.

use std::sync::Arc;

use fs_matdb::{ClaimSelection, MaterialAnswer};
use fs_material::state_point::{
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_phs::circuit::capacitance::CapacitorChargeCurve;
use fs_phs::circuit::converter::commutation::{
    capacitance::{battery_motor::PreparedGateBatteryMotorCycle, gate::GateCommutationParameters},
    gate::{GateChannelLaw, GateChannelResponse, GateConductanceCurve, GateDriver},
};
use fs_qty::{Capacitance, Conductance, ElectricCharge};
use fs_system::commutated_drive::{
    ChannelLaw, CommutatedBatteryMotorDrive, GateDriverResistances, GateRail, GateStorage,
    SwitchPosition,
};

use super::graph::{
    Card, CardResistance, Case, Coolant, Drivers, Lane, Row, held_to_owner, kernel_rows, leak,
    owner_drive, text,
};
use super::nine_state::coolant;
use super::*;

#[path = "current/graph.rs"]
mod current;

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
fn gate_fixture(miller_f: f64, feedback: bool, resistance_slope: f64) -> GateFixture {
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
        // A scalar is a constant law: LinearInside has no knot span to interpolate inside on it.
        let interpolation = if matches!(value, PropertyValue::Scalar { .. }) {
            InterpolationPolicy::ConstantWithinValidity
        } else {
            InterpolationPolicy::LinearInside
        };
        claims
            .insert_claim(PropertyClaim {
                key: key.clone(),
                value,
                validity,
                uncertainty: UncertaintyModel::Unstated,
                interpolation,
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
        let domain = thermal_domain().with_quantity(
            "Vgs",
            QuantitySpec::dimensional(Voltage::DIMS),
            -4.0,
            4.0,
        );
        add(&conductance[i], value, domain);
    }
    for (key, value) in capacitance.iter().zip([1.0, miller_f]) {
        let law = PropertyValue::Scalar {
            value,
            dims: Capacitance::DIMS,
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
    /// All eleven bodies born at one common temperature (the caloric curve starts at 299 K), time 0.
    fn thermal_at(&self, kelvin: f64) -> DepositionSimulation<'_> {
        let mut events: Vec<_> = self
            .base
            .thermal_at(kelvin)
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
                    specific_h_j_kg: kelvin - 299.0,
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

/// The timeline's cells with ports: the seven, then the drivers'.
const CELLS: usize = 11;

impl Lane for GateStorage {
    /// Each gate's on and off resistances follow the card at its driver
    /// cell, on the case's two driver claims.
    fn add_drivers(
        mut drive: CommutatedBatteryMotorDrive<Self>,
        case: &Case<'_, Self>,
    ) -> CommutatedBatteryMotorDrive<Self> {
        let Some(drivers) = case.drivers else {
            return drive;
        };
        let card = Card {
            material: &case.fixture.card,
            curve: &case.fixture.caloric,
            context: case.context,
        };
        let switches = [
            SwitchPosition::InputHigh,
            SwitchPosition::InputLow,
            SwitchPosition::OutputHigh,
            SwitchPosition::OutputLow,
        ];
        let rails = switches
            .into_iter()
            .flat_map(|switch| [(switch, GateRail::On), (switch, GateRail::Off)]);
        let keys = drivers.keys.iter().cycle();
        for ((tag, (switch, rail)), key) in (5..).zip(rails).zip(keys) {
            let law = CardResistance { card, key, tag };
            drive = drive.driver_following(switch, rail, Arc::new(law)).unwrap();
        }
        drive
    }
}

/// A switch's channel from its package's card, as the module resolved it:
/// the card's G line over the module's gate voltages, -4 V to 4 V, at the
/// package's state, under the claim its lower end selects.
#[derive(Debug)]
struct CardChannel {
    card: Card,
    key: &'static PropertyKey,
    tag: u64,
}

impl CardChannel {
    /// The line's knots at `kelvin`: gate voltage \[V\], conductance \[S\].
    fn knots(&self, kelvin: f64) -> Result<Vec<(f64, f64)>, String> {
        let query = self.card.query(kelvin)?;
        let at = |volts| query.clone().with("Vgs", volts).map_err(text);
        let lower = at(-4.0)?;
        let rule = ScalarAdmissibility::StrictlyPositive;
        let requirement = ScalarPropertyRequirement::try_with_key(self.key, rule).map_err(text)?;
        let snapshot = resolve_material_state_point(
            self.card.material,
            &lower,
            &[requirement],
            MaterialPropertySelection::SingleClaimOnly,
        )
        .map_err(text)?;
        let property = snapshot.property_by_key(self.key).ok_or("no channel")?;
        let claim = ClaimSelection::Pinned(property.answer().receipt.selected);
        let claims = self.card.material.claims();
        let line = claims
            .query_line_typed(self.key, "Vgs", &lower, &at(4.0)?, claim)
            .map_err(text)?;
        let knot = |(volts, answer): &(f64, MaterialAnswer)| (*volts, answer.evidence.value.value);
        Ok(line.knots.iter().map(knot).collect())
    }
}

impl ChannelLaw for CardChannel {
    fn channel(&self, kelvin: f64) -> Result<GateChannelLaw, String> {
        let curve = GateConductanceCurve::try_new(self.knots(kelvin)?).map_err(text)?;
        GateChannelLaw::new(curve).map_err(text)
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        words.push(self.tag);
    }
}

/// The lane's storage from the package card at 300 K: each package's Q(Vds)
/// line over the module's -100 V to 100 V, its Cgs and its Cgd. The module
/// required the storage constant, so one reading stands for every step. The
/// channels follow [`CardChannel`]; the drivers' resistances are the
/// module's authored ones until their laws answer.
fn card_storage(f: &'static GateFixture, card: Card) -> GateStorage {
    use ScalarAdmissibility::{NonNegative, StrictlyPositive};
    let query = card.query(300.0).unwrap();
    let at = |volts| query.clone().with("Vds", volts).unwrap();
    let charge = |key| {
        let claim = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
        let claims = f.base.card.claims();
        let line = claims
            .query_line_typed(key, "Vds", &at(-100.0), &at(100.0), claim)
            .unwrap();
        let knots = line.knots.iter();
        let knots = knots.map(|(volts, answer)| (answer.evidence.value.value, *volts));
        CapacitorChargeCurve::try_new(knots.collect()).unwrap()
    };
    let [source, drain] = f.capacitance.each_ref();
    let wanted = [(source, StrictlyPositive), (drain, NonNegative)];
    let c = card.values(300.0, &wanted).unwrap();
    let mut tags = 0..;
    GateStorage {
        intrinsic_drain_charge: f.charge.each_ref().map(charge),
        gate_source_capacitance_f: [c[0]; 4],
        gate_drain_capacitance_f: [c[1]; 4],
        maximum_cells_per_chord: 4096,
        channels: f.conductance.each_ref().map(|key| {
            let tag = tags.next().unwrap();
            Arc::new(CardChannel { card, key, tag }) as Arc<dyn ChannelLaw>
        }),
        drivers: drivers().map(|driver| GateDriverResistances {
            on_resistance_ohm: driver.on_resistance_ohm,
            off_resistance_ohm: driver.off_resistance_ohm,
            maximum_current_a: driver.maximum_current_a,
        }),
    }
}

/// A gate case of the module's fixture `f` authored as `p`: the drive
/// started as the module started it, every cell at `start_k`, one step of
/// `dt_s` per program. The card's storage must be the authored storage, as
/// the module's bind required.
fn gate_case<'s>(
    f: &'static GateFixture,
    steady: &'s SteadyDrive,
    p: &GateCommutationParameters,
    start_k: f64,
    dt_s: f64,
) -> Case<'s, GateStorage> {
    let context = leak(gate_context());
    let card = Card {
        material: &f.base.card,
        curve: &f.base.caloric,
        context,
    };
    let storage = card_storage(f, card);
    assert_eq!(
        storage.intrinsic_drain_charge,
        p.intrinsic_drain_charge_curves
    );
    for (read, authored) in [
        (
            storage.gate_source_capacitance_f,
            p.gate_source_capacitance_f,
        ),
        (storage.gate_drain_capacitance_f, p.gate_drain_capacitance_f),
    ] {
        assert_eq!(read.map(f64::to_bits), authored.map(f64::to_bits));
    }
    let rails = p
        .gate_drivers
        .map(|driver| [driver.on_voltage_v, driver.off_voltage_v]);
    let drivers = Drivers {
        cells: DRIVER,
        keys: f.driver_resistance.clone(),
        rails: rails.as_flattened().try_into().unwrap(),
        start: f.thermal_at(start_k),
    };
    Case {
        fixture: &f.base,
        context,
        storage_keys: &f.charge,
        storage,
        steady,
        start: initial(steady, p),
        dt_s,
        substeps: 1,
        load_cell: None,
        schedule: Vec::new(),
        linked: false,
        coolant: None,
        drivers: Some(leak(drivers)),
    }
}

/// The owner's step built by hand for a gate case (`held_to_owner`), with
/// the parameters `at` gives at the cells' end temperatures. The new state,
/// then the seven cells' heats and the four drivers'.
fn gate_owner<C: GateChannelResponse + Clone + 'static>(
    case: &Case<'_, GateStorage>,
    at: impl Fn(&[f64]) -> GateCommutationParameters<C>,
) -> impl Fn(usize, &[f64; 15], &[f64]) -> ([f64; 15], Vec<f64>) {
    move |k: usize, old: &[f64; 15], kelvin: &[f64]| {
        let p = at(kelvin);
        let intervals = [CommutationInterval {
            gates: case.gates(k),
            duration_s: case.dt_s,
        }];
        let (e, c) = (power_controls(), cycle_controls());
        let owner = PreparedGateBatteryMotorCycle::<C>::new_gate_controlled(
            owner_drive(),
            p,
            &intervals,
            case.substeps,
            e,
            c,
            || false,
        );
        let mut next = [0.0; 15];
        let step = owner
            .unwrap()
            .step_interval_commands_controlled(old, &[case.steady.command], &mut next, || false)
            .unwrap();
        let mut heat = per_body(step.heat_joules(), 0.0).to_vec();
        heat.extend_from_slice(&step.all_heat_joules()[14..]);
        (next, heat)
    }
}

/// The module's authored parameters with the card's `laws` (Cgd, whether
/// the channels feed back, the on resistance's slope) in closed form at the
/// cells' temperatures. Each channel is `g` at its package's temperature on
/// the module's three gate voltages, each gate's on resistance 0.25 ohm +
/// slope (T - 300 K) at its driver cell.
fn closed_form(
    (miller, feedback, slope): (f64, bool, f64),
) -> impl Fn(&[f64]) -> GateCommutationParameters {
    move |kelvin: &[f64]| {
        let mut p = authored(miller, feedback);
        p.channel_conductance = core::array::from_fn(|i| {
            let knot = |volts| (volts, g(i, volts, kelvin[i + 1], feedback));
            GateConductanceCurve::try_new([-4.0, 0.0, 4.0].map(knot).to_vec()).unwrap()
        });
        for (driver, t) in p.gate_drivers.iter_mut().zip(&kelvin[IDS.len()..]) {
            driver.on_resistance_ohm += slope * (t - 300.0);
        }
        p
    }
}

#[test]
fn a_gate_graph_drives_its_gate_on_and_off_through_a_finite_driver() {
    // The module's finite-driver fixture: no Miller charge, the input high
    // gate's on resistance 0.25 ohm + 0.125 ohm/K (T - 300 K) at its driver
    // cell. A step of 1 s with the steady gates, then one with every gate
    // off.
    let laws = (0.0, false, 0.125);
    let f = leak(gate_fixture(0.0, false, 0.125));
    let steady = steady_drive(false);
    let p = authored(0.0, false);
    let off = DeviceGateState { on: [false; 4] };
    let case = Case {
        schedule: vec![steady.gates, off],
        ..gate_case(f, &steady, &p, 300.0, 1.0)
    };
    let graph: Vec<Row<_, [f64; CELLS]>> = kernel_rows(&case, 2);
    // The owner by hand, and its closed form.
    held_to_owner(
        &case,
        &graph,
        gate_owner(&case, closed_form(laws)),
        [1e-14, 5e-12],
    );
    // The module's closed form: the gate charges to 2 V through 0.5 ohm,
    // the on resistance at 302 K, its driver's cell taking 2 J; then from
    // the 1 V rail through 2 ohm it falls to 1.6 V and 0.32 J more.
    near(graph[0].state[8], 2.0);
    near(graph[0].kelvin[7], 302.0);
    near(graph[1].state[8], 1.6);
    near(graph[1].kelvin[7], 302.32);
    for kelvin in &graph[1].kelvin[8..] {
        near(*kelvin, 300.0);
    }
}

#[test]
fn a_gate_graph_carries_the_thermal_channel_surface_and_the_miller_charge() {
    // The module's Miller fixture: Cgd 0.05 F, each channel a G(Vgs, T)
    // surface rising 0.1 per volt and 0.05 per kelvin, the input high gate's
    // on resistance rising 0.125 ohm/K. Steps of 1/8 s.
    let laws = (0.05, true, 0.125);
    let f = leak(gate_fixture(0.05, true, 0.125));
    let steady = steady_drive(false);
    let p = authored(0.05, true);
    let case = gate_case(f, &steady, &p, 300.0, 0.125);
    let graph: Vec<Row<_, [f64; CELLS]>> = kernel_rows(&case, 3);
    // The owner by hand, its surface `g` on three knots where the card's
    // line has more, so not bit for bit.
    held_to_owner(
        &case,
        &graph,
        gate_owner(&case, closed_form(laws)),
        [2e-13, 1e-12],
    );
    // The Miller charge is real, as the module checked: a gate-drain
    // voltage moves by more than 1 uV in the first step.
    let (v0, g0) = voltages(&case.start, 0.05);
    let (v1, g1) = voltages(&graph[0].state, 0.05);
    let moved = |i: usize| ((g1[i] - v1[i]) - (g0[i] - v0[i])).abs();
    assert!((0..4).any(|i| moved(i) > 1e-6));
}

#[test]
fn a_gate_graph_shares_a_finite_liquid_coolant_with_a_driver_cell() {
    // The module's driver-cooling fixture: the fixture's claims as a liquid
    // (the nine-state coolant's card), 2 kg at 300 K, joined to the input
    // high gate's driver cell at 0.5 W/K. The on resistance rises
    // 0.175 ohm/K; steps of 1 s.
    let laws = (0.0, false, 0.175);
    let f = leak(gate_fixture(0.0, false, 0.175));
    let steady = steady_drive(false);
    let p = authored(0.0, false);
    let base = coolant(&f.base);
    let liquid = leak(Coolant {
        context: common::point(300.0),
        cell: DepositedCell {
            id: CellId(120),
            segment_id: 11,
            body_id: 11,
            specific_h_j_kg: 1.0,
            geometry: CellGeometry {
                start_m: [30.0, 0.0, 0.0],
                end_m: [32.0, 0.0, 0.0],
                ..base.cell.geometry
            },
            ..base.cell
        },
        partner: (DRIVER[0], 0.5),
        ..base
    });
    let case = Case {
        coolant: Some(liquid),
        ..gate_case(f, &steady, &p, 300.0, 1.0)
    };
    let graph: Vec<Row<_, [f64; CELLS + 1]>> = kernel_rows(&case, 2);
    // The owner by hand, the coolant's energy by its own caloric curve, and
    // the closed form.
    held_to_owner(
        &case,
        &graph,
        gate_owner(&case, closed_form(laws)),
        [2e-13, 5e-12],
    );
    // The module's closed form for the first step, backward Euler from
    // 300 K through the link with cells of 1 and 2 J/K: the driver's 2 J
    // warm its cell 10/7 K and the coolant 2/7 K, where the on resistance
    // is 0.5 ohm and the gate reaches 2 V.
    near(graph[0].kelvin[7], 300.0 + 10.0 / 7.0);
    near(graph[0].kelvin[CELLS], 300.0 + 2.0 / 7.0);
    near(graph[0].state[8], 2.0);
}

/// The module's drive, bound at one accepted state, refused a step from
/// another (`gate_drive_bound_at_one_state_refuses_a_step_from_another_accepted_state`
/// and `gate_conductance_bound_at_one_state_refuses_a_step_from_another_accepted_state`).
/// The graph binds nothing: from that other state it takes the step at the
/// laws' answers there.
#[test]
fn a_gate_graph_steps_from_another_accepted_state_at_the_answers_there() {
    // The two refusals' fixtures from 320 K: the on resistance rising
    // 0.125 ohm/K (1 s), then only the channels' surface varying, with the
    // Miller charge (1/8 s).
    let fixtures = [((0.0, false, 0.125), 1.0), ((0.05, true, 0.0), 0.125)];
    for (laws @ (miller, feedback, slope), dt_s) in fixtures {
        let f = leak(gate_fixture(miller, feedback, slope));
        let steady = steady_drive(false);
        let p = authored(miller, feedback);
        let case = gate_case(f, &steady, &p, 320.0, dt_s);
        let graph: Vec<Row<_, [f64; CELLS]>> = kernel_rows(&case, 1);
        // The owner by hand at the laws' closed forms there, its new state
        // bit for bit.
        let owner = gate_owner(&case, closed_form(laws));
        let (next, _) = owner(0, &case.start, &graph[0].kelvin);
        assert_eq!(next.map(f64::to_bits), graph[0].state.map(f64::to_bits));
        held_to_owner(&case, &graph, owner, [0.0, 5e-13]);
        // The answers at the step's end differ in their bits from those at
        // 300 K, the state the module's drive was bound at.
        let end: [f64; CELLS] = graph[0].kelvin;
        assert!(end.iter().all(|&kelvin| kelvin >= 320.0), "{end:?}");
        let card = Card {
            material: &f.base.card,
            curve: &f.base.caloric,
            context: case.context,
        };
        let bits = |values: Vec<f64>| values.into_iter().map(f64::to_bits).collect::<Vec<_>>();
        if feedback {
            for (i, key) in f.conductance.iter().enumerate() {
                let law = CardChannel { card, key, tag: 0 };
                let line = |kelvin| {
                    law.knots(kelvin)
                        .unwrap()
                        .into_iter()
                        .flat_map(<[f64; 2]>::from)
                };
                assert_ne!(
                    bits(line(end[i + 1]).collect()),
                    bits(line(300.0).collect())
                );
            }
        } else {
            let wanted = [(
                &f.driver_resistance[0],
                ScalarAdmissibility::StrictlyPositive,
            )];
            let ohms = |kelvin| card.values(kelvin, &wanted).unwrap();
            assert_ne!(bits(ohms(end[7])), bits(ohms(300.0)));
        }
    }
}
