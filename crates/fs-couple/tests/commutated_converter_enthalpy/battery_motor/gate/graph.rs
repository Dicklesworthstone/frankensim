//! The gate fixtures as declared graphs (`battery_motor/graph.rs`).
//!
//! The drive runs in the `GateStorage` lane. Its storage is the package
//! card's, read once ([`card_storage`]), since the module holds it constant.
//! Each switch's channel follows the card's G(Vgs, T) line at its package
//! ([`CardChannel`]); each gate's on and off resistances follow the card at
//! its driver cell. The four driver cells join the timeline after the seven,
//! each receiving its gate's driver heat, and eight `Force` sources hold the
//! rails at the module's driver voltages.

use std::sync::Arc;

use fs_couple::circuit_enthalpy::converter::commutation::GateChannelBinding;
use fs_matdb::{ClaimSelection, MaterialAnswer};
use fs_material::state_point::{
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_phs::circuit::converter::commutation::gate::{GateChannelLaw, GateChannelResponse};
use fs_system::commutated_drive::{
    ChannelLaw, CommutatedBatteryMotorDrive, GateDriverResistances, GateRail, GateStorage,
    SwitchPosition,
};

use super::super::graph::{
    Card, CardResistance, Case, Coolant, Drivers, Lane, Row, cell_ids, cells_start, compare,
    heat_assembly, held_to_owner, kernel_run, leak, owner_drive, text,
};
use super::super::nine_state::coolant;
use super::*;

/// The timeline's cells with ports: the seven, then the drivers'.
pub(super) const CELLS: usize = 11;

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

/// A switch's channel from its package's card, as the module resolves it:
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
/// requires the storage constant, so one reading stands for every step. The
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
/// started as the module starts it, every cell at `start_k`, one step of
/// `dt_s` per program. The card's storage must be the authored storage, as
/// the module's bind requires.
pub(super) fn gate_case<'s>(
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

/// The temperatures of the case's timeline cells in `thermal`, in
/// `cell_ids` order.
fn kelvin<const N: usize>(
    case: &Case<'_, GateStorage>,
    thermal: &DepositionSimulation<'_>,
) -> [f64; N] {
    let ids = cell_ids(case);
    assert_eq!(ids.len(), N);
    core::array::from_fn(|i| thermal.temperature_k(ids[i]).unwrap())
}

/// The module's `steps` steps of the case from its drive bound at the case's
/// start, authored as `p`.
fn module_run(
    f: &GateFixture,
    case: &Case<'_, GateStorage>,
    p: &GateCommutationParameters,
    steps: usize,
) -> Vec<Row<[f64; 15], [f64; CELLS]>> {
    let context = gate_context();
    let cells = f.cells(&context);
    let (resistors, devices) = (drive_resistors(&f.base), f.base.devices());
    let (q, cgs, cgd) = (f.charges(), f.caps(false), f.caps(true));
    let (conductance, driver_r) = (f.conductances(), f.driver_resistors());
    let bindings = gate_bindings(&q, &cgs, &cgd, &conductance, &driver_r);
    module_steps(
        &coupled(&cells, &resistors, &devices, bindings),
        case,
        p,
        steps,
    )
}

/// The module's `steps` steps of the case on `network`, its drive bound at
/// the case's start, authored as `p`.
pub(super) fn module_steps<B: GateChannelBinding>(
    network: &GateBatteryMotorThermalNetwork<'_, B>,
    case: &Case<'_, GateStorage>,
    p: &GateCommutationParameters,
    steps: usize,
) -> Vec<Row<[f64; 15], [f64; CELLS]>> {
    let mut thermal = cells_start(case);
    let native = &case.steady.native;
    let mut drive = network
        .bind(&thermal, native, p, case.start, || false)
        .unwrap();
    let mut rows = Vec::new();
    for k in 0..steps {
        let intervals = [CommutationInterval {
            gates: case.gates(k),
            duration_s: case.dt_s,
        }];
        let next = network
            .step_candidate(
                &thermal,
                thermal.time_s() + case.dt_s,
                &drive,
                &intervals,
                &[case.steady.command],
                case.substeps,
                case.load_cell,
                power_controls(),
                cycle_controls(),
                drive_controls(CELLS),
                &[0.0; CELLS],
                || false,
            )
            .unwrap();
        rows.push(Row {
            state: *next.drive().state(),
            kelvin: kelvin(case, next.thermal()),
        });
        (thermal, drive, _) = next.into_parts();
    }
    rows
}

/// The module's `steps` steps of the case through its heat assembly, which
/// carries the case's coolant and link (`FixedCooling`).
fn module_linked_run(
    f: &GateFixture,
    case: &Case<'_, GateStorage>,
    p: &GateCommutationParameters,
    steps: usize,
) -> Vec<Row<[f64; 15], [f64; CELLS + 1]>> {
    let start = &cells_start(case);
    heat_assembly(case, start, |cx, assembly| {
        let (resistors, devices) = (drive_resistors(&f.base), f.base.devices());
        let (q, cgs, cgd) = (f.charges(), f.caps(false), f.caps(true));
        let (conductance, driver_r) = (f.conductances(), f.driver_resistors());
        let network = GateBatteryMotorHeatNetwork {
            drive: CommutatedBatteryMotorHeatNetwork {
                assembly,
                resistors: &resistors,
                devices: &devices,
                bearing_cell: IDS[0],
            },
            gate: gate_bindings(&q, &cgs, &cgd, &conductance, &driver_r),
            maximum_cells_per_chord: 4096,
        };
        let native = &case.steady.native;
        let mut drive = network.bind(cx, start, native, p, case.start).unwrap();
        let mut thermal = start.clone();
        let mut rows = Vec::new();
        for k in 0..steps {
            let intervals = [CommutationInterval {
                gates: case.gates(k),
                duration_s: case.dt_s,
            }];
            let next = network
                .step_candidate(
                    cx,
                    &thermal,
                    thermal.time_s() + case.dt_s,
                    &drive,
                    &intervals,
                    &[case.steady.command],
                    case.substeps,
                    case.load_cell,
                    power_controls(),
                    cycle_controls(),
                    drive_controls(CELLS + 1),
                    &[0.0; CELLS + 1],
                )
                .unwrap();
            rows.push(Row {
                state: *next.drive().state(),
                kelvin: kelvin(case, next.thermal()),
            });
            (thermal, drive, _) = next.into_parts();
        }
        rows
    })
}

/// The owner's step built by hand for a gate case (`held_to_owner`), with
/// the parameters `at` gives at the cells' end temperatures. The new state,
/// then the seven cells' heats and the four drivers'.
pub(super) fn gate_owner<C: GateChannelResponse + Clone + 'static>(
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
    let f = leak(gate_fixture(0.0, false, 0.125, false, false));
    let steady = steady_drive(false);
    let p = authored(0.0, false);
    let off = DeviceGateState { on: [false; 4] };
    let case = Case {
        schedule: vec![steady.gates, off],
        ..gate_case(f, &steady, &p, 300.0, 1.0)
    };
    let module = module_run(f, &case, &p, 2);
    let graph = kernel_run(&case, 2);
    compare(&module, &graph, [5e-13, 5e-12]);
    // Without the module: the owner by hand, and its closed form.
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
    let f = leak(gate_fixture(0.05, true, 0.125, false, false));
    let steady = steady_drive(false);
    let p = authored(0.05, true);
    let case = gate_case(f, &steady, &p, 300.0, 0.125);
    let module = module_run(f, &case, &p, 3);
    let graph = kernel_run(&case, 3);
    compare(&module, &graph, [5e-13, 5e-12]);
    // Without the module: the owner by hand, its surface `g` on three
    // knots where the card's line has more, so not bit for bit.
    held_to_owner(
        &case,
        &graph,
        gate_owner(&case, closed_form(laws)),
        [2e-13, 1e-12],
    );
    // The Miller charge is real, as the module checks: a gate-drain
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
    let f = leak(gate_fixture(0.0, false, 0.175, false, false));
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
    let module = module_linked_run(f, &case, &p, 2);
    let graph = kernel_run(&case, 2);
    compare(&module, &graph, [5e-13, 5e-12]);
    // Without the module: the owner by hand, the coolant's energy by its
    // own caloric curve, and the closed form.
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

/// The module's drive bound at one accepted state refuses a step from
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
        let f = leak(gate_fixture(miller, feedback, slope, false, false));
        let steady = steady_drive(false);
        let p = authored(miller, feedback);
        let case = gate_case(f, &steady, &p, 320.0, dt_s);
        let graph = kernel_run(&case, 1);
        // The module bound at that state takes the same step.
        let module = module_run(f, &case, &p, 1);
        compare(&module, &graph, [1e-14, 1e-12]);
        // Without the module: the owner by hand at the laws' closed forms
        // there, its new state bit for bit.
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
