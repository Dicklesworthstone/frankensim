//! The module's battery fixtures as declared graphs on the composition kernel.
//!
//! `circuit_enthalpy::battery` couples fs-phs's equivalent-circuit battery to
//! the enthalpy timeline by hand: backward Euler for the cells, the battery's
//! discrete-gradient step with each resistance at the material resistance of
//! its receiving cell's END temperature, Newton on the cell energies. The
//! graph says the same thing with no code for the pair:
//!
//! - `fs_system::battery::EquivalentCircuitBattery`: the same curve and
//!   parameters, each resistance following its receiving cell's card through
//!   [`CardResistance`] when the module resolves it from the card;
//! - an `elements::FlowSource` drawing the module's terminal current;
//! - the timeline of `kernel`, whose cell ports return end-of-step
//!   temperatures;
//! - one heat route per loss channel into the module's receiving cell.
//!
//! The two sides solve the same equations with different solvers, so they
//! agree only to within the solves' tolerances.

use std::sync::Arc;

use fs_couple::cavity_phs::piston::thermal::source::ThermalMaterialPoint;
use fs_couple::circuit_enthalpy::battery::{
    BatteryResistance, BatteryThermalChannel, BatteryThermalNetwork,
};
use fs_matdb::{MaterialCard, PropertyKey};
use fs_material::phase::ResolvedHeatCapacityCurve;
use fs_material::state_point::{
    MaterialPropertySelection, ResolvedScalarProperty, ScalarAdmissibility,
    ScalarPropertyRequirement, THERMAL_CONDUCTIVITY_DIMS, resolve_material_state_point,
};
use fs_ports::PortKind;
use fs_system::battery::{BatteryLoss, EquivalentCircuitBattery};
use fs_system::circuit::ResistanceLaw;
use fs_system::elements::{FlowSource, SourceAccount, Waveform};
use fs_system::{NodeId, RouteTemperature, System, SystemBuilder, heat_sink};
use fs_thermal::component::TimelineState;
use fs_thermal::enthalpy::CellId;

use super::{
    IDS, LINK, R, bindings, controls, electrical, key, materials, native, parameters, point,
};
use crate::common;
use crate::kernel::{CELLS, enthalpy_at, timeline, with_cx};

/// A loss channel's resistance, resolved as the module resolves it: the
/// receiving cell's card at the phase state of the temperature the channel's
/// heat port receives.
#[derive(Debug)]
struct CardResistance {
    card: &'static MaterialCard,
    curve: &'static ResolvedHeatCapacityCurve,
    key: &'static PropertyKey,
    /// Names the card in the battery's identity: the fixture's index.
    tag: u64,
}

impl ResistanceLaw for CardResistance {
    fn ohms(&self, kelvin: f64) -> Result<f64, String> {
        let h = enthalpy_at(self.curve, kelvin)?;
        let phase = self
            .curve
            .curve()
            .state_at_specific_enthalpy(h)
            .map_err(|error| format!("{error:?}"))?;
        let context = point(300.0);
        let material = ThermalMaterialPoint {
            material: self.card,
            context: &context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        };
        let query = material.at(phase).map_err(|error| format!("{error:?}"))?;
        let requirement = ScalarPropertyRequirement::try_with_key(
            self.key,
            ScalarAdmissibility::StrictlyPositive,
        )
        .map_err(|error| format!("{error:?}"))?;
        let snapshot = resolve_material_state_point(
            self.card,
            &query,
            &[requirement],
            MaterialPropertySelection::SingleClaimOnly,
        )
        .map_err(|error| format!("{error:?}"))?;
        snapshot
            .property_by_key(self.key)
            .map(ResolvedScalarProperty::value_si)
            .ok_or_else(|| "the card has no battery resistance".to_owned())
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        words.push(self.tag);
    }
}

/// One of the module's battery fixtures.
struct Fixture {
    capacity_c: f64,
    knots: &'static [(f64, f64)],
    /// State of charge, polarization voltage \[V\] and bus voltage \[V\].
    start: [f64; 3],
    /// Terminal current \[A\], positive charging.
    current_a: f64,
    /// The cells' specific heat \[J/(kg K)\].
    cp: [f64; 2],
    /// The cards' resistance slope \[1/K\].
    alpha: f64,
    /// Whether both resistances come from the receiving cells' cards.
    from_cards: bool,
    /// The receiving cell of the series, then the polarization, heat.
    receivers: [CellId; 2],
}

/// One accepted interval: the battery's charges \[C\], the two cells'
/// temperatures \[K\] and the heat of both resistances \[J\].
#[derive(Debug, Clone)]
struct Row {
    state: [f64; 3],
    kelvin: [f64; 2],
    joule_j: f64,
}

type Materials = (
    &'static [MaterialCard; 2],
    &'static [ResolvedHeatCapacityCurve; 2],
);

fn leaked(fixture: &Fixture) -> Materials {
    let (cards, curves) = materials(fixture.cp, fixture.alpha, 2000.0);
    (Box::leak(Box::new(cards)), Box::leak(Box::new(curves)))
}

/// The module's run: `steps` intervals of `dt_s` with `step_candidate`.
fn module_run(fixture: &Fixture, (cards, curves): Materials, dt_s: f64, steps: usize) -> Vec<Row> {
    let context = point(300.0);
    let thermal_key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = bindings(cards, curves, &context, &thermal_key);
    let resistance = key();
    let channels: [BatteryThermalChannel<'_>; 2] =
        core::array::from_fn(|i| BatteryThermalChannel {
            cell: fixture.receivers[i],
            resistance: if fixture.from_cards {
                BatteryResistance::Source {
                    key: &resistance,
                    selection: MaterialPropertySelection::SingleClaimOnly,
                }
            } else {
                BatteryResistance::Fixed
            },
        });
    let network = BatteryThermalNetwork {
        cells: &cells,
        links: &LINK,
        channels: &channels,
    };
    let native = native(fixture.capacity_c, fixture.knots);
    let [soc, vp, vb] = fixture.start;
    let mut thermal = common::initial(curves, 0.0, None);
    let mut battery = network
        .bind(
            &thermal,
            &native,
            native.state(soc, vp, vb).unwrap(),
            || false,
        )
        .unwrap();
    let mut rows = Vec::new();
    for k in 1..=steps {
        let next = network
            .step_candidate(
                &thermal,
                k as f64 * dt_s,
                &battery,
                fixture.current_a,
                electrical(),
                16,
                controls(),
                &[0.0; 2],
                || false,
            )
            .unwrap();
        rows.push(Row {
            state: *next.battery().state(),
            kelvin: CELLS.map(|id| next.thermal().temperature_k(id).unwrap()),
            joule_j: next.report().deposited_joule_j,
        });
        let (t, b, _) = next.into_parts();
        (thermal, battery) = (t, b);
    }
    rows
}

/// The same fixture as a graph: the battery with its laws, a current source
/// on its terminal, the cells of `kernel::timeline`, one heat route per loss
/// channel. Returns the system, the battery and the cells.
fn graph(fixture: &Fixture, (cards, curves): Materials, dt_s: f64) -> (System, NodeId, NodeId) {
    let (cells, cells_start) = timeline(curves, 0.0, [0.0; 2], dt_s);
    let mut battery = EquivalentCircuitBattery::new(
        fixture.capacity_c,
        fixture.knots,
        parameters(),
        electrical(),
        16,
        1e-8,
    )
    .unwrap();
    if fixture.from_cards {
        let key: &'static PropertyKey = Box::leak(Box::new(key()));
        for (loss, receiver) in [BatteryLoss::Series, BatteryLoss::Polarization]
            .into_iter()
            .zip(fixture.receivers)
        {
            let index = CELLS.iter().position(|&cell| cell == receiver).unwrap();
            let law = CardResistance {
                card: &cards[index],
                curve: &curves[index],
                key,
                tag: index as u64,
            };
            battery = battery.following(loss, Arc::new(law)).unwrap();
        }
    }
    let [soc, vp, vb] = fixture.start;
    let start = battery.state(soc, vp, vb).unwrap();

    let mut builder = SystemBuilder::new();
    let electrical = builder.add("battery", battery, start).unwrap();
    // Positive flow enters the source, so it draws minus the terminal current.
    let load = builder
        .add(
            "load",
            FlowSource::in_domain(
                PortKind::ElectricalVoltageCurrent,
                Waveform::constant(-fixture.current_a).unwrap(),
            ),
            SourceAccount::default(),
        )
        .unwrap();
    builder
        .common_flow(
            builder.port(load, "p").unwrap(),
            &[builder.port(electrical, "terminal").unwrap()],
        )
        .unwrap();
    let cells = builder.add("cells", cells, cells_start).unwrap();
    for (heat, receiver) in ["heat_series", "heat_polarization"]
        .into_iter()
        .zip(fixture.receivers)
    {
        builder
            .heat_route(
                builder.port(electrical, heat).unwrap(),
                &[(
                    builder.port(cells, &format!("cell{}", receiver.0)).unwrap(),
                    1.0,
                )],
                None,
                RouteTemperature::SinkEquivalent,
            )
            .unwrap();
    }
    // A cell that receives no loss is insulated, as in the module.
    for cell in CELLS {
        if !fixture.receivers.contains(&cell) {
            builder
                .boundary(
                    builder.port(cells, &format!("cell{}", cell.0)).unwrap(),
                    Box::new(|_, _, values| values[0] = 0.0),
                )
                .unwrap();
        }
    }
    builder
        .heat_route(
            builder.port(cells, "environment").unwrap(),
            &[],
            Some(heat_sink("environment", 300.0).unwrap()),
            RouteTemperature::Fixed(300.0),
        )
        .unwrap();
    let mut system = builder.seal().unwrap();
    for heat in ["heat_series", "heat_polarization"] {
        let port = system.port(electrical, heat).unwrap();
        assert!(system.guess(port, &[300.0]));
    }
    (system, electrical, cells)
}

fn graph_run(
    system: &mut System,
    electrical: NodeId,
    cells: NodeId,
    dt_s: f64,
    steps: usize,
) -> Vec<Row> {
    let mut rows = Vec::new();
    let mut joule_before = 0.0;
    with_cx(|cx| {
        for _ in 0..steps {
            let report = system.step(dt_s, cx).unwrap();
            assert_eq!(
                report.dt_s.to_bits(),
                dt_s.to_bits(),
                "the kernel took the whole interval"
            );
            let joule = system.ledger().nodes[electrical.index()]
                .dissipated_j
                .value();
            let timeline = system.state::<TimelineState<'static>>(cells).unwrap();
            rows.push(Row {
                state: *system.state::<[f64; 3]>(electrical).unwrap(),
                kelvin: CELLS.map(|id| timeline.simulation().temperature_k(id).unwrap()),
                joule_j: joule - joule_before,
            });
            joule_before = joule;
        }
    });
    rows
}

/// Bounds of one comparison: temperature \[K\], charge \[C\] and heat per
/// interval \[J\].
struct Bounds {
    kelvin: f64,
    charge_c: f64,
    joule_j: f64,
}

/// Runs the fixture both ways and holds every interval to `bounds`. Returns
/// the module's rows. The kernel's ledger closes, and the heat that left the
/// battery is in the cells.
fn compare(fixture: &Fixture, dt_s: f64, steps: usize, bounds: &Bounds) -> Vec<Row> {
    let materials = leaked(fixture);
    let module = module_run(fixture, materials, dt_s, steps);
    let (mut system, electrical, cells) = graph(fixture, materials, dt_s);
    let kernel = graph_run(&mut system, electrical, cells, dt_s, steps);
    assert_eq!((kernel.len(), module.len()), (steps, steps));
    for (k, (ours, theirs)) in kernel.iter().zip(&module).enumerate() {
        for cell in 0..2 {
            let (x, y) = (ours.kelvin[cell], theirs.kelvin[cell]);
            assert!(
                (x - y).abs() <= bounds.kelvin,
                "interval {k}, cell {cell}: {x} K vs {y} K"
            );
        }
        for (i, (x, y)) in ours.state.iter().zip(&theirs.state).enumerate() {
            assert!(
                (x - y).abs() <= bounds.charge_c,
                "interval {k}, charge {i}: {x} vs {y}"
            );
        }
        assert!(
            (ours.joule_j - theirs.joule_j).abs() <= bounds.joule_j,
            "interval {k}, heat: {} vs {}",
            ours.joule_j,
            theirs.joule_j
        );
    }
    let ledger = system.ledger();
    assert!(
        ledger.energy_residual_j().abs() <= 1e-8,
        "{}",
        ledger.energy_residual_j()
    );
    let heat = ledger.nodes[electrical.index()].dissipated_j.value();
    let stored = ledger.nodes[cells.index()].stored_change_j.value();
    assert!(
        heat > 0.0 && (stored - heat).abs() <= 1e-8,
        "{stored} vs {heat}"
    );
    module
}

#[test]
fn the_declared_graph_crosses_the_ocv_knot_as_the_module_does() {
    // The module's knot fixture: no terminal current, the bus charging from
    // the cell across the 0.5 SOC knot in one step of 30 ms, each loss
    // heating its own cell.
    let fixture = Fixture {
        capacity_c: 1.0,
        knots: &[(0.0, 3.0), (0.5, 3.0), (1.0, 4.2)],
        start: [0.6, 0.2, 0.4],
        current_a: 0.0,
        cp: [2.0, 3.0],
        alpha: 0.0,
        from_cards: false,
        receivers: IDS,
    };
    // On the export of 1e611fdf33 the two agree bit for bit. The bounds
    // allow for rounding on another platform.
    let bounds = Bounds {
        kelvin: 1e-11,
        charge_c: 1e-13,
        joule_j: 1e-13,
    };
    let module = compare(&fixture, 0.03, 1, &bounds);
    assert!(module[0].state[0] < 0.5, "{module:?}");
}

#[test]
fn the_declared_graph_reproduces_the_endpoint_material_resistances() {
    // The module's card fixture: both resistances rise 0.4 % per kelvin of
    // their cell's end temperature, and the small heat capacities warm by
    // kelvins an interval; 1 A drawn for four intervals of 0.2 s.
    let fixture = Fixture {
        capacity_c: 10.0,
        knots: &[(0.0, 4.0), (1.0, 4.0)],
        start: [0.8, 0.3, 3.0],
        current_a: -1.0,
        cp: [0.02, 0.03],
        alpha: 0.004,
        from_cards: true,
        receivers: IDS,
    };
    // On the export of 1e611fdf33 the first interval agrees bit for bit;
    // after that the two agree within 3.9e-12 K, 4.7e-15 C and 5.6e-16 J.
    // The bounds sit about fifty times above that.
    let bounds = Bounds {
        kelvin: 2e-10,
        charge_c: 2.5e-13,
        joule_j: 3e-14,
    };
    let module = compare(&fixture, 0.2, 4, &bounds);
    // The feedback is real: cell 30 warms by more than 15 K, so its
    // resistance rises by more than 6 %.
    assert!(module[3].kelvin[0] > 315.0, "{module:?}");
}

#[test]
fn the_declared_graph_routes_both_losses_to_one_receiver_in_discharge_and_regeneration() {
    // The module's common-receiver fixture: a flat curve at its stationary
    // state, both losses into cell 30, for both signs of the current.
    for current in [-1.0, 1.0] {
        let fixture = Fixture {
            capacity_c: 10.0,
            knots: &[(0.0, 4.0), (1.0, 4.0)],
            start: [0.5, -current * R[1], 4.0 + current * (R[0] + R[1])],
            current_a: current,
            cp: [2.0, 3.0],
            alpha: 0.0,
            from_cards: false,
            receivers: [IDS[0]; 2],
        };
        // On the export of 1e611fdf33 the two agree bit for bit but for one
        // rounding of cell 30's temperature (5.7e-14 K) in the last interval.
        let bounds = Bounds {
            kelvin: 1e-11,
            charge_c: 1e-13,
            joule_j: 1e-13,
        };
        let module = compare(&fixture, 0.2, 3, &bounds);
        assert!((module[2].kelvin[1] - 300.0).abs() < 1e-12, "{module:?}");
    }
}
