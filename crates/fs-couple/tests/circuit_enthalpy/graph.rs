//! The module's coupling as a declared graph on the composition kernel.
//!
//! `circuit_enthalpy::step_candidate` couples fs-phs's Kirchhoff circuit to
//! the enthalpy timeline by hand: backward Euler for the cells, the circuit's
//! midpoint step with each resistor at the material resistance of its cell's
//! END temperature, Newton on the cell energies. The graph says the same
//! thing with no code for the pair:
//!
//! - `fs_system::circuit::KirchhoffCircuit`: the same netlist, one heat port
//!   per resistor, each resistor following its card through
//!   [`CardResistance`] at the temperature its port receives;
//! - `fs_thermal::component::ThermalTimeline::implicit` with backward Euler,
//!   returning end-of-step temperatures, so that temperature is the cell's
//!   end temperature;
//! - one heat route per resistor into the module's destination cell.
//!
//! The two sides solve the same equations with different solvers (the
//! module's Newton--Krylov on cell energies and its own backward-Euler
//! transport; the kernel's interface iteration and the implicit conduction
//! graph), so they agree to the solves' tolerances, not bit for bit.

use std::sync::Arc;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_couple::cavity_phs::piston::thermal::power::joule::MaterialOhmicBody;
use fs_couple::cavity_phs::piston::thermal::source::ThermalMaterialPoint;
use fs_couple::circuit_enthalpy::step_candidate;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::MaterialCard;
use fs_material::phase::ResolvedHeatCapacityCurve;
use fs_material::state_point::MaterialPropertySelection;
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_system::circuit::{KirchhoffCircuit, ResistanceLaw};
use fs_system::{NodeId, RouteTemperature, System, SystemBuilder, heat_sink};
use fs_thermal::component::{ImplicitScheme, ThermalTimeline, TimelineState};
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionControls, ImplicitConductionGraph, ImplicitConductionLaw,
    ImplicitConductionLinearControls, ImplicitConductionState,
};
use fs_thermal::enthalpy::conduction::{ConductionObserver, ConductiveLink, ConductivityLaw};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation, MaterialId, ThermalLaw};

use crate::common::*;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let key = StreamKey {
            seed: 11,
            kernel_id: 29,
            tile: 0,
            iteration: 0,
        };
        f(&Cx::new(
            &gate,
            arena,
            key,
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

/// The specific enthalpy [J/kg] at `kelvin` on a caloric curve whose knots
/// rise in temperature, by the curve's own linear interpolation.
fn enthalpy_at(curve: &ResolvedHeatCapacityCurve, kelvin: f64) -> Result<f64, String> {
    let knots = curve.curve().knots();
    for pair in knots.windows(2) {
        let (low, high) = (&pair[0], &pair[1]);
        if (low.temperature_k..=high.temperature_k).contains(&kelvin)
            && high.temperature_k > low.temperature_k
        {
            let fraction = (kelvin - low.temperature_k) / (high.temperature_k - low.temperature_k);
            return Ok(low.specific_enthalpy_j_kg
                + fraction * (high.specific_enthalpy_j_kg - low.specific_enthalpy_j_kg));
        }
    }
    Err(format!("{kelvin} K is outside the caloric curve"))
}

/// The fixture's resistor, resolved as the module's `coefficients` resolves
/// it: the card's conductor at the phase state of the cell, here the state
/// at the temperature the resistor's heat port receives.
#[derive(Debug)]
struct CardResistance {
    card: &'static MaterialCard,
    curve: &'static ResolvedHeatCapacityCurve,
    /// Names the card in the circuit's identity: the fixture's index.
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
        let body = MaterialOhmicBody {
            material: ThermalMaterialPoint {
                material: self.card,
                context: &context,
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            selection: MaterialPropertySelection::SingleClaimOnly,
            length_m: 1.0,
            area_m2: 1.0,
        };
        body.resolve_conductor(phase, 1.0)
            .map(|conductor| conductor.resistance_ohm())
            .map_err(|error| format!("{error:?}"))
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        words.push(self.tag);
    }
}

/// The fixture's caloric curve with the card's constant conductivity `k`,
/// for the implicit conduction graph. The temperature is the curve's own.
struct CardConduction {
    curve: &'static ResolvedHeatCapacityCurve,
    k_w_m_k: f64,
}

impl ConductivityLaw for CardConduction {
    fn conductivity_w_m_k(&self, _: f64) -> Result<f64, DepositionError> {
        Ok(self.k_w_m_k)
    }
}

impl ImplicitConductionLaw for CardConduction {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self.curve.curve()
    }

    fn conductivity_law(&self) -> &dyn ConductivityLaw {
        self
    }

    fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        let knots = self.curve.curve().knots();
        [
            knots[0].specific_enthalpy_j_kg,
            knots[knots.len() - 1].specific_enthalpy_j_kg,
        ]
    }

    fn state_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<ImplicitConductionState, DepositionError> {
        let knots = self.curve.curve().knots();
        let (first, last) = (&knots[0], &knots[knots.len() - 1]);
        // The fixture's curve is one linear segment: dT/dh is its inverse slope.
        let slope = (last.temperature_k - first.temperature_k)
            / (last.specific_enthalpy_j_kg - first.specific_enthalpy_j_kg);
        Ok(ImplicitConductionState {
            temperature_k: self.curve.curve().temperature_k(h)?,
            temperature_slope_kg_k_j: slope,
            conductivity_w_m_k: self.k_w_m_k,
            conductivity_slope_kg_w_m_k_j: 0.0,
        })
    }
}

/// One accepted interval: the capacitor charge [C], the two cells'
/// temperatures [K] and the Joule heat [J].
#[derive(Debug, Clone, Copy)]
struct Row {
    charge_c: f64,
    kelvin: [f64; 2],
    joule_j: f64,
}

const CELLS: [CellId; 2] = [CellId(30), CellId(10)];

/// The module's run: `steps` intervals of `dt_s` with `step_candidate`.
fn module_run(
    cards: &[MaterialCard; 2],
    curves: &[ResolvedHeatCapacityCurve; 2],
    conductance_w_k: f64,
    dt_s: f64,
    steps: usize,
) -> Vec<Row> {
    let mut thermal = initial(curves, conductance_w_k, None);
    let mut electrical = circuit(&thermal, cards, curves);
    let mut rows = Vec::new();
    for k in 1..=steps {
        let next = step_candidate(
            &thermal,
            k as f64 * dt_s,
            input(&electrical),
            config(),
            &[0.0; 2],
            |s| coefficients(s, cards, curves),
            || false,
        )
        .unwrap();
        rows.push(Row {
            charge_c: next.circuit().charges_c()[0],
            kelvin: CELLS.map(|id| next.thermal().temperature_k(id).unwrap()),
            joule_j: next.report().deposited_joule_j,
        });
        let (t, c, _) = next.into_parts();
        thermal = t;
        electrical = c;
    }
    rows
}

/// The fixture's two cells as one backward-Euler timeline whose cell ports
/// return end-of-step temperatures, and its start. The link's conductance
/// comes from the cards' conductivity `k` across two half-cells of 0.5 m and
/// the link's 1 m^2: `1 / (0.5 / k + 0.5 / k)`, the fixture's authored
/// conductance when `k` is that conductance.
fn timeline(
    curves: &'static [ResolvedHeatCapacityCurve; 2],
    conductance_w_k: f64,
    dt_s: f64,
) -> (ThermalTimeline<'static>, TimelineState<'static>) {
    let start: DepositionSimulation<'static> = initial(curves, conductance_w_k, None);
    let laws: &'static [CardConduction; 2] =
        Box::leak(Box::new(core::array::from_fn(|i| CardConduction {
            curve: &curves[i],
            k_w_m_k: conductance_w_k,
        })));
    let observer = ConductionObserver::new(
        &start,
        &[
            (MaterialId(0), &laws[0] as &dyn ConductivityLaw),
            (MaterialId(1), &laws[1]),
        ],
        &[ConductiveLink {
            link_id: 7,
            half_distances_m: [0.5; 2],
            contact_resistance_m2_k_w: 0.0,
        }],
        &[],
        || Ok(()),
    )
    .unwrap();
    let conduction = ImplicitConductionGraph::new(
        &start,
        &observer,
        &[
            (MaterialId(0), &laws[0] as &dyn ImplicitConductionLaw),
            (MaterialId(1), &laws[1]),
        ],
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap();
    let controls = ImplicitConductionControls {
        maximum_dt_s: dt_s,
        cell_tolerance_j: 1e-10,
        total_tolerance_j: 1e-9,
        ..Default::default()
    };
    let timeline = ThermalTimeline::implicit(
        CELLS.to_vec(),
        conduction,
        controls,
        ImplicitScheme::BackwardEuler,
    )
    .unwrap()
    .returning_end_temperatures()
    .unwrap();
    let cells_start = timeline.start(&start).unwrap();
    (timeline, cells_start)
}

/// The same device as a graph: the module's netlist with each resistor
/// following its card, the cells of [`timeline`], one heat route per resistor
/// into its destination cell.
fn graph(
    cards: &'static [MaterialCard; 2],
    curves: &'static [ResolvedHeatCapacityCurve; 2],
    conductance_w_k: f64,
    dt_s: f64,
) -> (System, NodeId, NodeId) {
    let (timeline, cells_start) = timeline(curves, conductance_w_k, dt_s);
    // The module's netlist; its resistances are placeholders the laws replace.
    let netlist = CircuitGraph {
        node_count: 3,
        branches: vec![
            (1, 0, Branch::Capacitor { farads: 1.0 }),
            (1, 2, Branch::Resistor { ohms: 999.0 }),
            (2, 0, Branch::Resistor { ohms: 999.0 }),
        ],
        transformers: vec![],
    };
    let mut circuit = KirchhoffCircuit::new(netlist, vec![], 1e-8).unwrap();
    for (tag, (destination, (card, curve))) in
        (0_u64..).zip(DESTINATIONS.iter().zip(cards.iter().zip(curves.iter())))
    {
        circuit = circuit
            .following(
                destination.branch,
                Arc::new(CardResistance { card, curve, tag }),
            )
            .unwrap();
    }
    let charge = circuit.state(&[], &[10.0]).unwrap();

    let mut builder = SystemBuilder::new();
    let electrical = builder.add("circuit", circuit, charge).unwrap();
    let cells = builder.add("cells", timeline, cells_start).unwrap();
    for destination in DESTINATIONS {
        builder
            .heat_route(
                builder
                    .port(electrical, &format!("heat{}", destination.branch))
                    .unwrap(),
                &[(
                    builder
                        .port(cells, &format!("cell{}", destination.cell.0))
                        .unwrap(),
                    1.0,
                )],
                None,
                RouteTemperature::SinkEquivalent,
            )
            .unwrap();
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
    for destination in DESTINATIONS {
        let port = system
            .port(electrical, &format!("heat{}", destination.branch))
            .unwrap();
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
                charge_c: system.state::<Vec<f64>>(electrical).unwrap()[0],
                kelvin: CELLS.map(|id| timeline.simulation().temperature_k(id).unwrap()),
                joule_j: joule - joule_before,
            });
            joule_before = joule;
        }
    });
    rows
}

#[test]
fn the_declared_graph_reproduces_the_endpoint_resistance_feedback() {
    // The module's own feedback fixture: resistor 0 rises 0.5 % per kelvin,
    // the cells share a 1 W/K link, eight quarter-second intervals.
    let (a, ca) = material(1.0, 1.0, 0.005);
    let (b, cb) = material(3.0, 2.0, 0.0);
    let cards: &'static [MaterialCard; 2] = Box::leak(Box::new([a, b]));
    let curves: &'static [ResolvedHeatCapacityCurve; 2] = Box::leak(Box::new([ca, cb]));
    let (dt, steps, conductance) = (0.25, 8, 1.0);
    let module = module_run(cards, curves, conductance, dt, steps);
    let (mut system, electrical, cells) = graph(cards, curves, conductance, dt);
    let kernel = graph_run(&mut system, electrical, cells, dt, steps);

    // The feedback is real: cell 30 warms by kelvins, so resistor 0 moves by
    // several percent of itself over the run.
    assert!(
        module[steps - 1].kelvin[0] > 302.0,
        "{:?}",
        module[steps - 1]
    );
    for (k, (ours, theirs)) in kernel.iter().zip(&module).enumerate() {
        // Temperatures and charge agree bit for bit over the first four
        // intervals and the Joule heat within 9e-16 J; after that all agree
        // within 2e-9 K, 1.4e-12 C and 4e-12 J (export of 1e611fdf33): the
        // module's Newton and the kernel's interface iteration stop at
        // different points inside their tolerances. The bounds sit about
        // fifty times above that, far below the module's own 1e-6 J stop.
        for cell in 0..2 {
            let (x, y) = (ours.kelvin[cell], theirs.kelvin[cell]);
            assert!(
                (x - y).abs() <= 1e-7,
                "interval {k}, cell {cell}: {x} vs {y}"
            );
        }
        assert!(
            (ours.charge_c - theirs.charge_c).abs() <= 1e-10,
            "interval {k}: {} vs {}",
            ours.charge_c,
            theirs.charge_c
        );
        assert!(
            (ours.joule_j - theirs.joule_j).abs() <= 2e-10,
            "interval {k}: {} vs {}",
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
    // The field energy that left the capacitor is in the cells.
    let charge = kernel[steps - 1].charge_c;
    let field = (100.0 - charge * charge) / 2.0;
    let stored = ledger.nodes[cells.index()].stored_change_j.value();
    assert!((stored - field).abs() <= 1e-8, "{stored} vs {field}");
}
