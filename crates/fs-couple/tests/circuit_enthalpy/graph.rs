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
//! graph), so they agree only to within the solves' tolerances.

use std::sync::Arc;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_couple::cavity_phs::piston::thermal::power::joule::MaterialOhmicBody;
use fs_couple::cavity_phs::piston::thermal::source::ThermalMaterialPoint;
use fs_couple::circuit_enthalpy::step_candidate;
use fs_couple::material_circuit::MaterialCircuit;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::MaterialCard;
use fs_material::phase::ResolvedHeatCapacityCurve;
use fs_material::state_point::MaterialPropertySelection;
use fs_phs::circuit::{Branch, CircuitGraph, MutualInductance};
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
/// at the temperature the resistor's heat port receives, its geometry
/// checked against the cell's volume.
#[derive(Debug)]
struct CardResistance {
    card: &'static MaterialCard,
    curve: &'static ResolvedHeatCapacityCurve,
    /// The destination cell's volume \[m^3\].
    volume_m3: f64,
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
        body.resolve_conductor(phase, self.volume_m3)
            .map(|conductor| conductor.resistance_ohm())
            .map_err(|error| format!("{error:?}"))
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        words.extend([self.tag, self.volume_m3.to_bits()]);
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

/// One accepted interval: the circuit's differential state (inductor flux
/// linkages \[Wb\], then capacitor charges \[C\]), the two cells' temperatures
/// \[K\] and the Joule heat \[J\].
#[derive(Debug, Clone)]
struct Row {
    state: Vec<f64>,
    kelvin: [f64; 2],
    joule_j: f64,
}

const CELLS: [CellId; 2] = [CellId(30), CellId(10)];

/// A device of the module's tests: its netlist, whose resistances are
/// placeholders the cards replace, its mutual inductances and its start.
struct Device {
    netlist: CircuitGraph,
    mutual: Vec<MutualInductance>,
    fluxes_wb: Vec<f64>,
    charges_c: Vec<f64>,
}

/// The module's run: `steps` intervals of `dt_s` with `step_candidate`, the
/// device built as the module's tests build it.
fn module_run(
    device: &Device,
    cards: &[MaterialCard; 2],
    curves: &[ResolvedHeatCapacityCurve; 2],
    conductance_w_k: f64,
    dt_s: f64,
    steps: usize,
) -> Vec<Row> {
    let mut thermal = initial(curves, conductance_w_k, None);
    let mut electrical = MaterialCircuit::new_with_mutual_inductance(
        device.netlist.clone(),
        coefficients(&thermal, cards, curves).unwrap().resistors,
        device.mutual.clone(),
        device.fluxes_wb.clone(),
        device.charges_c.clone(),
        &[],
    )
    .unwrap();
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
        let circuit = next.circuit();
        rows.push(Row {
            state: circuit
                .fluxes_wb()
                .iter()
                .chain(circuit.charges_c())
                .copied()
                .collect(),
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

/// The same device as a graph: its netlist with each resistor following its
/// card, the cells of [`timeline`], one heat route per resistor into its
/// destination cell.
fn graph(
    device: &Device,
    cards: &'static [MaterialCard; 2],
    curves: &'static [ResolvedHeatCapacityCurve; 2],
    conductance_w_k: f64,
    dt_s: f64,
) -> (System, NodeId, NodeId) {
    let (timeline, cells_start) = timeline(curves, conductance_w_k, dt_s);
    let start = initial(curves, conductance_w_k, None);
    let mut circuit =
        KirchhoffCircuit::new(device.netlist.clone(), device.mutual.clone(), 1e-8).unwrap();
    for (tag, (destination, (card, curve))) in
        (0_u64..).zip(DESTINATIONS.iter().zip(cards.iter().zip(curves.iter())))
    {
        circuit = circuit
            .following(
                destination.branch,
                Arc::new(CardResistance {
                    card,
                    curve,
                    volume_m3: start.cells()[&destination.cell].geometry.volume_m3,
                    tag,
                }),
            )
            .unwrap();
    }
    let state = circuit.state(&device.fluxes_wb, &device.charges_c).unwrap();

    let mut builder = SystemBuilder::new();
    let electrical = builder.add("circuit", circuit, state).unwrap();
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
                state: system.state::<Vec<f64>>(electrical).unwrap().clone(),
                kelvin: CELLS.map(|id| timeline.simulation().temperature_k(id).unwrap()),
                joule_j: joule - joule_before,
            });
            joule_before = joule;
        }
    });
    rows
}

/// Bounds of one comparison: temperature \[K\], state coordinate \[Wb or C\]
/// and Joule heat per interval \[J\].
struct Bounds {
    kelvin: f64,
    state: f64,
    joule_j: f64,
}

/// Runs the device both ways and holds every interval to `bounds`. Returns
/// the module's rows. The kernel's ledger closes, and the energy that left
/// the circuit is in the cells.
fn compare(
    device: &Device,
    cards: &'static [MaterialCard; 2],
    curves: &'static [ResolvedHeatCapacityCurve; 2],
    (dt, steps): (f64, usize),
    bounds: &Bounds,
) -> Vec<Row> {
    let conductance = 1.0;
    let module = module_run(device, cards, curves, conductance, dt, steps);
    let (mut system, electrical, cells) = graph(device, cards, curves, conductance, dt);
    let kernel = graph_run(&mut system, electrical, cells, dt, steps);
    assert_eq!((kernel.len(), module.len()), (steps, steps));
    for (k, (ours, theirs)) in kernel.iter().zip(&module).enumerate() {
        for cell in 0..2 {
            let (x, y) = (ours.kelvin[cell], theirs.kelvin[cell]);
            assert!(
                (x - y).abs() <= bounds.kelvin,
                "interval {k}, cell {cell}: {x} vs {y}"
            );
        }
        assert_eq!(ours.state.len(), theirs.state.len());
        for (i, (x, y)) in ours.state.iter().zip(&theirs.state).enumerate() {
            assert!(
                (x - y).abs() <= bounds.state,
                "interval {k}, coordinate {i}: {x} vs {y}"
            );
        }
        assert!(
            (ours.joule_j - theirs.joule_j).abs() <= bounds.joule_j,
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
    let released = -ledger.nodes[electrical.index()].stored_change_j.value();
    let stored = ledger.nodes[cells.index()].stored_change_j.value();
    assert!(
        released > 1.0 && (stored - released).abs() <= 1e-8,
        "{stored} vs {released}"
    );
    module
}

fn cards(
    alphas: [f64; 2],
) -> (
    &'static [MaterialCard; 2],
    &'static [ResolvedHeatCapacityCurve; 2],
) {
    let (a, ca) = material(1.0, 1.0, alphas[0]);
    let (b, cb) = material(3.0, 2.0, alphas[1]);
    (Box::leak(Box::new([a, b])), Box::leak(Box::new([ca, cb])))
}

#[test]
fn the_declared_graph_reproduces_the_endpoint_resistance_feedback() {
    // The module's own feedback fixture: a 1 F capacitor at 10 C through the
    // two card resistors, resistor 0 rising 0.5 % per kelvin, the cells
    // sharing a 1 W/K link, eight quarter-second intervals.
    let (cards, curves) = cards([0.005, 0.0]);
    let device = Device {
        netlist: CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, Branch::Capacitor { farads: 1.0 }),
                (1, 2, Branch::Resistor { ohms: 999.0 }),
                (2, 0, Branch::Resistor { ohms: 999.0 }),
            ],
            transformers: vec![],
        },
        mutual: vec![],
        fluxes_wb: vec![],
        charges_c: vec![10.0],
    };
    // Temperatures and charge agree bit for bit over the first four
    // intervals and the Joule heat within 9e-16 J; after that all agree
    // within 2e-9 K, 1.4e-12 C and 4e-12 J (export of 1e611fdf33): the
    // module's Newton and the kernel's interface iteration stop at
    // different points inside their tolerances. The bounds sit about
    // fifty times above that, far below the module's own 1e-6 J stop.
    let bounds = Bounds {
        kelvin: 1e-7,
        state: 1e-10,
        joule_j: 2e-10,
    };
    let module = compare(&device, cards, curves, (0.25, 8), &bounds);
    // The feedback is real: cell 30 warms by kelvins, so resistor 0 moves by
    // several percent of itself over the run.
    assert!(module[7].kelvin[0] > 302.0, "{:?}", module[7]);
}

#[test]
fn the_declared_graph_reproduces_mutual_coils_with_both_resistors_following_their_cards() {
    // The module's mutual-coil fixture: 2 H and 3 H coupled by 1 H, each
    // shunted by a card resistor (0.5 % and 0.2 % per kelvin), linkages
    // 3.5 and 0.5 Wb, six quarter-second intervals.
    let (cards, curves) = cards([0.005, 0.002]);
    let device = Device {
        netlist: CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, Branch::Inductor { henries: 2.0 }),
                (1, 0, Branch::Resistor { ohms: 999.0 }),
                (2, 0, Branch::Resistor { ohms: 999.0 }),
                (2, 0, Branch::Inductor { henries: 3.0 }),
            ],
            transformers: vec![],
        },
        mutual: vec![MutualInductance {
            branches: (0, 3),
            henries: 1.0,
        }],
        fluxes_wb: vec![3.5, 0.5],
        charges_c: vec![],
    };
    // On the export of 1e611fdf33 the two agree within 2e-9 K, 1.1e-11 Wb
    // and 7.1e-12 J at every interval; the bounds sit about fifty times
    // above that.
    let bounds = Bounds {
        kelvin: 1e-7,
        state: 5e-10,
        joule_j: 4e-10,
    };
    let module = compare(&device, cards, curves, (0.25, 6), &bounds);
    assert!(module[5].kelvin[0] > 300.5, "{:?}", module[5]);
}
