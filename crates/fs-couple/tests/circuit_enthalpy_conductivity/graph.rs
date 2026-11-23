//! The conductivity coupling as a declared graph on the composition kernel.
//!
//! `circuit_enthalpy::conductivity::step_candidate` couples the circuit to the
//! enthalpy timeline by hand: backward Euler for the cells, the circuit's
//! midpoint step with each resistor at the material resistance of its cell's
//! END temperature and, through the material circuit's constitutive
//! iteration, of the field or current density the step solves. The graph
//! says the same thing with no code for the pair:
//!
//! - `fs_system::circuit::KirchhoffCircuit` with the module's conductivity
//!   controls as its `OperatingPointIteration`, each resistor following its
//!   card through [`CardConductor`] at the temperature its heat port receives
//!   and, for the axes the card binds, at its operating point;
//! - the backward-Euler timeline of `kernel::timeline`, returning end-of-step
//!   temperatures;
//! - one heat route per resistor into the module's destination cell.
//!
//! The two sides solve the same equations with different solvers (the
//! module's Newton--Krylov on cell energies around the material circuit's
//! iteration; the kernel's interface iteration around the circuit's), so
//! they agree only to within the solves' tolerances.

use std::sync::Arc;

use fs_couple::circuit_enthalpy::CircuitEnthalpyError;
use fs_couple::circuit_enthalpy::conductivity::step_candidate;
use fs_couple::material_circuit::conductivity::{
    ConductionAxis, ConductionCoordinate, ConductivityIteration, ConductorResponse,
};
use fs_couple::material_circuit::{MaterialCircuit, MaterialCircuitError, MaterialResistor};
use fs_matdb::{ClaimSet, MaterialCard, PropertyValue, QueryPoint};
use fs_material::conductor::{
    ELECTRICAL_RESISTIVITY_PROPERTY, ResolvedConductor, ResolvedConductorState,
    resolve_conductor_state,
};
use fs_material::phase::ResolvedHeatCapacityCurve;
use fs_material::state_point::MaterialPropertySelection;
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_qty::QuantitySpec;
use fs_system::circuit::{
    CircuitComponentError, KirchhoffCircuit, OperatingPointIteration, ResistanceLaw,
};
use fs_system::{NodeId, RouteTemperature, System, SystemBuilder, SystemError, heat_sink};
use fs_thermal::component::TimelineState;

use crate::common::{self, DESTINATIONS, config, initial, input, point};
use crate::kernel::{CELLS, enthalpy_at, timeline, with_cx};
use crate::{AXES, CURRENT_AXIS, CURRENT_DENSITY, coefficients, field_caloric, field_material};

/// A resistor of the module's fixtures as its conductivity step resolves it:
/// the seed conductor of the original circuit at the phase state of the
/// temperature its port receives, as the fixtures' `coefficients` resolve it,
/// then at the field `V / L` or current density `I / A` of its operating
/// point on each axis the card binds, as the material circuit refreshes it.
#[derive(Debug)]
struct CardConductor {
    card: &'static MaterialCard,
    curve: &'static ResolvedHeatCapacityCurve,
    seed: &'static ResolvedConductor,
    axes: &'static [ConductionAxis<'static>],
    /// Names the card in the circuit's identity: the fixture's index.
    tag: u64,
}

impl CardConductor {
    fn resolve(&self, kelvin: f64, volts: f64, amperes: f64) -> Result<f64, String> {
        let h = enthalpy_at(self.curve, kelvin)?;
        let phase = self
            .curve
            .curve()
            .state_at_specific_enthalpy(h)
            .map_err(|error| format!("{error:?}"))?;
        let context = self.seed.material().query();
        let physical = ThermalMaterialPoint {
            material: self.card,
            context: &context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        }
        .at(phase)
        .map_err(|error| format!("{error:?}"))?;
        let bundle = self
            .seed
            .material()
            .at_query_point(self.card, &physical)
            .map_err(|error| format!("{error:?}"))?;
        let conductor =
            ResolvedConductorState::try_from_material_with_basis(bundle, self.seed.source_basis())
                .map_err(|error| format!("{error:?}"))?
                .into_uniform_conductor(self.seed.length_m(), self.seed.area_m2())
                .map_err(|error| format!("{error:?}"))?;
        if self.axes.is_empty() {
            return Ok(conductor.resistance_ohm());
        }
        let field = volts / self.seed.length_m();
        let density = amperes / self.seed.area_m2();
        let mut query = conductor.material().query();
        for axis in self.axes {
            let value = match axis.coordinate {
                ConductionCoordinate::ElectricField => field,
                ConductionCoordinate::ElectricFieldMagnitude => field.abs(),
                ConductionCoordinate::CurrentDensity => density,
                ConductionCoordinate::CurrentDensityMagnitude => density.abs(),
            };
            query = query
                .with(axis.name, value)
                .map_err(|error| format!("{error:?}"))?;
        }
        let material = conductor
            .material()
            .at_query_point(self.card, &query)
            .map_err(|error| format!("{error:?}"))?;
        conductor
            .with_material(material)
            .map(|refreshed| refreshed.resistance_ohm())
            .map_err(|error| format!("{error:?}"))
    }
}

impl ResistanceLaw for CardConductor {
    fn ohms(&self, kelvin: f64) -> Result<f64, String> {
        self.resolve(kelvin, 0.0, 0.0)
    }

    fn follows_operating_point(&self) -> bool {
        !self.axes.is_empty()
    }

    fn ohms_at(&self, kelvin: f64, volts: f64, amperes: f64) -> Result<f64, String> {
        self.resolve(kelvin, volts, amperes)
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        words.push(self.tag);
    }
}

/// A fixture of the module's conductivity tests: the cards and caloric curves
/// of resistors 1 and 2 (cells 30 and 10), the axes each binds, the seed
/// circuit whose resistors carry the cards' start queries (a 1 F capacitor at
/// 10 C), and the module's conductivity controls.
struct Fixture {
    cards: &'static [MaterialCard; 2],
    curves: &'static [ResolvedHeatCapacityCurve; 2],
    axes: [&'static [ConductionAxis<'static>]; 2],
    seed: &'static MaterialCircuit,
    iteration: ConductivityIteration,
}

impl Fixture {
    /// The fixture with resistor `k` resolved from `cards[k]` at `queries[k]`.
    fn new(
        cards: [MaterialCard; 2],
        curves: [ResolvedHeatCapacityCurve; 2],
        queries: &[QueryPoint; 2],
        axes: [&'static [ConductionAxis<'static>]; 2],
        iteration: ConductivityIteration,
    ) -> Self {
        let cards: &'static [MaterialCard; 2] = Box::leak(Box::new(cards));
        let resistors = DESTINATIONS
            .iter()
            .zip(cards.iter().zip(queries))
            .map(|(destination, (card, query))| MaterialResistor {
                branch: destination.branch,
                conductor: resolve_conductor_state(
                    card,
                    query,
                    MaterialPropertySelection::SingleClaimOnly,
                )
                .unwrap()
                .into_uniform_conductor(1.0, 1.0)
                .unwrap(),
            })
            .collect();
        let seed = MaterialCircuit::new(netlist(), resistors, vec![], vec![10.0], &[]).unwrap();
        Self {
            cards,
            curves: Box::leak(Box::new(curves)),
            axes,
            seed: Box::leak(Box::new(seed)),
            iteration,
        }
    }

    /// The module's controls as the kernel circuit's iteration, at
    /// `relaxation`, with the circuit's power tolerance of `config`.
    fn kernel_iteration(&self, relaxation: f64) -> OperatingPointIteration {
        let solve = self.iteration.controls();
        OperatingPointIteration {
            maximum_trials: solve.maximum_iterations,
            absolute_resistance_tolerance_ohm: solve.absolute_resistance_tolerance_ohm,
            relative_resistance_tolerance: solve.relative_resistance_tolerance,
            power_tolerance_w: config().circuit.power_tolerance_w,
            relaxation,
        }
    }
}

/// The fixtures' netlist: the capacitor across node 1, resistor 1 from node
/// 1 to node 2 and resistor 2 to ground, at placeholder resistances.
fn netlist() -> CircuitGraph {
    CircuitGraph {
        node_count: 3,
        branches: vec![
            (1, 0, Branch::Capacitor { farads: 1.0 }),
            (1, 2, Branch::Resistor { ohms: 999.0 }),
            (2, 0, Branch::Resistor { ohms: 999.0 }),
        ],
        transformers: vec![],
    }
}

/// One accepted interval: the capacitor's charge \[C\], the two cells'
/// temperatures \[K\] and the Joule heat \[J\].
#[derive(Debug, Clone)]
struct Row {
    charge_c: f64,
    kelvin: [f64; 2],
    joule_j: f64,
}

/// The module's run: `steps` intervals of `dt_s`, every trial's resistors
/// resolved from the seed circuit by the fixtures' `coefficients`.
fn module_run(
    fixture: &Fixture,
    dt_s: f64,
    steps: usize,
) -> Result<Vec<Row>, CircuitEnthalpyError<String>> {
    let mut thermal = initial(fixture.curves, 1.0, None);
    let mut electrical = fixture.seed.clone();
    let responses: Vec<ConductorResponse<'_>> = (0..2)
        .map(|k| ConductorResponse {
            branch: DESTINATIONS[k].branch,
            card: &fixture.cards[k],
            axes: fixture.axes[k],
        })
        .collect();
    let mut rows = Vec::new();
    for k in 1..=steps {
        let next = step_candidate(
            &thermal,
            k as f64 * dt_s,
            input(&electrical),
            &responses,
            config(),
            fixture.iteration,
            &[0.0; 2],
            |s| coefficients(s, fixture.seed, fixture.cards, fixture.curves),
            || false,
        )?;
        rows.push(Row {
            charge_c: next.circuit().charges_c()[0],
            kelvin: CELLS.map(|id| next.thermal().temperature_k(id).unwrap()),
            joule_j: next.report().deposited_joule_j,
        });
        let (t, c, _) = next.into_parts();
        thermal = t;
        electrical = c;
    }
    Ok(rows)
}

/// The fixture as a graph: the netlist with each resistor following its card
/// by `iteration`, the cells of [`timeline`] joined by 1 W/K, one heat route
/// per resistor into its destination cell, and the cells' environment port
/// routed to a sink at 300 K.
fn graph(
    fixture: &Fixture,
    iteration: OperatingPointIteration,
    dt_s: f64,
) -> (System, NodeId, NodeId) {
    let (timeline, cells_start) = timeline(fixture.curves, 1.0, [0.0; 2], dt_s);
    let mut circuit = KirchhoffCircuit::new(netlist(), vec![], 1e-8)
        .unwrap()
        .iterating(iteration)
        .unwrap();
    for (tag, (k, destination)) in (0_u64..).zip(DESTINATIONS.iter().enumerate()) {
        let seed = fixture
            .seed
            .resistors()
            .iter()
            .find(|resistor| resistor.branch == destination.branch)
            .unwrap();
        circuit = circuit
            .following(
                destination.branch,
                Arc::new(CardConductor {
                    card: &fixture.cards[k],
                    curve: &fixture.curves[k],
                    seed: &seed.conductor,
                    axes: fixture.axes[k],
                    tag,
                }),
            )
            .unwrap();
    }
    let state = circuit.state(&[], &[10.0]).unwrap();

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
) -> Result<Vec<Row>, SystemError> {
    let mut rows = Vec::new();
    let mut joule_before = 0.0;
    with_cx(|cx| {
        for _ in 0..steps {
            let report = system.step(dt_s, cx)?;
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
        Ok(rows)
    })
}

/// Bounds of one comparison: temperature \[K\], charge \[C\] and Joule heat
/// per interval \[J\].
struct Bounds {
    kelvin: f64,
    charge_c: f64,
    joule_j: f64,
}

/// Runs the fixture both ways over `steps` intervals of `dt` and holds every
/// interval to `bounds`. Returns the module's rows. The kernel's ledger
/// closes and the energy that left the capacitor is in the cells.
fn compare(fixture: &Fixture, (dt, steps): (f64, usize), bounds: &Bounds) -> Vec<Row> {
    let module = module_run(fixture, dt, steps).unwrap();
    let iteration = fixture.kernel_iteration(fixture.iteration.relaxation());
    let (mut system, electrical, cells) = graph(fixture, iteration, dt);
    let kernel = graph_run(&mut system, electrical, cells, dt, steps).unwrap();
    assert_eq!((kernel.len(), module.len()), (steps, steps));
    for (k, (ours, theirs)) in kernel.iter().zip(&module).enumerate() {
        for cell in 0..2 {
            let (x, y) = (ours.kelvin[cell], theirs.kelvin[cell]);
            assert!(
                (x - y).abs() <= bounds.kelvin,
                "interval {k}, cell {cell}: {x} vs {y}"
            );
        }
        assert!(
            (ours.charge_c - theirs.charge_c).abs() <= bounds.charge_c,
            "interval {k}: {} vs {} C",
            ours.charge_c,
            theirs.charge_c
        );
        assert!(
            (ours.joule_j - theirs.joule_j).abs() <= bounds.joule_j,
            "interval {k}: {} vs {} J",
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

#[test]
fn the_declared_graph_reproduces_field_and_temperature_feedback() {
    // The module's g1 fixture: resistor 1's resistivity 1 + 0.05 E with E
    // its field magnitude [V/m], resistor 2's 3 (1 + 0.005 (T - 300 K)), the
    // 1 F capacitor at 10 C, a 1 W/K link, unrelaxed; eight 0.2 s intervals.
    let (a, ca) = field_material(0.05, 2000.0);
    let (b, cb) = common::material(3.0, 2.0, 0.005);
    let fixture = Fixture::new(
        [a, b],
        [ca, cb],
        &[crate::query(300.0), point(300.0)],
        [&AXES, &[]],
        crate::conductivity_controls().into(),
    );
    // On the export of 1e611fdf33 with series 1-8 the two agree within
    // 3.5e-9 K, 9.2e-12 C and 9.9e-12 J at every interval: the module's
    // Newton and the kernel's interface iteration stop at different points
    // inside their tolerances. The bounds sit about fifty times above that.
    let bounds = Bounds {
        kelvin: 2e-7,
        charge_c: 5e-10,
        joule_j: 5e-10,
    };
    let module = compare(&fixture, (0.2, 8), &bounds);
    // The feedback is real: cell 10 warms by kelvins, so resistor 2 moves by
    // several percent of itself over the run.
    assert!(module[7].kelvin[1] > 304.0, "{:?}", module[7]);
}

/// The module's steep current-density card: resistivity 0.1 + 1000 J with J
/// in [0, 16] A/m^2, at the base card's constant heat capacity (the module's
/// relaxed fixture also varies the heat capacity, which the kernel's
/// timeline of one caloric segment does not hold).
fn current_material() -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let (base, _) = common::material(1.0, 1.0, 0.0);
    let mut claims = ClaimSet::new();
    for (_, original) in base.claims().claims_ordered() {
        let mut claim = original.clone();
        if claim.key.name() == ELECTRICAL_RESISTIVITY_PROPERTY {
            claim.value = PropertyValue::Curve {
                abscissa: "J".into(),
                abscissa_dims: CURRENT_DENSITY,
                knots: vec![(0.0, 0.1), (16.0, 16000.1)],
                dims: claim.key.quantity().dims(),
            };
            claim.validity = claim.validity.with_quantity(
                "J",
                QuantitySpec::dimensional(CURRENT_DENSITY),
                0.0,
                16.0,
            );
        }
        claims.insert_claim(claim).unwrap();
    }
    let card = MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap();
    let curve = field_caloric(&card);
    (card, curve)
}

#[test]
fn the_declared_graph_reproduces_a_steep_current_law_only_relaxed() {
    // The module's relaxed fixture with a constant heat capacity: resistor
    // 1's resistivity 0.1 + 1000 J, resistor 2's 3 (1 + 0.005 (T - 300 K)),
    // the iteration halfway to each solved point; six 0.2 s intervals.
    let (a, ca) = current_material();
    let (b, cb) = common::material(3.0, 2.0, 0.005);
    let start = point(300.0)
        .with_quantity("J", QuantitySpec::dimensional(CURRENT_DENSITY), 0.0)
        .unwrap();
    let relaxed = crate::conductivity_controls().with_relaxation(0.5).unwrap();
    let fixture = Fixture::new(
        [a, b],
        [ca, cb],
        &[start, point(300.0)],
        [&[CURRENT_AXIS], &[]],
        relaxed,
    );
    // On the same export the charge agrees to the bit at every interval, the
    // temperatures within 1.8e-12 K and the Joule heat within 1.7e-16 J. The
    // temperature and heat bounds sit about fifty times above that; the
    // charge bound is the field fixture's.
    let bounds = Bounds {
        kelvin: 1e-10,
        charge_c: 5e-10,
        joule_j: 1e-14,
    };
    compare(&fixture, (0.2, 6), &bounds);

    // Unrelaxed, neither side meets the resistance tolerance in 64 trials of
    // the first interval: the module refuses its step, the kernel's circuit
    // its proposal.
    let unrelaxed = Fixture {
        iteration: crate::conductivity_controls().into(),
        ..fixture
    };
    assert!(matches!(
        module_run(&unrelaxed, 0.2, 1),
        Err(CircuitEnthalpyError::Circuit(
            MaterialCircuitError::ConductivityConvergence { .. }
        ))
    ));
    let (mut system, electrical, cells) = graph(&unrelaxed, unrelaxed.kernel_iteration(1.0), 0.2);
    let refused = graph_run(&mut system, electrical, cells, 0.2, 1).unwrap_err();
    let SystemError::Component { node, error, .. } = &refused else {
        panic!("{refused:?}");
    };
    assert_eq!(node, "circuit");
    assert!(
        matches!(
            error.downcast_ref::<CircuitComponentError>(),
            Some(CircuitComponentError::OperatingPoint { trials: 64, .. })
        ),
        "{refused:?}"
    );
}
