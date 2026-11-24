//! H1's fixture (`material_contact_radiation_and_magnetic_entropy_close_same_voltage_driven_endpoint`)
//! as a declared graph on the composition kernel.
//!
//! `MaterialHeatNetwork::step_nonlinear_winding_candidate` couples by hand
//! the winding's field and conductors, its material graph resistor, link
//! 17's interface contact and the two cells' radiation. The graph says the
//! same with the kernel's pieces:
//! - fs-em's `PlanarWinding::magnetic_at_port_temperatures` with the
//!   module's own magnetic laws, the winding's conductors
//!   ([`CellConductors`]) at the winding cell's port temperature and the
//!   graph resistor ([`BodyResistor`]) at the sink's;
//! - kernel.rs's conduction graph of the two cells with link 17's contact
//!   (fs-conduction's `PreparedMaterialContact`, through `with_contacts`)
//!   and the enclosure (fs-conduction's `RadiationExchange`, through
//!   `with_coupled_heat_source`), as one backward-Euler timeline.
//!
//! The ports return `-T1 (s1 - s0)` beside the Joule heat where the module
//! deposits `-P - (T1 s1 - T0 s0)`; the Component books the difference,
//! `-(P + (T1 - T0) s0)`, as numerical loss, as in graph_caloric.rs.

use std::collections::BTreeMap;
use std::sync::Arc;

use fs_conduction::transient::deposition::contact::implicit::PreparedMaterialContact;
use fs_em::winding::component::{
    ConductorLaw, HeatReceiver, MagneticLaw, MagneticLaws, PlanarWinding, PlanarWindingState,
    ResistorLaw,
};
use fs_system::{NodeId, System};
use fs_thermal::component::TimelineState;
use fs_thermal::enthalpy::conduction::ConductiveLink;
use fs_thermal::enthalpy::conduction::implicit::contact::ImplicitContactLaw;

use super::*;
use crate::graph::{BodyResistor, Row, STEPS, TOLERANCE_J, graph_run, module_row};
use crate::graph_caloric::closed;
use crate::kernel::{CellConductors, conduction, graph, ledger_closes, timeline, with_cx};

/// Bounds of the kernel's states against heat.rs's closed forms on them,
/// over every variant: the current \[A\], the sink's heat balance, the total
/// energy balance with the term, and the booked loss against the owner's
/// dissipation and the term \[J\]. They measure 3.3e-11, 1.0e-9, 8.2e-10
/// and 7.3e-11 over the six variants; the bounds sit about fifty times
/// above.
const CLOSED: [f64; 4] = [1.7e-9, 5e-8, 4.1e-8, 3.7e-9];

/// At slope 0, the kernel against the module: temperature \[K\], current
/// \[A\], the end field potential \[Wb/m\] and the two cells' energy gap
/// \[J\]. They measure 8.2e-10, 4.4e-11, 1.3e-11 and 7.5e-10; the bounds
/// sit about fifty times above.
const CONTROL: [f64; 4] = [4.1e-8, 2.2e-9, 6.7e-10, 3.8e-8];

/// At slope 0.004, the same, with the gap's departure from the summed
/// terms as a fraction of them. They measure 1.04e-2, 2.3e-4, 5.3e-5 and
/// 1.26e-2, all from the cells' feedback through their resistivities,
/// contact, radiation and permeability, as in graph_caloric.rs; the bounds
/// sit about a third above.
const FEEDBACK: [f64; 4] = [1.4e-2, 3e-4, 7e-5, 1.7e-2];

/// The module's run of `parts` over graph.rs's intervals: H1's step from
/// `circuit`, its winding conductors and graph resistor held to `law` and
/// `resistor` at its own end temperatures, bit for bit. Also its field
/// potential \[Wb/m\] at the end.
fn heat_module(
    parts: &Parts,
    circuit: &NonlinearWindingCircuit<'_>,
    law: &CellConductors,
    resistor: &BodyResistor,
) -> (Vec<Row>, Vec<f64>) {
    let mut thermal = parts.thermal.clone();
    let mut circuit = circuit.clone();
    let mut rows = Vec::new();
    with_cx(|cx| {
        let prepared = parts.assembly.prepare(cx, &thermal).unwrap();
        for step in 1..=STEPS {
            let next = prepared
                .step_nonlinear_winding_candidate(
                    cx,
                    &thermal,
                    step as f64 * DT,
                    &circuit,
                    &[VOLTAGE],
                    &routes(),
                    &DESTINATIONS,
                    1e-12,
                    Some((parts.bindings.as_slice(), linear::limits())),
                    config(),
                    field::controls(),
                    &[0.0; 2],
                )
                .unwrap();
            let t = [CELL, SINK].map(|id| next.thermal().temperature_k(id).unwrap());
            let electrical = &next.report().electrical;
            for loss in &electrical.cells {
                let triangle = loss.source.triangle;
                assert_eq!(
                    law.conductor(triangle, t[0]).unwrap(),
                    loss.source.conductor,
                    "step {step}, triangle {triangle}"
                );
            }
            assert_eq!(
                resistor.conductor(t[1]).unwrap(),
                electrical.graph_materials[0].conductor,
                "step {step}"
            );
            rows.push(module_row(t.to_vec(), next.circuit(), electrical, [0.0; 2]));
            (thermal, circuit, _) = next.into_parts();
        }
    });
    (rows, circuit.field().potential_wb_m().to_vec())
}

/// The winding's conductors and its graph resistor for `parts`, from
/// `start`'s conductors and the module's resistor.
fn port_laws(
    parts: &Parts,
    start: &NonlinearWindingCircuit<'static>,
) -> (Arc<CellConductors>, Arc<BodyResistor>) {
    let law = CellConductors::new(
        &parts.cards[0],
        parts.context.clone(),
        start.field().conductors(),
        start.conductors(),
    );
    let resistor = BodyResistor::new(
        &parts.cards[1],
        parts.context.clone(),
        parts.key,
        &parts.assembly.electrical.bulk.resistors[0],
    );
    (Arc::new(law), Arc::new(resistor))
}

/// The graph of `parts` from `start`: the winding with `law` and
/// `resistor` at its ports' temperatures beside the two cells' conduction
/// graph with link 17's contact and the enclosure, from the cells' start
/// temperatures.
fn heat_graph(
    parts: &Parts,
    start: &NonlinearWindingCircuit<'static>,
    law: &Arc<CellConductors>,
    resistor: &Arc<BodyResistor>,
) -> (System, NodeId, NodeId) {
    let magnetic = MagneticLaws {
        triangles: parts
            .bindings
            .iter()
            .map(|binding| MagneticLaw {
                triangle: binding.triangle,
                law: binding.law.clone(),
            })
            .collect(),
        assembly: linear::limits(),
    };
    let winding = PlanarWinding::magnetic_at_port_temperatures(
        start,
        vec![law.clone(), law.clone()],
        BTreeMap::from([(2, resistor.clone() as Arc<dyn ResistorLaw>)]),
        magnetic,
        vec![
            HeatReceiver {
                triangles: (0..8).collect(),
                resistors: vec![],
            },
            HeatReceiver {
                triangles: vec![],
                resistors: vec![2],
            },
        ],
        field::controls(),
        TOLERANCE_J,
    )
    .unwrap();

    // Link 17 conducts through the cards' 1 W/(m K) over 0.5 m to
    // each cell's centre, in series with the interface's contact;
    // the enclosure exchanges the cells' gray radiation.
    let bounds = parts.curves.each_ref().map(|curve| {
        let knots = curve.curve().knots();
        [&knots[0], &knots[knots.len() - 1]].map(|knot| knot.specific_enthalpy_j_kg)
    });
    let contacts = leak([(
        17,
        leak(
            PreparedMaterialContact::new(
                &parts.thermal,
                parts.assembly.thermal_contacts[0],
                bounds,
                || Ok(()),
            )
            .unwrap(),
        ) as &dyn ImplicitContactLaw,
    )]);
    let exchange = leak(
        parts
            .assembly
            .radiation
            .unwrap()
            .exchange(&parts.thermal, || Ok(()))
            .unwrap(),
    );
    let cells_graph = conduction(
        &parts.thermal,
        &[
            (MATERIAL, &parts.curves[0]),
            (SINK_MATERIAL, &parts.curves[1]),
        ],
        &[ConductiveLink {
            link_id: 17,
            half_distances_m: [0.5; 2],
            contact_resistance_m2_k_w: 0.0,
        }],
    )
    .with_contacts(&parts.thermal, contacts, || Ok(()))
    .unwrap()
    .with_coupled_heat_source(&parts.thermal, exchange.cells(), exchange, 4, || Ok(()))
    .unwrap();
    let t0 = [CELL, SINK].map(|id| parts.thermal.temperature_k(id).unwrap());
    graph(
        winding,
        PlanarWindingState::at_temperatures(start.clone(), t0.to_vec()),
        timeline(&parts.thermal, cells_graph),
        &[CELL, SINK],
        Some(|_| VOLTAGE),
    )
}

/// Each interval of `kernel`, from the cells at `t0` and no current, against
/// heat.rs's closed forms on its own states, in [`CLOSED`]'s order: the
/// current, the sink's heat balance, the total balance less the term the
/// module deposits and the kernel books, and the booked loss against the
/// owner's dissipation and the term. Also the terms summed \[J\].
fn closed_forms(
    kernel: &[Row],
    t0: [f64; 2],
    (e, contact, slope): (f64, f64, f64),
) -> ([f64; 4], f64) {
    let (mut from, mut most, mut summed) = ((t0, 0.0), [0.0_f64; 4], 0.0);
    for row in kernel {
        let to = ([row.kelvin[0], row.kelvin[1]], row.current_a);
        let [numerical, term] = closed((from.1, from.0[0]), (to.1, to.0[0]), slope);
        let [sink, total] = balances(from.0, from.1, to.0, to.1, (e, contact, slope));
        let values = [
            (to.1 - current(from.0, from.1, to.0, slope)).abs(),
            sink.abs(),
            (total + term).abs(),
            (row.numerical_j - numerical - term).abs(),
        ];
        for (most, value) in most.iter_mut().zip(values) {
            *most = most.max(value);
        }
        summed += term;
        from = to;
    }
    (most, summed)
}

/// The largest differences of `kernel` and `module` over the intervals:
/// temperature \[K\], current \[A\] and the end field potentials `ours`
/// and `theirs` \[Wb/m\].
fn farthest(kernel: &[Row], module: &[Row], ours: &[f64], theirs: &[f64]) -> [f64; 3] {
    let mut far = [0.0_f64; 3];
    for (ours, theirs) in kernel.iter().zip(module) {
        for (x, y) in ours.kelvin.iter().zip(&theirs.kelvin) {
            far[0] = far[0].max((x - y).abs());
        }
        far[1] = far[1].max((ours.current_a - theirs.current_a).abs());
    }
    for (x, y) in ours.iter().zip(theirs) {
        far[2] = far[2].max((x - y).abs());
    }
    far
}

/// H1's three variants at slope 0, the control, and at H1's 0.004.
#[test]
fn the_declared_graph_and_the_heat_assembly_differ_by_the_booked_term() {
    for slope in [0.0, SLOPE] {
        let mut ends = Vec::new();
        for (e, contact) in [(0.2, 0.2), (0.8, 0.2), (0.2, 2.0)] {
            let parts = parts(e, contact, 500.0, slope);
            let start = circuit(parts.field, 0.0);
            let (law, resistor) = port_laws(&parts, &start);
            let (mut system, electrical, cells) = heat_graph(&parts, &start, &law, &resistor);
            let kernel = graph_run(&mut system, electrical, cells, &[CELL, SINK]);
            ledger_closes(&system, electrical, cells, 0.0);
            let node = &system.ledger().nodes[electrical.index()];
            assert_eq!(node.entropy.as_ref().unwrap().audited_steps, STEPS as u64);

            let t0 = [CELL, SINK].map(|id| parts.thermal.temperature_k(id).unwrap());
            let (most, summed) = closed_forms(&kernel, t0, (e, contact, slope));
            let end = &kernel[STEPS - 1];
            ends.push([end.kelvin[0], end.kelvin[1], end.current_a]);
            println!(
                "graph H1 (slope {slope}, e {e}, contact {contact}): closed forms {:.3e} A, {:.3e} J, \
                 {:.3e} J, {:.3e} J; terms {summed:.4e} J; end {:.4} K, {:.4} K, {:.4} A",
                most[0], most[1], most[2], most[3], end.kelvin[0], end.kelvin[1], end.current_a
            );
            for (value, bound) in most.into_iter().zip(CLOSED) {
                assert!(value <= bound, "{value:e} above {bound:e}");
            }

            let (module, potential) = heat_module(&parts, &start, &law, &resistor);
            let state = system
                .state::<PlanarWindingState<'static>>(electrical)
                .unwrap();
            let far = farthest(
                &kernel,
                &module,
                state.circuit().field().potential_wb_m(),
                &potential,
            );
            let gap: f64 = (0..2)
                .map(|cell| HEAT_CAPACITY * (module[STEPS - 1].kelvin[cell] - end.kelvin[cell]))
                .sum();
            println!(
                "graph H1 (slope {slope}, e {e}, contact {contact}): kernel against module: \
                 {:.3e} K, {:.3e} A, {:.3e} Wb/m; cells' gap {gap:.4e} J",
                far[0], far[1], far[2]
            );
            if slope == 0.0 {
                assert_eq!(summed.to_bits(), 0.0_f64.to_bits());
                let values = [far[0], far[1], far[2], gap.abs()];
                for (value, bound) in values.into_iter().zip(CONTROL) {
                    assert!(value <= bound, "slope 0: {values:?} against {CONTROL:?}");
                }
            } else {
                // The term is resolved, two orders above the bound of its
                // booking; the module's cells end cooler by it, less the
                // feedback of their resistivities, contact, radiation and
                // permeability.
                assert!(-summed > 100.0 * CLOSED[3], "{summed:e}");
                let fraction = (gap - summed).abs() / summed.abs();
                println!(
                    "graph H1 (e {e}, contact {contact}): the gap departs from the terms by {fraction:.3e} of them"
                );
                let values = [far[0], far[1], far[2], fraction];
                for (value, bound) in values.into_iter().zip(FEEDBACK) {
                    assert!(
                        value <= bound,
                        "slope {slope}: {values:?} against {FEEDBACK:?}"
                    );
                }
            }
        }
        // Emissivity and contact each move both cells and the current, as
        // H1 asserts of the module.
        for changed in &ends[1..] {
            for (x, y) in changed.iter().zip(&ends[0]) {
                assert!((x - y).abs() > 1e-5, "slope {slope}: {ends:?}");
            }
        }
    }
}

/// H3's fixture (`missing_endpoint_emissivity_refuses_without_advancing_any_native_history`):
/// H1's at 0.004 per kelvin, with the emissivity claim ending at 300.01 K,
/// which the sink would pass within one interval. The module refuses the
/// whole interval. The kernel shortens it instead: each step ends with both
/// cells inside the claim, and the circuit, the field and the cells all at
/// the step's end time.
#[test]
fn the_declared_graph_shortens_the_interval_at_the_end_of_the_emissivity_claim() {
    let parts = parts(0.2, 0.2, 300.01, SLOPE);
    let start = circuit(parts.field, 0.0);
    let (law, resistor) = port_laws(&parts, &start);
    let (mut system, electrical, cells) = heat_graph(&parts, &start, &law, &resistor);
    with_cx(|cx| {
        for _ in 0..4 {
            let report = system.step(DT, cx).unwrap();
            let timeline = system.state::<TimelineState<'static>>(cells).unwrap();
            let kelvin = [CELL, SINK].map(|id| timeline.simulation().temperature_k(id).unwrap());
            println!("graph H3: {report:?}; cells {kelvin:?} K");
            assert!(report.shortenings > 0 && report.dt_s < DT, "{report:?}");
            assert!(kelvin.iter().all(|t| *t <= 300.01), "{kelvin:?}");
            let circuit = system
                .state::<PlanarWindingState<'static>>(electrical)
                .unwrap()
                .circuit();
            let times = [
                circuit.time_s(),
                circuit.field().time_s(),
                timeline.simulation().time_s(),
            ];
            for time in times {
                assert_eq!(time.to_bits(), report.time_s.to_bits(), "{times:?}");
            }
        }
    });
}
