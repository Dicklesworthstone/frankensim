//! Graphs 5 and 6: tests 5 and 6 of `hysteretic.rs`, the play-memory owner
//! (`SpatialHystereticWindingCircuit`) as fs-em's `SpatialWinding`, built
//! with `at_port_temperatures`. That owner refreshes no field material with
//! temperature, so it has no magnetic lane: its conductors and its graph
//! resistor follow the ports, its fits do not. The cube, card and one-mode
//! reference are `hysteretic.rs`'s.
//! - One heat port receives the heat of all 48 tetrahedra, eddy, strand,
//!   pinning and viscous alike, and, in graph 5, the graph resistor's; in
//!   graph 6 a second port receives that resistor's heat in its own body.
//! - The conductors follow the card at the port's temperature
//!   ([`CardConductors`]), and graph 6's resistor its body's card
//!   ([`BodyResistor`]).
//! - The source port holds test 5's 180 V, then 0 V; in graph 6, test 6's
//!   120 V, then 0 V.

use fs_couple::circuit_enthalpy::spatial_induction::SpatialMagneticThermalNetwork;
use fs_couple::circuit_enthalpy::spatial_induction::hysteretic::winding::SpatialWindingThermalNetwork;
use fs_em::spatial::vector_potential::induction::hysteretic::winding::{
    SpatialHystereticWindingCircuit, SpatialHystereticWindingStep,
};
use fs_em::spatial::vector_potential::induction::hysteretic::{
    SpatialHystereticInductionControls, SpatialHystereticProblem,
};
use fs_em::winding::nonlinear::resistors::GraphResistorMaterial;

use super::*;
use crate::hysteretic::data::{
    CAPACITY, DT, GRAPH_R, Mode, assert_field, circuit, controls, core_heat, cube, material,
    pinning_heat, voltage_endpoint,
};
use crate::hysteretic::{caloric, cells as module_cells, regions, thermal};

type PlayWinding = SpatialWinding<
    'static,
    ResolvedConductorState,
    SpatialHystereticWindingCircuit<'static>,
    SpatialHystereticInductionControls,
>;
type PlayState =
    SpatialWindingState<'static, ResolvedConductorState, SpatialHystereticWindingCircuit<'static>>;

/// Each interval's source voltage \[V\]: test 5's, then a release.
const DRIVE_5: [f64; STEPS] = [180.0, 0.0];
/// Test 6's, held. After it, the owner's circuit Newton refuses a release
/// to 0 V, in the module and the Component alike.
const DRIVE_6: [f64; STEPS] = [120.0, 120.0];

/// Test 5's and 6's card, its caloric curve and the cube, for the life of
/// the test.
fn fixture() -> (
    &'static Arc<MaterialCard>,
    &'static ResolvedHeatCapacityCurve,
    &'static SpatialHystereticProblem,
) {
    let card: &'static Arc<MaterialCard> =
        Box::leak(Box::new(material::card_with_heat_capacity(500.0, CAPACITY)));
    let curve: &'static ResolvedHeatCapacityCurve = Box::leak(Box::new(caloric(card)));
    (card, curve, Box::leak(Box::new(cube(card))))
}

/// The card's conductor law over `circuit`'s conductors at the module
/// cells' query point.
fn card_law(
    card: &MaterialCard,
    circuit: &SpatialHystereticWindingCircuit<'_>,
    trip: Option<(usize, &'static CancelGate)>,
) -> Arc<CardConductors> {
    let (core, strands) = (circuit.field().conductors(), circuit.conductors());
    let point = material::point(300.0);
    Arc::new(CardConductors::on(card, point, core, strands, trip))
}

/// Each cell's branch polarizations \[T\]: the play memory.
fn memory(circuit: &SpatialHystereticWindingCircuit<'_>) -> Vec<f64> {
    circuit
        .field()
        .cell_states()
        .iter()
        .flat_map(|cell| cell.polarizations_t().iter().flatten().copied())
        .collect()
}

/// One accepted interval, as [`Row`] holds it, with the play memory.
fn play_row(
    kelvin: Vec<f64>,
    circuit: &SpatialHystereticWindingCircuit<'_>,
    step: &SpatialHystereticWindingStep,
) -> (Row, Vec<f64>) {
    let strands = circuit.conductors().iter().map(|s| &s.conductor);
    let row = Row {
        kelvin,
        currents_a: circuit.currents_a().to_vec(),
        state: circuit.state().to_vec(),
        potential_wb: circuit.field().potential_wb().to_vec(),
        resistivity: strands
            .chain(circuit.field().conductors().iter().flatten())
            .map(ResolvedConductorState::resistivity_ohm_m)
            .collect(),
        resistance_ohm: step
            .graph_materials
            .iter()
            .map(|m| m.conductor.resistance_ohm())
            .collect(),
        flux_density_t: step.field.endpoint_flux_density_t.clone(),
        energies_j: [
            step.cells.iter().map(|cell| cell.heat_j).sum(),
            step.field.heat_j,
            step.graph_resistors.iter().map(|r| r.heat_j).sum(),
            step.heat_j,
            step.field.numerical_dissipation_j,
            step.numerical_dissipation_j,
            step.stored_change_j,
            step.supplied_work_j,
        ],
    };
    (row, memory(circuit))
}

/// The fixture's winding as the play-memory Component beside `bodies`'
/// timeline, its source port held at `drive` of the interval's start. The
/// owner closes its account only to its controls' 4e-7 J, so the Component
/// takes that allowance.
fn kernel(
    bodies: Bodies,
    circuit: &SpatialHystereticWindingCircuit<'static>,
    law: Arc<CardConductors>,
    drive: fn(f64) -> f64,
) -> (System, NodeId, NodeId) {
    let (receivers, ports, resistors) = layout(bodies.resistor);
    // One law per receiver, as in graph 4.
    let law = law as Arc<dyn SpatialConductorLaw<ResolvedConductorState>>;
    let winding: PlayWinding = SpatialWinding::at_port_temperatures(
        circuit,
        vec![law; receivers.len()],
        resistors,
        receivers,
        controls(),
        controls().energy_tolerance_j,
    )
    .unwrap();
    let start = &bodies.start;
    graph(
        winding,
        PlayState::new(circuit.clone()).unwrap(),
        timeline(start, conduction(start, &bodies.materials, &[])),
        &ports,
        Some(drive),
    )
}

/// The accepted interval the system holds.
fn row(system: &System, winding: NodeId, cells: NodeId) -> (Row, Vec<f64>) {
    let state = system.state::<PlayState>(winding).unwrap();
    let timeline = system.state::<TimelineState<'static>>(cells).unwrap();
    play_row(
        kelvins(timeline.simulation()),
        state.circuit(),
        state.step().unwrap(),
    )
}

fn kernel_run(system: &mut System, winding: NodeId, cells: NodeId) -> (Vec<Row>, Vec<Vec<f64>>) {
    with_cx(|cx| {
        (0..STEPS)
            .map(|_| {
                let report = system.step(DT, cx).unwrap();
                assert_eq!(
                    report.dt_s.to_bits(),
                    DT.to_bits(),
                    "the kernel took the whole interval"
                );
                row(system, winding, cells)
            })
            .unzip()
    })
}

/// The module's run: test 5's `step_hysteretic_candidate`, or with
/// `resistor`, test 6's `step_hysteretic_material_candidate` with the
/// resistor in its own body, from `circuit` under `drive`; its conductors
/// held to `law`, and its graph resistor to `resistor`, at its own end
/// temperatures, bit for bit.
fn module_run(
    card: &MaterialCard,
    curve: &ResolvedHeatCapacityCurve,
    circuit: &SpatialHystereticWindingCircuit<'_>,
    (law, resistor): (&CardConductors, Option<&BodyResistor>),
    drive: [f64; STEPS],
) -> (Vec<Row>, Vec<Vec<f64>>) {
    let query = material::point(300.0);
    let separate = resistor.is_some();
    let cells = module_cells(card, curve, &query, separate);
    let regions = regions();
    let routes = [ResistorCell {
        branch: ROUTES[0].branch,
        cell: if separate { RES_CELL } else { CELL },
    }];
    let resistors = resistor.map(|r| CellResistor {
        destination: routes[0],
        selection: MaterialPropertySelection::SingleClaimOnly,
        length_m: r.length_m,
        area_m2: r.area_m2,
    });
    let network = SpatialWindingThermalNetwork {
        magnetic: SpatialMagneticThermalNetwork {
            cells: &cells,
            links: &[],
            magnetic_regions: &regions,
            volume_tolerance_m3: 1e-12,
        },
        graph_resistors: &routes,
        fill_tolerance: 1e-12,
    };
    let (mut heat, mut accepted, mut target) = (thermal(curve, separate), circuit.clone(), 0.0);
    let mut rows = (Vec::new(), Vec::new());
    for (interval, input) in drive.into_iter().enumerate() {
        target += DT;
        let next = match &resistors {
            None => network.step_hysteretic_candidate(
                &heat,
                target,
                &accepted,
                &[input],
                thermal_controls(),
                controls(),
                &[0.0],
                || false,
            ),
            Some(cell) => network.step_hysteretic_material_candidate(
                &heat,
                target,
                &accepted,
                &[input],
                std::slice::from_ref(cell),
                &[(ROUTES[0].branch, BASIS)],
                thermal_controls(),
                controls(),
                &[0.0; 2],
                || false,
            ),
        }
        .unwrap();
        let (t, c) = (kelvins(next.thermal()), next.circuit());
        let conductors = (c.conductors(), c.field().conductors());
        follows(law, conductors, t[0], interval);
        let step = &next.report().coupled.electrical;
        if let Some(resistor) = resistor {
            assert_eq!(
                resistor.conductor(t[1]).unwrap(),
                step.graph_materials[0].conductor,
                "interval {interval}"
            );
        }
        let (row, memory) = play_row(t, c, step);
        rows.0.push(row);
        rows.1.push(memory);
        (heat, accepted, _) = next.into_parts();
    }
    rows
}

/// The largest difference between two runs' play memories \[T\].
fn memory_apart(kernel: &[Vec<f64>], module: &[Vec<f64>]) -> f64 {
    let pairs = kernel.iter().flatten().zip(module.iter().flatten());
    assert_eq!(kernel.concat().len(), module.concat().len());
    pairs.map(|(a, b)| (a - b).abs()).fold(0.0, f64::max)
}

/// Bounds of [`largest`] and [`memory_apart`], in that order, for graphs 5
/// and 6. On d94e1ed746's export graph 5's kernel and module agree within
/// 1.327e-10 K, 1.502e-11 A, 0 (the circuit state to the bit), 1.531e-14 Wb,
/// 1.651e-10 J and 1.099e-13 T; graph 6's within 4.604e-12 K, 1.080e-12 A,
/// 2.842e-14, 2.165e-15 Wb, 1.285e-11 J and 1.643e-14 T. The bounds sit 48
/// to 53 times above the larger of the two.
const PLAY: [f64; 6] = [7e-9, 8e-10, 1.5e-12, 8e-13, 8e-9, 5.5e-12];

/// The kernel against the module: each quantity within its bound.
fn agrees(name: &str, kernel: &(Vec<Row>, Vec<Vec<f64>>), module: &(Vec<Row>, Vec<Vec<f64>>)) {
    let [k, a, s, w, j] = largest(&kernel.0, &module.0);
    let most = [k, a, s, w, j, memory_apart(&kernel.1, &module.1)];
    println!(
        "{name}: kernel against module: {k:.3e} K, {a:.3e} A, {s:.3e}, {w:.3e} Wb, {j:.3e} J, memory {:.3e} T",
        most[5]
    );
    for (i, (value, bound)) in most.iter().zip(PLAY).enumerate() {
        assert!(*value <= bound, "{name}, quantity {i}: {most:?}");
    }
}

/// Graph 5: test 5's intervals at the graph's own drive, the graph
/// resistor heating the coil's cell. Test 5's scalar reference holds on the
/// kernel's rows: the cell's temperature solves `C dT` = the core's eddy
/// heat, the pinning heat and the strand and resistor heat; the current is
/// the voltage play's; the field carries the play's memory.
#[test]
fn the_declared_play_memory_graph_follows_the_scalar_play_and_its_four_losses() {
    let (card, curve, p) = fixture();
    let start = circuit(p, &material::conductor(card, 300.0));
    let law = card_law(card, &start, None);
    let bodies = Bodies {
        start: thermal(curve, false),
        materials: vec![(MAT, curve)],
        resistor: None,
    };
    let drive = |t0: f64| DRIVE_5[usize::from(t0 >= 0.5 * DT)];
    let (mut system, winding, cells) = kernel(bodies, &start, law.clone(), drive);
    let rows = kernel_run(&mut system, winding, cells);
    ledger_closes(&system, winding, cells, 0.0);
    let (mut old, mut before) = (Mode::default(), 300.0);
    for (n, (row, input)) in rows.0.iter().zip(DRIVE_5).enumerate() {
        let (t, i) = (row.kelvin[0], row.currents_a[0]);
        let physical = |t| {
            let next = voltage_endpoint(old, t, input, GRAPH_R);
            core_heat(old, next, t)
                + pinning_heat(old, next)
                + DT * (8.0 * material::resistivity(t) + GRAPH_R) * next.i.powi(2)
        };
        let reference = material::bracket(before, 480.0, |t| CAPACITY * (t - before) - physical(t));
        let expected = voltage_endpoint(old, t, input, GRAPH_R);
        let [strand, field, resistor, heat, ..] = row.energies_j;
        println!("graph 5, interval {n}: {t:.9} K against {reference:.9} K, {i:.9} A");
        near(t, reference, 3e-7);
        near(i, expected.i, 4e-7);
        let losses = core_heat(old, expected, t) + pinning_heat(old, expected);
        near(field, losses, 6e-8);
        near(strand, DT * 8.0 * material::resistivity(t) * i * i, 4e-7);
        near(resistor, DT * GRAPH_R * i * i, 4e-7);
        near(heat, physical(t), 5e-7);
        (old, before) = (expected, t);
    }
    // The second interval releases the drive, and the play keeps its memory.
    let state = system.state::<PlayState>(winding).unwrap();
    assert_field(state.circuit().field(), old);
    assert!(old.z > material::play(0.0, old.h) + 0.1);
    let module = module_run(card, curve, &start, (&law, None), DRIVE_5);
    agrees("graph 5", &rows, &module);
}

/// Graph 6: test 6's graph resistor in its own body, its law the body's card
/// at test 6's length and area. Test 6's equations hold on the kernel's
/// rows; the module refuses a source at 301 K, which the kernel has no input
/// to express; and a late cancellation publishes nothing.
#[test]
fn the_declared_play_memory_graph_heats_a_separate_resistor_body_and_refuses_a_late_cancellation() {
    let (card, curve, p) = fixture();
    let source = material::conductor(card, 300.0);
    let length = (GRAPH_R / material::resistivity(300.0)).sqrt();
    let area = 1.0 / length;
    let original = circuit(p, &source)
        .with_graph_resistors(
            vec![GraphResistorMaterial {
                branch: ROUTES[0].branch,
                conductor: source.into_uniform_conductor(length, area).unwrap(),
            }],
            &mut || false,
        )
        .unwrap();
    let resistor = Arc::new(BodyResistor::sized(card, material::point, length, area));
    let bodies = || Bodies {
        start: thermal(curve, true),
        materials: vec![(MAT, curve)],
        resistor: Some(resistor.clone() as Arc<dyn ResistorLaw>),
    };
    let drive = |t0: f64| DRIVE_6[usize::from(t0 >= 0.5 * DT)];
    let law = card_law(card, &original, None);
    let (mut system, winding, cells) = kernel(bodies(), &original, law.clone(), drive);
    let rows = kernel_run(&mut system, winding, cells);
    ledger_closes(&system, winding, cells, 0.0);
    let (mut old, mut before) = (Mode::default(), [300.0; 2]);
    for (n, (row, input)) in rows.0.iter().zip(DRIVE_6).enumerate() {
        let ([t, tr], i) = ([row.kelvin[0], row.kelvin[1]], row.currents_a[0]);
        let r = material::resistivity(tr) * length / area;
        let expected = voltage_endpoint(old, t, input, r);
        let [strand, field, graph, ..] = row.energies_j;
        println!("graph 6, interval {n}: coil {t:.9} K, body {tr:.9} K, {i:.9} A");
        near(i, expected.i, 4e-7);
        near(strand, DT * 8.0 * material::resistivity(t) * i * i, 4e-7);
        near(CAPACITY * (t - before[0]), field + strand, 5e-7);
        near(CAPACITY * (tr - before[1]), DT * r * i * i, 5e-7);
        near(graph, DT * r * i * i, 4e-7);
        assert_eq!(row.resistance_ohm.len(), 1, "interval {n}");
        near(row.resistance_ohm[0], r, 1e-12);
        assert!(n > 0 || tr > t && t > 300.0);
        (old, before) = (expected, [t, tr]);
    }
    let module = module_run(card, curve, &original, (&law, Some(&*resistor)), DRIVE_6);
    agrees("graph 6", &rows, &module);

    // A construction circuit whose resistor holds the card at 301 K: the
    // module's `admit_initial_graph_sources` refuses it; the kernel rebinds
    // every graph resistor from its law before the solve, so it steps to the
    // original's bits.
    let wrong = original
        .with_graph_resistors(
            vec![GraphResistorMaterial {
                branch: ROUTES[0].branch,
                conductor: material::conductor(card, 301.0)
                    .into_uniform_conductor(length, area)
                    .unwrap(),
            }],
            &mut || false,
        )
        .unwrap();
    let law = card_law(card, &wrong, None);
    let (mut system, winding, cells) = kernel(bodies(), &wrong, law, drive);
    with_cx(|cx| system.step(DT, cx)).unwrap();
    let (first, memory) = row(&system, winding, cells);
    assert_eq!(first.bits(), rows.0[0].bits());
    assert_eq!(memory, rows.1[0]);
    cancelled(card, &original, |law| {
        kernel(bodies(), &original, law, drive)
    });
}

/// The kernel's late cancellation of the first interval: a good step counts
/// its law's answers; the same step, cancelled at the last answer, publishes
/// nothing, and its retry is the good step, bit for bit.
fn cancelled(
    card: &MaterialCard,
    original: &SpatialHystereticWindingCircuit<'static>,
    built: impl Fn(Arc<CardConductors>) -> (System, NodeId, NodeId),
) {
    let counted = card_law(card, original, None);
    let (mut system, winding, cells) = built(counted.clone());
    with_cx(|cx| system.step(DT, cx)).unwrap();
    let good = row(&system, winding, cells);
    let answers = counted.answers.load(Ordering::Relaxed);
    let gate: &'static CancelGate = Box::leak(Box::new(CancelGate::new_clock_free()));
    let tripping = card_law(card, original, Some((answers, gate)));
    let (mut system, winding, cells) = built(tripping.clone());
    let error = with_gate(gate, |cx| system.step(DT, cx)).unwrap_err();
    assert!(matches!(error, SystemError::Cancelled), "{error:?}");
    assert_eq!(tripping.answers.load(Ordering::Relaxed), answers);
    assert_eq!(system.time_s().to_bits(), 0.0_f64.to_bits());
    let state = system.state::<PlayState>(winding).unwrap();
    assert!(*state == PlayState::new(original.clone()).unwrap() && state.step().is_none());
    let timeline = system.state::<TimelineState<'static>>(cells).unwrap();
    for kelvin in kelvins(timeline.simulation()) {
        assert_eq!(kelvin.to_bits(), 300.0_f64.to_bits());
    }
    with_cx(|cx| system.step(DT, cx)).unwrap();
    let retry = row(&system, winding, cells);
    assert_eq!((retry.0.bits(), retry.1), (good.0.bits(), good.1));
    println!("graph 6: {answers} conductor answers in the step; the cancelled one was the last");
}
