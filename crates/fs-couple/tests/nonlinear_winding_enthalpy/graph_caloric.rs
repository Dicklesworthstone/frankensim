//! The module's magnetic thermal trial (`magnetic_thermal`) as a declared
//! graph on the composition kernel.
//!
//! `WindingThermalNetwork::step_nonlinear_thermal_candidate` couples the
//! winding's field to its cell by hand: each law triangle's permeability at
//! the cell's END temperature, and the cell heated by `-P - (T1 s1 - T0 s0)`
//! beside the Joule heat, with `P` the owner's parameter work and `s` the
//! triangles' entropy. The graph is graph.rs's with fs-em's
//! `PlanarWinding::magnetic_at_port_temperatures`: the same laws at the
//! temperature its heat port receives, and a port heat of `-T1 (s1 - s0)`
//! beside the Joule heat. The two heats differ by `-(P + (T1 - T0) s0)`,
//! which the Component books as numerical loss and the module deposits.
//! With the same states the module's cell would end below the graph's by
//! that term summed over the intervals; the cooler cell's resistivity and
//! permeability then move the module's circuit by the term times their
//! temperature slopes.

use std::collections::BTreeMap;
use std::sync::Arc;

use fs_em::planar::nonlinear::PlanarNonlinearMagnetostatics;
use fs_em::winding::component::{
    ConductorLaw, HeatReceiver, MagneticLaw, MagneticLaws, PlanarWinding, PlanarWindingState,
};
use fs_em::winding::nonlinear::NonlinearWindingCircuit;

use super::graph::{CellConductors, Row, STEPS, TOLERANCE_J, graph_run, hold, largest, module_row};
use super::kernel::{graph, ledger_closes};
use super::magnetic_thermal::{self as trial, END_H, HEAT_CAPACITY, VOLTAGE, linear, mu};
use super::*;

/// The trial's closed forms (`reference` in thermal.rs) over an interval
/// from `(i0, t0)` to `(i1, t1)` \[A, K\], with `b = mu(t) i` in every
/// triangle: the owner's backward-Euler dissipation and the term
/// `-(P + (t1 - t0) s0)` \[J\], `s0 = -dpsi/dT` at the start.
fn closed((i0, t0): (f64, f64), (i1, t1): (f64, f64), slope: f64) -> [f64; 2] {
    let (m0, m1) = (mu(t0, slope), mu(t1, slope));
    let b0 = m0 * i0;
    let numerical = VOLUME * (m1 * i1 - b0).powi(2) / (2.0 * m1) + 0.5 * END_H * (i1 - i0).powi(2);
    let parameter = VOLUME * b0 * b0 * (1.0 / m1 - 1.0 / m0) / 2.0;
    let s0 = VOLUME * slope * b0 * b0 / (2.0 * m0 * m0);
    [numerical, -(parameter + (t1 - t0) * s0)]
}

/// Bounds at slope 0, beside graph.rs's `BOUNDS`: the end field potential
/// \[Wb/m\], the booked loss against the term \[J\] and the cell gap against
/// the summed term \[J\]. On the export of aa4d699b1f they measure 1.1e-11,
/// 1.5e-10 and 3.3e-10; the bounds sit about fifty times above.
const CONTROL: [f64; 3] = [5e-10, 7.5e-9, 1.6e-8];

/// Bounds at slope 0.004: temperature \[K\], current \[A\] and end field
/// potential \[Wb/m\] against the module, the booked loss against the term
/// \[J\], and the cell gap's departure from the summed term as a fraction of
/// it. They measure 4.0e-2, 7.9e-4, 1.9e-4, 1.7e-11 and 1.8 %; all but the
/// booked loss come from the cooler cell's feedback through its resistivity
/// and permeability.
const FEEDBACK: [f64; 5] = [5e-2, 1e-3, 2.5e-4, 7.5e-9, 0.025];

/// The largest difference of two field potentials \[Wb/m\].
fn farthest(ours: &[f64], theirs: &[f64]) -> f64 {
    ours.iter()
        .zip(theirs)
        .fold(0.0, |most, (x, y)| most.max((x - y).abs()))
}

/// The module's run of the trial: `step_nonlinear_thermal_candidate` from
/// `circuit`, its conductors held to `law` at its own end temperature, bit
/// for bit. Also its field potential \[Wb/m\] at the end.
fn caloric_module(
    card: &MaterialCard,
    curve: &ResolvedHeatCapacityCurve,
    problem: &PlanarNonlinearMagnetostatics,
    circuit: &NonlinearWindingCircuit<'_>,
    law: &CellConductors,
) -> (Vec<Row>, Vec<f64>) {
    let mut thermal = trial::thermal(curve);
    let mut circuit = circuit.clone();
    let bindings = trial::laws(card, problem);
    let context = linear::absolute_point(300.0);
    let cells = cells(card, curve, &context);
    let regions = routes();
    let layout = network(&cells, &regions);
    let mut rows = Vec::new();
    for step in 1..=STEPS {
        let candidate = layout
            .step_nonlinear_thermal_candidate(
                &thermal,
                step as f64 * DT,
                &circuit,
                &[VOLTAGE],
                &bindings,
                None,
                controls(),
                field::controls(),
                linear::limits(),
                &[0.0],
                || false,
            )
            .unwrap();
        let t = candidate.thermal().temperature_k(CELL).unwrap();
        let electrical = &candidate.report().electrical;
        for loss in &electrical.cells {
            let triangle = loss.source.triangle;
            assert_eq!(
                law.conductor(triangle, t).unwrap(),
                loss.source.conductor,
                "step {step}, triangle {triangle}"
            );
        }
        rows.push(module_row(
            vec![t],
            candidate.circuit(),
            electrical,
            [electrical.heat_j / t, 0.0],
        ));
        (thermal, circuit, _) = candidate.into_parts();
    }
    (rows, circuit.field().potential_wb_m().to_vec())
}

#[test]
fn the_declared_graph_and_the_magnetic_trial_differ_by_the_booked_term() {
    for slope in [0.0, 0.004] {
        let card: &'static MaterialCard = Box::leak(Box::new(trial::card(slope)));
        let curve: &'static ResolvedHeatCapacityCurve = Box::leak(Box::new(trial::caloric(card)));
        let problem: &'static PlanarNonlinearMagnetostatics =
            Box::leak(Box::new(trial::problem(card)));
        let circuit = trial::circuit(problem, 0.0);
        let mut law = CellConductors::new(card, &circuit);
        law.context = linear::absolute_point(300.0);
        let law = Arc::new(law);
        let (module, potential) = caloric_module(card, curve, problem, &circuit, &law);

        let magnetic = MagneticLaws {
            triangles: trial::laws(card, problem)
                .into_iter()
                .map(|binding| MagneticLaw {
                    triangle: binding.triangle,
                    law: binding.law,
                })
                .collect(),
            assembly: linear::limits(),
        };
        let winding = PlanarWinding::magnetic_at_port_temperatures(
            &circuit,
            vec![law],
            BTreeMap::new(),
            magnetic,
            vec![HeatReceiver {
                triangles: (0..8).collect(),
                resistors: vec![2],
            }],
            field::controls(),
            TOLERANCE_J,
        )
        .unwrap();
        let start = trial::thermal(curve);
        let (mut system, electrical, cells) = graph(
            winding,
            PlanarWindingState::at_temperatures(circuit, vec![300.0]),
            (&start, &[(MATERIAL, curve)], &[]),
            &[CELL],
            Some(|_| VOLTAGE),
        );
        let kernel = graph_run(&mut system, electrical, cells, &[CELL]);
        ledger_closes(&system, electrical, cells);
        let entropy = system.ledger().nodes[electrical.index()]
            .entropy
            .clone()
            .unwrap();
        assert_eq!(entropy.audited_steps, STEPS as u64);

        // Each interval from the kernel's own start: its numerical loss less
        // the owner's dissipation against the term, and the module's cell
        // energy less the kernel's against the terms summed.
        let (mut from, mut summed, mut terms) = ((0.0, 300.0), 0.0, Vec::new());
        let mut most = [0.0_f64; 2];
        for (ours, theirs) in kernel.iter().zip(&module) {
            let to = (ours.current_a, ours.kelvin[0]);
            let [numerical, term] = closed(from, to, slope);
            most[0] = most[0].max((ours.numerical_j - numerical - term).abs());
            summed += term;
            let gap = HEAT_CAPACITY * (theirs.kelvin[0] - ours.kelvin[0]);
            most[1] = most[1].max((gap - summed).abs());
            terms.push(term);
            from = to;
        }
        let differences = largest(&kernel, &module);
        let state = system
            .state::<PlanarWindingState<'static>>(electrical)
            .unwrap();
        let field = farthest(state.circuit().field().potential_wb_m(), &potential);
        println!(
            "graph_caloric (slope {slope}): kernel against module: {:.3e} K, {:.3e} A, {:.3e} Wb/m at the end; \
             booked against the term {:.3e} J; cell gap against the summed term {:.3e} J; \
             terms {terms:?} J, summed {summed:e} J; end {:.6} K, module {:.6} K",
            differences[0],
            differences[1],
            field,
            most[0],
            most[1],
            kernel[STEPS - 1].kelvin[0],
            module[STEPS - 1].kelvin[0]
        );
        if slope == 0.0 {
            hold(differences, "slope 0");
            for (value, bound) in [field, most[0], most[1]].into_iter().zip(CONTROL) {
                assert!(value <= bound, "slope 0: {value:e} above {bound:e}");
            }
        } else {
            let values = [
                differences[0],
                differences[1],
                field,
                most[0],
                most[1] / summed.abs(),
            ];
            for (value, bound) in values.into_iter().zip(FEEDBACK) {
                assert!(value <= bound, "slope {slope}: {value:e} above {bound:e}");
            }
            // The term is negative from the second interval (the first
            // starts with no field), and the module's cell ends cooler.
            assert!(terms[1..].iter().all(|term| *term < -1e-3), "{terms:?}");
            assert!(module[STEPS - 1].kelvin[0] < kernel[STEPS - 1].kelvin[0]);
        }
    }
}
