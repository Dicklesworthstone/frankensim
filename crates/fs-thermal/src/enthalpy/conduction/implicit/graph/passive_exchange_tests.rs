//! A coupled source that declares itself a closed passive exchange
//! (`ImplicitCoupledHeatSourceLaw::passive_exchange_tolerance_w`).
use super::*;
use crate::enthalpy::advection::{Advected, MassFlow};

/// `G (T_b - T_a)` moved into cell `a` out of cell `b` for each pair, plus
/// fixed powers `extra_w`, one per cell, on the original law's calorics.
struct Exchange<'a> {
    law: &'a Law,
    calorics: Vec<&'a dyn ThermalLaw>,
    pairs: Vec<(usize, usize, f64)>,
    extra_w: Vec<f64>,
    dependencies: Vec<(usize, usize)>,
    tolerance_w: Option<f64>,
}

fn exchange<'a>(
    law: &'a Law,
    pairs: &[(usize, usize, f64)],
    extra_w: &[f64],
    tolerance_w: Option<f64>,
) -> Exchange<'a> {
    let n = extra_w.len();
    let dependencies = if pairs.is_empty() {
        Vec::new()
    } else {
        (0..n)
            .flat_map(|row| (0..n).map(move |col| (row, col)))
            .collect()
    };
    Exchange {
        law,
        calorics: vec![law as &dyn ThermalLaw; n],
        pairs: pairs.to_vec(),
        extra_w: extra_w.to_vec(),
        dependencies,
        tolerance_w,
    }
}

impl ImplicitCoupledHeatSourceLaw for Exchange<'_> {
    fn caloric_laws(&self) -> &[&dyn ThermalLaw] {
        &self.calorics
    }
    fn power_dependencies(&self) -> &[(usize, usize)] {
        &self.dependencies
    }
    fn state_at_specific_enthalpies_and_time(
        &self,
        h: &[f64],
        _time_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitCoupledHeatSourceState, DepositionError> {
        checkpoint()?;
        let n = self.extra_w.len();
        let mut powers = self.extra_w.clone();
        let mut slopes = vec![vec![0.0; n]; n];
        for &(a, b, g) in &self.pairs {
            let states = [
                self.law.state_at_specific_enthalpy(h[a])?,
                self.law.state_at_specific_enthalpy(h[b])?,
            ];
            let q = g * (states[1].temperature_k - states[0].temperature_k);
            powers[a] += q;
            powers[b] -= q;
            for (row, sign) in [(a, 1.0), (b, -1.0)] {
                slopes[row][a] -= sign * g * states[0].temperature_slope_kg_k_j;
                slopes[row][b] += sign * g * states[1].temperature_slope_kg_k_j;
            }
        }
        Ok(ImplicitCoupledHeatSourceState {
            powers_w: powers,
            power_enthalpy_derivatives_w_kg_j: self
                .dependencies
                .iter()
                .map(|&(row, col)| slopes[row][col])
                .collect(),
        })
    }
    fn passive_exchange_tolerance_w(&self) -> Option<f64> {
        self.tolerance_w
    }
}

fn attached<'a>(
    thermal: &DepositionSimulation<'_>,
    observer: &ConductionObserver<'a>,
    law: &'a Law,
    source: &'a Exchange<'a>,
) -> Result<ImplicitConductionGraph<'a>, ImplicitConductionError> {
    let cells: Vec<CellId> = (0..source.extra_w.len()).map(id).collect();
    graph(thermal, observer, law).with_coupled_heat_source(
        thermal,
        &cells,
        source,
        source.dependencies.len(),
        || Ok(()),
    )
}

fn step<'m>(
    graph: &ImplicitConductionGraph<'_>,
    thermal: &DepositionSimulation<'m>,
    dt: f64,
    sdirk: bool,
) -> Result<(DepositionSimulation<'m>, ImplicitConductionReport), ImplicitConductionError> {
    let mut workspace = ImplicitConductionWorkspace::new(graph.cell_ids().len()).unwrap();
    if sdirk {
        graph.advance_sdirk2_to(thermal, dt, controls(dt), &mut workspace, || Ok(()))
    } else {
        graph.advance_to(thermal, dt, controls(dt), &mut workspace, || Ok(()))
    }
}

fn refused<T>(result: Result<T, ImplicitConductionError>, words: &str) {
    match result {
        Err(
            ImplicitConductionError::Invalid(message)
            | ImplicitConductionError::Thermal(DepositionError::InvalidInput(message)),
        ) => assert!(message.contains(words), "{message}"),
        Err(other) => panic!("expected a refusal naming {words}, got {other:?}"),
        Ok(_) => panic!("expected a refusal naming {words}"),
    }
}

#[test]
fn a_declared_exchange_takes_exactly_its_own_entropy_out_of_the_unclassified_sources() {
    // Cells at 300, 340 and 380 K, cell 0 linked to cell 1, 1.5 W authored
    // into cell 2; the exchange moves heat from cell 2 into cells 0 and 1.
    // The same law steps once declared and once not.
    let law = Law {
        plateau: false,
        nonlinear: false,
    };
    let (mut thermal, observer) = setup(
        &law,
        &[300.0, 340.0, 380.0],
        &[2.0, 3.0, 4.0],
        &[(0, 1, 0.0)],
    );
    thermal.cells.get_mut(&id(2)).unwrap().boundary.source_w = 1.5;
    let pairs = [(0, 2, 0.7), (1, 2, 0.4)];
    let declared = exchange(&law, &pairs, &[0.0; 3], Some(1e-12));
    let plain = exchange(&law, &pairs, &[0.0; 3], None);
    let dt = 0.5;
    for sdirk in [false, true] {
        let graph = attached(&thermal, &observer, &law, &declared).unwrap();
        let (next, report) = step(&graph, &thermal, dt, sdirk).unwrap();
        let graph = attached(&thermal, &observer, &law, &plain).unwrap();
        let (plain_next, plain_report) = step(&graph, &thermal, dt, sdirk).unwrap();
        // The classification moves no heat and changes no other term.
        for i in 0..3 {
            let [a, b] = [&next, &plain_next].map(|state| state.cells()[&id(i)].specific_h_j_kg);
            assert_eq!(a.to_bits(), b.to_bits());
        }
        assert_eq!(report.energy, plain_report.energy);
        let (a, b) = (report.entropy, plain_report.entropy);
        for (x, y) in [
            (a.thermal_storage_change_j_k, b.thermal_storage_change_j_k),
            (a.reservoir_entropy_input_j_k, b.reservoir_entropy_input_j_k),
            (a.conduction_generation_j_k, b.conduction_generation_j_k),
        ] {
            assert_eq!(x.to_bits(), y.to_bits());
        }
        assert_eq!(
            b.passive_exchange_generation_j_k.to_bits(),
            0.0_f64.to_bits()
        );
        assert!(a.passive_exchange_generation_j_k > 0.0);
        assert!(a.has_unclassified_sources && b.has_unclassified_sources);
        // Its two parts make the plain unclassified term, and the identity's
        // defect is the plain one, to rounding.
        let parts = a.unclassified_source_at_cell_j_k + a.passive_exchange_generation_j_k
            - b.unclassified_source_at_cell_j_k;
        let defect = a.temporal_defect_j_k - b.temporal_defect_j_k;
        println!(
            "declared exchange, {}: production {:.6e} J/K, unclassified {:.6e} J/K (plain {:.6e}); parts less plain {parts:.1e}, defect change {defect:.1e} (defect {:.1e})",
            if sdirk { "SDIRK2" } else { "BE" },
            a.passive_exchange_generation_j_k,
            a.unclassified_source_at_cell_j_k,
            b.unclassified_source_at_cell_j_k,
            b.temporal_defect_j_k,
        );
        assert!(parts.abs() <= 1e-15 && defect.abs() <= 1e-15);
        if !sdirk {
            // Backward Euler samples at the end temperatures: by hand from
            // the law's own powers there.
            let h: Vec<f64> = (0..3)
                .map(|i| next.cells()[&id(i)].specific_h_j_kg)
                .collect();
            let t: Vec<f64> = (0..3).map(|i| next.temperature_k(id(i)).unwrap()).collect();
            let powers = declared
                .state_at_specific_enthalpies_and_time(&h, dt, &mut || Ok(()))
                .unwrap()
                .powers_w;
            let by_hand = dt * (powers[0] / t[0] + powers[1] / t[1] + powers[2] / t[2]);
            let authored = dt * 1.5 / t[2];
            let errors = [
                a.passive_exchange_generation_j_k - by_hand,
                a.unclassified_source_at_cell_j_k - authored,
            ];
            println!(
                "  by hand at the end temperatures: {:.1e}, {:.1e}",
                errors[0], errors[1]
            );
            assert!(errors[0].abs() <= 1e-15 * by_hand && errors[1].abs() <= 1e-15 * authored);
        }
    }
}

#[test]
fn a_declared_exchange_is_refused_beyond_its_closure_and_beside_an_implicit_flow() {
    let law = Law {
        plateau: false,
        nonlinear: false,
    };
    let (thermal, observer) = setup(&law, &[300.0, 340.0], &[2.0, 3.0], &[(0, 1, 0.0)]);
    let pairs = [(0, 1, 0.5)];
    for tolerance_w in [-1e-12, f64::NAN, f64::INFINITY] {
        let source = exchange(&law, &pairs, &[0.0; 2], Some(tolerance_w));
        refused(
            attached(&thermal, &observer, &law, &source),
            "finite nonnegative",
        );
    }
    // 2e-9 W left over against a 1e-9 W tolerance is refused at the query
    // the attachment makes; 5e-10 W is inside it, and the step goes ahead.
    let open = exchange(&law, &pairs, &[2e-9, 0.0], Some(1e-9));
    refused(
        attached(&thermal, &observer, &law, &open),
        "exceed their closure tolerance",
    );
    let inside = exchange(&law, &pairs, &[5e-10, 0.0], Some(1e-9));
    step(
        &attached(&thermal, &observer, &law, &inside).unwrap(),
        &thermal,
        0.5,
        false,
    )
    .unwrap();
    // Beside an implicit flow, in either order; undeclared, the same source
    // attaches.
    let flows = [(10, 0, 1), (11, 1, 0)].map(|(flow, a, b)| MassFlow {
        id: flow,
        a: id(a),
        b: id(b),
        kg_s: 0.5,
    });
    let flow = Advected::new((), Vec::from(flows), 2).unwrap();
    let declared = exchange(&law, &pairs, &[0.0; 2], Some(1e-12));
    let plain = exchange(&law, &pairs, &[0.0; 2], None);
    let flowing = || {
        graph(&thermal, &observer, &law)
            .with_implicit_advection(&thermal, &flow, || Ok(()))
            .unwrap()
    };
    let cells = [id(0), id(1)];
    refused(
        flowing().with_coupled_heat_source(&thermal, &cells, &declared, 4, || Ok(())),
        "implicit flow",
    );
    refused(
        attached(&thermal, &observer, &law, &declared)
            .unwrap()
            .with_implicit_advection(&thermal, &flow, || Ok(())),
        "implicit flow",
    );
    flowing()
        .with_coupled_heat_source(&thermal, &cells, &plain, 4, || Ok(()))
        .unwrap();
}

#[test]
fn the_entropy_gate_keeps_rounding_below_zero_and_refuses_beyond_its_allowance() {
    // An exactly isothermal enclosure: three cells on the latent plateau at
    // 300 K whose 1, 3 and -4 W close exactly. Their quotients round, and the
    // signed production is that rounding, below zero, inside the band
    // dt * 3 * EPSILON * 8 W / 300 K.
    let plateau = Law {
        plateau: true,
        nonlinear: false,
    };
    let (thermal, observer) = setup(&plateau, &[400.0; 3], &[1.0; 3], &[]);
    let source = exchange(&plateau, &[], &[1.0, 3.0, -4.0], Some(0.0));
    let graph = attached(&thermal, &observer, &plateau, &source).unwrap();
    let dt = 0.5;
    let band = dt * 3.0 * f64::EPSILON * 8.0 / 300.0;
    for sdirk in [false, true] {
        let (next, report) = step(&graph, &thermal, dt, sdirk).unwrap();
        let production = report.entropy.passive_exchange_generation_j_k;
        println!(
            "isothermal exchange, {}: production {production:.3e} J/K, band {band:.3e} J/K",
            if sdirk { "SDIRK2" } else { "BE" },
        );
        assert!(production < 0.0 && -production <= band);
        for i in 0..3 {
            assert_eq!(
                next.temperature_k(id(i)).unwrap().to_bits(),
                300.0_f64.to_bits()
            );
        }
    }
    // 1 W moved from a 300 K cell into a 400 K one closes exactly and loses
    // about dt (1/300 - 1/400) = dt/1200 J/K (the cells move 0.5 K). A 0.2 W
    // tolerance allows about dt/1500 J/K: refused. A 0.3125 W one allows about
    // dt/960 J/K: admitted, the loss kept.
    let law = Law {
        plateau: false,
        nonlinear: false,
    };
    let (thermal, observer) = setup(&law, &[300.0, 400.0], &[1.0; 2], &[]);
    for (tolerance_w, admitted) in [(0.2, false), (0.3125, true)] {
        let source = exchange(&law, &[], &[-1.0, 1.0], Some(tolerance_w));
        let graph = attached(&thermal, &observer, &law, &source).unwrap();
        for sdirk in [false, true] {
            let result = step(&graph, &thermal, dt, sdirk);
            if admitted {
                let (next, report) = result.unwrap();
                let production = report.entropy.passive_exchange_generation_j_k;
                let t = [0, 1].map(|i| next.temperature_k(id(i)).unwrap());
                let by_hand = dt * (-1.0 / t[0] + 1.0 / t[1]);
                println!(
                    "reversed exchange admitted, {}: production {production:.6e} J/K, at the end temperatures {by_hand:.6e} J/K",
                    if sdirk { "SDIRK2" } else { "BE" },
                );
                assert!(production < 0.0 && (production - by_hand).abs() <= 0.01 * -by_hand);
                if !sdirk {
                    assert!((production - by_hand).abs() <= 1e-15 * -by_hand);
                }
            } else {
                refused(result, "beyond its closure and rounding allowance");
            }
        }
    }
}
