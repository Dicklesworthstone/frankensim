//! G0/G1/G3/G4 checks of the original liquid carrier in native implicit solves.
use super::*;
use crate::enthalpy::advection::liquid::{LiquidCarrierKnot, LiquidEnthalpyCarrier};

fn phase(datum: f64, nonlinear_k: bool) -> ShiftedLaw {
    ShiftedLaw {
        base: Law {
            nonlinear_k,
            plateau: true,
        },
        datum,
    }
}

fn carrier(law: &ShiftedLaw) -> LiquidEnthalpyCarrier {
    let knots = [
        (100.0, 0.0),
        (300.0, 0.0),
        (400.0, 0.5),
        (500.0, 1.0),
        (1000.0, 1.0),
    ]
    .map(|(h, liquid_mass_fraction)| LiquidCarrierKnot {
        specific_enthalpy_j_kg: h + law.datum,
        temperature_k: law.temperature_k(h + law.datum).unwrap(),
        liquid_mass_fraction,
        bulk_density_kg_m3: 1.0,
    });
    LiquidEnthalpyCarrier::from_source_chart(MaterialId(0), law, &knots, || Ok(())).unwrap()
}

fn solver<'a>(
    thermal: &DepositionSimulation<'_>,
    conduction: &ConductionObserver<'a>,
    law: &'a ShiftedLaw,
    flow: &Advected<()>,
    graph: bool,
) -> ImplicitConductionChain<'a> {
    ImplicitConductionChain::prepare(
        thermal,
        conduction,
        &[(MaterialId(0), law)],
        if graph {
            LinearTopology::Graph(ImplicitConductionLinearControls::default())
        } else {
            LinearTopology::Chain
        },
        || Ok(()),
    )
    .unwrap()
    .with_implicit_advection(thermal, flow, || Ok(()))
    .unwrap()
}

#[test]
fn implicit_liquid_isothermal_circulation_preserves_stationary_solid() {
    for datum in [0.0, 1_099_511_627_776.0] {
        let law = phase(datum, true);
        let (thermal, conduction) = setup(&law, &[350.0 + datum, 450.0 + datum], None, vec![]);
        // A Courant number of 1000 cannot force artificial mixing or thousands
        // of explicit substeps when the moving liquid is already isothermal.
        let flow = circulating((), 1000.0)
            .with_liquid_carrier(carrier(&law))
            .unwrap();
        let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
        for graph in [false, true] {
            let system = solver(&thermal, &conduction, &law, &flow, graph);
            for second_order in [false, true] {
                let (next, report) = if second_order {
                    system.advance_sdirk2_to(
                        &thermal,
                        1.0,
                        controls(1.0),
                        &mut workspace,
                        || Ok(()),
                    )
                } else {
                    system.advance_to(&thermal, 1.0, controls(1.0), &mut workspace, || Ok(()))
                }
                .unwrap();
                assert_eq!(next.cells(), thermal.cells());
                assert_eq!(report.energy.substeps, 1);
                assert_eq!(report.energy.external_heat_j, 0.0);
                assert_eq!(report.energy.stored_energy_change_j, 0.0);
            }
        }
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(flow.audit().substeps, 0);
    }
}

#[test]
fn implicit_liquid_superheat_melts_solid_and_closes_the_same_conduction_balance() {
    for datum in [0.0, 1_073_741_824.0] {
        let law = phase(datum, false);
        let (thermal, conduction) = setup(&law, &[600.0 + datum, 350.0 + datum], None, vec![]);
        let flow = circulating((), 3.0)
            .with_liquid_carrier(carrier(&law))
            .unwrap();
        let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
        for graph in [false, true] {
            let system = solver(&thermal, &conduction, &law, &flow, graph);
            let (next, report) = system
                .advance_to(&thermal, 1.0, controls(1.0), &mut workspace, || Ok(()))
                .unwrap();
            // Cp=m=G=1, mdot=3. The cold cell stays at Tm, so hot superheat
            // solves x_new=100-(1+3)*x_new. The SAME 80 J melts the cold cell.
            near(
                next.cells()[&CellId(0)].specific_h_j_kg - datum,
                520.0,
                1e-9,
            );
            near(
                next.cells()[&CellId(1)].specific_h_j_kg - datum,
                430.0,
                1e-9,
            );
            near(next.temperature_k(CellId(1)).unwrap(), 300.0, 0.0);
            near(report.energy.external_heat_j, 0.0, 0.0);
            near(report.energy.stored_energy_change_j, 0.0, 1e-9);
            near(report.energy.residual_j, 0.0, 1e-9);
            assert_eq!(report.energy.substeps, 1);
            assert!(report.entropy.has_unclassified_sources);
        }
    }
}

#[test]
fn implicit_liquid_open_heating_crosses_the_liquidus_and_counts_only_stream_heat() {
    let law = phase(0.0, false);
    let (thermal, conduction) = setup(&law, &[400.0], None, vec![]);
    let flow = open_flow(1, 2.0, 600.0)
        .with_liquid_carrier(carrier(&law))
        .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    for graph in [false, true] {
        let system = solver(&thermal, &conduction, &law, &flow, graph);
        let (next, report) = system
            .advance_to(&thermal, 1.0, controls(1.0), &mut workspace, || Ok(()))
            .unwrap();
        let h = (400.0 + 2.0 * 600.0) / 3.0;
        near(next.cells()[&CellId(0)].specific_h_j_kg, h, 1e-9);
        near(report.energy.external_heat_j, 2.0 * (600.0 - h), 1e-9);
        near(report.energy.stored_energy_change_j, h - 400.0, 1e-9);
        near(report.energy.residual_j, 0.0, 1e-9);
        assert_eq!(report.energy.substeps, 1);
    }
    assert_eq!(flow.audit().substeps, 0);
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn implicit_liquid_tangent_tracks_trial_phase_and_unequal_cell_mass() {
    let law = phase(0.0, true);
    let (mut thermal, _) = setup(&law, &[350.0, 650.0, 540.0], None, vec![]);
    for (i, mass) in [1.0, 2.0, 1.5].into_iter().enumerate() {
        let cell = thermal.cells.get_mut(&CellId(i as u64)).unwrap();
        cell.mass_kg = mass;
        cell.geometry.volume_m3 = mass;
    }
    let links = [0, 1].map(|link_id| ConductiveLink {
        link_id,
        half_distances_m: [0.5; 2],
        contact_resistance_m2_k_w: 0.0,
    });
    let conduction =
        ConductionObserver::new(&thermal, &[(MaterialId(0), &law)], &links, &[], || Ok(()))
            .unwrap();
    for rate in [-2.5, 2.5] {
        let flow = Advected::new(
            (),
            (0..3)
                .map(|i| MassFlow {
                    id: i,
                    a: CellId(i),
                    b: CellId((i + 1) % 3),
                    kg_s: rate,
                })
                .collect(),
            3,
        )
        .unwrap()
        .with_liquid_carrier(carrier(&law))
        .unwrap();
        let system = solver(&thermal, &conduction, &law, &flow, true);
        let step = BackwardEulerStep::new(&thermal, 0.7, || Ok(())).unwrap();
        // Unlike their initial phases, the reference cell is now pure liquid
        // and cell 1 lies on the plateau. A frozen initial tangent is wrong.
        let energy = [190.0, -360.0, -45.0];
        let mut report = empty_report();
        let mut workspace = ImplicitConductionWorkspace::new(3).unwrap();
        let evaluate = |energy: &[f64], report: &mut ImplicitConductionReport| {
            system.evaluate(
                &step,
                energy,
                &BTreeMap::new(),
                &[],
                false,
                controls(1.0),
                report,
                &mut || Ok(()),
            )
        };
        let trial = evaluate(&energy, &mut report).unwrap();
        system
            .jacobian(&thermal, trial.payload(), 0.7, &mut workspace, &mut || {
                Ok(())
            })
            .unwrap();
        let coefficient = |row, col| {
            if row == col {
                workspace.diagonal[row]
            } else {
                workspace
                    .graph_entries
                    .as_ref()
                    .unwrap()
                    .iter()
                    .filter(|&&(r, c, _)| r == row && c == col)
                    .map(|entry| entry.2)
                    .sum::<f64>()
            }
        };
        near(coefficient(1, 1), 1.0, 0.0);
        for col in 0..3 {
            let mut plus = energy;
            let mut minus = energy;
            plus[col] += 1e-3;
            minus[col] -= 1e-3;
            let a = evaluate(&plus, &mut report).unwrap();
            let b = evaluate(&minus, &mut report).unwrap();
            for row in 0..3 {
                near(
                    coefficient(row, col),
                    (a.residual_j()[row] - b.residual_j()[row]) / 0.002,
                    1e-7,
                );
            }
        }
    }
}

#[test]
fn implicit_liquid_sdirk_preserves_the_single_phase_limit_and_second_order() {
    let law = phase(0.0, false);
    let (thermal, conduction) = setup(&law, &[550.0], None, vec![]);
    let original = open_flow(1, 2.0, 600.0);
    let flow = original.clone().with_liquid_carrier(carrier(&law)).unwrap();
    let system = solver(&thermal, &conduction, &law, &flow, false);
    let homogeneous = solver(&thermal, &conduction, &law, &original, false);
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    let exact = 600.0 - 50.0 * (-2.0_f64).exp();
    let mut errors = Vec::new();
    for dt in [0.2, 0.1, 0.05] {
        let (next, report) = system
            .advance_sdirk2_to(&thermal, 1.0, controls(dt), &mut workspace, || Ok(()))
            .unwrap();
        let (bulk, bulk_report) = homogeneous
            .advance_sdirk2_to(&thermal, 1.0, controls(dt), &mut workspace, || Ok(()))
            .unwrap();
        assert_eq!(next.cells(), bulk.cells());
        assert_eq!(report, bulk_report);
        let h = next.cells()[&CellId(0)].specific_h_j_kg;
        errors.push((h - exact).abs());
        near(report.energy.external_heat_j, h - 550.0, 1e-9);
    }
    assert!(
        errors[1] < 0.3 * errors[0] && errors[2] < 0.3 * errors[1],
        "{errors:?}"
    );
}

#[test]
fn implicit_liquid_refusals_and_cancellation_preserve_original_histories() {
    let law = phase(0.0, false);
    let (thermal, conduction) = setup(&law, &[400.0], None, vec![]);
    let flow = open_flow(1, 2.0, 600.0)
        .with_liquid_carrier(carrier(&law))
        .unwrap();
    for h in [250.0, 300.0] {
        let mut solid = thermal.clone();
        solid.cells.get_mut(&CellId(0)).unwrap().specific_h_j_kg = h;
        assert!(
            ImplicitConductionChain::new(&solid, &conduction, &[(MaterialId(0), &law)], || Ok(()))
                .unwrap()
                .with_implicit_advection(&solid, &flow, || Ok(()))
                .is_err()
        );
    }
    let system = solver(&thermal, &conduction, &law, &flow, false);
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    let mut calls = 0;
    let (expected, expected_report) = system
        .advance_to(&thermal, 1.0, controls(1.0), &mut workspace, || {
            calls += 1;
            Ok(())
        })
        .unwrap();
    let mut cancelled_calls = 0;
    assert!(matches!(
        system.advance_to(&thermal, 1.0, controls(1.0), &mut workspace, || {
            cancelled_calls += 1;
            if cancelled_calls + 1 == calls {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        },),
        Err(ImplicitConductionError::Thermal(DepositionError::Cancelled))
    ));
    let (repeated, report) = system
        .advance_to(&thermal, 1.0, controls(1.0), &mut workspace, || Ok(()))
        .unwrap();
    assert_eq!(repeated.cells(), expected.cells());
    assert_eq!(report, expected_report);
    assert_eq!(thermal.cells()[&CellId(0)].specific_h_j_kg, 400.0);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(flow.audit().substeps, 0);

    // The original residual returns a domain exit for a trial with no liquid,
    // permitting backtracking without ever accepting flow through solid.
    let step = BackwardEulerStep::new(&thermal, 1.0, || Ok(())).unwrap();
    let mut report = empty_report();
    assert!(matches!(
        system.evaluate(
            &step,
            &[-101.0],
            &BTreeMap::new(),
            &[],
            false,
            controls(1.0),
            &mut report,
            &mut || Ok(())
        ),
        Err(ImplicitConductionError::Thermal(
            DepositionError::MaterialDomain
        ))
    ));
}
