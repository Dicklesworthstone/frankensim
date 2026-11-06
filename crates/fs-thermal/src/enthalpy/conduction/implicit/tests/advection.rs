use super::*;
use crate::enthalpy::advection::reconstruction::{AxialCellStencil, LimitedEnthalpyReconstruction};
use crate::enthalpy::advection::{Advected, MassBoundary, MassFlow};

#[path = "advection/liquid.rs"]
mod liquid;
#[path = "advection/species.rs"]
mod species;

fn circulating<O>(inner: O, rate: f64) -> Advected<O> {
    Advected::new(
        inner,
        vec![
            MassFlow {
                id: 10,
                a: CellId(0),
                b: CellId(1),
                kg_s: rate,
            },
            MassFlow {
                id: 11,
                a: CellId(1),
                b: CellId(0),
                kg_s: rate,
            },
        ],
        2,
    )
    .unwrap()
}

struct ShiftedLaw {
    base: Law,
    datum: f64,
}
impl ThermalLaw for ShiftedLaw {
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        Some([7; 32])
    }
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        self.base.temperature_k(h - self.datum)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
    fn specific_entropy_change_j_kg_k(
        &self,
        a: f64,
        b: f64,
    ) -> Result<Option<f64>, DepositionError> {
        self.base
            .specific_entropy_change_j_kg_k(a - self.datum, b - self.datum)
    }
}
impl ConductivityLaw for ShiftedLaw {
    fn conductivity_w_m_k(&self, h: f64) -> Result<f64, DepositionError> {
        self.base.conductivity_w_m_k(h - self.datum)
    }
}
impl ImplicitConductionLaw for ShiftedLaw {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self
    }
    fn conductivity_law(&self) -> &dyn ConductivityLaw {
        self
    }
    fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        [100.0 + self.datum, 1000.0 + self.datum]
    }
    fn state_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<ImplicitConductionState, DepositionError> {
        self.base.state_at_specific_enthalpy(h - self.datum)
    }
}

#[test]
fn observed_advection_stiff_chain_and_graph_preserve_energy_and_datum() {
    for datum in [0.0, 1e9] {
        let law = ShiftedLaw {
            base: Law {
                nonlinear_k: false,
                plateau: false,
            },
            datum,
        };
        let (thermal, conduction) = setup(&law, &[400.0 + datum, 200.0 + datum], None, vec![]);
        let chain = ImplicitConductionChain::new(
            &thermal,
            &conduction,
            &[(MaterialId(0), &law)],
            || Ok(()),
        )
        .unwrap();
        let graph = ImplicitConductionGraph::new(
            &thermal,
            &conduction,
            &[(MaterialId(0), &law)],
            ImplicitConductionLinearControls::default(),
            || Ok(()),
        )
        .unwrap();
        let flow = circulating((), 0.02).with_entropy_audit();
        let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
        for graph_path in [false, true] {
            let (next, accepted, report) = if graph_path {
                graph.advance_observed_to(
                    &thermal,
                    &flow,
                    8.0,
                    controls(10.0),
                    &mut workspace,
                    || Ok(()),
                )
            } else {
                chain.advance_observed_to(
                    &thermal,
                    &flow,
                    8.0,
                    controls(10.0),
                    &mut workspace,
                    || Ok(()),
                )
            }
            .unwrap();
            // Unit masses/Cp/G. Explicit reciprocal flow changes the initial
            // contrast by (1-2*r*dt/m); implicit conduction divides by
            // (1+2*G*dt/(m*Cp)). The dt=8 interval is ten diffusion CFL limits.
            let contrast = 200.0 * (1.0 - 2.0 * 0.02 * 8.0) / (1.0 + 2.0 * 8.0);
            near(
                next.temperature_k(CellId(0)).unwrap(),
                300.0 + contrast * 0.5,
                1e-10,
            );
            near(
                next.temperature_k(CellId(1)).unwrap(),
                300.0 - contrast * 0.5,
                1e-10,
            );
            near(report.energy.external_heat_j, 0.0, 0.0);
            near(report.energy.stored_energy_change_j, 0.0, 1e-10);
            near(report.energy.residual_j, 0.0, 1e-10);
            assert_eq!(report.energy.substeps, 1);
            assert!(report.entropy.has_unclassified_sources);
            near(accepted.audit().transported_mass_kg, 0.32, 1e-15);
            near(accepted.audit().external_enthalpy_j, 0.0, 0.0);
            near(accepted.audit().maximum_courant, 0.16, 1e-15);
            near(
                accepted
                    .entropy_audit()
                    .unwrap()
                    .advective_balance_residual_j_k,
                0.0,
                1e-12,
            );
        }
        assert_eq!(flow.audit().substeps, 0);
        assert_eq!(thermal.time_s(), 0.0);
    }
}

struct PhaseContact<'a>(&'a Law);
impl contact::ImplicitContactLaw for PhaseContact<'_> {
    fn caloric_laws(&self) -> [&dyn ThermalLaw; 2] {
        [self.0, self.0]
    }
    fn enthalpy_bounds_j_kg(&self) -> [[f64; 2]; 2] {
        [[100.0, 1000.0]; 2]
    }
    fn state_at_specific_enthalpies(
        &self,
        h: [f64; 2],
    ) -> Result<contact::ImplicitContactState, DepositionError> {
        Ok(contact::ImplicitContactState {
            resistance_m2_k_w: 0.1 + 0.0005 * (h[0] - 100.0) + 0.0002 * (h[1] - 100.0),
            resistance_enthalpy_derivatives_m2_kg_k_w_j: [0.0005, 0.0002],
        })
    }
}

#[test]
fn observed_advection_open_channel_keeps_latent_heat_and_nonlinear_contact() {
    let law = Law {
        nonlinear_k: true,
        plateau: true,
    };
    let (thermal, conduction) = setup(&law, &[400.0, 550.0], None, vec![]);
    let contact = PhaseContact(&law);
    let chain =
        ImplicitConductionChain::new(&thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap()
            .with_contacts(&thermal, &[(0, &contact)], || Ok(()))
            .unwrap();
    let flow = Advected::new(
        (),
        vec![MassFlow {
            id: 10,
            a: CellId(0),
            b: CellId(1),
            kg_s: 0.2,
        }],
        3,
    )
    .unwrap()
    .with_boundaries(
        vec![
            MassBoundary::Inlet {
                id: 11,
                cell: CellId(0),
                kg_s: 0.2,
                material_id: MaterialId(0),
                specific_h_j_kg: 600.0,
            },
            MassBoundary::Outlet {
                id: 12,
                cell: CellId(1),
                kg_s: 0.2,
            },
        ],
        3,
    )
    .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
    let (next, accepted, report) = chain
        .advance_observed_to(&thermal, &flow, 1.0, controls(2.0), &mut workspace, || {
            Ok(())
        })
        .unwrap();

    // Independent endpoint balance: the stream adds 40 J to A and removes
    // 30 J from B. Let q be conductive heat B -> A. A stays on its latent
    // plateau; B has T=320-q, k=T/100, and the original two-sided resistance.
    let (mut lo, mut hi) = (0.0, 20.0);
    for _ in 0..80 {
        let q = (lo + hi) * 0.5;
        let a = 440.0 + q;
        let b = 520.0 - q;
        let resistance = 0.1 + 0.0005 * (a - 100.0) + 0.0002 * (b - 100.0);
        let g = 1.0 / (0.5 / 3.0 + resistance + 0.5 / ((320.0 - q) / 100.0));
        if q < g * (20.0 - q) {
            lo = q;
        } else {
            hi = q;
        }
    }
    let q = (lo + hi) * 0.5;
    near(next.cells()[&CellId(0)].specific_h_j_kg, 440.0 + q, 1e-9);
    near(next.cells()[&CellId(1)].specific_h_j_kg, 520.0 - q, 1e-9);
    near(next.temperature_k(CellId(0)).unwrap(), 300.0, 0.0);
    near(next.temperature_k(CellId(1)).unwrap(), 320.0 - q, 1e-9);
    near(report.energy.external_heat_j, 10.0, 1e-12);
    near(report.energy.stored_energy_change_j, 10.0, 1e-9);
    near(accepted.audit().external_enthalpy_j, 10.0, 1e-12);
    near(accepted.audit().inlet_mass_kg, 0.2, 1e-15);
    near(accepted.audit().outlet_mass_kg, 0.2, 1e-15);
    assert_eq!(report.energy.substeps, 1);
    assert!(report.entropy.conduction_generation_j_k > 0.0);
}

#[test]
fn observed_advection_retries_prepare_flow_once_and_refusals_keep_history() {
    let law = FallingConductivity;
    let (thermal, conduction) = setup(&law, &[400.0, 420.0], Some(250.0), vec![]);
    let chain =
        ImplicitConductionChain::new(&thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap();
    let flow = circulating(Forcing::default(), 0.1);
    let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
    let (next, accepted, report) = chain
        .advance_observed_to(&thermal, &flow, 4.0, controls(8.0), &mut workspace, || {
            Ok(())
        })
        .unwrap();
    // The original flow limit reduces 8 -> 4 s; the starting latent-branch
    // reservoir Jacobian is singular at dt=4. Retrying it at dt=2 must reuse
    // the pending flow, then prepare exactly one new epoch for [2,4].
    assert_eq!(report.rejected_intervals, 1);
    assert_eq!(report.energy.substeps, 2);
    assert_eq!(accepted.inner().transfer_evaluations, 2);
    assert_eq!(accepted.inner().evaluations, 3);
    assert_eq!(accepted.inner().begins, 1);
    assert_eq!(
        accepted.inner().accepted,
        [(0.0, 2.0, 0.0), (2.0, 4.0, 0.0)]
    );
    assert_eq!(accepted.audit().substeps, 2);
    near(accepted.audit().transported_mass_kg, 0.8, 1e-14);
    near(accepted.audit().maximum_courant, 0.2, 1e-14);
    near(
        report.energy.stored_energy_change_j,
        report.energy.external_heat_j,
        1e-9,
    );
    assert_eq!(next.time_s(), 4.0);

    let refusing = circulating(
        Forcing {
            refuse_at_s: 4.0,
            ..Forcing::default()
        },
        0.1,
    );
    assert!(matches!(
        chain.advance_observed_to(
            &thermal,
            &refusing,
            4.0,
            controls(8.0),
            &mut workspace,
            || Ok(()),
        ),
        Err(ImplicitConductionObservedError::Observer(
            DepositionError::MaterialDomain
        ))
    ));
    let unbalanced = Advected::new(
        Forcing::default(),
        vec![MassFlow {
            id: 10,
            a: CellId(0),
            b: CellId(1),
            kg_s: 0.1,
        }],
        1,
    )
    .unwrap();
    assert!(matches!(
        chain.advance_observed_to(
            &thermal,
            &unbalanced,
            4.0,
            controls(8.0),
            &mut workspace,
            || Ok(()),
        ),
        Err(ImplicitConductionObservedError::Observer(
            DepositionError::InvalidInput(_)
        ))
    ));
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.cells()[&CellId(0)].specific_h_j_kg, 400.0);
    assert_eq!(thermal.cells()[&CellId(1)].specific_h_j_kg, 420.0);
    assert_eq!(refusing.audit().substeps, 0);
    assert_eq!(refusing.inner().time_s, 0.0);
    assert_eq!(refusing.inner().transfer_evaluations, 0);
    assert_eq!(refusing.inner().begins, 0);
}

#[test]
fn implicit_upwind_stiff_pair_uses_endpoint_enthalpy_and_retains_the_datum() {
    for datum in [0.0, 1e9] {
        let law = ShiftedLaw {
            base: Law {
                nonlinear_k: false,
                plateau: false,
            },
            datum,
        };
        let (thermal, conduction) = setup(&law, &[400.0 + datum, 200.0 + datum], None, vec![]);
        let flow = circulating((), 30.5);
        let chain = ImplicitConductionChain::new(
            &thermal,
            &conduction,
            &[(MaterialId(0), &law)],
            || Ok(()),
        )
        .unwrap()
        .with_implicit_advection(&thermal, &flow, || Ok(()))
        .unwrap();
        let graph = ImplicitConductionGraph::new(
            &thermal,
            &conduction,
            &[(MaterialId(0), &law)],
            ImplicitConductionLinearControls::default(),
            || Ok(()),
        )
        .unwrap()
        .with_implicit_advection(&thermal, &flow, || Ok(()))
        .unwrap();
        let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
        for graph_path in [false, true] {
            let (next, report) = if graph_path {
                graph.advance_to(&thermal, 1.0, controls(1.0), &mut workspace, || Ok(()))
            } else {
                chain.advance_to(&thermal, 1.0, controls(1.0), &mut workspace, || Ok(()))
            }
            .unwrap();
            // m=Cp=G=1, reciprocal mdot=30.5: the exact BE contrast is
            // 200/[1+2*(G/Cp+mdot)*dt/m] = 3.125. Courant=30.5 is well
            // above the explicit bound, without freezing the advected state.
            near(next.temperature_k(CellId(0)).unwrap(), 301.5625, 1e-9);
            near(next.temperature_k(CellId(1)).unwrap(), 298.4375, 1e-9);
            near(report.energy.external_heat_j, 0.0, 0.0);
            near(report.energy.stored_energy_change_j, 0.0, 1e-9);
            assert_eq!(report.energy.substeps, 1);
            assert!(report.entropy.has_unclassified_sources);
        }
        assert_eq!(flow.audit().substeps, 0);
        assert_eq!(thermal.time_s(), 0.0);
    }
}

fn open_flow(cells: usize, rate: f64, incoming_h: f64) -> Advected<()> {
    Advected::new(
        (),
        (0..cells - 1)
            .map(|i| MassFlow {
                id: i as u64,
                a: CellId(i as u64),
                b: CellId(i as u64 + 1),
                kg_s: rate,
            })
            .collect(),
        cells,
    )
    .unwrap()
    .with_boundaries(
        vec![
            MassBoundary::Inlet {
                id: cells as u64,
                cell: CellId(if rate > 0.0 { 0 } else { cells as u64 - 1 }),
                kg_s: rate.abs(),
                material_id: MaterialId(0),
                specific_h_j_kg: incoming_h,
            },
            MassBoundary::Outlet {
                id: cells as u64 + 1,
                cell: CellId(if rate > 0.0 { cells as u64 - 1 } else { 0 }),
                kg_s: rate.abs(),
            },
        ],
        cells + 1,
    )
    .unwrap()
}

#[test]
fn implicit_open_stream_moves_latent_heat_at_identical_temperatures() {
    let law = Law {
        nonlinear_k: true,
        plateau: true,
    };
    let (mut thermal, conduction) = setup(&law, &[350.0, 400.0], None, vec![]);
    thermal.cells.get_mut(&CellId(1)).unwrap().mass_kg = 2.0;
    let flow = open_flow(2, 4.0, 450.0);
    let chain =
        ImplicitConductionChain::new(&thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap()
            .with_implicit_advection(&thermal, &flow, || Ok(()))
            .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
    let (next, report) = chain
        .advance_to(&thermal, 2.0, controls(2.0), &mut workspace, || Ok(()))
        .unwrap();
    let first = (350.0 + 8.0 * 450.0) / 9.0;
    let second = (2.0 * 400.0 + 8.0 * first) / 10.0;
    near(next.cells()[&CellId(0)].specific_h_j_kg, first, 1e-10);
    near(next.cells()[&CellId(1)].specific_h_j_kg, second, 1e-10);
    for id in [CellId(0), CellId(1)] {
        near(next.temperature_k(id).unwrap(), 300.0, 0.0);
    }
    let external = 8.0 * (450.0 - second);
    near(report.energy.external_heat_j, external, 1e-9);
    near(
        report.energy.stored_energy_change_j,
        first - 350.0 + 2.0 * (second - 400.0),
        1e-9,
    );
    near(report.energy.residual_j, 0.0, 1e-9);
    assert_eq!(report.energy.substeps, 1);
    assert_eq!(next.cells()[&CellId(1)].mass_kg, 2.0);
}

#[test]
fn implicit_directed_cycle_has_the_original_nonlinear_jacobian_and_balance() {
    let law = Law {
        nonlinear_k: true,
        plateau: true,
    };
    let (mut thermal, conduction) = setup(&law, &[400.0, 350.0, 250.0], None, vec![]);
    for (i, mass) in [1.0, 2.0, 3.0].into_iter().enumerate() {
        thermal.cells.get_mut(&CellId(i as u64)).unwrap().mass_kg = mass;
    }
    let flow = Advected::new(
        (),
        vec![
            MassFlow {
                id: 1,
                a: CellId(0),
                b: CellId(2),
                kg_s: 2.5,
            },
            MassFlow {
                id: 2,
                a: CellId(1),
                b: CellId(2),
                kg_s: -2.5,
            },
            MassFlow {
                id: 3,
                a: CellId(1),
                b: CellId(0),
                kg_s: 2.5,
            },
        ],
        3,
    )
    .unwrap();
    assert!(matches!(
        ImplicitConductionChain::new(&thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()),)
            .unwrap()
            .with_implicit_advection(&thermal, &flow, || Ok(())),
        Err(ImplicitConductionError::Invalid(
            "nontridiagonal enthalpy transport requires the graph solver"
        ))
    ));
    let graph = ImplicitConductionChain::prepare(
        &thermal,
        &conduction,
        &[(MaterialId(0), &law)],
        LinearTopology::Graph(ImplicitConductionLinearControls::default()),
        || Ok(()),
    )
    .unwrap()
    .with_implicit_advection(&thermal, &flow, || Ok(()))
    .unwrap();
    let step = BackwardEulerStep::new(&thermal, 0.7, || Ok(())).unwrap();
    let energy = [3.0, -4.0, 2.0];
    let mut report = empty_report();
    let mut workspace = ImplicitConductionWorkspace::new(3).unwrap();
    let trial = graph
        .evaluate(
            &step,
            &energy,
            &BTreeMap::new(),
            &[],
            false,
            controls(1.0),
            &mut report,
            &mut || Ok(()),
        )
        .unwrap();
    graph
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
                .sum()
        }
    };
    assert!((coefficient(0, 1) - coefficient(1, 0)).abs() > 0.1);
    for col in 0..3 {
        let mut plus = energy;
        let mut minus = energy;
        plus[col] += 1e-3;
        minus[col] -= 1e-3;
        let a = graph
            .evaluate(
                &step,
                &plus,
                &BTreeMap::new(),
                &[],
                false,
                controls(1.0),
                &mut report,
                &mut || Ok(()),
            )
            .unwrap();
        let b = graph
            .evaluate(
                &step,
                &minus,
                &BTreeMap::new(),
                &[],
                false,
                controls(1.0),
                &mut report,
                &mut || Ok(()),
            )
            .unwrap();
        for row in 0..3 {
            near(
                coefficient(row, col),
                (a.residual_j()[row] - b.residual_j()[row]) / 0.002,
                1e-8,
            );
        }
    }
    let (next, report) = graph
        .advance_to(&thermal, 0.7, controls(0.7), &mut workspace, || Ok(()))
        .unwrap();
    let h: Vec<_> = (0..3)
        .map(|i| next.cells()[&CellId(i)].specific_h_j_kg)
        .collect();
    let t: Vec<_> = h.iter().map(|&h| law.temperature_k(h).unwrap()).collect();
    let k: Vec<_> = h
        .iter()
        .map(|&h| law.conductivity_w_m_k(h).unwrap())
        .collect();
    let q01 = (t[0] - t[1]) / (0.5 / k[0] + 0.5 / k[1]);
    let q12 = (t[1] - t[2]) / (0.5 / k[1] + 0.5 / k[2]);
    let powers = [
        -q01 + 2.5 * (h[1] - h[0]),
        q01 - q12 + 2.5 * (h[2] - h[1]),
        q12 + 2.5 * (h[0] - h[2]),
    ];
    for (i, old) in [400.0, 350.0, 250.0].into_iter().enumerate() {
        near((i + 1) as f64 * (h[i] - old), 0.7 * powers[i], 1e-9);
    }
    near(report.energy.external_heat_j, 0.0, 0.0);
    near(report.energy.stored_energy_change_j, 0.0, 1e-9);
}

#[test]
fn implicit_upwind_time_refinement_and_adaptation_follow_the_real_stream_balance() {
    let law = Law {
        nonlinear_k: false,
        plateau: false,
    };
    let (thermal, conduction) = setup(&law, &[400.0], None, vec![]);
    let flow = open_flow(1, 2.0, 300.0);
    let chain =
        ImplicitConductionChain::new(&thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap()
            .with_implicit_advection(&thermal, &flow, || Ok(()))
            .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    let exact = 300.0 + 100.0 * (-2.0_f64).exp();
    let mut errors = Vec::new();
    for dt in [0.1, 0.05, 0.025] {
        let (next, report) = chain
            .advance_to(&thermal, 1.0, controls(dt), &mut workspace, || Ok(()))
            .unwrap();
        let h = next.cells()[&CellId(0)].specific_h_j_kg;
        errors.push((h - exact).abs());
        near(report.energy.external_heat_j, h - 400.0, 1e-9);
    }
    assert!(
        errors[1] < 0.6 * errors[0] && errors[2] < 0.6 * errors[1],
        "{errors:?}"
    );
    let policy = |tolerance| adaptive::AdaptiveImplicitConductionControls {
        step: controls(1.0),
        initial_dt_s: 1.0,
        maximum_trials: 4096,
        enthalpy_absolute_tolerance_j_kg: tolerance,
        temperature_absolute_tolerance_k: tolerance,
        phase_fraction_absolute_tolerance: None,
        relative_tolerance: 0.0,
    };
    let (loose, a) = chain
        .advance_adaptive_to(&thermal, 1.0, policy(1.0), &mut workspace, || Ok(()))
        .unwrap();
    let (tight, b) = chain
        .advance_adaptive_to(&thermal, 1.0, policy(0.001), &mut workspace, || Ok(()))
        .unwrap();
    assert!(
        (tight.temperature_k(CellId(0)).unwrap() - exact).abs()
            < 0.1 * (loose.temperature_k(CellId(0)).unwrap() - exact).abs()
    );
    assert!(b.accepted_comparisons > a.accepted_comparisons);
    assert!(b.accuracy_rejections > 0 && b.maximum_accepted_normalized_error <= 1.0);
    near(
        b.solve.energy.external_heat_j,
        tight.cells()[&CellId(0)].specific_h_j_kg - 400.0,
        1e-9,
    );
}

#[test]
fn implicit_upwind_preserves_failure_rollback_and_explicit_source_ownership() {
    let law = Law {
        nonlinear_k: false,
        plateau: false,
    };
    let (thermal, conduction) = setup(&law, &[400.0], None, vec![]);
    let flow = open_flow(1, 2.0, 300.0);
    let bind = |flow: &Advected<()>| {
        ImplicitConductionChain::new(&thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap()
            .with_implicit_advection(&thermal, flow, || Ok(()))
    };
    assert!(bind(&flow.clone().with_entropy_audit()).is_err());
    assert!(bind(&open_flow(1, 2.0, 2000.0)).is_err());
    let chain = bind(&flow).unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    let source = Forcing {
        constant_w: 1.0,
        ..Forcing::default()
    };
    let mut calls = 0;
    let (next, accepted, report) = chain
        .advance_observed_to(
            &thermal,
            &source,
            2.0,
            controls(2.0),
            &mut workspace,
            || {
                calls += 1;
                Ok(())
            },
        )
        .unwrap();
    near(
        next.temperature_k(CellId(0)).unwrap(),
        (400.0 + 2.0 * (600.0 + 1.0)) / 5.0,
        1e-10,
    );
    near(accepted.energy_j, 2.0, 1e-12);
    near(
        report.energy.external_heat_j,
        next.temperature_k(CellId(0)).unwrap() - 400.0,
        1e-9,
    );
    let mut cancelled_calls = 0;
    assert!(matches!(
        chain.advance_observed_to(
            &thermal,
            &source,
            2.0,
            controls(2.0),
            &mut workspace,
            || {
                cancelled_calls += 1;
                if cancelled_calls + 1 == calls {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            }
        ),
        Err(ImplicitConductionObservedError::Implicit(
            ImplicitConductionError::Thermal(DepositionError::Cancelled)
        ))
    ));
    let mut exhausted = controls(2.0);
    exhausted.maximum_residual_evaluations = 1;
    assert!(matches!(
        chain.advance_to(&thermal, 2.0, exhausted, &mut workspace, || Ok(())),
        Err(ImplicitConductionError::Thermal(DepositionError::Budget(_)))
    ));
    let (retried, retried_source, retried_report) = chain
        .advance_observed_to(
            &thermal,
            &source,
            2.0,
            controls(2.0),
            &mut workspace,
            || Ok(()),
        )
        .unwrap();
    assert_eq!(retried.cells(), next.cells());
    assert_eq!(retried_source.energy_j, accepted.energy_j);
    assert_eq!(retried_report, report);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.cells()[&CellId(0)].specific_h_j_kg, 400.0);
    assert_eq!(flow.audit().substeps, 0);
    assert_eq!(source.energy_j, 0.0);
}

fn axial_plan(thermal: &DepositionSimulation<'_>) -> LimitedEnthalpyReconstruction {
    let ids: Vec<_> = thermal.cells().keys().copied().collect();
    let stencils: Vec<_> = ids
        .windows(3)
        .map(|ids| AxialCellStencil {
            previous: ids[0],
            cell: ids[1],
            next: ids[2],
        })
        .collect();
    LimitedEnthalpyReconstruction::try_new(thermal, &stencils, ids.len(), || Ok(())).unwrap()
}

#[test]
fn implicit_limited_advection_differentiates_the_actual_endpoint_and_flow_direction() {
    let law = Law {
        nonlinear_k: true,
        plateau: true,
    };
    for rate in [7.0, -7.0] {
        let (mut thermal, conduction) =
            setup(&law, &[200.0, 250.0, 400.0, 500.0, 700.0], None, vec![]);
        for (i, mass) in [1.0, 2.0, 1.5, 0.75, 3.0].into_iter().enumerate() {
            thermal.cells.get_mut(&CellId(i as u64)).unwrap().mass_kg = mass;
        }
        let flow = open_flow(5, rate, 450.0);
        let plan = axial_plan(&thermal);
        let graph = ImplicitConductionChain::prepare(
            &thermal,
            &conduction,
            &[(MaterialId(0), &law)],
            LinearTopology::Graph(ImplicitConductionLinearControls::default()),
            || Ok(()),
        )
        .unwrap()
        .with_implicit_advection_reconstruction(&thermal, &flow, &plan, || Ok(()))
        .unwrap();
        let step = BackwardEulerStep::new(&thermal, 0.7, || Ok(())).unwrap();
        // These candidate increments change the active minmod branches at
        // cells 1 and 2. A tangent frozen at the starting state is incorrect.
        let energy = [-1.0, 102.0, -80.0, 7.0, -3.0];
        let mut report = empty_report();
        let mut workspace = ImplicitConductionWorkspace::new(5).unwrap();
        let trial = graph
            .evaluate(
                &step,
                &energy,
                &BTreeMap::new(),
                &[],
                false,
                controls(1.0),
                &mut report,
                &mut || Ok(()),
            )
            .unwrap();
        graph
            .jacobian(&thermal, trial.payload(), 0.7, &mut workspace, &mut || {
                Ok(())
            })
            .unwrap();
        assert!(
            workspace
                .graph_entries
                .as_ref()
                .unwrap()
                .iter()
                .any(|&(row, col, value)| row.abs_diff(col) > 1 && value != 0.0)
        );
        for col in 0..5 {
            let mut plus = energy;
            let mut minus = energy;
            plus[col] += 1e-3;
            minus[col] -= 1e-3;
            let a = graph
                .evaluate(
                    &step,
                    &plus,
                    &BTreeMap::new(),
                    &[],
                    false,
                    controls(1.0),
                    &mut report,
                    &mut || Ok(()),
                )
                .unwrap();
            let b = graph
                .evaluate(
                    &step,
                    &minus,
                    &BTreeMap::new(),
                    &[],
                    false,
                    controls(1.0),
                    &mut report,
                    &mut || Ok(()),
                )
                .unwrap();
            for row in 0..5 {
                let analytic = if row == col {
                    workspace.diagonal[row]
                } else {
                    workspace
                        .graph_entries
                        .as_ref()
                        .unwrap()
                        .iter()
                        .filter(|&&(r, c, _)| r == row && c == col)
                        .map(|entry| entry.2)
                        .sum()
                };
                near(
                    analytic,
                    (a.residual_j()[row] - b.residual_j()[row]) / 0.002,
                    1e-8,
                );
            }
        }
    }
}

#[test]
fn implicit_limited_latent_channel_closes_the_large_step_affine_balance() {
    for datum in [0.0, 1e9] {
        let law = ShiftedLaw {
            base: Law {
                nonlinear_k: true,
                plateau: true,
            },
            datum,
        };
        for rate in [2.0, -2.0] {
            let mut initial = [350.0, 420.0, 410.0, 420.0, 410.0];
            let mut expected = [350.0, 360.0, 370.0, 380.0, 390.0];
            if rate < 0.0 {
                initial.reverse();
                expected.reverse();
            }
            let initial = initial.map(|h| h + datum);
            let (thermal, conduction) = setup(&law, &initial, None, vec![]);
            let flow = open_flow(5, rate, 350.0 + datum);
            let plan = axial_plan(&thermal);
            let graph = ImplicitConductionGraph::new(
                &thermal,
                &conduction,
                &[(MaterialId(0), &law)],
                ImplicitConductionLinearControls::default(),
                || Ok(()),
            )
            .unwrap()
            .with_implicit_advection_reconstruction(&thermal, &flow, &plan, || Ok(()))
            .unwrap();
            let mut workspace = ImplicitConductionWorkspace::new(5).unwrap();
            let mut policy = controls(2.0);
            policy.cell_tolerance_j = 1e-6;
            policy.total_tolerance_j = 5e-6;
            let mut calls = 0;
            let (next, report) = graph
                .advance_to(&thermal, 2.0, policy, &mut workspace, || {
                    calls += 1;
                    Ok(())
                })
                .unwrap();
            // At the affine endpoint, the interior outlet traces are
            // 365, 375 and 385 J/kg. The boundary donors remain 350 and 390.
            // With mdot*dt=4 kg, the five independently balanced cell stores
            // are exactly the authored initial values above. T stays 300 K:
            // this transports latent enthalpy, with zero conductive power.
            for (i, h) in expected.into_iter().enumerate() {
                near(
                    next.cells()[&CellId(i as u64)].specific_h_j_kg - datum,
                    h,
                    2e-6,
                );
                near(next.temperature_k(CellId(i as u64)).unwrap(), 300.0, 0.0);
            }
            near(report.energy.external_heat_j, -160.0, 5e-6);
            near(report.energy.stored_energy_change_j, -160.0, 5e-6);
            near(report.energy.residual_j, 0.0, 5e-6);
            assert_eq!(report.energy.substeps, 1); // Courant = 4.
            assert_eq!(flow.audit().substeps, 0);
            if datum == 0.0 && rate > 0.0 {
                let mut cancel_calls = 0;
                assert!(matches!(
                    graph.advance_to(&thermal, 2.0, policy, &mut workspace, || {
                        cancel_calls += 1;
                        if cancel_calls + 1 == calls {
                            Err(DepositionError::Cancelled)
                        } else {
                            Ok(())
                        }
                    }),
                    Err(ImplicitConductionError::Thermal(DepositionError::Cancelled))
                ));
                assert_eq!(thermal.time_s(), 0.0);
                for (i, h) in initial.into_iter().enumerate() {
                    assert_eq!(thermal.cells()[&CellId(i as u64)].specific_h_j_kg, h);
                }
            }
        }
    }
}

#[test]
fn implicit_limited_loop_refines_a_transported_wave_with_less_spatial_diffusion() {
    use std::f64::consts::{PI, TAU};
    let law = Law {
        nonlinear_k: false,
        plateau: false,
    };
    let mut limited_errors = Vec::new();
    for n in [16, 32] {
        let nf = n as f64;
        // Unit-perimeter polygon, rho*A=1: mdot=1 gives axial speed 1.
        let radius = 1.0 / (2.0 * nf * (PI / nf).sin());
        let vertices: Vec<_> = (0..n)
            .map(|i| {
                let angle = TAU * i as f64 / nf;
                [radius * angle.cos(), radius * angle.sin(), 0.0]
            })
            .collect();
        let amplitude = 20.0 * (PI / nf).sin() / (PI / nf);
        let average =
            |i: usize, time: f64| 300.0 + amplitude * (TAU * ((i as f64 + 0.5) / nf - time)).sin();
        let events = (0..n)
            .map(|i| TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id: CellId(i as u64),
                    material_id: MaterialId(0),
                    body_id: 0,
                    segment_id: i as u64,
                    mass_kg: 1.0 / nf,
                    specific_h_j_kg: average(i, 0.0),
                    geometry: CellGeometry {
                        start_m: vertices[i],
                        end_m: vertices[(i + 1) % n],
                        volume_m3: 1.0 / nf,
                        surface_area_m2: 1.0,
                    },
                    boundary: boundary(None),
                }),
            })
            .collect();
        let thermal = DepositionSimulation::new(
            events,
            &[(MaterialId(0), &law)],
            DepositionBudget {
                max_cells: n,
                max_links: n,
                max_events: n,
                max_substeps: 4096,
            },
        )
        .unwrap()
        .advance_to(
            0.0,
            StepConfig {
                max_dt_s: 1.0,
                energy_tolerance_j: 1e-9,
            },
            || Ok(()),
        )
        .unwrap()
        .0;
        let conduction =
            ConductionObserver::new(&thermal, &[(MaterialId(0), &law)], &[], &[], || Ok(()))
                .unwrap();
        let flow = Advected::new(
            (),
            (0..n)
                .map(|i| MassFlow {
                    id: i as u64,
                    a: CellId(i as u64),
                    b: CellId(((i + 1) % n) as u64),
                    kg_s: 1.0,
                })
                .collect(),
            n,
        )
        .unwrap();
        let stencils: Vec<_> = (0..n)
            .map(|i| AxialCellStencil {
                previous: CellId(((i + n - 1) % n) as u64),
                cell: CellId(i as u64),
                next: CellId(((i + 1) % n) as u64),
            })
            .collect();
        let plan =
            LimitedEnthalpyReconstruction::try_new(&thermal, &stencils, n, || Ok(())).unwrap();
        let bind = || {
            ImplicitConductionGraph::new(
                &thermal,
                &conduction,
                &[(MaterialId(0), &law)],
                ImplicitConductionLinearControls::default(),
                || Ok(()),
            )
            .unwrap()
        };
        let upwind = bind()
            .with_implicit_advection(&thermal, &flow, || Ok(()))
            .unwrap();
        let limited = bind()
            .with_implicit_advection_reconstruction(&thermal, &flow, &plan, || Ok(()))
            .unwrap();
        let mut errors = Vec::new();
        for graph in [&upwind, &limited] {
            let mut workspace = ImplicitConductionWorkspace::new(n).unwrap();
            // dt scales as dx² so first-order BE time error cannot mask
            // the spatial comparison; this is not a large-step accuracy claim.
            let (next, report) = graph
                .advance_to(
                    &thermal,
                    0.125,
                    controls(0.5 / nf.powi(2)),
                    &mut workspace,
                    || Ok(()),
                )
                .unwrap();
            errors.push(
                (0..n)
                    .map(|i| {
                        (next.cells()[&CellId(i as u64)].specific_h_j_kg - average(i, 0.125)).abs()
                    })
                    .sum::<f64>()
                    / nf,
            );
            near(report.energy.external_heat_j, 0.0, 0.0);
            near(report.energy.stored_energy_change_j, 0.0, 1e-8);
        }
        assert!(
            errors[1] < 0.5 * errors[0],
            "n={n}: upwind={}, limited={}",
            errors[0],
            errors[1]
        );
        limited_errors.push(errors[1]);
    }
    assert!(
        limited_errors[1] < 0.5 * limited_errors[0],
        "{limited_errors:?}"
    );
}

#[path = "advection_entropy.rs"]
mod entropy;


#[test]
fn implicit_entropy_keeps_cancelled_local_and_observed_sources_unclassified() {
    use crate::enthalpy::conduction::implicit::source::{
        ImplicitHeatSourceLaw, ImplicitHeatSourceState,
    };

    struct LocalSource<'a>(&'a Law);
    impl ImplicitHeatSourceLaw for LocalSource<'_> {
        fn caloric_law(&self) -> &dyn ThermalLaw { self.0 }
        fn state_at_specific_enthalpy(
            &self,
            h: f64,
            check: &mut dyn FnMut() -> Result<(), DepositionError>,
        ) -> Result<ImplicitHeatSourceState, DepositionError> {
            check()?;
            self.0.temperature_k(h)?;
            Ok(ImplicitHeatSourceState {
                power_w: 20.0,
                power_enthalpy_derivative_w_kg_j: 0.0,
            })
        }
    }
    #[derive(Clone)]
    struct ObservedSink;
    impl DepositionObserver for ObservedSink {
        type Error = DepositionError;
        fn sources(
            &mut self,
            _: &DepositionSimulation<'_>,
            _: &mut dyn FnMut() -> Result<(), DepositionError>,
        ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
            Ok(BTreeMap::from([(CellId(0), -20.0)]))
        }
    }

    let law = Law { nonlinear_k: false, plateau: false };
    let (thermal, conduction) = setup(&law, &[400.0], None, vec![]);
    let source = LocalSource(&law);
    let flow = open_flow(1, 2.0, 300.0);
    let chain = ImplicitConductionChain::new(
        &thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()),
    ).unwrap().with_heat_sources(&thermal, &[(CellId(0), &source)], || Ok(())).unwrap()
        .with_implicit_advection(&thermal, &flow, || Ok(())).unwrap()
        .with_implicit_advection_entropy(&thermal, || Ok(())).unwrap();
    let (next, _, report) = chain.advance_observed_to(
        &thermal, &ObservedSink, 0.5, controls(0.5),
        &mut ImplicitConductionWorkspace::new(1).unwrap(), || Ok(()),
    ).unwrap();
    near(next.temperature_k(CellId(0)).unwrap(), 350.0, 1e-10);
    near(report.energy.external_heat_j, -50.0, 1e-9);
    near(report.entropy.unclassified_source_at_cell_j_k, 0.0, 0.0);
    assert!(report.entropy.has_unclassified_sources);
    assert!(report.entropy.advection.is_some());
    near(report.entropy.thermal_storage_change_j_k, (350.0_f64 / 400.0).ln(), 1e-12);
}
