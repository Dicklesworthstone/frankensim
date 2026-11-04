use super::*;
use crate::enthalpy::advection::{Advected, MassBoundary, MassFlow};

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
