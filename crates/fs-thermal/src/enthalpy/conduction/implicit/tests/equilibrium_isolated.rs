//! Closed component energy and original native heat-rate acceptance.
use super::super::equilibrium::isolated::{
    CaloricEquilibriumError, IsolatedThermalControls, IsolatedThermalError,
};
use super::*;

fn disconnected<'a>(law: &'a Law) -> (DepositionSimulation<'a>, ConductionObserver<'a>) {
    let mut events = Vec::new();
    for (i, (energy, mass)) in [
        (500.0, 1.0),
        (100.0, 3.0),
        (900.0, 1.0),
        (300.0, 1.0),
        (444.125, 2.0),
    ]
    .into_iter()
    .enumerate()
    {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64),
                material_id: MaterialId(0),
                body_id: 0,
                segment_id: 0,
                mass_kg: mass,
                specific_h_j_kg: energy,
                geometry: CellGeometry {
                    start_m: [i as f64, 0.0, 0.0],
                    end_m: [i as f64 + 1.0, 0.0, 0.0],
                    volume_m3: 1.0,
                    surface_area_m2: 2.0,
                },
                boundary: boundary(None),
            }),
        });
    }
    let mut links = Vec::new();
    for (id, a, b) in [(0, 0, 1), (1, 2, 3)] {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id,
                a: CellId(a),
                b: CellId(b),
                area_m2: 1.0,
                conductance_w_k: 0.0,
            }),
        });
        links.push(ConductiveLink {
            link_id: id,
            half_distances_m: [0.5; 2],
            contact_resistance_m2_k_w: 0.0,
        });
    }
    let mut future = boundary(None);
    future.source_w = 10.0;
    events.push(TimedEvent {
        time_s: 8.0,
        event: Event::Boundary(CellId(0), future),
    });
    let thermal = DepositionSimulation::new(
        events,
        &[(MaterialId(0), law)],
        DepositionBudget {
            max_cells: 10,
            max_links: 10,
            max_events: 20,
            max_substeps: 100,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 0.0,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    let conduction =
        ConductionObserver::new(&thermal, &[(MaterialId(0), law)], &links, &[], || Ok(())).unwrap();
    (thermal, conduction)
}

#[test]
fn isolated_graph_conserves_each_disconnected_mass_inventory_and_keeps_future_events() {
    let law = Law {
        nonlinear_k: false,
        plateau: false,
    };
    let (thermal, conduction) = disconnected(&law);
    let before = thermal.cells().clone();
    let graph = ImplicitConductionGraph::new(
        &thermal,
        &conduction,
        &[(MaterialId(0), &law)],
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap();
    let result = graph
        .solve_isolated_equilibrium(&thermal, IsolatedThermalControls::default(), || Ok(()))
        .unwrap();
    assert_eq!(result.components().len(), 3);
    for (actual, expected) in result
        .field()
        .temperatures_k()
        .iter()
        .zip([200.0, 200.0, 600.0, 600.0, 444.125])
    {
        near(*actual, expected, 1e-9);
    }
    assert_eq!(result.components()[0].cells, [CellId(0), CellId(1)]);
    assert_eq!(result.components()[1].cells, [CellId(2), CellId(3)]);
    assert_eq!(result.components()[2].cells, [CellId(4)]);
    for component in result.components() {
        assert!(component.energy_residual_j.abs() <= 1e-9);
    }
    assert!(result.components()[2].retained_initial_state);
    assert!(
        result
            .field()
            .net_powers_w()
            .iter()
            .all(|p| p.abs() <= 1e-9)
    );
    assert_eq!(result.field().operating_time_s(), 0.0);
    assert_eq!(thermal.cells(), &before);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.cursor, 7);
}

#[test]
fn isolated_chain_retains_original_latent_state_but_refuses_unselected_allocation() {
    let law = Law {
        nonlinear_k: false,
        plateau: true,
    };
    for (energies, retained) in [([350.125, 499.75], true), ([100.0, 700.0], false)] {
        let (thermal, conduction) = setup(&law, &energies, None, Vec::new());
        let chain = ImplicitConductionChain::new(
            &thermal,
            &conduction,
            &[(MaterialId(0), &law)],
            || Ok(()),
        )
        .unwrap();
        let result =
            chain.solve_isolated_equilibrium(&thermal, IsolatedThermalControls::default(), || {
                Ok(())
            });
        if retained {
            let result = result.unwrap();
            assert_eq!(result.field().specific_enthalpies_j_kg(), &energies);
            assert!(result.components()[0].retained_initial_state);
        } else {
            assert!(matches!(
                result,
                Err(IsolatedThermalError::Caloric(
                    CaloricEquilibriumError::UnderdeterminedPhaseAllocation { cells: 2 }
                ))
            ));
        }
        assert_eq!(thermal.cells()[&CellId(0)].specific_h_j_kg, energies[0]);
        assert_eq!(thermal.cells()[&CellId(1)].specific_h_j_kg, energies[1]);
    }
}

#[test]
fn isolated_final_native_watt_gate_cannot_be_bypassed_by_a_loose_temperature_tolerance() {
    let law = Law {
        nonlinear_k: false,
        plateau: false,
    };
    let (thermal, conduction) = setup(&law, &[400.0, 200.0], None, Vec::new());
    let chain =
        ImplicitConductionChain::new(&thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap();
    let mut controls = IsolatedThermalControls::default();
    controls.caloric.temperature_tolerance_k = 300.0;
    assert!(matches!(
        chain.solve_isolated_equilibrium(&thermal, controls, || Ok(())),
        Err(IsolatedThermalError::Residual { .. })
    ));
    assert_eq!(thermal.cells()[&CellId(0)].specific_h_j_kg, 400.0);
}

#[test]
fn isolated_sources_baths_and_cancellation_refuse_without_publishing_history() {
    let law = Law {
        nonlinear_k: false,
        plateau: false,
    };
    for reservoir in [None, Some(300.0)] {
        let (mut thermal, conduction) = setup(&law, &[400.0], reservoir, Vec::new());
        if reservoir.is_none() {
            thermal.cells.get_mut(&CellId(0)).unwrap().boundary.source_w = 1.0;
        }
        let chain = ImplicitConductionChain::new(
            &thermal,
            &conduction,
            &[(MaterialId(0), &law)],
            || Ok(()),
        )
        .unwrap();
        assert!(matches!(
            chain.solve_isolated_equilibrium(&thermal, IsolatedThermalControls::default(), || Ok(
                ()
            )),
            Err(IsolatedThermalError::Unsupported(_))
        ));
        assert!(matches!(
            chain.solve_isolated_equilibrium(&thermal, IsolatedThermalControls::default(), || Err(
                DepositionError::Cancelled
            )),
            Err(IsolatedThermalError::Original(
                ImplicitConductionError::Thermal(DepositionError::Cancelled)
            ))
        ));
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(thermal.cells()[&CellId(0)].specific_h_j_kg, 400.0);
    }
}

struct ClosingConductance;
impl ThermalLaw for ClosingConductance {
    fn temperature_k(&self, energy: f64) -> Result<f64, DepositionError> {
        if !(100.0..=1000.0).contains(&energy) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(energy)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
}
impl ConductivityLaw for ClosingConductance {
    fn conductivity_w_m_k(&self, energy: f64) -> Result<f64, DepositionError> {
        Ok(((self.temperature_k(energy)? - 300.0) / 100.0).powi(2))
    }
}
impl ImplicitConductionLaw for ClosingConductance {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self
    }
    fn conductivity_law(&self) -> &dyn ConductivityLaw {
        self
    }
    fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        [100.0, 1000.0]
    }
    fn state_at_specific_enthalpy(
        &self,
        energy: f64,
    ) -> Result<ImplicitConductionState, DepositionError> {
        Ok(ImplicitConductionState {
            temperature_k: self.temperature_k(energy)?,
            temperature_slope_kg_k_j: 1.0,
            conductivity_w_m_k: self.conductivity_w_m_k(energy)?,
            conductivity_slope_kg_w_m_k_j: (energy - 300.0) / 5000.0,
        })
    }
}

#[test]
fn isolated_component_connectivity_is_rechecked_on_original_final_materials() {
    let law = ClosingConductance;
    let (thermal, conduction) = setup(&law, &[200.0, 400.0], None, Vec::new());
    let chain =
        ImplicitConductionChain::new(&thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap();
    assert!(matches!(
        chain.solve_isolated_equilibrium(&thermal, IsolatedThermalControls::default(), || Ok(())),
        Err(IsolatedThermalError::ConductingTopology { link: 0 })
    ));
    assert_eq!(thermal.cells()[&CellId(0)].specific_h_j_kg, 200.0);
    assert_eq!(thermal.cells()[&CellId(1)].specific_h_j_kg, 400.0);
}
