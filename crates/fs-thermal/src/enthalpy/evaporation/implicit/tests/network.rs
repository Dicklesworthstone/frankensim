//! G0 shared work/atomicity and G1 simultaneous heterogeneous phase surfaces.
use super::*;
use crate::enthalpy::evaporation::implicit::{
    EvaporatingCellBinding, EvaporatingNetworkBackwardEulerStep, EvaporationImplicitControls,
    EvaporationImplicitError, EvaporationNetworkImplicitControls,
};
use std::cell::Cell;

// A second concrete law with different molar mass, latent heat, source range
// and liquid density. These analytic fixtures make no empirical liquid claim.
struct OtherLiquid;
impl EvaporationLaw for OtherLiquid {
    fn temperature_range_k(&self) -> [f64; 2] {
        [285.0, 325.0]
    }
    fn properties(&self, t: f64) -> Result<EvaporationProperties, EvaporationError> {
        Ok(EvaporationProperties {
            saturation_pressure_pa: 6000.0 + 150.0 * (t - 300.0),
            saturation_pressure_derivative_pa_k: 150.0,
            latent_heat_j_kg: 1.1e6 + 1000.0 * (t - 300.0),
            liquid_density_kg_m3: 700.0,
            molar_mass_kg_mol: 0.044,
        })
    }
}
fn other_film() -> GasFilm {
    GasFilm {
        total_pressure_pa: 120_000.0,
        vapor_partial_pressure_pa: 9000.0,
        temperature_k: 305.0,
        mass_transfer_coefficient_m_s: 0.007,
    }
}
fn other_transfer(t: f64, area: f64, dt: f64) -> f64 {
    let film = other_film();
    film.mass_transfer_coefficient_m_s * 0.044 * film.total_pressure_pa
        / (MOLAR_GAS_CONSTANT * film.temperature_k)
        * ((film.total_pressure_pa - film.vapor_partial_pressure_pa)
            / (film.total_pressure_pa - 6000.0 - 150.0 * (t - 300.0)))
            .ln()
        * area
        * dt
}
fn binding<'a>(id: u64, law: &'a dyn EvaporationLaw, film: GasFilm) -> EvaporatingCellBinding<'a> {
    EvaporatingCellBinding {
        id: CellId(id),
        law,
        film,
        pure_liquid_enthalpy_range_j_kg: range(0.0),
    }
}
fn empty_network(
    _: &DepositionSimulation<'_>,
    _: &BTreeMap<CellId, LiquidSurfaceState>,
) -> Result<BoundaryExchangeStepCoefficients<()>, ThermalEvaporationError> {
    Ok(coefficients(()))
}

#[test]
fn heterogeneous_surfaces_share_one_solve_with_cross_inventory_heat_and_canonical_order() {
    let dt = 1000.0;
    let target = [305.0, 306.0, 300.0];
    let dm = [
        transfer(&VARYING, target[0], film(0.0), 0.007, dt),
        other_transfer(target[1], 0.007, dt),
    ];
    assert!(dm[0] > 0.0 && dm[1] < 0.0);
    let latent = [LATENT - 2000.0 * 5.0, 1.1e6 + 1000.0 * 6.0];
    let required = [
        0.1 * CP * 5.0 + dm[0] * latent[0],
        -0.2 * CP * 4.0 + dm[1] * latent[1],
        0.15 * CP * 10.0,
    ];
    let coupled_sources = [-100.0 * dm[1], 80.0 * dm[0]];
    let source = [
        required[0] / dt - (0.07 * -15.0 + 0.15 * 15.0 + 0.4) - coupled_sources[0],
        required[1] / dt - (0.14 * -6.0 - 0.4 - 0.2 * 6.0) - coupled_sources[1],
        required[2] / dt - 0.2 * 6.0,
    ];
    let shifted = Caloric {
        datum: 1e9,
        fixed_mass: None,
    };
    let mut a = cell(1, 0.1, 300.0, 0.0, source[0]);
    a.boundary.ambient_k = 290.0;
    a.boundary.htc_w_m2_k = 10.0;
    a.boundary.reservoir = Some(Reservoir {
        area_m2: 0.001,
        temperature_k: 320.0,
        conductance_w_k: 0.15,
    });
    let mut b = cell(2, 0.2, 310.0, 1e9, source[1]);
    b.material_id = MaterialId(2);
    b.boundary.htc_w_m2_k = 20.0;
    let c = cell(3, 0.15, 290.0, 0.0, source[2]);
    let mut events: Vec<_> = [a, b, c]
        .into_iter()
        .map(|cell| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell),
        })
        .collect();
    for (id, a, b, area, conductance) in [(9, 1, 2, 0.002, 0.4), (10, 2, 3, 0.001, 0.2)] {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id,
                a: CellId(a),
                b: CellId(b),
                area_m2: area,
                conductance_w_k: conductance,
            }),
        });
    }
    let old = DepositionSimulation::new(
        events,
        &[(MaterialId(1), &CALORIC), (MaterialId(2), &shifted)],
        DepositionBudget {
            max_cells: 3,
            max_links: 2,
            max_events: 5,
            max_substeps: 10,
        },
    )
    .unwrap()
    .advance_to(0.0, config(), || Ok(()))
    .unwrap()
    .0;
    let second = EvaporatingCellBinding {
        pure_liquid_enthalpy_range_j_kg: [1e9 - 10.0 * CP, 1e9 + 20.0 * CP],
        ..binding(2, &OtherLiquid, other_film())
    };
    let bindings = [binding(1, &VARYING, film(0.0)), second];
    let mut accepted = Vec::new();
    for order in [bindings, [bindings[1], bindings[0]]] {
        let calls = Cell::new(0usize);
        let (next, payload, report) = old
            .solve_implicit_evaporating_network_step(
                &order,
                dt,
                EvaporationNetworkImplicitControls::default(),
                |thermal, liquids| {
                    calls.set(calls.get() + 1);
                    assert_eq!(
                        liquids.keys().copied().collect::<Vec<_>>(),
                        [CellId(1), CellId(2)]
                    );
                    assert_eq!(thermal.cells()[&CellId(1)].mass_kg, 0.1);
                    assert_eq!(thermal.cells()[&CellId(2)].mass_kg, 0.2);
                    let mut result = coefficients(calls.get());
                    result.transport.heat.bulk.sources_w = BTreeMap::from([
                        (CellId(1), 100.0 * (liquids[&CellId(2)].mass_kg - 0.2)),
                        (CellId(2), -80.0 * (liquids[&CellId(1)].mass_kg - 0.1)),
                    ]);
                    Ok(result)
                },
                || Ok(()),
            )
            .unwrap();
        for (i, &temperature) in target.iter().enumerate() {
            let id = CellId(i as u64 + 1);
            close(next.temperature_k(id).unwrap(), temperature, 2e-9);
            assert_eq!(next.cells()[&id].boundary, old.cells()[&id].boundary);
            assert_eq!(next.cells()[&id].geometry, old.cells()[&id].geometry);
        }
        assert_eq!(next.cells()[&CellId(3)].mass_kg, 0.15);
        let audit = &report.physical;
        for i in 0..2 {
            let receipt = &audit.surfaces[&CellId(i as u64 + 1)];
            close(receipt.transferred_mass_kg, dm[i], 1e-12);
            close(receipt.latent_energy_j, dm[i] * latent[i], 1e-6);
            close(receipt.liquid.area_m2, 0.007, 1e-17);
        }
        let vapor = audit
            .surfaces
            .values()
            .map(|receipt| receipt.vapor_energy_j)
            .sum::<f64>();
        let expected_heat = required.iter().sum::<f64>();
        close(audit.external_heat_j, expected_heat, 1e-6);
        close(
            next.energy_j() - old.energy_j() + vapor,
            expected_heat,
            1e-5,
        );
        close(
            audit.stored_energy_change_j,
            next.energy_j() - old.energy_j(),
            1e-5,
        );
        close(
            next.mass_kg() + dm.iter().sum::<f64>(),
            old.mass_kg(),
            1e-12,
        );
        assert_eq!(payload, calls.get());
        assert_eq!(report.residual_evaluations, calls.get());
        assert_eq!(next.links(), old.links());
        accepted.push((next.cells().clone(), report));
    }
    assert_eq!(accepted[0], accepted[1]);
    assert_eq!(old.time_s(), 0.0);
}

#[test]
fn single_unsized_api_delegates_without_losing_its_explicit_liquid_domain() {
    let law: &dyn EvaporationLaw = &VARYING;
    let old = simulation(&CALORIC, vec![cell(1, 0.1, 300.0, 0.0, 15.0)], None, vec![]);
    let single = EvaporationImplicitControls::new(range(0.0));
    let (a, (), single_report) = old
        .solve_implicit_evaporating_step(
            CellId(1),
            law,
            film(0.0),
            1000.0,
            single,
            empty,
            || Ok(()),
        )
        .unwrap();
    let (b, (), network_report) = old
        .solve_implicit_evaporating_network_step(
            &[binding(1, law, film(0.0))],
            1000.0,
            single.network_controls(),
            empty_network,
            || Ok(()),
        )
        .unwrap();
    assert_eq!(a.cells(), b.cells());
    assert_eq!(
        single_report.physical.surface,
        network_report.physical.surfaces[&CellId(1)]
    );
    assert_eq!(
        single_report.physical.external_heat_j,
        network_report.physical.external_heat_j
    );
    assert_eq!(
        single_report.physical.energy_residual_j,
        network_report.physical.energy_residual_j
    );
    assert_eq!(
        single_report.newton_iterations,
        network_report.newton_iterations
    );
    assert_eq!(
        single_report.residual_evaluations,
        network_report.residual_evaluations
    );
    let restricted = EvaporationImplicitControls {
        pure_liquid_enthalpy_range_j_kg: [0.0, CP],
        ..single
    };
    assert!(matches!(
        old.solve_implicit_evaporating_step(
            CellId(1),
            law,
            film(0.0),
            1000.0,
            restricted,
            empty,
            || Ok(()),
        ),
        Err(EvaporationImplicitError::Physical(
            ThermalEvaporationError::Thermal(DepositionError::MaterialDomain)
        ))
    ));
}

#[test]
fn opposite_surface_defects_do_not_cancel_the_local_energy_or_mass_gates() {
    let sealed = GasFilm {
        mass_transfer_coefficient_m_s: 0.0,
        ..film(0.0)
    };
    let old = simulation(
        &CALORIC,
        vec![cell(1, 0.1, 300.0, 0.0, 0.0), cell(2, 0.1, 300.0, 0.0, 0.0)],
        None,
        vec![],
    );
    let step = EvaporatingNetworkBackwardEulerStep::new(
        &old,
        &[binding(1, &CONSTANT, sealed), binding(2, &CONSTANT, sealed)],
        1.0,
        || Ok(()),
    )
    .unwrap();
    let trial = step
        .evaluate(&[1.0, -1.0], empty_network, || Ok(()))
        .unwrap();
    assert_eq!(trial.audit().energy_residual_j, 0.0);
    assert_eq!(trial.maximum_cell_residual_j(), 1.0);
    assert!(matches!(
        trial.accept(1e-9, 1e-9, 1e-14, || Ok(())),
        Err(ThermalEvaporationError::Thermal(
            DepositionError::EnergyBalance(_)
        ))
    ));

    // Reciprocal carrier-pressure ratios give exactly opposing source-owned
    // fluxes. Their represented-rate defects cancel globally but still refuse.
    let power = transfer(&CONSTANT, 300.0, film(0.0), 0.01, 1.0) * LATENT;
    let old = simulation(
        &CALORIC,
        vec![
            cell(1, 0.1, 300.0, 0.0, power),
            cell(2, 0.1, 300.0, 0.0, -power),
        ],
        None,
        vec![],
    );
    let step = EvaporatingNetworkBackwardEulerStep::new(
        &old,
        &[
            binding(1, &CONSTANT, film(0.0)),
            binding(2, &CONSTANT, film(7840.0)),
        ],
        1.0,
        || Ok(()),
    )
    .unwrap();
    let trial = step
        .evaluate(&[0.0, 0.0], empty_network, || Ok(()))
        .unwrap();
    let defects: Vec<_> = trial
        .audit()
        .surfaces
        .values()
        .map(|s| s.mass_rate_residual_kg)
        .collect();
    assert!(defects.iter().all(|d| *d != 0.0));
    close(defects.iter().sum::<f64>(), 0.0, 1e-20);
    assert_eq!(trial.audit().mass_balance_residual_kg, 0.0);
    assert!(matches!(
        trial.accept(1e-7, 1e-7, 0.0, || Ok(())),
        Err(ThermalEvaporationError::MassBalance { .. })
    ));
}

#[test]
fn batch_rate_work_is_linear_in_cells_and_final_cancellation_discards_all_masses() {
    struct CountingCaloric<'a>(&'a Cell<usize>);
    impl ThermalLaw for CountingCaloric<'_> {
        fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
            self.0.set(self.0.get() + 1);
            CALORIC.temperature_k(h)
        }
        fn maximum_temperature_slope(&self) -> f64 {
            1.0 / CP
        }
    }
    let queries = Cell::new(0);
    let caloric = CountingCaloric(&queries);
    let power = transfer(&CONSTANT, 300.0, film(0.0), 0.01, 1.0) * LATENT;
    let old = simulation(
        &caloric,
        (1..=4).map(|id| cell(id, 0.1, 300.0, 0.0, power)).collect(),
        None,
        vec![],
    );
    let bindings: Vec<_> = (1..=4)
        .map(|id| binding(id, &CONSTANT, film(0.0)))
        .collect();
    let step = EvaporatingNetworkBackwardEulerStep::new(&old, &bindings, 1.0, || Ok(())).unwrap();
    queries.set(0);
    let calls = Cell::new(0);
    let trial = step
        .evaluate(
            &[0.0; 4],
            |_, liquids| {
                calls.set(calls.get() + 1);
                assert_eq!(liquids.len(), 4);
                assert!(liquids.values().all(|liquid| liquid.mass_kg < 0.1));
                Ok::<_, ThermalEvaporationError>(coefficients(()))
            },
            || Ok(()),
        )
        .unwrap();
    assert_eq!(calls.get(), 1);
    assert!(
        queries.get() <= 4 * old.cells().len(),
        "repeated full heat passes: {} queries",
        queries.get()
    );
    let mut polls = 0;
    let (accepted, (), _) = trial
        .accept(1e-7, 1e-7, 1e-12, || {
            polls += 1;
            Ok(())
        })
        .unwrap();
    assert!(accepted.cells().values().all(|cell| cell.mass_kg < 0.1));
    let trial = step.evaluate(&[0.0; 4], empty_network, || Ok(())).unwrap();
    let mut count = 0;
    assert!(matches!(
        trial.accept(1e-7, 1e-7, 1e-12, || {
            count += 1;
            if count == polls {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(ThermalEvaporationError::Thermal(DepositionError::Cancelled))
    ));
    assert!(old.cells().values().all(|cell| cell.mass_kg == 0.1));
    assert_eq!(old.time_s(), 0.0);
}

#[test]
fn duplicate_inactive_or_material_refused_binding_cannot_publish_other_surfaces() {
    let fixed = Caloric {
        datum: 0.0,
        fixed_mass: Some(0.1),
    };
    let a = cell(1, 0.1, 300.0, 0.0, 0.0);
    let mut b = cell(2, 0.1, 300.0, 0.0, 0.0);
    b.material_id = MaterialId(2);
    let old = DepositionSimulation::new(
        vec![
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(a),
            },
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(b),
            },
        ],
        &[(MaterialId(1), &CALORIC), (MaterialId(2), &fixed)],
        DepositionBudget {
            max_cells: 2,
            max_links: 0,
            max_events: 2,
            max_substeps: 4,
        },
    )
    .unwrap()
    .advance_to(0.0, config(), || Ok(()))
    .unwrap()
    .0;
    let a = binding(1, &VARYING, film(0.0));
    assert!(matches!(
        EvaporatingNetworkBackwardEulerStep::new(&old, &[a, a], 1.0, || Ok(())),
        Err(ThermalEvaporationError::Thermal(
            DepositionError::InvalidInput("duplicate evaporating cell binding")
        ))
    ));
    assert!(matches!(
        EvaporatingNetworkBackwardEulerStep::new(
            &old,
            &[binding(3, &VARYING, film(0.0))],
            1.0,
            || Ok(()),
        ),
        Err(ThermalEvaporationError::Thermal(
            DepositionError::InvalidInput("evaporation requires an active liquid cell")
        ))
    ));
    let step = EvaporatingNetworkBackwardEulerStep::new(
        &old,
        &[a, binding(2, &VARYING, film(0.0))],
        1.0,
        || Ok(()),
    )
    .unwrap();
    let called = Cell::new(false);
    let result = step.evaluate(
        &[0.0, 0.0],
        |_, _| {
            called.set(true);
            Ok::<_, ThermalEvaporationError>(coefficients(()))
        },
        || Ok(()),
    );
    assert!(matches!(
        result,
        Err(ThermalEvaporationError::Thermal(
            DepositionError::InvalidInput("original fixed inventory changed")
        ))
    ));
    assert!(!called.get());
    assert!(old.cells().values().all(|cell| cell.mass_kg == 0.1));
    assert_eq!(old.time_s(), 0.0);
}

#[test]
fn empty_phase_selection_keeps_the_original_implicit_heat_network() {
    let old = simulation(
        &CALORIC,
        vec![cell(1, 0.1, 300.0, 0.0, 0.0), cell(2, 0.1, 310.0, 0.0, 0.0)],
        Some(2.0),
        vec![],
    );
    let bindings: &[EvaporatingCellBinding<'_>] = &[];
    let (next, (), report) = old
        .solve_implicit_evaporating_network_step(
            bindings,
            1000.0,
            EvaporationNetworkImplicitControls::default(),
            |_, liquids| {
                assert!(liquids.is_empty());
                Ok(coefficients(()))
            },
            || Ok(()),
        )
        .unwrap();
    let exchange = 20_000.0 / 11.0;
    close(
        next.temperature_k(CellId(1)).unwrap(),
        300.0 + exchange / 400.0,
        1e-9,
    );
    close(
        next.temperature_k(CellId(2)).unwrap(),
        310.0 - exchange / 400.0,
        1e-9,
    );
    assert!(report.physical.surfaces.is_empty());
    assert_eq!(report.physical.external_heat_j, 0.0);
    close(next.energy_j(), old.energy_j(), 1e-8);
    assert_eq!(next.mass_kg(), old.mass_kg());
}

// Independent normal-scale pressure bisection with standard-library ln. The
// production implementation instead solves a scaled pressure-drop fraction.
fn kinetic_flux_reference(ps: f64, t: f64, alpha: f64, bath: GasFilm) -> f64 {
    let k = alpha * (0.018 / (2.0 * core::f64::consts::PI * MOLAR_GAS_CONSTANT * t)).sqrt();
    let c = bath.mass_transfer_coefficient_m_s * 0.018 * bath.total_pressure_pa
        / (MOLAR_GAS_CONSTANT * bath.temperature_k);
    let mut lo = ps.min(bath.vapor_partial_pressure_pa);
    let mut hi = ps.max(bath.vapor_partial_pressure_pa);
    for _ in 0..100 {
        let pi = 0.5 * (lo + hi);
        let residual = c
            * ((bath.total_pressure_pa - bath.vapor_partial_pressure_pa)
                / (bath.total_pressure_pa - pi))
                .ln()
            - k * (ps - pi);
        if residual < 0.0 {
            lo = pi;
        } else {
            hi = pi;
        }
    }
    k * (ps - 0.5 * (lo + hi))
}

#[test]
fn finite_interface_kinetics_reaches_the_original_shared_heat_network() {
    use crate::evaporation::{InterfaceKinetics, KineticEvaporation};
    let alpha = 1e-6;
    let law = KineticEvaporation::try_new(
        &VARYING,
        InterfaceKinetics::HertzKnudsen {
            accommodation_coefficient: alpha,
        },
    )
    .unwrap();
    let sealed = KineticEvaporation::try_new(
        &VARYING,
        InterfaceKinetics::Schrage {
            accommodation_coefficient: 0.0,
        },
    )
    .unwrap();
    let dt = 1000.0;
    let expected_dm = kinetic_flux_reference(4500.0, 305.0, alpha, film(0.0)) * 0.01 * dt;
    let required_heat = 0.1 * CP * 5.0 + expected_dm * (LATENT - 2000.0 * 5.0);
    let old = simulation(
        &CALORIC,
        vec![
            cell(1, 0.1, 300.0, 0.0, required_heat / dt),
            cell(2, 0.1, 310.0, 0.0, 0.0),
        ],
        None,
        vec![],
    );
    let bindings = [binding(1, &law, film(0.0)), binding(2, &sealed, film(0.0))];
    let (next, (), report) = old
        .solve_implicit_evaporating_network_step(
            &bindings,
            dt,
            EvaporationNetworkImplicitControls::default(),
            empty_network,
            || Ok(()),
        )
        .unwrap();
    close(next.temperature_k(CellId(1)).unwrap(), 305.0, 2e-8);
    close(
        report.physical.surfaces[&CellId(1)].transferred_mass_kg,
        expected_dm,
        2e-12,
    );
    assert_eq!(next.cells()[&CellId(2)].mass_kg, 0.1);
    assert_eq!(
        report.physical.surfaces[&CellId(2)].transferred_mass_kg,
        0.0
    );
    assert!(report.physical.energy_residual_j.abs() < 1e-6);
    assert!(report.physical.mass_balance_residual_kg.abs() < 1e-14);
    assert_eq!(old.time_s(), 0.0);
    assert_eq!(old.cells()[&CellId(1)].mass_kg, 0.1);
}

#[test]
fn inventory_retry_uses_the_selected_initial_and_endpoint_kinetic_rates() {
    use crate::evaporation::{InterfaceKinetics, KineticEvaporation};
    let alpha = 1e-6;
    let law = KineticEvaporation::try_new(
        &CONSTANT,
        InterfaceKinetics::HertzKnudsen {
            accommodation_coefficient: alpha,
        },
    )
    .unwrap();
    let rate = kinetic_flux_reference(4000.0, 300.0, alpha, film(0.0)) * 0.01;
    let expected_limit = 0.1 / rate;
    let old = simulation(&CALORIC, vec![cell(1, 0.1, 300.0, 0.0, 0.0)], None, vec![]);
    let bindings = [binding(1, &law, film(0.0))];
    let step =
        EvaporatingNetworkBackwardEulerStep::new(&old, &bindings, 2.0 * expected_limit, || Ok(()))
            .unwrap();
    let result = step.evaluate(&[0.0], empty_network, || Ok(()));
    let Err(ThermalEvaporationError::Surface(EvaporationError::InventoryLimit {
        maximum_duration_s,
    })) = result
    else {
        panic!("the physical endpoint must refuse exhausted inventory");
    };
    close(maximum_duration_s, expected_limit, expected_limit * 2e-13);
    assert_eq!(old.cells()[&CellId(1)].mass_kg, 0.1);
}

#[test]
fn parallel_film_hook_reaches_the_original_network_endpoint_after_heating() {
    use crate::evaporation::{
        EvaporationPatch, InterfaceKinetics, KineticEvaporation, ParallelFilmEvaporation,
    };
    let patches = [
        EvaporationPatch {
            area_m2: 0.002,
            mass_transfer_coefficient_m_s: 1e-4,
        },
        EvaporationPatch {
            area_m2: 0.008,
            mass_transfer_coefficient_m_s: 0.2,
        },
    ];
    let parallel = ParallelFilmEvaporation::try_new(&VARYING, &patches).unwrap();
    let alpha = 1e-5;
    let law = KineticEvaporation::try_new(
        &parallel,
        InterfaceKinetics::HertzKnudsen {
            accommodation_coefficient: alpha,
        },
    )
    .unwrap();
    let duration = 1000.0;
    let expected_dm: f64 = patches
        .iter()
        .map(|patch| {
            let bath = GasFilm {
                mass_transfer_coefficient_m_s: patch.mass_transfer_coefficient_m_s,
                ..film(0.0)
            };
            kinetic_flux_reference(4500.0, 305.0, alpha, bath) * patch.area_m2 * duration
        })
        .sum();
    let supplied_heat = 0.1 * CP * 5.0 + expected_dm * (LATENT - 2000.0 * 5.0);
    let old = simulation(
        &CALORIC,
        vec![cell(1, 0.1, 300.0, 0.0, supplied_heat / duration)],
        None,
        vec![],
    );
    let bindings = [binding(1, &law, film(0.0))];
    let (next, (), report) = old
        .solve_implicit_evaporating_network_step(
            &bindings,
            duration,
            EvaporationNetworkImplicitControls::default(),
            empty_network,
            || Ok(()),
        )
        .unwrap();
    close(next.temperature_k(CellId(1)).unwrap(), 305.0, 2e-8);
    close(
        report.physical.surfaces[&CellId(1)].transferred_mass_kg,
        expected_dm,
        2e-12,
    );
    assert!(report.physical.energy_residual_j.abs() < 1e-6);
    assert!(report.physical.mass_balance_residual_kg.abs() < 1e-14);
    assert_eq!(old.cells()[&CellId(1)].mass_kg, 0.1);
}

#[test]
fn parallel_film_hook_also_controls_the_cached_initial_inventory_retry() {
    use crate::evaporation::{
        EvaporationPatch, InterfaceKinetics, KineticEvaporation, ParallelFilmEvaporation,
    };
    let patches = [
        EvaporationPatch {
            area_m2: 0.002,
            mass_transfer_coefficient_m_s: 1e-4,
        },
        EvaporationPatch {
            area_m2: 0.008,
            mass_transfer_coefficient_m_s: 0.2,
        },
    ];
    let alpha = 1e-5;
    let selected = KineticEvaporation::try_new(
        &CONSTANT,
        InterfaceKinetics::HertzKnudsen {
            accommodation_coefficient: alpha,
        },
    )
    .unwrap();
    let parallel = ParallelFilmEvaporation::try_new(&selected, &patches).unwrap();
    let rate: f64 = patches
        .iter()
        .map(|patch| {
            let bath = GasFilm {
                mass_transfer_coefficient_m_s: patch.mass_transfer_coefficient_m_s,
                ..film(0.0)
            };
            kinetic_flux_reference(4000.0, 300.0, alpha, bath) * patch.area_m2
        })
        .sum();
    let expected_limit = 0.1 / rate;
    let old = simulation(&CALORIC, vec![cell(1, 0.1, 300.0, 0.0, 0.0)], None, vec![]);
    let bindings = [binding(1, &parallel, film(0.0))];
    let step =
        EvaporatingNetworkBackwardEulerStep::new(&old, &bindings, 2.0 * expected_limit, || Ok(()))
            .unwrap();
    let result = step.evaluate(&[0.0], empty_network, || Ok(()));
    let Err(ThermalEvaporationError::Surface(EvaporationError::InventoryLimit {
        maximum_duration_s,
    })) = result
    else {
        panic!("parallel patch inventory must refuse before any endpoint publication");
    };
    close(maximum_duration_s, expected_limit, expected_limit * 2e-11);
    assert_eq!(old.cells()[&CellId(1)].mass_kg, 0.1);
}
