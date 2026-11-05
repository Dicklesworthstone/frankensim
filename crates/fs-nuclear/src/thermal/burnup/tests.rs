// Child of thermal::tests: these manufactured coefficients test physical
// coupling equations, not evaluated isotope data or reactor calibration.
use super::*;
use crate::burnup::{
    BurnupBudget, BurnupData, BurnupInventory, BurnupLimits, BurnupModel, BurnupNuclide,
    BurnupProduct, BurnupRate, BurnupReaction,
};
use crate::kinetics::{
    CompositionFeedback, KineticsBudget, KineticsError, PointKineticsModel, PointKineticsState,
};
use fs_la::decay::DecayLimits;
use fs_thermochem::decay::DecayLawV1;

fn burnup_budget() -> BurnupHeatingBudget {
    BurnupHeatingBudget {
        source: BurnupBudget {
            limits: BurnupLimits {
                max_neutron_groups: 4,
                max_reactions: 16,
                max_products: 32,
                propagation: DecayLimits::default(),
            },
            nuclei_tolerance: 1e-8,
            energy_tolerance_j: 1e-8,
        },
        max_composition_feedback_terms: 4,
    }
}

fn reaction(
    parent: usize,
    rate: BurnupRate,
    products: &[(usize, f64)],
    energy: EnergyChannels,
) -> BurnupReaction {
    BurnupReaction {
        data_id: format!(
            "manufactured-product-reaction-{parent}-{}",
            energy.deposited_j
        ),
        parent,
        rate,
        products: products
            .iter()
            .map(|&(nuclide, yield_per_event)| BurnupProduct {
                nuclide,
                yield_per_event,
            })
            .collect(),
        energy_per_event: energy,
    }
}

fn deposited(value: f64) -> EnergyChannels {
    EnergyChannels {
        deposited_j: value,
        ..EnergyChannels::default()
    }
}

fn product_source(shutdown: Option<f64>, product_decay: bool, capture: bool) -> BurnupCellSource {
    let groups = Arc::new(
        NeutronGroups::new("manufactured-product-grid-v1".into(), vec![0.0, 1.0]).unwrap(),
    );
    let flux = GroupFlux::new(groups.clone(), vec![1.0]).unwrap();
    // These identifiers label manufactured laws; no physical half-life or
    // cross section for the named nuclides is asserted by this fixture.
    let identities = [
        (92, 235),
        (54, 135),
        (38, 98),
        (55, 135),
        (54, 136),
        (92, 236),
    ];
    let nuclides = identities
        .into_iter()
        .enumerate()
        .map(|(i, (z, a))| BurnupNuclide {
            id: NuclideId::new(z, a, 0).unwrap(),
            law: if i == 1 && product_decay {
                DecayLawV1::from_decay_constant(fs_qty::Frequency::new(0.7)).unwrap()
            } else {
                DecayLawV1::stable()
            },
        })
        .collect();
    let mut reactions = vec![reaction(
        0,
        BurnupRate::Fission {
            cross_sections_m2: vec![0.2],
        },
        &[(1, 1.0), (2, 1.0)],
        deposited(8.0),
    )];
    if product_decay {
        reactions.push(reaction(
            1,
            BurnupRate::Decay {
                branch_fraction: 1.0,
            },
            &[(3, 1.0)],
            EnergyChannels {
                deposited_j: 4.0,
                radiation_j: 1.0,
                escaped_j: 2.0,
            },
        ));
    }
    if capture {
        reactions.push(reaction(
            0,
            BurnupRate::Capture {
                cross_sections_m2: vec![0.1],
            },
            &[(5, 1.0)],
            deposited(3.0),
        ));
        reactions.push(reaction(
            1,
            BurnupRate::Capture {
                cross_sections_m2: vec![0.3],
            },
            &[(4, 1.0)],
            deposited(2.0),
        ));
    }
    let model = Arc::new(
        BurnupModel::new(
            BurnupData {
                data_id: "manufactured-product-network-v1".into(),
                groups: groups.clone(),
                temperature_min_k: 290.0,
                temperature_max_k: 1000.0,
                nuclides,
                reactions,
            },
            burnup_budget().source.limits,
        )
        .unwrap(),
    );
    let changes = shutdown
        .map(|time_s| {
            vec![NeutronFluxChange {
                time_s,
                flux: GroupFlux::new(groups, vec![0.0]).unwrap(),
            }]
        })
        .unwrap_or_default();
    BurnupCellSource {
        cell_id: CellId(1),
        material_id: MaterialId(1),
        inventory: BurnupInventory::new(model, vec![100.0, 0.0, 0.0, 0.0, 0.0, 0.0]).unwrap(),
        fluxes: NeutronFluxSchedule::new(flux, changes).unwrap(),
        composition_feedback: vec![],
    }
}

#[test]
fn products_react_and_decay_inside_the_same_thermal_interval() {
    let initial = thermal(
        vec![cell(1, 1, 2.0, 1200.0, 3.0)],
        &[(MaterialId(1), &HOT)],
        vec![],
    );
    let heating = BurnupHeating::new(
        &initial,
        vec![product_source(None, true, true)],
        burnup_budget(),
    )
    .unwrap();
    let (heated, final_state, audit) = heating
        .advance_to(&initial, 1.0, controls(10.0), || Ok(()))
        .unwrap();
    let a = 0.3_f64;
    let b = 1.0_f64;
    let fissions = 100.0 * 0.2 / a * -(-a).exp_m1();
    let parent_captures = 100.0 * 0.1 / a * -(-a).exp_m1();
    let product = 100.0 * 0.2 * ((-a).exp() - (-b).exp()) / (b - a);
    let product_time = 100.0 * 0.2 / (b - a) * (-(-a).exp_m1() / a + (-b).exp_m1() / b);
    let decays = 0.7 * product_time;
    let product_captures = 0.3 * product_time;
    let release = 8.0 * fissions + 3.0 * parent_captures + 4.0 * decays + 2.0 * product_captures;
    let source = final_state.source(CellId(1)).unwrap();
    let n = source.inventory.populations();
    for (&actual, expected) in n.iter().zip([
        100.0 * (-a).exp(),
        product,
        fissions,
        decays,
        product_captures,
        parent_captures,
    ]) {
        close(actual, expected, 3e-10);
    }
    close(final_state.audit().fissions, fissions, 3e-10);
    close(
        final_state.audit().captures,
        parent_captures + product_captures,
        3e-10,
    );
    close(final_state.audit().released.deposited_j, release, 3e-9);
    close(final_state.audit().released.radiation_j, decays, 3e-10);
    close(final_state.audit().decays, decays, 3e-10);
    close(final_state.audit().released.escaped_j, 2.0 * decays, 3e-10);
    close(audit.external_heat_j, release + 3.0, 3e-9);
    close(heated.energy_j() - initial.energy_j(), release + 3.0, 3e-9);
    assert_eq!(heated.cells()[&CellId(1)].mass_kg, 2.0);
    assert_eq!(heated.cells()[&CellId(1)].boundary.source_w, 3.0);
    assert_eq!(audit.substeps, 1);
    // Fission multiplicity changes tracked nucleus count; the solver must
    // not normalize the two fission products into a single conserved nucleus.
    close(n.iter().sum::<f64>(), 100.0 + fissions, 3e-10);
}

#[test]
fn shutdown_stops_neutron_reactions_while_created_products_keep_heating() {
    let initial = thermal(
        vec![cell(1, 1, 2.0, 1200.0, 0.0)],
        &[(MaterialId(1), &HOT)],
        vec![],
    );
    let source = BurnupHeating::new(
        &initial,
        vec![product_source(Some(0.4), true, true)],
        burnup_budget(),
    )
    .unwrap();
    let (at_shutdown, stopped, _) = source
        .advance_to(&initial, 0.4, controls(10.0), || Ok(()))
        .unwrap();
    let p0 = stopped.source(CellId(1)).unwrap().inventory.populations()[1];
    let f0 = stopped.audit().fissions;
    let c0 = stopped.audit().captures;
    let (cooled, ended, _) = stopped
        .advance_to(&at_shutdown, 1.0, controls(10.0), || Ok(()))
        .unwrap();
    let decay_count = p0 * -(-0.7_f64 * 0.6).exp_m1();
    close(ended.audit().fissions, 0.0, 3e-10);
    close(ended.audit().captures, 0.0, 3e-10);
    let cumulative = ended
        .source(CellId(1))
        .unwrap()
        .inventory
        .cumulative_reaction_counts();
    close(cumulative[0], f0, 3e-10);
    close(cumulative[2] + cumulative[3], c0, 3e-10);
    close(
        ended.last_step(CellId(1)).unwrap().released.deposited_j,
        4.0 * decay_count,
        3e-10,
    );
    close(
        cooled.energy_j() - at_shutdown.energy_j(),
        4.0 * decay_count,
        3e-9,
    );
    let (direct, direct_source, _) = source
        .advance_to(&initial, 1.0, controls(10.0), || Ok(()))
        .unwrap();
    assert_eq!(direct_source.audit().substeps, 2);
    close(direct.energy_j(), cooled.energy_j(), 3e-9);
    for (a, b) in direct_source
        .source(CellId(1))
        .unwrap()
        .inventory
        .populations()
        .iter()
        .zip(ended.source(CellId(1)).unwrap().inventory.populations())
    {
        close(*a, *b, 3e-10);
    }
}

fn poison_fit(maximum_nuclei: f64) -> CompositionFeedback {
    CompositionFeedback::new(
        "manufactured-poison-worth-v1".into(),
        NuclideId::new(54, 135, 0).unwrap(),
        0.0,
        0.0,
        maximum_nuclei,
        -0.02,
    )
    .unwrap()
}

fn burnup_kinetics_budget() -> KineticHeatingBudget<BurnupHeatingBudget> {
    KineticHeatingBudget {
        irradiation: burnup_budget(),
        kinetics: KineticsBudget {
            max_delayed_groups: 4,
            max_squarings: 32,
            max_taylor_terms: 24,
            population_tolerance: 1e-8,
            integral_tolerance_population_s: 1e-8,
        },
        max_feedback_terms: 4,
    }
}

fn poison_source(with_feedback: bool) -> KineticCellSource<BurnupCellSource> {
    let mut irradiation = product_source(None, false, false);
    if with_feedback {
        irradiation.composition_feedback.push(poison_fit(100.0));
    }
    KineticCellSource {
        irradiation,
        kinetics: PointKineticsState::new(
            Arc::new(
                PointKineticsModel::new("manufactured-poison-region-v1".into(), 1.0, vec![])
                    .unwrap(),
            ),
            1.0,
            vec![],
        )
        .unwrap(),
        controls: KineticsControlSchedule::new(
            KineticsControl {
                reactivity: 0.05,
                external_source_per_s: 0.0,
            },
            vec![],
        )
        .unwrap(),
        feedback: vec![],
    }
}

#[test]
fn product_worth_closes_neutron_feedback_and_refines_to_an_independent_invariant() {
    let initial = thermal(
        vec![cell(1, 1, 2.0, 1200.0, 0.0)],
        &[(MaterialId(1), &HOT)],
        vec![],
    );
    let heating = KineticHeating::new(
        &initial,
        vec![poison_source(true)],
        burnup_kinetics_budget(),
    )
    .unwrap();
    let unpoisoned = KineticHeating::new(
        &initial,
        vec![poison_source(false)],
        burnup_kinetics_budget(),
    )
    .unwrap();
    let (_, unpoisoned, _) = unpoisoned
        .advance_to(&initial, 1.0, controls(0.025), || Ok(()))
        .unwrap();
    let mut errors = Vec::new();
    for dt in [0.1, 0.05, 0.025] {
        let (heated, state, _) = heating
            .advance_to(&initial, 1.0, controls(dt), || Ok(()))
            .unwrap();
        let source = state.source(CellId(1)).unwrap();
        let fuel = source.irradiation.inventory.populations()[0];
        let fluence = -(fuel / 100.0).ln() / 0.2;
        // Eliminate t from F'=-0.2*n*F, P=100-F,
        // n'=(0.05-0.02*P)*n. This invariant is independent of our solver.
        let expected_n = 1.0 + 0.05 * fluence - 2.0 * (fluence + (-0.2 * fluence).exp_m1() / 0.2);
        errors.push((source.kinetics.neutron_population() - expected_n).abs());
        assert!(source.kinetics.neutron_population() < 1.0);
        assert!(state.audit().fissions < unpoisoned.audit().fissions);
        assert!(state.last_step(CellId(1)).unwrap().applied_reactivity < 0.05);
        close(
            heated.energy_j() - initial.energy_j(),
            state.audit().released.deposited_j,
            3e-8,
        );
    }
    for pair in errors.windows(2) {
        assert!(
            (1.8..2.2).contains(&(pair[0] / pair[1])),
            "composition refinement: {errors:?}"
        );
    }
    assert!(errors[2] < 0.006, "composition refinement: {errors:?}");
}

#[test]
fn composition_domain_failure_and_cancellation_publish_no_heat_or_products() {
    let initial = thermal(
        vec![cell(1, 1, 2.0, 1200.0, 0.0)],
        &[(MaterialId(1), &HOT)],
        vec![],
    );
    let mut limited = product_source(None, true, true);
    limited.composition_feedback.push(poison_fit(0.1));
    let heating = BurnupHeating::new(&initial, vec![limited], burnup_budget()).unwrap();
    assert!(matches!(
        heating.advance_to(&initial, 1.0, controls(1.0), || Ok(())),
        Err(FissionThermalError::Kinetics(KineticsError::InvalidInput(
            _
        )))
    ));
    assert_eq!(initial.energy_j(), 2400.0);
    assert_eq!(
        heating.source(CellId(1)).unwrap().inventory.populations(),
        &[100.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    );
    assert_eq!(heating.audit().released, EnergyChannels::default());
    let heating = BurnupHeating::new(
        &initial,
        vec![product_source(None, true, true)],
        burnup_budget(),
    )
    .unwrap();
    let mut polls = 0;
    let result = heating.advance_to(&initial, 1.0, controls(0.1), || {
        polls += 1;
        if polls == 20 {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    });
    assert!(matches!(
        result,
        Err(FissionThermalError::Thermal(DepositionError::Cancelled))
    ));
    assert_eq!(heating.audit().fissions, 0.0);
    assert_eq!(
        heating.source(CellId(1)).unwrap().inventory.populations()[0],
        100.0
    );
}

#[test]
fn composition_fits_require_a_tracked_identity_and_bounded_resources() {
    let initial = thermal(
        vec![cell(1, 1, 2.0, 1200.0, 0.0)],
        &[(MaterialId(1), &HOT)],
        vec![],
    );
    let mut source = product_source(None, false, false);
    source.composition_feedback.push(
        CompositionFeedback::new(
            "absent-nuclide-v1".into(),
            NuclideId::new(54, 136, 1).unwrap(),
            0.0,
            0.0,
            100.0,
            -0.01,
        )
        .unwrap(),
    );
    assert!(matches!(
        BurnupHeating::new(&initial, vec![source], burnup_budget()),
        Err(FissionThermalError::Invalid(_))
    ));
    let source = poison_source(true);
    let mut denied = burnup_kinetics_budget();
    denied.irradiation.max_composition_feedback_terms = 0;
    assert!(matches!(
        KineticHeating::new(&initial, vec![source], denied),
        Err(FissionThermalError::Kinetics(KineticsError::Budget(_)))
    ));
    assert!(
        CompositionFeedback::new(
            "overflow-v1".into(),
            NuclideId::new(54, 135, 0).unwrap(),
            0.0,
            0.0,
            1e308,
            2.0
        )
        .is_err()
    );
}
