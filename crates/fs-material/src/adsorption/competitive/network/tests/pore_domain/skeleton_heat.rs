//! Native pore-history feedback and atomic late contact-source refusal.
use super::super::skeleton_heat::{
    contact, contact_card, contact_point, fixture_with_heat, heat_link, layer,
};
use super::*;

fn setup(medium: &DustyGasMedium) -> (MaterialCard, [CompetitiveSolid; 2], AdsorbingPoreDomain) {
    let (card, solids, _, catalog) = fixture_with_heat(200.0, 0.4);
    let pure = catalog
        .with_composition(Composition::new(CompositionBasis::MoleFraction, vec![1.0, 0.0]).unwrap())
        .unwrap();
    let gas = |pressure| ClosedGasParcel::new(pure.clone(), point(300.0, pressure), 0.001).unwrap();
    let domain = PoreDomain::try_new(
        0.25,
        vec![PoreCellGeometry::try_new(medium.clone()).unwrap(); 2],
        gas(105000.0),
        vec![gas(100000.0), gas(99000.0)],
        gas(95000.0),
    )
    .unwrap();
    let domain = AdsorbingPoreDomain::try_new(
        domain,
        BTreeMap::from([(0, solids[0].clone()), (1, solids[1].clone())]),
        catalog,
    )
    .unwrap();
    (card, solids, domain)
}

#[test]
fn skeleton_heat_changes_later_uptake_and_actual_porous_flow() {
    let medium = medium();
    let (card, retained, original) = setup(&medium);
    let saved = original.clone();
    let context = point(300.0, 100000.0);
    let layers = retained
        .each_ref()
        .map(|solid| layer(solid, &card, &context, [200.0, 800.0]));
    let links = [heat_link(layers.each_ref(), None)];
    let faces = [faces(&medium, &context, 0.001)[0]; 3];
    let endpoints = [pair_controls().gas.gas; 4];
    let local = BTreeMap::from([(0, pair_controls()), (1, pair_controls())]);
    let dt = 0.00025;
    let duration = 8.0 * dt;
    let first = original
        .step_candidate_with_conduction(dt, &faces, &endpoints, window(dt), &local, &links, || {
            false
        })
        .unwrap();
    let isolated_first = original
        .step_candidate(dt, &faces, &endpoints, window(dt), &local, || false)
        .unwrap();
    // Kinetics/flow sample the SAME initial snapshot. Conduction changes this
    // interval's finite T, and reaches gas uptake/transport in later intervals.
    for cell in 0..2 {
        let species = SpeciesId::new("A").unwrap();
        near(
            first.domain().solids()[&cell]
                .adsorbed_moles(&species)
                .unwrap(),
            isolated_first.domain().solids()[&cell]
                .adsorbed_moles(&species)
                .unwrap(),
            1e-16,
        );
        near(
            first
                .domain()
                .pore_domain()
                .pore_gas(cell)
                .unwrap()
                .mass_kg(),
            isolated_first
                .domain()
                .pore_domain()
                .pore_gas(cell)
                .unwrap()
                .mass_kg(),
            1e-13,
        );
    }
    assert!(
        first.domain().solids()[&0].temperature_k()
            < isolated_first.domain().solids()[&0].temperature_k()
    );
    assert!(
        first.domain().solids()[&1].temperature_k()
            > isolated_first.domain().solids()[&1].temperature_k()
    );
    let result = original
        .step_candidate_with_conduction(
            duration,
            &faces,
            &endpoints,
            window(dt),
            &local,
            &links,
            || false,
        )
        .unwrap();
    let isolated = original
        .step_candidate(duration, &faces, &endpoints, window(dt), &local, || false)
        .unwrap();
    let species = SpeciesId::new("A").unwrap();
    assert!(
        result.domain().solids()[&0]
            .adsorbed_moles(&species)
            .unwrap()
            > isolated.domain().solids()[&0]
                .adsorbed_moles(&species)
                .unwrap()
    );
    assert!(
        result.domain().solids()[&1]
            .adsorbed_moles(&species)
            .unwrap()
            < isolated.domain().solids()[&1]
                .adsorbed_moles(&species)
                .unwrap()
    );
    let rates = &result.step().pore.last_substep.transport.dusty_gas;
    let old_rates = &isolated.step().pore.last_substep.transport.dusty_gas;
    assert!(
        rates
            .iter()
            .zip(old_rates)
            .any(|(a, b)| a.flux.net_molar_flux_mol_m2_s != b.flux.net_molar_flux_mol_m2_s)
    );
    assert_eq!(result.step().conduction.len(), 1);
    let heat = &result.step().conduction[0];
    assert_eq!(heat.cells, [0, 1]);
    assert_eq!(heat.id, 41);
    assert!(heat.gross_a_to_b_j > 0.0);
    assert_eq!(heat.gross_b_to_a_j, 0.0);
    assert!(heat.gross_a_to_b_j > result.step().pore.last_substep.conduction[0].energy_a_to_b_j);
    assert_eq!(
        result.step().pore.last_substep.conduction[0].caloric_curve_identities,
        layers
            .each_ref()
            .map(|layer| layer.caloric_curve().identity())
    );
    near(
        result.step().pore.total_energy_residual_j,
        0.0,
        controls().exchange.energy_tolerance_j,
    );
    near(
        result.step().pore.maximum_species_mass_residual_kg,
        0.0,
        controls().exchange.mass_tolerance_kg,
    );
    assert!(result.step().pore.entropy_change_j_k >= -controls().entropy_tolerance_j_k);
    assert_eq!(original, saved);

    // Continuing accepted windows preserve the same finite thermal history and
    // account for every applied heat increment, not only the final rate sample.
    let mut continued = original.clone();
    let mut total = 0.0;
    for _ in 0..8 {
        let step = continued
            .advance_with_conduction(dt, &faces, &endpoints, window(dt), &local, &links, || false)
            .unwrap();
        total += step.conduction[0].gross_a_to_b_j;
    }
    near(total, heat.gross_a_to_b_j, 1e-11);
    for cell in 0..2 {
        near(
            continued.solids()[&cell].temperature_k(),
            result.domain().solids()[&cell].temperature_k(),
            2e-10,
        );
        near(
            continued.solids()[&cell].adsorbed_moles(&species).unwrap(),
            result.domain().solids()[&cell]
                .adsorbed_moles(&species)
                .unwrap(),
            2e-15,
        );
        near(
            continued
                .pore_domain()
                .pore_gas(cell)
                .unwrap()
                .state()
                .state()
                .pressure,
            result
                .domain()
                .pore_domain()
                .pore_gas(cell)
                .unwrap()
                .state()
                .state()
                .pressure,
            1e-6,
        );
    }
    near(continued.pore_domain().time_s(), duration, 1e-16);
}

#[test]
fn late_ordered_contact_refusal_preserves_all_accepted_input_histories() {
    let medium = medium();
    let (card, retained, mut original) = setup(&medium);
    let saved = original.clone();
    let context = point(300.0, 100000.0);
    let layers = retained
        .each_ref()
        .map(|solid| layer(solid, &card, &context, [200.0, 800.0]));
    // Initial cold skeleton is 299 K. The first coupled interval reaches about
    // 299.0019 K; the next leaves this original contact's declared support.
    let source = contact_card(&card, 0.0001, 299.003);
    let contact_context = contact_point();
    let contact = contact(&retained, &source, &card, &contact_context);
    let links = [heat_link(layers.each_ref(), Some(&contact))];
    let faces = [faces(&medium, &context, 0.001)[0]; 3];
    let endpoints = [pair_controls().gas.gas; 4];
    let local = BTreeMap::from([(0, pair_controls()), (1, pair_controls())]);
    let dt = 0.00025;
    original
        .step_candidate_with_conduction(dt, &faces, &endpoints, window(dt), &local, &links, || {
            false
        })
        .unwrap();
    let error = original
        .advance_with_conduction(
            8.0 * dt,
            &faces,
            &endpoints,
            window(dt),
            &local,
            &links,
            || false,
        )
        .unwrap_err();
    assert!(matches!(error,PoreDomainError::Adsorption { elapsed_s,
        error:AdsorptionNetworkError::Conduction { id:41,.. },.. } if elapsed_s>0.0));
    assert_eq!(original, saved);
}
