//! Native transient pore flow, finite site storage and finite skeleton heat.
use super::*;
use crate::porous::permeability::{
    DarcyGasBarrier, IntrinsicPermeability, PorousGasLayer, PorousGasSupport,
};
use crate::thermochemical_gas::energy::parcel::transfer::species::diffusion::dusty_gas::{
    DustyGasFaceSolve, DustyGasMedium, DustyGasPores,
};
use fs_thermochem::{Composition, CompositionBasis};
use transport::pore_domain::{
    PoreCellGeometry, PoreDomain, PoreDomainError, PoreDomainSolve, PoreFace,
    adsorption::AdsorbingPoreDomain,
};

fn medium() -> DustyGasMedium {
    let layer = PorousGasLayer::try_new(
        IntrinsicPermeability::isotropic(2e-14).unwrap(),
        0.01,
        0.0,
        "manufactured finite adsorbing pore geometry",
    )
    .unwrap();
    let barrier = DarcyGasBarrier::try_new(
        vec![layer],
        [1.0, 0.0, 0.0],
        PorousGasSupport {
            temperature_k: (200.0, 800.0),
            pressure_pa: (1000.0, 1e7),
            maximum_superficial_speed_m_s: 100.0,
        },
    )
    .unwrap();
    DustyGasMedium::try_new(
        barrier,
        DustyGasPores::try_new(
            1e-7,
            0.4,
            2.0,
            "manufactured rigid pore source; no empirical validation",
        )
        .unwrap(),
    )
    .unwrap()
}
fn faces<'a>(medium: &'a DustyGasMedium, context: &'a QueryPoint, area: f64) -> [PoreFace<'a>; 2] {
    [PoreFace {
        medium,
        area_m2: area,
        pairs: &[],
        gas_context: context,
        solve: DustyGasFaceSolve {
            maximum_temperature_difference_k: 1.0,
            maximum_relative_pressure_difference: 0.25,
            gradient_tolerance_per_m: 1e-8,
        },
    }; 2]
}
fn window(dt: f64) -> PoreDomainSolve {
    PoreDomainSolve {
        maximum_step_s: dt,
        minimum_step_s: 1e-9,
        maximum_substeps: 1024,
        maximum_attempts: 4096,
        maximum_throat_volume_fraction: 0.02,
        transport: controls(),
    }
}
fn domain(
    medium: &DustyGasMedium,
    left: ClosedGasParcel,
    gas: ClosedGasParcel,
    right: ClosedGasParcel,
) -> PoreDomain {
    // The original fixture's 0.001 m³ gas equals phi * A * L exactly.
    PoreDomain::try_new(
        0.25,
        vec![PoreCellGeometry::try_new(medium.clone()).unwrap()],
        left,
        vec![gas],
        right,
    )
    .unwrap()
}
fn pure_fixture(medium: &DustyGasMedium) -> AdsorbingPoreDomain {
    let (fixture, _, catalog) = fixture([0.0, 0.0]);
    let solid = CompetitiveSolid::ordered(
        Arc::new(CompetitiveAdsorbent::try_new(vec![fixture.model.channels[0].clone()]).unwrap()),
        0.1,
        vec![0.0],
        300.0,
    )
    .unwrap();
    let pure = catalog
        .with_composition(Composition::new(CompositionBasis::MoleFraction, vec![1.0, 0.0]).unwrap())
        .unwrap();
    let gas = |pressure| ClosedGasParcel::new(pure.clone(), point(300.0, pressure), 0.001).unwrap();
    AdsorbingPoreDomain::try_new(
        domain(medium, gas(105000.0), gas(100000.0), gas(95000.0)),
        BTreeMap::from([(0, solid)]),
        catalog,
    )
    .unwrap()
}

#[test]
fn closed_pore_reuses_exact_site_balance_and_original_finite_heat() {
    let (solid, gas, catalog) = fixture([0.0, 0.0]);
    let medium = medium();
    let domain = domain(&medium, gas.clone(), gas.clone(), gas.clone());
    let original = AdsorbingPoreDomain::try_new(
        domain,
        BTreeMap::from([(0, solid.clone())]),
        catalog.clone(),
    )
    .unwrap();
    let saved = original.clone();
    let context = point(300.0, 100000.0);
    let step = 0.001;
    let result = original
        .step_candidate(
            step,
            &faces(&medium, &context, 0.0),
            &[pair_controls().gas.gas; 3],
            window(step),
            &BTreeMap::from([(0, pair_controls())]),
            || false,
        )
        .unwrap();
    let actual = &result.domain().solids()[&0];
    let expected_moles = [
        0.05 * step * 0.4 / (1.0 + step * 2.4),
        0.05 * step * 1.6 / (1.0 + step * 2.4),
    ];
    for (i, law) in solid.model.channels.iter().enumerate() {
        near(
            actual.adsorbed_moles(law.species()).unwrap(),
            expected_moles[i],
            2e-16,
        );
    }
    let heat = 6000.0 * expected_moles[0] + 9000.0 * expected_moles[1];
    let capacity = 0.1 * 500.0 + R * (3.5 * expected_moles[0] + 2.5 * expected_moles[1]);
    near(actual.temperature_k(), 300.0 + heat / capacity, 2e-10);
    let pair = super::super::super::advance_candidate(
        &solid,
        &gas,
        &catalog,
        step,
        pair_controls(),
        || false,
    )
    .unwrap();
    assert_eq!(actual, pair.solid());
    let pore_gas = result.domain().pore_domain().pore_gas(0).unwrap();
    near(pore_gas.mass_kg(), pair.gas().mass_kg(), 1e-13);
    near(
        pore_gas.state().state().temperature,
        pair.gas().state().state().temperature,
        2e-7,
    );
    assert_eq!(result.domain().pore_domain().left_reservoir(), &gas);
    assert_eq!(result.domain().pore_domain().right_reservoir(), &gas);
    assert!(result.step().pore.faces.is_empty());
    near(
        result.step().pore.pore_mass_change_kg + result.step().adsorbed_mass_change_kg,
        0.0,
        1e-13,
    );
    near(
        result.step().pore.total_energy_residual_j,
        0.0,
        controls().exchange.energy_tolerance_j,
    );
    assert!(result.step().pore.entropy_change_j_k >= -controls().entropy_tolerance_j_k);
    assert_eq!(original, saved);
}

#[test]
fn evolving_loading_and_heat_change_subsequent_real_porous_throughput() {
    let medium = medium();
    let original = pure_fixture(&medium);
    let context = point(300.0, 100000.0);
    let faces = faces(&medium, &context, 0.001);
    let dt = 0.00025;
    let duration = 8.0 * dt;
    let endpoints = [pair_controls().gas.gas; 3];
    let local = BTreeMap::from([(0, pair_controls())]);
    let inert = original
        .pore_domain()
        .step_candidate(duration, &faces, &endpoints, window(dt), || false)
        .unwrap();
    let result = original
        .step_candidate(duration, &faces, &endpoints, window(dt), &local, || false)
        .unwrap();
    let gas = result.domain().pore_domain().pore_gas(0).unwrap();
    assert!(result.domain().solids()[&0].temperature_k() > 300.0);
    assert!(
        gas.state().state().temperature
            < inert
                .domain()
                .pore_gas(0)
                .unwrap()
                .state()
                .state()
                .temperature
    );
    assert!(
        gas.state().state().pressure < inert.domain().pore_gas(0).unwrap().state().state().pressure
    );
    let flow = &result.step().pore.last_substep.transport.dusty_gas;
    let inert_flow = &inert.step().last_substep.dusty_gas;
    assert!(flow[0].flux.net_molar_flux_mol_m2_s > inert_flow[0].flux.net_molar_flux_mol_m2_s);
    assert!(flow[1].flux.net_molar_flux_mol_m2_s < inert_flow[1].flux.net_molar_flux_mol_m2_s);
    assert!(
        result
            .step()
            .pore
            .faces
            .iter()
            .any(|f| f.mass_left_to_right_kg > 0.0)
    );
    assert!(result.step().adsorbed_mass_change_kg > 0.0);
    assert!(
        result
            .step()
            .adsorption
            .iter()
            .all(|r| r.gross_adsorption_mol > 0.0 && r.gross_desorption_mol > 0.0)
    );
    assert_eq!(
        result.domain().pore_domain().cells(),
        original.pore_domain().cells()
    );
    assert_eq!(
        result.domain().solids()[&0].model(),
        original.solids()[&0].model()
    );
    near(
        result.step().pore.maximum_species_mass_residual_kg,
        0.0,
        controls().exchange.mass_tolerance_kg,
    );
    near(
        result.step().pore.total_energy_residual_j,
        0.0,
        controls().exchange.energy_tolerance_j,
    );
    near(result.domain().pore_domain().time_s(), duration, 1e-16);

    // Continuing windows preserve the same accepted site and temperature history.
    let mut continued = original.clone();
    for _ in 0..8 {
        continued
            .advance(dt, &faces, &endpoints, window(dt), &local, || false)
            .unwrap();
    }
    near(
        continued.solids()[&0].temperature_k(),
        result.domain().solids()[&0].temperature_k(),
        2e-10,
    );
    for (id, amount) in continued.solids()[&0].loadings() {
        near(
            amount,
            result.domain().solids()[&0].adsorbed_moles(id).unwrap(),
            2e-15,
        );
    }
    for (&id, gas) in continued.pore_domain().parcels() {
        let other = &result.domain().pore_domain().parcels()[&id];
        near(gas.mass_kg(), other.mass_kg(), 1e-13);
        near(
            gas.state().state().temperature,
            other.state().state().temperature,
            2e-7,
        );
    }
}

#[test]
fn typed_site_inventory_limits_refine_the_same_pore_window() {
    let medium = medium();
    let original = pure_fixture(&medium);
    let context = point(300.0, 100000.0);
    let mut local = pair_controls();
    local.maximum_gas_gross_fraction = 1e-4;
    let result = original
        .step_candidate(
            0.001,
            &faces(&medium, &context, 0.0),
            &[pair_controls().gas.gas; 3],
            window(0.001),
            &BTreeMap::from([(0, local)]),
            || false,
        )
        .unwrap();
    assert!(result.step().pore.attempts > result.step().pore.accepted_substeps);
    assert!(result.step().pore.accepted_substeps > 1);
    assert!(result.step().pore.maximum_accepted_step_s < 0.001);
    near(result.domain().pore_domain().time_s(), 0.001, 1e-16);
    assert!(result.step().adsorbed_mass_change_kg > 0.0);
}

#[test]
fn late_window_failure_or_cancellation_preserves_all_gas_sites_and_clock() {
    let medium = medium();
    let mut original = pure_fixture(&medium);
    let saved = original.clone();
    let context = point(300.0, 100000.0);
    let faces = faces(&medium, &context, 0.0);
    let endpoints = [pair_controls().gas.gas; 3];
    let local = BTreeMap::from([(0, pair_controls())]);
    let mut work = window(0.00025);
    work.maximum_substeps = 1;
    assert!(matches!(
        original.advance(0.0005, &faces, &endpoints, work, &local, || false),
        Err(PoreDomainError::WorkLimit {
            accepted_substeps: 1,
            ..
        })
    ));
    assert_eq!(original, saved);
    let mut polls = 0;
    assert!(matches!(
        original.advance(0.0005, &faces, &endpoints, window(0.00025), &local, || {
            polls += 1;
            polls > 20
        }),
        Err(PoreDomainError::Cancelled)
    ));
    assert_eq!(original, saved);
    assert!(
        original
            .step_candidate(
                0.00025,
                &faces,
                &endpoints,
                window(0.00025),
                &BTreeMap::new(),
                || false
            )
            .is_err()
    );
    let mut inconsistent = local;
    inconsistent
        .get_mut(&0)
        .unwrap()
        .gas
        .gas
        .maximum_temperature_k = 700.0;
    assert!(
        original
            .step_candidate(
                0.00025,
                &faces,
                &endpoints,
                window(0.00025),
                &inconsistent,
                || false
            )
            .is_err()
    );
    assert_eq!(original, saved);
}
