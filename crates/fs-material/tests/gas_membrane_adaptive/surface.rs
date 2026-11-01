//! Manufactured surface-transfer laws exercise the original film/gas owners.
use super::*;
use fs_matdb::{InterfaceSystemCard, SurfaceSpec, SystemContext};
use fs_material::thermochemical_gas::energy::parcel::transfer::species::network::transport::membrane::{
    MembraneInterfaceMaterial, MOLAR_INTERFACE_PERMEANCE_PROPERTY, MOLAR_PERMEANCE_DIMS,
    transient::ResolvedMembraneDiffusion,
};

fn interface(
    left: &MaterialCard,
    right: &MaterialCard,
    permeance: f64,
    pressures: [f64; 2],
) -> InterfaceSystemCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(claim(
            MOLAR_INTERFACE_PERMEANCE_PROPERTY,
            permeance,
            MOLAR_PERMEANCE_DIMS,
            pressures,
        ))
        .unwrap();
    InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: left.id().clone(),
            texture_frame: "left trace".into(),
        },
        SurfaceSpec {
            material: right.id().clone(),
            texture_frame: "right trace".into(),
        },
        SystemContext {
            medium: "explicit independent species transfer".into(),
            third_body: None,
            environment: "manufactured isothermal surface".into(),
            history: "fixed".into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn surface<'a>(
    name: &str,
    card: &'a InterfaceSystemCard,
    context: &'a QueryPoint,
) -> MembraneInterfaceMaterial<'a> {
    MembraneInterfaceMaterial {
        species: id(name),
        interface: card,
        context,
        selection: MaterialPropertySelection::SingleClaimOnly,
    }
}
fn narrow_polymer() -> MaterialCard {
    card(
        "surface-limited polymer",
        &[
            (MEMBRANE_DIFFUSIVITY_PROPERTY, D, MEMBRANE_DIFFUSIVITY_DIMS),
            (
                MOLAR_HENRY_SOLUBILITY_PROPERTY,
                S,
                MOLAR_HENRY_SOLUBILITY_DIMS,
            ),
        ],
        [99000.0, 101000.0],
    )
}

#[test]
fn source_backed_surface_limited_steady_flux_and_finite_relaxation_match_independent_laws() {
    let polymer = narrow_polymer();
    let bank = model(&[1.0], 0.0);
    let gas_material = &bank.transport_materials()[0].material;
    let half = 2.0 * D * S / THICKNESS;
    let h = 0.01 * half;
    let left = interface(gas_material, &polymer, h, [0.0, 2e5]);
    let right = interface(&polymer, gas_material, h, [0.0, 2e5]);
    let context = point(1e5);
    let sources = [channel("A", &polymer, &context, 1, 1e5)];
    let surfaces = [
        Some(surface("A", &left, &context)),
        Some(surface("A", &right, &context)),
    ];
    // The bulk law is qualified only near 100 kPa, while the gas boundaries
    // are 140/60 kPa. Actual polymer-side surface traces must govern D and S.
    let fixed = ResolvedMembraneDiffusion::new_with_surfaces(
        T,
        AREA,
        &sources[0].layers,
        &[],
        &surfaces,
        &[1e5],
        [1.4e5, 0.6e5],
        || false,
    )
    .unwrap();
    let sample = fixed
        .state()
        .sample_fluxes([1.4e5, 0.6e5], || false)
        .unwrap();
    let flux = 80000.0 / (2.0 / h + THICKNESS / (D * S));
    near(sample.face_molar_fluxes_mol_m2_s[0], flux, flux * 2e-14);
    near(sample.face_molar_fluxes_mol_m2_s[1], flux, flux * 2e-14);
    near(sample.face_partial_pressures_pa[0][0].unwrap(), 1.4e5, 1e-9);
    near(
        sample.face_partial_pressures_pa[0][1].unwrap(),
        1.4e5 - flux / h,
        1e-9,
    );
    near(
        sample.face_partial_pressures_pa[1][0].unwrap(),
        0.6e5 + flux / h,
        1e-9,
    );
    near(sample.face_partial_pressures_pa[1][1].unwrap(), 0.6e5, 1e-9);
    for receipt in &fixed.source_support().layers[0] {
        assert!((99000.0..=101000.0).contains(&receipt.state_coordinate("pressure").unwrap()));
    }
    assert_eq!(
        fixed.source_support().surfaces[0].as_ref().unwrap()[1].state_coordinate("pressure"),
        Some(1.4e5)
    );
    let mut native = fs_material::thermochemical_gas::energy::parcel::transfer::species::network::transport::membrane::transient::MembraneDiffusionWorkspace::new(1).unwrap();
    let steady = fixed
        .step_candidate(
            [1.4e5, 0.6e5],
            100.0,
            native_controls().membrane,
            &mut native,
            || false,
        )
        .unwrap();
    near(steady.state().state().partial_pressures_pa()[0], 1e5, 1e-8);
    near(
        steady.step().left_reservoir_to_membrane_mol,
        AREA * flux * 100.0,
        1e-18,
    );

    let volume = equal_capacity_volume();
    let gases = [
        gas(&[1.0], 0.0, 1.4e5, volume),
        gas(&[1.0], 0.0, 0.6e5, volume),
    ];
    let state = IsothermalGasMembraneMixture::new_with_surfaces(
        T,
        AREA,
        &sources,
        &[surfaces],
        gases.clone(),
        native_controls(),
        || false,
    )
    .unwrap();
    // Equal gas/film capacities give an exact antisymmetric mode. Surface
    // resistance changes its rate to Gface/Cgas, not the bulk-only 2D/L².
    let face_permeance = 1.0 / (1.0 / h + 1.0 / half);
    let rate = AREA * face_permeance / (volume / (R * T));
    let duration = 1.0 / rate;
    let mut work = state.adaptive_workspace().unwrap();
    let result = state
        .advance_candidate(duration, controls(duration, 0.005), &mut work, || false)
        .unwrap();
    let expected = 40000.0 * (-1.0_f64).exp();
    near(
        result.state().gases()[0].state().state().pressure,
        1e5 + expected,
        120.0,
    );
    near(
        result.state().gases()[1].state().state().pressure,
        1e5 - expected,
        120.0,
    );
    assert!(rate < 0.02 * (2.0 * D / THICKNESS.powi(2)));
    assert_eq!(
        result.report().source_relative_gas_energy_residual_j,
        [0.0; 2]
    );
    near(
        result.report().species[0].closed_species_residual_mol,
        0.0,
        1e-17,
    );
    assert_eq!(state.gases(), &gases);
}

#[test]
fn surface_selectivity_changes_finite_mixture_composition_with_distinct_actual_traces() {
    let polymer = narrow_polymer();
    let bank = model(&[0.5, 0.5], 0.0);
    let materials = bank.transport_materials();
    let ha = 4e-13;
    let hb = 4e-14;
    let left_a = interface(&materials[0].material, &polymer, ha, [0.0, 2e5]);
    let right_a = interface(&polymer, &materials[0].material, ha, [0.0, 2e5]);
    let left_b = interface(&materials[1].material, &polymer, hb, [0.0, 2e5]);
    let right_b = interface(&polymer, &materials[1].material, hb, [0.0, 2e5]);
    let context = point(1e5);
    let mut expected = None;
    for reverse in [false, true] {
        let mut sources = vec![
            channel("A", &polymer, &context, 1, 5e4),
            channel("B", &polymer, &context, 1, 5e4),
        ];
        let mut surfaces = vec![
            [
                Some(surface("A", &left_a, &context)),
                Some(surface("A", &right_a, &context)),
            ],
            [
                Some(surface("B", &left_b, &context)),
                Some(surface("B", &right_b, &context)),
            ],
        ];
        if reverse {
            sources.reverse();
            surfaces.reverse();
        }
        let volume = equal_capacity_volume();
        let state = IsothermalGasMembraneMixture::new_with_surfaces(
            T,
            AREA,
            &sources,
            &surfaces,
            [
                gas(&[0.9, 0.1], 0.0, 1e5, volume),
                gas(&[0.1, 0.9], 0.0, 1e5, volume),
            ],
            native_controls(),
            || false,
        )
        .unwrap();
        let mut work = state.adaptive_workspace().unwrap();
        let result = state
            .advance_candidate(1000.0, controls(1000.0, 0.01), &mut work, || false)
            .unwrap();
        let actual = result.state();
        assert!(actual.gases()[0].state().state().pressure < 99000.0);
        assert!(actual.gases()[1].state().state().pressure > 101000.0);
        assert!(
            actual.gases()[1]
                .species_partial_pressure_pa(&id("A"))
                .unwrap()
                / actual.gases()[1].state().state().pressure
                > 0.1
        );
        for name in ["A", "B"] {
            let support = actual.source_support(&id(name)).unwrap();
            for receipt in &support.layers[0] {
                assert!(
                    (99000.0..=101000.0).contains(&receipt.state_coordinate("pressure").unwrap())
                );
            }
            let contact = support.surfaces[0].as_ref().unwrap();
            assert!(contact[0].state_coordinate("pressure").unwrap() < 99000.0);
            assert!(contact[1].state_coordinate("pressure").unwrap() > 99000.0);
        }
        let rows = &result.report().species;
        assert!(rows[0].amount_into_gas_mol[0] < 0.0);
        assert!(rows[1].amount_into_gas_mol[0] > 0.0);
        assert_eq!(
            result.report().source_relative_gas_energy_residual_j,
            [0.0; 2]
        );
        let values = (
            actual.gases().each_ref().map(|g| g.mass_kg()),
            actual
                .gases()
                .each_ref()
                .map(|g| g.state().state().pressure),
            rows.iter()
                .map(|r| r.amount_into_gas_mol)
                .collect::<Vec<_>>(),
        );
        if let Some(expected) = &expected {
            assert_eq!(&values, expected);
        } else {
            expected = Some(values);
        }
    }
}

#[test]
fn surface_seals_and_original_order_identity_and_two_sided_domains_are_not_inferred() {
    let polymer = membrane("polymer", D, S, 2e5);
    let bank = model(&[1.0], 0.0);
    let gas_material = &bank.transport_materials()[0].material;
    let left = interface(gas_material, &polymer, 0.0, [0.0, 2e5]);
    let context = point(1e5);
    let sources = [channel("A", &polymer, &context, 1, 1e5)];
    let surfaces = [[Some(surface("A", &left, &context)), None]];
    let volume = equal_capacity_volume();
    let gases = [
        gas(&[1.0], 0.0, 1.4e5, volume),
        gas(&[1.0], 0.0, 0.6e5, volume),
    ];
    let state = IsothermalGasMembraneMixture::new_with_surfaces(
        T,
        AREA,
        &sources,
        &surfaces,
        gases.clone(),
        native_controls(),
        || false,
    )
    .unwrap();
    let mut work = state.adaptive_workspace().unwrap();
    let result = state
        .advance_candidate(1.0, controls(1.0, 0.1), &mut work, || false)
        .unwrap();
    assert_eq!(result.state().gases()[0], gases[0]);
    assert_eq!(result.report().species[0].amount_into_gas_mol[0], 0.0);
    assert!(result.report().species[0].amount_into_gas_mol[1] > 0.0);
    assert!(film(result.state(), "A")[0] < 1e5);

    let foreign_bank = model(&[0.5, 0.5], 0.0);
    let foreign = &foreign_bank.transport_materials()[1].material;
    let wrong_gas = interface(foreign, &polymer, 4e-13, [0.0, 2e5]);
    let wrong_order = interface(&polymer, gas_material, 4e-13, [0.0, 2e5]);
    let narrow = interface(gas_material, &polymer, 4e-13, [99000.0, 101000.0]);
    for contact in [&wrong_gas, &wrong_order, &narrow] {
        assert!(
            IsothermalGasMembraneMixture::new_with_surfaces(
                T,
                AREA,
                &sources,
                &[[Some(surface("A", contact, &context)), None]],
                gases.clone(),
                native_controls(),
                || false
            )
            .is_err()
        );
    }
    // With both the surface and adjacent bulk sealed, the intervening trace
    // has no constitutive determination. Do not invent pressure for its source.
    let blocked = membrane("blocked polymer", 0.0, S, 2e5);
    let double_seal = interface(gas_material, &blocked, 0.0, [0.0, 2e5]);
    let source = [channel("A", &blocked, &context, 1, 1e5)];
    assert!(
        IsothermalGasMembraneMixture::new_with_surfaces(
            T,
            AREA,
            &source,
            &[[Some(surface("A", &double_seal, &context)), None]],
            gases.clone(),
            native_controls(),
            || false
        )
        .is_err()
    );
    let saved = state.finite_state(&id("A")).unwrap().clone();
    assert!(matches!(
        state.advance_candidate(1.0, controls(1.0, 0.1), &mut work, || true),
        Err(MembraneMixtureAdvanceError::Cancelled)
    ));
    assert_eq!(state.gases(), &gases);
    assert_eq!(state.finite_state(&id("A")).unwrap(), &saved);
}
