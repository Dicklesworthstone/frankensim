//! Local diffuse-interface material states retain their original sources.

use super::*;
use fs_material::liquid::solution::gradient::BinaryLiquidSquareGradientLaw;
use fs_qty::AmountConcentration;
use fs_thermochem::solution::gradient::BinarySquareGradientV1;
use fs_thermochem::solution::gradient::interface::BinaryPlanarInterfaceOptionsV1;

fn field(law: &BinaryLiquidSolutionLaw<'_>) -> BinarySquareGradientV1 {
    BinarySquareGradientV1::new(
        law.thermodynamics().clone(),
        Temperature::new(300.0),
        AmountConcentration::new(10_000.0),
        1e-6,
    )
    .unwrap()
}

#[test]
fn g1_homogeneous_field_cell_retains_spinodal_composition_and_actual_sources() {
    // Explicit equal partial molar volumes: rho=c*(x*M_A+(1-x)*M_B).
    let card = card_with_density(CompositionBasis::MoleFraction, 1.0, [1000.0, 590.0, 180.0]);
    let liquid = law(&card);
    let adapter = BinaryLiquidSquareGradientLaw::try_new(&liquid, field(&liquid), 1e-12).unwrap();
    for amounts in [[0.0, 5.0], [2.0, 3.0], [2.5, 2.5], [4.0, 1.0], [5.0, 0.0]] {
        let local = adapter.resolve_cell(amounts, || false).unwrap();
        let phase = local.phase();
        assert_eq!(phase.component_moles(), amounts);
        near(phase.volume_m3(), 5.0 / 10_000.0);
        near(phase.molar_concentration_mol_m3(), 10_000.0);
        near(phase.mass_kg(), amounts[0] * 0.018 + amounts[1] * 0.10);
        assert_eq!(
            phase.liquid().material().state_coordinate("x_a"),
            Some(amounts[0] / 5.0)
        );
        assert_eq!(
            phase.thermal_transport_state().unwrap().temperature_k(),
            300.0
        );
        for property in phase.liquid().material().properties() {
            assert_eq!(
                property.answer().receipt.selected,
                liquid
                    .reference_state()
                    .material()
                    .property(property.requirement().name())
                    .unwrap()
                    .answer()
                    .receipt
                    .selected
            );
        }
    }
    let local = adapter.resolve_cell([2.5, 2.5], || false).unwrap();
    assert!(local.thermodynamics().thermodynamic_factor() < 0.0);
    assert_eq!(local.phase().mole_fractions(), [0.5, 0.5]);
    let bulk = liquid
        .equilibrate(
            Temperature::new(300.0),
            [2.5, 2.5],
            BinaryCoexistenceOptionsV1::default(),
            || false,
        )
        .unwrap();
    assert_eq!(bulk.phases().len(), 2);
    assert_ne!(
        bulk.phases()[0].mole_fractions(),
        local.phase().mole_fractions()
    );
}

#[test]
fn g0_field_adapter_refuses_incompatible_density_and_another_pair() {
    let incompatible = card(CompositionBasis::MoleFraction, 1.0);
    let liquid = law(&incompatible);
    let adapter = BinaryLiquidSquareGradientLaw::try_new(&liquid, field(&liquid), 1e-12).unwrap();
    assert!(adapter.resolve_cell([2.0, 3.0], || false).is_err());
    let other_pair = BinaryRegularSolutionV1::new(
        [SpeciesId::new("A").unwrap(), SpeciesId::new("B").unwrap()],
        0.0,
        1000.0,
        Temperature::new(290.0),
        Temperature::new(600.0),
    )
    .unwrap();
    let other = BinarySquareGradientV1::new(
        other_pair,
        Temperature::new(300.0),
        AmountConcentration::new(10_000.0),
        1e-6,
    )
    .unwrap();
    assert!(BinaryLiquidSquareGradientLaw::try_new(&liquid, other, 1e-12).is_err());
    assert!(BinaryLiquidSquareGradientLaw::try_new(&liquid, field(&liquid), 0.01).is_err());
}

#[test]
fn g4_field_source_failure_and_cancellation_publish_no_partial_cell() {
    let card = card_with_density(CompositionBasis::MoleFraction, 0.8, [1000.0, 590.0, 180.0]);
    let liquid = law(&card);
    let adapter = BinaryLiquidSquareGradientLaw::try_new(&liquid, field(&liquid), 1e-12).unwrap();
    let input = [2.0, 3.0];
    let mut calls = 0;
    let accepted = adapter
        .resolve_cell(input, || {
            calls += 1;
            false
        })
        .unwrap();
    for stop in 1..=calls {
        let mut count = 0;
        assert!(
            adapter
                .resolve_cell(input, || {
                    count += 1;
                    count == stop
                })
                .is_err()
        );
    }
    assert_eq!(accepted.phase().component_moles(), input);
    assert!(adapter.resolve_cell([4.5, 0.5], || false).is_err());
    assert!(adapter.resolve_cell([1.0, 1e-20], || false).is_err());
    assert!(adapter.resolve_cell([0.0, 0.0], || false).is_err());
}

#[test]
fn g1_planar_interface_retains_actual_inventory_and_coexisting_source_properties() {
    let card = card_with_density(CompositionBasis::MoleFraction, 1.0, [1000.0, 590.0, 180.0]);
    let liquid = law(&card);
    let adapter = BinaryLiquidSquareGradientLaw::try_new(&liquid, field(&liquid), 1e-12).unwrap();
    let options = BinaryPlanarInterfaceOptionsV1::default();
    let resolved = adapter
        .equilibrium_interface([4.0, 6.0], options, || false)
        .unwrap()
        .unwrap();
    assert_eq!(resolved.equilibrium().total_component_moles(), [4.0, 6.0]);
    assert_eq!(resolved.equilibrium().phases().len(), 2);
    assert_eq!(resolved.interface().model(), adapter.gradient_model());
    assert_eq!(
        resolved.interface(),
        &adapter
            .gradient_model()
            .equilibrium_interface(options, || false)
            .unwrap()
            .unwrap()
    );
    assert!(resolved.interface().surface_tension_n_m() > 0.0);
    assert!(resolved.interface().central_slope_width_m() > 0.0);
    for phase in resolved.equilibrium().phases() {
        near(phase.molar_concentration_mol_m3(), 10_000.0);
        assert_eq!(
            phase.liquid().material().card_identity(),
            liquid.reference_state().material().card_identity()
        );
    }
    near(resolved.equilibrium().total_volume_m3(), 0.001);
    near(
        resolved.equilibrium().total_mass_kg(),
        4.0 * 0.018 + 6.0 * 0.10,
    );
    assert!(
        adapter
            .equilibrium_interface([0.0, 10.0], options, || false)
            .unwrap()
            .is_none()
    );
}

#[test]
fn g4_planar_interface_never_bypasses_source_domain_or_volume_admission() {
    for (upper, density) in [
        (0.8, [1000.0, 590.0, 180.0]),
        (1.0, [1000.0, 1100.0, 800.0]),
    ] {
        let card = card_with_density(CompositionBasis::MoleFraction, upper, density);
        let liquid = law(&card);
        let adapter =
            BinaryLiquidSquareGradientLaw::try_new(&liquid, field(&liquid), 1e-12).unwrap();
        assert!(
            adapter
                .equilibrium_interface(
                    [4.0, 6.0],
                    BinaryPlanarInterfaceOptionsV1::default(),
                    || false
                )
                .is_err()
        );
    }
    let card = card_with_density(CompositionBasis::MoleFraction, 1.0, [1000.0, 590.0, 180.0]);
    let liquid = law(&card);
    let adapter = BinaryLiquidSquareGradientLaw::try_new(&liquid, field(&liquid), 1e-12).unwrap();
    assert!(
        adapter
            .equilibrium_interface(
                [4.0, 6.0],
                BinaryPlanarInterfaceOptionsV1::default(),
                || true
            )
            .is_err()
    );
}
