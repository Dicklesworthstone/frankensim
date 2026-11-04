use super::*;
use super::super::tests::{
    AREA, LENGTH, P, POROSITY, RADIUS, T, TORTUOSITY, binary, card, controls, knudsen, near,
    pair, pure, support,
};
use crate::porous::permeability::{DarcyGasBarrier, IntrinsicPermeability, PorousGasLayer};
use crate::thermochemical_gas::tests::point;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;

fn pores(radius: f64, porosity: f64, tortuosity: f64) -> DustyGasPores {
    DustyGasPores::try_new(radius, porosity, tortuosity, "synthetic layered pore geometry").unwrap()
}

fn layered(layers: &[(f64, f64, DustyGasPores)]) -> LayeredDustyGasMedium {
    let porous = layers.iter().map(|(k, length, _)| {
        PorousGasLayer::try_new(
            IntrinsicPermeability::isotropic(*k).unwrap(), *length, 0.0,
            "synthetic layered intrinsic permeability",
        ).unwrap()
    }).collect();
    let barrier = DarcyGasBarrier::try_new(porous, [1.0, 0.0, 0.0], support()).unwrap();
    LayeredDustyGasMedium::try_new(
        barrier, layers.iter().map(|(_, _, p)| p.clone()).collect(),
    ).unwrap()
}

#[test]
fn one_layer_is_the_original_homogeneous_dusty_gas_law() {
    let (left, right) = (binary(0.75, P), binary(0.25, P));
    let d = card(2e-5);
    let context = point(T, P);
    let bindings = [pair("A", "B", &d, &context)];
    let original = super::super::tests::medium(2e-14, LENGTH);
    let stack = layered(&[(
        2e-14, LENGTH, pores(RADIUS, POROSITY, TORTUOSITY),
    )]);
    let a = left.dusty_gas_flux_to(
        &right, &original, AREA, &bindings, &context, controls(), || false,
    ).unwrap();
    let b = left.layered_dusty_gas_flux_to(
        &right, &stack, AREA, &bindings, &context, controls(), || false,
    ).unwrap();
    assert_eq!(a.species.len(), b.species.len());
    for (a, b) in a.species.iter().zip(&b.species) {
        assert_eq!(a.species, b.species);
        near(b.face_mole_fraction, a.face_mole_fraction);
        near(b.effective_knudsen_diffusivity_m2_s, a.knudsen_diffusivity_m2_s);
        near(b.diffusive_molar_flux_mol_m2_s, a.diffusive_molar_flux_mol_m2_s);
        near(b.viscous_molar_flux_mol_m2_s, a.viscous_molar_flux_mol_m2_s);
        near(b.molar_flux_mol_m2_s, a.molar_flux_mol_m2_s);
        near(b.mass_rate_left_to_right_kg_s, a.mass_rate_left_to_right_kg_s);
    }
    near(b.net_molar_flux_mol_m2_s, a.net_molar_flux_mol_m2_s);
    near(
        b.entropy_production_w_m2_k,
        a.local_entropy_production_w_m3_k * LENGTH,
    );
    assert!(b.maximum_integrated_composition_residual
        <= controls().gradient_tolerance_per_m * LENGTH);
}

#[test]
fn splitting_one_material_into_identical_layers_does_not_change_transport() {
    let (left, right) = (binary(0.68, 106_000.0), binary(0.32, 94_000.0));
    let d = card(1.7e-5);
    let context = point(T, P);
    let bindings = [pair("A", "B", &d, &context)];
    let pore = pores(RADIUS, POROSITY, TORTUOSITY);
    let single = layered(&[(3e-14, LENGTH, pore.clone())]);
    let split = layered(&[
        (3e-14, 0.25 * LENGTH, pore.clone()),
        (3e-14, 0.75 * LENGTH, pore),
    ]);
    let a = left.layered_dusty_gas_flux_to(
        &right, &single, AREA, &bindings, &context, controls(), || false,
    ).unwrap();
    let b = left.layered_dusty_gas_flux_to(
        &right, &split, AREA, &bindings, &context, controls(), || false,
    ).unwrap();
    for (a, b) in a.species.iter().zip(&b.species) {
        near(b.effective_knudsen_diffusivity_m2_s, a.effective_knudsen_diffusivity_m2_s);
        near(b.knudsen_resistance_s_m, a.knudsen_resistance_s_m);
        near(b.diffusive_molar_flux_mol_m2_s, a.diffusive_molar_flux_mol_m2_s);
        near(b.viscous_molar_flux_mol_m2_s, a.viscous_molar_flux_mol_m2_s);
        near(b.molar_flux_mol_m2_s, a.molar_flux_mol_m2_s);
    }
    near(b.entropy_production_w_m2_k, a.entropy_production_w_m2_k);
    near(b.net_molar_flux_mol_m2_s, a.net_molar_flux_mol_m2_s);
}

#[test]
fn heterogeneous_pure_gas_uses_series_knudsen_resistance_not_average_geometry() {
    let (left, right) = (pure("A", 110_000.0), pure("A", 90_000.0));
    let context = point(T, P);
    let lengths = [0.003, 0.007];
    let radii = [0.5 * RADIUS, 2.0 * RADIUS];
    let porosities = [0.25, 0.55];
    let tortuosities = [1.5, 3.0];
    let stack = layered(&[
        (0.0, lengths[0], pores(radii[0], porosities[0], tortuosities[0])),
        (0.0, lengths[1], pores(radii[1], porosities[1], tortuosities[1])),
    ]);
    let flux = left.layered_dusty_gas_flux_to(
        &right, &stack, AREA, &[], &context, controls(), || false,
    ).unwrap();

    let speed = (8.0 * UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * T
        / (core::f64::consts::PI * 0.028)).sqrt();
    let dk = [
        porosities[0] / tortuosities[0] * (2.0 * radii[0] / 3.0) * speed,
        porosities[1] / tortuosities[1] * (2.0 * radii[1] / 3.0) * speed,
    ];
    let resistance = lengths[0] / dk[0] + lengths[1] / dk[1];
    let effective = LENGTH / resistance;
    let gradient = (90_000.0 - 110_000.0) / LENGTH;
    let expected = -effective * gradient / (UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * T);
    near(flux.species[0].effective_knudsen_diffusivity_m2_s, effective);
    near(flux.species[0].knudsen_resistance_s_m, resistance);
    near(flux.species[0].diffusive_molar_flux_mol_m2_s, expected);
    assert_eq!(flux.species[0].viscous_molar_flux_mol_m2_s, 0.0);
    near(flux.species[0].molar_flux_mol_m2_s, expected);
}

#[test]
fn layered_constructor_refuses_klinkenberg_double_counting() {
    let slip_layer = PorousGasLayer::try_new(
        IntrinsicPermeability::isotropic(2e-14).unwrap(), LENGTH, 20_000.0,
        "synthetic slip layer",
    ).unwrap();
    let gas = binary(0.5, P);
    let barrier = DarcyGasBarrier::try_new_for_gas_model(
        vec![slip_layer], [1.0, 0.0, 0.0], support(), gas.model(),
    ).unwrap();
    let result = LayeredDustyGasMedium::try_new(
        barrier, vec![pores(RADIUS, POROSITY, TORTUOSITY)],
    );
    assert!(result.is_err());
}
