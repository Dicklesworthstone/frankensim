//! Wet constitutive snapshots consumed by the actual Darcy and DGM owners.
use super::*;
use crate::porous::permeability::INTRINSIC_PERMEABILITY_DIMS;
use crate::porous::wet::{GasDiffusivityLaw, WetPoreGasLaw, WetPoreGasState, WetPoreSupport};

fn wet_law(residual_liquid: f64, residual_gas: f64) -> WetPoreGasLaw {
    WetPoreGasLaw::try_new(
        residual_liquid,
        residual_gas,
        0.5,
        GasDiffusivityLaw::ConstantTransportTortuosity { factor: TORTUOSITY },
        WetPoreSupport {
            liquid_saturation: (0.0, 1.0),
            total_porosity: (0.1, 0.8),
            mean_gas_pore_radius_m: (RADIUS / 4.0, 4.0 * RADIUS),
        },
        "synthetic prescribed wet geometry; no coupled liquid dynamics",
    )
    .unwrap()
}
fn sample(law: &WetPoreGasLaw, saturation: f64) -> WetPoreGasState {
    law.sample(saturation, POROSITY, RADIUS).unwrap()
}
fn sample_pure(slab: &DustyGasMedium) -> DustyGasFlux {
    pure("A", 110_000.0)
        .dusty_gas_flux_to(
            &pure("A", 90_000.0),
            slab,
            AREA,
            &[],
            &point(T, P),
            controls(),
            || false,
        )
        .unwrap()
}

#[test]
fn wet_pure_flux_applies_independent_viscous_and_diffusive_responses_once() {
    let k = 2e-14;
    let original = medium(k, LENGTH);
    let law = wet_law(0.0, 0.0);
    let state = sample(&law, 0.5);
    let slab = original.with_wet_state(state.clone()).unwrap();
    let dry = sample_pure(&original);
    let wet = sample_pure(&slab);
    // Independent m=1/2 polynomial: sqrt(1-Se)*(1-Se^2).
    let kr = 0.5_f64.sqrt() * 0.75;
    let gradient = -20_000.0 / LENGTH;
    let viscosity = 1.7e-5 * (1.0 + (T - 200.0) / 1800.0);
    let expected_diffusion =
        -0.5 * knudsen(0.028) * gradient / (UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * T);
    let expected_viscous =
        -k * kr * P * gradient / (viscosity * UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * T);
    near(
        wet.species[0].diffusive_molar_flux_mol_m2_s,
        expected_diffusion,
    );
    near(wet.species[0].viscous_molar_flux_mol_m2_s, expected_viscous);
    near(
        wet.species[0].molar_flux_mol_m2_s,
        expected_diffusion + expected_viscous,
    );
    near(
        wet.species[0].knudsen_diffusivity_m2_s,
        0.5 * dry.species[0].knudsen_diffusivity_m2_s,
    );
    near(
        wet.darcy.as_ref().unwrap().mass_rate_left_to_right_kg_s,
        kr * dry.darcy.as_ref().unwrap().mass_rate_left_to_right_kg_s,
    );
    assert!(wet.local_entropy_production_w_m3_k > 0.0);
    assert_eq!(slab.wet_state(), Some(&state));
    assert_eq!(slab.reference_pores(), original.reference_pores());
    assert_eq!(
        slab.barrier().layers()[0].permeability(),
        original.barrier().layers()[0].permeability()
    );
}

#[test]
fn wet_binary_friction_scales_all_diffusivities_once_and_retains_original_pair_source() {
    let original = medium(2e-14, LENGTH);
    let slab = original
        .with_wet_state(sample(&wet_law(0.0, 0.0), 0.5))
        .unwrap();
    let (left, right) = (binary(0.75, P), binary(0.25, P));
    let (material, context) = (card(2e-5), point(T, P));
    let pairs = [pair("A", "B", &material, &context)];
    let dry = left
        .dusty_gas_flux_to(
            &right,
            &original,
            AREA,
            &pairs,
            &context,
            controls(),
            || false,
        )
        .unwrap();
    let wet = left
        .dusty_gas_flux_to(&right, &slab, AREA, &pairs, &context, controls(), || false)
        .unwrap();
    // At unchanged radius every wall and binary D is half its dry value;
    // the inverse friction operator therefore has exactly half the conductance.
    for (actual, reference) in wet.species.iter().zip(&dry.species) {
        near(
            actual.diffusive_molar_flux_mol_m2_s,
            0.5 * reference.diffusive_molar_flux_mol_m2_s,
        );
        assert_eq!(actual.viscous_molar_flux_mol_m2_s, 0.0);
    }
    assert_eq!(
        wet.pairs, dry.pairs,
        "wet effective diffusivity must not rewrite the material's original D_ij"
    );
    near(wet.pairs[0].diffusivity_m2_s, 2e-5);
    near(
        wet.local_entropy_production_w_m3_k,
        0.5 * dry.local_entropy_production_w_m3_k,
    );
    assert!(wet.maximum_gradient_residual_per_m < controls().gradient_tolerance_per_m);
}

#[test]
fn residual_gas_can_have_zero_darcy_mobility_and_nonzero_independent_pore_diffusion() {
    let original = medium(2e-14, LENGTH);
    let law = wet_law(0.1, 0.2);
    let state = sample(&law, 0.8);
    let slab = original.with_wet_state(state).unwrap();
    assert_eq!(slab.barrier().geometric_resistance_per_m(), None);
    assert!(slab.gas_filled_porosity() > 0.0);
    let dry = sample_pure(&original);
    let residual = sample_pure(&slab);
    assert_eq!(residual.species[0].viscous_molar_flux_mol_m2_s, 0.0);
    assert_eq!(
        residual
            .darcy
            .as_ref()
            .unwrap()
            .mass_rate_left_to_right_kg_s,
        0.0
    );
    near(
        residual.species[0].diffusive_molar_flux_mol_m2_s,
        0.2 * dry.species[0].diffusive_molar_flux_mol_m2_s,
    );
    assert!(residual.species[0].molar_flux_mol_m2_s > 0.0);
    near(
        residual.species[0].molar_flux_mol_m2_s,
        residual.species[0].diffusive_molar_flux_mol_m2_s,
    );
}

#[test]
fn fully_flooded_pores_close_exactly_without_querying_missing_or_invalid_binary_sources() {
    let original = medium(2e-14, LENGTH);
    let slab = original
        .with_wet_state(sample(&wet_law(0.0, 0.0), 1.0))
        .unwrap();
    let (left, right) = (binary(0.75, 110_000.0), binary(0.25, 90_000.0));
    let no_context = QueryPoint::new();
    let closed = left
        .dusty_gas_flux_to(&right, &slab, AREA, &[], &no_context, controls(), || false)
        .unwrap();
    assert_eq!(slab.gas_filled_porosity(), 0.0);
    assert_eq!(slab.effective_diffusivity_factor(), 0.0);
    assert_eq!(slab.transport_tortuosity_factor(), None);
    assert!(closed.species.is_empty() && closed.pairs.is_empty());
    assert!(closed.face_state.is_none() && closed.darcy.is_none());
    assert_eq!(closed.net_molar_flux_mol_m2_s, 0.0);
    assert_eq!(closed.local_entropy_production_w_m3_k, 0.0);
    assert_eq!(closed.maximum_gradient_residual_per_m, 0.0);
    let wrong = source(
        PropertyValue::Scalar {
            value: 2e-5,
            dims: INTRINSIC_PERMEABILITY_DIMS,
        },
        INTRINSIC_PERMEABILITY_DIMS,
        InterpolationPolicy::ConstantWithinValidity,
        [95_000.0, 105_000.0],
    );
    let invalid_pairs = [pair("A", "B", &wrong, &no_context)];
    let closed_with_invalid = left
        .dusty_gas_flux_to(
            &right,
            &slab,
            AREA,
            &invalid_pairs,
            &no_context,
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(closed_with_invalid, closed);
    assert!(
        left.dusty_gas_flux_to(
            &right,
            &original,
            AREA,
            &invalid_pairs,
            &no_context,
            controls(),
            || false
        )
        .is_err(),
        "the deliberately invalid source must fail when the pore path is active"
    );
    assert!(matches!(
        left.dusty_gas_flux_to(&right, &slab, AREA, &[], &no_context, controls(), || true),
        Err(GasEnergyError::Cancelled)
    ));
}

#[test]
fn wet_updates_replace_prior_factors_and_dry_recovery_restores_original_transport() {
    let original = medium(2e-14, LENGTH);
    let saved = original.clone();
    let law = wet_law(0.0, 0.0);
    let first = original.with_wet_state(sample(&law, 0.25)).unwrap();
    let through_history = first.with_wet_state(sample(&law, 0.65)).unwrap();
    let direct = original.with_wet_state(sample(&law, 0.65)).unwrap();
    assert_eq!(through_history, direct);
    assert_eq!(sample_pure(&through_history), sample_pure(&direct));
    let flooded = through_history.with_wet_state(sample(&law, 1.0)).unwrap();
    let recovered = flooded.with_wet_state(sample(&law, 0.0)).unwrap();
    let dry_flux = sample_pure(&original);
    let recovered_flux = sample_pure(&recovered);
    assert_eq!(recovered_flux.species, dry_flux.species);
    assert_eq!(recovered_flux.darcy, dry_flux.darcy);
    assert_eq!(
        recovered_flux.local_entropy_production_w_m3_k,
        dry_flux.local_entropy_production_w_m3_k
    );
    assert_eq!(recovered.reference_pores(), original.reference_pores());
    assert_eq!(
        original, saved,
        "sampling updates must not mutate the reference medium"
    );
}

#[test]
fn prescribed_gas_pore_radius_changes_knudsen_transport_without_changing_darcy_mobility() {
    let original = medium(2e-14, LENGTH);
    let law = wet_law(0.0, 0.0);
    let narrow = original
        .with_wet_state(law.sample(0.4, POROSITY, RADIUS).unwrap())
        .unwrap();
    let wide = narrow
        .with_wet_state(law.sample(0.4, POROSITY, 2.0 * RADIUS).unwrap())
        .unwrap();
    let a = sample_pure(&narrow);
    let b = sample_pure(&wide);
    near(
        b.species[0].knudsen_diffusivity_m2_s,
        2.0 * a.species[0].knudsen_diffusivity_m2_s,
    );
    near(
        b.species[0].diffusive_molar_flux_mol_m2_s,
        2.0 * a.species[0].diffusive_molar_flux_mol_m2_s,
    );
    assert_eq!(
        b.species[0].viscous_molar_flux_mol_m2_s,
        a.species[0].viscous_molar_flux_mol_m2_s
    );
    assert_eq!(
        wide.barrier().geometric_resistance_per_m(),
        narrow.barrier().geometric_resistance_per_m()
    );
    assert_eq!(wide.reference_pores(), original.reference_pores());
}

fn source_bound_medium() -> (DustyGasMedium, ResolvedMaterialStatePoint) {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("intrinsic_permeability", INTRINSIC_PERMEABILITY_DIMS),
            value: PropertyValue::Scalar {
                value: 2e-14,
                dims: INTRINSIC_PERMEABILITY_DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "temperature",
                    typed(QuantityKind::AbsoluteTemperature),
                    300.0,
                    500.0,
                )
                .with_quantity(
                    "pressure",
                    typed(QuantityKind::Pressure),
                    90_000.0,
                    110_000.0,
                ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic original dry skeleton".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic porous solid".into(),
            phase: "solid".into(),
            process: "source-preservation fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let requirement = ScalarPropertyRequirement::try_new(
        "intrinsic_permeability",
        INTRINSIC_PERMEABILITY_DIMS,
        ScalarAdmissibility::NonNegative,
    )
    .unwrap();
    let source = resolve_material_state_point(
        &card,
        &point(T, P),
        &[requirement],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let layer = PorousGasLayer::from_material_isotropic(
        source.clone(),
        "intrinsic_permeability",
        LENGTH,
        None,
    )
    .unwrap();
    let barrier = DarcyGasBarrier::try_new(vec![layer], [1.0, 0.0, 0.0], support()).unwrap();
    let pores = DustyGasPores::try_new(RADIUS, POROSITY, TORTUOSITY, "reference geometry").unwrap();
    (DustyGasMedium::try_new(barrier, pores).unwrap(), source)
}

#[test]
fn wet_snapshots_preserve_material_identity_support_and_refuse_dry_slip_calibration() {
    let (original, source) = source_bound_medium();
    let law = wet_law(0.0, 0.0);
    let state = sample(&law, 0.5);
    let wet = original.with_wet_state(state.clone()).unwrap();
    assert_eq!(wet.barrier().layers()[0].source(), Some(&source));
    assert_eq!(
        wet.barrier().layers()[0].permeability(),
        original.barrier().layers()[0].permeability()
    );
    assert_eq!(
        wet.barrier().layers()[0].provenance(),
        original.barrier().layers()[0].provenance()
    );
    assert_eq!(wet.barrier().support(), original.barrier().support());
    assert_eq!(wet.wet_state().unwrap().law(), &law);
    sample_pure(&wet);
    let outside = pure("A", 89_999.0);
    assert!(
        pure("A", 110_000.0)
            .dusty_gas_flux_to(&outside, &wet, AREA, &[], &point(T, P), controls(), || {
                false
            })
            .is_err()
    );
    let slip = PorousGasLayer::try_new(
        IntrinsicPermeability::isotropic(2e-14).unwrap(),
        LENGTH,
        5000.0,
        "explicit dry gas slip calibration",
    )
    .unwrap();
    assert!(slip.with_wet_state(state.clone()).is_err());
    assert!(slip.wet_state().is_none());
    let gas = pure("A", P);
    let slip_barrier =
        DarcyGasBarrier::try_new_for_gas_model(vec![slip], [1.0, 0.0, 0.0], support(), gas.model())
            .unwrap();
    assert!(slip_barrier.with_wet_states(&[state]).is_err());
}

#[test]
fn layered_anisotropic_darcy_uses_each_wet_mobility_once_and_closes_with_one_flooded_layer() {
    let q = core::f64::consts::FRAC_1_SQRT_2;
    let tensor = IntrinsicPermeability::orthotropic(
        [2e-14, 8e-14, 0.0],
        [[q, q, 0.0], [-q, q, 0.0], [0.0, 0.0, 1.0]],
    )
    .unwrap();
    let layers = vec![
        PorousGasLayer::try_new(tensor, LENGTH, 0.0, "rotated dry layer").unwrap(),
        PorousGasLayer::try_new(
            IntrinsicPermeability::isotropic(1e-14).unwrap(),
            2.0 * LENGTH,
            0.0,
            "second dry layer",
        )
        .unwrap(),
    ];
    let original = DarcyGasBarrier::try_new(layers.clone(), [1.0, 0.0, 0.0], support()).unwrap();
    let law = wet_law(0.0, 0.0);
    let states = [sample(&law, 0.25), sample(&law, 0.5)];
    let wet = original.with_wet_states(&states).unwrap();
    let kr1 = 0.75_f64.sqrt() * (1.0 - 0.25 * 0.25);
    let kr2 = 0.5_f64.sqrt() * 0.75;
    let r1 = LENGTH / (5e-14 * kr1);
    let r2 = 2.0 * LENGTH / (1e-14 * kr2);
    near(wet.geometric_resistance_per_m().unwrap(), r1 + r2);
    let parcel = pure("A", P);
    let gas = parcel.state().state();
    let profile = wet
        .mass_flow_profile(gas, 110_000.0, 90_000.0, AREA)
        .unwrap();
    let expected =
        AREA * 20_000.0 * P / (gas.dynamic_viscosity * gas.specific_gas_constant * T * (r1 + r2));
    near(profile.flow.mass_rate_left_to_right_kg_s, expected);
    let interface = (110_000.0_f64.powi(2)
        - r1 / (r1 + r2) * (110_000.0_f64.powi(2) - 90_000.0_f64.powi(2)))
    .sqrt();
    near(
        profile.interface_pressures_pa.as_ref().unwrap()[1],
        interface,
    );
    assert_eq!(profile.iterations, 0);
    assert_eq!(wet.with_wet_states(&states).unwrap(), wet);
    for (updated, dry) in wet.layers().iter().zip(&layers) {
        assert_eq!(updated.permeability(), dry.permeability());
    }
    assert!(original.with_wet_states(&states[..1]).is_err());
    let sealed = wet
        .with_wet_states(&[states[0].clone(), sample(&law, 1.0)])
        .unwrap();
    let stopped = sealed
        .mass_flow_profile(gas, 110_000.0, 90_000.0, AREA)
        .unwrap();
    assert_eq!(stopped.flow.mass_rate_left_to_right_kg_s, 0.0);
    assert!(stopped.interface_pressures_pa.is_none());
}
