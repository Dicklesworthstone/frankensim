//! G0/G1/G3: thermal weakening keeps fracture history and explicit parameter work.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::brittle::source::*;
use fs_material::brittle::{
    BrittleCrackBand, BrittleCrackBandParameters, BrittleCrackBandSettings,
};
use fs_material::cohesive::source::FRACTURE_ENERGY_PER_AREA_DIMS;
use fs_material::elastic::IsotropicElastic;
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, POISSON_RATIO_PROPERTY, ScalarAdmissibility,
    ScalarPropertyRequirement, YOUNG_MODULUS_PROPERTY, resolve_material_state_point,
};
use fs_material::tensor::{rotate, rotation};
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature};

fn settings() -> BrittleCrackBandSettings {
    BrittleCrackBandSettings {
        band_width_m: 0.01,
        normal: [1.0, 0.0, 0.0],
        maximum_strain: 0.05,
        energy_tolerance_j_m3: 1e-6,
    }
}

fn law(e: f64, strength: f64, energy: f64) -> BrittleCrackBand {
    BrittleCrackBand::new(
        IsotropicElastic::new(1e9 * e, 0.2, 0.05).unwrap(),
        BrittleCrackBandParameters {
            tensile_strength_pa: 1e6 * strength,
            shear_strength_pa: 0.8e6 * strength,
            mode_i_fracture_energy_j_m2: 100.0 * energy,
            mode_ii_fracture_energy_j_m2: 200.0 * energy,
            bk_exponent: 1.5,
        },
        settings(),
    )
    .unwrap()
}

fn near(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 1e-9 * a.abs().max(b.abs()).max(1.0),
        "{a:e} != {b:e}"
    );
}

#[test]
fn modulus_change_accounts_for_all_elastic_storage_at_held_strain() {
    let cold = law(1.0, 1.0, 1.0);
    let warm = law(0.5, 1.0, 1.0);
    let strain = [0.0003, 0.0002, -0.0001, 0.0001, 0.0002, -0.00015];
    let old = cold.trial(&strain, &cold.initial_state()).unwrap();
    assert_eq!(old.state.damage(), 0.0);
    let transferred = old.state.transfer_temperature(&warm, 1e-6).unwrap();
    near(
        transferred.parameter_work_j_m3,
        -0.5 * old.state.stored_energy_j_m3(),
    );
    near(
        transferred.stored_change_j_m3,
        transferred.parameter_work_j_m3,
    );
    assert_eq!(transferred.damage_dissipated_j_m3, 0.0);
    assert_eq!(transferred.state.strain(), &strain);
    let held = warm.trial(&strain, &transferred.state).unwrap();
    assert_eq!(held.energy.work, 0.0);
    for (new, old) in held.stress_pa.into_iter().zip(old.stress_pa) {
        near(new, 0.5 * old);
    }
    near(transferred.defect_j_m3, 0.0);
}

#[test]
fn thermal_failure_keeps_past_fracture_loss_and_does_not_heal_on_cooling() {
    let cold = law(1.0, 1.0, 1.0);
    let hot = law(0.5, 0.5, 0.03);
    let strain = [0.004, 0.0002, -0.0001, 0.0003, 0.0004, -0.0002];
    let before = cold.trial(&strain, &cold.initial_state()).unwrap().state;
    assert!(before.damage() > 0.0 && !before.failed());
    let saved = before.clone();
    let heated = before.transfer_temperature(&hot, 1e-6).unwrap();
    assert!(heated.state.failed());
    assert!(heated.damage_dissipated_j_m3 > 0.0);
    near(
        heated.parameter_work_j_m3,
        -0.5 * before.stored_energy_j_m3(),
    );
    near(
        heated.state.dissipated_j_m3(),
        before.dissipated_j_m3() + heated.damage_dissipated_j_m3,
    );
    near(
        heated.parameter_work_j_m3,
        heated.stored_change_j_m3 + heated.damage_dissipated_j_m3,
    );
    let cooled = heated.state.transfer_temperature(&cold, 1e-6).unwrap();
    assert!(cooled.state.failed());
    assert_eq!(cooled.damage_dissipated_j_m3, 0.0);
    near(
        cooled.state.dissipated_j_m3(),
        heated.state.dissipated_j_m3(),
    );
    assert_eq!(before, saved);
    // Reclosure restores the original compressive stiffness, not cohesion.
    let compressed = cold
        .trial(&[-0.001, 0.0, 0.0, 0.0, 0.0, 0.0], &cooled.state)
        .unwrap();
    assert!(compressed.state.failed());
    near(
        compressed.stress_pa[0],
        -1e9 * (1.0 - 0.2) / (1.2 * 0.6) * 0.001,
    );
    assert_eq!(compressed.energy.dissipated, 0.0);
}

#[test]
fn toughness_increase_preserves_damage_and_accumulated_fracture_energy() {
    let old = law(1.0, 1.0, 1.0);
    let tough = law(1.0, 1.0, 3.0);
    let loaded = old
        .trial(&[0.004, 0.0, 0.0, 0.001, 0.0, 0.0], &old.initial_state())
        .unwrap();
    let unloaded = old
        .trial(&[0.001, 0.0, 0.0, 0.00025, 0.0, 0.0], &loaded.state)
        .unwrap();
    let step = unloaded.state.transfer_temperature(&tough, 1e-6).unwrap();
    assert!(step.state.damage() >= unloaded.state.damage());
    near(step.state.damage(), unloaded.state.damage());
    near(
        step.state.dissipated_j_m3(),
        unloaded.state.dissipated_j_m3(),
    );
    near(step.parameter_work_j_m3, 0.0);
    near(step.damage_dissipated_j_m3, 0.0);
    let held = tough.trial(step.state.strain(), &step.state).unwrap();
    for (a, b) in held.stress_pa.into_iter().zip(unloaded.stress_pa) {
        near(a, b);
    }
    let later = tough
        .trial(&[0.012, 0.0, 0.0, 0.003, 0.0, 0.0], &step.state)
        .unwrap();
    assert!(later.state.damage() >= step.state.damage());
    assert!(later.state.dissipated_j_m3() >= step.state.dissipated_j_m3());
}

#[test]
fn compression_parameter_work_and_rotated_mixed_history_are_objective() {
    let old = law(1.0, 1.0, 1.0);
    let next = law(0.7, 0.8, 0.4);
    let axis_scale = fs_math::det::sqrt(6.0);
    let r = rotation(
        [1.0 / axis_scale, 2.0 / axis_scale, -1.0 / axis_scale],
        0.73,
    );
    let rotate_law = |base: &BrittleCrackBand| {
        BrittleCrackBand::new(
            base.elastic().clone(),
            base.parameters(),
            BrittleCrackBandSettings {
                normal: [r[0][0], r[1][0], r[2][0]],
                ..base.settings()
            },
        )
        .unwrap()
    };
    let rotated_old = rotate_law(&old);
    let rotated_next = rotate_law(&next);
    for strain in [
        [-0.002, 0.0001, -0.0002, 0.0, 0.0, 0.0],
        [0.004, 0.0003, -0.0001, 0.0007, 0.0002, -0.0004],
    ] {
        let a = old.trial(&strain, &old.initial_state()).unwrap().state;
        let rotated_strain = rotate(&strain, &r);
        let b = rotated_old
            .trial(&rotated_strain, &rotated_old.initial_state())
            .unwrap()
            .state;
        let x = a.transfer_temperature(&next, 1e-6).unwrap();
        let y = b.transfer_temperature(&rotated_next, 1e-6).unwrap();
        near(x.parameter_work_j_m3, y.parameter_work_j_m3);
        near(x.damage_dissipated_j_m3, y.damage_dissipated_j_m3);
        near(x.state.damage(), y.state.damage());
        near(x.state.dissipated_j_m3(), y.state.dissipated_j_m3());
        let xs = next.trial(&strain, &x.state).unwrap().stress_pa;
        let ys = rotated_next
            .trial(&rotated_strain, &y.state)
            .unwrap()
            .stress_pa;
        for (expected, actual) in rotate(&xs, &r).into_iter().zip(ys) {
            near(expected, actual);
        }
    }
}

#[test]
fn changed_plane_width_poisson_and_inadmissible_strain_refuse_transactionally() {
    let old = law(1.0, 1.0, 1.0);
    let state = old
        .trial(&[0.004, 0.0, 0.0, 0.001, 0.0, 0.0], &old.initial_state())
        .unwrap()
        .state;
    let saved = state.clone();
    for (poisson, settings) in [
        (
            0.2,
            BrittleCrackBandSettings {
                normal: [0.0, 1.0, 0.0],
                ..settings()
            },
        ),
        (
            0.2,
            BrittleCrackBandSettings {
                band_width_m: 0.02,
                ..settings()
            },
        ),
        (0.25, settings()),
        (
            0.2,
            BrittleCrackBandSettings {
                maximum_strain: 0.003,
                ..settings()
            },
        ),
    ] {
        let target = BrittleCrackBand::new(
            IsotropicElastic::new(1e9, poisson, 0.05).unwrap(),
            old.parameters(),
            settings,
        )
        .unwrap();
        assert!(state.transfer_temperature(&target, 1e-6).is_err());
    }
    for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(state.transfer_temperature(&old, invalid).is_err());
    }
    assert_eq!(state, saved);
    old.trial(&[0.005, 0.0, 0.0, 0.001, 0.0, 0.0], &state)
        .unwrap();
}

fn source() -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, cold, hot) in [
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, 1e9, 0.5e9),
        (POISSON_RATIO_PROPERTY, Dims::NONE, 0.2, 0.2),
        (
            BRITTLE_TENSILE_STRENGTH_PROPERTY,
            Pressure::DIMS,
            1e6,
            0.5e6,
        ),
        (
            BRITTLE_SHEAR_STRENGTH_PROPERTY,
            Pressure::DIMS,
            0.8e6,
            0.4e6,
        ),
        (
            BRITTLE_MODE_I_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
            100.0,
            3.0,
        ),
        (
            BRITTLE_MODE_II_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
            200.0,
            6.0,
        ),
        (BRITTLE_BK_EXPONENT_PROPERTY, Dims::NONE, 1.5, 1.5),
        (DENSITY_PROPERTY, Density::DIMS, 2300.0, 2250.0),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(300.0, cold), (400.0, hot)],
                    dims,
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    300.0,
                    400.0,
                ),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured thermal brittle test; no measured material claim".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic brittle thermal material".into(),
            phase: "solid".into(),
            process: "numerical verification".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn resolved(card: &MaterialCard) -> ResolvedBrittleCrackBand {
    let mut requirements = brittle_crack_band_requirements().unwrap();
    requirements.push(
        ScalarPropertyRequirement::try_new(
            DENSITY_PROPERTY,
            Density::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap(),
    );
    let query = QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 300.0)
        .unwrap();
    ResolvedBrittleCrackBand::try_from_material(
        resolve_material_state_point(
            card,
            &query,
            &requirements,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap(),
        settings(),
    )
    .unwrap()
}

#[test]
fn original_source_claims_and_extra_properties_survive_heating_and_cooling() {
    let card = source();
    let cold = resolved(&card);
    let hot = cold.at_temperature(&card, "T", 400.0).unwrap();
    assert_eq!(
        hot.material().card_identity(),
        cold.material().card_identity()
    );
    assert_eq!(hot.material().properties().len(), 8);
    near(
        hot.material()
            .property(DENSITY_PROPERTY)
            .unwrap()
            .value_si(),
        2250.0,
    );
    for (a, b) in cold
        .material()
        .properties()
        .iter()
        .zip(hot.material().properties())
    {
        assert_eq!(a.requirement(), b.requirement());
        assert_eq!(a.answer().receipt.selected, b.answer().receipt.selected);
    }
    let before = cold
        .law()
        .trial(
            &[0.004, 0.0002, -0.0001, 0.0003, 0.0004, -0.0002],
            &cold.law().initial_state(),
        )
        .unwrap()
        .state;
    // A coefficient refresh alone never authorizes a raw state swap.
    assert!(hot.law().trial(before.strain(), &before).is_err());
    let heated = cold.transfer_temperature(&before, &hot, "T", 1e-6).unwrap();
    assert!(heated.state.failed());
    let refreshed_cold = hot.at_temperature(&card, "T", 300.0).unwrap();
    let cooled = hot
        .transfer_temperature(&heated.state, &refreshed_cold, "T", 1e-6)
        .unwrap();
    assert!(cooled.state.failed());
    near(
        cooled.state.dissipated_j_m3(),
        heated.state.dissipated_j_m3(),
    );
    for invalid in [250.0, 401.0, 0.0, f64::NAN] {
        assert!(cold.at_temperature(&card, "T", invalid).is_err());
    }
    assert!(
        cold.at_temperature(&card, "missing-temperature", 350.0)
            .is_err()
    );
    assert!(hot.transfer_temperature(&before, &cold, "T", 1e-6).is_err());
    let clone = cold.clone();
    assert!(
        clone
            .transfer_temperature(&before, &hot, "T", 1e-6)
            .is_err()
    );
    let oblique = ResolvedBrittleCrackBand::try_from_material(
        cold.material().clone(),
        BrittleCrackBandSettings {
            normal: [0.5631578947368421, -0.18461538461538463, 0.8648648648648649],
            ..settings()
        },
    )
    .unwrap();
    let oblique_hot = oblique.at_temperature(&card, "T", 400.0).unwrap();
    assert_eq!(oblique.law().local_axes(), oblique_hot.law().local_axes());
    let state = oblique
        .law()
        .trial(
            &[0.004, 0.0, 0.0, 0.001, 0.0, 0.0],
            &oblique.law().initial_state(),
        )
        .unwrap()
        .state;
    let heated = oblique
        .transfer_temperature(&state, &oblique_hot, "T", 1e-6)
        .unwrap();
    assert_eq!(
        heated.state.cohesive_point().jump_m(),
        state.cohesive_point().jump_m()
    );
}
