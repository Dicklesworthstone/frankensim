//! G0/G1/G3: independently sourced mixed Hill hardening and cyclic sheet memory.
//! Authored coefficients verify the law; they do not calibrate a rolled alloy.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::hill::oriented::OrientedHill48;
use fs_material::plastic::hill::source::{
    HILL_COEFFICIENT_PROPERTIES, HILL_HARDENING_MODULUS_PROPERTY, HILL_KINEMATIC_MODULUS_PROPERTY,
    HILL_REFERENCE_YIELD_STRESS_PROPERTY, ResolvedHill48, hill48_mixed_requirements,
    hill48_requirements, resolve_hill48, resolve_hill48_mixed,
};
use fs_material::plastic::hill::{Hill48, Hill48Coefficients};
use fs_material::plastic::section::energy::NativeSectionEnergyLaw;
use fs_material::plastic::source::KINEMATIC_HARDENING_MODULUS_PROPERTY;
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, POISSON_RATIO_PROPERTY, ScalarAdmissibility,
    ScalarPropertyRequirement, YOUNG_MODULUS_PROPERTY, resolve_material_state_point,
};
use fs_material::tensor::{contract, rotation};
use fs_material::{IsotropicElastic, SmallStrainLaw};
use fs_qty::{Density, Dims, Pressure};

const E: f64 = 1200.0;
const NU: f64 = 0.2;
const Y: f64 = 3.0;
const H: f64 = 40.0;
const C: f64 = 160.0;
const VALUES: [f64; 6] = [0.25, 0.9, 0.45, 1.4, 2.1, 1.8];
const IDENTITY: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

fn scalar(name: &str, dims: Dims, value: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained(),
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "authored mixed Hill coefficients".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}

fn properties(modulus: f64) -> Vec<PropertyClaim> {
    let mut properties = vec![
        scalar(YOUNG_MODULUS_PROPERTY, Pressure::DIMS, E),
        scalar(POISSON_RATIO_PROPERTY, Dims([0; 6]), NU),
        scalar(HILL_REFERENCE_YIELD_STRESS_PROPERTY, Pressure::DIMS, Y),
        scalar(HILL_HARDENING_MODULUS_PROPERTY, Pressure::DIMS, H),
        scalar(DENSITY_PROPERTY, Density::DIMS, 7800.0),
    ];
    properties.extend(
        HILL_COEFFICIENT_PROPERTIES
            .into_iter()
            .zip(VALUES)
            .map(|(name, value)| scalar(name, Dims([0; 6]), value)),
    );
    properties.push(scalar(
        HILL_KINEMATIC_MODULUS_PROPERTY,
        Pressure::DIMS,
        modulus,
    ));
    properties
}

fn card(properties: Vec<PropertyClaim>) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for property in properties {
        claims.insert_claim(property).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "authored rolled metal".into(),
            phase: "solid".into(),
            process: "mixed-hardening coupon fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn native(modulus: f64) -> Hill48 {
    let [f, g, h, l, m, n] = VALUES;
    Hill48::new(
        IsotropicElastic::new(E, NU, 0.05).unwrap(),
        Y,
        H,
        Hill48Coefficients { f, g, h, l, m, n },
    )
    .unwrap()
    .with_linear_kinematic_hardening(modulus)
    .unwrap()
}

fn gates() -> EnergyTolerances {
    EnergyTolerances {
        energy_j_m3: 1e-10,
        yield_pa: 1e-9,
    }
}

fn near(actual: f64, expected: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= 2e-11 * expected.abs().max(1.0),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn explicit_mixed_opt_in_preserves_all_original_receipts_and_native_cyclic_energy() {
    let card = card(properties(C));
    let point = QueryPoint::new().with("T", 350.0).unwrap();
    let frame = rotation([0.0, 0.0, 1.0], 0.43);
    let mut requirements = hill48_mixed_requirements().unwrap();
    assert_eq!(requirements.len(), 11);
    requirements.push(
        ScalarPropertyRequirement::try_new(
            DENSITY_PROPERTY,
            Density::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap(),
    );
    let bundle = resolve_material_state_point(
        &card,
        &point,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let base = ResolvedHill48::try_from_material(bundle.clone(), 0.05, frame).unwrap();
    assert_eq!(base.local_law().kinematic_modulus_pa(), 0.0);
    let resolved = base.with_linear_kinematic_hardening().unwrap();
    assert_eq!(resolved.material(), &bundle);
    assert_eq!(resolved.material().card_identity(), card.content_hash());
    assert_eq!(resolved.material().properties().len(), 12);
    for property in resolved.material().properties() {
        card.claims()
            .verify_receipt(&property.answer().receipt)
            .unwrap();
    }
    let expected = OrientedHill48::new(native(C), frame).unwrap();
    assert_eq!(resolved.law(), &expected);
    let mut old = resolved.initial_state();
    let mut previous = [0.0; 6];
    for strain in [
        [0.006, -0.002, -0.004, 0.002, 0.0, 0.0],
        [0.002, 0.003, -0.001, -0.003, 0.001, -0.001],
        [-0.005, 0.001, 0.002, -0.002, -0.001, 0.001],
    ] {
        let original = expected
            .advance_with_energy(&previous, &old, &strain, gates())
            .unwrap();
        let (state, stress, energy) = resolved
            .native_energy_increment(&previous, &old, &strain, gates())
            .unwrap();
        assert_eq!(resolved.response(&strain, &old).unwrap(), original.response);
        assert_eq!(state, original.response.state);
        assert_eq!(stress, original.response.stress_pa);
        assert_eq!(energy, original.energy);
        near(
            energy.final_state.kinematic_hardening_j_m3,
            0.5 * C * contract(&state.plastic_strain, &state.plastic_strain),
        );
        near(
            energy.physical_dissipation_j_m3,
            Y * (state.alpha - old.alpha),
        );
        assert!(energy.numerical_dissipation_j_m3 >= 0.0);
        old = state;
        previous = strain;
    }
    assert_eq!(resolved.material(), &bundle);
}

#[test]
fn sourced_directional_shear_has_the_independent_bauschinger_reverse_threshold() {
    let card = card(properties(C));
    let law = resolve_hill48_mixed(
        &card,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        0.05,
        IDENTITY,
    )
    .unwrap();
    let twice_shear = E / (1.0 + NU);
    for (index, a) in [(3, VALUES[5]), (4, VALUES[3]), (5, VALUES[4])] {
        let radius_factor = (2.0 * a).sqrt();
        let flow_factor = (a / 2.0).sqrt();
        let mut peak = [0.0; 6];
        peak[index] = 0.004;
        let alpha = (radius_factor * twice_shear * peak[index] - Y) / (H + (twice_shear + C) * a);
        let plastic = alpha * flow_factor;
        let loaded = law.response(&peak, &law.initial_state()).unwrap();
        near(loaded.state.alpha, alpha);
        near(loaded.state.plastic_strain[index], plastic);
        let radius = Y + H * alpha;
        let reverse_yield = C * plastic - radius / radius_factor;
        assert!(reverse_yield > -radius / radius_factor);
        let mut strain = [0.0; 6];
        strain[index] = plastic + (reverse_yield + 0.01 * Y / radius_factor) / twice_shear;
        let unloading = law.response(&strain, &loaded.state).unwrap();
        assert_eq!(unloading.state, loaded.state);
        near(
            unloading.stress_pa[index],
            reverse_yield + 0.01 * Y / radius_factor,
        );
        strain[index] = plastic + (C * plastic - 1.3 * radius / radius_factor) / twice_shear;
        let reverse = law.response(&strain, &loaded.state).unwrap();
        let increment = 0.3 * radius / (H + (twice_shear + C) * a);
        near(reverse.plastic_increment, increment);
        near(reverse.state.alpha, alpha + increment);
        near(
            reverse.state.plastic_strain[index],
            plastic - flow_factor * increment,
        );
        near(
            reverse.tangent_pa[index][index],
            twice_shear * (H + C * a) / (H + (twice_shear + C) * a),
        );
    }
}

#[test]
fn explicitly_sourced_zero_c_preserves_the_existing_law_and_work_exactly() {
    let card = card(properties(0.0));
    let point = QueryPoint::new();
    let frame = rotation([0.0, 0.0, 1.0], 0.31);
    let original = resolve_hill48(
        &card,
        &point,
        MaterialPropertySelection::SingleClaimOnly,
        0.05,
        frame,
    )
    .unwrap();
    let mixed = resolve_hill48_mixed(
        &card,
        &point,
        MaterialPropertySelection::SingleClaimOnly,
        0.05,
        frame,
    )
    .unwrap();
    assert_eq!(mixed.law(), original.law());
    assert_eq!(mixed.material().properties().len(), 11);
    assert_eq!(original.material().properties().len(), 10);
    let mut old = original.initial_state();
    let mut previous = [0.0; 6];
    for strain in [
        [0.004, -0.001, -0.003, 0.003, 0.0, 0.0],
        [-0.002, 0.003, -0.001, -0.002, 0.001, -0.001],
    ] {
        let a = original
            .native_energy_increment(&previous, &old, &strain, gates())
            .unwrap();
        let b = mixed
            .native_energy_increment(&previous, &old, &strain, gates())
            .unwrap();
        assert_eq!(a, b);
        assert_eq!(b.2.final_state.kinematic_hardening_j_m3, 0.0);
        old = a.0;
        previous = strain;
    }
}

#[test]
fn mixed_hill_requires_its_own_admissible_prager_source_without_j2_aliases() {
    let mut absent = properties(C);
    absent.pop();
    let absent = card(absent);
    let base = resolve_hill48(
        &absent,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        0.05,
        IDENTITY,
    )
    .unwrap();
    assert!(base.with_linear_kinematic_hardening().is_err());
    assert!(
        resolve_hill48_mixed(
            &absent,
            &QueryPoint::new(),
            MaterialPropertySelection::SingleClaimOnly,
            0.05,
            IDENTITY,
        )
        .is_err()
    );
    for (name, dims, value) in [
        (KINEMATIC_HARDENING_MODULUS_PROPERTY, Pressure::DIMS, C),
        (HILL_KINEMATIC_MODULUS_PROPERTY, Dims([0; 6]), C),
        (HILL_KINEMATIC_MODULUS_PROPERTY, Pressure::DIMS, -C),
    ] {
        let mut data = properties(C);
        *data.last_mut().unwrap() = scalar(name, dims, value);
        let card = card(data);
        assert!(
            resolve_hill48_mixed(
                &card,
                &QueryPoint::new(),
                MaterialPropertySelection::SingleClaimOnly,
                0.05,
                IDENTITY,
            )
            .is_err()
        );
        // A pre-resolved bundle with permissive scalar admission must not
        // bypass the opt-in constructor's own units, role or native C checks.
        let mut requirements = hill48_requirements().unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap(),
        );
        let bundle = resolve_material_state_point(
            &card,
            &QueryPoint::new(),
            &requirements,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        assert!(
            ResolvedHill48::try_from_material(bundle, 0.05, IDENTITY)
                .unwrap()
                .with_linear_kinematic_hardening()
                .is_err()
        );
    }
    let mut data = properties(C);
    let evolving = data.last_mut().unwrap();
    evolving.value = PropertyValue::Curve {
        abscissa: "strain".into(),
        abscissa_dims: Dims([0; 6]),
        knots: vec![(0.0, C), (0.02, 2.0 * C)],
        dims: Pressure::DIMS,
    };
    evolving.validity = ValidityDomain::unconstrained().with("strain", 0.0, 0.02);
    evolving.interpolation = InterpolationPolicy::LinearInside;
    assert!(
        resolve_hill48_mixed(
            &card(data),
            &QueryPoint::new().with("strain", 0.01).unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
            0.05,
            IDENTITY,
        )
        .is_err()
    );
}
