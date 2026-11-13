//! Joint original exposed-depth/native-history hardness, with independent
//! numerical references. Synthetic fixtures are not material qualification.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, HardnessLoadStep, HardnessTestContext, InterpolationPolicy, MaterialCard,
    MaterialStateId, ObservationDataset, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    QueryPoint, UncertaintyModel,
};
use fs_material::SmallStrainLaw;
use fs_material::hardness::INDENTATION_HARDNESS_QUANTITY;
use fs_material::hardness::depth::history::DepthHardnessHistoryError;
use fs_material::hardness::depth::{
    DepthHardnessError, EXPOSED_MATERIAL_DEPTH_AXIS as Z, ResolvedDepthHardness,
};
use fs_material::hardness::plastic_history::ACCUMULATED_EQUIVALENT_PLASTIC_STRAIN_AXIS as ALPHA;
use fs_material::plastic::chaboche::source::{BackstressProperties, resolve_chaboche_branches};
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::source::{
    HardeningPropertyBasis, ISOTROPIC_HARDENING_MODULUS_PROPERTY,
    KINEMATIC_HARDENING_MODULUS_PROPERTY, ResolvedJ2Plasticity, resolve_j2_plasticity,
};
use fs_material::plastic::voce::source::{
    VOCE_SATURATION_INCREMENT_PROPERTY, VOCE_SATURATION_RATE_PROPERTY,
};
use fs_material::state_point::{
    MaterialPropertySelection, POISSON_RATIO_PROPERTY, YIELD_STRESS_PROPERTY,
    YOUNG_MODULUS_PROPERTY,
};
use fs_qty::{Dims, Pressure, QtyAny, QuantitySpec};

const SELECT: MaterialPropertySelection = MaterialPropertySelection::SingleClaimOnly;
const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);
const MAX_DEPTH: f64 = 0.002;
const SLOPE_ALPHA: f64 = 10000.0;
const GATES: EnergyTolerances = EnergyTolerances {
    energy_j_m3: 1e-8,
    yield_pa: 1e-8,
};

fn query() -> QueryPoint {
    QueryPoint::new().with("T", 300.0).unwrap()
}
fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-11 * expected.abs().max(1e-8),
        "{actual:e} != {expected:e}"
    );
}
fn surface(spans: &[&[(f64, f64)]]) -> PropertyValue {
    let mut samples = Vec::new();
    let mut triangles = Vec::new();
    for knots in spans {
        let base = samples.len() as u32;
        for &(z, h) in *knots {
            for alpha in [0.0, 0.04] {
                samples.push(([z, alpha], h + SLOPE_ALPHA * alpha));
            }
        }
        for interval in 0..knots.len() - 1 {
            let i = base + 2 * interval as u32;
            triangles.extend([[i, i + 2, i + 3], [i, i + 3, i + 1]]);
        }
    }
    PropertyValue::Surface2D {
        axes: [Z.into(), ALPHA.into()],
        axis_dims: [LENGTH, Dims::NONE],
        samples,
        triangles,
        dims: Pressure::DIMS,
    }
}
fn fixture(value: PropertyValue) -> (MaterialCard, PropertyKey) {
    let provenance = Provenance {
        source: "synthetic joint depth/history hardness".into(),
        license: "CC0-1.0".into(),
        artifact: None,
    };
    let environment = ValidityDomain::unconstrained().with("T", 280.0, 320.0);
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, 1000.0),
        (POISSON_RATIO_PROPERTY, Dims::NONE, 0.3),
        (YIELD_STRESS_PROPERTY, Pressure::DIMS, 1.0),
        (ISOTROPIC_HARDENING_MODULUS_PROPERTY, Pressure::DIMS, 40.0),
        (KINEMATIC_HARDENING_MODULUS_PROPERTY, Pressure::DIMS, 20.0),
        (VOCE_SATURATION_INCREMENT_PROPERTY, Pressure::DIMS, 0.6),
        (VOCE_SATURATION_RATE_PROPERTY, Dims::NONE, 80.0),
        ("fast_C", Pressure::DIMS, 300.0),
        ("fast_gamma", Dims::NONE, 250.0),
        ("slow_C", Pressure::DIMS, 80.0),
        ("slow_gamma", Dims::NONE, 30.0),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: environment.clone(),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: provenance.clone(),
            })
            .unwrap();
    }
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: "synthetic uniformly preconditioned graded specimen".into(),
            method: "synthetic force/projected-area indentation".into(),
            artifact: fs_blake3::hash_domain("fixture", b"native history and depth"),
            caveats: "numerical regression only".into(),
            provenance: provenance.clone(),
        })
        .unwrap();
    let key = PropertyKey::with_quantity("H", INDENTATION_HARDNESS_QUANTITY)
        .with_hardness_test(
            HardnessTestContext::new(
                "synthetic sphere radius 1 mm",
                vec![
                    HardnessLoadStep::new(
                        QtyAny::new(10.0, Dims([1, 1, -2, 0, 0, 0])),
                        QtyAny::new(1.0, Dims([0, 0, 1, 0, 0, 0])),
                    )
                    .unwrap(),
                ],
                "synthetic depth/history protocol",
                observation,
            )
            .unwrap(),
        )
        .unwrap();
    claims
        .insert_claim(PropertyClaim {
            key: key.clone(),
            value,
            validity: environment
                .with_quantity(Z, QuantitySpec::dimensional(LENGTH), 0.0, MAX_DEPTH)
                .with_quantity(ALPHA, QuantitySpec::dimensional(Dims::NONE), 0.0, 0.04),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![observation],
            provenance,
        })
        .unwrap();
    (
        MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic graded solid".into(),
                phase: "solid".into(),
                process: "synthetic uniform native preconditioning".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap(),
        key,
    )
}
fn ordinary_fixture() -> (MaterialCard, PropertyKey) {
    fixture(surface(&[&[
        (0.0, 300.0),
        (0.001, 180.0),
        (MAX_DEPTH, 100.0),
    ]]))
}
fn plastic(card: &MaterialCard) -> ResolvedJ2Plasticity {
    resolve_j2_plasticity(
        card,
        &query(),
        SELECT,
        HardeningPropertyBasis::PlasticStrainModulus,
        0.05,
    )
    .unwrap()
}
fn base_hardness(z: f64) -> f64 {
    if z <= 0.001 {
        300.0 - 120000.0 * z
    } else {
        180.0 - 80000.0 * (z - 0.001)
    }
}

#[test]
fn g1_native_loading_unloading_and_reversal_set_a_joint_positive_depth_profile() {
    let (card, key) = ordinary_fixture();
    let plastic = plastic(&card);
    let mut old_strain = [0.0; 6];
    let mut state = plastic.initial_state();
    let mut previous: Option<ResolvedDepthHardness> = None;
    for (step, axial) in [0.006, 0.00599, -0.006].into_iter().enumerate() {
        let strain = [axial, 0.0, 0.0, 0.0, 0.0, 0.0];
        let native = plastic
            .law()
            .advance_with_energy(&old_strain, &state, &strain, GATES)
            .unwrap();
        let before = native.clone();
        let profile = ResolvedDepthHardness::from_j2_accepted(
            &plastic,
            &strain,
            &native.state,
            GATES.yield_pa,
            &card,
            key.clone(),
            SELECT,
            MAX_DEPTH,
        )
        .unwrap();
        let history = profile.accepted_history().unwrap();
        assert_eq!(history.mechanical_source(), plastic.material());
        assert_eq!(history.alpha(), native.state.alpha);
        if step == 0 {
            let mu = 1000.0 / 2.6;
            near(
                history.alpha(),
                (2.0 * mu * axial - 1.0) / (3.0 * mu + 40.0),
            );
        }
        for z in [0.0, 0.0003, 0.001, 0.0017, MAX_DEPTH] {
            near(
                profile.hardness_pa_at(z).unwrap(),
                base_hardness(z) + SLOPE_ALPHA * history.alpha(),
            );
            let endpoint = profile.resolve_at_depth(&card, z).unwrap();
            near(endpoint.hardness_pa(), profile.hardness_pa_at(z).unwrap());
            assert_eq!(endpoint.key(), &key);
            assert_eq!(
                endpoint.resolved().state_coordinate(ALPHA),
                Some(history.alpha())
            );
            assert_eq!(endpoint.resolved().state_coordinate(Z), Some(z));
            assert_eq!(
                endpoint.property().answer().receipt.selected,
                profile
                    .reference_source()
                    .property()
                    .answer()
                    .receipt
                    .selected
            );
        }
        let integral = 0.38 + SLOPE_ALPHA * history.alpha() * MAX_DEPTH;
        near(
            profile.hardness_integral_between(0.0, MAX_DEPTH).unwrap(),
            integral,
        );
        let interior_integral = 0.3464 + SLOPE_ALPHA * history.alpha() * 0.0017;
        near(
            profile
                .advance_by_integral(0.0, interior_integral)
                .unwrap()
                .depth_m,
            0.0017,
        );
        if let Some(old) = previous {
            assert_eq!(
                history.native_source_identity(),
                old.accepted_history().unwrap().native_source_identity()
            );
            if step == 1 {
                assert_eq!(profile, old);
            } else {
                assert!(
                    profile.hardness_pa_at(0.001).unwrap() > old.hardness_pa_at(0.001).unwrap()
                );
            }
        }
        assert_eq!(profile.at_temperature(&card, 300.0).unwrap(), profile);
        assert!(profile.at_temperature(&card, 301.0).is_err());
        assert_eq!(
            native, before,
            "property preparation must not alter native work or history"
        );
        previous = Some(profile);
        old_strain = strain;
        state = native.state;
    }
}

#[test]
fn g0_scalar_alpha_invalid_history_and_interior_source_failures_are_not_admitted() {
    let (card, key) = ordinary_fixture();
    let native = plastic(&card);
    let state = native.initial_state();
    let arbitrary = query()
        .with_quantity(ALPHA, QuantitySpec::dimensional(Dims::NONE), 0.01)
        .unwrap();
    assert!(ResolvedDepthHardness::new(&card, &arbitrary, key.clone(), SELECT, MAX_DEPTH).is_err());
    let strain = [0.006, 0.0, 0.0, 0.0, 0.0, 0.0];
    assert!(matches!(
        ResolvedDepthHardness::from_j2_accepted(
            &native,
            &strain,
            &state,
            GATES.yield_pa,
            &card,
            key.clone(),
            SELECT,
            MAX_DEPTH,
        ),
        Err(DepthHardnessHistoryError::NativeJ2(_))
    ));
    // The law's small-strain bound is 0.05 (Frobenius norm of the strain), so the history
    // leaves the surface's alpha range (0 to 0.04) by a reversal at 0.045, not one 0.09 step.
    let loaded = native
        .law()
        .advance_with_energy(&[0.0; 6], &state, &[0.045, 0.0, 0.0, 0.0, 0.0, 0.0], GATES)
        .unwrap();
    let outside = native
        .law()
        .advance_with_energy(
            &loaded.strain,
            &loaded.state,
            &[-0.045, 0.0, 0.0, 0.0, 0.0, 0.0],
            GATES,
        )
        .unwrap();
    assert!(outside.state.alpha > 0.04);
    assert!(
        ResolvedDepthHardness::from_j2_accepted(
            &native,
            &outside.strain,
            &outside.state,
            GATES.yield_pa,
            &card,
            key.clone(),
            SELECT,
            MAX_DEPTH,
        )
        .is_err()
    );
    for value in [
        surface(&[
            &[(0.0, 300.0), (0.0008, 250.0)],
            &[(0.0012, 180.0), (MAX_DEPTH, 100.0)],
        ]),
        surface(&[&[(0.0, 300.0), (0.001, -1000.0), (MAX_DEPTH, 100.0)]]),
    ] {
        let (other, other_key) = fixture(value);
        let other_native = plastic(&other);
        assert!(matches!(
            ResolvedDepthHardness::from_j2_accepted(
                &other_native,
                &[0.0; 6],
                &other_native.initial_state(),
                GATES.yield_pa,
                &other,
                other_key,
                SELECT,
                MAX_DEPTH,
            ),
            Err(DepthHardnessHistoryError::Profile(
                DepthHardnessError::Support(_) | DepthHardnessError::Invalid(_)
            ))
        ));
        assert!(
            ResolvedDepthHardness::from_j2_accepted(
                &native,
                &[0.0; 6],
                &state,
                GATES.yield_pa,
                &other,
                key.clone(),
                SELECT,
                MAX_DEPTH,
            )
            .is_err()
        );
    }
}

#[test]
fn g3_full_cyclic_backstresses_are_admitted_before_the_depth_profile() {
    let (card, key) = ordinary_fixture();
    let law = resolve_chaboche_branches(
        &card,
        &query(),
        SELECT,
        HardeningPropertyBasis::MixedPlasticStrainModuli,
        0.05,
        [
            BackstressProperties::new("fast_C", "fast_gamma"),
            BackstressProperties::new("slow_C", "slow_gamma"),
        ],
    )
    .unwrap();
    let mut strain = [0.0; 6];
    let mut state = law.initial_state();
    for axial in [0.006, -0.006] {
        let next_strain = [axial, 0.0, 0.0, 0.0, 0.0, 0.0];
        let native = law
            .law()
            .advance_with_energy(&strain, &state, &next_strain, GATES)
            .unwrap();
        assert_ne!(native.state.backstress_pa[0], native.state.backstress_pa[1]);
        let profile = ResolvedDepthHardness::from_cyclic_accepted(
            &law,
            &next_strain,
            &native.state,
            GATES.yield_pa,
            &card,
            key.clone(),
            SELECT,
            MAX_DEPTH,
        )
        .unwrap();
        near(
            profile.hardness_pa_at(0.001).unwrap(),
            180.0 + SLOPE_ALPHA * native.state.plastic.alpha,
        );
        let mut invalid = native.state.clone();
        invalid.backstress_pa[1][0] += 1.0;
        assert!(matches!(
            ResolvedDepthHardness::from_cyclic_accepted(
                &law,
                &next_strain,
                &invalid,
                GATES.yield_pa,
                &card,
                key.clone(),
                SELECT,
                MAX_DEPTH,
            ),
            Err(DepthHardnessHistoryError::NativeCyclic(_))
        ));
        strain = next_strain;
        state = native.state;
    }
}
