//! G0/G1/G3: independent 3-D strength polynomial and original-source admission.
//! Manufactured coefficients verify implementation, not a physical material.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::MaterialPropertySelection;
use fs_material::strength::OrthotropicTsaiWu;
use fs_qty::{Dims, Pressure, Temperature};

const STRENGTH: [&str; 9] = [
    "x1t", "x1c", "x2t", "x2c", "x3t", "x3c", "s12", "s23", "s31",
];
const INTERACTION: [&str; 3] = ["f12", "f23", "f31"];
const SELECT: MaterialPropertySelection = MaterialPropertySelection::SingleClaimOnly;

fn card(
    strength: [f64; 9],
    interaction: [f64; 3],
    omit: Option<&str>,
    bad_units: Option<&str>,
    rate: bool,
) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (i, name) in STRENGTH.into_iter().chain(INTERACTION).enumerate() {
        if omit == Some(name) {
            continue;
        }
        let value = if i < 9 {
            strength[i]
        } else {
            interaction[i - 9]
        };
        let dims = if bad_units == Some(name) {
            Dims::NONE
        } else if i < 9 {
            Pressure::DIMS
        } else {
            Dims::NONE
        };
        let mut validity = ValidityDomain::unconstrained().with("T", 300.0, 500.0);
        if rate {
            validity = validity.with("strain_rate", 1.0, 2.0);
        }
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![
                        (300.0, value),
                        (500.0, if i < 9 { value * 0.5 } else { value }),
                    ],
                    dims,
                },
                validity,
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured 3-D Tsai-Wu fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "solid".into(),
            process: "declared orthotropic directions".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn point(t: f64) -> QueryPoint {
    QueryPoint::new().with("T", t).unwrap()
}
fn law(strength: [f64; 9], interaction: [f64; 3]) -> OrthotropicTsaiWu {
    OrthotropicTsaiWu::resolve(
        &card(strength, interaction, None, None, false),
        &point(300.0),
        STRENGTH,
        INTERACTION,
        SELECT,
    )
    .unwrap()
}
fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-13 * expected.abs().max(1.0),
        "{actual:e} != {expected:e}"
    );
}

// Deliberately evaluate the dimensional polynomial, independently of the
// production Cholesky factor and scaled radial root. Test magnitudes are small.
fn polynomial(s: [f64; 6], x: [f64; 9], f: [f64; 3]) -> f64 {
    let mut value = 0.0;
    for i in 0..3 {
        value += s[i] * (1.0 / x[2 * i] - 1.0 / x[2 * i + 1]);
        value += s[i] * s[i] / (x[2 * i] * x[2 * i + 1]);
        value += s[i + 3] * s[i + 3] / (x[i + 6] * x[i + 6]);
    }
    for (k, (i, j)) in [(0, 1), (1, 2), (2, 0)].into_iter().enumerate() {
        value +=
            2.0 * f[k] * s[i] * s[j] / (x[2 * i] * x[2 * i + 1] * x[2 * j] * x[2 * j + 1]).sqrt();
    }
    value
}

#[test]
fn g1_every_original_tension_compression_and_shear_intercept_is_recovered() {
    let strengths = [100.0, 180.0, 15.0, 35.0, 8.0, 19.0, 6.0, 4.0, 5.0];
    let law = law(strengths, [0.2, -0.3, 0.1]);
    for axis in 0..3 {
        for sign in [1.0, -1.0] {
            let mut stress = [0.0; 6];
            stress[axis] = sign * strengths[2 * axis + usize::from(sign < 0.0)];
            near(
                law.assess_material_stress(&stress).unwrap().utilization(),
                1.0,
            );
            assert!(
                !law.assess_material_stress(&stress.map(|v| v * 0.99))
                    .unwrap()
                    .reaches_limit()
            );
            assert!(
                law.assess_material_stress(&stress.map(|v| v * 1.01))
                    .unwrap()
                    .reaches_limit()
            );
            stress = [0.0; 6];
            stress[axis + 3] = sign * strengths[axis + 6];
            near(
                law.assess_material_stress(&stress).unwrap().utilization(),
                1.0,
            );
        }
    }
}

#[test]
fn g1_complete_stress_surface_matches_independent_dimensional_polynomial() {
    let strengths = [100.0, 180.0, 15.0, 35.0, 8.0, 19.0, 6.0, 4.0, 5.0];
    let interaction = [0.2, -0.3, 0.1];
    let law = law(strengths, interaction);
    for stress in [
        [31.0, 4.0, 3.0, -2.0, 1.0, 2.0],
        [-55.0, -17.0, -4.0, 1.0, 0.5, -3.0],
        [83.0, -14.0, 7.0, 4.0, -3.0, 1.0],
        [0.0, 0.0, 0.0, -3.0, 2.0, 4.0],
    ] {
        let assessment = law.assess_material_stress(&stress).unwrap();
        let u = assessment.utilization();
        assert_eq!(assessment.material_stress_pa(), stress);
        near(
            polynomial(stress.map(|v| v / u), strengths, interaction),
            1.0,
        );
        assert_eq!(
            assessment.reaches_limit(),
            polynomial(stress, strengths, interaction) >= 1.0
        );
        near(
            law.assess_material_stress(&stress.map(|v| 7.0 * v))
                .unwrap()
                .utilization(),
            7.0 * u,
        );
    }
    let zero = law.assess_material_stress(&[0.0; 6]).unwrap();
    assert_eq!(zero.utilization(), 0.0);
    assert!(!zero.reaches_limit());
}

#[test]
fn g1_biaxial_tension_and_out_of_plane_shear_can_govern_before_uniaxial_limits() {
    let law = law(
        [100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 50.0, 50.0, 50.0],
        [0.0; 3],
    );
    for stress in [
        [80.0, 80.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 80.0, 0.0, 40.0, 0.0],
    ] {
        let result = law.assess_material_stress(&stress).unwrap();
        near(result.utilization(), 1.28_f64.sqrt());
        assert!(result.reaches_limit());
    }
    let positive = law
        .assess_material_stress(&[80.0, 0.0, 0.0, 20.0, -10.0, 5.0])
        .unwrap();
    let reverse_shear = law
        .assess_material_stress(&[80.0, 0.0, 0.0, -20.0, 10.0, -5.0])
        .unwrap();
    assert_eq!(positive.utilization(), reverse_shear.utilization());
}

#[test]
fn g0_full_convexity_is_checked_instead_of_three_independent_pair_bounds() {
    for interactions in [[-0.75; 3], [-0.5; 3], [0.9, 0.9, -0.9], [1.0, 0.0, 0.0]] {
        let card = card([100.0; 9], interactions, None, None, false);
        assert!(
            OrthotropicTsaiWu::resolve(&card, &point(300.0), STRENGTH, INTERACTION, SELECT)
                .is_err()
        );
    }
    // Nonzero, unequal interactions in all three planes are admitted together.
    let law = law([100.0; 9], [0.2, -0.3, 0.1]);
    assert_eq!(law.normalized_interactions(), [0.2, -0.3, 0.1]);
}

#[test]
fn g3_extreme_units_and_compressive_asymmetry_keep_representable_results() {
    let strengths = [100.0, 180.0, 15.0, 35.0, 8.0, 19.0, 6.0, 4.0, 5.0];
    let stress = [-51.0, -9.0, 4.0, 2.0, -1.0, 3.0];
    let reference = law(strengths, [0.2, -0.3, 0.1])
        .assess_material_stress(&stress)
        .unwrap()
        .utilization();
    for scale in [1e-280, 1e280] {
        near(
            law(strengths.map(|v| v * scale), [0.2, -0.3, 0.1])
                .assess_material_stress(&stress.map(|v| v * scale))
                .unwrap()
                .utilization(),
            reference,
        );
    }
    // Naive subtraction of sqrt(B^2+4A)+B loses the positive root here.
    let asymmetric = law([1e-100, 1e100, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], [0.0; 3]);
    near(
        asymmetric
            .assess_material_stress(&[-1e100, 0.0, 0.0, 0.0, 0.0, 0.0])
            .unwrap()
            .utilization(),
        1.0,
    );
}

#[test]
fn g0_refresh_retains_all_original_sources_and_refuses_unowned_or_missing_data() {
    let strengths = [100.0; 9];
    let original = card(strengths, [0.2, -0.3, 0.1], None, None, false);
    let cold = OrthotropicTsaiWu::resolve(&original, &point(300.0), STRENGTH, INTERACTION, SELECT)
        .unwrap();
    let warm = cold.at_query_point(&original, &point(400.0)).unwrap();
    assert_eq!(warm.strengths_pa(), [75.0; 9]);
    assert_eq!(
        warm.normalized_interactions(),
        cold.normalized_interactions()
    );
    assert_eq!(warm.strength_properties(), STRENGTH);
    assert_eq!(warm.interaction_properties(), INTERACTION);
    assert_eq!(
        warm.resolved().card_identity(),
        cold.resolved().card_identity()
    );
    assert_ne!(warm.resolved().identity(), cold.resolved().identity());
    assert!(cold.at_query_point(&original, &point(501.0)).is_err());
    assert!(
        cold.at_query_point(
            &card([90.0; 9], [0.2, -0.3, 0.1], None, None, false),
            &point(300.0)
        )
        .is_err()
    );
    for missing in ["x3c", "s23", "f31"] {
        assert!(
            OrthotropicTsaiWu::resolve(
                &card(strengths, [0.0; 3], Some(missing), None, false),
                &point(300.0),
                STRENGTH,
                INTERACTION,
                SELECT
            )
            .is_err()
        );
    }
    assert!(
        OrthotropicTsaiWu::resolve(
            &card(strengths, [0.0; 3], None, Some("x2t"), false),
            &point(300.0),
            STRENGTH,
            INTERACTION,
            SELECT
        )
        .is_err()
    );
    let rate_point = point(300.0).with("strain_rate", 1.5).unwrap();
    assert!(
        OrthotropicTsaiWu::resolve(
            &card(strengths, [0.0; 3], None, None, true),
            &rate_point,
            STRENGTH,
            INTERACTION,
            SELECT
        )
        .is_err()
    );
    let mut duplicate = STRENGTH;
    duplicate[1] = duplicate[0];
    assert!(
        OrthotropicTsaiWu::try_from_resolved(cold.resolved().clone(), duplicate, INTERACTION)
            .is_err()
    );
}

#[test]
fn g0_nonfinite_stress_and_unrepresentable_ratios_refuse() {
    let ordinary = law([1.0; 9], [0.0; 3]);
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        for i in 0..6 {
            let mut stress = [0.0; 6];
            stress[i] = value;
            assert!(ordinary.assess_material_stress(&stress).is_err());
        }
    }
    let tiny = law([1e-300; 9], [0.0; 3]);
    assert!(
        tiny.assess_material_stress(&[1e300, 0.0, 0.0, 0.0, 0.0, 0.0])
            .is_err()
    );
    let huge = law([1e300; 9], [0.0; 3]);
    assert!(
        huge.assess_material_stress(&[1e-300, 0.0, 0.0, 0.0, 0.0, 0.0])
            .is_err()
    );
}
