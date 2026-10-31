//! Independent K/G invariants and source-domain checks; synthetic, not measured.
use super::*;
use crate::state_point::{
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, QueryPoint, UncertaintyModel,
};
use fs_qty::{
    Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

const K: f64 = 1200.0;
const KT: f64 = -1.5;
const G: f64 = 400.0;
const GT: f64 = 0.8;
const ALPHA: f64 = 0.0001;
const ALPHAT: f64 = 0.000001;

fn quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new().with_quantity("T", quantity(), t).unwrap()
}
fn card(nonaffine: bool) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value, slope) in [
        (BULK_MODULUS_PROPERTY, Pressure::DIMS, K, KT),
        (SHEAR_MODULUS_PROPERTY, Pressure::DIMS, G, GT),
        (
            LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
            Dims([0, 0, 0, -1, 0, 0]),
            ALPHA,
            ALPHAT,
        ),
        (
            ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES[0],
            Dims([0, 0, 0, -1, 0, 0]),
            ALPHA,
            ALPHAT,
        ),
        (
            ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES[1],
            Dims([0, 0, 0, -1, 0, 0]),
            -0.5 * ALPHA,
            0.2 * ALPHAT,
        ),
        (
            ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES[2],
            Dims([0, 0, 0, -1, 0, 0]),
            2.0 * ALPHA,
            -0.5 * ALPHAT,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![
                        (280.0, value - 20.0 * slope),
                        (
                            300.0,
                            value
                                + if nonaffine && name == BULK_MODULUS_PROPERTY {
                                    1.0
                                } else {
                                    0.0
                                },
                        ),
                        (360.0, value + 60.0 * slope),
                    ],
                    dims,
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    quantity(),
                    280.0,
                    360.0,
                ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "authored K/G thermoelastic fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "independent bulk/shear fixture".into(),
            phase: "solid".into(),
            process: "numerical derivative test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn sources(
    card: &MaterialCard,
    t: f64,
    orthotropic: bool,
) -> (ResolvedMaterialStatePoint, ResolvedMaterialStatePoint) {
    let requirements = [BULK_MODULUS_PROPERTY, SHEAR_MODULUS_PROPERTY].map(|name| {
        ScalarPropertyRequirement::try_new(
            name,
            Pressure::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap()
    });
    let names = if orthotropic {
        ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES.to_vec()
    } else {
        vec![LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY]
    };
    let expansion: Vec<_> = names
        .into_iter()
        .map(|name| {
            ScalarPropertyRequirement::try_new(
                name,
                Dims([0, 0, 0, -1, 0, 0]),
                ScalarAdmissibility::Finite,
            )
            .unwrap()
        })
        .collect();
    (
        resolve_material_state_point(
            card,
            &point(t),
            &requirements,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap(),
        resolve_material_state_point(
            card,
            &point(t),
            &expansion,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap(),
    )
}
fn law(card: &MaterialCard, expansion: ThermoelasticSymmetry) -> AffineThermoelastic {
    let (elastic, thermal) = sources(card, 300.0, expansion != ThermoelasticSymmetry::Isotropic);
    AffineThermoelastic::new(
        card,
        &elastic,
        &thermal,
        ThermoelasticSymmetry::IsotropicBulkShear,
        expansion,
        [280.0, 360.0],
        0.1,
    )
    .unwrap()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
}

// Trace/deviatoric invariants, independent of the production Mandel matrix.
fn oracle(t: f64, total: [f64; 6], volume: f64) -> (f64, f64, f64, [f64; 6], [f64; 6]) {
    let dt = t - 300.0;
    let (k, g) = (K + KT * dt, G + GT * dt);
    let free = ALPHA * dt + 0.5 * ALPHAT * dt * dt;
    let alpha = ALPHA + ALPHAT * dt;
    let total_trace = total[..3].iter().sum::<f64>();
    let trace = total_trace - 3.0 * free;
    let dev: [f64; 6] =
        core::array::from_fn(|i| total[i] - if i < 3 { total_trace / 3.0 } else { 0.0 });
    let dev2 = dev.iter().map(|x| x * x).sum::<f64>();
    let psi = volume * (0.5 * k * trace * trace + g * dev2);
    let entropy = volume * (3.0 * k * alpha * trace - 0.5 * KT * trace * trace - GT * dev2);
    let capacity = -t
        * volume
        * (-6.0 * KT * alpha * trace + 9.0 * k * alpha * alpha - 3.0 * k * ALPHAT * trace);
    let stress = core::array::from_fn(|i| 2.0 * g * dev[i] + if i < 3 { k * trace } else { 0.0 });
    let derivative = core::array::from_fn(|i| {
        2.0 * GT * dev[i]
            + if i < 3 {
                KT * trace - 3.0 * k * alpha
            } else {
                0.0
            }
    });
    (psi, entropy, capacity, stress, derivative)
}

#[test]
fn g1_bulk_shear_sources_drive_exact_thermal_potentials_and_derivatives() {
    let card = card(false);
    let law = law(&card, ThermoelasticSymmetry::Isotropic);
    let total = [0.013, -0.004, 0.002, 0.007, -0.003, 0.005];
    let volume = 0.7;
    for t in [281.0, 300.0, 321.0, 359.0] {
        let at = law.at_temperature(t).unwrap();
        let actual = at.state(total, volume).unwrap();
        let (psi, entropy, capacity, stress, derivative) = oracle(t, total, volume);
        close(actual.free_energy_j, psi, 1e-13);
        close(actual.entropy_j_k, entropy, 1e-14);
        close(actual.internal_energy_j, psi + t * entropy, 1e-12);
        close(actual.heat_capacity_j_k, capacity, 1e-13);
        for i in 0..6 {
            close(actual.stress_mandel_pa[i], stress[i], 2e-12);
            close(
                actual.stress_temperature_derivative_mandel_pa_k[i],
                derivative[i],
                2e-13,
            );
        }
        let h = 0.002;
        let low = oracle(t - h, total, volume);
        let high = oracle(t + h, total, volume);
        close(actual.entropy_j_k, -(high.0 - low.0) / (2.0 * h), 2e-10);
        close(
            actual.heat_capacity_j_k,
            t * (high.1 - low.1) / (2.0 * h),
            2e-8,
        );
        let relaxed = at.state(at.free_strain_mandel(), volume).unwrap();
        assert_eq!(relaxed.free_energy_j, 0.0);
        assert_eq!(relaxed.stress_mandel_pa, [0.0; 6]);
        assert!(at.elastic.property(YOUNG_MODULUS_PROPERTY).is_none());
        assert!(at.elastic.property(POISSON_RATIO_PROPERTY).is_none());
        let (elastic, thermal) = sources(&card, t, false);
        assert_eq!(at, law.at_materials(&elastic, &thermal).unwrap());
    }
}

#[test]
fn g3_bulk_shear_retains_independently_oriented_expansion_and_capacity_bound() {
    let card = card(false);
    let q = [[0.8, -0.6, 0.0], [0.6, 0.8, 0.0], [0.0, 0.0, 1.0]];
    let principal = law(
        &card,
        ThermoelasticSymmetry::OrientedOrthotropic {
            principal_to_target: IDENTITY,
        },
    );
    let rotated = law(
        &card,
        ThermoelasticSymmetry::OrientedOrthotropic {
            principal_to_target: q,
        },
    );
    let total = [0.013, -0.004, 0.002, 0.007, -0.003, 0.005];
    let transform = mandel_rotation(q);
    let rotate =
        |x: [f64; 6]| core::array::from_fn(|i| (0..6).map(|j| transform[i][j] * x[j]).sum::<f64>());
    let bound = rotated.minimum_heat_capacity_j_k(0.7).unwrap();
    for t in [281.0, 300.0, 321.0, 359.0] {
        let a = principal
            .at_temperature(t)
            .unwrap()
            .state(total, 0.7)
            .unwrap();
        let b = rotated
            .at_temperature(t)
            .unwrap()
            .state(rotate(total), 0.7)
            .unwrap();
        for (actual, expected) in [
            b.free_energy_j,
            b.entropy_j_k,
            b.internal_energy_j,
            b.heat_capacity_j_k,
        ]
        .into_iter()
        .zip([
            a.free_energy_j,
            a.entropy_j_k,
            a.internal_energy_j,
            a.heat_capacity_j_k,
        ]) {
            close(actual, expected, 3e-12);
        }
        for (actual, expected) in b
            .stress_mandel_pa
            .into_iter()
            .zip(rotate(a.stress_mandel_pa))
        {
            close(actual, expected, 3e-12);
        }
        assert!(b.heat_capacity_j_k >= bound);
    }
    for i in 0..=128 {
        let t = 280.0 + 80.0 * i as f64 / 128.0;
        for sign in [-1.0, 1.0] {
            let strain = [0.06 * sign; 6];
            assert!(
                rotated
                    .at_temperature(t)
                    .unwrap()
                    .state(strain, 0.7)
                    .unwrap()
                    .heat_capacity_j_k
                    >= bound
            );
        }
    }
}

#[test]
fn g0_bulk_shear_refuses_nonaffine_sources_wrong_expansion_schema_and_domain_exit() {
    let good = card(false);
    let law = law(&good, ThermoelasticSymmetry::Isotropic);
    assert!(law.at_temperature(360.01).is_err());
    let (elastic, thermal) = sources(&good, 300.0, false);
    assert!(
        AffineThermoelastic::new(
            &good,
            &elastic,
            &thermal,
            ThermoelasticSymmetry::IsotropicBulkShear,
            ThermoelasticSymmetry::IsotropicBulkShear,
            [280.0, 360.0],
            0.1
        )
        .is_err()
    );
    let bad = card(true);
    let (elastic, thermal) = sources(&bad, 300.0, false);
    assert!(
        AffineThermoelastic::new(
            &bad,
            &elastic,
            &thermal,
            ThermoelasticSymmetry::IsotropicBulkShear,
            ThermoelasticSymmetry::Isotropic,
            [280.0, 360.0],
            0.1
        )
        .is_err()
    );
    assert!(law.at_materials(&elastic, &thermal).is_err());
}
