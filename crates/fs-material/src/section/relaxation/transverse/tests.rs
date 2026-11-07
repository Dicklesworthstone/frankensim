//! G0/G1 source admission and independent annular transverse-shear response.

use super::*;
use crate::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use crate::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, YOUNG_MODULUS_PROPERTY, resolve_material_state_point,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_qty::{Density, Pressure, Time};

fn material(
    shear: f64,
    equilibrium: f64,
    arms: &[(&str, &str, f64, f64)],
    band: Option<(f64, f64)>,
    shear_band: Option<(f64, f64)>,
) -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    let mut values = vec![
        (DENSITY_PROPERTY, Density::DIMS, 7800.0),
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, 2e9),
        ("shear", Pressure::DIMS, shear),
        ("equilibrium", Pressure::DIMS, equilibrium),
    ];
    for &(delta, time, value, tau) in arms {
        values.extend([(delta, Pressure::DIMS, value), (time, Time::DIMS, tau)]);
    }
    for (name, dims, value) in values {
        let support = if name == "shear" {
            shear_band.or(band)
        } else {
            band
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: support.map_or_else(ValidityDomain::unconstrained, |(lo, hi)| {
                    ValidityDomain::unconstrained().with("omega", lo, hi)
                }),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic layered transverse relaxation fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "solid".into(),
            process: "shear spectrum test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let (lo, hi) = match (band, shear_band) {
        (Some(a), Some(b)) => (a.0.max(b.0), a.1.min(b.1)),
        (Some(a), None) | (None, Some(a)) => a,
        (None, None) => (0.0, 200.0),
    };
    let query = QueryPoint::new().with("omega", 0.5 * (lo + hi)).unwrap();
    resolve_material_state_point(
        &card,
        &query,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn selection(pairs: &[(&str, &str)]) -> LayerRelaxationProperties {
    LayerRelaxationProperties {
        equilibrium_modulus: "equilibrium".into(),
        branches: pairs
            .iter()
            .map(|(a, b)| ((*a).into(), (*b).into()))
            .collect(),
    }
}

fn section(a: ResolvedMaterialStatePoint, b: ResolvedMaterialStatePoint) -> ConcentricShearSection {
    ConcentricElasticSection::try_new(
        0.0004,
        vec![
            ConcentricMaterialLayer {
                outer_radius_m: 0.001,
                material: a,
            },
            ConcentricMaterialLayer {
                outer_radius_m: 0.002,
                material: b,
            },
        ],
    )
    .unwrap()
    .with_shear_modulus("shear", 0.71)
    .unwrap()
}

fn near(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 3e-13 * a.abs().max(b.abs()).max(1e-30),
        "{a} != {b}"
    );
}

#[test]
fn g1_corrected_layer_areas_keep_each_original_shear_spectrum() {
    let core = material(
        0.8e9,
        0.8e9,
        &[("core_delta", "core_time", 1.3e9, 0.002)],
        Some((0.0, 2000.0)),
        None,
    );
    let sleeve = material(
        1.4e9,
        1.4e9,
        &[
            ("sleeve_delta", "sleeve_time", 0.5e9, 0.025),
            ("zero_delta", "zero_time", 0.0, 0.9),
        ],
        Some((10.0, 1000.0)),
        None,
    );
    let shear = section(core, sleeve);
    let properties = [
        selection(&[("core_delta", "core_time")]),
        selection(&[("sleeve_delta", "sleeve_time"), ("zero_delta", "zero_time")]),
    ];
    let law = shear.clone().with_prony_transverse(&properties).unwrap();
    assert_eq!(law.shear(), &shear);
    assert_eq!(law.properties(), properties);
    assert_eq!(law.omega_band_rad_s(), (10.0, 1000.0));
    assert_eq!(law.spectra()[0].terms, vec![(1.3e9, 0.002)]);
    assert_eq!(law.spectra()[1].terms, vec![(0.5e9, 0.025), (0.0, 0.9)]);
    assert_eq!(
        law.branches()
            .iter()
            .map(|b| b.layer_index)
            .collect::<Vec<_>>(),
        [0, 1, 1]
    );
    assert_eq!(
        law.branches()
            .iter()
            .map(|b| b.relaxation_time_s)
            .collect::<Vec<_>>(),
        [0.002, 0.025, 0.9]
    );
    let areas = [
        core::f64::consts::PI * (0.001_f64.powi(2) - 0.0004_f64.powi(2)),
        core::f64::consts::PI * (0.002_f64.powi(2) - 0.001_f64.powi(2)),
    ];
    near(
        law.branches()[0].relaxing_rigidity_n,
        0.71 * areas[0] * 1.3e9,
    );
    near(
        law.branches()[1].relaxing_rigidity_n,
        0.71 * areas[1] * 0.5e9,
    );
    assert_eq!(
        law.branches()[2].relaxing_rigidity_n.to_bits(),
        0.0_f64.to_bits()
    );
    let aggregate = GeneralizedMaxwell::new(
        shear.shear_stiffness_n(),
        law.branches()
            .iter()
            .map(|b| (b.relaxing_rigidity_n, b.relaxation_time_s))
            .collect(),
    )
    .unwrap();
    for omega in [10.0, 50.0, 1000.0] {
        let mut storage = 0.71 * (areas[0] * 0.8e9 + areas[1] * 1.4e9);
        let mut loss = 0.0;
        for (rigidity, tau) in [
            (0.71 * areas[0] * 1.3e9, 0.002),
            (0.71 * areas[1] * 0.5e9, 0.025),
        ] {
            let x = omega * tau;
            storage += rigidity * x * x / (1.0 + x * x);
            loss += rigidity * x / (1.0 + x * x);
        }
        let actual = aggregate.modulus(omega);
        near(actual.0, storage);
        near(actual.1, loss);
    }
}

#[test]
fn g0_explicit_elastic_and_zero_arms_retain_sources_but_missing_or_wrong_g_refuses() {
    let a = material(0.8e9, 0.8e9, &[], Some((0.0, 2000.0)), None);
    let b = material(
        1.4e9,
        1.4e9,
        &[("delta", "time", 0.0, 0.003)],
        Some((0.0, 1000.0)),
        None,
    );
    let original = section(a.clone(), b.clone());
    let plain = selection(&[]);
    let zero = selection(&[("delta", "time")]);
    let law = original
        .clone()
        .with_prony_transverse(&[plain.clone(), zero.clone()])
        .unwrap();
    assert_eq!(law.shear(), &original);
    assert!(law.spectra()[0].terms.is_empty());
    assert_eq!(law.spectra()[1].terms, vec![(0.0, 0.003)]);
    assert_eq!(law.branches().len(), 1);
    assert_eq!(law.branches()[0].layer_index, 1);
    assert_eq!(law.branches()[0].relaxing_rigidity_n, 0.0);
    assert!(
        original
            .clone()
            .with_prony_transverse(&[plain.clone()])
            .is_err()
    );
    assert!(
        original
            .clone()
            .with_prony_transverse(&[zero.clone(), zero])
            .is_err()
    );
    assert!(
        original
            .with_prony_transverse(&[plain.clone(), selection(&[("time", "delta")])])
            .is_err()
    );
    // Young's E cannot stand in for the selected long-time shear modulus G.
    let wrong = material(0.8e9, 2e9, &[], Some((0.0, 2000.0)), None);
    assert!(
        section(wrong, b)
            .with_prony_transverse(&[plain.clone(), selection(&[("delta", "time")])])
            .is_err()
    );
    for (delta, tau) in [(-1.0, 0.003), (1.0, 0.0)] {
        let invalid = material(
            1.4e9,
            1.4e9,
            &[("delta", "time", delta, tau)],
            Some((0.0, 1000.0)),
            None,
        );
        assert!(
            section(a.clone(), invalid)
                .with_prony_transverse(&[plain.clone(), selection(&[("delta", "time")])])
                .is_err()
        );
    }
}

#[test]
fn g0_frequency_support_intersects_actual_shear_parents_and_rejects_missing_or_disjoint_bands() {
    let a = material(
        0.8e9,
        0.8e9,
        &[("delta", "time", 1e9, 0.002)],
        Some((0.0, 2000.0)),
        Some((30.0, 700.0)),
    );
    let b = material(1.4e9, 1.4e9, &[], Some((20.0, 1000.0)), None);
    let binding = [selection(&[("delta", "time")]), selection(&[])];
    assert_eq!(
        section(a.clone(), b)
            .with_prony_transverse(&binding)
            .unwrap()
            .omega_band_rad_s(),
        (30.0, 700.0)
    );
    let disjoint = material(1.4e9, 1.4e9, &[], Some((800.0, 1000.0)), None);
    assert!(
        section(a.clone(), disjoint)
            .with_prony_transverse(&binding)
            .is_err()
    );
    // A band's presence on G does not supply a missing explicit law/time band.
    let missing_empty = material(1.4e9, 1.4e9, &[], None, Some((20.0, 1000.0)));
    assert!(
        section(a.clone(), missing_empty)
            .with_prony_transverse(&binding)
            .is_err()
    );
    let missing_time = material(
        1.4e9,
        1.4e9,
        &[("delta", "time", 1e9, 0.003)],
        None,
        Some((20.0, 1000.0)),
    );
    assert!(
        section(a, missing_time)
            .with_prony_transverse(&[binding[0].clone(), binding[0].clone()])
            .is_err()
    );
}
