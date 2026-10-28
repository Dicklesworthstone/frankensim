//! G0/G1/G3 source-to-specimen tests. Synthetic fixtures verify mechanics and
//! source binding; they are not experimental validation of a real foam.
#![deny(unsafe_code)]

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::visco::source::{
    PronyReferenceLimits, PronyReferenceRequest, ResolvedPronyReference,
};
use fs_qty::{Density, Dims, Pressure, QuantitySpec};

const TIME: Dims = Dims([0, 0, 1, 0, 0, 0]);
const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);
const AREA: Dims = Dims([2, 0, 0, 0, 0, 0]);
const SINGLE: SelectionPolicy = SelectionPolicy::SingleClaimOnly;

fn context() -> ValidityDomain {
    ValidityDomain::unconstrained()
        .with_quantity("sample", QuantitySpec::dimensional(Dims::NONE), 1.0, 1.0)
        .with_quantity(
            "original-length",
            QuantitySpec::dimensional(LENGTH),
            0.02,
            0.02,
        )
        .with_quantity(
            "original-area",
            QuantitySpec::dimensional(AREA),
            0.001,
            0.001,
        )
}
fn claim(name: &str, dims: Dims, value: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: context(),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic standard linear solid; numerical test only".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}
fn rows() -> Vec<PropertyClaim> {
    vec![
        claim("prony-branch-count", Dims::NONE, 1.0),
        claim("prony-equilibrium-modulus", Pressure::DIMS, 2e6),
        claim("prony-modulus-1", Pressure::DIMS, 3e6),
        claim("prony-relaxation-time-1", TIME, 4.0),
    ]
}
fn card(rows: Vec<PropertyClaim>) -> MaterialCard {
    let mut set = ClaimSet::new();
    for row in rows {
        set.insert_claim(row).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic SLS".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        set,
        vec![],
    )
    .unwrap()
}
fn point() -> QueryPoint {
    QueryPoint::new()
        .with_quantity("sample", QuantitySpec::dimensional(Dims::NONE), 1.0)
        .unwrap()
        .with_quantity("original-length", QuantitySpec::dimensional(LENGTH), 0.02)
        .unwrap()
        .with_quantity("original-area", QuantitySpec::dimensional(AREA), 0.001)
        .unwrap()
}
fn request() -> PronyReferenceRequest {
    PronyReferenceRequest {
        limits: PronyReferenceLimits {
            minimum_incremental_strain: -0.1,
            maximum_incremental_strain: 0.1,
            maximum_absolute_strain_rate_s: 1.0,
            maximum_elapsed_s: 1000.0,
        },
        original_length_axis: Some("original-length".into()),
        original_area_axis: Some("original-area".into()),
    }
}
fn resolve(source: &MaterialCard) -> ResolvedPronyReference {
    ResolvedPronyReference::resolve(source, &point(), &request(), SINGLE).unwrap()
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 5e-11 * expected.abs().max(1e-8),
        "{actual:.17e} != {expected:.17e}"
    );
}
fn balanced(r: &fs_material::visco::source::PronyReferenceResponse) {
    assert!(
        r.energy_balance_residual_j.abs()
            <= 5e-11
                * (r.model_work_j.abs() + r.stored_energy_j + r.dissipated_energy_j).max(1e-12)
    );
    assert!(r.incremental_dissipated_energy_j >= 0.0);
}

#[test]
fn g1_strain_hold_relaxes_with_exact_source_time_and_releases_storage() {
    let source = card(rows());
    let mut coupon = resolve(&source).coupon(0.02, 0.001).unwrap();
    let ramp = coupon.advance_strain(0.02, 1.0).unwrap();
    let q = 0.02 * 4.0 * (1.0 - (-0.25_f64).exp());
    close(ramp.incremental_stress_pa, 2e6 * 0.02 + 3e6 * q);
    close(ramp.incremental_force_n, ramp.incremental_stress_pa * 0.001);
    close(ramp.incremental_displacement_m, 0.0004);
    let hold = coupon.advance_strain(0.02, 10.0).unwrap();
    close(
        hold.incremental_stress_pa,
        2e6 * 0.02 + 3e6 * q * (-2.5_f64).exp(),
    );
    close(hold.incremental_model_work_j, 0.0);
    close(
        hold.incremental_dissipated_energy_j,
        ramp.stored_energy_j - hold.stored_energy_j,
    );
    balanced(&ramp);
    balanced(&hold);
}

#[test]
fn g1_force_ramp_hold_and_zero_load_recovery_match_sls_creep() {
    let mut coupon = resolve(&card(rows())).coupon(0.02, 0.001).unwrap();
    let ramp = coupon.advance_force(10.0, 2.0).unwrap();
    // J(t)=1/Einf - E1/(Einf*E0)*exp(-t/Tr), Tr=tau*E0/Einf=10s.
    let equilibrium = 10000.0 / 2e6;
    let strain = equilibrium - (3e6 / (2e6 * 5e6)) * 5000.0 * 10.0 * (1.0 - (-0.2_f64).exp());
    close(ramp.incremental_strain, strain);
    let hold = coupon.advance_force(10.0, 30.0).unwrap();
    close(
        hold.incremental_strain,
        equilibrium + (strain - equilibrium) * (-3.0_f64).exp(),
    );
    assert!(hold.incremental_strain > ramp.incremental_strain);
    let unload = coupon.advance_force(0.0, 2.0).unwrap();
    let recover = coupon.advance_force(0.0, 30.0).unwrap();
    close(
        recover.incremental_strain,
        unload.incremental_strain * (-3.0_f64).exp(),
    );
    close(recover.incremental_force_n, 0.0);
    close(recover.incremental_model_work_j, 0.0);
    for response in [&ramp, &hold, &unload, &recover] {
        balanced(response);
    }
}

#[test]
fn g3_control_switching_keeps_memory_energy_and_deterministic_history() {
    let law = resolve(&card(rows()));
    let mut a = law.coupon(0.02, 0.001).unwrap();
    let mut b = a.clone();
    for coupon in [&mut a, &mut b] {
        let load = coupon.advance_strain(0.03, 0.5).unwrap();
        let creep = coupon.advance_force(load.incremental_force_n, 3.0).unwrap();
        assert!(creep.incremental_strain > load.incremental_strain);
        let held = coupon
            .advance_strain(creep.incremental_strain, 2.0)
            .unwrap();
        assert!(held.incremental_force_n < creep.incremental_force_n);
        let unload = coupon.advance_strain(0.0, 0.5).unwrap();
        assert!(unload.incremental_model_work_j < 0.0);
        assert!(unload.stored_energy_j > 0.0); // Loop work is not all dissipative heat.
        balanced(&unload);
        coupon.advance_force(0.0, 1.0).unwrap();
        let recovered = coupon.advance_force(0.0, 100.0).unwrap();
        balanced(&recovered);
    }
    assert_eq!(a.identity(), b.identity());
    assert_eq!(a.last_response(), b.last_response());
    let fresh = law.coupon(0.02, 0.001).unwrap();
    assert_ne!(fresh.identity(), a.identity());
}

#[test]
fn g0_original_receipts_uncertainty_geometry_and_missing_density_are_preserved() {
    let source = card(rows());
    let law = resolve(&source);
    assert_eq!(law.card_identity(), source.content_hash());
    assert_eq!(law.properties().len(), 4);
    assert_eq!(law.branches(), [(3e6, 4.0)]);
    for (_, answer) in law.properties() {
        source.claims().verify_receipt(&answer.receipt).unwrap();
        assert_eq!(
            answer.evidence.value.uncertainty,
            UncertaintyModel::Unstated
        );
    }
    assert_eq!(law.coupon(0.02, 0.001).unwrap().mass_kg(), None);
    assert!(law.coupon(0.015, 0.001).is_err());
    assert!(law.coupon(0.02, 0.002).is_err());
    let mut with_density = rows();
    with_density.push(claim("density", Density::DIMS, 65.0));
    close(
        resolve(&card(with_density))
            .coupon(0.02, 0.001)
            .unwrap()
            .mass_kg()
            .unwrap(),
        0.0013,
    );
    let extra = point().with("invented-rate", 100.0).unwrap();
    assert!(ResolvedPronyReference::resolve(&source, &extra, &request(), SINGLE).is_err());
    let mut bad = request();
    bad.original_area_axis = Some("original-length".into());
    assert!(ResolvedPronyReference::resolve(&source, &point(), &bad, SINGLE).is_err());
}

#[test]
fn g0_bad_coefficients_and_source_conditions_never_create_a_law() {
    let resolve_bad =
        |rows| ResolvedPronyReference::resolve(&card(rows), &point(), &request(), SINGLE);
    let mut missing = rows();
    missing.pop();
    assert!(resolve_bad(missing).is_err());
    let mut fractional = rows();
    fractional[0] = claim("prony-branch-count", Dims::NONE, 1.5);
    assert!(resolve_bad(fractional).is_err());
    let mut zero_arm = rows();
    zero_arm[2] = claim("prony-modulus-1", Pressure::DIMS, 0.0);
    assert!(resolve_bad(zero_arm).is_err());
    let mut wrong_units = rows();
    wrong_units[3] = claim("prony-relaxation-time-1", Pressure::DIMS, 4.0);
    assert!(resolve_bad(wrong_units).is_err());
    let mut curve = rows();
    curve[2].value = PropertyValue::Curve {
        abscissa: "sample".into(),
        abscissa_dims: Dims::NONE,
        knots: vec![(1.0, 3e6), (2.0, 3e6)],
        dims: Pressure::DIMS,
    };
    curve[2].interpolation = InterpolationPolicy::TabulatedOnly;
    assert!(resolve_bad(curve).is_err());
    let source = card(rows());
    assert!(
        ResolvedPronyReference::resolve(
            &source,
            &point().with("sample", 2.0).unwrap(),
            &request(),
            SINGLE
        )
        .is_err()
    );
}

#[test]
fn g0_history_limit_refusals_are_atomic_and_replay_stays_usable() {
    let source = card(rows());
    let mut scope = request();
    scope.limits.maximum_absolute_strain_rate_s = 0.01;
    scope.limits.maximum_elapsed_s = 10.0;
    let mut coupon = ResolvedPronyReference::resolve(&source, &point(), &scope, SINGLE)
        .unwrap()
        .coupon(0.02, 0.001)
        .unwrap();
    coupon.advance_strain(0.01, 2.0).unwrap();
    let checkpoint = coupon.clone();
    for (strain, dt) in [
        (0.02, 0.1),
        (0.2, 5.0),
        (0.01, 9.0),
        (0.01, 0.0),
        (f64::NAN, 1.0),
        (0.0, f64::INFINITY),
    ] {
        assert!(coupon.advance_strain(strain, dt).is_err());
        assert_eq!(coupon.last_response(), checkpoint.last_response());
        assert_eq!(coupon.identity(), checkpoint.identity());
    }
    for force in [1e6, f64::INFINITY, f64::NAN] {
        assert!(coupon.advance_force(force, 1.0).is_err());
        assert_eq!(coupon.identity(), checkpoint.identity());
    }
    let mut replay = checkpoint;
    assert_eq!(
        coupon.advance_force(0.0, 3.0).unwrap(),
        replay.advance_force(0.0, 3.0).unwrap()
    );
}

#[test]
fn g1_empty_spectrum_is_elastic_under_both_controls() {
    let source = card(vec![
        claim("prony-branch-count", Dims::NONE, 0.0),
        claim("prony-equilibrium-modulus", Pressure::DIMS, 2e6),
    ]);
    let mut coupon = resolve(&source).coupon(0.02, 0.001).unwrap();
    let load = coupon.advance_force(20.0, 2.0).unwrap();
    close(load.incremental_strain, 0.01);
    close(load.model_work_j, 0.5 * 2e6 * 0.01 * 0.01 * 0.02 * 0.001);
    close(load.dissipated_energy_j, 0.0);
    let unload = coupon.advance_strain(0.0, 2.0).unwrap();
    close(unload.stored_energy_j, 0.0);
    close(unload.model_work_j, 0.0);
    close(unload.dissipated_energy_j, 0.0);
}
