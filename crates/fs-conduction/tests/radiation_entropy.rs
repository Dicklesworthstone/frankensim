//! Original material-card/ambient-law integration, not experimental validation.
use fs_conduction::{AmbientRadiationPatch, SurfaceEmissivity,
    EMISSIVITY_DIMS, SURFACE_EMISSIVITY_PROPERTY, STEFAN_BOLTZMANN_W_M2_K4};
use fs_conduction::radiation_entropy::{GrayPhotonEmission, PhotonEntropyQuadrature};
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, SelectionPolicy, UncertaintyModel};

fn patch(emissivity: f64, ambient: f64) -> AmbientRadiationPatch {
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
        value: PropertyValue::Scalar { value: emissivity, dims: EMISSIVITY_DIMS },
        validity: ValidityDomain::unconstrained().with("T", 250.0, 500.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: Vec::new(),
        provenance: Provenance { source: "analytic entropy fixture".into(),
            license: "internal-test-use".into(), artifact: None },
    }).unwrap();
    let card = MaterialCard::assemble(MaterialStateId { chemistry: "declared gray surface".into(),
        phase: "solid".into(), process: "test".into(), revision: 0 }, claims, Vec::new()).unwrap();
    let source = SurfaceEmissivity::from_card("wall", &card, 350.0,
        SelectionPolicy::SingleClaimOnly).unwrap();
    AmbientRadiationPatch::new("wall", source, ambient).unwrap()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a-b).abs() <= tolerance, "{a:e} != {b:e}");
}

#[test]
fn original_signed_heat_closes_surface_and_bath_entropy_in_both_directions() {
    for e in [0.1, 0.5, 1.0] {
        for (temperature, bath) in [(400.0_f64, 300.0_f64), (300.0, 400.0)] {
            let p = patch(e, bath);
            let saved = p.clone();
            let a = p.material_entropy_exchange(temperature).unwrap();
            let expected_q = e*STEFAN_BOLTZMANN_W_M2_K4*(temperature.powi(4)-bath.powi(4));
            assert_eq!(a.heat_out_w_m2, p.heat_flux_w_m2(temperature).unwrap());
            close(a.heat_out_w_m2, expected_q, 2e-12);
            close(a.surface_storage_rate_w_m2_k, -expected_q/temperature, 1e-14);
            close(a.ambient_storage_rate_w_m2_k, expected_q/bath, 1e-14);
            close(a.generation_w_m2_k, expected_q*(1.0/bath-1.0/temperature), 1e-14);
            close(a.surface_storage_rate_w_m2_k+a.ambient_storage_rate_w_m2_k,
                a.generation_w_m2_k+a.arithmetic_residual_w_m2_k, 1e-15);
            assert!(a.generation_w_m2_k > 0.0);
            assert_eq!(p, saved);
        }
    }
}

#[test]
fn adjacent_temperatures_retain_the_quadratic_entropy_signal() {
    let t0 = 300.0_f64;
    let p = patch(0.7, t0);
    for t in [f64::from_bits(t0.to_bits()-1), f64::from_bits(t0.to_bits()+1)] {
        let a = p.material_entropy_exchange(t).unwrap();
        let delta = t-t0;
        let leading = 4.0*0.7*STEFAN_BOLTZMANN_W_M2_K4*t0*delta*delta;
        assert!(a.generation_w_m2_k > 0.0);
        close(a.generation_w_m2_k/leading, 1.0, 2e-14);
    }
    let a = p.material_entropy_exchange(t0).unwrap();
    assert_eq!(a.heat_out_w_m2, 0.0);
    assert_eq!(a.generation_w_m2_k, 0.0);
    assert_eq!(a.arithmetic_residual_w_m2_k, 0.0);
    let emission = GrayPhotonEmission::new(t0, p.emissivity().value()).unwrap()
        .bolometric(PhotonEntropyQuadrature::default(), || false).unwrap();
    assert!(emission.entropy_w_m2_k > 0.0, "equilibrium net exchange is not zero emission");
}

#[test]
fn source_domains_are_retained_and_bath_temperature_is_not_a_surface_query() {
    let p = patch(0.5, 1000.0);
    let saved = p.clone();
    let a = p.material_entropy_exchange(350.0).unwrap();
    assert!(a.heat_out_w_m2 < 0.0 && a.generation_w_m2_k > 0.0);
    for invalid in [249.0, 501.0, 0.0, f64::NAN, f64::INFINITY] {
        assert!(p.material_entropy_exchange(invalid).is_err());
    }
    assert_eq!(p, saved);
}
