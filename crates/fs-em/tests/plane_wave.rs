#![cfg(feature = "plane-wave")]
//! Synthetic source-bound plane-wave checks, not material qualification.
use fs_em::plane_wave::{DielectricResponse, PlanarDielectricInterface, Polarization};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    field_storage::{
        ELECTRIC_PERMITTIVITY_DIMS as EPS, ELECTRIC_PERMITTIVITY_PROPERTY as EPS_NAME,
        MAGNETIC_PERMEABILITY_DIMS as MU, MAGNETIC_PERMEABILITY_PROPERTY as MU_NAME,
        dielectric::{DebyePoleProperties, ResolvedDebyeCapacitor},
    },
    state_point::{
        MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
        ScalarPropertyRequirement, resolve_material_state_point,
    },
};
use fs_math::c64::C64;
use fs_qty::{Dims, Time};

const E: f64 = 1e-11;
const M: f64 = 1e-6;

fn source(
    properties: &[(&str, Dims, f64)],
    validity: ValidityDomain,
    point: QueryPoint,
) -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for &(name, dims, value) in properties {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: validity.clone(),
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic plane-wave constitutive fixture".into(),
                license: "test-only".into(),
                artifact: None,
            },
        }).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(
            name, dims, ScalarAdmissibility::StrictlyPositive,
        ).unwrap());
    }
    let card = MaterialCard::assemble(MaterialStateId {
        chemistry: "synthetic dielectric".into(),
        phase: "solid".into(),
        process: "analytic scattering fixture".into(),
        revision: 0,
    }, claims, vec![]).unwrap();
    resolve_material_state_point(
        &card, &point, &requirements, MaterialPropertySelection::SingleClaimOnly,
    ).unwrap()
}
fn material(epsilon: f64, mu: f64) -> ResolvedMaterialStatePoint {
    source(
        &[(EPS_NAME, EPS, epsilon), (MU_NAME, MU, mu)],
        ValidityDomain::unconstrained().with("T", 280.0, 320.0),
        QueryPoint::new().with("T", 300.0).unwrap(),
    )
}
fn near(a: f64, b: f64, tol: f64) {
    assert!((a-b).abs() <= tol, "{a} != {b}");
}
fn complex_near(a: C64, b: C64, tol: f64) {
    assert!((a-b).abs() <= tol, "{a:?} != {b:?}");
}

#[test]
fn normal_incidence_and_magnetic_impedance_use_the_original_sources() {
    let incident = material(E, M);
    let transmitted = material(4.0*E, M);
    let interface = PlanarDielectricInterface::try_new(
        &incident, DielectricResponse::Lossless(&transmitted),
    ).unwrap();
    assert!(std::ptr::eq(interface.incident_material(), &incident));
    assert!(std::ptr::eq(interface.transmitted_material(), &transmitted));
    for p in [Polarization::Te, Polarization::Tm] {
        let r = interface.response(1e9, 1.0, p).unwrap();
        complex_near(r.reflected_tangential_e, C64::from_re(-1.0/3.0), 1e-14);
        complex_near(r.transmitted_tangential_e, C64::from_re(2.0/3.0), 1e-14);
        near(r.reflected_power_fraction, 1.0/9.0, 1e-14);
        near(r.transmitted_power_fraction, 8.0/9.0, 1e-14);
    }
    // Small positive index must survive the normal-incidence subtraction.
    let unit_medium = material(1.0, 1.0);
    let low_index = material(1e-20, 1.0);
    let low = PlanarDielectricInterface::try_new(
        &unit_medium, DielectricResponse::Lossless(&low_index),
    ).unwrap();
    for p in [Polarization::Te, Polarization::Tm] {
        let r = low.response(1.0, 1.0, p).unwrap();
        near(r.transmitted_normal_wavenumber_rad_m.re, 1e-10, 1e-24);
        near(r.transmitted_power_fraction, 4e-10, 1e-18);
    }

    // Equal impedance, different phase velocity. Assuming mu=mu0 would fail.
    let magnetic = material(4.0*E, 4.0*M);
    let matched = PlanarDielectricInterface::try_new(
        &incident, DielectricResponse::Lossless(&magnetic),
    ).unwrap();
    for p in [Polarization::Te, Polarization::Tm] {
        let r = matched.response(1e9, 1.0, p).unwrap();
        complex_near(r.reflected_tangential_e, C64::ZERO, 1e-14);
        near(r.transmitted_power_fraction, 1.0, 1e-14);
        near(r.transmitted_normal_wavenumber_rad_m.re
            / r.incident_normal_wavenumber_rad_m, 4.0, 1e-14);
    }
}

#[test]
fn brewster_zero_and_lossless_poynting_flux_are_polarization_specific() {
    let incident = material(E, M);
    let transmitted = material(4.0*E, M);
    let interface = PlanarDielectricInterface::try_new(
        &incident, DielectricResponse::Lossless(&transmitted),
    ).unwrap();
    let cosine = 1.0/5.0_f64.sqrt(); // tan(theta_B)=n_t/n_i=2
    let tm = interface.response(1e9, cosine, Polarization::Tm).unwrap();
    let te = interface.response(1e9, cosine, Polarization::Te).unwrap();
    assert!(tm.reflected_tangential_e.abs() < 2e-15);
    near(tm.transmitted_power_fraction, 1.0, 1e-14);
    near(te.reflected_power_fraction, 0.36, 1e-14);
    for c in [0.01, 0.3, 0.8, 1.0] {
        for p in [Polarization::Te, Polarization::Tm] {
            let r = interface.response(1e9, c, p).unwrap();
            let ki = r.incident_normal_wavenumber_rad_m;
            let kt = r.transmitted_normal_wavenumber_rad_m.re;
            let ratio = match p {
                Polarization::Te => kt/ki,
                Polarization::Tm => 4.0*ki/kt,
            };
            near(r.transmitted_power_fraction,
                r.transmitted_tangential_e.norm_sq()*ratio, 2e-14);
            near(r.reflected_power_fraction+r.transmitted_power_fraction, 1.0, 2e-14);
        }
    }
}

#[test]
fn total_internal_reflection_keeps_phase_and_a_decaying_evanescent_field() {
    let incident = material(4.0*E, M);
    let transmitted = material(E, M);
    let interface = PlanarDielectricInterface::try_new(
        &incident, DielectricResponse::Lossless(&transmitted),
    ).unwrap();
    for p in [Polarization::Te, Polarization::Tm] {
        let r = interface.response(1e9, 0.6, p).unwrap();
        near(r.reflected_power_fraction, 1.0, 2e-14);
        near(r.transmitted_power_fraction, 0.0, 1e-14);
        assert_eq!(r.transmitted_normal_wavenumber_rad_m.re, 0.0);
        assert!(r.transmitted_normal_wavenumber_rad_m.im < 0.0);
        assert!(r.reflected_tangential_e.im.abs() > 0.4,
            "unit magnitude must not erase the TIR phase");
        complex_near(r.transmitted_tangential_e,
            C64::ONE+r.reflected_tangential_e, 1e-14);
        match p {
            Polarization::Te => assert!(r.reflected_tangential_e.im > 0.0),
            Polarization::Tm => assert!(r.reflected_tangential_e.im < 0.0),
        }
    }
    // Exact dyadic critical angle: (n_t/n_i)^2=3/4 and cos(theta)=1/2.
    let critical_medium = material(3.0*E, M);
    let critical = PlanarDielectricInterface::try_new(
        &incident, DielectricResponse::Lossless(&critical_medium),
    ).unwrap();
    for (p, expected) in [(Polarization::Te, 1.0), (Polarization::Tm, -1.0)] {
        let r = critical.response(1e9, 0.5, p).unwrap();
        complex_near(r.reflected_tangential_e, C64::from_re(expected), 1e-14);
        near(r.transmitted_power_fraction, 0.0, 1e-14);
    }
}

#[test]
fn identical_media_remain_transparent_arbitrarily_close_to_grazing() {
    let medium = material(E, M);
    let interface = PlanarDielectricInterface::try_new(
        &medium, DielectricResponse::Lossless(&medium),
    ).unwrap();
    for cosine in [1.0, 0.5, 1e-10, 1e-200] {
        for p in [Polarization::Te, Polarization::Tm] {
            let r = interface.response(1e9, cosine, p).unwrap();
            assert_eq!(r.reflected_tangential_e, C64::ZERO);
            assert_eq!(r.transmitted_tangential_e, C64::ONE);
            near(r.transmitted_power_fraction, 1.0, 1e-14);
        }
    }
}

#[test]
fn original_debye_response_supplies_phase_flux_and_geometry_independent_bulk_epsilon() {
    let incident = material(E, M);
    let state = source(
        &[(EPS_NAME, EPS, 2.0*E), (MU_NAME, MU, M),
          ("delta", EPS, 2.0*E), ("tau", Time::DIMS, 1e-9)],
        ValidityDomain::unconstrained().with("T", 280.0, 320.0),
        QueryPoint::new().with("T", 300.0).unwrap(),
    );
    let poles = [DebyePoleProperties {
        permittivity_increment: "delta".into(), relaxation_time: "tau".into(),
    }];
    let a = ResolvedDebyeCapacitor::try_from_material(
        state.clone(), EPS_NAME, &poles, 0.02, 0.001,
    ).unwrap();
    let b = ResolvedDebyeCapacitor::try_from_material(
        state, EPS_NAME, &poles, 0.3, 0.02,
    ).unwrap();
    let ia = PlanarDielectricInterface::try_new(&incident, DielectricResponse::Debye(&a)).unwrap();
    let ib = PlanarDielectricInterface::try_new(&incident, DielectricResponse::Debye(&b)).unwrap();
    assert!(std::ptr::eq(ia.transmitted_material(), a.material()));
    for p in [Polarization::Te, Polarization::Tm] {
        for c in [0.2, 0.7, 1.0] {
            let r = ia.response(1e9, c, p).unwrap();
            let other = ib.response(1e9, c, p).unwrap();
            complex_near(r.transmitted_permittivity_f_m, C64::new(3.0*E, -E), 1e-24);
            complex_near(r.reflected_tangential_e, other.reflected_tangential_e, 2e-14);
            if c == 1.0 {
                // Independent 80-digit Cartesian evaluation of sqrt(3-i).
                complex_near(r.reflected_tangential_e,
                    C64::new(-0.28180664470320127, 0.07424789305525832), 2e-14);
            }
            assert!(r.transmitted_normal_wavenumber_rad_m.re > 0.0);
            assert!(r.transmitted_normal_wavenumber_rad_m.im < 0.0);
            let ki = r.incident_normal_wavenumber_rad_m;
            let kt = r.transmitted_normal_wavenumber_rad_m;
            let (yi, yt) = match p {
                Polarization::Te => (ki/(1e9*M), kt.scale(1.0/(1e9*M))),
                Polarization::Tm => (1e9*E/ki, r.transmitted_permittivity_f_m.scale(1e9)/kt),
            };
            let flux = r.transmitted_tangential_e.norm_sq()*yt.re/yi;
            near(r.transmitted_power_fraction, flux, 2e-14);
            near(r.reflected_power_fraction+flux, 1.0, 2e-14);
            assert!(flux > 0.0 && flux < 1.0);
        }
    }
}

#[test]
fn missing_magnetic_source_and_invalid_incidence_refuse() {
    let medium = material(E, M);
    let missing_mu = source(
        &[(EPS_NAME, EPS, E)], ValidityDomain::unconstrained(), QueryPoint::new(),
    );
    assert!(PlanarDielectricInterface::try_new(
        &medium, DielectricResponse::Lossless(&missing_mu),
    ).is_err());
    // A frequency-limited static scalar is not promoted into this explicit
    // nondispersive law just because its current query happens to be in-band.
    let bounded = source(
        &[(EPS_NAME, EPS, E), (MU_NAME, MU, M)],
        ValidityDomain::unconstrained().with("omega", 0.0, 1e9),
        QueryPoint::new().with("omega", 5e8).unwrap(),
    );
    assert!(PlanarDielectricInterface::try_new(
        &medium, DielectricResponse::Lossless(&bounded),
    ).is_err());
    // Positive finite source values must not underflow into a zero-index medium.
    let extreme_incident = material(1e300, 1e-300);
    let extreme_transmitted = material(1e-300, 1e-300);
    let extreme = PlanarDielectricInterface::try_new(
        &extreme_incident, DielectricResponse::Lossless(&extreme_transmitted),
    ).unwrap();
    assert!(extreme.response(1e9, 1.0, Polarization::Te).is_err());
    assert!(extreme.response(1e9, 1.0, Polarization::Tm).is_err());
    let interface = PlanarDielectricInterface::try_new(
        &medium, DielectricResponse::Lossless(&medium),
    ).unwrap();
    for (omega, cosine) in [(0.0, 1.0), (f64::NAN, 1.0), (1e9, 0.0),
        (1e9, -0.1), (1e9, 1.1), (1e9, f64::NAN)] {
        assert!(interface.response(omega, cosine, Polarization::Te).is_err());
    }
}
