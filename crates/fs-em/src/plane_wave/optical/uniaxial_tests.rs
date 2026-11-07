//! Actual source-to-interface controls. These cards are authored numerical
//! fixtures, not measured optical data or crystal qualification.
use super::*;
use crate::plane_wave::{
    DielectricResponse, PlanarDielectricInterface, Polarization,
    uniaxial::{PrincipalPlaneUniaxialInterface, UniaxialWaveKind},
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    field_storage::{
        ELECTRIC_PERMITTIVITY_DIMS, ELECTRIC_PERMITTIVITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS,
        MAGNETIC_PERMEABILITY_PROPERTY,
    },
    state_point::{
        MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
        SpectralComplexIndexKnot, resolve_material_state_point,
    },
};
use fs_math::c64::C64;
use fs_qty::Dims;

const MU: f64 = 1.256_637_062_12e-6;
const C0: f64 = 299_792_458.0;
const SPAN: f64 = 410e-9;

fn snapshot(properties: &[(&str, Dims, f64)], revision: u32) -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for &(name, dims, value) in properties {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained(),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "authored two-axis optical Maxwell control".into(),
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
            chemistry: "explicit principal-plane optical control".into(),
            phase: "caller-owned transparent homogeneous state".into(),
            process: "source-role and Maxwell admission control".into(),
            revision,
        },
        claims,
        vec![],
    )
    .unwrap();
    resolve_material_state_point(
        &card,
        &QueryPoint::new(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn axes(
    ordinary: [f64; 2],
    extraordinary: [f64; 2],
    extinction: [f64; 2],
    mu: Option<f64>,
    extraordinary_band: [f64; 2],
    revision: u32,
) -> (
    SpectralComplexIndexStatePoint,
    SpectralComplexIndexStatePoint,
) {
    let mut properties = vec![
        ("ordinary_n_0", Dims::NONE, ordinary[0]),
        ("ordinary_n_1", Dims::NONE, ordinary[1]),
        ("ordinary_k_0", Dims::NONE, extinction[0]),
        ("ordinary_k_1", Dims::NONE, extinction[0]),
        ("extraordinary_n_0", Dims::NONE, extraordinary[0]),
        ("extraordinary_n_1", Dims::NONE, extraordinary[1]),
        ("extraordinary_k_0", Dims::NONE, extinction[1]),
        ("extraordinary_k_1", Dims::NONE, extinction[1]),
    ];
    if let Some(mu) = mu {
        properties.push((
            MAGNETIC_PERMEABILITY_PROPERTY,
            MAGNETIC_PERMEABILITY_DIMS,
            mu,
        ));
    }
    let original = snapshot(&properties, revision);
    let bindings = |prefix: &str, wavelengths: [f64; 2]| {
        [0, 1].map(|i| SpectralComplexIndexKnot {
            wavelength_nm: wavelengths[i],
            real_index_property: PropertyKey::new(format!("{prefix}_n_{i}"), Dims::NONE),
            extinction_index_property: PropertyKey::new(format!("{prefix}_k_{i}"), Dims::NONE),
        })
    };
    let ordinary = SpectralComplexIndexStatePoint::try_from_resolved(
        original.clone(),
        &bindings("ordinary", [400.0, 800.0]),
    )
    .unwrap();
    let extraordinary = SpectralComplexIndexStatePoint::try_from_resolved(
        original,
        &bindings("extraordinary", extraordinary_band),
    )
    .unwrap();
    (ordinary, extraordinary)
}

fn incident(index: f64, mu: f64) -> ResolvedMaterialStatePoint {
    snapshot(
        &[
            (
                ELECTRIC_PERMITTIVITY_PROPERTY,
                ELECTRIC_PERMITTIVITY_DIMS,
                index * index / (mu * C0 * C0),
            ),
            (
                MAGNETIC_PERMEABILITY_PROPERTY,
                MAGNETIC_PERMEABILITY_DIMS,
                mu,
            ),
        ],
        0,
    )
}

fn frequency(wavelength_nm: f64) -> f64 {
    (core::f64::consts::TAU * C0 * 1e9) / wavelength_nm
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}; tolerance={tolerance:.3e}"
    );
}

fn complex_near(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:?} != {expected:?}"
    );
}

#[test]
fn original_axis_roles_complete_snapshot_and_dispersion_survive_into_refraction() {
    let (o, e) = axes(
        [1.4, 1.6],
        [1.8, 2.0],
        [0.0; 2],
        Some(MU),
        [400.0, 800.0],
        0,
    );
    let law = UniaxialOpticalResponse::try_new(&o, &e, SPAN).unwrap();
    assert!(core::ptr::eq(law.ordinary_source(), &o));
    assert!(core::ptr::eq(law.extraordinary_source(), &e));
    assert_eq!(law.material().identity(), o.resolved().identity());
    assert_eq!(law.material().identity(), e.resolved().identity());
    assert_ne!(o.identity(), e.identity());
    assert_eq!(law.permeability_h_m(), MU);
    assert_eq!(law.wavelength_bounds_nm(), [400.0, 800.0]);
    let outside = incident(1.0, MU);
    let interface = PrincipalPlaneUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside),
        law,
        [1.0, 0.0],
    )
    .unwrap();
    for (wavelength, indices) in [
        (400.0, [1.4, 1.8]),
        (600.0, [1.5, 1.9]),
        (800.0, [1.6, 2.0]),
    ] {
        let omega = frequency(wavelength);
        let actual = law.phase_indices(omega).unwrap();
        for i in [0, 1] {
            near(actual[i], indices[i], 4e-15);
        }
        for (polarization, n) in [
            (Polarization::Te, indices[0]),
            (Polarization::Tm, indices[1]),
        ] {
            let result = interface.response(omega, 1.0, polarization).unwrap();
            assert_eq!(result.wave_kind, UniaxialWaveKind::Propagating);
            near(
                result.reflected_tangential_e.re,
                (1.0 - n) / (1.0 + n),
                4e-15,
            );
            near(
                result.transmitted_power_fraction,
                4.0 * n / ((1.0 + n) * (1.0 + n)),
                4e-15,
            );
        }
    }
}

#[test]
fn same_state_is_required_even_when_axis_values_and_permeability_match() {
    let (o, e) = axes([1.5; 2], [1.7; 2], [0.0; 2], Some(MU), [400.0, 800.0], 0);
    let (_, other) = axes([1.5; 2], [1.7; 2], [0.0; 2], Some(MU), [400.0, 800.0], 1);
    assert!(UniaxialOpticalResponse::try_new(&o, &e, SPAN).is_ok());
    assert!(UniaxialOpticalResponse::try_new(&o, &other, SPAN).is_err());
}

#[test]
fn missing_permeability_interpolation_gaps_and_disjoint_support_refuse() {
    let (o, e) = axes([1.5; 2], [1.7; 2], [0.0; 2], None, [400.0, 800.0], 0);
    assert!(UniaxialOpticalResponse::try_new(&o, &e, SPAN).is_err());
    let (o, e) = axes([1.5; 2], [1.7; 2], [0.0; 2], Some(MU), [900.0, 1000.0], 0);
    assert!(UniaxialOpticalResponse::try_new(&o, &e, SPAN).is_err());
    let (o, e) = axes([1.5; 2], [1.7; 2], [0.0; 2], Some(MU), [400.0, 800.0], 0);
    for span in [0.0, f64::NAN, 399e-9] {
        assert!(UniaxialOpticalResponse::try_new(&o, &e, span).is_err());
    }
}

#[test]
fn each_original_spectral_support_is_enforced_before_endpoint_rounding() {
    let (o, e) = axes([1.5; 2], [1.7; 2], [0.0; 2], Some(MU), [500.0, 700.0], 0);
    let law = UniaxialOpticalResponse::try_new(&o, &e, SPAN).unwrap();
    assert_eq!(law.wavelength_bounds_nm(), [500.0, 700.0]);
    assert!(law.phase_indices(frequency(500.0)).is_ok());
    assert!(law.phase_indices(frequency(700.0)).is_ok());
    assert!(law.phase_indices(frequency(500.0).next_up()).is_err());
    assert!(law.phase_indices(frequency(700.0).next_down()).is_err());
    assert!(law.phase_indices(frequency(450.0)).is_err());
    assert!(law.phase_indices(frequency(750.0)).is_err());
}

#[test]
fn extinction_in_either_original_axis_is_not_silently_dropped() {
    for loss in [[0.0, 1e-12], [1e-12, 0.0]] {
        let (o, e) = axes([1.5; 2], [1.7; 2], loss, Some(MU), [400.0, 800.0], 0);
        let law = UniaxialOpticalResponse::try_new(&o, &e, SPAN).unwrap();
        assert!(law.phase_indices(frequency(550.0)).is_err());
        let outside = incident(1.0, MU);
        let interface = PrincipalPlaneUniaxialInterface::try_new(
            DielectricResponse::Lossless(&outside),
            law,
            [0.6, 0.8],
        )
        .unwrap();
        for mode in [Polarization::Te, Polarization::Tm] {
            assert!(interface.response(frequency(550.0), 0.8, mode).is_err());
        }
    }
}

#[test]
fn isotropic_limit_matches_the_original_source_bound_interface_owner() {
    let (o, _) = axes(
        [1.3, 1.7],
        [1.3, 1.7],
        [0.0; 2],
        Some(1.6 * MU),
        [400.0, 800.0],
        0,
    );
    let law = UniaxialOpticalResponse::try_new(&o, &o, SPAN).unwrap();
    let outside = incident(1.5, MU);
    let native = PlanarDielectricInterface::try_new(
        &outside,
        DielectricResponse::SpectralOptical(
            SpectralOpticalResponse::try_new(SpectralOpticalSource::ComplexIndex(&o), SPAN)
                .unwrap(),
        ),
    )
    .unwrap();
    let uniaxial = PrincipalPlaneUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside),
        law,
        [0.37, -0.81],
    )
    .unwrap();
    for wavelength in [400.0, 550.0, 800.0] {
        for cosine in [0.1, 0.5, 0.9, 1.0] {
            for mode in [Polarization::Te, Polarization::Tm] {
                let reference = native
                    .response(frequency(wavelength), cosine, mode)
                    .unwrap();
                let actual = uniaxial
                    .response(frequency(wavelength), cosine, mode)
                    .unwrap();
                complex_near(
                    actual.reflected_tangential_e,
                    reference.reflected_tangential_e,
                    3e-13,
                );
                complex_near(
                    actual.transmitted_tangential_e,
                    reference.transmitted_tangential_e,
                    3e-13,
                );
                near(
                    actual.transmitted_power_fraction,
                    reference.transmitted_power_fraction,
                    3e-13,
                );
            }
        }
    }
}

#[test]
fn sourced_tilted_crystal_keeps_negative_phase_normal_and_positive_energy() {
    let (o, e) = axes([1.0; 2], [2.0; 2], [0.0; 2], Some(MU), [400.0, 800.0], 0);
    let law = UniaxialOpticalResponse::try_new(&o, &e, SPAN).unwrap();
    let outside = incident(2.0, MU);
    let interface = PrincipalPlaneUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside),
        law,
        [1.0, 1.0],
    )
    .unwrap();
    let response = interface
        .response_at_tangential_index(frequency(550.0), 1.5, Polarization::Tm)
        .unwrap();
    assert!(response.transmitted_normal_wavenumber_rad_m.re < 0.0);
    assert!(response.phase_direction.unwrap()[1] < 0.0);
    assert!(response.energy_direction.unwrap()[1] > 0.0);
    assert!(response.transmitted_power_fraction > 0.0);
    near(
        response.reflected_power_fraction + response.transmitted_power_fraction,
        1.0,
        2e-14,
    );
}

#[test]
fn invalid_axis_and_absorbing_incident_port_are_not_admitted_as_lossless() {
    let (o, e) = axes([1.4; 2], [1.8; 2], [0.0; 2], Some(MU), [400.0, 800.0], 0);
    let law = UniaxialOpticalResponse::try_new(&o, &e, SPAN).unwrap();
    let outside = incident(1.0, MU);
    for axis in [[0.0, 0.0], [f64::NAN, 1.0], [1.0, f64::INFINITY]] {
        assert!(
            PrincipalPlaneUniaxialInterface::try_new(
                DielectricResponse::Lossless(&outside),
                law,
                axis,
            )
            .is_err()
        );
    }
    let (lossy, _) = axes([1.0; 2], [1.0; 2], [0.01, 0.0], Some(MU), [400.0, 800.0], 0);
    let incoming =
        SpectralOpticalResponse::try_new(SpectralOpticalSource::ComplexIndex(&lossy), SPAN)
            .unwrap();
    let interface = PrincipalPlaneUniaxialInterface::try_new(
        DielectricResponse::SpectralOptical(incoming),
        law,
        [0.0, 1.0],
    )
    .unwrap();
    assert!(
        interface
            .response(frequency(550.0), 1.0, Polarization::Te)
            .is_err()
    );
}

#[test]
fn oriented_crystal_retains_original_sources_dispersion_and_refusal_domains() {
    use crate::plane_wave::uniaxial::oriented::OrientedUniaxialInterface;
    let (ordinary, extraordinary) = axes(
        [1.4, 1.6],
        [1.8, 2.0],
        [0.0; 2],
        Some(MU),
        [400.0, 800.0],
        0,
    );
    let crystal = UniaxialOpticalResponse::try_new(&ordinary, &extraordinary, SPAN).unwrap();
    let outside = incident(1.0, MU);
    let interface = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside),
        crystal,
        [1.0, 2.0, 3.0],
    )
    .unwrap();
    assert!(core::ptr::eq(
        interface.incident_response().material(),
        &outside
    ));
    assert!(core::ptr::eq(
        interface.crystal_response().ordinary_source(),
        &ordinary
    ));
    assert!(core::ptr::eq(
        interface.crystal_response().extraordinary_source(),
        &extraordinary
    ));
    let illumination = [C64::new(0.3, -0.5), C64::new(0.7, 0.2)];
    for wavelength in [400.0, 550.0, 800.0] {
        let omega = frequency(wavelength);
        let [no, _] = crystal.phase_indices(omega).unwrap();
        let r = interface.response(omega, 0.8, illumination).unwrap();
        let qo = r.transmitted[0].normal_wavenumber_rad_m.re * C0 / omega;
        let p = r.tangential_wavenumber_rad_m * C0 / omega;
        near(qo * qo + p * p, no * no, 4e-15);
        near(
            r.reflected_power_fraction + r.transmitted_power_fraction,
            1.0,
            2e-14,
        );
        for j in 0..2 {
            complex_near(r.incident_electric_v_m[j], illumination[j], 1e-16);
        }
        let signed = interface
            .response_at_tangential_index(omega, -0.6, illumination)
            .unwrap();
        assert!(signed.tangential_wavenumber_rad_m < 0.0);
        near(
            signed.reflected_power_fraction + signed.transmitted_power_fraction,
            1.0,
            2e-14,
        );
    }
    for wavelength in [399.0, 801.0] {
        assert!(
            interface
                .response(frequency(wavelength), 0.8, illumination)
                .is_err()
        );
    }
    for loss in [[1e-12, 0.0], [0.0, 1e-12]] {
        let (o, e) = axes([1.5; 2], [1.9; 2], loss, Some(MU), [400.0, 800.0], 0);
        let lossy = UniaxialOpticalResponse::try_new(&o, &e, SPAN).unwrap();
        let interface = OrientedUniaxialInterface::try_new(
            DielectricResponse::Lossless(&outside),
            lossy,
            [1.0, 2.0, 3.0],
        )
        .unwrap();
        assert!(
            interface
                .response(frequency(550.0), 0.8, illumination)
                .is_err()
        );
    }
}

#[test]
fn prepared_crystal_spectrum_binds_original_laws_and_coherent_spatial_fields() {
    use crate::plane_wave::{
        angular_spectrum::{AngularSpectrumError, AngularSpectrumSide},
        uniaxial::oriented::{
            OrientedUniaxialInterface,
            spectrum::{
                PreparedUniaxialAngularSpectrum, UniaxialIncidentMode, UniaxialSpectrumLimits,
            },
        },
    };
    let (ordinary, extraordinary) = axes(
        [1.4, 1.6],
        [1.8, 2.0],
        [0.0; 2],
        Some(MU),
        [400.0, 800.0],
        0,
    );
    let crystal = UniaxialOpticalResponse::try_new(&ordinary, &extraordinary, SPAN).unwrap();
    let outside = incident(1.0, MU);
    let interface = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside),
        crystal,
        [1.0, 2.0, 3.0],
    )
    .unwrap();
    let input = [
        UniaxialIncidentMode {
            tangential_phase_index: 0.6,
            incidence_plane_axis: [1.0, 0.0],
            electric_v_m: [C64::ONE, C64::ZERO],
        },
        UniaxialIncidentMode {
            tangential_phase_index: -0.4,
            incidence_plane_axis: [3.0, -4.0],
            electric_v_m: [C64::new(0.2, 0.7), C64::new(-0.3, 0.1)],
        },
    ];
    let limits = UniaxialSpectrumLimits {
        maximum_modes: 2,
        maximum_workspace_bytes: 10000,
    };
    let mut depths = Vec::new();
    for wavelength in [400.0, 550.0, 800.0] {
        let omega = frequency(wavelength);
        let prepared =
            PreparedUniaxialAngularSpectrum::prepare(&interface, omega, &input, limits, || false)
                .unwrap();
        assert!(core::ptr::eq(prepared.source(), &interface));
        assert_eq!(prepared.angular_frequency_rad_s(), omega);
        assert_eq!(prepared.retained_mode_count(), 2);
        assert!(prepared.workspace_bytes() <= limits.maximum_workspace_bytes);
        let point = [0.3e-6, -0.5e-6, 0.0];
        let left = prepared
            .field_at(AngularSpectrumSide::TotalIncidentSide, point, || false)
            .unwrap();
        let right = prepared
            .field_at(AngularSpectrumSide::Transmitted, point, || false)
            .unwrap();
        for j in 0..2 {
            complex_near(left.electric_v_m[j], right.electric_v_m[j], 3e-14);
            complex_near(left.magnetic_a_m[j], right.magnetic_a_m[j], 3e-16);
        }
        near(left.poynting_w_m2[2], right.poynting_w_m2[2], 3e-16);
        depths.push(
            prepared
                .field_at(AngularSpectrumSide::Transmitted, [0.0, 0.0, 1e-6], || false)
                .unwrap(),
        );
    }
    assert!((depths[0].electric_v_m[0] - depths[2].electric_v_m[0]).abs() > 0.1);
    assert!(matches!(
        PreparedUniaxialAngularSpectrum::prepare(
            &interface,
            frequency(399.0),
            &input,
            limits,
            || false,
        ),
        Err(AngularSpectrumError::Source(_))
    ));
    assert!(matches!(
        PreparedUniaxialAngularSpectrum::prepare(
            &interface,
            frequency(399.0),
            &input,
            UniaxialSpectrumLimits {
                maximum_modes: 0,
                ..limits
            },
            || false,
        ),
        Err(AngularSpectrumError::Invalid(_))
    ));
    assert!(matches!(
        PreparedUniaxialAngularSpectrum::prepare(
            &interface,
            frequency(399.0),
            &input,
            limits,
            || true,
        ),
        Err(AngularSpectrumError::Cancelled)
    ));
}
