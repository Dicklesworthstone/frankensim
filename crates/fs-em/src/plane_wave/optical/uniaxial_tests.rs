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
    let mut card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "explicit principal-plane optical control".into(),
            phase: "caller-owned transparent homogeneous state".into(),
            process: "source-role and Maxwell admission control".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    // Revision is part of the original source identity. Preserve the
    // requested revision through MatDb's immutable predecessor chain.
    for _ in 0..revision {
        card = MaterialCard::supersede(&card, card.claims().clone(), vec![]).unwrap();
    }
    assert_eq!(card.id().revision, revision);
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
fn original_absorbing_axes_feed_coupled_fields_and_physical_heating() {
    use crate::plane_wave::uniaxial::oriented::{
        OrientedCrystalFieldSide, OrientedUniaxialInterface,
    };
    let source_incident = incident(1.0, MU);
    let (o, e) = axes(
        [1.4, 1.6],
        [1.8, 2.0],
        [0.08, 0.2],
        Some(MU),
        [400.0, 800.0],
        4,
    );
    let law = UniaxialOpticalResponse::try_new(&o, &e, SPAN).unwrap();
    let interface = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&source_incident),
        law,
        [1.0, 2.0, -0.5],
    )
    .unwrap();
    let electric = [C64::new(0.6, -0.2), C64::new(0.1, 0.7)];
    for wavelength in [400.0, 550.0, 800.0] {
        let omega = frequency(wavelength);
        assert!(interface.response(omega, 0.8, electric).is_err());
        let prepared = interface
            .prepare_absorbing_field(omega, 0.8, electric)
            .unwrap();
        assert!(core::ptr::eq(prepared.source(), &interface));
        assert!(core::ptr::eq(
            prepared.source().crystal_response().ordinary_source(),
            &o
        ));
        assert!(core::ptr::eq(
            prepared.source().crystal_response().extraordinary_source(),
            &e
        ));
        assert_eq!(prepared.omega_rad_s(), omega);
        let expected_loss = law
            .principal_permittivities_f_m(omega)
            .unwrap()
            .map(|v| -v.im);
        assert_eq!(
            prepared
                .boundary_response()
                .principal_loss_permittivities_f_m(),
            expected_loss
        );
        let boundary = prepared
            .field(OrientedCrystalFieldSide::Transmitted, [0.0; 3])
            .unwrap();
        let depth = prepared
            .field(
                OrientedCrystalFieldSide::Transmitted,
                [0.0, 0.0, wavelength * 1e-9],
            )
            .unwrap();
        assert!(boundary.absorbed_power_density_w_m3 > 0.0);
        assert!(depth.absorbed_power_density_w_m3 > 0.0);
        assert!(depth.field.poynting_w_m2[2] < boundary.field.poynting_w_m2[2]);
        let area = 2e-8;
        let end = wavelength * 1e-9;
        let absorbed = prepared.slab_absorption(area, 0.0, end).unwrap();
        assert!(core::ptr::eq(absorbed.prepared_field(), &prepared));
        assert!(core::ptr::eq(absorbed.source(), &interface));
        assert_eq!(
            absorbed.source().crystal_response().material().identity(),
            law.material().identity()
        );
        assert_eq!(absorbed.area_m2(), area);
        assert_eq!(absorbed.start_depth_m(), 0.0);
        assert_eq!(absorbed.end_depth_m(), end);
        assert_eq!(absorbed.volume_m3(), area * end);
        near(
            absorbed.absorbed_power_w(),
            area * (boundary.field.poynting_w_m2[2] - depth.field.poynting_w_m2[2]),
            3e-12 * absorbed.absorbed_power_w(),
        );
        for (area, start, end) in [
            (0.0, 0.0, end),
            (-1.0, 0.0, end),
            (f64::INFINITY, 0.0, end),
            (area, -end, end),
            (area, end, end),
            (area, end, 0.0),
            (area, 0.0, f64::INFINITY),
            (area, f64::NAN, end),
            (f64::MAX, 0.0, f64::MAX),
        ] {
            assert!(prepared.slab_absorption(area, start, end).is_err());
        }
        let by_index = interface
            .prepare_absorbing_field_at_tangential_index(omega, 0.6, electric)
            .unwrap();
        near(
            prepared.boundary_response().reflected_power_fraction,
            by_index.boundary_response().reflected_power_fraction,
            2e-14,
        );
    }
    for omega in [frequency(800.0).next_down(), frequency(400.0).next_up()] {
        assert!(
            interface
                .prepare_absorbing_field(omega, 0.8, electric)
                .is_err()
        );
    }
    // The original spectral owner rejects the unsupported knot gap during
    // source admission, before an absorbing interface can be prepared.
    assert!(UniaxialOpticalResponse::try_new(&o, &e, 10e-9).is_err());
}

#[test]
fn crystal_statistical_illumination_retains_original_spectra_volume_and_source_identity() {
    use crate::plane_wave::uniaxial::oriented::{
        OrientedCrystalFieldSide, OrientedUniaxialInterface,
        absorbing::{AbsorbingCrystalIllumination, CrystalIlluminationError},
    };
    let outside = incident(1.0, MU);
    let (o, e) = axes(
        [1.4, 1.6],
        [1.8, 2.0],
        [0.08, 0.2],
        Some(MU),
        [400.0, 800.0],
        4,
    );
    let law = UniaxialOpticalResponse::try_new(&o, &e, SPAN).unwrap();
    let interface = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside),
        law,
        [1.0, 2.0, -0.5],
    )
    .unwrap();
    let factors = [
        interface
            .prepare_absorbing_field_in_plane(
                frequency(400.0),
                0.8,
                [1.0, 0.0],
                [C64::new(0.6, -0.2), C64::new(0.1, 0.7)],
            )
            .unwrap(),
        interface
            .prepare_absorbing_field_at_tangential_index_in_plane(
                frequency(800.0),
                -0.7,
                [3.0, 4.0],
                [C64::new(-0.2, 0.4), C64::new(0.3, 0.1)],
            )
            .unwrap(),
    ];
    let illumination =
        AbsorbingCrystalIllumination::try_new(&interface, &factors, 2, || false).unwrap();
    assert!(core::ptr::eq(illumination.source(), &interface));
    assert!(core::ptr::eq(
        illumination.factors().as_ptr(),
        factors.as_ptr()
    ));
    assert_eq!(illumination.factor_count(), 2);
    assert_eq!(illumination.factors()[0].omega_rad_s(), frequency(400.0));
    assert_eq!(illumination.factors()[1].omega_rad_s(), frequency(800.0));
    assert_eq!(
        illumination.factors()[1]
            .boundary_response()
            .incidence_plane_axis(),
        [0.6, 0.8]
    );
    assert_eq!(
        factors[0].boundary_response(),
        interface
            .prepare_absorbing_field(
                frequency(400.0),
                0.8,
                [C64::new(0.6, -0.2), C64::new(0.1, 0.7)]
            )
            .unwrap()
            .boundary_response()
    );
    let side = OrientedCrystalFieldSide::Transmitted;
    let point = [0.0, 0.0, 300e-9];
    let statistics = illumination
        .field_statistics(side, point, || false)
        .unwrap();
    let expected_density = factors
        .iter()
        .map(|f| f.field(side, point).unwrap().absorbed_power_density_w_m3)
        .sum::<f64>();
    near(
        statistics.absorbed_power_density_w_m3,
        expected_density,
        3e-15 * expected_density,
    );
    assert_eq!(
        illumination.poynting_w_m2(side, point, || false).unwrap(),
        statistics.poynting_w_m2
    );
    let absorbed = illumination
        .slab_absorption(2e-8, 100e-9, 900e-9, || false)
        .unwrap();
    assert!(core::ptr::eq(absorbed.illumination(), &illumination));
    assert!(core::ptr::eq(absorbed.source(), &interface));
    assert_eq!(
        absorbed.source().crystal_response().material().identity(),
        law.material().identity()
    );
    assert_eq!(absorbed.area_m2(), 2e-8);
    assert_eq!(absorbed.start_depth_m(), 100e-9);
    assert_eq!(absorbed.end_depth_m(), 900e-9);
    assert_eq!(absorbed.volume_m3(), 2e-8 * (900e-9 - 100e-9));
    let expected_power = factors
        .iter()
        .map(|f| {
            f.slab_absorption(2e-8, 100e-9, 900e-9)
                .unwrap()
                .absorbed_power_w()
        })
        .sum::<f64>();
    near(
        absorbed.absorbed_power_w(),
        expected_power,
        3e-15 * expected_power,
    );
    let power = illumination.power();
    near(
        power.incident_w_m2,
        power.reflected_w_m2 + power.transmitted_w_m2,
        3e-14 * power.incident_w_m2,
    );
    assert!(matches!(
        AbsorbingCrystalIllumination::try_new(&interface, &factors, 1, || false),
        Err(CrystalIlluminationError::Electromagnetic(
            ReflectionError::Invalid(_)
        ))
    ));
    let equivalent_interface = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside),
        law,
        [1.0, 2.0, -0.5],
    )
    .unwrap();
    assert!(matches!(
        AbsorbingCrystalIllumination::try_new(&equivalent_interface, &factors, 2, || false),
        Err(CrystalIlluminationError::Electromagnetic(
            ReflectionError::Invalid(_)
        ))
    ));
    let dark = AbsorbingCrystalIllumination::try_new(&interface, &[], 0, || false).unwrap();
    assert_eq!(dark.power().incident_w_m2, 0.0);
    assert_eq!(
        dark.slab_absorption(2e-8, 0.0, 900e-9, || false)
            .unwrap()
            .absorbed_power_w(),
        0.0
    );
    assert!(dark.slab_absorption(0.0, 0.0, 900e-9, || false).is_err());
    assert!(
        dark.field_statistics(side, [0.0, 0.0, -1.0], || false)
            .is_err()
    );
    assert!(matches!(
        AbsorbingCrystalIllumination::try_new(&interface, &factors, 2, || true),
        Err(CrystalIlluminationError::Cancelled)
    ));
    assert!(matches!(
        illumination.field_statistics(side, point, || true),
        Err(CrystalIlluminationError::Cancelled)
    ));
    assert!(matches!(
        illumination.poynting_w_m2(side, point, || true),
        Err(CrystalIlluminationError::Cancelled)
    ));
    assert!(matches!(
        illumination.slab_absorption(2e-8, 0.0, 900e-9, || true),
        Err(CrystalIlluminationError::Cancelled)
    ));
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
fn prepared_oriented_fields_keep_original_sources_frequency_and_boundary_traces() {
    use crate::plane_wave::uniaxial::oriented::{
        OrientedCrystalFieldSide, OrientedUniaxialInterface,
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
    let omega = frequency(550.0);
    let illumination = [C64::new(0.3, -0.5), C64::new(0.7, 0.2)];
    let prepared = interface.prepare_field(omega, 0.8, illumination).unwrap();
    assert!(core::ptr::eq(prepared.source(), &interface));
    assert!(core::ptr::eq(
        prepared.source().crystal_response().ordinary_source(),
        &ordinary,
    ));
    assert_eq!(prepared.omega_rad_s(), omega);
    assert_eq!(
        *prepared.boundary_response(),
        interface.response(omega, 0.8, illumination).unwrap(),
    );
    let outside_trace = prepared
        .field(OrientedCrystalFieldSide::TotalIncidentSide, [0.0; 3])
        .unwrap();
    let inside_trace = prepared
        .field(OrientedCrystalFieldSide::Transmitted, [0.0; 3])
        .unwrap();
    for j in 0..2 {
        complex_near(
            outside_trace.electric_v_m[j],
            inside_trace.electric_v_m[j],
            3e-15,
        );
        complex_near(
            outside_trace.magnetic_a_m[j],
            inside_trace.magnetic_a_m[j],
            3e-17,
        );
    }
    let point = [0.12e-6, 0.7, 0.37e-6];
    let propagated = prepared
        .field(OrientedCrystalFieldSide::Transmitted, point)
        .unwrap();
    near(
        propagated.poynting_w_m2[2],
        prepared.boundary_response().incident_normal_power_w_m2
            * prepared.boundary_response().transmitted_power_fraction,
        2e-17,
    );
    assert!(
        prepared
            .field(OrientedCrystalFieldSide::Incident, point)
            .is_err()
    );
    assert_eq!(
        propagated,
        prepared
            .field(OrientedCrystalFieldSide::Transmitted, point)
            .unwrap()
    );
    let signed = interface
        .prepare_field_at_tangential_index(omega, -0.6, illumination)
        .unwrap();
    assert!(signed.boundary_response().tangential_wavenumber_rad_m < 0.0);
    assert!(
        signed
            .field(OrientedCrystalFieldSide::Transmitted, point)
            .is_ok()
    );
    for wavelength in [399.0, 801.0] {
        assert!(
            interface
                .prepare_field(frequency(wavelength), 0.8, illumination)
                .is_err()
        );
        assert!(
            interface
                .prepare_field_at_tangential_index(frequency(wavelength), -0.6, illumination)
                .is_err()
        );
    }
    for loss in [[1e-12, 0.0], [0.0, 1e-12]] {
        let (o, e) = axes([1.5; 2], [1.9; 2], loss, Some(MU), [400.0, 800.0], 0);
        let interface = OrientedUniaxialInterface::try_new(
            DielectricResponse::Lossless(&outside),
            UniaxialOpticalResponse::try_new(&o, &e, SPAN).unwrap(),
            [1.0, 2.0, 3.0],
        )
        .unwrap();
        assert!(interface.prepare_field(omega, 0.8, illumination).is_err());
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

#[test]
fn finite_crystal_slab_and_spectrum_retain_original_ports_and_both_interfaces() {
    use crate::plane_wave::{
        angular_spectrum::AngularSpectrumError,
        uniaxial::oriented::{
            OrientedUniaxialInterface,
            slab::{OrientedUniaxialSlab, PreparedUniaxialSlabSpectrum, SlabRegion},
            spectrum::{UniaxialIncidentMode, UniaxialSpectrumLimits},
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
    let exit = incident(1.2, 0.8 * MU);
    let entrance = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside),
        crystal,
        [1.0, 2.0, 3.0],
    )
    .unwrap();
    let slab =
        OrientedUniaxialSlab::try_new(&entrance, DielectricResponse::Lossless(&exit), 0.5e-6)
            .unwrap();
    assert!(core::ptr::eq(slab.entrance_interface(), &entrance));
    assert!(core::ptr::eq(slab.exit_response().material(), &exit));
    let input = [
        UniaxialIncidentMode {
            tangential_phase_index: 0.4,
            incidence_plane_axis: [1.0, 0.0],
            electric_v_m: [C64::new(0.7, 0.2), C64::new(-0.3, 0.5)],
        },
        UniaxialIncidentMode {
            tangential_phase_index: -0.6,
            incidence_plane_axis: [3.0, 4.0],
            electric_v_m: [C64::ONE, C64::ZERO],
        },
    ];
    let limits = UniaxialSpectrumLimits {
        maximum_modes: 2,
        maximum_workspace_bytes: 20000,
    };
    let mut transmitted = Vec::new();
    for wavelength in [400.0, 550.0, 800.0] {
        let omega = frequency(wavelength);
        let one = slab
            .response_at_tangential_index(omega, 0.4, input[0].electric_v_m)
            .unwrap();
        near(
            one.reflected_power_fraction + one.transmitted_power_fraction,
            1.0,
            3e-13,
        );
        transmitted.push(one.transmitted_electric_v_m[0]);
        let prepared =
            PreparedUniaxialSlabSpectrum::prepare(&slab, omega, &input, limits, || false).unwrap();
        assert!(core::ptr::eq(prepared.source(), &slab));
        assert_eq!(prepared.angular_frequency_rad_s(), omega);
        assert_eq!(prepared.retained_mode_count(), 2);
        assert!(prepared.workspace_bytes() <= limits.maximum_workspace_bytes);
        for (z, left, right) in [
            (0.0, SlabRegion::TotalIncidentSide, SlabRegion::Crystal),
            (
                slab.thickness_m(),
                SlabRegion::Crystal,
                SlabRegion::Transmitted,
            ),
        ] {
            let a = prepared
                .field_at(left, [0.3e-6, -0.2e-6, z], || false)
                .unwrap();
            let b = prepared
                .field_at(right, [0.3e-6, -0.2e-6, z], || false)
                .unwrap();
            for j in 0..2 {
                complex_near(a.electric_v_m[j], b.electric_v_m[j], 3e-13);
                complex_near(a.magnetic_a_m[j], b.magnetic_a_m[j], 3e-15);
            }
            near(a.poynting_w_m2[2], b.poynting_w_m2[2], 3e-15);
        }
    }
    assert!((transmitted[0] - transmitted[2]).abs() > 0.1);
    assert!(
        slab.response_at_tangential_index(frequency(399.0), 0.4, input[0].electric_v_m)
            .is_err()
    );
    assert!(matches!(
        PreparedUniaxialSlabSpectrum::prepare(
            &slab,
            frequency(399.0),
            &input,
            UniaxialSpectrumLimits {
                maximum_modes: 0,
                ..limits
            },
            || false
        ),
        Err(AngularSpectrumError::Invalid(_))
    ));
    assert!(matches!(
        PreparedUniaxialSlabSpectrum::prepare(&slab, frequency(399.0), &input, limits, || true),
        Err(AngularSpectrumError::Cancelled)
    ));
    let (lossy, _) = axes([1.2; 2], [1.2; 2], [0.01, 0.0], Some(MU), [400.0, 800.0], 0);
    let absorbing =
        SpectralOpticalResponse::try_new(SpectralOpticalSource::ComplexIndex(&lossy), SPAN)
            .unwrap();
    let refused = OrientedUniaxialSlab::try_new(
        &entrance,
        DielectricResponse::SpectralOptical(absorbing),
        0.5e-6,
    )
    .unwrap();
    assert!(
        refused
            .response(frequency(550.0), 0.8, input[0].electric_v_m)
            .is_err()
    );
}

#[test]
fn coherent_crystal_stack_retains_original_layers_and_matches_the_single_slab_owner() {
    use crate::plane_wave::{
        angular_spectrum::AngularSpectrumError,
        uniaxial::oriented::{
            OrientedUniaxialInterface,
            slab::OrientedUniaxialSlab,
            spectrum::{UniaxialIncidentMode, UniaxialSpectrumLimits},
            stack::{
                OrientedUniaxialStack, PreparedUniaxialStackSpectrum, StackRegion,
                UniaxialStackLayer, UniaxialStackLimits,
            },
        },
    };
    let (ordinary_a, extraordinary_a) = axes(
        [1.4, 1.6],
        [1.8, 2.0],
        [0.0; 2],
        Some(MU),
        [400.0, 800.0],
        0,
    );
    let (ordinary_b, extraordinary_b) = axes(
        [1.3, 1.45],
        [1.65, 1.85],
        [0.0; 2],
        Some(1.2 * MU),
        [400.0, 800.0],
        1,
    );
    let a = UniaxialOpticalResponse::try_new(&ordinary_a, &extraordinary_a, SPAN).unwrap();
    let b = UniaxialOpticalResponse::try_new(&ordinary_b, &extraordinary_b, SPAN).unwrap();
    let outside = incident(1.0, MU);
    let exit = incident(1.2, 0.8 * MU);
    let layers = [
        UniaxialStackLayer::try_new(a, [1.0, 2.0, 3.0], 0.4e-6).unwrap(),
        UniaxialStackLayer::try_new(b, [2.0, -1.0, 1.0], 0.3e-6).unwrap(),
    ];
    let limits = UniaxialStackLimits {
        maximum_layers: 2,
        maximum_workspace_bytes: 50_000,
    };
    let stack = OrientedUniaxialStack::try_new(
        DielectricResponse::Lossless(&outside),
        DielectricResponse::Lossless(&exit),
        &layers,
        limits,
    )
    .unwrap();
    assert!(core::ptr::eq(stack.layers(), layers.as_slice()));
    assert!(core::ptr::eq(
        stack.incident_response().material(),
        &outside
    ));
    assert!(core::ptr::eq(stack.exit_response().material(), &exit));
    let one_layer = [layers[0]];
    let one_stack = OrientedUniaxialStack::try_new(
        DielectricResponse::Lossless(&outside),
        DielectricResponse::Lossless(&exit),
        &one_layer,
        limits,
    )
    .unwrap();
    let entrance = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside),
        a,
        layers[0].optic_axis(),
    )
    .unwrap();
    let slab = OrientedUniaxialSlab::try_new(
        &entrance,
        DielectricResponse::Lossless(&exit),
        layers[0].thickness_m(),
    )
    .unwrap();
    let input = [
        UniaxialIncidentMode {
            tangential_phase_index: 0.4,
            incidence_plane_axis: [1.0, 0.0],
            electric_v_m: [C64::new(0.7, 0.2), C64::new(-0.3, 0.5)],
        },
        UniaxialIncidentMode {
            tangential_phase_index: -0.5,
            incidence_plane_axis: [3.0, 4.0],
            electric_v_m: [C64::ONE, C64::ZERO],
        },
    ];
    let spectrum_limits = UniaxialSpectrumLimits {
        maximum_modes: 2,
        maximum_workspace_bytes: 100_000,
    };
    for wavelength in [400.0, 550.0, 800.0] {
        let omega = frequency(wavelength);
        let x = one_stack
            .response_at_tangential_index(omega, 0.4, input[0].electric_v_m, || false)
            .unwrap();
        let y = slab
            .response_at_tangential_index(omega, 0.4, input[0].electric_v_m)
            .unwrap();
        for j in 0..3 {
            complex_near(
                x.reflected_electric_v_m[j],
                y.reflected_electric_v_m[j],
                4e-13,
            );
            complex_near(
                x.transmitted_electric_v_m[j],
                y.transmitted_electric_v_m[j],
                4e-13,
            );
        }
        let prepared =
            PreparedUniaxialStackSpectrum::prepare(&stack, omega, &input, spectrum_limits, || {
                false
            })
            .unwrap();
        assert!(core::ptr::eq(prepared.source(), &stack));
        assert_eq!(prepared.retained_mode_count(), 2);
        assert!(prepared.workspace_bytes() <= spectrum_limits.maximum_workspace_bytes);
        for (z, left, right) in [
            (0.0, StackRegion::TotalIncidentSide, StackRegion::Layer(0)),
            (
                layers[0].thickness_m(),
                StackRegion::Layer(0),
                StackRegion::Layer(1),
            ),
            (
                stack.thickness_m(),
                StackRegion::Layer(1),
                StackRegion::Transmitted,
            ),
        ] {
            let x = prepared
                .field_at(left, [0.2e-6, -0.1e-6, z], || false)
                .unwrap();
            let y = prepared
                .field_at(right, [0.2e-6, -0.1e-6, z], || false)
                .unwrap();
            for j in 0..2 {
                complex_near(x.electric_v_m[j], y.electric_v_m[j], 5e-13);
                complex_near(x.magnetic_a_m[j], y.magnetic_a_m[j], 5e-15);
            }
            near(x.poynting_w_m2[2], y.poynting_w_m2[2], 5e-15);
        }
    }
    assert!(
        stack
            .response_at_tangential_index(frequency(399.0), 0.4, input[0].electric_v_m, || false)
            .is_err()
    );
    assert!(
        OrientedUniaxialStack::try_new(
            DielectricResponse::Lossless(&outside),
            DielectricResponse::Lossless(&exit),
            &layers,
            UniaxialStackLimits {
                maximum_layers: 1,
                ..limits
            }
        )
        .is_err()
    );
    assert!(matches!(
        PreparedUniaxialStackSpectrum::prepare(
            &stack,
            frequency(399.0),
            &input,
            UniaxialSpectrumLimits {
                maximum_modes: 0,
                ..spectrum_limits
            },
            || false
        ),
        Err(AngularSpectrumError::Invalid(_))
    ));
    assert!(matches!(
        PreparedUniaxialStackSpectrum::prepare(
            &stack,
            frequency(399.0),
            &input,
            spectrum_limits,
            || true
        ),
        Err(AngularSpectrumError::Cancelled)
    ));
}

#[test]
fn coherent_absorbing_crystal_illumination_retains_source_frequency_and_global_fields() {
    use crate::plane_wave::uniaxial::oriented::{
        OrientedCrystalFieldSide, OrientedUniaxialInterface,
        absorbing::{AbsorbingCrystalCoherentIllumination, CrystalIlluminationError},
    };
    let outside = incident(1.0, MU);
    let (o, e) = axes(
        [1.4, 1.6],
        [1.8, 2.0],
        [0.08, 0.2],
        Some(MU),
        [400.0, 800.0],
        4,
    );
    let law = UniaxialOpticalResponse::try_new(&o, &e, SPAN).unwrap();
    let interface = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside),
        law,
        [1.0, 2.0, -0.5],
    )
    .unwrap();
    let omega = frequency(550.0);
    let electric = [C64::new(0.6, -0.2), C64::new(0.1, 0.7)];
    let components = [
        interface
            .prepare_absorbing_field_at_tangential_index_in_plane(
                omega,
                0.5,
                [1.0, 0.0],
                electric,
            )
            .unwrap(),
        interface
            .prepare_absorbing_field_at_tangential_index_in_plane(
                omega,
                -0.4,
                [3.0, 4.0],
                electric,
            )
            .unwrap(),
    ];
    let illumination =
        AbsorbingCrystalCoherentIllumination::try_new(&interface, &components, 2, || false)
            .unwrap();
    assert!(core::ptr::eq(illumination.source(), &interface));
    assert!(core::ptr::eq(
        illumination.components().as_ptr(),
        components.as_ptr()
    ));
    assert_eq!(illumination.component_count(), 2);
    assert_eq!(illumination.omega_rad_s(), Some(omega));
    assert_eq!(
        illumination.source().crystal_response().material().identity(),
        law.material().identity()
    );
    let side = OrientedCrystalFieldSide::Transmitted;
    let point = [0.11e-6, -0.07e-6, 0.23e-6];
    let actual = illumination.field(side, point, || false).unwrap();
    let first = components[0].field(side, point).unwrap();
    let second = components[1].field(side, point).unwrap();
    for j in 0..3 {
        complex_near(
            actual.field.electric_v_m[j],
            first.field.electric_v_m[j] + second.field.electric_v_m[j],
            5e-15,
        );
        complex_near(
            actual.field.magnetic_a_m[j],
            first.field.magnetic_a_m[j] + second.field.magnetic_a_m[j],
            5e-17,
        );
    }
    assert_eq!(
        illumination.poynting_w_m2(side, point, || false).unwrap(),
        actual.field.poynting_w_m2
    );
    let single =
        AbsorbingCrystalCoherentIllumination::try_new(&interface, &components[..1], 1, || false)
            .unwrap();
    assert_eq!(single.field(side, point, || false).unwrap(), first);
    assert!(matches!(
        AbsorbingCrystalCoherentIllumination::try_new(&interface, &components, 1, || false),
        Err(CrystalIlluminationError::Electromagnetic(
            ReflectionError::Invalid(_)
        ))
    ));
    let equivalent_interface = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside),
        law,
        [1.0, 2.0, -0.5],
    )
    .unwrap();
    assert!(matches!(
        AbsorbingCrystalCoherentIllumination::try_new(
            &equivalent_interface,
            &components,
            2,
            || false,
        ),
        Err(CrystalIlluminationError::Electromagnetic(
            ReflectionError::Invalid(_)
        ))
    ));
    // Even one frequency ULP needs a temporal beat/exposure model; proximity
    // is not permission to invent a common monochromatic carrier.
    let mixed_frequencies = [
        interface
            .prepare_absorbing_field(omega, 0.8, electric)
            .unwrap(),
        interface
            .prepare_absorbing_field(omega.next_up(), 0.8, electric)
            .unwrap(),
    ];
    assert!(matches!(
        AbsorbingCrystalCoherentIllumination::try_new(
            &interface,
            &mixed_frequencies,
            2,
            || false,
        ),
        Err(CrystalIlluminationError::Electromagnetic(
            ReflectionError::Invalid(_)
        ))
    ));
    let dark = AbsorbingCrystalCoherentIllumination::try_new(&interface, &[], 0, || false)
        .unwrap();
    assert_eq!(dark.omega_rad_s(), None);
    assert_eq!(
        dark.field(side, point, || false).unwrap().absorbed_power_density_w_m3,
        0.0
    );
    assert!(dark.field(side, [0.0, 0.0, -1.0], || false).is_err());
    for cancel_at in 1..=4 {
        let mut polls = 0;
        let result =
            AbsorbingCrystalCoherentIllumination::try_new(&interface, &components, 2, || {
                polls += 1;
                polls == cancel_at
            });
        assert!(matches!(result, Err(CrystalIlluminationError::Cancelled)));
        assert_eq!(polls, cancel_at);
    }
    assert!(matches!(
        illumination.field(side, point, || true),
        Err(CrystalIlluminationError::Cancelled)
    ));
    assert!(matches!(
        illumination.poynting_w_m2(side, point, || true),
        Err(CrystalIlluminationError::Cancelled)
    ));
}

#[test]
fn coherent_crystal_box_retains_original_material_volume_and_resource_admission() {
    use crate::plane_wave::uniaxial::oriented::{
        OrientedUniaxialInterface,
        absorbing::{
            AbsorbingCrystalBox, AbsorbingCrystalBoxLimits,
            AbsorbingCrystalCoherentBoxAbsorption, AbsorbingCrystalCoherentIllumination,
            CrystalIlluminationError,
        },
    };
    let outside = incident(1.0, MU);
    let (o, e) = axes(
        [1.5; 2], [1.9; 2], [0.08, 0.2], Some(MU), [400.0, 800.0], 3,
    );
    let law = UniaxialOpticalResponse::try_new(&o, &e, SPAN).unwrap();
    let interface = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside), law, [0.0, 0.0, 1.0],
    ).unwrap();
    let components = [
        interface.prepare_absorbing_field_at_tangential_index(
            frequency(550.0), 0.25, [C64::ZERO, C64::ONE],
        ).unwrap(),
        interface.prepare_absorbing_field_at_tangential_index(
            frequency(550.0), 0.75, [C64::ZERO, C64::new(0.4, 0.9)],
        ).unwrap(),
    ];
    let illumination = AbsorbingCrystalCoherentIllumination::try_new(
        &interface, &components, 2, || false,
    ).unwrap();
    let region = AbsorbingCrystalBox::try_new(
        [0.13e-6, -0.17e-6, 0.1e-6], [0.31e-6, 0.4e-6, 0.9e-6],
    ).unwrap();
    let mut limits = AbsorbingCrystalBoxLimits {
        maximum_mode_pairs: 16,
        maximum_workspace_bytes: 1 << 20,
    };
    let result = illumination.box_absorption(region, limits, || false).unwrap();
    assert!(core::ptr::eq(result.illumination(), &illumination));
    assert!(core::ptr::eq(result.source(), &interface));
    assert_eq!(result.source().crystal_response().material().identity(), law.material().identity());
    assert_eq!(result.region(), region);
    assert_eq!(result.volume_m3(), region.volume_m3());
    assert_eq!(result.mode_pair_count(), 16);
    assert!(result.absorbed_power_w() > 0.0);
    assert!(
        result.workspace_bytes()
            >= core::mem::size_of::<AbsorbingCrystalCoherentBoxAbsorption<'_, '_, '_, '_>>()
    );
    let faces = result.outward_face_power_w();
    let scale = faces.iter().map(|v| v.abs()).sum::<f64>();
    near(result.residual_w(), 0.0, 5e-12 * scale);
    assert!(faces[0] < 0.0);
    assert!(faces[1] > 0.0);
    assert!(faces[4] < 0.0);
    assert!(faces[5] > 0.0);
    limits.maximum_workspace_bytes = result.workspace_bytes();
    let repeated = illumination.box_absorption(region, limits, || false).unwrap();
    assert_eq!(repeated.absorbed_power_w(), result.absorbed_power_w());
    assert_eq!(repeated.outward_face_power_w(), faces);
    limits.maximum_workspace_bytes -= 1;
    assert!(illumination.box_absorption(region, limits, || false).is_err());
    limits.maximum_workspace_bytes = 1 << 20;
    limits.maximum_mode_pairs = 15;
    assert!(illumination.box_absorption(region, limits, || false).is_err());
    limits.maximum_mode_pairs = 16;
    assert!(matches!(
        illumination.box_absorption(region, limits, || true),
        Err(CrystalIlluminationError::Cancelled)
    ));
    let dark = AbsorbingCrystalCoherentIllumination::try_new(&interface, &[], 0, || false).unwrap();
    limits.maximum_mode_pairs = 0;
    let result = dark.box_absorption(region, limits, || false).unwrap();
    assert_eq!(result.absorbed_power_w(), 0.0);
    assert_eq!(result.outward_face_power_w(), [0.0; 6]);
    assert_eq!(result.mode_pair_count(), 0);
}

#[path = "uniaxial_tests/partial_coherence.rs"]
mod partial_coherence;

#[path = "uniaxial_tests/absorbing_slab_box.rs"]
mod absorbing_slab_box;


#[test]
fn g1_original_crystal_group_slopes_reach_native_carrier_and_independent_dispersion_gradient() {
    let source_incident = incident(1.3, 0.8 * MU);
    let (o, e) = axes([2.8, 1.2], [3.1, 1.4], [0.0; 2], Some(1.7 * MU), [400.0, 800.0], 2);
    let law = UniaxialOpticalResponse::try_new(&o, &e, SPAN).unwrap();
    let omega = frequency(600.0);
    // Original n(lambda) lines have intercepts 4.4 and 4.8. Their group
    // indices are these intercepts, not their carrier phase-index samples.
    let groups = law.group_indices(omega).unwrap();
    near(groups[0], 4.4, 3e-14);
    near(groups[1], 4.8, 3e-14);
    fn slope(f: impl Fn(f64) -> f64, x: f64, h: f64) -> f64 {
        (f(x - 2.0 * h) - 8.0 * f(x - h) + 8.0 * f(x + h) - f(x + 2.0 * h))
            / (12.0 * h)
    }
    for authored_axis in [[1.0_f64, 2.0], [-1.0, 2.0]] {
        let length = authored_axis[0].hypot(authored_axis[1]);
        let [ax, az] = authored_axis.map(|v| v / length);
        let interface = PrincipalPlaneUniaxialInterface::try_new(
            DielectricResponse::Lossless(&source_incident), law, authored_axis,
        ).unwrap();
        for p in [-0.65, 0.65] {
            for polarization in [Polarization::Te, Polarization::Tm] {
                let prepared = interface.group_velocity_response_at_tangential_index(
                    omega, p, polarization,
                ).unwrap();
                assert!(core::ptr::eq(prepared.source(), &interface));
                assert!(core::ptr::eq(prepared.source().crystal_response().ordinary_source(), &o));
                assert_eq!(prepared.angular_frequency_rad_s(), omega);
                assert_eq!(prepared.carrier(),
                    interface.response_at_tangential_index(omega, p, polarization).unwrap());
                assert_eq!(prepared.principal_group_indices(), groups);
                let carrier = prepared.carrier();
                let k = [carrier.tangential_wavenumber_rad_m, carrier.transmitted_normal_wavenumber_rad_m.re];
                let velocity = prepared.group_velocity_m_s().unwrap();
                // Independently perturb the full k/omega dispersion relation,
                // evaluating the authored wavelength lines without a production jet.
                let dispersion = |w: f64, k: [f64; 2]| {
                    let lambda = core::f64::consts::TAU * C0 * 1e9 / w;
                    let no = 2.8 + (1.2 - 2.8) * ((lambda - 400.0) / 400.0);
                    let ne = 3.1 + (1.4 - 3.1) * ((lambda - 400.0) / 400.0);
                    let transverse = match polarization {
                        Polarization::Te => (k[0] * k[0] + k[1] * k[1]) / (no * no),
                        Polarization::Tm => {
                            let parallel = ax * k[0] + az * k[1];
                            let perpendicular = -az * k[0] + ax * k[1];
                            parallel * parallel / (no * no) + perpendicular * perpendicular / (ne * ne)
                        }
                    };
                    transverse - (w / C0) * (w / C0)
                };
                let d_omega = slope(|w| dispersion(w, k), omega, omega * 1e-5);
                let step_k = k[0].hypot(k[1]) * 1e-5;
                for j in 0..2 {
                    let d_k = slope(|value| {
                        let mut changed = k;
                        changed[j] = value;
                        dispersion(omega, changed)
                    }, k[j], step_k);
                    near(velocity[j], -d_k / d_omega, 5e-10 * C0);
                }
            }
        }
        let by_angle = interface.group_velocity_response(omega, 0.8, Polarization::Tm).unwrap();
        let by_tangent = interface.group_velocity_response_at_tangential_index(
            omega, 1.3 * 0.6, Polarization::Tm,
        ).unwrap();
        for j in 0..2 {
            near(by_angle.group_velocity_m_s().unwrap()[j],
                by_tangent.group_velocity_m_s().unwrap()[j], 2e-13 * C0);
        }
    }
}

#[test]
fn g0_crystal_group_sources_keep_loss_band_and_zero_or_anomalous_slope_semantics() {
    let source_incident = incident(1.0, MU);
    let (o, e) = axes([1.4, 1.6], [1.8, 2.0], [0.0; 2], Some(MU), [400.0, 800.0], 0);
    let law = UniaxialOpticalResponse::try_new(&o, &e, SPAN).unwrap();
    let interface = PrincipalPlaneUniaxialInterface::try_new(
        DielectricResponse::Lossless(&source_incident), law, [0.6, 0.8],
    ).unwrap();
    for wavelength in [400.0, 800.0] {
        // Inclusive monochromatic support is not a two-sided derivative domain.
        assert!(interface.response(frequency(wavelength), 0.85, Polarization::Tm).is_ok());
        assert!(law.group_indices(frequency(wavelength)).is_err());
        assert!(interface.group_velocity_response(frequency(wavelength), 0.85, Polarization::Tm).is_err());
    }
    let (loss_o, loss_e) = axes([1.4, 1.6], [1.8, 2.0], [0.01, 0.0], Some(MU), [400.0, 800.0], 0);
    let lossy = UniaxialOpticalResponse::try_new(&loss_o, &loss_e, SPAN).unwrap();
    assert!(lossy.group_indices(frequency(600.0)).is_err());
    for (samples, expected_group) in [([1.0, 2.0], 0.0), ([0.5, 1.5], -0.5)] {
        let (o, e) = axes(samples, samples, [0.0; 2], Some(MU), [400.0, 800.0], 0);
        let law = UniaxialOpticalResponse::try_new(&o, &e, SPAN).unwrap();
        let interface = PrincipalPlaneUniaxialInterface::try_new(
            DielectricResponse::Lossless(&source_incident), law, [0.6, 0.8],
        ).unwrap();
        for polarization in [Polarization::Te, Polarization::Tm] {
            let prepared = interface.group_velocity_response(frequency(600.0), 0.85, polarization).unwrap();
            for value in prepared.principal_group_indices() {
                near(value, expected_group, 2e-15);
            }
            if expected_group == 0.0 {
                assert_eq!(prepared.group_velocity_m_s(), None);
            } else {
                // Preserve the signed local spectral slope without clipping it
                // into a positive ray-energy or subluminal signal-speed claim.
                assert!(prepared.group_velocity_m_s().unwrap()[1] < 0.0);
                assert!(prepared.carrier().energy_direction.unwrap()[1] > 0.0);
            }
        }
    }
}

#[test]
fn g1_oriented_group_transport_retains_original_sources_and_the_exact_prepared_carrier() {
    use crate::plane_wave::uniaxial::oriented::OrientedUniaxialInterface;

    let source_incident = incident(1.3, 0.8 * MU);
    let (o, e) = axes(
        [2.8, 1.2],
        [3.1, 1.4],
        [0.0; 2],
        Some(1.7 * MU),
        [400.0, 800.0],
        2,
    );
    let law = UniaxialOpticalResponse::try_new(&o, &e, SPAN).unwrap();
    let interface = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&source_incident),
        law,
        [1.0, 2.0, 3.0],
    )
    .unwrap();
    let omega = frequency(600.0);
    let electric = [C64::new(1.0, 0.2), C64::new(0.4, -0.3)];
    for p in [-0.65, 0.65] {
        let prepared = interface
            .prepare_field_at_tangential_index(omega, p, electric)
            .unwrap();
        let boundary = *prepared.boundary_response();
        let from_prepared = prepared.group_velocity_response().unwrap();
        let from_interface = interface
            .group_velocity_response_at_tangential_index(omega, p, electric)
            .unwrap();
        assert!(core::ptr::eq(from_prepared.source(), &interface));
        assert!(core::ptr::eq(
            from_prepared.source().crystal_response().ordinary_source(),
            &o,
        ));
        assert!(core::ptr::eq(
            from_prepared.source().crystal_response().extraordinary_source(),
            &e,
        ));
        assert_eq!(from_prepared.omega_rad_s(), omega);
        assert_eq!(*from_prepared.boundary_response(), boundary);
        assert_eq!(*from_interface.boundary_response(), boundary);
        assert_eq!(*prepared.boundary_response(), boundary);
        assert_eq!(
            from_prepared.group_velocities_m_s(),
            from_interface.group_velocities_m_s(),
        );
        near(from_prepared.principal_group_indices()[0], 4.4, 3e-14);
        near(from_prepared.principal_group_indices()[1], 4.8, 3e-14);
        let [ordinary, extraordinary] = from_prepared.group_velocities_m_s();
        let ordinary = ordinary.unwrap();
        let extraordinary = extraordinary.unwrap();
        assert!(ordinary[2] > 0.0 && extraordinary[2] > 0.0);
        near(ordinary[1], 0.0, 2e-14 * C0);
        assert!(extraordinary[1].abs() > 1e-3 * C0);
        assert_eq!(boundary.transmitted[1].phase_direction.unwrap()[1], 0.0);
    }
    let by_angle = interface.group_velocity_response(omega, 0.8, electric).unwrap();
    let by_tangent = interface
        .group_velocity_response_at_tangential_index(omega, 1.3 * 0.6, electric)
        .unwrap();
    for m in 0..2 {
        for j in 0..3 {
            near(
                by_angle.group_velocities_m_s()[m].unwrap()[j],
                by_tangent.group_velocities_m_s()[m].unwrap()[j],
                2e-13 * C0,
            );
        }
    }
    // Preparing monochromatic fields at an original source endpoint does
    // not grant a two-sided derivative domain to the added group query.
    for wavelength in [400.0, 800.0] {
        let prepared = interface.prepare_field(frequency(wavelength), 0.8, electric).unwrap();
        assert!(prepared.group_velocity_response().is_err());
        assert!(interface.group_velocity_response(frequency(wavelength), 0.8, electric).is_err());
    }
}

#[test]
fn original_sources_preserve_matched_extraordinary_power_at_grazing_incidence() {
    use crate::plane_wave::uniaxial::oriented::OrientedUniaxialInterface;
    use fs_math::det;

    let outside = incident(1.5, MU);
    let omega = frequency(550.0);
    // Match the ACTUAL incident law, including its native permittivity-to-index
    // rounding. A separately rounded nominal 1.5 is a different cutoff.
    let epsilon = DielectricResponse::Lossless(&outside)
        .lossless_permittivity(omega)
        .unwrap();
    let ni = (C0 * det::sqrt(MU)) * det::sqrt(epsilon);
    let (ordinary, extraordinary) = axes([1.0; 2], [ni; 2], [0.0; 2], Some(MU), [400.0, 800.0], 0);
    let crystal = UniaxialOpticalResponse::try_new(&ordinary, &extraordinary, SPAN).unwrap();
    let interface = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside),
        crystal,
        [0.0, 0.0, 1.0],
    )
    .unwrap();
    assert!(core::ptr::eq(
        interface.crystal_response().extraordinary_source(),
        &extraordinary,
    ));
    for cosine in [1e-10, 1e-16, 1e-20] {
        let result = interface
            .response(omega, cosine, [C64::from_re(cosine), C64::ZERO])
            .unwrap();
        let transmitted = result.transmitted[1];
        let kz = (omega / C0) * cosine;
        assert_eq!(transmitted.wave_kind, UniaxialWaveKind::Propagating);
        near(transmitted.normal_wavenumber_rad_m.re, kz, 4e-14 * kz);
        assert_eq!(transmitted.normal_wavenumber_rad_m.im, 0.0);
        assert!(transmitted.normal_power_w_m2 > 0.0);
        near(
            result.transmitted_power_fraction,
            4.0 * ni / ((ni + 1.0) * (ni + 1.0)),
            4e-14,
        );
        near(
            result.reflected_power_fraction + result.transmitted_power_fraction,
            1.0,
            4e-14,
        );
    }
}

#[test]
fn original_homogeneous_sources_keep_grazing_phase_through_slabs_and_stacks() {
    use crate::plane_wave::uniaxial::oriented::{
        OrientedUniaxialInterface,
        slab::{OrientedUniaxialSlab, SlabRegion},
        stack::{OrientedUniaxialStack, StackRegion, UniaxialStackLayer, UniaxialStackLimits},
    };
    use fs_math::det;

    let outside = incident(1.5, MU);
    let epsilon = DielectricResponse::Lossless(&outside)
        .lossless_permittivity(frequency(550.0))
        .unwrap();
    let ni = (C0 * det::sqrt(MU)) * det::sqrt(epsilon);
    let (ordinary, extraordinary) = axes([ni; 2], [ni; 2], [0.0; 2], Some(MU), [400.0, 800.0], 0);
    let crystal = UniaxialOpticalResponse::try_new(&ordinary, &extraordinary, SPAN).unwrap();
    let entrance = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside),
        crystal,
        [1.0, 2.0, 3.0],
    )
    .unwrap();
    for wavelength in [450.0, 650.0] {
        let omega = frequency(wavelength);
        for cosine in [1e-8, 1e-10, 1e-12] {
            let kz = (omega / C0) * ni * cosine;
            // Keep a finite, measurable normal phase at every grazing cosine.
            let first_depth = 0.23 / kz;
            let second_depth = 0.37 / kz;
            let depth = first_depth + second_depth;
            let electric = [C64::new(0.7 * cosine, -0.2 * cosine), C64::new(0.3, 0.4)];
            let slab = OrientedUniaxialSlab::try_new(
                &entrance,
                DielectricResponse::Lossless(&outside),
                depth,
            )
            .unwrap();
            let layers = [
                UniaxialStackLayer::try_new(crystal, [1.0, 2.0, 3.0], first_depth).unwrap(),
                UniaxialStackLayer::try_new(crystal, [2.0, -1.0, 1.0], second_depth).unwrap(),
            ];
            let stack = OrientedUniaxialStack::try_new(
                DielectricResponse::Lossless(&outside),
                DielectricResponse::Lossless(&outside),
                &layers,
                UniaxialStackLimits {
                    maximum_layers: 2,
                    maximum_workspace_bytes: 50_000,
                },
            )
            .unwrap();
            let a = slab.response(omega, cosine, electric).unwrap();
            let b = stack.response(omega, cosine, electric, || false).unwrap();
            for (kind, normal, reflected, transmitted) in [
                (
                    a.exit_wave_kind,
                    a.exit_normal_wavenumber_rad_m,
                    a.reflected_power_fraction,
                    a.transmitted_power_fraction,
                ),
                (
                    b.exit_wave_kind,
                    b.exit_normal_wavenumber_rad_m,
                    b.reflected_power_fraction,
                    b.transmitted_power_fraction,
                ),
            ] {
                assert_eq!(kind, UniaxialWaveKind::Propagating);
                near(normal.re, kz, 4e-14 * kz);
                assert_eq!(normal.im, 0.0);
                near(reflected, 0.0, 4e-12);
                near(transmitted, 1.0, 4e-12);
            }
            for (index, layer) in b.layers.iter().enumerate() {
                let z = 0.5 * (layer.entrance_z_m + layer.exit_z_m);
                let field_a = a.field_at(SlabRegion::Crystal, [0.0, 0.0, z]).unwrap();
                let field_b = b
                    .field_at(StackRegion::Layer(index), [0.0, 0.0, z])
                    .unwrap();
                let phase = C64::new(det::cos(kz * z), -det::sin(kz * z));
                for component in 0..3 {
                    let expected_e = a.incident_electric_v_m[component] * phase;
                    let expected_h = a.incident_magnetic_a_m[component] * phase;
                    complex_near(field_a.electric_v_m[component], expected_e, 4e-12);
                    complex_near(field_b.electric_v_m[component], expected_e, 4e-12);
                    complex_near(field_a.magnetic_a_m[component], expected_h, 4e-14);
                    complex_near(field_b.magnetic_a_m[component], expected_h, 4e-14);
                }
            }
            // A caller-authored tangential index exactly at ni is critical;
            // it does not carry the positive cosine retained above.
            assert!(
                slab.response_at_tangential_index(omega, ni, electric)
                    .is_err()
            );
            assert!(
                stack
                    .response_at_tangential_index(omega, ni, electric, || false)
                    .is_err()
            );
            assert!(slab.response(omega, 0.0, electric).is_err());
            assert!(stack.response(omega, 0.0, electric, || false).is_err());
        }
    }
}

#[test]
fn crystal_normal_frequency_jets_retain_sources_carriers_and_explicit_frequency_paths() {
    use crate::plane_wave::dispersion::IncidenceFrequencyConstraint;
    let outside = incident(1.1, MU);
    let (ordinary, extraordinary) = axes([1.6; 2], [1.9; 2], [0.0; 2], Some(MU), [400.0, 800.0], 0);
    let crystal = UniaxialOpticalResponse::try_new(&ordinary, &extraordinary, SPAN).unwrap();
    let interface = PrincipalPlaneUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside),
        crystal,
        [0.8, 0.6],
    )
    .unwrap();
    let omega = frequency(550.0);
    for p in [-0.4, 0.4] {
        for polarization in [Polarization::Te, Polarization::Tm] {
            let carrier = interface
                .response_at_tangential_index(omega, p, polarization)
                .unwrap();
            let angle = interface
                .frequency_dispersion_at_tangential_index(
                    omega,
                    p,
                    polarization,
                    IncidenceFrequencyConstraint::FixedAngle,
                )
                .unwrap();
            assert!(core::ptr::eq(angle.source(), &interface));
            assert_eq!(angle.angular_frequency_rad_s(), omega);
            assert_eq!(angle.carrier(), carrier);
            assert_eq!(
                angle.frequency_constraint(),
                IncidenceFrequencyConstraint::FixedAngle
            );
            complex_near(
                angle.transmitted_normal_wavenumber_derivative_s_m(),
                carrier
                    .transmitted_normal_wavenumber_rad_m
                    .scale(1.0 / omega),
                3e-14 / C0,
            );
            assert_eq!(
                angle.transmitted_normal_wavenumber_second_derivative_s2_m(),
                C64::ZERO
            );
            let distance = 1e-6;
            let propagation = angle.transmitted_propagation(distance).unwrap();
            let phase = carrier.transmitted_normal_wavenumber_rad_m.re * distance;
            complex_near(
                propagation.transfer,
                C64::new(fs_math::det::cos(phase), -fs_math::det::sin(phase)),
                3e-14,
            );
            near(
                propagation.group_delay_s,
                angle.transmitted_normal_wavenumber_derivative_s_m().re * distance,
                1e-30,
            );
            assert_eq!(propagation.group_delay_dispersion_s2, 0.0);
            assert!(angle.transmitted_propagation(-distance).is_err());
        }
        let normal = interface
            .frequency_dispersion_at_tangential_index(
                omega,
                p,
                Polarization::Te,
                IncidenceFrequencyConstraint::FixedTransverseWavenumber,
            )
            .unwrap();
        let q = normal.carrier().transmitted_normal_wavenumber_rad_m.re;
        let kx = normal.carrier().tangential_wavenumber_rad_m;
        let k = omega * 1.6 / C0;
        let first = (k / q) * (1.6 / C0);
        let second = -(1.6 / C0).powi(2) * (kx / q).powi(2) / q;
        near(
            normal.transmitted_normal_wavenumber_derivative_s_m().re,
            first,
            4e-14 * first,
        );
        near(
            normal
                .transmitted_normal_wavenumber_second_derivative_s2_m()
                .re,
            second,
            5e-14 * second.abs(),
        );
        assert!(second < 0.0);
    }
    // The carrier remains available at a source endpoint. A local frequency
    // jet cannot invent the missing side of that original spectral domain.
    for wavelength in [400.0, 800.0] {
        assert!(
            interface
                .response(frequency(wavelength), 0.8, Polarization::Tm)
                .is_ok()
        );
        assert!(
            interface
                .frequency_dispersion(
                    frequency(wavelength),
                    0.8,
                    Polarization::Tm,
                    IncidenceFrequencyConstraint::FixedTransverseWavenumber,
                )
                .is_err()
        );
    }
}

#[test]
fn matched_grazing_crystal_dispersion_keeps_the_original_ordinary_bulk_slope() {
    use crate::plane_wave::dispersion::IncidenceFrequencyConstraint;
    let outside = incident(1.5, MU);
    let omega = frequency(550.0);
    let epsilon = DielectricResponse::Lossless(&outside)
        .lossless_permittivity(omega)
        .unwrap();
    let ni = (C0 * fs_math::det::sqrt(MU)) * fs_math::det::sqrt(epsilon);
    let (ordinary, extraordinary) =
        axes([1.2, 1.4], [ni; 2], [0.0; 2], Some(MU), [400.0, 800.0], 0);
    let crystal = UniaxialOpticalResponse::try_new(&ordinary, &extraordinary, SPAN).unwrap();
    let interface = PrincipalPlaneUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside),
        crystal,
        [0.0, 1.0],
    )
    .unwrap();
    let ordinary_law =
        SpectralOpticalResponse::try_new(SpectralOpticalSource::ComplexIndex(&ordinary), SPAN)
            .unwrap();
    let ordinary_bulk = DielectricResponse::SpectralOptical(ordinary_law)
        .bulk_dispersion(omega)
        .unwrap();
    for cosine in [1e-10, 1e-16, 1e-20] {
        let carrier = interface.response(omega, cosine, Polarization::Tm).unwrap();
        let jet = interface
            .frequency_dispersion(
                omega,
                cosine,
                Polarization::Tm,
                IncidenceFrequencyConstraint::FixedAngle,
            )
            .unwrap();
        assert_eq!(jet.carrier(), carrier);
        assert_eq!(carrier.wave_kind, UniaxialWaveKind::Propagating);
        let expected = ordinary_bulk.wavenumber_derivative_s_m().scale(cosine);
        complex_near(
            jet.transmitted_normal_wavenumber_derivative_s_m(),
            expected,
            3e-14 * expected.abs(),
        );
        // The original ordinary n(lambda) is linear. Its exact zero bulk
        // curvature survives the grazing normal-coordinate projection.
        assert_eq!(
            jet.transmitted_normal_wavenumber_second_derivative_s2_m(),
            C64::ZERO
        );
        let propagated = jet.transmitted_propagation(0.003).unwrap();
        near(
            propagated.group_delay_s,
            expected.re * 0.003,
            expected.abs() * 1e-15,
        );
        assert_eq!(propagated.log_amplitude_derivative_s, 0.0);
    }
}
#[test]
fn oriented_frequency_jets_retain_vector_carriers_and_separate_modal_propagation() {
    use crate::plane_wave::{
        dispersion::IncidenceFrequencyConstraint,
        uniaxial::{UniaxialMode, oriented::OrientedUniaxialInterface},
    };
    let outside = incident(1.1, MU);
    let (ordinary, extraordinary) = axes([1.6; 2], [1.9; 2], [0.0; 2], Some(MU), [400.0, 800.0], 0);
    let crystal = UniaxialOpticalResponse::try_new(&ordinary, &extraordinary, SPAN).unwrap();
    let interface = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside),
        crystal,
        [1.0, 2.0, 3.0],
    )
    .unwrap();
    let omega = frequency(550.0);
    let electric = [C64::new(0.7, 0.2), C64::new(0.3, -0.4)];
    for p in [-0.4, 0.4] {
        let carrier = interface
            .response_at_tangential_index(omega, p, electric)
            .unwrap();
        for constraint in [
            IncidenceFrequencyConstraint::FixedAngle,
            IncidenceFrequencyConstraint::FixedTransverseWavenumber,
        ] {
            let jet = interface
                .frequency_dispersion_at_tangential_index(omega, p, electric, constraint)
                .unwrap();
            assert!(core::ptr::eq(jet.source(), &interface));
            assert!(core::ptr::eq(
                jet.source().crystal_response().ordinary_source(),
                &ordinary
            ));
            assert!(core::ptr::eq(
                jet.source().crystal_response().extraordinary_source(),
                &extraordinary
            ));
            assert_eq!(jet.omega_rad_s(), omega);
            assert_eq!(jet.frequency_constraint(), constraint);
            assert_eq!(*jet.boundary_response(), carrier);
            let first = jet.transmitted_normal_wavenumber_derivatives_s_m();
            let second = jet.transmitted_normal_wavenumber_second_derivatives_s2_m();
            assert_ne!(first[0], first[1]);
            for (index, mode) in [UniaxialMode::Ordinary, UniaxialMode::Extraordinary]
                .into_iter()
                .enumerate()
            {
                let distance = 3e-6;
                let propagation = jet.transmitted_propagation(mode, distance).unwrap();
                let phase = carrier.transmitted[index].normal_wavenumber_rad_m.re * distance;
                complex_near(
                    propagation.transfer,
                    C64::new(fs_math::det::cos(phase), -fs_math::det::sin(phase)),
                    4e-14,
                );
                near(propagation.group_delay_s, first[index].re * distance, 1e-30);
                near(
                    propagation.group_delay_dispersion_s2,
                    second[index].re * distance,
                    1e-44,
                );
                assert!(jet.transmitted_propagation(mode, -distance).is_err());
                if constraint == IncidenceFrequencyConstraint::FixedAngle {
                    complex_near(
                        first[index],
                        carrier.transmitted[index]
                            .normal_wavenumber_rad_m
                            .scale(1.0 / omega),
                        4e-14 / C0,
                    );
                    assert_eq!(second[index], C64::ZERO);
                } else {
                    assert!(second[index].re < 0.0);
                }
            }
        }
    }
    // The derivative owner retains the original two-sided spectral domain.
    for wavelength in [400.0, 800.0] {
        assert!(
            interface
                .response(frequency(wavelength), 0.8, electric)
                .is_ok()
        );
        assert!(
            interface
                .frequency_dispersion(
                    frequency(wavelength),
                    0.8,
                    electric,
                    IncidenceFrequencyConstraint::FixedTransverseWavenumber,
                )
                .is_err()
        );
    }
}

#[test]
fn oriented_frequency_jets_require_a_spectral_basis_at_equal_phase_indices() {
    use crate::plane_wave::{
        dispersion::IncidenceFrequencyConstraint, uniaxial::oriented::OrientedUniaxialInterface,
    };
    let outside = incident(1.1, MU);
    let omega = frequency(600.0);
    let electric = [C64::new(0.7, 0.2), C64::new(0.3, -0.4)];
    // Bind an exactly equal carrier phase index from the original table,
    // while retaining the ordinary table's distinct (here zero) bulk slope.
    let (probe_o, probe_e) = axes([1.0, 2.0], [1.5; 2], [0.0; 2], Some(MU), [400.0, 800.0], 0);
    let phase = UniaxialOpticalResponse::try_new(&probe_o, &probe_e, SPAN)
        .unwrap()
        .phase_indices(omega)
        .unwrap()[0];
    let (ordinary, extraordinary) = axes(
        [1.0, 2.0],
        [phase; 2],
        [0.0; 2],
        Some(MU),
        [400.0, 800.0],
        0,
    );
    let crystal = UniaxialOpticalResponse::try_new(&ordinary, &extraordinary, SPAN).unwrap();
    let phases = crystal.phase_indices(omega).unwrap();
    assert_eq!(phases[0], phases[1]);
    let mixed = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside),
        crystal,
        [1.0, 2.0, 3.0],
    )
    .unwrap();
    let in_plane = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside),
        crystal,
        [1.0, 0.0, 3.0],
    )
    .unwrap();
    let principal = PrincipalPlaneUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside),
        crystal,
        [1.0, 3.0],
    )
    .unwrap();
    for constraint in [
        IncidenceFrequencyConstraint::FixedAngle,
        IncidenceFrequencyConstraint::FixedTransverseWavenumber,
    ] {
        assert!(mixed.response(omega, 0.8, electric).is_ok());
        assert!(
            mixed
                .frequency_dispersion(omega, 0.8, electric, constraint)
                .is_err()
        );
        let jet = in_plane
            .frequency_dispersion(omega, 0.8, electric, constraint)
            .unwrap();
        for (index, polarization) in [Polarization::Te, Polarization::Tm].into_iter().enumerate() {
            let expected = principal
                .frequency_dispersion(omega, 0.8, polarization, constraint)
                .unwrap();
            complex_near(
                jet.transmitted_normal_wavenumber_derivatives_s_m()[index],
                expected.transmitted_normal_wavenumber_derivative_s_m(),
                5e-14 / C0,
            );
            complex_near(
                jet.transmitted_normal_wavenumber_second_derivatives_s2_m()[index],
                expected.transmitted_normal_wavenumber_second_derivative_s2_m(),
                5e-14 / (C0 * omega),
            );
        }
    }
    // An identical original local law has an admitted isotropic spectral
    // basis even when the authored optic axis has an out-of-plane component.
    let (isotropic_o, isotropic_e) = axes(
        [1.4, 1.6],
        [1.4, 1.6],
        [0.0; 2],
        Some(MU),
        [400.0, 800.0],
        0,
    );
    let isotropic = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside),
        UniaxialOpticalResponse::try_new(&isotropic_o, &isotropic_e, SPAN).unwrap(),
        [1.0, 2.0, 3.0],
    )
    .unwrap();
    for constraint in [
        IncidenceFrequencyConstraint::FixedAngle,
        IncidenceFrequencyConstraint::FixedTransverseWavenumber,
    ] {
        let jet = isotropic
            .frequency_dispersion(omega, 0.8, electric, constraint)
            .unwrap();
        let first = jet.transmitted_normal_wavenumber_derivatives_s_m();
        let second = jet.transmitted_normal_wavenumber_second_derivatives_s2_m();
        complex_near(first[0], first[1], 5e-14 / C0);
        complex_near(second[0], second[1], 5e-14 / (C0 * omega));
    }
}
#[test]
fn prepared_oriented_frequency_jets_reuse_exact_signed_and_grazing_coordinates() {
    use crate::plane_wave::{
        dispersion::IncidenceFrequencyConstraint,
        uniaxial::oriented::{OrientedCrystalFieldSide, OrientedUniaxialInterface},
    };
    let outside = incident(1.5, MU);
    let omega = frequency(550.0);
    let epsilon = DielectricResponse::Lossless(&outside)
        .lossless_permittivity(omega)
        .unwrap();
    let ni = (C0 * fs_math::det::sqrt(MU)) * fs_math::det::sqrt(epsilon);
    let (ordinary, extraordinary) = axes([ni; 2], [ni; 2], [0.0; 2], Some(MU), [400.0, 800.0], 0);
    let crystal = UniaxialOpticalResponse::try_new(&ordinary, &extraordinary, SPAN).unwrap();
    let interface = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside),
        crystal,
        [1.0, 2.0, 3.0],
    )
    .unwrap();
    for p in [-0.4, 0.4] {
        let electric = [C64::new(0.7, 0.2), C64::new(0.3, -0.4)];
        let prepared = interface
            .prepare_field_at_tangential_index(omega, p, electric)
            .unwrap();
        for constraint in [
            IncidenceFrequencyConstraint::FixedAngle,
            IncidenceFrequencyConstraint::FixedTransverseWavenumber,
        ] {
            let direct = interface
                .frequency_dispersion_at_tangential_index(omega, p, electric, constraint)
                .unwrap();
            let added = prepared.frequency_dispersion(constraint).unwrap();
            assert!(core::ptr::eq(added.source(), &interface));
            assert_eq!(*added.boundary_response(), *prepared.boundary_response());
            assert_eq!(*added.boundary_response(), *direct.boundary_response());
            assert_eq!(
                added.transmitted_normal_wavenumber_derivatives_s_m(),
                direct.transmitted_normal_wavenumber_derivatives_s_m(),
            );
            assert_eq!(
                added.transmitted_normal_wavenumber_second_derivatives_s2_m(),
                direct.transmitted_normal_wavenumber_second_derivatives_s2_m(),
            );
        }
    }
    for cosine in [1e-10, 1e-16, 1e-20] {
        let electric = [C64::new(0.7 * cosine, 0.2 * cosine), C64::new(0.3, -0.4)];
        let prepared = interface.prepare_field(omega, cosine, electric).unwrap();
        let position = [0.0, 0.0, 0.13 / ((omega / C0) * ni * cosine)];
        let field_before = prepared
            .field(OrientedCrystalFieldSide::Transmitted, position)
            .unwrap();
        for constraint in [
            IncidenceFrequencyConstraint::FixedAngle,
            IncidenceFrequencyConstraint::FixedTransverseWavenumber,
        ] {
            let direct = interface
                .frequency_dispersion(omega, cosine, electric, constraint)
                .unwrap();
            let added = prepared.frequency_dispersion(constraint).unwrap();
            assert_eq!(*added.boundary_response(), *prepared.boundary_response());
            assert_eq!(*added.boundary_response(), *direct.boundary_response());
            let first = added.transmitted_normal_wavenumber_derivatives_s_m();
            let second = added.transmitted_normal_wavenumber_second_derivatives_s2_m();
            assert_eq!(
                first,
                direct.transmitted_normal_wavenumber_derivatives_s_m()
            );
            assert_eq!(
                second,
                direct.transmitted_normal_wavenumber_second_derivatives_s2_m()
            );
            for mode in 0..2 {
                assert_eq!(
                    added.boundary_response().transmitted[mode].wave_kind,
                    UniaxialWaveKind::Propagating
                );
                let expected = match constraint {
                    IncidenceFrequencyConstraint::FixedAngle => ni * cosine / C0,
                    IncidenceFrequencyConstraint::FixedTransverseWavenumber => ni / (C0 * cosine),
                };
                near(first[mode].re, expected, 8e-14 * expected);
                assert_eq!(first[mode].im, 0.0);
                if constraint == IncidenceFrequencyConstraint::FixedAngle {
                    assert_eq!(second[mode], C64::ZERO);
                } else {
                    let expected =
                        -ni * (1.0 - cosine) * (1.0 + cosine) / ((C0 * omega) * cosine.powi(3));
                    near(second[mode].re, expected, 1e-13 * expected.abs());
                    assert_eq!(second[mode].im, 0.0);
                }
            }
        }
        // Adding spectral derivatives neither replaces the prepared carrier
        // nor changes its native coherent field propagation.
        assert_eq!(
            prepared
                .field(OrientedCrystalFieldSide::Transmitted, position)
                .unwrap(),
            field_before,
        );
    }
}

#[path = "uniaxial_tests/boundary_dispersion.rs"]
mod boundary_dispersion;

#[path = "uniaxial_tests/oriented_boundary_dispersion.rs"]
mod oriented_boundary_dispersion;
#[path = "uniaxial_tests/slab_thickness.rs"]
mod slab_thickness;
