//! Original-source arbitrary-axis vector boundary frequency controls.
use super::*;
use crate::plane_wave::{
    dispersion::IncidenceFrequencyConstraint,
    uniaxial::{
        UniaxialMode,
        oriented::{OrientedCrystalFieldSide, OrientedUniaxialInterface},
    },
};

#[test]
fn prepared_vector_boundary_jets_retain_original_carrier_sources_and_maxwell_traces() {
    let outside = incident(1.2, MU);
    let (ordinary, extraordinary) = axes(
        [1.4, 1.7],
        [1.9, 2.1],
        [0.0; 2],
        Some(1.3 * MU),
        [400.0, 800.0],
        3,
    );
    let source = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside),
        UniaxialOpticalResponse::try_new(&ordinary, &extraordinary, SPAN).unwrap(),
        [1.0, 2.0, -0.7],
    )
    .unwrap();
    let omega = frequency(550.0);
    let electric = [C64::new(0.6, -0.2), C64::new(0.1, 0.7)];
    for p in [-0.4, 0.4] {
        let prepared = source
            .prepare_field_at_tangential_index(omega, p, electric)
            .unwrap();
        let field_before = prepared
            .field(OrientedCrystalFieldSide::Transmitted, [1e-7, 2e-7, 3e-7])
            .unwrap();
        for constraint in [
            IncidenceFrequencyConstraint::FixedAngle,
            IncidenceFrequencyConstraint::FixedTransverseWavenumber,
        ] {
            let direct = source
                .boundary_frequency_dispersion_at_tangential_index(omega, p, electric, constraint)
                .unwrap();
            let added = prepared.boundary_frequency_dispersion(constraint).unwrap();
            assert!(core::ptr::eq(added.source(), &source));
            assert!(core::ptr::eq(
                added.source().crystal_response().ordinary_source(),
                &ordinary,
            ));
            assert!(core::ptr::eq(
                added.source().crystal_response().extraordinary_source(),
                &extraordinary,
            ));
            assert_eq!(added.omega_rad_s(), omega);
            assert_eq!(added.frequency_constraint(), constraint);
            assert_eq!(*added.boundary_response(), *prepared.boundary_response());
            assert_eq!(*added.boundary_response(), *direct.boundary_response());
            assert_eq!(added.incident_field(), direct.incident_field());
            assert_eq!(added.reflected_field(), direct.reflected_field());
            for (index, mode) in [UniaxialMode::Ordinary, UniaxialMode::Extraordinary]
                .into_iter()
                .enumerate()
            {
                assert_eq!(
                    added.transmitted_field(mode),
                    direct.transmitted_field(mode)
                );
                assert_eq!(
                    added.transmitted_field(mode).electric_v_m(),
                    prepared.boundary_response().transmitted[index].electric_v_m,
                );
                assert_eq!(
                    added.transmitted_field(mode).magnetic_a_m(),
                    prepared.boundary_response().transmitted[index].magnetic_a_m,
                );
            }
            let original_normal = prepared.frequency_dispersion(constraint).unwrap();
            assert_eq!(
                added
                    .normal_dispersion()
                    .transmitted_normal_wavenumber_derivatives_s_m(),
                original_normal.transmitted_normal_wavenumber_derivatives_s_m(),
            );
            assert_eq!(
                added
                    .normal_dispersion()
                    .transmitted_normal_wavenumber_second_derivatives_s2_m(),
                original_normal.transmitted_normal_wavenumber_second_derivatives_s2_m(),
            );
            let incident = added.incident_field();
            let reflected = added.reflected_field();
            let transmitted = added.transmitted_boundary_field().unwrap();
            for order in 0..2 {
                for component in 0..2 {
                    assert_eq!(
                        incident.scaled_electric_derivatives_v_m()[order][component],
                        C64::ZERO
                    );
                    for (a, b, c) in [
                        (
                            incident.scaled_electric_derivatives_v_m(),
                            reflected.scaled_electric_derivatives_v_m(),
                            transmitted.scaled_electric_derivatives_v_m(),
                        ),
                        (
                            incident.scaled_magnetic_derivatives_a_m(),
                            reflected.scaled_magnetic_derivatives_a_m(),
                            transmitted.scaled_magnetic_derivatives_a_m(),
                        ),
                    ] {
                        complex_near(
                            a[order][component] + b[order][component],
                            c[order][component],
                            2e-11 * c[order][component].abs().max(1e-3),
                        );
                    }
                }
            }
        }
        assert_eq!(
            prepared
                .field(OrientedCrystalFieldSide::Transmitted, [1e-7, 2e-7, 3e-7],)
                .unwrap(),
            field_before
        );
    }
}

#[test]
fn vector_principal_sections_match_original_scalar_fresnel_jets_including_tir() {
    let outside = incident(2.4, MU);
    let (ordinary, extraordinary) = axes(
        [1.4, 1.6],
        [1.8, 2.1],
        [0.0; 2],
        Some(1.2 * MU),
        [400.0, 800.0],
        0,
    );
    let law = UniaxialOpticalResponse::try_new(&ordinary, &extraordinary, SPAN).unwrap();
    let scalar = PrincipalPlaneUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside),
        law,
        [1.0, 0.7],
    )
    .unwrap();
    let vector = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside),
        law,
        [1.0, 0.0, 0.7],
    )
    .unwrap();
    let omega = frequency(550.0);
    for cosine in [0.95, 0.6] {
        for (polarization, mode, component, input) in [
            (
                Polarization::Te,
                UniaxialMode::Ordinary,
                1,
                [C64::ZERO, C64::ONE],
            ),
            (
                Polarization::Tm,
                UniaxialMode::Extraordinary,
                0,
                [C64::ONE, C64::ZERO],
            ),
        ] {
            for constraint in [
                IncidenceFrequencyConstraint::FixedAngle,
                IncidenceFrequencyConstraint::FixedTransverseWavenumber,
            ] {
                let expected = scalar
                    .boundary_frequency_dispersion(omega, cosine, polarization, constraint)
                    .unwrap();
                let actual = vector
                    .boundary_frequency_dispersion(omega, cosine, input, constraint)
                    .unwrap();
                complex_near(
                    actual.reflected_field().electric_v_m()[component],
                    expected.reflected_tangential_e().amplitude(),
                    2e-12,
                );
                for order in 0..2 {
                    complex_near(
                        actual.reflected_field().scaled_electric_derivatives_v_m()[order]
                            [component],
                        expected.reflected_tangential_e().scaled_derivatives()[order],
                        2e-10,
                    );
                    complex_near(
                        actual
                            .transmitted_field(mode)
                            .scaled_electric_derivatives_v_m()[order][component],
                        expected.transmitted_tangential_e().scaled_derivatives()[order],
                        2e-10,
                    );
                }
                if expected.carrier().wave_kind == UniaxialWaveKind::Evanescent {
                    let phase = actual
                        .reflected_field()
                        .electric_component_logarithmic_dispersion()[component]
                        .unwrap();
                    let reference = expected
                        .reflected_tangential_e()
                        .logarithmic_dispersion()
                        .unwrap();
                    near(
                        phase.group_delay_s * omega,
                        reference.group_delay_s * omega,
                        2e-10,
                    );
                    near(
                        (phase.group_delay_dispersion_s2 * omega) * omega,
                        (reference.group_delay_dispersion_s2 * omega) * omega,
                        2e-10,
                    );
                }
            }
        }
    }
}

#[test]
fn boundary_jet_admission_preserves_normal_only_fixed_kx_and_original_endpoints() {
    let (_, incoming) = axes([1.0; 2], [1.1, 1.3], [0.0; 2], Some(MU), [500.0, 700.0], 0);
    let (ordinary, extraordinary) = axes(
        [1.4, 1.6],
        [1.8, 2.0],
        [0.0; 2],
        Some(MU),
        [400.0, 800.0],
        0,
    );
    let incident = DielectricResponse::SpectralOptical(
        SpectralOpticalResponse::try_new(SpectralOpticalSource::ComplexIndex(&incoming), SPAN)
            .unwrap(),
    );
    let source = OrientedUniaxialInterface::try_new(
        incident,
        UniaxialOpticalResponse::try_new(&ordinary, &extraordinary, SPAN).unwrap(),
        [1.0, 2.0, 3.0],
    )
    .unwrap();
    let input = [C64::new(0.7, 0.2), C64::new(0.3, -0.4)];
    for wavelength in [500.0, 700.0] {
        let omega = frequency(wavelength);
        let prepared = source.prepare_field(omega, 0.8, input).unwrap();
        assert!(
            prepared
                .frequency_dispersion(IncidenceFrequencyConstraint::FixedTransverseWavenumber,)
                .is_ok()
        );
        assert!(
            prepared
                .boundary_frequency_dispersion(
                    IncidenceFrequencyConstraint::FixedTransverseWavenumber,
                )
                .is_err()
        );
    }
    let outside = super::incident(1.1, MU);
    let source = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside),
        UniaxialOpticalResponse::try_new(&ordinary, &extraordinary, SPAN).unwrap(),
        [1.0, 2.0, 3.0],
    )
    .unwrap();
    for wavelength in [400.0, 800.0] {
        let prepared = source
            .prepare_field(frequency(wavelength), 0.8, input)
            .unwrap();
        for constraint in [
            IncidenceFrequencyConstraint::FixedAngle,
            IncidenceFrequencyConstraint::FixedTransverseWavenumber,
        ] {
            assert!(prepared.boundary_frequency_dispersion(constraint).is_err());
        }
    }
}

#[test]
fn matched_prepared_grazing_boundary_preserves_zero_fields_and_incident_admittance_jets() {
    let outside = incident(1.5, MU);
    let omega = frequency(550.0);
    let epsilon = DielectricResponse::Lossless(&outside)
        .lossless_permittivity(omega)
        .unwrap();
    let index = (C0 * fs_math::det::sqrt(MU)) * fs_math::det::sqrt(epsilon);
    let (ordinary, extraordinary) = axes(
        [index; 2],
        [index; 2],
        [0.0; 2],
        Some(MU),
        [400.0, 800.0],
        0,
    );
    let source = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside),
        UniaxialOpticalResponse::try_new(&ordinary, &extraordinary, SPAN).unwrap(),
        [1.0, 2.0, 3.0],
    )
    .unwrap();
    for cosine in [0.8, 1e-10, 1e-20] {
        let prepared = source
            .prepare_field(omega, cosine, [C64::ZERO, C64::ONE])
            .unwrap();
        let response = prepared
            .boundary_frequency_dispersion(IncidenceFrequencyConstraint::FixedTransverseWavenumber)
            .unwrap();
        let reflection = response.reflected_field();
        assert_eq!(reflection.electric_v_m(), [C64::ZERO; 3]);
        assert_eq!(
            reflection.scaled_electric_derivatives_v_m(),
            [[C64::ZERO; 3]; 2]
        );
        assert_eq!(
            reflection.electric_component_logarithmic_dispersion(),
            [None; 3]
        );
        let transmitted = response.transmitted_field(UniaxialMode::Ordinary);
        assert_eq!(
            transmitted.scaled_electric_derivatives_v_m(),
            [[C64::ZERO; 3]; 2]
        );
        let h = response.incident_field().magnetic_a_m();
        let first = response.incident_field().scaled_magnetic_derivatives_a_m()[0];
        let second = response.incident_field().scaled_magnetic_derivatives_a_m()[1];
        let ratio = ((1.0 - cosine) * (1.0 + cosine) / cosine) / cosine;
        complex_near(first[0], h[0].scale(ratio), 2e-13 * first[0].abs());
        complex_near(
            second[0],
            h[0].scale(-3.0 * ratio - ratio * ratio),
            3e-13 * second[0].abs(),
        );
        complex_near(first[2], -h[2], 2e-14 * h[2].abs());
        complex_near(second[2], h[2].scale(2.0), 2e-14 * h[2].abs());
    }
}

#[path = "oriented_boundary_dispersion/spatial.rs"]
mod spatial;
