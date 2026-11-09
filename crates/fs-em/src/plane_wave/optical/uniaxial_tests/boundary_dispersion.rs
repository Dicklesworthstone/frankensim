//! Source-owned crystal boundary dispersion controls; authored numerical laws.
use super::*;
use crate::plane_wave::dispersion::{FresnelAmplitudeDispersion, IncidenceFrequencyConstraint};

fn amplitude_near(actual: FresnelAmplitudeDispersion, expected: FresnelAmplitudeDispersion) {
    complex_near(actual.amplitude(), expected.amplitude(), 5e-13);
    for (a, b) in actual
        .scaled_derivatives()
        .into_iter()
        .zip(expected.scaled_derivatives())
    {
        complex_near(a, b, 4e-11 * b.abs().max(1.0));
    }
}

#[test]
fn original_isotropic_boundary_owner_matches_crystal_fields_phase_and_power_jets() {
    let outside = incident(1.8, MU);
    let (ordinary, extraordinary) = axes(
        [1.3, 1.7],
        [1.3, 1.7],
        [0.0; 2],
        Some(1.6 * MU),
        [400.0, 800.0],
        0,
    );
    let crystal = UniaxialOpticalResponse::try_new(&ordinary, &extraordinary, SPAN).unwrap();
    let original = PlanarDielectricInterface::try_new(
        &outside,
        DielectricResponse::SpectralOptical(
            SpectralOpticalResponse::try_new(SpectralOpticalSource::ComplexIndex(&ordinary), SPAN)
                .unwrap(),
        ),
    )
    .unwrap();
    let interface = PrincipalPlaneUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside),
        crystal,
        [0.37, -0.81],
    )
    .unwrap();
    for wavelength in [550.0, 650.0] {
        let omega = frequency(wavelength);
        for cosine in [0.2, 0.9] {
            for polarization in [Polarization::Te, Polarization::Tm] {
                for constraint in [
                    IncidenceFrequencyConstraint::FixedAngle,
                    IncidenceFrequencyConstraint::FixedTransverseWavenumber,
                ] {
                    let native = interface.response(omega, cosine, polarization).unwrap();
                    let expected = original
                        .frequency_dispersion(omega, cosine, polarization, constraint)
                        .unwrap();
                    let actual = interface
                        .boundary_frequency_dispersion(omega, cosine, polarization, constraint)
                        .unwrap();
                    assert!(core::ptr::eq(actual.source(), &interface));
                    assert_eq!(actual.angular_frequency_rad_s(), omega);
                    assert_eq!(actual.frequency_constraint(), constraint);
                    assert_eq!(actual.carrier(), native);
                    amplitude_near(
                        actual.reflected_tangential_e(),
                        expected.reflected_tangential_e(),
                    );
                    amplitude_near(
                        actual.transmitted_tangential_e(),
                        expected.transmitted_tangential_e(),
                    );
                    amplitude_near(
                        actual.transmitted_conjugate_h(),
                        expected.transmitted_conjugate_h(),
                    );
                    let upgraded = interface
                        .frequency_dispersion(omega, cosine, polarization, constraint)
                        .unwrap()
                        .boundary_frequency_dispersion()
                        .unwrap();
                    assert_eq!(upgraded.carrier(), native);
                    assert_eq!(
                        upgraded.reflected_tangential_e(),
                        actual.reflected_tangential_e()
                    );
                    assert_eq!(
                        upgraded.transmitted_tangential_e(),
                        actual.transmitted_tangential_e()
                    );
                    assert_eq!(
                        upgraded.transmitted_conjugate_h(),
                        actual.transmitted_conjugate_h()
                    );
                    let reflected = actual.reflected_power_fraction().unwrap();
                    let transmitted = actual.transmitted_power_fraction().unwrap();
                    assert_eq!(reflected.fraction(), native.reflected_power_fraction);
                    assert_eq!(transmitted.fraction(), native.transmitted_power_fraction);
                    near(
                        reflected.scaled_derivative() + transmitted.scaled_derivative(),
                        0.0,
                        5e-12,
                    );
                    if native.wave_kind == UniaxialWaveKind::Evanescent {
                        assert_eq!(reflected.scaled_derivative(), 0.0);
                        assert_eq!(transmitted.scaled_derivative(), 0.0);
                        let phase = actual
                            .reflected_tangential_e()
                            .logarithmic_dispersion()
                            .unwrap();
                        let reference = expected
                            .reflected_tangential_e()
                            .logarithmic_dispersion()
                            .unwrap();
                        near(
                            phase.group_delay_s * omega,
                            reference.group_delay_s * omega,
                            3e-11,
                        );
                        near(
                            (phase.group_delay_dispersion_s2 * omega) * omega,
                            (reference.group_delay_dispersion_s2 * omega) * omega,
                            4e-10,
                        );
                        assert_ne!(phase.group_delay_s, 0.0);
                    }
                }
            }
        }
    }
}

#[test]
fn source_owned_exact_fresnel_zero_keeps_raw_derivatives_without_inventing_phase() {
    let outside = incident(1.5, MU);
    let omega = frequency(600.0);
    let epsilon = DielectricResponse::Lossless(&outside)
        .lossless_permittivity(omega)
        .unwrap();
    let ni = (C0 * fs_math::det::sqrt(MU)) * fs_math::det::sqrt(epsilon);
    let (ordinary, extraordinary) = axes(
        [ni - 0.25, ni + 0.25],
        [ni - 0.25, ni + 0.25],
        [0.0; 2],
        Some(MU),
        [400.0, 800.0],
        0,
    );
    let crystal = UniaxialOpticalResponse::try_new(&ordinary, &extraordinary, SPAN).unwrap();
    assert_eq!(crystal.phase_indices(omega).unwrap(), [ni; 2]);
    let interface = PrincipalPlaneUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside),
        crystal,
        [0.6, 0.8],
    )
    .unwrap();
    for polarization in [Polarization::Te, Polarization::Tm] {
        let normal = interface
            .frequency_dispersion(
                omega,
                1.0,
                polarization,
                IncidenceFrequencyConstraint::FixedTransverseWavenumber,
            )
            .unwrap();
        let result = normal.boundary_frequency_dispersion().unwrap();
        let reflected = result.reflected_tangential_e();
        assert_eq!(result.carrier(), normal.carrier());
        assert_eq!(reflected.amplitude(), C64::ZERO);
        assert_ne!(reflected.scaled_derivatives()[0], C64::ZERO);
        assert!(reflected.logarithmic_dispersion().is_none());
        assert_eq!(
            result
                .reflected_power_fraction()
                .unwrap()
                .scaled_derivative(),
            0.0
        );
        assert_eq!(
            result.transmitted_tangential_e().scaled_derivatives(),
            reflected.scaled_derivatives(),
        );
        assert_eq!(
            result.transmitted_conjugate_h().scaled_derivatives(),
            reflected.scaled_derivatives().map(core::ops::Neg::neg),
        );
    }
}

#[test]
fn fixed_kx_boundary_requires_original_incident_jet_without_narrowing_normal_propagation() {
    // The incident carrier is supported at 500nm, but its original table
    // supplies no two-sided frequency neighbourhood at that endpoint.
    let (_, incident_axis) = axes([1.2; 2], [1.2; 2], [0.0; 2], Some(MU), [500.0, 700.0], 0);
    let incident_law =
        SpectralOpticalResponse::try_new(SpectralOpticalSource::ComplexIndex(&incident_axis), SPAN)
            .unwrap();
    let (ordinary, extraordinary) = axes([1.6; 2], [1.9; 2], [0.0; 2], Some(MU), [400.0, 800.0], 0);
    let crystal = UniaxialOpticalResponse::try_new(&ordinary, &extraordinary, SPAN).unwrap();
    let interface = PrincipalPlaneUniaxialInterface::try_new(
        DielectricResponse::SpectralOptical(incident_law),
        crystal,
        [0.8, 0.6],
    )
    .unwrap();
    let omega = frequency(500.0);
    for polarization in [Polarization::Te, Polarization::Tm] {
        let normal = interface
            .frequency_dispersion_at_tangential_index(
                omega,
                -0.4,
                polarization,
                IncidenceFrequencyConstraint::FixedTransverseWavenumber,
            )
            .unwrap();
        assert!(normal.transmitted_propagation(1e-6).is_ok());
        assert!(normal.boundary_frequency_dispersion().is_err());
        assert!(
            interface
                .boundary_frequency_dispersion_at_tangential_index(
                    omega,
                    -0.4,
                    polarization,
                    IncidenceFrequencyConstraint::FixedTransverseWavenumber,
                )
                .is_err()
        );
        // The same original sources have an admitted interior frequency.
        let interior = interface
            .boundary_frequency_dispersion_at_tangential_index(
                frequency(600.0),
                -0.4,
                polarization,
                IncidenceFrequencyConstraint::FixedTransverseWavenumber,
            )
            .unwrap();
        assert!(core::ptr::eq(
            interior.source().incident_material(),
            incident_axis.resolved()
        ));
    }
}
