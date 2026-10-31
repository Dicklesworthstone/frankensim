//! Independent continuous/analytic authority and separate live spatial axes.
use core::f64::consts::TAU;
use fs_duct::modal::time_domain::{ModalEnd, ModalTimeError, assemble_lossless};
use fs_duct::modal::{ModeCount, MouthLoad, mm_profile_impedance, modal_reflection_from_impedance};
use fs_duct::{DuctError, LossModel, RadiusProfile, RiccatiSteps, Termination};
use fs_material::gas::{GasSpec, GasState};
use fs_math::{c64::C64, det};
use fs_phs::harmonic::{HarmonicBudget, HarmonicModel};

fn air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
}

fn harmonic_budget() -> HarmonicBudget {
    HarmonicBudget {
        max_dimension: 124,
        max_dense_bytes: 16 * 1024 * 1024,
        max_dense_work: 4_000_000,
        relative_residual_tolerance: 2.0e-10,
    }
}

fn reflection(z: &[C64], reference: f64, n: usize) -> Vec<C64> {
    modal_reflection_from_impedance(z, &vec![C64::from_re(reference); n], n).unwrap()
}

fn distance(a: &[C64], b: &[C64]) -> f64 {
    det::sqrt(a.iter().zip(b).map(|(&x, &y)| (x - y).abs().powi(2)).sum())
}

// Exercise the actual canonical G/Q/J system through the independent shared
// harmonic owner, not an algebraic duplicate of the assembler's formulas.
fn spatial_reflection(
    profile: &RadiusProfile,
    gas: &GasState,
    n: usize,
    elements: usize,
    end: ModalEnd,
    omega: f64,
) -> Vec<C64> {
    let image =
        assemble_lossless(profile, gas, ModeCount(n), elements, end, 124, || false).unwrap();
    let system = image.system();
    assert_eq!(system.state_dim(), 2 * image.source().coordinates());
    assert_eq!(system.port_dim(), n);
    assert!(system.structure().1.iter().all(|&r| r == 0.0));
    let harmonic = HarmonicModel::from_ode(system, harmonic_budget(), || false).unwrap();
    let factor = harmonic.factor(omega, || false).unwrap();
    let mut impedance = vec![C64::ZERO; n * n];
    for column in 0..n {
        let mut input = vec![C64::ZERO; n];
        input[column] = C64::ONE;
        let response = factor.response(&input, || false).unwrap();
        for row in 0..n {
            impedance[row * n + column] = response.output[row];
        }
        assert_eq!(response.dissipated_power_w, 0.0);
    }
    for row in 0..n {
        for column in 0..n {
            assert!(
                (impedance[row * n + column] - impedance[column * n + row]).abs()
                    < 1.0e-8 * image.source().inlet_reference_impedance()
            );
        }
    }
    reflection(&impedance, image.source().inlet_reference_impedance(), n)
}

fn authority(
    profile: &RadiusProfile,
    gas: &GasState,
    n: usize,
    end: ModalEnd,
    omega: f64,
) -> Vec<C64> {
    let termination = match end {
        ModalEnd::Rigid => Termination::Closed,
        ModalEnd::PressureRelease => Termination::IdealOpen,
    };
    let exact = mm_profile_impedance(
        profile,
        gas,
        omega,
        LossModel::Lossless,
        MouthLoad::Analytic(termination),
        ModeCount(n),
        RiccatiSteps(512),
    )
    .unwrap();
    let z0 =
        gas.density * gas.sound_speed / (core::f64::consts::PI * profile.inlet_radius().powi(2));
    reflection(&exact.impedance_matrix, z0, n)
}

#[test]
fn uniform_all_mode_ideal_ends_follow_propagating_and_evanescent_formulas() {
    let gas = air();
    let (radius, length) = (0.05, 0.04);
    let profile = RadiusProfile::exponential(radius, 0.0, length).unwrap();
    for end in [ModalEnd::Rigid, ModalEnd::PressureRelease] {
        let image =
            assemble_lossless(&profile, &gas, ModeCount(3), 16, end, 124, || false).unwrap();
        assert_eq!(
            image.source().nodes_m().len(),
            16 + usize::from(end == ModalEnd::Rigid)
        );
        let z0 = image.source().inlet_reference_impedance();
        for hz in [700.0, 5000.0] {
            let k0 = TAU * hz / gas.sound_speed;
            let mut exact = vec![C64::ZERO; 9];
            for (mode, &gamma) in image.source().radial_eigenvalues().iter().enumerate() {
                let squared = k0 * k0 - (gamma / radius).powi(2);
                let imaginary = if squared > 0.0 {
                    let k = det::sqrt(squared);
                    let tangent = det::tan(k * length);
                    z0 * k0 / k
                        * match end {
                            ModalEnd::Rigid => 1.0 / tangent,
                            ModalEnd::PressureRelease => -tangent,
                        }
                } else {
                    let alpha = det::sqrt(-squared);
                    let e = det::exp(-2.0 * alpha * length);
                    -z0 * k0 / alpha
                        * match end {
                            ModalEnd::Rigid => (1.0 + e) / (1.0 - e),
                            ModalEnd::PressureRelease => (1.0 - e) / (1.0 + e),
                        }
                };
                exact[mode * 3 + mode] = C64::new(0.0, imaginary);
            }
            let actual = spatial_reflection(&profile, &gas, 3, 16, end, TAU * hz);
            let expected = reflection(&exact, z0, 3);
            assert!(distance(&actual, &expected) < 0.025, "end={end:?}, hz={hz}");
            // No accidental coupling survives a truly constant radius.
            for i in 0..3 {
                for j in 0..3 {
                    if i != j {
                        assert!(actual[i * 3 + j].abs() < 1.0e-10);
                    }
                }
            }
        }
    }
}

#[test]
fn axial_ladder_converges_every_matrix_entry_to_original_continuous_riccati() {
    let gas = air();
    let profile = RadiusProfile::exponential(0.025, 4.0, 0.12).unwrap();
    for end in [ModalEnd::Rigid, ModalEnd::PressureRelease] {
        let expected = authority(&profile, &gas, 2, end, TAU * 600.0);
        let errors: Vec<_> = [4, 8, 16]
            .into_iter()
            .map(|elements| {
                distance(
                    &spatial_reflection(&profile, &gas, 2, elements, end, TAU * 600.0),
                    &expected,
                )
            })
            .collect();
        assert!(
            errors[0] > 0.03,
            "a live spatial defect is required: {errors:?}"
        );
        assert!(
            errors[1] < 0.6 * errors[0] && errors[2] < 0.4 * errors[1],
            "{errors:?}"
        );
        assert!(errors[2] < 0.01, "{errors:?}");
        assert!(
            expected[1].abs() > 1.0e-3,
            "the continuous cross-mode reflection must be live"
        );
    }
}

#[test]
fn separate_radial_ladder_changes_and_converges_plane_reflection() {
    let gas = air();
    let profile = RadiusProfile::exponential(0.025, 10.0, 0.12).unwrap();
    let expected = authority(&profile, &gas, 4, ModalEnd::Rigid, TAU * 600.0)[0];
    let values: Vec<_> = [1, 2, 4]
        .into_iter()
        .map(|n| spatial_reflection(&profile, &gas, n, 14, ModalEnd::Rigid, TAU * 600.0)[0])
        .collect();
    let errors: Vec<_> = values.iter().map(|&v| (v - expected).abs()).collect();
    assert!(
        errors[0] > 0.008,
        "a live omitted-mode defect is required: {errors:?}"
    );
    assert!(
        errors[1] < 0.6 * errors[0] && errors[2] < 0.6 * errors[1],
        "{errors:?}"
    );
    assert!(errors[2] < 0.003, "{errors:?}");
    assert!((values[1] - values[0]).abs() > 0.004);
    assert!((values[2] - values[1]).abs() > 0.001);
}

#[test]
fn source_preserves_all_coordinates_and_refuses_capacity_or_cancellation() {
    let gas = air();
    let profile = RadiusProfile::exponential(0.025, 4.0, 0.12).unwrap();
    // All 136 states exceed the shared full exact-flow capacity; never keep
    // only the low eigenmodes, even if the caller supplies a larger ceiling.
    assert!(matches!(
        assemble_lossless(
            &profile,
            &gas,
            ModeCount(4),
            16,
            ModalEnd::Rigid,
            1000,
            || false
        ),
        Err(ModalTimeError::Invalid(_))
    ));
    assert!(matches!(
        assemble_lossless(&profile, &gas, ModeCount(2), 4, ModalEnd::Rigid, 19, || {
            false
        }),
        Err(ModalTimeError::Invalid(_))
    ));
    assert!(matches!(
        assemble_lossless(
            &profile,
            &gas,
            ModeCount(0),
            4,
            ModalEnd::Rigid,
            124,
            || false
        ),
        Err(ModalTimeError::Duct(DuctError::BadParameter { .. }))
    ));
    let mut calls = 0;
    assert!(matches!(
        assemble_lossless(
            &profile,
            &gas,
            ModeCount(2),
            8,
            ModalEnd::Rigid,
            124,
            || {
                calls += 1;
                calls == 3
            }
        ),
        Err(ModalTimeError::Duct(DuctError::Cancelled {
            completed_steps: 1
        }))
    ));
}
