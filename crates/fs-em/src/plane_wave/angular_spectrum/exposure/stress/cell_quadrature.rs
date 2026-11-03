//! G1/G3/G4: continuous source-bound stress, including unpaired Nyquist bins.

use super::*;
use crate::plane_wave::DielectricResponse;
use fs_evidence::ValidityDomain;
use fs_fft::harmonic_window::HarmonicWindowProfile;
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
        MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
        ScalarPropertyRequirement, resolve_material_state_point,
    },
};

fn material() -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, dims) in [
        (ELECTRIC_PERMITTIVITY_PROPERTY, ELECTRIC_PERMITTIVITY_DIMS),
        (MAGNETIC_PERMEABILITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value: 1.0, dims },
                validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic modal-stress constitutive fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
                .unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic dielectric".into(),
            phase: "solid".into(),
            process: "analytic modal stress fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    resolve_material_state_point(
        &card,
        &QueryPoint::new().with("T", 300.0).unwrap(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn exposure_limits() -> AngularSpectrumExposureLimits {
    AngularSpectrumExposureLimits {
        maximum_lines: 8,
        maximum_point_line_pairs: 64,
        maximum_point_frequency_pairs: 64,
        maximum_workspace_bytes: 2 << 20,
        maximum_fft_work: 1 << 20,
    }
}

fn ensemble_limits() -> AngularSpectrumExposureEnsembleLimits {
    AngularSpectrumExposureEnsembleLimits {
        maximum_factors: 4,
        maximum_point_line_pairs: 256,
        maximum_point_frequency_pairs: 256,
        maximum_workspace_bytes: 4 << 20,
        maximum_fft_work: 1 << 20,
    }
}

fn prepare<'s, 'm>(
    source: &'s PlanarDielectricLayers<'m>,
    grid: AngularSpectrumGrid,
    samples: &[[C64; 2]],
) -> PlanarAngularSpectrum<'s, 'm> {
    PlanarAngularSpectrum::prepare(
        source,
        grid,
        samples,
        AngularSpectrumLimits {
            maximum_grid_points: 64,
            maximum_workspace_bytes: 1 << 20,
            maximum_fft_work: 1 << 20,
        },
        || false,
    )
    .unwrap()
}

fn close(actual: [[f64; 3]; 3], expected: [[f64; 3]; 3]) {
    for i in 0..3 {
        for j in 0..3 {
            assert!(
                (actual[i][j] - expected[i][j]).abs() < 3e-12 * expected[i][j].abs().max(1.0),
                "stress[{i},{j}]: {} != {}",
                actual[i][j],
                expected[i][j]
            );
        }
    }
}

#[test]
fn g1_source_nyquist_stress_avoids_sample_alias_for_mean_real_and_proper_banks() {
    let medium = material();
    let omega: f64 = 5.0;
    let stack =
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), omega, &[])
            .unwrap();
    let distance = 0.31;
    let transverse = -core::f64::consts::PI;
    let kz = (omega * omega - transverse * transverse).sqrt();
    for along_y in [false, true] {
        let grid = AngularSpectrumGrid {
            width: if along_y { 1 } else { 2 },
            height: if along_y { 2 } else { 1 },
            spacing_m: [1.0; 2],
        };
        let samples: Vec<_> = [C64::ONE, -C64::ONE]
            .into_iter()
            .map(|a| if along_y { [a, C64::ZERO] } else { [C64::ZERO, a] })
            .collect();
        let lines = [prepare(&stack, grid, &samples)];
        let mean =
            AngularSpectrumExposure::try_new(&lines[0], &lines, exposure_limits(), || false)
                .unwrap();
        let factors = [
            AngularSpectrumExposure::try_new(&lines[0], &lines, exposure_limits(), || false)
                .unwrap(),
            AngularSpectrumExposure::try_new(&lines[0], &lines, exposure_limits(), || false)
                .unwrap(),
        ];
        let statistics = [
            AngularSpectrumFactorStatistics::Real,
            AngularSpectrumFactorStatistics::ProperComplex,
        ];
        let mixed = AngularSpectrumExposureEnsemble::try_new_with_statistics(
            &lines[0],
            Some(&mean),
            &factors,
            &statistics,
            ensemble_limits(),
            || false,
        )
        .unwrap();
        let proper = AngularSpectrumExposureEnsemble::try_new(
            &lines[0],
            None,
            &factors[..1],
            ensemble_limits(),
            || false,
        )
        .unwrap();
        let mut direction = [0.0, 0.0, kz / omega];
        direction[usize::from(along_y)] = transverse / omega;
        // E has unit peak amplitude and epsilon=mu=1. At every time the
        // continuous transverse mean of cos²(k.x-omega*t) is exactly 1/2.
        // Thus integral(T_ij) = -area*n_i*n_j/2, with area=2.
        let expected: [[f64; 3]; 3] =
            std::array::from_fn(|i| std::array::from_fn(|j| -direction[i] * direction[j]));
        for profile in [
            HarmonicWindowProfile::Uniform,
            HarmonicWindowProfile::Triangular,
        ] {
            let window = HarmonicWindow::with_profile(0.17, 0.12, profile).unwrap();
            let query_limits = AngularSpectrumExposureStressLimits {
                maximum_workspace_bytes: mixed.stress_workspace_bytes().unwrap(),
            };
            let actual = mean
                .stress_field(
                    AngularSpectrumSide::Transmitted,
                    distance,
                    window,
                    query_limits,
                    || false,
                )
                .unwrap();
            close(actual.integrated_stress_n(), expected);
            let sampled = actual
                .points()
                .iter()
                .map(|p| p.maxwell_stress_pa[2][2])
                .sum::<f64>();
            // Squaring the negative Nyquist samples aliases their carrier
            // term to DC: this sampled result is almost twice the true one.
            assert!(sampled.abs() > 1.9 * expected[2][2].abs());
            if matches!(profile, HarmonicWindowProfile::Uniform) {
                let a = omega * window.start_s() - kz * distance;
                let b = a + omega * window.duration_s();
                let mean_square = 0.5 + ((2.0 * b).sin() - (2.0 * a).sin())
                    / (4.0 * omega * window.duration_s());
                for point in actual.points() {
                    close(
                        point.maxwell_stress_pa,
                        expected.map(|row| row.map(|v| v * mean_square)),
                    );
                }
            }
            let proper_field = proper
                .stress_field(
                    AngularSpectrumSide::Transmitted,
                    distance,
                    window,
                    query_limits,
                    || false,
                )
                .unwrap();
            close(proper_field.integrated_stress_n(), expected);
            let mixed_field = mixed
                .stress_field(
                    AngularSpectrumSide::Transmitted,
                    distance,
                    window,
                    query_limits,
                    || false,
                )
                .unwrap();
            close(
                mixed_field.integrated_stress_n(),
                expected.map(|row| row.map(|v| 3.0 * v)),
            );
        }
    }
}

#[test]
fn g3_true_opposite_modes_keep_real_interference_and_authored_properness() {
    let medium = material();
    let omega: f64 = 5.0;
    let stack =
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), omega, &[])
            .unwrap();
    let grid = AngularSpectrumGrid {
        width: 4,
        height: 1,
        spacing_m: [1.0; 2],
    };
    // Ey=2i*cos(kx*x). Maxwell gives Hx=-2i*(kz/omega)*cos(kx*x)
    // and Hz=2*(kx/omega)*sin(kx*x), so their real carrier phases differ.
    let lines = [prepare(
        &stack,
        grid,
        &[
            [C64::ZERO, C64::new(0.0, 2.0)],
            [C64::ZERO, C64::ZERO],
            [C64::ZERO, C64::new(0.0, -2.0)],
            [C64::ZERO, C64::ZERO],
        ],
    )];
    let real =
        AngularSpectrumExposure::try_new(&lines[0], &lines, exposure_limits(), || false).unwrap();
    let factors = [
        AngularSpectrumExposure::try_new(&lines[0], &lines, exposure_limits(), || false).unwrap(),
    ];
    let proper = AngularSpectrumExposureEnsemble::try_new(
        &lines[0],
        None,
        &factors,
        ensemble_limits(),
        || false,
    )
    .unwrap();
    let window = HarmonicWindow::new(0.17, 0.12).unwrap();
    let phase0 = omega * window.start_s();
    let phase1 = phase0 + omega * window.duration_s();
    let sin2 = 0.5 - ((2.0 * phase1).sin() - (2.0 * phase0).sin())
        / (4.0 * omega * window.duration_s());
    let transverse2 = (core::f64::consts::FRAC_PI_2 / omega).powi(2);
    let limits = AngularSpectrumExposureStressLimits {
        maximum_workspace_bytes: proper.stress_workspace_bytes().unwrap(),
    };
    for (field, square) in [
        (
            real.stress_field(AngularSpectrumSide::Incident, 0.0, window, limits, || false)
                .unwrap(),
            sin2,
        ),
        (
            proper
                .stress_field(AngularSpectrumSide::Incident, 0.0, window, limits, || false)
                .unwrap(),
            0.5,
        ),
    ] {
        close(
            field.integrated_stress_n(),
            [
                [-4.0 * transverse2, 0.0, 0.0],
                [0.0, 4.0 * transverse2 * (2.0 * square - 1.0), 0.0],
                [0.0, 0.0, 4.0 * (transverse2 - 2.0 * square)],
            ],
        );
    }
}

#[test]
fn g4_modal_pass_preserves_work_admission_workspace_cancellation_and_retry() {
    let medium = material();
    let stack =
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), 5.0, &[])
            .unwrap();
    let grid = AngularSpectrumGrid {
        width: 2,
        height: 1,
        spacing_m: [1.0; 2],
    };
    let lines = [prepare(
        &stack,
        grid,
        &[[C64::ZERO, C64::ONE], [C64::ZERO, -C64::ONE]],
    )];
    let bank =
        AngularSpectrumExposure::try_new(&lines[0], &lines, exposure_limits(), || false).unwrap();
    assert_eq!(bank.point_frequency_pairs(), 2);
    assert!(
        AngularSpectrumExposure::try_new(
            &lines[0],
            &lines,
            AngularSpectrumExposureLimits {
                maximum_point_frequency_pairs: 1,
                ..exposure_limits()
            },
            || false,
        )
        .is_err()
    );
    let window = HarmonicWindow::new(0.17, 0.12).unwrap();
    let exact = AngularSpectrumExposureStressLimits {
        maximum_workspace_bytes: bank.stress_workspace_bytes().unwrap(),
    };
    assert!(stress_workspace(usize::MAX, 0).is_err());
    assert!(
        bank.stress_field(
            AngularSpectrumSide::Incident,
            0.0,
            window,
            AngularSpectrumExposureStressLimits {
                maximum_workspace_bytes: exact.maximum_workspace_bytes - 1,
            },
            || false,
        )
        .is_err()
    );
    let modal = bank
        .frequency_modal_field(
            0,
            ExposureSide::OneSided(AngularSpectrumSide::Incident),
            0.0,
            &mut || false,
        )
        .unwrap();
    let fields = [modal];
    let kernels = [window.product_kernel(5.0, 5.0).unwrap()];
    let mut polls = 0;
    assert!(matches!(
        integrate_cell_dyad(
            grid,
            &fields,
            &kernels,
            1.0,
            1.0,
            HarmonicProduct::mean_real_product_parts_dd,
            &mut || {
                polls += 1;
                polls == 3
            },
        ),
        Err(AngularSpectrumError::Cancelled)
    ));
    let first = bank
        .stress_field(AngularSpectrumSide::Incident, 0.0, window, exact, || false)
        .unwrap();
    assert!(matches!(
        bank.stress_field(AngularSpectrumSide::Incident, 0.0, window, exact, || true),
        Err(AngularSpectrumError::Cancelled)
    ));
    let retry = bank
        .stress_field(AngularSpectrumSide::Incident, 0.0, window, exact, || false)
        .unwrap();
    assert_eq!(first.points(), retry.points());
    assert_eq!(first.integrated_stress_n(), retry.integrated_stress_n());
}
