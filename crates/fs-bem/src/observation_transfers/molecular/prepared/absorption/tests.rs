use super::{MolecularAbsorptionLimits, MolecularAbsorptionReading};
use crate::{
    helmholtz::laplace::complex::ComplexOptions,
    observation_transfers::{
        HarmonicField, HarmonicVolumeQuadrature, MolecularBemObservation, MolecularResponseLimits,
        PreparedMolecularResponse, ReceiverSpec, SourceId,
    },
    panel3d::SpherePanels,
};
use fs_material::{
    gas::{GasSpec, GasState},
    relaxation::{RelaxationState, RelaxingMedium, input::read_relaxation},
};
use fs_math::c64::C64;

const LIMITS: MolecularAbsorptionLimits = MolecularAbsorptionLimits {
    maximum_quadrature_points: 16,
    maximum_coherence_modes: 4,
    maximum_products: 1_000_000,
    maximum_workspace_bytes: 1_000_000,
};

fn material() -> RelaxationState {
    // The original prepared-response fixture's declared positive reservoir,
    // not an experimental atmospheric absorption fit.
    let text = "declared 1 synthetic 1 vibration 1 900 1 power 101.325 293.15 0 200 500";
    let law = read_relaxation(&mut text.split_ascii_whitespace()).unwrap();
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap();
    let medium = RelaxingMedium::new(gas, vec![(law.species()[0].species.clone(), 1.0)]).unwrap();
    law.evaluate(&medium).unwrap()
}

fn response(duplicate_traces: bool) -> PreparedMolecularResponse {
    let surface = SpherePanels::icosphere(0.1, 0).unwrap();
    let panels = surface.areas().len();
    let ids = [
        SourceId::StructuralMode("radial".into()),
        SourceId::StructuralMode("second drive".into()),
    ];
    let traces = [
        vec![1.0; panels],
        if duplicate_traces {
            vec![1.0; panels]
        } else {
            surface.normals().iter().map(|normal| normal[2]).collect()
        },
    ];
    let receivers = [
        ReceiverSpec { id: "axial".into(), position_m: [0.0, 0.0, 1.0] },
        ReceiverSpec { id: "oblique".into(), position_m: [0.4, -0.2, 0.8] },
    ];
    let observer = MolecularBemObservation::new(
        &surface,
        material(),
        &ids,
        &traces,
        &receivers,
        0.01,
        [1.0, 3_000.0],
        ComplexOptions {
            max_panels: panels,
            max_bytes: 32 * panels * panels,
            max_work: 512 * (panels * panels) as u64,
        },
    )
    .unwrap();
    observer
        .prepare_response(
            600.0,
            MolecularResponseLimits {
                maximum_retained_entries: 3 * 2 * 2 + 2 * 2 * panels + 8 * 2 * 2,
                maximum_panel_products: (3 * 2 * 2 * panels) as u64,
            },
            || false,
        )
        .unwrap()
}

fn quadrature() -> [HarmonicVolumeQuadrature; 3] {
    // Authored physical volume weights; these tests qualify the discrete
    // observation and units, not spatial cubature convergence.
    [
        HarmonicVolumeQuadrature { receiver: 0, volume_m3: 0.125 },
        HarmonicVolumeQuadrature { receiver: 1, volume_m3: 0.5 },
        HarmonicVolumeQuadrature { receiver: 0, volume_m3: 0.375 },
    ]
}

fn values(reading: MolecularAbsorptionReading<'_>) -> (Vec<f64>, f64) {
    (reading.absorbed_power_density_w_m3.to_vec(), reading.absorbed_power_w)
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs()
            <= 5e-12 * actual.abs().max(expected.abs()).max(f64::MIN_POSITIVE),
        "actual {actual:e}, expected {expected:e}"
    );
}

#[test]
fn g1_original_caloric_heat_matches_real_pressure_cycles_and_physical_volume() {
    let response = response(false);
    let q = quadrature();
    let rates = [C64::new(0.5, -0.2), C64::new(-0.3, 0.4)];
    let mut field = [HarmonicField::default(); 2];
    response.harmonic().evaluate_coherent(&rates, &mut field, || false).unwrap();
    let mut absorption = response.prepare_absorption(&q, LIMITS, || false).unwrap();
    assert!(core::ptr::eq(absorption.source(), response.source()));
    assert!(core::ptr::eq(absorption.response(), &response));
    assert_eq!(absorption.sources(), response.sources());
    assert_eq!(absorption.receivers(), response.harmonic().receivers());
    assert_eq!(absorption.producer_hash(), response.harmonic().producer_hash());
    assert_eq!(absorption.angular_frequency_rad_s(), 600.0);
    assert_eq!(absorption.volume_m3(), 1.0);

    let omega = response.harmonic().angular_frequency_rad_s();
    let compliance = response.source().caloric_response(omega).unwrap().bulk_modulus_pa.recip();
    assert!(compliance.im > 0.0);
    let mut expected = [0.0; 2];
    for (r, value) in field.iter().enumerate() {
        // Independently time-average the REAL pressure, avoiding a duplicated
        // peak-phasor factor of two in the new observation's expected value.
        let mut mean_square = 0.0;
        for sample in 0..256 {
            let phase = core::f64::consts::TAU * (sample as f64 + 0.25) / 256.0;
            let pressure = value.pressure_pa.re * phase.cos() + value.pressure_pa.im * phase.sin();
            mean_square += pressure * pressure / 256.0;
        }
        expected[r] = omega * compliance.im * mean_square;
    }
    let (density, watts) = values(absorption.evaluate_coherent(&rates, || false).unwrap());
    for r in 0..2 {
        close(density[r], expected[r]);
        assert!(density[r] > 0.0);
    }
    let integrated: f64 = q.iter().map(|sample| sample.volume_m3 * expected[sample.receiver]).sum();
    close(watts, integrated);
}

#[test]
fn g3_coherent_cancellation_and_independent_sources_have_distinct_physical_heat() {
    let response = response(true);
    let q = quadrature();
    let mut absorption = response.prepare_absorption(&q, LIMITS, || false).unwrap();
    let a = C64::new(0.6, -0.25);
    let (_, solo) = values(absorption.evaluate_coherent(&[a, C64::ZERO], || false).unwrap());
    assert!(solo > 0.0);
    let (_, cancelling) = values(absorption.evaluate_coherent(&[a, -a], || false).unwrap());
    let (_, reinforcing) = values(absorption.evaluate_coherent(&[a, a], || false).unwrap());
    let (independent_density, independent) = values(
        absorption.evaluate_independent_sources(&[a.norm_sq(), a.norm_sq()], || false).unwrap(),
    );
    assert!(cancelling <= 1e-24 * solo);
    close(reinforcing, 4.0 * solo);
    close(independent, 2.0 * solo);

    let factors = [a, C64::ZERO, C64::ZERO, a];
    let (factor_density, factor_power) = values(
        absorption.evaluate_coherence_modes(&factors, 2, || false).unwrap(),
    );
    close(factor_power, independent);
    for (factor, diagonal) in factor_density.iter().zip(&independent_density) {
        close(*factor, *diagonal);
    }
    let (_, correlated) = values(
        absorption.evaluate_coherence_modes(&[a, -a], 1, || false).unwrap(),
    );
    close(correlated, cancelling);
}

#[test]
fn g3_factor_rotations_and_volume_partitions_preserve_their_original_semantics() {
    let response = response(false);
    let q = quadrature();
    let merged = [
        HarmonicVolumeQuadrature { receiver: 1, volume_m3: 0.5 },
        HarmonicVolumeQuadrature { receiver: 0, volume_m3: 0.5 },
    ];
    let enlarged = [
        HarmonicVolumeQuadrature { receiver: 0, volume_m3: 1.5 },
        HarmonicVolumeQuadrature { receiver: 1, volume_m3: 1.5 },
    ];
    let mut split = response.prepare_absorption(&q, LIMITS, || false).unwrap();
    let mut merged = response.prepare_absorption(&merged, LIMITS, || false).unwrap();
    let mut enlarged = response.prepare_absorption(&enlarged, LIMITS, || false).unwrap();
    let factors = [
        C64::new(0.3, -0.2), C64::new(-0.1, 0.4),
        C64::new(0.5, 0.1), C64::new(0.2, -0.3),
    ];
    let h = core::f64::consts::FRAC_1_SQRT_2;
    let rotated = [
        (factors[0] + factors[2]).scale(h), (factors[1] + factors[3]).scale(h),
        (factors[0] - factors[2]).scale(h), (factors[1] - factors[3]).scale(h),
    ];
    // Unitary mixing changes the factorization but keeps B B* unchanged.
    let (density, power) = values(split.evaluate_coherence_modes(&factors, 2, || false).unwrap());
    let (rotated_density, rotated_power) =
        values(split.evaluate_coherence_modes(&rotated, 2, || false).unwrap());
    close(rotated_power, power);
    for (a, b) in rotated_density.iter().zip(&density) {
        close(*a, *b);
    }
    let (merged_density, merged_power) =
        values(merged.evaluate_coherence_modes(&factors, 2, || false).unwrap());
    assert_eq!(merged_density, density);
    close(merged_power, power);
    let (enlarged_density, enlarged_power) =
        values(enlarged.evaluate_coherence_modes(&factors, 2, || false).unwrap());
    assert_eq!(enlarged_density, density);
    assert_eq!(enlarged.volume_m3(), 3.0);
    close(enlarged_power, 3.0 * power);
}

#[test]
fn g0_volume_work_storage_and_rank_limits_refuse_before_successful_observation() {
    let response = response(false);
    let q = quadrature();
    for bad in [
        HarmonicVolumeQuadrature { receiver: 2, volume_m3: 1.0 },
        HarmonicVolumeQuadrature { receiver: usize::MAX, volume_m3: 0.0 },
        HarmonicVolumeQuadrature { receiver: 0, volume_m3: -1.0 },
        HarmonicVolumeQuadrature { receiver: 0, volume_m3: f64::NAN },
        HarmonicVolumeQuadrature { receiver: 0, volume_m3: f64::INFINITY },
    ] {
        assert!(response.prepare_absorption(&[bad], LIMITS, || false).is_err());
    }
    let workspace = response.prepare_absorption(&q, LIMITS, || false).unwrap().workspace_bytes();
    assert!(workspace > 0);
    assert!(response.prepare_absorption(
        &q, MolecularAbsorptionLimits { maximum_workspace_bytes: workspace - 1, ..LIMITS },
        || false,
    ).is_err());
    assert!(response.prepare_absorption(
        &q, MolecularAbsorptionLimits { maximum_quadrature_points: q.len() - 1, ..LIMITS },
        || false,
    ).is_err());
    let rates = [C64::ONE, C64::ONE];
    // A fixed allowance may be rejected at preparation or on a query whose
    // actual work exceeds it; neither route may successfully compute this field.
    if let Ok(mut limited) = response.prepare_absorption(
        &q, MolecularAbsorptionLimits { maximum_products: 1, ..LIMITS }, || false,
    ) {
        assert!(limited.evaluate_coherent(&rates, || false).is_err());
    }
    let mut single_mode = response.prepare_absorption(
        &q, MolecularAbsorptionLimits { maximum_coherence_modes: 1, ..LIMITS }, || false,
    ).unwrap();
    assert!(single_mode.evaluate_coherence_modes(&[C64::ONE; 4], 2, || false).is_err());
    assert!(single_mode.evaluate_coherent(&rates, || false).is_ok());
}

#[test]
fn g4_refusal_and_cancellation_leave_preparation_reusable_and_silence_exact() {
    let response = response(false);
    let q = quadrature();
    assert!(response.prepare_absorption(&q, LIMITS, || true).is_err());
    let mut absorption = response.prepare_absorption(&q, LIMITS, || false).unwrap();
    let rates = [C64::new(0.2, -0.1), C64::new(0.4, 0.3)];
    let baseline = values(absorption.evaluate_coherent(&rates, || false).unwrap());
    assert!(absorption.evaluate_coherent(&rates[..1], || false).is_err());
    assert!(absorption.evaluate_coherent(&[C64::new(f64::NAN, 0.0), C64::ONE], || false).is_err());
    assert!(absorption.evaluate_coherence_modes(&rates, 2, || false).is_err());
    for bad in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(absorption.evaluate_independent_sources(&[1.0, bad], || false).is_err());
    }
    assert!(absorption.evaluate_independent_sources(&[1.0], || false).is_err());
    let mut calls = 0;
    assert!(absorption.evaluate_coherent(&rates, || {
        calls += 1;
        calls == 2
    }).is_err());
    assert_eq!(calls, 2);
    assert!(absorption.evaluate_coherence_modes(&rates, 1, || true).is_err());
    assert_eq!(values(absorption.evaluate_coherent(&rates, || false).unwrap()), baseline);

    let (empty_modes, empty_power) =
        values(absorption.evaluate_coherence_modes(&[], 0, || false).unwrap());
    assert_eq!(empty_modes, vec![0.0; 2]);
    assert_eq!(empty_power, 0.0);
    let (silent, silent_power) =
        values(absorption.evaluate_independent_sources(&[0.0; 2], || false).unwrap());
    assert_eq!(silent, vec![0.0; 2]);
    assert_eq!(silent_power, 0.0);
    assert_eq!(values(absorption.evaluate_coherent(&rates, || false).unwrap()), baseline);
}
