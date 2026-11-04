//! Measured coherence through the original 3D pressure/velocity propagators.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::{
    CoherenceBudget, CoherentCompactField, CompactFieldStatistics, CompactSourcePhasor,
    ConvectedMonopoleLimits, IdealAcousticPlane, IdealPlaneBoundary, PreparedCompactSourceObserver,
    PreparedConvectedMonopoleObserver,
};
use fs_fft::welch::{WelchBudget, WelchConfig, WelchCrossSpectrum, WelchDetrend, WelchWindow};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;
use std::f64::consts::TAU;

const N: usize = 64;
const BIN: usize = 7;
const DT: f64 = 1.0 / 8192.0;
const POSITIONS: [[f64; 3]; 2] = [[0.05, 0.08, 0.12], [0.13, 0.02, 0.21]];
const RECEIVER: [f64; 3] = [0.39, 0.28, 0.17];
const WORK: CoherenceBudget = CoherenceBudget {
    max_modes: 3,
    max_source_mode_pairs: 6,
};

fn medium() -> FluidAcousticState {
    FluidAcousticState::try_new(1.17, 347.0).unwrap()
}

fn record<const D: usize>(segments: &[[C64; D]]) -> Vec<[f64; D]> {
    segments
        .iter()
        .flat_map(|amplitudes| {
            (0..N).map(move |j| {
                let theta = TAU * BIN as f64 * j as f64 / N as f64;
                // Real exp(-i omega t) record, not synthesized estimator factors.
                amplitudes.map(|a| a.re * theta.cos() + a.im * theta.sin())
            })
        })
        .collect()
}

fn estimate<const D: usize>(segments: &[[C64; D]], window: WelchWindow) -> WelchCrossSpectrum {
    WelchCrossSpectrum::estimate_interleaved(
        &record(segments),
        DT,
        WelchConfig {
            segment_samples: N,
            hop_samples: N,
            window,
            detrend: WelchDetrend::None,
        },
        WelchBudget {
            maximum_channels: 8,
            maximum_samples_per_channel: 3 * N,
            maximum_transforms: 24,
            maximum_workspace_bytes: 1024 * 1024,
        },
        || false,
    )
    .unwrap()
}

fn segments() -> [[CompactSourcePhasor; 2]; 3] {
    let a = CompactSourcePhasor {
        volume_flow_m3_s: C64::new(2e-5, -1e-5),
        force_on_fluid_n: [
            C64::new(0.012, -0.006),
            C64::new(-0.003, 0.008),
            C64::new(0.002, 0.005),
        ],
    };
    let b = CompactSourcePhasor {
        volume_flow_m3_s: C64::new(-1e-5, 3e-5),
        force_on_fluid_n: [
            C64::new(0.004, 0.007),
            C64::new(0.009, -0.002),
            C64::new(-0.006, 0.003),
        ],
    };
    [C64::ONE, C64::new(0.0, 1.0), C64::new(-0.5, 0.2)].map(|phase| {
        [
            a,
            CompactSourcePhasor {
                volume_flow_m3_s: b.volume_flow_m3_s * phase,
                force_on_fluid_n: b.force_on_fluid_n.map(|f| f * phase),
            },
        ]
    })
}

fn channels(sources: [CompactSourcePhasor; 2]) -> [C64; 8] {
    core::array::from_fn(|i| {
        let source = sources[i / 4];
        if i % 4 == 0 {
            source.volume_flow_m3_s
        } else {
            source.force_on_fluid_n[i % 4 - 1]
        }
    })
}

fn observer(frequency: f64) -> PreparedCompactSourceObserver {
    PreparedCompactSourceObserver::new(medium(), frequency, RECEIVER, &POSITIONS, 2).unwrap()
}

fn close(actual: C64, expected: C64) {
    assert!(
        (actual - expected).abs() <= 3e-11 * expected.abs().max(1e-13),
        "actual={actual:?}, expected={expected:?}"
    );
}

// Average known physical source phasors after original deterministic propagation.
// No Welch factor or cross-spectrum value enters this reference calculation.
fn check_moments(actual: &CompactFieldStatistics, fields: &[CoherentCompactField], weight: f64) {
    let moment = |i: usize, j: usize| {
        let component = |field: &CoherentCompactField, axis| {
            if axis == 0 {
                field.pressure_pa
            } else {
                field.particle_velocity_m_s[axis - 1]
            }
        };
        fields
            .iter()
            .fold(C64::ZERO, |sum, field| {
                sum + component(field, i) * component(field, j).conj()
            })
            .scale(0.5 * weight / fields.len() as f64)
    };
    close(C64::from_re(actual.mean_square_pressure_pa2), moment(0, 0));
    for axis in 0..3 {
        close(
            C64::new(
                actual.active_intensity_w_m2[axis],
                actual.reactive_intensity_w_m2[axis],
            ),
            moment(0, axis + 1),
        );
        for other in 0..3 {
            close(
                actual.velocity_covariance_m2_s2[axis][other],
                moment(axis + 1, other + 1),
            );
        }
    }
    assert_eq!(actual.scope, fields[0].scope);
}

#[test]
fn measured_mixed_sources_retain_partial_coherence_and_complex_velocity_cross_terms() {
    let inputs = segments();
    let spectrum = estimate(&inputs.map(channels), WelchWindow::Rectangular);
    let prepared = observer(spectrum.frequency_hz(BIN).unwrap());
    let fields: Vec<_> = inputs
        .iter()
        .map(|s| prepared.evaluate(s).unwrap())
        .collect();
    let actual = prepared
        .evaluate_welch_bin(&spectrum, BIN, WORK, || false)
        .unwrap();
    check_moments(&actual, &fields, 1.0);
    assert!(actual.velocity_covariance_m2_s2[0][1].im.abs() > 1e-10);
    assert_eq!(prepared.acoustic(), medium());
    let repeat = prepared
        .evaluate_welch_bin(&spectrum, BIN, WORK, || false)
        .unwrap();
    assert_eq!(
        actual.mean_square_pressure_pa2,
        repeat.mean_square_pressure_pa2
    );
    assert_eq!(
        actual.velocity_covariance_m2_s2,
        repeat.velocity_covariance_m2_s2
    );
}

#[test]
fn periodic_hann_retains_physical_bin_energy_in_center_and_adjacent_bins() {
    let inputs = segments();
    let spectrum = estimate(&inputs.map(channels), WelchWindow::PeriodicHann);
    // A bin-centered real sinusoid contributes 2/3 of its power at the center
    // and 1/6 in each adjacent bin under a periodic Hann modified periodogram.
    // Each bin still uses its OWN physical propagation frequency.
    for (bin, fraction) in [(BIN - 1, 1.0 / 6.0), (BIN, 2.0 / 3.0), (BIN + 1, 1.0 / 6.0)] {
        let prepared = observer(spectrum.frequency_hz(bin).unwrap());
        let fields: Vec<_> = inputs
            .iter()
            .map(|s| prepared.evaluate(s).unwrap())
            .collect();
        let actual = prepared
            .evaluate_welch_bin(&spectrum, bin, WORK, || false)
            .unwrap();
        check_moments(&actual, &fields, fraction);
    }
}

#[test]
fn measured_volume_force_cancellation_keeps_velocity_and_changes_with_correlation() {
    let f = BIN as f64 / (N as f64 * DT);
    let r = 0.31;
    let acoustic = medium();
    let q = C64::from_re(2e-5);
    // From p_Q=-i omega rho G Q and p_F=(1/r-i k) G Fx on the x axis.
    let force = C64::new(0.0, TAU * f * acoustic.density_kg_m3()) * q
        / C64::new(1.0 / r, -TAU * f / acoustic.sound_speed_m_s());
    let destructive = [q, force, C64::ZERO, C64::ZERO];
    let constructive = [q, -force, C64::ZERO, C64::ZERO];
    let prepared =
        PreparedCompactSourceObserver::new(acoustic, f, [r, 0.0, 0.0], &[[0.0; 3]], 1).unwrap();
    let mono = prepared
        .evaluate(&[CompactSourcePhasor::monopole(q)])
        .unwrap();
    let single_power = 0.5 * mono.pressure_pa.norm_sq();
    let coherent = estimate(&[destructive; 2], WelchWindow::Rectangular);
    let quiet = prepared
        .evaluate_welch_bin(&coherent, BIN, WORK, || false)
        .unwrap();
    assert!(quiet.mean_square_pressure_pa2 < single_power * 1e-24);
    assert!(quiet.velocity_covariance_m2_s2[0][0].re > 1e-12);

    // Reversing the force in half the measured segments removes Q/F cross
    // covariance. It must restore the sum of their individual powers.
    let uncorrelated = estimate(&[destructive, constructive], WelchWindow::Rectangular);
    let restored = prepared
        .evaluate_welch_bin(&uncorrelated, BIN, WORK, || false)
        .unwrap();
    close(
        C64::from_re(restored.mean_square_pressure_pa2),
        C64::from_re(2.0 * single_power),
    );
}

#[test]
fn measured_boundary_paths_remain_coherent_at_rigid_and_pressure_release_walls() {
    let inputs = segments();
    let spectrum = estimate(&inputs.map(channels), WelchWindow::Rectangular);
    let frequency = spectrum.frequency_hz(BIN).unwrap();
    let wall_receiver = [0.39, 0.28, 0.0];
    let boundary_work = CoherenceBudget {
        max_source_mode_pairs: 12,
        ..WORK
    };
    for boundary in [
        IdealPlaneBoundary::Rigid,
        IdealPlaneBoundary::PressureRelease,
    ] {
        let plane = IdealAcousticPlane::new([0.0; 3], [0.0, 0.0, 1.0], boundary).unwrap();
        let direct =
            PreparedCompactSourceObserver::new(medium(), frequency, wall_receiver, &POSITIONS, 2)
                .unwrap();
        let half = direct.clone().with_plane(plane, 4).unwrap();
        let folded = direct.with_orthogonal_planes(&[plane], 4).unwrap();
        let actual = half
            .evaluate_welch_bin(&spectrum, BIN, boundary_work, || false)
            .unwrap();
        let fields: Vec<_> = inputs.iter().map(|s| half.evaluate(s).unwrap()).collect();
        check_moments(&actual, &fields, 1.0);
        let folded_actual = folded
            .evaluate_welch_bin(&spectrum, BIN, WORK, || false)
            .unwrap();
        let folded_fields: Vec<_> = inputs.iter().map(|s| folded.evaluate(s).unwrap()).collect();
        check_moments(&folded_actual, &folded_fields, 1.0);
        close(
            C64::from_re(folded_actual.mean_square_pressure_pa2),
            C64::from_re(actual.mean_square_pressure_pa2),
        );
        if boundary == IdealPlaneBoundary::Rigid {
            assert!(actual.mean_square_pressure_pa2 > 1e-3);
            assert_eq!(actual.velocity_covariance_m2_s2[2][2], C64::ZERO);
        } else {
            assert!(actual.mean_square_pressure_pa2 < 1e-28);
            assert!(actual.velocity_covariance_m2_s2[2][2].re > 1e-9);
        }
        assert!(
            half.evaluate_welch_bin(&spectrum, BIN, WORK, || false)
                .is_err()
        );
    }
}

#[test]
fn measured_uniform_flow_preserves_correlated_fields_and_advected_energy() {
    let flows = [
        [C64::new(2e-5, -1e-5), C64::new(-1e-5, 3e-5)],
        [C64::new(2e-5, -1e-5), C64::new(-3e-5, -1e-5)],
        [C64::new(2e-5, -1e-5), C64::new(1e-5, -1.5e-5)],
    ];
    let spectrum = estimate(&flows, WelchWindow::Rectangular);
    for velocity in [[0.0; 3], [61.0, -28.0, 13.0]] {
        let prepared = PreparedConvectedMonopoleObserver::new(
            medium(),
            velocity,
            spectrum.frequency_hz(BIN).unwrap(),
            RECEIVER,
            &POSITIONS,
            2,
            ConvectedMonopoleLimits {
                maximum_mach: 0.8,
                minimum_distance_m: 1e-5,
            },
        )
        .unwrap();
        let fields: Vec<_> = flows
            .iter()
            .map(|q| prepared.evaluate(q).unwrap())
            .collect();
        let actual = prepared
            .evaluate_welch_bin(&spectrum, BIN, WORK, || false)
            .unwrap();
        let reference: Vec<_> = fields.iter().map(|f| f.field.clone()).collect();
        check_moments(&actual.field, &reference, 1.0);
        let energy = fields
            .iter()
            .map(|f| f.mean_perturbation_energy_density_j_m3)
            .sum::<f64>()
            / 3.0;
        close(
            C64::from_re(actual.mean_perturbation_energy_density_j_m3),
            C64::from_re(energy),
        );
        for axis in 0..3 {
            let flux = fields
                .iter()
                .map(|f| f.mean_perturbation_energy_flux_w_m2[axis])
                .sum::<f64>()
                / 3.0;
            close(
                C64::from_re(actual.mean_perturbation_energy_flux_w_m2[axis]),
                C64::from_re(flux),
            );
            close(
                C64::from_re(actual.mean_perturbation_energy_flux_w_m2[axis]),
                C64::from_re(actual.field.active_intensity_w_m2[axis] + velocity[axis] * energy),
            );
        }
    }
}

#[test]
fn recorded_3d_layout_frequency_work_and_cancellation_refuse_without_poisoning_replay() {
    let spectrum = estimate(&segments().map(channels), WelchWindow::Rectangular);
    let prepared = observer(spectrum.frequency_hz(BIN).unwrap());
    for bin in [0, N / 2, BIN + 1, usize::MAX] {
        assert!(matches!(
            prepared.evaluate_welch_bin(&spectrum, bin, WORK, || false),
            Err(AeroacError::InvalidParameter { .. })
        ));
    }
    let wrong_channels = estimate(&[[C64::ONE; 2]; 3], WelchWindow::Rectangular);
    assert!(
        prepared
            .evaluate_welch_bin(&wrong_channels, BIN, WORK, || false)
            .is_err()
    );
    for limited in [
        CoherenceBudget {
            max_modes: 2,
            ..WORK
        },
        CoherenceBudget {
            max_source_mode_pairs: 5,
            ..WORK
        },
    ] {
        assert!(
            prepared
                .evaluate_welch_bin(&spectrum, BIN, limited, || false)
                .is_err()
        );
    }
    let mut polls = 0;
    assert!(matches!(
        prepared.evaluate_welch_bin(&spectrum, BIN, WORK, || {
            polls += 1;
            polls == 5
        }),
        Err(AeroacError::Cancelled { .. })
    ));
    let result = prepared
        .evaluate_welch_bin(&spectrum, BIN, WORK, || false)
        .unwrap();
    let fields: Vec<_> = segments()
        .iter()
        .map(|s| prepared.evaluate(s).unwrap())
        .collect();
    check_moments(&result, &fields, 1.0);

    let flow = PreparedConvectedMonopoleObserver::new(
        medium(),
        [10.0, 0.0, 0.0],
        spectrum.frequency_hz(BIN).unwrap(),
        RECEIVER,
        &POSITIONS,
        2,
        ConvectedMonopoleLimits {
            maximum_mach: 0.5,
            minimum_distance_m: 1e-5,
        },
    )
    .unwrap();
    assert!(
        flow.evaluate_welch_bin(&spectrum, BIN, WORK, || false)
            .is_err()
    );
    for bin in [0, N / 2, BIN + 1] {
        assert!(
            flow.evaluate_welch_bin(&wrong_channels, bin, WORK, || false)
                .is_err()
        );
    }
    assert!(matches!(
        flow.evaluate_welch_bin(&wrong_channels, BIN, WORK, || true),
        Err(AeroacError::Cancelled { .. })
    ));
}

#[test]
fn measured_molecular_propagation_retains_original_relaxation_and_reflected_phase() {
    use fs_material::gas::{GasSpec, GasState};
    use fs_material::relaxation::{
        InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw,
        RelaxingMedium, SpeciesRelaxation,
    };
    use fs_thermochem::SpeciesId;

    // The existing synthetic caloric-reservoir fixture tests the original
    // constitutive law; it does not claim measured parameters for a real gas.
    let gas = GasState::try_new(
        &GasSpec {
            gamma: 1.25,
            ..GasSpec::dry_air_ussa1976()
        },
        293.15,
        101325.0,
    )
    .unwrap();
    let species = SpeciesId::new("synthetic").unwrap();
    let relaxation = RelaxationAbsorption::new(vec![SpeciesRelaxation {
        species: species.clone(),
        processes: vec![RelaxationProcess {
            reservoir: InternalReservoir::Rotational {
                degrees_of_freedom: 2,
            },
            rate: RelaxationRateLaw::PartnerIndependent {
                p_tau_ref_pa_s: gas.pressure * 3e-4,
                t_ref_k: gas.temperature,
                exponent: 0.0,
                t_min_k: 200.0,
                t_max_k: 500.0,
            },
        }],
    }])
    .unwrap()
    .evaluate(&RelaxingMedium::new(gas, vec![(species, 1.0)]).unwrap())
    .unwrap();
    let inputs = segments();
    let spectrum = estimate(&inputs.map(channels), WelchWindow::Rectangular);
    let molecular = PreparedCompactSourceObserver::new_molecular(
        relaxation,
        spectrum.frequency_hz(BIN).unwrap(),
        RECEIVER,
        &POSITIONS,
        2,
    )
    .unwrap();
    assert!(molecular.wavenumber_per_m().im > 0.0);
    let actual = molecular
        .evaluate_welch_bin(&spectrum, BIN, WORK, || false)
        .unwrap();
    let fields: Vec<_> = inputs
        .iter()
        .map(|s| molecular.evaluate(s).unwrap())
        .collect();
    check_moments(&actual, &fields, 1.0);
    let plane =
        IdealAcousticPlane::new([0.0; 3], [0.0, 0.0, 1.0], IdealPlaneBoundary::Rigid).unwrap();
    let reflected = molecular.with_plane(plane, 4).unwrap();
    let actual = reflected
        .evaluate_welch_bin(
            &spectrum,
            BIN,
            CoherenceBudget {
                max_source_mode_pairs: 12,
                ..WORK
            },
            || false,
        )
        .unwrap();
    let fields: Vec<_> = inputs
        .iter()
        .map(|s| reflected.evaluate(s).unwrap())
        .collect();
    check_moments(&actual, &fields, 1.0);
}
