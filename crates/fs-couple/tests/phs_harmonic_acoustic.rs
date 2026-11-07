//! G1/G3/G4: nonlinear physical harmonics reach actual BEM fields and PCM.
use fs_aeroac::observation::{
    AcousticFieldComponent, HarmonicProductSpectrumLimits, HarmonicSampleBudget,
    HarmonicSampleClock,
};
use fs_bem::{
    helmholtz::Medium,
    near_field::Options,
    observation_transfers::{
        CompactSourceSpec, ForceDipoleSpec, HarmonicField, PointSourceSpec,
        PreparedHarmonicObservation, ReceiverSpec, ScatteringObservation,
    },
    panel3d::SpherePanels,
};
use fs_couple::{
    pcm_wav::{
        encode_pcm16_wav,
        observation::harmonic::HarmonicPressureRenderer,
        stream::{Pcm16WavStream, render_pressure_pcm16},
    },
    phs_orbit::{
        AutonomousPortFeedback, PeriodicPortInput, PhsOrbitError,
        acoustic::{
            PortAcousticBinding, PortAcousticError, PortAcousticLimits, PortAcousticSpectrum,
        },
        audit_autonomous_feedback_orbit, audit_forced_orbit,
        sideband::PhsPortSidebandResponse,
        spectrum::{PortSpectrum, PortSpectrumControls, observe_port_spectrum},
    },
};
use fs_exec::CancelGate;
use fs_fft::harmonic_window::HarmonicWindow;
use fs_math::{c64::C64, dd::Dd, det};
use fs_orbit::{HbOrbit, sideband::SidebandResponse};
use fs_phs::{PortHamiltonian, QuadraticStorage, Storage};
use std::{f64::consts::TAU, io::Cursor};

const LIMITS: PortAcousticLimits = PortAcousticLimits {
    maximum_ports: 8,
    maximum_lines: 32,
    maximum_receivers: 8,
    maximum_products: 1_000_000,
    maximum_workspace_bytes: 1 << 20,
};
// Port zero is force [N]; port one is velocity [m/s]. The second source is
// explicitly a piston of area 0.002 m². No unit/source conversion is inferred.
const MAP: [C64; 4] = [C64::ONE, C64::ZERO, C64::ZERO, C64::new(0.002, 0.0)];

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "actual {actual:e}, expected {expected:e}"
    );
}
fn close_complex(actual: C64, expected: C64) {
    close(actual.re, expected.re, 4e-11);
    close(actual.im, expected.im, 4e-11);
}
fn receivers() -> [ReceiverSpec; 2] {
    [
        ReceiverSpec {
            id: "front".into(),
            position_m: [0.12, 0.08, 0.03],
        },
        ReceiverSpec {
            id: "side".into(),
            position_m: [-0.08, 0.09, 0.02],
        },
    ]
}
fn banks(frequencies: &[f64], receivers: &[ReceiverSpec]) -> Vec<PreparedHarmonicObservation> {
    let surface = SpherePanels::icosphere(0.02, 0).unwrap();
    let sources = [
        CompactSourceSpec::ForceDipole(ForceDipoleSpec {
            id: "applied fluid force".into(),
            position_m: [0.07, 0.0, 0.0],
            direction: [1.0, 0.0, 0.0],
        }),
        CompactSourceSpec::VolumeFlow(PointSourceSpec {
            id: "piston volume flow".into(),
            position_m: [0.0, 0.06, 0.0],
        }),
    ];
    let owner = ScatteringObservation::new_mixed(
        &surface,
        Medium {
            density: 1.2,
            sound_speed: 343.0,
        },
        &sources,
        receivers,
        10_000.0,
        0.001,
        Options::default(),
        || false,
    )
    .unwrap();
    frequencies
        .iter()
        .map(|&omega| owner.prepare_harmonic(omega, || false).unwrap())
        .collect()
}
fn bindings<'a>(
    banks: &'a [PreparedHarmonicObservation],
    map: &'a [C64],
) -> Vec<PortAcousticBinding<'a>> {
    banks
        .iter()
        .map(|observation| PortAcousticBinding {
            observation,
            source_from_ports: map,
        })
        .collect()
}
fn native_fields(
    bank: &PreparedHarmonicObservation,
    force: C64,
    velocity: C64,
) -> Vec<HarmonicField> {
    let mut fields = vec![HarmonicField::default(); bank.receivers().len()];
    bank.evaluate_coherent(&[force, velocity.scale(0.002)], &mut fields, || false)
        .unwrap();
    fields
}
fn sample(amplitude: C64, omega: f64, time: f64) -> f64 {
    let (sin, cos) = (omega * time).sin_cos();
    amplitude.re * cos + amplitude.im * sin
}

struct Spring;
impl Storage for Spring {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        1.5 * x[0].powi(2) + 1.75 * x[0].powi(4) + 0.5 * x[1].powi(2)
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        out[0] = 3.0 * x[0] + 7.0 * x[0].powi(3);
        out[1] = x[1];
    }
}
fn nonlinear_output() -> PortSpectrum {
    nonlinear_output_at(TAU * 83.0)
}
fn nonlinear_output_at(omega: f64) -> PortSpectrum {
    let source = PortHamiltonian::new(
        2,
        2,
        vec![0.0, 1.0, -1.0, 0.0],
        vec![0.0; 4],
        vec![1.0, 0.0, 0.0, 1.0],
        Box::new(Spring),
    )
    .unwrap();
    // Manufactured full-state motion isolates the actual nonlinear output law;
    // it does not assert this waveform solves the unforced source dynamics.
    let angle: f64 = 0.37;
    let waveform = HbOrbit {
        omega,
        coeffs: vec![
            vec![
                C64::from_re(0.4),
                C64::new(angle.cos(), angle.sin()).scale(0.4),
            ],
            vec![C64::from_re(-0.2), C64::new(0.0, -0.15)],
        ],
        residual_trace: vec![],
        residual: 0.0,
        aft_samples: 0,
        sampling_defect: 0.0,
    };
    observe_port_spectrum(
        &source,
        &waveform,
        &PortSpectrumControls {
            harmonics: 3,
            initial_samples: 32,
            max_samples: 256,
            max_source_evaluations: 480,
            max_workspace_bytes: 1 << 20,
            rms_tolerances: vec![1e-10; 2],
        },
        || false,
    )
    .unwrap()
}

#[test]
fn nonlinear_even_and_odd_harmonics_reach_physical_bem_receivers_and_pcm() {
    let output = nonlinear_output();
    let spectrum = PortAcousticSpectrum::from_harmonics(&output, LIMITS, || false).unwrap();
    let receivers = receivers();
    let banks = banks(spectrum.angular_frequencies_rad_s(), &receivers);
    let bindings = bindings(&banks, &MAP);
    let acoustic = spectrum
        .propagate(&receivers, &bindings, LIMITS, || false)
        .unwrap();
    assert_eq!(acoustic.receivers(), &receivers);
    assert_eq!(acoustic.dc_port_outputs(), output.dc);
    for (observation, bank) in acoustic.observations().iter().zip(&banks) {
        assert_eq!(observation.producer_hash(), bank.producer_hash());
        assert_eq!(
            observation.producer_uncertainty(),
            bank.producer_uncertainty()
        );
    }
    // Independent expansion of k*q + beta*q³ for q=b+a*cos(wt+angle).
    let (b, a, angle) = (0.4_f64, 0.8_f64, 0.37_f64);
    let force_amplitudes = [
        3.0 * a + 7.0 * (3.0 * b * b * a + 0.75 * a.powi(3)),
        10.5 * b * a * a,
        1.75 * a.powi(3),
    ];
    let expected: Vec<_> = banks
        .iter()
        .enumerate()
        .map(|(index, bank)| {
            let phase = -((index + 1) as f64) * angle;
            let force = C64::new(phase.cos(), phase.sin()).scale(force_amplitudes[index]);
            let velocity = if index == 0 {
                C64::new(0.0, 0.3)
            } else {
                C64::ZERO
            };
            native_fields(bank, force, velocity)
        })
        .collect();
    for (receiver, field) in acoustic.fields().iter().enumerate() {
        assert_eq!(field.frequency_count(), 3);
        for time in [0.0, 0.00019, 0.0023, 0.02] {
            let observed = field.field_at_clock_dd(Dd::from_f64(time)).unwrap();
            let pressure: f64 = expected
                .iter()
                .zip(spectrum.angular_frequencies_rad_s())
                .map(|(f, &w)| sample(f[receiver].pressure_pa, w, time))
                .sum();
            close(observed.pressure_pa, pressure, 8e-11);
            for axis in 0..3 {
                let velocity: f64 = expected
                    .iter()
                    .zip(spectrum.angular_frequencies_rad_s())
                    .map(|(f, &w)| sample(f[receiver].particle_velocity_m_s[axis], w, time))
                    .sum();
                close(observed.particle_velocity_m_s[axis], velocity, 8e-11);
            }
        }
    }
    let clock = HarmonicSampleClock::new(Dd::from_f64(0.00019), 8000).unwrap();
    let mut renderer = HarmonicPressureRenderer::new(
        &acoustic.fields()[0],
        clock,
        HarmonicSampleBudget {
            max_samples: 64,
            max_tone_samples: 192,
        },
    )
    .unwrap();
    let mut stream = Pcm16WavStream::new(Cursor::new(Vec::new()), 8000, 10_000.0, 64).unwrap();
    render_pressure_pcm16(
        &mut renderer,
        &mut stream,
        &CancelGate::new_clock_free(),
        &mut [0.0; 19],
        129,
    )
    .unwrap();
    let (wav, summary) = stream.finish().unwrap();
    let pressure: Vec<_> = (0..129)
        .map(|index| {
            let time = 0.00019 + f64::from(index) / 8000.0;
            expected
                .iter()
                .zip(spectrum.angular_frequencies_rad_s())
                .map(|(f, &w)| sample(f[0].pressure_pa, w, time))
                .sum::<f64>()
        })
        .collect();
    let (expected_wav, clips) = encode_pcm16_wav(&pressure, 8000, 10_000.0).unwrap();
    assert_eq!(wav.into_inner(), expected_wav);
    assert_eq!(summary.samples, 129);
    assert_eq!(summary.clipped_samples, clips as u64);
}

#[test]
fn original_angular_clock_survives_long_epoch_observation_products_and_pcm() {
    let output = nonlinear_output_at(10_000.0 / 3.0);
    let spectrum = PortAcousticSpectrum::from_harmonics(&output, LIMITS, || false).unwrap();
    let receivers = receivers();
    let banks = banks(spectrum.angular_frequencies_rad_s(), &receivers);
    let acoustic = spectrum
        .propagate(&receivers, &bindings(&banks, &MAP), LIMITS, || false)
        .unwrap();
    let native: Vec<_> = banks
        .iter()
        .enumerate()
        .map(|(index, bank)| {
            native_fields(
                bank,
                output.peak_phasors[0][index],
                output.peak_phasors[1][index],
            )
        })
        .collect();
    let field = &acoustic.fields()[0];
    let native_pressure = |time: Dd| {
        native
            .iter()
            .zip(spectrum.angular_frequencies_rad_s())
            .fold(Dd::ZERO, |sum, (fields, &omega)| {
                // Original angular source, with no Hz conversion. At this
                // epoch that round trip changes the phase by order one.
                let (sin, cos) = det::sin_cos_dd(Dd::from_f64(omega) * time);
                let pressure = fields[0].pressure_pa;
                sum + Dd::from_f64(pressure.re) * cos + Dd::from_f64(pressure.im) * sin
            })
            .to_f64()
    };
    let epoch = 1e12;
    let clock = HarmonicSampleClock::new(Dd::from_pair(epoch, 0.00019), 8000).unwrap();
    let pressure: Vec<_> = (0..129)
        .map(|index| native_pressure(clock.time_at(index).unwrap()))
        .collect();
    for (index, &expected) in pressure.iter().enumerate() {
        close(
            field
                .pressure_at_clock_dd(clock.time_at(index as u64).unwrap())
                .unwrap(),
            expected,
            2e-14,
        );
    }
    let product = field
        .product_spectrum(
            AcousticFieldComponent::Pressure,
            AcousticFieldComponent::Pressure,
            HarmonicProductSpectrumLimits {
                maximum_products: 9,
                maximum_generated_terms: 18,
                maximum_workspace_bytes: 16_384,
            },
        )
        .unwrap();
    close(
        product.evaluate(epoch).unwrap(),
        native_pressure(Dd::from_f64(epoch)).powi(2),
        2e-13,
    );
    let window = HarmonicWindow::new(epoch, 0.0013).unwrap();
    let mut expected_square = 0.0;
    let mut expected_intensity = [0.0; 3];
    for (a, &wa) in native.iter().zip(spectrum.angular_frequencies_rad_s()) {
        for (b, &wb) in native.iter().zip(spectrum.angular_frequencies_rad_s()) {
            let kernel = window.product_kernel(wa, wb).unwrap();
            expected_square += kernel
                .mean_real_product(a[0].pressure_pa, b[0].pressure_pa)
                .unwrap();
            for (axis, value) in expected_intensity.iter_mut().enumerate() {
                *value += kernel
                    .mean_real_product(a[0].pressure_pa, b[0].particle_velocity_m_s[axis])
                    .unwrap();
            }
        }
    }
    let observed = field.observe_window(window).unwrap();
    close(observed.mean_square_pressure_pa2, expected_square, 2e-13);
    close(
        product.observe_window(window).unwrap(),
        expected_square,
        2e-13,
    );
    for (actual, expected) in observed
        .mean_intensity_w_m2
        .into_iter()
        .zip(expected_intensity)
    {
        close(actual, expected, 2e-13);
    }
    let mut renderer = HarmonicPressureRenderer::new(
        field,
        clock,
        HarmonicSampleBudget {
            max_samples: 64,
            max_tone_samples: 192,
        },
    )
    .unwrap();
    let mut stream = Pcm16WavStream::new(Cursor::new(Vec::new()), 8000, 10_000.0, 64).unwrap();
    render_pressure_pcm16(
        &mut renderer,
        &mut stream,
        &CancelGate::new_clock_free(),
        &mut [0.0; 19],
        129,
    )
    .unwrap();
    let (wav, _) = stream.finish().unwrap();
    let (expected_wav, _) = encode_pcm16_wav(&pressure, 8000, 10_000.0).unwrap();
    assert_eq!(wav.into_inner(), expected_wav);
}

fn sidebands(carrier: f64) -> PhsPortSidebandResponse {
    // A manufactured observed probe in the native positive-time envelope
    // chart; the bridge consumes physical ports, not the state coefficients.
    PhsPortSidebandResponse {
        omega_rad_s: 1600.0,
        carrier_angular_frequency: carrier,
        state_harmonics: 1,
        output_harmonics: 1,
        state_dimension: 1,
        port_count: 2,
        state_response: SidebandResponse {
            coefficients: vec![C64::ZERO; 3],
            relative_residual: 0.0,
            omitted_residual: 0.0,
        },
        port_coefficients: vec![
            C64::new(1.0, 2.0),
            C64::new(0.3, -0.4),
            C64::new(-0.25, 1.0),
            C64::new(-0.7, 0.2),
            C64::new(0.8, -0.5),
            C64::new(0.1, 0.6),
        ],
        base_residual: 0.0,
        tangent_sampling_defect: 0.0,
        observation_trace: vec![],
        hessian_actions: 0,
    }
}

#[test]
fn half_integer_carrier_combines_conjugate_collisions_before_propagation() {
    let output = sidebands(800.0);
    let spectrum = PortAcousticSpectrum::from_sidebands(&output, LIMITS, || false).unwrap();
    assert_eq!(spectrum.angular_frequencies_rad_s(), &[800.0, 2400.0]);
    let low = [C64::new(0.75, 1.0), C64::new(-0.4, -0.6)];
    let high = [C64::new(0.8, 0.5), C64::new(0.1, -0.6)];
    for (line, expected) in [low, high].iter().enumerate() {
        for (&actual, &expected) in spectrum
            .peak_port_outputs(line)
            .unwrap()
            .iter()
            .zip(expected)
        {
            close_complex(actual, expected);
        }
    }
    let receivers = receivers();
    let banks = banks(spectrum.angular_frequencies_rad_s(), &receivers);
    let bindings = bindings(&banks, &MAP);
    let acoustic = spectrum
        .propagate(&receivers, &bindings, LIMITS, || false)
        .unwrap();
    let expected = [
        native_fields(&banks[0], low[0], low[1]),
        native_fields(&banks[1], high[0], high[1]),
    ];
    for time in [0.0, 0.0007, 0.01] {
        for receiver in 0..2 {
            let actual = acoustic.fields()[receiver]
                .field_at_clock_dd(Dd::from_f64(time))
                .unwrap();
            let pressure = sample(expected[0][receiver].pressure_pa, 800.0, time)
                + sample(expected[1][receiver].pressure_pa, 2400.0, time);
            close(actual.pressure_pa, pressure, 2e-11);
        }
    }
}

#[test]
fn inharmonic_sidebands_retain_real_beats_and_vector_intensity() {
    let output = sidebands(592.1);
    let spectrum = PortAcousticSpectrum::from_sidebands(&output, LIMITS, || false).unwrap();
    assert_eq!(spectrum.angular_frequencies_rad_s().len(), 3);
    let receivers = receivers();
    let banks = banks(spectrum.angular_frequencies_rad_s(), &receivers);
    // The first source has an explicitly authored quadrature phase. This is
    // a complex positive-frequency law applied AFTER folding the real probe.
    let phased_map = [C64::new(0.0, 1.0), C64::ZERO, C64::ZERO, MAP[3]];
    let bindings = bindings(&banks, &phased_map);
    let acoustic = spectrum
        .propagate(&receivers, &bindings, LIMITS, || false)
        .unwrap();
    let expected: Vec<_> = banks
        .iter()
        .enumerate()
        .map(|(i, bank)| {
            let row = spectrum.peak_port_outputs(i).unwrap();
            native_fields(bank, C64::new(-row[0].im, row[0].re), row[1])[0]
        })
        .collect();
    let (start, duration) = (0.00017, 0.00122);
    let observed = acoustic.fields()[0]
        .observe_window(HarmonicWindow::new(start, duration).unwrap())
        .unwrap();
    // Independent time integration of actual propagated p(t) and u(t), not
    // a second spectral-product implementation or pressure-only plane wave.
    let n = 16_384;
    let (mut mean, mut square, mut separate) = (0.0, 0.0, 0.0);
    let mut intensity = [0.0; 3];
    for index in 0..n {
        let t = start + (f64::from(index) + 0.5) * duration / f64::from(n);
        let (mut p, mut u) = (0.0, [0.0; 3]);
        for (field, &w) in expected.iter().zip(spectrum.angular_frequencies_rad_s()) {
            let line = sample(field.pressure_pa, w, t);
            separate += line * line / f64::from(n);
            p += line;
            for axis in 0..3 {
                u[axis] += sample(field.particle_velocity_m_s[axis], w, t);
            }
        }
        mean += p / f64::from(n);
        square += p * p / f64::from(n);
        for axis in 0..3 {
            intensity[axis] += p * u[axis] / f64::from(n);
        }
    }
    close(observed.mean_pressure_pa, mean, 2e-8);
    close(observed.mean_square_pressure_pa2, square, 2e-8);
    for axis in 0..3 {
        close(observed.mean_intensity_w_m2[axis], intensity[axis], 2e-8);
    }
    assert!(
        (square - separate).abs() > 0.01 * separate,
        "independent-tone powers must not replace this coherent finite window"
    );
}

#[test]
fn static_ports_are_preserved_without_inventing_acoustic_dc() {
    let shifted = sidebands(1600.0);
    let spectrum = PortAcousticSpectrum::from_sidebands(&shifted, LIMITS, || false).unwrap();
    assert_eq!(spectrum.dc_port_outputs(), &[1.0, 0.3]);
    assert_eq!(spectrum.angular_frequencies_rad_s(), &[1600.0, 3200.0]);
    // The imaginary zero-frequency coefficient is absent from the real signal.
    let constant = PortSpectrum {
        omega_rad_s: 1600.0,
        dc: vec![2.0, -1.0],
        peak_phasors: vec![vec![], vec![]],
        refinement_trace: vec![],
        source_evaluations: 0,
    };
    let spectrum = PortAcousticSpectrum::from_harmonics(&constant, LIMITS, || false).unwrap();
    let receivers = receivers();
    let acoustic = spectrum
        .propagate(&receivers, &[], LIMITS, || false)
        .unwrap();
    assert_eq!(acoustic.dc_port_outputs(), &[2.0, -1.0]);
    assert!(acoustic.observations().is_empty());
    for field in acoustic.fields() {
        assert_eq!(field.frequency_count(), 0);
        assert_eq!(field.pressure_at_clock_dd(Dd::from_f64(7.0)).unwrap(), 0.0);
    }
}

#[test]
fn exact_clocks_maps_geometry_budgets_and_cancellation_are_required() {
    let output = sidebands(800.0);
    let mut malformed = sidebands(800.0);
    malformed.output_harmonics = usize::MAX / 2;
    malformed.port_count = usize::MAX;
    assert!(matches!(
        PortAcousticSpectrum::from_sidebands(&malformed, LIMITS, || false),
        Err(PortAcousticError::Budget(_))
    ));
    let mut limits = LIMITS;
    limits.maximum_lines = 2; // All three signed bands must be admitted first.
    assert!(matches!(
        PortAcousticSpectrum::from_sidebands(&output, limits, || false),
        Err(PortAcousticError::Budget(_))
    ));
    assert!(matches!(
        PortAcousticSpectrum::from_sidebands(&output, LIMITS, || true),
        Err(PortAcousticError::Cancelled)
    ));
    let spectrum = PortAcousticSpectrum::from_sidebands(&output, LIMITS, || false).unwrap();
    let receivers = receivers();
    let banks = banks(spectrum.angular_frequencies_rad_s(), &receivers);
    let mut bindings = bindings(&banks, &MAP);
    bindings.swap(0, 1);
    assert!(matches!(
        spectrum.propagate(&receivers, &bindings, LIMITS, || false),
        Err(PortAcousticError::Invalid(_))
    ));
    bindings.swap(0, 1);
    let mut moved = receivers.clone();
    moved[0].position_m[0] += 0.01;
    assert!(matches!(
        spectrum.propagate(&moved, &bindings, LIMITS, || false),
        Err(PortAcousticError::Invalid(_))
    ));
    assert!(matches!(
        spectrum.propagate(&receivers, &bindings[..1], LIMITS, || false),
        Err(PortAcousticError::Invalid(_))
    ));
    let bad_map = [C64::new(f64::NAN, 0.0); 4];
    let bad_bindings = [
        PortAcousticBinding {
            observation: &banks[0],
            source_from_ports: &bad_map,
        },
        bindings[1],
    ];
    assert!(matches!(
        spectrum.propagate(&receivers, &bad_bindings, LIMITS, || false),
        Err(PortAcousticError::Invalid(_))
    ));
    for limits in [
        PortAcousticLimits {
            maximum_products: 1,
            ..LIMITS
        },
        PortAcousticLimits {
            maximum_workspace_bytes: 1,
            ..LIMITS
        },
    ] {
        assert!(matches!(
            spectrum.propagate(&receivers, &bindings, limits, || false),
            Err(PortAcousticError::Budget(_))
        ));
    }
    let mut calls = 0;
    assert!(matches!(
        spectrum.propagate(&receivers, &bindings, LIMITS, || {
            calls += 1;
            calls == 10
        }),
        Err(PortAcousticError::Cancelled)
    ));
    assert_eq!(
        spectrum
            .propagate(&receivers, &bindings, LIMITS, || false)
            .unwrap()
            .fields()
            .len(),
        2
    );
}

struct FollowFlow;
impl AutonomousPortFeedback for FollowFlow {
    fn port_count(&self) -> usize {
        1
    }
    fn input(&self, _state: &[f64], flows: &[f64], out: &mut [f64]) -> bool {
        out.copy_from_slice(flows);
        true
    }
}
fn audit_wave(order: usize) -> HbOrbit {
    let mut row = vec![C64::ZERO; order + 1];
    row[order] = C64::from_re(0.5);
    HbOrbit {
        omega: 2.0,
        coeffs: vec![row],
        residual_trace: vec![],
        residual: 0.0,
        aft_samples: 0,
        sampling_defect: 0.0,
    }
}

#[test]
fn cycle_work_cannot_sample_unresolved_drive_or_state_harmonics_as_dc() {
    let system = PortHamiltonian::new(
        1,
        1,
        vec![0.0],
        vec![1.0],
        vec![1.0],
        Box::new(QuadraticStorage::new(vec![1.0], 1).unwrap()),
    )
    .unwrap();
    let constant = PeriodicPortInput::new(2.0, vec![vec![C64::ZERO]]).unwrap();
    let high = audit_wave(32);
    assert!(matches!(
        audit_forced_orbit(&system, &constant, &high, 32),
        Err(PhsOrbitError::Invalid { .. })
    ));
    assert!(matches!(
        audit_autonomous_feedback_orbit(&system, &FollowFlow, &high, 32),
        Err(PhsOrbitError::Invalid { .. })
    ));
    let input = PeriodicPortInput::new(2.0, high.coeffs.clone()).unwrap();
    assert!(matches!(
        audit_forced_orbit(&system, &input, &audit_wave(1), 32),
        Err(PhsOrbitError::Invalid { .. })
    ));
    let resolved = audit_autonomous_feedback_orbit(&system, &FollowFlow, &high, 256).unwrap();
    close(resolved.mean_stored_energy_j, 0.25, 1e-12);
    close(resolved.mean_resistance_loss_w, 0.5, 1e-12);
    close(resolved.mean_supplied_power_w, 0.5, 1e-12);
    close(resolved.cycle_power_residual_w, 0.0, 1e-12);
}
