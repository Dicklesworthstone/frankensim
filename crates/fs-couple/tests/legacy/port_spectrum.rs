//! G1/G3/G4: physical nonlinear outputs, phase conventions, and bounded work.
use fs_couple::phs_orbit::spectrum::{
    PortSpectrumControls, PortSpectrumError, observe_port_spectrum,
};
use fs_math::c64::C64;
use fs_orbit::HbOrbit;
use fs_phs::{PortHamiltonian, Storage};
use std::f64::consts::{PI, TAU};

struct Spring {
    k: f64,
    beta: f64,
}
impl Storage for Spring {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        0.5 * self.k * x[0].powi(2) + 0.25 * self.beta * x[0].powi(4) + 0.5 * x[1].powi(2)
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        out[0] = self.k * x[0] + self.beta * x[0].powi(3);
        out[1] = x[1];
    }
}
fn system(k: f64, beta: f64) -> PortHamiltonian {
    PortHamiltonian::new(
        2,
        2,
        vec![0.0, 1.0, -1.0, 0.0],
        vec![0.0; 4],
        vec![1.0, 0.0, 0.0, 1.0],
        Box::new(Spring { k, beta }),
    )
    .unwrap()
}
fn phase(angle: f64) -> C64 {
    C64::new(angle.cos(), angle.sin())
}

// A manufactured full-state waveform, not a claim that this forced trajectory
// solves the unforced oscillator. This isolates the physical observation law.
fn waveform(bias: f64, amplitude: f64, angle: f64) -> HbOrbit {
    HbOrbit {
        omega: TAU * 83.0,
        coeffs: vec![
            vec![C64::from_re(bias), phase(angle).scale(0.5 * amplitude)],
            vec![C64::from_re(-0.2), C64::new(0.0, -0.15)],
        ],
        residual_trace: Vec::new(),
        residual: f64::NAN,
        aft_samples: 0,
        sampling_defect: f64::NAN,
    }
}
fn controls(harmonics: usize) -> PortSpectrumControls {
    PortSpectrumControls {
        harmonics,
        initial_samples: 32,
        max_samples: 256,
        max_source_evaluations: 480,
        max_workspace_bytes: 1 << 20,
        rms_tolerances: vec![1e-10; 2],
    }
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 3e-11 * (1.0 + expected.abs()),
        "actual {actual:e}, expected {expected:e}"
    );
}

#[test]
fn original_nonlinear_storage_creates_dc_even_and_third_harmonics_with_correct_phase() {
    let (k, beta, bias, amplitude, angle) = (3.0, 7.0, 0.4, 0.8, 0.37);
    let source = system(k, beta);
    let candidate = waveform(bias, amplitude, angle);
    let result = observe_port_spectrum(&source, &candidate, &controls(3), || false).unwrap();
    close(
        result.dc[0],
        k * bias + beta * (bias.powi(3) + 1.5 * bias * amplitude.powi(2)),
    );
    close(result.dc[1], -0.2);
    let amplitudes = [
        k * amplitude + beta * (3.0 * bias.powi(2) * amplitude + 0.75 * amplitude.powi(3)),
        1.5 * beta * bias * amplitude.powi(2),
        0.25 * beta * amplitude.powi(3),
    ];
    for (index, amplitude) in amplitudes.into_iter().enumerate() {
        let expected = phase(-((index + 1) as f64) * angle).scale(amplitude);
        close(result.peak_phasors[0][index].re, expected.re);
        close(result.peak_phasors[0][index].im, expected.im);
    }
    close(result.peak_phasors[1][0].re, 0.0);
    close(result.peak_phasors[1][0].im, 0.3);
    assert_eq!(result.refinement_trace.len(), 3);
    assert_eq!(result.source_evaluations, 32 + 64 + 128);
    assert!(result.refinement_trace[0].sampling_change_rms.is_none());
    for check in &result.refinement_trace[1..] {
        assert!(check.omitted_rms.iter().all(|&r| r < 1e-10));
        assert!(
            check
                .sampling_change_rms
                .as_ref()
                .unwrap()
                .iter()
                .all(|&r| r < 1e-10)
        );
    }
    assert_eq!(
        candidate.coeffs[0].len(),
        2,
        "observation must not alter state truncation"
    );
}

#[test]
fn omitted_physical_harmonics_refuse_instead_of_following_the_state_order() {
    let source = system(3.0, 7.0);
    let candidate = waveform(0.0, 0.8, 0.37);
    let error = observe_port_spectrum(&source, &candidate, &controls(1), || false).unwrap_err();
    let PortSpectrumError::Unresolved { last } = error else {
        panic!("{error:?}")
    };
    assert_eq!(last.samples, 256);
    close(
        last.omitted_rms[0],
        0.25 * 7.0 * 0.8_f64.powi(3) / 2.0_f64.sqrt(),
    );
    assert!(last.sampling_change_rms.unwrap()[0] < 1e-10);
}

#[test]
fn physical_scale_does_not_set_an_implicit_harmonic_cutoff() {
    let candidate = waveform(0.4, 0.8, 0.37);
    let ordinary =
        observe_port_spectrum(&system(3.0, 7.0), &candidate, &controls(3), || false).unwrap();
    let mut tiny = controls(3);
    tiny.rms_tolerances[0] = 1e-160;
    let scaled =
        observe_port_spectrum(&system(3e-150, 7e-150), &candidate, &tiny, || false).unwrap();
    close(scaled.dc[0] / 1e-150, ordinary.dc[0]);
    for (a, b) in scaled.peak_phasors[0].iter().zip(&ordinary.peak_phasors[0]) {
        close(a.re / 1e-150, b.re);
        close(a.im / 1e-150, b.im);
    }
}

#[test]
fn cancellation_of_strong_modes_does_not_erase_nonlinear_output_harmonics() {
    struct CoupledOutput;
    impl Storage for CoupledOutput {
        fn hamiltonian(&self, x: &[f64]) -> f64 {
            0.5 * (x[0].powi(2) + x[2].powi(2)) + 0.25 * x[1].powi(4)
        }
        fn gradient(&self, x: &[f64], out: &mut [f64]) {
            out.copy_from_slice(&[x[0], x[1].powi(3), x[2]]);
        }
    }
    let source = PortHamiltonian::new(3, 1, vec![0.0; 9], vec![0.0; 9],
        vec![1.0; 3], Box::new(CoupledOutput)).unwrap();
    let angle = 0.37;
    let half = phase(angle).scale(0.5);
    let large = 2.0_f64.powi(54);
    // A manufactured physical waveform isolates the observation law. The two
    // linear efforts cancel exactly; the remaining cubic effort is cos^3(theta).
    let candidate = HbOrbit {
        omega: TAU * 83.0,
        coeffs: vec![
            vec![C64::ZERO, half.scale(large)],
            vec![C64::ZERO, half],
            vec![C64::ZERO, half.scale(-large)],
        ],
        residual_trace: Vec::new(), residual: f64::NAN,
        aft_samples: 0, sampling_defect: f64::NAN,
    };
    let mut limits = controls(3);
    limits.rms_tolerances = vec![1e-10];
    let result = observe_port_spectrum(&source, &candidate, &limits, || false).unwrap();
    close(result.dc[0], 0.0);
    for (index, amplitude) in [0.75, 0.0, 0.25].into_iter().enumerate() {
        let expected = phase(-((index + 1) as f64) * angle).scale(amplitude);
        close(result.peak_phasors[0][index].re, expected.re);
        close(result.peak_phasors[0][index].im, expected.im);
    }
    assert_eq!(result.source_evaluations, 224);
}

#[test]
fn dc_only_output_retains_the_physical_static_load_without_doubling() {
    let candidate = waveform(0.4, 0.0, 0.0);
    let mut candidate = candidate;
    candidate.coeffs[1][1] = C64::ZERO;
    let result =
        observe_port_spectrum(&system(3.0, 7.0), &candidate, &controls(0), || false).unwrap();
    close(result.dc[0], 3.0 * 0.4 + 7.0 * 0.4_f64.powi(3));
    assert!(result.peak_phasors.iter().all(Vec::is_empty));
}

#[test]
fn malformed_rows_budgets_source_failure_and_cancellation_publish_no_spectrum() {
    let source = system(3.0, 7.0);
    let mut candidate = waveform(0.4, 0.8, 0.37);
    let mut c = controls(3);
    c.max_workspace_bytes = 1;
    assert!(matches!(
        observe_port_spectrum(&source, &candidate, &c, || false),
        Err(PortSpectrumError::Budget(_))
    ));
    c = controls(3);
    c.max_source_evaluations = 223;
    assert!(matches!(
        observe_port_spectrum(&source, &candidate, &c, || false),
        Err(PortSpectrumError::Budget(_))
    ));
    for at in [0, 2, 20, 45, 100] {
        let mut calls = 0;
        assert_eq!(
            observe_port_spectrum(&source, &candidate, &controls(3), || {
                let stop = calls == at;
                calls += 1;
                stop
            })
            .unwrap_err(),
            PortSpectrumError::Cancelled
        );
    }
    candidate.coeffs[1].push(C64::ZERO);
    assert!(matches!(
        observe_port_spectrum(&source, &candidate, &controls(3), || false),
        Err(PortSpectrumError::Invalid(_))
    ));
    candidate = waveform(0.4, 0.8, 0.37);
    candidate.coeffs[0][0].im = 1.0;
    assert!(matches!(
        observe_port_spectrum(&source, &candidate, &controls(3), || false),
        Err(PortSpectrumError::Invalid(_))
    ));
    assert_eq!(
        observe_port_spectrum(
            &system(3.0, f64::NAN),
            &waveform(0.4, 0.8, 0.37),
            &controls(3),
            || false
        )
        .unwrap_err(),
        PortSpectrumError::SourceRefused
    );
}

#[test]
fn physical_harmonics_reach_original_acoustic_propagation_and_pressure_renderer() {
    use fs_aeroac::curle3d::{CompactSourcePhasor, PreparedCompactSourceObserver};
    use fs_aeroac::observation::{
        HarmonicObservationBudget, HarmonicSampleBudget, HarmonicSampleClock, PreparedHarmonicField,
    };
    use fs_couple::pcm_wav::observation::{PressureRenderer, harmonic::HarmonicPressureRenderer};
    use fs_material::fluid_acoustics::FluidAcousticState;
    use fs_math::dd::Dd;

    let (k, beta, bias, amplitude, angle) = (3.0, 7.0, 0.4, 0.8, 0.37);
    let candidate = waveform(bias, amplitude, angle);
    let result =
        observe_port_spectrum(&system(k, beta), &candidate, &controls(3), || false).unwrap();
    let (distance, sound_speed) = (2.3, 340.0);
    let medium = FluidAcousticState::try_new(1.2, sound_speed).unwrap();
    // Explicit one-way physical map: port 0 supplies the force on the fluid
    // along x. Static loading is retained in result.dc and not sent to an
    // acoustic spectral line. No reciprocal source loading is claimed here.
    let lines: Vec<_> = result.peak_phasors[0]
        .iter()
        .enumerate()
        .map(|(h, &force)| {
            let observer = PreparedCompactSourceObserver::new(
                medium,
                (h + 1) as f64 * 83.0,
                [distance, 0.0, 0.0],
                &[[0.0; 3]],
                1,
            )
            .unwrap();
            observer
                .evaluate_harmonic_line(&[CompactSourcePhasor::force_dipole([
                    force,
                    C64::ZERO,
                    C64::ZERO,
                ])])
                .unwrap()
        })
        .collect();
    let field = PreparedHarmonicField::new(
        &lines,
        HarmonicObservationBudget {
            max_lines: 3,
            max_line_pairs: 9,
        },
    )
    .unwrap();
    let start = 0.003;
    let clock = HarmonicSampleClock::new(Dd::from_f64(start), 8000).unwrap();
    let mut renderer = HarmonicPressureRenderer::new(
        &field,
        clock,
        HarmonicSampleBudget {
            max_samples: 128,
            max_tone_samples: 384,
        },
    )
    .unwrap();
    let mut pressure = [0.0; 128];
    renderer.block(&mut pressure).unwrap();
    for (j, &actual) in pressure.iter().enumerate() {
        let time = start + j as f64 / 8000.0 - distance / sound_speed;
        let theta = candidate.omega * time + angle;
        let q = bias + amplitude * theta.cos();
        let qdot = -candidate.omega * amplitude * theta.sin();
        let force = k * q + beta * q.powi(3) - result.dc[0];
        let force_rate = (k + 3.0 * beta * q * q) * qdot;
        let expected =
            (force / (distance * distance) + force_rate / (sound_speed * distance)) / (4.0 * PI);
        close(actual, expected);
    }
}
