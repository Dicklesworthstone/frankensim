//! Original boundary solve -> measured source covariance -> physical receivers
//! and radiation. This is separate from the analytic retained-column tests.
use fs_bem::{
    helmholtz::{self, DirectivityProjectionOptions, Formulation, Medium},
    near_field::FirstOrder,
    observation_transfers::{
        BemObservation, HarmonicField, HarmonicFieldStatistics, HarmonicWorkBudget,
        PreparedRadiationOptions, ReceiverSpec, SourceId,
    },
    panel3d::SpherePanels,
};
use fs_exec::CancelGate;
use fs_fft::welch::{WelchBudget, WelchConfig, WelchCrossSpectrum, WelchDetrend, WelchWindow};
use fs_math::c64::C64;

fn close(actual: C64, expected: C64) {
    assert!(
        (actual - expected).abs() <= 2e-9 * expected.abs().max(1e-16),
        "actual {actual:?}, expected {expected:?}"
    );
}

#[test]
fn recorded_source_coherence_reaches_original_boundary_fields_microphones_and_radiation() {
    let n = 32;
    let bin = 3;
    let dt = 1.0 / 2000.0;
    let physical = [
        [C64::new(0.001, 0.0004), C64::new(-0.0007, 0.0009)],
        [C64::new(-0.0003, 0.0008), C64::new(0.0011, -0.0002)],
    ];
    let record: Vec<[f64; 2]> = physical
        .iter()
        .flat_map(|a| {
            (0..n).map(move |j| {
                let theta = core::f64::consts::TAU * bin as f64 * j as f64 / n as f64;
                a.map(|v| v.re * theta.cos() + v.im * theta.sin())
            })
        })
        .collect();
    let spectrum = WelchCrossSpectrum::estimate_interleaved(
        &record,
        dt,
        WelchConfig {
            segment_samples: n,
            hop_samples: n,
            window: WelchWindow::Rectangular,
            detrend: WelchDetrend::None,
        },
        WelchBudget {
            maximum_channels: 2,
            maximum_samples_per_channel: 64,
            maximum_transforms: 4,
            maximum_workspace_bytes: 100_000,
        },
        || false,
    )
    .unwrap();
    let omega = core::f64::consts::TAU * spectrum.frequency_hz(bin).unwrap();
    let medium = Medium {
        density: 1.7,
        sound_speed: 298.0,
    };
    let surface = SpherePanels::icosphere(0.025, 0).unwrap();
    let sources = [
        SourceId::StructuralMode("radial".into()),
        SourceId::StructuralMode("signed-z".into()),
    ];
    let traces = [
        vec![1.0; surface.areas().len()],
        surface.normals().iter().map(|v| v[2]).collect::<Vec<_>>(),
    ];
    let receivers = [
        ReceiverSpec {
            id: "near".into(),
            position_m: [0.06, 0.02, 0.08],
        },
        ReceiverSpec {
            id: "side".into(),
            position_m: [0.1, -0.04, 0.09],
        },
    ];
    let owner = BemObservation::new(&surface, medium, &sources, &traces, &receivers, 1e-3).unwrap();
    let solutions = helmholtz::solve_radiation_batch(
        &surface,
        omega / medium.sound_speed,
        medium,
        &traces.iter().map(Vec::as_slice).collect::<Vec<_>>(),
        Formulation::BurtonMiller,
    )
    .unwrap();
    let field = owner
        .prepare_harmonic_from_solutions(omega, &solutions, || false)
        .unwrap();
    let work = HarmonicWorkBudget {
        maximum_products: 64_000_000,
    };
    let mut observed = [HarmonicFieldStatistics::default(); 2];
    field
        .evaluate_welch_bin(&spectrum, bin, &sources, work, &mut observed, || false)
        .unwrap();
    assert_eq!(field.sources(), &sources);
    assert_eq!(field.receivers(), &receivers);
    let patterns = [
        FirstOrder::default(),
        FirstOrder::new(0.5, [-1.0, 0.0, 0.0]).unwrap(),
    ];
    let mut microphones = field.microphones(&patterns).unwrap();
    let cross = microphones
        .evaluate_welch_bin(&spectrum, bin, &sources, work, || false)
        .unwrap();
    let mut expected_pressure = [0.0; 2];
    let mut expected_intensity = [[C64::ZERO; 3]; 2];
    let mut expected_cross = [C64::ZERO; 4];
    for drive in physical {
        let mut coherent = [HarmonicField::default(); 2];
        field
            .evaluate_coherent(&drive, &mut coherent, || false)
            .unwrap();
        let mut m = [C64::ZERO; 2];
        for receiver in 0..2 {
            let f = coherent[receiver];
            expected_pressure[receiver] += 0.25 * f.pressure_pa.norm_sq();
            for axis in 0..3 {
                expected_intensity[receiver][axis] = expected_intensity[receiver][axis]
                    + (f.pressure_pa * f.particle_velocity_m_s[axis].conj()).scale(0.25);
            }
            m[receiver] = patterns[receiver]
                .observe(f.pressure_pa, f.particle_velocity_m_s, medium)
                .unwrap();
        }
        for i in 0..2 {
            for j in 0..2 {
                expected_cross[i * 2 + j] =
                    expected_cross[i * 2 + j] + (m[i] * m[j].conj()).scale(0.25);
            }
        }
    }
    for receiver in 0..2 {
        close(
            C64::from_re(observed[receiver].mean_square_pressure_pa2),
            C64::from_re(expected_pressure[receiver]),
        );
        for axis in 0..3 {
            close(
                observed[receiver].complex_intensity_w_m2[axis],
                expected_intensity[receiver][axis],
            );
        }
    }
    for i in 0..4 {
        close(cross[i], expected_cross[i]);
    }
    let directions = [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 0.0, -1.0]];
    let radiation = owner
        .prepare_radiation_from_solutions(
            omega,
            &solutions,
            &directions,
            PreparedRadiationOptions {
                angular: DirectivityProjectionOptions::default(),
                initial_latitudes: 8,
                refinements: 4,
                maximum_kernel_evaluations: 50_000_000,
                maximum_source_entries: 1_000_000,
            },
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    let mut density = [0.0; 3];
    let measured = radiation
        .evaluate_welch_bin(&spectrum, bin, &sources, work, &mut density, || false)
        .unwrap();
    let mut total = 0.0;
    let mut first = [0.0; 3];
    let mut reference_density = [0.0; 3];
    for drive in physical {
        let mut far = [C64::ZERO; 3];
        let power = radiation
            .evaluate_coherent(&drive, work, &mut far, || false)
            .unwrap();
        total += 0.5 * power.total_power_w;
        for axis in 0..3 {
            first[axis] += 0.5 * power.first_directional_moment_w[axis];
            reference_density[axis] +=
                0.25 * far[axis].norm_sq() / (medium.density * medium.sound_speed);
        }
    }
    close(C64::from_re(measured.total_power_w), C64::from_re(total));
    for i in 0..3 {
        close(
            C64::from_re(measured.first_directional_moment_w[i]),
            C64::from_re(first[i]),
        );
        close(C64::from_re(density[i]), C64::from_re(reference_density[i]));
    }
}
