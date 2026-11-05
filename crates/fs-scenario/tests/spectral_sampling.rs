//! G1/G3: physical frequency quadrature is independent of the observation clock.
use fs_qty::{Dims, QtyAny};
use fs_scenario::{
    RealizationBudget, SpectralDiscretization, SpectrumModel, StochasticEnsemble,
};

const TIME: Dims = Dims([0, 0, 1, 0, 0, 0]);
const RATE: Dims = Dims([0, 0, -1, 0, 0, 0]);

fn ensemble(dt: f64, samples: usize) -> StochasticEnsemble {
    StochasticEnsemble {
        name: "fixed-physical-spectrum".to_string(),
        seed: 0x5eed_1977,
        members: 4096,
        duration: QtyAny::new(samples as f64 * dt, TIME),
        dt: QtyAny::new(dt, TIME),
        model: SpectrumModel::KanaiTajimi {
            s0: 0.02,
            omega_g: QtyAny::new(12.0, RATE),
            zeta_g: 0.55,
        },
    }
}

fn grid() -> SpectralDiscretization {
    SpectralDiscretization {
        omega_max_rad_s: 30.0,
        frequency_lines: 31,
    }
}

fn budget() -> RealizationBudget {
    RealizationBudget {
        max_samples: 1024,
        max_work: 100_000,
    }
}

#[test]
fn g3_clock_refinement_and_record_extension_preserve_the_same_excitation() {
    let base = ensemble(0.04, 80)
        .realize_with_spectrum(7, grid(), budget())
        .unwrap();
    let refined = ensemble(0.02, 160)
        .realize_with_spectrum(7, grid(), budget())
        .unwrap();
    let extended = ensemble(0.04, 160)
        .realize_with_spectrum(7, grid(), budget())
        .unwrap();
    for i in 0..base.values.len() {
        assert_eq!(base.times[i].to_bits(), refined.times[2 * i].to_bits());
        assert_eq!(base.values[i].to_bits(), refined.values[2 * i].to_bits());
        assert_eq!(base.values[i].to_bits(), extended.values[i].to_bits());
    }
    // Preserve the same guarantee for the other actual spectral model.
    let mut wind = ensemble(0.04, 80);
    wind.model = SpectrumModel::Dryden {
        sigma: QtyAny::new(2.0, Dims([1, 0, -1, 0, 0, 0])),
        length_scale: QtyAny::new(8.0, Dims([1, 0, 0, 0, 0, 0])),
        mean_speed: QtyAny::new(10.0, Dims([1, 0, -1, 0, 0, 0])),
    };
    let coarse = wind.realize_with_spectrum(11, grid(), budget()).unwrap();
    wind.dt.value *= 0.5;
    let fine = wind.realize_with_spectrum(11, grid(), budget()).unwrap();
    for (i, value) in coarse.values.iter().enumerate() {
        assert_eq!(value.to_bits(), fine.values[2 * i].to_bits());
    }
}

#[test]
fn g1_gaussian_variance_and_lag_covariance_follow_angular_psd_quadrature() {
    let model = ensemble(0.025, 4);
    let spectrum = SpectralDiscretization {
        omega_max_rad_s: 30.0,
        frequency_lines: 5,
    };
    let dw = spectrum.omega_max_rad_s / spectrum.frequency_lines as f64;
    let lag = 3.0 * model.dt.value;
    let mut expected_variance = 0.0;
    let mut expected_covariance = 0.0;
    for k in 0..spectrum.frequency_lines {
        let omega = (k as f64 + 0.5) * dw;
        let variance = model.model.try_psd(omega).unwrap() * dw;
        expected_variance += variance;
        expected_covariance += variance * (omega * lag).cos();
    }
    let (mut mean0, mut mean1, mut second0, mut second1, mut cross) =
        (0.0, 0.0, 0.0, 0.0, 0.0);
    for member in 0..model.members {
        let record = model
            .realize_with_spectrum(member, spectrum, budget())
            .unwrap();
        let (x, y) = (record.values[0], record.values[3]);
        mean0 += x;
        mean1 += y;
        second0 += x * x;
        second1 += y * y;
        cross += x * y;
    }
    let n = f64::from(model.members);
    let variance0 = second0 / n - (mean0 / n).powi(2);
    let variance1 = second1 / n - (mean1 / n).powi(2);
    let covariance = cross / n - (mean0 / n) * (mean1 / n);
    // Bounded statistical oracle over independent Gaussian members, not a
    // per-record periodogram or an assertion that each trace has unit RMS.
    assert!((mean0 / n).abs() < 0.08 * expected_variance.sqrt());
    assert!((mean1 / n).abs() < 0.08 * expected_variance.sqrt());
    assert!((variance0 - expected_variance).abs() < 0.1 * expected_variance);
    assert!((variance1 - expected_variance).abs() < 0.1 * expected_variance);
    assert!((covariance - expected_covariance).abs() < 0.1 * expected_variance);
}

#[test]
fn g3_time_frequency_and_density_rescaling_preserves_physical_samples() {
    let original = ensemble(0.04, 80);
    let a = original.realize_with_spectrum(17, grid(), budget()).unwrap();
    let scale = 16.0;
    let mut changed = original.clone();
    changed.dt.value *= scale;
    changed.duration.value *= scale;
    if let SpectrumModel::KanaiTajimi { s0, omega_g, .. } = &mut changed.model {
        *s0 *= scale;
        omega_g.value /= scale;
    }
    let spectrum = SpectralDiscretization {
        omega_max_rad_s: grid().omega_max_rad_s / scale,
        ..grid()
    };
    let b = changed.realize_with_spectrum(17, spectrum, budget()).unwrap();
    for (x, y) in a.values.iter().zip(&b.values) {
        assert!((x - y).abs() < 1e-12 * x.abs().max(y.abs()).max(1.0));
    }
}

#[test]
fn g0_frequency_aliasing_parameter_bands_and_resource_exhaustion_refuse() {
    let model = ensemble(0.04, 80);
    for invalid in [
        SpectralDiscretization { omega_max_rad_s: f64::NAN, ..grid() },
        SpectralDiscretization { omega_max_rad_s: 0.0, ..grid() },
        SpectralDiscretization { omega_max_rad_s: 80.0, ..grid() },
        SpectralDiscretization { frequency_lines: 0, ..grid() },
        SpectralDiscretization { frequency_lines: usize::MAX, ..grid() },
    ] {
        assert!(model.realize_with_spectrum(0, invalid, budget()).is_err());
    }
    let nyquist = SpectralDiscretization {
        omega_max_rad_s: core::f64::consts::PI / model.dt.value,
        frequency_lines: 8,
    };
    assert!(model.realize_with_spectrum(0, nyquist, budget()).is_ok());
    assert!(model.realize_with_spectrum(model.members, grid(), budget()).is_err());
    let required_work = 80 * 31 + 80 + 31;
    assert!(model.realize_with_spectrum(
        0, grid(), RealizationBudget { max_samples: 80, max_work: required_work }
    ).is_ok());
    assert!(model.realize_with_spectrum(
        0, grid(), RealizationBudget { max_samples: 80, max_work: required_work - 1 }
    ).is_err());
    assert!(model.realize_with_spectrum(
        0, grid(), RealizationBudget { max_samples: 79, ..budget() }
    ).is_err());
    let mut band = model.clone();
    let viscosity = QtyAny::new(1.0, Dims([-1, 1, -1, 0, 0, 0]));
    band.model = SpectrumModel::CarreauBand {
        eta_zero: [viscosity; 2],
        eta_inf: [viscosity; 2],
        lambda: [QtyAny::new(1.0, TIME); 2],
        n: [1.0; 2],
    };
    assert!(band.realize(0).is_ok());
    assert!(band.realize_with_spectrum(0, grid(), budget()).is_err());
}
