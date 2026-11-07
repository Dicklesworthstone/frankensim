//! Exact control-volume deposition retains Maxwell power and signed pools.
use super::*;

#[test]
fn finite_slab_heat_matches_flux_loss_local_fields_and_disjoint_cells() {
    let gas = gas()
        .with_coulomb_drag(fixture::coulomb_model(gas().state()), || false)
        .unwrap();
    let harmonic = gas.prepare_harmonic(|| false).unwrap();
    let background = background(&gas);
    let omega = 1e12;
    let transport = harmonic
        .prepare_magnetized(
            omega,
            [0.0, 0.0, 0.6 * omega * fixture::ME / fixture::E],
            magnetic_limits(),
            || false,
        )
        .unwrap();
    let interface = MagnetizedIonizedGasInterface::try_new(
        DielectricResponse::Lossless(&background),
        &transport,
        &background,
        [0.0, 0.0, 1.0],
        || false,
    )
    .unwrap();
    let wave = interface
        .prepare_wave(
            [C64::new(0.02, 0.003), C64::new(-0.004, 0.01), C64::ZERO],
            || false,
        )
        .unwrap();
    let decay = interface
        .circular_modes()
        .iter()
        .map(|m| -m.transmitted_normal_wavenumber_rad_m.im)
        .fold(0.0_f64, f64::max);
    let (area, start, end) = (0.7, 0.2 / decay, 0.9 / decay);
    let slab = wave.slab_absorption(area, start, end).unwrap();
    assert!(std::ptr::eq(slab.wave(), &wave));
    assert_eq!(slab.volume_m3(), area * (end - start));
    near(
        slab.absorbed_power_w(),
        area * (wave.field(start).unwrap().complex_poynting_w_m2[2].re
            - wave.field(end).unwrap().complex_poynting_w_m2[2].re),
        2e-11,
    );
    let split = 0.43 * start + 0.57 * end;
    let left = wave.slab_absorption(area, start, split).unwrap();
    let right = wave.slab_absorption(area, split, end).unwrap();
    near(
        slab.absorbed_power_w(),
        left.absorbed_power_w() + right.absorbed_power_w(),
        2e-13,
    );
    near(
        slab.electron_power_w(),
        left.electron_power_w() + right.electron_power_w(),
        2e-13,
    );
    near(
        slab.heavy_power_w(),
        left.heavy_power_w() + right.heavy_power_w(),
        2e-13,
    );
    near(
        slab.exchange_correction_to_electrons_w(),
        left.exchange_correction_to_electrons_w() + right.exchange_correction_to_electrons_w(),
        2e-13,
    );
    near(
        slab.electron_power_w(),
        slab.electron_friction_power_w() + slab.exchange_correction_to_electrons_w(),
        2e-13,
    );
    near(
        slab.heavy_power_w(),
        slab.heavy_friction_power_w() - slab.exchange_correction_to_electrons_w(),
        2e-13,
    );
    assert_ne!(slab.exchange_correction_to_electrons_w(), 0.0);
    assert!(slab.power_balance_relative_residual() < 1e-10);
    // Independent composite Simpson integration of actual local carrier fields.
    let n = 512;
    let step = (end - start) / n as f64;
    let mut integral = [0.0; 6];
    for i in 0..=n {
        let p = wave.field(start + i as f64 * step).unwrap();
        let w = if i == 0 || i == n {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        let values = [
            p.mean_heating_w_m3(),
            p.mean_electron_heating_w_m3(),
            p.mean_heavy_heating_w_m3(),
            p.mean_exchange_correction_to_electrons_w_m3(),
            p.mean_drift_kinetic_energy_j_m3(),
            p.mean_background_field_energy_j_m3(),
        ];
        for j in 0..6 {
            integral[j] += w * values[j];
        }
    }
    let actual = [
        slab.absorbed_power_w(),
        slab.electron_power_w(),
        slab.heavy_power_w(),
        slab.exchange_correction_to_electrons_w(),
        slab.drift_kinetic_energy_j(),
        slab.background_field_energy_j(),
    ];
    for i in 0..6 {
        near(actual[i], integral[i] * area * step / 3.0, 2e-11);
    }
}

#[test]
fn weak_opaque_deep_and_zero_amplitude_limits_preserve_physical_heat() {
    let gas = gas();
    let harmonic = gas.prepare_harmonic(|| false).unwrap();
    let background = background(&gas);
    let omega = 1e12;
    let transport = harmonic
        .prepare_magnetized(omega, [0.0; 3], magnetic_limits(), || false)
        .unwrap();
    let interface = MagnetizedIonizedGasInterface::try_new(
        DielectricResponse::Lossless(&background),
        &transport,
        &background,
        [0.0, 0.0, 1.0],
        || false,
    )
    .unwrap();
    let wave = interface
        .prepare_wave([C64::from_re(0.01), C64::ZERO, C64::ZERO], || false)
        .unwrap();
    let alpha = -interface.circular_modes()[0]
        .transmitted_normal_wavenumber_rad_m
        .im;
    let width = 1e-18 / alpha;
    let thin = wave.slab_absorption(3.0, 0.0, width).unwrap();
    assert!(thin.absorbed_power_w() > 0.0);
    near(
        thin.absorbed_power_w(),
        3.0 * width * wave.field(0.0).unwrap().mean_heating_w_m3(),
        3e-13,
    );
    let opaque = wave.slab_absorption(3.0, 0.0, 1000.0 / alpha).unwrap();
    near(
        opaque.absorbed_power_w(),
        3.0 * wave.transmitted_power_w_m2(),
        2e-11,
    );
    assert_eq!(
        wave.slab_absorption(3.0, 900.0 / alpha, 1000.0 / alpha)
            .unwrap()
            .absorbed_power_w(),
        0.0
    );
    let zero = interface.prepare_wave([C64::ZERO; 3], || false).unwrap();
    let empty = zero.slab_absorption(3.0, 0.0, 1.0).unwrap();
    assert_eq!(empty.absorbed_power_w(), 0.0);
    assert_eq!(empty.electron_power_w(), 0.0);
    assert_eq!(empty.heavy_power_w(), 0.0);
    assert_eq!(empty.drift_kinetic_energy_j(), 0.0);
    for (area, start, end) in [
        (0.0, 0.0, 1.0),
        (-1.0, 0.0, 1.0),
        (1.0, -1.0, 1.0),
        (1.0, 1.0, 1.0),
        (1.0, 2.0, 1.0),
        (1.0, 0.0, f64::INFINITY),
        (f64::NAN, 0.0, 1.0),
        (f64::MAX, 0.0, f64::MAX),
    ] {
        assert!(wave.slab_absorption(area, start, end).is_err());
    }
}
