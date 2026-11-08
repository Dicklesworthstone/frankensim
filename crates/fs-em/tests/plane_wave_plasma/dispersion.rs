//! Original carrier inertia and mutual drag in bulk group/GVD propagation.
use super::*;

#[test]
fn g1_original_and_coulomb_carriers_supply_native_maxwell_frequency_derivatives() {
    let independent = gas();
    let coupled = independent.clone()
        .with_coulomb_drag(fixture::coulomb_model(independent.state()), || false).unwrap();
    for gas in [&independent, &coupled] {
        let transport = gas.prepare_harmonic(|| false).unwrap();
        let background = background(gas);
        for omega in [1e8_f64, 1e10, 1e12, 1e14] {
            let law = IonizedGasWaveResponse::try_new(
                &transport, &background, omega, || false,
            ).unwrap();
            let derivatives = transport.permittivity_frequency_derivatives(omega, || false).unwrap();
            assert!(std::ptr::eq(derivatives.source, gas));
            assert_eq!(derivatives.omega_rad_s, omega);
            let jet = DielectricResponse::IonizedGas(law).bulk_dispersion(omega).unwrap();
            assert_eq!(jet.permittivity_f_m(), law.effective_permittivity_f_m());
            assert_eq!(jet.scaled_permittivity_derivatives_f_m(),
                derivatives.scaled_permittivity_derivatives_f_m);
            // Five separate original value solves; no new derivative helper
            // participates in this oracle, including the mutual-drag case.
            let step = omega * 3e-4;
            let background_slowness = (EPSILON * MU).sqrt();
            let samples: [[C64; 2]; 5] = core::array::from_fn(|i| {
                let w = omega + (i as f64 - 2.0) * step;
                let sample = IonizedGasWaveResponse::try_new(
                    &transport, &background, w, || false,
                ).unwrap();
                let k = PlanarDielectricInterface::try_new(
                    &background, DielectricResponse::IonizedGas(sample),
                ).unwrap().response(w, 1.0, Polarization::Te).unwrap()
                    .transmitted_normal_wavenumber_rad_m;
                // Rationalize the original carrier excess instead of
                // subtracting the much larger linear background wave. This
                // oracle uses only native conductivity and Maxwell values.
                let sigma = sample.conductivity_s_m();
                let chi = C64::new(sigma.im / w, -sigma.re / w);
                let slowness = C64::new(k.re / w, k.im / w);
                let excess = (chi / (slowness + C64::from_re(background_slowness)))
                    .scale(w * MU);
                [k, excess]
            });
            let values = samples.map(|sample| sample[0]);
            let excess = samples.map(|sample| sample[1]);
            let first = (values[0] - values[1].scale(8.0)
                + values[3].scale(8.0) - values[4]).scale(1.0 / (12.0 * step));
            // The constant background has exactly zero curvature.
            let second = (-excess[0] + excess[1].scale(16.0)
                - excess[2].scale(30.0) + excess[3].scale(16.0) - excess[4])
                .scale(1.0 / (12.0 * step * step));
            complex_near(jet.wavenumber_rad_m(), values[2], 1e-13);
            complex_near(jet.wavenumber_derivative_s_m(), first, 3e-6);
            complex_near(jet.wavenumber_second_derivative_s2_m(), second, 2e-4);
            let distance = 0.05 / jet.wavenumber_rad_m().abs();
            let propagation = jet.propagation(distance).unwrap();
            near(propagation.group_delay_s, first.re * distance, 3e-6);
            near(propagation.group_delay_dispersion_s2, second.re * distance, 2e-4);
        }
    }
}

#[test]
fn g0_carrier_jets_retain_cancellation_and_prepared_frequency_admission() {
    let gas = gas();
    let transport = gas.prepare_harmonic(|| false).unwrap();
    for omega in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(transport.permittivity_frequency_derivatives(omega, || false).is_err());
    }
    assert!(matches!(transport.permittivity_frequency_derivatives(1e10, || true),
        Err(IonizedGasError::Cancelled)));
    let background = background(&gas);
    let law = IonizedGasWaveResponse::try_new(
        &transport, &background, 1e10, || false,
    ).unwrap();
    assert!(DielectricResponse::IonizedGas(law).bulk_dispersion(1.01e10).is_err());
    let jet = DielectricResponse::IonizedGas(law).bulk_dispersion(1e10).unwrap();
    assert!(std::ptr::eq(jet.source().material(), &background));
}
