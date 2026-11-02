use super::*;
use crate::compose::ComposeOptions;
use crate::porous::PorousLayer;
use fs_material::gas::{GasSpec, GasState};
use fs_material::porous::{PorousCard, PorousFluidLaw};

fn matrix(a: f64, b: f64, cross: C64) -> [C64; 4] {
    [C64::from_re(a), cross, cross.conj(), C64::from_re(b)]
}
fn close(a: f64, b: f64, relative: f64) {
    assert!(
        (a - b).abs() <= relative * a.abs().max(b.abs()),
        "{a:.17e} != {b:.17e}"
    );
}
fn reconstructed(source: PressureCovariance) -> [C64; 4] {
    let mut result = [C64::ZERO; 4];
    for p in source.peak_components() {
        for i in 0..2 {
            for j in 0..2 {
                // Form the Hermitian outer product explicitly. The native
                // complex multiply uses one FMA and one rounded product; for
                // p*conj(p) that can leave a tiny imaginary diagonal residual.
                // A covariance diagonal is real by definition, independently
                // of that generic complex-product rounding artifact.
                let h = p[i].scale(0.5);
                result[2 * i + j].re += h.re * p[j].re + h.im * p[j].im;
                result[2 * i + j].im += h.im * p[j].re - h.re * p[j].im;
            }
        }
    }
    result
}

#[test]
fn g0_covariance_refuses_nonhermitian_indefinite_and_nonfinite_data_without_projection() {
    let good = matrix(4.0, 1.0, C64::new(1.0, 0.5));
    assert!(PressureCovariance::try_new(good).is_ok());
    let mut invalid = vec![
        matrix(-1.0, 1.0, C64::ZERO),
        matrix(0.0, 1.0, C64::from_re(1e-200)),
        matrix(
            1.0,
            1.0,
            C64::from_re(f64::from_bits(1.0_f64.to_bits() + 1)),
        ),
        matrix(f64::INFINITY, 1.0, C64::ZERO),
        matrix(1.0, 1.0, C64::new(1.0, 1e-170)),
        matrix(1.0, 1.0, C64::new(1.0, 1e-200)),
        matrix(
            2.0_f64.powi(900),
            2.0_f64.powi(900),
            C64::new(2.0_f64.powi(900), 1e-200),
        ),
        matrix(1.0, 1.0, C64::new(f64::NAN, 0.0)),
    ];
    let mut asymmetric = good;
    asymmetric[2] = good[1];
    invalid.push(asymmetric);
    let mut imaginary_diagonal = good;
    imaginary_diagonal[0].im = 1e-300;
    invalid.push(imaginary_diagonal);
    for bad in invalid {
        assert!(PressureCovariance::try_new(bad).is_err());
    }
    for weight in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(PressureCovariance::from_angular_density(good, weight).is_err());
    }
    let zero = PressureCovariance::try_new([C64::ZERO; 4]).unwrap();
    assert_eq!(zero.peak_components(), &[[C64::ZERO; 2]]);
    for source in [matrix(0.0, 7.0, C64::ZERO), matrix(7.0, 0.0, C64::ZERO)] {
        let covariance = PressureCovariance::try_new(source).unwrap();
        for (actual, expected) in reconstructed(covariance).into_iter().zip(source) {
            close(actual.re, expected.re, 4e-16);
            close(actual.im, expected.im, 4e-16);
        }
    }
}

#[test]
fn g1_scaled_covariance_factor_keeps_weak_incoherence_and_extreme_variance_ratios() {
    let c = f64::from_bits(1.0_f64.to_bits() - 1);
    let source = PressureCovariance::try_new(matrix(1.0, 1.0, C64::from_re(c))).unwrap();
    assert_eq!(source.peak_components().len(), 2);
    let weak = source.peak_components()[1][1].norm_sq() * 0.5;
    close(weak, (1.0 - c) * (1.0 + c), 4e-16);
    assert!(weak > 0.0);
    let edge = 2.0_f64.powi(-26);
    let imaginary = f64::from_bits(edge.to_bits() - 1);
    let b = 1.0 + f64::EPSILON;
    let almost = PressureCovariance::try_new(matrix(1.0, b, C64::new(1.0, imaginary))).unwrap();
    let weak = almost.peak_components()[1][0].norm_sq() * 0.5;
    close(weak, (edge - imaginary) * (edge + imaginary) / b, 8e-16);
    assert!(weak > 0.0);
    for authored in [
        matrix(1e300, 1e-300, C64::ZERO),
        matrix(2.0_f64.powi(900), 2.0_f64.powi(-900), C64::new(0.3, 0.4)),
        matrix(2.0_f64.powi(-900), 2.0_f64.powi(900), C64::new(-0.3, 0.4)),
        matrix(4.0, 1.0, C64::new(1.0, 0.5)),
        matrix(f64::from_bits(1), f64::from_bits(1), C64::ZERO),
    ] {
        let actual = reconstructed(PressureCovariance::try_new(authored).unwrap());
        for (a, b) in actual.into_iter().zip(authored) {
            assert!(
                (a.re - b.re).abs() <= f64::from_bits(2) + 8e-16 * b.re.abs(),
                "real: {a:?} != {b:?}"
            );
            assert!(
                (a.im - b.im).abs() <= f64::from_bits(2) + 8e-16 * b.im.abs(),
                "imaginary: {a:?} != {b:?}"
            );
        }
    }
    let p = [C64::new(3.1, 4.7), C64::new(-5.2, 12.3)];
    assert_eq!(
        PressureCovariance::from_peak_pressure(p)
            .unwrap()
            .peak_components(),
        &[p]
    );
}

#[test]
fn g1_same_frequency_covariance_fields_match_an_independent_four_phase_ensemble() {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 101325.0).unwrap();
    let material = PorousFluidLaw::try_new(PorousCard::melamine_uf_nguyen2024(), 1.81e-5).unwrap();
    let layers = [0.03125, 0.046875].map(|length| {
        PorousLayer::try_new(
            &material,
            &gas,
            core::f64::consts::TAU * 1100.0,
            length,
            0.004,
            gas.characteristic_impedance / 0.004,
        )
        .unwrap()
    });
    let stack = PorousLayerStack::try_new(&layers, ComposeOptions::default(), || false).unwrap();
    // These two RMS vectors produce the explicitly authored covariance below.
    // Averaging four coherent phases removes cross terms exactly; neither the
    // ensemble nor its original field evaluation uses the covariance factor.
    let u = [C64::new(3.0, 4.0), C64::new(1.0, -2.0)];
    let v = [C64::new(-2.0, 1.0), C64::from_re(4.0)];
    let authored = matrix(30.0, 21.0, C64::new(-13.0, 14.0));
    let source = PressureCovariance::try_new(authored).unwrap();
    let components = stack.prepare_covariance(source, || false).unwrap();
    assert_eq!(components.len(), 2);
    let phases = [C64::ONE, C64::new(0.0, 1.0), -C64::ONE, C64::new(0.0, -1.0)];
    let ensemble: Vec<_> = phases
        .into_iter()
        .map(|phase| {
            let p = std::array::from_fn(|i| (u[i] + phase * v[i]).scale(core::f64::consts::SQRT_2));
            stack.prepare_excitation(p, || false).unwrap()
        })
        .collect();
    let cov_heat: f64 = components
        .iter()
        .map(|e| e.dissipation().absorbed_power_w)
        .sum();
    close(
        cov_heat,
        ensemble
            .iter()
            .map(|e| 0.25 * e.dissipation().absorbed_power_w)
            .sum(),
        2e-14,
    );
    close(
        components
            .iter()
            .map(|e| e.boundary_power().outgoing_power_w)
            .sum(),
        ensemble
            .iter()
            .map(|e| 0.25 * e.boundary_power().outgoing_power_w)
            .sum(),
        2e-14,
    );
    close(
        components
            .iter()
            .map(|e| e.boundary_power().incident_power_w)
            .sum(),
        51.0 / stack.two_port().reference_impedance,
        8e-16,
    );
    for depth in [0.0, 0.017, 0.03125, 0.057, 0.078125] {
        let moments = |fields: &[PorousStackExcitation<'_, '_, '_>], weight: f64| {
            let mut result = [0.0; 4];
            for e in fields {
                let field = e.at_depth(depth).unwrap().field;
                for (sum, value) in result.iter_mut().zip([
                    0.5 * field.pressure_pa.norm_sq(),
                    0.5 * field.superficial_velocity_m_s.norm_sq(),
                    field.viscous_dissipation_w_m3,
                    field.thermal_dissipation_w_m3,
                ]) {
                    *sum += weight * value;
                }
            }
            result
        };
        for (a, b) in moments(&components, 1.0)
            .into_iter()
            .zip(moments(&ensemble, 0.25))
        {
            close(a, b, 3e-14);
        }
    }
    let weighted = PressureCovariance::from_angular_density(authored, 12.5).unwrap();
    let weighted = stack.prepare_covariance(weighted, || false).unwrap();
    close(
        weighted
            .iter()
            .map(|e| e.dissipation().absorbed_power_w)
            .sum(),
        12.5 * cov_heat,
        2e-14,
    );
    assert!(stack.prepare_covariance(source, || true).is_err());
    let original = stack.prepare_excitation(u, || false).unwrap();
    let coherent = stack
        .prepare_covariance(PressureCovariance::from_peak_pressure(u).unwrap(), || false)
        .unwrap();
    assert_eq!(original.dissipation(), coherent[0].dissipation());
    assert_eq!(
        original.boundary_power().outgoing_pressure_pa,
        coherent[0].boundary_power().outgoing_pressure_pa
    );
}
