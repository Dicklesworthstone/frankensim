//! Full flanged radiation reaches both existing multimodal line authorities.
use super::*;
use fs_material::gas::GasSpec;
use fs_math::det;

fn air(temperature: f64, pressure: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), temperature, pressure).unwrap()
}

fn multiply(matrix: &[C64], vector: &[C64]) -> Vec<C64> {
    let n = vector.len();
    (0..n)
        .map(|i| (0..n).fold(C64::ZERO, |sum, j| sum + matrix[i * n + j] * vector[j]))
        .collect()
}

fn power(pressure: &[C64], flow: &[C64]) -> f64 {
    pressure
        .iter()
        .zip(flow)
        .map(|(p, u)| (*p * u.conj()).re * 0.5)
        .sum()
}

#[test]
fn flanged_line_selector_retains_full_physical_port_work() {
    let gas = air(343.15, 101_325.0);
    let exterior = air(293.15, 101_325.0);
    let (radius, length, n) = (0.04, 0.03, 3);
    let duct = Duct {
        segments: vec![Segment::Cylinder { radius, length }],
    };
    let config = MmLineConfig {
        sample_rate_hz: 48_000,
        n_modes: n,
        extra_slices: 1,
        loss: LossModel::Lossless,
    };
    let load = MmLoad::FlangedModal(&exterior);
    validate_mouth_load(&duct, &gas, &load, config.sample_rate_hz, 65_536).unwrap();
    assert!(load.is_strict());
    let z0 = gas.density * gas.sound_speed / (core::f64::consts::PI * radius * radius);
    let outgoing = [C64::new(1.0, 0.2), C64::new(-0.3, 0.4), C64::new(0.1, -0.2)];
    let mut roots = vec![0.0];
    roots.extend(fs_duct::modal::j1_roots(n - 1).unwrap());
    for ka in [1.2, 4.6, 7.8] {
        let omega = ka * gas.sound_speed / radius;
        let reflection = matrix_reflection_at(&duct, &gas, omega, &config, &load, z0).unwrap();
        let incoming = multiply(&reflection, &outgoing);
        let pressure: Vec<_> = outgoing
            .iter()
            .zip(&incoming)
            .map(|(a, b)| *a + *b)
            .collect();
        let flow: Vec<_> = outgoing
            .iter()
            .zip(&incoming)
            .map(|(a, b)| (*a - *b).scale(1.0 / z0))
            .collect();
        // Carry physical P/U through the exact uniform-cylinder equations,
        // independently of the production reflection recursion.
        let mut mouth_pressure = Vec::new();
        let mut mouth_flow = Vec::new();
        for i in 0..n {
            let k0 = omega / gas.sound_speed;
            let k = C64::from_re(k0 * k0 - (roots[i] / radius).powi(2)).sqrt();
            let zc = C64::from_re(z0 * k0) / k;
            let phase = k.scale(length);
            let (sin, cos) = det::sin_cos(phase.re);
            let forward = C64::new(cos, sin).scale(det::exp(-phase.im));
            let backward = C64::new(cos, -sin).scale(det::exp(phase.im));
            let c = (forward + backward).scale(0.5);
            let s = (forward - backward) * C64::new(0.0, -0.5);
            mouth_pressure.push(c * pressure[i] + C64::new(0.0, 1.0) * zc * s * flow[i]);
            mouth_flow.push(C64::new(0.0, 1.0) * s / zc * pressure[i] + c * flow[i]);
        }
        let radiation =
            fs_duct::flanged::flanged_modal_impedance(&exterior, radius, omega, n).unwrap();
        let expected = multiply(&radiation, &mouth_flow);
        let scale = mouth_pressure.iter().map(|p| p.abs()).fold(1.0, f64::max);
        for (p, expected) in mouth_pressure.iter().zip(expected) {
            assert!((*p - expected).abs() < 5.0e-10 * scale);
        }
        let input_power = power(&pressure, &flow);
        assert!(input_power > 0.0);
        assert!((input_power - power(&mouth_pressure, &mouth_flow)).abs() < 5.0e-10 * input_power);
    }
    let dc = dc_response(&duct, &gas, &config, &load).unwrap();
    assert!(
        dc.impedance_matrix
            .iter()
            .all(|z| z.re.is_finite() && z.im.is_finite())
    );
}

#[test]
fn flanged_line_preflights_exterior_band_and_pressure() {
    let gas = air(293.15, 101_325.0);
    let different_pressure = air(293.15, 90_000.0);
    let duct = Duct {
        segments: vec![Segment::Cylinder {
            radius: 0.06,
            length: 0.1,
        }],
    };
    assert!(
        validate_mouth_load(
            &duct,
            &gas,
            &MmLoad::FlangedModal(&different_pressure),
            8000,
            65_536,
        )
        .is_err()
    );
    assert!(matches!(
        validate_mouth_load(&duct, &gas, &MmLoad::FlangedModal(&gas), 48_000, 65_536,),
        Err(MmLineError::Duct(DuctError::RadiationKaTooLarge { .. }))
    ));
}
