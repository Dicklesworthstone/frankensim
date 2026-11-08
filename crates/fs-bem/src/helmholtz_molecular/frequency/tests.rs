use super::*;
use fs_material::gas::{GasSpec, GasState};
use fs_material::relaxation::{RelaxingMedium, input::read_relaxation};

pub(super) fn model(relaxing: bool) -> MolecularExterior {
    // Synthetic reservoir tests constitutive dispersion, not measured dry air.
    let data = if relaxing {
        "declared 1 synthetic 1 vibration 1 900 1 power 101.325 293.15 0 200 500"
    } else {
        "declared 1 synthetic 0"
    };
    let law = read_relaxation(&mut data.split_ascii_whitespace()).unwrap();
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    let medium = RelaxingMedium::new(gas, vec![(law.species()[0].species.clone(), 1.0)]).unwrap();
    MolecularExterior::new(law.evaluate(&medium).unwrap(), [100.0, 2000.0]).unwrap()
}

pub(super) fn limits(n: usize) -> ComplexOptions {
    ComplexOptions {
        max_panels: n,
        max_bytes: 64 * n * n,
        max_work: 512 * (n * n) as u64,
    }
}

pub(super) fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-12),
        "{actual:?} != {expected:?}"
    );
}

#[test]
fn original_value_parity_and_complete_frequency_slope_match_independent_solves() {
    let surface = SpherePanels::icosphere(0.2, 0).unwrap();
    let n = surface.areas().len();
    let velocity: Vec<_> = (0..n).map(|i| C64::new(0.2 + 0.01 * i as f64, -0.04)).collect();
    let slope: Vec<_> = (0..n).map(|i| C64::new(-0.0001, 0.00002 * i as f64)).collect();
    let omega = 700.0;
    let h = 0.7;
    for relaxing in [false, true] {
        let material = model(relaxing);
        for formulation in [Formulation::PlainCbie, Formulation::BurtonMiller] {
            let response = material.solve_batch_with_frequency_derivative(
                &surface, omega, &[&velocity], Some(&[&slope]), formulation, &limits(n), || false,
            ).unwrap().remove(0);
            let original = material.solve_batch(
                &surface, omega, &[&velocity], formulation, &limits(n), || false,
            ).unwrap().remove(0);
            assert_eq!(response.response.pressure_pa(), original.pressure_pa());
            assert_eq!(response.response.boundary_work_w(), original.boundary_work_w());
            assert_eq!(response.velocity_derivative_m_s_per_rad_s, slope);
            assert_eq!(response.response.source(), material.source());
            let samples: Vec<_> = [-2.0, -1.0, 1.0, 2.0].into_iter().map(|offset| {
                let drive: Vec<_> = velocity.iter().zip(&slope)
                    .map(|(&v, &s)| v + s.scale(offset * h)).collect();
                material.solve_batch(
                    &surface, omega + offset * h, &[&drive], formulation, &limits(n), || false,
                ).unwrap().remove(0)
            }).collect();
            for i in 0..n {
                let expected = (samples[0].pressure_pa()[i] - samples[1].pressure_pa()[i].scale(8.0)
                    + samples[2].pressure_pa()[i].scale(8.0) - samples[3].pressure_pa()[i])
                    .scale(1.0 / (12.0 * h));
                close(response.pressure_derivative_pa_per_rad_s[i], expected, 3e-8);
            }
            let expected_work = (samples[0].boundary_work_w() - 8.0 * samples[1].boundary_work_w()
                + 8.0 * samples[2].boundary_work_w() - samples[3].boundary_work_w()) / (12.0 * h);
            assert!((response.boundary_work_derivative_w_per_rad_s - expected_work).abs()
                < 3e-8 * expected_work.abs());
        }
    }
}

#[test]
fn molecular_sphere_frequency_slope_converges_to_analytic_radial_impedance() {
    let material = model(true);
    let (omega, radius) = (600.0, 0.2);
    let gas = material.source().gas();
    let process = &material.source().processes()[0];
    let frozen = material.source().frozen_cv_j_kg_k();
    let cv0 = frozen + process.heat_capacity_j_kg_k;
    let wt = omega * process.relaxation_time_s;
    let n = C64::new(cv0, -wt * frozen);
    let d = C64::new(cv0 + gas.specific_gas_constant, -wt * (frozen + gas.specific_gas_constant));
    // Independent one-reservoir rational caloric elimination.
    let slowness = (n / d).scale(gas.density / gas.pressure).sqrt();
    let k = slowness.scale(omega);
    let dk = slowness * (C64::ONE + C64::new(0.0,
        0.5 * wt * gas.specific_gas_constant * process.heat_capacity_j_kg_k) / (n * d));
    let denominator = C64::new(-1.0 / radius - k.im, k.re);
    let q = C64::new(0.0, gas.density * omega);
    let expected = C64::new(0.0, gas.density) / denominator
        - q * C64::new(-dk.im, dk.re) / (denominator * denominator);
    let mut errors = Vec::new();
    for level in 0..=2 {
        let surface = SpherePanels::icosphere(radius, level).unwrap();
        let count = surface.areas().len();
        let velocity = vec![C64::ONE; count];
        let response = material.solve_batch_with_frequency_derivative(
            &surface, omega, &[&velocity], None, Formulation::BurtonMiller, &limits(count), || false,
        ).unwrap().remove(0);
        let area: f64 = surface.areas().iter().sum();
        let average = response.pressure_derivative_pa_per_rad_s.iter().zip(surface.areas())
            .fold(C64::ZERO, |sum, (&p, &a)| sum + p.scale(a / area));
        errors.push((average - expected).abs() / expected.abs());
        close(response.wavenumber_derivative_s_m, dk, 5e-15);
        assert!(response.response.boundary_work_w() > 0.0);
    }
    assert!(errors[1] < errors[0] * 0.8 && errors[2] < errors[1] * 0.8, "{errors:?}");
    assert!(errors[2] < 0.06, "{errors:?}");
}

#[test]
fn one_factor_batch_preserves_superposition_and_a_zero_crossing_drive_slope() {
    let material = model(true);
    let surface = SpherePanels::icosphere(0.2, 0).unwrap();
    let n = surface.areas().len();
    let velocity = vec![C64::new(0.2, -0.05); n];
    let doubled: Vec<_> = velocity.iter().map(|v| v.scale(2.0)).collect();
    let zero = vec![C64::ZERO; n];
    let responses = material.solve_batch_with_frequency_derivative(
        &surface, 700.0, &[&velocity, &doubled, &zero], Some(&[&zero, &zero, &velocity]),
        Formulation::BurtonMiller, &limits(n), || false,
    ).unwrap();
    for i in 0..n {
        close(responses[1].pressure_derivative_pa_per_rad_s[i],
            responses[0].pressure_derivative_pa_per_rad_s[i].scale(2.0), 1e-13);
        close(responses[2].pressure_derivative_pa_per_rad_s[i],
            responses[0].response.pressure_pa()[i], 1e-13);
    }
    assert_eq!(responses[2].response.boundary_work_w(), 0.0);
    assert_eq!(responses[2].boundary_work_derivative_w_per_rad_s, 0.0);
    let fixed = material.solve_batch_with_frequency_derivative(
        &surface, 700.0, &[&velocity], None, Formulation::BurtonMiller, &limits(n), || false,
    ).unwrap();
    assert_eq!(fixed[0].pressure_derivative_pa_per_rad_s, responses[0].pressure_derivative_pa_per_rad_s);
}

#[test]
fn derivatives_admit_full_storage_and_complete_drives_and_never_publish_cancelled_batches() {
    let material = model(true);
    let surface = SpherePanels::icosphere(0.2, 0).unwrap();
    let n = surface.areas().len();
    let velocity = vec![C64::ONE; n];
    let formulation = Formulation::BurtonMiller;
    let mut budget = limits(n);
    budget.max_bytes -= 1;
    assert!(matches!(material.solve_batch_with_frequency_derivative(
        &surface, 700.0, &[&velocity], None, formulation, &budget, || false,
    ), Err(MolecularError::Operator(ComplexError::Budget { .. }))));
    for slopes in [Vec::<&[C64]>::new(), vec![&velocity[..n-1]]] {
        assert!(material.solve_batch_with_frequency_derivative(
            &surface, 700.0, &[&velocity], Some(&slopes), formulation, &limits(n), || false,
        ).is_err());
    }
    let mut invalid = velocity.clone();
    invalid[0].im = f64::NAN;
    assert!(material.solve_batch_with_frequency_derivative(
        &surface, 700.0, &[&velocity], Some(&[&invalid]), formulation, &limits(n), || false,
    ).is_err());
    assert!(material.solve_batch_with_frequency_derivative(
        &surface, 99.0, &[&velocity], None, formulation, &limits(n), || false,
    ).is_err());
    let mut polls = 0;
    material.solve_batch_with_frequency_derivative(
        &surface, 700.0, &[&velocity], None, formulation, &limits(n), || { polls += 1; false },
    ).unwrap();
    for stop in [1, polls / 2, polls] {
        let mut observed = 0;
        assert!(material.solve_batch_with_frequency_derivative(
            &surface, 700.0, &[&velocity], None, formulation, &limits(n),
            || { observed += 1; observed >= stop },
        ).is_err());
    }
}
