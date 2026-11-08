//! Exact-volume thermal work, physical rate closure and pressure-gauge invariance.

use fs_material::hyper::incompressible_plane_stress::IncompressiblePlaneStressPotential;
use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_material::hyper::reference::{ReferenceDeformation, ReferencedHyperelastic};
use fs_material::hyper::thermoelastic::ThermalReferenceJet;
use fs_material::hyper::thermoelastic::incompressible_plane_stress::IncompressibleThermoelasticPlaneStressState;
use fs_material::hyper::{Hyperelastic, HyperelasticModel};

const A: [f64; 4] = [1.22, 0.16, -0.05, 0.91];
const G: [f64; 9] = [1.09, 0.07, -0.025, 0.03, 0.97, 0.04, -0.015, 0.02, 1.04];
const G1: [f64; 9] = [
    0.002, 0.001, -0.0004, -0.0003, -0.001, 0.0006, 0.0002, 0.0008, 0.0015,
];
const G2: [f64; 9] = [2e-5, -3e-5, 1e-5, 1e-5, 4e-5, 0.0, -2e-5, 1e-5, -1e-5];

fn settings() -> PlaneStressSettings {
    PlaneStressSettings {
        relative_tolerance: 2e-12,
        ..Default::default()
    }
}
fn law() -> Hyperelastic {
    Hyperelastic::new(
        HyperelasticModel::MooneyRivlin {
            c10: 4.0,
            c01: 2.0,
            kappa: 40.0,
        },
        4.0,
    )
    .unwrap()
}
fn nh(lambda: f64) -> Hyperelastic {
    Hyperelastic::new(HyperelasticModel::NeoHookean { mu: 6.0, lambda }, 4.0).unwrap()
}
fn near(a: f64, b: f64, absolute: f64, relative: f64) {
    assert!(a.is_finite() && b.is_finite());
    let tolerance = absolute + relative * b.abs();
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}; tolerance {tolerance:e}"
    );
}
fn diagonal(value: f64) -> [f64; 9] {
    [value, 0.0, 0.0, 0.0, value, 0.0, 0.0, 0.0, value]
}
fn multiply(a: [f64; 9], b: [f64; 9]) -> [f64; 9] {
    core::array::from_fn(|ij| {
        (0..3)
            .map(|k| a[3 * (ij / 3) + k] * b[3 * k + ij % 3])
            .sum()
    })
}
fn cof(f: [f64; 9]) -> [f64; 9] {
    [
        f[4] * f[8] - f[5] * f[7],
        f[5] * f[6] - f[3] * f[8],
        f[3] * f[7] - f[4] * f[6],
        f[2] * f[7] - f[1] * f[8],
        f[0] * f[8] - f[2] * f[6],
        f[1] * f[6] - f[0] * f[7],
        f[1] * f[5] - f[2] * f[4],
        f[2] * f[3] - f[0] * f[5],
        f[0] * f[4] - f[1] * f[3],
    ]
}
fn jet(delta: f64) -> ThermalReferenceJet {
    ThermalReferenceJet::try_new(
        ReferenceDeformation::try_new(core::array::from_fn(|i| {
            G[i] + delta * G1[i] + 0.5 * delta * delta * G2[i]
        }))
        .unwrap(),
        core::array::from_fn(|i| G1[i] + delta * G2[i]),
        G2,
    )
    .unwrap()
}
fn state(
    material: &Hyperelastic,
    a: &[f64; 4],
    delta: f64,
) -> IncompressibleThermoelasticPlaneStressState {
    IncompressibleThermoelasticPlaneStressState::try_new(
        material,
        a,
        310.0 + delta,
        jet(delta),
        None,
        settings(),
    )
    .unwrap()
}
fn derivative(values: [f64; 4], h: f64) -> f64 {
    (values[0] - 8.0 * values[1] + 8.0 * values[2] - values[3]) / (12.0 * h)
}

#[test]
fn changing_natural_volume_matches_independent_scalar_energy_and_pressure() {
    let (mu, s, alpha, temperature) = (6.0, 1.1_f64, 0.002_f64, 310.0);
    let thermal = ThermalReferenceJet::try_new(
        ReferenceDeformation::try_new(diagonal(s)).unwrap(),
        diagonal(alpha * s),
        diagonal(alpha * alpha * s),
    )
    .unwrap();
    let d = A[0] * A[3] - A[1] * A[2];
    let da = [A[3], -A[2], -A[1], A[0]];
    let ia: f64 = A.iter().map(|x| x * x).sum();
    // Substituting F22=s^3/d into the NH energy gives this scalar
    // potential directly; no constrained derivative API enters the oracle.
    let psi = 0.5 * mu * (s * ia + s.powi(7) / (d * d) - 3.0 * s.powi(3));
    let entropy = -0.5 * mu * alpha * (s * ia + 7.0 * s.powi(7) / (d * d) - 9.0 * s.powi(3));
    let capacity = -0.5
        * temperature
        * mu
        * alpha
        * alpha
        * (s * ia + 49.0 * s.powi(7) / (d * d) - 27.0 * s.powi(3));
    let pressure = mu * (s.powi(4) / (d * d) - 1.0);
    let pressure_t = 4.0 * mu * alpha * s.powi(4) / (d * d);
    let mut previous = None;
    for bulk in [14.0, 1e250] {
        let state = IncompressibleThermoelasticPlaneStressState::try_new(
            &nh(bulk),
            &A,
            temperature,
            thermal,
            None,
            settings(),
        )
        .unwrap();
        near(state.free_energy_j_m3(), psi, 2e-12, 2e-12);
        near(state.entropy_j_m3_k(), entropy, 2e-13, 2e-11);
        near(
            state.internal_energy_j_m3(),
            psi + temperature * entropy,
            2e-11,
            2e-11,
        );
        near(
            state.mechanical_heat_capacity_j_m3_k(),
            capacity,
            2e-12,
            2e-11,
        );
        near(state.response().pressure_pa, pressure, 2e-12, 2e-12);
        near(
            state.pressure_temperature_derivative_pa_k(),
            pressure_t,
            2e-13,
            2e-11,
        );
        near(
            state.transverse_temperature_derivative_per_k()[2],
            3.0 * alpha * s.powi(3) / d,
            2e-13,
            2e-12,
        );
        near(
            state.natural_volume_derivative_per_k(),
            3.0 * alpha * s.powi(3),
            2e-14,
            2e-12,
        );
        near(
            state.natural_volume_second_derivative_per_k2(),
            9.0 * alpha * alpha * s.powi(3),
            2e-15,
            2e-12,
        );
        for i in 0..4 {
            let p = mu * (s * A[i] - s.powi(7) * da[i] / d.powi(3));
            let beta = mu * alpha * (s * A[i] - 7.0 * s.powi(7) * da[i] / d.powi(3));
            near(state.response().piola[i], p, 2e-12, 2e-12);
            near(
                state.stress_temperature_derivative_pa_k()[i],
                beta,
                2e-13,
                2e-11,
            );
        }
        let observed = (
            state.free_energy_j_m3(),
            state.entropy_j_m3_k(),
            state.mechanical_heat_capacity_j_m3_k(),
            state.response().full_piola,
            *state.stress_temperature_derivative_pa_k(),
        );
        if let Some(previous) = previous {
            assert_eq!(observed, previous);
        }
        previous = Some(observed);
    }
    assert!(pressure.abs() > 0.1);
    assert!(capacity < 0.0); // The mechanical contribution is not clipped.
}

#[test]
fn noncommuting_thermal_path_closes_volume_tractions_pressure_and_energy_rates() {
    let material = law();
    let center = state(&material, &A, 0.0);
    let adot = [0.08, -0.03, 0.01, -0.05];
    let tdot = 2.3;
    let rate = center.rate(adot, tdot).unwrap();
    let h = 0.005;
    let samples = [-2.0 * h, -h, h, 2.0 * h].map(|time| {
        state(
            &material,
            &core::array::from_fn(|i| A[i] + time * adot[i]),
            tdot * time,
        )
    });
    for i in 0..9 {
        near(
            rate.material().deformation_gradient_rate_per_s()[i],
            derivative(
                samples
                    .each_ref()
                    .map(|s| s.response().deformation_gradient[i]),
                h,
            ),
            4e-8,
            3e-7,
        );
        near(
            rate.material().elastic_deformation_gradient_rate_per_s()[i],
            derivative(
                samples
                    .each_ref()
                    .map(|s| s.elastic_deformation_gradient()[i]),
                h,
            ),
            4e-8,
            3e-7,
        );
        near(
            rate.material().first_piola_rate_pa_s()[i],
            derivative(samples.each_ref().map(|s| s.response().full_piola[i]), h),
            3e-7,
            3e-7,
        );
    }
    for (actual, values) in [
        (
            rate.pressure_rate_pa_s(),
            samples.each_ref().map(|s| s.response().pressure_pa),
        ),
        (
            rate.material().free_energy_rate_w_m3(),
            samples.each_ref().map(|s| s.free_energy_j_m3()),
        ),
        (
            rate.material().entropy_rate_w_m3_k(),
            samples.each_ref().map(|s| s.entropy_j_m3_k()),
        ),
        (
            rate.material().internal_energy_rate_w_m3(),
            samples.each_ref().map(|s| s.internal_energy_j_m3()),
        ),
    ] {
        near(actual, derivative(values, h), 3e-7, 3e-7);
    }
    let fdot = rate.material().deformation_gradient_rate_per_s();
    let fedot = rate.material().elastic_deformation_gradient_rate_per_s();
    near(
        (0..9)
            .map(|i| cof(center.response().deformation_gradient)[i] * fdot[i])
            .sum(),
        rate.natural_volume_ratio_rate_per_s(),
        2e-13,
        2e-12,
    );
    near(
        (0..9)
            .map(|i| cof(*center.elastic_deformation_gradient())[i] * fedot[i])
            .sum(),
        0.0,
        2e-13,
        0.0,
    );
    for i in [2, 5, 8] {
        near(rate.material().first_piola_rate_pa_s()[i], 0.0, 3e-12, 0.0);
    }
    // Pressure work is a real difference between remainder and physical power.
    let remainder = ReferencedHyperelastic::new(&material, *jet(0.0).reference())
        .incompressible_remainder_response(&center.response().deformation_gradient)
        .unwrap();
    let remainder_power: f64 = (0..9).map(|i| remainder.first_piola_pa[i] * fdot[i]).sum();
    let physical_power: f64 = (0..9)
        .map(|i| center.response().full_piola[i] * fdot[i])
        .sum();
    near(
        remainder_power - physical_power,
        center.response().pressure_pa * rate.natural_volume_ratio_rate_per_s(),
        3e-12,
        3e-12,
    );
    assert!((remainder_power - physical_power).abs() > 1e-5);
    assert!(center.transverse_temperature_derivative_per_k()[0].abs() > 1e-6);
    assert!(center.transverse_temperature_derivative_per_k()[1].abs() > 1e-6);
}

#[test]
fn constrained_thermodynamics_are_independent_of_native_remainder_pressure_gauge() {
    let original = nh(14.0);
    let isochoric = Hyperelastic::new(
        HyperelasticModel::MooneyRivlin {
            c10: 3.0,
            c01: 0.0,
            kappa: 1e250,
        },
        4.0,
    )
    .unwrap();
    let first = state(&original, &A, 0.0);
    let second = state(&isochoric, &A, 0.0);
    assert!((first.response().pressure_pa - second.response().pressure_pa).abs() > 1e-3);
    for (a, b) in [
        (first.free_energy_j_m3(), second.free_energy_j_m3()),
        (first.entropy_j_m3_k(), second.entropy_j_m3_k()),
        (first.internal_energy_j_m3(), second.internal_energy_j_m3()),
        (
            first.mechanical_heat_capacity_j_m3_k(),
            second.mechanical_heat_capacity_j_m3_k(),
        ),
    ] {
        near(a, b, 2e-10, 2e-10);
    }
    for (a, b) in first
        .response()
        .full_piola
        .iter()
        .zip(second.response().full_piola)
    {
        near(*a, b, 2e-10, 2e-10);
    }
    for (a, b) in first
        .stress_temperature_derivative_pa_k()
        .iter()
        .zip(second.stress_temperature_derivative_pa_k())
    {
        near(*a, *b, 2e-11, 2e-10);
    }
    let first_rate = first.rate([0.08, -0.03, 0.01, -0.05], 2.3).unwrap();
    let second_rate = second.rate([0.08, -0.03, 0.01, -0.05], 2.3).unwrap();
    for (a, b) in first_rate
        .material()
        .first_piola_rate_pa_s()
        .iter()
        .zip(second_rate.material().first_piola_rate_pa_s())
    {
        near(*a, *b, 2e-10, 2e-10);
    }
    for (a, b) in [
        (
            first_rate.material().free_energy_rate_w_m3(),
            second_rate.material().free_energy_rate_w_m3(),
        ),
        (
            first_rate.material().entropy_rate_w_m3_k(),
            second_rate.material().entropy_rate_w_m3_k(),
        ),
        (
            first_rate.material().internal_energy_rate_w_m3(),
            second_rate.material().internal_energy_rate_w_m3(),
        ),
    ] {
        near(a, b, 2e-10, 2e-10);
    }
}

#[test]
fn physical_rates_transform_with_a_changing_planar_observer() {
    let material = law();
    let center = state(&material, &A, 0.0);
    let adot = [0.08, -0.03, 0.01, -0.05];
    let rate = center.rate(adot, -2.1).unwrap();
    let q = [0.6, -0.8, 0.0, 0.8, 0.6, 0.0, 0.0, 0.0, 1.0];
    let spin = [0.0, -0.7, 0.0, 0.7, 0.0, 0.0, 0.0, 0.0, 0.0];
    let f = multiply(q, center.response().deformation_gradient);
    let q_fdot = multiply(q, *rate.material().deformation_gradient_rate_per_s());
    let spin_f = multiply(spin, f);
    let expected_fdot: [f64; 9] = core::array::from_fn(|i| q_fdot[i] + spin_f[i]);
    let rotated = IncompressibleThermoelasticPlaneStressState::try_new(
        &material,
        &[f[0], f[1], f[3], f[4]],
        310.0,
        jet(0.0),
        Some([f[2], f[5]]),
        settings(),
    )
    .unwrap();
    let observed = rotated
        .rate(
            [
                expected_fdot[0],
                expected_fdot[1],
                expected_fdot[3],
                expected_fdot[4],
            ],
            -2.1,
        )
        .unwrap();
    let q_pdot = multiply(q, *rate.material().first_piola_rate_pa_s());
    let spin_p = multiply(spin, multiply(q, center.response().full_piola));
    for i in 0..9 {
        near(
            observed.material().deformation_gradient_rate_per_s()[i],
            expected_fdot[i],
            2e-12,
            2e-11,
        );
        near(
            observed.material().first_piola_rate_pa_s()[i],
            q_pdot[i] + spin_p[i],
            2e-10,
            2e-10,
        );
    }
    near(
        observed.pressure_rate_pa_s(),
        rate.pressure_rate_pa_s(),
        2e-11,
        2e-10,
    );
    for (a, b) in [
        (
            observed.material().free_energy_rate_w_m3(),
            rate.material().free_energy_rate_w_m3(),
        ),
        (
            observed.material().entropy_rate_w_m3_k(),
            rate.material().entropy_rate_w_m3_k(),
        ),
        (
            observed.material().internal_energy_rate_w_m3(),
            rate.material().internal_energy_rate_w_m3(),
        ),
    ] {
        near(a, b, 2e-10, 2e-10);
    }
}

#[test]
fn free_expansion_has_zero_elastic_rate_and_invalid_queries_leave_state_unchanged() {
    let (s, alpha, tdot) = (1.1, 0.002, 3.0);
    let material = nh(14.0);
    let thermal = ThermalReferenceJet::try_new(
        ReferenceDeformation::try_new(diagonal(s)).unwrap(),
        diagonal(alpha * s),
        diagonal(alpha * alpha * s),
    )
    .unwrap();
    let center = IncompressibleThermoelasticPlaneStressState::try_new(
        &material,
        &[s, 0.0, 0.0, s],
        310.0,
        thermal,
        None,
        settings(),
    )
    .unwrap();
    let before = center.response().deformation_gradient;
    let rate = center
        .rate([alpha * s * tdot, 0.0, 0.0, alpha * s * tdot], tdot)
        .unwrap();
    near(rate.pressure_rate_pa_s(), 0.0, 2e-12, 0.0);
    for (actual, expected) in rate
        .material()
        .deformation_gradient_rate_per_s()
        .iter()
        .zip(diagonal(alpha * s * tdot))
    {
        near(*actual, expected, 2e-12, 2e-12);
    }
    for actual in rate
        .material()
        .elastic_deformation_gradient_rate_per_s()
        .iter()
        .chain(rate.material().first_piola_rate_pa_s())
    {
        near(*actual, 0.0, 2e-12, 0.0);
    }
    for actual in [
        rate.material().free_energy_rate_w_m3(),
        rate.material().entropy_rate_w_m3_k(),
        rate.material().internal_energy_rate_w_m3(),
    ] {
        near(actual, 0.0, 2e-11, 0.0);
    }
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(center.rate([bad; 4], 0.0).is_err());
        assert!(center.rate([0.0; 4], bad).is_err());
    }
    for temperature in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            IncompressibleThermoelasticPlaneStressState::try_new(
                &material,
                &A,
                temperature,
                thermal,
                None,
                settings(),
            )
            .is_err()
        );
    }
    for a in [[1.0, 2.0, 2.0, 1.0], [10.0, 0.0, 0.0, 10.0], [f64::NAN; 4]] {
        assert!(
            IncompressibleThermoelasticPlaneStressState::try_new(
                &material,
                &a,
                310.0,
                thermal,
                None,
                settings(),
            )
            .is_err()
        );
    }
    assert_eq!(center.response().deformation_gradient, before);
}

#[test]
fn native_exact_thermal_increments_retain_changes_below_absolute_storage_resolution() {
    let material = nh(1e250);
    let constant =
        ThermalReferenceJet::try_new(ReferenceDeformation::identity(), [0.0; 9], [0.0; 9]).unwrap();
    let before = IncompressibleThermoelasticPlaneStressState::try_new(
        &material,
        &[2.0, 0.0, 0.0, 2.0],
        300.0,
        constant,
        None,
        settings(),
    )
    .unwrap();
    let shear = 1e-8;
    let after = IncompressibleThermoelasticPlaneStressState::try_new(
        &material,
        &[2.0, shear, 0.0, 2.0],
        300.0,
        constant,
        None,
        settings(),
    )
    .unwrap();
    // The area and exact thickness are unchanged. The NH remainder has the
    // independent increment mu*shear^2/2 despite a much larger prestress energy.
    let expected = 3.0 * shear * shear;
    assert_eq!(after.free_energy_j_m3() - before.free_energy_j_m3(), 0.0);
    near(
        before.free_energy_increment_to(&after).unwrap(),
        expected,
        1e-29,
        1e-12,
    );
    near(
        before.internal_energy_increment_to(&after).unwrap(),
        expected,
        1e-29,
        1e-12,
    );
    near(
        after.free_energy_increment_to(&before).unwrap(),
        -expected,
        1e-29,
        1e-12,
    );
    assert_eq!(before.free_energy_increment_to(&before).unwrap(), 0.0);
    assert_eq!(before.internal_energy_increment_to(&before).unwrap(), 0.0);
}
