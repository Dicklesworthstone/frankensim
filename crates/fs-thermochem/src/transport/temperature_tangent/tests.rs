use super::*;
use crate::transport::{WilkeGasViscosityComponentV1, WilkeGasViscosityMixtureV1};
use crate::{Composition, CompositionBasis, SpeciesId};
use fs_qty::{DynViscosity, MolarMass, Temperature};

fn evaluation(
    viscosities: &[f64],
    masses: &[f64],
    fractions: &[f64],
) -> WilkeGasViscosityEvaluationV1 {
    let components = viscosities
        .iter()
        .zip(masses)
        .enumerate()
        .map(|(i, (&mu, &mass))| {
            WilkeGasViscosityComponentV1::new(
                SpeciesId::new(format!("S{i}")).unwrap(),
                MolarMass::new(mass),
                DynViscosity::new(mu),
            )
            .unwrap()
        })
        .collect();
    WilkeGasViscosityMixtureV1::new(
        Temperature::new(300.0),
        components,
        Composition::new(CompositionBasis::MoleFraction, fractions.to_vec()).unwrap(),
    )
    .unwrap()
    .evaluate()
    .unwrap()
}

fn near(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 2e-7 * a.abs().max(b.abs()),
        "{a:e} != {b:e}"
    );
}

#[test]
fn full_wilke_temperature_response_matches_independent_perturbed_original_values() {
    let viscosity = [1.7e-5, 2.0e-5, 3.2e-5];
    let derivatives = [4e-8, -1e-8, 2e-8];
    let masses = [0.028, 0.004, 0.044];
    let fractions = [0.5, 0.25, 0.25];
    let original = evaluation(&viscosity, &masses, &fractions);
    let before = original.clone();
    let tangent = original.temperature_tangent(&derivatives).unwrap();
    let h = 0.01;
    let values = [-2.0, -1.0, 1.0, 2.0].map(|factor| {
        evaluation(
            &[0, 1, 2].map(|i| viscosity[i] + factor * h * derivatives[i]),
            &masses,
            &fractions,
        )
    });
    let fd = |v: [f64; 4]| (v[0] - 8.0 * v[1] + 8.0 * v[2] - v[3]) / (12.0 * h);
    near(
        tangent.dynamic_viscosity_derivative_pa_s_k,
        fd(values
            .each_ref()
            .map(|value| value.dynamic_viscosity().value())),
    );
    assert!(
        tangent
            .denominator_derivatives_per_k()
            .iter()
            .any(|value| *value != 0.0)
    );
    for i in 0..3 {
        near(
            tangent.denominator_derivatives_per_k()[i],
            fd(values
                .each_ref()
                .map(|value| f64::from_bits(value.receipt().denominator_bits()[i]))),
        );
    }
    let stationary = original.temperature_tangent(&[0.0; 3]).unwrap();
    assert_eq!(stationary.dynamic_viscosity_derivative_pa_s_k, 0.0);
    assert_eq!(stationary.denominator_derivatives_per_k(), &[0.0; 3]);
    assert_eq!(original, before);
}

#[test]
fn exact_pure_and_common_scaling_limits_retain_signed_slopes_and_refuse_malformed_inputs() {
    let pure = evaluation(&[1e-5], &[0.028], &[1.0]);
    for derivative in [0.0, -1e-8, 3e-8, f64::MAX] {
        let response = pure.temperature_tangent(&[derivative]).unwrap();
        assert_eq!(response.dynamic_viscosity_derivative_pa_s_k, derivative);
        assert_eq!(response.denominator_derivatives_per_k(), &[0.0]);
    }
    let mixed = evaluation(&[2e-5, 4e-5], &[0.028, 0.004], &[0.5, 0.5]);
    let response = mixed.temperature_tangent(&[2e-8, 4e-8]).unwrap();
    assert_eq!(response.denominator_derivatives_per_k(), &[0.0; 2]);
    near(
        response.dynamic_viscosity_derivative_pa_s_k,
        mixed.dynamic_viscosity().value() * 0.001,
    );
    assert!(matches!(
        mixed.temperature_tangent(&[1e-8]),
        Err(WilkeTemperatureTangentErrorV1::Length { .. })
    ));
    for value in [f64::NAN, f64::INFINITY] {
        assert!(matches!(
            mixed.temperature_tangent(&[value, 0.0]),
            Err(WilkeTemperatureTangentErrorV1::Arithmetic { .. })
        ));
    }
}
