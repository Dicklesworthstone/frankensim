//! Complete physical direction covariance under a change of displacement units.
//! Original power ports and material strains stay fixed; every state survives.

use fs_math::{c64::C64, det};
use fs_phs::RelaxationBranch;
use fs_phs::harmonic::{HarmonicBudget, HarmonicFrequencyParameterResponse, HarmonicModel};
use fs_phs::mechanical::{
    LinearMechanicalSystem, MechanicalHarmonicDirection, MechanicalHarmonicDirectionControls,
    MechanicalRelaxationDirection, PreparedMechanicalHarmonicDirection,
};

fn near(actual: C64, expected: C64, name: &str, scale: f64, order: usize) {
    assert!(
        (actual - expected).abs() < 3e-11 * expected.abs().max(1.0),
        "{name}, scale={scale:e}, order={order}: {actual:?} != {expected:?}"
    );
}

fn physical_model(scale: f64) -> (HarmonicModel, PreparedMechanicalHarmonicDirection) {
    // q_new=a*q, p_new=p/a, M_new=M/a², P_new=P/a and b_new=b/a.
    // This is a fixed chart change, not an extra parameter-dependent basis.
    let inverse = 1.0 / scale;
    let square = inverse * inverse;
    let mass = [2.0 * square];
    let damping = [0.125 * square];
    let stiffness = [3.0 * square];
    let force_map = [0.75 * inverse];
    let system = LinearMechanicalSystem::new(
        &mass,
        &damping,
        &stiffness,
        &force_map,
        1,
        1,
        vec![RelaxationBranch {
            projection: vec![0.5 * inverse],
            stiffness: 1.25,
            relaxation_time_s: 0.875,
        }],
    )
    .unwrap();
    let d_mass = [0.25 * square];
    let d_damping = [-0.03125 * square];
    let d_stiffness = [0.5 * square];
    let d_force_map = [0.125 * inverse];
    let d_projection = [-0.0625 * inverse];
    let d_branch = [MechanicalRelaxationDirection {
        projection: Some(&d_projection),
        stiffness: 0.375,
        relaxation_time_s: 0.125,
    }];
    let direction = system
        .prepare_harmonic_parameter_direction(
            MechanicalHarmonicDirection {
                mass: Some(&d_mass),
                damping: Some(&d_damping),
                relaxed_stiffness: Some(&d_stiffness),
                force_map: Some(&d_force_map),
                branches: Some(&d_branch),
            },
            MechanicalHarmonicDirectionControls {
                max_state_dimension: 16,
                max_dense_bytes: 1_000_000,
                max_multiply_terms: 1_000_000,
            },
            || false,
        )
        .unwrap();
    assert_eq!(direction.coordinates(), 1);
    assert_eq!(direction.state_dim(), 3);
    assert_eq!(direction.port_dim(), 1);
    near(
        C64::from_re(direction.branch_root_derivatives()[0]),
        C64::from_re(0.375 / (2.0 * det::sqrt(1.25))),
        "physical branch root direction",
        scale,
        0,
    );
    near(
        C64::from_re(direction.branch_rate_derivatives()[0]),
        C64::from_re(-0.125 / (0.875 * 0.875)),
        "physical branch rate direction",
        scale,
        0,
    );
    let model = HarmonicModel::from_ode(
        system.system(),
        HarmonicBudget {
            max_dimension: 16,
            max_dense_bytes: 2_000_000,
            max_dense_work: 8_000_000,
            relative_residual_tolerance: 1e-12,
        },
        || false,
    )
    .unwrap();
    (model, direction)
}

fn solve(
    model: &HarmonicModel,
    direction: &PreparedMechanicalHarmonicDirection,
    omega: f64,
) -> HarmonicFrequencyParameterResponse {
    let delay = 0.2;
    let delay_parameter = 0.05;
    let phase = C64::new(det::cos(omega * delay), det::sin(omega * delay));
    let input = C64::new(0.75, -0.25) * phase;
    let frequency = C64::new(0.0, delay) * input;
    let second = input.scale(-delay * delay);
    let parameter =
        C64::new(0.125, 0.0625) * phase + C64::new(0.0, omega * delay_parameter) * input;
    let mixed = C64::new(0.0, delay) * parameter + C64::new(0.0, delay_parameter) * input;
    model
        .factor(omega, || false)
        .unwrap()
        .response_with_frequency_parameter_derivative(
            &[input],
            direction.as_harmonic_direction(),
            Some(&[frequency]),
            Some(&[second]),
            Some(&[parameter]),
            Some(&[mixed]),
            || false,
        )
        .unwrap()
}

fn fields(response: &HarmonicFrequencyParameterResponse) -> [[&[C64]; 5]; 3] {
    let second = &response.frequency_response;
    let first = &second.frequency_response;
    let base = &first.response;
    let parameter = &response.parameter_response;
    [
        [
            &base.state,
            &first.state_derivative,
            &second.state_second_derivative,
            &parameter.state_derivative,
            &response.state_mixed_derivative,
        ],
        [
            &base.effort,
            &first.effort_derivative,
            &second.effort_second_derivative,
            &parameter.effort_derivative,
            &response.effort_mixed_derivative,
        ],
        [
            &base.output,
            &first.output_derivative,
            &second.output_second_derivative,
            &parameter.output_derivative,
            &response.output_mixed_derivative,
        ],
    ]
}

fn powers(response: &HarmonicFrequencyParameterResponse) -> [[f64; 5]; 4] {
    let second = &response.frequency_response;
    let first = &second.frequency_response;
    let base = &first.response;
    let parameter = &response.parameter_response;
    [
        [
            base.mean_stored_energy_j,
            first.mean_stored_energy_derivative,
            second.mean_stored_energy_second_derivative,
            parameter.mean_stored_energy_derivative,
            response.mean_stored_energy_mixed_derivative,
        ],
        [
            base.supplied_power_w,
            first.supplied_power_derivative,
            second.supplied_power_second_derivative,
            parameter.supplied_power_derivative,
            response.supplied_power_mixed_derivative,
        ],
        [
            base.dissipated_power_w,
            first.dissipated_power_derivative,
            second.dissipated_power_second_derivative,
            parameter.dissipated_power_derivative,
            response.dissipated_power_mixed_derivative,
        ],
        [
            base.port_power_w[0],
            first.port_power_derivative[0],
            second.port_power_second_derivative[0],
            parameter.port_power_derivative[0],
            response.port_power_mixed_derivative[0],
        ],
    ]
}

#[test]
fn full_mechanical_and_maxwell_directions_preserve_physical_jets_across_coordinate_units() {
    let frequencies = [0.4, 1.3, 2.5];
    let (reference_model, reference_direction) = physical_model(1.0);
    let reference: Vec<_> = frequencies
        .iter()
        .map(|&omega| solve(&reference_model, &reference_direction, omega))
        .collect();
    for scale in [2.0_f64.powi(-400), 2.0_f64.powi(400)] {
        let (model, direction) = physical_model(scale);
        for (&omega, reference) in frequencies.iter().zip(&reference) {
            let response = solve(&model, &direction, omega);
            let actual_fields = fields(&response);
            let reference_fields = fields(reference);
            for order in 0..5 {
                // Normalized memory strain has unchanged physical meaning.
                for (coordinate, conversion) in [1.0 / scale, scale, 1.0].into_iter().enumerate() {
                    near(
                        actual_fields[0][order][coordinate].scale(conversion),
                        reference_fields[0][order][coordinate],
                        "original state",
                        scale,
                        order,
                    );
                }
                for (coordinate, conversion) in [scale, 1.0 / scale, 1.0].into_iter().enumerate() {
                    near(
                        actual_fields[1][order][coordinate].scale(conversion),
                        reference_fields[1][order][coordinate],
                        "power-conjugate effort",
                        scale,
                        order,
                    );
                }
                near(
                    actual_fields[2][order][0],
                    reference_fields[2][order][0],
                    "original port output",
                    scale,
                    order,
                );
            }
            for (actual, expected) in powers(&response).into_iter().zip(powers(reference)) {
                for order in 0..5 {
                    near(
                        C64::from_re(actual[order]),
                        C64::from_re(expected[order]),
                        "original energy or work",
                        scale,
                        order,
                    );
                }
            }
            assert!(response.relative_mixed_derivative_residual < 1e-12);
            assert!(
                response
                    .frequency_response
                    .relative_second_derivative_residual
                    < 1e-12
            );
            // Frequency and material work are independent observations, not
            // inferred heat. Their signed balances remain physical after scaling.
            near(
                C64::from_re(response.power_residual_mixed_derivative),
                C64::ZERO,
                "mixed power balance",
                scale,
                4,
            );
        }
    }
}
