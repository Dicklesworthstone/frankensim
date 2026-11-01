use super::*;
use fs_bem::observation_transfers::{PreparedRadiationOptions, RadiatedPowerMoments};

const DIRECTIONS: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]];
const ANGULAR_BUDGET: HarmonicWorkBudget = HarmonicWorkBudget {
    maximum_products: 1_000_000,
};

fn options() -> PreparedRadiationOptions {
    PreparedRadiationOptions {
        angular: helmholtz::DirectivityProjectionOptions::default(),
        initial_latitudes: 8,
        refinements: 3,
        maximum_kernel_evaluations: 4_000_000,
        maximum_source_entries: 4096,
    }
}

fn mixed_sources() -> [PortSource; 2] {
    [
        PortSource {
            label: "inlet".into(),
            excitation: PortExcitation::VolumeFlow(C64::new(1e-6, -0.4e-6)),
        },
        PortSource {
            label: "aux".into(),
            excitation: PortExcitation::Pressure {
                open_circuit_pressure_pa: C64::new(0.6, 0.2),
                impedance_pa_s_m3: C64::new(2e6, -1e6),
            },
        },
    ]
}

fn asymmetric_model() -> RadiatingNetwork {
    let mut model = two_source_model();
    // Attach the second source at the left opening: its angular pattern no
    // longer equals the inlet's, so phase changes affect directional moments.
    model.driving.network.edges.last_mut().unwrap().to = 2;
    model
}

fn direct_solution(
    model: &RadiatingNetwork,
    omega: f64,
    flows: &[C64],
) -> helmholtz::RadiationSolution {
    let mut velocity = vec![C64::ZERO; model.surface().areas().len()];
    for ((patch, &area), &flow) in model
        .patches()
        .iter()
        .zip(model.opening_areas_m2())
        .zip(flows)
    {
        for &panel in &patch.panels {
            velocity[panel] = flow.scale(area.recip());
        }
    }
    let medium = Medium::from(&model.gas());
    helmholtz::solve_radiation(
        model.surface(),
        omega / medium.sound_speed,
        medium,
        &velocity,
        Formulation::BurtonMiller,
    )
    .unwrap()
}

fn direct_far_field(
    model: &RadiatingNetwork,
    solution: &helmholtz::RadiationSolution,
    directions: &[[f64; 3]],
) -> Vec<C64> {
    helmholtz::far_field_integrated_cancellable(
        model.surface(),
        solution,
        Medium::from(&model.gas()),
        directions,
        10_000_000,
        || false,
    )
    .unwrap()
}

// Independent composite Simpson integration in mu=cos(theta), with uniform
// azimuths, not the cached Gauss-Legendre ladder. The positive 65-node mu rule
// resolves a pure dipole's fourth-order directional moment to 3.82e-7 of its
// power; a 64-slice midpoint rule errs by 4.88e-4. A newly assembled physical
// velocity BEM solve separately checks the source and field composition.
fn simpson_moments(
    model: &RadiatingNetwork,
    solution: &helmholtz::RadiationSolution,
) -> RadiatedPowerMoments {
    let mut directions = Vec::new();
    for latitude in 0..=64 {
        let z = -1.0 + latitude as f64 / 32.0;
        let radius = (1.0 - z * z).sqrt();
        for longitude in 0..64 {
            let phase = core::f64::consts::TAU * (longitude as f64 + 0.5) / 64.0;
            directions.push([radius * phase.cos(), radius * phase.sin(), z]);
        }
    }
    let amplitudes = direct_far_field(model, solution, &directions);
    let weight = core::f64::consts::TAU
        / (32.0 * 3.0 * 64.0)
        / (2.0 * model.gas().density * model.gas().sound_speed);
    let mut result = RadiatedPowerMoments::default();
    for (index, (direction, amplitude)) in directions.iter().zip(amplitudes).enumerate() {
        let latitude = index / 64;
        let simpson_weight = if latitude == 0 || latitude == 64 {
            1.0
        } else if latitude % 2 == 1 {
            4.0
        } else {
            2.0
        };
        let power = weight * simpson_weight * amplitude.norm_sq();
        result.total_power_w += power;
        for i in 0..3 {
            result.first_directional_moment_w[i] += power * direction[i];
            for j in 0..3 {
                result.second_directional_moment_w[i][j] += power * direction[i] * direction[j];
            }
        }
    }
    result
}

fn moments_close(actual: RadiatedPowerMoments, expected: RadiatedPowerMoments, tolerance: f64) {
    let close = |a: f64, b: f64| {
        assert!(
            (a - b).abs() <= tolerance * expected.total_power_w.max(1e-30),
            "{a:e} != {b:e}"
        )
    };
    close(actual.total_power_w, expected.total_power_w);
    for i in 0..3 {
        close(
            actual.first_directional_moment_w[i],
            expected.first_directional_moment_w[i],
        );
        for j in 0..3 {
            close(
                actual.second_directional_moment_w[i][j],
                expected.second_directional_moment_w[i][j],
            );
        }
    }
}

#[test]
fn g1_loaded_mixed_units_reach_independent_far_field_and_positive_power_integral() {
    let model = asymmetric_model();
    let sources = mixed_sources();
    let omega = core::f64::consts::TAU * 1200.0;
    let gate = CancelGate::new_clock_free();
    let plain = model
        .prepare_sources(omega, LossModel::Bessel, &sources, None, &[], &gate)
        .unwrap();
    assert!(plain.far_field().is_none());
    let prepared = model
        .prepare_sources_with_radiation(
            omega,
            LossModel::Bessel,
            &sources,
            None,
            &[],
            &DIRECTIONS,
            options(),
            &gate,
        )
        .unwrap();
    assert_eq!(prepared.boundary().z, plain.boundary().z);
    assert!(prepared.boundary().z[1].abs() > 0.0);
    let far = prepared.far_field().unwrap();
    assert_eq!(far.source_coordinates(), prepared.source_coordinates());
    assert_eq!(
        far.source_coordinates()[0].quantity,
        NetworkSourceQuantity::OpenCircuitPressure
    );
    assert_eq!(
        far.source_coordinates()[1].quantity,
        NetworkSourceQuantity::VolumeFlow
    );
    assert_eq!(far.directions(), DIRECTIONS);
    assert_eq!(far.cached_angular_directions(), 136 + 351 + 820);
    let kernel_work = far.kernel_evaluations();
    let mut response = prepared.response_buffer();
    let mut plain_response = plain.response_buffer();
    let drive = prepared.default_excitations();
    prepared.drive_into(drive, &mut response, &gate).unwrap();
    plain.drive_into(drive, &mut plain_response, &gate).unwrap();
    assert_eq!(response.receiver_fields(), plain_response.receiver_fields());
    assert_eq!(
        response.exterior_surface_power_w(),
        plain_response.exterior_surface_power_w()
    );
    let solution = direct_solution(&model, omega, response.outward_flows_m3_s().unwrap());
    let expected = direct_far_field(&model, &solution, &DIRECTIONS);
    let mut actual = [C64::ZERO; 3];
    let power = far
        .evaluate_coherent(drive, ANGULAR_BUDGET, &mut actual, &gate)
        .unwrap();
    for (&a, &b) in actual.iter().zip(&expected) {
        close(a, b);
    }
    assert!(power.total_power_w > 0.0);
    moments_close(power, simpson_moments(&model, &solution), 3e-4);
    real_close(
        power
            .second_directional_moment_w
            .iter()
            .enumerate()
            .map(|(i, row)| row[i])
            .sum(),
        power.total_power_w,
    );

    // A common complex phase and amplitude change must propagate in physical SI.
    let phase = C64::new(0.6, 0.8).scale(3.0);
    let scaled: Vec<_> = drive.iter().map(|&x| x * phase).collect();
    let scaled_power = far
        .evaluate_coherent(&scaled, ANGULAR_BUDGET, &mut actual, &gate)
        .unwrap();
    for (&a, &b) in actual.iter().zip(&expected) {
        close(a, b * phase);
    }
    real_close(scaled_power.total_power_w, 9.0 * power.total_power_w);
    assert_eq!(far.kernel_evaluations(), kernel_work);
}

#[test]
fn g3_dark_loaded_generators_and_partial_coherence_keep_radiation_and_internal_work_distinct() {
    let model = two_source_model();
    let omega = core::f64::consts::TAU * 1200.0;
    let bank = bank(omega);
    let gate = CancelGate::new_clock_free();
    let prepared = model
        .prepare_sources_with_radiation(
            omega,
            LossModel::Lossless,
            &[],
            Some(&bank),
            &[],
            &DIRECTIONS,
            options(),
            &gate,
        )
        .unwrap();
    let far = prepared.far_field().unwrap();
    let mut bright = [C64::ZERO; 3];
    let bright_power = far
        .evaluate_coherent(&[C64::ONE; 2], ANGULAR_BUDGET, &mut bright, &gate)
        .unwrap();
    let mut dark = [C64::ZERO; 3];
    let dark_power = far
        .evaluate_coherent(&[C64::ONE, -C64::ONE], ANGULAR_BUDGET, &mut dark, &gate)
        .unwrap();
    assert!(dark_power.total_power_w < 1e-18 * bright_power.total_power_w);
    let mut response = prepared.response_buffer();
    prepared
        .drive_into(&[C64::ONE, -C64::ONE], &mut response, &gate)
        .unwrap();
    assert!(response.network().unwrap().source_dissipation_w > 0.0);
    let mut independent = [0.0; 3];
    let independent_power = far
        .evaluate_independent_sources(&[1.0; 2], ANGULAR_BUDGET, &mut independent, &gate)
        .unwrap();
    real_close(
        independent_power.total_power_w,
        0.5 * (bright_power.total_power_w + dark_power.total_power_w),
    );
    let mut factors = [0.0; 3];
    let scale = core::f64::consts::FRAC_1_SQRT_2;
    let rotated = [
        C64::from_re(scale),
        C64::from_re(scale),
        C64::from_re(scale),
        C64::from_re(-scale),
    ];
    let factor_power = far
        .evaluate_coherence_modes(&rotated, 2, ANGULAR_BUDGET, &mut factors, &gate)
        .unwrap();
    moments_close(factor_power, independent_power, 3e-9);
    for (&a, &b) in factors.iter().zip(&independent) {
        real_close(a, b);
    }
}

#[test]
fn g3_factor_phases_and_independent_generators_retain_correlated_opening_fields() {
    let model = asymmetric_model();
    let sources = mixed_sources();
    let omega = core::f64::consts::TAU * 1200.0;
    let gate = CancelGate::new_clock_free();
    let prepared = model
        .prepare_sources_with_radiation(
            omega,
            LossModel::Bessel,
            &sources,
            None,
            &[],
            &DIRECTIONS,
            options(),
            &gate,
        )
        .unwrap();
    let far = prepared.far_field().unwrap();
    let factors = [
        C64::new(0.7, -0.2),
        C64::new(1e-6, 0.5e-6),
        C64::new(-0.1, 0.4),
        C64::new(-0.3e-6, 0.2e-6),
    ];
    let mut first = [C64::ZERO; 3];
    let mut second = [C64::ZERO; 3];
    let a = far
        .evaluate_coherent(&factors[..2], ANGULAR_BUDGET, &mut first, &gate)
        .unwrap();
    let b = far
        .evaluate_coherent(&factors[2..], ANGULAR_BUDGET, &mut second, &gate)
        .unwrap();
    let mut density = [0.0; 3];
    let partial = far
        .evaluate_coherence_modes(&factors, 2, ANGULAR_BUDGET, &mut density, &gate)
        .unwrap();
    let mut expected = a;
    expected.total_power_w += b.total_power_w;
    for i in 0..3 {
        real_close(
            density[i],
            (first[i].norm_sq() + second[i].norm_sq())
                / (2.0 * model.gas().density * model.gas().sound_speed),
        );
        expected.first_directional_moment_w[i] += b.first_directional_moment_w[i];
        for j in 0..3 {
            expected.second_directional_moment_w[i][j] += b.second_directional_moment_w[i][j];
        }
    }
    moments_close(partial, expected, 3e-9);

    // One independent generator still drives two mutually coherent openings.
    // Compare independently solved individual opening fields, and explicitly
    // demonstrate that dropping their cross term gives a different result.
    let mut response = prepared.response_buffer();
    prepared
        .drive_into(&[C64::ONE, C64::ZERO], &mut response, &gate)
        .unwrap();
    let flows = response.outward_flows_m3_s().unwrap();
    let left = direct_far_field(
        &model,
        &direct_solution(&model, omega, &[flows[0], C64::ZERO]),
        &DIRECTIONS,
    );
    let right = direct_far_field(
        &model,
        &direct_solution(&model, omega, &[C64::ZERO, flows[1]]),
        &DIRECTIONS,
    );
    far.observe_independent_sources(&[1.0, 0.0], ANGULAR_BUDGET, &mut density, &gate)
        .unwrap();
    let denominator = 2.0 * model.gas().density * model.gas().sound_speed;
    let mut differs = false;
    for i in 0..3 {
        let correct = (left[i] + right[i]).norm_sq() / denominator;
        let diagonal_openings = (left[i].norm_sq() + right[i].norm_sq()) / denominator;
        real_close(density[i], correct);
        differs |= (correct - diagonal_openings).abs() > 1e-3 * correct;
    }
    assert!(differs, "fixture must resolve opening interference");
}

#[test]
fn g0_angular_preparation_and_application_keep_budgets_cancellation_and_silence_explicit() {
    let model = asymmetric_model();
    let sources = mixed_sources();
    let omega = core::f64::consts::TAU * 1200.0;
    let gate = CancelGate::new_clock_free();
    let mut too_small = options();
    too_small.maximum_source_entries = 1;
    assert!(
        model
            .prepare_sources_with_radiation(
                omega,
                LossModel::Lossless,
                &sources,
                None,
                &[],
                &DIRECTIONS,
                too_small,
                &gate,
            )
            .is_err()
    );
    let prepared = model
        .prepare_sources_with_radiation(
            omega,
            LossModel::Lossless,
            &sources,
            None,
            &[],
            &DIRECTIONS,
            options(),
            &gate,
        )
        .unwrap();
    let far = prepared.far_field().unwrap();
    let mut output = [C64::ZERO; 3];
    // Direction-only BEM work is 3*(2+1+10)=39. Four loaded mapping products
    // also belong to this call; a budget that omits them must refuse.
    assert!(
        far.observe_coherent(
            prepared.default_excitations(),
            HarmonicWorkBudget {
                maximum_products: 39
            },
            &mut output,
            &gate
        )
        .is_err()
    );
    far.observe_coherent(
        prepared.default_excitations(),
        HarmonicWorkBudget {
            maximum_products: 43,
        },
        &mut output,
        &gate,
    )
    .unwrap();
    assert!(
        far.evaluate_coherent(
            prepared.default_excitations(),
            HarmonicWorkBudget {
                maximum_products: 43
            },
            &mut output,
            &gate
        )
        .is_err()
    );
    assert!(
        far.observe_coherent(
            &[C64::new(f64::NAN, 0.0); 2],
            ANGULAR_BUDGET,
            &mut output,
            &gate
        )
        .is_err()
    );
    assert!(
        far.observe_coherent(
            prepared.default_excitations(),
            ANGULAR_BUDGET,
            &mut [],
            &gate
        )
        .is_err()
    );
    let mut density = [1.0; 3];
    assert!(
        far.observe_coherence_modes(&[], 33, ANGULAR_BUDGET, &mut density, &gate)
            .is_err()
    );
    assert!(
        far.observe_independent_sources(&[-1.0, 1.0], ANGULAR_BUDGET, &mut density, &gate)
            .is_err()
    );
    let silence = far
        .evaluate_coherence_modes(&[], 0, ANGULAR_BUDGET, &mut density, &gate)
        .unwrap();
    assert_eq!(silence.total_power_w, 0.0);
    assert_eq!(density, [0.0; 3]);
    let stopped = CancelGate::new_clock_free();
    stopped.request();
    assert!(matches!(
        far.evaluate_coherent(
            prepared.default_excitations(),
            ANGULAR_BUDGET,
            &mut output,
            &stopped
        ),
        Err(ExteriorNetworkError::Cancelled)
    ));
    far.observe_coherent(
        prepared.default_excitations(),
        ANGULAR_BUDGET,
        &mut output,
        &gate,
    )
    .unwrap();
}
