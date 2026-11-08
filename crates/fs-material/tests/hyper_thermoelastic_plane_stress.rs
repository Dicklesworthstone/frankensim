//! Native transverse relaxation retains full finite-strain thermodynamics.
use fs_material::hyper::mixed::VolumetricTreatment;
use fs_material::hyper::plane_stress::{PlaneStressSettings, plane_stress_with_seed};
use fs_material::hyper::reference::{ReferenceDeformation, ReferencedHyperelastic};
use fs_material::hyper::thermoelastic::plane_stress::FiniteThermoelasticPlaneStressState;
use fs_material::hyper::thermoelastic::{FiniteThermoelasticState, ThermalReferenceJet};
use fs_material::hyper::{Hyperelastic, HyperelasticModel};
use fs_material::visco_hyper::{BranchEnergy, FiniteStrainBranch, ViscoHyperelastic};

const G: [f64; 9] = [1.09, 0.07, -0.025, 0.03, 0.97, 0.04, -0.015, 0.02, 1.04];
const G1: [f64; 9] = [
    0.002, 0.001, -0.0004, -0.0003, -0.001, 0.0006, 0.0002, 0.0008, 0.0015,
];
const G2: [f64; 9] = [2e-5, -3e-5, 1e-5, 1e-5, 4e-5, 0.0, -2e-5, 1e-5, -1e-5];
const A: [f64; 4] = [1.22, 0.16, -0.05, 0.91];

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

fn settings() -> PlaneStressSettings {
    PlaneStressSettings {
        relative_tolerance: 2e-12,
        ..Default::default()
    }
}

fn near(actual: f64, expected: f64, absolute: f64, relative: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    let allowance = absolute + relative * expected.abs();
    assert!(
        (actual - expected).abs() <= allowance,
        "{actual:e} != {expected:e}, allowance {allowance:e}"
    );
}

fn diagonal(value: f64) -> [f64; 9] {
    [value, 0.0, 0.0, 0.0, value, 0.0, 0.0, 0.0, value]
}

fn jet(delta_k: f64) -> ThermalReferenceJet {
    let reference = ReferenceDeformation::try_new(core::array::from_fn(|i| {
        G[i] + delta_k * G1[i] + 0.5 * delta_k * delta_k * G2[i]
    }))
    .unwrap();
    ThermalReferenceJet::try_new(
        reference,
        core::array::from_fn(|i| G1[i] + delta_k * G2[i]),
        G2,
    )
    .unwrap()
}

fn derivative(values: [f64; 4], step: f64) -> f64 {
    (values[0] - 8.0 * values[1] + 8.0 * values[2] - values[3]) / (12.0 * step)
}

#[test]
fn isotropic_natural_state_matches_closed_form_thermal_schur_elimination() {
    let (mu, lambda, stretch, rate, temperature) = (6.0, 12.0, 1.1, 0.002, 300.0);
    let law = Hyperelastic::new(HyperelasticModel::NeoHookean { mu, lambda }, 4.0).unwrap();
    let jet = ThermalReferenceJet::try_new(
        ReferenceDeformation::try_new(diagonal(stretch)).unwrap(),
        diagonal(rate),
        [0.0; 9],
    )
    .unwrap();
    let state = FiniteThermoelasticPlaneStressState::try_new(
        &law,
        &[stretch, 0.0, 0.0, stretch],
        temperature,
        jet,
        settings(),
    )
    .unwrap();
    let bulk = lambda + 2.0 * mu / 3.0;
    let beta = -3.0 * bulk * stretch * rate;
    let normal = stretch * (lambda + 2.0 * mu);
    let beta_reduced = beta * (2.0 * mu / (lambda + 2.0 * mu));
    let full_capacity = -9.0 * temperature * bulk * stretch * rate * rate;
    let capacity = full_capacity + temperature * beta * beta / normal;
    near(
        state.full_state().mechanical_heat_capacity_j_m3_k(),
        full_capacity,
        2e-12,
        2e-11,
    );
    near(
        state.mechanical_heat_capacity_j_m3_k(),
        capacity,
        2e-12,
        2e-11,
    );
    assert!(state.mechanical_heat_capacity_j_m3_k() > full_capacity);
    for (actual, expected) in state.stress_temperature_derivative_pa_k().iter().zip([
        beta_reduced,
        0.0,
        0.0,
        beta_reduced,
    ]) {
        near(*actual, expected, 2e-12, 2e-11);
    }
    for (actual, expected) in
        state
            .transverse_temperature_derivative_per_k()
            .iter()
            .zip([0.0, 0.0, -beta / normal])
    {
        near(*actual, expected, 2e-12, 2e-11);
    }
    let plane_lambda = 2.0 * mu * lambda / (lambda + 2.0 * mu);
    let expected = [
        [plane_lambda + 2.0 * mu, 0.0, 0.0, plane_lambda],
        [0.0, mu, mu, 0.0],
        [0.0, mu, mu, 0.0],
        [plane_lambda, 0.0, 0.0, plane_lambda + 2.0 * mu],
    ];
    for (actual, expected) in state
        .response()
        .tangent
        .iter()
        .flatten()
        .zip(expected.iter().flatten())
    {
        near(*actual, stretch * expected, 2e-11, 2e-11);
    }
    for (actual, expected) in state
        .full_state()
        .deformation_gradient()
        .iter()
        .zip(diagonal(stretch))
    {
        near(*actual, expected, 2e-12, 2e-12);
    }
    near(state.free_energy_j_m3(), 0.0, 2e-12, 0.0);
    near(state.entropy_j_m3_k(), 0.0, 2e-12, 0.0);
    near(state.internal_energy_j_m3(), 0.0, 2e-10, 0.0);
}

#[test]
fn general_thermal_map_matches_relaxed_energy_stress_and_entropy_derivatives() {
    let law = law();
    let temperature = 310.0;
    let state =
        FiniteThermoelasticPlaneStressState::try_new(&law, &A, temperature, jet(0.0), settings())
            .unwrap();
    let step = 0.05;
    // Only the original elastic potential and native transverse solve enter
    // these temperature differences; no thermal derivative API supplies them.
    let samples = [-2.0 * step, -step, step, 2.0 * step].map(|delta| {
        let reference = *jet(delta).reference();
        plane_stress_with_seed(
            &ReferencedHyperelastic::new(&law, reference),
            &A,
            [
                0.0,
                0.0,
                reference.volume_ratio() / (A[0] * A[3] - A[1] * A[2]),
            ],
            settings(),
        )
        .unwrap()
    });
    let energies = samples.map(|p| p.energy_density_pa);
    let psi_tt = (-energies[3] + 16.0 * energies[2] - 30.0 * state.free_energy_j_m3()
        + 16.0 * energies[1]
        - energies[0])
        / (12.0 * step * step);
    near(
        state.entropy_j_m3_k(),
        -derivative(energies, step),
        2e-9,
        2e-7,
    );
    near(
        state.mechanical_heat_capacity_j_m3_k(),
        -temperature * psi_tt,
        3e-7,
        3e-6,
    );
    for i in 0..4 {
        near(
            state.stress_temperature_derivative_pa_k()[i],
            derivative(samples.map(|p| p.piola[i]), step),
            2e-9,
            2e-6,
        );
        // Independent Maxwell relation for the reduced potential: beta_bar
        // equals -dS_bar/da while the transverse state is allowed to relax.
        let h = 2e-5;
        let entropy = [-2.0 * h, -h, h, 2.0 * h].map(|delta| {
            let mut a = A;
            a[i] += delta;
            FiniteThermoelasticPlaneStressState::try_new(
                &law,
                &a,
                temperature,
                jet(0.0),
                settings(),
            )
            .unwrap()
            .entropy_j_m3_k()
        });
        near(
            state.stress_temperature_derivative_pa_k()[i],
            -derivative(entropy, h),
            2e-8,
            2e-6,
        );
    }
    for (i, dof) in [2, 5, 8].into_iter().enumerate() {
        near(
            state.transverse_temperature_derivative_per_k()[i],
            derivative(samples.map(|p| p.deformation_gradient[dof]), step),
            2e-9,
            2e-6,
        );
    }
    assert!(state.transverse_temperature_derivative_per_k()[0].abs() > 1e-6);
    assert!(state.transverse_temperature_derivative_per_k()[1].abs() > 1e-6);
    assert!(
        state.mechanical_heat_capacity_j_m3_k()
            - state.full_state().mechanical_heat_capacity_j_m3_k()
            > 1e-6
    );
    assert_eq!(state.entropy_j_m3_k(), state.full_state().entropy_j_m3_k());
    assert_eq!(
        state.internal_energy_j_m3(),
        state.full_state().internal_energy_j_m3()
    );
}

#[test]
fn snapshot_reduction_retains_its_owner_and_refuses_unrelaxed_or_constrained_states() {
    let law = law();
    let state = FiniteThermoelasticPlaneStressState::try_new(&law, &A, 310.0, jet(0.0), settings())
        .unwrap();
    let full = state.full_state().clone();
    let reduced = FiniteThermoelasticPlaneStressState::try_from_state(
        full.clone(),
        VolumetricTreatment::OriginalMaterial,
        settings(),
    )
    .unwrap();
    assert_eq!(reduced.full_state(), &full);
    assert_eq!(reduced.response().iterations, 0);
    assert_eq!(
        reduced.response().deformation_gradient,
        state.response().deformation_gradient
    );
    assert_eq!(reduced.response().tangent, state.response().tangent);
    assert_eq!(
        reduced.stress_temperature_derivative_pa_k(),
        state.stress_temperature_derivative_pa_k()
    );
    assert_eq!(
        reduced.mechanical_heat_capacity_j_m3_k(),
        state.mechanical_heat_capacity_j_m3_k()
    );
    assert!(
        FiniteThermoelasticPlaneStressState::try_from_state(
            full.clone(),
            VolumetricTreatment::ExactIncompressible,
            settings()
        )
        .is_err()
    );
    let unrelaxed =
        FiniteThermoelasticState::try_new(&law, &diagonal(1.4), 310.0, jet(0.0)).unwrap();
    assert!(
        FiniteThermoelasticPlaneStressState::try_from_state(
            unrelaxed,
            VolumetricTreatment::OriginalMaterial,
            settings()
        )
        .is_err()
    );
    let mut tilted = *full.deformation_gradient();
    tilted[6] = 0.02;
    let tilted = FiniteThermoelasticState::try_new(&law, &tilted, 310.0, jet(0.0)).unwrap();
    assert!(
        FiniteThermoelasticPlaneStressState::try_from_state(
            tilted,
            VolumetricTreatment::OriginalMaterial,
            settings()
        )
        .is_err()
    );
    assert!(
        FiniteThermoelasticPlaneStressState::try_new(
            &law,
            &[1.0, 2.0, 2.0, 1.0],
            310.0,
            jet(0.0),
            settings()
        )
        .is_err()
    );
    assert!(
        FiniteThermoelasticPlaneStressState::try_new(&law, &A, 0.0, jet(0.0), settings()).is_err()
    );
    let constant = ThermalReferenceJet::try_new(*jet(0.0).reference(), [0.0; 9], [0.0; 9]).unwrap();
    let no_thermal_work =
        FiniteThermoelasticPlaneStressState::try_new(&law, &A, 310.0, constant, settings())
            .unwrap();
    assert_eq!(
        no_thermal_work.stress_temperature_derivative_pa_k(),
        &[0.0; 4]
    );
    assert_eq!(no_thermal_work.mechanical_heat_capacity_j_m3_k(), 0.0);
    assert_eq!(
        no_thermal_work.internal_energy_j_m3(),
        no_thermal_work.free_energy_j_m3()
    );
}

#[test]
fn aged_viscoelastic_snapshot_uses_the_same_reduction_without_changing_history() {
    let mut material = ViscoHyperelastic::new(
        law(),
        vec![FiniteStrainBranch {
            energy: BranchEnergy::NeoHookean { mu_pa: 7.0 },
            tau_s: 0.2,
        }],
    )
    .unwrap()
    .with_reference_deformation(*jet(0.0).reference());
    let accepted = material
        .step_plane_stress_response(
            &material.initial_state(),
            &[1.18, 0.1, -0.02, 0.94],
            0.06,
            None,
            settings(),
        )
        .unwrap()
        .material
        .state;
    assert!(accepted.dissipated > 0.0);
    let saved = accepted.clone();
    let state = material
        .frozen_thermoelastic_plane_stress_state(&accepted, &A, 310.0, jet(0.0), None, settings())
        .unwrap();
    let frozen = material.freeze(&accepted, 0.0).unwrap();
    let native = plane_stress_with_seed(
        &frozen,
        &A,
        [
            accepted.deformation_gradient[2],
            accepted.deformation_gradient[5],
            jet(0.0).reference().volume_ratio() / (A[0] * A[3] - A[1] * A[2]),
        ],
        settings(),
    )
    .unwrap();
    near(
        state.free_energy_j_m3(),
        native.energy_density_pa,
        2e-11,
        2e-11,
    );
    for i in 0..4 {
        near(state.response().piola[i], native.piola[i], 2e-11, 2e-11);
        for j in 0..4 {
            near(
                state.response().tangent[i][j],
                native.tangent[i][j],
                2e-10,
                2e-11,
            );
        }
    }
    material.branches[0].tau_s *= 10_000.0;
    let slow = material
        .frozen_thermoelastic_plane_stress_state(&accepted, &A, 310.0, jet(0.0), None, settings())
        .unwrap();
    assert_eq!(slow.full_state(), state.full_state());
    assert_eq!(
        slow.stress_temperature_derivative_pa_k(),
        state.stress_temperature_derivative_pa_k()
    );
    assert_eq!(
        slow.mechanical_heat_capacity_j_m3_k(),
        state.mechanical_heat_capacity_j_m3_k()
    );
    let rates = state.rate([0.08, -0.03, 0.01, -0.05], 3.0).unwrap();
    let slow_rates = slow.rate([0.08, -0.03, 0.01, -0.05], 3.0).unwrap();
    assert_eq!(rates, slow_rates);
    for i in [2, 5, 8] {
        near(rates.first_piola_rate_pa_s()[i], 0.0, 3e-11, 0.0);
    }
    assert_eq!(accepted, saved);
}

#[test]
fn coupled_plane_rates_follow_the_relaxed_thermal_path_in_all_three_transverse_components() {
    let material = law();
    let temperature = 310.0;
    let adot = [0.08, -0.03, 0.01, -0.05];
    let tdot = 3.0;
    let state = FiniteThermoelasticPlaneStressState::try_new(
        &material,
        &A,
        temperature,
        jet(0.0),
        settings(),
    )
    .unwrap();
    let rate = state.rate(adot, tdot).unwrap();
    let h = 0.005;
    let samples = [-2.0 * h, -h, h, 2.0 * h].map(|time| {
        let a = core::array::from_fn(|i| A[i] + time * adot[i]);
        FiniteThermoelasticPlaneStressState::try_new(
            &material,
            &a,
            temperature + tdot * time,
            jet(tdot * time),
            settings(),
        )
        .unwrap()
    });
    for i in 0..9 {
        near(
            rate.deformation_gradient_rate_per_s()[i],
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
            rate.elastic_deformation_gradient_rate_per_s()[i],
            derivative(
                samples
                    .each_ref()
                    .map(|s| s.full_state().elastic_deformation_gradient()[i]),
                h,
            ),
            4e-8,
            3e-7,
        );
        near(
            rate.first_piola_rate_pa_s()[i],
            derivative(
                samples
                    .each_ref()
                    .map(|s| s.full_state().first_piola_pa()[i]),
                h,
            ),
            3e-7,
            3e-7,
        );
    }
    for (actual, values) in [
        (
            rate.free_energy_rate_w_m3(),
            samples.each_ref().map(|s| s.free_energy_j_m3()),
        ),
        (
            rate.entropy_rate_w_m3_k(),
            samples.each_ref().map(|s| s.entropy_j_m3_k()),
        ),
        (
            rate.internal_energy_rate_w_m3(),
            samples.each_ref().map(|s| s.internal_energy_j_m3()),
        ),
    ] {
        near(actual, derivative(values, h), 3e-7, 3e-7);
    }
    for i in [2, 5, 8] {
        near(rate.first_piola_rate_pa_s()[i], 0.0, 3e-12, 0.0);
    }
    for (i, component) in [0, 1, 3, 4].into_iter().enumerate() {
        let expected = state.stress_temperature_derivative_pa_k()[i] * tdot
            + (0..4)
                .map(|j| state.response().tangent[i][j] * adot[j])
                .sum::<f64>();
        near(
            rate.first_piola_rate_pa_s()[component],
            expected,
            3e-12,
            3e-12,
        );
    }
    let expected_entropy = state.mechanical_heat_capacity_j_m3_k() / temperature * tdot
        - (0..4)
            .map(|i| state.stress_temperature_derivative_pa_k()[i] * adot[i])
            .sum::<f64>();
    near(rate.entropy_rate_w_m3_k(), expected_entropy, 3e-13, 3e-12);
    let held = state.kinematics().deformation_gradient_rate(adot).unwrap();
    let isothermal = state.rate(adot, 0.0).unwrap();
    assert_eq!(isothermal.deformation_gradient_rate_per_s(), &held);
    // A frozen-T lift would violate transverse equilibrium on this path.
    let frozen_temperature = state.full_state().rate(held, tdot).unwrap();
    assert!(
        [2, 5, 8]
            .iter()
            .any(|&i| frozen_temperature.first_piola_rate_pa_s()[i].abs() > 1e-3)
    );
}

#[test]
fn plane_thermal_free_expansion_recovers_thickness_without_elastic_strain_or_heat() {
    let (stretch, expansion, tdot) = (1.1, 0.002, 3.0);
    let material = law();
    let thermal = ThermalReferenceJet::try_new(
        ReferenceDeformation::try_new(diagonal(stretch)).unwrap(),
        diagonal(expansion),
        [0.0; 9],
    )
    .unwrap();
    let state = FiniteThermoelasticPlaneStressState::try_new(
        &material,
        &[stretch, 0.0, 0.0, stretch],
        300.0,
        thermal,
        settings(),
    )
    .unwrap();
    let before = state.full_state().clone();
    let rate = state
        .rate([expansion * tdot, 0.0, 0.0, expansion * tdot], tdot)
        .unwrap();
    for (actual, expected) in rate
        .deformation_gradient_rate_per_s()
        .iter()
        .zip(diagonal(expansion * tdot))
    {
        near(*actual, expected, 2e-12, 2e-12);
    }
    for actual in rate
        .elastic_deformation_gradient_rate_per_s()
        .iter()
        .chain(rate.first_piola_rate_pa_s())
    {
        near(*actual, 0.0, 3e-12, 0.0);
    }
    for actual in [
        rate.free_energy_rate_w_m3(),
        rate.entropy_rate_w_m3_k(),
        rate.internal_energy_rate_w_m3(),
    ] {
        near(actual, 0.0, 3e-11, 0.0);
    }
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(state.rate([bad; 4], 0.0).is_err());
        assert!(state.rate([0.0; 4], bad).is_err());
    }
    assert!(state.rate([f64::MAX; 4], 0.0).is_err());
    assert_eq!(state.full_state(), &before);
}
