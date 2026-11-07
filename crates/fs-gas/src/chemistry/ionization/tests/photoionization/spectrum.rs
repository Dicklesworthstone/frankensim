//! G1/G3 spectral attenuation, conservative heating and time/order refinement.
mod column;

use super::*;
use crate::chemistry::photoionization::spectrum::{
    SpectralPhotoionizationError, SpectralPhotoionizationResult, SpectralPhotoionizationSolve,
    advance_spectral_photoionization,
};

fn solve_controls() -> SpectralPhotoionizationSolve {
    let mut controls =
        SpectralPhotoionizationSolve::new(Temperature::new(1000.0), Temperature::new(50000.0));
    controls.relative_tolerance = 1e-7;
    controls.absolute_amount_tolerance_mol = 1e-13;
    controls.relative_charge_tolerance = 1e-12;
    controls.initial_step_s = 0.1;
    controls.maximum_step_s = 1.0;
    controls
}

fn spectral_source(channel: AtomicIonizationV1) -> AtomicPhotoionizationV1 {
    AtomicPhotoionizationV1::new(
        channel,
        Provenance {
            source: "manufactured two-band optical-depth and photon-energy oracle".into(),
            license: "analytic test fixture".into(),
            artifact: None,
        },
        vec![
            PhotoionizationCrossSectionPointV1 {
                photon_energy: Energy::new(14.0 * EV),
                cross_section: Area::new(SIGMA),
            },
            PhotoionizationCrossSectionPointV1 {
                photon_energy: Energy::new(16.0 * EV),
                cross_section: Area::new(SIGMA),
            },
            PhotoionizationCrossSectionPointV1 {
                photon_energy: Energy::new(20.0 * EV),
                cross_section: Area::new(0.2 * SIGMA),
            },
        ],
    )
    .unwrap()
}
fn groups() -> [PhotonNumberFluxV1; 2] {
    [
        PhotonNumberFluxV1 {
            photon_energy: Energy::new(16.0 * EV),
            number_flux_per_m2_s: 1.5 / SIGMA,
        },
        PhotonNumberFluxV1 {
            photon_energy: Energy::new(20.0 * EV),
            number_flux_per_m2_s: 0.8 / (0.2 * SIGMA),
        },
    ]
}
fn area() -> Area {
    Area::new(SIGMA * NA / 3.0)
}
fn spectral_run(
    initial: &FrozenIonizedGasStateV1,
    source: &AtomicPhotoionizationV1,
    groups: &[PhotonNumberFluxV1],
    controls: SpectralPhotoionizationSolve,
) -> SpectralPhotoionizationResult {
    advance_spectral_photoionization(
        initial,
        source,
        groups,
        area(),
        Time::new(1.0),
        controls,
        || false,
    )
    .unwrap()
}

// Independent explicit RK4 of d[n_A,abs_0,abs_1]/dt. This never calls either
// production depletion kernel or the split integrator, and uses std expm1.
fn independent_reference() -> [f64; 3] {
    let rate = |state: [f64; 3]| {
        let soft = 0.5 * -(-3.0 * state[0]).exp_m1();
        let hard = (4.0 / 3.0) * -(-0.6 * state[0]).exp_m1();
        [-soft - hard, soft, hard]
    };
    let add = |a: [f64; 3], b: [f64; 3], h: f64| std::array::from_fn(|i| a[i] + h * b[i]);
    let mut state = [1.0, 0.0, 0.0];
    let dt = 1.0 / 8192.0;
    for _ in 0..8192 {
        let k1 = rate(state);
        let k2 = rate(add(state, k1, dt * 0.5));
        let k3 = rate(add(state, k2, dt * 0.5));
        let k4 = rate(add(state, k3, dt));
        state = std::array::from_fn(|i| {
            state[i] + dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0
        });
    }
    state
}

#[test]
fn one_band_matches_exact_thin_thick_and_neutral_tail_updates() {
    for seed in [0.0, 0.1] {
        let (initial, atomic) = initial(seed, 12000.0, 4000.0, 0.0, 0.0);
        let source = photon_source(atomic);
        for (tau, dose) in [(1.0, 1.0), (1e-12, 50.0), (1e6, 100.0)] {
            let beam = exposure(1.0 - seed, tau, dose);
            let exact = advance_monochromatic_photoionization(
                &initial,
                &source,
                beam,
                photo_controls(),
                || false,
            )
            .unwrap();
            let mut controls = solve_controls();
            controls.initial_step_s = 1.0;
            let result = advance_spectral_photoionization(
                &initial,
                &source,
                &[beam.beam],
                beam.illuminated_area,
                beam.duration,
                controls,
                || false,
            )
            .unwrap();
            relative(
                result.state().amounts()[0].value(),
                exact.state().amounts()[0].value(),
                1e-11,
            );
            relative(
                result.reaction_extent().value(),
                exact.reaction_extent().value(),
                1e-13,
            );
            relative(
                result.state().electron_temperature().value(),
                exact.state().electron_temperature().value(),
                1e-13,
            );
            relative(
                result.balance().absorbed_photons,
                exact.balance().absorbed_photons,
                1e-13,
            );
            if exact.balance().transmitted_photons > 0.0 {
                relative(
                    result.balance().transmitted_photons,
                    exact.balance().transmitted_photons,
                    1e-11,
                );
            }
            assert_eq!(result.accepted_steps(), 1);
            assert_eq!(result.rejected_steps(), 0);
        }
    }
}

#[test]
fn concurrent_bands_match_independent_ode_and_conserve_each_photon_energy() {
    let (initial, atomic) = initial(0.0, 12000.0, 4000.0, 0.0, 0.0);
    let source = spectral_source(atomic);
    let groups = groups();
    let result = spectral_run(&initial, &source, &groups, solve_controls());
    let expected = independent_reference();
    relative(result.state().amounts()[0].value(), expected[0], 2e-7);
    let mut deposited = 0.0;
    for (g, balance) in result.group_balances().iter().enumerate() {
        relative(balance.absorbed_photons / NA, expected[g + 1], 2e-7);
        relative(
            balance.incident_photons,
            groups[g].number_flux_per_m2_s * area().value(),
            1e-14,
        );
        relative(
            balance.absorbed_photons + balance.transmitted_photons,
            balance.incident_photons,
            2e-14,
        );
        relative(
            balance.absorbed_energy_j,
            balance.absorbed_photons * groups[g].photon_energy.value(),
            1e-14,
        );
        deposited += expected[g + 1] * (groups[g].photon_energy.value() * NA - CHI);
    }
    let expected_temperature = deposited / (1.5 * R * (expected[1] + expected[2]));
    relative(
        result.state().electron_temperature().value(),
        expected_temperature,
        2e-7,
    );
    assert_eq!(result.state().amounts()[1], result.state().amounts()[2]);
    relative(result.state().mass().value(), initial.mass().value(), 1e-14);
    assert_eq!(
        result.state().heavy_temperature(),
        initial.heavy_temperature()
    );
    assert_eq!(result.state().thermodynamics(), initial.thermodynamics());
    assert!(result.balance().energy_residual_j.abs() <= solve_controls().energy_tolerance_j);
    assert!(result.maximum_normalized_local_error() <= 1.0);
    // More hard photons penetrate. A flux-weighted mean photon energy would
    // therefore give the wrong absorbed energy and newborn electron temperature.
    let budgets = result.group_balances();
    assert!(
        budgets[1].transmitted_photons / budgets[1].incident_photons
            > budgets[0].transmitted_photons / budgets[0].incident_photons
    );
    let incident_mean = (groups[0].photon_energy.value() * groups[0].number_flux_per_m2_s
        + groups[1].photon_energy.value() * groups[1].number_flux_per_m2_s)
        / (groups[0].number_flux_per_m2_s + groups[1].number_flux_per_m2_s);
    let false_temperature = (incident_mean * NA - CHI) / (1.5 * R);
    assert!((false_temperature / expected_temperature - 1.0).abs() > 0.1);
}

#[test]
fn second_order_time_refinement_removes_group_order_bias() {
    let (initial, atomic) = initial(0.0, 12000.0, 4000.0, 0.0, 0.0);
    let source = spectral_source(atomic);
    let forward = groups();
    let reverse = [forward[1], forward[0]];
    let expected = independent_reference();
    let mut errors = Vec::new();
    let mut order_bias = Vec::new();
    for h in [0.25, 0.125, 0.0625] {
        let mut controls = solve_controls();
        controls.relative_tolerance = 0.1;
        controls.initial_step_s = h;
        controls.maximum_step_s = h;
        let a = spectral_run(&initial, &source, &forward, controls);
        let b = spectral_run(&initial, &source, &reverse, controls);
        assert_eq!(a.rejected_steps(), 0);
        assert_eq!(b.rejected_steps(), 0);
        errors.push((a.state().amounts()[0].value() - expected[0]).abs());
        order_bias.push((a.state().amounts()[0].value() - b.state().amounts()[0].value()).abs());
    }
    for values in [&errors, &order_bias] {
        for pair in values.windows(2) {
            assert!(
                pair[0] / pair[1] > 3.8 && pair[0] / pair[1] < 4.2,
                "{values:?}"
            );
        }
    }
}

#[test]
fn empty_dark_and_subthreshold_groups_preserve_gas_and_pass_photons() {
    let (initial, atomic) = initial(0.1, 12000.0, 4000.0, 0.0, 0.0);
    let source = spectral_source(atomic);
    let dark = PhotonNumberFluxV1 {
        photon_energy: Energy::new(16.0 * EV),
        number_flux_per_m2_s: 0.0,
    };
    let subthreshold = PhotonNumberFluxV1 {
        photon_energy: Energy::new(10.0 * EV),
        number_flux_per_m2_s: 1e20,
    };
    for groups in [vec![], vec![dark], vec![dark, subthreshold]] {
        let result = spectral_run(&initial, &source, &groups, solve_controls());
        assert_eq!(result.state(), &initial);
        assert_eq!(result.reaction_extent().value(), 0.0);
        assert_eq!(result.accepted_steps(), 0);
        for balance in result.group_balances() {
            assert_eq!(balance.absorbed_photons, 0.0);
            assert_eq!(balance.incident_photons, balance.transmitted_photons);
        }
    }
}

#[test]
fn mixed_birth_temperatures_use_the_combined_physical_electron_pool() {
    let (initial, atomic) = initial(0.0, 12000.0, 4000.0, 0.0, 0.0);
    let source = photon_source(atomic);
    let groups = [14.0, 22.0].map(|ev| PhotonNumberFluxV1 {
        photon_energy: Energy::new(ev * EV),
        number_flux_per_m2_s: 1.0 / SIGMA,
    });
    let mut controls = solve_controls();
    controls.minimum_electron_temperature_k = 5000.0;
    controls.maximum_electron_temperature_k = 45000.0;
    let result = spectral_run(&initial, &source, &groups, controls);
    // The two equally absorbed bands individually create electrons below/above
    // the chosen support. Their physical mixture is supported throughout t>0.
    let low = (14.0 * EV * NA - CHI) / (1.5 * R);
    let high = (22.0 * EV * NA - CHI) / (1.5 * R);
    assert!(low < controls.minimum_electron_temperature_k);
    assert!(high > controls.maximum_electron_temperature_k);
    relative(
        result.state().electron_temperature().value(),
        0.5 * (low + high),
        3e-7,
    );
}

#[test]
fn original_formation_references_cancel_before_spectral_energy_differences() {
    let mut outputs = Vec::new();
    for (elemental, charge) in [(0.0, 0.0), (1e12, 1e8)] {
        let (initial, atomic) = initial(0.05, 20000.0, 4000.0, elemental, charge);
        let source = spectral_source(atomic);
        let mut controls = solve_controls();
        controls.relative_tolerance = 0.1;
        controls.initial_step_s = 0.125;
        controls.maximum_step_s = 0.125;
        let result = spectral_run(&initial, &source, &groups(), controls);
        assert!(result.balance().energy_residual_j.abs() <= controls.energy_tolerance_j);
        outputs.push(result);
    }
    assert_eq!(outputs[0].state().amounts(), outputs[1].state().amounts());
    relative(
        outputs[0].state().electron_temperature().value(),
        outputs[1].state().electron_temperature().value(),
        1e-9,
    );
}

#[test]
fn source_budgets_and_late_cancellation_never_publish_partial_ionization() {
    let (initial, atomic) = initial(0.0, 12000.0, 4000.0, 0.0, 0.0);
    let source = spectral_source(atomic);
    let before = initial.clone();
    let groups = groups();
    let mut controls = solve_controls();
    controls.maximum_groups = 1;
    assert!(matches!(
        advance_spectral_photoionization(
            &initial,
            &source,
            &groups,
            area(),
            Time::new(1.0),
            controls,
            || false
        ),
        Err(SpectralPhotoionizationError::GroupLimit { .. })
    ));
    controls = solve_controls();
    controls.maximum_steps = 1;
    controls.initial_step_s = 1e-4;
    assert_eq!(
        advance_spectral_photoionization(
            &initial,
            &source,
            &groups,
            area(),
            Time::new(1.0),
            controls,
            || false
        ),
        Err(SpectralPhotoionizationError::StepBudget)
    );
    let mut calls = 0;
    assert_eq!(
        advance_spectral_photoionization(
            &initial,
            &source,
            &groups,
            area(),
            Time::new(1.0),
            solve_controls(),
            || {
                calls += 1;
                calls >= 40
            }
        ),
        Err(SpectralPhotoionizationError::Cancelled)
    );
    let unsupported_dark = [PhotonNumberFluxV1 {
        photon_energy: Energy::new(21.0 * EV),
        number_flux_per_m2_s: 0.0,
    }];
    assert!(matches!(
        advance_spectral_photoionization(
            &initial,
            &source,
            &unsupported_dark,
            area(),
            Time::new(1.0),
            solve_controls(),
            || false
        ),
        Err(SpectralPhotoionizationError::Group { .. })
    ));
    let mut cold = solve_controls();
    cold.maximum_electron_temperature_k = 13000.0;
    assert_eq!(
        advance_spectral_photoionization(
            &initial,
            &source,
            &groups,
            area(),
            Time::new(1.0),
            cold,
            || false
        ),
        Err(SpectralPhotoionizationError::OutsideDomain)
    );
    assert_eq!(initial, before);
}
