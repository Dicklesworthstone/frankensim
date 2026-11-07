//! Independent simultaneous reaction/attenuation references and physical limits.
use super::collisional_radiative::fixture;
use super::{CHI, EV, KB, ME, MP, NA, R, id};
use crate::chemistry::collisional_radiative::{
    CollisionalRadiativeSolve, CollisionalRadiativeThermal, advance_collisional_radiative,
};
use crate::chemistry::energy::ElectronHeavyCollisionSources;
use crate::chemistry::radiation_transport::*;
use fs_matdb::Provenance;
use fs_qty::{Area, Energy, Power, Temperature, Time, Volume};
use fs_thermochem::electron_collision::{
    ElectronCollisionKindV1, ElectronCrossSectionPointV1, ElectronCrossSectionV1,
    MaxwellianRateControlsV1,
};
use fs_thermochem::equilibrium::composition::FrozenIonizedGasStateV1;
use fs_thermochem::kinetics::atomic_ionization::photoionization::PhotonNumberFluxV1;

const SIGMA: f64 = 1e-22;
const CV: f64 = 1.5 * R;

fn controls(relative: f64, duration: f64) -> RadiationTransportSolve {
    let mut local =
        CollisionalRadiativeSolve::new(Temperature::new(1000.0), Temperature::new(50000.0));
    local.relative_tolerance = relative;
    local.initial_step_s = duration / 8.0;
    local.maximum_step_s = duration;
    local.energy_tolerance_j = 1e-6;
    let mut result = RadiationTransportSolve::new(local);
    result.maximum_cells = 4;
    result.maximum_groups = 4;
    result.maximum_cell_groups = 16;
    result
}
fn area() -> Area {
    Area::new(SIGMA * NA / 2.0)
}
fn beams(energies: [f64; 2], rates: [f64; 2]) -> [PhotonNumberFluxV1; 2] {
    std::array::from_fn(|g| PhotonNumberFluxV1 {
        photon_energy: Energy::new(energies[g]),
        number_flux_per_m2_s: rates[g] / SIGMA,
    })
}
fn resize(state: &FrozenIonizedGasStateV1, volume: f64) -> FrozenIonizedGasStateV1 {
    FrozenIonizedGasStateV1::new(
        state.thermodynamics(),
        state.amounts(),
        Volume::new(volume),
        state.heavy_temperature(),
        state.electron_temperature(),
        1e-12,
        || false,
    )
    .unwrap()
}
fn close(actual: f64, expected: f64, absolute: f64, relative: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    assert!(
        (actual - expected).abs() <= absolute + relative * expected.abs(),
        "{actual:.16e} != {expected:.16e}"
    );
}

// Each cell has [ions, Ue, absorbed_band0, absorbed_band1, net_collision,
// radiative_capture, emitted_energy, Uh]. Independent RK4 uses closed-form
// ramp-impact Maxwellian and constant-cross-section Milne moments, statistical
// detailed balance, neutral elastic exchange and instantaneous Beer-Lambert attenuation. No production
// rate evaluator, subflow, splitting or time-integration routine is called.
fn reference(
    initial: [&FrozenIonizedGasStateV1; 2],
    duration: f64,
    impact: f64,
    energies: [f64; 2],
    rates: [f64; 2],
    powers: [[f64; 2]; 2],
    momentum_cross_section: f64,
) -> [[f64; 8]; 2] {
    let gap = CHI / NA;
    let rate = |state: [[f64; 8]; 2]| {
        let mut result = [[0.0; 8]; 2];
        let mut incident = rates.map(|gamma| gamma / 2.0); // A F / NA [mol/s]
        for c in 0..2 {
            let x = state[c][0];
            let te = if x > 0.0 {
                state[c][1] / (CV * x)
            } else {
                let photons = rates[0] + rates[1];
                if photons > 0.0 {
                    ((rates[0] * energies[0] + rates[1] * energies[1]) * NA / photons - CHI) / CV
                } else {
                    initial[c].electron_temperature().value()
                }
            };
            let theta = KB * te;
            let volume = initial[c].volume().value();
            let slope = impact * 1e-18 / (1000.0 * EV);
            let kf = (8.0 / (core::f64::consts::PI * ME)).sqrt()
                * slope
                * theta.sqrt()
                * (gap + 2.0 * theta)
                * (-gap / theta).exp();
            let h = 6.626_070_15e-34;
            let kc = (2.0 * core::f64::consts::PI * theta * (MP * ME / (MP + ME)) / (h * h))
                .powf(1.5)
                * (-gap / theta).exp()
                / NA;
            let collision = NA * kf * (x * (1.0 - x) / volume - x * x * x / (kc * volume * volume));
            let p = gap * gap + 2.0 * gap * theta + 2.0 * theta * theta;
            let alpha = (8.0 / (core::f64::consts::PI * ME)).sqrt() * SIGMA
                / (ME * 299_792_458.0_f64.powi(2))
                * p
                / theta.sqrt();
            let mean_capture = theta * (gap * gap + 4.0 * gap * theta + 6.0 * theta * theta) / p;
            let capture = alpha * NA * x * x / volume;
            // Constant source: energy-weighted momentum moment is 4/3 times
            // sigma*mean_speed, with the stationary-heavy Banks mass factor.
            let km = (4.0 / 3.0)
                * momentum_cross_section
                * (8.0 * theta / (core::f64::consts::PI * ME)).sqrt();
            let conductance = 3.0 * R * x * ((1.0 - x) / volume * NA * (ME / (MP + ME)) * km);
            let exchange = conductance * (te - state[c][7] / CV);
            let tau = 2.0 * (1.0 - x);
            let absorption = incident.map(|photons| photons * -(-tau).exp_m1());
            for photons in &mut incident {
                *photons *= (-tau).exp();
            }
            result[c][0] = absorption[0] + absorption[1] + collision - capture;
            result[c][1] = absorption[0] * (energies[0] * NA - CHI)
                + absorption[1] * (energies[1] * NA - CHI)
                - CHI * collision
                - capture * mean_capture * NA
                + powers[c][0]
                - exchange;
            result[c][2] = absorption[0];
            result[c][3] = absorption[1];
            result[c][4] = collision;
            result[c][5] = capture;
            result[c][6] = capture * (CHI + mean_capture * NA);
            result[c][7] = powers[c][1] + exchange;
        }
        result
    };
    let add = |a: [[f64; 8]; 2], b: [[f64; 8]; 2], dt: f64| {
        std::array::from_fn(|c| std::array::from_fn(|i| a[c][i] + dt * b[c][i]))
    };
    let mut state = std::array::from_fn(|c| {
        let x = initial[c].amounts()[1].value();
        [
            x,
            CV * x * initial[c].electron_temperature().value(),
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            CV * initial[c].heavy_temperature().value(),
        ]
    });
    let dt = duration / 16384.0;
    for _ in 0..16384 {
        let a = rate(state);
        let b = rate(add(state, a, dt * 0.5));
        let c = rate(add(state, b, dt * 0.5));
        let d = rate(add(state, c, dt));
        state = std::array::from_fn(|j| {
            std::array::from_fn(|k| {
                state[j][k] + dt * (a[j][k] + 2.0 * b[j][k] + 2.0 * c[j][k] + d[j][k]) / 6.0
            })
        });
    }
    state
}

fn check_energy_and_photons(result: &RadiationTransportResult, tolerance: f64) {
    let mut energy = 0.0;
    for (c, cell) in result.cells().iter().enumerate() {
        energy +=
            cell.state().internal_energy().value() - cell.initial_state().internal_energy().value();
        assert_eq!(cell.state().amounts()[1], cell.state().amounts()[2]);
        close(
            cell.state().mass().value(),
            cell.initial_state().mass().value(),
            0.0,
            2e-13,
        );
        assert_eq!(
            cell.state().thermodynamics(),
            cell.initial_state().thermodynamics()
        );
        for (g, b) in cell.group_balances().iter().enumerate() {
            close(
                b.absorbed_photons + b.transmitted_photons,
                b.incident_photons,
                0.0,
                2e-12,
            );
            if c > 0 {
                assert_eq!(
                    b.incident_photons,
                    result.cells()[c - 1].group_balances()[g].transmitted_photons
                );
            }
        }
    }
    for b in result.group_balances() {
        close(
            b.absorbed_photons + b.transmitted_photons,
            b.incident_photons,
            0.0,
            2e-12,
        );
    }
    let b = result.balance();
    assert!(b.energy_residual_j.abs() <= tolerance);
    assert!(
        (energy + b.emitted_radiation_energy_j
            - b.absorbed_radiation_energy_j
            - b.electron_deposited_energy_j
            - b.heavy_deposited_energy_j)
            .abs()
            < 4.0 * tolerance
    );
    close(
        b.photoionization_extent_mol * NA,
        b.absorbed_photons,
        0.0,
        2e-13,
    );
    assert!(b.maximum_normalized_local_error <= 1.0);
}

#[test]
fn competing_reactions_and_evolving_shielding_match_independent_two_cell_ode() {
    let (front, source) = fixture(0.1, 10000.0, 1.0, 0.0, 0.0);
    let (back, _) = fixture(0.12, 12000.0, 1.0, 0.0, 0.0);
    let back = resize(&back, 0.35);
    let energies = [16.0 * EV, 20.0 * EV];
    let rates = [400.0, 600.0];
    let duration = 4e-5;
    let result = advance_radiation_transport(
        &[
            RadiationTransportCell::new(&front, &source),
            RadiationTransportCell::new(&back, &source),
        ],
        &beams(energies, rates),
        area(),
        Time::new(duration),
        controls(1e-7, duration),
        || false,
    )
    .unwrap();
    let expected = reference(
        [&front, &back],
        duration,
        1.0,
        energies,
        rates,
        [[0.0; 2]; 2],
        0.0,
    );
    for (c, cell) in result.cells().iter().enumerate() {
        close(cell.state().amounts()[1].value(), expected[c][0], 3e-7, 0.0);
        close(
            cell.state().electron_temperature().value(),
            expected[c][1] / (CV * expected[c][0]),
            0.05,
            0.0,
        );
        for g in 0..2 {
            close(
                cell.group_balances()[g].absorbed_photons / NA,
                expected[c][2 + g],
                3e-7,
                0.0,
            );
        }
        close(
            cell.balance().collisional_extent_mol,
            expected[c][4],
            3e-7,
            0.0,
        );
        close(
            cell.balance().radiative_recombination_extent_mol,
            expected[c][5],
            3e-7,
            0.0,
        );
        close(
            cell.balance().emitted_radiation_energy_j,
            expected[c][6],
            0.5,
            0.0,
        );
        assert_eq!(
            cell.state().heavy_temperature(),
            cell.initial_state().heavy_temperature()
        );
    }
    // Restored neutral population increases shielding during illumination;
    // a photoionization-only column cannot reproduce this direction of change.
    assert!(result.cells()[0].state().amounts()[0].value() > front.amounts()[0].value());
    assert!(
        result.balance().radiative_recombination_extent_mol
            > result.balance().photoionization_extent_mol
    );
    check_energy_and_photons(&result, 1e-6);
}

#[test]
fn photon_budget_survives_repeated_creation_and_capture_from_zero_electrons() {
    let gap = CHI / NA;
    let theta = gap / (1.0 + 7.0_f64.sqrt());
    let energy = gap + 1.5 * theta;
    let (initial, source) = fixture(0.0, 10000.0, 0.0, 0.0, 0.0);
    let duration = 5e-3;
    let rates = [400.0, 600.0];
    let result = advance_radiation_transport(
        &[RadiationTransportCell::new(&initial, &source); 2],
        &beams([energy; 2], rates),
        area(),
        Time::new(duration),
        controls(1e-6, duration),
        || false,
    )
    .unwrap();
    let expected = reference(
        [&initial; 2],
        duration,
        0.0,
        [energy; 2],
        rates,
        [[0.0; 2]; 2],
        0.0,
    );
    for (c, cell) in result.cells().iter().enumerate() {
        close(cell.state().amounts()[1].value(), expected[c][0], 3e-6, 0.0);
        close(
            cell.state().electron_temperature().value(),
            theta / KB,
            3e-5,
            0.0,
        );
        close(
            cell.balance().radiative_recombination_extent_mol,
            expected[c][5],
            2e-5,
            0.0,
        );
    }
    assert!(result.balance().photoionization_extent_mol > 2.0);
    assert!(result.balance().radiative_recombination_extent_mol > 2.0);
    assert!(result.balance().collisional_extent_mol.abs() < 1e-12);
    check_energy_and_photons(&result, 1e-6);
}

#[test]
fn coupled_time_refinement_reduces_the_continuous_ode_error_at_second_order() {
    let (front, source) = fixture(0.1, 10000.0, 1.0, 0.0, 0.0);
    let (back, _) = fixture(0.12, 12000.0, 1.0, 0.0, 0.0);
    // Keep all four channels, but resolve the fastest density-dependent mode
    // with these modest fixed step counts. The dense, stiff case is covered
    // separately by the adaptive continuous-ODE and repeated-turnover tests.
    let front = resize(&front, 20.0);
    let back = resize(&back, 35.0);
    let energies = [16.0 * EV, 20.0 * EV];
    let rates = [400.0, 600.0];
    let duration = 2e-5;
    let expected = reference(
        [&front, &back],
        duration,
        1.0,
        energies,
        rates,
        [[0.0; 2]; 2],
        0.0,
    );
    let mut errors = Vec::new();
    for count in [16.0, 32.0, 64.0] {
        let mut solve = controls(0.1, duration);
        solve.chemistry.initial_step_s = duration / count;
        solve.chemistry.maximum_step_s = duration / count;
        let result = advance_radiation_transport(
            &[
                RadiationTransportCell::new(&front, &source),
                RadiationTransportCell::new(&back, &source),
            ],
            &beams(energies, rates),
            area(),
            Time::new(duration),
            solve,
            || false,
        )
        .unwrap();
        assert_eq!(result.balance().rejected_steps, 0);
        errors.push(
            result
                .cells()
                .iter()
                .enumerate()
                .map(|(c, cell)| {
                    (cell.state().amounts()[1].value() - expected[c][0]).abs()
                        + (cell.state().electron_temperature().value()
                            - expected[c][1] / (CV * expected[c][0]))
                            .abs()
                            / 100000.0
                })
                .sum::<f64>(),
        );
    }
    for pair in errors.windows(2) {
        assert!(
            pair[0] / pair[1] > 3.6 && pair[0] / pair[1] < 4.4,
            "{errors:?}"
        );
    }
}

#[test]
fn held_pool_powers_and_neutral_exchange_share_state_with_transport() {
    let (initial, source) = fixture(0.1, 10000.0, 1.0, 0.0, 0.0);
    let mut thermal =
        CollisionalRadiativeThermal::new(Temperature::new(1000.0), Temperature::new(50000.0));
    thermal.electron_power = Power::new(2e7);
    thermal.heavy_power = Power::new(1e7);
    let momentum_cross_section = 1e-24;
    let momentum = ElectronCrossSectionV1::tabulated(
        id("H"),
        ElectronCollisionKindV1::ElasticMomentumTransfer,
        Provenance {
            source: "manufactured constant neutral momentum source for reacting transport oracle"
                .into(),
            license: "analytic test fixture".into(),
            artifact: None,
        },
        vec![
            ElectronCrossSectionPointV1 {
                energy: Energy::new(0.0),
                cross_section: Area::new(momentum_cross_section),
            },
            ElectronCrossSectionPointV1 {
                energy: Energy::new(500.0 * EV),
                cross_section: Area::new(momentum_cross_section),
            },
        ],
        Area::new(momentum_cross_section),
    )
    .unwrap();
    thermal.collisions = Some(ElectronHeavyCollisionSources {
        momentum_sources: vec![momentum],
        coulomb: None,
        // This manufactured neutral-only closure covers this trajectory's
        // charged fraction (<0.2); an unrestricted plasma requires Coulomb data.
        maximum_charged_mole_fraction: 0.5,
        maximum_target_thermal_speed_ratio: 0.1,
        rate_controls: MaxwellianRateControlsV1 {
            relative_tolerance: 1e-10,
            maximum_subdivisions: 4096,
        },
    });
    let energies = [16.0 * EV, 20.0 * EV];
    let rates = [400.0, 600.0];
    let duration = 2e-5;
    let result = advance_radiation_transport(
        &[RadiationTransportCell::new(&initial, &source).with_thermal(&thermal); 2],
        &beams(energies, rates),
        area(),
        Time::new(duration),
        controls(1e-7, duration),
        || false,
    )
    .unwrap();
    let expected = reference(
        [&initial; 2],
        duration,
        1.0,
        energies,
        rates,
        [[2e7, 1e7]; 2],
        momentum_cross_section,
    );
    for (c, cell) in result.cells().iter().enumerate() {
        close(cell.state().amounts()[1].value(), expected[c][0], 3e-7, 0.0);
        close(
            cell.state().electron_temperature().value(),
            expected[c][1] / (CV * expected[c][0]),
            0.05,
            0.0,
        );
        close(
            cell.state().heavy_temperature().value(),
            expected[c][7] / CV,
            0.002,
            0.0,
        );
        let expected_exchange =
            expected[c][7] - CV * initial.heavy_temperature().value() - 1e7 * duration;
        close(
            cell.balance().electron_to_heavy_neutral_energy_j,
            expected_exchange,
            0.025,
            0.0,
        );
        assert!(cell.balance().electron_to_heavy_neutral_energy_j > 100.0);
        assert_eq!(cell.balance().electron_to_heavy_coulomb_energy_j, 0.0);
        assert_eq!(cell.thermal_inputs(), Some(&thermal));
    }
    close(
        result.balance().electron_deposited_energy_j,
        2.0 * 2e7 * duration,
        1e-9,
        0.0,
    );
    close(
        result.balance().heavy_deposited_energy_j,
        2.0 * 1e7 * duration,
        1e-9,
        0.0,
    );
    check_energy_and_photons(&result, 1e-6);
}

#[test]
fn dark_and_identity_limits_keep_source_admission_and_atomic_budget_refusals() {
    let (initial, source) = fixture(0.1, 10000.0, 1.0, 0.0, 0.0);
    let before = initial.clone();
    let duration = 2e-5;
    let solve = controls(1e-7, duration);
    let cell = RadiationTransportCell::new(&initial, &source);
    let dark =
        advance_radiation_transport(&[cell], &[], area(), Time::new(duration), solve, || false)
            .unwrap();
    let expected = advance_collisional_radiative(
        &initial,
        &source,
        &[],
        Time::new(duration),
        solve.chemistry,
        || false,
    )
    .unwrap();
    close(
        dark.cells()[0].state().amounts()[1].value(),
        expected.state().amounts()[1].value(),
        3e-7,
        0.0,
    );
    close(
        dark.cells()[0].state().electron_temperature().value(),
        expected.state().electron_temperature().value(),
        0.05,
        0.0,
    );
    assert!(dark.balance().emitted_photons > 0.0);
    let photons = beams([16.0 * EV, 20.0 * EV], [400.0, 600.0]);
    let clear =
        advance_radiation_transport(&[], &photons, area(), Time::new(duration), solve, || false)
            .unwrap();
    assert!(clear.cells().is_empty());
    for b in clear.group_balances() {
        assert_eq!(b.incident_photons, b.transmitted_photons);
    }
    let zero =
        advance_radiation_transport(&[cell], &photons, area(), Time::new(0.0), solve, || false)
            .unwrap();
    assert_eq!(zero.cells()[0].state(), &initial);
    assert_eq!(zero.balance().accepted_steps, 0);
    let (neutral, _) = fixture(0.0, 10000.0, 1.0, 0.0, 0.0);
    let unseeded_dark = advance_radiation_transport(
        &[RadiationTransportCell::new(&neutral, &source)],
        &[],
        area(),
        Time::new(duration),
        solve,
        || false,
    )
    .unwrap();
    assert_eq!(unseeded_dark.cells()[0].state(), &neutral);
    let mut limited = solve;
    limited.maximum_cell_groups = 3;
    let mut polls = 0;
    assert!(matches!(
        advance_radiation_transport(
            &[cell; 2],
            &photons,
            area(),
            Time::new(duration),
            limited,
            || {
                polls += 1;
                false
            }
        ),
        Err(RadiationTransportError::CellGroupLimit { .. })
    ));
    assert_eq!(polls, 1);
    limited = solve;
    limited.chemistry.maximum_steps = 1;
    limited.chemistry.initial_step_s = 1e-10;
    assert_eq!(
        advance_radiation_transport(
            &[cell; 2],
            &photons,
            area(),
            Time::new(duration),
            limited,
            || false
        ),
        Err(RadiationTransportError::StepBudget)
    );
    let mut completed_polls = 0;
    advance_radiation_transport(
        &[cell; 2],
        &photons,
        area(),
        Time::new(duration),
        solve,
        || {
            completed_polls += 1;
            false
        },
    )
    .unwrap();
    polls = 0;
    assert_eq!(
        advance_radiation_transport(
            &[cell; 2],
            &photons,
            area(),
            Time::new(duration),
            solve,
            || {
                polls += 1;
                polls >= completed_polls - 1
            }
        ),
        Err(RadiationTransportError::Cancelled)
    );
    assert_eq!(initial, before);
}
