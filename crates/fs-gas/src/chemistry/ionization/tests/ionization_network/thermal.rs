use super::*;
use crate::chemistry::collisional_radiative::CollisionalRadiativeThermal;
use crate::chemistry::energy::{ElectronHeavyCollisionSources, ElectronHeavyEnergyError};
use fs_qty::Power;
use fs_thermochem::coulomb::{CoulombLimitsV1, CoulombModelV1, CoulombScreeningV1};

const MOMENTUM: [f64; 2] = [4e-19, 7e-19];

fn parcel(
    neutral: bool,
    shifts: [f64; 3],
) -> (FrozenIonizedGasStateV1, Vec<AtomicRadiativeRecombinationV1>) {
    let (state, sources) = fixture(neutral, shifts);
    (
        FrozenIonizedGasStateV1::new(
            state.thermodynamics(),
            state.amounts(),
            Volume::new(10000.0),
            state.heavy_temperature(),
            state.electron_temperature(),
            1e-12,
            || false,
        )
        .unwrap(),
        sources,
    )
}

fn thermal(state: &FrozenIonizedGasStateV1, coulomb: bool) -> CollisionalRadiativeThermal {
    let mut inputs =
        CollisionalRadiativeThermal::new(Temperature::new(1000.0), Temperature::new(50000.0));
    inputs.collisions = Some(ElectronHeavyCollisionSources {
        momentum_sources: [0, 2]
            .into_iter()
            .zip(MOMENTUM)
            .map(|(index, sigma)| {
                ElectronCrossSectionV1::tabulated(
                    id(NAMES[index]),
                    ElectronCollisionKindV1::ElasticMomentumTransfer,
                    Provenance {
                        source: "manufactured constant H/He momentum sections".into(),
                        license: "analytic test fixture".into(),
                        artifact: None,
                    },
                    [0.0, 1000.0 * EV]
                        .map(|e| ElectronCrossSectionPointV1 {
                            energy: Energy::new(e),
                            cross_section: Area::new(sigma),
                        })
                        .to_vec(),
                    Area::new(sigma),
                )
                .unwrap()
            })
            .collect(),
        coulomb: coulomb.then(|| {
            CoulombModelV1::new(
                state,
                CoulombScreeningV1::StaticAllChargedDebye,
                CoulombLimitsV1 {
                    minimum_coulomb_logarithm: 4.0,
                    minimum_particles_per_debye_sphere: 10.0,
                    maximum_coulomb_coupling: 0.1,
                    maximum_thermal_phase_space_density: 0.1,
                    maximum_thermal_energy_rest_energy_ratio: 0.001,
                    maximum_quantum_to_classical_ratio: 0.5,
                    maximum_relative_drift_thermal_ratio: 0.1,
                },
                || false,
            )
            .unwrap()
        }),
        maximum_charged_mole_fraction: if coulomb { 1.0 } else { 0.9 },
        maximum_target_thermal_speed_ratio: 0.1,
        rate_controls: MaxwellianRateControlsV1 {
            relative_tolerance: 1e-10,
            maximum_subdivisions: 4096,
        },
    });
    inputs.electron_power = Power::new(2e8);
    inputs.heavy_power = Power::new(3e7);
    inputs
}

// Independent simultaneous ODE: the parent analytic impact/Milne equations
// are augmented by two original constant-section Maxwellian neutral moments
// and a direct SI finite-mass Landau electron/ion sum, including He++ charge².
// No production rate, Coulomb pair, quadrature or thermal time step is called.
// Banks (1966): https://ntrs.nasa.gov/citations/19660020172
// Finite-mass moments: https://farside.ph.utexas.edu/teaching/plasma/lectures/node42.html
// Slots 0..25 retain the parent chemical oracle; 25 is heavy sensible energy,
// 26 and 27 independently integrate neutral and Coulomb energy transfers.
fn coupled_rhs(
    y: &[f64; 28],
    volume: f64,
    photons: &[PhotonNumberFluxV1],
    inputs: &CollisionalRadiativeThermal,
) -> [f64; 28] {
    let chemistry = rhs(&core::array::from_fn(|i| y[i]), volume, photons);
    let mut out = [0.0; 28];
    out[..25].copy_from_slice(&chemistry);
    let te = y[6] / (CV * y[5]);
    let th = y[25] / (CV * y[..5].iter().sum::<f64>());
    let neutral_g = [0, 2]
        .into_iter()
        .zip(MOMENTUM)
        .map(|(i, sigma)| {
            3.0 * R
                * y[5]
                * (y[i] * NA / volume)
                * (ME / MASS[i])
                * (4.0 / 3.0)
                * sigma
                * (8.0 * KB * te / (core::f64::consts::PI * ME)).sqrt()
        })
        .sum::<f64>();
    let e = 1.602_176_634e-19_f64;
    let eps = 8.854_187_818_8e-12_f64;
    let ne = y[5] * NA / volume;
    let ionic_screen = (y[1] + y[3] + 4.0 * y[4]) * NA / (volume * th);
    let debye = (eps * KB / (e * e * (ne / te + ionic_screen))).sqrt();
    let coulomb_g = [(1, 1.0), (3, 1.0), (4, 2.0)]
        .into_iter()
        .map(|(i, charge)| {
            let variance = KB * (te / ME + th / MASS[i]);
            let reduced_mass = ME * MASS[i] / (ME + MASS[i]);
            let c = charge * e * e / (4.0 * core::f64::consts::PI * eps);
            let b90 = c / (reduced_mass * variance);
            4.0 * (2.0 * core::f64::consts::PI).sqrt()
                * KB
                * ne
                * (y[i] * NA)
                * c
                * c
                * (debye / b90).ln()
                / (ME * MASS[i] * variance.powf(1.5))
        })
        .sum::<f64>();
    let qn = neutral_g * (te - th);
    let qc = coulomb_g * (te - th);
    out[6] += inputs.electron_power.value() - qn - qc;
    out[25] = inputs.heavy_power.value() + qn + qc;
    out[26] = qn;
    out[27] = qc;
    out
}

fn coupled_reference(
    initial: &FrozenIonizedGasStateV1,
    photons: &[PhotonNumberFluxV1],
    inputs: &CollisionalRadiativeThermal,
    duration: f64,
) -> [f64; 28] {
    let mut y = [0.0; 28];
    for (target, amount) in y.iter_mut().zip(initial.amounts()) {
        *target = amount.value();
    }
    y[6] = CV * y[5] * initial.electron_temperature().value();
    y[25] = CV * y[..5].iter().sum::<f64>() * initial.heavy_temperature().value();
    let dt = duration / 8192.0;
    for _ in 0..8192 {
        let f = |v| coupled_rhs(&v, initial.volume().value(), photons, inputs);
        let a = f(y);
        let b = f(core::array::from_fn(|i| y[i] + 0.5 * dt * a[i]));
        let c = f(core::array::from_fn(|i| y[i] + 0.5 * dt * b[i]));
        let d = f(core::array::from_fn(|i| y[i] + dt * c[i]));
        for i in 0..y.len() {
            y[i] += (dt / 6.0) * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]);
        }
    }
    y
}

fn assert_thermal_accounts(result: &IonizationNetworkResult) {
    let b = result.balance();
    let inputs = result.thermal_inputs().unwrap();
    let time = result.elapsed_time().value();
    assert_eq!(
        b.electron_deposited_energy_j,
        inputs.electron_power.value() * time
    );
    assert_eq!(
        b.heavy_deposited_energy_j,
        inputs.heavy_power.value() * time
    );
    assert_eq!(result.state().volume(), result.initial_state().volume());
    assert_eq!(
        result.state().thermodynamics(),
        result.initial_state().thermodynamics()
    );
    let n = result.state().amounts();
    assert!((n[0].value() + n[1].value() - 1.0).abs() < 2e-13);
    assert!((n[2].value() + n[3].value() + n[4].value() - 0.7).abs() < 2e-13);
    assert!(b.relative_charge_residual < 1e-12);
    for residual in [
        b.energy_residual_j,
        b.electron_energy_residual_j,
        b.heavy_energy_residual_j,
        b.electron_to_heavy_energy_j
            - b.electron_to_heavy_neutral_energy_j
            - b.electron_to_heavy_coulomb_energy_j,
    ] {
        assert!(
            residual.abs() <= result.controls().local.energy_tolerance_j,
            "{b:?}"
        );
    }
    assert!(b.maximum_normalized_local_error <= 1.0);
}

#[test]
fn g1_network_thermal_pools_and_charge_dependent_exchange_match_simultaneous_ode() {
    let (state, sources) = parcel(false, [0.0; 3]);
    let inputs = thermal(&state, true);
    let photons = field();
    let time = 2e-5;
    let reference = coupled_reference(&state, &photons, &inputs, time);
    let result = advance_ionization_network_with_thermal(
        &state,
        &sources,
        &photons,
        &inputs,
        Time::new(time),
        controls(2e-8),
        || false,
    )
    .unwrap();
    for (actual, expected) in result.state().amounts().iter().zip(&reference[..6]) {
        assert!(
            (actual.value() - expected).abs() < 1e-7,
            "{actual:?} versus {expected}"
        );
    }
    let te = reference[6] / (CV * reference[5]);
    let th = reference[25] / (CV * reference[..5].iter().sum::<f64>());
    assert!(
        (result.state().electron_temperature().value() - te).abs() < 0.015,
        "Te {te} {:?}",
        result.state().electron_temperature()
    );
    assert!(
        (result.state().heavy_temperature().value() - th).abs() < 0.006,
        "Th {th} {:?}",
        result.state().heavy_temperature()
    );
    assert!((result.balance().electron_to_heavy_neutral_energy_j - reference[26]).abs() < 0.05);
    assert!((result.balance().electron_to_heavy_coulomb_energy_j - reference[27]).abs() < 0.05);
    assert!(result.state().heavy_temperature().value() > state.heavy_temperature().value());
    assert!(result.balance().electron_to_heavy_neutral_energy_j > 0.0);
    assert!(result.balance().electron_to_heavy_coulomb_energy_j > 0.0);
    assert!(
        result
            .channel_balances()
            .iter()
            .all(|b| b.radiative_recombination_extent_mol > 0.0
                && b.photoionization_extent_mol > 0.0)
    );
    assert_eq!(result.thermal_inputs(), Some(&inputs));
    assert_thermal_accounts(&result);
    let mut reversed = sources.clone();
    reversed.reverse();
    assert_eq!(
        result,
        advance_ionization_network_with_thermal(
            &state,
            &reversed,
            &photons,
            &inputs,
            Time::new(time),
            controls(2e-8),
            || false,
        )
        .unwrap()
    );
}

#[test]
fn g3_network_thermal_refinement_and_formation_reference_invariance() {
    let (state, sources) = parcel(false, [0.0; 3]);
    let inputs = thermal(&state, true);
    let photons = field();
    let time = 2e-5;
    let reference = coupled_reference(&state, &photons, &inputs, time);
    let expected = [
        reference[6] / (CV * reference[5]),
        reference[25] / (CV * 1.7),
    ];
    let mut errors = Vec::new();
    for steps in [16, 32, 64] {
        let mut solve = controls(0.1);
        solve.local.temperature_tolerance_k = 1e6;
        solve.local.absolute_amount_tolerance_mol = 1.0;
        solve.local.initial_step_s = time / steps as f64;
        solve.local.maximum_step_s = solve.local.initial_step_s;
        let result = advance_ionization_network_with_thermal(
            &state,
            &sources,
            &photons,
            &inputs,
            Time::new(time),
            solve,
            || false,
        )
        .unwrap();
        assert_eq!(result.balance().rejected_steps, 0);
        errors.push(
            [
                result.state().electron_temperature().value(),
                result.state().heavy_temperature().value(),
            ]
            .into_iter()
            .zip(expected)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max),
        );
        assert_thermal_accounts(&result);
    }
    for pair in errors.windows(2) {
        assert!((3.6..4.4).contains(&(pair[0] / pair[1])), "{errors:?}");
    }
    let base = advance_ionization_network_with_thermal(
        &state,
        &sources,
        &photons,
        &inputs,
        Time::new(time),
        controls(1e-6),
        || false,
    )
    .unwrap();
    let (shifted, sources) = parcel(false, [1e10, -2e10, 5e9]);
    let inputs = thermal(&shifted, true);
    let shifted = advance_ionization_network_with_thermal(
        &shifted,
        &sources,
        &photons,
        &inputs,
        Time::new(time),
        controls(1e-6),
        || false,
    )
    .unwrap();
    for (a, b) in base.state().amounts().iter().zip(shifted.state().amounts()) {
        assert!((a.value() - b.value()).abs() < 1e-9);
    }
    assert!(
        (base.state().electron_temperature().value()
            - shifted.state().electron_temperature().value())
        .abs()
            < 2e-5
    );
    assert!(
        (base.state().heavy_temperature().value() - shifted.state().heavy_temperature().value())
            .abs()
            < 2e-5
    );
    assert_thermal_accounts(&shifted);
}

#[test]
fn g0_network_thermal_birth_complete_neutral_support_and_atomic_refusals() {
    let (neutral, sources) = parcel(true, [0.0; 3]);
    let mut inputs = thermal(&neutral, true);
    inputs.electron_power = Power::new(0.0);
    inputs.heavy_power = Power::new(0.0);
    let photons = [PhotonNumberFluxV1 {
        photon_energy: Energy::new(14.0 * EV),
        number_flux_per_m2_s: 2e24,
    }];
    let born = advance_ionization_network_with_thermal(
        &neutral,
        &sources,
        &photons,
        &inputs,
        Time::new(1e-4),
        controls(1e-7),
        || false,
    )
    .unwrap();
    assert!(born.state().amounts()[5].value() > 0.0);
    assert!(born.balance().electron_to_heavy_energy_j < 0.0);
    assert!(born.state().heavy_temperature().value() < neutral.heavy_temperature().value());
    assert_thermal_accounts(&born);
    let dark = advance_ionization_network_with_thermal(
        &neutral,
        &sources,
        &[],
        &inputs,
        Time::new(1000.0),
        controls(1e-6),
        || false,
    )
    .unwrap();
    assert_eq!(dark.state(), &neutral);
    assert_eq!(dark.balance().accepted_steps, 0);
    inputs.electron_power = Power::new(1.0);
    assert!(matches!(
        advance_ionization_network_with_thermal(
            &neutral,
            &sources,
            &photons,
            &inputs,
            Time::new(1.0),
            controls(1e-6),
            || false,
        ),
        Err(IonizationNetworkError::Invalid(_))
    ));

    let (state, sources) = parcel(false, [0.0; 3]);
    let state = FrozenIonizedGasStateV1::new(
        state.thermodynamics(),
        &[0.7, 0.3, 0.0, 0.6, 0.1, 1.1].map(Amount::new),
        state.volume(),
        state.heavy_temperature(),
        state.electron_temperature(),
        1e-12,
        || false,
    )
    .unwrap();
    let mut inputs = thermal(&state, true);
    let original = state.clone();
    inputs.collisions.as_mut().unwrap().momentum_sources.pop();
    assert!(matches!(advance_ionization_network_with_thermal(
        &state, &sources, &[], &inputs, Time::new(0.0), controls(1e-6), || false,
    ), Err(IonizationNetworkError::Thermal(ElectronHeavyEnergyError::MissingSource(id))) if id.as_str() == "He"));
    let inputs = thermal(&state, true);
    let mut limited = controls(1e-6);
    limited.maximum_rate_evaluations = 1;
    assert!(matches!(
        advance_ionization_network_with_thermal(
            &state,
            &sources,
            &[],
            &inputs,
            Time::new(1e-4),
            limited,
            || false,
        ),
        Err(IonizationNetworkError::StepBudget)
    ));
    let mut calls = 0;
    assert!(matches!(
        advance_ionization_network_with_thermal(
            &state,
            &sources,
            &[],
            &inputs,
            Time::new(1e-4),
            controls(1e-6),
            || {
                calls += 1;
                calls > 500
            },
        ),
        Err(IonizationNetworkError::Cancelled)
    ));
    assert_eq!(state, original);
}
