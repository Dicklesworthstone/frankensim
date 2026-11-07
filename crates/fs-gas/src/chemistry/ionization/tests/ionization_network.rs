use super::{EV, KB, ME, MP, NA, R, id, species};
use crate::chemistry::collisional_radiative::CollisionalRadiativeSolve;
use crate::chemistry::ionization_network::*;
use fs_matdb::Provenance;
use fs_qty::{Amount, Area, Energy, Temperature, Time, Volume};
use fs_thermochem::electron_collision::{
    ElectronCollisionKindV1, ElectronCrossSectionPointV1, ElectronCrossSectionV1,
    MaxwellianRateControlsV1,
};
use fs_thermochem::equilibrium::composition::{FrozenIonizedGasStateV1, IdealGasEquilibriumV1};
use fs_thermochem::kinetics::ElectronImpactRateV1;
use fs_thermochem::kinetics::atomic_ionization::photoionization::{
    AtomicPhotoionizationV1, PhotoionizationCrossSectionPointV1, PhotonNumberFluxV1,
};
use fs_thermochem::kinetics::atomic_ionization::radiative_recombination::AtomicRadiativeRecombinationV1;
use fs_thermochem::kinetics::atomic_ionization::{AtomicIonizationV1, AtomicStatisticalWeightsV1};
use fs_thermochem::{ChargeVector, ElementId, ElementalMatrix, ReactionId};

const CV: f64 = 1.5 * R;
const ALPHA_MASS: f64 = 6.644_657_335_7e-27;
// NIST ASD 5.12 H I, He I and He II ionization energies; these level gaps
// do not make the manufactured impact/photo cross sections empirical data.
// https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
// https://physics.nist.gov/cgi-bin/ASD/ie.pl?spectra=He&units=1
const GAPS_EV: [f64; 3] = [13.598_434_599_702, 24.587_389_011, 54.417_765_528_2];
const SIGMA: [f64; 3] = [1e-22, 0.8e-22, 0.5e-22];
const IMPACT: [f64; 3] = [1.0e-18, 0.8e-18, 0.5e-18];
const EDGES: [[usize; 2]; 3] = [[0, 1], [2, 3], [3, 4]];
const WEIGHTS: [[u32; 2]; 3] = [[2, 1], [1, 2], [2, 1]];
const NAMES: [&str; 6] = ["H", "H+", "He", "He+", "He++", "e"];
const MASS: [f64; 6] = [
    MP + ME,
    MP,
    ALPHA_MASS + 2.0 * ME,
    ALPHA_MASS + ME,
    ALPHA_MASS,
    ME,
];

fn fixture(
    neutral: bool,
    shifts: [f64; 3],
) -> (FrozenIonizedGasStateV1, Vec<AtomicRadiativeRecombinationV1>) {
    let ids = NAMES.map(id).to_vec();
    let [h_shift, he_shift, q_shift] = shifts;
    let gaps = GAPS_EV.map(|e| e * EV * NA);
    let formations = [
        h_shift,
        h_shift + gaps[0] + q_shift,
        he_shift,
        he_shift + gaps[1] + q_shift,
        he_shift + gaps[1] + gaps[2] + 2.0 * q_shift,
        -q_shift,
    ];
    let bank = IdealGasEquilibriumV1::new(
        ElementalMatrix::new(
            vec![ElementId::new("H").unwrap(), ElementId::new("He").unwrap()],
            ids.clone(),
            vec![vec![1, 1, 0, 0, 0, 0], vec![0, 0, 1, 1, 1, 0]],
        )
        .unwrap(),
        ChargeVector::new(ids, vec![0, 1, 0, 1, 2, -1]).unwrap(),
        (0..6)
            .map(|i| {
                species(
                    NAMES[i],
                    MASS[i],
                    [2.0, 1.0, 1.0, 2.0, 1.0, 2.0][i],
                    formations[i],
                )
            })
            .collect(),
    )
    .unwrap();
    let provenance = Provenance {
        source: "Manufactured ramp-impact/constant-photo cross sections on source-bound H/He ground-state levels; not measured cross sections".into(),
        license: "analytic test fixture".into(), artifact: None,
    };
    let mut sources = Vec::new();
    for (edge, [lower, upper]) in EDGES.into_iter().enumerate() {
        // Use the actual retained source gap, including authored reference
        // roundoff, rather than substituting the nominal tabulated threshold.
        let formed = bank.models()[upper].regions()[0].coefficients()[7] * R
            - bank.models()[lower].regions()[0].coefficients()[7] * R
            + bank.models()[5].regions()[0].coefficients()[7] * R;
        let gap = formed / NA;
        let impact = ElectronCrossSectionV1::tabulated(
            id(NAMES[lower]),
            ElectronCollisionKindV1::Ionization,
            provenance.clone(),
            vec![
                ElectronCrossSectionPointV1 {
                    energy: Energy::new(gap),
                    cross_section: Area::new(0.0),
                },
                ElectronCrossSectionPointV1 {
                    energy: Energy::new(gap + 1000.0 * EV),
                    cross_section: Area::new(IMPACT[edge]),
                },
            ],
            Area::new(IMPACT[edge]),
        )
        .unwrap();
        let channel = AtomicIonizationV1::new(
            bank.clone(),
            ReactionId::new(format!("network-edge-{edge}")).unwrap(),
            id(NAMES[upper]),
            ElectronImpactRateV1::new(
                id("e"),
                impact,
                MaxwellianRateControlsV1 {
                    relative_tolerance: 1e-10,
                    maximum_subdivisions: 4096,
                },
                0.1,
            )
            .unwrap(),
            AtomicStatisticalWeightsV1 {
                neutral: WEIGHTS[edge][0],
                ion: WEIGHTS[edge][1],
            },
            1e-12,
        )
        .unwrap();
        let photon = AtomicPhotoionizationV1::new(
            channel,
            provenance.clone(),
            vec![
                PhotoionizationCrossSectionPointV1 {
                    photon_energy: Energy::new(gap),
                    cross_section: Area::new(SIGMA[edge]),
                },
                PhotoionizationCrossSectionPointV1 {
                    photon_energy: Energy::new(gap + 1000.0 * EV),
                    cross_section: Area::new(SIGMA[edge]),
                },
            ],
        )
        .unwrap();
        sources.push(
            AtomicRadiativeRecombinationV1::new(photon, Area::new(SIGMA[edge]), provenance.clone())
                .unwrap(),
        );
    }
    let amounts = if neutral {
        [1.0, 0.0, 0.7, 0.0, 0.0, 0.0]
    } else {
        [0.7, 0.3, 0.4, 0.2, 0.1, 0.7]
    };
    let state = FrozenIonizedGasStateV1::new(
        &bank,
        &amounts.map(Amount::new),
        Volume::new(20.0),
        Temperature::new(6000.0),
        Temperature::new(20000.0),
        1e-12,
        || false,
    )
    .unwrap();
    (state, sources)
}

fn field() -> [PhotonNumberFluxV1; 3] {
    [14.0, 26.0, 60.0].map(|e| PhotonNumberFluxV1 {
        photon_energy: Energy::new(e * EV),
        number_flux_per_m2_s: 2e24,
    })
}
fn controls(tolerance: f64) -> IonizationNetworkSolve {
    let mut local =
        CollisionalRadiativeSolve::new(Temperature::new(1000.0), Temperature::new(50000.0));
    local.relative_tolerance = tolerance;
    local.initial_step_s = 1e-6;
    local.energy_tolerance_j = 2e-6;
    IonizationNetworkSolve::new(local)
}

// Independent simultaneous ODE. It uses analytic infinite-domain Maxwellian
// integrals, not any production rate evaluator, quadrature or split stage.
// The finite production tables end more than 230 kT above threshold here;
// their omitted moments are negligible compared with the test tolerance.
// State: six populations, electron sensible energy, six accounts per edge.
fn rhs(y: &[f64; 25], volume: f64, photons: &[PhotonNumberFluxV1]) -> [f64; 25] {
    let mut out = [0.0; 25];
    let mut birth = 0.0;
    let mut birth_energy = 0.0;
    for (edge, [lower, _]) in EDGES.into_iter().enumerate() {
        for photon in photons {
            if photon.photon_energy.value() >= GAPS_EV[edge] * EV {
                let q = y[lower] * SIGMA[edge] * photon.number_flux_per_m2_s;
                birth += q;
                birth_energy += q * NA * (photon.photon_energy.value() - GAPS_EV[edge] * EV);
            }
        }
    }
    let te = if y[5] > 0.0 {
        y[6] / (CV * y[5])
    } else if birth > 0.0 {
        birth_energy / (CV * birth)
    } else {
        20000.0
    };
    let theta = KB * te;
    let h = 6.626_070_15e-34_f64;
    let light = 299_792_458.0_f64;
    for (edge, [lower, upper]) in EDGES.into_iter().enumerate() {
        let gap = GAPS_EV[edge] * EV;
        let g = f64::from(WEIGHTS[edge][1]) / f64::from(WEIGHTS[edge][0]);
        let coefficient = (8.0 / (core::f64::consts::PI * ME)).sqrt()
            * (IMPACT[edge] / (1000.0 * EV))
            * theta.sqrt()
            * (gap + 2.0 * theta)
            * (-gap / theta).exp();
        let saha = 2.0
            * g
            * (2.0 * core::f64::consts::PI * ME * theta / (h * h)).powf(1.5)
            * (MASS[upper] / MASS[lower]).powf(1.5)
            / NA
            * (-gap / theta).exp();
        let collision = NA * coefficient * y[lower] * y[5] / volume
            - (NA * coefficient / saha) * y[upper] * y[5] * y[5] / (volume * volume);
        let polynomial = gap * gap + 2.0 * gap * theta + 2.0 * theta * theta;
        let alpha = (1.0 / g) * (2.0 / core::f64::consts::PI).sqrt()
            / (ME.powf(1.5) * light * light)
            * SIGMA[edge]
            * polynomial
            / theta.sqrt();
        let captured_energy =
            theta * (gap * gap + 4.0 * gap * theta + 6.0 * theta * theta) / polynomial;
        let capture = NA * alpha * y[upper] * y[5] / volume;
        let mut photo = 0.0;
        let mut absorbed = 0.0;
        for photon in photons {
            if photon.photon_energy.value() >= gap {
                let q = y[lower] * SIGMA[edge] * photon.number_flux_per_m2_s;
                photo += q;
                absorbed += q * NA * photon.photon_energy.value();
            }
        }
        let net = collision + photo - capture;
        let captured = NA * captured_energy * capture;
        out[lower] -= net;
        out[upper] += net;
        out[5] += net;
        out[6] += absorbed - NA * gap * (photo + collision) - captured;
        for (slot, value) in [net, photo, collision, capture, absorbed, captured]
            .into_iter()
            .enumerate()
        {
            out[7 + 6 * edge + slot] = value;
        }
    }
    out
}

fn reference(
    initial: &FrozenIonizedGasStateV1,
    photons: &[PhotonNumberFluxV1],
    duration: f64,
) -> [f64; 25] {
    let mut y = [0.0; 25];
    for (target, amount) in y.iter_mut().zip(initial.amounts()) {
        *target = amount.value();
    }
    y[6] = CV * y[5] * initial.electron_temperature().value();
    let steps = 4096;
    let dt = duration / steps as f64;
    for _ in 0..steps {
        let k1 = rhs(&y, initial.volume().value(), photons);
        let k2 = rhs(
            &core::array::from_fn(|i| y[i] + 0.5 * dt * k1[i]),
            initial.volume().value(),
            photons,
        );
        let k3 = rhs(
            &core::array::from_fn(|i| y[i] + 0.5 * dt * k2[i]),
            initial.volume().value(),
            photons,
        );
        let k4 = rhs(
            &core::array::from_fn(|i| y[i] + dt * k3[i]),
            initial.volume().value(),
            photons,
        );
        for i in 0..25 {
            y[i] += (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }
    }
    y
}

fn assert_accounts(result: &IonizationNetworkResult) {
    let amounts = result.state().amounts();
    assert!((amounts[0].value() + amounts[1].value() - 1.0).abs() < 3e-13);
    assert!((amounts[2].value() + amounts[3].value() + amounts[4].value() - 0.7).abs() < 3e-13);
    assert!(
        (amounts[1].value() + amounts[3].value() + 2.0 * amounts[4].value() - amounts[5].value())
            .abs()
            < 3e-13
    );
    assert_eq!(
        result.state().heavy_temperature(),
        result.initial_state().heavy_temperature()
    );
    assert_eq!(result.state().volume(), result.initial_state().volume());
    assert_eq!(
        result.state().thermodynamics(),
        result.initial_state().thermodynamics()
    );
    assert!(result.balance().energy_residual_j.abs() <= result.controls().local.energy_tolerance_j);
    assert!(result.balance().relative_charge_residual < 1e-12);
    let sum = |f: fn(&IonizationChannelBalance) -> f64| {
        result.channel_balances().iter().map(f).sum::<f64>()
    };
    assert!(
        (sum(|c| c.photoionization_extent_mol) * NA / result.balance().absorbed_photons - 1.0)
            .abs()
            < 1e-13
    );
    assert!(
        (sum(|c| c.radiative_recombination_extent_mol) * NA / result.balance().emitted_photons
            - 1.0)
            .abs()
            < 1e-13
    );
}

#[test]
fn g1_shared_electrons_couple_adjacent_stages_and_elements_against_simultaneous_ode() {
    let (initial, sources) = fixture(false, [0.0; 3]);
    let photons = field();
    let time = 1e-4;
    let expected = reference(&initial, &photons, time);
    let result = advance_ionization_network(
        &initial,
        &sources,
        &photons,
        Time::new(time),
        controls(2e-8),
        || false,
    )
    .unwrap();
    for (actual, expected) in result.state().amounts().iter().zip(&expected[..6]) {
        assert!(
            (actual.value() - expected).abs() < 2e-7,
            "{} versus {expected}",
            actual.value()
        );
    }
    let te = expected[6] / (CV * expected[5]);
    assert!(
        (result.state().electron_temperature().value() - te).abs() < 0.03,
        "{} versus {te}",
        result.state().electron_temperature().value()
    );
    for (edge, account) in result.channel_balances().iter().enumerate() {
        for (actual, expected) in [
            account.net_extent_mol,
            account.photoionization_extent_mol,
            account.collisional_extent_mol,
            account.radiative_recombination_extent_mol,
        ]
        .into_iter()
        .zip(&expected[7 + 6 * edge..11 + 6 * edge])
        {
            assert!(
                (actual - expected).abs() < 3e-7,
                "edge {edge}: {actual} versus {expected}"
            );
        }
        assert!((account.absorbed_radiation_energy_j - expected[11 + 6 * edge]).abs() < 0.8);
        let emitted = GAPS_EV[edge] * EV * NA * expected[10 + 6 * edge] + expected[12 + 6 * edge];
        assert!((account.emitted_radiation_energy_j - emitted).abs() < 0.8);
    }
    assert!(result.state().amounts()[4].value() != initial.amounts()[4].value());
    assert!(
        result
            .channel_balances()
            .iter()
            .all(|c| c.photoionization_extent_mol > 0.0
                && c.radiative_recombination_extent_mol > 0.0)
    );
    assert_accounts(&result);
    let mut reversed = sources.clone();
    reversed.reverse();
    assert_eq!(
        result,
        advance_ionization_network(
            &initial,
            &reversed,
            &photons,
            Time::new(time),
            controls(2e-8),
            || false
        )
        .unwrap()
    );
}

#[test]
fn g3_network_refinement_is_second_order_and_reference_shifts_preserve_dynamics() {
    let (initial, sources) = fixture(false, [0.0; 3]);
    let photons = field();
    let time = 1e-4;
    // Initial H impact has a 1.18e5/s birth coefficient and a roughly 1.4us
    // electron cooling timescale. Resolve that same dense initial transient
    // on a 1us interval; the separate 100us adaptive oracle above is unchanged.
    let resolved_time = 1e-6;
    let expected = reference(&initial, &photons, resolved_time);
    let mut errors = Vec::new();
    for steps in [16, 32, 64] {
        let mut solve = controls(0.1);
        solve.local.absolute_amount_tolerance_mol = 1.0;
        solve.local.temperature_tolerance_k = 1e6;
        solve.local.initial_step_s = resolved_time / steps as f64;
        solve.local.maximum_step_s = solve.local.initial_step_s;
        let result = advance_ionization_network(
            &initial,
            &sources,
            &photons,
            Time::new(resolved_time),
            solve,
            || false,
        )
        .unwrap();
        assert_eq!(result.balance().rejected_steps, 0);
        errors.push(
            result
                .state()
                .amounts()
                .iter()
                .zip(&expected[..6])
                .map(|(n, y)| (n.value() - y).abs())
                .fold(0.0_f64, f64::max),
        );
        assert_accounts(&result);
    }
    for pair in errors.windows(2) {
        assert!(
            (3.6..4.4).contains(&(pair[0] / pair[1])),
            "errors {errors:?}"
        );
    }
    let base = advance_ionization_network(
        &initial,
        &sources,
        &photons,
        Time::new(time),
        controls(2e-6),
        || false,
    )
    .unwrap();
    let (shifted, shifted_sources) = fixture(false, [1e10, -2e10, 5e9]);
    let shifted = advance_ionization_network(
        &shifted,
        &shifted_sources,
        &photons,
        Time::new(time),
        controls(2e-6),
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
    assert_accounts(&shifted);
}

#[test]
fn g1_exact_neutral_photo_birth_feeds_the_shared_pool_without_a_seed() {
    let (initial, sources) = fixture(true, [0.0; 3]);
    let photons = [PhotonNumberFluxV1 {
        photon_energy: Energy::new(14.0 * EV),
        number_flux_per_m2_s: 2e24,
    }];
    let time = 2e-5;
    let expected = reference(&initial, &photons, time);
    let result = advance_ionization_network(
        &initial,
        &sources,
        &photons,
        Time::new(time),
        controls(1e-6),
        || false,
    )
    .unwrap();
    assert_eq!(initial.amounts()[5].value(), 0.0);
    assert!(result.state().amounts()[5].value() > 0.0);
    for (actual, expected) in result.state().amounts().iter().zip(&expected[..6]) {
        assert!((actual.value() - expected).abs() < 5e-8);
    }
    assert!(
        (result.state().electron_temperature().value() - expected[6] / (CV * expected[5])).abs()
            < 0.02
    );
    assert_accounts(&result);
    let dark = advance_ionization_network(
        &initial,
        &sources,
        &[],
        Time::new(1e6),
        controls(1e-6),
        || false,
    )
    .unwrap();
    assert_eq!(dark.state(), &initial);
    assert_eq!(dark.balance().accepted_steps, 0);
    assert_eq!(dark.balance().rate_evaluations, 0);
}

#[test]
fn network_admission_work_cancellation_and_zero_time_are_transactional() {
    let (initial, sources) = fixture(false, [0.0; 3]);
    let unchanged = initial.clone();
    let photons = field();
    let zero = advance_ionization_network(
        &initial,
        &sources,
        &photons,
        Time::new(0.0),
        controls(1e-6),
        || false,
    )
    .unwrap();
    assert_eq!(zero.state(), &initial);
    assert_eq!(zero.balance().rate_evaluations, 0);
    let empty = advance_ionization_network(
        &initial,
        &[],
        &photons,
        Time::new(1.0),
        controls(1e-6),
        || false,
    )
    .unwrap();
    assert_eq!(empty.state(), &initial);
    assert!(empty.sources().is_empty());
    let mut duplicate = sources.clone();
    duplicate.push(sources[0].clone());
    assert!(matches!(
        advance_ionization_network(
            &initial,
            &duplicate,
            &photons,
            Time::new(1e-4),
            controls(1e-6),
            || false
        ),
        Err(IonizationNetworkError::Invalid(_))
    ));
    let (_, foreign) = fixture(false, [1e9, 0.0, 0.0]);
    assert!(matches!(
        advance_ionization_network(
            &initial,
            &foreign,
            &photons,
            Time::new(1e-4),
            controls(1e-6),
            || false
        ),
        Err(IonizationNetworkError::Invalid(_))
    ));
    let mut limited = controls(1e-6);
    limited.maximum_channels = 2;
    assert!(matches!(
        advance_ionization_network(
            &initial,
            &sources,
            &photons,
            Time::new(1e-4),
            limited,
            || false
        ),
        Err(IonizationNetworkError::Invalid(_))
    ));
    limited.maximum_channels = 3;
    limited.maximum_rate_evaluations = 1;
    assert!(matches!(
        advance_ionization_network(
            &initial,
            &sources,
            &photons,
            Time::new(1e-4),
            limited,
            || false
        ),
        Err(IonizationNetworkError::StepBudget)
    ));
    for stop in [1, 150] {
        let mut polls = 0;
        assert!(matches!(
            advance_ionization_network(
                &initial,
                &sources,
                &photons,
                Time::new(1e-4),
                controls(1e-6),
                || {
                    polls += 1;
                    polls >= stop
                }
            ),
            Err(IonizationNetworkError::Cancelled)
        ));
    }
    assert_eq!(initial, unchanged);
}
