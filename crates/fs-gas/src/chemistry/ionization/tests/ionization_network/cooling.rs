use super::*;
use crate::chemistry::diffuse_recombination::DiffuseRecombinationSolve;
use crate::chemistry::ionization_network::transport::diffuse::{
    NetworkDiffuseSolve, advance_network_diffuse_recombination,
};
use crate::chemistry::ionization_network::transport::{
    NetworkRadiationCell, NetworkRadiationError, NetworkRadiationSolve,
    advance_network_radiation_transport,
};
use fs_thermochem::plasma_radiation::{
    BremsstrahlungControlsV1, BremsstrahlungGauntFactorV1, ElectronIonBremsstrahlungV1,
};

fn cooling(
    state: &FrozenIonizedGasStateV1,
    law: BremsstrahlungGauntFactorV1,
) -> ElectronIonBremsstrahlungV1 {
    ElectronIonBremsstrahlungV1::new(
        state.thermodynamics().clone(),
        BremsstrahlungControlsV1 {
            frequency_integrated_gaunt_factor: law,
            minimum_electron_temperature_k: 1000.0,
            maximum_electron_temperature_k: 50000.0,
            maximum_electron_thermal_rest_energy_ratio: 0.01,
            maximum_ion_to_electron_thermal_variance_ratio: 0.01,
            maximum_electron_thermal_phase_space_density: 0.01,
        },
    )
    .unwrap()
}

// Original admitted level/caloric banks with genuinely zero authored impact
// and photon cross sections isolate the continuum cooling integration.
fn silent_sources(
    state: &FrozenIonizedGasStateV1,
    originals: &[AtomicRadiativeRecombinationV1],
) -> Vec<AtomicRadiativeRecombinationV1> {
    originals
        .iter()
        .enumerate()
        .map(|(edge, original)| {
            let gap = original.source().channel().ionization_energy_j_mol() / NA;
            let [lower, upper] = EDGES[edge];
            let provenance = Provenance {
                source: "Manufactured zero atomic rates isolating original free-free cooling"
                    .into(),
                license: "analytic test fixture".into(),
                artifact: None,
            };
            let impact = ElectronCrossSectionV1::tabulated(
                id(NAMES[lower]),
                ElectronCollisionKindV1::Ionization,
                provenance.clone(),
                [gap, gap + 1000.0 * EV]
                    .map(|e| ElectronCrossSectionPointV1 {
                        energy: Energy::new(e),
                        cross_section: Area::new(0.0),
                    })
                    .to_vec(),
                Area::new(0.0),
            )
            .unwrap();
            let channel = AtomicIonizationV1::new(
                state.thermodynamics().clone(),
                ReactionId::new(format!("silent-edge-{edge}")).unwrap(),
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
            let photo = AtomicPhotoionizationV1::new(
                channel,
                provenance.clone(),
                [gap, gap + 1000.0 * EV]
                    .map(|e| PhotoionizationCrossSectionPointV1 {
                        photon_energy: Energy::new(e),
                        cross_section: Area::new(0.0),
                    })
                    .to_vec(),
            )
            .unwrap();
            AtomicRadiativeRecombinationV1::new(photo, Area::new(0.0), provenance).unwrap()
        })
        .collect()
}

fn near(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * expected.abs(),
        "{actual} versus {expected}"
    );
}

// Independent published bolometric normalization and Table2 polynomial from
// Itoh et al., arXiv:astro-ph/0011044, Eq.11–14. No production cooling stage,
// rate evaluator or source slope participates in this integration oracle.
fn continuum_power(amounts: &[f64], volume: f64, te: f64, law: BremsstrahlungGauntFactorV1) -> f64 {
    let gaunt = |z: f64| match law {
        BremsstrahlungGauntFactorV1::Constant(g) => g,
        BremsstrahlungGauntFactorV1::ItohNonrelativistic => {
            let x = ((z * z * 1.579e5 / te).log10() + 0.5) / 2.5;
            [
                1.42752, 0.170677, -0.813616, -0.488059, 1.33141, 0.727788, -1.50917, -0.599317,
                0.986146, 0.199586, -0.272392,
            ]
            .iter()
            .enumerate()
            .map(|(i, b)| b * x.powi(i as i32))
            .sum::<f64>()
        }
    };
    let ions = (amounts[1] + amounts[3]) * gaunt(1.0) + 4.0 * amounts[4] * gaunt(2.0);
    1.426e-40 * NA * NA * amounts[5] * ions / volume * te.sqrt()
}

fn rk4<const N: usize>(
    mut y: [f64; N],
    duration: f64,
    f: impl Fn(&[f64; N]) -> [f64; N],
) -> [f64; N] {
    let dt = duration / 4096.0;
    for _ in 0..4096 {
        let a = f(&y);
        let b = f(&core::array::from_fn(|i| y[i] + 0.5 * dt * a[i]));
        let c = f(&core::array::from_fn(|i| y[i] + 0.5 * dt * b[i]));
        let d = f(&core::array::from_fn(|i| y[i] + dt * c[i]));
        for i in 0..N {
            y[i] += dt / 6.0 * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]);
        }
    }
    y
}

#[test]
fn constant_gaunt_cooling_matches_exact_solution_once_per_cell() {
    let (initial, original) = fixture(false, [0.0; 3]);
    let sources = silent_sources(&initial, &original);
    let law = BremsstrahlungGauntFactorV1::Constant(1.2);
    let model = cooling(&initial, law);
    let duration = 2e-4;
    let amounts: Vec<_> = initial.amounts().iter().map(|n| n.value()).collect();
    let te = initial.electron_temperature().value();
    let capacity = CV * amounts[5];
    let coefficient = continuum_power(&amounts, initial.volume().value(), te, law) / te.sqrt();
    let expected = (te.sqrt() - coefficient * duration / (2.0 * capacity)).powi(2);
    let mut solve = controls(1e-8);
    solve.local.initial_step_s = duration / 4.0;
    solve.local.maximum_step_s = duration / 2.0;
    for edges in [&sources[..1], &sources[..]] {
        let result = advance_ionization_network_with_optically_thin_bremsstrahlung(
            &initial,
            edges,
            &[],
            None,
            &model,
            Time::new(duration),
            solve,
            || false,
        )
        .unwrap();
        assert_eq!(result.state().amounts(), initial.amounts());
        assert_eq!(
            result.state().heavy_temperature(),
            initial.heavy_temperature()
        );
        assert_eq!(result.free_free_cooling_source(), Some(&model));
        near(
            result.state().electron_temperature().value(),
            expected,
            5e-13,
            2e-9,
        );
        let b = result.balance();
        near(
            b.free_free_escaped_energy_j,
            capacity * (te - expected),
            5e-12,
            2e-8,
        );
        assert_eq!(b.absorbed_photons, 0.0);
        assert_eq!(b.emitted_photons, 0.0);
        assert_eq!(b.emitted_radiation_energy_j, 0.0);
        assert_eq!(b.net_radiation_energy_j, 0.0);
        assert_eq!(b.net_gas_radiation_energy_j, -b.free_free_escaped_energy_j);
        assert!(b.free_free_evaluations > 1 && b.rate_evaluations >= b.free_free_evaluations);
        assert!(b.energy_residual_j.abs() <= solve.local.energy_tolerance_j);
    }
}

#[test]
fn evolving_gaunt_factor_follows_independent_cooling_ode() {
    let (initial, original) = fixture(false, [0.0; 3]);
    let sources = silent_sources(&initial, &original);
    let law = BremsstrahlungGauntFactorV1::ItohNonrelativistic;
    let model = cooling(&initial, law);
    let duration = 2e-4;
    let amounts: Vec<_> = initial.amounts().iter().map(|n| n.value()).collect();
    let capacity = CV * amounts[5];
    let te = initial.electron_temperature().value();
    let expected = rk4([te], duration, |y| {
        [-continuum_power(&amounts, initial.volume().value(), y[0], law) / capacity]
    })[0];
    let coefficient = continuum_power(&amounts, initial.volume().value(), te, law) / te.sqrt();
    let frozen = (te.sqrt() - coefficient * duration / (2.0 * capacity)).powi(2);
    assert!((expected - frozen).abs() > 1.0);
    let mut solve = controls(1e-9);
    solve.local.initial_step_s = duration / 4.0;
    solve.local.maximum_step_s = duration / 4.0;
    let result = advance_ionization_network_with_optically_thin_bremsstrahlung(
        &initial,
        &sources,
        &[],
        None,
        &model,
        Time::new(duration),
        solve,
        || false,
    )
    .unwrap();
    near(
        result.state().electron_temperature().value(),
        expected,
        0.0,
        5e-3,
    );
    near(
        result.balance().free_free_escaped_energy_j,
        capacity * (te - expected),
        0.0,
        5e-2,
    );
    assert!(result.balance().energy_residual_j.abs() <= solve.local.energy_tolerance_j);
}

#[test]
fn cooling_changes_multistage_kinetics_against_simultaneous_ode() {
    let (initial, sources) = fixture(false, [0.0; 3]);
    let photons = field();
    let law = BremsstrahlungGauntFactorV1::Constant(1.2);
    let model = cooling(&initial, law);
    let duration = 5e-6;
    let mut y = [0.0; 26];
    for (target, n) in y.iter_mut().zip(initial.amounts()) {
        *target = n.value();
    }
    y[6] = CV * y[5] * initial.electron_temperature().value();
    let expected = rk4(y, duration, |y| {
        let base = rhs(
            &core::array::from_fn(|i| y[i]),
            initial.volume().value(),
            &photons,
        );
        let mut rates = [0.0; 26];
        rates[..25].copy_from_slice(&base);
        let power = continuum_power(y, initial.volume().value(), y[6] / (CV * y[5]), law);
        rates[6] -= power;
        rates[25] = power;
        rates
    });
    let solve = controls(2e-8);
    let result = advance_ionization_network_with_optically_thin_bremsstrahlung(
        &initial,
        &sources,
        &photons,
        None,
        &model,
        Time::new(duration),
        solve,
        || false,
    )
    .unwrap();
    for (n, expected) in result.state().amounts().iter().zip(&expected[..6]) {
        near(n.value(), *expected, 0.0, 3e-7);
    }
    near(
        result.state().electron_temperature().value(),
        expected[6] / (CV * expected[5]),
        0.0,
        0.03,
    );
    near(
        result.balance().free_free_escaped_energy_j,
        expected[25],
        2e-6,
        1e-2,
    );
    let uncooled = advance_ionization_network(
        &initial,
        &sources,
        &photons,
        Time::new(duration),
        solve,
        || false,
    )
    .unwrap();
    assert!(
        (result.state().amounts()[5].value() - uncooled.state().amounts()[5].value()).abs() > 1e-7
    );
    let b = result.balance();
    let gas = b.chemical_storage_change_j
        + b.electron_sensible_energy_change_j
        + b.heavy_sensible_energy_change_j;
    near(
        gas,
        b.net_gas_radiation_energy_j,
        0.0,
        solve.local.energy_tolerance_j,
    );
    near(
        b.net_gas_radiation_energy_j,
        b.net_radiation_energy_j - b.free_free_escaped_energy_j,
        0.0,
        solve.local.energy_tolerance_j,
    );
    assert_accounts(&result);
}

#[test]
fn finite_and_diffuse_owners_keep_continuum_escape_out_of_bound_free_photons() {
    let (initial, original) = fixture(false, [0.0; 3]);
    let sources = silent_sources(&initial, &original);
    let model = cooling(&initial, BremsstrahlungGauntFactorV1::Constant(1.2));
    let duration = 2e-4;
    let mut local = controls(1e-8);
    local.local.initial_step_s = duration / 4.0;
    local.local.maximum_step_s = duration / 4.0;
    let control = NetworkRadiationSolve::new(local);
    let cell =
        NetworkRadiationCell::new(&initial, &sources).with_optically_thin_bremsstrahlung(&model);
    let finite = advance_network_radiation_transport(
        &[cell, cell],
        Area::new(100.0),
        &[],
        Time::new(duration),
        control,
        || false,
    )
    .unwrap();
    let diffuse = advance_network_diffuse_recombination(
        &[cell, cell],
        Area::new(100.0),
        &[],
        Time::new(duration),
        NetworkDiffuseSolve::new(control, DiffuseRecombinationSolve::default()),
        || false,
    )
    .unwrap();
    for result in [&finite, &diffuse] {
        let b = result.balance();
        assert!(b.free_free_escaped_energy_j > 0.0);
        assert_eq!(b.absorbed_photons, 0.0);
        assert_eq!(b.emitted_photons, 0.0);
        assert_eq!(b.net_radiation_energy_j, 0.0);
        near(
            b.gas_energy_change_j,
            -b.free_free_escaped_energy_j,
            0.0,
            local.local.energy_tolerance_j,
        );
        assert_eq!(b.net_gas_radiation_energy_j, -b.free_free_escaped_energy_j);
        assert!(b.rate_evaluations >= b.free_free_evaluations && b.free_free_evaluations >= 2);
        for c in result.cells() {
            assert_eq!(c.free_free_cooling_source(), Some(&model));
        }
    }
    for (a, b) in finite.cells().iter().zip(diffuse.cells()) {
        near(
            a.state().electron_temperature().value(),
            b.state().electron_temperature().value(),
            0.0,
            2e-8,
        );
    }
    let radiation = diffuse.diffuse_radiation().unwrap();
    assert_eq!(radiation.balance().escaped_photons, 0.0);
    assert_eq!(radiation.balance().escaped_energy_j, 0.0);
    assert_eq!(
        radiation.balance().free_free_escaped_energy_j,
        diffuse.balance().free_free_escaped_energy_j
    );
    assert!(radiation.balance().energy_residual_j.abs() <= local.local.energy_tolerance_j);
    assert!(
        radiation
            .interfaces()
            .iter()
            .all(|f| f.right_energy_j == 0.0 && f.left_energy_j == 0.0)
    );
}

#[test]
fn source_binding_and_shared_cooling_budget_refuse_atomically() {
    let (initial, original) = fixture(false, [0.0; 3]);
    let sources = silent_sources(&initial, &original);
    let model = cooling(&initial, BremsstrahlungGauntFactorV1::Constant(1.2));
    let before = initial.clone();
    let mut bounded = controls(1e-7);
    bounded.maximum_rate_evaluations = 1;
    assert!(matches!(
        advance_ionization_network_with_optically_thin_bremsstrahlung(
            &initial,
            &sources,
            &[],
            None,
            &model,
            Time::new(1e-5),
            bounded,
            || false,
        ),
        Err(IonizationNetworkError::StepBudget)
    ));
    let cell =
        NetworkRadiationCell::new(&initial, &sources).with_optically_thin_bremsstrahlung(&model);
    assert!(matches!(
        advance_network_radiation_transport(
            &[cell, cell],
            Area::new(100.0),
            &[],
            Time::new(1e-5),
            NetworkRadiationSolve::new(bounded),
            || false,
        ),
        Err(NetworkRadiationError::WorkBudget)
    ));
    let identity = advance_ionization_network_with_optically_thin_bremsstrahlung(
        &initial,
        &sources,
        &[],
        None,
        &model,
        Time::new(0.0),
        bounded,
        || false,
    )
    .unwrap();
    assert_eq!(identity.state(), &initial);
    assert_eq!(identity.balance().rate_evaluations, 1);
    assert_eq!(identity.balance().free_free_evaluations, 1);
    let (other, _) = fixture(false, [2e6, 0.0, 0.0]);
    let wrong = cooling(&other, BremsstrahlungGauntFactorV1::Constant(1.2));
    assert!(
        advance_ionization_network_with_optically_thin_bremsstrahlung(
            &initial,
            &sources,
            &[],
            None,
            &wrong,
            Time::new(0.0),
            bounded,
            || false,
        )
        .is_err()
    );
    assert!(matches!(
        advance_ionization_network_with_optically_thin_bremsstrahlung(
            &initial,
            &sources,
            &[],
            None,
            &model,
            Time::new(1.0),
            controls(1e-7),
            || true,
        ),
        Err(IonizationNetworkError::Cancelled)
    ));
    let mut restricted = model.controls();
    restricted.minimum_electron_temperature_k = 19999.0;
    let restricted =
        ElectronIonBremsstrahlungV1::new(initial.thermodynamics().clone(), restricted).unwrap();
    assert!(
        advance_ionization_network_with_optically_thin_bremsstrahlung(
            &initial,
            &sources,
            &[],
            None,
            &restricted,
            Time::new(2e-4),
            controls(1e-7),
            || false,
        )
        .is_err()
    );
    assert_eq!(initial, before);
}
