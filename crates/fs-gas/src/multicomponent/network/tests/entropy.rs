//! Independent polynomial-caloric entropy, work and rollback oracles.
use super::*;
use crate::multicomponent::{MixtureEntropyControl, Reconstruction};
use fs_thermochem::{Nasa9RegionV1, Nasa9StandardStateModelV1};

mod reservoir;

fn first_order() -> StepControl {
    StepControl {
        reconstruction: Reconstruction::FirstOrder,
        ..StepControl::default()
    }
}

// Original NASA source, evaluated independently of its entropy methods and of
// the numerical account: Cp_i/R = a_i + b_i*T, p0=1, R/M_i = [1, 1/2].
fn source_entropy(network: &MixtureNetwork1d) -> f64 {
    let density = |q: Cell<'_>| {
        let t = network.gas().primitive(q, || false).unwrap().temperature;
        [1.0, 0.5]
            .into_iter()
            .zip([3.5, 5.0])
            .zip([0.25, 0.7])
            .zip(q.partial_densities)
            .filter(|(_, rho)| **rho > 0.0)
            .map(|(((r, a), b), &rho)| rho * r * ((a - 1.0) * t.ln() + b * t - (rho * r).ln()))
            .sum::<f64>()
    };
    network
        .branches()
        .iter()
        .map(|branch| {
            branch
                .mesh()
                .volumes()
                .iter()
                .enumerate()
                .map(|(i, volume)| volume * density(branch.cell(i).unwrap()))
                .sum::<f64>()
        })
        .sum::<f64>()
        + network
            .plenums()
            .iter()
            .map(|plenum| plenum.volume_m3() * density(plenum.state()))
            .sum::<f64>()
}

fn isolated(gas: &MixtureGas) -> MixtureNetwork1d {
    let branch = MixtureEuler1d::new(
        gas.clone(),
        Length::new(10.0),
        0.01,
        vec![gas.conserved(&[0.7, 0.6], 0.0, 1.0).unwrap(); 2],
        Boundaries::Periodic,
        || false,
    )
    .unwrap();
    MixtureNetwork1d::new(
        vec![branch],
        vec![Plenum::new(1.0, gas.conserved(&[0.7, 0.6], 0.0, 1.0).unwrap()).unwrap()],
        vec![],
        || false,
    )
    .unwrap()
}

fn geometry(network: &MixtureNetwork1d) -> Vec<f64> {
    network.plenums().iter().map(Plenum::volume_m3).collect()
}

#[test]
fn g1_closed_network_entropy_matches_independent_variable_cp_inventory_and_native_work() {
    let gas = caloric_gas(true, 0.0);
    for rate in [0.0, -0.2, 0.2] {
        let mut network = connected(&gas, 4, false);
        let mut native = connected(&gas, 4, false);
        let before_entropy = source_entropy(&network);
        let before = integrated(&network);
        let (step, entropy) = network
            .step_with_volume_rates_and_entropy(
                Time::new(0.002),
                first_order(),
                &[rate],
                MixtureEntropyControl::default(),
                || false,
            )
            .unwrap();
        let original = native
            .step_with_volume_rates(Time::new(0.002), first_order(), &[rate], || false)
            .unwrap();
        assert_eq!(
            step, original,
            "entropy admission changed the native accepted flow"
        );
        assert_eq!(snapshot(&network), snapshot(&native));
        assert_eq!(geometry(&network), geometry(&native));
        assert!(
            entropy.entropy_change_j_k > 0.0,
            "finite mixing should generate entropy: {entropy:?}"
        );
        close(
            entropy.source_entropy_change_j_k,
            source_entropy(&network) - before_entropy,
            1.0,
            3e-13,
        );
        close(
            entropy.source_entropy_change_j_k - entropy.reference_entropy_residual_j_k,
            entropy.entropy_change_j_k,
            1.0,
            2e-16,
        );
        assert_eq!(entropy.rejected_entropy_steps, 0);
        assert!(entropy.maximum_relative_species_residual < 2e-14);
        let after = integrated(&network);
        for k in 0..2 {
            close(after.0[k], before.0[k], before.0[k], 2e-14);
        }
        close(after.1 - before.1, step.boundary_work_on_gas_j, 1.0, 3e-14);
        assert!(
            network
                .accepted_port_exchanges()
                .unwrap()
                .iter()
                .any(|port| port.energy_to_plenum != 0.0)
        );
    }
}

#[test]
fn g1_moving_mixture_entropy_and_pressure_work_converge_to_the_reversible_isentrope() {
    let gas = caloric_gas(true, 0.0);
    for rate in [-2.0, 2.0] {
        let mut errors = Vec::new();
        for steps in [20, 40] {
            let mut network = isolated(&gas);
            let before = integrated(&network);
            let initial_entropy = source_entropy(&network);
            let mut entropy_sum = 0.0;
            let mut work = 0.0;
            for _ in 0..steps {
                let (step, entropy) = network
                    .step_with_volume_rates_and_entropy(
                        Time::new(0.25 / f64::from(steps)),
                        first_order(),
                        &[rate],
                        MixtureEntropyControl {
                            absolute_tolerance_j_k: 1e-4,
                            ..MixtureEntropyControl::default()
                        },
                        || false,
                    )
                    .unwrap();
                assert_eq!(step.halvings, 0);
                assert_eq!(step.boundary_change.energy, 0.0);
                entropy_sum += entropy.source_entropy_change_j_k;
                work += step.boundary_work_on_gas_j;
            }
            let chamber = &network.plenums()[0];
            let volume = chamber.volume_m3();
            let temperature = gas
                .primitive(chamber.state(), || false)
                .unwrap()
                .temperature;
            // Fixed masses m_A=0.7, m_B=0.6 give total Cv(T)=2.95+0.385*T
            // and total R=1. Integrate Cv/T*dT + R*dV/V = 0 independently.
            let (mut lo, mut hi) = (0.5_f64, 2.0_f64);
            for _ in 0..80 {
                let t = lo.midpoint(hi);
                if 2.95 * t.ln() + 0.385 * (t - 1.0) + volume.ln() > 0.0 {
                    hi = t;
                } else {
                    lo = t;
                }
            }
            errors.push((temperature - lo.midpoint(hi)).abs());
            let independent_change =
                2.95 * temperature.ln() + 0.385 * (temperature - 1.0) + volume.ln();
            close(entropy_sum, independent_change, 1.0, 5e-13);
            close(
                entropy_sum,
                source_entropy(&network) - initial_entropy,
                1.0,
                5e-13,
            );
            let after = integrated(&network);
            for k in 0..2 {
                close(after.0[k], before.0[k], before.0[k], 3e-14);
            }
            close(after.1 - before.1, work, work.abs(), 2e-13);
            close(volume, 1.0 + rate * 0.25, 1.0, 2e-14);
            assert_eq!(work.is_sign_positive(), rate.is_sign_negative());
            assert_eq!(entropy_sum.is_sign_positive(), rate.is_sign_negative());
        }
        let order = (errors[0] / errors[1]).log2();
        assert!(
            (1.8..2.2).contains(&order),
            "order={order}, errors={errors:?}"
        );
    }
}

#[test]
fn g3_network_entropy_admission_is_independent_of_species_entropy_and_formation_datums() {
    let gas = caloric_gas(true, 0.0);
    let shifted = caloric_gas(true, 1e18);
    let models = shifted
        .source()
        .models()
        .iter()
        .enumerate()
        .map(|(k, model)| {
            let mut card = model.regions()[0].card().clone();
            card.parameters.get_mut("a8").unwrap().value = if k == 0 {
                2.0_f64.powi(50)
            } else {
                -2.0_f64.powi(51)
            };
            Nasa9StandardStateModelV1::new(
                model.species().clone(),
                model.molar_mass(),
                model.convention().clone(),
                vec![Nasa9RegionV1::from_card(card).unwrap()],
            )
            .unwrap()
        })
        .collect();
    let shifted = MixtureGas::new(
        crate::multicomponent::tests::bank(models),
        gas.domain(),
        || false,
    )
    .unwrap();
    let mut baseline = connected(&gas, 4, false);
    let mut changed = connected(&shifted, 4, false);
    let (original_step, original_entropy) = baseline
        .step_with_volume_rates_and_entropy(
            Time::new(0.002),
            first_order(),
            &[-0.2],
            MixtureEntropyControl::default(),
            || false,
        )
        .unwrap();
    let (changed_step, changed_entropy) = changed
        .step_with_volume_rates_and_entropy(
            Time::new(0.002),
            first_order(),
            &[-0.2],
            MixtureEntropyControl::default(),
            || false,
        )
        .unwrap();
    assert_eq!(original_step, changed_step);
    assert_eq!(snapshot(&baseline), snapshot(&changed));
    assert_eq!(geometry(&baseline), geometry(&changed));
    assert_eq!(
        original_entropy.entropy_change_j_k,
        changed_entropy.entropy_change_j_k
    );
    assert_eq!(
        original_entropy.maximum_relative_species_residual,
        changed_entropy.maximum_relative_species_residual
    );
    assert!(
        changed_entropy.reference_entropy_residual_j_k.abs() > 1e-6,
        "fixture must expose the numerical source-reference residue: {changed_entropy:?}"
    );
    close(
        changed_entropy.source_entropy_change_j_k - changed_entropy.reference_entropy_residual_j_k,
        original_entropy.entropy_change_j_k,
        1.0,
        2e-15,
    );
}

#[test]
fn g4_entropy_refusal_retries_the_whole_moving_network_and_exhaustion_publishes_nothing() {
    let gas = caloric_gas(true, 0.0);
    let entropy_controls = MixtureEntropyControl {
        absolute_tolerance_j_k: 1e-8,
        ..MixtureEntropyControl::default()
    };
    let mut network = isolated(&gas);
    let (step, entropy) = network
        .step_with_volume_rates_and_entropy(
            Time::new(0.05),
            first_order(),
            &[2.0],
            entropy_controls,
            || false,
        )
        .unwrap();
    assert!(step.elapsed.value() > 0.0 && step.elapsed.value() < 0.05);
    assert_eq!(entropy.rejected_entropy_steps, step.halvings);
    assert!(entropy.rejected_entropy_steps >= 3);
    assert!(entropy.entropy_change_j_k >= -entropy_controls.absolute_tolerance_j_k);
    close(
        network.plenums()[0].volume_m3(),
        1.0 + 2.0 * step.elapsed.value(),
        1.0,
        1e-15,
    );
    let mut direct = isolated(&gas);
    direct
        .step_with_volume_rates(step.elapsed, first_order(), &[2.0], || false)
        .unwrap();
    assert_eq!(network.plenums()[0].state(), direct.plenums()[0].state());
    assert_eq!(geometry(&network), geometry(&direct));
    close(
        entropy.source_entropy_change_j_k,
        source_entropy(&network) - source_entropy(&isolated(&gas)),
        1.0,
        3e-13,
    );

    let mut failed = isolated(&gas);
    let before = snapshot(&failed);
    let volumes = geometry(&failed);
    let result = failed.step_with_volume_rates_and_entropy(
        Time::new(0.05),
        StepControl {
            maximum_halvings: 0,
            ..first_order()
        },
        &[2.0],
        entropy_controls,
        || false,
    );
    assert!(
        matches!(result, Err(Error::StepRetries { attempts: 1, last })
        if matches!(*last, Error::EntropyDecrease { .. }))
    );
    assert_eq!(snapshot(&failed), before);
    assert_eq!(geometry(&failed), volumes);
    assert!(failed.accepted_plenum_boundary_work().is_none());
    failed
        .step_with_volume_rates_and_entropy(
            Time::new(0.05),
            first_order(),
            &[2.0],
            entropy_controls,
            || false,
        )
        .unwrap();
    assert_eq!(snapshot(&failed), snapshot(&network));
    assert_eq!(geometry(&failed), geometry(&network));
}

#[test]
fn g4_late_entropy_cancellation_retains_moving_geometry_work_reports_and_restart() {
    let gas = caloric_gas(true, 0.0);
    let advance = |network: &mut MixtureNetwork1d| {
        network
            .step_with_volume_rates_and_entropy(
                Time::new(0.01),
                first_order(),
                &[-2.0],
                MixtureEntropyControl::default(),
                || false,
            )
            .unwrap()
    };
    let mut probe = isolated(&gas);
    advance(&mut probe);
    let mut polls = 0;
    probe
        .step_with_volume_rates_and_entropy(
            Time::new(0.01),
            first_order(),
            &[-2.0],
            MixtureEntropyControl::default(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut network = isolated(&gas);
    advance(&mut network);
    let before = snapshot(&network);
    let volumes = geometry(&network);
    let work = network.accepted_plenum_boundary_work().unwrap().to_vec();
    let mut observed = 0;
    assert!(matches!(
        network.step_with_volume_rates_and_entropy(
            Time::new(0.01),
            first_order(),
            &[-2.0],
            MixtureEntropyControl::default(),
            || {
                observed += 1;
                observed == polls
            },
        ),
        Err(Error::Cancelled)
    ));
    assert_eq!(snapshot(&network), before);
    assert_eq!(geometry(&network), volumes);
    assert_eq!(network.accepted_plenum_boundary_work().unwrap(), work);
    advance(&mut network);
    assert_eq!(snapshot(&network), snapshot(&probe));
    assert_eq!(geometry(&network), geometry(&probe));
}

#[test]
fn g0_closed_network_entropy_requires_closed_exterior_and_zero_rates_preserve_rigid_identity() {
    let gas = caloric_gas(true, 0.0);
    let mut rigid = connected(&gas, 4, false);
    let mut zero = connected(&gas, 4, false);
    assert_eq!(
        rigid
            .step_with_entropy(
                Time::new(0.002),
                first_order(),
                MixtureEntropyControl::default(),
                || false
            )
            .unwrap(),
        zero.step_with_volume_rates_and_entropy(
            Time::new(0.002),
            first_order(),
            &[0.0],
            MixtureEntropyControl::default(),
            || false
        )
        .unwrap(),
    );
    assert_eq!(snapshot(&rigid), snapshot(&zero));
    let branch = MixtureEuler1d::new(
        gas.clone(),
        Length::new(1.0),
        1.0,
        vec![gas.conserved(&[0.7, 0.6], 0.0, 1.0).unwrap(); 2],
        Boundaries::Ends {
            left: Boundary::Transmissive,
            right: Boundary::Reflecting,
        },
        || false,
    )
    .unwrap();
    let mut open = MixtureNetwork1d::new(vec![branch], vec![], vec![], || false).unwrap();
    let before = snapshot(&open);
    assert!(matches!(
        open.step_with_entropy(
            Time::new(0.002),
            first_order(),
            MixtureEntropyControl::default(),
            || false,
        ),
        Err(Error::EntropyDomain(_))
    ));
    assert_eq!(snapshot(&open), before);
}
