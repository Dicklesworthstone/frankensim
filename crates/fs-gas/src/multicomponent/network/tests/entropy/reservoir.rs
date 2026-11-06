//! Independent open-system polynomial-Cp entropy and native pressure-work checks.
use super::*;
use crate::multicomponent::Totals;

// Cp_i/R = a_i + b_i*T, specific R_i = r_i, source pressure p0=1.
const ROWS: [[f64; 3]; 2] = [[3.5, 0.25, 1.0], [5.0, 0.7, 0.5]];
const T_MIN: f64 = 0.001;

fn independent_temperature(q: Cell<'_>) -> f64 {
    let rho: f64 = q.partial_densities.iter().sum();
    let a: f64 = ROWS
        .iter()
        .zip(q.partial_densities)
        .map(|(row, density)| density * row[2] * (row[0] - 1.0))
        .sum();
    let b: f64 = ROWS
        .iter()
        .zip(q.partial_densities)
        .map(|(row, density)| 0.5 * density * row[2] * row[1])
        .sum();
    let target = q.energy - 0.5 * q.momentum * q.momentum / rho + a * T_MIN + b * T_MIN * T_MIN;
    2.0 * target / (a + (a * a + 4.0 * b * target).sqrt())
}

fn independent_entropy(network: &MixtureNetwork1d) -> f64 {
    let density = |q: Cell<'_>| {
        let t = independent_temperature(q);
        ROWS.iter()
            .zip(q.partial_densities)
            .filter(|(_, rho)| **rho > 0.0)
            .map(|([a, b, r], rho)| rho * r * ((a - 1.0) * t.ln() + b * t - (rho * r).ln()))
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

// A stationary prescribed bath receives the opposite conservative transfer.
// Its source-referenced chemical effort is h_i,ref/T - s_i; advective pressure
// work is already inside the native energy transfer. No source entropy methods
// or production entropy variables/reducers are used by this oracle.
fn reservoir_input(source: Cell<'_>, transfer: &Totals) -> f64 {
    assert_eq!(source.momentum, 0.0);
    let t = independent_temperature(source);
    let mut input = transfer.energy / t;
    for (k, [a, b, r]) in ROWS.into_iter().enumerate() {
        if transfer.species_masses[k] == 0.0 {
            continue;
        }
        let rho = source.partial_densities[k];
        assert!(rho > 0.0);
        let u = r * ((a - 1.0) * (t - T_MIN) + 0.5 * b * (t * t - T_MIN * T_MIN));
        let s = r * ((a - 1.0) * t.ln() + b * t - (rho * r).ln());
        input -= (u / t + r - s) * transfer.species_masses[k];
    }
    input
}

fn boundary_input(network: &MixtureNetwork1d, sources: &[[Option<Cell<'_>>; 2]]) -> [f64; 2] {
    let mut input = [0.0; 2];
    for (report, ends) in network
        .accepted_branch_reports()
        .unwrap()
        .iter()
        .zip(sources)
    {
        for (end, transfer) in [report.left_boundary_change, report.right_boundary_change]
            .iter()
            .enumerate()
        {
            if let Some(source) = ends[end] {
                input[end] += reservoir_input(source, transfer);
            }
        }
    }
    input
}

fn driven(gas: &MixtureGas, reversed: bool) -> MixtureNetwork1d {
    let inlet = branch(gas, 4, 0.1, |x| {
        let x = if reversed { 1.0 - x } else { x };
        ([1.2, 0.0], 0.0, if x > 0.75 { 1.6 } else { 1.0 })
    });
    let periodic = MixtureEuler1d::new(
        gas.clone(),
        Length::new(0.25),
        0.02,
        vec![gas.conserved(&[0.7, 0.6], 0.3, 1.0).unwrap(); 4],
        Boundaries::Periodic,
        || false,
    )
    .unwrap();
    MixtureNetwork1d::new(
        vec![
            inlet,
            branch(gas, 4, 0.04, |_| ([0.0, 2.0], 0.0, 1.0)),
            branch(gas, 4, 0.06, |_| ([0.0, 2.0], 0.0, 1.0)),
            periodic,
        ],
        vec![Plenum::new(0.03, gas.conserved(&[0.0, 2.0], 0.0, 1.0).unwrap()).unwrap()],
        (0..3)
            .map(|index| PlenumConnection {
                branch: index,
                end: if (index == 0) != reversed {
                    BranchEnd::Right
                } else {
                    BranchEnd::Left
                },
                plenum: 0,
            })
            .collect(),
        || false,
    )
    .unwrap()
}

fn drives(states: &[Conserved; 3], reversed: bool) -> [[Option<Cell<'_>>; 2]; 4] {
    let mut ends = [[None; 2]; 4];
    ends[0][usize::from(reversed)] = Some(states[0].as_cell());
    ends[1][usize::from(!reversed)] = Some(states[1].as_cell());
    ends[2][usize::from(!reversed)] = Some(states[2].as_cell());
    ends
}

fn reservoirs(gas: &MixtureGas) -> [Conserved; 3] {
    [[0.4, 3.2], [0.2, 1.0], [0.8, 1.2]]
        .map(|rho| gas.conserved(&rho, 0.0, rho[0] + 0.5 * rho[1]).unwrap())
}

#[test]
fn g1_open_moving_network_matches_independent_reservoir_entropy_and_extensive_work() {
    let gas = caloric_gas(true, 0.0);
    let states = reservoirs(&gas);
    for rate in [0.0, -0.2, 0.2] {
        let controls = if rate < 0.0 {
            StepControl::default()
        } else {
            first_order()
        };
        for reversed in [false, true] {
            let sources = drives(&states, reversed);
            let mut network = driven(&gas, reversed);
            let mut native = driven(&gas, reversed);
            let before_s = independent_entropy(&network);
            let before = integrated(&network);
            let (step, entropy) = network
                .step_with_reservoirs_and_volume_rates_and_entropy(
                    Time::new(0.002),
                    controls,
                    &sources,
                    &[rate],
                    MixtureEntropyControl::default(),
                    || false,
                )
                .unwrap();
            let ordinary = native
                .step_with_reservoirs_and_volume_rates(
                    Time::new(0.002),
                    controls,
                    &sources,
                    &[rate],
                    || false,
                )
                .unwrap();
            assert_eq!(step, ordinary);
            assert_eq!(snapshot(&network), snapshot(&native));
            assert_eq!(geometry(&network), geometry(&native));
            let input = boundary_input(&network, &sources);
            for end in 0..2 {
                close(
                    entropy.boundary_entropy_input_j_k[end],
                    input[end],
                    1.0,
                    3e-13,
                );
            }
            let measured = independent_entropy(&network) - before_s;
            close(entropy.source_entropy_change_j_k, measured, 1.0, 4e-13);
            close(
                entropy.source_entropy_generation_j_k,
                measured - input.iter().sum::<f64>(),
                1.0,
                4e-13,
            );
            assert!(entropy.entropy_generation_j_k > 0.0);
            assert_eq!(entropy.rejected_entropy_steps, 0);
            let after = integrated(&network);
            close(
                after.1 - before.1,
                step.boundary_change.absolute_energy(&gas).unwrap() + step.boundary_work_on_gas_j,
                1.0,
                3e-14,
            );
            for k in 0..2 {
                close(
                    after.0[k] - before.0[k],
                    step.boundary_change.species_masses[k],
                    1.0,
                    3e-14,
                );
            }
            assert!(
                network
                    .accepted_port_exchanges()
                    .unwrap()
                    .iter()
                    .any(|port| port.energy_to_plenum != 0.0)
            );
            // These nonzero periodic advective faces are internal, not exterior throughput.
            for i in 0..4 {
                assert_eq!(network.branches()[3].cell(i), native.branches()[3].cell(i));
            }
        }
    }
}

#[test]
fn g1_reservoir_outflow_can_decrease_network_entropy_with_positive_total_generation() {
    let gas = caloric_gas(true, 0.0);
    let mut network = MixtureNetwork1d::new(
        vec![branch(&gas, 4, 0.1, |_| ([0.7, 0.6], 0.0, 1.0))],
        vec![Plenum::new(0.1, gas.conserved(&[0.7, 0.6], 0.0, 1.0).unwrap()).unwrap()],
        vec![PlenumConnection {
            branch: 0,
            end: BranchEnd::Right,
            plenum: 0,
        }],
        || false,
    )
    .unwrap();
    let sink = gas.conserved(&[0.35, 0.3], 0.0, 0.5).unwrap();
    let sources = [[Some(sink.as_cell()), None]];
    let before = independent_entropy(&network);
    let (step, entropy) = network
        .step_with_reservoirs_and_entropy(
            Time::new(0.002),
            first_order(),
            &sources,
            MixtureEntropyControl::default(),
            || false,
        )
        .unwrap();
    assert!(step.boundary_change.mass() < 0.0);
    assert!(entropy.source_entropy_change_j_k < 0.0);
    assert!(entropy.boundary_entropy_input_j_k[0] < entropy.source_entropy_change_j_k);
    assert!(entropy.entropy_generation_j_k > 0.0);
    close(
        entropy.source_entropy_change_j_k,
        independent_entropy(&network) - before,
        1.0,
        3e-13,
    );
    close(
        entropy.boundary_entropy_input_j_k[0],
        boundary_input(&network, &sources)[0],
        1.0,
        3e-13,
    );
    assert_eq!(entropy.boundary_entropy_input_j_k[1], 0.0);
}

fn datum_gas(formation: f64, datum: f64) -> MixtureGas {
    let original = caloric_gas(true, formation);
    let models = original
        .source()
        .models()
        .iter()
        .enumerate()
        .map(|(k, model)| {
            let mut card = model.regions()[0].card().clone();
            card.parameters.get_mut("a8").unwrap().value = if k == 1 { datum } else { 0.0 };
            Nasa9StandardStateModelV1::new(
                model.species().clone(),
                model.molar_mass(),
                model.convention().clone(),
                vec![Nasa9RegionV1::from_card(card).unwrap()],
            )
            .unwrap()
        })
        .collect();
    MixtureGas::new(
        crate::multicomponent::tests::bank(models),
        original.domain(),
        || false,
    )
    .unwrap()
}

fn introducing(gas: &MixtureGas, reversed: bool) -> MixtureNetwork1d {
    MixtureNetwork1d::new(
        vec![branch(gas, 4, 0.1, |_| ([1.0, 0.0], 0.0, 1.0))],
        vec![Plenum::new(0.1, gas.conserved(&[1.0, 0.0], 0.0, 1.0).unwrap()).unwrap()],
        vec![PlenumConnection {
            branch: 0,
            end: if reversed {
                BranchEnd::Left
            } else {
                BranchEnd::Right
            },
            plenum: 0,
        }],
        || false,
    )
    .unwrap()
}

#[test]
fn g3_incoming_absent_species_retains_physical_entropy_datums_without_changing_admission() {
    let baseline_gas = caloric_gas(true, 0.0);
    let datum = 2.0_f64.powi(50);
    let shifted_gas = datum_gas(1e18, datum);
    for reversed in [false, true] {
        let source = baseline_gas.conserved(&[0.0, 4.0], 0.0, 2.0).unwrap();
        let shifted_source = shifted_gas.conserved(&[0.0, 4.0], 0.0, 2.0).unwrap();
        let mut sources = [[None; 2]];
        let mut shifted_sources = [[None; 2]];
        let end = usize::from(reversed);
        sources[0][end] = Some(source.as_cell());
        shifted_sources[0][end] = Some(shifted_source.as_cell());
        let mut baseline = introducing(&baseline_gas, reversed);
        let mut changed = introducing(&shifted_gas, reversed);
        assert_eq!(integrated(&baseline).0[1], 0.0);
        let (step, original) = baseline
            .step_with_reservoirs_and_entropy(
                Time::new(0.002),
                first_order(),
                &sources,
                MixtureEntropyControl::default(),
                || false,
            )
            .unwrap();
        let (shifted_step, shifted) = changed
            .step_with_reservoirs_and_entropy(
                Time::new(0.002),
                first_order(),
                &shifted_sources,
                MixtureEntropyControl::default(),
                || false,
            )
            .unwrap();
        assert_eq!(step, shifted_step);
        assert_eq!(snapshot(&baseline), snapshot(&changed));
        assert!(step.boundary_change.species_masses[1] > 0.0);
        assert!(shifted.maximum_relative_species_residual < 1e-12);
        assert_eq!(
            original.entropy_generation_j_k,
            shifted.entropy_generation_j_k
        );
        let physical_change = 0.5 * datum * integrated(&changed).0[1];
        let physical_input = 0.5 * datum * shifted_step.boundary_change.species_masses[1];
        assert!(physical_input.abs() > 1e7);
        close(
            shifted.source_entropy_change_j_k - original.source_entropy_change_j_k,
            physical_change,
            physical_change.abs(),
            3e-13,
        );
        close(
            shifted.boundary_entropy_input_j_k[end] - original.boundary_entropy_input_j_k[end],
            physical_input,
            physical_input.abs(),
            3e-13,
        );
        close(
            shifted.source_entropy_generation_j_k - original.source_entropy_generation_j_k,
            shifted.reference_entropy_residual_j_k - original.reference_entropy_residual_j_k,
            1.0,
            3e-13,
        );
    }
}

fn expanding(gas: &MixtureGas) -> MixtureNetwork1d {
    let branch = MixtureEuler1d::new(
        gas.clone(),
        Length::new(10.0),
        1e-6,
        vec![gas.conserved(&[0.7, 0.6], 0.0, 1.0).unwrap(); 2],
        walls(),
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

#[test]
fn g4_open_entropy_retries_discard_exterior_transfers_and_cancel_atomically() {
    let gas = caloric_gas(true, 0.0);
    let source = gas.conserved(&[0.707, 0.606], 0.0, 1.01).unwrap();
    let sources = [[Some(source.as_cell()), None]];
    let allowance = MixtureEntropyControl {
        absolute_tolerance_j_k: 1e-8,
        ..MixtureEntropyControl::default()
    };
    let prime = |network: &mut MixtureNetwork1d| {
        network
            .step_with_reservoirs_and_volume_rates_and_entropy(
                Time::new(0.005),
                first_order(),
                &sources,
                &[-2.0],
                MixtureEntropyControl::default(),
                || false,
            )
            .unwrap()
    };
    let mut expected = expanding(&gas);
    prime(&mut expected);
    let before_s = independent_entropy(&expected);
    let volume = expected.plenums()[0].volume_m3();
    let mut polls = 0;
    let (step, entropy) = expected
        .step_with_reservoirs_and_volume_rates_and_entropy(
            Time::new(0.05),
            first_order(),
            &sources,
            &[2.0],
            allowance,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert!(step.halvings >= 3);
    assert_eq!(step.halvings, entropy.rejected_entropy_steps);
    assert!(step.elapsed.value() > 0.0 && step.elapsed.value() < 0.05);
    assert!(step.boundary_change.mass() > 0.0);
    close(
        expected.plenums()[0].volume_m3(),
        volume + 2.0 * step.elapsed.value(),
        1.0,
        2e-15,
    );
    let input = boundary_input(&expected, &sources);
    close(entropy.boundary_entropy_input_j_k[0], input[0], 1.0, 1e-14);
    close(
        entropy.source_entropy_generation_j_k,
        independent_entropy(&expected) - before_s - input.iter().sum::<f64>(),
        1.0,
        5e-13,
    );
    assert!(entropy.entropy_generation_j_k >= -allowance.absolute_tolerance_j_k);

    for cancel in [false, true] {
        let mut network = expanding(&gas);
        prime(&mut network);
        let before = snapshot(&network);
        let volumes = geometry(&network);
        let work = network.accepted_plenum_boundary_work().unwrap().to_vec();
        let mut visits = 0;
        let result = network.step_with_reservoirs_and_volume_rates_and_entropy(
            Time::new(0.05),
            StepControl {
                maximum_halvings: if cancel { 16 } else { 0 },
                ..first_order()
            },
            &sources,
            &[2.0],
            allowance,
            || {
                visits += 1;
                cancel && visits == polls
            },
        );
        if cancel {
            assert_eq!(result, Err(Error::Cancelled));
        } else {
            assert!(
                matches!(result, Err(Error::StepRetries { attempts: 1, last })
            if matches!(*last, Error::EntropyDecrease { .. }))
            );
        }
        assert_eq!(snapshot(&network), before);
        assert_eq!(geometry(&network), volumes);
        assert_eq!(network.accepted_plenum_boundary_work().unwrap(), work);
        let initial_s = independent_entropy(&network);
        let (recovered, entropy) = network
            .step_with_reservoirs_and_volume_rates_and_entropy(
                Time::new(0.05),
                first_order(),
                &sources,
                &[2.0],
                allowance,
                || false,
            )
            .unwrap();
        assert_eq!(recovered.elapsed, step.elapsed);
        let input = boundary_input(&network, &sources);
        close(
            entropy.source_entropy_generation_j_k,
            independent_entropy(&network) - initial_s - input.iter().sum::<f64>(),
            1.0,
            5e-13,
        );
    }
}

#[test]
fn g0_open_entropy_preserves_closed_zero_motion_and_refuses_unowned_or_zero_activity_flux() {
    let gas = caloric_gas(true, 0.0);
    let mut closed = driven(&gas, false);
    let mut absent = driven(&gas, false);
    let (native, entropy) = closed
        .step_with_volume_rates_and_entropy(
            Time::new(0.002),
            first_order(),
            &[-0.2],
            MixtureEntropyControl::default(),
            || false,
        )
        .unwrap();
    let (same, open) = absent
        .step_with_reservoirs_and_volume_rates_and_entropy(
            Time::new(0.002),
            first_order(),
            &[[None; 2]; 4],
            &[-0.2],
            MixtureEntropyControl::default(),
            || false,
        )
        .unwrap();
    assert_eq!(native, same);
    assert_eq!(snapshot(&closed), snapshot(&absent));
    assert_eq!(open.entropy_generation_j_k, entropy.entropy_change_j_k);
    assert_eq!(
        open.source_entropy_change_j_k,
        entropy.source_entropy_change_j_k
    );
    assert_eq!(open.boundary_entropy_input_j_k, [0.0; 2]);
    let states = reservoirs(&gas);
    let sources = drives(&states, false);
    let mut rigid = driven(&gas, false);
    let mut zero = driven(&gas, false);
    assert_eq!(
        rigid
            .step_with_reservoirs_and_entropy(
                Time::new(0.002),
                first_order(),
                &sources,
                MixtureEntropyControl::default(),
                || false,
            )
            .unwrap(),
        zero.step_with_reservoirs_and_volume_rates_and_entropy(
            Time::new(0.002),
            first_order(),
            &sources,
            &[0.0],
            MixtureEntropyControl::default(),
            || false,
        )
        .unwrap()
    );
    assert_eq!(snapshot(&rigid), snapshot(&zero));

    // Pure B cannot receive a nonzero A transfer at finite chemical activity.
    let mut refused = introducing(&gas, false);
    let sink = gas.conserved(&[0.0, 1.0], 0.0, 0.5).unwrap();
    let before = snapshot(&refused);
    assert!(matches!(refused.step_with_reservoirs_and_entropy(
        Time::new(0.002), StepControl { maximum_halvings: 0, ..first_order() },
        &[[Some(sink.as_cell()), None]], MixtureEntropyControl::default(), || false,
    ), Err(Error::StepRetries { last, .. }) if matches!(*last, Error::EntropyDomain(_))));
    assert_eq!(snapshot(&refused), before);

    let branch = MixtureEuler1d::new(
        gas.clone(),
        Length::new(0.5),
        0.1,
        vec![gas.conserved(&[0.7, 0.6], 0.0, 1.0).unwrap(); 2],
        Boundaries::Ends {
            left: Boundary::Transmissive,
            right: Boundary::Transmissive,
        },
        || false,
    )
    .unwrap();
    let mut unowned = MixtureNetwork1d::new(vec![branch], vec![], vec![], || false).unwrap();
    let before = snapshot(&unowned);
    assert!(matches!(
        unowned.step_with_reservoirs_and_entropy(
            Time::new(0.002),
            first_order(),
            &[[Some(states[0].as_cell()), None]],
            MixtureEntropyControl::default(),
            || false,
        ),
        Err(Error::EntropyDomain(_))
    ));
    assert_eq!(snapshot(&unowned), before);
}
