//! Physical storage/flux checks on the original NASA gas and hydrostatic owners.

use super::*;
use crate::dynamics::hydrostatic::{
    AtmosphereQuadrature, HydrostaticReference, IsentropicAtmosphere, IsothermalAtmosphere,
};
use crate::dynamics::{Boundaries, Boundary, FrozenGas, Mesh1d, Reconstruction, StepControl};
use fs_qty::{Acceleration, Length, Pressure, Temperature, Time};

const WALLS: Boundaries = Boundaries::Ends {
    left: Boundary::Reflecting,
    right: Boundary::Reflecting,
};

fn column(
    gas: &FrozenGas,
    isentropic: bool,
    n: usize,
    acceleration: f64,
    datum_position: f64,
    reference: (f64, f64),
) -> Euler1d {
    let positions: Vec<_> = (0..=n).map(|i| (i as f64 / n as f64).powi(2)).collect();
    let volumes = positions.windows(2).map(|x| 0.7 * (x[1] - x[0])).collect();
    let mesh = Mesh1d::new(positions, vec![0.7; n + 1], volumes).unwrap();
    if isentropic {
        Euler1d::new_isentropic_atmosphere(
            gas.clone(),
            mesh,
            IsentropicAtmosphere {
                acceleration: Acceleration::new(acceleration),
                reference_temperature: Temperature::new(reference.0),
                reference_pressure: Pressure::new(reference.1),
                reference_position: Length::new(datum_position),
            },
            AtmosphereQuadrature::default(),
            WALLS,
            || false,
        )
        .unwrap()
    } else {
        Euler1d::new_isothermal_atmosphere(
            gas.clone(),
            mesh,
            IsothermalAtmosphere {
                acceleration: Acceleration::new(acceleration),
                temperature: Temperature::new(reference.0),
                reference_pressure: Pressure::new(reference.1),
                reference_position: Length::new(datum_position),
            },
            WALLS,
            || false,
        )
        .unwrap()
    }
}

fn connect(branch: usize, end: BranchEnd, plenum: usize) -> PlenumConnection {
    PlenumConnection {
        branch,
        end,
        plenum,
    }
}

fn set_velocity(branch: &mut Euler1d, velocity: f64) {
    for state in &mut branch.states {
        state.momentum = state.density * velocity;
        state.energy += 0.5 * state.density * velocity * velocity;
    }
}

fn moving_chain(datum: f64, formation: f64, isentropic: bool) -> GasNetwork1d {
    let gas = crate::dynamics::tests::gas(0.4, formation);
    let mut lower = column(&gas, isentropic, 9, -0.5, 0.0, (2.0, 3.0));
    let base = lower.hydrostatic_reference_face(0).unwrap();
    let middle = lower.hydrostatic_reference_face(9).unwrap();
    let mut upper = column(
        &gas,
        isentropic,
        11,
        -0.5,
        0.0,
        (middle.1.temperature, middle.1.pressure),
    );
    let top = upper.hydrostatic_reference_face(11).unwrap();
    set_velocity(&mut lower, 0.06);
    set_velocity(&mut upper, 0.06);
    GasNetwork1d::new_with_gravity(
        vec![lower, upper],
        vec![
            Plenum::new(0.2, base.0).unwrap(),
            Plenum::new(0.3, middle.0).unwrap(),
            Plenum::new(0.2, top.0).unwrap(),
        ],
        vec![
            connect(0, BranchEnd::Left, 0),
            connect(0, BranchEnd::Right, 1),
            connect(1, BranchEnd::Left, 1),
            connect(1, BranchEnd::Right, 2),
        ],
        NetworkGravity {
            branch_offsets_j_kg: vec![datum, datum + 0.5],
            plenum_potentials_j_kg: vec![datum, datum + 0.5, datum + 1.0],
        },
        || false,
    )
    .unwrap()
}

// Independent direct quadrature of actual conserved fields. Do not use the
// network's total/potential helpers or its reported work to reconstruct storage.
fn storage(network: &GasNetwork1d) -> [f64; 3] {
    let gravity = network.gravity().unwrap();
    let mut sum = [0.0; 3];
    for (b, branch) in network.branches().iter().enumerate() {
        let (g, datum) = match branch.hydrostatic_reference() {
            Some(HydrostaticReference::Isothermal(r)) => {
                (r.acceleration.value(), r.reference_position.value())
            }
            Some(HydrostaticReference::Isentropic(r)) => {
                (r.acceleration.value(), r.reference_position.value())
            }
            None => (0.0, 0.0),
        };
        for (i, q) in branch.states().iter().enumerate() {
            let volume = branch.mesh().volumes()[i];
            let x = 0.5 * (branch.mesh().positions()[i] + branch.mesh().positions()[i + 1]);
            let potential = gravity.branch_offsets_j_kg[b] - g * (x - datum);
            sum[0] += volume * q.density;
            sum[1] += volume * q.energy;
            sum[2] += volume * q.density * potential;
        }
    }
    for (p, plenum) in network.plenums().iter().enumerate() {
        let q = plenum.state();
        sum[0] += plenum.volume_m3() * q.density;
        sum[1] += plenum.volume_m3() * q.energy;
        sum[2] += plenum.volume_m3() * q.density * gravity.plenum_potentials_j_kg[p];
    }
    sum
}

#[test]
fn g0_joined_isothermal_and_isentropic_columns_remain_exactly_at_rest() {
    for isentropic in [false, true] {
        for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
            let gas = crate::dynamics::tests::gas(0.4, 1e18);
            // Both right ends meet at one height; the second local axis is
            // physically reversed, so its signed acceleration is reversed too.
            let first = column(&gas, isentropic, 9, -0.5, 1.0, (2.0, 3.0));
            let second = column(&gas, isentropic, 13, 0.5, 1.0, (2.0, 3.0));
            let at_node = first.hydrostatic_reference_face(9).unwrap().0;
            let initial = vec![first.states().to_vec(), second.states().to_vec()];
            let mut network = GasNetwork1d::new_with_gravity(
                vec![first, second],
                vec![Plenum::new(0.3, at_node).unwrap()],
                vec![
                    connect(0, BranchEnd::Right, 0),
                    connect(1, BranchEnd::Right, 0),
                ],
                NetworkGravity {
                    branch_offsets_j_kg: vec![2.0, 2.0],
                    plenum_potentials_j_kg: vec![2.0],
                },
                || false,
            )
            .unwrap();
            for _ in 0..30 {
                let report = network
                    .step(
                        Time::new(0.003),
                        StepControl {
                            reconstruction,
                            ..StepControl::default()
                        },
                        || false,
                    )
                    .unwrap();
                for (branch, before) in network.branches().iter().zip(&initial) {
                    assert_eq!(branch.states(), before);
                }
                assert_eq!(network.plenums()[0].state(), at_node);
                for exchange in network.accepted_port_exchanges().unwrap() {
                    assert_eq!(exchange.mass_to_plenum, 0.0);
                    assert_eq!(exchange.energy_to_plenum, 0.0);
                }
                let gravity = report.gravity.unwrap();
                assert_eq!(gravity.gravity_work, 0.0);
                assert_eq!(
                    gravity.before_centered_potential_energy,
                    gravity.after_centered_potential_energy
                );
                assert_eq!(gravity.conservation_residual, 0.0);
            }
        }
    }
}

#[test]
fn g1_multilevel_transient_closes_independent_gas_and_potential_storage() {
    for isentropic in [false, true] {
        let mut network = moving_chain(0.0, 0.0, isentropic);
        let initial = storage(&network);
        let mut work = 0.0;
        let mut transported = 0.0;
        for _ in 0..60 {
            let report = network
                .step(Time::new(0.001), StepControl::default(), || false)
                .unwrap();
            assert!(report.conservation_residual.mass.abs() < 2e-13);
            assert!(report.conservation_residual.energy.abs() < 2e-12);
            assert_eq!(report.boundary_change.mass, 0.0);
            assert_eq!(report.boundary_change.energy, 0.0);
            let gravity = report.gravity.unwrap();
            assert!(gravity.conservation_residual.abs() < 3e-12);
            work += gravity.gravity_work;
            transported += network
                .accepted_port_exchanges()
                .unwrap()
                .iter()
                .map(|p| p.mass_to_plenum.abs())
                .sum::<f64>();
            for branch in network.accepted_branch_reports().unwrap() {
                assert!(branch.conservation_residual.energy.abs() < 2e-12);
                assert!(branch.conservation_residual.momentum.abs() < 2e-13);
            }
        }
        let final_state = storage(&network);
        assert!(transported > 1e-4);
        assert!(work < -1e-5);
        assert!((final_state[0] - initial[0]).abs() < 2e-12);
        assert!((final_state[1] - initial[1] - work).abs() < 2e-11);
        assert!((final_state[2] - initial[2] + work).abs() < 2e-11);
        assert!((final_state[1] + final_state[2] - initial[1] - initial[2]).abs() < 2e-11);
    }
}

#[test]
fn g1_external_transport_includes_the_actual_exit_potential_once() {
    let gas = crate::dynamics::tests::gas(0.4, 0.0);
    let mut branch = column(&gas, true, 11, -0.5, 0.0, (2.0, 3.0));
    let plenum = Plenum::new(0.3, branch.hydrostatic_reference_face(0).unwrap().0).unwrap();
    set_velocity(&mut branch, 0.06);
    branch
        .set_boundaries(
            Boundaries::Ends {
                left: Boundary::Reflecting,
                right: Boundary::Transmissive,
            },
            || false,
        )
        .unwrap();
    let mut network = GasNetwork1d::new_with_gravity(
        vec![branch],
        vec![plenum],
        vec![connect(0, BranchEnd::Left, 0)],
        NetworkGravity {
            branch_offsets_j_kg: vec![5.0],
            plenum_potentials_j_kg: vec![5.0],
        },
        || false,
    )
    .unwrap();
    let before = storage(&network);
    let report = network
        .step(Time::new(0.001), StepControl::default(), || false)
        .unwrap();
    let after = storage(&network);
    let gravity = report.gravity.unwrap();
    assert!(report.boundary_change.mass < 0.0);
    assert!(
        (gravity.boundary_centered_potential_energy_change - 0.5 * report.boundary_change.mass)
            .abs()
            < 1e-16
    );
    assert!((after[0] - before[0] - report.boundary_change.mass).abs() < 2e-13);
    assert!(
        (after[1] + after[2]
            - before[1]
            - before[2]
            - report.boundary_change.energy
            - gravity.boundary_centered_potential_energy_change
            - gravity.reference_potential_j_kg * report.boundary_change.mass)
            .abs()
            < 2e-12
    );
    assert!(gravity.conservation_residual.abs() < 2e-12);
}

#[test]
fn g3_common_potential_and_original_formation_datums_do_not_change_flow() {
    let mut base = moving_chain(0.0, 0.0, true);
    let shift = 2.0_f64.powi(40);
    let mut shifted = moving_chain(shift, 1e18, true);
    for _ in 0..35 {
        let a = base
            .step(Time::new(0.001), StepControl::default(), || false)
            .unwrap();
        let b = shifted
            .step(Time::new(0.001), StepControl::default(), || false)
            .unwrap();
        assert_eq!(a.elapsed, b.elapsed);
        assert_eq!(a.before, b.before);
        assert_eq!(a.after, b.after);
        assert_eq!(a.boundary_change, b.boundary_change);
        assert_eq!(a.conservation_residual, b.conservation_residual);
        for (a, b) in base.branches().iter().zip(shifted.branches()) {
            assert_eq!(a.states(), b.states());
        }
        assert_eq!(base.plenums(), shifted.plenums());
        assert_eq!(
            base.accepted_port_exchanges(),
            shifted.accepted_port_exchanges()
        );
        let (ga, gb) = (a.gravity.unwrap(), b.gravity.unwrap());
        assert_eq!(ga.gravity_work, gb.gravity_work);
        assert_eq!(gb.reference_potential_j_kg, shift);
        assert_eq!(
            ga.before_centered_potential_energy,
            gb.before_centered_potential_energy
        );
        assert_eq!(
            ga.after_centered_potential_energy,
            gb.after_centered_potential_energy
        );
        assert_eq!(
            ga.boundary_centered_potential_energy_change,
            gb.boundary_centered_potential_energy_change
        );
        assert_eq!(ga.conservation_residual, gb.conservation_residual);
        assert_eq!(
            gb.absolute_conservation_residual,
            ga.conservation_residual + shift * a.conservation_residual.mass
        );
    }
}

#[test]
fn g0_potential_mismatch_and_unrepresentable_common_datum_are_refused() {
    let gas = crate::dynamics::tests::gas(0.4, 0.0);
    let build = |offset: f64, node: f64| {
        let branch = column(&gas, false, 5, -0.5, 0.0, (2.0, 3.0));
        let plenum = Plenum::new(0.2, branch.hydrostatic_reference_face(5).unwrap().0).unwrap();
        GasNetwork1d::new_with_gravity(
            vec![branch],
            vec![plenum],
            vec![connect(0, BranchEnd::Right, 0)],
            NetworkGravity {
                branch_offsets_j_kg: vec![offset],
                plenum_potentials_j_kg: vec![node],
            },
            || false,
        )
    };
    assert!(build(0.0, 0.5).is_ok());
    assert!(build(0.0, 0.6).is_err());
    assert!(build(f64::NAN, 0.5).is_err());
    // Absolute-datum-scaled tolerances would silently erase this real half-J/kg
    // head: both naive absolute face/node values round to the same huge number.
    let huge = 2.0_f64.powi(60);
    assert_eq!(huge + 0.5, huge);
    assert!(build(huge, huge).is_err());
    let branch = column(&gas, false, 5, -0.5, 0.0, (2.0, 3.0));
    assert!(GasNetwork1d::new(vec![branch], vec![], vec![], || false).is_err());
    let branch = column(&gas, false, 5, -0.5, 0.0, (2.0, 3.0));
    assert!(
        GasNetwork1d::new_with_gravity(
            vec![branch],
            vec![],
            vec![],
            NetworkGravity {
                branch_offsets_j_kg: vec![],
                plenum_potentials_j_kg: vec![]
            },
            || false
        )
        .is_err()
    );
}

#[test]
fn g4_final_gravity_report_cancellation_preserves_every_published_owner() {
    let mut probe = moving_chain(0.0, 0.0, true);
    let mut cancelled_network = moving_chain(0.0, 0.0, true);
    for network in [&mut probe, &mut cancelled_network] {
        network
            .step(Time::new(0.001), StepControl::default(), || false)
            .unwrap();
    }
    let states: Vec<_> = cancelled_network
        .branches()
        .iter()
        .map(|b| b.states().to_vec())
        .collect();
    let plenums = cancelled_network.plenums().to_vec();
    let report = cancelled_network.last_report();
    let ports = cancelled_network
        .accepted_port_exchanges()
        .unwrap()
        .to_vec();
    let branches = cancelled_network
        .accepted_branch_reports()
        .unwrap()
        .to_vec();
    let time = cancelled_network.time();
    let potential = cancelled_network.gravitational_potential_energy();
    let mut polls = 0;
    probe
        .step(Time::new(0.001), StepControl::default(), || {
            polls += 1;
            false
        })
        .unwrap();
    let mut observed = 0;
    let result = cancelled_network.step(Time::new(0.001), StepControl::default(), || {
        observed += 1;
        observed == polls
    });
    assert!(matches!(result, Err(Error::Cancelled)));
    assert_eq!(observed, polls);
    for (branch, old) in cancelled_network.branches().iter().zip(states) {
        assert_eq!(branch.states(), old);
    }
    assert_eq!(cancelled_network.plenums(), plenums);
    assert_eq!(cancelled_network.last_report(), report);
    assert_eq!(cancelled_network.accepted_port_exchanges().unwrap(), ports);
    assert_eq!(
        cancelled_network.accepted_branch_reports().unwrap(),
        branches
    );
    assert_eq!(cancelled_network.time(), time);
    assert_eq!(
        cancelled_network.gravitational_potential_energy(),
        potential
    );
}
