//! Original finite ducts exchange continuously at pressure-matched junctions.
use super::*;
use crate::dg_port::DgBranch;
use fs_duct::dg::{DgBoundary, DgDuct, DgMesh, DgSidePort, DgState, PortId};
use fs_duct::{Duct, Segment};
use fs_material::gas::GasSpec;
use fs_vfit::mimo::{
    MimoModel,
    load::{MimoLoadLimits, PreparedMimoLoad},
};
use fs_vfit::model::PoleTerm;

const INLET: PortId = PortId(7);
const OUTLET: PortId = PortId(9);
const FIRST: PortId = PortId(11);
const SECOND: PortId = PortId(19);
const PORTS: [PortId; 3] = [SECOND, OUTLET, FIRST];
const DIAGONAL: [f64; 3] = [2e7, 3e7, 4e7];
const MUTUAL: f64 = 2e6;
const RESIDUE_FACTOR: f64 = 100.0;
const DECAY: f64 = 400.0;

fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}
fn cylinder(radius: f64, length: f64, elements: usize, inlet: PortId, outlet: PortId) -> DgDuct {
    DgDuct::new(
        &Duct {
            segments: vec![Segment::Cylinder { radius, length }],
        },
        &gas(),
        DgMesh { elements, order: 2 },
        DgBoundary::Port(inlet),
        DgBoundary::Port(outlet),
    )
    .unwrap()
}
fn branched_load(rate: u32, wall: bool, reverse_branches: bool) -> DgPortLoad {
    let dt = 1.0 / f64::from(rate);
    let main = cylinder(0.001, 0.05, 3, INLET, OUTLET)
        .with_side_ports(&[
            DgSidePort {
                id: FIRST,
                downstream_element: 1,
            },
            DgSidePort {
                id: SECOND,
                downstream_element: 2,
            },
        ])
        .unwrap();
    let first = cylinder(0.001, 0.012, 1, PortId(0), PortId(1));
    let first = if wall {
        first.with_viscothermal(8, [50.0, 2000.0], 0.01).unwrap()
    } else {
        first
    };
    let mut branches = vec![
        DgBranch {
            port: FIRST,
            duct: first,
        },
        DgBranch {
            port: SECOND,
            duct: cylinder(0.0012, 0.018, 1, PortId(0), PortId(1)),
        },
    ];
    if reverse_branches {
        branches.reverse();
    }
    // D = diag(d_i) + mu 11^T is SPD. The original mutual memory obeys
    // z'=-a*z+u, p=D*u+k*D*z, with known physical storage k*z^T*D*z/2.
    let d: Vec<_> = (0..3)
        .flat_map(|i| (0..3).map(move |j| MUTUAL + if i == j { DIAGONAL[i] } else { 0.0 }))
        .collect();
    let model = MimoModel {
        ports: 3,
        poles: vec![PoleTerm::Real {
            pole: -DECAY,
            residue: 0.0,
        }],
        residues: vec![d.iter().map(|x| C64::from_re(RESIDUE_FACTOR * x)).collect()],
        d,
        e: vec![0.0; 9],
        reciprocal: true,
    };
    let terminal = PreparedMimoLoad::new(
        &model,
        dt,
        MimoLoadLimits {
            band_rad_s: [1.0, core::f64::consts::PI / dt],
            sampled_passivity_points: 129,
            max_states: 3,
            max_memory_bytes: 1 << 20,
        },
    )
    .unwrap();
    DgPortLoad::new_with_terminal_branches(main, dt, 1e-11, 256, PORTS.to_vec(), terminal, branches)
        .unwrap()
}
fn branched(rate: u32, wall: bool, reverse_branches: bool) -> SharedCavityDg {
    cavity(rate, wall, 1e6)
        .with_exact_downstream_duct(
            branched_load(rate, wall, reverse_branches),
            budget(),
            &CancelGate::new_clock_free(),
        )
        .unwrap()
}
fn pressure(duct: &DgDuct, state: &DgState, port: PortId) -> f64 {
    let original = duct.continuous_operator();
    let (_, column) = original.ports().iter().find(|(id, _)| *id == port).unwrap();
    column
        .iter()
        .zip(original.energy_coordinates(state).unwrap())
        .map(|(b, x)| b * x)
        .sum()
}
fn same_network(a: &SharedCavityDg, b: &SharedCavityDg) {
    same_state(a, b);
    for original in a.load().branches() {
        let other = b
            .load()
            .branches()
            .iter()
            .find(|x| x.port() == original.port())
            .unwrap();
        let av = original.duct().continuous_operator();
        let bv = other.duct().continuous_operator();
        let ax = av.energy_coordinates(original.state()).unwrap();
        let bx = bv.energy_coordinates(other.state()).unwrap();
        assert_eq!(ax.len(), bx.len());
        for (&x, &y) in ax.iter().zip(bx) {
            close(x, y, 2e-16);
        }
        close(
            original.state().stored_energy_j(),
            other.state().stored_energy_j(),
            2e-18,
        );
    }
    for (x, y) in a
        .load()
        .terminal()
        .unwrap()
        .pole_states()
        .iter()
        .zip(b.load().terminal().unwrap().pole_states())
    {
        close(x.re, y.re, 2e-21);
        close(x.im, y.im, 2e-21);
    }
}

#[test]
fn finite_branch_trajectories_and_true_junction_work_survive_three_clocks_and_owner_order() {
    let gate = CancelGate::new_clock_free();
    let mut whole = branched(24000, true, false);
    let mut halves = branched(48000, true, true);
    let mut quarters = branched(96000, true, false);
    let mut covariance = 0.0_f64;
    for _ in 0..16 {
        let frame = whole.step(&[0.001], 0.03, &gate).unwrap().clone();
        let junctions = whole.load().last_branch_frames().to_vec();
        for (system, count) in [(&mut halves, 2), (&mut quarters, 4)] {
            let mut sums = [0.0; 4];
            let mut work = [0.0; 2];
            let mut volume = [0.0; 2];
            for _ in 0..count {
                let step = system.step(&[0.001], 0.03, &gate).unwrap();
                for (total, value) in sums.iter_mut().zip([
                    step.authored_external_work_j,
                    step.inferred_total_dissipation_j,
                    step.inlet_work_j,
                    step.terminal_work_j,
                ]) {
                    *total += value;
                }
                for (i, reference) in junctions.iter().enumerate() {
                    let j = system
                        .load()
                        .last_branch_frames()
                        .iter()
                        .find(|j| j.port == reference.port)
                        .unwrap();
                    work[i] += j.bore_work_j;
                    volume[i] += j.junction_flow_m3_s * system.load().propagator().dt();
                }
            }
            same_network(&whole, system);
            for (actual, expected) in sums.into_iter().zip([
                frame.authored_external_work_j,
                frame.inferred_total_dissipation_j,
                frame.inlet_work_j,
                frame.terminal_work_j,
            ]) {
                close(actual, expected, 2e-16);
            }
            for (i, reference) in junctions.iter().enumerate() {
                close(work[i], reference.bore_work_j, 2e-18);
                close(volume[i], reference.junction_flow_m3_s / 24000.0, 2e-21);
            }
        }
        for junction in &junctions {
            close(junction.bore_work_j, junction.branch_work_j, 1e-22);
            close(
                junction.bore_port.conjugate_rate,
                junction.branch_inlet.conjugate_rate,
                2e-8,
            );
            assert!(junction.duct.ledger.balance_checked);
            covariance = covariance.max(
                (junction.bore_work_j
                    - junction.bore_port.conjugate_rate * junction.junction_flow_m3_s / 24000.0)
                    .abs(),
            );
            let branch = whole
                .load()
                .branches()
                .iter()
                .find(|b| b.port() == junction.port)
                .unwrap();
            close(
                pressure(whole.load().duct(), whole.load().state(), junction.port),
                pressure(branch.duct(), branch.state(), PortId(0)),
                2e-8,
            );
        }
        assert!(frame.exact_joint_evolution && frame.duct.ledger.balance_checked);
        assert!(frame.balance_residual_j.abs() <= frame.balance_tolerance_j);
    }
    assert!(
        covariance > 1e-18,
        "junction product-of-means error {covariance:e}"
    );
    let wall = whole
        .load()
        .branches()
        .iter()
        .find(|b| b.port() == FIRST)
        .unwrap();
    assert!(wall.state().wall_stored_energy_j() > 0.0);
    assert!(wall.state().cumulative_dissipated_j() > 0.0);
    assert!(
        whole
            .load()
            .terminal()
            .unwrap()
            .pole_states()
            .iter()
            .any(|z| z.re != 0.0)
    );
}

// Independent physical terminal reference. The Sherman-Morrison formula here
// inverts the authored diag(d_i)+mu*11^T; it does not use TerminalSolve, the
// prepared generator, or the runtime's integrated work/loss matrices.
fn terminal_storage_and_heat(system: &SharedCavityDg) -> (f64, f64) {
    let load = system.load();
    let z: Vec<_> = load
        .terminal()
        .unwrap()
        .pole_states()
        .iter()
        .map(|z| z.re)
        .collect();
    let zsum: f64 = z.iter().sum();
    let dz: Vec<_> = (0..3).map(|i| DIAGONAL[i] * z[i] + MUTUAL * zsum).collect();
    let pressures: Vec<_> = PORTS
        .iter()
        .map(|&id| {
            if let Some(branch) = load.branches().iter().find(|b| b.port() == id) {
                pressure(branch.duct(), branch.state(), PortId(1))
            } else {
                pressure(load.duct(), load.state(), id)
            }
        })
        .collect();
    let rhs: Vec<_> = pressures
        .iter()
        .zip(&dz)
        .map(|(p, history)| p - RESIDUE_FACTOR * history)
        .collect();
    let correction = MUTUAL * (0..3).map(|i| rhs[i] / DIAGONAL[i]).sum::<f64>()
        / (1.0 + MUTUAL * DIAGONAL.iter().map(|d| 1.0 / d).sum::<f64>());
    let u: Vec<_> = (0..3)
        .map(|i| (rhs[i] - correction) / DIAGONAL[i])
        .collect();
    let zdotdz: f64 = z.iter().zip(&dz).map(|(x, y)| x * y).sum();
    let resistance = (0..3).map(|i| DIAGONAL[i] * u[i] * u[i]).sum::<f64>()
        + MUTUAL * u.iter().sum::<f64>().powi(2);
    (
        0.5 * RESIDUE_FACTOR * zdotdz,
        resistance + DECAY * RESIDUE_FACTOR * zdotdz,
    )
}

#[test]
fn original_mutual_terminal_heat_and_each_branch_balance_close_before_atomic_publication() {
    let gate = CancelGate::new_clock_free();
    let rate = 24000;
    let subdivisions = 128;
    let mut coarse = branched(rate, true, false);
    let mut fine = branched(rate * subdivisions, true, false);
    let initial = coarse.stored_energy_j().unwrap();
    let (mut authored, mut physical_loss, mut terminal_work) = (0.0, 0.0, 0.0);
    let mut heats = [0.0; 3];
    let mut fine_terminal_work = 0.0;
    for _ in 0..8 {
        let step = coarse.step(&[0.001], 0.03, &gate).unwrap();
        authored += step.authored_external_work_j;
        physical_loss += step.inferred_total_dissipation_j;
        terminal_work += step.terminal_work_j;
        let mut quadrature = [terminal_storage_and_heat(&fine).1; 3];
        for substep in 1..=subdivisions {
            fine_terminal_work += fine.step(&[0.001], 0.03, &gate).unwrap().terminal_work_j;
            let power = terminal_storage_and_heat(&fine).1;
            for (rung, stride) in [4, 2, 1].into_iter().enumerate() {
                if substep % stride == 0 {
                    let i = substep / stride;
                    let weight = if substep == subdivisions {
                        1.0
                    } else if i % 2 == 0 {
                        2.0
                    } else {
                        4.0
                    };
                    quadrature[rung] += weight * power;
                }
            }
        }
        for (rung, stride) in [4, 2, 1].into_iter().enumerate() {
            heats[rung] +=
                quadrature[rung] * f64::from(stride) / (3.0 * f64::from(rate * subdivisions));
        }
        same_network(&coarse, &fine);
    }
    let energy = terminal_storage_and_heat(&coarse).0;
    let errors = heats.map(|heat| (terminal_work - energy - heat).abs());
    eprintln!(
        "branched terminal physical heat errors32/64/128 {errors:?}; work {terminal_work:e}, stored {energy:e}, heat {:e}",
        heats[2]
    );
    assert!(errors[1] < errors[0] && errors[2] < errors[1], "{errors:?}");
    assert!(errors[2] < 1e-18 + 2e-6 * terminal_work.abs().max(energy + heats[2]));
    close(terminal_work, fine_terminal_work, 2e-18);
    assert!(
        (coarse.stored_energy_j().unwrap() + energy + physical_loss + heats[2]
            - initial
            - authored)
            .abs()
            < 2e-14 + errors[2]
    );
    assert!(energy > 0.0 && heats[2] > 0.0);

    let frame = coarse.last_frame().unwrap().clone();
    let duct = coarse.load().state().clone();
    let branches: Vec<_> = coarse
        .load()
        .branches()
        .iter()
        .map(|b| b.state().clone())
        .collect();
    let junctions = coarse.load().last_branch_frames().to_vec();
    let poles = coarse.load().terminal().unwrap().pole_states().to_vec();
    let state = endpoint(&coarse);
    let memory = coarse.cavity().network().material_memory().to_vec();
    assert!(coarse.step(&[1e20], 0.0, &gate).is_err());
    let cancel = CancelGate::new_clock_free();
    cancel.request();
    assert!(coarse.step(&[0.0], 0.0, &cancel).is_err());
    assert_eq!(coarse.last_frame().unwrap(), &frame);
    assert_eq!(coarse.load().state(), &duct);
    for (branch, original) in coarse.load().branches().iter().zip(branches) {
        assert_eq!(branch.state(), &original);
    }
    assert_eq!(coarse.load().last_branch_frames(), junctions);
    assert_eq!(coarse.load().terminal().unwrap().pole_states(), poles);
    assert_eq!(endpoint(&coarse), state);
    assert_eq!(coarse.cavity().network().material_memory(), memory);
    // The untouched candidate can continue identically after late refusal.
    coarse.step(&[0.001], 0.03, &gate).unwrap();
    for _ in 0..subdivisions {
        fine.step(&[0.001], 0.03, &gate).unwrap();
    }
    same_network(&coarse, &fine);
}
