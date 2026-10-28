use super::*;
use crate::dynamics::tests::gas;
use crate::dynamics::{Boundary, Mesh1d, StepControl};
use fs_qty::Length;

fn walls() -> Boundaries {
    Boundaries::Ends {
        left: Boundary::Reflecting,
        right: Boundary::Reflecting,
    }
}
fn branch(
    gas: &FrozenGas,
    n: usize,
    boundaries: Boundaries,
    f: impl Fn(f64) -> [f64; 3],
) -> Euler1d {
    let states = (0..n)
        .map(|i| {
            let [rho, u, p] = f((i as f64 + 0.5) / n as f64);
            gas.conserved(rho, u, p).unwrap()
        })
        .collect();
    Euler1d::new(
        gas.clone(),
        Length::new(1.0 / n as f64),
        0.1,
        states,
        boundaries,
        || false,
    )
    .unwrap()
}
fn connected() -> GasNetwork1d {
    let gas = gas(0.2, 0.0);
    let a = branch(&gas, 12, walls(), |x| {
        [1.0, 0.0, if x > 0.7 { 1.1 } else { 1.0 }]
    });
    let b = branch(&gas, 12, walls(), |_| [1.0, 0.0, 1.0]);
    GasNetwork1d::new(
        vec![a, b],
        vec![Plenum::new(0.02, gas.conserved(1.0, 0.0, 1.0).unwrap()).unwrap()],
        vec![
            PlenumConnection {
                branch: 0,
                end: BranchEnd::Right,
                plenum: 0,
            },
            PlenumConnection {
                branch: 1,
                end: BranchEnd::Left,
                plenum: 0,
            },
        ],
        || false,
    )
    .unwrap()
}

#[test]
fn an_isolated_branch_retains_the_original_euler_step_and_field() {
    let gas = gas(0.2, -5.0);
    let make = || {
        branch(&gas, 32, Boundaries::Periodic, |x| {
            [
                1.0 + 0.1 * fs_math::det::sin(2.0 * core::f64::consts::PI * x),
                0.3,
                1.0,
            ]
        })
    };
    let mut direct = make();
    let mut network = GasNetwork1d::new(vec![make()], Vec::new(), Vec::new(), || false).unwrap();
    assert!(network.accepted_branch_reports().is_none());
    assert!(network.accepted_port_exchanges().is_none());
    for _ in 0..30 {
        let original = direct
            .step(Time::new(0.003), StepControl::default(), || false)
            .unwrap();
        let common = network
            .step(Time::new(0.003), StepControl::default(), || false)
            .unwrap();
        assert_eq!(network.branches()[0].states(), direct.states());
        assert_eq!(network.time(), direct.time());
        assert_eq!(network.accepted_branch_reports().unwrap(), &[original]);
        assert_eq!(common.elapsed, original.elapsed);
        assert_eq!(common.maximum_cfl, original.maximum_cfl);
        assert_eq!(common.before.mass, original.before.mass);
        assert_eq!(common.after.energy, original.after.energy);
        assert_eq!(common.boundary_change, NetworkTotals::default());
        assert!(network.accepted_port_exchanges().unwrap().is_empty());
    }
}

#[test]
fn compact_network_preserves_rest_across_density_contacts_and_tapered_branches() {
    let source = gas(0.3, 0.0);
    let mut domain = source.domain();
    // This roundoff-scale equilibrium check needs a tighter caloric solve than
    // the general fixture's 1e-10 J/kg absolute inversion tolerance.
    domain.absolute_energy_tolerance = 1e-14;
    domain.relative_energy_tolerance = 1e-14;
    let gas = FrozenGas::new(source.source().clone(), domain, || false).unwrap();
    let a = Euler1d::new_on_mesh(
        gas.clone(),
        Mesh1d::new(
            vec![0.0, 0.2, 0.4, 1.0],
            vec![0.1, 0.2, 0.15, 0.1],
            vec![0.03, 0.035, 0.075],
        )
        .unwrap(),
        [1.2, 1.3, 1.4]
            .map(|rho| gas.conserved(rho, 0.0, 1.0).unwrap())
            .to_vec(),
        walls(),
        || false,
    )
    .unwrap();
    let b = branch(&gas, 8, walls(), |_| [0.9, 0.0, 1.0]);
    let plenum = Plenum::new(0.03, gas.conserved(0.8, 0.0, 1.0).unwrap()).unwrap();
    let initial_plenum = plenum.state();
    let mut network = GasNetwork1d::new(
        vec![a, b],
        vec![plenum],
        vec![
            PlenumConnection {
                branch: 0,
                end: BranchEnd::Right,
                plenum: 0,
            },
            PlenumConnection {
                branch: 1,
                end: BranchEnd::Left,
                plenum: 0,
            },
        ],
        || false,
    )
    .unwrap();
    let before = network.totals();
    for _ in 0..20 {
        let report = network
            .step(Time::new(0.003), StepControl::default(), || false)
            .unwrap();
        assert_eq!(report.boundary_change.mass, 0.0);
        assert_eq!(report.boundary_change.energy, 0.0);
        assert!(report.conservation_residual.mass.abs() < 2e-13);
        assert!(report.conservation_residual.energy.abs() < 2e-12);
        for branch in network.branches() {
            for &state in branch.states() {
                let primitive = gas.primitive(state, || false).unwrap();
                assert!(primitive.velocity.abs() < 2e-12);
                assert!(
                    (primitive.pressure - 1.0).abs() < 2e-12,
                    "rest pressure {}, velocity {}, state {state:?}",
                    primitive.pressure,
                    primitive.velocity
                );
            }
        }
    }
    let after = network.totals();
    assert!((after.mass - before.mass).abs() < 2e-13);
    assert!((after.energy - before.energy).abs() < 2e-12);
    assert!((network.plenums()[0].state().density - initial_plenum.density).abs() < 2e-12);
    assert!((network.plenums()[0].state().energy - initial_plenum.energy).abs() < 2e-12);
}

#[test]
fn a_small_plenum_limits_the_common_step_and_every_branch_clock() {
    let gas = gas(0.0, 0.0);
    let make = || branch(&gas, 4, walls(), |_| [1.0, 0.0, 1.0]);
    let mut network = GasNetwork1d::new(
        vec![make(), make()],
        vec![Plenum::new(1e-5, gas.conserved(1.0, 0.0, 1.01).unwrap()).unwrap()],
        vec![
            PlenumConnection {
                branch: 0,
                end: BranchEnd::Right,
                plenum: 0,
            },
            PlenumConnection {
                branch: 1,
                end: BranchEnd::Left,
                plenum: 0,
            },
        ],
        || false,
    )
    .unwrap();
    let report = network
        .step(Time::new(0.1), StepControl::default(), || false)
        .unwrap();
    // The branch-only acoustic limit is O(0.07 s), while the compact
    // volume's summed port rate requires an O(3e-5 s) common increment.
    assert!(report.elapsed.value() > 0.0 && report.elapsed.value() < 1e-4);
    assert!(report.maximum_cfl <= 0.400_000_000_000_01);
    assert_eq!(network.time(), report.elapsed);
    assert_eq!(report.boundary_change.mass, 0.0);
    assert_eq!(report.boundary_change.energy, 0.0);
    for branch in network.branches() {
        assert_eq!(branch.time(), network.time());
    }
    for branch_report in network.accepted_branch_reports().unwrap() {
        assert_eq!(branch_report.elapsed, report.elapsed);
        assert!(branch_report.maximum_cfl < report.maximum_cfl);
    }
    assert!(
        network
            .accepted_port_exchanges()
            .unwrap()
            .iter()
            .all(|exchange| exchange.mass_to_plenum < 0.0)
    );
}

#[test]
fn network_admission_rejects_inconsistent_physical_owners_and_topology() {
    let gas = gas(0.2, 0.0);
    let make = || branch(&gas, 4, walls(), |_| [1.0, 0.0, 1.0]);
    let p = || Plenum::new(0.1, gas.conserved(1.0, 0.0, 1.0).unwrap()).unwrap();
    let connection = PlenumConnection {
        branch: 0,
        end: BranchEnd::Left,
        plenum: 0,
    };
    assert!(matches!(
        GasNetwork1d::new(Vec::new(), Vec::new(), Vec::new(), || false),
        Err(Error::Invalid(_))
    ));
    assert!(matches!(
        Plenum::new(0.0, p().state()),
        Err(Error::Invalid(_))
    ));
    assert!(matches!(
        Plenum::new(0.1, gas.conserved(1.0, 0.1, 1.0).unwrap()),
        Err(Error::Invalid(_))
    ));
    assert!(matches!(
        GasNetwork1d::new(
            vec![make()],
            vec![p()],
            vec![connection, connection],
            || false
        ),
        Err(Error::Invalid(_))
    ));
    assert!(matches!(
        GasNetwork1d::new(
            vec![make()],
            vec![p()],
            vec![PlenumConnection {
                branch: 1,
                ..connection
            }],
            || false
        ),
        Err(Error::Invalid(_))
    ));
    assert!(matches!(
        GasNetwork1d::new(
            vec![make()],
            vec![p()],
            vec![PlenumConnection {
                plenum: 1,
                ..connection
            }],
            || false
        ),
        Err(Error::Invalid(_))
    ));
    let periodic = branch(&gas, 4, Boundaries::Periodic, |_| [1.0, 0.0, 1.0]);
    assert!(matches!(
        GasNetwork1d::new(vec![periodic], vec![p()], vec![connection], || false),
        Err(Error::Invalid(_))
    ));
    let other = crate::dynamics::tests::gas(0.3, 0.0);
    let incompatible = branch(&other, 4, walls(), |_| [1.0, 0.0, 1.0]);
    assert!(matches!(
        GasNetwork1d::new(vec![make(), incompatible], Vec::new(), Vec::new(), || false),
        Err(Error::Invalid(_))
    ));
    let mut advanced = make();
    advanced
        .step(Time::new(0.001), StepControl::default(), || false)
        .unwrap();
    assert!(matches!(
        GasNetwork1d::new(vec![make(), advanced], Vec::new(), Vec::new(), || false),
        Err(Error::Invalid(_))
    ));
    let unsupported = Plenum::new(
        0.1,
        Conserved {
            density: 1.0,
            momentum: 0.0,
            energy: 1e10,
        },
    )
    .unwrap();
    assert!(matches!(
        GasNetwork1d::new(vec![make()], vec![unsupported], vec![connection], || false),
        Err(Error::EnergyDomain { .. })
    ));
    assert!(matches!(
        GasNetwork1d::new(vec![make()], vec![p()], vec![connection], || true),
        Err(Error::Cancelled)
    ));
}

type Snapshot = (
    Vec<Vec<Conserved>>,
    Vec<Conserved>,
    f64,
    Option<NetworkStepReport>,
    Option<Vec<PortExchange>>,
    Option<Vec<StepReport>>,
);
fn snapshot(network: &GasNetwork1d) -> Snapshot {
    (
        network
            .branches()
            .iter()
            .map(|b| b.states().to_vec())
            .collect(),
        network.plenums().iter().map(Plenum::state).collect(),
        network.time().value(),
        network.last_report(),
        network.accepted_port_exchanges().map(<[_]>::to_vec),
        network.accepted_branch_reports().map(<[_]>::to_vec),
    )
}

#[test]
fn cancellation_preserves_every_owner_clock_and_previous_accepted_exchange() {
    let mut probe = connected();
    probe
        .step(Time::new(0.002), StepControl::default(), || false)
        .unwrap();
    let mut count = 0;
    probe
        .step(Time::new(0.002), StepControl::default(), || {
            count += 1;
            false
        })
        .unwrap();
    assert!(count > 100);
    for stop in [1, count / 2, count] {
        let mut network = connected();
        network
            .step(Time::new(0.002), StepControl::default(), || false)
            .unwrap();
        let before = snapshot(&network);
        let mut visits = 0;
        let result = network.step(Time::new(0.002), StepControl::default(), || {
            visits += 1;
            visits >= stop
        });
        assert!(
            matches!(result, Err(Error::Cancelled)),
            "stop={stop}, visits={visits}"
        );
        assert_eq!(snapshot(&network), before);
        for branch in network.branches() {
            assert_eq!(branch.time().value(), before.2);
        }
    }
}

#[test]
fn exhausted_source_retry_and_invalid_controls_publish_no_partial_network() {
    let gas = gas(0.0, 0.0);
    // Incoming Mach~42 gas requires a stopped wall state above the admitted
    // 50 K caloric path. A smaller dt cannot manufacture that missing source.
    let hot_wall = branch(
        &gas,
        4,
        Boundaries::Ends {
            left: Boundary::Transmissive,
            right: Boundary::Reflecting,
        },
        |_| [1.0, 50.0, 1.0],
    );
    let mut network = GasNetwork1d::new(vec![hot_wall], Vec::new(), Vec::new(), || false).unwrap();
    let before = snapshot(&network);
    assert!(matches!(
        network.step(
            Time::new(0.002),
            StepControl {
                maximum_halvings: 0,
                ..StepControl::default()
            },
            || false
        ),
        Err(Error::StepRetries { attempts: 1, .. })
    ));
    assert_eq!(snapshot(&network), before);
    assert!(matches!(
        network.step(Time::new(0.0), StepControl::default(), || false),
        Err(Error::Invalid(_))
    ));
    assert_eq!(snapshot(&network), before);
}
