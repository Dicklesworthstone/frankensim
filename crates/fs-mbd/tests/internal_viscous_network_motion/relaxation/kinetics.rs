//! Held kinetics reach both native moving source halves without rebinding memory.
use super::*;
use fs_mbd::spatial_loads::internal_motion::viscous::network::relaxation::{
    HeldRelaxationKinetics, step_relaxation_network_with_loads_and_kinetics,
};

#[test]
fn g1_moving_held_kinetics_matches_original_effective_laws_and_is_atomic() {
    let (states, properties, viscous) = full();
    let links = memory_links(viscous);
    let original = links.clone();
    let history = memory();
    let preserved = history.clone();
    let d = drives();
    let efforts = [
        HeldInternalEffort::World(d[0]),
        HeldInternalEffort::MaterialFirst(d[1]),
        HeldInternalEffort::World(d[2]),
    ];
    let viscosity = [0.3, 1.6, 0.0];
    let times = vec![vec![0.4], vec![2.0], vec![4.0]];
    let kinetics = HeldRelaxationKinetics {
        viscosity_scales: &viscosity,
        relaxation_time_scales: &times,
    };
    let effective = core::array::from_fn::<_, 3, _>(|e| {
        let mut link = links[e].clone();
        link.viscous.damping = link.viscous.damping.scaled(viscosity[e]).unwrap();
        link.branches[0] = SpatialRelaxationBranch::from_gram_rows(
            *links[e].branches[0].gram_rows(),
            links[e].branches[0].relaxation_time_s() * times[e][0],
        )
        .unwrap();
        link
    });
    let wrenches = [
        Wrench {
            force_world: Vec3::new(0.1, 0.2, -0.1),
            torque_body: Vec3::new(0.01, 0.0, 0.02),
        },
        Wrench {
            force_world: Vec3::new(-0.1, 0.05, 0.2),
            torque_body: Vec3::new(0.0, 0.03, -0.01),
        },
        Wrench {
            force_world: Vec3::new(0.02, -0.03, 0.04),
            torque_body: Vec3::new(-0.01, 0.02, 0.01),
        },
    ];
    let loads = HeldExternalLoads {
        wrenches: &wrenches,
        gravity: Gravity::new(Vec3::new(0.0, 0.0, -0.2)).unwrap(),
    };
    let expected = step_relaxation_network_with_loads(
        &states,
        &properties,
        &effective,
        &history,
        &efforts,
        loads,
        0.002,
        memory_controls(),
        || false,
    )
    .unwrap();
    let mut polls = 0;
    let actual = step_relaxation_network_with_loads_and_kinetics(
        &states,
        &properties,
        &links,
        &history,
        &efforts,
        kinetics,
        loads,
        0.002,
        memory_controls(),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    assert_eq!(actual, expected);
    near(actual.work_heat_residual_j, 0.0, 3e-10);
    assert!(actual.connector_dissipation_j.iter().all(|q| *q > 0.0));
    let mut seen = 0;
    assert_eq!(
        step_relaxation_network_with_loads_and_kinetics(
            &states,
            &properties,
            &links,
            &history,
            &efforts,
            kinetics,
            loads,
            0.002,
            memory_controls(),
            || {
                seen += 1;
                seen == polls
            }
        ),
        Err(RelaxationNetworkMotionError::Cancelled)
    );
    let invalid = vec![vec![0.4], vec![], vec![4.0]];
    assert!(
        step_relaxation_network_with_loads_and_kinetics(
            &states,
            &properties,
            &links,
            &history,
            &efforts,
            HeldRelaxationKinetics {
                relaxation_time_scales: &invalid,
                ..kinetics
            },
            loads,
            0.002,
            memory_controls(),
            || false
        )
        .is_err()
    );
    assert_eq!(
        step_relaxation_network_with_loads_and_kinetics(
            &states,
            &properties,
            &links,
            &history,
            &efforts,
            kinetics,
            loads,
            0.002,
            memory_controls(),
            || false
        )
        .unwrap(),
        expected
    );
    assert_eq!(links, original);
    assert_eq!(history, preserved);
}
