//! G3 adapter qualification against the shared receiver on the actual outlet flow.
use super::*;
use fs_exec::CancelGate;

const SOURCE: &str = include_str!("../../../examples/plate-valve-exterior.performance");

fn fixture() -> file::PlateValvePerformance {
    file::PlateValvePerformance::from_bytes(SOURCE.as_bytes(), 113, &CancelGate::new_clock_free())
        .unwrap()
}

fn accepted() -> AperturePerformance {
    let loaded = fixture();
    let load = loaded.info().radiation_load.unwrap();
    let old = loaded.into_renderer();
    AperturePerformance::new_with_observation_policy(
        old.system,
        old.observation,
        old.schedule,
        old.config,
        old.exterior_gas,
        Some(load),
        ObservationPolicy::Acceptance,
    )
    .unwrap()
}

#[test]
fn aperture_performance_observer_matches_owner() {
    let mut source = accepted();
    let (CoupledAperture::Network(network), ApertureObservation::NetworkBaffled { node, receiver }) =
        (&source.system, source.observation)
    else {
        panic!("fixture must have a physical receiver")
    };
    let radius = network
        .spec()
        .sections
        .iter()
        .find(|s| s.nodes.contains(&node))
        .unwrap()
        .radius_m;
    let mut owner = BaffledPressure::circular_outlet(
        radius,
        source.config.sample_rate_hz,
        receiver,
        source.baffled_medium().unwrap(),
    )
    .unwrap();
    let mut peak = 0.0_f64;
    let start = std::time::Instant::now();
    for _ in 0..2048 {
        let mut observed = [0.0];
        source.block(&mut observed).unwrap();
        let CoupledAperture::Network(network) = &source.system else {
            unreachable!()
        };
        let flow = network.node_frame(node).unwrap().net_flow_into_node_m3_s;
        assert_eq!(
            observed[0].to_bits(),
            owner.step(&[flow]).unwrap().to_bits()
        );
        peak = peak.max(observed[0].abs());
    }
    assert!(peak > 0.0);
    println!(
        "{{\"case\":\"aperture_performance_observer_matches_owner\",\"inputs\":\"authored plate-valve-exterior.performance; accepted terminal volume flow\",\"observed\":0,\"err\":0,\"tol\":0,\"verdict\":\"pass\",\"elapsed_s\":{}}}",
        start.elapsed().as_secs_f64()
    );
}

#[test]
fn reaction_receiver_gas_mismatch_refuses() {
    let loaded = fixture();
    let load = loaded.info().radiation_load.unwrap();
    let old = loaded.into_renderer();
    let external = *old.exterior_gas().unwrap();
    let changed = GasState::try_new_moist_air(293.15, external.pressure, 0.0).unwrap();
    let error = AperturePerformance::new_with_observation_policy(
        old.system,
        old.observation,
        old.schedule,
        old.config,
        Some(changed),
        Some(load),
        ObservationPolicy::Acceptance,
    )
    .err()
    .unwrap();
    assert!(
        error.to_string().contains("reaction/receiver gas mismatch"),
        "{error}"
    );

    let loaded = fixture();
    let old = loaded.into_renderer();
    let error = AperturePerformance::new_with_observation_policy(
        old.system,
        old.observation,
        old.schedule,
        old.config,
        old.exterior_gas,
        None,
        ObservationPolicy::Acceptance,
    )
    .err()
    .unwrap();
    assert!(
        error
            .to_string()
            .contains("actual matched reaction-load binding")
    );
}

#[test]
fn observer_bitwise_repeat() {
    let (mut a, mut b) = (accepted(), accepted());
    let (mut x, mut y) = ([0.0; 2048], [0.0; 2048]);
    for chunk in x.chunks_mut(113) {
        a.block(chunk).unwrap();
    }
    for chunk in y.chunks_mut(37) {
        b.block(chunk).unwrap();
    }
    assert!(x.iter().any(|v| *v != 0.0));
    assert!(x.iter().zip(y).all(|(a, b)| a.to_bits() == b.to_bits()));
    assert_eq!(a.observation_policy(), ObservationPolicy::Acceptance);
}
