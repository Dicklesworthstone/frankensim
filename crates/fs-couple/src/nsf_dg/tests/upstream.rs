//! G0/G3: both original acoustic loads participate in one resolved interval.
use super::*;

fn coupled(high: bool, upstream_boundary: DgBoundary) -> NsfDg<D3q27> {
    let gate = CancelGate::new_clock_free();
    let mut owner = system(100000.0, high, 1, DgBoundary::Matched);
    owner.limits.incoming_bracket_m_s = [-0.5, 0.5];
    owner.limits.maximum_evaluations = 128;
    owner.limits.maximum_cell_steps = CELLS * 128;
    let other = load(
        owner.zone(),
        *owner.reference_gas(),
        1.0,
        1.0,
        upstream_boundary,
    );
    owner.with_upstream_load(other, &gate).unwrap()
}

#[test]
fn g0_two_sided_pressure_and_source_work_replay_both_original_owners() {
    let gate = CancelGate::new_clock_free();
    for high in [false, true] {
        let mut coupled = coupled(high, DgBoundary::Closed);
        let reference = *coupled.reference_gas();
        let mut fluid = zone(GasSpec::dry_air_ussa1976(), reference, high, 1);
        let mut main = load(&fluid, reference, 1.0, 1.0, DgBoundary::Matched);
        let mut upstream = load(&fluid, reference, 1.0, 1.0, DgBoundary::Closed);
        for rise in [2.0, 2.0, -1.0, 0.0] {
            let frame = coupled
                .step_with_upstream_pressure_rise(rise, &gate)
                .unwrap();
            let source = frame.upstream.as_ref().unwrap();
            assert_eq!(source.face, usize::from(!high));
            assert_eq!(frame.drive, NsfDgDrive::UpstreamPressureRise(rise));
            assert!(frame.evaluations <= 128);
            let mut incoming = [source.incoming_characteristic_m_s; 2];
            incoming[usize::from(high)] = frame.incoming_characteristic_m_s;
            let raw = fluid
                .preview_characteristic_interval(incoming, 1, CELLS, || false)
                .unwrap();
            assert_eq!(raw.boundary_exchange(), &frame.faces);
            raw.accept();
            main.preview_step(frame.inlet_flow_m3_s).unwrap().accept();
            upstream
                .preview_step(source.inlet_flow_m3_s)
                .unwrap()
                .accept();
            assert_eq!(fluid.state_hash(), coupled.zone().state_hash());
            assert_eq!(main.state(), coupled.load().state());
            assert_eq!(upstream.state(), coupled.upstream_load().unwrap().state());
            let face = frame.faces[source.face];
            assert!(
                (face.mean_acoustic_pressure_pa - source.inlet_pressure_pa - rise).abs() < 1e-6
            );
            assert!(frame.pressure_residual_pa.abs() < 1e-6);
            let work = -rise * face.normal_volume_flow_m3_s / f64::from(RATE);
            assert!((source.source_work_j - work).abs() < 1e-12 * work.abs().max(1e-24));
            // The imposed series effort has its own independent energy port.
            let reconstructed = face.acoustic_pressure_work_out_j + source.source_work_j
                - source.inlet_pressure_pa * source.inlet_flow_m3_s / f64::from(RATE);
            let work_scale = face.acoustic_pressure_work_out_j.abs()
                + source.source_work_j.abs()
                + (source.inlet_pressure_pa * source.inlet_flow_m3_s / f64::from(RATE)).abs();
            assert!(
                (reconstructed - source.acoustic_work_residual_j - source.port_work_residual_j)
                    .abs()
                    <= 32.0 * f64::EPSILON * work_scale.max(1e-30)
            );
        }
    }
}

#[test]
fn g3_returning_upstream_waves_change_source_flow_and_the_bore() {
    let gate = CancelGate::new_clock_free();
    let mut closed = coupled(true, DgBoundary::Closed);
    let mut matched = coupled(true, DgBoundary::Matched);
    let mut source_difference: f64 = 0.0;
    let mut bore_difference: f64 = 0.0;
    for _ in 0..96 {
        let a = closed.step_with_upstream_pressure_rise(2.0, &gate).unwrap();
        let b = matched
            .step_with_upstream_pressure_rise(2.0, &gate)
            .unwrap();
        source_difference = source_difference
            .max((a.faces[0].normal_volume_flow_m3_s - b.faces[0].normal_volume_flow_m3_s).abs());
        bore_difference = bore_difference.max((a.inlet_pressure_pa - b.inlet_pressure_pa).abs());
    }
    assert!(source_difference > 1e-12);
    assert!(bore_difference > 1e-9);
    assert_ne!(closed.zone().state_hash(), matched.zone().state_hash());
    assert_ne!(closed.load().state(), matched.load().state());
}

#[test]
fn g4_two_sided_budget_domain_and_cancellation_refuse_before_three_way_publication() {
    let gate = CancelGate::new_clock_free();
    let mut owner = coupled(true, DgBoundary::Closed);
    owner.step_with_upstream_pressure_rise(2.0, &gate).unwrap();
    let before = (
        owner.zone().state_hash(),
        owner.load().state().clone(),
        owner.upstream_load().unwrap().state().clone(),
        owner.accepted_steps(),
    );
    for rise in [f64::NAN, 1e10] {
        assert!(owner.step_with_upstream_pressure_rise(rise, &gate).is_err());
        assert_eq!(
            before,
            (
                owner.zone().state_hash(),
                owner.load().state().clone(),
                owner.upstream_load().unwrap().state().clone(),
                owner.accepted_steps()
            )
        );
    }
    assert!(
        owner
            .step_with_boundary(NsfBoundaryDrive::PressureOffset(2.0), &gate)
            .is_err()
    );
    let budget = owner.limits.maximum_evaluations;
    owner.limits.maximum_evaluations = 3;
    assert!(owner.step_with_upstream_pressure_rise(2.0, &gate).is_err());
    owner.limits.maximum_evaluations = budget;
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        owner.step_with_upstream_pressure_rise(2.0, &cancelled),
        Err(NsfDgError::Cancelled)
    ));
    assert_eq!(
        before,
        (
            owner.zone().state_hash(),
            owner.load().state().clone(),
            owner.upstream_load().unwrap().state().clone(),
            owner.accepted_steps()
        )
    );
    let mut replay = coupled(true, DgBoundary::Closed);
    replay.step_with_upstream_pressure_rise(2.0, &gate).unwrap();
    let a = owner.step_with_upstream_pressure_rise(2.0, &gate).unwrap();
    let b = replay.step_with_upstream_pressure_rise(2.0, &gate).unwrap();
    assert_eq!(a.inlet_pressure_pa, b.inlet_pressure_pa);
    assert_eq!(owner.zone().state_hash(), replay.zone().state_hash());
    assert_eq!(
        owner.upstream_load().unwrap().state(),
        replay.upstream_load().unwrap().state()
    );
}

#[test]
fn g0_upstream_requires_original_gas_area_clock_and_cold_topology() {
    let gate = CancelGate::new_clock_free();
    for (area, clock) in [(2.0, 1.0), (1.0, 2.0)] {
        let owner = system(100000.0, true, 1, DgBoundary::Matched);
        let other = load(
            owner.zone(),
            *owner.reference_gas(),
            area,
            clock,
            DgBoundary::Closed,
        );
        assert!(owner.with_upstream_load(other, &gate).is_err());
    }
    let mut owner = system(100000.0, true, 1, DgBoundary::Matched);
    assert!(owner.step_with_upstream_pressure_rise(2.0, &gate).is_err());
    owner.step(0.0, &gate).unwrap();
    let other = load(
        owner.zone(),
        *owner.reference_gas(),
        1.0,
        1.0,
        DgBoundary::Closed,
    );
    assert!(owner.with_upstream_load(other, &gate).is_err());
}
