//! Remote reservoir reflections act on the original nonlinear paired blades.
use super::*;
use fs_couple::upstream_airway::StoredPressureSource;
use fs_vfit::impedance::{ImpedanceState, SeriesImpedance, SeriesImpedanceSpec};

const COMPLIANCE: f64 = 2.0e-9;
const INITIAL_PRESSURE: f64 = 50.0;

fn cold(up_delay: usize, bore_reflection: f64, steps: u64, driven: bool) -> TwoSidedAperture {
    let mut upstream = duct(gas(), up_delay, -1.0);
    if driven {
        upstream.graph.nodes[1] = NetworkNode::Driven;
    }
    let sides = [upstream, duct(gas(), 8, bore_reflection)];
    let (g, z) = combined_inlet_state(&sides).unwrap();
    // This slow, source-authored pair lies well inside the original midpoint
    // one-cent band. The original nonlinear slit, swept flow and pressure forces
    // remain live; there is no prescribed opening or substitute valve law.
    let b = ApertureBlade {
        mass_kg: 1.0e-4,
        stiffness_n_m: 5.0,
        damping_n_s_m: 0.005,
        pressure_area_m2: 1.0e-5,
    };
    let valve = PairedAperture::new(
        PairedApertureSpec {
            rest_opening_m: 4.0e-4,
            width_m: 0.008,
            blades: [b; 2],
            density_kg_m3: g.density,
            impedance_pa_s_m3: z,
            time_step_s: DT,
            max_steps: steps,
            limits: PairedApertureLimits {
                maximum_abs_displacement_m: 0.01,
                maximum_abs_velocity_m_s: 100.0,
                maximum_energy_j: 1.0,
                max_evaluations: 256,
                velocity_tolerance_m_s: 1.0e-14,
                relative_tolerance: 1.0e-12,
            },
        },
        [BladeState::default(); 2],
        contact(),
    )
    .unwrap();
    TwoSidedAperture::new(valve.into(), sides, limits(), &CancelGate::new_clock_free()).unwrap()
}

fn charge(
    model: &TwoSidedAperture,
    compliance: Option<f64>,
    pressure: f64,
) -> StoredPressureSource {
    StoredPressureSource::new(
        SeriesImpedance::new_with_state(
            SeriesImpedanceSpec {
                resistance_pa_s_m3: 0.0,
                inertance_pa_s2_m3: 0.0,
                compliance_m3_pa: compliance,
            },
            model.represented_sections(ValveSide::Upstream)[0].impedance_pa_s_m3,
            DT,
            ImpedanceState {
                inertive_flow_m3_s: 0.0,
                compliance_pressure_pa: pressure,
            },
        )
        .unwrap(),
    )
}
fn charged(up_delay: usize, reflection: f64, steps: u64) -> TwoSidedAperture {
    let model = cold(up_delay, reflection, steps, true);
    let source = charge(&model, Some(COMPLIANCE), INITIAL_PRESSURE);
    model
        .with_stored_pressure_source(1, source, &CancelGate::new_clock_free())
        .unwrap()
}
fn source_pressure(model: &TwoSidedAperture) -> f64 {
    model
        .stored_pressure_source()
        .unwrap()
        .1
        .impedance()
        .state()
        .compliance_pressure_pa
}
fn close(a: f64, b: f64, scale: f64) {
    assert!(
        (a - b).abs() <= 2.0e-11 * scale.max(a.abs()).max(b.abs()),
        "{a:.17e} != {b:.17e}"
    );
}

#[test]
fn g1_initial_reservoir_pressure_moves_real_blades_only_after_duct_transit() {
    let mut short = charged(6, -0.7, 256);
    let mut long = charged(12, -0.7, 256);
    let gate = CancelGate::new_clock_free();
    let initial = 0.5 * COMPLIANCE * INITIAL_PRESSURE * INITIAL_PRESSURE;
    close(short.stored_energy_j(), initial, initial);
    let mut losses = [0.0; 2];
    let mut moved = [false; 2];
    let mut nonuniform_flow = false;
    for n in 0..256 {
        for (index, (model, delay)) in [(&mut short, 6), (&mut long, 12)].into_iter().enumerate() {
            let f = model
                .step_with_stored_source_under_gate(drive(0.0), 0.0, &gate)
                .unwrap();
            let source = f.stored_source.unwrap();
            if n < delay {
                assert_eq!(f.aperture.bore_flow_m3_s, 0.0);
                assert_eq!(f.aperture.state.opening_m, 4.0e-4);
            } else {
                moved[index] |= f.aperture.state.opening_m != 4.0e-4;
            }
            nonuniform_flow |= (source.flow_out_m3_s - f.aperture.bore_flow_m3_s).abs() > 1.0e-8;
            close(f.upstream.driven_work_j, source.outgoing_work_j, initial);
            close(f.balance_residual_j(), 0.0, initial);
            close(
                f.aperture.pressure_work_j,
                DT * f.pressure_drop_pa * f.aperture.bore_flow_m3_s,
                initial,
            );
            assert_eq!(f.source_work_j, 0.0);
            losses[index] += f.dissipated_energy_j;
            close(f.stored_energy_j + losses[index], initial, initial);
        }
    }
    assert!(moved.into_iter().all(|m| m) && nonuniform_flow);
    assert!(source_pressure(&short) != INITIAL_PRESSURE);
}

#[test]
fn g1_bore_reflection_reaches_stored_source_after_complete_return_path() {
    let (up_delay, down_delay) = (6, 8);
    let mut matched = charged(up_delay, 0.0, 256);
    let mut reflected = charged(up_delay, -0.8, 256);
    let gate = CancelGate::new_clock_free();
    let first_source_return = 2 * up_delay + 2 * down_delay;
    let mut valve_changed = false;
    let mut source_changed = false;
    for n in 0..256 {
        let a = matched
            .step_with_stored_source_under_gate(drive(0.0), 0.0, &gate)
            .unwrap();
        let b = reflected
            .step_with_stored_source_under_gate(drive(0.0), 0.0, &gate)
            .unwrap();
        if n < up_delay + 2 * down_delay {
            assert_eq!(a.aperture.state, b.aperture.state);
        } else {
            valve_changed |= a.aperture.state != b.aperture.state;
        }
        if n < first_source_return {
            assert_eq!(source_pressure(&matched), source_pressure(&reflected));
        } else {
            source_changed |= source_pressure(&matched) != source_pressure(&reflected);
        }
    }
    assert!(valve_changed && source_changed);
}

#[test]
fn g0_remote_and_local_drives_keep_separate_flow_work_in_one_balance() {
    let model = cold(6, -0.7, 384, true);
    let z = model.represented_sections(ValveSide::Upstream)[0].impedance_pa_s_m3;
    let resistance = 0.2 * z;
    let source = StoredPressureSource::new(
        SeriesImpedance::new_with_state(
            SeriesImpedanceSpec {
                resistance_pa_s_m3: resistance,
                inertance_pa_s2_m3: 0.7 * z * DT,
                compliance_m3_pa: Some(COMPLIANCE),
            },
            z,
            DT,
            ImpedanceState {
                inertive_flow_m3_s: -2.0e-6,
                compliance_pressure_pa: INITIAL_PRESSURE,
            },
        )
        .unwrap(),
    );
    let mut model = model
        .with_stored_pressure_source(1, source, &CancelGate::new_clock_free())
        .unwrap();
    let initial = model.stored_energy_j();
    let mut work = 0.0;
    let mut loss = 0.0;
    let gate = CancelGate::new_clock_free();
    let mut distinct = false;
    for n in 0..384 {
        let local = if n < 96 { 2.0 } else { 0.0 };
        let remote = if n < 192 { 10.0 } else { 0.0 };
        let f = model
            .step_with_stored_source_under_gate(drive(local), remote, &gate)
            .unwrap();
        let source = f.stored_source.unwrap();
        let expected = local * f.aperture.bore_flow_m3_s * DT + remote * source.flow_out_m3_s * DT;
        close(f.source_work_j, expected, 1.0e-5);
        close(
            source.load.dissipated_energy_j,
            resistance * source.flow_out_m3_s * source.flow_out_m3_s * DT,
            1.0e-5,
        );
        close(source.balance_residual_j(), 0.0, initial);
        distinct |= (source.flow_out_m3_s - f.aperture.bore_flow_m3_s).abs() > 1.0e-8;
        work += f.source_work_j;
        loss += f.dissipated_energy_j;
        close(f.stored_energy_j - initial + loss, work, 1.0e-5);
    }
    assert!(distinct);
}

#[test]
fn g3_zero_storage_terminal_matches_original_pressure_release_duct() {
    let mut original = cold(6, -0.7, 256, false);
    let cold_source = cold(6, -0.7, 256, true);
    let source = charge(&cold_source, None, 0.0);
    let gate = CancelGate::new_clock_free();
    let mut explicit = cold_source
        .with_stored_pressure_source(1, source, &gate)
        .unwrap();
    for n in 0..256 {
        let local = if n < 64 { 5.0 } else { 0.0 };
        let a = original.step_under_gate(drive(local), &gate).unwrap();
        let b = explicit
            .step_with_stored_source_under_gate(drive(local), 0.0, &gate)
            .unwrap();
        assert!(a.stored_source.is_none());
        close(a.aperture.bore_flow_m3_s, b.aperture.bore_flow_m3_s, 1.0e-6);
        close(
            a.aperture.state.opening_m,
            b.aperture.state.opening_m,
            4.0e-4,
        );
        close(a.stored_energy_j, b.stored_energy_j, 1.0e-7);
        assert_eq!(b.stored_source.unwrap().load.stored_energy_j, 0.0);
    }
}

#[test]
fn g4_remote_source_and_all_physical_histories_survive_refusal_and_retry() {
    let gate = CancelGate::new_clock_free();
    let mut actual = charged(6, -0.7, 32);
    let mut reference = charged(6, -0.7, 32);
    for _ in 0..16 {
        assert_eq!(
            actual
                .step_with_stored_source_under_gate(drive(0.0), 0.0, &gate)
                .unwrap(),
            reference
                .step_with_stored_source_under_gate(drive(0.0), 0.0, &gate)
                .unwrap()
        );
    }
    let state = actual
        .stored_pressure_source()
        .unwrap()
        .1
        .impedance()
        .state();
    let energy = actual.stored_energy_j();
    let node = *actual.node_frame(ValveSide::Upstream, 1).unwrap();
    let accepted = actual.valve().accepted_steps();
    // A large, finite remote drive passes the RLC solve and graph previews
    // but fails the joint energy gate. Neither network nor source may publish.
    let refused = actual
        .step_with_stored_source_under_gate(drive(0.0), 1.0e9, &gate)
        .unwrap_err();
    assert!(format!("{refused:?}").contains("complete storage"));
    assert!(actual.step_under_gate(drive(0.0), &gate).is_err());
    assert!(
        actual
            .step_with_stored_source_under_gate(drive(0.0), f64::NAN, &gate)
            .is_err()
    );
    assert!(
        actual
            .step_with_stored_source_under_gate(
                TwoSidedDrive {
                    source_pressure_pa: 0.0,
                    blade_forces_n: [f64::INFINITY, 0.0]
                },
                0.0,
                &gate
            )
            .is_err()
    );
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(
        actual
            .step_with_stored_source_under_gate(drive(0.0), 0.0, &cancelled)
            .is_err()
    );
    assert_eq!(
        actual
            .stored_pressure_source()
            .unwrap()
            .1
            .impedance()
            .state(),
        state
    );
    assert_eq!(actual.stored_energy_j(), energy);
    assert_eq!(*actual.node_frame(ValveSide::Upstream, 1).unwrap(), node);
    assert_eq!(actual.valve().accepted_steps(), accepted);
    for _ in 16..32 {
        assert_eq!(
            actual
                .step_with_stored_source_under_gate(drive(0.0), 0.0, &gate)
                .unwrap(),
            reference
                .step_with_stored_source_under_gate(drive(0.0), 0.0, &gate)
                .unwrap()
        );
    }
    let state = actual
        .stored_pressure_source()
        .unwrap()
        .1
        .impedance()
        .state();
    assert!(
        actual
            .step_with_stored_source_under_gate(drive(0.0), 0.0, &gate)
            .is_err()
    );
    assert_eq!(
        actual
            .stored_pressure_source()
            .unwrap()
            .1
            .impedance()
            .state(),
        state
    );
}

#[test]
fn g0_source_attachment_refuses_duplicate_termination_clock_and_storage() {
    let gate = CancelGate::new_clock_free();
    let passive = cold(6, -0.7, 16, false);
    let source = charge(&passive, Some(COMPLIANCE), INITIAL_PRESSURE);
    assert!(
        passive
            .with_stored_pressure_source(1, source, &gate)
            .is_err()
    );
    let model = cold(6, -0.7, 16, true);
    let wrong_clock = StoredPressureSource::new(
        SeriesImpedance::new(
            *source.impedance().spec(),
            source.impedance().port_impedance(),
            2.0 * DT,
        )
        .unwrap(),
    );
    assert!(
        model
            .with_stored_pressure_source(1, wrong_clock, &gate)
            .is_err()
    );
    let model = cold(6, -0.7, 16, true);
    let large = charge(&model, Some(1.0), 100.0);
    assert!(model.with_stored_pressure_source(1, large, &gate).is_err());
    let model = charged(6, -0.7, 16);
    assert!(model.with_stored_pressure_source(1, source, &gate).is_err());
}
