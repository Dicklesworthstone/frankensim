use super::*;
use crate::flue_loop::FlueControl;
use fs_duct::dg::{DgBoundary, DgDuct, DgMesh, PortId, MAX_DG_PROPAGATOR_DOFS};
use fs_duct::{Duct, Segment};
use fs_material::gas::{GasSpec, GasState};

fn gas(pressure: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, pressure).unwrap()
}

fn source(gas: GasState) -> FlueDgSpec {
    let duct = DgDuct::new(
        &Duct { segments: vec![Segment::Cylinder { radius: 0.001, length: 0.05 }] },
        &gas,
        DgMesh { elements: 1, order: 2 },
        DgBoundary::Port(PortId(0)),
        DgBoundary::IdealOpen,
    ).unwrap().with_viscothermal(8, [50.0, 2000.0], 0.01).unwrap();
    FlueDgSpec {
        load: DgPortLoad::new(duct, 1.0 / 96000.0, 1e-11, MAX_DG_PROPAGATOR_DOFS).unwrap(),
        minimum_jet_speed_m_s: 1.0,
        maximum_jet_mach: 0.3,
        max_jet_history_samples: 8192,
    }
}

fn geometry(source: &FlueDgSpec) -> FlueGeometry {
    FlueGeometry {
        flue_width_m: 0.010,
        flue_height_m: 0.001,
        cut_up_m: 0.004,
        labium_offset_m: 0.0,
        jet_angle_rad: 0.0,
        bore_area_m2: source.load.duct().inlet_area_m2(),
        mouth_loss_rel: 0.02,
    }
}

fn voice(gas: GasState) -> FlueVoice {
    let source = source(gas);
    let card = fs_aeroac::jetcard::mint_tonal_interim_card().unwrap();
    FlueVoice::from_dg(
        Some(&card), 8.0, 1e-6, geometry(&source), source, 96000, 3, 0,
        &CancelGate::new_clock_free(),
    ).unwrap()
}

#[test]
fn g0_original_jet_and_dg_exchange_the_same_mean_pressure_flow_and_work() {
    let mut voice = voice(gas(101325.0));
    voice.apply(FlueControl::SetBlowingPressure(20.0)).unwrap();
    let mut output = vec![0.0; 256];
    let mut integrated = 0.0;
    let mut last_step = 0;
    voice.step_block_observed(&mut output, |resonator| {
        let FlueResonator::Dg(duct) = resonator else { panic!("original DG bore") };
        let frame = duct.last_frame().unwrap();
        let port = frame.duct.port_outputs[0];
        assert_eq!(duct.load().state().step_count(), last_step + 1);
        last_step += 1;
        assert!((port.conjugate_rate - frame.inlet_pressure_pa).abs()
            <= 1e-9 * (1.0 + frame.inlet_pressure_pa.abs()));
        let work = frame.inlet_pressure_pa * frame.inlet_flow_m3_s / 96000.0;
        assert!(frame.port_work_residual_j.abs() <= 1e-10 * work.abs().max(1e-25));
        integrated += frame.duct.ledger.port_work_step_j;
        assert!(!frame.duct.ledger.balance_checked, "inferred loss stays disclosed");
        Ok(())
    }).unwrap();
    let diag = voice.diagnostics().last().unwrap();
    let scale = integrated.abs() + diag.pipe_work_j.abs();
    assert!((integrated - diag.pipe_work_j).abs() < 1e-10 * scale.max(1e-25));
    assert!(output.iter().any(|pressure| pressure.abs() > 1e-6));
    assert_eq!(last_step, 256);
    assert!(voice.dg().unwrap().stored_energy_j() > 0.0);
}

#[test]
fn g3_gas_substitution_changes_original_excitation_and_complete_dg_response() {
    let mut reference = voice(gas(101325.0));
    let mut changed = voice(gas(60000.0));
    for voice in [&mut reference, &mut changed] {
        voice.apply(FlueControl::SetBlowingPressure(20.0)).unwrap();
    }
    assert_ne!(reference.jet_speed_m_s(), changed.jet_speed_m_s());
    assert_ne!(reference.reynolds(), changed.reynolds());
    assert_ne!(reference.dg().unwrap().load().impedance_pa_s_m3(),
        changed.dg().unwrap().load().impedance_pa_s_m3());
    let mut a = vec![0.0; 256];
    let mut b = vec![0.0; 256];
    reference.step_block(&mut a).unwrap();
    changed.step_block(&mut b).unwrap();
    assert!(a.iter().zip(&b).any(|(x, y)| (x - y).abs() > 1e-8));
    assert_ne!(reference.dg().unwrap().load().state(), changed.dg().unwrap().load().state());
}

#[test]
fn g3_callback_partition_preserves_the_original_jet_and_all_bore_memory() {
    let mut whole = voice(gas(101325.0));
    let mut split = voice(gas(101325.0));
    for voice in [&mut whole, &mut split] {
        voice.apply(FlueControl::SetBlowingPressure(20.0)).unwrap();
    }
    let mut a = vec![0.0; 256];
    let mut b = vec![0.0; 256];
    whole.step_block(&mut a).unwrap();
    for chunk in b.chunks_mut(37) { split.step_block(chunk).unwrap(); }
    assert_eq!(a, b);
    assert_eq!(whole.dg().unwrap().load().state(), split.dg().unwrap().load().state());
    assert_eq!(whole.dg().unwrap().last_frame(), split.dg().unwrap().last_frame());
}

#[test]
fn g0_mismatched_clock_area_and_nonphysical_jet_controls_refuse_before_running() {
    let card = fs_aeroac::jetcard::mint_tonal_interim_card().unwrap();
    for case in 0..3 {
        let mut source = source(gas(101325.0));
        let mut geometry = geometry(&source);
        let mut rate = 96000;
        if case == 0 { rate = 48000; }
        if case == 1 { geometry.bore_area_m2 *= 1.1; }
        if case == 2 { source.maximum_jet_mach = 0.4; }
        assert!(FlueVoice::from_dg(
            Some(&card), 8.0, 1e-6, geometry, source, rate, 3, 0,
            &CancelGate::new_clock_free(),
        ).is_err());
    }
}

#[test]
fn g0_port_work_audit_keeps_the_pressure_scale_where_the_inlet_pressure_crosses_zero() {
    // music_render_wind refused here at 96 kHz: H = -1.39493 Pa and Z Q = +1.39493 Pa
    // left -2.1e-6 Pa, and the DG mean pressure matched it to 4.75e-16 Pa, the
    // rounding of H and Z Q, which is 2.3e-10 of the cancelled work.
    let spec = source(gas(101325.0));
    let inlet = *spec.load.duct().inlet_gas();
    let zc = inlet.density * inlet.sound_speed / spec.load.duct().inlet_area_m2();
    let mut dg = FlueDg { load: spec.load, characteristic_impedance_pa_s_m3: zc,
        maximum_jet_mach: spec.maximum_jet_mach, last: None };
    for _ in 0..64 { dg.load.preview_step(2e-6).unwrap().accept(); }
    let relation = dg.load.port_response().unwrap();
    let (h, z) = (relation.history_pressure_pa, relation.impedance_pa_s_m3);
    assert!(h.abs() > 1e-3, "the bore carries a history pressure to cancel");
    let flow = -h / z * (1.0 + 1e-9);
    let pressure = h + z * flow;
    assert!(pressure != 0.0 && pressure.abs() < 1e-8 * h.abs());
    dg.push(MouthStep { flow_m3_s: flow, outgoing_pa: 0.5 * (pressure + zc * flow),
        incoming_pa: 0.5 * (pressure - zc * flow), pressure_pa: pressure }).unwrap();
    let frame = dg.last_frame().unwrap();
    let scale = (h.abs() + (z * flow).abs() + pressure.abs()) * (flow / 96000.0).abs();
    assert!(frame.port_work_residual_j.abs() <= 1e-10 * scale);
    assert_eq!(dg.load.state().step_count(), 65);
}
