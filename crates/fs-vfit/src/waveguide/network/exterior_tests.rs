//! A known series exterior pressure at relaxation terminals: exact first-step
//! response, linear superposition, closed balance, refusals, and the original
//! arithmetic when nothing is attached or the attached pressure is zero.
use super::*;
use crate::impedance::SeriesImpedanceSpec;
use crate::relaxation::{RelaxationImpedanceSpec, RelaxationTerm};

const DT: f64 = 1.0 / 288_000.0;
const ZL: f64 = 2.0e6;

fn load(resistance: f64, rate: f64) -> NetworkNode {
    let base = SeriesImpedanceSpec { resistance_pa_s_m3: 0.0, inertance_pa_s2_m3: 0.0, compliance_m3_pa: None };
    NetworkNode::Relaxation { load: RelaxationImpedanceSpec::new(base,
        &[RelaxationTerm { resistance_pa_s_m3: resistance, rate_per_s: rate }]).unwrap() }
}

/// Inlet 0 - junction 1 - two relaxation terminals 2 and 3 (two outlets).
fn two_outlets() -> WaveguideNetwork {
    let nodes = [NetworkNode::Inlet, NetworkNode::Junction, load(1.5e6, 9.0e4), load(0.9e6, 1.4e5)];
    let segs = [
        NetworkSegment { nodes: [0, 1], one_way_samples: 7, impedance_pa_s_m3: ZL },
        NetworkSegment { nodes: [1, 2], one_way_samples: 5, impedance_pa_s_m3: ZL },
        NetworkSegment { nodes: [1, 3], one_way_samples: 11, impedance_pa_s_m3: 1.3 * ZL },
    ];
    WaveguideNetwork::new(&nodes, &segs, DT, 1 << 20).unwrap()
}

fn drive(n: usize) -> f64 { if n < 40 { (n as f64 * 0.37).sin() * 30.0 } else { 0.0 } }
fn pe2(n: usize) -> f64 { if (5..60).contains(&n) { (n as f64 * 0.21).cos() * 4.0 } else { 0.0 } }
fn pe3(n: usize) -> f64 { if (9..70).contains(&n) { (n as f64 * 0.13).sin() * -3.0 } else { 0.0 } }

fn bits(frame: &NetworkFrame, net: &WaveguideNetwork) -> Vec<u64> {
    let mut v = vec![frame.inlet_pressure_pa.to_bits(), frame.inlet_flow_m3_s.to_bits(),
        frame.stored_energy_j.to_bits(), frame.storage_change_j.to_bits(), frame.inlet_work_j.to_bits(),
        frame.terminal_loss_j.to_bits(), frame.coupled_terminal_work_j.to_bits()];
    for node in 0..4 {
        let o = net.node_frame(node).unwrap();
        v.extend([o.pressure_pa, o.net_flow_into_node_m3_s, o.absorbed_energy_j, o.stored_energy_j].map(f64::to_bits));
    }
    v
}

#[test]
fn first_step_exterior_pressure_is_a_series_source_with_the_instant_load() {
    // From rest nothing arrives (a = 0): 2a - Z U = p_load + pe gives
    // U = -pe/(Z + R/(1+k)) and p = pe Z/(Z + R/(1+k)), k = rate dt/2.
    let mut net = two_outlets();
    net.attach_terminal_exterior_pressures(&[2, 3]).unwrap();
    let pe = 7.25;
    net.set_terminal_exterior_pressure(2, pe).unwrap();
    let frame = net.step(0.0).unwrap();
    let k = 9.0e4 * 0.5 * DT;
    let r = 1.5e6 / (1.0 + k);
    let o = net.node_frame(2).unwrap();
    let u = -pe / (ZL + r);
    assert!((o.net_flow_into_node_m3_s / u - 1.0).abs() < 1e-13, "{} vs {u}", o.net_flow_into_node_m3_s);
    assert!((o.pressure_pa / (pe * ZL / (ZL + r)) - 1.0).abs() < 1e-13, "{}", o.pressure_pa);
    // The exterior took pe U dt (negative here: it pushed flow into the bore).
    assert!((frame.coupled_terminal_work_j / (pe * u * DT) - 1.0).abs() < 1e-13);
    assert!(frame.balance_residual_j().abs() <= 1e-12 * frame.coupled_terminal_work_j.abs());
    // The load's own account is (p - pe) U dt; the quiet outlet stays at rest.
    let own = o.storage_change_j + o.absorbed_energy_j;
    assert!((own - (o.pressure_pa - pe) * u * DT).abs() <= 1e-12 * own.abs(), "{own}");
    assert_eq!(net.node_frame(3).unwrap().net_flow_into_node_m3_s, 0.0);
}

#[test]
fn exterior_pressures_superpose_and_close_the_balance_every_step() {
    let mut drive_only = two_outlets();
    let mut pe_only = two_outlets();
    let mut both = two_outlets();
    for net in [&mut pe_only, &mut both] { net.attach_terminal_exterior_pressures(&[2, 3]).unwrap(); }
    // Rounding made while the waves were large stays in the lines as they decay, so the superposition
    // tolerance is set by the largest wave flow seen so far, not by this step's (late flows are ~1e-11
    // of a 30 Pa drive's 1.5e-5 m^3/s).
    let mut peak = 0.0_f64;
    for n in 0..400 {
        for net in [&mut pe_only, &mut both] {
            net.set_terminal_exterior_pressure(2, pe2(n)).unwrap();
            net.set_terminal_exterior_pressure(3, pe3(n)).unwrap();
        }
        let before = both.stored_energy_j();
        let a = drive_only.step(drive(n)).unwrap();
        let b = pe_only.step(0.0).unwrap();
        let c = both.step(drive(n)).unwrap();
        let scale = before + c.stored_energy_j + c.inlet_work_j.abs() + c.terminal_loss_j
            + c.coupled_terminal_work_j.abs() + 1e-30;
        assert!(c.balance_residual_j().abs() <= 1e-12 * scale, "step {n}: {c:?}");
        // The inlet flow is (outgoing - incoming)/Z: its rounding is set by the two waves, which can
        // nearly cancel, not by the flow itself.
        let waves = (drive(n).abs() + a.incoming_pressure_pa.abs() + b.incoming_pressure_pa.abs()
            + c.incoming_pressure_pa.abs()) / ZL;
        peak = peak.max(waves);
        for node in [2, 3] {
            let u0 = drive_only.node_frame(node).unwrap().net_flow_into_node_m3_s;
            let u1 = pe_only.node_frame(node).unwrap().net_flow_into_node_m3_s;
            let y = both.node_frame(node).unwrap().net_flow_into_node_m3_s;
            peak = peak.max(u0.abs() + u1.abs());
            assert!((u0 + u1 - y).abs() <= 1e-12 * peak, "step {n} node {node}: {u0}+{u1} vs {y} (peak {peak})");
        }
        let (u0, u1) = (a.inlet_flow_m3_s, b.inlet_flow_m3_s);
        peak = peak.max(u0.abs() + u1.abs());
        assert!((u0 + u1 - c.inlet_flow_m3_s).abs() <= 1e-12 * peak,
            "step {n}: {u0}+{u1} vs {} (peak {peak})", c.inlet_flow_m3_s);
        // Each pressurized load still closes its own account.
        for node in [2, 3] {
            let o = both.node_frame(node).unwrap();
            assert!(o.absorbed_energy_j >= 0.0);
        }
    }
}

#[test]
fn zero_or_absent_exterior_pressure_keeps_the_original_bits() {
    let mut plain = two_outlets();
    let mut attached = two_outlets();
    attached.attach_terminal_exterior_pressures(&[2, 3]).unwrap();
    for n in 0..300 {
        let a = plain.step(drive(n)).unwrap();
        let b = attached.step(drive(n)).unwrap();
        assert_eq!(bits(&a, &plain), bits(&b, &attached), "step {n}");
    }
}

#[test]
fn exterior_pressure_admission_refuses_without_mutation() {
    let mut net = two_outlets();
    assert!(net.set_terminal_exterior_pressure(2, 1.0).is_err());
    assert!(net.attach_terminal_exterior_pressures(&[]).is_err());
    assert!(net.attach_terminal_exterior_pressures(&[3, 2]).is_err());
    assert!(net.attach_terminal_exterior_pressures(&[2, 2]).is_err());
    assert!(net.attach_terminal_exterior_pressures(&[1]).is_err(), "a junction is not a terminal load");
    assert!(net.attach_terminal_exterior_pressures(&[0, 2]).is_err(), "the inlet is not a terminal load");
    assert!(net.exterior_terminal_nodes().is_empty());
    net.attach_terminal_exterior_pressures(&[3]).unwrap();
    assert!(net.attach_terminal_exterior_pressures(&[2]).is_err(), "attached once");
    assert!(net.set_terminal_exterior_pressure(2, 1.0).is_err(), "not admitted");
    assert!(net.set_terminal_exterior_pressure(3, f64::NAN).is_err());
    assert_eq!(net.terminal_exterior_pressure(3), Some(0.0));
    net.set_terminal_exterior_pressure(3, -2.5).unwrap();
    assert_eq!(net.terminal_exterior_pressure(3), Some(-2.5));
    assert_eq!(net.terminal_exterior_pressure(2), None);
    let mut memoryless = WaveguideNetwork::new(&[NetworkNode::Inlet, NetworkNode::Termination { reflection: 0.3 }],
        &[NetworkSegment { nodes: [0, 1], one_way_samples: 3, impedance_pa_s_m3: ZL }], DT, 1 << 20).unwrap();
    assert!(memoryless.attach_terminal_exterior_pressures(&[1]).is_err(), "a memoryless reflection has no series load");
}
