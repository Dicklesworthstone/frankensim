//! N2 gates for the network port-flow lines: (a) one walled edge is the series duct's line and mouth flow, (b) a closed
//! side-branch edge is the series duct's closed tone hole, (c) a branched graph with two nostrils and a closed sinus is
//! passive. The dimensions are test geometry, not anatomy. The tolerance of (a) and (b), 1e-9 of the largest tap, is
//! the one-edge network test's (`network_line_matches_the_duct_line_for_one_edge`): the two sides solve the same chain
//! by different arithmetic, so equal bits are not expected.

use super::super::NetworkDrivingPoint;
use super::super::tests::{gas, schema};
use crate::driving_point::outlet::realize_mouth_flow_ir;
use crate::driving_point::{LineLadder, PASSIVITY_TOL, analyze_reflectance_ir, realize_line};
use fs_duct::network::{
    DuctNetwork, JunctionCorrections, NetworkEdge, NetworkLoad, NetworkMedia, NetworkNode, NodeBoundary,
    PortTermination, network_drive_sources_regional,
};
use fs_duct::{Duct, HoleState, LossModel, Segment, Termination};
use fs_material::gas::GasState;
use fs_material::gas_species::input::GasRecord;
use fs_phs::WallPin;
use fs_vfit::discretize::DelayedFilter;

const RATE_HZ: u32 = 48_000;
const TAP_TOL: f64 = 1.0e-9;

fn tissue() -> WallPin {
    WallPin { surface_density: 15.0, stiffness_per_area: 3.0e5, resistance: 8000.0 }
}
fn cyl(radius: f64, length: f64) -> Segment {
    Segment::Cylinder { radius, length }
}
fn node(id: &str, boundary: NodeBoundary) -> NetworkNode {
    NetworkNode { id: id.into(), boundary }
}
fn port(label: &str) -> NetworkNode {
    node(label, NodeBoundary::Port { label: label.into() })
}
fn edge(id: &str, from: usize, to: usize, segments: Vec<Segment>) -> NetworkEdge {
    NetworkEdge { id: id.into(), from, to, duct: Duct { segments }, corrections: JunctionCorrections::Disabled }
}
fn flanged(label: &str, radius: f64) -> PortTermination {
    PortTermination { label: label.into(), load: NetworkLoad::Radiation { radius, termination: Termination::FlangedOpen } }
}
fn media(net: &DuctNetwork, gas: GasState, wall: Option<WallPin>) -> NetworkMedia {
    let bound = NetworkMedia::new(
        net, net.edges.iter().map(|e| (e.id.clone(), GasRecord::Legacy(gas))).collect(), GasRecord::Legacy(gas),
    ).unwrap();
    match wall {
        Some(pin) => bound.with_edge_walls(net, net.edges.iter().map(|e| (e.id.clone(), pin)).collect()).unwrap(),
        None => bound,
    }
}
fn line_taps(line: &DelayedFilter, n: usize) -> Vec<f64> {
    let mut probe = line.clone();
    (0..n).map(|k| probe.push(if k == 0 { 1.0 } else { 0.0 }).expect("push")).collect()
}
fn relative_difference(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len());
    let scale = a.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
    a.iter().zip(b).fold(0.0_f64, |m, (x, y)| m.max((x - y).abs())) / scale
}
fn zc_of(gas: &GasState, radius: f64) -> f64 {
    gas.characteristic_impedance / (core::f64::consts::PI * radius * radius)
}

/// Series line and mouth flow of `duct` against the network's glottal line and lip flow, all on the same wall.
fn compare_with_series(case: &str, duct: &Duct, net: DuctNetwork, mouth_radius: f64, wall: WallPin) -> (Vec<f64>, Vec<f64>) {
    let gas = gas();
    let zc = zc_of(&gas, duct.segments[0].inlet_radius());
    let ladder = LineLadder::default();
    let (line, receipt) = realize_line(duct, &gas, Termination::FlangedOpen, RATE_HZ, zc, Some(&wall), &ladder).unwrap();
    let (flow, flow_receipt) =
        realize_mouth_flow_ir(duct, &gas, Termination::FlangedOpen, RATE_HZ, zc, Some(&wall), &ladder).unwrap();
    let bound = media(&net, gas, Some(wall));
    let point = NetworkDrivingPoint::new(net, "glottis".into(), vec![flanged("lips", mouth_radius)], schema()).unwrap();
    let (network_line, network_receipt) =
        point.realize_line_regional(&bound, LossModel::Bessel, RATE_HZ, zc, &ladder, None).unwrap();
    let (network_flow, network_flow_receipt) =
        point.realize_port_flow_ir_regional(&bound, LossModel::Bessel, "lips", RATE_HZ, zc, &ladder).unwrap();
    assert_eq!(network_receipt.n_fft, receipt.n_fft);
    assert_eq!(network_receipt.rungs_tried.len(), receipt.rungs_tried.len());
    assert_eq!(network_flow_receipt.n_fft, flow_receipt.n_fft);
    assert_eq!(network_flow_receipt.rungs.len(), flow_receipt.rungs.len());
    let taps = line_taps(&line, receipt.n_fft);
    let line_error = relative_difference(&taps, &line_taps(&network_line, receipt.n_fft));
    let flow_error = relative_difference(&flow, &network_flow);
    eprintln!(
        "{{\"suite\":\"fs-couple\",\"case\":\"{case}\",\"n_fft\":{},\"flow_n_fft\":{},\"line_max_relative_tap_difference\":{line_error:.3e},\"flow_max_relative_tap_difference\":{flow_error:.3e},\"tol\":{TAP_TOL:e}}}",
        receipt.n_fft, flow_receipt.n_fft
    );
    assert!(line_error <= TAP_TOL, "{case}: network line taps differ from the series line by {line_error:.3e}");
    assert!(flow_error <= TAP_TOL, "{case}: network lip flow differs from the series mouth flow by {flow_error:.3e}");
    (taps, flow)
}

/// (a) A graph with no branch is the series tract: the glottal line and the lip flow of one walled edge are the series
/// duct's on the same pin, and the pin reaches the network solve (the rigid graph differs).
#[test]
fn one_walled_edge_is_the_series_line_and_mouth_flow() {
    let segments = vec![cyl(0.008, 0.07), cyl(0.011, 0.05), cyl(0.009, 0.045)];
    let net = DuctNetwork { nodes: vec![port("glottis"), port("lips")], edges: vec![edge("tract", 0, 1, segments.clone())] };
    let duct = Duct { segments };
    let (_, walled_flow) = compare_with_series("network-port-flow-one-walled-edge", &duct, net.clone(), 0.009, tissue());
    let gas = gas();
    let point = NetworkDrivingPoint::new(net.clone(), "glottis".into(), vec![flanged("lips", 0.009)], schema()).unwrap();
    let (rigid_flow, _) = point.realize_port_flow_ir_regional(&media(&net, gas, None), LossModel::Bessel, "lips",
        RATE_HZ, zc_of(&gas, 0.008), &LineLadder::default()).unwrap();
    let n = walled_flow.len().min(rigid_flow.len());
    let change = relative_difference(&walled_flow[..n], &rigid_flow[..n]);
    assert!(change > 1.0e-6, "the edge wall did not reach the lip flow: {change:.3e}");
}

/// (b) A closed cylindrical side edge at a junction is the series duct's `ToneHole { Closed }` at the same boundary:
/// the series shunt is the closed chimney's line impedance on the same wall, the network solves the chimney as an edge.
#[test]
fn closed_side_edge_is_the_series_closed_tone_hole() {
    let (bore, hole, chimney) = (0.008, 0.004, 0.025);
    let duct = Duct { segments: vec![
        cyl(bore, 0.06),
        Segment::ToneHole { hole_radius: hole, chimney_height: chimney, bore_radius: bore, state: HoleState::Closed },
        cyl(bore, 0.08),
    ] };
    let net = DuctNetwork {
        nodes: vec![port("glottis"), node("branch", NodeBoundary::Junction), port("lips"),
            node("pouch", NodeBoundary::Termination(NetworkLoad::Closed))],
        edges: vec![edge("front", 0, 1, vec![cyl(bore, 0.06)]), edge("back", 1, 2, vec![cyl(bore, 0.08)]),
            edge("side", 1, 3, vec![cyl(hole, chimney)])],
    };
    let (taps, _) = compare_with_series("network-port-flow-closed-side-edge", &duct, net, bore, tissue());
    let gas = gas();
    let plain = Duct { segments: vec![cyl(bore, 0.06), cyl(bore, 0.08)] };
    let (line, receipt) = realize_line(&plain, &gas, Termination::FlangedOpen, RATE_HZ, zc_of(&gas, bore),
        Some(&tissue()), &LineLadder::default()).unwrap();
    let n = taps.len().min(receipt.n_fft);
    let change = relative_difference(&taps[..n], &line_taps(&line, receipt.n_fft)[..n]);
    assert!(change > 1.0e-3, "the closed branch did not change the line: {change:.3e}");
}

fn nasal_like() -> (DuctNetwork, Vec<PortTermination>) {
    let net = DuctNetwork {
        nodes: vec![port("glottis"), node("velum", NodeBoundary::Junction), port("lips"),
            node("choana", NodeBoundary::Junction), port("nostril-left"), port("nostril-right"),
            node("sinus", NodeBoundary::Termination(NetworkLoad::Closed))],
        edges: vec![
            edge("pharynx", 0, 1, vec![cyl(0.008, 0.09)]),
            edge("oral", 1, 2, vec![cyl(0.010, 0.05), cyl(0.009, 0.03)]),
            edge("port", 1, 3, vec![cyl(0.003, 0.012)]),
            edge("left", 3, 4, vec![cyl(0.004, 0.09)]),
            edge("right", 3, 5, vec![cyl(0.004, 0.09)]),
            edge("sinus", 3, 6, vec![cyl(0.0015, 0.006), cyl(0.012, 0.03)]),
        ],
    };
    (net, vec![flanged("lips", 0.009), flanged("nostril-left", 0.004), flanged("nostril-right", 0.004)])
}

/// Per unit incident wave at the glottis: the power into the graph, `(1 - |R|^2) / (2 zc)`, against the power each
/// load port radiates, `Re(p_k conj(T_k)) / 2` with `T_k` the port flow. Also checks `T_k` is the owner's outward flow.
fn power_balance(point: &NetworkDrivingPoint, bound: &NetworkMedia, loss: LossModel, zc: f64, omega: f64) -> (f64, Vec<f64>) {
    let response = network_drive_sources_regional(&point.network, bound, omega, loss, &point.reflection_source(zc).unwrap(),
        &point.terminations, None).unwrap();
    let r = point.reflected_wave(&response, zc).unwrap();
    let fields = &response.fields;
    let radiated = point.terminations.iter().map(|load| {
        let k = fields.ports.iter().position(|label| *label == load.label).unwrap();
        let flow = point.port_flow_per_incident_regional(bound, omega, loss, zc, &load.label).unwrap();
        assert_eq!(flow, fields.port_inflows[k].scale(-1.0), "port flow is not the owner's outward flow");
        0.5 * (fields.port_pressures[k] * flow.conj()).re
    }).collect();
    ((1.0 - r.norm_sq()) / (2.0 * zc), radiated)
}

/// (c) The branched graph's glottal line, walled and lossy, keeps its 4x-oversampled `|R|` under 1, every port flow
/// realizes, and at the line's own bins the ports radiate no more than enters the graph (exactly as much when the graph
/// is lossless and rigid).
#[test]
fn branched_nasal_like_graph_is_passive() {
    let gas = gas();
    let (net, loads) = nasal_like();
    let zc = zc_of(&gas, 0.008);
    // The voice realizes its glottal and airway lines one rung over the default, tolerances unchanged
    // (`GLOTTAL_LINE_MAX_TAPS`). On the default 65536 taps this walled graph's lip flow refuses with a tail energy of
    // 1.131e-10 against 1e-10 (ts2, 2026-11-08), so the gate uses the voice's ladder.
    let ladder = LineLadder { max_taps: crate::glottis::GLOTTAL_LINE_MAX_TAPS, ..LineLadder::default() };
    let point = NetworkDrivingPoint::new(net.clone(), "glottis".into(), loads.clone(), schema()).unwrap();
    let walled = media(&net, gas, Some(tissue()));
    let (line, receipt) = point.realize_line_regional(&walled, LossModel::Bessel, RATE_HZ, zc, &ladder, None).unwrap();
    let dt = 1.0 / f64::from(RATE_HZ);
    let geo_delay = point.regional_round_trip_seconds(&walled).unwrap() / dt;
    let (peak, _) = analyze_reflectance_ir(&line_taps(&line, receipt.n_fft), dt, geo_delay);
    for load in &loads {
        let (flow, flow_receipt) =
            point.realize_port_flow_ir_regional(&walled, LossModel::Bessel, &load.label, RATE_HZ, zc, &ladder).unwrap();
        assert!(flow.iter().all(|h| h.is_finite()) && flow.iter().any(|h| *h != 0.0), "{}: empty flow line", load.label);
        eprintln!("{{\"suite\":\"fs-couple\",\"case\":\"network-port-flow-nasal-like\",\"port\":\"{}\",\"n_fft\":{}}}",
            load.label, flow_receipt.n_fft);
    }
    let rigid = media(&net, gas, None);
    let mut worst_excess = f64::NEG_INFINITY;
    // One lossless bin: k, frequency, ka at the lips (radius 0.009 m in `nasal_like`), the transmitted share 1 - |R|^2, the
    // balance relative to the net input and the balance per unit incident wave.
    #[derive(Clone, Copy, Default)]
    struct LosslessBin {
        k: usize,
        hz: f64,
        ka_lips: f64,
        transmission: f64,
        relative: f64,
        incident: f64,
    }
    let (mut worst_incident, mut worst_relative) = (LosslessBin::default(), LosslessBin::default());
    let mut totals = vec![0.0; loads.len()];
    for k in (1..receipt.n_fft / 2).step_by(4) {
        let omega = core::f64::consts::TAU * k as f64 / (receipt.n_fft as f64 * dt);
        let (into, radiated) = power_balance(&point, &walled, LossModel::Bessel, zc, omega);
        for (total, &p) in totals.iter_mut().zip(&radiated) {
            assert!(p >= -PASSIVITY_TOL * into, "bin {k}: a load port supplies power {p:e}");
            *total += p;
        }
        worst_excess = worst_excess.max((radiated.iter().sum::<f64>() - into) / into);
        let (into, radiated) = power_balance(&point, &rigid, LossModel::Lossless, zc, omega);
        // Energy conservation is incident = reflected + transmitted, per unit incident wave (incident power 1 / (2 zc)). At
        // ka ~1e-4 the transmitted share 1 - |R|^2 is a difference of nearly equal numbers, so a balance relative to the
        // net input carries that difference's roundoff; the gate reads the balance on the incident scale. On the net-input
        // scale the worst bin read 7.874e-9 (ts2, 2026-11-09).
        let imbalance = radiated.iter().sum::<f64>() - into;
        let bin = LosslessBin {
            k,
            hz: omega / core::f64::consts::TAU,
            ka_lips: omega / gas.sound_speed * 0.009,
            transmission: into * 2.0 * zc,
            relative: (imbalance / into).abs(),
            incident: (imbalance * 2.0 * zc).abs(),
        };
        if bin.incident > worst_incident.incident {
            worst_incident = bin;
        }
        if bin.relative > worst_relative.relative {
            worst_relative = bin;
        }
    }
    for (measure, bin) in [("incident", worst_incident), ("relative", worst_relative)] {
        eprintln!(
            "{{\"suite\":\"fs-couple\",\"case\":\"network-port-flow-nasal-like-lossless\",\"worst_by\":\"{measure}\",\"k\":{},\"hz\":{:.4},\"ka_lips\":{:.3e},\"transmission\":{:.3e},\"relative_balance\":{:.3e},\"incident_balance\":{:.3e}}}",
            bin.k, bin.hz, bin.ka_lips, bin.transmission, bin.relative, bin.incident
        );
    }
    let worst_lossless = worst_incident.incident;
    eprintln!(
        "{{\"suite\":\"fs-couple\",\"case\":\"network-port-flow-nasal-like\",\"n_fft\":{},\"peak_reflection\":{peak:.12},\"worst_radiated_over_input_minus_1\":{worst_excess:.3e},\"worst_lossless_incident_balance\":{worst_lossless:.3e}}}",
        receipt.n_fft
    );
    assert!(peak <= 1.0 + PASSIVITY_TOL, "branched line is active on the 4x grid: {peak}");
    assert!(worst_excess <= PASSIVITY_TOL, "the ports radiate more than enters the graph: {worst_excess:.3e}");
    assert!(worst_lossless <= TAP_TOL, "the lossless rigid graph loses or gains power per incident wave: {worst_lossless:.3e}");
    assert!(totals.iter().all(|t| *t > 0.0), "a load port radiates nothing: {totals:?}");
}
