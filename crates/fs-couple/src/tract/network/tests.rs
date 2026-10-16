use super::*;
use crate::driving_point::network::NetworkDrivingPoint;
use fs_duct::{
    LossModel,
    network::{NetworkLoad, NetworkResponse, PortTermination, network_port_impedance},
};
use fs_material::gas::{GasSpec, GasState};
use std::time::Instant;

const SOURCE: &str = "synthetic diagnostic nasal branch: uniform cylinder L=0.13m, r=0.004m; no human-specimen equivalence claimed";

fn schema(kind: fs_ports::PortKind) -> fs_ports::PortSchema {
    use fs_ports::{CoordinateBinding, PortOrientation, PortTimestamp, StableId};
    let id = |s: &str| StableId::new(s).unwrap();
    kind.scalar_seed_schema(
        id("network-tract-test"),
        CoordinateBinding::new(
            id("pressure-volume-flow"),
            id("glottis-reference-plane"),
            PortOrientation::OutwardFromOwner,
        ),
        PortTimestamp::new(id("declared-frequency-sweep"), 0),
    )
    .unwrap()
}

fn branch() -> TractBranch {
    TractBranch {
        id: "nasal-synthetic".into(),
        after_section: 2,
        duct: Duct {
            segments: vec![Segment::Cylinder {
                radius: 0.004,
                length: 0.13,
            }],
        },
        corrections: JunctionCorrections::Disabled,
        port_label: "nasal".into(),
        source_id: SOURCE.into(),
    }
}
fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
}
fn loads(branched: bool) -> Vec<PortTermination> {
    let mut result = vec![PortTermination {
        label: "oral".into(),
        load: NetworkLoad::PressureRelease,
    }];
    if branched {
        result.push(PortTermination {
            label: "nasal".into(),
            load: NetworkLoad::PressureRelease,
        });
    }
    result
}

fn invariant_error(response: &NetworkResponse) -> f64 {
    let n = response.ports.len();
    let z = &response.z;
    let norm = z.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
    let mut error = 0.0_f64;
    for i in 0..n {
        for j in 0..n {
            error = error.max((z[n * i + j] - z[n * j + i]).abs() / norm);
        }
    }
    let mut h = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            h[n * i + j] = 0.5 * (z[n * i + j].re + z[n * j + i].re) / norm;
        }
    }
    let mut minors = Vec::new();
    for i in 0..n {
        minors.push(h[n * i + i]);
        for j in i + 1..n {
            minors.push(h[n * i + i] * h[n * j + j] - h[n * i + j] * h[n * j + i]);
        }
    }
    if n == 3 {
        minors.push(
            h[0] * (h[4] * h[8] - h[5] * h[7]) - h[1] * (h[3] * h[8] - h[5] * h[6])
                + h[2] * (h[3] * h[7] - h[4] * h[6]),
        );
    }
    error.max(-minors.into_iter().fold(0.0_f64, f64::min))
}

struct Sweep {
    rows: Vec<(f64, f64, f64)>,
    formants: Vec<f64>,
    worst_invariant: f64,
}
fn sweep(net: &DuctNetwork, gas: &GasState, branched: bool, step: f64) -> Sweep {
    let mut rows = Vec::new();
    let mut f = 80.0;
    let mut worst = 0.0_f64;
    while f <= 1600.0 {
        let response =
            network_port_impedance(net, gas, core::f64::consts::TAU * f, LossModel::Bessel)
                .unwrap();
        worst = worst.max(invariant_error(&response));
        let oral = response.ports.iter().position(|p| p == "oral").unwrap();
        let drive = response.drive("glottis", &loads(branched), gas).unwrap();
        rows.push((f, drive.impedance.abs(), drive.port_inflows[oral].abs()));
        f += step;
    }
    let formants = rows
        .windows(3)
        .filter(|w| w[1].1 > w[0].1 && w[1].1 > w[2].1)
        .map(|w| w[1].0)
        .collect();
    Sweep {
        rows,
        formants,
        worst_invariant: worst,
    }
}

#[test]
fn network_branched_tract_nasal() {
    // Qualitative nasal/side-cavity pole-zero behavior: Dang, Honda & Suzuki,
    // JASA 96,2088 (1994), https://pubmed.ncbi.nlm.nih.gov/7963023/ .
    // The branch geometry is synthetic; no MRI or measured-nose claim is made.
    let chart = TractChart::assaneo_a();
    let gas = gas();
    let bare = chart.network_with_branches(&[]).unwrap();
    let nasal = chart.network_with_branches(&[branch()]).unwrap();
    let anti_reference = gas.sound_speed / (2.0 * 0.13);
    let mut formant_ladder = Vec::new();
    let mut anti_ladder = Vec::new();
    for (rung, step) in [8.0, 4.0, 2.0].into_iter().enumerate() {
        let start = Instant::now();
        let a = sweep(&bare, &gas, false, step);
        let b = sweep(&nasal, &gas, true, step);
        assert!(!a.formants.is_empty() && !b.formants.is_empty());
        let anti = b
            .rows
            .iter()
            .filter(|r| r.0 > 0.85 * anti_reference && r.0 < 1.15 * anti_reference)
            .min_by(|a, b| a.2.total_cmp(&b.2))
            .unwrap();
        let baseline = a.rows.iter().find(|r| r.0 == anti.0).unwrap();
        let ratio = anti.2 / baseline.2;
        let invariant = a.worst_invariant.max(b.worst_invariant);
        let pass = b.formants[0] < a.formants[0]
            && ratio < 0.2
            && (anti.0 / anti_reference - 1.0).abs() < 0.05
            && invariant <= 1.0e-12;
        eprintln!(
            "{{\"suite\":\"fs-couple/tract\",\"case\":\"network-branched-tract-nasal\",\"inputs\":{{\"chart\":\"Assaneo2011/table2/a CC-BY\",\"branch\":\"{SOURCE}\",\"mouths\":\"explicit ideal pressure release\"}},\"axis\":\"frequency_grid_hz\",\"rung\":{rung},\"step_hz\":{step},\"observed\":{{\"bare_formants_hz\":{:?},\"branched_formants_hz\":{:?},\"antiformant_hz\":{},\"oral_flow_ratio_at_antiformant\":{ratio}}},\"oracle\":\"side-branch pole-zero mechanism; nasal geometry is synthetic\",\"err\":{invariant},\"tol\":1e-12,\"verdict\":\"{}\",\"elapsed_us\":{}}}",
            a.formants,
            b.formants,
            anti.0,
            if pass { "pass" } else { "fail" },
            start.elapsed().as_micros()
        );
        assert!(
            pass,
            "nasal fixture: first formants {:?}/{:?}; anti={} ratio={ratio}; invariant={invariant}",
            a.formants, b.formants, anti.0
        );
        formant_ladder.push(b.formants[0]);
        anti_ladder.push(anti.0);
    }
    assert!((formant_ladder[2] - formant_ladder[1]).abs() <= 4.0);
    assert!((anti_ladder[2] - anti_ladder[1]).abs() <= 4.0);
}

#[test]
fn network_player_airway_port_uses_declared_graph_without_defaults() {
    let chart = TractChart::assaneo_a();
    let net = chart.network_with_branches(&[branch()]).unwrap();
    let gas = gas();
    let w = core::f64::consts::TAU * 731.0;
    assert!(
        NetworkDrivingPoint::new(
            net.clone(),
            "glottis".into(),
            loads(false),
            schema(fs_ports::PortKind::FluidPressureFlux)
        )
        .is_err()
    );
    assert!(
        NetworkDrivingPoint::new(
            net.clone(),
            "unknown".into(),
            loads(true),
            schema(fs_ports::PortKind::FluidPressureFlux)
        )
        .is_err()
    );
    assert!(
        NetworkDrivingPoint::new(
            net.clone(),
            "glottis".into(),
            loads(true),
            schema(fs_ports::PortKind::MechanicalForceVelocity)
        )
        .is_err()
    );
    let adapter = NetworkDrivingPoint::new(
        net.clone(),
        "glottis".into(),
        loads(true),
        schema(fs_ports::PortKind::FluidPressureFlux),
    )
    .unwrap();
    let direct = network_port_impedance(&net, &gas, w, LossModel::Bessel)
        .unwrap()
        .driving_point("glottis", &loads(true), &gas)
        .unwrap();
    let consumed = adapter.impedance(&gas, w, LossModel::Bessel).unwrap();
    assert_eq!(consumed.re.to_bits(), direct.re.to_bits());
    assert_eq!(consumed.im.to_bits(), direct.im.to_bits());
    let response = adapter.response(&gas, w, LossModel::Bessel).unwrap();
    assert_eq!(response.port_inflows.len(), 3);
    assert_eq!(adapter.network().nodes.len(), net.nodes.len());
    let flow = fs_math::c64::C64::new(0.001, -0.0004);
    let port = adapter
        .harmonic_port(&gas, w, LossModel::Bessel, flow)
        .unwrap();
    assert_eq!(port.outward_volume_flow, -flow);
    assert_eq!(port.pressure, consumed * flow);
    assert_eq!(port.schema.kind(), fs_ports::PortKind::FluidPressureFlux);
    assert!(
        adapter
            .harmonic_port(
                &gas,
                w,
                LossModel::Bessel,
                fs_math::c64::C64::new(f64::NAN, 0.0)
            )
            .is_err()
    );
}

#[test]
fn network_tract_branch_admission_and_unbranched_chain_equivalence() {
    let chart = TractChart::assaneo_a();
    let gas = gas();
    let net = chart.network_with_branches(&[]).unwrap();
    let adapter = NetworkDrivingPoint::new(
        net,
        "glottis".into(),
        loads(false),
        schema(fs_ports::PortKind::FluidPressureFlux),
    )
    .unwrap();
    let w = core::f64::consts::TAU * 653.0;
    let graph = adapter.impedance(&gas, w, LossModel::Bessel).unwrap();
    let chain = fs_duct::input_impedance(
        &chart.to_duct_for_test(),
        &gas,
        w,
        LossModel::Bessel,
        fs_duct::Termination::IdealOpen,
    )
    .unwrap()
    .impedance;
    assert!((graph - chain).abs() / chain.abs() < 1.0e-10);
    let mut missing = branch();
    missing.source_id.clear();
    assert!(chart.network_with_branches(&[missing]).is_err());
    let mut beyond = branch();
    beyond.after_section = chart.sections().len() + 1;
    assert!(chart.network_with_branches(&[beyond]).is_err());
    let mut duplicate = branch();
    duplicate.port_label = "oral".into();
    assert!(chart.network_with_branches(&[duplicate]).is_err());
}
