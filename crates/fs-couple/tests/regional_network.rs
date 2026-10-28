//! Regional material gases reach the actual acoustic source, delay and observer.
use fs_bem::panel3d::SpherePanels;
use fs_blake3::ContentHash;
use fs_couple::driving_point::network::{NetworkDrivingPoint, radiation::{
    ExteriorNetworkBudget, OpeningPatch, RadiatingNetwork,
}};
use fs_duct::{Duct, LossModel, Segment, network::{
    DuctNetwork, JunctionCorrections, NetworkEdge, NetworkMedia, NetworkNode,
    NodeBoundary, network_port_impedance_regional,
}};
use fs_exec::CancelGate;
use fs_material::gas_species::input::{GasRecord, read_gas_record};
use fs_math::c64::C64;
use fs_ports::{CoordinateBinding, PortKind, PortOrientation, PortTimestamp, StableId};

fn record(text: &str) -> GasRecord {
    let mut words = text.split_ascii_whitespace();
    let result = read_gas_record(&mut words, None, None).unwrap();
    assert!(words.next().is_none());
    result
}
fn air() -> GasRecord { record("mixture 293.15 100000 mole 2 N2 0.79 O2 0.21") }
fn graph() -> DuctNetwork {
    let port = |id: &str| NetworkNode { id: id.into(), boundary: NodeBoundary::Port { label: id.into() } };
    let edge = |id: &str, from, to, length| NetworkEdge {
        id: id.into(), from, to,
        duct: Duct { segments: vec![Segment::Cylinder { radius: 0.006, length }] },
        corrections: JunctionCorrections::Disabled,
    };
    DuctNetwork {
        nodes: vec![port("inlet"), NetworkNode { id: "junction".into(), boundary: NodeBoundary::Junction },
            port("left"), port("right")],
        edges: vec![edge("stem", 0, 1, 0.13), edge("left", 1, 2, 0.09), edge("right", 1, 3, 0.09)],
    }
}
fn driver(net: DuctNetwork) -> NetworkDrivingPoint {
    let id = |s: &str| StableId::new(s).unwrap();
    let schema = PortKind::FluidPressureFlux.scalar_seed_schema(id("regional-network-test"),
        CoordinateBinding::new(id("pressure-volume-flow"), id("inlet-plane"), PortOrientation::OutwardFromOwner),
        PortTimestamp::new(id("harmonic"), 0)).unwrap();
    NetworkDrivingPoint::new_coupled(net, "inlet".into(), vec![], vec!["left".into(), "right".into()], schema).unwrap()
}
fn media(left: GasRecord, exterior: GasRecord) -> NetworkMedia {
    NetworkMedia::new(&graph(), vec![("stem".into(), air()), ("left".into(), left),
        ("right".into(), air())], exterior).unwrap()
}
fn model(media: NetworkMedia) -> RadiatingNetwork {
    let surface = SpherePanels::icosphere(0.03, 1).unwrap();
    let normal = surface.normals()[0];
    let dot = |v: &[f64;3]| normal.iter().zip(v).map(|(a,b)| a*b).sum::<f64>();
    let opposite = surface.normals().iter().enumerate()
        .min_by(|(_,a),(_,b)| dot(a).total_cmp(&dot(b))).unwrap().0;
    let patches = vec![OpeningPatch { port: "left".into(), panels: vec![0] },
        OpeningPatch { port: "right".into(), panels: vec![opposite] }];
    RadiatingNetwork::new_regional(driver(graph()), surface, patches, media,
        vec![[0.13,0.05,0.11]], 0.001, ContentHash([13;32]),
        ExteriorNetworkBudget { max_panels: 128, max_field_entries: 4096, max_solve_terms: 2_000_000 },
        &CancelGate::new_clock_free()).unwrap()
}
fn close(a: C64, b: C64) { assert!((a-b).abs() < 2e-9*b.abs().max(1e-20), "{a:?} != {b:?}"); }

#[test]
fn uniform_regional_media_preserve_the_original_coupled_response() {
    let m = media(air(), air());
    let point = driver(graph());
    let omega = core::f64::consts::TAU*721.0;
    let boundary = fs_duct::network::CoupledTermination { omega,
        labels: vec!["left".into(), "right".into()],
        z: vec![C64::new(2e5,-3e5), C64::new(2e4,-1e4), C64::new(2e4,-1e4), C64::new(3e5,-4e5)],
        source_id: "synthetic passive two-port".into() };
    for loss in [LossModel::Lossless, LossModel::Bessel] {
        let old = point.response_coupled(air().state(), omega, loss, &boundary).unwrap();
        let regional = point.response_coupled_regional(&m, omega, loss, &boundary).unwrap();
        assert_eq!(old.impedance, regional.impedance);
        assert_eq!(old.port_inflows, regional.port_inflows);
        assert_eq!(old.port_pressures, regional.port_pressures);
    }
    let mut wrong = boundary.clone(); wrong.labels.swap(0,1);
    assert!(point.response_coupled_regional(&m,omega,LossModel::Lossless,&wrong).is_err());
    assert!(point.response_regional(&m,omega,LossModel::Lossless).is_err());
}

#[test]
fn actual_branch_sound_speeds_determine_bulk_delay_and_keep_edge_identity() {
    let carbon = record("mixture 280 100000 mole 1 CO2 1");
    let gas = air(); let left_speed = carbon.state().sound_speed;
    let m = media(carbon, gas.clone());
    let expected = 2.0*(0.13/gas.state().sound_speed + 0.09/left_speed);
    let time = driver(graph()).regional_round_trip_seconds(&m).unwrap();
    assert_eq!(time, expected);
    let uniform = driver(graph()).regional_round_trip_seconds(&media(air(),air())).unwrap();
    assert!(time > uniform);
    let mut reordered = graph(); reordered.edges.reverse();
    assert_eq!(driver(reordered).regional_round_trip_seconds(&m).unwrap(),time);
}

#[test]
fn changing_only_an_interior_gas_changes_actual_sound_and_preserves_exterior_work() {
    let base = model(media(air(),air()));
    let helium = record("mixture 310 100000 mole 1 He 1");
    let changed = model(media(helium.clone(),air()));
    assert_eq!(base.boundary_hash(),changed.boundary_hash());
    assert_eq!(changed.media().unwrap().edge_record("left"),Some(&helium));
    let gate = CancelGate::new_clock_free();
    let omega = core::f64::consts::TAU*1200.0;
    let flow = C64::new(1e-6,-2e-6);
    let a = base.solve(omega,LossModel::Lossless,flow,&gate).unwrap();
    let b = changed.solve(omega,LossModel::Lossless,flow,&gate).unwrap();
    assert_eq!(a.boundary.z,b.boundary.z);
    assert!((a.response.impedance-b.response.impedance).abs() > 1e-5*a.response.impedance.abs());
    assert!((a.receiver_fields[0].pressure_pa-b.receiver_fields[0].pressure_pa).abs()
        > 1e-5*a.receiver_fields[0].pressure_pa.abs());
    let direct = network_port_impedance_regional(changed.driving_point().network(),
        changed.media().unwrap(),omega,LossModel::Lossless).unwrap()
        .drive_coupled_bound("inlet",&[],&b.boundary).unwrap();
    assert_eq!(b.response.impedance,direct.impedance);
    let input_power = 0.5*(b.inlet.pressure*flow.conj()).re;
    assert!(input_power > 0.0);
    assert!((input_power-b.exterior_surface_power_w).abs() < 2e-9*input_power);
    for row in 0..2 {
        let pressure = (0..2).fold(C64::ZERO,|sum,j|sum+b.boundary.z[2*row+j]*b.outward_flows_m3_s[j]);
        close(pressure,b.response.port_pressures[row+1]*flow);
    }
}

#[test]
fn exterior_composition_changes_reaction_and_observation_on_the_same_interior() {
    let gas = air();
    let a = model(media(gas.clone(),gas));
    let b = model(media(air(),record("mixture 293.15 100000 mole 1 CO2 1")));
    assert_ne!(a.boundary_hash(),b.boundary_hash());
    for edge in ["stem","left","right"] {
        assert_eq!(a.media().unwrap().edge_record(edge),b.media().unwrap().edge_record(edge));
    }
    let gate = CancelGate::new_clock_free(); let omega = core::f64::consts::TAU*1200.0;
    let x = a.solve(omega,LossModel::Bessel,C64::from_re(1e-6),&gate).unwrap();
    let y = b.solve(omega,LossModel::Bessel,C64::from_re(1e-6),&gate).unwrap();
    assert_ne!(x.boundary.z,y.boundary.z);
    assert_ne!(x.response.impedance,y.response.impedance);
    assert_ne!(x.receiver_fields[0].pressure_pa,y.receiver_fields[0].pressure_pa);
    let input_power = 0.5*y.inlet.pressure.re*1e-6;
    assert!(input_power > y.exterior_surface_power_w);
}
