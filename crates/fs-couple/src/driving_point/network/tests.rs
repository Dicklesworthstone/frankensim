use super::*;
use fs_duct::{Duct, Segment, network::{JunctionCorrections, NetworkEdge, NetworkLoad, NetworkNode}};
use fs_material::gas::GasSpec;
use fs_ports::{CoordinateBinding, PortTimestamp, StableId};

pub(super) fn schema() -> PortSchema {
    let id = |s: &str| StableId::new(s).unwrap();
    PortKind::FluidPressureFlux.scalar_seed_schema(id("network-radiation-test"),
        CoordinateBinding::new(id("pressure-volume-flow"), id("inlet-plane"), PortOrientation::OutwardFromOwner),
        PortTimestamp::new(id("declared-harmonic"), 0)).unwrap()
}
pub(super) fn graph() -> DuctNetwork {
    let node = |id: &str, boundary| NetworkNode { id: id.into(), boundary };
    let port = |s: &str| NodeBoundary::Port { label: s.into() };
    let edge = |id: &str, from, to, length| NetworkEdge { id: id.into(), from, to,
        duct: Duct { segments: vec![Segment::Cylinder { radius: 0.006, length }] },
        corrections: JunctionCorrections::Disabled };
    DuctNetwork { nodes: vec![node("in", port("inlet")), node("junction", NodeBoundary::Junction),
        node("left", port("left")), node("right", port("right"))],
        edges: vec![edge("stem",0,1,0.13), edge("left",1,2,0.09), edge("right",1,3,0.09)] }
}
pub(super) fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
}
fn boundary(w: f64, diagonal: C64, mutual: C64) -> CoupledTermination {
    CoupledTermination { omega: w, labels: vec!["right".into(), "left".into()],
        z: vec![diagonal, mutual, mutual, diagonal], source_id: "synthetic full boundary".into() }
}
fn coupled() -> NetworkDrivingPoint {
    NetworkDrivingPoint::new_coupled(graph(), "inlet".into(), vec![],
        vec!["right".into(), "left".into()], schema()).unwrap()
}
#[test]
fn coupled_even_mode_changes_inlet_and_retains_the_original_owner() {
    let gas = gas(); let w = 2.0*core::f64::consts::PI*731.0;
    let z = C64::new(100_000.0, -400_000.0); let m = C64::new(50_000.0,-100_000.0);
    let supplied = boundary(w,z,m); let model = coupled();
    let actual = model.response_coupled(&gas,w,LossModel::Lossless,&supplied).unwrap();
    let direct = network_port_impedance(model.network(),&gas,w,LossModel::Lossless).unwrap()
        .drive_coupled("inlet",&[],&supplied,&gas).unwrap();
    assert_eq!(actual.impedance,direct.impedance);
    assert_eq!(actual.port_pressures,direct.port_pressures);
    assert_eq!(actual.port_inflows,direct.port_inflows);
    let independent = |value| NetworkDrivingPoint::new(graph(),"inlet".into(),
        ["left","right"].map(|s| PortTermination { label:s.into(), load:NetworkLoad::Impedance(value) }).to_vec(),schema()).unwrap();
    let even = independent(z+m).response(&gas,w,LossModel::Lossless).unwrap();
    assert!((actual.impedance-even.impedance).abs() < 1e-10*even.impedance.abs());
    let dropped = independent(z).impedance(&gas,w,LossModel::Lossless).unwrap();
    assert!((actual.impedance-dropped).abs() > 1e-4*actual.impedance.abs());
    let diagonal = model.response_coupled(&gas,w,LossModel::Lossless,&boundary(w,z,C64::ZERO)).unwrap();
    assert_eq!(diagonal.impedance,dropped);
    let flow = C64::new(1e-6,-2e-6);
    let port = model.harmonic_port_coupled(&gas,w,LossModel::Lossless,&supplied,flow).unwrap();
    assert_eq!(port.pressure,actual.impedance*flow);
    assert_eq!(port.outward_volume_flow,-flow);
}
#[test]
fn coupling_partition_is_complete_and_cannot_silently_fall_back() {
    for labels in [vec![],vec!["left"],vec!["left","left"],vec!["inlet","right"],vec!["left","absent"]] {
        assert!(NetworkDrivingPoint::new_coupled(graph(),"inlet".into(),vec![],
            labels.into_iter().map(String::from).collect(),schema()).is_err());
    }
    let model = coupled(); let gas=gas(); let w=2000.0;
    assert!(model.response(&gas,w,LossModel::Lossless).is_err());
    let mut load=boundary(w,C64::new(1e5,-2e5),C64::ZERO);
    load.labels.swap(0,1);
    assert!(model.response_coupled(&gas,w,LossModel::Lossless,&load).is_err());
    load.labels.swap(0,1); load.omega=w+1.0;
    assert!(model.response_coupled(&gas,w,LossModel::Lossless,&load).is_err());
    load.omega=w; load.z[1]=C64::new(3e5,0.0); load.z[2]=load.z[1];
    assert!(model.response_coupled(&gas,w,LossModel::Lossless,&load).is_err());
}
#[test]
fn mixed_scalar_and_coupled_boundaries_keep_original_units_and_signs() {
    let other=vec![PortTermination{ label:"right".into(),load:NetworkLoad::Closed }];
    let model=NetworkDrivingPoint::new_coupled(graph(),"inlet".into(),other.clone(),vec!["left".into()],schema()).unwrap();
    let gas=gas(); let w=2000.0; let z=C64::new(1e5,-2e5);
    let load=CoupledTermination{omega:w,labels:vec!["left".into()],z:vec![z],source_id:"synthetic one-port".into()};
    let a=model.response_coupled(&gas,w,LossModel::Bessel,&load).unwrap();
    let mut all=other; all.push(PortTermination{label:"left".into(),load:NetworkLoad::Impedance(z)});
    let b=NetworkDrivingPoint::new(graph(),"inlet".into(),all,schema()).unwrap().response(&gas,w,LossModel::Bessel).unwrap();
    assert_eq!(a.impedance,b.impedance); assert_eq!(a.port_inflows,b.port_inflows);
    assert!((a.port_pressures[1]+z*a.port_inflows[1]).abs()<1e-8*a.impedance.abs());
}


// G1: the source and load close the bore before requiring an open-circuit
// N-port impedance. The characteristic-line and tract adapters share this path.
#[test]
fn loaded_half_wave_response_does_not_require_open_port_impedance() {
    let gas = gas();
    let radius = 0.006;
    let length = 0.2;
    let omega = core::f64::consts::PI * gas.sound_speed / length;
    let zc = gas.density * gas.sound_speed / (core::f64::consts::PI * radius * radius);
    let net = DuctNetwork {
        nodes: ["input", "end"].map(|id| NetworkNode {
            id: id.into(),
            boundary: NodeBoundary::Port { label: id.into() },
        }).to_vec(),
        edges: vec![NetworkEdge {
            id: "bore".into(), from: 0, to: 1,
            duct: Duct { segments: vec![Segment::Cylinder { radius, length }] },
            corrections: JunctionCorrections::Disabled,
        }],
    };
    assert!(network_port_impedance(&net, &gas, omega, LossModel::Lossless).is_err());
    let model = NetworkDrivingPoint::new(net, "input".into(), vec![PortTermination {
        label: "end".into(), load: NetworkLoad::Impedance(C64::from_re(zc)),
    }], schema()).unwrap();
    let response = model.response(&gas, omega, LossModel::Lossless).unwrap();
    assert!((response.impedance - C64::from_re(zc)).abs() / zc < 2.0e-12);
    assert!((response.port_inflows[1] - C64::ONE).abs() < 2.0e-12);
}
