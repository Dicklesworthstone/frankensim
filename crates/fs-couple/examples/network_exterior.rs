//! Harmonic mutual radiation: two capped openings, one duct graph and one BEM.
//! Synthetic geometry; this is NOT a flute/clarinet anatomy or a rendered note.
//! cargo run --release -p fs-couple --example network_exterior > exterior.csv
use fs_bem::panel3d::SpherePanels;
use fs_blake3::{ContentHash, DomainHasher};
use fs_couple::driving_point::network::{NetworkDrivingPoint, radiation::{
    ExteriorNetworkBudget, OpeningPatch, RadiatingNetwork,
}};
use fs_duct::{Duct, LossModel, Segment, network::{
    DuctNetwork, JunctionCorrections, NetworkEdge, NetworkNode, NodeBoundary,
}};
use fs_exec::CancelGate;
use fs_material::gas::{GasSpec, GasState};
use fs_math::c64::C64;
use fs_ports::{CoordinateBinding, PortKind, PortOrientation, PortTimestamp, StableId};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let gas=GasState::try_new(&GasSpec::dry_air_ussa1976(),293.15,101325.0)?;
    // Every triangular cap represents one scalar opening. Interior/exterior
    // reference-plane correspondence is declared, not inferred from a radius.
    let surface=SpherePanels::icosphere(0.03,1)?;
    let normal=surface.normals()[0];
    let dot=|v:&[f64;3]|normal.iter().zip(v).map(|(a,b)|a*b).sum::<f64>();
    let opposite=surface.normals().iter().enumerate()
        .min_by(|(_,a),(_,b)|dot(a).total_cmp(&dot(b))).unwrap().0;
    let patches=vec![OpeningPatch{port:"left".into(),panels:vec![0]},
        OpeningPatch{port:"right".into(),panels:vec![opposite]}];
    let radius=(surface.areas()[0]/core::f64::consts::PI).sqrt();
    let port=|name:&str|NetworkNode{id:name.into(),boundary:NodeBoundary::Port{label:name.into()}};
    let edge=|name:&str,from,to,length|NetworkEdge{id:name.into(),from,to,
        duct:Duct{segments:vec![Segment::Cylinder{radius,length}]},corrections:JunctionCorrections::Disabled};
    let graph=DuctNetwork{nodes:vec![port("inlet"),NetworkNode{id:"junction".into(),boundary:NodeBoundary::Junction},
        port("left"),port("right")],edges:vec![edge("stem",0,1,0.13),edge("left",1,2,0.09),edge("right",1,3,0.09)]};
    let id=|s:&str|StableId::new(s).unwrap();
    let schema=PortKind::FluidPressureFlux.scalar_seed_schema(id("synthetic-radiating-graph"),
        CoordinateBinding::new(id("pressure-volume-flow"),id("inlet-plane"),PortOrientation::OutwardFromOwner),
        PortTimestamp::new(id("frequency-sweep"),0))?;
    let driver=NetworkDrivingPoint::new_with_coupled_ports(graph,"inlet".into(),vec![],
        patches.iter().map(|p|p.port.clone()).collect(),schema)?;
    let mut hash=DomainHasher::new("fs-couple:network-exterior-example:v1");
    hash.update(b"synthetic equal-area cap/cylinder junction; no measured specimen equivalence");
    let provenance:ContentHash=hash.finalize();
    let gate=CancelGate::new_clock_free();
    let model=RadiatingNetwork::new(driver,surface,patches,gas,vec![[0.13,0.05,0.11]],0.001,
        provenance,ExteriorNetworkBudget{max_panels:128,max_field_entries:4096,max_solve_terms:2_000_000},&gate)?;
    println!("frequency_hz,inlet_z_re,inlet_z_im,left_flow_re,left_flow_im,microphone_pa_re,microphone_pa_im,exterior_power_w");
    for hz in [400.0,800.0,1200.0,1600.0] {
        let row=model.solve(core::f64::consts::TAU*hz,LossModel::Bessel,C64::from_re(1e-6),&gate)?;
        let z=row.response.impedance;let q=row.outward_flows_m3_s[0];let p=row.receiver_fields[0].pressure_pa;
        println!("{hz},{},{},{},{},{},{},{}",z.re,z.im,q.re,q.im,p.re,p.im,row.exterior_surface_power_w);
    }
    Ok(())
}
