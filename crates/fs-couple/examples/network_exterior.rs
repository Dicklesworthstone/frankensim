//! Harmonic mutual radiation: two capped openings, one duct graph and one BEM.
//! Synthetic geometry; this is NOT a flute/clarinet anatomy or a rendered note.
//! cargo run --release -p fs-couple --example network_exterior > exterior.csv
//!
//! Bind different gases without changing physics code:
//! --edge-gas stem stem.gas --edge-gas left left.gas --edge-gas right right.gas
//! --exterior-gas outside.gas
//! Each file uses the shared T/P/RH, mixture or thermochemical gas grammar.
//! Thermochemical records take the same --gas-component SPECIES NASA9 TRANSPORT,
//! --gas-water and --gas-phase options as material_strings/music_render.
#[path = "../src/bin/music_render/thermochemical.rs"]
mod thermochemical;
#[path = "../src/bin/music_render/json.rs"]
mod json;
use json::json_string;
use fs_bem::panel3d::SpherePanels;
use fs_blake3::{ContentHash, DomainHasher};
use fs_couple::driving_point::network::{NetworkDrivingPoint, radiation::{
    ExteriorNetworkBudget, OpeningPatch, RadiatingNetwork,
}};
use fs_duct::{Duct, LossModel, Segment, network::{
    DuctNetwork, JunctionCorrections, NetworkEdge, NetworkMedia, NetworkNode, NodeBoundary,
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
    let media=read_media(&graph,&std::env::args().skip(1).collect::<Vec<_>>())?;
    let id=|s:&str|StableId::new(s).unwrap();
    let schema=PortKind::FluidPressureFlux.scalar_seed_schema(id("synthetic-radiating-graph"),
        CoordinateBinding::new(id("pressure-volume-flow"),id("inlet-plane"),PortOrientation::OutwardFromOwner),
        PortTimestamp::new(id("frequency-sweep"),0))
        .map_err(|e|format!("port schema refused: {e:?}"))?;
    let driver=NetworkDrivingPoint::new_coupled(graph,"inlet".into(),vec![],
        patches.iter().map(|p|p.port.clone()).collect(),schema)?;
    let mut hash=DomainHasher::new("fs-couple:network-exterior-example:v1");
    hash.update(b"synthetic equal-area cap/cylinder junction; no measured specimen equivalence");
    let provenance:ContentHash=hash.finalize();
    let gate=CancelGate::new_clock_free();
    let budget=ExteriorNetworkBudget{max_panels:128,max_field_entries:4096,max_solve_terms:2_000_000};
    let receivers=vec![[0.13,0.05,0.11]];
    let model=match media {
        Some(media)=>RadiatingNetwork::new_regional(driver,surface,patches,media,receivers,0.001,provenance,budget,&gate)?,
        None=>RadiatingNetwork::new(driver,surface,patches,gas,receivers,0.001,provenance,budget,&gate)?,
    };
    println!("frequency_hz,inlet_z_re,inlet_z_im,left_flow_re,left_flow_im,microphone_pa_re,microphone_pa_im,exterior_power_w");
    for hz in [400.0,800.0,1200.0,1600.0] {
        let row=model.solve(core::f64::consts::TAU*hz,LossModel::Bessel,C64::from_re(1e-6),&gate)?;
        let z=row.response.impedance;let q=row.outward_flows_m3_s[0];let p=row.receiver_fields[0].pressure_pa;
        println!("{hz},{},{},{},{},{},{},{}",z.re,z.im,q.re,q.im,p.re,p.im,row.exterior_surface_power_w);
    }
    Ok(())
}

fn read_media(graph:&DuctNetwork,args:&[String])->Result<Option<NetworkMedia>,String> {
    let (args,options)=thermochemical::Options::extract(args)?;
    if options.record.is_some() {
        return Err("use --edge-gas EDGE FILE and --exterior-gas FILE to bind every physical region".into());
    }
    if args.is_empty() && !options.uses_external_files() { return Ok(None); }
    let mut paths=Vec::new(); let mut exterior=None; let mut words=args.iter();
    while let Some(option)=words.next() {
        match option.as_str() {
            "--edge-gas"=>{
                let edge=words.next().ok_or("--edge-gas needs EDGE FILE")?;
                let file=words.next().ok_or("--edge-gas needs EDGE FILE")?;
                if paths.len()>=graph.edges.len() { return Err("too many network gas bindings".into()); }
                paths.push((edge.clone(),file));
            }
            "--exterior-gas"=>{
                if exterior.is_some() { return Err("--exterior-gas may only appear once".into()); }
                exterior=Some(words.next().ok_or("--exterior-gas needs FILE")?);
            }
            _=>return Err(format!("unknown network exterior option {option:?}")),
        }
    }
    let exterior=exterior.ok_or("regional media need an explicit --exterior-gas FILE")?;
    if paths.len()!=graph.edges.len() { return Err("regional media need one --edge-gas for every edge".into()); }
    let loaded=options.load()?; let bindings=loaded.bindings(); let registry=loaded.registry(&bindings);
    let mut edges=Vec::with_capacity(paths.len());
    for (edge,path) in paths {
        edges.push((edge,thermochemical::read_record(std::path::Path::new(path),registry.as_ref())?));
    }
    let outside=thermochemical::read_record(std::path::Path::new(exterior),registry.as_ref())?;
    NetworkMedia::new(graph,edges,outside).map(Some).map_err(|e|e.to_string())
}
