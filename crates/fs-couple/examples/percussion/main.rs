//! Sourced-dimension percussion reference, not a calibrated or real-time audio plugin.
//! See README.md for measured anchors/estimates and AUDIO.md for pressure export.
//! cargo run -p fs-couple --example percussion -- splash 4096 > splash.csv
//! cargo run -p fs-couple --example percussion -- drum 4096 > drum.csv
use fs_couple::render::plate::impact::{BodyPotential,ImpactBody,ImpactSystem,ImpactConfig,VolumeSpring};
use fs_couple::render::plate::impact::felt::{FeltPad,KelvinBranch};
use fs_couple::render::plate::impact::striker::{RadiusStation,StrikerProperties};
use fs_couple::modal_acoustic_time::ModalAcousticState;
use fs_material::fiber::WoolFelt;
use fs_plate::shell::profile::{ProfileStation,ProfileBudget,revolve};
use fs_plate::shell::reduction::{ShellReduction,ReductionBudget};
use fs_plate::shell::head::{TensionedDisk,TensionedDiskSpec};
use fs_plate::shell::{ShellSupport,modes_shell};
use fs_plate::{ModePair,SliceOptions};
use fs_exec::CancelGate;
use fs_dcontact::Obstacle;
use std::io::Write;

mod acoustics;

type Error=Box<dyn std::error::Error>;
fn mesh_budget()->ProfileBudget {ProfileBudget{max_nodes:10000,max_triangles:20000,max_feature_evaluations:100000}}
fn config(steps:u64,dt_s:f64)->ImpactConfig {ImpactConfig{dt_s,max_steps:steps,maximum_energy_j:20.0,
    energy_absolute_tolerance_j:1e-9,energy_relative_tolerance:1e-6,maximum_generalized_force:1e6}}
fn zero_body(potential:BodyPotential,omegas:&[f64])->ImpactBody {
    ImpactBody{potential,initial:vec![ModalAcousticState::default();omegas.len()],
        // Explicit research loss; not identified from a Zildjian or Remo sample.
        damping_per_s:omegas.iter().map(|w|0.002*w).collect()}
}
fn stick()->Result<(ImpactBody,f64),Error> {
    // Only total length/maximum shaft diameter, wood category and oval/medium
    // taper designation are published Z5A data. The stations and density below
    // are EDITABLE ESTIMATES, not a manufacturer CAD file or batch measurement.
    let profile=[(0.0,0.0068),(0.015,0.007112),(0.280,0.007112),(0.365,0.0035),
        (0.380,0.0025),(0.393,0.0045),(0.403,0.0035),(0.4064,0.0)]
        .map(|(position_m,radius_m)|RadiusStation{position_m,radius_m});
    let mass=StrikerProperties::from_profile(&profile,800.0,0.12,0.400)?;
    eprintln!("estimated_Z5A_shape: full_mass_kg={:.8}, contact_effective_mass_kg={:.8}, grip_m=0.12",
        mass.mass_kg,mass.contact_effective_mass_kg);
    // A 0.2mm approach gap and 0.8m/s downward free stroke are authored inputs.
    Ok(ImpactBody::free_mass(mass.contact_effective_mass_kg,-0.0002,0.8)?)
}
fn elastic_contact(weights:Vec<f64>)->Result<Obstacle,Error> {
    let n=weights.len();
    // Hertz research approximation: wood transverse effective modulus and tip
    // curvature are NOT measured; elastic K=4 E_eff sqrt(R_tip)/3.
    let stiffness=4.0/3.0*0.8e9*0.003_f64.sqrt();
    Ok(Obstacle::new(weights,1,n,vec![0.0],vec![1.0],stiffness,1.5,
        "estimated isotropic Hertz tip: E_eff=0.8GPa, R=3mm; not identified hickory-shell contact".into())?)
}
struct Experiment {system:ImpactSystem,force:Vec<f64>,observer_a:Vec<f64>,observer_b:Vec<f64>,pressure:Option<VolumeSpring>,acoustics:Option<acoustics::Boundary>}
fn splash(steps:u64,dt_s:f64,audio:bool)->Result<Experiment,Error> {
    // Published anchors: diameter 203.2mm; bell diameter78mm; hole diameter12.3mm;
    // edge thickness0.5mm; literature B20 E112.6GPa,nu.342,rho8607.
    // ALL interior heights and thicknesses below are explicit estimates.
    let stations=[(0.00615,0.022,0.0016),(0.012,0.0215,0.0016),(0.020,0.020,0.0015),
        (0.039,0.008,0.0011),(0.050,0.0055,0.0009),(0.070,0.003,0.00075),
        (0.085,0.0015,0.0006),(0.1016,0.0,0.0005)]
        .map(|(radius_m,height_m,thickness_m)|ProfileStation{radius_m,height_m,thickness_m});
    // Unknown proprietary hammer pattern is NOT fabricated as measured data.
    // profile::revolve accepts explicit dents/lathe relief when surveyed.
    let shell=revolve(&stations,32,112.6e9,0.342,8607.0,&[],&[],mesh_budget())?;
    let model=shell.assemble(&[],ShellSupport::Free)?;
    let pi=std::f64::consts::PI;
    let report=modes_shell(&model,((2.0*pi*50.0).powi(2),(2.0*pi*1200.0).powi(2)),&SliceOptions::default())?;
    // Keep a true vertical free-translation coordinate for felt mounting.
    // Other rigid rotations/translations are omitted in this bounded example;
    // this is NOT a fully rocking 6-DOF cymbal stand model.
    let mut phi=vec![0.0;model.free];
    for node in 0..shell.mesh.nodes.len() {if let Some(i)=model.dof_map[6*node+2] {phi[i]=1.0/shell.mass_kg.sqrt();}}
    let mut defect=vec![0.0;model.free];model.k.spmv(&phi,&mut defect);
    let residual=defect.iter().enumerate().map(|(i,r)|r*r/model.m.get(i,i)).sum::<f64>().sqrt();
    let mut modes=vec![ModePair{lambda:0.0,phi,residual,interval:(-residual,residual)}];
    modes.extend(report.modes);
    if modes.len()>32 {return Err("splash retains too many modes for this declared reference; narrow the explicit window or increase the host budget".into());}
    let reduction=ShellReduction::new(&shell.mesh,&shell.sections,&model,&modes,
        ReductionBudget{max_modes:32,max_facet_modes:20000,relative_tolerance:1e-5})?;
    let acoustics=if audio {
        let surface=reduction.radiation_surface(&shell.nodal_thickness_m,
            fs_plate::shell::reduction::radiation::RadiationSurfaceBudget{max_panels:2048,max_panel_modes:65536})?;
        Some(acoustics::Boundary::shell(&surface,1)?)
    }else{None};
    let nearest=|x:f64,y:f64|shell.mesh.tris.iter().enumerate().min_by(|(_,a),(_,b)| {
        let distance=|t:&&[usize;3]| {let cx=t.iter().map(|i|shell.mesh.nodes[*i][0]/3.0).sum::<f64>();
            let cy=t.iter().map(|i|shell.mesh.nodes[*i][1]/3.0).sum::<f64>();(cx-x).hypot(cy-y)};
        distance(a).total_cmp(&distance(b))
    }).map(|(i,_)|i).expect("nonempty admitted shell");
    let port=reduction.point_port(nearest(0.1016*2.0/3.0,0.0),[1.0/3.0;3],[0.0,0.0,-1.0])?;
    let (stick,stick_weight)=stick()?;let n=1+modes.len();let mut contact=vec![stick_weight];
    contact.extend(port.iter().map(|b|-b));
    let mut pads=Vec::new();
    // Estimated felt annulus: OD30mm/ID13mm,6mm thickness,three loaded patches
    // on each face. NOT published Zildjian dimensions or material coefficients.
    let area=pi*(0.015_f64.powi(2)-0.0065_f64.powi(2))/3.0;
    for i in 0..3 {let angle=2.0*pi*i as f64/3.0;
        let weights=reduction.point_port(nearest(0.012*angle.cos(),0.012*angle.sin()),[1.0/3.0;3],[0.0,0.0,1.0])?;
        for sign in [-1.0,1.0] {let mut b=vec![0.0];b.extend(weights.iter().map(|b|sign*b));
            pads.push(FeltPad{area_m2:area,thickness_m:0.006,precompression_m:0.0003,weights:b,
                law:WoolFelt::new(30000.0,0.2,2.2,3.0,0.15,0.7)?,prior_maximum_strain:0.05,
                creep:vec![KelvinBranch{stiffness_n_m:1500.0,viscosity_n_s_m:8.0}]});
        }
    }
    eprintln!("estimated splash reconstruction: mass_kg={},modes={},facets={},max_edge_m={}",shell.mass_kg,modes.len(),shell.mesh.tris.len(),shell.max_edge_m);
    eprintln!("modal frequencies_hz={:?}",reduction.omegas().iter().map(|w|w/(2.0*pi)).collect::<Vec<_>>());
    let omegas=reduction.omegas().to_vec();let body=zero_body(BodyPotential::Shell(reduction),&omegas);
    let system=ImpactSystem::new(vec![stick,body],vec![elastic_contact(contact)?],pads,vec![],config(steps,dt_s))?;
    let mut a=vec![0.0];a.extend(port);let mut b=vec![0.0;n];b[0]=stick_weight;
    Ok(Experiment{system,force:vec![0.0;n],observer_a:a,observer_b:b,pressure:None,acoustics})
}
fn drum(steps:u64,dt_s:f64,audio:bool)->Result<Experiment,Error> {
    // Pearl MM6 published 14x6.5in,7.5mm maple shell. Rigid cylindrical cavity
    // and clear-span radius below are geometric approximations of that shell;
    // maple elasticity, bearing-edge shape, hoops and snare wires are NOT solved.
    let radius=0.1778-0.0075;let depth=0.1651;let pi=std::f64::consts::PI;
    let mut films=Vec::new();let mut mode_sets=Vec::new();
    for (thickness,tension) in [(0.000254,3000.0),(0.0000762,1500.0)] {
        let film=TensionedDisk::new(TensionedDiskSpec{radius_m:radius,thickness_m:thickness,
            young_pa:4e9,poisson:0.38,density_kg_m3:1390.0,tension_n_m:tension,radial_intervals:5,azimuths:32},mesh_budget())?;
        let modes=fs_modal::slice_window(&film.model.k,&film.model.m,
            ((2.0*pi*80.0).powi(2),(2.0*pi*500.0).powi(2)),&SliceOptions::default())?.modes;
        if modes.is_empty() {return Err("head frequency window is empty".into());}
        mode_sets.push(modes);films.push(film);
    }
    let acoustics=if audio {Some(acoustics::Boundary::drum(&films,&mode_sets,depth,0.1778)?)}else{None};
    let (stick,stick_weight)=stick()?;let mut bodies=vec![stick];let n=1+mode_sets.iter().map(Vec::len).sum::<usize>();
    let mut contact=vec![0.0;n];contact[0]=stick_weight;let mut area=vec![0.0;n];let mut top=vec![0.0;n];let mut bottom=vec![0.0;n];let mut offset=1;
    for (head,(film,modes)) in films.iter().zip(&mode_sets).enumerate() {
        let point=film.mesh.nodes.iter().enumerate().min_by(|(_,a),(_,b)|
            (a.0-0.06).hypot(a.1).total_cmp(&(b.0-0.06).hypot(b.1))).unwrap().0;
        for (i,mode) in modes.iter().enumerate() {
            let shape=film.model.dof_map[3*point].map_or(0.0,|k|mode.phi[k]);
            // Both head coordinates are positive downward: this signed area
            // integrates COMPRESSION (negative exterior swept volume).
            area[offset+i]=(if head==0 {1.0}else{-1.0})*film.modal_area(&mode.phi)?;
            if head==0 {contact[offset+i]=-shape;top[offset+i]=shape;}else{bottom[offset+i]=shape;}
        }
        let omegas:Vec<_>=modes.iter().map(|m|m.lambda.sqrt()).collect();
        bodies.push(zero_body(BodyPotential::Linear(omegas.clone()),&omegas));offset+=modes.len();
        eprintln!("head {head}: film_mass_kg={},frequencies_hz={:?}; PET constants and tension are estimates",film.mass_kg,omegas.iter().map(|w|w/(2.0*pi)).collect::<Vec<_>>());
    }
    let volume=VolumeSpring{bulk_modulus_pa:1.2*343.0*343.0,volume_m3:pi*radius*radius*depth,areas:area};
    let system=ImpactSystem::new(bodies,vec![elastic_contact(contact)?],vec![],vec![volume.clone()],config(steps,dt_s))?;
    Ok(Experiment{system,force:vec![0.0;n],observer_a:top,observer_b:bottom,pressure:Some(volume),acoustics})
}
// The stored drum areas encode compression, so positive contraction means
// positive internal pressure. The volume-spring Hamiltonian is unchanged.
fn cavity_pressure(volume:&VolumeSpring,state:&[f64])->f64 {
    (volume.bulk_modulus_pa/volume.volume_m3)*volume.areas.iter().enumerate()
        .map(|(i,a)|a*state[2*i]).sum::<f64>()
}
fn run()->Result<(),Error> {
    let args:Vec<_>=std::env::args().skip(1).collect();
    if args.is_empty() || args.len()>3 {return Err("usage: percussion splash|drum [mechanics_steps]; or splash-wav|drum-wav [audio_frames] [full_scale_pa]; see AUDIO.md".into());}
    let audio=matches!(args[0].as_str(),"splash-wav"|"drum-wav");
    if !audio && args.len()>2 {return Err("mechanics CSV accepts only a step count".into());}
    let count=if args.len()>=2 {args[1].parse::<u64>()?}else if audio {48000}else{4096};
    let maximum=if audio {480000}else{1_000_000};
    if count==0 || count>maximum {return Err(format!("requested count must be 1..={maximum}").into());}
    let full_scale_pa=if args.len()==3 {args[2].parse::<f64>()?}else{1.0};
    if !full_scale_pa.is_finite() || full_scale_pa<=0.0 {return Err("full_scale_pa must be positive and finite".into());}
    let steps=if audio {count.checked_mul(acoustics::SUBSTEPS as u64).ok_or("sample budget overflow")?}else{count};
    let dt_s=if audio {acoustics::MECHANICAL_DT}else{2e-6};
    let mut experiment=match args[0].as_str(){
        "splash"|"splash-wav"=>splash(steps,dt_s,audio)?,
        "drum"|"drum-wav"=>drum(steps,dt_s,audio)?,
        _=>return Err("unknown experiment".into()),
    };
    let stdout=std::io::stdout();let mut out=std::io::BufWriter::new(stdout.lock());
    if audio {
        // Render and admit the complete candidate before writing a WAV header.
        let wav=acoustics::render(&mut experiment,usize::try_from(count)?,full_scale_pa)?;
        out.write_all(&wav)?;out.flush()?;return Ok(());
    }
    let gate=CancelGate::new_clock_free();
    writeln!(out,"time_s,point_a_displacement_m,point_a_velocity_m_s,point_b_displacement_m,cavity_internal_pa,total_energy_j,felt_crush_j,loss_j,balance_j")?;
    for _ in 0..steps {
        let f=experiment.system.step(&experiment.force,&gate)?;let x=experiment.system.state();
        let displacement=|weights:&[f64]|weights.iter().enumerate().map(|(i,b)|b*x[2*i]).sum::<f64>();
        let velocity=experiment.observer_a.iter().enumerate().map(|(i,b)|b*x[2*i+1]).sum::<f64>();
        let pressure=experiment.pressure.as_ref().map_or(0.0,|v|cavity_pressure(v,x));
        writeln!(out,"{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
            f.time_s,displacement(&experiment.observer_a),velocity,displacement(&experiment.observer_b),pressure,
            f.stored_energy_j,f.felt_crush_loss_j,f.dissipated_energy_j,f.balance_residual_j)?;
    }
    out.flush()?;Ok(())
}
fn main(){if let Err(e)=run(){eprintln!("percussion reference refused: {e}");std::process::exit(1);}}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn downward_batter_motion_compresses_air_and_bottom_motion_releases_it() {
        let v=VolumeSpring{bulk_modulus_pa:100.0,volume_m3:2.0,areas:vec![0.0,3.0,-3.0]};
        assert!((cavity_pressure(&v,&[0.,0.,0.2,0.,0.,0.])-30.0).abs()<1e-13);
        assert!((cavity_pressure(&v,&[0.,0.,0.,0.,0.2,0.])+30.0).abs()<1e-13);
        assert_eq!(cavity_pressure(&v,&[0.,0.,0.2,0.,0.2,0.]),0.0);
    }
}
