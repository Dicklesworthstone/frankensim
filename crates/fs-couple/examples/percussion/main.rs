//! Sourced-dimension percussion reference, not a calibrated or real-time audio plugin.
//! See README.md for measured anchors/estimates and AUDIO.md for pressure export.
//! cargo run -p fs-couple --example percussion -- splash 4096 > splash.csv
//! cargo run -p fs-couple --example percussion -- drum 4096 > drum.csv
use fs_couple::render::plate::impact::{BodyPotential,ImpactBody,ImpactSystem,ImpactConfig,VolumeSpring};
use fs_couple::render::plate::impact::felt::{FeltPad,KelvinBranch};
use fs_couple::render::plate::impact::striker::{RadiusStation,StrikerProperties};
use fs_couple::modal_acoustic_time::ModalAcousticState;
use fs_material::fiber::WoolFelt;
use fs_plate::shell::profile::ProfileBudget;
use fs_plate::shell::reduction::{ShellReduction,ReductionBudget};
use fs_plate::shell::head::{TensionedDisk,TensionedDiskSpec};
use fs_plate::shell::{ShellSupport,modes_shell};
use fs_plate::{ModePair,SliceOptions};
use fs_exec::CancelGate;
use fs_dcontact::Obstacle;
use std::io::Write;

mod acoustics;
mod cavity;
mod snare;
mod nonlinear_snare;
mod head_relaxation;
mod mallets;
mod vented;
mod mechanics;
mod playing;
mod specimen;
mod shell_prepare;
mod hihat;
#[path="hihat_sticks.rs"]
mod flexible_sticks;
mod shaft_playing;
mod drum_spec;
mod sticks;
mod muffling;
mod compliant_mute;
use playing::Stroke;
use mechanics::Mechanics;

type Error=Box<dyn std::error::Error>;
fn mesh_budget()->ProfileBudget {ProfileBudget{max_nodes:10000,max_triangles:20000,max_feature_evaluations:100000}}
fn config(steps:u64,dt_s:f64)->ImpactConfig {ImpactConfig{dt_s,max_steps:steps,maximum_energy_j:20.0,
    energy_absolute_tolerance_j:1e-9,energy_relative_tolerance:1e-6,maximum_generalized_force:1e6}}
fn zero_body(potential:BodyPotential,omegas:&[f64])->ImpactBody {
    ImpactBody{potential,initial:vec![ModalAcousticState::default();omegas.len()],
        // Explicit research loss; not identified from a Zildjian or Remo sample.
        damping_per_s:omegas.iter().map(|w|0.002*w).collect()}
}
fn stick()->Result<(ImpactBody,f64),Error> { stick_with_speed(0.8) }
fn stick_with_speed(speed_m_s:f64)->Result<(ImpactBody,f64),Error> {
    if !speed_m_s.is_finite() || !(0.0..=20.0).contains(&speed_m_s) {return Err("strike speed must be finite in 0..=20 m/s".into());}
    // Only total length/maximum shaft diameter, wood category and oval/medium
    // taper designation are published Z5A data. The stations and density below
    // are EDITABLE ESTIMATES, not a manufacturer CAD file or batch measurement.
    let profile=[(0.0,0.0068),(0.015,0.007112),(0.280,0.007112),(0.365,0.0035),
        (0.380,0.0025),(0.393,0.0045),(0.403,0.0035),(0.4064,0.0)]
        .map(|(position_m,radius_m)|RadiusStation{position_m,radius_m});
    let mass=StrikerProperties::from_profile(&profile,800.0,0.12,0.400)?;
    eprintln!("estimated_Z5A_shape: full_mass_kg={:.8}, contact_effective_mass_kg={:.8}, grip_m=0.12",
        mass.mass_kg,mass.contact_effective_mass_kg);
    // The 0.2mm gap and declared downward launch speed are mechanical inputs.
    // Changing speed changes kinetic energy, not the observer gain.
    Ok(ImpactBody::free_mass(mass.contact_effective_mass_kg,-0.0002,speed_m_s)?)
}
fn elastic_contact(weights:Vec<f64>)->Result<Obstacle,Error> {
    let n=weights.len();
    // Hertz research approximation: wood transverse effective modulus and tip
    // curvature are NOT measured; elastic K=4 E_eff sqrt(R_tip)/3.
    let stiffness=4.0/3.0*0.8e9*0.003_f64.sqrt();
    Ok(Obstacle::new(weights,1,n,vec![0.0],vec![1.0],stiffness,1.5,
        "estimated isotropic Hertz tip: E_eff=0.8GPa, R=3mm; not identified hickory-shell contact".into())?)
}
struct Experiment {flexible_sticks:[Option<fs_couple::render::plate::impact::striker::flexible::StrikerPorts>;2],mute:Option<compliant_mute::Attachment>,system:Mechanics,force:Vec<f64>,stick_weight:f64,second_stick:Option<sticks::Port>,observer_a:Vec<f64>,observer_b:Vec<f64>,pressure:Option<VolumeSpring>,acoustics:Option<acoustics::Boundary>,air:Option<cavity::InteriorPressure>}
fn splash(steps:u64,dt_s:f64,audio:bool)->Result<Experiment,Error> {
    splash_with_stroke(steps,dt_s,audio,Stroke::default())
}
fn splash_with_stroke(steps:u64,dt_s:f64,audio:bool,stroke:Stroke)->Result<Experiment,Error> {
    splash_with_specimen(steps,dt_s,audio,stroke,None)
}
fn splash_with_specimen(steps:u64,dt_s:f64,audio:bool,stroke:Stroke,supplied:Option<specimen::Specimen>)->Result<Experiment,Error> {
    splash_with_mufflers(steps,dt_s,audio,stroke,supplied,&[])
}
fn splash_with_mufflers(steps:u64,dt_s:f64,audio:bool,stroke:Stroke,supplied:Option<specimen::Specimen>,mufflers:&[muffling::Muffler])->Result<Experiment,Error> {
    splash_with_sticks(steps,dt_s,audio,stroke,supplied,mufflers,None)
}
fn splash_with_sticks(steps:u64,dt_s:f64,audio:bool,stroke:Stroke,supplied:Option<specimen::Specimen>,mufflers:&[muffling::Muffler],second:Option<Stroke>)->Result<Experiment,Error> {
    splash_with_compliant_mute(steps,dt_s,audio,stroke,supplied,mufflers,second,None)
}
#[allow(clippy::too_many_arguments)]
fn splash_with_compliant_mute(steps:u64,dt_s:f64,audio:bool,stroke:Stroke,supplied:Option<specimen::Specimen>,mufflers:&[muffling::Muffler],second:Option<Stroke>,mute:Option<&compliant_mute::Spec>)->Result<Experiment,Error> {
    splash_with_mallets(steps,dt_s,audio,stroke,supplied,mufflers,second,mute,&mallets::Selection::default())
}
#[allow(clippy::too_many_arguments)]
fn splash_with_mallets(steps:u64,dt_s:f64,audio:bool,stroke:Stroke,supplied:Option<specimen::Specimen>,mufflers:&[muffling::Muffler],second:Option<Stroke>,mute:Option<&compliant_mute::Spec>,mallets:&mallets::Selection)->Result<Experiment,Error> {
    splash_with_shafts(steps,dt_s,audio,stroke,supplied,mufflers,second,mute,mallets,&shaft_playing::Selection::default())
}
#[allow(clippy::too_many_arguments)]
fn splash_with_shafts(steps:u64,dt_s:f64,audio:bool,stroke:Stroke,supplied:Option<specimen::Specimen>,mufflers:&[muffling::Muffler],second:Option<Stroke>,mute:Option<&compliant_mute::Spec>,mallets:&mallets::Selection,shafts:&shaft_playing::Selection)->Result<Experiment,Error> {
    shafts.admit("splash",second,mallets)?;
    mallets.admit("splash",stroke,second)?;
    muffling::admit_command(mufflers,"splash")?;
    if let Some(spec)=mute {spec.admit_command("splash")?;}
    let imported=supplied.is_some();let specimen=supplied.unwrap_or_else(specimen::Specimen::reference);
    let (shell,reduction)=shell_prepare::prepare(&specimen,dt_s)?;
    let pi=std::f64::consts::PI;let [lower,upper]=specimen.band_hz;
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
    let (triangle,barycentric)=match stroke.position_m {
        Some(p)=>playing::shell_location(&shell.mesh.nodes,&shell.mesh.tris,p)?,
        None if imported=>{
            playing::shell_location(&shell.mesh.nodes,&shell.mesh.tris,specimen.default_strike_position()?)?
        }
        None=>(nearest(0.1016*2.0/3.0,0.0),[1.0/3.0;3]),
    };
    let port=reduction.point_port(triangle,barycentric,[0.0,0.0,-1.0])?;
    // Preserve the original stick/shell prefix. A second striker follows the
    // SAME shell modes, before the time owner's private Kelvin coordinates.
    let second_coordinate=1+reduction.mode_count();let structural=second_coordinate+usize::from(second.is_some());
    let attachment=mute.map(|spec|spec.shell(&reduction,&shell.mesh.nodes,&shell.mesh.tris,structural)).transpose()?;
    let base=structural+attachment.as_ref().map_or(0,|a|a.bodies.len());
    let mut shaft=shafts.build_with_mallets(base,second_coordinate,stroke,second,dt_s,mallets)?;
    let n=shaft.total;
    let mut mallet_pads=Vec::new();
    let (stick,stick_weight)=match &mallets.first {
        Some(spec)=>{let tip=spec.compile_shell(&reduction,&shell,stroke,0,1,n)?;
                let tip=spec.bind_shaft(tip,0,&mut shaft)?;
            mallet_pads.extend(tip.pads);(tip.body,tip.port.inverse_sqrt_mass)},
        None=>match shaft.bodies[0].take() {
            Some(body)=>body,None=>stick_with_speed(stroke.speed_m_s)?,
        },
    };
    let second=second.map(|stroke|->Result<_,Error>{
        match &mallets.second {
            Some(spec)=>{let tip=spec.compile_shell(&reduction,&shell,stroke,second_coordinate,1,n)?;
                let tip=spec.bind_shaft(tip,1,&mut shaft)?;
                mallet_pads.extend(tip.pads);
                Ok((tip.body,None,sticks::Port{coordinate:second_coordinate,weight:tip.port.inverse_sqrt_mass}))},
            None=>{
                if let Some((body,weight))=shaft.bodies[1].take() {
                    let p=stroke.position_m.ok_or("second shaft needs a physical XY station")?;
                    let (face,bary)=playing::shell_location(&shell.mesh.nodes,&shell.mesh.tris,p)?;
                    let shapes=reduction.point_port(face,bary,[0.0,0.0,-1.0])?;
                    if shapes.iter().all(|b|*b==0.0){return Err("second shaft has no retained shell participation".into());}
                    let mut row=shaft.tip_row(1,weight)?;
                    for (i,b) in shapes.iter().enumerate(){row[1+i]=-b;}
                    Ok((body,Some(elastic_contact(row)?),sticks::Port{coordinate:second_coordinate,weight}))
                }else{let (body,contact,port)=sticks::build_shell(stroke,&reduction,&shell.mesh.nodes,
                    &shell.mesh.tris,second_coordinate,n)?;Ok((body,Some(contact),port))}
            },
        }
    }).transpose()?;
    let mut contact=shaft.tip_row(0,stick_weight)?;
    for (i,b) in port.iter().enumerate(){contact[1+i]=-b;}
    let mut pads=Vec::new();
    // Estimated felt annulus: OD30mm/ID13mm,6mm thickness,three loaded patches
    // on each face. NOT published Zildjian dimensions or material coefficients.
    let area=pi*(0.015_f64.powi(2)-0.0065_f64.powi(2))/3.0;
    for i in 0..3 {let angle=2.0*pi*i as f64/3.0;
        let p=[0.012*angle.cos(),0.012*angle.sin()];
        // Supplied geometry cannot snap a stand pad across a mounting hole.
        // Hardware dimensions/materials remain explicit estimates of this host.
        let (face,bary)=if imported {playing::shell_location(&shell.mesh.nodes,&shell.mesh.tris,p)?}
            else {(nearest(p[0],p[1]),[1.0/3.0;3])};
        let weights=reduction.point_port(face,bary,[0.0,0.0,1.0])?;
        for sign in [-1.0,1.0] {let mut b=vec![0.0];b.extend(weights.iter().map(|b|sign*b));
            b.resize(n,0.0); // Neither stick directly compresses the stand felt.
            pads.push(FeltPad{area_m2:area,thickness_m:0.006,precompression_m:0.0003,weights:b,
                law:WoolFelt::new(30000.0,0.2,2.2,3.0,0.15,0.7)?,prior_maximum_strain:0.05,
                creep:vec![KelvinBranch{stiffness_n_m:1500.0,viscosity_n_s_m:8.0}]});
        }
    }
    eprintln!("shell input={}, mass_kg={},modes={},facets={},max_edge_m={},band_hz={lower}..{upper}; stand felt, stick and contact remain estimated; no calibration or full-band claim",
        if imported {specimen.input_label()}else{"estimated splash"},shell.mass_kg,reduction.mode_count(),shell.mesh.tris.len(),shell.max_edge_m);
    eprintln!("modal frequencies_hz={:?}",reduction.omegas().iter().map(|w|w/(2.0*pi)).collect::<Vec<_>>());
    let mut dampers=muffling::shell_ports(mufflers,&reduction,&shell.mesh.nodes,&shell.mesh.tris)?;
    for damper in &mut dampers {damper.weights.resize(n,0.0);}
    let omegas=reduction.omegas().to_vec();let body=zero_body(BodyPotential::Shell(reduction),&omegas);
    let mut bodies=vec![stick,body];
    // A selected felt face REPLACES its Hertz contact. Preserve all six stand
    // histories first; each following mallet site has its own material history.
    let mut contacts=if mallets.first.is_some(){Vec::new()}else{vec![elastic_contact(contact)?]};
    let second_stick=second.map(|(body,contact,port)|{bodies.push(body);if let Some(c)=contact{contacts.push(c);}port});
    pads.extend(mallet_pads);
    let mute=attachment.map(|a| {
        let observation=mute.expect("compiled mute specification").observation(a.ports,pads.len());
        bodies.extend(a.bodies);pads.extend(a.pads);observation
    });
    // Stand/mute/other-hand mallet rows do not act directly on shaft flexure.
    for pad in &mut pads {pad.weights.resize(n,0.0);}
    bodies.append(&mut shaft.elastic);
    let system=ImpactSystem::new_with_dampers(bodies,contacts,pads,vec![],dampers,config(steps,dt_s))?;
    let mut a=vec![0.0];a.extend(port);a.resize(n,0.0);let b=shaft.tip_row(0,stick_weight)?;
    Ok(Experiment{flexible_sticks:shaft.ports,mute,system:Mechanics::Reference(system),force:vec![0.0;n],stick_weight,second_stick,observer_a:a,observer_b:b,pressure:None,acoustics,air:None})
}
fn drum(steps:u64,dt_s:f64,audio:bool,prepared:bool)->Result<Experiment,Error> {
    drum_with_wires(steps,dt_s,audio,prepared,None)
}
fn drum_with_wires(steps:u64,dt_s:f64,audio:bool,prepared:bool,snares:Option<snare::SnareSet>)->Result<Experiment,Error> {
    drum_with_playing(steps,dt_s,audio,prepared,snares,false,Stroke::default())
}
fn drum_with_playing(steps:u64,dt_s:f64,audio:bool,prepared:bool,snares:Option<snare::SnareSet>,stretching:bool,stroke:Stroke)->Result<Experiment,Error> {
    drum_with_air(steps,dt_s,audio,prepared,snares,stretching,stroke,false,None)
}
#[allow(clippy::too_many_arguments)]
fn drum_with_air(steps:u64,dt_s:f64,audio:bool,prepared:bool,snares:Option<snare::SnareSet>,stretching:bool,stroke:Stroke,distributed_cavity:bool,neck:Option<cavity::NeckOptions>)->Result<Experiment,Error> {
    drum_with_spec(steps,dt_s,audio,prepared,snares,stretching,stroke,distributed_cavity,neck,None)
}
#[allow(clippy::too_many_arguments)]
fn drum_with_spec(steps:u64,dt_s:f64,audio:bool,prepared:bool,snares:Option<snare::SnareSet>,stretching:bool,stroke:Stroke,distributed_cavity:bool,neck:Option<cavity::NeckOptions>,supplied:Option<drum_spec::Spec>)->Result<Experiment,Error> {
    drum_with_sticks(steps,dt_s,audio,prepared,snares,stretching,stroke,distributed_cavity,neck,supplied,None)
}
#[allow(clippy::too_many_arguments)]
fn drum_with_sticks(steps:u64,dt_s:f64,audio:bool,prepared:bool,snares:Option<snare::SnareSet>,stretching:bool,stroke:Stroke,distributed_cavity:bool,neck:Option<cavity::NeckOptions>,supplied:Option<drum_spec::Spec>,second:Option<Stroke>)->Result<Experiment,Error> {
    drum_with_mufflers(steps,dt_s,audio,prepared,snares,stretching,stroke,distributed_cavity,neck,supplied,second,&[])
}
#[allow(clippy::too_many_arguments)]
fn drum_with_mufflers(steps:u64,dt_s:f64,audio:bool,prepared:bool,snares:Option<snare::SnareSet>,stretching:bool,stroke:Stroke,distributed_cavity:bool,neck:Option<cavity::NeckOptions>,supplied:Option<drum_spec::Spec>,second:Option<Stroke>,mufflers:&[muffling::Muffler])->Result<Experiment,Error> {
    drum_with_cavity_loss(steps,dt_s,audio,prepared,snares,stretching,stroke,distributed_cavity,
        neck,supplied,second,mufflers,0.0)
}
#[allow(clippy::too_many_arguments)]
fn drum_with_cavity_loss(steps:u64,dt_s:f64,audio:bool,prepared:bool,snares:Option<snare::SnareSet>,stretching:bool,stroke:Stroke,distributed_cavity:bool,neck:Option<cavity::NeckOptions>,supplied:Option<drum_spec::Spec>,second:Option<Stroke>,mufflers:&[muffling::Muffler],drag_per_s:f64)->Result<Experiment,Error> {
    drum_with_compliant_mute(steps,dt_s,audio,prepared,snares,stretching,stroke,distributed_cavity,
        neck,supplied,second,mufflers,drag_per_s,None)
}
#[allow(clippy::too_many_arguments)]
fn drum_with_compliant_mute(steps:u64,dt_s:f64,audio:bool,prepared:bool,snares:Option<snare::SnareSet>,stretching:bool,stroke:Stroke,distributed_cavity:bool,neck:Option<cavity::NeckOptions>,supplied:Option<drum_spec::Spec>,second:Option<Stroke>,mufflers:&[muffling::Muffler],drag_per_s:f64,mute:Option<&compliant_mute::Spec>)->Result<Experiment,Error> {
    drum_with_radiation(steps,dt_s,audio,prepared,snares,stretching,stroke,distributed_cavity,
        neck,supplied,second,mufflers,drag_per_s,mute,false)
}
#[allow(clippy::too_many_arguments)]
fn drum_with_radiation(steps:u64,dt_s:f64,audio:bool,prepared:bool,snares:Option<snare::SnareSet>,stretching:bool,stroke:Stroke,distributed_cavity:bool,neck:Option<cavity::NeckOptions>,supplied:Option<drum_spec::Spec>,second:Option<Stroke>,mufflers:&[muffling::Muffler],drag_per_s:f64,mute:Option<&compliant_mute::Spec>,prescribed_vent:bool)->Result<Experiment,Error> {
    drum_with_material(steps,dt_s,audio,prepared,snares,stretching,stroke,distributed_cavity,neck,
        supplied,second,mufflers,drag_per_s,mute,prescribed_vent,None)
}
#[allow(clippy::too_many_arguments)]
fn drum_with_material(steps:u64,dt_s:f64,audio:bool,prepared:bool,snares:Option<snare::SnareSet>,stretching:bool,stroke:Stroke,distributed_cavity:bool,neck:Option<cavity::NeckOptions>,supplied:Option<drum_spec::Spec>,second:Option<Stroke>,mufflers:&[muffling::Muffler],drag_per_s:f64,mute:Option<&compliant_mute::Spec>,prescribed_vent:bool,relaxation:Option<&head_relaxation::Spec>)->Result<Experiment,Error> {
    drum_with_mallets(steps,dt_s,audio,prepared,snares,stretching,stroke,distributed_cavity,neck,
        supplied,second,mufflers,drag_per_s,mute,prescribed_vent,relaxation,&mallets::Selection::default())
}
#[allow(clippy::too_many_arguments)]
fn drum_with_mallets(steps:u64,dt_s:f64,audio:bool,prepared:bool,snares:Option<snare::SnareSet>,stretching:bool,stroke:Stroke,distributed_cavity:bool,neck:Option<cavity::NeckOptions>,supplied:Option<drum_spec::Spec>,second:Option<Stroke>,mufflers:&[muffling::Muffler],drag_per_s:f64,mute:Option<&compliant_mute::Spec>,prescribed_vent:bool,relaxation:Option<&head_relaxation::Spec>,mallets:&mallets::Selection)->Result<Experiment,Error> {
    drum_with_shafts(steps,dt_s,audio,prepared,snares,stretching,stroke,distributed_cavity,neck,supplied,
        second,mufflers,drag_per_s,mute,prescribed_vent,relaxation,mallets,&shaft_playing::Selection::default())
}
#[allow(clippy::too_many_arguments)]
fn drum_with_shafts(steps:u64,dt_s:f64,audio:bool,prepared:bool,snares:Option<snare::SnareSet>,stretching:bool,stroke:Stroke,distributed_cavity:bool,neck:Option<cavity::NeckOptions>,supplied:Option<drum_spec::Spec>,second:Option<Stroke>,mufflers:&[muffling::Muffler],drag_per_s:f64,mute:Option<&compliant_mute::Spec>,prescribed_vent:bool,relaxation:Option<&head_relaxation::Spec>,mallets:&mallets::Selection,shafts:&shaft_playing::Selection)->Result<Experiment,Error> {
    shafts.admit(if snares.is_some(){"snare"}else{"drum"},second,mallets)?;
    if prepared && mallets.enabled() {return Err("felt mallets require the coupled felt/history owner".into());}
    mallets.admit(if snares.is_some() {"snare"}else{"drum"},stroke,second)?;
    if prescribed_vent && (!audio || !distributed_cavity || neck.is_none()) {
        return Err("prescribed vent radiation requires an audio observer and a distributed-cavity neck".into());
    }
    if let Some(spec)=mute {
        spec.admit_command("drum")?;
        if prepared || snares.is_some() {return Err("compliant drum pad needs the nonlinear contact/felt owner".into());}
    }
    muffling::admit_command(mufflers,"drum")?;
    cavity::validate_drag(drag_per_s)?;
    if drag_per_s!=0.0 && !distributed_cavity {return Err("acoustic drag requires distributed cavity inertia".into());}
    if neck.is_some() && (!distributed_cavity || audio && !prescribed_vent) {return Err("vented audio requires explicit --prescribed-vent-radiation; fully coupled radiation loading is not implemented".into());}
    nonlinear_snare::admit_image(prepared,snares.is_some(),
        stretching || relaxation.is_some() || mallets.enabled() || snares.is_some_and(|s|s.stretching.is_some() || s.carrier.is_some()))?;
    let extra_modes=match snares {Some(spec)=>spec.mode_count()?,None=>0};
    // One declaration supplies BOTH head pencils, the air volume and the
    // closed exterior. No independently retuned oscillator or stock drum mesh.
    let imported=supplied.is_some();let spec=supplied.unwrap_or_else(drum_spec::Spec::reference);
    spec.admit_clock(dt_s,audio)?;
    if let Some(material)=relaxation {material.admit(&spec)?;}
    if let Some(wires)=snares {spec.admit_snare(wires)?;}
    let radius=spec.radius_m;let depth=spec.depth_m;let pi=std::f64::consts::PI;
    let (films,mode_sets)=spec.prepare(dt_s,audio)?;
    let acoustics=if audio {Some(acoustics::Boundary::drum(&films,&mode_sets,depth,spec.outer_radius_m)?)}else{None};
    // Keep the original nearest-node strike for the no-file reference. A new
    // geometry instead gets an interior relative location unless explicitly set.
    let position=stroke.position_m.or_else(||imported.then_some([0.35*radius,0.0]));
    eprintln!("drum input={}; clear_radius_m={radius}, depth_m={depth}, outer_radius_m={}, head_window_hz={:?}, radial_intervals={}, azimuths={}, strike_xy_m={position:?}; rigid shell/rim, unchanged declared gas and estimated stick/contact; no specimen certification",
        if imported {"supplied SI specification"}else{"estimated 14x6.5in reference"},
        spec.outer_radius_m,spec.band_hz,spec.radial_intervals,spec.azimuths);
    let structural=1+mode_sets.iter().map(Vec::len).sum::<usize>()+extra_modes+usize::from(second.is_some());
    let attachment=mute.map(|spec|spec.head(&films,&mode_sets,structural)).transpose()?;
    let base=structural+attachment.as_ref().map_or(0,|a|a.bodies.len());
    let second_coordinate=1+mode_sets.iter().map(Vec::len).sum::<usize>();
    let mut shaft=shafts.build_with_mallets(base,second_coordinate,stroke,second,dt_s,mallets)?;
    let n=shaft.total;
    let mut mallet_pads=Vec::new();
    let (stick,stick_weight)=match &mallets.first {
        Some(spec)=>{let tip=spec.compile(&films[0],&mode_sets[0],stroke,0,n)?;
                let tip=spec.bind_shaft(tip,0,&mut shaft)?;
            mallet_pads.extend(tip.pads);(tip.body,tip.port.inverse_sqrt_mass)},
        None=>match shaft.bodies[0].take() {
            Some(body)=>body,None=>stick_with_speed(stroke.speed_m_s)?,
        },
    };
    let mut bodies=vec![stick];
    let mut contact=shaft.tip_row(0,stick_weight)?;let mut area=vec![0.0;n];let mut top=vec![0.0;n];let mut bottom=vec![0.0;n];let mut offset=1;
    for (head,(film,modes)) in films.iter().zip(&mode_sets).enumerate() {
        let point=film.mesh.nodes.iter().enumerate().min_by(|(_,a),(_,b)|
            (a.0-0.06).hypot(a.1).total_cmp(&(b.0-0.06).hypot(b.1))).unwrap().0;
        let explicit_shapes=match position {
            Some(position)=>Some(fs_couple::render::plate::impact::linear::wire::film_shapes(film,modes,&[position])?.remove(0)),
            None=>None,
        };
        for (i,mode) in modes.iter().enumerate() {
            let shape=explicit_shapes.as_ref().map_or_else(
                ||film.model.dof_map[3*point].map_or(0.0,|k|mode.phi[k]),|row|row[i]);
            // Both head coordinates are positive downward: this signed area
            // integrates COMPRESSION (negative exterior swept volume).
            area[offset+i]=(if head==0 {1.0}else{-1.0})*film.modal_area(&mode.phi)?;
            if head==0 {contact[offset+i]=-shape;top[offset+i]=shape;}else{bottom[offset+i]=shape;}
        }
        let omegas:Vec<_>=modes.iter().map(|m|m.lambda.sqrt()).collect();
        let potential=if stretching {
            let rim:Vec<_>=(0..film.mesh.nodes.len()).filter(|&i|film.model.dof_map[3*i].is_none()).collect();
            let law=fs_couple::render::plate::impact::membrane::MembranePotential::from_pencil(
                &film.mesh,&film.section,&film.model,modes,&rim,
                fs_plate::shell::head::nonlinear::MembraneReductionBudget {
                    max_modes:31,max_nodes:512,max_facet_pairs:1_000_000,
                    max_solve_entries:2_000_000,relative_tolerance:1e-5,
                },0.2)?;
            eprintln!("head {head}: geometric stretching enabled, in_plane_residual={}, maximum_accepted_slope={}; fixed in-plane rim, relaxed interior",
                law.reduction().solve_residual(),law.slope_limit());
            BodyPotential::Membrane(law)
        }else{BodyPotential::Linear(omegas.clone())};
        let mut body=zero_body(potential,&omegas);
        body.damping_per_s=omegas.iter().map(|w|2.0*spec.heads[head].damping_ratio*w).collect();
        bodies.push(body);offset+=modes.len();
        eprintln!("head {head}: film_mass_kg={},frequencies_hz={:?}, tension_n_m={}, damping_ratio={}; material authority belongs to the input, not a brand name",
            film.mass_kg,omegas.iter().map(|w|w/(2.0*pi)).collect::<Vec<_>>(),
            spec.heads[head].tension_n_m,spec.heads[head].damping_ratio);
    }
    // The selected felt face REPLACES its hard-tip contact, never adds to it.
    let mut contacts=if mallets.first.is_some() {Vec::new()}else{vec![elastic_contact(contact)?]};
    // Keep BOTH original head ranges fixed. The second striker follows them,
    // before any wires/air, and has zero direct volume/radiation participation.
    let head_end=offset;
    let second_stick=if let Some(stroke)=second {
        let port=match &mallets.second {
            Some(spec)=>{let tip=spec.compile(&films[0],&mode_sets[0],stroke,offset,n)?;
                let tip=spec.bind_shaft(tip,1,&mut shaft)?;
                let port=sticks::Port {coordinate:offset,weight:tip.port.inverse_sqrt_mass};
                bodies.push(tip.body);mallet_pads.extend(tip.pads);port},
            None=>{
                if let Some((body,weight))=shaft.bodies[1].take() {
                    let position=stroke.position_m.ok_or("second shaft needs a physical XY station")?;
                    let shapes=fs_couple::render::plate::impact::linear::wire::film_shapes(
                        &films[0],&mode_sets[0],&[position])?.remove(0);
                    if shapes.iter().all(|b|*b==0.0){return Err("second shaft has no retained head participation; rim strikes are not implemented".into());}
                    let mut row=shaft.tip_row(1,weight)?;
                    for (i,b) in shapes.iter().enumerate(){row[1+i]=-b;}
                    bodies.push(body);contacts.push(elastic_contact(row)?);
                    sticks::Port{coordinate:offset,weight}
                }else{let (body,contact,port)=sticks::build(stroke,&films[0],&mode_sets[0],offset,n)?;
                    bodies.push(body);contacts.push(contact);port}
            },
        };
        offset+=1;Some(port)
    }else{None};
    if let Some(spec)=snares {
        let bottom_start=1+mode_sets[0].len();
        let (wire_bodies,wire_contacts)=spec.assemble(&films[1],&mode_sets[1],
            bottom_start..head_end,offset,n)?;
        bodies.extend(wire_bodies);contacts.extend(wire_contacts);
    }
    let mut pads=mallet_pads;
    let mute=attachment.map(|a| {
        let observation=mute.expect("compiled mute specification").observation(a.ports,pads.len());
        bodies.extend(a.bodies);pads.extend(a.pads);observation
    });
    for pad in &mut pads {pad.weights.resize(n,0.0);}
    bodies.append(&mut shaft.elastic);
    let volume=VolumeSpring{bulk_modulus_pa:1.2*343.0*343.0,volume_m3:spec.volume_m3(),areas:area};
    // Both images consume the identical geometric reduction, strike port,
    // constitutive contact, loss coefficients and air volume. Only the discrete
    // realization changes. Neither path modifies the microphone/BEM boundary.
    let dampers=muffling::head_ports(mufflers,&films,&mode_sets,n)?;
    let mut air=None;
    let system=if distributed_cavity {
        if prepared {
            let mut configuration=mechanics::coupled_config(steps,dt_s,snares.is_some())?;
            if drag_per_s>0.0 || neck.is_some() {
                // Explicit work envelope for this new acoustic-loss image:
                // up to 8 pressure springs + 7 momentum drags + 1 vent drag
                // + 16 solid mufflers. Modal/contact/energy limits are unchanged.
                configuration.coupling.max_connections=32;
                configuration.coupling.max_setup_terms=1_000_000;
            }
            let (system,probe)=cavity::build_prepared_with_losses(&films,&mode_sets,bodies,contacts,
                dampers,radius,depth,configuration,neck,drag_per_s)?;
            air=Some(probe);Mechanics::Prepared(system)
        } else {
            let (system,probe)=cavity::build_with_pads(&films,&mode_sets,bodies,contacts,pads,
                dampers,radius,depth,steps,dt_s,neck,drag_per_s)?;
            air=Some(probe);Mechanics::Reference(system)
        }
    }else if prepared {
        Mechanics::prepared_with_dampers(bodies,contacts,volume.clone(),pi*radius*radius,dampers,
            mechanics::coupled_config(steps,dt_s,snares.is_some())?)?
    }else{
        Mechanics::Reference(ImpactSystem::new_with_dampers(bodies,contacts,pads,vec![volume.clone()],dampers,config(steps,dt_s))?)
    };
    let system=match (system,relaxation) {
        (Mechanics::Reference(system),Some(material))=>Mechanics::Reference(material.attach(system,&films,&mode_sets,dt_s)?),
        (system,None)=>system,
        (_,Some(_))=>return Err("material memory cannot be discarded by the linear-only modal image".into()),
    };
    // Acoustic inertia never receives an external strike or a direct solid-radiation projection.
    let count=air.as_ref().map_or(n,|a|a.coupling.total_modes());
    top.resize(count,0.0);bottom.resize(count,0.0);
    let acoustics=if prescribed_vent {
        let neck=neck.ok_or("missing prescribed neck geometry")?;
        let port=air.as_ref().ok_or("missing neck mechanics")?.coupling.neck_radiation_port(0)?;
        Some(acoustics.ok_or("missing exterior boundary")?.with_sidewall_aperture(
            port,neck.azimuth_rad,neck.axial_position_m,depth)?)
    } else {acoustics};
    Ok(Experiment{flexible_sticks:shaft.ports,mute,system,force:vec![0.0;count],stick_weight,second_stick,observer_a:top,observer_b:bottom,pressure:Some(volume),acoustics,air})
}
// The stored drum areas encode compression, so positive contraction means
// positive internal pressure. The volume-spring Hamiltonian is unchanged.
fn cavity_pressure(volume:&VolumeSpring,state:&[f64])->f64 {
    (volume.bulk_modulus_pa/volume.volume_m3)*volume.areas.iter().enumerate()
        .map(|(i,a)|a*state[2*i]).sum::<f64>()
}
fn run()->Result<(),Error> {
    run_args(std::env::args().skip(1).collect())
}
fn run_args(mut raw_args:Vec<String>)->Result<(),Error> {
    if hihat::is_command(raw_args.first().map(String::as_str)) {return hihat::run(raw_args);}
    if specimen::export_command(&raw_args)? {return Ok(());}
    let head_stretching=nonlinear_snare::option(&mut raw_args)?;
    let mallet_paths=mallets::options(&mut raw_args)?;
    let head_relaxation_path=head_relaxation::option(&mut raw_args)?;
    let prescribed_vent=acoustics::aperture::option(&mut raw_args)?;
    let radiation_spec=acoustics::stereo::radiation_spec::option(&mut raw_args)?;
    let radiation_feedback=acoustics::stereo::feedback::option(&mut raw_args)?;
    let right_microphone=acoustics::stereo::option(&mut raw_args)?;
    let microphone_spec=acoustics::receivers::input::option(&mut raw_args)?;
    let analytic_newton=mechanics::analytic_option(&mut raw_args)?;
    let impact_substeps=mechanics::substeps_option(&mut raw_args)?;
    let prepared_nonlinear=mechanics::prepared_option(&mut raw_args)? || analytic_newton || impact_substeps.is_some();
    let distributed_cavity=cavity::option(&mut raw_args)?;
    let neck=cavity::neck_option(&mut raw_args)?;
    let cavity_drag=cavity::drag_option(&mut raw_args)?;
    let shell_path=specimen::option(&mut raw_args)?;
    let shell_mesh_path=specimen::mesh_option(&mut raw_args)?;
    let drum_path=drum_spec::option(&mut raw_args)?;
    let snare_path=snare::spec::option(&mut raw_args)?;
    let carrier_path=snare::carrier::option(&mut raw_args)?;
    let second=sticks::option(&mut raw_args)?;
    let shafts=shaft_playing::Selection::options(&mut raw_args)?;
    let mufflers=muffling::options(&mut raw_args)?;
    let compliant_mute=compliant_mute::option(&mut raw_args)?;
    let playing_force=mechanics::drive::option(&mut raw_args)?;
    let second_force=mechanics::drive::second_option(&mut raw_args)?;
    if second_force.is_some() && second.is_none() {
        return Err("--second-stick-force-file requires --second-stick-position-m X Y".into());
    }
    let driven=playing_force.is_some() || second_force.is_some() || compliant_mute.is_some() || carrier_path.is_some();
    let (args,stroke)=playing::parse(raw_args)?;
    if args.is_empty() || args.len()>6 {return Err("usage: percussion splash|drum [mechanics_steps]; splash-wav|drum-wav [audio_frames] [full_scale_pa]; splash-mic|drum-mic [audio_frames] [full_scale_pa] [x_m y_m z_m]; prepared drum: drum-modal[-wav|-mic] with the same arguments; see AUDIO.md, PREPARED.md and SNARES.md; snare[-off][-wav|-mic] adds explicit wire coupling; drum-stretch[-wav|-mic] adds geometric stretching; --head-stretching enables both nonlinear heads on snare[-off][-wav|-mic] without dropping wires or loss (see NONLINEAR_SNARE.md); --snare-spec wires.fsn supplies bank geometry, tension, damping, contact and optional wire stretching (see SNARE_SPEC.md); --snare-carrier input.fsc adds force-driven moving supports without resetting the wires (see SNARE_CARRIER.md); --strike-speed-m-s V and --strike-position-m X Y set physical launch inputs; --prepared-nonlinear prepares the unchanged splash/drum/drum-stretch model; --analytic-newton selects its analytic storage tangents (see ANALYTIC.md); --impact-substeps DEPTH ATTEMPTS adds bounded hard-impact recovery without changing the output clock (see SUBSTEPS.md); --cavity-modes adds distributed enclosed air to all drum/snare commands; --cavity-drag-per-s D supplies nonuniform acoustic momentum drag, and --cavity-neck radius_m length_eff_m resistance_Pa_s_m3 azimuth_rad z_m adds a vent to any drum/snare mechanics CSV (see CAVITY.md and SNARE_CAVITY.md); --drum-spec instrument.fsd supplies geometry, independent head materials/tensions/losses and the mesh/window (see DRUM_SPEC.md); --head-relaxation material.fshr supplies hereditary bending with zero separate head damping (see HEAD_RELAXATION.md); --microphone-right X,Y,Z adds a physical stereo receiver to -mic commands (see STEREO.md); --microphone-spec input.frm supplies geometry-aware close/directional microphones (see MICROPHONES.md); --compliant-mute file.fsm adds moving felt-pad squeeze/retract mechanics (see COMPLIANT_MUTE.md); --prescribed-vent-radiation adds explicit one-way neck-flow BEM radiation to a vented drum/snare audio command (see VENT_RADIATION.md); --shell-mesh input.fss supplies an explicit 3D shell and thickness/material fields; export-shell-mesh OUTPUT.fss [INPUT.profile] exports a profile before modal preparation (see SHELL_MESH.md); --radiation-spec input.fra selects the acoustic band, source-preserving boundary refinement and fitting/work limits (see RADIATION_BAND.md); --mallet-spec and --second-mallet-spec replace hard tips with supplied finite-area felt faces (see MALLETS.md); --radiation-feedback couples a passive BEM load into nonlinear-capable pressure playback (see RADIATION_FEEDBACK.md); --flexible-stick and --second-flexible-stick supply physical shaft geometry and hand/tip stations (see FLEXIBLE_STICKS.md)".into());}
    if let Some(spec)=&microphone_spec {spec.admit_command(&args[0],args.len()>3,right_microphone.is_some())?;}
    acoustics::stereo::feedback::admit_command(radiation_feedback,&args[0],neck.is_some())?;
    head_relaxation::admit_command(head_relaxation_path.is_some(),&args[0])?;
    let head_relaxation=head_relaxation_path.as_deref().map(head_relaxation::Spec::load).transpose()?;
    let mallets=mallets::Selection::load(mallet_paths)?;
    mallets.admit(&args[0],stroke,second)?;
    shafts.admit(&args[0],second,&mallets)?;
    let selected_snare=snare::spec::select(snare_path.as_deref(),&args[0])?;
    let (selected_snare,carrier)=snare::carrier::select(carrier_path.as_deref(),selected_snare)?;
    let has_carrier=carrier.is_some();
    let carrier_body=3+usize::from(second.is_some());
    let nonlinear_wires=selected_snare.is_some_and(|s|s.stretching.is_some());
    let nonlinear_instrument=head_stretching || nonlinear_wires || has_carrier
        || ((head_relaxation.is_some() || mallets.enabled() || shafts.enabled()) && selected_snare.is_some());
    nonlinear_snare::admit_command(head_stretching,&args[0])?;
    nonlinear_snare::admit_prepared_command(prepared_nonlinear,nonlinear_instrument,&args[0])?;
    if let Some(spec)=&compliant_mute {spec.admit_command(&args[0])?;}
    acoustics::stereo::admit_command(right_microphone,&args[0])?;
    if let Some(spec)=radiation_spec {spec.admit_command(&args[0])?;}
    cavity::admit_command(distributed_cavity,&args[0])?;
    vented::admit(neck,distributed_cavity,&args[0],prescribed_vent)?;
    cavity::admit_drag_command(cavity_drag,distributed_cavity,&args[0])?;
    specimen::admit_selection(shell_path.as_deref(),shell_mesh_path.as_deref(),&args[0])?;
    drum_spec::admit_command(drum_path.as_deref(),&args[0])?;
    sticks::admit_command(second.is_some(),&args[0])?;
    muffling::admit_command(&mufflers,&args[0])?;
    for spec in &mufflers {
        eprintln!("fixed viscous muffler: {:?}, xy_m={:?}, resistance_Ns_m={}; mechanical attachment, not a measured finger/gel or a timed choke",spec.surface,spec.position_m,spec.resistance_n_s_m);
    }
    let microphone=matches!(args[0].as_str(),"splash-mic"|"drum-mic"|"drum-modal-mic"|"snare-mic"|"snare-off-mic"|"drum-stretch-mic");
    let audio=microphone || matches!(args[0].as_str(),"splash-wav"|"drum-wav"|"drum-modal-wav"|"snare-wav"|"snare-off-wav"|"drum-stretch-wav");
    let stretching=head_stretching || matches!(args[0].as_str(),"drum-stretch"|"drum-stretch-wav"|"drum-stretch-mic");
    if args.len()>3 && (!microphone || args.len()!=6) {return Err("microphone position needs exactly x_m y_m z_m after frames and full-scale".into());}
    if !audio && args.len()>2 {return Err("mechanics CSV accepts only a step count".into());}
    let count=if args.len()>=2 {args[1].parse::<u64>()?}else if audio {48000}else{4096};
    let maximum=if audio {480000}else{1_000_000};
    if count==0 || count>maximum {return Err(format!("requested count must be 1..={maximum}").into());}
    let full_scale_pa=if args.len()>=3 {args[2].parse::<f64>()?}else{1.0};
    if !full_scale_pa.is_finite() || full_scale_pa<=0.0 {return Err("full_scale_pa must be positive and finite".into());}
    let receiver=if microphone {
        let position=if args.len()==6 {[args[3].parse::<f64>()?,args[4].parse::<f64>()?,args[5].parse::<f64>()?]}else{[0.08,0.05,0.35]};
        acoustics::Receiver::FinitePoint(position)
    }else{acoustics::Receiver::FarField([1.5,0.7,1.5])};
    let steps=if audio {count.checked_mul(acoustics::SUBSTEPS as u64).ok_or("sample budget overflow")?}else{count};
    let dt_s=if audio {acoustics::MECHANICAL_DT}else{2e-6};
    let supplied_shell=specimen::load_selection(shell_path.as_deref(),shell_mesh_path.as_deref())?;
    let supplied_drum=drum_path.as_deref().map(drum_spec::Spec::load).transpose()?;
    let drag_per_s=cavity_drag.unwrap_or(0.0);
    let experiment=match args[0].as_str(){
        "splash"|"splash-wav"|"splash-mic"=>splash_with_shafts(steps,dt_s,audio,stroke,supplied_shell,&mufflers,second,compliant_mute.as_ref(),&mallets,&shafts)?,
        "drum"|"drum-wav"|"drum-mic"=>drum_with_shafts(steps,dt_s,audio,false,None,false,stroke,distributed_cavity,neck,supplied_drum,second,&mufflers,drag_per_s,compliant_mute.as_ref(),prescribed_vent,head_relaxation.as_ref(),&mallets,&shafts)?,
        "drum-stretch"|"drum-stretch-wav"|"drum-stretch-mic"=>drum_with_shafts(steps,dt_s,audio,false,None,true,stroke,distributed_cavity,neck,supplied_drum,second,&mufflers,drag_per_s,compliant_mute.as_ref(),prescribed_vent,head_relaxation.as_ref(),&mallets,&shafts)?,
        "drum-modal"|"drum-modal-wav"|"drum-modal-mic"=>drum_with_shafts(steps,dt_s,audio,true,None,false,stroke,distributed_cavity,neck,supplied_drum,second,&mufflers,drag_per_s,None,prescribed_vent,head_relaxation.as_ref(),&mallets,&shafts)?,
        "snare"|"snare-wav"|"snare-mic"=>drum_with_shafts(steps,dt_s,audio,!nonlinear_instrument,selected_snare,head_stretching,stroke,distributed_cavity,neck,supplied_drum,second,&mufflers,drag_per_s,None,prescribed_vent,head_relaxation.as_ref(),&mallets,&shafts)?,
        "snare-off"|"snare-off-wav"|"snare-off-mic"=>drum_with_shafts(steps,dt_s,audio,!nonlinear_instrument,selected_snare,head_stretching,stroke,distributed_cavity,neck,supplied_drum,second,&mufflers,drag_per_s,None,prescribed_vent,head_relaxation.as_ref(),&mallets,&shafts)?,
        _=>return Err("unknown experiment".into()),
    };
    let receivers=match microphone_spec {
        Some(spec)=>spec.into_receivers(),
        None=>match right_microphone {
            Some(right)=>vec![receiver,acoustics::Receiver::FinitePoint(right)],None=>vec![receiver],
        },
    };
    let (mut experiment,loaded)=if radiation_feedback {
        let (e,bake)=acoustics::stereo::feedback::prepare(experiment,usize::try_from(count)?,full_scale_pa,
            &receivers,radiation_spec.unwrap_or_default(),&CancelGate::new_clock_free())?;
        (e,Some(bake))
    }else{(experiment,None)};
    if prepared_nonlinear {
        experiment.system=if analytic_newton {experiment.system.into_analytic_nonlinear()?}
            else {experiment.system.into_prepared_nonlinear()?};
        eprintln!("mechanical image: prepared nonlinear Gonzalez; analytic_newton={analytic_newton}; unchanged geometry, materials, felt history and clocks; real-time performance unqualified");
    }
    if let Some(bounds)=impact_substeps {
        experiment.system=experiment.system.with_impact_substeps(bounds)?;
        eprintln!("internal impact refinement: depth={}, maximum solve attempts={} per mechanical output tick; unchanged output/force clocks, transactional state/history; no temporal accuracy or real-time claim",bounds.max_depth,bounds.max_attempts);
    }
    let (mut inputs,hands)=shaft_playing::player_inputs(&experiment,[playing_force,second_force])?;
    if let Some(spec)=carrier {
        inputs.push(spec.into_input(&experiment.system,carrier_body)?);
    }
    if let Some(spec)=compliant_mute {
        inputs.extend(spec.into_inputs(experiment.mute.as_ref().ok_or("missing moving mute assembly")?)?);
    }
    if !inputs.is_empty() || !hands.is_empty() {
        eprintln!("external physical performance: {} independent SI force programs, one accepted mechanical clock; physical hand rows for flexible shafts, no state reset; see DRIVE.md and FLEXIBLE_STICKS.md",inputs.len()+hands.len());
        experiment.system=experiment.system.with_player_drives(inputs,hands,dt_s,steps,experiment.force.len())?;
    }
    for (hand,p) in experiment.flexible_sticks.iter().enumerate() {
        if let Some(p)=p {eprintln!("physical shaft {}: {} elastic modes, rigid-only launch, tip contact and distinct hand-force port; no direct acoustic source",hand+1,p.elastic_modes());}
    }
    eprintln!("physical stroke: speed_m_s={}, explicit_xy_m={:?}; no output normalization or pitch control",stroke.speed_m_s,stroke.position_m);
    let stdout=std::io::stdout();let mut out=std::io::BufWriter::new(stdout.lock());
    if audio {
        // Render and admit the complete candidate before writing a WAV header.
        let wav=match loaded {
            Some(bake)=>bake.render(&mut experiment,usize::try_from(count)?,full_scale_pa,&CancelGate::new_clock_free())?,
            None=>acoustics::stereo::render_receivers_with_spec(&mut experiment,usize::try_from(count)?,
                full_scale_pa,&receivers,radiation_spec.unwrap_or_default(),&CancelGate::new_clock_free())?,
        };
        out.write_all(&wav)?;out.flush()?;return Ok(());
    }
    let gate=CancelGate::new_clock_free();
    let extra=if stretching {",batter_slope,resonant_slope,head_stretching_energy_j"}else{""};
    let wire_columns=if nonlinear_wires || has_carrier {",snare_max_slope_bound,snare_max_tension_n,snare_stretching_energy_j"}else{""};
    let carrier_columns=if has_carrier {",snare_carrier_position_m,snare_carrier_velocity_m_s"}else{""};
    let air_columns=if distributed_cavity {",cavity_point_a_pa,cavity_point_b_pa"}else{""};
    let neck_columns=if neck.is_some() {",neck_volume_m3,neck_flow_m3_s,neck_pressure_pa,neck_loss_power_w"}else{""};
    let drive_columns=if driven {",player_work_j"}else{""};
    let memory_columns=if head_relaxation.is_some() {",head_memory_energy_j,head_relaxation_power_w"}else{""};
    let stick_columns=if second.is_some() {",stick_1_displacement_m,stick_1_velocity_m_s,stick_2_displacement_m,stick_2_velocity_m_s"}else{""};
    write!(out,"time_s,point_a_displacement_m,point_a_velocity_m_s,point_b_displacement_m,cavity_internal_pa,total_energy_j,felt_crush_j,loss_j,balance_j{extra}{wire_columns}{carrier_columns}{air_columns}{neck_columns}{drive_columns}{stick_columns}{memory_columns}")?;
    if let Some(mute)=&experiment.mute {mute.header(&mut out)?;}
    shaft_playing::header(&experiment,&mut out)?;
    writeln!(out)?;
    for _ in 0..steps {
        let f=experiment.system.step(&experiment.force,&gate)?;let x=experiment.system.state();
        let displacement=|weights:&[f64]|weights.iter().enumerate().map(|(i,b)|b*x[2*i]).sum::<f64>();
        let velocity=experiment.observer_a.iter().enumerate().map(|(i,b)|b*x[2*i+1]).sum::<f64>();
        let pressure=match &experiment.air {
            Some(air)=>air.uniform_pressure(x)?,
            None=>experiment.pressure.as_ref().map_or(0.0,|v|cavity_pressure(v,x)),
        };
        write!(out,"{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
            f.time_s,displacement(&experiment.observer_a),velocity,displacement(&experiment.observer_b),pressure,
            f.stored_energy_j,f.felt_crush_loss_j,f.dissipated_energy_j,f.balance_residual_j)?;
        if stretching {
            let a=experiment.system.membrane_observation(1).ok_or("missing batter stretching state")?;
            let b=experiment.system.membrane_observation(2).ok_or("missing resonant-head stretching state")?;
            write!(out,",{:.17e},{:.17e},{:.17e}",a.maximum_slope,b.maximum_slope,a.stretching_energy_j+b.stretching_energy_j)?;
        }
        if has_carrier {
            let o=snare::carrier::observe(&experiment.system,carrier_body).ok_or("missing accepted snare carrier state")?;
            write!(out,",{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",o.maximum_wire_slope,
                o.maximum_wire_tension_n,o.wire_stretching_energy_j,o.position_m,o.velocity_m_s)?;
        } else if nonlinear_wires {
            let spec=selected_snare.ok_or("missing selected wire bank")?;
            let first_body=3+usize::from(second.is_some());
            let (mut slope,mut tension,mut energy)=(0.0_f64,0.0_f64,0.0_f64);
            for body in first_body..first_body+spec.strands {
                let value=snare::observe(&experiment.system,body).ok_or("missing stretching wire state")?;
                slope=slope.max(value.slope_bound);tension=tension.max(value.tension_n);
                energy+=value.stretching_energy_j;
            }
            write!(out,",{slope:.17e},{tension:.17e},{energy:.17e}")?;
        }
        if let Some(air)=&experiment.air {
            let (a,b)=air.points(x)?;write!(out,",{a:.17e},{b:.17e}")?;
            if air.coupling.neck_count()>0 {
                let n=air.coupling.neck_observation(x,0)?;
                write!(out,",{:.17e},{:.17e},{:.17e},{:.17e}",n.displaced_volume_m3,
                    n.volume_flow_m3_s,n.driving_pressure_pa,n.dissipated_power_w)?;
            }
        }
        if driven {write!(out,",{:.17e}",f.supplied_work_j)?;}
        if experiment.second_stick.is_some() {
            let (a,av)=shaft_playing::tip_motion(&experiment,0)?;
            let (b,bv)=shaft_playing::tip_motion(&experiment,1)?;
            write!(out,",{a:.17e},{av:.17e},{b:.17e},{bv:.17e}")?;
        }
        if head_relaxation.is_some() {
            let memory=head_relaxation::observation(&experiment.system);
            write!(out,",{:.17e},{:.17e}",memory.stored_energy_j,memory.dissipated_power_w)?;
        }
        if let Some(mute)=&experiment.mute {mute.row(&experiment.system,&mut out)?;}
        shaft_playing::row(&experiment,&mut out)?;
        writeln!(out)?;
    }
    out.flush()?;Ok(())
}
fn main(){if let Err(e)=run(){eprintln!("percussion reference refused: {e}");std::process::exit(1);}}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn a_later_player_force_drives_real_drum_contact_without_resetting_ringdown() {
        let gate=CancelGate::new_clock_free();
        for prepared in [false,true] {
            let mut manual=drum(256,2e-6,false,prepared).unwrap();
            let mut driven=drum(256,2e-6,false,prepared).unwrap();
            let initial=driven.system.state().to_vec();
            let text="0,0\n0.000256,0\n0.000384,2\n0.000512,0";
            let mut staging=mechanics::drive::StickDrive::new(
                mechanics::drive::Program::parse(text).unwrap(),2e-6,256,
                manual.stick_weight,manual.force.len()).unwrap();
            driven.system=driven.system.with_stick_drive(
                mechanics::drive::Program::parse(text).unwrap(),2e-6,256,
                driven.stick_weight,driven.force.len()).unwrap();
            assert_eq!(driven.system.state(),initial);
            let mut head_motion=0.0_f64;let mut player_work=0.0_f64;
            for step in 0..256 {
                let f=driven.system.step(&driven.force,&gate).unwrap();
                manual.system.step(staging.forces(&manual.force).unwrap(),&gate).unwrap();
                staging.accept();
                assert_eq!(driven.system.state(),manual.system.state());
                if step<128 {assert_eq!(f.supplied_work_j,0.0);}
                player_work+=f.supplied_work_j.abs();
                head_motion=head_motion.max(driven.system.state()[2..].iter().map(|x|x.abs()).fold(0.0_f64,f64::max));
                assert!(f.balance_residual_j.abs()<1e-7);
            }
            assert!(head_motion>0.0 && player_work>0.0);
        }
    }
    #[test]
    fn downward_batter_motion_compresses_air_and_bottom_motion_releases_it() {
        let v=VolumeSpring{bulk_modulus_pa:100.0,volume_m3:2.0,areas:vec![0.0,3.0,-3.0]};
        assert!((cavity_pressure(&v,&[0.,0.,0.2,0.,0.,0.])-30.0).abs()<1e-13);
        assert!((cavity_pressure(&v,&[0.,0.,0.,0.,0.2,0.])+30.0).abs()<1e-13);
        assert_eq!(cavity_pressure(&v,&[0.,0.,0.2,0.,0.2,0.]),0.0);
    }
}

#[cfg(test)]
mod stretching_tests {
    use super::*;

    #[test]
    fn stronger_stick_launch_changes_kinetic_energy_not_mass_geometry_or_gap() {
        let (a,wa)=stick().unwrap(); let (b,wb)=stick_with_speed(1.6).unwrap();
        assert_eq!(wa.to_bits(),wb.to_bits());
        assert_eq!(a.initial[0].displacement_m_sqrt_kg.to_bits(),b.initial[0].displacement_m_sqrt_kg.to_bits());
        assert_eq!(b.initial[0].velocity_m_sqrt_kg_per_s,2.0*a.initial[0].velocity_m_sqrt_kg_per_s);
        assert_eq!(a.damping_per_s,b.damping_per_s);
        assert!(stick_with_speed(f64::NAN).is_err());
    }

    #[test]
    fn stretching_drum_changes_actual_head_motion_without_changing_the_basis_or_air() {
        let stroke=Stroke{speed_m_s:4.0,position_m:Some([0.06,0.01])};
        let mut linear=drum_with_playing(64,2e-6,false,false,None,false,stroke).unwrap();
        let mut nonlinear=drum_with_playing(64,2e-6,false,false,None,true,stroke).unwrap();
        assert_eq!(linear.system.state(),nonlinear.system.state());
        assert_eq!(linear.observer_a,nonlinear.observer_a);
        assert_eq!(linear.observer_b,nonlinear.observer_b);
        assert_eq!(linear.pressure.as_ref().unwrap().areas,nonlinear.pressure.as_ref().unwrap().areas);
        assert!(linear.system.membrane_observation(1).is_none());
        assert_eq!(nonlinear.system.membrane_observation(1).unwrap().stretching_energy_j,0.0);
        let gate=CancelGate::new_clock_free(); let mut stretch=0.0_f64; let mut changed=0.0_f64;
        for _ in 0..64 {
            linear.system.step(&linear.force,&gate).unwrap();
            nonlinear.system.step(&nonlinear.force,&gate).unwrap();
            let state=nonlinear.system.membrane_observation(1).unwrap();
            assert!(state.maximum_slope<=0.2); stretch=stretch.max(state.stretching_energy_j);
            for (&a,&b) in linear.system.state().iter().zip(nonlinear.system.state()) { changed=changed.max((a-b).abs()); }
        }
        assert!(stretch>0.0 && changed>1e-16,"stretching must enter the actual contact/air/mechanics solve");
        assert!(drum_with_playing(64,2e-6,false,true,None,true,stroke).is_err());
    }
}
