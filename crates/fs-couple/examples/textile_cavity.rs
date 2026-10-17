//! Two source-backed textile walls exchanging work through one material gas.
//! No constitutive or integration law lives in this example. Gas parsing,
//! membrane equilibrium/tangents, eigenanalysis and the coupled step all use
//! existing owners. Inputs define a fixed prestretched reference enclosure;
//! motion is small-signal about it, not a moving-mesh or wrinkled-cloth model.
use std::{collections::BTreeSet,error::Error,io::Read,str::SplitAsciiWhitespace};
use fs_couple::{cavity_modal::{HelmholtzModalCavity,CavityNeckLoss,shared::{SharedCavity,CavitySurface}},
    modal_acoustic_time::{ModalAcousticMode,ModalAcousticState,ModalAcousticTimeBudget,ModalAcousticTimeModel},
    render::schedule::force::coupled::ModalCouplingConfig};
use fs_exec::CancelGate;
use fs_material::{gas::GasState,gas_species::input::read_gas_state,
    sheet::{SheetCapabilities,resolve_sheet},state_point::MaterialPropertySelection};
use fs_matdb::{MaterialCard,MaterialStateId,NormalizedPack,QueryPoint};
use fs_math::c64::C64;
use fs_scenario::HelmholtzCavity;
use fs_solid::{membrane::{MembranePatch,MembraneBudget,MembraneDomain,MembraneSolveSettings,
    vibration::{MembraneVibration,VibrationSettings}},mesh2::{Mesh2,Patch}};
use fs_solver::nonlinear::NewtonKrylovConfig;
type Result<T> = std::result::Result<T,Box<dyn Error>>;

struct Input {
    gas:GasState,geometry:HelmholtzCavity,resistance:f64,temperature:f64,stretch:f64,grain:f64,strain:f64,
    rate:u32,frames:usize,pulse_frames:usize,force:f64,max_displacement:f64,max_volume_fraction:f64,
}
fn word<'a>(f:&mut SplitAsciiWhitespace<'a>)->Result<&'a str> {f.next().ok_or_else(||"missing cavity field".into())}
fn scalar(f:&mut SplitAsciiWhitespace<'_>)->Result<f64> {
    let x:f64=word(f)?.parse()?;if !x.is_finite() {return Err("finite cavity scalar required".into());}Ok(x)
}
fn label(f:&mut SplitAsciiWhitespace<'_>,expected:&str)->Result<()> {
    if word(f)?!=expected {return Err(format!("expected {expected} record").into());}Ok(())
}
impl Input {
    fn parse(text:&str)->Result<Self> {
        if text.len()>16384 {return Err("cavity input exceeds 16 KiB".into());}
        let mut f=text.split_ascii_whitespace();label(&mut f,"textile-cavity-v1")?;label(&mut f,"gas")?;
        let gas=read_gas_state(&mut f,None)?;label(&mut f,"cavity")?;
        // The input predates DEF-09: this fixture's neck is a hole through the
        // cavity's own wall, so both ends are declared flanged here.
        let geometry=HelmholtzCavity {volume_m3:scalar(&mut f)?,neck_radius_m:scalar(&mut f)?,neck_length_m:scalar(&mut f)?,
            inner_end:fs_scenario::NeckEnd::Flanged,outer_end:fs_scenario::NeckEnd::Flanged};
        let resistance=scalar(&mut f)?;label(&mut f,"surface")?;
        let temperature=scalar(&mut f)?;let stretch=scalar(&mut f)?;let grain=scalar(&mut f)?;let strain=scalar(&mut f)?;
        label(&mut f,"drive")?;let rate:u32=word(&mut f)?.parse()?;
        let frames:usize=word(&mut f)?.parse()?;let pulse_frames:usize=word(&mut f)?.parse()?;let force=scalar(&mut f)?;
        label(&mut f,"limits")?;let max_displacement=scalar(&mut f)?;let max_volume_fraction=scalar(&mut f)?;
        if f.next().is_some() {return Err("extra cavity input fields".into());}
        if !(8000..=192000).contains(&rate) || frames==0 || frames>192000 || pulse_frames==0 || pulse_frames>frames
            || temperature<=0.0 || stretch<=1.0 || strain<=0.0 || max_displacement<=0.0
            || max_displacement>0.001 || max_volume_fraction<=0.0 || max_volume_fraction>0.01 {
            return Err("invalid bounded clock, tensile prestretch or small-signal limits".into());
        }
        // Admission belongs to the same physical source later used by the run.
        HelmholtzModalCavity::new(geometry,gas,CavityNeckLoss::Resistance(resistance))?;
        Ok(Self {gas,geometry,resistance,temperature,stretch,grain,strain,rate,frames,pulse_frames,force,
            max_displacement,max_volume_fraction})
    }
}

fn patch(card:&MaterialCard,i:&Input)->Result<(MembranePatch,Vec<(usize,f64)>)> {
    let point=QueryPoint::new().with("T",i.temperature)?.with("ftx_characterization",0.0)?;
    let sheet=resolve_sheet(card,&point,SheetCapabilities {membrane:true,bending:false},MaterialPropertySelection::SingleClaimOnly)?;
    let mesh=Mesh2::triangles(0.1,0.1,2,2);
    let boundary:BTreeSet<_>=[Patch::Left,Patch::Right,Patch::Top,Patch::Bottom].into_iter().flat_map(|p|mesh.patch_nodes(p)).collect();
    let mut fixed=Vec::new();
    for node in boundary {
        let [x,y]=mesh.nodes[node];fixed.extend([(3*node,(i.stretch-1.0)*(x-0.05)),
            (3*node+1,(i.stretch-1.0)*(y-0.05)),(3*node+2,0.0)]);
    }
    Ok((MembranePatch::try_new(mesh,vec![sheet;8],vec![i.grain;8],
        MembraneDomain {minimum_area_ratio:0.01,maximum_abs_strain:i.strain},
        MembraneBudget {max_nodes:9,max_triangles:8},||false)?,fixed))
}
fn vibration<'a>(p:&'a MembranePatch,bc:&[(usize,f64)])->Result<MembraneVibration<'a>> {
    Ok(p.solve_vibration(&[[0.0;3];9],bc,&[[0.0;3];9],MembraneSolveSettings {
        newton:NewtonKrylovConfig {linear_restart:12,absolute_tolerance:1e-9,relative_tolerance:1e-10,..Default::default()},
        max_iterations:80,max_evaluations:2000,max_krylov_restart:12,
    },VibrationSettings {max_free_dofs:3,relative_residual_tolerance:1e-8},||false)?)
}
fn surface(v:&MembraneVibration<'_>,sign:f64,i:&Input)->Result<CavitySurface> {
    // Bound each nodal Cartesian increment by ||phi_row||_1 * max|q|.
    // The original modal owner enforces this conservative coordinate limit
    // BEFORE publication; there is no post-step clamp or trajectory rescaling.
    let row_norm=(0..v.free_dofs().len()).map(|r|v.modes().iter().map(|m|m.phi[r].abs()).sum::<f64>()).fold(0.0_f64,f64::max);
    let modal_limit=i.max_displacement/row_norm;
    if !modal_limit.is_finite() || modal_limit<=0.0 {return Err("unrepresentable modal displacement bound".into());}
    let mut budget=ModalAcousticTimeBudget::audible_reference();
    budget.maximum_abs_displacement_m_sqrt_kg=modal_limit;
    let model=ModalAcousticTimeModel::try_new(i.rate,v.modes().iter().map(|m|ModalAcousticMode {
        angular_frequency_rad_s:m.lambda.sqrt(),damping_ratio:0.0,pressure_per_modal_velocity:C64::ZERO,
    }).collect(),budget)?;
    Ok(CavitySurface {model,outward_volume_shapes:v.volume_coupling().iter().map(|a|sign*a).collect()})
}
#[derive(Debug,Clone,PartialEq)]
struct Row {time:f64,a:f64,b:f64,pressure:f64,flow:f64,volume_fraction:f64,energy:f64,work:f64,loss:f64,residual:f64}
struct Run {rows:Vec<Row>,frequencies:[Vec<f64>;2],masses:[f64;2]}
fn simulate(a:&MaterialCard,b:&MaterialCard,i:&Input)->Result<Run> {
    let (pa,bca)=patch(a,i)?;let (pb,bcb)=patch(b,i)?;
    let va=vibration(&pa,&bca)?;let vb=vibration(&pb,&bcb)?;
    let source=HelmholtzModalCavity::new(i.geometry,i.gas,CavityNeckLoss::Resistance(i.resistance))?;
    // E_total >= E_compression. This conservative total-energy ceiling enforces
    // the supplied linearized gas-volume limit inside the original transaction.
    let extension=i.max_volume_fraction*i.geometry.volume_m3/source.neck_area_m2();
    let max_energy=0.5*source.compression_stiffness_n_m()*extension*extension;
    if !max_energy.is_finite() || max_energy<=0.0 {return Err("unrepresentable small-signal energy ceiling".into());}
    let config=ModalCouplingConfig {max_modes:7,max_connections:2,max_setup_terms:128,nyquist_guard_fraction:0.9,
        maximum_total_energy_j:max_energy,maximum_abs_pressure_pa:1e6,maximum_abs_connection_force_n:1e5,
        solve_relative_tolerance:1e-10,energy_absolute_tolerance_j:1e-12,energy_relative_tolerance:1e-8};
    let gate=CancelGate::new_clock_free();
    // The two reference sheets are opposite walls, so local +z has opposite
    // outward signs. All frequencies/masses/shapes come from their own sources.
    let mut coupled=SharedCavity::new(source,vec![surface(&va,1.0,i)?,surface(&vb,-1.0,i)?],i.rate,
        ModalAcousticState::default(),ModalAcousticTimeBudget::audible_reference(),config,&gate)?;
    let mut nodal=[[0.0;3];9];nodal[4][2]=i.force;
    let mut drive=va.project_nodal_force(&nodal)?;drive.extend(vec![0.0;vb.modes().len()]);
    let zero=vec![0.0;drive.len()];let mut work=0.0;let mut loss=0.0;
    let stride=(i.frames/200).max(1);let mut rows=Vec::new();
    for sample in 0..i.frames {
        let frame=coupled.step(if sample<i.pulse_frames {&drive}else{&zero},0.0,&gate)?;
        work+=frame.external_work_j;loss+=frame.component_dissipation_j+frame.connection_dissipation_j;
        if (sample+1)%stride==0 || sample+1==i.frames {
            let observation=coupled.observe()?;let energy=coupled.network().total_energy_j()?;
            let q=|component:usize|coupled.network().components()[component].states().iter()
                .map(|s|s.displacement_m_sqrt_kg).collect::<Vec<_>>();
            rows.push(Row {time:(sample+1) as f64/i.rate as f64,a:va.reconstruct(&q(0))?[4][2],b:vb.reconstruct(&q(1))?[4][2],
                pressure:observation.pressure_increment_pa,flow:observation.neck_outflow_m3_s,
                volume_fraction:observation.fractional_volume_change,energy,work,loss,residual:energy+loss-work});
        }
    }
    Ok(Run {rows,frequencies:[va.modes().iter().map(|m|m.lambda.sqrt()/std::f64::consts::TAU).collect(),
        vb.modes().iter().map(|m|m.lambda.sqrt()/std::f64::consts::TAU).collect()],masses:[pa.reference_mass_kg(),pb.reference_mass_kg()]})
}
fn bytes(path:&str,limit:u64)->Result<Vec<u8>> {
    let mut b=Vec::new();std::fs::File::open(path)?.take(limit+1).read_to_end(&mut b)?;
    if b.len() as u64>limit {return Err("input file exceeds its size limit".into());}Ok(b)
}
fn card(pack:&NormalizedPack)->Result<MaterialCard> {
    Ok(MaterialCard::assemble(MaterialStateId {chemistry:pack.pack_id().into(),phase:"characterized sheet".into(),
        process:"imported fixed specimen/conditioning".into(),revision:0},pack.claims().clone(),vec![])?)
}
fn main()->Result<()> {
    let args:Vec<_>=std::env::args().collect();if args.len()!=4 {return Err("usage: textile_cavity PACK_A PACK_B CAVITY.txt".into());}
    let a=NormalizedPack::from_bytes(&bytes(&args[1],64*1024*1024)?)?;
    let b=NormalizedPack::from_bytes(&bytes(&args[2],64*1024*1024)?)?;
    let i=Input::parse(std::str::from_utf8(&bytes(&args[3],16384)?)?)?;
    // Retain BOTH original packs, cards and derived structural owners throughout
    // simulation. Emit no partial CSV on a refused run.
    let run=simulate(&card(&a)?,&card(&b)?,&i)?;
    println!("# packs={:?},{:?}; gas={:?}",a.pack_id(),b.pack_id(),i.gas);
    println!("# reference_masses_kg={:?}; uncoupled_prestressed_modes_Hz={:?}",run.masses,run.frequencies);
    println!("# fixed-basis small-signal membrane/cavity; zero structural loss; no radiation or thermal evolution");
    println!("time_s,wall_a_increment_m,wall_b_increment_m,cavity_pressure_Pa,neck_outflow_m3_s,volume_fraction,stored_J,input_work_J,loss_J,residual_J");
    for r in run.rows {println!("{:.9},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e}",
        r.time,r.a,r.b,r.pressure,r.flow,r.volume_fraction,r.energy,r.work,r.loss,r.residual);}
    Ok(())
}

#[cfg(test)]
#[path="textile_cavity/tests.rs"]
mod tests;
