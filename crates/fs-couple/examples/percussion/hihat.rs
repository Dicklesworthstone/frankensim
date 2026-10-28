//! Paired nonlinear shells with axial pedal mechanics and distributed contact.
//! This is a fixed-reference, small-slope hi-hat research image, not a sample mix.
use super::*;
use fs_plate::shell::survey::MeshShell;
use fs_plate::shell::reduction::radiation::{ShellFace,ShellRadiationSurface,RadiationSurfaceBudget};

#[path="hihat_input.rs"]
mod input;
use input::{Spec,Mount};
#[path="hihat_play.rs"]
mod play;
pub use play::{is_command,run};
#[path="hihat_squeeze.rs"]
mod squeeze;
#[path="hihat_targets.rs"]
mod targets;
use fs_couple::render::plate::impact::striker::flexible::{FlexibleStriker,StrikerPorts};

// One independently prepared physical shell and its exact inner-skin chart.
struct Shell {
    mesh:MeshShell,
    reduction:ShellReduction,
    skin:ShellRadiationSurface,
    inner_nodes:Vec<[f64;3]>,
    outer_nodes:Vec<[f64;3]>,
}
impl Shell {
    fn new(spec:&specimen::Specimen,dt:f64)->Result<Self,Error> {
        Self::new_with_rocking(spec,dt,None)
    }
    fn new_with_rocking(spec:&specimen::Specimen,dt:f64,rocking:Option<f64>)->Result<Self,Error> {
        let (mesh,reduction)=match rocking {
            Some(limit)=>shell_prepare::prepare_rocking(spec,dt,limit)?,
            None=>shell_prepare::prepare(spec,dt)?,
        };
        let skin=reduction.radiation_surface(&mesh.nodal_thickness_m,
            RadiationSurfaceBudget{max_panels:2048,max_panel_modes:65536})?;
        // The negative skin reverses vertex order; restore SOURCE connectivity
        // before locating contact barycentrics. Do not project onto midsurface XY.
        let mut inner_nodes=vec![[0.;3];mesh.mesh.nodes.len()];
        let mut outer_nodes=vec![[0.;3];mesh.mesh.nodes.len()];
        for (f,t) in mesh.mesh.tris.iter().enumerate() {
            for a in 0..3 {
                inner_nodes[t[a]]=skin.triangles()[2*f+1][2-a];
                outer_nodes[t[a]]=skin.triangles()[2*f][a];
            }
        }
        Ok(Self{mesh,reduction,skin,inner_nodes,outer_nodes})
    }
    fn port(&self,p:[f64;2],face:ShellFace)->Result<fs_plate::shell::reduction::radiation::ShellSurfacePort,Error> {
        let nodes=if face==ShellFace::Positive {&self.outer_nodes}else{&self.inner_nodes};
        let (f,b)=playing::shell_location(nodes,&self.mesh.mesh.tris,p)?;
        Ok(self.reduction.surface_point_port(&self.mesh.nodal_thickness_m,f,b,face,[0.,0.,1.])?)
    }
}
struct Pair {
    experiment:Experiment,
    pedal:sticks::Port,
    collision:Obstacle,
    upper_modes:std::ops::Range<usize>,
    lower_modes:std::ops::Range<usize>,
    flexible_sticks:[Option<StrikerPorts>;2],
    targets:targets::Selection,
}

// Mounts are opposed compression-only washers, not an imposed shell position.
// Upper local +z is up; pedal travel is down, so relative offset is z_top + x_p.
fn washers(s:&Shell,m:&Mount,start:usize,total:usize,pedal:Option<sticks::Port>)->Result<Vec<FeltPad>,Error> {
    let mut pads=Vec::with_capacity(6);
    for i in 0..3 {
        let a=2.*std::f64::consts::PI*i as f64/3.;
        let (f,b)=playing::shell_location(&s.mesh.mesh.nodes,&s.mesh.mesh.tris,
            [m.radius*a.cos(),m.radius*a.sin()])?;
        let row=s.reduction.point_port(f,b,[0.,0.,1.])?;
        for sign in [-1.,1.] {
            let mut weights=vec![0.;total];
            for (k,w) in row.iter().enumerate(){weights[start+k]=sign*w;}
            if let Some(p)=pedal {weights[p.coordinate]=sign*p.weight;}
            pads.push(FeltPad{weights,area_m2:m.area/3.,thickness_m:m.thickness,
                precompression_m:m.precompression,law:m.law.clone(),prior_maximum_strain:m.prior,
                // Whole-face parallel area split: both K and eta scale, keeping tau.
                creep:vec![KelvinBranch{stiffness_n_m:m.k/3.,viscosity_n_s_m:m.eta/3.}]});
        }
    }
    Ok(pads)
}

fn collision(spec:&Spec,upper:&Shell,lower:&Shell,lo:usize,total:usize)->Result<Obstacle,Error> {
    let mut rows=Vec::with_capacity(spec.sites.len()*total);let mut gaps=Vec::new();
    for &(x,y,_) in &spec.sites {
        // Lower is placed by a proper pi rotation about x: (x,y,z)->(x,-y,-z).
        let a=upper.port([x,y],ShellFace::Negative)?;
        let b=lower.port([x,-y],ShellFace::Negative)?;
        let gap=spec.separation+a.position_m[2]+b.position_m[2];
        if !gap.is_finite() || gap<=0. {return Err("paired skins must begin separated at every contact site".into());}
        let mut row=vec![0.;total];
        for (k,w) in a.weights.iter().enumerate(){row[1+k]=-w;}
        for (k,w) in b.weights.iter().enumerate(){row[lo+k]=-w;}
        rows.extend(row);gaps.push(gap);
    }
    Ok(Obstacle::new(rows,spec.sites.len(),total,gaps,spec.sites.iter().map(|s|s.2).collect(),
        spec.contact[0],spec.contact[1],"supplied effective normal inter-cymbal contact; fixed reference sites, not measured or collision-detected".into())?
        .with_internal_loss(spec.contact[2])?)
}

fn build(spec:&Spec,upper:&specimen::Specimen,lower:&specimen::Specimen,stroke:Stroke,
    second:Option<Stroke>,steps:u64,dt:f64,audio:bool)->Result<Pair,Error> {
    build_with_squeeze(spec,upper,lower,stroke,second,steps,dt,audio,None)
}
#[allow(clippy::too_many_arguments)]
fn build_with_squeeze(spec:&Spec,upper:&specimen::Specimen,lower:&specimen::Specimen,stroke:Stroke,
    second:Option<Stroke>,steps:u64,dt:f64,audio:bool,film:Option<&squeeze::Config>)->Result<Pair,Error> {
    build_with_strikers(spec,upper,lower,stroke,second,steps,dt,audio,film,[None,None])
}
#[allow(clippy::too_many_arguments)]
fn build_with_strikers(spec:&Spec,upper:&specimen::Specimen,lower:&specimen::Specimen,stroke:Stroke,
    second:Option<Stroke>,steps:u64,dt:f64,audio:bool,film:Option<&squeeze::Config>,
    flexible:[Option<&FlexibleStriker>;2])->Result<Pair,Error> {
    let shafts=shaft_playing::Selection{first:flexible[0].cloned(),second:flexible[1].cloned()};
    build_with_mallets(spec,upper,lower,stroke,second,steps,dt,audio,film,&shafts,
        &mallets::Selection::default())
}
#[allow(clippy::too_many_arguments)]
fn build_with_mallets(spec:&Spec,upper:&specimen::Specimen,lower:&specimen::Specimen,stroke:Stroke,
    second:Option<Stroke>,steps:u64,dt:f64,audio:bool,film:Option<&squeeze::Config>,
    shafts:&shaft_playing::Selection,mallets:&mallets::Selection)->Result<Pair,Error> {
    build_with_material(spec,upper,lower,stroke,second,steps,dt,audio,film,shafts,mallets,None)
}
#[allow(clippy::too_many_arguments)]
fn build_with_material(spec:&Spec,upper:&specimen::Specimen,lower:&specimen::Specimen,stroke:Stroke,
    second:Option<Stroke>,steps:u64,dt:f64,audio:bool,film:Option<&squeeze::Config>,
    shafts:&shaft_playing::Selection,mallets:&mallets::Selection,material:Option<&shell_relaxation::Spec>)->Result<Pair,Error> {
    build_with_mute(spec,upper,lower,stroke,second,steps,dt,audio,film,shafts,mallets,material,None)
}
#[allow(clippy::too_many_arguments)]
fn build_with_mute(spec:&Spec,upper:&specimen::Specimen,lower:&specimen::Specimen,stroke:Stroke,
    second:Option<Stroke>,steps:u64,dt:f64,audio:bool,film:Option<&squeeze::Config>,
    shafts:&shaft_playing::Selection,mallets:&mallets::Selection,material:Option<&shell_relaxation::Spec>,
    mute:Option<&compliant_mute::Spec>)->Result<Pair,Error> {
    build_with_targets(spec,upper,lower,stroke,second,steps,dt,audio,film,shafts,mallets,material,
        mute,targets::Selection::default())
}
#[allow(clippy::too_many_arguments)]
fn build_with_targets(spec:&Spec,upper:&specimen::Specimen,lower:&specimen::Specimen,stroke:Stroke,
    second:Option<Stroke>,steps:u64,dt:f64,audio:bool,film:Option<&squeeze::Config>,
    shafts:&shaft_playing::Selection,mallets:&mallets::Selection,material:Option<&shell_relaxation::Spec>,
    mute:Option<&compliant_mute::Spec>,targets:targets::Selection)->Result<Pair,Error> {
    build_with_rocking(spec,upper,lower,stroke,second,steps,dt,audio,film,shafts,mallets,material,
        mute,targets,None)
}
#[allow(clippy::too_many_arguments)]
fn build_with_rocking(spec:&Spec,upper:&specimen::Specimen,lower:&specimen::Specimen,stroke:Stroke,
    second:Option<Stroke>,steps:u64,dt:f64,audio:bool,film:Option<&squeeze::Config>,
    shafts:&shaft_playing::Selection,mallets:&mallets::Selection,material:Option<&shell_relaxation::Spec>,
    mute:Option<&compliant_mute::Spec>,targets:targets::Selection,rocking:Option<f64>)->Result<Pair,Error> {
    targets.admit(second)?;
    // The hi-hat's own row owns the air; a shell that declares one must agree.
    if [upper,lower].iter().any(|s|s.ambient.as_ref().is_some_and(|a|a.gas()!=spec.ambient.gas())) {
        return Err("hi-hat shells and the hi-hat description declare different ambient gases".into());
    }
    let ambient=super::ambient::effective(&spec.ambient);
    eprintln!("{}",ambient.receipt());
    mallets.admit("hihat",stroke,second)?;
    if let Some(mute)=mute {mute.admit_command("hihat")?;}
    if shafts.second.is_some() && second.is_none() {
        return Err("a flexible second stick requires its own strike position".into());
    }
    for (shaft,mallet) in [shafts.first.as_ref(),shafts.second.as_ref()].into_iter()
        .zip([mallets.first.as_ref(),mallets.second.as_ref()]) {
        if let Some(mallet)=mallet {mallet.admit_shaft(shaft.is_some())?;}
    }
    // Both complete opposed mounts already own twelve independent felt sites.
    // One four-site head fits the existing sixteen-pad owner; never remove
    // stand contacts or collapse a footprint to admit a second felt head.
    if mallets.first.is_some() && mallets.second.is_some() {
        return Err("paired cymbals admit one four-site felt mallet: twelve mount sites plus two mallets exceed the sixteen-pad limit".into());
    }
    let mallet_sites=4*(usize::from(mallets.first.is_some())+usize::from(mallets.second.is_some()));
    if 12+mallet_sites+mute.map_or(0,compliant_mute::Spec::pad_count)>16 {
        return Err("paired mute and mallet exceed the sixteen-pad limit with all twelve stand contacts retained".into());
    }
    spec.validate()?;
    if let Some(material)=material {
        material.admit_windows(&[upper.band_hz,lower.band_hz])?;
        for side in 0..2 {
            if material.selected(side) && spec.damping[side]!=0. {
                return Err("selected hi-hat shell relaxation requires explicit zero intrinsic damping in the pair input; no double-counted or silently overwritten material loss".into());
            }
        }
    }
    let upper=Shell::new_with_rocking(upper,dt,rocking)?;
    let lower=Shell::new_with_rocking(lower,dt,rocking)?;
    // Conservative reference separation for the ENTIRE finite skins, not just
    // selected rim sites. No ambiguous overlapping initial acoustic bodies.
    let min_z=|s:&Shell|s.skin.triangles().iter().flatten().map(|p|p[2]).fold(f64::INFINITY,f64::min);
    if spec.separation+min_z(&upper)+min_z(&lower)<=0. {
        return Err("paired reference skins fail the separating-plane test".into());
    }
    let hi=1..1+upper.reduction.mode_count();
    let lo=hi.end..hi.end+lower.reduction.mode_count();
    let pedal_coord=lo.end;let second_coord=pedal_coord+1;
    // Jaws follow the unchanged shell/pedal/stick prefix; shaft modes follow
    // jaws. Neither acoustic shell nor existing rigid force addresses move.
    let base=second_coord+usize::from(second.is_some());
    let attachment=mute.map(|mute|{
        let (shell,start)=if mute.lower_shell()? {(&lower,lo.start)}else{(&upper,hi.start)};
        mute.paired_shell(&shell.reduction,&shell.mesh.nodal_thickness_m,&shell.outer_nodes,
            &shell.mesh.mesh.tris,start,base)
    }).transpose()?;
    let mut shaft=shafts.build_with_mallets(base+mute.map_or(0,compliant_mute::Spec::jaw_count),
        second_coord,stroke,second,dt,mallets)?;
    let total=shaft.total;
    if total>fs_couple::render::plate::impact::MAX_IMPACT_MODES {
        return Err("paired cymbals exceed the original complete-state mode ceiling".into());
    }
    let (mut carriage,weight)=ImpactBody::free_mass(spec.carriage[0],0.,0.)?;
    carriage.potential=BodyPotential::Linear(vec![(spec.carriage[1]/spec.carriage[0]).sqrt()]);
    carriage.damping_per_s[0]=spec.carriage[2]/spec.carriage[0];
    if (spec.carriage[1]/spec.carriage[0]).sqrt()*dt>=0.9*std::f64::consts::PI {
        return Err("carriage return frequency exceeds the mechanical Nyquist guard".into());
    }
    let pedal=sticks::Port{coordinate:pedal_coord,weight};
    let mut pads=washers(&upper,&spec.mounts[0],hi.start,total,Some(pedal))?;
    pads.extend(washers(&lower,&spec.mounts[1],lo.start,total,None)?);
    let selected=|target|match target {
        targets::Target::Upper=>(&upper,hi.start),targets::Target::Lower=>(&lower,lo.start),
    };
    let (first_shell,first_start)=selected(targets.first);
    let position=targets.first.local_point(stroke.position_m.unwrap_or(spec.strike));
    let p=first_shell.port(position,ShellFace::Positive)?;
    let (first,stick_weight)=match &mallets.first {
        Some(spec)=>{
            let tip=targets.first.mallet(spec,first_shell,stroke,0,first_start,total)?;
            let tip=spec.bind_shaft(tip,0,&mut shaft)?;
            pads.extend(tip.pads);(tip.body,tip.port.inverse_sqrt_mass)
        }
        None=>match shaft.bodies[0].take() {
            Some(body)=>body,None=>stick_with_speed(stroke.speed_m_s)?,
        },
    };
    let inter=collision(spec,&upper,&lower,lo.start,total)?;
    let mut contacts=Vec::new();
    if mallets.first.is_none() {
        let mut hit=shaft.tip_row(0,stick_weight)?;
        // For either exterior, stick travel is inward and local shell travel
        // is outward: penetration = stick + shell. The lower world's force
        // direction changes through its proper rotation, not a row sign flip.
        for (k,w) in p.weights.iter().enumerate(){hit[first_start+k]=*w;}
        contacts.push(elastic_contact(hit)?);
    }
    contacts.push(inter.clone());
    let second=second.map(|stroke|->Result<_,Error>{
        let (shell,start)=selected(targets.second);
        if let Some(spec)=&mallets.second {
            let tip=targets.second.mallet(spec,shell,stroke,second_coord,start,total)?;
            let tip=spec.bind_shaft(tip,1,&mut shaft)?;
            pads.extend(tip.pads);
            return Ok((tip.body,None,sticks::Port{coordinate:second_coord,weight:tip.port.inverse_sqrt_mass}));
        }
        let stroke=targets.second.local_stroke(stroke);
        let p=shell.port(stroke.position_m.ok_or("second hi-hat stick requires a station")?,ShellFace::Positive)?;
        let (body,weight)=match shaft.bodies[1].take() {
            Some(body)=>body,None=>stick_with_speed(stroke.speed_m_s)?,
        };
        let mut b=shaft.tip_row(1,weight)?;
        for (k,w) in p.weights.iter().enumerate(){b[start+k]=*w;}
        let port=sticks::Port{coordinate:second_coord,weight};
        Ok((body,Some(elastic_contact(b)?),port))
    }).transpose()?;
    let acoustics=if audio {Some(acoustics::Boundary::shell_pair(&upper.skin,hi.start,
        &lower.skin,lo.start,spec.separation)?)}else{None};
    let mut a=vec![0.;total];let mut b=vec![0.;total];
    // Keep the existing upper probe when hand one plays the upper shell. A
    // lower strike need not lie over that shell: use the already admitted
    // first inter-cymbal site for the upper probe instead.
    let p=if targets.first==targets::Target::Upper {p}else{
        upper.port([spec.sites[0].0,spec.sites[0].1],ShellFace::Negative)?
    };
    for (k,w) in p.weights.iter().enumerate(){a[hi.start+k]=-w;}
    let p=lower.port([spec.sites[0].0,-spec.sites[0].1],ShellFace::Negative)?;
    for (k,w) in p.weights.iter().enumerate(){b[lo.start+k]=*w;}
    // Compile before moving either basis; cell/face gaps use both real skins.
    // Stick and pedal columns are zero: pressure loads the shells reciprocally.
    let gas=film.and_then(squeeze::Config::gas);
    let material_gas=film.and_then(squeeze::Config::material_gas).cloned();
    let film=film.map(|f|f.compile(spec,&upper,&lower,hi.start,lo.start,total)).transpose()?;
    let upper_omega=upper.reduction.omegas().to_vec();let lower_omega=lower.reduction.omegas().to_vec();
    let mut up=zero_body(BodyPotential::Shell(upper.reduction),&upper_omega);
    let mut down=zero_body(BodyPotential::Shell(lower.reduction),&lower_omega);
    up.damping_per_s=upper_omega.iter().map(|w|2.*spec.damping[0]*w).collect();
    down.damping_per_s=lower_omega.iter().map(|w|2.*spec.damping[1]*w).collect();
    let mut bodies=vec![first,up,down,carriage];
    let second_stick=second.map(|(body,contact,port)|{
        bodies.push(body);if let Some(contact)=contact {contacts.push(contact);}port
    });
    let mute=attachment.map(|a| {
        let observation=mute.expect("compiled paired mute").observation(a.ports,pads.len());
        bodies.extend(a.bodies);pads.extend(a.pads);observation
    });
    for pad in &mut pads {pad.weights.resize(total,0.);}
    bodies.append(&mut shaft.elastic);
    let flexible_sticks=shaft.ports;
    let system=ImpactSystem::new(bodies,contacts,pads,vec![],config(steps,dt))?;
    let system=match (film,material_gas,gas) {
        (Some(film),Some(source),_)=>system.with_material_compressible_squeeze_film(film,source)?,
        (Some(film),None,Some(gas))=>system.with_compressible_squeeze_film(film,gas)?,
        (Some(film),None,None)=>system.with_squeeze_film(film)?,
        (None,_,_)=>system,
    };
    let system=match material {Some(material)=>material.attach(system)?,None=>system};
    Ok(Pair{experiment:Experiment{flexible_sticks:flexible_sticks.clone(),mute,system:Mechanics::Reference(system),force:vec![0.;total],
        stick_weight,second_stick,observer_a:a,observer_b:b,pressure:None,acoustics,air:None,ambient,head_materials:[None,None]},
        pedal,collision:inter,upper_modes:hi,lower_modes:lo,flexible_sticks,targets})
}

#[cfg(test)]
#[path="hihat_tests.rs"]
mod tests;

#[cfg(test)]
#[path="hihat_flexible_tests.rs"]
mod flexible_tests;

#[cfg(test)]
#[path="hihat_gas_tests.rs"]
mod gas_tests;

#[cfg(test)]
#[path="hihat_thermochemical_tests.rs"]
mod thermochemical_tests;

#[cfg(test)]
#[path="hihat_mallet_tests.rs"]
mod mallet_tests;

#[cfg(test)]
#[path="hihat_target_tests.rs"]
mod target_tests;

#[cfg(test)]
#[path="hihat_rocking_tests.rs"]
mod rocking_tests;

#[cfg(test)]
#[path="hihat_mute_tests.rs"]
mod mute_tests;

#[cfg(test)]
#[path="shell_relaxation_pair_tests.rs"]
mod relaxation_tests;
