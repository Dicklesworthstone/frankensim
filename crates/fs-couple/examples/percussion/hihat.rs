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

// One independently prepared physical shell and its exact inner-skin chart.
struct Shell {
    mesh:MeshShell,
    reduction:ShellReduction,
    skin:ShellRadiationSurface,
    inner_nodes:Vec<[f64;3]>,
}
impl Shell {
    fn new(spec:&specimen::Specimen,dt:f64)->Result<Self,Error> {
        let (mesh,reduction)=shell_prepare::prepare(spec,dt)?;
        let skin=reduction.radiation_surface(&mesh.nodal_thickness_m,
            RadiationSurfaceBudget{max_panels:2048,max_panel_modes:65536})?;
        // The negative skin reverses vertex order; restore SOURCE connectivity
        // before locating contact barycentrics. Do not project onto midsurface XY.
        let mut inner_nodes=vec![[0.;3];mesh.mesh.nodes.len()];
        for (f,t) in mesh.mesh.tris.iter().enumerate() {
            for a in 0..3 {inner_nodes[t[a]]=skin.triangles()[2*f+1][2-a];}
        }
        Ok(Self{mesh,reduction,skin,inner_nodes})
    }
    fn port(&self,p:[f64;2],face:ShellFace)->Result<fs_plate::shell::reduction::radiation::ShellSurfacePort,Error> {
        let (f,b)=playing::shell_location(&self.inner_nodes,&self.mesh.mesh.tris,p)?;
        Ok(self.reduction.surface_point_port(&self.mesh.nodal_thickness_m,f,b,face,[0.,0.,1.])?)
    }
}
struct Pair {
    experiment:Experiment,
    pedal:sticks::Port,
    collision:Obstacle,
    upper_modes:std::ops::Range<usize>,
    lower_modes:std::ops::Range<usize>,
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
    spec.validate()?;
    let upper=Shell::new(upper,dt)?;let lower=Shell::new(lower,dt)?;
    // Conservative reference separation for the ENTIRE finite skins, not just
    // selected rim sites. No ambiguous overlapping initial acoustic bodies.
    let min_z=|s:&Shell|s.skin.triangles().iter().flatten().map(|p|p[2]).fold(f64::INFINITY,f64::min);
    if spec.separation+min_z(&upper)+min_z(&lower)<=0. {
        return Err("paired reference skins fail the separating-plane test".into());
    }
    let hi=1..1+upper.reduction.mode_count();
    let lo=hi.end..hi.end+lower.reduction.mode_count();
    let pedal_coord=lo.end;let second_coord=pedal_coord+1;
    let total=second_coord+usize::from(second.is_some());
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
    let position=stroke.position_m.unwrap_or(spec.strike);
    let p=upper.port(position,ShellFace::Positive)?;
    let (stick,stick_weight)=stick_with_speed(stroke.speed_m_s)?;
    let mut hit=vec![0.;total];hit[0]=stick_weight;
    for (k,w) in p.weights.iter().enumerate(){hit[hi.start+k]=*w;}
    let inter=collision(spec,&upper,&lower,lo.start,total)?;
    let mut contacts=vec![elastic_contact(hit)?,inter.clone()];
    let second=second.map(|stroke|->Result<_,Error>{
        let p=upper.port(stroke.position_m.ok_or("second hi-hat stick requires a station")?,ShellFace::Positive)?;
        let (body,weight)=stick_with_speed(stroke.speed_m_s)?;
        let mut b=vec![0.;total];b[second_coord]=weight;
        for (k,w) in p.weights.iter().enumerate(){b[hi.start+k]=*w;}
        Ok((body,elastic_contact(b)?,sticks::Port{coordinate:second_coord,weight}))
    }).transpose()?;
    let acoustics=if audio {Some(acoustics::Boundary::shell_pair(&upper.skin,hi.start,
        &lower.skin,lo.start,spec.separation)?)}else{None};
    let mut a=vec![0.;total];let mut b=vec![0.;total];
    for (k,w) in p.weights.iter().enumerate(){a[hi.start+k]=-w;}
    let p=lower.port([spec.sites[0].0,-spec.sites[0].1],ShellFace::Negative)?;
    for (k,w) in p.weights.iter().enumerate(){b[lo.start+k]=*w;}
    let upper_omega=upper.reduction.omegas().to_vec();let lower_omega=lower.reduction.omegas().to_vec();
    let mut up=zero_body(BodyPotential::Shell(upper.reduction),&upper_omega);
    let mut down=zero_body(BodyPotential::Shell(lower.reduction),&lower_omega);
    up.damping_per_s=upper_omega.iter().map(|w|2.*spec.damping[0]*w).collect();
    down.damping_per_s=lower_omega.iter().map(|w|2.*spec.damping[1]*w).collect();
    let mut bodies=vec![stick,up,down,carriage];
    let second_stick=second.map(|(body,contact,port)|{bodies.push(body);contacts.push(contact);port});
    let system=ImpactSystem::new(bodies,contacts,pads,vec![],config(steps,dt))?;
    Ok(Pair{experiment:Experiment{mute:None,system:Mechanics::Reference(system),force:vec![0.;total],
        stick_weight,second_stick,observer_a:a,observer_b:b,pressure:None,acoustics,air:None},
        pedal,collision:inter,upper_modes:hi,lower_modes:lo})
}

#[cfg(test)]
#[path="hihat_tests.rs"]
mod tests;
