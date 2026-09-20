//! Distributed enclosed air from the same drum dimensions and actual head modes.
use super::{Error,ImpactBody,ImpactSystem,ModePair,Obstacle,TensionedDisk,config};
use fs_couple::render::plate::impact::cavity::{CavityCoupling,cylinder::{CylinderSpec,CylindricalCavity}};
use fs_couple::vibroacoustic::{AcousticMedium,StructuralModes,assemble_coupling};
use fs_exec::CancelGate;

pub fn option(args:&mut Vec<String>)->Result<bool,Error> {
    let count=args.iter().filter(|a|a.as_str()=="--cavity-modes").count();
    if count>1 {return Err("--cavity-modes may be supplied only once".into());}
    args.retain(|a|a!="--cavity-modes");Ok(count==1)
}

/// Cached interior observations, not an exterior microphone or extra force.
pub struct InteriorPressure {pub coupling:CavityCoupling,first:Vec<f64>,second:Vec<f64>}
impl InteriorPressure {
    pub fn points(&self,state:&[f64])->Result<(f64,f64),Error> {
        Ok((self.coupling.pressure_at(state,&self.first)?,
            self.coupling.pressure_at(state,&self.second)?))
    }
}

fn basis(radius:f64,depth:f64,gate:&CancelGate)->Result<CylindricalCavity,Error> {
    Ok(CylindricalCavity::new(CylinderSpec {radius_m:radius,depth_m:depth,
        radial_intervals:32,maximum_azimuthal_order:4,maximum_axial_order:2,
        maximum_frequency_hz:1100.0,maximum_modes:8,eigen_residual_tolerance:1e-7},
        AcousticMedium {rho0:1.2,c0:343.0},gate)?)
}

/// Same piecewise-linear head displacement, sampled at three triangle points.
/// Positive normal motion is OUTWARD: downward top motion has negative sign;
/// downward bottom motion has positive sign. Cavity and solid shapes use the
/// exact same quadrature points and positive area measures.
fn interface(films:&[TensionedDisk],modes:&[Vec<ModePair>],air:&CylindricalCavity,
    gate:&CancelGate)->Result<Vec<f64>,Error> {
    if films.len()!=2 || modes.len()!=2 {return Err("cavity needs the two original heads".into());}
    let total=1+modes.iter().map(Vec::len).sum::<usize>();
    let points=films.iter().map(|f|f.mesh.tris.len()*3).sum::<usize>();
    if points>12000 || total>64 {return Err("cavity/head quadrature budget exceeded".into());}
    let mut positions=Vec::with_capacity(points);let mut areas=Vec::with_capacity(points);
    let mut shapes=vec![vec![0.0;points];total];let mut offset=1;let mut point=0;
    for (head,(film,pairs)) in films.iter().zip(modes).enumerate() {
        let sign=if head==0 {-1.0}else{1.0};let z=head as f64*air.spec().depth_m;
        for triangle in &film.mesh.tris {
            if gate.is_requested() {return Err("cavity interface assembly cancelled".into());}
            let nodes=triangle.map(|i|film.mesh.nodes[i]);
            let area=((nodes[1].0-nodes[0].0)*(nodes[2].1-nodes[0].1)
                -(nodes[2].0-nodes[0].0)*(nodes[1].1-nodes[0].1)).abs()/2.0;
            if !area.is_finite() || area<=0.0 {return Err("cavity interface has an invalid triangle".into());}
            for bary in [[2.0/3.0,1.0/6.0,1.0/6.0],[1.0/6.0,2.0/3.0,1.0/6.0],[1.0/6.0,1.0/6.0,2.0/3.0]] {
                positions.push([nodes.iter().zip(bary).map(|(n,b)|n.0*b).sum(),
                    nodes.iter().zip(bary).map(|(n,b)|n.1*b).sum(),z]);areas.push(area/3.0);
                for (i,mode) in pairs.iter().enumerate() {
                    shapes[offset+i][point]=sign*triangle.iter().zip(bary).map(|(&node,b)|
                        b*film.model.dof_map[3*node].map_or(0.0,|d|mode.phi[d])).sum::<f64>();
                }
                point+=1;
            }
        }
        offset+=pairs.len();
    }
    let sampled=air.sample(&positions,100000)?;
    let structure=StructuralModes {omegas:std::iter::once(0.0).chain(modes.iter().flatten().map(|m|m.lambda.sqrt())).collect(),
        shapes,loss_factor:0.0};
    Ok(assemble_coupling(&structure,&sampled,&areas)?)
}

/// Compile distributed air IN PLACE OF the uniform spring. The zero pressure
/// mode already supplies the full bulk compliance; adding the old spring again
/// would double it. Acoustic state only reaches exterior sound through head motion.
#[allow(clippy::too_many_arguments)]
pub fn build(films:&[TensionedDisk],modes:&[Vec<ModePair>],bodies:Vec<ImpactBody>,
    contacts:Vec<Obstacle>,radius:f64,depth:f64,steps:u64,dt_s:f64)
    ->Result<(ImpactSystem,InteriorPressure),Error> {
    let gate=CancelGate::new_clock_free();let air=basis(radius,depth,&gate)?;
    let coupling=interface(films,modes,&air,&gate)?;
    let structural=1+modes.iter().map(Vec::len).sum::<usize>();
    let sampled=air.sample(&[[0.0,0.0,0.0]],8)?;
    eprintln!("distributed cavity: R={radius}m, depth={depth}m; radial_intervals=32, acoustic_hz={:?}; explicit zero acoustic drag; rigid cylindrical sidewall, not calibrated losses",
        sampled.omegas.iter().map(|w|w/core::f64::consts::TAU).collect::<Vec<_>>());
    let compiled=CavityCoupling::new(&sampled,structural,&coupling,&vec![0.0;air.modes().len()])?;
    let (system,compiled)=compiled.build(bodies,contacts,vec![],config(steps,dt_s),&gate)?;
    let first=air.values_at([0.4*radius,0.2*radius,0.0])?;
    let second=air.values_at([-0.4*radius,-0.2*radius,depth])?;
    Ok((system,InteriorPressure {coupling:compiled,first,second}))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn cavity_flag_preserves_playing_controls_and_rejects_duplicates() {
        let mut args=vec!["drum-mic".into(),"--cavity-modes".into(),"128".into(),
            "--prepared-nonlinear".into(),"--strike-speed-m-s".into(),"4".into()];
        assert!(option(&mut args).unwrap());
        assert!(super::super::mechanics::prepared_option(&mut args).unwrap());
        let (positional,stroke)=super::super::playing::parse(args).unwrap();
        assert_eq!(positional,["drum-mic","128"]);assert_eq!(stroke.speed_m_s,4.0);
        assert!(option(&mut vec!["--cavity-modes".into();2]).is_err());
    }
    #[test]
    fn actual_head_quadrature_preserves_uniform_area_and_both_outward_signs() {
        let gate=CancelGate::new_clock_free();let radius=0.1703;let depth=0.1651;
        let make=||TensionedDisk::new(super::super::TensionedDiskSpec {radius_m:radius,thickness_m:0.000254,
            young_pa:4e9,poisson:0.38,density_kg_m3:1390.0,tension_n_m:3000.0,
            radial_intervals:2,azimuths:16},super::super::mesh_budget()).unwrap();
        let films=vec![make(),make()];let mut modes=Vec::new();
        for film in &films {
            let n=film.model.free;let mut k=vec![0.0;n*n];let mut m=k.clone();
            for i in 0..n {for j in 0..n {k[i*n+j]=film.model.k.get(i,j);m[i*n+j]=film.model.m.get(i,j);}}
            let mut pairs=fs_modal::eigh_gen_dense(&k,&m,n).unwrap();pairs.truncate(2);modes.push(pairs);
        }
        let air=basis(radius,depth,&gate).unwrap();let c=interface(&films,&modes,&air,&gate).unwrap();
        let count=air.modes().len();let axial=air.modes().iter().position(|m|m.azimuthal_order==0 && m.axial_order==1).unwrap();
        let mut offset=1;
        for head in 0..2 {for mode in &modes[head] {
            let area=films[head].modal_area(&mode.phi).unwrap();
            let expected=area*if head==0 {-1.0}else{1.0};
            assert!((c[offset*count]-expected).abs()<1e-12);
            assert!((c[offset*count+axial]+area).abs()<1e-12);
            offset+=1;
        }}
        assert!(c[..count].iter().all(|v|*v==0.0),"stick does not displace enclosed air directly");
    }
    #[test]
    fn distributed_air_changes_real_struck_head_motion_without_direct_acoustic_drive() {
        let stroke=super::super::Stroke {speed_m_s:4.0,position_m:Some([0.06,0.01])};
        let mut compact=super::super::drum_with_air(128,2e-6,false,false,None,true,stroke,false).unwrap();
        let mut distributed=super::super::drum_with_air(128,2e-6,false,false,None,true,stroke,true).unwrap();
        let solid=compact.force.len();
        assert_eq!(compact.system.state(),&distributed.system.state()[..2*solid]);
        assert_eq!(&compact.observer_a[..],&distributed.observer_a[..solid]);
        assert!(distributed.observer_a[solid..].iter().all(|v|*v==0.0));
        assert!(distributed.force[solid..].iter().all(|v|*v==0.0));
        assert_eq!(distributed.air.as_ref().unwrap().coupling.cavity_modes(),6);
        compact.system=compact.system.into_prepared_nonlinear().unwrap();
        distributed.system=distributed.system.into_prepared_nonlinear().unwrap();
        let gate=CancelGate::new_clock_free();let mut changed=0.0_f64;let mut nonuniform=0.0_f64;
        for _ in 0..128 {
            compact.system.step(&compact.force,&gate).unwrap();
            let frame=distributed.system.step(&distributed.force,&gate).unwrap();
            for (&a,&b) in compact.system.state().iter().zip(distributed.system.state()) {changed=changed.max((a-b).abs());}
            let (a,b)=distributed.air.as_ref().unwrap().points(distributed.system.state()).unwrap();
            nonuniform=nonuniform.max((a-b).abs());assert!(frame.balance_residual_j.abs()<1e-7);
        }
        assert!(changed>1e-12);assert!(nonuniform>1e-5);
        assert!(distributed.system.membrane_observation(1).unwrap().stretching_energy_j>0.0);
    }
}
