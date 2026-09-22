//! Distributed enclosed air from the same drum dimensions and actual head modes.
use super::{Error,ImpactBody,ImpactSystem,ModePair,Obstacle,TensionedDisk,config};
use fs_couple::render::plate::impact::cavity::{CavityCoupling,
    cylinder::{CylinderSpec,CylindricalCavity,SidewallAperture},neck::CavityNeck};
use fs_couple::vibroacoustic::{AcousticMedium,StructuralModes,assemble_coupling};
use fs_exec::CancelGate;
use fs_couple::render::plate::impact::linear::{LinearImpactConfig,LinearImpactSystem};

pub fn option(args:&mut Vec<String>)->Result<bool,Error> {
    let count=args.iter().filter(|a|a.as_str()=="--cavity-modes").count();
    if count>1 {return Err("--cavity-modes may be supplied only once".into());}
    args.retain(|a|a!="--cavity-modes");Ok(count==1)
}

/// Both existing mechanical images can now consume the same enclosed air.
pub fn admit_command(distributed:bool,command:&str)->Result<(),Error> {
    if distributed && !matches!(command,"drum"|"drum-wav"|"drum-mic"|
        "drum-stretch"|"drum-stretch-wav"|"drum-stretch-mic"|
        "drum-modal"|"drum-modal-wav"|"drum-modal-mic"|
        "snare"|"snare-wav"|"snare-mic"|"snare-off"|"snare-off-wav"|"snare-off-mic") {
        return Err("--cavity-modes requires a drum or snare command, including WAV/microphone variants".into());
    }
    Ok(())
}

/// All neck physics and position are explicit; this is not a measured drum card.
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct NeckOptions {
    pub radius_m:f64,
    pub effective_length_m:f64,
    pub resistance_pa_s_m3:f64,
    pub azimuth_rad:f64,
    pub axial_position_m:f64,
}

pub fn neck_option(args:&mut Vec<String>)->Result<Option<NeckOptions>,Error> {
    let mut positions=args.iter().enumerate().filter(|(_,v)|v.as_str()=="--cavity-neck");
    let Some((start,_))=positions.next() else {return Ok(None);};
    if positions.next().is_some() || args.len()-start<6 {
        return Err("--cavity-neck needs one radius_m length_eff_m resistance_Pa_s_m3 azimuth_rad z_m declaration".into());
    }
    let mut values=[0.0_f64;5];
    for (value,text) in values.iter_mut().zip(&args[start+1..start+6]) {*value=text.parse()?;}
    if values.iter().any(|v|!v.is_finite()) || values[0]<=0.0 || values[1]<=0.0 || values[2]<0.0 {
        return Err("neck needs positive finite radius/length, passive resistance and finite wall coordinates".into());
    }
    let neck=NeckOptions {radius_m:values[0],effective_length_m:values[1],
        resistance_pa_s_m3:values[2],azimuth_rad:values[3],axial_position_m:values[4]};
    args.drain(start..start+6);
    Ok(Some(neck))
}

pub fn admit_neck_command(neck:Option<NeckOptions>,distributed:bool,command:&str)->Result<(),Error> {
    if neck.is_some() && (!distributed || !matches!(command,"drum"|"drum-stretch")) {
        return Err("--cavity-neck requires --cavity-modes and drum/drum-stretch CSV; vented exterior radiation is not implemented, so WAV/microphone export is refused".into());
    }
    Ok(())
}

/// Cached interior observations, not an exterior microphone or extra force.
pub struct InteriorPressure {pub coupling:CavityCoupling,first:Vec<f64>,second:Vec<f64>}
impl InteriorPressure {
    pub fn uniform_pressure(&self,state:&[f64])->Result<f64,Error> {
        let mut pressure=[0.0;8];
        self.coupling.pressures_into(state,&mut pressure[..self.coupling.cavity_modes()])?;
        // The cylindrical basis owns an exact constant first member. Include
        // displaced neck volume here, not only the old solid-head contraction.
        Ok(pressure[0])
    }
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
    contacts:Vec<Obstacle>,radius:f64,depth:f64,steps:u64,dt_s:f64,neck:Option<NeckOptions>)
    ->Result<(ImpactSystem,InteriorPressure),Error> {
    let gate=CancelGate::new_clock_free();
    // Both head ranges remain the original prefix. Include every appended
    // striker before allocating cavity inertia; its coupling row stays zero.
    let structural=bodies.iter().try_fold(0usize,|n,b|n.checked_add(b.initial.len()))
        .ok_or("cavity body-count overflow")?;
    let InteriorPressure {coupling,first,second}=compile(films,modes,radius,depth,structural,
        fs_couple::render::plate::impact::MAX_IMPACT_MODES,neck,&gate)?;
    let (system,coupling)=coupling.build(bodies,contacts,vec![],config(steps,dt_s),&gate)?;
    Ok((system,InteriorPressure {coupling,first,second}))
}

/// Same cavity basis and surface integrals, with wires retained after the heads.
/// No direct air coupling for the striker or filaments; wire reactions act on
/// the resonant head. Appended inertia follows EVERY original solid coordinate.
#[allow(clippy::too_many_arguments)]
pub fn build_prepared(films:&[TensionedDisk],modes:&[Vec<ModePair>],bodies:Vec<ImpactBody>,
    contacts:Vec<Obstacle>,radius:f64,depth:f64,configuration:LinearImpactConfig)
    ->Result<(LinearImpactSystem,InteriorPressure),Error> {
    let gate=CancelGate::new_clock_free();
    let structural=bodies.iter().try_fold(0usize,|n,b|n.checked_add(b.initial.len()))
        .ok_or("prepared cavity body-count overflow")?;
    let InteriorPressure {coupling,first,second}=compile(films,modes,radius,depth,structural,
        configuration.coupling.max_modes,None,&gate)?;
    let (system,coupling)=coupling.build_linear(bodies,contacts,
        core::f64::consts::PI*radius*radius,configuration,&gate)?;
    eprintln!("mechanical image: prepared modal heads/wires plus simultaneous distributed-air/contact reactions; no wire homogenization or direct gas audio; real-time performance unqualified");
    Ok((system,InteriorPressure {coupling,first,second}))
}

#[allow(clippy::too_many_arguments)]
fn compile(films:&[TensionedDisk],modes:&[Vec<ModePair>],radius:f64,depth:f64,
    structural:usize,maximum_modes:usize,neck:Option<NeckOptions>,gate:&CancelGate)
    ->Result<InteriorPressure,Error> {
    let heads=1+modes.iter().map(Vec::len).sum::<usize>();
    if structural<heads || structural>maximum_modes || maximum_modes>4096 {
        return Err("cavity needs the original head prefix and a bounded structural layout".into());
    }
    let air=basis(radius,depth,gate)?;
    let head_coupling=interface(films,modes,&air,gate)?;
    let mut coupling=vec![0.0;structural*air.modes().len()];
    coupling[..head_coupling.len()].copy_from_slice(&head_coupling);
    let sampled=air.sample(&[[0.0,0.0,0.0]],8)?;
    eprintln!("distributed cavity: R={radius}m, depth={depth}m; radial_intervals=32, acoustic_hz={:?}; explicit zero acoustic drag; rigid cylindrical sidewall, not calibrated losses",
        sampled.omegas.iter().map(|w|w/core::f64::consts::TAU).collect::<Vec<_>>());
    let mut compiled=CavityCoupling::new_with_mode_budget(&sampled,structural,&coupling,
        &vec![0.0;air.modes().len()],maximum_modes)?;
    if let Some(neck)=neck {
        let averages=air.sidewall_averages(SidewallAperture {radius_m:neck.radius_m,
            azimuth_rad:neck.azimuth_rad,axial_position_m:neck.axial_position_m,
            radial_rings:8,angular_points:32,maximum_terms:2048},gate)?;
        compiled=compiled.with_necks(vec![CavityNeck {
            area_m2:core::f64::consts::PI*neck.radius_m*neck.radius_m,
            effective_length_m:neck.effective_length_m,resistance_pa_s_m3:neck.resistance_pa_s_m3,
            pressure_shape_averages:averages,initial_volume_m3:0.0,initial_flow_m3_s:0.0,
        }],gate)?;
        eprintln!("compact neck: radius={}m, effective_length={}m, resistance={}Pa*s/m^3, azimuth={}rad, z={}m; 8x32 wall-area quadrature; zero-gauge reservoir; parameters are declarations, not calibrated losses; exterior audio refused",
            neck.radius_m,neck.effective_length_m,neck.resistance_pa_s_m3,neck.azimuth_rad,neck.axial_position_m);
    }
    let first=air.values_at([0.4*radius,0.2*radius,0.0])?;
    let second=air.values_at([-0.4*radius,-0.2*radius,depth])?;
    Ok(InteriorPressure {coupling:compiled,first,second})
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn neck_arguments_preserve_playing_controls_and_refuse_unimplemented_audio() {
        let mut args:Vec<String>=["drum-stretch","128","--cavity-neck","0.005","0.012","1000","0.4","0.08",
            "--cavity-modes","--strike-speed-m-s","4","--prepared-nonlinear"].map(String::from).to_vec();
        let neck=neck_option(&mut args).unwrap().unwrap();
        assert_eq!(neck.radius_m,0.005);assert_eq!(neck.effective_length_m,0.012);
        assert_eq!(neck.resistance_pa_s_m3,1000.0);assert_eq!(neck.azimuth_rad,0.4);assert_eq!(neck.axial_position_m,0.08);
        assert!(option(&mut args).unwrap());
        assert!(super::super::mechanics::prepared_option(&mut args).unwrap());
        let (positional,stroke)=super::super::playing::parse(args).unwrap();
        assert_eq!(positional,["drum-stretch","128"]);assert_eq!(stroke.speed_m_s,4.0);
        assert!(admit_neck_command(Some(neck),true,"drum-stretch").is_ok());
        for command in ["drum-mic","drum-wav","drum-stretch-mic","drum-stretch-wav","snare","splash"] {
            assert!(admit_neck_command(Some(neck),true,command).is_err());
        }
        assert!(admit_neck_command(Some(neck),false,"drum").is_err());
        assert!(admit_neck_command(None,false,"splash").is_ok());
        for declaration in ["--cavity-neck", "--cavity-neck 0.005 0.012 -1 0 0.08",
            "--cavity-neck NaN 0.012 1000 0 0.08", "--cavity-neck 0 0.012 1000 0 0.08",
            "--cavity-neck 0.005 0.012 1000 0 0.08 --cavity-neck 0.005 0.012 1000 0 0.08"] {
            let mut args:Vec<String>=declaration.split_whitespace().map(String::from).collect();
            let original=args.clone();assert!(neck_option(&mut args).is_err());assert_eq!(args,original);
        }
    }
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
        let mut compact=super::super::drum_with_air(128,2e-6,false,false,None,true,stroke,false,None).unwrap();
        let mut distributed=super::super::drum_with_air(128,2e-6,false,false,None,true,stroke,true,None).unwrap();
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

    #[test]
    fn finite_neck_changes_struck_stretching_heads_and_pressure_contains_outward_slug_volume() {
        let stroke=super::super::Stroke {speed_m_s:4.0,position_m:Some([0.06,0.01])};
        let neck=NeckOptions {radius_m:0.005,effective_length_m:0.012,resistance_pa_s_m3:1000.0,
            azimuth_rad:0.4,axial_position_m:0.08};
        let mut sealed=super::super::drum_with_air(256,2e-6,false,false,None,true,stroke,true,None).unwrap();
        let mut vented=super::super::drum_with_air(256,2e-6,false,false,None,true,stroke,true,Some(neck)).unwrap();
        let original=sealed.force.len();
        assert_eq!(vented.force.len(),original+1);
        assert_eq!(sealed.system.state(),&vented.system.state()[..2*original]);
        assert_eq!(vented.observer_a[original],0.0);assert_eq!(vented.observer_b[original],0.0);
        assert_eq!(vented.force[original],0.0);assert!(vented.acoustics.is_none());
        sealed.system=sealed.system.into_prepared_nonlinear().unwrap();
        vented.system=vented.system.into_prepared_nonlinear().unwrap();
        let gate=CancelGate::new_clock_free();let mut changed=0.0_f64;let mut flow=0.0_f64;let mut loss=0.0_f64;
        for _ in 0..256 {
            sealed.system.step(&sealed.force,&gate).unwrap();
            let f=vented.system.step(&vented.force,&gate).unwrap();
            let x=vented.system.state();let air=vented.air.as_ref().unwrap();
            let n=air.coupling.neck_observation(x,0).unwrap();
            // The old volume column knows only head displacement and would
            // falsely report sealed compression. This checks the CSV observable.
            let volume=vented.pressure.as_ref().unwrap();
            let heads=super::super::cavity_pressure(volume,x);
            let expected=heads-volume.bulk_modulus_pa/volume.volume_m3*n.displaced_volume_m3;
            assert!((air.uniform_pressure(x).unwrap()-expected).abs()<1e-7*(1.0+expected.abs()));
            let solid=air.coupling.structural_modes();
            changed=changed.max(sealed.system.state()[..2*solid].iter().zip(&x[..2*solid])
                .map(|(a,b)|(a-b).abs()).fold(0.0_f64,f64::max));
            flow=flow.max(n.volume_flow_m3_s.abs());loss+=n.dissipated_power_w*2e-6;
            assert!(f.balance_residual_j.abs()<1e-7);
        }
        assert!(changed>1e-12);assert!(flow>1e-10);assert!(loss>0.0);
        assert!(vented.system.membrane_observation(1).unwrap().stretching_energy_j>0.0);
        // Reject before any geometry or exterior transfer is constructed.
        assert!(super::super::drum_with_air(1,2e-6,true,false,None,true,stroke,true,Some(neck)).is_err());
    }
}

#[cfg(test)]
#[path = "cavity_prepared_tests.rs"]
mod prepared_tests;
