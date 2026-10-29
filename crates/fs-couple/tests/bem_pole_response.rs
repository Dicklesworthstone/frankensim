//! Original aluminium shell/BEM/material equations, not a fitted modal response.
#![cfg(not(target_arch = "wasm32"))]
use fs_bem::helmholtz::{Formulation,Medium};
use fs_bem::helmholtz::laplace::complex::ComplexOptions;
use fs_couple::bem_loaded_modes::{BemLoadedModes,BemModePole,BemModesOptions};
use fs_couple::bem_pole_response::{BemParameterDirection,BemPoleResponseOptions};
use fs_la::eigen_complex::lu_complex;
use fs_la::nonlinear_eig::{NonlinearEigError,PoleOptions,newton_pole_with_cancel};
use fs_la::pole_response::PoleResponseOptions;
use fs_math::c64::C64;
use fs_phs::RelaxationBranch;
use fs_plate::shell::reduction::radiation::RadiationSurfaceBudget;
use std::cell::Cell;
#[path="support/material_packs.rs"] mod material_packs;
#[path="bem_loaded_plate_modes/source.rs"] mod source;

fn numeric(w:f64)->PoleOptions{
    PoleOptions{spectral_scale:w,step_tolerance:1e-12,residual_tolerance:1e-12,
        max_iterations:20,max_dimension:512,max_dense_bytes:64*1024*1024,max_work:4_000_000_000}
}
fn analysis(w:f64)->BemPoleResponseOptions{
    BemPoleResponseOptions{numerical:PoleResponseOptions{
        spectral_scale:w,pole_residual_tolerance:1e-11,pole_step_tolerance:1e-11,
        linear_residual_tolerance:1e-11,max_dimension:512,max_vectors:8,
        max_dense_bytes:64*1024*1024,max_work:4_000_000_000},
        physical_residual_tolerance:1e-10,max_operator_bytes:64*1024*1024,max_action_work:1_000_000_000}
}
fn make(source:&source::Plate,medium:Medium,strength:f64,time:f64)->BemLoadedModes{
    let w=source.reduction.omegas()[0];
    let settings=BemModesOptions{formulation:Formulation::PlainCbie,
        boundary:ComplexOptions{max_panels:512,max_bytes:32*1024*1024,max_work:512*512*512},
        panel_refinements:0,reference_angular_frequency_rad_s:w,physical_residual_tolerance:1e-10,
        max_dimension:512,max_dense_bytes:64*1024*1024,max_operator_evaluations:128,max_projection_work:100_000_000};
    BemLoadedModes::from_shell(&source.reduction,&source.thickness,
        RadiationSurfaceBudget{max_panels:512,max_panel_modes:4096},medium,None,settings,||false).unwrap()
        .with_relaxation_branches(&[RelaxationBranch{projection:vec![1.0,-0.3,0.2],
            stiffness:strength,relaxation_time_s:time}],1,||false).unwrap()
}
fn resolve(model:&BemLoadedModes,w:C64,seed:&[C64])->BemModePole{
    model.resonance_pole(w,seed,&numeric(w.abs()),||false).unwrap()
}
fn norm(v:&[C64])->f64{v.iter().map(|z|z.norm_sq()).sum::<f64>().sqrt()}

#[test]
fn g1_original_loaded_response_converges_to_force_residue_with_full_material_history(){
    let plate=source::plate(3,3,3);let w=plate.reduction.omegas()[0];
    let model=make(&plate,source::air(1.0),0.04*w*w,0.8/w);
    let pole=resolve(&model,C64::new(w,-0.01*w),&[C64::ONE,C64::ZERO,C64::ZERO]);
    let prepared=model.prepare_pole_response(&pole,analysis(w),||false).unwrap();
    let force=[C64::new(1.0,0.3),C64::new(-0.1,0.2),C64::new(0.2,-0.4)];
    let residue=prepared.force_residues(&[&force],||false).unwrap().remove(0);
    let mut expected=residue.modal_displacement_residue.clone();
    expected.extend_from_slice(&residue.surface_pressure_residue);
    expected.extend_from_slice(&residue.elastic_memory_residue);
    let mut errors=Vec::new();
    for fraction in [2e-4,1e-4,5e-5]{
        let delta=C64::new(fraction*w,0.4*fraction*w);
        let raw=model.original_pencil(pole.angular_frequency_rad_s+delta,||false).unwrap();
        let mut x=vec![C64::ZERO;model.dimension()];x[..3].copy_from_slice(&force);
        lu_complex(&raw.matrix,model.dimension()).unwrap().solve(&mut x);
        let difference:Vec<_>=x.iter().zip(&expected).map(|(&a,&b)|delta*a-b).collect();
        errors.push(norm(&difference)/norm(&expected));
    }
    assert!(errors[0]>0.0);
    for pair in errors.windows(2){assert!(pair[0]/pair[1]>1.8&&pair[0]/pair[1]<2.2,"{errors:?}");}
    assert!(residue.linear_residual<1e-10);
    assert!(prepared.physical_residuals().iter().all(|&r|r<1e-10));
    for i in 0..model.panel_count(){
        let expected=(0..3).fold(C64::ZERO,|sum,j|sum+
            (C64::new(pole.angular_frequency_rad_s.im,-pole.angular_frequency_rad_s.re)
            *residue.modal_displacement_residue[j]).scale(model.normal_velocity_weights()[j][i]));
        assert!((expected-residue.surface_normal_velocity_residue[i]).abs()<1e-10*(1.0+expected.abs()));
    }
}

#[test]
fn g3_gas_and_material_derivatives_match_independently_rebuilt_loaded_models(){
    let plate=source::plate(3,3,3);let w=plate.reduction.omegas()[0];let air=source::air(1.0);
    let strength=0.04*w*w;let time=0.8/w;let model=make(&plate,air,strength,time);
    let pole=resolve(&model,C64::new(w,-0.01*w),&[C64::ONE,C64::ZERO,C64::ZERO]);
    let prepared=model.prepare_pole_response(&pole,analysis(w),||false).unwrap();
    let dc:Vec<_>=model.relaxation_rows_per_s()[0].iter().map(|c|0.5*c).collect();
    let dt=[time];
    let directions=[
        BemParameterDirection{density_kg_m3:air.density,..Default::default()},
        BemParameterDirection{sound_speed_m_s:air.sound_speed,..Default::default()},
        BemParameterDirection{relaxation_rows_per_s:&dc,..Default::default()},
        BemParameterDirection{relaxation_times_s:&dt,..Default::default()},
    ];
    let derivatives=prepared.parameter_derivatives(&directions,||false).unwrap();
    let h=2.5e-4;
    for (index,actual) in derivatives.iter().enumerate(){
        let mut frequencies=Vec::new();
        for sign in [-1.0,1.0]{
            let scale=1.0+sign*h;let mut medium=air;
            if index==0{medium.density*=scale;}if index==1{medium.sound_speed*=scale;}
            let changed=make(&plate,medium,if index==2{strength*scale}else{strength},
                if index==3{time*scale}else{time});
            frequencies.push(resolve(&changed,pole.angular_frequency_rad_s,&pole.modal_displacement).angular_frequency_rad_s);
        }
        let fd=(frequencies[1]-frequencies[0]).scale(0.5/h);
        let error=(fd-actual.angular_frequency_rad_s).abs();
        assert!(error<2e-5*actual.angular_frequency_rad_s.abs().max(1e-5),"direction={index}: {fd:?} vs {actual:?}");
        assert!(actual.linear_residual<1e-10);
    }
}

#[test]
fn g1_full_mass_stiffness_and_viscosity_directions_match_original_pencil_changes(){
    let plate=source::plate(3,3,3);let w=plate.reduction.omegas()[0];
    let model=make(&plate,source::air(1.0),0.04*w*w,0.8/w);
    let pole=resolve(&model,C64::new(w,-0.01*w),&[C64::ONE,C64::ZERO,C64::ZERO]);
    let prepared=model.prepare_pole_response(&pole,analysis(w),||false).unwrap();
    let mass=[1.0,0.2,0.0,0.2,0.4,0.1,0.0,0.1,0.3];
    let k:Vec<_>=mass.iter().map(|x|0.3*w*w*x).collect();
    let c:Vec<_>=mass.iter().map(|x|0.02*w*x).collect();
    let dirs=[BemParameterDirection{mass:&mass,..Default::default()},
        BemParameterDirection{stiffness_per_s2:&k,..Default::default()},
        BemParameterDirection{viscosity_per_s:&c,..Default::default()}];
    let derivatives=prepared.parameter_derivatives(&dirs,||false).unwrap();
    let n=model.dimension();let h=1e-4;
    for (which,derivative) in derivatives.iter().enumerate(){
        let mut values=Vec::new();
        for sign in [-1.0,1.0]{
            let root=newton_pole_with_cancel(n,pole.angular_frequency_rad_s,prepared.reference_mode(),&numeric(w),
                |omega,cancel|{
                    let mut raw=model.original_pencil(omega,cancel)
                        .map_err(|_|NonlinearEigError::FunctionRefused{at:omega})?;
                    for i in 0..3{for j in 0..3{
                        let (value,slope)=match which{
                            0=>(-(omega*omega).scale(mass[3*i+j]),-omega.scale(2.0*mass[3*i+j])),
                            1=>(C64::from_re(k[3*i+j]),C64::ZERO),
                            _=>(C64::new(omega.im,-omega.re).scale(c[3*i+j]),C64::new(0.0,-c[3*i+j])),
                        };
                        raw.matrix[i*n+j]=raw.matrix[i*n+j]+value.scale(sign*h);
                        raw.derivative[i*n+j]=raw.derivative[i*n+j]+slope.scale(sign*h);
                        raw.residual_scale+=h*value.abs();
                    }}
                    Ok(raw)
                },||false).unwrap();
            values.push(root.value);
        }
        let fd=(values[1]-values[0]).scale(0.5/h);
        assert!((fd-derivative.angular_frequency_rad_s).abs()<2e-6*derivative.angular_frequency_rad_s.abs().max(1e-6));
    }
}

#[test]
fn g0_g4_preparation_refuses_cross_wired_modes_and_batches_replay_after_cancellation(){
    let plate=source::plate(3,3,3);let w=plate.reduction.omegas()[0];
    let model=make(&plate,source::air(1.0),0.04*w*w,0.8/w);
    let pole=resolve(&model,C64::new(w,-0.01*w),&[C64::ONE,C64::ZERO,C64::ZERO]);
    let mut corrupt=pole.clone();corrupt.surface_pressure_pa[0]=corrupt.surface_pressure_pa[0]+C64::from_re(1e6);
    assert!(model.prepare_pole_response(&corrupt,analysis(w),||false).is_err());
    let mut limits=analysis(w);limits.max_operator_bytes=1;
    assert!(model.prepare_pole_response(&pole,limits,||false).is_err());
    let prepared=model.prepare_pole_response(&pole,analysis(w),||false).unwrap();
    let force=[C64::ONE,C64::ZERO,C64::ZERO];let count=Cell::new(0);
    let before=prepared.force_residues(&[&force,&force],||{count.set(count.get()+1);false}).unwrap();
    let total=count.get();
    for stop in [1,total/2,total]{
        let calls=Cell::new(0);
        assert!(prepared.force_residues(&[&force,&force],||{calls.set(calls.get()+1);calls.get()==stop}).is_err());
    }
    let after=prepared.force_residues(&[&force],||false).unwrap();
    assert_eq!(before[0].surface_pressure_residue,after[0].surface_pressure_residue);
    assert!(prepared.parameter_derivatives(&[BemParameterDirection{mass:&[1.0],..Default::default()}],||false).is_err());
    assert!(prepared.parameter_derivatives(&[BemParameterDirection{density_kg_m3:f64::NAN,..Default::default()}],||false).is_err());
    assert!(prepared.force_residues(&[],||false).is_err());
}
