use super::*;
use super::super::tests::oscillator;
use fs_phs::RelaxationBranch;
use std::cell::Cell;

fn options() -> BemHarmonicOptions {
    BemHarmonicOptions {
        max_load_cases: 8, max_response_bytes: 1_000_000,
        max_work: 1_000_000, power_relative_tolerance: 1e-10,
    }
}
fn damped() -> BemLoadedModes {
    let mut base = oscillator(0.2);
    base.settings.max_dimension = 12;
    base.with_relaxation_branches(&[RelaxationBranch {
        projection: vec![0.5,-0.3], stiffness: 9_000.0, relaxation_time_s: 0.003,
    }],1,||false).unwrap()
}
fn solve(model: &BemLoadedModes, w: f64, f: &[C64]) -> BemHarmonicResponse {
    model.harmonic_response_batch(w,&[f],&options(),||false).unwrap().remove(0)
}
fn close(a: C64, b: C64, tol: f64) {
    assert!((a-b).abs()<=tol*(1.0+b.abs()),"{a:?} vs {b:?}");
}

#[test]
fn g1_driven_maxwell_response_matches_exact_modulus_and_material_heat() {
    let mut base = oscillator(0.0);
    base.settings.max_dimension=12;
    base.stiffness[3]=60_000.0;
    base.viscosity[3]=0.0;
    let model=base.with_relaxation_branches(&[RelaxationBranch {
        projection:vec![0.0,1.0],stiffness:50_000.0,relaxation_time_s:1.0/300.0,
    }],1,||false).unwrap();
    let f=C64::new(0.7,-0.4);
    for w in [10.0,100.0,295.8,600.0,1200.0] {
        let r=solve(&model,w,&[C64::ZERO,f]);
        let s=C64::new(0.0,-w);
        let factor=s.scale(1.0/300.0)/(C64::ONE+s.scale(1.0/300.0));
        let x=f/(C64::from_re(60_000.0-w*w)+factor.scale(50_000.0));
        let e=factor*x.scale(50_000.0_f64.sqrt());
        close(r.modal_displacement[1],x,1e-13);
        close(r.elastic_memory_sqrt_j[0],e,1e-12);
        let heat=150.0*e.norm_sq();
        assert!((r.power.relaxation_w-heat).abs()<1e-10*heat);
        assert!((r.power.input_w-heat).abs()<1e-10*heat);
        assert!(r.power.kelvin_w.abs()<1e-20 && r.power.radiated_w.abs()<1e-20);
        assert!(r.power.relative_balance<1e-11);
        assert!(r.mechanical_residual<1e-10 && r.memory_residual<1e-10);
    }
}

#[test]
fn g3_driven_surface_is_original_real_bem_and_all_power_channels_close() {
    let model=damped();
    for w in [100.0,299.0,400.0,510.0,700.0] {
        let f=[C64::new(1.0,0.3),C64::new(-0.4,0.7)];
        let r=solve(&model,w,&f);
        let boundary=fs_bem::helmholtz::solve_radiation(model.surface(),
            w/model.medium().sound_speed,model.medium(),&r.surface_normal_velocity_m_s,
            Formulation::BurtonMiller).unwrap();
        let scale=r.surface_pressure_pa.iter().map(|v|v.abs()).fold(0.0,f64::max);
        for (&actual,&expected) in r.surface_pressure_pa.iter().zip(&boundary.pressure) {
            assert!((actual-expected).abs()<2e-10*scale);
        }
        assert!(r.power.input_w>0.0 && r.power.kelvin_w>0.0 && r.power.relaxation_w>0.0);
        assert!(r.power.radiated_w>0.0,"authored breathing mode must radiate");
        assert!((r.power.radiated_w-boundary.radiated_power).abs()<1e-9*r.power.radiated_w);
        let sum=r.power.kelvin_w+r.power.relaxation_w+r.power.radiated_w;
        assert!((r.power.input_w-sum).abs()<1e-9*r.power.input_w);
        // Drop the signed material cross term: the independently computed
        // material response must change, so a diagonal-loss surrogate cannot pass.
        let c=model.relaxation_rows_per_s()[0].as_slice();
        let s=C64::new(0.0,-w);
        let factor=s.scale(0.003)/(C64::ONE+s.scale(0.003));
        let e=factor*(r.modal_displacement[0].scale(c[0])+r.modal_displacement[1].scale(c[1]));
        close(r.elastic_memory_sqrt_j[0],e,1e-12);
    }
}

#[test]
fn g5_batch_reuses_one_assembly_and_preserves_complex_linearity_and_zero_drive() {
    let mut model=damped();
    model.settings.max_operator_evaluations=1;
    let f=[C64::new(1.0,0.3),C64::ZERO];
    let g=[C64::ZERO,C64::new(-0.4,0.7)];
    let h=[f[0]+g[0],f[1]+g[1]];
    let zero=[C64::ZERO;2];
    let batch=model.harmonic_response_batch(350.0,&[&f,&g,&h,&zero],&options(),||false).unwrap();
    for (i,force) in [&f[..],&g[..],&h[..],&zero[..]].into_iter().enumerate() {
        let one=solve(&model,350.0,force);
        assert_eq!(one.modal_displacement,batch[i].modal_displacement);
        assert_eq!(one.surface_pressure_pa,batch[i].surface_pressure_pa);
        assert_eq!(one.elastic_memory_sqrt_j,batch[i].elastic_memory_sqrt_j);
    }
    for i in 0..2 { close(batch[2].modal_displacement[i],batch[0].modal_displacement[i]+batch[1].modal_displacement[i],1e-12); }
    for i in 0..8 { close(batch[2].surface_pressure_pa[i],batch[0].surface_pressure_pa[i]+batch[1].surface_pressure_pa[i],1e-12); }
    let z=&batch[3];
    assert!(z.modal_displacement.iter().chain(&z.surface_pressure_pa).chain(&z.elastic_memory_sqrt_j).all(|v|v.abs()==0.0));
    assert_eq!(z.power.input_w,0.0);
    assert_eq!(z.power.relative_balance,0.0);
}

#[test]
fn g3_physical_drive_scaling_changes_response_and_power_not_loss_law() {
    let model=damped();
    let f=[C64::new(1.0,0.3),C64::new(-0.4,0.7)];
    let reference=solve(&model,350.0,&f);
    for gain in [1e-100,1e100] {
        let result=solve(&model,350.0,&f.map(|z|z.scale(gain)));
        for i in 0..2 { close(divide(result.modal_displacement[i],gain),reference.modal_displacement[i],1e-12); }
        for (actual,expected) in [
            (result.power.input_w,reference.power.input_w),
            (result.power.kelvin_w,reference.power.kelvin_w),
            (result.power.relaxation_w,reference.power.relaxation_w),
            (result.power.radiated_w,reference.power.radiated_w),
        ] { assert!(((actual/gain)/gain-expected).abs()<1e-8*expected.abs()); }
    }
}

#[test]
fn g0_g4_batch_limits_and_late_cancellation_leave_model_reusable() {
    let model=damped();
    let f=[C64::ONE,C64::ZERO];
    let before=solve(&model,350.0,&f);
    let polls=Cell::new(0);
    let reference=model.harmonic_response_batch(350.0,&[&f,&f],&options(),||{
        polls.set(polls.get()+1);false
    }).unwrap();
    let total=polls.get();
    for cut in [1,4,total/2,total-1,total] {
        let count=Cell::new(0);
        let result=model.harmonic_response_batch(350.0,&[&f,&f],&options(),||{
            count.set(count.get()+1);count.get()>=cut
        });
        assert!(result.is_err(),"cancellation at {cut}/{total}");
    }
    let after=solve(&model,350.0,&f);
    assert_eq!(before.modal_displacement,after.modal_displacement);
    assert_eq!(reference[0].surface_pressure_pa,after.surface_pressure_pa);
    for defect in 0..4 {
        let mut settings=options();
        match defect {
            0=>settings.max_load_cases=0,
            1=>settings.max_response_bytes=0,
            2=>settings.max_work=0,
            _=>settings.power_relative_tolerance=f64::NAN,
        }
        assert!(model.harmonic_response_batch(350.0,&[&f],&settings,||false).is_err());
    }
    for w in [0.0,-1.0,f64::NAN,f64::INFINITY] {
        assert!(model.harmonic_response_batch(w,&[&f],&options(),||false).is_err());
    }
    assert!(model.harmonic_response_batch(350.0,&[],&options(),||false).is_err());
    assert!(model.harmonic_response_batch(350.0,&[&f[..1]],&options(),||false).is_err());
}
