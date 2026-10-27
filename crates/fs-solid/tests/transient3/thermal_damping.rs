//! Changing material stiffness shares one passive damping/motion/energy path.
use super::*;
use fs_blake3::ContentHash;

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a-b).abs() <= tolerance, "{a} != {b}, tolerance {tolerance}");
}

#[test]
fn material_softening_and_stiffening_keep_damping_passive_and_acceleration_physical() {
    with_cx(false, |cx| {
        let fixed: Vec<_> = (0..12).filter(|d| *d != 3).collect();
        let frame = ContentHash([74;32]);
        let original = TetElasticMaterial::try_new(6.0,60.0,0.0,[73;32]).unwrap();
        let old_problem = problem(&original, &fixed);
        let geometry = geometry(cx);
        let f0 = 0.0002;
        let f1 = -0.0001;
        let load0 = geometry.assemble(&[],VolumeLoad3::None,&[(3,f0)],cx).unwrap();
        let load1 = geometry.assemble(&[],VolumeLoad3::None,&[(3,f1)],cx).unwrap();
        let cfg = TransientConfig3 { mass_damping_per_s:0.4, stiffness_damping_s:0.1, ..config() };
        let (q0,v0,eps0,eps1,h,m,k0) = (0.002,-0.001,0.0003,0.0012,0.02,0.1,10.0);
        let mut q = ZERO; q[1][0] = q0;
        let mut v = ZERO; v[1][0] = v0;
        let old_strain = TetThermalStrainState::try_new(300.0,290.0,
            [eps0,0.0,0.0,0.0,0.0,0.0],original.material_state_identity,frame,0.01).unwrap();
        let accepted = load0.begin_transient(&old_problem,
            Some(TetThermalStrainField::Uniform(&old_strain)),0.0,&q,&v,cfg,cx).unwrap();
        let saved = accepted.state().clone();
        for young in [6.0,36.0,300.0] {
            let replacement = TetElasticMaterial::try_new(6.0,young,0.0,[75;32]).unwrap();
            let next_problem = problem(&replacement,&fixed);
            let next_strain = TetThermalStrainState::try_new(301.0,290.0,
                [eps1,0.0,0.0,0.0,0.0,0.0],replacement.material_state_identity,frame,0.01).unwrap();
            let (mut next,report) = accepted.advance_thermal_candidate(&next_problem,
                TetThermalStrainField::Uniform(&next_strain),&load1,h,frame,0.01,cx).unwrap();
            let k1 = young/6.0;
            let c0 = cfg.mass_damping_per_s*m + cfg.stiffness_damping_s*k0;
            let c1 = cfg.mass_damping_per_s*m + cfg.stiffness_damping_s*k1;
            let cb = (c0+c1)/2.0;
            let force_mid = (f0+k0*eps0+f1+k1*eps1)/2.0;
            // Solve the independent displacement/velocity two-equation system.
            let expected_v = ((m/h-cb/2.0)*v0-k0*q0/2.0+force_mid
                -k1*(q0+h*v0/2.0)/2.0)/(m/h+cb/2.0+h*k1/4.0);
            let expected_q = q0+h*(v0+expected_v)/2.0;
            let expected_a = (f1+k1*(eps1-expected_q)-c1*expected_v)/m;
            close(next.state().displacement[0],expected_q,1e-14);
            close(next.state().velocity[0],expected_v,1e-13);
            close(next.state().acceleration[0],expected_a,1e-11);
            close(report.native.damping_loss_j,h*cb*((v0+expected_v)/2.0).powi(2),1e-14);
            assert!(report.native.damping_loss_j >= 0.0);
            close(report.native.stiffness_parameter_work_j,
                0.5*q0*(k1-k0)*expected_q,1e-13);
            close(report.stiffness_material_work_j,
                0.5*(q0-eps0)*(k1-k0)*(expected_q-eps1),1e-13);
            close(report.kinetic_plus_elastic_after_j,
                0.5*m*expected_v*expected_v+0.5*k1*(expected_q-eps1).powi(2),1e-13);
            assert!(report.native.balance_error_j.abs() < 1e-12);
            assert!(report.physical_balance_error_j.abs() < 1e-12);
            assert!(report.native.iterations <= cfg.solve.maximum_iterations);
            close(next.recover_response(frame,0.01,cx).unwrap().stress.elements()[0]
                .stress_mandel_pa()[0],young*(expected_q-eps1),1e-11);
            // Return to the existing fixed-operator Newmark path without reset.
            let qp = expected_q+h*expected_v+h*h*expected_a/4.0;
            let vp = expected_v+h*expected_a/2.0;
            let a2 = (f1+k1*(eps1-qp)-c1*vp)/(m+h*c1/2.0+h*h*k1/4.0);
            let second = next.advance(&load1,h,cx).unwrap();
            close(next.state().displacement[0],qp+h*h*a2/4.0,1e-13);
            close(next.state().velocity[0],vp+h*a2/2.0,1e-12);
            assert!(second.balance_error_j.abs() < 1e-12);
            assert_eq!(next.state().step,2);
            assert_eq!(accepted.state(),&saved);
        }
    });
}

#[test]
fn smooth_material_damping_evolution_has_second_order_temporal_convergence() {
    with_cx(false, |cx| {
        let fixed: Vec<_> = (0..12).filter(|d| *d != 3).collect();
        let geometry = geometry(cx);
        let frame = ContentHash([74;32]);
        let cfg = TransientConfig3 { mass_damping_per_s:0.4, stiffness_damping_s:0.04, ..config() };
        let load_value = |t: f64| {
            let k = 10.0*(1.0+t);
            (k-0.1)*0.001*t.cos()-(0.4*0.1+0.04*k)*0.001*t.sin()
        };
        let mut errors = Vec::new();
        for count in [16,32,64] {
            let materials: Vec<_> = (0..=count).map(|i|
                TetElasticMaterial::try_new(6.0,60.0*(1.0+i as f64/count as f64),
                    0.0,[73;32]).unwrap()).collect();
            let problems: Vec<_> = materials.iter().map(|m| problem(m,&fixed)).collect();
            let first = geometry.assemble(&[],VolumeLoad3::None,&[(3,load_value(0.0))],cx).unwrap();
            let mut initial = ZERO; initial[1][0] = 0.001;
            let mut run = first.begin_transient(&problems[0],None,0.0,&initial,&ZERO,cfg,cx).unwrap();
            for i in 1..=count {
                let t = i as f64/count as f64;
                let load = geometry.assemble(&[],VolumeLoad3::None,&[(3,load_value(t))],cx).unwrap();
                let zero = TetThermalStrainState::try_new(300.0+t,300.0,[0.0;6],
                    materials[i].material_state_identity,frame,0.01).unwrap();
                let (next,report) = run.advance_thermal_candidate(&problems[i],
                    TetThermalStrainField::Uniform(&zero),&load,1.0/count as f64,frame,0.01,cx).unwrap();
                assert!(report.native.damping_loss_j >= 0.0);
                assert!(report.physical_balance_error_j.abs() < 1e-12);
                run = next;
            }
            errors.push((run.state().displacement[0]-0.001*1.0_f64.cos()).abs()
                +(run.state().velocity[0]+0.001*1.0_f64.sin()).abs());
            assert_eq!(run.state().step,count);
        }
        for i in 0..2 {
            assert!((3.6..4.4).contains(&(errors[i]/errors[i+1])),"{errors:?}");
        }
    });
}

#[test]
fn changing_stiffness_damping_does_not_damp_rigid_translation() {
    with_cx(false, |cx| {
        let before = TetElasticMaterial::try_new(6.0,60.0,0.25,[73;32]).unwrap();
        let after = TetElasticMaterial::try_new(6.0,300.0,0.25,[75;32]).unwrap();
        let p0 = problem(&before,&[]); let p1 = problem(&after,&[]);
        let geometry = geometry(cx);
        let load = geometry.assemble(&[],VolumeLoad3::None,&[],cx).unwrap();
        let velocity = [[0.01,-0.02,0.03];4];
        let frame = ContentHash([74;32]);
        let zero = TetThermalStrainState::try_new(301.0,300.0,[0.0;6],
            after.material_state_identity,frame,0.01).unwrap();
        let run = load.begin_transient(&p0,None,0.0,&ZERO,&velocity,
            TransientConfig3 { stiffness_damping_s:0.1, ..config() },cx).unwrap();
        let (next,report) = run.advance_thermal_candidate(&p1,
            TetThermalStrainField::Uniform(&zero),&load,0.02,frame,0.01,cx).unwrap();
        let motion = next.nodal_motion(cx).unwrap();
        for node in 0..4 { for axis in 0..3 {
            close(motion[0][node][axis],0.02*velocity[node][axis],1e-11);
            close(motion[1][node][axis],velocity[node][axis],1e-10);
            close(motion[2][node][axis],0.0,1e-9);
        }}
        assert!(report.native.damping_loss_j.abs() < 1e-12);
        assert!(report.physical_balance_error_j.abs() < 1e-12);
    });
}

#[test]
fn fully_clamped_material_change_still_accounts_for_elastic_storage() {
    with_cx(false, |cx| {
        let fixed: Vec<_> = (0..12).collect();
        let before = TetElasticMaterial::try_new(6.0,60.0,0.0,[73;32]).unwrap();
        let after = TetElasticMaterial::try_new(6.0,30.0,0.0,[75;32]).unwrap();
        let p0=problem(&before,&fixed); let p1=problem(&after,&fixed);
        let geometry=geometry(cx);
        let load=geometry.assemble(&[],VolumeLoad3::None,&[],cx).unwrap();
        let frame=ContentHash([74;32]);
        let e0=TetThermalStrainState::try_new(310.0,300.0,[0.001,0.0,0.0,0.0,0.0,0.0],
            before.material_state_identity,frame,0.01).unwrap();
        let e1=TetThermalStrainState::try_new(320.0,300.0,[0.002,0.0,0.0,0.0,0.0,0.0],
            after.material_state_identity,frame,0.01).unwrap();
        let run=load.begin_transient(&p0,Some(TetThermalStrainField::Uniform(&e0)),
            0.0,&ZERO,&ZERO,TransientConfig3 { stiffness_damping_s:0.1,..config() },cx).unwrap();
        let (next,report)=run.advance_thermal_candidate(&p1,TetThermalStrainField::Uniform(&e1),
            &load,0.01,frame,0.01,cx).unwrap();
        assert!(next.state().displacement.is_empty());
        assert_eq!(report.native.damping_loss_j,0.0);
        close(report.kinetic_plus_elastic_before_j,5e-6,1e-15);
        close(report.kinetic_plus_elastic_after_j,1e-5,1e-15);
        close(report.thermal_material_work_j,1e-5,1e-15);
        close(report.stiffness_material_work_j,-5e-6,1e-15);
        assert!(report.physical_balance_error_j.abs() < 1e-15);
    });
}

#[test]
fn damped_material_trials_refuse_late_strain_and_exhausted_work_without_publication() {
    with_cx(false, |cx| {
        let fixed: Vec<_>=(0..12).filter(|d| *d!=3).collect();
        let before=TetElasticMaterial::try_new(6.0,60.0,0.0,[73;32]).unwrap();
        let after=TetElasticMaterial::try_new(6.0,30.0,0.0,[75;32]).unwrap();
        let p0=problem(&before,&fixed); let p1=problem(&after,&fixed);
        let geometry=geometry(cx);
        let load=geometry.assemble(&[],VolumeLoad3::None,&[],cx).unwrap();
        let frame=ContentHash([74;32]);
        let strain=TetThermalStrainState::try_new(310.0,300.0,[0.001,0.0,0.0,0.0,0.0,0.0],
            after.material_state_identity,frame,0.01).unwrap();
        let cfg=TransientConfig3 { stiffness_damping_s:0.1,..config() };
        let run=load.begin_transient(&p0,None,0.0,&ZERO,&ZERO,cfg,cx).unwrap();
        let saved=run.state().clone();
        // Initial stress is zero. The rejected endpoint has nonzero elastic
        // strain, after the midpoint and endpoint-acceleration solves succeed.
        assert!(run.advance_thermal_candidate(&p1,TetThermalStrainField::Uniform(&strain),
            &load,0.01,frame,1e-6,cx).is_err());
        assert_eq!(run.state(),&saved);
        let preview=run.advance_thermal_candidate(&p1,TetThermalStrainField::Uniform(&strain),
            &load,0.01,frame,0.01,cx).unwrap();
        assert_eq!(run.state(),&saved);
        let mut limited=cfg; limited.solve.maximum_iterations=1;
        let capped=load.begin_transient(&p0,None,0.0,&ZERO,&ZERO,limited,cx).unwrap();
        assert!(matches!(capped.advance_thermal_candidate(&p1,TetThermalStrainField::Uniform(&strain),
            &load,0.01,frame,0.01,cx),Err(TransientError3::Solve {..})|Err(TransientError3::Budget)));
        assert_eq!(capped.state(),&saved);
        with_cx(true, |cancelled| {
            assert!(matches!(run.advance_thermal_candidate(&p1,TetThermalStrainField::Uniform(&strain),
                &load,0.01,frame,0.01,cancelled),Err(TransientError3::Cancelled)));
        });
        let again=run.advance_thermal_candidate(&p1,TetThermalStrainField::Uniform(&strain),
            &load,0.01,frame,0.01,cx).unwrap();
        assert_eq!(preview.0.state(),again.0.state());
        assert_eq!(preview.1,again.1);
    });
}
