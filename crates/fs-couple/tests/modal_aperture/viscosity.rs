use super::*;
use fs_exec::CancelGate;

fn undamped() -> ModalApertureSpec {
    let mut s = spec();
    for m in &mut s.modes { m.damping_n_s_m = 0.0; }
    s
}
fn matrix() -> Vec<f64> { vec![0.03,0.005,-0.002, 0.005,0.02,0.003, -0.002,0.003,0.01] }
fn coupled(s: ModalApertureSpec, c: Vec<f64>) -> ModalAperture {
    model(s).with_viscosity(c, &CancelGate::new_clock_free()).unwrap()
}

#[test]
fn signed_modes_full_material_damping_and_channel_viscosity_close_one_work_balance() {
    let mut s = undamped(); s.rest_opening_m = 4e-5;
    let mut aperture = coupled(s, matrix());
    let law = fs_couple::bernoulli_aperture::viscous::ViscousSlit::new(1.81e-5, 0.001).unwrap();
    aperture.set_viscous_slit(law).unwrap();
    for sample in 0..256 {
        let input = drive(if sample<192 {2000.0} else {-100.0});
        let frame = aperture.step(input,[0.001,-0.0003]).unwrap();
        let expected = law.evaluate(aperture.spec().width_m,frame.aperture.midpoint_opening_m,
            aperture.spec().density_kg_m3,frame.pressure_drop_pa).unwrap();
        assert_eq!(frame.aperture.jet_flow_m3_s,expected.flow_m3_s);
        assert!(frame.momentum_residual_n.iter().all(|r|r.abs()<1e-8));
        let scale = frame.aperture.storage_change_j.abs()+frame.aperture.dissipated_energy_j
            +frame.aperture.pressure_work_j.abs()+frame.aperture.mechanical_work_j.abs();
        assert!(frame.aperture.balance_residual_j().abs()<1e-8*scale.max(1e-30));
    }
    assert!(aperture.modal_states().iter().all(|m|m.displacement_m!=0.0));
}

#[test]
fn full_viscosity_matches_manufactured_simultaneous_momentum_flow_and_work() {
    let mut s = undamped(); s.modes.truncate(2);
    s.modes[0].force_projection = [1.0,0.0];
    s.modes[1].force_projection = [0.0,1.0];
    let old = [BladeState { displacement_m:1e-6, velocity_m_s:0.01 },
        BladeState { displacement_m:-2e-6, velocity_m_s:0.005 }];
    let midpoint = [0.02,-0.006];
    for c in [[0.03,0.012,0.012,0.02], [0.03,-0.03,-0.03,0.03]] {
        for dp in [-100.0_f64,0.0,100.0] {
            let mut a = ModalAperture::new(s.clone(), old.to_vec(), contact()).unwrap()
                .with_viscosity(c.to_vec(), &CancelGate::new_clock_free()).unwrap();
            let dt = s.time_step_s;
            let opening = s.rest_opening_m + s.modes.iter().enumerate()
                .map(|(i,m)| m.gap_projection*(old[i].displacement_m + 0.5*dt*midpoint[i])).sum::<f64>();
            let jet = s.width_m*opening*dp.signum()*(2.0*dp.abs()/s.density_kg_m3).sqrt();
            let swept = -s.modes.iter().zip(midpoint).map(|(m,v)|m.pressure_area_m2*v).sum::<f64>();
            let incoming = 13.0;
            let body = -1e-7;
            let upstream = dp + 2.0*incoming + s.impedance_pa_s_m3*(jet + swept + body);
            let forces = std::array::from_fn(|i| {
                let m = s.modes[i];
                2.0*m.mass_kg*(midpoint[i]-old[i].velocity_m_s)/dt
                    + m.stiffness_n_m*(old[i].displacement_m+0.5*dt*midpoint[i])
                    + c[2*i]*midpoint[0]+c[2*i+1]*midpoint[1]+m.pressure_area_m2*dp
            });
            let frame = a.step(ApertureDrive { upstream_pressure_pa:upstream,
                incoming_pressure_pa:incoming, body_flow_m3_s:body }, forces).unwrap();
            for i in 0..2 {
                assert!((frame.modes[i].displacement_m-old[i].displacement_m-dt*midpoint[i]).abs()<1e-14);
                assert!((frame.modes[i].velocity_m_s-(2.0*midpoint[i]-old[i].velocity_m_s)).abs()<1e-9);
                assert!(frame.momentum_residual_n[i].abs()<1e-9);
            }
            assert!((frame.pressure_drop_pa-dp).abs()<1e-6);
            let power: f64 = (0..2).map(|i| (0..2).map(|j|midpoint[i]*c[2*i+j]*midpoint[j]).sum::<f64>()).sum();
            assert!((frame.aperture.dissipated_energy_j-dt*(dp*jet+power)).abs()<1e-12);
            assert!(frame.aperture.balance_residual_j().abs()<1e-12);
        }
    }
}

#[test]
fn relative_dashpot_preserves_common_motion_and_resolves_differential_decay() {
    use fs_couple::bernoulli_aperture::paired::BladePropagator;
    let mut s = undamped(); s.modes.truncate(2);
    s.modes[1].mass_kg = s.modes[0].mass_kg;
    s.modes[1].stiffness_n_m = s.modes[0].stiffness_n_m;
    for m in &mut s.modes { m.pressure_area_m2 = 0.0; }
    s.modes[1].gap_projection = -s.modes[0].gap_projection;
    let c = vec![0.03,-0.03,-0.03,0.03];
    for propagator in [BladePropagator::Midpoint, BladePropagator::ExactZoh] {
        for sign in [1.0,-1.0] {
            let old = vec![BladeState { displacement_m:0.0,velocity_m_s:0.1 },
                BladeState { displacement_m:0.0,velocity_m_s:sign*0.1 }];
            let mut a = ModalAperture::new_with_propagator(s.clone(),old.clone(),contact(),propagator)
                .unwrap().with_viscosity(c.clone(),&CancelGate::new_clock_free()).unwrap();
            let mut lossless = ModalAperture::new_with_propagator(s.clone(),old,contact(),propagator).unwrap();
            assert_eq!(a.viscosity_matrix_n_s_m(),Some(c.as_slice()));
            let initial = a.stored_energy_j(); let mut loss = 0.0;
            for _ in 0..512 {
                let frame = a.step(drive(0.0),[0.0;2]).unwrap();
                loss += frame.aperture.dissipated_energy_j;
                if sign == 1.0 {
                    let reference = lossless.step(drive(0.0),[0.0;2]).unwrap();
                    for (actual,expected) in frame.modes.iter().zip(reference.modes) {
                        assert!((actual.displacement_m-expected.displacement_m).abs()<1e-13);
                        assert!((actual.velocity_m_s-expected.velocity_m_s).abs()<1e-9);
                    }
                } else {
                    assert!(frame.aperture.dissipated_energy_j>=0.0);
                }
            }
            assert!((initial-a.stored_energy_j()-loss).abs()<1e-17);
            if sign == 1.0 { assert!(loss.abs()<1e-24); }
            else { assert!(loss>0.5*initial); }
            if propagator == BladePropagator::ExactZoh {
                let t = 512.0*s.time_step_s;
                let alpha = if sign == 1.0 {0.0} else {0.03/s.modes[0].mass_kg};
                let omega = (s.modes[0].stiffness_n_m/s.modes[0].mass_kg-alpha*alpha).sqrt();
                let q = 0.1*(-alpha*t).exp()*(omega*t).sin()/omega;
                let v = 0.1*(-alpha*t).exp()*((omega*t).cos()-alpha*(omega*t).sin()/omega);
                assert!((a.modal_states()[0].displacement_m-q).abs()<1e-13);
                assert!((a.modal_states()[0].velocity_m_s-v).abs()<1e-9);
            }
        }
    }
}

#[test]
fn singular_viscosity_retains_signed_physical_coordinate_scaling() {
    let mut s = undamped(); s.modes.truncate(2);
    let mut transformed = s.clone(); let scales = [-2.0,0.5];
    for (m,scale) in transformed.modes.iter_mut().zip(scales) {
        m.mass_kg *= scale*scale; m.stiffness_n_m *= scale*scale;
        m.pressure_area_m2 *= scale; m.gap_projection *= scale;
        for p in &mut m.force_projection { *p *= scale; }
    }
    let c = vec![0.03,-0.03,-0.03,0.03];
    let scaled = (0..4).map(|ij|c[ij]*scales[ij/2]*scales[ij%2]).collect();
    let (mut a,mut b) = (coupled(s,c),coupled(transformed,scaled));
    for n in 0..256 {
        let input = drive(if n<128 {2000.0} else {0.0});
        let f = [if n<50 {0.001} else {0.0},0.0];
        let x = a.step(input,f).unwrap(); let y = b.step(input,f).unwrap();
        for i in 0..2 {
            assert!((x.modes[i].displacement_m-y.modes[i].displacement_m*scales[i]).abs()<1e-11);
            assert!((x.modes[i].velocity_m_s-y.modes[i].velocity_m_s*scales[i]).abs()<1e-7);
        }
        assert!((x.aperture.dissipated_energy_j-y.aperture.dissipated_energy_j).abs()<1e-12);
        assert!((x.aperture.outgoing_pressure_pa-y.aperture.outgoing_pressure_pa).abs()<1e-5);
    }
}

#[test]
fn off_diagonal_material_loss_transfers_motion_without_adding_energy() {
    let mut s = undamped();
    for m in &mut s.modes { m.pressure_area_m2 = 0.0; }
    let old = vec![BladeState { displacement_m:0.0,velocity_m_s:0.1 },
        BladeState::default(),BladeState::default()];
    let mut full = ModalAperture::new(s.clone(),old.clone(),contact()).unwrap()
        .with_viscosity(matrix(),&CancelGate::new_clock_free()).unwrap();
    let mut diagonal = s;
    for (i,m) in diagonal.modes.iter_mut().enumerate() { m.damping_n_s_m = matrix()[3*i+i]; }
    let mut diagonal = ModalAperture::new(diagonal,old,contact()).unwrap();
    let initial = full.stored_energy_j();
    let mut loss = 0.0;
    let first = full.step(drive(0.0),[0.0;2]).unwrap();
    loss += first.aperture.dissipated_energy_j;
    assert!(first.modes[1].velocity_m_s.abs()>1e-6);
    assert_eq!(diagonal.step(drive(0.0),[0.0;2]).unwrap().modes[1],BladeState::default());
    for _ in 1..256 {
        let before = full.stored_energy_j();
        let frame = full.step(drive(0.0),[0.0;2]).unwrap();
        assert!(frame.aperture.dissipated_energy_j>=0.0);
        assert!(full.stored_energy_j()<=before+1e-18);
        loss += frame.aperture.dissipated_energy_j;
    }
    assert!((initial-full.stored_energy_j()-loss).abs()<1e-17);
}

#[test]
fn zero_tensor_preserves_legacy_arithmetic_and_signed_congruence_preserves_physics() {
    let s = undamped();
    let (mut legacy,mut zero) = (model(s.clone()),coupled(s.clone(),vec![0.0;9]));
    let mut transformed = s.clone(); let scales = [-2.0,0.5,-1.5];
    for (m,scale) in transformed.modes.iter_mut().zip(scales) {
        m.mass_kg *= scale*scale; m.stiffness_n_m *= scale*scale;
        m.pressure_area_m2 *= scale; m.gap_projection *= scale;
        for p in &mut m.force_projection { *p *= scale; }
    }
    let mut c = matrix();
    for i in 0..3 {for j in 0..3 {c[3*i+j] *= scales[i]*scales[j];}}
    let (mut a,mut b) = (coupled(s,matrix()),coupled(transformed,c));
    for n in 0..512 {
        let input = drive(if n<200 {2000.0} else {0.0});
        let f = [if n<50 {0.001} else {0.0},0.0];
        assert_eq!(legacy.step(input,f).unwrap(),zero.step(input,f).unwrap());
        let x=a.step(input,f).unwrap(); let y=b.step(input,f).unwrap();
        for i in 0..3 {
            assert!((x.modes[i].displacement_m-y.modes[i].displacement_m*scales[i]).abs()<1e-11);
            assert!((x.modes[i].velocity_m_s-y.modes[i].velocity_m_s*scales[i]).abs()<1e-7);
        }
        assert!((x.aperture.outgoing_pressure_pa-y.aperture.outgoing_pressure_pa).abs()<1e-5);
        assert!((x.aperture.stored_energy_j-y.aperture.stored_energy_j).abs()<1e-12);
    }
}

#[test]
fn full_viscosity_refuses_invalid_laws_and_failed_steps_preserve_accepted_state() {
    for c in [vec![], vec![f64::NAN;9], vec![-1.0;9],
        vec![1.0,2.0,0.0,2.0,1.0,0.0,0.0,0.0,1.0],
        vec![1.0,0.1,0.0,0.2,1.0,0.0,0.0,0.0,1.0],
        vec![0.0,0.1,0.0,0.1,1.0,0.0,0.0,0.0,1.0]] {
        assert!(model(undamped()).with_viscosity(c,&CancelGate::new_clock_free()).is_err());
    }
    assert!(model(spec()).with_viscosity(matrix(),&CancelGate::new_clock_free()).is_err());
    let gate=CancelGate::new_clock_free();gate.request();
    assert!(model(undamped()).with_viscosity(matrix(),&gate).is_err());
    let mut a=coupled(undamped(),matrix());a.step(drive(50.0),[0.0;2]).unwrap();
    let before=a.modal_states().to_vec();let energy=a.stored_energy_j();
    for p in [f64::NAN,f64::INFINITY,f64::MAX] {
        assert!(a.step(drive(p),[0.0;2]).is_err());
        assert_eq!(a.modal_states(),before);assert_eq!(a.stored_energy_j(),energy);
        assert_eq!(a.accepted_steps(),1);
    }
    assert!(a.with_viscosity(matrix(),&CancelGate::new_clock_free()).is_err());
    assert!(coupled(undamped(),matrix()).with_viscosity(matrix(),&CancelGate::new_clock_free()).is_err());
    let c=vec![0.03,0.005,0.0,0.005,0.02,0.0,0.0,0.0,0.0];
    let mut mixed=coupled(undamped(),c.clone());
    assert_eq!(mixed.viscosity_matrix_n_s_m(),Some(c.as_slice()));
    mixed.step(drive(30.0),[0.0;2]).unwrap();
}

#[path = "viscosity/material.rs"]
mod material;
