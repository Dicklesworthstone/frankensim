//! Independent manufactured momentum/work checks for spatially distinct contacts.
use super::*;
fn distributed() -> Obstacle {
    Obstacle::new(vec![-0.4,-1.0,-1.7,0.3],4,1,
        vec![-2e-5,2e-5,1e-5,-1e-5],vec![0.1,0.2,0.3,0.4],1e8,2.0,
        "synthetic four-point affine lay".into()).unwrap().with_internal_loss(0.5).unwrap()
}
fn potential(law: &Obstacle, q: f64) -> f64 {
    law.collocation().iter().zip(law.gaps()).zip(law.weights()).map(|((&b,&gap),&w)|
        w*law.stiffness()*(b*q-gap).max(0.0).powi(3)/3.0).sum()
}
// Independent exact polynomial secant for alpha=2, not the adapter under test.
fn reaction(law: &Obstacle, q0: f64, q1: f64, vm: f64) -> (f64,f64) {
    let mut force = 0.0;
    let mut loss = 0.0;
    for ((&b,&gap),&weight) in law.collocation().iter().zip(law.gaps()).zip(law.weights()) {
        let a = (b*q0-gap).max(0.0);
        let z = (b*q1-gap).max(0.0);
        let elastic = if q0 == q1 || b == 0.0 { weight*law.stiffness()*a*a }
            else if a > 0.0 && z > 0.0 { weight*law.stiffness()*(a*a+a*z+z*z)/3.0 }
            else { weight*law.stiffness()*(z.powi(3)-a.powi(3))/(3.0*b*(q1-q0)) };
        let local_velocity = -b*vm;
        let f = (elastic*(1.0-law.internal_loss()*local_velocity)).max(0.0);
        force -= b*f;
        loss += (elastic-f)*local_velocity;
    }
    (force,loss)
}
#[test]
fn multiple_lay_points_solve_actual_momentum_pressure_and_contact_work_together() {
    let law = distributed();
    // The pole-rate bound is about 12100 rad/s: midpoint runs inside its
    // band, the exact block at the original 1e-5 s clock. Both hold the same
    // secant contact reaction over the step.
    for (propagator,dt) in [(BladePropagator::Midpoint,3e-6),(BladePropagator::ExactZoh,1e-5)] {
    let mut s = spec(1); s.time_step_s = dt;
    let (m,k,rest) = (s.mass_kg,s.stiffness_n_m,s.aperture.rest_opening_m);
    let area = k*rest/s.aperture.closing_pressure_pa;
    let damping = 2.0*s.damping_ratio*(m*k).sqrt();
    for (y,next_y,v) in [(-1e-4,-9e-5,0.2),(2e-5,-1e-5,-0.1),(-1e-5,2e-5,0.1)] {
        let vm = (next_y-y)/dt;
        let (force, contact_loss) = reaction(&law,y,next_y,vm);
        // The total force the propagator held, its end velocity and viscous loss.
        let (held,v1,viscous) = match propagator {
            BladePropagator::Midpoint => {
                let v1 = 2.0*vm-v;
                (m*(v1-v)/dt+k*(f64::midpoint(y,next_y)-rest)+damping*vm,v1,dt*damping*vm*vm)
            }
            BladePropagator::ExactZoh => {
                let (held,v1) = exact_held_force(&s,dt,(y,next_y,v));
                // Energy theorem under a constant force: H dy = dKE + dPE + viscous.
                (held,v1,held*(next_y-y)-0.5*m*(v1*v1-v*v)-0.5*k*((next_y-rest).powi(2)-(y-rest).powi(2)))
            }
        };
        let dp = (force-held)/area;
        let jet = s.aperture.width_m*f64::midpoint(y,next_y).max(0.0)*dp.signum()
            *(2.0*dp.abs()/s.density_kg_m3).sqrt();
        let incoming = 75.0;
        let body = 2e-7;
        let outgoing = incoming+s.impedance_pa_s_m3*(jet-area*vm+body);
        let mut valve = DynamicAperture::new_with_propagator(s,ApertureState { opening_m:y, opening_velocity_m_s:v },
            law.clone(),propagator).unwrap();
        let f = valve.step(ApertureDrive { upstream_pressure_pa:dp+outgoing+incoming,
            incoming_pressure_pa:incoming,body_flow_m3_s:body }).unwrap();
        assert!((f.state.opening_m-next_y).abs() < 1e-13,"{propagator:?}");
        assert!((f.state.opening_velocity_m_s-v1).abs() < 1e-8,"{propagator:?}");
        assert!((f.outgoing_pressure_pa-outgoing).abs() < 1e-6,"{propagator:?}");
        let storage = 0.5*m*v1*v1+0.5*k*(next_y-rest).powi(2)+potential(&law,next_y);
        assert!((f.stored_energy_j-storage).abs() < 1e-12*storage.max(1e-20),"{propagator:?}");
        assert!(viscous > 0.0,"{propagator:?}");
        let expected_loss = dt*(dp*jet+contact_loss)+viscous;
        assert!((f.dissipated_energy_j-expected_loss).abs() < 1e-9*expected_loss.max(1e-20),"{propagator:?}");
        let scale = f.stored_energy_j+f.storage_change_j.abs()+f.dissipated_energy_j+f.pressure_work_j.abs();
        assert!(f.balance_residual_j().abs() <= 1e-10*scale,"{propagator:?}");
    }
    }
}
#[test]
fn physical_gap_distribution_changes_motion_without_collapsing_to_a_mean_lay() {
    let initial = ApertureState { opening_m:-1e-4,opening_velocity_m_s:-0.2 };
    let law = distributed();
    let altered = Obstacle::new(law.collocation().to_vec(),4,1,vec![4e-4;4],law.weights().to_vec(),
        law.stiffness(),law.alpha(),law.provenance().into()).unwrap().with_internal_loss(law.internal_loss()).unwrap();
    let mut a = DynamicAperture::new(spec(32),initial,law.clone()).unwrap();
    let mut b = DynamicAperture::new(spec(32),initial,altered).unwrap();
    let f = a.step(ApertureDrive::default()).unwrap();
    let g = b.step(ApertureDrive::default()).unwrap();
    assert!((f.state.opening_velocity_m_s-g.state.opening_velocity_m_s).abs() > 1e-4);
    assert!((f.swept_flow_m3_s-g.swept_flow_m3_s).abs() > 1e-10);
    assert_eq!(a.contact_law().collocation(),law.collocation());
    assert_eq!(a.contact_law().gaps(),law.gaps());
    assert_eq!(a.contact_law().provenance(),law.provenance());
}
#[test]
fn distributed_state_survives_refusal_and_cancelled_callback_for_exact_retry() {
    let initial = ApertureState { opening_m:-1e-4,opening_velocity_m_s:0.2 };
    let mut a = DynamicAperture::new(spec(32),initial,distributed()).unwrap();
    let mut b = DynamicAperture::new(spec(32),initial,distributed()).unwrap();
    let valid = ApertureDrive { upstream_pressure_pa:500.0,..ApertureDrive::default() };
    assert_eq!(bits(a.step(valid).unwrap()),bits(b.step(valid).unwrap()));
    let state = a.state();
    assert!(a.step(ApertureDrive { body_flow_m3_s:f64::NAN,..valid }).is_err());
    let gate = CancelGate::new_clock_free(); gate.request();
    let sentinel = ApertureFrame { step:u64::MAX,..ApertureFrame::default() };
    let mut output = [sentinel;8];
    let progress = a.advance_block(&[valid;8],&mut output,&gate).unwrap();
    assert_eq!(progress.completed,0); assert_eq!(a.state(),state);
    assert!(output.iter().all(|f| f.step==u64::MAX));
    for _ in 0..20 { assert_eq!(bits(a.step(valid).unwrap()),bits(b.step(valid).unwrap())); }
}
