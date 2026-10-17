use super::*;
use super::super::dynamic::{DynamicAperture, DynamicApertureSpec};
use super::super::BernoulliAperture;

fn spec() -> PairedApertureSpec {
    PairedApertureSpec { rest_opening_m: 4e-4, width_m: 0.013,
        blades: [ApertureBlade { mass_kg: 1e-5, stiffness_n_m: 350.0,
            damping_n_s_m: 0.03, pressure_area_m2: 1e-4 }; 2],
        density_kg_m3: 1.2, impedance_pa_s_m3: 2.7e7, time_step_s: 1.0/3072000.0,
        max_steps: 2000*64, limits: PairedApertureLimits {
            maximum_abs_displacement_m: 0.1, maximum_abs_velocity_m_s: 100.0,
            maximum_energy_j: 10.0, max_evaluations: 256,
            velocity_tolerance_m_s: 1e-15, relative_tolerance: 1e-13,
        } }
}
fn contact(k: f64) -> Obstacle {
    Obstacle::new(vec![-1.0],1,1,vec![0.0],vec![1.0],k,2.0,
        "synthetic shared blade interface; not measured cane".into()).unwrap()
        .with_internal_loss(0.15).unwrap()
}

#[test]
fn equal_blades_match_the_existing_symmetric_midpoint_image() {
    let s = spec(); let b=s.blades[0];
    let mut pair=PairedAperture::new(s,[BladeState::default();2],contact(2e9)).unwrap();
    let mut scalar=DynamicAperture::new(DynamicApertureSpec {
        aperture: BernoulliAperture { rest_opening_m:s.rest_opening_m,width_m:s.width_m,
            closing_pressure_pa:0.5*b.stiffness_n_m*s.rest_opening_m/b.pressure_area_m2 },
        mass_kg:0.5*b.mass_kg,stiffness_n_m:0.5*b.stiffness_n_m,
        damping_ratio:b.damping_n_s_m/(2.0*(b.mass_kg*b.stiffness_n_m).sqrt()),
        density_kg_m3:s.density_kg_m3,impedance_pa_s_m3:s.impedance_pa_s_m3,
        time_step_s:s.time_step_s,max_steps:s.max_steps,
    }, ApertureState {opening_m:s.rest_opening_m,opening_velocity_m_s:0.0},contact(2e9)).unwrap();
    for i in 0..240 * 64 {
        let drive=ApertureDrive {upstream_pressure_pa:if i<120*64 {1200.0} else {-300.0},
            incoming_pressure_pa:75.0,body_flow_m3_s:1e-7};
        let a = pair.step(drive).unwrap();
        let r = scalar.step(drive).unwrap();
        assert!((a.aperture.state.opening_m-r.state.opening_m).abs()<2e-10);
        assert!((a.aperture.state.opening_velocity_m_s-r.state.opening_velocity_m_s).abs()<2e-6);
        assert!((a.aperture.bore_pressure_pa-r.bore_pressure_pa).abs()<2e-3);
        assert!((a.aperture.stored_energy_j-r.stored_energy_j).abs()<2e-9);
        assert_eq!(a.blades[0],a.blades[1]);
    }
}

#[test]
fn independent_motion_is_not_collapsed_into_the_total_gap() {
    let s=spec();let initial=[BladeState {displacement_m:1e-5,velocity_m_s:0.1},
        BladeState {displacement_m:-1e-5,velocity_m_s:-0.1}];
    let mut pair=PairedAperture::new(s,initial,contact(0.0)).unwrap();
    let before=pair.stored_energy_j();
    let frame=pair.step(ApertureDrive::default()).unwrap();
    assert!((frame.aperture.state.opening_m-s.rest_opening_m).abs()<1e-18);
    assert!(frame.aperture.state.opening_velocity_m_s.abs()<1e-13);
    assert!(frame.aperture.bore_flow_m3_s.abs()<1e-14);
    assert!(frame.blades[0].velocity_m_s.abs()>0.01);
    assert!(frame.aperture.stored_energy_j>0.0 && frame.aperture.stored_energy_j<before);
    assert!(frame.aperture.balance_residual_j().abs()<1e-16);
}

#[test]
fn unequal_blades_satisfy_both_momenta_and_the_same_fluid_equation() {
    // Manufacture individual midpoint velocities and pressure. Required forces
    // follow each original momentum equation, independent of condensation/root.
    for dp in [-1000.0_f64,0.0,1000.0] {
        for k in [0.0,2e9] {
            let mut s=spec();s.blades[1]=ApertureBlade {mass_kg:2.5e-5,stiffness_n_m:700.0,
                damping_n_s_m:0.08,pressure_area_m2:1.6e-4};
            let old=[BladeState {displacement_m:-2.6e-4,velocity_m_s:0.08},
                BladeState {displacement_m:-1.6e-4,velocity_m_s:-0.02}];
            let velocities=[-0.05,0.01];let dt=s.time_step_s;
            let gap=s.rest_opening_m+old[0].displacement_m+old[1].displacement_m;
            let next_gap=gap+dt*(velocities[0]+velocities[1]);
            let law=contact(k);
            let fc=ApertureContactStep::new(&law,gap).unwrap()
                .response(next_gap,velocities[0]+velocities[1]).unwrap().force;
            let forces=std::array::from_fn(|i| {let b=s.blades[i];
                b.mass_kg*2.0*(velocities[i]-old[i].velocity_m_s)/dt
                    +b.stiffness_n_m*(old[i].displacement_m+0.5*dt*velocities[i])
                    +b.damping_n_s_m*velocities[i]+b.pressure_area_m2*dp-fc});
            let jet=s.width_m*f64::midpoint(gap,next_gap).max(0.0)*dp.signum()
                *(2.0*dp.abs()/s.density_kg_m3).sqrt();
            let swept=-(0..2).map(|i|s.blades[i].pressure_area_m2*velocities[i]).sum::<f64>();
            let incoming=75.0;let body=2e-7;
            let upstream=dp+2.0*incoming+s.impedance_pa_s_m3*(jet+swept+body);
            let mut pair=PairedAperture::new(s,old,law).unwrap();
            let f=pair.step_with_forces(ApertureDrive {upstream_pressure_pa:upstream,
                incoming_pressure_pa:incoming,body_flow_m3_s:body},forces).unwrap();
            for i in 0..2 {
                assert!((f.blades[i].velocity_m_s-(2.0*velocities[i]-old[i].velocity_m_s)).abs()<1e-10);
                assert!(f.momentum_residual_n[i].abs()<1e-9);
            }
            assert!((f.pressure_drop_pa-dp).abs()<1e-7);
            assert!(f.aperture.balance_residual_j().abs()<1e-13);
            assert!(f.aperture.flow_residual_m3_s.abs()<1e-14);
        }
    }
}

#[test]
fn a_mass_only_data_change_alters_dynamics_without_retuning_elastic_closure() {
    let s = spec();
    let mut heavier = s;
    heavier.blades[1].mass_kg *= 4.0;
    let mut a=PairedAperture::new(s,[BladeState::default();2],contact(0.0)).unwrap();
    let mut b=PairedAperture::new(heavier,[BladeState::default();2],contact(0.0)).unwrap();
    let drive = ApertureDrive {
        upstream_pressure_pa: 1000.0,
        ..ApertureDrive::default()
    };
    // Observe 0.25 ms of actual motion. One original 48 kHz tick is only
    // 20.8 us: after resolving its pressure feedback the mass-only velocity
    // separation has not yet reached the unchanged 0.01 m/s observable.
    for _ in 1..768 {
        a.step(drive).unwrap();
        b.step(drive).unwrap();
    }
    let fa = a.step(drive).unwrap();
    let fb = b.step(drive).unwrap();
    assert!((fa.blades[1].velocity_m_s-fb.blades[1].velocity_m_s).abs()>0.01);
    assert!((fa.aperture.bore_pressure_pa-fb.aperture.bore_pressure_pa).abs()>0.01);
    // The source's static compliance is unchanged; no frequency was prescribed.
    assert_eq!(s.blades[1].pressure_area_m2/s.blades[1].stiffness_n_m,
        heavier.blades[1].pressure_area_m2/heavier.blades[1].stiffness_n_m);
}

#[test]
fn refusal_and_budget_exhaustion_publish_neither_blade() {
    let mut s=spec();s.max_steps=1;
    let mut pair=PairedAperture::new(s,[BladeState::default();2],contact(0.0)).unwrap();
    assert!(pair.step(ApertureDrive {upstream_pressure_pa:f64::NAN,..ApertureDrive::default()}).is_err());
    assert_eq!(pair.accepted_steps(),0);
    assert_eq!(pair.blade_states(),&[BladeState::default();2]);
    pair.step(ApertureDrive {upstream_pressure_pa:1000.0,..ApertureDrive::default()}).unwrap();
    let old=*pair.blade_states();let energy=pair.stored_energy_j();
    assert!(pair.step(ApertureDrive::default()).is_err());
    assert_eq!(pair.blade_states(),&old);assert_eq!(pair.stored_energy_j(),energy);
    pair.extend_step_budget(2).unwrap();pair.step(ApertureDrive::default()).unwrap();
    assert_eq!(pair.accepted_steps(),2);
    let mut s=spec();s.limits.max_evaluations=1;
    let mut pair=PairedAperture::new(s,[BladeState::default();2],contact(0.0)).unwrap();
    assert!(pair.step(ApertureDrive {upstream_pressure_pa:1000.0,..ApertureDrive::default()}).is_err());
    assert_eq!(pair.accepted_steps(),0);assert_eq!(pair.stored_energy_j(),0.0);
}

/// G1: at omega*dt = 1, where Midpoint refuses (it would be 130.7 c flat), the
/// exact blocks follow the closed-form damped oscillator of each blade.
#[test]
fn exact_blades_follow_the_damped_oscillator_beyond_the_midpoint_band() {
    use fs_math::det;
    // A negligible pressure face isolates the mechanical map (A*dp ~ 1e-26 N).
    let blade=ApertureBlade {mass_kg:1e-5,stiffness_n_m:500.0,damping_n_s_m:2e-3,pressure_area_m2:1e-12};
    let omega=det::sqrt(blade.stiffness_n_m/blade.mass_kg);
    let mut s=spec();
    s.blades=[blade;2];s.impedance_pa_s_m3=1.0;s.time_step_s=1.0/omega;
    let initial=[BladeState {displacement_m:1e-5,velocity_m_s:0.0},
        BladeState {displacement_m:0.0,velocity_m_s:0.05}];
    assert!(matches!(PairedAperture::new(s,initial,contact(0.0)),Err(AcousticRealizeError::TimeResolution(_))));
    let mut pair=PairedAperture::new_with_propagator(s,initial,contact(0.0),BladePropagator::ExactZoh).unwrap();
    assert_eq!(pair.blade_propagator(),BladePropagator::ExactZoh);
    let zeta=blade.damping_n_s_m/(2.0*det::sqrt(blade.stiffness_n_m*blade.mass_kg));
    let wd=omega*det::sqrt(1.0-zeta*zeta);
    let closed=|b:BladeState,t:f64| {
        let decay=det::exp(-zeta*omega*t);
        let (c,sn)=(det::cos(wd*t),det::sin(wd*t));
        (decay*(b.displacement_m*c+(b.velocity_m_s+zeta*omega*b.displacement_m)/wd*sn),
         decay*(b.velocity_m_s*c-(omega*omega*b.displacement_m+zeta*omega*b.velocity_m_s)/wd*sn))
    };
    let energy=pair.stored_energy_j();
    for step in 1..=100 {
        let f=pair.step(ApertureDrive::default()).unwrap();
        let t=f64::from(step)*s.time_step_s;
        for i in 0..2 {
            let (x,v)=closed(initial[i],t);
            let amplitude=initial[i].displacement_m.abs()+initial[i].velocity_m_s.abs()/omega;
            assert!((f.blades[i].displacement_m-x).abs()<=1e-11*amplitude,"step {step} blade {i}");
            assert!((f.blades[i].velocity_m_s-v).abs()<=1e-11*omega*amplitude,"step {step} blade {i}");
            assert!(f.momentum_residual_n[i].abs()<=1e-10*blade.stiffness_n_m*amplitude);
        }
        assert!(f.aperture.balance_residual_j().abs()<=1e-12*energy);
        assert!(f.aperture.dissipated_energy_j>=-1e-15*energy);
    }
}

/// G0/G1: the authored reed blades (about 1.1 kHz, heavy pressure coupling)
/// refuse under Midpoint at 48 kHz. The exact blocks run the shared contact
/// and Bernoulli junction there with closed work, flow and momentum ledgers.
#[test]
fn exact_blades_run_the_reed_band_with_closed_ledgers() {
    let b=ApertureBlade {mass_kg:1e-5,stiffness_n_m:500.0,damping_n_s_m:0.03,pressure_area_m2:1.25e-4};
    let mut s=spec();
    s.blades=[b,ApertureBlade {mass_kg:1.3e-5,..b}];s.time_step_s=1.0/48_000.0;s.max_steps=4_800;
    assert!(matches!(PairedAperture::new(s,[BladeState::default();2],contact(2e9)),
        Err(AcousticRealizeError::TimeResolution(_))));
    let mut pair=PairedAperture::new_with_propagator(s,[BladeState::default();2],contact(2e9),
        BladePropagator::ExactZoh).unwrap();
    let mut moved=0.0_f64;
    for i in 0..4_800 {
        let drive=ApertureDrive {upstream_pressure_pa:if i<2_400 {1200.0} else {-300.0},
            incoming_pressure_pa:75.0,body_flow_m3_s:1e-7};
        // Each accepted step already passed the internal work/flow balance gate.
        let f=pair.step(drive).unwrap();
        for k in 0..2 {
            let (blade,state)=(s.blades[k],f.blades[k]);
            let scale=blade.mass_kg*state.velocity_m_s.abs()/s.time_step_s
                +blade.stiffness_n_m*state.displacement_m.abs()+blade.damping_n_s_m*state.velocity_m_s.abs()
                +blade.pressure_area_m2*f.pressure_drop_pa.abs()+f.blade_contact_forces_n[k].abs()+1e-12;
            assert!(f.momentum_residual_n[k].abs()<=1e-9*scale,"step {i}: {:?}",f.momentum_residual_n);
        }
        moved=moved.max((f.blades[0].displacement_m-f.blades[1].displacement_m).abs());
    }
    // The unequal masses keep a real antisymmetric state.
    assert!(moved>0.0);
}

/// G1: the exact single-coordinate valve follows the closed-form damped
/// oscillator at omega*dt = 1, where its Midpoint constructor refuses.
#[test]
fn exact_dynamic_aperture_follows_the_damped_oscillator() {
    use fs_math::det;
    // A 1e12 Pa closing pressure leaves a 2e-13 m^2 face: the mechanics are isolated.
    let (m,k,zeta)=(1e-5,500.0,0.02);
    let omega=det::sqrt(k/m);
    let spec=DynamicApertureSpec {
        aperture: BernoulliAperture {rest_opening_m:4e-4,width_m:0.013,closing_pressure_pa:1e12},
        mass_kg:m,stiffness_n_m:k,damping_ratio:zeta,density_kg_m3:1.2,impedance_pa_s_m3:1.0,
        time_step_s:1.0/omega,max_steps:100,
    };
    let x0=1e-5;
    let start=ApertureState {opening_m:4e-4+x0,opening_velocity_m_s:0.0};
    assert!(matches!(DynamicAperture::new(spec,start,contact(0.0)),Err(AcousticRealizeError::TimeResolution(_))));
    let mut valve=DynamicAperture::new_with_propagator(spec,start,contact(0.0),BladePropagator::ExactZoh).unwrap();
    let wd=omega*det::sqrt(1.0-zeta*zeta);
    let energy=valve.stored_energy_j();
    for step in 1..=100 {
        let f=valve.step(ApertureDrive::default()).unwrap();
        let t=f64::from(step)*spec.time_step_s;
        let decay=det::exp(-zeta*omega*t);
        let (c,s)=(det::cos(wd*t),det::sin(wd*t));
        let x=decay*(x0*c+zeta*omega*x0/wd*s);
        let v=-decay*omega*omega*x0/wd*s;
        assert!((f.state.opening_m-4e-4-x).abs()<=1e-11*x0,"step {step}");
        assert!((f.state.opening_velocity_m_s-v).abs()<=1e-11*omega*x0,"step {step}");
        assert!(f.balance_residual_j().abs()<=1e-12*energy,"step {step}: {}",f.balance_residual_j());
    }
}

/// G1: two equal exact blades are the symmetric image of one exact opening
/// coordinate (blade mass, stiffness and damping twice the opening's), through
/// two independent junction solvers, where both Midpoint constructors refuse.
#[test]
fn exact_dynamic_aperture_matches_the_exact_symmetric_paired_image() {
    let mut s=spec();s.time_step_s=1.0/48_000.0;s.max_steps=2_400;
    let b=s.blades[0];
    let dynamic_spec=DynamicApertureSpec {
        aperture: BernoulliAperture {rest_opening_m:s.rest_opening_m,width_m:s.width_m,
            closing_pressure_pa:0.5*b.stiffness_n_m*s.rest_opening_m/b.pressure_area_m2},
        mass_kg:0.5*b.mass_kg,stiffness_n_m:0.5*b.stiffness_n_m,
        damping_ratio:b.damping_n_s_m/(2.0*(b.mass_kg*b.stiffness_n_m).sqrt()),
        density_kg_m3:s.density_kg_m3,impedance_pa_s_m3:s.impedance_pa_s_m3,
        time_step_s:s.time_step_s,max_steps:s.max_steps,
    };
    let rest=ApertureState {opening_m:s.rest_opening_m,opening_velocity_m_s:0.0};
    assert!(PairedAperture::new(s,[BladeState::default();2],contact(2e9)).is_err());
    assert!(DynamicAperture::new(dynamic_spec,rest,contact(2e9)).is_err());
    let mut pair=PairedAperture::new_with_propagator(s,[BladeState::default();2],contact(2e9),
        BladePropagator::ExactZoh).unwrap();
    let mut scalar=DynamicAperture::new_with_propagator(dynamic_spec,rest,contact(2e9),
        BladePropagator::ExactZoh).unwrap();
    for i in 0..2_400 {
        let drive=ApertureDrive {upstream_pressure_pa:if i<1_200 {1200.0} else {-300.0},
            incoming_pressure_pa:75.0,body_flow_m3_s:1e-7};
        let a=pair.step(drive).unwrap();
        let r=scalar.step(drive).unwrap();
        assert!((a.aperture.state.opening_m-r.state.opening_m).abs()<2e-10,"step {i}");
        assert!((a.aperture.state.opening_velocity_m_s-r.state.opening_velocity_m_s).abs()<2e-6,"step {i}");
        assert!((a.aperture.bore_pressure_pa-r.bore_pressure_pa).abs()<2e-3,"step {i}");
        assert!((a.aperture.stored_energy_j-r.stored_energy_j).abs()<2e-9,"step {i}");
        assert!(r.balance_residual_j().abs()<=1e-9*r.stored_energy_j.abs().max(1e-12),"step {i}");
        assert_eq!(a.blades[0],a.blades[1]);
    }
}
