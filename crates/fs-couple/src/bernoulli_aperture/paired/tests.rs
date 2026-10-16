use super::*;
use super::super::dynamic::{DynamicAperture, DynamicApertureSpec};
use super::super::BernoulliAperture;

fn spec() -> PairedApertureSpec {
    PairedApertureSpec { rest_opening_m: 4e-4, width_m: 0.013,
        blades: [ApertureBlade { mass_kg: 1e-5, stiffness_n_m: 350.0,
            damping_n_s_m: 0.03, pressure_area_m2: 1e-4 }; 2],
        density_kg_m3: 1.2, impedance_pa_s_m3: 2.7e7, time_step_s: 1.0/48000.0,
        max_steps: 2000, limits: PairedApertureLimits {
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
    for i in 0..240 {
        let drive=ApertureDrive {upstream_pressure_pa:if i<120 {1200.0} else {-300.0},
            incoming_pressure_pa:75.0,body_flow_m3_s:1e-7};
        let a=pair.step(drive).unwrap();let r=scalar.step(drive).unwrap();
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
    let s=spec();let mut heavier=s;heavier.blades[1].mass_kg*=4.0;
    let mut a=PairedAperture::new(s,[BladeState::default();2],contact(0.0)).unwrap();
    let mut b=PairedAperture::new(heavier,[BladeState::default();2],contact(0.0)).unwrap();
    let drive=ApertureDrive {upstream_pressure_pa:1000.0,..ApertureDrive::default()};
    let fa=a.step(drive).unwrap();let fb=b.step(drive).unwrap();
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
