//! Both existing valve images retain the same gas-derived fractional loss owner.
use super::*;
use fs_couple::bernoulli_aperture::{BernoulliAperture,
    dynamic::{ApertureState, DynamicAperture, DynamicApertureSpec},
    paired::{ApertureBlade, BladeState, PairedAperture, PairedApertureLimits, PairedApertureSpec},
    two_sided::{ApertureDuctSpec, LoadedAperture, TwoSidedAperture, TwoSidedDrive,
        TwoSidedLimits, ValveSide, combined_inlet_state}};
use fs_dcontact::Obstacle;

fn model(upstream_length: f64, paired: bool) -> TwoSidedAperture {
    let g=gas(293.15);
    let sources=[upstream_length,0.4001].map(|length| ApertureDuctSpec {
        graph: graph(g,length), gases: vec![g], viscothermal: vec![selection(0)] });
    let (_, z)=combined_inlet_state(&sources).unwrap();
    let b=ApertureBlade {mass_kg:1e-5,stiffness_n_m:500.0,damping_n_s_m:0.03,
        pressure_area_m2:1.25e-4};
    let contact=Obstacle::new(vec![-1.0],1,1,vec![0.0],vec![1.0],1e8,2.0,
        "synthetic fractional loss valve fixture".into()).unwrap();
    let valve:LoadedAperture=if paired {
        let mut blades=[b;2];blades[1].mass_kg*=1.7;
        PairedAperture::new(PairedApertureSpec {
            rest_opening_m:4e-4,width_m:0.013,blades,density_kg_m3:g.density,
            impedance_pa_s_m3:z,time_step_s:DT,max_steps:1024,
            limits:PairedApertureLimits {maximum_abs_displacement_m:0.01,
                maximum_abs_velocity_m_s:100.0,maximum_energy_j:1.0,max_evaluations:256,
                velocity_tolerance_m_s:1e-13,relative_tolerance:1e-11},
        },[BladeState::default();2],contact).unwrap().into()
    } else {
        DynamicAperture::new(DynamicApertureSpec {
            aperture:BernoulliAperture {rest_opening_m:4e-4,width_m:0.013,
                closing_pressure_pa:0.5*b.stiffness_n_m*4e-4/b.pressure_area_m2},
            mass_kg:0.5*b.mass_kg,stiffness_n_m:0.5*b.stiffness_n_m,
            damping_ratio:b.damping_n_s_m/(2.0*(b.mass_kg*b.stiffness_n_m).sqrt()),
            density_kg_m3:g.density,impedance_pa_s_m3:z,time_step_s:DT,max_steps:1024,
        },ApertureState {opening_m:4e-4,opening_velocity_m_s:0.0},contact).unwrap().into()
    };
    TwoSidedAperture::with_fractional_transits(valve,sources,TwoSidedLimits {
        maximum_energy_j:1.0,energy_absolute_tolerance_j:1e-16,energy_relative_tolerance:1e-7,
        flow_absolute_tolerance_m3_s:1e-16,flow_relative_tolerance:1e-8,
    },[policy();2],&CancelGate::new_clock_free()).unwrap()
}

#[test]
fn original_single_and_independent_blade_sources_use_unrounded_lossy_tract_feedback() {
    let gate=CancelGate::new_clock_free();
    for paired in [false,true] {
        let mut traces=Vec::new();
        for length in [0.4001,0.4002] {
            let mut m=model(length,paired);
            for (side,length) in [(ValveSide::Upstream,length),(ValveSide::Downstream,0.4001)] {
                assert_eq!(m.source(side).graph.sections[0].length_m,length);
                assert!((m.viscothermal_sections(side)[0].represented_length_m-length).abs()<1e-12);
                assert_eq!(m.fractional_sections(side).unwrap().1.len(),32);
            }
            let mut trace=Vec::new();let mut work=0.0;let mut loss=0.0;
            for k in 0..512 {
                let f=m.step_under_gate(TwoSidedDrive {source_pressure_pa:if k<128 {100.0} else {0.0},
                    blade_forces_n:[0.0;2]},&gate).unwrap();
                assert!(f.balance_residual_j().abs()<1e-12);
                assert!((f.upstream.inlet_flow_m3_s+f.downstream.inlet_flow_m3_s).abs()<1e-14);
                work+=f.source_work_j;loss+=f.dissipated_energy_j;
                trace.push((f.aperture.state.opening_m,
                    m.node_frame(ValveSide::Downstream,1).unwrap().pressure_pa));
            }
            assert!(loss>0.0);
            assert!((m.stored_energy_j()+loss-work).abs()<1e-11);
            traces.push(trace);
        }
        assert_ne!(traces[0],traces[1],"upstream geometry must reach mechanics and downstream sound");
        assert!(traces[0].iter().zip(&traces[1]).any(|(a,b)|a.0!=b.0));
        assert!(traces[0].iter().zip(&traces[1]).any(|(a,b)|a.1!=b.1));
    }
}

#[test]
fn rejected_fractional_loss_samples_preserve_both_acoustic_memories_and_the_valve() {
    let gate=CancelGate::new_clock_free();
    let mut m=model(0.4001,true);
    let drive=TwoSidedDrive {source_pressure_pa:100.0,blade_forces_n:[0.0;2]};
    for _ in 0..128 {m.step_under_gate(drive,&gate).unwrap();}
    let energy=m.stored_energy_j();let accepted=m.valve().accepted_steps();
    let frames=[ValveSide::Upstream,ValveSide::Downstream]
        .map(|side|*m.node_frame(side,1).unwrap());
    assert!(m.step_under_gate(TwoSidedDrive {source_pressure_pa:f64::NAN,..drive},&gate).is_err());
    gate.request();
    assert!(m.step_under_gate(drive,&gate).is_err());
    assert_eq!(m.stored_energy_j(),energy);
    assert_eq!(m.valve().accepted_steps(),accepted);
    assert_eq!(frames,[ValveSide::Upstream,ValveSide::Downstream]
        .map(|side|*m.node_frame(side,1).unwrap()));
}
