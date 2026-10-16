//! Synthetic reactive-cavity consumer; independent geometry and coupled twins.
use fs_couple::bernoulli_aperture::BernoulliAperture;
use fs_couple::bernoulli_aperture::cavity::HelmholtzLoadSpec;
use fs_couple::bernoulli_aperture::dynamic::{ApertureState, ApertureTerminal, DynamicAperture, DynamicApertureSpec};
use fs_couple::bernoulli_aperture::network::{ApertureNetwork, ApertureNetworkFrame, NetworkNode, TubeNetworkSpec, TubeSection};
use fs_couple::bernoulli_aperture::tube::TubeDrive;
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;

fn cavity(volume: f64) -> HelmholtzLoadSpec {
    HelmholtzLoadSpec { volume_m3: volume, neck_radius_m: 0.003,
        effective_neck_length_m: 0.02, resistance_pa_s_m3: 2e5 }
}
fn model(volume: f64, budget: u64, state: ApertureState) -> ApertureNetwork {
    let dt=1e-5;
    let section=|a,b,n,radius| TubeSection { nodes:[a,b],length_m:f64::from(n)*(343.0*dt),
        radius_m:radius,max_length_error_m:1e-15 };
    let spec=TubeNetworkSpec {
        nodes:vec![NetworkNode::Inlet,NetworkNode::Junction,NetworkNode::Termination { reflection:-0.8 },
            cavity(volume).termination(1.2,343.0).unwrap()],
        sections:vec![section(0,1,8,0.007),section(1,2,12,0.009),section(1,3,5,0.003)],
        sound_speed_m_s:343.0,max_wave_memory_bytes:1<<20,
    };
    let mechanics=DynamicApertureSpec {
        aperture:BernoulliAperture { rest_opening_m:4e-4,width_m:0.013,closing_pressure_pa:6000.0 },
        mass_kg:1e-5,stiffness_n_m:500.0,damping_ratio:0.35,density_kg_m3:1.2,
        impedance_pa_s_m3:spec.inlet_impedance(1.2).unwrap(),time_step_s:dt,max_steps:budget,
    };
    let lay=Obstacle::new(vec![-1.0],1,1,vec![0.0],vec![1.0],1e8,2.0,
        "synthetic contact; not identified material data".into()).unwrap().with_internal_loss(5.0).unwrap();
    ApertureNetwork::new(DynamicAperture::new(mechanics,state,lay).unwrap(),spec).unwrap()
}
fn rest() -> ApertureState { ApertureState { opening_m:4e-4,opening_velocity_m_s:0.0 } }

#[test]
fn cavity_coefficients_and_resonance_come_from_geometry_and_fluid() {
    let c=cavity(1e-4);
    let a=c.impedance(1.2,343.0).unwrap();
    let area=core::f64::consts::PI*c.neck_radius_m.powi(2);
    assert!((a.inertance_pa_s2_m3-1.2*c.effective_neck_length_m/area).abs()<1e-10);
    assert!((a.compliance_m3_pa.unwrap()-c.volume_m3/(1.2*343.0*343.0)).abs()<1e-23);
    let expected=343.0/(2.0*core::f64::consts::PI)*(area/(c.effective_neck_length_m*c.volume_m3)).sqrt();
    assert!((c.resonance_hz(1.2,343.0).unwrap()-expected).abs()<1e-10);
    assert!((cavity(4e-4).resonance_hz(1.2,343.0).unwrap()-0.5*expected).abs()<1e-10);
    let dense=c.impedance(2.4,343.0).unwrap();
    assert!((dense.inertance_pa_s2_m3-2.0*a.inertance_pa_s2_m3).abs()<1e-10);
    assert!((dense.compliance_m3_pa.unwrap()-0.5*a.compliance_m3_pa.unwrap()).abs()<1e-23);
}

#[test]
fn driven_contact_and_cavity_storage_close_one_energy_window() {
    let mut omitted_storage_error=0.0_f64;
    for volume in [1e-4,4e-4] {
        for state in [rest(),ApertureState { opening_m:-1e-4,opening_velocity_m_s:-0.2 }] {
            let mut net=model(volume,2048,state);
            let mut saw_storage=false;
            for n in 0..2048 {
                let before=net.stored_energy_j();
                let old_load=net.node_frame(3).unwrap().stored_energy_j;
                let f=net.step(TubeDrive { upstream_pressure_pa:if n<512 {1200.0} else {0.0},
                    body_flow_m3_s:if n<800 {2e-7*(0.1*f64::from(n)).sin()} else {0.0} }).unwrap();
                let scale=(before+f.stored_energy_j+f.dissipated_energy_j+f.upstream_work_j.abs()+f.body_work_j.abs())
                    .max(f64::MIN_POSITIVE);
                assert!(f.balance_residual_j().abs()<3e-10*scale);
                assert!(f.dissipated_energy_j>=0.0);
                let load=net.node_frame(3).unwrap();
                saw_storage|=load.stored_energy_j>0.0;
                omitted_storage_error=omitted_storage_error.max((f.balance_residual_j()
                    -(load.stored_energy_j-old_load)).abs());
                if n>=800 { assert!(f.stored_energy_j<=before+3e-10*scale); }
            }
            assert!(saw_storage);
        }
    }
    assert!(omitted_storage_error>1e-10,"wave storage alone must not pass this balance");
}

#[test]
fn changing_cavity_volume_changes_the_valve_only_after_the_return_path() {
    let (mut a,mut b)=(model(1e-4,1024,rest()),model(4e-4,1024,rest()));
    let (mut dp,mut dy)=(0.0_f64,0.0_f64);
    for n in 0..1024 {
        let drive=TubeDrive { upstream_pressure_pa:800.0,body_flow_m3_s:0.0 };
        let (x,y)=(a.step(drive).unwrap(),b.step(drive).unwrap());
        if n<26 {
            assert_eq!(x.aperture.bore_pressure_pa.to_bits(),y.aperture.bore_pressure_pa.to_bits());
            assert_eq!(x.aperture.state.opening_m.to_bits(),y.aperture.state.opening_m.to_bits());
        }
        dp=dp.max((x.aperture.bore_pressure_pa-y.aperture.bore_pressure_pa).abs());
        dy=dy.max((x.aperture.state.opening_m-y.aperture.state.opening_m).abs());
    }
    assert!(dp>1.0 && dy>1e-9,"cavity must affect real feedback, not only an independent oscillator");
}

#[test]
fn cancellation_and_budget_resume_retain_reactive_energy_and_state() {
    let drives:Vec<_>=(0..512).map(|n| TubeDrive {
        upstream_pressure_pa:if n<256 {800.0} else {0.0},body_flow_m3_s:2e-7*(0.1*f64::from(n)).sin(),
    }).collect();
    let mut uninterrupted=model(1e-4,512,rest());
    let mut expected=vec![ApertureNetworkFrame::default();512];
    let gate=CancelGate::new_clock_free();
    uninterrupted.advance_block(&drives,&mut expected,&gate).unwrap();
    let mut resumed=model(1e-4,137,rest());
    let mut actual=vec![ApertureNetworkFrame::default();512];
    resumed.advance_block(&drives[..100],&mut actual[..100],&gate).unwrap();
    let before=resumed.stored_energy_j();
    let cancel=CancelGate::new_clock_free();cancel.request();
    let progress=resumed.advance_block(&drives[100..],&mut actual[100..],&cancel).unwrap();
    assert_eq!(progress.completed,0);assert_eq!(progress.terminal,ApertureTerminal::Cancelled);
    assert_eq!(resumed.stored_energy_j(),before);
    assert_eq!(actual[100..],vec![ApertureNetworkFrame::default();412]);
    let progress=resumed.advance_block(&drives[100..],&mut actual[100..],&gate).unwrap();
    assert_eq!(progress.completed,37);assert_eq!(progress.terminal,ApertureTerminal::BudgetExhausted);
    assert_eq!(actual[137..],vec![ApertureNetworkFrame::default();375]);
    resumed.extend_step_budget(512).unwrap();
    resumed.advance_block(&drives[137..],&mut actual[137..],&gate).unwrap();
    assert_eq!(actual,expected);
    assert_eq!(resumed.node_frame(3),uninterrupted.node_frame(3));
}

#[test]
fn invalid_samples_and_memoryless_controls_cannot_clear_cavity_storage() {
    let (mut a,mut b)=(model(1e-4,256,rest()),model(1e-4,256,rest()));
    let drive=TubeDrive { upstream_pressure_pa:800.0,body_flow_m3_s:2e-7 };
    for _ in 0..100 { a.step(drive).unwrap();b.step(drive).unwrap(); }
    let before=*a.node_frame(3).unwrap();
    assert!(before.stored_energy_j>0.0);
    assert!(a.set_terminal_reflection(3,0.0).is_err());
    for bad in [f64::NAN,f64::INFINITY,f64::MAX] {
        assert!(a.step(TubeDrive {body_flow_m3_s:bad,..drive}).is_err());
        assert_eq!(*a.node_frame(3).unwrap(),before);
        assert_eq!(a.aperture().accepted_steps(),100);
    }
    for _ in 0..128 { assert_eq!(a.step(drive).unwrap(),b.step(drive).unwrap()); }
}

#[test]
fn invalid_cavity_geometry_or_fluid_is_refused_without_substitution() {
    for field in 0..4 {
        let mut c=cavity(1e-4);
        match field {0=>c.volume_m3=0.0,1=>c.neck_radius_m=f64::NAN,
            2=>c.effective_neck_length_m=-1.0,_=>c.resistance_pa_s_m3=-1.0}
        assert!(c.termination(1.2,343.0).is_err());
    }
    for value in [0.0,-1.0,f64::NAN,f64::INFINITY] {
        assert!(cavity(1e-4).termination(value,343.0).is_err());
        assert!(cavity(1e-4).termination(1.2,value).is_err());
    }
}
