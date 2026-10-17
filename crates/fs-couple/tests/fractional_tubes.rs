#![allow(unused_imports)]
//! Tests exercise the actual valve/network owners, not a stand-in output filter.
use fs_couple::bernoulli_aperture::dynamic::{ApertureState, DynamicAperture, DynamicApertureSpec};
use fs_couple::bernoulli_aperture::network::{NetworkNode, TubeNetworkSpec, TubeSection};
use fs_couple::bernoulli_aperture::paired::{ApertureBlade, BladeState, PairedAperture,
    PairedApertureLimits, PairedApertureSpec};
use fs_couple::bernoulli_aperture::two_sided::{ApertureDuctSpec, LoadedAperture,
    TwoSidedAperture, TwoSidedDrive, TwoSidedFrame, TwoSidedLimits, ValveSide, combined_inlet_state};
use fs_couple::bernoulli_aperture::{BernoulliAperture, dynamic::ApertureTerminal};
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;
use fs_material::gas::{GasSpec, GasState};

const DT: f64 = 1.0 / 48_000.0;
fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
}
fn duct(g: GasState, delay: usize, reflection: f64) -> ApertureDuctSpec {
    ApertureDuctSpec {
        graph: TubeNetworkSpec {
            nodes: vec![NetworkNode::Inlet, NetworkNode::Termination { reflection }],
            sections: vec![TubeSection { nodes: [0, 1], length_m: delay as f64 * g.sound_speed * DT,
                radius_m: 0.007, max_length_error_m: 1e-12 }],
            sound_speed_m_s: g.sound_speed, max_wave_memory_bytes: 1 << 20,
        }, gases: vec![g], viscothermal: vec![],
    }
}
fn contact() -> Obstacle {
    Obstacle::new(vec![-1.0], 1, 1, vec![0.0], vec![1.0], 1e8, 2.0,
        "numerical two-sided valve contact fixture".into()).unwrap()
}
fn limits() -> TwoSidedLimits {
    TwoSidedLimits { maximum_energy_j: 1.0,
        energy_absolute_tolerance_j: 1e-16, energy_relative_tolerance: 1e-7,
        flow_absolute_tolerance_m3_s: 1e-16, flow_relative_tolerance: 1e-8 }
}
fn blade() -> ApertureBlade {
    ApertureBlade { mass_kg: 1e-5, stiffness_n_m: 500.0,
        damping_n_s_m: 0.03, pressure_area_m2: 1.25e-4 }
}
// About 24900 rad/s of pole rate at 48 kHz is outside the midpoint band.
fn paired(z: f64, g: GasState, steps: u64, second_mass_ratio: f64) -> PairedAperture {
    let mut b = [blade(); 2];
    b[1].mass_kg *= second_mass_ratio;
    PairedAperture::new_with_propagator(PairedApertureSpec {
        rest_opening_m: 4e-4, width_m: 0.013, blades: b,
        density_kg_m3: g.density, impedance_pa_s_m3: z, time_step_s: DT, max_steps: steps,
        limits: PairedApertureLimits {
            maximum_abs_displacement_m: 0.01, maximum_abs_velocity_m_s: 100.0,
            maximum_energy_j: 1.0, max_evaluations: 256,
            velocity_tolerance_m_s: 1e-13, relative_tolerance: 1e-11,
        },
    }, [BladeState::default(); 2], contact(),
        fs_couple::bernoulli_aperture::paired::BladePropagator::ExactZoh).unwrap()
}
fn drive(p: f64) -> TwoSidedDrive {
    TwoSidedDrive { source_pressure_pa: p, blade_forces_n: [0.0; 2] }
}

use fs_couple::bernoulli_aperture::fractional::FractionalTubePolicy;
fn policy() -> FractionalTubePolicy {
    FractionalTubePolicy {maximum_frequency_hz:3000.0,maximum_summed_phase_error_rad:0.01}
}
fn fractional_system(extra: f64) -> TwoSidedAperture {
    let g=gas(); let mut sides=[duct(g,12,0.7),duct(g,32,-0.8)];
    sides[0].graph.sections[0].length_m+=(0.2+extra)*g.sound_speed*DT;
    sides[1].graph.sections[0].length_m+=0.3*g.sound_speed*DT;
    let (_,z)=combined_inlet_state(&sides).unwrap();
    TwoSidedAperture::with_fractional_transits(paired(z,g,1024,1.7).into(),sides,
        limits(),[policy();2],&CancelGate::new_clock_free()).unwrap()
}

#[test]
fn subcell_geometry_reaches_actual_reed_feedback_without_rounding_the_length() {
    let mut a=fractional_system(0.0); let mut b=fractional_system(0.1);
    let gate=CancelGate::new_clock_free();let mut different=false;
    for n in 0..1024 {
        let x=drive(if n<256 {20.0} else {0.0});
        let fa=a.step_under_gate(x,&gate).unwrap();let fb=b.step_under_gate(x,&gate).unwrap();
        assert!(fa.balance_residual_j().abs()<1e-13);
        assert!(fb.balance_residual_j().abs()<1e-13);
        if n<24 {assert_eq!(fa.aperture.state,fb.aperture.state);}
        else {different|=(fa.aperture.state.opening_m-fb.aperture.state.opening_m).abs()>1e-13;}
    }
    assert!(different,"these two lengths formerly had the SAME rounded transit");
    assert!(a.stored_energy_j()>0.0,"release must retain physical vibration");
    let (_,sections)=a.fractional_sections(ValveSide::Upstream).unwrap();
    assert_eq!(sections[0].integer_samples,12);
    assert!((sections[0].transit.samples()-0.2).abs()<1e-13);
    assert!((sections[0].realized_dc_length_m-sections[0].requested_length_m).abs()<1e-12);
}

#[test]
fn gas_temperature_changes_the_original_load_and_subsample_transit() {
    let a=gas();let b=GasState::try_new(&GasSpec::dry_air_ussa1976(),293.25,101325.0).unwrap();
    let mut spec=duct(a,42,0.0).graph;
    spec.sections[0].length_m+=0.3*a.sound_speed*DT;
    let (_,ra,da)=spec.realize_fractional_waveguide(&[a],DT,policy(),&CancelGate::new_clock_free()).unwrap();
    spec.sound_speed_m_s=b.sound_speed;
    let (_,rb,db)=spec.realize_fractional_waveguide(&[b],DT,policy(),&CancelGate::new_clock_free()).unwrap();
    assert_eq!(da[0].integer_samples,db[0].integer_samples);
    assert_ne!(da[0].transit,db[0].transit);
    assert_ne!(ra[0].impedance_pa_s_m3,rb[0].impedance_pa_s_m3);
    assert!((db[0].requested_samples-spec.sections[0].length_m/(b.sound_speed*DT)).abs()<1e-12);
}

#[test]
fn geometry_phase_and_work_budgets_refuse_instead_of_silently_retuning() {
    let g=gas();let mut graph=duct(g,12,0.7).graph;
    graph.sections[0].length_m+=0.37*g.sound_speed*DT;
    let gate=CancelGate::new_clock_free();
    assert!(graph.realize_waveguide(&[g],DT,&gate).is_err(),"old integer length does not fit this geometry tolerance");
    assert!(graph.realize_fractional_waveguide(&[g],DT,policy(),&gate).is_ok());
    let mut p=policy();p.maximum_summed_phase_error_rad=0.0;
    assert!(graph.realize_fractional_waveguide(&[g],DT,p,&gate).is_err());
    p=policy();p.maximum_frequency_hz=24000.0;
    assert!(graph.realize_fractional_waveguide(&[g],DT,p,&gate).is_err());
    graph.max_wave_memory_bytes=1;
    assert!(graph.realize_fractional_waveguide(&[g],DT,policy(),&gate).is_err());
    graph.max_wave_memory_bytes=1<<20;graph.sections[0].length_m=0.9*g.sound_speed*DT;
    assert!(graph.realize_fractional_waveguide(&[g],DT,policy(),&gate).is_err());
    gate.request();assert!(graph.realize_fractional_waveguide(&[g],DT,policy(),&gate).is_err());
}
