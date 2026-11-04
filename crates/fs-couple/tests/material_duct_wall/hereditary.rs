//! Original material spectrum -> shared pressure memory -> actual jet/outlet PCM.
use super::*;
use fs_couple::bernoulli_aperture::wall::material::{MaterialWallInput, WallRelaxation};
use fs_material::visco::GeneralizedMaxwell;
use fs_math::c64::C64;
use fs_couple::flue_loop::{FlueGeometry, FlueVoice};
use fs_couple::flue_loop::network::FlueNetworkSpec;
use fs_couple::flue_loop::network::performance::{FlueOutletConfig, FlueOutletPerformance};
use fs_couple::pcm_wav::baffled::CircularOutletReceiver;
use fs_couple::pcm_wav::observation::PressureRenderer;
use fs_scenario::gesture::{GestureSchedule, GestureTrack, GestureTarget, GestureValue};

fn spectrum(scale: f64) -> WallRelaxation {
    WallRelaxation::new(GeneralizedMaxwell::new(1e6,vec![(8e5,0.003*scale),(4e5,0.0005*scale)]).unwrap(),
        "synthetic homogeneous viscoelastic wall at 293.15 K; not measured cane or metal".into(),
        BAND,1.0,0.8).unwrap()
}
fn specimen(scale: f64) -> MaterialWallAssignment {
    MaterialWallAssignment { loss:WallMaterialLoss::Maxwell(spectrum(scale)), ..assignment(1000.0,1e6,2.0) }
}

#[test]
fn source_spectrum_maps_exact_modulus_and_geometry_without_reidentification() {
    let a=specimen(1.0); let w=bound(a.clone());
    assert_eq!(w.source(),&a);
    assert_eq!(w.load().compliance_terms().len(),2);
    assert!(w.load().terms().is_empty());
    assert_eq!(w.patch().wall.resistance,0.0); // The unused KV property is NOT added twice.
    let WallMaterialLoss::Maxwell(source)=&a.loss else {unreachable!()};
    assert!(source.provenance().contains("synthetic"));
    for f in [20.0,73.0,400.0,1000.0,4000.0] {
        let omega=core::f64::consts::TAU*f; let s=C64::new(0.0,omega);
        let base=w.load().base();
        let mut z=s.scale(base.inertance_pa_s2_m3)+s.scale(base.compliance_m3_pa.unwrap()).recip();
        for t in w.load().compliance_terms() {
            z=z+(s+C64::from_re(t.rate_per_s)).recip().scale(t.resistance_pa_s_m3*t.rate_per_s);
        }
        // The original material owns E*(omega), including both storage and loss.
        let (storage,loss)=source.material().modulus(omega);
        let factor=w.patch().wall.stiffness_per_area/source.material().e_inf/w.patch().area_m2;
        let expected=C64::new(factor*loss/omega,base.inertance_pa_s2_m3*omega-factor*storage/omega);
        assert!((z-expected).abs()<3e-12*expected.abs());
    }
    let changed=bound(specimen(3.0));
    assert_eq!(changed.patch(),w.patch()); // Only the history times changed, not equilibrium mass/stiffness.
    assert_ne!(changed.load(),w.load());
}

#[test]
fn empty_or_zero_material_arms_reproduce_the_elastic_wall_exactly() {
    let mut a=assignment(1000.0,1e6,2.0); a.loss=WallMaterialLoss::Elastic;
    let elastic=bound(a.clone());
    for terms in [vec![],vec![(0.0,0.001)]] {
        a.loss=WallMaterialLoss::Maxwell(WallRelaxation::new(GeneralizedMaxwell::new(1e6,terms).unwrap(),
            "explicit zero-memory control".into(),BAND,1.0,0.8).unwrap());
        let zero=bound(a.clone());
        assert_eq!(elastic.shunt(),zero.shunt()); assert_eq!(elastic.patch(),zero.patch());
        zero.validate_time_step(1.0/f64::from(RATE)).unwrap();
    }
}

#[test]
fn equilibrium_source_band_and_memory_resolution_cannot_be_silently_changed() {
    let gate=CancelGate::new_clock_free();
    let wrong=|model,band| {
        let mut a=specimen(1.0);
        a.loss=WallMaterialLoss::Maxwell(WallRelaxation::new(model,"supplied source".into(),band,1.0,0.8).unwrap());
        assert!(bind_material_walls(graph(),vec![a],BAND,1,&gate).is_err());
    };
    wrong(GeneralizedMaxwell::new(2e6,vec![(8e5,0.001)]).unwrap(),BAND);
    wrong(GeneralizedMaxwell::new(1e6,vec![(8e5,0.001)]).unwrap(),[100.0,2000.0]);
    for (e,tau) in [(-1.0,0.001),(8e5,0.0),(8e5,f64::NAN)] {
        assert!(WallRelaxation::new(GeneralizedMaxwell{e_inf:1e6,terms:vec![(e,tau)]},
            "invalid source".into(),BAND,1.0,0.8).is_err());
    }
    assert!(WallRelaxation::new(GeneralizedMaxwell::new(1e6,vec![(8e5,0.001);9]).unwrap(),
        "over-budget source".into(),BAND,1.0,0.8).is_err());
    assert!(WallRelaxation::new(GeneralizedMaxwell::new(1e6,vec![]).unwrap(),
        "".into(),BAND,1.0,0.8).is_err());
    let w=bound(specimen(1.0));
    w.validate_time_step(1.0/f64::from(RATE)).unwrap();
    assert!(w.validate_time_step(0.01).is_err()); // Unresolved material pole.
    assert!(w.validate_time_step(0.0005).is_err()); // Poles admitted, instantaneous wall frequency unresolved.
}

#[test]
fn hereditary_wall_closes_physical_momentum_and_all_spring_energy_in_the_graph() {
    let dt=1.0/f64::from(RATE); let gate=CancelGate::new_clock_free();
    let (g,w)=bind_material_walls(graph(),vec![specimen(1.0)],BAND,1,&gate).unwrap();
    let w=&w[0]; w.validate_time_step(dt).unwrap();
    let (mut n,_)=g.realize_waveguide(&[gas(),gas()],dt,&gate).unwrap();
    let mut old=w.observe_waveguide(&n).unwrap(); let mut stored=false;
    for i in 0..2048 {
        let old_pressures=n.load_compliance_pressures(1).unwrap().to_vec();
        let frame=n.step(if i<20 {1.0} else {0.0}).unwrap();
        let new=w.observe_waveguide(&n).unwrap();
        let energy=w.relaxation_stored_energy_j(&n).unwrap();
        let node=n.node_frame(1).unwrap();
        close(new.stored_energy_j+energy,node.stored_energy_j);
        let inertia=w.patch().wall.surface_density*(new.velocity_m_s-old.velocity_m_s)/dt;
        let spring=w.patch().wall.stiffness_per_area*0.5*(old.displacement_m+new.displacement_m);
        let pressures=n.load_compliance_pressures(1).unwrap();
        let mut memory_pressure=0.0;
        let mut loss=0.0;
        for (j,term) in w.load().compliance_terms().iter().enumerate() {
            let mid=0.5*(old_pressures[j]+pressures[j]); memory_pressure+=mid;
            let change=(pressures[j]-old_pressures[j])/dt;
            let rhs=term.rate_per_s*(term.resistance_pa_s_m3*node.load_flow_m3_s-mid);
            assert!((change-rhs).abs()<2e-9*(change.abs()+rhs.abs()+term.rate_per_s*mid.abs()).max(1e-30));
            loss+=mid*mid/term.resistance_pa_s_m3*dt;
        }
        let scale=(inertia.abs()+spring.abs()+memory_pressure.abs()+node.pressure_pa.abs()).max(1e-30);
        assert!((inertia+spring+memory_pressure-node.pressure_pa).abs()<2e-10*scale);
        close(loss,node.absorbed_energy_j);
        close(w.patch().area_m2*0.5*(old.velocity_m_s+new.velocity_m_s),node.load_flow_m3_s);
        assert!(frame.balance_residual_j().abs()<1e-18);
        stored|=energy>0.0;
        old=new;
    }
    assert!(stored);
    let (rigid,_)=graph().realize_waveguide(&[gas(),gas()],dt,&gate).unwrap();
    assert!(w.observe_waveguide(&rigid).is_err());
}

fn voice(scale: f64, rate: u32) -> Result<FlueVoice, fs_couple::flue_loop::FlueError> {
    let card=fs_aeroac::jetcard::mint_tonal_interim_card().unwrap();
    let source=FlueNetworkSpec { graph:graph(),gases:vec![gas(),
        GasState::try_new(&GasSpec::dry_air_ussa1976(),320.0,gas().pressure).unwrap()],
        viscothermal:vec![],minimum_jet_speed_m_s:1.0,maximum_jet_mach:0.3,max_jet_history_samples:4096 };
    let geometry=FlueGeometry {flue_width_m:0.01,flue_height_m:0.001,cut_up_m:0.004,
        labium_offset_m:0.0,jet_angle_rad:0.0,bore_area_m2:core::f64::consts::PI*0.006*0.006,mouth_loss_rel:0.02};
    FlueVoice::from_material_network(Some(&card),90.0,1e-6,geometry,source,
        MaterialWallInput {assignments:vec![specimen(scale)],band_hz:BAND,max_walls:1},
        rate,3,0,&CancelGate::new_clock_free())
}

#[test]
fn material_memory_changes_actual_microphone_pressure_and_survives_callbacks() {
    let mut signals=Vec::new();
    for scale in [1.0,3.0] {
        let mut partitioned=None;
        for block in [1,37,512] {
            let schedule=GestureSchedule::try_new(700,vec![GestureTrack{id:"blow".into(),
                target:GestureTarget::BlowingPressure,initial:GestureValue::PressurePa(300.0),events:vec![]}]).unwrap();
            let config=FlueOutletConfig {sample_rate_hz:RATE,samples:1024,max_block:512,max_controls:64,
                max_compile_work:10000,outlet_node:2,receiver:CircularOutletReceiver{
                    position_m:[0.0,0.0,1.0],radial_rings:4,angular_points:16,maximum_frequency_hz:1000.0},
                maximum_abs_pressure_pa:1e6};
            let mut p=FlueOutletPerformance::new(voice(scale,RATE).unwrap(),schedule,config).unwrap();
            let mut output=vec![0.0;1024];
            for chunk in output.chunks_mut(block) {p.block(chunk).unwrap();}
            let n=p.voice().network().unwrap();
            assert_eq!(n.material_walls()[0].source(),&specimen(scale));
            let memory=n.wall_relaxation_energy_j(1).unwrap(); assert!(memory>0.0);
            close(memory+n.wall_state(1).unwrap().stored_energy_j,n.node_frame(1).unwrap().stored_energy_j);
            if let Some(previous)=&partitioned {assert_eq!(&output,previous);} else {partitioned=Some(output);}
        }
        let output=partitioned.unwrap(); assert!(output.iter().any(|p|p.abs()>1e-12)); signals.push(output);
    }
    assert_ne!(signals[0],signals[1]);
    assert!(voice(0.001,RATE).is_err()); // Source construction enforces its actual material-memory clock.
}

/// G3: changing only a sourced relaxation time changes the native DG pressure
/// and material memory. The equilibrium mass/stiffness and geometry stay fixed.
#[test]
fn sourced_relaxation_reaches_native_dg_motion_pressure_and_loss() {
    use fs_couple::bernoulli_aperture::wall::material::LocalMaterialWall;
    use fs_duct::dg::{DgBoundary, DgDuct, DgMesh};
    let radius = 0.006;
    let geometry = fs_duct::Duct { segments: vec![fs_duct::Segment::Cylinder {
        radius, length: 0.3,
    }] };
    let mut histories = Vec::new();
    let mut base_pin = None;
    for tau in [0.0005, 0.002] {
        let mut values = properties(1000.0, 1e6, 2.0).to_vec();
        values.push(("wall_arm_modulus", Pressure::DIMS, 8e5));
        values.push(("wall_arm_time", Dims([0, 0, 1, 0, 0, 0]), tau));
        let source = state(&values, None, None);
        let spectrum = WallRelaxation::from_state_point(
            &source,
            &[("wall_arm_modulus".into(), "wall_arm_time".into())],
            BAND, 1.0, 0.8,
        ).unwrap();
        assert_eq!(spectrum.source_identity(), Some(source.identity()));
        let wall = LocalMaterialWall::new(
            source.clone(), 0.0002, AxialConstraint::StressFree,
            WallMaterialLoss::Maxwell(spectrum), BAND,
        ).unwrap();
        let dg = DgDuct::new(
            &geometry, &gas(), DgMesh { elements: 2, order: 2 },
            DgBoundary::Closed, DgBoundary::Closed,
        ).unwrap();
        let bound = wall.bind_dg(dg).unwrap();
        assert_eq!(bound.wall().material().identity(), source.identity());
        let dg = bound.duct();
        let pin = dg.compliant_wall().unwrap();
        assert_eq!(pin.resistance, 0.0); // The unused KV property is not added.
        if let Some(previous) = base_pin { assert_eq!(pin, previous); }
        base_pin = Some(pin);
        assert_eq!(dg.total_dofs(), 5 * dg.field_dofs());
        let propagator = dg.propagator(1e-5, 1e-13).unwrap();
        let mut current = dg.project_initial_state(|_| 1.0, |_| 0.0);
        let initial = current.stored_energy_j();
        let mut previous = initial;
        let mut trace = Vec::new();
        let mut maximum_memory = 0.0_f64;
        for _ in 0..100 {
            propagator.step(&mut current, &[]).unwrap();
            propagator.sync_nodal(&mut current).unwrap();
            let motion = dg.compliant_wall_state(&current).unwrap();
            assert_eq!(motion.relaxation_pressure_pa.len(), 1);
            close(motion.stored_energy_j, current.wall_stored_energy_j());
            close(motion.dissipated_power_w, dg.wall_dissipated_power_w(&current).unwrap());
            assert!(current.stored_energy_j() <= previous + 1e-12 * initial);
            previous = current.stored_energy_j();
            maximum_memory = maximum_memory.max(motion.relaxation_stored_energy_j);
            trace.push(current.pressure()[0]);
        }
        assert!(maximum_memory > 0.0);
        assert!(current.stored_energy_j() < initial);
        histories.push(trace);
    }
    assert_ne!(histories[0], histories[1]);
}
