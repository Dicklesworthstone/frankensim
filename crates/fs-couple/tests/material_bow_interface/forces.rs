//! Physical force programs reuse the same source-bound strings, twist and body.
use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::performance::{MultiBowPerformance,
    force::{StringForceCoupling, StringForcePort, StringForceProgram}};
use fs_couple::bowed_string::runtime::reciprocal::plate::MaterialBodyCavity;
use fs_couple::cavity_modal::{CavityNeckLoss, HelmholtzModalCavity};
use fs_couple::render::schedule::force::{ForceInitialization, ForceRenderConfig, ModalForceEvent};
use fs_couple::render::RenderError;
use fs_scenario::HelmholtzCavity;

fn released() -> FrictionDrive { FrictionDrive { speed_m_s: 0.0, normal_force_n: 0.0 } }
fn point(coupling: StringForceCoupling, station: f64) -> StringForcePort {
    StringForcePort { location: StringBowLocation { string: 0, station_fraction: station }, coupling }
}
fn event(sample: u64, port: usize, force_n: f64) -> ModalForceEvent {
    ModalForceEvent { sample, voice: 0, port, force_n }
}
fn forces(initialization: ForceInitialization) -> StringForceProgram {
    StringForceProgram {
        ports: vec![point(StringForceCoupling::Surface, 0.2), point(StringForceCoupling::Centreline, 0.31)],
        initial_force_n: vec![0.002, 0.001],
        events: vec![event(64,1,0.0), event(32,0,0.0)], initialization,
        limits: ForceRenderConfig { sample_rate_hz:48_000, max_block:128,
            max_events:64, max_controls:1024, max_projection_terms:32_768 },
    }
}
fn bank(cx:&Cx<'_>, gate:&CancelGate, cavity:bool) -> MultiBowStringBank {
    let binding=interface(0.4,10.0);
    let base=if cavity {
        let gas=GasState::try_new(&GasSpec::dry_air_ussa1976(),293.15,101_325.0).unwrap();
        let cavity=MaterialBodyCavity {
            cavity:HelmholtzModalCavity::new(HelmholtzCavity {volume_m3:0.003,neck_radius_m:0.01,neck_length_m:0.02,
                inner_end:fs_scenario::NeckEnd::Flanged,outer_end:fs_scenario::NeckEnd::Flanged},
                gas,CavityNeckLoss::Resistance(100.0)).unwrap(),
            modal_budget:ModalAcousticTimeBudget::audible_reference(),listener_m:1.0,
        };
        ReciprocalStringBank::from_material_plate_with_cavity(cx,vec![source()],gas,body(),cavity,
            binding.card().friction_law().clone(),options(),gate).unwrap().with_bow_interface(binding).unwrap()
    } else {build(cx,gate,binding)};
    base.with_additional_bows(vec![],MultiFrictionConfig{max_ports:4,max_sweeps:128,max_setup_terms:16_384},gate).unwrap()
}
fn schedule(drive:FrictionDrive) -> GestureSchedule {
    GestureSchedule::try_new(700,vec![GestureTrack{id:"bow".into(),target:GestureTarget::BowStroke{string:0},
        initial:GestureValue::Bow {velocity_m_per_s:drive.speed_m_s,normal_force_n:drive.normal_force_n,station:0.2},
        events:vec![]}]).unwrap()
}
fn performance(b:MultiBowStringBank,d:FrictionDrive)->MultiBowPerformance {
    b.into_scheduled(schedule(d),vec!["bow".into()],
        MultiBowPerformanceConfig{samples:192,max_block:128,max_work:10_000,max_events:64}).unwrap()
}
fn render(mut p:MultiBowPerformance,chunk:usize)->(MultiBowPerformance,Vec<f64>) {
    let mut out=vec![0.0;192];for block in out.chunks_mut(chunk){p.block(block).unwrap();}(p,out)
}

#[test]
fn physical_port_has_exact_transverse_and_surface_torque_work() {
    with_cx(|cx,gate| {
        let b=bank(cx,gate,true);
        let ports=[point(StringForceCoupling::Centreline,0.2),point(StringForceCoupling::Surface,0.2)];
        let columns=b.string_force_columns(&ports,1024,gate).unwrap();
        let card=&b.cards()[0]; let modes=b.system().components()[0].modes().len();
        assert_eq!(modes,7); assert_eq!(&columns[1][..modes],b.system().port(0).unwrap().left.shapes.as_slice());
        for n in 1..=4 {
            let independent=(n as f64*core::f64::consts::PI*0.2).sin()
                /(card.linear_density_kg_m*card.length_m*0.5).sqrt();
            assert!((columns[0][n-1]-independent).abs()<1e-12);
        }
        assert!(columns[0][4..].iter().all(|x|*x==0.0));
        assert!(columns[1][4..modes].iter().any(|x|*x!=0.0));
        assert!(columns[1][modes..].iter().all(|x|*x==0.0));
        let dq:Vec<_>=(0..b.system().mode_count()).map(|i|1e-6*(i+1) as f64).collect();
        let f=[0.002,-0.001];
        let virtual_work:f64=columns.iter().zip(f).map(|(col,f)|f*col.iter().zip(&dq).map(|(p,q)|p*q).sum::<f64>()).sum();
        let modal_work:f64=(0..dq.len()).map(|k|(columns[0][k]*f[0]+columns[1][k]*f[1])*dq[k]).sum();
        assert!((virtual_work-modal_work).abs()<1e-18);
    });
}

#[test]
fn settled_pluck_moves_the_same_body_and_cavity_then_rings_on_release() {
    with_cx(|cx,gate| {
        let mut program=forces(ForceInitialization::StaticPreload);
        // A single settled load: no excitation remains after the release.
        program.initial_force_n[1]=0.0;program.events=vec![event(32,0,0.0)];
        let mut p=performance(bank(cx,gate,true),released()).with_string_forces(program,gate).unwrap();
        assert_eq!(p.samples_rendered(),0);
        let e=p.preload_energy_j().unwrap();assert!(e>0.0);
        let components=p.bank().system().components();
        assert!(components.iter().flat_map(|m|m.states()).all(|s|s.velocity_m_sqrt_kg_per_s==0.0));
        assert!(components[1].states().iter().any(|s|s.displacement_m_sqrt_kg!=0.0));
        assert!(components[2].states().iter().any(|s|s.displacement_m_sqrt_kg!=0.0));
        assert!(p.bank().source().bank().material_cavity().is_some());
        assert_eq!(p.bank().bow_interface(0),Some(&interface(0.4,10.0)));
        let mut hold=vec![0.0;32];p.block(&mut hold).unwrap();
        assert!(hold.iter().all(|x|x.abs()<1e-8),"force-balanced preload emitted pressure: {hold:?}");
        let mut after=vec![0.0;128];p.block(&mut after).unwrap();
        assert!(after.iter().any(|x|x.abs()>1e-10));
        let final_energy=p.bank().system().total_energy_j().unwrap();
        assert!(final_energy>0.0 && final_energy<=e+1e-12);
        let frame=p.bank().system().last_frame().unwrap();
        assert_eq!(frame.external_work_j,0.0); assert_eq!(frame.drive_work_j,0.0);
        assert!(frame.energy_residual_j.abs()<=frame.energy_tolerance_j);
    });
}

#[test]
fn actual_load_history_matches_direct_joint_physics_at_every_callback_partition() {
    with_cx(|cx,gate| {
        let d=FrictionDrive{speed_m_s:0.01,normal_force_n:0.001};
        let mut program=forces(ForceInitialization::RetainState);
        // Input order at equal time: the last assignment wins, the other port persists.
        program.events.extend([event(10,0,0.02),event(10,0,-0.003),event(10,1,0.004)]);
        let mut direct=bank(cx,gate,false);
        let cols=direct.string_force_columns(&program.ports,1024,gate).unwrap();
        let mut expected=Vec::new();
        for sample in 0..192 {
            let (f0,f1)=if sample<10{(0.002,0.001)}else if sample<32{(-0.003,0.004)}
                else if sample<64{(0.0,0.004)}else{(0.0,0.0)};
            let external:Vec<_>=(0..direct.system().mode_count()).map(|i|cols[0][i]*f0+cols[1][i]*f1).collect();
            let frame=direct.step(&external,&[d]).unwrap();
            assert!(frame.energy_residual_j.abs()<=frame.energy_tolerance_j);
            expected.push(frame.observer_pressure_pa.to_bits());
        }
        for partition in [1,37,128] {
            let p=performance(bank(cx,gate,false),d).with_string_forces(program.clone(),gate).unwrap();
            assert!(p.preload_energy_j().is_none());
            assert_eq!(p.string_force_program().unwrap().events,program.events);
            let (_,out)=render(p,partition);
            assert_eq!(out.iter().map(|x|x.to_bits()).collect::<Vec<_>>(),expected);
        }
    });
}

#[test]
fn changing_only_force_location_changes_twist_and_real_pressure_without_retuning() {
    with_cx(|cx,gate| {
        let a=bank(cx,gate,false);let b=bank(cx,gate,false);
        for (a,b) in a.system().components().iter().zip(b.system().components()){assert_eq!(a.modes(),b.modes());}
        let mut pa=forces(ForceInitialization::StaticPreload);pa.initial_force_n[1]=0.0;pa.events=vec![event(0,0,0.0)];
        let mut pb=pa.clone();pb.ports[0].location.station_fraction=0.31;
        let (a,ap)=render(performance(a,released()).with_string_forces(pa,gate).unwrap(),37);
        let (b,bp)=render(performance(b,released()).with_string_forces(pb,gate).unwrap(),128);
        assert!(ap.iter().zip(&bp).any(|(a,b)|(a-b).abs()>1e-10));
        assert_ne!(a.bank().system().components()[0].states(),b.bank().system().components()[0].states());
        assert!(a.bank().system().components()[0].states()[4..].iter().any(|s|s.velocity_m_sqrt_kg_per_s!=0.0));
        assert_eq!(a.bank().source().bank().gas(),b.bank().source().bank().gas());
    });
}

#[test]
fn no_load_program_keeps_legacy_bits_and_retained_motion_is_never_reset() {
    with_cx(|cx,gate| {
        let d=FrictionDrive{speed_m_s:0.01,normal_force_n:0.001};
        let mut original=bank(cx,gate,false);let mut loaded=bank(cx,gate,false);
        let zero=vec![0.0;original.system().mode_count()];
        for _ in 0..19 {original.step(&zero,&[d]).unwrap();loaded.step(&zero,&[d]).unwrap();}
        let mut program=forces(ForceInitialization::RetainState);program.initial_force_n.fill(0.0);program.events.clear();
        let loaded=performance(loaded,d).with_string_forces(program,gate).unwrap();
        assert_eq!(loaded.bank().system().samples_rendered(),19);
        assert_eq!(loaded.samples_rendered(),0);
        let (_,a)=render(performance(original,d),37);let(_,b)=render(loaded,128);assert_eq!(a,b);
    });
}

#[test]
fn bad_force_programs_refuse_before_preload_and_cannot_replace_a_live_binding() {
    with_cx(|cx,gate| {
        for which in 0..9 {
            let mut program=forces(ForceInitialization::StaticPreload);
            match which {
                0=>program.ports[0].location.station_fraction=1.0,
                1=>program.ports[0].location.string=1,
                2=>program.events.push(event(192,0,0.0)),
                3=>program.events[0].force_n=f64::NAN,
                4=>program.limits.max_projection_terms=0,
                5=>program.limits.max_controls=0,
                6=>program.limits.sample_rate_hz=44_100,
                7=>program.initial_force_n.pop().map(|_|()).unwrap(),
                _=>program.events[0].voice=1,
            }
            assert!(performance(bank(cx,gate,false),released()).with_string_forces(program,gate).is_err());
        }
        let p=performance(bank(cx,gate,false),released()).with_string_forces(forces(ForceInitialization::RetainState),gate).unwrap();
        assert!(p.with_string_forces(forces(ForceInitialization::RetainState),gate).is_err());
        let mut p=performance(bank(cx,gate,false),released());p.block(&mut[0.0]).unwrap();
        assert!(p.with_string_forces(forces(ForceInitialization::RetainState),gate).is_err());
        let stop=CancelGate::new_clock_free();stop.request();
        assert!(matches!(performance(bank(cx,gate,false),released()).with_string_forces(forces(ForceInitialization::StaticPreload),&stop),
            Err(RenderError::Coupled(_))));
        let mut moving=bank(cx,gate,false);let zero=vec![0.0;moving.system().mode_count()];moving.step(&zero,&[drivers()[0]]).unwrap();
        assert!(performance(moving,released()).with_string_forces(forces(ForceInitialization::StaticPreload),gate).is_err());
    });
}

#[test]
fn material_density_reaches_force_maps_and_plucked_sound_at_fixed_load_and_geometry() {
    with_cx(|cx,gate| {
        let make=|rho| {
            let binding=interface(0.4,10.0);
            let gas=GasState::try_new(&GasSpec::dry_air_ussa1976(),293.15,101_325.0).unwrap();
            ReciprocalStringBank::from_material_plate(cx,vec![source_at_density(rho)],gas,body(),
                binding.card().friction_law().clone(),options(),gate).unwrap().with_bow_interface(binding).unwrap()
                .with_additional_bows(vec![],MultiFrictionConfig{max_ports:4,max_sweeps:128,max_setup_terms:16_384},gate).unwrap()
        };
        let a=make(7800.0);let b=make(15600.0);
        let ports=[point(StringForceCoupling::Surface,0.2)];
        let ca=a.string_force_columns(&ports,1024,gate).unwrap();let cb=b.string_force_columns(&ports,1024,gate).unwrap();
        for i in 0..7 { assert!((cb[0][i]/ca[0][i]-0.5_f64.sqrt()).abs()<1e-12); }
        assert_eq!(a.system().components()[1].modes(),b.system().components()[1].modes());
        let mut p=forces(ForceInitialization::StaticPreload);p.initial_force_n[1]=0.0;p.events=vec![event(0,0,0.0)];
        let pa=performance(a,released()).with_string_forces(p.clone(),gate).unwrap();
        let pb=performance(b,released()).with_string_forces(p,gate).unwrap();
        // Density changes inertia, not the same-force static physical stiffness.
        let ea=pa.preload_energy_j().unwrap();let eb=pb.preload_energy_j().unwrap();
        assert!((ea-eb).abs()<1e-9*ea.max(eb));
        let(_,oa)=render(pa,37);let(_,ob)=render(pb,37);
        assert!(oa.iter().zip(&ob).any(|(a,b)|(a-b).abs()>1e-12));
    });
}

#[test]
fn refused_physical_force_sample_keeps_joint_state_and_poisons_only_the_callback_wrapper() {
    with_cx(|cx,gate| {
        let mut program=forces(ForceInitialization::RetainState);
        program.initial_force_n.fill(0.0);program.events=vec![event(1,0,1e6)];
        let mut p=performance(bank(cx,gate,false),drivers()[0]).with_string_forces(program,gate).unwrap();
        p.block(&mut[0.0]).unwrap();
        let before:Vec<_>=p.bank().system().components().iter().map(|m|m.states().to_vec()).collect();
        let frame=p.bank().system().last_frame().unwrap().clone();
        let mut out=[123.0;2];assert!(p.block(&mut out).is_err());
        assert_eq!(out,[123.0;2]);assert_eq!(p.samples_rendered(),1);
        assert_eq!(p.bank().system().last_frame(),Some(&frame));
        for (m,q) in p.bank().system().components().iter().zip(before){assert_eq!(m.states(),q);}
        assert!(matches!(p.validate_sample_count(0),Err(RenderError::Poisoned)));
        assert!(matches!(p.block(&mut out),Err(RenderError::Poisoned)));
    });
}
