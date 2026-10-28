//! Material string + normal contact + original friction + body pressure.
use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::polarization::{TransverseDirection,TransverseDirections};
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::contact_load::BowContactLoad;
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::stops::{MaterialStringStop,StringStopSide};
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::performance::motion::{
    StopMotionKnot,StopMotionProgram,ObstacleMotionLimits,
};
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::performance::bow_stations::{
    BowStationKnot,BowStationProgram,BowStationMotionLimits,
};
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::performance::stations::{
    StopStationKnot,StopStationProgram,StopStationMotionLimits,
};
use fs_couple::render::schedule::force::coupled::contact::ModalContactConfig;
use fs_couple::render::schedule::force::coupled::contact::friction::multiple::obstacles::FrictionObstacleConfig;
use fs_couple::render::schedule::force::coupled::contact::friction::ModalFrictionRegime;
use fs_dcontact::Obstacle;
use fs_scenario::gesture::GestureEvent;

fn binding(mu: f64, maximum_pressure: f64) -> DryFrictionBinding {
    DryFrictionBinding::new(DryInterfaceSystemCard::new(InterfaceSystemRef::new(
        "synthetic/string-surface", "synthetic/history", "synthetic/contact", InputAuthority::SyntheticFixture,
        InterfaceMedium::Dry).unwrap(), FrictionLaw::Coulomb {static_mu:0.8,kinetic_mu:mu},
        DryFrictionApplicability::new(Range::new(280.0,320.0).unwrap(),
            Range::new(0.0,maximum_pressure).unwrap(),Range::new(0.0,10.0).unwrap()).unwrap()).unwrap(),
        293.15,0.0001).unwrap()
}
fn contact() -> MaterialStringStop {
    MaterialStringStop { location:options().bow,
        side:StringStopSide::Along(TransverseDirection::try_new([0.0,-1.0]).unwrap()),body:None,
        law:Obstacle::new(vec![-1.0],1,1,vec![-1e-7],vec![1.0],5000.0,1.0,
            "synthetic normal compliance, not measured hair/rosin".into()).unwrap(),
        limits:ModalContactConfig {max_iterations:96,maximum_force_n:10.0,maximum_penetration_m:0.001,
            force_absolute_tolerance_n:1e-11,force_relative_tolerance:1e-9} }
}
fn raw(cx:&Cx<'_>,gate:&CancelGate,rho:f64,interface:DryFrictionBinding,normal:MaterialStringStop)->MultiBowStringBank {
    let source=source_at_density(rho).with_transverse_polarizations(
        TransverseDirections::try_new([1.0,0.0],[0.6,0.8]).unwrap()).unwrap();
    ReciprocalStringBank::from_material_plate(cx,vec![source],
        GasState::try_new(&GasSpec::dry_air_ussa1976(),293.15,101325.0).unwrap(),body(),
        interface.card().friction_law().clone(),options(),gate).unwrap().with_bow_interface(interface).unwrap()
        .with_additional_bows(vec![],MultiFrictionConfig {max_ports:4,max_sweeps:128,max_setup_terms:100000},gate).unwrap()
        .with_string_stops(vec![normal],FrictionObstacleConfig {max_contacts:4,max_setup_terms:100000},gate).unwrap()
}
fn loaded(cx:&Cx<'_>,gate:&CancelGate,rho:f64,interface:DryFrictionBinding)->MultiBowStringBank {
    raw(cx,gate,rho,interface,contact()).with_bow_contact_loads(vec![BowContactLoad {bow:0,stop:0}],gate).unwrap()
}

#[test]
fn joint_motion_admission_retains_the_original_contact_binding_order() {
    with_cx(|cx,gate| {
        let mut moving=loaded(cx,gate,7800.0,binding(0.4,1e6))
            .with_bow_station_motion(bow_motion().limits,gate).unwrap();
        assert_eq!(moving.bow_station_fractions(),Some([0.2].as_slice()));
        assert_eq!(moving.stop_station_fractions(),Some([0.2].as_slice()));
        let mut fixed=loaded(cx,gate,7800.0,binding(0.4,1e6));
        let zero=vec![0.0;fixed.system().mode_count()];
        for _ in 0..16 {
            assert_eq!(moving.step_with_bow_stations(&zero,&[drive()],&[0.2],None).unwrap(),
                fixed.step(&zero,&[drive()]).unwrap());
        }
        let moving=raw(cx,gate,7800.0,binding(0.4,1e6),contact())
            .with_bow_station_motion(bow_motion().limits,gate).unwrap();
        assert!(moving.with_bow_contact_loads(vec![BowContactLoad {bow:0,stop:0}],gate).is_err());
    });
}

fn drive()->FrictionDrive { FrictionDrive {speed_m_s:0.3,normal_force_n:0.0} }
fn motion()->StopMotionProgram {
    StopMotionProgram {knots:vec![StopMotionKnot {sample:64,offsets_m:vec![-3e-6]},
        StopMotionKnot {sample:128,offsets_m:vec![-3e-6]},
        StopMotionKnot {sample:256,offsets_m:vec![3e-5]}],
        limits:ObstacleMotionLimits {maximum_abs_offset_m:0.001,maximum_speed_m_s:0.1},
        max_knots:8,max_coordinate_values:64}
}
fn offset(sample:u64)->f64 {
    if sample<=64 {-3e-6*(sample as f64/64.0)}
    else if sample<=128 {-3e-6}
    else if sample<256 {-3e-6+(3e-5-(-3e-6))*((sample-128) as f64/128.0)}
    else {3e-5}
}
fn schedule(normal:f64)->GestureSchedule {
    GestureSchedule::try_new(1000,vec![GestureTrack {id:"bow".into(),target:GestureTarget::BowStroke {string:0},
        initial:GestureValue::Bow {velocity_m_per_s:0.3,normal_force_n:normal,station:0.2},events:vec![]}]).unwrap()
}
fn performance_config()->MultiBowPerformanceConfig {
    MultiBowPerformanceConfig {samples:320,max_block:64,max_work:10000,max_events:64}
}
fn bow_motion()->BowStationProgram {
    BowStationProgram {knots:vec![
        BowStationKnot {sample:64,station_fractions:vec![0.202]},
        BowStationKnot {sample:128,station_fractions:vec![0.201]}],
        limits:BowStationMotionLimits {maximum_speed_m_s:2.0,maximum_step_m:0.00003,max_update_terms:1_048_576},
        max_knots:8,max_coordinate_values:64}
}
fn bow_station(sample:u64)->f64 {
    if sample<=64 {0.2+(0.202-0.2)*(sample as f64/64.0)}
    else if sample<128 {0.202+(0.201-0.202)*((sample-64) as f64/64.0)}
    else {0.201}
}
fn displacement(bank:&MultiBowStringBank)->Vec<f64> {
    bank.system().components().iter().flat_map(|m|m.states().iter()
        .map(|s|s.displacement_m_sqrt_kg)).collect()
}

#[test]
fn moving_loaded_bow_keeps_the_material_basis_and_accounts_for_normal_patch_work() {
    with_cx(|cx,gate| {
        let mut bank=loaded(cx,gate,7800.0,binding(0.4,1e6))
            .with_stop_motion(motion().limits,gate).unwrap()
            .with_bow_station_motion(bow_motion().limits,gate).unwrap();
        let modes:Vec<_>=bank.system().components().iter().map(|m|m.modes().to_vec()).collect();
        let zero=vec![0.0;bank.system().mode_count()];
        let initial=bank.system().total_energy_j().unwrap();
        let (mut supply,mut loss,mut tolerance,mut old_offset)=(0.0,0.0,0.0,0.0);
        let (mut migration_work,mut sound,mut released)=(0.0,0.0_f64,false);
        for sample in 1..=320 {
            let old_q=displacement(&bank);
            let old_normal=bank.system().obstacle_normal_column(0).unwrap().to_vec();
            let next_offset=offset(sample);
            let station=bow_station(sample);
            let frame=bank.step_with_bow_stations(&zero,&[drive()],&[station],Some(&[next_offset]))
                .unwrap().clone();
            let q=displacement(&bank);
            let normal=bank.system().obstacle_normal_column(0).unwrap();
            let contact=bank.system().last_obstacle_frames().unwrap()[0];
            // Independent discrete product rule: moving the normal patch is an
            // actuator, not a free change of modal shape or stored potential.
            let work=contact.normal_force_n*normal.iter().zip(&old_normal).zip(q.iter().zip(&old_q))
                .map(|((end,start),(q,q0))|(end-start)*0.5*(q+q0)).sum::<f64>();
            let actuator_work=work-contact.normal_force_n*(next_offset-old_offset);
            assert!((contact.attachment_work_j-work).abs()<1e-15);
            assert!((frame.obstacle_work_j-actuator_work).abs()<1e-15);
            assert_eq!(frame.ports[0].normal_force_n,contact.normal_force_n);
            assert_eq!(bank.bow_station_fractions(),Some([station].as_slice()));
            assert_eq!(bank.stop_station_fractions(),Some([station].as_slice()));
            if contact.normal_force_n==0.0 {released=true;assert_eq!(frame.ports[0].traction_n,0.0);}
            supply+=frame.external_work_j+frame.drive_work_j+actuator_work;
            loss+=frame.network_dissipation_j+frame.friction_dissipation_j+contact.dissipation_j;
            tolerance+=frame.energy_tolerance_j;
            assert!((bank.system().total_energy_j().unwrap()-initial+loss-supply).abs()<=tolerance);
            assert!(frame.energy_residual_j.abs()<=frame.energy_tolerance_j);
            migration_work+=work.abs();sound=sound.max(frame.observer_pressure_pa.abs());
            old_offset=next_offset;
        }
        assert!(migration_work>0.0 && sound>1e-12 && released);
        for (model,original) in bank.system().components().iter().zip(&modes) {assert_eq!(model.modes(),original);}
        assert_eq!(bank.source().bows()[0].location.station_fraction,0.2);
        assert_eq!(bank.source().stops()[0].location.station_fraction,0.2);
        // Settling the held tangential midpoint to its accepted endpoint must
        // not move the normal again or create another geometry-work pulse.
        assert_eq!(bank.step(&zero,&[drive()]).unwrap().obstacle_work_j,0.0);
    });
}

#[test]
fn refused_joint_patch_motion_keeps_both_maps_histories_and_the_next_hold() {
    with_cx(|cx,gate| {
        let build=||loaded(cx,gate,7800.0,binding(0.4,1e6))
            .with_bow_station_motion(bow_motion().limits,gate).unwrap();
        let mut actual=build();let mut reference=build();
        let zero=vec![0.0;actual.system().mode_count()];
        actual.step_with_bow_stations(&zero,&[drive()],&[0.20001],None).unwrap();
        reference.step_with_bow_stations(&zero,&[drive()],&[0.20001],None).unwrap();
        let states:Vec<_>=actual.system().components().iter().map(|m|m.states().to_vec()).collect();
        let frame=actual.system().last_frame().unwrap().clone();
        let contacts=actual.system().last_obstacle_frames().unwrap().to_vec();
        let tangent=actual.system().port(0).unwrap().clone();
        let normal=actual.system().obstacle_normal_column(0).unwrap().to_vec();
        // This fails inside the coupled solve, after new tangential maps have
        // been staged: a resolved normal cannot also accept a prescribed load.
        let invalid_drive=FrictionDrive {normal_force_n:1.0,..drive()};
        assert!(actual.step_with_bow_stations(&zero,&[invalid_drive],&[0.20002],None).is_err());
        assert!(actual.step_with_bow_and_stop_stations_under_gate(&zero,&[drive()],&[0.20002],
            Some(&[0.20003]),None,gate).is_err());
        assert!(actual.step_with_stop_stations(&zero,&[drive()],&[0.20002],None).is_err());
        let cancelled=CancelGate::new_clock_free();cancelled.request();
        assert!(actual.step_with_bow_stations_under_gate(&zero,&[drive()],&[0.20002],None,&cancelled).is_err());
        assert_eq!(actual.system().last_frame(),Some(&frame));
        assert_eq!(actual.system().last_obstacle_frames().unwrap(),contacts);
        assert_eq!(actual.system().port(0).unwrap().left,tangent.left);
        assert_eq!(actual.system().port(0).unwrap().right,tangent.right);
        assert_eq!(actual.system().obstacle_normal_column(0).unwrap(),normal);
        assert_eq!(actual.bow_station_fractions(),Some([0.20001].as_slice()));
        assert_eq!(actual.stop_station_fractions(),Some([0.20001].as_slice()));
        for (model,state) in actual.system().components().iter().zip(&states) {assert_eq!(model.states(),state);}
        for _ in 0..16 {
            assert_eq!(actual.step(&zero,&[drive()]).unwrap(),reference.step(&zero,&[drive()]).unwrap());
        }
    });
}

#[test]
fn joint_bow_and_normal_paths_share_the_original_performance_clock() {
    with_cx(|cx,gate| {
        let mut direct=loaded(cx,gate,7800.0,binding(0.4,1e6))
            .with_stop_motion(motion().limits,gate).unwrap()
            .with_bow_station_motion(bow_motion().limits,gate).unwrap();
        let zero=vec![0.0;direct.system().mode_count()];
        let expected:Vec<_>=(1..=320).map(|sample|direct.step_with_bow_stations(&zero,&[drive()],
            &[bow_station(sample)],Some(&[offset(sample)])).unwrap().observer_pressure_pa.to_bits()).collect();
        for separate_stop_path in [false,true] {
            for partition in [1,37,64] {
                let mut performance=loaded(cx,gate,7800.0,binding(0.4,1e6))
                    .into_scheduled(schedule(0.0),vec!["bow".into()],performance_config()).unwrap()
                    .with_stop_motion(motion(),gate).unwrap();
                if separate_stop_path {
                    performance=performance.with_stop_station_motion(StopStationProgram {
                        knots:vec![StopStationKnot {sample:128,station_fractions:vec![0.2]}],
                        limits:StopStationMotionLimits {maximum_speed_m_s:2.0,maximum_step_m:0.00003,
                            max_update_terms:1_048_576},max_knots:8,max_coordinate_values:64,
                    },gate).unwrap();
                }
                let mut performance=performance.with_bow_station_motion(bow_motion(),gate).unwrap();
                let mut output=vec![0.0;320];
                for block in output.chunks_mut(partition) {performance.block(block).unwrap();}
                assert_eq!(output.iter().map(|p|p.to_bits()).collect::<Vec<_>>(),expected);
                assert_eq!(performance.bank().stop_station_fractions(),Some([0.201].as_slice()));
            }
        }
        let slower_normals=loaded(cx,gate,7800.0,binding(0.4,1e6))
            .into_scheduled(schedule(0.0),vec!["bow".into()],performance_config()).unwrap()
            .with_stop_station_motion(StopStationProgram {
                knots:vec![StopStationKnot {sample:128,station_fractions:vec![0.2]}],
                limits:StopStationMotionLimits {maximum_speed_m_s:0.1,maximum_step_m:0.00003,
                    max_update_terms:1_048_576},max_knots:8,max_coordinate_values:64,
            },gate).unwrap();
        assert!(slower_normals.with_bow_station_motion(bow_motion(),gate).is_err());
    });
}

#[test]
fn material_contact_motion_sets_bow_load_and_releases_without_a_force_command() {
    with_cx(|cx,gate| {
        let interface=binding(0.4,1e6);
        let mut bank=loaded(cx,gate,7800.0,interface.clone()).with_stop_motion(motion().limits,gate).unwrap();
        assert_eq!(bank.source().bow_contact_loads(),&[BowContactLoad {bow:0,stop:0}]);
        assert_eq!(bank.bow_interface(0),Some(&interface));
        let tangent=&bank.system().port(0).unwrap().left.shapes;
        let normal=&bank.system().obstacle(0).unwrap().left.shapes;
        assert!(tangent[..4].iter().any(|v|*v!=0.0));
        assert!(tangent[4..8].iter().all(|v|*v==0.0));
        assert!(tangent[8..].iter().any(|v|*v!=0.0));
        assert!(normal[..4].iter().all(|v|*v==0.0));
        assert!(normal[4..8].iter().any(|v|*v!=0.0));
        assert!(normal[8..].iter().all(|v|*v==0.0));
        let zero=vec![0.0;bank.system().mode_count()];
        let (mut loaded_steps,mut released_steps,mut sound)=(0,0,0.0_f64);
        for sample in 1..=320 {
            let f=bank.step_with_stop_offsets(&zero,&[drive()],&[offset(sample)]).unwrap().clone();
            let n=bank.system().last_obstacle_frames().unwrap()[0].normal_force_n;
            assert_eq!(f.ports[0].normal_force_n,n);
            if n==0.0 {released_steps+=1;assert_eq!(f.ports[0].traction_n,0.0);
                assert_eq!(f.ports[0].regime,ModalFrictionRegime::Released);}
            else {loaded_steps+=1;}
            sound=sound.max(f.observer_pressure_pa.abs());
            assert!(f.energy_residual_j.abs()<=f.energy_tolerance_j);
        }
        assert!(loaded_steps>0 && released_steps>0 && sound>1e-12);
        assert!(bank.system().total_energy_j().unwrap()>0.0);
    });
}

#[test]
fn material_and_interface_changes_reach_resolved_load_and_actual_body_sound() {
    with_cx(|cx,gate| {
        let make=|rho,mu|loaded(cx,gate,rho,binding(mu,1e6)).with_stop_motion(motion().limits,gate).unwrap()
            .with_bow_station_motion(bow_motion().limits,gate).unwrap();
        let mut a=make(7800.0,0.4);let mut b=make(15600.0,0.4);let mut c=make(7800.0,0.2);
        for (a,c) in a.system().components().iter().zip(c.system().components()) {assert_eq!(a.modes(),c.modes());}
        assert_eq!(a.system().components()[1].modes(),b.system().components()[1].modes());
        let an=&a.system().obstacle(0).unwrap().left.shapes;
        let bn=&b.system().obstacle(0).unwrap().left.shapes;
        for (x,y) in an[4..8].iter().zip(&bn[4..8]) {assert!((y/x-0.5_f64.sqrt()).abs()<1e-12);}
        let zero=vec![0.0;a.system().mode_count()];
        let (mut material_sound,mut interface_sound,mut normal_change)=(0.0,0.0,0.0);
        for sample in 1..=128 {
            let station=[bow_station(sample)];let offsets=[offset(sample)];
            let fa=a.step_with_bow_stations(&zero,&[drive()],&station,Some(&offsets)).unwrap().clone();
            let fb=b.step_with_bow_stations(&zero,&[drive()],&station,Some(&offsets)).unwrap().clone();
            let fc=c.step_with_bow_stations(&zero,&[drive()],&station,Some(&offsets)).unwrap().clone();
            material_sound+=(fa.observer_pressure_pa-fb.observer_pressure_pa).abs();
            interface_sound+=(fa.observer_pressure_pa-fc.observer_pressure_pa).abs();
            normal_change+=(fa.ports[0].normal_force_n-fb.ports[0].normal_force_n).abs();
        }
        assert!(material_sound>1e-12 && interface_sound>1e-12 && normal_change>1e-10);
    });
}

#[test]
fn original_schedule_and_support_sampler_match_direct_contact_loaded_pressure() {
    with_cx(|cx,gate| {
        let mut direct=loaded(cx,gate,7800.0,binding(0.4,1e6)).with_stop_motion(motion().limits,gate).unwrap();
        let zero=vec![0.0;direct.system().mode_count()];
        let expected:Vec<_>=(1..=320).map(|sample|direct.step_with_stop_offsets(&zero,&[drive()],&[offset(sample)])
            .unwrap().observer_pressure_pa.to_bits()).collect();
        for partition in [1,37,64] {
            let mut performance=loaded(cx,gate,7800.0,binding(0.4,1e6))
                .into_scheduled(schedule(0.0),vec!["bow".into()],performance_config()).unwrap()
                .with_stop_motion(motion(),gate).unwrap();
            let mut pressure=vec![0.0;320];for block in pressure.chunks_mut(partition){performance.block(block).unwrap();}
            assert_eq!(pressure.iter().map(|v|v.to_bits()).collect::<Vec<_>>(),expected);
            assert_eq!(performance.bank().source().bow_contact_loads(),&[BowContactLoad {bow:0,stop:0}]);
            assert_eq!(performance.bank().bow_interface(0),Some(&binding(0.4,1e6)));
        }
        for late in [false,true] {
            let mut s=schedule(if late {0.0} else {0.0001});
            if late {s=GestureSchedule::try_new(1000,vec![GestureTrack {id:"bow".into(),target:GestureTarget::BowStroke {string:0},
                initial:GestureValue::Bow {velocity_m_per_s:0.3,normal_force_n:0.0,station:0.2},
                events:vec![GestureEvent {time_s:0.001,transition_s:0.0,value:GestureValue::Bow {
                    velocity_m_per_s:0.3,normal_force_n:0.0001,station:0.2}}]}]).unwrap();}
            assert!(loaded(cx,gate,7800.0,binding(0.4,1e6)).into_scheduled(s,vec!["bow".into()],performance_config()).is_err());
        }
    });
}

#[test]
fn source_domain_checks_resolved_pressure_and_refusal_keeps_the_whole_instrument() {
    with_cx(|cx,gate| {
        // The authored drive is zero-normal; the actual overlap still has load.
        let mut bank=loaded(cx,gate,7800.0,binding(0.4,0.001));
        let zero=vec![0.0;bank.system().mode_count()];
        let before:Vec<_>=bank.system().components().iter().map(|m|m.states().to_vec()).collect();
        assert!(bank.step(&zero,&[drive()]).is_err());
        assert_eq!(bank.system().samples_rendered(),0);assert!(bank.system().last_frame().is_none());
        for (m,s) in bank.system().components().iter().zip(&before){assert_eq!(m.states(),s);}
        assert!(bank.system().last_obstacle_frames().is_none());
        // A genuinely separated contact bypasses its positive-load applicability.
        let mut separated=contact();separated.law=Obstacle::new(vec![-1.0],1,1,vec![0.001],vec![1.0],5000.0,1.0,
            "synthetic separated surface".into()).unwrap();
        let mut separated=raw(cx,gate,7800.0,binding(0.4,0.001),separated)
            .with_bow_contact_loads(vec![BowContactLoad {bow:0,stop:0}],gate).unwrap();
        assert_eq!(separated.step(&zero,&[drive()]).unwrap().ports[0].normal_force_n,0.0);
    });
}

#[test]
fn geometry_binding_cannot_relabel_distant_oblique_or_opposite_torque_contacts() {
    with_cx(|cx,gate| {
        for variant in 0..4 {
            let mut stop=contact();
            match variant {
                0=>stop.location.station_fraction=0.31,
                1=>stop.side=StringStopSide::Positive,
                2=>stop.side=StringStopSide::Along(TransverseDirection::try_new([0.6,-0.8]).unwrap()),
                _=>stop.side=StringStopSide::Along(TransverseDirection::try_new([0.0,1.0]).unwrap()),
            }
            assert!(raw(cx,gate,7800.0,binding(0.4,1e6),stop)
                .with_bow_contact_loads(vec![BowContactLoad {bow:0,stop:0}],gate).is_err());
        }
        let b=BowContactLoad {bow:0,stop:0};
        assert!(loaded(cx,gate,7800.0,binding(0.4,1e6)).with_bow_contact_loads(vec![b],gate).is_err());
        assert!(raw(cx,gate,7800.0,binding(0.4,1e6),contact())
            .with_bow_contact_loads(vec![b,b],gate).is_err());
    });
}

#[test]
fn resolved_tangent_and_normal_cannot_name_different_receiving_surfaces() {
    use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::surfaces::BowBodySurface;
    with_cx(|cx,gate| {
        for same_body in [false,true] {
            let source=source_at_density(7800.0).with_transverse_polarizations(
                TransverseDirections::try_new([1.0,0.0],[0.6,0.8]).unwrap()).unwrap();
            let bank=ReciprocalStringBank::from_material_plate(cx,vec![source],
                GasState::try_new(&GasSpec::dry_air_ussa1976(),293.15,101325.0).unwrap(),body(),
                binding(0.4,1e6).card().friction_law().clone(),options(),gate).unwrap()
                .with_additional_bows(vec![],MultiFrictionConfig {max_ports:4,max_sweeps:128,max_setup_terms:100000},gate).unwrap();
            let footprint=bank.source().bank().strings()[0].connections[0].body.clone();
            let bank=bank.with_bow_body_surfaces(vec![BowBodySurface {bow:0,body:footprint.clone()}],gate).unwrap();
            let mut stop=contact();
            if same_body {stop.body=Some(footprint);}
            let bank=bank.with_string_stops(vec![stop],
                FrictionObstacleConfig {max_contacts:4,max_setup_terms:100000},gate).unwrap();
            let result=bank.with_bow_contact_loads(vec![BowContactLoad {bow:0,stop:0}],gate);
            assert_eq!(result.is_ok(),same_body);
        }
    });
}
