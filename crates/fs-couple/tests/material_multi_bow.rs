//! Synthetic material/geometry -> shared strings/body/cavity -> joint bow sound.
//! These source coefficients are wiring fixtures, not measured instrument data.
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::{MaterialBowContact, MultiBowStringBank, performance::MultiBowPerformanceConfig};
use fs_couple::bowed_string::runtime::reciprocal::{bank::{MaterialStringSource, ReciprocalStringBank,
    StringBankOptions, StringBowLocation, plate::{SharedMaterialPlateBody, SharedPlateBridge}}, plate::MaterialBodyCavity};
use fs_couple::cavity_modal::{CavityNeckLoss, HelmholtzModalCavity};
use fs_couple::modal_acoustic_time::{ModalAcousticState, ModalAcousticTimeBudget};
use fs_couple::pcm_wav::observation::PressureRenderer;
use fs_couple::render::{RenderError, schedule::force::coupled::{ModalCouplingConfig,
    contact::friction::{FrictionDrive, ModalFrictionConfig, multiple::MultiFrictionConfig},
    render::friction::FrictionGestureConfig}};
use fs_couple::string_specimen::{KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY as ETA, with_uniform_circular_material_state};
use fs_couple::thin_plate::{PlateMaterialModel, PlateThicknessConstraint, with_uniform_plate_material_state};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::{gas::{GasSpec, GasState}, state_point::{MaterialPropertySelection,
    ResolvedMaterialStatePoint, ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point}};
use fs_qty::{Density, Dims, DynViscosity, Pressure};
use fs_scenario::{HelmholtzCavity, PrestressedString, ThinPlate,
    gesture::{GestureEvent, GestureSchedule, GestureTarget, GestureTrack, GestureValue}};
use fs_tribo::FrictionLaw;

#[path = "material_multi_bow/fiber_ribbon.rs"]
mod fiber_ribbon;

#[path = "material_multi_bow/shear_layer.rs"]
mod shear_layer;

#[path = "material_multi_bow/shear_thermal.rs"]
mod shear_thermal;

#[path = "material_multi_bow/stations.rs"]
mod stations;

fn with_cx<T>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> T) -> T {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(&gate, arena, StreamKey { seed: 6, kernel_id: 5, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic), &gate)
    })
}
fn material(rho: f64) -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new(); let mut requirements = Vec::new();
    for (key, dims, value) in [("density", Density::DIMS, rho), ("young_modulus", Pressure::DIMS, 200e9),
        ("poisson_ratio", Dims::NONE, 0.3), (ETA, DynViscosity::DIMS, 10.0)] {
        claims.insert_claim(PropertyClaim { key: PropertyKey::new(key,dims), value: PropertyValue::Scalar {value,dims},
            validity: ValidityDomain::unconstrained().with("T",280.0,310.0).with("omega",0.0,100000.0),
            uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![], provenance: Provenance { source:"synthetic multi-bow".into(), license:"CC0-1.0".into(),artifact:None },
        }).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(key,dims,ScalarAdmissibility::Finite).unwrap());
    }
    let card=MaterialCard::assemble(MaterialStateId {chemistry:"synthetic".into(),phase:"solid".into(),process:"fixture".into(),revision:0},claims,vec![]).unwrap();
    resolve_material_state_point(&card,&QueryPoint::new().with("T",293.15).unwrap().with("omega",1000.0).unwrap(),
        &requirements,MaterialPropertySelection::SingleClaimOnly).unwrap()
}
fn source(rho:f64)->MaterialStringSource {
    let string=PrestressedString {length_m:0.5,tension_n:20.0,lin_density_kg_m:1.0,axial_stiffness_n:1.0,
        width_m:1.0,n_modes:4,damping_ratio:0.0,rayleigh:None,bending_stiffness_n_m2:0.0,
        kelvin_voigt_bending:None,relaxation_bending:None,polarization_detune:0.0,moving_end:false,
        tension_loss_factor:None};
    MaterialStringSource::Homogeneous(with_uniform_circular_material_state(string,0.0003,&material(rho))
        .unwrap().with_kelvin_voigt_bending_loss().unwrap())
}
fn law()->FrictionLaw {FrictionLaw::Coulomb {static_mu:0.8,kinetic_mu:0.4}}
fn options()->StringBankOptions {
    StringBankOptions {sample_rate_hz:48000,bow:StringBowLocation {string:0,station_fraction:0.2},
        string_budget:ModalAcousticTimeBudget::audible_reference(),
        network:ModalCouplingConfig {max_modes:32,max_connections:8,max_setup_terms:16384,nyquist_guard_fraction:0.9,
            maximum_total_energy_j:100.0,maximum_abs_pressure_pa:1e6,maximum_abs_connection_force_n:100.0,
            solve_relative_tolerance:1e-10,energy_absolute_tolerance_j:1e-12,energy_relative_tolerance:1e-8},
        friction:ModalFrictionConfig {max_iterations:96,maximum_force_n:100.0,maximum_slip_speed_m_s:10.0,
            force_absolute_tolerance_n:1e-10,force_relative_tolerance:1e-9,sticking_tolerance_m:1e-12}}
}
fn joint_limits()->MultiFrictionConfig {MultiFrictionConfig {max_ports:4,max_sweeps:128,max_setup_terms:100000}}
fn extra()->MaterialBowContact {MaterialBowContact {location:StringBowLocation {string:1,station_fraction:0.3},law:law(),friction:options().friction}}
fn base(cx:&Cx<'_>,gate:&CancelGate,rho:f64)->ReciprocalStringBank {
    base_with_body_density(cx,gate,rho,7800.0)
}
fn base_with_body_density(cx:&Cx<'_>,gate:&CancelGate,rho:f64,body_rho:f64)->ReciprocalStringBank {
    let plate=with_uniform_plate_material_state(ThinPlate {length_m:0.4,width_m:0.3,thickness_m:1.0,density_kg_m3:1.0,
        e1_pa:1.0,e2_pa:1.0,nu12:0.0,g12_pa:1.0,material_angle_rad:0.0,damping_ratio:0.0,thermoelastic:None,
        kelvin_voigt_bending:None,n_modes:2,geometric_nonlinearity:false,pretension_n_m:0.0,clamped:false},
        &material(body_rho),PlateMaterialModel::Isotropic,PlateThicknessConstraint::FixedThickness(0.001))
        .unwrap().with_kelvin_voigt_bending_loss().unwrap();
    let body=SharedMaterialPlateBody {specimen:plate,bridges:vec![SharedPlateBridge {station_fraction:0.8,
        stiffness_n_m:3000.0,damping_n_s_m:0.1,rest_extension_m:0.0};2],
        modal_budget:ModalAcousticTimeBudget::audible_reference(),listener_m:1.0};
    let gas=GasState::try_new(&GasSpec::dry_air_ussa1976(),293.15,101325.0).unwrap();
    let cavity=MaterialBodyCavity {cavity:HelmholtzModalCavity::new(HelmholtzCavity {volume_m3:0.005,
        neck_radius_m:0.006,neck_length_m:0.02,inner_end:fs_scenario::NeckEnd::Flanged,
        outer_end:fs_scenario::NeckEnd::Flanged},gas,CavityNeckLoss::Resistance(1000.0)).unwrap(),
        modal_budget:ModalAcousticTimeBudget::audible_reference(),listener_m:1.0};
    ReciprocalStringBank::from_material_plate_with_cavity(cx,vec![source(7800.0),source(rho)],gas,body,cavity,law(),options(),gate).unwrap()
}
fn joint(cx:&Cx<'_>,gate:&CancelGate,rho:f64)->MultiBowStringBank {
    base(cx,gate,rho).with_additional_bows(vec![extra()],joint_limits(),gate).unwrap()
}
fn states(bank:&MultiBowStringBank)->Vec<ModalAcousticState> {
    bank.system().components().iter().flat_map(|m|m.states().iter().copied()).collect()
}
fn value(speed:f64,load:f64,station:f64)->GestureValue {GestureValue::Bow {velocity_m_per_s:speed,normal_force_n:load,station}}
fn schedule()->GestureSchedule {
    GestureSchedule::try_new(700,(0..2).map(|i| {
        let station=if i==0 {0.2} else {0.3};
        GestureTrack {id:format!("bow-{i}"),target:GestureTarget::BowStroke {string:i},
            initial:value(0.01,0.1,station),events:vec![
                GestureEvent {time_s:0.002,transition_s:0.005,value:value(0.02,0.2,station)},
                GestureEvent {time_s:0.004,transition_s:0.002,value:value(-0.01,0.1,station)},
                GestureEvent {time_s:0.007,transition_s:0.0,value:value(0.0,0.0,station)},
            ]}
    }).collect()).unwrap()
}
fn tracks()->Vec<String> {vec!["bow-0".into(),"bow-1".into()]}
fn performance_config()->MultiBowPerformanceConfig {
    MultiBowPerformanceConfig {samples:513,max_block:128,max_work:10000,max_events:128}
}
fn drive(speed:f64,load:f64)->FrictionDrive {FrictionDrive {speed_m_s:speed,normal_force_n:load}}

fn reacting_surface(mut bank:MultiBowStringBank,gate:&CancelGate)->MultiBowStringBank {
    use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::surfaces::BowBodySurface;
    // A synthetic string rubbing a plate at a retained physical footprint.
    // These are the source-projected rows, not invented bow-hair parameters.
    let surfaces=bank.source().bank().strings().iter().enumerate().map(|(bow,string)|
        BowBodySurface {bow,body:string.connections[0].body.clone()}).collect();
    let before=states(&bank);
    bank=bank.with_bow_body_surfaces(surfaces,gate).unwrap();
    assert_eq!(states(&bank),before);
    bank
}

#[test]
fn a_cavity_neck_cannot_be_relabeled_as_a_solid_friction_receiver() {
    use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::surfaces::BowBodySurface;
    use fs_couple::render::schedule::force::coupled::ModalAttachment;
    with_cx(|cx, gate| {
        let bank = joint(cx, gate, 7800.0);
        assert!(bank.source().bank().material_cavity().is_some());
        let cavity_modes = bank.system().components()[3].modes().len();
        assert!(bank.with_bow_body_surfaces(vec![BowBodySurface {
            bow: 0, body: ModalAttachment { component: 1, shapes: vec![1.0; cavity_modes] },
        }], gate).is_err());
    });
}

#[test]
fn receiving_material_changes_contact_motion_and_shared_body_pressure() {
    with_cx(|cx,gate| {
        let make=|density| reacting_surface(base_with_body_density(cx,gate,7800.0,density)
            .with_additional_bows(vec![extra()],joint_limits(),gate).unwrap(),gate);
        let mut a=make(7800.0);let mut b=make(15600.0);
        assert_eq!(a.source().body_surfaces().len(),2);
        for port in 0..2 {
            let pa=a.system().port(port).unwrap();let pb=b.system().port(port).unwrap();
            assert_eq!(pa.left,pb.left);
            let ra=pa.right.as_ref().unwrap();let rb=pb.right.as_ref().unwrap();
            assert_eq!(ra.component,2);assert_eq!(rb.component,2);
            assert_eq!(ra.shapes,a.source().body_surfaces()[port].body.shapes);
            for (x,y) in ra.shapes.iter().zip(&rb.shapes) {
                assert!((y-x*0.5_f64.sqrt()).abs()<1e-12*x.abs().max(1.0));
            }
        }
        let zero=vec![0.0;a.system().mode_count()];
        let (mut pressure_difference,mut traction_difference)=(0.0,0.0);
        for _ in 0..256 {
            let drives=[drive(0.001,0.1),drive(-0.0005,0.1)];
            let fa=a.step(&zero,&drives).unwrap().clone();
            let fb=b.step(&zero,&drives).unwrap().clone();
            pressure_difference+=(fa.observer_pressure_pa-fb.observer_pressure_pa).abs();
            traction_difference+=(fa.ports[0].traction_n-fb.ports[0].traction_n).abs();
            assert!(fa.energy_residual_j.abs()<=fa.energy_tolerance_j);
            assert!(fb.energy_residual_j.abs()<=fb.energy_tolerance_j);
        }
        assert!(pressure_difference>1e-12 && traction_difference>1e-12);
        assert!(a.system().components()[2].states().iter().any(|s|s.velocity_m_sqrt_kg_per_s!=0.0));
        let surface=a.source().body_surfaces()[0].clone();
        assert!(a.with_bow_body_surfaces(vec![surface.clone()],gate).is_err());
        let mut absent=surface;absent.body.component=2; // after the plate and cavity
        assert!(joint(cx,gate,7800.0).with_bow_body_surfaces(vec![absent],gate).is_err());
    });
}

fn ribbon_limits()->MultiFrictionConfig {
    MultiFrictionConfig {max_ports:15,max_sweeps:128,max_setup_terms:1_000_000}
}

#[test]
fn finite_ribbon_quadrature_converges_to_the_uniform_physical_load_integral() {
    // G1: integrate the actual sine modes independently over a fixed interval.
    // This tests the mechanical force map, not a spectral filter on the output.
    with_cx(|cx,gate| {
        let (width,centre)=(0.08,0.2);
        let mut previous_error=f64::INFINITY;
        for cells in [1,3,9] {
            let bank=base(cx,gate,7800.0).with_uniform_bow_ribbon(width,cells,ribbon_limits(),gate).unwrap();
            let ribbon=bank.source().ribbon().unwrap();
            let card=&bank.cards()[0];
            assert_eq!(ribbon.centre(),options().bow);
            assert_eq!(ribbon.width_m(),width);
            assert_eq!(ribbon.stations().len(),cells);
            let mut error=0.0;
            for k in 1..=card.mode_count {
                let wavenumber=k as f64*core::f64::consts::PI;
                let half_phase=wavenumber*width/(2.0*card.length_m);
                let exact=(wavenumber*centre).sin()*half_phase.sin()/half_phase
                    /(card.linear_density_kg_m*card.length_m/2.0).sqrt();
                let actual=(0..cells).map(|port|bank.system().port(port).unwrap().left.shapes[k-1]
                    *ribbon.load_fraction()).sum::<f64>();
                error+=(actual-exact).abs();
            }
            assert!(error<previous_error/4.0,"midpoint load integration did not refine: {error} >= {previous_error}");
            previous_error=error;
            let mut drives=vec![drive(0.0,0.0);cells];
            ribbon.project_drive(drive(-0.04,0.7),&mut drives).unwrap();
            assert!(drives.iter().all(|d|d.speed_m_s == -0.04));
            assert!((drives.iter().map(|d|d.normal_force_n).sum::<f64>()-0.7).abs()<2e-15);
        }
    });
}

#[test]
fn finite_ribbon_solves_separate_slip_with_shared_body_and_conserved_work() {
    with_cx(|cx,gate| {
        let mut original=base(cx,gate,15600.0);
        let zero=vec![0.0;original.system().mode_count()];
        for _ in 0..7 { original.step(&zero,drive(0.1,0.0001)).unwrap(); }
        let initial:Vec<_>=original.system().components().iter().flat_map(|c|c.states().iter().copied()).collect();
        let mut bank=original.with_uniform_bow_ribbon(0.03,3,ribbon_limits(),gate).unwrap();
        assert_eq!(states(&bank),initial);
        assert_eq!(bank.system().samples_rendered(),7);
        assert!(bank.source().bank().material_cavity().is_some());
        assert_eq!(bank.source().bank().strings().len(),2);
        let ribbon=bank.source().ribbon().unwrap().clone();
        let mut drives=vec![drive(0.0,0.0);3];
        let (mut nonuniform_slip,mut peak)=(0.0_f64,0.0_f64);
        for sample in 0..256 {
            ribbon.project_drive(drive(0.1,if sample<192 {0.0003} else {0.0}),&mut drives).unwrap();
            let frame=bank.step(&zero,&drives).unwrap();
            assert!(frame.energy_residual_j.abs()<=frame.energy_tolerance_j);
            assert!(frame.ports.iter().all(|p|p.dissipation_j>=0.0));
            nonuniform_slip=nonuniform_slip.max((frame.ports[1].slip_distance_m-frame.ports[2].slip_distance_m).abs());
            peak=peak.max(frame.observer_pressure_pa.abs());
            if sample>=192 { assert_eq!(frame.drive_work_j,0.0);assert_eq!(frame.friction_dissipation_j,0.0); }
        }
        assert!(peak>1e-12);
        assert!(nonuniform_slip>1e-14,"spatial cells were collapsed to one averaged slip");
        assert!(bank.system().total_energy_j().unwrap()>0.0,"ribbon release erased sympathetic ringing");
    });
}

#[test]
fn ribbon_gestures_keep_original_total_load_ramps_and_block_invariance() {
    with_cx(|cx,gate| {
        let source=GestureSchedule::try_new(700,vec![GestureTrack {
            id:"ribbon".into(),target:GestureTarget::BowStroke {string:0},
            initial:value(0.1,0.0003,0.2),events:vec![
                GestureEvent {time_s:0.002,transition_s:0.006,value:value(0.2,0.0006,0.2)},
                GestureEvent {time_s:0.004,transition_s:0.0,value:value(0.1,0.00015,0.2)},
                GestureEvent {time_s:0.007,transition_s:0.0,value:value(0.0,0.0,0.2)},
            ],
        }]).unwrap();
        let build=||base(cx,gate,7800.0).with_uniform_bow_ribbon(0.03,3,ribbon_limits(),gate).unwrap()
            .into_scheduled_ribbon(source.clone(),"ribbon".into(),performance_config()).unwrap();
        let mut a=build();let mut b=build();
        assert_eq!(a.schedule().to_canonical_bytes(),source.to_canonical_bytes());
        let compiled=fs_couple::render::schedule::force::coupled::render::friction::compile_friction_track(
            &source,"ribbon",FrictionGestureConfig {sample_rate_hz:48000,samples:513,max_block:128,
                max_work:10000,max_events:128,station_fraction:0.2}).unwrap();
        assert_eq!(a.events().len(),compiled.events.len()*3);
        for (chunk,event) in a.events().chunks_exact(3).zip(&compiled.events) {
            assert!(chunk.iter().all(|e|e.sample==event.sample && e.drive.speed_m_s==event.drive.speed_m_s));
            assert!((chunk.iter().map(|e|e.drive.normal_force_n).sum::<f64>()-event.drive.normal_force_n).abs()<1e-18);
        }
        let mut first=vec![0.0;513];let mut second=first.clone();
        for block in first.chunks_mut(1) { a.block(block).unwrap(); }
        for block in second.chunks_mut(127) { b.block(block).unwrap(); }
        assert_eq!(first,second);assert_eq!(states(a.bank()),states(b.bank()));
        assert!(first.iter().any(|p|p.abs()>1e-12));
    });
}

#[test]
fn finite_ribbon_refuses_invalid_geometry_budgets_and_partial_drives() {
    with_cx(|cx,gate| {
        for (width,cells) in [(0.0,3),(-0.01,3),(f64::NAN,3),(0.25,3),(0.01,0),(0.01,2),(0.01,65)] {
            assert!(base(cx,gate,7800.0).with_uniform_bow_ribbon(width,cells,ribbon_limits(),gate).is_err());
        }
        assert!(base(cx,gate,7800.0).with_uniform_bow_ribbon(0.01,5,joint_limits(),gate).is_err());
        let bank=base(cx,gate,7800.0).with_uniform_bow_ribbon(0.01,3,ribbon_limits(),gate).unwrap();
        let ribbon=bank.source().ribbon().unwrap();
        let mut partial=[drive(1.0,2.0);2];
        assert!(ribbon.project_drive(drive(0.1,0.1),&mut partial).is_err());
        assert!(partial.iter().all(|d|d.speed_m_s==1.0 && d.normal_force_n==2.0));
        let mut complete=[drive(1.0,2.0);3];
        assert!(ribbon.project_drive(drive(0.1,-0.1),&mut complete).is_err());
        assert!(complete.iter().all(|d|d.speed_m_s==1.0 && d.normal_force_n==2.0));
        assert!(bank.into_scheduled_ribbon(schedule(),"bow-0".into(),performance_config()).is_err());
    });
}

#[test]
fn adding_material_contacts_preserves_vibration_sources_and_physical_pressure_basis() {
    with_cx(|cx,gate| {
        let mut original=base(cx,gate,15600.0);
        let zero=vec![0.0;original.system().mode_count()];
        for _ in 0..19 {original.step(&zero,drive(0.01,0.1)).unwrap();}
        let old:Vec<_>=original.system().components().iter().flat_map(|m|m.states().iter().copied()).collect();
        let weights=original.system().acceleration_pressure_weights().unwrap().to_vec();
        let card=original.cards()[1].clone();
        let body_id=original.source().material_plate().unwrap().specimen.material().identity();
        let bank=original.with_additional_bows(vec![extra()],joint_limits(),gate).unwrap();
        assert_eq!(states(&bank),old); assert_eq!(bank.system().samples_rendered(),19);
        assert_eq!(bank.system().acceleration_pressure_weights(),Some(weights.as_slice()));
        assert_eq!(bank.source().bank().material_plate().unwrap().specimen.material().identity(),body_id);
        assert!(bank.source().bank().material_cavity().is_some());
        assert_eq!(bank.source().bank().strings().len(),2);
        let port=bank.system().port(1).unwrap(); assert_eq!(port.left.component,1);
        for (k,&shape) in port.left.shapes.iter().enumerate() {
            let expected=((k+1) as f64*core::f64::consts::PI*0.3).sin()/(card.linear_density_kg_m*card.length_m/2.0).sqrt();
            assert!((shape-expected).abs()<1e-12*expected.abs().max(1.0));
        }
    });
}

#[test]
fn second_material_and_contact_change_joint_reactions_and_radiated_sound() {
    with_cx(|cx,gate| {
        let mut a=joint(cx,gate,7800.0);let mut b=joint(cx,gate,15600.0);let mut released=joint(cx,gate,7800.0);
        let zero=vec![0.0;a.system().mode_count()];
        for (&x,&y) in a.system().port(1).unwrap().left.shapes.iter().zip(&b.system().port(1).unwrap().left.shapes) {
            assert!((y/x-0.5_f64.sqrt()).abs()<1e-12);
        }
        let (mut material_difference,mut contact_difference,mut peak)=(0.0_f64,0.0_f64,0.0_f64);
        for sample in 0..512 {
            let load=if sample<384 {0.1} else {0.0};
            let drives=[drive(0.01,load),drive(0.01,load)];
            let fa=a.step(&zero,&drives).unwrap().clone();let fb=b.step(&zero,&drives).unwrap().clone();
            let fr=released.step(&zero,&[drives[0],drive(0.01,0.0)]).unwrap().clone();
            for f in [&fa,&fb,&fr] {assert!(f.energy_residual_j.abs()<=f.energy_tolerance_j);}
            material_difference+=(fa.observer_pressure_pa-fb.observer_pressure_pa).abs();
            contact_difference+=(fa.ports[0].traction_n-fr.ports[0].traction_n).abs();
            peak=peak.max(fa.observer_pressure_pa.abs());
            if sample>=384 {assert_eq!(fa.drive_work_j,0.0);assert_eq!(fa.friction_dissipation_j,0.0);}
        }
        assert!(peak>1e-12);assert!(material_difference>1e-12);assert!(contact_difference>1e-12);
        assert!(a.system().total_energy_j().unwrap()>0.0,"release erased body/cavity vibration");
    });
}

#[test]
fn all_tracks_use_source_ramps_before_one_joint_step_independent_of_callbacks() {
    with_cx(|cx,gate| {
        for reacted in [false,true] {
        let make=|| {
            let bank=joint(cx,gate,7800.0);
            if reacted {reacting_surface(bank,gate)} else {bank}
        };
        let source=schedule();let mut direct=make();
        let zero=vec![0.0;direct.system().mode_count()];let mut expected=Vec::new();
        for sample in 0..513_u64 {
            let tick=sample*700/48000;
            let drives:Vec<_>=tracks().iter().map(|id| {
                let GestureValue::Bow {velocity_m_per_s,normal_force_n,..}=source.sample_value(id,tick).unwrap() else {panic!("bow")};
                drive(velocity_m_per_s,normal_force_n)
            }).collect();
            expected.push(direct.step(&zero,&drives).unwrap().observer_pressure_pa.to_bits());
        }
        for partition in [1,37,128] {
            let mut performance=make().into_scheduled(source.clone(),tracks(),performance_config()).unwrap();
            assert!(performance.events().windows(2).any(|e|e[0].sample==e[1].sample && e[0].port!=e[1].port));
            let mut output=vec![0.0;513];for block in output.chunks_mut(partition) {performance.block(block).unwrap();}
            assert_eq!(output.iter().map(|p|p.to_bits()).collect::<Vec<_>>(),expected);
            assert_eq!(states(performance.bank()),states(&direct));
            assert_eq!(performance.bank().source().body_surfaces(),direct.source().body_surfaces());
            assert_eq!(performance.samples_rendered(),513);assert_eq!(performance.remaining_samples(),0);
            let mut sentinel=[123.0];assert!(performance.block(&mut sentinel).is_err());assert_eq!(sentinel,[123.0]);
        }
        }
    });
}

#[test]
fn singleton_schedule_keeps_original_renderer_and_prior_physical_clock_exact() {
    with_cx(|cx,gate| {
        let mut a=base(cx,gate,7800.0);let mut b=base(cx,gate,7800.0);
        let zero=vec![0.0;a.system().mode_count()];
        for _ in 0..17 {a.step(&zero,drive(0.01,0.1)).unwrap();b.step(&zero,drive(0.01,0.1)).unwrap();}
        let source=schedule();
        let mut legacy=a.into_scheduled(&source,"bow-0",FrictionGestureConfig {sample_rate_hz:48000,samples:513,
            max_block:128,max_work:10000,max_events:128,station_fraction:0.2}).unwrap();
        let mut joint=b.with_additional_bows(vec![],joint_limits(),gate).unwrap()
            .into_scheduled(source,vec!["bow-0".into()],performance_config()).unwrap();
        let mut expected=vec![0.0;513];let mut actual=expected.clone();
        for block in expected.chunks_mut(128) {legacy.renderer_mut().block(block).unwrap();}
        for block in actual.chunks_mut(37) {joint.block(block).unwrap();}
        assert_eq!(expected.iter().map(|p|p.to_bits()).collect::<Vec<_>>(),actual.iter().map(|p|p.to_bits()).collect::<Vec<_>>());
        assert_eq!(joint.samples_rendered(),513);assert_eq!(joint.bank().system().samples_rendered(),530);
    });
}

#[test]
fn wrong_contacts_tracks_and_aggregate_control_budgets_refuse() {
    with_cx(|cx,gate| {
        for location in [StringBowLocation {string:2,station_fraction:0.3},StringBowLocation {string:1,station_fraction:1.0},StringBowLocation {string:1,station_fraction:f64::NAN}] {
            let mut bow=extra();bow.location=location;
            assert!(base(cx,gate,7800.0).with_additional_bows(vec![bow],joint_limits(),gate).is_err());
        }
        for bindings in [vec!["bow-0".into()],vec!["bow-1".into(),"bow-0".into()],vec!["bow-0".into(),"bow-0".into()]] {
            assert!(joint(cx,gate,7800.0).into_scheduled(schedule(),bindings,performance_config()).is_err());
        }
        for config in [MultiBowPerformanceConfig {max_work:60,..performance_config()},
            MultiBowPerformanceConfig {max_events:5,..performance_config()},MultiBowPerformanceConfig {samples:0,..performance_config()}] {
            assert!(joint(cx,gate,7800.0).into_scheduled(schedule(),tracks(),config).is_err());
        }
        let mut performance=joint(cx,gate,7800.0).into_scheduled(schedule(),tracks(),performance_config()).unwrap();
        let before=states(performance.bank());let mut oversized=[123.0;129];
        assert!(performance.block(&mut oversized).is_err());assert_eq!(oversized,[123.0;129]);
        assert_eq!(states(performance.bank()),before);assert_eq!(performance.samples_rendered(),0);
        assert!(performance.validate_sample_rate(44100).is_err());performance.block(&mut [0.0]).unwrap();
        // Physical refusal is different from sizing: no continued audio after a failed callback.
        let mut limited=extra();limited.friction.maximum_force_n=1e-30;
        let mut failing=base(cx,gate,7800.0).with_additional_bows(vec![limited],joint_limits(),gate).unwrap()
            .into_scheduled(schedule(),tracks(),performance_config()).unwrap();
        let before=states(failing.bank());let mut sample=[123.0];
        assert!(failing.block(&mut sample).is_err());assert_eq!(sample,[123.0]);
        assert_eq!(states(failing.bank()),before);assert_eq!(failing.samples_rendered(),0);
        assert!(matches!(failing.block(&mut sample),Err(RenderError::Poisoned)));
    });
}
