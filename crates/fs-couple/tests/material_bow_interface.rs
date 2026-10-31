//! Synthetic source-bound friction on actual material string/body/twist modes.
use fs_couple::bowed_string::runtime::reciprocal::bank::{MaterialStringSource,ReciprocalStringBank,StringBankOptions,StringBowLocation};
use fs_couple::bowed_string::runtime::reciprocal::bank::plate::{SharedMaterialPlateBody,SharedPlateBridge};
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::{MaterialBowContact,MultiBowStringBank};
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::performance::MultiBowPerformanceConfig;
use fs_couple::modal_acoustic_time::ModalAcousticTimeBudget;
use fs_couple::pcm_wav::observation::PressureRenderer;
use fs_couple::render::schedule::force::coupled::ModalCouplingConfig;
use fs_couple::render::schedule::force::coupled::contact::friction::{FrictionDrive,ModalFrictionConfig};
use fs_couple::render::schedule::force::coupled::contact::friction::interface::DryFrictionBinding;
use fs_couple::render::schedule::force::coupled::contact::friction::multiple::MultiFrictionConfig;
use fs_couple::string_specimen::{KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY as ETA,with_uniform_circular_material_state};
use fs_couple::thin_plate::{PlateMaterialModel,PlateThicknessConstraint,with_uniform_plate_material_state};
use fs_evidence::ValidityDomain;
use fs_exec::{CancelGate,Cx};
use fs_matdb::{ClaimSet,InterpolationPolicy,MaterialCard,MaterialStateId,PropertyClaim,PropertyKey,PropertyValue,Provenance,QueryPoint,UncertaintyModel};
use fs_material::gas::{GasSpec,GasState};
use fs_material::state_point::{MaterialPropertySelection,ResolvedMaterialStatePoint,ScalarAdmissibility,ScalarPropertyRequirement,resolve_material_state_point};
use fs_qty::{Density,Dims,DynViscosity,Pressure};
use fs_scenario::{PrestressedString,ThinPlate,gesture::{GestureSchedule,GestureTarget,GestureTrack,GestureValue}};
use fs_tribo::{ApplicabilityRange as Range,DryFrictionApplicability,DryInterfaceSystemCard,FrictionLaw,InputAuthority,InterfaceMedium,InterfaceSystemRef};

fn with_cx<R>(f:impl FnOnce(&Cx<'_>,&CancelGate)->R)->R {
    let gate=CancelGate::new_clock_free();let pool=fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena|f(&Cx::new(&gate,arena,fs_exec::StreamKey {seed:6,kernel_id:5,tile:0,iteration:0},
        fs_exec::Budget::INFINITE,fs_exec::ExecMode::Deterministic),&gate))
}
fn material()->ResolvedMaterialStatePoint { material_at_density(7800.0) }
fn material_at_density(rho:f64)->ResolvedMaterialStatePoint {
    let mut claims=ClaimSet::new();let mut requirements=Vec::new();
    for (name,dims,value) in [("density",Density::DIMS,rho),("young_modulus",Pressure::DIMS,200e9),
        ("poisson_ratio",Dims::NONE,0.3),(ETA,DynViscosity::DIMS,10.0),
        ("G",Pressure::DIMS,80e9),("eta_s",DynViscosity::DIMS,1000.0)] {
        claims.insert_claim(PropertyClaim {key:PropertyKey::new(name,dims),value:PropertyValue::Scalar{value,dims},
            validity:ValidityDomain::unconstrained().with("T",280.0,320.0).with("omega",0.0,200_000.0),
            uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::ConstantWithinValidity,
            observations:vec![],provenance:Provenance{source:"synthetic interface fixture".into(),license:"CC0-1.0".into(),artifact:None},
        }).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(name,dims,ScalarAdmissibility::Finite).unwrap());
    }
    let card=MaterialCard::assemble(MaterialStateId{chemistry:"synthetic".into(),phase:"solid".into(),process:"test".into(),revision:0},claims,vec![]).unwrap();
    resolve_material_state_point(&card,&QueryPoint::new().with("T",293.15).unwrap().with("omega",1000.0).unwrap(),
        &requirements,MaterialPropertySelection::SingleClaimOnly).unwrap()
}
fn source()->MaterialStringSource { source_at_density(7800.0) }
fn source_at_density(rho:f64)->MaterialStringSource {
    MaterialStringSource::Homogeneous(with_uniform_circular_material_state(PrestressedString {
        length_m:0.5,tension_n:20.0,lin_density_kg_m:1.0,axial_stiffness_n:1.0,width_m:1.0,
        n_modes:4,damping_ratio:0.0,rayleigh:None,bending_stiffness_n_m2:0.0,kelvin_voigt_bending:None,
        relaxation_bending:None,polarization_detune:0.0,moving_end:false,
        tension_loss_factor:None,
        bending_loss_factor:None,
    },0.0003,&material_at_density(rho)).unwrap().with_kelvin_voigt_bending_loss().unwrap()).with_torsion("G","eta_s",3).unwrap()
}
fn options()->StringBankOptions {
    StringBankOptions{sample_rate_hz:48_000,bow:StringBowLocation{string:0,station_fraction:0.2},
        string_budget:ModalAcousticTimeBudget::audible_reference(),
        network:ModalCouplingConfig{max_modes:32,max_connections:8,max_setup_terms:16_384,
            nyquist_guard_fraction:0.9,maximum_total_energy_j:100.0,maximum_abs_pressure_pa:1e6,
            maximum_abs_connection_force_n:100.0,solve_relative_tolerance:1e-10,
            energy_absolute_tolerance_j:1e-12,energy_relative_tolerance:1e-8},
        friction:ModalFrictionConfig{max_iterations:96,maximum_force_n:100.0,maximum_slip_speed_m_s:10.0,
            force_absolute_tolerance_n:1e-10,force_relative_tolerance:1e-9,sticking_tolerance_m:1e-12}}
}
fn body()->SharedMaterialPlateBody {
    let specimen=with_uniform_plate_material_state(ThinPlate{length_m:0.4,width_m:0.3,thickness_m:1.0,
        density_kg_m3:1.0,e1_pa:1.0,e2_pa:1.0,nu12:0.0,g12_pa:1.0,material_angle_rad:0.0,
        damping_ratio:0.0,thermoelastic:None,kelvin_voigt_bending:None,n_modes:2,geometric_nonlinearity:false,
        pretension_n_m:0.0,clamped:false},&material(),PlateMaterialModel::Isotropic,PlateThicknessConstraint::FixedThickness(0.001))
        .unwrap().with_kelvin_voigt_bending_loss().unwrap();
    SharedMaterialPlateBody{specimen,bridges:vec![SharedPlateBridge{station_fraction:0.8,stiffness_n_m:3000.0,
        damping_n_s_m:0.1,rest_extension_m:0.0}],modal_budget:ModalAcousticTimeBudget::audible_reference(),listener_m:1.0}
}
fn interface(mu:f64,maximum_speed:f64)->DryFrictionBinding {
    DryFrictionBinding::new(DryInterfaceSystemCard::new(InterfaceSystemRef::new(
        "synthetic/string->rosined-bow","synthetic/history","synthetic/interface",InputAuthority::SyntheticFixture,InterfaceMedium::Dry).unwrap(),
        FrictionLaw::Coulomb{static_mu:0.8,kinetic_mu:mu},DryFrictionApplicability::new(
            Range::new(280.0,320.0).unwrap(),Range::new(0.1,1e6).unwrap(),Range::new(0.0,maximum_speed).unwrap()).unwrap()).unwrap(),293.15,0.0001).unwrap()
}
fn build(cx:&Cx<'_>,gate:&CancelGate,binding:DryFrictionBinding)->ReciprocalStringBank {
    let gas=GasState::try_new(&GasSpec::dry_air_ussa1976(),293.15,101_325.0).unwrap();
    ReciprocalStringBank::from_material_plate(cx,vec![source()],gas,body(),binding.card().friction_law().clone(),options(),gate)
        .unwrap().with_bow_interface(binding).unwrap()
}
fn drivers()->[FrictionDrive;2] {[FrictionDrive{speed_m_s:0.3,normal_force_n:0.0001},FrictionDrive{speed_m_s:-0.2,normal_force_n:0.0001}]}
fn multi(cx:&Cx<'_>,gate:&CancelGate)->MultiBowStringBank {
    let binding=interface(0.4,10.0);
    build(cx,gate,binding.clone()).with_additional_bows(vec![MaterialBowContact{
        location:StringBowLocation{string:0,station_fraction:0.31},law:binding.card().friction_law().clone(),friction:options().friction,
    }],MultiFrictionConfig{max_ports:4,max_sweeps:128,max_setup_terms:16_384},gate).unwrap().with_bow_interface(1,binding).unwrap()
}

#[test]
fn interface_data_alone_changes_actual_traction_twist_body_motion_and_pressure() {
    with_cx(|cx,gate| {
        let mut a=build(cx,gate,interface(0.4,10.0));let mut b=build(cx,gate,interface(0.2,10.0));
        for (a,b) in a.system().components().iter().zip(b.system().components()) {assert_eq!(a.modes(),b.modes());}
        let zero=vec![0.0;a.system().mode_count()];let mut difference=0.0;
        for _ in 0..128 {
            let fa=*a.step(&zero,drivers()[0]).unwrap();let fb=*b.step(&zero,drivers()[0]).unwrap();
            difference+=(fa.observer_pressure_pa-fb.observer_pressure_pa).abs();
            assert!(fa.energy_residual_j.abs()<=fa.energy_tolerance_j);
            assert!(fb.energy_residual_j.abs()<=fb.energy_tolerance_j);
        }
        assert!(difference>1e-15,"interface change never reached body pressure");
        assert_ne!(a.system().components()[0].states(),b.system().components()[0].states());
        assert_eq!(a.bow_interface(),Some(&interface(0.4,10.0)));
        assert_eq!(a.source().material_plate().unwrap().specimen.material().identity(),b.source().material_plate().unwrap().specimen.material().identity());
    });
}

#[test]
fn source_domain_uses_torsional_surface_motion_not_just_prescribed_bow_speed() {
    with_cx(|cx,gate| {
        let mut bank=build(cx,gate,interface(0.4,0.001));let mut force=vec![0.0;bank.system().mode_count()];
        force[4]=10.0; // First twist coordinate; the prescribed bow itself is still.
        let before:Vec<_>=bank.system().components().iter().map(|m|m.states().to_vec()).collect();
        assert!(bank.step(&force,FrictionDrive{speed_m_s:0.0,normal_force_n:0.0001}).is_err());
        assert_eq!(bank.system().samples_rendered(),0);
        for (m,s) in bank.system().components().iter().zip(before) {assert_eq!(m.states(),s);}
        // The identical mechanical input is allowed when there is no contact.
        bank.step(&force,FrictionDrive{speed_m_s:0.0,normal_force_n:0.0}).unwrap();
        assert!(bank.system().components()[0].states()[4].velocity_m_sqrt_kg_per_s.abs()>0.0);
    });
}

#[test]
fn both_interfaces_survive_joint_scheduling_and_callback_partitions() {
    with_cx(|cx,gate| {
        let tracks:Vec<_>= ["first","second"].into_iter().zip(drivers()).zip([0.2,0.31]).map(|((id,d),station)|GestureTrack{
            id:id.into(),target:GestureTarget::BowStroke{string:0},initial:GestureValue::Bow{
                velocity_m_per_s:d.speed_m_s,normal_force_n:d.normal_force_n,station},events:vec![],
        }).collect();
        let schedule=GestureSchedule::try_new(1000,tracks).unwrap();
        let mut direct=multi(cx,gate);let zero=vec![0.0;direct.system().mode_count()];
        let expected:Vec<_>=(0..96).map(|_|direct.step(&zero,&drivers()).unwrap().observer_pressure_pa.to_bits()).collect();
        for partition in [1,37,64] {
            let mut performance=multi(cx,gate).into_scheduled(schedule.clone(),vec!["first".into(),"second".into()],
                MultiBowPerformanceConfig{samples:96,max_block:64,max_work:10_000,max_events:64}).unwrap();
            for port in 0..2 {assert_eq!(performance.bank().bow_interface(port),Some(&interface(0.4,10.0)));}
            let mut out=vec![0.0;96];for block in out.chunks_mut(partition){performance.block(block).unwrap();}
            assert_eq!(out.iter().map(|v|v.to_bits()).collect::<Vec<_>>(),expected);
        }
    });
}

#[path = "material_bow_interface/forces.rs"]
mod forces;

#[path = "material_bow_interface/stops.rs"]
mod stops;

#[path = "material_bow_interface/directions.rs"]
mod directions;

#[path = "material_bow_interface/contact_load.rs"]
mod contact_load;
