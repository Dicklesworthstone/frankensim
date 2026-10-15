//! Synthetic sources test physical wiring; no measured instrument timbre claim.
use fs_couple::bowed_string::runtime::reciprocal::{InteriorStringConnection, ReciprocalStringError};
use fs_couple::bowed_string::runtime::reciprocal::bank::{
    MaterialStringMember, MaterialStringSource, ReciprocalStringBank, StringBankOptions, StringBowLocation,
};
use fs_couple::bowed_string::runtime::reciprocal::bank::plate::{SharedMaterialPlateBody, SharedPlateBridge};
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::{MaterialBowContact, MultiBowStringBank};
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::performance::MultiBowPerformanceConfig;
use fs_couple::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeModel, ModalAcousticTimeBudget};
use fs_couple::pcm_wav::observation::PressureRenderer;
use fs_couple::render::schedule::force::coupled::{ModalAttachment, ModalCouplingConfig};
use fs_couple::render::schedule::force::coupled::contact::friction::{FrictionDrive, ModalFrictionConfig, ModalFrictionRegime};
use fs_couple::render::schedule::force::coupled::contact::friction::multiple::MultiFrictionConfig;
use fs_couple::string_specimen::{KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY as ETA,
    StringPrestress, with_uniform_circular_material_state};
use fs_couple::thin_plate::{PlateMaterialModel, PlateThicknessConstraint, with_uniform_plate_material_state};
use fs_evidence::ValidityDomain;
use fs_exec::{CancelGate, Cx};
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::gas::{GasSpec, GasState};
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::state_point::{MaterialPropertySelection, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point};
use fs_qty::{Density, Dims, DynViscosity, Pressure};
use fs_scenario::{PrestressedString, ThinPlate, gesture::{GestureSchedule, GestureTarget, GestureTrack, GestureValue}};
use fs_tribo::FrictionLaw;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate=CancelGate::new_clock_free();
    let pool=fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| f(&Cx::new(&gate,arena,fs_exec::StreamKey { seed:6,kernel_id:5,tile:0,iteration:0 },
        fs_exec::Budget::INFINITE,fs_exec::ExecMode::Deterministic),&gate))
}
fn material(rho:f64,g:f64,eta:f64,hi:f64)->ResolvedMaterialStatePoint {
    let mut claims=ClaimSet::new(); let mut requirements=Vec::new();
    for (name,dims,value) in [("density",Density::DIMS,rho),("young_modulus",Pressure::DIMS,200e9),
        ("poisson_ratio",Dims::NONE,0.3),(ETA,DynViscosity::DIMS,10.0),
        ("G",Pressure::DIMS,g),("eta_s",DynViscosity::DIMS,eta)] {
        claims.insert_claim(PropertyClaim {key:PropertyKey::new(name,dims),value:PropertyValue::Scalar{value,dims},
            validity:ValidityDomain::unconstrained().with("T",280.0,320.0).with("omega",0.0,hi),
            uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::ConstantWithinValidity,
            observations:vec![],provenance:Provenance{source:"synthetic material twist".into(),license:"CC0-1.0".into(),artifact:None},
        }).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(name,dims,ScalarAdmissibility::Finite).unwrap());
    }
    let card=MaterialCard::assemble(MaterialStateId{chemistry:"synthetic".into(),phase:"solid".into(),process:"test".into(),revision:0},claims,vec![]).unwrap();
    resolve_material_state_point(&card,&QueryPoint::new().with("T",293.15).unwrap().with("omega",1000.0).unwrap(),
        &requirements,MaterialPropertySelection::SingleClaimOnly).unwrap()
}
fn template()->PrestressedString {
    PrestressedString {length_m:0.5,tension_n:20.0,lin_density_kg_m:1.0,axial_stiffness_n:1.0,width_m:1.0,
        n_modes:4,damping_ratio:0.0,rayleigh:None,bending_stiffness_n_m2:0.0,kelvin_voigt_bending:None,
        relaxation_bending:None,polarization_detune:0.0,moving_end:false}
}
fn source(g:f64,eta:f64)->MaterialStringSource {
    MaterialStringSource::Homogeneous(with_uniform_circular_material_state(template(),0.0003,
        &material(7800.0,g,eta,200_000.0)).unwrap().with_kelvin_voigt_bending_loss().unwrap())
}
fn gas()->GasState {GasState::try_new(&GasSpec::dry_air_ussa1976(),293.15,101_325.0).unwrap()}
fn law()->FrictionLaw {FrictionLaw::Coulomb{static_mu:0.8,kinetic_mu:0.4}}
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
fn raw(cx:&Cx<'_>,gate:&CancelGate,src:MaterialStringSource,o:StringBankOptions)->Result<ReciprocalStringBank,ReciprocalStringError> {
    let body=ModalAcousticTimeModel::try_new(o.sample_rate_hz,vec![ModalAcousticMode{
        angular_frequency_rad_s:450.0,damping_ratio:0.01,pressure_per_modal_velocity:fs_math::c64::C64::ZERO,
    }],ModalAcousticTimeBudget::audible_reference()).unwrap();
    ReciprocalStringBank::new(cx,vec![MaterialStringMember{source:src,connections:vec![InteriorStringConnection{
        station_fraction:0.8,body:ModalAttachment{component:0,shapes:vec![2.0]},stiffness_n_m:3000.0,
        damping_n_s_m:0.1,rest_extension_m:0.0}]}],gas(),vec![body],law(),o,gate)
}
fn body(count:usize)->SharedMaterialPlateBody {
    let m=material(7800.0,80e9,1000.0,200_000.0);
    let specimen=with_uniform_plate_material_state(ThinPlate{length_m:0.4,width_m:0.3,thickness_m:1.0,
        density_kg_m3:1.0,e1_pa:1.0,e2_pa:1.0,nu12:0.0,g12_pa:1.0,material_angle_rad:0.0,
        damping_ratio:0.0,thermoelastic:None,kelvin_voigt_bending:None,n_modes:2,geometric_nonlinearity:false,
        pretension_n_m:0.0,clamped:false},&m,PlateMaterialModel::Isotropic,PlateThicknessConstraint::FixedThickness(0.001))
        .unwrap().with_kelvin_voigt_bending_loss().unwrap();
    SharedMaterialPlateBody{specimen,bridges:vec![SharedPlateBridge{station_fraction:0.8,stiffness_n_m:3000.0,
        damping_n_s_m:0.1,rest_extension_m:0.0};count],modal_budget:ModalAcousticTimeBudget::audible_reference(),listener_m:1.0}
}
fn radiating(cx:&Cx<'_>,gate:&CancelGate,g:f64)->ReciprocalStringBank {
    ReciprocalStringBank::from_material_plate(cx,vec![source(g,1000.0).with_torsion("G","eta_s",3).unwrap()],
        gas(),body(1),law(),options(),gate).unwrap()
}
fn drive()->FrictionDrive {FrictionDrive{speed_m_s:0.01,normal_force_n:0.1}}

#[test]
fn sourced_twist_has_physical_frequency_loss_and_surface_torque_normalization() {
    with_cx(|cx,gate| {
        let bank=raw(cx,gate,source(80e9,1000.0).with_torsion("G","eta_s",3).unwrap(),options()).unwrap();
        let MaterialStringSource::WithTorsion(t)=&bank.source().strings()[0].source else {panic!("lost torsional source")};
        let rotary=t.section().rotary_inertia_kg_m(); let radius=t.section().shear().elastic().outer_radius_m();
        assert_eq!(t.mode_count(),3); assert_eq!(bank.system().components()[0].modes().len(),7);
        let modes=bank.system().components()[0].modes(); let phi=&bank.system().port().left.shapes;
        for n in 1..=3 {
            let k=n as f64*core::f64::consts::PI/0.5;
            let omega=k*(80e9/7800.0_f64).sqrt();
            assert!((modes[3+n].angular_frequency_rad_s/omega-1.0).abs()<2e-14);
            let expected_zeta=(1000.0/7800.0)*k*k/(2.0*omega);
            assert!((modes[3+n].damping_ratio/expected_zeta-1.0).abs()<2e-14);
            let angle=(n as f64*core::f64::consts::PI*0.2).sin()/(rotary*0.5*0.5).sqrt();
            assert!((phi[3+n]/(radius*angle)-1.0).abs()<2e-14);
        }
        assert!(modes.iter().all(|m|m.pressure_per_modal_velocity==fs_math::c64::C64::ZERO));
    });
}

#[test]
fn sticking_constrains_surface_motion_not_transverse_centre_motion_and_release_keeps_twist() {
    with_cx(|cx,gate| {
        let mut bank=raw(cx,gate,source(80e9,1000.0).with_torsion("G","eta_s",3).unwrap(),options()).unwrap();
        let zero=vec![0.0;bank.system().mode_count()]; let phi=bank.system().port().left.shapes.clone();
        let mut twist_distance=0.0_f64;
        for _ in 0..32 {
            let before=bank.system().components()[0].states().to_vec();
            let f=*bank.step(&zero,FrictionDrive{speed_m_s:0.01,normal_force_n:1.0}).unwrap();
            assert_eq!(f.regime,ModalFrictionRegime::Sticking);
            let after=bank.system().components()[0].states();
            let delta:Vec<_>=after.iter().zip(&before).map(|(a,b)|a.displacement_m_sqrt_kg-b.displacement_m_sqrt_kg).collect();
            let centre:f64=phi[..4].iter().zip(&delta).map(|(p,d)|p*d).sum();
            let twist:f64=phi[4..].iter().zip(&delta[4..]).map(|(p,d)|p*d).sum();
            assert!((centre+twist-0.01/48_000.0).abs()<1e-12);
            assert!(f.energy_residual_j.abs()<=f.energy_tolerance_j);
            twist_distance=twist_distance.max(twist.abs());
        }
        assert!(twist_distance>1e-10,"surface rotation was discarded");
        let initial=bank.system().total_energy_j().unwrap();
        for _ in 0..64 {
            let f=bank.step(&zero,FrictionDrive{speed_m_s:0.0,normal_force_n:0.0}).unwrap();
            assert_eq!(f.regime,ModalFrictionRegime::Released); assert_eq!(f.drive_work_j,0.0);
        }
        assert!(bank.system().total_energy_j().unwrap()<=initial+1e-12);
        assert!(bank.system().components()[0].states()[4..].iter().any(|s|s.velocity_m_sqrt_kg_per_s!=0.0));
    });
}

#[test]
fn shear_only_material_swap_changes_string_motion_and_real_body_pressure_without_retuning() {
    with_cx(|cx,gate| {
        let mut a=radiating(cx,gate,80e9); let mut b=radiating(cx,gate,20e9);
        assert_eq!(a.cards()[0].fundamental_hz(),b.cards()[0].fundamental_hz());
        assert_eq!(a.cards()[0].zetas,b.cards()[0].zetas);
        assert_eq!(a.system().components()[1].modes(),b.system().components()[1].modes());
        let zero=vec![0.0;a.system().mode_count()];let mut dp=0.0_f64;let mut dq=0.0_f64;
        for _ in 0..512 {
            let pa=a.step(&zero,drive()).unwrap().observer_pressure_pa;
            let pb=b.step(&zero,drive()).unwrap().observer_pressure_pa;
            dp=dp.max((pa-pb).abs());
            dq=dq.max((a.system().components()[0].states()[0].displacement_m_sqrt_kg
                -b.system().components()[0].states()[0].displacement_m_sqrt_kg).abs());
        }
        assert!(dp>1e-12 && dq>1e-12,"twist only changed a diagnostic: dp={dp:e}, dq={dq:e}");
        let weights=a.system().acceleration_pressure_weights().unwrap();
        assert!(weights[..7].iter().all(|w|*w==0.0));
        assert!(weights[7..].iter().any(|w|*w!=0.0));
    });
}

#[test]
fn layered_twist_uses_all_original_annuli_and_the_same_selected_shear_data() {
    with_cx(|cx,gate| {
        let section=ConcentricElasticSection::try_new(0.0,vec![
            ConcentricMaterialLayer{outer_radius_m:0.0002,material:material(7800.0,80e9,1000.0,200_000.0)},
            ConcentricMaterialLayer{outer_radius_m:0.0003,material:material(15_600.0,20e9,2000.0,200_000.0)},
        ]).unwrap().with_kelvin_voigt_bending(ETA).unwrap();
        let original=section.elastic().clone();
        let s=MaterialStringSource::Concentric{template:template(),section,prestress:StringPrestress::FixedTension(20.0)}
            .with_torsion("G","eta_s",2).unwrap();
        let b=raw(cx,gate,s,options()).unwrap();
        let MaterialStringSource::WithTorsion(t)=&b.source().strings()[0].source else {panic!("missing source")};
        assert_eq!(t.section().shear().elastic(),&original);
        let MaterialStringSource::Concentric{section,..}=t.bending() else {panic!("lost layered bending")};
        assert_eq!(section.elastic(),&original);
        assert_eq!(b.cards()[0].linear_density_kg_m,original.line_density_kg_m());
        assert_eq!(b.system().components()[0].modes().len(),6);
    });
}

fn multi(cx:&Cx<'_>,gate:&CancelGate)->MultiBowStringBank {
    let bank=ReciprocalStringBank::from_material_plate(cx,vec![
        source(80e9,1000.0).with_torsion("G","eta_s",3).unwrap(),
        source(20e9,1000.0).with_torsion("G","eta_s",2).unwrap(),
    ],gas(),body(2),law(),options(),gate).unwrap();
    bank.with_additional_bows(vec![MaterialBowContact{location:StringBowLocation{string:1,station_fraction:0.31},
        law:law(),friction:options().friction}],MultiFrictionConfig{max_ports:4,max_sweeps:128,max_setup_terms:16_384},gate).unwrap()
}
#[test]
fn additional_bows_and_existing_scheduler_keep_complete_torsional_basis_and_partition_equivalence() {
    with_cx(|cx,gate| {
        let schedule=GestureSchedule::try_new(700,vec![
            GestureTrack{id:"a".into(),target:GestureTarget::BowStroke{string:0},initial:GestureValue::Bow{
                velocity_m_per_s:0.01,normal_force_n:0.1,station:0.2},events:vec![]},
            GestureTrack{id:"b".into(),target:GestureTarget::BowStroke{string:1},initial:GestureValue::Bow{
                velocity_m_per_s:-0.01,normal_force_n:0.08,station:0.31},events:vec![]},
        ]).unwrap();
        let mut direct=multi(cx,gate);
        assert_eq!(direct.system().port(0).unwrap().left.shapes.len(),7);
        assert_eq!(direct.system().port(1).unwrap().left.shapes.len(),6);
        assert!(direct.system().port(1).unwrap().left.shapes[4..].iter().any(|x|*x!=0.0));
        let zero=vec![0.0;direct.system().mode_count()];
        let expected:Vec<_>=(0..129).map(|_|direct.step(&zero,&[drive(),FrictionDrive{speed_m_s:-0.01,normal_force_n:0.08}])
            .unwrap().observer_pressure_pa.to_bits()).collect();
        for partition in [1,37,128] {
            let mut p=multi(cx,gate).into_scheduled(schedule.clone(),vec!["a".into(),"b".into()],
                MultiBowPerformanceConfig{samples:129,max_block:128,max_work:1000,max_events:16}).unwrap();
            let mut out=vec![0.0;129];
            for block in out.chunks_mut(partition) {p.block(block).unwrap();}
            assert_eq!(out.iter().map(|p|p.to_bits()).collect::<Vec<_>>(),expected);
            assert_eq!(p.samples_rendered(),129);
        }
    });
}

#[test]
fn source_band_nyquist_complete_dimension_and_invalid_selectors_refuse_without_fallback() {
    with_cx(|cx,gate| {
        let base=source(80e9,1000.0);
        for n in [0,4097] {assert!(base.clone().with_torsion("G","eta_s",n).is_err());}
        assert!(base.clone().with_torsion("G","missing",1).is_err());
        assert!(base.clone().with_torsion("G","G",1).is_err());
        assert!(base.clone().with_torsion("G","eta_s",1).unwrap().with_torsion("G","eta_s",1).is_err());
        let band=MaterialStringSource::Homogeneous(with_uniform_circular_material_state(template(),0.0003,
            &material(7800.0,80e9,1000.0,10_000.0)).unwrap().with_kelvin_voigt_bending_loss().unwrap());
        assert!(raw(cx,gate,band.with_torsion("G","eta_s",1).unwrap(),options()).is_err());
        let mut too_small=options();too_small.network.max_modes=5;
        assert!(raw(cx,gate,base.clone().with_torsion("G","eta_s",1).unwrap(),too_small).is_err());
        let mut low_rate=options();low_rate.sample_rate_hz=8000;
        assert!(raw(cx,gate,base.with_torsion("G","eta_s",3).unwrap(),low_rate).is_err());
    });
}

#[test]
fn refused_steps_keep_both_transverse_and_twist_history() {
    with_cx(|cx,gate| {
        let make=||raw(cx,gate,source(80e9,1000.0).with_torsion("G","eta_s",3).unwrap(),options()).unwrap();
        let mut a=make();let mut b=make();let zero=vec![0.0;a.system().mode_count()];
        for _ in 0..16 {assert_eq!(a.step(&zero,drive()).unwrap(),b.step(&zero,drive()).unwrap());}
        let before=a.system().components()[0].states().to_vec();
        assert!(a.step(&zero,FrictionDrive{speed_m_s:f64::NAN,normal_force_n:0.1}).is_err());
        assert_eq!(a.system().samples_rendered(),16);
        assert_eq!(a.system().components()[0].states(),before);
        assert_eq!(a.step(&zero,drive()).unwrap(),b.step(&zero,drive()).unwrap());
    });
}
