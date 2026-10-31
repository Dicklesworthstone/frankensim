//! Synthetic sources test physical wiring, and the measured Dominant string packs test their data path; no measured
//! instrument timbre claim.
use fs_couple::bowed_string::runtime::reciprocal::{InteriorStringConnection, ReciprocalStringError};
use fs_couple::bowed_string::runtime::reciprocal::bank::{
    MaterialStringMember, MaterialStringSource, ReciprocalStringBank, StringBankOptions, StringBowLocation,
};
use fs_couple::bowed_string::runtime::reciprocal::bank::plate::{SharedMaterialPlateBody, SharedPlateBridge};
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::{MaterialBowContact, MultiBowStringBank};
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::performance::MultiBowPerformanceConfig;
use fs_couple::bowed_string::runtime::reciprocal::bank::torsion::{
    ModalRotationAttachment, StringBodyRotationConnection, StringBridgeLocation,
    SHEAR_VISCOSITY_REFERENCE_OMEGA as TORSION_REFERENCE, SHEAR_VISCOSITY_REFERENCE_OMEGA_DIMS as TORSION_REFERENCE_DIMS,
    STRING_LOSS_SPECTRUM_BAND_RAD_S,
};
use fs_couple::shear_shell::constant_loss::ConstantLossSpectrum;
use fs_couple::bowed_string::runtime::reciprocal::bank::chart::{
    SharedMaterialPlateChartBody, PlateBridgeFootprints,
};
use fs_couple::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeModel, ModalAcousticTimeBudget};
use fs_couple::pcm_wav::observation::PressureRenderer;
use fs_couple::render::schedule::force::coupled::{ModalAttachment, ModalCouplingConfig};
use fs_couple::render::schedule::force::coupled::contact::friction::{FrictionDrive, ModalFrictionConfig, ModalFrictionRegime};
use fs_couple::render::schedule::force::coupled::contact::friction::multiple::MultiFrictionConfig;
use fs_couple::string_specimen::{KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY as ETA,
    StringPrestress, with_uniform_circular_material_state};
use fs_couple::thin_plate::{PlateMaterialModel, PlateThicknessConstraint, with_uniform_plate_material_state,
    PlateRegionMaterial, PlateChartRadiation, compile_plate_material_chart};
use fs_evidence::ValidityDomain;
use fs_exec::{CancelGate, Cx};
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, NormalizedPack, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::gas::{GasSpec, GasState};
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::state_point::{MaterialPropertySelection, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point};
use fs_qty::{Density, Dims, DynViscosity, Pressure};
use fs_plate::{AssemblyOptions, EdgeSupport, PlateMesh, PlateRegion};
use fs_scenario::{PrestressedString, ThinPlate, gesture::{GestureSchedule, GestureTarget, GestureTrack, GestureValue}};
use fs_tribo::FrictionLaw;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate=CancelGate::new_clock_free();
    let pool=fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| f(&Cx::new(&gate,arena,fs_exec::StreamKey { seed:6,kernel_id:5,tile:0,iteration:0 },
        fs_exec::Budget::INFINITE,fs_exec::ExecMode::Deterministic),&gate))
}
fn material(rho:f64,g:f64,eta:f64,hi:f64)->ResolvedMaterialStatePoint { material_with(rho,g,eta,hi,&[]) }
fn material_with(rho:f64,g:f64,eta:f64,hi:f64,extra:&[(&str,Dims,f64)])->ResolvedMaterialStatePoint {
    let mut claims=ClaimSet::new(); let mut requirements=Vec::new();
    for (name,dims,value) in [("density",Density::DIMS,rho),("young_modulus",Pressure::DIMS,200e9),
        ("poisson_ratio",Dims::NONE,0.3),(ETA,DynViscosity::DIMS,10.0),
        ("G",Pressure::DIMS,g),("eta_s",DynViscosity::DIMS,eta)].into_iter().chain(extra.iter().copied()) {
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
        relaxation_bending:None,polarization_detune:0.0,moving_end:false,tension_loss_factor:None,bending_loss_factor:None}
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

/// A source stating SHEAR_VISCOSITY_REFERENCE_OMEGA [rad/s] beside its shear viscosity.
fn anchored(g:f64,eta:f64,reference:f64)->MaterialStringSource {
    MaterialStringSource::Homogeneous(with_uniform_circular_material_state(template(),0.0003,
        &material_with(7800.0,g,eta,200_000.0,&[(TORSION_REFERENCE,TORSION_REFERENCE_DIMS,reference)])).unwrap()
        .with_kelvin_voigt_bending_loss().unwrap())
}

// The receipt for the string band: every weight is positive, and G stays within 2 % of one from 100 Hz to 30 kHz
// (measured 0.9849..1.0169 in the numpy probe) and within 4.5 % over the whole 50 Hz-40 kHz fit band (0.9624..1.0169).
#[test]
fn the_string_loss_spectrum_is_flat_to_thirty_kilohertz() {
    let spectrum=ConstantLossSpectrum::fit_over(STRING_LOSS_SPECTRUM_BAND_RAD_S);
    assert!(spectrum.branches().all(|(_,g)|g>0.0),"{spectrum:?}");
    let extreme=|lo:f64,hi:f64| (0..2001).map(|k| core::f64::consts::TAU*lo*(hi/lo).powf(k as f64/2000.0))
        .map(|w| spectrum.shape(w)).fold((f64::MAX,f64::MIN),|(a,b),g|(a.min(g),b.max(g)));
    let inner=extreme(100.0,30_000.0); let band=extreme(50.0,40_000.0);
    eprintln!("G over 100 Hz-30 kHz {:.4}..{:.4}; over 50 Hz-40 kHz {:.4}..{:.4}",inner.0,inner.1,band.0,band.1);
    assert!(inner.0>0.98 && inner.1<1.02,"{inner:?}");
    assert!(band.0>0.955 && band.1<1.045,"{band:?}");
}

// A shear viscosity written as eta_s = eta_t G / omega_ref with omega_ref stated gives every twist mode the loss
// factor eta_t G(omega_n) (damping ratio half that); the same number read as Kelvin-Voigt gives eta_t omega_n / omega_ref.
// Frequencies and normalization do not change.
#[test]
fn a_stated_shear_viscosity_reference_gives_every_twist_mode_its_constant_loss_factor() {
    with_cx(|cx,gate| {
        let (g,eta_t,reference)=(80e9,0.0585,core::f64::consts::TAU*1000.0);
        let eta_s=eta_t*g/reference;
        let spectrum=ConstantLossSpectrum::fit_over(STRING_LOSS_SPECTRUM_BAND_RAD_S);
        let kelvin_voigt=raw(cx,gate,source(g,eta_s).with_torsion("G","eta_s",3).unwrap(),options()).unwrap();
        let constant=raw(cx,gate,anchored(g,eta_s,reference).with_torsion("G","eta_s",3).unwrap(),options()).unwrap();
        let (kv,cl)=(kelvin_voigt.system().components()[0].modes(),constant.system().components()[0].modes());
        assert_eq!((kv.len(),cl.len()),(7,7));
        for i in 0..4 { assert_eq!(kv[i],cl[i]); }
        assert_eq!(kelvin_voigt.system().port().left.shapes,constant.system().port().left.shapes);
        for n in 1..=3 {
            let omega=cl[3+n].angular_frequency_rad_s;
            assert_eq!(omega.to_bits(),kv[3+n].angular_frequency_rad_s.to_bits());
            let loss=2.0*cl[3+n].damping_ratio;
            assert!((loss/(eta_t*spectrum.shape(omega))-1.0).abs()<1e-12,"mode {n}: {loss}");
            assert!((2.0*kv[3+n].damping_ratio/(eta_t*omega/reference)-1.0).abs()<1e-12,"mode {n}");
        }
    });
}

// The law is claimed only inside its spectrum's fit band: a twist mode below 50 Hz refuses under constant loss, while
// the same string read as Kelvin-Voigt is admitted.
#[test]
fn constant_loss_torsion_refuses_a_twist_mode_outside_its_spectrum_band() {
    with_cx(|cx,gate| {
        // c_t = sqrt(1e7 / 7800) = 35.8 m/s puts the first twist mode of the 0.5 m span at 35.8 Hz.
        let (g,reference)=(1e7,core::f64::consts::TAU*1000.0);
        assert!(raw(cx,gate,source(g,100.0).with_torsion("G","eta_s",1).unwrap(),options()).is_ok());
        let refused=raw(cx,gate,anchored(g,100.0,reference).with_torsion("G","eta_s",1).unwrap(),options()).err()
            .expect("a 35.8 Hz constant-loss twist mode must refuse");
        assert!(refused.to_string().contains("exceed material support"),"{refused}");
    });
}

/// Pickering (1985) as committed: the scalar `<key>_<row>` of his whole-string table.
fn pickering(key:&str,row:&str)->f64 {
    let name=format!("{key}_{row}");
    include_str!("../../../data/matdb/seed-v1/violin-strings-whole-pickering-jcas1985/properties.tsv").lines()
        .map(|l| l.split('\t').collect::<Vec<_>>()).find(|f| f[0]=="scalar" && f[3]==name)
        .map(|f| f[4].parse().unwrap()).unwrap()
}
/// A material_strings string pack resolved as the example resolves a twisting string with the constant-loss anchor.
fn string_pack(bytes:&[u8])->ResolvedMaterialStatePoint {
    let pack=NormalizedPack::from_bytes(bytes).unwrap();
    let card=MaterialCard::assemble(MaterialStateId{chemistry:format!("source-pack:{}",pack.pack_id()),
        phase:"declared solid model".into(),process:"frozen source-pack realization".into(),revision:0},
        pack.claims().clone(),vec![]).unwrap();
    let requirements=[("density",Density::DIMS),("young_modulus",Pressure::DIMS),(ETA,DynViscosity::DIMS),
        ("shear_modulus",Pressure::DIMS),("kelvin_voigt_shear_viscosity",DynViscosity::DIMS),
        (TORSION_REFERENCE,TORSION_REFERENCE_DIMS)]
        .map(|(name,dims)| ScalarPropertyRequirement::try_new(name,dims,ScalarAdmissibility::Finite).unwrap());
    resolve_material_state_point(&card,&QueryPoint::new().with("T",293.15).unwrap().with("omega",1000.0).unwrap(),
        &requirements,MaterialPropertySelection::SingleClaimOnly).unwrap()
}

// The four Thomastik Dominant violin strings of Pickering (1985) as material_strings ships them
// (examples/material_strings/gen_dominant_string_packs.py). At his mean radius each pack holds his linear density and
// his string torsional stiffness on the 329 mm bed (both halves twist), and its Young's modulus is the least-squares
// inharmonicity of his whole-cent partials. On a violin span every twist mode carries the string's own loss factor at
// its own frequency: the pendulum-and-string decrement with the pendulum's 30.355 g cm/rad of gravity stiffness taken
// out, at the midpoint of the bracket the pendulum's unstated own loss leaves (0 up to all of the steel E's loss).
#[test]
fn the_dominant_packs_hold_pickerings_strings_and_their_measured_twist_loss() {
    use core::f64::consts::PI;
    let packs:[(&str,f64,&[u8]);4]=[
        ("s01",659.26,include_bytes!("../examples/material_strings/packs/violin-string-dominant-e-pickering1985-section.fsmatpk")),
        ("s12",440.0,include_bytes!("../examples/material_strings/packs/violin-string-dominant-a-pickering1985-section.fsmatpk")),
        ("s22",293.66,include_bytes!("../examples/material_strings/packs/violin-string-dominant-d-pickering1985-section.fsmatpk")),
        ("s32",196.0,include_bytes!("../examples/material_strings/packs/violin-string-dominant-g-pickering1985-section.fsmatpk"))];
    let (bed,gram_cm,pendulum,speaking,span)=(0.329,9.80665e-5,30.355,0.328,0.383);
    let system=|row:&str| pickering("damping_coefficient_printed_x100",row)/100.0/PI;
    let k_e=pickering("torsion_stiffness_printed_g_cm_per_rad","s01");
    let pendulum_loss=system("s01")*(k_e+pendulum)/pendulum;
    let spectrum=ConstantLossSpectrum::fit_over(STRING_LOSS_SPECTRUM_BAND_RAD_S);
    with_cx(|cx,gate| for (row,pitch,bytes) in packs {
        let m=string_pack(bytes);
        let v=|name:&str| m.property(name).unwrap().value_si();
        let r=0.25*(pickering("diameter_max",row)+pickering("diameter_min",row));
        let (mu,tension,k)=(pickering("linear_density",row),pickering("tension",row),
            pickering("torsion_stiffness_printed_g_cm_per_rad",row));
        assert!((v("density")*PI*r*r/mu-1.0).abs()<1e-8,"{row}");
        assert!((4.0*v("shear_modulus")*(PI*r.powi(4)/2.0)/bed/(k*gram_cm)-1.0).abs()<1e-8,"{row}");
        let squares=|b:f64| (2..=8).map(|n| {
            let cents=1200.0*((1.0+b*(n*n) as f64).sqrt()/(1.0+b).sqrt()).log2();
            (cents-pickering(&format!("partial_{n}_deviation_cents"),row)).powi(2)
        }).sum::<f64>();
        let b=PI*PI*v("young_modulus")*(PI*r.powi(4)/4.0)/(tension*bed*bed);
        assert!(squares(b)<squares(0.99*b) && squares(b)<squares(1.01*b),"{row}: B {b:e}");
        let eta_hi=system(row)*(k+pendulum)/k;
        let eta_lo=((system(row)*(k+pendulum)-pendulum_loss*pendulum)/k).max(0.0);
        let eta_t=v("kelvin_voigt_shear_viscosity")*v(TORSION_REFERENCE)/v("shear_modulus");
        assert!((eta_t/(0.5*(eta_lo+eta_hi))-1.0).abs()<1e-8,"{row}: {eta_t} in {eta_lo}..{eta_hi}");
        assert_eq!(v(ETA),0.0,"{row}: the bending loss is stated unknown, written 0");
        eprintln!("{row}: B {b:.4e} (rms {:.2} cents), twist loss factor {eta_t:.5} ({eta_lo:.5}..{eta_hi:.5})",
            (squares(b)/7.0).sqrt());
        // The il Cannone span at the open string's pitch from its speaking length.
        let string=PrestressedString{length_m:span,tension_n:mu*(2.0*speaking*pitch).powi(2),..template()};
        let source=MaterialStringSource::Homogeneous(with_uniform_circular_material_state(string,r,&m).unwrap()
            .with_kelvin_voigt_bending_loss().unwrap()).with_torsion("shear_modulus","kelvin_voigt_shear_viscosity",3).unwrap();
        let bank=raw(cx,gate,source,options()).unwrap();
        let modes=bank.system().components()[0].modes();
        assert_eq!(modes.len(),7);
        for n in 1..=3 {
            let omega=n as f64*PI/span*(v("shear_modulus")/v("density")).sqrt();
            assert!((modes[3+n].angular_frequency_rad_s/omega-1.0).abs()<2e-14,"{row} twist {n}");
            let loss=2.0*modes[3+n].damping_ratio;
            assert!((loss/(eta_t*spectrum.shape(omega))-1.0).abs()<1e-12,"{row} twist {n}: {loss}");
        }
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

fn rotational_link(body: ModalRotationAttachment) -> StringBodyRotationConnection {
    StringBodyRotationConnection {
        string: 0, station_fraction: 0.8, body,
        stiffness_n_m_per_rad: 0.001, damping_n_m_s_per_rad: 1e-7,
        rest_angle_rad: 0.0,
    }
}

fn radiating_rotation(cx: &Cx<'_>, gate: &CancelGate, g: f64) -> ReciprocalStringBank {
    let bank = radiating(cx, gate, g);
    let moving = &bank.source().strings()[0].connections[0].body;
    // Small-angle chord from an explicitly fixed body support to its actual
    // material-derived moving footprint, over a declared 0.1 m baseline.
    let fixed = ModalAttachment { component: moving.component, shapes: vec![0.0; moving.shapes.len()] };
    let rotation = ModalRotationAttachment::from_displacement_footprints(moving, &fixed, 0.1).unwrap();
    bank.with_rotational_body_connections(vec![rotational_link(rotation)], gate).unwrap()
}

#[test]
fn g3_torsional_bridge_radiates_after_release_and_shear_swap_changes_that_physical_response() {
    with_cx(|cx, gate| {
        let mut unconnected = radiating(cx, gate, 80e9);
        let mut original = radiating_rotation(cx, gate, 80e9);
        let mut changed = radiating_rotation(cx, gate, 20e9);
        assert_eq!(original.cards()[0].fundamental_hz(), changed.cards()[0].fundamental_hz());
        assert_eq!(original.cards()[0].zetas, changed.cards()[0].zetas);
        assert_eq!(original.system().components()[1].modes(), changed.system().components()[1].modes());
        assert_eq!(original.system().acceleration_pressure_weights(), unconnected.system().acceleration_pressure_weights());
        assert_eq!(original.source().rotational_body_connections(), changed.source().rotational_body_connections());
        let mut force = vec![0.0; original.system().mode_count()];
        let released = FrictionDrive { speed_m_s: 0.0, normal_force_n: 0.0 };
        let mut peak_pressure = 0.0_f64;
        let mut changed_pressure = 0.0_f64;
        let mut previous_energy = [0.0_f64; 2];
        for sample in 0..512 {
            // Directly excite ONLY one twist coordinate with a finite held
            // force, then release it. No bow can relay torque to the body here.
            force[4] = if sample < 32 { 0.02 } else { 0.0 };
            let baseline = *unconnected.step(&force, released).unwrap();
            let a = *original.step(&force, released).unwrap();
            let b = *changed.step(&force, released).unwrap();
            assert_eq!(baseline.observer_pressure_pa, 0.0);
            for (i, frame) in [a, b].iter().enumerate() {
                assert_eq!(frame.regime, ModalFrictionRegime::Released);
                assert_eq!(frame.drive_work_j, 0.0);
                assert_eq!(frame.friction_dissipation_j, 0.0);
                assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
                assert!(frame.network_dissipation_j >= 0.0);
                if sample >= 32 {
                    assert_eq!(frame.external_work_j, 0.0);
                    assert!(frame.network_energy_j <= previous_energy[i] + frame.energy_tolerance_j);
                }
                previous_energy[i] = frame.network_energy_j;
            }
            if sample >= 32 {
                peak_pressure = peak_pressure.max(a.observer_pressure_pa.abs());
                changed_pressure = changed_pressure.max((a.observer_pressure_pa - b.observer_pressure_pa).abs());
            }
        }
        assert!(peak_pressure > 1e-12, "released torque never reached physical body pressure: {peak_pressure:e}");
        assert!(changed_pressure > 1e-12, "shear swap never reached the released body: {changed_pressure:e}");
        assert!(original.system().components()[1].states().iter().any(|s| s.displacement_m_sqrt_kg != 0.0));
        assert!(unconnected.system().components()[1].states().iter().all(|s| s.displacement_m_sqrt_kg == 0.0));
    });
}

#[test]
fn g2_reciprocal_torque_matches_independent_two_oscillator_solution() {
    with_cx(|cx, gate| {
        let o = options();
        let receiver = |omega| ModalAcousticTimeModel::try_new(o.sample_rate_hz, vec![ModalAcousticMode {
            angular_frequency_rad_s: omega, damping_ratio: 0.0,
            pressure_per_modal_velocity: fs_math::c64::C64::ZERO,
        }], ModalAcousticTimeBudget::audible_reference()).unwrap();
        let src = source(80e9, 0.0).with_torsion("G", "eta_s", 1).unwrap();
        let mut bank = ReciprocalStringBank::new(cx, vec![MaterialStringMember {
            source: src, connections: vec![InteriorStringConnection {
                station_fraction: 0.8, body: ModalAttachment { component: 0, shapes: vec![2.0] },
                stiffness_n_m: 3000.0, damping_n_s_m: 0.0, rest_extension_m: 0.0,
            }],
        }], gas(), vec![receiver(450.0), receiver(300.0)], law(), o, gate).unwrap();
        let mut link = rotational_link(ModalRotationAttachment { component: 1, angular_shapes: vec![100.0] });
        link.damping_n_m_s_per_rad = 0.0;
        bank = bank.with_rotational_body_connections(vec![link.clone()], gate).unwrap();
        let MaterialStringSource::WithTorsion(twist) = &bank.source().strings()[0].source else { panic!("missing source") };
        let angular = (core::f64::consts::PI * link.station_fraction).sin()
            / (twist.section().rotary_inertia_kg_m() * bank.cards()[0].length_m * 0.5).sqrt();
        let w = bank.system().components()[0].modes()[4].angular_frequency_rad_s;
        // Independent physical angular spring: K += k*[phi_s,-phi_b]^T*[phi_s,-phi_b].
        // Its only driven subspace is the one twist and one receiving body mode.
        let aa = w * w + link.stiffness_n_m_per_rad * angular * angular;
        let dd = 300.0_f64.powi(2) + link.stiffness_n_m_per_rad * 100.0_f64.powi(2);
        let off = -link.stiffness_n_m_per_rad * angular * 100.0;
        let centre = 0.5 * (aa + dd);
        let spread = (0.5 * (aa - dd)).hypot(off);
        let hi = centre + spread;
        let lo = (aa * dd - off * off) / hi;
        let mut force = vec![0.0; bank.system().mode_count()];
        force[4] = 1.0;
        let released = FrictionDrive { speed_m_s: 0.0, normal_force_n: 0.0 };
        for sample in 1..=128 {
            let frame = *bank.step(&force, released).unwrap();
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
            let time = f64::from(sample) / 48_000.0;
            let g = |lambda: f64| (1.0 - (lambda.sqrt() * time).cos()) / lambda;
            let h = |lambda: f64| (lambda.sqrt() * time).sin() / lambda.sqrt();
            let apply = |f: &dyn Fn(f64) -> f64| {
                let gradient = (f(hi) - f(lo)) / (hi - lo);
                [f(lo) + gradient * (aa - lo), gradient * off]
            };
            let q = apply(&g);
            let v = apply(&h);
            let actual = [bank.system().components()[0].states()[4], bank.system().components()[2].states()[0]];
            for i in 0..2 {
                assert!((actual[i].displacement_m_sqrt_kg - q[i]).abs() <= 1e-16 + 1e-7 * q[i].abs(),
                    "sample {sample}, coordinate {i}: q={:?}, expected={q:?}", actual);
                assert!((actual[i].velocity_m_sqrt_kg_per_s - v[i]).abs() <= 1e-12 + 1e-7 * v[i].abs());
            }
            assert!(bank.system().components()[0].states()[..4].iter().all(|s| s.displacement_m_sqrt_kg == 0.0));
            assert_eq!(bank.system().components()[1].states()[0].displacement_m_sqrt_kg, 0.0);
        }
    });
}

#[test]
fn g0_rotational_bridge_rejects_missing_sources_bad_units_dimensions_and_live_rebinding() {
    with_cx(|cx, gate| {
        let make = || raw(cx, gate, source(80e9, 1000.0).with_torsion("G", "eta_s", 1).unwrap(), options()).unwrap();
        let good = rotational_link(ModalRotationAttachment { component: 0, angular_shapes: vec![20.0] });
        for case in 0..11 {
            let mut bad = good.clone();
            match case {
                0 => bad.string = 1,
                1 => bad.station_fraction = 0.0,
                2 => bad.station_fraction = 1.0,
                3 => bad.body.component = 1,
                4 => bad.body.angular_shapes = vec![],
                5 => bad.body.angular_shapes[0] = 0.0,
                6 => bad.body.angular_shapes[0] = f64::NAN,
                7 => bad.stiffness_n_m_per_rad = -1.0,
                8 => bad.damping_n_m_s_per_rad = f64::INFINITY,
                9 => bad.rest_angle_rad = f64::NAN,
                _ => { bad.stiffness_n_m_per_rad = 0.0; bad.damping_n_m_s_per_rad = 0.0; }
            }
            assert!(make().with_rotational_body_connections(vec![bad], gate).is_err(), "case {case}");
        }
        assert!(raw(cx, gate, source(80e9, 1000.0), options()).unwrap()
            .with_rotational_body_connections(vec![good.clone()], gate).is_err());
        let mut small = options();
        small.network.max_connections = 1;
        assert!(raw(cx, gate, source(80e9, 1000.0).with_torsion("G", "eta_s", 1).unwrap(), small).unwrap()
            .with_rotational_body_connections(vec![good.clone()], gate).is_err());
        assert!(make().with_rotational_body_connections(vec![good.clone()], gate).unwrap()
            .with_rotational_body_connections(vec![good.clone()], gate).is_err());
        let mut playing = make();
        playing.step(&vec![0.0; playing.system().mode_count()], drive()).unwrap();
        assert!(playing.with_rotational_body_connections(vec![good.clone()], gate).is_err());
        let pending = make();
        gate.request();
        assert!(pending.with_rotational_body_connections(vec![good], gate).is_err());
    });
}

fn spatial_rotation_bank(cx: &Cx<'_>, gate: &CancelGate, body_density: f64) -> ReciprocalStringBank {
    let mesh = PlateMesh::rectangle(0.4, 0.3, 5, 4);
    let region = PlateRegionMaterial {
        region: PlateRegion { name: "body".into(), triangle_indices: (0..mesh.tris.len()).collect() },
        material: material(body_density, 80e9, 1000.0, 200_000.0),
        model: PlateMaterialModel::Isotropic,
        thickness_constraint: PlateThicknessConstraint::FixedThickness(0.001), material_angle_rad: 0.0,
    };
    let point = |at: usize| { let mut row = vec![0.0; mesh.nodes.len()]; row[at] = 1.0; row };
    let weights = vec![point(14), point(15)];
    let body = SharedMaterialPlateChartBody {
        specimen: compile_plate_material_chart(mesh, PlateMesh::rectangle_boundary(5, 4), vec![region]).unwrap(),
        reduction: PlateChartRadiation {
            assembly: AssemblyOptions { pretension: 0.0, support: EdgeSupport::SimplySupported },
            eigenvalue_window: (1.0, 1e8), n_modes: 4, damping_ratio: 0.0,
            unit_force_weights: weights[0].clone(),
        },
        bridges: vec![SharedPlateBridge { station_fraction: 0.8, stiffness_n_m: 3000.0,
            damping_n_s_m: 0.1, rest_extension_m: 0.0 }; 2],
        max_plate_dofs: 90, max_plate_elements: 40,
        modal_budget: ModalAcousticTimeBudget::audible_reference(), listener_m: 1.0,
    };
    ReciprocalStringBank::from_material_plate_chart_with_force_ports(cx,
        vec![source(80e9, 1000.0).with_torsion("G", "eta_s", 1).unwrap(); 2],
        gas(), body, PlateBridgeFootprints { weights, max_projection_terms: 100_000 }, None,
        law(), options(), gate).unwrap()
}

#[test]
fn g3_material_body_swap_rebuilds_rotation_from_original_spatial_footprints_in_joint_bowing() {
    with_cx(|cx, gate| {
        let mut rows = Vec::new();
        let mut signals = Vec::new();
        for density in [7800.0, 15_600.0] {
            let bank = spatial_rotation_bank(cx, gate, density);
            let rotation = bank.source().rotation_between_bridges(
                StringBridgeLocation { string: 0, connection: 0 },
                StringBridgeLocation { string: 1, connection: 0 }, 0.08,
            ).unwrap();
            let bridge = |string: usize| &bank.source().strings()[string].connections[0].body.shapes;
            let expected: Vec<_> = bridge(0).iter().zip(bridge(1)).map(|(p, n)| (p - n) / 0.08).collect();
            assert_eq!(rotation.angular_shapes, expected);
            rows.push(rotation.angular_shapes.clone());
            let old_weights = bank.system().acceleration_pressure_weights().unwrap().to_vec();
            let old_body = bank.system().components()[2].modes().to_vec();
            let link = rotational_link(rotation);
            let mut bank = bank.with_rotational_body_connections(vec![link.clone()], gate).unwrap()
                .with_additional_bows(vec![MaterialBowContact {
                    location: StringBowLocation { string: 1, station_fraction: 0.31 },
                    law: law(), friction: options().friction,
                }], MultiFrictionConfig { max_ports: 4, max_sweeps: 128, max_setup_terms: 16_384 }, gate).unwrap();
            assert_eq!(bank.source().bank().rotational_body_connections(), &[link]);
            assert!(bank.source().bank().material_chart_bridges().is_some());
            assert_eq!(bank.system().acceleration_pressure_weights().unwrap(), old_weights);
            assert_eq!(bank.system().components()[2].modes(), old_body);
            let zero = vec![0.0; bank.system().mode_count()];
            let mut signal = Vec::new();
            for sample in 0..128 {
                let drives = [
                    FrictionDrive { speed_m_s: 0.01, normal_force_n: if sample < 32 { 0.001 } else { 0.0 } },
                    FrictionDrive { speed_m_s: -0.01, normal_force_n: if sample < 32 { 0.001 } else { 0.0 } },
                ];
                let frame = bank.step(&zero, &drives).unwrap();
                assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
                signal.push(frame.observer_pressure_pa);
            }
            signals.push(signal);
        }
        assert!(rows[0].iter().any(|x| x.abs() > 1e-6));
        for (a, b) in rows[0].iter().zip(&rows[1]) {
            // Doubling only body density changes mass normalization by sqrt(2).
            // Mode sign remains arbitrary; the physical row scales in magnitude.
            assert!((a.abs() / 2.0_f64.sqrt() - b.abs()).abs() < 1e-8 * a.abs().max(1.0));
        }
        assert!(signals[0][32..].iter().any(|x| x.abs() > 1e-12));
        assert!(signals[0][32..].iter().zip(&signals[1][32..]).any(|(a, b)| (a - b).abs() > 1e-12));
    });
}

#[test]
fn g0_rotation_geometry_refusals_and_stress_free_angle_have_physical_units() {
    let positive = ModalAttachment { component: 0, shapes: vec![2.0, -1.0] };
    let negative = ModalAttachment { component: 0, shapes: vec![-1.0, 1.0] };
    let rotation = ModalRotationAttachment::from_displacement_footprints(&positive, &negative, 0.2).unwrap();
    assert_eq!(rotation.angular_shapes, vec![15.0, -10.0]);
    for separation in [0.0, -0.1, f64::NAN, f64::INFINITY] {
        assert!(ModalRotationAttachment::from_displacement_footprints(&positive, &negative, separation).is_err());
    }
    assert!(ModalRotationAttachment::from_displacement_footprints(&positive, &positive, 0.2).is_err());
    let mut other = negative.clone(); other.component = 1;
    assert!(ModalRotationAttachment::from_displacement_footprints(&positive, &other, 0.2).is_err());
    other = negative.clone(); other.shapes.pop();
    assert!(ModalRotationAttachment::from_displacement_footprints(&positive, &other, 0.2).is_err());
    with_cx(|cx, gate| {
        let mut link = rotational_link(ModalRotationAttachment { component: 0, angular_shapes: vec![20.0] });
        link.rest_angle_rad = 0.01;
        let mut bank = raw(cx, gate, source(80e9, 1000.0).with_torsion("G", "eta_s", 1).unwrap(), options()).unwrap()
            .with_rotational_body_connections(vec![link.clone()], gate).unwrap();
        let initial = bank.system().total_energy_j().unwrap();
        let expected = 0.5 * link.stiffness_n_m_per_rad * link.rest_angle_rad.powi(2);
        assert!((initial / expected - 1.0).abs() < 1e-14);
        let zero = vec![0.0; bank.system().mode_count()];
        let released = FrictionDrive { speed_m_s: 0.0, normal_force_n: 0.0 };
        let mut loss = 0.0;
        for _ in 0..128 {
            let frame = bank.step(&zero, released).unwrap();
            loss += frame.network_dissipation_j;
            assert_eq!(frame.external_work_j, 0.0);
            assert_eq!(frame.drive_work_j, 0.0);
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        }
        assert!(loss > 0.0);
        assert!((bank.system().total_energy_j().unwrap() + loss - initial).abs() < 1e-13);
    });
}

#[path = "material_string_torsion/polarization.rs"]
mod polarization;

#[path = "material_string_torsion/fluid.rs"]
mod fluid;
