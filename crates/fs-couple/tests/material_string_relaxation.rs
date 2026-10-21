//! Source substitutions drive the actual reciprocal string/body/friction state.
#[path = "material_string_nonlinear/support.rs"]
#[allow(dead_code)]
mod support;
use support::{with_cx,template,gas};
use fs_couple::bowed_string::runtime::reciprocal::{InteriorStringConnection,bank::{
    MaterialStringMember,MaterialStringSource,ReciprocalStringBank,StringBankOptions,StringBowLocation,
    polarization::TransverseDirections,relaxation::StringBendingRelaxationSelection,
}};
use fs_couple::modal_acoustic_time::{ModalAcousticMode,ModalAcousticTimeBudget,ModalAcousticTimeModel};
use fs_couple::render::schedule::force::coupled::{ModalAttachment,ModalCouplingConfig};
use fs_couple::render::schedule::force::coupled::contact::friction::{FrictionDrive,ModalFrictionConfig};
use fs_couple::string_specimen::{StringPrestress,KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY as ETA};
use fs_exec::{Cx,CancelGate};
use fs_material::section::{ConcentricElasticSection,ConcentricMaterialLayer};
use fs_material::state_point::{ResolvedMaterialStatePoint,ScalarPropertyRequirement,
    ScalarAdmissibility,MaterialPropertySelection,resolve_material_state_point};
use fs_matdb::{ClaimSet,PropertyClaim,PropertyKey,PropertyValue,Provenance,MaterialCard,
    MaterialStateId,QueryPoint,UncertaintyModel,InterpolationPolicy};
use fs_evidence::ValidityDomain;
use fs_math::c64::C64;
use fs_qty::{Dims,Density,Pressure,DynViscosity};
use fs_tribo::FrictionLaw;

fn material(delta: f64,tau: f64,band: f64) -> ResolvedMaterialStatePoint {
    let mut claims=ClaimSet::new(); let mut requirements=Vec::new();
    for (name,dims,value) in [("density",Density::DIMS,7800.0),
        ("young_modulus",Pressure::DIMS,200e9),(ETA,DynViscosity::DIMS,0.0),
        ("G",Pressure::DIMS,60e9),("eta_s",DynViscosity::DIMS,1000.0),
        ("bending_delta",Pressure::DIMS,delta),("bending_time",Dims([0,0,1,0,0,0]),tau)] {
        let high=if name.starts_with("bending_") { band } else {200_000.0};
        claims.insert_claim(PropertyClaim {key:PropertyKey::new(name,dims),
            value:PropertyValue::Scalar {value,dims},
            validity:ValidityDomain::unconstrained().with("T",280.0,320.0).with("omega",0.0,high),
            uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::ConstantWithinValidity,
            observations:vec![],provenance:Provenance {source:"manufactured layered bending spectrum".into(),
                license:"CC0-1.0".into(),artifact:None}}).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(name,dims,ScalarAdmissibility::Finite).unwrap());
    }
    let card=MaterialCard::assemble(MaterialStateId {chemistry:"synthetic".into(),phase:"solid".into(),
        process:"test".into(),revision:0},claims,vec![]).unwrap();
    resolve_material_state_point(&card,&QueryPoint::new().with("T",293.15).unwrap().with("omega",1000.0).unwrap(),
        &requirements,MaterialPropertySelection::SingleClaimOnly).unwrap()
}
fn bank(cx:&Cx<'_>,gate:&CancelGate,delta:f64,tau:f64,band:f64,modes:usize,count:usize) -> ReciprocalStringBank {
    let section=ConcentricElasticSection::try_new(0.0,vec![
        ConcentricMaterialLayer {outer_radius_m:0.0002,material:material(delta,0.002,band)},
        ConcentricMaterialLayer {outer_radius_m:0.0003,material:material(delta,tau,band)},
    ]).unwrap().with_kelvin_voigt_bending(ETA).unwrap();
    let source=MaterialStringSource::Concentric {template:template(modes),section,
        prestress:StringPrestress::FixedTension(20.0)}.with_torsion("G","eta_s",1).unwrap()
        .with_transverse_polarizations(TransverseDirections::try_new([0.6,0.8],[1.0,0.0]).unwrap()).unwrap();
    let members=(0..count).map(|_|MaterialStringMember {source:source.clone(),connections:vec![
        InteriorStringConnection {station_fraction:0.8,body:ModalAttachment {component:0,shapes:vec![1.0]},
            stiffness_n_m:3000.0,damping_n_s_m:0.1,rest_extension_m:0.0},
    ]}).collect();
    let body=ModalAcousticTimeModel::try_new(48_000,vec![ModalAcousticMode {
        angular_frequency_rad_s:500.0,damping_ratio:0.02,
        pressure_per_modal_velocity:C64::new(1.0,0.0),
    }],ModalAcousticTimeBudget::audible_reference()).unwrap();
    ReciprocalStringBank::new(cx,members,gas(),vec![body],FrictionLaw::Coulomb {static_mu:0.8,kinetic_mu:0.4},
        StringBankOptions {sample_rate_hz:48_000,bow:StringBowLocation {string:0,station_fraction:0.2},
            string_budget:ModalAcousticTimeBudget::audible_reference(),
            network:ModalCouplingConfig {max_modes:64,max_connections:64,max_setup_terms:1_048_576,
                nyquist_guard_fraction:0.9,maximum_total_energy_j:100.0,maximum_abs_pressure_pa:1e6,
                maximum_abs_connection_force_n:1000.0,solve_relative_tolerance:1e-10,
                energy_absolute_tolerance_j:1e-12,energy_relative_tolerance:1e-8},
            friction:ModalFrictionConfig {max_iterations:96,maximum_force_n:100.0,
                maximum_slip_speed_m_s:10.0,force_absolute_tolerance_n:1e-10,
                force_relative_tolerance:1e-9,sticking_tolerance_m:1e-12}},gate).unwrap()
}
fn selection() -> Vec<StringBendingRelaxationSelection> {
    vec![StringBendingRelaxationSelection {string:0,
        properties:vec![("bending_delta".into(),"bending_time".into())]}]
}
fn render(bank:&mut ReciprocalStringBank) -> Vec<f64> {
    let forces=vec![0.0;bank.system().mode_count()];
    (0..128).map(|_| {
        let frame=bank.step(&forces,FrictionDrive {speed_m_s:0.3,normal_force_n:0.001}).unwrap();
        assert!(frame.energy_residual_j.abs()<=frame.energy_tolerance_j);
        frame.observer_pressure_pa
    }).collect()
}

#[test]
fn independent_layer_memory_reaches_both_planes_shared_body_and_silent_string() {
    with_cx(|cx,gate| {
        let mut first=bank(cx,gate,40e9,0.004,200_000.0,2,2)
            .with_material_string_bending_relaxation(selection(),None,16,gate).unwrap();
        let mut changed=bank(cx,gate,40e9,0.012,200_000.0,2,2)
            .with_material_string_bending_relaxation(selection(),None,16,gate).unwrap();
        let image=first.source().string_bending_relaxation().unwrap();
        assert_eq!(image.channels().len(),8);
        assert_eq!(image.spectra()[0][0].terms[0].1,0.002);
        assert_eq!(image.spectra()[0][1].terms[0].1,0.004);
        assert_eq!(changed.source().string_bending_relaxation().unwrap().spectra()[0][1].terms[0].1,0.012);
        let channels=image.channels();
        let c=channels[0].coupling_per_s;
        assert_eq!(channels[2].coupling_per_s,c);
        assert!((channels[1].coupling_per_s/c-4.0).abs()<1e-13);
        let card=&first.cards()[0];
        let expected=40e9*image.sections()[0].layer_integrals()[0].second_moment_m4()
            *(core::f64::consts::PI/card.length_m).powi(4)/card.linear_density_kg_m;
        assert!((c*c-expected).abs()<1e-12*expected);
        for channel in channels { assert!(channel.mode<4); } // Twist remains independent.
        for (a,b) in first.system().components().iter().zip(changed.system().components()) {
            assert_eq!(a.modes(),b.modes()); // No elastic retuning or changed scalar loss.
        }
        let p=render(&mut first); let q=render(&mut changed);
        assert!(p.iter().any(|x|x.abs()>1e-15));
        assert!(p.iter().zip(q).any(|(p,q)|(p-q).abs()>1e-18));
        assert_ne!(first.system().components()[1].states(),changed.system().components()[1].states());
        assert!(first.system().material_energy_j().unwrap()>0.0);
        assert!(first.system().coupled_frame().unwrap().material_dissipation_j>0.0);
    });
}

#[test]
fn explicit_zero_is_elastic_and_invalid_or_stale_selections_refuse() {
    with_cx(|cx,gate| {
        let make=||bank(cx,gate,0.0,0.004,200_000.0,1,1);
        let mut original=make();
        let mut zero=make().with_material_string_bending_relaxation(selection(),None,16,gate).unwrap();
        assert!(zero.source().string_bending_relaxation().is_some());
        assert!(zero.system().material_memory().is_empty());
        assert_eq!(render(&mut original),render(&mut zero));
        for (delta,tau,band,limit) in [(40e9,0.004,200_000.0,1),
            (-1.0,0.004,200_000.0,16),(40e9,0.0,200_000.0,16),
            (40e9,1e-7,200_000.0,16),(40e9,0.004,8_000.0,16)] {
            assert!(bank(cx,gate,delta,tau,band,1,1)
                .with_material_string_bending_relaxation(selection(),None,limit,gate).is_err());
        }
        let mut missing=selection();missing[0].properties[0].0="missing".into();
        assert!(make().with_material_string_bending_relaxation(missing,None,16,gate).is_err());
        let mut duplicate=selection();duplicate.push(duplicate[0].clone());
        assert!(make().with_material_string_bending_relaxation(duplicate,None,16,gate).is_err());
        assert!(zero.with_material_string_bending_relaxation(selection(),None,16,gate).is_err());
    });
}

#[test]
fn gas_mass_enters_projection_and_axial_extension_retains_bending_history() {
    use fs_couple::air_path::cylinder::CylinderFlowResolution;
    use fs_couple::bowed_string::runtime::reciprocal::bank::nonlinear::longitudinal::LongitudinalStringResolution;
    use fs_couple::render::schedule::force::coupled::contact::friction::nonlinear::NonlinearModalConfig;
    with_cx(|cx,gate| {
        let make=||bank(cx,gate,40e9,0.004,200_000.0,1,1);
        let gas_resolution=CylinderFlowResolution {minimum_wavenumber_radius:0.001,
            maximum_wavenumber_radius:30.0,intervals:4};
        let dry=make().with_material_string_bending_relaxation(selection(),None,16,gate).unwrap();
        let loaded=make().with_cylinder_gas_loading(gas_resolution,gate).unwrap()
            .with_material_string_bending_relaxation(selection(),None,16,gate).unwrap();
        let a=dry.source().string_bending_relaxation().unwrap().channels()[0].coupling_per_s;
        let b=loaded.source().string_bending_relaxation().unwrap().channels()[0].coupling_per_s;
        let ratio=(dry.cards()[0].linear_density_kg_m/loaded.cards()[0].linear_density_kg_m).sqrt();
        assert!((b/a-ratio).abs()<1e-12);
        assert!(dry.with_cylinder_gas_loading(gas_resolution,gate).is_err());
        let memory=loaded.system().material_memory().to_vec();
        let mut axial=loaded.with_longitudinal_waves(vec![LongitudinalStringResolution {
            axial_modes:vec![1],maximum_slope:0.2,maximum_strain:0.01}],NonlinearModalConfig {
            max_iterations:64,maximum_generalized_force:1e5,
            force_absolute_tolerance:1e-11,force_relative_tolerance:1e-10},gate).unwrap();
        assert_eq!(axial.system().material_memory(),memory);
        assert_eq!(axial.source().string_bending_relaxation().unwrap().channels().len(),4);
        assert!(axial.system().material_relaxation_branches().iter()
            .all(|branch|branch.projection.len()==2*axial.system().mode_count()));
        assert!(render(&mut axial).iter().all(|p|p.is_finite()));
    });
}
