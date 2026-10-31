//! Source substitutions drive the actual reciprocal string/body/friction state.
#[path = "material_string_nonlinear/support.rs"]
#[allow(dead_code)]
mod support;
use support::{with_cx,template,gas};
use fs_couple::bowed_string::runtime::reciprocal::{InteriorStringConnection,bank::{
    MaterialStringMember,MaterialStringSource,ReciprocalStringBank,StringBankOptions,StringBowLocation,
    polarization::TransverseDirections,relaxation::{StringBendingRelaxationSelection,
        LayeredStringBendingRelaxationSelection},
}};
use fs_couple::modal_acoustic_time::{ModalAcousticMode,ModalAcousticTimeBudget,ModalAcousticTimeModel};
use fs_couple::render::schedule::force::coupled::{ModalAttachment,ModalCouplingConfig};
use fs_couple::render::schedule::force::coupled::contact::friction::{FrictionDrive,ModalFrictionConfig};
use fs_couple::string_specimen::{StringPrestress,KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY as ETA};
use fs_exec::{Cx,CancelGate};
use fs_material::section::{ConcentricElasticSection,ConcentricMaterialLayer};
use fs_material::section::relaxation::{LayerBendingRelaxationProperties,LayerRelaxationProperties};
use fs_material::state_point::{ResolvedMaterialStatePoint,ScalarPropertyRequirement,
    ScalarAdmissibility,MaterialPropertySelection,resolve_material_state_point};
use fs_matdb::{ClaimSet,PropertyClaim,PropertyKey,PropertyValue,Provenance,MaterialCard,
    MaterialStateId,QueryPoint,UncertaintyModel,InterpolationPolicy};
use fs_evidence::ValidityDomain;
use fs_math::c64::C64;
use fs_qty::{Dims,Density,Pressure,DynViscosity};
use fs_tribo::FrictionLaw;

fn material(delta: f64,tau: f64,band: f64) -> ResolvedMaterialStatePoint {
    material_with_arms(&[("bending_delta","bending_time",delta,tau)],band)
}
fn material_with_arms(arms: &[(&str,&str,f64,f64)],band: f64) -> ResolvedMaterialStatePoint {
    let mut claims=ClaimSet::new(); let mut requirements=Vec::new();
    let mut values=vec![("density",Density::DIMS,7800.0),
        ("young_modulus",Pressure::DIMS,200e9),
        ("equilibrium_young_modulus",Pressure::DIMS,200e9),(ETA,DynViscosity::DIMS,0.0),
        ("G",Pressure::DIMS,60e9),("eta_s",DynViscosity::DIMS,1000.0)];
    for &(modulus,time,delta,tau) in arms {
        values.extend([(modulus,Pressure::DIMS,delta),(time,Dims([0,0,1,0,0,0]),tau)]);
    }
    for (name,dims,value) in values {
        let high=if arms.iter().any(|(modulus,time,_,_)|name==*modulus || name==*time) { band } else {200_000.0};
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
    bank_from_materials(cx,gate,[material(delta,0.002,band),material(delta,tau,band)],modes,count)
}
fn bank_from_materials(cx:&Cx<'_>,gate:&CancelGate,materials:[ResolvedMaterialStatePoint;2],
    modes:usize,count:usize) -> ReciprocalStringBank {
    bank_with_sources(cx,gate,vec![source_from_materials(materials,modes);count]).unwrap()
}
fn source_from_materials(materials:[ResolvedMaterialStatePoint;2],modes:usize) -> MaterialStringSource {
    let [core,sleeve]=materials;
    let section=ConcentricElasticSection::try_new(0.0,vec![
        ConcentricMaterialLayer {outer_radius_m:0.0002,material:core},
        ConcentricMaterialLayer {outer_radius_m:0.0003,material:sleeve},
    ]).unwrap().with_kelvin_voigt_bending(ETA).unwrap();
    MaterialStringSource::Concentric {template:template(modes),section,
        prestress:StringPrestress::FixedTension(20.0)}.with_torsion("G","eta_s",1).unwrap()
        .with_transverse_polarizations(TransverseDirections::try_new([0.6,0.8],[1.0,0.0]).unwrap()).unwrap()
}
fn bank_with_sources(cx:&Cx<'_>,gate:&CancelGate,sources:Vec<MaterialStringSource>)
    -> Result<ReciprocalStringBank,fs_couple::bowed_string::runtime::reciprocal::ReciprocalStringError>
{
    let members=sources.into_iter().map(|source|MaterialStringMember {source,connections:vec![
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
                force_relative_tolerance:1e-9,sticking_tolerance_m:1e-12}},gate)
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
            (40e9,0.004,8_000.0,16)] {
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

fn original_source(delta:f64,tau:f64,band:f64,prony:bool) -> MaterialStringSource {
    use fs_couple::string_specimen::with_uniform_circular_material_state;
    use fs_scenario::acoustic::BendingRelaxationProperties;
    let specimen=with_uniform_circular_material_state(template(2),0.0003,
        &material(delta,tau,band)).unwrap();
    let specimen=if prony {
        specimen.with_prony_bending_loss(&[BendingRelaxationProperties {
            modulus:"bending_delta".into(),relaxation_time:"bending_time".into(),
        }])
    } else { specimen.with_kelvin_voigt_bending_loss() }.unwrap();
    MaterialStringSource::Homogeneous(specimen).with_torsion("G","eta_s",1).unwrap()
        .with_transverse_polarizations(TransverseDirections::try_new([0.6,0.8],[1.0,0.0]).unwrap()).unwrap()
}

#[test]
fn original_prony_material_reaches_pressure_and_sympathetic_string_without_reselection() {
    with_cx(|cx,gate| {
        let make=|tau|bank_with_sources(cx,gate,vec![
            original_source(40e9,tau,200_000.0,true),
            original_source(40e9,0.004,200_000.0,false),
        ]).unwrap();
        let mut a=make(0.004); let mut b=make(0.012);
        assert_eq!(a.system().material_memory().len(),4);
        assert!(a.source().string_bending_relaxation().is_none());
        for (a,b) in a.system().components().iter().zip(b.system().components()) {
            assert_eq!(a.modes(),b.modes());
        }
        let p=render(&mut a); let q=render(&mut b);
        assert!(p.iter().any(|x|x.abs()>1e-15));
        assert!(p.iter().zip(q).any(|(a,b)|(a-b).abs()>1e-18));
        assert_ne!(a.system().components()[1].states(),b.system().components()[1].states());
        assert!(a.system().material_memory().iter().any(|z|*z!=0.0));
        assert!(a.system().coupled_frame().unwrap().material_dissipation_j>0.0);
    });
}

#[test]
fn original_and_explicit_string_memory_share_gas_axial_and_friction_owners() {
    use fs_couple::air_path::cylinder::CylinderFlowResolution;
    use fs_couple::bowed_string::runtime::reciprocal::bank::nonlinear::longitudinal::LongitudinalStringResolution;
    use fs_couple::render::schedule::force::coupled::contact::friction::nonlinear::NonlinearModalConfig;
    with_cx(|cx,gate| {
        let make=||bank_with_sources(cx,gate,vec![
            original_source(40e9,0.004,200_000.0,true),
            original_source(20e9,0.008,200_000.0,false),
        ]).unwrap();
        let dry=make();
        let loaded=make().with_cylinder_gas_loading(CylinderFlowResolution {
            minimum_wavenumber_radius:0.001,maximum_wavenumber_radius:30.0,intervals:4,
        },gate).unwrap();
        let dry_k=dry.system().material_relaxation_branches()[0].stiffness;
        let wet_k=loaded.system().material_relaxation_branches()[0].stiffness;
        assert!((wet_k/dry_k-dry.cards()[0].linear_density_kg_m/loaded.cards()[0].linear_density_kg_m).abs()<1e-12);
        let before=loaded.system().material_memory().to_vec();
        let mut other=selection(); other[0].string=1;
        let joined=loaded.with_material_string_bending_relaxation(other,None,8,gate).unwrap();
        assert_eq!(&joined.system().material_memory()[..4],before);
        assert_eq!(joined.system().material_memory().len(),8);
        assert_eq!(joined.source().string_bending_relaxation().unwrap().channels().len(),4);
        let states=joined.system().material_memory().to_vec();
        let resolution=LongitudinalStringResolution {
            axial_modes:vec![1],maximum_slope:0.2,maximum_strain:0.01,
        };
        let mut axial=joined.with_longitudinal_waves(vec![resolution.clone(),resolution],
            NonlinearModalConfig {max_iterations:64,maximum_generalized_force:1e5,
                force_absolute_tolerance:1e-11,force_relative_tolerance:1e-10},gate).unwrap();
        assert_eq!(axial.system().material_memory(),states);
        assert!(axial.system().material_relaxation_branches().iter()
            .all(|b|b.projection.len()==2*axial.system().mode_count()));
        assert!(render(&mut axial).iter().all(|p|p.is_finite()));
        assert!(axial.system().material_memory()[..4].iter().any(|z|*z!=0.0));
        assert!(axial.system().material_memory()[4..].iter().any(|z|*z!=0.0));
        assert!(make().with_material_string_bending_relaxation(selection(),None,8,gate).is_err());
    });
}

#[test]
fn original_zero_memory_preserves_elastic_trajectory_and_source_limits_refuse() {
    with_cx(|cx,gate| {
        let make=|prony|bank_with_sources(cx,gate,
            vec![original_source(0.0,0.004,200_000.0,prony)]).unwrap();
        let mut elastic=make(false); let mut zero=make(true);
        assert!(zero.system().material_memory().is_empty());
        assert_eq!(render(&mut elastic),render(&mut zero));
        assert!(bank_with_sources(cx,gate,
            vec![original_source(40e9,0.004,8_000.0,true)]).is_err());
    });
}

#[test]
fn uniform_selectors_lower_to_the_same_layered_operator_and_pressure_bits() {
    with_cx(|cx,gate| {
        let make=||bank(cx,gate,40e9,0.004,200_000.0,2,2);
        let layers=vec![LayerBendingRelaxationProperties {
            equilibrium_modulus:"young_modulus".into(),branches:selection()[0].properties.clone(),
        };2];
        let mut uniform=make().with_material_string_bending_relaxation(selection(),None,16,gate).unwrap();
        let mut explicit=make().with_material_layered_string_bending_relaxation(vec![
            LayeredStringBendingRelaxationSelection {string:0,layers:layers.clone()},
        ],None,16,gate).unwrap();
        assert_eq!(uniform.source().string_bending_relaxation().unwrap().channels(),
            explicit.source().string_bending_relaxation().unwrap().channels());
        assert_eq!(explicit.source().string_bending_relaxation().unwrap().layered_selections()[0].layers,layers);
        let bits=|pressure:Vec<f64>|pressure.into_iter().map(f64::to_bits).collect::<Vec<_>>();
        assert_eq!(bits(render(&mut uniform)),bits(render(&mut explicit)));
    });
}

#[test]
fn elastic_core_and_independent_two_arm_coating_reach_shared_body_and_pressure() {
    with_cx(|cx,gate| {
        let layers=vec![
            LayerBendingRelaxationProperties {equilibrium_modulus:"young_modulus".into(),branches:vec![]},
            LayerBendingRelaxationProperties {equilibrium_modulus:"young_modulus".into(),branches:vec![
                ("coating_fast".into(),"coating_fast_time".into()),
                ("coating_slow".into(),"coating_slow_time".into()),
            ]},
        ];
        let make=|tau|bank_from_materials(cx,gate,[
            material_with_arms(&[],200_000.0),
            material_with_arms(&[("coating_fast","coating_fast_time",20e9,0.002),
                ("coating_slow","coating_slow_time",40e9,tau)],200_000.0),
        ],1,2);
        let bind=|tau|make(tau).with_material_layered_string_bending_relaxation(vec![
            LayeredStringBendingRelaxationSelection {string:0,layers:layers.clone()},
        ],None,16,gate).unwrap();
        let mut first=bind(0.004); let mut changed=bind(0.012);
        let image=first.source().string_bending_relaxation().unwrap();
        assert!(image.spectra()[0][0].terms.is_empty());
        assert_eq!(image.spectra()[0][1].terms,vec![(20e9,0.002),(40e9,0.004)]);
        assert_eq!(image.channels().len(),4); // Two coating arms, two bending planes.
        for (a,b) in first.system().components().iter().zip(changed.system().components()) {
            assert_eq!(a.modes(),b.modes());
        }
        let pressure=render(&mut first); let altered=render(&mut changed);
        assert!(pressure.iter().any(|p|p.abs()>1e-15));
        assert!(pressure.iter().zip(altered).any(|(a,b)|(a-b).abs()>1e-18));
        assert_ne!(first.system().components()[1].states(),changed.system().components()[1].states());
        assert!(first.system().coupled_frame().unwrap().material_dissipation_j>0.0);
        assert!(make(0.004).with_material_string_bending_relaxation(selection(),None,16,gate).is_err());
        assert!(make(0.004).with_material_layered_string_bending_relaxation(vec![
            LayeredStringBendingRelaxationSelection {string:0,layers:vec![layers[1].clone()]},
        ],None,16,gate).is_err());
    });
}


fn shear_layers() -> Vec<LayerRelaxationProperties> {
    vec![
        LayerRelaxationProperties {equilibrium_modulus:"G".into(),branches:vec![]},
        LayerRelaxationProperties {equilibrium_modulus:"G".into(),branches:vec![
            ("shear_fast".into(),"shear_fast_time".into()),
            ("shear_slow".into(),"shear_slow_time".into()),
        ]},
    ]
}
fn shear_source(delta:f64,tau:f64,band:f64,bind:bool) -> MaterialStringSource {
    let source=source_from_materials([
        material_with_arms(&[("bending_delta","bending_time",10e9,0.004)],200_000.0),
        material_with_arms(&[("bending_delta","bending_time",10e9,0.004),
            ("shear_fast","shear_fast_time",0.5*delta,0.002),
            ("shear_slow","shear_slow_time",delta,tau)],band),
    ],2);
    if bind {source.with_torsional_relaxation(&shear_layers()).unwrap()} else {source}
}

#[test]
fn torsional_material_memory_drives_actual_surface_friction_bridge_and_pressure() {
    use fs_couple::bowed_string::runtime::reciprocal::bank::torsion::{
        ModalRotationAttachment,StringBodyRotationConnection,
    };
    with_cx(|cx,gate| {
        let make=|tau|bank_with_sources(cx,gate,vec![
            shear_source(40e9,tau,200_000.0,true),
            original_source(20e9,0.004,200_000.0,false),
        ]).unwrap().with_rotational_body_connections(vec![StringBodyRotationConnection {
            string:0,station_fraction:0.75,
            body:ModalRotationAttachment {component:0,angular_shapes:vec![12.0]},
            stiffness_n_m_per_rad:0.002,damping_n_m_s_per_rad:0.0,rest_angle_rad:0.0,
        }],gate).unwrap();
        let mut first=make(0.004); let mut changed=make(0.012);
        let arms=first.system().material_relaxation_branches();
        assert_eq!(arms.len(),2); // Two coating arms, one actual twist mode.
        let j=core::f64::consts::PI*0.0003_f64.powi(4)/2.0;
        let coating_j=core::f64::consts::PI*(0.0003_f64.powi(4)-0.0002_f64.powi(4))/2.0;
        let wave_number=core::f64::consts::PI/first.cards()[0].length_m;
        for (arm,delta) in arms.iter().zip([20e9,40e9]) {
            let expected=delta*coating_j*wave_number.powi(2)/(7800.0*j);
            assert!((arm.stiffness-expected).abs()<1e-12*expected);
            assert_eq!(arm.projection[8],1.0); // After both two-mode bending planes.
            assert!(arm.projection.iter().enumerate().all(|(i,p)|i==8 || *p==0.0));
        }
        for (a,b) in first.system().components().iter().zip(changed.system().components()) {
            assert_eq!(a.modes(),b.modes());
        }
        let a=render(&mut first); let b=render(&mut changed);
        assert!(a.iter().any(|p|p.abs()>1e-15));
        assert!(a.iter().zip(b).any(|(p,q)|(p-q).abs()>1e-18));
        assert_ne!(first.system().components()[1].states(),changed.system().components()[1].states());
        assert!(first.system().components()[0].states()[4].velocity_m_sqrt_kg_per_s.abs()>0.0);
        assert!(first.system().material_memory().iter().any(|z|*z!=0.0));
        assert!(first.system().coupled_frame().unwrap().material_dissipation_j>0.0);
    });
}

#[test]
fn source_shear_history_keeps_bending_prefix_and_survives_gas_and_axial_rebuilds() {
    use fs_couple::air_path::cylinder::CylinderFlowResolution;
    use fs_couple::bowed_string::runtime::reciprocal::bank::nonlinear::longitudinal::LongitudinalStringResolution;
    use fs_couple::render::schedule::force::coupled::contact::friction::nonlinear::NonlinearModalConfig;
    with_cx(|cx,gate| {
        let make=||bank_with_sources(cx,gate,vec![
            original_source(40e9,0.004,200_000.0,true),
            shear_source(40e9,0.012,200_000.0,true),
        ]).unwrap();
        let dry=make();
        let loaded=make().with_cylinder_gas_loading(CylinderFlowResolution {
            minimum_wavenumber_radius:0.001,maximum_wavenumber_radius:30.0,intervals:4,
        },gate).unwrap();
        assert_eq!(loaded.system().material_memory().len(),6); // Four bending, two torsion.
        for (arm,mode) in loaded.system().material_relaxation_branches()[..4].iter().zip(0..4) {
            assert_eq!(arm.projection[2*mode],1.0);
        }
        for (a,b) in dry.system().material_relaxation_branches()[4..].iter()
            .zip(&loaded.system().material_relaxation_branches()[4..]) {
            assert_eq!(a.stiffness,b.stiffness); // Gas changes transverse, not solid rotary inertia.
            assert_eq!(b.projection[18],1.0); // Second member's actual twist coordinate.
        }
        let before=loaded.system().material_memory().to_vec();
        let mut explicit=selection(); explicit[0].string=1;
        let joined=loaded.with_material_string_bending_relaxation(explicit,None,16,gate).unwrap();
        assert_eq!(&joined.system().material_memory()[..6],before);
        assert_eq!(joined.system().material_memory().len(),14);
        let states=joined.system().material_memory().to_vec();
        let resolution=LongitudinalStringResolution {
            axial_modes:vec![1],maximum_slope:0.2,maximum_strain:0.01,
        };
        let mut axial=joined.with_longitudinal_waves(vec![resolution.clone(),resolution],
            NonlinearModalConfig {max_iterations:64,maximum_generalized_force:1e5,
                force_absolute_tolerance:1e-11,force_relative_tolerance:1e-10},gate).unwrap();
        assert_eq!(axial.system().material_memory(),states);
        assert!(axial.system().material_relaxation_branches().iter()
            .all(|branch|branch.projection.len()==2*axial.system().mode_count()));
        assert!(render(&mut axial).iter().all(|p|p.is_finite()));
    });
}

#[test]
fn zero_shear_arms_preserve_pressure_and_unsupported_history_consumers_refuse() {
    with_cx(|cx,gate| {
        let make=|bind|bank_with_sources(cx,gate,vec![shear_source(0.0,0.004,200_000.0,bind)]).unwrap();
        let mut ordinary=make(false); let mut zero=make(true);
        assert!(zero.system().material_memory().is_empty());
        assert_eq!(render(&mut ordinary),render(&mut zero));
        assert!(bank_with_sources(cx,gate,vec![shear_source(40e9,0.004,8_000.0,true)]).is_err());
        assert!(shear_source(40e9,0.004,200_000.0,true)
            .with_torsional_relaxation(&shear_layers()).is_err());
        let mut wrong=shear_layers(); wrong[0].equilibrium_modulus="young_modulus".into();
        assert!(shear_source(40e9,0.004,200_000.0,false).with_torsional_relaxation(&wrong).is_err());
        let initial=vec![fs_couple::modal_acoustic_time::ModalAcousticState {
            displacement_m_sqrt_kg:0.0,velocity_m_sqrt_kg_per_s:0.0,
        };5];
        assert!(shear_source(40e9,0.004,200_000.0,true)
            .into_nonlinear(cx,gas(),48_000,0.2,initial,gate).is_err());
    });
}

#[test]
fn fast_layered_bending_original_bending_and_torsion_share_the_exact_material_owner() {
    with_cx(|cx,gate| {
        let tau = 1e-7;
        let layered = bank(cx,gate,40e9,tau,200_000.0,1,1)
            .with_material_string_bending_relaxation(selection(),None,16,gate).unwrap();
        let original = bank_with_sources(cx,gate,
            vec![original_source(40e9,tau,200_000.0,true)]).unwrap();
        let torsion = bank_with_sources(cx,gate,
            vec![shear_source(40e9,tau,200_000.0,true)]).unwrap();
        for mut bank in [layered, original, torsion] {
            assert!(bank.system().sample_period_s() / tau > 200.0);
            assert!(bank.system().material_relaxation_branches().iter()
                .any(|branch| branch.relaxation_time_s == tau));
            let pressure = render(&mut bank);
            assert!(pressure.iter().all(|p|p.is_finite()));
            assert!(pressure.iter().any(|p|p.abs()>1e-15));
            assert!(bank.system().material_memory().iter().any(|z|*z!=0.0));
            assert!(bank.system().coupled_frame().unwrap().material_dissipation_j>0.0);
        }
    });
}

#[path = "material_string_relaxation/exclusive.rs"]
mod exclusive;
