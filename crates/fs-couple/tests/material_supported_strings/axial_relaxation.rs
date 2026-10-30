//! Data-only layer spectra reach the original moving-support/impact owner.
use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::{polarization::TransverseDirections,
    nonlinear::{longitudinal::{LongitudinalStringResolution,relaxation::LongitudinalRelaxationSelection},
        supported::longitudinal::MaterialLongitudinalSupport}};
use fs_couple::modal_acoustic_time::ModalAcousticState;
use fs_couple::render::plate::impact::{relaxation::InitialMemory,
    supported::{StringEndMotion,longitudinal::StringSupportGeometry}};
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet,InterpolationPolicy,MaterialCard,MaterialStateId,PropertyClaim,PropertyKey,
    PropertyValue,Provenance,QueryPoint,UncertaintyModel};
use fs_material::{section::relaxation::LayerRelaxationProperties,
    state_point::{MaterialPropertySelection,ResolvedMaterialStatePoint,ScalarAdmissibility,
        ScalarPropertyRequirement,resolve_material_state_point}};
use fs_qty::{Density,Dims,DynViscosity,Pressure,Time};

const RATE:u32=6_144_000;
fn material(young:f64,delta:f64,tau:Option<f64>,high:f64,loss_high:f64)->ResolvedMaterialStatePoint {
    let mut values=vec![("density",Density::DIMS,7800.0),("young_modulus",Pressure::DIMS,young),
        ("poisson_ratio",Dims::NONE,0.3),(ETA,DynViscosity::DIMS,0.0),("delta",Pressure::DIMS,delta)];
    if let Some(tau)=tau {values.push(("tau",Time::DIMS,tau));}
    let mut claims=ClaimSet::new();let mut requirements=Vec::new();
    for (name,dims,value) in values {
        claims.insert_claim(PropertyClaim {key:PropertyKey::new(name,dims),value:PropertyValue::Scalar {value,dims},
            validity:ValidityDomain::unconstrained().with("T",280.0,320.0).with("omega",0.0,
                if name=="tau" {high} else if name==ETA {loss_high} else {1e8}),
            uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::ConstantWithinValidity,
            observations:vec![],provenance:Provenance {source:"synthetic supported hereditary specimen".into(),
                license:"CC0-1.0".into(),artifact:None},
        }).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(name,dims,ScalarAdmissibility::Finite).unwrap());
    }
    let card=MaterialCard::assemble(MaterialStateId {chemistry:"synthetic".into(),phase:"solid".into(),
        process:"hereditary-test".into(),revision:0},claims,vec![]).unwrap();
    resolve_material_state_point(&card,&QueryPoint::new().with("T",293.15).unwrap().with("omega",1000.0).unwrap(),
        &requirements,MaterialPropertySelection::SingleClaimOnly).unwrap()
}
fn source(cx:&Cx<'_>,gate:&CancelGate,delta:f64,tau:Option<f64>,high:f64,loss_high:f64)->NonlinearMaterialString {
    let section=ConcentricElasticSection::try_new(0.0,vec![
        ConcentricMaterialLayer {outer_radius_m:0.0002,material:material(200e9,delta,tau,high,loss_high)},
        ConcentricMaterialLayer {outer_radius_m:0.0003,material:material(50e9,0.5*delta,tau,high,loss_high)},
    ]).unwrap().with_kelvin_voigt_bending(ETA).unwrap();
    MaterialStringSource::Concentric {template:template(1),section,prestress:StringPrestress::FixedTension(20.0)}
        .with_transverse_polarizations(TransverseDirections::try_new([1.0,0.0],[1.0,0.0]).unwrap()).unwrap()
        .into_nonlinear(cx,gas(),RATE,0.2,vec![ModalAcousticState {displacement_m_sqrt_kg:1e-6,
            velocity_m_sqrt_kg_per_s:0.0},ModalAcousticState::default()],gate).unwrap()
}
fn image(cx:&Cx<'_>,gate:&CancelGate,delta:f64,tau:Option<f64>,high:f64)->SupportedMaterialStrings {
    image_with_loss_band(cx,gate,delta,tau,high,1e8)
}
fn image_with_loss_band(cx:&Cx<'_>,gate:&CancelGate,delta:f64,tau:Option<f64>,high:f64,loss_high:f64)->SupportedMaterialStrings {
    let geometry=StringSupportGeometry::try_new([[0.0;3],[0.5,0.0,0.0]],
        [0.0,1.0,0.0],[0.6,0.8,0.0]).unwrap();
    SupportedMaterialStrings::in_stationary_gas_with_longitudinal_waves(
        vec![source(cx,gate,delta,tau,high,loss_high)],TranslatingSupport {initial_position_m:0.0,
            initial_velocity_m_s:0.0,..rail()},vec![MaterialLongitudinalSupport {
                moving_end:StringEndMotion::Second,geometry,resolution:LongitudinalStringResolution {
                    axial_modes:vec![1,2],maximum_slope:0.2,maximum_strain:0.01},
                initial_axial:vec![ModalAcousticState::default();2],
            }],300.0,gate).unwrap()
}
fn selection()->Vec<LongitudinalRelaxationSelection> {
    vec![LongitudinalRelaxationSelection {string:0,layers:vec![LayerRelaxationProperties {
        equilibrium_modulus:"young_modulus".into(),branches:vec![("delta".into(),"tau".into())],
    };2]}]
}
#[test]
fn g3_source_modulus_changes_support_dynamics_without_changing_geometry_or_equilibrium() {
    with_cx(|cx,gate| {
        let mut states=Vec::new();let mut ea=Vec::new();
        for delta in [10e9,30e9] {
            let original=image(cx,gate,delta,Some(0.001),1e8);ea.push(original.sources()[0].axial_rigidity_n());
            let bound=original.supported().clock_bounds().to_vec();
            let loaded=original.with_material_longitudinal_relaxation(selection(),InitialMemory::Relaxed,6,gate).unwrap();
            assert_eq!(loaded.axial_material_relaxation().unwrap().sections().len(),1);
            assert!(loaded.supported().clock_bounds().iter().zip(bound).all(|(a,b)|*a>b));
            // Copying the source body into a generic impact assembly preserves
            // all six material coordinates and the source physical clock.
            let mut direct=ImpactSystem::new(loaded.bodies().to_vec(),vec![],vec![],vec![],config_at(RATE)).unwrap();
            assert_eq!(direct.relaxation_branch_count(),6);
            let force=loaded.support_shape().iter().map(|b|b*0.01).collect::<Vec<_>>();
            let mut run=loaded.prepare(config_at(RATE),gate).unwrap();let before=run.system().stored_energy_j();
            let(mut work,mut loss)=(0.0,0.0);
            for _ in 0..32 {
                let frame=run.step(&force,gate).unwrap();direct.step(&force,gate).unwrap();
                work+=frame.supplied_work_j;loss+=frame.dissipated_energy_j;
                assert!((frame.stored_energy_j+loss-before-work).abs()<1e-9);
            }
            for (a,b) in run.system().state().iter().zip(direct.state()) {assert!((a-b).abs()<1e-7);}
            assert!(loss>0.0);states.push(run.system().state().to_vec());
        }
        assert_eq!(ea[0],ea[1]);
        // Mechanical support displacement must change, not merely the energy
        // normalization of the appended material-history coordinates.
        assert!((states[0][8]-states[1][8]).abs()>1e-18);
    });
}
#[test]
fn g0_source_missing_data_band_clock_and_complete_history_budget_are_not_fallbacks() {
    with_cx(|cx,gate| {
        for (tau,band,budget) in [(None,1e8,6),(Some(0.001),4000.0,6),
            (Some(1e-12),1e8,6),(Some(0.001),1e8,5)] {
            assert!(image(cx,gate,10e9,tau,band).with_material_longitudinal_relaxation(
                selection(),InitialMemory::Relaxed,budget,gate).is_err());
        }
        assert!(image_with_loss_band(cx,gate,10e9,Some(0.001),1e8,4000.0)
            .with_material_longitudinal_relaxation(selection(),InitialMemory::Relaxed,6,gate).is_err());
        let loaded=image(cx,gate,10e9,Some(0.001),1e8).with_material_longitudinal_relaxation(
            selection(),InitialMemory::Relaxed,6,gate).unwrap();
        let mut wrong=config_at(RATE);wrong.dt_s*=0.5;
        assert!(ImpactSystem::new(loaded.bodies().to_vec(),vec![],vec![],vec![],wrong).is_err());
        assert!(loaded.with_material_longitudinal_relaxation(selection(),InitialMemory::Relaxed,6,gate).is_err());
    });
}
