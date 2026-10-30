//! The SAME material constitutive receipts must reach axial waves, the actual
//! moving termination and unplayed strings, while retaining gas and torsion.
use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::{
    polarization::TransverseDirections,
    nonlinear::{longitudinal::{LongitudinalMaterialLoss,LongitudinalStringResolution},
        supported::longitudinal::MaterialLongitudinalSupport},
};
use fs_couple::modal_acoustic_time::ModalAcousticState;
use fs_couple::render::plate::impact::supported::{StringEndMotion,longitudinal::StringSupportGeometry};
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet,InterpolationPolicy,MaterialCard,MaterialStateId,PropertyClaim,
    PropertyKey,PropertyValue,Provenance,QueryPoint,UncertaintyModel};
use fs_material::{section::axial::YOUNG_VISCOSITY_PROPERTY,
    state_point::{MaterialPropertySelection,ResolvedMaterialStatePoint,ScalarAdmissibility,
        ScalarPropertyRequirement,resolve_material_state_point}};
use fs_qty::{Density,Dims,DynViscosity,Pressure};

const RATE:u32=1_536_000;

fn material(rho:f64,young:f64,eta:Option<f64>,band:(f64,f64))->ResolvedMaterialStatePoint {
    let mut claims=ClaimSet::new();let mut requirements=Vec::new();
    let mut values=vec![("density",Density::DIMS,rho),("young_modulus",Pressure::DIMS,young),
        ("poisson_ratio",Dims::NONE,0.3),(ETA,DynViscosity::DIMS,10.0),
        ("G",Pressure::DIMS,60e9),("eta_s",DynViscosity::DIMS,1000.0)];
    if let Some(eta)=eta {values.push((YOUNG_VISCOSITY_PROPERTY,DynViscosity::DIMS,eta));}
    for (name,dims,value) in values {
        let (lo,hi)=if name==YOUNG_VISCOSITY_PROPERTY {band} else {(0.0,200_000.0)};
        claims.insert_claim(PropertyClaim {key:PropertyKey::new(name,dims),value:PropertyValue::Scalar {value,dims},
            validity:ValidityDomain::unconstrained().with("T",280.0,320.0).with("omega",lo,hi),
            uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::ConstantWithinValidity,
            observations:vec![],provenance:Provenance {source:"synthetic supported axial viscosity".into(),
                license:"CC0-1.0".into(),artifact:None},
        }).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(name,dims,ScalarAdmissibility::Finite).unwrap());
    }
    let card=MaterialCard::assemble(MaterialStateId {chemistry:"synthetic".into(),phase:"solid".into(),
        process:"axial-test".into(),revision:0},claims,vec![]).unwrap();
    resolve_material_state_point(&card,&QueryPoint::new().with("T",293.15).unwrap().with("omega",1000.0).unwrap(),
        &requirements,MaterialPropertySelection::SingleClaimOnly).unwrap()
}
fn sources(cx:&Cx<'_>,gate:&CancelGate,eta:Option<f64>,band:(f64,f64))->Vec<NonlinearMaterialString> {
    let section=ConcentricElasticSection::try_new(0.0,vec![
        ConcentricMaterialLayer {outer_radius_m:0.0002,material:material(7800.0,200e9,Some(0.0),band)},
        ConcentricMaterialLayer {outer_radius_m:0.0003,material:material(15600.0,50e9,eta,band)},
    ]).unwrap().with_kelvin_voigt_bending(ETA).unwrap();
    let directions=TransverseDirections::try_new([0.6,0.8],[1.0,0.0]).unwrap();
    let mut first_initial=initial(3,0.0002);first_initial[2].velocity_m_sqrt_kg_per_s=0.001;
    let first=MaterialStringSource::Concentric {template:template(1),section,prestress:StringPrestress::FixedTension(20.0)}
        .with_torsion("G","eta_s",1).unwrap().with_transverse_polarizations(directions).unwrap()
        .into_nonlinear(cx,gas(),RATE,0.2,first_initial,gate).unwrap();
    let second=MaterialStringSource::Homogeneous(fs_couple::string_specimen::with_uniform_circular_material_state(
        template(1),0.0003,&material(7800.0,200e9,Some(0.0),(0.0,200_000.0)))
        .unwrap().with_kelvin_voigt_bending_loss().unwrap())
        .with_transverse_polarizations(directions).unwrap()
        .into_nonlinear(cx,gas(),RATE,0.2,initial(2,0.0),gate).unwrap();
    vec![first,second]
}
fn attachments()->Vec<MaterialLongitudinalSupport> {
    let geometry=StringSupportGeometry::try_new([[0.0;3],[0.5,0.0,0.0]],
        [0.0,1.0,0.0],[0.6,0.8,0.0]).unwrap();
    [StringEndMotion::Second,StringEndMotion::First].into_iter().enumerate()
        .map(|(i,moving_end)|MaterialLongitudinalSupport {moving_end,geometry,
            resolution:LongitudinalStringResolution {axial_modes:vec![1,2],maximum_slope:0.2,maximum_strain:0.01},
            initial_axial:vec![ModalAcousticState {displacement_m_sqrt_kg:if i==0 {1e-8} else {0.0},
                velocity_m_sqrt_kg_per_s:0.0},ModalAcousticState::default()],
        }).collect()
}
fn support()->TranslatingSupport {
    TranslatingSupport {initial_position_m:0.0,initial_velocity_m_s:0.0,..rail()}
}
fn selected()->LongitudinalMaterialLoss {
    LongitudinalMaterialLoss::KelvinVoigt {viscosity_property:YOUNG_VISCOSITY_PROPERTY.into()}
}
fn build(cx:&Cx<'_>,gate:&CancelGate,eta:Option<f64>,band:(f64,f64),loss:LongitudinalMaterialLoss)
    ->Result<SupportedMaterialStrings,fs_couple::render::plate::impact::ImpactError>
{
    SupportedMaterialStrings::in_stationary_gas_with_longitudinal_waves_and_loss(
        sources(cx,gate,eta,band),support(),attachments(),loss,300.0,gate)
}

#[test]
fn g3_layer_material_loss_reaches_axial_waves_support_and_unplayed_string() {
    with_cx(|cx,gate| {
        let mut histories=Vec::new();let mut initial_energies=Vec::new();
        for eta in [1e5,2e6] {
            let sources=sources(cx,gate,Some(eta),(0.0,200_000.0));
            let old=sources.clone();
            let elastic=SupportedMaterialStrings::in_stationary_gas_with_longitudinal_waves(
                sources.clone(),support(),attachments(),300.0,gate).unwrap();
            let image=SupportedMaterialStrings::in_stationary_gas_with_longitudinal_waves_and_loss(
                sources,support(),attachments(),selected(),300.0,gate).unwrap();
            let laws=image.axial_material_loss().unwrap();
            let expected=eta*core::f64::consts::PI*(0.0003_f64.powi(2)-0.0002_f64.powi(2));
            assert!((laws[0].viscous_rigidity_n_s()-expected).abs()<1e-14*expected);
            assert_eq!(laws[0].layer_viscous_rigidity_n_s()[0],0.0);
            assert_eq!(laws[1].viscous_rigidity_n_s(),0.0);
            assert_eq!(image.supported().axial_viscosity_n_s().unwrap(),
                &[laws[0].viscous_rigidity_n_s(),0.0]);
            assert_eq!(image.support().damping_n_s_m,support().damping_n_s_m);
            assert_eq!(image.supported().total_mass_kg(),elastic.supported().total_mass_kg());
            assert_eq!(image.supported().absolute_drag(),elastic.supported().absolute_drag());
            assert_eq!(image.torsion_ranges(),elastic.torsion_ranges());
            for (actual,reference) in image.bodies().iter().zip(elastic.bodies()) {
                assert_eq!(actual.initial,reference.initial);
                assert_eq!(actual.damping_per_s,reference.damping_per_s);
            }
            assert_eq!(image.bodies()[1].initial,old[0].bodies()[1].initial);
            let mu=laws[0].elastic().line_density_kg_m();let ea=laws[0].elastic().axial_stiffness_n();
            let w=image.supported().axial_angular_frequencies(0).unwrap()[0];
            let k=core::f64::consts::PI/old[0].card().length_m;
            assert!((w*w/(k*k)-ea/mu).abs()<1e-13*ea/mu);
            let support_index=image.supported().support_coordinate();
            let mut run=image.prepare(config_at(RATE),gate).unwrap();
            let before=run.system().stored_energy_j();initial_energies.push(before);
            let (mut loss,mut bridge,mut neighbor)=(0.0,0.0_f64,0.0_f64);
            for _ in 0..512 {
                let frame=run.step(&[0.0;10],gate).unwrap();
                assert!(frame.dissipated_energy_j>=0.0);loss+=frame.dissipated_energy_j;
                assert!((frame.stored_energy_j+loss-before).abs()<1e-9);
                let state=run.system().state();
                bridge=bridge.max(state[2*support_index].abs());neighbor=neighbor.max(state[12].abs());
            }
            assert!(loss>0.0 && bridge>1e-14 && neighbor>1e-16);
            histories.push(run.system().state().to_vec());
        }
        assert_eq!(initial_energies[0],initial_energies[1],"viscosity does not alter mass or storage");
        assert_ne!(histories[0][16],histories[1][16],"source viscosity must change the physical support");
        assert_ne!(histories[0][12],histories[1][12],"source viscosity must reach the unplayed string");
    });
}

#[test]
fn g0_explicit_zero_axial_viscosity_is_bitwise_elastic_without_dropping_receipts() {
    with_cx(|cx,gate| {
        let damped=build(cx,gate,Some(0.0),(0.0,200_000.0),selected()).unwrap();
        let elastic=build(cx,gate,Some(0.0),(0.0,200_000.0),LongitudinalMaterialLoss::Elastic).unwrap();
        assert!(damped.axial_material_loss().is_some());assert!(elastic.axial_material_loss().is_none());
        let mut damped=damped.prepare(config_at(RATE),gate).unwrap();
        let mut elastic=elastic.prepare(config_at(RATE),gate).unwrap();
        for _ in 0..64 {
            let a=damped.step(&[0.0;10],gate).unwrap();let b=elastic.step(&[0.0;10],gate).unwrap();
            assert_eq!(a.stored_energy_j,b.stored_energy_j);assert_eq!(a.dissipated_energy_j,b.dissipated_energy_j);
            assert_eq!(damped.system().state(),elastic.system().state());
        }
    });
}

#[test]
fn g0_missing_layer_wrong_loss_or_uncovered_support_frequency_refuses_before_playback() {
    with_cx(|cx,gate| {
        for (eta,band) in [(None,(0.0,200_000.0)),(Some(-1.0),(0.0,200_000.0)),
            (Some(1e5),(1.0,200_000.0)),(Some(1e5),(0.0,10_000.0))] {
            assert!(build(cx,gate,eta,band,selected()).is_err());
        }
        assert!(build(cx,gate,Some(1e5),(0.0,200_000.0),LongitudinalMaterialLoss::KelvinVoigt {
            viscosity_property:"density".into(),
        }).is_err());
        let elastic=build(cx,gate,Some(1e5),(0.0,200_000.0),LongitudinalMaterialLoss::Elastic).unwrap();
        let fixed=elastic.supported().axial_angular_frequencies(0).unwrap().iter().copied().fold(0.0_f64,f64::max);
        let coupled=elastic.supported().clock_bounds()[0];assert!(coupled>fixed);
        // The fixed-end poles fit, but their moving support does not. Checking
        // only omega_n as the fixed-end bank does would accept this source.
        assert!(build(cx,gate,Some(1e5),(0.0,f64::midpoint(fixed,coupled)),selected()).is_err());
        let source=sources(cx,gate,Some(1e5),(0.0,200_000.0));
        let stop=CancelGate::new_clock_free();stop.request();
        assert!(SupportedMaterialStrings::in_stationary_gas_with_longitudinal_waves_and_loss(
            source,support(),attachments(),selected(),300.0,&stop).is_err());
    });
}

#[test]
fn g3_full_material_strain_loss_drives_the_same_support_and_closes_mechanical_work() {
    with_cx(|cx,gate| {
        let full=||LongitudinalMaterialLoss::FiniteStrainKelvinVoigt {
            viscosity_property:YOUNG_VISCOSITY_PROPERTY.into(),
        };
        let mut trajectories=Vec::new();
        for eta in [1e5,2e6] {
            let image=build(cx,gate,Some(eta),(0.0,200_000.0),full()).unwrap();
            assert!(image.supported().has_finite_strain_axial_viscosity());
            let small=build(cx,gate,Some(eta),(0.0,200_000.0),selected()).unwrap();
            assert_eq!(image.bodies()[0].initial,small.bodies()[0].initial);
            assert_eq!(image.supported().absolute_drag(),small.supported().absolute_drag());
            assert_eq!(image.torsion_ranges(),small.torsion_ranges());
            let index=image.supported().support_coordinate();
            let row=image.support_shape();
            let mut run=image.prepare(config_at(RATE),gate).unwrap();
            let initial=run.system().stored_energy_j();
            let (mut loss,mut work)=(0.0,0.0);
            for tick in 0..128 {
                let force=if tick<16 {0.01} else {0.0};
                let before=run.system().state().to_vec();
                let frame=run.step(&row.iter().map(|b|force*b).collect::<Vec<_>>(),gate).unwrap();
                let state=run.system().state();
                let actual_work=force*row.iter().enumerate().map(|(i,b)|b*(state[2*i]-before[2*i])).sum::<f64>();
                assert!((frame.supplied_work_j-actual_work).abs()<1e-12);
                assert!(frame.dissipated_energy_j>=0.0);
                loss+=frame.dissipated_energy_j;work+=frame.supplied_work_j;
                assert!((frame.stored_energy_j+loss-initial-work).abs()<1e-9);
            }
            assert!(loss>0.0);
            trajectories.push(run.system().state().to_vec());
            assert_ne!(trajectories.last().unwrap()[2*index],0.0);
            // A failed/cancelled tick must leave the original mechanical clock
            // and every coupled state untouched, including torsion and support.
            let before=run.system().state().to_vec();let samples=run.system().samples();
            let stop=CancelGate::new_clock_free();stop.request();
            assert!(run.step(&vec![0.0;row.len()],&stop).is_err());
            assert_eq!(before,run.system().state());assert_eq!(samples,run.system().samples());
        }
        assert_ne!(trajectories[0],trajectories[1],"only the source viscosity changed");
        let image=build(cx,gate,Some(0.0),(0.0,200_000.0),full()).unwrap();
        assert!(image.axial_material_loss().is_some());
        assert!(!image.supported().has_finite_strain_axial_viscosity());
        let mut zero=image.prepare(config_at(RATE),gate).unwrap();
        let mut elastic=build(cx,gate,Some(0.0),(0.0,200_000.0),LongitudinalMaterialLoss::Elastic)
            .unwrap().prepare(config_at(RATE),gate).unwrap();
        for _ in 0..16 {
            let a=zero.step(&[0.0;10],gate).unwrap();let b=elastic.step(&[0.0;10],gate).unwrap();
            assert_eq!(a.dissipated_energy_j,b.dissipated_energy_j);
            assert_eq!(zero.system().state(),elastic.system().state());
        }
        for (eta,band) in [(None,(0.0,200_000.0)),(Some(1e5),(1.0,200_000.0)),
            (Some(1e5),(0.0,10_000.0))] {
            assert!(build(cx,gate,eta,band,full()).is_err());
        }
    });
}
