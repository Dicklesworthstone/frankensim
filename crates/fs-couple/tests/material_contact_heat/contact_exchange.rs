//! Original mechanical footprint drives a source-backed finite thermal contact.
use super::*;
use fs_conduction::interface::{InterfaceResistance,AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
    AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY};
use fs_conduction::interface::lumped_pair::{advance_contact_pair,ContactPairReport};
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet,InterfaceSystemCard,InterpolationPolicy,PropertyClaim,PropertyKey,
    PropertyValue,Provenance,SelectionPolicy,SurfaceSpec,SystemContext,UncertaintyModel};
use fs_couple::render::plate::impact::material::heating::transport::ImpactThermalTrial;

fn source_resistance(value:f64)->InterfaceResistance {
    let a=solid("A",1000.0,1e7);let b=solid("B",2000.0,2e7);
    let dims=AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS;let mut claims=ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key:PropertyKey::new(AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,dims),
        value:PropertyValue::Scalar {value,dims},validity:ValidityDomain::unconstrained(),
        uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::ConstantWithinValidity,
        observations:vec![],provenance:Provenance {source:"synthetic thermal interface".into(),
            license:"CC0-1.0".into(),artifact:None},
    }).unwrap();
    let card=InterfaceSystemCard::assemble(
        SurfaceSpec {material:a.id().clone(),texture_frame:"surface-a".into()},
        SurfaceSpec {material:b.id().clone(),texture_frame:"surface-b".into()},
        SystemContext {medium:"dry".into(),third_body:None,environment:"test".into(),history:"initial".into()},
        claims,vec![]).unwrap();
    InterfaceResistance::from_card("touch",&card,&query(),SelectionPolicy::SingleClaimOnly).unwrap()
}
fn proposal(trial:ImpactThermalTrial<'_>,cx:&Cx<'_>,source:&InterfaceResistance)
    ->Result<CoupledThermalStep<ContactPairReport>,ImpactError> {
    let area=trial.contact_patch_area_m2("touch")?;
    let transport=LumpedThermalTransport::try_declared(100.0,0.0).unwrap();
    let inputs=trial.reservoirs();assert_eq!(inputs.len(),2);
    let transfers=advance_contact_pair(cx,[inputs[0].thermal,inputs[1].thermal],
        [&transport,&transport],source,area,1e-13).map_err(|e|ImpactError::Owner(e.to_string()))?;
    Ok(CoupledThermalStep {transfers:Vec::from(transfers),boundary_heat_j:0.0})
}

#[test]
fn actual_collision_opens_and_closes_source_thermal_exchange_without_changing_mechanics() {
    let c=curve(0.02,1.0);let low=source_resistance(1e-4);let high=source_resistance(4e-4);
    context(|gate,cx| {
        let make=||HeatedImpactSystem::new(collision(Some(0.4),0.0,false).prepare_analytic().unwrap(),
            warm_stores(&c),vec![route("touch")],1e-10).unwrap();
        let mut a=make();let mut b=make();let initial=a.mechanics().stored_energy_j();
        let (mut transferred,mut slow)=(0.0,0.0);let (mut dry,mut active)=(0,0);
        for sample in 1..=1800 {
            let fa=a.step_with_coupled_transport(&[0.0;2],gate,|trial|proposal(trial,cx,&low)).unwrap();
            let fb=b.step_with_coupled_transport(&[0.0;2],gate,|trial|proposal(trial,cx,&high)).unwrap();
            assert_eq!(a.mechanics().state(),b.mechanics().state());
            let physical=a.mechanics().material_contact_observation(0).unwrap();
            let report=&fa.transfers[0].report;
            near(report.area_m2,std::f64::consts::PI*physical.patch_radius_m.powi(2),1e-18);
            assert_eq!(report.interface,low);assert_eq!(fb.transfers[0].report.interface,high);
            assert_eq!(fa.boundary_heat_j,0.0);assert_eq!(fb.boundary_heat_j,0.0);
            if report.area_m2==0.0 {
                dry+=1;assert!(report.contact.is_none() && report.resistance.is_none());
                assert_eq!(fa.transfers[0].external_heat_j,0.0);
                assert_eq!(fa.transfers[1].external_heat_j,0.0);
                if sample>1500 {assert_eq!(fa.impact.deposited_heat_j,0.0);}
            } else {
                active+=1;assert!(report.contact.is_some());
                near(report.resistance.as_ref().unwrap().budget().value_k_per_w,
                    low.value_m2_k_per_w()/report.area_m2,1e-10);
                transferred+=fa.transfers[1].external_heat_j;slow+=fb.transfers[1].external_heat_j;
            }
            assert!(fa.impact.combined_energy_residual_j.abs()<1e-9);
        }
        assert!(dry>0 && active>0);assert!(transferred>slow+1e-7 && slow>0.0);
        assert!(a.reservoirs()[0].state().temperature_k()<b.reservoirs()[0].state().temperature_k());
        near(a.mechanics().stored_energy_j()+a.reservoirs().iter().map(|s|s.storage_change_j()).sum::<f64>(),initial,1e-9);
        let before=a.reservoirs().to_vec();
        a.step_with_coupled_transport(&[0.0;2],gate,|trial|proposal(trial,cx,&low)).unwrap();
        assert_eq!(a.reservoirs(),before);
    });
}
fn input(c:&EquilibriumEnthalpyPhaseCurve,h:f64,heat:f64,dt:f64)->UniformEnthalpyStepInput<'_> {
    UniformEnthalpyStepInput {curve:c,initial:c.state_at_specific_enthalpy(h).unwrap(),
        mass_kg:0.001,volume_m3:1e-6,surface_area_m2:0.0006,internal_heat_j:heat,duration_s:dt}
}

#[test]
fn simultaneous_exchange_admits_a_valid_joint_endpoint_without_an_invalid_heat_only_stage() {
    let c=curve(0.02,1.0);let source=source_resistance(1e-4);
    let transport=LumpedThermalTransport::try_declared(100.0,0.0).unwrap();
    context(|_,cx| {
        let inputs=[input(&c,0.9,0.0002,0.01),input(&c,0.1,0.0,0.01)];
        // Heating A alone would leave its chart, but the simultaneous finite
        // exchange has both endpoints strictly inside the original charts.
        assert!(c.state_at_specific_enthalpy(1.1).is_err());
        let next=advance_contact_pair(cx,inputs,[&transport,&transport],&source,1e-6,1e-12).unwrap();
        assert!(next.iter().all(|p|p.state.specific_enthalpy_j_kg()<1.0));
        assert_eq!(next[0].external_heat_j,-next[1].external_heat_j);
        assert!(advance_contact_pair(cx,inputs,[&transport,&transport],&source,0.0,1e-12).is_err());
        let mut opened=[inputs[0],inputs[1]];
        for i in 0..2 {opened[i].initial=next[i].state;opened[i].internal_heat_j=0.0;}
        let dry=advance_contact_pair(cx,opened,[&transport,&transport],&source,0.0,1e-12).unwrap();
        for i in 0..2 {assert_eq!(dry[i].state,next[i].state);assert_eq!(dry[i].external_heat_j,0.0);}
        let resumed=advance_contact_pair(cx,opened,[&transport,&transport],&source,1e-6,1e-12).unwrap();
        assert!(resumed[0].state.temperature_k()<next[0].state.temperature_k());
        assert!(resumed[1].state.temperature_k()>next[1].state.temperature_k());
    });
}

#[test]
fn source_latent_heat_changes_fraction_on_the_same_finite_contact_exchange() {
    let latent=|width|EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([9;32]),vec![
        EnthalpyPhaseKnot {specific_enthalpy_j_kg:0.0,temperature_k:310.0,liquid_mass_fraction:0.0,bulk_density_kg_m3:1000.0},
        EnthalpyPhaseKnot {specific_enthalpy_j_kg:width,temperature_k:310.0,liquid_mass_fraction:1.0,bulk_density_kg_m3:900.0},
    ]).unwrap();
    let small=latent(0.4);let large=latent(0.8);let hot=curve(0.02,1.0);
    let source=source_resistance(1e-4);let transport=LumpedThermalTransport::try_declared(100.0,0.0).unwrap();
    context(|_,cx| {
        let advance=|c|advance_contact_pair(cx,[input(c,0.0,0.0,0.001),input(&hot,0.5,0.0,0.001)],
            [&transport,&transport],&source,1e-6,1e-13).unwrap();
        let a=advance(&small);let b=advance(&large);
        near(a[0].external_heat_j,b[0].external_heat_j,2e-13);assert!(a[0].external_heat_j>0.0);
        assert_eq!(a[0].state.temperature_k(),310.0);assert_eq!(b[0].state.temperature_k(),310.0);
        near(a[0].state.liquid_mass_fraction(),2.0*b[0].state.liquid_mass_fraction(),1e-8);
    });
}

#[test]
fn open_and_closed_pair_refuse_foreign_states_invalid_clocks_areas_and_cancelled_work() {
    let c=curve(0.02,1.0);let foreign=curve(0.04,1.0);let source=source_resistance(1e-4);
    let transport=LumpedThermalTransport::try_declared(100.0,0.0).unwrap();
    context(|gate,cx| {
        let inputs=[input(&c,0.4,0.0,0.001),input(&c,0.1,0.0,0.001)];
        for area in [0.0,1e-6] {
            let mut bad=inputs;bad[1].duration_s*=2.0;
            assert!(advance_contact_pair(cx,bad,[&transport,&transport],&source,area,1e-12).is_err());
            bad=inputs;bad[1].initial=foreign.state_at_specific_enthalpy(0.1).unwrap();
            assert!(advance_contact_pair(cx,bad,[&transport,&transport],&source,area,1e-12).is_err());
            bad=inputs;bad[0].internal_heat_j=f64::NAN;
            assert!(advance_contact_pair(cx,bad,[&transport,&transport],&source,area,1e-12).is_err());
        }
        for area in [-1.0,0.001,f64::NAN,f64::INFINITY] {
            assert!(advance_contact_pair(cx,inputs,[&transport,&transport],&source,area,1e-12).is_err());
        }
        gate.request();
        for area in [0.0,1e-6] {
            assert!(matches!(advance_contact_pair(cx,inputs,[&transport,&transport],&source,area,1e-12),
                Err(fs_conduction::ConductionError::Cancelled {..})));
        }
    });
}
