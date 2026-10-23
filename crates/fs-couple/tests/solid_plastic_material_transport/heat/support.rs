//! G0/G1/G3: solved support pressure and bulk heat transport retain their owners.
#[path = "support/electrical.rs"]
mod electrical;
use super::*;
use fs_couple::solid_enthalpy::plastic::contact::{ContactError, ContactSolve, ThermalSupportContact};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use std::convert::Infallible;

fn support_controls() -> ContactSolve {
    ContactSolve { maximum_cells:2, maximum_mechanical_trials:1000,
        maximum_newton_iterations:32, difference_step_j:1e-7,
        cell_energy_tolerance_j:1e-10, energy_tolerance_j:1e-8, yield_tolerance_pa:1e-3,
        newton:NewtonKrylovConfig { absolute_tolerance:1e-12, relative_tolerance:1e-12,
            linear_restart:2, ..NewtonKrylovConfig::default() } }
}
fn displacements(compression:f64, shear:f64) -> Vec<(usize,f64)> {
    (0..12).filter(|&d|d!=9).map(|d|(d,-compression*NODES[d/3][d%3]
        + if d%3==0 {2.0*shear*NODES[d/3][1]} else {0.0})).collect()
}

fn with_support<R>(extra:Option<PropertyClaim>, maximum_pressure:f64,
    run:impl FnOnce(&Cx<'_>,&CancelGate,&Fixture,&DepositionSimulation<'_>,
        &ThermoplasticSolid<'_>,&ThermalSupportContact<'_>,&MaterialHeatNetwork<'_>)->R)->R {
    with_heat_fixture(fixture_with_property(0.7,extra),0.1,|cx,gate,f,t,_,network| {
        // The lower thermal cell is a finite support, NOT another mechanical tet.
        let s=ThermoplasticSolid::try_new(t,&NODES[..4],&TETS[..1],
            f.sources().into_iter().take(1).collect(),TetAssemblyBudget::standard(),cx).unwrap();
        let mut claim=scalar(CONTACT_R,CONTACT_DIMS,0.5);
        claim.value=PropertyValue::Curve {abscissa:"pressure".into(),abscissa_dims:Pressure::DIMS,
            knots:vec![(0.0,0.5),(maximum_pressure,0.1)],dims:CONTACT_DIMS};
        claim.validity=ValidityDomain::unconstrained().with("pressure",0.0,maximum_pressure)
            .with("Ta",280.0,320.0).with("Tb",280.0,320.0);
        claim.interpolation=InterpolationPolicy::LinearInside;
        let mut claims=ClaimSet::new();let pin=claims.insert_claim(claim).unwrap();
        let card=InterfaceSystemCard::assemble(
            SurfaceSpec {material:f.cards[0].id().clone(),texture_frame:"specimen".into()},
            SurfaceSpec {material:f.cards[1].id().clone(),texture_frame:"support".into()},
            SystemContext {medium:"synthetic dry contact".into(),third_body:None,
                environment:"test".into(),history:"fixed footprint".into()},claims,vec![]).unwrap();
        // Never the mechanical pressure: the coupling must replace this value.
        let point=QueryPoint::new().with("Ta",300.0).unwrap().with("Tb",300.0).unwrap()
            .with("pressure",7e6).unwrap();
        let contact=ThermalSupportContact {link:91,face_nodes:[2,0,1],
            point:ThermalInterfacePoint {interface:&card,materials:[&f.cards[0],&f.cards[1]],
                context:&point,temperature_axes:["Ta","Tb"],liquid_fraction_axes:[None,None],
                mean_temperature_axis:None},caloric:[f.caloric[0].curve(),f.caloric[1].curve()],
            claim:pin,pressure_axis:"pressure"};
        run(cx,gate,f,t,&s,&contact,network)
    })
}

fn coefficients(bulk:Option<f64>)->Result<StepCoefficients<()>,Infallible> {
    Ok(StepCoefficients {sources_w:BTreeMap::from([(CELLS[0],0.001)]),
        conductances_w_k:bulk.map(|g|BTreeMap::from([(91,g)])).unwrap_or_default(),payload:()})
}

#[test]
fn finite_bulk_and_solved_pressure_match_two_cell_backward_euler() {
    with_support(None,8e6,|cx,_,_,t,s,contact,_| {
        let mut outcomes=Vec::new();
        for (compression,bulk) in [(0.001,1e-4),(0.002,1e-4),(0.001,2e-4)] {
            let next=s.advance_with_series_thermal_contacts(t,std::slice::from_ref(contact),1.0,
                &displacements(compression,0.0),&[],IncrementSettings::default(),support_controls(),
                &[0.0;2],|_|coefficients(Some(bulk)),cx).unwrap();
            let sample=&next.report().contacts[0];
            let pressure=2e9*compression;
            let resistance=0.5-0.4*pressure/8e6;
            let g=1.0/(1.0/bulk+resistance/LINK_AREA);
            let capacity=MASS*5.0;
            let ta=300.0+0.001/(capacity+g*capacity/(capacity+g));
            let tb=300.0+g/(capacity+g)*(ta-300.0);
            close(sample.compression_pa,pressure,1e-5);
            close(sample.resistance.value_m2_k_per_w(),resistance,1e-12);
            assert_eq!(sample.bulk_conductance_w_k,Some(bulk));
            close(sample.conductance_w_k,g,1e-15);
            close(next.thermal().temperature_k(CELLS[0]).unwrap(),ta,1e-7);
            close(next.thermal().temperature_k(CELLS[1]).unwrap(),tb,1e-7);
            close(sample.heat_a_to_b_j,capacity*(tb-300.0),1e-10);
            close(next.report().balance.external_heat_j,0.001,1e-10);
            assert!(next.report().balance.combined_residual_j.abs()<1e-8);
            assert_eq!(next.solid().solution().element_states[0].alpha,0.0);
            outcomes.push(tb);
        }
        assert!(outcomes[1]>outcomes[0]);assert!(outcomes[2]>outcomes[0]);
        assert_eq!(s.time_s(),0.0);assert_eq!(t.time_s(),0.0);
    });
}

#[test]
fn original_material_conductivity_and_contact_sources_both_survive_series_composition() {
    with_support(None,8e6,|cx,_,_,t,s,contact,network| {
        let prepared=network.electrical.bulk.prepare(t,||false).unwrap();
        let next=s.advance_with_series_thermal_contacts(t,std::slice::from_ref(contact),1.0,
            &displacements(0.001,0.0),&[],IncrementSettings::default(),support_controls(),&[0.0;2],
            |endpoint|prepared.resolve(endpoint,||false).map(|c|StepCoefficients {
                sources_w:BTreeMap::from([(CELLS[0],0.001)]),
                conductances_w_k:c.conductances_w_k,payload:c.source}),cx).unwrap();
        let bulk=next.report().source.links[0].conductance_w_k;
        let sample=&next.report().contacts[0];
        assert_eq!(sample.bulk_conductance_w_k,Some(bulk));
        close(bulk,LINK_AREA/(L/(2.0*0.01)),1e-15);
        close(sample.conductance_w_k,1.0/(1.0/bulk+sample.resistance.value_m2_k_per_w()/LINK_AREA),1e-15);
        assert_eq!(sample.resistance.receipt().selected,contact.claim);
        assert!(sample.conductance_w_k<bulk);
    });
}

#[test]
fn isolation_missing_bulk_and_old_interface_only_ownership_are_distinct() {
    with_support(None,8e6,|cx,gate,_,t,s,contact,_| {
        let supports=displacements(0.001,0.0);
        for bulk in [None,Some(-1.0),Some(f64::NAN)] {
            assert!(s.advance_with_series_thermal_contacts(t,std::slice::from_ref(contact),1.0,
                &supports,&[],IncrementSettings::default(),support_controls(),&[0.0;2],
                |_|coefficients(bulk),cx).is_err());
        }
        let isolated=s.advance_with_series_thermal_contacts(t,std::slice::from_ref(contact),1.0,
            &supports,&[],IncrementSettings::default(),support_controls(),&[0.0;2],
            |_|coefficients(Some(0.0)),cx).unwrap();
        assert_eq!(isolated.report().contacts[0].conductance_w_k,0.0);
        assert_eq!(isolated.report().contacts[0].heat_a_to_b_j,0.0);
        close(isolated.thermal().temperature_k(CELLS[1]).unwrap(),300.0,1e-10);
        let old=s.advance_with_thermal_contacts(t,std::slice::from_ref(contact),1.0,
            &supports,&[],IncrementSettings::default(),support_controls(),&[0.0;2],
            |_|coefficients(None),cx).unwrap();
        assert_eq!(old.report().contacts[0].bulk_conductance_w_k,None);
        close(old.report().contacts[0].conductance_w_k,LINK_AREA/0.4,1e-15);
        assert!(s.advance_with_thermal_contacts(t,std::slice::from_ref(contact),1.0,
            &supports,&[],IncrementSettings::default(),support_controls(),&[0.0;2],
            |_|coefficients(Some(1e-4)),cx).is_err());
        gate.request();
        assert!(matches!(s.advance_with_series_thermal_contacts(t,std::slice::from_ref(contact),1.0,
            &supports,&[],IncrementSettings::default(),support_controls(),&[0.0;2],
            |_|coefficients(Some(1e-4)),cx),Err(ContactError::Cancelled)));
        assert_eq!(s.time_s(),0.0);assert_eq!(t.time_s(),0.0);
    });
}
