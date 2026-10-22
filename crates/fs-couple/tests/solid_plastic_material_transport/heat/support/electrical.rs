//! G0/G1/G3: pressure-dependent cooling changes finite current and loaded J2 response.
use super::*;
use fs_couple::circuit_enthalpy::{CircuitEnthalpyError,CircuitEnthalpySolve,ResistorCell};
use fs_couple::circuit_enthalpy::material::CellResistor;
use fs_couple::material_circuit::{MaterialCircuit,MaterialCircuitSolve};
use fs_couple::solid_enthalpy::plastic::electrical::*;
use fs_material::conductor::{ElectricalPropertyBasis,ELECTRICAL_RESISTIVITY_DIMS};
use fs_phs::circuit::{Branch,CircuitGraph};

const CAP:f64=0.01;
const CHARGE:f64=0.01;
const CURRENT_AREA:f64=VOLUME/L;
fn powered_controls()->ElectricalSolidSolve {
    ElectricalSolidSolve {coupled:CircuitEnthalpySolve {
        maximum_cells:2,maximum_circuit_trials:1000,maximum_newton_iterations:32,
        difference_step_j:1e-7,cell_energy_tolerance_j:1e-10,thermal_energy_tolerance_j:1e-10,
        total_energy_tolerance_j:1e-8,circuit:MaterialCircuitSolve {
            energy_tolerance_j:1e-10,power_tolerance_w:1e-10},
        newton:NewtonKrylovConfig {absolute_tolerance:1e-12,relative_tolerance:1e-12,
            forcing_minimum:1e-4,linear_restart:2,..NewtonKrylovConfig::default()}},
        maximum_mechanical_trials:1000,yield_tolerance_pa:1e-3}
}
fn powered<R>(rho:f64,maximum_pressure:f64,duplicate:bool,charge:f64,
    run:impl FnOnce(&Cx<'_>,&CancelGate,&DepositionSimulation<'_>,&ThermoplasticSolid<'_>,
        &ThermalSupportContact<'_>,&PreparedHeatNetwork<'_,'_>,&MaterialCircuit)->R)->R {
    let mut property=scalar("electrical_resistivity",ELECTRICAL_RESISTIVITY_DIMS,rho);
    property.value=PropertyValue::Curve {abscissa:"T".into(),abscissa_dims:Temperature::DIMS,
        knots:vec![(280.0,0.8*rho),(320.0,1.2*rho)],dims:ELECTRICAL_RESISTIVITY_DIMS};
    property.validity=ValidityDomain::unconstrained().with("T",280.0,320.0);
    property.interpolation=InterpolationPolicy::LinearInside;
    with_support(Some(property),maximum_pressure,|cx,gate,_,t,s,support,old| {
        let resistors=[CellResistor {destination:ResistorCell {branch:1,cell:CELLS[0]},
            selection:MaterialPropertySelection::SingleClaimOnly,length_m:L,area_m2:CURRENT_AREA}];
        let bulk=MaterialNetwork {cells:old.electrical.bulk.cells,
            links:old.electrical.bulk.links,resistors:&resistors};
        let bases=[(1,ElectricalPropertyBasis::Resistivity)];
        let electrical=ElectricalContactNetwork {bulk:&bulk,contacts:&[],electrical_bases:&bases};
        let heat=MaterialHeatNetwork {electrical:&electrical,
            thermal_contacts:if duplicate {old.thermal_contacts}else{&[]},radiation:old.radiation};
        let prepared=heat.prepare(cx,t).unwrap();
        let initial=prepared.resolve(cx,t).unwrap();
        let circuit=MaterialCircuit::new(CircuitGraph {node_count:2,branches:vec![
            (1,0,Branch::Capacitor {farads:CAP}),(1,0,Branch::Resistor {ohms:0.0})],
            transformers:vec![]},initial.bulk.resistors,vec![],vec![charge],&[]).unwrap();
        run(cx,gate,t,s,support,&prepared,&circuit)
    })
}

// Independent RC midpoint, two-cell BE and scalar J2 traction equations.
fn oracle(rho:f64,compression:f64,maximum_pressure:f64,charge:f64)->[f64;4] {
    let capacity=MASS*5.0;
    let pressure=2e9*compression;
    let interface=0.5-0.4*pressure/maximum_pressure;
    let bulk=LINK_AREA/(L/(2.0*0.01));
    let g=1.0/(1.0/bulk+interface/LINK_AREA);
    let radiation=5.670374419e-8*L*L/(2.0/0.7-1.0);
    let exchange=|a:f64,b:f64|g*(a-b)+radiation*(a.powi(4)-b.powi(4));
    let lower=|a| {
        let (mut lo,mut hi)=(300.0,a);
        for _ in 0..80 {let b=0.5*(lo+hi);
            if capacity*(b-300.0)<exchange(a,b) {lo=b;}else{hi=b;}}
        0.5*(lo+hi)
    };
    let q=|a| {
        let r=rho*(1.0+0.01*(a-300.0))*L/CURRENT_AREA;
        let z=1.0/(2.0*r*CAP);charge*(1.0-z)/(1.0+z)
    };
    let residual=|a| {
        let (_,alpha,fy)=apex_oracle(a);
        capacity*(a-300.0)+exchange(a,lower(a))-VOLUME*fy*alpha
            -(charge*charge-q(a)*q(a))/(2.0*CAP)
    };
    let (mut lo,mut hi)=(300.0,319.0);
    assert!(residual(lo)<0.0 && residual(hi)>0.0);
    for _ in 0..80 {let mid=0.5*(lo+hi);if residual(mid)<0.0 {lo=mid;}else{hi=mid;}}
    let a=0.5*(lo+hi);[a,lower(a),2.0*L*apex_oracle(a).0,q(a)]
}

#[test]
fn compression_changes_current_and_plastic_displacement_through_the_original_heat_path() {
    let mut outcomes=Vec::new();
    for (compression,maximum_pressure) in [(0.001,8e6),(0.003,8e6),(0.001,16e6)] {
        powered(5.0,maximum_pressure,false,CHARGE,|cx,_,t,s,support,p,c| {
            let next=s.advance_electrical_with_supports_candidate(t,c,p,std::slice::from_ref(support),
                &[],1.0,&displacements(compression,SHEAR),&applied(),IncrementSettings::default(),
                powered_controls(),&[0.0;2],cx).unwrap();
            let expected=oracle(5.0,compression,maximum_pressure,CHARGE);
            let actual=[next.thermal().temperature_k(CELLS[0]).unwrap(),
                next.thermal().temperature_k(CELLS[1]).unwrap(),next.solid().solution().displacement_m[3][0],
                next.circuit().charges_c()[0]];
            for i in 0..2 {close(actual[i],expected[i],3e-6);}
            close(actual[2],expected[2],3e-11);close(actual[3],expected[3],1e-10);
            let r=next.report();let contact=&r.support_contacts[0];
            close(contact.compression_pa,2e9*compression,1e-5);
            assert_eq!(contact.resistance.receipt().selected,support.claim);
            let bulk=r.circuit_thermal.source.material.material.links[0].conductance_w_k;
            assert_eq!(contact.bulk_conductance_w_k,Some(bulk));
            close(contact.conductance_w_k,1.0/(1.0/bulk+contact.resistance.value_m2_k_per_w()/LINK_AREA),1e-15);
            assert!(r.circuit_thermal.source.material.thermal_contacts.is_empty());
            assert!(r.physical_plastic_heat_j>0.0 && r.circuit_thermal.deposited_joule_j>0.0);
            assert!(r.combined_residual_j.abs()<=1e-8 && r.thermal_residual_j.abs()<=1e-8);
            close(r.external_heat_j,0.0,1e-12);
            assert_eq!(t.time_s(),0.0);assert_eq!(s.time_s(),0.0);assert_eq!(c.charges_c(),&[CHARGE]);
            outcomes.push(actual);
        });
    }
    assert!(outcomes[1][0]<outcomes[0][0]-1e-4);
    assert!(outcomes[1][1]>outcomes[0][1]+1e-4);
    assert!((outcomes[1][2]-outcomes[0][2]).abs()>1e-10);
    assert!((outcomes[1][3]-outcomes[0][3]).abs()>1e-12);
    // The third run changes only the interface source curve, not the load,
    // bulk material, geometry or electrical supply.
    assert!(outcomes[2][0]>outcomes[0][0]+1e-4);
    assert!(outcomes[2][1]<outcomes[0][1]-1e-4);
    assert!((outcomes[2][2]-outcomes[0][2]).abs()>1e-10);
    assert!((outcomes[2][3]-outcomes[0][3]).abs()>1e-12);
}

#[test]
fn both_charge_signs_replay_the_native_circuit_and_series_contact_solver() {
    for charge in [CHARGE,-CHARGE] {
        powered(5.0,8e6,false,charge,|cx,_,t,s,support,p,c| {
            let u=displacements(0.002,SHEAR);
            let next=s.advance_electrical_with_supports_candidate(t,c,p,std::slice::from_ref(support),
                &[],1.0,&u,&applied(),IncrementSettings::default(),powered_controls(),&[0.0;2],cx).unwrap();
            let endpoint=p.resolve(cx,next.thermal()).unwrap();
            let native=c.step_candidate(endpoint.bulk.resistors,&[],1.0,powered_controls().coupled.circuit,
                ||false).unwrap();
            assert_eq!(next.circuit(),native.circuit());
            let heat=next.report().circuit_thermal.deposited_joule_j;
            let replay=s.advance_with_series_thermal_contacts(t,std::slice::from_ref(support),1.0,&u,
                &applied(),IncrementSettings::default(),support_controls(),&[0.0;2],|state| {
                    let mut k=p.resolve(cx,state)?.bulk;
                    *k.sources_w.entry(CELLS[0]).or_insert(0.0)+=heat;
                    Ok::<_,HeatNetworkError>(StepCoefficients {sources_w:k.sources_w,
                        conductances_w_k:k.conductances_w_k,payload:k.source})
                },cx).unwrap();
            for id in CELLS {close(next.thermal().temperature_k(id).unwrap(),
                replay.thermal().temperature_k(id).unwrap(),3e-6);}
            close(next.solid().solution().displacement_m[3][0],replay.solid().solution().displacement_m[3][0],3e-11);
            close(next.report().support_contacts[0].heat_a_to_b_j,replay.report().contacts[0].heat_a_to_b_j,1e-9);
            assert!(next.circuit().charges_c()[0].abs()<charge.abs());
        });
    }
}

#[test]
fn pressure_unloading_preserves_finite_charge_and_plastic_history() {
    powered(5.0,8e6,false,CHARGE,|cx,_,t,s,support,p,c| {
        let first=s.advance_electrical_with_supports_candidate(t,c,p,std::slice::from_ref(support),
            &[],1.0,&displacements(0.002,SHEAR),&applied(),IncrementSettings::default(),
            powered_controls(),&[0.0;2],cx).unwrap();
        let old=first.solid().solution().element_states[0].alpha;
        let next=first.solid().advance_electrical_with_supports_candidate(first.thermal(),first.circuit(),p,
            std::slice::from_ref(support),&[],2.0,&displacements(0.001,SHEAR-0.0005),&[],
            IncrementSettings::default(),powered_controls(),&[0.0;2],cx).unwrap();
        assert!(next.solid().solution().element_states[0].alpha>=old && old>0.0);
        assert!(next.report().support_contacts[0].compression_pa<first.report().support_contacts[0].compression_pa);
        assert!(next.circuit().charges_c()[0]<first.circuit().charges_c()[0]);
        assert!(next.report().combined_residual_j.abs()<=1e-8);
        assert_eq!(first.solid().time_s(),1.0);assert_eq!(c.charges_c(),&[CHARGE]);
    });
}

#[test]
fn late_pressure_failure_duplicate_correction_and_work_exhaustion_publish_no_history() {
    powered(5.0,1e6,false,CHARGE,|cx,_,t,s,support,p,c| {
        let saved=t.cells().clone();let history=s.solution().clone();let charge=c.clone();
        assert!(s.advance_electrical_with_supports_candidate(t,c,p,std::slice::from_ref(support),
            &[],1.0,&displacements(0.002,SHEAR),&applied(),IncrementSettings::default(),
            powered_controls(),&[0.0;2],cx).is_err());
        assert_eq!(t.cells(),&saved);assert_eq!(s.solution(),&history);assert_eq!(c,&charge);
    });
    powered(5.0,8e6,true,CHARGE,|cx,_,t,s,support,p,c| {
        assert!(matches!(s.advance_electrical_with_supports_candidate(t,c,p,std::slice::from_ref(support),
            &[],1.0,&displacements(0.001,SHEAR),&applied(),IncrementSettings::default(),
            powered_controls(),&[0.0;2],cx),Err(ElectricalSolidError::Invalid(_))));
        assert_eq!(s.time_s(),0.0);assert_eq!(c.charges_c(),&[CHARGE]);
    });
    powered(5.0,8e6,false,CHARGE,|cx,gate,t,s,support,p,c| {
        let u=displacements(0.001,SHEAR);
        let good=s.advance_electrical_with_supports_candidate(t,c,p,std::slice::from_ref(support),
            &[],1.0,&u,&applied(),IncrementSettings::default(),powered_controls(),&[0.0;2],cx).unwrap();
        let mut short=powered_controls();short.maximum_mechanical_trials=good.report().mechanical_trials-1;
        assert!(matches!(s.advance_electrical_with_supports_candidate(t,c,p,std::slice::from_ref(support),
            &[],1.0,&u,&applied(),IncrementSettings::default(),short,&[0.0;2],cx),
            Err(ElectricalSolidError::Coupled(CircuitEnthalpyError::Source(ElectricalSolidSourceError::MechanicalBudget)))));
        gate.request();
        assert!(matches!(s.advance_electrical_with_supports_candidate(t,c,p,std::slice::from_ref(support),
            &[],1.0,&u,&applied(),IncrementSettings::default(),powered_controls(),&[0.0;2],cx),
            Err(ElectricalSolidError::Cancelled)));
        assert_eq!(t.time_s(),0.0);assert_eq!(s.time_s(),0.0);assert_eq!(c.charges_c(),&[CHARGE]);
    });
}

#[test]
fn empty_support_list_retains_the_original_powered_arithmetic() {
    powered(5.0,8e6,false,CHARGE,|cx,_,t,s,_,p,c| {
        let u=displacements(0.001,SHEAR);
        let old=s.advance_electrical_candidate(t,c,p,&[],1.0,&u,&applied(),
            IncrementSettings::default(),powered_controls(),&[0.0;2],cx).unwrap();
        let new=s.advance_electrical_with_supports_candidate(t,c,p,&[],&[],1.0,&u,&applied(),
            IncrementSettings::default(),powered_controls(),&[0.0;2],cx).unwrap();
        assert_eq!(new.circuit(),old.circuit());assert_eq!(new.thermal().cells(),old.thermal().cells());
        assert_eq!(new.solid().solution(),old.solid().solution());assert!(new.report().support_contacts.is_empty());
    });
}
