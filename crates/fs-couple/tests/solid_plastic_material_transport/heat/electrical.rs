//! G0/G3: finite supply and actual native tetrahedra, never a surrogate spring.
#[path = "electrical/dynamics.rs"]
mod dynamics;
#[path = "electrical/geometry.rs"]
mod geometry;
use super::*;
use fs_couple::circuit_enthalpy::{CircuitEnthalpyError, CircuitEnthalpySolve, ResistorCell};
use fs_couple::circuit_enthalpy::material::CellResistor;
use fs_couple::circuit_enthalpy::contacts::material::MaterialElectricalContact;
use fs_couple::material_circuit::{MaterialCircuit, MaterialCircuitSolve};
use fs_couple::solid_enthalpy::plastic::electrical::*;
use fs_material::conductor::{ElectricalPropertyBasis, ELECTRICAL_RESISTIVITY_DIMS};
use fs_material::electrical_interface::{ELECTRICAL_CONTACT_RESISTANCE_PROPERTY as ER,
    ELECTRICAL_CONTACT_RESISTANCE_DIMS as ER_DIMS};
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_solver::NewtonKrylovConfig;
use fs_solid::plastic3::{SmallStrainTetProblem, SmallStrainLawField};

const CAPACITANCE: f64 = 0.01;
const CHARGE: f64 = 0.01;
const ELECTRICAL_AREA: f64 = VOLUME / L;
const EMISSIVITY: f64 = 0.7;
const THERMAL_CONTACT: f64 = 0.1;

fn controls() -> ElectricalSolidSolve {
    ElectricalSolidSolve {
        coupled: CircuitEnthalpySolve {
            maximum_cells: 2, maximum_circuit_trials: 1000, maximum_newton_iterations: 32,
            difference_step_j: 1e-7, cell_energy_tolerance_j: 1e-10,
            thermal_energy_tolerance_j: 1e-10, total_energy_tolerance_j: 1e-8,
            circuit: MaterialCircuitSolve { energy_tolerance_j: 1e-10, power_tolerance_w: 1e-10 },
            newton: NewtonKrylovConfig { absolute_tolerance: 1e-12, relative_tolerance: 1e-12,
                forcing_minimum: 1e-4, linear_restart: 4, ..NewtonKrylovConfig::default() },
        }, maximum_mechanical_trials: 1000, yield_tolerance_pa: 1e-3,
    }
}

fn with_powered<R>(rho: f64, contact: bool, short_domain: bool, charge: f64,
    run: impl FnOnce(&Cx<'_>, &CancelGate, &DepositionSimulation<'_>, &ThermoplasticSolid<'_>,
        &MaterialCircuit, &PreparedHeatNetwork<'_, '_>) -> R) -> R {
    let mut property = scalar("electrical_resistivity", ELECTRICAL_RESISTIVITY_DIMS, rho);
    property.value = PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
        knots: vec![(280.0, 0.8*rho), (320.0, 1.2*rho)], dims: ELECTRICAL_RESISTIVITY_DIMS };
    property.interpolation = InterpolationPolicy::LinearInside;
    property.validity = ValidityDomain::unconstrained().with("T", 280.0, if short_domain {300.05}else{320.0});
    with_heat_fixture(fixture_with_property(EMISSIVITY, Some(property)), THERMAL_CONTACT,
        |cx, gate, _, thermal, solid, original| {
        let resistors = [CellResistor { destination: ResistorCell { branch: 1, cell: CELLS[0] },
            selection: MaterialPropertySelection::SingleClaimOnly, length_m: L, area_m2: ELECTRICAL_AREA }];
        let bulk = MaterialNetwork { cells: original.electrical.bulk.cells,
            links: original.electrical.bulk.links, resistors: &resistors };
        let bases = [(1, ElectricalPropertyBasis::Resistivity)];
        let old_interface = original.thermal_contacts[0].point;
        let mut claims = old_interface.interface.claims().clone();
        claims.insert_claim(scalar(ER, ER_DIMS, 0.005)).unwrap();
        let interface = InterfaceSystemCard::assemble(
            old_interface.interface.surface_a().clone(), old_interface.interface.surface_b().clone(),
            SystemContext { medium: "fixed synthetic interface".into(), third_body: None,
                environment: "test".into(), history: "fixed".into() }, claims, vec![]).unwrap();
        let point = ThermalInterfacePoint { interface: &interface, ..old_interface };
        let contacts = [MaterialElectricalContact { link: 91, branch: 1, point,
            selection: MaterialPropertySelection::SingleClaimOnly, heat_fraction_to_first: 0.3 }];
        let electrical = ElectricalContactNetwork { bulk: &bulk,
            contacts: if contact {&contacts}else{&[]}, electrical_bases: &bases };
        let thermal_contact = [MaterialContactLink { point, ..original.thermal_contacts[0] }];
        let network = MaterialHeatNetwork { electrical: &electrical,
            thermal_contacts: &thermal_contact, radiation: original.radiation };
        let prepared = network.prepare(cx, thermal).unwrap();
        let initial = prepared.resolve(cx, thermal).unwrap();
        let circuit = MaterialCircuit::new_with_contacts(CircuitGraph { node_count: 2,
            branches: vec![(1,0,Branch::Capacitor { farads: CAPACITANCE }),
                (1,0,Branch::Resistor { ohms: 0.0 })], transformers: vec![] },
            initial.bulk.resistors, initial.contacts, vec![], vec![charge], &[]).unwrap();
        run(cx, gate, thermal, solid, &circuit, &prepared)
    })
}

// Independent closed-form J2 reduction + capacitor midpoint + two-cell BE heat.
// It uses neither the native solid solve nor the circuit/thermal residual owner.
fn reference(rho: f64, charge: f64) -> [f64; 4] {
    let capacity = MASS * 5.0;
    let bulk = LINK_AREA / (L / (2.0 * 0.01));
    let conductance = 1.0 / (1.0/bulk + THERMAL_CONTACT/LINK_AREA);
    let radiation = 5.670374419e-8 * L*L / (2.0/EMISSIVITY - 1.0);
    let exchange = |ta: f64, tb: f64| conductance*(ta-tb) + radiation*(ta.powi(4)-tb.powi(4));
    let lower_temperature = |ta| {
        let (mut lo, mut hi) = (300.0, ta);
        for _ in 0..80 {
            let tb = 0.5*(lo+hi);
            if capacity*(tb-300.0) < exchange(ta,tb) {lo=tb;} else {hi=tb;}
        }
        0.5*(lo+hi)
    };
    let q_at = |ta| {
        let r = rho*(1.0+0.01*(ta-300.0))*L/ELECTRICAL_AREA;
        let z = 1.0/(2.0*r*CAPACITANCE);
        charge*(1.0-z)/(1.0+z)
    };
    let residual = |ta| {
        let (_, alpha, fy) = apex_oracle(ta);
        let q = q_at(ta);
        let joule = (charge*charge-q*q)/(2.0*CAPACITANCE);
        capacity*(ta-300.0)-VOLUME*fy*alpha-joule+exchange(ta,lower_temperature(ta))
    };
    let (mut lo, mut hi) = (300.0,319.0);
    assert!(residual(lo)<0.0 && residual(hi)>0.0);
    for _ in 0..80 {
        let mid=0.5*(lo+hi); if residual(mid)<0.0 {lo=mid;} else {hi=mid;}
    }
    let ta=0.5*(lo+hi);
    [ta, lower_temperature(ta), 2.0*L*apex_oracle(ta).0, q_at(ta)]
}

#[test]
fn material_resistivity_changes_actual_plastic_displacement_and_matches_independent_equations() {
    let mut observations = Vec::new();
    for rho in [5.0,20.0] {
        with_powered(rho,false,false,CHARGE,|cx,_,t,s,c,p| {
            let next=s.advance_electrical_candidate(t,c,p,&[],1.0,&supports(SHEAR),&applied(),
                IncrementSettings::default(),controls(),&[0.0;2],cx).unwrap();
            let expected=reference(rho,CHARGE);
            let actual=[next.thermal().temperature_k(CELLS[0]).unwrap(),
                next.thermal().temperature_k(CELLS[1]).unwrap(),next.solid().solution().displacement_m[3][0],
                next.circuit().charges_c()[0]];
            for i in 0..2 {close(actual[i],expected[i],3e-6);}
            close(actual[2],expected[2],3e-11); close(actual[3],expected[3],1e-10);
            let report=next.report(); let electrical=&report.circuit_thermal;
            assert!(report.physical_plastic_heat_j>0.0 && electrical.deposited_joule_j>0.0);
            assert!(report.numerical_loss_j>report.physical_plastic_heat_j);
            close(electrical.thermal.stored_energy_change_j,
                report.physical_plastic_heat_j+electrical.deposited_joule_j,1e-8);
            assert_eq!(report.reversible_heat_j,0.0); close(report.external_heat_j,0.0,1e-12);
            assert!(report.combined_residual_j.abs()<=controls().coupled.total_energy_tolerance_j);
            assert!(report.thermal_residual_j.abs()<=controls().coupled.total_energy_tolerance_j);
            assert_eq!(t.time_s(),0.0); assert_eq!(s.time_s(),0.0); assert_eq!(c.charges_c(),&[CHARGE]);
            observations.push(actual);
        });
    }
    assert!(observations[0][0]>observations[1][0]+0.1);
    assert!(observations[0][2]>observations[1][2]+1e-8);
    assert!(observations[0][3]<observations[1][3]);
}

#[test]
fn contact_heating_replays_native_circuit_and_loaded_tetrahedra_at_the_same_endpoint() {
    for charge in [CHARGE,-CHARGE] {
        with_powered(10.0,true,false,charge,|cx,_,t,s,c,p| {
            let supports=supports(SHEAR); let forces=applied();
            let next=s.advance_electrical_candidate(t,c,p,&[],1.0,&supports,&forces,
                IncrementSettings::default(),controls(),&[0.0;2],cx).unwrap();
            let exact=p.resolve(cx,next.thermal()).unwrap();
            let original=c.step_with_contacts_candidate(exact.bulk.resistors,exact.contacts,&[],1.0,
                controls().coupled.circuit,||false).unwrap();
            assert_eq!(next.circuit(),original.circuit());
            assert_eq!(next.report().circuit_thermal.source.destinations[0].heat_fraction_to_first,0.3);
            assert_eq!(next.report().circuit_thermal.electrical.contacts.len(),1);
            assert!(next.report().circuit_thermal.electrical.contacts[0].heat_j>0.0);
            let problem=SmallStrainTetProblem { nodes_m:&NODES,tetrahedra:&TETS,
                laws:SmallStrainLawField::PerElement(next.solid().materials()),prescribed_m:&supports,
                nodal_forces_n:&forces,load_path:&[1.0],budget:TetAssemblyBudget::standard(),
                settings:IncrementSettings::default() };
            let native=problem.solve_from(s.solution().displacement_m.as_flattened(),
                &s.solution().element_states,cx).unwrap();
            for (a,b) in native.displacement_m.as_flattened().iter()
                .zip(next.solid().solution().displacement_m.as_flattened()) {close(*a,*b,1e-11);}
            for (a,b) in native.element_states.iter().zip(&next.solid().solution().element_states) {
                close(a.alpha,b.alpha,1e-10);
            }
            assert_eq!(next.report().circuit_thermal.source.material.radiation,exact.bulk.source.radiation);
            assert!(next.circuit().charges_c()[0].abs()<charge.abs());
        });
    }
}

#[test]
fn complete_heat_trait_object_preserves_original_contact_and_solid_trajectory() {
    use fs_couple::circuit_enthalpy::contacts::heat::{HeatAssembly, HeatNetworkSample};
    with_powered(10.0, true, false, CHARGE, |cx, _, t, s, c, p| {
        let assembly: &dyn HeatAssembly<Sample = HeatNetworkSample, Error = HeatNetworkError> = p;
        let original = s.advance_electrical_candidate(t, c, p, &[], 1.0,
            &supports(SHEAR), &applied(), IncrementSettings::default(), controls(), &[0.0; 2], cx).unwrap();
        let complete = s.advance_electrical_candidate(t, c, assembly, &[], 1.0,
            &supports(SHEAR), &applied(), IncrementSettings::default(), controls(), &[0.0; 2], cx).unwrap();
        assert_eq!(complete.thermal().cells(), original.thermal().cells());
        assert_eq!(complete.solid().solution(), original.solid().solution());
        assert_eq!(complete.circuit(), original.circuit());
        assert_eq!(complete.report().circuit_thermal.source.destinations,
            original.report().circuit_thermal.source.destinations);
        assert_eq!(complete.report().combined_residual_j, original.report().combined_residual_j);
    });
}

#[test]
fn unloading_keeps_plastic_and_charge_history_without_thermalizing_algorithmic_loss() {
    with_powered(10.0,true,false,CHARGE,|cx,_,t,s,c,p| {
        let first=s.advance_electrical_candidate(t,c,p,&[],1.0,&supports(SHEAR),&applied(),
            IncrementSettings::default(),controls(),&[0.0;2],cx).unwrap();
        let previous=first.solid().solution().element_states.clone();
        let next=first.solid().advance_electrical_candidate(first.thermal(),first.circuit(),p,&[],2.0,
            &supports(SHEAR-0.0005),&[],IncrementSettings::default(),controls(),&[0.0;2],cx).unwrap();
        assert!(next.circuit().charges_c()[0]<first.circuit().charges_c()[0]);
        assert!(next.solid().solution().element_states[0].alpha>=previous[0].alpha);
        assert!(next.solid().solution().element_states[0].alpha>0.0);
        assert!(next.report().numerical_loss_j>0.0);
        close(next.report().circuit_thermal.thermal.stored_energy_change_j,
            next.report().deposited_mechanical_heat_j+next.report().circuit_thermal.deposited_joule_j,1e-8);
        assert!(next.report().combined_residual_j.abs()<=controls().coupled.total_energy_tolerance_j);
        assert_eq!(next.solid().time_s(),next.thermal().time_s());
        assert_eq!(first.solid().time_s(),1.0); assert_eq!(c.charges_c(),&[CHARGE]);
    });
}

#[test]
fn final_trial_budget_and_cancellation_leave_all_original_histories_replayable() {
    with_powered(10.0,true,false,CHARGE,|cx,gate,t,s,c,p| {
        let good=s.advance_electrical_candidate(t,c,p,&[],1.0,&supports(SHEAR),&applied(),
            IncrementSettings::default(),controls(),&[0.0;2],cx).unwrap();
        let mut short=controls();short.maximum_mechanical_trials=good.report().mechanical_trials-1;
        assert!(matches!(s.advance_electrical_candidate(t,c,p,&[],1.0,&supports(SHEAR),&applied(),
            IncrementSettings::default(),short,&[0.0;2],cx),
            Err(ElectricalSolidError::Coupled(CircuitEnthalpyError::Source(ElectricalSolidSourceError::MechanicalBudget)))));
        let again=s.advance_electrical_candidate(t,c,p,&[],1.0,&supports(SHEAR),&applied(),
            IncrementSettings::default(),controls(),&[0.0;2],cx).unwrap();
        assert_eq!(again.circuit(),good.circuit());assert_eq!(again.thermal().cells(),good.thermal().cells());
        assert_eq!(again.solid().solution(),good.solid().solution());
        gate.request();
        assert!(matches!(s.advance_electrical_candidate(t,c,p,&[],1.0,&supports(SHEAR),&applied(),
            IncrementSettings::default(),controls(),&[0.0;2],cx),Err(ElectricalSolidError::Cancelled)));
        assert_eq!(t.time_s(),0.0);assert_eq!(s.time_s(),0.0);assert_eq!(c.charges_c(),&[CHARGE]);
    });
}

#[test]
fn electrical_source_expiration_rejects_the_thermal_and_plastic_candidate() {
    with_powered(10.0,false,true,CHARGE,|cx,_,t,s,c,p| {
        assert!(p.resolve(cx,t).is_ok());
        let original=t.cells().clone();let history=s.solution().clone();let charge=c.clone();
        assert!(s.advance_electrical_candidate(t,c,p,&[],1.0,&supports(SHEAR),&applied(),
            IncrementSettings::default(),controls(),&[0.0;2],cx).is_err());
        assert_eq!(t.cells(),&original);assert_eq!(s.solution(),&history);assert_eq!(c,&charge);
    });
}
