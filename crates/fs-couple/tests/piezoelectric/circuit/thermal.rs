//! Synthetic numerical ceramic: independent constitutive and midpoint references.
//! These intentionally scaled numbers are not measured ceramic calibration.

mod heat;
use super::{near, thermal_controls};
use fs_couple::circuit_enthalpy::{CircuitEnthalpyError, ResistorCell};
use fs_couple::circuit_enthalpy::material::{CellResistor, MaterialCell, MaterialNetwork};
use fs_couple::circuit_enthalpy::piezoelectric::thermal::{
    self as heated, PiezoelectricThermalBinding,
};
use fs_couple::material_circuit::piezoelectric::PiezoelectricCircuit;
use fs_couple::piezoelectric::{PiezoelectricRod, SmallSignalLimits};
use fs_evidence::ValidityDomain;
use fs_material::conductor::ElectricalPropertyBasis;
use fs_material::field_storage::ELECTRIC_PERMITTIVITY_DIMS;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use fs_material::piezoelectric::{
    COMPLIANCE_33_E_PROPERTY, COMPLIANCE_DIMS, LongitudinalPiezoelectric,
    PERMITTIVITY_33_T_PROPERTY, STRAIN_CHARGE_33_PROPERTY, STRAIN_CHARGE_DIMS,
    thermal::ThermalLongitudinalPiezoelectric,
};
use fs_material::state_point::MaterialPropertySelection;
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_qty::{Density, QuantitySpec, Temperature, semantic::{QuantityKind, SemanticType, ValueForm}};
use fs_thermal::enthalpy::*;

const LENGTH: f64 = 2.0;
const AREA: f64 = 0.5;
const DENSITY: f64 = 3.0;
const CAPACITY: f64 = 0.05;
const CELL: CellId = CellId(7);
const MATERIAL: MaterialId = MaterialId(3);
const SLOPES: [f64; 3] = [0.0004, 0.0015, 0.006];

fn temperature_quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(QuantityKind::AbsoluteTemperature, ValueForm::Static))
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new().with_quantity("T", temperature_quantity(), t).unwrap()
}
fn source_card(d: f64, slopes: [f64; 3], density_slope: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    let basis = ElectricalPropertyBasis::Resistivity;
    for (name, dims, value, slope) in [
        (COMPLIANCE_33_E_PROPERTY, COMPLIANCE_DIMS, 0.1, slopes[0]),
        (STRAIN_CHARGE_33_PROPERTY, STRAIN_CHARGE_DIMS, d, slopes[1]),
        (PERMITTIVITY_33_T_PROPERTY, ELECTRIC_PERMITTIVITY_DIMS, 2.0, slopes[2]),
        ("density", Density::DIMS, DENSITY, density_slope),
        (basis.property_name(), basis.dimensions(), 0.25, 0.0025),
    ] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Curve {
                abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                knots: vec![(280.0, value-20.0*slope), (360.0, value+60.0*slope)], dims,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity("T", temperature_quantity(), 280.0, 360.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside, observations: vec![],
            provenance: Provenance { source: "authored affine coupled ceramic fixture".into(),
                license: "test-only".into(), artifact: None },
        }).unwrap();
    }
    MaterialCard::assemble(MaterialStateId {
        chemistry: "synthetic thermal piezoelectric".into(), phase: "solid".into(),
        process: "numerical reciprocal reduction".into(), revision: 0,
    }, claims, vec![]).unwrap()
}
fn ceramic(card: &MaterialCard) -> LongitudinalPiezoelectric {
    LongitudinalPiezoelectric::resolve(card, &point(300.0),
        MaterialPropertySelection::SingleClaimOnly).unwrap()
}
fn thermal_binding(card: &MaterialCard, range: [f64; 2]) -> PiezoelectricThermalBinding {
    PiezoelectricThermalBinding {
        cell: CELL,
        law: ThermalLongitudinalPiezoelectric::new(card, &ceramic(card), "T", range).unwrap(),
        maximum_abs_displacement_c_m2: 0.22, volume_tolerance_m3: 1e-12,
    }
}
fn caloric_curve(card: &MaterialCard, phase: bool) -> EquilibriumEnthalpyPhaseCurve {
    let knots = vec![
        EnthalpyPhaseKnot { specific_enthalpy_j_kg: 0.0, temperature_k: 280.0,
            liquid_mass_fraction: 0.0, bulk_density_kg_m3: DENSITY },
        EnthalpyPhaseKnot { specific_enthalpy_j_kg: 80.0*CAPACITY/DENSITY, temperature_k: 360.0,
            liquid_mass_fraction: if phase {1.0} else {0.0}, bulk_density_kg_m3: DENSITY },
    ];
    if phase { EquilibriumEnthalpyPhaseCurve::try_new(card.content_hash(), knots).unwrap() }
    else { EquilibriumEnthalpyPhaseCurve::try_single_phase(
        card.content_hash(), SolidLiquidPhase::Solid, knots).unwrap() }
}
fn initial_thermal(curve: &EquilibriumEnthalpyPhaseCurve, power: f64) -> DepositionSimulation<'_> {
    DepositionSimulation::new(vec![TimedEvent { time_s: 0.0, event: Event::Birth(DepositedCell {
        id: CELL, segment_id: 0, body_id: 0, material_id: MATERIAL, mass_kg: DENSITY,
        specific_h_j_kg: 20.0*CAPACITY/DENSITY,
        geometry: CellGeometry { start_m: [0.0;3], end_m: [LENGTH,0.0,0.0],
            volume_m3: LENGTH*AREA, surface_area_m2: 6.0 },
        boundary: Boundary { ambient_k: 300.0, htc_w_m2_k: 0.0, reservoir: None, source_w: power },
    }) }], &[(MATERIAL, curve as &dyn ThermalLaw)], DepositionBudget {
        max_cells: 1, max_links: 0, max_events: 1, max_substeps: 100,
    }).unwrap().advance_to(0.0, StepConfig { max_dt_s: 0.1, energy_tolerance_j: 1e-9 },
        || Ok(())).unwrap().0
}
fn cell<'a>(card: &'a MaterialCard, curve: &'a EquilibriumEnthalpyPhaseCurve,
    context: &'a QueryPoint) -> MaterialCell<'a> {
    MaterialCell { cell: CELL, material_id: MATERIAL, curve,
        point: ThermalMaterialPoint { material: card, context,
            temperature_axis: "T", liquid_fraction_axis: None }, conductivity: None }
}
fn resistor() -> CellResistor {
    CellResistor { destination: ResistorCell { branch: 1, cell: CELL },
        selection: MaterialPropertySelection::SingleClaimOnly, length_m: LENGTH, area_m2: AREA }
}
fn model(card: &MaterialCard, material: &MaterialNetwork<'_>, thermal: &DepositionSimulation<'_>,
    x: f64, v: f64, q: f64) -> PiezoelectricCircuit {
    let mut branches = vec![(1,0,Branch::Capacitor { farads: f64::NAN })];
    if !material.resistors.is_empty() { branches.push((1,0,Branch::Resistor { ohms: f64::NAN })); }
    PiezoelectricCircuit::new(CircuitGraph { node_count: 2, branches, transformers: vec![] },
        material.resolve(thermal, || false).unwrap().resistors, vec![], vec![q], 0,
        PiezoelectricRod::new(ceramic(card), LENGTH, AREA, SmallSignalLimits {
            max_abs_strain: 0.03, max_abs_electric_field_v_m: 10.0,
        }).unwrap(), x, v, &[], || false).unwrap()
}

// Independent inverse of A=[[s,d],[d,epsilon]], evaluated only from authored numbers.
fn coefficients(t: f64, d: f64, slopes: [f64; 3]) -> [f64; 3] {
    [0.1+slopes[0]*(t-300.0), d+slopes[1]*(t-300.0), 2.0+slopes[2]*(t-300.0)]
}
fn potentials(t: f64, strain: f64, displacement: f64, d: f64, slopes: [f64; 3],
    volume: f64) -> [f64; 4] {
    let [s,d,e] = coefficients(t,d,slopes);
    let determinant = s*e-d*d;
    let stress = (e*strain-d*displacement)/determinant;
    let field = (s*displacement-d*strain)/determinant;
    let a = slopes[0]*stress+slopes[1]*field;
    let b = slopes[1]*stress+slopes[2]*field;
    let psi = 0.5*volume*(strain*stress+displacement*field);
    let entropy = 0.5*volume*(stress*a+field*b);
    [psi, entropy, psi+t*entropy,
        -t*volume*(e*a*a-2.0*d*a*b+s*b*b)/determinant]
}
// Hand elimination of the native three midpoint equations at the proposed T.
// No runtime matrices, source values, material response or native state helper
// define this reference. Voltage is the held-T constitutive midpoint voltage.
fn motion(t: f64, old: [f64; 3], d: f64, slopes: [f64; 3],
    force: f64, dt: f64) -> [f64; 5] {
    let [x,v,q] = old;
    let [s,d,e] = coefficients(t,d,slopes);
    let mass = DENSITY*LENGTH*AREA/3.0;
    let stiffness = AREA/(s*LENGTH);
    let capacitance = (e-d*d/s)*AREA/LENGTH;
    let theta = (d/s)*AREA/LENGTH;
    let resistance = 1.0+0.01*(t-300.0);
    let denominator = 1.0+dt*dt*stiffness/(4.0*mass);
    let a = (x+0.5*dt*v+dt*dt*force/(4.0*mass))/denominator;
    let b = dt*dt*theta/(4.0*mass*denominator);
    let voltage = (q-theta*a)/(capacitance+dt/(2.0*resistance)+theta*b);
    let middle_x = a+b*voltage;
    [2.0*middle_x-x, 4.0*(middle_x-x)/dt-v, q-dt*voltage/resistance,
        voltage, dt*voltage*voltage/resistance]
}
fn root(mut low: f64, mut high: f64, residual: impl Fn(f64) -> f64) -> f64 {
    assert!(residual(low)<0.0 && residual(high)>0.0);
    for _ in 0..100 {
        let middle = 0.5*(low+high);
        if residual(middle)<0.0 { low=middle; } else { high=middle; }
    }
    0.5*(low+high)
}
fn internal(t: f64, x: f64, q: f64, d: f64, slopes: [f64;3]) -> f64 {
    potentials(t,x/LENGTH,q/AREA,d,slopes,LENGTH*AREA)[2]
}
fn endpoint(t0: f64, old: [f64;3], d: f64, slopes: [f64;3],
    force: f64, dt: f64, heat: f64) -> (f64,[f64;5]) {
    let mass = DENSITY*LENGTH*AREA/3.0;
    let before = internal(t0,old[0],old[2],d,slopes)+0.5*mass*old[1]*old[1];
    let t = root(280.0,360.0, |t| {
        let [x,v,q,_,_] = motion(t,old,d,slopes,force,dt);
        CAPACITY*(t-t0)+internal(t,x,q,d,slopes)+0.5*mass*v*v-before-force*(x-old[0])-heat
    });
    (t,motion(t,old,d,slopes,force,dt))
}

#[test]
fn coupled_compliance_thermal_derivatives_keep_signed_and_zero_cross_coefficients() {
    for (d,slopes) in [(0.2,SLOPES),(-0.2,[SLOPES[0],-SLOPES[1],SLOPES[2]]),
        (0.0,[SLOPES[0],0.0,SLOPES[2]]), (0.0,SLOPES), (0.2,[0.0;3])] {
        let card = source_card(d,slopes,0.0); let law = thermal_binding(&card,[280.0,360.0]).law;
        let (t,strain,displacement,volume) = (321.0,0.013,-0.19,0.7);
        let actual = law.at_temperature(t).unwrap().state(strain,displacement,volume).unwrap();
        let expected = potentials(t,strain,displacement,d,slopes,volume);
        for (a,b) in [actual.free_energy_j,actual.entropy_j_k,
            actual.internal_energy_j,actual.heat_capacity_j_k].into_iter().zip(expected) {
            near(a,b,2e-14);
        }
        let low = law.at_temperature(t-0.001).unwrap().state(strain,displacement,volume).unwrap();
        let high = law.at_temperature(t+0.001).unwrap().state(strain,displacement,volume).unwrap();
        near(-(high.free_energy_j-low.free_energy_j)/0.002,actual.entropy_j_k,2e-11);
        near((high.internal_energy_j-low.internal_energy_j)/0.002,actual.heat_capacity_j_k,2e-11);
        assert!(actual.heat_capacity_j_k>=law.minimum_heat_capacity_j_k(0.03,0.22,volume).unwrap());
        assert!(actual.heat_capacity_j_k<=0.0);
        if slopes==[0.0;3] {
            assert_eq!(actual.entropy_j_k,0.0);
            assert_eq!(actual.heat_capacity_j_k,0.0);
            assert_eq!(actual.internal_energy_j,actual.free_energy_j);
        }
        // The affine d path passes through zero without a relative-error singularity.
        law.at_temperature(300.0).unwrap();
        assert_eq!(law.at_temperature(t).unwrap().material().material().card_identity(),card.content_hash());
    }
}

#[test]
fn ceramic_temperature_updates_native_load_motion_and_total_internal_energy() {
    let mut final_states = Vec::new();
    for (d,slopes) in [(0.2,SLOPES),(-0.2,[SLOPES[0],-SLOPES[1],SLOPES[2]]),
        (0.0,[SLOPES[0],0.0,SLOPES[2]]), (0.2,[0.0;3])] {
        let card=source_card(d,slopes,0.0); let curve=caloric_curve(&card,false); let point=point(300.0);
        let cells=[cell(&card,&curve,&point)]; let resistors=[resistor()];
        let network=MaterialNetwork { cells:&cells,resistors:&resistors,links:&[] };
        let binding=thermal_binding(&card,[280.0,360.0]);
        let mut thermal=initial_thermal(&curve,0.04);
        let mut old=[0.02,0.01,if d<0.0 {-0.1} else {0.1}]; let mut t=300.0;
        let mut body=model(&card,&network,&thermal,old[0],old[1],old[2]);
        let start=thermal.energy_j()+internal(t,old[0],old[2],d,slopes)+0.5*old[1]*old[1];
        let mut work=0.0; let mut heat=0.0;
        for interval in 1..=6 {
            let target=f64::from(interval)*0.1; let dt=target-thermal.time_s();
            let expected=endpoint(t,old,d,slopes,0.03,dt,0.04*dt);
            let before_body=body.clone(); let before_cells=thermal.cells().clone();
            let next=heated::step_candidate(&network,&thermal,target,&body,&[],0.03,&binding,
                &[(1,ElectricalPropertyBasis::Resistivity)],thermal_controls(),&[0.0],||false).unwrap();
            let report=next.report(); let coupled=&report.coupled;
            let (expected_t,[x,v,q,voltage,joule])=expected;
            near(next.thermal().temperature_k(CELL).unwrap(),expected_t,2e-8);
            near(next.model().extension_m(),x,2e-10);
            near(next.model().velocity_m_s(),v,2e-9);
            near(next.model().charge_c(),q,2e-10);
            near(coupled.electrical.resistors[0].voltage_v,voltage,2e-9);
            near(coupled.deposited_joule_j,joule,2e-10);
            near(coupled.electrical.mechanical_supply_j,0.03*(x-old[0]),2e-10);
            near(coupled.thermal.stored_energy_change_j,
                joule+coupled.source.reversible_heat_j+0.04*dt,2e-10);
            near(report.external_heat_j,0.04*dt,2e-10);
            near(report.total_energy_residual_j,0.0,2e-10);
            let before=potentials(t,old[0]/LENGTH,old[2]/AREA,d,slopes,1.0);
            let frozen=potentials(expected_t,old[0]/LENGTH,old[2]/AREA,d,slopes,1.0);
            let after=potentials(expected_t,x/LENGTH,q/AREA,d,slopes,1.0);
            near(coupled.electrical.parameter_work_j,frozen[0]-before[0],2e-10);
            near(coupled.source.reversible_heat_j,
                -(frozen[0]-before[0])-(expected_t*after[1]-t*before[1]),2e-10);
            let actual_t=coupled.source.temperature_k[1];
            assert_eq!(actual_t.to_bits(),next.thermal().temperature_k(CELL).unwrap().to_bits());
            assert_eq!(next.model().rod().material(),&coupled.source.material);
            assert_eq!(coupled.source.material.material().card_identity(),card.content_hash());
            assert_eq!(body,before_body); assert_eq!(thermal.cells(),&before_cells);
            if slopes==[0.0;3] {
                // With the SAME endpoint resistor data, fixed ceramic sources
                // take the exact existing native trajectory: no artificial work.
                let fixed=body.step_candidate(network.resolve(next.thermal(),||false).unwrap().resistors,
                    &[],0.03,dt,thermal_controls().circuit,||false).unwrap();
                assert_eq!(fixed.model().charge_c().to_bits(),next.model().charge_c().to_bits());
                assert_eq!(fixed.model().extension_m().to_bits(),next.model().extension_m().to_bits());
                assert_eq!(fixed.model().momentum_kg_m_s().to_bits(),next.model().momentum_kg_m_s().to_bits());
                assert_eq!(coupled.electrical.parameter_work_j,0.0);
                assert_eq!(coupled.source.reversible_heat_j,0.0);
            }
            work+=coupled.electrical.mechanical_supply_j; heat+=report.external_heat_j;
            near(next.thermal().energy_j()+internal(expected_t,x,q,d,slopes)+0.5*v*v,
                start+work+heat,3e-9);
            (thermal,body,_) = next.into_parts();
            old=[x,v,q]; t=expected_t;
        }
        final_states.push((t,body.extension_m(),body.charge_c()));
    }
    near(final_states[0].0,final_states[1].0,2e-8);
    near(final_states[0].1,final_states[1].1,2e-10);
    near(final_states[0].2,-final_states[1].2,2e-10);
    assert!((final_states[0].0-final_states[3].0).abs()>0.05);
    assert!((final_states[0].1-final_states[3].1).abs()>2e-5);
    assert!((final_states[0].2-final_states[3].2).abs()>4e-5);
}

#[test]
fn isolated_charge_uses_internal_energy_and_the_original_zero_coupling_limit() {
    let slopes=[0.0,0.0,0.006]; let card=source_card(0.0,slopes,0.0);
    let curve=caloric_curve(&card,false); let point=point(300.0); let thermal=initial_thermal(&curve,0.04);
    let cells=[cell(&card,&curve,&point)];
    let network=MaterialNetwork { cells:&cells,resistors:&[],links:&[] };
    let binding=thermal_binding(&card,[280.0,360.0]);
    let body=model(&card,&network,&thermal,0.0,0.0,0.1);
    let expected=root(280.0,360.0,|t| CAPACITY*(t-300.0)
        +internal(t,0.0,0.1,0.0,slopes)-internal(300.0,0.0,0.1,0.0,slopes)-0.004);
    let next=heated::step_candidate(&network,&thermal,0.1,&body,&[],0.0,&binding,
        &[],thermal_controls(),&[0.0],||false).unwrap();
    near(next.thermal().temperature_k(CELL).unwrap(),expected,2e-8);
    near(next.model().charge_c(),0.1,1e-14);
    near(next.model().extension_m(),0.0,1e-14);
    near(next.model().momentum_kg_m_s(),0.0,1e-14);
    assert_eq!(next.report().coupled.deposited_joule_j,0.0);
    assert!(next.report().coupled.source.reversible_heat_j>0.0);
    assert!(expected>300.0+0.004/CAPACITY);
}

#[test]
fn ceramic_thermal_refusals_preserve_native_and_caloric_history() {
    let card=source_card(0.2,SLOPES,0.0); let curve=caloric_curve(&card,false); let point=point(300.0);
    let thermal=initial_thermal(&curve,0.04); let cells=[cell(&card,&curve,&point)];
    let resistors=[resistor()]; let network=MaterialNetwork {cells:&cells,resistors:&resistors,links:&[]};
    let binding=thermal_binding(&card,[280.0,360.0]); let body=model(&card,&network,&thermal,0.02,0.01,0.1);
    let before=body.clone(); let before_cells=thermal.cells().clone();
    let bases=[(1,ElectricalPropertyBasis::Resistivity)];
    let mut polls=0;
    let accepted=heated::step_candidate(&network,&thermal,0.1,&body,&[],0.03,&binding,
        &bases,thermal_controls(),&[0.0],||{polls+=1;false}).unwrap();
    let mut calls=0;
    assert!(matches!(heated::step_candidate(&network,&thermal,0.1,&body,&[],0.03,&binding,
        &bases,thermal_controls(),&[0.0],||{calls+=1;calls==polls}),Err(CircuitEnthalpyError::Cancelled)));
    let mut budget=thermal_controls(); budget.maximum_circuit_trials=1;
    assert!(heated::step_candidate(&network,&thermal,0.1,&body,&[],0.03,&binding,
        &bases,budget,&[0.0],||false).is_err());
    let narrow=PiezoelectricThermalBinding { law:ThermalLongitudinalPiezoelectric::new(
        &card,&ceramic(&card),"T",[280.0,300.00001]).unwrap(),..binding };
    assert!(heated::step_candidate(&network,&thermal,0.1,&body,&[],0.03,&narrow,
        &bases,thermal_controls(),&[0.0],||false).is_err());
    let binding=thermal_binding(&card,[280.0,360.0]);
    let excessive=PiezoelectricThermalBinding {maximum_abs_displacement_c_m2:100.0,
        law:binding.law.clone(),cell:CELL,volume_tolerance_m3:1e-12};
    assert!(heated::step_candidate(&network,&thermal,0.1,&body,&[],0.03,&excessive,
        &bases,thermal_controls(),&[0.0],||false).is_err());
    // A different card cannot silently replace either accepted material history.
    let other_card=source_card(0.21,SLOPES,0.0);
    let foreign=thermal_binding(&other_card,[280.0,360.0]);
    assert!(heated::step_candidate(&network,&thermal,0.1,&body,&[],0.03,&foreign,
        &bases,thermal_controls(),&[0.0],||false).is_err());
    assert_eq!(body,before); assert_eq!(thermal.cells(),&before_cells); assert_eq!(thermal.time_s(),0.0);
    let retry=heated::step_candidate(&network,&thermal,0.1,&body,&[],0.03,&binding,
        &bases,thermal_controls(),&[0.0],||false).unwrap();
    assert_eq!(accepted.model(),retry.model());
    assert_eq!(accepted.thermal().cells(),retry.thermal().cells());
    // Unsupported density/phase response remains a refusal, not hidden mass/work.
    let density_card=source_card(0.2,SLOPES,0.001);
    assert!(ThermalLongitudinalPiezoelectric::new(&density_card,&ceramic(&density_card),
        "T",[280.0,360.0]).is_err());
    let unstable=source_card(0.2,[0.0,0.02,0.0],0.0);
    assert!(ThermalLongitudinalPiezoelectric::new(&unstable,&ceramic(&unstable),
        "T",[280.0,360.0]).is_err());
    let phase_curve=caloric_curve(&card,true); let phase_thermal=initial_thermal(&phase_curve,0.04);
    let phase_cells=[cell(&card,&phase_curve,&point)];
    let phase_network=MaterialNetwork {cells:&phase_cells,resistors:&resistors,links:&[]};
    assert!(heated::step_candidate(&phase_network,&phase_thermal,0.1,&body,&[],0.03,&binding,
        &bases,thermal_controls(),&[0.0],||false).is_err());
}
