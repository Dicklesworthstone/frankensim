//! Material changes drive the actual reciprocal field/body/support descriptor.
//! Synthetic source coefficients test causality and accounting, not alloy fidelity.
#[path = "electrical_contacts/source.rs"]
mod source;
#[path = "circuit_enthalpy/common.rs"]
mod thermal;

use fs_couple::material_circuit::motion::{FieldMotionInput, MaterialFieldMotion};
use fs_couple::material_circuit::motion::support::AxialSupport;
use fs_couple::material_circuit::storage::{MaterialStorage, MaterialStorageElement};
use fs_couple::material_circuit::{MaterialCircuit, MaterialCircuitError, MaterialCircuitSolve};
use fs_evidence::ValidityDomain;
use fs_ga::Vec3;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::field_storage::*;
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::section::thermal::ConcentricThermoelasticSection;
use fs_material::state_point::{INVERSE_TEMPERATURE_DIMS, MaterialPropertySelection,
    ScalarAdmissibility, ScalarPropertyRequirement, integrate_isotropic_thermal_expansion,
    resolve_material_state_point};
use fs_mbd_material::{RegionBinding, ResolvedBody, VolumeMoments, resolve_body};
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_qty::{Density, Pressure, Temperature};

const GAP: f64 = 1.01;
const LENGTH: f64 = 2.0;
const REFERENCE: f64 = 1.0;

fn section(outer_e: f64, outer_rho: f64, temperature: f64, variable_e: bool)
    -> ConcentricThermoelasticSection
{
    let mut layers = Vec::new();
    let mut paths = Vec::new();
    for (index, (area, rho, young, alpha)) in [(0.01_f64, 3.0, 1000.0, 1e-4),
        (0.03_f64, outer_rho, outer_e, 2e-4)].into_iter().enumerate()
    {
        let mut claims = ClaimSet::new();
        for (key, value, dims) in [("density", rho, Density::DIMS),
            ("young_modulus", young, Pressure::DIMS),
            ("linear_thermal_expansion_coefficient", alpha, INVERSE_TEMPERATURE_DIMS)]
        {
            let variable = variable_e && key == "young_modulus";
            claims.insert_claim(PropertyClaim { key: PropertyKey::new(key, dims),
                value: if variable { PropertyValue::Curve { abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS, knots: vec![(300.0,value),(400.0,2.0*value)], dims } }
                    else { PropertyValue::Scalar { value, dims } },
                validity: ValidityDomain::unconstrained().with("T",300.0,400.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: if variable { InterpolationPolicy::LinearInside }
                    else { InterpolationPolicy::ConstantWithinValidity }, observations: vec![],
                provenance: Provenance { source: "synthetic layered actuator support".into(),
                    license: "test-only".into(), artifact: None } }).unwrap();
        }
        let card = MaterialCard::assemble(MaterialStateId { chemistry: format!("layer-{index}"),
            phase: "solid".into(), process: "synthetic".into(), revision: 0 }, claims, vec![]).unwrap();
        let point = QueryPoint::new().with("T",temperature).unwrap();
        let material = resolve_material_state_point(&card, &point, &[
            ScalarPropertyRequirement::try_new("density",Density::DIMS,ScalarAdmissibility::StrictlyPositive).unwrap(),
            ScalarPropertyRequirement::try_new("young_modulus",Pressure::DIMS,ScalarAdmissibility::StrictlyPositive).unwrap(),
        ], MaterialPropertySelection::SingleClaimOnly).unwrap();
        paths.push(integrate_isotropic_thermal_expansion(&card,
            &QueryPoint::new().with("T",300.0).unwrap(), &point,
            MaterialPropertySelection::SingleClaimOnly).unwrap());
        layers.push(ConcentricMaterialLayer {
            outer_radius_m: (area/core::f64::consts::PI).sqrt(), material });
    }
    ConcentricThermoelasticSection::try_new(
        ConcentricElasticSection::try_new(0.0,layers).unwrap(),paths).unwrap()
}
fn support(young: f64, density: f64, t: f64) -> AxialSupport {
    AxialSupport::try_new(section(young,density,t,false),LENGTH,REFERENCE,0.05).unwrap()
}
fn body() -> ResolvedBody {
    let card = source::scalar_card("density",Density::DIMS,1.0);
    resolve_body(&[RegionBinding { key: "guided body", material: &card, point: &QueryPoint::new(),
        selection: MaterialPropertySelection::SingleClaimOnly,
        geometry: VolumeMoments::cuboid(Vec3::new(1.0,1.0,1.0)).unwrap() }],1,||false).unwrap()
}
fn field(magnetic: bool) -> MaterialStorage {
    let element = if magnetic {
        let card = source::scalar_card(MAGNETIC_PERMEABILITY_PROPERTY,MAGNETIC_PERMEABILITY_DIMS,1.0);
        MaterialStorageElement::Inductor(ResolvedInductor::resolve(&card,&QueryPoint::new(),
            MaterialPropertySelection::SingleClaimOnly,1.0,GAP,1).unwrap())
    } else {
        let card = source::scalar_card(ELECTRIC_PERMITTIVITY_PROPERTY,ELECTRIC_PERMITTIVITY_DIMS,1.0);
        MaterialStorageElement::Capacitor(ResolvedCapacitor::resolve(&card,&QueryPoint::new(),
            MaterialPropertySelection::SingleClaimOnly,1.0,GAP).unwrap())
    };
    MaterialStorage { branch: 0, element }
}
fn isolated(magnetic: bool, member: AxialSupport, velocity: f64) -> MaterialFieldMotion {
    let (branches,fluxes,charges,ports) = if magnetic {
        (vec![(1,0,Branch::Inductor { henries: f64::NAN }),
            (1,0,Branch::VoltageSource { port: 0 })],vec![0.2],vec![],vec![0.0])
    } else { (vec![(1,0,Branch::Capacitor { farads: f64::NAN })],vec![],vec![0.2],vec![]) };
    let circuit = MaterialCircuit::new_with_storage(CircuitGraph { node_count: 2, branches,
        transformers: vec![] },vec![],vec![field(magnetic)],fluxes,charges,&ports).unwrap();
    MaterialFieldMotion::new_with_support(circuit,0,body(),member,velocity,[0.5,1.5]).unwrap()
}
fn controls() -> MaterialCircuitSolve {
    MaterialCircuitSolve { energy_tolerance_j: 1e-9, power_tolerance_w: 1e-9 }
}
fn near(a: f64,b: f64) { assert!((a-b).abs()<3e-9,"{a} != {b}"); }

#[test]
fn both_field_types_follow_the_same_layered_mass_stiffness_and_preload() {
    let (dt,velocity,external) = (0.1,0.025,0.02);
    for magnetic in [false,true] {
        for (e,rho,t) in [(2000.0,6.0,300.0),(8000.0,6.0,300.0),
            (2000.0,24.0,300.0),(2000.0,6.0,320.0)]
        {
            let member = support(e,rho,t);
            let mass = 1.0+(3.0*0.01+rho*0.02)*LENGTH/3.0;
            let k = (1000.0*0.01+e*0.02)/LENGTH;
            let preload = (1000.0*0.01*1e-4+e*0.02*2e-4)*(t-300.0);
            let force = external-0.5*0.2*0.2+preload;
            let a = dt*dt*k/(4.0*mass);
            let v1 = ((1.0-a)*velocity+dt*(force-k*(GAP-REFERENCE))/mass)/(1.0+a);
            let gap1 = GAP+dt*(velocity+v1)/2.0;
            let model = isolated(magnetic,member.clone(),velocity);
            let before = model.clone();
            let ports = if magnetic { vec![0.0] } else { vec![] };
            let next = model.step_candidate(vec![],FieldMotionInput { electrical_ports: &ports,
                external_force_n: external, duration_s: dt },controls(),||false).unwrap();
            near(model.generalized_mass_kg(),mass);
            near(next.model().gap_m(),gap1);
            near(next.model().velocity_m_s(),v1);
            near(next.step().kinetic_energy_j[1],0.5*mass*v1*v1);
            let elastic = [(1000.0,0.01,1e-4),(e,0.02,2e-4)].into_iter()
                .map(|(young,area,alpha)| 0.5*young*area*LENGTH
                    *((gap1-REFERENCE)/LENGTH-alpha*(t-300.0)).powi(2)).sum::<f64>();
            near(next.step().elastic_energy_j[1],elastic);
            near(next.step().field_work_on_body_j,-0.5*0.2*0.2*(gap1-GAP));
            near(next.step().energy_residual_j,0.0);
            assert_eq!(next.step().heat_j,0.0);
            assert_eq!(next.model().support(),Some(&member));
            assert_eq!(model,before);
        }
    }
}

#[test]
fn changing_support_data_changes_charge_motion_and_finite_joule_heating() {
    use fs_couple::circuit_enthalpy::motion::{self,FieldMotionEnthalpyInput};
    let run = |young,density| {
        let (a,ca) = thermal::material(1.0,1.0,0.005);
        let (b,cb) = thermal::material(2.0,2.0,0.0);
        let cards = [a,b]; let curves = [ca,cb];
        let mut heat = thermal::initial(&curves,0.3,None);
        let circuit = MaterialCircuit::new_with_storage(CircuitGraph { node_count: 3, branches: vec![
            (1,0,Branch::Capacitor { farads: f64::NAN }),
            (1,2,Branch::Resistor { ohms: f64::NAN }),
            (2,0,Branch::Resistor { ohms: f64::NAN })],transformers:vec![] },
            thermal::coefficients(&heat,&cards,&curves).unwrap().resistors,
            vec![field(false)],vec![],vec![0.6],&[]).unwrap();
        let mut moving = MaterialFieldMotion::new_with_support(circuit,0,body(),
            support(young,density,300.0),0.02,[0.5,1.5]).unwrap();
        let storage = |m: &MaterialFieldMotion| 0.5*m.circuit().charges_c()[0].powi(2)*m.gap_m()
            +0.5*m.generalized_mass_kg()*m.velocity_m_s().powi(2)+m.elastic_energy_j().unwrap();
        let initial = heat.energy_j()+storage(&moving);
        for index in 1..=8 {
            let next = motion::step_candidate(&heat,f64::from(index)*0.025,
                FieldMotionEnthalpyInput { motion:&moving,ports:&[],external_force_n:0.0,
                    destinations:&thermal::DESTINATIONS },thermal::config(),&[0.0;2],
                |state| thermal::coefficients(state,&cards,&curves),||false).unwrap();
            let report = next.report();
            near(report.electrical.elastic_energy_j[1],next.motion().elastic_energy_j().unwrap());
            thermal::near(next.thermal().energy_j()+storage(next.motion()),initial);
            thermal::near(report.total_energy_residual_j,0.0);
            thermal::near(report.deposited_joule_j,report.electrical.heat_j);
            let (next_heat,next_motion,_) = next.into_parts();
            heat=next_heat; moving=next_motion;
        }
        (moving.gap_m(),moving.circuit().charges_c()[0],heat.energy_j())
    };
    let original = run(2000.0,6.0);
    for changed in [run(8000.0,6.0),run(2000.0,24.0)] {
        assert!((changed.0-original.0).abs()>1e-5);
        assert!((changed.1-original.1).abs()>1e-7,"support must change the actual circuit history");
        assert!((changed.2-original.2).abs()>1e-8,"changed circuit work must reach thermal stores");
    }
}

#[test]
fn cancellation_and_layer_strain_refusal_preserve_the_whole_assembly() {
    assert!(AxialSupport::try_new(section(2000.0,6.0,320.0,true),LENGTH,REFERENCE,0.05).is_err());
    assert!(AxialSupport::try_new(section(2000.0,6.0,320.0,false),LENGTH,REFERENCE,1e-4).is_err());
    let member = support(2000.0,6.0,300.0);
    assert!(member.force_and_energy(1.2).is_err());
    let model = isolated(false,member,0.025);
    let before = model.clone();
    let input = FieldMotionInput { electrical_ports:&[],external_force_n:1e6,duration_s:1.0 };
    assert!(model.step_candidate(vec![],input,controls(),||false).is_err());
    assert_eq!(model,before);
    assert!(matches!(model.step_candidate(vec![],input,controls(),||true),
        Err(MaterialCircuitError::Cancelled)));
    assert_eq!(model,before);
}
