//! G1/G3/G4 full-field temperature feedback, reversible heat and physical shunts.
use super::*;
use fs_couple::piezoelectric::spatial::transient::load::thermal::tensor::{
    ThermalTensorPiezoelectricBody, TensorThermalReport,
};
use fs_material::piezoelectric::thermal::ThermalTensorPiezoelectric;

const RHO: f64 = 7500.0;
const CP: f64 = 1.0;
const E33: f64 = 12.0;
const C_PRIME: f64 = -1e8;
const E_PRIME: f64 = 0.01;
const P_PRIME: f64 = 2e-11;
const D_INITIAL: f64 = 1e-3;

fn absolute_temperature() -> QuantitySpec {
    use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
    QuantitySpec::semantic(SemanticType::new(QuantityKind::AbsoluteTemperature, ValueForm::Static))
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new().with_quantity("T", absolute_temperature(), t).unwrap()
}
fn resistor_source() -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let mut claims=ClaimSet::new();
    for (name,dims,value) in [
        ("density",Density::DIMS,1.0),
        ("specific_heat_capacity",SPECIFIC_HEAT_CAPACITY_DIMS,2.0),
        (ELECTRICAL_RESISTIVITY_PROPERTY,ELECTRICAL_RESISTIVITY_DIMS,RESISTANCE*AREA/LENGTH),
    ] {
        claims.insert_claim(PropertyClaim {
            key:PropertyKey::new(name,dims),value:PropertyValue::Scalar{value,dims},
            validity:ValidityDomain::unconstrained().with_quantity(
                "T",absolute_temperature(),280.0,400.0),
            uncertainty:UncertaintyModel::Unstated,
            interpolation:InterpolationPolicy::ConstantWithinValidity,observations:vec![],
            provenance:Provenance{source:"authored isolated resistor".into(),
                license:"CC0-1.0".into(),artifact:None},
        }).unwrap();
    }
    let card=MaterialCard::assemble(MaterialStateId {
        chemistry:"synthetic physical shunt".into(),phase:"solid".into(),
        process:"uniform conductor".into(),revision:0,
    },claims,vec![]).unwrap();
    let curve=ResolvedHeatCapacityCurve::try_new(&card,HeatCapacitySourceQuery{
        heat_capacity:&PropertyKey::new("specific_heat_capacity",SPECIFIC_HEAT_CAPACITY_DIMS),
        density:&PropertyKey::new("density",Density::DIMS),temperature_axis:"T",
        lower:&point(280.0),upper:&point(400.0),selection:SelectionPolicy::SingleClaimOnly,
    },SolidLiquidPhase::Solid,0.0,1e-8).unwrap();
    (card,curve)
}

fn ceramic_source() -> (MaterialCard, ResolvedHeatCapacityCurve, TensorPiezoelectric) {
    let source = material([0.0, 0.0, E33, 0.0, 0.0, 0.0]);
    let mut claims = ClaimSet::new();
    for property in source.material().properties() {
        let name = property.requirement().name();
        let dims = property.requirement().dims();
        let slope = match name {
            "piezoelectric_stiffness_e_mandel_33" => C_PRIME,
            "piezoelectric_stress_charge_mandel_33" => E_PRIME,
            "piezoelectric_permittivity_s_33" => P_PRIME,
            _ => 0.0,
        };
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Curve {
                abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                knots: vec![(280.0, property.value_si()-20.0*slope),
                    (400.0, property.value_si()+100.0*slope)], dims,
            },
            validity: ValidityDomain::unconstrained().with_quantity(
                "T", absolute_temperature(), 280.0, 400.0,
            ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "authored full-tensor finite-caloric fixture".into(),
                license: "CC0-1.0".into(), artifact: None,
            },
        }).unwrap();
    }
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
        value: PropertyValue::Scalar { value: CP, dims: SPECIFIC_HEAT_CAPACITY_DIMS },
        validity: ValidityDomain::unconstrained().with_quantity(
            "T", absolute_temperature(), 280.0, 400.0,
        ),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: Provenance { source: "authored caloric fixture".into(),
            license: "CC0-1.0".into(), artifact: None },
    }).unwrap();
    let card = MaterialCard::assemble(MaterialStateId {
        chemistry: "synthetic reciprocal temperature-dependent ceramic".into(),
        phase: "solid".into(), process: "fixed reference geometry".into(), revision: 0,
    }, claims, vec![]).unwrap();
    let curve = ResolvedHeatCapacityCurve::try_new(&card, HeatCapacitySourceQuery {
        heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
        density: &PropertyKey::new("density", Density::DIMS),
        temperature_axis: "T", lower: &point(280.0), upper: &point(400.0),
        selection: SelectionPolicy::SingleClaimOnly,
    }, SolidLiquidPhase::Solid, 0.0, 1e-8).unwrap();
    let tensor = TensorPiezoelectric::resolve(&card, &point(300.0),
        MaterialPropertySelection::SingleClaimOnly, FRAME).unwrap();
    (card, curve, tensor)
}

fn caloric(curves: &[ResolvedHeatCapacityCurve], ceramic_power: f64) -> DepositionSimulation<'_> {
    let laws = curves.iter().enumerate().map(|(i, curve)|
        (MaterialId(i as u64), curve.curve() as &dyn ThermalLaw)
    ).collect::<Vec<_>>();
    let events = curves.iter().enumerate().map(|(i, _)| TimedEvent {
        time_s: 0.0, event: Event::Birth(DepositedCell {
            id: CELLS[i], segment_id: i as u64, body_id: i as u64,
            material_id: MaterialId(i as u64),
            mass_kg: VOLUME * if i == 0 { RHO } else { 1.0 },
            specific_h_j_kg: if i == 0 { 20.0*CP } else { 40.0 },
            geometry: CellGeometry {
                start_m: [0.0; 3], end_m: [LENGTH; 3], volume_m3: VOLUME,
                surface_area_m2: 6.0*AREA,
            },
            boundary: Boundary { ambient_k: 300.0, htc_w_m2_k: 0.0,
                reservoir: None, source_w: if i == 0 { ceramic_power } else { 0.0 } },
        }),
    }).collect();
    DepositionSimulation::new(events, &laws, DepositionBudget {
        max_cells: 2, max_links: 0, max_events: 2, max_substeps: 8,
    }).unwrap().advance_to(0.0, StepConfig {
        max_dt_s: 1e-6, energy_tolerance_j: 1e-16,
    }, || Ok(())).unwrap().0
}
fn cells<'a>(
    cards: &'a [MaterialCard], curves: &'a [ResolvedHeatCapacityCurve], context: &'a QueryPoint,
) -> Vec<MaterialCell<'a>> {
    cards.iter().zip(curves).enumerate().map(|(i, (card, curve))| MaterialCell {
        cell: CELLS[i], material_id: MaterialId(i as u64), curve: curve.curve(),
        point: ThermalMaterialPoint {
            material: card, context, temperature_axis: "T", liquid_fraction_axis: None,
        },
        conductivity: None,
    }).collect()
}
fn law(card: &MaterialCard, tensor: &TensorPiezoelectric) -> ThermalTensorPiezoelectric {
    ThermalTensorPiezoelectric::new(card, tensor, "T", [280.0, 400.0]).unwrap()
}

// Independent closed-form full stress-charge potential for the uniform axial
// fixture, including the exact natural-coordinate thermal derivative.
fn potential(t: f64, strain: f64, displacement: f64) -> [f64; 3] {
    let c = STIFFNESS[2] + C_PRIME*(t-300.0);
    let e = E33 + E_PRIME*(t-300.0);
    let p = EPSILON[2] + P_PRIME*(t-300.0);
    let z = displacement-e*strain;
    let field = z/p;
    let free = 0.5*VOLUME*(c*strain*strain+z*field);
    let entropy = VOLUME*(-0.5*C_PRIME*strain*strain
        + field*E_PRIME*strain + 0.5*P_PRIME*field*field);
    [free, entropy, free+t*entropy]
}
fn stress(t: f64, strain: f64, displacement: f64) -> f64 {
    let c = STIFFNESS[2]+C_PRIME*(t-300.0);
    let e = E33+E_PRIME*(t-300.0);
    c*strain-e*(displacement-e*strain)/(EPSILON[2]+P_PRIME*(t-300.0))
}
fn root(mut residual: impl FnMut(f64) -> f64) -> f64 {
    let (mut low, mut high) = (280.001, 399.999);
    assert!(residual(low)<0.0 && residual(high)>0.0);
    for _ in 0..80 {
        let mid = 0.5*(low+high);
        if residual(mid)>0.0 { high=mid; } else { low=mid; }
    }
    0.5*(low+high)
}
fn temperatures(thermal: &DepositionSimulation<'_>, curves: &[ResolvedHeatCapacityCurve]) -> Vec<f64> {
    curves.iter().enumerate().map(|(i, curve)|
        curve.curve().state_at_specific_enthalpy(thermal.cells()[&CELLS[i]].specific_h_j_kg)
            .unwrap().temperature_k()
    ).collect()
}
fn check_report(report: &TensorThermalReport) {
    near(report.total_energy_residual_j, 0.0, 1e-16, 0.0);
    near(report.coupled.total_energy_residual_j, 0.0, 1e-16, 0.0);
    near(report.coupled.electrical.loads.energy_residual_j, 0.0, 1e-16, 0.0);
    near(report.coupled.non_electrical_heat_j,
        report.external_heat_j+report.coupled.source.reversible_heat_j, 1e-16, 0.0);
}

#[test]
fn g1_g3_full_tensor_support_motion_and_reversible_heat_follow_the_same_physical_state() {
    with_cx(|cx, _| {
        for rotate in [false, true] {
            let mut geometry = cube(1);
            let (card, curve, original) = ceramic_source();
            let q = if rotate { rotation([1.0/3.0_f64.sqrt(); 3], 0.63) }
                else { [[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]] };
            let old_nodes = geometry.nodes.clone();
            for node in &mut geometry.nodes { *node=apply(&q, node); }
            let tensor = if rotate {
                original.rotated(ContentHash([32;32]), q).unwrap()
            } else { original };
            let materials = vec![tensor.clone()];
            let model = Model::assemble_owned_materials(&geometry.nodes, &geometry.tets,
                materials, tensor.frame_identity(), budget(), cx).unwrap();
            let mut strain = 1e-5;
            let initial: Vec<_> = old_nodes.iter().map(|p| apply(&q, &[0.0,0.0,strain*p[2]])).collect();
            let supports = initial.iter().enumerate().flat_map(|(i, u)|
                (0..3).map(move |a| (3*i+a,u[a]))).collect::<Vec<_>>();
            let zero = vec![[0.0;3]; geometry.nodes.len()];
            let charges = vec![0.0; geometry.nodes.len()];
            let ordinary = Body::new(&model, &supports,
                &electrodes(&geometry,0.0,Condition::Charge(-D_INITIAL*AREA)),
                &initial,&zero,&charges,field_controls(),cx).unwrap();
            let cards = [card];
            let curves = [curve];
            let context = point(300.0);
            let material_cells = cells(&cards,&curves,&context);
            let network = MaterialNetwork { cells:&material_cells,resistors:&[],links:&[] };
            let mut thermal = caloric(&curves,0.0);
            let prepared = network.prepare(&thermal,||false).unwrap();
            let mut body = ThermalTensorPiezoelectricBody::new(&ordinary,
                vec![law(&cards[0],&tensor)],CELLS[0],0.002,cx).unwrap();
            let dt=1e-6;
            let mut temperature=300.0;
            let mut rate=0.0;
            for next_rate in [20.0,20.0,-20.0,0.0] {
                let next_strain=strain+0.5*dt*(rate+next_rate);
                let before=potential(temperature,strain,D_INITIAL);
                let expected=root(|t| potential(t,next_strain,D_INITIAL)[2]-before[2]
                    +RHO*VOLUME*CP*(t-temperature)
                    -VOLUME*stress(t,0.5*(strain+next_strain),D_INITIAL)*(next_strain-strain));
                let support_velocity=old_nodes.iter().enumerate().flat_map(|(i,p)| {
                    let v=apply(&q,&[0.0,0.0,next_rate*p[2]]);
                    (0..3).map(move |a| (3*i+a,v[a]))
                }).collect::<Vec<_>>();
                let candidate=body.step_candidate(&thermal,thermal.time_s()+dt,&prepared,&[],
                    &[Drive::Voltage(0.0),Drive::Current(0.0)],&zero,&charges,
                    Some(&support_velocity),field_controls(),thermal_controls(),&[0.0],cx).unwrap();
                check_report(candidate.report());
                let actual_t=temperatures(candidate.thermal(),&curves)[0];
                near(actual_t,expected,2e-8,0.0);
                let report=&candidate.report().coupled;
                near(report.source.field_internal_energy_j[1],
                    potential(expected,next_strain,D_INITIAL)[2],2e-16,3e-8);
                near(report.electrical.native.kinetic_energy_change_j,
                    RHO*VOLUME*LENGTH*LENGTH*(next_rate*next_rate-rate*rate)/6.0,1e-16,3e-8);
                assert_eq!(report.deposited_joule_j,0.0);
                assert_eq!(report.electrical.native.material_heat_j,0.0);
                near(candidate.body().state().electrodes()[1].charge_c,-D_INITIAL*AREA,1e-18,3e-8);
                for (p,u) in old_nodes.iter().zip(candidate.body().state().displacement_m()) {
                    let expected_u=apply(&q,&[0.0,0.0,next_strain*p[2]]);
                    for a in 0..3 { near(u[a],expected_u[a],2e-15,3e-8); }
                }
                let (next_thermal,next_body,_)=candidate.into_parts();
                thermal=next_thermal;body=next_body;
                temperature=expected;strain=next_strain;rate=next_rate;
            }
            assert!(body.state().time_s()>0.0);
            assert_ne!(body.materials()[0].identity(),tensor.identity());
            assert_eq!(body.materials()[0].material().card_identity(),cards[0].content_hash());
        }
    });
}

#[test]
fn g1_heating_updates_spatial_capacitance_at_conserved_charge_and_keeps_external_heat_separate() {
    with_cx(|cx,_| {
        for refinement in [1,2] {
            let geometry=cube(refinement);
            let (card,curve,tensor)=ceramic_source();
            let materials=[tensor.clone()];
            let model=Model::assemble(&geometry.nodes,&geometry.tets,&materials,FRAME,budget(),cx).unwrap();
            let n=geometry.nodes.len();
            let zero=vec![[0.0;3];n];let charges=vec![0.0;n];
            let supports=(0..3*n).map(|i|(i,0.0)).collect::<Vec<_>>();
            let ordinary=Body::new(&model,&supports,
                &electrodes(&geometry,0.0,Condition::Charge(-D_INITIAL*AREA)),
                &zero,&zero,&charges,field_controls(),cx).unwrap();
            let cards=[card];let curves=[curve];let context=point(300.0);
            let material_cells=cells(&cards,&curves,&context);
            let network=MaterialNetwork{cells:&material_cells,resistors:&[],links:&[]};
            let power=0.015;let dt=2e-4;
            let thermal=caloric(&curves,power);
            let prepared=network.prepare(&thermal,||false).unwrap();
            let body=ThermalTensorPiezoelectricBody::new(&ordinary,
                vec![law(&cards[0],&tensor)],CELLS[0],0.002,cx).unwrap();
            let before=potential(300.0,0.0,D_INITIAL)[2];
            let expected=root(|t| potential(t,0.0,D_INITIAL)[2]-before
                +RHO*VOLUME*CP*(t-300.0)-power*dt);
            let candidate=body.step_candidate(&thermal,dt,&prepared,&[],
                &[Drive::Voltage(0.0),Drive::Current(0.0)],&zero,&charges,None,
                field_controls(),thermal_controls(),&[0.0],cx).unwrap();
            check_report(candidate.report());
            near(temperatures(candidate.thermal(),&curves)[0],expected,2e-8,0.0);
            near(candidate.report().external_heat_j,power*dt,1e-16,0.0);
            near(candidate.body().state().electrodes()[1].voltage_v,
                -LENGTH*D_INITIAL/(EPSILON[2]+P_PRIME*(expected-300.0)),2e-6,3e-8);
            assert_eq!(candidate.body().state().velocity_m_s(),zero);
            near(body.state().electrodes()[1].voltage_v,-LENGTH*D_INITIAL/EPSILON[2],2e-6,3e-8);
            assert_eq!(thermal.time_s(),0.0);
        }
    });
}

#[test]
fn g1_temperature_dependent_spatial_storage_and_physical_resistor_heat_close_together() {
    with_cx(|cx,_| {
        let geometry=cube(1);
        let (card,curve,tensor)=ceramic_source();
        let (resistor_card,resistor_curve)=resistor_source();
        let materials=[tensor.clone()];
        let model=Model::assemble(&geometry.nodes,&geometry.tets,&materials,FRAME,budget(),cx).unwrap();
        let n=geometry.nodes.len();
        let zero=vec![[0.0;3];n];let charges=vec![0.0;n];
        let supports=(0..3*n).map(|i|(i,0.0)).collect::<Vec<_>>();
        let ordinary=Body::new(&model,&supports,
            &electrodes(&geometry,0.0,Condition::Charge(-D_INITIAL*AREA)),
            &zero,&zero,&charges,field_controls(),cx).unwrap();
        let cards=[card,resistor_card];let curves=[curve,resistor_curve];let context=point(300.0);
        let material_cells=cells(&cards,&curves,&context);
        let resistors=[CellResistor {
            destination:ResistorCell{branch:LOAD.resistor_branch,cell:CELLS[1]},
            selection:MaterialPropertySelection::SingleClaimOnly,length_m:LENGTH,area_m2:AREA,
        }];
        let network=MaterialNetwork{cells:&material_cells,resistors:&resistors,links:&[]};
        let thermal=caloric(&curves,0.0);
        let prepared=network.prepare(&thermal,||false).unwrap();
        let body=ThermalTensorPiezoelectricBody::new(&ordinary,
            vec![law(&cards[0],&tensor)],CELLS[0],0.002,cx).unwrap();
        let dt=1e-7;
        let q0=-D_INITIAL*AREA;
        let electrical=|t| {
            let c=(EPSILON[2]+P_PRIME*(t-300.0))*LENGTH;
            let v=q0/(c+dt/(2.0*RESISTANCE));
            let q=q0-dt*v/RESISTANCE;
            [q,dt*v*v/RESISTANCE]
        };
        let before=potential(300.0,0.0,D_INITIAL)[2];
        let expected=root(|t| {
            let [q,heat]=electrical(t);
            potential(t,0.0,-q/AREA)[2]-before+RHO*VOLUME*CP*(t-300.0)+heat
        });
        let [q,heat]=electrical(expected);
        let candidate=body.step_candidate(&thermal,dt,&prepared,&[LOAD],
            &[Drive::Voltage(0.0),Drive::Current(0.0)],&zero,&charges,None,
            field_controls(),thermal_controls(),&[0.0,0.0],cx).unwrap();
        check_report(candidate.report());
        let t=temperatures(candidate.thermal(),&curves);
        near(t[0],expected,2e-8,0.0);
        near(t[1],300.0+heat/(2.0*VOLUME),2e-7,0.0);
        near(candidate.body().state().electrodes()[1].charge_c,q,1e-18,3e-8);
        near(candidate.report().coupled.deposited_joule_j,heat,1e-16,3e-8);
        near(candidate.report().external_heat_j,0.0,1e-16,0.0);
        assert!(candidate.report().coupled.source.reversible_heat_j.abs()>1e-13);
        assert_eq!(candidate.report().coupled.electrical.destinations[0].cell,CELLS[1]);
    });
}

#[test]
fn g0_g4_tensor_thermal_refuses_unpaid_clamps_bad_capacity_and_interrupted_publication() {
    with_cx(|cx,gate| {
        let geometry=cube(1);let (card,curve,tensor)=ceramic_source();
        let materials=[tensor.clone()];
        let model=Model::assemble(&geometry.nodes,&geometry.tets,&materials,FRAME,budget(),cx).unwrap();
        let n=geometry.nodes.len();let zero=vec![[0.0;3];n];let charges=vec![0.0;n];
        let supports=(0..3*n).map(|i|(i,0.0)).collect::<Vec<_>>();
        let clamped=Body::new(&model,&supports,
            &electrodes(&geometry,0.0,Condition::Voltage(-100.0)),
            &zero,&zero,&charges,field_controls(),cx).unwrap();
        assert!(ThermalTensorPiezoelectricBody::new(&clamped,
            vec![law(&card,&tensor)],CELLS[0],0.002,cx).is_err());
        let ordinary=Body::new(&model,&supports,
            &electrodes(&geometry,0.0,Condition::Charge(-D_INITIAL*AREA)),
            &zero,&zero,&charges,field_controls(),cx).unwrap();
        let cards=[card];let curves=[curve];let context=point(300.0);
        let material_cells=cells(&cards,&curves,&context);
        let network=MaterialNetwork{cells:&material_cells,resistors:&[],links:&[]};
        let thermal=caloric(&curves,0.015);let prepared=network.prepare(&thermal,||false).unwrap();
        let bad=ThermalTensorPiezoelectricBody::new(&ordinary,
            vec![law(&cards[0],&tensor)],CELLS[0],1e3,cx).unwrap();
        let drives=[Drive::Voltage(0.0),Drive::Current(0.0)];
        assert!(bad.step_candidate(&thermal,1e-4,&prepared,&[],&drives,&zero,&charges,None,
            field_controls(),thermal_controls(),&[0.0],cx).is_err());
        let body=ThermalTensorPiezoelectricBody::new(&ordinary,
            vec![law(&cards[0],&tensor)],CELLS[0],0.002,cx).unwrap();
        let mut exhausted=thermal_controls();exhausted.maximum_circuit_trials=1;
        assert!(body.step_candidate(&thermal,1e-4,&prepared,&[],&drives,&zero,&charges,None,
            field_controls(),exhausted,&[0.0],cx).is_err());
        let accepted=body.step_candidate(&thermal,1e-4,&prepared,&[],&drives,&zero,&charges,None,
            field_controls(),thermal_controls(),&[0.0],cx).unwrap();
        let retry=body.step_candidate(&thermal,1e-4,&prepared,&[],&drives,&zero,&charges,None,
            field_controls(),thermal_controls(),&[0.0],cx).unwrap();
        assert_eq!(accepted.body().state().voltage_v(),retry.body().state().voltage_v());
        assert_eq!(accepted.thermal().cells()[&CELLS[0]].specific_h_j_kg,
            retry.thermal().cells()[&CELLS[0]].specific_h_j_kg);
        assert_eq!(body.state().voltage_v(),ordinary.state().voltage_v());
        assert_eq!(thermal.time_s(),0.0);
        gate.request();
        assert!(matches!(body.step_candidate(&thermal,1e-4,&prepared,&[],&drives,&zero,&charges,None,
            field_controls(),thermal_controls(),&[0.0],cx),Err(CircuitEnthalpyError::Cancelled)));
        assert_eq!(body.state().voltage_v(),ordinary.state().voltage_v());
    });
}
