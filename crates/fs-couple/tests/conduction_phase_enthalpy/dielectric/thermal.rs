//! Independent thermodynamic and native spatial electrical/thermal references.
//! These constitutive numbers are synthetic, not qualified dielectric data.
use super::*;
use fs_material::field_storage::{
    planar::{PlanarDielectricLaw, PRINCIPAL_ELECTRIC_PERMITTIVITY_PROPERTIES},
    thermal::electric::ThermalPlanarDielectricMaterial,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

const CAPACITY: f64 = 4.0;
fn absolute_temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(QuantityKind::AbsoluteTemperature, ValueForm::Static))
}
fn context(t: f64) -> QueryPoint {
    QueryPoint::new().with_quantity("T", absolute_temperature(), t).unwrap()
}
fn source(slope: f64, resistance: f64) -> MaterialCard {
    let basis = ElectricalPropertyBasis::Resistivity;
    let mut claims = ClaimSet::new();
    for (name, dims, value, derivative) in [
        (basis.property_name(), basis.dimensions(), resistance, 0.01 * resistance),
        (ELECTRIC_PERMITTIVITY_PROPERTY, ELECTRIC_PERMITTIVITY_DIMS, 2.0, slope),
        (PRINCIPAL_ELECTRIC_PERMITTIVITY_PROPERTIES[0], ELECTRIC_PERMITTIVITY_DIMS, 2.0, slope),
        (PRINCIPAL_ELECTRIC_PERMITTIVITY_PROPERTIES[1], ELECTRIC_PERMITTIVITY_DIMS, 4.0, -0.5 * slope),
    ] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                knots: vec![(280.0, value - 20.0 * derivative), (360.0, value + 60.0 * derivative)], dims },
            validity: ValidityDomain::unconstrained().with_quantity("T", absolute_temperature(), 280.0, 360.0),
            interpolation: InterpolationPolicy::LinearInside, uncertainty: UncertaintyModel::Unstated,
            observations: vec![], provenance: Provenance { source: "authored affine electrocaloric fixture".into(),
                license: "test-only".into(), artifact: None },
        }).unwrap();
    }
    MaterialCard::assemble(MaterialStateId { chemistry: "synthetic thermodielectric".into(),
        phase: "solid".into(), process: "numerical fixture".into(), revision: 0 }, claims, vec![]).unwrap()
}
fn source_conductor(card: &MaterialCard, point: &QueryPoint) -> ResolvedConductorState {
    let basis = ElectricalPropertyBasis::Resistivity;
    let mut requirements = vec![ScalarPropertyRequirement::try_new(
        basis.property_name(), basis.dimensions(), ScalarAdmissibility::StrictlyPositive,
    ).unwrap()];
    for name in [ELECTRIC_PERMITTIVITY_PROPERTY, PRINCIPAL_ELECTRIC_PERMITTIVITY_PROPERTIES[0],
        PRINCIPAL_ELECTRIC_PERMITTIVITY_PROPERTIES[1]] {
        requirements.push(ScalarPropertyRequirement::try_new(name, ELECTRIC_PERMITTIVITY_DIMS,
            ScalarAdmissibility::StrictlyPositive).unwrap());
    }
    ResolvedConductorState::try_from_material_with_basis(resolve_material_state_point(
        card, point, &requirements, MaterialPropertySelection::SingleClaimOnly,
    ).unwrap(), basis).unwrap()
}
fn source_law(card: &MaterialCard, angle: Option<f64>, range: [f64; 2]) -> ThermalPlanarDielectricMaterial {
    let native = source_conductor(card, &context(300.0));
    let law = angle.map_or(PlanarDielectricLaw::Isotropic, |angle| {
        PlanarDielectricLaw::Orthotropic(PlanarElectricalFrame::new(
            ContentHash([19;32]), ContentHash([2;32]), angle,
        ).unwrap())
    });
    ThermalPlanarDielectricMaterial::new(card, native.material(), law, "T", range).unwrap()
}
fn caloric(card: &MaterialCard, phase: bool) -> EquilibriumEnthalpyPhaseCurve {
    let knots = vec![
        EnthalpyPhaseKnot { specific_enthalpy_j_kg: 0.0, temperature_k: 280.0,
            liquid_mass_fraction: 0.0, bulk_density_kg_m3: 1.0 },
        EnthalpyPhaseKnot { specific_enthalpy_j_kg: 80.0 * CAPACITY, temperature_k: 360.0,
            liquid_mass_fraction: if phase {1.0} else {0.0}, bulk_density_kg_m3: 1.0 },
    ];
    if phase {
        EquilibriumEnthalpyPhaseCurve::try_new(card.content_hash(), knots).unwrap()
    } else {
        EquilibriumEnthalpyPhaseCurve::try_single_phase(card.content_hash(),
            fs_material::phase::SolidLiquidPhase::Solid, knots).unwrap()
    }
}
fn initial_thermal(curve: &EquilibriumEnthalpyPhaseCurve, power: f64) -> DepositionSimulation<'_> {
    DepositionSimulation::new(vec![TimedEvent { time_s: 0.0, event: Event::Birth(DepositedCell {
        id: CellId(7), segment_id: 0, body_id: 0, material_id: MaterialId(3), mass_kg: 1.0,
        specific_h_j_kg: 20.0 * CAPACITY,
        geometry: CellGeometry { start_m: [0.0;3], end_m: [1.0,0.0,0.0], volume_m3: 1.0, surface_area_m2: 6.0 },
        boundary: Boundary { ambient_k: 300.0, htc_w_m2_k: 0.0, reservoir: None, source_w: power },
    }) }], &[(MaterialId(3), curve as &dyn ThermalLaw)], DepositionBudget {
        max_cells: 1, max_links: 0, max_events: 1, max_substeps: 100,
    }).unwrap().advance_to(0.0, StepConfig { max_dt_s: 0.1, energy_tolerance_j: 1e-8 }, || Ok(())).unwrap().0
}
fn cell<'a>(card: &'a MaterialCard, curve: &'a EquilibriumEnthalpyPhaseCurve,
    point: &'a QueryPoint) -> MaterialCell<'a> {
    MaterialCell { cell: CellId(7), material_id: MaterialId(3), curve,
        point: ThermalMaterialPoint { material: card, context: point,
            temperature_axis: "T", liquid_fraction_axis: None }, conductivity: None }
}
fn thermal_field(card: &MaterialCard, angle: Option<f64>, driven: bool, range: [f64; 2])
    -> PlanarLeakyDielectric<ResolvedConductorState> {
    let conduction = PlanarConduction::new(mesh(), vec![source_conductor(card, &context(300.0));4],
        limits(), &mut || false).unwrap();
    PlanarLeakyDielectric::new_with_thermal_dielectric_laws(conduction,
        vec![source_law(card, angle, range);4], vec![0.0,4.0,4.0,0.0,2.0],
        if driven { &[0,1,2,3] } else { &[0,3] }, 0.0,
        limits(), electric_controls(), &mut || false).unwrap()
}
fn epsilon(t: f64, e: f64, slope: f64) -> f64 { e + slope * (t - 300.0) }
fn field_u(t: f64, voltage: f64, e: f64, slope: f64) -> f64 {
    0.5 * (epsilon(t,e,slope) + t * slope) * voltage * voltage
}
fn bisect(mut low: f64, mut high: f64, residual: impl Fn(f64) -> f64) -> f64 {
    assert!(residual(low) < 0.0 && residual(high) > 0.0);
    for _ in 0..100 { let mid = 0.5 * (low + high);
        if residual(mid) < 0.0 { low = mid; } else { high = mid; }
    }
    0.5 * (low + high)
}
// Independent scalar elimination for the unit-volume affine field. No runtime
// matrices, source evaluations or thermal-state helpers define this reference.
fn endpoint(t0: f64, v0: f64, e: f64, slope: f64, resistance: f64,
    dt: f64, external_heat: f64, driven: Option<f64>) -> (f64,f64,f64,f64) {
    let sigma = |t: f64| 1.0 / (resistance * (1.0 + 0.01 * (t - 300.0)));
    let voltage = |t: f64| driven.unwrap_or_else(|| {
        v0 * (epsilon(t0,e,slope) - 0.5 * dt * sigma(t))
            / (epsilon(t,e,slope) + 0.5 * dt * sigma(t))
    });
    let work = |t: f64| if driven.is_some() {
        let v1 = voltage(t); let middle = 0.5 * (v0+v1);
        middle * (epsilon(t,e,slope)*v1-epsilon(t0,e,slope)*v0 + dt*sigma(t)*middle)
    } else {0.0};
    let t = bisect(280.0,360.0,|t| CAPACITY*(t-t0)
        + field_u(t,voltage(t),e,slope) - field_u(t0,v0,e,slope) - work(t) - external_heat);
    let v = voltage(t);
    (t,v,dt*sigma(t)*(0.5*(v0+v)).powi(2),work(t))
}

#[test]
fn oriented_dielectric_thermodynamics_uses_fixed_displacement_source_derivatives() {
    let card = source(0.004,1.0);
    let angle: f64 = 0.37;
    let law = source_law(&card,Some(angle),[280.0,360.0]);
    let (t,volume) = (323.0,0.7);
    let d = [3.0,-1.0];
    let (s,c) = angle.sin_cos();
    let principal = [c*d[0]+s*d[1],-s*d[0]+c*d[1]];
    let mut expected = [0.0;4];
    for (i,(e,b)) in [(2.0,0.004),(4.0,-0.002)].into_iter().enumerate() {
        let eps = epsilon(t,e,b);
        let psi = volume*principal[i]*principal[i]/(2.0*eps);
        let entropy = volume*principal[i]*principal[i]*b/(2.0*eps*eps);
        expected[0] += psi; expected[1] += entropy; expected[2] += psi+t*entropy;
        expected[3] -= t*volume*principal[i]*principal[i]*b*b/(eps*eps*eps);
    }
    let at = law.at_temperature(t).unwrap();
    let actual = at.state(d,volume).unwrap();
    for (a,b) in [actual.free_energy_j,actual.entropy_j_k,actual.internal_energy_j,actual.heat_capacity_j_k]
        .into_iter().zip(expected) { assert!((a-b).abs()<1e-12,"{a} != {b}"); }
    let low = law.at_temperature(t-0.001).unwrap().state(d,volume).unwrap();
    let high = law.at_temperature(t+0.001).unwrap().state(d,volume).unwrap();
    assert!((-(high.free_energy_j-low.free_energy_j)/0.002-actual.entropy_j_k).abs()<1e-9);
    assert!(((high.internal_energy_j-low.internal_energy_j)/0.002-actual.heat_capacity_j_k).abs()<1e-9);
    assert!(actual.heat_capacity_j_k<0.0);
    assert!(actual.heat_capacity_j_k>=law.minimum_heat_capacity_j_k(d[0].hypot(d[1]),volume).unwrap());
    assert_eq!(at.material().card_identity(),card.content_hash());
}

#[test]
fn native_thermal_dielectric_closes_charge_joule_reversible_heat_and_electrode_work() {
    for slope in [-0.004,0.0,0.004] {
        for (angle,e,b,driven) in [(None,2.0,slope,false),
            (Some(core::f64::consts::FRAC_PI_4),2.0,slope,true),
            (Some(3.0*core::f64::consts::FRAC_PI_4),4.0,-0.5*slope,false)] {
            let card = source(slope,1.0); let curve = caloric(&card,false);
            let point = context(300.0); let thermal = initial_thermal(&curve,0.4);
            let cells = [cell(&card,&curve,&point)]; let layout = network(&cells);
            let original = thermal_field(&card,angle,driven,[280.0,360.0]);
            let mut old = original.clone(); let mut heat = thermal;
            let (mut t0,mut v0) = (300.0,4.0);
            let initial_total = 20.0*CAPACITY+field_u(t0,v0,e,b);
            let mut work_total = 0.0;
            for step in 1..=2 {
                let voltage = if driven {Some(4.0+0.5*f64::from(step))} else {None};
                let boundary = if let Some(v) = voltage {vec![(0,0.0),(3,0.0),(1,v),(2,v)]}
                    else {vec![(0,0.0),(3,0.0)]};
                let expected = endpoint(t0,v0,e,b,1.0,0.5,0.2,voltage);
                let candidate = layout.step_dielectric_thermal_candidate(&heat,0.5*f64::from(step),
                    &old,&boundary,&[24.0;4],limits(),electric_controls(),config(),&[2.0],||false).unwrap();
                near(candidate.thermal().temperature_k(CellId(7)).unwrap(),expected.0);
                near(candidate.field().voltage_v()[1],expected.1);
                let report = candidate.report(); let native = &report.coupled.electrical;
                near(native.heat_j,expected.2); near(native.supplied_work_j,expected.3);
                near(report.external_heat_j,0.2);
                near(report.coupled.source.field_internal_energy_j[1],field_u(expected.0,expected.1,e,b));
                near(report.coupled.thermal.stored_energy_change_j,
                    native.heat_j+report.coupled.source.reversible_heat_j+0.2);
                let p = -0.5*(epsilon(expected.0,e,b)-epsilon(t0,e,b))*v0*expected.1;
                near(native.parameter_work_j,p);
                near(native.stored_after_j-native.stored_before_j+native.heat_j-native.supplied_work_j,p);
                for field in &report.coupled.source.fields {
                    near(field.reversible_heat_j,-field.parameter_work_j
                        -(field.temperature_k[1]*field.state[1].entropy_j_k-field.temperature_k[0]*field.state[0].entropy_j_k));
                }
                work_total += native.supplied_work_j;
                near(candidate.thermal().energy_j()+report.coupled.source.field_internal_energy_j[1],
                    initial_total+work_total+0.2*f64::from(step));
                if slope==0.0 { assert_eq!(native.parameter_work_j,0.0);
                    assert_eq!(report.coupled.source.reversible_heat_j,0.0); }
                for source in candidate.field().conduction().conductors() {
                    near(source.material().query().axes()["T"],expected.0);
                    assert_eq!(source.material().card_identity(),card.content_hash());
                }
                (t0,v0)=(expected.0,expected.1);
                let (next_heat,next_field,_) = candidate.into_parts(); heat=next_heat; old=next_field;
            }
            near(original.voltage_v()[1],4.0); assert_eq!(original.time_s(),0.0);
        }
    }
}

#[test]
fn isolated_charge_limit_and_source_capacity_phase_and_late_cancel_refuse_atomically() {
    let slope = 0.004; let card = source(slope,1e12); let curve = caloric(&card,false);
    let point = context(300.0); let thermal = initial_thermal(&curve,4.0);
    let cells = [cell(&card,&curve,&point)]; let layout = network(&cells);
    let original = thermal_field(&card,None,false,[280.0,360.0]);
    let charge = original.nodal_charge_c().to_vec();
    // The native model has positive finite conductivity. Its 10^12-ohm m limit
    // differs negligibly from independently conserved D=epsilon0*V0=8 C/m2.
    let isolated_t = bisect(300.0,310.0,|t| CAPACITY*(t-300.0)
        + field_u(t,8.0/epsilon(t,2.0,slope),2.0,slope)-field_u(300.0,4.0,2.0,slope)-2.0);
    let mut polls=0;
    let accepted=layout.step_dielectric_thermal_candidate(&thermal,0.5,&original,
        &[(0,0.0),(3,0.0)],&[10.0;4],limits(),electric_controls(),config(),&[2.0],
        ||{polls+=1;false}).unwrap();
    near(accepted.thermal().temperature_k(CellId(7)).unwrap(),isolated_t);
    near(accepted.field().voltage_v()[1],8.0/epsilon(isolated_t,2.0,slope));
    assert!(accepted.report().coupled.electrical.heat_j<1e-9);
    assert!(accepted.report().coupled.source.reversible_heat_j>0.0);
    assert!(isolated_t>300.5);
    let mut calls=0;
    assert!(matches!(layout.step_dielectric_thermal_candidate(&thermal,0.5,&original,
        &[(0,0.0),(3,0.0)],&[10.0;4],limits(),electric_controls(),config(),&[2.0],
        ||{calls+=1;calls>=polls}),Err(CircuitEnthalpyError::Cancelled)));
    assert!(layout.step_dielectric_thermal_candidate(&thermal,0.5,&original,
        &[(0,0.0),(3,0.0)],&[1e4;4],limits(),electric_controls(),config(),&[2.0],||false).is_err());
    let restricted=thermal_field(&card,None,false,[280.0,300.1]);
    assert!(layout.step_dielectric_thermal_candidate(&thermal,0.5,&restricted,
        &[(0,0.0),(3,0.0)],&[10.0;4],limits(),electric_controls(),config(),&[2.0],||false).is_err());
    assert!(layout.step_dielectric_candidate(&thermal,0.5,&original,
        &[(0,0.0),(3,0.0)],limits(),electric_controls(),config(),&[2.0],||false).is_err());
    let mut budget=config(); budget.maximum_circuit_trials=1;
    assert!(layout.step_dielectric_thermal_candidate(&thermal,0.5,&original,
        &[(0,0.0),(3,0.0)],&[10.0;4],limits(),electric_controls(),budget,&[2.0],||false).is_err());
    let phase=caloric(&card,true); let phase_thermal=initial_thermal(&phase,4.0);
    let phase_cells=[cell(&card,&phase,&point)];
    assert!(network(&phase_cells).step_dielectric_thermal_candidate(&phase_thermal,0.5,&original,
        &[(0,0.0),(3,0.0)],&[10.0;4],limits(),electric_controls(),config(),&[2.0],||false).is_err());
    assert_eq!(original.nodal_charge_c(),charge.as_slice()); assert_eq!(original.time_s(),0.0);
    near(thermal.energy_j(),80.0); near(original.stored_energy_j(),16.0);
    // A retry starts from the same charge/thermal state and reproduces acceptance.
    let retry=layout.step_dielectric_thermal_candidate(&thermal,0.5,&original,
        &[(0,0.0),(3,0.0)],&[10.0;4],limits(),electric_controls(),config(),&[2.0],||false).unwrap();
    assert_eq!(retry.field().voltage_v(),accepted.field().voltage_v());
    assert_eq!(retry.thermal().energy_j().to_bits(),accepted.thermal().energy_j().to_bits());
}
