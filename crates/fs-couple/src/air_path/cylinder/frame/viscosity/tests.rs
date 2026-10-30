//! Native solver/source regressions. All claims are synthetic, not material measurements.
use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::state_point::{MaterialPropertySelection, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point, SPECIFIC_HEAT_CAPACITY_DIMS,
    THERMAL_CONDUCTIVITY_DIMS};
use fs_qty::{Density, DynViscosity, Pressure, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm}};
use fs_solid::frame2d::{FrameMassFormulation, FrameViscousSection, Geometry};

pub(crate) fn point(t: f64) -> QueryPoint {
    QueryPoint::new().with("T", t).unwrap().with("omega", 1.0).unwrap()
}
fn claim(key: PropertyKey, value: PropertyValue, viscous_solid: bool) -> PropertyClaim {
    let mut validity = ValidityDomain::unconstrained().with("T", 290.0, 500.0);
    if viscous_solid { validity = validity.with("omega", 0.0, 1e7); }
    let interpolation = match &value {
        PropertyValue::Scalar {..} => InterpolationPolicy::ConstantWithinValidity,
        _ => InterpolationPolicy::LinearInside,
    };
    PropertyClaim {key, value, validity, interpolation, uncertainty: UncertaintyModel::Unstated,
        observations: vec![], provenance: Provenance {source: "manufactured frame thermal test".into(),
            license: "CC0-1.0".into(), artifact: None}}
}
fn scalar(name: &str, dims: fs_qty::Dims, value: f64) -> PropertyClaim {
    claim(PropertyKey::new(name, dims), PropertyValue::Scalar {value, dims}, false)
}
fn assemble(claims: ClaimSet, phase: &str) -> MaterialCard {
    MaterialCard::assemble(MaterialStateId {chemistry: "manufactured frame medium".into(),
        phase: phase.into(), process: "analytic test, not measurements".into(), revision: 0},
        claims, vec![]).unwrap()
}
fn curve(name: &str, dims: fs_qty::Dims, low: f64, high: f64, solid: bool) -> PropertyClaim {
    claim(PropertyKey::new(name, dims), PropertyValue::Curve {abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS, knots: vec![(290.0, low), (500.0, high)], dims}, solid)
}
pub(crate) fn solid_card(alpha: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for c in [scalar("density", Density::DIMS, 1000.0),
        scalar("young_modulus", Pressure::DIMS, 2e6),
        scalar("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS, 10.0),
        curve("axial_eta", DynViscosity::DIMS, 500.0*(1.0-10.0*alpha), 500.0*(1.0+200.0*alpha), true),
        curve("bending_eta", DynViscosity::DIMS, 500.0*(1.0-10.0*alpha), 500.0*(1.0+200.0*alpha), true)] {
        claims.insert_claim(c).unwrap();
    }
    assemble(claims, "solid")
}
pub(crate) fn frame(card: &MaterialCard) -> Frame2d {
    let requirements: Vec<_> = [("density", Density::DIMS), ("young_modulus", Pressure::DIMS),
        ("axial_eta", DynViscosity::DIMS), ("bending_eta", DynViscosity::DIMS)]
        .into_iter().map(|(n,d)| ScalarPropertyRequirement::try_new(n,d,ScalarAdmissibility::Finite).unwrap()).collect();
    let material = resolve_material_state_point(card,&point(300.0),&requirements,
        MaterialPropertySelection::SingleClaimOnly).unwrap();
    let section = ConcentricElasticSection::try_new(0.0,vec![ConcentricMaterialLayer {
        outer_radius_m: 0.02, material}]).unwrap();
    let section = FrameViscousSection::try_new(section,"axial_eta","bending_eta").unwrap();
    let mut f = Frame2d::with_mass_formulation(Geometry::Linear,FrameMassFormulation::ConsistentReference);
    f.add_node(0.0,0.0); f.add_node(1.0,0.0);
    f.fix(0,[true;3]); f.fix(1,[true,false,true]);
    f.add_material_viscous_element(0,1,3,section).unwrap();
    f
}
fn liquid_card(variable_density: bool, basis: LiquidViscosityBasis) -> MaterialCard {
    let mut claims = ClaimSet::new();
    let density = if variable_density {curve("density",Density::DIMS,400.0,500.0,false)}
        else {scalar("density",Density::DIMS,400.0)};
    claims.insert_claim(density).unwrap();
    claims.insert_claim(scalar("specific_heat_capacity",SPECIFIC_HEAT_CAPACITY_DIMS,4000.0)).unwrap();
    let quantity = QuantitySpec::semantic(SemanticType::new(QuantityKind::ThermalConductivity,ValueForm::Static));
    claims.insert_claim(claim(PropertyKey::with_quantity("thermal_conductivity",quantity),
        PropertyValue::Scalar {value:0.5,dims:THERMAL_CONDUCTIVITY_DIMS},false)).unwrap();
    let (name,dims,scale) = match basis {
        LiquidViscosityBasis::Dynamic => (DYNAMIC_VISCOSITY_PROPERTY,DynViscosity::DIMS,1.0),
        LiquidViscosityBasis::Kinematic => (KINEMATIC_VISCOSITY_PROPERTY,
            fs_material::liquid::KINEMATIC_VISCOSITY_DIMS,1.0/400.0),
    };
    claims.insert_claim(curve(name,dims,0.07*scale,0.7*scale,false)).unwrap();
    assemble(claims,"liquid")
}
fn law(card: &MaterialCard,basis:LiquidViscosityBasis) -> LiquidTemperatureLaw<'_> {
    LiquidTemperatureLaw::try_new(card,&point(300.0),"T",basis,
        MaterialPropertySelection::SingleClaimOnly).unwrap()
}
pub(crate) fn resolution() -> CylinderFlowResolution {
    CylinderFlowResolution {minimum_wavenumber_radius:0.01,maximum_wavenumber_radius:10.0,intervals:3}
}
pub(crate) fn model<F: FluidTransport>(frame:&Frame2d, fluid:F) -> MaterialFluidFrame<'_,F> {
    MaterialFluidFrame::new(frame,1,vec![FrameFluidBinding {element:0,fluid}],resolution(),
        FrameVibrationBudget {max_frame_dofs:6,max_dense_entries:36},
        FrameFluidBudget {max_states:14,max_dense_entries:196,max_projection_terms:1024,
            mass_inverse_tolerance:1e-10},&CancelGate::new_clock_free()).unwrap()
}
pub(crate) fn near(a:f64,b:f64,tol:f64) {
    assert!((a-b).abs()<=tol,"{a:.16e} != {b:.16e}, tolerance={tol:e}");
}

#[test]
fn prepared_physical_factors_match_original_dense_owner_with_all_memory() {
    let solid=solid_card(0.0); let f=frame(&solid);
    let card=liquid_card(false,LiquidViscosityBasis::Dynamic); let law=law(&card,LiquidViscosityBasis::Dynamic);
    let model=model(&f,law.reference_state().clone()); let gate=CancelGate::new_clock_free();
    let mut prepared=model.prepare(&gate).unwrap();
    assert_eq!(prepared.condensed_dimension(),Some(3));
    let mut a=model.state_from_motion(&[ModalAcousticState {displacement_m_sqrt_kg:0.002,
        velocity_m_sqrt_kg_per_s:-0.003}],&[0.001;12]).unwrap();
    let mut b=a.clone();
    let mut condensed=0;
    for _ in 0..20 {
        let dense=model.step_candidate(&a,&[0.0,0.0,0.0,0.0,0.03,0.0],0.002,1e-9,&gate).unwrap();
        let next=prepared.step_candidate(&b,&[0.0,0.0,0.0,0.0,0.03,0.0],0.002,1e-9,&gate).unwrap();
        for (x,y) in dense.native.x.iter().zip(&next.state) {near(*x,*y,2e-9);}
        near(dense.nodal_work_j,next.nodal_work_j,1e-10);
        near(dense.fluid_member_dissipation_j[0],next.fluid_member_dissipation_j[0],1e-10);
        condensed+=next.linear_solve_counts.0;
        a=dense.native.x; b=next.state;
    }
    assert!(condensed>0);
}

#[test]
fn prescribed_temperature_changes_original_source_loss_without_resetting_inertia() {
    let solid=solid_card(0.0); let f=frame(&solid); let gate=CancelGate::new_clock_free();
    for basis in [LiquidViscosityBasis::Dynamic,LiquidViscosityBasis::Kinematic] {
        let card=liquid_card(false,basis); let source=law(&card,basis);
        let model=model(&f,source.reference_state().clone());
        let state=model.state_from_motion(&[ModalAcousticState {displacement_m_sqrt_kg:0.002,
            velocity_m_sqrt_kg_per_s:0.004}],&[0.003;12]).unwrap();
        let saved=state.clone();
        let inertia=model.modal_mass_matrix().to_vec(); let masses=model.memory_masses_kg().to_vec();
        let mut prepared=LiquidViscosityFrame::new(&model,vec![source.clone()],&gate).unwrap();
        let a=prepared.step_at_temperatures(&state,&[0.0;6],&[300.0],0.001,1e-9,&gate).unwrap();
        let b=prepared.step_at_temperatures(&state,&[0.0;6],&[400.0],0.001,1e-9,&gate).unwrap();
        assert!((a.mechanical.state[0]-b.mechanical.state[0]).abs()>1e-12);
        near(b.liquids[0].dynamic_viscosity_pa_s()/a.liquids[0].dynamic_viscosity_pa_s(),4.0,2e-14);
        assert_eq!(state,saved); assert_eq!(model.modal_mass_matrix(),inertia);
        assert_eq!(model.memory_masses_kg(),masses);
        for property in b.liquids[0].material().properties() {
            assert_eq!(property.answer().receipt.selected,
                source.reference_state().material().property(property.requirement().name()).unwrap().answer().receipt.selected);
        }
        assert!(b.mechanical.state[2..].iter().all(|p| *p!=0.0));
    }
}

#[test]
fn updated_source_step_matches_cold_native_hot_fluid_without_history_transfer() {
    let solid=solid_card(0.0); let f=frame(&solid); let gate=CancelGate::new_clock_free();
    let card=liquid_card(false,LiquidViscosityBasis::Dynamic); let law=law(&card,LiquidViscosityBasis::Dynamic);
    let cold=model(&f,law.reference_state().clone());
    let hot=model(&f,law.at_temperature(Temperature::new(400.0)).unwrap());
    let motion=[ModalAcousticState {displacement_m_sqrt_kg:0.001,velocity_m_sqrt_kg_per_s:-0.003}];
    let xc=cold.state_from_motion(&motion,&[0.002;12]).unwrap();
    let xh=hot.state_from_motion(&motion,&[0.002;12]).unwrap();
    let mut prepared=LiquidViscosityFrame::new(&cold,vec![law],&gate).unwrap();
    let next=prepared.step_at_temperatures(&xc,&[0.0,0.0,0.0,0.0,0.1,0.0],&[400.0],0.001,1e-9,&gate).unwrap();
    let reference=hot.step_candidate(&xh,&[0.0,0.0,0.0,0.0,0.1,0.0],0.001,1e-9,&gate).unwrap();
    for (a,b) in next.mechanical.state.iter().zip(reference.native.x) {near(*a,b,2e-9);}
}

#[test]
fn changing_density_source_or_incomplete_binding_is_not_a_viscosity_update() {
    let solid=solid_card(0.0); let f=frame(&solid); let gate=CancelGate::new_clock_free();
    let card=liquid_card(true,LiquidViscosityBasis::Dynamic); let source=law(&card,LiquidViscosityBasis::Dynamic);
    let model=model(&f,source.reference_state().clone());
    assert!(LiquidViscosityFrame::new(&model,vec![source],&gate).is_err());
    assert!(LiquidViscosityFrame::new(&model,vec![],&gate).is_err());
}

#[test]
fn source_domain_and_native_refusals_leave_history_available_for_retry() {
    let solid=solid_card(0.0); let f=frame(&solid);
    let card=liquid_card(false,LiquidViscosityBasis::Dynamic); let law=law(&card,LiquidViscosityBasis::Dynamic);
    let model=model(&f,law.reference_state().clone()); let gate=CancelGate::new_clock_free();
    let mut prepared=LiquidViscosityFrame::new(&model,vec![law],&gate).unwrap();
    let state=model.state_from_motion(&[ModalAcousticState::default()],&[0.001;12]).unwrap();
    let saved=state.clone();
    for temperatures in [&[501.0][..],&[f64::NAN][..],&[][..]] {
        assert!(prepared.step_at_temperatures(&state,&[0.0;6],temperatures,0.001,1e-9,&gate).is_err());
    }
    assert!(prepared.step_at_temperatures(&state,&[f64::NAN;6],&[300.0],0.001,1e-9,&gate).is_err());
    let retry=prepared.step_at_temperatures(&state,&[0.0;6],&[300.0],0.001,1e-9,&gate).unwrap();
    assert!(retry.mechanical.fluid_member_dissipation_j[0]>0.0);
    gate.request();
    assert!(prepared.step_at_temperatures(&state,&[0.0;6],&[300.0],0.001,1e-9,&gate).is_err());
    assert_eq!(state,saved);
}
