//! G0/G1/G3: source geometry partitions, original materials and native motion.
use super::*;
use fs_solid::frame2d::{FrameMassFormulation, FrameViscousSection, Geometry, SectionModel};
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::gas::{GasSpec, GasState};
use fs_material::liquid::{ResolvedLiquidState, resolve_liquid_state,
    DYNAMIC_VISCOSITY_DIMS, DYNAMIC_VISCOSITY_PROPERTY};
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::state_point::{MaterialPropertySelection, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point, SPECIFIC_HEAT_CAPACITY_DIMS,
    THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY};
use fs_qty::{Density, DynViscosity, Pressure, QuantitySpec};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

fn geometry(points: [[f64; 2]; 2], reverse: bool, kind: Geometry) -> Frame2d {
    let mut frame = Frame2d::new(kind);
    for p in points { frame.add_node(p[0], p[1]); }
    let (i, j) = if reverse { (1, 0) } else { (0, 1) };
    frame.add_element(i, j, 3, 1.0, &|| SectionModel::Elastic { ea: 100.0, ei: 2.0 }).unwrap();
    frame
}
fn near(a: f64, b: f64) { assert!((a-b).abs() < 2e-11*(1.0+a.abs()+b.abs()), "{a:e} != {b:e}"); }
fn split(frame: &Frame2d, normal: [f64; 2], offset: f64) -> Result<Vec<FrameFluidRegion<&'static str>>> {
    FrameFluidRegion::split_at_reference_plane(frame, 0, normal, offset, "negative", "positive",
        &CancelGate::new_clock_free())
}

#[test]
fn crossing_uses_original_orientation_translation_and_signed_material_sides() {
    let frame = geometry([[2.0,-1.0],[5.0,3.0]], false, Geometry::Linear);
    let regions = split(&frame, [0.6,0.8], 2.4).unwrap();
    assert_eq!(regions.len(), 2);
    near(regions[0].span[1], 0.4);
    assert_eq!(regions[0].span[1], regions[1].span[0]);
    assert_eq!(regions[0].fluid, "negative");
    assert_eq!(regions[1].fluid, "positive");
    let reversed = geometry([[2.0,-1.0],[5.0,3.0]], true, Geometry::Linear);
    let r = split(&reversed, [0.6,0.8], 2.4).unwrap();
    near(r[0].span[1], 0.6);
    assert_eq!(r[0].fluid, "positive");
    assert_eq!(r[1].fluid, "negative");
    let moved = geometry([[6.0,6.0],[9.0,10.0]], false, Geometry::Linear);
    let r = split(&moved, [0.6,0.8], 10.4).unwrap();
    near(r[0].span[1], 0.4);
    let r = split(&frame, [-0.6,-0.8], -2.4).unwrap();
    near(r[0].span[1], 0.4);
    assert_eq!(r[0].fluid, "positive");
    assert_eq!(frame.displacements(), &[0.0; 6]);
}

#[test]
fn touching_endpoint_has_no_empty_region_and_coincident_member_refuses() {
    let frame = geometry([[0.0,0.0],[1.0,0.0]], false, Geometry::Linear);
    for (offset, source) in [(0.0,"positive"),(1.0,"negative"),(-1.0,"positive"),(2.0,"negative")] {
        let r = split(&frame, [1.0,0.0], offset).unwrap();
        assert_eq!(r.len(), 1);
        assert_eq!(r[0].span, [0.0,1.0]);
        assert_eq!(r[0].fluid, source);
    }
    assert!(split(&frame, [0.0,1.0], 0.0).is_err());
    assert_eq!(split(&frame, [0.0,1.0], 1.0).unwrap()[0].fluid, "negative");
}

#[test]
fn bad_planes_geometry_and_cancellation_preserve_inputs() {
    let frame = geometry([[0.0,0.0],[1.0,0.0]], false, Geometry::Linear);
    for normal in [[0.0,0.0],[2.0,0.0],[f64::NAN,1.0],[0.0,f64::INFINITY]] {
        assert!(split(&frame, normal, 0.5).is_err());
    }
    assert!(split(&frame, [1.0,0.0], f64::NAN).is_err());
    assert!(split(&geometry([[0.0,0.0],[1.0,0.0]],false,Geometry::Corotational),[1.0,0.0],0.5).is_err());
    let gate = CancelGate::new_clock_free();
    assert!(FrameFluidRegion::split_at_reference_plane(&frame, 9, [1.0,0.0], 0.5, (), (), &gate).is_err());
    gate.request();
    assert!(matches!(FrameFluidRegion::split_at_reference_plane(&frame, 0, [1.0,0.0], 0.5, (), (), &gate),
        Err(FrameFluidError::Cancelled)));
    assert_eq!(frame.displacements(), &[0.0; 6]);
}

fn source_card(phase: &str, values: Vec<(&str, fs_qty::Dims, f64)>) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in values {
        let key = if name == THERMAL_CONDUCTIVITY_PROPERTY {
            PropertyKey::with_quantity(name, QuantitySpec::semantic(SemanticType::new(
                QuantityKind::ThermalConductivity, ValueForm::Static)))
        } else { PropertyKey::new(name,dims) };
        let mut validity = ValidityDomain::unconstrained().with("T",300.0,400.0);
        if name == "eta" { validity = validity.with("omega",0.0,1e7); }
        claims.insert_claim(PropertyClaim { key, value: PropertyValue::Scalar {value,dims}, validity,
            uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![], provenance: Provenance {source:"synthetic mixed-fluid plane fixture".into(),
                license:"CC0-1.0".into(),artifact:None} }).unwrap();
    }
    MaterialCard::assemble(MaterialStateId {chemistry:"manufactured transport/solid".into(),
        phase:phase.into(),process:"reference data only".into(),revision:0},claims,vec![]).unwrap()
}
fn material_frame() -> Frame2d {
    let values = vec![("density",Density::DIMS,1000.0),("young_modulus",Pressure::DIMS,2e6),
        ("eta",DynViscosity::DIMS,500.0)];
    let requirements = values.iter().map(|(n,d,_)|
        ScalarPropertyRequirement::try_new(n,*d,ScalarAdmissibility::Finite).unwrap()).collect::<Vec<_>>();
    let card = source_card("solid",values);
    let point = QueryPoint::new().with("T",300.0).unwrap().with("omega",1.0).unwrap();
    let material = resolve_material_state_point(&card,&point,&requirements,MaterialPropertySelection::SingleClaimOnly).unwrap();
    let section = ConcentricElasticSection::try_new(0.0,vec![ConcentricMaterialLayer {
        outer_radius_m:0.02,material}]).unwrap();
    let mut frame = Frame2d::with_mass_formulation(Geometry::Linear,FrameMassFormulation::ConsistentReference);
    frame.add_node(0.0,0.0);frame.add_node(1.0,0.0);
    frame.fix(0,[true;3]);frame.fix(1,[true,false,true]);
    frame.add_material_viscous_element(0,1,3,FrameViscousSection::try_new(section,"eta","eta").unwrap()).unwrap();
    frame
}
fn liquid() -> ResolvedLiquidState {
    let card = source_card("liquid",vec![("density",Density::DIMS,1000.0),
        ("specific_heat_capacity",SPECIFIC_HEAT_CAPACITY_DIMS,4000.0),
        (THERMAL_CONDUCTIVITY_PROPERTY,THERMAL_CONDUCTIVITY_DIMS,0.6),
        (DYNAMIC_VISCOSITY_PROPERTY,DYNAMIC_VISCOSITY_DIMS,0.001)]);
    resolve_liquid_state(&card,&QueryPoint::new().with("T",300.0).unwrap(),MaterialPropertySelection::SingleClaimOnly).unwrap()
}

#[test]
fn original_liquid_and_gas_sources_share_native_motion_without_transport_erasure() {
    let frame = material_frame(); let liquid = liquid();
    let gate = CancelGate::new_clock_free(); let mut masses = Vec::new();
    for pressure in [100_000.0,200_000.0] {
        let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(),300.0,pressure).unwrap();
        let negative: &dyn FluidTransport = &liquid;
        let positive: &dyn FluidTransport = &gas;
        let regions = FrameFluidRegion::split_at_reference_plane(&frame,0,[1.0,0.0],0.4,
            negative,positive,&gate).unwrap();
        let sim = MaterialFluidFrame::new_with_regions(&frame,1,regions,
            CylinderFlowResolution {minimum_wavenumber_radius:0.01,maximum_wavenumber_radius:10.0,intervals:3},
            FrameVibrationBudget {max_frame_dofs:6,max_dense_entries:36},
            FrameFluidBudget {max_states:26,max_dense_entries:676,max_projection_terms:4096,mass_inverse_tolerance:1e-10},
            &gate).unwrap();
        assert!(core::ptr::eq(*sim.members()[0].flow().fluid(),negative));
        assert!(core::ptr::eq(*sim.members()[1].flow().fluid(),positive));
        assert_eq!(*sim.members()[0].flow().transport(),liquid.transport_state().unwrap());
        assert_eq!(*sim.members()[1].flow().transport(),gas.transport_state().unwrap());
        let primitive = |x:f64|9.0*x.powi(5)/5.0-2.0*x.powi(6)+4.0*x.powi(7)/7.0;
        let shape = sim.basis().modes()[0].displacement[4];
        let expected = 1.0 + core::f64::consts::PI*0.02_f64.powi(2)*shape*shape*
            (1000.0*primitive(0.4)+gas.density*(primitive(1.0)-primitive(0.4)));
        near(sim.modal_mass_matrix()[0],expected);
        masses.push(sim.modal_mass_matrix()[0]);
        let mut state = sim.state_from_motion(&[ModalAcousticState::default()],&[0.0;24]).unwrap();
        let mut prepared = sim.prepare(&gate).unwrap(); let mut heat = [0.0;2];
        for _ in 0..16 {
            let step = prepared.step_candidate(&state,&[0.0,0.0,0.0,0.0,0.01,0.0],0.001,1e-9,&gate).unwrap();
            assert_eq!(step.fluid_member_dissipation_j.len(),2);
            for (sum,loss) in heat.iter_mut().zip(&step.fluid_member_dissipation_j) { *sum += loss; }
            near(step.energy_residual_j,0.0); state = step.state;
        }
        assert!(state[0]*shape>0.0 && heat.iter().all(|q|*q>0.0));
        assert_eq!(frame.ndof(),6);
    }
    assert!(masses[1]>masses[0]);
}
