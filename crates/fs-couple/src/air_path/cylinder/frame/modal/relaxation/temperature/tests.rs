//! Original source curves, independent layer temperatures and continuing histories.
use super::*;
use crate::air_path::cylinder::{CylinderFlowResolution, frame::{FrameFluidBinding, FrameFluidBudget}};
use crate::frame_vibration::FrameVibrationBudget;
use crate::modal_acoustic_time::{ModalAcousticState, ModalAcousticTimeBudget};
use crate::render::schedule::force::coupled::ModalCouplingConfig;
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, UncertaintyModel};
use fs_material::fluid_transport::FluidTransportState;
use fs_material::state_point::{MaterialPropertySelection, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point};
use fs_qty::{Density, Pressure, Time};
use fs_solid::frame2d::{Frame2d, FrameMassFormulation, Geometry};

const RATE: u32 = 256;
const RADII: [f64; 2] = [0.06, 0.1];
const TIMES: [f64; 2] = [0.025, 0.12];
fn query(t: f64) -> QueryPoint {
    QueryPoint::new().with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t).unwrap()
        .with("omega", 1.0).unwrap()
}
fn card(layer: usize, constant_time: bool, curved_storage: bool, foreign: bool) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        ("density", Density::DIMS, [3.0, 5.0][layer]),
        ("young_modulus", Pressure::DIMS, [1000.0, 1500.0][layer]),
        ("equilibrium", Pressure::DIMS, [1000.0, 1500.0][layer]),
        ("increment", Pressure::DIMS, [600.0, 900.0][layer]),
        ("time", Time::DIMS, TIMES[layer]),
    ] {
        let curve = name == "time" && !constant_time || name == "increment" && curved_storage;
        let values = if name == "time" { [value*1.2, value, value*0.2] }
            else { [value, value, value] };
        let mut validity = ValidityDomain::unconstrained().with_quantity("T",
            QuantitySpec::dimensional(Temperature::DIMS), 290.0, 340.0);
        if name == "time" { validity = validity.with("omega", 0.0, 10_000.0); }
        claims.insert_claim(PropertyClaim { key: PropertyKey::new(name, dims),
            value: if curve { PropertyValue::Curve { abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS, knots: vec![(290.0,values[0]),(300.0,values[1]),(340.0,values[2])], dims }
                } else { PropertyValue::Scalar { value, dims } },
            validity, uncertainty: UncertaintyModel::Unstated,
            interpolation: if curve { InterpolationPolicy::LinearInside }
                else { InterpolationPolicy::ConstantWithinValidity },
            observations: vec![], provenance: Provenance { source: "manufactured thermal mobility".into(),
                license: "test-only".into(), artifact: None } }).unwrap();
    }
    MaterialCard::assemble(MaterialStateId { chemistry: format!("layer {layer}"),
        phase: "solid".into(), process: String::from(if foreign {"foreign reference"} else {"synthetic reference"}), revision: 0 }, claims, vec![]).unwrap()
}
fn make_frame(cards: &[MaterialCard; 2]) -> Frame2d {
    let requirements = [("density",Density::DIMS),("young_modulus",Pressure::DIMS),
        ("equilibrium",Pressure::DIMS),("increment",Pressure::DIMS),("time",Time::DIMS)]
        .map(|(name,dims)| ScalarPropertyRequirement::try_new(name,dims,ScalarAdmissibility::Finite).unwrap());
    let layers = cards.iter().enumerate().map(|(layer,card)| ConcentricMaterialLayer {
        outer_radius_m: RADII[layer], material: resolve_material_state_point(card,&query(300.0),
            &requirements,MaterialPropertySelection::SingleClaimOnly).unwrap() }).collect();
    let section = ConcentricElasticSection::try_new(0.0,layers).unwrap();
    let mut frame = Frame2d::with_mass_formulation(Geometry::Linear,FrameMassFormulation::ConsistentReference);
    frame.add_node(0.0,0.0); frame.add_node(1.0,0.0);
    frame.fix(0,[true;3]); frame.fix(1,[true,false,true]);
    frame.add_material_element(0,1,3,section).unwrap();
    frame
}
fn make_model(frame: &Frame2d) -> MaterialFluidFrame<'_,FluidTransportState> {
    MaterialFluidFrame::new(frame,1,vec![FrameFluidBinding {element:0,
        fluid:FluidTransportState::try_new(2.0,4000.0,0.5,0.02).unwrap()}],
        CylinderFlowResolution {minimum_wavenumber_radius:0.01,maximum_wavenumber_radius:2.0,intervals:2},
        FrameVibrationBudget {max_frame_dofs:6,max_dense_entries:36},
        FrameFluidBudget {max_states:16,max_dense_entries:256,max_projection_terms:4096,mass_inverse_tolerance:1e-10},
        &CancelGate::new_clock_free()).unwrap()
}
fn native<'m,'f>(model: &'m MaterialFluidFrame<'f,FluidTransportState>)
    -> (FrameFluidPlacement<'m,'f,FluidTransportState>,FrameRelaxationBinding<'m,'f,FluidTransportState>,CoupledModalSystem) {
    let gate = CancelGate::new_clock_free();
    let state = model.state_from_motion(&[ModalAcousticState {displacement_m_sqrt_kg:0.0003,
        velocity_m_sqrt_kg_per_s:-0.001}],&[0.0001;8]).unwrap();
    let config = ModalCouplingConfig {max_modes:16,max_connections:32,max_setup_terms:1_000_000,
        nyquist_guard_fraction:0.9,maximum_total_energy_j:100.0,maximum_abs_pressure_pa:1e6,
        maximum_abs_connection_force_n:1e6,solve_relative_tolerance:1e-10,
        energy_absolute_tolerance_j:1e-11,energy_relative_tolerance:1e-9};
    let (placement,system) = model.coupled_component(&state,RATE,ModalAcousticTimeBudget::audible_reference(),
        config,&gate).unwrap().into_system(None,&gate).unwrap();
    // A foreign, constrained material history must not be refreshed or reset.
    let extent = 2*system.mode_count();
    let system = system.with_material_relaxation(vec![RelaxationBranch {
        projection:vec![0.0;extent],stiffness:2.0,relaxation_time_s:0.7}],
        Some(vec![0.007]),&gate).unwrap();
    let selected = [FrameLongitudinalRelaxation {element:0,layers:vec![LayerRelaxationProperties {
        equilibrium_modulus:"equilibrium".into(),branches:vec![("increment".into(),"time".into())]};2]}];
    let (binding,system) = placement.with_longitudinal_relaxation(system,&selected,
        Some(vec![0.0002,0.003,-0.001,0.0005,-0.002,0.0007]),
        FrameRelaxationBudget {max_members:1,max_source_arms:2,max_memory_coordinates:6,max_setup_terms:100_000},&gate).unwrap();
    (placement,binding,system)
}
fn temperatures(cards: &[MaterialCard;2], t: [f64;2]) -> [FrameRelaxationLayerTemperature<'_>;2] {
    core::array::from_fn(|layer| FrameRelaxationLayerTemperature {element:0,layer,card:&cards[layer],
        temperature_axis:"T",temperature_k:t[layer]})
}
fn near(a:f64,b:f64) {assert!((a-b).abs()<2e-10*(1.0+a.abs()+b.abs()),"{a:e} != {b:e}");}
fn snapshot(binding:&FrameRelaxationBinding<'_,'_,FluidTransportState>) -> Vec<ConcentricPronySection> {
    binding.members().iter().map(|m|m.source.clone()).collect()
}

#[test]
fn independent_layer_temperatures_preserve_complete_history_and_change_continuing_wet_motion() {
    let cards = [card(0,false,false,false),card(1,false,false,false)];
    let frame = make_frame(&cards); let model = make_model(&frame); let gate = CancelGate::new_clock_free();
    let (placement,mut binding,mut system) = native(&model);
    let (_,_,mut control) = native(&model);
    let mut force = placement.project_nodal_force(&[0.0,0.0,0.0,0.0,0.002,0.0],&gate).unwrap();
    force.resize(system.mode_count(),0.0);
    for _ in 0..5 { system.step_under_gate(&force,&gate).unwrap();control.step_under_gate(&force,&gate).unwrap(); }
    let physical = placement.physical_state(system.components(),&gate).unwrap();
    let memories = system.material_memory().to_vec(); let clock = system.samples_rendered();
    let energy = system.total_energy_j().unwrap(); let old = snapshot(&binding);
    let report = binding.set_layer_temperatures(&mut system,&temperatures(&cards,[320.0,330.0]),&gate).unwrap();
    assert_eq!(report.changes.len(),6);
    assert_eq!(placement.physical_state(system.components(),&gate).unwrap(),physical);
    assert_eq!(system.material_memory(),memories);
    assert_eq!(system.total_energy_j().unwrap(),energy);
    assert_eq!(system.samples_rendered(),clock);
    assert_eq!(system.material_relaxation_branches()[0].relaxation_time_s,0.7);
    for (layer,factor) in [0.6,0.4].into_iter().enumerate() {
        near(binding.members()[0].source.spectra()[layer].terms[0].1,TIMES[layer]*factor);
        let original = &old[0].elastic().layers()[layer].material;
        let current = &binding.members()[0].source.elastic().layers()[layer].material;
        assert_eq!(original.card_identity(),current.card_identity());
        assert_eq!(current.query_point().iter().find(|(name,_)|name=="T").unwrap().1,[320.0,330.0][layer]);
        for property in original.properties() {
            assert_eq!(property.answer().receipt.selected,
                current.property(property.requirement().name()).unwrap().answer().receipt.selected);
        }
    }
    binding.states(&system,&gate).unwrap();
    let mut work = 0.0; let mut loss = 0.0;
    for _ in 0..20 {
        let step = system.step_under_gate(&force,&gate).unwrap().clone();
        work += step.external_work_j; loss += step.total_dissipation_j();
        control.step_under_gate(&force,&gate).unwrap();
        near(system.total_energy_j().unwrap()-energy+loss,work);
    }
    let changed = placement.physical_state(system.components(),&gate).unwrap();
    let unchanged = placement.physical_state(control.components(),&gate).unwrap();
    assert!((changed[0]-unchanged[0]).abs()>1e-10);
    assert!(changed[2..].iter().any(|p|*p!=0.0));
    assert_eq!(model.frame().material_section(0).unwrap().layers()[0].material.query(),query(300.0));
}

#[test]
fn late_source_refusal_foreign_card_and_cancel_preserve_both_owners_for_retry() {
    let cards = [card(0,false,false,false),card(1,false,false,false)];
    let frame = make_frame(&cards); let model = make_model(&frame); let gate = CancelGate::new_clock_free();
    let (placement,mut binding,mut system) = native(&model);
    let force = vec![0.0;system.mode_count()]; system.step_under_gate(&force,&gate).unwrap();
    let original = snapshot(&binding);let physical = placement.physical_state(system.components(),&gate).unwrap();
    let memories = system.material_memory().to_vec();let clock = system.samples_rendered();
    // First layer succeeds, second fails after some complete candidate queries.
    assert!(binding.set_layer_temperatures(&mut system,&temperatures(&cards,[320.0,350.0]),&gate).is_err());
    let foreign = card(1,false,false,true);let mut inputs = temperatures(&cards,[320.0,330.0]);inputs[1].card = &foreign;
    assert!(binding.set_layer_temperatures(&mut system,&inputs,&gate).is_err());
    let cancelled = CancelGate::new_clock_free();cancelled.request();
    assert!(binding.set_layer_temperatures(&mut system,&temperatures(&cards,[320.0,330.0]),&cancelled).is_err());
    assert!(binding.set_layer_temperatures(&mut system,&[],&gate).is_err());
    let mut reversed = temperatures(&cards,[320.0,330.0]);reversed.swap(0,1);
    assert!(binding.set_layer_temperatures(&mut system,&reversed,&gate).is_err());
    assert_eq!(snapshot(&binding),original);
    assert_eq!(placement.physical_state(system.components(),&gate).unwrap(),physical);
    assert_eq!(system.material_memory(),memories);assert_eq!(system.samples_rendered(),clock);
    binding.set_layer_temperatures(&mut system,&temperatures(&cards,[320.0,330.0]),&gate).unwrap();
    system.step_under_gate(&force,&gate).unwrap();assert_eq!(system.samples_rendered(),clock+1);
}

#[test]
fn unchanged_rates_refresh_original_receipts_but_curve_storage_and_stale_native_edits_refuse() {
    let gate = CancelGate::new_clock_free();
    let cards = [card(0,true,false,false),card(1,true,false,false)];
    let frame = make_frame(&cards); let model = make_model(&frame);
    let (_,mut binding,mut system) = native(&model);
    let initial = system.material_memory().to_vec();
    assert!(binding.set_layer_temperatures(&mut system,&temperatures(&cards,[320.0,330.0]),&gate).unwrap().changes.is_empty());
    assert_eq!(system.material_memory(),initial);
    assert_eq!(binding.members()[0].source.elastic().layers()[1].material.query(),query(330.0));
    let before = snapshot(&binding);
    system.set_material_relaxation_times(&[MaterialRelaxationTime {branch:binding.memory_range().start,
        relaxation_time_s:0.009}],&gate).unwrap();
    assert!(binding.set_layer_temperatures(&mut system,&temperatures(&cards,[320.0,330.0]),&gate).is_err());
    assert_eq!(snapshot(&binding),before);
    // Even a flat curve at all sampled knots is not a validity-wide scalar
    // storage claim. The rate-only owner must not guess its thermal potential.
    let curved = [card(0,false,true,false),card(1,false,false,false)];
    let frame2 = make_frame(&curved); let model2 = make_model(&frame2);
    let (_,mut binding2,mut system2) = native(&model2);
    assert!(binding2.set_layer_temperatures(&mut system2,&temperatures(&curved,[300.0,300.0]),&gate).is_err());
    assert_eq!(system2.samples_rendered(),0);
}
