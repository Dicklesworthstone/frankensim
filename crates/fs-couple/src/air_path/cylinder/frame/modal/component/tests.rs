//! G0/G1/G3: native complete equations, continuing histories and source-only swaps.
mod ports;
use super::*;
use crate::air_path::cylinder::{CylinderFlowResolution, frame::{FrameFluidBudget, FrameFluidRegion}};
use crate::frame_vibration::FrameVibrationBudget;
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::fluid_transport::FluidTransportState;
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::state_point::{MaterialPropertySelection, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point};
use fs_phs::linear_exact::{LinearSplitBudget, PreparedExactZoh};
use fs_qty::{Density, DynViscosity, Pressure};
use fs_solid::frame2d::{Frame2d, FrameMassFormulation, FrameViscousSection, Geometry};

fn frame(two: bool, density: f64, young: f64, eta: f64) -> Frame2d {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, dims, value) in [("density", Density::DIMS, density),
        ("young_modulus", Pressure::DIMS, young), ("eta", DynViscosity::DIMS, eta)] {
        let mut validity = ValidityDomain::unconstrained().with("T", 300.0, 400.0);
        if name == "eta" { validity = validity.with("omega", 0.0, 1e6); }
        claims.insert_claim(PropertyClaim { key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims }, validity,
            uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![], provenance: Provenance { source: "manufactured complete-fluid assembly".into(),
                license: "CC0-1.0".into(), artifact: None } }).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap());
    }
    let card = MaterialCard::assemble(MaterialStateId { chemistry: "synthetic".into(), phase: "solid".into(),
        process: "test-only constants".into(), revision: 0 }, claims, vec![]).unwrap();
    let point = QueryPoint::new().with("T", 300.0).unwrap().with("omega", 1.0).unwrap();
    let material = resolve_material_state_point(&card, &point, &requirements,
        MaterialPropertySelection::SingleClaimOnly).unwrap();
    let elastic = ConcentricElasticSection::try_new(0.0, vec![ConcentricMaterialLayer {
        outer_radius_m: 0.02, material }]).unwrap();
    let section = FrameViscousSection::try_new(elastic, "eta", "eta").unwrap();
    let mut result = Frame2d::with_mass_formulation(Geometry::Linear, FrameMassFormulation::ConsistentReference);
    result.add_node(0.0, 0.0); result.add_node(if two { 0.5 } else { 1.0 }, 0.0);
    if two { result.add_node(1.0, 0.0); }
    result.fix(0, [true; 3]); result.fix(1, [true, false, !two]);
    if two { result.fix(2, [true, false, false]); }
    for i in 0..(if two { 2 } else { 1 }) {
        result.add_material_viscous_element(i, i+1, 3, section.clone()).unwrap();
    }
    result
}
fn model(frame: &Frame2d, two: bool, density: f64, viscosity: f64)
    -> MaterialFluidFrame<'_, FluidTransportState> {
    let fluid = FluidTransportState::try_new(density, 4000.0, 0.5, viscosity).unwrap();
    let regions = if two { vec![
        FrameFluidRegion { element: 1, span: [0.0, 0.4], fluid },
        FrameFluidRegion { element: 1, span: [0.4, 1.0],
            fluid: FluidTransportState::try_new(3.0, 4000.0, 0.5, 0.002).unwrap() },
    ] } else { vec![FrameFluidRegion { element: 0, span: [0.0, 1.0], fluid }] };
    MaterialFluidFrame::new_with_regions(frame, if two {4} else {1}, regions,
        CylinderFlowResolution { minimum_wavenumber_radius: 0.01, maximum_wavenumber_radius: 2.0, intervals: 2 },
        FrameVibrationBudget { max_frame_dofs: 12, max_dense_entries: 144 },
        FrameFluidBudget { max_states: 64, max_dense_entries: 4096, max_projection_terms: 4096,
            mass_inverse_tolerance: 1e-10 }, &CancelGate::new_clock_free()).unwrap()
}
fn config() -> ModalCouplingConfig {
    ModalCouplingConfig { max_modes: 48, max_connections: 96, max_setup_terms: 1_000_000,
        nyquist_guard_fraction: 0.9, maximum_total_energy_j: 100.0, maximum_abs_pressure_pa: 100.0,
        maximum_abs_connection_force_n: 1e6, solve_relative_tolerance: 1e-10,
        energy_absolute_tolerance_j: 1e-9, energy_relative_tolerance: 1e-8 }
}
fn prepared() -> PreparedModalBudget {
    PreparedModalBudget { max_modes: 48, max_dense_bytes: 128*1024*1024, max_multiply_terms: 20_000_000_000 }
}
fn initial(model: &MaterialFluidFrame<'_, FluidTransportState>) -> Vec<f64> {
    let motion: Vec<_> = (0..model.basis().modes().len()).map(|i| ModalAcousticState {
        displacement_m_sqrt_kg: 2e-5*(i+1) as f64,
        velocity_m_sqrt_kg_per_s: -0.003*(i+1) as f64 }).collect();
    let velocity: Vec<_> = (0..model.memory_masses_kg().len()).map(|i| 1e-4*(i%3+1) as f64).collect();
    model.state_from_motion(&motion, &velocity).unwrap()
}
fn near(a: f64, b: f64) {
    assert!((a-b).abs() <= 2e-8*(1.0+a.abs()+b.abs()), "{a:e} != {b:e}");
}
fn complete<'m, 'f>(model: &'m MaterialFluidFrame<'f, FluidTransportState>, state: &[f64], cached: bool)
    -> (FrameFluidPlacement<'m, 'f, FluidTransportState>, CoupledModalSystem) {
    model.coupled_component(state, 2000, ModalAcousticTimeBudget::audible_reference(),
        config(), &CancelGate::new_clock_free()).unwrap()
        .into_system(cached.then(prepared), &CancelGate::new_clock_free()).unwrap()
}

#[test]
fn original_complete_phs_trajectory_work_and_loss_match_shared_native_assembly() {
    for two in [false, true] { for cached in [false, true] {
        let frame = frame(two, 1000.0, 2e6, 500.0);
        let model = model(&frame, two, 400.0, 0.02);
        let mut state = initial(&model);
        let (placement, mut network) = complete(&model, &state, cached);
        let gate = CancelGate::new_clock_free();
        for (&a, b) in state.iter().zip(placement.physical_state(network.components(), &gate).unwrap()) { near(a,b); }
        near(network.total_energy_j().unwrap(), model.system().hamiltonian(&state));
        let mut exact = PreparedExactZoh::new(model.system(), 0.0005,
            LinearSplitBudget { max_modes: model.system().state_dim(), max_omega_dt: core::f64::consts::PI }).unwrap();
        let m = model.basis().modes().len();
        let dofs = model.basis().reference_displacement().len();
        let mut nodal = vec![0.0; dofs];
        let mut dry_force = vec![0.0; m];
        let mut force = vec![0.0; network.mode_count()];
        for step in 0..12 {
            nodal[dofs-2] = if step<5 { 0.01 } else { -0.004 };
            model.basis().project_force_into(&nodal, [0.0;2], &mut dry_force, ||false).unwrap();
            force[..m].copy_from_slice(&placement.project_nodal_force(&nodal, &gate).unwrap());
            let reference = exact.step(&state, &dry_force).unwrap();
            let actual = network.step(&force, &gate).unwrap();
            for (&a, b) in reference.x.iter().zip(placement.physical_state(network.components(), &gate).unwrap()) { near(a,b); }
            near(actual.external_work_j, reference.supplied);
            near(actual.total_dissipation_j(), reference.dissipated);
            near(network.total_energy_j().unwrap(), model.system().hamiltonian(&reference.x));
            state = reference.x;
        }
        assert_eq!(network.samples_rendered(), 12);
        assert!(state[2*m..].iter().any(|p| *p != 0.0));
    }}
}

#[test]
fn complete_component_maps_every_original_loss_range_and_relocates_all_addresses() {
    let frame = frame(true, 1000.0, 2e6, 500.0);
    let model = model(&frame, true, 400.0, 0.02);
    let gate = CancelGate::new_clock_free();
    let component = model.coupled_component(&initial(&model), 2000,
        ModalAcousticTimeBudget::audible_reference(), config(), &gate).unwrap();
    assert!(core::ptr::eq(component.source(), &model));
    let count = component.connection_count();
    let (map, components, connections) = component.into_parts(3, 7, &gate).unwrap();
    assert_eq!(map.component_range(), 3..5);
    assert_eq!(map.connections().elastic, 7..11);
    assert_eq!(map.connections().material, 11..17);
    assert_eq!(map.connections().fluid_regions, vec![17..29, 29..41]);
    assert_eq!(count, 34);
    assert_eq!(components.len(), 2);
    assert_eq!(connections.len(), count);
    for (i, c) in connections.iter().enumerate() {
        assert_eq!(c.left.component, 3);
        assert!([3,4].contains(&c.right.component));
        if i < 4 { assert!(c.stiffness_n_m > 0.0); assert_eq!(c.damping_n_s_m, 0.0); }
        else { assert_eq!(c.stiffness_n_m, 0.0); assert!(c.damping_n_s_m >= 0.0); }
    }
    assert!(map.physical_state(&components, &gate).is_err());
}

#[test]
fn dry_component_is_the_same_source_damped_frame_without_artificial_memory() {
    let frame = frame(false, 1000.0, 2e6, 500.0);
    let gate = CancelGate::new_clock_free();
    let model = MaterialFluidFrame::<FluidTransportState>::new_with_regions(&frame,1,vec![],
        CylinderFlowResolution { minimum_wavenumber_radius:0.01,maximum_wavenumber_radius:2.0,intervals:2 },
        FrameVibrationBudget {max_frame_dofs:6,max_dense_entries:36},
        FrameFluidBudget {max_states:2,max_dense_entries:4,max_projection_terms:100,mass_inverse_tolerance:1e-10},&gate).unwrap();
    let state = initial(&model);
    let mut dry = model.basis().coupled_component(2000, &[0.0], ModalAcousticTimeBudget::audible_reference(),100,||false).unwrap();
    dry.restore_states(&model.modal_states(&state).unwrap()).unwrap();
    let (_, mut original) = dry.into_system(config(), None, &gate).unwrap();
    let (map, mut network) = complete(&model, &state, false);
    assert_eq!(map.fluid_component(),None);
    assert!(map.connections().fluid_regions.is_empty());
    for _ in 0..8 {
        network.step(&[0.01], &gate).unwrap(); original.step(&[0.01], &gate).unwrap();
        let actual = map.physical_state(network.components(), &gate).unwrap();
        near(actual[0], original.components()[0].states()[0].displacement_m_sqrt_kg);
        near(actual[1], original.components()[0].states()[0].velocity_m_sqrt_kg_per_s);
    }
}

#[test]
fn data_only_solid_and_exterior_substitutions_change_native_forced_motion() {
    let mut outputs = Vec::new();
    for (rho, young, fluid_rho, mu) in [(1000.0,2e6,400.0,0.02),(2000.0,2e6,400.0,0.02),
        (1000.0,8e6,400.0,0.02),(1000.0,2e6,1600.0,0.02),(1000.0,2e6,400.0,0.2)] {
        let frame = frame(false,rho,young,500.0);
        let model = model(&frame,false,fluid_rho,mu);
        let zero = vec![0.0;model.system().state_dim()];
        let (map,mut network)=complete(&model,&zero,false);
        let gate=CancelGate::new_clock_free();
        let mut force=vec![0.0;network.mode_count()];
        force[0]=map.project_nodal_force(&[0.0,0.0,0.0,0.0,0.01,0.0],&gate).unwrap()[0];
        for _ in 0..24 { network.step(&force,&gate).unwrap(); }
        let (mut u,mut v)=([0.0;6],[0.0;6]);
        map.reconstruct_into(network.components(),&mut u,&mut v,&gate).unwrap();
        outputs.push(u[4]);
    }
    for &changed in &outputs[1..] { assert!((changed-outputs[0]).abs()>1e-12); }
}

#[test]
fn refusal_cancellation_and_retry_do_not_reset_source_or_any_memory() {
    let frame=frame(false,1000.0,2e6,500.0);
    let model=model(&frame,false,400.0,0.02);
    let state=initial(&model);let saved=state.clone();let gate=CancelGate::new_clock_free();
    let mut too_small=config();too_small.max_modes=8;
    assert!(model.coupled_component(&state,2000,ModalAcousticTimeBudget::audible_reference(),too_small,&gate).is_err());
    let mut too_small=config();too_small.max_connections=15;
    assert!(model.coupled_component(&state,2000,ModalAcousticTimeBudget::audible_reference(),too_small,&gate).is_err());
    assert!(model.coupled_component(&state,1_000_000,ModalAcousticTimeBudget::audible_reference(),config(),&gate).is_err());
    let (map,mut network)=complete(&model,&state,false);
    let (control_map,mut control)=complete(&model,&state,false);
    let before=map.physical_state(network.components(),&gate).unwrap();
    let force=vec![0.0;network.mode_count()];
    let mut bad=force.clone();bad[0]=1e100;
    assert!(network.step(&bad,&gate).is_err());
    assert_eq!(map.physical_state(network.components(),&gate).unwrap(),before);
    assert_eq!(network.samples_rendered(),0);
    let cancelled=CancelGate::new_clock_free();cancelled.request();
    assert!(network.step(&force,&cancelled).is_err());
    assert!(map.reconstruct_into(network.components(),&mut [0.0;6],&mut [0.0;6],&cancelled).is_err());
    let a=network.step(&force,&gate).unwrap();let b=control.step(&force,&gate).unwrap();
    assert_eq!(a,b);
    assert_eq!(map.physical_state(network.components(),&gate).unwrap(),control_map.physical_state(control.components(),&gate).unwrap());
    assert_eq!(state,saved);
}
