//! Whole-assembly condition changes use one native replacement, regardless of
//! declaration order, and cannot publish a successful prefix of source queries.
use super::*;
use super::super::shared::set_layer_temperatures;
use crate::render::schedule::force::coupled::{ModalAttachment, ModalConnection};

fn pair<'m, 'f>(a: &'m MaterialFluidFrame<'f, FluidTransportState>,
    b: &'m MaterialFluidFrame<'f, FluidTransportState>) -> (
        Vec<FrameFluidPlacement<'m, 'f, FluidTransportState>>,
        Vec<FrameRelaxationBinding<'m, 'f, FluidTransportState>>, CoupledModalSystem) {
    let gate = CancelGate::new_clock_free();
    let config = ModalCouplingConfig { max_modes: 32, max_connections: 64,
        max_setup_terms: 1_000_000, nyquist_guard_fraction: 0.9,
        maximum_total_energy_j: 100.0, maximum_abs_pressure_pa: 1e6,
        maximum_abs_connection_force_n: 1e6, solve_relative_tolerance: 1e-10,
        energy_absolute_tolerance_j: 1e-11, energy_relative_tolerance: 1e-9 };
    let mut models = Vec::new(); let mut links = Vec::new(); let mut placements = Vec::new();
    for (i, model) in [a, b].into_iter().enumerate() {
        let state = model.state_from_motion(&[ModalAcousticState {
            displacement_m_sqrt_kg: 0.0003*(i+1) as f64,
            velocity_m_sqrt_kg_per_s: -0.001 }], &[0.0001; 8]).unwrap();
        let (placement, components, connections) = model.coupled_component(&state, RATE,
            ModalAcousticTimeBudget::audible_reference(), config, &gate).unwrap()
            .into_parts(models.len(), links.len(), &gate).unwrap();
        models.extend(components); links.extend(connections); placements.push(placement);
    }
    let point = [0.0, 0.0, 0.0, 0.0, 1.0, 0.0];
    links.push(ModalConnection {
        left: ModalAttachment { component: placements[0].structural_component(),
            shapes: placements[0].project_nodal_force(&point, &gate).unwrap() },
        right: ModalAttachment { component: placements[1].structural_component(),
            shapes: placements[1].project_nodal_force(&point, &gate).unwrap() },
        stiffness_n_m: 4.0, damping_n_s_m: 0.02, rest_extension_m: 0.0 });
    let mut mechanical = CoupledModalSystem::new(models, links, config, &gate).unwrap();
    let extent = 2*mechanical.mode_count();
    let foreign = || vec![RelaxationBranch { projection: vec![0.0; extent],
        stiffness: 2.0, relaxation_time_s: 0.7 }];
    let selected = [FrameLongitudinalRelaxation { element: 0,
        layers: vec![LayerRelaxationProperties { equilibrium_modulus: "equilibrium".into(),
            branches: vec![("increment".into(), "time".into())] }; 2] }];
    let mut bindings = Vec::new();
    for placement in &placements {
        // Foreign histories occupy the prefix AND the gap between source frames.
        mechanical = mechanical.with_additional_material_relaxation(foreign(), Some(vec![0.007]), &gate).unwrap();
        let (binding, next) = placement.with_longitudinal_relaxation(mechanical, &selected,
            Some(vec![0.0002, 0.003, -0.001, 0.0005, -0.002, 0.0007]),
            FrameRelaxationBudget { max_members: 1, max_source_arms: 2,
                max_memory_coordinates: 6, max_setup_terms: 1_000_000 }, &gate).unwrap();
        bindings.push(binding); mechanical = next;
    }
    mechanical = mechanical.with_additional_material_relaxation(foreign(), Some(vec![0.007]), &gate).unwrap();
    (placements, bindings, mechanical)
}
fn whole_state(system: &CoupledModalSystem) -> Vec<ModalAcousticState> {
    system.components().iter().flat_map(|model| model.states().iter().copied()).collect()
}
fn sources(bindings: &[FrameRelaxationBinding<'_, '_, FluidTransportState>]) -> Vec<Vec<ConcentricPronySection>> {
    bindings.iter().map(snapshot).collect()
}

#[test]
fn reversed_shared_sources_match_sequential_conditions_without_a_history_reset() {
    let cards = [card(0,false,false,false), card(1,false,false,false)];
    let a = make_frame(&cards); let b = make_frame(&cards);
    let a = make_model(&a); let b = make_model(&b); let gate = CancelGate::new_clock_free();
    let (_, mut bindings, mut system) = pair(&a, &b);
    let (_, mut sequential, mut reference) = pair(&a, &b);
    let mut force = vec![0.0; system.mode_count()]; force[0] = 0.001;
    for _ in 0..3 {system.step_under_gate(&force,&gate).unwrap(); reference.step_under_gate(&force,&gate).unwrap();}
    let accepted = whole_state(&system); let memory = system.material_memory().to_vec();
    let energy = system.total_energy_j().unwrap(); let sample = system.samples_rendered();
    let first = temperatures(&cards,[320.0,330.0]); let second = temperatures(&cards,[310.0,315.0]);
    bindings.reverse();
    let report = set_layer_temperatures(&mut bindings,&mut system,&[&second,&first],&gate).unwrap();
    assert_eq!(report.changes.len(),12);
    assert!(report.changes.windows(2).all(|p|p[0].branch < p[1].branch));
    assert_eq!(whole_state(&system),accepted); assert_eq!(system.material_memory(),memory);
    assert_eq!(system.total_energy_j().unwrap(),energy); assert_eq!(system.samples_rendered(),sample);
    sequential[0].set_layer_temperatures(&mut reference,&first,&gate).unwrap();
    sequential[1].set_layer_temperatures(&mut reference,&second,&gate).unwrap();
    for _ in 0..8 {
        system.step_under_gate(&force,&gate).unwrap(); reference.step_under_gate(&force,&gate).unwrap();
        for (a,b) in whole_state(&system).iter().zip(whole_state(&reference)) {
            near(a.displacement_m_sqrt_kg,b.displacement_m_sqrt_kg);
            near(a.velocity_m_sqrt_kg_per_s,b.velocity_m_sqrt_kg_per_s);
        }
        for (&a,&b) in system.material_memory().iter().zip(reference.material_memory()) {near(a,b);}
    }
    for index in [0,7,14] {assert_eq!(system.material_relaxation_branches()[index].relaxation_time_s,0.7);}
    for binding in &bindings {binding.states(&system,&gate).unwrap();}
}

#[test]
fn late_body_refusal_and_cancellation_preserve_all_sources_and_retry() {
    let cards = [card(0,false,false,false),card(1,false,false,false)];
    let a = make_frame(&cards); let b = make_frame(&cards);
    let a = make_model(&a); let b = make_model(&b); let gate = CancelGate::new_clock_free();
    let (_,mut bindings,mut system) = pair(&a,&b);
    let force = vec![0.0;system.mode_count()];system.step_under_gate(&force,&gate).unwrap();
    let snapshots = sources(&bindings);let state=whole_state(&system);let memory=system.material_memory().to_vec();
    let first=temperatures(&cards,[320.0,330.0]);let invalid=temperatures(&cards,[315.0,350.0]);
    assert!(set_layer_temperatures(&mut bindings,&mut system,&[&first,&invalid],&gate).is_err());
    assert_eq!(sources(&bindings),snapshots);assert_eq!(whole_state(&system),state);assert_eq!(system.material_memory(),memory);
    let second=temperatures(&cards,[310.0,315.0]);let cancelled=CancelGate::new_clock_free();cancelled.request();
    assert!(set_layer_temperatures(&mut bindings,&mut system,&[&first,&second],&cancelled).is_err());
    assert!(set_layer_temperatures(&mut bindings,&mut system,&[&first],&gate).is_err());
    assert_eq!(sources(&bindings),snapshots);assert_eq!(system.samples_rendered(),1);
    set_layer_temperatures(&mut bindings,&mut system,&[&first,&second],&gate).unwrap();
    system.step_under_gate(&force,&gate).unwrap();assert_eq!(system.samples_rendered(),2);
}

#[test]
fn constant_rates_refresh_each_source_and_overlapping_memory_blocks_refuse() {
    let cards=[card(0,true,false,false),card(1,true,false,false)];
    let a=make_frame(&cards);let b=make_frame(&cards);let a=make_model(&a);let b=make_model(&b);
    let gate=CancelGate::new_clock_free();let (_,mut bindings,mut system)=pair(&a,&b);
    let first=temperatures(&cards,[320.0,330.0]);let second=temperatures(&cards,[310.0,315.0]);
    let memory=system.material_memory().to_vec();
    assert!(set_layer_temperatures(&mut bindings,&mut system,&[&first,&second],&gate).unwrap().changes.is_empty());
    assert_eq!(bindings[0].members()[0].source.elastic().layers()[1].material.query(),query(330.0));
    assert_eq!(bindings[1].members()[0].source.elastic().layers()[1].material.query(),query(315.0));
    assert_eq!(system.material_memory(),memory);
    let snapshots=sources(&bindings);
    let overlapping=bindings[0].memory_range(); bindings[1].memory=overlapping;
    assert!(set_layer_temperatures(&mut bindings,&mut system,&[&first,&second],&gate).is_err());
    assert_eq!(sources(&bindings),snapshots);assert_eq!(system.material_memory(),memory);
}
