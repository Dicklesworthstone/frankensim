//! G1/G3: source thermoelasticity, retained exterior histories and reciprocal bodies.
use super::*;
use crate::air_path::cylinder::frame::harmonic_ports::FrameHarmonicPortBudget;
use crate::air_path::cylinder::frame::modal::FrameFluidPlacement;
use crate::render::schedule::force::coupled::ModalConnection;

const OTHER: CellId = CellId(23);
const JOINT_K: f64 = 7.0;
const JOINT_C: f64 = 0.03;
fn shared_config() -> FrameFluidThermoelasticConfig {
    let mut c = config();
    c.coupling.max_modes = 32;
    c.coupling.max_connections = 64;
    c.coupling.max_setup_terms = 1_000_000;
    c
}
fn port_budget() -> FrameHarmonicPortBudget {
    FrameHarmonicPortBudget { max_ports: 1, max_projection_terms: 128,
        max_mobility_entries: 1, power_absolute_tolerance_w: 1e-10,
        power_relative_tolerance: 1e-9 }
}
fn diagonal_column() -> Vec<f64> {
    let mut b = vec![0.0; 6];
    b[3] = core::f64::consts::FRAC_1_SQRT_2;
    b[4] = core::f64::consts::FRAC_1_SQRT_2;
    b
}
fn two_frames<R>(run: impl FnOnce(&MaterialFluidFrame<'_, FluidTransportState>,
    &MaterialFluidFrame<'_, FluidTransportState>, DepositionSimulation<'_>, &MaterialNetwork<'_>) -> R) -> R {
    let materials = [sources(-2000.0, 1e-4), sources(1000.0, 2e-4)];
    let frames = materials.each_ref().map(|(card, _)| frame(card));
    let models = frames.each_ref().map(|frame| model(frame, true));
    let cells_ids = [CELL, OTHER];
    let events = (0..2).map(|i| TimedEvent { time_s: 0.0,
        event: Event::Birth(DepositedCell { id: cells_ids[i], segment_id: i as u64,
            body_id: i as u64, material_id: MaterialId(i as u64), mass_kg: RHO*area(),
            specific_h_j_kg: 100.0*CP,
            geometry: CellGeometry { start_m: [0.0, i as f64, 0.0], end_m: [1.0, i as f64, 0.0],
                volume_m3: area(), surface_area_m2: 2.0*core::f64::consts::PI*RADIUS },
            boundary: Boundary { ambient_k: 300.0, htc_w_m2_k: 0.0, reservoir: None,
                source_w: if i == 0 { 10.0 } else { 0.0 } },
        }) }).collect();
    let laws = [(MaterialId(0), materials[0].1.curve() as &dyn ThermalLaw),
        (MaterialId(1), materials[1].1.curve() as &dyn ThermalLaw)];
    let thermal = DepositionSimulation::new(events, &laws,
        DepositionBudget { max_cells: 2, max_links: 1, max_events: 2, max_substeps: 1000 }).unwrap()
        .advance_to(0.0, StepConfig { max_dt_s: 0.001, energy_tolerance_j: 1e-10 }, || Ok(())).unwrap().0;
    let context = point(300.0);
    let cells: Vec<_> = materials.iter().enumerate().map(|(i,(card,curve))| MaterialCell {
        cell: cells_ids[i], material_id: MaterialId(i as u64), curve: curve.curve(),
        point: ThermalMaterialPoint { material: card, context: &context, temperature_axis: "T",
            liquid_fraction_axis: None }, conductivity: None,
    }).collect();
    run(&models[0], &models[1], thermal, &MaterialNetwork { cells: &cells, resistors: &[], links: &[] })
}
fn layer_total(map: &FrameFluidPlacement<'_, '_, FluidTransportState>,
    components: &[ModalAcousticTimeModel], thermal: &DepositionSimulation<'_>,
    cell: CellId, slope: f64, alpha: f64) -> f64 {
    let gate = CancelGate::new_clock_free();
    let state = map.physical_state(components, &gate).unwrap();
    let motion = map.source().modal_states(&state).unwrap();
    let (mut u, mut v) = ([0.0; 6], [0.0; 6]);
    map.reconstruct_into(components, &mut u, &mut v, &gate).unwrap();
    let m = motion.len();
    let kinetic = 0.5*(0..m).map(|i| state[m+i]*motion[i].velocity_m_sqrt_kg_per_s).sum::<f64>()
        + 0.5*state[2*m..].iter().zip(map.source().memory_masses_kg()).map(|(p,m)| p*p/m).sum::<f64>();
    let t = thermal.temperature_k(cell).unwrap();
    let e = YOUNG+slope*(t-300.0);
    let elastic = u[3]-alpha*(t-300.0);
    kinetic + area()*(0.5*(e-t*slope)*elastic*elastic+t*e*alpha*elastic)
        + 0.5*(e-t*slope)*12.0*moment()*u[4]*u[4] + RHO*area()*CP*(t-200.0)
}

#[test]
fn shared_parts_preserve_standalone_temperature_motion_and_exterior_history() {
    fixture(-2000.0, 1e-4, 10.0, true, |model, mut thermal, source| {
        let gate = CancelGate::new_clock_free();
        let initial = initial(model);
        let mut standalone_heat = thermal.clone();
        let mut standalone = FrameFluidThermoelasticSystem::new(model, &initial, &standalone_heat,
            source, &DESTINATIONS, config(), &gate).unwrap();
        let p = model.thermoelastic_parts(&initial, &thermal, source, &DESTINATIONS,
            config(), 0, 0, &gate).unwrap();
        let mechanical = CoupledModalSystem::new(p.components, p.connections, config().coupling, &gate).unwrap();
        let mut shared = p.binding.into_system(mechanical, &thermal, &gate).unwrap();
        for _ in 0..12 {
            let nodal = [0.0,0.0,0.0,0.0,0.001,0.0];
            let mut force = vec![0.0; shared.mechanical().mode_count()];
            force[..2].copy_from_slice(&p.placement.project_nodal_force(&nodal, &gate).unwrap());
            let a = shared.step(&mut thermal, &force, controls(), &[0.0], &gate).unwrap();
            let b = standalone.step(&mut standalone_heat, &nodal, controls(), &[0.0], &gate).unwrap();
            assert_eq!(p.placement.physical_state(shared.mechanical().components(), &gate).unwrap(),
                standalone.physical_state().unwrap());
            assert_eq!(thermal.cells(), standalone_heat.cells());
            near(a.total_energy_residual_j, b.total_energy_residual_j, 1e-12);
            assert_eq!(a.external_mechanical_loss_j, b.external_mechanical_loss_j);
        }
    });
}

#[test]
fn heating_one_wet_body_changes_the_other_on_one_shared_material_clock() {
    two_frames(|a, b, mut thermal, source| {
        let gate = CancelGate::new_clock_free();
        let ap = a.harmonic_ports(&[diagonal_column()], port_budget(), &gate).unwrap();
        let bp = b.harmonic_ports(&[diagonal_column()], port_budget(), &gate).unwrap();
        let pa = a.thermoelastic_parts(&vec![0.0; a.system().state_dim()], &thermal, source,
            &DESTINATIONS, shared_config(), 0, 0, &gate).unwrap();
        let pb = b.thermoelastic_parts(&vec![0.0; b.system().state_dim()], &thermal, source,
            &[FrameLayerCell { element: 0, layer: 0, cell: OTHER }], shared_config(),
            pa.components.len(), pa.connections.len(), &gate).unwrap();
        let amap = pa.placement.bind_ports(&ap, 128, &gate).unwrap();
        let bmap = pb.placement.bind_ports(&bp, 128, &gate).unwrap();
        let mut components = pa.components; components.extend(pb.components);
        let mut connections = pa.connections; connections.extend(pb.connections);
        connections.push(ModalConnection { left: amap.attachment(0).unwrap(), right: bmap.attachment(0).unwrap(),
            stiffness_n_m: JOINT_K, damping_n_s_m: JOINT_C, rest_extension_m: 0.0 });
        let mechanical = CoupledModalSystem::new(components, connections, shared_config().coupling, &gate).unwrap();
        // Deliberately reverse declarations: source reporting remains native-order.
        let mut owner = FrameFluidThermoelasticBinding::into_shared_system(vec![pb.binding, pa.binding],
            mechanical, &thermal, &gate).unwrap();
        assert_eq!(owner.frames()[0].component_index(), 0);
        assert_eq!(owner.frames()[1].component_index(), 2);
        assert_eq!(owner.members().iter().map(|m| m.component).collect::<Vec<_>>(), vec![0,2]);
        let start = layer_total(&pa.placement, owner.mechanical().components(), &thermal, CELL, -2000.0, 1e-4)
            + layer_total(&pb.placement, owner.mechanical().components(), &thermal, OTHER, 1000.0, 2e-4);
        let force = vec![0.0; owner.mechanical().mode_count()];
        let (mut heat, mut outside) = (0.0, 0.0);
        let (mut aq, mut av, mut bq, mut bv) = ([0.0], [0.0], [0.0], [0.0]);
        for _ in 0..32 {
            let report = owner.step(&mut thermal, &force, controls(), &[0.0; 2], &gate).unwrap();
            amap.observe_into(owner.mechanical().components(), &mut aq, &mut av, &gate).unwrap();
            bmap.observe_into(owner.mechanical().components(), &mut bq, &mut bv, &gate).unwrap();
            heat += report.external_thermal_heat_j;
            outside += report.external_mechanical_loss_j;
            let total = layer_total(&pa.placement, owner.mechanical().components(), &thermal, CELL, -2000.0, 1e-4)
                + layer_total(&pb.placement, owner.mechanical().components(), &thermal, OTHER, 1000.0, 2e-4)
                + 0.5*JOINT_K*(aq[0]-bq[0]).powi(2);
            near(total-start+outside-heat, 0.0, 2e-7);
            near(report.mechanical.external_work_j, 0.0, 1e-14);
            near(report.total_energy_residual_j, 0.0, 2e-8);
        }
        assert!(aq[0] > 0.0 && bq[0] > 1e-9);
        assert!(outside > 0.0);
        assert!(thermal.temperature_k(CELL).unwrap() > 300.0);
        assert!(owner.members()[0].source.layers[0].young_modulus_pa < YOUNG);
        assert!(pb.placement.physical_state(owner.mechanical().components(), &gate).unwrap()[4..]
            .iter().any(|p| *p != 0.0));
        assert_eq!(owner.mechanical().samples_rendered(), 32);
    });
}

#[test]
fn shared_admission_protects_fluid_history_coefficients_and_reciprocal_columns() {
    fixture(-2000.0, 1e-4, 10.0, true, |model, thermal, source| {
        let gate = CancelGate::new_clock_free();
        let state = initial(model);
        for fault in 0..4 {
            let mut p = model.thermoelastic_parts(&state, &thermal, source, &DESTINATIONS,
                config(), 0, 0, &gate).unwrap();
            let fluid_start = p.placement.connections().fluid_regions[0].start;
            match fault {
                0 => p.connections[fluid_start].damping_n_s_m *= 2.0,
                1 => p.connections[fluid_start+1].right.shapes.fill(0.0),
                2 => { let mut changed = p.components[1].states().to_vec();
                    changed[0].velocity_m_sqrt_kg_per_s *= 2.0;
                    p.components[1].restore_states(&changed).unwrap(); },
                _ => { p.connections.pop(); },
            }
            let mechanical = CoupledModalSystem::new(p.components, p.connections, config().coupling, &gate).unwrap();
            assert!(p.binding.into_system(mechanical, &thermal, &gate).is_err());
            assert_eq!(thermal.time_s(), 0.0);
        }
    });
}

#[test]
fn a_late_shared_heat_refusal_preserves_both_complete_bodies_and_caloric_cells() {
    two_frames(|a, b, mut thermal, source| {
        let gate = CancelGate::new_clock_free();
        let pa = a.thermoelastic_parts(&initial(a), &thermal, source, &DESTINATIONS,
            shared_config(), 0, 0, &gate).unwrap();
        let pb = b.thermoelastic_parts(&initial(b), &thermal, source,
            &[FrameLayerCell { element: 0, layer: 0, cell: OTHER }], shared_config(),
            pa.components.len(), pa.connections.len(), &gate).unwrap();
        let mut components = pa.components; components.extend(pb.components);
        let mut connections = pa.connections; connections.extend(pb.connections);
        let mechanical = CoupledModalSystem::new(components, connections, shared_config().coupling, &gate).unwrap();
        let mut owner = FrameFluidThermoelasticBinding::into_shared_system(vec![pa.binding,pb.binding],
            mechanical, &thermal, &gate).unwrap();
        let before_a = pa.placement.physical_state(owner.mechanical().components(), &gate).unwrap();
        let before_b = pb.placement.physical_state(owner.mechanical().components(), &gate).unwrap();
        let before_heat = thermal.cells().clone();
        let force = vec![0.0; owner.mechanical().mode_count()];
        let mut exhausted = controls(); exhausted.maximum_material_trials = 1;
        assert!(owner.step(&mut thermal, &force, exhausted, &[0.0; 2], &gate).is_err());
        assert_eq!(pa.placement.physical_state(owner.mechanical().components(), &gate).unwrap(), before_a);
        assert_eq!(pb.placement.physical_state(owner.mechanical().components(), &gate).unwrap(), before_b);
        assert_eq!(thermal.cells(), &before_heat);
        assert_eq!(owner.mechanical().samples_rendered(), 0);
        assert_eq!(thermal.time_s(), 0.0);
        let cancelled = CancelGate::new_clock_free(); cancelled.request();
        assert!(owner.step(&mut thermal, &force, controls(), &[0.0; 2], &cancelled).is_err());
        owner.step(&mut thermal, &force, controls(), &[0.0; 2], &gate).unwrap();
        assert_eq!(owner.mechanical().samples_rendered(), 1);
        assert!(thermal.temperature_k(CELL).unwrap() > 300.0);
    });
}
