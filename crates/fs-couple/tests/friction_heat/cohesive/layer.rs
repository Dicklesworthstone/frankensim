//! Finite adhesive storage from the real native joint, with analytic references.
use super::*;
use fs_couple::friction_heat::third_body::cohesive::*;
use fs_couple::friction_heat::third_body::{ThirdBodyHeatRoute, ThirdBodyHeatTolerances};

const LAYER: CellId = CellId(2);
const LINKS: [u64; 2] = [17, 18];
const CAPACITY: [f64; 3] = [2.0, 3.0, 0.025];
const K: [f64; 3] = [4.0, 1.0, 2.0];
const R: [f64; 2] = [0.2, 0.8];
const ADHESIVE: &str = "original finite adhesive";

fn layer_network() -> DepositionSimulation<'static> {
    let mut events = Vec::new();
    for i in 0..3 {
        let [start, end] = [[-1.0, 0.0], [0.1, 1.1], [0.0, 0.1]][i];
        let length = end - start;
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64),
                segment_id: i as u64,
                body_id: i as u64,
                material_id: MaterialId(i as u64),
                mass_kg: CAPACITY[i],
                specific_h_j_kg: 0.0,
                geometry: CellGeometry {
                    start_m: [start, 0.0, 0.0],
                    end_m: [end, 0.0, 0.0],
                    volume_m3: length,
                    surface_area_m2: 2.0 + 4.0 * length,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        });
    }
    for i in 0..2 {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: LINKS[i],
                a: CellId(i as u64),
                b: LAYER,
                area_m2: 1.0,
                conductance_w_k: series_contact_conductance_w_k(
                    1.0,
                    [0.5, 0.05],
                    [K[i], K[2]],
                    R[i],
                )
                .unwrap(),
            }),
        });
    }
    DepositionSimulation::new(
        events,
        &[
            (MaterialId(0), &SENSIBLE as &dyn ThermalLaw),
            (MaterialId(1), &SENSIBLE),
            (MaterialId(2), &SENSIBLE),
        ],
        DepositionBudget {
            max_cells: 3,
            max_links: 2,
            max_events: 5,
            max_substeps: 100,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.01,
            energy_tolerance_j: 1e-12,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

fn binding<'a>(
    route: &'a ThirdBodyHeatRoute,
    card: &'a InterfaceSystemCard,
) -> CohesiveLayerHeatBinding<'a> {
    CohesiveLayerHeatBinding {
        patch_id: 7,
        route,
        card,
        layer_material: MaterialId(2),
        layer_name: ADHESIVE,
        bodies: [
            CohesiveLayerHeatBody {
                body_id: 0,
                thermal_material: MaterialId(0),
                mechanical_material: &card.surface_a().material,
            },
            CohesiveLayerHeatBody {
                body_id: 1,
                thermal_material: MaterialId(1),
                mechanical_material: &card.surface_b().material,
            },
        ],
    }
}
fn layer_controls() -> CohesiveLayerHeatControls {
    CohesiveLayerHeatControls {
        maximum_cells: 3,
        maximum_trials: 128,
        maximum_newton_iterations: 20,
        difference_step_j: 1e-6,
        heat: ThirdBodyHeatTolerances {
            cell_j: 1e-11,
            thermal_j: 1e-11,
            layer_j: 1e-11,
            total_j: 1e-11,
        },
        combined_energy_tolerance_j: 2e-10,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-12,
            relative_tolerance: 1e-12,
            linear_restart: 3,
            forcing_minimum: 1e-6,
            ..Default::default()
        },
        audit_entropy: true,
    }
}
fn transport() -> CohesiveLayerHeatSource<()> {
    CohesiveLayerHeatSource {
        conductivity_w_m_k: K,
        area_resistance_m2_k_w: R,
        transport: StepCoefficients {
            sources_w: BTreeMap::new(),
            conductances_w_k: BTreeMap::new(),
            payload: (),
        },
    }
}
fn reference(dt: f64, work: f64, before: [f64; 3]) -> [f64; 3] {
    let g: [f64; 2] = std::array::from_fn(|i| 1.0 / (0.5 / K[i] + 0.05 / K[2] + R[i]));
    let weights: [f64; 2] = std::array::from_fn(|i| dt * g[i] / (CAPACITY[i] + dt * g[i]));
    let layer = (work
        + CAPACITY[2] * before[2]
        + CAPACITY[0] * weights[0] * before[0]
        + CAPACITY[1] * weights[1] * before[1])
        / (CAPACITY[2] + CAPACITY[0] * weights[0] + CAPACITY[1] * weights[1]);
    [
        (1.0 - weights[0]) * before[0] + weights[0] * layer,
        (1.0 - weights[1]) * before[1] + weights[1] * layer,
        layer,
    ]
}

#[test]
fn native_bond_viscosity_retains_real_layer_heat_and_conducts_with_two_resistances() {
    let (source, card) = material_and_card(true, Some(ADHESIVE));
    let elastic = Scalar(1000.0);
    let mass = Scalar(10.0);
    let joint = system(&elastic, &source, &[7], true);
    let motion = compressed(&joint, 0.0)
        .into_dynamics(&mass, vec![0.02], 0.0, dynamics())
        .unwrap();
    let thermal = layer_network();
    let route = ThirdBodyHeatRoute::from_links(&thermal, LAYER, LINKS, K, R).unwrap();
    let dt = 0.01;
    let u = (2.0 * mass.0 * 0.02 / dt) / (2.0 * mass.0 / (dt * dt) + 1000.0 + 4.0 / dt);
    let viscous = 4.0 * u * u / dt;
    let next = propose_deformable_cohesive_layer_heat(
        &motion,
        &[0.0],
        &[[-0.001, 0.0, 0.0]],
        dt,
        binding(&route, &card),
        &thermal,
        layer_controls(),
        &[0.0; 3],
        |state, t, mechanics| {
            close(state.temperature_k(LAYER).unwrap(), t, 1e-12);
            close(mechanics.viscous_dissipated_j, viscous, 1e-12);
            Ok::<_, &'static str>(transport())
        },
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    close(
        next.mechanical.motion.state().displacements_m()[0],
        u,
        1e-12,
    );
    close(next.heat.dissipated_j, viscous, 1e-12);
    let temperatures = reference(dt, viscous, [300.0; 3]);
    for i in 0..3 {
        close(
            next.thermal.temperature_k(CellId(i as u64)).unwrap(),
            temperatures[i],
            2e-10,
        );
    }
    close(
        next.heat.layer_stored_change_j,
        CAPACITY[2] * (temperatures[2] - 300.0),
        1e-11,
    );
    assert!(next.heat.layer_stored_change_j > 0.0);
    assert!(next.heat.layer_stored_change_j < viscous);
    assert!(next.heat.body_heat_j[0] > next.heat.body_heat_j[1]);
    close(
        next.heat.layer_stored_change_j + next.heat.body_heat_j.iter().sum::<f64>(),
        viscous,
        1e-11,
    );
    let entropy = next.heat.entropy.unwrap();
    close(
        entropy.friction_generation_j_k,
        viscous / temperatures[2],
        1e-14,
    );
    assert!(entropy.network_conduction_generation_j_k > 0.0);
    assert_eq!(entropy.interface_conduction_generation_j_k, 0.0);
    assert!(!entropy.has_unclassified_sources);
    assert_eq!(next.fracture_work_j, 0.0);
    close(next.combined_energy_residual_j, 0.0, 2e-10);
    close(next.summed_energy_residual_j, 0.0, 2e-10);
    assert_eq!(motion.time_s(), 0.0);
    assert_eq!(thermal.time_s(), 0.0);

    // The accepted native state and existing finite heat store continue together.
    let again = propose_deformable_cohesive_layer_heat(
        &next.mechanical.motion,
        &[0.0],
        &[[-0.001, 0.0, 0.0]],
        dt,
        binding(&route, &card),
        &next.thermal,
        layer_controls(),
        &[0.0; 3],
        |_, _, _| Ok::<_, &'static str>(transport()),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    let expected = reference(dt, again.mechanical.viscous_dissipated_j, temperatures);
    for i in 0..3 {
        close(
            again.thermal.temperature_k(CellId(i as u64)).unwrap(),
            expected[i],
            2e-10,
        );
    }
    assert_eq!(again.thermal.time_s(), again.mechanical.motion.time_s());
    assert_eq!(again.mechanical.motion.steps(), 2);
}

#[test]
fn progressive_adhesive_damage_heats_only_by_surviving_bond_viscosity() {
    let (source, card) = material_and_card(true, Some(ADHESIVE));
    let elastic = Scalar(1000.0);
    let mass = Scalar(10.0);
    let joint = system(&elastic, &source, &[7], false);
    let initial = joint
        .bonded_state()
        .solve_over_time_controlled(
            &[0.0],
            &[[-0.001, 0.0005, 0.0]],
            1.0,
            dynamics().solve,
            &mut || false,
        )
        .unwrap()
        .state;
    let motion = initial
        .into_dynamics(&mass, vec![0.0], 0.0, dynamics())
        .unwrap();
    let thermal = layer_network();
    let route = ThirdBodyHeatRoute::from_links(&thermal, LAYER, LINKS, K, R).unwrap();
    let next = propose_deformable_cohesive_layer_heat(
        &motion,
        &[0.0],
        &[[-0.001, 0.003, 0.0]],
        0.01,
        binding(&route, &card),
        &thermal,
        layer_controls(),
        &[0.0; 3],
        |_, _, mechanics| {
            close(
                mechanics.motion.state().points()[0].damage(),
                8.0 / 9.0,
                1e-12,
            );
            Ok::<_, &'static str>(transport())
        },
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    let survival_integral = 0.0005 + (4.0 / 3.0) * 0.001 * 3.0_f64.ln() - 0.002 / 3.0;
    let viscous = 4.0 * (0.003 - 0.0005) / 0.01 * survival_integral;
    let bond_work =
        0.5 * 1000.0 * (0.001_f64.powi(2) - 0.0005_f64.powi(2)) + 0.5 * (1.0 + 1.0 / 3.0) * 0.002;
    let storage_change = 0.5 * 1000.0 * ((1.0 / 9.0) * 0.003_f64.powi(2) - 0.0005_f64.powi(2));
    close(next.heat.dissipated_j, viscous, 1e-11);
    close(next.fracture_work_j, bond_work - storage_change, 1e-11);
    close(next.thermal.energy_j() - thermal.energy_j(), viscous, 1e-10);
    close(
        next.mechanical.prescribed_jump_work_j,
        bond_work + viscous,
        1e-10,
    );
    close(next.combined_energy_residual_j, 0.0, 2e-10);
    assert!(next.fracture_work_j > viscous);
    assert_eq!(next.mechanical.motion.velocity_m_s(), &[0.0]);
}

#[test]
fn unlocated_surface_friction_and_material_mismatch_refuse_before_transport() {
    let (source, card) = material_and_card(true, Some(ADHESIVE));
    let elastic = Scalar(1000.0);
    let mass = Scalar(10.0);
    let joint = system(&elastic, &source, &[7], true);
    let motion = compressed(&joint, 0.002)
        .into_dynamics(&mass, vec![0.02], 0.0, dynamics())
        .unwrap()
        .with_damaged_area_friction_controlled(&[Some(1e-4)], &mut || false)
        .unwrap();
    let thermal = layer_network();
    let route = ThirdBodyHeatRoute::from_links(&thermal, LAYER, LINKS, K, R).unwrap();
    let original_point = motion.state().points()[0];
    let cells = thermal.cells().clone();
    let calls = Cell::new(0);
    for failure in 0..5 {
        let mut contact = binding(&route, &card);
        match failure {
            1 => contact.layer_material = MaterialId(0),
            2 => contact.layer_name = "unrelated adhesive",
            3 => contact.bodies.swap(0, 1),
            4 => contact.patch_id = 99,
            _ => (),
        }
        let refused = propose_deformable_cohesive_layer_heat(
            &motion,
            &[0.0],
            &[[-0.001, 0.0, 0.0]],
            0.01,
            contact,
            &thermal,
            layer_controls(),
            &[0.0; 3],
            |_, _, _| {
                calls.set(calls.get() + 1);
                Ok::<_, &'static str>(transport())
            },
            &CancelGate::new_clock_free(),
        );
        assert!(matches!(refused, Err(CohesiveHeatError::Invalid(_))));
    }
    assert_eq!(calls.get(), 0);
    assert_eq!(motion.state().points()[0], original_point);
    assert_eq!(motion.velocity_m_s(), &[0.02]);
    assert_eq!(motion.time_s(), 0.0);
    assert_eq!(thermal.cells(), &cells);
}

#[test]
fn late_transport_failure_and_cancellation_preserve_native_joint_and_layer() {
    let (source, card) = material_and_card(true, Some(ADHESIVE));
    let elastic = Scalar(1000.0);
    let mass = Scalar(10.0);
    let joint = system(&elastic, &source, &[7], true);
    let motion = compressed(&joint, 0.0)
        .into_dynamics(&mass, vec![0.02], 0.0, dynamics())
        .unwrap();
    let thermal = layer_network();
    let route = ThirdBodyHeatRoute::from_links(&thermal, LAYER, LINKS, K, R).unwrap();
    let original_point = motion.state().points()[0];
    let cells = thermal.cells().clone();
    for cancel in [false, true] {
        let calls = Cell::new(0);
        let gate = CancelGate::new_clock_free();
        let refused = propose_deformable_cohesive_layer_heat(
            &motion,
            &[0.0],
            &[[-0.001, 0.0, 0.0]],
            0.01,
            binding(&route, &card),
            &thermal,
            layer_controls(),
            &[0.0; 3],
            |_, _, _| {
                calls.set(calls.get() + 1);
                if calls.get() == 2 {
                    if cancel {
                        gate.request();
                    } else {
                        return Err("outside conductivity source");
                    }
                }
                Ok(transport())
            },
            &gate,
        );
        if cancel {
            assert!(refused.is_err());
        } else {
            assert!(matches!(
                refused,
                Err(CohesiveHeatError::Source("outside conductivity source"))
            ));
        }
        assert!(calls.get() >= 2);
        assert_eq!(motion.state().points()[0], original_point);
        assert_eq!(motion.velocity_m_s(), &[0.02]);
        assert_eq!(motion.time_s(), 0.0);
        assert_eq!(thermal.cells(), &cells);
    }
}
