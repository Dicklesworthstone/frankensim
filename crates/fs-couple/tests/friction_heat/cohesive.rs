//! Native joint loss and finite-body heat, with independent momentum/caloric
//! references. Constants below are analytical fixtures, not material estimates.
use super::*;
use fs_couple::friction_heat::cohesive::*;
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SurfaceSpec, SystemContext,
    UncertaintyModel,
};
use fs_material::cohesive::source::*;
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    interface_stiffness_per_area_dims, resolve_interface_state_point,
};
use fs_qty::{Dims, Pressure, Temperature};
use fs_solid::cohesive::dynamics::DynamicsConfig;
use fs_solid::cohesive::dynamics::coulomb::{CoulombActivation, CoulombConfig, CoulombRegime};
use fs_solid::cohesive::{Budget, CohesivePatch, CohesiveState, CohesiveSystem, Config, JumpTerm};
use fs_solver::{LinearOp, NewtonKrylovConfig};
use std::cell::Cell;

#[path = "cohesive/layer.rs"]
mod layer;

#[path = "cohesive/interval.rs"]
mod interval;

fn material(viscous: bool) -> ResolvedBilinearCohesive {
    material_and_card(viscous, None).0
}

fn material_and_card(
    viscous: bool,
    third_body: Option<&str>,
) -> (ResolvedBilinearCohesive, InterfaceSystemCard) {
    let mut rows = vec![
        (
            COHESIVE_STIFFNESS_PROPERTY,
            interface_stiffness_per_area_dims(),
            1000.0,
        ),
        (COHESIVE_STRENGTH_PROPERTY, Pressure::DIMS, 1.0),
        (
            COHESIVE_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
            0.002,
        ),
        (
            COHESIVE_COMPRESSION_STIFFNESS_PROPERTY,
            interface_stiffness_per_area_dims(),
            2000.0,
        ),
        (RESIDUAL_KINETIC_FRICTION_PROPERTY, Dims::NONE, 0.3),
        (RESIDUAL_STATIC_FRICTION_PROPERTY, Dims::NONE, 0.5),
    ];
    if viscous {
        rows.extend([
            (
                COHESIVE_NORMAL_VISCOSITY_PROPERTY,
                COHESIVE_VISCOSITY_DIMS,
                2.0,
            ),
            (
                COHESIVE_TANGENTIAL_VISCOSITY_PROPERTY,
                COHESIVE_VISCOSITY_DIMS,
                4.0,
            ),
        ]);
    }
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, dims, value) in rows {
        requirements.push(
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap(),
        );
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(280.0, value), (400.0, value)],
                    dims,
                },
                validity: ValidityDomain::unconstrained().with("T", 280.0, 400.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured deformable heat fixture; no specimen claim".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let surface = |name: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "original finite member".into(),
            revision: 1,
        },
        texture_frame: "fixed normal/shear axes".into(),
    };
    let card = InterfaceSystemCard::assemble(
        surface("A"),
        surface("B"),
        SystemContext {
            medium: "original bond and residual surfaces".into(),
            third_body: third_body.map(str::to_owned),
            environment: "manufactured constant coefficients over admitted temperature".into(),
            history: "initially bonded".into(),
        },
        claims,
        vec![],
    )
    .unwrap();
    let resolved = resolve_interface_state_point(
        &card,
        &QueryPoint::new().with("T", 300.0).unwrap(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let source = if viscous {
        ResolvedBilinearCohesive::try_from_viscous_resolved(
            &card,
            resolved,
            CohesiveEnvelope::CommonMode,
            0.02,
            1.0,
        )
        .unwrap()
    } else {
        ResolvedBilinearCohesive::try_from_resolved(&card, resolved, 0.02).unwrap()
    };
    (source, card)
}

struct Scalar(f64);
impl LinearOp for Scalar {
    fn n(&self) -> usize {
        1
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        y[0] = self.0 * x[0];
    }
}
fn system<'a>(
    elastic: &'a dyn LinearOp,
    source: &'a ResolvedBilinearCohesive,
    ids: &[u64],
    free: bool,
) -> CohesiveSystem<'a> {
    let patches = ids
        .iter()
        .map(|&id| {
            CohesivePatch::from_interface(
                id,
                1.0,
                source,
                [
                    &source.surfaces()[0].material,
                    &source.surfaces()[1].material,
                ],
                if free {
                    vec![JumpTerm {
                        dof: 0,
                        local: [0.0, 1.0, 0.0],
                    }]
                } else {
                    vec![]
                },
            )
            .unwrap()
        })
        .collect();
    CohesiveSystem::new(
        elastic,
        patches,
        Budget {
            max_dofs: 1,
            max_patches: ids.len(),
            max_jump_terms: ids.len(),
            max_krylov_basis_entries: 128,
        },
    )
    .unwrap()
}
fn dynamics() -> DynamicsConfig {
    DynamicsConfig {
        solve: Config {
            newton: NewtonKrylovConfig {
                absolute_tolerance: 1e-10,
                relative_tolerance: 1e-12,
                linear_restart: 3,
                forcing_maximum: 0.01,
                ..Default::default()
            },
            max_newton_iterations: 40,
            max_operator_evaluations: 2000,
            force_tolerance_n: 1e-8,
            local_energy_tolerance_j_m2: 1e-10,
            total_cohesive_energy_tolerance_j: 1e-10,
            maximum_displacement_m: 0.02,
        },
        difference_step_m: 1e-8,
        maximum_step_s: 0.02,
        maximum_velocity_m_s: 2.0,
        maximum_steps: 50,
        absolute_energy_tolerance_j: 1e-10,
        relative_energy_tolerance: 1e-8,
    }
}
fn compressed<'s, 'a>(system: &'s CohesiveSystem<'a>, opening: f64) -> CohesiveState<'s, 'a> {
    let n = system.patches().len();
    system
        .bonded_state()
        .solve_over_time_controlled(
            &[0.0],
            &vec![[opening, 0.0, 0.0]; n],
            1.0,
            dynamics().solve,
            &mut || false,
        )
        .unwrap()
        .state
        .solve_over_time_controlled(
            &[0.0],
            &vec![[-0.001, 0.0, 0.0]; n],
            1.0,
            dynamics().solve,
            &mut || false,
        )
        .unwrap()
        .state
}
fn controls(cells: usize, interfaces: usize) -> CohesiveHeatControls {
    CohesiveHeatControls {
        maximum_cells: cells,
        maximum_interfaces: interfaces,
        maximum_unknowns: cells + interfaces,
        maximum_trials: 512,
        maximum_newton_iterations: 30,
        difference_step_j: 1e-6,
        temperature_scale_j_per_k: 1.0,
        heat: tolerances(),
        combined_energy_tolerance_j: 2e-10,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-12,
            relative_tolerance: 1e-12,
            linear_restart: cells + interfaces,
            forcing_minimum: 1e-6,
            ..Default::default()
        },
        audit_entropy: true,
    }
}
fn inputs(ids: &[u64], k: impl Fn(u64) -> [f64; 2]) -> CohesiveHeatSource<()> {
    CohesiveHeatSource {
        interfaces: ids
            .iter()
            .map(|&link_id| CohesiveHeatConductivity {
                link_id,
                conductivity_w_m_k: k(link_id),
            })
            .collect(),
        transport: StepCoefficients {
            sources_w: BTreeMap::new(),
            conductances_w_k: BTreeMap::new(),
            payload: (),
        },
    }
}

#[test]
fn partial_joint_deceleration_heats_finite_bodies_by_actual_viscosity_and_friction() {
    let source = material(true);
    let elastic = Scalar(1000.0);
    let mass = Scalar(10.0);
    let joint = system(&elastic, &source, &[7], true);
    let initial = compressed(&joint, 0.002);
    close(initial.points()[0].damage(), 2.0 / 3.0, 1e-13);
    let motion = initial
        .into_dynamics(&mass, vec![0.02], 0.0, dynamics())
        .unwrap()
        .with_damaged_area_friction_controlled(&[Some(1e-4)], &mut || false)
        .unwrap();
    let thermal = network(
        1,
        [1.0, 3.0],
        [2.0, 9.0],
        [0.5, 2.0],
        [0.0; 2],
        [40.0, 10.0],
    );
    let route = FrictionHeatRoute::from_link(&thermal, CONTACT, [2.0, 9.0]).unwrap();
    let dt = 0.01;
    let resistance = 0.3 * (2.0 / 3.0) * 2000.0 * 0.001;
    let viscosity = 4.0 / 3.0;
    let u = (2.0 * mass.0 * 0.02 / dt - resistance)
        / (2.0 * mass.0 / (dt * dt) + 0.5 * (1000.0 + 1000.0 / 3.0) + viscosity / dt);
    let heat = resistance * u + viscosity * u * u / dt;
    let mut reference = factor(&thermal, dt, [4.0, 6.0], 0.0);
    let (energy, temperature) = endpoint(&thermal, &mut reference, dt, heat / dt);
    let next = propose_deformable_cohesive_heat(
        &motion,
        &[0.0],
        &[[-0.001, 0.0, 0.0]],
        dt,
        &[CohesiveHeatContact {
            patch_id: 7,
            route: &route,
            temperature_guess_k: 320.0,
        }],
        &thermal,
        controls(2, 1),
        &[0.0; 2],
        |_, ids, _| Ok::<_, &'static str>(inputs(ids, |_| [2.0, 9.0])),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    close(
        next.mechanical.motion.state().displacements_m()[0],
        u,
        1e-12,
    );
    close(
        next.mechanical.motion.velocity_m_s()[0],
        2.0 * u / dt - 0.02,
        1e-10,
    );
    close(next.mechanical.friction_dissipated_j, resistance * u, 1e-12);
    close(
        next.mechanical.viscous_dissipated_j,
        viscosity * u * u / dt,
        1e-12,
    );
    close(next.heat.dissipated_j, heat, 1e-12);
    for (i, cell) in thermal.cells().values().enumerate() {
        close(
            next.thermal.cells()[&cell.id].specific_h_j_kg,
            cell.specific_h_j_kg + energy[i] / cell.mass_kg,
            1e-10,
        );
    }
    close(
        next.heat.interfaces[0].interface_temperature_k,
        temperature,
        1e-10,
    );
    assert!(next.heat.interfaces[0].body_heat_j[0] < 0.0);
    assert!(next.heat.interfaces[0].body_heat_j[1] > heat);
    assert!(
        next.heat
            .entropy
            .as_ref()
            .unwrap()
            .interface_conduction_generation_j_k
            > 0.0
    );
    assert_eq!(next.fracture_work_j, 0.0);
    close(next.combined_energy_residual_j, 0.0, 2e-10);
    close(next.summed_energy_residual_j, 0.0, 2e-10);
    assert_eq!(motion.time_s(), 0.0);
    assert_eq!(thermal.time_s(), 0.0);
    // Continue with both actual endpoints, retaining original damage and clocks.
    let again = propose_deformable_cohesive_heat(
        &next.mechanical.motion,
        &[0.0],
        &[[-0.001, 0.0, 0.0]],
        dt,
        &[CohesiveHeatContact {
            patch_id: 7,
            route: &route,
            temperature_guess_k: temperature,
        }],
        &next.thermal,
        controls(2, 1),
        &[0.0; 2],
        |_, ids, _| Ok::<_, &'static str>(inputs(ids, |_| [2.0, 9.0])),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    assert_eq!(again.mechanical.motion.time_s(), again.thermal.time_s());
    assert_eq!(again.mechanical.motion.steps(), 2);
    close(
        again.mechanical.motion.state().points()[0].damage(),
        2.0 / 3.0,
        1e-13,
    );
    close(again.combined_energy_residual_j, 0.0, 2e-10);
}

#[test]
fn prescribed_damage_growth_keeps_fracture_work_out_of_heat() {
    let source = material(false);
    let elastic = Scalar(1000.0);
    let mass = Scalar(10.0);
    let joint = system(&elastic, &source, &[7], false);
    let start = [-0.001, 0.0005, 0.0];
    let initial = joint
        .bonded_state()
        .solve_controlled(&[0.0], &[start], dynamics().solve, &mut || false)
        .unwrap()
        .state;
    let motion = initial
        .into_dynamics(&mass, vec![0.0], 0.0, dynamics())
        .unwrap()
        .with_damaged_area_friction_controlled(&[Some(1e-4)], &mut || false)
        .unwrap();
    let thermal = network(1, [1.0; 2], [1.0; 2], [1.0; 2], [0.0; 2], [0.0; 2]);
    let route = FrictionHeatRoute::from_link(&thermal, CONTACT, [1.0; 2]).unwrap();
    let next = propose_deformable_cohesive_heat(
        &motion,
        &[0.0],
        &[[-0.001, 0.003, 0.0]],
        0.01,
        &[CohesiveHeatContact {
            patch_id: 7,
            route: &route,
            temperature_guess_k: 300.0,
        }],
        &thermal,
        controls(2, 1),
        &[0.0; 2],
        |_, ids, _| Ok::<_, &'static str>(inputs(ids, |_| [1.0; 2])),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    // Analytic integral of d(q)=(4/3)(1-.001/q), q=.001..003.
    let heat = 0.3 * 2000.0 * 0.001 * (4.0 / 3.0) * (0.002 - 0.001 * 3.0_f64.ln());
    let bond_work =
        0.5 * 1000.0 * (0.001_f64.powi(2) - 0.0005_f64.powi(2)) + 0.5 * (1.0 + 1.0 / 3.0) * 0.002;
    let storage_change = 0.5 * 1000.0 * ((1.0 / 9.0) * 0.003_f64.powi(2) - 0.0005_f64.powi(2));
    close(next.fracture_work_j, bond_work - storage_change, 1e-12);
    assert!(next.fracture_work_j > heat);
    close(next.thermal.energy_j() - thermal.energy_j(), heat, 1e-10);
    close(next.heat.dissipated_j, heat, 1e-12);
    close(
        next.mechanical.prescribed_jump_work_j,
        bond_work + heat,
        1e-12,
    );
    close(next.combined_energy_residual_j, 0.0, 2e-10);
    assert_eq!(next.mechanical.motion.velocity_m_s(), &[0.0]);
}

fn three_cells() -> DepositionSimulation<'static> {
    let mut events = Vec::new();
    for i in 0..3 {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64 + 1),
                segment_id: i as u64,
                body_id: i as u64,
                material_id: MaterialId(1),
                mass_kg: [0.5, 2.0, 1.5][i],
                specific_h_j_kg: [10.0, 20.0, 30.0][i],
                geometry: CellGeometry {
                    start_m: [i as f64, 0.0, 0.0],
                    end_m: [i as f64 + 1.0, 0.0, 0.0],
                    volume_m3: 1.0,
                    surface_area_m2: 6.0,
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
    for (id, a, b, k) in [(20, 1, 2, [1.0, 2.0]), (10, 2, 3, [2.0, 3.0])] {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id,
                a: CellId(a),
                b: CellId(b),
                area_m2: 1.0,
                conductance_w_k: series_conductance_w_k(1.0, [0.5; 2], k).unwrap(),
            }),
        });
    }
    DepositionSimulation::new(
        events,
        &[(MaterialId(1), &SENSIBLE as &dyn ThermalLaw)],
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
            energy_tolerance_j: 1e-10,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

#[test]
fn two_original_patch_losses_share_one_finite_cell_without_duplication_or_order_dependence() {
    let source = material(false);
    let elastic = Scalar(1000.0);
    let mass = Scalar(10.0);
    let thermal = three_cells();
    let routes = [
        FrictionHeatRoute::from_link(&thermal, 20, [1.0, 2.0]).unwrap(),
        FrictionHeatRoute::from_link(&thermal, 10, [2.0, 3.0]).unwrap(),
    ];
    let slips = [0.0002, 0.0007];
    let heats = slips.map(|s| 0.3 * 2000.0 * 0.001 * s);
    let dt = 0.01;
    // Independent 3-body/2-face temperature balance with deliberately reversed
    // face IDs. No call to the production heat residual or its reported losses.
    let mut matrix = vec![0.0; 25];
    let mut rhs = vec![0.0; 5];
    for (i, cell) in thermal.cells().values().enumerate() {
        matrix[i * 5 + i] = cell.mass_kg / dt;
        rhs[i] = cell.mass_kg / dt * cell.specific_h_j_kg;
    }
    for (p, (a, b, g)) in [(0usize, 1usize, [2.0, 4.0]), (1, 2, [4.0, 6.0])]
        .into_iter()
        .enumerate()
    {
        let f = 3 + p;
        rhs[f] = heats[p] / dt;
        for (i, conductance) in [a, b].into_iter().zip(g) {
            matrix[i * 5 + i] += conductance;
            matrix[i * 5 + f] -= conductance;
            matrix[f * 5 + i] -= conductance;
            matrix[f * 5 + f] += conductance;
        }
    }
    let mut lu = LuWorkspace::new(5).unwrap();
    lu.factor(&matrix).unwrap();
    let mut reference = vec![0.0; 5];
    lu.solve_factored_into(&rhs, &mut reference).unwrap();
    let mut endpoints = Vec::new();
    for reversed in [false, true] {
        let ids = if reversed { [9, 7] } else { [7, 9] };
        let joint = system(&elastic, &source, &ids, false);
        let motion = compressed(&joint, 0.005)
            .into_dynamics(&mass, vec![0.0], 0.0, dynamics())
            .unwrap()
            .with_damaged_area_friction_controlled(&[Some(1e-4); 2], &mut || false)
            .unwrap();
        let offsets: Vec<_> = ids
            .iter()
            .map(|id| [-0.001, slips[usize::from(*id == 9)], 0.0])
            .collect();
        let mut contacts = vec![
            CohesiveHeatContact {
                patch_id: 7,
                route: &routes[0],
                temperature_guess_k: 315.0,
            },
            CohesiveHeatContact {
                patch_id: 9,
                route: &routes[1],
                temperature_guess_k: 325.0,
            },
        ];
        if reversed {
            contacts.reverse();
        }
        let next = propose_deformable_cohesive_heat(
            &motion,
            &[0.0],
            &offsets,
            dt,
            &contacts,
            &thermal,
            controls(3, 2),
            &[0.0; 3],
            |_, ids, _| {
                Ok::<_, &'static str>(inputs(
                    ids,
                    |id| if id == 10 { [2.0, 3.0] } else { [1.0, 2.0] },
                ))
            },
            &CancelGate::new_clock_free(),
        )
        .unwrap();
        assert_eq!(next.patch_ids, [9, 7]);
        close(next.heat.interfaces[0].dissipated_j, heats[1], 1e-12);
        close(next.heat.interfaces[1].dissipated_j, heats[0], 1e-12);
        close(next.heat.dissipated_j, heats.iter().sum(), 1e-12);
        for (i, cell) in next.thermal.cells().values().enumerate() {
            close(cell.specific_h_j_kg, reference[i], 1e-10);
        }
        close(
            next.heat.interfaces[0].interface_temperature_k,
            300.0 + reference[4],
            1e-10,
        );
        close(
            next.heat.interfaces[1].interface_temperature_k,
            300.0 + reference[3],
            1e-10,
        );
        close(next.combined_energy_residual_j, 0.0, 2e-10);
        endpoints.push(next.thermal.cells().clone());
    }
    assert_eq!(endpoints[0], endpoints[1]);
}

#[test]
fn static_joint_has_no_artificial_heat_and_late_refusal_keeps_both_original_states() {
    let source = material(false);
    let elastic = Scalar(1000.0);
    let mass = Scalar(10.0);
    let joint = system(&elastic, &source, &[7], true);
    let motion = compressed(&joint, 0.005)
        .into_dynamics(&mass, vec![0.0], 0.0, dynamics())
        .unwrap()
        .with_damaged_area_coulomb_controlled(
            &[Some(CoulombActivation::AtRest)],
            CoulombConfig {
                traction_scale_pa_per_m: 1000.0,
                slip_tolerance_m: 1e-11,
                velocity_tolerance_m_s: 1e-8,
            },
            &mut || false,
        )
        .unwrap();
    let thermal = network(1, [1.0; 2], [1.0; 2], [1.0; 2], [0.0; 2], [10.0; 2]);
    let route = FrictionHeatRoute::from_link(&thermal, CONTACT, [1.0; 2]).unwrap();
    let contacts = [CohesiveHeatContact {
        patch_id: 7,
        route: &route,
        temperature_guess_k: 310.0,
    }];
    let next = propose_deformable_cohesive_heat(
        &motion,
        &[0.4],
        &[[-0.001, 0.0, 0.0]],
        0.01,
        &contacts,
        &thermal,
        controls(2, 1),
        &[1e-30, -1e-30],
        |_, ids, _| Ok::<_, &'static str>(inputs(ids, |_| [1.0; 2])),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    assert_eq!(
        next.mechanical
            .motion
            .damaged_area_coulomb_states()
            .unwrap()[0]
            .occupancy,
        CoulombRegime::Sticking
    );
    assert_eq!(next.heat.dissipated_j, 0.0);
    assert_eq!(next.mechanical.friction_dissipated_j, 0.0);
    assert_eq!(next.thermal.cells(), thermal.cells());
    let original_point = motion.state().points()[0];
    let original_cells = thermal.cells().clone();
    let calls = Cell::new(0usize);
    let gate = CancelGate::new_clock_free();
    let result = propose_deformable_cohesive_heat(
        &motion,
        &[0.4],
        &[[-0.001, 0.0, 0.0]],
        0.01,
        &contacts,
        &thermal,
        controls(2, 1),
        &[0.0; 2],
        |_, ids, _| {
            calls.set(calls.get() + 1);
            gate.request();
            Ok::<_, &'static str>(inputs(ids, |_| [1.0; 2]))
        },
        &gate,
    );
    assert!(result.is_err());
    assert!(calls.get() > 0); // cancellation after the real mechanical solve
    let missing = propose_deformable_cohesive_heat(
        &motion,
        &[0.4],
        &[[-0.001, 0.0, 0.0]],
        0.01,
        &contacts,
        &thermal,
        controls(2, 1),
        &[0.0; 2],
        |_, _, _| Err::<CohesiveHeatSource<()>, _>("source domain"),
        &CancelGate::new_clock_free(),
    );
    assert!(matches!(
        missing,
        Err(CohesiveHeatError::Source("source domain"))
    ));
    let opened = propose_deformable_cohesive_heat(
        &motion,
        &[0.4],
        &[[0.0001, 0.0, 0.0]],
        0.01,
        &contacts,
        &thermal,
        controls(2, 1),
        &[0.0; 2],
        |_, ids, _| Ok::<_, &'static str>(inputs(ids, |_| [1.0; 2])),
        &CancelGate::new_clock_free(),
    );
    assert!(opened.is_err());
    assert_eq!(motion.state().points()[0], original_point);
    assert_eq!(motion.state().displacements_m(), &[0.0]);
    assert_eq!(motion.velocity_m_s(), &[0.0]);
    assert_eq!(motion.time_s(), 0.0);
    assert_eq!(thermal.cells(), &original_cells);
    assert_eq!(thermal.time_s(), 0.0);
}
