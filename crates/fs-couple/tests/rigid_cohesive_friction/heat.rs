//! G1/G3/G4: real rigid joint work heats original finite cells once.
use super::*;
use fs_couple::friction_heat::cohesive::CohesiveHeatError;
use fs_couple::friction_heat::rigid_cohesive::*;
use fs_couple::friction_heat::{FrictionHeatRoute, FrictionHeatTolerances};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::*;
use std::cell::Cell;
use std::collections::BTreeMap;

struct Caloric;
impl ThermalLaw for Caloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = 300.0 + h;
        if t.is_finite() && t > 0.0 {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
    fn specific_entropy_change_j_kg_k(
        &self,
        a: f64,
        b: f64,
    ) -> Result<Option<f64>, DepositionError> {
        let ta = self.temperature_k(a)?;
        self.temperature_k(b)?;
        Ok(Some(((b - a) / ta).ln_1p()))
    }
}
static CALORIC: Caloric = Caloric;
const K: [f64; 2] = [2.0, 9.0];
const G: [f64; 2] = [0.004, 0.006];
const CAPACITY: [f64; 2] = [0.01, 0.03];

fn thermal(initial: [f64; 2]) -> DepositionSimulation<'static> {
    let lengths = [0.001, 0.003];
    let mut events = Vec::new();
    for i in 0..2 {
        let (start, end) = if i == 0 {
            (-lengths[i], 0.0)
        } else {
            (0.0, lengths[i])
        };
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64 + 1),
                segment_id: i as u64,
                body_id: [10, 20][i],
                material_id: MaterialId(101 + i as u64),
                mass_kg: CAPACITY[i],
                specific_h_j_kg: initial[i],
                geometry: CellGeometry {
                    start_m: [start, 0.0, 0.0],
                    end_m: [end, 0.0, 0.0],
                    volume_m3: AREA * lengths[i],
                    surface_area_m2: 2.0 * AREA + 4.0 * AREA.sqrt() * lengths[i],
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
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 7,
            a: CellId(1),
            b: CellId(2),
            area_m2: AREA,
            conductance_w_k: series_conductance_w_k(AREA, [0.0005, 0.0015], K).unwrap(),
        }),
    });
    DepositionSimulation::new(
        events,
        &[(MaterialId(101), &CALORIC), (MaterialId(102), &CALORIC)],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
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
fn bodies(source: &ResolvedBilinearCohesive) -> [RigidCohesiveHeatBody<'_>; 2] {
    core::array::from_fn(|i| RigidCohesiveHeatBody {
        body_id: [10, 20][i],
        thermal_material: MaterialId(101 + i as u64),
        mechanical_material: &source.surfaces()[i].material,
    })
}
fn controls() -> CohesiveHeatControls {
    CohesiveHeatControls {
        maximum_cells: 2,
        maximum_interfaces: 1,
        maximum_unknowns: 3,
        maximum_trials: 256,
        maximum_newton_iterations: 24,
        difference_step_j: 1e-6,
        temperature_scale_j_per_k: 1.0,
        heat: FrictionHeatTolerances {
            cell_j: 1e-10,
            thermal_j: 1e-10,
            interface_j: 1e-10,
            total_j: 1e-10,
        },
        combined_energy_tolerance_j: 3e-10,
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
fn inputs(ids: &[u64]) -> CohesiveHeatSource<()> {
    CohesiveHeatSource {
        interfaces: ids
            .iter()
            .map(|&link_id| CohesiveHeatConductivity {
                link_id,
                conductivity_w_m_k: K,
            })
            .collect(),
        transport: StepCoefficients {
            sources_w: BTreeMap::new(),
            conductances_w_k: BTreeMap::new(),
            payload: (),
        },
    }
}
fn contact(route: &FrictionHeatRoute) -> [RigidCohesiveHeatContact<'_>; 1] {
    [RigidCohesiveHeatContact {
        site_index: 0,
        route,
        temperature_guess_k: 305.0,
    }]
}

/// Independent elimination of the two C*dT/dt = g*(Tface-T) equations.
fn thermal_reference(initial: [f64; 2], dt: f64, heat: f64) -> ([f64; 2], f64) {
    let effective =
        core::array::from_fn::<_, 2, _>(|i| G[i] * CAPACITY[i] / (CAPACITY[i] + dt * G[i]));
    let face = (heat / dt + effective[0] * initial[0] + effective[1] * initial[1])
        / (effective[0] + effective[1]);
    (
        core::array::from_fn(|i| {
            (CAPACITY[i] * initial[i] + dt * G[i] * face) / (CAPACITY[i] + dt * G[i])
        }),
        300.0 + face,
    )
}

#[test]
fn native_rigid_viscosity_and_sliding_heat_finite_bodies_with_signed_conduction() {
    let (source, friction) = sources([0.6, 0.3], true);
    let point = partial(&source);
    let dt = 0.01;
    let delta = Vec3::new(0.0, 6e-5, 8e-5);
    let mut owner = system(
        &source,
        Some(&friction),
        point,
        &[[0.0; 3]],
        UnitQuaternion::IDENTITY,
        [Vec3::ZERO, delta.scale(1.0 / dt)],
        zero(),
        dt,
        config(),
    );
    let mut reference = system(
        &source,
        Some(&friction),
        point,
        &[[0.0; 3]],
        UnitQuaternion::IDENTITY,
        [Vec3::ZERO, delta.scale(1.0 / dt)],
        zero(),
        dt,
        config(),
    );
    let initial = [10.0, 0.0];
    let mut heat = thermal(initial);
    let route = FrictionHeatRoute::from_link(&heat, 7, K).unwrap();
    let viscous = 3000.0 * (1.0 - point.damage()) * AREA;
    let load = balanced_load(
        Vec3::new(-0.4, 0.05625, 0.075).add(delta.scale(viscous / dt)),
        v(point.jump_m()).add(delta.scale(0.5)),
    );
    let mechanical_reference = *reference.step(load).unwrap();
    let report = apply_rigid_cohesive_heat(
        &mut owner,
        load,
        &mut heat,
        bodies(&source),
        &contact(&route),
        controls(),
        &[0.0; 2],
        |_, ids, _| Ok::<_, &'static str>(inputs(ids)),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    assert_eq!(owner.bodies(), reference.bodies());
    assert_eq!(report.mechanical, mechanical_reference);
    close_vec(
        v(owner.site_frames()[0].unwrap().jump_m),
        v(point.jump_m()).add(delta),
        1e-10,
    );
    let friction_heat = 0.3 * 0.625 * 0.4 * 1e-4;
    let viscosity_heat = viscous * 1e-8 / dt;
    close(
        report.mechanical.friction_dissipation_j,
        friction_heat,
        1e-11,
    );
    close(
        report.mechanical.viscous_dissipation_j,
        viscosity_heat,
        1e-13,
    );
    close(
        report.heat.dissipated_j,
        friction_heat + viscosity_heat,
        1e-11,
    );
    let (enthalpy, face) = thermal_reference(initial, dt, friction_heat + viscosity_heat);
    for i in 0..2 {
        close(
            heat.cells()[&CellId(i as u64 + 1)].specific_h_j_kg,
            enthalpy[i],
            2e-9,
        );
    }
    close(
        report.heat.interfaces[0].interface_temperature_k,
        face,
        2e-8,
    );
    assert!(report.heat.interfaces[0].body_heat_j[0] < 0.0);
    assert!(report.heat.interfaces[0].body_heat_j[1] > report.heat.dissipated_j);
    assert!(
        report
            .heat
            .entropy
            .as_ref()
            .unwrap()
            .interface_conduction_generation_j_k
            > 0.0
    );
    assert_eq!(report.fracture_work_j, 0.0);
    close(report.combined_energy_residual_j, 0.0, 3e-10);
    close(report.summed_energy_residual_j, 0.0, 3e-10);
    assert_eq!(owner.time_s(), heat.time_s());
    // A second real step retains original source, damage, momenta and clock.
    let original = owner.sites()[0].point();
    let next = v(original.numerical().material().jump_m()).add(delta);
    let (_, mean) = original
        .numerical()
        .material()
        .trial_over_time_with_mean_traction([next.x, next.y, next.z], dt, 2e-7)
        .unwrap();
    let load = balanced_load(
        v(mean.map(|t| t * AREA)).add(Vec3::new(0.0, 0.045, 0.060)),
        v(original.numerical().material().jump_m()).add(delta.scale(0.5)),
    );
    let expected = *reference.step(load).unwrap();
    let report = apply_rigid_cohesive_heat(
        &mut owner,
        load,
        &mut heat,
        bodies(&source),
        &contact(&route),
        controls(),
        &[0.0; 2],
        |_, ids, _| Ok::<_, &'static str>(inputs(ids)),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    assert_eq!(report.mechanical, expected);
    assert_eq!(owner.bodies(), reference.bodies());
    assert_eq!(owner.time_s(), 2.0 * dt);
    assert_eq!(owner.time_s(), heat.time_s());
}

#[test]
fn rigid_surface_creation_remains_separate_from_actual_friction_and_viscous_heat() {
    let (source, friction) = sources([0.6, 0.3], true);
    let point = partial(&source);
    let dt = 0.01;
    let end = [-2e-4, 7.2e-4, 9.6e-4];
    let delta = v(end).sub(v(point.jump_m()));
    let (material, mean) = point
        .trial_over_time_with_mean_traction(end, dt, 2e-7)
        .unwrap();
    // Exact integral of the bilinear damage curve along this shear path:
    // d=.625 until r=.0004, d=1.25-.00025/r to r=.001, then d=1.
    let mean_damage = 1.0 - (0.00025 / 0.0012) * 2.5_f64.ln();
    let capacity = 0.3 * 0.4 * mean_damage;
    let force = v(mean.map(|t| AREA * t)).add(Vec3::new(0.0, 0.6 * capacity, 0.8 * capacity));
    let load = balanced_load(force, v(point.jump_m()).add(delta.scale(0.5)));
    let mut owner = system(
        &source,
        Some(&friction),
        point,
        &[[0.0; 3]],
        UnitQuaternion::IDENTITY,
        [Vec3::ZERO, delta.scale(1.0 / dt)],
        zero(),
        dt,
        config(),
    );
    let mut heat = thermal([0.0; 2]);
    let route = FrictionHeatRoute::from_link(&heat, 7, K).unwrap();
    let report = apply_rigid_cohesive_heat(
        &mut owner,
        load,
        &mut heat,
        bodies(&source),
        &contact(&route),
        controls(),
        &[0.0; 2],
        |_, ids, _| Ok::<_, &'static str>(inputs(ids)),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    assert!(owner.sites()[0].point().numerical().material().failed());
    close(report.fracture_work_j, 75.0 * AREA, 3e-11);
    close(
        report.mechanical.friction_dissipation_j,
        capacity * 0.0012,
        3e-11,
    );
    let physical_heat = capacity * 0.0012 + AREA * material.viscous_dissipated_j_m2;
    close(report.heat.dissipated_j, physical_heat, 3e-11);
    let enthalpy: f64 = heat
        .cells()
        .values()
        .map(|c| c.mass_kg * c.specific_h_j_kg)
        .sum();
    close(enthalpy, physical_heat, 3e-10);
    assert!((enthalpy - physical_heat - report.fracture_work_j).abs() > 5e-5);
    close(report.combined_energy_residual_j, 0.0, 3e-10);
}

#[test]
fn late_thermal_source_and_cancellation_refuse_without_advancing_either_native_owner() {
    let (source, friction) = sources([0.6, 0.3], false);
    let point = failed(&source);
    let mut owner = system(
        &source,
        Some(&friction),
        point,
        &[[0.0; 3]],
        UnitQuaternion::IDENTITY,
        zero(),
        zero(),
        0.002,
        config(),
    );
    let load = balanced_load(Vec3::new(-0.4, 0.16, 0.0), v(point.jump_m()));
    let mut heat = thermal([0.0; 2]);
    let route = FrictionHeatRoute::from_link(&heat, 7, K).unwrap();
    let old_bodies = owner.bodies().clone();
    let old_sites = owner.sites().to_vec();
    let old_friction = owner.damaged_area_friction_states().unwrap().to_vec();
    let old_cells = heat.cells().clone();
    let failed = apply_rigid_cohesive_heat(
        &mut owner,
        load,
        &mut heat,
        bodies(&source),
        &contact(&route),
        controls(),
        &[0.0; 2],
        |_, _, _| Err::<CohesiveHeatSource<()>, _>("conductivity domain"),
        &CancelGate::new_clock_free(),
    );
    assert!(matches!(
        failed,
        Err(RigidCohesiveHeatError::Heat(CohesiveHeatError::Source(
            "conductivity domain"
        )))
    ));
    let gate = CancelGate::new_clock_free();
    let calls = Cell::new(0);
    let cancelled = apply_rigid_cohesive_heat(
        &mut owner,
        load,
        &mut heat,
        bodies(&source),
        &contact(&route),
        controls(),
        &[0.0; 2],
        |_, ids, _| {
            calls.set(calls.get() + 1);
            gate.request();
            Ok::<_, &'static str>(inputs(ids))
        },
        &gate,
    );
    assert!(cancelled.is_err());
    assert!(calls.get() > 0);
    assert_eq!(owner.bodies(), &old_bodies);
    assert_eq!(owner.sites(), old_sites);
    assert_eq!(owner.damaged_area_friction_states().unwrap(), old_friction);
    assert!(owner.last_frame().is_none());
    assert!(owner.site_frames().iter().all(Option::is_none));
    assert_eq!(owner.time_s(), 0.0);
    assert_eq!(heat.time_s(), 0.0);
    assert_eq!(heat.cells(), &old_cells);
    // Stale staged buffers cannot leak into a subsequent successful proposal.
    let accepted = apply_rigid_cohesive_heat(
        &mut owner,
        load,
        &mut heat,
        bodies(&source),
        &contact(&route),
        controls(),
        &[1e-30, -1e-30],
        |_, ids, _| Ok::<_, &'static str>(inputs(ids)),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    assert_eq!(accepted.heat.dissipated_j, 0.0);
    // Successful finite-tolerance thermal solves can retain a numerical energy
    // tail. Refused proposals above must still preserve the exact original cells.
    for (id, cell) in heat.cells() {
        let change_j = cell.mass_kg * (cell.specific_h_j_kg - old_cells[id].specific_h_j_kg);
        close(change_j, 0.0, controls().newton.absolute_tolerance);
    }
    close(
        accepted.heat.thermal.stored_energy_change_j,
        0.0,
        controls().newton.absolute_tolerance,
    );
    close(
        accepted.heat.total_energy_residual_j,
        0.0,
        controls().newton.absolute_tolerance,
    );
    assert_eq!(owner.time_s(), heat.time_s());
}

fn cold_sources(
    card: &InterfaceSystemCard,
    viscous: bool,
) -> (ResolvedBilinearCohesive, ResidualCoulombSource) {
    let query = QueryPoint::new().with("T", 280.0).unwrap();
    let selection = MaterialPropertySelection::SingleClaimOnly;
    let bond = if viscous {
        resolve_viscous_bilinear_cohesive(
            card,
            &query,
            selection.clone(),
            CohesiveEnvelope::CommonMode,
            0.05,
            10.0,
        )
    } else {
        resolve_bilinear_cohesive(card, &query, selection.clone(), 0.05)
    }
    .unwrap();
    let friction = resolve_residual_coulomb(card, &query, selection).unwrap();
    (bond, friction)
}

#[test]
fn prescribed_hot_sources_and_finite_heat_accept_the_same_native_increment() {
    let card = source_card([0.6, 0.3], true, Some([0.3, 0.1]));
    let (bond, residual) = cold_sources(&card, true);
    let hot_bond = bond.at_temperature(&card, "T", 320.0).unwrap();
    let hot_residual = residual.at_temperature(&card, "T", 320.0).unwrap();
    let point = partial(&bond);
    let dt = 0.01;
    let delta = Vec3::new(0.0, 6e-5, 8e-5);
    let make = || {
        system(
            &bond,
            Some(&residual),
            point,
            &[[0.0; 3]],
            UnitQuaternion::IDENTITY,
            [Vec3::ZERO, delta.scale(1.0 / dt)],
            zero(),
            dt,
            config(),
        )
    };
    let mut owner = make();
    let mut reference = make();
    let viscosity = 6000.0 * (1.0 - point.damage()) * AREA;
    let load = balanced_load(
        Vec3::new(-0.6, 0.028125, 0.0375).add(delta.scale(viscosity / dt)),
        v(point.jump_m()).add(delta.scale(0.5)),
    );
    let expected = *reference
        .step_with_temperature_and_friction(load, &[&hot_bond], &[&hot_residual], "T")
        .unwrap();
    let initial = [10.0, 0.0];
    let mut heat = thermal(initial);
    let route = FrictionHeatRoute::from_link(&heat, 7, K).unwrap();
    let report = apply_rigid_cohesive_heat_with_temperature(
        &mut owner,
        load,
        RigidCohesiveHeatTemperature {
            cohesive: &[&hot_bond],
            friction: Some(&[&hot_residual]),
            axis: "T",
        },
        &mut heat,
        bodies(&bond),
        &contact(&route),
        controls(),
        &[0.0; 2],
        |_, ids, _| Ok::<_, &'static str>(inputs(ids)),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    assert_eq!(report.mechanical, expected);
    assert_eq!(owner.bodies(), reference.bodies());
    assert_eq!(owner.site_frames(), reference.site_frames());
    assert_eq!(owner.steps(), 1);
    close(report.mechanical.parameter_work_j, 20e-6, 3e-12);
    let actual_heat = 0.1 * 0.625 * 0.6 * 1e-4 + viscosity * 1e-8 / dt;
    close(report.heat.dissipated_j, actual_heat, 3e-11);
    let (expected_h, expected_face) = thermal_reference(initial, dt, actual_heat);
    for (cell, expected) in heat.cells().values().zip(expected_h) {
        close(cell.specific_h_j_kg, expected, 2e-8);
    }
    close(
        report.heat.interfaces[0].interface_temperature_k,
        expected_face,
        2e-8,
    );
    close(
        report.heat.thermal.stored_energy_change_j,
        actual_heat,
        3e-10,
    );
    close(report.combined_energy_residual_j, 0.0, 3e-10);
    assert!(core::ptr::eq(owner.sites()[0].point().source(), &hot_bond));
    assert!(core::ptr::eq(
        owner.damaged_area_friction_sources().unwrap()[0].source,
        &hot_residual,
    ));
    assert_eq!(owner.time_s(), heat.time_s());
}

#[test]
fn thermal_damage_and_source_work_stay_outside_heat_and_late_refusals_preserve_sources() {
    let card = source_card([0.6, 0.3], false, Some([0.3, 0.1]));
    let (bond, residual) = cold_sources(&card, false);
    let hot_bond = bond.at_temperature(&card, "T", 320.0).unwrap();
    let hot_residual = residual.at_temperature(&card, "T", 320.0).unwrap();
    let jump = [-2e-4, 1.8e-4, 0.0];
    let point = bond.law().bonded_point().trial(jump, 1e-8).unwrap().point;
    let mut owner = system(
        &bond,
        Some(&residual),
        point,
        &[[0.0; 3]],
        UnitQuaternion::IDENTITY,
        zero(),
        zero(),
        0.002,
        config(),
    );
    let mut heat = thermal([0.0; 2]);
    let route = FrictionHeatRoute::from_link(&heat, 7, K).unwrap();
    let damage = (0.00125 / (0.00125 - 0.00016)) * (1.0 - 0.00016 / 0.00018);
    let load = balanced_load(Vec3::new(-0.6, (1.0 - damage) * 0.09 + 0.02, 0.0), v(jump));
    let before_bodies = owner.bodies().clone();
    let before_sites = owner.sites().to_vec();
    let before_friction = owner.damaged_area_friction_states().unwrap().to_vec();
    let before_cells = heat.cells().clone();
    for cancel in [false, true] {
        let gate = CancelGate::new_clock_free();
        let calls = Cell::new(0);
        let refused = apply_rigid_cohesive_heat_with_temperature(
            &mut owner,
            load,
            RigidCohesiveHeatTemperature {
                cohesive: &[&hot_bond],
                friction: Some(&[&hot_residual]),
                axis: "T",
            },
            &mut heat,
            bodies(&bond),
            &contact(&route),
            controls(),
            &[0.0; 2],
            |_, ids, _| {
                calls.set(calls.get() + 1);
                if cancel {
                    gate.request();
                    Ok(inputs(ids))
                } else {
                    Err("outside caloric conductivity domain")
                }
            },
            &gate,
        );
        assert!(refused.is_err());
        assert!(calls.get() > 0);
        assert_eq!(owner.bodies(), &before_bodies);
        assert_eq!(owner.sites(), before_sites);
        assert_eq!(
            owner.damaged_area_friction_states().unwrap(),
            before_friction
        );
        assert!(core::ptr::eq(owner.sites()[0].point().source(), &bond));
        assert!(core::ptr::eq(
            owner.damaged_area_friction_sources().unwrap()[0].source,
            &residual
        ));
        assert_eq!(owner.steps(), 0);
        assert!(owner.last_frame().is_none());
        assert!(owner.site_frames().iter().all(Option::is_none));
        assert_eq!(heat.cells(), &before_cells);
        assert_eq!(heat.time_s(), 0.0);
    }
    // Both successful source directions are still real mechanical increments;
    // the imposed equilibrium loads keep geometry fixed without a projection.
    for (next_bond, next_friction, force, parameter_work, fracture) in [
        (
            &hot_bond,
            &hot_residual,
            Vec3::new(-0.6, (1.0 - damage) * 0.09 + 0.02, 0.0),
            11.9e-6,
            8.1 * damage * AREA,
        ),
        (
            &bond,
            &residual,
            Vec3::new(-0.4, (1.0 - damage) * 0.18 + 0.02, 0.0),
            (-20.0 + 8.1 * (1.0 - damage)) * AREA,
            0.0,
        ),
    ] {
        let report = apply_rigid_cohesive_heat_with_temperature(
            &mut owner,
            balanced_load(force, v(jump)),
            RigidCohesiveHeatTemperature {
                cohesive: &[next_bond],
                friction: Some(&[next_friction]),
                axis: "T",
            },
            &mut heat,
            bodies(&bond),
            &contact(&route),
            controls(),
            &[0.0; 2],
            |_, ids, _| Ok::<_, &'static str>(inputs(ids)),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
        close(report.mechanical.parameter_work_j, parameter_work, 3e-11);
        close(report.fracture_work_j, fracture, 3e-11);
        close(
            report.mechanical.thermal_damage_dissipation_j,
            fracture,
            3e-11,
        );
        close(report.heat.dissipated_j, 0.0, 1e-15);
        close(report.heat.thermal.stored_energy_change_j, 0.0, 1e-12);
        close(report.combined_energy_residual_j, 0.0, 3e-10);
        close(
            owner.sites()[0].point().numerical().material().damage(),
            damage,
            3e-10,
        );
        assert!(core::ptr::eq(owner.sites()[0].point().source(), next_bond));
        assert_eq!(owner.time_s(), heat.time_s());
    }
    assert_eq!(owner.steps(), 2);
}

#[test]
fn ordinary_native_heat_continues_from_an_interval_between_nominal_ticks() {
    let (source, friction) = sources([0.6, 0.3], false);
    let point = failed(&source);
    let mut owner = system(
        &source,
        Some(&friction),
        point,
        &[[0.0; 3]],
        UnitQuaternion::IDENTITY,
        zero(),
        zero(),
        0.002,
        config(),
    );
    let load = balanced_load(Vec3::new(-0.4, 0.16, 0.0), v(point.jump_m()));
    owner
        .advance_interval(
            load,
            0.0037,
            RigidCohesiveIntervalConfig {
                minimum_step_s: 1e-12,
                reversal_time_tolerance_s: 1e-10,
                maximum_trials: 16,
                maximum_substeps: 16,
                maximum_total_evaluations: 1_000_000,
                maximum_total_material_evaluations: 2_000_000,
                maximum_total_iterations: 1024,
            },
        )
        .unwrap();
    let mut heat = thermal([0.0; 2])
        .advance_to(
            owner.time_s(),
            StepConfig {
                max_dt_s: 0.002,
                energy_tolerance_j: 1e-10,
            },
            || Ok(()),
        )
        .unwrap()
        .0;
    let route = FrictionHeatRoute::from_link(&heat, 7, K).unwrap();
    let target = owner.next_step_time_s().unwrap();
    assert_eq!(target, 0.0037 + 0.002);
    assert_ne!(target, (owner.steps() + 1) as f64 * owner.step_duration_s());
    let report = apply_rigid_cohesive_heat(
        &mut owner,
        load,
        &mut heat,
        bodies(&source),
        &contact(&route),
        controls(),
        &[0.0; 2],
        |_, ids, _| Ok::<_, &'static str>(inputs(ids)),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    assert_eq!(owner.time_s(), target);
    assert_eq!(heat.time_s(), target);
    assert_eq!(owner.steps(), 3);
    assert_eq!(report.heat.dissipated_j, 0.0);
    close(report.combined_energy_residual_j, 0.0, 3e-10);
}
