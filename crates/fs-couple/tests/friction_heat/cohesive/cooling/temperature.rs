//! G1/G3/G4: actual heat changes the SAME deformable and thermal bond history.
use super::*;
use fs_couple::friction_heat::cohesive::cooling::temperature::*;

fn thermal_material(mixed: bool, viscous: bool) -> (ResolvedBilinearCohesive, InterfaceSystemCard) {
    let (_, base) = material_and_card(viscous, None);
    let stiffness = interface_stiffness_per_area_dims();
    let mut rows = if mixed {
        vec![
            (
                COHESIVE_NORMAL_STIFFNESS_PROPERTY,
                stiffness,
                [1000.0, 400.0],
            ),
            (
                COHESIVE_TANGENTIAL_STIFFNESS_PROPERTY,
                stiffness,
                [500.0, 200.0],
            ),
            (
                COHESIVE_NORMAL_STRENGTH_PROPERTY,
                Pressure::DIMS,
                [1.0, 0.1],
            ),
            (
                COHESIVE_TANGENTIAL_STRENGTH_PROPERTY,
                Pressure::DIMS,
                [1.5, 0.1],
            ),
            (
                COHESIVE_MODE_I_FRACTURE_ENERGY_PROPERTY,
                FRACTURE_ENERGY_PER_AREA_DIMS,
                [0.002, 0.00002],
            ),
            (
                COHESIVE_MODE_II_FRACTURE_ENERGY_PROPERTY,
                FRACTURE_ENERGY_PER_AREA_DIMS,
                [0.006, 0.00006],
            ),
            (COHESIVE_BK_EXPONENT_PROPERTY, Dims::NONE, [2.0; 2]),
        ]
    } else {
        vec![
            (COHESIVE_STIFFNESS_PROPERTY, stiffness, [1000.0, 400.0]),
            (COHESIVE_STRENGTH_PROPERTY, Pressure::DIMS, [1.0, 0.1]),
            (
                COHESIVE_FRACTURE_ENERGY_PROPERTY,
                FRACTURE_ENERGY_PER_AREA_DIMS,
                [0.002, 0.00002],
            ),
        ]
    };
    rows.push((
        COHESIVE_COMPRESSION_STIFFNESS_PROPERTY,
        stiffness,
        [2000.0; 2],
    ));
    if viscous {
        rows.extend([
            (
                COHESIVE_NORMAL_VISCOSITY_PROPERTY,
                COHESIVE_VISCOSITY_DIMS,
                [2.0, 0.5],
            ),
            (
                COHESIVE_TANGENTIAL_VISCOSITY_PROPERTY,
                COHESIVE_VISCOSITY_DIMS,
                [4.0, 1.0],
            ),
        ]);
    }
    let mut claims = ClaimSet::new();
    for (name, dims, values) in rows {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(300.0, values[0]), (400.0, values[1])],
                    dims,
                },
                validity: ValidityDomain::unconstrained().with("T", 300.0, 400.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "analytical thermal weakening; no specimen calibration".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = InterfaceSystemCard::assemble(
        base.surface_a().clone(),
        base.surface_b().clone(),
        base.context().clone(),
        claims,
        vec![],
    )
    .unwrap();
    let query = QueryPoint::new().with("T", 300.0).unwrap();
    let select = MaterialPropertySelection::SingleClaimOnly;
    let source = if viscous {
        resolve_viscous_bilinear_cohesive(
            &card,
            &query,
            select,
            if mixed {
                CohesiveEnvelope::FrozenOnsetMixedMode
            } else {
                CohesiveEnvelope::CommonMode
            },
            0.02,
            1.0,
        )
        .unwrap()
    } else if mixed {
        resolve_mixed_mode_bilinear_cohesive(&card, &query, select, 0.02).unwrap()
    } else {
        resolve_bilinear_cohesive(&card, &query, select, 0.02).unwrap()
    };
    (source, card)
}

fn binding(card: &InterfaceSystemCard, weight: f64) -> [CohesiveCoolingTemperatureBinding<'_>; 1] {
    [CohesiveCoolingTemperatureBinding {
        patch_id: CONTACT,
        card,
        temperature_weight_a: weight,
    }]
}
fn controls() -> CohesiveCoolingTemperatureControls {
    CohesiveCoolingTemperatureControls {
        local_energy_tolerance_j_m2: 1e-10,
        total_energy_tolerance_j: 1e-10,
    }
}
fn tracked<'a>(
    thermal: DepositionSimulation<'static>,
    point: CohesivePoint<'a>,
    open_gap: f64,
) -> FractureCooling<'static, 'a> {
    let reference = thermal.links()[&CONTACT];
    FractureCooling::new(
        thermal,
        vec![Interface {
            reference,
            point,
            closed: true,
            policy: ContactPolicy {
                open_gap_m: open_gap,
                reclose_gap_m: 0.0,
                reclosed_conductance_w_k: 0.25,
                maximum_slip_m: 0.02,
            },
        }],
        1,
        || Ok(()),
    )
    .unwrap()
}

fn heat_cycle() -> DepositionSimulation<'static> {
    let initial = network(1, [1.0; 2], [1.0; 2], [1.0; 2], [0.0; 2], [0.0; 2]);
    let mut events: Vec<_> = initial
        .cells()
        .values()
        .cloned()
        .map(|mut cell| {
            cell.boundary.source_w = 10_000.0;
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(cell),
            }
        })
        .collect();
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(initial.links()[&CONTACT]),
    });
    for (time_s, source_w) in [(0.01, 0.0), (0.02, -10_000.0), (0.03, 0.0)] {
        for cell in initial.cells().values() {
            let mut boundary = cell.boundary;
            boundary.source_w = source_w;
            events.push(TimedEvent {
                time_s,
                event: Event::Boundary(cell.id, boundary),
            });
        }
    }
    DepositionSimulation::new(
        events,
        &[(MaterialId(1), &SENSIBLE)],
        DepositionBudget {
            max_cells: 2,
            max_links: 2,
            max_events: 10,
            max_substeps: 100,
        },
    )
    .unwrap()
    .advance_to(0.0, cooling_controls().thermal, || Ok(()))
    .unwrap()
    .0
}

#[test]
fn actual_heating_weakens_native_motion_and_heat_contact_then_cooling_cannot_rebond() {
    let (cold, card) = thermal_material(false, false);
    let elastic = Scalar(500.0);
    let mass = Scalar(1.0);
    let joint = normal_system(&elastic, &cold);
    let force = 0.75;
    let state = joint
        .bonded_state()
        .solve_controlled(&[force], &[[0.0; 3]], dynamics().solve, &mut || false)
        .unwrap()
        .state;
    let motion = state
        .into_dynamics(&mass, vec![0.0], 0.0, dynamics())
        .unwrap();
    let original = tracked(heat_cycle(), motion.state().points()[0], 0.00045);
    let gate = CancelGate::new_clock_free();
    let heated = propose_deformable_cohesive_cooling(
        &motion,
        &original,
        &[force],
        &[[0.0; 3]],
        0.01,
        &partition(0.0, 0.0),
        cooling_controls(),
        &gate,
    )
    .unwrap();
    close(heated.transport.external_heat_j, 200.0, 1e-10);
    assert_eq!(
        heated.cooling.thermal().temperature_k(CellId(1)).unwrap(),
        400.0
    );
    assert!(!heated.mechanical.motion.state().points()[0].failed());
    let inputs = binding(&card, 0.021);
    let sources = prepare_deformable_cohesive_temperature_sources(
        &heated.mechanical.motion,
        &heated.cooling,
        &inputs,
        "T",
        &gate,
    )
    .unwrap();
    assert_eq!(sources[0].resolved().state_coordinate("T"), Some(400.0));
    let target = joint
        .with_temperature_sources_controlled(&[&sources[0]], "T", &mut || false)
        .unwrap();
    let updated = propose_deformable_cohesive_temperature(
        &heated.mechanical.motion,
        &heated.cooling,
        &target,
        &inputs,
        "T",
        controls(),
        &gate,
    )
    .unwrap();

    let gap = heated.mechanical.motion.state().points()[0].jump_m()[0];
    close(gap, 0.0005, 1e-12);
    let old_storage = 0.5 * 1000.0 * gap * gap;
    let changed_storage = 0.5 * 400.0 * gap * gap;
    close(
        updated.mechanical.energy.parameter_work_j,
        changed_storage - old_storage,
        1e-12,
    );
    close(
        updated.mechanical.energy.damage_dissipated_j,
        changed_storage,
        1e-12,
    );
    close(
        updated.mechanical.energy.stored_change_j,
        -old_storage,
        1e-12,
    );
    close(updated.combined_energy_residual_j, 0.0, 1e-12);
    assert!(updated.mechanical.motion.state().points()[0].failed());
    assert_eq!(
        updated.mechanical.motion.state().points()[0],
        updated.cooling.interfaces()[&CONTACT].point
    );
    assert_eq!(updated.contact.contact.opened, 1);
    assert_eq!(updated.contact.contact.viscous_dissipation_j, 0.0);
    assert!(!updated.cooling.thermal().links().contains_key(&CONTACT));
    assert_eq!(
        updated.cooling.thermal().cells(),
        heated.cooling.thermal().cells()
    );
    assert_eq!(
        updated.mechanical.motion.state().displacements_m(),
        heated.mechanical.motion.state().displacements_m()
    );
    assert_eq!(
        updated.mechanical.motion.velocity_m_s(),
        heated.mechanical.motion.velocity_m_s()
    );
    assert_eq!(updated.mechanical.motion.time_s(), 0.01);
    assert_eq!(
        updated.mechanical.motion.steps(),
        heated.mechanical.motion.steps()
    );
    assert!(core::ptr::eq(
        updated.mechanical.motion.mass_operator(),
        motion.mass_operator()
    ));

    // The SAME original mass now accelerates under the reduced restoring force.
    // Midpoint inertia + the unchanged bulk spring has this independent root.
    let continued = propose_deformable_cohesive_cooling(
        &updated.mechanical.motion,
        &updated.cooling,
        &[force],
        &[[0.0; 3]],
        0.02,
        &partition(1.0, 1.0),
        cooling_controls(),
        &gate,
    )
    .unwrap();
    let dt = 0.02 - updated.mechanical.motion.time_s();
    let c = 2.0 * mass.0 / (dt * dt);
    let expected = ((c - 0.5 * elastic.0) * gap
        + 2.0 * mass.0 * updated.mechanical.motion.velocity_m_s()[0] / dt
        + force)
        / (c + 0.5 * elastic.0);
    close(
        continued.mechanical.motion.state().displacements_m()[0],
        expected,
        1e-12,
    );
    assert!(expected > gap + 1e-5);
    assert_eq!(continued.heat.contact.fracture_dissipation_j, 0.0);
    assert_eq!(continued.heat.contact.viscous_dissipation_j, 0.0);
    for cell in continued.cooling.thermal().cells().values() {
        assert_eq!(cell.specific_h_j_kg, 100.0);
    }

    let cooled = propose_deformable_cohesive_cooling(
        &continued.mechanical.motion,
        &continued.cooling,
        &[force],
        &[[0.0; 3]],
        0.03,
        &partition(1.0, 1.0),
        cooling_controls(),
        &gate,
    )
    .unwrap();
    let sources_cold = prepare_deformable_cohesive_temperature_sources(
        &cooled.mechanical.motion,
        &cooled.cooling,
        &inputs,
        "T",
        &gate,
    )
    .unwrap();
    close(
        sources_cold[0].resolved().state_coordinate("T").unwrap(),
        300.0,
        1e-12,
    );
    let target_cold = target
        .with_temperature_sources_controlled(&[&sources_cold[0]], "T", &mut || false)
        .unwrap();
    let returned = propose_deformable_cohesive_temperature(
        &cooled.mechanical.motion,
        &cooled.cooling,
        &target_cold,
        &inputs,
        "T",
        controls(),
        &gate,
    )
    .unwrap();
    assert!(returned.mechanical.motion.state().points()[0].failed());
    close(
        returned.mechanical.motion.state().points()[0].dissipated_j_m2(),
        changed_storage,
        1e-12,
    );
    assert_eq!(returned.mechanical.energy.parameter_work_j, 0.0);
    assert_eq!(returned.mechanical.energy.damage_dissipated_j, 0.0);
    assert!(!returned.cooling.interfaces()[&CONTACT].closed);
    assert_eq!(original.thermal().time_s(), 0.0);
    assert!(!motion.state().points()[0].failed());
}

#[test]
fn weighted_mixed_viscous_feedback_retains_onset_direction_and_exact_no_op() {
    let (cold, card) = thermal_material(true, true);
    let elastic = Scalar(1000.0);
    let mass = Scalar(1.0);
    let joint = system(&elastic, &cold, &[CONTACT], false);
    let state = joint
        .bonded_state()
        .solve_over_time_controlled(
            &[0.0],
            &[[0.0015, 0.0005, -0.0002]],
            1.0,
            dynamics().solve,
            &mut || false,
        )
        .unwrap()
        .state;
    let motion = state
        .into_dynamics(&mass, vec![0.02], 0.0, dynamics())
        .unwrap();
    let original = cooling(
        network(1, [1.0; 2], [1.0; 2], [1.0; 2], [0.0; 2], [100.0, 0.0]),
        motion.state().points()[0],
    );
    let gate = CancelGate::new_clock_free();
    let inputs = binding(&card, 0.25);
    let sources =
        prepare_deformable_cohesive_temperature_sources(&motion, &original, &inputs, "T", &gate)
            .unwrap();
    assert_eq!(sources[0].resolved().state_coordinate("T"), Some(325.0));
    let target = joint
        .with_temperature_sources_controlled(&[&sources[0]], "T", &mut || false)
        .unwrap();
    let next = propose_deformable_cohesive_temperature(
        &motion,
        &original,
        &target,
        &inputs,
        "T",
        controls(),
        &gate,
    )
    .unwrap();
    let before = motion.state().points()[0];
    let after = next.mechanical.motion.state().points()[0];
    assert!(after.damage() >= before.damage());
    assert_eq!(
        after.mixed_mode_history().unwrap().onset_direction,
        before.mixed_mode_history().unwrap().onset_direction
    );
    assert_eq!(after.jump_m(), before.jump_m());
    assert_eq!(next.mechanical.motion.velocity_m_s(), motion.velocity_m_s());
    assert_eq!(next.contact.contact.viscous_dissipation_j, 0.0);
    assert_eq!(next.cooling.thermal().cells(), original.thermal().cells());
    for (old, new) in cold
        .resolved()
        .properties()
        .iter()
        .zip(sources[0].resolved().properties())
    {
        assert_eq!(old.answer().receipt.selected, new.answer().receipt.selected);
    }
    let again = prepare_deformable_cohesive_temperature_sources(
        &next.mechanical.motion,
        &next.cooling,
        &inputs,
        "T",
        &gate,
    )
    .unwrap();
    assert_eq!(
        again[0], sources[0],
        "unchanged actual T retains the exact source clone"
    );
    let same = propose_deformable_cohesive_temperature(
        &next.mechanical.motion,
        &next.cooling,
        &target,
        &inputs,
        "T",
        controls(),
        &gate,
    )
    .unwrap();
    assert_eq!(
        same.mechanical.motion.state().points(),
        next.mechanical.motion.state().points()
    );
    assert_eq!(same.mechanical.energy.parameter_work_j, 0.0);
    assert_eq!(same.mechanical.energy.damage_dissipated_j, 0.0);
    assert_eq!(same.mechanical.energy.stored_change_j, 0.0);
    assert_eq!(same.combined_energy_residual_j, 0.0);
}

#[test]
fn foreign_source_stale_temperature_clock_and_geometry_refuse_both_owners() {
    let (cold, card) = thermal_material(false, false);
    let elastic = Scalar(1000.0);
    let mass = Scalar(1.0);
    let joint = system(&elastic, &cold, &[CONTACT], false);
    let state = joint
        .bonded_state()
        .solve_controlled(&[0.0], &[[0.0005, 0.0, 0.0]], dynamics().solve, &mut || {
            false
        })
        .unwrap()
        .state;
    let motion = state
        .into_dynamics(&mass, vec![0.02], 0.0, dynamics())
        .unwrap();
    let original = cooling(
        network(1, [1.0; 2], [1.0; 2], [1.0; 2], [0.0; 2], [100.0; 2]),
        motion.state().points()[0],
    );
    let gate = CancelGate::new_clock_free();
    let inputs = binding(&card, 0.5);
    let (_, foreign_card) = material_and_card(false, None);
    assert!(
        prepare_deformable_cohesive_temperature_sources(
            &motion,
            &original,
            &binding(&foreign_card, 0.5),
            "T",
            &gate
        )
        .is_err()
    );
    let mut wrong_id = inputs;
    wrong_id[0].patch_id += 1;
    assert!(
        prepare_deformable_cohesive_temperature_sources(&motion, &original, &wrong_id, "T", &gate)
            .is_err()
    );
    let (stale, _) = original
        .advance_to(0.001, cooling_controls().thermal, || Ok(()))
        .unwrap();
    assert!(
        prepare_deformable_cohesive_temperature_sources(&motion, &stale, &inputs, "T", &gate)
            .is_err()
    );
    let expired = cooling(
        network(1, [1.0; 2], [1.0; 2], [1.0; 2], [0.0; 2], [150.0; 2]),
        motion.state().points()[0],
    );
    assert!(matches!(
        prepare_deformable_cohesive_temperature_sources(&motion, &expired, &inputs, "T", &gate),
        Err(CohesiveCoolingError::Source { .. })
    ));

    let wrong_temperature = cold.at_temperature(&card, "T", 350.0).unwrap();
    let wrong = joint
        .with_temperature_sources_controlled(&[&wrong_temperature], "T", &mut || false)
        .unwrap();
    assert!(matches!(
        propose_deformable_cohesive_temperature(
            &motion,
            &original,
            &wrong,
            &inputs,
            "T",
            controls(),
            &gate
        ),
        Err(CohesiveCoolingError::Temperature(_))
    ));
    let hot = cold.at_temperature(&card, "T", 400.0).unwrap();
    let target = joint
        .with_temperature_sources_controlled(&[&hot], "T", &mut || false)
        .unwrap();
    let foreign_operator = Scalar(1000.0);
    let replaced = system(&foreign_operator, &hot, &[CONTACT], false);
    assert!(matches!(
        propose_deformable_cohesive_temperature(
            &motion,
            &original,
            &replaced,
            &inputs,
            "T",
            controls(),
            &gate
        ),
        Err(CohesiveCoolingError::Mechanical(_))
    ));
    gate.request();
    assert!(
        propose_deformable_cohesive_temperature(
            &motion,
            &original,
            &target,
            &inputs,
            "T",
            controls(),
            &gate
        )
        .is_err()
    );
    let retry = propose_deformable_cohesive_temperature(
        &motion,
        &original,
        &target,
        &inputs,
        "T",
        controls(),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    assert!(retry.mechanical.motion.state().points()[0].failed());
    assert!(!motion.state().points()[0].failed());
    assert_eq!(
        original.interfaces()[&CONTACT].point,
        motion.state().points()[0]
    );
    assert_eq!(original.thermal().time_s(), 0.0);
    assert_eq!(motion.time_s(), 0.0);
    assert_eq!(motion.velocity_m_s(), &[0.02]);
}

struct CancelOnTemperature {
    gate: &'static CancelGate,
    armed: Cell<bool>,
}
impl ThermalLaw for CancelOnTemperature {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if self.armed.get() {
            self.gate.request();
        }
        Ok(300.0 + h)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
}

#[test]
fn cancellation_after_mechanical_transfer_keeps_the_original_pair_retryable() {
    let (cold, card) = thermal_material(false, false);
    let elastic = Scalar(1000.0);
    let mass = Scalar(1.0);
    let joint = system(&elastic, &cold, &[CONTACT], false);
    let state = joint
        .bonded_state()
        .solve_controlled(&[0.0], &[[0.0005, 0.0, 0.0]], dynamics().solve, &mut || {
            false
        })
        .unwrap()
        .state;
    let motion = state
        .into_dynamics(&mass, vec![0.02], 0.0, dynamics())
        .unwrap();
    let gate = Box::leak(Box::new(CancelGate::new_clock_free()));
    let law = Box::leak(Box::new(CancelOnTemperature {
        gate,
        armed: Cell::new(false),
    }));
    let original = cooling(
        network_with_law(1, [1.0; 2], [1.0; 2], [1.0; 2], [0.0; 2], [100.0; 2], law),
        motion.state().points()[0],
    );
    let inputs = binding(&card, 0.5);
    let sources =
        prepare_deformable_cohesive_temperature_sources(&motion, &original, &inputs, "T", gate)
            .unwrap();
    let target = joint
        .with_temperature_sources_controlled(&[&sources[0]], "T", &mut || false)
        .unwrap();
    law.armed.set(true);
    assert!(
        propose_deformable_cohesive_temperature(
            &motion,
            &original,
            &target,
            &inputs,
            "T",
            controls(),
            gate
        )
        .is_err()
    );
    assert!(gate.is_requested());
    assert_eq!(
        motion.state().points()[0],
        original.interfaces()[&CONTACT].point
    );
    assert!(!motion.state().points()[0].failed());
    law.armed.set(false);
    let retry = propose_deformable_cohesive_temperature(
        &motion,
        &original,
        &target,
        &inputs,
        "T",
        controls(),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    close(
        retry.contact.contact.fracture_dissipation_j,
        0.5 * 400.0 * 0.0005_f64.powi(2),
        1e-12,
    );
    assert_eq!(retry.cooling.thermal().cells(), original.thermal().cells());
    assert_eq!(retry.mechanical.motion.time_s(), motion.time_s());
}
