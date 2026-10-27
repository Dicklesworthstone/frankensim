//! Live source-temperature changes, not freshly initialized comparison voices.

use super::*;
use fs_qty::Temperature;

fn thermal_card() -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, values) in [
        ("density", Density::DIMS, [1000.0, 1000.0]),
        ("young_modulus", Pressure::DIMS, [1e6, 5e5]),
        ("poisson_ratio", Dims([0; 6]), [0.3, 0.3]),
        (YOUNG_VISCOSITY_PROPERTY, DynViscosity::DIMS, [2.0, 5.0]),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(300.0, values[0]), (360.0, values[1])],
                    dims,
                },
                validity: ValidityDomain::unconstrained().with("T", 300.0, 360.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic temperature-dependent wall; not measured qualification"
                        .into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic thermal membrane".into(),
            phase: "solid".into(),
            process: "fixed geometry test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn wall_at(
    card: &MaterialCard,
    temperature: f64,
) -> Result<LocalMaterialWall, AcousticRealizeError> {
    let requirements: Vec<_> = properties(1000.0, 1e6, 2.0)
        .into_iter()
        .map(|(name, dims, _)| {
            ScalarPropertyRequirement::try_with_quantity(
                name,
                QuantitySpec::dimensional(dims),
                ScalarAdmissibility::Finite,
            )
            .unwrap()
        })
        .collect();
    let material = resolve_material_state_point(
        card,
        &QueryPoint::new().with("T", temperature).unwrap(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
    LocalMaterialWall::new(
        material,
        0.0002,
        AxialConstraint::StressFree,
        WallMaterialLoss::KelvinVoigt,
        BAND,
    )
}

#[test]
fn same_source_at_new_temperature_changes_live_wall_and_pressure_without_restart() {
    let card = thermal_card();
    let cold = wall_at(&card, 300.0).unwrap();
    let hot = wall_at(&card, 360.0).unwrap();
    assert_eq!(
        cold.material().card_identity(),
        hot.material().card_identity()
    );
    assert!(wall_at(&card, 361.0).is_err());
    let field = gas_field(0.2, 350.0).unwrap();
    let mut bound = cold.clone().bind_dg(duct(&field)).unwrap();
    let mut propagation = bound.duct().propagator(2e-6, 1e-13).unwrap();
    let mut state = bound.duct().project_initial_state(|_| 10.0, |_| 0.0);
    for _ in 0..30 {
        propagation.step(&mut state, &[]).unwrap();
    }
    propagation.sync_nodal(&mut state).unwrap();
    let before = state.clone();
    let original_motion = bound.duct().compliant_wall_state(&state).unwrap();
    let mut control = state.clone();
    let control_propagation = propagation.clone();
    let original_initial_energy = state.initial_energy_j();
    let update = bound
        .rebind_controlled(hot.clone(), &mut state, &mut propagation, 1.0, &mut || {
            false
        })
        .unwrap();
    assert_eq!(bound.wall(), &hot);
    assert_eq!(bound.duct().gas_field(), Some(&field));
    assert_eq!(state.step_count(), 30);
    assert_eq!(state.initial_energy_j(), original_initial_energy);
    assert_eq!(state.pressure(), before.pressure());
    assert_eq!(state.volume_flow(), before.volume_flow());
    assert_eq!(
        state.cumulative_dissipated_j(),
        before.cumulative_dissipated_j()
    );
    let transferred = bound.duct().compliant_wall_state(&state).unwrap();
    for i in 0..bound.duct().field_dofs() {
        close(
            transferred.displacement_m[i],
            original_motion.displacement_m[i],
        );
        close(transferred.velocity_m_s[i], original_motion.velocity_m_s[i]);
    }
    assert!(update.parameter_work_j < 0.0);
    close(
        update.energy_after_j - update.energy_before_j,
        update.parameter_work_j,
    );
    assert_eq!(state.cumulative_parameter_work_j(), update.parameter_work_j);
    let mut signal = 0.0;
    let mut difference = 0.0;
    for _ in 0..30 {
        let frame = propagation.step(&mut state, &[]).unwrap();
        control_propagation.step(&mut control, &[]).unwrap();
        propagation.sync_nodal(&mut state).unwrap();
        control_propagation.sync_nodal(&mut control).unwrap();
        signal += control.pressure()[0].powi(2);
        difference += (state.pressure()[0] - control.pressure()[0]).powi(2);
        assert_eq!(
            frame.ledger.cumulative_parameter_work_j,
            update.parameter_work_j
        );
        assert!(frame.ledger.energy_balance_residual_j.abs() < 1e-10 * original_initial_energy);
        close(
            state.stored_energy_j() + state.cumulative_dissipated_j(),
            original_initial_energy + update.parameter_work_j,
        );
    }
    assert!(difference > 1e-8 * signal);
    assert_eq!(state.step_count(), 60);
    assert_eq!(propagation.dt(), control_propagation.dt());
    assert_eq!(propagation.tolerance(), control_propagation.tolerance());
}

#[test]
fn regional_live_sources_work_budget_and_late_cancellation_are_atomic() {
    let field = gas_field(0.2, 350.0).unwrap();
    let left = kelvin_wall(1000.0, 1e6, 2.0, 0.0002);
    let right = kelvin_wall(1800.0, 4e6, 5.0, 0.00015);
    let original =
        RegionalMaterialWallDg::new(vec![left.clone(), right.clone()], duct(&field)).unwrap();
    let original_propagation = original.duct().propagator(2e-6, 1e-13).unwrap();
    let mut initial = original.duct().project_initial_state(|_| 10.0, |_| 0.0);
    for _ in 0..20 {
        original_propagation.step(&mut initial, &[]).unwrap();
    }
    original_propagation.sync_nodal(&mut initial).unwrap();
    let stiffer = kelvin_wall(1000.0, 2e6, 2.0, 0.0002);
    let mut bound = original.clone();
    let mut state = initial.clone();
    let mut propagation = original_propagation.clone();
    let replacement = vec![stiffer.clone(), right.clone()];
    let update = bound
        .rebind_controlled(
            replacement.clone(),
            &mut state,
            &mut propagation,
            1.0,
            &mut || false,
        )
        .unwrap();
    assert_eq!(bound.walls(), replacement.as_slice());
    assert_eq!(bound.band_hz(), original.band_hz());
    assert!(update.parameter_work_j > 0.0);
    assert_eq!(state.step_count(), initial.step_count());
    assert_eq!(state.pressure(), initial.pressure());
    let native_motion = bound.duct().compliant_wall_state(&state).unwrap();
    let prior_motion = original.duct().compliant_wall_state(&initial).unwrap();
    for (a, b) in native_motion
        .displacement_m
        .iter()
        .zip(&prior_motion.displacement_m)
    {
        close(*a, *b);
    }

    let narrower = LocalMaterialWall::new(
        left.material().clone(),
        left.thickness_m(),
        left.axial_constraint(),
        WallMaterialLoss::KelvinVoigt,
        [100.0, 1000.0],
    )
    .unwrap();
    for (walls, budget) in [
        (replacement.clone(), 0.0),
        (vec![narrower, right.clone()], 1.0),
        (vec![source_wall(1000.0, 1e6), right.clone()], 1.0),
    ] {
        let mut bound = original.clone();
        let mut state = initial.clone();
        let mut propagation = original_propagation.clone();
        assert!(
            bound
                .rebind_controlled(walls, &mut state, &mut propagation, budget, &mut || false)
                .is_err()
        );
        assert_eq!(state, initial);
        assert_eq!(bound.walls(), original.walls());
        assert_eq!(
            bound.duct().operator_identity(),
            original.duct().operator_identity()
        );
        let mut control = initial.clone();
        assert_eq!(
            propagation.step(&mut state, &[]).unwrap(),
            original_propagation.step(&mut control, &[]).unwrap()
        );
    }
    let mut polls = 0;
    let mut probe = original.clone();
    probe
        .rebind_controlled(
            replacement.clone(),
            &mut initial.clone(),
            &mut original_propagation.clone(),
            1.0,
            &mut || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut calls = 0;
    let mut bound = original.clone();
    let mut state = initial.clone();
    let mut propagation = original_propagation.clone();
    assert!(matches!(
        bound.rebind_controlled(
            replacement.clone(),
            &mut state,
            &mut propagation,
            1.0,
            &mut || {
                calls += 1;
                calls == polls
            }
        ),
        Err(AcousticRealizeError::Cancelled)
    ));
    assert_eq!(state, initial);
    assert_eq!(bound.walls(), original.walls());
    assert_eq!(
        bound.duct().operator_identity(),
        original.duct().operator_identity()
    );
    let mut foreign = duct(&field).propagator(2e-6, 1e-13).unwrap();
    assert!(
        bound
            .rebind_controlled(replacement, &mut state, &mut foreign, 1.0, &mut || false)
            .is_err()
    );
    assert_eq!(state, initial);
    assert_eq!(bound.walls(), original.walls());

    // Upstream regional binding supports partial Maxwell histories. A live
    // elastic replacement must still refuse instead of truncating that arm.
    let mut mixed =
        RegionalMaterialWallDg::new(vec![source_wall(1000.0, 1e6), right.clone()], duct(&field))
            .unwrap();
    let mut mixed_propagation = mixed.duct().propagator(2e-6, 1e-13).unwrap();
    let mut memory = mixed.duct().project_initial_state(|_| 10.0, |_| 0.0);
    mixed_propagation.step(&mut memory, &[]).unwrap();
    assert!(
        mixed
            .duct()
            .compliant_wall_state(&memory)
            .unwrap()
            .relaxation_stored_energy_j
            > 0.0
    );
    let saved_memory = memory.clone();
    let saved_sources = mixed.walls().to_vec();
    let saved_identity = mixed.duct().operator_identity();
    assert!(
        mixed
            .rebind_controlled(
                vec![left, right],
                &mut memory,
                &mut mixed_propagation,
                1.0,
                &mut || false
            )
            .is_err()
    );
    assert_eq!(memory, saved_memory);
    assert_eq!(mixed.walls(), saved_sources.as_slice());
    assert_eq!(mixed.duct().operator_identity(), saved_identity);
}
