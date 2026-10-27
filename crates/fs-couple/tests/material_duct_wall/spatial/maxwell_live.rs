//! Original source claims refresh a continuing native hereditary wall.
use super::*;
use fs_qty::Temperature;

fn source() -> (MaterialCard, ResolvedMaterialStatePoint) {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, dims, values) in [
        ("density", Density::DIMS, [1000.0; 3]),
        ("young_modulus", Pressure::DIMS, [1e6; 3]),
        ("poisson_ratio", Dims([0; 6]), [0.3; 3]),
        ("arm_a_modulus", Pressure::DIMS, [8e5, 1.6e6, 0.0]),
        (
            "arm_a_time",
            Dims([0, 0, 1, 0, 0, 0]),
            [0.0005, 0.0015, 0.0015],
        ),
        ("arm_b_modulus", Pressure::DIMS, [3e5; 3]),
        (
            "arm_b_time",
            Dims([0, 0, 1, 0, 0, 0]),
            [0.0002, 0.0008, 0.0008],
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(300.0, values[0]), (360.0, values[1]), (420.0, values[2])],
                    dims,
                },
                validity: ValidityDomain::unconstrained().with("T", 300.0, 420.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source:
                        "synthetic temperature-dependent Maxwell arms; not measured qualification"
                            .into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_with_quantity(
                name,
                QuantitySpec::dimensional(dims),
                ScalarAdmissibility::Finite,
            )
            .unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic hereditary membrane".into(),
            phase: "solid".into(),
            process: "fixed geometry test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let resolved = resolve_material_state_point(
        &card,
        &QueryPoint::new().with("T", 300.0).unwrap(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    (card, resolved)
}

fn wall(material: &ResolvedMaterialStatePoint, reversed: bool) -> LocalMaterialWall {
    let mut pairs = vec![
        ("arm_a_modulus".into(), "arm_a_time".into()),
        ("arm_b_modulus".into(), "arm_b_time".into()),
    ];
    if reversed {
        pairs.reverse();
    }
    let relaxation = WallRelaxation::from_state_point(material, &pairs, BAND, 1.0, 0.8).unwrap();
    LocalMaterialWall::new(
        material.clone(),
        0.0002,
        AxialConstraint::StressFree,
        WallMaterialLoss::Maxwell(relaxation),
        BAND,
    )
    .unwrap()
}

#[test]
fn source_temperature_and_reordered_arms_change_the_continuing_gas_wall_trajectory() {
    let (card, cold_source) = source();
    let hot_source = cold_source
        .at_query_point(&card, &cold_source.query().with("T", 360.0).unwrap())
        .unwrap();
    let cold = wall(&cold_source, false);
    let hot = wall(&hot_source, true);
    let field = gas_field(0.2, 350.0).unwrap();
    let mut bound = RegionalMaterialWallDg::new(vec![cold.clone(); 2], duct(&field)).unwrap();
    let mut propagation = bound.duct().propagator(2e-6, 1e-13).unwrap();
    let mut state = bound.duct().project_initial_state(|_| 10.0, |_| 0.0);
    for _ in 0..30 {
        propagation.step(&mut state, &[]).unwrap();
    }
    propagation.sync_nodal(&mut state).unwrap();
    let before = state.clone();
    let original_motion = bound.duct().compliant_wall_state(&state).unwrap();
    assert_eq!(original_motion.relaxation_pressure_pa.len(), 2);
    assert!(original_motion.relaxation_stored_energy_j > 0.0);
    let mut control = state.clone();
    let control_propagation = propagation.clone();
    let account = bound
        .rebind_controlled(
            vec![cold.clone(), hot.clone()],
            &mut state,
            &mut propagation,
            1.0,
            &mut || false,
        )
        .unwrap();
    assert_eq!(bound.walls(), &[cold, hot]);
    assert_eq!(bound.duct().gas_field(), Some(&field));
    assert_eq!(state.step_count(), before.step_count());
    assert_eq!(state.initial_energy_j(), before.initial_energy_j());
    assert_eq!(
        state.cumulative_dissipated_j(),
        before.cumulative_dissipated_j()
    );
    assert_eq!(state.pressure(), before.pressure());
    assert_eq!(state.volume_flow(), before.volume_flow());
    assert!(account.parameter_work_j > 0.0);
    close(
        account.energy_after_j - account.energy_before_j,
        account.parameter_work_j,
    );
    let after = bound.duct().compliant_wall_state(&state).unwrap();
    assert_eq!(after.relaxation_pressure_pa.len(), 4);
    let local = bound.duct().mesh().order + 1;
    for node in 0..bound.duct().field_dofs() {
        close(
            after.displacement_m[node],
            original_motion.displacement_m[node],
        );
        close(after.velocity_m_s[node], original_motion.velocity_m_s[node]);
        let right = node >= local;
        // A doubles its spring modulus only on the hot right element. Its
        // elastic extension and viscous offset survive; B only changes tau.
        close(
            after.relaxation_pressure_pa[if right { 3 } else { 0 }][node],
            (if right { 2.0 } else { 1.0 }) * original_motion.relaxation_pressure_pa[0][node],
        );
        close(
            after.relaxation_pressure_pa[if right { 2 } else { 1 }][node],
            original_motion.relaxation_pressure_pa[1][node],
        );
    }
    let mut signal = 0.0;
    let mut difference = 0.0;
    for _ in 0..40 {
        let frame = propagation.step(&mut state, &[]).unwrap();
        control_propagation.step(&mut control, &[]).unwrap();
        propagation.sync_nodal(&mut state).unwrap();
        control_propagation.sync_nodal(&mut control).unwrap();
        for (&a, &b) in state.pressure().iter().zip(control.pressure()) {
            signal += b * b;
            difference += (a - b) * (a - b);
        }
        assert_eq!(
            frame.ledger.cumulative_parameter_work_j,
            account.parameter_work_j
        );
        close(
            state.stored_energy_j() + state.cumulative_dissipated_j(),
            before.initial_energy_j() + account.parameter_work_j,
        );
    }
    assert!(difference > 1e-10 * signal);
    assert_eq!(state.step_count(), 70);
    assert_eq!(propagation.dt(), control_propagation.dt());
}

#[test]
fn mismatched_source_history_geometry_work_and_cancellation_are_atomic() {
    let (card, cold_source) = source();
    let hot_source = cold_source
        .at_query_point(&card, &cold_source.query().with("T", 360.0).unwrap())
        .unwrap();
    let vanished_source = cold_source
        .at_query_point(&card, &cold_source.query().with("T", 420.0).unwrap())
        .unwrap();
    let cold = wall(&cold_source, false);
    let hot = wall(&hot_source, true);
    let field = gas_field(0.2, 350.0).unwrap();
    let original = cold.clone().bind_dg(duct(&field)).unwrap();
    let original_propagation = original.duct().propagator(2e-6, 1e-13).unwrap();
    let mut initial = original.duct().project_initial_state(|_| 10.0, |_| 0.0);
    for _ in 0..20 {
        original_propagation.step(&mut initial, &[]).unwrap();
    }
    original_propagation.sync_nodal(&mut initial).unwrap();
    assert!(
        original
            .duct()
            .compliant_wall_state(&initial)
            .unwrap()
            .relaxation_stored_energy_j
            > 0.0
    );
    let WallMaterialLoss::Maxwell(cold_law) = cold.loss() else {
        unreachable!()
    };
    let authored = WallRelaxation::new(
        cold_law.material().clone(),
        "same coefficients, no identified source arms".into(),
        BAND,
        1.0,
        0.8,
    )
    .unwrap();
    let relinked = WallRelaxation::from_state_point(
        &cold_source,
        &[
            ("arm_a_modulus".into(), "arm_b_time".into()),
            ("arm_b_modulus".into(), "arm_a_time".into()),
        ],
        BAND,
        1.0,
        0.8,
    )
    .unwrap();
    let with_loss = |loss| {
        LocalMaterialWall::new(
            cold_source.clone(),
            0.0002,
            AxialConstraint::StressFree,
            loss,
            BAND,
        )
        .unwrap()
    };
    let thicker = LocalMaterialWall::new(
        cold_source.clone(),
        0.00025,
        AxialConstraint::StressFree,
        cold.loss().clone(),
        BAND,
    )
    .unwrap();
    for (replacement, budget) in [
        (hot.clone(), 0.0),
        (wall(&vanished_source, false), 1.0),
        (with_loss(WallMaterialLoss::Maxwell(authored)), 1.0),
        (with_loss(WallMaterialLoss::Maxwell(relinked)), 1.0),
        (source_wall(1000.0, 1e6), 1.0),
        (thicker, 1.0),
    ] {
        let mut bound = original.clone();
        let mut state = initial.clone();
        let mut propagation = original_propagation.clone();
        assert!(
            bound
                .rebind_controlled(
                    replacement,
                    &mut state,
                    &mut propagation,
                    budget,
                    &mut || false
                )
                .is_err()
        );
        assert_eq!(state, initial);
        assert_eq!(bound.wall(), original.wall());
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
    // Refuse after the final native propagator has been prepared, so even this
    // latest cancellation cannot publish a source with stale/no memory state.
    let mut polls = 0;
    original
        .clone()
        .rebind_controlled(
            hot.clone(),
            &mut initial.clone(),
            &mut original_propagation.clone(),
            1.0,
            &mut || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut count = 0;
    let mut bound = original.clone();
    let mut state = initial.clone();
    let mut propagation = original_propagation.clone();
    assert!(matches!(
        bound.rebind_controlled(hot, &mut state, &mut propagation, 1.0, &mut || {
            count += 1;
            count == polls
        }),
        Err(AcousticRealizeError::Cancelled)
    ));
    assert_eq!(state, initial);
    assert_eq!(bound.wall(), original.wall());
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
