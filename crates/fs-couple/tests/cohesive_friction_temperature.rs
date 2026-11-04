//! Paired-source temperature transfer of the original cohesive friction owners.
use fs_couple::cohesive_contact::frictional_damage::{
    DamageFrictionConfig, DamageFrictionError, DamageFrictionLoadLine, DamageFrictionRegime,
};
use fs_couple::cohesive_contact::reclosure::{
    ReclosureActivation, ReclosureConfig, ReclosureError, ReclosureLoadLine, ReclosureRegime,
    ResidualCoulombSource,
    source::{
        RESIDUAL_KINETIC_FRICTION_PROPERTY, RESIDUAL_STATIC_FRICTION_PROPERTY,
        resolve_residual_coulomb,
    },
};
use fs_couple::cohesive_contact::CohesiveContactConfig;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SurfaceSpec, SystemContext,
    UncertaintyModel,
};
use fs_material::cohesive::CohesivePoint;
use fs_material::cohesive::source::*;
use fs_material::state_point::{MaterialPropertySelection, interface_stiffness_per_area_dims};
use fs_qty::{Dims, Pressure, Temperature};

const AREA: f64 = 1e-6;
const JUMP: [f64; 3] = [-2e-4, 0.0, 0.0];

fn card(hot_mu: [f64; 2]) -> InterfaceSystemCard {
    let stiffness = interface_stiffness_per_area_dims();
    let mut claims = ClaimSet::new();
    for (name, dims, cold, hot) in [
        (COHESIVE_STIFFNESS_PROPERTY, stiffness, 1e9, 5e8),
        (COHESIVE_STRENGTH_PROPERTY, Pressure::DIMS, 2e5, 8e4),
        (
            COHESIVE_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
            100.0,
            50.0,
        ),
        (COHESIVE_COMPRESSION_STIFFNESS_PROPERTY, stiffness, 2e9, 3e9),
        (RESIDUAL_STATIC_FRICTION_PROPERTY, Dims::NONE, 0.6, hot_mu[0]),
        (RESIDUAL_KINETIC_FRICTION_PROPERTY, Dims::NONE, 0.3, hot_mu[1]),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(280.0, cold), (320.0, hot)],
                    dims,
                },
                validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured temperature-friction fixture, not measured data".into(),
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
            process: "original temperature-friction fixture".into(),
            revision: 1,
        },
        texture_frame: name.into(),
    };
    InterfaceSystemCard::assemble(
        surface("A"),
        surface("B"),
        SystemContext {
            medium: "explicit interface".into(),
            third_body: None,
            environment: "temperature query".into(),
            history: "original bond and residual friction".into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn sources(card: &InterfaceSystemCard) -> (ResolvedBilinearCohesive, ResidualCoulombSource) {
    let query = QueryPoint::new().with("T", 280.0).unwrap();
    (
        resolve_bilinear_cohesive(
            card,
            &query,
            MaterialPropertySelection::SingleClaimOnly,
            0.02,
        )
        .unwrap(),
        resolve_residual_coulomb(card, &query, MaterialPropertySelection::SingleClaimOnly).unwrap(),
    )
}

fn receiving(source: &ResolvedBilinearCohesive) -> [&MaterialStateId; 2] {
    [&source.surfaces()[0].material, &source.surfaces()[1].material]
}

fn compliance() -> [[f64; 3]; 3] {
    [[2e-5, 2e-6, -1e-6], [2e-6, 3e-5, 4e-6], [-1e-6, 4e-6, 4e-5]]
}

fn free_for(jump: [f64; 3], force: [f64; 3]) -> [f64; 3] {
    core::array::from_fn(|i| {
        jump[i] + (0..3).map(|j| compliance()[i][j] * force[j]).sum::<f64>()
    })
}

fn reclosure_config() -> ReclosureConfig {
    ReclosureConfig {
        max_iterations: 128,
        max_evaluations: 4096,
        jump_absolute_tolerance_m: 1e-14,
        jump_relative_tolerance: 1e-10,
        maximum_jump_increment_m: 0.01,
        maximum_force_n: 100.0,
        force_absolute_tolerance_n: 1e-11,
        force_relative_tolerance: 1e-9,
        material_energy_tolerance_j_m2: 1e-8,
        work_absolute_tolerance_j: 1e-13,
        work_relative_tolerance: 1e-8,
        minimum_monotonicity_margin: 0.02,
    }
}

fn damage_config() -> DamageFrictionConfig {
    DamageFrictionConfig {
        contact: CohesiveContactConfig {
            max_iterations: 48,
            max_material_evaluations: 8192,
            max_backtracks: 16,
            jump_absolute_tolerance_m: 1e-14,
            jump_relative_tolerance: 1e-10,
            maximum_jump_increment_m: 0.01,
            maximum_force_n: 100.0,
            material_energy_tolerance_j_m2: 1e-8,
        },
        max_friction_iterations: 96,
        force_absolute_tolerance_n: 1e-11,
        force_relative_tolerance: 1e-9,
        work_absolute_tolerance_j: 1e-13,
        work_relative_tolerance: 1e-8,
    }
}

fn history(source: &ResolvedBilinearCohesive, maximum: f64) -> CohesivePoint<'_> {
    source
        .law()
        .bonded_point()
        .trial([maximum, 0.0, 0.0], 1e-8)
        .unwrap()
        .point
        .trial(JUMP, 1e-8)
        .unwrap()
        .point
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}, tolerance={tolerance:e}");
}

#[test]
fn reclosure_temperature_keeps_failed_history_and_occupied_kinetic_branch() {
    let card = card([0.3, 0.1]);
    let (cold, friction) = sources(&card);
    let hot = cold.at_temperature(&card, "T", 320.0).unwrap();
    let hot_friction = friction.at_temperature(&card, "T", 320.0).unwrap();
    let owner =
        ReclosureLoadLine::new(&cold, &friction, AREA, compliance(), reclosure_config()).unwrap();
    let point = history(&cold, 0.002);
    let original = owner
        .adopt(point, receiving(&cold), ReclosureActivation::Sliding)
        .unwrap();
    let (heated, retained, transfer) = owner
        .transfer_temperature(original, &hot, &hot_friction, "T")
        .unwrap();
    assert_eq!(retained.point(), transfer.point);
    assert_eq!(retained.point().jump_m(), JUMP);
    assert_eq!(retained.point().damage(), 1.0);
    assert_eq!(retained.point().dissipated_j_m2(), point.dissipated_j_m2());
    assert_eq!(retained.occupancy(), ReclosureRegime::Sliding);
    assert!(core::ptr::eq(retained.point().law(), hot.law()));
    assert_eq!(heated.config(), owner.config());
    assert_eq!(heated.area_m2(), owner.area_m2());
    assert_eq!(heated.compliance_m_per_n(), owner.compliance_m_per_n());
    close(transfer.parameter_work_j_m2, 20.0, 1e-12);
    close(transfer.stored_change_j_m2, 20.0, 1e-12);
    assert_eq!(transfer.damage_dissipated_j_m2, 0.0);

    // This load fits the hot static disk but exceeds its kinetic disk. A
    // mistaken readoption at rest would suppress the actual ongoing slip.
    let end = [JUMP[0], 1e-6, 0.0];
    let expected = [-0.6, 0.06, 0.0];
    let free = free_for(end, expected);
    let sliding = heated.trial(retained, free, 0.01, || true).unwrap();
    assert_eq!(sliding.audit.regime, ReclosureRegime::Sliding);
    for i in 0..3 {
        close(sliding.mean_force_n[i], expected[i], 2e-9);
        close(sliding.audit.state.point().jump_m()[i], end[i], 2e-12);
    }
    close(sliding.audit.friction_dissipation_j, 0.06e-6, 2e-13);
    let at_rest = heated
        .adopt(retained.point(), receiving(&hot), ReclosureActivation::AtRest)
        .unwrap();
    let held = heated.trial(at_rest, free, 0.01, || true).unwrap();
    assert_eq!(held.audit.regime, ReclosureRegime::Sticking);
    assert_eq!(held.audit.friction_dissipation_j, 0.0);

    let (cooled, cold_state, cooling) = heated
        .transfer_temperature(retained, &cold, &friction, "T")
        .unwrap();
    assert_eq!(cold_state.occupancy(), original.occupancy());
    assert_eq!(cold_state.point().jump_m(), original.point().jump_m());
    assert_eq!(cold_state.point().damage(), 1.0);
    assert_eq!(cold_state.point().dissipated_j_m2(), point.dissipated_j_m2());
    close(cooling.parameter_work_j_m2, -20.0, 1e-12);
    assert_eq!(cooling.damage_dissipated_j_m2, 0.0);
    assert!(core::ptr::eq(cooled.cohesive_source(), &cold));
    assert!(core::ptr::eq(cooled.friction_source(), &friction));
    assert_eq!(original.point(), point);
    assert!(core::ptr::eq(owner.cohesive_source(), &cold));
}

#[test]
fn temperature_damage_is_counted_once_and_cooling_cannot_heal_it() {
    let card = card([0.3, 0.1]);
    let (cold, friction) = sources(&card);
    let hot = cold.at_temperature(&card, "T", 320.0).unwrap();
    let hot_friction = friction.at_temperature(&card, "T", 320.0).unwrap();
    let owner =
        DamageFrictionLoadLine::new(&cold, &friction, AREA, compliance(), damage_config()).unwrap();
    let jump = [-2e-4, 1.8e-4, 0.0];
    let point = cold.law().bonded_point().trial(jump, 1e-8).unwrap().point;
    assert_eq!(point.damage(), 0.0);
    let original = owner
        .adopt(point, receiving(&cold), ReclosureActivation::AtRest)
        .unwrap();
    let (mut heated, retained, transfer) = owner
        .transfer_temperature(original, &hot, &hot_friction, "T")
        .unwrap();
    let damage = (0.00125 / (0.00125 - 0.00016)) * (1.0 - 0.00016 / 0.00018);
    close(retained.point().damage(), damage, 1e-13);
    assert_eq!(retained.point().jump_m(), jump);
    assert_eq!(retained.point(), transfer.point);
    assert_eq!(retained.occupancy(), original.occupancy());
    assert_eq!(retained.occupancy(), DamageFrictionRegime::Inactive);
    close(transfer.parameter_work_j_m2, 11.9, 1e-12);
    close(transfer.damage_dissipated_j_m2, 8.1 * damage, 1e-12);
    close(
        transfer.stored_change_j_m2,
        transfer.parameter_work_j_m2 - transfer.damage_dissipated_j_m2,
        1e-12,
    );
    close(retained.point().dissipated_j_m2(), 8.1 * damage, 1e-12);
    assert_eq!(heated.config(), owner.config());

    // The next real mechanical audit can select sticking, without repeating
    // the thermal damage or turning the coefficient work into friction heat.
    let force = [-0.6, (1.0 - damage) * 0.09 + 0.02, 0.0];
    let held = heated
        .trial(retained, free_for(jump, force), 0.01, || true)
        .unwrap();
    assert_eq!(held.audit.regime, DamageFrictionRegime::Sticking);
    assert_eq!(held.audit.material.equilibrium.energy.dissipated, 0.0);
    assert_eq!(held.audit.material.viscous_dissipated_j_m2, 0.0);
    assert_eq!(held.audit.friction_dissipation_j, 0.0);
    close(held.audit.damaged_compression_force_n, 0.6 * damage, 2e-9);

    let (_, cooled, cooling) = heated
        .transfer_temperature(retained, &cold, &friction, "T")
        .unwrap();
    assert!(cooled.point().damage() >= retained.point().damage());
    close(cooled.point().damage(), damage, 1e-13);
    close(cooled.point().dissipated_j_m2(), retained.point().dissipated_j_m2(), 1e-12);
    close(cooling.damage_dissipated_j_m2, 0.0, 1e-12);
    assert_eq!(cooled.point().jump_m(), jump);
    assert_eq!(cooled.occupancy(), retained.occupancy());
    assert_eq!(original.point(), point);
}

#[test]
fn partially_damaged_temperature_transfer_preserves_sliding_occupancy_and_new_force() {
    let card = card([0.3, 0.1]);
    let (cold, friction) = sources(&card);
    let hot = cold.at_temperature(&card, "T", 320.0).unwrap();
    let hot_friction = friction.at_temperature(&card, "T", 320.0).unwrap();
    let owner =
        DamageFrictionLoadLine::new(&cold, &friction, AREA, compliance(), damage_config()).unwrap();
    let point = history(&cold, 4e-4);
    close(point.damage(), 0.625, 1e-14);
    let original = owner
        .adopt(point, receiving(&cold), ReclosureActivation::Sliding)
        .unwrap();
    let (mut heated, retained, transfer) = owner
        .transfer_temperature(original, &hot, &hot_friction, "T")
        .unwrap();
    assert_eq!(retained.occupancy(), DamageFrictionRegime::Sliding);
    close(retained.point().damage(), 0.625, 1e-14);
    close(retained.point().dissipated_j_m2(), point.dissipated_j_m2(), 1e-12);
    assert_eq!(transfer.damage_dissipated_j_m2, 0.0);
    assert!(core::ptr::eq(heated.friction_source(), &hot_friction));
    assert!(core::ptr::eq(retained.point().law(), hot.law()));

    let slip = 1e-6;
    let cap = 0.1 * 0.625 * 0.6;
    let end = [JUMP[0], slip, 0.0];
    let expected = [-0.6, cap + 0.5 * AREA * 5e8 * (1.0 - 0.625) * slip, 0.0];
    let free = free_for(end, expected);
    let sliding = heated.trial(retained, free, 0.01, || true).unwrap();
    assert_eq!(sliding.audit.regime, DamageFrictionRegime::Sliding);
    for i in 0..3 {
        close(sliding.mean_force_n[i], expected[i], 2e-9);
        close(sliding.audit.state.point().jump_m()[i], end[i], 2e-12);
    }
    close(sliding.audit.friction_dissipation_j, cap * slip, 2e-13);
    assert_eq!(sliding.audit.material.equilibrium.energy.dissipated, 0.0);
    let at_rest = heated
        .adopt(retained.point(), receiving(&hot), ReclosureActivation::AtRest)
        .unwrap();
    let held = heated.trial(at_rest, free, 0.01, || true).unwrap();
    assert_eq!(held.audit.regime, DamageFrictionRegime::Sticking);
    assert_eq!(held.audit.friction_dissipation_j, 0.0);
}

#[test]
fn temperature_refusals_preserve_both_original_source_owners_and_states() {
    let card = card([0.3, 0.1]);
    let (cold, friction) = sources(&card);
    let hot = cold.at_temperature(&card, "T", 320.0).unwrap();
    let hot_friction = friction.at_temperature(&card, "T", 320.0).unwrap();
    let foreign_card = self::card([0.4, 0.1]);
    let (foreign_cold, foreign_friction) = sources(&foreign_card);
    let foreign_hot = foreign_cold.at_temperature(&foreign_card, "T", 320.0).unwrap();
    let foreign_hot_friction = foreign_friction.at_temperature(&foreign_card, "T", 320.0).unwrap();
    let reclosure =
        ReclosureLoadLine::new(&cold, &friction, AREA, compliance(), reclosure_config()).unwrap();
    let mut damage =
        DamageFrictionLoadLine::new(&cold, &friction, AREA, compliance(), damage_config()).unwrap();
    let failed = reclosure
        .adopt(history(&cold, 0.002), receiving(&cold), ReclosureActivation::Sliding)
        .unwrap();
    let partial = damage
        .adopt(history(&cold, 4e-4), receiving(&cold), ReclosureActivation::Sliding)
        .unwrap();
    for (next, next_friction, axis) in [
        (&hot, &friction, "T"),
        (&cold, &hot_friction, "T"),
        (&foreign_hot, &foreign_hot_friction, "T"),
        (&hot, &hot_friction, "pressure"),
    ] {
        assert!(matches!(
            reclosure.transfer_temperature(failed, next, next_friction, axis),
            Err(ReclosureError::Source(_))
        ));
        assert!(matches!(
            damage.transfer_temperature(partial, next, next_friction, axis),
            Err(DamageFrictionError::Source(_))
        ));
    }
    let other_reclosure =
        ReclosureLoadLine::new(&cold, &friction, 2.0 * AREA, compliance(), reclosure_config()).unwrap();
    assert!(matches!(
        other_reclosure.transfer_temperature(failed, &hot, &hot_friction, "T"),
        Err(ReclosureError::StateMismatch)
    ));
    let other_damage =
        DamageFrictionLoadLine::new(&cold, &friction, 2.0 * AREA, compliance(), damage_config()).unwrap();
    assert!(matches!(
        other_damage.transfer_temperature(partial, &hot, &hot_friction, "T"),
        Err(DamageFrictionError::StateMismatch)
    ));

    // Unchanged source queries are a zero-work rebind, not a new history.
    let (_, same_failed, same) = reclosure
        .transfer_temperature(failed, &cold, &friction, "T")
        .unwrap();
    assert_eq!(same_failed, failed);
    assert_eq!(same.parameter_work_j_m2, 0.0);
    assert_eq!(same.damage_dissipated_j_m2, 0.0);
    let (_, same_partial, same) = damage
        .transfer_temperature(partial, &cold, &friction, "T")
        .unwrap();
    assert_eq!(same_partial, partial);
    assert_eq!(same.parameter_work_j_m2, 0.0);
    assert_eq!(same.damage_dissipated_j_m2, 0.0);

    // Both old owners remain mechanically usable after every rejected stage.
    let free = free_for(JUMP, [-0.4, 0.0, 0.0]);
    assert_eq!(
        reclosure.trial(failed, free, 0.01, || true).unwrap().audit.friction_dissipation_j,
        0.0
    );
    assert_eq!(
        damage.trial(partial, free, 0.01, || true).unwrap().audit.friction_dissipation_j,
        0.0
    );
}

#[test]
fn reclosure_temperature_rechecks_original_cross_axis_friction_admission() {
    let card = card([10.0, 0.1]);
    let (cold, friction) = sources(&card);
    let hot = cold.at_temperature(&card, "T", 320.0).unwrap();
    let hot_friction = friction.at_temperature(&card, "T", 320.0).unwrap();
    let owner =
        ReclosureLoadLine::new(&cold, &friction, AREA, compliance(), reclosure_config()).unwrap();
    let point = history(&cold, 0.002);
    let original = owner
        .adopt(point, receiving(&cold), ReclosureActivation::Sliding)
        .unwrap();
    assert!(matches!(
        owner.transfer_temperature(original, &hot, &hot_friction, "T"),
        Err(ReclosureError::CouplingUnestablished { .. })
    ));
    assert!(core::ptr::eq(owner.friction_source(), &friction));
    assert_eq!(original.point(), point);
}
