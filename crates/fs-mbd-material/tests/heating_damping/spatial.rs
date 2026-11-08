//! Full native spatial loss enters the same finite source-bound thermal history.

use super::*;
use fs_mbd::Vec3 as RigidVec3;
use fs_mbd::spatial_loads::internal_impulse::InternalWrenchImpulseControls;
use fs_mbd::spatial_loads::internal_impulse::viscous::ViscousConnectorControls;
use fs_mbd_material::thermal::heating::motion::damping::spatial::{
    SpatialDampedRigidPair, SpatialDampingSolve, SpatialResistanceShape,
};

fn controls() -> SpatialDampingSolve {
    SpatialDampingSolve {
        thermal: solve(),
        connector: ViscousConnectorControls {
            impulse: InternalWrenchImpulseControls {
                angular_origin_world_m: RigidVec3::ZERO,
                linear_tolerance_ns: 1e-12,
                angular_tolerance_nms: 1e-12,
                work_tolerance_j: 1e-11,
            },
            linear_speed_tolerance_m_s: 1e-12,
            angular_speed_tolerance_rad_s: 1e-12,
        },
    }
}

fn shape() -> SpatialResistanceShape {
    let mut normalized = [[0.0; 6]; 6];
    normalized[0][0] = 1.0;
    normalized[3][3] = 1.0;
    SpatialResistanceShape::new(normalized, 0.5).unwrap()
}

fn spinning_pair<'a>(
    reference: &'a ResolvedBody,
    solid: &'a MaterialCard,
    chart: &'a ResolvedHeatCapacityCurve,
    interface: &'a InterfaceSystemCard,
    watts: [f64; 2],
) -> DampedRigidPair<'a> {
    let (_, mut damper, _) = pair(reference, solid, chart, interface, 0.0, 0.0, watts).into_parts();
    damper.first_heat_fraction = 0.5;
    let bodies = [0, 1].map(|i| {
        let sign = if i == 0 { 1.0 } else { -1.0 };
        HeatedRigidBody::new(
            heated(reference, solid, chart, watts[i], vec![]),
            Se3::from_parts(So3::identity(), Vec3::new(-sign, 0.0, 0.0)).unwrap(),
            Twist::new(Vec3::new(sign, 0.0, 0.0), Vec3::new(sign, 0.0, 0.0)),
            RigidBodyIntegrator::new(Gravity::ZERO),
            HeatedMotionControls {
                maximum_step_s: 0.125,
                relative_tensor_tolerance: 1e-12,
                transfer: ConfigurationTransferControls {
                    relative_momentum_tolerance: 1e-12,
                    maximum_parameter_work_j: 1.0,
                },
                maximum_numerical_energy_change_j: 1e-7,
            },
            || Ok(()),
        )
        .unwrap()
    });
    DampedRigidPair::new(bodies, damper).unwrap()
}

fn world_spin(state: fs_mbd::RigidBodyState) -> RigidVec3 {
    state
        .pose()
        .orientation()
        .rotate_body_to_world(state.angular_momentum_body())
}

fn assert_unchanged(old: &SpatialDampedRigidPair<'_>, saved: &SpatialDampedRigidPair<'_>) {
    let (old_pair, old_shape) = old.clone().into_parts();
    let (saved_pair, saved_shape) = saved.clone().into_parts();
    unchanged(&old_pair, &saved_pair);
    assert_eq!(old_shape, saved_shape);
}

#[test]
fn independent_translation_and_spin_heat_enter_the_original_finite_cells_once() {
    let solid = card(10.0, 0.0, 500.0);
    let chart = caloric(&solid, 1e5, SolidLiquidPhase::Solid);
    let reference = reference(&solid, 1);
    let interface = interface(&solid, 0.0, 500.0, 0.0);
    let (bodies, mut damper, _) =
        spinning_pair(&reference, &solid, &chart, &interface, [2.0, 3.0]).into_parts();
    damper.first_heat_fraction = 0.25;
    let original = DampedRigidPair::new(bodies, damper)
        .unwrap()
        .with_spatial_resistance(shape());
    let saved = original.clone();
    let (next, audit) = original
        .step_to(0.1, controls(), &[0.3, 0.6], || false)
        .unwrap();

    // Two original 2 kg cubes have I_x=1/3 kg m². With ell=1/2 m,
    // linear slip decays at 1/s and spin slip at 1.5/s.
    let linear = (-0.1_f64).exp();
    let spin = (-0.15_f64).exp();
    let loss = 2.0 * (1.0 - linear * linear) + (1.0 - spin * spin) / 3.0;
    for (i, body) in next.bodies().iter().enumerate() {
        let sign = if i == 0 { 1.0 } else { -1.0 };
        close(body.state().linear_momentum_world().x, sign * 2.0 * linear);
        close(world_spin(body.state()).x, sign * spin / 3.0);
        let deposited = [0.25, 0.75][i] * loss;
        close(audit.heat_by_body_j[i], deposited);
        close(body.thermal().energy_j(), [0.2, 0.3][i] + deposited);
        close(
            body.thermal().temperature_k(CellId(0)).unwrap(),
            300.0 + ([0.2, 0.3][i] + deposited) / 20.0,
        );
        close(body.thermal().time_s(), 0.1);
    }
    close(audit.dissipated_heat_j, loss);
    close(audit.impulse.impulse.kinetic_energy_change_j, -loss);
    close(audit.impulse.kinetic_heat_residual_j, 0.0);
    close(audit.external_heat_j, 0.5);
    close(audit.combined_residual_j, 0.0);
    assert_eq!(audit.parameter_work_j, [0.0; 2]);
    assert_unchanged(&original, &saved);
}

#[test]
fn finite_heat_capacity_changes_source_viscosity_and_both_momentum_channels() {
    let mut observed = Vec::new();
    for cp in [1.0, 10.0] {
        let solid = card(cp, 0.0, 500.0);
        let chart = caloric(&solid, 1e5, SolidLiquidPhase::Solid);
        let reference = reference(&solid, 1);
        let interface = interface(&solid, 1.0, 500.0, 0.0);
        let old = spinning_pair(&reference, &solid, &chart, &interface, [0.0; 2])
            .with_spatial_resistance(shape());
        let (next, audit) = old.step_to(0.1, controls(), &[0.2; 2], || false).unwrap();

        // Independent scalar caloric root: equal heat routing gives
        // T_A=300+Q/(4Cp), hence c(T_A)=1+Q/(4Cp).
        let loss_at = |heat: f64| {
            let c = 1.0 + heat / (4.0 * cp);
            2.0 * -(-0.2 * c).exp_m1() + -(-0.3 * c).exp_m1() / 3.0
        };
        let (mut low, mut high) = (0.0, 7.0 / 3.0);
        for _ in 0..80 {
            let middle = 0.5 * (low + high);
            if middle < loss_at(middle) {
                low = middle;
            } else {
                high = middle;
            }
        }
        let expected_heat = 0.5 * (low + high);
        let c = 1.0 + expected_heat / (4.0 * cp);
        close(audit.dissipated_heat_j, expected_heat);
        close(audit.source.damping_n_s_m(), c);
        close(
            next.bodies()[0].state().linear_momentum_world().x,
            2.0 * (-0.1 * c).exp(),
        );
        close(
            world_spin(next.bodies()[0].state()).x,
            (-0.15 * c).exp() / 3.0,
        );
        close(audit.combined_residual_j, 0.0);
        assert_eq!(next.damper().card.content_hash(), interface.content_hash());
        observed.push((
            audit.dissipated_heat_j,
            next.bodies()[0].state().linear_momentum_world().x,
        ));
    }
    assert!(observed[0].0 > observed[1].0);
    assert!(observed[0].1 < observed[1].1);
}

#[test]
fn expansion_preserves_spatial_recoil_and_separates_rotational_parameter_work() {
    let solid = card(10.0, 0.01, 500.0);
    let chart = caloric(&solid, 1e5, SolidLiquidPhase::Solid);
    let reference = reference(&solid, 1);
    let interface = interface(&solid, 0.0, 500.0, 0.0);
    let old = spinning_pair(&reference, &solid, &chart, &interface, [0.0; 2])
        .with_spatial_resistance(shape());
    let (next, audit) = old.step_to(0.1, controls(), &[0.2; 2], || false).unwrap();
    for i in 0..2 {
        let original = old.bodies()[i].principal_frame().properties();
        let expanded = next.bodies()[i].principal_frame().properties();
        assert_eq!(original.mass(), expanded.mass());
        assert!(expanded.principal_inertia_body().x > original.principal_inertia_body().x);
        let impulse_state = audit.impulse.impulse.bodies[i].state_after;
        let physical_spin = world_spin(impulse_state);
        close(world_spin(next.bodies()[i].state()).x, physical_spin.x);
        close(
            next.bodies()[i].state().linear_momentum_world().x,
            impulse_state.linear_momentum_world().x,
        );
        let expected_work = 0.5
            * physical_spin.x
            * physical_spin.x
            * (1.0 / expanded.principal_inertia_body().x
                - 1.0 / original.principal_inertia_body().x);
        close(audit.parameter_work_j[i], expected_work);
        assert!(audit.parameter_work_j[i] < 0.0);
        close(
            next.bodies()[i].thermal().energy_j(),
            audit.heat_by_body_j[i],
        );
    }
    close(audit.combined_residual_j, 0.0);

    // The next original source uses the accepted expanded inertia rather than
    // reconstructing a new state at the previous angular velocity.
    let inertia = next.bodies()[0]
        .principal_frame()
        .properties()
        .principal_inertia_body()
        .x;
    let spin_before = world_spin(next.bodies()[0].state()).x;
    let (continued, _) = next.step_to(0.2, controls(), &[0.2; 2], || false).unwrap();
    close(
        world_spin(continued.bodies()[0].state()).x,
        spin_before * (-0.1 * 0.5 / inertia).exp(),
    );
}

#[test]
fn zero_shape_and_late_source_budget_or_cancellation_refusal_preserve_both_histories() {
    let solid = card(10.0, 0.0, 500.0);
    let chart = caloric(&solid, 1e5, SolidLiquidPhase::Solid);
    let reference = reference(&solid, 1);
    let interface = interface(&solid, 0.0, 500.0, 0.0);
    let original = spinning_pair(&reference, &solid, &chart, &interface, [2.0, 3.0]);
    let zero = original
        .clone()
        .with_spatial_resistance(SpatialResistanceShape::new([[0.0; 6]; 6], 0.5).unwrap());
    let (idle, receipt) = zero
        .step_to(0.1, controls(), &[0.2, 0.3], || false)
        .unwrap();
    assert_eq!(receipt.dissipated_heat_j, 0.0);
    assert_eq!(receipt.heat_by_body_j, [0.0; 2]);
    for i in 0..2 {
        close(
            idle.bodies()[i].state().linear_momentum_world().x,
            original.bodies()[i].state().linear_momentum_world().x,
        );
        close(
            world_spin(idle.bodies()[i].state()).x,
            world_spin(original.bodies()[i].state()).x,
        );
        close(idle.bodies()[i].thermal().energy_j(), [0.2, 0.3][i]);
    }
    let old = original.with_spatial_resistance(shape());
    let saved = old.clone();
    let mut calls = 0;
    let (expected, _) = old
        .step_to(0.1, controls(), &[0.3, 0.6], || {
            calls += 1;
            false
        })
        .unwrap();
    let mut observed = 0;
    assert!(matches!(
        old.step_to(0.1, controls(), &[0.3, 0.6], || {
            observed += 1;
            observed == calls
        }),
        Err(DampingError::Thermal(DepositionError::Cancelled))
    ));
    let mut limited = controls();
    limited.thermal.maximum_relaxation_step = 1e-8;
    assert!(matches!(
        old.step_to(0.1, limited, &[0.3, 0.6], || false),
        Err(DampingError::Budget(_))
    ));
    let mut invalid = controls();
    invalid.connector.impulse.work_tolerance_j = f64::NAN;
    assert!(old.step_to(0.1, invalid, &[0.3, 0.6], || false).is_err());
    assert_unchanged(&old, &saved);
    let (retry, _) = old.step_to(0.1, controls(), &[0.3, 0.6], || false).unwrap();
    for i in 0..2 {
        assert_eq!(retry.bodies()[i].state(), expected.bodies()[i].state());
        assert_eq!(
            retry.bodies()[i].thermal().cells(),
            expected.bodies()[i].thermal().cells()
        );
    }
}
