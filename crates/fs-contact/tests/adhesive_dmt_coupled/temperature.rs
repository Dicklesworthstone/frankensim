//! G1/G3/G4: original joint histories survive prescribed source temperatures.
use super::*;
#[path = "../adhesive_temperature/common.rs"]
mod thermal;
use fs_contact::interface_binding::adhesive_shear::SourceAdhesiveShear;
use fs_contact::normal_patch::adhesive::temperature::AdhesiveTemperatures;
use fs_material::state_point::MaterialPropertySelection;

fn temperatures() -> AdhesiveTemperatures {
    AdhesiveTemperatures {
        surface_a_k: 310.0,
        surface_b_k: 290.0,
        interface_k: 315.0,
    }
}

#[test]
fn source_temperature_preserves_static_sliding_and_released_states_and_original_compliance() {
    let fixture = thermal::Fixture::with_options(1e5, 0.01, true, 320.0);
    let normal = fixture.dmt(Z);
    let shear = SourceAdhesiveShear::resolve(
        normal.normal_state(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let owner = DmtCoupledLoadLine::new(&normal, &shear, compliance(), config()).unwrap();
    for (delta, velocity, regime) in [
        (0.3 * Z, [0.0; 2], DmtCoupledRegime::Sticking),
        (0.3 * Z, [1e-3, -2e-3], DmtCoupledRegime::Sliding),
        (-0.4 * Z, [1e-3, -2e-3], DmtCoupledRegime::Released),
        (-2.0 * Z, [0.0; 2], DmtCoupledRegime::Released),
    ] {
        let state = owner
            .initialize_with_velocity([delta, 0.2 * Z, -0.7 * Z], velocity)
            .unwrap();
        let refreshed = normal
            .refresh_at_held_indentation(delta, fixture.cards(), temperatures())
            .unwrap();
        let hot_shear = shear.rebind_dmt_temperature(&normal, &refreshed).unwrap();
        let (next, successor) = owner
            .with_temperature_refresh(state, &refreshed, &hot_shear, &CancelGate::new_clock_free())
            .unwrap();
        assert_eq!(successor.position_m(), state.position_m());
        assert_eq!(successor.regime(), regime);
        assert_eq!(next.compliance_m_per_n(), owner.compliance_m_per_n());
        assert_eq!(next.config(), owner.config());
        assert!(core::ptr::eq(next.normal(), refreshed.source()));
        assert!(core::ptr::eq(next.shear(), &hot_shear));
        let path = next.path(delta, delta).unwrap();
        close(
            path.mean_real_area_m2,
            core::f64::consts::PI * 1e-6 * delta.max(0.0),
            1e-30,
            1e-14,
        );
        close(
            path.static_capacity_n,
            2550.0 * path.mean_real_area_m2,
            1e-30,
            1e-14,
        );
        close(
            path.kinetic_capacity_n,
            1700.0 * path.mean_real_area_m2,
            1e-30,
            1e-14,
        );
        assert_eq!(path.normal_work_j, 0.0);
        if regime == DmtCoupledRegime::Released {
            assert_eq!(path.static_capacity_n, 0.0);
            assert_eq!(path.kinetic_capacity_n, 0.0);
        }
        // Both old and new state tokens are tied to their actual source owner.
        assert!(owner.with_endpoint_velocity(successor, velocity).is_err());
        assert!(next.with_endpoint_velocity(state, velocity).is_err());
        assert_eq!(
            owner.with_endpoint_velocity(state, velocity).unwrap(),
            state
        );
    }
    let state = owner
        .initialize_with_velocity([0.2 * Z, 0.0, 0.0], [1e-3, 0.0])
        .unwrap();
    let unchanged = normal
        .refresh_at_held_indentation(
            0.2 * Z,
            fixture.cards(),
            AdhesiveTemperatures {
                surface_a_k: 300.0,
                surface_b_k: 300.0,
                interface_k: 300.0,
            },
        )
        .unwrap();
    let same_shear = shear.rebind_dmt_temperature(&normal, &unchanged).unwrap();
    let (_, same_state) = owner
        .with_temperature_refresh(
            state,
            &unchanged,
            &same_shear,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(same_state, state);
    assert_eq!(unchanged.parameter_work_j(), 0.0);
}

#[test]
fn thermal_joint_transfer_rejects_stale_indentation_shear_control_and_cancellation() {
    let fixture = thermal::Fixture::with_options(1e5, 0.01, false, 320.0);
    let normal = fixture.dmt(Z);
    let shear = SourceAdhesiveShear::resolve(
        normal.normal_state(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let owner = DmtCoupledLoadLine::new(&normal, &shear, compliance(), config()).unwrap();
    let state = owner
        .initialize_with_velocity([0.2 * Z, 0.0, 0.0], [1e-3, 0.0])
        .unwrap();
    let valid = normal
        .refresh_at_held_indentation(0.2 * Z, fixture.cards(), temperatures())
        .unwrap();
    let hot_shear = shear.rebind_dmt_temperature(&normal, &valid).unwrap();
    let stale = normal
        .refresh_at_held_indentation(0.3 * Z, fixture.cards(), temperatures())
        .unwrap();
    let force = normal.evaluate(0.2 * Z).unwrap().normal_force_n;
    let held_load = normal
        .refresh_at_held_load(force, fixture.cards(), temperatures())
        .unwrap();
    for (refresh, next_shear) in [
        (&stale, &hot_shear),
        (&valid, &shear),
        (&held_load, &hot_shear),
    ] {
        assert!(matches!(
            owner.with_temperature_refresh(
                state,
                refresh,
                next_shear,
                &CancelGate::new_clock_free()
            ),
            Err(DmtCoupledError::StateMismatch)
        ));
    }
    let gate = CancelGate::new_clock_free();
    gate.request();
    assert!(matches!(
        owner.with_temperature_refresh(state, &valid, &hot_shear, &gate),
        Err(DmtCoupledError::Cancelled)
    ));
    assert_eq!(
        owner.with_endpoint_velocity(state, [1e-3, 0.0]).unwrap(),
        state
    );
    assert!(
        owner
            .with_temperature_refresh(state, &valid, &hot_shear, &CancelGate::new_clock_free())
            .is_ok()
    );
}
