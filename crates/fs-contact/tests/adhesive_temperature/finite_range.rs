//! Finite-range normal energy and source-owned junction shear through temperature changes.
use super::{close, common, modulus, temperatures};
use core::f64::consts::PI;
use fs_contact::interface_binding::adhesive_shear::SourceAdhesiveShear;
use fs_contact::normal_patch::adhesive::temperature::AdhesiveTemperatureError;
use fs_contact::normal_patch::dmt::{DmtEnvelope, DmtError, DmtRegime, SourceDmtSphereContact};
use fs_contact::normal_patch::maugis::continuation::{MaugisContinuationError, MaugisPhase};
use fs_contact::normal_patch::maugis::{MaugisEnvelope, MaugisError, SourceMaugisSphereContact};
use fs_material::state_point::MaterialPropertySelection;

fn dmt(fixture: &common::Fixture, range: f64) -> SourceDmtSphereContact {
    SourceDmtSphereContact::new(
        fixture.normal.clone(),
        1e-6,
        range,
        DmtEnvelope {
            maximum_indentation_m: 1e-9,
            half_space_depth_m: 1e-3,
            layer_thickness_m: 1e-3,
            maximum_interaction_radius_to_radius: 0.1,
            maximum_interaction_radius_to_depth: 0.1,
            maximum_interaction_radius_to_layer: 0.1,
            maximum_nominal_strain: 0.1,
            maximum_repulsive_pressure_to_yield: 0.3,
            maximum_adhesive_stress_to_yield: 0.1,
            maximum_interaction_range_to_radius: 0.01,
            maximum_tabor_parameter: 0.1,
        },
    )
    .unwrap()
}

fn maugis(fixture: &common::Fixture, lambda: f64) -> SourceMaugisSphereContact {
    let radius: f64 = 0.01;
    let e = fixture.normal.elastic().reduced_modulus_pa();
    let w = fixture.normal.adhesion_energy_j_per_m2();
    let a0 = (PI * w * radius * radius / (4.0 * e / 3.0)).cbrt();
    let d0 = a0 * a0 / radius;
    SourceMaugisSphereContact::new(
        fixture.normal.clone(),
        radius,
        2.0 * d0 / (PI * lambda),
        MaugisEnvelope {
            maximum_indentation_m: 20.0 * d0,
            half_space_depth_m: 0.1,
            layer_thickness_m: 0.1,
            maximum_interaction_radius_to_radius: 0.1,
            maximum_interaction_radius_to_depth: 0.1,
            maximum_interaction_radius_to_layer: 0.1,
            maximum_nominal_strain: 0.1,
            maximum_repulsive_pressure_to_yield: 0.5,
            maximum_adhesive_stress_to_yield: 0.9,
            maximum_interaction_range_to_radius: 0.01,
        },
    )
    .unwrap()
}

#[test]
fn dmt_temperature_retains_contact_gap_and_separated_energy_without_annular_shear() {
    let fixture = common::Fixture::with_options(1e5, 0.01, false, 320.0);
    let law = dmt(&fixture, 1e-9);
    let original = law.clone();
    let shear = SourceAdhesiveShear::resolve(
        law.normal_state(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let e = 1e5 * modulus(310.0, 290.0);
    let w = 0.0055;
    let r: f64 = 1e-6;
    let z: f64 = 1e-9;
    let pull = 2.0 * PI * r * w;
    for delta in [-2.0 * z, -z, -0.5 * z, 0.0, 0.2 * z, z] {
        let update = law
            .refresh_at_held_indentation(delta, fixture.cards(), temperatures(310.0, 290.0, 315.0))
            .unwrap();
        let (force, energy, area) = if delta >= 0.0 {
            let repulsive = (4.0 / 3.0) * e * r.sqrt() * delta.powf(1.5);
            (
                repulsive - pull,
                0.4 * repulsive * delta - pull * (delta + 0.5 * z),
                PI * r * delta,
            )
        } else if delta > -z {
            (
                -pull * (1.0 + delta / z),
                -0.5 * pull * (delta + z).powi(2) / z,
                0.0,
            )
        } else {
            (0.0, 0.0, 0.0)
        };
        close(update.response().normal_force_n, force, 3e-13);
        close(update.response().reversible_energy_j, energy, 3e-13);
        close(update.response().contact_area_m2, area, 3e-13);
        assert_eq!(update.response().regime, update.before().regime);
        assert_eq!(update.response().indentation_m, delta);
        close(
            update.parameter_work_j(),
            energy - update.before().reversible_energy_j,
            3e-13,
        );
        assert_eq!(update.mechanical_work_j(), 0.0);
        let next_shear = shear.rebind_dmt_temperature(&law, &update).unwrap();
        let patch = next_shear
            .sliding_patch_dmt(update.source().normal_state(), update.response())
            .unwrap();
        close(patch.patch().parameters().shear_strength_pa, 1700.0, 3e-15);
        close(
            patch.static_patch().parameters().shear_strength_pa,
            2550.0,
            3e-15,
        );
        let response = patch.evaluate([0.0; 2], 10.0).unwrap();
        close(
            response.spin_moment_nm,
            -(2.0 / 3.0) * 1700.0 * area * update.response().contact_radius_m,
            3e-13,
        );
        if area == 0.0 {
            let mixed = patch.evaluate([0.1, -0.2], 30.0).unwrap();
            assert_eq!(mixed.force_n, [0.0; 2]);
            assert_eq!(mixed.spin_moment_nm, 0.0);
            assert_eq!(mixed.dissipated_power_w, 0.0);
        }
        assert!(
            shear
                .sliding_patch_dmt(update.source().normal_state(), update.response())
                .is_err()
        );
        assert!(
            shear
                .rebind_dmt_temperature(&dmt(&fixture, 1.1e-9), &update)
                .is_err()
        );
    }
    assert_eq!(law, original);
}

#[test]
fn dmt_held_load_changes_real_area_and_keeps_actuator_work_pull_off_and_tabor_controls() {
    let fixture = common::Fixture::with_options(1e5, 0.01, false, 320.0);
    let law = dmt(&fixture, 1e-9);
    let hot = temperatures(310.0, 290.0, 315.0);
    let e = 1e5 * modulus(310.0, 290.0);
    let w = 0.0055;
    let r: f64 = 1e-6;
    for force in [0.0, 1e-8, -1e-8] {
        let update = law
            .refresh_at_held_load(force, fixture.cards(), hot)
            .unwrap();
        let radius = ((3.0 * r / (4.0 * e)) * (force + 2.0 * PI * r * w)).cbrt();
        let delta = radius * radius / r;
        close(update.response().contact_radius_m, radius, 3e-13);
        close(update.response().indentation_m, delta, 3e-13);
        assert_eq!(update.response().regime, DmtRegime::Contact);
        assert!(update.response().contact_area_m2 > 0.0);
        close(
            update.mechanical_work_j(),
            force * (delta - update.before().indentation_m),
            3e-13,
        );
        close(
            update.parameter_work_j(),
            update.potential_change_j() - update.mechanical_work_j(),
            3e-15,
        );
    }
    assert!(matches!(
        law.refresh_at_held_load(
            0.9 * law.force_controlled_pull_off_n(),
            fixture.cards(),
            hot
        ),
        Err(AdhesiveTemperatureError::Dmt(
            DmtError::NoContactEquilibrium { .. }
        ))
    ));
    let e0 = law.normal_state().elastic().reduced_modulus_pa();
    let range = (r * (0.01 / e0).powi(2)).cbrt() / 0.08;
    let near_limit = dmt(&fixture, range);
    assert!(matches!(
        near_limit.at_temperatures(fixture.cards(), temperatures(315.0, 315.0, 280.0)),
        Err(AdhesiveTemperatureError::Dmt(
            DmtError::OutsideApplicability { .. }
        ))
    ));
}

#[test]
fn maugis_temperature_updates_finite_zone_equations_absolute_potential_and_original_shear() {
    let fixture = common::Fixture::with_options(1e5, 1e-7, false, 320.0);
    let law = maugis(&fixture, 1.0);
    let shear = SourceAdhesiveShear::resolve(
        law.normal_state(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let before = law.evaluate_contact_energy(0.0).unwrap();
    let update = law
        .refresh_at_held_indentation(before, fixture.cards(), temperatures(310.0, 290.0, 315.0))
        .unwrap();
    let e = 1e5 * modulus(310.0, 290.0);
    let w = 0.55e-7;
    let r = law.effective_radius_m();
    let a0 = (PI * w * r * r / (4.0 * e / 3.0)).cbrt();
    let d0 = a0 * a0 / r;
    let p0 = PI * w * r;
    let lambda = 2.0 * d0 / (PI * law.interaction_range_m());
    let contact = update.response().contact_response().unwrap();
    let a = contact.contact_radius_m / a0;
    let ratio = contact.cohesive_radius_m / contact.contact_radius_m;
    let u = (ratio * ratio - 1.0).sqrt();
    let theta = (1.0 / ratio).acos();
    let closure = 0.5 * lambda * a * a * ((ratio * ratio - 2.0) * theta + u)
        + (4.0 / 3.0) * lambda * lambda * a * (u * theta - ratio + 1.0);
    close(closure, 1.0, 2e-8);
    close(
        contact.normal_force_n / p0,
        a.powi(3) - lambda * a * a * (u + ratio * ratio * theta),
        2e-8,
    );
    assert!((contact.indentation_m / d0 - (a * a - (4.0 / 3.0) * lambda * a * u)).abs() < 2e-8);
    close(update.source().maugis_parameter(), lambda, 3e-14);
    assert_eq!(update.response().phase(), MaugisPhase::Contact);
    assert_eq!(update.before(), &before);
    assert_eq!(
        update.parameter_work_j(),
        update.response().potential_j() - before.potential_j()
    );
    let next_shear = shear.rebind_maugis_temperature(&law, &update).unwrap();
    let patch = next_shear
        .sliding_patch_maugis_equilibrium(update.source(), update.response())
        .unwrap();
    close(
        patch.contact_area_m2(),
        PI * contact.contact_radius_m.powi(2),
        3e-15,
    );
    assert!(patch.contact_area_m2() < PI * contact.cohesive_radius_m.powi(2));
    let spin = patch.evaluate([0.0; 2], 5.0).unwrap();
    close(
        spin.spin_moment_nm,
        -(2.0 / 3.0) * 1700.0 * patch.contact_area_m2() * contact.contact_radius_m,
        3e-13,
    );
    assert!(
        shear
            .sliding_patch_maugis_equilibrium(update.source(), update.response())
            .is_err()
    );

    // At fixed displacement and elastic modulus, stationarity gives
    // dU/dW = U_surface/W. This checks the signed thermal parameter work against
    // the physical surface term, independently of subtracting two output fields.
    let h = 1e-3;
    let lower = law
        .refresh_at_held_indentation(
            before,
            fixture.cards(),
            temperatures(300.0, 300.0, 300.0 - h),
        )
        .unwrap();
    let upper = law
        .refresh_at_held_indentation(
            before,
            fixture.cards(),
            temperatures(300.0, 300.0, 300.0 + h),
        )
        .unwrap();
    let observed = (upper.parameter_work_j() - lower.parameter_work_j()) / (2.0 * h);
    let expected =
        before.surface_energy_j() / law.normal_state().adhesion_energy_j_per_m2() * -3e-9;
    assert!((observed - expected).abs() < 5e-6 * expected.abs());

    let force = 0.2 * PI * law.normal_state().adhesion_energy_j_per_m2() * r;
    let load = law
        .refresh_at_held_load(force, fixture.cards(), temperatures(310.0, 290.0, 315.0))
        .unwrap();
    assert!((load.response().normal_force_n() - force).abs() < 1e-10 * p0);
    close(
        load.mechanical_work_j(),
        force * (load.response().indentation_m() - load.before().indentation_m()),
        3e-15,
    );
}

#[test]
fn maugis_noncontact_phase_keeps_deformed_gap_energy_zero_shear_and_explicit_refusals() {
    let fixture = common::Fixture::with_options(1e5, 1e-7, false, 320.0);
    let law = maugis(&fixture, 0.2);
    let shear = SourceAdhesiveShear::resolve(
        law.normal_state(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let limits = law.noncontact_limits().unwrap().unwrap();
    let delta = 0.25 * limits.minimum_indentation_m + 0.75 * limits.first_contact_indentation_m;
    let previous = law.evaluate_noncontact(delta).unwrap();
    let update = law
        .refresh_at_held_indentation(previous, fixture.cards(), temperatures(300.0, 300.0, 315.0))
        .unwrap();
    let e = law.normal_state().elastic().reduced_modulus_pa();
    let r = law.effective_radius_m();
    let z = law.interaction_range_m();
    let sigma = 0.55e-7 / z;
    let offset = 4.0 * sigma * r / (PI * e);
    let c = offset + (offset * offset + 2.0 * r * (delta + z)).sqrt();
    let elastic = 8.0 * sigma * sigma * c.powi(3) / (3.0 * e);
    let potential = -PI * sigma * c.powi(4) / (4.0 * r) + 0.5 * elastic;
    assert_eq!(update.response().phase(), MaugisPhase::Noncontact);
    assert!(update.response().contact_response().is_none());
    close(update.response().cohesive_radius_m(), c, 2e-12);
    close(
        update.response().normal_force_n(),
        -PI * sigma * c * c,
        2e-12,
    );
    close(
        update.response().central_gap_m(),
        -delta - 2.0 * sigma * c / e,
        2e-12,
    );
    close(update.response().elastic_energy_j(), elastic, 2e-12);
    close(update.response().potential_j(), potential, 2e-12);
    let next_shear = shear.rebind_maugis_temperature(&law, &update).unwrap();
    let patch = next_shear
        .sliding_patch_maugis_equilibrium(update.source(), update.response())
        .unwrap();
    let response = patch.evaluate([0.1, -0.2], 30.0).unwrap();
    assert_eq!(patch.contact_area_m2(), 0.0);
    assert_eq!(response.force_n, [0.0; 2]);
    assert_eq!(response.spin_moment_nm, 0.0);
    assert_eq!(response.dissipated_power_w, 0.0);
    let noop = law
        .refresh_at_held_indentation(previous, fixture.cards(), temperatures(300.0, 300.0, 300.0))
        .unwrap();
    assert_eq!(noop.source(), &law);
    assert_eq!(noop.before(), &previous);
    assert_eq!(noop.response(), &previous);
    assert_eq!(noop.parameter_work_j(), 0.0);
    let detached = law.evaluate_detached(-2.0 * z).unwrap();
    let distant = law
        .refresh_at_held_indentation(detached, fixture.cards(), temperatures(310.0, 290.0, 315.0))
        .unwrap();
    assert_eq!(distant.response().phase(), MaugisPhase::Detached);
    assert_eq!(distant.parameter_work_j(), 0.0);

    let foreign = maugis(&fixture, 0.3);
    assert!(matches!(
        foreign.refresh_at_held_indentation(
            previous,
            fixture.cards(),
            temperatures(300.0, 300.0, 315.0)
        ),
        Err(AdhesiveTemperatureError::MaugisContinuation(
            MaugisContinuationError::HistoryMismatch
        ))
    ));
    assert!(shear.rebind_maugis_temperature(&foreign, &update).is_err());
    assert!(
        next_shear
            .sliding_patch_maugis_equilibrium(&foreign, update.response())
            .is_err()
    );
    let strong = maugis(&fixture, 5.0);
    let attached = strong
        .evaluate_contact_energy(0.98 * strong.minimum_attached_indentation_m())
        .unwrap();
    assert!(matches!(
        strong.refresh_at_held_indentation(
            attached,
            fixture.cards(),
            temperatures(300.0, 300.0, 320.0)
        ),
        Err(AdhesiveTemperatureError::MaugisContinuation(
            MaugisContinuationError::Contact(MaugisError::NoContactEquilibrium { .. })
        ))
    ));
    assert!(matches!(
        strong.refresh_at_held_load(
            0.9 * strong.force_controlled_pull_off_n(),
            fixture.cards(),
            temperatures(300.0, 300.0, 320.0)
        ),
        Err(AdhesiveTemperatureError::Maugis(
            MaugisError::BelowPullOff { .. }
        ))
    ));
}
