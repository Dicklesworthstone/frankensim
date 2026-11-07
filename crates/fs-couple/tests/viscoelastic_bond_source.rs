//! Original ordered-interface spectra reach the native coupled load line.
use fs_couple::cohesive_contact::{CohesiveContactConfig, CohesiveLoadLine};
use fs_matdb::QueryPoint;
use fs_material::cohesive::relaxation::source::*;
use fs_material::state_point::MaterialPropertySelection;
#[path = "viscoelastic_bond_source/common.rs"]
mod common;
use common::{card, limits};

#[test]
fn sourced_relaxation_drives_native_bonded_motion_and_retains_domain_refusals() {
    let card = card(false, false);
    let query = QueryPoint::new().with("T", 300.0).unwrap();
    let source = resolve_viscoelastic_bond(
        &card,
        &query,
        MaterialPropertySelection::SingleClaimOnly,
        limits(),
    )
    .unwrap();
    assert_eq!(source.resolved().card_identity(), card.content_hash());
    assert_eq!(source.resolved().properties().len(), 9);
    let pins = source
        .resolved()
        .properties()
        .iter()
        .map(|p| {
            (
                p.requirement().name().to_owned(),
                p.answer().receipt.selected,
            )
        })
        .collect::<Vec<_>>();
    let pinned = resolve_viscoelastic_bond(
        &card,
        &query,
        MaterialPropertySelection::PinnedByProperty(pins.clone()),
        limits(),
    )
    .unwrap();
    assert_eq!(pinned.law(), source.law());
    for (a, b) in pinned
        .resolved()
        .properties()
        .iter()
        .zip(source.resolved().properties())
    {
        assert_eq!(a.answer().receipt.selected, b.answer().receipt.selected);
    }
    let mut foreign_pins = pins;
    foreign_pins.push(("unowned_relaxation_coefficient".into(), foreign_pins[0].1));
    assert!(
        resolve_viscoelastic_bond(
            &card,
            &query,
            MaterialPropertySelection::PinnedByProperty(foreign_pins),
            limits()
        )
        .is_err()
    );
    assert_eq!(source.surfaces()[0], *card.surface_a());
    let dt = 0.001;
    let area = 1e-4;
    let compliance = 0.5 * dt * dt;
    let config = CohesiveContactConfig {
        max_iterations: 8,
        max_material_evaluations: 16,
        max_backtracks: 8,
        jump_absolute_tolerance_m: 1e-14,
        jump_relative_tolerance: 1e-11,
        maximum_jump_increment_m: 1e-4,
        maximum_force_n: 100.0,
        material_energy_tolerance_j_m2: 1e-9,
    };
    let mut solver = CohesiveLoadLine::new(
        area,
        [
            [compliance, 0.0, 0.0],
            [0.0, compliance, 0.0],
            [0.0, 0.0, compliance],
        ],
        config,
    )
    .unwrap();
    let prepared = source.law().prepare(dt).unwrap();
    let mut point = source.law().bonded_point();
    let mut velocity = [0.0; 3];
    let load = [0.5, -0.2, 0.1];
    let mut work = 0.0;
    let mut heat = 0.0;
    for _ in 0..1000 {
        let old = point.jump_m();
        let free = core::array::from_fn(|i| old[i] + dt * velocity[i] + compliance * load[i]);
        let step = solver
            .trial_relaxation(&point, &prepared, free, || true)
            .unwrap();
        let new = step.material.point.jump_m();
        for i in 0..3 {
            work += load[i] * (new[i] - old[i]);
            velocity[i] = 2.0 * (new[i] - old[i]) / dt - velocity[i];
        }
        heat += area * step.material.viscous_dissipated_j_m2;
        point = step.material.point;
    }
    assert!(heat > 0.0);
    let kinetic = 0.5 * velocity.iter().map(|v| v * v).sum::<f64>();
    assert!((work - kinetic - area * point.stored_energy_j_m2() - heat).abs() < 1e-10);
    assert!((point.jump_m()[0] - 5e-6).abs() < 1e-9);
    assert!(
        point
            .trial(&source.law().prepare(1.1).unwrap(), point.jump_m(), 1e-9)
            .is_err()
    );
    assert!(
        resolve_viscoelastic_bond(
            &card,
            &QueryPoint::new().with("T", 321.0).unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
            limits()
        )
        .is_err()
    );
}
#[test]
fn incomplete_and_dimensionally_wrong_interface_spectra_do_not_become_bulk_defaults() {
    let query = QueryPoint::new().with("T", 300.0).unwrap();
    for interface in [card(true, false), card(false, true)] {
        assert!(
            resolve_viscoelastic_bond(
                &interface,
                &query,
                MaterialPropertySelection::SingleClaimOnly,
                limits()
            )
            .is_err()
        );
    }
}

#[test]
fn pinned_temperature_refresh_preserves_memory_before_native_motion_uses_hot_spectra() {
    let card = common::temperature_card();
    let source = resolve_viscoelastic_bond(
        &card,
        &QueryPoint::new().with("T", 280.0).unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
        limits(),
    )
    .unwrap();
    let hot = source.at_temperature(&card, "T", 320.0).unwrap();
    source.require_temperature_successor(&hot, "T").unwrap();
    assert_eq!(source.surfaces(), hot.surfaces());
    assert_eq!(
        source.resolved().card_identity(),
        hot.resolved().card_identity()
    );
    for (a, b) in source
        .resolved()
        .properties()
        .iter()
        .zip(hot.resolved().properties())
    {
        assert_eq!(a.requirement(), b.requirement());
        assert_eq!(a.answer().receipt.selected, b.answer().receipt.selected);
    }
    let config = CohesiveContactConfig {
        max_iterations: 8,
        max_material_evaluations: 16,
        max_backtracks: 8,
        jump_absolute_tolerance_m: 1e-14,
        jump_relative_tolerance: 1e-11,
        maximum_jump_increment_m: 1e-4,
        maximum_force_n: 100.0,
        material_energy_tolerance_j_m2: 1e-9,
    };
    let mut native = CohesiveLoadLine::new(
        1e-4,
        [[1e-8, 2e-9, 0.0], [2e-9, 8e-9, 1e-9], [0.0, 1e-9, 7e-9]],
        config,
    )
    .unwrap();
    let dt = 0.001;
    let cold_ramp = source.law().prepare(dt).unwrap();
    let mut before = source.law().bonded_point();
    for free in [
        [1e-5, -5e-6, 3e-6],
        [-3e-6, 2e-6, 8e-6],
        [4e-6, -6e-6, 2e-6],
    ] {
        before = native
            .trial_relaxation(&before, &cold_ramp, free, || true)
            .unwrap()
            .material
            .point;
    }
    let transfer = before.transfer_temperature(hot.law(), 1e-10).unwrap();
    assert_eq!(transfer.point.jump_m(), before.jump_m());
    assert_eq!(transfer.point.elapsed_s(), before.elapsed_s());
    assert!(transfer.parameter_work_j_m2 < 0.0);
    assert!(
        (transfer.stored_change_j_m2
            - (transfer.point.stored_energy_j_m2() - before.stored_energy_j_m2()))
        .abs()
            < 1e-12
    );
    assert!(transfer.defect_j_m2.abs() < 1e-12);
    let hot_ramp = hot.law().prepare(dt).unwrap();
    let free = [5e-6, -4e-6, 1e-6];
    let continued = native
        .trial_relaxation(&transfer.point, &hot_ramp, free, || true)
        .unwrap();
    let cold = native
        .trial_relaxation(&before, &cold_ramp, free, || true)
        .unwrap();
    let reset = native
        .trial_relaxation(&hot.law().bonded_point(), &hot_ramp, free, || true)
        .unwrap();
    assert_ne!(continued.material.traction_pa, cold.material.traction_pa);
    assert_ne!(continued.material.traction_pa, reset.material.traction_pa);
    assert_eq!(
        continued.material.point.elapsed_s(),
        before.elapsed_s() + dt
    );
    assert!(continued.material.viscous_dissipated_j_m2 > 0.0);
    assert!(
        (continued.material.work_j_m2
            - continued.material.stored_change_j_m2
            - continued.material.viscous_dissipated_j_m2)
            .abs()
            < 1e-10
    );
}

#[test]
fn temperature_refresh_refuses_untyped_axis_source_replacement_and_geometry_change() {
    let query = QueryPoint::new().with("T", 280.0).unwrap();
    let interface = common::temperature_card();
    let source = resolve_viscoelastic_bond(
        &interface,
        &query,
        MaterialPropertySelection::SingleClaimOnly,
        limits(),
    )
    .unwrap();
    for temperature in [279.0, 321.0, 0.0, f64::NAN] {
        assert!(source.at_temperature(&interface, "T", temperature).is_err());
    }
    assert!(
        source
            .at_temperature(&interface, "missing_axis", 300.0)
            .is_err()
    );
    let untyped = resolve_viscoelastic_bond(
        &card(false, false),
        &query,
        MaterialPropertySelection::SingleClaimOnly,
        limits(),
    )
    .unwrap();
    assert!(
        untyped
            .at_temperature(&card(false, false), "T", 300.0)
            .is_err()
    );
    assert!(source.require_temperature_successor(&untyped, "T").is_err());
    let mut other_limits = limits();
    other_limits.maximum_elapsed_s = 1.0;
    let restricted = resolve_viscoelastic_bond(
        &interface,
        &query,
        MaterialPropertySelection::SingleClaimOnly,
        other_limits,
    )
    .unwrap();
    assert!(
        source
            .require_temperature_successor(&restricted, "T")
            .is_err()
    );
    let expanding = common::expanding_temperature_card();
    let changing_geometry = resolve_viscoelastic_bond(
        &expanding,
        &query,
        MaterialPropertySelection::SingleClaimOnly,
        limits(),
    )
    .unwrap();
    assert!(
        changing_geometry
            .at_temperature(&expanding, "T", 320.0)
            .is_err()
    );
    assert!(source.at_temperature(&expanding, "T", 300.0).is_err());
    // Every failed requery leaves the original claim plan and point retryable.
    let retry = source.at_temperature(&interface, "T", 320.0).unwrap();
    assert_eq!(retry.law().thickness_m(), source.law().thickness_m());
    assert_eq!(retry.law().normal_spectrum().terms[0].1, 0.002);
}
