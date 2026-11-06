//! Native actual slip, capped pressure and original per-body wear stay coupled.

#[path = "sliding/temperature.rs"]
mod temperature;

#[path = "sliding/coefficients.rs"]
mod coefficients;

use super::*;
use fs_contact::interface_binding::thermal_friction::{
    SourceThermalCoulomb, SourceThermalCoulombError, ThermalFrictionAxes,
};
use fs_contact::normal_patch::plastic_source::wear::sliding::{
    PlasticSlidingInput, SourcePlasticSlidingError, SourcePlasticSlidingWear,
};
use fs_contact::normal_patch::{
    FiniteGapPlasticContinuationKind, FiniteGapPlasticDrive, FiniteGapPlasticLoadControls,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Pressure, QuantitySpec};
use fs_tribo::{FrictionRegime, TangentialSolveConfig};

fn friction_controls() -> TangentialSolveConfig {
    TangentialSolveConfig {
        max_iterations: 64,
        maximum_force_n: 1.0,
        force_absolute_tolerance_n: 1.0e-14,
        force_relative_tolerance: 1.0e-12,
    }
}

fn normal_controls() -> FiniteGapPlasticLoadControls {
    FiniteGapPlasticLoadControls {
        minimum_approach_m: -1.0e-4,
        maximum_approach_m: 1.0e-4,
        absolute_force_tolerance_n: 1.0e-14,
        relative_force_tolerance: 1.0e-12,
        approach_tolerance_m: 1.0e-18,
        maximum_equilibrium_evaluations: 128,
    }
}

fn motion() -> Vec<PlasticSlidingInput> {
    vec![
        PlasticSlidingInput {
            released_increment_m: [0.0; 2],
            compliance_m_per_n: 0.0,
        };
        25
    ]
}

fn source(
    hardness_pa: [f64; 2],
    maximum_recession_m: [f64; 2],
    two_cells: bool,
    pressure_domain_pa: Option<[f64; 2]>,
    cx: &Cx<'_>,
) -> SourcePlasticSlidingWear {
    // Reuse the existing real material-card, hardness-test and Archard fixture;
    // only its immutable interface is rebuilt with the two friction claims.
    let original = fixture(
        [0.01, 0.03],
        hardness_pa,
        maximum_recession_m,
        two_cells,
        cx,
    );
    let normal = original.contact().normal_source();
    let old_card = normal.elastic().interface().card();
    let mut claims = old_card.claims().clone();
    let pressure_quantity =
        QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static));
    for (name, value) in [("mu_static", 0.6), ("mu_kinetic", 0.4)] {
        let mut property = claim(PropertyKey::new(name, Dims::NONE), value);
        if let Some([lo, hi]) = pressure_domain_pa {
            property.validity = property
                .validity
                .with_quantity("p", pressure_quantity, lo, hi);
            property.value = PropertyValue::Curve {
                abscissa: "p".into(),
                abscissa_dims: Pressure::DIMS,
                knots: vec![(lo, value), (hi, 0.5 * value)],
                dims: Dims::NONE,
            };
            property.interpolation = InterpolationPolicy::LinearInside;
        }
        claims.insert_claim(property).unwrap();
    }
    let card = InterfaceSystemCard::assemble(
        old_card.surface_a().clone(),
        old_card.surface_b().clone(),
        old_card.context().clone(),
        claims,
        old_card.models().to_vec(),
    )
    .unwrap();
    let dry = bind_dry_interface_system_card(&card, InputAuthority::SyntheticFixture).unwrap();
    let elastic = bind_isotropic_elastic_interface(
        &dry,
        normal.elastic().surface_a_state(),
        normal.elastic().surface_b_state(),
    )
    .unwrap();
    let query = common::query();
    let select = MaterialPropertySelection::SingleClaimOnly;
    let adhesion = resolve_interface_state_point(
        &card,
        &query,
        &[ScalarPropertyRequirement::try_new(
            ADHESION_ENERGY_PROPERTY,
            ADHESION_ENERGY_DIMS,
            ScalarAdmissibility::NonNegative,
        )
        .unwrap()],
        select.clone(),
    )
    .unwrap();
    let normal = bind_normal_interface_state(&elastic, &adhesion).unwrap();
    let config = original.contact().model().config();
    let contact = SourcePlasticContact::new_with_linear_solver(
        normal,
        original.contact().hardness_sources().clone(),
        config.grid,
        config.undeformed_gap_m.clone(),
        PlasticContactControls {
            maximum_active_set_iterations: config.maximum_active_set_iterations,
            complementarity_tolerance_m: config.complementarity_tolerance_m,
            boundary_clearance_cells: config.boundary_clearance_cells,
        },
        PlasticContactEnvelope {
            maximum_approach_m: 1.0e-4,
            maximum_plastic_depth_m: 1.0e-4,
            half_space_depth_m: [1.0e-2; 2],
            maximum_patch_to_depth: [0.1; 2],
        },
        original.contact().model().linear_solver(),
        cx,
    )
    .unwrap();
    let keys = original.coefficient_keys().clone();
    let coefficients = resolve_interface_state_point(
        &card,
        &query,
        &keys.each_ref().map(|key| {
            ScalarPropertyRequirement::try_with_key(key, ScalarAdmissibility::NonNegative).unwrap()
        }),
        select.clone(),
    )
    .unwrap();
    let wear = SourcePlasticArchardWear::new(
        contact,
        coefficients,
        keys,
        original.hardness_sources().clone(),
        maximum_recession_m,
    )
    .unwrap();
    let friction_query = if let Some([lo, hi]) = pressure_domain_pa {
        query
            .with_quantity("p", pressure_quantity, 0.5 * lo + 0.5 * hi)
            .unwrap()
    } else {
        query
    };
    let requirement = |name| {
        ScalarPropertyRequirement::try_new(name, Dims::NONE, ScalarAdmissibility::NonNegative)
            .unwrap()
    };
    let friction = SourceThermalCoulomb::resolve(
        dry,
        friction_query,
        requirement("mu_static"),
        requirement("mu_kinetic"),
        select,
        ThermalFrictionAxes {
            temperature: "T".into(),
            nominal_pressure: pressure_domain_pa.map(|_| "p".into()),
        },
    )
    .unwrap();
    SourcePlasticSlidingWear::new(wear, friction).unwrap()
}

#[test]
fn capped_pressure_and_actual_slip_drive_traction_wear_and_neighbor_indentation() {
    with_cx(|cx, _| {
        let law = source([1.0e6, 2.0e6], [1.0e-4; 2], true, None, cx);
        let loaded = law
            .wear_source()
            .advance_normal(&law.wear_source().initial_state(), 5.0e-5, cx)
            .unwrap();
        let before = law.state_from_wear(loaded.candidate).unwrap();
        let saved = before.clone();
        let mut input = motion();
        input[12] = PlasticSlidingInput {
            released_increment_m: [1.0e-5, 0.0],
            compliance_m_per_n: 0.01,
        };
        input[13] = PlasticSlidingInput {
            released_increment_m: [2.0e-4, 0.0],
            compliance_m_per_n: 0.01,
        };
        let next = law
            .advance(&before, &input, 0.1, friction_controls(), cx)
            .unwrap();
        let stuck = &next.cells[12].response;
        let sliding = &next.cells[13].response;
        assert_eq!(stuck.regime, FrictionRegime::Sticking);
        assert_eq!(stuck.slip_increment_m, [0.0; 2]);
        assert_eq!(next.recession_increment_m[0][12], 0.0);
        assert_eq!(next.recession_increment_m[1][12], 0.0);
        close(stuck.traction_n[0], 1.0e-3);
        let pressure = 1.0e6;
        let force = pressure * 1.0e-8;
        let traction = 0.4 * force;
        let distance = 2.0e-4 - 0.01 * traction;
        close(next.cells[13].pressure_pa, pressure);
        assert_eq!(sliding.regime, FrictionRegime::Sliding);
        close(sliding.traction_n[0], traction);
        close(sliding.slip_increment_m[0], distance);
        close(next.frictional_work_j, traction * distance);
        close(next.pressure_weighted_sliding_j, force * distance);
        for side in 0..2 {
            let depth = [0.01, 0.03][side] * pressure * distance / [1.0e6, 2.0e6][side];
            close(next.recession_increment_m[side][13], depth);
            close(next.wear[side].increment_m3, depth * 1.0e-8);
            close(
                next.candidate.wear_state().volumes_m3()[side],
                depth * 1.0e-8,
            );
        }
        assert_eq!(next.friction_sources.len(), 1);
        assert_eq!(next.cells[12].source_index, next.cells[13].source_index);
        assert_eq!(next.source_identity, law.identity());
        assert!(next.normal.response.plastic_increment_m[12] > 0.0);
        assert!(
            next.candidate
                .wear_state()
                .contact_state()
                .contact_state()
                .pressure_pa()[13]
                < pressure
        );
        assert_eq!(before, saved);
    });
}

#[test]
fn original_hardness_changes_traction_and_sliding_without_a_yield_conversion() {
    with_cx(|cx, _| {
        let mut outcomes = Vec::new();
        for hardness in [0.5e6, 1.0e6] {
            let law = source([hardness, 2.0e6], [1.0e-4; 2], false, None, cx);
            let loaded = law
                .wear_source()
                .advance_normal(&law.wear_source().initial_state(), 5.0e-5, cx)
                .unwrap();
            let before = law.state_from_wear(loaded.candidate).unwrap();
            let mut input = motion();
            input[12] = PlasticSlidingInput {
                released_increment_m: [2.0e-4, 0.0],
                compliance_m_per_n: 0.01,
            };
            let next = law
                .advance(&before, &input, 0.1, friction_controls(), cx)
                .unwrap();
            close(next.cells[12].pressure_pa, hardness);
            close(
                next.cells[12].response.traction_n[0],
                0.4 * hardness * 1.0e-8,
            );
            outcomes.push((
                next.cells[12].response.traction_n[0],
                next.cells[12].response.slip_increment_m[0],
            ));
        }
        assert!(outcomes[0].0 < outcomes[1].0);
        assert!(outcomes[0].1 > outcomes[1].1);
    });
}

#[test]
fn all_normal_boundaries_use_the_same_sliding_and_wear_step() {
    with_cx(|cx, _| {
        let law = source([1.0e6, 2.0e6], [1.0e-4; 2], false, None, cx);
        let force = 0.005;
        let loaded = law
            .wear_source()
            .advance_normal_at_load(
                &law.wear_source().initial_state(),
                force,
                normal_controls(),
                cx,
            )
            .unwrap();
        let before = law.state_from_wear(loaded.candidate).unwrap();
        let approach = before
            .wear_state()
            .contact_state()
            .contact_state()
            .approach_m();
        let mut input = motion();
        input[12].released_increment_m = [1.0e-4, 0.0];
        let held = law
            .advance(&before, &input, 0.1, friction_controls(), cx)
            .unwrap();
        let loaded = law
            .advance_at_load(
                &before,
                &input,
                0.1,
                friction_controls(),
                force,
                normal_controls(),
                cx,
            )
            .unwrap();
        let device = FiniteGapPlasticDrive {
            stiffness_n_per_m: 1000.0,
            position_m: approach + force / 1000.0,
        };
        let driven = law
            .advance_with_drive(
                &before,
                &input,
                0.1,
                friction_controls(),
                device,
                normal_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(held.cells, loaded.cells);
        assert_eq!(held.cells, driven.cells);
        assert_eq!(held.recession_increment_m, loaded.recession_increment_m);
        assert_eq!(held.recession_increment_m, driven.recession_increment_m);
        assert_eq!(
            held.normal.kind,
            FiniteGapPlasticContinuationKind::SurfaceRecession
        );
        assert_eq!(
            loaded.normal.kind,
            FiniteGapPlasticContinuationKind::SurfaceRecessionAtLoad
        );
        assert_eq!(
            driven.normal.kind,
            FiniteGapPlasticContinuationKind::SurfaceRecessionWithDrive
        );
        close(held.normal.approach_m, approach);
        assert!(held.normal.response.normal.normal_force_n < force);
        close(loaded.normal.response.normal.normal_force_n, force);
        assert!(loaded.normal.approach_m > approach);
        assert!(driven.normal.approach_m > approach);
        assert!(driven.normal.response.normal.normal_force_n < force);
        close(
            driven.boundary.drive_force_n,
            driven.normal.response.normal.normal_force_n,
        );
        assert!(driven.normal.normal_work_j > 0.0);
    });
}

#[test]
fn reversals_remove_material_using_each_actual_path_and_retain_native_shear_history() {
    with_cx(|cx, _| {
        let law = source([1.0e6, 2.0e6], [1.0e-4; 2], false, None, cx);
        let loaded = law
            .advance_normal(&law.initial_state(), 5.0e-5, cx)
            .unwrap();
        let before = loaded.candidate;
        let mut input = motion();
        input[12].released_increment_m = [1.0e-4, 0.0];
        let forward = law
            .advance(&before, &input, 0.1, friction_controls(), cx)
            .unwrap();
        let new_force = 0.003;
        let reloaded = law
            .advance_normal_at_load(&forward.candidate, new_force, normal_controls(), cx)
            .unwrap();
        close(reloaded.receipt.normal.normal_force_n, new_force);
        assert_eq!(
            reloaded.candidate.tangential_states(),
            forward.candidate.tangential_states()
        );
        assert_eq!(
            reloaded.candidate.wear_state().recession_m(),
            forward.candidate.wear_state().recession_m()
        );
        assert_eq!(
            reloaded.candidate.wear_state().volumes_m3(),
            forward.candidate.wear_state().volumes_m3()
        );
        let device = FiniteGapPlasticDrive {
            stiffness_n_per_m: 1000.0,
            position_m: reloaded
                .candidate
                .wear_state()
                .contact_state()
                .contact_state()
                .approach_m()
                + new_force / 1000.0,
        };
        let held = law
            .advance_normal_with_drive(&reloaded.candidate, device, normal_controls(), cx)
            .unwrap();
        assert_eq!(held.candidate, reloaded.candidate);
        input[12].released_increment_m[0] = -1.0e-4;
        let reverse = law
            .advance(&held.candidate, &input, 0.1, friction_controls(), cx)
            .unwrap();
        assert_eq!(reverse.candidate.tangential_states()[12].slip_m(), [0.0; 2]);
        assert!(reverse.cells[12].response.traction_n[0] < 0.0);
        assert!(forward.frictional_work_j > 0.0 && reverse.frictional_work_j > 0.0);
        for side in 0..2 {
            assert!(forward.wear[side].increment_m3 > 0.0 && reverse.wear[side].increment_m3 > 0.0);
            close(
                reverse.candidate.wear_state().volumes_m3()[side],
                forward.wear[side].increment_m3 + reverse.wear[side].increment_m3,
            );
        }
    });
}

#[test]
fn unloaded_cells_skip_positive_pressure_sources_but_validate_native_controls() {
    with_cx(|cx, _| {
        let law = source([1.0e6, 2.0e6], [1.0e-4; 2], true, Some([1.0e5, 2.0e6]), cx);
        let before = law.initial_state();
        let mut input = motion();
        input[0].released_increment_m = [0.001, -0.002];
        let unloaded = law
            .advance(&before, &input, 0.1, friction_controls(), cx)
            .unwrap();
        assert!(unloaded.friction_sources.is_empty());
        assert_eq!(unloaded.frictional_work_j, 0.0);
        assert_eq!(unloaded.candidate.wear_state().volumes_m3(), [0.0; 2]);
        assert!(unloaded.cells.iter().all(|cell| cell.source_index.is_none()
            && cell.response.traction_n == [0.0; 2]
            && cell.response.dissipated_j == 0.0));
        assert!(
            law.advance(&before, &input, 0.0, friction_controls(), cx)
                .is_err()
        );
        assert!(
            law.advance(
                &before,
                &input,
                0.1,
                TangentialSolveConfig {
                    max_iterations: 0,
                    ..friction_controls()
                },
                cx
            )
            .is_err()
        );
        input[24].compliance_m_per_n = -1.0;
        assert!(matches!(
            law.advance(&before, &input, 0.1, friction_controls(), cx),
            Err(SourcePlasticSlidingError::Cell { index: 24, .. })
        ));
        let loaded = law
            .wear_source()
            .advance_normal(before.wear_state(), 5.0e-5, cx)
            .unwrap();
        let loaded = law.state_from_wear(loaded.candidate).unwrap();
        let mut input = motion();
        input[12].released_increment_m = [1.0e-5, 0.0];
        input[13].released_increment_m = [0.0, 1.0e-5];
        let next = law
            .advance(&loaded, &input, 0.1, friction_controls(), cx)
            .unwrap();
        assert_eq!(next.friction_sources.len(), 1);
        let selected = &next.friction_sources[0];
        assert_eq!(selected.state_coordinate("p"), Some(1.0e6));
        let mu = selected.property("mu_kinetic").unwrap().value_si();
        close(mu, 0.4 * (1.0 - 0.5 * (1.0e6 - 1.0e5) / (2.0e6 - 1.0e5)));
        close(next.cells[12].response.traction_n[0], mu * 1.0e6 * 1.0e-8);
        assert_eq!(
            selected.card_identity(),
            law.friction_source().initial_source().card_identity()
        );
        assert_eq!(
            selected
                .property("mu_kinetic")
                .unwrap()
                .answer()
                .receipt
                .selected,
            law.friction_source()
                .initial_source()
                .property("mu_kinetic")
                .unwrap()
                .answer()
                .receipt
                .selected
        );
    });
}

#[test]
fn source_mismatch_late_cell_wear_and_normal_failures_leave_every_history_unchanged() {
    with_cx(|cx, gate| {
        let law = source([1.0e6, 2.0e6], [1.0e-8; 2], true, None, cx);
        let foreign = source([1.0e6, 2.0e6], [1.0e-8; 2], true, Some([1.0e5, 2.0e6]), cx);
        assert!(
            SourcePlasticSlidingWear::new(
                law.wear_source().clone(),
                foreign.friction_source().clone()
            )
            .is_err()
        );
        let loaded = law
            .wear_source()
            .advance_normal(&law.wear_source().initial_state(), 5.0e-5, cx)
            .unwrap();
        let before = law.state_from_wear(loaded.candidate).unwrap();
        let saved = before.clone();
        let mut input = motion();
        input[12].released_increment_m = [1.0e-4, 0.0];
        input[13].released_increment_m = [f64::NAN, 0.0];
        assert!(matches!(
            law.advance(&before, &input, 0.1, friction_controls(), cx),
            Err(SourcePlasticSlidingError::Cell {
                index: 13,
                error: SourceThermalCoulombError::Tribo(_)
            })
        ));
        assert_eq!(before, saved);
        input[13].released_increment_m = [0.0; 2];
        assert!(matches!(
            law.advance(&before, &input, 0.1, friction_controls(), cx),
            Err(SourcePlasticSlidingError::Wear(
                SourcePlasticWearError::Contact(_)
            ))
        ));
        assert_eq!(before, saved);
        input[12].released_increment_m = [1.0e-8, 0.0];
        assert!(
            law.advance_at_load(
                &before,
                &input,
                0.1,
                friction_controls(),
                0.001,
                normal_controls(),
                cx
            )
            .is_err()
        );
        assert_eq!(before, saved);
        assert!(
            law.advance(&before, &input, 0.1, friction_controls(), cx)
                .is_ok()
        );
        gate.request();
        assert!(
            law.advance(&before, &input, 0.1, friction_controls(), cx)
                .is_err()
        );
        assert_eq!(before, saved);
    });
}
