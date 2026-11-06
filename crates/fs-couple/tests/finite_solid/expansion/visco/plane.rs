//! Source-bound thin thermal events keep quadrature memory and accept atomically.
mod dynamics;

use super::*;
use fs_couple::finite_solid::thermal::visco::plane::{
    ExpansionViscoPlaneError, ExpansionViscoPlaneField, ExpansionViscoPlaneSolid,
};
use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_solid::hyper2d::visco::{
    FiniteViscoPlaneError, FiniteViscoPlaneMode, FiniteViscoStepControls,
};
use fs_solid::{Mesh2, NewtonSettings, Patch};

fn settings() -> NewtonSettings {
    NewtonSettings {
        tol: 2e-9,
        max_iters: 40,
        load_steps: 1,
        max_backtracks: 30,
    }
}

fn allowances() -> FiniteViscoStepControls {
    FiniteViscoStepControls {
        equilibrium_work_tolerance_j: 1e-8,
        maximum_endpoint_work_defect_j: 1e-3,
    }
}

fn held(x: f64, y: f64) -> [f64; 2] {
    [0.08 * x + 0.025 * x * y, -0.04 * y + 0.015 * x * y]
}

fn walls<'a>(
    value: &'a dyn Fn(f64, f64) -> [f64; 2],
) -> [(Patch, &'a dyn Fn(f64, f64) -> [f64; 2]); 4] {
    [
        (Patch::Left, value),
        (Patch::Right, value),
        (Patch::Bottom, value),
        (Patch::Top, value),
    ]
}

#[test]
fn thermal_plane_field_retains_nonuniform_memory_and_separates_both_reference_inputs() {
    with_cx(|_, cx| {
        // These source and mechanical calibrations are the existing authored
        // test fixtures; the adapter does not infer Prony data from expansion.
        let source = card(STRETCHES, "T", false);
        let datum = endpoint(&source, Some(shift()));
        let field = ExpansionViscoPlaneField::uniform(datum.clone());
        let element_field =
            ExpansionViscoPlaneField::per_element(vec![datum.clone(), datum.clone()], cx).unwrap();
        let temperatures = [300.2, 299.8];
        let target = element_field
            .at_temperatures(&[&source, &source], &temperatures, cx)
            .unwrap();
        for (element, temperature) in temperatures.into_iter().enumerate() {
            let expected = datum.at_temperature(&source, temperature).unwrap();
            assert_eq!(target.endpoint_at(element).unwrap(), &expected);
            near(field.endpoints()[0].temperature_k(), 300.0, 0.0);
        }
        let mesh = Mesh2::quads(1.0, 0.7, 2, 1);
        let initial: Vec<_> = mesh.nodes.iter().map(|[x, y]| held(*x, *y)).collect();
        let boundary = walls(&held);
        for mode in [
            FiniteViscoPlaneMode::FiniteBulk,
            FiniteViscoPlaneMode::Incompressible,
        ] {
            let virgin = ExpansionViscoPlaneSolid::try_new_from(
                &mesh,
                &field,
                mode,
                0.12,
                &initial,
                PlaneStressSettings::default(),
                cx,
            )
            .unwrap();
            let accepted = virgin
                .advance(0.1, &boundary, &[], settings(), allowances(), cx)
                .unwrap()
                .solid;
            let before_full = accepted
                .native()
                .point_states()
                .map(|states| states.to_vec());
            let before_remainder = accepted
                .native()
                .point_remainder_states()
                .map(|states| states.to_vec());
            let before_u = accepted.native().displacement_m().to_vec();
            assert_eq!(accepted.native().material_points().len(), 8);
            let changed = accepted
                .change_reference(&target, allowances(), cx)
                .unwrap();
            assert_eq!(changed.solid.native().time_s(), accepted.native().time_s());
            assert_eq!(changed.solid.native().displacement_m(), before_u);
            assert_eq!(
                changed.solid.native().material_points(),
                accepted.native().material_points()
            );
            if let Some(before) = &before_full {
                let after = changed.solid.native().point_states().unwrap();
                assert!(before.iter().any(|state| state.dissipated > 0.0));
                assert_ne!(
                    before[0].inverse_viscous_cauchy_green,
                    before[1].inverse_viscous_cauchy_green
                );
                for (old, new) in before.iter().zip(after) {
                    assert_eq!(
                        old.inverse_viscous_cauchy_green,
                        new.inverse_viscous_cauchy_green
                    );
                    assert_eq!(old.dissipated, new.dissipated);
                }
            }
            if let Some(before) = &before_remainder {
                let after = changed.solid.native().point_remainder_states().unwrap();
                assert!(before.iter().any(|state| state.dissipated > 0.0));
                assert_ne!(
                    before[0].inverse_viscous_cauchy_green,
                    before[1].inverse_viscous_cauchy_green
                );
                for (old, new) in before.iter().zip(after) {
                    assert_eq!(
                        old.inverse_viscous_cauchy_green,
                        new.inverse_viscous_cauchy_green
                    );
                    assert_eq!(old.dissipated, new.dissipated);
                }
            }
            let reference = &changed.energy;
            near(
                reference.stored_after_j - reference.stored_before_j,
                reference.reference_work_j
                    + reference.constraint_reference_work_j
                    + reference.transverse_endpoint_work_j
                    - reference.endpoint_work_defect_j
                    - reference.material_balance_error_j,
                1e-11,
            );
            if mode == FiniteViscoPlaneMode::Incompressible {
                // Changing the source's natural volume produces physical
                // pressure actuation, not a spurious dissipation increment.
                assert!(reference.constraint_reference_work_j.abs() > 1e-9);
            } else {
                near(reference.constraint_reference_work_j, 0.0, 0.0);
            }
            let sequential = changed
                .solid
                .advance(0.03, &boundary, &[], settings(), allowances(), cx)
                .unwrap();
            let combined = accepted
                .change_reference_and_advance(
                    &target,
                    0.03,
                    &boundary,
                    &[],
                    settings(),
                    allowances(),
                    allowances(),
                    cx,
                )
                .unwrap();
            assert_eq!(combined.reference, changed.energy);
            assert_eq!(combined.mechanical, sequential.energy);
            assert_eq!(combined.reactions_n, sequential.reactions_n);
            assert_eq!(
                combined.solid.native().point_states(),
                sequential.solid.native().point_states()
            );
            assert_eq!(
                combined.solid.native().point_remainder_states(),
                sequential.solid.native().point_remainder_states()
            );
            near(
                combined.reference_input_j(),
                combined.reference.reference_work_j
                    + combined.reference.constraint_reference_work_j,
                0.0,
            );
            assert!(std::ptr::eq(combined.solid.source_field(), &target));
            assert!(std::ptr::eq(
                combined.solid.native().reference_mesh(),
                &mesh
            ));
            near(combined.solid.native().reference_thickness_m(), 0.12, 0.0);
            near(combined.solid.native().time_s(), 0.13, 1e-15);
            assert!(std::ptr::eq(accepted.source_field(), &field));
            assert_eq!(accepted.native().point_states(), before_full.as_deref());
            assert_eq!(
                accepted.native().point_remainder_states(),
                before_remainder.as_deref()
            );
            assert_eq!(accepted.native().displacement_m(), before_u);
        }
    });
}

#[test]
fn thermal_plane_source_admission_and_late_mechanical_failure_preserve_accepted_owner() {
    with_cx(|_, cx| {
        let source = card(STRETCHES, "T", false);
        let datum = endpoint(&source, Some(shift()));
        let field = ExpansionViscoPlaneField::uniform(datum.clone());
        let target = field.at_temperatures(&[&source], &[300.1], cx).unwrap();
        assert!(field.at_temperatures(&[&source], &[], cx).is_err());
        assert!(
            field
                .at_temperatures(&[&source, &source], &[300.0], cx)
                .is_err()
        );
        assert!(field.at_temperatures(&[&source], &[320.0], cx).is_err());
        assert!(ExpansionViscoPlaneField::per_element(vec![], cx).is_err());

        let mesh = Mesh2::quads(1.0, 0.7, 2, 1);
        let initial: Vec<_> = mesh.nodes.iter().map(|[x, y]| held(*x, *y)).collect();
        let different_source = card([1.12, 0.96, 1.04], "T", false);
        let incompatible = ExpansionViscoPlaneField::per_element(
            vec![datum.clone(), endpoint(&different_source, Some(shift()))],
            cx,
        )
        .unwrap();
        let wrong_count = ExpansionViscoPlaneField::per_element(vec![datum.clone()], cx).unwrap();
        let wrong_calibration =
            ExpansionViscoPlaneField::per_element(vec![datum.clone(), endpoint(&source, None)], cx)
                .unwrap();
        let next_held = |x: f64, y: f64| {
            let mut u = held(x, y);
            u[0] += 0.03 * x;
            u
        };
        let boundary = walls(&next_held);
        for mode in [
            FiniteViscoPlaneMode::FiniteBulk,
            FiniteViscoPlaneMode::Incompressible,
        ] {
            let accepted = ExpansionViscoPlaneSolid::try_new_from(
                &mesh,
                &field,
                mode,
                0.12,
                &initial,
                PlaneStressSettings::default(),
                cx,
            )
            .unwrap();
            let before_full = accepted
                .native()
                .point_states()
                .map(|states| states.to_vec());
            let before_remainder = accepted
                .native()
                .point_remainder_states()
                .map(|states| states.to_vec());
            for invalid in [&incompatible, &wrong_calibration] {
                assert!(matches!(
                    accepted.change_reference(invalid, allowances(), cx),
                    Err(ExpansionViscoPlaneError::Source { element: 1, .. })
                ));
            }
            assert!(matches!(
                accepted.change_reference(&wrong_count, allowances(), cx),
                Err(ExpansionViscoPlaneError::Invalid(_))
            ));
            // The source event is admitted on its own. A deliberately zero
            // endpoint-work allowance then rejects the actual loading step,
            // after source and native local preparation have both succeeded.
            accepted
                .change_reference(&target, allowances(), cx)
                .unwrap();
            let strict = FiniteViscoStepControls {
                maximum_endpoint_work_defect_j: 0.0,
                ..allowances()
            };
            assert!(matches!(
                accepted.change_reference_and_advance(
                    &target,
                    0.03,
                    &boundary,
                    &[],
                    settings(),
                    allowances(),
                    strict,
                    cx,
                ),
                Err(ExpansionViscoPlaneError::Native(
                    FiniteViscoPlaneError::EnergyGate { .. }
                ))
            ));
            assert!(std::ptr::eq(accepted.source_field(), &field));
            assert_eq!(accepted.native().point_states(), before_full.as_deref());
            assert_eq!(
                accepted.native().point_remainder_states(),
                before_remainder.as_deref()
            );
            assert_eq!(accepted.native().displacement_m(), initial);
            near(accepted.native().time_s(), 0.0, 0.0);
            near(field.endpoints()[0].temperature_k(), 300.0, 0.0);
        }
    });
}
