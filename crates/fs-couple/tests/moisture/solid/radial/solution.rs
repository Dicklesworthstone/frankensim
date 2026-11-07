//! Solution activity drives actual native geometry and wet mechanical inertia.
use super::*;
use crate::{atmosphere as air, solution as bath};

#[test]
fn finite_solution_radial_uptake_drives_native_shape_and_wet_inertia() {
    with_cx(|cx| {
        let volume = 1e-4;
        let reservoir = bath::reservoir(volume, vec![0.3]);
        let saved = reservoir.solution().clone();
        let original = radial_specimen(card([0.01; 2], 0.3), vec![0.0, 0.04], vec![0.1], cx);
        let binding = radial_binding(&original, vec![0, 1, 2, 4, 5, 11], cx);
        let faces = [None, Some(air::sorption(0.25))];
        for second_order in [false, true] {
            let (next_reservoir, next, solid, transfer) = if second_order {
                let (candidate, solid) = binding
                    .advance_finite_solution_sdirk2(
                        &reservoir,
                        &original,
                        0.025,
                        &faces,
                        config(),
                        nonlinear(),
                        bath::controls(volume),
                        cx,
                    )
                    .unwrap();
                let (r, m, step) = candidate.into_parts();
                (r, m, solid, step.water_to_solid_kg)
            } else {
                let (candidate, solid) = binding
                    .advance_finite_solution(
                        &reservoir,
                        &original,
                        0.025,
                        &faces,
                        config(),
                        nonlinear(),
                        bath::controls(volume),
                        cx,
                    )
                    .unwrap();
                let (r, m, step) = candidate.into_parts();
                (r, m, solid, step.water_to_solid_kg)
            };
            assert!(transfer > 0.0);
            assert!(
                next_reservoir.solution().water_activity().unwrap()
                    < reservoir.solution().water_activity().unwrap()
            );
            near(
                next_reservoir.solvent_mass_kg() + solid.wet_mass_kg(),
                reservoir.solvent_mass_kg() + original.wet_mass_kg(),
                3e-12,
            );
            near(
                solid.tangent().unwrap().pencil().total_mass_kg,
                next.wet_mass_kg(),
                1e-12,
            );
            let increment = next.moisture()[0] - original.moisture()[0];
            for (x, u) in binding
                .reference_nodes_m()
                .iter()
                .zip(solid.solution().displacement_m())
            {
                for k in 0..3 {
                    near(u[k], x[k] * [0.01, 0.04, 0.02][k] * increment, 1e-10);
                }
            }
            assert_eq!(solid.time_s(), next.time_s());
            assert_eq!(
                next_reservoir.solution().solute_moles(),
                saved.solute_moles()
            );
        }
        assert_eq!(reservoir.solution(), &saved);
        assert_eq!(original.moisture(), &[0.1]);
    });
}

#[test]
fn late_radial_solution_mechanical_refusal_preserves_both_original_inventories() {
    with_cx(|cx| {
        let volume = 1e-4;
        let reservoir = bath::reservoir(volume, vec![0.3]);
        let saved = reservoir.solution().clone();
        let limited = radial_specimen(card([0.01; 2], 0.12), vec![0.0, 0.04], vec![0.1], cx);
        let binding = radial_binding(&limited, vec![0, 1, 2, 4, 5, 11], cx);
        let faces = [None, Some(air::sorption(0.25))];
        let native = reservoir
            .advance_radial_sdirk2(
                &limited,
                0.025,
                &faces,
                config(),
                nonlinear(),
                bath::controls(volume),
                cx,
            )
            .unwrap();
        assert!(native.moisture().moisture()[0] > 0.12);
        assert!(
            binding
                .advance_finite_solution_sdirk2(
                    &reservoir,
                    &limited,
                    0.025,
                    &faces,
                    config(),
                    nonlinear(),
                    bath::controls(volume),
                    cx
                )
                .is_err()
        );
        assert!(
            binding
                .advance_finite_solution(
                    &reservoir,
                    &limited,
                    0.025,
                    &faces,
                    config(),
                    nonlinear(),
                    bath::controls(volume),
                    cx
                )
                .is_err()
        );
        assert_eq!(reservoir.solution(), &saved);
        assert_eq!(limited.moisture(), &[0.1]);
        assert_eq!(limited.time_s(), 0.0);
    });
}
