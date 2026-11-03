//! G1 independent moving-length/area flux equations and G3 transactional limits.
use super::*;
use fs_couple::moisture::layered::MoistureActivityTransfer;

fn root(mut lo: f64, mut hi: f64, residual: impl Fn(f64) -> f64) -> f64 {
    assert!(residual(lo) <= 0.0 && residual(hi) >= 0.0);
    for _ in 0..100 {
        let mid = lo + 0.5 * (hi - lo);
        if residual(mid) < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    lo + 0.5 * (hi - lo)
}

fn one_cell(beta: [f64; 3], cx: &Cx<'_>) -> MoistureLayeredSlab {
    MoistureLayeredSlab::new(
        vec![layer(material(500.0, beta, 0.3, 0.3), 0.1, vec![0.1], cx)],
        16,
        cx,
    )
    .unwrap()
}

#[test]
fn g1_actual_swelling_length_enters_the_implicit_water_balance() {
    with_cx(|cx| {
        let beta = 0.1;
        let stack = one_cell([beta, 0.0, 0.0], cx);
        let solid = bind_stack(&stack, free_supports(), cx);
        let faces = [MoistureLayeredFace::PrescribedRelativeHumidity(2.0 / 3.0); 2];
        let solved = solid
            .advance_spatial_diffusion(&stack, 0.25, &faces, config(), nonlinear(), 1e-9, cx)
            .unwrap();
        // One cell, two half-cell Dirichlet faces, D=.01, L=.1, dt=.25:
        // (M-.1) * (1 + beta*(M-.1))² = .2-M. The native source-integrated
        // solid supplies the stretch; this scalar polynomial is independent.
        let expected = root(0.1, 0.2, |m| {
            (m - 0.1) * (1.0 + beta * (m - 0.1)).powi(2) - (0.2 - m)
        });
        let m = solved.moisture.layers()[0].slab().moisture()[0];
        near(m, expected, 2e-11);
        assert!(m < 0.15 - 1e-5);
        near(
            solved.transport.axial_stretch[0],
            1.0 + beta * (m - 0.1),
            2e-11,
        );
        near(solved.transport.current_cross_section_m2, 0.001, 2e-12);
        let water = solved.transport.transport.boundary_water_kg;
        near(water[0] + water[1], 0.05 * (m - 0.1), 1e-12);
        near(solved.moisture.dry_mass_kg(), stack.dry_mass_kg(), 0.0);
        near(
            solved.solid.wet_mass_kg(),
            solved.moisture.wet_mass_kg(),
            1e-13,
        );
        assert!(solved.transport.native_solid_solves > 2);
        assert_eq!(
            solved.transport.native_solid_solves,
            solved.transport.transport.nonlinear.iterations + 1,
        );
        assert_eq!(
            solved.transport.geometry_identity,
            solved.solid.geometry().identity()
        );
        assert_eq!(stack.layers()[0].slab().moisture(), &[0.1]);
        assert_eq!(solved.moisture.time_s(), 0.25);
    });
}

#[test]
fn g1_finite_surface_exchange_uses_actual_area_and_conserved_dry_mass() {
    with_cx(|cx| {
        let beta = [0.1, 0.08, 0.04];
        let h = 0.4;
        let stack = one_cell(beta, cx)
            .with_transfer(
                [MoistureActivityTransfer::ConductanceKgM2S(h); 2],
                vec![],
                cx,
            )
            .unwrap();
        let solid = bind_stack(&stack, free_supports(), cx);
        let dt = 20.0;
        let solved = solid
            .advance_spatial_diffusion(
                &stack,
                dt,
                &[MoistureLayeredFace::PrescribedRelativeHumidity(2.0 / 3.0); 2],
                config(),
                nonlinear(),
                1e-9,
                cx,
            )
            .unwrap();
        let expected = root(0.1, 0.2, |m| {
            let stretch = beta.map(|b| 1.0 + b * (m - 0.1));
            let area = 0.001 * stretch[1] * stretch[2];
            let half = 2.0 * 0.05 * 0.01 * 0.3 / (0.1 * stretch[0]).powi(2);
            let g = 1.0 / (1.0 / half + 1.0 / (h * area));
            0.05 * (m - 0.1) - 2.0 * dt * g * (2.0 / 3.0 - m / 0.3)
        });
        let m = solved.moisture.layers()[0].slab().moisture()[0];
        near(m, expected, 3e-11);
        let stretch = beta.map(|b| 1.0 + b * (m - 0.1));
        near(
            solved.transport.current_cross_section_m2,
            0.001 * stretch[1] * stretch[2],
            1e-12,
        );
        for side in 0..2 {
            let a_surface = solved.transport.transport.surface_activities[side].unwrap();
            near(
                solved.transport.transport.boundary_water_kg[side],
                dt * h * solved.transport.current_cross_section_m2 * (2.0 / 3.0 - a_surface),
                1e-12,
            );
        }
        near(
            solved.moisture.water_mass_kg() - stack.water_mass_kg(),
            solved.transport.transport.boundary_water_kg.iter().sum(),
            1e-12,
        );
    });
}

#[test]
fn g1_unlike_stretched_layers_exchange_one_conserved_interface_inventory() {
    with_cx(|cx| {
        let stack = MoistureLayeredSlab::new(
            vec![
                layer(
                    material(500.0, [0.1, 0.0, 0.0], 0.3, 0.3),
                    0.1,
                    vec![0.2],
                    cx,
                ),
                layer(
                    material(800.0, [0.2, 0.0, 0.0], 0.3, 0.2),
                    0.2,
                    vec![0.1],
                    cx,
                ),
            ],
            16,
            cx,
        )
        .unwrap();
        let solid = bind_stack(&stack, free_supports(), cx);
        let dt = 0.2;
        let solved = solid
            .advance_spatial_diffusion(
                &stack,
                dt,
                &[MoistureLayeredFace::Sealed; 2],
                config(),
                nonlinear(),
                1e-9,
                cx,
            )
            .unwrap();
        let maximum_q = (0.2 / 0.3 - 0.1 / 0.2) / (1.0 / (0.05 * 0.3) + 1.0 / (0.16 * 0.2));
        let q = root(0.0, maximum_q, |q| {
            let m = [0.2 - q / 0.05, 0.1 + q / 0.16];
            let stretch = [1.0 + 0.1 * (m[0] - 0.1), 1.0 + 0.2 * (m[1] - 0.1)];
            let g = [
                2.0 * 0.05 * 0.01 * 0.3 / (0.1 * stretch[0]).powi(2),
                2.0 * 0.16 * 0.01 * 0.2 / (0.2 * stretch[1]).powi(2),
            ];
            q - dt * (m[0] / 0.3 - m[1] / 0.2) / (1.0 / g[0] + 1.0 / g[1])
        });
        near(
            solved.moisture.layers()[0].slab().moisture()[0],
            0.2 - q / 0.05,
            3e-11,
        );
        near(
            solved.moisture.layers()[1].slab().moisture()[0],
            0.1 + q / 0.16,
            3e-11,
        );
        near(
            solved.transport.transport.interfaces[0].water_left_to_right_kg,
            q,
            2e-12,
        );
        near(
            solved.moisture.water_mass_kg(),
            stack.water_mass_kg(),
            2e-12,
        );
        assert_eq!(solved.transport.transport.boundary_water_kg, [0.0; 2]);
        assert_ne!(
            solved.transport.axial_stretch[0],
            solved.transport.axial_stretch[1]
        );
    });
}

#[test]
fn g3_undeformed_limit_and_material_split_keep_the_original_transport() {
    with_cx(|cx| {
        let stack = one_cell([0.0; 3], cx);
        let solid = bind_stack(&stack, free_supports(), cx);
        let faces = [MoistureLayeredFace::PrescribedRelativeHumidity(2.0 / 3.0); 2];
        let original = stack
            .advance(0.25, &faces, config(), nonlinear(), cx)
            .unwrap();
        let spatial = solid
            .advance_spatial_diffusion(&stack, 0.25, &faces, config(), nonlinear(), 1e-9, cx)
            .unwrap();
        assert_eq!(original.0.activity(), spatial.moisture.activity());
        assert_eq!(
            original.1.boundary_water_kg,
            spatial.transport.transport.boundary_water_kg
        );
        assert_eq!(spatial.transport.axial_stretch, vec![1.0]);

        let card = material(500.0, [0.1, 0.0, 0.0], 0.3, 0.3);
        let whole =
            MoistureLayeredSlab::new(vec![layer(card.clone(), 0.2, vec![0.2, 0.1], cx)], 16, cx)
                .unwrap();
        let split = MoistureLayeredSlab::new(
            vec![
                layer(card.clone(), 0.1, vec![0.2], cx),
                layer(card, 0.1, vec![0.1], cx),
            ],
            16,
            cx,
        )
        .unwrap();
        let a = bind_stack(&whole, free_supports(), cx)
            .advance_spatial_diffusion(
                &whole,
                0.2,
                &[MoistureLayeredFace::Sealed; 2],
                config(),
                nonlinear(),
                1e-9,
                cx,
            )
            .unwrap();
        let b = bind_stack(&split, free_supports(), cx)
            .advance_spatial_diffusion(
                &split,
                0.2,
                &[MoistureLayeredFace::Sealed; 2],
                config(),
                nonlinear(),
                1e-9,
                cx,
            )
            .unwrap();
        for (&left, &right) in a.moisture.activity().iter().zip(b.moisture.activity()) {
            near(left, right, 2e-12);
        }
        assert_eq!(a.transport.transport.interfaces.len(), 0);
        assert_eq!(b.transport.transport.interfaces.len(), 1);
    });
}

#[test]
fn g3_work_source_and_incompatible_geometry_refusals_leave_original_state() {
    with_cx(|cx| {
        let stack = one_cell([0.1, 0.0, 0.0], cx);
        let solid = bind_stack(&stack, free_supports(), cx);
        let faces = [MoistureLayeredFace::PrescribedRelativeHumidity(2.0 / 3.0); 2];
        let mut bounded = nonlinear();
        bounded.maximum_iterations = 1;
        assert!(matches!(
            solid.advance_spatial_diffusion(&stack, 0.25, &faces, config(), bounded, 1e-9, cx,),
            Err(MoistureError::NonlinearNotConverged(_))
        ));
        assert_eq!(stack.layers()[0].slab().moisture(), &[0.1]);
        assert_eq!(stack.time_s(), 0.0);

        let restricted = MoistureLayeredSlab::new(
            vec![layer(
                material(500.0, [0.1, 0.0, 0.0], 0.12, 0.3),
                0.1,
                vec![0.1],
                cx,
            )],
            16,
            cx,
        )
        .unwrap();
        let binding = bind_stack(&restricted, free_supports(), cx);
        assert!(
            binding
                .advance_spatial_diffusion(
                    &restricted,
                    0.25,
                    &faces,
                    config(),
                    nonlinear(),
                    1e-9,
                    cx,
                )
                .is_err()
        );
        assert_eq!(restricted.layers()[0].slab().moisture(), &[0.1]);

        // A laterally clamped face with transverse free swelling bends/shears
        // the prism. Native mechanics accepts it; the 1-D spatial reduction
        // must refuse it rather than invent an axial-only metric.
        let warped = one_cell([0.01, 0.08, 0.04], cx);
        let clamped = bind_stack(&warped, (0..12).collect(), cx);
        assert!(
            clamped
                .advance(&warped, 0.25, &faces, config(), nonlinear(), cx)
                .is_ok()
        );
        assert!(matches!(clamped.advance_spatial_diffusion(
            &warped, 0.25, &faces, config(), nonlinear(), 1e-9, cx,
        ), Err(MoistureError::Invalid(message)) if message.contains("axis-aligned")));
        assert_eq!(warped.layers()[0].slab().moisture(), &[0.1]);
        assert_eq!(warped.time_s(), 0.0);
    });
}
