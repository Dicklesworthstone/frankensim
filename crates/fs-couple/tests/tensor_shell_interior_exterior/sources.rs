//! Original volume/P2 injection drives the material shell and exterior reaction.

use super::*;
use fs_feec::acoustics::helmholtz::heterogeneous::HeterogeneousError;

// Independent simplex moments for normalized uniform volume injection and
// g=4*l0*l1-4*l1*l2 on face [0,1,2]. The latter has ZERO integrated velocity,
// but nonzero vertex and edge loads. Basis order is original vertex + edge.
fn exact_load(volume_flow: C64, boundary_velocity: C64, linear_velocity: C64) -> Vec<C64> {
    let area = 0.4 * 0.3 / 2.0;
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let mut load = vec![volume_flow.scale(0.25); 4];
    load[0] = load[0] + boundary_velocity.scale(area / 15.0);
    load[2] = load[2] - boundary_velocity.scale(area / 15.0);
    for value in &mut load[..3] {
        *value = *value + linear_velocity.scale(area / 3.0);
    }
    for edge in &mesh.edges {
        let signed = match *edge {
            [0, 1] => area / 45.0,
            [1, 2] => -area / 45.0,
            _ => 0.0,
        };
        let linear = if edge.iter().all(|&v| v < 3) {
            area / 12.0
        } else {
            0.0
        };
        load.push(
            volume_flow.scale(1.0 / 20.0)
                + boundary_velocity.scale(signed)
                + linear_velocity.scale(linear),
        );
    }
    load
}

#[test]
fn distributed_injection_drives_original_moving_shell_and_full_mutual_radiation() {
    fixture(|cx, _, shell, interior, geometry| {
        let owner = TensorShellInterior::new(
            interior,
            shell,
            &interfaces(),
            &wall(),
            OMEGA,
            interior_controls(),
            cx,
        )
        .unwrap();
        let boundary = geometry.prepare(120.0, boundary_budget(), cx).unwrap();
        let prepared = owner
            .prepare_with_exterior(&boundary, exterior_controls(), cx)
            .unwrap();
        let volume = interior
            .volume_velocity_source(
                4,
                &mut |_, _| Ok::<_, HeterogeneousError>(100.0),
                &mut || false,
            )
            .unwrap();
        let signed = interior
            .quadratic_velocity_source(&[([0, 1, 2], [0.0, 0.0, 0.0, 1.0, -1.0, 0.0])], &mut || {
                false
            })
            .unwrap();
        assert!((volume.integrated_shape() - 1.0).abs() < 3e-14);
        assert!(signed.integrated_shape_m2().abs() < 3e-16);
        assert!(signed.has_nonzero_projection());
        let q = C64::new(1e-8, -4e-9);
        let u = C64::new(-3e-7, 2e-7);
        let linear = C64::new(4e-8, 1e-8);
        let response = prepared
            .solve_with_sources(
                &[C64::ZERO; 6],
                &[([0, 1, 2], [linear; 3])],
                &[(&signed, u)],
                &[(&volume, q)],
                cx,
            )
            .unwrap();
        let field = &response.coupled.field;
        let load = exact_load(q, u, linear);
        let wall_mass = interior.geometry().boundary_mass(&wall()[0].faces).unwrap();
        for (row, &expected) in load.iter().enumerate() {
            let obtained = product(interior.stiffness(), &field.pressure_potential_pa_s, row)
                - product(interior.mass(), &field.pressure_potential_pa_s, row)
                    .scale(OMEGA * OMEGA)
                + product(&wall_mass, &field.pressure_pa, row).scale(WALL_Y)
                + field.displaced_flow_m3_s[row];
            assert!((obtained - expected).abs() < 3e-8 * q.abs());
        }
        let source_work = work(&load, &field.pressure_pa);
        assert!(
            (source_work - field.average_acoustic_input_power_w).abs()
                < 2e-10 * source_work.abs().max(1e-25)
        );
        assert_eq!(field.average_structural_input_power_w, 0.0);
        assert!(response.radiation.outgoing_power_w > 0.0);
        assert!(field.average_structural_dissipation_w > 0.0);
        assert!(field.pressure_pa.iter().skip(4).any(|p| p.abs() > 1e-8));
        assert!(field.structural_velocity.iter().any(|v| v.abs() > 1e-14));
        let native_work = field.average_structural_dissipation_w
            + field.average_wall_dissipation_w
            + response.radiation.outgoing_power_w;
        assert!((source_work - native_work).abs() < 3e-8 * source_work.abs());

        // All pressure coordinates and all mutual exterior loads must give
        // the same answer as the independently derived original-coordinate RHS.
        let explicit = prepared.solve(&[C64::ZERO; 6], &load, cx).unwrap();
        for (&actual, &expected) in field
            .pressure_pa
            .iter()
            .zip(&explicit.coupled.field.pressure_pa)
        {
            near(actual, expected);
        }
        for (&actual, &expected) in field
            .structural_velocity
            .iter()
            .zip(&explicit.coupled.field.structural_velocity)
        {
            near(actual, expected);
        }
        for (&actual, &expected) in response
            .exterior
            .surface_pressure_pa()
            .iter()
            .zip(explicit.exterior.surface_pressure_pa())
        {
            near(actual, expected);
        }
    });
}

#[test]
fn physical_source_bridge_preserves_coherent_shell_work_and_rejects_foreign_owners() {
    fixture(|cx, gate, shell, interior, geometry| {
        let owner = TensorShellInterior::new(
            interior,
            shell,
            &interfaces(),
            &wall(),
            OMEGA,
            interior_controls(),
            cx,
        )
        .unwrap();
        let volume = interior
            .volume_velocity_source(
                4,
                &mut |_, _| Ok::<_, HeterogeneousError>(100.0),
                &mut || false,
            )
            .unwrap();
        let boundary = geometry.prepare(120.0, boundary_budget(), cx).unwrap();
        let prepared = owner
            .prepare_with_exterior(&boundary, exterior_controls(), cx)
            .unwrap();
        let q = C64::new(1e-8, 4e-9);
        let driven = prepared
            .solve_with_sources(&force(), &[], &[], &[(&volume, q)], cx)
            .unwrap();
        let mechanical = prepared
            .solve_with_sources(&force(), &[], &[], &[], cx)
            .unwrap();
        let injection = prepared
            .solve_with_sources(&[C64::ZERO; 6], &[], &[], &[(&volume, q)], cx)
            .unwrap();
        let scale = driven
            .coupled
            .field
            .pressure_pa
            .iter()
            .map(|p| p.abs())
            .fold(1e-25_f64, f64::max);
        for ((&all, &a), &b) in driven
            .coupled
            .field
            .pressure_pa
            .iter()
            .zip(&mechanical.coupled.field.pressure_pa)
            .zip(&injection.coupled.field.pressure_pa)
        {
            assert!((all - a - b).abs() < 2e-8 * scale);
        }
        let expected_work = work(
            &exact_load(q, C64::ZERO, C64::ZERO),
            &driven.coupled.field.pressure_pa,
        );
        assert!(
            (expected_work - driven.coupled.field.average_acoustic_input_power_w).abs()
                < 2e-10 * expected_work.abs().max(1e-25)
        );

        // The material-only path uses exactly the same native source loading;
        // it rechecks full-cycle material strains and its own reciprocal work.
        let interior_only = owner.prepare(cx).unwrap();
        let inner = interior_only
            .solve_with_sources(&force(), &[], &[], &[(&volume, q)], cx)
            .unwrap();
        let original = interior_only
            .solve(&force(), &exact_load(q, C64::ZERO, C64::ZERO), cx)
            .unwrap();
        for (&actual, &expected) in inner
            .field
            .pressure_pa
            .iter()
            .zip(&original.field.pressure_pa)
        {
            near(actual, expected);
        }
        let other = HeterogeneousHelmholtzFe::new(
            interior.geometry(),
            4,
            &mut |_, _| {
                Ok::<_, AcousticFeError>(Medium {
                    density_kg_m3: 12.0,
                    sound_speed_m_s: 200.0,
                })
            },
            &mut || false,
        )
        .unwrap();
        let foreign = other
            .volume_velocity_source(
                4,
                &mut |_, _| Ok::<_, HeterogeneousError>(100.0),
                &mut || false,
            )
            .unwrap();
        assert!(
            prepared
                .solve_with_sources(&[C64::ZERO; 6], &[], &[], &[(&foreign, C64::ZERO)], cx)
                .is_err()
        );
        assert!(
            interior_only
                .solve_with_sources(&[C64::ZERO; 6], &[], &[], &[(&foreign, C64::ZERO)], cx)
                .is_err()
        );
        gate.request();
        assert!(
            prepared
                .solve_with_sources(&[C64::ZERO; 6], &[], &[], &[(&volume, q)], cx)
                .is_err()
        );
    });
}
