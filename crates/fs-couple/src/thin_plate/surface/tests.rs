use super::*;

fn close(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() < 3e-12 * scale,
        "{actual:e} != {expected:e} (scale {scale:e})"
    );
}

fn rectangle() -> PlateMesh {
    PlateMesh {
        nodes: vec![(0.0, 0.0), (0.4, 0.0), (0.4, 0.3), (0.0, 0.3)],
        tris: vec![[0, 1, 2], [0, 2, 3]],
    }
}

fn samples(
    mesh: &PlateMesh,
    levels: usize,
    modes: usize,
) -> (Vec<[f64; 3]>, Vec<SurfaceSample>, f64) {
    let (divisions, points, total) = surface_work(mesh.tris.len(), modes, levels).unwrap();
    let weights = triangle_weights(divisions, points);
    assert_eq!(weights.len(), points);
    let (surface, diameter) = sample_surface(mesh, &weights, divisions, modes, total).unwrap();
    assert_eq!(surface.len(), total);
    (weights, surface, diameter)
}

#[test]
fn three_rungs_preserve_physical_area_signed_p1_fields_and_geometric_resolution() {
    let mesh = rectangle();
    let map: Vec<_> = (0..3 * mesh.nodes.len()).map(Some).collect();
    // Deliberately unrelated slope DOFs must not leak into the declared P1 trace.
    let phi: Vec<_> = mesh
        .nodes
        .iter()
        .flat_map(|&(x, _)| [x - 0.2, 777.0, -51.0])
        .collect();
    let scaled: Vec<_> = phi.iter().map(|v| -4.0 * v).collect();
    for level in 0..=2 {
        let (weights, mut surface, diameter) = samples(&mesh, level, 2);
        append_mode(&mesh, &map, &phi, &weights, &mut surface).unwrap();
        append_mode(&mesh, &map, &scaled, &weights, &mut surface).unwrap();
        close(diameter, 0.5 / (1usize << level) as f64, 0.5);
        close(surface.iter().map(|p| p.area_m2).sum(), 0.12, 0.12);
        close(
            surface.iter().map(|p| p.area_m2 * p.mode_shape[0]).sum(),
            0.0,
            0.024,
        );
        // Exact affine kinetic/overlap moment, not only cancellation of a mean.
        close(
            surface
                .iter()
                .map(|p| p.area_m2 * p.mode_shape[0].powi(2))
                .sum(),
            0.3 * 0.4_f64.powi(3) / 12.0,
            0.002,
        );
        for p in &surface {
            assert!(p.area_m2 > 0.0);
            assert_eq!(p.position_m[2], 0.0);
            close(p.mode_shape[0], p.position_m[0] - 0.2, 0.4);
            close(p.mode_shape[1], -4.0 * p.mode_shape[0], 1.6);
        }
        // Real coordinate changes must reach points, area and the true largest
        // edge. An area-derived diameter or origin-centred surrogate fails this.
        let mut moved = mesh.clone();
        for p in &mut moved.nodes {
            *p = (3.0 * p.0 + 1.2, 3.0 * p.1 - 0.4);
        }
        let (_, mut transformed, moved_diameter) = samples(&moved, level, 1);
        append_mode(&moved, &map, &phi, &weights, &mut transformed).unwrap();
        close(moved_diameter, 3.0 * diameter, 1.5);
        for (a, b) in surface.iter().zip(&transformed) {
            close(b.position_m[0], 3.0 * a.position_m[0] + 1.2, 3.0);
            close(b.position_m[1], 3.0 * a.position_m[1] - 0.4, 1.0);
            close(b.area_m2, 9.0 * a.area_m2, 9.0 * a.area_m2);
            assert_eq!(a.mode_shape[0], b.mode_shape[0]);
        }
        // A supported w coordinate contributes zero, even if a supplied phi
        // has a nonzero value at the coincident unreduced index.
        let mut supported = map.clone();
        supported[0] = None;
        let (_, mut fixed, _) = samples(&mesh, level, 1);
        append_mode(
            &mesh,
            &supported,
            &vec![1.0; map.len()],
            &weights,
            &mut fixed,
        )
        .unwrap();
        close(
            fixed.iter().map(|p| p.area_m2 * p.mode_shape[0]).sum(),
            0.08,
            0.12,
        );
    }
}

#[test]
fn zero_monopole_surface_retains_its_nonzero_finite_wavelength_directivity() {
    let mesh = rectangle();
    let map: Vec<_> = (0..3 * mesh.nodes.len()).map(Some).collect();
    let phi: Vec<_> = mesh
        .nodes
        .iter()
        .flat_map(|&(x, _)| [x - 0.2, 0.0, 0.0])
        .collect();
    // Exact far-field spatial Fourier factor for an antisymmetric affine
    // rectangle: int[-a,a] x*exp(i*k*x) dx times its physical height.
    // Net volume flow is zero; finite tangential wavenumber is not silent.
    let k = 9.0_f64;
    let a = 0.2;
    let expected_im = 2.0 * 0.3 * ((k * a).sin() / (k * k) - a * (k * a).cos() / k);
    assert!(expected_im.abs() > 1e-3);
    let mut errors = Vec::new();
    for level in 0..=2 {
        let (weights, mut surface, _) = samples(&mesh, level, 1);
        append_mode(&mesh, &map, &phi, &weights, &mut surface).unwrap();
        let mut real = 0.0;
        let mut imag = 0.0;
        for p in &surface {
            let phase = k * (p.position_m[0] - a);
            let amplitude = p.area_m2 * p.mode_shape[0];
            real += amplitude * phase.cos();
            imag += amplitude * phase.sin();
        }
        close(real, 0.0, expected_im.abs());
        errors.push((imag / expected_im - 1.0).abs());
    }
    assert!(errors[1] < 0.2 * errors[0], "{errors:?}");
    assert!(errors[2] < 0.2 * errors[1], "{errors:?}");
    assert!(errors[2] < 1e-3, "{errors:?}");
}

#[test]
fn complete_surface_work_is_bounded_before_rule_or_shape_allocation() {
    assert_eq!(surface_work(40, 64, 2).unwrap(), (4, 48, 1920));
    assert_eq!(
        surface_work(1, MAX_RAYLEIGH_INPUTS, 6).unwrap(),
        (64, 12_288, 12_288)
    );
    assert_eq!(
        surface_work(MAX_PLATE_SURFACE_SAMPLES / 3, 1, 0).unwrap().2,
        16_383
    );
    for (triangles, modes, level) in [
        (0, 1, 0),
        (1, 0, 0),
        (1, MAX_RAYLEIGH_INPUTS + 1, 0),
        (1, 1, 7),
        (40, 64, 4),
        (MAX_PLATE_SURFACE_SAMPLES / 3 + 1, 1, 0),
        (usize::MAX, 1, 0),
        (1, 1, usize::MAX),
    ] {
        assert!(surface_work(triangles, modes, level).is_err());
    }
}

fn physical_coefficients(body: &CompactBody) -> [u64; 8] {
    [
        body.area_m2,
        body.mass_kg,
        body.drive_participation,
        body.omega,
        body.zeta,
        body.y,
        body.v,
        body.piston_area_m2,
    ]
    .map(f64::to_bits)
}

#[test]
fn one_shared_harvest_preserves_every_body_coefficient_and_signed_surface_integral() {
    use fs_plate::{AssemblyOptions, EdgeSupport, PlateSection};
    let mesh = PlateMesh::rectangle(0.4, 0.3, 4, 3);
    let section = PlateSection::isotropic(12e9, 0.25, 0.003, 600.0).unwrap();
    let mut force = vec![0.0; mesh.nodes.len()];
    force[6] = 1.0;
    let chart = PlateChart::from_mesh(mesh, section).unwrap();
    let options = PlateChartRadiation {
        assembly: AssemblyOptions {
            support: EdgeSupport::SimplySupported,
            pretension: 0.0,
        },
        eigenvalue_window: (1.0, 6e6),
        n_modes: 3,
        damping_ratio: 0.02,
        unit_force_weights: force,
    };
    for mass in [PlateMassModel::Lumped, PlateMassModel::ConsistentDkt] {
        let compact =
            super::super::certified_chart_radiators_with_mass(&chart, &options, mass).unwrap();
        assert_eq!(compact.len(), 3);
        for level in 0..=2 {
            let image = certified_chart_surface_radiators(&chart, &options, mass, level).unwrap();
            assert_eq!(
                image
                    .bodies()
                    .iter()
                    .map(physical_coefficients)
                    .collect::<Vec<_>>(),
                compact
                    .iter()
                    .map(physical_coefficients)
                    .collect::<Vec<_>>()
            );
            assert_eq!(
                image.surface().len(),
                3 * 4usize.pow(level as u32) * chart.mesh.tris.len()
            );
            for (mode, body) in image.bodies().iter().enumerate() {
                let signed = image
                    .surface()
                    .iter()
                    .map(|s| s.area_m2 * s.mode_shape[mode])
                    .sum();
                let magnitude: f64 = image
                    .surface()
                    .iter()
                    .map(|s| s.area_m2 * s.mode_shape[mode].abs())
                    .sum();
                close(signed, body.area_m2, magnitude);
            }
            let (bodies, surface) = image.into_parts();
            assert!(surface.iter().all(|p| p.mode_shape.len() == bodies.len()));
        }
    }
}

#[test]
fn modal_velocity_stays_defined_for_a_mechanically_excited_dark_monopole() {
    let mut body = CompactBody::from_mass_normalized_projection(1000.0, 0.0, 2.0, 0.0).unwrap();
    let dt = 1e-5;
    body.drive(2.0, dt).unwrap();
    assert_eq!(body.volume_velocity(), 0.0);
    let expected = 2.0 / 1000.0 * (1000.0_f64 * dt).sin();
    close(body.modal_velocity_m_s(), expected, expected);
}
