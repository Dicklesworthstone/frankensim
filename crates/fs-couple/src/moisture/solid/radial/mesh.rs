//! Conforming cylindrical shell tessellation for the native hygroscopic solid.
//! Each triangular prism uses the same globally ordered face diagonal, including
//! the angular seam. Curved/polygonal volume discrepancy is admitted explicitly.

#[derive(Clone, Debug)]
pub(super) struct CylinderMesh {
    pub nodes_m: Vec<[f64; 3]>,
    pub tetrahedra: Vec<[usize; 4]>,
    pub shell_for_element: Vec<usize>,
    pub element_volumes_m3: Vec<f64>,
    pub shell_volumes_m3: Vec<f64>,
    pub maximum_relative_volume_error: f64,
}

#[derive(Debug)]
pub(super) enum MeshError<E> {
    Invalid(&'static str),
    Capacity,
    Checkpoint(E),
}

fn reserved<T, E>(count: usize) -> Result<Vec<T>, MeshError<E>> {
    let mut result = Vec::new();
    result
        .try_reserve_exact(count)
        .map_err(|_| MeshError::Capacity)?;
    Ok(result)
}

fn volume(nodes: &[[f64; 3]], tet: [usize; 4]) -> f64 {
    let a: [f64; 3] = core::array::from_fn(|i| nodes[tet[1]][i] - nodes[tet[0]][i]);
    let b: [f64; 3] = core::array::from_fn(|i| nodes[tet[2]][i] - nodes[tet[0]][i]);
    let c: [f64; 3] = core::array::from_fn(|i| nodes[tet[3]][i] - nodes[tet[0]][i]);
    (a[0] * (b[1] * c[2] - b[2] * c[1])
        + a[1] * (b[2] * c[0] - b[0] * c[2])
        + a[2] * (b[0] * c[1] - b[1] * c[0]))
        / 6.0
}

#[allow(clippy::too_many_arguments)]
pub(super) fn cylinder<E>(
    radii: &[f64],
    length: f64,
    angular: usize,
    axial: usize,
    maximum_nodes: usize,
    maximum_tetrahedra: usize,
    maximum_relative_volume_error: f64,
    mut checkpoint: impl FnMut() -> Result<(), E>,
) -> Result<CylinderMesh, MeshError<E>> {
    checkpoint().map_err(MeshError::Checkpoint)?;
    if radii.len() < 2
        || angular < 3
        || axial == 0
        || !length.is_finite()
        || length <= 0.0
        || radii.iter().any(|r| !r.is_finite() || *r < 0.0)
        || radii.windows(2).any(|r| r[0] >= r[1])
        || !maximum_relative_volume_error.is_finite()
        || maximum_relative_volume_error <= 0.0
        || maximum_relative_volume_error >= 1.0
    {
        return Err(MeshError::Invalid(
            "ordered radial faces, positive cylinder and bounded geometric error required",
        ));
    }
    let center = usize::from(radii[0] == 0.0);
    let rings = radii.len() - center;
    let per_plane = rings
        .checked_mul(angular)
        .and_then(|n| n.checked_add(center))
        .ok_or(MeshError::Capacity)?;
    let nodes = axial
        .checked_add(1)
        .and_then(|n| n.checked_mul(per_plane))
        .ok_or(MeshError::Capacity)?;
    let triangles = (radii.len() - 1)
        .checked_mul(2)
        .and_then(|n| n.checked_sub(center))
        .and_then(|n| n.checked_mul(angular))
        .ok_or(MeshError::Capacity)?;
    let elements = triangles
        .checked_mul(axial)
        .and_then(|n| n.checked_mul(3))
        .ok_or(MeshError::Capacity)?;
    if nodes > maximum_nodes || elements > maximum_tetrahedra {
        return Err(MeshError::Invalid(
            "radial solid exceeds native mesh budget",
        ));
    }
    let mut mesh = CylinderMesh {
        nodes_m: reserved(nodes)?,
        tetrahedra: reserved(elements)?,
        shell_for_element: reserved(elements)?,
        element_volumes_m3: reserved(elements)?,
        shell_volumes_m3: reserved(radii.len() - 1)?,
        maximum_relative_volume_error: 0.0,
    };
    mesh.shell_volumes_m3.resize(radii.len() - 1, 0.0);
    for z in 0..=axial {
        checkpoint().map_err(MeshError::Checkpoint)?;
        let height = length * (z as f64 / axial as f64);
        if center == 1 {
            mesh.nodes_m.push([0.0, 0.0, height]);
        }
        for &r in &radii[center..] {
            for angle in 0..angular {
                if angle % 128 == 0 {
                    checkpoint().map_err(MeshError::Checkpoint)?;
                }
                let theta = core::f64::consts::TAU * (angle as f64 / angular as f64);
                let (s, c) = fs_math::det::sin_cos(theta);
                mesh.nodes_m.push([r * c, r * s, height]);
            }
        }
    }
    let ring = |face: usize, angle: usize| center + (face - center) * angular + angle;
    for z in 0..axial {
        for shell in 0..radii.len() - 1 {
            for angle in 0..angular {
                if angle % 64 == 0 {
                    checkpoint().map_err(MeshError::Checkpoint)?;
                }
                let next = (angle + 1) % angular;
                let outer = [ring(shell + 1, angle), ring(shell + 1, next)];
                let (triangles, count) = if center == 1 && shell == 0 {
                    ([[0, outer[0], outer[1]], [0; 3]], 1)
                } else {
                    let inner = [ring(shell, angle), ring(shell, next)];
                    (
                        [
                            [inner[0], outer[0], outer[1]],
                            [inner[0], outer[1], inner[1]],
                        ],
                        2,
                    )
                };
                for mut triangle in triangles.into_iter().take(count) {
                    triangle.sort_unstable();
                    let a = triangle.map(|i| i + z * per_plane);
                    let b = a.map(|i| i + per_plane);
                    for mut tet in [
                        [a[0], a[1], a[2], b[2]],
                        [a[0], a[1], b[1], b[2]],
                        [a[0], b[0], b[1], b[2]],
                    ] {
                        let signed = volume(&mesh.nodes_m, tet);
                        if !signed.is_finite() || signed == 0.0 {
                            return Err(MeshError::Invalid(
                                "radial solid tetrahedron is unrepresentable",
                            ));
                        }
                        if signed < 0.0 {
                            tet.swap(1, 2);
                        }
                        let v = signed.abs();
                        mesh.tetrahedra.push(tet);
                        mesh.shell_for_element.push(shell);
                        mesh.element_volumes_m3.push(v);
                        mesh.shell_volumes_m3[shell] += v;
                    }
                }
            }
        }
    }
    for (shell, &actual) in mesh.shell_volumes_m3.iter().enumerate() {
        checkpoint().map_err(MeshError::Checkpoint)?;
        let [a, b] = [radii[shell], radii[shell + 1]];
        let exact = core::f64::consts::PI * length * (b - a) * (b + a);
        if !actual.is_finite() || actual <= 0.0 || !exact.is_finite() || exact <= 0.0 {
            return Err(MeshError::Invalid(
                "radial solid shell volume is unrepresentable",
            ));
        }
        let relative = (actual / exact - 1.0).abs();
        if !relative.is_finite() || relative > maximum_relative_volume_error {
            return Err(MeshError::Invalid(
                "refine angular mesh to meet curved shell volume tolerance",
            ));
        }
        mesh.maximum_relative_volume_error = mesh.maximum_relative_volume_error.max(relative);
    }
    checkpoint().map_err(MeshError::Checkpoint)?;
    Ok(mesh)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeMap;

    #[test]
    fn cylindrical_shell_mesh_closes_seams_and_retains_the_analytic_volume_error() {
        for radii in [vec![0.0, 0.3, 0.7, 1.0], vec![0.2, 0.35, 1.0]] {
            let mut errors = Vec::new();
            for angular in [8, 16] {
                let mesh = cylinder(&radii, 2.0, angular, 2, 10000, 50000, 0.11, || {
                    Ok::<(), ()>(())
                })
                .unwrap();
                let fraction = (core::f64::consts::TAU / angular as f64).sin() * angular as f64
                    / core::f64::consts::TAU;
                for (shell, &v) in mesh.shell_volumes_m3.iter().enumerate() {
                    let expected = core::f64::consts::PI
                        * 2.0
                        * (radii[shell + 1].powi(2) - radii[shell].powi(2))
                        * fraction;
                    assert!((v / expected - 1.0).abs() < 2e-14);
                }
                let mut faces = BTreeMap::new();
                for &tet in &mesh.tetrahedra {
                    assert!(volume(&mesh.nodes_m, tet) > 0.0);
                    for omitted in 0..4 {
                        let mut f = [0; 3];
                        let mut j = 0;
                        for (i, &vertex) in tet.iter().enumerate() {
                            if i != omitted {
                                f[j] = vertex;
                                j += 1;
                            }
                        }
                        f.sort_unstable();
                        *faces.entry(f).or_insert(0) += 1;
                    }
                }
                assert!(faces.values().all(|&count| count == 1 || count == 2));
                let rings = radii.len() - usize::from(radii[0] == 0.0);
                let end_triangles =
                    angular * (2 * (radii.len() - 1) - usize::from(radii[0] == 0.0));
                let walls = if radii[0] == 0.0 { 1 } else { 2 };
                let boundary = 2 * end_triangles + walls * 2 * angular * 2;
                assert_eq!(
                    faces.values().filter(|&&count| count == 1).count(),
                    boundary
                );
                assert_eq!(
                    mesh.nodes_m.len(),
                    3 * (rings * angular + usize::from(radii[0] == 0.0))
                );
                errors.push(mesh.maximum_relative_volume_error);
            }
            assert!(errors[0] > 3.8 * errors[1]);
        }
    }

    #[test]
    fn curved_geometry_work_and_accuracy_limits_refuse_before_publication() {
        assert!(cylinder(&[0.0, 1.0], 1.0, 8, 1, 100, 100, 0.01, || Ok::<(), ()>(())).is_err());
        assert!(
            cylinder(&[0.0, 1.0], 1.0, usize::MAX, 1, 100, 100, 0.1, || Ok::<
                (),
                (),
            >(
                ()
            ))
            .is_err()
        );
        let mut calls = 0;
        assert!(matches!(
            cylinder(&[0.0, 0.5, 1.0], 1.0, 64, 4, 10000, 50000, 0.1, || {
                calls += 1;
                if calls == 5 { Err("stop") } else { Ok(()) }
            }),
            Err(MeshError::Checkpoint("stop"))
        ));
    }
}
