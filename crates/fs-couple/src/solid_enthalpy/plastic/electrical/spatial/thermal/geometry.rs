//! Native principal-normal stretch of original isotropic planar heat paths.
//! Geometry changes transport, never the reference caloric mass or stored heat.

use super::*;
use crate::circuit_enthalpy::contacts::heat::ThermalLinkMetric;
use crate::circuit_enthalpy::material::ConductionLink;
use crate::solid_enthalpy::same_measure;
use fs_math::det::sqrt;
use fs_solid::loads3::{LoadError3, TetLoadGeometry3};
use fs_solid::plastic3::{SmallStrainLawField, SmallStrainTetError, SmallStrainTetProblem};

/// Bounds for the explicit planar-half-cell thermal geometry reduction.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpatialThermalGeometrySolve {
    /// Maximum Frobenius norm of grad(u), in (0,1). Native material
    /// small-strain and rotation bounds remain independently authoritative.
    pub maximum_displacement_gradient: f64,
    /// Relative principal-direction defect |F^T F N - lambda^2 N|/lambda^2,
    /// in (0,1). This admits numerical roundoff, not a sheared Fourier model.
    pub orthogonality_tolerance: f64,
}

pub(in crate::solid_enthalpy::plastic::electrical::spatial) struct MetricRoute {
    link: u64,
    elements: [usize; 2],
    vertices: [[usize; 3]; 2],
    normals: [[f64; 3]; 2],
    area_m2: f64,
    distances_m: [f64; 2],
    controls: SpatialThermalGeometrySolve,
}

pub(in crate::solid_enthalpy::plastic::electrical::spatial) fn routes<E>(
    solid: &ThermoplasticSolid<'_>,
    problem: &SpatialConduction,
    thermal_routes: &[ThermalPressureRoute<'_>],
    paths: &[ConductionLink],
    controls: Option<SpatialThermalGeometrySolve>,
    cx: &Cx<'_>,
) -> Result<Vec<MetricRoute>, ElectricalSolidSourceError<E>> {
    let Some(controls) = controls else {
        return Ok(Vec::new());
    };
    let invalid = ElectricalSolidSourceError::Invalid;
    if thermal_routes.is_empty()
        || !positive(controls.maximum_displacement_gradient)
        || controls.maximum_displacement_gradient >= 1.0
        || !positive(controls.orthogonality_tolerance)
        || controls.orthogonality_tolerance >= 1.0
    {
        return Err(invalid(
            "deforming heat paths require mapped interfaces and bounded positive geometry controls",
        ));
    }
    let mut routes = Vec::with_capacity(thermal_routes.len());
    for route in thermal_routes {
        cx.checkpoint().map_err(|_| DepositionError::Cancelled)?;
        let contact = &problem.contacts()[route.contact];
        let elements = problem.contact_cells()[route.contact];
        let vertices = [contact.first_vertices, contact.second_vertices];
        let path = paths
            .iter()
            .find(|path| path.link == route.link)
            .ok_or(invalid(
                "deforming contact must retain its original two-sided Fourier path",
            ))?;
        let distances_m = [path.distance_a_m, path.distance_b_m];
        let mut normals = [[0.0; 3]; 2];
        let mut areas = [0.0; 2];
        for side in 0..2 {
            let mut nodes = vertices[side];
            nodes.sort_unstable();
            let face = solid
                .geometry
                .faces()
                .iter()
                .find(|face| face.nodes == nodes)
                .ok_or(invalid(
                    "deforming thermal path lost its original exterior face",
                ))?;
            let distance = 0.75 * solid.cells[elements[side]].volume_m3 / face.area_m2;
            if !same_measure(distances_m[side], distance) {
                return Err(invalid(
                    "thermal half-path must be the original tetrahedron centroid-to-face normal distance",
                ));
            }
            normals[side] = face.outward_normal;
            areas[side] = face.area_m2;
        }
        if !same_measure(areas[0], areas[1]) {
            return Err(invalid(
                "paired thermal paths require the same original matching face measure",
            ));
        }
        routes.push(MetricRoute {
            link: route.link,
            elements,
            vertices,
            normals,
            area_m2: areas[0],
            distances_m,
            controls,
        });
    }
    Ok(routes)
}

/// Extract accepted gradients from stored u with the SAME native P1 owner.
/// Target-load mechanics must never substitute for accepted geometry.
pub(in crate::solid_enthalpy::plastic::electrical::spatial) fn accepted<E>(
    solid: &ThermoplasticSolid<'_>,
    routes: &[MetricRoute],
    positions: &[[f64; 3]],
    prescribed_m: &[(usize, f64)],
    settings: IncrementSettings,
    cx: &Cx<'_>,
) -> Result<Vec<ThermalLinkMetric>, ElectricalSolidSourceError<E>> {
    if routes.is_empty() {
        return Ok(Vec::new());
    }
    let native = SmallStrainTetProblem {
        nodes_m: solid.nodes_m,
        tetrahedra: solid.tetrahedra,
        laws: SmallStrainLawField::PerElement(&solid.materials),
        prescribed_m,
        nodal_forces_n: &[],
        load_path: &[1.0],
        budget: solid.budget,
        settings,
    };
    let gradients = native
        .displacement_gradients(solid.solution.displacement_m.as_flattened(), cx)
        .map_err(|error| {
            ElectricalSolidSourceError::Mechanical(match error {
                SmallStrainTetError::Cancelled => ThermoplasticError::Cancelled,
                other => ThermoplasticError::Solid(other),
            })
        })?;
    sample(solid, routes, positions, &gradients, cx)
}

pub(in crate::solid_enthalpy::plastic::electrical::spatial) fn sample<E>(
    solid: &ThermoplasticSolid<'_>,
    routes: &[MetricRoute],
    positions: &[[f64; 3]],
    gradients: &[[[f64; 3]; 3]],
    cx: &Cx<'_>,
) -> Result<Vec<ThermalLinkMetric>, ElectricalSolidSourceError<E>> {
    if routes.is_empty() {
        return Ok(Vec::new());
    }
    let invalid = ElectricalSolidSourceError::Invalid;
    let maximum_faces = solid
        .tetrahedra
        .len()
        .checked_mul(4)
        .ok_or(invalid("thermal metric face budget overflow"))?;
    let geometry =
        TetLoadGeometry3::new(positions, solid.tetrahedra, solid.budget, maximum_faces, cx)
            .map_err(|error| {
                ElectricalSolidSourceError::Mechanical(match error {
                    LoadError3::Cancelled => ThermoplasticError::Cancelled,
                    other => ThermoplasticError::Geometry(other),
                })
            })?;
    let mut metrics = Vec::with_capacity(routes.len());
    for route in routes {
        cx.checkpoint().map_err(|_| DepositionError::Cancelled)?;
        if route.vertices[0]
            .iter()
            .zip(route.vertices[1])
            .any(|(&a, b)| positions[a] != positions[b])
        {
            return Err(invalid(
                "thermal metric requires the original exactly matching material traces",
            ));
        }
        let mut areas = [0.0; 2];
        let mut distances_m = [0.0; 2];
        let mut normals = [[0.0; 3]; 2];
        for side in 0..2 {
            let gradient = gradients.get(route.elements[side]).ok_or(invalid(
                "thermal path lost its native displacement gradient",
            ))?;
            let norm2: f64 = gradient.iter().flatten().map(|g| g * g).sum();
            if !norm2.is_finite() || norm2 > route.controls.maximum_displacement_gradient.powi(2) {
                return Err(invalid(
                    "thermal displacement gradient exceeds its admitted small-strain geometry bound",
                ));
            }
            let mut f = *gradient;
            for (i, row) in f.iter_mut().enumerate() {
                row[i] += 1.0;
            }
            let j = f[0][0] * (f[1][1] * f[2][2] - f[1][2] * f[2][1])
                - f[0][1] * (f[1][0] * f[2][2] - f[1][2] * f[2][0])
                + f[0][2] * (f[1][0] * f[2][1] - f[1][1] * f[2][0]);
            let n = route.normals[side];
            let tangent = f.map(|row| dot(row, n));
            let stretch2 = dot(tangent, tangent);
            let residual = std::array::from_fn(|a| {
                f[0][a] * tangent[0] + f[1][a] * tangent[1] + f[2][a] * tangent[2] - stretch2 * n[a]
            });
            if !positive(j)
                || !positive(stretch2)
                || !dot(residual, residual).is_finite()
                || sqrt(dot(residual, residual)) > route.controls.orthogonality_tolerance * stretch2
            {
                return Err(invalid(
                    "deforming scalar heat path requires positive volume and a principal reference normal",
                ));
            }
            let stretch = sqrt(stretch2);
            let mut nodes = route.vertices[side];
            nodes.sort_unstable();
            let face = geometry
                .faces()
                .iter()
                .find(|face| face.nodes == nodes)
                .ok_or(invalid(
                    "deformed thermal trace lost its original exterior ownership",
                ))?;
            let distance = route.distances_m[side] * stretch;
            if !same_measure(face.area_m2, route.area_m2 * j / stretch)
                || !same_measure(
                    distance,
                    0.75 * geometry.element_volumes_m3()[route.elements[side]] / face.area_m2,
                )
            {
                return Err(invalid(
                    "principal thermal path must agree with actual native face area and centroid distance",
                ));
            }
            areas[side] = face.area_m2;
            distances_m[side] = distance;
            normals[side] = face.outward_normal;
        }
        if !same_measure(areas[0], areas[1]) || dot(normals[0], normals[1]) > -1.0 + 1e-12 {
            return Err(invalid(
                "deformed thermal traces require one matching area and opposite adjacent volumes",
            ));
        }
        metrics.push(ThermalLinkMetric {
            link: route.link,
            area_m2: areas[0],
            distances_m,
        });
    }
    Ok(metrics)
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
