//! Admit one original tetrahedral approximation of a complete cylinder or tube.
//!
//! This checks its fixed occupancy, not moving facets or field convergence.
//! Radial bounds cover entire facets and cells, including projected chords that
//! could otherwise bridge an annular hole despite admissible vertex radii.

use super::{CylindricalConductor, SpatialRotorError};
use fs_em::spatial::vector_potential::induction::SpatialInduction;
use fs_ga::Vec3;
use fs_math::dd::Dd;
use std::collections::{BTreeMap, BTreeSet};

type Result<T> = core::result::Result<T, SpatialRotorError>;

pub(super) fn admit_cylinder(
    field: &SpatialInduction<'_>,
    cells: &[usize],
    origin: Vec3,
    shape: CylindricalConductor,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64> {
    let CylindricalConductor {
        inner_radius_m: inner,
        outer_radius_m: outer,
        half_length_m: half,
        maximum_boundary_error_m: tolerance,
    } = shape;
    if [inner, outer, half, tolerance, origin.x, origin.y, origin.z]
        .into_iter()
        .any(|v| !v.is_finite())
        || inner < 0.0
        || outer <= inner
        || half <= 0.0
        || tolerance <= 0.0
        || tolerance >= 0.125 * (outer - inner).min(half)
        || inner > 0.0 && tolerance >= 0.125 * inner
    {
        return Err(invalid(
            "bounded positive cylindrical approximation required",
        ));
    }
    let problem = field.problem();
    let mesh = problem.mesh();
    if cells.is_empty() || cells.len() > mesh.tets.len() {
        return Err(invalid("nonempty original cylindrical conductor required"));
    }
    let point = |v: u32| {
        let p = problem.positions_m()[v as usize];
        Vec3::new(p[0] - origin.x, p[1] - origin.y, p[2] - origin.z)
    };
    let mut unique = BTreeSet::new();
    let mut faces = BTreeMap::<[u32; 3], Vec<(usize, [u32; 3])>>::new();
    for (local, &cell) in cells.iter().enumerate() {
        poll(cancelled)?;
        let tet = *mesh
            .tets
            .get(cell)
            .ok_or(invalid("rotor cell outside mesh"))?;
        if !unique.insert(cell) || !field.conductors()[cell].is_some() {
            return Err(invalid("distinct original conducting rotor cells required"));
        }
        let p = tet.map(point);
        let determinant = finite((p[1] - p[0]).dot((p[2] - p[0]).cross(p[3] - p[0])))?;
        if determinant <= 0.0
            || p.iter().any(|p| {
                let r = p.x.hypot(p.y);
                !r.is_finite()
                    || !p.z.is_finite()
                    || r > outer + tolerance
                    || p.z.abs() > half + tolerance
            })
        {
            return Err(invalid(
                "rotor cells must be positively oriented inside cylinder",
            ));
        }
        // Outward faces of the positively oriented original tetrahedron.
        for local_face in [[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]] {
            let face = local_face.map(|i| tet[i]);
            // In two dimensions every point of a four-vertex convex hull lies
            // in some three-vertex hull. All four projections therefore also
            // bound the minimum radius throughout the original tetrahedron.
            if inner > 0.0 && minimum_radius(local_face.map(|i| p[i]))? < inner - tolerance {
                return Err(invalid(
                    "rotor tetrahedron crosses the declared cylindrical hole",
                ));
            }
            let mut key = face;
            key.sort_unstable();
            let entries = faces.entry(key).or_default();
            entries.push((local, face));
            if entries.len() > 2 {
                return Err(invalid("nonmanifold rotor face"));
            }
        }
    }
    let mut adjacency = vec![Vec::new(); cells.len()];
    let mut boundary_edges = BTreeMap::<[u32; 2], (usize, i32)>::new();
    let mut links = BTreeMap::<u32, BTreeMap<u32, u32>>::new();
    let mut kinds = [0usize; 4]; // outer, inner, bottom, top
    let mut areas = [Dd::ZERO; 2];
    let mut maximum_error: f64 = 0.0;
    let mut boundary_faces = 0usize;
    for entries in faces.values() {
        poll(cancelled)?;
        let face = match entries.as_slice() {
            &[(a, u), (b, v)]
                if (0..3)
                    .any(|i| u[0] == v[i] && u[1] == v[(i + 2) % 3] && u[2] == v[(i + 1) % 3]) =>
            {
                adjacency[a].push(b);
                adjacency[b].push(a);
                continue;
            }
            &[(_, face)] => face,
            _ => {
                return Err(invalid(
                    "multiply covered or inconsistently oriented rotor face",
                ));
            }
        };
        let p = face.map(point);
        let normal = (p[1] - p[0]).cross(p[2] - p[0]);
        let bottom = p.iter().map(|p| (p.z + half).abs()).fold(0.0, f64::max);
        let top = p.iter().map(|p| (p.z - half).abs()).fold(0.0, f64::max);
        let (kind, error) = if bottom <= tolerance && finite(normal.z)? < 0.0 {
            areas[0] = areas[0] + Dd::from_f64(-0.5 * normal.z);
            (2, bottom)
        } else if top <= tolerance && finite(normal.z)? > 0.0 {
            areas[1] = areas[1] + Dd::from_f64(0.5 * normal.z);
            (3, top)
        } else {
            let least = minimum_radius(p)?;
            let greatest = p.iter().map(|p| p.x.hypot(p.y)).fold(0.0, f64::max);
            let radial = finite(
                normal.x * (p[0].x + p[1].x + p[2].x) + normal.y * (p[0].y + p[1].y + p[2].y),
            )?;
            let outside = (least - outer).abs().max((greatest - outer).abs());
            let inside = (least - inner).abs().max((greatest - inner).abs());
            if outside <= tolerance && radial > 0.0 {
                (0, outside)
            } else if inner > 0.0 && inside <= tolerance && radial < 0.0 {
                (1, inside)
            } else {
                return Err(invalid(
                    "rotor boundary is not a declared cylindrical facet or end",
                ));
            }
        };
        kinds[kind] += 1;
        boundary_faces += 1;
        maximum_error = maximum_error.max(error);
        for i in 0..3 {
            let (a, b, c) = (face[i], face[(i + 1) % 3], face[(i + 2) % 3]);
            let edge = boundary_edges.entry([a.min(b), a.max(b)]).or_default();
            edge.0 += 1;
            edge.1 += if a < b { 1 } else { -1 };
            if links.entry(a).or_default().insert(b, c).is_some() {
                return Err(invalid("nonmanifold cylindrical boundary vertex"));
            }
        }
    }
    if kinds[0] == 0
        || kinds[2] == 0
        || kinds[3] == 0
        || (kinds[1] > 0) != (inner > 0.0)
        || boundary_edges
            .values()
            .any(|&(n, sign)| n != 2 || sign != 0)
    {
        return Err(invalid(
            "complete closed oriented cylindrical boundary required",
        ));
    }
    // Each boundary-vertex link must be ONE circle; edge pairing alone admits
    // pinched vertices. The admitted solid/tube surfaces have genus zero/one.
    for link in links.values() {
        let start = *link
            .keys()
            .next()
            .ok_or(invalid("empty boundary vertex link"))?;
        let mut vertex = start;
        for count in 0..link.len() {
            poll(cancelled)?;
            vertex = *link
                .get(&vertex)
                .ok_or(invalid("open boundary vertex link"))?;
            if vertex == start && count + 1 != link.len() {
                return Err(invalid("disconnected boundary vertex link"));
            }
        }
        if vertex != start {
            return Err(invalid("unclosed boundary vertex link"));
        }
    }
    let euler = links.len() as i128 - boundary_edges.len() as i128 + boundary_faces as i128;
    if euler != if inner == 0.0 { 2 } else { 0 } {
        return Err(invalid("rotor boundary has undeclared topology"));
    }
    // A double covering can pass local incidence checks. Each complete cap's
    // projected area must lie in the analytic annulus's radial error envelope.
    let lower = finite(
        core::f64::consts::PI
            * ((outer - tolerance).powi(2)
                - if inner == 0.0 {
                    0.0
                } else {
                    (inner + tolerance).powi(2)
                }),
    )?;
    let upper = finite(
        core::f64::consts::PI
            * ((outer + tolerance).powi(2) - (inner - tolerance).max(0.0).powi(2)),
    )?;
    for area in areas {
        let area = finite(area.to_f64())?;
        if area < lower || area > upper {
            return Err(invalid(
                "cylindrical cap does not cover its declared section once",
            ));
        }
    }
    let mut visited = vec![false; cells.len()];
    let mut pending = vec![0];
    visited[0] = true;
    while let Some(cell) = pending.pop() {
        poll(cancelled)?;
        for &other in &adjacency[cell] {
            if !visited[other] {
                visited[other] = true;
                pending.push(other);
            }
        }
    }
    if visited.iter().any(|&v| !v) {
        return Err(invalid(
            "rotor conductor must be one face-connected cylinder",
        ));
    }
    poll(cancelled)?;
    Ok(maximum_error)
}

fn minimum_radius(p: [Vec3; 3]) -> Result<f64> {
    let cross = |a: Vec3, b: Vec3| finite(a.x * b.y - a.y * b.x);
    let sides = [cross(p[0], p[1])?, cross(p[1], p[2])?, cross(p[2], p[0])?];
    if cross(p[1] - p[0], p[2] - p[0])? != 0.0
        && (sides.iter().all(|&s| s >= 0.0) || sides.iter().all(|&s| s <= 0.0))
    {
        return Ok(0.0);
    }
    let mut nearest = f64::INFINITY;
    for i in 0..3 {
        let a = p[i];
        let d = p[(i + 1) % 3] - a;
        let square = finite(d.x * d.x + d.y * d.y)?;
        let t = if square == 0.0 {
            0.0
        } else {
            (finite(-a.x * d.x - a.y * d.y)? / square).clamp(0.0, 1.0)
        };
        nearest = nearest.min(finite((a.x + t * d.x).hypot(a.y + t * d.y))?);
    }
    Ok(nearest)
}

fn invalid(reason: &'static str) -> SpatialRotorError {
    SpatialRotorError::Invalid(reason)
}
fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("nonfinite cylindrical geometry"))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(SpatialRotorError::Cancelled)
    } else {
        Ok(())
    }
}
