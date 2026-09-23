//! An actual two-sided acoustic boundary from a reduced shell, not a monopole.
//!
//! The reference midsurface is extruded along area-weighted vertex directors.
//! Both faces and every boundary loop (including mounting holes) are retained.
//! Panel velocities include physical axial rotations: u_skin = u + theta x d.
//! The same weights transpose pressure work back into the mechanical basis.
//!
//! This is a cold, piecewise-linear, finite-thickness geometric approximation.
//! Refinement and close-panel quadrature belong to the downstream acoustic
//! solve. Local inversion and nonmanifold topology are refused; global offset
//! self-intersection is NOT certified. No radiation law, damping, time stepper,
//! fluid loading or moving-boundary acoustics is manufactured here.
use super::{PlateError, ShellReduction, bad, dot};
use std::collections::BTreeMap;

#[path = "surface.rs"]
mod surface;
pub use surface::{ShellFace, ShellSurfacePort};

/// Explicit cold-work ceilings, independent of an acoustic solver's limits.
#[derive(Debug, Clone, Copy)]
pub struct RadiationSurfaceBudget {
    /// Maximum triangles, including both faces and boundary walls.
    pub max_panels: usize,
    /// Maximum stored panel-by-mode normal-velocity weights.
    pub max_panel_modes: usize,
}

/// Immutable boundary geometry and mode-major normal-velocity projection.
/// Mode ordering and normalization are exactly those of the source reduction.
#[derive(Debug, Clone)]
pub struct ShellRadiationSurface {
    triangles: Vec<[[f64; 3]; 3]>,
    areas: Vec<f64>,
    weights: Vec<Vec<f64>>,
    boundary_panels: usize,
    volume_m3: f64,
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    core::array::from_fn(|c| a[c] - b[c])
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
}
fn reserve<T>(n: usize) -> Result<Vec<T>, PlateError> {
    let mut out = Vec::new();
    out.try_reserve_exact(n).map_err(|_| bad("acoustic skin allocation refused"))?;
    Ok(out)
}

impl ShellReduction {
    /// Extrude this reduction's own geometry and project its own modal basis.
    /// Nodal thickness must average to the section thickness on EVERY facet;
    /// a different mesh, thickness scale or unrelated modal list is not accepted.
    /// Normals are smoothed geometrically, not inferred from an instrument name.
    ///
    /// # Errors
    /// Invalid thickness, inconsistent/nonmanifold topology, local offset fold,
    /// nonfinite geometry/projection, allocation failure or exceeded work budget.
    pub fn radiation_surface(
        &self,
        nodal_thickness_m: &[f64],
        budget: RadiationSurfaceBudget,
    ) -> Result<ShellRadiationSurface, PlateError> {
        let nodes = self.nodes;
        let modes = self.mode_count();
        let faces = self.triangles.len();
        if nodal_thickness_m.len() != nodes
            || nodal_thickness_m.iter().any(|h| !h.is_finite() || *h <= 0.0)
            || faces.checked_mul(2).is_none_or(|n| n > budget.max_panels)
        {
            return Err(bad("acoustic skin thickness or face budget is invalid"));
        }
        let mut edges = BTreeMap::<(usize, usize), Vec<(usize, usize)>>::new();
        let mut links = vec![Vec::<(usize, usize)>::new(); nodes];
        let directors = surface::directors(self, nodal_thickness_m)?;
        for tri in &self.triangles {
            for a in 0..3 {
                let (i, j, k) = (tri[a], tri[(a+1)%3], tri[(a+2)%3]);
                edges.entry((i.min(j), i.max(j))).or_default().push((i, j));
                links[i].push((j, k));
            }
        }
        // An edge may belong to one boundary face or two oppositely wound faces.
        let mut boundary = Vec::new();
        for uses in edges.values() {
            match uses.as_slice() {
                &[edge] => boundary.push(edge),
                &[(a, b), (c, d)] if a == d && b == c => {}
                _ => return Err(bad("acoustic skin requires consistently oriented manifold edges")),
            }
        }
        // Edge manifoldness alone does not reject two fans touching at a vertex.
        for link in &links { validate_link(link)?; }
        let wall_count = boundary.len().checked_mul(2)
            .ok_or_else(|| bad("acoustic skin panel count overflow"))?;
        let count = faces.checked_mul(2).and_then(|n| n.checked_add(wall_count))
            .ok_or_else(|| bad("acoustic skin panel count overflow"))?;
        if count > budget.max_panels
            || count.checked_mul(modes).is_none_or(|n| n > budget.max_panel_modes)
        {
            return Err(bad("acoustic skin exceeds its panel/projection budget"));
        }
        // (midsurface node, sign of director offset), rather than an independent
        // acoustic mode table: side-wall motion cannot lose its mechanical origin.
        let mut panels = reserve::<[(usize, f64); 3]>(count)?;
        for tri in &self.triangles {
            panels.push(tri.map(|i| (i, 1.0)));
            panels.push([(tri[2], -1.0), (tri[1], -1.0), (tri[0], -1.0)]);
        }
        for (a, b) in boundary {
            panels.push([(a, 1.0), (a, -1.0), (b, -1.0)]);
            panels.push([(a, 1.0), (b, -1.0), (b, 1.0)]);
        }
        let mut triangles = reserve(count)?;
        let mut areas = reserve(count)?;
        let mut weights = reserve(modes)?;
        for _ in 0..modes { weights.push(reserve::<f64>(count)?); }
        let origin = self.reference_positions[0];
        let mut volume_m3 = 0.0;
        for (index, panel) in panels.iter().enumerate() {
            let triangle: [[f64; 3]; 3] = panel.map(|(node, sign)| {
                core::array::from_fn(|c| self.reference_positions[node][c] + sign * directors[node][c])
            });
            let area_normal = cross(sub(triangle[1], triangle[0]), sub(triangle[2], triangle[0]));
            let norm = dot(area_normal, area_normal).sqrt();
            if triangle.iter().flatten().any(|v| !v.is_finite()) || !norm.is_finite() || norm <= 0.0 {
                return Err(bad("acoustic skin has nonfinite or collapsed panels"));
            }
            if index < 2 * faces {
                let sign = if index % 2 == 0 { 1.0 } else { -1.0 };
                if sign * dot(area_normal, self.area_normals[index/2]) <= 0.0 {
                    return Err(bad("acoustic skin offset inverted a face"));
                }
            }
            let normal = area_normal.map(|v| v / norm);
            volume_m3 += dot(sub(triangle[0], origin), area_normal) / 6.0;
            for (mode, row) in weights.iter_mut().enumerate() {
                let mut velocity = [0.0; 3];
                for &(node, sign) in panel {
                    let offset = directors[node].map(|d| d * sign);
                    let rotation = cross(self.rotations[mode*nodes+node], offset);
                    for c in 0..3 {
                        velocity[c] += (self.translations[mode*nodes+node][c] + rotation[c]) / 3.0;
                    }
                }
                let value = dot(velocity, normal);
                if !value.is_finite() { return Err(bad("acoustic skin modal projection overflow")); }
                row.push(value);
            }
            triangles.push(triangle);
            areas.push(0.5 * norm);
        }
        if !volume_m3.is_finite() || volume_m3 <= 0.0 {
            return Err(bad("acoustic skin does not enclose positive finite material volume"));
        }
        Ok(ShellRadiationSurface { triangles, areas, weights, boundary_panels: wall_count, volume_m3 })
    }
}

fn validate_link(link: &[(usize, usize)]) -> Result<(), PlateError> {
    let mut adjacent = BTreeMap::<usize, Vec<usize>>::new();
    for &(a, b) in link {
        adjacent.entry(a).or_default().push(b);
        adjacent.entry(b).or_default().push(a);
    }
    let ends = adjacent.values().filter(|v| v.len() == 1).count();
    if adjacent.is_empty() || (ends != 0 && ends != 2)
        || adjacent.values().any(|v| v.len() > 2)
    {
        return Err(bad("acoustic skin vertex link is not a manifold chain or cycle"));
    }
    let mut pending = vec![*adjacent.keys().next().expect("nonempty checked")];
    let mut visited = std::collections::BTreeSet::new();
    while let Some(i) = pending.pop() {
        if visited.insert(i) { pending.extend(adjacent[&i].iter().copied()); }
    }
    if visited.len() != adjacent.len() { return Err(bad("acoustic skin has disconnected vertex fans")); }
    Ok(())
}

impl ShellRadiationSurface {
    /// Exact oriented triangles for the downstream closed-surface BEM constructor.
    #[must_use]
    pub fn triangles(&self) -> &[[[f64; 3]; 3]] { &self.triangles }
    /// Triangle areas [m^2] in the identical panel order.
    #[must_use]
    pub fn areas_m2(&self) -> &[f64] { &self.areas }
    /// Mode-major weights [1/sqrt(kg)]. Multiply each row by that mode's velocity
    /// [m sqrt(kg)/s] to obtain its contribution to outward panel velocity [m/s].
    #[must_use]
    pub fn normal_velocity_weights(&self) -> &[Vec<f64>] { &self.weights }
    /// Number of panels on all outer/inner boundary walls.
    #[must_use]
    pub fn boundary_panel_count(&self) -> usize { self.boundary_panels }
    /// Signed enclosed material volume, not the volume inside a drum or bowl.
    #[must_use]
    pub fn material_volume_m3(&self) -> f64 { self.volume_m3 }

    /// Allocation-free modal-to-panel velocity projection. Invalid input or an
    /// overflowing contraction leaves output unchanged. No modes are dropped.
    pub fn normal_velocity(&self, modal: &[f64], out: &mut [f64]) -> Result<(), PlateError> {
        if modal.len() != self.weights.len() || out.len() != self.areas.len()
            || modal.iter().any(|v| !v.is_finite())
        { return Err(bad("acoustic skin velocity shape/finite-set mismatch")); }
        let value = |p: usize| self.weights.iter().zip(modal).map(|(row, v)| row[p]*v).sum::<f64>();
        if (0..out.len()).any(|p| !value(p).is_finite()) {
            return Err(bad("acoustic skin velocity contraction overflow"));
        }
        for (p, out) in out.iter_mut().enumerate() { *out = value(p); }
        Ok(())
    }

    /// Exterior pressure [Pa] produces INWARD force on the solid. This is the
    /// negative area-weighted transpose of `normal_velocity`: solid port power
    /// is minus the acoustic pressure-times-outward-volume-flow power. It is a
    /// projection only, NOT a coupled acoustic load or a passivity certificate.
    /// Allocation-free; refusal leaves output unchanged.
    pub fn pressure_force(&self, pressure: &[f64], out: &mut [f64]) -> Result<(), PlateError> {
        if pressure.len() != self.areas.len() || out.len() != self.weights.len()
            || pressure.iter().any(|v| !v.is_finite())
        { return Err(bad("acoustic skin pressure shape/finite-set mismatch")); }
        let value = |k: usize| -self.weights[k].iter().zip(&self.areas).zip(pressure)
            .map(|((b, a), p)| b*a*p).sum::<f64>();
        if (0..out.len()).any(|k| !value(k).is_finite()) {
            return Err(bad("acoustic skin pressure contraction overflow"));
        }
        for (k, out) in out.iter_mut().enumerate() { *out = value(k); }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::shell::ShellMesh;

    fn budget() -> RadiationSurfaceBudget {
        RadiationSurfaceBudget { max_panels: 1000, max_panel_modes: 10000 }
    }
    // Geometric/kinematic fixture, NOT an eigensolve or mass-normalization proof.
    // Production can construct a ShellReduction only through the owning reducer.
    fn fixture(nodes: Vec<[f64; 3]>, triangles: Vec<[usize; 3]>, h: Vec<f64>) -> ShellReduction {
        let mesh = ShellMesh::new(nodes, triangles).unwrap();
        let count = mesh.nodes.len();
        let mut translations = Vec::new();
        let mut rotations = Vec::new();
        // Vertical translation and a physical rigid rotation about the y axis.
        for _ in &mesh.nodes { translations.push([0.0, 0.0, 1.0]); rotations.push([0.0; 3]); }
        for &p in &mesh.nodes { translations.push(cross([0.0, 1.0, 0.0], p)); rotations.push([0.0, 1.0, 0.0]); }
        ShellReduction {
            omegas: vec![0.0, 0.0], remainder: vec![0.0; 4], facets: vec![],
            translations, rotations, reference_positions: mesh.nodes.clone(),
            section_thicknesses: mesh.tris.iter().map(|t| t.iter().map(|&i| h[i]/3.0).sum()).collect(),
            nodes: count, triangles: mesh.tris.clone(),
            area_normals: (0..mesh.tris.len()).map(|f| {
                let g = mesh.facet(f).unwrap(); g.frame[2].map(|v| v*g.area_m2)
            }).collect(),
        }
    }
    fn square(h: Vec<f64>) -> ShellReduction {
        fixture(vec![[0.,0.,0.], [1.,0.,0.], [1.,1.,0.], [0.,1.,0.]], vec![[0,1,2], [0,2,3]], h)
    }
    #[test]
    fn skin_point_contact_and_acoustic_panel_share_geometry_motion_and_moments() {
        use super::surface::ShellFace;
        let h = vec![0.01, 0.02, 0.03, 0.015];
        let r = fixture(vec![[0.,0.,0.],[1.,0.,0.2],[1.,1.,0.3],[0.,1.,0.1]],
            vec![[0,1,2],[0,2,3]], h.clone());
        let skin = r.radiation_surface(&h, budget()).unwrap();
        for f in 0..2 { for (side, face) in [ShellFace::Positive, ShellFace::Negative].into_iter().enumerate() {
            let panel = skin.triangles[2*f+side];
            let area = cross(sub(panel[1],panel[0]), sub(panel[2],panel[0]));
            let normal = area.map(|x|x/dot(area,area).sqrt());
            let p = r.surface_point_port(&h,f,[1./3.;3],face,normal).unwrap();
            for c in 0..3 {
                assert!((p.position_m[c]-panel.iter().map(|v|v[c]/3.).sum::<f64>()).abs()<1e-14);
            }
            for m in 0..2 { assert!((p.weights[m]-skin.weights[m][2*f+side]).abs()<1e-14); }
            // Independent rigid-motion oracle, INCLUDING the offset lever arm.
            let qdot = [0.3,-0.7]; let force = 1.9;
            let velocity = dot(normal,[0.,0.,qdot[0]])
                + qdot[1]*dot(normal,cross([0.,1.,0.],p.position_m));
            let modal_work: f64 = qdot.iter().zip(&p.weights).map(|(v,b)|v*b*force).sum();
            assert!((modal_work-force*velocity).abs()<1e-14);
        }}
        let top = r.surface_point_port(&h,0,[0.2,0.3,0.5],ShellFace::Positive,[1.,0.,0.]).unwrap();
        let mid = r.point_port(0,[0.2,0.3,0.5],[1.,0.,0.]).unwrap();
        assert!((top.weights[1]-mid[1]).abs()>0.001, "offset moment must survive");
    }

    #[test]
    fn skin_point_preserves_source_barycentrics_and_refuses_incompatible_inputs() {
        use super::surface::ShellFace;
        let r = square(vec![0.02;4]); let h = [0.02;4];
        for face in [ShellFace::Positive,ShellFace::Negative] {
            let p = r.surface_point_port(&h,0,[0.,1.,0.],face,[0.,0.,1.]).unwrap();
            assert_eq!(p.position_m[0],1.); assert_eq!(p.position_m[1],0.);
            assert_eq!(p.weights,vec![1.,-1.]);
        }
        for bary in [[-0.1,0.5,0.6],[0.,0.,0.],[f64::NAN,0.,1.]] {
            assert!(r.surface_point_port(&h,0,bary,ShellFace::Positive,[0.,0.,1.]).is_err());
        }
        assert!(r.surface_point_port(&h,2,[1.,0.,0.],ShellFace::Positive,[0.,0.,1.]).is_err());
        assert!(r.surface_point_port(&[0.03;4],0,[1.,0.,0.],ShellFace::Positive,[0.,0.,1.]).is_err());
        assert!(r.surface_point_port(&h,0,[1.,0.,0.],ShellFace::Positive,[0.,0.,2.]).is_err());
    }

    #[test]
    fn finite_skin_closes_both_faces_and_preserves_opposite_velocities() {
        let s = square(vec![0.02; 4]).radiation_surface(&[0.02; 4], budget()).unwrap();
        assert_eq!(s.triangles().len(), 12);
        assert_eq!(s.boundary_panel_count(), 8);
        assert!((s.material_volume_m3()-0.02).abs() < 1e-14);
        assert_eq!(&s.weights[0][..4], &[1., -1., 1., -1.]);
        for row in &s.weights {
            let flux: f64 = row.iter().zip(&s.areas).map(|(v,a)| v*a).sum();
            assert!(flux.abs() < 1e-14, "rigid motion cannot create a monopole: {flux}");
        }
        // Offset rotation is load-bearing on the boundary walls.
        assert!(s.weights[1][4..].iter().any(|v| v.abs() > 0.001));
    }
    #[test]
    fn pressure_velocity_power_is_reciprocal_and_pressure_pushes_inward() {
        let h = vec![0.01, 0.02, 0.03, 0.015];
        let s = square(h.clone()).radiation_surface(&h, budget()).unwrap();
        let velocity = [0.37, -0.81];
        let pressure: Vec<_> = (0..s.areas.len()).map(|i| (i as f64 + 1.0)*17.0).collect();
        let mut vn = vec![0.0; s.areas.len()]; let mut force = [0.0; 2];
        s.normal_velocity(&velocity, &mut vn).unwrap();
        s.pressure_force(&pressure, &mut force).unwrap();
        let air: f64 = pressure.iter().zip(&s.areas).zip(&vn).map(|((p,a),v)| p*a*v).sum();
        let solid: f64 = velocity.iter().zip(force).map(|(v,f)| v*f).sum();
        assert!((air + solid).abs() < 1e-12);
        assert!(s.triangles[0][0][2] != s.triangles[0][1][2], "taper changes actual geometry");
        let mut p = vec![0.0; s.areas.len()]; p[0] = 1.0;
        s.pressure_force(&p, &mut force).unwrap();
        assert!(force[0] < 0.0);
    }
    #[test]
    fn mounting_hole_is_walled_not_capped() {
        let mut nodes = Vec::new();
        for r in [0.5, 1.0] { for p in [[-1.,-1.], [1.,-1.], [1.,1.], [-1.,1.]] {
            nodes.push([r*p[0], r*p[1], 0.0]);
        }}
        let mut tris = Vec::new();
        for i in 0..4 { let j = (i+1)%4; tris.push([i,i+4,j+4]); tris.push([i,j+4,j]); }
        let s = fixture(nodes, tris, vec![0.01; 8]).radiation_surface(&[0.01; 8], budget()).unwrap();
        assert_eq!(s.boundary_panel_count(), 16);
        assert_eq!(s.triangles.len(), 32);
        assert!((s.material_volume_m3() - 0.03).abs() < 1e-14);
        for row in &s.weights { assert!(row.iter().zip(&s.areas).map(|(b,a)| b*a).sum::<f64>().abs() < 1e-14); }
    }
    #[test]
    fn skin_rejects_inconsistent_thickness_winding_and_work_limits() {
        let s = square(vec![0.02; 4]);
        assert!(s.radiation_surface(&[0.03; 4], budget()).is_err());
        assert!(s.radiation_surface(&[f64::NAN; 4], budget()).is_err());
        assert!(s.radiation_surface(&[0.02; 4], RadiationSurfaceBudget { max_panels: 11, max_panel_modes: 1000 }).is_err());
        assert!(s.radiation_surface(&[0.02; 4], RadiationSurfaceBudget { max_panels: 12, max_panel_modes: 23 }).is_err());
        let wrong = fixture(s.reference_positions.clone(), vec![[0,1,2], [0,3,2]], vec![0.02; 4]);
        assert!(wrong.radiation_surface(&[0.02; 4], budget()).is_err());
        assert!(validate_link(&[(1,2), (3,4)]).is_err());
    }
    #[test]
    fn projection_refusal_is_transactional() {
        let s = square(vec![0.02; 4]).radiation_surface(&[0.02; 4], budget()).unwrap();
        let mut out = vec![42.0; 12];
        assert!(s.normal_velocity(&[f64::NAN, 1.0], &mut out).is_err());
        assert_eq!(out, vec![42.0; 12]);
        let mut force = [42.0; 2];
        assert!(s.pressure_force(&[f64::INFINITY; 12], &mut force).is_err());
        assert_eq!(force, [42.0; 2]);
    }
}
