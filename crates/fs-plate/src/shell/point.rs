//! Work-conjugate point kinematics on six-DOF shell nodes.
//!
//! These are the P1 interface rows used by the existing shell point and skin
//! ports, not a reconstruction of an element's enriched interior rotation.
//! For a reference offset `r`, material motion is `u + theta cross r`.
//! Its transpose carries both the force and its moment on the SAME nodes.
//! Offsets are global Cartesian metres; rotations are physical axial radians.
//! No contact search, finite-rotation update or mesh-node snapping occurs here.

/// Invalid physical direction, interpolation weights or reference offset.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ShellPointError(&'static str);

impl ShellPointError {
    /// Stable explanation suitable for an enclosing mechanical owner's error.
    #[must_use]
    pub const fn reason(self) -> &'static str { self.0 }
}
impl core::fmt::Display for ShellPointError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-SHELL-POINT: {}", self.0)
    }
}
impl std::error::Error for ShellPointError {}

/// Unit-force row `[d, r cross d]` on `(u,v,w,theta_x,theta_y,theta_z)`.
/// Dot with nodal motion to observe the displacement/velocity at the offset
/// point; multiply by signed force to return its conjugate force and moment.
///
/// # Errors
/// Refuses a nonunit/nonfinite direction, nonfinite offset or moment overflow.
pub fn nodal_point_row(direction: [f64; 3], offset_m: [f64; 3])
    -> Result<[f64; 6], ShellPointError>
{
    let norm2 = direction.iter().map(|v| v * v).sum::<f64>();
    if !norm2.is_finite() || (norm2 - 1.0).abs() > 1e-12
        || offset_m.iter().any(|v| !v.is_finite())
    {
        return Err(ShellPointError("shell point requires a finite unit direction and finite offset"));
    }
    let [x, y, z] = direction;
    let [rx, ry, rz] = offset_m;
    let row = [x, y, z, ry * z - rz * y, rz * x - rx * z, rx * y - ry * x];
    if row.iter().any(|v| !v.is_finite()) {
        return Err(ShellPointError("shell point moment overflows"));
    }
    Ok(row)
}

/// Three original-node rows for one explicit triangle material point.
/// Nonnegative barycentrics sum to one; zero weights are admitted on edges and
/// vertices. The caller retains the original triangle's node order and admits
/// its mesh geometry. Each vertex may have its own reference offset, as on a
/// curved finite-thickness skin. A common offset represents a point displaced
/// from the interpolated midsurface location.
///
/// The point position is `sum_a weight[a] * (node[a] + offset[a])` and its
/// velocity is `sum_a weight[a] * (v[a] + omega[a] cross offset[a])`.
/// These SAME rows preserve total force, moment and rigid-motion virtual work.
///
/// # Errors
/// Refuses invalid barycentrics or any invalid nodal point row. It never clamps
/// weights, normalizes a direction or renormalizes a footprint.
pub fn triangle_point_rows(barycentric: [f64; 3], direction: [f64; 3],
    offsets_m: [[f64; 3]; 3]) -> Result<[[f64; 6]; 3], ShellPointError>
{
    if barycentric.iter().any(|v| !v.is_finite() || *v < 0.0)
        || (barycentric.iter().sum::<f64>() - 1.0).abs() > 1e-12
    {
        return Err(ShellPointError("shell point requires nonnegative unit-sum barycentrics"));
    }
    let mut rows = [[0.0; 6]; 3];
    for a in 0..3 {
        rows[a] = nodal_point_row(direction, offsets_m[a])?.map(|v| barycentric[a] * v);
        if rows[a].iter().any(|v| !v.is_finite()) {
            return Err(ShellPointError("shell point weighted row overflows"));
        }
    }
    Ok(rows)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
        [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
    }
    fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
        a.iter().zip(b).map(|(x, y)| x*y).sum()
    }
    fn near(a: f64, b: f64) {
        assert!((a-b).abs() < 2e-13*a.abs().max(b.abs()).max(1.0), "{a} != {b}");
    }

    #[test]
    fn offset_force_has_physical_moment_and_work_conjugate_rotation() {
        assert_eq!(nodal_point_row([1.0,0.0,0.0], [0.0,0.0,0.02]).unwrap(),
            [1.0,0.0,0.0,0.0,0.02,0.0]);
        let direction = [0.6,0.0,0.8];
        let offset = [0.02,-0.03,0.005];
        let translation = [0.4,-0.2,0.1];
        let rotation = [-0.1,0.7,0.3];
        let row = nodal_point_row(direction, offset).unwrap();
        let v = cross(rotation, offset);
        let point = core::array::from_fn(|i| translation[i]+v[i]);
        near(row[..3].iter().zip(translation).map(|(a,b)|a*b).sum::<f64>()
            + row[3..].iter().zip(rotation).map(|(a,b)|a*b).sum::<f64>(), dot(direction,point));
    }

    #[test]
    fn curved_skin_preserves_resultants_rigid_work_and_arbitrary_nodal_power() {
        let nodes = [[0.1,-0.2,0.02],[0.4,0.1,0.07],[-0.2,0.3,-0.01]];
        let offsets = [[0.001,0.002,0.003],[-0.002,0.001,0.004],[0.0,-0.001,0.002]];
        let weights = [0.2,0.3,0.5];
        let direction = [0.6,0.0,0.8];
        let rows = triangle_point_rows(weights,direction,offsets).unwrap();
        let point: [f64;3] = core::array::from_fn(|i| (0..3)
            .map(|a|weights[a]*(nodes[a][i]+offsets[a][i])).sum());
        let force_n = 7.3;
        let expected_force = direction.map(|v|force_n*v);
        let expected_moment = cross(point,expected_force);
        let mut force = [0.0;3]; let mut moment = [0.0;3];
        let rigid_v = [0.3,-0.1,0.2]; let rigid_omega = [-0.8,0.4,0.1];
        let mut rigid_power = 0.0; let mut nodal_power = 0.0; let mut point_v = [0.0;3];
        for a in 0..3 {
            let f = core::array::from_fn(|i|force_n*rows[a][i]);
            let m = cross(nodes[a],f);
            let rv = cross(rigid_omega,nodes[a]);
            let v = [0.3*(a as f64+1.0),-0.7,0.2];
            let omega = [0.1,-0.2*a as f64,0.6];
            let arm_v = cross(omega,offsets[a]);
            for i in 0..3 {
                force[i] += f[i];
                moment[i] += m[i]+force_n*rows[a][3+i];
                rigid_power += force_n*(rows[a][i]*(rigid_v[i]+rv[i])
                    + rows[a][3+i]*rigid_omega[i]);
                nodal_power += force_n*(rows[a][i]*v[i]+rows[a][3+i]*omega[i]);
                point_v[i] += weights[a]*(v[i]+arm_v[i]);
            }
        }
        for i in 0..3 { near(force[i],expected_force[i]); near(moment[i],expected_moment[i]); }
        near(rigid_power,dot(expected_force,rigid_v)+dot(expected_moment,rigid_omega));
        near(nodal_power,dot(expected_force,point_v));
    }

    #[test]
    fn same_point_in_refined_triangle_keeps_affine_motion_without_snapping() {
        let coarse = [[0.0,0.0,0.0],[1.0,0.0,0.0],[0.0,1.0,0.0]];
        let fine = [[0.0,0.0,0.0],[0.5,0.0,0.0],[0.0,0.5,0.0]];
        let direction = [0.0,0.0,1.0]; let offset = [0.02,0.0,0.01];
        let sample = |nodes: [[f64;3];3], weights| {
            let rows = triangle_point_rows(weights,direction,[offset;3]).unwrap();
            (0..3).map(|a| {
                let [x,y,z] = nodes[a];
                let state = [x-y,2.0*y,z+3.0*x-y,0.2+y,0.3-x,0.1];
                rows[a].iter().zip(state).map(|(b,q)|b*q).sum::<f64>()
            }).sum::<f64>()
        };
        near(sample(coarse,[0.6,0.25,0.15]),sample(fine,[0.2,0.5,0.3]));
    }

    #[test]
    fn invalid_footprints_refuse_instead_of_normalizing_or_clipping() {
        for weights in [[-0.1,0.5,0.6],[0.2,0.3,0.4],[f64::NAN,0.0,1.0],
            [f64::MAX,f64::MAX,0.0]] {
            assert!(triangle_point_rows(weights,[0.0,0.0,1.0],[[0.0;3];3]).is_err());
        }
        for direction in [[0.0;3],[0.0,0.0,2.0],[f64::NAN,0.0,1.0],[f64::MAX,0.0,0.0]] {
            assert!(nodal_point_row(direction,[0.0;3]).is_err());
        }
        assert!(nodal_point_row([1.0,0.0,0.0],[0.0,f64::INFINITY,0.0]).is_err());
        let q = core::f64::consts::FRAC_1_SQRT_2;
        assert!(nodal_point_row([q,q,0.0],[f64::MAX,-f64::MAX,0.0]).is_err());
        let vertex = triangle_point_rows([0.0,1.0,0.0],[1.0,0.0,0.0],[[0.0;3];3]).unwrap();
        assert_eq!(vertex,[[0.0;6],[1.0,0.0,0.0,0.0,0.0,0.0],[0.0;6]]);
    }
}
