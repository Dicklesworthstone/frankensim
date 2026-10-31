//! Work-conjugate point kinematics on six-DOF shell nodes.
//!
//! These are the P1 interface rows used by the existing shell point and skin
//! ports, not a reconstruction of an element's enriched interior rotation.
//! Director-aligned skin offsets use `u + theta cross r`; arbitrary rigid arms
//! use the physical orientation recovered from directors and membrane spin.
//! Their transposes carry force and moment on the SAME original nodes, without
//! treating a numerical drilling coordinate as physical material rotation.
//! Offsets are global Cartesian metres; physical angles are axial radians.
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
/// This is kinematic algebra, not admission of all three shell rotations as
/// physical. A shell node with a drilling gauge must use
/// [`nodal_director_point_row`] or a triangle's physical orientation.
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

/// Nodal offset force that is blind to the actual director-axis drilling gauge.
/// A zero or director-aligned offset is admitted. Other arms are admitted only
/// when their moment has no director-axis component, to geometric roundoff.
/// A node alone does not select a facet's membrane spin; an arbitrary offset
/// requiring that spin must use [`triangle_offset_point_rows`] instead.
///
/// # Errors
/// Invalid point/director or an offset moment along the numerical drilling axis.
pub fn nodal_director_point_row(direction: [f64;3], offset_m: [f64;3], director: [f64;3])
    -> Result<[f64;6], ShellPointError>
{
    let row = nodal_point_row(direction,offset_m)?;
    nodal_point_row(director,[0.0;3])?;
    let scale = row[3..].iter().map(|v|v.abs()).fold(0.0_f64,f64::max);
    if scale > 0.0 && (0..3).map(|i|(row[3+i]/scale)*director[i]).sum::<f64>().abs()>1e-12 {
        return Err(ShellPointError("shell node offset requires a physical triangle orientation, not director-axis drilling"));
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
/// For shells with numerical drilling, the offsets must follow the actual
/// nodal directors (or each moment must otherwise annihilate its nodal gauge).
/// A common arbitrary rigid arm instead uses [`triangle_offset_point_rows`].
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

/// Force/velocity rows for a rigid reference arm on a physical triangle point.
/// The point is `sum(w*x) + offset_m`; its velocity is `sum(w*v) + Omega cross
/// offset_m`, where `Omega` is exactly [`triangle_rotation_rows`]'s gauge-free
/// P1 orientation. The transpose adds the conjugate torque `offset_m cross F`
/// through that SAME owner, including its membrane forces. It never loads the
/// independent director-axis drilling coordinates.
///
/// Frame, gradients and directors belong to the original triangle. This does
/// not add arm inertia/compliance or reconstruct the enriched element interior.
///
/// # Errors
/// Invalid point geometry, invalid physical orientation needed by its torque,
/// or an unrepresentable force/moment row.
pub fn triangle_offset_point_rows(barycentric: [f64;3], direction: [f64;3], offset_m: [f64;3],
    frame: [[f64;3];3], gradient: [[f64;2];3], directors: [[f64;3];3])
    -> Result<[[f64;6];3], ShellPointError>
{
    let mut rows = triangle_point_rows(barycentric,direction,[[0.0;3];3])?;
    let moment = nodal_point_row(direction,offset_m)?;
    // Unit Cartesian axes avoid normalizing a zero/tiny/overflowing torque.
    // Linearity gives the same physical orientation map for every arm.
    for i in 0..3 {
        if moment[3+i] == 0.0 { continue; }
        let mut axis = [0.0;3]; axis[i] = 1.0;
        let torque = triangle_rotation_rows(barycentric,axis,frame,gradient,directors)?;
        for a in 0..3 { for dof in 0..6 { rows[a][dof] += moment[3+i]*torque[a][dof]; } }
    }
    if rows.iter().flatten().any(|v|!v.is_finite()) {
        return Err(ShellPointError("shell offset point pullback overflows"));
    }
    Ok(rows)
}

/// Gauge-free physical angle/torque rows at explicit triangle barycentrics.
/// `frame` and `gradient` are the ORIGINAL `ShellMesh::facet` outputs; directors
/// are the original unit nodal directors in the same triangle order.
///
/// Let d=sum(w*n), delta_d=sum(w*(theta cross n)). Director motion determines
/// t=(d cross delta_d)/|d|²; the original P1 membrane determines facet-normal
/// spin chi=0.5 curl(u).n_f. Physical rotation is
/// Omega=t+d*(chi-n_f.t)/(n_f.d). These rows observe axis.Omega and return its
/// conjugate torque. Every independent theta_a parallel to n_a is annihilated:
/// numerical drilling stabilization never supplies physical torsion stiffness.
///
/// In the flat-director limit, tangent rotation comes from director motion and
/// normal rotation from membrane spin. Nodal forces sum to zero and forces plus
/// moments carry the requested torque. This is the shared P1 interface motion,
/// not a reconstruction of the condensed enriched element interior.
///
/// # Errors
/// Refuses invalid interpolation/axis/frame/gradients, incompatible original
/// directors or unrepresentable orientation recovery. No input is repaired.
pub fn triangle_rotation_rows(barycentric: [f64; 3], axis: [f64; 3],
    frame: [[f64;3];3], gradient: [[f64;2];3], directors: [[f64;3];3])
    -> Result<[[f64;6];3], ShellPointError>
{
    triangle_point_rows(barycentric,axis,[[0.0;3];3])?;
    let dot = |a: [f64;3],b: [f64;3]| a.iter().zip(b).map(|(x,y)|x*y).sum::<f64>();
    let cross = |a: [f64;3],b: [f64;3]| [a[1]*b[2]-a[2]*b[1],
        a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
    for v in frame.iter().chain(&directors) { nodal_point_row(*v,[0.0;3])?; }
    if (0..3).any(|i|(0..i).any(|j|dot(frame[i],frame[j]).abs()>1e-12))
        || (dot(cross(frame[0],frame[1]),frame[2])-1.0).abs()>1e-12
        || directors.iter().any(|n|dot(*n,frame[2])<=0.0)
        || gradient.iter().flatten().any(|v|!v.is_finite())
    {
        return Err(ShellPointError("shell rotation needs original orthonormal facet geometry and aligned directors"));
    }
    let scale = gradient.iter().flatten().map(|v|v.abs()).fold(0.0_f64,f64::max);
    if scale == 0.0 || (0..2).any(|i|gradient.iter().map(|g|g[i]/scale).sum::<f64>().abs()>1e-12) {
        return Err(ShellPointError("shell rotation needs original partition-preserving facet gradients"));
    }
    let d: [f64;3] = core::array::from_fn(|i|(0..3).map(|a|barycentric[a]*directors[a][i]).sum());
    let square = dot(d,d); let alignment = dot(frame[2],d);
    if !square.is_finite() || square<=0.0 || !alignment.is_finite() || alignment<=0.0 {
        return Err(ShellPointError("shell rotation director interpolation is singular"));
    }
    let alpha = dot(axis,d)/alignment;
    let b = core::array::from_fn(|i|axis[i]-alpha*frame[2][i]);
    let c = cross(b,d).map(|v|v/square);
    let mut rows = [[0.0;6];3];
    for a in 0..3 {
        let moment = cross(directors[a],c);
        for i in 0..3 {
            // Same facet membrane-spin row as the original MITC owner.
            rows[a][i] = 0.5*alpha*(gradient[a][0]*frame[1][i]-gradient[a][1]*frame[0][i]);
            rows[a][3+i] = barycentric[a]*moment[i];
        }
    }
    if rows.iter().flatten().any(|v|!v.is_finite()) {
        return Err(ShellPointError("shell physical rotation pullback overflows"));
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

    #[test]
    fn physical_rotation_carries_a_pure_couple_with_rigid_work_and_no_drill_gauge() {
        let frame = [[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]];
        let gradient = [[-1.0,-1.0],[1.0,0.0],[0.0,1.0]];
        let nodes = [[0.0,0.0,0.0],[1.0,0.0,0.0],[0.0,1.0,0.0]];
        let barycentric = [0.2,0.3,0.5]; let rigid = [0.7,-0.2,0.3];
        for directors in [[[0.0,0.0,1.0];3],
            [[0.2,0.0,0.96_f64.sqrt()],[0.0,0.3,0.91_f64.sqrt()],[-0.1,-0.15,0.9675_f64.sqrt()]]]
        {
            for axis in [[1.0,0.0,0.0],[0.0,0.0,1.0],[0.6,0.0,0.8]] {
                let rows = triangle_rotation_rows(barycentric,axis,frame,gradient,directors).unwrap();
                let mut force = [0.0;3]; let mut moment = [0.0;3]; let mut power = 0.0;
                for a in 0..3 {
                    let f = [rows[a][0],rows[a][1],rows[a][2]];
                    let m = cross(nodes[a],f); let v = cross(rigid,nodes[a]);
                    // Each node has its OWN independent director-axis gauge.
                    near((0..3).map(|i|rows[a][3+i]*directors[a][i]).sum(),0.0);
                    for i in 0..3 {
                        force[i] += f[i]; moment[i] += m[i]+rows[a][3+i];
                        power += f[i]*(v[i]+0.19) + rows[a][3+i]
                            *(rigid[i]+(a as f64+1.0)*0.73*directors[a][i]);
                    }
                }
                for i in 0..3 { near(force[i],0.0); near(moment[i],axis[i]); }
                near(power,dot(axis,rigid));
            }
        }
    }

    #[test]
    fn physical_rotation_power_matches_deforming_director_and_membrane_motion() {
        let frame = [[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]];
        let gradient = [[-1.0,-1.0],[1.0,0.0],[0.0,1.0]];
        let directors = [[0.2,0.0,0.96_f64.sqrt()],[0.0,0.3,0.91_f64.sqrt()],
            [-0.1,-0.15,0.9675_f64.sqrt()]];
        let weights = [0.2,0.3,0.5]; let axis = [0.6,0.0,0.8];
        let velocity = [[0.2,-0.1,0.3],[0.7,0.4,-0.2],[-0.3,0.5,0.6]];
        let omega = [[0.1,-0.3,0.4],[0.7,0.2,-0.5],[-0.4,0.9,0.3]];
        let mut d = [0.0;3]; let mut delta_d = [0.0;3];
        for a in 0..3 {
            let motion = cross(omega[a],directors[a]);
            for i in 0..3 { d[i]+=weights[a]*directors[a][i]; delta_d[i]+=weights[a]*motion[i]; }
        }
        let tangent = cross(d,delta_d).map(|v|v/dot(d,d));
        let spin = 0.5*((velocity[1][1]-velocity[0][1])-(velocity[2][0]-velocity[0][0]));
        let alpha = (spin-tangent[2])/d[2];
        let physical = core::array::from_fn(|i|tangent[i]+alpha*d[i]);
        let rows = triangle_rotation_rows(weights,axis,frame,gradient,directors).unwrap();
        let torque = 0.031;
        let power: f64 = (0..3).map(|a|(0..3).map(|i|torque*
            (rows[a][i]*velocity[a][i]+rows[a][3+i]*omega[a][i])).sum::<f64>()).sum();
        near(power,torque*dot(axis,physical));
    }

    #[test]
    fn arbitrary_offset_force_uses_physical_orientation_work_and_annihilates_drilling() {
        let frame = [[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]];
        let gradient = [[-1.0,-1.0],[1.0,0.0],[0.0,1.0]];
        let nodes = [[0.0,0.0,0.0],[1.0,0.0,0.0],[0.0,1.0,0.0]];
        let weights = [0.2,0.3,0.5];
        let velocity = [[0.2,-0.1,0.3],[0.7,0.4,-0.2],[-0.3,0.5,0.6]];
        let omega = [[0.1,-0.3,0.4],[0.7,0.2,-0.5],[-0.4,0.9,0.3]];
        for directors in [[[0.0,0.0,1.0];3],
            [[0.2,0.0,0.96_f64.sqrt()],[0.0,0.3,0.91_f64.sqrt()],[-0.1,-0.15,0.9675_f64.sqrt()]]]
        {
            for (direction,offset) in [([0.0,1.0,0.0],[0.03,0.0,0.0]),
                ([0.6,0.0,0.8],[0.02,-0.03,0.005])]
            {
                let rows = triangle_offset_point_rows(weights,direction,offset,frame,gradient,directors).unwrap();
                let point = core::array::from_fn(|i|offset[i]+(0..3).map(|a|weights[a]*nodes[a][i]).sum::<f64>());
                let mut resultant = [0.0;3]; let mut moment = [0.0;3]; let mut power = 0.0;
                let mut d = [0.0;3]; let mut delta_d = [0.0;3]; let mut mean_v = [0.0;3];
                for a in 0..3 {
                    let f = [rows[a][0],rows[a][1],rows[a][2]];
                    let m = cross(nodes[a],f); let motion = cross(omega[a],directors[a]);
                    near((0..3).map(|i|rows[a][3+i]*directors[a][i]).sum(),0.0);
                    for i in 0..3 {
                        resultant[i] += f[i]; moment[i] += m[i]+rows[a][3+i];
                        power += f[i]*velocity[a][i]+rows[a][3+i]
                            *(omega[a][i]+0.73*(a as f64+1.0)*directors[a][i]);
                        d[i] += weights[a]*directors[a][i];
                        delta_d[i] += weights[a]*motion[i]; mean_v[i] += weights[a]*velocity[a][i];
                    }
                }
                let tangent = cross(d,delta_d).map(|v|v/dot(d,d));
                let spin = 0.5*((velocity[1][1]-velocity[0][1])-(velocity[2][0]-velocity[0][0]));
                let physical = core::array::from_fn(|i|tangent[i]+d[i]*(spin-tangent[2])/d[2]);
                let arm_v = cross(physical,offset);
                let expected_moment = cross(point,direction);
                for i in 0..3 { near(resultant[i],direction[i]); near(moment[i],expected_moment[i]); }
                near(power,dot(direction,core::array::from_fn(|i|mean_v[i]+arm_v[i])));
            }
        }
    }

    #[test]
    fn node_only_offsets_admit_bending_and_refuse_ambiguous_drilling_moments() {
        for director in [[0.0,0.0,1.0],[0.6,0.0,0.8]] {
            let direction = [0.0,1.0,0.0];
            for offset in [[0.0;3],director.map(|v|0.02*v),[0.0,0.03,0.0]] {
                let row = nodal_director_point_row(direction,offset,director).unwrap();
                assert_eq!(row,nodal_point_row(direction,offset).unwrap());
                near((0..3).map(|i|row[3+i]*director[i]).sum(),0.0);
            }
            for offset in [[0.03,0.0,0.0],[-0.03,0.0,0.0]] {
                assert!(nodal_director_point_row(direction,offset,director).is_err());
            }
        }
    }

    #[test]
    fn torque_ports_refuse_invalid_axes_and_interpolation_without_repair() {
        let frame = [[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]];
        let gradient = [[-1.0,-1.0],[1.0,0.0],[0.0,1.0]];
        let directors = [[0.0,0.0,1.0];3];
        for axis in [[0.0;3],[0.0,0.0,2.0],[f64::NAN,0.0,1.0]] {
            assert!(triangle_rotation_rows([0.2,0.3,0.5],axis,frame,gradient,directors).is_err());
        }
        for barycentric in [[-0.1,0.5,0.6],[0.2,0.3,0.4],[f64::NAN,0.0,1.0]] {
            assert!(triangle_rotation_rows(barycentric,[0.0,1.0,0.0],frame,gradient,directors).is_err());
        }
        assert!(triangle_rotation_rows([0.2,0.3,0.5],[0.0,1.0,0.0],frame,gradient,
            [[0.0,0.0,-1.0];3]).is_err());
        assert!(triangle_rotation_rows([0.2,0.3,0.5],[0.0,1.0,0.0],frame,[[0.0;2];3],directors).is_err());
    }
}
