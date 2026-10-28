//! Exact body/right linearization of the existing geodesic-midpoint rod.
//!
//! Kinematics, section law, global equilibrium and LU remain their original
//! owners. This module differentiates the same strain map without probing a
//! displaced/material state. With Q=Ra^T Rb, phi=Log(Q), H=Exp(phi/2), M=Ra H:
//! dphi=Ba theta_a+Bb theta_b; eta=Ta theta_a+Tb theta_b;
//! ds=M^T d(rb-ra)/L+s cross eta. Virtual work supplies both nodal torques.
//! Differentiating these BODY torques includes their changing basis. The result
//! is the residual Jacobian, not a forcibly symmetrized coordinate Hessian.
//!
//! SO(3) Exp/Log and all Jacobian differentials come from fs-ga. The principal
//! logarithm cut at pi is explicitly refused for derivatives, not differentiated
//! through a quaternion-sign switch. Refine the rotation field across that cut.

use super::super::{Rod, SolidError, group_error};
use fs_ga::{Mat3, So3, So3Tangent, Vec3, so3_left_jacobian_inverse,
    so3_right_jacobian, so3_right_jacobian_inverse};
use fs_ga::lie_differential::{so3_left_jacobian_inverse_directional,
    so3_right_jacobian_directional, so3_right_jacobian_inverse_directional};

pub(in crate::rod) fn add(a: Mat3, b: Mat3) -> Mat3 {
    Mat3 { m: std::array::from_fn(|i| a.m[i]+b.m[i]) }
}
pub(in crate::rod) fn scale(a: Mat3, s: f64) -> Mat3 { Mat3 { m: a.m.map(|v| s*v) } }
pub(in crate::rod) fn hat(v: Vec3) -> Mat3 {
    Mat3 { m: [0.0,-v.z,v.y,v.z,0.0,-v.x,-v.y,v.x,0.0] }
}
pub(in crate::rod) fn vadd(a: Vec3, b: Vec3) -> Vec3 { Vec3::new(a.x+b.x,a.y+b.y,a.z+b.z) }
pub(in crate::rod) fn vsub(a: Vec3, b: Vec3) -> Vec3 { Vec3::new(a.x-b.x,a.y-b.y,a.z-b.z) }
pub(in crate::rod) fn array(v: Vec3) -> [f64;3] { [v.x,v.y,v.z] }
pub(in crate::rod) fn vector(v: [f64;3]) -> Vec3 { Vec3::new(v[0],v[1],v[2]) }
fn diagonal(d: [f64;3], v: Vec3) -> Vec3 { Vec3::new(d[0]*v.x,d[1]*v.y,d[2]*v.z) }
fn invalid(what: &str) -> SolidError { SolidError::InvalidInput { what: what.into() } }

/// One native element's admitted configuration, not a second pose representation.
pub(in crate::rod) struct Kinematics {
    length: f64,
    relative: So3Tangent,
    half: So3,
    midpoint: So3,
    local: Vec3,
}
impl Kinematics {
    pub(in crate::rod) fn new(rod: &Rod, segment: usize) -> Result<Self, SolidError> {
        if rod.positions.len() < 2 || rod.frames.len() != rod.positions.len()
            || segment >= rod.positions.len()-1 || !rod.l0.is_finite() || rod.l0 <= 0.0 {
            return Err(invalid("rod kinematics require matching nodes/frames, a valid segment and positive finite reference length"));
        }
        let relative = rod.frames[segment+1].body_minus(rod.frames[segment]).map_err(group_error)?;
        let half = So3::exp(relative.scale(0.5)).map_err(group_error)?;
        let midpoint = rod.frames[segment].compose(half).map_err(group_error)?;
        let delta = vsub(vector(rod.positions[segment+1]),vector(rod.positions[segment]));
        let d = Vec3::new(delta.x/rod.l0,delta.y/rod.l0,delta.z/rod.l0);
        let local = midpoint.inverse().rotate(d).map_err(group_error)?;
        if array(local).iter().chain(&array(relative.angular)).any(|v| !v.is_finite()) {
            return Err(invalid("rod strain kinematics are not representable"));
        }
        Ok(Self { length: rod.l0, relative, half, midpoint, local })
    }

    pub(in crate::rod) fn strains(&self) -> ([f64;3],[f64;3]) {
        ([self.local.x-1.0,self.local.y,self.local.z], [self.relative.angular.x/self.length,self.relative.angular.y/self.length,self.relative.angular.z/self.length])
    }

    /// n and m are current section resultants; the two diagonal arrays are
    /// dn/dgamma and dm/dkappa at this same material state. Eigenstrains change
    /// resultants, not this kinematic pullback. Units remain N, N m and metres.
    pub(in crate::rod) fn linearize(
        &self, n: Vec3, m: Vec3, force_moduli: [f64;3], moment_moduli: [f64;3],
        tangent: bool,
    ) -> Result<PairResponse, SolidError> {
        let angle = self.relative.angular.norm();
        if !angle.is_finite() || angle >= core::f64::consts::PI*(1.0-64.0*f64::EPSILON) {
            return Err(invalid("rod relative rotation is on the SO(3) logarithm cut; refine segments before differentiating"));
        }
        let b = [scale(so3_left_jacobian_inverse(self.relative).map_err(group_error)?, -1.0),
            so3_right_jacobian_inverse(self.relative).map_err(group_error)?];
        let j = so3_right_jacobian(self.relative.scale(0.5)).map_err(group_error)?;
        let half_t = self.half.matrix().transpose();
        let tmap = [add(half_t,scale(j.compose(b[0]),0.5)),scale(j.compose(b[1]),0.5)];
        let rotation = self.midpoint.matrix();
        let t = n.cross(self.local).scale(self.length);
        let world = rotation.apply(n);
        let moments: [Vec3;2] = std::array::from_fn(|i|
            vadd(tmap[i].transpose().apply(t), b[i].transpose().apply(m)));
        let gradient = [join(world.scale(-1.0),moments[0]),join(world,moments[1])];
        let matrix = if tangent {
            let mut matrix = [[0.0;12];12];
            for column in 0..12 {
                let side = column/6;
                let component = column%6;
                let mut direction = [0.0;3];
                direction[component%3] = 1.0;
                let direction = vector(direction);
                let (dp,eta,dr) = if component < 3 {
                    (Vec3::new(0.0,0.0,0.0),Vec3::new(0.0,0.0,0.0),
                        direction.scale(if side==0 {-1.0} else {1.0}))
                } else {
                    (b[side].apply(direction),tmap[side].apply(direction),Vec3::new(0.0,0.0,0.0))
                };
                let ds = vadd(rotation.transpose().apply(dr).scale(1.0/self.length),self.local.cross(eta));
                let dn = diagonal(force_moduli,ds);
                let dm = diagonal(moment_moduli,dp.scale(1.0/self.length));
                let dt = vadd(dn.cross(self.local),n.cross(ds)).scale(self.length);
                let dworld = rotation.apply(vadd(eta.cross(n),dn));
                let (db,dtmap) = if component < 3 {
                    ([Mat3::zero();2],[Mat3::zero();2])
                } else {
                    let differential = So3Tangent::new(dp);
                    let db = [scale(so3_left_jacobian_inverse_directional(self.relative,differential)
                        .map_err(group_error)?, -1.0),
                        so3_right_jacobian_inverse_directional(self.relative,differential).map_err(group_error)?];
                    let dj = so3_right_jacobian_directional(self.relative.scale(0.5),differential.scale(0.5))
                        .map_err(group_error)?;
                    let dhalf_t = scale(hat(j.apply(dp.scale(0.5))).compose(half_t),-1.0);
                    let dtmap = [add(dhalf_t,scale(add(dj.compose(b[0]),j.compose(db[0])),0.5)),
                        scale(add(dj.compose(b[1]),j.compose(db[1])),0.5)];
                    (db,dtmap)
                };
                for i in 0..2 {
                    let moment = vadd(vadd(dtmap[i].transpose().apply(t),tmap[i].transpose().apply(dt)),
                        vadd(db[i].transpose().apply(m),b[i].transpose().apply(dm)));
                    let column_values = join(dworld.scale(if i==0 {-1.0} else {1.0}),moment);
                    for row in 0..6 { matrix[6*i+row][column]=column_values[row]; }
                }
            }
            Some(matrix)
        } else { None };
        let result = PairResponse { gradient, tangent: matrix };
        result.check()?;
        Ok(result)
    }
}

pub(in crate::rod) fn join(force: Vec3, moment: Vec3) -> [f64;6] {
    [force.x,force.y,force.z,moment.x,moment.y,moment.z]
}
pub(in crate::rod) struct PairResponse {
    pub gradient: [[f64;6];2],
    pub tangent: Option<[[f64;12];12]>,
}
impl PairResponse {
    pub(in crate::rod) fn check(&self) -> Result<(), SolidError> {
        if self.gradient.iter().flatten().any(|v| !v.is_finite())
            || self.tangent.as_ref().is_some_and(|k| k.iter().flatten().any(|v| !v.is_finite())) {
            return Err(invalid("analytical rod element derivative is not representable"));
        }
        Ok(())
    }
}

/// Full nodal gradient, including root, plus an optional reduced dense Jacobian.
/// The original solver still owns factorization and the free-coordinate policy.
pub(in crate::rod) struct Linearization {
    pub gradient: Vec<[f64;6]>,
    pub tangent: Option<Vec<f64>>,
}
pub(in crate::rod) struct Assembly {
    result: Linearization,
    row_of: Vec<usize>,
    dimension: usize,
}
impl Assembly {
    pub(in crate::rod) fn new(nodes: usize, free: &[usize], tangent: bool) -> Result<Self,SolidError> {
        let size = nodes.checked_mul(6).ok_or_else(|| invalid("rod derivative size overflow"))?;
        let matrix_size = free.len().checked_mul(free.len()).ok_or_else(|| invalid("rod tangent size overflow"))?;
        let mut row_of = vec![usize::MAX;size];
        for (row,&coordinate) in free.iter().enumerate() {
            let index = coordinate.checked_add(6).ok_or_else(|| invalid("rod coordinate overflow"))?;
            if index >= size || row_of[index] != usize::MAX { return Err(invalid("invalid free rod coordinate map")); }
            row_of[index]=row;
        }
        Ok(Self { result: Linearization { gradient: vec![[0.0;6];nodes],
            tangent: tangent.then(|| vec![0.0;matrix_size]) }, row_of, dimension:free.len() })
    }
    pub(in crate::rod) fn add_pair(&mut self,nodes:[usize;2], response:&PairResponse) {
        for (side,&node) in nodes.iter().enumerate() {
            for component in 0..6 { self.result.gradient[node][component]+=response.gradient[side][component]; }
        }
        if let (Some(global),Some(local))=(&mut self.result.tangent,&response.tangent) {
            for row in 0..12 {
                let a=self.row_of[6*nodes[row/6]+row%6];
                if a==usize::MAX {continue;}
                for col in 0..12 {
                    let b=self.row_of[6*nodes[col/6]+col%6];
                    if b!=usize::MAX {global[a*self.dimension+b]+=local[row][col];}
                }
            }
        }
    }
    pub(in crate::rod) fn add_node(&mut self,node:usize,gradient:[f64;6],tangent:[[f64;6];6]) {
        for c in 0..6 {self.result.gradient[node][c]+=gradient[c];}
        if let Some(global)=&mut self.result.tangent {
            for row in 0..6 {
                let a=self.row_of[6*node+row];
                if a==usize::MAX {continue;}
                for col in 0..6 {
                    let b=self.row_of[6*node+col];
                    if b!=usize::MAX {global[a*self.dimension+b]+=tangent[row][col];}
                }
            }
        }
    }
    pub(in crate::rod) fn finish(self)->Result<Linearization,SolidError> {
        if self.result.gradient.iter().flatten().any(|v| !v.is_finite())
            || self.result.tangent.as_ref().is_some_and(|k| k.iter().any(|v| !v.is_finite())) {
            return Err(invalid("assembled analytical rod derivative is not representable"));
        }
        Ok(self.result)
    }
}

#[cfg(test)]
mod tests;
