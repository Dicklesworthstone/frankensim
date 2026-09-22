//! Analytic tangent of the existing distributed unilateral potential.
use super::ContactStorage;
use fs_math::det;

impl ContactStorage {
    /// Compose an inner-storage Hessian action with these exact contact rows.
    /// The callback must describe `inner_storage()` with the same frozen state;
    /// momentum and internal coordinates need not be unit-mass or quadratic.
    /// Contact stiffness acts only on the interleaved displacement prefix.
    /// At zero penetration the inactive branch is selected, as in `gradient`.
    /// This differentiates elastic storage, not the Hunt--Crossley loss port.
    ///
    /// No allocation. Returns false for incompatible/nonfinite data or an
    /// unavailable inner action; `out` is trial scratch and may then be partial.
    pub fn hessian_vector_with<F>(&self, x: &[f64], direction: &[f64],
        out: &mut [f64], inner_action: F) -> bool
    where F: FnOnce(&[f64], &[f64], &mut [f64]) -> bool {
        if self.n_modes.checked_mul(2).is_none_or(|n|x.len()<n)
            || x.len()!=direction.len() || out.len()!=x.len()
            || x.iter().chain(direction).any(|v|!v.is_finite()) { return false; }
        if !inner_action(x,direction,out) { return false; }
        if (0..self.n_modes).any(|k|direction[2*k]!=0.0) {
            for ob in &self.obstacles { for i in 0..ob.n_points {
                let p=ob.penetration_at(self.n_modes,x,i);
                if p<=0.0 { continue; }
                let row=&ob.collocation[i*self.n_modes..(i+1)*self.n_modes];
                let rate=row.iter().enumerate().map(|(k,b)|b*direction[2*k]).sum::<f64>();
                let df=ob.weights[i]*ob.stiffness*ob.alpha*det::pow(p,ob.alpha-1.0)*rate;
                for (k,b) in row.iter().enumerate() { out[2*k]+=df*b; }
            }}
        }
        out.iter().all(|v|v.is_finite())
    }
}
