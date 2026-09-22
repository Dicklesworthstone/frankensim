//! Tangent of the same statically relaxed, positive facet stretching energy.
use super::{MembraneReduction, dot, mul};

impl MembraneReduction {
    /// Apply the exact reduced potential Hessian without a quartic tensor or
    /// heap scratch. Installed tension/bending is included exactly once; the
    /// condensed in-plane strains retain all mixed transverse-mode terms.
    /// Invalid dimensions/nonfinite inputs fill `out` with NaN.
    pub fn hessian_vector(&self, q: &[f64], direction: &[f64], out: &mut [f64]) {
        let n=self.mode_count();
        if q.len()!=n || direction.len()!=n || out.len()!=n
            || q.iter().chain(direction).any(|v|!v.is_finite()) {
            out.fill(f64::NAN); return;
        }
        out.fill(0.0);
        if direction.iter().all(|d|*d==0.0) { return; }
        for i in 0..n { out[i]=(0..n).map(|j|self.linear[i*n+j]*direction[j]).sum(); }
        for f in &self.facets {
            let stress=mul(&f.a,self.strain(f,q));
            let mut ds=[0.0;3];
            for (&(i,j),s) in self.pairs.iter().zip(&f.strains) { for c in 0..3 {
                ds[c]+=(direction[i]*q[j]+q[i]*direction[j])*s[c];
            }}
            let dstress=mul(&f.a,ds);
            for (&(i,j),s) in self.pairs.iter().zip(&f.strains) {
                let force=dot(*s,stress); let df=dot(*s,dstress);
                out[i]+=direction[j]*force+q[j]*df;
                out[j]+=direction[i]*force+q[i]*df;
            }
        }
    }
}
