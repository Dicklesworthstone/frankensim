//! Exact directional derivative of this shell's own metric-strain energy.
use super::{ShellReduction, dot, mul};

impl ShellReduction {
    /// Apply the potential Hessian to a direction in the original modal basis.
    /// Includes material AND geometric stiffness from the curved-shell strain;
    /// it is not the small-signal pencil or a positive-only approximation.
    /// No allocation. Invalid dimensions/nonfinite inputs fill `out` with NaN.
    pub fn hessian_vector(&self, q: &[f64], direction: &[f64], out: &mut [f64]) {
        let n = self.mode_count();
        if q.len()!=n || direction.len()!=n || out.len()!=n
            || q.iter().chain(direction).any(|v|!v.is_finite()) {
            out.fill(f64::NAN); return;
        }
        out.fill(0.0);
        if direction.iter().all(|d|*d==0.0) { return; }
        for i in 0..n {
            out[i]=(0..n).map(|j|self.remainder[i*n+j]*direction[j]).sum();
        }
        for f in &self.facets {
            let (strain,dx,dy)=Self::strain(f,q);
            let (mut ddx,mut ddy,mut ds)=([0.0;3],[0.0;3],[0.0;3]);
            for (&d,s) in direction.iter().zip(&f.strains) { for c in 0..3 {
                ddx[c]+=s.dx[c]*d; ddy[c]+=s.dy[c]*d; ds[c]+=s.linear[c]*d;
            }}
            ds[0]+=dot(dx,ddx); ds[1]+=dot(dy,ddy);
            ds[2]+=dot(ddx,dy)+dot(dx,ddy);
            let stress=mul(&f.membrane,strain);
            let dstress=mul(&f.membrane,ds);
            for (i,s) in f.strains.iter().enumerate() {
                let b=[s.linear[0]+dot(dx,s.dx),s.linear[1]+dot(dy,s.dy),
                    s.linear[2]+dot(dy,s.dx)+dot(dx,s.dy)];
                let db=[dot(ddx,s.dx),dot(ddy,s.dy),dot(ddy,s.dx)+dot(ddx,s.dy)];
                out[i]+=dot(b,dstress)+dot(db,stress);
            }
        }
    }
}
