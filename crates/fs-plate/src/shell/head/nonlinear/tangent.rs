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

    /// Assemble all columns of the SAME potential Hessian at one displacement.
    ///
    /// `out` is row-major n*n; `strain_gradients` is caller-owned n-by-3 scratch.
    /// Each facet's strain/stress and modal strain gradients are evaluated once,
    /// rather than once per direction. Both B^T A B and stress.d2(strain) remain:
    /// this is not a frozen linear pencil or a positive-only tangent. The cold
    /// reduction's mixed modal-pair coefficients are used without truncation.
    ///
    /// Work is O(facets*n^2), versus O(facets*n^3) for n separate actions of this
    /// condensed film. No heap allocation or dense fourth-order tensor is added.
    /// This operation count is not a measured end-to-end speedup. The existing
    /// action remains useful for individual directions. Rounding may differ.
    /// Invalid dimensions/nonfinite input fills `out` with NaN. Scratch and
    /// matrix are numerical trial data; finite-output admission belongs to caller.
    pub fn hessian_matrix(&self, q: &[f64], out: &mut [f64], strain_gradients: &mut [[f64;3]]) {
        let n=self.mode_count();
        if q.len()!=n || n.checked_mul(n)!=Some(out.len()) || strain_gradients.len()!=n
            || q.iter().any(|v|!v.is_finite()) {
            out.fill(f64::NAN); return;
        }
        out.copy_from_slice(&self.linear);
        for f in &self.facets {
            let stress=mul(&f.a,self.strain(f,q));
            strain_gradients.fill([0.0;3]);
            for (&(i,j),s) in self.pairs.iter().zip(&f.strains) {
                for c in 0..3 {
                    strain_gradients[i][c]+=q[j]*s[c];
                    strain_gradients[j][c]+=q[i]*s[c];
                }
                // Diagonal pairs contribute twice; off-diagonal pairs once to
                // each symmetric entry, exactly as differentiating the gradient.
                let geometric=dot(*s,stress);
                out[i*n+j]+=geometric;
                out[j*n+i]+=geometric;
            }
            for i in 0..n {
                let stress_direction=mul(&f.a,strain_gradients[i]);
                for j in 0..=i {
                    let material=dot(strain_gradients[j],stress_direction);
                    out[i*n+j]+=material;
                    if i!=j {out[j*n+i]+=material;}
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::Facet;

    fn specimen(n:usize)->MembraneReduction {
        let pairs:Vec<_>=(0..n).flat_map(|i|(i..n).map(move|j|(i,j))).collect();
        let facets=(0..3).map(|f|Facet {
            a:[4.0,0.4,0.0, 0.4,3.0,0.0, 0.0,0.0,1.2],
            gradients:vec![[0.0;2];n],
            strains:pairs.iter().map(|&(i,j)| {
                let t=(1+i+3*j+7*f) as f64;
                [0.2*t.sin(),0.3*(t*0.7).cos(),0.1*(t*1.3).sin()]
            }).collect(),
        }).collect();
        let mut linear=vec![0.0;n*n];
        for i in 0..n {linear[i*n+i]=1.0+i as f64;}
        MembraneReduction {omegas:vec![1.0;n],linear,pairs,facets,solve_residual:0.0}
    }

    #[test]
    fn g0_batched_film_tangent_retains_every_mixed_pair_and_matches_gradient_derivatives() {
        for n in [1,3,8] {
            let law=specimen(n);let mut matrix=vec![0.0;n*n];let mut scratch=vec![[0.0;3];n];
            for scale in [0.0,0.01,0.4] {
                let q:Vec<_>=(0..n).map(|i|scale*(i as f64+0.3).sin()).collect();
                law.hessian_matrix(&q,&mut matrix,&mut scratch);
                for col in 0..n {
                    let mut d=vec![0.0;n];d[col]=1.0;let mut product=vec![0.0;n];
                    law.hessian_vector(&q,&d,&mut product);
                    let mut plus=q.clone();let mut minus=q.clone();plus[col]+=1e-6;minus[col]-=1e-6;
                    let mut gp=vec![0.0;n];let mut gm=vec![0.0;n];
                    law.gradient(&plus,&mut gp);law.gradient(&minus,&mut gm);
                    for row in 0..n {
                        let actual=matrix[row*n+col];let size=1.0+actual.abs();
                        assert!((actual-product[row]).abs()<2e-13*size);
                        assert!((actual-(gp[row]-gm[row])/2e-6).abs()<2e-8*size);
                        assert_eq!(actual.to_bits(),matrix[col*n+row].to_bits());
                    }
                }
            }
        }
    }

    #[test]
    fn g0_film_matrix_reuses_scratch_and_refuses_bad_shapes_without_panics() {
        let law=specimen(3);let mut out=[0.0;9];let mut scratch=[[0.0;3];3];
        law.hessian_matrix(&[0.1,-0.2,0.3],&mut out,&mut scratch);let expected=out;
        out.fill(123.0);scratch.fill([456.0;3]);
        law.hessian_matrix(&[0.1,-0.2,0.3],&mut out,&mut scratch);assert_eq!(out,expected);
        law.hessian_matrix(&[f64::NAN,0.0,0.0],&mut out,&mut scratch);assert!(out.iter().all(|v|v.is_nan()));
        law.hessian_matrix(&[0.0;3],&mut out,&mut []);assert!(out.iter().all(|v|v.is_nan()));
        law.hessian_matrix(&[0.0;3],&mut out[..8],&mut scratch);assert!(out[..8].iter().all(|v|v.is_nan()));
        law.hessian_matrix(&[0.0;3],&mut out,&mut scratch);assert_eq!(out,law.linear.as_slice());
    }
}
