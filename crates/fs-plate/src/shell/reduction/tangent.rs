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

    /// Assemble the exact row-major potential Hessian using caller-owned n-by-3
    /// strain-gradient scratch. Each facet's strain/stress is evaluated once for
    /// all columns. Material and geometric stiffness are both retained; negative
    /// tangent curvature is not clipped. The linear DKT remainder is counted once.
    /// No heap allocation, dense quartic tensor, altered modal basis or new law.
    /// This O(facets*n^2) kernel may round differently from repeated actions.
    /// Invalid dimensions/nonfinite input fills `out` with NaN; callers must
    /// reject nonfinite trial results before numerical or physical publication.
    pub fn hessian_matrix(&self, q: &[f64], out: &mut [f64], strain_gradients: &mut [[f64; 3]]) {
        let _ = self.hessian_matrix_with_cancel(q, out, strain_gradients, || false);
    }

    /// The same original Hessian arithmetic with cancellation at bounded rows
    /// and facet boundaries. A failure leaves caller-owned scratch partially
    /// written; it must not be published as a physical tangent.
    pub fn hessian_matrix_with_cancel(
        &self,
        q: &[f64],
        out: &mut [f64],
        strain_gradients: &mut [[f64; 3]],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), super::PlateError> {
        if cancelled() {
            return Err(super::PlateError::Cancelled);
        }
        let n = self.mode_count();
        if q.len() != n
            || n.checked_mul(n) != Some(out.len())
            || strain_gradients.len() != n
            || q.iter().any(|v| !v.is_finite())
        {
            out.fill(f64::NAN);
            return Err(super::bad(
                "shell Hessian dimensions or finite set are invalid",
            ));
        }
        out.copy_from_slice(&self.remainder);
        for f in &self.facets {
            if cancelled() {
                return Err(super::PlateError::Cancelled);
            }
            let (strain, dx, dy) = Self::strain(f, q);
            let stress = mul(&f.membrane, strain);
            for (b, s) in strain_gradients.iter_mut().zip(&f.strains) {
                *b = [
                    s.linear[0] + dot(dx, s.dx),
                    s.linear[1] + dot(dy, s.dy),
                    s.linear[2] + dot(dy, s.dx) + dot(dx, s.dy),
                ];
            }
            for i in 0..n {
                if i % 64 == 0 && cancelled() {
                    return Err(super::PlateError::Cancelled);
                }
                let stress_direction = mul(&f.membrane, strain_gradients[i]);
                let a = &f.strains[i];
                for j in 0..=i {
                    let b = &f.strains[j];
                    let second = [
                        dot(a.dx, b.dx),
                        dot(a.dy, b.dy),
                        dot(a.dx, b.dy) + dot(a.dy, b.dx),
                    ];
                    let tangent = dot(strain_gradients[j], stress_direction) + dot(second, stress);
                    out[i * n + j] += tangent;
                    if i != j {
                        out[j * n + i] += tangent;
                    }
                }
            }
        }
        if cancelled() {
            return Err(super::PlateError::Cancelled);
        }
        Ok(())
    }
}

// Material-memory consumers need the physical flexural part of this tangent.
#[path = "bending.rs"]
mod bending;

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::{Facet,ModeStrain};

    fn specimen(n:usize)->ShellReduction {
        let facets=(0..3).map(|f|Facet {
            membrane:[4.0,0.4,0.0, 0.4,3.0,0.0, 0.0,0.0,1.2], bending:[0.;9],
            strains:(0..n).map(|i| {
                let t=(1+i+5*f) as f64;
                let dx=[0.2*t.sin(),0.1*(0.8*t).cos(),0.3*(0.6*t).sin()];
                let dy=[0.1*(1.1*t).cos(),0.3*t.sin(),0.2*(0.4*t).cos()];
                ModeStrain {dx,dy,linear:[dx[0],dy[1],dx[1]+dy[0]]}
            }).collect(),
        }).collect();
        let mut remainder=vec![0.0;n*n];for i in 0..n {remainder[i*n+i]=1.0+i as f64;}
        ShellReduction {omegas:vec![1.0;n],remainder,facets,translations:vec![],rotations:vec![],
            reference_positions:vec![],section_thicknesses:vec![],nodes:0,triangles:vec![],
            area_normals:vec![],rigid_motion:None}
    }

    #[test]
    fn g0_batched_shell_tangent_matches_metric_energy_gradient_and_original_actions() {
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
    fn g0_shell_geometric_stiffness_is_not_replaced_with_a_positive_surrogate() {
        let mut law=specimen(1);law.remainder[0]=0.0;
        law.facets=vec![Facet {membrane:[1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0], bending:[0.;9],
            strains:vec![ModeStrain {dx:[1.0,0.0,0.0],dy:[0.0;3],linear:[1.0,0.0,0.0]}]}];
        let mut matrix=[0.0];let mut scratch=[[0.0;3]];
        law.hessian_matrix(&[-0.5],&mut matrix,&mut scratch);
        assert_eq!(matrix[0],-0.125); // (1+q)^2 + (q + q^2/2)
        law.hessian_matrix(&[0.0],&mut matrix,&mut scratch);assert_eq!(matrix[0],1.0);
        law.hessian_matrix(&[f64::NAN],&mut matrix,&mut scratch);assert!(matrix[0].is_nan());
        law.hessian_matrix(&[0.0],&mut matrix,&mut []);assert!(matrix[0].is_nan());
        law.hessian_matrix(&[],&mut matrix,&mut scratch);assert!(matrix[0].is_nan());
        law.hessian_matrix(&[0.0],&mut [],&mut scratch);
    }

    #[test]
    fn cancellable_hessian_reuses_original_arithmetic_and_refuses_partial_work() {
        let law = specimen(8);
        let q = vec![0.03; 8];
        let mut expected = vec![0.0; 64];
        let mut scratch = vec![[0.0; 3]; 8];
        law.hessian_matrix(&q, &mut expected, &mut scratch);
        let mut actual = vec![0.0; 64];
        law.hessian_matrix_with_cancel(&q, &mut actual, &mut scratch, || false)
            .unwrap();
        assert_eq!(actual, expected);
        let mut calls = 0;
        assert!(matches!(
            law.hessian_matrix_with_cancel(&q, &mut actual, &mut scratch, || {
                calls += 1;
                calls > 2
            }),
            Err(super::super::PlateError::Cancelled)
        ));
        law.hessian_matrix_with_cancel(&q, &mut actual, &mut scratch, || false)
            .unwrap();
        assert_eq!(actual, expected);
    }
}
