//! Physical chain rule on the original descriptor model, in a fixed state chart.
//! With S=J-R and e=W*v: dA*v=(dJ-dR)*e+S*(dQ*v)-dT*v. T is the optional
//! total nonlinear resisting-flow derivative. No storage inverse is required.
use super::*;

/// Derivatives per unit of ONE caller-declared real parameter, evaluated at the
/// prepared model. None means exactly zero. Matrices are row-major in original
/// state/port order. The differential/algebraic split and state chart stay fixed.
/// Storage and dissipation tangents are symmetric, interconnection is skew.
/// Signed tangents need NOT be positive semidefinite: decreasing mass, stiffness
/// or loss is a valid local direction. No finite parameter interval is admitted.
#[derive(Clone, Copy, Debug, Default)]
pub struct HarmonicParameterDirection<'a> {
    pub interconnection: Option<&'a [f64]>,
    pub dissipation: Option<&'a [f64]>,
    /// Full n*n dQ with exactly zero algebraic rows/columns.
    pub storage: Option<&'a [f64]>,
    /// n*ports dG, used in BOTH source and observation semantics. Free-mode
    /// eigenvalues do not depend on G, but their observed output does.
    pub ports: Option<&'a [f64]>,
    /// Full n*n TOTAL derivative of T=D_x+D_e W, not an extra R matrix.
    /// No symmetry/PSD assumption. Include every constitutive, Hessian and
    /// operating-point chain-rule term; None holds T fixed independently of Q.
    /// A changing bias is not solved or differentiated automatically here.
    pub resisting_flow: Option<&'a [f64]>,
}

/// Derivatives per unit of the declared parameter, in v[pivot]=1 gauge.
/// Energy/loss derivatives are gauge-dependent observations, not derivatives
/// of a forced response or a positive modal-energy normalization.
#[derive(Clone, Debug)]
pub struct HarmonicPoleDerivative {
    pub laplace_pole_per_s: C64,
    pub angular_frequency_rad_s: C64,
    pub state: Vec<C64>,
    pub effort: Vec<C64>,
    pub output: Vec<C64>,
    pub stored_energy_j: f64,
    /// Includes the derivative of signed incremental work for nonlinear images.
    pub dissipated_power_w: f64,
    /// d(2 Re(s) H + D), independently formed, not used to invent loss.
    pub decay_power_residual_w: f64,
    pub relative_residual: f64,
    pub linear_residual: f64,
}
fn at(matrix: Option<&[f64]>, index: usize) -> f64 { matrix.map_or(0.0, |a| a[index]) }

impl HarmonicParameterDirection<'_> {
    fn validate(&self, model: &HarmonicModel, cancelled: &mut impl FnMut() -> bool)
        -> Result<(), HarmonicPoleResponseError>
    {
        let n = model.n;
        for matrix in [self.interconnection, self.dissipation, self.storage, self.resisting_flow].into_iter().flatten() {
            check_cancel(cancelled)?;
            if matrix.len() != n*n || matrix.iter().any(|v| !v.is_finite()) {
                return Err(invalid("complete finite n*n parameter matrices required"));
            }
        }
        if let Some(g) = self.ports {
            if g.len() != n*model.m || g.iter().any(|v| !v.is_finite()) {
                return Err(invalid("complete finite n*ports port derivative required"));
            }
        }
        // Reuse the source owner's symmetry admission, never repair a tangent.
        if let Some(j) = self.interconnection {
            crate::require_skew(j,n,"interconnection derivative")
                .map_err(|_| invalid("interconnection derivative must be skew"))?;
        }
        for matrix in [self.dissipation,self.storage].into_iter().flatten() {
            crate::require_symmetric(matrix,n,"material derivative")
                .map_err(|_| invalid("storage and dissipation derivatives must be symmetric"))?;
        }
        if let Some(q) = self.storage {
            for i in 0..n {
                check_cancel(cancelled)?;
                for j in 0..n {
                    if (i >= model.n_diff || j >= model.n_diff) && q[i*n+j] != 0.0 {
                        return Err(invalid("algebraic efforts cannot acquire storage through a derivative"));
                    }
                }
            }
        }
        check_cancel(cancelled)
    }
}

impl PreparedHarmonicPole<'_> {
    /// Differentiate free poles, full eigenstates, efforts and port outputs.
    /// Reuses this factor, with O(n^2+n*ports) work per direction. No perturbed
    /// model, finite-difference step, mode rematching or storage inverse.
    /// The caller supplies the complete constitutive/circuit chain rule,
    /// including dT on nonlinear local models. No change of coordinates,
    /// topology or differential rank is inferred. Batches are transactional.
    pub fn parameter_derivatives(
        &self, directions: &[HarmonicParameterDirection<'_>], mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<HarmonicPoleDerivative>, HarmonicPoleResponseError> {
        check_cancel(&mut cancelled)?;
        if directions.is_empty() || directions.len() > self.options.numerical.max_vectors {
            return Err(invalid("complete nonempty parameter batch within its ceiling required"));
        }
        let model = self.model;
        let (n,m) = (model.n,model.m);
        let v = self.right_state();
        let s = self.laplace_pole_per_s();
        let h = 0.5*(0..model.n_diff).map(|i|(v[i].conj()*self.effort[i]).re).sum::<f64>();
        let mut results = Vec::with_capacity(directions.len());
        for direction in directions {
            direction.validate(model,&mut cancelled)?;
            let mut dqv = vec![C64::ZERO;n];
            let mut dq_envelope = vec![0.0;n];
            for i in 0..model.n_diff {
                check_cancel(&mut cancelled)?;
                let mut correction = C64::ZERO;
                for j in 0..model.n_diff {
                    let q = at(direction.storage,i*n+j);
                    accumulate(&mut dqv[i],&mut correction,v[j].scale(q));
                    dq_envelope[i] += q.abs()*v[j].abs();
                }
            }
            let mut action = vec![C64::ZERO;n];
            let mut envelope = vec![0.0;n];
            for i in 0..n {
                check_cancel(&mut cancelled)?;
                let mut correction = C64::ZERO;
                for j in 0..n {
                    let dj = at(direction.interconnection,i*n+j);
                    let dr = at(direction.dissipation,i*n+j);
                    let a = self.effort[j].scale(dj-dr)+dqv[j].scale(model.structure[i*n+j]);
                    accumulate(&mut action[i],&mut correction,-a);
                    envelope[i] += (dj.abs()+dr.abs())*self.effort[j].abs()
                        +model.structure[i*n+j].abs()*dq_envelope[j];
                    if let Some(dt)=direction.resisting_flow {
                        accumulate(&mut action[i],&mut correction,v[j].scale(dt[i*n+j]));
                        envelope[i]+=dt[i*n+j].abs()*v[j].abs();
                    }
                }
            }
            if action.iter().chain(&dqv).any(|&z| !finite(z))
                || envelope.iter().chain(&dq_envelope).any(|x| !x.is_finite()) {
                return Err(invalid("physical parameter action is not representable"));
            }
            let scaled: Vec<_> = action.iter().zip(&self.row_scales).map(|(&z,&r)|divide(z,r)).collect();
            let mut solved = self.numerical.parameter_derivatives(&[&scaled],&mut cancelled)?;
            let tangent = solved.pop().expect("one admitted parameter");
            let ds = tangent.value;
            let dv = tangent.vector;
            let mut residual = 0.0_f64;
            for i in 0..n {
                check_cancel(&mut cancelled)?;
                let scale = self.row_scales[i];
                let mut sum = divide(action[i],scale);
                let mut correction = C64::ZERO;
                let mut denominator = envelope[i]/scale;
                if i < model.n_diff {
                    let a = divide(s,scale)*dv[i];
                    let b = divide(ds,scale)*v[i];
                    accumulate(&mut sum,&mut correction,a);
                    accumulate(&mut sum,&mut correction,b);
                    denominator += a.abs()+b.abs();
                }
                for j in 0..n {
                    let a = dv[j].scale(model.generator[i*n+j]/scale);
                    accumulate(&mut sum,&mut correction,-a);
                    denominator += a.abs();
                }
                if !finite(sum) || !denominator.is_finite() {
                    return Err(invalid("original tangent residual overflow"));
                }
                let row = if denominator == 0.0 { if sum.abs()==0.0 {0.0} else {f64::INFINITY} }
                    else {sum.abs()/denominator};
                residual = residual.max(row);
            }
            if !residual.is_finite() || residual > model.budget.relative_residual_tolerance {
                return Err(HarmonicPoleResponseError::Physical(HarmonicPoleError::Residual {
                    relative:residual,tolerance:model.budget.relative_residual_tolerance,
                }));
            }
            let mut effort = dv.clone();
            let mut dh = 0.0;
            for i in 0..model.n_diff {
                check_cancel(&mut cancelled)?;
                effort[i]=dqv[i];
                let mut correction=C64::ZERO;
                for j in 0..model.n_diff {
                    accumulate(&mut effort[i],&mut correction,dv[j].scale(model.q[i*n+j]));
                }
                dh += (dv[i].conj()*self.effort[i]).re+0.5*(v[i].conj()*dqv[i]).re;
            }
            let mut output=vec![C64::ZERO;m];
            for p in 0..m {
                check_cancel(&mut cancelled)?;
                let mut correction=C64::ZERO;
                for i in 0..n {
                    accumulate(&mut output[p],&mut correction,effort[i].scale(model.g[i*m+p])
                        +self.effort[i].scale(at(direction.ports,i*m+p)));
                }
            }
            let mut dd=0.0;
            for i in 0..n {
                check_cancel(&mut cancelled)?;
                for j in 0..n {
                    dd += 2.0*(effort[i].conj()*self.effort[j]).re*model.r[i*n+j]
                        +(self.effort[i].conj()*self.effort[j]).re*at(direction.dissipation,i*n+j);
                    if let Some(t)=&model.resisting_tangent {
                        dd += (effort[i].conj()*v[j]+self.effort[i].conj()*dv[j]).re*t[i*n+j];
                    }
                    if let Some(dt)=direction.resisting_flow {
                        dd += (self.effort[i].conj()*v[j]).re*dt[i*n+j];
                    }
                }
            }
            let balance=2.0*ds.re*h+2.0*s.re*dh+dd;
            if effort.iter().chain(&output).any(|&z| !finite(z))
                || [h,dh,dd,balance].iter().any(|x| !x.is_finite()) {
                return Err(invalid("physical tangent observation overflow"));
            }
            results.push(HarmonicPoleDerivative {
                laplace_pole_per_s:ds,angular_frequency_rad_s:C64::new(-ds.im,ds.re),
                state:dv,effort,output,stored_energy_j:dh,dissipated_power_w:dd,
                decay_power_residual_w:balance,relative_residual:residual,
                linear_residual:tangent.linear_residual,
            });
        }
        check_cancel(&mut cancelled)?;
        Ok(results)
    }
}
