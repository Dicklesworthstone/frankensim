//! Original Maxwell branch work, without a global correction of kinetic effort.
//!
//! For e=s(x)-z, H_arm=e^2/2 and globally quadratic s, the symmetric rule is
//! DG_x H=mean(e_a,e_b)*mean(grad s_a,grad s_b), DG_z H=-mean(e_a,e_b).
//! Its contraction is the original branch energy change. Differentiating this
//! rule retains prestress curvature and mechanical/memory cross terms. Neither
//! the strain law nor its energy is fitted, linearized or replaced.

use super::{RelaxationCoupling, RelaxationStrainMap};
use crate::interconnection::accumulation::Sum;
use crate::{ConstitutiveDiscreteGradient, RelaxationStorage};

/// A constant analytic derivative image of an explicitly declared polynomial.
/// The original callback still evaluates strain, including its source domain.
pub(crate) struct QuadraticStrainImage {
    linear: Box<[f64]>,
    hessian: Box<[f64]>,
}

impl QuadraticStrainImage {
    pub(crate) fn prepare(strain: &dyn RelaxationStrainMap, scale: f64) -> Option<Self> {
        let n = strain.state_dimension();
        if !strain.is_quadratic_strain()
            || n > crate::linear_exact::MAX_LINEAR_EXPONENTIAL_BLOCK_DIM
        {
            return None;
        }
        // A global quadratic declaration, not sampled agreement, establishes
        // the constant Hessian. Basis-vector analytic actions only materialize
        // its original coefficients. At most 256^2+256 doubles per mapped arm.
        let origin = vec![0.0; n];
        let mut linear = vec![f64::NAN; n];
        let mut hessian = vec![f64::NAN; n * n];
        let mut direction = vec![0.0; n];
        let mut action = vec![f64::NAN; n];
        strain.gradient(&origin, &mut linear);
        for value in &mut linear {
            let scaled = *value * scale;
            if scaled == 0.0 && *value != 0.0 {
                return None;
            }
            *value = scaled;
        }
        for j in 0..n {
            direction[j] = 1.0;
            action.fill(f64::NAN);
            strain.hessian_vector(&origin, &direction, &mut action);
            direction[j] = 0.0;
            for i in 0..n {
                let scaled = scale * action[i];
                if scaled == 0.0 && action[i] != 0.0 {
                    return None;
                }
                hessian[i * n + j] = scaled;
            }
        }
        if linear.iter().chain(&hessian).any(|v| !v.is_finite()) {
            // A polynomial can have a physical source domain excluding the
            // coordinate origin. A missing image only declines this cache;
            // it must not reject a previously valid original material law.
            return None;
        }
        for i in 0..n {
            for j in 0..i {
                if hessian[i * n + j] != hessian[j * n + i] {
                    return None;
                }
            }
        }
        Some(Self {
            linear: linear.into_boxed_slice(),
            hessian: hessian.into_boxed_slice(),
        })
    }

    pub(crate) fn gradient(&self, x: &[f64], i: usize) -> f64 {
        let mut sum = Sum::default();
        let result = (|| {
            sum.product(1.0, self.linear[i])?;
            for j in 0..x.len() {
                sum.product(self.hessian[i * x.len() + j], x[j])?;
            }
            sum.value()
        })();
        result.unwrap_or(f64::NAN)
    }

    // Average endpoint gradients, retaining represented-midpoint rounding.
    fn mean_gradient(&self, a: &[f64], b: &[f64], i: usize) -> f64 {
        let mut sum = Sum::default();
        let result = (|| {
            sum.product(1.0, self.linear[i])?;
            for j in 0..a.len() {
                let mid = f64::midpoint(a[j], b[j]);
                let h = self.hessian[i * a.len() + j];
                sum.product(h, mid)?;
                sum.product(0.5 * h, (a[j] - mid) + (b[j] - mid))?;
            }
            sum.value()
        })();
        result.unwrap_or(f64::NAN)
    }

    pub(crate) fn curvature(&self, direction: &[f64], i: usize) -> f64 {
        let n = direction.len();
        Sum::dot(&self.hessian[i * n..(i + 1) * n], direction)
            .and_then(Sum::value)
            .unwrap_or(f64::NAN)
    }

    /// Same independent polynomial remainder as the original mapped owner,
    /// contracted from coordinate increments without allocating derivative rows.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn remainder(
        &self,
        a: &[f64],
        b: &[f64],
        mid: &[f64],
        za: f64,
        zb: f64,
        zm: f64,
        strain_mid: f64,
    ) -> f64 {
        let mut change = Sum::default();
        let mut offset = Sum::default();
        let mut cross = Sum::default();
        let result = (|| {
            change.product(-1.0, zb - za)?;
            offset.product(-1.0, (za - zm) + (zb - zm))?;
            for i in 0..a.len() {
                let (mut hd, mut ha, mut hb) = (Sum::default(), Sum::default(), Sum::default());
                for j in 0..a.len() {
                    let h = self.hessian[i * a.len() + j];
                    hd.product(h, b[j] - a[j])?;
                    ha.product(h, a[j] - mid[j])?;
                    hb.product(h, b[j] - mid[j])?;
                }
                let round = (a[i] - mid[i]) + (b[i] - mid[i]);
                let gradient = self.gradient(mid, i);
                cross.add_scaled(hd, 0.5 * round)?;
                change.product(gradient, b[i] - a[i])?;
                offset.product(gradient, round)?;
                offset.add_scaled(ha, 0.5 * (a[i] - mid[i]))?;
                offset.add_scaled(hb, 0.5 * (b[i] - mid[i]))?;
            }
            change.add_scaled(cross, 1.0)?;
            let mut result = Sum::default();
            let (leading, tail) = change.components()?;
            result.add_half_scaled(offset, leading)?;
            if tail != 0.0 {
                result.add_half_scaled(offset, tail)?;
            }
            result.add_scaled(cross, strain_mid - zm)?;
            result.value()
        })();
        result.unwrap_or(f64::NAN)
    }
}

impl RelaxationStorage {
    pub(crate) fn has_constitutive_rule(&self) -> bool {
        // Preserve existing unsupported/global-Gonzalez behavior when an arm
        // lacks a declared image. Nested native constitutive laws remain native.
        (self.base.discrete_gradient_rule().is_some()
            || self.base.exact_quadratic_matrix().is_some()
            || self.polynomial.iter().any(Option::is_some))
            && self
                .couplings
                .iter()
                .zip(&self.polynomial)
                .all(|(c, p)| matches!(c, RelaxationCoupling::Linear(_)) || p.is_some())
    }

    fn mean_strain_gradient(&self, arm: usize, a: &[f64], b: &[f64], i: usize) -> f64 {
        match &self.couplings[arm] {
            RelaxationCoupling::Linear(row) => row[i],
            RelaxationCoupling::Mapped { .. } => self.polynomial[arm]
                .as_ref()
                .map_or(f64::NAN, |image| image.mean_gradient(a, b, i)),
        }
    }

    fn strain_gradient(&self, arm: usize, x: &[f64], i: usize) -> f64 {
        match &self.couplings[arm] {
            RelaxationCoupling::Linear(row) => row[i],
            RelaxationCoupling::Mapped { .. } => self.polynomial[arm]
                .as_ref()
                .map_or(f64::NAN, |image| image.gradient(x, i)),
        }
    }
}

impl ConstitutiveDiscreteGradient for RelaxationStorage {
    fn evaluate(&self, a: &[f64], b: &[f64], mid: &[f64], out: &mut [f64]) {
        let n = self.base_dim + self.couplings.len();
        out.fill(f64::NAN);
        if a.len() != n
            || b.len() != n
            || mid.len() != n
            || out.len() != n
            || a.iter().chain(b).chain(mid).any(|v| !v.is_finite())
        {
            return;
        }
        let d = self.base_dim;
        crate::discrete_gradient_from_midpoint_unchecked(
            self.base.as_ref(),
            &a[..d],
            &b[..d],
            &mid[..d],
            &mut out[..d],
            d,
        );
        for (arm, coupling) in self.couplings.iter().enumerate() {
            let z = d + arm;
            let effort = f64::midpoint(
                coupling.strain(&a[..d]) - a[z],
                coupling.strain(&b[..d]) - b[z],
            );
            for (i, value) in out[..d].iter_mut().enumerate() {
                *value += effort * self.mean_strain_gradient(arm, &a[..d], &b[..d], i);
            }
            out[z] = -effort;
        }
    }

    fn end_jacobian_action(
        &self,
        a: &[f64],
        b: &[f64],
        direction: &[f64],
        out: &mut [f64],
    ) -> bool {
        let n = self.base_dim + self.couplings.len();
        if a.len() != n
            || b.len() != n
            || direction.len() != n
            || out.len() != n
            || a.iter().chain(b).chain(direction).any(|v| !v.is_finite())
        {
            return false;
        }
        let d = self.base_dim;
        if !crate::component_end_jacobian_action(
            self.base.as_ref(),
            &a[..d],
            &b[..d],
            &direction[..d],
            &mut out[..d],
        ) {
            return false;
        }
        for (arm, coupling) in self.couplings.iter().enumerate() {
            let z = d + arm;
            let effort = f64::midpoint(
                coupling.strain(&a[..d]) - a[z],
                coupling.strain(&b[..d]) - b[z],
            );
            let mut change = Sum::default();
            if change.product(-0.5, direction[z]).is_none() {
                return false;
            }
            for (i, &v) in direction[..d].iter().enumerate() {
                if change
                    .product(0.5 * self.strain_gradient(arm, &b[..d], i), v)
                    .is_none()
                {
                    return false;
                }
            }
            let Some(change) = change.value() else {
                return false;
            };
            for (i, value) in out[..d].iter_mut().enumerate() {
                let curvature = self.polynomial[arm]
                    .as_ref()
                    .map_or(0.0, |image| image.curvature(&direction[..d], i));
                *value += change * self.mean_strain_gradient(arm, &a[..d], &b[..d], i)
                    + (0.5 * effort) * curvature;
            }
            out[z] = -change;
        }
        out.iter().all(|v| v.is_finite())
    }
}
