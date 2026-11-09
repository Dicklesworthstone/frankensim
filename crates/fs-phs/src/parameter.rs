//! Original physical coefficient directions without a linearized storage model.

use crate::{PhsError, PortHamiltonian, interconnection::accumulation::Sum};

/// One real parameter's complete original J/R/G derivatives, row-major.
/// None holds that matrix fixed. J's direction is skew and R's is symmetric;
/// signed directions need NOT be positive semidefinite. The physical state and
/// port charts stay fixed. This admits a local derivative, not a finite update.
#[derive(Clone, Copy, Debug, Default)]
pub struct PhsStructureParameterDirection<'a> {
    pub interconnection: Option<&'a [f64]>,
    pub dissipation: Option<&'a [f64]>,
    /// The SAME n*ports derivative changes both excitation and observation.
    pub ports: Option<&'a [f64]>,
}

/// Validated coefficient direction borrowed together with its original source.
/// No storage, physical matrix or perturbed model is copied or reconstructed.
pub struct PreparedPhsParameterDirection<'system, 'direction> {
    system: &'system PortHamiltonian,
    direction: PhsStructureParameterDirection<'direction>,
}

fn invalid(what: &'static str) -> PhsError {
    PhsError::Unadmitted { what }
}

impl PortHamiltonian {
    /// Gradient of the original resistive power `e^T R e` with respect to
    /// the complete original storage effort. Computes `(R + R^T) e` with
    /// one compensated sum per output, retaining the source's actual R.
    ///
    /// A trusted raw-parts source may be nonsymmetric or active; neither a
    /// symmetry shortcut nor positivity repair is applied. Composing this
    /// covector with the original storage Hessian gives its state gradient.
    /// Work is O(n²), with no allocation and constant scratch.
    ///
    /// # Errors
    /// Incomplete/nonfinite source or efforts, or unrepresentable contraction.
    /// Extent/input failures precede writes. Arithmetic refusal leaves scratch.
    pub fn resistance_power_effort_gradient_into(
        &self,
        effort: &[f64],
        out: &mut [f64],
    ) -> Result<(), PhsError> {
        if effort.len() != self.n
            || out.len() != self.n
            || self.n.checked_mul(self.n) != Some(self.r.len())
            || effort.iter().chain(&self.r).any(|v| !v.is_finite())
        {
            return Err(invalid(
                "complete finite resistance power gradient inputs required",
            ));
        }
        for (i, value) in out.iter_mut().enumerate() {
            let mut sum = Sum::default();
            for (j, &e) in effort.iter().enumerate() {
                for r in [self.r[i * self.n + j], self.r[j * self.n + i]] {
                    if r != 0.0 && e != 0.0 {
                        sum.product(r, e).ok_or_else(|| {
                            invalid("resistance power gradient is not representable")
                        })?;
                    }
                }
            }
            *value = sum
                .value()
                .ok_or_else(|| invalid("resistance power gradient is not representable"))?;
        }
        Ok(())
    }

    /// Bind complete finite coefficient derivatives to this original system.
    /// Reuses the source's skew/symmetry admission without symmetrizing or
    /// imposing positivity on a signed derivative. No storage inverse,
    /// eigenproblem or source reassembly is needed; work is O(n²+n*ports).
    ///
    /// # Errors
    /// Incomplete/nonfinite source or direction matrices, or wrong derivative
    /// symmetry. Trusted raw sources retain their original admission contract.
    pub fn prepare_parameter_direction<'s, 'd>(
        &'s self,
        direction: PhsStructureParameterDirection<'d>,
    ) -> Result<PreparedPhsParameterDirection<'s, 'd>, PhsError> {
        let square = self.n.checked_mul(self.n);
        let ports = self.n.checked_mul(self.m);
        if square != Some(self.j.len())
            || square != Some(self.r.len())
            || ports != Some(self.g.len())
            || self
                .j
                .iter()
                .chain(&self.r)
                .chain(&self.g)
                .any(|v| !v.is_finite())
        {
            return Err(invalid(
                "complete finite original parameter-source matrices required",
            ));
        }
        for matrix in [direction.interconnection, direction.dissipation]
            .into_iter()
            .flatten()
        {
            if square != Some(matrix.len()) || matrix.iter().any(|v| !v.is_finite()) {
                return Err(invalid(
                    "complete finite n*n coefficient derivative required",
                ));
            }
        }
        if let Some(g) = direction.ports {
            if ports != Some(g.len()) || g.iter().any(|v| !v.is_finite()) {
                return Err(invalid(
                    "complete finite n*ports coefficient derivative required",
                ));
            }
        }
        if let Some(j) = direction.interconnection {
            crate::require_skew(j, self.n, "interconnection derivative")?;
        }
        if let Some(r) = direction.dissipation {
            crate::require_symmetric(r, self.n, "dissipation derivative")?;
        }
        Ok(PreparedPhsParameterDirection {
            system: self,
            direction,
        })
    }
}

impl PreparedPhsParameterDirection<'_, '_> {
    /// Apply `(J-R)e_p + (J_p-R_p)e + G_p u + G u_p` in original units.
    /// The supplied effort derivative can be a direct constitutive derivative
    /// at fixed state, or a total derivative including Hxx*x_p. The caller owns
    /// that distinction and the complete original storage/drive chain rule.
    /// Every original product shares one compensated sum per equation, so a
    /// weak source is not lost by rounding opposing contributions separately.
    /// No allocation; O(n²+n*ports) work and constant scratch.
    ///
    /// # Errors
    /// Incomplete/nonfinite inputs or unrepresentable contraction. Shape/input
    /// failures precede output writes; arithmetic refusal leaves scratch output.
    pub fn vector_field_from_effort_into(
        &self,
        effort: &[f64],
        effort_derivative: &[f64],
        input: &[f64],
        input_derivative: &[f64],
        out: &mut [f64],
    ) -> Result<(), PhsError> {
        let (n, m) = (self.system.n, self.system.m);
        if effort.len() != n
            || effort_derivative.len() != n
            || out.len() != n
            || input.len() != m
            || input_derivative.len() != m
            || effort
                .iter()
                .chain(effort_derivative)
                .chain(input)
                .chain(input_derivative)
                .any(|v| !v.is_finite())
        {
            return Err(invalid(
                "complete finite physical parameter-flow inputs required",
            ));
        }
        for (row, value) in out.iter_mut().enumerate() {
            let mut sum = Sum::default();
            let mut product = |a: f64, b: f64| {
                if a != 0.0 {
                    sum.product(a, b)
                        .ok_or_else(|| invalid("parameter flow is not representable"))?;
                }
                Ok::<(), PhsError>(())
            };
            for col in 0..n {
                let index = row * n + col;
                product(self.system.j[index], effort_derivative[col])?;
                product(-self.system.r[index], effort_derivative[col])?;
                if let Some(j) = self.direction.interconnection {
                    product(j[index], effort[col])?;
                }
                if let Some(r) = self.direction.dissipation {
                    product(-r[index], effort[col])?;
                }
            }
            for port in 0..m {
                let index = row * m + port;
                if let Some(g) = self.direction.ports {
                    product(g[index], input[port])?;
                }
                product(self.system.g[index], input_derivative[port])?;
            }
            *value = sum
                .value()
                .ok_or_else(|| invalid("parameter flow is not representable"))?;
        }
        Ok(())
    }

    /// Observe `G^T e_p + G_p^T e` using the same reciprocal port direction.
    /// Both terms share the original compensated contraction before rounding.
    /// No allocation; O(n*ports) work. Effort derivatives are not state slopes.
    /// # Errors
    /// Incomplete/nonfinite efforts or unrepresentable physical output.
    pub fn output_from_effort_into(
        &self,
        effort: &[f64],
        effort_derivative: &[f64],
        out: &mut [f64],
    ) -> Result<(), PhsError> {
        let (n, m) = (self.system.n, self.system.m);
        if effort.len() != n
            || effort_derivative.len() != n
            || out.len() != m
            || effort
                .iter()
                .chain(effort_derivative)
                .any(|v| !v.is_finite())
        {
            return Err(invalid(
                "complete finite physical parameter-output efforts required",
            ));
        }
        for (port, value) in out.iter_mut().enumerate() {
            let mut sum = Sum::default();
            for row in 0..n {
                let index = row * m + port;
                if self.system.g[index] != 0.0 {
                    sum.product(self.system.g[index], effort_derivative[row])
                        .ok_or_else(|| invalid("parameter output is not representable"))?;
                }
                if let Some(g) = self.direction.ports {
                    if g[index] != 0.0 {
                        sum.product(g[index], effort[row])
                            .ok_or_else(|| invalid("parameter output is not representable"))?;
                    }
                }
            }
            *value = sum
                .value()
                .ok_or_else(|| invalid("parameter output is not representable"))?;
        }
        Ok(())
    }
}
