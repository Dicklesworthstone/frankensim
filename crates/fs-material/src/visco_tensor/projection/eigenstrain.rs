//! Affine free strain in the original quadrature and factor coordinates.
use super::*;

impl ProjectedTensorProny {
    /// Bind an explicitly admitted total stress-free strain in this point's
    /// Mandel frame. Geometry, reference mass and every original strain column
    /// remain unchanged. The caller owns its source and temperature history.
    /// A changing Kelvin eigenstrain needs a resolved thermal strain rate and
    /// is deliberately outside the endpoint Maxwell continuation model.
    pub fn with_eigenstrain_mandel(&self, free: Voigt) -> Result<Self> {
        self.model
            .validate_mandel_strain(&free)
            .map_err(|_| invalid("free strain exceeds original small-strain domain"))?;
        if free.iter().any(|x| *x != 0.0) && self.model.kv_pa_s().is_some() {
            return Err(invalid(
                "Kelvin thermal expansion requires a resolved strain-rate coupling",
            ));
        }
        let mut next = self.clone();
        next.free_strain_mandel = free;
        Ok(next)
    }

    /// Total stress-free strain in the same Mandel frame as the material.
    pub const fn free_strain_mandel(&self) -> &Voigt {
        &self.free_strain_mandel
    }

    /// Equilibrium affine load and constant free energy: psi_eq =
    /// .5 u^T K u - u^T load + constant. These are internal material storage,
    /// not an external force port or a source of dissipative heat.
    pub fn equilibrium_eigenstrain_terms(&self) -> Result<(Vec<f64>, f64)> {
        let stress = apply(&self.model.c_eq, &self.free_strain_mandel);
        let load: Vec<_> = self
            .columns
            .iter()
            .map(|column| self.volume_m3 * dot(&mandel(column), &stress))
            .collect();
        let constant = 0.5 * self.volume_m3 * dot(&self.free_strain_mandel, &stress);
        if load.iter().any(|v| !v.is_finite()) || !constant.is_finite() || constant < 0.0 {
            return Err(invalid("equilibrium eigenstrain storage overflow"));
        }
        Ok((load, constant))
    }

    /// Original volume-scaled branch factors applied to the total free strain.
    /// Subtract these from A*u only when explicitly initializing instantaneous
    /// memory. Existing histories must use transfer_factor_memory instead.
    pub fn eigenstrain_factor_coordinates(&self) -> Result<Vec<f64>> {
        let root = det::sqrt(self.volume_m3);
        let mut values = Vec::with_capacity(6 * self.model.branches.len());
        for (branch, factor) in self.model.branches.iter().zip(&self.model.factors) {
            values.extend(
                apply(&branch.resolved_gram_rows(factor), &self.free_strain_mandel)
                    .map(|v| root * v),
            );
        }
        if values.iter().any(|v| !v.is_finite()) {
            return Err(invalid("eigenstrain factor coordinates overflow"));
        }
        Ok(values)
    }

    /// Positive magnitude of -dU/dT / T at fixed total and viscous strains,
    /// for CONSTANT stored tensors and constant free-strain derivative alpha.
    /// The complete equilibrium and Maxwell tensors contribute, including
    /// constrained directions absent from a reduced displacement basis.
    pub fn expansion_capacity_coefficient_j_k2(&self, alpha_mandel: Voigt) -> Result<f64> {
        if !finite(&alpha_mandel) {
            return Err(invalid("nonfinite expansion derivative"));
        }
        let mut value = dot(&alpha_mandel, &apply(&self.model.c_eq, &alpha_mandel));
        for branch in &self.model.branches {
            value += dot(&alpha_mandel, &apply(&branch.c_pa, &alpha_mandel));
        }
        value *= self.volume_m3;
        if !value.is_finite() || value < 0.0 {
            return Err(invalid("thermal capacity contribution overflow"));
        }
        Ok(value)
    }
}
