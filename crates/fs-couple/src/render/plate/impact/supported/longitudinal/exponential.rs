//! Numerical split of the original moving-support strain storage.
//!
//! All physical coefficients, strain maps, signed endpoint reactions and
//! mass transformations remain in SupportedStrings. Only its declared
//! quadratic polynomial and full source tangent envelope are prepared here.
use super::{relaxation::SupportedStrainMap, *};
use fs_math::dd::Dd;
use fs_phs::{
    QuadraticStorage, RelaxationStrainMap, Storage, quadratic_strain_midpoint_energy_remainder,
};

pub(in crate::render::plate::impact) struct SupportedLinearSplit {
    stiffness: Vec<f64>,
    strain_quadratic: QuadraticStorage,
    strains: Vec<(SupportedStrainMap, f64)>,
    modes: usize,
}

impl SupportedStrings {
    /// Numerical metadata only. No state, physical rigidity, original band
    /// or source clock is added/replaced. Repreparation is idempotent.
    pub(in crate::render::plate::impact) fn with_exponential_envelope(
        mut self,
    ) -> Result<Self, ImpactError> {
        self.require_exponential_strain()?;
        if self.0.exponential_equilibrium_bound.is_some() {
            return Ok(self);
        }
        let mut bound = 0.0;
        for (map, stiffness) in self.equilibrium_strains()? {
            // SAME canonical map and declared slope/strain domain as storage.
            // The old small-signal bound and every original Maxwell stress
            // envelope remain. Adding the full equilibrium envelope is
            // conservative; no small-signal norm is subtracted from a bound.
            bound += stiffness * map.stiffness_bound(map.potential().maximum_strain)?;
        }
        if !bound.is_finite() || bound < 0.0 {
            return Err(invalid("supported equilibrium tangent envelope overflows"));
        }
        let data = Rc::make_mut(&mut self.0);
        for omega in &mut data.clock_bounds {
            *omega = det::sqrt(*omega * *omega + bound);
            if !omega.is_finite() {
                return Err(invalid("complete supported exponential envelope overflows"));
            }
        }
        data.exponential_equilibrium_bound = Some(bound);
        Ok(self)
    }

    pub(in crate::render::plate::impact) fn has_exponential_envelope(&self) -> bool {
        self.0.exponential_equilibrium_bound.is_some()
    }

    fn require_exponential_strain(&self) -> Result<(), ImpactError> {
        if self
            .0
            .potentials
            .iter()
            .any(|p| !matches!(p, SupportedPotential::Longitudinal(_)))
        {
            return Err(invalid(
                "exponential supports require every original longitudinal strain and explicit geometry",
            ));
        }
        if self.0.finite_strain_axial_viscosity {
            return Err(invalid(
                "exponential supports cannot replace configuration-dependent axial viscosity by constant resistance",
            ));
        }
        Ok(())
    }

    fn equilibrium_strains(&self) -> Result<Vec<(SupportedStrainMap, f64)>, ImpactError> {
        self.require_exponential_strain()?;
        let mut strains = Vec::new();
        for (string, potential) in self.0.potentials.iter().enumerate() {
            let SupportedPotential::Longitudinal(p) = potential else {
                unreachable!()
            };
            for channel in 0..=p.model.axial_modes().len() {
                // Original fs-nlmodal coefficient, not a recomputed material EA
                // or a rounded square-root energy normalization.
                let stiffness =
                    p.model.strain_energy_coefficient() * if channel == 0 { 2.0 } else { 1.0 };
                if !stiffness.is_finite() || stiffness <= 0.0 {
                    return Err(invalid(
                        "supported equilibrium strain normalization overflows",
                    ));
                }
                strains.push((
                    SupportedStrainMap {
                        support: self.clone(),
                        string,
                        channel,
                        offset: 0,
                        base_dim: 2 * self.mode_count(),
                    },
                    stiffness,
                ));
            }
        }
        Ok(strains)
    }

    pub(in crate::render::plate::impact) fn exponential_split(
        &self,
    ) -> Result<SupportedLinearSplit, ImpactError> {
        let strains = self.equilibrium_strains()?;
        let n = self.mode_count();
        let support = self.support_coordinate();
        let mut stiffness = vec![0.0; n * n];
        let mut strain_matrix = vec![0.0; n * n];
        let weight = self.0.weight;
        stiffness[support * n + support] = self.0.support.stiffness_n_m * weight * weight;
        for (wire, potential) in self.0.potentials.iter().enumerate() {
            let SupportedPotential::Longitudinal(p) = potential else {
                unreachable!()
            };
            stiffness[support * n + support] += p.lift_stiffness * weight * weight;
            for (mode, &omega) in self.0.ranges[wire].clone().zip(&p.omegas) {
                // q_mode = z_mode - translation_mode * weight * z_support.
                // Retain the full reciprocal dyad, including support inertia.
                let k = omega * omega;
                let a = self.0.translation[mode] * weight;
                stiffness[mode * n + mode] += k;
                stiffness[mode * n + support] -= k * a;
                stiffness[support * n + mode] = stiffness[mode * n + support];
                stiffness[support * n + support] += k * a * a;
            }
        }
        for (map, k) in &strains {
            // The balanced source has epsilon(0)=0. Its declared quadratic
            // strain gives an exact polynomial quadratic term k*B0^T*B0;
            // no finite-state tangent replaces the original nonlinear energy.
            let gradient = map
                .zero_constant_linear_strain()
                .ok_or_else(|| invalid("supported original linear strain declaration overflows"))?;
            for i in 0..n {
                for j in 0..=i {
                    strain_matrix[i * n + j] += k * gradient[2 * i] * gradient[2 * j];
                    strain_matrix[j * n + i] = strain_matrix[i * n + j];
                }
            }
        }
        for (total, strain) in stiffness.iter_mut().zip(&strain_matrix) {
            *total += strain;
        }
        if stiffness.iter().any(|v| !v.is_finite()) {
            return Err(invalid("supported original quadratic stiffness overflows"));
        }
        Ok(SupportedLinearSplit {
            stiffness,
            strain_quadratic: QuadraticStorage::new(strain_matrix, n)
                .map_err(|error| ImpactError::Owner(error.to_string()))?,
            strains,
            modes: n,
        })
    }
}

impl SupportedLinearSplit {
    pub(in crate::render::plate::impact) fn mode_count(&self) -> usize {
        self.modes
    }

    pub(in crate::render::plate::impact) fn stiffness(&self) -> &[f64] {
        &self.stiffness
    }

    /// Nonquadratic work only. Original quadratic work is already in the
    /// complete mechanical Q. Both terms below are stable small remainders;
    /// no difference of large total energies recovers weak material coupling.
    pub(in crate::render::plate::impact) fn midpoint_nonlinear_remainder(
        &self,
        a: &[f64],
        b: &[f64],
        mid: &[f64],
    ) -> Option<f64> {
        let dim = 2 * self.modes;
        if a.len() != dim || b.len() != dim || mid.len() != dim {
            return None;
        }
        let mut result = Dd::from_f64(0.0);
        for (map, stiffness) in &self.strains {
            result = result
                + Dd::from_f64(*stiffness)
                    * Dd::from_f64(quadratic_strain_midpoint_energy_remainder(
                        map, 1.0, a, b, mid, [0.0; 3],
                    )?);
        }
        let (mut qa, mut qb, mut qm) = (
            [0.0; MAX_IMPACT_MODES],
            [0.0; MAX_IMPACT_MODES],
            [0.0; MAX_IMPACT_MODES],
        );
        for i in 0..self.modes {
            qa[i] = a[2 * i];
            qb[i] = b[2 * i];
            qm[i] = mid[2 * i];
        }
        result = result
            - Dd::from_f64(self.strain_quadratic.midpoint_energy_remainder(
                &qa[..self.modes],
                &qb[..self.modes],
                &qm[..self.modes],
            )?);
        let value = result.to_f64();
        value.is_finite().then_some(value)
    }
}

#[cfg(test)]
mod tests;
