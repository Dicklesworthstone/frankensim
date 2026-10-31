//! Exact original Prony ramp and work-conjugate quadratic-strain pullback.
use super::*;
use crate::render::schedule::force::coupled::contact::friction::nonlinear::history::ModalMaterialHistory;

pub(super) struct MeanStrainHistory {
    pub potential: StringPotential,
    pub length_m: f64,
    pub material: StrainRelaxation,
}
fn refused(what: &'static str) -> ModalCouplingError {
    ModalCouplingError::Invalid { what }
}

impl MeanStrainHistory {
    fn strain(&self, q: &[f64]) -> Result<f64, ModalCouplingError> {
        self.potential
            .observe(q)
            .map(|o| o.additional_strain)
            .map_err(|_| {
                refused("mean-strain history exceeds original string geometry or slope support")
            })
    }
    fn admit(&self, q: &[f64], memory: &[f64]) -> Result<f64, ModalCouplingError> {
        let observation = self
            .potential
            .observe(q)
            .map_err(|_| refused("mean-strain history lost original string coordinates"))?;
        let strain = observation.additional_strain;
        let mut stress = [0.0];
        self.material
            .stress(&[strain], memory, &mut stress)
            .map_err(|_| {
                refused("mean-strain history requires every original finite relaxed strain")
            })?;
        let tension = observation.tension_n + stress[0] / self.length_m;
        if !tension.is_finite() || tension <= 0.0 {
            return Err(refused(
                "hereditary core and winding stress loses positive string tension",
            ));
        }
        Ok(strain)
    }
    fn pull(&self, q: &[f64], stress_j: f64, out: &mut [f64]) -> Result<(), ModalCouplingError> {
        self.potential.mean_strain_gradient(q, out).map_err(|_| {
            refused("mean-strain material force requires its original coordinate gradient")
        })?;
        for value in out.iter_mut() {
            *value *= -stress_j;
        }
        if out.iter().any(|v| !v.is_finite()) {
            return Err(refused("mean-strain material force overflows"));
        }
        Ok(())
    }
    #[allow(clippy::too_many_arguments)]
    fn stage_inner(
        &self,
        before: &[f64],
        after: &[f64],
        memory: &[f64],
        dt: f64,
        candidate: &mut [f64],
        reaction: &mut [f64],
        channels: Option<&mut [f64]>,
    ) -> Result<f64, ModalCouplingError> {
        let a = [self.admit(before, memory)?];
        let b = [self.strain(after)?];
        let mut stress = [0.0];
        let loss = match channels {
            Some(channels) => self.material.stage_with_arm_dissipation(
                &a,
                &b,
                memory,
                dt,
                candidate,
                &mut stress,
                channels,
            ),
            None => self
                .material
                .stage(&a, &b, memory, dt, candidate, &mut stress),
        }
        .map_err(|_| refused("original Prony mean-strain trial is outside its finite domain"))?;
        self.admit(after, candidate)?;
        // Exact identity for epsilon=q^T E q/4:
        // grad(epsilon)(midpoint(q)) dot Delta(q) = Delta(epsilon).
        let midpoint: Vec<_> = before
            .iter()
            .zip(after)
            .map(|(&a, &b)| f64::midpoint(a, b))
            .collect();
        self.pull(&midpoint, stress[0], reaction)?;
        Ok(loss)
    }
}

impl ModalMaterialHistory for MeanStrainHistory {
    fn coordinate_count(&self) -> usize {
        self.potential.mode_count()
    }
    fn memory_count(&self) -> usize {
        self.material.memory_count()
    }
    fn stiffness_upper_bound_per_s2(
        &self,
        initial_memory: &[f64],
    ) -> Result<f64, ModalCouplingError> {
        if initial_memory.len() != self.memory_count()
            || initial_memory.iter().any(|v| !v.is_finite())
        {
            return Err(refused(
                "mean-strain stiffness needs all original initial memories",
            ));
        }
        // epsilon <= integral(|slope|^2)/2 <= slope_limit^2/2.
        // For the original positive quadratic strain with Hessian B,
        // |grad epsilon|^2 <= 2*||B||*epsilon. Exact Prony ramp memory stays
        // in the convex hull of its accepted value and endpoint strains.
        let strain_bound = 0.5 * self.potential.slope_limit().powi(2);
        let hessian = self.potential.mean_strain_hessian_norm();
        let bound: f64 = self
            .material
            .arms()
            .iter()
            .zip(initial_memory)
            .map(|(arm, z)| {
                arm.stiffness_j * hessian * (3.0 * strain_bound + z.abs().max(strain_bound))
            })
            .sum();
        if !bound.is_finite() || bound < 0.0 {
            return Err(refused("mean-strain hereditary stiffness overflows"));
        }
        Ok(bound)
    }
    fn relaxed_memory(&self, q: &[f64]) -> Result<Vec<f64>, ModalCouplingError> {
        self.material
            .relaxed_memory(&[self.strain(q)?])
            .map_err(|_| refused("initial mean-strain memory is unrepresentable"))
    }
    fn energy(&self, q: &[f64], memory: &[f64]) -> Result<f64, ModalCouplingError> {
        self.material
            .energy(&[self.admit(q, memory)?], memory)
            .map_err(|_| refused("mean-strain material storage is unrepresentable"))
    }
    fn reaction(
        &self,
        q: &[f64],
        memory: &[f64],
        out: &mut [f64],
    ) -> Result<(), ModalCouplingError> {
        let strain = [self.admit(q, memory)?];
        let mut stress = [0.0];
        self.material
            .stress(&strain, memory, &mut stress)
            .map_err(|_| refused("mean-strain material stress is unrepresentable"))?;
        self.pull(q, stress[0], out)
    }
    fn stage(
        &self,
        before: &[f64],
        after: &[f64],
        memory: &[f64],
        dt: f64,
        candidate: &mut [f64],
        reaction: &mut [f64],
    ) -> Result<f64, ModalCouplingError> {
        self.stage_inner(before, after, memory, dt, candidate, reaction, None)
    }
    fn dissipation_channel_count(&self) -> usize {
        self.material.arms().len()
    }
    fn stage_with_dissipation(
        &self,
        before: &[f64],
        after: &[f64],
        memory: &[f64],
        dt: f64,
        candidate: &mut [f64],
        reaction: &mut [f64],
        channels_j: &mut [f64],
    ) -> Result<f64, ModalCouplingError> {
        self.stage_inner(
            before,
            after,
            memory,
            dt,
            candidate,
            reaction,
            Some(channels_j),
        )
    }
}

#[cfg(test)]
mod tests;
