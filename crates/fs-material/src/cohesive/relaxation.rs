//! Intact, small-strain adhesive layers with genuine relaxation memory.
//!
//! For original thickness h, the three strains are [gn,gt1,gt2]/h. Independent
//! normal and shear Maxwell spectra describe this particular constrained joint;
//! a bulk Young modulus, Poisson ratio, hardness or Kelvin viscosity does not
//! identify them. Bondline inertia and spatial lateral fields are neglected.
//! The two shear coordinates are isotropic. Compression remains
//! the same intact material. There is no fracture, contact opening, damage,
//! healing, temperature transfer or finite-rotation interpretation.
//!
//! The existing PronyState and exact PronyRamp own branch memory and physical
//! heat. Surface work/storage/heat are h times their volume densities. Prepared
//! ramp coefficients serve every point with the same law and duration.
//! Endpoint and work-mean tangents are distinct; both are exact for a straight
//! jump ramp. Held jumps relax, and the existing mechanical load line can use
//! the affine mean force to solve coupled motion and creep.
//!
//! Constitutive background: Abaqus, Time Domain Viscoelasticity, traction-
//! separation elasticity (independent normal/shear relaxation).
//! https://docs.software.vt.edu/abaqusv2026/English/SIMACAEMATRefMap/simamat-c-timevisco.htm

use super::{CohesiveError, finite, norm, positive};
use crate::visco::{GeneralizedMaxwell, PronyRamp, PronyState};

pub mod source;

/// Explicit local model-use bounds, not inferred experimental qualification.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ViscoelasticBondLimits {
    /// Euclidean separation norm, including compressive displacement [m].
    pub maximum_jump_m: f64,
    /// Euclidean rate of the declared straight jump ramp [m/s].
    pub maximum_rate_m_s: f64,
    pub maximum_elapsed_s: f64,
}

/// Independent effective normal/shear moduli [Pa] and original bond thickness.
/// The fixed frame and original bonded area belong to the mechanical caller.
#[derive(Debug, Clone, PartialEq)]
pub struct IntactViscoelasticBond {
    spectra: [GeneralizedMaxwell; 2],
    thickness: f64,
    limits: ViscoelasticBondLimits,
}

impl IntactViscoelasticBond {
    pub fn new(
        normal: GeneralizedMaxwell,
        shear: GeneralizedMaxwell,
        thickness_m: f64,
        limits: ViscoelasticBondLimits,
    ) -> Result<Self, CohesiveError> {
        for x in [
            thickness_m,
            limits.maximum_jump_m,
            limits.maximum_rate_m_s,
            limits.maximum_elapsed_s,
        ] {
            positive(x)?;
        }
        let strain_bound = positive(limits.maximum_jump_m / thickness_m)?;
        for spectrum in [&normal, &shear] {
            positive(spectrum.e_inf)?;
            if spectrum.terms.len() > 64 {
                return Err(CohesiveError::Invalid(
                    "at most 64 bond relaxation branches",
                ));
            }
            let mut instantaneous = spectrum.e_inf;
            for &(modulus, time) in &spectrum.terms {
                if !modulus.is_finite() || modulus < 0.0 {
                    return Err(CohesiveError::Invalid(
                        "nonnegative branch modulus required",
                    ));
                }
                positive(time)?;
                instantaneous = positive(instantaneous + modulus)?;
            }
            positive(instantaneous / thickness_m)?;
            positive((instantaneous * strain_bound) * strain_bound)?;
        }
        Ok(Self {
            spectra: [normal, shear],
            thickness: thickness_m,
            limits,
        })
    }
    pub const fn thickness_m(&self) -> f64 {
        self.thickness
    }
    pub const fn limits(&self) -> ViscoelasticBondLimits {
        self.limits
    }
    pub fn normal_spectrum(&self) -> &GeneralizedMaxwell {
        &self.spectra[0]
    }
    pub fn shear_spectrum(&self) -> &GeneralizedMaxwell {
        &self.spectra[1]
    }
    /// Original stress-free state; source identity stays with this borrowed law.
    pub fn bonded_point(&self) -> ViscoelasticBondPoint<'_> {
        ViscoelasticBondPoint {
            law: self,
            states: core::array::from_fn(|i| self.spectra[usize::from(i > 0)].state()),
            elapsed: 0.0,
        }
    }
    /// Allocate and evaluate exponentials once per law/duration; reuse for every
    /// trial and every point. Ordinary trial work is linear in branch count.
    pub fn prepare(&self, dt_s: f64) -> Result<ViscoelasticBondRamp<'_>, CohesiveError> {
        positive(dt_s)?;
        if dt_s > self.limits.maximum_elapsed_s {
            return Err(CohesiveError::Invalid("bond duration outside domain"));
        }
        let ramps = core::array::from_fn(|mode| {
            self.spectra[mode]
                .terms
                .iter()
                .map(|&(_, tau)| PronyRamp::new(dt_s / tau))
                .collect::<Vec<_>>()
        });
        let mut endpoint = [0.0; 2];
        let mut mean = [0.0; 2];
        for mode in 0..2 {
            endpoint[mode] = self.spectra[mode].e_inf;
            mean[mode] = 0.5 * self.spectra[mode].e_inf;
            for (&(e, _), ramp) in self.spectra[mode].terms.iter().zip(&ramps[mode]) {
                let (end, average, _) = ramp.apply(0.0, 1.0);
                endpoint[mode] = finite(endpoint[mode] + e * end)?;
                mean[mode] = finite(mean[mode] + e * average)?;
            }
            endpoint[mode] = positive(endpoint[mode] / self.thickness)?;
            mean[mode] = positive(mean[mode] / self.thickness)?;
        }
        Ok(ViscoelasticBondRamp {
            law: self,
            dt: dt_s,
            ramps,
            endpoint_stiffness: [endpoint[0], endpoint[1], endpoint[1]],
            mean_stiffness: [mean[0], mean[1], mean[1]],
        })
    }
}

/// Exact ramp data bound to one original law.
#[derive(Debug, Clone)]
pub struct ViscoelasticBondRamp<'a> {
    law: &'a IntactViscoelasticBond,
    dt: f64,
    ramps: [Vec<PronyRamp>; 2],
    endpoint_stiffness: [f64; 3],
    mean_stiffness: [f64; 3],
}
impl ViscoelasticBondRamp<'_> {
    pub const fn duration_s(&self) -> f64 {
        self.dt
    }
    pub const fn endpoint_stiffness_pa_m(&self) -> [f64; 3] {
        self.endpoint_stiffness
    }
    pub const fn mean_stiffness_pa_m(&self) -> [f64; 3] {
        self.mean_stiffness
    }
}

/// The original scalar Prony histories, never publicly replaceable or relabeled.
#[derive(Debug, Clone)]
pub struct ViscoelasticBondPoint<'a> {
    law: &'a IntactViscoelasticBond,
    states: [PronyState; 3],
    elapsed: f64,
}

/// Surface energy is per ORIGINAL bonded area; viscosity is physical heat.
#[derive(Debug, Clone)]
pub struct ViscoelasticBondStep<'a> {
    pub point: ViscoelasticBondPoint<'a>,
    pub traction_pa: [f64; 3],
    pub mean_traction_pa: [f64; 3],
    pub work_j_m2: f64,
    pub stored_change_j_m2: f64,
    pub viscous_dissipated_j_m2: f64,
    pub defect_j_m2: f64,
}
impl<'a> ViscoelasticBondPoint<'a> {
    pub const fn law(&self) -> &'a IntactViscoelasticBond {
        self.law
    }
    pub const fn elapsed_s(&self) -> f64 {
        self.elapsed
    }
    pub fn jump_m(&self) -> [f64; 3] {
        core::array::from_fn(|i| self.states[i].strain * self.law.thickness)
    }
    pub fn traction_pa(&self) -> [f64; 3] {
        core::array::from_fn(|i| self.states[i].stress)
    }
    pub fn stored_energy_j_m2(&self) -> f64 {
        self.states
            .iter()
            .enumerate()
            .map(|(i, s)| {
                let law = &self.law.spectra[usize::from(i > 0)];
                0.5 * (law.e_inf * s.strain) * s.strain
                    + law
                        .terms
                        .iter()
                        .zip(&s.q)
                        .map(|(&(e, _), &q)| 0.5 * (e * q) * q)
                        .sum::<f64>()
            })
            .sum::<f64>()
            * self.law.thickness
    }
    fn admit(&self, prepared: &ViscoelasticBondRamp<'_>) -> Result<f64, CohesiveError> {
        if !core::ptr::eq(self.law, prepared.law) {
            return Err(CohesiveError::Invalid(
                "relaxation ramp belongs to a different bond law",
            ));
        }
        let time = finite(self.elapsed + prepared.dt)?;
        if time <= self.elapsed || time > self.law.limits.maximum_elapsed_s {
            return Err(CohesiveError::Invalid("bond elapsed time outside domain"));
        }
        Ok(time)
    }
    /// Held-jump mean traction. Together with prepared.mean_stiffness this
    /// gives the exact affine mean map: held_mean + Kmean*(new_jump-old_jump).
    pub fn held_mean_traction(
        &self,
        prepared: &ViscoelasticBondRamp<'_>,
    ) -> Result<[f64; 3], CohesiveError> {
        self.admit(prepared)?;
        let mut result = [0.0; 3];
        for i in 0..3 {
            let mode = usize::from(i > 0);
            let law = &self.law.spectra[mode];
            let s = &self.states[i];
            result[i] = finite(law.e_inf * s.strain)?;
            for ((&(e, _), &q), ramp) in law.terms.iter().zip(&s.q).zip(&prepared.ramps[mode]) {
                result[i] = finite(result[i] + e * ramp.apply(q, 0.0).1)?;
            }
        }
        Ok(result)
    }
    /// Retryable immutable trial. No endpoint trapezoid or subtraction of
    /// accumulated work/heat inventories replaces the exact local increments.
    pub fn trial(
        &self,
        prepared: &ViscoelasticBondRamp<'_>,
        jump_m: [f64; 3],
        energy_tolerance_j_m2: f64,
    ) -> Result<ViscoelasticBondStep<'a>, CohesiveError> {
        positive(energy_tolerance_j_m2)?;
        let elapsed = self.admit(prepared)?;
        if norm(jump_m)? > self.law.limits.maximum_jump_m {
            return Err(CohesiveError::SeparationLimit);
        }
        let old_jump = self.jump_m();
        let delta = core::array::from_fn::<_, 3, _>(|i| jump_m[i] - old_jump[i]);
        if norm(delta.map(|x| x / prepared.dt))? > self.law.limits.maximum_rate_m_s {
            return Err(CohesiveError::Invalid("bond jump rate outside domain"));
        }
        let mut point = self.clone();
        point.elapsed = elapsed;
        let mut mean = [0.0; 3];
        let mut work = 0.0;
        let mut storage = 0.0;
        let mut heat = 0.0;
        for i in 0..3 {
            let mode = usize::from(i > 0);
            let law = &self.law.spectra[mode];
            let old = &self.states[i];
            let next = &mut point.states[i];
            let increment = finite(delta[i] / self.law.thickness)?;
            next.strain = finite(jump_m[i] / self.law.thickness)?;
            let mut end = finite(law.e_inf * next.strain)?;
            mean[i] = finite(law.e_inf * (old.strain + 0.5 * increment))?;
            let mut stored = finite((0.5 * law.e_inf * increment) * (old.strain + next.strain))?;
            let mut lost = 0.0;
            for j in 0..old.q.len() {
                let e = law.terms[j].0;
                let (q, average, loss) = prepared.ramps[mode][j].apply(old.q[j], increment);
                next.q[j] = finite(q)?;
                mean[i] = finite(mean[i] + e * average)?;
                end = finite(end + e * q)?;
                stored = finite(stored + (0.5 * e * (q - old.q[j])) * (q + old.q[j]))?;
                lost = finite(lost + e * loss)?;
            }
            if lost < 0.0 {
                return Err(CohesiveError::Invalid("negative bond relaxation heat"));
            }
            next.stress = end;
            next.work = finite(old.work + mean[i] * increment)?;
            next.dissipated = finite(old.dissipated + lost)?;
            work = finite(work + mean[i] * delta[i])?;
            storage = finite(storage + self.law.thickness * stored)?;
            heat = finite(heat + self.law.thickness * lost)?;
        }
        finite(point.stored_energy_j_m2())?;
        let defect = finite(work - storage - heat)?;
        if defect.abs() > energy_tolerance_j_m2 {
            return Err(CohesiveError::EnergyDefect(defect));
        }
        Ok(ViscoelasticBondStep {
            traction_pa: point.traction_pa(),
            point,
            mean_traction_pa: mean,
            work_j_m2: work,
            stored_change_j_m2: storage,
            viscous_dissipated_j_m2: heat,
            defect_j_m2: defect,
        })
    }
}

#[cfg(test)]
mod tests;
