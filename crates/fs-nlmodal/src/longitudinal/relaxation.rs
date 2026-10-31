//! Source Maxwell stress on retained mean/cosine material strain coordinates.
//!
//! Strings and fiber ribbons use this same constitutive algebra. The caller
//! owns geometry, coordinate pullback, admitted material domain and publication.
//! A trial never mutates accepted memory or independently advances mechanics.
use super::{MAX_LONGITUDINAL_MODES, invalid};
use crate::NlModalError;

/// One original source arm; no averaging over layers or relaxation times.
#[derive(Clone, Debug)]
pub struct StrainRelaxationArm {
    /// Delta(EA) times the original material span [J].
    pub stiffness_j: f64,
    pub relaxation_time_s: f64,
}

/// Passive Maxwell storage over the SAME retained material strain coordinates.
/// Memory is dimensionless relaxed strain, arm-major, then mean/cosine order.
#[derive(Clone, Debug)]
pub struct StrainRelaxation {
    width: usize,
    arms: Vec<StrainRelaxationArm>,
}

impl StrainRelaxation {
    /// Retain the declared axial harmonics plus the geometric mean strain.
    pub fn new(harmonics: usize, arms: Vec<StrainRelaxationArm>) -> Result<Self, NlModalError> {
        if harmonics > MAX_LONGITUDINAL_MODES
            || arms.is_empty()
            || arms.iter().any(|arm| {
                !arm.stiffness_j.is_finite()
                    || arm.stiffness_j <= 0.0
                    || !arm.relaxation_time_s.is_finite()
                    || arm.relaxation_time_s <= 0.0
            })
            || arms.len().checked_mul(harmonics + 1).is_none()
        {
            return Err(invalid(
                "strain memory needs retained harmonics and positive finite source arms",
            ));
        }
        Ok(Self {
            width: harmonics + 1,
            arms,
        })
    }

    pub const fn width(&self) -> usize {
        self.width
    }
    pub fn arms(&self) -> &[StrainRelaxationArm] {
        &self.arms
    }
    pub fn memory_count(&self) -> usize {
        self.width * self.arms.len()
    }

    fn check(&self, strain: &[f64], memory: &[f64]) -> Result<(), NlModalError> {
        if strain.len() != self.width
            || memory.len() != self.memory_count()
            || strain.iter().chain(memory).any(|v| !v.is_finite())
        {
            return Err(invalid(
                "complete finite mean/cosine strain and material memory required",
            ));
        }
        Ok(())
    }

    /// Relaxed at the actual strain, rather than assuming zero displacement.
    pub fn relaxed_memory(&self, strain: &[f64]) -> Result<Vec<f64>, NlModalError> {
        if strain.len() != self.width || strain.iter().any(|v| !v.is_finite()) {
            return Err(invalid("complete finite initial material strain required"));
        }
        Ok(strain.repeat(self.arms.len()))
    }

    /// Mean storage is C r_0²/2 and each retained cosine stores C r_j²/4.
    pub fn energy(&self, strain: &[f64], memory: &[f64]) -> Result<f64, NlModalError> {
        self.check(strain, memory)?;
        let mut energy = 0.0;
        for (arm, z) in self.arms.iter().zip(memory.chunks_exact(self.width)) {
            for (j, (&epsilon, &z)) in strain.iter().zip(z).enumerate() {
                energy += 0.5 * arm.stiffness_j * weight(j) * (epsilon - z).powi(2);
            }
        }
        if !energy.is_finite() {
            return Err(invalid("material strain storage overflows"));
        }
        Ok(energy)
    }

    /// Instantaneous work-conjugate stress coefficients [J]. The mean has
    /// weight one and each cosine weight one half, from spatial orthogonality.
    /// Physical axial-force coefficients are out[0]/L and 2*out[j]/L (j>0).
    pub fn stress(
        &self,
        strain: &[f64],
        memory: &[f64],
        out: &mut [f64],
    ) -> Result<(), NlModalError> {
        self.check(strain, memory)?;
        if out.len() != self.width {
            return Err(invalid("material strain stress shape"));
        }
        out.fill(0.0);
        for (arm, z) in self.arms.iter().zip(memory.chunks_exact(self.width)) {
            for (j, ((force, epsilon), z)) in out.iter_mut().zip(strain).zip(z).enumerate() {
                *force += arm.stiffness_j * weight(j) * (epsilon - z);
            }
        }
        if out.iter().any(|v| !v.is_finite()) {
            return Err(invalid("material strain stress overflows"));
        }
        Ok(())
    }

    /// Propose memory and stress for one interval of the original joint solve.
    /// The midpoint of ENDPOINT strains preserves work for quadratic geometry.
    /// Returns physical dissipation [J]; accepted memory remains untouched.
    pub fn stage(
        &self,
        before: &[f64],
        after: &[f64],
        memory: &[f64],
        dt: f64,
        candidate: &mut [f64],
        stress: &mut [f64],
    ) -> Result<f64, NlModalError> {
        self.check(before, memory)?;
        self.check(after, memory)?;
        if candidate.len() != memory.len()
            || stress.len() != self.width
            || !dt.is_finite()
            || dt <= 0.0
        {
            return Err(invalid(
                "material history trial needs complete state and positive physical duration",
            ));
        }
        stress.fill(0.0);
        let mut loss = 0.0;
        for ((arm, z0), z1) in self
            .arms
            .iter()
            .zip(memory.chunks_exact(self.width))
            .zip(candidate.chunks_exact_mut(self.width))
        {
            let rate = dt / arm.relaxation_time_s;
            for j in 0..self.width {
                let (next, effort, heat) = maxwell_step(
                    before[j],
                    after[j],
                    z0[j],
                    rate,
                    arm.stiffness_j * weight(j),
                )
                .ok_or_else(|| {
                    invalid("material clock is unresolved or strain-history step overflows")
                })?;
                z1[j] = next;
                stress[j] += effort;
                loss += heat;
            }
        }
        if !loss.is_finite() || stress.iter().any(|v| !v.is_finite()) {
            return Err(invalid("material strain dissipation or stress overflows"));
        }
        Ok(loss)
    }
}

fn weight(index: usize) -> f64 {
    if index == 0 { 1.0 } else { 0.5 }
}

/// One Maxwell strain coordinate, eliminated inside the existing joint trial.
/// Returns relaxed strain, midpoint stress coefficient, and heat. dt/tau <= 2
/// keeps the update a convex combination of accepted memory and endpoint strain.
pub fn maxwell_step(
    before: f64,
    after: f64,
    memory: f64,
    rate: f64,
    stiffness: f64,
) -> Option<(f64, f64, f64)> {
    if [before, after, memory, rate, stiffness]
        .iter()
        .any(|v| !v.is_finite())
        || rate <= 0.0
        || rate > 2.0
        || stiffness <= 0.0
    {
        return None;
    }
    let residual = (f64::midpoint(before, after) - memory) / (1.0 + 0.5 * rate);
    let next = memory + rate * residual;
    let force = stiffness * residual;
    let loss = stiffness * rate * residual * residual;
    [next, force, loss]
        .iter()
        .all(|v| v.is_finite())
        .then_some((next, force, loss))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn g1_strain_memory_matches_independent_spatial_storage_force_and_heat() {
        let law = StrainRelaxation::new(
            2,
            vec![
                StrainRelaxationArm {
                    stiffness_j: 1200.0,
                    relaxation_time_s: 0.001,
                },
                StrainRelaxationArm {
                    stiffness_j: 450.0,
                    relaxation_time_s: 0.003,
                },
            ],
        )
        .unwrap();
        let a = [0.002, 0.0003, -0.0001];
        let b = [0.0021, 0.00025, -0.00008];
        let memory = [0.0018, 0.0001, -0.00005, 0.0019, 0.0002, -0.00004];
        let mut next = vec![0.0; memory.len()];
        let mut stress = [0.0; 3];
        let dt = 0.0002;
        let loss = law
            .stage(&a, &b, &memory, dt, &mut next, &mut stress)
            .unwrap();
        let (mut initial, mut final_energy, mut heat) = (0.0, 0.0, 0.0);
        let mut spatial_stress = [0.0; 3];
        for cell in 0..512 {
            let x = (cell as f64 + 0.5) / 512.0;
            let shapes = [
                1.0,
                (core::f64::consts::PI * x).cos(),
                (2.0 * core::f64::consts::PI * x).cos(),
            ];
            let field = |q: &[f64]| q.iter().zip(shapes).map(|(q, s)| q * s).sum::<f64>();
            for (arm, z) in law.arms().iter().zip(memory.chunks_exact(3)) {
                let (ea, eb, z) = (field(&a), field(&b), field(z));
                let rate = dt / arm.relaxation_time_s;
                let z1 = ((1.0 - rate / 2.0) * z + rate * (ea + eb) / 2.0) / (1.0 + rate / 2.0);
                let residual = (ea + eb - z - z1) / 2.0;
                initial += 0.5 * arm.stiffness_j * (ea - z).powi(2) / 512.0;
                final_energy += 0.5 * arm.stiffness_j * (eb - z1).powi(2) / 512.0;
                heat += arm.stiffness_j * rate * residual.powi(2) / 512.0;
                for (f, s) in spatial_stress.iter_mut().zip(shapes) {
                    *f += arm.stiffness_j * residual * s / 512.0;
                }
            }
        }
        assert!((law.energy(&a, &memory).unwrap() - initial).abs() < 1e-15);
        assert!((law.energy(&b, &next).unwrap() - final_energy).abs() < 1e-15);
        assert!((loss - heat).abs() < 1e-15);
        for (a, b) in stress.iter().zip(spatial_stress) {
            assert!((a - b).abs() < 1e-13);
        }
        let work: f64 = stress
            .iter()
            .zip(a.iter().zip(b))
            .map(|(s, (a, b))| s * (b - a))
            .sum();
        assert!((final_energy - initial + heat - work).abs() < 1e-15);
    }

    #[test]
    fn g1_strain_memory_recovers_source_relaxation_time() {
        let law = StrainRelaxation::new(
            1,
            vec![StrainRelaxationArm {
                stiffness_j: 1700.0,
                relaxation_time_s: 0.001,
            }],
        )
        .unwrap();
        let epsilon = [0.002, 0.0003];
        let original = [0.0003, -0.0001];
        let mut errors = Vec::new();
        for steps in [16, 32, 64] {
            let dt = 0.002 / steps as f64;
            let mut memory = original.to_vec();
            let mut candidate = memory.clone();
            let mut stress = [0.0; 2];
            let mut heat = 0.0;
            for _ in 0..steps {
                heat += law
                    .stage(&epsilon, &epsilon, &memory, dt, &mut candidate, &mut stress)
                    .unwrap();
                core::mem::swap(&mut memory, &mut candidate);
            }
            let exact = epsilon[0] + (original[0] - epsilon[0]) * (-2.0_f64).exp();
            errors.push((memory[0] - exact).abs());
            assert!(
                (law.energy(&epsilon, &memory).unwrap() + heat
                    - law.energy(&epsilon, &original).unwrap())
                .abs()
                    < 1e-14
            );
        }
        for pair in errors.windows(2) {
            assert!((pair[0] / pair[1] - 4.0).abs() < 0.02);
        }
    }

    #[test]
    fn g0_strain_trials_preserve_accepted_memory_and_reject_unresolved_clock() {
        let law = StrainRelaxation::new(
            1,
            vec![StrainRelaxationArm {
                stiffness_j: 100.0,
                relaxation_time_s: 0.001,
            }],
        )
        .unwrap();
        let a = [0.001, 0.0001];
        let b = [0.0012, 0.0002];
        let memory = law.relaxed_memory(&a).unwrap();
        let original = memory.clone();
        let mut next = vec![0.0; memory.len()];
        let mut force = vec![0.0; a.len()];
        let loss = law
            .stage(&a, &b, &memory, 0.0002, &mut next, &mut force)
            .unwrap();
        let mut repeat = next.clone();
        let mut again = force.clone();
        assert_eq!(
            loss,
            law.stage(&a, &b, &memory, 0.0002, &mut repeat, &mut again)
                .unwrap()
        );
        assert_eq!(next, repeat);
        assert_eq!(force, again);
        assert_eq!(memory, original);
        assert!(
            law.stage(&a, &b, &memory, 0.0021, &mut repeat, &mut again)
                .is_err()
        );
        assert_eq!(memory, original);
    }
}
