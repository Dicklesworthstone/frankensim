//! Source Maxwell stress on retained mean/cosine material strain coordinates.
//!
//! Strings and fiber ribbons use this same constitutive algebra. The caller
//! owns geometry, coordinate pullback, admitted material domain and publication.
//! A trial never mutates accepted memory or independently advances mechanics.
use super::{MAX_LONGITUDINAL_MODES, invalid};
use crate::NlModalError;
use fs_material::visco::PronyRamp;

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
    /// The material owner integrates a linear ramp between ENDPOINT strains
    /// exactly, returning the interval-mean work-conjugate stress and positive
    /// dashpot heat. Pulling this stress back at midpoint coordinates preserves
    /// work for quadratic strain geometry. Mechanical motion still requires
    /// time refinement; exact constitutive integration is not an exact joint
    /// nonlinear trajectory. Returns heat [J]; accepted memory is untouched.
    pub fn stage(
        &self,
        before: &[f64],
        after: &[f64],
        memory: &[f64],
        dt: f64,
        candidate: &mut [f64],
        stress: &mut [f64],
    ) -> Result<f64, NlModalError> {
        self.stage_inner(before, after, memory, dt, candidate, stress, None)
    }

    /// The same constitutive update, retaining each original arm's physical
    /// loss separately. Source geometry may route these channels to different
    /// finite bodies; their sum is the returned total. No second history update
    /// or residual-energy partition is used.
    #[allow(clippy::too_many_arguments)]
    pub fn stage_with_arm_dissipation(
        &self,
        before: &[f64],
        after: &[f64],
        memory: &[f64],
        dt: f64,
        candidate: &mut [f64],
        stress: &mut [f64],
        arm_dissipation_j: &mut [f64],
    ) -> Result<f64, NlModalError> {
        if arm_dissipation_j.len() != self.arms.len() {
            return Err(invalid(
                "one dissipation output per original Maxwell arm required",
            ));
        }
        self.stage_inner(
            before,
            after,
            memory,
            dt,
            candidate,
            stress,
            Some(arm_dissipation_j),
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn stage_inner(
        &self,
        before: &[f64],
        after: &[f64],
        memory: &[f64],
        dt: f64,
        candidate: &mut [f64],
        stress: &mut [f64],
        mut arm_dissipation_j: Option<&mut [f64]>,
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
        for (index, ((arm, z0), z1)) in self
            .arms
            .iter()
            .zip(memory.chunks_exact(self.width))
            .zip(candidate.chunks_exact_mut(self.width))
            .enumerate()
        {
            let rate = dt / arm.relaxation_time_s;
            // Every harmonic in this source arm shares one relaxation clock.
            // Prepare the exponential coefficients once per trial, then reuse
            // the positive scalar ramp without per-coordinate transcendental work.
            let ramp = PronyRamp::new(rate);
            let mut arm_loss = 0.0;
            for j in 0..self.width {
                let (next, effort, heat) = maxwell_step_prepared(
                    before[j],
                    after[j],
                    z0[j],
                    rate,
                    arm.stiffness_j * weight(j),
                    &ramp,
                )
                .ok_or_else(|| invalid("material strain-history step overflows"))?;
                z1[j] = next;
                stress[j] += effort;
                loss += heat;
                arm_loss += heat;
            }
            if let Some(channels) = &mut arm_dissipation_j {
                channels[index] = arm_loss;
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
/// Returns relaxed strain, interval-mean stress coefficient, and physical heat
/// for a linear ramp between the two supplied total strains.
///
/// The shared material owner retains exponential relaxation at every dt/tau,
/// including stiff arms; it does not use midpoint decay or extrapolate a
/// damping ratio. Memory remains a convex combination of its initial value
/// and the endpoint strains. Zero/infinite rate are the elastic/fully-relaxed
/// limits, including underflow/overflow of positive finite dt/tau.
pub fn maxwell_step(
    before: f64,
    after: f64,
    memory: f64,
    rate: f64,
    stiffness: f64,
) -> Option<(f64, f64, f64)> {
    if rate.is_nan() || rate < 0.0 {
        return None;
    }
    let ramp = PronyRamp::new(rate);
    maxwell_step_prepared(before, after, memory, rate, stiffness, &ramp)
}

fn maxwell_step_prepared(
    before: f64,
    after: f64,
    memory: f64,
    rate: f64,
    stiffness: f64,
    ramp: &PronyRamp,
) -> Option<(f64, f64, f64)> {
    if [before, after, memory, stiffness]
        .iter()
        .any(|v| !v.is_finite())
        || stiffness <= 0.0
    {
        return None;
    }
    let elastic = before - memory;
    let increment = after - before;
    if !elastic.is_finite() || !increment.is_finite() {
        return None;
    }
    let (endpoint, mean, heat) = ramp.apply(elastic, increment);
    // z' = elastic_strain/tau. Integrating its mean retains tiny memory
    // changes without subtracting two large endpoint strains. For stiff
    // histories use the decayed elastic strain, avoiding infinity * zero.
    let next = if rate <= 1.0 {
        memory + rate * mean
    } else {
        after - endpoint
    };
    let force = stiffness * mean;
    let loss = stiffness * heat;
    ([next, force, loss].iter().all(|v| v.is_finite()) && loss >= 0.0)
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
                // Independent continuum solution q(t)=slope+(q0-slope)e^-t/tau,
                // integrated by composite Simpson quadrature, not the owner.
                let slope = (eb - ea) / rate;
                let transient = ea - z - slope;
                let endpoint = slope + transient * (-rate).exp();
                let (mut mean, mut squared) = (0.0, 0.0);
                for sample in 0..=128 {
                    let t = sample as f64 / 128.0;
                    let q = slope + transient * (-rate * t).exp();
                    let w = if sample == 0 || sample == 128 {
                        1.0
                    } else if sample % 2 == 0 {
                        2.0
                    } else {
                        4.0
                    };
                    mean += w * q / 384.0;
                    squared += w * q * q / 384.0;
                }
                initial += 0.5 * arm.stiffness_j * (ea - z).powi(2) / 512.0;
                final_energy += 0.5 * arm.stiffness_j * endpoint.powi(2) / 512.0;
                heat += arm.stiffness_j * rate * squared / 512.0;
                for (f, s) in spatial_stress.iter_mut().zip(shapes) {
                    *f += arm.stiffness_j * mean * s / 512.0;
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
    fn g1_strain_memory_recovers_source_relaxation_time_without_clock_ringing() {
        for tau in [0.001, 1e-6] {
            let law = StrainRelaxation::new(
                1,
                vec![StrainRelaxationArm {
                    stiffness_j: 1700.0,
                    relaxation_time_s: tau,
                }],
            )
            .unwrap();
            let epsilon = [0.002, 0.0003];
            let original = [0.0003, -0.0001];
            for steps in [1, 4, 16, 64] {
                let dt = 0.002 / steps as f64;
                let mut memory = original.to_vec();
                let mut candidate = memory.clone();
                let mut stress = [0.0; 2];
                let mut heat = 0.0;
                for _ in 0..steps {
                    let loss = law
                        .stage(&epsilon, &epsilon, &memory, dt, &mut candidate, &mut stress)
                        .unwrap();
                    assert!(loss >= 0.0);
                    for ((&next, &old), &target) in candidate.iter().zip(&memory).zip(&epsilon) {
                        assert!(next >= old && next <= target, "relaxation cannot overshoot");
                    }
                    core::mem::swap(&mut memory, &mut candidate);
                    heat += loss;
                }
                for ((&actual, &initial), &target) in memory.iter().zip(&original).zip(&epsilon) {
                    let exact = target + (initial - target) * (-0.002 / tau).exp();
                    assert!((actual - exact).abs() < 8e-18);
                }
                assert!(
                    (law.energy(&epsilon, &memory).unwrap() + heat
                        - law.energy(&epsilon, &original).unwrap())
                    .abs()
                        < 1e-14
                );
            }
        }
    }

    #[test]
    fn g3_linear_strain_partition_preserves_memory_work_and_positive_heat() {
        let (before, after, memory, stiffness) = (0.004, -0.001, -0.002, 1300.0);
        for rate in [1e-12, 1e-5, 0.7, 2.0, 1000.0] {
            let (end, mean, heat) = maxwell_step(before, after, memory, rate, stiffness).unwrap();
            let work = mean * (after - before);
            for steps in [2, 16, 64] {
                let mut z = memory;
                let (mut split_work, mut split_heat) = (0.0, 0.0);
                for n in 0..steps {
                    let a = before + (after - before) * n as f64 / steps as f64;
                    let b = before + (after - before) * (n + 1) as f64 / steps as f64;
                    let (next, stress, loss) =
                        maxwell_step(a, b, z, rate / steps as f64, stiffness).unwrap();
                    assert!(loss >= 0.0);
                    split_work += stress * (b - a);
                    split_heat += loss;
                    z = next;
                }
                assert!((z - end).abs() < 5e-16, "rate={rate}, steps={steps}");
                assert!((split_work - work).abs() < 2e-14);
                assert!((split_heat - heat).abs() < 2e-14);
            }
            let old = 0.5 * stiffness * (before - memory).powi(2);
            let new = 0.5 * stiffness * (after - end).powi(2);
            assert!((new - old + heat - work).abs() < 2e-14);
        }
    }

    #[test]
    fn g0_strain_trials_preserve_accepted_memory_and_reject_invalid_clock() {
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
        // A stiff source arm is now integrated directly, without midpoint
        // ringing or the old dt/tau <= 2 refusal.
        assert!(
            law.stage(&a, &b, &memory, 0.021, &mut repeat, &mut again)
                .is_ok()
        );
        for dt in [0.0, -0.1, f64::INFINITY, f64::NAN] {
            assert!(
                law.stage(&a, &b, &memory, dt, &mut repeat, &mut again)
                    .is_err()
            );
        }
        assert_eq!(memory, original);
    }

    #[test]
    fn g1_original_arm_heat_matches_spatial_integrals_and_preserves_the_joint_update() {
        let law = StrainRelaxation::new(
            2,
            vec![
                StrainRelaxationArm {
                    stiffness_j: 1200.0,
                    relaxation_time_s: 0.001,
                },
                StrainRelaxationArm {
                    stiffness_j: 310.0,
                    relaxation_time_s: 0.0004,
                },
            ],
        )
        .unwrap();
        let before = [0.004, -0.002, 0.001];
        let after = [0.002, 0.003, -0.0005];
        let memory = [0.0006, 0.0002, -0.0003, -0.0008, 0.0005, 0.0001];
        let accepted = memory;
        let dt = 0.0002;
        let (mut ordinary, mut candidate) = ([0.0; 6], [0.0; 6]);
        let (mut ordinary_stress, mut stress) = ([0.0; 3], [0.0; 3]);
        let total = law
            .stage(
                &before,
                &after,
                &memory,
                dt,
                &mut ordinary,
                &mut ordinary_stress,
            )
            .unwrap();
        let mut channels = [f64::NAN; 2];
        assert_eq!(
            total,
            law.stage_with_arm_dissipation(
                &before,
                &after,
                &memory,
                dt,
                &mut candidate,
                &mut stress,
                &mut channels
            )
            .unwrap()
        );
        assert_eq!(candidate, ordinary);
        assert_eq!(stress, ordinary_stress);
        let field = |coeff: &[f64], x: f64| {
            coeff[0]
                + coeff[1..]
                    .iter()
                    .enumerate()
                    .map(|(j, v)| v * ((j + 1) as f64 * core::f64::consts::PI * x).cos())
                    .sum::<f64>()
        };
        let mut spatial = [0.0; 2];
        for point in 0..512 {
            let x = (point as f64 + 0.5) / 512.0;
            for (arm, source) in law.arms().iter().enumerate() {
                let rate = dt / source.relaxation_time_s;
                let ea = field(&before, x);
                let eb = field(&after, x);
                let old = field(&memory[arm * 3..(arm + 1) * 3], x);
                // Independent continuous Maxwell solution, integrated in time
                // at each spatial quadrature point. The tight source-channel
                // comparison uses a resolved 2048-panel Simpson rule.
                let slope = (eb - ea) / rate;
                let transient = ea - old - slope;
                let mut squared = 0.0;
                for sample in 0..=2048 {
                    let t = sample as f64 / 2048.0;
                    let elastic = slope + transient * (-rate * t).exp();
                    let w = if sample == 0 || sample == 2048 {
                        1.0
                    } else if sample % 2 == 0 {
                        2.0
                    } else {
                        4.0
                    };
                    squared += w * elastic * elastic / 6144.0;
                }
                spatial[arm] += source.stiffness_j * rate * squared / 512.0;
            }
        }
        for (observed, expected) in channels.iter().zip(spatial) {
            assert!(
                (*observed - expected).abs() < 2e-14 * expected,
                "arm heat {observed:e} vs {expected:e}, relative error {}",
                (*observed - expected).abs() / expected
            );
        }
        assert!(channels[0] > 0.0 && channels[1] > 0.0 && channels[0] != channels[1]);
        assert!((channels.iter().sum::<f64>() - total).abs() < 2e-15 * total);
        assert_eq!(memory, accepted);
        assert!(
            law.stage_with_arm_dissipation(
                &before,
                &after,
                &memory,
                0.0,
                &mut candidate,
                &mut stress,
                &mut channels
            )
            .is_err()
        );
        assert_eq!(
            law.stage_with_arm_dissipation(
                &before,
                &after,
                &memory,
                dt,
                &mut candidate,
                &mut stress,
                &mut channels
            )
            .unwrap(),
            total
        );
        assert_eq!(candidate, ordinary);
        assert_eq!(stress, ordinary_stress);
        assert_eq!(memory, accepted);
        for (observed, expected) in channels.iter().zip(spatial) {
            assert!(
                (*observed - expected).abs() < 2e-14 * expected,
                "arm heat {observed:e} vs {expected:e}, relative error {}",
                (*observed - expected).abs() / expected
            );
        }
        assert!(
            law.stage_with_arm_dissipation(
                &before,
                &after,
                &memory,
                dt,
                &mut candidate,
                &mut stress,
                &mut channels[..1]
            )
            .is_err()
        );
    }
}
