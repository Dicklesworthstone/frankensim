//! Original fiber geometry pulling back the shared Maxwell strain material.
use super::*;
use crate::render::schedule::force::coupled::contact::friction::nonlinear::history::ModalMaterialHistory;
use fs_material::section::relaxation::ConcentricPronySection;
use fs_nlmodal::longitudinal::relaxation::{StrainRelaxation, StrainRelaxationArm};

#[derive(Clone, Debug)]
pub(in crate::fiber_ribbon) struct FiberStrainHistory {
    strain: FiberStrain,
    length_m: f64,
    material: StrainRelaxation,
}

impl FiberStrainHistory {
    pub(in crate::fiber_ribbon) fn from_source(
        strain: &FiberStrain,
        source: &ConcentricPronySection,
        length_m: f64,
        share: f64,
    ) -> Result<Option<Self>, FiberRibbonError> {
        let mut arms = Vec::new();
        for (spectrum, geometry) in source
            .spectra()
            .iter()
            .zip(source.elastic().layer_integrals())
        {
            for &(delta, relaxation_time_s) in &spectrum.terms {
                if delta == 0.0 {
                    continue;
                }
                arms.push(StrainRelaxationArm {
                    stiffness_j: delta * geometry.area_m2() * share * length_m,
                    relaxation_time_s,
                });
            }
        }
        if arms.is_empty() {
            return Ok(None);
        }
        let material = StrainRelaxation::new(strain.model.axial_modes().len(), arms)
            .map_err(FiberRibbonError::Strain)?;
        Ok(Some(Self {
            strain: strain.clone(),
            length_m,
            material,
        }))
    }

    pub(in crate::fiber_ribbon) fn stiffness_bound(&self, memory: Option<&[f64]>) -> f64 {
        let axial = 0.5
            * self
                .strain
                .model
                .axial_slope_scales()
                .iter()
                .map(|v| v * v)
                .sum::<f64>();
        let normal = 0.5 * self.strain.slope_scales.iter().map(|v| v * v).sum::<f64>();
        let epsilon = self.strain.limits.maximum_strain;
        self.material
            .arms()
            .iter()
            .enumerate()
            .map(|(i, arm)| {
                // dt/tau <= 2 retains this continuous-span envelope for every
                // accepted future memory, including geometric-stress stiffness.
                let z = memory.map_or(epsilon, |m| {
                    m[i * self.material.width()..(i + 1) * self.material.width()]
                        .iter()
                        .map(|v| v.abs())
                        .sum::<f64>()
                        .max(epsilon)
                });
                arm.stiffness_j
                    * (axial + normal * (self.strain.limits.maximum_slope.powi(2) + epsilon + z))
            })
            .sum()
    }

    fn strain(&self, q: &[f64]) -> Result<Vec<f64>, ModalCouplingError> {
        self.strain.validate_coordinates(q)?;
        let mut strain = vec![0.0; self.material.width()];
        strain[0] = self
            .strain
            .model
            .geometric_strain(&self.strain.slopes(q), &mut strain[1..])
            .map_err(|_| invalid())?;
        for (j, &scale) in self.strain.model.axial_slope_scales().iter().enumerate() {
            strain[j + 1] += scale * q[j];
        }
        Ok(strain)
    }

    fn validate_memory(&self, strain: &[f64], memory: &[f64]) -> Result<(), ModalCouplingError> {
        let mut tension = vec![0.0; self.material.width()];
        self.material
            .stress(strain, memory, &mut tension)
            .map_err(|_| invalid())?;
        for (j, (force, epsilon)) in tension.iter_mut().zip(strain).enumerate() {
            *force = *force / self.length_m * if j == 0 { 1.0 } else { 2.0 }
                + self.strain.axial_rigidity_n * epsilon;
        }
        tension[0] += self.strain.tension_n;
        if tension.iter().any(|v| !v.is_finite())
            || tension[0] - tension[1..].iter().map(|v| v.abs()).sum::<f64>() <= 0.0
        {
            return Err(invalid());
        }
        Ok(())
    }

    fn pull(&self, q: &[f64], stress: &[f64], out: &mut [f64]) -> Result<(), ModalCouplingError> {
        if out.len() != self.coordinate_count() {
            return Err(invalid());
        }
        let slopes = self.strain.slopes(q);
        let mut gradient = vec![0.0; slopes.len()];
        self.strain
            .model
            .geometric_strain_gradient(&slopes, stress[0], &stress[1..], &mut gradient)
            .map_err(|_| invalid())?;
        let (axial, normal) = out.split_at_mut(self.strain.model.axial_modes().len());
        for ((force, &scale), &stress) in axial
            .iter_mut()
            .zip(self.strain.model.axial_slope_scales())
            .zip(&stress[1..])
        {
            *force = -scale * stress;
        }
        for ((force, &scale), &gradient) in normal
            .iter_mut()
            .zip(&self.strain.slope_scales)
            .zip(&gradient[1..])
        {
            *force = -scale * gradient;
        }
        if out.iter().any(|v| !v.is_finite()) {
            return Err(invalid());
        }
        Ok(())
    }
}

impl ModalMaterialHistory for FiberStrainHistory {
    fn coordinate_count(&self) -> usize {
        self.strain.coordinate_count()
    }
    fn memory_count(&self) -> usize {
        self.material.memory_count()
    }
    fn stiffness_upper_bound_per_s2(&self, memory: &[f64]) -> Result<f64, ModalCouplingError> {
        if memory.len() != self.memory_count() || memory.iter().any(|v| !v.is_finite()) {
            return Err(invalid());
        }
        let bound = self.stiffness_bound(Some(memory));
        if !bound.is_finite() {
            return Err(invalid());
        }
        Ok(bound)
    }
    fn relaxed_memory(&self, q: &[f64]) -> Result<Vec<f64>, ModalCouplingError> {
        self.material
            .relaxed_memory(&self.strain(q)?)
            .map_err(|_| invalid())
    }
    fn energy(&self, q: &[f64], memory: &[f64]) -> Result<f64, ModalCouplingError> {
        let strain = self.strain(q)?;
        self.validate_memory(&strain, memory)?;
        self.material.energy(&strain, memory).map_err(|_| invalid())
    }
    fn reaction(
        &self,
        q: &[f64],
        memory: &[f64],
        out: &mut [f64],
    ) -> Result<(), ModalCouplingError> {
        let strain = self.strain(q)?;
        self.validate_memory(&strain, memory)?;
        let mut stress = vec![0.0; self.material.width()];
        self.material
            .stress(&strain, memory, &mut stress)
            .map_err(|_| invalid())?;
        self.pull(q, &stress, out)
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
        let a = self.strain(before)?;
        let b = self.strain(after)?;
        self.validate_memory(&a, memory)?;
        let mut stress = vec![0.0; self.material.width()];
        let loss = self
            .material
            .stage(&a, &b, memory, dt, candidate, &mut stress)
            .map_err(|_| invalid())?;
        self.validate_memory(&b, candidate)?;
        let midpoint: Vec<_> = before
            .iter()
            .zip(after)
            .map(|(&a, &b)| f64::midpoint(a, b))
            .collect();
        self.pull(&midpoint, &stress, reaction)?;
        Ok(loss)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn law(share: f64) -> FiberStrainHistory {
        let length = 0.6;
        let mu = 0.0005 * share;
        let ea = 2000.0 * share;
        let tension = 10.0 * share;
        let model = LongitudinalString::new(
            KcStringParams {
                length,
                lin_density: mu,
                tension,
                ea,
            },
            2,
            1,
            &[1, 2, 3, 4],
        )
        .unwrap();
        FiberStrainHistory {
            strain: FiberStrain {
                model,
                slope_scales: (1..=2)
                    .map(|n| {
                        (2.0 / (mu * length)).sqrt() * n as f64 * core::f64::consts::PI / length
                    })
                    .collect(),
                tension_n: tension,
                axial_rigidity_n: ea,
                viscous_coefficient_j_s: 0.0,
                limits: FiberStrainLimits {
                    maximum_slope: 0.2,
                    maximum_strain: 0.01,
                },
                config: NonlinearModalConfig {
                    max_iterations: 128,
                    maximum_generalized_force: 1e6,
                    force_absolute_tolerance: 1e-10,
                    force_relative_tolerance: 1e-9,
                },
            },
            length_m: length,
            material: StrainRelaxation::new(
                4,
                vec![
                    StrainRelaxationArm {
                        stiffness_j: 250.0 * share,
                        relaxation_time_s: 0.001,
                    },
                    StrainRelaxationArm {
                        stiffness_j: 150.0 * share,
                        relaxation_time_s: 0.003,
                    },
                ],
            )
            .unwrap(),
        }
    }

    fn spatial(q: &[f64], fraction: f64, share: f64) -> (f64, [f64; 6]) {
        let scale = (2.0 / (0.0005 * share * 0.6)).sqrt() * core::f64::consts::PI / 0.6;
        let shape =
            |n: usize| scale * n as f64 * (n as f64 * core::f64::consts::PI * fraction).cos();
        let mut gradient = [0.0; 6];
        let mut strain = 0.0;
        for mode in 0..4 {
            gradient[mode] = shape(mode + 1);
            strain += gradient[mode] * q[mode];
        }
        let slope = shape(1) * q[4] + shape(2) * q[5];
        strain += 0.5 * slope * slope;
        for mode in 0..2 {
            gradient[4 + mode] = shape(mode + 1) * slope;
        }
        (strain, gradient)
    }

    #[test]
    fn g1_fiber_history_matches_independent_spatial_work_with_normal_axial_coupling() {
        let law = law(1.0);
        let a = [1e-8, -2e-8, 0.0, 1e-8, 2e-5, -1e-5];
        let delta = [2e-9, -1e-9, 1e-9, -2e-9, 2e-6, -1e-6];
        let b: Vec<_> = a.iter().zip(delta).map(|(a, d)| a + d).collect();
        let midpoint: Vec<_> = a
            .iter()
            .zip(&b)
            .map(|(&a, &b)| f64::midpoint(a, b))
            .collect();
        let mut memory = law.relaxed_memory(&a).unwrap();
        memory[0] -= 2e-6;
        memory[6] += 1e-7;
        let original = memory.clone();
        let mut next = vec![0.0; memory.len()];
        let mut force = [0.0; 6];
        let dt = 0.0001;
        let loss = law
            .stage(&a, &b, &memory, dt, &mut next, &mut force)
            .unwrap();
        let (mut heat, mut initial, mut stored) = (0.0, 0.0, 0.0);
        let mut spatial_force = [0.0; 6];
        for cell in 0..512 {
            let x = (cell as f64 + 0.5) / 512.0;
            let ea = spatial(&a, x, 1.0).0;
            let eb = spatial(&b, x, 1.0).0;
            let gradient = spatial(&midpoint, x, 1.0).1;
            for (arm, z) in law.material.arms().iter().zip(memory.chunks_exact(5)) {
                let z = z[0]
                    + z[1..]
                        .iter()
                        .enumerate()
                        .map(|(j, z)| z * ((j + 1) as f64 * core::f64::consts::PI * x).cos())
                        .sum::<f64>();
                let rate = dt / arm.relaxation_time_s;
                let z1 = ((1.0 - rate / 2.0) * z + rate * (ea + eb) / 2.0) / (1.0 + rate / 2.0);
                let residual = (ea + eb - z - z1) / 2.0;
                initial += 0.5 * arm.stiffness_j * (ea - z).powi(2) / 512.0;
                stored += 0.5 * arm.stiffness_j * (eb - z1).powi(2) / 512.0;
                heat += arm.stiffness_j * rate * residual.powi(2) / 512.0;
                for (force, g) in spatial_force.iter_mut().zip(gradient) {
                    *force -= arm.stiffness_j * residual * g / 512.0;
                }
            }
        }
        assert!((loss - heat).abs() < 1e-15);
        assert!((law.energy(&a, &memory).unwrap() - initial).abs() < 1e-15);
        assert!((law.energy(&b, &next).unwrap() - stored).abs() < 1e-15);
        for (a, b) in force.iter().zip(spatial_force) {
            assert!((a - b).abs() < 1e-10 * a.abs().max(1e-12));
        }
        assert!(
            (stored - initial + heat + force.iter().zip(delta).map(|(f, d)| f * d).sum::<f64>())
                .abs()
                < 1e-15
        );
        assert_eq!(memory, original);
        let mut repeat = next.clone();
        let mut again = force;
        assert_eq!(
            loss,
            law.stage(&a, &b, &memory, dt, &mut repeat, &mut again)
                .unwrap()
        );
        assert_eq!(repeat, next);
        assert_eq!(again, force);
    }

    #[test]
    fn g1_fiber_history_preserves_physical_width_partition_and_stiffness_admission() {
        let whole = law(1.0);
        let part = law(1.0 / 3.0);
        let q = [1e-8, -2e-8, 0.0, 1e-8, 2e-5, -1e-5];
        let qp: Vec<_> = q.iter().map(|q| q / 3.0_f64.sqrt()).collect();
        let mut memory = whole.relaxed_memory(&q).unwrap();
        memory[0] -= 2e-5;
        let mut fw = [0.0; 6];
        let mut fp = [0.0; 6];
        whole.reaction(&q, &memory, &mut fw).unwrap();
        part.reaction(&qp, &memory, &mut fp).unwrap();
        assert!(
            (whole.energy(&q, &memory).unwrap() - 3.0 * part.energy(&qp, &memory).unwrap()).abs()
                < 1e-15
        );
        for (a, b) in fw.iter().zip(fp) {
            assert!((a - 3.0_f64.sqrt() * b).abs() < 1e-10 * a.abs().max(1e-12));
        }
        let bound = whole.stiffness_upper_bound_per_s2(&memory).unwrap();
        let direction = [0.2, -0.3, 0.15, -0.25, 0.1, 0.4];
        let h = 1e-10;
        let plus: Vec<_> = q.iter().zip(direction).map(|(q, d)| q + h * d).collect();
        let minus: Vec<_> = q.iter().zip(direction).map(|(q, d)| q - h * d).collect();
        whole.reaction(&plus, &memory, &mut fw).unwrap();
        whole.reaction(&minus, &memory, &mut fp).unwrap();
        let tangent: f64 = fw
            .iter()
            .zip(fp)
            .zip(direction)
            .map(|((p, m), d)| -d * (p - m) / (2.0 * h))
            .sum();
        assert!(tangent.abs() <= bound * direction.iter().map(|d| d * d).sum::<f64>());
    }
}
