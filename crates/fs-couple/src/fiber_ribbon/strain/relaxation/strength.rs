//! Total original-layer stress from the SAME native axial and bending histories.
//! Memory is borrowed, never restaged or reset by this admission gate.
use super::super::strength::{FiberTensileAssessment, invalid_strength, product, response, sum};
use super::*;
use crate::fiber_ribbon::relaxation::FiberRelaxationChannel;
use fs_nlmodal::longitudinal::{MAX_LONGITUDINAL_MODES, MAX_TRANSVERSE_MODES};

#[derive(Clone, Debug)]
pub(super) struct FiberLayerArm {
    pub layer: usize,
    pub delta_pa: f64,
}

#[derive(Clone, Debug)]
pub(super) struct FiberBendingHistory {
    index: usize,
    channel: FiberRelaxationChannel,
    /// Delta E / c, with z stored by the original exact owner as c*q.
    stress_per_memory: f64,
}

impl FiberStrainHistory {
    pub(in crate::fiber_ribbon) fn bind_bending_history(
        &mut self,
        first: usize,
        channels: &[FiberRelaxationChannel],
        source: &ConcentricPronySection,
    ) -> Result<(), FiberRibbonError> {
        if self.strain.tensile.is_none() {
            return Ok(());
        }
        self.bending_history = channels
            .iter()
            .enumerate()
            .map(|(i, channel)| {
                let delta = source
                    .spectra()
                    .get(channel.layer)
                    .and_then(|s| s.terms.get(channel.arm))
                    .map(|arm| arm.0)
                    .ok_or(FiberRibbonError::Invalid(
                        "fiber tensile history lost its original source arm",
                    ))?;
                let stress_per_memory = delta / channel.coupling_per_s;
                if !stress_per_memory.is_finite() || stress_per_memory <= 0.0 {
                    return Err(FiberRibbonError::Invalid(
                        "fiber tensile history stress map is unrepresentable",
                    ));
                }
                Ok(FiberBendingHistory {
                    index: first.checked_add(i).ok_or(FiberRibbonError::Invalid(
                        "fiber tensile memory offset overflows",
                    ))?,
                    channel: *channel,
                    stress_per_memory,
                })
            })
            .collect::<Result<_, _>>()?;
        Ok(())
    }

    pub(super) fn validate_tensile_binding(
        &self,
        coordinates: &[usize],
        branches: &[fs_phs::RelaxationBranch],
    ) -> Result<(), ModalCouplingError> {
        if self.strain.tensile.is_none() {
            return Ok(());
        }
        let normal = self.strain.slope_scales.len();
        if coordinates.len() != self.coordinate_count()
            || self.bending_history.len() != self.layer_arms.len() * normal
        {
            return Err(invalid_strength());
        }
        for original in &self.bending_history {
            let branch = branches.get(original.index).ok_or_else(invalid_strength)?;
            let mode = coordinates
                .get(original.channel.mode)
                .ok_or_else(invalid_strength)?;
            let displacement = mode.checked_mul(2).ok_or_else(invalid_strength)?;
            if branch.stiffness != 1.0
                || branch.relaxation_time_s != original.channel.relaxation_time_s
                || branch.projection.get(displacement) != Some(&original.channel.coupling_per_s)
                || branch
                    .projection
                    .iter()
                    .enumerate()
                    .any(|(i, &value)| i != displacement && value != 0.0)
            {
                return Err(ModalCouplingError::Invalid {
                    what: "fiber tensile admission must retain its original bending material branches",
                });
            }
        }
        Ok(())
    }

    pub(super) fn tensile_assessment(
        &self,
        q: &[f64],
        velocity: &[f64],
        memory: &[f64],
        linear_memory: &[f64],
    ) -> Result<Option<FiberTensileAssessment>, ModalCouplingError> {
        let Some(guard) = &self.strain.tensile else {
            return Ok(None);
        };
        let axial = self.strain.model.axial_modes().len();
        let normal = self.strain.slope_scales.len();
        if memory.len() != self.material.memory_count()
            || memory.iter().any(|z| !z.is_finite())
            || self.bending_history.len() != self.layer_arms.len() * normal
        {
            return Err(invalid_strength());
        }
        let mut epsilon = [0.0; MAX_LONGITUDINAL_MODES + 1];
        let mut rate = [0.0; MAX_LONGITUDINAL_MODES + 1];
        let (mean, mean_rate) = self.strain.tensile_kinematics(
            q,
            velocity,
            &mut epsilon[1..=axial],
            &mut rate[1..=axial],
        )?;
        epsilon[0] = mean;
        rate[0] = mean_rate;
        let width = axial + 1;
        let capacity = guard.capacity.elastic_capacity();
        let preload = self.strain.tension_n / self.strain.axial_rigidity_n;
        let mut force = [0.0; MAX_LONGITUDINAL_MODES + 1];
        let mut result = FiberTensileAssessment {
            per_fiber_tension_bound_n: 0.0,
            curvature_bound_per_m: 0.0,
            curvature_rate_bound_per_m_s: 0.0,
            governing_layer: 0,
            governing_tensile_stress_bound_pa: 0.0,
            utilization: 0.0,
        };
        for (mode, &scale) in guard.curvature_scales.iter().enumerate() {
            result.curvature_bound_per_m = sum(
                result.curvature_bound_per_m,
                product(scale, q[axial + mode])?.abs(),
            )?;
            result.curvature_rate_bound_per_m_s = sum(
                result.curvature_rate_bound_per_m_s,
                product(scale, velocity[axial + mode])?.abs(),
            )?;
        }
        let mut stress = [0.0; MAX_LONGITUDINAL_MODES + 1];
        let mut bending = [0.0; MAX_TRANSVERSE_MODES];
        // Both original source sequences are layer-major. Walk each memory
        // once, retaining linear cost in the admitted state count.
        let mut arm = 0;
        let mut bend = 0;
        for (layer, coefficients) in guard.capacity.layers().iter().enumerate() {
            for j in 0..width {
                stress[j] = response(
                    coefficients.young_modulus_pa(),
                    if j == 0 {
                        sum(preload, epsilon[j])?
                    } else {
                        epsilon[j]
                    },
                    coefficients.axial_viscosity_pa_s(),
                    rate[j],
                )?;
            }
            while let Some(original) = self.layer_arms.get(arm).filter(|arm| arm.layer == layer) {
                for j in 0..width {
                    stress[j] = sum(
                        stress[j],
                        product(
                            original.delta_pa,
                            sum(epsilon[j], -memory[arm * width + j])?,
                        )?,
                    )?;
                }
                arm += 1;
            }
            let geometry = &capacity.section().layer_integrals()[layer];
            let mut maximum = stress[0];
            for j in 0..width {
                force[j] = sum(force[j], product(geometry.area_m2(), stress[j])?)?;
                if j > 0 {
                    maximum = sum(maximum, stress[j].abs())?;
                }
            }
            for (mode, stress) in bending[..normal].iter_mut().enumerate() {
                *stress = response(
                    coefficients.young_modulus_pa(),
                    q[axial + mode],
                    coefficients.bending_viscosity_pa_s(),
                    velocity[axial + mode],
                )?;
            }
            while let Some(original) = self
                .bending_history
                .get(bend)
                .filter(|history| history.channel.layer == layer)
            {
                let z = *linear_memory
                    .get(original.index)
                    .ok_or_else(invalid_strength)?;
                let mode = original
                    .channel
                    .mode
                    .checked_sub(axial)
                    .ok_or_else(invalid_strength)?;
                let residual = sum(
                    product(original.channel.coupling_per_s, q[axial + mode])?,
                    -z,
                )?;
                bending[mode] = sum(
                    bending[mode],
                    product(original.stress_per_memory, residual)?,
                )?;
                bend += 1;
            }
            let mut curvature_stress = 0.0;
            for (&scale, &stress) in guard.curvature_scales.iter().zip(&bending[..normal]) {
                curvature_stress = sum(curvature_stress, product(scale, stress)?.abs())?;
            }
            let radius = capacity.section().layers()[layer].outer_radius_m;
            let amplitude = product(radius, curvature_stress)?;
            let assessment = capacity
                .assess_layer_normal_stress(layer, maximum, amplitude)
                .map_err(|_| invalid_strength())?;
            if assessment.tensile_utilization > result.utilization {
                result.governing_layer = layer;
                result.governing_tensile_stress_bound_pa =
                    assessment.maximum_axial_stress_pa.max(0.0);
                result.utilization = assessment.tensile_utilization;
            }
        }
        result.per_fiber_tension_bound_n = force[0];
        for value in &force[1..width] {
            result.per_fiber_tension_bound_n = sum(result.per_fiber_tension_bound_n, value.abs())?;
        }
        Ok(Some(result))
    }
}
