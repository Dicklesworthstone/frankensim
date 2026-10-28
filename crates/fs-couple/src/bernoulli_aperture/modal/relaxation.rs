//! Proportional material spectra in the original plate modes, with real memory.
//!
//! If every region has C_j=(E_j/E_inf) C_inf with the same ratios and times,
//! and no prestress contributes to K, the retained modal branch stiffness is
//! exactly k_ij=(E_j/E_inf) k_i. This is a constitutive restriction, not weak
//! damping or a fitted modal decay. Nonproportional regions need full matrix
//! projection in `with_regional_material_relaxation`. All source claims remain
//! on the original chart; the original proportional entry point is unchanged.

use super::super::paired::ExactArm;
use super::{AcousticRealizeError, BladePropagator, BladeState, ModalAperture, invalid};
use crate::thin_plate::PlateMaterialModel;
use fs_exec::CancelGate;
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_material::visco::GeneralizedMaxwell;
use fs_math::det;
use fs_qty::{Density, Dims, Pressure, QuantitySpec};

mod regional;
mod tensor;
pub use tensor::{TensorModalMaterialPoint, TensorModalRelaxationLimits};

/// Initial internal strain of each physical Maxwell arm.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum InitialModalMemory {
    /// Each arm starts equilibrated at the supplied modal displacement.
    Relaxed,
    /// Zero prior viscous deformation; the initial displacement loads each arm.
    Unrelaxed,
}

/// Source spectra by plate/region and the accepted energy-normalized memories.
pub struct ModalRelaxation {
    spectra: Vec<Vec<GeneralizedMaxwell>>,
    initial: InitialModalMemory,
    pub(super) arms: Vec<Vec<ExactArm>>,
    pub(super) states: Vec<Vec<f64>>,
    pub(super) coupled_arms: Option<Vec<super::exact::CoupledArm>>,
    tensor_points: Vec<TensorModalMaterialPoint>,
}
impl ModalRelaxation {
    /// Original reference Young moduli and branches, in plate/region order.
    /// Orthotropic branches scale the whole original tensor proportionally.
    pub fn material_spectra(&self) -> &[Vec<GeneralizedMaxwell>] {
        &self.spectra
    }
    /// Original initial-strain declaration; subsequent accepted states evolve.
    pub const fn initial_memory(&self) -> InitialModalMemory {
        self.initial
    }
    /// Accepted energy-normalized viscous coordinates. Proportional selection
    /// uses one row per original mode; regional selection uses one row with
    /// every coupled material channel. No positive channel is truncated.
    pub fn memory_sqrt_j(&self) -> &[Vec<f64>] {
        &self.states
    }
    /// Whether separate original regions retain full cross-mode memory.
    pub fn is_regional(&self) -> bool { self.coupled_arms.is_some() }
    /// Full 3D material points, present only for exact tensor selection.
    /// Their source laws, normal constraints and histories remain observable.
    pub fn tensor_material_points(&self) -> &[TensorModalMaterialPoint] { &self.tensor_points }
    pub(super) fn validate_material_state(&self, modes: &[BladeState], states: &[Vec<f64>])
        -> Result<(), AcousticRealizeError>
    {
        if self.tensor_points.is_empty() { return Ok(()); }
        let memory = states.first().ok_or_else(|| invalid("tensor material memory is absent"))?;
        for point in &self.tensor_points { point.surface_strain(modes,memory)?; }
        Ok(())
    }
    /// A regional channel's original physical-coordinate row c, in sqrt(N/m).
    /// Its positive energy is `(c^T q - z)^2/2`; exact-zero rows are absent.
    pub fn regional_coupling_row(&self, channel: usize) -> Option<&[f64]> {
        self.coupled_arms.as_ref()?.get(channel).map(|arm| arm.coupling.as_slice())
    }
    /// Relaxation time [s] of that original regional material channel.
    pub fn regional_relaxation_time_s(&self, channel: usize) -> Option<f64> {
        self.coupled_arms.as_ref()?.get(channel).map(|arm| arm.relaxation_time_s)
    }
    pub(super) fn energy(&self, modes: &[BladeState], states: &[Vec<f64>]) -> f64 {
        if let Some(arms) = &self.coupled_arms {
            let memory = states.first().map_or(&[][..], Vec::as_slice);
            if memory.len() != arms.len() { return f64::NAN; }
            return arms.iter().zip(memory).map(|(arm, z)| {
                let strain = arm.coupling.iter().zip(modes)
                    .map(|(c, mode)| c * mode.displacement_m).sum::<f64>() - z;
                0.5 * strain * strain
            }).sum();
        }
        self.arms
            .iter()
            .zip(modes)
            .zip(states)
            .map(|((arms, mode), states)| {
                arms.iter()
                    .zip(states)
                    .map(|(arm, z)| {
                        let strain = arm.coupling * mode.displacement_m - z;
                        0.5 * strain * strain
                    })
                    .sum::<f64>()
            })
            .sum()
    }
}

fn property(
    material: &ResolvedMaterialStatePoint,
    key: &str,
    dims: Dims,
    band: [f64; 2],
) -> Result<f64, AcousticRealizeError> {
    let p = material
        .property(key)
        .ok_or_else(|| invalid("modal spectrum lacks a selected material property"))?;
    if p.requirement().quantity() != QuantitySpec::dimensional(dims)
        || !p.is_constant_at_fixed_temperature()
        || !p.value_si().is_finite()
    {
        return Err(invalid(
            "modal spectrum requires dimensioned constant coefficients at the original material state",
        ));
    }
    if let Some((lo, hi)) = p.angular_frequency_band_rad_s().map_err(super::nonlinear)?
        && (band[0] < lo || band[1] > hi)
    {
        return Err(invalid(
            "modal spectrum source does not cover the retained equilibrium and instantaneous frequencies",
        ));
    }
    Ok(p.value_si())
}

impl ModalAperture {
    /// Attach source-bound proportional Maxwell branches before any motion.
    /// Each pair names (additional Young modulus [Pa], relaxation time [s])
    /// properties already resolved on EVERY original plate region. The card's
    /// Young modulus is the equilibrium modulus; for orthotropic plates its
    /// axis-1 modulus defines the ratio that multiplies the COMPLETE tensor.
    /// All regions of each blade must share exactly the same normalized law.
    /// Geometry, thickness, orientation and equilibrium elasticity may differ.
    ///
    /// Uses the same ExactBlade exponential and independently integrated loss
    /// as the scalar aperture. All requested modes, contact, pressure forces and
    /// swept volume remain in the same coupled solve. Initial memory is explicit.
    /// No frequency-independent damping approximation or output filter is added.
    ///
    /// # Errors
    /// Requires ExactZoh, original material plates, zero pretension and no other
    /// material damping. At most eight source terms per region, 256 live arms
    /// total, dt/tau <= 2 and the original instantaneous-frequency clock cap.
    /// Nonproportional regions, unsupported source domains, nonfinite states,
    /// cancellation and late/repeated attachment refuse without advancing time.
    pub fn with_material_relaxation(
        mut self,
        properties: &[(String, String)],
        initial: InitialModalMemory,
        gate: &CancelGate,
    ) -> Result<Self, AcousticRealizeError> {
        let poll = || {
            if gate.is_requested() {
                Err(AcousticRealizeError::Cancelled)
            } else {
                Ok(())
            }
        };
        poll()?;
        if self.samples != 0
            || self.relaxation.is_some()
            || self.viscosity.is_some()
            || self.blade_propagator() != BladePropagator::ExactZoh
            || self.material_plates.is_empty()
            || properties.is_empty()
            || properties.len() > 8
            || self.spec.modes.iter().any(|m| m.damping_n_s_m != 0.0)
        {
            return Err(invalid(
                "material modal relaxation requires one cold exact-ZOH attachment, material plates, 1..=8 terms and no separate damping",
            ));
        }
        let mut spectra = Vec::new();
        let mut arms = vec![Vec::<ExactArm>::new(); self.state.len()];
        for plate in &self.material_plates {
            poll()?;
            if plate.options().reduction.assembly.pretension != 0.0
                || plate.source().bending_viscosity().is_some()
            {
                return Err(invalid(
                    "proportional plate relaxation cannot relax geometric prestress or duplicate material viscosity",
                ));
            }
            let range = plate.mode_range();
            let low = self.spec.modes[range.clone()]
                .iter()
                .map(|m| det::sqrt(m.stiffness_n_m / m.mass_kg))
                .fold(f64::INFINITY, f64::min);
            let high = self.spec.modes[range.clone()]
                .iter()
                .map(|m| det::sqrt(m.stiffness_n_m / m.mass_kg))
                .fold(0.0, f64::max);
            let mut common = None::<Vec<(f64, f64)>>;
            let mut regional = Vec::new();
            for region in plate.source().regions() {
                poll()?;
                let material = &region.input.material;
                let key = match region.input.model {
                    PlateMaterialModel::Isotropic => "young_modulus",
                    _ => "young_modulus_1",
                };
                let e_inf = property(material, key, Pressure::DIMS, [low, high])?;
                let terms = properties
                    .iter()
                    .map(|(modulus, time)| {
                        Ok((
                            property(material, modulus, Pressure::DIMS, [low, high])?,
                            property(material, time, Dims([0, 0, 1, 0, 0, 0]), [low, high])?,
                        ))
                    })
                    .collect::<Result<Vec<_>, AcousticRealizeError>>()?;
                let law = GeneralizedMaxwell::new(e_inf, terms).map_err(super::nonlinear)?;
                let normalized: Vec<_> =
                    law.terms.iter().map(|(e, tau)| (e / e_inf, *tau)).collect();
                if common.as_ref().is_some_and(|c| c != &normalized) {
                    return Err(invalid(
                        "nonproportional regional relaxation requires full matrix projection; no modal averaging is admitted",
                    ));
                }
                let factor = 1.0 + normalized.iter().map(|(ratio, _)| ratio).sum::<f64>();
                let instantaneous = high * det::sqrt(factor);
                if !instantaneous.is_finite()
                    || instantaneous * self.spec.time_step_s > plate.options().maximum_angular_step
                    || normalized
                        .iter()
                        .any(|(ratio, tau)| *ratio > 0.0 && self.spec.time_step_s / tau > 2.0)
                {
                    return Err(invalid(
                        "material spectrum exceeds its instantaneous-frequency or relaxation-time resolution",
                    ));
                }
                for (modulus, time) in properties {
                    property(material, modulus, Pressure::DIMS, [low, instantaneous])?;
                    property(
                        material,
                        time,
                        Dims([0, 0, 1, 0, 0, 0]),
                        [low, instantaneous],
                    )?;
                }
                // The instantaneous branch stiffness can reach beyond the
                // equilibrium eigenwindow. Keep every coefficient of the
                // original mass/elastic pencil inside its own source domain.
                let elastic: &[(&str, Dims)] = match region.input.model {
                    PlateMaterialModel::Isotropic => &[
                        ("density", Density::DIMS),
                        ("young_modulus", Pressure::DIMS),
                        ("poisson_ratio", Dims::NONE),
                    ],
                    _ => &[
                        ("density", Density::DIMS),
                        ("young_modulus_1", Pressure::DIMS),
                        ("young_modulus_2", Pressure::DIMS),
                        ("shear_modulus_12", Pressure::DIMS),
                        ("poisson_ratio_12", Dims::NONE),
                    ],
                };
                for &(key, dims) in elastic {
                    property(material, key, dims, [low, instantaneous])?;
                }
                common = Some(normalized);
                regional.push(law);
            }
            let common = common.ok_or_else(|| invalid("material plate has no original regions"))?;
            for i in range {
                for &(ratio, relaxation_time_s) in &common {
                    if ratio == 0.0 {
                        continue;
                    }
                    let coupling = det::sqrt(ratio * self.spec.modes[i].stiffness_n_m);
                    if !coupling.is_finite() || coupling <= 0.0 {
                        return Err(invalid(
                            "material modal projection lost a positive Maxwell arm",
                        ));
                    }
                    arms[i].push(ExactArm {
                        coupling,
                        relaxation_time_s,
                    });
                }
            }
            spectra.push(regional);
        }
        let states: Vec<Vec<f64>> = arms
            .iter()
            .zip(&self.state)
            .map(|(arms, state)| {
                arms.iter()
                    .map(|arm| match initial {
                        InitialModalMemory::Relaxed => arm.coupling * state.displacement_m,
                        InitialModalMemory::Unrelaxed => 0.0,
                    })
                    .collect()
            })
            .collect();
        if states.iter().flatten().any(|z| !z.is_finite()) {
            return Err(invalid(
                "initial modal material memory is not representable",
            ));
        }
        let exact = super::exact::ExactModes::with_arms(&self.spec, arms.clone())?;
        poll()?;
        let relaxation = ModalRelaxation {
            spectra,
            initial,
            arms,
            states,
            coupled_arms: None,
            tensor_points: Vec::new(),
        };
        let energy = self.energy(&self.state) + relaxation.energy(&self.state, &relaxation.states);
        if !energy.is_finite() || energy > self.spec.limits.maximum_energy_j {
            return Err(invalid(
                "initial material modal memory exceeds the original energy budget",
            ));
        }
        self.exact = Some(exact);
        self.relaxation = Some(relaxation);
        Ok(self)
    }

    /// Original source spectra and accepted material histories, never a fitted decay.
    pub fn material_relaxation(&self) -> Option<&ModalRelaxation> {
        self.relaxation.as_ref()
    }
}
