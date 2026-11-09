//! Source-bound rigid-frame porous propagation in the original DG owner.
//!
//! The high-frequency effective inertia/compliance belong to the native DG
//! mass matrices. Positive JCA/JCAL viscous and thermal memories are appended
//! by the SAME uniform diffusion stamp used by circular wall loss. Thus one
//! continuous operator and one native exponential advance pressure, superficial
//! flow and every memory coordinate; no impedance is frozen at a drive tone.
//! All declared frequency-domain queries retain the producer's source band.

use super::*;
use crate::fluid::FluidDuctError;
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_material::fluid_acoustics::thermoviscous::ThermoviscousFluid;
use fs_material::porous::diffusive::PorousDiffusive;

mod heat;
pub use heat::{DgPorousElementHeat, DgPorousHeat, DgPorousHeating};

fn invalid(what: &'static str) -> DuctError {
    DuctError::BadParameter { what }
}

/// Original saturating-fluid source plus its admitted rigid porous realization.
/// The frame's source/calibration remains on the realization; the actual fluid
/// owner retains its EOS, phase choices, properties and source receipts.
#[derive(Clone, Debug)]
pub struct DgPorousMaterial<M> {
    source: M,
    model: PorousDiffusive,
    acoustic: FluidAcousticState,
}

impl<M> DgPorousMaterial<M> {
    pub const fn source(&self) -> &M {
        &self.source
    }
    pub const fn model(&self) -> &PorousDiffusive {
        &self.model
    }
    /// High-frequency effective storage, already including porosity once.
    pub const fn acoustic(&self) -> FluidAcousticState {
        self.acoustic
    }
}
impl<M: Clone> fluid::sealed::Medium for DgPorousMaterial<M> {}
impl<M: Clone> DgMedium for DgPorousMaterial<M> {
    fn density_kg_m3(&self) -> f64 {
        self.acoustic.density_kg_m3()
    }
    fn sound_speed_m_s(&self) -> f64 {
        self.acoustic.sound_speed_m_s()
    }
    fn admit(&self) -> Result<(), DuctError> {
        Ok(())
    }
    fn porous_diffusive_model(&self) -> Option<&PorousDiffusive> {
        Some(&self.model)
    }
}

impl DgDuct {
    /// Bind a prepared JCA/JCAL realization and its ORIGINAL saturating fluid
    /// to the existing homogeneous, constant-area DG propagation owner.
    ///
    /// The full cross-sectional area carries superficial volume flow, exactly
    /// as in `PorousLayer`. No circular-pore wall loss is added to this porous
    /// law. High-frequency effective rho and K provide the acoustic masses;
    /// drag and both sets of positive arms add their own physical storage/loss.
    /// Generic liquids keep their actual Ks and Cp/Cv, without a gas EOS.
    ///
    /// `model` must have been admitted for this same fluid snapshot. This is a
    /// frozen small-signal rigid-frame specimen, not Biot poroelasticity, mean
    /// flow, acoustic startup qualification outside the declared band, or live
    /// temperature feedback. The continuous harmonic/scattering/pole and time
    /// owners all consume the same completed operator.
    ///
    /// # Errors
    /// Changed fluid tuple, nonuniform section, native geometric/DOF limits,
    /// or an unrepresentable augmented operator. No existing state is mutated.
    pub fn from_porous_material<M: ThermoviscousFluid + Clone>(
        duct: &Duct,
        model: &PorousDiffusive,
        fluid: &M,
        mesh: DgMesh,
        inlet: DgBoundary,
        outlet: DgBoundary,
    ) -> Result<DgDuct<DgPorousMaterial<M>>, FluidDuctError> {
        let actual = fluid.thermoviscous_state()?;
        if actual != model.fluid() {
            return Err(invalid(
                "porous DG requires the original fluid tuple admitted by its model",
            )
            .into());
        }
        if viscothermal::uniform_radius(duct).is_none() {
            return Err(invalid(
                "porous DG currently requires one homogeneous constant-area specimen",
            )
            .into());
        }
        let acoustic = FluidAcousticState::from_adiabatic_bulk_modulus(
            model.density_infinite_kg_m3(),
            model.bulk_infinite_pa(),
        )?;
        let medium = DgPorousMaterial {
            source: fluid.clone(),
            model: model.clone(),
            acoustic,
        };
        let mut result = DgDuct::assemble(
            duct,
            &medium,
            &medium,
            mesh,
            inlet,
            outlet,
            None,
            &mut |_| Ok::<_, DuctError>(acoustic),
        )?;
        result.stamp_uniform_diffusion(
            model.series_drag_per_s(),
            model
                .series_terms()
                .iter()
                .map(|t| (t.gain_per_s, t.rate_per_s)),
            model
                .shunt_terms()
                .iter()
                .map(|t| (t.gain_per_s, t.rate_per_s)),
        )?;
        Ok(result)
    }
}

impl<M: Clone> DgDuct<DgPorousMaterial<M>> {
    /// Original saturating-fluid source, porous realization and admitted storage.
    pub const fn porous_material(&self) -> &DgPorousMaterial<M> {
        &self.gas
    }
}

/// Separate instantaneous intrinsic material source powers. These exclude
/// terminals, numerical damping and any attached structural dissipation.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct DgPorousPower {
    pub viscous_w: f64,
    pub thermal_w: f64,
}

impl<M: DgMedium> DgDuct<M> {
    /// Original porous law, fluid snapshot, positive arms and sampled admission.
    /// None means this duct has no intrinsic porous material.
    pub fn porous_diffusive_model(&self) -> Option<&PorousDiffusive> {
        self.gas.porous_diffusive_model()
    }

    /// Observe the original static Darcy drag and every material memory's
    /// positive source square. This uses neither a total energy deficit nor
    /// the rounded symmetric part of the augmented generator.
    pub fn porous_dissipated_power_w(&self, state: &DgState) -> Result<DgPorousPower, DuctError> {
        if state.operator_identity != self.operator_identity
            || state.tilde.len() != self.total_dofs()
            || state.tilde.iter().any(|value| !value.is_finite())
        {
            return Err(invalid(
                "porous source power needs a finite state of its complete DG operator",
            ));
        }
        let (viscous_w, thermal_w) = self.porous_loss_pair::<false>(&state.tilde, &state.tilde);
        if !viscous_w.is_finite() || !thermal_w.is_finite() {
            return Err(invalid("porous source power is unrepresentable"));
        }
        Ok(DgPorousPower {
            viscous_w,
            thermal_w,
        })
    }

    pub(super) fn porous_loss_pair<const CROSS: bool>(
        &self,
        state: &[f64],
        tangent: &[f64],
    ) -> (f64, f64) {
        let Some(model) = self.gas.porous_diffusive_model() else {
            return (0.0, 0.0);
        };
        let m = self.field_dofs;
        let mut viscous = model.series_drag_per_s()
            * (0..m)
                .map(|i| state[m + i] * if CROSS { tangent[m + i] } else { state[m + i] })
                .sum::<f64>();
        let mut thermal = 0.0;
        for (branch, (offset, term, is_thermal)) in model
            .series_terms()
            .iter()
            .map(|t| (m, t, false))
            .chain(model.shunt_terms().iter().map(|t| (0, t, true)))
            .enumerate()
        {
            let root_gain = det::sqrt(term.gain_per_s);
            let root_rate = det::sqrt(term.rate_per_s);
            for i in 0..m {
                let memory = 2 * m + branch * m + i;
                let slip = root_gain * state[offset + i] - root_rate * state[memory];
                let other = if CROSS {
                    root_gain * tangent[offset + i] - root_rate * tangent[memory]
                } else {
                    slip
                };
                if is_thermal {
                    thermal += slip * other;
                } else {
                    viscous += slip * other;
                }
            }
        }
        (viscous, thermal)
    }
}

#[cfg(test)]
mod tests;
