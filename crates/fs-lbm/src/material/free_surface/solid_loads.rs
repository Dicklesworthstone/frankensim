//! Source-scaled immersed-body reactions on the original free-surface stream.

use super::*;
pub use crate::d3q19::freesurface3::SubmergedSolidLoadControls3;
use crate::d3q19::freesurface3::{FreeSurfaceError3, SubmergedSolid3, SubmergedSolidLoads3};

/// One original stationary solid tile assembly with a physical moment origin.
#[derive(Clone, Copy, Debug)]
pub struct MaterialSubmergedSolid3<'a> {
    pub solid_tiles: &'a [(u32, u32, u32)],
    pub moment_origin_m: [f64; 3],
}

/// Actual hydrodynamic reaction ON one completely submerged stationary object.
/// Fluid impulse and angular impulse are the negatives of these values.
/// Uniform pressure cancels on the complete closed surface; no absolute
/// pressure datum is needed. Stationary boundary work is zero.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialSubmergedSolidLoad3 {
    pub impulse_kg_m_s: [f64; 3],
    pub angular_impulse_kg_m2_s: [f64; 3],
    /// Mean force over the original physical timestep.
    pub force_n: [f64; 3],
    /// Mean torque about this object's declared physical origin.
    pub torque_n_m: [f64; 3],
    pub measured_links: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub struct MaterialSubmergedSolidLoads3 {
    pub bodies: Vec<MaterialSubmergedSolidLoad3>,
    pub geometry_epoch: u64,
    pub mass_before_kg: f64,
    pub mass_after_kg: f64,
    pub relative_mass_defect: f64,
    pub maximum_speed_m_s: f64,
}

fn vector(value: [f64; 3], unit: f64) -> Result<[f64; 3], MaterialLbmError> {
    let mut result = [0.0; 3];
    for a in 0..3 {
        result[a] = finite(
            value[a] * unit,
            "submerged-solid SI reaction is unrepresentable",
        )?;
    }
    Ok(result)
}

impl MaterialFreeSurface3<'_> {
    /// Advance the original fluid once and return SI force/torque on all
    /// selected stationary solid objects, in input-object order.
    ///
    /// The original source density, viscosity, dx/dt, gravity, capillarity,
    /// Guo force history and VOF evolution remain authoritative. Loads are
    /// actual reflected-population impulses, not a hydrostatic approximation.
    /// Their complete closed surfaces must stay fully submerged. Constant
    /// reference pressure has zero resultant and moment on each such surface.
    /// No moving-body dynamics, partial-wetting area or thermal work is inferred.
    ///
    /// # Errors
    /// The original tile/link/population/mass/native/cancellation gates apply.
    /// Physical endpoint speeds must also obey the original material scenario,
    /// and every SI quantity must be representable before publication.
    pub fn step_with_submerged_solid_loads_si(
        &mut self,
        bodies: &[MaterialSubmergedSolid3<'_>],
        controls: SubmergedSolidLoadControls3,
        gate: &CancelGate,
    ) -> Result<MaterialSubmergedSolidLoads3, MaterialLbmError> {
        let (candidate, report) = self.submerged_solid_load_candidate_si(bodies, controls, gate)?;
        if gate.is_requested() {
            return Err(FreeSurfaceError3::Cancelled.into());
        }
        *self = candidate;
        Ok(report)
    }

    /// Propose the same fluid/load result without publishing this liquid.
    /// Coupled mechanics can finish its own acceptance before installing the
    /// returned candidate. Dropping it preserves the original fluid state.
    pub fn submerged_solid_load_candidate_si(
        &self,
        bodies: &[MaterialSubmergedSolid3<'_>],
        controls: SubmergedSolidLoadControls3,
        gate: &CancelGate,
    ) -> Result<(Self, MaterialSubmergedSolidLoads3), MaterialLbmError> {
        if gate.is_requested() {
            return Err(FreeSurfaceError3::Cancelled.into());
        }
        if bodies.len() > controls.maximum_solid_tiles {
            return Err(MaterialLbmError::Invalid(
                "submerged-solid object budget exceeded",
            ));
        }
        let s = &self.scaling;
        let impulse_unit = positive(
            s.mass_unit_kg * s.velocity_scale,
            "submerged-solid impulse unit is unrepresentable",
        )?;
        let angular_unit = positive(
            impulse_unit * s.cell_size_m,
            "submerged-solid angular impulse unit is unrepresentable",
        )?;
        let force_unit = positive(
            impulse_unit / s.timestep_s,
            "submerged-solid force unit is unrepresentable",
        )?;
        let torque_unit = positive(
            force_unit * s.cell_size_m,
            "submerged-solid torque unit is unrepresentable",
        )?;
        let mut selection = Vec::with_capacity(bodies.len());
        for body in bodies {
            let mut origin = [0.0; 3];
            for a in 0..3 {
                origin[a] = finite(
                    body.moment_origin_m[a] / s.cell_size_m,
                    "submerged-solid moment origin is unrepresentable",
                )?;
            }
            selection.push(SubmergedSolid3 {
                solid_tiles: body.solid_tiles,
                moment_origin: origin,
            });
        }
        let (candidate, raw): (_, SubmergedSolidLoads3) = self
            .lattice
            .submerged_solid_load_candidate(&selection, controls, gate)?;
        if raw.maximum_speed > s.flow.u_lattice {
            return Err(MaterialLbmError::Invalid(
                "submerged-solid liquid exceeds its original material scenario speed",
            ));
        }
        let mut reactions = Vec::with_capacity(raw.bodies.len());
        for body in raw.bodies {
            reactions.push(MaterialSubmergedSolidLoad3 {
                impulse_kg_m_s: vector(body.wall_impulse, impulse_unit)?,
                angular_impulse_kg_m2_s: vector(body.wall_angular_impulse, angular_unit)?,
                force_n: vector(body.wall_impulse, force_unit)?,
                torque_n_m: vector(body.wall_angular_impulse, torque_unit)?,
                measured_links: body.measured_links,
            });
        }
        let report = MaterialSubmergedSolidLoads3 {
            bodies: reactions,
            geometry_epoch: raw.geometry_epoch,
            mass_before_kg: positive(
                raw.mass_before * s.mass_unit_kg,
                "submerged-solid initial physical mass is unrepresentable",
            )?,
            mass_after_kg: positive(
                raw.mass_after * s.mass_unit_kg,
                "submerged-solid final physical mass is unrepresentable",
            )?,
            relative_mass_defect: raw.relative_mass_defect,
            maximum_speed_m_s: finite(
                raw.maximum_speed * s.velocity_scale,
                "submerged-solid physical speed is unrepresentable",
            )?,
        };
        if gate.is_requested() {
            return Err(FreeSurfaceError3::Cancelled.into());
        }
        Ok((
            Self {
                lattice: candidate,
                scaling: self.scaling,
            },
            report,
        ))
    }
}
