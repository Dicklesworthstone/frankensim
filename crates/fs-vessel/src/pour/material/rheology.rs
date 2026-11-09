//! Source-driven shear viscosity on the existing SI weir pour.
//!
//! The actual force-corrected D2Q9 strain re-queries the original material
//! claims before every collision. All geometry, rotating-frame motion, native
//! mass transport, contact choice, SI reporting and endpoint limits are shared
//! with MaterialPour. Temperature, density and capillarity remain fixed.

use fs_exec::CancelGate;
use fs_lbm::freesurface::ContactModel;
use fs_lbm::material::free_surface::MaterialFreeSurfaceScaling;
use fs_lbm::material::free_surface::rheology::{
    IsothermalRheologyLaw, MaterialFreeSurfaceRheology2, RheologySolve,
};

use super::{MaterialPour, MaterialPourError, MaterialPourOutcome, MaterialPourSpec};

/// An admitted physical pour with one immutable, evolving-shear source law.
///
/// The source may thin or thicken with physical rate; no analytic fit or frozen
/// apparent viscosity replaces it. Source and root refusals return no partial
/// pour. Absolute-pressure dependence, temperature evolution, viscoelastic
/// memory and yield history require other owners and are not inferred here.
pub struct MaterialRheologyPour<'s, 'm> {
    physical: MaterialPour<'s>,
    rheology: MaterialFreeSurfaceRheology2<'s, 'm>,
}

/// Native field and physical report, plus the last collision's source receipts.
pub struct MaterialRheologyPourOutcome<'s, 'm> {
    pour: MaterialPourOutcome<'s>,
    rheology: MaterialFreeSurfaceRheology2<'s, 'm>,
}

impl<'s, 'm> MaterialRheologyPour<'s, 'm> {
    pub fn try_new(
        law: &'s IsothermalRheologyLaw<'m>,
        spec: MaterialPourSpec,
        solve: RheologySolve,
    ) -> Result<Self, MaterialPourError> {
        let rheology = MaterialFreeSurfaceRheology2::try_new(
            law,
            spec.cells,
            spec.cell_size_m,
            spec.timestep_s,
            spec.maximum_speed_m_s,
            solve,
        )?;
        let physical = MaterialPour::prepare(spec, || Ok(*rheology.scaling()))?;
        Ok(Self { physical, rheology })
    }

    pub const fn scaling(&self) -> &MaterialFreeSurfaceScaling<'s> {
        self.rheology.scaling()
    }

    pub const fn spec(&self) -> &MaterialPourSpec {
        self.physical.spec()
    }

    pub const fn law(&self) -> &'s IsothermalRheologyLaw<'m> {
        self.rheology.law()
    }

    /// Advance one source-resolved collision per native pour timestep.
    ///
    /// The last accepted collision samples precede streaming/conversion:
    /// they describe the actual viscosity used at (steps-1)*dt, not a second
    /// source sweep of the final field. Rates are objective physical shear
    /// in the rotating vessel; the speed envelope is vessel-frame speed.
    /// Cancellation/source/native refusal publishes no partial outcome and
    /// leaves this immutable preparation available for an exact retry.
    pub fn run(
        &self,
        contact: ContactModel,
        gate: &CancelGate,
    ) -> Result<MaterialRheologyPourOutcome<'s, 'm>, MaterialPourError> {
        // The preparation has no sampled epochs; this clone starts with empty
        // buffers, then reuses the same two source buffers throughout the run.
        let mut rheology = self.rheology.clone();
        let pour = self.physical.run_with(contact, gate, |liquid, _step| {
            rheology.step(liquid, gate)?;
            Ok(())
        })?;
        Ok(MaterialRheologyPourOutcome { pour, rheology })
    }

    /// Balance initial liquid pressure against gravity, then evolve source viscosity.
    ///
    /// Uses [`MaterialPour::run_hydrostatic`]'s complete below-crest layers,
    /// neutral contact and zero surface tension. The returned `pour()` retains
    /// the positive free-surface pressure datum and reports the actual initial
    /// mass, including bounded numerical compression. This changes neither
    /// the source's fixed density nor its original query coordinates.
    ///
    /// Every collision still closes physical shear against the original law;
    /// no Newtonian approximation or frozen apparent viscosity is introduced.
    /// `last_collision()` retains the same source receipts as [`Self::run`].
    /// The pressure datum does not add pressure-dependent rheology, and all
    /// original source, root, speed, density and cancellation refusals apply.
    pub fn run_hydrostatic(
        &self,
        free_surface_pressure_pa: f64,
        gate: &CancelGate,
    ) -> Result<MaterialRheologyPourOutcome<'s, 'm>, MaterialPourError> {
        let mut rheology = self.rheology.clone();
        let pour = self.physical.run_hydrostatic_with(
            free_surface_pressure_pa,
            gate,
            |liquid, _step| {
                rheology.step(liquid, gate)?;
                Ok(())
            },
        )?;
        Ok(MaterialRheologyPourOutcome { pour, rheology })
    }
}

impl<'s, 'm> MaterialRheologyPourOutcome<'s, 'm> {
    pub const fn pour(&self) -> &MaterialPourOutcome<'s> {
        &self.pour
    }

    /// Complete source states/rates used by the final native collision.
    /// Cell wetness may subsequently change during streaming and conversion.
    pub const fn last_collision(&self) -> &MaterialFreeSurfaceRheology2<'s, 'm> {
        &self.rheology
    }
}
