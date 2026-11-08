//! Mutually coherent incident beams in the same absorbing crystal.
//!
//! All components have the exact same original interface and frequency, with
//! their authored complex peak amplitudes anchored at the same global origin.
//! Directions and polarizations may differ. Every incident, reflected and
//! transmitted mode interferes before Poynting flux or dielectric work is
//! formed. This preserves spatial absorption fringes, including their dark
//! nodes, without assigning a separate positive heat source to each beam.

use super::super::OrientedCrystalFieldSide;
use super::{
    AbsorbingCrystalField, CrystalIlluminationError, OrientedUniaxialInterface,
    PreparedAbsorbingCrystalField, ReflectionError, kernel,
};

mod volume;
pub use volume::{
    AbsorbingCrystalBox, AbsorbingCrystalBoxLimits, AbsorbingCrystalCoherentBoxAbsorption,
};

type Result<T> = core::result::Result<T, CrystalIlluminationError>;

/// A borrowed monochromatic coherent sum of the original prepared fields.
/// This explicitly declares mutual coherence; use
/// [`super::AbsorbingCrystalIllumination`] for independent statistical factors.
/// Construction and each observation take O(components) work and allocate
/// nothing. Material laws and Maxwell boundary solves are never repeated.
pub struct AbsorbingCrystalCoherentIllumination<'p, 's, 'm> {
    source: &'s OrientedUniaxialInterface<'m>,
    components: &'p [PreparedAbsorbingCrystalField<'s, 'm>],
}

impl<'p, 's, 'm> AbsorbingCrystalCoherentIllumination<'p, 's, 'm> {
    /// Bind a coherent angular spectrum to its exact original interface.
    /// Each component already passed source-spectrum and incident-light-cone
    /// admission. Frequencies must agree exactly: different frequencies need
    /// an explicit temporal coherence/exposure model, not a phasor sum.
    /// Empty components mean darkness and carry no implicit field or carrier.
    pub fn try_new(
        source: &'s OrientedUniaxialInterface<'m>,
        components: &'p [PreparedAbsorbingCrystalField<'s, 'm>],
        maximum_components: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        kernel::poll_illumination(&mut cancelled)?;
        if components.len() > maximum_components {
            return Err(ReflectionError::Invalid(
                "coherent crystal illumination exceeds its component budget",
            )
            .into());
        }
        let omega = components.first().map(|field| field.omega_rad_s());
        for component in components {
            kernel::poll_illumination(&mut cancelled)?;
            if !core::ptr::eq(source, component.source()) {
                return Err(ReflectionError::Invalid(
                    "coherent crystal components require the exact original interface",
                )
                .into());
            }
            if Some(component.omega_rad_s()) != omega {
                return Err(ReflectionError::Invalid(
                    "coherent crystal components require exactly the same frequency",
                )
                .into());
            }
        }
        kernel::poll_illumination(&mut cancelled)?;
        Ok(Self { source, components })
    }

    pub const fn source(&self) -> &'s OrientedUniaxialInterface<'m> {
        self.source
    }

    pub const fn components(&self) -> &'p [PreparedAbsorbingCrystalField<'s, 'm>] {
        self.components
    }

    pub const fn component_count(&self) -> usize {
        self.components.len()
    }

    /// The common original positive carrier; an empty dark bank has none.
    pub fn omega_rad_s(&self) -> Option<f64> {
        self.components.first().map(|field| field.omega_rad_s())
    }

    /// Full coherent Cartesian peak E/H, cycle-mean flux, and local heat.
    /// Coordinates are physical global metres, with the interface at z=0;
    /// the same half-space admission applies even to an empty dark bank.
    /// Transmitted heat is omega/2 E*[-Im(epsilon)]E evaluated after ALL
    /// components and their original eigenmodes have interfered. Ordinary or
    /// extraordinary selection retains only that eigenmode from each beam.
    ///
    /// Crossed beams produce transverse variation and lateral energy flow.
    /// A receiving thermal owner must integrate this local work over its
    /// actual volume; multiplying one point by an area or adding independent
    /// single-beam slab integrals would discard these interference terms.
    pub fn field(
        &self,
        side: OrientedCrystalFieldSide,
        position_m: [f64; 3],
        cancelled: impl FnMut() -> bool,
    ) -> Result<AbsorbingCrystalField> {
        kernel::coherent_illumination_field(
            self.components.iter().map(|field| field.boundary_response()),
            side,
            position_m,
            cancelled,
        )
        .map_err(CrystalIlluminationError::from)
    }

    /// Coherent signed energy flux without requiring a representable local
    /// dielectric-work magnitude. No per-beam power fractions are added.
    pub fn poynting_w_m2(
        &self,
        side: OrientedCrystalFieldSide,
        position_m: [f64; 3],
        cancelled: impl FnMut() -> bool,
    ) -> Result<[f64; 3]> {
        kernel::coherent_illumination_poynting(
            self.components.iter().map(|field| field.boundary_response()),
            side,
            position_m,
            cancelled,
        )
        .map_err(CrystalIlluminationError::from)
    }
}
