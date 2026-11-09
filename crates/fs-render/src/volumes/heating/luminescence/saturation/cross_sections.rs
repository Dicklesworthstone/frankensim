//! One molecular absorption source for external and recycled photons.
//!
//! Rows follow the caller's actual vacuum-wavelength nodes; columns are
//! transitions `[ground -> singlet, ground -> triplet]`. Both transitions
//! compete for the same photons and original finite ground inventory.
//! Cross sections are square metres per emitter, distinct from the host's
//! passive mass absorption. The source carries no populations or optical path.

#[cfg(feature = "material-state")]
use super::finite;
use super::{ContentHash, OpticalHeatingError, invalid, nonnegative};
#[cfg(feature = "material-state")]
use std::sync::Arc;

#[derive(Clone, Debug)]
pub(super) enum CrossSections<const B: usize> {
    Authored {
        identity: ContentHash,
        values: [[f64; 2]; B],
    },
    #[cfg(feature = "material-state")]
    Material {
        card: Arc<fs_matdb::MaterialCard>,
        source: fs_material::phosphorescence::pump::GroundStatePumpStatePoint,
    },
}

/// Molecular absorption at the actual photon wavelengths of one original cell.
///
/// The same carrier serves external pump beams and reabsorbed emission. Source
/// claims describe a molecule at a physical state, regardless of photon origin.
/// An authored row belongs to its caller-declared spectral node; material rows
/// are resolved at every requested node and accepted cell temperature.
#[derive(Clone, Debug)]
pub struct EmitterAbsorptionCrossSections<const B: usize> {
    pub(super) grid_cell: usize,
    pub(super) source: CrossSections<B>,
}

impl<const B: usize> EmitterAbsorptionCrossSections<B> {
    /// Bind explicit, nonnegative square metres per ground-state emitter.
    /// The source identity attributes authored data; it does not certify it.
    pub fn authored(
        grid_cell: usize,
        source_identity: ContentHash,
        cross_sections_m2: [[f64; 2]; B],
    ) -> Result<Self, OpticalHeatingError> {
        if B == 0
            || source_identity == ContentHash([0; 32])
            || !cross_sections_m2.into_iter().flatten().all(nonnegative)
        {
            return Err(invalid("invalid molecular absorption cross sections"));
        }
        Ok(Self {
            grid_cell,
            source: CrossSections::Authored {
                identity: source_identity,
                values: cross_sections_m2,
            },
        })
    }

    /// Bind the original selected molecular absorption claims. Both vacuum
    /// wavelength and absolute temperature must be explicit owned coordinates;
    /// every other condition remains fixed by the native query. Unsupported
    /// wavelengths and temperatures refuse without substituting a coefficient.
    #[cfg(feature = "material-state")]
    pub fn from_material(
        grid_cell: usize,
        card: Arc<fs_matdb::MaterialCard>,
        source: fs_material::phosphorescence::pump::GroundStatePumpStatePoint,
    ) -> Result<Self, OpticalHeatingError> {
        if B == 0 || source.temperature_axis().is_none() {
            return Err(invalid(
                "thermal molecular absorption requires spectral nodes and an owned absolute-temperature axis",
            ));
        }
        source
            .at_query_point(&card, &source.resolved().query())
            .map_err(|_| invalid("absorption cross sections differ from original material card"))?;
        Ok(Self {
            grid_cell,
            source: CrossSections::Material { card, source },
        })
    }

    pub const fn grid_cell(&self) -> usize {
        self.grid_cell
    }

    pub(super) fn identity(&self) -> ContentHash {
        match &self.source {
            CrossSections::Authored { identity, .. } => *identity,
            #[cfg(feature = "material-state")]
            CrossSections::Material { source, .. } => source.identity(),
        }
    }

    /// Resolve actual photon nodes, sharing the same physical source for a
    /// pump or recycled emission. No spectral row is borrowed from a neighbor.
    pub(super) fn evaluate(
        &self,
        wavelengths_nm: [f64; B],
        temperature_k: f64,
    ) -> Result<[[f64; 2]; B], OpticalHeatingError> {
        if !temperature_k.is_finite()
            || temperature_k <= 0.0
            || wavelengths_nm
                .iter()
                .any(|wavelength| !wavelength.is_finite() || *wavelength <= 0.0)
        {
            return Err(invalid(
                "invalid molecular absorption wavelength or temperature",
            ));
        }
        match &self.source {
            CrossSections::Authored { values, .. } => Ok(*values),
            #[cfg(feature = "material-state")]
            CrossSections::Material { card, source } => {
                let mut values = [[0.0; 2]; B];
                for lane in 0..B {
                    values[lane] = source
                        .at_wavelength_temperature(
                            card,
                            finite(wavelengths_nm[lane] / 1e9)?,
                            Some(temperature_k),
                        )
                        .map_err(|_| {
                            invalid("molecular absorption wavelength/temperature is outside original material support")
                        })?
                        .cross_sections_m2();
                }
                Ok(values)
            }
        }
    }
}
