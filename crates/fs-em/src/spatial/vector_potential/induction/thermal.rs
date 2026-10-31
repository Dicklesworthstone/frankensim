//! Source-pinned temperature refresh of a retained three-dimensional A–phi state.
//!
//! The original geometry, fixed tangential trace, relative tree gauge and
//! electrical components are retained. Endpoint conductivity and permeability
//! come from the same source selections and temperature coordinate; no field is
//! re-equilibrated or reset. The original midpoint owner supplies signed
//! magnetic parameter work, which a finite thermal owner must account for with
//! its material entropy/internal-energy law. This entry point alone does not
//! choose a thermal temperature, deposit heat or infer a heat-capacity law.

use super::super::super::{SpatialMagneticMaterial, invalid, poll};
use super::{
    MagneticOwner, Prepared, Result, SpatialInduction, SpatialInductionCandidate,
    SpatialInductionControls, coupling::Geometry,
};
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_material::{
    conductor::spatial::SpatialConductorMaterial,
    field_storage::admit_magnetic_temperature_refresh,
};
use std::sync::Arc;

impl<'m, C: SpatialConductorMaterial> SpatialInduction<'m, C> {
    /// Advance the same accepted history with complete endpoint source bundles.
    ///
    /// Magnetic cells must remain original induced-linear materials; their
    /// positive permeability may change only through the source-pinned absolute
    /// temperature coordinate. Conducting support, conductivity/resistivity
    /// basis and selected electrical claims remain fixed. The returned native
    /// report independently closes free energy, Ohmic heat, Lorentz work,
    /// impressed work and signed magnetic parameter work. A caloric consumer
    /// must own that parameter transfer and the field entropy correction.
    ///
    /// Material-weighted magnetic/conductivity operators are reassembled on
    /// shared immutable geometric allocations. Every candidate starts from the
    /// same retained A and source history; failure publishes no partial refresh.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_materials_candidate(
        &self,
        endpoint_current_a_m2: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        duration_s: f64,
        materials: Vec<SpatialMagneticMaterial>,
        conductors: Vec<Option<C>>,
        assembly: WeightedAssemblyLimits,
        controls: SpatialInductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialInductionCandidate<'m, C>> {
        poll(cancelled)?;
        let old = &self.prepared;
        if materials.len() != old.problem.materials().len() || conductors.len() != materials.len() {
            return Err(invalid(
                "complete original spatial thermal source map required",
            ));
        }
        for (before, after) in old.problem.materials().iter().zip(&materials) {
            poll(cancelled)?;
            let (SpatialMagneticMaterial::Linear(before), SpatialMagneticMaterial::Linear(after)) =
                (before, after)
            else {
                return Err(invalid(
                    "temperature refresh requires original induced-linear magnetic cells",
                ));
            };
            admit_magnetic_temperature_refresh(before, after)?;
        }
        let problem = Arc::new(old.problem.with_materials(materials, assembly, cancelled)?);
        let coupling = old.coupling.refresh_sources(
            Geometry {
                geometry: &problem.geometry,
                cells: &problem.cells,
                curl: &problem.curl,
            },
            conductors,
            |cell, source| {
                let Some(original) = old.sources[cell].as_ref() else {
                    return Err(invalid("thermal refresh changed conducting support"));
                };
                let SpatialMagneticMaterial::Linear(material) = &problem.materials()[cell] else {
                    unreachable!("linear source map admitted above")
                };
                if source.material() != material
                    || !original.has_same_electrical_source(source)
                    || original
                        .with_material(source.material().clone())
                        .map_err(|_| {
                            invalid("thermal refresh left its original electrical source")
                        })?
                        != *source
                {
                    return Err(invalid(
                        "thermal conductor changed its original material, basis or source",
                    ));
                }
                Ok(())
            },
            cancelled,
        )?;
        let magnetic = coupling.restrict(&problem.stiffness, cancelled)?;
        let prepared = Arc::new(Prepared {
            problem: MagneticOwner::Owned(problem),
            coupling,
            magnetic,
        });
        self.step_prepared_candidate(
            prepared,
            endpoint_current_a_m2,
            velocity_m_s,
            duration_s,
            controls,
            cancelled,
        )
    }
}
