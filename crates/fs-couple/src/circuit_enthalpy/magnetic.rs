//! Shared geometry, caloric and source binding for magnetic thermal consumers.
//!
//! Stranded windings and resolved eddy currents use the same original cell
//! identities, fixed solid volumes, complete material snapshots and permeability
//! admission. This module does not own another field or thermal integrator.

/// Shared field entropy, heat-capacity admission and reversible heat transfer.
pub mod caloric;

use super::{
    CircuitEnthalpyError, finite,
    material::{ConductionLink, MaterialCell, MaterialNetworkError, MaterialNetworkSample},
};
use fs_em::{
    planar::{
        PlanarMagnetostatics,
        eddy::{EddyError, nonlinear::NonlinearEddyError},
    },
    winding::ohmic::OhmicWindingError,
};
use fs_material::{
    conductor::{ConductorError, ResolvedConductorState},
    field_storage::{
        FieldStorageError, MAGNETIC_PERMEABILITY_PROPERTY, admit_fixed_magnetic_refresh,
        admit_magnetic_temperature_refresh,
    },
    state_point::ResolvedMaterialStatePoint,
    thermal::source::ThermalMaterialError,
};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use std::{cell::RefCell, collections::BTreeMap};

/// Exactly one thermal or explicitly fixed correspondence per original triangle.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MagneticRegion {
    Cell {
        triangle: usize,
        cell: CellId,
    },
    /// Unheated region outside the thermal network, at its prepared state.
    Fixed {
        triangle: usize,
    },
}
impl MagneticRegion {
    fn triangle(self) -> usize {
        match self {
            Self::Cell { triangle, .. } | Self::Fixed { triangle } => triangle,
        }
    }
}
#[derive(Clone, Debug, PartialEq)]
pub struct MagneticThermalSample {
    pub triangle: usize,
    pub cell: Option<CellId>,
    pub material: ResolvedMaterialStatePoint,
}
#[derive(Debug)]
pub enum MagneticThermalError {
    Material(MaterialNetworkError),
    Heat(super::contacts::heat::HeatNetworkError),
    Conductor(ConductorError),
    Winding(OhmicWindingError),
    NonlinearWinding(fs_em::winding::nonlinear::NonlinearWindingError),
    Induction(EddyError),
    SpatialInduction(fs_em::spatial::SpatialMagneticError),
    SpatialNonlinearInduction(
        fs_em::spatial::vector_potential::induction::nonlinear::SpatialNonlinearInductionError,
    ),
    NonlinearInduction(NonlinearEddyError),
    Magnetic(FieldStorageError),
    Invalid(&'static str),
}
impl core::fmt::Display for MagneticThermalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material magnetic/enthalpy: {self:?}")
    }
}
impl std::error::Error for MagneticThermalError {}
pub(super) type Error = CircuitEnthalpyError<MagneticThermalError>;
pub(super) fn invalid(what: &'static str) -> Error {
    CircuitEnthalpyError::Source(MagneticThermalError::Invalid(what))
}
pub(super) fn material(error: MaterialNetworkError) -> Error {
    match error {
        MaterialNetworkError::Thermal(DepositionError::Cancelled) => {
            CircuitEnthalpyError::Cancelled
        }
        other => CircuitEnthalpyError::Source(MagneticThermalError::Material(other)),
    }
}
pub(super) fn check(cancelled: &RefCell<impl FnMut() -> bool>) -> Result<(), Error> {
    if (*cancelled.borrow_mut())() {
        Err(CircuitEnthalpyError::Cancelled)
    } else {
        Ok(())
    }
}

/// One data-only layout shared by retained field and circuit field consumers.
pub struct MagneticThermalNetwork<'a> {
    /// Complete original caloric/card/conductivity binding for active cells.
    pub cells: &'a [MaterialCell<'a>],
    pub links: &'a [ConductionLink],
    pub magnetic_regions: &'a [MagneticRegion],
    /// Absolute volume agreement allowance [m3], never an inferred fill factor.
    pub volume_tolerance_m3: f64,
}
pub(super) struct MagneticRoutes {
    pub magnetic: BTreeMap<usize, Option<CellId>>,
    pub conductors: BTreeMap<usize, CellId>,
}
impl MagneticThermalNetwork<'_> {
    pub(super) fn routes(
        &self,
        thermal: &DepositionSimulation<'_>,
        problem: &PlanarMagnetostatics,
        conducting_triangles: &[usize],
        cancelled: &RefCell<impl FnMut() -> bool>,
    ) -> Result<MagneticRoutes, Error> {
        self.routes_for_mesh(thermal, problem.mesh(), conducting_triangles, cancelled)
    }

    // All magnetic laws use the same original triangle/thermal volume map.
    pub(super) fn routes_for_mesh(
        &self,
        thermal: &DepositionSimulation<'_>,
        mesh: &fs_rep_mesh::TriComplex2,
        conducting_triangles: &[usize],
        cancelled: &RefCell<impl FnMut() -> bool>,
    ) -> Result<MagneticRoutes, Error> {
        let volumes = (0..mesh.faces().len())
            .map(|cell| {
                check(cancelled)?;
                mesh.face_measure(cell)
                    .ok_or_else(|| invalid("missing original magnetic triangle geometry"))
            })
            .collect::<Result<Vec<_>, _>>()?;
        self.routes_for_volumes(thermal, &volumes, conducting_triangles, cancelled)
    }

    // Original cell indices and physical volumes are sufficient for thermal
    // routing; planar and tetrahedral owners supply their own admitted geometry.
    pub(super) fn routes_for_volumes(
        &self,
        thermal: &DepositionSimulation<'_>,
        cell_volumes_m3: &[f64],
        conducting_triangles: &[usize],
        cancelled: &RefCell<impl FnMut() -> bool>,
    ) -> Result<MagneticRoutes, Error> {
        if !self.volume_tolerance_m3.is_finite()
            || self.volume_tolerance_m3 < 0.0
            || self.magnetic_regions.len() != cell_volumes_m3.len()
        {
            return Err(invalid(
                "every magnetic triangle needs an explicit thermal/fixed region and a finite volume tolerance",
            ));
        }
        let mut magnetic = BTreeMap::new();
        let mut volumes = BTreeMap::new();
        for &region in self.magnetic_regions {
            check(cancelled)?;
            let triangle = region.triangle();
            let cell = match region {
                MagneticRegion::Cell { cell, .. } => Some(cell),
                MagneticRegion::Fixed { .. } => None,
            };
            if triangle >= cell_volumes_m3.len()
                || cell.is_some_and(|id| !thermal.cells().contains_key(&id))
                || magnetic.insert(triangle, cell).is_some()
            {
                return Err(invalid(
                    "duplicate, foreign or inactive magnetic region correspondence",
                ));
            }
            if let Some(cell) = cell {
                let volume = cell_volumes_m3[triangle];
                if !volume.is_finite() || volume <= 0.0 {
                    return Err(invalid(
                        "original magnetic cell volume must be finite positive",
                    ));
                }
                let sum = volumes.entry(cell).or_insert(0.0);
                *sum = finite(*sum + volume)?;
            }
        }
        let mut conductors = BTreeMap::new();
        for &triangle in conducting_triangles {
            check(cancelled)?;
            let cell = magnetic.get(&triangle).copied().flatten().ok_or_else(|| {
                invalid("every conducting triangle must have an actual thermal heat destination")
            })?;
            if conductors.insert(triangle, cell).is_some() {
                return Err(invalid("repeated conducting triangle"));
            }
        }
        for (id, volume) in volumes {
            check(cancelled)?;
            let cell = &thermal.cells()[&id];
            if (volume - cell.geometry.volume_m3).abs() > self.volume_tolerance_m3 {
                return Err(invalid(
                    "mapped magnetic triangle volumes must cover the whole receiving thermal cell",
                ));
            }
            let binding = self
                .cells
                .iter()
                .find(|s| s.cell == id)
                .ok_or_else(|| invalid("magnetic thermal cell lacks caloric binding"))?;
            for knot in binding.curve.knots() {
                check(cancelled)?;
                if knot.liquid_mass_fraction != 0.0
                    || (finite(cell.mass_kg / knot.bulk_density_kg_m3)? - cell.geometry.volume_m3)
                        .abs()
                        > self.volume_tolerance_m3
                {
                    return Err(invalid(
                        "fixed magnetic geometry requires a fully solid caloric chart with matching constant density",
                    ));
                }
            }
        }
        Ok(MagneticRoutes {
            magnetic,
            conductors,
        })
    }

    pub(super) fn materials_at(
        &self,
        thermal: &DepositionSimulation<'_>,
        problem: &PlanarMagnetostatics,
        conductors: &BTreeMap<usize, &ResolvedConductorState>,
        sample: &MaterialNetworkSample,
        routes: &MagneticRoutes,
        cancelled: &RefCell<impl FnMut() -> bool>,
    ) -> Result<
        (
            Vec<(usize, ResolvedConductorState)>,
            Vec<MagneticThermalSample>,
        ),
        Error,
    > {
        self.materials_at_impl(
            thermal, problem, conductors, sample, routes, cancelled, false,
        )
    }

    // Used only by the thermal-permeability owner, which independently binds
    // every changing triangle to its source derivative and caloric receiver.
    pub(super) fn materials_at_thermal(
        &self,
        thermal: &DepositionSimulation<'_>,
        problem: &PlanarMagnetostatics,
        conductors: &BTreeMap<usize, &ResolvedConductorState>,
        sample: &MaterialNetworkSample,
        routes: &MagneticRoutes,
        cancelled: &RefCell<impl FnMut() -> bool>,
    ) -> Result<
        (
            Vec<(usize, ResolvedConductorState)>,
            Vec<MagneticThermalSample>,
        ),
        Error,
    > {
        self.materials_at_impl(
            thermal, problem, conductors, sample, routes, cancelled, true,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn materials_at_impl(
        &self,
        thermal: &DepositionSimulation<'_>,
        problem: &PlanarMagnetostatics,
        conductors: &BTreeMap<usize, &ResolvedConductorState>,
        sample: &MaterialNetworkSample,
        routes: &MagneticRoutes,
        cancelled: &RefCell<impl FnMut() -> bool>,
        thermal_permeability: bool,
    ) -> Result<
        (
            Vec<(usize, ResolvedConductorState)>,
            Vec<MagneticThermalSample>,
        ),
        Error,
    > {
        let originals: Vec<_> = problem.materials().iter().collect();
        self.materials_at_using(
            thermal,
            &originals,
            conductors,
            sample,
            routes,
            cancelled,
            |_, original, snapshot, binding| {
                let admission = if thermal_permeability {
                    admit_magnetic_temperature_refresh(original, snapshot)
                } else {
                    admit_fixed_magnetic_refresh(original, snapshot)
                };
                admission.map_err(|e| match e {
                    FieldStorageError::Invalid(message) => invalid(message),
                    other => CircuitEnthalpyError::Source(MagneticThermalError::Magnetic(other)),
                })?;
                let property = snapshot
                    .property(MAGNETIC_PERMEABILITY_PROPERTY)
                    .expect("admitted permeability");
                if !thermal_permeability {
                    binding
                        .point
                        .admit_static(property)
                        .map_err(|e| material(MaterialNetworkError::Material(e)))?;
                }
                Ok(())
            },
        )
    }

    // Conductivity, caloric correspondence, fixed geometry and source refresh
    // are common. Only the magnetic law's supported refresh admission differs.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn materials_at_using<C: fs_material::conductor::spatial::SpatialConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'_>,
        originals: &[&ResolvedMaterialStatePoint],
        conductors: &BTreeMap<usize, &C>,
        sample: &MaterialNetworkSample,
        routes: &MagneticRoutes,
        cancelled: &RefCell<impl FnMut() -> bool>,
        admit: impl Fn(
            usize,
            &ResolvedMaterialStatePoint,
            &ResolvedMaterialStatePoint,
            &MaterialCell<'_>,
        ) -> Result<(), Error>,
    ) -> Result<
        (
            Vec<(usize, C)>,
            Vec<MagneticThermalSample>,
        ),
        Error,
    > {
        let bindings: BTreeMap<_, _> = self.cells.iter().map(|c| (c.cell, c)).collect();
        let phases: BTreeMap<_, _> = sample.cells.iter().map(|c| (c.cell, c.phase)).collect();
        let mut result = Vec::with_capacity(conductors.len());
        let mut magnetic = Vec::with_capacity(routes.magnetic.len());
        for (&triangle, &cell) in &routes.magnetic {
            check(cancelled)?;
            let original = *originals
                .get(triangle)
                .ok_or_else(|| invalid("magnetic material changed triangle correspondence"))?;
            let source = conductors.get(&triangle);
            let snapshot = if let Some(id) = cell {
                let binding = bindings
                    .get(&id)
                    .ok_or_else(|| invalid("magnetic cell lacks original material binding"))?;
                let phase = *phases
                    .get(&id)
                    .ok_or_else(|| invalid("magnetic cell lacks actual endpoint phase"))?;
                if !phase.is_fully_solid()
                    || (finite(thermal.cells()[&id].mass_kg / phase.bulk_density_kg_m3())?
                        - thermal.cells()[&id].geometry.volume_m3)
                        .abs()
                        > self.volume_tolerance_m3
                {
                    return Err(invalid("fixed solid magnetic density/geometry changed"));
                }
                let query = binding
                    .point
                    .at(phase)
                    .map_err(|e| material(MaterialNetworkError::Material(e)))?;
                let base = source.map_or(original, |s| s.material());
                let snapshot = base
                    .at_query_point(binding.point.material, &query)
                    .map_err(|e| {
                        material(MaterialNetworkError::Material(
                            ThermalMaterialError::Material(e),
                        ))
                    })?;
                admit(triangle, original, &snapshot, binding)?;
                snapshot
            } else {
                original.clone()
            };
            if let Some(source) = source {
                let conductor = source.with_material(snapshot.clone()).map_err(|e| {
                    CircuitEnthalpyError::Source(MagneticThermalError::Conductor(e))
                })?;
                result.push((triangle, conductor));
            }
            magnetic.push(MagneticThermalSample {
                triangle,
                cell,
                material: snapshot,
            });
        }
        Ok((result, magnetic))
    }
}

// A pinned usage receipt differs from an initial single-claim receipt without
// changing its selected physical source. Compare states without dropping either.
pub(super) fn same_material(
    x: &ResolvedMaterialStatePoint,
    y: &ResolvedMaterialStatePoint,
) -> bool {
    x.card_identity() == y.card_identity()
        && x.material() == y.material()
        && x.axis_quantities() == y.axis_quantities()
        && x.query_point().len() == y.query_point().len()
        && x.query_point()
            .iter()
            .zip(y.query_point())
            .all(|((n, v), (m, w))| n == m && v.to_bits() == w.to_bits())
        && x.properties().len() == y.properties().len()
        && x.properties().iter().zip(y.properties()).all(|(p, q)| {
            p.requirement() == q.requirement()
                && p.source_axes() == q.source_axes()
                && p.answer().receipt.selected == q.answer().receipt.selected
                && p.value_si().to_bits() == q.value_si().to_bits()
        })
}
