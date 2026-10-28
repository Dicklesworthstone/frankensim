//! Sourced sigma(T), resolved surface stress and original thermal flow.
//!
//! The ordinary path uses local caloric cell temperature. A complete solved
//! heat boundary can instead supply source sigma at its interface temperature,
//! extended normally within each admitted surface column. Bulk material states
//! retain their actual caloric temperature. Normal pressure and tangential
//! traction share that surface sigma; neither supplies inferred capillary heat.

use super::*;
use crate::d3q19::freesurface3::thermocapillary::{
    HeightGraphSurface3, LocalSurface3, PlanarSurface3,
};
use crate::material::free_surface::{
    MaterialFreeSurfaceDomain3, MaterialFreeSurfaceInitialCell3, MaterialFreeSurfaceScaling,
    capillary_scaling,
};
use crate::material::transport::boussinesq::BoussinesqConfig;
use fs_exec::Cx;

/// Explicit validity envelope for a nearly flat liquid/gas surface in a tank.
///
/// The liquid occupies every column from one domain wall to its free surface.
/// Internal walls, overhangs, holes, multiple partial cells per column,
/// crossing the declared surface-cell layer and an imposed contact-angle law
/// are outside this reduction. Tangential stress uses the declared Cartesian
/// normal and each column's projected area dx²;
/// a checked height departure admits only the stated flat-surface approximation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlanarThermocapillarity {
    /// Cartesian normal axis, 0=x, 1=y, 2=z.
    pub normal_axis: usize,
    /// +1 when liquid is on the lower-coordinate side; -1 on the upper side.
    pub outward_sign: i8,
    /// Declared surface coordinate in the original grid frame [m].
    pub surface_position_m: f64,
    /// Maximum actual column-height departure [m], nonnegative and < dx/2.
    /// This is a model admission, not an estimated error or a stability proof.
    pub maximum_departure_m: f64,
    /// Positive lower bound (at most one) for the surface cell's liquid fill.
    /// Traction acts on its actual occupied volume; vanishing films refuse.
    pub minimum_surface_fill: f64,
}

/// Finite-slope, single-valued liquid/gas height graph over a Cartesian tank.
///
/// Unlike the planar reduction, the surface may cross lattice layers and
/// contain several cut cells per column. Complete monotone liquid-to-gas
/// columns determine its height; their volume-weighted interface sigma drives
/// both Laplace pressure and the full two-dimensional surface gradient.
///
/// This is a resolved graph reduction, not an arbitrary-drop or overhang model.
/// Internal solid walls, imposed contact angles and a nonzero mass carry refuse.
/// Raw interface fractions must remain in [0,1] and monotone along each column.
/// The native conversion hysteresis can temporarily retain over/underfilled
/// cells; such a candidate is refused by this stricter geometric model.
/// The original Guo force owner distributes each curved patch's traction over
/// its interface band's occupied liquid volume. Bulk heat and material states
/// retain their own temperatures; no artificial capillary heat is introduced.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HeightGraphThermocapillarity {
    /// Cartesian column axis, 0=x, 1=y, 2=z.
    pub normal_axis: usize,
    /// +1 for liquid below the graph, -1 for liquid above it.
    pub outward_sign: i8,
    /// Maximum admitted |grad h|, finite and in [0,1]. A bound of one keeps
    /// the chosen column axis dominant; it is a resolution/model admission.
    pub maximum_slope: f64,
    /// Minimum total occupied volume of the interface band in a column,
    /// measured in cell volumes, positive and at most one. Individual cut
    /// cells may be arbitrarily small because traction is shared by the band.
    pub minimum_surface_fill: f64,
}

#[derive(Clone, Copy)]
pub(super) enum ThermocapillarySurface {
    Planar(PlanarThermocapillarity),
    HeightGraph(HeightGraphThermocapillarity),
    Local(LocalThermocapillarity),
}

/// Locally resolved interfaces, including detached drops and changing graph axes.
///
/// Complete raw-VOF height stencils reconstruct the tangent plane. A bounded
/// plane/cube cut supplies actual patch area; a local affine liquid-source fit
/// extrapolates sigma to that patch and removes its normal gradient from the
/// tangential stress. The reconstructed sigma also drives Laplace pressure.
/// Walls/contact lines, unresolved thin structures, unlocated carry and poorly
/// conditioned source fits refuse. This is not a surfactant or surface-energy
/// evolution model, and cell-source extrapolation requires spatial refinement.
/// Raw interface fractions outside [0,1], including transient native conversion
/// hysteresis states, also refuse; no geometric clamp discards their volume.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LocalThermocapillarity {
    /// Lower bound on a positive cut-cell fill, finite and strictly in (0,1).
    /// Smaller films refuse instead of receiving unbounded area/volume loading.
    pub minimum_surface_fill: f64,
    /// Maximum infinity-norm condition number of the local source-fit matrix.
    /// Finite and at least one; admission is a numerical resolution constraint.
    pub maximum_gradient_condition: f64,
}

impl LocalThermocapillarity {
    fn native(self) -> LocalSurface3 {
        LocalSurface3 {
            minimum_fill: self.minimum_surface_fill,
            maximum_gradient_condition: self.maximum_gradient_condition,
        }
    }
}

impl HeightGraphThermocapillarity {
    pub(crate) fn native(self) -> HeightGraphSurface3 {
        HeightGraphSurface3 {
            axis: self.normal_axis,
            outward_sign: self.outward_sign,
            maximum_slope: self.maximum_slope,
            minimum_fill: self.minimum_surface_fill,
        }
    }
}

impl PlanarThermocapillarity {
    fn native(self, dx: f64) -> Result<PlanarSurface3, Error> {
        let position = self.surface_position_m / dx;
        let maximum_departure = self.maximum_departure_m / dx;
        if !position.is_finite()
            || !maximum_departure.is_finite()
            || self.maximum_departure_m != 0.0 && maximum_departure == 0.0
        {
            return Err(Error::Invalid(
                "planar surface geometry cannot be represented in native cells",
            ));
        }
        Ok(PlanarSurface3 {
            axis: self.normal_axis,
            outward_sign: self.outward_sign,
            position,
            maximum_departure,
            minimum_fill: self.minimum_surface_fill,
        })
    }
}

impl<'s> MaterialFreeSurfaceScaling<'s> {
    /// Evolve source-bound thermocapillary stress on locally resolved interfaces.
    ///
    /// This uses the same native flow, material refresh, mass routes and thermal
    /// clock as the planar and height-graph entries. Every force/pressure field
    /// is reconstructed before and after the native step; geometric or material
    /// refusal publishes none of the candidate state. Complete local stencils
    /// replace the global tank-column requirement, so detached drops and smooth
    /// surfaces spanning different Cartesian orientations are admitted.
    ///
    /// A solved-interface-temperature heat boundary is not supplied by this
    /// entry; sigma is extrapolated from actual bulk material temperatures.
    #[allow(clippy::too_many_arguments)]
    pub fn free_surface3_from_cells_with_local_thermocapillarity<'m>(
        self,
        domain: MaterialFreeSurfaceDomain3<'_>,
        initial: &[MaterialFreeSurfaceInitialCell3],
        law: &'m LiquidTemperatureLaw<'m>,
        caloric: &'m EquilibriumEnthalpyPhaseCurve,
        controls: FreeSurfaceEnthalpyControls,
        surface: LocalThermocapillarity,
        buoyancy: Option<BoussinesqConfig<3>>,
        rheology: Option<rheology::RheologySolve>,
        initial_h: impl Fn([usize; 3]) -> f64,
        cx: &Cx<'_>,
        gate: &CancelGate,
    ) -> Result<MaterialFreeSurfaceEnthalpy3<'s, 'm>, Error> {
        poll(gate)?;
        let binding = rheology
            .map(|solve| rheology::RheologyBinding::try_new(self.material(), law, solve))
            .transpose()?;
        let flow = self.free_surface3_from_cells_inner(domain, initial, cx)?;
        let result = flow.with_enthalpy_model(
            law,
            caloric,
            controls,
            initial_h,
            buoyancy,
            binding,
            Some(ThermocapillarySurface::Local(surface)),
            gate,
        )?;
        cx.checkpoint().map_err(|_| FreeSurfaceError3::Cancelled)?;
        Ok(result)
    }

    /// Bind actual initial liquid, source-varying sigma and optional existing
    /// buoyancy/rheology owners to the same native enthalpy trajectory.
    ///
    /// The complete original material bundle MUST include surface_tension in
    /// N/m, with either a source-declared constant or only the law's temperature
    /// response axis. Every evaluation pins that original claim. The same local
    /// coefficient drives the native Laplace pressure and the surface gradient.
    /// Tangential force is grad_s(sigma) times dx² / (fill*dx³), converted with
    /// the original reference mass and clock. It composes with original gravity,
    /// Boussinesq buoyancy and any later electrical Lorentz loading.
    ///
    /// Actual source temperature is cell-centered, not an independently solved
    /// interface temperature. A first-order held-force split and planar geometry
    /// approximation require refinement. Normal capillary mechanics retain the
    /// native curvature reconstruction. No surface-energy, total-mechanical-energy,
    /// entropy or thermocapillary heat claim is added. Geometry, source, budget
    /// and cancellation refusal publish none of the candidate states.
    #[allow(clippy::too_many_arguments)]
    pub fn free_surface3_from_cells_with_thermocapillarity<'m>(
        self,
        domain: MaterialFreeSurfaceDomain3<'_>,
        initial: &[MaterialFreeSurfaceInitialCell3],
        law: &'m LiquidTemperatureLaw<'m>,
        caloric: &'m EquilibriumEnthalpyPhaseCurve,
        controls: FreeSurfaceEnthalpyControls,
        surface: PlanarThermocapillarity,
        buoyancy: Option<BoussinesqConfig<3>>,
        rheology: Option<rheology::RheologySolve>,
        initial_h: impl Fn([usize; 3]) -> f64,
        cx: &Cx<'_>,
        gate: &CancelGate,
    ) -> Result<MaterialFreeSurfaceEnthalpy3<'s, 'm>, Error> {
        poll(gate)?;
        let binding = rheology
            .map(|solve| rheology::RheologyBinding::try_new(self.material(), law, solve))
            .transpose()?;
        let flow = self.free_surface3_from_cells_inner(domain, initial, cx)?;
        let result = flow.with_enthalpy_model(
            law,
            caloric,
            controls,
            initial_h,
            buoyancy,
            binding,
            Some(ThermocapillarySurface::Planar(surface)),
            gate,
        )?;
        cx.checkpoint().map_err(|_| FreeSurfaceError3::Cancelled)?;
        Ok(result)
    }

    /// Bind sourced sigma(T) to a finite-slope height graph on the original
    /// thermal VOF trajectory. Geometry and surface gradients are reconstructed
    /// at every accepted endpoint; layer crossings need no new declaration.
    ///
    /// The source, caloric, rheology, buoyancy, time-split and no-energy-claim
    /// obligations of the planar constructor apply. Surface sigma is averaged
    /// over the actual liquid interface band in each column. With the
    /// `enthalpy-coupling` feature, `height_graph_interface_patches` exposes the
    /// same normals, centroids and actual areas for solved surface exchange.
    /// That boundary supplies the original sigma at its solved surface
    /// temperature to both traction and pressure; normal conduction distances
    /// remain explicit physical inputs.
    #[allow(clippy::too_many_arguments)]
    pub fn free_surface3_from_cells_with_height_graph_thermocapillarity<'m>(
        self,
        domain: MaterialFreeSurfaceDomain3<'_>,
        initial: &[MaterialFreeSurfaceInitialCell3],
        law: &'m LiquidTemperatureLaw<'m>,
        caloric: &'m EquilibriumEnthalpyPhaseCurve,
        controls: FreeSurfaceEnthalpyControls,
        surface: HeightGraphThermocapillarity,
        buoyancy: Option<BoussinesqConfig<3>>,
        rheology: Option<rheology::RheologySolve>,
        initial_h: impl Fn([usize; 3]) -> f64,
        cx: &Cx<'_>,
        gate: &CancelGate,
    ) -> Result<MaterialFreeSurfaceEnthalpy3<'s, 'm>, Error> {
        poll(gate)?;
        let binding = rheology
            .map(|solve| rheology::RheologyBinding::try_new(self.material(), law, solve))
            .transpose()?;
        let flow = self.free_surface3_from_cells_inner(domain, initial, cx)?;
        let result = flow.with_enthalpy_model(
            law,
            caloric,
            controls,
            initial_h,
            buoyancy,
            binding,
            Some(ThermocapillarySurface::HeightGraph(surface)),
            gate,
        )?;
        cx.checkpoint().map_err(|_| FreeSurfaceError3::Cancelled)?;
        Ok(result)
    }
}

impl MaterialFreeSurfaceEnthalpy3<'_, '_> {
    /// Retained local interface/source reconstruction controls, if selected.
    #[must_use]
    pub const fn local_thermocapillary_surface(&self) -> Option<LocalThermocapillarity> {
        match self.thermocapillary {
            Some(ThermocapillarySurface::Local(surface)) => Some(surface),
            _ => None,
        }
    }

    /// Retained declared planar-interface approximation, if selected.
    #[must_use]
    pub const fn thermocapillary_surface(&self) -> Option<PlanarThermocapillarity> {
        match self.thermocapillary {
            Some(ThermocapillarySurface::Planar(surface)) => Some(surface),
            _ => None,
        }
    }

    /// Retained finite-slope height-graph admission, if selected.
    #[must_use]
    pub const fn height_graph_thermocapillary_surface(
        &self,
    ) -> Option<HeightGraphThermocapillarity> {
        match self.thermocapillary {
            Some(ThermocapillarySurface::HeightGraph(surface)) => Some(surface),
            _ => None,
        }
    }

    #[cfg(feature = "enthalpy-coupling")]
    pub(super) fn temperature_dependent_thermocapillarity(&self) -> bool {
        self.thermocapillary.is_some()
            && self
                .law
                .reference_state()
                .material()
                .property(SURFACE_TENSION_PROPERTY)
                .is_some_and(|sigma| !sigma.source_axes().is_empty())
    }

    pub(super) fn resolve_thermocapillarity(
        &self,
        lattice: &FreeSurface3,
        states: &BTreeMap<MassLocation3, ResolvedLiquidState>,
        forces: &mut Option<BTreeMap<(u64, usize), [f64; 3]>>,
        solved_tensions: Option<BTreeMap<(u64, usize), f64>>,
        gate: &CancelGate,
    ) -> Result<Option<BTreeMap<(u64, usize), f64>>, Error> {
        let Some(surface) = self.thermocapillary else {
            return Ok(None);
        };
        let scaling = &self.flow.scaling;
        let mut tensions =
            if let Some(tensions) = solved_tensions {
                tensions
            } else {
                let mut tensions = BTreeMap::new();
                for (&location, state) in states {
                    poll(gate)?;
                    let MassLocation3::Cell(key, lane) = location else {
                        continue;
                    };
                    let sigma = state.material().property(SURFACE_TENSION_PROPERTY).ok_or(
                        Error::Invalid("evolving surface tension left the shared material bundle"),
                    )?;
                    let native = capillary_scaling(
                        Some(sigma.value_si()),
                        scaling.pressure_unit_pa,
                        scaling.cell_size_m(),
                    )?;
                    tensions.insert((key, lane), native);
                }
                tensions
            };
        let extra = match surface {
            ThermocapillarySurface::Planar(surface) => lattice.planar_marangoni_forces(
                surface.native(scaling.cell_size_m())?,
                &tensions,
                gate,
            )?,
            ThermocapillarySurface::HeightGraph(surface) => {
                lattice.height_graph_marangoni_forces(surface.native(), &mut tensions, gate)?
            }
            ThermocapillarySurface::Local(surface) => {
                lattice.local_marangoni_forces(surface.native(), &mut tensions, gate)?
            }
        };
        if !extra.is_empty() {
            let forces = forces.get_or_insert_with(|| {
                tensions
                    .keys()
                    .map(|&location| (location, lattice.reference_force()))
                    .collect()
            });
            for (location, traction) in extra {
                poll(gate)?;
                let force = forces.get_mut(&location).ok_or(Error::Invalid(
                    "Marangoni traction requires a resolved wet cell",
                ))?;
                for axis in 0..3 {
                    force[axis] += traction[axis];
                    if !force[axis].is_finite() {
                        return Err(Error::Invalid(
                            "combined thermocapillary force is not finite",
                        ));
                    }
                }
            }
        } else if forces.is_none() && self.thermocapillary.is_some() {
            // A former gradient may vanish. Restore the original force frame
            // at every endpoint without adding a second half-force impulse.
            *forces = Some(
                tensions
                    .keys()
                    .map(|&location| (location, lattice.reference_force()))
                    .collect(),
            );
        }
        Ok(Some(tensions))
    }
}
