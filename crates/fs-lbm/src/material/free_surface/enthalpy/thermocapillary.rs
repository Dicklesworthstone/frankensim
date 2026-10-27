//! Sourced sigma(T), planar tangential stress and the original thermal flow.
//!
//! Source temperature is the actual local caloric cell temperature, held over
//! one native flow interval. No surface temperature, surfactant transport,
//! surface entropy or capillary heat model is inferred. Normal pressure and
//! tangential traction have distinct owners and neither is added twice.

use super::*;
use crate::d3q19::freesurface3::thermocapillary::PlanarSurface3;
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
    /// native smoothed-fill accuracy. No surface-energy, total-mechanical-energy,
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
            Some(surface),
            gate,
        )?;
        cx.checkpoint().map_err(|_| FreeSurfaceError3::Cancelled)?;
        Ok(result)
    }
}

impl MaterialFreeSurfaceEnthalpy3<'_, '_> {
    /// Retained declared planar-interface approximation, if selected.
    #[must_use]
    pub const fn thermocapillary_surface(&self) -> Option<PlanarThermocapillarity> {
        self.thermocapillary
    }

    pub(super) fn resolve_thermocapillarity(
        &self,
        lattice: &FreeSurface3,
        states: &BTreeMap<MassLocation3, ResolvedLiquidState>,
        forces: &mut Option<BTreeMap<(u64, usize), [f64; 3]>>,
        gate: &CancelGate,
    ) -> Result<Option<BTreeMap<(u64, usize), f64>>, Error> {
        let Some(surface) = self.thermocapillary else {
            return Ok(None);
        };
        let scaling = &self.flow.scaling;
        let mut tensions = BTreeMap::new();
        for (&location, state) in states {
            poll(gate)?;
            let MassLocation3::Cell(key, lane) = location else {
                continue;
            };
            let sigma =
                state
                    .material()
                    .property(SURFACE_TENSION_PROPERTY)
                    .ok_or(Error::Invalid(
                        "evolving surface tension left the shared material bundle",
                    ))?;
            let native = capillary_scaling(
                Some(sigma.value_si()),
                scaling.pressure_unit_pa,
                scaling.cell_size_m(),
            )?;
            tensions.insert((key, lane), native);
        }
        let extra = lattice.planar_marangoni_forces(
            surface.native(scaling.cell_size_m())?,
            &tensions,
            gate,
        )?;
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
