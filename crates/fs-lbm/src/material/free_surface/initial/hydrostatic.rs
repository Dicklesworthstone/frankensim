//! Rest states of the original density-dependent Guo/VOF operator.

use super::*;
use std::collections::BTreeSet;

fn checkpoint(cx: &Cx<'_>) -> Result<(), MaterialLbmError> {
    cx.checkpoint()
        .map_err(|_| crate::d3q19::freesurface3::FreeSurfaceError3::Cancelled.into())
}

/// Shared D2Q9/D3Q19 rest profile, ordered from the free surface into the liquid.
/// The caller admits its domain, pressure datum and density-excursion bound.
pub(in crate::material::free_surface) fn hydrostatic_density_layers(
    scaling: &MaterialFreeSurfaceScaling<'_>,
    acceleration_lattice: f64,
    depth: usize,
    maximum_density_excursion: f64,
    checkpoint: impl FnMut() -> Result<(), MaterialLbmError>,
) -> Result<Vec<f64>, MaterialLbmError> {
    crate::freesurface::hydrostatic::density_layers(
        acceleration_lattice,
        depth,
        maximum_density_excursion,
        checkpoint,
        |density| scaling.absolute_pressure_si(density).map(|_| ()),
        MaterialLbmError::Invalid,
    )
}

impl<'a> MaterialFreeSurfaceScaling<'a> {
    /// Initialize a planar liquid column at rest under its actual `rho * g` force.
    ///
    /// `supporting_wall` is the domain face toward which gravity points. The
    /// liquid occupies `liquid_depth_cells` complete layers next to that face,
    /// except for the original permanent solid tiles. Its flat surface lies
    /// between the last wet cell center and the first gas cell center, at the
    /// physical depth `liquid_depth_cells * cell_size_m` from the support.
    /// There must be at least one gas layer. Gravity may be zero; otherwise it
    /// must be parallel to, and point toward, the selected Cartesian wall.
    ///
    /// Requires [`Self::with_pressure_reference`]. That positive absolute
    /// pressure [Pa] is the gas pressure at the free surface. Ordinary
    /// [`MaterialFreeSurface3::step`] retains this same pressure boundary;
    /// no pressure callback or relaxation run is needed. A change of this
    /// datum alone leaves populations, pressure differences and loads unchanged.
    /// Capillarity and contact-angle models must be explicitly absent. Solid
    /// tiles remain stationary obstacles; this initializer supplies no moving
    /// geometry, partially wetted body loads or contact-line force model.
    ///
    /// The profile balances the existing collision, streaming and gas-link
    /// reconstruction to population roundoff. With `a = |g| dt^2/dx`,
    /// `h = a/(2 cs^2)` and `n` full layers below the nearest wet layer,
    /// its lattice density is `rho(n) = (1+h)^n / (1-h)^(n+1)`.
    /// The factor `1/(1-h)` at the free surface is essential: a continuum
    /// exponential evaluated at the half-cell position does not satisfy the
    /// native pressure reconstruction. The original initializer supplies the
    /// physical-zero-velocity Guo frame; no second collision is introduced.
    ///
    /// This is a well-balanced *weakly compressible numerical state*, not a
    /// liquid equation of state. Pressure is still
    /// `(rho-1) rho_source (dx/dt)^2 cs^2` above the surface datum, and native
    /// mass includes this density perturbation. Every layer must satisfy the
    /// caller's `maximum_density_excursion` in `(0,1)`. Reducing `dt/dx` reduces
    /// the compressibility error relative to `rho_source |g| depth`; reducing
    /// `dx` alone at fixed `dt/dx` does not remove that physical-model error.
    /// The existing viscosity, Mach and allocation bounds remain necessary.
    ///
    /// # Errors
    /// Refuses incompatible gravity, an empty or fully filled column, missing
    /// pressure datum, capillarity/wetting, invalid physical extents, excessive
    /// or unrepresentable densities/pressures, allocation failure, native input
    /// refusal or cancellation. Construction returns only a complete candidate
    /// and does not change any existing fluid or material source.
    pub fn free_surface3_hydrostatic(
        self,
        domain: MaterialFreeSurfaceDomain3<'_>,
        supporting_wall: PlanarWall3,
        liquid_depth_cells: usize,
        maximum_density_excursion: f64,
        cx: &Cx<'_>,
    ) -> Result<MaterialFreeSurface3<'a>, MaterialLbmError> {
        checkpoint(cx)?;
        self.require_frozen_source()?;
        admit_domain(&domain, &|_, _, _| false)?;
        if self.reference_pressure_pa.is_none() {
            return Err(MaterialLbmError::Invalid(
                "hydrostatic initialization requires an explicit free-surface pressure reference",
            ));
        }
        if !self.density_dependent_gravity {
            return Err(MaterialLbmError::Invalid(
                "hydrostatic initialization requires density-dependent material gravity",
            ));
        }
        if self.surface_tension_lattice != 0.0
            || self.planar_wetting.is_some()
            || domain.contact != ContactModel3::Neutral
        {
            return Err(MaterialLbmError::Invalid(
                "planar hydrostatic initialization requires absent capillarity and neutral walls",
            ));
        }
        if !maximum_density_excursion.is_finite()
            || maximum_density_excursion <= 0.0
            || maximum_density_excursion >= 1.0
        {
            return Err(MaterialLbmError::Invalid(
                "hydrostatic density-excursion bound must be finite and in (0,1)",
            ));
        }
        let (axis, upper) = match supporting_wall {
            PlanarWall3::XMin => (0, false),
            PlanarWall3::XMax => (0, true),
            PlanarWall3::YMin => (1, false),
            PlanarWall3::YMax => (1, true),
            PlanarWall3::ZMin => (2, false),
            PlanarWall3::ZMax => (2, true),
        };
        let extent = domain.cells[axis];
        if liquid_depth_cells == 0 || liquid_depth_cells >= extent {
            return Err(MaterialLbmError::Invalid(
                "hydrostatic column needs positive liquid depth and at least one gas layer",
            ));
        }
        for n in domain.cells {
            positive(
                n as f64 * self.cell_size_m,
                "hydrostatic physical domain extent is unrepresentable",
            )?;
        }
        let gravity = self.gravity_vector_lattice(domain.gravity_m_s2)?;
        if gravity
            .iter()
            .enumerate()
            .any(|(a, &g)| a != axis && g != 0.0)
            || (upper && gravity[axis] < 0.0)
            || (!upper && gravity[axis] > 0.0)
        {
            return Err(MaterialLbmError::Invalid(
                "hydrostatic gravity must point toward its Cartesian supporting wall",
            ));
        }

        let layers = hydrostatic_density_layers(
            &self,
            gravity[axis].abs(),
            liquid_depth_cells,
            maximum_density_excursion,
            || checkpoint(cx),
        )?;

        // Resolve the one-dimensional profile once, then reuse it for every
        // cell in a layer. The original native constructor still owns sparse
        // activation, full-interface fill, mass closure and Guo reframing.
        let mut lower = [0; 3];
        let mut higher = domain.cells;
        if upper {
            lower[axis] = extent - liquid_depth_cells;
        } else {
            higher[axis] = liquid_depth_cells;
        }
        let capacity = (higher[0] - lower[0]) * (higher[1] - lower[1]) * (higher[2] - lower[2]);
        let mut initial = Vec::new();
        initial
            .try_reserve_exact(capacity)
            .map_err(|_| MaterialLbmError::Invalid("hydrostatic initial-cell allocation failed"))?;
        let solids: BTreeSet<_> = domain.solid_tiles.iter().copied().collect();
        for z in lower[2]..higher[2] {
            for y in lower[1]..higher[1] {
                for x in lower[0]..higher[0] {
                    checkpoint(cx)?;
                    if solids.contains(&((x / TILE) as u32, (y / TILE) as u32, (z / TILE) as u32)) {
                        continue;
                    }
                    initial.push(MaterialFreeSurfaceInitialCell3 {
                        coordinates: [x, y, z],
                        fill_fraction: 1.0,
                        velocity_m_s: [0.0; 3],
                    });
                }
            }
        }
        self.free_surface3_from_cells_with_density(domain, &initial, cx, |coordinates| {
            let from_wall = if upper {
                extent - 1 - coordinates[axis]
            } else {
                coordinates[axis]
            };
            Ok(layers[liquid_depth_cells - 1 - from_wall])
        })
    }
}
