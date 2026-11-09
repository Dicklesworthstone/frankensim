//! Planar rest states in the original density-dependent Guo force frame.

use super::{ContactModel, FreeSurface, FreeSurfacePressureError2, neighbor, pressure};
use crate::core2::{Cell, CollisionModel2, Grid};
use crate::{CS2, Q, equilibrium};
use fs_exec::CancelGate;

fn invalid(reason: &'static str) -> FreeSurfacePressureError2 {
    FreeSurfacePressureError2::Invalid { coordinates: None, reason }
}

/// One recurrence for native D2Q9 and material-bound D2Q9/D3Q19.
/// Callers retain their own cancellation, physical-pressure and error owners.
pub(crate) fn density_layers<Error>(
    acceleration_lattice: f64,
    depth: usize,
    maximum_density_excursion: f64,
    mut checkpoint: impl FnMut() -> Result<(), Error>,
    mut admit_density: impl FnMut(f64) -> Result<(), Error>,
    invalid: impl Fn(&'static str) -> Error,
) -> Result<Vec<f64>, Error> {
    // For u=0, the original collision maps w*rho-S/2 to w*rho+S/2.
    // Streaming one cell downwards therefore requires
    // rho_next*(1-h) = rho_previous*(1+h). At a missing gas link,
    // f_in = 2*w*rho_boundary - f_out instead gives
    // rho_nearest*(1-h) = rho_boundary = 1. This is not sqrt(ratio).
    let half_head = acceleration_lattice / (2.0 * CS2);
    if !half_head.is_finite() || half_head >= 1.0 {
        return Err(invalid(
            "hydrostatic lattice acceleration cannot retain positive rest populations",
        ));
    }
    let ratio = (1.0 + half_head) / (1.0 - half_head);
    let mut density = 1.0 / (1.0 - half_head);
    let mut layers = Vec::new();
    layers.try_reserve_exact(depth).map_err(|_| {
        invalid("hydrostatic density-profile allocation failed")
    })?;
    for layer in 0..depth {
        checkpoint()?;
        if !density.is_finite() || (density - 1.0).abs() > maximum_density_excursion {
            return Err(invalid(
                "hydrostatic profile exceeds its declared density-excursion bound",
            ));
        }
        admit_density(density)?;
        layers.push(density);
        if layer + 1 < depth {
            density *= ratio;
        }
    }
    Ok(layers)
}

impl FreeSurface {
    /// Initialize complete liquid columns in an existing closed wall/weir mask.
    ///
    /// Supply a fresh [`Grid::uniform`] with only Fluid/Gas/Wall flags changed,
    /// nonperiodic boundaries, and downward `grid.g = [0, -acceleration]`.
    /// Each occupied column must contain Fluid from row 1 through
    /// `liquid_top_y`, with gas immediately above. Gas links may cross only
    /// that common horizontal free surface; exposed vertical sides refuse.
    /// Exterior wall rows and the original interior walls remain unchanged.
    ///
    /// The constructor uses the same density recurrence as the material-bound
    /// hydrostatic owners, then the original Guo force-frame transform to
    /// establish physical zero velocity. It uses BGK, neutral walls and zero
    /// capillarity; it does not settle or advance the liquid. The ordinary
    /// pressure boundary retains gas density one.
    ///
    /// Interface cells retain FULL mass. The represented surface lies halfway
    /// between the last wet and first gas cell centres. Initial inventory thus
    /// includes every declared layer and its bounded numerical compression;
    /// it is not renormalized to the generic constructor's half-filled state.
    /// This initializes a new liquid, not a conservative remap of an old one.
    ///
    /// # Errors
    /// Refuses incomplete/periodic grids, incompatible geometry or forcing,
    /// existing interface/evolved populations, non-BGK collision, invalid
    /// relaxation, an empty column or insufficient gas space, a density
    /// excursion outside the finite bound in `(0,1)`, or cancellation. No
    /// partially initialized flow is returned. Zero gravity is supported.
    pub fn new_hydrostatic_column(
        mut grid: Grid,
        liquid_top_y: usize,
        maximum_density_excursion: f64,
        gate: &CancelGate,
    ) -> Result<Self, FreeSurfacePressureError2> {
        pressure::check_cancel(Some(gate))?;
        let (nx, ny) = (grid.nx, grid.ny);
        let n = nx.checked_mul(ny).ok_or_else(|| invalid("invalid hydrostatic grid extent"))?;
        if nx < 3 || ny < 4 || grid.periodic_x || grid.periodic_y
            || grid.flags.len() != n || grid.f.len() != n
            || grid.tau.len() != n || grid.fext.len() != n
        {
            return Err(invalid("hydrostatic initialization requires a complete closed nonperiodic grid"));
        }
        if liquid_top_y == 0 || liquid_top_y >= ny - 2 {
            return Err(invalid("hydrostatic column needs positive liquid depth and at least one interior gas layer"));
        }
        if !maximum_density_excursion.is_finite()
            || maximum_density_excursion <= 0.0 || maximum_density_excursion >= 1.0
        {
            return Err(invalid("hydrostatic density-excursion bound must be finite and in (0,1)"));
        }
        if grid.collision != CollisionModel2::Bgk {
            return Err(invalid("masked hydrostatic initialization requires the original BGK owner"));
        }
        if grid.g[0] != 0.0 || !grid.g[1].is_finite() || grid.g[1] > 0.0 {
            return Err(invalid("hydrostatic gravity must point toward the lower supporting wall"));
        }
        let mut occupied = false;
        let fresh = equilibrium(1.0, 0.0, 0.0);
        for y in 0..ny {
            pressure::check_cancel(Some(gate))?;
            for x in 0..nx {
                let i = grid.idx(x, y);
                let flag = grid.flags[i];
                if grid.fext[i] != [0.0; 2] {
                    return Err(invalid("hydrostatic initialization requires zero external force"));
                }
                if x == 0 || y == 0 || x + 1 == nx || y + 1 == ny {
                    if flag != Cell::Wall {
                        return Err(invalid("hydrostatic exterior boundary must consist of wall cells"));
                    }
                    continue;
                }
                if flag == Cell::Interface {
                    return Err(invalid("hydrostatic initialization requires a fresh Fluid/Gas/Wall mask"));
                }
                if flag != Cell::Wall && (!grid.tau[i].is_finite() || grid.tau[i] <= 0.5) {
                    return Err(invalid("hydrostatic relaxation must be finite and greater than one half"));
                }
                let column_wet = grid.flags[grid.idx(x, liquid_top_y)] == Cell::Fluid;
                if column_wet && y <= liquid_top_y && flag != Cell::Fluid {
                    return Err(invalid("hydrostatic wet columns must contain complete layers above the floor"));
                }
                if flag != Cell::Fluid {
                    continue;
                }
                occupied = true;
                if !column_wet || y > liquid_top_y || grid.f[i] != fresh {
                    return Err(invalid("hydrostatic initialization requires fresh unit-density columns below one common surface"));
                }
                if y == liquid_top_y && grid.flags[grid.idx(x, y + 1)] != Cell::Gas {
                    return Err(invalid("every hydrostatic wet column must expose the declared upper gas face"));
                }
                for q in 1..Q {
                    if let Some(j) = neighbor(&grid, x, y, q)
                        && grid.flags[j] == Cell::Gas
                        && (y != liquid_top_y || j / nx != liquid_top_y + 1)
                    {
                        return Err(invalid("hydrostatic liquid cannot expose a vertical or submerged gas face"));
                    }
                }
            }
        }
        if !occupied {
            return Err(invalid("hydrostatic initialization requires a nonempty liquid column"));
        }
        let gravity = grid.g;
        let layers = density_layers(
            gravity[1].abs(), liquid_top_y, maximum_density_excursion,
            || pressure::check_cancel(Some(gate)), |_| Ok(()), invalid,
        )?;
        // This is a fresh state: install its rest density before binding the
        // held force. Reusing the original reframe avoids a second Guo owner.
        grid.g = [0.0; 2];
        for y in 1..=liquid_top_y {
            pressure::check_cancel(Some(gate))?;
            let density = layers[liquid_top_y - y];
            for x in 1..nx - 1 {
                let i = grid.idx(x, y);
                if grid.flags[i] == Cell::Fluid && density != 1.0 {
                    grid.f[i] = equilibrium(density, 0.0, 0.0);
                }
            }
        }
        grid.rebind_gravity(gravity, &mut Vec::new()).map_err(|_| {
            invalid("D2Q9 gravity cannot represent the physical rest state")
        })?;
        let mut flow = Self::new(grid, 0.0, ContactModel::Neutral);
        for x in 1..nx - 1 {
            pressure::check_cancel(Some(gate))?;
            let i = flow.grid.idx(x, liquid_top_y);
            if flow.grid.flags[i] == Cell::Interface {
                flow.mass[i] = flow.grid.f[i].iter().sum();
            }
        }
        pressure::check_cancel(Some(gate))?;
        Ok(flow)
    }
}
