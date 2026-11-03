//! The same sourced rho/mu/sigma and wetting state drive the dense D2Q9 owner.
//! Two-dimensional mass is reported per unit out-of-plane depth (kg/m).
//! The flow is translationally invariant in that direction; no second surface
//! curvature, finite-depth wall friction or invented 3-D droplet is included.

use super::{MaterialFreeSurfaceScaling, MaterialLbmError, PlanarWall3, finite, positive};
use crate::core2::{Cell, Grid};
use crate::freesurface::{ContactModel, FreeSurface, FreeSurfacePressureError2, GasPressureLink2, PlanarWall2};
use fs_exec::CancelGate;

mod wetting;
pub use wetting::MaterialFreeSurfaceWetting2;

/// Physical coordinates of one actual native D2Q9 liquid/gas link.
/// The domain origin is zero. The midpoint is halfway between cell centres,
/// not a reconstructed interface intersection or geometric surface normal.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialGasPressureLink2 {
    /// Original native cell indices and incoming population direction.
    pub link: GasPressureLink2,
    /// Liquid cell centre [m].
    pub liquid_center_m: [f64; 2],
    /// Gas cell centre [m].
    pub gas_center_m: [f64; 2],
    /// Geometric midpoint of these two cell centres [m].
    pub midpoint_m: [f64; 2],
}

/// Closed stationary rectangular D2Q9 domain, including its wall rows.
#[derive(Clone, Copy, Debug)]
pub struct MaterialFreeSurfaceDomain2 {
    /// Cell counts including one wall cell at each end of both axes.
    pub cells: [usize; 2],
    /// Gravity in m/s²; the native dense owner uses density-dependent rho*g.
    pub gravity_m_s2: [f64; 2],
    /// Explicit legacy wall model; sourced angles require Neutral elsewhere.
    pub contact: ContactModel,
    /// Bound on all allocated cells, including walls and gas.
    pub max_cells: usize,
}

/// Native wet-cell observation in the original source's SI units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialFreeSurfaceSample2 {
    /// Native fluid/interface classification.
    pub cell: Cell,
    /// Native liquid fill fraction.
    pub fill_fraction: f64,
    /// Numerical liquid density around the source reference, kg/m³.
    pub density_kg_m3: f64,
    /// Guo half-force-corrected velocity, m/s.
    pub velocity_m_s: [f64; 2],
    /// Incompressible perturbation about reference density one, Pa.
    pub pressure_perturbation_pa: f64,
}

/// Original dense free-surface solver, with immutable source and unit binding.
pub struct MaterialFreeSurface2<'a> {
    lattice: FreeSurface,
    scaling: MaterialFreeSurfaceScaling<'a>,
    mass_per_depth_unit_kg_m: f64,
}

impl<'a> MaterialFreeSurfaceScaling<'a> {
    /// Construct a closed D2Q9 free surface with the same liquid and interface
    /// sources as the sparse path. Use `with_planar_wetting` with an X or Y
    /// face to supply a measured or Young–Dupré-derived equilibrium angle.
    /// Z faces refuse: no silent dimensional remapping is performed.
    ///
    /// The classifier receives only interior cell coordinates. Interface
    /// insertion follows the native half-filled-cell convention. This is an
    /// isothermal equilibrium-angle closure, not a dynamic wetting model.
    pub fn free_surface2(
        self,
        domain: MaterialFreeSurfaceDomain2,
        fluid: impl Fn(usize, usize) -> bool,
    ) -> Result<MaterialFreeSurface2<'a>, MaterialLbmError> {
        self.require_frozen_source()?;
        let planar = match self.planar_wetting {
            None => None,
            Some((wall, wetting)) => {
                if domain.contact != ContactModel::Neutral {
                    return Err(MaterialLbmError::Invalid(
                        "a sourced planar angle requires Neutral on all other walls",
                    ));
                }
                let wall = match wall {
                    PlanarWall3::XMin => PlanarWall2::XMin,
                    PlanarWall3::XMax => PlanarWall2::XMax,
                    PlanarWall3::YMin => PlanarWall2::YMin,
                    PlanarWall3::YMax => PlanarWall2::YMax,
                    PlanarWall3::ZMin | PlanarWall3::ZMax => {
                        return Err(MaterialLbmError::Invalid(
                            "a 2-D free surface has no Z-normal wetting wall",
                        ));
                    }
                };
                Some((wall, wetting.contact_angle_radians()))
            }
        };
        let [nx, ny] = domain.cells;
        super::super::admit_cells(nx, ny, domain.max_cells)?;
        if nx < 4 || ny < 4 {
            return Err(MaterialLbmError::Invalid(
                "closed 2-D free surface needs at least two interior cells per axis",
            ));
        }
        let gravity = [
            self.gravity_lattice(domain.gravity_m_s2[0])?,
            self.gravity_lattice(domain.gravity_m_s2[1])?,
        ];
        let mass_per_depth_unit_kg_m = positive(
            (self.material.density_kg_m3() * self.cell_size_m) * self.cell_size_m,
            "2-D free-surface mass per depth is unrepresentable",
        )?;
        let mut grid = Grid::uniform(nx, ny, self.flow.tau);
        grid.periodic_x = false;
        grid.periodic_y = false;
        grid.g = gravity;
        for y in 0..ny {
            for x in 0..nx {
                let i = grid.idx(x, y);
                grid.flags[i] = if x == 0 || y == 0 || x + 1 == nx || y + 1 == ny {
                    Cell::Wall
                } else if fluid(x, y) {
                    Cell::Fluid
                } else {
                    Cell::Gas
                };
            }
        }
        let mut lattice = FreeSurface::new(grid, self.surface_tension_lattice, domain.contact);
        if let Some((wall, angle)) = planar {
            lattice.set_planar_contact_angle(wall, angle).map_err(|_| {
                MaterialLbmError::Invalid("2-D planar contact-angle stencil refused")
            })?;
        }
        Ok(MaterialFreeSurface2 {
            lattice,
            scaling: self,
            mass_per_depth_unit_kg_m,
        })
    }
}

impl<'a> MaterialFreeSurface2<'a> {
    /// Read-only native populations, geometry, mass and conversion accounting.
    #[must_use]
    pub const fn lattice(&self) -> &FreeSurface {
        &self.lattice
    }

    /// The original liquid, interface and physical unit owner.
    #[must_use]
    pub const fn scaling(&self) -> &MaterialFreeSurfaceScaling<'a> {
        &self.scaling
    }

    /// One native dense step at the admitted physical timestep. Preserves the
    /// dense solver's original panic and mass-conversion semantics; unlike the
    /// sparse owner this entry does not add cancellation or rollback support.
    pub fn step(&mut self) {
        self.lattice.step();
    }

    /// Advance the original dense owner under absolute external pressure [Pa]
    /// sampled at exposed liquid cell centres [m]. Pressure is frozen for this
    /// timestep and converted about the explicit construction-time reference.
    /// The combined gas/capillary density departure must fit the supplied bound
    /// in (0,1); this is not a physical equation of state or error certificate.
    ///
    /// The original material coefficients, mass-per-depth ledger and conversion
    /// owner are retained. Pressure/geometry/cancellation refusal publishes
    /// nothing. Gas storage and pressure-work heating are outside this frozen
    /// isothermal liquid model, as for the sparse 3-D pressure adapter.
    pub fn step_with_gas_pressure(
        &mut self,
        gate: &CancelGate,
        maximum_density_excursion: f64,
        pressure_pa: impl Fn([f64; 2]) -> f64,
    ) -> Result<(), MaterialLbmError> {
        let reference = self.scaling.reference_pressure_pa.ok_or(MaterialLbmError::Invalid(
            "absolute gas loading requires an explicit fixed pressure reference",
        ))?;
        let dx = self.scaling.cell_size_m;
        let unit = self.scaling.pressure_unit_pa;
        self.lattice.step_with_gas_pressure(gate, maximum_density_excursion, |coordinates| {
            let position = coordinates.map(|c| (c as f64 + 0.5) * dx);
            pressure_perturbation(reference, unit, pressure_pa(position), coordinates, &position)
        })?;
        Ok(())
    }

    /// Prescribe distinct absolute gas pressures at the actual missing links.
    /// Uniform loads reproduce `step_with_gas_pressure` exactly. Both centres
    /// and the midpoint are supplied in metres; the source decides its spatial
    /// sampling. The same source binding, pressure bound and atomicity apply.
    pub fn step_with_gas_link_pressure(
        &mut self,
        gate: &CancelGate,
        maximum_density_excursion: f64,
        pressure_pa: impl Fn(MaterialGasPressureLink2) -> f64,
    ) -> Result<(), MaterialLbmError> {
        let reference = self.scaling.reference_pressure_pa.ok_or(MaterialLbmError::Invalid(
            "absolute gas loading requires an explicit fixed pressure reference",
        ))?;
        let dx = self.scaling.cell_size_m;
        let unit = self.scaling.pressure_unit_pa;
        self.lattice.step_with_gas_link_pressure(gate, maximum_density_excursion, |link| {
            let liquid_center_m = link.liquid_cell.map(|c| (c as f64 + 0.5) * dx);
            let gas_center_m = link.gas_cell.map(|c| (c as f64 + 0.5) * dx);
            let midpoint_m = core::array::from_fn(|d| f64::midpoint(liquid_center_m[d], gas_center_m[d]));
            let physical = MaterialGasPressureLink2 { link, liquid_center_m, gas_center_m, midpoint_m };
            pressure_perturbation(
                reference, unit, pressure_pa(physical), link.liquid_cell,
                &[liquid_center_m[0], liquid_center_m[1], gas_center_m[0], gas_center_m[1], midpoint_m[0], midpoint_m[1]],
            )
        })?;
        Ok(())
    }

    /// Conserved mass per unit out-of-plane depth, kg/m, including carry.
    pub fn ledger_mass_kg_per_m(&self) -> Result<f64, MaterialLbmError> {
        finite(
            self.lattice.ledger_mass() * self.mass_per_depth_unit_kg_m,
            "2-D free-surface mass ledger is unrepresentable",
        )
    }

    /// Unregioned conversion carry per unit out-of-plane depth, kg/m.
    pub fn carry_mass_kg_per_m(&self) -> Result<f64, MaterialLbmError> {
        finite(
            self.lattice.carry * self.mass_per_depth_unit_kg_m,
            "2-D free-surface carry is unrepresentable",
        )
    }

    /// Physical observation of a liquid/interface cell; gas, walls and outside
    /// the domain remain None. This does not certify the entire flow field.
    pub fn sample_at_si(
        &self,
        x: usize,
        y: usize,
    ) -> Result<Option<MaterialFreeSurfaceSample2>, MaterialLbmError> {
        let grid = &self.lattice.grid;
        if x >= grid.nx || y >= grid.ny {
            return Ok(None);
        }
        let i = grid.idx(x, y);
        let cell = grid.flags[i];
        if !matches!(cell, Cell::Fluid | Cell::Interface) {
            return Ok(None);
        }
        let moments = grid.moments(i);
        Ok(Some(MaterialFreeSurfaceSample2 {
            cell,
            fill_fraction: finite(self.lattice.fill(i), "2-D free-surface fill must be finite")?,
            density_kg_m3: positive(
                moments.rho * self.scaling.material.density_kg_m3(),
                "2-D free-surface density must be positive and finite",
            )?,
            velocity_m_s: [
                self.scaling.velocity_si(moments.u[0])?,
                self.scaling.velocity_si(moments.u[1])?,
            ],
            pressure_perturbation_pa: self.scaling.pressure_perturbation_si(moments.rho)?,
        }))
    }
}

fn pressure_perturbation(
    reference: f64,
    unit: f64,
    absolute: f64,
    coordinates: [usize; 2],
    positions: &[f64],
) -> Result<f64, FreeSurfacePressureError2> {
    let perturbation = (absolute - reference) / unit;
    if !absolute.is_finite() || absolute <= 0.0 || !perturbation.is_finite()
        || positions.iter().any(|x| !x.is_finite())
    {
        return Err(FreeSurfacePressureError2::Invalid {
            coordinates: Some(coordinates),
            reason: "absolute gas pressure and physical coordinates must admit a finite positive load",
        });
    }
    Ok(perturbation)
}
