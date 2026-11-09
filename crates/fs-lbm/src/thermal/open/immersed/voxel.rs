//! Cartesian solid conduction on the original immersed voxel geometry.
//!
//! Each selected voxel owns one thermal node. Neighboring selected cells share
//! the harmonic conductance of their two half cells. Surface resistance adds
//! the solid half-cell path to the explicitly supplied fluid-contact term.
//! The prepared network remains the only temperature and energy owner.

use super::conduction::ImmersedConductionLink3;
use super::{FiniteThermalWall, ImmersedThermalBody3};
use crate::thermal::boundary::ThermalWall;
use crate::thermal::boundary::conduction::WallConductionControls;
use crate::thermal::open::{ThermalOpenError3, ThermalOpenFace3, ThermalOpenField3};

/// Fixed isotropic properties of one ORIGINAL solid voxel, in lattice units.
/// For cubic spacing dx, time step dt and reference fluid volumetric capacity
/// rho_f Cp_f: capacity = rho_s Cp_s/(rho_f Cp_f), conductivity =
/// k_s dt/(rho_f Cp_f dx^2), contact resistance = R_contact rho_f Cp_f dx/dt.
/// A supplied heat rate density maps to q''' dt/(rho_f Cp_f delta_T).
/// All temperatures use the thermal field's same affine coordinate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ImmersedThermalVoxel3 {
    /// Original x-fast native grid index.
    pub cell: usize,
    pub temperature: f64,
    /// TOTAL capacity of this one voxel, relative to one fluid-cell capacity.
    pub heat_capacity: f64,
    /// Positive isotropic conductivity in the units defined above.
    pub conductivity: f64,
    /// Additional nonnegative surface-to-fluid resistance, excluding the solid
    /// half cell. It does not add a resistance between neighboring solid cells.
    pub contact_resistance: f64,
    /// Signed heat supplied to this one voxel per original lattice step.
    pub heat_per_step: f64,
}

/// Fixed grid-aligned orthotropic solid on one original cubic native voxel.
///
/// The diagonal conductivity tensor is `[k_x, k_y, k_z]` in the ORIGINAL grid
/// frame, in the same lattice units as `ImmersedThermalVoxel3`. Both reciprocal
/// solid conduction and the solid half-cell path to each fluid face use that
/// face's normal component. No scalar average replaces the tensor.
/// Rotated tensors with off-diagonal components require a different stencil.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OrthotropicThermalVoxel3 {
    pub cell: usize,
    pub temperature: f64,
    pub heat_capacity: f64,
    /// Three positive finite principal conductivities aligned with x/y/z.
    pub conductivity: [f64; 3],
    /// Additional nonnegative bulk-to-fluid contact resistance by normal axis.
    /// Opposing faces share a value; solid-to-solid contact remains ideal.
    pub contact_resistance: [f64; 3],
    pub heat_per_step: f64,
}

impl From<ImmersedThermalVoxel3> for OrthotropicThermalVoxel3 {
    fn from(voxel: ImmersedThermalVoxel3) -> Self {
        Self {
            cell: voxel.cell,
            temperature: voxel.temperature,
            heat_capacity: voxel.heat_capacity,
            conductivity: [voxel.conductivity; 3],
            contact_resistance: [voxel.contact_resistance; 3],
            heat_per_step: voxel.heat_per_step,
        }
    }
}

/// Separate bounds on temporary geometry, retained bodies, and sparse factors.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ImmersedVoxelControls3 {
    /// Payload of sorted input, returned mapping, declarations and raw faces.
    /// Input storage and allocator metadata are outside this preparation budget.
    pub maximum_geometry_bytes: usize,
    pub maximum_body_storage_bytes: usize,
    pub conduction: WallConductionControls,
}

impl Default for ImmersedVoxelControls3 {
    fn default() -> Self {
        Self {
            maximum_geometry_bytes: 16 * 1024 * 1024,
            maximum_body_storage_bytes: 32 * 1024 * 1024,
            conduction: WallConductionControls::default(),
        }
    }
}

/// Native-cell mapping for the installed network, in canonical x-fast order.
/// It contains no copied temperatures or thermal history. It remains valid
/// until another immersed-body/network installation replaces this network.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ImmersedVoxelMap3 {
    cells: Vec<usize>,
    conduction_links: usize,
}

impl ImmersedVoxelMap3 {
    /// Original cell index for each installed node, in declaration order.
    pub fn cells(&self) -> &[usize] {
        &self.cells
    }
    /// Node index to use with the original `immersed_body` and source setter.
    pub fn node_for_cell(&self, cell: usize) -> Option<usize> {
        self.cells.binary_search(&cell).ok()
    }
    /// Distinct undirected paths after combining periodic parallel faces.
    pub const fn conduction_links(&self) -> usize {
        self.conduction_links
    }
}

fn invalid(message: &'static str) -> ThermalOpenError3 {
    ThermalOpenError3::Invalid(message)
}

fn positive_neighbor(field: &ThermalOpenField3, cell: usize, axis: usize) -> Option<usize> {
    let stride = [1, field.dims[0], field.dims[0] * field.dims[1]][axis];
    let width = field.dims[axis];
    if cell / stride % width + 1 < width {
        Some(cell + stride)
    } else if field.boundaries.faces[2 * axis + 1] == ThermalOpenFace3::Periodic {
        Some(cell - (width - 1) * stride)
    } else {
        None
    }
}

fn harmonic_conductance(a: f64, b: f64) -> f64 {
    let small = a.min(b);
    let large = a.max(b);
    // 1/(0.5/a + 0.5/b), without forming a*b or a+b.
    small * (2.0 / (1.0 + small / large))
}

impl OrthotropicThermalVoxel3 {
    fn surface_resistance(self) -> [f64; 3] {
        std::array::from_fn(|axis| 0.5 / self.conductivity[axis] + self.contact_resistance[axis])
    }

    fn wall(self) -> Result<FiniteThermalWall, ThermalOpenError3> {
        if self
            .conductivity
            .iter()
            .any(|k| !k.is_finite() || *k <= 0.0)
            || self
                .contact_resistance
                .iter()
                .any(|r| !r.is_finite() || *r < 0.0)
        {
            return Err(invalid(
                "voxel conductivity must be positive and contact resistance nonnegative",
            ));
        }
        let resistance = self.surface_resistance();
        if resistance.iter().any(|r| !r.is_finite()) {
            return Err(invalid(
                "immersed voxel half-cell resistance is unrepresentable",
            ));
        }
        let solid = FiniteThermalWall {
            temperature: self.temperature,
            heat_capacity: self.heat_capacity,
            // For uniform faces preserve the original scalar declaration.
            // Otherwise the node has no common path: its authoritative actual
            // face paths are retained in the body's `surface_resistance`.
            resistance: if resistance.iter().all(|r| *r == resistance[0]) {
                resistance[0]
            } else {
                0.0
            },
            heat_per_step: self.heat_per_step,
        };
        ThermalWall::FiniteSolid(solid)
            .validate()
            .map_err(ThermalOpenError3::Wall)?;
        Ok(solid)
    }
}

impl ThermalOpenField3 {
    /// Install spatial solid conduction and fluid exchange on selected native
    /// solid voxels, with one capacity/source per voxel and no invented faces.
    ///
    /// Selected neighbors have ideal solid contact; heterogeneous conductivities
    /// combine through their physical half-cell series resistance. Periodic
    /// seams use the original thermal geometry. Two distinct periodic faces
    /// between the same nodes add in parallel, rather than being dropped as a
    /// duplicate. Degenerate self-neighbors carry no heat and are omitted.
    /// Unselected solid neighbors and nonperiodic exterior solid faces insulate.
    ///
    /// Coefficients are fixed over the prepared network's lifetime. This is a
    /// cubic finite-volume solid with a native stair-stepped fluid interface;
    /// it does not resolve rotated anisotropy, latent heat or subcell geometry.
    /// Interior nodes require a conducting path to an exposed node. All checks
    /// and factor preparation precede replacement of the accepted network.
    pub fn set_immersed_thermal_voxels(
        &mut self,
        voxels: &[ImmersedThermalVoxel3],
        controls: ImmersedVoxelControls3,
    ) -> Result<ImmersedVoxelMap3, ThermalOpenError3> {
        self.prepare_immersed_voxels(voxels.iter().copied().map(Into::into), controls)
    }

    /// Resolve grid-aligned orthotropic solid conduction and actual convective
    /// interfaces on the ORIGINAL stationary voxel mask.
    ///
    /// For an axis-normal solid face the reciprocal conductance is
    /// `1/(0.5/k_a,axis + 0.5/k_b,axis)`. The path to the fluid uses
    /// `0.5/k_axis + R_contact,axis`, followed ONCE by the original D3Q7 fluid
    /// half-cell resistance. Every exposed face retains its own conductance;
    /// one simultaneous native solve advances the voxel's single inventory.
    ///
    /// Same topology, budgets, periodic-face handling and atomic replacement as
    /// the isotropic entry. Fixed diagonal tensors only; no rotated anisotropy,
    /// material-source inference, variable heat capacity or nonlinear k(T).
    pub fn set_immersed_orthotropic_thermal_voxels(
        &mut self,
        voxels: &[OrthotropicThermalVoxel3],
        controls: ImmersedVoxelControls3,
    ) -> Result<ImmersedVoxelMap3, ThermalOpenError3> {
        self.prepare_immersed_voxels(voxels.iter().copied(), controls)
    }

    fn prepare_immersed_voxels(
        &mut self,
        voxels: impl ExactSizeIterator<Item = OrthotropicThermalVoxel3>,
        controls: ImmersedVoxelControls3,
    ) -> Result<ImmersedVoxelMap3, ThermalOpenError3> {
        let n = voxels.len();
        if n == 0 || n > controls.conduction.max_nodes || n > self.solid.len() {
            return Err(invalid(
                "immersed voxel count exceeds the admitted native/node bounds",
            ));
        }
        let record_bytes = core::mem::size_of::<OrthotropicThermalVoxel3>()
            + core::mem::size_of::<usize>()
            + core::mem::size_of::<ImmersedThermalBody3<'_>>();
        let base_bytes = n
            .checked_mul(record_bytes)
            .ok_or_else(|| invalid("immersed voxel geometry size overflow"))?;
        if base_bytes > controls.maximum_geometry_bytes {
            return Err(invalid("immersed voxel geometry storage budget exceeded"));
        }
        let mut ordered = Vec::new();
        ordered
            .try_reserve_exact(n)
            .map_err(|_| ThermalOpenError3::Allocation)?;
        for voxel in voxels {
            if voxel.cell >= self.solid.len() || !self.solid[voxel.cell] {
                return Err(invalid("immersed voxel must name an original solid cell"));
            }
            voxel.wall()?;
            ordered.push(voxel);
        }
        ordered.sort_unstable_by_key(|voxel| voxel.cell);
        if ordered.windows(2).any(|pair| pair[0].cell == pair[1].cell) {
            return Err(invalid(
                "a native solid voxel may own only one thermal node",
            ));
        }
        let mut faces = 0usize;
        for voxel in &ordered {
            for axis in 0..3 {
                if let Some(next) = positive_neighbor(self, voxel.cell, axis)
                    && next != voxel.cell
                    && ordered
                        .binary_search_by_key(&next, |other| other.cell)
                        .is_ok()
                {
                    faces = faces
                        .checked_add(1)
                        .ok_or_else(|| invalid("immersed voxel face count overflow"))?;
                }
            }
        }
        let bytes = faces
            .checked_mul(core::mem::size_of::<ImmersedConductionLink3>())
            .and_then(|value| base_bytes.checked_add(value))
            .ok_or_else(|| invalid("immersed voxel geometry size overflow"))?;
        if bytes > controls.maximum_geometry_bytes {
            return Err(invalid("immersed voxel geometry storage budget exceeded"));
        }
        let mut links = Vec::new();
        links
            .try_reserve_exact(faces)
            .map_err(|_| ThermalOpenError3::Allocation)?;
        for (a, voxel) in ordered.iter().enumerate() {
            for axis in 0..3 {
                if let Some(next) = positive_neighbor(self, voxel.cell, axis)
                    && next != voxel.cell
                    && let Ok(b) = ordered.binary_search_by_key(&next, |other| other.cell)
                {
                    let conductance = harmonic_conductance(
                        voxel.conductivity[axis],
                        ordered[b].conductivity[axis],
                    );
                    if !conductance.is_finite() || conductance <= 0.0 {
                        return Err(invalid("immersed voxel conductance is unrepresentable"));
                    }
                    links.push(ImmersedConductionLink3 {
                        a: a.min(b),
                        b: a.max(b),
                        conductance,
                    });
                }
            }
        }
        links.sort_unstable_by_key(|link| (link.a, link.b));
        let mut distinct = 0;
        for i in 0..links.len() {
            let link = links[i];
            if distinct > 0 && (links[distinct - 1].a, links[distinct - 1].b) == (link.a, link.b) {
                let total = links[distinct - 1].conductance + link.conductance;
                if !total.is_finite() {
                    return Err(invalid(
                        "parallel immersed voxel conductance is unrepresentable",
                    ));
                }
                links[distinct - 1].conductance = total;
            } else {
                links[distinct] = link;
                distinct += 1;
            }
        }
        links.truncate(distinct);
        if distinct > controls.conduction.max_links {
            return Err(invalid("immersed voxel link budget exceeded"));
        }
        let mut cells = Vec::new();
        cells
            .try_reserve_exact(n)
            .map_err(|_| ThermalOpenError3::Allocation)?;
        cells.extend(ordered.iter().map(|voxel| voxel.cell));
        let mut declarations = Vec::new();
        declarations
            .try_reserve_exact(n)
            .map_err(|_| ThermalOpenError3::Allocation)?;
        for (voxel, cell) in ordered.iter().zip(&cells) {
            declarations.push(ImmersedThermalBody3 {
                solid_cells: core::slice::from_ref(cell),
                solid: voxel.wall()?,
            });
        }
        let mut prepared = self.prepare_immersed_bodies_with_surface_resistance(
            &declarations,
            controls.maximum_body_storage_bytes,
            true,
            |node| ordered[node].surface_resistance(),
        )?;
        prepared.conduction = Some(super::conduction::ImmersedConduction3::try_new(
            &prepared.bodies,
            &links,
            controls.conduction,
        )?);
        self.immersed = Some(prepared);
        Ok(ImmersedVoxelMap3 {
            cells,
            conduction_links: distinct,
        })
    }
}

#[cfg(test)]
mod tests;
