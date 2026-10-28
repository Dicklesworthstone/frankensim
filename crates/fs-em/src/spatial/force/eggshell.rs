//! Distributed virtual-displacement (eggshell) magnetic force and torque.
//!
//! A continuous original-mesh P1 weight w equals one on the rigid target and
//! zero on the exterior boundary. Its transition lies in one homogeneous,
//! induced-linear, current-free exterior. With the original matched B/H,
//! F=-integral(T grad(w)). For torque the virtual motion is the ORIGINAL P1
//! interpolant sum_i N_i*w_i*(omega cross (x_i-o)), which is exactly rigid
//! on the weight-one plateau and zero outside. Integrating its complete
//! gradient avoids a quadratic rotation field inconsistent with the P1 mesh.
//!
//! Original cell fields and the virtual-motion gradient are constant; original
//! volumes integrate both functionals exactly, without field smoothing or new solves.
//! This spreads the functional through a collar instead of selecting one
//! discontinuous face trace. Different weights remain a discretization check,
//! not an error bound. No source/control functional shape adjoint, equilibrated
//! stress, moving-domain energy exchange or material striction is claimed.
//!
//! Henrotte, Deliége, Hameyer, "The eggshell approach for the computation of
//! electromagnetic forces in 2D and 3D", Eq. (24):
//! https://orbi.uliege.be/bitstream/2268/29071/1/pub_1132.pdf .

use super::super::{
    Face, SpatialMagneticError, SpatialMagneticField, SpatialMagneticMaterial,
    SpatialMagnetostatics, cross, dot, finite, invalid,
    nonlinear::{SpatialNonlinearField, SpatialNonlinearMagnetostatics, SpatialNonlinearMaterial},
    poll, sub,
    vector_potential::{SpatialVectorField, SpatialVectorMagnetostatics},
};
use super::SpatialMaxwellForceLimits;
use fs_feec::whitney::ElementGeometry;
use fs_ga::{Se3, Vec3, Wrench};
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_rep_mesh::TetComplex;

type Result<T> = core::result::Result<T, SpatialMagneticError>;

struct Cell {
    index: usize,
    /// Original cell volume times the P1 weight gradient [m2].
    volume_gradient_m2: [f64; 3],
    centroid_from_reference_m: [f64; 3],
    /// V sum_i (w_i-w_0)(x_i-centroid) tensor grad(N_i) [m3].
    /// The omitted constant-weight identity contracts to zero against the
    /// symmetric stress in a rigid rotation.
    moment_gradient_m3: [[f64; 3]; 3],
}

/// Cached geometric functional tied to one immutable original magnetic mesh.
///
/// A weight-one plateau declares the target (possibly disconnected), including
/// its attached exterior region; every exterior vertex must have weight zero.
/// Intermediate plateaus and all changing weights are admitted as exterior
/// material at application time. A magnetic/current-carrying body cannot hide
/// in a fractional plateau and contribute a fractional unreported force.
///
/// The preparation retains original geometry identity, so it cannot be applied
/// to a different problem merely because its topology or dimensions match.
/// It can be reused across every accepted excitation of that same problem.
pub struct PreparedSpatialEggshell<'m> {
    mesh: &'m TetComplex,
    positions_m: &'m [[f64; 3]],
    weights: Vec<f64>,
    cells: Vec<Cell>,
    reference_m: [f64; 3],
    plateau_cells: usize,
    collar_volume_m3: f64,
}
impl PreparedSpatialEggshell<'_> {
    pub fn nodal_weights(&self) -> &[f64] {
        &self.weights
    }
    /// Includes intermediate constant plateaus, which also require exterior admission.
    pub fn collar_cell_count(&self) -> usize {
        self.cells.len()
    }
    pub const fn plateau_cell_count(&self) -> usize {
        self.plateau_cells
    }
    pub const fn collar_volume_m3(&self) -> f64 {
        self.collar_volume_m3
    }
}

fn prepare<'m>(
    mesh: &'m TetComplex,
    positions_m: &'m [[f64; 3]],
    geometry: &ElementGeometry,
    faces: &[Face],
    weights: &[f64],
    limits: SpatialMaxwellForceLimits,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<PreparedSpatialEggshell<'m>> {
    poll(cancelled)?;
    if mesh.tets.len() > limits.max_cells
        || faces.len() > limits.max_faces
        || weights.len() != positions_m.len()
    {
        return Err(invalid(
            "complete original-mesh eggshell weights and mesh budgets required",
        ));
    }
    for &weight in weights {
        poll(cancelled)?;
        if !weight.is_finite() || !(0.0..=1.0).contains(&weight) {
            return Err(invalid(
                "eggshell weights must be finite and between zero and one",
            ));
        }
    }
    for (index, face) in faces.iter().enumerate() {
        poll(cancelled)?;
        if face.cells.len() == 1
            && mesh.faces[index]
                .iter()
                .any(|&v| weights[v as usize] != 0.0)
        {
            return Err(invalid(
                "eggshell virtual displacement must vanish on the entire exterior",
            ));
        }
    }
    let mut first_plateau = None;
    let mut plateau_cells = 0;
    for (cell, tet) in mesh.tets.iter().enumerate() {
        poll(cancelled)?;
        if tet.iter().all(|&v| weights[v as usize] == 1.0) {
            first_plateau.get_or_insert(cell);
            plateau_cells += 1;
        }
    }
    let first = first_plateau
        .ok_or_else(|| invalid("eggshell requires a complete weight-one target cell"))?;
    let reference_m = positions_m[mesh.tets[first][0] as usize];
    let mut cells = Vec::new();
    let mut collar_volume_m3 = 0.0;
    for (cell, tet) in mesh.tets.iter().enumerate() {
        poll(cancelled)?;
        let w = tet.map(|v| weights[v as usize]);
        if w.iter().all(|&x| x == 0.0) || w.iter().all(|&x| x == 1.0) {
            continue;
        }
        let mut gradient = [0.0; 3];
        for local in 1..4 {
            let difference = w[local] - w[0];
            for axis in 0..3 {
                gradient[axis] =
                    finite(difference.mul_add(geometry.grads[cell][local][axis], gradient[axis]))?;
            }
        }
        let volume = geometry.vol_signed[cell].abs();
        let volume_gradient_m2 = gradient.map(|x| x * volume);
        if w.iter().any(|&x| x != w[0]) && volume_gradient_m2.iter().all(|&x| x == 0.0) {
            return Err(invalid("eggshell gradient is not representable"));
        }
        let offset = tet.map(|v| sub(positions_m[v as usize], reference_m));
        let centroid_from_reference_m = std::array::from_fn(|axis| {
            offset[0][axis] / 4.0
                + offset[1][axis] / 4.0
                + offset[2][axis] / 4.0
                + offset[3][axis] / 4.0
        });
        let mut moment_gradient_m3 = [[0.0; 3]; 3];
        for local in 1..4 {
            for position in 0..3 {
                let arm = offset[local][position] - centroid_from_reference_m[position];
                let weight = finite(volume * (w[local] - w[0]) * arm)?;
                for gradient in 0..3 {
                    moment_gradient_m3[position][gradient] = finite(weight.mul_add(
                        geometry.grads[cell][local][gradient],
                        moment_gradient_m3[position][gradient],
                    ))?;
                }
            }
        }
        for x in volume_gradient_m2
            .into_iter()
            .chain(centroid_from_reference_m)
        {
            finite(x)?;
        }
        collar_volume_m3 = finite(collar_volume_m3 + volume)?;
        cells.push(Cell {
            index: cell,
            volume_gradient_m2,
            centroid_from_reference_m,
            moment_gradient_m3,
        });
    }
    if cells.is_empty() {
        return Err(invalid("eggshell requires a nonempty exterior transition"));
    }
    poll(cancelled)?;
    Ok(PreparedSpatialEggshell {
        mesh,
        positions_m,
        weights: weights.to_vec(),
        cells,
        reference_m,
        plateau_cells,
        collar_volume_m3,
    })
}

impl SpatialMagnetostatics {
    /// Prepare a reusable P1 virtual-displacement collar; field/source admission
    /// remains mandatory when its wrench is evaluated.
    pub fn prepare_eggshell(
        &self,
        weights: &[f64],
        limits: SpatialMaxwellForceLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedSpatialEggshell<'_>> {
        prepare(
            &self.mesh,
            &self.positions_m,
            &self.geometry,
            &self.faces,
            weights,
            limits,
            cancelled,
        )
    }
}
impl SpatialNonlinearMagnetostatics {
    /// Same weight, exterior and geometry contract as the linear preparation.
    pub fn prepare_eggshell(
        &self,
        weights: &[f64],
        limits: SpatialMaxwellForceLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedSpatialEggshell<'_>> {
        prepare(
            &self.geometry.mesh,
            &self.geometry.positions_m,
            &self.geometry.geometry,
            &self.geometry.faces,
            weights,
            limits,
            cancelled,
        )
    }
}
impl SpatialVectorMagnetostatics {
    /// Same geometric preparation; impressed current must vanish throughout
    /// the entire transition when the field is applied.
    pub fn prepare_eggshell(
        &self,
        weights: &[f64],
        limits: SpatialMaxwellForceLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedSpatialEggshell<'_>> {
        prepare(
            &self.geometry.mesh,
            &self.geometry.positions_m,
            &self.geometry.geometry,
            &self.geometry.faces,
            weights,
            limits,
            cancelled,
        )
    }
}

impl super::super::vector_potential::nonlinear::SpatialNonlinearVectorMagnetostatics {
    /// Same geometric preparation; impressed current must vanish throughout
    /// the entire transition when the field is applied.
    pub fn prepare_eggshell(
        &self,
        weights: &[f64],
        limits: SpatialMaxwellForceLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedSpatialEggshell<'_>> {
        prepare(
            &self.geometry.mesh,
            &self.geometry.positions_m,
            &self.geometry.geometry,
            &self.geometry.faces,
            weights,
            limits,
            cancelled,
        )
    }
}

fn admit_source<'a>(
    source: &'a ResolvedMaterialStatePoint,
    exterior: &mut Option<&'a ResolvedMaterialStatePoint>,
) -> Result<()> {
    let reference = *exterior.get_or_insert(source);
    if source != reference {
        return Err(invalid(
            "eggshell requires one complete homogeneous original exterior source",
        ));
    }
    Ok(())
}
impl SpatialMagneticField<'_> {
    /// Raw volume-functional load in the receiving frame and about its declared
    /// origin. No finite-mesh equality with a surface integral is assumed.
    pub fn eggshell_wrench(
        &self,
        shell: &PreparedSpatialEggshell<'_>,
        placement: &Se3,
        moment_origin_m: [f64; 3],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Wrench> {
        let p = self.problem();
        let mut exterior = None;
        apply(
            shell,
            p.mesh(),
            p.positions_m(),
            placement,
            moment_origin_m,
            |cell| {
                let SpatialMagneticMaterial::Linear(source) = &p.materials()[cell] else {
                    return Err(invalid(
                        "eggshell transition must be induced-linear exterior material",
                    ));
                };
                admit_source(source, &mut exterior)?;
                Ok((self.flux_density_t()[cell], self.field_strength_a_m()[cell]))
            },
            cancelled,
        )
    }
}
impl SpatialNonlinearField<'_> {
    /// Saturation/recoil may occur inside the weight-one plateau, never inside
    /// the exterior transition or an intermediate constant plateau.
    pub fn eggshell_wrench(
        &self,
        shell: &PreparedSpatialEggshell<'_>,
        placement: &Se3,
        moment_origin_m: [f64; 3],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Wrench> {
        let p = self.problem();
        let mut exterior = None;
        apply(
            shell,
            p.mesh(),
            p.positions_m(),
            placement,
            moment_origin_m,
            |cell| {
                let SpatialNonlinearMaterial::Linear(source) = &p.materials()[cell] else {
                    return Err(invalid(
                        "eggshell transition must be induced-linear exterior material",
                    ));
                };
                admit_source(source, &mut exterior)?;
                Ok((self.flux_density_t()[cell], self.field_strength_a_m()[cell]))
            },
            cancelled,
        )
    }
}
impl SpatialVectorField<'_> {
    /// Source-identical induced-linear and current-free transition only.
    pub fn eggshell_wrench(
        &self,
        shell: &PreparedSpatialEggshell<'_>,
        placement: &Se3,
        moment_origin_m: [f64; 3],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Wrench> {
        let p = self.problem();
        let mut exterior = None;
        apply(
            shell,
            p.mesh(),
            p.positions_m(),
            placement,
            moment_origin_m,
            |cell| {
                if self.impressed_current_density_a_m2()[cell]
                    .iter()
                    .any(|&j| j != 0.0)
                {
                    return Err(invalid(
                        "eggshell transition must be free of impressed current",
                    ));
                }
                let SpatialMagneticMaterial::Linear(source) = &p.materials()[cell] else {
                    return Err(invalid(
                        "eggshell transition must be induced-linear exterior material",
                    ));
                };
                admit_source(source, &mut exterior)?;
                Ok((self.flux_density_t()[cell], self.field_strength_a_m()[cell]))
            },
            cancelled,
        )
    }
}

impl super::super::vector_potential::nonlinear::SpatialNonlinearVectorField<'_> {
    /// Source-identical induced-linear and current-free transition only.
    pub fn eggshell_wrench(
        &self,
        shell: &PreparedSpatialEggshell<'_>,
        placement: &Se3,
        moment_origin_m: [f64; 3],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Wrench> {
        let p = self.problem();
        let mut exterior = None;
        apply(
            shell,
            p.mesh(),
            p.positions_m(),
            placement,
            moment_origin_m,
            |cell| {
                if self.impressed_current_density_a_m2()[cell]
                    .iter()
                    .any(|&j| j != 0.0)
                {
                    return Err(invalid(
                        "eggshell transition must be free of impressed current",
                    ));
                }
                let SpatialNonlinearMaterial::Linear(source) = &p.materials()[cell] else {
                    return Err(invalid(
                        "eggshell transition must be induced-linear exterior material",
                    ));
                };
                admit_source(source, &mut exterior)?;
                Ok((self.flux_density_t()[cell], self.field_strength_a_m()[cell]))
            },
            cancelled,
        )
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn apply(
    shell: &PreparedSpatialEggshell<'_>,
    mesh: &TetComplex,
    positions_m: &[[f64; 3]],
    placement: &Se3,
    moment_origin_m: [f64; 3],
    mut fields: impl FnMut(usize) -> Result<([f64; 3], [f64; 3])>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Wrench> {
    poll(cancelled)?;
    if !std::ptr::eq(shell.mesh, mesh)
        || !std::ptr::eq(shell.positions_m, positions_m)
        || moment_origin_m.iter().any(|x| !x.is_finite())
    {
        return Err(invalid(
            "eggshell belongs to another original geometry or has an invalid origin",
        ));
    }
    let mut sum = [0.0; 6];
    let mut correction = [0.0; 6];
    for cell in &shell.cells {
        poll(cancelled)?;
        let (b, h) = fields(cell.index)?;
        for x in b.into_iter().chain(h) {
            finite(x)?;
        }
        let g = cell.volume_gradient_m2;
        let pressure = finite(0.5 * dot(b, h))?;
        let normal_h = finite(dot(h, g))?;
        let force = std::array::from_fn(|axis| (-b[axis]).mul_add(normal_h, pressure * g[axis]));
        let mut torque = cross(cell.centroid_from_reference_m, force);
        // The product w*(omega cross x) is quadratic. Its nodal P1 interpolant
        // has an additional gradient moment; omitting it destroys the actual
        // discrete virtual-work torque even though the continuum limits agree.
        let mut moment_stress = [[0.0; 3]; 3];
        for (position, row) in cell.moment_gradient_m3.iter().enumerate() {
            let projection = finite(dot(h, *row))?;
            for force in 0..3 {
                moment_stress[position][force] =
                    finite(b[force].mul_add(projection, -pressure * row[force]))?;
            }
        }
        torque[0] = finite(torque[0] + moment_stress[2][1] - moment_stress[1][2])?;
        torque[1] = finite(torque[1] + moment_stress[0][2] - moment_stress[2][0])?;
        torque[2] = finite(torque[2] + moment_stress[1][0] - moment_stress[0][1])?;
        for (index, value) in torque.into_iter().chain(force).enumerate() {
            let adjusted = finite(value - correction[index])?;
            let next = finite(sum[index] + adjusted)?;
            correction[index] = finite((next - sum[index]) - adjusted)?;
            sum[index] = next;
        }
    }
    let frame_error = |_| invalid("eggshell wrench frame transformation is unrepresentable");
    let force = placement
        .transform_vector(Vec3::new(sum[3], sum[4], sum[5]))
        .map_err(frame_error)?;
    let torque = placement
        .transform_vector(Vec3::new(sum[0], sum[1], sum[2]))
        .map_err(frame_error)?;
    let reference = placement
        .transform_point(Vec3::new(
            shell.reference_m[0],
            shell.reference_m[1],
            shell.reference_m[2],
        ))
        .map_err(frame_error)?;
    let shift = Vec3::new(
        reference.x - moment_origin_m[0],
        reference.y - moment_origin_m[1],
        reference.z - moment_origin_m[2],
    )
    .cross(force);
    let result = Wrench::new(
        Vec3::new(torque.x + shift.x, torque.y + shift.y, torque.z + shift.z),
        force,
    );
    for x in result.to_array() {
        finite(x)?;
    }
    poll(cancelled)?;
    Ok(result)
}
