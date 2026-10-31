//! Distributed stranded windings on a fixed planar magnetic problem.
//!
//! Each profile is signed J_z / terminal_current [1/m2] in triangle order.
//! It already includes turns and the outgoing/return orientation. Area times
//! this density integrates to a signed turn count; volume times A_z times the
//! density is flux LINKAGE [Wb], including the problem's actual axial length.
//! Each profile must balance outgoing and return current on EVERY connected
//! cross-section component (each has its own constant-potential gauge).
//! End closure is outside this 2-D model: callers explicitly supply nonnegative
//! series end/leakage self-inductances [H], including zero when neglected.
//!
//! Unit-current solves reuse the existing FEEC operator. The symmetric magnetic
//! energy Gram matrix supplies active-section self/mutual inductance. Every
//! entry is compared with the independently evaluated A.J_unit pairing before
//! publication; the residual is an estimate, not a discretization certificate.
//! Geometry/material state stays fixed. No moving winding, nonlinear magnetic
//! law, skin/proximity loss, end mutual coupling or temperature-work law is
//! inferred. The optional `ohmic` owner derives intrinsic series resistance
//! from the same spatial material/profile data and retains explicit graph losses.

pub mod circuit;
/// Differential self/mutual inductance at a solved nonlinear operating point.
pub mod incremental;
/// Circuit-driven reversible nonlinear fields with spatial conductor losses.
pub mod nonlinear;
/// Spatial material-derived winding resistance and native-circuit Joule heat.
pub mod ohmic;

use crate::{
    planar::{MagneticError, PlanarField, PlanarMagnetostatics, poll},
    transfer::{CellSample, PlanarFrame},
};
use fs_la::factor::cholesky;
use fs_solver::CheckedCgConfig;

/// Bound the small dense coupled magnetic block independently of mesh size.
pub const MAX_WINDINGS: usize = 16;

/// Allocation and measured consistency tolerances for winding extraction.
#[derive(Clone, Copy, Debug)]
pub struct WindingConfig {
    /// Maximum profiles, in 1..=MAX_WINDINGS.
    pub max_windings: usize,
    /// Ceiling on retained scalar basis values: windings*(2*vertices+7*cells),
    /// plus 2*vertices+12*cells+windings for a retained fixed-remanence bias.
    /// These include each solved field's original impressed-current snapshot
    /// and the bias field's recoil operating-point margins.
    pub max_basis_values: usize,
    /// Absolute source/energy flux-pairing agreement tolerance [H].
    pub pairing_tolerance_h: f64,
    /// Absolute outgoing-minus-return current tolerance [turns].
    pub current_balance_tolerance_turns: f64,
    /// Maximum row-sum residual of L*L_inverse-I; finite, in (0,1).
    pub inverse_residual_tolerance: f64,
}

/// Source, physical-model or bounded numerical refusal.
#[derive(Debug)]
pub enum WindingError {
    /// Existing material/FEEC/field-solve refusal.
    Magnetic(MagneticError),
    /// Invalid winding, correspondence, limits or unrepresentable arithmetic.
    Invalid(&'static str),
    /// Unit-field energy and original source pairing do not agree.
    Pairing {
        /// Source-current basis index.
        source: usize,
        /// Receiving winding index.
        receiver: usize,
        /// Direct integral of A_source.J_receiver [H].
        flux_h: f64,
        /// Magnetic energy Gram entry [H].
        energy_h: f64,
    },
    /// The magnetic terminal block is singular or its inverse was not admitted.
    SingularInductance,
    /// Existing Kirchhoff/consistent-initial-state refusal.
    Circuit(fs_phs::circuit::CircuitError),
    /// Existing descriptor/quadratic-storage admission refusal.
    Phs(fs_phs::PhsError),
}
impl core::fmt::Display for WindingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "planar winding: {self:?}")
    }
}
impl std::error::Error for WindingError {}
impl From<MagneticError> for WindingError {
    fn from(error: MagneticError) -> Self {
        Self::Magnetic(error)
    }
}
impl From<fs_phs::circuit::CircuitError> for WindingError {
    fn from(error: fs_phs::circuit::CircuitError) -> Self {
        Self::Circuit(error)
    }
}
impl From<fs_phs::PhsError> for WindingError {
    fn from(error: fs_phs::PhsError) -> Self {
        Self::Phs(error)
    }
}

fn finite(value: f64) -> Result<f64, WindingError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(WindingError::Invalid(
            "winding arithmetic is not representable",
        ))
    }
}

/// Fixed-mesh/material unit-current fields and their reciprocal magnetic block.
/// Original problem, complete material snapshots, signed winding data and field
/// residual reports remain reachable.
/// Incremental circuit storage is delta_lambda^T L^-1 delta_lambda/2.
/// Fixed remanence additionally supplies total lambda=lambda_0+delta_lambda,
/// a source-admitted zero-current field, and its constant storage reference W0.
pub struct PlanarWindingBasis<'a> {
    problem: &'a PlanarMagnetostatics,
    prescribed_vertices: Vec<usize>,
    profiles_per_m2: Vec<Vec<f64>>,
    fields: Vec<PlanarField<'a>>,
    bias_field: Option<PlanarField<'a>>,
    bias_flux_linkages_wb: Vec<f64>,
    maximum_bias_stationarity_defect_wb: f64,
    active_inductance_h: Vec<f64>,
    inductance_h: Vec<f64>,
    inverse_inductance_per_h: Vec<f64>,
    end_inductance_h: Vec<f64>,
    net_turns_by_component: Vec<Vec<f64>>,
    maximum_pairing_defect_h: f64,
    inverse_residual: f64,
}
impl<'a> PlanarWindingBasis<'a> {
    /// Solve one 1-A excitation per profile with HOMOGENEOUS prescribed A_z.
    /// Boundary vertex selection uses the planar owner's complete-edge and
    /// connected-component admission. Background fields/permanent magnets and
    /// nonzero boundary bias are separate source-energy problems, not an L.
    ///
    /// End/leakage entries are caller-declared series self-inductances; zero
    /// explicitly neglects unrepresented ends. No extra leakage is invented to
    /// regularize dependent winding modes. A singular total block refuses.
    pub fn new(
        problem: &'a PlanarMagnetostatics,
        profiles_per_m2: Vec<Vec<f64>>,
        prescribed_vertices: &[usize],
        end_inductance_h: Vec<f64>,
        config: WindingConfig,
        linear: CheckedCgConfig,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, WindingError> {
        Self::new_impl(
            problem,
            profiles_per_m2,
            prescribed_vertices,
            end_inductance_h,
            config,
            None,
            linear,
            cancelled,
        )
    }

    /// Extract an affine fixed-remanence winding response. A homogeneous
    /// zero-current field must first satisfy every original recoil envelope.
    /// Unit-current derivatives need not themselves be admissible physical
    /// magnet states; actual reconstructed fields are checked as a whole.
    ///
    /// The explicit Wb allowance bounds the measured stationarity defect
    /// integral H0.B_unit dV. It is distinct from the H-valued inductance
    /// pairing allowance. Geometry, Br, direction and material state stay fixed.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_fixed_remanence(
        problem: &'a PlanarMagnetostatics,
        profiles_per_m2: Vec<Vec<f64>>,
        prescribed_vertices: &[usize],
        end_inductance_h: Vec<f64>,
        config: WindingConfig,
        bias_stationarity_tolerance_wb: f64,
        linear: CheckedCgConfig,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, WindingError> {
        if !bias_stationarity_tolerance_wb.is_finite() || bias_stationarity_tolerance_wb < 0.0 {
            return Err(WindingError::Invalid(
                "fixed-remanence stationarity allowance must be finite and nonnegative",
            ));
        }
        Self::new_impl(
            problem,
            profiles_per_m2,
            prescribed_vertices,
            end_inductance_h,
            config,
            Some(bias_stationarity_tolerance_wb),
            linear,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn new_impl(
        problem: &'a PlanarMagnetostatics,
        profiles_per_m2: Vec<Vec<f64>>,
        prescribed_vertices: &[usize],
        end_inductance_h: Vec<f64>,
        config: WindingConfig,
        bias_stationarity_tolerance_wb: Option<f64>,
        linear: CheckedCgConfig,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, WindingError> {
        poll(cancelled)?;
        let has_bias = problem.has_permanent_magnets();
        if has_bias && bias_stationarity_tolerance_wb.is_none() {
            return Err(WindingError::Invalid(
                "permanent-magnet basis requires explicit fixed-remanence stationarity admission",
            ));
        }
        let count = profiles_per_m2.len();
        let cells = problem.mesh().faces().len();
        let vertices = problem.mesh().vertices().len();
        let thickness_m = match problem.mesh().metric() {
            fs_rep_mesh::Metric2::Planar { thickness } => thickness,
            fs_rep_mesh::Metric2::Axisymmetric { .. } => unreachable!("admitted planar problem"),
        };
        let values = cells
            .checked_mul(7)
            .and_then(|n| vertices.checked_mul(2).and_then(|v| n.checked_add(v)))
            .and_then(|n| n.checked_mul(count))
            .and_then(|n| {
                if has_bias {
                    cells
                        .checked_mul(12)
                        .and_then(|c| vertices.checked_mul(2).and_then(|v| c.checked_add(v)))
                        .and_then(|bias| bias.checked_add(count))
                        .and_then(|bias| n.checked_add(bias))
                } else {
                    Some(n)
                }
            });
        if count == 0
            || count > config.max_windings
            || config.max_windings == 0
            || config.max_windings > MAX_WINDINGS
            || values.is_none_or(|n| n > config.max_basis_values)
            || end_inductance_h.len() != count
            || prescribed_vertices.len() > vertices
            || !config.pairing_tolerance_h.is_finite()
            || config.pairing_tolerance_h < 0.0
            || !config.current_balance_tolerance_turns.is_finite()
            || config.current_balance_tolerance_turns < 0.0
            || !(config.inverse_residual_tolerance.is_finite()
                && config.inverse_residual_tolerance > 0.0
                && config.inverse_residual_tolerance < 1.0)
            || end_inductance_h.iter().any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(WindingError::Invalid(
                "winding count, storage, end model or tolerances are invalid",
            ));
        }
        let net_turns_by_component = admit_profiles(
            problem.mesh(), &profiles_per_m2, thickness_m,
            config.current_balance_tolerance_turns, cancelled,
        )?;
        let boundary: Vec<_> = prescribed_vertices.iter().map(|&v| (v, 0.0)).collect();
        let bias_field = if has_bias {
            Some(problem.solve(&vec![0.0; cells], &boundary, linear, cancelled)?)
        } else {
            None
        };
        let mut fields = Vec::with_capacity(count);
        for profile in &profiles_per_m2 {
            poll(cancelled)?;
            fields.push(problem.solve_response(profile, &boundary, linear, cancelled)?);
        }
        let mut bias_flux_linkages_wb = vec![0.0; count];
        let mut maximum_bias_stationarity_defect_wb = 0.0_f64;
        if let Some(bias) = &bias_field {
            for i in 0..count {
                bias_flux_linkages_wb[i] = pairing(
                    problem,
                    bias.potential_wb_m(),
                    &profiles_per_m2[i],
                    cancelled,
                )?;
                let mut stationarity = 0.0;
                for cell in 0..cells {
                    poll(cancelled)?;
                    let h = bias.magnetic_field_a_m()[cell];
                    let b = fields[i].magnetic_flux_density_t()[cell];
                    let derivative = finite(h[0].mul_add(b[0], h[1] * b[1]))?;
                    stationarity = finite(
                        problem
                            .mesh()
                            .face_measure(cell)
                            .unwrap()
                            .mul_add(derivative, stationarity),
                    )?;
                }
                maximum_bias_stationarity_defect_wb =
                    maximum_bias_stationarity_defect_wb.max(stationarity.abs());
                if stationarity.abs()
                    > bias_stationarity_tolerance_wb.expect("admitted fixed remanence")
                {
                    return Err(WindingError::Invalid(
                        "zero-current magnet field exceeds the fixed-remanence stationarity allowance",
                    ));
                }
            }
        }
        let mut active = vec![0.0; count * count];
        for i in 0..count {
            for j in i..count {
                poll(cancelled)?;
                let mut value = 0.0;
                for cell in 0..cells {
                    poll(cancelled)?;
                    let a = fields[i].magnetic_flux_density_t()[cell];
                    let b = fields[j].magnetic_flux_density_t()[cell];
                    let volume = problem.mesh().face_measure(cell).unwrap();
                    let density =
                        finite(a[0].mul_add(b[0], a[1] * b[1]) / problem.permeability_h_m()[cell])?;
                    value = finite(volume.mul_add(density, value))?;
                }
                active[i * count + j] = value;
                active[j * count + i] = value;
            }
            if active[i * count + i] <= 0.0 {
                return Err(WindingError::Invalid(
                    "winding has no representable resolved magnetic energy",
                ));
            }
        }
        let mut maximum_pairing_defect_h = 0.0_f64;
        for i in 0..count {
            for j in 0..count {
                poll(cancelled)?;
                let flux = pairing(
                    problem,
                    fields[i].potential_wb_m(),
                    &profiles_per_m2[j],
                    cancelled,
                )?;
                let defect = finite(flux - active[i * count + j])?.abs();
                if defect > config.pairing_tolerance_h {
                    return Err(WindingError::Pairing {
                        source: i,
                        receiver: j,
                        flux_h: flux,
                        energy_h: active[i * count + j],
                    });
                }
                maximum_pairing_defect_h = maximum_pairing_defect_h.max(defect);
            }
        }
        let mut inductance_h = active.clone();
        for i in 0..count {
            inductance_h[i * count + i] =
                finite(inductance_h[i * count + i] + end_inductance_h[i])?;
        }
        let (inverse_inductance_per_h, inverse_residual) = inverse(
            &inductance_h,
            count,
            config.inverse_residual_tolerance,
            cancelled,
        )?;
        let mut prescribed_vertices = prescribed_vertices.to_vec();
        prescribed_vertices.sort_unstable();
        poll(cancelled)?;
        Ok(Self {
            problem,
            prescribed_vertices,
            profiles_per_m2,
            fields,
            bias_field,
            bias_flux_linkages_wb,
            maximum_bias_stationarity_defect_wb,
            active_inductance_h: active,
            inductance_h,
            inverse_inductance_per_h,
            end_inductance_h,
            net_turns_by_component,
            maximum_pairing_defect_h,
            inverse_residual,
        })
    }

    /// Exact original spatial material/geometry owner.
    pub fn problem(&self) -> &PlanarMagnetostatics {
        self.problem
    }
    /// Canonical vertices carrying the original homogeneous A_z boundary.
    /// Material replacement retains this physical boundary, not just the mesh.
    pub fn prescribed_vertices(&self) -> &[usize] {
        &self.prescribed_vertices
    }
    /// Number of terminal-current coordinates.
    pub fn winding_count(&self) -> usize {
        self.fields.len()
    }
    /// Signed J_z/terminal-current profiles [1/m2], in winding and cell order.
    pub fn profiles_per_m2(&self) -> &[Vec<f64>] {
        &self.profiles_per_m2
    }
    /// Original incremental unit-current responses and their solver reports.
    /// Remanence is omitted from these derivatives, never from the actual field.
    pub fn unit_fields(&self) -> &[PlanarField<'a>] {
        &self.fields
    }
    /// Source-admitted zero-current field, including fixed-remanence storage.
    pub fn bias_field(&self) -> Option<&PlanarField<'a>> {
        self.bias_field.as_ref()
    }
    /// Constant total-flux offset [Wb], in original winding order.
    pub fn bias_flux_linkages_wb(&self) -> &[f64] {
        &self.bias_flux_linkages_wb
    }
    /// Constant W0 [J] of the zero-current field under its declared recoil
    /// reference. The native circuit evolves energy relative to this state.
    pub fn bias_energy_j(&self) -> f64 {
        self.bias_field.as_ref().map_or(0.0, PlanarField::energy_j)
    }
    /// Measured largest |integral H0.B_unit dV| [Wb]. Exact stationarity is zero.
    pub fn maximum_bias_stationarity_defect_wb(&self) -> f64 {
        self.maximum_bias_stationarity_defect_wb
    }
    /// Symmetric row-major active-section inductance [H].
    pub fn active_inductance_h(&self) -> &[f64] {
        &self.active_inductance_h
    }
    /// Total row-major inductance [H], including explicit end/leakage diagonals.
    pub fn inductance_h(&self) -> &[f64] {
        &self.inductance_h
    }
    /// Symmetric inverse used by the original circuit's Hamiltonian [1/H].
    pub fn inverse_inductance_per_h(&self) -> &[f64] {
        &self.inverse_inductance_per_h
    }
    /// Declared series end/leakage self-inductance [H], not a fitted field result.
    pub fn end_inductance_h(&self) -> &[f64] {
        &self.end_inductance_h
    }
    /// Measured signed closure residual [turns] for each winding and connected
    /// component. Components are ordered by their first mesh vertex.
    pub fn net_turns_by_component(&self) -> &[Vec<f64>] {
        &self.net_turns_by_component
    }
    /// Largest original source-pairing versus energy-Gram discrepancy [H].
    pub fn maximum_pairing_defect_h(&self) -> f64 {
        self.maximum_pairing_defect_h
    }
    /// Measured infinity-norm residual of the symmetric inverse.
    pub fn inverse_residual(&self) -> f64 {
        self.inverse_residual
    }

    /// Total physical flux linkage lambda=lambda0+L*i [Wb], including declared
    /// end storage. Native circuit coordinates are the INCREMENT L*i instead.
    pub fn flux_linkages_wb(&self, currents_a: &[f64]) -> Result<Vec<f64>, WindingError> {
        let mut flux = self.incremental_flux_linkages_wb(currents_a)?;
        for (value, &bias) in flux.iter_mut().zip(&self.bias_flux_linkages_wb) {
            *value = finite(*value + bias)?;
        }
        Ok(flux)
    }
    /// Incremental linkage delta_lambda=L*i [Wb], the native circuit coordinate.
    pub fn incremental_flux_linkages_wb(
        &self,
        currents_a: &[f64],
    ) -> Result<Vec<f64>, WindingError> {
        multiply(&self.inductance_h, currents_a, self.winding_count())
    }
    /// Terminal current i=L^-1*(lambda-lambda0) [A] from TOTAL physical linkage.
    pub fn currents_a(&self, flux_linkages_wb: &[f64]) -> Result<Vec<f64>, WindingError> {
        if flux_linkages_wb.len() != self.winding_count() {
            return Err(WindingError::Invalid(
                "one total flux linkage per winding is required",
            ));
        }
        let incremental = flux_linkages_wb
            .iter()
            .zip(&self.bias_flux_linkages_wb)
            .map(|(&total, &bias)| finite(total - bias))
            .collect::<Result<Vec<_>, _>>()?;
        self.currents_from_incremental_flux_wb(&incremental)
    }
    /// Recover currents from the native incremental magnetic coordinates,
    /// without adding and subtracting a potentially much larger permanent flux.
    pub fn currents_from_incremental_flux_wb(
        &self,
        incremental_flux_wb: &[f64],
    ) -> Result<Vec<f64>, WindingError> {
        multiply(
            &self.inverse_inductance_per_h,
            incremental_flux_wb,
            self.winding_count(),
        )
    }
    /// Evaluate the ORIGINAL distributed source pairing on any finite nodal A.
    /// This exposes constant-gauge and independently solved-field cross-checks.
    /// It includes only the active section; declared end storage is separate.
    pub fn pair_flux_linkages_wb(
        &self,
        potential_wb_m: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<f64>, WindingError> {
        self.profiles_per_m2
            .iter()
            .map(|profile| pairing(self.problem, potential_wb_m, profile, cancelled))
            .collect()
    }
    /// Reconstruct the same fixed affine field from terminal currents, then
    /// transfer it under the existing explicit triangle/position/frame contract.
    /// End/leakage energy has no field in this cross-section. This is still
    /// point sampling, with no inferred conservative voxel projection.
    pub fn sample_for_currents(
        &self,
        currents_a: &[f64],
        frame: &PlanarFrame,
        samples: &[CellSample],
        position_tolerance_m: f64,
        max_samples: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<[f64; 3]>, WindingError> {
        if self.bias_field.is_some() {
            return Ok(self
                .field_for_currents(currents_a, cancelled)?
                .sample_cells(frame, samples, position_tolerance_m, max_samples, cancelled)?);
        }
        // Preserve the sparse point-query cost of the original unbiased lane.
        // No recoil envelope exists there that could require an unseen-cell scan.
        self.admit_currents(currents_a, cancelled)?;
        let mut result = Vec::new();
        for (index, field) in self.fields.iter().enumerate() {
            poll(cancelled)?;
            let contribution =
                field.sample_cells(frame, samples, position_tolerance_m, max_samples, cancelled)?;
            if index == 0 {
                result = vec![[0.0; 3]; contribution.len()];
            }
            for (out, value) in result.iter_mut().zip(contribution) {
                for axis in 0..3 {
                    out[axis] = finite(currents_a[index].mul_add(value[axis], out[axis]))?;
                }
            }
        }
        poll(cancelled)?;
        Ok(result)
    }

    /// Reconstruct A=A0+sum(i_w*A_unit,w) on the original mesh and invoke the
    /// original physical recovery. EVERY magnet cell must satisfy its source
    /// envelope, including cells absent from a subsequent sampling request.
    /// The result has no new solver report; original bias/response reports stay
    /// attached here. Declared end inductance has no cross-section field.
    /// Original static current snapshots are superposed with the potentials,
    /// so the recovered actual field can use the existing static-force owner.
    pub fn field_for_currents(
        &self,
        currents_a: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PlanarField<'a>, WindingError> {
        poll(cancelled)?;
        if currents_a.len() != self.winding_count() || currents_a.iter().any(|v| !v.is_finite()) {
            return Err(WindingError::Invalid(
                "one finite terminal current per winding is required",
            ));
        }
        let mut potential = self.bias_field.as_ref().map_or_else(
            || vec![0.0; self.problem.mesh().vertices().len()],
            |bias| bias.potential_wb_m().to_vec(),
        );
        let mut current_density = vec![0.0; self.problem.mesh().faces().len()];
        if let Some(bias) = &self.bias_field {
            current_density.copy_from_slice(bias.impressed_current_density_a_m2().ok_or(
                WindingError::Invalid(
                    "winding reconstruction requires the original static bias current snapshot",
                ),
            )?);
        }
        for (index, field) in self.fields.iter().enumerate() {
            poll(cancelled)?;
            for (out, &value) in potential.iter_mut().zip(field.potential_wb_m()) {
                poll(cancelled)?;
                *out = finite(currents_a[index].mul_add(value, *out))?;
            }
            let source = field
                .impressed_current_density_a_m2()
                .ok_or(WindingError::Invalid(
                    "winding reconstruction requires original static unit-current snapshots",
                ))?;
            for (out, &value) in current_density.iter_mut().zip(source) {
                poll(cancelled)?;
                *out = finite(currents_a[index].mul_add(value, *out))?;
            }
        }
        Ok(self
            .problem
            .recover_static(potential, current_density, cancelled)?)
    }

    /// Admit actual currents against the complete fixed material field. Circuit
    /// and Ohmic candidate owners share this validator with field sampling.
    pub fn admit_currents(
        &self,
        currents_a: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), WindingError> {
        poll(cancelled)?;
        if currents_a.len() != self.winding_count() || currents_a.iter().any(|v| !v.is_finite()) {
            return Err(WindingError::Invalid(
                "one finite terminal current per winding is required",
            ));
        }
        if self.bias_field.is_none() {
            return Ok(());
        }
        self.field_for_currents(currents_a, cancelled)?;
        Ok(())
    }
}

fn multiply(matrix: &[f64], vector: &[f64], n: usize) -> Result<Vec<f64>, WindingError> {
    if vector.len() != n || vector.iter().any(|v| !v.is_finite()) {
        return Err(WindingError::Invalid(
            "complete finite winding coordinates are required",
        ));
    }
    let mut result = vec![0.0; n];
    for i in 0..n {
        for j in 0..n {
            result[i] = finite(matrix[i * n + j].mul_add(vector[j], result[i]))?;
        }
    }
    Ok(result)
}
fn pairing(
    problem: &PlanarMagnetostatics,
    potential: &[f64],
    profile: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, WindingError> {
    pairing_on_mesh(problem.mesh(), potential, profile, cancelled)
}

fn pairing_on_mesh(
    mesh: &fs_rep_mesh::TriComplex2,
    potential: &[f64],
    profile: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, WindingError> {
    if potential.len() != mesh.vertices().len()
        || profile.len() != mesh.faces().len()
        || profile.iter().any(|v| !v.is_finite())
        || potential.iter().any(|v| !v.is_finite())
    {
        return Err(WindingError::Invalid(
            "flux pairing needs a complete finite nodal potential",
        ));
    }
    let mut flux = 0.0;
    for (cell, face) in mesh.faces().iter().enumerate() {
        poll(cancelled)?;
        // P1 exact cell average; divide before summing to limit intermediate scale.
        let mean = finite(face.iter().map(|&v| potential[v as usize] / 3.0).sum())?;
        let load = finite(profile[cell] * mesh.face_measure(cell).unwrap())?;
        flux = finite(load.mul_add(mean, flux))?;
    }
    Ok(flux)
}
fn inverse(
    matrix: &[f64],
    n: usize,
    tolerance: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<f64>, f64), WindingError> {
    poll(cancelled)?;
    let factor = cholesky(matrix, n).map_err(|_| WindingError::SingularInductance)?;
    let mut inverse = vec![0.0; n * n];
    for column in 0..n {
        poll(cancelled)?;
        let mut rhs = vec![0.0; n];
        rhs[column] = 1.0;
        factor.solve(&mut rhs);
        for row in 0..n {
            inverse[row * n + column] = finite(rhs[row])?;
        }
    }
    // Only floating-point inverse roundoff is averaged. The physical L was
    // already an exactly symmetric energy Gram matrix, checked against J.A.
    for i in 0..n {
        for j in 0..i {
            let value = 0.5 * inverse[i * n + j] + 0.5 * inverse[j * n + i];
            inverse[i * n + j] = value;
            inverse[j * n + i] = value;
        }
    }
    cholesky(&inverse, n).map_err(|_| WindingError::SingularInductance)?;
    let mut residual = 0.0_f64;
    for i in 0..n {
        poll(cancelled)?;
        let mut row = 0.0;
        for j in 0..n {
            let mut entry = if i == j { -1.0 } else { 0.0 };
            for k in 0..n {
                entry = finite(matrix[i * n + k].mul_add(inverse[k * n + j], entry))?;
            }
            row = finite(row + entry.abs())?;
        }
        residual = residual.max(row);
    }
    if residual > tolerance {
        return Err(WindingError::SingularInductance);
    }
    Ok((inverse, residual))
}

fn admit_profiles(
    mesh: &fs_rep_mesh::TriComplex2,
    profiles: &[Vec<f64>],
    thickness_m: f64,
    current_balance_tolerance_turns: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<Vec<f64>>, WindingError> {
    let components = component_labels(mesh, cancelled)?;
    let component_count = components.iter().copied().max().unwrap_or(0) + 1;
    let mut result = Vec::with_capacity(profiles.len());
    for profile in profiles {
        poll(cancelled)?;
        if profile.len() != mesh.faces().len() || profile.iter().any(|v| !v.is_finite()) {
            return Err(WindingError::Invalid("one finite current density per triangle is required"));
        }
        let mut net = vec![0.0; component_count];
        let mut positive = vec![0.0; component_count];
        let mut negative = vec![0.0; component_count];
        for (cell, &density) in profile.iter().enumerate() {
            poll(cancelled)?;
            let component = components[mesh.faces()[cell][0] as usize];
            let turns = finite(density * (mesh.face_measure(cell).unwrap() / thickness_m))?;
            net[component] = finite(net[component] + turns)?;
            if turns > 0.0 { positive[component] = finite(positive[component] + turns)?; }
            if turns < 0.0 { negative[component] = finite(negative[component] - turns)?; }
        }
        if positive.iter().all(|&v| v == 0.0) {
            return Err(WindingError::Invalid("winding needs nonzero outgoing and return current"));
        }
        for component in 0..component_count {
            if (positive[component] > 0.0) != (negative[component] > 0.0)
                || net[component].abs() > current_balance_tolerance_turns
            {
                return Err(WindingError::Invalid("unit current must close separately on every connected mesh component"));
            }
        }
        result.push(net);
    }
    Ok(result)
}

fn component_labels(
    mesh: &fs_rep_mesh::TriComplex2,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<usize>, WindingError> {
    let n = mesh.vertices().len();
    let mut neighbors = vec![Vec::new(); n];
    for &[a, b] in mesh.edges() {
        poll(cancelled)?;
        neighbors[a as usize].push(b as usize);
        neighbors[b as usize].push(a as usize);
    }
    let mut labels = vec![usize::MAX; n];
    let mut stack = Vec::new();
    let mut component = 0;
    for start in 0..n {
        if labels[start] != usize::MAX {
            continue;
        }
        labels[start] = component;
        stack.push(start);
        while let Some(v) = stack.pop() {
            poll(cancelled)?;
            for &w in &neighbors[v] {
                if labels[w] == usize::MAX {
                    labels[w] = component;
                    stack.push(w);
                }
            }
        }
        component += 1;
    }
    Ok(labels)
}
