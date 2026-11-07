//! Dispersive forced motion of the original CST/DKT shell.
//!
//! Each retained homogeneous or laminate tensor law is reduced by the material owner's
//! COMPLEX plane-stress constraint at the requested frequency. The original
//! CST membrane and DKT curvature operators integrate its complete anisotropic
//! storage and loss. The shared harmonic owner solves the full physical pencil
//! K'(omega) - omega^2 M - i K''(omega), including optional physical B/Z loads.
//! No mode fit, scalar loss factor, or hereditary time-step approximation enters.

use super::{
    ShellViscosity, ShellViscosityBudget, assemble_constitutive_operator,
    assemble_signed_constitutive_operator, bad, compensated_products, projection_terms,
    strain_rows,
};
use crate::shell::{
    CstDktMass, DRILLING_ALPHA, ReferenceDensityMoments, ShellMesh, ShellModel, ShellSection,
    ShellSupport, add_drilling_stiffness, assemble_shell_materials_with_mass,
    assemble_shell_materials_with_reference_density_moments,
};
use crate::{PlateError, PlateSection};
use fs_exec::Cx;
use fs_material::plane_stress_viscosity::MaterialPlane;
use fs_material::sheet::SheetMatrix;
use fs_material::visco::ViscoError;
use fs_material::visco::tensor::{PlaneStressHarmonic, TensorProny};
use fs_math::{c64::C64, dd::Dd};
use fs_modal::harmonic::{
    HarmonicParameterControls, HarmonicParameterResponse, HarmonicPortControls,
    HarmonicPortParameterDirection, HarmonicResponse, HarmonicResponseControls,
    HarmonicResponseError, PreparedHarmonicResponse,
};
use fs_modal::quadratic_modes::QuadraticParameterDirection;
use fs_sparse::{Coo, Csr};

type Resultants = ([f64; 9], [f64; 9], [f64; 9]);

mod section;
use section::PreparedSection;
pub use section::TensorPronyShellSection;

/// The original material, spatial or shared harmonic refusal.
#[derive(Debug)]
pub enum TensorShellHarmonicError {
    /// Native geometry, reference section, work budget or cancellation.
    Plate(PlateError),
    /// Original tensor frequency, passivity, strain domain or representability.
    Material { element: usize, source: ViscoError },
    /// Original layer, orientation, section integration or source strain domain.
    Laminate {
        element: usize,
        source: fs_material::laminate::LaminateError,
    },
    /// Shared dense factor, physical port, tangent or numerical audit refusal.
    Harmonic(HarmonicResponseError),
}
impl core::fmt::Display for TensorShellHarmonicError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-PLATE-TENSOR-HARMONIC: {self:?}")
    }
}
impl std::error::Error for TensorShellHarmonicError {}
impl From<PlateError> for TensorShellHarmonicError {
    fn from(value: PlateError) -> Self {
        Self::Plate(value)
    }
}
impl From<HarmonicResponseError> for TensorShellHarmonicError {
    fn from(value: HarmonicResponseError) -> Self {
        Self::Harmonic(value)
    }
}
fn checkpoint(cx: &Cx<'_>) -> Result<(), TensorShellHarmonicError> {
    cx.checkpoint().map_err(|_| PlateError::Cancelled.into())
}
fn invalid(what: &'static str) -> TensorShellHarmonicError {
    HarmonicResponseError::Invalid(what).into()
}

/// Native shell and the exact material/geometry selections that assembled it.
/// The model's original reference drilling stiffness and all chosen mass terms
/// remain fixed. Only physical membrane/bending storage disperses; stabilization
/// never acquires material heat. This front door establishes source identity
/// by assembling through the existing shell owner, not by accepting an unrelated
/// public `ShellModel` with matching dimensions.
pub struct TensorPronyShell {
    mesh: ShellMesh,
    sections: Vec<TensorPronyShellSection>,
    drilling: Csr,
    model: ShellModel,
    budget: ShellViscosityBudget,
    reference_omega_rad_s: f64,
}

impl TensorPronyShell {
    /// Assemble the original CST/DKT pencil at an explicitly admitted reference
    /// frequency. Supply one original homogeneous or laminate section per facet.
    /// Native supports retain their full coupling. Laminates require consistent
    /// inertia, including their signed first and physical second density moments.
    /// Material loss is absent from the reference pencil and added exactly once
    /// by each prepared harmonic operator. No modal projection is performed.
    ///
    /// The spatial budget bounds elements, free/full coordinates and cold
    /// assembly work before native matrix allocation. `max_modes` is unused.
    /// Material branch storage belongs to the supplied retained source laws.
    #[allow(clippy::too_many_arguments)] // original source, geometry, inertia and explicit budgets
    pub fn assemble(
        mesh: ShellMesh,
        sections: Vec<TensorPronyShellSection>,
        boundary_nodes: &[usize],
        support: ShellSupport,
        mass: CstDktMass,
        reference_omega_rad_s: f64,
        budget: ShellViscosityBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, TensorShellHarmonicError> {
        checkpoint(cx)?;
        admit_construction(&mesh, &sections, boundary_nodes, support, budget, cx)?;
        let layered = sections.iter().any(|s| s.laminate().is_some());
        if layered && mass != CstDktMass::Consistent {
            return Err(invalid(
                "layered harmonic shells require consistent reference-density inertia",
            ));
        }
        let mut native = Vec::with_capacity(sections.len());
        let mut storage_sections = Vec::with_capacity(sections.len());
        for (element, section) in sections.iter().enumerate() {
            checkpoint(cx)?;
            let (_, (a, b, d), _) = section.prepare(element, reference_omega_rad_s, cx)?;
            let membrane = SheetMatrix::from_coefficients(a)
                .map_err(|_| bad("tensor shell reference membrane is unrepresentable"))?;
            native.push(ShellSection::new(
                PlateSection {
                    d,
                    thickness: section.thickness_m(),
                    density: section.density_kg_m3(),
                },
                membrane,
            )?);
            storage_sections.push(ShellViscosity::new_coupled(a, b, d)?);
        }
        // The original source assembly owns inertia, physical axial rotations,
        // homogeneous support elimination and relative-spin drilling stiffness.
        let mut model = if layered {
            let moments: Result<Vec<_>, PlateError> = sections
                .iter()
                .map(|s| s.density_moments().map(|m| [m; 3]))
                .collect();
            assemble_shell_materials_with_reference_density_moments(
                &mesh,
                &native,
                &moments?,
                boundary_nodes,
                support,
            )?
        } else {
            assemble_shell_materials_with_mass(&mesh, &native, boundary_nodes, support, mass)?
        };
        let drilling = assemble_drilling(&mesh, &native, &model, cx)?;
        if layered {
            // Basic native sections have no B block. Use the SAME original
            // coupled strain owner for reference K and all prepared frequencies;
            // never subtract large reference terms to recover physical storage.
            let physical =
                assemble_constitutive_operator(&mesh, &storage_sections, &model, budget, || {
                    cx.checkpoint().is_err()
                })?;
            model.k = add_operators(&physical, &drilling, model.free, cx)?;
        }
        checkpoint(cx)?;
        Ok(Self {
            mesh,
            sections,
            drilling,
            model,
            budget,
            reference_omega_rad_s,
        })
    }

    /// Original physical geometry. No later load or frequency mutates it.
    pub const fn mesh(&self) -> &ShellMesh {
        &self.mesh
    }
    /// Original supported reference M/K and full-to-free coordinate map.
    pub const fn model(&self) -> &ShellModel {
        &self.model
    }
    /// Immutable material selection per original facet.
    pub fn sections(&self) -> &[TensorPronyShellSection] {
        &self.sections
    }
    /// Reference frequency used only for the fixed native pencil [rad/s].
    pub const fn reference_omega_rad_s(&self) -> f64 {
        self.reference_omega_rad_s
    }

    /// Prepare complete physical forced response at a positive admitted frequency.
    /// The dense controls include the shared factor, three dense source matrices,
    /// and retained material reductions / resultant arrays. Native sparse/COO
    /// storage retains the explicit spatial assembly budget above. A different
    /// frequency needs a new material evaluation and factorization.
    pub fn prepare_harmonic_response(
        &self,
        omega_rad_s: f64,
        controls: HarmonicResponseControls,
        cx: &Cx<'_>,
    ) -> Result<PreparedTensorPronyShell<'_>, TensorShellHarmonicError> {
        self.prepare(omega_rad_s, controls, None, cx)
    }

    /// Attach a full physical boundary impedance in the SAME free coordinate
    /// order. The real row-major B maps conjugate port efforts to forces/torques;
    /// flow is B^T v. Mutual and reactive Z entries participate in the same
    /// factorization and can regularize a bare resonance. Boundary absorbed work
    /// remains separate from intrinsic material heat. The caller owns Z's source
    /// frequency domain and causality; pointwise passivity is checked centrally.
    #[allow(clippy::too_many_arguments)] // physical load and independent source/load budgets
    pub fn prepare_harmonic_response_with_port_impedance(
        &self,
        omega_rad_s: f64,
        force_map: &[f64],
        ports: usize,
        impedance: &[C64],
        controls: HarmonicResponseControls,
        port_controls: HarmonicPortControls,
        cx: &Cx<'_>,
    ) -> Result<PreparedTensorPronyShell<'_>, TensorShellHarmonicError> {
        self.prepare(
            omega_rad_s,
            controls,
            Some((force_map, ports, impedance, port_controls)),
            cx,
        )
    }

    fn prepare(
        &self,
        omega: f64,
        controls: HarmonicResponseControls,
        boundary: Option<(&[f64], usize, &[C64], HarmonicPortControls)>,
        cx: &Cx<'_>,
    ) -> Result<PreparedTensorPronyShell<'_>, TensorShellHarmonicError> {
        checkpoint(cx)?;
        if !omega.is_finite() || omega <= 0.0 {
            return Err(invalid(
                "tensor shell harmonic response requires finite positive frequency",
            ));
        }
        let n = self.model.free;
        if n > controls.maximum_coordinates {
            return Err(HarmonicResponseError::Budget {
                quantity: "shell harmonic coordinates",
                requested: n,
                cap: controls.maximum_coordinates,
            }
            .into());
        }
        if !controls.maximum_relative_residual.is_finite()
            || controls.maximum_relative_residual <= 0.0
            || controls.maximum_relative_residual >= 1.0
        {
            return Err(invalid("shell harmonic residual budget must lie in (0,1)"));
        }
        let section_bytes = self
            .sections
            .iter()
            .try_fold(0usize, |sum, section| {
                PreparedSection::payload_bytes(section.layer_count())
                    .and_then(|v| v.checked_add(2 * core::mem::size_of::<ShellViscosity>()))
                    .and_then(|v| sum.checked_add(v))
            })
            .ok_or_else(|| invalid("shell harmonic layer payload overflows"))?;
        let source_bytes = n
            .checked_mul(n)
            .and_then(|v| v.checked_mul(24))
            .and_then(|v| v.checked_add(section_bytes))
            .ok_or_else(|| invalid("shell harmonic source payload overflows"))?;
        let total = n
            .checked_mul(n)
            .and_then(|v| v.checked_mul(144))
            .and_then(|v| n.checked_mul(128).and_then(|w| v.checked_add(w)))
            .and_then(|v| v.checked_add(source_bytes))
            .ok_or_else(|| invalid("shell harmonic total payload overflows"))?;
        if total > controls.maximum_payload_bytes {
            return Err(HarmonicResponseError::Budget {
                quantity: "shell harmonic numerical payload bytes",
                requested: total,
                cap: controls.maximum_payload_bytes,
            }
            .into());
        }
        let mut materials = Vec::with_capacity(self.sections.len());
        let mut storage_sections = Vec::with_capacity(self.sections.len());
        let mut losses = Vec::with_capacity(self.sections.len());
        for (element, section) in self.sections.iter().enumerate() {
            checkpoint(cx)?;
            let (material, (a, b, d), mut viscosity) = section.prepare(element, omega, cx)?;
            storage_sections.push(ShellViscosity::new_coupled(a, b, d)?);
            for value in viscosity
                .0
                .iter_mut()
                .chain(&mut viscosity.1)
                .chain(&mut viscosity.2)
            {
                let loss = *value;
                *value /= omega;
                if !value.is_finite() || loss != 0.0 && *value == 0.0 {
                    return Err(invalid(
                        "effective shell material viscosity is not representable",
                    ));
                }
            }
            let (a, b, d) = viscosity;
            losses.push(ShellViscosity::new_coupled(a, b, d)?);
            materials.push(material);
        }
        let mass = dense(&self.model.m, n, cx)?;
        let physical_storage = assemble_constitutive_operator(
            &self.mesh,
            &storage_sections,
            &self.model,
            self.budget,
            || cx.checkpoint().is_err(),
        )?;
        let mut storage = dense(&physical_storage, n, cx)?;
        add_sparse(&mut storage, &self.drilling, n, cx)?;
        let loss =
            assemble_constitutive_operator(&self.mesh, &losses, &self.model, self.budget, || {
                cx.checkpoint().is_err()
            })?;
        let damping = dense(&loss, n, cx)?;
        let solve_controls = HarmonicResponseControls {
            maximum_payload_bytes: controls.maximum_payload_bytes - source_bytes,
            ..controls
        };
        let response = if let Some((b, ports, z, load_controls)) = boundary {
            PreparedHarmonicResponse::try_with_port_impedance(
                &mass,
                &damping,
                &storage,
                n,
                b,
                ports,
                z,
                omega,
                solve_controls,
                load_controls,
                cx,
            )?
        } else {
            PreparedHarmonicResponse::try_new(
                &mass,
                &damping,
                &storage,
                n,
                omega,
                solve_controls,
                cx,
            )?
        };
        checkpoint(cx)?;
        Ok(PreparedTensorPronyShell {
            source: self,
            omega_rad_s: omega,
            materials,
            response,
        })
    }
}

/// A cached full-coordinate factor and full complex strain recovery per facet.
/// Repeated forcing cases reuse both; no material inverse is formed during solve.
pub struct PreparedTensorPronyShell<'a> {
    source: &'a TensorPronyShell,
    omega_rad_s: f64,
    materials: Vec<PreparedSection>,
    response: PreparedHarmonicResponse,
}
impl<'a> PreparedTensorPronyShell<'a> {
    /// Unchanged source shell, including original physical free-coordinate map.
    pub const fn source(&self) -> &'a TensorPronyShell {
        self.source
    }
    /// Prepared angular frequency [rad/s].
    pub const fn omega_rad_s(&self) -> f64 {
        self.omega_rad_s
    }

    /// Solve peak force/torque phasors under exp(-i omega t), in original free
    /// displacement/physical-rotation order. `average_viscous_power_w` in the
    /// shared result is exactly intrinsic tensor-Prony plus Kelvin heat; port
    /// absorbed work is separate. No numerical drilling loss is introduced.
    ///
    /// Before publication, recover all six complex Mandel strains and admit the
    /// maximum norm over the ENTIRE harmonic cycle at every original layer face.
    /// DKT curvature is affine over each facet: its vertex extrema bound this
    /// finite-element strain field by convexity, including membrane/bending
    /// interference. This is not a spatial-discretization or shell-model bound.
    pub fn solve(
        &self,
        force: &[C64],
        cx: &Cx<'_>,
    ) -> Result<HarmonicResponse, TensorShellHarmonicError> {
        let result = self.response.solve(force, cx)?;
        self.admit_displacement(&result.displacement, cx)?;
        Ok(result)
    }

    fn admit_displacement(
        &self,
        displacement: &[C64],
        cx: &Cx<'_>,
    ) -> Result<(), TensorShellHarmonicError> {
        for (element, (tri, material)) in self
            .source
            .mesh
            .tris
            .iter()
            .zip(&self.materials)
            .enumerate()
        {
            checkpoint(cx)?;
            let g = self.source.mesh.facet(element)?;
            let local: [C64; 18] = core::array::from_fn(|i| {
                let base = 6 * tri[i / 6] + 3 * ((i % 6) / 3);
                let terms: [(f64, C64); 3] = core::array::from_fn(|p| {
                    (
                        g.frame[i % 3][p],
                        self.source.model.dof_map[base + p].map_or(C64::ZERO, |r| displacement[r]),
                    )
                });
                C64::new(
                    compensated_products(terms.iter().map(|(a, b)| (*a, b.re))),
                    compensated_products(terms.iter().map(|(a, b)| (*a, b.im))),
                )
            });
            let rows = strain_rows(&g, element)?;
            let midpoint: [[C64; 6]; 3] = rows.map(|point| {
                core::array::from_fn(|r| {
                    C64::new(
                        compensated_products((0..18).map(|i| (point[r][i], local[i].re))),
                        compensated_products((0..18).map(|i| (point[r][i], local[i].im))),
                    )
                })
            });
            // Original midpoints are 01, 12, 20; affine extrapolation gives all
            // three vertices without introducing another curvature operator.
            for vertex in 0..3 {
                checkpoint(cx)?;
                let strain: [C64; 6] = core::array::from_fn(|r| {
                    midpoint[vertex][r] + midpoint[(vertex + 2) % 3][r]
                        - midpoint[(vertex + 1) % 3][r]
                });
                material.admit(&self.source.sections[element], element, &strain, cx)?;
            }
        }
        checkpoint(cx)
    }

    /// Cache the exact material-frequency direction beside the same factor.
    /// The added payload limit includes three n*n source arrays plus two signed
    /// resultant arrays. Native sparse temporaries retain the spatial budget.
    /// Reference inertia, drilling and supports have exactly zero derivative.
    pub fn with_frequency_derivative(
        self,
        maximum_payload_bytes: usize,
        cx: &Cx<'_>,
    ) -> Result<PreparedTensorPronyShellFrequency<'a>, TensorShellHarmonicError> {
        checkpoint(cx)?;
        let n = self.source.model.free;
        let layer_bytes = self
            .source
            .sections
            .iter()
            .try_fold(0usize, |sum, section| {
                section
                    .layer_count()
                    .checked_mul(2 * core::mem::size_of::<[f64; 9]>())
                    .and_then(|v| v.checked_add(2 * core::mem::size_of::<Resultants>()))
                    .and_then(|v| sum.checked_add(v))
            })
            .ok_or_else(|| invalid("shell material frequency layer payload overflows"))?;
        let bytes = n
            .checked_mul(n)
            .and_then(|v| v.checked_mul(24))
            .and_then(|v| v.checked_add(layer_bytes))
            .ok_or_else(|| invalid("shell material frequency payload overflows"))?;
        if bytes > maximum_payload_bytes {
            return Err(HarmonicResponseError::Budget {
                quantity: "shell material frequency payload bytes",
                requested: bytes,
                cap: maximum_payload_bytes,
            }
            .into());
        }
        let mut storage = Vec::with_capacity(self.source.sections.len());
        let mut viscosity = Vec::with_capacity(self.source.sections.len());
        for (element, section) in self.source.sections.iter().enumerate() {
            checkpoint(cx)?;
            let (re, slope) = section.frequency_derivative(
                &self.materials[element],
                element,
                self.omega_rad_s,
                cx,
            )?;
            storage.push(re);
            viscosity.push(slope);
        }
        let model = &self.source.model;
        let k = assemble_signed_constitutive_operator(
            &self.source.mesh,
            &storage,
            model,
            self.source.budget,
            || cx.checkpoint().is_err(),
        )?;
        let c = assemble_signed_constitutive_operator(
            &self.source.mesh,
            &viscosity,
            model,
            self.source.budget,
            || cx.checkpoint().is_err(),
        )?;
        let stiffness_derivative = dense(&k, n, cx)?;
        let damping_derivative = dense(&c, n, cx)?;
        let mass_derivative = vec![0.0; n * n];
        checkpoint(cx)?;
        Ok(PreparedTensorPronyShellFrequency {
            base: self,
            mass_derivative,
            damping_derivative,
            stiffness_derivative,
        })
    }
}

/// One factor and its full material-frequency direction, reusable across loads.
pub struct PreparedTensorPronyShellFrequency<'a> {
    base: PreparedTensorPronyShell<'a>,
    mass_derivative: Vec<f64>,
    damping_derivative: Vec<f64>,
    stiffness_derivative: Vec<f64>,
}
impl<'a> PreparedTensorPronyShellFrequency<'a> {
    /// Original prepared response remains available without differentiation.
    pub const fn base(&self) -> &PreparedTensorPronyShell<'a> {
        &self.base
    }
    /// Differentiate actual dispersive motion, heat and boundary work per rad/s.
    /// Authored force and full B/Z slopes are optional; None holds them fixed.
    /// Explicit frequency differentiation still includes inertia, velocity and
    /// load kinematics. No neighboring-frequency solve or eigensystem is used.
    /// Primal full-cycle strain admission applies before returning either result.
    pub fn solve(
        &self,
        force: &[C64],
        force_derivative: Option<&[C64]>,
        port_direction: Option<HarmonicPortParameterDirection<'_>>,
        controls: HarmonicParameterControls,
        cx: &Cx<'_>,
    ) -> Result<HarmonicParameterResponse, TensorShellHarmonicError> {
        let result = self.base.response.solve_with_total_parameter_derivative(
            force,
            QuadraticParameterDirection {
                mass: &self.mass_derivative,
                damping: &self.damping_derivative,
                stiffness: &self.stiffness_derivative,
            },
            force_derivative,
            port_direction,
            1.0,
            controls,
            cx,
        )?;
        self.base
            .admit_displacement(&result.response.displacement, cx)?;
        Ok(result)
    }
}

fn resultants(q: [f64; 9], h: f64) -> Result<Resultants, PlateError> {
    let a = q.map(|v| v * h);
    let d = a.map(|v| (v * (h / 12.0)) * h);
    if a.iter().chain(&d).any(|v| !v.is_finite()) {
        return Err(bad("tensor shell resultants overflow"));
    }
    Ok((a, [0.0; 9], d))
}
fn dense(matrix: &Csr, n: usize, cx: &Cx<'_>) -> Result<Vec<f64>, TensorShellHarmonicError> {
    let mut out = vec![0.0; n * n];
    add_sparse(&mut out, matrix, n, cx)?;
    Ok(out)
}
fn add_sparse(
    out: &mut [f64],
    matrix: &Csr,
    n: usize,
    cx: &Cx<'_>,
) -> Result<(), TensorShellHarmonicError> {
    for i in 0..n {
        checkpoint(cx)?;
        let (columns, values) = matrix.row(i);
        for (&j, &value) in columns.iter().zip(values) {
            out[n * i + j] += value;
            if !out[n * i + j].is_finite() {
                return Err(invalid("shell spectral matrix overflow"));
            }
        }
    }
    Ok(())
}

fn add_operators(
    physical: &Csr,
    drilling: &Csr,
    n: usize,
    cx: &Cx<'_>,
) -> Result<Csr, TensorShellHarmonicError> {
    let mut out = Coo::new(n, n);
    for i in 0..n {
        checkpoint(cx)?;
        for matrix in [physical, drilling] {
            let (columns, values) = matrix.row(i);
            for (&j, &value) in columns.iter().zip(values) {
                if value != 0.0 {
                    out.push(i, j, value);
                }
            }
        }
    }
    let out = out.assemble();
    if (0..n).any(|r| out.row(r).1.iter().any(|v| !v.is_finite())) {
        return Err(bad("layered shell reference stiffness overflow").into());
    }
    checkpoint(cx)?;
    Ok(out)
}

fn assemble_drilling(
    mesh: &ShellMesh,
    sections: &[ShellSection],
    model: &ShellModel,
    cx: &Cx<'_>,
) -> Result<Csr, TensorShellHarmonicError> {
    let mut out = Coo::new(model.free, model.free);
    for (element, (tri, section)) in mesh.tris.iter().zip(sections).enumerate() {
        checkpoint(cx)?;
        let g = mesh.facet(element)?;
        let mut local = [0.0; 324];
        add_drilling_stiffness(
            &g,
            section.membrane().coefficients()[8],
            DRILLING_ALPHA,
            &mut local,
        );
        for i in 0..18 {
            let Some(row) = model.dof_map[6 * tri[i / 6] + i % 6] else {
                continue;
            };
            let bi = 6 * (i / 6) + 3 * ((i % 6) / 3);
            for j in i..18 {
                let Some(column) = model.dof_map[6 * tri[j / 6] + j % 6] else {
                    continue;
                };
                let bj = 6 * (j / 6) + 3 * ((j % 6) / 3);
                let mut value = 0.0;
                for p in 0..3 {
                    for q in 0..3 {
                        value +=
                            g.frame[p][i % 3] * local[(bi + p) * 18 + bj + q] * g.frame[q][j % 3];
                    }
                }
                if !value.is_finite() {
                    return Err(bad("shell drilling projection overflow").into());
                }
                if value != 0.0 {
                    out.push(row, column, value);
                    if i != j {
                        out.push(column, row, value);
                    }
                }
            }
        }
    }
    let out = out.assemble();
    if (0..model.free).any(|r| out.row(r).1.iter().any(|v| !v.is_finite())) {
        return Err(bad("shell drilling assembly overflow").into());
    }
    checkpoint(cx)?;
    Ok(out)
}

fn admit_construction(
    mesh: &ShellMesh,
    sections: &[TensorPronyShellSection],
    boundary: &[usize],
    support: ShellSupport,
    budget: ShellViscosityBudget,
    cx: &Cx<'_>,
) -> Result<(), TensorShellHarmonicError> {
    if sections.len() != mesh.tris.len() {
        return Err(PlateError::SectionCount {
            expected: mesh.tris.len(),
            actual: sections.len(),
        }
        .into());
    }
    if mesh.tris.len() > budget.max_elements {
        return Err(PlateError::Budget {
            what: "tensor shell elements",
            requested: mesh.tris.len(),
            limit: budget.max_elements,
        }
        .into());
    }
    let full = mesh
        .nodes
        .len()
        .checked_mul(6)
        .ok_or_else(|| bad("shell coordinate overflow"))?;
    let mut material_work = 0usize;
    for section in sections {
        checkpoint(cx)?;
        material_work = section
            .material_work()
            .and_then(|v| material_work.checked_add(v))
            .ok_or_else(|| bad("shell material work overflow"))?;
    }
    let work = projection_terms(full, mesh.tris.len(), 18)
        .and_then(|v| v.checked_add(material_work))
        .ok_or_else(|| bad("shell assembly work overflow"))?;
    for (what, requested, limit) in [
        ("tensor shell full coordinates", full, budget.max_full_dofs),
        ("tensor shell assembly work", work, budget.max_work),
    ] {
        if requested > limit {
            return Err(PlateError::Budget {
                what,
                requested,
                limit,
            }
            .into());
        }
    }
    let mut fixed = vec![false; mesh.nodes.len()];
    for &node in boundary {
        checkpoint(cx)?;
        let n = fixed.len();
        *fixed.get_mut(node).ok_or(PlateError::BadBoundary {
            node,
            node_count: n,
        })? = true;
    }
    let per_node = match support {
        ShellSupport::Free => 0,
        ShellSupport::Pinned => 3,
        ShellSupport::Clamped => 6,
    };
    let free = full - per_node * fixed.iter().filter(|v| **v).count();
    if free > budget.max_dofs {
        return Err(PlateError::Budget {
            what: "tensor shell free coordinates",
            requested: free,
            limit: budget.max_dofs,
        }
        .into());
    }
    if free == 0 {
        return Err(bad("tensor shell requires free physical coordinates").into());
    }
    checkpoint(cx)
}
