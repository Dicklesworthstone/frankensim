//! Physical in-plane vibration about a native finite-strain equilibrium.
//!
//! The original plane-stress or mean-dilatation owner solves the loaded state
//! and supplies its complete condensed Hessian, including initial-stress
//! effects. Native P1/Q1 inertia integrates the authored reference density and
//! thickness. Neither the relaxed/current thickness nor a prescribed natural
//! deformation changes that original mass inventory.
//!
//! The existing fs-modal dense owner returns the complete signed spectrum and
//! residual intervals on at most 128 free coordinates. Negative and null
//! eigenvalues are retained; neither an absolute value nor a frequency clamp
//! turns an unstable or unresolved mode into a stable vibration.
//!
//! Transverse relaxation remains quasistatic, as in the original plane owner:
//! no transverse inertia, bending, shell or wrinkling model is added. These are
//! small vibrations about one held elastic equilibrium, not evolving plastic
//! or viscoelastic modes. Frozen incremental potentials are explicitly refused.

use super::{HyperProblem, NewtonReport};
use crate::SolidError;
use crate::inertia::{InertiaError, ReferenceMass};
use fs_exec::Cx;
use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_modal::{ModalError, ModePair, eigh_gen_dense};
use fs_sparse::Csr;

pub use crate::hyper3::vibration::{MAX_DENSE_FREE_DOFS, VibrationSettings};

/// The SAME constitutive reduction governs equilibrium and vibration.
#[derive(Clone, Copy, Debug)]
pub enum PlaneVibrationFormulation {
    PlaneStress(PlaneStressSettings),
    IncompressiblePlaneStress(PlaneStressSettings),
    /// Original Q1 mean dilatation, with its existing native-law restrictions.
    MeanDilatation,
}

/// Density in the original reference mesh, retained without averaging.
#[derive(Clone, Copy, Debug)]
pub enum PlaneReferenceDensity<'a> {
    PerElement(&'a [f64]),
    /// Positive original P1/Q1 nodal volumetric density [kg/m³].
    Nodal(&'a [f64]),
}

/// Explicit native mass-assembly and bounded dense-spectrum controls.
#[derive(Clone, Copy, Debug)]
pub struct PlaneVibrationSettings {
    pub spectrum: VibrationSettings,
    pub maximum_mass_assembly_bytes: usize,
}

#[derive(Debug)]
pub enum PlaneVibrationError {
    Solid(SolidError),
    Modal(ModalError),
    UnresolvedMode {
        index: usize,
        eigenvalue_per_s2: f64,
        residual_per_s2: f64,
    },
}

impl From<SolidError> for PlaneVibrationError {
    fn from(error: SolidError) -> Self {
        Self::Solid(error)
    }
}

impl core::fmt::Display for PlaneVibrationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite-strain plane vibration: {self:?}")
    }
}

impl std::error::Error for PlaneVibrationError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Solid(error) => Some(error),
            Self::Modal(error) => Some(error),
            Self::UnresolvedMode { .. } => None,
        }
    }
}

/// Immutable loaded state, original mass and signed tangent spectrum.
///
/// The original mesh, material/reference owners, dead tractions and support
/// positions remain borrowed together. Modes use the retained free-DOF order
/// and zero increments on prescribed coordinates. K and M are row-major.
pub struct PlaneVibration<'p, 'm> {
    problem: &'p HyperProblem<'m>,
    formulation: PlaneVibrationFormulation,
    reference_thickness_m: f64,
    density: PlaneReferenceDensity<'p>,
    reference_mass: ReferenceMass,
    displacement_m: Vec<[f64; 2]>,
    equilibrium: NewtonReport,
    reactions_n: Vec<(usize, f64)>,
    free_dofs: Vec<usize>,
    stiffness_n_m: Vec<f64>,
    mass_kg: Vec<f64>,
    modes: Vec<ModePair>,
}

fn invalid(what: &str) -> SolidError {
    SolidError::InvalidInput { what: what.into() }
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), SolidError> {
    cx.checkpoint().map_err(|_| SolidError::Cancelled)
}

fn inertia_error(error: InertiaError) -> SolidError {
    match error {
        InertiaError::Cancelled => SolidError::Cancelled,
        InertiaError::InvalidInput { what } => invalid(&what),
    }
}

fn restrict(matrix: &Csr, free: &[usize], scale: f64, cx: &Cx<'_>) -> Result<Vec<f64>, SolidError> {
    let n = free.len();
    let mut dense = vec![0.0; n * n];
    for (row, &i) in free.iter().enumerate() {
        checkpoint(cx)?;
        for (col, &j) in free.iter().enumerate() {
            let value = matrix.get(i, j);
            let scaled = value * scale;
            if !scaled.is_finite() || (value != 0.0 && scaled == 0.0) {
                return Err(invalid("unrepresentable physical plane vibration pencil"));
            }
            dense[row * n + col] = scaled;
        }
    }
    let largest = dense.iter().map(|v| v.abs()).fold(0.0, f64::max);
    let tolerance = (64.0 * f64::EPSILON * n as f64) * largest;
    // Admit the physical pencil before the shared modal owner symmetrizes
    // whitening roundoff. Never silently repair a nonsymmetric material.
    for row in 0..n {
        for col in 0..row {
            if (dense[row * n + col] - dense[col * n + row]).abs() > tolerance {
                return Err(invalid(
                    "plane vibration pencil is not symmetric to roundoff",
                ));
            }
        }
    }
    Ok(dense)
}

impl<'m> HyperProblem<'m> {
    /// Solve original loaded equilibrium, then its complete in-plane spectrum.
    ///
    /// Density [kg/m³] and thickness [m] belong to the original reference.
    /// Dead tractions retain the native stress units; thickness multiplies both
    /// their force resultants and the consistent stiffness. Material domains,
    /// exact-volume/thickness relaxation, support lifting and Newton controls
    /// remain with the selected existing owner.
    ///
    /// Only original elastic stored-energy owners are admitted. No supplied
    /// displacement can bypass equilibrium or be reused with foreign sources.
    pub fn solve_vibration<'p>(
        &'p self,
        formulation: PlaneVibrationFormulation,
        reference_thickness_m: f64,
        density: PlaneReferenceDensity<'p>,
        settings: PlaneVibrationSettings,
        cx: &Cx<'_>,
    ) -> Result<PlaneVibration<'p, 'm>, PlaneVibrationError> {
        self.solve_vibration_impl(
            formulation,
            reference_thickness_m,
            density,
            None,
            settings,
            cx,
        )
    }

    /// The same physical solve with an admitted Newton starting field.
    /// The starting field is re-equilibrated under the original absolute
    /// support targets and load ramps; it is not an accepted spectral state.
    pub fn solve_vibration_from<'p>(
        &'p self,
        initial_displacement: &[[f64; 2]],
        formulation: PlaneVibrationFormulation,
        reference_thickness_m: f64,
        density: PlaneReferenceDensity<'p>,
        settings: PlaneVibrationSettings,
        cx: &Cx<'_>,
    ) -> Result<PlaneVibration<'p, 'm>, PlaneVibrationError> {
        self.solve_vibration_impl(
            formulation,
            reference_thickness_m,
            density,
            Some(initial_displacement),
            settings,
            cx,
        )
    }

    fn solve_vibration_impl<'p>(
        &'p self,
        formulation: PlaneVibrationFormulation,
        reference_thickness_m: f64,
        density: PlaneReferenceDensity<'p>,
        initial: Option<&[[f64; 2]]>,
        settings: PlaneVibrationSettings,
        cx: &Cx<'_>,
    ) -> Result<PlaneVibration<'p, 'm>, PlaneVibrationError> {
        checkpoint(cx)?;
        let tolerance = settings.spectrum.relative_residual_tolerance;
        if !tolerance.is_finite() || tolerance <= 0.0 || tolerance >= 1.0 {
            return Err(invalid("modal residual tolerance must be finite and in (0,1)").into());
        }
        self.validate_materials_with(|| checkpoint(cx))?;
        for (element, _) in self.mesh.elems.iter().enumerate() {
            checkpoint(cx)?;
            for point in 0..self.material.material_point_count(element).unwrap_or(1) {
                if self
                    .point_material(element, point)?
                    .elastic_potential()
                    .is_none()
                {
                    return Err(invalid(
                        "physical vibration requires elastic stored energy, not a frozen incremental potential",
                    ).into());
                }
            }
        }
        // Native fixed-DOF selection indexes the original boundary nodes.
        if self
            .mesh
            .patches
            .iter()
            .flat_map(|(_, edges)| edges)
            .any(|&(a, b)| a >= self.mesh.nodes.len() || b >= self.mesh.nodes.len())
        {
            return Err(invalid("plane vibration boundary references an absent mesh node").into());
        }
        let fixed = self.fixed_dofs(1.0)?;
        if fixed.values().any(|value| !value.is_finite()) {
            return Err(invalid("finite prescribed plane displacement required").into());
        }
        let ndofs = self
            .mesh
            .nodes
            .len()
            .checked_mul(2)
            .ok_or_else(|| invalid("plane vibration displacement count overflows"))?;
        let free_dofs: Vec<_> = (0..ndofs).filter(|d| !fixed.contains_key(d)).collect();
        let n = free_dofs.len();
        if n == 0 || n > settings.spectrum.maximum_free_dofs || n > MAX_DENSE_FREE_DOFS {
            return Err(invalid("plane dense vibration requires 1..=128 bounded free DOFs").into());
        }
        // Admit and assemble the original consistent inertia before Newton.
        let reference_mass = match density {
            PlaneReferenceDensity::PerElement(values) => ReferenceMass::from_plane(
                self.mesh,
                reference_thickness_m,
                values,
                settings.maximum_mass_assembly_bytes,
                cx,
            ),
            PlaneReferenceDensity::Nodal(values) => ReferenceMass::from_plane_nodal_density(
                self.mesh,
                reference_thickness_m,
                values,
                settings.maximum_mass_assembly_bytes,
                cx,
            )
            .map(|(mass, _)| mass),
        }
        .map_err(inertia_error)?;
        let mass_kg = restrict(reference_mass.matrix(), &free_dofs, 1.0, cx)?;
        let (displacement_m, equilibrium) = match formulation {
            PlaneVibrationFormulation::PlaneStress(local) => match initial {
                Some(u) => self.solve_plane_stress_from(u, local, cx),
                None => self.solve_plane_stress(local, cx),
            },
            PlaneVibrationFormulation::IncompressiblePlaneStress(local) => match initial {
                Some(u) => self.solve_incompressible_plane_stress_from(u, local, cx),
                None => self.solve_incompressible_plane_stress(local, cx),
            },
            PlaneVibrationFormulation::MeanDilatation => match initial {
                Some(u) => self.solve_mean_dilatation_from(u, cx),
                None => self.solve_mean_dilatation(cx),
            },
        }?;
        let u: Vec<_> = displacement_m.iter().flatten().copied().collect();
        let (residual, tangent) = match formulation {
            PlaneVibrationFormulation::PlaneStress(local) => {
                self.plane_stress_residual_and_tangent(&u, 1.0, local, cx)
            }
            PlaneVibrationFormulation::IncompressiblePlaneStress(local) => {
                self.incompressible_plane_stress_residual_and_tangent(&u, 1.0, local, cx)
            }
            PlaneVibrationFormulation::MeanDilatation => {
                self.mean_dilatation_residual_and_tangent(&u, 1.0, cx)
            }
        }?;
        let stiffness_n_m = restrict(&tangent, &free_dofs, reference_thickness_m, cx)?;
        let reactions_n: Vec<_> = fixed
            .keys()
            .map(|&d| (d, residual[d] * reference_thickness_m))
            .collect();
        if reactions_n.iter().any(|(_, value)| !value.is_finite()) {
            return Err(invalid("unrepresentable plane support reaction").into());
        }
        checkpoint(cx)?;
        let modes =
            eigh_gen_dense(&stiffness_n_m, &mass_kg, n).map_err(PlaneVibrationError::Modal)?;
        checkpoint(cx)?;
        if modes.len() != n {
            return Err(invalid("modal owner returned an incomplete plane basis").into());
        }
        let scale = modes
            .iter()
            .map(|mode| mode.lambda.abs())
            .fold(0.0, f64::max);
        for (index, mode) in modes.iter().enumerate() {
            if mode.phi.len() != n
                || mode.phi.iter().any(|v| !v.is_finite())
                || !mode.lambda.is_finite()
                || !mode.residual.is_finite()
                || mode.residual < 0.0
                || !mode.interval.0.is_finite()
                || !mode.interval.1.is_finite()
                || mode.residual > tolerance * scale
            {
                return Err(PlaneVibrationError::UnresolvedMode {
                    index,
                    eigenvalue_per_s2: mode.lambda,
                    residual_per_s2: mode.residual,
                });
            }
        }
        Ok(PlaneVibration {
            problem: self,
            formulation,
            reference_thickness_m,
            density,
            reference_mass,
            displacement_m,
            equilibrium,
            reactions_n,
            free_dofs,
            stiffness_n_m,
            mass_kg,
            modes,
        })
    }
}

impl<'p, 'm> PlaneVibration<'p, 'm> {
    pub fn problem(&self) -> &HyperProblem<'m> {
        self.problem
    }
    pub const fn formulation(&self) -> PlaneVibrationFormulation {
        self.formulation
    }
    pub const fn reference_thickness_m(&self) -> f64 {
        self.reference_thickness_m
    }
    pub const fn reference_density(&self) -> PlaneReferenceDensity<'p> {
        self.density
    }
    pub fn reference_mass_kg(&self) -> f64 {
        self.reference_mass.reference_mass_kg()
    }
    pub fn displacement_m(&self) -> &[[f64; 2]] {
        &self.displacement_m
    }
    pub fn equilibrium_report(&self) -> &NewtonReport {
        &self.equilibrium
    }
    pub fn reactions_n(&self) -> &[(usize, f64)] {
        &self.reactions_n
    }
    pub fn free_dofs(&self) -> &[usize] {
        &self.free_dofs
    }
    pub fn stiffness_n_m(&self) -> &[f64] {
        &self.stiffness_n_m
    }
    pub fn mass_kg(&self) -> &[f64] {
        &self.mass_kg
    }
    /// Signed eigenvalues [1/s²], original mass-normalized shapes and residual intervals.
    pub fn modes(&self) -> &[ModePair] {
        &self.modes
    }

    /// Actual nodal force from authored acceleration in the original frame.
    /// Applies the SAME full consistent mass, including prescribed coordinates.
    pub fn reference_acceleration_load_n(
        &self,
        acceleration_m_s2: &[[f64; 2]],
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, PlaneVibrationError> {
        self.reference_mass
            .action(acceleration_m_s2.as_flattened(), cx)
            .map_err(inertia_error)
            .map_err(Into::into)
    }
}
