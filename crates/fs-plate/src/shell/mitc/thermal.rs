//! Prescribed thermal eigenstrain on the original shear-shell operator.
//!
//! Each facet receives an already integrated free strain
//! `epsilon_th(z) = membrane + z * gradient`, in its local tangent frame.
//! The engineering ordering is `[epsilon_xx, epsilon_yy, gamma_xy]`.
//! Each coefficient may be constant or P1-interpolated over the facet. The
//! gradient is the coefficient of material distance along the director, not
//! a DKT slope-curvature convention. The original A/B/D coefficients yield
//! `N_th = A membrane + B gradient`, `M_th = B membrane + D gradient`.
//!
//! The native enriched strain field, elastic bubble condensation, support
//! coordinates and static solver are reused. Both the condensed thermal force
//! and the affine thermal bubble displacement enter recovery. Inertia and the
//! reference geometry remain those of the original immutable model.
//!
//! This is one-way prescribed-temperature, infinitesimal-strain elasticity.
//! Section stiffness stays fixed at the caller's explicitly supplied values.
//! Callers retain source receipts and integrate coefficients before lowering;
//! neither endpoint CTE multiplication nor an isotropic material default is
//! inferred. Layer-dependent free strains use the original laminate owner's
//! thickness integration through [`LaminateThermalField`]. Each layer may have
//! its own affine profile; a nonlinear profile needs thickness refinement.
//! No heat equation, buckling tangent, physical
//! drilling stress, interlaminar stress or thermoelastic heating is implied.

mod laminate;
pub use laminate::LaminateThermalField;

use super::{MitcShellModel, PlateError, element, solve_static};
use element::thermal::{PreparedThermalFacet, ThermalLoading};
use fs_material::state_point::{
    IntegratedIsotropicThermalExpansion, IntegratedOrthotropicThermalExpansion,
};
use fs_math::det;

/// Integrated thermal strain, affine through thickness at a spatial point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ShellThermalStrain {
    membrane: [f64; 3],
    gradient_per_m: [f64; 3],
}

/// Spatial interpolation of an already integrated thermal strain profile.
/// Nodal values are in this facet's SAME ex/ey frame and original node order.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FacetThermalStrain {
    /// Constant over the triangular midsurface.
    Uniform(ShellThermalStrain),
    /// Affine interpolation of three independently integrated node profiles.
    /// This interpolates strain, not CTE or temperature: nonlinear alpha(T)
    /// must be integrated before selecting this explicit spatial approximation.
    AtNodes([ShellThermalStrain; 3]),
}

impl FacetThermalStrain {
    pub(super) fn at(&self, lambda: [f64; 3]) -> ShellThermalStrain {
        match self {
            Self::Uniform(strain) => *strain,
            Self::AtNodes(nodes) => ShellThermalStrain {
                membrane: core::array::from_fn(|i| {
                    (0..3).map(|n| lambda[n] * nodes[n].membrane[i]).sum()
                }),
                gradient_per_m: core::array::from_fn(|i| {
                    (0..3).map(|n| lambda[n] * nodes[n].gradient_per_m[i]).sum()
                }),
            },
        }
    }
}

impl From<ShellThermalStrain> for FacetThermalStrain {
    fn from(strain: ShellThermalStrain) -> Self {
        Self::Uniform(strain)
    }
}

impl ShellThermalStrain {
    /// Admit already integrated engineering free strain and its signed
    /// thickness gradient in the facet's ex/ey axes. Negative expansion and
    /// cooling are supported. These are strains, never CTEs in 1/K.
    pub fn try_new(membrane: [f64; 3], gradient_per_m: [f64; 3]) -> Result<Self, PlateError> {
        if membrane
            .iter()
            .chain(&gradient_per_m)
            .any(|v| !v.is_finite())
        {
            return Err(PlateError::BadSection {
                what: "finite integrated shell thermal strains required",
            });
        }
        Ok(Self {
            membrane,
            gradient_per_m,
        })
    }

    /// Lower the shared material owner's source-integrated small-strain path
    /// to a uniform-thickness plane-stress profile. Retain that original path
    /// for its material identity, temperature domain and source receipts.
    pub fn from_isotropic(path: &IntegratedIsotropicThermalExpansion) -> Result<Self, PlateError> {
        let strain = path.free_linear_strain();
        Self::try_new([strain, strain, 0.0], [0.0; 3])
    }

    /// Lower source-integrated principal strains with axis 3 explicitly normal
    /// to the shell and axis 1 rotated by `angle_rad` from facet ex. The plane-
    /// stress shell uses axes 1/2; it does not constrain thickness strain or
    /// convert normal expansion into an in-plane stress. Retain the source path.
    pub fn from_orthotropic_in_plane(
        path: &IntegratedOrthotropicThermalExpansion,
        angle_rad: f64,
    ) -> Result<Self, PlateError> {
        if !angle_rad.is_finite() {
            return Err(PlateError::BadSection {
                what: "finite in-plane thermal material angle required",
            });
        }
        let [a, b, _] = *path.free_principal_strains();
        let c = det::cos(angle_rad);
        let s = det::sin(angle_rad);
        Self::try_new(
            [
                c * c * a + s * s * b,
                s * s * a + c * c * b,
                2.0 * c * s * (a - b),
            ],
            [0.0; 3],
        )
    }

    /// An explicitly affine free-strain profile between z=-h/2 and z=+h/2.
    /// The mesh reference surface must be the same midplane. Integrate the
    /// actual material law at both faces before calling; a nonlinear profile
    /// is not silently reconstructed from temperatures or endpoint CTEs.
    pub fn from_faces(
        bottom: [f64; 3],
        top: [f64; 3],
        thickness_m: f64,
    ) -> Result<Self, PlateError> {
        if !thickness_m.is_finite() || thickness_m <= 0.0 {
            return Err(PlateError::BadSection {
                what: "positive finite thermal shell thickness required",
            });
        }
        Self::try_new(
            core::array::from_fn(|i| 0.5 * bottom[i] + 0.5 * top[i]),
            core::array::from_fn(|i| (top[i] - bottom[i]) / thickness_m),
        )
    }

    /// Engineering free strain at the mesh reference surface.
    #[must_use]
    pub const fn membrane(&self) -> [f64; 3] {
        self.membrane
    }

    /// Signed engineering strain gradient along the material director [1/m].
    #[must_use]
    pub const fn gradient_per_m(&self) -> [f64; 3] {
        self.gradient_per_m
    }
}

/// Area-averaged physical stress resultants in one original facet frame.
/// Shear resultants use actual stress components, without engineering doubling.
#[derive(Debug, Clone, PartialEq)]
pub struct ThermalFacetResponse {
    /// Membrane `[Nxx, Nyy, Nxy]` [N/m], after subtracting thermal free strain.
    pub membrane_force_n_m: [f64; 3],
    /// Signed first stress moments `[Mxx, Myy, Mxy]` about the mesh surface [N].
    pub bending_moment_n: [f64; 3],
    /// Transverse shear `[Qx, Qy]` [N/m].
    pub transverse_shear_n_m: [f64; 2],
    /// Integrated physical elastic storage [J], with thermal strain removed.
    pub strain_energy_j: f64,
    /// World nodal forces/moments, in original triangle order, excluding the
    /// numerical drilling gauge. They are not an independent applied load.
    pub material_nodal_force: [f64; 18],
}

/// Compatible nodal state and physical material response from the same solve.
#[derive(Debug, Clone, PartialEq)]
pub struct ThermalShellResponse {
    /// World displacement/physical rotation, six entries per original node.
    pub displacement_6n: Vec<f64>,
    /// Original-facet-order stress resultants and physical strain energy.
    pub facets: Vec<ThermalFacetResponse>,
    /// Sum of facet material energies, excluding numerical drilling storage.
    pub strain_energy_j: f64,
    /// Assembled world material forces/moments, excluding numerical drilling.
    /// Support reactions may additionally contain the model's gauge force.
    pub material_nodal_force_6n: Vec<f64>,
}

/// Immutable thermal loading prepared against the exact original shell pencil.
/// Borrowing prevents geometry, supports or material changes after admission.
pub struct PreparedThermalShell<'a> {
    model: &'a MitcShellModel,
    facets: Vec<PreparedThermalFacet>,
    load_6n: Vec<f64>,
}

impl MitcShellModel {
    /// Bind one integrated affine thickness profile per original facet.
    /// Cold preparation is O(facets), polls cancellation per quadrature point,
    /// and never publishes partially assembled loads.
    pub fn prepare_thermal(
        &self,
        strains: &[FacetThermalStrain],
        cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedThermalShell<'_>, PlateError> {
        self.prepare_thermal_fields(
            strains.len(),
            strains.iter().copied().map(ThermalLoading::Strain),
            cancelled,
        )
    }

    /// Bind independently integrated layer expansion to each original facet.
    /// The laminate, reference surface, section angle and all physical section
    /// coefficients must match the section used to assemble this model. A
    /// shared prepared laminate may be reused across facets without copying its
    /// layers. Spatial variation is explicit through the per-facet fields.
    pub fn prepare_laminate_thermal(
        &self,
        fields: &[LaminateThermalField],
        cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedThermalShell<'_>, PlateError> {
        self.prepare_thermal_fields(
            fields.len(),
            fields.iter().cloned().map(ThermalLoading::Laminate),
            cancelled,
        )
    }

    fn prepare_thermal_fields(
        &self,
        count: usize,
        fields: impl Iterator<Item = ThermalLoading>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedThermalShell<'_>, PlateError> {
        if cancelled() {
            return Err(PlateError::Cancelled);
        }
        if count != self.mesh.tris.len() {
            return Err(PlateError::SectionCount {
                expected: self.mesh.tris.len(),
                actual: count,
            });
        }
        let mut facets = Vec::with_capacity(count);
        let mut load_6n = vec![0.0; 6 * self.mesh.nodes.len()];
        for (e, (loading, tri)) in fields.zip(&self.mesh.tris).enumerate() {
            if cancelled() {
                return Err(PlateError::Cancelled);
            }
            if let ThermalLoading::Laminate(laminate) = &loading {
                laminate.admit_section(&self.sections[e])?;
            }
            let field = element::PreparedStrainField::new(
                &self.mesh,
                tri,
                self.mesh.facet(e)?,
                tri.map(|n| self.directors.vectors[n]),
                &self.sections[e],
                &mut cancelled,
            )?;
            let facet =
                PreparedThermalFacet::new(field, &self.sections[e], loading, &mut cancelled)?;
            for i in 0..18 {
                load_6n[6 * tri[i / 6] + i % 6] += facet.load[i];
            }
            facets.push(facet);
        }
        if load_6n.iter().any(|v| !v.is_finite()) {
            return Err(PlateError::BadSection {
                what: "thermal shell load accumulation overflow",
            });
        }
        if cancelled() {
            return Err(PlateError::Cancelled);
        }
        Ok(PreparedThermalShell {
            model: self,
            facets,
            load_6n,
        })
    }
}

impl PreparedThermalShell<'_> {
    /// Exact work-conjugate world force/moment vector, including constrained DOFs.
    /// Existing source displacement maps and modal projections can consume it.
    #[must_use]
    pub fn load_6n(&self) -> &[f64] {
        &self.load_6n
    }

    /// Solve with these thermal forces plus explicit world dead nodal loads.
    /// Uses the original static factorization and true-residual gate. Cancellation
    /// is polled around that existing factorization and throughout recovery;
    /// this adapter does not add interruption inside the native direct solve.
    pub fn solve_static(
        &self,
        external_load_6n: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ThermalShellResponse, PlateError> {
        if cancelled() {
            return Err(PlateError::Cancelled);
        }
        self.model.project_load(external_load_6n)?;
        let load: Vec<_> = self
            .load_6n
            .iter()
            .zip(external_load_6n)
            .map(|(a, b)| a + b)
            .collect();
        let displacement = solve_static(self.model, &load)?;
        self.recover(&displacement, &mut cancelled)
    }

    /// Recover actual material stress resultants and energy at an explicitly
    /// supplied physical displacement snapshot. This performs no solve and
    /// makes no equilibrium claim. It includes thermally driven internal
    /// director motion even when all physical nodal displacements are zero.
    pub fn recover(
        &self,
        displacement_6n: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ThermalShellResponse, PlateError> {
        if cancelled() {
            return Err(PlateError::Cancelled);
        }
        self.model.project_load(displacement_6n)?;
        let mut facets = Vec::with_capacity(self.facets.len());
        let mut material_nodal_force_6n = vec![0.0; displacement_6n.len()];
        let mut strain_energy_j = 0.0;
        for (prepared, tri) in self.facets.iter().zip(&self.model.mesh.tris) {
            let nodal = core::array::from_fn(|i| displacement_6n[6 * tri[i / 6] + i % 6]);
            let facet = prepared.recover(&nodal, &mut cancelled)?;
            strain_energy_j += facet.strain_energy_j;
            for i in 0..18 {
                material_nodal_force_6n[6 * tri[i / 6] + i % 6] += facet.material_nodal_force[i];
            }
            facets.push(facet);
        }
        if !strain_energy_j.is_finite() || material_nodal_force_6n.iter().any(|v| !v.is_finite()) {
            return Err(PlateError::BadSection {
                what: "thermal shell response accumulation overflow",
            });
        }
        if cancelled() {
            return Err(PlateError::Cancelled);
        }
        Ok(ThermalShellResponse {
            displacement_6n: displacement_6n.to_vec(),
            facets,
            strain_energy_j,
            material_nodal_force_6n,
        })
    }
}
