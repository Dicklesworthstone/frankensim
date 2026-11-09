//! Work-conjugate binding of the original layer thermal law to a facet frame.

use super::super::ShearShellSection;
use super::PlateError;
use fs_material::laminate::thermal::nodal::PreparedNodalLaminateThermal;
use fs_material::laminate::thermal::{
    LaminateThermalError, LaminateThermalResponse, LayerThermalStress, PreparedLaminateThermal,
};
use fs_math::det;
use std::sync::Arc;

/// Thermal mismatch of the original bonded layers, expressed in a facet frame.
///
/// Layer free strains are integrated by the material owner before binding.
/// `section_angle_rad` rotates the entire laminate from its section x/y axes
/// into facet ex/ey, exactly as `ShearShellSection::from_laminate`. The original
/// section reference surface must coincide with the mesh surface. Preparation
/// caches the strain transform and thermal force/moment; clones share all
/// source layers and thermal samples through `Arc`.
#[derive(Debug, Clone)]
pub struct LaminateThermalField {
    source: Arc<PreparedLaminateThermal>,
    section: ShearShellSection,
    section_angle_rad: f64,
    strain_to_section: [[f64; 3]; 3],
    membrane_force_n_m: [f64; 3],
    bending_moment_n: [f64; 3],
}

impl LaminateThermalField {
    /// Bind one fixed layer state and its declared section-to-facet angle.
    pub fn try_new(
        source: Arc<PreparedLaminateThermal>,
        section_angle_rad: f64,
    ) -> Result<Self, PlateError> {
        let section = ShearShellSection::from_laminate(source.section(), section_angle_rad)?;
        let c = det::cos(section_angle_rad);
        let s = det::sin(section_angle_rad);
        let strain_to_section = [
            [c * c, s * s, c * s],
            [s * s, c * c, -c * s],
            [-2.0 * c * s, 2.0 * c * s, c * c - s * s],
        ];
        let membrane_force_n_m = to_facet(&strain_to_section, source.membrane_force_n_m());
        let bending_moment_n = to_facet(&strain_to_section, source.bending_moment_n());
        if membrane_force_n_m
            .iter()
            .chain(&bending_moment_n)
            .any(|x| !x.is_finite())
        {
            return Err(PlateError::BadSection {
                what: "oriented laminate thermal load overflow",
            });
        }
        Ok(Self {
            source,
            section,
            section_angle_rad,
            strain_to_section,
            membrane_force_n_m,
            bending_moment_n,
        })
    }

    /// Exact original laminate, source layers and integrated thermal profiles.
    pub fn source(&self) -> &PreparedLaminateThermal {
        &self.source
    }

    /// Exact oriented section required by this field's native shell assembly.
    /// Reuse it when rebuilding a temperature-dependent section to avoid a
    /// second rotation and keep stiffness, inertia and thermal admission aligned.
    pub const fn section(&self) -> &ShearShellSection {
        &self.section
    }

    /// Whole-section angle relative to the facet's ex/ey axes [rad].
    pub const fn section_angle_rad(&self) -> f64 {
        self.section_angle_rad
    }

    /// Positive equivalent membrane thermal loading `[Nxx,Nyy,Nxy]` [N/m].
    pub const fn membrane_force_n_m(&self) -> [f64; 3] {
        self.membrane_force_n_m
    }

    /// Positive equivalent thermal first moment `[Mxx,Myy,Mxy]` [N].
    pub const fn bending_moment_n(&self) -> [f64; 3] {
        self.bending_moment_n
    }

    /// Conservative whole-facet/layer bound on represented elastic strain,
    /// including transverse shear but not an inferred normal thickness strain.
    /// Supply a field from the SAME prepared thermal shell and nodal snapshot.
    pub fn elastic_strain_bound(
        &self,
        layer: usize,
        field: &super::super::strain::RealElementStrain,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, PlateError> {
        let faces = self
            .source
            .layer_free_strain_in_section(layer)
            .map_err(material_error)?;
        let integral = &self.source.section().layer_integrals()[layer];
        let t = &self.strain_to_section;
        // Inverse engineering-strain rotation; the work-conjugate stress
        // transpose is deliberately not used to transform a free strain.
        let free = faces.map(|e| {
            [
                t[0][0] * e[0] + t[0][1] * e[1] - t[0][2] * e[2],
                t[1][0] * e[0] + t[1][1] * e[1] - t[1][2] * e[2],
                -t[2][0] * e[0] - t[2][1] * e[1] + t[2][2] * e[2],
            ]
        });
        field.reduced_elastic_layer_bound([integral.z_bottom_m, integral.z_top_m], free, cancelled)
    }

    /// Actual layer-face plane stress at supplied compatible facet strains.
    /// Both returned tensors are in facet ex/ey, with undoubled xy shear [Pa].
    /// This is a material evaluation; it does not claim nodal equilibrium or
    /// interlaminar stress recovery.
    pub fn layer_stress(
        &self,
        layer: usize,
        membrane: [f64; 3],
        gradient_per_m: [f64; 3],
    ) -> Result<LayerThermalStress, PlateError> {
        let stress = self
            .source
            .layer_stress(
                layer,
                self.to_section(membrane),
                self.to_section(gradient_per_m),
            )
            .map_err(material_error)?;
        Ok(LayerThermalStress {
            bottom_pa: to_facet(&self.strain_to_section, stress.bottom_pa),
            top_pa: to_facet(&self.strain_to_section, stress.top_pa),
        })
    }

    pub(super) fn admit_section(&self, actual: &ShearShellSection) -> Result<(), PlateError> {
        let expected = &self.section;
        if actual.coupled_matrices() != expected.coupled_matrices()
            || actual.areal_mass_kg_m2() != expected.areal_mass_kg_m2()
            || actual.mass_first_moment_kg_m() != expected.mass_first_moment_kg_m()
            || actual.rotary_inertia_kg() != expected.rotary_inertia_kg()
        {
            return Err(PlateError::BadSection {
                what: "thermal laminate must match original shell section, reference and orientation",
            });
        }
        Ok(())
    }

    pub(in crate::shell::mitc) fn evaluate(
        &self,
        membrane: [f64; 3],
        gradient_per_m: [f64; 3],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<LaminateThermalResponse, PlateError> {
        let response = self
            .source
            .evaluate(
                self.to_section(membrane),
                self.to_section(gradient_per_m),
                cancelled,
            )
            .map_err(material_error)?;
        Ok(LaminateThermalResponse {
            membrane_force_n_m: to_facet(&self.strain_to_section, response.membrane_force_n_m),
            bending_moment_n: to_facet(&self.strain_to_section, response.bending_moment_n),
            strain_energy_j_m2: response.strain_energy_j_m2,
        })
    }

    pub(super) fn to_section(&self, strain: [f64; 3]) -> [f64; 3] {
        std::array::from_fn(|i| {
            (0..3)
                .map(|j| self.strain_to_section[i][j] * strain[j])
                .sum()
        })
    }
}

/// Nodal layer expansion in one fixed laminate and facet material frame.
///
/// Each original triangle node retains its own integrated layer profiles.
/// The existing section transform acts after the material owner interpolates
/// free strain; no section stiffness, source or layer ordering is mixed.
#[derive(Debug, Clone)]
pub struct NodalLaminateThermalField {
    source: Arc<PreparedNodalLaminateThermal>,
    frame: LaminateThermalField,
}

impl NodalLaminateThermalField {
    pub fn try_new(
        source: Arc<PreparedNodalLaminateThermal>,
        section_angle_rad: f64,
    ) -> Result<Self, PlateError> {
        let frame = LaminateThermalField::try_new(source.nodes()[0].clone(), section_angle_rad)?;
        Ok(Self { source, frame })
    }

    pub fn source(&self) -> &PreparedNodalLaminateThermal {
        &self.source
    }

    pub const fn section(&self) -> &ShearShellSection {
        self.frame.section()
    }

    pub const fn section_angle_rad(&self) -> f64 {
        self.frame.section_angle_rad()
    }

    /// Positive equivalent membrane load and bending moment at this point.
    pub fn thermal_resultants_at(
        &self,
        barycentric: [f64; 3],
    ) -> Result<([f64; 3], [f64; 3]), PlateError> {
        let (force, moment) = self
            .source
            .thermal_resultants_at(barycentric)
            .map_err(material_error)?;
        let force = to_facet(&self.frame.strain_to_section, force);
        let moment = to_facet(&self.frame.strain_to_section, moment);
        if force.iter().chain(&moment).any(|value| !value.is_finite()) {
            return Err(PlateError::BadSection {
                what: "oriented nodal laminate thermal load overflow",
            });
        }
        Ok((force, moment))
    }

    /// Actual layer-face stress at a supplied compatible spatial point.
    /// Returned shear is physical tau_xy, not doubled engineering shear.
    pub fn layer_stress(
        &self,
        layer: usize,
        barycentric: [f64; 3],
        membrane: [f64; 3],
        gradient_per_m: [f64; 3],
    ) -> Result<LayerThermalStress, PlateError> {
        let stress = self
            .source
            .layer_stress(
                layer,
                barycentric,
                self.frame.to_section(membrane),
                self.frame.to_section(gradient_per_m),
            )
            .map_err(material_error)?;
        Ok(LayerThermalStress {
            bottom_pa: to_facet(&self.frame.strain_to_section, stress.bottom_pa),
            top_pa: to_facet(&self.frame.strain_to_section, stress.top_pa),
        })
    }

    pub(super) fn admit_section(&self, actual: &ShearShellSection) -> Result<(), PlateError> {
        self.frame.admit_section(actual)
    }

    pub(in crate::shell::mitc) fn evaluate(
        &self,
        barycentric: [f64; 3],
        membrane: [f64; 3],
        gradient_per_m: [f64; 3],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<LaminateThermalResponse, PlateError> {
        let response = self
            .source
            .evaluate(
                barycentric,
                self.frame.to_section(membrane),
                self.frame.to_section(gradient_per_m),
                cancelled,
            )
            .map_err(material_error)?;
        Ok(LaminateThermalResponse {
            membrane_force_n_m: to_facet(
                &self.frame.strain_to_section,
                response.membrane_force_n_m,
            ),
            bending_moment_n: to_facet(&self.frame.strain_to_section, response.bending_moment_n),
            strain_energy_j_m2: response.strain_energy_j_m2,
        })
    }
}

// The transpose of the engineering-strain transform maps work-conjugate
// resultants and actual stress components. It is not the strain transform.
fn to_facet(transform: &[[f64; 3]; 3], stress: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| (0..3).map(|j| transform[j][i] * stress[j]).sum())
}

fn material_error(error: LaminateThermalError) -> PlateError {
    match error {
        LaminateThermalError::Cancelled => PlateError::Cancelled,
        _ => PlateError::BadSection {
            what: "laminate thermal constitutive response refused",
        },
    }
}
