//! Physical ply shear on the existing condensed MITC trial space.
//!
//! The material owner prepares the equilibrium profile once per source stack.
//! Each facet caches only its original strain observation and frame rotation.
//! A displacement snapshot is observed once; subsequent layer/point queries
//! allocate nothing and retain the complete native assumed-shear polynomial.
//! The result inherits the source's declared two-probe bending approximation;
//! it is not a free-edge, delamination or transverse-normal-stress solution.

use super::strain::{PreparedElementStrain, RealElementStrain};
use super::{MitcShellModel, PlateError, ShearShellSection};
use fs_material::laminate::{LaminateShearRecovery, LayerShearResponse};
use fs_math::det;

fn invalid(what: &'static str) -> PlateError {
    PlateError::BadSection { what }
}

/// Original facet observation bound to the exact oriented laminate operators.
pub struct PreparedLaminateShearFacet<'a> {
    model: &'a MitcShellModel,
    element: usize,
    source: &'a LaminateShearRecovery,
    observation: PreparedElementStrain,
    sin: f64,
    cos: f64,
}

/// Actual native strain at one instant, with the original layer recovery.
/// It can be queried repeatedly without re-observing the nodal displacement.
pub struct LaminateShearField<'a> {
    source: &'a LaminateShearRecovery,
    native: RealElementStrain,
    sin: f64,
    cos: f64,
}

/// The same local response in original ply/section axes and facet ex/ey axes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FacetLayerShearResponse {
    /// Original material response; its xz/yz components use SECTION axes.
    pub layer: LayerShearResponse,
    /// Physical [tau_ex,n, tau_ey,n] in facet axes [Pa].
    pub facet_stress_pa: [f64; 2],
    /// Physical engineering [gamma_ex,n, gamma_ey,n] in facet axes.
    pub facet_engineering_strain: [f64; 2],
}

impl MitcShellModel {
    /// Prepare local ply shear for one original facet and laminate source.
    /// `section_angle_rad` has exactly the meaning used by
    /// [`ShearShellSection::from_laminate`]. All elastic blocks and inertia
    /// moments must match the retained assembly; a stiffness-only lookalike or
    /// changed reference surface refuses. No solve or motion update occurs.
    pub fn prepare_laminate_shear<'a>(
        &'a self,
        element: usize,
        source: &'a LaminateShearRecovery,
        section_angle_rad: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedLaminateShearFacet<'a>, PlateError> {
        if cancelled() {
            return Err(PlateError::Cancelled);
        }
        let actual = self
            .sections
            .get(element)
            .ok_or(invalid("laminate shear facet out of range"))?;
        admit(actual, source, section_angle_rad)?;
        let observation = PreparedElementStrain::new(
            &self.mesh,
            element,
            &self.directors,
            actual,
            &mut cancelled,
        )?;
        Ok(PreparedLaminateShearFacet {
            model: self,
            element,
            source,
            observation,
            sin: det::sin(section_angle_rad),
            cos: det::cos(section_angle_rad),
        })
    }
}

fn admit(
    actual: &ShearShellSection,
    source: &LaminateShearRecovery,
    angle: f64,
) -> Result<(), PlateError> {
    let expected = ShearShellSection::from_laminate(source.section(), angle)?;
    if actual.coupled_matrices() != expected.coupled_matrices()
        || actual.areal_mass_kg_m2() != expected.areal_mass_kg_m2()
        || actual.mass_first_moment_kg_m() != expected.mass_first_moment_kg_m()
        || actual.rotary_inertia_kg() != expected.rotary_inertia_kg()
    {
        return Err(invalid(
            "shear recovery must match the original laminate, reference and facet orientation",
        ));
    }
    Ok(())
}

impl<'a> PreparedLaminateShearFacet<'a> {
    /// Observe a full physical world-coordinate displacement/rotation snapshot.
    /// Validation reads only this facet's nodes, so E facet observations cost
    /// O(E), not O(E times all nodes). Support coordinates must be expanded by
    /// the original model before this call. Thermal director offsets require
    /// [`super::thermal::PreparedThermalShell::laminate_shear_field`] instead.
    pub fn evaluate(
        &self,
        displacement_6n: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<LaminateShearField<'a>, PlateError> {
        if cancelled() {
            return Err(PlateError::Cancelled);
        }
        if displacement_6n.len() != 6 * self.model.mesh.nodes.len() {
            return Err(invalid(
                "one physical displacement/rotation per original shell DOF required",
            ));
        }
        let tri = self.model.mesh.tris[self.element];
        let nodal = core::array::from_fn(|i| displacement_6n[6 * tri[i / 6] + i % 6]);
        let native = self.observation.real_field(&nodal, &mut cancelled)?;
        Ok(LaminateShearField {
            source: self.source,
            native,
            sin: self.sin,
            cos: self.cos,
        })
    }
}

impl<'a> LaminateShearField<'a> {
    pub(super) fn from_thermal(
        section: &ShearShellSection,
        source: &'a LaminateShearRecovery,
        angle: f64,
        native: RealElementStrain,
    ) -> Result<Self, PlateError> {
        admit(section, source, angle)?;
        Ok(Self {
            source,
            native,
            sin: det::sin(angle),
            cos: det::cos(angle),
        })
    }

    /// Complete original native strain polynomial, including any thermal
    /// director offset supplied by the native thermal owner.
    pub const fn native_strain(&self) -> &RealElementStrain {
        &self.native
    }

    /// Physical shear at an interior barycentric location and a fraction
    /// through an original ply. Stress is continuous across bonded interfaces;
    /// strains follow each ply's original transverse compliance.
    pub fn at_barycentric(
        &self,
        barycentric: [f64; 3],
        layer: usize,
        fraction: f64,
    ) -> Result<FacetLayerShearResponse, PlateError> {
        let gamma = self.native.at_barycentric(barycentric)?.transverse_shear;
        let section_gamma = [
            self.cos * gamma[0] + self.sin * gamma[1],
            -self.sin * gamma[0] + self.cos * gamma[1],
        ];
        let layer = self
            .source
            .from_section_strain(layer, fraction, section_gamma)
            .map_err(|_| invalid("original laminate shear recovery refused"))?;
        let to_facet = |v: [f64; 2]| {
            [
                self.cos * v[0] - self.sin * v[1],
                self.sin * v[0] + self.cos * v[1],
            ]
        };
        let facet_stress_pa = to_facet(layer.stress_pa);
        let facet_engineering_strain = to_facet(layer.engineering_strain);
        if facet_stress_pa
            .iter()
            .chain(&facet_engineering_strain)
            .any(|x| !x.is_finite())
        {
            return Err(invalid("unrepresentable rotated laminate shear response"));
        }
        Ok(FacetLayerShearResponse {
            layer,
            facet_stress_pa,
            facet_engineering_strain,
        })
    }
}
