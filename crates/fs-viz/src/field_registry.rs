//! Versioned field registry for scientific visualization and field export.
//!
//! Bead: `frankensim-extreal-program-f85xj.6.8`
//!
//! Binds requested field names to typed semantic descriptors, component counts,
//! physical units, mesh element association (PointData vs CellData), admissible
//! container formats (VTU, XDMF), and maturity no-claim contracts.

use super::vtu::DataAssociation;

/// Semantic kind of a field.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FieldKind {
    /// Temperature field (scalar).
    Temperature,
    /// Heat flux vector (3 components).
    HeatFlux,
    /// Thermal conductivity (scalar).
    ThermalConductivity,
    /// Fluid pressure (scalar).
    Pressure,
    /// Fluid velocity (vector).
    Velocity,
    /// Geometric region index (discrete integer scalar).
    RegionIndex,
    /// Heat transfer coefficient (scalar).
    HeatTransferCoefficient,
    /// Volumetric heat generation source (scalar).
    VolumetricHeatSource,
}

/// Output format container.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExportFormat {
    /// VTK Unstructured Grid XML format (.vtu).
    Vtu,
    /// eXtensible Data Model and Format (.xdmf).
    Xdmf,
}

/// Specification descriptor for a registrable simulation field.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FieldDescriptor {
    /// Unique semantic identifier.
    pub semantic_id: &'static str,
    /// Canonical display / export name.
    pub canonical_name: &'static str,
    /// Physical field category.
    pub kind: FieldKind,
    /// Association with mesh points or cells.
    pub association: DataAssociation,
    /// Number of scalar components (1 for scalar, 3 for 3D vector).
    pub components: usize,
    /// Default standard unit.
    pub default_unit: &'static str,
    /// Set of allowed export container formats.
    pub allowed_formats: &'static [ExportFormat],
    /// Development maturity level.
    pub maturity_level: &'static str,
    /// Explicit scope and no-claim boundary.
    pub no_claim_boundary: &'static str,
}

/// The authoritative registry of Cooling and Flow solution fields.
pub const FIELD_REGISTRY: &[FieldDescriptor] = &[
    FieldDescriptor {
        semantic_id: "field.thermal.temperature",
        canonical_name: "Temperature",
        kind: FieldKind::Temperature,
        association: DataAssociation::PointData,
        components: 1,
        default_unit: "K",
        allowed_formats: &[ExportFormat::Vtu, ExportFormat::Xdmf],
        maturity_level: "L2",
        no_claim_boundary: "Scalar nodal field derived from steady-state conduction PCG solve; does not attest transient or radiation effects",
    },
    FieldDescriptor {
        semantic_id: "field.thermal.heat_flux",
        canonical_name: "HeatFlux",
        kind: FieldKind::HeatFlux,
        association: DataAssociation::CellData,
        components: 3,
        default_unit: "W/m^2",
        allowed_formats: &[ExportFormat::Vtu, ExportFormat::Xdmf],
        maturity_level: "L2",
        no_claim_boundary: "3D vector element flux derived from -k grad(T); piecewise-constant across tetrahedra",
    },
    FieldDescriptor {
        semantic_id: "field.material.conductivity",
        canonical_name: "ThermalConductivity",
        kind: FieldKind::ThermalConductivity,
        association: DataAssociation::CellData,
        components: 1,
        default_unit: "W/(m·K)",
        allowed_formats: &[ExportFormat::Vtu, ExportFormat::Xdmf],
        maturity_level: "L2",
        no_claim_boundary: "Isotropic conductivity assigned from admitted material card packs",
    },
    FieldDescriptor {
        semantic_id: "field.geometry.region_id",
        canonical_name: "RegionId",
        kind: FieldKind::RegionIndex,
        association: DataAssociation::CellData,
        components: 1,
        default_unit: "1",
        allowed_formats: &[ExportFormat::Vtu, ExportFormat::Xdmf],
        maturity_level: "L2",
        no_claim_boundary: "Categorical integer tag distinguishing die, TIM, spreader, heatsink, and chassis domains",
    },
    FieldDescriptor {
        semantic_id: "field.boundary.htc",
        canonical_name: "HeatTransferCoefficient",
        kind: FieldKind::HeatTransferCoefficient,
        association: DataAssociation::PointData,
        components: 1,
        default_unit: "W/(m^2·K)",
        allowed_formats: &[ExportFormat::Vtu, ExportFormat::Xdmf],
        maturity_level: "L2",
        no_claim_boundary: "Surface Robin conductance derived from flow-network channel velocities; ideal correlation, not local 3D boundary layer",
    },
    FieldDescriptor {
        semantic_id: "field.source.volumetric_heat",
        canonical_name: "VolumetricHeatSource",
        kind: FieldKind::VolumetricHeatSource,
        association: DataAssociation::CellData,
        components: 1,
        default_unit: "W/m^3",
        allowed_formats: &[ExportFormat::Vtu, ExportFormat::Xdmf],
        maturity_level: "L2",
        no_claim_boundary: "Prescribed active die heat generation distributed uniformly over target tetrahedral volume",
    },
];

/// Look up a field descriptor by its canonical name (case-insensitive).
#[must_use]
pub fn find_field_by_name(name: &str) -> Option<&'static FieldDescriptor> {
    FIELD_REGISTRY.iter().find(|f| {
        f.canonical_name.eq_ignore_ascii_case(name) || f.semantic_id.eq_ignore_ascii_case(name)
    })
}

/// Validate whether an output request matches the field registry.
pub fn validate_output_request(
    name: &str,
    format: ExportFormat,
) -> Result<&'static FieldDescriptor, String> {
    let desc = find_field_by_name(name)
        .ok_or_else(|| format!("unknown output field `{name}`; see FIELD_REGISTRY"))?;
    if !desc.allowed_formats.contains(&format) {
        return Err(format!(
            "field `{}` does not support requested format `{:?}`",
            desc.canonical_name, format
        ));
    }
    Ok(desc)
}
