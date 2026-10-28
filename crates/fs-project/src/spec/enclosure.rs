//! Declared closed-enclosure irradiation, separate from ambient reservoirs.

use super::{
    BTreeSet, ConductionSetup, QtyAny, ThermalBoundaryCondition, Violation, check_dims, dims,
    is_canonical_binding_text, is_hex_digest, violation,
};

/// One material-backed patch in an explicitly ordered closed enclosure.
#[derive(Debug, Clone, PartialEq)]
pub struct EnclosureRadiatingSurface {
    /// Unique result identity. Matrix rows and columns follow surface order.
    pub name: String,
    /// Existing complete uniform convective or prescribed-flux boundary target.
    pub target: String,
    /// Immutable surface-finish material-card content hash.
    pub card: String,
    /// Exact emissivity claim pin, when selection would otherwise be ambiguous.
    pub claim: Option<String>,
    /// Explicit temperature for the immutable emissivity query, K.
    pub query_temperature: QtyAny,
}

/// Declared origin of supplied factors, not proof of the imported geometry.
#[derive(Debug, Clone, PartialEq)]
pub enum EnclosureViewFactorEvidence {
    /// Named analytic geometry or formula revision supplied by the author.
    Analytic {
        /// Geometry or closed-form formula identifier and revision.
        geometry: String,
    },
    /// External randomized generation with explicit reproducibility coordinates.
    ExternalQmc {
        /// Generator seed, never inferred from the project's physics seed.
        seed: u64,
        /// Positive sample count representable as an exact wire integer.
        samples: u64,
        /// Generator and version identifier.
        generator: String,
    },
}

/// Closed gray-diffuse surface exchange, with no external radiation reservoir.
///
/// View factors and source evidence are authored inputs. Execution binds every
/// patch to actual boundary faces and checks area-weighted reciprocity there.
/// Row closure alone cannot validate geometric visibility or occlusion.
#[derive(Debug, Clone, PartialEq)]
pub struct ConductionEnclosure {
    /// Two to 64 uniquely named, nonoverlapping boundary targets.
    pub surfaces: Vec<EnclosureRadiatingSurface>,
    /// Complete square matrix in the exact declared surface order.
    pub view_factors: Vec<Vec<f64>>,
    /// Explicit analytic or external-generation origin.
    pub evidence: EnclosureViewFactorEvidence,
    /// Absolute allowed row-sum discrepancy, finite and in [0, 1).
    pub closure_tolerance: f64,
    /// Allowed relative area-weighted reciprocity discrepancy, in [0, 1).
    pub reciprocity_relative_tolerance: f64,
    /// Positive outer-iteration budget.
    pub max_iterations: u32,
    /// Positive surface-temperature convergence tolerance, K.
    pub temperature_tolerance: QtyAny,
    /// Positive patch-heat convergence tolerance, W.
    pub heat_tolerance: QtyAny,
    /// Temperature under-relaxation in (0, 1].
    pub relaxation: f64,
}

pub(super) fn validate(
    setup: &ConductionSetup,
    enclosure: &ConductionEnclosure,
    out: &mut Vec<Violation>,
) {
    if setup.radiation.is_some() {
        out.push(violation(
            "project-enclosure-radiation-overlap",
            "ambient reservoirs and a closed enclosure are both declared",
            "select one radiation model; no reservoir declaration may be silently ignored",
        ));
    }
    let n = enclosure.surfaces.len();
    if !(2..=64).contains(&n) {
        out.push(violation(
            "project-enclosure-radiation-surfaces",
            format!("closed enclosure has {n} surfaces"),
            "declare between two and 64 complete boundary patches",
        ));
        return;
    }
    if enclosure.max_iterations == 0 {
        out.push(violation(
            "project-enclosure-radiation-iterations",
            "closed-enclosure max-iterations is zero",
            "declare a positive outer-iteration budget",
        ));
    }
    if !(enclosure.relaxation.is_finite()
        && enclosure.relaxation > 0.0
        && enclosure.relaxation <= 1.0)
    {
        out.push(violation(
            "project-enclosure-radiation-relaxation",
            "closed-enclosure relaxation is outside (0, 1]",
            "declare finite under-relaxation in (0, 1]",
        ));
    }
    for (name, value, expected) in [
        (
            "temperature-tolerance",
            enclosure.temperature_tolerance,
            dims::TEMPERATURE,
        ),
        ("heat-tolerance", enclosure.heat_tolerance, dims::POWER),
    ] {
        check_dims(
            out,
            "project-enclosure-radiation-dims",
            name,
            value,
            expected,
        );
        if !(value.value.is_finite() && value.value > 0.0) {
            out.push(violation(
                "project-enclosure-radiation-tolerance",
                format!("closed-enclosure {name} is not finite and positive"),
                "declare a positive convergence tolerance with the named physical unit",
            ));
        }
    }
    for (name, tolerance) in [
        ("closure", enclosure.closure_tolerance),
        ("reciprocity", enclosure.reciprocity_relative_tolerance),
    ] {
        if !tolerance.is_finite() || !(0.0..1.0).contains(&tolerance) {
            out.push(violation(
                "project-enclosure-radiation-factor-tolerance",
                format!("{name} tolerance must be finite and in [0, 1)"),
                "declare explicit bounded tolerances; no normalization or factor repair is implied",
            ));
        }
    }
    match &enclosure.evidence {
        EnclosureViewFactorEvidence::Analytic { geometry }
            if !is_canonical_binding_text(geometry) =>
        {
            out.push(violation(
                "project-enclosure-radiation-evidence",
                "analytic view-factor evidence has no canonical geometry/formula identifier",
                "name the analytic formula and its geometry assumptions",
            ));
        }
        EnclosureViewFactorEvidence::ExternalQmc {
            samples, generator, ..
        } if *samples == 0
            || i64::try_from(*samples).is_err()
            || !is_canonical_binding_text(generator) =>
        {
            out.push(violation(
                "project-enclosure-radiation-evidence",
                "external-QMC evidence needs a generator and positive exact sample count",
                "retain the actual seed, sample count fitting i64, and generator/version identifier",
            ));
        }
        _ => {}
    }
    if enclosure.view_factors.len() != n || enclosure.view_factors.iter().any(|row| row.len() != n)
    {
        out.push(violation(
            "project-enclosure-radiation-matrix",
            "view-factor dimensions differ from the ordered surfaces",
            "supply exactly one row and column per surface",
        ));
    }
    for (i, row) in enclosure.view_factors.iter().enumerate() {
        if row
            .iter()
            .any(|value| !value.is_finite() || !(0.0..=1.0).contains(value))
        {
            out.push(violation(
                "project-enclosure-radiation-factor",
                format!("view-factor row {i} contains a value outside [0, 1]"),
                "supply finite physical factors without clipping or implicit normalization",
            ));
        }
        let sum: f64 = row.iter().sum();
        if !sum.is_finite() || (sum - 1.0).abs() > enclosure.closure_tolerance {
            out.push(violation(
                "project-enclosure-radiation-closure",
                format!("view-factor row {i} does not close within its declared tolerance"),
                "close the declared enclosure; an omitted external environment is not an enclosure surface",
            ));
        }
    }
    let mut names = BTreeSet::new();
    let mut targets = BTreeSet::new();
    for surface in &enclosure.surfaces {
        if !is_canonical_binding_text(&surface.name)
            || !names.insert(&surface.name)
            || !is_canonical_binding_text(&surface.target)
            || !targets.insert(&surface.target)
        {
            out.push(violation(
                "project-enclosure-radiation-identity",
                "enclosure surface names and boundary targets must be unique and canonical",
                "name each complete physical patch once; matrix order follows this declaration",
            ));
        }
        if !setup.boundaries.iter().any(|row| {
            row.target == surface.target
                && matches!(
                    row.condition,
                    ThermalBoundaryCondition::Convection { .. }
                        | ThermalBoundaryCondition::AirflowConvection { .. }
                        | ThermalBoundaryCondition::NaturalConvection { .. }
                        | ThermalBoundaryCondition::HeatFlux { .. }
                )
        }) {
            out.push(violation(
                "project-enclosure-radiation-boundary",
                format!("enclosure surface {} lacks a supported original boundary", surface.name),
                "bind each patch to complete uniform convection or prescribed outward heat flux; zero flux declares vacuum, and fixed-temperature patches remain separate",
            ));
        }
        if !is_hex_digest(&surface.card)
            || surface
                .claim
                .as_ref()
                .is_some_and(|claim| !is_hex_digest(claim))
        {
            out.push(violation(
                "project-enclosure-radiation-card",
                format!(
                    "surface {} has an invalid material or claim identity",
                    surface.name
                ),
                "supply the exact 64-hex immutable emissivity card and any explicit claim pin",
            ));
        }
        check_dims(
            out,
            "project-enclosure-radiation-dims",
            "query-temperature",
            surface.query_temperature,
            dims::TEMPERATURE,
        );
        if !surface.query_temperature.value.is_finite() || surface.query_temperature.value <= 0.0 {
            out.push(violation(
                "project-enclosure-radiation-temperature",
                format!("surface {} has an invalid emissivity query temperature", surface.name),
                "declare a positive absolute temperature; execution separately checks card validity",
            ));
        }
    }
}
