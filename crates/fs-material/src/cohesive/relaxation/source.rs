//! Original ordered-interface spectra; no bulk-to-joint isotropy defaults.
use super::{IntactViscoelasticBond, ViscoelasticBondLimits};
use crate::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedInterfaceStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_interface_state_point,
};
use crate::visco::GeneralizedMaxwell;
use fs_matdb::{InterfaceSystemCard, QueryPoint, SurfaceSpec};
use fs_qty::{Dims, Length, Pressure, QuantitySpec, Time};

mod temperature;

pub const BONDLINE_THICKNESS_PROPERTY: &str = "cohesive_bondline_thickness";
pub const NORMAL_RELAXATION_PREFIX: &str = "cohesive_normal_relaxation";
pub const SHEAR_RELAXATION_PREFIX: &str = "cohesive_shear_relaxation";

fn invalid(reason: &'static str) -> MaterialStatePointError {
    MaterialStatePointError::InvalidRequirement {
        property: "intact viscoelastic bond".into(),
        reason,
    }
}
fn coefficient(
    source: &ResolvedInterfaceStatePoint,
    name: &str,
    dims: Dims,
) -> Result<f64, MaterialStatePointError> {
    let p = source
        .property(name)
        .ok_or_else(|| MaterialStatePointError::InvalidRequirement {
            property: name.into(),
            reason: "missing or ambiguous original interface relaxation coefficient",
        })?;
    if p.requirement().quantity() != QuantitySpec::dimensional(dims)
        || !p.is_static_at_fixed_external_conditions()
    {
        return Err(invalid(
            "relaxation coefficients require fixed external conditions and correct units",
        ));
    }
    let v = p.value_si();
    if !v.is_finite() || v < 0.0 {
        return Err(invalid(
            "nonnegative finite relaxation coefficient required",
        ));
    }
    Ok(v)
}
fn branch_count(
    source: &ResolvedInterfaceStatePoint,
    prefix: &str,
) -> Result<usize, MaterialStatePointError> {
    let n = coefficient(source, &format!("{prefix}_branch_count"), Dims::NONE)?;
    if n.fract() != 0.0 || n > 64.0 {
        return Err(invalid(
            "relaxation branch count must be an integer in 0..64",
        ));
    }
    Ok(n as usize)
}
fn requirement(
    name: impl Into<String>,
    dims: Dims,
    admissibility: ScalarAdmissibility,
) -> Result<ScalarPropertyRequirement, MaterialStatePointError> {
    ScalarPropertyRequirement::try_new(name, dims, admissibility)
}

/// Immutable original source bundle plus the explicit intact bond law.
#[derive(Debug, Clone)]
pub struct ResolvedViscoelasticBond {
    source: ResolvedInterfaceStatePoint,
    surfaces: [SurfaceSpec; 2],
    law: IntactViscoelasticBond,
}
impl ResolvedViscoelasticBond {
    /// Require thickness [m] plus independent normal/shear spectra. For either
    /// prefix: _branch_count [-], _equilibrium_modulus [Pa], and indexed
    /// _modulus_N [Pa]/_time_N [s], N=1..count. The normal modulus is the
    /// measured/declared JOINT normal response under its lateral constraint;
    /// it is not silently replaced by bulk E or bulk K.
    pub fn try_from_resolved(
        card: &InterfaceSystemCard,
        source: ResolvedInterfaceStatePoint,
        limits: ViscoelasticBondLimits,
    ) -> Result<Self, MaterialStatePointError> {
        if source.card_identity() != card.content_hash() {
            return Err(invalid("original ordered interface card required"));
        }
        let spectrum = |prefix: &str| -> Result<GeneralizedMaxwell, MaterialStatePointError> {
            let n = branch_count(&source, prefix)?;
            let equilibrium = coefficient(
                &source,
                &format!("{prefix}_equilibrium_modulus"),
                Pressure::DIMS,
            )?;
            let mut terms = Vec::with_capacity(n);
            for j in 1..=n {
                terms.push((
                    coefficient(&source, &format!("{prefix}_modulus_{j}"), Pressure::DIMS)?,
                    coefficient(&source, &format!("{prefix}_time_{j}"), Time::DIMS)?,
                ));
            }
            GeneralizedMaxwell::new(equilibrium, terms)
                .map_err(|_| invalid("inadmissible original relaxation spectrum"))
        };
        let thickness = coefficient(&source, BONDLINE_THICKNESS_PROPERTY, Length::DIMS)?;
        let law = IntactViscoelasticBond::new(
            spectrum(NORMAL_RELAXATION_PREFIX)?,
            spectrum(SHEAR_RELAXATION_PREFIX)?,
            thickness,
            limits,
        )
        .map_err(|_| invalid("bond thickness, spectra or model-use limits are inadmissible"))?;
        Ok(Self {
            source,
            surfaces: [card.surface_a().clone(), card.surface_b().clone()],
            law,
        })
    }
    pub fn law(&self) -> &IntactViscoelasticBond {
        &self.law
    }
    pub fn resolved(&self) -> &ResolvedInterfaceStatePoint {
        &self.source
    }
    pub fn surfaces(&self) -> &[SurfaceSpec; 2] {
        &self.surfaces
    }
}

/// Select every branch declared by the original source; no truncation or
/// missing-mode fallback. Runtime jump/rate/time limits are caller-authored
/// applicability gates, not new experimental validity inferred from this query.
pub fn resolve_viscoelastic_bond(
    card: &InterfaceSystemCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
    limits: ViscoelasticBondLimits,
) -> Result<ResolvedViscoelasticBond, MaterialStatePointError> {
    let mut requirements = Vec::new();
    for prefix in [NORMAL_RELAXATION_PREFIX, SHEAR_RELAXATION_PREFIX] {
        requirements.push(requirement(
            format!("{prefix}_branch_count"),
            Dims::NONE,
            ScalarAdmissibility::NonNegative,
        )?);
    }
    // The count discovery query owns only its two pins. The final query still
    // validates the caller's COMPLETE plan, including missing or foreign pins.
    let count_selection = match &selection {
        MaterialPropertySelection::PinnedByProperty(pins) => {
            MaterialPropertySelection::PinnedByProperty(
                pins.iter()
                    .filter(|(name, _)| requirements.iter().any(|r| r.name() == name))
                    .cloned()
                    .collect(),
            )
        }
        other => other.clone(),
    };
    let counts = resolve_interface_state_point(card, point, &requirements, count_selection)?;
    requirements.push(requirement(
        BONDLINE_THICKNESS_PROPERTY,
        Length::DIMS,
        ScalarAdmissibility::StrictlyPositive,
    )?);
    for prefix in [NORMAL_RELAXATION_PREFIX, SHEAR_RELAXATION_PREFIX] {
        requirements.push(requirement(
            format!("{prefix}_equilibrium_modulus"),
            Pressure::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?);
        for j in 1..=branch_count(&counts, prefix)? {
            requirements.push(requirement(
                format!("{prefix}_modulus_{j}"),
                Pressure::DIMS,
                ScalarAdmissibility::NonNegative,
            )?);
            requirements.push(requirement(
                format!("{prefix}_time_{j}"),
                Time::DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )?);
        }
    }
    let resolved = resolve_interface_state_point(card, point, &requirements, selection)?;
    ResolvedViscoelasticBond::try_from_resolved(card, resolved, limits)
}
