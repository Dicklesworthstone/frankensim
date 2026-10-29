//! Free geometry and optional constant-mass density from engineering length changes.
//!
//! Source epsilon is `(L(T)-L_datum)/L_datum`, with a fixed datum for each
//! direction. It is neither a thermal expansion coefficient nor an arbitrary
//! infinitesimal strain tensor. Normalize each source at the declared reference:
//! `s_i(T) = (1+epsilon_i(T))/(1+epsilon_i(T_ref))`, then `J = s_1*s_2*s_3`
//! and, when a compatible reference density exists, `rho(T) = rho_ref/J`.
//! Dimensional change itself requires no density. A fitted datum residual is preserved.
//!
//! The caller declares isotropy or three fixed orthogonal principal directions
//! and unconstrained free expansion. Density composition also requires fixed mass.
//! This adapter does not infer
//! matched specimens, rotating directions, shear, phase changes, mass loss,
//! moisture evolution or constrained thermoelastic response. Source uncertainty
//! survives unchanged; no covariance or propagated uncertainty is invented.

use super::{ThermalTransportError, admit_temperature_point, positive, source_line, source_value};
use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::{
    ClaimSelection, MaterialAnswer, MaterialCard, PropertyKey, PropertyLineSlice, QueryPoint,
};
use fs_qty::{Dims, QuantitySpec};

/// One or three canonical source lines, with no adaptive resampling.
pub const MAX_FREE_EXPANSION_KNOTS: usize = 3 * fs_matdb::MAX_PROPERTY_LINE_KNOTS;

/// Explicit geometric assumption; a scalar source does not establish isotropy.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FreeExpansionBasis {
    Isotropic,
    /// Nonzero identity of fixed, mutually orthogonal material directions.
    FixedPrincipal {
        frame: ContentHash,
    },
}

/// Caller-declared engineering length-change sources and independent selections.
#[derive(Clone, Copy, Debug)]
pub enum FreeExpansionStrainSources<'a> {
    Isotropic {
        strain: &'a PropertyKey,
        selection: ClaimSelection,
    },
    FixedPrincipal {
        frame: ContentHash,
        /// Three engineering length changes in the declared frame's axis order.
        strains: [&'a PropertyKey; 3],
        selections: [ClaimSelection; 3],
    },
}

/// Exact density reference and continuous, fixed-context free-expansion path.
#[derive(Clone, Copy, Debug)]
pub struct FreeExpansionSourceQuery<'a> {
    pub density: &'a PropertyKey,
    pub density_selection: ClaimSelection,
    pub reference: &'a QueryPoint,
    pub strains: FreeExpansionStrainSources<'a>,
    pub temperature_axis: &'a str,
    pub lower: &'a QueryPoint,
    pub upper: &'a QueryPoint,
}

/// Continuous free dimensional change with an explicit reference and source basis.
/// Density is not an input to length ratios or their volume ratio.
#[derive(Clone, Copy, Debug)]
pub struct FreeExpansionGeometryQuery<'a> {
    pub reference: &'a QueryPoint,
    pub strains: FreeExpansionStrainSources<'a>,
    pub temperature_axis: &'a str,
    pub lower: &'a QueryPoint,
    pub upper: &'a QueryPoint,
}

/// Derived free geometry, not an independently measured property claim.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreeExpansionGeometrySample {
    path_identity: ContentHash,
    temperature_k: f64,
    principal_stretches: [f64; 3],
    relative_volume: f64,
}

impl FreeExpansionGeometrySample {
    #[must_use]
    pub const fn path_identity(&self) -> ContentHash {
        self.path_identity
    }
    #[must_use]
    pub const fn temperature_k(&self) -> f64 {
        self.temperature_k
    }
    /// Length ratios relative to the declared reference, in principal-axis order.
    #[must_use]
    pub const fn principal_stretches(&self) -> [f64; 3] {
        self.principal_stretches
    }
    /// Finite volume ratio `V(T)/V(T_ref)` of the declared free expansion.
    #[must_use]
    pub const fn relative_volume(&self) -> f64 {
        self.relative_volume
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct StrainKnot {
    temperature: f64,
    strains: [f64; 3],
}

/// Immutable source-backed reduction. Only engineering strains interpolate;
/// stretches and their product are derived at evaluation time. Missing density
/// does not block the geometry supported by these original strain measurements.
#[derive(Clone, Debug, PartialEq)]
pub struct PreparedFreeExpansion {
    material_card: ContentHash,
    basis: FreeExpansionBasis,
    temperature_axis: String,
    lower: QueryPoint,
    upper: QueryPoint,
    reference: QueryPoint,
    sources: Vec<PropertyLineSlice>,
    reference_lengths: [f64; 3],
    knots: Vec<StrainKnot>,
    identity: ContentHash,
}

impl PreparedFreeExpansion {
    /// Admit each original strain source through canonical typed selection over
    /// the full interval. The reference shares every fixed condition and lies
    /// inside that interval. Unsupported spans and changing selections refuse.
    /// Each source explicitly bounds temperature. No datum value is guessed.
    pub fn try_new(
        card: &MaterialCard,
        query: FreeExpansionGeometryQuery<'_>,
    ) -> Result<Self, ThermalTransportError> {
        let invalid = ThermalTransportError::Invalid;
        let (basis, selected) = match query.strains {
            FreeExpansionStrainSources::Isotropic { strain, selection } => {
                (FreeExpansionBasis::Isotropic, vec![(strain, selection)])
            }
            FreeExpansionStrainSources::FixedPrincipal {
                frame,
                strains,
                selections,
            } => {
                if frame == ContentHash([0; 32]) {
                    return Err(invalid(
                        "fixed principal directions require a nonzero frame identity",
                    ));
                }
                (
                    FreeExpansionBasis::FixedPrincipal { frame },
                    strains.into_iter().zip(selections).collect(),
                )
            }
        };
        if selected
            .iter()
            .any(|(key, _)| !scalar_role(key, Dims::NONE))
        {
            return Err(invalid(
                "engineering length changes require dimensionless scalar source roles",
            ));
        }
        for point in [query.lower, query.upper, query.reference] {
            admit_temperature_point(point, query.temperature_axis)?;
            if point.axis_quantities() != query.lower.axis_quantities()
                || !point.axes().keys().eq(query.lower.axes().keys())
                || point.axes().iter().any(|(axis, value)| {
                    axis != query.temperature_axis && query.lower.axes()[axis] != *value
                })
            {
                return Err(invalid(
                    "reference and interval must share every non-temperature coordinate and descriptor",
                ));
            }
        }
        let low = query.lower.axes()[query.temperature_axis];
        let high = query.upper.axes()[query.temperature_axis];
        let reference_t = query.reference.axes()[query.temperature_axis];
        if low >= high || reference_t < low || reference_t > high {
            return Err(invalid(
                "free expansion reference must lie inside a positive-width temperature interval",
            ));
        }
        let sources: Vec<_> = selected
            .into_iter()
            .map(|(key, selection)| {
                source_line(
                    card,
                    key,
                    query.temperature_axis,
                    query.lower,
                    query.upper,
                    selection,
                )
            })
            .collect::<Result<_, _>>()?;
        let strains_at = |t| std::array::from_fn(|i| source_value(&sources[i % sources.len()], t));
        let reference_lengths = strains_at(reference_t).map(|strain| 1.0 + strain);
        if !reference_lengths.into_iter().all(positive) {
            return Err(invalid(
                "source reference lengths must be finite and positive relative to their datum",
            ));
        }
        let mut temperatures: Vec<_> = sources
            .iter()
            .flat_map(|source| source.knots.iter().map(|(t, _)| *t))
            .chain([reference_t])
            .collect();
        temperatures.sort_by(f64::total_cmp);
        temperatures.dedup();
        if temperatures.len() > MAX_FREE_EXPANSION_KNOTS {
            return Err(invalid(
                "free expansion source union exceeds its breakpoint budget",
            ));
        }
        let material_card = card.content_hash();
        let mut identity = identity_prefix(
            "org.frankensim.fs-material.prepared-free-expansion.v1",
            material_card,
            query.temperature_axis,
            basis,
        );
        identity.update(&reference_t.to_bits().to_le_bytes());
        hash_sources(&mut identity, &sources);
        let identity = identity.finalize();
        let mut knots = Vec::with_capacity(temperatures.len());
        for temperature in temperatures {
            let strains = strains_at(temperature);
            checked_geometry(identity, temperature, strains, reference_lengths)?;
            knots.push(StrainKnot {
                temperature,
                strains,
            });
        }
        Ok(Self {
            material_card,
            basis,
            temperature_axis: query.temperature_axis.to_owned(),
            lower: query.lower.clone(),
            upper: query.upper.clone(),
            reference: query.reference.clone(),
            sources,
            reference_lengths,
            knots,
            identity,
        })
    }

    /// One binary search; no allocation, hashing, query, mutation or extrapolation.
    /// Interpolate the independent engineering strains, then form normalized
    /// stretches and their product. The reference is a union knot, so its unit
    /// stretches are exact even when the original fitted datum residual is nonzero.
    /// Preparation checks all knots, but extreme interior products or ratios can
    /// still overflow or underflow; that sample refuses without modifying the path.
    pub fn at_temperature(
        &self,
        temperature_k: f64,
    ) -> Result<FreeExpansionGeometrySample, ThermalTransportError> {
        let span = self.temperature_span_k();
        if !temperature_k.is_finite() || temperature_k < span[0] || temperature_k > span[1] {
            return Err(ThermalTransportError::Invalid(
                "temperature is outside the prepared free expansion interval",
            ));
        }
        let index = self
            .knots
            .partition_point(|knot| knot.temperature < temperature_k);
        let upper = self.knots[index];
        let strains = if upper.temperature == temperature_k {
            upper.strains
        } else {
            let lower = self.knots[index - 1];
            let weight =
                (temperature_k - lower.temperature) / (upper.temperature - lower.temperature);
            std::array::from_fn(|i| {
                lower.strains[i] + weight * (upper.strains[i] - lower.strains[i])
            })
        };
        checked_geometry(
            self.identity,
            temperature_k,
            strains,
            self.reference_lengths,
        )
    }

    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    #[must_use]
    pub const fn material_card_identity(&self) -> ContentHash {
        self.material_card
    }
    #[must_use]
    pub const fn basis(&self) -> FreeExpansionBasis {
        self.basis
    }
    #[must_use]
    pub fn temperature_axis(&self) -> &str {
        &self.temperature_axis
    }
    #[must_use]
    pub fn query_endpoints(&self) -> [&QueryPoint; 2] {
        [&self.lower, &self.upper]
    }
    #[must_use]
    pub const fn reference_query_point(&self) -> &QueryPoint {
        &self.reference
    }
    /// Original isotropic line, or three principal lines in declared axis order.
    #[must_use]
    pub fn source_lines(&self) -> &[PropertyLineSlice] {
        &self.sources
    }
    #[must_use]
    pub fn temperature_span_k(&self) -> [f64; 2] {
        [
            self.knots[0].temperature,
            self.knots[self.knots.len() - 1].temperature,
        ]
    }
    /// Union of original source breakpoints plus the exact geometry reference.
    pub fn temperatures_k(&self) -> impl ExactSizeIterator<Item = f64> + '_ {
        self.knots.iter().map(|knot| knot.temperature)
    }
}

/// Constant-mass density composed with the original free geometry.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreeExpansionDensitySample {
    path_identity: ContentHash,
    geometry: FreeExpansionGeometrySample,
    density_kg_m3: f64,
}

impl FreeExpansionDensitySample {
    #[must_use]
    pub const fn path_identity(&self) -> ContentHash {
        self.path_identity
    }
    #[must_use]
    pub const fn temperature_k(&self) -> f64 {
        self.geometry.temperature_k
    }
    #[must_use]
    pub const fn principal_stretches(&self) -> [f64; 3] {
        self.geometry.principal_stretches
    }
    #[must_use]
    pub const fn relative_volume(&self) -> f64 {
        self.geometry.relative_volume
    }
    #[must_use]
    pub const fn density_kg_m3(&self) -> f64 {
        self.density_kg_m3
    }
    /// Original geometry result retains its own source-path identity.
    #[must_use]
    pub const fn geometry(&self) -> FreeExpansionGeometrySample {
        self.geometry
    }
}

/// Optional density composition; geometric expansion remains independently usable.
#[derive(Clone, Debug, PartialEq)]
pub struct PreparedFreeExpansionDensity {
    expansion: PreparedFreeExpansion,
    density: MaterialAnswer,
    identity: ContentHash,
}

impl PreparedFreeExpansionDensity {
    /// Prepare the original free geometry and resolve density only at its exact
    /// reference. This preserves the existing constructor and density identity.
    pub fn try_new(
        card: &MaterialCard,
        query: FreeExpansionSourceQuery<'_>,
    ) -> Result<Self, ThermalTransportError> {
        let expansion = PreparedFreeExpansion::try_new(
            card,
            FreeExpansionGeometryQuery {
                reference: query.reference,
                strains: query.strains,
                temperature_axis: query.temperature_axis,
                lower: query.lower,
                upper: query.upper,
            },
        )?;
        Self::try_from_expansion(card, expansion, query.density, query.density_selection)
    }

    /// Reuse already prepared geometry without re-querying any strain source.
    /// Density must belong to the same immutable card and support the geometry's
    /// complete reference context. This does not assert matched specimens or
    /// manufacture a missing density, phase law, mass loss or covariance model.
    pub fn try_from_expansion(
        card: &MaterialCard,
        expansion: PreparedFreeExpansion,
        density_key: &PropertyKey,
        selection: ClaimSelection,
    ) -> Result<Self, ThermalTransportError> {
        if card.content_hash() != expansion.material_card_identity()
            || !scalar_role(density_key, super::DENSITY)
        {
            return Err(ThermalTransportError::Invalid(
                "density composition requires the original card and a scalar density source",
            ));
        }
        let reference = expansion.reference_query_point();
        let density = source_line(
            card,
            density_key,
            expansion.temperature_axis(),
            reference,
            reference,
            selection,
        )?
        .knots
        .pop()
        .expect("canonical point line is nonempty")
        .1;
        if !positive(density.evidence.value.value) {
            return Err(ThermalTransportError::Invalid(
                "reference density must be finite and positive",
            ));
        }
        // Preserve the established density identity byte order and domain.
        let mut identity = identity_prefix(
            "org.frankensim.fs-material.prepared-free-expansion-density.v1",
            expansion.material_card_identity(),
            expansion.temperature_axis(),
            expansion.basis(),
        );
        identity.update(density.receipt.content_hash().as_bytes());
        hash_sources(&mut identity, expansion.source_lines());
        for knot in &expansion.knots {
            let geometry = checked_geometry(
                expansion.identity,
                knot.temperature,
                knot.strains,
                expansion.reference_lengths,
            )?;
            density_from_geometry(&geometry, density.evidence.value.value)?;
        }
        Ok(Self {
            expansion,
            density,
            identity: identity.finalize(),
        })
    }

    /// One original geometry evaluation followed by rho_ref/J. No allocation,
    /// database query or hashing occurs during successful evaluation.
    pub fn at_temperature(
        &self,
        temperature_k: f64,
    ) -> Result<FreeExpansionDensitySample, ThermalTransportError> {
        let geometry = self.expansion.at_temperature(temperature_k)?;
        let density_kg_m3 = density_from_geometry(&geometry, self.density.evidence.value.value)?;
        Ok(FreeExpansionDensitySample {
            path_identity: self.identity,
            geometry,
            density_kg_m3,
        })
    }

    #[must_use]
    pub const fn expansion(&self) -> &PreparedFreeExpansion {
        &self.expansion
    }
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    #[must_use]
    pub const fn material_card_identity(&self) -> ContentHash {
        self.expansion.material_card_identity()
    }
    #[must_use]
    pub const fn basis(&self) -> FreeExpansionBasis {
        self.expansion.basis()
    }
    #[must_use]
    pub fn temperature_axis(&self) -> &str {
        self.expansion.temperature_axis()
    }
    #[must_use]
    pub fn query_endpoints(&self) -> [&QueryPoint; 2] {
        self.expansion.query_endpoints()
    }
    #[must_use]
    pub const fn reference_query_point(&self) -> &QueryPoint {
        self.expansion.reference_query_point()
    }
    #[must_use]
    pub const fn reference_density_answer(&self) -> &MaterialAnswer {
        &self.density
    }
    #[must_use]
    pub fn source_lines(&self) -> &[PropertyLineSlice] {
        self.expansion.source_lines()
    }
    #[must_use]
    pub fn temperature_span_k(&self) -> [f64; 2] {
        self.expansion.temperature_span_k()
    }
    pub fn temperatures_k(&self) -> impl ExactSizeIterator<Item = f64> + '_ {
        self.expansion.temperatures_k()
    }
}

fn identity_prefix(
    domain: &str,
    card: ContentHash,
    axis: &str,
    basis: FreeExpansionBasis,
) -> DomainHasher {
    let mut identity = DomainHasher::new(domain);
    identity.update(card.as_bytes());
    identity.update(&(axis.len() as u64).to_le_bytes());
    identity.update(axis.as_bytes());
    match basis {
        FreeExpansionBasis::Isotropic => {
            identity.update(&[0]);
        }
        FreeExpansionBasis::FixedPrincipal { frame } => {
            identity.update(&[1]);
            identity.update(frame.as_bytes());
        }
    }
    identity
}

fn hash_sources(identity: &mut DomainHasher, sources: &[PropertyLineSlice]) {
    for source in sources {
        identity.update(&(source.knots.len() as u64).to_le_bytes());
        for (_, answer) in &source.knots {
            identity.update(answer.receipt.content_hash().as_bytes());
        }
    }
}

fn density_from_geometry(
    geometry: &FreeExpansionGeometrySample,
    reference_density: f64,
) -> Result<f64, ThermalTransportError> {
    let density = reference_density / geometry.relative_volume;
    if !positive(density) {
        return Err(ThermalTransportError::Invalid(
            "derived free expansion density must be finite and positive",
        ));
    }
    Ok(density)
}

fn scalar_role(key: &PropertyKey, dims: Dims) -> bool {
    key.quantity() == QuantitySpec::dimensional(dims)
        && key.hardness_test().is_none()
        && key.elastic_component().is_none()
        && key.strain_component().is_none()
        && key.stress_component().is_none()
}

fn checked_geometry(
    identity: ContentHash,
    temperature_k: f64,
    strains: [f64; 3],
    reference_lengths: [f64; 3],
) -> Result<FreeExpansionGeometrySample, ThermalTransportError> {
    let principal_stretches = std::array::from_fn(|i| (1.0 + strains[i]) / reference_lengths[i]);
    // Multiply the smallest and largest first: principal-axis order must not
    // create intermediate overflow/underflow when the final J is representable.
    let mut ordered = principal_stretches;
    ordered.sort_unstable_by(f64::total_cmp);
    let relative_volume = (ordered[0] * ordered[2]) * ordered[1];
    if !principal_stretches
        .into_iter()
        .chain([relative_volume])
        .all(positive)
    {
        return Err(ThermalTransportError::Invalid(
            "free expansion stretches and volume must be finite and positive",
        ));
    }
    Ok(FreeExpansionGeometrySample {
        path_identity: identity,
        temperature_k,
        principal_stretches,
        relative_volume,
    })
}
