//! Source-bound singlet/triplet kinetics for phosphorescence and delayed fluorescence.
//!
//! Original material claims provide both state energies and all six competing
//! first-order rates. The numerical law is the existing
//! [`SingletTripletDecay`]; this adapter neither fits lifetimes nor infers a
//! quantum yield from a material name. Rates may be source curves or surfaces
//! in temperature and explicitly typed fixed external conditions, including
//! quencher concentration. No oxygen-solubility or partial-pressure conversion
//! is assumed, and a nonradiative rate already containing quenching must not
//! receive the same channel again. Nonradiative rates retain the native model's
//! local-thermalization meaning; no excited-quencher store or photochemistry is
//! inferred. Thermally admissible reverse rates remain the source's responsibility.
//!
//! A snapshot holds its complete physical context fixed during an interval.
//! Refreshing reuses the same immutable card and selected claims. An owner of
//! excited populations must keep state energies fixed or account explicitly
//! for their parameter work; changing a snapshot does not update stored energy.
//! Pump absorption, photon spectra and photon transport remain separate owners.
//! This dilute two-state model does not describe trap distributions, saturation,
//! triplet annihilation, photobleaching or coherent spin dynamics.

use crate::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::{Dims, QuantitySpec};
use fs_thermal::luminescence::transfer::{SingletTripletDecay, SingletTripletError};

/// Source-bound triplet pair encounters and productive singlet branching.
pub mod annihilation;
/// Source-bound molecular absorption into the singlet and triplet states.
pub mod pump;

/// Original `[singlet, triplet]` excitation energies relative to one ground-state
/// datum, in joules per excitation, not joules per mole or activation barriers.
pub const PHOSPHORESCENCE_ENERGY_PROPERTIES: [&str; 2] = [
    "phosphorescence_singlet_energy",
    "phosphorescence_triplet_energy",
];
/// Original radiative ground-state rates, ordered `[singlet, triplet]`, in s^-1.
pub const PHOSPHORESCENCE_RADIATIVE_RATE_PROPERTIES: [&str; 2] = [
    "phosphorescence_singlet_radiative_rate",
    "phosphorescence_triplet_radiative_rate",
];
/// Original nonradiative ground-state rates, ordered `[singlet, triplet]`, in s^-1.
pub const PHOSPHORESCENCE_NONRADIATIVE_RATE_PROPERTIES: [&str; 2] = [
    "phosphorescence_singlet_nonradiative_rate",
    "phosphorescence_triplet_nonradiative_rate",
];
/// Original `[singlet -> triplet, triplet -> singlet]` rates, in s^-1.
pub const PHOSPHORESCENCE_CROSSING_RATE_PROPERTIES: [&str; 2] = [
    "phosphorescence_intersystem_crossing_rate",
    "phosphorescence_reverse_crossing_rate",
];
/// SI dimensions of an individual excitation's energy, joules.
pub const PHOSPHORESCENCE_ENERGY_DIMS: Dims = Dims([2, 1, -2, 0, 0, 0]);
/// SI dimensions of a first-order transition rate, inverse seconds.
pub const PHOSPHORESCENCE_RATE_DIMS: Dims = Dims([0, 0, -1, 0, 0, 0]);

/// Original source, scalar-context or kinetic admission failure.
#[derive(Clone, Debug, PartialEq)]
pub enum PhosphorescenceError {
    /// Source selection, source support or query-coordinate failure.
    Material(MaterialStatePointError),
    /// A canonical role is absent, ambiguous or not a supported scalar snapshot.
    Property {
        /// Canonical role requiring attention.
        name: &'static str,
        /// Stable reason; no default coefficient is substituted.
        reason: &'static str,
    },
    /// Energies or competing rates violate the shared kinetic model's domain.
    Kinetics(SingletTripletError),
}

impl core::fmt::Display for PhosphorescenceError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(formatter, "phosphorescence material: {self:?}")
    }
}
impl std::error::Error for PhosphorescenceError {}
impl From<MaterialStatePointError> for PhosphorescenceError {
    fn from(value: MaterialStatePointError) -> Self {
        Self::Material(value)
    }
}
impl From<SingletTripletError> for PhosphorescenceError {
    fn from(value: SingletTripletError) -> Self {
        Self::Kinetics(value)
    }
}

/// A complete original material snapshot lowered into the native two-state law.
/// Identities bind sources and conditions; they are not experimental validation.
#[derive(Clone, Debug, PartialEq)]
pub struct PhosphorescenceStatePoint {
    resolved: ResolvedMaterialStatePoint,
    decay: SingletTripletDecay,
    identity: ContentHash,
}

impl PhosphorescenceStatePoint {
    /// Resolve all eight required source roles at one complete physical point.
    /// Energies are strictly positive; an explicitly zero rate is permitted.
    pub fn resolve(
        card: &MaterialCard,
        point: &QueryPoint,
        selection: MaterialPropertySelection,
    ) -> Result<Self, PhosphorescenceError> {
        let mut requirements = Vec::with_capacity(8);
        for name in PHOSPHORESCENCE_ENERGY_PROPERTIES {
            requirements.push(ScalarPropertyRequirement::try_new(
                name,
                PHOSPHORESCENCE_ENERGY_DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )?);
        }
        for name in PHOSPHORESCENCE_RADIATIVE_RATE_PROPERTIES
            .into_iter()
            .chain(PHOSPHORESCENCE_NONRADIATIVE_RATE_PROPERTIES)
            .chain(PHOSPHORESCENCE_CROSSING_RATE_PROPERTIES)
        {
            requirements.push(ScalarPropertyRequirement::try_new(
                name,
                PHOSPHORESCENCE_RATE_DIMS,
                ScalarAdmissibility::NonNegative,
            )?);
        }
        Self::try_from_material(resolve_material_state_point(
            card,
            point,
            &requirements,
            selection,
        )?)
    }

    /// Adopt a complete original bundle, retaining unrelated optical, thermal
    /// and mechanical properties. Exact scalar quantities and source contexts
    /// are rechecked even if the original bundle used permissive requirements.
    pub fn try_from_material(
        resolved: ResolvedMaterialStatePoint,
    ) -> Result<Self, PhosphorescenceError> {
        let pair = |names: [&'static str; 2], dims| {
            Ok::<_, PhosphorescenceError>([
                coefficient(&resolved, names[0], dims)?,
                coefficient(&resolved, names[1], dims)?,
            ])
        };
        let decay = SingletTripletDecay::try_new(
            pair(
                PHOSPHORESCENCE_ENERGY_PROPERTIES,
                PHOSPHORESCENCE_ENERGY_DIMS,
            )?,
            pair(
                PHOSPHORESCENCE_RADIATIVE_RATE_PROPERTIES,
                PHOSPHORESCENCE_RATE_DIMS,
            )?,
            pair(
                PHOSPHORESCENCE_NONRADIATIVE_RATE_PROPERTIES,
                PHOSPHORESCENCE_RATE_DIMS,
            )?,
            pair(
                PHOSPHORESCENCE_CROSSING_RATE_PROPERTIES,
                PHOSPHORESCENCE_RATE_DIMS,
            )?,
        )?;
        let mut hash = DomainHasher::new("org.frankensim.material.phosphorescence-state.v1");
        hash.update(resolved.identity().as_bytes());
        Ok(Self {
            resolved,
            decay,
            identity: hash.finalize(),
        })
    }

    /// Original immutable card, complete physical point and selected receipts.
    #[must_use]
    pub const fn resolved(&self) -> &ResolvedMaterialStatePoint {
        &self.resolved
    }

    /// Native fixed-rate model; no source query or allocation occurs here.
    #[must_use]
    pub const fn decay(&self) -> &SingletTripletDecay {
        &self.decay
    }

    /// Identity of this two-state model family and its entire original bundle.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Refresh exact original claims and coordinate meanings atomically.
    /// Unsupported temperatures or quencher conditions refuse. The returned
    /// snapshot owns no excitation history and supplies no unaccounted work
    /// when source energies change; an evolving physical owner must admit that
    /// change separately before replacing its coefficients.
    pub fn at_query_point(
        &self,
        card: &MaterialCard,
        point: &QueryPoint,
    ) -> Result<Self, PhosphorescenceError> {
        Self::try_from_material(self.resolved.at_query_point(card, point)?)
    }
}

fn coefficient(
    material: &ResolvedMaterialStatePoint,
    name: &'static str,
    dims: Dims,
) -> Result<f64, PhosphorescenceError> {
    let invalid = |reason| PhosphorescenceError::Property { name, reason };
    let property = material
        .property(name)
        .ok_or_else(|| invalid("missing or ambiguous original property"))?;
    let requirement = property.requirement();
    if requirement.quantity() != QuantitySpec::dimensional(dims)
        || requirement.hardness_test().is_some()
        || requirement.elastic_component().is_some()
        || requirement.strain_component().is_some()
        || requirement.stress_component().is_some()
        || !property.is_static_at_fixed_external_conditions()
    {
        return Err(invalid(
            "requires an exact SI scalar at fixed external conditions, without frequency, field, strain, stress or history dependence or bounds",
        ));
    }
    Ok(property.value_si())
}
