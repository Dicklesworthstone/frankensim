//! Energy-dependent deposition of an explicitly supplied radiation spectrum.
//!
//! Each physical source has one independently closed deposited-energy response
//! per energy group. Different groups may heat different receivers, so changing
//! radioactive composition can change deposition without replacing a response.
//! Coefficients must describe the declared radiation, geometry, materials and
//! within-group spectra; this module does not derive them from attenuation or
//! thermal emissivity. Grouping itself does not establish transport accuracy.
//!
//! Canonical source/group pairs are mapped to the original [`RadiationResponse`]
//! owner. Its sparse receiver reduction, closure checks, arithmetic accounting
//! and cancellation behavior remain authoritative. Input values are integrated
//! emitted joules, not particle counts or densities. Applying this map creates
//! no new energy allowance: coupling must reconcile its input with the original
//! outward-radiation source account before crediting material heat.

use std::sync::Arc;

use super::{
    RadiationDeposit, RadiationDeposition, RadiationError, RadiationResponse,
    RadiationResponseLimits, RadiationSourceResponse, nonnegative,
};
use crate::endf::decay::spectra::DecayRadiation;

/// Supplied energy fractions for one physical source and one energy group.
#[derive(Clone, Debug, PartialEq)]
pub struct RadiationGroupResponse {
    /// Energy thermalized at each receiver, as a fraction of this group's input.
    pub deposits: Vec<RadiationDeposit>,
    /// Group energy escaping the modeled receiving domain.
    pub escape_fraction: f64,
}

/// One physical emitter's responses in the declared energy-group order.
#[derive(Clone, Debug, PartialEq)]
pub struct GroupedRadiationSourceResponse {
    /// Original emitting identity, independent of the internal group-row indices.
    pub source_id: u64,
    /// Exactly one response for every interval of the supplied group grid.
    pub groups: Vec<RadiationGroupResponse>,
}

/// Explicit admission bounds for physical sources and the expanded sparse map.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct GroupedRadiationResponseLimits {
    /// Maximum number of distinct physical emitting sources.
    pub max_sources: usize,
    /// Maximum number of energy intervals per physical source.
    pub max_groups: usize,
    /// Maximum number of distinct receiver identities across every group.
    pub max_receivers: usize,
    /// Maximum `sources * groups`; checked before allocating internal rows.
    pub max_source_groups: usize,
    /// Maximum deposited-energy entries across all source/group responses.
    pub max_transfers: usize,
}

/// Immutable grouped response. Clones share coefficients, group grid and mapping.
#[derive(Clone, Debug)]
pub struct GroupedRadiationResponse(Arc<GroupedResponseData>);

#[derive(Clone, Debug)]
struct GroupedResponseData {
    radiation: DecayRadiation,
    edges_j: Vec<f64>,
    source_ids: Vec<u64>,
    response: RadiationResponse,
    preparation_roundoff: Vec<f64>,
}

impl GroupedRadiationResponse {
    /// Admit a finite, nonnegative, strictly increasing SI energy grid and a
    /// complete response for every source/group pair. Fractions are never
    /// renormalized. Their physical applicability remains the caller's claim.
    pub fn new(
        data_id: String,
        radiation: DecayRadiation,
        group_edges_j: Vec<f64>,
        mut rows: Vec<GroupedRadiationSourceResponse>,
        limits: GroupedRadiationResponseLimits,
    ) -> Result<Self, RadiationError> {
        let groups = group_edges_j
            .len()
            .checked_sub(1)
            .filter(|&count| count > 0)
            .ok_or(RadiationError::InvalidInput(
                "at least two grouped radiation energy edges are required",
            ))?;
        if groups > limits.max_groups {
            return Err(RadiationError::Budget("energy group count"));
        }
        if !group_edges_j.iter().copied().all(nonnegative)
            || !group_edges_j.windows(2).all(|pair| pair[0] < pair[1])
        {
            return Err(RadiationError::InvalidInput(
                "grouped radiation energies must be finite, nonnegative and ascending",
            ));
        }
        if rows.len() > limits.max_sources {
            return Err(RadiationError::Budget("physical source count"));
        }
        let row_count = rows
            .len()
            .checked_mul(groups)
            .ok_or(RadiationError::Budget("source/group count"))?;
        if row_count > limits.max_source_groups || u64::try_from(row_count).is_err() {
            return Err(RadiationError::Budget("source/group count"));
        }
        let mut transfers = 0_usize;
        for row in &rows {
            if row.groups.len() != groups {
                return Err(RadiationError::InvalidInput(
                    "each radiation source requires exactly one response per energy group",
                ));
            }
            for group in &row.groups {
                transfers = transfers
                    .checked_add(group.deposits.len())
                    .ok_or(RadiationError::Budget("transfer count"))?;
                if transfers > limits.max_transfers {
                    return Err(RadiationError::Budget("transfer count"));
                }
            }
        }
        rows.sort_by_key(|row| row.source_id);
        if rows
            .windows(2)
            .any(|pair| pair[0].source_id == pair[1].source_id)
        {
            return Err(RadiationError::InvalidInput(
                "duplicate physical source identity",
            ));
        }
        let mut source_ids = Vec::with_capacity(rows.len());
        let mut expanded = Vec::with_capacity(row_count);
        for row in rows {
            source_ids.push(row.source_id);
            for group in row.groups {
                // row_count was checked to fit both usize and u64 above.
                expanded.push(RadiationSourceResponse {
                    source_id: expanded.len() as u64,
                    deposits: group.deposits,
                    escape_fraction: group.escape_fraction,
                });
            }
        }
        let response = RadiationResponse::new(
            data_id,
            expanded,
            RadiationResponseLimits {
                max_sources: limits.max_source_groups,
                max_receivers: limits.max_receivers,
                max_transfers: limits.max_transfers,
            },
        )?;
        Ok(Self(Arc::new(GroupedResponseData {
            radiation,
            edges_j: group_edges_j,
            source_ids,
            response,
            preparation_roundoff: vec![0.0; row_count],
        })))
    }

    /// Retain generated-response arithmetic defects without changing its rows.
    pub(crate) fn with_preparation_roundoff(
        mut self,
        errors: Vec<(u64, Vec<f64>)>,
    ) -> Result<Self, RadiationError> {
        let groups = self.0.edges_j.len() - 1;
        if errors.len() != self.0.source_ids.len() {
            return Err(RadiationError::InvalidInput("exact preparation source coverage required"));
        }
        let mut expanded = vec![0.0; self.0.preparation_roundoff.len()];
        let mut seen = vec![false; self.0.source_ids.len()];
        for (id, values) in errors {
            let source = self.0.source_ids.binary_search(&id)
                .map_err(|_| RadiationError::InvalidInput("unknown preparation source"))?;
            if seen[source] || values.len() != groups || !values.iter().copied().all(nonnegative) {
                return Err(RadiationError::InvalidInput("exact nonnegative preparation group coverage required"));
            }
            seen[source] = true;
            expanded[source * groups..(source + 1) * groups].copy_from_slice(&values);
        }
        let data = Arc::make_mut(&mut self.0);
        data.response.add_preparation_roundoff(&expanded)?;
        for (total, extra) in data.preparation_roundoff.iter_mut().zip(expanded) {
            let mut error = super::ErrorBudget::new(f64::MAX);
            error.charge(*total)?;
            error.charge(extra)?;
            *total = error.used;
        }
        Ok(self)
    }

    pub(crate) fn preparation_roundoff(&self, source_id: u64, group: usize) -> Option<f64> {
        let groups = self.0.edges_j.len() - 1;
        if group >= groups { return None; }
        let source = self.0.source_ids.binary_search(&source_id).ok()?;
        Some(self.0.preparation_roundoff[source * groups + group])
    }

    /// Borrow one canonical coefficient row without applying unit-energy inputs.
    pub(crate) fn source_group_row(
        &self,
        source_id: u64,
        group: usize,
    ) -> Option<&RadiationSourceResponse> {
        let groups = self.0.edges_j.len() - 1;
        if group >= groups {
            return None;
        }
        let source = self.0.source_ids.binary_search(&source_id).ok()?;
        self.0.response.rows().get(source * groups + group)
    }

    /// Supplied dataset identity, including its intended physical applicability.
    #[must_use]
    pub fn data_id(&self) -> &str {
        self.0.response.data_id()
    }

    /// Radiation kind to which all supplied response coefficients apply.
    #[must_use]
    pub fn radiation(&self) -> DecayRadiation {
        self.0.radiation
    }

    /// Original SI energy edges; group order is never sorted or reassigned.
    #[must_use]
    pub fn group_edges_j(&self) -> &[f64] {
        &self.0.edges_j
    }

    /// Original physical emitting identities in ascending order.
    pub fn source_ids(&self) -> impl ExactSizeIterator<Item = u64> + '_ {
        self.0.source_ids.iter().copied()
    }

    /// Distinct physical receivers across all energy groups, in ascending order.
    pub fn receiver_ids(&self) -> impl ExactSizeIterator<Item = u64> + '_ {
        self.0.response.receiver_ids()
    }

    /// Distribute integrated emitted energies on this response's exact radiation
    /// kind and group grid. Every physical source must appear exactly once,
    /// including zero emitters, with one finite nonnegative energy per group.
    ///
    /// The input order is irrelevant. All sums use canonical physical-source,
    /// energy-group and receiver order in the original scalar response owner.
    /// The caller must establish spectral-grid identity and original source
    /// energy coverage; this method neither regroups nor rescales supplied joules.
    /// Refusal and cancellation return no partial receiver credit or mutable state.
    pub fn apply(
        &self,
        emitted_by_source_group_j: &[(u64, &[f64])],
        energy_tolerance_j: f64,
        mut checkpoint: impl FnMut() -> Result<(), RadiationError>,
    ) -> Result<RadiationDeposition, RadiationError> {
        if !nonnegative(energy_tolerance_j) {
            return Err(RadiationError::InvalidInput(
                "energy tolerance must be finite and nonnegative",
            ));
        }
        if emitted_by_source_group_j.len() != self.0.source_ids.len() {
            return Err(RadiationError::InvalidInput(
                "physical radiation source coverage must be exact",
            ));
        }
        checkpoint()?;
        let groups = self.0.edges_j.len() - 1;
        let mut seen = vec![false; self.0.source_ids.len()];
        let mut expanded = Vec::with_capacity(self.0.response.rows().len());
        for &(id, energies) in emitted_by_source_group_j {
            checkpoint()?;
            let source = self.0.source_ids.binary_search(&id).map_err(|_| {
                RadiationError::InvalidInput("unknown physical radiation source identity")
            })?;
            if seen[source] {
                return Err(RadiationError::InvalidInput(
                    "duplicate physical emitted source identity",
                ));
            }
            seen[source] = true;
            if energies.len() != groups {
                return Err(RadiationError::InvalidInput(
                    "each radiation source requires exactly one emitted energy per group",
                ));
            }
            for (group, &energy) in energies.iter().enumerate() {
                checkpoint()?;
                expanded.push(((source * groups + group) as u64, energy));
            }
        }
        self.0
            .response
            .apply(&expanded, energy_tolerance_j, checkpoint)
    }
}

#[cfg(test)]
mod tests;
