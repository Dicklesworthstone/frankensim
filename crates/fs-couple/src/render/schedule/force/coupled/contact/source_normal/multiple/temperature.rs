//! Atomic mixed JKR/DMT/Maugis source changes at held native coordinates.
//!
//! Only selected contact sources change. Every original modal, connection,
//! material and radiation state, the clock, and the last mechanical frame stay
//! fixed. Signed parameter work belongs to the explicit temperature driver; it
//! does not become caloric heat or a mechanical impulse. Source laws remain
//! frozen within each subsequent mechanical step.
//!
//! Requests are evaluated in contact-index order. Each sealed source change
//! passes its own small-work gate before the complete successor set passes the
//! same coupled negative-curvature screen used by mechanical stepping. A failed
//! source, receiver or cancellation publishes none of the replacements.

use super::super::temperature::{
    SourceNormalThermalChange, prepare_dmt_contact_refresh, prepare_jkr_contact_refresh,
    prepare_maugis_contact_refresh,
};
use super::{
    BoundNormalContact, CancelGate, MultipleSourceNormalModalSystem, SourceNormalModalError,
    check_contact_work, check_coupled_stability, finite, invalid, limit, negative_slope_root,
    observation, poll,
};
use fs_blake3::ContentHash;
use fs_contact::normal_patch::adhesive::temperature::{
    SourceDmtThermalRefresh, SourceJkrThermalRefresh, SourceMaugisThermalRefresh,
};
use fs_math::dd::Dd;

/// One original sealed change for one contact. Indices must be distinct and
/// in range; request order does not change evaluation or publication order.
#[derive(Clone, Copy)]
pub enum NormalSourceRefresh<'a> {
    Jkr {
        contact: usize,
        refreshed: &'a SourceJkrThermalRefresh,
    },
    Dmt {
        contact: usize,
        refreshed: &'a SourceDmtThermalRefresh,
    },
    Maugis {
        contact: usize,
        refreshed: &'a SourceMaugisThermalRefresh,
    },
}
impl NormalSourceRefresh<'_> {
    pub const fn contact_index(&self) -> usize {
        match self {
            Self::Jkr { contact, .. }
            | Self::Dmt { contact, .. }
            | Self::Maugis { contact, .. } => *contact,
        }
    }
}

/// One source's independently checked contribution to the accepted batch.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IndexedSourceNormalThermalChange {
    pub contact: usize,
    pub original_model_identity: ContentHash,
    pub model_identity: ContentHash,
    pub change: SourceNormalThermalChange,
}

/// One accepted batch at unchanged native coordinates and clock. Each selected
/// source retains its own identity and signed-work contribution below.
#[derive(Clone, Debug, PartialEq)]
pub struct MultipleSourceNormalRefreshReport {
    pub sample: u64,
    /// Unchanged original native storage, including all retained memories.
    pub network_energy_j: f64,
    /// All original contact potentials, including contacts outside the batch.
    pub contact_energy_before_j: f64,
    pub contact_energy_after_j: f64,
    pub energy_before_j: f64,
    pub energy_after_j: f64,
    /// Sum of the original signed parameter works; no mechanical motion occurs.
    pub parameter_work_j: f64,
    /// Sum of changed-source potential increments minus their signed work.
    /// This gate never subtracts the large total native energy observations.
    pub energy_residual_j: f64,
    /// Sum of the individually enforced contact-sized tolerances.
    pub energy_tolerance_j: f64,
    /// Changed contacts only, in canonical contact-index order.
    pub contacts: Vec<IndexedSourceNormalThermalChange>,
}

type Replacements = Vec<(usize, BoundNormalContact)>;

impl MultipleSourceNormalModalSystem {
    /// Apply a nonempty bounded set of original held-indentation refreshes as
    /// one source transition. Neither geometry nor any native history advances.
    /// The explicit driver must account for the returned signed parameter work.
    ///
    /// Request ordering, replacement bindings and reports allocate in proportion
    /// to the selected contacts. The native network and its histories are never
    /// cloned; the complete stability screen reuses the existing contact scratch.
    pub fn apply_source_refreshes(
        &mut self,
        refreshes: &[NormalSourceRefresh<'_>],
        gate: &CancelGate,
    ) -> Result<MultipleSourceNormalRefreshReport, SourceNormalModalError> {
        let (replacements, report) = self.prepare_source_refreshes(refreshes, gate)?;
        poll(Some(gate))?;
        self.publish_source_refreshes(replacements);
        Ok(report)
    }

    /// Stage a receiving history against the complete candidate batch report and
    /// the replacement bindings in contact-index order. The callback must return
    /// a candidate whose later installation is infallible, as in native stepping.
    /// A late refusal or cancellation preserves every accepted source and frame.
    pub fn apply_source_refreshes_with_admission<R, E>(
        &mut self,
        refreshes: &[NormalSourceRefresh<'_>],
        gate: &CancelGate,
        admit: impl FnOnce(
            &MultipleSourceNormalRefreshReport,
            &[(usize, BoundNormalContact)],
        ) -> Result<R, E>,
    ) -> Result<R, E>
    where
        E: From<SourceNormalModalError>,
    {
        let (replacements, report) = self
            .prepare_source_refreshes(refreshes, gate)
            .map_err(E::from)?;
        let receiving = admit(&report, &replacements)?;
        poll(Some(gate))
            .map_err(SourceNormalModalError::from)
            .map_err(E::from)?;
        self.publish_source_refreshes(replacements);
        Ok(receiving)
    }

    fn publish_source_refreshes(&mut self, replacements: Replacements) {
        // All indices, laws, energy gates and receiver state have passed. This
        // bounded move contains no fallible operation or partial publication.
        for (index, replacement) in replacements {
            self.points[index].contact = replacement;
        }
    }

    fn prepare_source_refreshes(
        &mut self,
        refreshes: &[NormalSourceRefresh<'_>],
        gate: &CancelGate,
    ) -> Result<(Replacements, MultipleSourceNormalRefreshReport), SourceNormalModalError> {
        poll(Some(gate))?;
        if refreshes.is_empty() || refreshes.len() > self.points.len() {
            return Err(invalid("nonempty bounded original-source refresh batch required").into());
        }
        // A canonical order makes signed-work accumulation and receiving
        // candidates independent of how the caller enumerates the batch.
        let mut ordered = refreshes.to_vec();
        ordered.sort_unstable_by_key(NormalSourceRefresh::contact_index);
        if ordered.last().unwrap().contact_index() >= self.points.len()
            || ordered
                .windows(2)
                .any(|pair| pair[0].contact_index() == pair[1].contact_index())
        {
            return Err(invalid("source refresh indices must be distinct and in range").into());
        }
        let mut requests = ordered.iter().peekable();
        let mut replacements = Vec::with_capacity(ordered.len());
        let mut changes = Vec::with_capacity(ordered.len());
        let native = self.network.total_energy_j()?;
        let mut total_before = Dd::from_f64(native);
        let mut total_after = Dd::from_f64(native);
        let mut contact_before = Dd::ZERO;
        let mut contact_after = Dd::ZERO;
        let mut increments = Dd::ZERO;
        let mut parameter_work = Dd::ZERO;
        let mut tolerance = Dd::ZERO;
        for (index, point) in self.points.iter().enumerate() {
            poll(Some(gate))?;
            let approach = point.approach(self.network.components())?;
            let (before, after, negative) = if requests
                .peek()
                .is_some_and(|next| next.contact_index() == index)
            {
                let request = requests.next().unwrap();
                let (candidate, change) = match *request {
                    NormalSourceRefresh::Jkr { refreshed, .. } => prepare_jkr_contact_refresh(
                        &point.contact,
                        point.config,
                        approach,
                        refreshed,
                        gate,
                    )?,
                    NormalSourceRefresh::Dmt { refreshed, .. } => prepare_dmt_contact_refresh(
                        &point.contact,
                        point.config,
                        approach,
                        refreshed,
                        gate,
                    )?,
                    NormalSourceRefresh::Maugis { refreshed, .. } => {
                        prepare_maugis_contact_refresh(
                            &point.contact,
                            point.config,
                            approach,
                            refreshed,
                            gate,
                        )?
                    }
                };
                let negative =
                    negative_slope_root(candidate.source(), point.config, approach, approach)?;
                increments = increments
                    + (Dd::from_f64(change.contact_energy_after_j)
                        - Dd::from_f64(change.contact_energy_before_j));
                parameter_work = parameter_work + Dd::from_f64(change.parameter_work_j);
                tolerance = tolerance + Dd::from_f64(change.energy_tolerance_j);
                changes.push(IndexedSourceNormalThermalChange {
                    contact: index,
                    original_model_identity: point.contact.source().identity(),
                    model_identity: candidate.source().identity(),
                    change,
                });
                replacements.push((index, candidate));
                (
                    change.contact_energy_before_j,
                    change.contact_energy_after_j,
                    negative,
                )
            } else {
                let energy = observation(point.contact.source(), approach)?.0;
                (
                    energy,
                    energy,
                    negative_slope_root(point.contact.source(), point.config, approach, approach)?,
                )
            };
            self.negative_slopes[index] = negative;
            contact_before = contact_before + Dd::from_f64(before);
            contact_after = contact_after + Dd::from_f64(after);
            total_before = total_before + Dd::from_f64(before);
            total_after = total_after + Dd::from_f64(after);
        }
        // All successor sources enter the same coupled guard together. An
        // individually safe warm contact can destabilize an adhesive neighbour.
        check_coupled_stability(&self.compliance, &self.negative_slopes, gate)?;
        let residual = finite((increments - parameter_work).to_f64())?;
        let tolerance = finite(tolerance.to_f64())?;
        check_contact_work(residual, tolerance)?;
        let energy_after = finite(total_after.to_f64())?;
        limit(
            "refreshed simultaneous native/source energy magnitude",
            energy_after.abs(),
            self.network.config.maximum_total_energy_j,
        )?;
        let report = MultipleSourceNormalRefreshReport {
            sample: self.samples_rendered(),
            network_energy_j: native,
            contact_energy_before_j: finite(contact_before.to_f64())?,
            contact_energy_after_j: finite(contact_after.to_f64())?,
            energy_before_j: finite(total_before.to_f64())?,
            energy_after_j: energy_after,
            parameter_work_j: finite(parameter_work.to_f64())?,
            energy_residual_j: residual,
            energy_tolerance_j: tolerance,
            contacts: changes,
        };
        poll(Some(gate))?;
        Ok((replacements, report))
    }
}
