//! Atomic original normal/shear temperature changes at held native coordinates.
//!
//! Each selected contact consumes its own sealed held-indentation normal refresh
//! and matching original successor shear. The complete successor set passes the
//! same conservative smooth-branch curvature screen as construction. Every
//! original native mode, connection, material memory, radiation state, clock,
//! contact position/occupancy and last mechanical frame stays fixed.
//!
//! Signed parameter work is the explicit temperature driver's responsibility;
//! it is not frictional heat, an impulse or an implicit thermal solve. Each
//! source passes a contact-sized arithmetic gate before double-double aggregate
//! work admission. Unrelated native storage cannot hide a source-work mismatch.

use super::{
    CancelGate, Dd, DmtCoupledLoadLine, DmtCoupledRegime, DmtCoupledState,
    MultipleAdhesiveDmtModalError, MultipleAdhesiveDmtModalSystem, SourceAdhesiveShear, check_work,
    coupled_stability, finite, invalid, limit, poll, source_error,
};
use fs_blake3::ContentHash;
use fs_contact::normal_patch::adhesive::temperature::SourceDmtThermalRefresh;

/// One original source transition. Indices must be distinct and in range;
/// request order does not change evaluation, reporting or publication order.
#[derive(Clone, Copy)]
pub struct AdhesiveDmtSourceRefresh<'a> {
    pub contact: usize,
    pub refreshed: &'a SourceDmtThermalRefresh,
    pub next_shear: &'a SourceAdhesiveShear,
}

/// Independently admitted original-source work at one unchanged contact.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IndexedAdhesiveDmtThermalChange {
    pub contact: usize,
    pub original_model_identity: ContentHash,
    pub model_identity: ContentHash,
    pub position_m: [f64; 3],
    pub regime: DmtCoupledRegime,
    pub surface_energy_before_j: f64,
    pub surface_energy_after_j: f64,
    pub parameter_work_j: f64,
    pub energy_residual_j: f64,
    pub energy_tolerance_j: f64,
}

/// One accepted batch without native motion, loss or clock advancement.
#[derive(Clone, Debug, PartialEq)]
pub struct MultipleAdhesiveDmtSourceRefreshReport {
    pub sample: u64,
    /// Unchanged storage of the complete original native network.
    pub network_energy_j: f64,
    /// All contact potentials, including contacts outside the selected batch.
    pub surface_energy_before_j: f64,
    pub surface_energy_after_j: f64,
    pub energy_before_j: f64,
    pub energy_after_j: f64,
    /// Sum of the original signed parameter works, accumulated in double-double.
    pub parameter_work_j: f64,
    /// Sum of changed contact potential increments minus their signed works.
    /// This gate never subtracts the large total native energy observations.
    pub energy_residual_j: f64,
    /// Sum of independently enforced contact-sized arithmetic tolerances.
    pub energy_tolerance_j: f64,
    /// Selected contacts only, in ascending contact-index order.
    pub contacts: Vec<IndexedAdhesiveDmtThermalChange>,
}

type Replacements<'a> = Vec<(usize, DmtCoupledLoadLine<'a>, DmtCoupledState)>;

impl<'a> MultipleAdhesiveDmtModalSystem<'a> {
    /// Apply one nonempty bounded selected batch of sealed original normal and
    /// shear source changes. The temperature driver owns the returned signed
    /// parameter work exactly once; geometry and native histories do not move.
    /// Stale/foreign receipts, mismatched shear, coupled instability or any work
    /// refusal publish none of the batch. The complete native network is never
    /// cloned; cold source-transition storage is bounded by the contact count.
    pub fn apply_source_refreshes(
        &mut self,
        refreshes: &[AdhesiveDmtSourceRefresh<'a>],
        gate: &CancelGate,
    ) -> Result<MultipleAdhesiveDmtSourceRefreshReport, MultipleAdhesiveDmtModalError> {
        let (replacements, report) = self.prepare_source_refreshes(refreshes, gate)?;
        poll(Some(gate))?;
        self.publish_source_refreshes(replacements);
        Ok(report)
    }

    /// Stage another owner's receiving candidate against the complete report
    /// and indexed original replacement load lines/states, in ascending index
    /// order. Install that candidate infallibly only after this call succeeds.
    /// Receiver refusal and late cancellation leave every accepted source,
    /// contact branch, native history, clock and mechanical frame unchanged.
    pub fn apply_source_refreshes_with_admission<R, E>(
        &mut self,
        refreshes: &[AdhesiveDmtSourceRefresh<'a>],
        gate: &CancelGate,
        admit: impl FnOnce(
            &MultipleAdhesiveDmtSourceRefreshReport,
            &[(usize, DmtCoupledLoadLine<'a>, DmtCoupledState)],
        ) -> Result<R, E>,
    ) -> Result<R, E>
    where
        E: From<MultipleAdhesiveDmtModalError>,
    {
        let (replacements, report) = self
            .prepare_source_refreshes(refreshes, gate)
            .map_err(E::from)?;
        let received = admit(&report, &replacements)?;
        poll(Some(gate))
            .map_err(MultipleAdhesiveDmtModalError::from)
            .map_err(E::from)?;
        self.publish_source_refreshes(replacements);
        Ok(received)
    }

    fn publish_source_refreshes(&mut self, replacements: Replacements<'a>) {
        // All original sources, indices, coupled/work gates and receiver state
        // passed. No fallible operation or partial publication remains.
        for (index, load_line, state) in replacements {
            self.points[index].load_line = load_line;
            self.points[index].state = state;
        }
    }

    fn prepare_source_refreshes(
        &self,
        refreshes: &[AdhesiveDmtSourceRefresh<'a>],
        gate: &CancelGate,
    ) -> Result<
        (Replacements<'a>, MultipleAdhesiveDmtSourceRefreshReport),
        MultipleAdhesiveDmtModalError,
    > {
        poll(Some(gate))?;
        if refreshes.is_empty() || refreshes.len() > self.points.len() {
            return Err(invalid("nonempty bounded original DMT refresh batch required").into());
        }
        let mut ordered = refreshes.to_vec();
        ordered.sort_unstable_by_key(|r| r.contact);
        if ordered.last().unwrap().contact >= self.points.len()
            || ordered
                .windows(2)
                .any(|pair| pair[0].contact == pair[1].contact)
        {
            return Err(invalid("DMT refresh indices must be distinct and in range").into());
        }
        let mut requests = ordered.iter().peekable();
        let mut replacements = Vec::with_capacity(ordered.len());
        let mut changes = Vec::with_capacity(ordered.len());
        let native = self.network.total_energy_j()?;
        let mut total_before = Dd::from_f64(native);
        let mut total_after = Dd::from_f64(native);
        let mut surface_before = Dd::ZERO;
        let mut surface_after = Dd::ZERO;
        let mut increments = Dd::ZERO;
        let mut work = Dd::ZERO;
        let mut tolerance = Dd::ZERO;
        for (index, point) in self.points.iter().enumerate() {
            poll(Some(gate))?;
            let before = finite(
                point
                    .load_line
                    .normal()
                    .evaluate(point.state.position_m()[0])
                    .map_err(|e| source_error(index, e.into()))?
                    .reversible_energy_j,
            )?;
            let after = if requests.peek().is_some_and(|r| r.contact == index) {
                let request = requests.next().unwrap();
                let (next, state) = point
                    .load_line
                    .with_temperature_refresh(
                        point.state,
                        request.refreshed,
                        request.next_shear,
                        gate,
                    )
                    .map_err(|e| source_error(index, e))?;
                let after = finite(
                    next.normal()
                        .evaluate(state.position_m()[0])
                        .map_err(|e| source_error(index, e.into()))?
                        .reversible_energy_j,
                )?;
                let parameter = finite(request.refreshed.parameter_work_j())?;
                let change = Dd::from_f64(after) - Dd::from_f64(before);
                let residual = finite((change - Dd::from_f64(parameter)).to_f64())?;
                // The original held-geometry transition has no nonlinear motion
                // error. Match the single owner's pure arithmetic gate, without
                // adding mechanical tolerances or unrelated contact/native work.
                let allowance = finite(
                    128.0 * f64::EPSILON * before.abs().max(after.abs()).max(parameter.abs()),
                )?;
                check_work(Some(index), residual, allowance)?;
                check_work(
                    Some(index),
                    finite(
                        (change - Dd::from_f64(request.refreshed.potential_change_j())).to_f64(),
                    )?,
                    allowance,
                )?;
                increments = increments + change;
                work = work + Dd::from_f64(parameter);
                tolerance = tolerance + Dd::from_f64(allowance);
                changes.push(IndexedAdhesiveDmtThermalChange {
                    contact: index,
                    original_model_identity: point.load_line.normal().identity(),
                    model_identity: next.normal().identity(),
                    position_m: state.position_m(),
                    regime: state.regime(),
                    surface_energy_before_j: before,
                    surface_energy_after_j: after,
                    parameter_work_j: parameter,
                    energy_residual_j: residual,
                    energy_tolerance_j: allowance,
                });
                replacements.push((index, next, state));
                after
            } else {
                before
            };
            surface_before = surface_before + Dd::from_f64(before);
            surface_after = surface_after + Dd::from_f64(after);
            total_before = total_before + Dd::from_f64(before);
            total_after = total_after + Dd::from_f64(after);
        }
        // Local refreshed load lines alone cannot screen the shared response.
        // Include untouched neighbours and all replacements together without
        // cloning the native owner or modifying its accepted/scratch histories.
        let mut changed = replacements.iter().peekable();
        coupled_stability(
            self.points.iter().enumerate().map(|(index, point)| {
                let line = if changed.peek().is_some_and(|r| r.0 == index) {
                    &changed.next().unwrap().1
                } else {
                    &point.load_line
                };
                (line.normal(), line.shear())
            }),
            &self.compliance,
            gate,
        )?;
        let residual = finite((increments - work).to_f64())?;
        let allowance = finite(tolerance.to_f64())?;
        check_work(None, residual, allowance)?;
        let after = finite(total_after.to_f64())?;
        limit(
            "refreshed simultaneous native/DMT energy magnitude",
            after.abs(),
            self.network.config.maximum_total_energy_j,
        )?;
        let report = MultipleAdhesiveDmtSourceRefreshReport {
            sample: self.samples_rendered(),
            network_energy_j: native,
            surface_energy_before_j: finite(surface_before.to_f64())?,
            surface_energy_after_j: finite(surface_after.to_f64())?,
            energy_before_j: finite(total_before.to_f64())?,
            energy_after_j: after,
            parameter_work_j: finite(work.to_f64())?,
            energy_residual_j: residual,
            energy_tolerance_j: allowance,
            contacts: changes,
        };
        poll(Some(gate))?;
        Ok((replacements, report))
    }
}
