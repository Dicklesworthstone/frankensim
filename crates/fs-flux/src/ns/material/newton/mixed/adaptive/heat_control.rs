//! Optional physical heat and signed temporal-kinetic-defect control.

use super::{
    MaterialNavierStokesAdaptiveAdvance, MaterialNavierStokesSdirk2Step,
    MaterialSdirk2AdaptiveError, poll,
};
use fs_exec::CancelGate;
use fs_math::dd::Dd;

/// Absolute physical heat allowance [J/m] and relative local allowance.
/// Relative scales use this comparison's generated heat, never cumulative heat.
#[derive(Clone, Copy, Debug)]
pub struct MaterialSdirk2HeatTolerance {
    pub absolute_j_per_m: f64,
    pub relative: f64,
}
impl MaterialSdirk2HeatTolerance {
    fn valid(self) -> bool {
        self.absolute_j_per_m.is_finite()
            && self.absolute_j_per_m >= 0.0
            && self.relative.is_finite()
            && self.relative >= 0.0
    }
    fn ratio(self, error: f64, scale: f64) -> Result<f64, MaterialSdirk2AdaptiveError> {
        Ok(ratio(
            error,
            finite(self.absolute_j_per_m + self.relative * scale)?,
        ))
    }
}

/// Heat-aware acceptance on the SAME original SDIRK2 trial solves.
///
/// Cell/edge indicators compare one coarse integral with the sum of two fine
/// integrals, divided by three for the second-order local-error estimate. The
/// total indicators sum absolute spatial discrepancies, preventing cancellation.
///
/// The separate time-defect allowance is A + r*H for the whole successful call,
/// where H is accepted original bulk PLUS slip heat and the charged cost is the
/// sum of absolute signed kinetic defects of every accepted fine half step.
/// Each comparison receives A*dt/requested_duration + r*H_comparison. This
/// explicit allocation prevents an early step spending later absolute credit;
/// it can cost extra steps during a short stiff transient in a long interval.
/// Unused local credit is not borrowed by later comparisons.
///
/// This controls temporal resolution, not mesh error, pressure/stress error,
/// caloric evolution or full energy balance. SDIRK2's signed time defect is
/// not nonnegative numerical damping and is never deposited as heat.
#[derive(Clone, Copy, Debug)]
pub struct MaterialSdirk2HeatControls {
    pub cell_bulk: MaterialSdirk2HeatTolerance,
    pub boundary_slip: MaterialSdirk2HeatTolerance,
    pub total_bulk: MaterialSdirk2HeatTolerance,
    pub total_slip: MaterialSdirk2HeatTolerance,
    pub absolute_temporal_defect_budget_j_per_m: f64,
    /// Finite in [0,1); refers to accepted physical bulk plus slip heat.
    pub relative_temporal_defect_budget: f64,
}
impl Default for MaterialSdirk2HeatControls {
    fn default() -> Self {
        let tolerance = MaterialSdirk2HeatTolerance {
            absolute_j_per_m: 1e-8,
            relative: 1e-3,
        };
        Self {
            cell_bulk: tolerance,
            boundary_slip: tolerance,
            total_bulk: tolerance,
            total_slip: tolerance,
            absolute_temporal_defect_budget_j_per_m: 1e-8,
            relative_temporal_defect_budget: 1e-3,
        }
    }
}

/// Accepted-comparison maxima and the complete accepted temporal energy account.
#[derive(Clone, Copy, Debug, Default)]
pub struct MaterialSdirk2HeatReport {
    pub maximum_cell_bulk_error_j_per_m: f64,
    pub maximum_boundary_slip_error_j_per_m: f64,
    pub maximum_bulk_error_l1_j_per_m: f64,
    pub maximum_slip_error_l1_j_per_m: f64,
    pub maximum_local_heat_error_ratio: f64,
    pub maximum_temporal_defect_ratio: f64,
    /// Signed sum, which can cancel and therefore does not define the budget.
    pub signed_temporal_kinetic_defect_j_per_m: f64,
    /// Sum of |D| for EACH accepted fine half step, with no sign cancellation.
    pub absolute_temporal_kinetic_defect_j_per_m: f64,
    pub temporal_defect_allowance_j_per_m: f64,
}

/// Whole-call flow/heat result plus opt-in temporal-accuracy diagnostics.
pub struct MaterialNavierStokesHeatAdaptiveAdvance<G> {
    pub flow: MaterialNavierStokesAdaptiveAdvance<G>,
    pub heat_control: MaterialSdirk2HeatReport,
}

pub(super) struct HeatControl {
    controls: MaterialSdirk2HeatControls,
    duration_s: f64,
    signed: Dd,
    absolute: Dd,
    heat: Dd,
    report: MaterialSdirk2HeatReport,
}

pub(super) struct HeatTrial {
    signed: Dd,
    absolute: Dd,
    heat: Dd,
    cell: f64,
    edge: f64,
    bulk_l1: f64,
    slip_l1: f64,
    local_ratio: f64,
    defect_ratio: f64,
}
impl HeatTrial {
    pub(super) fn error_ratio(&self) -> f64 {
        self.local_ratio.max(self.defect_ratio)
    }
}
impl HeatControl {
    pub(super) fn new(
        controls: MaterialSdirk2HeatControls,
        duration_s: f64,
    ) -> Result<Self, MaterialSdirk2AdaptiveError> {
        if ![
            controls.cell_bulk,
            controls.boundary_slip,
            controls.total_bulk,
            controls.total_slip,
        ]
        .iter()
        .all(|x| x.valid())
            || !controls.absolute_temporal_defect_budget_j_per_m.is_finite()
            || controls.absolute_temporal_defect_budget_j_per_m < 0.0
            || !controls.relative_temporal_defect_budget.is_finite()
            || !(0.0..1.0).contains(&controls.relative_temporal_defect_budget)
            || !duration_s.is_finite()
            || duration_s <= 0.0
        {
            return Err(invalid(
                "finite heat tolerances and temporal-defect budgets required",
            ));
        }
        Ok(Self {
            controls,
            duration_s,
            signed: Dd::ZERO,
            absolute: Dd::ZERO,
            heat: Dd::ZERO,
            report: MaterialSdirk2HeatReport::default(),
        })
    }

    pub(super) fn trial<G>(
        &self,
        coarse: &MaterialNavierStokesSdirk2Step<G>,
        first: &MaterialNavierStokesSdirk2Step<G>,
        second: &MaterialNavierStokesSdirk2Step<G>,
        dt: f64,
        gate: &CancelGate,
    ) -> Result<HeatTrial, MaterialSdirk2AdaptiveError> {
        poll(gate)?;
        let c = &coarse.dissipated_heat;
        let a = &first.dissipated_heat;
        let b = &second.dissipated_heat;
        let cells = c.cell_bulk_energy_j_per_m();
        if cells.len() != a.cell_bulk_energy_j_per_m().len()
            || cells.len() != b.cell_bulk_energy_j_per_m().len()
            || c.boundary_slip_energy().len() != a.boundary_slip_energy().len()
            || c.boundary_slip_energy().len() != b.boundary_slip_energy().len()
        {
            return Err(invalid("original heat comparison topology changed"));
        }
        let mut bulk_errors = SpatialHeatError::default();
        let mut slip_errors = SpatialHeatError::default();
        for (i, ((&coarse, &first), &second)) in cells
            .iter()
            .zip(a.cell_bulk_energy_j_per_m())
            .zip(b.cell_bulk_energy_j_per_m())
            .enumerate()
        {
            if i % 256 == 0 {
                poll(gate)?;
            }
            bulk_errors.add(coarse, first, second, self.controls.cell_bulk)?;
        }
        for (i, ((coarse, first), second)) in c
            .boundary_slip_energy()
            .iter()
            .zip(a.boundary_slip_energy())
            .zip(b.boundary_slip_energy())
            .enumerate()
        {
            if i % 256 == 0 {
                poll(gate)?;
            }
            let identity = (coarse.edge, coarse.component);
            if identity != (first.edge, first.component)
                || identity != (second.edge, second.component)
            {
                return Err(invalid(
                    "original heat comparison boundary identity changed",
                ));
            }
            slip_errors.add(
                coarse.energy_j_per_m,
                first.energy_j_per_m,
                second.energy_j_per_m,
                self.controls.boundary_slip,
            )?;
        }
        let bulk_l1 = finite(bulk_errors.l1.to_f64())?;
        let slip_l1 = finite(slip_errors.l1.to_f64())?;
        let cell = bulk_errors.maximum;
        let edge = slip_errors.maximum;
        let (_, bulk_scale) = heat_difference(
            c.bulk_energy_j_per_m(),
            a.bulk_energy_j_per_m(),
            b.bulk_energy_j_per_m(),
        )?;
        let (_, slip_scale) = heat_difference(
            c.slip_energy_j_per_m(),
            a.slip_energy_j_per_m(),
            b.slip_energy_j_per_m(),
        )?;
        let local_ratio = bulk_errors
            .ratio
            .max(slip_errors.ratio)
            .max(self.controls.total_bulk.ratio(bulk_l1, bulk_scale)?)
            .max(self.controls.total_slip.ratio(slip_l1, slip_scale)?);
        let d1 = first.temporal_kinetic_defect_j_per_m;
        let d2 = second.temporal_kinetic_defect_j_per_m;
        let signed = Dd::from_f64(finite(d1)?) + Dd::from_f64(finite(d2)?);
        let absolute = Dd::from_f64(d1.abs()) + Dd::from_f64(d2.abs());
        let heat = Dd::from_f64(a.bulk_energy_j_per_m())
            + Dd::from_f64(a.slip_energy_j_per_m())
            + Dd::from_f64(b.bulk_energy_j_per_m())
            + Dd::from_f64(b.slip_energy_j_per_m());
        let allowance = Dd::from_f64(self.controls.absolute_temporal_defect_budget_j_per_m)
            * (Dd::from_f64(dt) / Dd::from_f64(self.duration_s))
            + Dd::from_f64(self.controls.relative_temporal_defect_budget) * heat;
        let mut defect_ratio = ratio(finite(absolute.to_f64())?, finite(allowance.to_f64())?);
        // An inequality can fail below one f64 ulp. Preserve that refusal even
        // when the reported quotient rounds to exactly one.
        if allowance.lt(absolute) && defect_ratio <= 1.0 {
            defect_ratio = above_one();
        }
        let total_absolute = self.absolute + absolute;
        let total_allowance = Dd::from_f64(self.controls.absolute_temporal_defect_budget_j_per_m)
            + Dd::from_f64(self.controls.relative_temporal_defect_budget) * (self.heat + heat);
        if total_allowance.lt(total_absolute) && defect_ratio <= 1.0 {
            defect_ratio = above_one();
        }
        Ok(HeatTrial {
            signed,
            absolute,
            heat,
            cell,
            edge,
            bulk_l1,
            slip_l1,
            local_ratio,
            defect_ratio,
        })
    }

    pub(super) fn accept(&mut self, trial: HeatTrial) -> Result<(), MaterialSdirk2AdaptiveError> {
        self.signed = self.signed + trial.signed;
        self.absolute = self.absolute + trial.absolute;
        self.heat = self.heat + trial.heat;
        finite(self.signed.to_f64())?;
        finite(self.absolute.to_f64())?;
        finite(self.heat.to_f64())?;
        self.report.maximum_cell_bulk_error_j_per_m =
            self.report.maximum_cell_bulk_error_j_per_m.max(trial.cell);
        self.report.maximum_boundary_slip_error_j_per_m = self
            .report
            .maximum_boundary_slip_error_j_per_m
            .max(trial.edge);
        self.report.maximum_bulk_error_l1_j_per_m =
            self.report.maximum_bulk_error_l1_j_per_m.max(trial.bulk_l1);
        self.report.maximum_slip_error_l1_j_per_m =
            self.report.maximum_slip_error_l1_j_per_m.max(trial.slip_l1);
        self.report.maximum_local_heat_error_ratio = self
            .report
            .maximum_local_heat_error_ratio
            .max(trial.local_ratio);
        self.report.maximum_temporal_defect_ratio = self
            .report
            .maximum_temporal_defect_ratio
            .max(trial.defect_ratio);
        Ok(())
    }

    pub(super) fn finish(
        mut self,
    ) -> Result<MaterialSdirk2HeatReport, MaterialSdirk2AdaptiveError> {
        self.report.signed_temporal_kinetic_defect_j_per_m = finite(self.signed.to_f64())?;
        self.report.absolute_temporal_kinetic_defect_j_per_m = finite(self.absolute.to_f64())?;
        let allowance = Dd::from_f64(self.controls.absolute_temporal_defect_budget_j_per_m)
            + Dd::from_f64(self.controls.relative_temporal_defect_budget) * self.heat;
        self.report.temporal_defect_allowance_j_per_m = finite(allowance.to_f64())?;
        if allowance.lt(self.absolute)
            || self.report.absolute_temporal_kinetic_defect_j_per_m
                > self.report.temporal_defect_allowance_j_per_m
        {
            return Err(invalid("accepted temporal-defect budget exceeded"));
        }
        Ok(self.report)
    }
}

#[derive(Default)]
struct SpatialHeatError {
    maximum: f64,
    l1: Dd,
    ratio: f64,
}
impl SpatialHeatError {
    fn add(
        &mut self,
        coarse: f64,
        first: f64,
        second: f64,
        tolerance: MaterialSdirk2HeatTolerance,
    ) -> Result<(), MaterialSdirk2AdaptiveError> {
        let (error, scale) = heat_difference(coarse, first, second)?;
        self.maximum = self.maximum.max(error);
        self.l1 = self.l1 + Dd::from_f64(error);
        self.ratio = self.ratio.max(tolerance.ratio(error, scale)?);
        Ok(())
    }
}

fn heat_difference(
    coarse: f64,
    first: f64,
    second: f64,
) -> Result<(f64, f64), MaterialSdirk2AdaptiveError> {
    let fine = Dd::from_f64(first) + Dd::from_f64(second);
    let error = finite(((fine - Dd::from_f64(coarse)) / Dd::from_f64(3.0)).to_f64())?.abs();
    let scale = finite(fine.to_f64())?.max(coarse);
    Ok((error, scale))
}
fn ratio(cost: f64, allowance: f64) -> f64 {
    if cost == 0.0 {
        0.0
    } else if allowance == 0.0 {
        f64::INFINITY
    } else {
        cost / allowance
    }
}
fn above_one() -> f64 {
    f64::from_bits(1.0_f64.to_bits() + 1)
}
fn finite(value: f64) -> Result<f64, MaterialSdirk2AdaptiveError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("unrepresentable heat-control account"))
    }
}
fn invalid(message: &'static str) -> MaterialSdirk2AdaptiveError {
    MaterialSdirk2AdaptiveError::Invalid(message)
}

#[cfg(test)]
#[path = "heat_control/tests.rs"]
mod tests;
