//! Transient porous absorption on the original acoustic and thermal owners.
//!
//! Each full DG element binds one original thermal prism. Original Darcy and
//! memory source squares give its separate viscous and thermal Joules over a
//! native acoustic step, including startup, transport and ringdown. No energy
//! deficit, boundary loss or numerical residual becomes a thermal source.
//!
//! Thermal substeps receive the exact interval's heat as average watts. Refine
//! the acoustic interval to resolve its thermal timing; refining thermal CFL
//! substeps alone does not resolve changes within this held-average source.
//! The material is frozen and coupling is one-way. Each cell's original storage
//! law must represent the intended fluid/skeleton thermal equilibrium; the
//! acoustic mechanism split does not infer separate constituent temperatures.

use super::*;
use fs_duct::dg::{
    DgDuct, DgMedium, DgPorousElementHeat, DgPorousHeating, DgPropagator, DgState, DgStepOutput,
    PortInput,
};
use fs_phs::linear_exact::DenseExponentialBudget;

// A large epoch must not silently change either native owner's physical step.
// Admit ordinary clock rounding only while retaining half of f64's significant
// bits in the interval. This is a time-resolution bound, not an error estimate
// for the acoustic or thermal model.
fn check_clock(native: f64, represented: f64) -> Result<(), PorousHeatingError> {
    if represented <= 0.0 || (represented - native).abs() > det::sqrt(f64::EPSILON) * native {
        Err(PorousHeatingError::Invalid(
            "thermal clock cannot resolve the original transient porous interval",
        ))
    } else {
        Ok(())
    }
}

/// Optional element-local source caches on the complete original DG operator.
/// Cold preparation shares one total dense byte/work allowance across all cells.
/// Warm advances reuse the caches and one source scratch buffer.
pub struct TransientPorousHeating<'d, M: DgMedium> {
    duct: &'d DgDuct<M>,
    propagator: &'d DgPropagator,
    elements: Vec<DgPorousElementHeat>,
    cells: Vec<PorousCellHeating>,
}

/// Physical acoustic heat and its actual native thermal credit for one cell.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransientPorousCellHeat {
    /// Original thermal destination, in DG element order.
    pub cell: CellId,
    /// Independently integrated original viscous and thermal source squares.
    pub generated: DgPorousHeating,
    /// Actual heat credited across accepted native thermal substeps [J].
    pub deposited_j: f64,
}

/// Both candidates publish together; inputs remain unchanged on any refusal.
pub struct TransientPorousHeatingStep<'t> {
    /// Unmodified original acoustic propagator's complete candidate state.
    pub acoustic: DgState,
    /// Original held-port work and acoustic ledger; never recast as cell heat.
    pub acoustic_step: DgStepOutput,
    /// Original finite-storage thermal candidate, retaining authored exchanges.
    pub thermal: DepositionSimulation<'t>,
    /// Native thermal audit, already including deposited porous heat once.
    pub heat_report: EnergyAudit,
    /// Physical source and deposition receipts in original DG element order.
    pub cells: Vec<TransientPorousCellHeat>,
    /// Sum of independent intrinsic porous sources over the complete specimen.
    pub generated: DgPorousHeating,
    /// Actual heat credited to all bound native cells [J].
    pub deposited_j: f64,
    /// Signed E_initial + port work - E_final - intrinsic porous heat [J].
    /// Includes any boundary/other loss AND numerical residual; it is not an
    /// independent loss measurement and is never deposited or clipped.
    pub other_acoustic_loss_and_residual_j: f64,
}

impl<'d, M: DgMedium> TransientPorousHeating<'d, M> {
    /// Bind exactly one active full-area thermal prism to every original DG
    /// element, in left-to-right order. Retain the original fluid/skeleton
    /// material law, all acoustic memories and the native held-input interval.
    /// The complete state must fit the existing 128-coordinate source owner;
    /// neither material poles nor spatial states are removed to fit a budget.
    pub fn try_new(
        duct: &'d DgDuct<M>,
        propagator: &'d DgPropagator,
        thermal: &DepositionSimulation<'_>,
        destinations: &[CellId],
        budget: DenseExponentialBudget,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, PorousHeatingError> {
        checkpoint()?;
        let count = duct.mesh().elements;
        if count == 0 || count > MAX_POROUS_HEATING_CELLS || destinations.len() != count {
            return Err(PorousHeatingError::Invalid(
                "transient porous heat needs one bounded destination per original DG element",
            ));
        }
        let element_budget = DenseExponentialBudget {
            max_dimension: budget.max_dimension,
            max_dense_bytes: budget.max_dense_bytes / count,
            max_multiply_terms: budget.max_multiply_terms / count as u64,
        };
        let mut elements = Vec::with_capacity(count);
        let mut cells = Vec::with_capacity(count);
        let mut seen = std::collections::BTreeSet::new();
        let mut axis = [0.0; 3];
        for (index, &destination) in destinations.iter().enumerate() {
            checkpoint()?;
            if !seen.insert(destination) {
                return Err(PorousHeatingError::Invalid(
                    "transient porous destinations must be distinct",
                ));
            }
            let mut refusal = None;
            let element =
                duct.prepare_porous_element_heat(propagator, index, element_budget, &mut || {
                    match checkpoint() {
                        Ok(()) => false,
                        Err(error) => {
                            refusal = Some(error);
                            true
                        }
                    }
                });
            if let Some(error) = refusal {
                return Err(error.into());
            }
            let element = element?;
            let [start, end] = element.interval_m();
            let cell = bind_cell(
                thermal,
                destination,
                start,
                end,
                element.area_m2(),
                cells.last(),
                &mut axis,
            )?;
            cells.push(cell);
            elements.push(element);
        }
        checkpoint()?;
        Ok(Self {
            duct,
            propagator,
            elements,
            cells,
        })
    }

    /// Original complete acoustic/material owner, without a replacement model.
    pub const fn duct(&self) -> &DgDuct<M> {
        self.duct
    }

    /// Advance one original acoustic interval and the same thermal interval as
    /// one transaction. `acoustic_time_s` is the caller's actual acoustic clock
    /// and must equal the native thermal clock; a step count is not a time.
    /// The held ports are identical for physical heat and native propagation.
    /// Whole and thermal-substep clocks must resolve their native intervals
    /// within `sqrt(f64::EPSILON)` relative error; use a local epoch otherwise.
    pub fn advance<'t>(
        &self,
        acoustic: &DgState,
        acoustic_time_s: f64,
        thermal: &DepositionSimulation<'t>,
        ports: &[PortInput],
        config: StepConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<TransientPorousHeatingStep<'t>, PorousHeatingError> {
        checkpoint()?;
        if !acoustic_time_s.is_finite() || acoustic_time_s != thermal.time_s() {
            return Err(PorousHeatingError::Invalid(
                "transient porous acoustic and thermal clocks must agree",
            ));
        }
        validate_cells(&self.cells, thermal, &mut checkpoint)?;
        let target = finite(thermal.time_s() + self.propagator.dt())?;
        let elapsed = target - thermal.time_s();
        check_clock(self.propagator.dt(), elapsed)?;
        let mut scratch = vec![0.0; self.elements[0].heat().scratch_len()];
        let mut generated = Vec::with_capacity(self.cells.len());
        let mut sources_w = BTreeMap::new();
        let (mut viscous, mut relaxation) = (Dd::ZERO, Dd::ZERO);
        for (element, cell) in self.elements.iter().zip(&self.cells) {
            let mut refusal = None;
            let heat =
                element
                    .heat()
                    .integrated_heat_into_j(
                        acoustic,
                        ports,
                        &mut scratch,
                        &mut || match checkpoint() {
                            Ok(()) => false,
                            Err(error) => {
                                refusal = Some(error);
                                true
                            }
                        },
                    );
            if let Some(error) = refusal {
                return Err(error.into());
            }
            let heat = heat?;
            sources_w.insert(cell.cell, finite(heat.total_j() / elapsed)?);
            viscous = viscous + Dd::from_f64(heat.viscous_j);
            relaxation = relaxation + Dd::from_f64(heat.thermal_j);
            generated.push(heat);
        }
        let mut acoustic_candidate = acoustic.clone();
        let acoustic_step = self.propagator.step(&mut acoustic_candidate, ports)?;
        checkpoint()?;
        let observer = IntervalHeat {
            cells: &self.cells,
            mean_sources_w: sources_w.values().copied().collect(),
            sources_w,
            deposited: vec![Dd::ZERO; self.cells.len()],
            steps: 0,
        };
        let (thermal, observer, heat_report) =
            thermal.advance_observed(&observer, target, config, &mut checkpoint)?;
        validate_cells(&self.cells, &thermal, &mut checkpoint)?;
        let mut deposited = Dd::ZERO;
        let mut cells = Vec::with_capacity(self.cells.len());
        for ((cell, generated), actual) in self.cells.iter().zip(generated).zip(observer.deposited)
        {
            checkpoint()?;
            let actual = finite(actual.to_f64())?;
            let expected = generated.total_j();
            let tolerance = finite(
                512.0
                    * f64::EPSILON
                    * (observer.steps + 1) as f64
                    * actual.abs().max(expected.abs()),
            )?;
            if (actual - expected).abs() > tolerance {
                return Err(PorousHeatingError::Invalid(
                    "native thermal interval did not receive its original porous source integral",
                ));
            }
            deposited = deposited + Dd::from_f64(actual);
            cells.push(TransientPorousCellHeat {
                cell: cell.cell,
                generated,
                deposited_j: actual,
            });
        }
        let generated = DgPorousHeating {
            viscous_j: finite(viscous.to_f64())?,
            thermal_j: finite(relaxation.to_f64())?,
        };
        let other_acoustic_loss_and_residual_j = finite(
            acoustic.stored_energy_j() + acoustic_step.ledger.port_work_step_j
                - acoustic_candidate.stored_energy_j()
                - generated.total_j(),
        )?;
        checkpoint()?;
        Ok(TransientPorousHeatingStep {
            acoustic: acoustic_candidate,
            acoustic_step,
            thermal,
            heat_report,
            cells,
            generated,
            deposited_j: finite(deposited.to_f64())?,
            other_acoustic_loss_and_residual_j,
        })
    }
}

#[derive(Clone)]
struct IntervalHeat<'c> {
    cells: &'c [PorousCellHeating],
    // Same stable sorted-cell order as sources_w, never updated by a substep.
    mean_sources_w: Vec<f64>,
    sources_w: BTreeMap<CellId, f64>,
    deposited: Vec<Dd>,
    steps: usize,
}
impl DepositionObserver for IntervalHeat<'_> {
    type Error = PorousHeatingError;
    fn requires_step_interval(&self) -> bool {
        true
    }
    fn step_sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        start_s: f64,
        dt_s: f64,
        end_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        validate_cells(self.cells, thermal, checkpoint)?;
        let represented = end_s - start_s;
        check_clock(dt_s, represented)?;
        // The native thermal owner integrates dt_s even when its rounded clock
        // difference differs slightly. Credit this represented share of the
        // original Q exactly once; contiguous shares telescope to the full Q.
        let clock_ratio = represented / dt_s;
        for ((_, watts), &mean) in self.sources_w.iter_mut().zip(&self.mean_sources_w) {
            checkpoint()?;
            *watts = finite(mean * clock_ratio)?;
        }
        Ok(self.sources_w.clone())
    }
    fn event(&mut self, _: &Event, thermal: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        validate_cells(self.cells, thermal, &mut || Ok(()))
    }
    fn step(
        &mut self,
        _: f64,
        dt: f64,
        _: &BTreeMap<CellId, f64>,
        _: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        for (cell, deposited) in self.cells.iter().zip(&mut self.deposited) {
            checkpoint()?;
            *deposited = *deposited + Dd::from_f64(finite(self.sources_w[&cell.cell] * dt)?);
        }
        self.steps += 1;
        Ok(())
    }
}
