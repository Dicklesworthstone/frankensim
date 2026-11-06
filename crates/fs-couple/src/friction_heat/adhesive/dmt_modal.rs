//! Original simultaneous DMT junction loss and shared finite-cell heat.
//!
//! The native modal/contact owner is solved once. Only its actual junction
//! dissipation enters the existing finite-volume heat solve; normal potential,
//! prescribed work and network damping keep their own mechanical accounts.
//! Every heat trial retains that same native candidate. Both owners publish
//! together only after the physical work gates and downstream admission pass.
//!
//! Mechanical sources and the incoming conducting footprints are frozen over
//! this split. True DMT area still evolves in the mechanical law; its fractional
//! change must fit an explicit allowance. Contact discovery, opening, moving
//! coverage and thermodynamic surface-entropy feedback are not inferred here.

use crate::friction_heat::FrictionHeatRoute;
use crate::friction_heat::cohesive::{CohesiveHeatError, solve_frozen_cohesive_heat};
use crate::friction_heat::network::{FrictionHeatNetworkReport, FrictionHeatNetworkStep};
use crate::render::schedule::force::coupled::contact::adhesive_dmt::multiple::{
    MultipleAdhesiveDmtModalError, MultipleAdhesiveDmtModalFrame, MultipleAdhesiveDmtModalSystem,
};
use fs_exec::CancelGate;
use fs_math::dd::Dd;
use fs_solver::NewtonReport;
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};

/// Reuse the original finite-cell transport inputs and nonlinear controls.
pub use crate::friction_heat::cohesive::{
    CohesiveHeatConductivity as DmtModalHeatConductivity,
    CohesiveHeatControls as DmtModalThermalControls, CohesiveHeatSource as DmtModalHeatSource,
};

/// One physical junction and its original ordered thermal face.
/// The caller declares the source A/B surfaces to be these A/B cells; modal
/// work columns alone cannot reconstruct their world placement or material IDs.
#[derive(Clone, Copy)]
pub struct DmtModalHeatContact<'a> {
    pub contact_index: usize,
    pub route: &'a FrictionHeatRoute,
    /// Positive numerical seed, never an imposed interface temperature.
    pub temperature_guess_k: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct DmtModalHeatControls {
    pub thermal: DmtModalThermalControls,
    /// Bound |A_end/A_start - 1| in [0,1); the affine DMT approach makes the
    /// endpoint areas the path extrema. No numerical area is substituted into
    /// the original mechanics. The thermal footprint remains A_start.
    pub maximum_fractional_area_change: f64,
    /// Bound the relative difference between the representable thermal interval
    /// and the native fixed sample period, in [0,1). Long clocks must not erase
    /// or materially stretch a physical sample while keeping its original loss.
    pub maximum_relative_clock_error: f64,
}

#[derive(Debug)]
pub enum DmtModalHeatError<E> {
    Invalid(&'static str),
    Mechanical(MultipleAdhesiveDmtModalError),
    Thermal(DepositionError),
    Transport(CohesiveHeatError<E>),
    AreaChange {
        contact: usize,
        fractional_change: f64,
    },
    Energy {
        residual_j: f64,
        absolute_stage_residual_j: f64,
    },
}
impl<E> From<MultipleAdhesiveDmtModalError> for DmtModalHeatError<E> {
    fn from(error: MultipleAdhesiveDmtModalError) -> Self {
        Self::Mechanical(error)
    }
}
impl<E> From<DepositionError> for DmtModalHeatError<E> {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl<E> From<CohesiveHeatError<E>> for DmtModalHeatError<E> {
    fn from(error: CohesiveHeatError<E>) -> Self {
        Self::Transport(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for DmtModalHeatError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native simultaneous DMT heat refused: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for DmtModalHeatError<E> {}

/// Mechanical construction order is independent of the canonical thermal order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DmtModalHeatFootprint {
    pub contact_index: usize,
    pub link_id: u64,
    pub incoming_area_m2: f64,
    pub endpoint_area_m2: f64,
    pub mean_real_area_m2: f64,
    pub fractional_area_change: f64,
}

pub struct HeatedDmtModalFrame<P> {
    pub mechanical: MultipleAdhesiveDmtModalFrame,
    pub heat: FrictionHeatNetworkReport,
    /// Same ascending link order as heat.interfaces, preserving contact indices.
    pub footprints: Vec<DmtModalHeatFootprint>,
    pub source: P,
    pub nonlinear: NewtonReport,
    pub thermal_trials: usize,
    pub thermal_duration_s: f64,
    /// Mechanical + caloric balance with friction eliminated as internal heat.
    /// Native network damping remains an explicitly exported mechanical loss.
    pub combined_energy_residual_j: f64,
    /// Absolute native and per-interface defects cannot cancel into acceptance.
    pub absolute_stage_energy_residual_j: f64,
}

impl MultipleAdhesiveDmtModalSystem<'_> {
    /// Advance all original junctions and all their shared thermal cells once.
    ///
    /// Every mechanical contact needs one distinct original thermal route.
    /// Bindings may be permuted; the source callback and thermal report use
    /// ascending link IDs. Shared cells are retained in one thermal solve.
    /// Conductivity comes from each original body's own proposed thermal state;
    /// the solved common face temperatures determine the signed heat exchange.
    /// The callback supplies OTHER heat only and must stage external history.
    ///
    /// Every straight approach path must remain in positive-area solid contact.
    /// The route must match its incoming true area. Rebind physical coverage
    /// before the next step when that area changes; do not reuse an obsolete
    /// face merely by relaxing the area-change control. Source-temperature
    /// refresh is a separate explicit boundary operation with parameter work.
    /// This is frozen-source/frozen-footprint splitting plus backward Euler
    /// transport, requiring temporal and spatial refinement, not strong feedback.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_friction_heat<'m, P, E>(
        &mut self,
        thermal: &mut DepositionSimulation<'m>,
        bindings: &[DmtModalHeatContact<'_>],
        external: &[f64],
        prescribed_rates_m_s: &[[f64; 3]],
        controls: DmtModalHeatControls,
        energy_guess_j: &[f64],
        resolve: impl Fn(&DepositionSimulation<'m>, &[u64], &[f64]) -> Result<DmtModalHeatSource<P>, E>,
        gate: &CancelGate,
    ) -> Result<HeatedDmtModalFrame<P>, DmtModalHeatError<E>> {
        self.step_with_friction_heat_admission(
            thermal,
            bindings,
            external,
            prescribed_rates_m_s,
            controls,
            energy_guess_j,
            resolve,
            gate,
            |_, _| Ok(()),
        )
        .map(|(report, ())| report)
    }

    /// Stage an additional receiver against both real candidate endpoints.
    /// Its returned value must be infallibly installable. A refusal or late
    /// cancellation preserves native modes, source history, last frame and heat.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn step_with_friction_heat_admission<'m, P, E, R>(
        &mut self,
        thermal: &mut DepositionSimulation<'m>,
        bindings: &[DmtModalHeatContact<'_>],
        external: &[f64],
        prescribed_rates_m_s: &[[f64; 3]],
        controls: DmtModalHeatControls,
        energy_guess_j: &[f64],
        resolve: impl Fn(&DepositionSimulation<'m>, &[u64], &[f64]) -> Result<DmtModalHeatSource<P>, E>,
        gate: &CancelGate,
        admit: impl FnOnce(
            &HeatedDmtModalFrame<P>,
            &DepositionSimulation<'m>,
        ) -> Result<R, DmtModalHeatError<E>>,
    ) -> Result<(HeatedDmtModalFrame<P>, R), DmtModalHeatError<E>> {
        let check = || {
            if gate.is_requested() {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        };
        check()?;
        let count = self.contact_count();
        let sample = self.samples_rendered();
        let next_sample = sample
            .checked_add(1)
            .ok_or(DmtModalHeatError::Invalid("native sample clock exhausted"))?;
        let dt = self.sample_period_s();
        let start = sample as f64 * dt;
        let target = next_sample as f64 * dt;
        let duration = target - start;
        if sample >= (1_u64 << 53)
            || thermal.time_s() != start
            || !duration.is_finite()
            || duration <= 0.0
            || !controls.maximum_relative_clock_error.is_finite()
            || !(0.0..1.0).contains(&controls.maximum_relative_clock_error)
            || ((duration - dt) / dt).abs() > controls.maximum_relative_clock_error
            || !controls.maximum_fractional_area_change.is_finite()
            || !(0.0..1.0).contains(&controls.maximum_fractional_area_change)
            || bindings.len() != count
            || count == 0
            || count > controls.thermal.maximum_interfaces
        {
            return Err(DmtModalHeatError::Invalid(
                "synchronized bounded contacts, sample clock and frozen-area allowances required",
            ));
        }
        let mut seen = vec![false; count];
        let mut ordered: Vec<_> = bindings.iter().collect();
        ordered.sort_by_key(|binding| binding.route.link_id());
        let mut footprints = Vec::with_capacity(count);
        let mut guesses = Vec::with_capacity(count);
        for (i, binding) in ordered.iter().enumerate() {
            check()?;
            let c = binding.contact_index;
            if c >= count
                || seen[c]
                || (i > 0 && ordered[i - 1].route.link_id() == binding.route.link_id())
                || !binding.temperature_guess_k.is_finite()
                || binding.temperature_guess_k <= 0.0
            {
                return Err(DmtModalHeatError::Invalid(
                    "one positive-temperature original heat route per contact required",
                ));
            }
            seen[c] = true;
            let response = self
                .normal(c)
                .unwrap()
                .evaluate(self.state(c).unwrap().position_m()[0])
                .map_err(|error| MultipleAdhesiveDmtModalError::Contact {
                    contact: c,
                    source: error.into(),
                })?;
            let area = response.contact_area_m2;
            let actual_face = binding.route.contact_area_m2();
            if area <= 0.0
                || (area - actual_face).abs() > 128.0 * f64::EPSILON * area.max(actual_face)
            {
                return Err(DmtModalHeatError::Invalid(
                    "thermal footprint must equal the incoming true attached DMT area",
                ));
            }
            footprints.push(DmtModalHeatFootprint {
                contact_index: c,
                link_id: binding.route.link_id(),
                incoming_area_m2: area,
                endpoint_area_m2: area,
                mean_real_area_m2: area,
                fractional_area_change: 0.0,
            });
            guesses.push(binding.temperature_guess_k);
        }
        let routes: Vec<_> = ordered.iter().map(|binding| binding.route).collect();
        let step = FrictionHeatNetworkStep::new(
            thermal,
            target,
            &routes,
            controls.thermal.maximum_interfaces,
            check,
        )?;
        let (next_thermal, report, admitted) =
            self.step_with_admission(external, prescribed_rates_m_s, gate, |mechanical, _| {
                let mut losses = Vec::with_capacity(count);
                for footprint in &mut footprints {
                    check()?;
                    let contact = &mechanical.contacts[footprint.contact_index];
                    let endpoint = contact.path.contact.contact_area_m2;
                    let change = ((endpoint - footprint.incoming_area_m2)
                        / footprint.incoming_area_m2)
                        .abs();
                    if endpoint <= 0.0
                        || !change.is_finite()
                        || change > controls.maximum_fractional_area_change
                    {
                        return Err(DmtModalHeatError::AreaChange {
                            contact: footprint.contact_index,
                            fractional_change: change,
                        });
                    }
                    footprint.endpoint_area_m2 = endpoint;
                    footprint.mean_real_area_m2 = contact.path.mean_real_area_m2;
                    footprint.fractional_area_change = change;
                    losses.push(contact.friction_dissipation_j);
                }
                let solved = solve_frozen_cohesive_heat(
                    thermal,
                    &step,
                    &guesses,
                    &losses,
                    controls.thermal,
                    energy_guess_j,
                    resolve,
                    gate,
                )?;
                let transfer = Dd::from_f64(solved.heat.dissipated_j)
                    - Dd::from_f64(mechanical.friction_dissipation_j);
                let residual = (Dd::from_f64(mechanical.energy_residual_j)
                    + Dd::from_f64(solved.heat.total_energy_residual_j)
                    + transfer)
                    .to_f64();
                let mut absolute = Dd::from_f64(mechanical.energy_residual_j.abs())
                    + Dd::from_f64(solved.heat.thermal.residual_j.abs())
                    + Dd::from_f64(transfer.to_f64().abs());
                for interface in &solved.heat.interfaces {
                    check()?;
                    absolute = absolute + Dd::from_f64(interface.interface_residual_j.abs());
                }
                let absolute = absolute.to_f64();
                if !residual.is_finite()
                    || !absolute.is_finite()
                    || residual.abs() > controls.thermal.combined_energy_tolerance_j
                    || absolute > controls.thermal.combined_energy_tolerance_j
                {
                    return Err(DmtModalHeatError::Energy {
                        residual_j: residual,
                        absolute_stage_residual_j: absolute,
                    });
                }
                let report = HeatedDmtModalFrame {
                    mechanical: mechanical.clone(),
                    heat: solved.heat,
                    footprints,
                    source: solved.source,
                    nonlinear: solved.nonlinear,
                    thermal_trials: solved.trials,
                    thermal_duration_s: duration,
                    combined_energy_residual_j: residual,
                    absolute_stage_energy_residual_j: absolute,
                };
                let admitted = admit(&report, &solved.thermal)?;
                // The original native transaction checks the gate once more
                // before publication. There is no fallible work after it returns.
                Ok((solved.thermal, report, admitted))
            })?;
        *thermal = next_thermal;
        Ok((report, admitted))
    }
}
