//! Adaptive molecular exchange with the original finite thermal bodies.
//!
//! Every trial repeats the native thermal-then-molecular transaction. The two
//! accepted half steps alone update gas composition, momentum and phase storage.
//! Local estimates include each gross molecular beam, its signed net, each
//! body's heat transfer and latent storage. They are numerical estimates, not
//! bounds on source uncertainty or on the omitted spatial/container dynamics.
//!
//! Body and prescribed gas heats describe the WHOLE requested window and are
//! applied as held powers. Radiation and state-resolved convection retain their
//! original native laws. Fixed window accounts use original-to-final gas U+K
//! and body storage, with one original source datum per conserved gas species.

use super::super::{
    MolecularNetworkError, MolecularNetworkLink, MovingGasBoundarySolve, MovingGasParcel,
    trajectory::{
        self as molecular, MolecularNetworkAdaptiveControls, MolecularNetworkLinkTotals,
        MolecularNetworkTrajectoryError, MolecularNetworkWindowBalance,
        accounts::{Sources, Totals, stream_error_from_states, window_balance},
    },
};
use super::{
    MAX_NETWORK_THERMAL_BODIES, MolecularGasThermalBoundary, MolecularNetworkThermalCandidate,
    MolecularNetworkThermalControls, MolecularNetworkThermalError, MolecularNetworkThermalStep,
    step_candidate_with_convection,
};
use crate::finite_gas_thermal::{
    ConvectionCoefficient,
    bodies::{FiniteGasBody, MAX_BODIES},
};
use fs_exec::Cx;
use fs_material::{
    phase::EquilibriumPhaseState, thermochemical_gas::ResolvedThermochemicalGasState,
};
use fs_math::dd::Dd;
use std::{collections::BTreeMap, convert::Infallible};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolecularNetworkThermalAdaptiveControls {
    pub molecular: MolecularNetworkAdaptiveControls,
    pub body_temperature_tolerance_k: f64,
    /// Actual mass times the specific-enthalpy increment; no absolute datum.
    pub body_energy_tolerance_j: f64,
    /// Latent evolution remains controlled on an isothermal phase plateau.
    pub body_liquid_fraction_tolerance: f64,
    /// Applied separately to each body's convection and radiation.
    pub heat_tolerance_j: f64,
}

#[derive(Debug)]
pub enum MolecularNetworkThermalTrajectoryError<E = Infallible> {
    Adaptive(MolecularNetworkTrajectoryError),
    Step {
        elapsed_s: f64,
        duration_s: f64,
        error: MolecularNetworkThermalError<E>,
    },
}
impl<E> From<MolecularNetworkTrajectoryError> for MolecularNetworkThermalTrajectoryError<E> {
    fn from(error: MolecularNetworkTrajectoryError) -> Self {
        Self::Adaptive(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for MolecularNetworkThermalTrajectoryError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive thermal molecular network refused: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for MolecularNetworkThermalTrajectoryError<E> {}
type Error<E = Infallible> = MolecularNetworkThermalTrajectoryError<E>;

#[derive(Clone, Debug, PartialEq)]
pub struct MolecularNetworkBodyHeatTotals {
    pub id: u64,
    pub gas: u64,
    pub solid_to_gas_heat_j: f64,
    pub radiation_heat_j: f64,
    pub internal_heat_j: f64,
    pub storage_change_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct MolecularNetworkThermalTrajectoryStep<S = ()> {
    pub duration_s: f64,
    pub energy_frame_velocity_m_s: [f64; 3],
    pub attempted_intervals: usize,
    pub accepted_intervals: usize,
    pub rejected_intervals: usize,
    pub maximum_accepted_error: f64,
    /// Exactly two original native transactions per accepted interval.
    pub accepted_steps: Vec<MolecularNetworkThermalStep<S>>,
    pub links: Vec<MolecularNetworkLinkTotals>,
    pub bodies: Vec<MolecularNetworkBodyHeatTotals>,
    /// The authored whole-window gas/body sources plus actual radiation.
    pub external_heat_j: f64,
    pub body_energy_change_j: f64,
    /// Actual rounded prescribed source sum minus the authored window request.
    pub source_partition_residual_j: f64,
    pub residual: MolecularNetworkWindowBalance,
    /// Absolute native defects, actual stage closure and source partition
    /// roundoff consume one fixed physical budget; signs never refund it.
    pub debit: MolecularNetworkWindowBalance,
}

#[derive(Clone, Debug, PartialEq)]
pub struct MolecularNetworkThermalTrajectoryCandidate<S = ()> {
    gases: BTreeMap<u64, MovingGasParcel>,
    bodies: BTreeMap<u64, EquilibriumPhaseState>,
    step: MolecularNetworkThermalTrajectoryStep<S>,
}
impl<S> MolecularNetworkThermalTrajectoryCandidate<S> {
    pub fn gases(&self) -> &BTreeMap<u64, MovingGasParcel> {
        &self.gases
    }
    pub fn bodies(&self) -> &BTreeMap<u64, EquilibriumPhaseState> {
        &self.bodies
    }
    pub const fn step(&self) -> &MolecularNetworkThermalTrajectoryStep<S> {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        BTreeMap<u64, MovingGasParcel>,
        BTreeMap<u64, EquilibriumPhaseState>,
        MolecularNetworkThermalTrajectoryStep<S>,
    ) {
        (self.gases, self.bodies, self.step)
    }
}

fn invalid(detail: &'static str) -> MolecularNetworkTrajectoryError {
    MolecularNetworkTrajectoryError::Invalid(detail)
}
fn unresolved(at: f64, detail: &'static str) -> MolecularNetworkTrajectoryError {
    MolecularNetworkTrajectoryError::Unrepresentable {
        elapsed_s: at,
        detail,
    }
}
fn poll(cx: &Cx<'_>, at: f64) -> Result<(), MolecularNetworkTrajectoryError> {
    cx.checkpoint()
        .map_err(|_| MolecularNetworkTrajectoryError::Cancelled { elapsed_s: at })
}
fn value(x: Dd, at: f64) -> Result<f64, MolecularNetworkTrajectoryError> {
    let rounded = x.to_f64();
    if x.hi.is_finite() && x.lo.is_finite() && rounded.is_finite() {
        Ok(rounded)
    } else {
        Err(unresolved(at, "thermal trajectory account"))
    }
}
fn scaled_heat(
    heat: f64,
    dt: f64,
    window: f64,
    at: f64,
) -> Result<f64, MolecularNetworkTrajectoryError> {
    let scaled = value(
        Dd::from_f64(heat) * Dd::from_f64(dt) / Dd::from_f64(window),
        at,
    )?;
    if heat != 0.0 && scaled == 0.0 {
        Err(unresolved(at, "prescribed thermal source underflow"))
    } else {
        Ok(scaled)
    }
}

#[derive(Clone, Copy)]
struct Body {
    gas: u64,
    initial: EquilibriumPhaseState,
    mass_kg: f64,
    internal_heat_j: f64,
}
struct Layout {
    bodies: BTreeMap<u64, Body>,
    gas_sources: BTreeMap<u64, f64>,
    source_heat_j: Dd,
}
impl Layout {
    fn new(
        cx: &Cx<'_>,
        gases: &BTreeMap<u64, MovingGasParcel>,
        boundaries: &[MolecularGasThermalBoundary<'_>],
        duration: f64,
        controls: MolecularNetworkThermalControls,
        adaptive: MolecularNetworkThermalAdaptiveControls,
    ) -> Result<Self, MolecularNetworkTrajectoryError> {
        if controls.maximum_bodies > MAX_NETWORK_THERMAL_BODIES
            || ![
                adaptive.body_temperature_tolerance_k,
                adaptive.body_energy_tolerance_j,
                adaptive.body_liquid_fraction_tolerance,
                adaptive.heat_tolerance_j,
            ]
            .into_iter()
            .all(|x| x.is_finite() && x > 0.0)
            || adaptive.body_liquid_fraction_tolerance > 1.0
        {
            return Err(invalid(
                "positive finite thermal state/heat estimates required",
            ));
        }
        let mut result = Self {
            bodies: BTreeMap::new(),
            gas_sources: BTreeMap::new(),
            source_heat_j: Dd::ZERO,
        };
        for boundary in boundaries {
            poll(cx, 0.0)?;
            if !gases.contains_key(&boundary.gas)
                || result
                    .gas_sources
                    .insert(boundary.gas, boundary.gas_heat_j)
                    .is_some()
                || boundary.bodies.is_empty()
                || boundary.bodies.len() > MAX_BODIES
                || !boundary.gas_heat_j.is_finite()
            {
                return Err(invalid(
                    "unique existing thermal gases and nonempty body groups required",
                ));
            }
            for body in boundary.bodies {
                poll(cx, 0.0)?;
                let input = body.input;
                if input.duration_s != duration
                    || !input.internal_heat_j.is_finite()
                    || !input.mass_kg.is_finite()
                    || input.mass_kg <= 0.0
                    || result
                        .bodies
                        .insert(
                            body.id,
                            Body {
                                gas: boundary.gas,
                                initial: input.initial,
                                mass_kg: input.mass_kg,
                                internal_heat_j: input.internal_heat_j,
                            },
                        )
                        .is_some()
                {
                    return Err(invalid(
                        "distinct bodies must retain positive mass and the whole-window duration/sources",
                    ));
                }
                if result.bodies.len() > controls.maximum_bodies {
                    return Err(invalid("thermal body count exceeds the native work limit"));
                }
            }
        }
        for &heat in result.gas_sources.values() {
            poll(cx, 0.0)?;
            result.source_heat_j = result.source_heat_j + Dd::from_f64(heat);
        }
        for body in result.bodies.values() {
            poll(cx, 0.0)?;
            result.source_heat_j = result.source_heat_j + Dd::from_f64(body.internal_heat_j);
        }
        let work = adaptive
            .molecular
            .states
            .maximum_accepted_intervals
            .checked_mul(2)
            .and_then(|n| n.checked_mul(result.bodies.len().max(1)));
        if work.is_none() {
            return Err(invalid("bounded thermal history storage required"));
        }
        value(result.source_heat_j, 0.0)?;
        Ok(result)
    }

    fn initial(&self) -> BTreeMap<u64, EquilibriumPhaseState> {
        self.bodies
            .iter()
            .map(|(&id, body)| (id, body.initial))
            .collect()
    }
    fn storage(
        &self,
        after: &BTreeMap<u64, EquilibriumPhaseState>,
        at: f64,
        cx: &Cx<'_>,
    ) -> Result<Dd, MolecularNetworkTrajectoryError> {
        if !self.bodies.keys().eq(after.keys()) {
            return Err(invalid("native thermal body identities changed"));
        }
        let mut sum = Dd::ZERO;
        for (&id, body) in &self.bodies {
            poll(cx, at)?;
            let state = after[&id];
            if state.phase_curve_identity() != body.initial.phase_curve_identity() {
                return Err(invalid(
                    "native thermal body changed its original caloric chart",
                ));
            }
            sum = sum
                + Dd::from_f64(body.mass_kg)
                    * (Dd::from_f64(state.specific_enthalpy_j_kg())
                        - Dd::from_f64(body.initial.specific_enthalpy_j_kg()));
        }
        value(sum, at)?;
        Ok(sum)
    }
}

/// Rebuild immutable descriptors only; the original native owners perform all
/// caloric inversions, convection/radiation solves and molecular transactions.
#[allow(clippy::too_many_arguments)]
fn trial<S, E>(
    cx: &Cx<'_>,
    gases: &BTreeMap<u64, MovingGasParcel>,
    bodies: &BTreeMap<u64, EquilibriumPhaseState>,
    links: &[MolecularNetworkLink<'_>],
    duration: f64,
    window: f64,
    at: f64,
    endpoints: &BTreeMap<u64, MovingGasBoundarySolve>,
    boundaries: &[MolecularGasThermalBoundary<'_>],
    controls: MolecularNetworkThermalControls,
    resolve: &mut impl FnMut(
        u64,
        &FiniteGasBody<'_>,
        &ResolvedThermochemicalGasState,
        [f64; 3],
    ) -> Result<ConvectionCoefficient<S>, E>,
) -> Result<MolecularNetworkThermalCandidate<S>, Error<E>> {
    poll(cx, at)?;
    let mut proposals = Vec::with_capacity(boundaries.len());
    for boundary in boundaries {
        poll(cx, at)?;
        let mut group = Vec::with_capacity(boundary.bodies.len());
        for body in boundary.bodies {
            poll(cx, at)?;
            let mut body = *body;
            body.input.initial = *bodies
                .get(&body.id)
                .ok_or_else(|| invalid("trial lost an original thermal body"))?;
            body.input.duration_s = duration;
            body.input.internal_heat_j =
                scaled_heat(body.input.internal_heat_j, duration, window, at)?;
            group.push(body);
        }
        proposals.push(group);
    }
    let mut sources = Vec::with_capacity(boundaries.len());
    for (boundary, bodies) in boundaries.iter().zip(&proposals) {
        sources.push(MolecularGasThermalBoundary {
            gas: boundary.gas,
            bodies,
            gas_heat_j: scaled_heat(boundary.gas_heat_j, duration, window, at)?,
            exchange: boundary.exchange,
            maximum_gas_temperature_increment_k: boundary.maximum_gas_temperature_increment_k,
        });
    }
    step_candidate_with_convection(
        cx, gases, links, duration, endpoints, &sources, controls, resolve,
    )
    .map_err(|error| match error {
        MolecularNetworkThermalError::Cancelled => {
            MolecularNetworkTrajectoryError::Cancelled { elapsed_s: at }.into()
        }
        error => Error::Step {
            elapsed_s: at,
            duration_s: duration,
            error,
        },
    })
}

fn gas_momentum_error(
    before: &BTreeMap<u64, MovingGasParcel>,
    coarse: &BTreeMap<u64, MovingGasParcel>,
    fine: &BTreeMap<u64, MovingGasParcel>,
    frame: [f64; 3],
    controls: MolecularNetworkAdaptiveControls,
    at: f64,
    cx: &Cx<'_>,
) -> Result<f64, MolecularNetworkTrajectoryError> {
    let mut error = 0.0_f64;
    for (&id, old) in before {
        poll(cx, at)?;
        for d in 0..3 {
            // Momentum increments in the declared translating energy frame.
            // Products and differences retain small recoil under common boosts.
            let increment = |gas: &MovingGasParcel| {
                value(
                    Dd::from_f64(gas.momentum_kg_m_s()[d])
                        - Dd::from_f64(old.momentum_kg_m_s()[d])
                        - Dd::from_f64(frame[d])
                            * (Dd::from_f64(gas.parcel().mass_kg())
                                - Dd::from_f64(old.parcel().mass_kg())),
                    at,
                )
            };
            let (a, b) = (increment(&coarse[&id])?, increment(&fine[&id])?);
            error = error.max(molecular::ratio(
                a,
                b,
                a.abs().max(b.abs()),
                controls.stream_momentum_tolerance_n_s,
                controls.states.relative_tolerance,
                at,
            )?);
        }
    }
    Ok(error)
}

fn body_error(
    layout: &Layout,
    before: &BTreeMap<u64, EquilibriumPhaseState>,
    coarse: &BTreeMap<u64, EquilibriumPhaseState>,
    fine: &BTreeMap<u64, EquilibriumPhaseState>,
    controls: MolecularNetworkThermalAdaptiveControls,
    at: f64,
    cx: &Cx<'_>,
) -> Result<f64, MolecularNetworkTrajectoryError> {
    if !layout.bodies.keys().eq(before.keys())
        || !before.keys().eq(coarse.keys())
        || !before.keys().eq(fine.keys())
    {
        return Err(invalid("native thermal state identities changed"));
    }
    let mut error = 0.0_f64;
    for (&id, body) in &layout.bodies {
        poll(cx, at)?;
        let [old, full, half] = [before[&id], coarse[&id], fine[&id]];
        let energy = |state: EquilibriumPhaseState| {
            value(
                Dd::from_f64(body.mass_kg)
                    * (Dd::from_f64(state.specific_enthalpy_j_kg())
                        - Dd::from_f64(old.specific_enthalpy_j_kg())),
                at,
            )
        };
        let [a, b] = [energy(full)?, energy(half)?];
        for (full, half, scale, absolute) in [
            (
                full.temperature_k(),
                half.temperature_k(),
                old.temperature_k()
                    .max(full.temperature_k())
                    .max(half.temperature_k()),
                controls.body_temperature_tolerance_k,
            ),
            (a, b, a.abs().max(b.abs()), controls.body_energy_tolerance_j),
            (
                full.liquid_mass_fraction(),
                half.liquid_mass_fraction(),
                1.0,
                controls.body_liquid_fraction_tolerance,
            ),
        ] {
            error = error.max(molecular::ratio(
                full,
                half,
                scale,
                absolute,
                controls.molecular.states.relative_tolerance,
                at,
            )?);
        }
    }
    Ok(error)
}

fn heat_error<S>(
    steps: [&MolecularNetworkThermalStep<S>; 3],
    controls: MolecularNetworkThermalAdaptiveControls,
    at: f64,
    cx: &Cx<'_>,
) -> Result<f64, MolecularNetworkTrajectoryError> {
    let mut rows = [BTreeMap::new(), BTreeMap::new(), BTreeMap::new()];
    for (map, step) in rows.iter_mut().zip(steps) {
        for (&gas, thermal) in &step.thermal {
            poll(cx, at)?;
            for body in thermal.bodies() {
                let report = &body.report;
                if map
                    .insert(
                        report.id,
                        (
                            gas,
                            [
                                report.solid_to_gas_heat_j,
                                report.radiation_heat_j,
                                report.enthalpy_change_j,
                            ],
                        ),
                    )
                    .is_some()
                {
                    return Err(invalid("native heat receipt duplicated a physical body"));
                }
            }
        }
    }
    if !rows[0].keys().eq(rows[1].keys()) || !rows[0].keys().eq(rows[2].keys()) {
        return Err(invalid("native heat receipt lost a physical body"));
    }
    let mut error = 0.0_f64;
    for (&id, (gas, full)) in &rows[0] {
        poll(cx, at)?;
        let (one_gas, one) = rows[1][&id];
        let (two_gas, two) = rows[2][&id];
        if *gas != one_gas || *gas != two_gas {
            return Err(invalid("native thermal contact changed its gas identity"));
        }
        for k in 0..3 {
            let absolute = if k == 2 {
                controls.body_energy_tolerance_j
            } else {
                controls.heat_tolerance_j
            };
            // Keep each direction of physical convection/radiation even if
            // its sign reverses during the accepted half steps.
            for direction in 0..if k == 2 { 1 } else { 3 } {
                let select = |q: f64| match direction {
                    1 => q.max(0.0),
                    2 => (-q).max(0.0),
                    _ => q,
                };
                let coarse = select(full[k]);
                let half = value(
                    Dd::from_f64(select(one[k])) + Dd::from_f64(select(two[k])),
                    at,
                )?;
                error = error.max(molecular::ratio(
                    coarse,
                    half,
                    coarse.abs().max(half.abs()),
                    absolute,
                    controls.molecular.states.relative_tolerance,
                    at,
                )?);
            }
        }
    }
    Ok(error)
}

#[derive(Default)]
struct BodyTotal {
    heat: Dd,
    radiation: Dd,
    internal: Dd,
}
struct ThermalTotals {
    bodies: BTreeMap<u64, BodyTotal>,
    source: Dd,
    radiation: Dd,
    source_roundoff: Dd,
    energy_debit: [Dd; 2],
    momentum_debit: [Dd; 3],
}
impl Default for ThermalTotals {
    fn default() -> Self {
        Self {
            bodies: BTreeMap::new(),
            source: Dd::ZERO,
            radiation: Dd::ZERO,
            source_roundoff: Dd::ZERO,
            energy_debit: [Dd::ZERO; 2],
            momentum_debit: [Dd::ZERO; 3],
        }
    }
}
impl ThermalTotals {
    fn accept<S>(
        &mut self,
        step: &MolecularNetworkThermalStep<S>,
        layout: &Layout,
        window: f64,
        at: f64,
        cx: &Cx<'_>,
    ) -> Result<(), MolecularNetworkTrajectoryError> {
        for (&gas, thermal) in &step.thermal {
            poll(cx, at)?;
            let report = thermal.report();
            self.source = self.source
                + Dd::from_f64(report.external_gas_heat_j)
                + Dd::from_f64(report.internal_body_heat_j);
            self.radiation = self.radiation + Dd::from_f64(report.radiation_heat_j);
            for body in thermal.bodies() {
                poll(cx, at)?;
                let report = &body.report;
                let original = layout
                    .bodies
                    .get(&report.id)
                    .ok_or_else(|| invalid("accepted heat lost its original body"))?;
                if original.gas != gas {
                    return Err(invalid("accepted heat changed its original gas"));
                }
                let internal = scaled_heat(original.internal_heat_j, step.duration_s, window, at)?;
                let total = self.bodies.entry(report.id).or_default();
                total.heat = total.heat + Dd::from_f64(report.solid_to_gas_heat_j);
                total.radiation = total.radiation + Dd::from_f64(report.radiation_heat_j);
                total.internal = total.internal + Dd::from_f64(internal);
            }
        }
        for heat in layout
            .gas_sources
            .values()
            .copied()
            .chain(layout.bodies.values().map(|b| b.internal_heat_j))
        {
            poll(cx, at)?;
            let exact = Dd::from_f64(heat) * Dd::from_f64(step.duration_s) / Dd::from_f64(window);
            let rounded = Dd::from_f64(scaled_heat(heat, step.duration_s, window, at)?);
            self.source_roundoff = self.source_roundoff + (rounded - exact).abs();
        }
        for debit in &mut self.energy_debit {
            *debit = *debit + Dd::from_f64(step.thermal_absolute_energy_defect_j);
        }
        for d in 0..3 {
            self.momentum_debit[d] =
                self.momentum_debit[d] + Dd::from_f64(step.thermal_absolute_momentum_defect_n_s[d]);
        }
        Ok(())
    }

    fn source_partition(&self, layout: &Layout, at: f64, window: f64) -> Dd {
        self.source - layout.source_heat_j * Dd::from_f64(at) / Dd::from_f64(window)
    }
    fn external(&self, layout: &Layout, at: f64, window: f64) -> Dd {
        layout.source_heat_j * Dd::from_f64(at) / Dd::from_f64(window) + self.radiation
    }
    fn debit(
        &self,
        mut molecular: MolecularNetworkWindowBalance,
        layout: &Layout,
        at: f64,
        window: f64,
    ) -> Result<MolecularNetworkWindowBalance, MolecularNetworkTrajectoryError> {
        let source = self.source_roundoff + self.source_partition(layout, at, window).abs();
        molecular.energy_j = value(
            Dd::from_f64(molecular.energy_j) + self.energy_debit[0] + source,
            at,
        )?;
        molecular.laboratory_energy_j = value(
            Dd::from_f64(molecular.laboratory_energy_j) + self.energy_debit[1] + source,
            at,
        )?;
        for d in 0..3 {
            molecular.momentum_n_s[d] = value(
                Dd::from_f64(molecular.momentum_n_s[d]) + self.momentum_debit[d],
                at,
            )?;
        }
        Ok(molecular)
    }
    fn reports(
        &self,
        layout: &Layout,
        bodies: &BTreeMap<u64, EquilibriumPhaseState>,
        at: f64,
        cx: &Cx<'_>,
    ) -> Result<Vec<MolecularNetworkBodyHeatTotals>, MolecularNetworkTrajectoryError> {
        layout
            .bodies
            .iter()
            .map(|(&id, original)| {
                poll(cx, at)?;
                let total = self
                    .bodies
                    .get(&id)
                    .ok_or_else(|| invalid("accepted history lost an original thermal body"))?;
                let storage = Dd::from_f64(original.mass_kg)
                    * (Dd::from_f64(bodies[&id].specific_enthalpy_j_kg())
                        - Dd::from_f64(original.initial.specific_enthalpy_j_kg()));
                Ok(MolecularNetworkBodyHeatTotals {
                    id,
                    gas: original.gas,
                    solid_to_gas_heat_j: value(total.heat, at)?,
                    radiation_heat_j: value(total.radiation, at)?,
                    internal_heat_j: value(total.internal, at)?,
                    storage_change_j: value(storage, at)?,
                })
            })
            .collect()
    }
}

/// Evolve with each original body's declared whole-surface film coefficient.
#[allow(clippy::too_many_arguments)]
pub fn advance_candidate(
    cx: &Cx<'_>,
    gases: &BTreeMap<u64, MovingGasParcel>,
    links: &[MolecularNetworkLink<'_>],
    duration_s: f64,
    endpoints: &BTreeMap<u64, MovingGasBoundarySolve>,
    boundaries: &[MolecularGasThermalBoundary<'_>],
    controls: MolecularNetworkThermalControls,
    adaptive: MolecularNetworkThermalAdaptiveControls,
) -> Result<MolecularNetworkThermalTrajectoryCandidate, Error> {
    advance_candidate_with_convection(
        cx,
        gases,
        links,
        duration_s,
        endpoints,
        boundaries,
        controls,
        adaptive,
        |_, body, _, _| {
            Ok(ConvectionCoefficient {
                coefficient_w_per_m2_k: body.environment.convection_w_per_m2_k,
                source: (),
            })
        },
    )
}

/// Evolve one complete window with original source and thermal owners.
///
/// Body input durations equal this window; prescribed heats are held powers.
/// Every trial uses accepted body phase states and gas composition/momentum.
/// The callback is repeatable and receives each native proposed gas state; its
/// external side effects cannot be rolled back. Only typed molecular gross
/// withdrawal permits a physical retry. Source, phase-domain, convection-law,
/// inversion and conservation refusals remain hard errors.
///
/// Storage is bounded by twice the accepted-interval limit times the admitted
/// native network/body counts. Cancellation or any refusal publishes no state.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn advance_candidate_with_convection<S, E>(
    cx: &Cx<'_>,
    gases: &BTreeMap<u64, MovingGasParcel>,
    links: &[MolecularNetworkLink<'_>],
    duration_s: f64,
    endpoints: &BTreeMap<u64, MovingGasBoundarySolve>,
    boundaries: &[MolecularGasThermalBoundary<'_>],
    controls: MolecularNetworkThermalControls,
    adaptive: MolecularNetworkThermalAdaptiveControls,
    mut resolve: impl FnMut(
        u64,
        &FiniteGasBody<'_>,
        &ResolvedThermochemicalGasState,
        [f64; 3],
    ) -> Result<ConvectionCoefficient<S>, E>,
) -> Result<MolecularNetworkThermalTrajectoryCandidate<S>, Error<E>> {
    poll(cx, 0.0)?;
    molecular::admit(
        gases,
        links,
        endpoints,
        duration_s,
        controls.exchange,
        adaptive.molecular,
    )?;
    let layout = Layout::new(cx, gases, boundaries, duration_s, controls, adaptive)?;
    let mut cancelled = || cx.checkpoint().is_err();
    let sources = Sources::with_thermal_gases(
        gases,
        links,
        boundaries.iter().map(|b| b.gas),
        &mut cancelled,
    )?;
    let a = adaptive.molecular.states;
    let mut current = gases.clone();
    let mut bodies = layout.initial();
    let (mut elapsed, mut next_step, mut maximum_error) = (0.0, a.initial_step_s, 0.0_f64);
    let (mut attempts, mut rejected, mut accepted) = (0, 0, 0);
    let mut history = Vec::new();
    let mut totals = Totals::default();
    let mut thermal = ThermalTotals::default();
    let mut residual = MolecularNetworkWindowBalance::default();
    while elapsed < duration_s {
        poll(cx, elapsed)?;
        if attempts == a.maximum_attempts || accepted == a.maximum_accepted_intervals {
            return Err(MolecularNetworkTrajectoryError::WorkBudget {
                elapsed_s: elapsed,
                attempted_intervals: attempts,
                accepted_intervals: accepted,
            }
            .into());
        }
        let remaining = duration_s - elapsed;
        let mut h = next_step
            .min(controls.exchange.maximum_duration_s)
            .min(remaining);
        let tail = remaining - h;
        if tail > 0.0 && tail < a.minimum_step_s && 0.5 * remaining >= a.minimum_step_s {
            h = 0.5 * remaining;
        }
        let (half, end, middle_time) = (0.5 * h, elapsed + h, elapsed + 0.5 * h);
        if half <= 0.0
            || half + half != h
            || end <= elapsed
            || end > duration_s
            || middle_time <= elapsed
            || middle_time >= end
            || (h == remaining && end != duration_s)
        {
            return Err(unresolved(
                elapsed,
                "trial halves or final physical time cannot advance exactly",
            )
            .into());
        }
        attempts += 1;
        let proposal = (|| {
            let mut step = |gas, body, dt, at| {
                trial(
                    cx,
                    gas,
                    body,
                    links,
                    dt,
                    duration_s,
                    at,
                    endpoints,
                    boundaries,
                    controls,
                    &mut resolve,
                )
            };
            let coarse = step(&current, &bodies, h, elapsed)?;
            let middle = step(&current, &bodies, half, elapsed)?;
            let fine = step(middle.gases(), middle.bodies(), half, middle_time)?;
            let error = molecular::state_error(
                &current,
                coarse.gases(),
                fine.gases(),
                &sources,
                a,
                elapsed,
                &mut cancelled,
            )?
            .max(stream_error_from_states(
                [
                    coarse.molecular_initial_gases(),
                    middle.molecular_initial_gases(),
                    fine.molecular_initial_gases(),
                ],
                [
                    &coarse.step().molecular,
                    &middle.step().molecular,
                    &fine.step().molecular,
                ],
                &sources,
                adaptive.molecular,
                elapsed,
                &mut cancelled,
            )?)
            .max(gas_momentum_error(
                &current,
                coarse.gases(),
                fine.gases(),
                controls.exchange.energy_frame_velocity_m_s,
                adaptive.molecular,
                elapsed,
                cx,
            )?)
            .max(body_error(
                &layout,
                &bodies,
                coarse.bodies(),
                fine.bodies(),
                adaptive,
                elapsed,
                cx,
            )?)
            .max(heat_error(
                [coarse.step(), middle.step(), fine.step()],
                adaptive,
                elapsed,
                cx,
            )?);
            Ok::<_, Error<E>>((middle, fine, error))
        })();
        let (middle, fine, error) = match proposal {
            Ok(result) => result,
            Err(Error::Step {
                error:
                    MolecularNetworkThermalError::Molecular(MolecularNetworkError::Withdrawal {
                        ..
                    }),
                ..
            }) if half >= a.minimum_step_s => {
                rejected += 1;
                next_step = half;
                continue;
            }
            Err(error) => return Err(error),
        };
        if error > 1.0 {
            rejected += 1;
            if half < a.minimum_step_s {
                return Err(MolecularNetworkTrajectoryError::MinimumStep {
                    elapsed_s: elapsed,
                    duration_s: h,
                    normalized_error: error,
                }
                .into());
            }
            next_step = half;
            continue;
        }
        for (candidate, at) in [(&middle, middle_time), (&fine, end)] {
            totals.accept(
                candidate.molecular_initial_gases(),
                candidate.gases(),
                &candidate.step().molecular,
                &sources,
                at,
                &mut cancelled,
            )?;
            thermal.accept(candidate.step(), &layout, duration_s, at, cx)?;
        }
        residual = window_balance(
            gases,
            fine.gases(),
            &sources,
            controls.exchange.energy_frame_velocity_m_s,
            end,
            &mut cancelled,
        )?;
        let body_change = layout.storage(fine.bodies(), end, cx)?;
        let external = thermal.external(&layout, end, duration_s);
        residual.energy_j = value(
            Dd::from_f64(residual.energy_j) + body_change - external,
            end,
        )?;
        residual.laboratory_energy_j = value(
            Dd::from_f64(residual.laboratory_energy_j) + body_change - external,
            end,
        )?;
        let debit = thermal.debit(totals.debit(end)?, &layout, end, duration_s)?;
        molecular::admit_window(residual, debit, adaptive.molecular, end)?;
        poll(cx, end)?;
        let (_, _, first) = middle.into_parts();
        let (next_gases, next_bodies, second) = fine.into_parts();
        history.extend([first, second]);
        current = next_gases;
        bodies = next_bodies;
        elapsed = end;
        accepted += 1;
        maximum_error = maximum_error.max(error);
        next_step = if error <= 0.25 {
            (2.0 * h).min(a.maximum_step_s)
        } else {
            h
        };
    }
    let debit = thermal.debit(totals.debit(elapsed)?, &layout, elapsed, duration_s)?;
    let links = totals.links(elapsed, &mut cancelled)?;
    let body_reports = thermal.reports(&layout, &bodies, elapsed, cx)?;
    let external_heat_j = value(thermal.external(&layout, elapsed, duration_s), elapsed)?;
    let body_energy_change_j = value(layout.storage(&bodies, elapsed, cx)?, elapsed)?;
    let source_partition_residual_j = value(
        thermal.source_partition(&layout, elapsed, duration_s),
        elapsed,
    )?;
    poll(cx, elapsed)?;
    Ok(MolecularNetworkThermalTrajectoryCandidate {
        gases: current,
        bodies,
        step: MolecularNetworkThermalTrajectoryStep {
            duration_s: elapsed,
            energy_frame_velocity_m_s: controls.exchange.energy_frame_velocity_m_s,
            attempted_intervals: attempts,
            accepted_intervals: accepted,
            rejected_intervals: rejected,
            maximum_accepted_error: maximum_error,
            accepted_steps: history,
            links,
            bodies: body_reports,
            external_heat_j,
            body_energy_change_j,
            source_partition_residual_j,
            residual,
            debit,
        },
    })
}
