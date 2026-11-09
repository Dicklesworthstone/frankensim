//! Simultaneous finite-reservoir contact and native slab pressure work.
//!
//! Backward Euler evaluates contact at BOTH mechanically corrected endpoints.
//! Original caloric functions, the existing bounded heat root and the original
//! pressure/remap retain ownership. The split reservoir entry remains available.

use super::*;
use fs_material::phase::PhaseStateError;

/// Separate allowances for original native work, represented reservoir storage,
/// joint energy and the endpoint contact equation. All are absolute joules.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreeSurfaceSlabCoupledReservoirHeatControls {
    pub heat: FreeSurfaceSlabHeatControls,
    pub maximum_reservoir_energy_defect_j: f64,
    pub maximum_total_energy_defect_j: f64,
    pub maximum_contact_heat_defect_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum FreeSurfaceSlabCoupledReservoirHeatError {
    Native(FreeSurfaceSlabHeatError),
    Caloric(PhaseStateError),
}
impl core::fmt::Display for FreeSurfaceSlabCoupledReservoirHeatError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "coupled finite-reservoir thermal slab: {self:?}")
    }
}
impl std::error::Error for FreeSurfaceSlabCoupledReservoirHeatError {}
impl From<FreeSurfaceSlabHeatError> for FreeSurfaceSlabCoupledReservoirHeatError {
    fn from(error: FreeSurfaceSlabHeatError) -> Self {
        Self::Native(error)
    }
}
impl From<Error> for FreeSurfaceSlabCoupledReservoirHeatError {
    fn from(error: Error) -> Self {
        Self::Native(error.into())
    }
}
impl From<PhaseStateError> for FreeSurfaceSlabCoupledReservoirHeatError {
    fn from(error: PhaseStateError) -> Self {
        Self::Caloric(error)
    }
}

/// An immutable joint proposal; the original liquid and reservoir remain intact.
pub struct FreeSurfaceSlabCoupledReservoirHeatPlan<'a, 's, 'm, 'r> {
    native: FreeSurfaceSlabHeatPlan<'a, 's, 'm>,
    reservoir: FinitePhaseBody<'r>,
    reservoir_endpoint: EquilibriumPhaseState,
    scalar_contact_heat_defect_j: f64,
    controls: FreeSurfaceSlabCoupledReservoirHeatControls,
    conductance_w_k: f64,
    duration_s: f64,
    exposure_j_k: Dd,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreeSurfaceSlabCoupledReservoirHeatReport {
    pub duration_s: f64,
    pub conductance_w_k: f64,
    /// Actual native liquid Delta H plus measured correction-pressure work [J].
    pub heat_to_liquid_j: f64,
    /// G*dt*(T_reservoir,end - T_liquid,end), using represented endpoints [J].
    pub endpoint_contact_heat_j: f64,
    /// Actual heat into the liquid minus the endpoint contact equation [J].
    pub contact_heat_defect_j: f64,
    /// Original scalar contact residual, before independent native readback [J].
    pub scalar_contact_heat_defect_j: f64,
    pub reservoir_enthalpy_change_j: f64,
    /// Represented reservoir Delta H plus the scalar's supplied heat [J].
    pub reservoir_energy_defect_j: f64,
    /// Actual liquid Delta H + measured native work + reservoir Delta H [J].
    pub total_energy_defect_j: f64,
}

/// Both accepted endpoints. The reservoir's fixed mass and original chart are
/// retained; its `initial` is the accepted state for a subsequent interval.
#[derive(Clone)]
pub struct FreeSurfaceSlabCoupledReservoirHeatStep<'s, 'm, 'r> {
    pub heat: FreeSurfaceSlabHeatStep<'s, 'm>,
    pub reservoir: FinitePhaseBody<'r>,
    pub report: FreeSurfaceSlabCoupledReservoirHeatReport,
}

impl<'a, 's, 'm, 'r> FreeSurfaceSlabCoupledReservoirHeatPlan<'a, 's, 'm, 'r> {
    #[must_use]
    pub const fn native(&self) -> &FreeSurfaceSlabHeatPlan<'a, 's, 'm> {
        &self.native
    }

    #[must_use]
    pub const fn reservoir_endpoint(&self) -> EquilibriumPhaseState {
        self.reservoir_endpoint
    }

    /// Execute ONE original moving-interface pressure/remap. No participant,
    /// clock or source debit is published unless every independent gate passes.
    pub fn advance<'next>(
        &'next self,
        gate: &CancelGate,
        cx: &Cx<'_>,
    ) -> Result<
        FreeSurfaceSlabCoupledReservoirHeatStep<'next, 'next, 'r>,
        FreeSurfaceSlabCoupledReservoirHeatError,
    >
    where
        'm: 'next,
    {
        let heat = self.native.advance(gate, cx)?;
        let actual = &heat.volume.report;
        let into_liquid = Dd::from_f64(actual.requested_enthalpy_increment_j)
            + Dd::from_f64(actual.enthalpy_defect_j)
            + Dd::from_f64(actual.required_mechanical_work_j);
        let reservoir_change = Dd::from_f64(self.reservoir.mass_kg)
            * (Dd::from_f64(self.reservoir_endpoint.specific_enthalpy_j_kg())
                - Dd::from_f64(self.reservoir.initial.specific_enthalpy_j_kg()));
        // This binding is installed only AFTER the original volume owner has
        // audited actual uniform native h, mass, geometry, momentum and energy.
        let mapped = heat.volume.liquid.volume_slab.as_ref().ok_or_else(|| {
            invalid("coupled reservoir lost its accepted native slab binding")
        })?;
        let liquid_temperature = self
            .native
            .original
            .caloric
            .temperature_at_specific_enthalpy(mapped.specific_enthalpy_j_kg)?;
        let contact = self.exposure_j_k
            * (Dd::from_f64(self.reservoir_endpoint.temperature_k())
                - Dd::from_f64(liquid_temperature));
        let reservoir_defect =
            value(reservoir_change + Dd::from_f64(self.native.supplied_heat_j));
        let total_defect = value(into_liquid + reservoir_change);
        let contact_defect = value(into_liquid - contact);
        if !reservoir_defect.is_finite()
            || !total_defect.is_finite()
            || !contact_defect.is_finite()
            || reservoir_defect.abs() > self.controls.maximum_reservoir_energy_defect_j
            || total_defect.abs() > self.controls.maximum_total_energy_defect_j
            || contact_defect.abs() > self.controls.maximum_contact_heat_defect_j
        {
            return Err(invalid(
                "represented coupled reservoir storage, joint energy or endpoint contact exceeds its gate",
            ));
        }
        let report = FreeSurfaceSlabCoupledReservoirHeatReport {
            duration_s: self.duration_s,
            conductance_w_k: self.conductance_w_k,
            heat_to_liquid_j: value(into_liquid),
            endpoint_contact_heat_j: value(contact),
            contact_heat_defect_j: contact_defect,
            scalar_contact_heat_defect_j: self.scalar_contact_heat_defect_j,
            reservoir_enthalpy_change_j: value(reservoir_change),
            reservoir_energy_defect_j: reservoir_defect,
            total_energy_defect_j: total_defect,
        };
        checkpoint(gate, cx)?;
        Ok(FreeSurfaceSlabCoupledReservoirHeatStep {
            heat,
            reservoir: FinitePhaseBody {
                initial: self.reservoir_endpoint,
                ..self.reservoir
            },
            report,
        })
    }
}

impl<'s, 'm> MaterialFreeSurfaceEnthalpy3<'s, 'm> {
    /// Simultaneous backward-Euler finite-capacity contact and native work.
    ///
    /// For the trial liquid increment d, the ORIGINAL pressure response gives
    /// W(d)=B*u*(u-u0). The one exchanged heat is Q=M*d+W, the reservoir endpoint
    /// is h_r1=h_r0-Q/m_r, and the scalar equation is
    /// Q=G*dt*(T_r(h_r1)-T_l(h_l0+d)). Thus mechanical cooling affects contact
    /// DURING this same interval, rather than following a provisional exchange.
    /// Original caloric temperature/derivative queries retain measured Cp,
    /// latent intervals and source identity. Their generalized knot tangents
    /// guide the safeguarded root; only the physical joule residual accepts it.
    ///
    /// This is a FIRST-ORDER endpoint discretization, not an exact continuous
    /// contact trajectory. Timestep refinement measures its time error. The
    /// caller's bracket selects the branch; every trial must remain in BOTH
    /// original source domains after paying work. A trial outside support
    /// refuses instead of clipping or extrapolating a reservoir temperature.
    ///
    /// The existing uniform-temperature pure-liquid slab, native incoming
    /// affine motion and pressure/remap admissions all apply. The reservoir
    /// may traverse its original equilibrium latent path. Both fixed masses,
    /// fixed G over the native dt and spatially uniform temperatures must be
    /// physically justified. The reservoir is a caloric body with no additional
    /// mechanical owner. Original fixed-pressure h already includes ambient
    /// p*dV; no second ambient work or invented pressure EOS is introduced.
    /// Projection/mapping losses remain separate and never become thermal heat.
    /// This adds no external power, mass transport, general ALE, viscous march
    /// or full coupled entropy claim. Both accepted endpoints publish together.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn prepare_uniform_slab_coupled_reservoir_heat<'a, 'r>(
        &'a self,
        anchored_face: Face3,
        reservoir: FinitePhaseBody<'r>,
        conductance_w_k: f64,
        controls: FreeSurfaceSlabCoupledReservoirHeatControls,
        gate: &CancelGate,
        cx: &Cx<'_>,
    ) -> Result<
        FreeSurfaceSlabCoupledReservoirHeatPlan<'a, 's, 'm, 'r>,
        FreeSurfaceSlabCoupledReservoirHeatError,
    > {
        checkpoint(gate, cx)?;
        if !reservoir.mass_kg.is_finite()
            || reservoir.mass_kg <= 0.0
            || [
                conductance_w_k,
                controls.maximum_reservoir_energy_defect_j,
                controls.maximum_total_energy_defect_j,
                controls.maximum_contact_heat_defect_j,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(invalid(
                "positive finite reservoir mass and finite nonnegative contact/gates required",
            ));
        }
        // The original owner checks identity, including the no-contact case.
        reservoir
            .curve
            .advance_specific_energy(reservoir.initial, 0.0)?;
        let prepared =
            self.prepare_uniform_slab_heat_state(anchored_face, controls.heat, gate, cx)?;
        let duration_s = prepared.dt;
        let exposure = Dd::from_f64(conductance_w_k) * Dd::from_f64(duration_s);
        if !value(exposure).is_finite()
            || (conductance_w_k > 0.0 && value(exposure) == 0.0)
        {
            return Err(invalid("coupled reservoir exposure is unrepresentable"));
        }
        let reservoir_mass = Dd::from_f64(reservoir.mass_kg);
        let mut accepted = None;
        let mut solution = balance::solve_residual(
            prepared.mass,
            0.0,
            controls.heat.net_specific_enthalpy_bracket_j_kg,
            controls
                .maximum_contact_heat_defect_j
                .min(controls.heat.maximum_heat_defect_j),
            controls.heat.maximum_source_evaluations,
            |increment| -> Result<balance::EnergyResidual, FreeSurfaceSlabCoupledReservoirHeatError> {
                let (liquid, mechanical) =
                    prepared.state_and_pressure_response(increment, gate, cx)?;
                let heat = prepared.mass * Dd::from_f64(increment) + mechanical.work_j;
                let heat_j = value(heat);
                let reservoir_increment = -value(heat / reservoir_mass);
                if !heat_j.is_finite()
                    || !reservoir_increment.is_finite()
                    || (heat_j != 0.0 && reservoir_increment == 0.0)
                {
                    return Err(invalid("coupled reservoir heat is unrepresentable"));
                }
                let endpoint = reservoir
                    .curve
                    .advance_specific_energy(reservoir.initial, reservoir_increment)?;
                let liquid_slope = self
                    .caloric
                    .temperature_derivative_at_specific_enthalpy(liquid.specific_enthalpy_j_kg())?;
                let reservoir_slope = reservoir
                    .curve
                    .temperature_derivative_at_specific_enthalpy(endpoint.specific_enthalpy_j_kg())?;
                let heat_slope = prepared.mass + mechanical.work_derivative_kg;
                let contact = exposure
                    * (Dd::from_f64(endpoint.temperature_k())
                        - Dd::from_f64(liquid.temperature_k()));
                accepted = Some((increment, endpoint, heat));
                Ok(balance::EnergyResidual {
                    residual_j: heat - contact,
                    derivative_kg: heat_slope
                        * (Dd::from_f64(1.0)
                            + exposure * Dd::from_f64(reservoir_slope) / reservoir_mass)
                        + exposure * Dd::from_f64(liquid_slope),
                    work_j: mechanical.work_j,
                })
            },
        )
        .map_err(|failure| match failure {
            balance::Failure::Invalid(what) => invalid(what),
            balance::Failure::Evaluation(error) => error,
            balance::Failure::Balance {
                evaluations,
                residual_j,
            } => FreeSurfaceSlabHeatError::Balance {
                evaluations,
                residual_j,
            }
            .into(),
        })?;
        let (increment, reservoir_endpoint, heat) = accepted
            .ok_or_else(|| invalid("coupled reservoir lost its accepted source response"))?;
        if increment.to_bits() != solution.increment_j_kg.to_bits() {
            return Err(invalid("coupled reservoir retained a different source trial"));
        }
        let scalar_contact_heat_defect_j = solution.residual_j;
        let heat_j = value(heat);
        // Preserve the native report's scalar HEAT meaning. The coupled root
        // solved contact; rounding its compensated M*d+W into supplied heat is
        // a distinct residual, retained before the native readback.
        solution.residual_j = value(heat - Dd::from_f64(heat_j));
        let native = prepared.finish(solution, heat_j, gate, cx)?;
        Ok(FreeSurfaceSlabCoupledReservoirHeatPlan {
            native,
            reservoir,
            reservoir_endpoint,
            scalar_contact_heat_defect_j,
            controls,
            conductance_w_k,
            duration_s,
            exposure_j_k: exposure,
        })
    }
}

fn invalid(what: &'static str) -> FreeSurfaceSlabCoupledReservoirHeatError {
    FreeSurfaceSlabHeatError::Invalid(what).into()
}
