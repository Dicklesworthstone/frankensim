//! Source-backed, isothermal Newtonian journal shear in articulated dynamics.
//!
//! The exact concentric-cylinder laminar resistance is
//! `c = 4 pi mu L ri² ro² / (ro²-ri²)`, with joint torque `-c qdot`.
//! Equal and opposite cylinder torques give loss `c qdot²`; the joint's native
//! work-conjugate coordinate accounts for both moving links. See Fitzpatrick,
//! *Taylor-Couette Flow*, equations 10.38–10.40:
//! <https://farside.ph.utexas.edu/teaching/336L/Fluidhtml/node137.html>.
//!
//! A film has negligible mass/thermal storage and equilibrates to its explicit
//! receiving wall cell. No pressure-supported bearing load, eccentricity,
//! cavitation, turbulence, end effects or finite-time fluid spin-up is modeled.
//! Geometry is prescribed even when link inertia thermally expands. Its
//! slenderness and quasisteady validity remain caller modeling obligations.
//! The creeping-flow limit is checked at every native mechanical RK stage.
//! Viscosity is frozen only within one explicitly capped thermal substep, then
//! refreshed from the same liquid claims at the accepted wall temperature.

use super::{ArticulatedModel, CellId, DepositionSimulation, Error, finite, positive};
use fs_matdb::ClaimSelection;
use fs_material::liquid::ResolvedLiquidState;
use fs_material::liquid_temperature::LiquidTemperatureLaw;
use fs_material::state_point::MaterialStatePointError;
use fs_mbd::articulated::JointKind;
use fs_qty::Temperature;

/// Prescribed concentric journal dimensions, in SI units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JournalGeometry {
    pub inner_radius_m: f64,
    pub outer_radius_m: f64,
    pub length_m: f64,
}

/// Explicit physical-domain and numerical-accuracy obligations for one film.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JournalLimits {
    /// `rho |qdot| ri (ro-ri) / mu`; positive and at most one.
    pub maximum_gap_reynolds: f64,
    /// Maximum fractional endpoint viscosity change within ONE thermal step.
    /// Positive and at most one. Refine motion `maximum_step_s` on refusal.
    pub maximum_relative_viscosity_change: f64,
    /// Absolute tolerance for this film's additional RK loss coordinate [J].
    pub absolute_loss_tolerance_j: f64,
}

/// One physical film; multiple distinct films may act on the same joint.
#[derive(Clone, Debug)]
pub struct JournalBearing<'a> {
    /// Child-link index, whose relative joint MUST be revolute.
    pub link: usize,
    /// Active thermal cell receiving all loss and setting film temperature.
    pub heat_cell: CellId,
    pub geometry: JournalGeometry,
    pub limits: JournalLimits,
    /// Immutable card, fixed nonthermal context and originally selected claims.
    pub liquid: LiquidTemperatureLaw<'a>,
}

/// Material and physical observations at the actual accepted joint endpoint.
#[derive(Clone, Debug, PartialEq)]
pub struct JournalBearingState {
    pub link: usize,
    pub heat_cell: CellId,
    pub liquid: ResolvedLiquidState,
    pub resistance_nm_s: f64,
    pub resisting_torque_nm: f64,
    pub dissipated_power_w: f64,
    pub gap_reynolds: f64,
}

#[derive(Clone)]
pub(super) struct BoundBearing<'a> {
    pub definition: JournalBearing<'a>,
    pub dof: usize,
    pub state: JournalBearingState,
    geometry_factor_m3: f64,
}

impl<'a> BoundBearing<'a> {
    pub fn new(
        definition: JournalBearing<'a>,
        model: &ArticulatedModel,
        thermal: &DepositionSimulation<'_>,
        velocity: &[f64],
    ) -> Result<Self, Error> {
        let g = definition.geometry;
        let limits = definition.limits;
        if ![g.inner_radius_m, g.outer_radius_m, g.length_m]
            .into_iter()
            .all(positive)
            || g.outer_radius_m <= g.inner_radius_m
            || !positive(limits.maximum_gap_reynolds)
            || limits.maximum_gap_reynolds > 1.0
            || !positive(limits.maximum_relative_viscosity_change)
            || limits.maximum_relative_viscosity_change > 1.0
            || !positive(limits.absolute_loss_tolerance_j)
        {
            return Err(Error::Invalid(
                "journal geometry and creeping-flow/loss bounds",
            ));
        }
        let link = model
            .links()
            .get(definition.link)
            .ok_or(Error::Invalid("journal child link does not exist"))?;
        if link.joint().kind() != JointKind::Revolute {
            return Err(Error::Invalid("journal shear requires a revolute joint"));
        }
        let dof = model
            .dof_index(definition.link)
            .ok_or(Error::Invalid("journal joint has no motion coordinate"))?;
        // Avoid subtracting squared nearly equal radii in a narrow gap.
        let factor = 4.0
            * core::f64::consts::PI
            * g.length_m
            * g.inner_radius_m.powi(2)
            * g.outer_radius_m.powi(2)
            / ((g.outer_radius_m - g.inner_radius_m) * (g.outer_radius_m + g.inner_radius_m));
        if !positive(factor) {
            return Err(Error::Invalid(
                "journal geometry resistance is unrepresentable",
            ));
        }
        let reference = definition.liquid.reference_state();
        if reference
            .material()
            .properties()
            .iter()
            .any(|p| !p.is_static_at_fixed_external_conditions())
        {
            return Err(Error::Invalid(
                "Newtonian film cannot freeze strain, rate or field dependent source claims",
            ));
        }
        let temperature = thermal.temperature_k(definition.heat_cell)?;
        let liquid = definition
            .liquid
            .at_temperature(Temperature::new(temperature))?;
        support(&definition.liquid, reference, &liquid)?;
        let state = JournalBearingState {
            link: definition.link,
            heat_cell: definition.heat_cell,
            resistance_nm_s: factor * liquid.dynamic_viscosity_pa_s(),
            liquid,
            resisting_torque_nm: 0.0,
            dissipated_power_w: 0.0,
            gap_reynolds: 0.0,
        };
        if !positive(state.resistance_nm_s) {
            return Err(Error::Invalid(
                "journal material resistance is unrepresentable",
            ));
        }
        let mut result = Self {
            definition,
            dof,
            state,
            geometry_factor_m3: factor,
        };
        result.update_rate(velocity[dof])?;
        Ok(result)
    }

    /// This exact resistance supplies both native torque and heat derivative.
    pub fn response(&self, rate: f64) -> Result<(f64, f64, f64), Error> {
        let g = self.definition.geometry;
        let reynolds = finite(
            rate.abs() * g.inner_radius_m * (g.outer_radius_m - g.inner_radius_m)
                / self.state.liquid.kinematic_viscosity_m2_s(),
        )?;
        if reynolds > self.definition.limits.maximum_gap_reynolds {
            return Err(Error::Budget("journal creeping-flow Reynolds limit"));
        }
        let torque = finite(-self.state.resistance_nm_s * rate)?;
        let power = finite(-torque * rate)?;
        if rate != 0.0 && (torque == 0.0 || power == 0.0) {
            return Err(Error::Invalid("nonzero journal work underflow"));
        }
        Ok((torque, power, reynolds))
    }

    fn update_rate(&mut self, rate: f64) -> Result<(), Error> {
        let (torque, power, reynolds) = self.response(rate)?;
        self.state.resisting_torque_nm = torque;
        self.state.dissipated_power_w = power;
        self.state.gap_reynolds = reynolds;
        Ok(())
    }

    pub fn refresh(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        velocity: &[f64],
    ) -> Result<(), Error> {
        let temperature = thermal.temperature_k(self.definition.heat_cell)?;
        let liquid = self
            .definition
            .liquid
            .at_temperature(Temperature::new(temperature))?;
        support(&self.definition.liquid, &self.state.liquid, &liquid)?;
        let old = self.state.liquid.dynamic_viscosity_pa_s();
        let change = finite((liquid.dynamic_viscosity_pa_s() - old).abs() / old)?;
        if change > self.definition.limits.maximum_relative_viscosity_change {
            return Err(Error::Budget(
                "journal viscosity change requires a smaller thermal substep",
            ));
        }
        let resistance = self.geometry_factor_m3 * liquid.dynamic_viscosity_pa_s();
        if !positive(resistance) {
            return Err(Error::Invalid(
                "journal material resistance is unrepresentable",
            ));
        }
        self.state.liquid = liquid;
        self.state.resistance_nm_s = resistance;
        self.update_rate(velocity[self.dof])
    }
}

// Ask the original source owner for continuous support of EVERY retained role;
// endpoint queries cannot jump over a hole or change the selected formulation.
fn support(
    law: &LiquidTemperatureLaw<'_>,
    from: &ResolvedLiquidState,
    to: &ResolvedLiquidState,
) -> Result<(), Error> {
    let from = from.material().query();
    let to = to.material().query();
    let axis = law.temperature_axis();
    let (low, high) = if from.axes()[axis] <= to.axes()[axis] {
        (&from, &to)
    } else {
        (&to, &from)
    };
    for property in law.reference_state().material().properties() {
        let selected = property.answer().receipt.selected;
        let claim = law
            .card()
            .claims()
            .claim(selected)
            .ok_or(Error::Invalid("journal liquid selected claim missing"))?;
        let path = law
            .card()
            .claims()
            .query_line_typed(
                &claim.key,
                axis,
                low,
                high,
                ClaimSelection::Pinned(selected),
            )
            .map_err(|source| MaterialStatePointError::Support {
                property: property.requirement().name().to_owned(),
                source: Box::new(source),
            })?;
        if path
            .knots
            .iter()
            .any(|(_, answer)| !positive(answer.evidence.value.value))
        {
            return Err(Error::Invalid(
                "journal liquid path contains an inadmissible coefficient",
            ));
        }
    }
    Ok(())
}
