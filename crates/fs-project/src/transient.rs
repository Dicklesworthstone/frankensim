//! Native thermal time grids and enthalpy charts shared by admission and execution.

use std::collections::BTreeSet;

use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use fs_qty::QtyAny;
use fs_scenario::Violation;

use crate::spec::{
    ConductionTransient, TransientInitialTemperature, TransientPhaseConductivity,
    TransientRegionEnthalpy, dims,
};

fn invalid(code: &'static str, what: impl Into<String>, fix: &str) -> Violation {
    Violation {
        code,
        what: what.into(),
        fix: fix.to_string(),
    }
}

fn work_limit() -> Violation {
    invalid(
        "project-conduction-transient-steps",
        "the event-aligned coarse grid and its nested half grid exceed the declared combined step cap",
        "increase max-steps to cover three times the coarse intervals, at most 10000; each workload switch consumes a time-grid boundary",
    )
}

impl TransientPhaseConductivity {
    /// Admit the declared mass-fraction law through the shared conduction owner.
    /// The multiplier acts on the existing base tensor at mean nodal temperature,
    /// evaluated at the element mean nodal liquid mass fraction. It does not
    /// change that tensor's source validity or create a material-query receipt.
    ///
    /// # Errors
    /// Refuses a missing source or nonpositive/nonfinite endpoint multipliers.
    pub fn model(
        &self,
    ) -> Result<fs_conduction::material::LiquidMassFractionConductivity, Violation> {
        if self.source.trim().is_empty() {
            return Err(invalid(
                "project-conduction-transient-phase-conductivity-source",
                "phase-dependent conductivity has no source or engineering basis",
                "state the basis for the declared phase multipliers separately from the base conductivity card",
            ));
        }
        fs_conduction::material::LiquidMassFractionConductivity::declared(
            self.solid_multiplier,
            self.liquid_multiplier,
        )
        .map_err(|error| {
            invalid(
                "project-conduction-transient-phase-conductivity",
                format!("phase-dependent conductivity refuses: {error}"),
                "declare finite strictly positive dimensionless solid and liquid multipliers",
            )
        })
    }
}

impl TransientRegionEnthalpy {
    /// Admit the explicit stationary-domain chart using fs-material's common
    /// equilibrium semantics. The reference density is fixed at every knot;
    /// this model transports heat without volume change or mass advection.
    ///
    /// The card identity is caller attribution. This method never authenticates
    /// a material card, resolves source claims or creates query receipts.
    ///
    /// # Errors
    /// Refuses missing source, wrong SI dimensions, nonpositive/nonfinite
    /// reference density, an unbounded knot count or invalid phase data.
    /// Initial-state consistency is checked separately by project admission.
    pub fn phase_curve(&self) -> Result<EquilibriumEnthalpyPhaseCurve, Violation> {
        if self.source.trim().is_empty() {
            return Err(invalid(
                "project-conduction-transient-enthalpy-source",
                format!("enthalpy model for `{}` has no source", self.region),
                "provide the engineering basis or source of this explicit enthalpy model",
            ));
        }
        if self.reference_density.dims != dims::DENSITY
            || !self.reference_density.value.is_finite()
            || self.reference_density.value <= 0.0
        {
            return Err(invalid(
                "project-conduction-transient-enthalpy-density",
                format!(
                    "enthalpy model for `{}` needs a positive finite reference density in kg/m³",
                    self.region
                ),
                "declare the stationary reference density; conductivity does not determine mass",
            ));
        }
        if !(2..=4096).contains(&self.knots.len()) {
            return Err(invalid(
                "project-conduction-transient-enthalpy-knots",
                format!(
                    "enthalpy model for `{}` carries {} knots; 2..=4096 are required",
                    self.region,
                    self.knots.len()
                ),
                "declare a bounded chart spanning the intended enthalpy range",
            ));
        }
        if self.knots.iter().any(|knot| {
            knot.specific_enthalpy.dims != dims::SPECIFIC_ENTHALPY
                || knot.temperature.dims != dims::TEMPERATURE
        }) {
            return Err(invalid(
                "project-conduction-transient-enthalpy-knot-dims",
                format!(
                    "enthalpy model for `{}` has a knot with incorrect units",
                    self.region
                ),
                "declare specific enthalpy in J/kg and absolute temperature in K",
            ));
        }
        let knots = self
            .knots
            .iter()
            .map(|knot| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: knot.specific_enthalpy.value,
                temperature_k: knot.temperature.value,
                liquid_mass_fraction: knot.liquid_mass_fraction,
                bulk_density_kg_m3: self.reference_density.value,
            })
            .collect();
        let curve = if self
            .knots
            .iter()
            .all(|knot| knot.liquid_mass_fraction == 0.0)
        {
            EquilibriumEnthalpyPhaseCurve::try_single_phase(
                self.material_card_identity,
                SolidLiquidPhase::Solid,
                knots,
            )
        } else if self
            .knots
            .iter()
            .all(|knot| knot.liquid_mass_fraction == 1.0)
        {
            EquilibriumEnthalpyPhaseCurve::try_single_phase(
                self.material_card_identity,
                SolidLiquidPhase::Liquid,
                knots,
            )
        } else {
            EquilibriumEnthalpyPhaseCurve::try_new(self.material_card_identity, knots)
        };
        curve.map_err(|error| invalid(
            "project-conduction-transient-enthalpy-curve",
            format!("enthalpy model for `{}` refuses: {error}", self.region),
            "supply finite monotone phase data and a nonzero card identity; a melting chart must span fully solid through fully liquid",
        ))
    }
}

impl TransientInitialTemperature {
    /// Return the authored uniform quantity, or `None` when every regional
    /// temperature must be resolved from its declared enthalpy state.
    #[must_use]
    pub const fn uniform_temperature(self) -> Option<QtyAny> {
        match self {
            Self::Uniform(temperature) => Some(temperature),
            Self::FromRegionalEnthalpy => None,
        }
    }
}

impl ConductionTransient {
    /// Build coarse endpoints including every declared workload switch.
    ///
    /// Each segment between switches is divided into
    /// `ceil(segment_duration / max_step)` intervals. The nested fine grid
    /// bisects each returned interval, so their combined work is three times
    /// this vector's length. With no schedules, this retains the original
    /// uniform grid's endpoint arithmetic. At most `max_steps / 3` endpoints
    /// and distinct switch times are allocated, regardless of input size.
    ///
    /// # Errors
    /// Refuses invalid time units, nonpositive or nonfinite controls, invalid
    /// power intervals, an insufficient combined work cap, or a coarse/fine
    /// interval whose positive duration cannot be represented. Region ownership
    /// and source coverage are additionally checked by `ProjectSpec::validate`.
    pub fn coarse_step_ends_s(&self) -> Result<Vec<f64>, Violation> {
        let horizon = self.horizon.value;
        let max_step = self.max_step.value;
        if self.horizon.dims != dims::TIME
            || self.max_step.dims != dims::TIME
            || !(horizon.is_finite() && horizon > 0.0 && max_step.is_finite() && max_step > 0.0)
        {
            return Err(invalid(
                "project-conduction-transient-quantity",
                "horizon and max-step must be positive finite quantities in seconds",
                "declare the physical horizon and maximum coarse step in coherent SI time units",
            ));
        }
        if !(3..=10_000).contains(&self.max_steps) {
            return Err(work_limit());
        }
        let cap = (self.max_steps / 3) as usize;
        // Positive finite f64 bit patterns have the same ordering as their
        // values. This set unifies simultaneous switches without tolerances
        // that could erase a short but explicitly declared pulse.
        let mut switches = BTreeSet::new();
        switches.insert(horizon.to_bits());
        for schedule in &self.power_schedules {
            if schedule.steps.is_empty() {
                return Err(invalid(
                    "project-conduction-transient-power-step",
                    format!("power schedule for `{}` has no intervals", schedule.region),
                    "declare one or more strictly increasing interval endpoints ending at the horizon",
                ));
            }
            // One schedule's strict endpoints are distinct, hence each needs
            // at least one coarse interval before any union is constructed.
            if schedule.steps.len() > cap {
                return Err(work_limit());
            }
            let mut previous = 0.0;
            for step in &schedule.steps {
                let until = step.until.value;
                if step.until.dims != dims::TIME
                    || !until.is_finite()
                    || until <= previous
                    || until > horizon
                    || step.watts.dims != dims::POWER
                    || !step.watts.value.is_finite()
                    || step.watts.value < 0.0
                {
                    return Err(invalid(
                        "project-conduction-transient-power-step",
                        format!(
                            "power schedule for `{}` needs increasing endpoints within the horizon and finite nonnegative delivered watts",
                            schedule.region
                        ),
                        "use seconds for until and watts for power; every interval must have positive duration",
                    ));
                }
                if switches.len() == cap && !switches.contains(&until.to_bits()) {
                    return Err(work_limit());
                }
                switches.insert(until.to_bits());
                previous = until;
            }
            if previous != horizon {
                return Err(invalid(
                    "project-conduction-transient-power-horizon",
                    format!(
                        "power schedule for `{}` ends at {previous} s, not the declared {horizon} s horizon",
                        schedule.region
                    ),
                    "end every regional workload exactly at the horizon; no tail power is inferred",
                ));
            }
        }
        let mut endpoints = Vec::with_capacity(cap);
        let mut start = 0.0;
        for bits in switches {
            let end = f64::from_bits(bits);
            let duration = end - start;
            let count = (duration / max_step).ceil().max(1.0);
            if !count.is_finite() || count > (cap - endpoints.len()) as f64 {
                return Err(work_limit());
            }
            let count = count as usize;
            let mut previous = start;
            for ordinal in 1..=count {
                let time = if ordinal == count {
                    end
                } else {
                    start + duration * (ordinal as f64 / count as f64)
                };
                let midpoint = previous + 0.5 * (time - previous);
                if !(time > previous && time <= end && midpoint > previous && midpoint < time) {
                    return Err(invalid(
                        "project-conduction-transient-time-resolution",
                        "a requested coarse interval or its half step has no representably positive duration",
                        "choose representable switch times and a time step that supports a nested half grid; no pulse can be skipped",
                    ));
                }
                endpoints.push(time);
                previous = time;
            }
            start = end;
        }
        Ok(endpoints)
    }
}
