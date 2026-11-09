//! Finite-rate exchange between a neutral dissolved component and its solid.
//!
//! With fixed temperature, solvent volume V and declared conductance G=kA,
//! `dn_d/dt = G (c_sat - n_d/V)`. The original pair law supplies c_sat;
//! its value is not a kinetic coefficient. The active trajectory is
//! `n_d(t) = n_d(0) + (V*c_sat-n_d(0))*(1-exp(-G*t/V))`.
//! Its saturation-deficit half-life is `ln(2)*V/G`: the magnitude of
//! `V*c_sat-n_d` halves during each interval of that length while exchange
//! remains active. This is not the time to dissolve half the initial solid.
//! If it exhausts the finite solid stock, the exact exhaustion event ends
//! dissolution. The default is a fixed-area, well-mixed film approximation.
//! [`SolubilityKineticsV1::with_particle_surface`] instead retains a fixed
//! population of self-similar particles: `G(s) = G_ref (s/s_ref)^(2/3)`.
//! Its finite-bath trajectory evolves the area throughout the step, including
//! crystal growth and disappearance, rather than freezing the initial area.
//!
//! Precipitation requires existing solid or an explicitly declared persistent
//! growth substrate. No component moles are supplied by that substrate. This
//! scalar law does not evolve temperature or infer heat transfer. The
//! fixed-area dissolution premise is the Noyes--Whitney model:
//! <https://doi.org/10.1021/ja02086a003>.
//! The self-similar surface premise is the Hixson--Crowell model:
//! <https://doi.org/10.1021/ie50260a018>. Neither mode supplies nucleation,
//! fragmentation, aggregation, a size distribution, or evolving film thickness.

mod particles;

use super::{IdealSolubilityV1, SolubilityErrorV1, SolubilityPartitionV1, partition_at_saturation};
use fs_math::decay::{attenuate, decayed_fraction};
use fs_math::det::{exp, expm1, ln, pow};
use fs_qty::{Amount, AmountConcentration, Temperature, Time, Volume, VolumetricFlowRate};

/// The caller, not supersaturation alone, supplies a surface for crystal growth.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PrecipitationPolicyV1 {
    /// No new solid can appear when the actual solid stock is zero.
    ExistingSolidOnly,
    /// A caller-declared persistent growth substrate remains after dissolution.
    /// It supplies no solute, and no nucleation induction time is modeled.
    SeededSurface,
}

/// Actual finite stocks; a supersaturated or undersaturated state is permitted.
/// Both phases remain explicit so a small donor survives a dominant receiver.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SolubilityKineticStateV1 {
    total_solute: Amount,
    dissolved: Amount,
    solid: Amount,
    balance_residual: Amount,
}
impl SolubilityKineticStateV1 {
    /// Declare both actual phases, retaining a trace solid or dissolved donor
    /// even when adding it to the dominant phase rounds the total unchanged.
    /// The rounded sum is the component total; the compensated residual keeps
    /// the discrepancy inspectable. This is also the transport boundary: only
    /// dissolved stock changes when a bare membrane retains the crystals.
    pub fn from_phases(dissolved: Amount, solid: Amount) -> Result<Self, SolubilityErrorV1> {
        Self::with_phases(
            Amount::new(dissolved.value() + solid.value()),
            dissolved.value(),
            solid.value(),
        )
    }

    /// Declare a complete component total and its initially dissolved part.
    pub fn new(total_solute: Amount, dissolved: Amount) -> Result<Self, SolubilityErrorV1> {
        if !nonnegative(total_solute.value())
            || !nonnegative(dissolved.value())
            || dissolved.value() > total_solute.value()
        {
            return Err(SolubilityErrorV1::Invalid("finite kinetic phase inventory"));
        }
        Self::with_phases(
            total_solute,
            dissolved.value(),
            total_solute.value() - dissolved.value(),
        )
    }

    fn with_phases(total: Amount, dissolved: f64, solid: f64) -> Result<Self, SolubilityErrorV1> {
        if !nonnegative(total.value()) || !nonnegative(dissolved) || !nonnegative(solid) {
            return Err(SolubilityErrorV1::Invalid(
                "finite nonnegative kinetic stocks",
            ));
        }
        let sum = dissolved + solid;
        let low = if dissolved >= solid {
            (dissolved - sum) + solid
        } else {
            (solid - sum) + dissolved
        };
        let residual = (sum - total.value()) + low;
        if !residual.is_finite() || residual.abs() > 8.0 * f64::EPSILON * total.value() {
            return Err(SolubilityErrorV1::Unrepresentable(
                "kinetic component balance",
            ));
        }
        Ok(Self {
            total_solute: total,
            dissolved: Amount::new(dissolved),
            solid: Amount::new(solid),
            balance_residual: Amount::new(residual),
        })
    }

    pub const fn total_solute(&self) -> Amount {
        self.total_solute
    }
    pub const fn dissolved(&self) -> Amount {
        self.dissolved
    }
    pub const fn solid(&self) -> Amount {
        self.solid
    }
    /// Compensated dissolved-plus-solid-minus-total [mol], not an error bound.
    pub const fn balance_residual(&self) -> Amount {
        self.balance_residual
    }
}

/// Original solubility plus an explicit transfer conductance [m³/s].
/// The default holds area fixed; an optional original particle stock defines
/// the reference area for self-similar dissolution and growth.
/// Numeric conductance does not establish experimental material applicability.
#[derive(Clone, Debug, PartialEq)]
pub struct SolubilityKineticsV1 {
    law: IdealSolubilityV1,
    conductance: VolumetricFlowRate,
    precipitation: PrecipitationPolicyV1,
    particle_reference_solid: Option<Amount>,
}
impl SolubilityKineticsV1 {
    pub fn new(
        law: IdealSolubilityV1,
        conductance: VolumetricFlowRate,
        precipitation: PrecipitationPolicyV1,
    ) -> Result<Self, SolubilityErrorV1> {
        if !nonnegative(conductance.value()) {
            return Err(SolubilityErrorV1::Invalid("nonnegative finite kA required"));
        }
        Ok(Self {
            law,
            conductance,
            precipitation,
            particle_reference_solid: None,
        })
    }

    pub const fn law(&self) -> &IdealSolubilityV1 {
        &self.law
    }
    /// Declared fixed-area conductance, or the reference conductance at
    /// [`Self::particle_reference_solid`] when particle evolution is selected.
    pub const fn conductance(&self) -> VolumetricFlowRate {
        self.conductance
    }
    pub const fn precipitation_policy(&self) -> PrecipitationPolicyV1 {
        self.precipitation
    }

    /// Select a fixed number of equally sized, self-similar particles with
    /// constant density and film mass-transfer coefficient. The declared
    /// conductance applies at `reference_solid`; actual conductance scales as
    /// `(solid/reference_solid)^(2/3)`. This original reference is retained
    /// across transport, temperature changes, cloning and split time steps.
    ///
    /// Only `ExistingSolidOnly` is admitted: when all particles disappear,
    /// their area is zero and they cannot regrow without a new nucleation law.
    /// A persistent inert seed has a different area/volume relation and cannot
    /// be represented by inventing additional component moles here. Solute
    /// volume remains negligible in the original incompressible solvent model.
    pub fn with_particle_surface(
        mut self,
        reference_solid: Amount,
    ) -> Result<Self, SolubilityErrorV1> {
        if !reference_solid.value().is_finite()
            || reference_solid.value() <= 0.0
            || self.precipitation != PrecipitationPolicyV1::ExistingSolidOnly
        {
            return Err(SolubilityErrorV1::Invalid(
                "particle surface requires a positive reference solid and existing-solid growth",
            ));
        }
        self.particle_reference_solid = Some(reference_solid);
        Ok(self)
    }

    /// Original stock at which the declared particle conductance applies.
    /// `None` denotes the unchanged fixed-area model.
    pub const fn particle_reference_solid(&self) -> Option<Amount> {
        self.particle_reference_solid
    }

    /// Instantaneous conductance at the actual solid stock, without allocation.
    /// Zero particle stock has zero exposed area. For fixed-area models this
    /// returns the original conductance; the growth policy is applied by
    /// [`Self::advance`], not by this coefficient query.
    pub fn conductance_at(&self, solid: Amount) -> Result<VolumetricFlowRate, SolubilityErrorV1> {
        if !nonnegative(solid.value()) {
            return Err(SolubilityErrorV1::Invalid("finite nonnegative solid stock"));
        }
        let Some(reference) = self.particle_reference_solid else {
            return Ok(self.conductance);
        };
        let conductance = self.conductance.value();
        if solid.value() == 0.0 || conductance == 0.0 {
            return Ok(VolumetricFlowRate::new(0.0));
        }
        if solid == reference {
            return Ok(self.conductance);
        }
        let ratio = solid.value() / reference.value();
        let direct = conductance * pow(ratio, 2.0 / 3.0);
        let value = if ratio.is_normal() && direct.is_normal() {
            direct
        } else {
            exp(ln(conductance) + (2.0 / 3.0) * (ln(solid.value()) - ln(reference.value())))
        };
        if !value.is_finite() || value <= 0.0 {
            return Err(SolubilityErrorV1::Unrepresentable(
                "positive particle conductance",
            ));
        }
        Ok(VolumetricFlowRate::new(value))
    }

    /// Time for an active saturation deficit or supersaturation to halve [s].
    ///
    /// At fixed temperature, volume and surface conductance, this is
    /// `ln(2)*V/G`. It describes the distance from saturation, not an arbitrary
    /// fraction of the initial solid stock. Finite solid exhaustion ends the
    /// trajectory, and a blocked precipitation branch does not relax at all.
    /// No finite half-life exists for zero conductance (`Ok(None)`). A positive
    /// conductance with an unrepresentable time scale returns an error.
    pub fn saturation_deficit_half_life(
        &self,
        solvent_volume: Volume,
    ) -> Result<Option<Time>, SolubilityErrorV1> {
        let volume = solvent_volume.value();
        if !volume.is_finite() || volume <= 0.0 {
            return Err(SolubilityErrorV1::Invalid(
                "finite positive solvent volume required",
            ));
        }
        if self.conductance.value() == 0.0 {
            return Ok(None);
        }
        if self.particle_reference_solid.is_some() {
            return Err(SolubilityErrorV1::Invalid(
                "evolving particle area has no constant saturation-deficit half-life",
            ));
        }
        let half_life = product_ratio(volume, core::f64::consts::LN_2, self.conductance.value())?;
        if !half_life.is_finite() {
            return Err(SolubilityErrorV1::Unrepresentable(
                "finite saturation-deficit half-life",
            ));
        }
        Ok(Some(Time::new(half_life)))
    }

    /// Allocation-free isothermal update at fixed T and V. Fixed-area kinetics
    /// use the exact exponential; particle kinetics solve the separable radius
    /// trajectory with bounded scalar work. The original reference area and
    /// film coefficient remain fixed, while the actual particle area evolves.
    /// The total is
    /// unchanged, no solid crosses an interface, and refusal changes no input.
    /// Every nonlinear source evaluation remains owned by the original law.
    pub fn advance(
        &self,
        temperature: Temperature,
        solvent_volume: Volume,
        initial: &SolubilityKineticStateV1,
        duration: Time,
    ) -> Result<SolubilityKineticStepV1, SolubilityErrorV1> {
        self.advance_at_saturation(
            temperature,
            solvent_volume,
            initial,
            duration,
            self.law.saturation_concentration(temperature)?,
        )
    }

    /// Advance the same original film/particle law toward an explicitly
    /// qualified conditional surface saturation [mol/m³ solvent].
    ///
    /// A coupled thermodynamic owner may compute this concentration from the
    /// original pure-pair calibration and the actual surrounding composition.
    /// It owns that activity calculation, its source support and chemical-work
    /// account. This numerical boundary supplies none of that qualification.
    /// The original law's temperature domain, rate coefficient, particle
    /// reference geometry and growth policy remain binding. The conditional
    /// saturation is held constant during this scalar stage; evolving other
    /// species requires an enclosing coupled-step accuracy controller.
    ///
    /// The exact fixed-area or bounded evolving-particle trajectory is reused
    /// without cloning or rewriting the original model. Returned equilibrium
    /// and exhaustion receipts refer to the supplied conditional saturation.
    /// No excess caloric model, temperature response or nucleation is inferred.
    pub fn advance_at_saturation(
        &self,
        temperature: Temperature,
        solvent_volume: Volume,
        initial: &SolubilityKineticStateV1,
        duration: Time,
        saturation: AmountConcentration,
    ) -> Result<SolubilityKineticStepV1, SolubilityErrorV1> {
        self.law.admit_temperature(temperature)?;
        if !nonnegative(duration.value()) {
            return Err(SolubilityErrorV1::Invalid(
                "finite nonnegative kinetic duration",
            ));
        }
        let equilibrium =
            partition_at_saturation(saturation, solvent_volume, initial.total_solute)?;
        let capacity = equilibrium.capacity.value();
        let dissolved = initial.dissolved.value();
        let solid = initial.solid.value();
        // Recover the saturation deficit from the small remaining solid when
        // the dissolved receiver dominates. In particular, capacity==total
        // gives deficit==solid even after dissolved rounds to capacity; later
        // steps must keep decaying that positive donor instead of freezing it.
        let excess_capacity = capacity - initial.total_solute.value();
        let driving = if excess_capacity >= 0.0 || dissolved > 0.5 * initial.total_solute.value() {
            excess_capacity + solid
        } else {
            capacity - dissolved
        };
        let blocked = driving < 0.0
            && solid == 0.0
            && self.precipitation == PrecipitationPolicyV1::ExistingSolidOnly;
        let mut result = SolubilityKineticStepV1 {
            state: *initial,
            equilibrium,
            transferred_to_solution: Amount::new(0.0),
            solid_exhaustion_time: None,
            precipitation_blocked: blocked,
        };
        if duration.value() == 0.0
            || self.conductance.value() == 0.0
            || driving == 0.0
            || blocked
            || (driving > 0.0 && solid == 0.0)
        {
            return Ok(result);
        }
        if let Some(reference_solid) = self.particle_reference_solid {
            return particles::advance(
                self,
                reference_solid,
                solvent_volume,
                initial,
                duration,
                result,
            );
        }
        let exponent = product_ratio(
            self.conductance.value(),
            duration.value(),
            solvent_volume.value(),
        )?;
        let fraction = decayed_fraction(exponent);
        if !(0.0..=1.0).contains(&fraction) {
            return Err(SolubilityErrorV1::Unrepresentable(
                "kinetic relaxation fraction",
            ));
        }
        let (next_dissolved, next_solid, transfer) = if driving > 0.0 {
            let gap = driving;
            let exhaustion = if excess_capacity > 0.0 {
                // ln(gap/(gap-solid)), without losing a small finite stock
                // against a large saturation capacity through 1-solid/gap.
                // The positive capacity-total difference determines the event;
                // gap-solid can round it away when the solid dominates.
                let ratio = solid / excess_capacity;
                let x = if ratio.is_finite() {
                    let shifted = 1.0 + ratio;
                    ln(shifted) - ((shifted - 1.0) - ratio) / shifted
                } else {
                    ln(solid) - ln(excess_capacity)
                };
                if !x.is_finite() || x <= 0.0 {
                    return Err(SolubilityErrorV1::Unrepresentable(
                        "solid exhaustion exponent",
                    ));
                }
                Some(x)
            } else {
                None
            };
            if let Some(event) = exhaustion.filter(|&event| exponent >= event) {
                let event_time =
                    product_ratio(solvent_volume.value(), event, self.conductance.value())?;
                if !event_time.is_finite() {
                    return Err(SolubilityErrorV1::Unrepresentable(
                        "finite solid exhaustion time",
                    ));
                }
                result.solid_exhaustion_time = Some(Time::new(event_time));
                (initial.total_solute.value(), 0.0, solid)
            } else {
                let transfer = gap * fraction;
                let remaining = if transfer <= 0.5 * solid {
                    solid - transfer
                } else if let Some(event) = exhaustion {
                    // Near the event, avoid subtracting a nearly exhausted
                    // donor. The event has not occurred, so this is positive.
                    excess_capacity * expm1(event - exponent)
                } else {
                    (solid - gap) + attenuate(gap, exponent)
                };
                (dissolved + transfer, remaining, transfer)
            }
        } else {
            let surplus = -driving;
            let transfer = surplus * fraction;
            let remaining = if transfer <= 0.5 * dissolved {
                dissolved - transfer
            } else {
                capacity + attenuate(surplus, exponent)
            };
            (remaining, solid + transfer, -transfer)
        };
        if transfer == 0.0 {
            return Err(SolubilityErrorV1::Unrepresentable(
                "nonzero kinetic phase transfer",
            ));
        }
        // A true asymptotic tail may round below the smallest positive float.
        // This numerical zero is not reported as a finite exhaustion event.
        result.state = SolubilityKineticStateV1::with_phases(
            initial.total_solute,
            next_dissolved,
            next_solid,
        )?;
        result.transferred_to_solution = Amount::new(transfer);
        Ok(result)
    }
}

/// Completed phase transfer and the separate instantaneous-equilibrium target.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SolubilityKineticStepV1 {
    state: SolubilityKineticStateV1,
    equilibrium: SolubilityPartitionV1,
    transferred_to_solution: Amount,
    solid_exhaustion_time: Option<Time>,
    precipitation_blocked: bool,
}
impl SolubilityKineticStepV1 {
    pub const fn state(&self) -> &SolubilityKineticStateV1 {
        &self.state
    }
    /// Target only: the actual state need not have reached this partition.
    pub const fn equilibrium(&self) -> &SolubilityPartitionV1 {
        &self.equilibrium
    }
    /// Integrated signed transfer [mol]; positive is solid to solution.
    pub const fn transferred_to_solution(&self) -> Amount {
        self.transferred_to_solution
    }
    /// Time from this step's start at which solid was exhausted, if reached.
    /// Surface conductance for dissolution becomes zero after that event.
    pub const fn solid_exhaustion_time(&self) -> Option<Time> {
        self.solid_exhaustion_time
    }
    pub const fn precipitation_blocked(&self) -> bool {
        self.precipitation_blocked
    }
}

fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}

// Prefer normal intermediate values, then fall back to representable subnormal
// arithmetic. No dimensionless relaxation time is silently replaced with zero.
fn product_ratio(a: f64, b: f64, c: f64) -> Result<f64, SolubilityErrorV1> {
    let candidates = [
        (a / c, (a / c) * b),
        (b / c, (b / c) * a),
        (a * b, (a * b) / c),
    ];
    candidates
        .iter()
        .find(|(first, result)| first.is_normal() && result.is_finite() && *result > 0.0)
        .or_else(|| {
            candidates
                .iter()
                .find(|(_, result)| result.is_finite() && *result > 0.0)
        })
        .map(|&(_, result)| result)
        .or_else(|| {
            candidates
                .iter()
                .any(|(_, result)| *result == f64::INFINITY)
                .then_some(f64::INFINITY)
        })
        .ok_or(SolubilityErrorV1::Unrepresentable(
            "positive kinetic time scale",
        ))
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod particle_tests;

#[cfg(test)]
mod conditional_tests;
