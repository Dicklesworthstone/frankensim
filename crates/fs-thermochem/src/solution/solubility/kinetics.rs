//! Finite-rate exchange between a neutral dissolved component and its solid.
//!
//! With fixed temperature, solvent volume V and declared conductance G=kA,
//! `dn_d/dt = G (c_sat - n_d/V)`. The original pair law supplies c_sat;
//! its value is not a kinetic coefficient. The active trajectory is
//! `n_d(t) = n_d(0) + (V*c_sat-n_d(0))*(1-exp(-G*t/V))`.
//! If it exhausts the finite solid stock, the exact exhaustion event ends
//! dissolution. This is a fixed-area, well-mixed film approximation, not a
//! shrinking-particle, nucleation, or crystal-population model.
//!
//! Precipitation requires existing solid or an explicitly declared persistent
//! growth substrate. No component moles are supplied by that substrate. This
//! scalar law does not evolve temperature or infer heat transfer. The
//! fixed-area dissolution premise is the Noyes--Whitney model:
//! <https://doi.org/10.1021/ja02086a003>.

use super::{IdealSolubilityV1, SolubilityErrorV1, SolubilityPartitionV1};
use fs_math::det::{exp, expm1, ln};
use fs_qty::{Amount, Temperature, Time, Volume, VolumetricFlowRate};

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

/// Original solubility plus an explicit fixed-area transfer conductance [m³/s].
/// Numeric conductance does not establish experimental material applicability.
#[derive(Clone, Debug, PartialEq)]
pub struct SolubilityKineticsV1 {
    law: IdealSolubilityV1,
    conductance: VolumetricFlowRate,
    precipitation: PrecipitationPolicyV1,
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
        })
    }

    pub const fn law(&self) -> &IdealSolubilityV1 {
        &self.law
    }
    pub const fn conductance(&self) -> VolumetricFlowRate {
        self.conductance
    }
    pub const fn precipitation_policy(&self) -> PrecipitationPolicyV1 {
        self.precipitation
    }

    /// Allocation-free analytical update at fixed T, V and kA. The total is
    /// unchanged, no solid crosses an interface, and refusal changes no input.
    /// Every nonlinear source evaluation remains owned by the original law.
    pub fn advance(
        &self,
        temperature: Temperature,
        solvent_volume: Volume,
        initial: &SolubilityKineticStateV1,
        duration: Time,
    ) -> Result<SolubilityKineticStepV1, SolubilityErrorV1> {
        if !nonnegative(duration.value()) {
            return Err(SolubilityErrorV1::Invalid(
                "finite nonnegative kinetic duration",
            ));
        }
        let equilibrium =
            self.law
                .equilibrate(temperature, solvent_volume, initial.total_solute)?;
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
        let exponent = product_ratio(
            self.conductance.value(),
            duration.value(),
            solvent_volume.value(),
        )?;
        let fraction = -expm1(-exponent);
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
                    (solid - gap) + decayed(gap, exponent)
                };
                (dissolved + transfer, remaining, transfer)
            }
        } else {
            let surplus = -driving;
            let transfer = surplus * fraction;
            let remaining = if transfer <= 0.5 * dissolved {
                dissolved - transfer
            } else {
                capacity + decayed(surplus, exponent)
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

// A subnormal exp(-x) loses relative accuracy even when its scaled donor is
// normal. Splitting the exponent avoids that loss and the cancellation in
// ln(amount)-x; each multiplication can only decrease the finite donor.
fn decayed(amount: f64, exponent: f64) -> f64 {
    let fraction = exp(-exponent);
    let direct = amount * fraction;
    if fraction.is_normal() && direct > 0.0 {
        direct
    } else {
        let half = exp(-0.5 * exponent);
        (amount * half) * half
    }
}

#[cfg(test)]
mod tests;
