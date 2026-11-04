//! Coherent work of the original finite radiators and passive wall patches.
//! Pressure, actual velocity and prescribed velocity are combined before any
//! quadratic product. No impedance matrix or additional acoustic solve is used.

use super::*;
use crate::ComplexSum;

/// Cycle-average work for peak phasors in the original source coordinates.
/// These are integrals of the retained discrete boundary fields, not error
/// bounds, a passivity certificate, or a replacement for far-field quadrature.
#[derive(Clone, Copy, Debug)]
pub struct IterativeBoundaryPower {
    pub omega_rad_s: f64,
    pub operator_hash: ContentHash,
    /// Work supplied through the original signed finite velocity traces [W].
    pub input_power_w: f64,
    /// Signed outward pressure/actual-normal-velocity work [W].
    /// Negative values are preserved, including discretization defects.
    pub outward_power_w: f64,
    /// Total work accepted by the original local passive wall laws [W].
    /// The historical name includes outgoing backing-fluid transmission;
    /// it is not an irreversible-heat account for a transmitting wall.
    pub absorbed_power_w: f64,
    /// Input minus outward work minus total accepted wall work [W].
    pub balance_w: f64,
    /// Integral of |p|(|v_drive|+|v_actual|)/2 plus accepted wall work [W].
    /// This remains meaningful when signed active powers nearly cancel.
    pub apparent_power_scale_w: f64,
    /// Absolute balance divided by the apparent scale; zero for a zero field.
    pub relative_balance: f64,
    /// Two source/panel visits, for scaling and coherent contraction.
    pub source_panel_products: u64,
}

impl IterativeBoundaryPower {
    /// Total power leaving the exterior through the passive wall laws [W].
    pub const fn accepted_power_w(&self) -> f64 {
        self.absorbed_power_w
    }
}

/// Coherent wall work split without replacing either the original boundary
/// flux or the material's independent positive-loss integrals.
#[derive(Clone, Copy, Debug)]
pub struct IterativeBoundaryPowerPartition {
    pub boundary_work: IterativeBoundaryPower,
    pub accepted_power_w: f64,
    /// Irreversible resistance/porous loss [W], not a transient temperature solve.
    pub irreversible_loss_power_w: f64,
    /// Power leaving into the explicitly lossless planar receiving fluid [W].
    /// This is not a second three-dimensional BEM radiation field.
    pub transmitted_power_w: f64,
    /// Accepted work minus independently evaluated loss and transmission [W].
    pub partition_residual_w: f64,
    /// Declared layers integrated across distinct transmitting wall laws.
    pub absorption_layers: usize,
}

struct BoundaryWorkEvaluation {
    work: IterativeBoundaryPower,
    partition: Option<[f64; 2]>,
    absorption_layers: usize,
}

impl IterativeBoundaryResponse<'_, '_> {
    /// Integrate coherent finite-source work in O(sources * panels), retaining
    /// the phase and cross terms that sums of per-source powers would discard.
    /// The original material owner supplies beta = rho*c/Z at this frequency;
    /// accepted work is integral Re(beta)*|p|²/(2*rho*c) dA. This cheap query
    /// does not integrate porous heat: use coherent_boundary_power_partition
    /// to distinguish irreversible loss from receiving-fluid transmission.
    /// Actual normal velocity
    /// is taken from the retained solve, so work closure is measured rather
    /// than imposed by reconstructing velocity from that wall law.
    ///
    /// `maximum_panel_products` admits both scaling and contraction passes
    /// before any material query. Temporary storage is one original material
    /// coefficient per panel plus one normalized rate per source. Cancellation
    /// returns no partial result. Compact-source banks refuse: their retained
    /// traces are scattered fields and they have no finite surface self port.
    ///
    /// Binary scaling protects large/tiny amplitudes and avoids squaring a
    /// source amplitude before restoring its physical units. Unrepresentable
    /// nonzero terms or final work refuse instead of silently becoming zero.
    pub fn coherent_boundary_power(
        &self,
        rates: &[C64],
        maximum_panel_products: u64,
        gate: &CancelGate,
    ) -> Result<IterativeBoundaryPower, IterativeObservationError> {
        Ok(self.evaluate_boundary_work(rates, maximum_panel_products, None, gate)?.work)
    }

    /// Separate outgoing transmission from irreversible loss using the same
    /// coherent TOTAL boundary pressure as the retained work calculation.
    /// The native material owner integrates positive viscous/thermal loss only
    /// on request. maximum_absorption_layers bounds declared layers summed over
    /// distinct transmitting laws, including identities; zero admits an empty
    /// matched-fluid backing. Source/panel and material budgets are checked
    /// before constitutive work. No power residual is assigned to either channel.
    pub fn coherent_boundary_power_partition(
        &self,
        rates: &[C64],
        maximum_panel_products: u64,
        maximum_absorption_layers: usize,
        gate: &CancelGate,
    ) -> Result<IterativeBoundaryPowerPartition, IterativeObservationError> {
        let evaluation = self.evaluate_boundary_work(
            rates, maximum_panel_products, Some(maximum_absorption_layers), gate,
        )?;
        let [irreversible_loss_power_w, transmitted_power_w] = evaluation.partition
            .ok_or_else(|| invalid("missing requested wall power partition"))?;
        let accepted_power_w = evaluation.work.accepted_power_w();
        let mut residual = ComplexSum::default();
        for value in [accepted_power_w, -irreversible_loss_power_w, -transmitted_power_w] {
            residual.add(C64::from_re(value));
        }
        let partition_residual_w = residual.finish().re;
        if !partition_residual_w.is_finite() {
            return Err(invalid("nonfinite wall power partition"));
        }
        poll(gate)?;
        Ok(IterativeBoundaryPowerPartition {
            boundary_work: evaluation.work,
            accepted_power_w,
            irreversible_loss_power_w,
            transmitted_power_w,
            partition_residual_w,
            absorption_layers: evaluation.absorption_layers,
        })
    }

    fn evaluate_boundary_work(
        &self,
        rates: &[C64],
        maximum_panel_products: u64,
        maximum_absorption_layers: Option<usize>,
        gate: &CancelGate,
    ) -> Result<BoundaryWorkEvaluation, IterativeObservationError> {
        poll(gate)?;
        let Excitation::Surface { owner, .. } = &self.observer.observation.excitation else {
            return Err(invalid(
                "coherent boundary work requires finite surface sources",
            ));
        };
        if rates.len() != self.sources().len() || rates.iter().any(|&v| !finite(v)) {
            return Err(invalid(
                "boundary work rates must match the finite source coordinates",
            ));
        }
        let products = (rates.len() as u64)
            .checked_mul(self.surface().areas().len() as u64)
            .and_then(|n| n.checked_mul(2))
            .ok_or_else(|| invalid("boundary work extent overflow"))?;
        if products > maximum_panel_products {
            return Err(invalid(
                "boundary work exceeds declared panel-product budget",
            ));
        }
        let rates: Vec<_> = rates
            .iter()
            .copied()
            .map(ScaledField::new)
            .collect::<Result<_, _>>()?;
        let boundary = &self.observer.observation.boundary;
        let (partition, absorption_layers) = if let Some(maximum) = maximum_absorption_layers {
            let (coefficients, layers) = boundary.power_coefficients(
                self.omega, maximum, || gate.is_requested(),
            )?;
            (Some(coefficients), layers)
        } else {
            (None, 0)
        };
        let beta = if let Some(coefficients) = &partition {
            coefficients.iter().map(|value| value.admittance).collect()
        } else {
            boundary.normalized_admittance(self.omega, || gate.is_requested())?
        };
        let rho_c = self.medium().density * self.medium().sound_speed;
        let mut sums = [ComplexSum::default(); 6];
        for (panel, &area) in self.surface().areas().iter().enumerate() {
            poll(gate)?;
            let mut exponents = [i32::MIN; 3];
            // Separate pressure/velocity scales do not mix their physical
            // units. Rate and trace are scaled independently, so a change of
            // source-coordinate magnitude does not create intermediate powers.
            for (source, rate) in rates.iter().enumerate() {
                if source % 256 == 0 {
                    poll(gate)?;
                }
                if rate.value == C64::ZERO {
                    continue;
                }
                let fields = [
                    self.solutions[source].pressure[panel],
                    self.solutions[source].velocity[panel],
                    owner.fields[source][panel],
                ];
                for (channel, field) in fields.into_iter().enumerate() {
                    if field != C64::ZERO {
                        exponents[channel] =
                            exponents[channel].max(exponent(field) + rate.exponent);
                    }
                }
            }
            let mut fields = [ComplexSum::default(); 3];
            for (source, rate) in rates.iter().enumerate() {
                if source % 256 == 0 {
                    poll(gate)?;
                }
                if rate.value == C64::ZERO {
                    continue;
                }
                let column = [
                    self.solutions[source].pressure[panel],
                    self.solutions[source].velocity[panel],
                    owner.fields[source][panel],
                ];
                for (channel, field) in column.into_iter().enumerate() {
                    if field == C64::ZERO {
                        continue;
                    }
                    let field = ScaledField::new(field)?;
                    let term = field.value * rate.value;
                    let shift = field.exponent + rate.exponent - exponents[channel];
                    fields[channel].add(C64::new(scale(term.re, shift)?, scale(term.im, shift)?));
                }
            }
            let mut combined = [ScaledField::default(); 3];
            for (channel, sum) in fields.into_iter().enumerate() {
                combined[channel] = ScaledField::new(sum.finish())?;
                if combined[channel].value != C64::ZERO {
                    combined[channel].exponent += exponents[channel];
                }
            }
            let [pressure, velocity, drive] = combined;
            let input = work(pressure, drive, area, false)?;
            let outward = work(pressure, velocity, area, false)?;
            let absorbed = product(
                pressure.value.norm_sq(),
                2 * pressure.exponent,
                &[0.5, area, beta[panel].re],
                rho_c,
            )?;
            sums[0].add(C64::from_re(input));
            sums[1].add(C64::from_re(outward));
            sums[2].add(C64::from_re(absorbed));
            sums[3].add(C64::from_re(work(pressure, drive, area, true)?));
            sums[3].add(C64::from_re(work(pressure, velocity, area, true)?));
            sums[3].add(C64::from_re(absorbed));
            if let Some(coefficients) = &partition {
                let coefficient = coefficients[panel];
                for (index, conductance) in [
                    coefficient.irreversible_conductance,
                    coefficient.transmitted_conductance,
                ].into_iter().enumerate() {
                    let value = product(
                        pressure.value.norm_sq(),
                        2 * pressure.exponent,
                        &[0.5, area, conductance],
                        rho_c,
                    )?;
                    sums[4 + index].add(C64::from_re(value));
                }
            }
        }
        let [
            input_power_w,
            outward_power_w,
            absorbed_power_w,
            apparent_power_scale_w,
            irreversible_loss_power_w,
            transmitted_power_w,
        ] = sums.map(|sum| sum.finish().re);
        let mut balance = ComplexSum::default();
        for value in [input_power_w, -outward_power_w, -absorbed_power_w] {
            balance.add(C64::from_re(value));
        }
        let balance_w = balance.finish().re;
        let relative_balance = if apparent_power_scale_w == 0.0 {
            0.0
        } else {
            balance_w.abs() / apparent_power_scale_w
        };
        if [
            input_power_w,
            outward_power_w,
            absorbed_power_w,
            apparent_power_scale_w,
            balance_w,
            relative_balance,
            irreversible_loss_power_w,
            transmitted_power_w,
        ]
        .iter()
        .any(|v| !v.is_finite())
        {
            return Err(invalid("nonfinite coherent boundary work"));
        }
        poll(gate)?;
        let work = IterativeBoundaryPower {
            omega_rad_s: self.omega,
            operator_hash: self.operator_hash(),
            input_power_w,
            outward_power_w,
            absorbed_power_w,
            balance_w,
            apparent_power_scale_w,
            relative_balance,
            source_panel_products: products,
        };
        Ok(BoundaryWorkEvaluation {
            work,
            partition: partition.map(|_| [irreversible_loss_power_w, transmitted_power_w]),
            absorption_layers,
        })
    }
}

fn invalid(message: &'static str) -> IterativeObservationError {
    TransferError::Invalid(message).into()
}

#[derive(Clone, Copy, Default)]
struct ScaledField {
    value: C64,
    exponent: i32,
}
impl ScaledField {
    fn new(value: C64) -> Result<Self, IterativeObservationError> {
        if !finite(value) {
            return Err(invalid("nonfinite coherent boundary field"));
        }
        if value == C64::ZERO {
            return Ok(Self::default());
        }
        let exponent = exponent(value);
        Ok(Self {
            value: C64::new(scale(value.re, -exponent)?, scale(value.im, -exponent)?),
            exponent,
        })
    }
}

// The exponent of the largest component, including subnormal binary64 values.
fn exponent(value: C64) -> i32 {
    let magnitude = value.re.abs().max(value.im.abs());
    let bits = magnitude.to_bits();
    let normal = ((bits >> 52) & 0x7ff) as i32;
    if normal == 0 {
        63 - bits.leading_zeros() as i32 - 1074
    } else {
        normal - 1023
    }
}

fn scale(mut value: f64, mut exponent: i32) -> Result<f64, IterativeObservationError> {
    if value == 0.0 {
        return Ok(value);
    }
    while exponent > 1023 {
        value *= f64::from_bits(2046_u64 << 52);
        exponent -= 1023;
    }
    while exponent < -1022 {
        value *= f64::MIN_POSITIVE;
        exponent += 1022;
    }
    value *= f64::from_bits(((exponent + 1023) as u64) << 52);
    if !value.is_finite() || value == 0.0 {
        Err(invalid(
            "coherent boundary work exceeds representable numerical range",
        ))
    } else {
        Ok(value)
    }
}

// Factor before restoring the binary exponent. Neither tiny/large areas nor
// rates are squared or multiplied into a physical intermediate prematurely.
fn product(
    value: f64,
    mut exponent: i32,
    factors: &[f64],
    divisor: f64,
) -> Result<f64, IterativeObservationError> {
    if value == 0.0 || factors.contains(&0.0) {
        return Ok(0.0);
    }
    let mut normalized = ScaledField::new(C64::from_re(value))?;
    exponent += normalized.exponent;
    for &factor in factors {
        let factor = ScaledField::new(C64::from_re(factor))?;
        normalized.value.re *= factor.value.re;
        exponent += factor.exponent;
    }
    let divisor = ScaledField::new(C64::from_re(divisor))?;
    normalized.value.re /= divisor.value.re;
    exponent -= divisor.exponent;
    scale(normalized.value.re, exponent)
}

fn work(
    pressure: ScaledField,
    velocity: ScaledField,
    area: f64,
    apparent: bool,
) -> Result<f64, IterativeObservationError> {
    let value = if apparent {
        pressure.value.abs() * velocity.value.abs()
    } else {
        let mut dot = ComplexSum::default();
        for (p, v) in [
            (pressure.value.re, velocity.value.re),
            (pressure.value.im, velocity.value.im),
        ] {
            let (high, low) = fs_math::eft::two_prod(p, v);
            dot.add(C64::from_re(high));
            dot.add(C64::from_re(low));
        }
        dot.finish().re
    };
    product(
        value,
        pressure.exponent + velocity.exponent,
        &[0.5, area],
        1.0,
    )
}
