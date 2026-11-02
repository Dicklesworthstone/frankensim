use super::*;
use fs_qty::{Length, Power, Pressure, VolumetricFlowRate};

/// Fully developed, incompressible, homogeneous circular-pipe flow with no slip.
/// The caller declares a laminar wall-rate ceiling; this type does not infer a
/// transition criterion or a non-settling regime from rheology alone. Pressure
/// drop is the signed frictional driving difference, excluding elevation,
/// entrance, fitting, acceleration and pump-efficiency contributions.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LaminarCircularPipe {
    radius: Length,
    length: Length,
    maximum_laminar_wall_rate: Frequency,
}

/// Pressure-driven state. Positive flow and pressure drop use the same direction.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PipeFlow {
    /// Signed frictional pressure drop over the declared pipe length, in pascals.
    pub pressure_drop: Pressure,
    /// Signed mixture volume flow in cubic metres per second.
    pub volumetric_flow: VolumetricFlowRate,
    /// Signed wall shear stress in pascals, oriented with the pressure drop.
    pub wall_stress: Stress,
    /// Signed wall shear rate in inverse seconds; zero in a stationary pipe.
    pub wall_shear_rate: Frequency,
    /// Stationary full cross section below yield; moving unsheared core above yield.
    pub plug_radius: Length,
    /// Pressure times volume flow, always nonnegative; excludes pump losses.
    pub dissipated_power: Power,
}

/// An inverse result with its achieved signed flow residual and bounded work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PressureSolution {
    /// Achieved pressure-driven state inside the source and laminar support.
    pub flow: PipeFlow,
    /// Achieved minus requested signed volume flow, in cubic metres per second.
    pub flow_residual: VolumetricFlowRate,
    /// Bisection iterations used; zero for an accepted endpoint or zero-flow request.
    pub iterations: u32,
}

impl LaminarCircularPipe {
    /// Admit positive inner radius/length in metres and a caller-declared laminar
    /// wall-rate ceiling in inverse seconds; the ceiling is not inferred from rheology.
    pub fn new(
        radius: Length,
        length: Length,
        maximum_laminar_wall_rate: Frequency,
    ) -> Result<Self, SlurryError> {
        positive(radius.value(), "pipe radius must be positive")?;
        positive(length.value(), "pipe length must be positive")?;
        positive(
            maximum_laminar_wall_rate.value(),
            "laminar wall-rate ceiling must be positive",
        )?;
        Ok(Self {
            radius,
            length,
            maximum_laminar_wall_rate,
        })
    }

    /// Exact analytic laminar integral, including the yield threshold and plug.
    pub fn flow_at_pressure(
        &self,
        card: &HerschelBulkleyCard,
        state: SlurryConditions,
        pressure_drop: Pressure,
    ) -> Result<PipeFlow, SlurryError> {
        self.admit_source(card, state)?;
        let dp = pressure_drop.value();
        require(dp.is_finite(), "pressure drop must be finite")?;
        let wall = dp.abs() * (self.radius.value() / self.length.value()) * 0.5;
        require(
            wall.is_finite() && (dp == 0.0 || wall > 0.0),
            "wall stress is unrepresentable",
        )?;
        let rate = card.shear_rate_at_stress(Stress::new(wall), state)?.value();
        admit(
            "laminar wall shear rate",
            rate,
            [0.0, self.maximum_laminar_wall_rate.value()],
        )?;
        self.response(card, dp, wall, rate)
    }

    /// Invert monotone Q versus wall rate inside the source and laminar ceilings.
    /// Absolute flow tolerance is in m³/s; at most 128 iterations are allowed.
    /// Zero requested flow chooses zero pressure from the static pressure interval.
    /// A positive requested flow never returns a stationary result, even with a
    /// loose tolerance. The residual is numerical, not a source uncertainty bound.
    pub fn pressure_for_flow(
        &self,
        card: &HerschelBulkleyCard,
        state: SlurryConditions,
        requested: VolumetricFlowRate,
        absolute_flow_tolerance: VolumetricFlowRate,
        maximum_iterations: u32,
    ) -> Result<PressureSolution, SlurryError> {
        self.admit_source(card, state)?;
        let target = requested.value();
        require(target.is_finite(), "requested flow must be finite")?;
        let tolerance = positive(
            absolute_flow_tolerance.value(),
            "flow tolerance must be positive",
        )?;
        require(
            (1..=128).contains(&maximum_iterations),
            "inverse iteration budget must be 1..=128",
        )?;
        if target == 0.0 {
            return Ok(PressureSolution {
                flow: self.flow_at_pressure(card, state, Pressure::new(0.0))?,
                flow_residual: VolumetricFlowRate::new(0.0),
                iterations: 0,
            });
        }
        let mut low = 0.0;
        let mut high = self
            .maximum_laminar_wall_rate
            .value()
            .min(card.domain.shear_rate[1].value());
        let ceiling = self.flow_at_wall_rate(card, high, target)?;
        admit(
            "flow magnitude",
            target.abs(),
            [0.0, ceiling.volumetric_flow.value().abs()],
        )?;
        let residual = ceiling.volumetric_flow.value() - target;
        if residual.abs() <= tolerance {
            return Ok(PressureSolution {
                flow: ceiling,
                flow_residual: VolumetricFlowRate::new(residual),
                iterations: 0,
            });
        }
        for iterations in 1..=maximum_iterations {
            let rate = low + (high - low) * 0.5;
            if rate == low || rate == high {
                return Err(SlurryError::NoConvergence { iterations });
            }
            let flow = self.flow_at_wall_rate(card, rate, target)?;
            let residual = flow.volumetric_flow.value() - target;
            if residual.abs() <= tolerance {
                return Ok(PressureSolution {
                    flow,
                    flow_residual: VolumetricFlowRate::new(residual),
                    iterations,
                });
            }
            if flow.volumetric_flow.value().abs() < target.abs() {
                low = rate;
            } else {
                high = rate;
            }
        }
        Err(SlurryError::NoConvergence {
            iterations: maximum_iterations,
        })
    }

    fn admit_source(
        &self,
        card: &HerschelBulkleyCard,
        state: SlurryConditions,
    ) -> Result<(), SlurryError> {
        card.domain.admit_conditions(state)?;
        require(
            card.domain.shear_rate[0].value() == 0.0,
            "pipe requires source support over every rate from zero to the wall",
        )
    }

    fn flow_at_wall_rate(
        &self,
        card: &HerschelBulkleyCard,
        rate: f64,
        direction: f64,
    ) -> Result<PipeFlow, SlurryError> {
        let wall = card.stress_magnitude(rate)?;
        let dp = positive(
            (wall * (self.length.value() / self.radius.value())) * 2.0,
            "pressure drop is unrepresentable",
        )?
        .copysign(direction);
        self.response(card, dp, wall, rate)
    }

    fn response(
        &self,
        card: &HerschelBulkleyCard,
        dp: f64,
        wall: f64,
        rate: f64,
    ) -> Result<PipeFlow, SlurryError> {
        let r = self.radius.value();
        let y = card.parameters.yield_stress.value();
        let (q, plug) = if rate == 0.0 {
            (0.0, r)
        } else {
            let a = y / wall;
            // Compute the thin yielded annulus without cancellation in 1-a.
            let z = (wall - y) / wall;
            let m = 1.0 / card.parameters.flow_index.value();
            // Q = pi R³ integral_0^1 gamma(tau_wall*x) x² dx.
            // Substituting x=a+z*t leaves three elementary positive integrals.
            let factor = z * (z * z / (m + 3.0) + 2.0 * a * z / (m + 2.0) + a * a / (m + 1.0));
            let q = positive(
                core::f64::consts::PI * r * r * r * rate * factor,
                "volume flow is unrepresentable",
            )?;
            (q.copysign(dp), r * a)
        };
        let power = dp * q;
        require(
            power.is_finite() && power >= 0.0 && (q == 0.0 || power > 0.0),
            "hydraulic dissipation is unrepresentable",
        )?;
        Ok(PipeFlow {
            pressure_drop: Pressure::new(dp),
            volumetric_flow: VolumetricFlowRate::new(q),
            wall_stress: Stress::new(wall.copysign(dp)),
            wall_shear_rate: Frequency::new(rate.copysign(dp)),
            plug_radius: Length::new(plug),
            dissipated_power: Power::new(power),
        })
    }
}
