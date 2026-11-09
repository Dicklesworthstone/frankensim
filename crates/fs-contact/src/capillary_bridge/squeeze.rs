//! Leading-order planar Reynolds squeeze flow on the original liquid bridge.
//!
//! For span b, gap h, invariant width w and signed separation rate v,
//! p_visc(x)=6 mu v (x²-b²/4)/h³, so F_visc=-mu w (b/h)³ v.
//! This is the parallel, impermeable, no-slip Newtonian Reynolds solution;
//! see Shim et al. (2019), section 3.5, equation (3.1):
//! https://pmc.ncbi.nlm.nih.gov/articles/PMC6528685/.
//!
//! The original fixed-volume circular geometry gives b=A/h+k h. Integrating
//! c(h)=mu w (A/h²+k)³ along the declared constant-rate path is analytic. The
//! positive divided differences below avoid endpoint subtraction even at rest.
//! F_mean=-integral(c dh)/dt, D=-F_mean*Delta_h, dF_mean/dh1=-c(h1)/dt.
//! Surface work and D therefore close without frozen-area or endpoint-force
//! approximations. These are exact identities of the selected reduced model,
//! not error bounds for the lubrication approximation.
//!
//! Meniscus-edge flow, contact-line wedge/slip loss, transient fluid inertia,
//! cavitation, non-Newtonian response and finite-width ends are not resolved.
//! The original mobile Young boundary remains an ideal zero-line-loss closure.
//! Bulk no slip applies away from its unresolved microscopic contact-line zone.
//! Thin-gap, slow-flow and pressure limits are caller-selected applicability
//! screens; topology transitions refuse instead of silently losing bulk loss.

use super::{
    BridgeError, BridgeEvent, BridgeState, SourcePlanarCapillaryBridge, finite, limit, positive,
};
use fs_blake3::DomainHasher;

/// Explicit no-slip lubrication and absolute-pressure operating envelope.
/// No pressure, cavitation threshold, wall slip or geometry scale is inferred.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlanarSqueezeFlowLimits {
    /// h/(b/2), with b the actual full wet span. Must be in (0,0.1].
    pub maximum_gap_to_half_span_ratio: f64,
    /// rho*U_lateral*h/mu, with U_lateral=|v|*b/(2h), in (0,0.1].
    pub maximum_reynolds_number: f64,
    /// rho*h²/(mu*dt), in (0,0.1]; screens unresolved momentum diffusion.
    pub maximum_unsteady_stokes_number: f64,
    /// Absolute pressure of the surrounding gas [Pa].
    pub ambient_pressure_pa: f64,
    /// Caller-supplied lower liquid pressure, above its cavitation limit [Pa].
    pub minimum_liquid_pressure_pa: f64,
    /// Upper absolute liquid pressure supported by the incompressible model [Pa].
    pub maximum_liquid_pressure_pa: f64,
}

#[derive(Clone, Copy, Default)]
pub(super) struct SqueezeStep {
    pub force_n: f64,
    pub tangent_n_per_m: f64,
    pub loss_j: f64,
}

impl SourcePlanarCapillaryBridge {
    /// Explicitly select leading-order no-slip bulk squeeze resistance.
    ///
    /// Call before initialization: this selection changes the model identity.
    /// Original source viscosity, volume, Young angles and surface energy remain
    /// authoritative. Source refresh preserves this selection and resolves the
    /// new viscosity. Pressure bounds include BOTH capillary and viscous loads.
    /// No-flow behavior remains the default; the separate dynamic-contact-line
    /// owner cannot silently consume a source with this bulk law enabled.
    pub fn with_no_slip_squeeze_flow(
        mut self,
        limits: PlanarSqueezeFlowLimits,
    ) -> Result<Self, BridgeError> {
        limits.validate()?;
        if self.squeeze_flow.is_some() {
            return Err(BridgeError::Invalid {
                field: "squeeze flow already selected",
            });
        }
        let mut hash = DomainHasher::new("org.frankensim.fs-contact.planar-reynolds-squeeze.v1");
        hash.update(self.identity.as_bytes());
        for value in limits.numbers() {
            hash.update(&value.to_bits().to_le_bytes());
        }
        self.identity = hash.finalize();
        self.squeeze_flow = Some(limits);
        Ok(self)
    }

    #[must_use]
    pub const fn squeeze_flow_limits(&self) -> Option<PlanarSqueezeFlowLimits> {
        self.squeeze_flow
    }

    pub(super) fn check_squeeze_geometry(&self, h: f64) -> Result<(), BridgeError> {
        let Some(limits) = self.squeeze_flow else {
            return Ok(());
        };
        let span_over_gap = (self.geometry.area / h) / h + self.geometry.segment;
        limit(
            2.0 / span_over_gap,
            limits.maximum_gap_to_half_span_ratio,
            "squeeze gap to half-span ratio",
        )?;
        self.check_squeeze_pressure(limits, h, h, 0.0)
    }

    pub(super) fn squeeze_step(
        &self,
        previous: BridgeState,
        h1: f64,
        event: BridgeEvent,
        dt: f64,
    ) -> Result<SqueezeStep, BridgeError> {
        let Some(limits) = self.squeeze_flow else {
            return Ok(SqueezeStep::default());
        };
        if event != BridgeEvent::None {
            return Err(BridgeError::Invalid {
                field: "squeeze flow requires a connected path without liquid topology change",
            });
        }
        if !previous.attached {
            return Ok(SqueezeStep::default());
        }
        let h0 = previous.gap_m;
        let dh = h1 - h0;
        let velocity = finite(dh / dt, "signed squeeze speed")?;
        let lower = h0.min(h1);
        let upper = h0.max(h1);
        let mu = self.liquid.dynamic_viscosity_pa_s();
        let rho = self.liquid.density_kg_m3();
        let g = self.geometry;
        // For b=A/h+k h, b/h is decreasing; b is convex. Thus endpoint
        // maxima bound the entire path, not only its sampled endpoints.
        let maximum_span = g
            .attached(lower)?
            .wetted_span_m
            .max(g.attached(upper)?.wetted_span_m);
        let lateral_speed = finite(
            0.5 * ((g.area / lower) / lower + g.segment) * velocity.abs(),
            "maximum squeeze lateral speed",
        )?;
        limit(
            rho * 0.5 * velocity.abs() * maximum_span / mu,
            limits.maximum_reynolds_number,
            "squeeze Reynolds number",
        )?;
        limit(
            (rho / mu) * upper * (upper / dt),
            limits.maximum_unsteady_stokes_number,
            "squeeze unsteady Stokes number",
        )?;
        limit(
            mu * lateral_speed / self.law.surface_tension_n_m(),
            self.limits.maximum_capillary_number,
            "squeeze lateral capillary number",
        )?;
        limit(
            rho * lateral_speed * lateral_speed * maximum_span.max(upper)
                / self.law.surface_tension_n_m(),
            self.limits.maximum_weber_number,
            "squeeze lateral Weber number",
        )?;
        self.check_squeeze_pressure(limits, lower, upper, velocity)?;
        let lower_aspect = (g.area / lower) / lower + g.segment;
        let maximum_resistance = positive(
            mu * g.width * lower_aspect * lower_aspect * lower_aspect,
            "maximum squeeze resistance",
        )?;
        // A conservative bound also covers an interior cancellation extremum
        // between capillary attraction and viscous compression on closing.
        limit(
            g.attached(lower)?.normal_force_n.abs() + maximum_resistance * velocity.abs(),
            self.limits.maximum_abs_force_n,
            "instantaneous capillary and squeeze force envelope",
        )?;
        let x = (g.area / h0) / h1;
        let k = g.segment;
        let ratio = h1 / h0;
        let inverse = h0 / h1;
        let mean_resistance = positive(
            mu * g.width
                * (k * k * k
                    + 3.0 * x * k * k
                    + x * x * k * (ratio + 1.0 + inverse)
                    + x * x * x * (ratio * ratio + ratio + 1.0 + inverse + inverse * inverse)
                        / 5.0),
            "positive integrated squeeze resistance",
        )?;
        let endpoint = (g.area / h1) / h1 + k;
        let endpoint_resistance = positive(
            mu * g.width * endpoint * endpoint * endpoint,
            "positive endpoint squeeze resistance",
        )?;
        let force = finite(-mean_resistance * velocity, "mean squeeze force")?;
        let loss = finite(-force * dh, "bulk squeeze dissipation")?;
        if loss < 0.0 || (dh != 0.0 && (force == 0.0 || loss == 0.0)) {
            return Err(BridgeError::Invalid {
                field: "representable positive squeeze dissipation",
            });
        }
        limit(
            loss,
            self.limits.maximum_abs_step_work_j,
            "bulk squeeze dissipation",
        )?;
        Ok(SqueezeStep {
            force_n: force,
            tangent_n_per_m: finite(-endpoint_resistance / dt, "squeeze tangent")?,
            loss_j: loss,
        })
    }

    fn check_squeeze_pressure(
        &self,
        limits: PlanarSqueezeFlowLimits,
        lower_gap: f64,
        upper_gap: f64,
        velocity: f64,
    ) -> Result<(), BridgeError> {
        // Along a constant-rate path, both suction and |p_visc(0)| reach
        // their maxima at the smallest gap. On closing, subtracting the two
        // pressures can have an interior extremum; independent extrema give
        // conservative whole-path bounds without missing that minimum.
        let aspect = (self.geometry.area / lower_gap) / lower_gap + self.geometry.segment;
        let viscous = finite(
            1.5 * self.liquid.dynamic_viscosity_pa_s() * velocity.abs() * aspect * aspect
                / lower_gap,
            "maximum central viscous pressure",
        )?;
        let minimum_edge = limits.ambient_pressure_pa
            + self
                .geometry
                .attached(lower_gap)?
                .liquid_minus_gas_pressure_pa;
        let maximum_edge = limits.ambient_pressure_pa
            + self
                .geometry
                .attached(upper_gap)?
                .liquid_minus_gas_pressure_pa;
        let minimum = finite(
            minimum_edge - if velocity > 0.0 { viscous } else { 0.0 },
            "minimum absolute squeeze pressure",
        )?;
        let maximum = finite(
            maximum_edge + if velocity < 0.0 { viscous } else { 0.0 },
            "maximum absolute squeeze pressure",
        )?;
        if minimum < limits.minimum_liquid_pressure_pa {
            return Err(BridgeError::Invalid {
                field: "squeeze cavitation pressure envelope",
            });
        }
        limit(
            maximum,
            limits.maximum_liquid_pressure_pa,
            "squeeze absolute pressure envelope",
        )
    }
}

impl PlanarSqueezeFlowLimits {
    fn numbers(self) -> [f64; 6] {
        [
            self.maximum_gap_to_half_span_ratio,
            self.maximum_reynolds_number,
            self.maximum_unsteady_stokes_number,
            self.ambient_pressure_pa,
            self.minimum_liquid_pressure_pa,
            self.maximum_liquid_pressure_pa,
        ]
    }
    fn validate(self) -> Result<(), BridgeError> {
        for value in self.numbers() {
            positive(value, "squeeze flow limit")?;
        }
        for value in [
            self.maximum_gap_to_half_span_ratio,
            self.maximum_reynolds_number,
            self.maximum_unsteady_stokes_number,
        ] {
            limit(value, 0.1, "thin creeping squeeze flow screen ceiling")?;
        }
        if self.minimum_liquid_pressure_pa >= self.ambient_pressure_pa
            || self.maximum_liquid_pressure_pa <= self.ambient_pressure_pa
        {
            return Err(BridgeError::Invalid {
                field: "ordered absolute squeeze pressure limits",
            });
        }
        Ok(())
    }
}
