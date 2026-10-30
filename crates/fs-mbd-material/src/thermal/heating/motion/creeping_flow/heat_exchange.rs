//! Conductive thermal reservoir on the original moving, expanding sphere.
//!
//! For an isothermal sphere in a constant-property unbounded medium, the
//! steady radial conduction solution gives G = 4*pi*k*r (Nu based on diameter
//! is 2). This is the infinite-exterior-radius limit of the spherical-shell
//! conductance 4*pi*k/(1/r_inner - 1/r_outer), independently checked in
//! fs-conduction's spherical-shell analytic test. No empirical convection
//! correlation, film-temperature rule or second thermal integrator is added.

use super::*;
use fs_material::fluid_transport::{ThermalFluidTransport, ThermalFluidTransportState};
use fs_material::thermal::transport::PreparedThermalTransport;

mod admission;
use admission::{admit_boundary, admit_solid, positive};

/// Explicit small-parameter allowances for the spherical conduction limit.
/// Each must be in (0,0.1]; these are model restrictions, not error certificates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SphereConductionLimits {
    /// Both translational and rotational Pe, plus expansion radial Pe.
    pub maximum_peclet: f64,
    /// Absolute diameter-based Rayleigh number, using actual gravity and
    /// source-owned expansion. Low translation speed alone does not admit it.
    pub maximum_rayleigh: f64,
    /// h*(V/A)/k_solid = k_fluid/(3*k_solid) for the uniform spherical cell.
    pub maximum_biot: f64,
    /// Exterior r²/alpha divided by the fastest caloric relaxation time;
    /// also bounds fractional temperature evolution over that diffusion time.
    pub maximum_diffusion_ratio: f64,
    /// |T_wall-T_fluid|/T_fluid for this frozen-property linearization.
    pub maximum_relative_temperature_contrast: f64,
}

/// Source conductivity for the SAME solid inventory and fixed source context.
/// Its density/Cp selections must agree with the original caloric sources.
/// Conductivity is explicitly treated as homogeneous isotropic bulk k; a
/// scalar directional measurement alone does not establish that assumption.
#[derive(Clone, Copy)]
pub struct SphereHeatExchange<'a> {
    pub solid: &'a PreparedThermalTransport,
    pub limits: SphereConductionLimits,
}

/// Instantaneous admitted exchange at the actual body temperature and radius.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SphereConductionSample {
    pub wall_temperature_k: f64,
    pub radius_m: f64,
    pub conductance_w_k: f64,
    /// Signed reservoir-to-body watts. Cooling is negative.
    pub heat_into_body_w: f64,
    pub translation_peclet: f64,
    pub rotation_peclet: f64,
    pub absolute_rayleigh: f64,
    pub biot: f64,
    pub thermal_diffusion_time_s: f64,
    pub thermal_relaxation_rate_per_s: f64,
    pub thermal_diffusion_ratio: f64,
}

/// One accepted reservoir exchange, thermal history and continuing trajectory.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HeatedCreepingThermalAudit {
    /// Existing storage, dynamics, fluid viscous loss and configuration work.
    pub coupled: HeatedCreepingAudit,
    /// Same source T, rho, cp, k and mu used for heat and motion in this window.
    pub fluid: ThermalFluidTransportState,
    pub before: SphereConductionSample,
    pub after: SphereConductionSample,
    /// Already included ONCE in coupled.thermal.external_heat_j. The prescribed
    /// fluid reservoir receives the opposite heat; its temperature is not solved.
    pub reservoir_heat_into_body_j: f64,
    /// Largest sampled values, including accepted expanded endpoints.
    pub maximum_peclet: f64,
    pub maximum_absolute_rayleigh: f64,
    pub maximum_biot: f64,
    pub maximum_thermal_diffusion_ratio: f64,
    /// Actual radius-change rate*r/alpha, separate from rigid translation/rotation.
    pub maximum_expansion_peclet: f64,
    /// r²/alpha * |dT_wall/dt|/T_fluid on each accepted interval.
    pub maximum_temperature_variation_ratio: f64,
}

impl<'a> HeatedRigidBody<'a> {
    /// Continue immersed-sphere motion while its actual material fluid heats or
    /// cools the SAME solid caloric cell. Source T and transport are bound once;
    /// no unrelated ambient temperature or hand-entered h is accepted.
    ///
    /// This uses the original observer source and additive decay-rate bound.
    /// Every accepted substep updates source expansion, momentum, current radius,
    /// drag and the next heat exchange. It admits only fully exposed spherical
    /// geometry, small Pe/Ra/Bi and slow thermal/geometry evolution. Authored
    /// Robin/reservoir/contact coverage refuses; independent heater watts remain.
    /// Ideal-gas beta=1/T is at fixed composition. A liquid with unknown beta
    /// requires zero gravity or exact thermal equilibrium, rechecked each step.
    ///
    /// The exterior temperature profile is assumed already quasi-steady; its
    /// initial diffusion transient, finite fluid heat storage, advective heat,
    /// radiation, phase change and reciprocal thermoelastic work are not modeled.
    /// Frozen fluid properties require source support at their stated point and
    /// the declared small temperature contrast; this is not a proof of accuracy
    /// for arbitrarily temperature-sensitive materials. Existing continuum,
    /// no-slip and hydrodynamic assumptions remain. Time refinement is required.
    pub fn advance_creeping_sphere_in_thermal_fluid<F: ThermalFluidTransport + ?Sized>(
        &self,
        sphere: &CreepingSphere<'_, F>,
        exchange: SphereHeatExchange<'_>,
        controls: HeatedCreepingControls,
        target_s: f64,
        config: StepConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, HeatedCreepingThermalAudit), HeatedMotionError> {
        let observer =
            self.thermal_creeping_observer(sphere, exchange, controls, &mut checkpoint)?;
        let before = observer.creeping.snapshot(&mut checkpoint)?;
        let (thermal, observer, audit) =
            self.thermal
                .advance_observed(&observer, target_s, config, &mut checkpoint)?;
        self.accept_thermal_creeping(before, thermal, observer, audit, &mut checkpoint)
    }

    /// Add an independent original heat/field owner on exactly these substeps.
    /// All proposed histories commit together. The added owner must not apply
    /// another heat-transfer law to the already occupied exterior surface.
    pub fn advance_creeping_sphere_in_thermal_fluid_with<
        O: DepositionObserver,
        F: ThermalFluidTransport + ?Sized,
    >(
        &self,
        observer: &O,
        sphere: &CreepingSphere<'_, F>,
        exchange: SphereHeatExchange<'_>,
        controls: HeatedCreepingControls,
        target_s: f64,
        config: StepConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<
        (Self, O, HeatedCreepingThermalAudit),
        ObserverCouplingError<O::Error, HeatedMotionError>,
    > {
        let motion = self
            .thermal_creeping_observer(sphere, exchange, controls, &mut checkpoint)
            .map_err(ObserverCouplingError::Second)?;
        let before = motion
            .creeping
            .snapshot(&mut checkpoint)
            .map_err(ObserverCouplingError::Second)?;
        let pair = CoupledObservers::new(observer.clone(), motion);
        let (thermal, pair, audit) =
            self.thermal
                .advance_observed(&pair, target_s, config, &mut checkpoint)?;
        let (observer, motion) = pair.into_parts();
        let (next, audit) = self
            .accept_thermal_creeping(before, thermal, motion, audit, &mut checkpoint)
            .map_err(ObserverCouplingError::Second)?;
        Ok((next, observer, audit))
    }

    fn thermal_creeping_observer<'m, 'f, 'k, F: ThermalFluidTransport + ?Sized>(
        &self,
        sphere: &'m CreepingSphere<'f, F>,
        exchange: SphereHeatExchange<'k>,
        controls: HeatedCreepingControls,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ThermalCreepingExpansion<'a, 'm, 'f, 'k, F>, HeatedMotionError> {
        let creeping = self.creeping_observer(sphere, controls, checkpoint)?;
        admit_solid(&creeping, exchange, checkpoint)?;
        let fluid = sphere
            .source()
            .thermal_transport_state()
            .map_err(crate::creeping_flow::CreepingFlowError::from)?;
        if fluid.transport() != sphere.transport() {
            return Err(HeatedMotionError::Invalid(
                "thermal and mechanical fluid snapshots disagree",
            ));
        }
        // Storage and mass are immutable on this path. Cache the original
        // caloric chart's conservative slope once, not on every source call.
        let region = &creeping.motion.expansion.regions[0];
        let inverse_capacity = positive(
            region.caloric.curve().maximum_temperature_slope()
                / self.thermal.cells()[&region.cell_id].mass_kg,
        )?;
        let initial = sample(
            &creeping,
            &self.thermal,
            exchange,
            fluid,
            inverse_capacity,
            checkpoint,
        )?;
        Ok(ThermalCreepingExpansion {
            creeping,
            exchange,
            fluid,
            inverse_capacity,
            initial,
            pending: initial,
            current: initial,
            heat_j: 0.0,
            maximum_peclet: initial.translation_peclet.max(initial.rotation_peclet),
            maximum_rayleigh: initial.absolute_rayleigh,
            maximum_biot: initial.biot,
            maximum_diffusion_ratio: initial.thermal_diffusion_ratio,
            maximum_expansion_peclet: 0.0,
            maximum_temperature_variation_ratio: 0.0,
        })
    }

    fn accept_thermal_creeping<F: ThermalFluidTransport + ?Sized>(
        &self,
        before: CreepingSnapshot,
        thermal: DepositionSimulation<'a>,
        mut observer: ThermalCreepingExpansion<'a, '_, '_, '_, F>,
        audit: EnergyAudit,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, HeatedCreepingThermalAudit), HeatedMotionError> {
        observer.update(&thermal, checkpoint)?;
        let (next, coupled) =
            self.accept_creeping(before, thermal, observer.creeping, audit, checkpoint)?;
        Ok((
            next,
            HeatedCreepingThermalAudit {
                coupled,
                fluid: observer.fluid,
                before: observer.initial,
                after: observer.current,
                reservoir_heat_into_body_j: observer.heat_j,
                maximum_peclet: observer.maximum_peclet,
                maximum_absolute_rayleigh: observer.maximum_rayleigh,
                maximum_biot: observer.maximum_biot,
                maximum_thermal_diffusion_ratio: observer.maximum_diffusion_ratio,
                maximum_expansion_peclet: observer.maximum_expansion_peclet,
                maximum_temperature_variation_ratio: observer.maximum_temperature_variation_ratio,
            },
        ))
    }
}

// The original CreepingExpansion owns every momentum/configuration update.
// This wrapper owns only its physical exterior heat source and exchange audit.
struct ThermalCreepingExpansion<'a, 'm, 'f, 'k, F: ThermalFluidTransport + ?Sized> {
    creeping: CreepingExpansion<'a, 'm, 'f, F>,
    exchange: SphereHeatExchange<'k>,
    fluid: ThermalFluidTransportState,
    inverse_capacity: f64,
    initial: SphereConductionSample,
    pending: SphereConductionSample,
    current: SphereConductionSample,
    heat_j: f64,
    maximum_peclet: f64,
    maximum_rayleigh: f64,
    maximum_biot: f64,
    maximum_diffusion_ratio: f64,
    maximum_expansion_peclet: f64,
    maximum_temperature_variation_ratio: f64,
}
impl<F: ThermalFluidTransport + ?Sized> Clone for ThermalCreepingExpansion<'_, '_, '_, '_, F> {
    fn clone(&self) -> Self {
        Self {
            creeping: self.creeping.clone(),
            ..*self
        }
    }
}
impl<F: ThermalFluidTransport + ?Sized> ThermalCreepingExpansion<'_, '_, '_, '_, F> {
    fn update(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), HeatedMotionError> {
        let next = sample(
            &self.creeping,
            thermal,
            self.exchange,
            self.fluid,
            self.inverse_capacity,
            checkpoint,
        )?;
        self.maximum_peclet = self
            .maximum_peclet
            .max(next.translation_peclet)
            .max(next.rotation_peclet);
        self.maximum_rayleigh = self.maximum_rayleigh.max(next.absolute_rayleigh);
        self.maximum_biot = self.maximum_biot.max(next.biot);
        self.maximum_diffusion_ratio = self
            .maximum_diffusion_ratio
            .max(next.thermal_diffusion_ratio);
        self.current = next;
        Ok(())
    }
}

impl<F: ThermalFluidTransport + ?Sized> DepositionObserver
    for ThermalCreepingExpansion<'_, '_, '_, '_, F>
{
    type Error = HeatedMotionError;
    fn begin(&mut self) -> Result<(), Self::Error> {
        self.creeping.begin()
    }
    fn event(
        &mut self,
        event: &Event,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<(), Self::Error> {
        self.creeping.event(event, thermal)?;
        self.update(thermal, &mut || Ok(()))
    }
    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.update(thermal, checkpoint)?;
        self.pending = self.current;
        let cell = self.creeping.motion.expansion.regions[0].cell_id;
        Ok(BTreeMap::from([(cell, self.pending.heat_into_body_w)]))
    }
    fn limit_timestep_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        proposed_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, Self::Error> {
        let proposed = self
            .creeping
            .limit_timestep_s(thermal, proposed_s, checkpoint)?;
        // The prior proposal already contains every earlier decay restriction.
        // Adding inverse stable times, rather than taking independent minima,
        // cannot double-spend that proposal's explicit thermal CFL allowance.
        positive(1.0 / (1.0 / proposed + self.pending.thermal_relaxation_rate_per_s / 0.8))
    }
    fn step(
        &mut self,
        start_s: f64,
        dt_s: f64,
        old: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.creeping
            .step(start_s, dt_s, old, thermal, checkpoint)?;
        self.update(thermal, checkpoint)?;
        let radius = self.pending.radius_m.max(self.current.radius_m);
        let alpha = self.fluid.transport().thermal_diffusivity_m2_s();
        let expansion_pe =
            finite((self.current.radius_m - self.pending.radius_m).abs() / dt_s * radius / alpha)?;
        let variation = finite(
            radius * radius / alpha
                * (self.current.wall_temperature_k - self.pending.wall_temperature_k).abs()
                / dt_s
                / self.fluid.temperature_k(),
        )?;
        if expansion_pe > self.exchange.limits.maximum_peclet
            || variation > self.exchange.limits.maximum_diffusion_ratio
        {
            return Err(HeatedMotionError::Invalid(
                "sphere expansion or temperature changes too rapidly for quasi-steady exterior conduction",
            ));
        }
        self.maximum_expansion_peclet = self.maximum_expansion_peclet.max(expansion_pe);
        self.maximum_temperature_variation_ratio =
            self.maximum_temperature_variation_ratio.max(variation);
        self.heat_j = finite(self.heat_j + finite(self.pending.heat_into_body_w * dt_s)?)?;
        checkpoint()?;
        Ok(())
    }
}

fn sample<F: FluidTransport + ?Sized>(
    creeping: &CreepingExpansion<'_, '_, '_, F>,
    thermal: &DepositionSimulation<'_>,
    exchange: SphereHeatExchange<'_>,
    fluid: ThermalFluidTransportState,
    inverse_capacity: f64,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<SphereConductionSample, HeatedMotionError> {
    checkpoint()?;
    let region = &creeping.motion.expansion.regions[0];
    admit_boundary(thermal, region.cell_id, creeping.sphere.radius_m())?;
    let wall = thermal.temperature_k(region.cell_id)?;
    let conductivity = exchange.solid.at_temperature(wall)?.conductivity_w_m_k();
    let mechanics = creeping.snapshot(checkpoint)?;
    let radius = mechanics.radius_m;
    let transport = fluid.transport();
    let contrast = finite(wall - fluid.temperature_k())?;
    let translation_peclet = finite(mechanics.load.translation_reynolds * transport.prandtl())?;
    let rotation_peclet = finite(mechanics.load.rotation_reynolds * transport.prandtl())?;
    let g = creeping.sphere.gravity().acceleration_world();
    let gravity = finite(g.x.hypot(g.y).hypot(g.z))?;
    let rayleigh = if gravity == 0.0 || contrast == 0.0 {
        0.0
    } else {
        let beta = fluid.volumetric_expansion_per_k().ok_or(HeatedMotionError::Invalid(
            "nonzero gravity and thermal contrast require the fluid source's volumetric expansion",
        ))?;
        transport
            .rayleigh(gravity, beta, contrast, 2.0 * radius)
            .map_err(crate::creeping_flow::CreepingFlowError::from)?
            .abs()
    };
    let biot = positive(transport.thermal_conductivity_w_m_k() / (3.0 * conductivity))?;
    let conductance =
        positive(4.0 * std::f64::consts::PI * transport.thermal_conductivity_w_m_k() * radius)?;
    let diffusion_time = positive(radius * radius / transport.thermal_diffusivity_m2_s())?;
    let relaxation_rate = positive(conductance * inverse_capacity)?;
    let diffusion_ratio = positive(diffusion_time * relaxation_rate)?;
    let limits = exchange.limits;
    if translation_peclet > limits.maximum_peclet
        || rotation_peclet > limits.maximum_peclet
        || rayleigh > limits.maximum_rayleigh
        || biot > limits.maximum_biot
        || diffusion_ratio > limits.maximum_diffusion_ratio
        || contrast.abs() / fluid.temperature_k() > limits.maximum_relative_temperature_contrast
    {
        return Err(HeatedMotionError::Invalid(
            "sphere heat exchange exceeds its Peclet, Rayleigh, Biot, thermal diffusion or temperature-contrast limit",
        ));
    }
    Ok(SphereConductionSample {
        wall_temperature_k: wall,
        radius_m: radius,
        conductance_w_k: conductance,
        heat_into_body_w: finite(-conductance * contrast)?,
        translation_peclet,
        rotation_peclet,
        absolute_rayleigh: rayleigh,
        biot,
        thermal_diffusion_time_s: diffusion_time,
        thermal_relaxation_rate_per_s: relaxation_rate,
        thermal_diffusion_ratio: diffusion_ratio,
    })
}
