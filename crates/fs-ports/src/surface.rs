//! Material momentum, stress reaction and energy through a moving planar cut.
//!
//! For outward unit normal n, area A, material velocity u and cut velocity w,
//! the same evaluated state supplies
//!
//! - mdot = rho A (u-w).n;
//! - momentum_out = mdot u + A (p n - tau n);
//! - energy_out = mdot (h + |u|²/2 + gz) + p A w.n - A (tau n).u.
//!
//! The swept pressure term is necessary even at an impermeable moving wall,
//! where mass flow is zero. It is separate from the deviatoric stress work
//! already owned by the existing enthalpy chart. Source calorics, composition
//! and stress are inputs; no EOS, viscosity, wall law or additional pressure
//! force is inferred. The returned StreamPort uses the existing chart,
//! constituent-axis admission and complete conservation accounting.
//!
//! This operator evaluates one uniform face or quadrature point. Mesh motion,
//! volume/GCL integration, face quadrature, material-boundary conditions and
//! evolution of finite bodies remain with their original owners. It supports
//! advected constituents and entropy; heat and diffusive constituent fluxes
//! require their separate owners. A stress reaction is not the total momentum
//! flow of an open stream: captured or redirected matter requires its actual
//! receiving mechanics.

use fs_qty::{
    Area, Density, Force, Length, MassFlowRate, Power, Pressure, Qty, Velocity, VolumetricFlowRate,
};

use crate::{
    ChemicalEnergyAccounting, CoupleError, DeviatoricStressWork, EntropyFlowRate,
    MovingStreamEnthalpyChart, SpecificEnergy, StableId, StreamChartBinding, StreamConstituentFlow,
    StreamConstituentId, StreamEnergyChart, StreamKinematics, StreamPort,
    ensure_finite_stream_value,
};

/// Source constituent amount per unit material mass, mol/kg.
pub type SpecificAmount = Qty<0, -1, 0, 0, 0, 1>;
/// Source specific entropy, J/(kg K), retaining its original reference.
pub type SpecificEntropy = Qty<2, 0, -2, -1, 0, 0>;
/// Outward angular-momentum rate, N m, about an explicitly supplied origin.
pub type AngularMomentumFlow = Qty<2, 1, -2, 0, 0, 0>;

/// One oriented uniform face or area-weighted quadrature point, in the same
/// inertial Cartesian frame as the material velocity and stream binding.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StreamControlSurface {
    /// Finite nonnegative physical area; zero closes this face.
    pub area: Area,
    /// Finite unit vector pointing out of the material owner. Never normalized.
    pub outward_normal: [f64; 3],
    /// Velocity of this control-surface point, not necessarily material velocity.
    pub velocity: [Velocity; 3],
    /// Centroid or quadrature-point position used for angular-momentum flux.
    pub position: [Length; 3],
}

/// One original source evaluation at the face. Enthalpy and entropy retain the
/// source's reference and chemical ownership. Its owner admits EOS validity.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StreamSurfaceState {
    /// Positive finite material density.
    pub density: Density,
    /// Signed isotropic pressure in sigma = -p I + tau. No pressure datum is added.
    pub pressure: Pressure,
    /// Source specific enthalpy h; pressure-flow work is already included.
    pub specific_enthalpy: SpecificEnergy,
    /// Actual inertial material velocity and specific gravitational energy.
    pub kinematics: StreamKinematics,
    /// Source specific entropy; reference-dependent negative values are allowed.
    pub specific_entropy: SpecificEntropy,
    /// Original outward deviatoric/extra traction tau n, in pressure units.
    /// It EXCLUDES isotropic pressure. No no-slip or constitutive law is assumed.
    pub deviatoric_traction: [Pressure; 3],
}

/// Original source composition carried by this material stream.
#[derive(Debug, Clone, PartialEq)]
pub struct AdvectedConstituent {
    /// Species or element coordinate from the existing bound constituent axis.
    pub id: StreamConstituentId,
    /// Nonnegative mol/kg, evaluated by the original material owner.
    pub amount_per_mass: SpecificAmount,
}

/// Geometry-derived flux retained inside the existing enthalpy chart.
/// Private rates cannot be substituted independently after construction.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SurfaceTransport {
    surface: StreamControlSurface,
    relative_volume_flow: VolumetricFlowRate,
    swept_volume_flow: VolumetricFlowRate,
    mass_flow: MassFlowRate,
    advective_momentum_flow: [Force; 3],
    stress_momentum_flow: [Force; 3],
    momentum_flow: [Force; 3],
    pressure_motion_work: Power,
    stress_power_out_of_material: Power,
    stress_power_on_surface: Power,
}

impl SurfaceTransport {
    /// Original point, orientation, area and control-surface velocity.
    #[must_use]
    pub const fn surface(&self) -> StreamControlSurface {
        self.surface
    }

    /// A (u-w).n, signed out of the material owner, in m³/s.
    #[must_use]
    pub const fn relative_volume_flow(&self) -> VolumetricFlowRate {
        self.relative_volume_flow
    }

    /// A w.n, the geometric swept-volume rate for the original moving mesh.
    #[must_use]
    pub const fn swept_volume_flow(&self) -> VolumetricFlowRate {
        self.swept_volume_flow
    }

    /// Momentum carried by the crossing material, mdot u, in N.
    #[must_use]
    pub const fn advective_momentum_flow(&self) -> [Force; 3] {
        self.advective_momentum_flow
    }

    /// Outward non-advective momentum A (p n - tau n), in N.
    /// Its opposite acts on the material; this vector acts on an adjoining
    /// body only when the original boundary law couples that body to this face.
    #[must_use]
    pub const fn stress_momentum_flow(&self) -> [Force; 3] {
        self.stress_momentum_flow
    }

    /// The pressure work p A w.n absent from relative-flow enthalpy, in W.
    #[must_use]
    pub const fn pressure_motion_work(&self) -> Power {
        self.pressure_motion_work
    }

    /// A (p n - tau n).u, mechanical power leaving the material, in W.
    #[must_use]
    pub const fn stress_power_out_of_material(&self) -> Power {
        self.stress_power_out_of_material
    }

    /// A (p n - tau n).w, work on the moving cut, in W.
    /// It equals material-side stress power at a no-slip impermeable wall.
    /// At a slipping boundary their difference requires the original
    /// interfacial law; it is not silently deposited as heat.
    #[must_use]
    pub const fn stress_power_on_surface(&self) -> Power {
        self.stress_power_on_surface
    }

    /// Complete orbital angular-momentum rate r cross momentum_out.
    /// The origin is in the same inertial frame; no intrinsic couple stress or
    /// spin carried by the continuum is inferred from this Cauchy model.
    pub fn angular_momentum_flow_about(
        &self,
        origin: [Length; 3],
    ) -> Result<[AngularMomentumFlow; 3], CoupleError> {
        self.moment_about(origin, self.momentum_flow)
    }

    /// Moment of the stress reaction alone, excluding transported matter.
    /// Together with stress_momentum_flow it is the force/couple pair at origin.
    pub fn stress_moment_about(
        &self,
        origin: [Length; 3],
    ) -> Result<[AngularMomentumFlow; 3], CoupleError> {
        self.moment_about(origin, self.stress_momentum_flow)
    }

    fn moment_about(
        &self,
        origin: [Length; 3],
        force: [Force; 3],
    ) -> Result<[AngularMomentumFlow; 3], CoupleError> {
        let arm = core::array::from_fn(|i| self.surface.position[i].value() - origin[i].value());
        check_vector("angular_momentum_arm", arm)?;
        let f = force.map(Force::value);
        let moment = [
            arm[1].mul_add(f[2], -arm[2] * f[1]),
            arm[2].mul_add(f[0], -arm[0] * f[2]),
            arm[0].mul_add(f[1], -arm[1] * f[0]),
        ];
        check_vector("angular_momentum_flow", moment)?;
        Ok(moment.map(AngularMomentumFlow::new))
    }

    pub(super) fn require_rates(
        &self,
        mass: MassFlowRate,
        momentum: [Force; 3],
    ) -> Result<(), CoupleError> {
        if mass.value().to_bits() != self.mass_flow.value().to_bits() {
            return Err(CoupleError::StreamSurfaceFlowMismatch {
                field: "mass_flow",
                index: None,
            });
        }
        for (i, (actual, expected)) in momentum.into_iter().zip(self.momentum_flow).enumerate() {
            if actual.value().to_bits() != expected.value().to_bits() {
                return Err(CoupleError::StreamSurfaceFlowMismatch {
                    field: "momentum_flow",
                    index: Some(i),
                });
            }
        }
        Ok(())
    }
}

impl StreamPort {
    /// Evaluate one source-bound material stream through its actual moving cut.
    ///
    /// The resulting existing enthalpy chart owns all transported kinetic,
    /// gravitational and chemical energy, original deviatoric work, and the
    /// additional swept pressure work. Momentum and mass rates remain bound
    /// to these same inputs even when the chart is reused in Self::try_new.
    ///
    /// Constituents and entropy are advected with the same signed relative
    /// mass flow. There is no heat conduction or species diffusion in this
    /// constructor. The supplied stress operator/evidence IDs retain exactly
    /// the meaning of DeviatoricStressWork; this constructor adds no authority.
    ///
    /// # Errors
    ///
    /// Refuses nonfinite/unrepresentable arithmetic, nonpositive density,
    /// negative area/concentration, a nonunit normal, and the existing
    /// stream orientation, basis, chemical ownership and chart mismatches.
    #[allow(clippy::too_many_arguments)]
    pub fn try_from_uniform_surface(
        binding: StreamChartBinding,
        surface: StreamControlSurface,
        state: StreamSurfaceState,
        constituents: impl IntoIterator<Item = AdvectedConstituent>,
        chemical_energy: ChemicalEnergyAccounting,
        stress_operator_id: StableId,
        stress_evidence_id: StableId,
    ) -> Result<Self, CoupleError> {
        let (transport, deviatoric_power) = evaluate_surface(surface, state)?;
        let mut flows = Vec::new();
        for constituent in constituents {
            let concentration = constituent.amount_per_mass.value();
            ensure_finite_stream_value("constituent_amount_per_mass", None, concentration)?;
            if concentration < 0.0 {
                return Err(CoupleError::InvalidStreamSurface {
                    what: "advected constituent amount per mass must be nonnegative",
                });
            }
            flows.push(StreamConstituentFlow::try_new(
                constituent.id,
                crate::AmountFlowRate::new(transport.mass_flow.value() * concentration),
            )?);
        }
        let deviatoric = DeviatoricStressWork::try_new(
            &binding,
            deviatoric_power,
            stress_operator_id,
            stress_evidence_id,
        )?;
        let mut chart = MovingStreamEnthalpyChart::try_new(
            binding.clone(),
            state.specific_enthalpy,
            state.kinematics,
            deviatoric,
            chemical_energy,
        )?;
        chart.surface_transport = Some(transport);
        let energy = chart.energy_rate(transport.mass_flow);
        let entropy =
            EntropyFlowRate::new(transport.mass_flow.value() * state.specific_entropy.value());
        Self::try_new(
            binding,
            transport.mass_flow,
            flows,
            transport.momentum_flow,
            energy,
            entropy,
            StreamEnergyChart::MovingStreamEnthalpy(Box::new(chart)),
        )
    }
}

fn evaluate_surface(
    surface: StreamControlSurface,
    state: StreamSurfaceState,
) -> Result<(SurfaceTransport, Power), CoupleError> {
    for (field, value) in [
        ("surface_area", surface.area.value()),
        ("surface_density", state.density.value()),
        ("surface_pressure", state.pressure.value()),
        ("surface_enthalpy", state.specific_enthalpy.value()),
        ("surface_entropy", state.specific_entropy.value()),
    ] {
        ensure_finite_stream_value(field, None, value)?;
    }
    if state.density.value() <= 0.0 {
        return Err(CoupleError::NonPositiveStreamDensity);
    }
    if surface.area.value() < 0.0 {
        return Err(CoupleError::InvalidStreamSurface {
            what: "surface area must be nonnegative",
        });
    }
    check_vector("surface_normal", surface.outward_normal)?;
    let normal_squared = dot(surface.outward_normal, surface.outward_normal);
    if (normal_squared - 1.0).abs() > 64.0 * f64::EPSILON {
        return Err(CoupleError::InvalidStreamSurface {
            what: "surface normal must be unit length",
        });
    }
    let u = state.kinematics.velocity.map(Velocity::value);
    let w = surface.velocity.map(Velocity::value);
    check_vector("surface_velocity", w)?;
    check_vector("surface_position", surface.position.map(Length::value))?;
    let tau = state.deviatoric_traction.map(Pressure::value);
    check_vector("deviatoric_traction", tau)?;
    let relative = core::array::from_fn(|i| u[i] - w[i]);
    check_vector("surface_relative_velocity", relative)?;
    let area_normal = surface.outward_normal.map(|n| surface.area.value() * n);
    let relative_volume = dot(area_normal, relative);
    let swept_volume = dot(area_normal, w);
    let mass = state.density.value() * relative_volume;
    let advective = u.map(|v| mass * v);
    let extra_force = tau.map(|t| surface.area.value() * t);
    let stress = core::array::from_fn(|i| state.pressure.value() * area_normal[i] - extra_force[i]);
    let momentum = core::array::from_fn(|i| advective[i] + stress[i]);
    let pressure_motion = state.pressure.value() * swept_volume;
    let deviatoric_power = -dot(extra_force, u);
    let material_power = dot(stress, u);
    let surface_power = dot(stress, w);
    for (field, value) in [
        ("relative_volume_flow", relative_volume),
        ("swept_volume_flow", swept_volume),
        ("surface_mass_flow", mass),
        ("pressure_motion_work", pressure_motion),
        ("deviatoric_work", deviatoric_power),
        ("material_stress_power", material_power),
        ("surface_stress_power", surface_power),
    ] {
        ensure_finite_stream_value(field, None, value)?;
    }
    for (field, values) in [
        ("advective_momentum_flow", advective),
        ("stress_momentum_flow", stress),
        ("surface_momentum_flow", momentum),
    ] {
        check_vector(field, values)?;
    }
    if relative_volume != 0.0 && mass == 0.0 {
        return Err(CoupleError::InvalidStreamSurface {
            what: "nonzero relative volume flow lost its represented mass flux",
        });
    }
    Ok((
        SurfaceTransport {
            surface,
            relative_volume_flow: VolumetricFlowRate::new(relative_volume),
            swept_volume_flow: VolumetricFlowRate::new(swept_volume),
            mass_flow: MassFlowRate::new(mass),
            advective_momentum_flow: advective.map(Force::new),
            stress_momentum_flow: stress.map(Force::new),
            momentum_flow: momentum.map(Force::new),
            pressure_motion_work: Power::new(pressure_motion),
            stress_power_out_of_material: Power::new(material_power),
            stress_power_on_surface: Power::new(surface_power),
        },
        Power::new(deviatoric_power),
    ))
}

fn check_vector(field: &'static str, values: [f64; 3]) -> Result<(), CoupleError> {
    for (i, value) in values.into_iter().enumerate() {
        ensure_finite_stream_value(field, Some(i), value)?;
    }
    Ok(())
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0].mul_add(b[0], a[1].mul_add(b[1], a[2] * b[2]))
}

#[cfg(test)]
mod tests;
