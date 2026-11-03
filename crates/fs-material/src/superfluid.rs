//! Local Landau two-fluid response with explicit thermodynamic inputs.
//!
//! Normal and superfluid densities share one material state. Entropy is carried
//! by the normal component; mass and kinetic energy involve both components.
//! Finite thermal expansion couples the two sound modes. This is neither an
//! ideal-gas helium model nor a zero-viscosity Newtonian-fluid adapter.
//!
//! References: Landau and Lifshitz, *Fluid Mechanics*, section 141 (two-fluid
//! sound); Donnelly and Barenghi, J. Phys. Chem. Ref. Data 27, 1217 (1998),
//! DOI: 10.1063/1.556028 (helium-II properties and conventions). The latter is
//! available at <https://srd.nist.gov/jpcrdreprint/1.556028.pdf>. No helium
//! property fit, phase boundary or experimental authority is inferred here.
//!
//! All properties are frozen at the caller's reference state. The model admits
//! positive normal and superfluid fractions; the zero-temperature and normal-
//! fluid endpoints need their own limiting models. Friction is a homogeneous
//! local source update with a declared coefficient, not vortex generation,
//! vortex-line evolution, wall friction, or a nonlinear thermal feedback solve.

use crate::MaterialError;
use crate::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_math::det;

/// Source-backed helium-II reference points with their original domain limits.
pub mod helium;

/// Explicit SI properties at ONE equilibrium state; no material-name defaults.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TwoFluidProperties {
    /// Absolute temperature [K], strictly positive.
    pub temperature_k: f64,
    /// Total mass density [kg/m³].
    pub density_kg_m3: f64,
    /// Superfluid mass fraction `rho_s / rho`, strictly between zero and one.
    pub superfluid_fraction: f64,
    /// Specific entropy per TOTAL mass [J/(kg K)].
    pub entropy_j_kg_k: f64,
    /// Isobaric heat capacity per total mass [J/(kg K)].
    pub specific_heat_cp_j_kg_k: f64,
    /// Adiabatic thermodynamic speed `sqrt((dp/drho)_s)` [m/s]. With nonzero
    /// expansion this is not either coupled two-fluid sound-mode speed.
    pub adiabatic_sound_speed_m_s: f64,
    /// Volumetric expansion `(1/V)(dV/dT)_p` [1/K], of either sign.
    pub thermal_expansion_per_k: f64,
    /// Normal-component dynamic shear viscosity [Pa s], nonnegative.
    pub normal_viscosity_pa_s: f64,
}

/// The two longitudinal, nondissipative Landau modes, sorted by speed.
/// Fast/slow labels avoid assuming which mode is density- or entropy-dominated.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SoundModes {
    /// Larger longitudinal phase speed [m/s].
    pub fast_m_s: f64,
    /// Smaller longitudinal phase speed [m/s].
    pub slow_m_s: f64,
}

/// Immutable checked state; no independently mutable derived coefficients.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TwoFluidState {
    properties: TwoFluidProperties,
    normal_fraction: f64,
    normal_density: f64,
    superfluid_density: f64,
    reduced_density: f64,
    kappa_s: f64,
    kappa_t: f64,
    cv: f64,
    isothermal_sound_speed: f64,
    sound: SoundModes,
}

/// Distinct normal and superfluid velocities, in the same inertial frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TwoFluidVelocity {
    /// Entropy-carrying component velocity [m/s].
    pub normal_m_s: [f64; 3],
    /// Superfluid component velocity [m/s].
    pub superfluid_m_s: [f64; 3],
}

/// Mass, entropy, and heat currents from the same pair of component velocities.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TwoFluidCurrents {
    /// Total mass current `rho_n v_n + rho_s v_s` [kg/(m² s)].
    pub mass_kg_m2_s: [f64; 3],
    /// Entropy current `rho s v_n` [W/(m² K)].
    pub entropy_w_m2_k: [f64; 3],
    /// Heat current RELATIVE TO the barycentric frame, `rho_s s T (v_n-v_s)`.
    /// This excludes advected enthalpy and pressure work [W/m²].
    pub heat_w_m2: [f64; 3],
}

/// Kinetic energy per total volume. Counterflow energy is already part of total.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct KineticEnergyDensity {
    /// Center-of-mass kinetic storage [J/m³].
    pub barycentric_j_m3: f64,
    /// Relative-motion kinetic storage [J/m³].
    pub counterflow_j_m3: f64,
    /// Sum of both kinetic stores [J/m³].
    pub total_j_m3: f64,
}

/// Linear, inviscid thermomechanical accelerations about the reference state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TwoFluidAcceleration {
    /// Normal-component acceleration [m/s²].
    pub normal_m_s2: [f64; 3],
    /// Superfluid-component acceleration [m/s²].
    pub superfluid_m_s2: [f64; 3],
}

/// Exact homogeneous friction source at frozen material properties.
/// The single dissipation account is removed from relative kinetic storage;
/// callers may transfer it to their thermal owner ONCE. Temperature is not
/// changed or a heat capacity integrated by this local velocity update.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MutualFrictionStep {
    /// Proposed endpoint velocities, in the input inertial frame.
    pub velocity: TwoFluidVelocity,
    /// Impulse per volume on the normal component [kg/(m² s)].
    pub normal_impulse_kg_m2_s: [f64; 3],
    /// Exact opposite of the normal-component impulse [kg/(m² s)].
    pub superfluid_impulse_kg_m2_s: [f64; 3],
    /// Physical heat released from relative kinetic energy [J/m³].
    pub dissipated_energy_j_m3: f64,
}

impl TwoFluidState {
    /// Admit a finite-temperature local two-fluid state. The thermodynamic
    /// identities are `kappa_T = kappa_s + alpha² T/(rho cp)` and
    /// `cv = cp kappa_s/kappa_T`. Neither entropy nor a phase fraction is
    /// inferred from temperature. Unrepresentable derived states refuse.
    pub fn try_new(p: TwoFluidProperties) -> Result<Self, MaterialError> {
        if ![
            p.temperature_k,
            p.density_kg_m3,
            p.entropy_j_kg_k,
            p.specific_heat_cp_j_kg_k,
            p.adiabatic_sound_speed_m_s,
        ]
        .into_iter()
        .all(positive)
            || !p.superfluid_fraction.is_finite()
            || p.superfluid_fraction <= 0.0
            || p.superfluid_fraction >= 1.0
            || !p.thermal_expansion_per_k.is_finite()
            || !nonnegative(p.normal_viscosity_pa_s)
        {
            return Err(invalid(
                "two-fluid T, rho, s, cp and adiabatic speed must be positive finite; require 0 < superfluid fraction < 1, finite expansion and nonnegative normal viscosity",
            ));
        }
        let normal_fraction = 1.0 - p.superfluid_fraction;
        let normal_density = p.density_kg_m3 * normal_fraction;
        let superfluid_density = p.density_kg_m3 * p.superfluid_fraction;
        let reduced_density = superfluid_density * normal_fraction;
        let c_ad2 = p.adiabatic_sound_speed_m_s * p.adiabatic_sound_speed_m_s;
        let kappa_s = (1.0 / p.density_kg_m3) / c_ad2;
        let expansion = p.thermal_expansion_per_k;
        let kappa_t = kappa_s
            + expansion * expansion * (p.temperature_k / p.density_kg_m3)
                / p.specific_heat_cp_j_kg_k;
        let cv = p.specific_heat_cp_j_kg_k * (kappa_s / kappa_t);
        let c_t2 = c_ad2 * (kappa_s / kappa_t);
        let b = (p.superfluid_fraction / normal_fraction)
            * (p.temperature_k / cv)
            * p.entropy_j_kg_k
            * p.entropy_j_kg_k;
        if ![
            normal_density,
            superfluid_density,
            reduced_density,
            c_ad2,
            kappa_s,
            kappa_t,
            cv,
            c_t2,
            b,
        ]
        .into_iter()
        .all(positive)
            || !nonnegative(p.normal_viscosity_pa_s / normal_density)
        {
            return Err(invalid(
                "derived two-fluid density, compressibility, heat capacity, viscosity or sound coefficient is unrepresentable",
            ));
        }
        // u^4 - (c_ad²+b)u² + c_T² b = 0. Scale before squaring, and recover
        // the small root from the product to avoid subtractive cancellation.
        let scale = c_ad2.max(b);
        let a_scaled = c_ad2 / scale;
        let b_scaled = b / scale;
        let c_scaled = c_t2 / scale;
        let discriminant = det::hypot(
            a_scaled - b_scaled,
            2.0 * (b_scaled * (a_scaled - c_scaled).max(0.0)).sqrt(),
        );
        let fast2 = (0.5 * (a_scaled + b_scaled + discriminant)) * scale;
        let slow2 = (c_t2 / fast2) * b;
        if !positive(fast2) || !positive(slow2) {
            return Err(invalid("two-fluid sound-mode roots are unrepresentable"));
        }
        Ok(Self {
            properties: p,
            normal_fraction,
            normal_density,
            superfluid_density,
            reduced_density,
            kappa_s,
            kappa_t,
            cv,
            isothermal_sound_speed: c_t2.sqrt(),
            sound: SoundModes {
                fast_m_s: fast2.sqrt(),
                slow_m_s: slow2.sqrt(),
            },
        })
    }

    /// Original supplied coefficients; changes to this copy do not alter the state.
    pub const fn properties(&self) -> TwoFluidProperties {
        self.properties
    }
    /// Absolute reference temperature [K].
    pub const fn temperature_k(&self) -> f64 {
        self.properties.temperature_k
    }
    /// Total reference density [kg/m³].
    pub const fn density_kg_m3(&self) -> f64 {
        self.properties.density_kg_m3
    }
    /// Superfluid mass fraction `rho_s/rho`.
    pub const fn superfluid_fraction(&self) -> f64 {
        self.properties.superfluid_fraction
    }
    /// Normal mass fraction `rho_n/rho`.
    pub const fn normal_fraction(&self) -> f64 {
        self.normal_fraction
    }
    /// Normal-component density per total volume [kg/m³].
    pub const fn normal_density_kg_m3(&self) -> f64 {
        self.normal_density
    }
    /// Superfluid density per total volume [kg/m³].
    pub const fn superfluid_density_kg_m3(&self) -> f64 {
        self.superfluid_density
    }
    /// Counterflow inertia `rho_n rho_s / rho` [kg/m³].
    pub const fn reduced_density_kg_m3(&self) -> f64 {
        self.reduced_density
    }
    /// Specific entropy per total mass [J/(kg K)].
    pub const fn entropy_j_kg_k(&self) -> f64 {
        self.properties.entropy_j_kg_k
    }
    /// Isobaric heat capacity per total mass [J/(kg K)].
    pub const fn specific_heat_cp_j_kg_k(&self) -> f64 {
        self.properties.specific_heat_cp_j_kg_k
    }
    /// Isochoric heat capacity from the same thermodynamics [J/(kg K)].
    pub const fn specific_heat_cv_j_kg_k(&self) -> f64 {
        self.cv
    }
    /// Thermodynamic `sqrt((dp/drho)_s)`, distinct from the coupled modes [m/s].
    pub const fn adiabatic_sound_speed_m_s(&self) -> f64 {
        self.properties.adiabatic_sound_speed_m_s
    }
    /// Thermodynamic `sqrt((dp/drho)_T)` [m/s].
    pub const fn isothermal_sound_speed_m_s(&self) -> f64 {
        self.isothermal_sound_speed
    }
    /// Adiabatic compressibility `1/(rho*c_ad²)` [1/Pa].
    pub const fn adiabatic_compressibility_per_pa(&self) -> f64 {
        self.kappa_s
    }
    /// Isothermal compressibility, including the supplied expansion [1/Pa].
    pub const fn isothermal_compressibility_per_pa(&self) -> f64 {
        self.kappa_t
    }
    /// Signed isobaric volumetric expansion [1/K].
    pub const fn thermal_expansion_per_k(&self) -> f64 {
        self.properties.thermal_expansion_per_k
    }
    /// Normal-component dynamic shear viscosity [Pa s].
    pub const fn normal_viscosity_pa_s(&self) -> f64 {
        self.properties.normal_viscosity_pa_s
    }
    /// Normal-component momentum diffusivity `mu_n/rho_n`, not `mu_n/rho`.
    pub fn normal_kinematic_viscosity_m2_s(&self) -> f64 {
        self.properties.normal_viscosity_pa_s / self.normal_density
    }
    /// Both finite-expansion longitudinal sound modes.
    pub const fn sound_modes(&self) -> SoundModes {
        self.sound
    }

    /// Component velocities from mass-weighted velocity and `w = v_n-v_s`.
    pub fn velocities_from_barycentric(
        &self,
        barycentric_m_s: [f64; 3],
        counterflow_m_s: [f64; 3],
    ) -> Result<TwoFluidVelocity, MaterialError> {
        finite_vector(barycentric_m_s)?;
        finite_vector(counterflow_m_s)?;
        let v = TwoFluidVelocity {
            normal_m_s: std::array::from_fn(|i| {
                barycentric_m_s[i] + self.superfluid_fraction() * counterflow_m_s[i]
            }),
            superfluid_m_s: std::array::from_fn(|i| {
                barycentric_m_s[i] - self.normal_fraction * counterflow_m_s[i]
            }),
        };
        validate_velocity(v)?;
        Ok(v)
    }

    /// Total mass current divided by total mass density [m/s].
    pub fn barycentric_velocity_m_s(
        &self,
        velocity: TwoFluidVelocity,
    ) -> Result<[f64; 3], MaterialError> {
        validate_velocity(velocity)?;
        finite_vector(std::array::from_fn(|i| {
            self.normal_fraction * velocity.normal_m_s[i]
                + self.superfluid_fraction() * velocity.superfluid_m_s[i]
        }))
    }

    /// Relative velocity `v_n-v_s`; positive heat flows in this direction.
    pub fn counterflow_velocity_m_s(
        &self,
        velocity: TwoFluidVelocity,
    ) -> Result<[f64; 3], MaterialError> {
        validate_velocity(velocity)?;
        finite_vector(std::array::from_fn(|i| {
            velocity.normal_m_s[i] - velocity.superfluid_m_s[i]
        }))
    }

    /// Evaluate mass, entropy and barycentric-relative heat flux in SI.
    pub fn currents(&self, velocity: TwoFluidVelocity) -> Result<TwoFluidCurrents, MaterialError> {
        let mean = self.barycentric_velocity_m_s(velocity)?;
        let relative = self.counterflow_velocity_m_s(velocity)?;
        let entropy_density = self.density_kg_m3() * self.entropy_j_kg_k();
        let heat_coefficient =
            self.superfluid_density * self.entropy_j_kg_k() * self.temperature_k();
        Ok(TwoFluidCurrents {
            mass_kg_m2_s: finite_vector(mean.map(|v| self.density_kg_m3() * v))?,
            entropy_w_m2_k: finite_vector(velocity.normal_m_s.map(|v| entropy_density * v))?,
            heat_w_m2: finite_vector(relative.map(|w| heat_coefficient * w))?,
        })
    }

    /// Separate mean and relative kinetic storage, with their total [J/m³].
    pub fn kinetic_energy_density(
        &self,
        velocity: TwoFluidVelocity,
    ) -> Result<KineticEnergyDensity, MaterialError> {
        let mean = speed(self.barycentric_velocity_m_s(velocity)?);
        let relative = speed(self.counterflow_velocity_m_s(velocity)?);
        let barycentric_j_m3 = 0.5 * self.density_kg_m3() * mean * mean;
        let counterflow_j_m3 = 0.5 * self.reduced_density * relative * relative;
        let total_j_m3 = barycentric_j_m3 + counterflow_j_m3;
        if !nonnegative(total_j_m3) {
            return Err(invalid(
                "two-fluid kinetic energy density is unrepresentable",
            ));
        }
        Ok(KineticEnergyDensity {
            barycentric_j_m3,
            counterflow_j_m3,
            total_j_m3,
        })
    }

    /// `a_s = -grad(p)/rho + s grad(T)` and
    /// `a_n = -grad(p)/rho - (rho_s/rho_n) s grad(T)`.
    /// These are local linear source terms; convection, normal viscous stress,
    /// vortices and body forces remain the spatial momentum owner's work.
    pub fn thermomechanical_acceleration(
        &self,
        pressure_gradient_pa_m: [f64; 3],
        temperature_gradient_k_m: [f64; 3],
    ) -> Result<TwoFluidAcceleration, MaterialError> {
        finite_vector(pressure_gradient_pa_m)?;
        finite_vector(temperature_gradient_k_m)?;
        let normal_entropy =
            (self.superfluid_fraction() / self.normal_fraction) * self.entropy_j_kg_k();
        Ok(TwoFluidAcceleration {
            normal_m_s2: finite_vector(std::array::from_fn(|i| {
                -pressure_gradient_pa_m[i] / self.density_kg_m3()
                    - normal_entropy * temperature_gradient_k_m[i]
            }))?,
            superfluid_m_s2: finite_vector(std::array::from_fn(|i| {
                -pressure_gradient_pa_m[i] / self.density_kg_m3()
                    + self.entropy_j_kg_k() * temperature_gradient_k_m[i]
            }))?,
        })
    }

    /// Exact `F_n = -gamma (v_n-v_s)` at frozen gamma [kg/(m³ s)].
    /// The total mass current is conserved, and relative kinetic energy falls
    /// exponentially. Gamma is explicit, nonnegative, and not inferred from mu.
    pub fn relax_linear_friction(
        &self,
        velocity: TwoFluidVelocity,
        gamma_kg_m3_s: f64,
        dt_s: f64,
    ) -> Result<MutualFrictionStep, MaterialError> {
        validate_velocity(velocity)?;
        friction_inputs(gamma_kg_m3_s, dt_s)?;
        if gamma_kg_m3_s == 0.0 || dt_s == 0.0 {
            return Ok(unchanged(velocity));
        }
        let rate_time = gamma_kg_m3_s * dt_s / self.reduced_density;
        if !nonnegative(rate_time) {
            return Err(invalid(
                "linear-friction rate times duration is unrepresentable",
            ));
        }
        let velocity_loss = -det::expm1(-rate_time);
        let energy_loss = -det::expm1(-2.0 * rate_time);
        self.friction_step(velocity, velocity_loss, energy_loss)
    }

    /// Exact isotropic Gorter-Mellink `F_n = -A rho_n rho_s |w|² w`,
    /// with explicit `A [m s/kg] >= 0` and `w = v_n-v_s`. This phenomenological
    /// turbulent-counterflow law has no implied onset or universal validity.
    /// `w(t) = w(0)/sqrt(1 + 2 A rho |w(0)|² t)`; no explicit-step stiffness.
    pub fn relax_gorter_mellink(
        &self,
        velocity: TwoFluidVelocity,
        a_m_s_kg: f64,
        dt_s: f64,
    ) -> Result<MutualFrictionStep, MaterialError> {
        validate_velocity(velocity)?;
        friction_inputs(a_m_s_kg, dt_s)?;
        if a_m_s_kg == 0.0 || dt_s == 0.0 {
            return Ok(unchanged(velocity));
        }
        let w = speed(self.counterflow_velocity_m_s(velocity)?);
        if w == 0.0 {
            return Ok(unchanged(velocity));
        }
        let x = 2.0 * a_m_s_kg * self.density_kg_m3() * w * w * dt_s;
        if !nonnegative(x) {
            return Err(invalid(
                "Gorter-Mellink rate times duration is unrepresentable",
            ));
        }
        let root = (1.0 + x).sqrt();
        let velocity_loss = if x < 1.0 {
            x / (root * (root + 1.0))
        } else {
            1.0 - 1.0 / root
        };
        let energy_loss = x / (1.0 + x);
        self.friction_step(velocity, velocity_loss, energy_loss)
    }

    fn friction_step(
        &self,
        velocity: TwoFluidVelocity,
        velocity_loss: f64,
        energy_loss: f64,
    ) -> Result<MutualFrictionStep, MaterialError> {
        let relative = self.counterflow_velocity_m_s(velocity)?;
        let kinetic = self.kinetic_energy_density(velocity)?;
        let impulse = finite_vector(relative.map(|w| -self.reduced_density * velocity_loss * w))?;
        let next = TwoFluidVelocity {
            normal_m_s: std::array::from_fn(|i| {
                velocity.normal_m_s[i] - self.superfluid_fraction() * velocity_loss * relative[i]
            }),
            superfluid_m_s: std::array::from_fn(|i| {
                velocity.superfluid_m_s[i] + self.normal_fraction * velocity_loss * relative[i]
            }),
        };
        validate_velocity(next)?;
        let dissipated_energy_j_m3 = kinetic.counterflow_j_m3 * energy_loss;
        if !nonnegative(dissipated_energy_j_m3) {
            return Err(invalid("mutual-friction dissipation is unrepresentable"));
        }
        Ok(MutualFrictionStep {
            velocity: next,
            normal_impulse_kg_m2_s: impulse,
            superfluid_impulse_kg_m2_s: impulse.map(|v| -v),
            dissipated_energy_j_m3,
        })
    }

    /// Explicit density-wave projection with `v_n = v_s` and zero expansion.
    /// This retains the original state but discards the independent entropy
    /// mode. Nonzero expansion refuses; ordinary thermoviscous wall or fluid
    /// transport traits are deliberately not implemented for either owner.
    pub fn first_sound_projection(&self) -> Result<FirstSoundProjection, MaterialError> {
        if self.thermal_expansion_per_k() != 0.0 {
            return Err(invalid(
                "a scalar in-phase first-sound projection requires explicitly zero thermal expansion; use the coupled two-fluid modes otherwise",
            ));
        }
        Ok(FirstSoundProjection {
            source: *self,
            acoustic: FluidAcousticState::try_new(
                self.density_kg_m3(),
                self.adiabatic_sound_speed_m_s(),
            )?,
        })
    }
}

/// Explicit, lossless density-wave approximation; never the full two-fluid law.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FirstSoundProjection {
    source: TwoFluidState,
    acoustic: FluidAcousticState,
}
impl FirstSoundProjection {
    /// Complete original admitted two-fluid state, including its discarded mode.
    pub const fn source(&self) -> &TwoFluidState {
        &self.source
    }
}
impl FluidAcoustics for FirstSoundProjection {
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        Ok(self.acoustic)
    }
}

fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn nonnegative(x: f64) -> bool {
    x.is_finite() && x >= 0.0
}
fn invalid(what: &str) -> MaterialError {
    MaterialError::Parameters { what: what.into() }
}
fn finite_vector(v: [f64; 3]) -> Result<[f64; 3], MaterialError> {
    if v.into_iter().all(f64::is_finite) {
        Ok(v)
    } else {
        Err(invalid(
            "two-fluid vector contains a nonfinite or unrepresentable component",
        ))
    }
}
fn validate_velocity(v: TwoFluidVelocity) -> Result<(), MaterialError> {
    finite_vector(v.normal_m_s)?;
    finite_vector(v.superfluid_m_s)?;
    Ok(())
}
fn speed(v: [f64; 3]) -> f64 {
    det::hypot(det::hypot(v[0], v[1]), v[2])
}
fn friction_inputs(coefficient: f64, time: f64) -> Result<(), MaterialError> {
    if nonnegative(coefficient) && nonnegative(time) {
        Ok(())
    } else {
        Err(invalid(
            "mutual-friction coefficient and physical duration must be finite and nonnegative",
        ))
    }
}
fn unchanged(velocity: TwoFluidVelocity) -> MutualFrictionStep {
    MutualFrictionStep {
        velocity,
        normal_impulse_kg_m2_s: [0.0; 3],
        superfluid_impulse_kg_m2_s: [0.0; 3],
        dissipated_energy_j_m3: 0.0,
    }
}
