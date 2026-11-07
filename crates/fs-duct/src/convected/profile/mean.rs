use super::super::{FIRST_TRANSVERSE_ROOT, bad, positive};
use crate::DuctError;
use crate::fluid::FluidDuctError;
use crate::profile::RadiusProfile;
use fs_material::fluid_acoustics::FluidAcoustics;
use fs_material::gas::GasState;
use fs_math::det;

/// Minimum admitted `1-M²` for numerically inverted isentropic backgrounds.
/// This is 64*sqrt(f64::EPSILON), a precision margin from the singular sonic
/// branch, not a physical choking correction. The exact uniform owner is
/// unaffected. An admitted mass residual alone cannot resolve a sonic limit.
pub const MINIMUM_SUBSONIC_MARGIN: f64 = 1.0 / 1_048_576.0;

/// One steady section of the authored calorically-perfect isentropic flow.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IsentropicMeanState {
    /// Local circular area [m²].
    pub area_m2: f64,
    /// Static temperature [K].
    pub temperature_k: f64,
    /// Static pressure [Pa].
    pub pressure_pa: f64,
    /// Mass density [kg/m³].
    pub density_kg_m3: f64,
    /// Adiabatic sound speed [m/s].
    pub sound_speed_m_s: f64,
    /// Signed mean velocity along +x [m/s].
    pub velocity_m_s: f64,
    /// Signed subsonic Mach number U/c.
    pub mach: f64,
    /// Conjugate characteristic impedance c/(rho A) [J s/kg²].
    pub flow_impedance_j_s_kg2: f64,
}

/// Maintained, homentropic subsonic flow through the original smooth profile.
///
/// `new` explicitly selects a calorically-perfect closure: the original
/// `GasState`'s gamma and gas constant remain fixed over the computed range.
/// Its coherent numerical state and the complete geometry remain borrowed.
/// This is not a generic-fluid or variable-Cp EOS adapter; original source
/// support and the validity of frozen composition/Cp remain the caller's duty.
/// The state is authored at x=0 even when flow travels toward negative x.
#[derive(Debug)]
pub struct IsentropicMeanFlow<'source> {
    pub(super) source: &'source GasState,
    pub(super) profile: &'source RadiusProfile,
    velocity_at_zero: f64,
    area_at_zero: f64,
    mass_parameter: f64,
    mach_squared_at_zero: f64,
    exponent: f64,
    half_gamma_minus_one: f64,
    cutoff_bound: f64,
    mass_flow: f64,
}

impl<'source> IsentropicMeanFlow<'source> {
    /// Resolve mass conservation, constant entropy and stagnation enthalpy.
    /// Every cubic span's control hull participates in admission: an unseen
    /// throat cannot be skipped by choosing fewer acoustic integration steps.
    pub fn new(
        source: &'source GasState,
        profile: &'source RadiusProfile,
        velocity_at_zero_m_s: f64,
    ) -> Result<Self, FluidDuctError> {
        let state = source.acoustic_state()?;
        let c = state.sound_speed_m_s();
        if !velocity_at_zero_m_s.is_finite() || velocity_at_zero_m_s.abs() >= c {
            return Err(bad("isentropic flow must be finite and strictly subsonic at x=0").into());
        }
        let mach = velocity_at_zero_m_s.abs() / c;
        let q = 0.5 * (source.gamma - 1.0);
        let exponent = (source.gamma + 1.0) / (2.0 * (source.gamma - 1.0));
        let factor = 1.0 + q * mach * mach;
        let area = core::f64::consts::PI * profile.inlet_radius() * profile.inlet_radius();
        let mass_parameter = area * mach * det::exp(-exponent * log1p(q * mach * mach));
        let mass_flow = source.density * area * velocity_at_zero_m_s;
        if ![area, factor, exponent].into_iter().all(positive)
            || !mass_flow.is_finite()
            || (velocity_at_zero_m_s != 0.0 && (!positive(mass_parameter) || mass_flow == 0.0))
        {
            return Err(bad("isentropic inlet invariants are unrepresentable").into());
        }
        let mut flow = Self {
            source,
            profile,
            velocity_at_zero: velocity_at_zero_m_s,
            area_at_zero: area,
            mass_parameter,
            mach_squared_at_zero: mach * mach,
            exponent,
            half_gamma_minus_one: q,
            cutoff_bound: 0.0,
            mass_flow,
        };
        let [r_min, r_max] = profile.radius_bounds_m();
        let fastest = flow.at_radius(r_min)?;
        let slowest = flow.at_radius(r_max)?;
        // c sqrt(1-M²) decreases with |M| on the isentropic subsonic
        // branch. Combining its minimum with the maximum radius gives a
        // conservative lower bound at EVERY x, not merely sampled nodes.
        flow.cutoff_bound = FIRST_TRANSVERSE_ROOT
            * fastest.sound_speed_m_s
            * det::sqrt((1.0 - fastest.mach) * (1.0 + fastest.mach))
            / r_max;
        if !positive(flow.cutoff_bound) || !positive(slowest.flow_impedance_j_s_kg2) {
            return Err(bad("isentropic global acoustic bounds are unrepresentable").into());
        }
        Ok(flow)
    }

    /// Original coherent material state, without an invented replacement gas.
    pub const fn source(&self) -> &'source GasState {
        self.source
    }
    /// Original complete smooth geometry.
    pub const fn radius_profile(&self) -> &'source RadiusProfile {
        self.profile
    }
    /// Signed conserved steady mass flow [kg/s].
    pub const fn mass_flow_kg_s(&self) -> f64 {
        self.mass_flow
    }
    /// Conservative whole-profile first non-plane cut-on [rad/s].
    /// Control-hull admission may refuse a profile which tighter bounds admit.
    pub const fn transverse_cutoff_bound_rad_s(&self) -> f64 {
        self.cutoff_bound
    }
    /// Resolve the actual steady state at x [m]. Cold, bounded scalar work;
    /// prepared acoustic observations cache their local state instead.
    pub fn at_position(&self, x_m: f64) -> Result<IsentropicMeanState, DuctError> {
        self.at_radius(self.profile.radius(x_m)?)
    }

    fn mass_function(&self, mach: f64) -> f64 {
        mach * det::exp(-self.exponent * log1p(self.half_gamma_minus_one * mach * mach))
    }

    fn mach_at_area(&self, area: f64) -> Result<f64, DuctError> {
        if self.velocity_at_zero == 0.0 {
            return Ok(0.0);
        }
        if area == self.area_at_zero {
            return Ok(self.velocity_at_zero.abs() / self.source.sound_speed);
        }
        let target = self.mass_parameter / area;
        if !positive(target) || target >= self.mass_function(1.0) {
            return Err(bad(
                "isentropic profile reaches choking or an unrepresentable mass flux",
            ));
        }
        let (mut lo, mut hi, mut mach) = (0.0, 1.0, target);
        for _ in 0..80 {
            let value = self.mass_function(mach);
            let residual = value - target;
            if residual.abs() <= 4.0 * f64::EPSILON * target {
                return Ok(mach);
            }
            if residual < 0.0 {
                lo = mach;
            } else {
                hi = mach;
            }
            let derivative = (1.0 - mach * mach)
                * det::exp(-(self.exponent + 1.0) * log1p(self.half_gamma_minus_one * mach * mach));
            let newton = mach - residual / derivative;
            mach = if newton > lo && newton < hi {
                newton
            } else {
                lo + 0.5 * (hi - lo)
            };
        }
        Err(bad("isentropic subsonic mass-flow solve exhausted"))
    }

    fn at_radius(&self, radius: f64) -> Result<IsentropicMeanState, DuctError> {
        let area = core::f64::consts::PI * radius * radius;
        if !positive(area) {
            return Err(bad("isentropic profile area is unrepresentable"));
        }
        let mach = self.mach_at_area(area)?.copysign(self.velocity_at_zero);
        if (1.0 - mach) * (1.0 + mach) <= MINIMUM_SUBSONIC_MARGIN {
            return Err(bad(
                "isentropic profile is too close to the ill-conditioned sonic branch",
            ));
        }
        let temperature_delta = self.half_gamma_minus_one
            * (self.mach_squared_at_zero - mach * mach)
            / (1.0 + self.half_gamma_minus_one * mach * mach);
        let ratio = 1.0 + temperature_delta;
        let rho_ratio = det::exp(log1p(temperature_delta) / (self.source.gamma - 1.0));
        let density = self.source.density * rho_ratio;
        let sound_speed = self.source.sound_speed * det::sqrt(ratio);
        let state = IsentropicMeanState {
            area_m2: area,
            temperature_k: self.source.temperature * ratio,
            pressure_pa: self.source.pressure * rho_ratio * ratio,
            density_kg_m3: density,
            sound_speed_m_s: sound_speed,
            velocity_m_s: mach * sound_speed,
            mach,
            flow_impedance_j_s_kg2: (sound_speed / density) / area,
        };
        if ![
            state.temperature_k,
            state.pressure_pa,
            density,
            sound_speed,
            state.flow_impedance_j_s_kg2,
        ]
        .into_iter()
        .all(positive)
            || !state.velocity_m_s.is_finite()
            || mach.abs() >= 1.0
        {
            return Err(bad("isentropic local mean state is unrepresentable"));
        }
        Ok(state)
    }
}

// Compensate the rounding in 1+x before multiplying by 1/(gamma-1).
// This retains the finite isothermal limit as admitted gamma approaches one.
fn log1p(x: f64) -> f64 {
    let sum = 1.0 + x;
    if sum == 1.0 {
        x
    } else {
        det::ln(sum) * (x / (sum - 1.0))
    }
}
