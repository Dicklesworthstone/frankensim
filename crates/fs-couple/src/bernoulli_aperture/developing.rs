//! Quasi-steady laminar developing flow through a straight parallel-plate slit.
//!
//! A slit of gap `h` and streamwise length `L`, fed through a uniform inlet from a
//! plenum at gauge pressure `P` and exhausting at the ambient pressure, holds the
//! bulk speed `U` of `P = rho U^2/2 [1 + 4 f_app L/D_h]` with `D_h = 2h`: the
//! inlet kinetic head plus the developing-flow drop from the uniform inlet to the
//! exit. The apparent Fanning friction factor is Shah's correlation for parallel
//! plates (R. K. Shah, J. Fluids Eng. 100 (1978) 177-179, doi 10.1115/1.3448626;
//! Shah and London 1978), whose form and constants McCreery et al. state openly
//! (NURETH-15 2013, INL/CON-12-27833, eq. 2, in the Darcy form 4 f_app):
//! `f_app Re = 3.44/sqrt(x+) + (24 + 0.674/(4 x+) - 3.44/sqrt(x+)) / (1 + 0.000029/x+^2)`,
//! `x+ = L/(D_h Re)`, `Re = U D_h/nu`. No constant is tuned. As `x+` grows the
//! drop tends to the fully developed `12 mu L U/h^2` of
//! [`super::viscous::ViscousSlit`] plus the incremental `0.674 rho U^2/2`.
//!
//! While the walls' layers have not met, the inviscid core leaves at
//! `u_c = sqrt(2P/rho)`, so each wall displaces `delta* = (h/2)(1 - U/u_c)`.
//!
//! The held state is refused outside the correlation's validity, never
//! extrapolated ([`DevelopingSlitState::validity`]):
//! - `x+ < 4e-4`: the correlation agrees with Bodoia and Osterle's solution within
//!   2.5 % for `x+ > 0.0004` (and with Liu's within 2.4 %), as Shah and London state
//!   (secondary, search snippet only; primary not read).
//! - low Reynolds number: the developing-flow solutions depend on Re there; Chen
//!   (1973) gives `K(inf) = 0.64 + 38/Re` (secondary, search snippet of Shah and
//!   London's text only; primary not read). Refused where that inlet term departs
//!   from Shah's 0.674 by enough to move the held drop by more than the stated 2.4 %.
//! - laminarity: an undeveloped exit is refused above `Re_delta* = 520`, where the
//!   flat-plate (Blasius) layer turns linearly unstable to Tollmien-Schlichting waves,
//!   with `Re_delta*` built on the local displacement thickness and the free-stream
//!   speed (Manneville 2016, arXiv 1604.00840, sec. 2.3.4); here the core speed. A
//!   developed exit is refused above `Re_Dh = 1760`: plane Poiseuille flow sustains
//!   turbulence above about `Re_cl = 660` on the laminar centreline speed and the
//!   half gap (Kashyap, Duguet and Dauchot 2020, Entropy 22, 1001, sec. 3), and a
//!   laminar slit has `Re_Dh = (8/3) Re_cl`.
//!
//! It is quasi-steady: an attack or release is not covered. A converging windway
//! declared as a straight slit of its exit gap most likely overstates the loss.

use crate::acoustic_realize::AcousticRealizeError;
use fs_material::gas::GasState;

/// Smallest admitted `x+ = L/(D_h Re)`.
pub const MINIMUM_X_PLUS: f64 = 4e-4;
/// The correlation's stated accuracy (2.4 %), the largest admitted change of the
/// held drop by Chen's low-Reynolds-number inlet term.
pub const STATED_ACCURACY: f64 = 0.024;
/// Shah's incremental pressure drop of the fully developed limit, `K(inf)`.
pub const SHAH_K_INFINITY: f64 = 0.674;
/// Largest admitted `Re_delta*` of an undeveloped exit (Blasius linear instability).
pub const LAYER_CRITICAL_REYNOLDS_DISPLACEMENT: f64 = 520.0;
/// Largest admitted `Re_Dh` of a developed exit: `Re_cl = 660` as `(8/3) Re_cl`.
pub const CHANNEL_CRITICAL_REYNOLDS_HYDRAULIC: f64 = 660.0 * 8.0 / 3.0;

fn invalid(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription { what }
}

/// Shah's apparent Fanning friction factor times `Re` at `x+ = L/(D_h Re)`.
#[must_use]
pub fn shah_apparent_friction_re(x_plus: f64) -> f64 {
    let short = 3.44 / x_plus.sqrt();
    short + (24.0 + 0.674 / (4.0 * x_plus) - short) / (1.0 + 0.000_029 / (x_plus * x_plus))
}

/// A straight parallel-plate slit of explicit length bound to one source gas.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DevelopingSlit {
    length_m: f64,
    gas: GasState,
}

/// The held state at one plenum pressure and its validity measures.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DevelopingSlitState {
    /// Plenum gauge pressure [Pa].
    pub pressure_pa: f64,
    /// Bulk (mean) exit speed `U` [m/s].
    pub bulk_speed_m_s: f64,
    /// Volume flow `U w h` [m^3/s].
    pub flow_m3_s: f64,
    /// `x+ = L/(D_h Re)` at the exit.
    pub x_plus: f64,
    /// `Re = U D_h/nu`.
    pub reynolds_hydraulic: f64,
    /// `Re_L = u_c L/nu`, the wall layers' length Reynolds number at the exit.
    pub reynolds_length: f64,
    /// Inviscid core speed `u_c = sqrt(2P/rho)` [m/s].
    pub core_speed_m_s: f64,
    /// Displacement thickness of each wall layer, `(h/2)(1 - U/u_c)` [m].
    pub displacement_thickness_m: f64,
    /// `Re_delta* = u_c delta*/nu` at the exit.
    pub reynolds_displacement: f64,
    /// `|0.64 + 38/Re - 0.674| (rho U^2/2)/P`: the share of the held drop that
    /// Chen's low-Re inlet term would move.
    pub chen_inlet_shift: f64,
    /// Whether the core estimate has reached the Poiseuille centre speed
    /// (`u_c >= 1.5 U`), i.e. the wall layers have met before the exit.
    pub developed: bool,
}

impl DevelopingSlit {
    /// Bind the length to one coherent source gas (density and viscosity).
    ///
    /// # Errors
    /// An incoherent gas or a nonpositive/nonfinite length.
    pub fn from_gas(gas: &GasState, length_m: f64) -> Result<Self, AcousticRealizeError> {
        gas.validate_coherence().map_err(|_| invalid("developing slit requires one coherent source gas"))?;
        if !(length_m.is_finite() && length_m > 0.0) {
            return Err(invalid("developing slit requires a positive finite streamwise length"));
        }
        Ok(Self { length_m, gas: *gas })
    }

    /// Supplied streamwise length [m].
    #[must_use]
    pub const fn length_m(&self) -> f64 { self.length_m }

    /// Source gas whose density and viscosity the law uses.
    #[must_use]
    pub const fn gas(&self) -> &GasState { &self.gas }

    /// Plenum pressure holding bulk speed `U` through a gap `h` [Pa].
    fn pressure_at(&self, height_m: f64, speed_m_s: f64) -> f64 {
        if speed_m_s == 0.0 { return 0.0; }
        let (rho, mu, l) = (self.gas.density, self.gas.dynamic_viscosity, self.length_m);
        let x_plus = mu * l / (4.0 * rho * speed_m_s * height_m * height_m);
        // rho U^2/2 * 4 f_app L/D_h = mu L U (f_app Re)/(2 h^2).
        0.5 * rho * speed_m_s * speed_m_s
            + mu * l * speed_m_s * shah_apparent_friction_re(x_plus) / (2.0 * height_m * height_m)
    }

    /// Held state at plenum pressure `pressure_pa` through a `width_m` by
    /// `height_m` exit. The held pressure rises strictly with speed, so the bulk
    /// speed is the unique root, bracketed by `[0, sqrt(2P/rho)]` and bisected to
    /// adjacent floating-point values. Validity is reported, not enforced here.
    ///
    /// # Errors
    /// Nonpositive/nonfinite geometry, a negative or nonfinite pressure, or an
    /// unrepresentable result.
    pub fn held(&self, width_m: f64, height_m: f64, pressure_pa: f64)
        -> Result<DevelopingSlitState, AcousticRealizeError>
    {
        if !(width_m.is_finite() && width_m > 0.0 && height_m.is_finite() && height_m > 0.0) {
            return Err(invalid("developing slit requires positive finite width and gap"));
        }
        if !(pressure_pa.is_finite() && pressure_pa >= 0.0) {
            return Err(invalid("developing slit requires a finite nonnegative plenum pressure"));
        }
        let (rho, mu) = (self.gas.density, self.gas.dynamic_viscosity);
        let core = (2.0 * pressure_pa / rho).sqrt();
        if pressure_pa == 0.0 {
            return Ok(DevelopingSlitState { pressure_pa, bulk_speed_m_s: 0.0, flow_m3_s: 0.0,
                x_plus: f64::INFINITY, reynolds_hydraulic: 0.0, reynolds_length: 0.0,
                core_speed_m_s: 0.0, displacement_thickness_m: 0.5 * height_m,
                reynolds_displacement: 0.0, chen_inlet_shift: f64::INFINITY, developed: true });
        }
        let (mut lo, mut hi) = (0.0_f64, core);
        loop {
            let mid = 0.5 * (lo + hi);
            if mid <= lo || mid >= hi { break; }
            if self.pressure_at(height_m, mid) > pressure_pa { hi = mid; } else { lo = mid; }
        }
        let speed = if pressure_pa - self.pressure_at(height_m, lo)
            <= self.pressure_at(height_m, hi) - pressure_pa { lo } else { hi };
        let nu = mu / rho;
        let reynolds = speed * 2.0 * height_m / nu;
        let x_plus = self.length_m / (2.0 * height_m * reynolds);
        let displacement = 0.5 * height_m * (1.0 - speed / core);
        let head = 0.5 * rho * speed * speed;
        let state = DevelopingSlitState {
            pressure_pa,
            bulk_speed_m_s: speed,
            flow_m3_s: speed * width_m * height_m,
            x_plus,
            reynolds_hydraulic: reynolds,
            reynolds_length: core * self.length_m / nu,
            core_speed_m_s: core,
            displacement_thickness_m: displacement,
            reynolds_displacement: core * displacement / nu,
            chen_inlet_shift: (0.64 + 38.0 / reynolds - SHAH_K_INFINITY).abs() * head / pressure_pa,
            developed: core >= 1.5 * speed,
        };
        if ![state.bulk_speed_m_s, state.flow_m3_s, state.x_plus, state.reynolds_hydraulic,
            state.reynolds_length, state.displacement_thickness_m, state.reynolds_displacement,
            state.chen_inlet_shift].iter().all(|v| v.is_finite()) || !(speed > 0.0)
        {
            return Err(invalid("developing slit held state cannot be represented as finite physical values"));
        }
        Ok(state)
    }
}

impl DevelopingSlitState {
    /// Refuse a held state outside the correlation's validity (module docs).
    ///
    /// # Errors
    /// The first violated bound, by name.
    pub fn validity(&self) -> Result<(), &'static str> {
        if !(self.x_plus >= MINIMUM_X_PLUS) {
            return Err("developing slit x+ below 4e-4, where Shah's correlation is not stated");
        }
        if !(self.chen_inlet_shift <= STATED_ACCURACY) {
            return Err("developing slit Reynolds number too low: Chen's inlet term moves the held drop by more than 2.4 %");
        }
        if self.developed {
            if !(self.reynolds_hydraulic <= CHANNEL_CRITICAL_REYNOLDS_HYDRAULIC) {
                return Err("developing slit developed above Re_Dh 1760, where plane-channel turbulence can be sustained");
            }
        } else if !(self.reynolds_displacement <= LAYER_CRITICAL_REYNOLDS_DISPLACEMENT) {
            return Err("developing slit wall layer above Re_delta* 520, where the Blasius layer is linearly unstable");
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests;
