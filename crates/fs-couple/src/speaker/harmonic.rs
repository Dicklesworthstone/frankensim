//! Continuous coil/cone/cabinet response with full baffled-piston radiation.
//!
//! This is the frequency-domain companion of `TsDriver`, using its original
//! mechanical assembly and the shared `fs_phs::baffled_piston_impedance` owner.
//! The piston impedance supplies BOTH radiation loss and reactive fluid loading
//! at the actual queried frequency; no fixed resistance, frequency pin or
//! additional compact air mass enters this path. The original pHS harmonic
//! solver owns factorization, equilibration, residuals and cancellation.
//!
//! At one frequency, `Z_rad = R_rad - i omega M_rad` is inserted into the
//! original mechanical equation. This pointwise equivalent mass is NOT a
//! causal time-domain realization or a stored-acoustic-energy model. Only
//! cycle-averaged powers, motion and phasors are published. Peak amplitudes use
//! exp(-i omega t), consistently with the acoustic owners. The load is one
//! uniform rigid circular face in an infinite baffle, not an unbaffled cabinet,
//! cone breakup, a panel radiation law or mutual exterior loading.
//!
//! Analytical reference: COMSOL, Baffled Piston Radiation, equations (2), (4):
//! <https://doc.comsol.com/6.3/doc/com.comsol.help.models.aco.baffled_piston_radiation/baffled_piston_radiation.html>.

use super::{DriverAssembly, PanelMode, SpeakerError, TsCard, assemble_driver_with_radiation};
use fs_material::gas::GasState;
use fs_math::{c64::C64, det};
use fs_phs::harmonic::{HarmonicBudget, HarmonicError, HarmonicModel};
use fs_phs::{PortHamiltonian, QuadraticStorage};

/// Source, harmonic solve, excursion or observer-domain refusal.
#[derive(Debug)]
pub enum BaffledDriverError {
    Source(SpeakerError),
    Harmonic(HarmonicError),
    Invalid(&'static str),
}

impl core::fmt::Display for BaffledDriverError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "baffled driver response refused: {self:?}")
    }
}
impl std::error::Error for BaffledDriverError {}
impl From<SpeakerError> for BaffledDriverError {
    fn from(error: SpeakerError) -> Self {
        Self::Source(error)
    }
}
impl From<fs_phs::PhsError> for BaffledDriverError {
    fn from(error: fs_phs::PhsError) -> Self {
        Self::Source(SpeakerError::Phs(error))
    }
}
impl From<HarmonicError> for BaffledDriverError {
    fn from(error: HarmonicError) -> Self {
        Self::Harmonic(error)
    }
}

/// Fixed original driver and media, reusable for independently chosen frequencies.
pub struct BaffledTsDriver {
    card: TsCard,
    exterior_gas: GasState,
    cavity_gas: Option<GasState>,
    panel: Option<PanelMode>,
    assembly: DriverAssembly,
    radius_m: f64,
}

/// Actual peak phasors and independent physical cycle-average power accounts.
#[derive(Clone, Debug, PartialEq)]
pub struct BaffledDriverResponse {
    pub angular_frequency_rad_s: f64,
    pub terminal_voltage_v: C64,
    pub coil_current_a: C64,
    pub cone_displacement_m: C64,
    pub cone_velocity_m_s: C64,
    /// The original retained panel coordinate; no panel radiation is inferred.
    pub panel_displacement_m: Option<C64>,
    pub panel_velocity_m_s: Option<C64>,
    /// Complete mechanical radiation impedance F/v [N s/m], exp(-i omega t).
    pub radiation_impedance_n_s_m: C64,
    pub supplied_power_w: f64,
    pub coil_heat_w: f64,
    pub suspension_heat_w: f64,
    pub panel_heat_w: f64,
    /// Exported acoustic power, not heat deposited in the driver or cavity.
    pub radiated_power_w: f64,
    /// Sum of independently computed sinks minus terminal electrical supply.
    pub power_residual_w: f64,
    /// Original equation residual from the shared harmonic solver.
    pub relative_residual: f64,
    radius_m: f64,
    exterior_gas: GasState,
    area_m2: f64,
}

/// Far-field cone pressure and geometric approximation limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BaffledFarFieldPressure {
    /// Peak exp(-i omega t) pressure, with the physical e^(ikr) propagation.
    pub pressure_pa: C64,
    /// Upper bound on the phase dropped by linearizing source-to-receiver range.
    pub maximum_path_phase_error_rad: f64,
    /// Bound a/(r-a) on replacing the Rayleigh kernel's 1/R by 1/r.
    /// This is a KERNEL bound, not a relative pressure bound near directivity zeros.
    pub maximum_kernel_amplitude_error: f64,
}

impl BaffledTsDriver {
    /// Bind the same coil, suspension, optional cavity and panel as `TsDriver`.
    /// Exterior and cavity gas admission is shared with the original constructor.
    /// `None` box means no rear-cavity stiffness; the front face remains baffled.
    pub fn new_with_media(
        card: TsCard,
        box_volume_m3: Option<f64>,
        panel: Option<PanelMode>,
        exterior_gas: GasState,
        cavity_gas: Option<GasState>,
    ) -> Result<Self, BaffledDriverError> {
        let assembly = assemble_driver_with_radiation(
            card,
            box_volume_m3,
            panel,
            &exterior_gas,
            cavity_gas.as_ref(),
            false,
        )?;
        let radius_m = det::sqrt(card.sd_m2 / core::f64::consts::PI);
        // The same storage/structure admission used by the time-domain owner.
        make_system(&assembly, card.mms_kg, card.rms_n_s_m)?;
        Ok(Self {
            card,
            exterior_gas,
            cavity_gas,
            panel,
            assembly,
            radius_m,
        })
    }

    pub fn exterior_gas(&self) -> &GasState {
        &self.exterior_gas
    }
    pub fn cavity_gas(&self) -> Option<&GasState> {
        self.cavity_gas.as_ref()
    }

    /// Solve at a positive continuous angular frequency for one peak terminal
    /// voltage. No sample clock, Nyquist cap, frequency interpolation or stepped
    /// sinusoid is involved. The 3- or 5-state system obeys the caller's original
    /// harmonic resources and residual tolerance. The linear travel limit still
    /// applies to the complete sinusoidal excursion, not only sampled endpoints.
    pub fn response(
        &self,
        omega_rad_s: f64,
        voltage_v: C64,
        budget: HarmonicBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BaffledDriverResponse, BaffledDriverError> {
        if cancelled() {
            return Err(HarmonicError::Cancelled.into());
        }
        if !omega_rad_s.is_finite() || omega_rad_s <= 0.0 || !finite(voltage_v) {
            return Err(BaffledDriverError::Invalid(
                "positive finite angular frequency and finite voltage required",
            ));
        }
        let (r, x) = fs_phs::baffled_piston_impedance(
            self.exterior_gas.density,
            self.exterior_gas.sound_speed,
            self.radius_m,
            omega_rad_s,
        )?;
        // p/U -> F/v uses the conjugate area at both sides of the port.
        let impedance = C64::new(r, x).scale(self.card.sd_m2).scale(self.card.sd_m2);
        let loaded_mass = self.card.mms_kg - impedance.im / omega_rad_s;
        let resistance = self.card.rms_n_s_m + impedance.re;
        if !finite(impedance)
            || !loaded_mass.is_finite()
            || loaded_mass <= 0.0
            || !resistance.is_finite()
            || resistance <= 0.0
        {
            return Err(BaffledDriverError::Invalid(
                "radiation load is not representable",
            ));
        }
        let system = make_system(&self.assembly, loaded_mass, resistance)?;
        let model = HarmonicModel::from_ode(&system, budget, &mut cancelled)?;
        let physical = model
            .factor(omega_rad_s, &mut cancelled)?
            .response(&[voltage_v], &mut cancelled)?;
        let displacement = physical.state[1];
        let excursion = displacement.abs();
        if excursion > self.card.x_max_m {
            return Err(SpeakerError::ExcursionExceeded {
                excursion_m: excursion,
                x_max_m: self.card.x_max_m,
            }
            .into());
        }
        let velocity = physical.effort[2];
        let current = physical.output[0];
        let panel_displacement_m = self.panel.map(|_| physical.state[3]);
        let panel_velocity_m_s = self.panel.map(|_| physical.effort[4]);
        let coil_heat_w = 0.5 * self.card.re_ohm * current.norm_sq();
        let suspension_heat_w = 0.5 * self.card.rms_n_s_m * velocity.norm_sq();
        let radiated_power_w = 0.5 * impedance.re * velocity.norm_sq();
        let panel_heat_w = self.panel.map_or(0.0, |p| {
            p.damping_ratio
                * p.mass_kg
                * core::f64::consts::TAU
                * p.frequency_hz
                * physical.effort[4].norm_sq()
        });
        let power_residual_w = coil_heat_w + suspension_heat_w + panel_heat_w + radiated_power_w
            - physical.supplied_power_w;
        if ![
            coil_heat_w,
            suspension_heat_w,
            panel_heat_w,
            radiated_power_w,
            power_residual_w,
        ]
        .iter()
        .all(|v| v.is_finite())
        {
            return Err(BaffledDriverError::Invalid(
                "physical cycle power is not representable",
            ));
        }
        if cancelled() {
            return Err(HarmonicError::Cancelled.into());
        }
        Ok(BaffledDriverResponse {
            angular_frequency_rad_s: omega_rad_s,
            terminal_voltage_v: voltage_v,
            coil_current_a: current,
            cone_displacement_m: displacement,
            cone_velocity_m_s: velocity,
            panel_displacement_m,
            panel_velocity_m_s,
            radiation_impedance_n_s_m: impedance,
            supplied_power_w: physical.supplied_power_w,
            coil_heat_w,
            suspension_heat_w,
            panel_heat_w,
            radiated_power_w,
            power_residual_w,
            relative_residual: physical.relative_residual,
            radius_m: self.radius_m,
            exterior_gas: self.exterior_gas,
            area_m2: self.card.sd_m2,
        })
    }
}

fn make_system(
    a: &DriverAssembly,
    mass: f64,
    resistance: f64,
) -> Result<PortHamiltonian, fs_phs::PhsError> {
    let mut q = a.q.clone();
    let mut r = a.r.clone();
    q[2 * a.n + 2] = 1.0 / mass;
    r[2 * a.n + 2] = resistance;
    PortHamiltonian::new(
        a.n,
        1,
        a.j.clone(),
        r,
        a.g.clone(),
        Box::new(QuadraticStorage::new(q, a.n)?),
    )
}

impl BaffledDriverResponse {
    /// Cone-only Rayleigh far field in the front hemisphere. The Bessel
    /// directivity keeps sign reversals at sidelobes, and the 2*pi half-space
    /// normalization agrees with the radiation load's exported power.
    ///
    /// The caller bounds both geometric approximations explicitly. For r>a,
    /// Taylor's theorem bounds the omitted path by a²/[2(r-a)]; multiplying by
    /// k bounds phase, while a/(r-a) bounds the inverse-distance kernel. These
    /// binary64 model estimates are not interval certificates or relative
    /// pressure guarantees at destructive-interference zeros. No panel output,
    /// atmospheric absorption, room, near-field or moving-source term is added.
    pub fn far_field_pressure(
        &self,
        range_m: f64,
        polar_angle_rad: f64,
        maximum_phase_error_rad: f64,
        maximum_kernel_amplitude_error: f64,
    ) -> Result<BaffledFarFieldPressure, BaffledDriverError> {
        if !range_m.is_finite()
            || range_m <= self.radius_m
            || !polar_angle_rad.is_finite()
            || !(0.0..=core::f64::consts::FRAC_PI_2).contains(&polar_angle_rad)
            || !maximum_phase_error_rad.is_finite()
            || maximum_phase_error_rad <= 0.0
            || !maximum_kernel_amplitude_error.is_finite()
            || maximum_kernel_amplitude_error <= 0.0
        {
            return Err(BaffledDriverError::Invalid(
                "receiver needs finite front-hemisphere geometry and positive error limits",
            ));
        }
        let k = self.angular_frequency_rad_s / self.exterior_gas.sound_speed;
        let amplitude_error = self.radius_m / (range_m - self.radius_m);
        let phase_error = 0.5 * (k * self.radius_m) * amplitude_error;
        let phase = k * range_m;
        if !phase_error.is_finite()
            || !phase.is_finite()
            || phase_error > maximum_phase_error_rad
            || amplitude_error > maximum_kernel_amplitude_error
        {
            return Err(BaffledDriverError::Invalid(
                "receiver lies outside the declared far-field accuracy",
            ));
        }
        let z = k * self.radius_m * det::sin(polar_angle_rad);
        let directivity = if z == 0.0 {
            1.0
        } else {
            2.0 * fs_math::special::bessel_j1(z) / z
        };
        let (sin, cos) = det::sin_cos(phase);
        let pressure = C64::new(sin, -cos)
            * self.cone_velocity_m_s.scale(
                self.exterior_gas.density * self.angular_frequency_rad_s * self.area_m2
                    / (core::f64::consts::TAU * range_m)
                    * directivity,
            );
        if !finite(pressure) {
            return Err(BaffledDriverError::Invalid(
                "receiver pressure is not representable",
            ));
        }
        Ok(BaffledFarFieldPressure {
            pressure_pa: pressure,
            maximum_path_phase_error_rad: phase_error,
            maximum_kernel_amplitude_error: amplitude_error,
        })
    }
}

fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}

#[cfg(test)]
mod tests;
