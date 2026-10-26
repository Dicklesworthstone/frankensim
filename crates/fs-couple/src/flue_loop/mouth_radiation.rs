//! Compact self reaction of the actual rectangular flue mouth in its exterior.
//!
//! The same uniform area cells as the pressure receiver enter the shared
//! Rayleigh work integral. The exact rectangular static inertance M and monopole
//! resistance coefficient C = rho/(2*pi*c) determine R = M*M/C and p = M/C in
//! the positive-real load R*s/(s+p). This is a two-term asymptotic realization,
//! not an equivalent-circle end correction or an authored damping constant.
//! The shared fs-vfit midpoint owner retains near-field energy and escaped
//! radiation separately. Both continuous and realized responses must agree with
//! the rectangle's Rayleigh image within 5% on the declared compact band.
//! The static integral independently checks the surface quadrature. Finite-band
//! sampled checks still retain quadrature error; they are not a continuum
//! certificate, a broadband load, or mouth-to-bell mutual radiation.

use super::{FlueError, FlueVoice};
use fs_bem::helmholtz::{
    Medium,
    rayleigh_modal::{AreaSample, PreparedSurface},
};
use fs_exec::CancelGate;
use fs_material::gas::GasState;
use fs_vfit::impedance::SeriesImpedanceSpec;
use fs_vfit::relaxation::{
    RelaxationFrame, RelaxationImpedance, RelaxationImpedanceSpec, RelaxationTerm,
};

/// Bounded cold construction: a complete 32-frequency rectangular self integral.
pub const MAX_MOUTH_REACTION_CELLS: usize = 1024;

fn invalid(what: &'static str) -> FlueError {
    FlueError::Invalid { what }
}
fn physics(error: impl core::fmt::Display) -> FlueError {
    FlueError::Physics(crate::acoustic_realize::AcousticRealizeError::Nonlinear(
        error.to_string(),
    ))
}

/// Immutable geometry, complete exterior gas and actual passive coefficients.
/// Constructed by the performance owner from its already admitted mouth and
/// microphone. No independently authored radiation radius or loss enters it.
#[derive(Clone, Copy, Debug)]
pub struct RectangularMouthLoad {
    size_m: [f64; 2],
    cells: [usize; 2],
    exterior: GasState,
    time_step_s: f64,
    maximum_frequency_hz: f64,
    static_quadrature_relative_error: f64,
    max_complex_relative_error: f64,
    max_resistance_relative_error: f64,
    load: RelaxationImpedanceSpec,
}

impl RectangularMouthLoad {
    /// Derive one compact load from the same area basis used by the receiver.
    /// `k*diagonal <= 0.5` and `dt*fmax <= 0.1` bound this approximation;
    /// independent 5% complex AND resistance checks can refuse sooner.
    /// Cell self reactance uses the shared Rayleigh owner's equivalent-cell
    /// integral; the entire aperture retains its rectangular shape.
    ///
    /// # Errors
    /// Invalid gas/geometry/clock, unresolved static quadrature, excessive
    /// compact-band error, cancellation or any original passive-owner refusal.
    pub fn new(
        size_m: [f64; 2],
        cells: [usize; 2],
        exterior: GasState,
        time_step_s: f64,
        maximum_frequency_hz: f64,
        gate: &CancelGate,
    ) -> Result<Self, FlueError> {
        exterior.validate_coherence().map_err(physics)?;
        let count = cells[0].checked_mul(cells[1]);
        if size_m.iter().any(|x| !x.is_finite() || *x <= 0.0)
            || cells.iter().any(|n| !(1..=512).contains(n))
            || count.is_none_or(|n| n > MAX_MOUTH_REACTION_CELLS)
            || !time_step_s.is_finite()
            || time_step_s <= 0.0
            || !maximum_frequency_hz.is_finite()
            || maximum_frequency_hz <= 0.0
        {
            return Err(invalid(
                "mouth reaction requires physical geometry, gas, clock and at most 1024 cells",
            ));
        }
        let omega_max = core::f64::consts::TAU * maximum_frequency_hz;
        let phase = omega_max * size_m[0].hypot(size_m[1]) / exterior.sound_speed;
        let area = size_m[0] * size_m[1];
        if !phase.is_finite()
            || phase > 0.5
            || time_step_s * maximum_frequency_hz > 0.1
            || !area.is_finite()
            || area <= 0.0
            || !(1.0 / area).is_finite()
        {
            return Err(invalid(
                "compact rectangular reaction requires k*diagonal <= 0.5 and dt*fmax <= 0.1",
            ));
        }
        if gate.is_requested() {
            return Err(invalid("mouth radiation construction cancelled"));
        }
        let count = count.expect("admitted cell product");
        let mut samples = Vec::with_capacity(count);
        for y in 0..cells[1] {
            for x in 0..cells[0] {
                samples.push(AreaSample {
                    position_m: [
                        (x as f64 + 0.5) * (size_m[0] / cells[0] as f64) - 0.5 * size_m[0],
                        (y as f64 + 0.5) * (size_m[1] / cells[1] as f64) - 0.5 * size_m[1],
                    ],
                    area_m2: area / count as f64,
                    mode_shape: vec![1.0 / area],
                });
            }
        }
        let surface = PreparedSurface::new(&samples).map_err(physics)?;
        let medium = Medium {
            density: exterior.density,
            sound_speed: exterior.sound_speed,
        };
        // The actual rectangle's exact static integral sets its near-field
        // inertia. Coarse equivalent-cell quadrature cannot silently replace
        // that shape with an inaccurate acoustic end correction.
        let inertance = rectangular_inertance(size_m, exterior.density)?;
        // This small phase gives O(1e-14) truncation error in the static limit.
        // Keep the original Rayleigh kernel, including its finite-cell self term.
        let omega_small = omega_max / 1_048_576.0;
        let quadrature_inertance = -surface
            .power_consistent_impedance(omega_small, medium)
            .map_err(physics)?[0]
            .im
            / omega_small;
        let static_quadrature_relative_error = (quadrature_inertance / inertance - 1.0).abs();
        if !static_quadrature_relative_error.is_finite() || static_quadrature_relative_error > 0.05
        {
            return Err(invalid(
                "mouth surface exceeds 5% exact rectangular static-inertia error; refine its cells",
            ));
        }
        let resistance_coefficient =
            exterior.density / (core::f64::consts::TAU * exterior.sound_speed);
        let rate = inertance / resistance_coefficient;
        let term = RelaxationTerm {
            resistance_pa_s_m3: inertance * rate,
            rate_per_s: rate,
        };
        let load = RelaxationImpedanceSpec::new(
            SeriesImpedanceSpec {
                resistance_pa_s_m3: 0.0,
                inertance_pa_s2_m3: 0.0,
                compliance_m3_pa: None,
            },
            &[term],
        )
        .map_err(physics)?;
        RelaxationImpedance::new(
            load,
            exterior.density * exterior.sound_speed / area,
            time_step_s,
        )
        .map_err(physics)?;
        let mut result = Self {
            size_m,
            cells,
            exterior,
            time_step_s,
            maximum_frequency_hz,
            static_quadrature_relative_error,
            max_complex_relative_error: 0.0,
            max_resistance_relative_error: 0.0,
            load,
        };
        for i in 1..=32 {
            if gate.is_requested() {
                return Err(invalid("mouth radiation construction cancelled"));
            }
            let frequency = maximum_frequency_hz * f64::from(i) / 32.0;
            let exact = surface
                .power_consistent_impedance(core::f64::consts::TAU * frequency, medium)
                .map_err(physics)?[0];
            for discrete in [false, true] {
                let (r, x) = result.impedance_at(frequency, discrete)?;
                let complex = (r - exact.re).hypot(x - exact.im) / exact.re.hypot(exact.im);
                let resistance = (r - exact.re).abs() / exact.re;
                if !complex.is_finite()
                    || !resistance.is_finite()
                    || complex > 0.05
                    || resistance > 0.05
                {
                    return Err(invalid(
                        "compact rectangular reaction exceeds the fixed 5% complex or resistance check",
                    ));
                }
                result.max_complex_relative_error = result.max_complex_relative_error.max(complex);
                result.max_resistance_relative_error =
                    result.max_resistance_relative_error.max(resistance);
            }
        }
        Ok(result)
    }

    /// Analogue or actually realized bilinear impedance under exp(-i*omega*t).
    ///
    /// # Errors
    /// Nonfinite, negative or out-of-band frequency, or unrepresentable response.
    pub fn impedance_at(&self, frequency_hz: f64, discrete: bool) -> Result<(f64, f64), FlueError> {
        if !frequency_hz.is_finite()
            || frequency_hz < 0.0
            || frequency_hz > self.maximum_frequency_hz
        {
            return Err(invalid(
                "mouth radiation query exceeds its admitted frequency band",
            ));
        }
        let mut omega = core::f64::consts::TAU * frequency_hz;
        if discrete {
            let angle = omega * (0.5 * self.time_step_s);
            omega = 2.0 * (fs_math::det::sin(angle) / fs_math::det::cos(angle)) / self.time_step_s;
        }
        let term = self.load.terms()[0];
        let ratio = omega / term.rate_per_s;
        let x = term.resistance_pa_s_m3 * ratio / (1.0 + ratio * ratio);
        let r = x * ratio;
        if !r.is_finite() || !x.is_finite() {
            return Err(invalid("mouth radiation impedance overflow"));
        }
        Ok((r, -x))
    }
    /// Complete passive coefficients, with no legacy end correction in series.
    #[must_use]
    pub const fn load(&self) -> RelaxationImpedanceSpec {
        self.load
    }
    /// Physical width and cut-up, not an equivalent radius.
    #[must_use]
    pub const fn size_m(&self) -> [f64; 2] {
        self.size_m
    }
    /// Shared source/receiver area quadrature.
    #[must_use]
    pub const fn cells(&self) -> [usize; 2] {
        self.cells
    }
    /// Complete gas used by both reaction and outgoing propagation.
    #[must_use]
    pub const fn exterior(&self) -> &GasState {
        &self.exterior
    }
    /// Declared compact-band upper frequency.
    #[must_use]
    pub const fn maximum_frequency_hz(&self) -> f64 {
        self.maximum_frequency_hz
    }
    /// Static surface quadrature discrepancy from the exact rectangular integral.
    #[must_use]
    pub const fn static_quadrature_relative_error(&self) -> f64 {
        self.static_quadrature_relative_error
    }
    /// Largest sampled continuous/discrete discrepancy from the shared surface.
    #[must_use]
    pub const fn max_complex_relative_error(&self) -> f64 {
        self.max_complex_relative_error
    }
    /// Independently checked weak radiation resistance discrepancy.
    #[must_use]
    pub const fn max_resistance_relative_error(&self) -> f64 {
        self.max_resistance_relative_error
    }
}

// Integrate rho/(2*pi*A^2) * 4*(a-x)*(b-y)/sqrt(x*x+y*y) over
// 0<x<a, 0<y<b. Dividing the closed form by t^2 first avoids cancellation
// in 1+t^3-(1+t^2)^(3/2) for a slender rectangle. This is an exact
// geometric integral, not a fit or an equivalent circular aperture.
fn rectangular_inertance(size: [f64; 2], density: f64) -> Result<f64, FlueError> {
    fs_bem::helmholtz::rayleigh_modal::rectangular_static_inertance(size, density)
        .map_err(physics)
}

/// A wave-port adapter only: all memory, stepping and energy remain in fs-vfit.
#[derive(Clone, Copy, Debug)]
pub(super) struct MouthRadiation {
    pub model: RectangularMouthLoad,
    state: RelaxationImpedance,
    port_impedance: f64,
    direct: f64,
}
impl MouthRadiation {
    pub fn new(model: RectangularMouthLoad, port_impedance: f64) -> Result<Self, FlueError> {
        let state = RelaxationImpedance::new(model.load, port_impedance, model.time_step_s)
            .map_err(physics)?;
        let term = model.load.terms()[0];
        let retained = 1.0 / (1.0 + term.rate_per_s * (0.5 * model.time_step_s));
        let direct = term.resistance_pa_s_m3 * retained;
        Ok(Self {
            model,
            state,
            port_impedance,
            direct,
        })
    }
    /// Midpoint pressure = direct*CURRENT mouth flow + retained history.
    /// This is the shared owner's algebraic elimination, not a lagged load.
    pub fn affine_pressure(&self) -> (f64, f64) {
        (
            self.direct,
            -self.direct * self.state.branch_flows_m3_s()[0],
        )
    }
    pub fn candidate(&self, flow: f64) -> Result<(Self, RelaxationFrame), FlueError> {
        let (direct, history) = self.affine_pressure();
        let pressure = direct * flow + history;
        let incident = 0.5 * (pressure + self.port_impedance * flow);
        let mut next = *self;
        let frame = next.state.step(incident).map_err(physics)?;
        let scale = flow
            .abs()
            .max(self.state.branch_flows_m3_s()[0].abs())
            .max(f64::MIN_POSITIVE);
        if (frame.port.flow_m3_s - flow).abs() > 64.0 * f64::EPSILON * scale {
            return Err(invalid(
                "mouth radiation wave adapter disagrees with the coupled port flow",
            ));
        }
        Ok((next, frame))
    }
}

impl FlueVoice {
    pub(super) fn attach_mouth_radiation(
        &mut self,
        model: RectangularMouthLoad,
        aperture_size_m: [f64; 2],
    ) -> Result<(), FlueError> {
        if self.poisoned || self.samples != 0 || self.mouth_radiation.is_some()
            || self.window_radiation.is_some()
            || self.dg().is_some_and(|duct| duct.load().radiation_record()
                .is_some_and(|record| record.driven_mouth().is_some()))
        {
            return Err(invalid("mouth radiation requires one cold attachment"));
        }
        if self.geometry.mouth_loss_rel != 0.0 {
            return Err(invalid(
                "physical mouth radiation replaces the legacy loss; declare mouth_loss_rel zero",
            ));
        }
        // The acoustic aperture is bound by the performer that owns its
        // receiver. A lip jet can be narrower than this opening.
        if model.size_m != aperture_size_m
            || model.time_step_s.to_bits() != self.dt.to_bits()
            || model.exterior.pressure.to_bits() != self.gas.pressure.to_bits()
        {
            return Err(invalid(
                "mouth radiation must share declared aperture, clock and static pressure",
            ));
        }
        let load = MouthRadiation::new(model, self.zc)?;
        self.mouth_radiation = Some(load);
        // Replace both old passive coefficients. Aerodynamic source law and
        // independently declared separation closure retain their own owners.
        self.inertance = 0.0;
        self.resistance = 0.0;
        Ok(())
    }
    /// Actual compact rectangular reaction, absent on the unchanged legacy path.
    #[must_use]
    pub fn mouth_radiation_load(&self) -> Option<&RectangularMouthLoad> {
        self.mouth_radiation.as_ref().map(|load| &load.model)
    }
}

#[cfg(test)]
mod tests;
