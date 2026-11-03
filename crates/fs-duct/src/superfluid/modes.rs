//! Energy-normalized singular modes of the linear Landau equations.

use super::{CellSample, ChannelError, FaceSample};
use fs_material::superfluid::TwoFluidState;

#[derive(Clone, Copy, Debug)]
pub(super) struct Modes {
    pub speed: [f64; 2],
    // Rows are orthonormal left/right singular vectors of G.
    thermodynamic: [[f64; 2]; 2],
    kinetic: [[f64; 2]; 2],
    pressure_scale: f64,
    temperature_scale: f64,
    pressure_temperature: f64,
    bulk_scale: f64,
    relative_scale: f64,
}

impl Modes {
    pub fn new(state: &TwoFluidState) -> Result<Self, ChannelError> {
        let rho = state.density_kg_m3();
        let t = state.temperature_k();
        let cp = state.specific_heat_cp_j_kg_k();
        let fs = state.superfluid_fraction();
        let fn_ = state.normal_fraction();
        let c = state.adiabatic_sound_speed_m_s();
        let u = state.entropy_j_kg_k() * (fs / fn_ * (t / cp)).sqrt();
        let g = c * state.thermal_expansion_per_k() * (t / cp).sqrt() * u;
        // Whiten thermodynamics by its exact Cholesky factor, and momentum
        // by the component kinetic masses. Then G = [[c,g],[0,u]], and the
        // complete finite-expansion system has Euclidean quadratic storage.
        // Scale before the 2x2 eigensolve so sound speeds are never squared
        // merely to discover the singular-vector rotation.
        let scale = c.max(u).max(g.abs());
        let (cs, us, gs) = (c / scale, u / scale, g / scale);
        let difference = cs * cs + gs * gs - us * us;
        let off = gs * us;
        let (cos, sin) = if off == 0.0 {
            if difference >= 0.0 {
                (1.0, 0.0)
            } else {
                (0.0, 1.0)
            }
        } else {
            let gap = difference.hypot(2.0 * off);
            if difference >= 0.0 {
                let cos = (0.5 * (1.0 + difference / gap)).sqrt();
                (cos, off / (gap * cos))
            } else {
                let sin = (0.5 * (1.0 - difference / gap)).sqrt().copysign(off);
                (off / (gap * sin), sin)
            }
        };
        let norm = cos.hypot(sin);
        let left = [cos / norm, sin / norm];
        let right = [cs * left[0], gs * left[0] + us * left[1]];
        let norm = right[0].hypot(right[1]);
        let right = [right[0] / norm, right[1] / norm];
        let sound = state.sound_modes();
        let modes = Self {
            speed: [sound.fast_m_s, sound.slow_m_s],
            thermodynamic: [left, [-left[1], left[0]]],
            kinetic: [right, [-right[1], right[0]]],
            pressure_scale: state.adiabatic_compressibility_per_pa().sqrt(),
            temperature_scale: (rho * cp / t).sqrt(),
            pressure_temperature: state.thermal_expansion_per_k() * t / (rho * cp),
            bulk_scale: rho.sqrt(),
            relative_scale: (rho * fs * fn_).sqrt(),
        };
        if [
            modes.pressure_scale,
            modes.temperature_scale,
            modes.bulk_scale,
            modes.relative_scale,
            norm,
            scale,
        ]
        .into_iter()
        .any(|x| !x.is_finite() || x <= 0.0)
            || !modes.pressure_temperature.is_finite()
        {
            return Err(ChannelError::Invalid(
                "two-fluid mode transformation is not representable",
            ));
        }
        // Refuse a poorly representable transformation instead of silently
        // replacing a coupled thermal mode with two unrelated wave copies.
        for mode in 0..2 {
            let l = modes.thermodynamic[mode];
            let r = modes.kinetic[mode];
            let expected = modes.speed[mode] / scale;
            let residual = (cs * r[0] + gs * r[1] - expected * l[0])
                .abs()
                .max((us * r[1] - expected * l[1]).abs());
            if !residual.is_finite() || residual > 512.0 * f64::EPSILON {
                return Err(ChannelError::Invalid(
                    "two-fluid singular modes are unresolved",
                ));
            }
        }
        Ok(modes)
    }

    pub fn cell_to_modes(self, pressure: f64, temperature: f64) -> [f64; 2] {
        let p = self.pressure_scale * pressure;
        let t = self.temperature_scale * (temperature - self.pressure_temperature * pressure);
        self.thermodynamic.map(|row| row[0] * p + row[1] * t)
    }

    pub fn face_to_modes(self, bulk: f64, relative: f64) -> [f64; 2] {
        let v = self.bulk_scale * bulk;
        let w = self.relative_scale * relative;
        self.kinetic.map(|row| row[0] * v + row[1] * w)
    }

    pub fn cell(self, state: &TwoFluidState, values: [f64; 2]) -> CellSample {
        let pressure = (self.thermodynamic[0][0] * values[0]
            + self.thermodynamic[1][0] * values[1])
            / self.pressure_scale;
        let temperature = (self.thermodynamic[0][1] * values[0]
            + self.thermodynamic[1][1] * values[1])
            / self.temperature_scale
            + self.pressure_temperature * pressure;
        CellSample {
            pressure_perturbation_pa: pressure,
            temperature_perturbation_k: temperature,
            density_perturbation_kg_m3: state.density_kg_m3()
                * (state.isothermal_compressibility_per_pa() * pressure
                    - state.thermal_expansion_per_k() * temperature),
            specific_entropy_perturbation_j_kg_k: state.specific_heat_cp_j_kg_k()
                / state.temperature_k()
                * temperature
                - state.thermal_expansion_per_k() / state.density_kg_m3() * pressure,
        }
    }

    pub fn face(self, state: &TwoFluidState, values: [f64; 2]) -> FaceSample {
        let v = (self.kinetic[0][0] * values[0] + self.kinetic[1][0] * values[1]) / self.bulk_scale;
        let w =
            (self.kinetic[0][1] * values[0] + self.kinetic[1][1] * values[1]) / self.relative_scale;
        let vn = v + state.superfluid_fraction() * w;
        FaceSample {
            normal_velocity_m_s: vn,
            superfluid_velocity_m_s: v - state.normal_fraction() * w,
            barycentric_velocity_m_s: v,
            counterflow_velocity_m_s: w,
            mass_flux_kg_m2_s: state.density_kg_m3() * v,
            entropy_flux_w_m2_k: state.density_kg_m3() * state.entropy_j_kg_k() * vn,
            counterflow_heat_flux_w_m2: state.superfluid_density_kg_m3()
                * state.entropy_j_kg_k()
                * state.temperature_k()
                * w,
        }
    }
}
