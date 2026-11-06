//! Fully developed Newtonian film with prescribed tangential surface traction.
//!
//! With a stationary no-slip wall at y=0 and free surface at y=h,
//! `mu*u'' = -rho*g_t`, `u(0)=0`, and `mu*u'(h)=tau_s`. Hence
//! `u(y) = g_t*(h*y-y*y/2)/nu + tau_s*y/(rho*nu)` and
//! `q = g_t*h^3/(3*nu) + tau_s*h^2/(2*rho*nu)`.
//! The surface traction is an explicit input, for example from gas drag or a
//! local surface-tension gradient. This owner does not infer it from gas speed.
//! The additive shear profile and flux agree with Hossain et al.,
//! <https://arxiv.org/abs/2405.12623>, equations (25)-(26), in the flat,
//! constant-property limit. Power follows directly by integrating `mu*u'^2`.
//!
//! Constant properties, thickness, acceleration and traction are held. No wave
//! stability, transient inertia, gas evolution or thermal evolution is implied.

use super::{FilmError, FilmMaterial, represented};

/// Mechanical power per wetted wall area for a stationary no-slip wall.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FilmPower {
    /// Signed body-force work `rho*g_t*q` [W/m²].
    pub body_work_w_m2: f64,
    /// Signed work of the imposed surface traction `tau_s*u(h)` [W/m²].
    pub surface_work_w_m2: f64,
    /// Nonnegative `integral_0^h mu*(du/dy)^2 dy` [W/m²]. It equals
    /// body plus surface work in the admitted steady model, up to roundoff.
    /// A thermal consumer may use this once as viscous heat; this owner does
    /// not itself change temperature or remove energy from another fluid.
    pub viscous_dissipation_w_m2: f64,
}

/// A signed Couette--Nusselt profile at a prescribed positive film thickness.
///
/// Positive velocity, acceleration and surface traction share one tangential
/// direction. Opposing shear can reverse part or all of the profile; zero net
/// flux may still contain flow and positive dissipation. Construction and all
/// queries allocate nothing. Nonzero unrepresentable coefficients refuse.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ShearedFilm {
    thickness: f64,
    body_velocity: f64,
    shear_velocity: f64,
    body_shear: f64,
    surface_shear: f64,
    mean_velocity: f64,
    flow: f64,
    power: FilmPower,
}

impl ShearedFilm {
    /// Admit a flat laminar film on a stationary no-slip wall. Acceleration
    /// and surface traction may have either sign or be zero. A positive
    /// thickness is required even when both driving terms vanish.
    pub fn new(
        material: FilmMaterial,
        thickness_m: f64,
        tangential_acceleration_m_s2: f64,
        surface_shear_pa: f64,
    ) -> Result<Self, FilmError> {
        let h = thickness_m;
        let rho = material.density_kg_m3;
        let nu = material.kinematic_viscosity_m2_s;
        let g = tangential_acceleration_m_s2;
        for (name, value) in [
            ("thickness", h),
            ("density", rho),
            ("kinematic viscosity", nu),
        ] {
            if !value.is_finite() || value <= 0.0 {
                return Err(FilmError::InvalidInput(name));
            }
        }
        if !material.surface_tension_n_m.is_finite() || material.surface_tension_n_m < 0.0 {
            return Err(FilmError::InvalidInput("surface tension"));
        }
        for (name, value) in [
            ("tangential acceleration", g),
            ("surface shear", surface_shear_pa),
        ] {
            if !value.is_finite() {
                return Err(FilmError::InvalidInput(name));
            }
        }
        let body_velocity = product_ratio([g, h, h], [nu])?;
        let shear_velocity = product_ratio([surface_shear_pa, h], [rho, nu])?;
        let body_shear = product_ratio([rho, g, h], [])?;
        let wall_shear = body_shear + surface_shear_pa;
        let surface_velocity = 0.5 * body_velocity + shear_velocity;
        let mean_velocity = body_velocity / 3.0 + shear_velocity / 2.0;
        for value in [wall_shear, surface_velocity, mean_velocity] {
            represented(value, false)?;
        }
        let flow = product_ratio([h, mean_velocity], [])?;

        // Integrate the linear shear stress in a scaled positive form. This
        // avoids cancellation of body and surface work during counterflow and
        // overflowing stress squared before viscosity or thickness scales it.
        let scale = wall_shear.abs().max(surface_shear_pa.abs());
        let dissipation = if scale == 0.0 {
            0.0
        } else {
            let a = wall_shear / scale;
            let b = surface_shear_pa / scale;
            let norm = (a * a + a * b + b * b) / 3.0;
            product_ratio([scale, scale, h, norm], [rho, nu])?
        };
        Ok(Self {
            thickness: h,
            body_velocity,
            shear_velocity,
            body_shear,
            surface_shear: surface_shear_pa,
            mean_velocity,
            flow,
            power: FilmPower {
                body_work_w_m2: product_ratio([rho, g, flow], [])?,
                surface_work_w_m2: product_ratio([surface_shear_pa, surface_velocity], [])?,
                viscous_dissipation_w_m2: dissipation,
            },
        })
    }

    #[must_use]
    pub const fn thickness_m(self) -> f64 {
        self.thickness
    }

    /// Signed volumetric flux per width [m²/s], including any counterflow.
    #[must_use]
    pub const fn flow_per_width_m2_s(self) -> f64 {
        self.flow
    }

    #[must_use]
    pub const fn mean_velocity_m_s(self) -> f64 {
        self.mean_velocity
    }

    #[must_use]
    pub fn surface_velocity_m_s(self) -> f64 {
        0.5 * self.body_velocity + self.shear_velocity
    }

    /// Signed fluid-on-wall tangential traction [Pa]. Wall-on-fluid traction
    /// has the opposite sign, since the outward normal there points into the wall.
    #[must_use]
    pub fn wall_shear_pa(self) -> f64 {
        self.body_shear + self.surface_shear
    }

    /// Prescribed traction on the liquid at its outward-facing free surface.
    #[must_use]
    pub const fn surface_shear_pa(self) -> f64 {
        self.surface_shear
    }

    #[must_use]
    pub const fn power(self) -> FilmPower {
        self.power
    }

    /// Signed tangential velocity [m/s] at wall distance y in [0,h].
    pub fn velocity_m_s(self, y_m: f64) -> Result<f64, FilmError> {
        let z = self.coordinate(y_m)?;
        let coefficient = self
            .body_velocity
            .mul_add(1.0 - 0.5 * z, self.shear_velocity);
        let value = (z * self.body_velocity).mul_add(1.0 - 0.5 * z, z * self.shear_velocity);
        represented(value, z != 0.0 && coefficient != 0.0)?;
        Ok(value)
    }

    /// Signed Cauchy shear component `mu*du/dy` [Pa] at y in [0,h].
    /// It is not a magnitude: an internal zero indicates a velocity extremum.
    pub fn shear_pa(self, y_m: f64) -> Result<f64, FilmError> {
        self.coordinate(y_m)?;
        let value = self
            .body_shear
            .mul_add((self.thickness - y_m) / self.thickness, self.surface_shear);
        represented(value, false)?;
        Ok(value)
    }

    fn coordinate(self, y: f64) -> Result<f64, FilmError> {
        if !y.is_finite() || y < 0.0 || y > self.thickness {
            return Err(FilmError::InvalidInput("wall-normal coordinate"));
        }
        let z = y / self.thickness;
        represented(z, y != 0.0)?;
        Ok(z)
    }
}

// Only construction uses this fallback. Ordinary coefficients retain direct
// IEEE arithmetic; deterministic logarithmic scaling rescues a representable
// result when intermediate multiplication/division leaves the normal range.
fn product_ratio<const N: usize, const D: usize>(
    numerator: [f64; N],
    denominator: [f64; D],
) -> Result<f64, FilmError> {
    if numerator.contains(&0.0) {
        return Ok(0.0);
    }
    let mut value = 1.0;
    let mut ordinary = true;
    for factor in numerator {
        value *= factor;
        ordinary &= value.is_normal();
    }
    for factor in denominator {
        value /= factor;
        ordinary &= value.is_normal();
    }
    if !ordinary {
        let mut sign = 1.0;
        let mut logarithm = 0.0;
        for factor in numerator {
            sign *= factor.signum();
            logarithm += fs_math::det::ln(factor.abs());
        }
        for factor in denominator {
            logarithm -= fs_math::det::ln(factor);
        }
        value = sign * fs_math::det::exp(logarithm);
    }
    represented(value, true)?;
    Ok(value)
}

#[cfg(test)]
mod tests;
