//! Conductive stationary walls coupled to the original gas thermal layer.
//!
//! Under exp(+i omega t), q=(1+i)/delta and Y=k*q. Continuity of
//! temperature and normal heat flux gives T_wall/T_bulk=Y_g/(Y_g+Y_w).
//! The gas temperature deficit multiplies the original thermal Wentzell term.
//! For a homogeneous planar slab, Y_w=k_w*q_w*coth(q_w*d) with an isothermal
//! backing, or k_w*q_w*tanh(q_w*d) with an insulated backing. These follow
//! directly from k*T''=i*omega*C_vol*T. The infinite-depth limit is Y_w=k_w*q_w.
//! See the conductive-wall BLI formulation and finite-wall boundary conditions:
//! <https://doc.comsol.com/6.3/doc/com.comsol.help.aco/aco_ug_pressure.05.032.html>.
//!
//! This is local normal Fourier conduction, linearized at one common gas/solid
//! equilibrium temperature. Source selection and caloric constraints belong to
//! the material owner. No thermoelastic heating, interface resistance, changing
//! mean temperature, lateral conduction, multilayer wall or time-step state is
//! inferred. The wall's physical normal and tangential length scales must admit
//! this one-dimensional thermal profile, independently of mesh resolution.

use super::{
    AcousticFeError, BoundaryLayerGeometry, BoundaryLayerTransport, Medium, ReducedHelmholtz,
    invalid,
};
use fs_math::det;

/// Local source-evaluated solid transport, under the caller's caloric constraint.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WallThermalTransport {
    /// Conductivity in the wall-normal direction [W/(m K)].
    pub conductivity_w_m_k: f64,
    /// Volumetric heat capacity [J/(m3 K)].
    pub volumetric_heat_capacity_j_m3_k: f64,
}

/// Homogeneous wall depth and its remote thermal boundary condition.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ThermalWallDepth {
    /// A semi-infinite solid: the oscillating temperature decays at infinity.
    HalfSpace,
    /// A planar finite-thickness solid with the stated backside condition.
    Slab {
        /// Physical wall thickness [m], positive and finite.
        thickness_m: f64,
        /// Oscillatory temperature or heat flux prescribed at the backside.
        backing: ThermalBacking,
    },
}

/// Remote boundary of a finite wall; these act on the temperature perturbation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ThermalBacking {
    /// Zero oscillating temperature at the remote face.
    Isothermal,
    /// Zero normal heat flux at the remote face.
    Adiabatic,
}

/// A stationary, perfectly contacting solid with a local normal heat profile.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConductiveThermalWall {
    /// Original material transport at the common equilibrium temperature.
    pub transport: WallThermalTransport,
    /// Physical wall depth and remote thermal condition.
    pub depth: ThermalWallDepth,
    /// Smallest physical tangential variation length [m] of pressure, wall
    /// properties, thickness or remote condition; infinity only when constant.
    /// This is not the face diameter. The penetrated depth must be smaller by
    /// the original geometry's maximum_layer_ratio to neglect lateral diffusion.
    pub minimum_tangential_scale_m: f64,
}

/// Thermal boundary attached to the same no-slip, stationary acoustic wall.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ThermalWall {
    /// Original zero-temperature-perturbation wall, retained bitwise.
    Isothermal,
    /// No thermal exchange. The original viscous boundary layer remains active.
    Adiabatic,
    /// Conjugate gas/solid temperature and normal heat-flux continuity.
    Conductive(ConductiveThermalWall),
}

/// Dimensionless response to one unit of the bulk adiabatic gas temperature.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalWallResponse {
    /// Complex T_wall/T_bulk, using exp(+i omega t).
    pub wall_temperature: [f64; 2],
    /// Complex (T_bulk-T_wall)/T_bulk; multiplies the thermal Wentzell form.
    pub gas_temperature_deficit: [f64; 2],
    /// Gas thermal entropy-production factor relative to the isothermal wall.
    /// The physical coefficient is k_g/(2*T0*delta_g)*|T_bulk|^2 per area.
    pub gas_loss_factor: f64,
    /// Solid thermal entropy-production factor in the same normalization.
    /// Total acoustic thermal loss is the sum, without clipping either term.
    pub solid_loss_factor: f64,
}

impl ConductiveThermalWall {
    fn validate(self) -> Result<(), AcousticFeError> {
        let t = self.transport;
        if ![t.conductivity_w_m_k, t.volumetric_heat_capacity_j_m3_k]
            .into_iter()
            .all(|x| x.is_finite() && x > 0.0)
            || self.minimum_tangential_scale_m.is_nan()
            || self.minimum_tangential_scale_m <= 0.0
            || matches!(self.depth, ThermalWallDepth::Slab { thickness_m, .. }
                if !thickness_m.is_finite() || thickness_m <= 0.0)
        {
            return Err(invalid("conductive thermal wall transport or geometry"));
        }
        let alpha = t.conductivity_w_m_k / t.volumetric_heat_capacity_j_m3_k;
        if !alpha.is_finite() || alpha <= 0.0 {
            return Err(invalid("conductive thermal wall diffusivity"));
        }
        Ok(())
    }
}

impl ThermalWall {
    /// Evaluate the exact local normal diffusion profile after both gas and
    /// solid thin-layer admissions. No matrix or material query is allocated.
    pub fn response(
        self,
        gas: BoundaryLayerTransport,
        medium: Medium,
        geometry: BoundaryLayerGeometry,
        frequency_hz: f64,
    ) -> Result<ThermalWallResponse, AcousticFeError> {
        gas.scales(medium, geometry, frequency_hz)?;
        let isothermal = ThermalWallResponse {
            wall_temperature: [0.0, 0.0],
            gas_temperature_deficit: [1.0, 0.0],
            gas_loss_factor: 1.0,
            solid_loss_factor: 0.0,
        };
        let wall = match self {
            Self::Isothermal => return Ok(isothermal),
            Self::Adiabatic => {
                return Ok(ThermalWallResponse {
                    wall_temperature: [1.0, 0.0],
                    gas_temperature_deficit: [0.0, 0.0],
                    gas_loss_factor: 0.0,
                    solid_loss_factor: 0.0,
                });
            }
            Self::Conductive(wall) => wall,
        };
        wall.validate()?;
        let t = wall.transport;
        let omega = 2.0 * core::f64::consts::PI * frequency_hz;
        let delta =
            det::sqrt(2.0 * (t.conductivity_w_m_k / t.volumetric_heat_capacity_j_m3_k) / omega);
        let penetrated = match wall.depth {
            ThermalWallDepth::HalfSpace => delta,
            ThermalWallDepth::Slab { thickness_m, .. } => delta.min(thickness_m),
        };
        let scale = geometry
            .minimum_curvature_radius_m
            .min(wall.minimum_tangential_scale_m);
        if !delta.is_finite()
            || delta <= 0.0
            || penetrated / scale > geometry.maximum_layer_ratio
            || omega / medium.sound_speed_m_s * penetrated > geometry.maximum_layer_ratio
        {
            return Err(invalid(
                "conductive wall normal diffusion violates local-layer geometry",
            ));
        }
        if gas.thermal_conductivity_w_m_k == 0.0 {
            // No gas heat flux can drive the stationary solid; the entire
            // acoustic thermal coefficient is zero in the original operator.
            return Ok(isothermal);
        }
        let eg = det::sqrt(gas.thermal_conductivity_w_m_k)
            * det::sqrt(medium.density_kg_m3)
            * det::sqrt(gas.cp_j_kg_k);
        let ew = det::sqrt(t.conductivity_w_m_k) * det::sqrt(t.volumetric_heat_capacity_j_m3_k);
        let ratio = ew / eg;
        if !ratio.is_finite() || ratio <= 0.0 {
            return Err(invalid("gas/solid thermal effusivity ratio"));
        }
        let y = match wall.depth {
            ThermalWallDepth::HalfSpace => [1.0, 1.0],
            ThermalWallDepth::Slab {
                thickness_m,
                backing,
            } => slab_admittance(thickness_m / delta, backing)?,
        };
        // Scale gas and wall admittances together before division. All real
        // components are passive. This avoids overflow for a very stiff
        // thermal sink and cancellation in thin insulated-wall dissipation.
        let rs = ratio.max(1.0);
        let ys = y[0].max(y[1]).max(1.0);
        let g = (1.0 / rs) / ys;
        let w = [(ratio / rs) * (y[0] / ys), (ratio / rs) * (y[1] / ys)];
        let d = [w[0] + g, w[1] + g];
        let norm = d[0] * d[0] + d[1] * d[1];
        let deficit = [
            (w[0] * d[0] + w[1] * d[1]) / norm,
            (w[1] * d[0] - w[0] * d[1]) / norm,
        ];
        let temperature = [g * (d[0] + d[1]) / norm, g * (d[0] - d[1]) / norm];
        let result = ThermalWallResponse {
            wall_temperature: temperature,
            gas_temperature_deficit: deficit,
            gas_loss_factor: (w[0] * w[0] + w[1] * w[1]) / norm,
            solid_loss_factor: 2.0 * w[0] * g / norm,
        };
        if result
            .wall_temperature
            .iter()
            .chain(&result.gas_temperature_deficit)
            .chain([&result.gas_loss_factor, &result.solid_loss_factor])
            .any(|x| !x.is_finite())
        {
            return Err(invalid("conductive thermal wall response overflow"));
        }
        Ok(result)
    }
}

// Y_w / (k_w/delta_w), with explicitly nonnegative real component. Evaluating
// complex tanh then subtracting its nearly equal parts loses thin-wall loss.
fn slab_admittance(x: f64, backing: ThermalBacking) -> Result<[f64; 2], AcousticFeError> {
    if x.is_nan() || x <= 0.0 {
        return Err(invalid("conductive wall dimensionless depth"));
    }
    let y = if x < 0.1 {
        // Series of sinh(2x) +/- sin(2x) and cosh(2x) +/- cos(2x).
        // Next relative term is below 9e-18 throughout this branch.
        let x2 = x * x;
        let x4 = x2 * x2;
        let plus = 1.0 + x4 * (2.0 / 15.0 + x4 * (2.0 / 2835.0));
        let minus = 1.0 + x4 * (2.0 / 105.0 + x4 * (2.0 / 51975.0));
        match backing {
            ThermalBacking::Adiabatic => {
                let den = 1.0 + x4 * (2.0 / 3.0 + x4 * (2.0 / 315.0));
                [(4.0 / 3.0) * x * x * x * minus / den, 2.0 * x * plus / den]
            }
            ThermalBacking::Isothermal => {
                let den = 1.0 + x4 * (2.0 / 45.0 + x4 * (2.0 / 14175.0));
                [plus / (x * den), (2.0 / 3.0) * x * minus / den]
            }
        }
    } else {
        let e = det::exp(-2.0 * x);
        if e == 0.0 {
            [1.0, 1.0]
        } else {
            let (sin, cos) = det::sin_cos(2.0 * x);
            let a = 1.0 - e * e;
            let b = 2.0 * e * sin;
            match backing {
                ThermalBacking::Adiabatic => {
                    let den = 1.0 + e * e + 2.0 * e * cos;
                    [(a - b) / den, (a + b) / den]
                }
                ThermalBacking::Isothermal => {
                    let den = 1.0 + e * e - 2.0 * e * cos;
                    [(a + b) / den, (a - b) / den]
                }
            }
        }
    };
    if y.into_iter().any(|v| !v.is_finite() || v < 0.0) {
        return Err(invalid(
            "conductive wall admittance overflow or nonpassivity",
        ));
    }
    Ok(y)
}

impl ReducedHelmholtz {
    /// Change the thermal condition on the already prepared viscous/thermal
    /// wall patch. The original projected matrices and pressure basis are
    /// retained; only frequency-local scalar coefficients change. A model
    /// without an attached boundary-layer wall refuses.
    pub fn with_thermal_wall(mut self, wall: ThermalWall) -> Result<Self, AcousticFeError> {
        if let ThermalWall::Conductive(conductive) = wall {
            conductive.validate()?;
        }
        let layer = self
            .boundary_layer
            .as_mut()
            .ok_or_else(|| invalid("thermal wall requires prepared boundary-layer faces"))?;
        layer.thermal_wall = wall;
        Ok(self)
    }

    /// Thermal boundary configuration retained by this original FE reduction.
    #[must_use]
    pub fn thermal_wall(&self) -> Option<ThermalWall> {
        self.boundary_layer.as_ref().map(|layer| layer.thermal_wall)
    }
}

#[cfg(test)]
mod tests;
