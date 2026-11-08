//! Equilibrium smooth-wall traction on the original conservative gas faces.
//!
//! The shared Spalding profile supplies the **secant** viscosity for traction;
//! its tangent is not a wall stress. The nonlinear traction Jacobian, derived
//! from the same profile, supplies the explicit diffusion bound. An explicit
//! thermal model integrates the same profile with the authored turbulent Pr.
//! No heat or species wall law is inferred from momentum controls alone.

use super::*;
use fs_gas::multicomponent::dissipation::DissipationError;
use fs_lbm::conjugate::{ChtError, FluidProperties};
use fs_lbm::conjugate::turbulence::law_of_the_wall_ratios;
pub use fs_lbm::conjugate::turbulence::thermal_wall::WallHeatConfig;
use fs_lbm::conjugate::turbulence::thermal_wall::equilibrium_wall_heat_transfer_cancellable;

/// Explicit smooth-wall equilibrium stress for an SGS gas calculation.
///
/// Uses the shared Spalding profile with `kappa = 0.41`, `E = 8.6`. The
/// physical half-cell distance and instantaneous wall-relative velocity fix
/// `Re_y = rho |u_t-u_wall| y / mu`. This profile connects the viscous sublayer
/// continuously to the logarithmic layer without fitting a new coefficient.
///
/// Admit only attached, approximately equilibrium, smooth-wall, low-Mach
/// flow. The caller must establish that the adjacent sample lies in the
/// inner layer; the numerical `y+` bound alone does not establish this.
/// Roughness, pressure-gradient/separation corrections, transition and
/// nonequilibrium wall history are absent. Used alone this admits adiabatic
/// no-slip walls. [`EquilibriumThermalWallConfig`] explicitly adds the thermal
/// law needed for an isothermal wall with an unresolved viscous sublayer.
///
/// Reference: Spalding (1961), *A single formula for the "law of the wall"*,
/// Journal of Applied Mechanics 28, 455--458, doi:10.1115/1.3641728.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquilibriumWallConfig {
    /// Inclusive admitted `y u_tau / nu`, positive and at most 1000.
    /// Larger values refuse the complete trial rather than clipping stress.
    pub maximum_wall_distance_plus: f64,
    /// Inclusive wall-relative Mach ceiling, positive and at most 0.3.
    /// Includes the normal velocity when checking the low-Mach assumption.
    pub maximum_mach_number: f64,
}

impl EquilibriumWallConfig {
    pub fn validate(self) -> Result<(), &'static str> {
        if !positive(self.maximum_wall_distance_plus)
            || self.maximum_wall_distance_plus > 1000.0
            || !positive(self.maximum_mach_number)
            || self.maximum_mach_number > 0.3
        {
            return Err("invalid equilibrium SGS wall applicability bounds");
        }
        Ok(())
    }

    pub(in crate::compressible::spatial) fn admit_wall(
        self,
        wall: MolecularWall3d,
    ) -> Result<(), &'static str> {
        self.validate()?;
        if !matches!(
            wall,
            MolecularWall3d::NoSlip {
                thermal: WallThermalCondition::Adiabatic,
                ..
            }
        ) {
            return Err("equilibrium SGS stress requires an adiabatic no-slip wall");
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub(in crate::compressible::spatial) fn viscosity(
        self,
        wall: MolecularWall3d,
        normal: usize,
        velocity: [f64; 3],
        density: f64,
        sound_speed: f64,
        distance: f64,
        molecular: f64,
    ) -> Result<WallViscosity, &'static str> {
        self.validate()?;
        let MolecularWall3d::NoSlip { velocity: wall, .. } = wall else {
            return Err("equilibrium SGS stress requires a no-slip wall");
        };
        if ![density, sound_speed, distance, molecular]
            .into_iter()
            .all(positive)
        {
            return Err("equilibrium SGS wall requires positive source density and viscosity");
        }
        let relative = std::array::from_fn(|axis| velocity[axis] - wall[axis]);
        if relative.iter().any(|value: &f64| !value.is_finite()) {
            return Err("unrepresentable equilibrium SGS wall relative velocity");
        }
        if norm(relative) / sound_speed > self.maximum_mach_number {
            return Err("equilibrium SGS wall relative Mach bound exceeded");
        }
        let speed = norm(std::array::from_fn(|axis| {
            if axis == normal { 0.0 } else { relative[axis] }
        }));
        if speed == 0.0 {
            return Ok(WallViscosity {
                secant: molecular,
                tangent: molecular,
            });
        }
        let reynolds = density * speed * distance / molecular;
        if !reynolds.is_finite() {
            return Err("unrepresentable equilibrium SGS wall Reynolds number");
        }
        // y+ >= u+, hence Re_y <= (y+)^2. This necessary domain bound also
        // keeps the shared inversion away from any extreme exponential range.
        if reynolds > self.maximum_wall_distance_plus.powi(2) {
            return Err("equilibrium SGS wall distance-plus bound exceeded");
        }
        let (eddy_tangent, secant) = law_of_the_wall_ratios(reynolds);
        let y_plus = (reynolds * secant).sqrt();
        if y_plus > self.maximum_wall_distance_plus {
            return Err("equilibrium SGS wall distance-plus bound exceeded");
        }
        let slope = 1.0 + eddy_tangent;
        // d(tau)/dU * y/mu = 2 r s/(r+s), where r=y+/u+ and s=dy+/du+.
        // This is the largest eigenvalue of the tangential traction Jacobian;
        // the orthogonal tangential direction uses the secant. Using r alone
        // underestimates explicit stiffness as the logarithmic layer develops.
        let tangent = (2.0 * secant / (1.0 + secant / slope)).max(secant);
        let result = WallViscosity {
            secant: molecular * secant,
            tangent: molecular * tangent,
        };
        if ![y_plus, result.secant, result.tangent]
            .into_iter()
            .all(f64::is_finite)
            || result.secant < molecular
        {
            return Err("unrepresentable equilibrium SGS wall traction");
        }
        Ok(result)
    }
}

/// Equilibrium smooth-wall stress and heat for a prescribed-temperature wall.
///
/// Integrates `dT+/dy+ = Pr / (1 + Pr/Pr_t * nu_t/nu)` using the same Spalding
/// profile as momentum, including the molecular sublayer and buffer layer.
/// Molecular viscosity/conductivity and source Cp are resampled at the
/// arithmetic fluid/wall temperature on every stage; density and velocity
/// retain the same adjacent sample as momentum. Properties are held constant
/// across the wall integration at the current composition. Use only for small
/// property variations and the low-Mach attached equilibrium regime admitted
/// by `momentum`; strong heating, compressibility and viscous-heating wall
/// corrections are absent. The interior SGS heat closure remains independent.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquilibriumThermalWallConfig {
    pub momentum: EquilibriumWallConfig,
    /// Explicit turbulent Prandtl number and bounded thermal quadrature.
    pub heat: WallHeatConfig,
}

impl EquilibriumThermalWallConfig {
    pub fn validate(self) -> Result<(), &'static str> {
        self.momentum.validate()?;
        self.heat
            .validate()
            .map_err(|_| "invalid equilibrium SGS thermal wall controls")
    }
}

/// One optional wall model; no new cell/face workspace is needed for heat.
#[derive(Clone, Copy)]
pub(in crate::compressible::spatial) struct EquilibriumWallModel {
    pub(in crate::compressible::spatial) momentum: EquilibriumWallConfig,
    pub(in crate::compressible::spatial) heat: Option<WallHeatConfig>,
}

impl From<EquilibriumWallConfig> for EquilibriumWallModel {
    fn from(momentum: EquilibriumWallConfig) -> Self {
        Self { momentum, heat: None }
    }
}

impl From<EquilibriumThermalWallConfig> for EquilibriumWallModel {
    fn from(config: EquilibriumThermalWallConfig) -> Self {
        Self { momentum: config.momentum, heat: Some(config.heat) }
    }
}

impl EquilibriumWallModel {
    pub(in crate::compressible::spatial) fn thermal_config(self) -> Option<EquilibriumThermalWallConfig> {
        Some(EquilibriumThermalWallConfig { momentum: self.momentum, heat: self.heat? })
    }

    pub(in crate::compressible::spatial) fn admit_wall(self, wall: MolecularWall3d) -> Result<(), &'static str> {
        if let Some(config) = self.thermal_config() {
            config.validate()?;
            if !matches!(wall, MolecularWall3d::NoSlip { thermal: WallThermalCondition::Isothermal(_), .. }) {
                return Err("equilibrium SGS thermal transfer requires an isothermal no-slip wall");
            }
            Ok(())
        } else {
            self.momentum.admit_wall(wall)
        }
    }

    pub(in crate::compressible::spatial) fn conductivity(
        self,
        wall: MolecularWall3d,
        normal: usize,
        velocity: [f64; 3],
        distance: f64,
        fluid: &FluidProperties,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, DissipationError> {
        let Some(heat) = self.heat else { return Ok(fluid.conductivity_w_m_k); };
        self.admit_wall(wall).map_err(DissipationError::Invalid)?;
        let MolecularWall3d::NoSlip { velocity: wall_velocity, .. } = wall else { unreachable!(); };
        let speed = norm(std::array::from_fn(|axis| {
            if axis == normal { 0.0 } else { velocity[axis] - wall_velocity[axis] }
        }));
        let transfer = equilibrium_wall_heat_transfer_cancellable(speed, distance, fluid, heat, cancelled)
            .map_err(|error| match error {
                ChtError::Cancelled => DissipationError::Cancelled,
                _ => DissipationError::Invalid("equilibrium SGS thermal wall properties or quadrature refused"),
            })?;
        if transfer.y_plus > self.momentum.maximum_wall_distance_plus {
            return Err(DissipationError::Invalid("equilibrium SGS thermal wall distance-plus bound exceeded"));
        }
        // Return the exact molecular coefficient at rest. Elsewhere constrain
        // only roundoff in the nonnegative eddy contribution; failed quadrature
        // has already refused rather than falling back to molecular heat.
        let conductivity = if speed == 0.0 { fluid.conductivity_w_m_k }
            else { (distance * transfer.heat_transfer_w_m2_k).max(fluid.conductivity_w_m_k) };
        if !positive(conductivity) {
            return Err(DissipationError::Invalid("unrepresentable equilibrium SGS wall conductivity"));
        }
        Ok(conductivity)
    }
}

fn norm(value: [f64; 3]) -> f64 {
    let scale = value.into_iter().map(f64::abs).fold(0.0_f64, f64::max);
    if scale == 0.0 {
        return 0.0;
    }
    scale
        * value
            .into_iter()
            .map(|v| (v / scale).powi(2))
            .sum::<f64>()
            .sqrt()
}

#[derive(Debug)]
pub(in crate::compressible::spatial) struct WallViscosity {
    /// Coefficient of the actual half-cell shear traction.
    pub(in crate::compressible::spatial) secant: f64,
    /// Upper eigenvalue of d(traction)/d(relative velocity), times distance.
    pub(in crate::compressible::spatial) tangent: f64,
}

impl<L: TransportLaw> MolecularTransport3d<L> {
    /// Add an explicitly bounded equilibrium traction to one already authored
    /// adiabatic no-slip wall. The original SGS interior and wall-work flux
    /// remain the owners; no additional wall workspace or heat source is added.
    pub fn with_equilibrium_sgs_wall(
        mut self,
        side: BoundarySide,
        config: EquilibriumWallConfig,
    ) -> Result<Self, MolecularError> {
        let index = 2 * side.axis() + usize::from(side.sign() > 0.0);
        config
            .admit_wall(self.walls[index])
            .map_err(MolecularError::Invalid)?;
        let sgs = self.sgs.as_mut().ok_or(MolecularError::Invalid(
            "equilibrium wall stress requires an explicit SGS operator",
        ))?;
        sgs.wall_models[index] = Some(config.into());
        Ok(self)
    }

    /// Model both traction and heat on one already authored isothermal no-slip
    /// wall. The resulting conductance enters the same conservative energy
    /// face and parabolic admission as molecular conduction; no heat source or
    /// new workspace is introduced. Zero molecular viscosity/conductivity,
    /// domain exit, quadrature exhaustion and cancellation refuse the trial.
    pub fn with_equilibrium_sgs_thermal_wall(
        mut self,
        side: BoundarySide,
        config: EquilibriumThermalWallConfig,
    ) -> Result<Self, MolecularError> {
        let index = 2 * side.axis() + usize::from(side.sign() > 0.0);
        let model = EquilibriumWallModel::from(config);
        model.admit_wall(self.walls[index]).map_err(MolecularError::Invalid)?;
        let sgs = self.sgs.as_mut().ok_or(MolecularError::Invalid(
            "equilibrium thermal wall requires an explicit SGS operator",
        ))?;
        sgs.wall_models[index] = Some(model);
        Ok(self)
    }

    /// Exact authored wall controls, independent of the interior SGS model.
    pub fn equilibrium_sgs_wall(&self, side: BoundarySide) -> Option<EquilibriumWallConfig> {
        Some(self.sgs.as_ref()?.wall_models[2 * side.axis() + usize::from(side.sign() > 0.0)]?.momentum)
    }

    pub fn equilibrium_sgs_thermal_wall(&self, side: BoundarySide) -> Option<EquilibriumThermalWallConfig> {
        self.sgs.as_ref()?.wall_models[2 * side.axis() + usize::from(side.sign() > 0.0)]?.thermal_config()
    }
}

impl SubgridTransport {
    pub(in crate::compressible::spatial::molecular) fn has_equilibrium_wall(&self) -> bool {
        self.wall_models.iter().any(Option::is_some)
    }

    pub(in crate::compressible::spatial::molecular) fn admit_wall(
        &self,
        side: BoundarySide,
        wall: MolecularWall3d,
    ) -> Result<(), MolecularError> {
        if let Some(config) = self.wall_models[2 * side.axis() + usize::from(side.sign() > 0.0)] {
            config.admit_wall(wall).map_err(MolecularError::Invalid)?;
        }
        Ok(())
    }

    pub(in crate::compressible::spatial::molecular) fn wall_viscosity(
        &self,
        side: usize,
        wall: MolecularWall3d,
        primitive: Primitive,
        distance: f64,
        molecular: f64,
    ) -> Result<WallViscosity, MolecularError> {
        let Some(config) = self.wall_models[side] else {
            return Ok(WallViscosity {
                secant: molecular,
                tangent: molecular,
            });
        };
        config.admit_wall(wall).map_err(MolecularError::Invalid)?;
        config
            .momentum
            .viscosity(
                wall,
                side / 2,
                primitive.velocity,
                primitive.density,
                primitive.sound_speed,
                distance,
                molecular,
            )
            .map_err(MolecularError::Invalid)
    }

    #[allow(clippy::too_many_arguments)]
    pub(in crate::compressible::spatial::molecular) fn wall_conductivity(
        &self,
        gas: &FrozenGas,
        side: usize,
        wall: MolecularWall3d,
        primitive: Primitive,
        distance: f64,
        molecular_viscosity: f64,
        molecular_conductivity: f64,
        temperature: f64,
        pressure: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, MolecularError> {
        let Some(model) = self.wall_models[side].filter(|model| model.heat.is_some()) else {
            return Ok(molecular_conductivity);
        };
        let cp = gas
            .source()
            .evaluate_properties(Temperature::new(temperature), Pressure::new(pressure))
            .map_err(Error::from)?
            .mass_specific_cp()
            .value();
        Ok(model.conductivity(
            wall,
            side / 2,
            primitive.velocity,
            distance,
            &FluidProperties {
                density_kg_m3: primitive.density,
                specific_heat_j_kg_k: cp,
                conductivity_w_m_k: molecular_conductivity,
                kinematic_viscosity_m2_s: molecular_viscosity / primitive.density,
            },
            cancelled,
        )?)
    }
}

#[cfg(test)]
mod tests;
