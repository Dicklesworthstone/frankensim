//! Equilibrium smooth-wall traction on the original conservative gas faces.
//!
//! The shared Spalding profile supplies the **secant** viscosity for traction;
//! its tangent is not a wall stress. The nonlinear traction Jacobian, derived
//! from the same profile, supplies the explicit diffusion bound. No heat or
//! species wall law is inferred from this momentum closure.

use super::*;
use fs_lbm::conjugate::turbulence::law_of_the_wall_ratios;

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
/// nonequilibrium wall history are absent. Only adiabatic no-slip walls are
/// supported; a molecular temperature jump is not a turbulent heat model.
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
        self.admit_wall(wall)?;
        let MolecularWall3d::NoSlip { velocity: wall, .. } = wall else {
            unreachable!();
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
        sgs.wall_models[index] = Some(config);
        Ok(self)
    }

    /// Exact authored wall controls, independent of the interior SGS model.
    pub fn equilibrium_sgs_wall(&self, side: BoundarySide) -> Option<EquilibriumWallConfig> {
        self.sgs.as_ref()?.wall_models[2 * side.axis() + usize::from(side.sign() > 0.0)]
    }
}

impl SubgridTransport {
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
        config
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
}

#[cfg(test)]
mod tests;
