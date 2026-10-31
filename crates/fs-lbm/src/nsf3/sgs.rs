//! Explicit low-Mach SGS transport for the conservative NSF owner.
//!
//! Only deviatoric momentum relaxation receives `nu_t`; molecular bulk
//! relaxation is unchanged. The same shear rate supplies the stress work in
//! the total-energy flux. Gradient diffusion adds `rho cp nu_t / Pr_t` to
//! interior heat conduction. No second viscous-heating source is introduced.
//! This low-Mach algebraic closure supplies no SGS kinetic-energy equation,
//! isotropic SGS pressure, pressure-dilatation or shocklet model. It does not
//! extend the molecular solver's Mach/stability qualification. The modeled
//! viscosity and Prandtl heat-flux form follows the limited closure in NASA
//! report 20200002378, equations (8)-(10), without its extra compressible terms.

use super::*;
use crate::sgs3::{SgsConfig, SgsError, SgsModel, eddy_viscosity};

/// Every modeling choice is explicit; there is no turbulence default.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NsfSgsConfig {
    pub closure: SgsConfig,
    /// Positive, finite gradient-diffusion Prandtl number. This is a declared
    /// modeling parameter, not the molecular Prandtl number of the gas source.
    pub turbulent_prandtl: f64,
}

impl NsfSgsConfig {
    pub fn validate(self) -> Result<(), NsfError> {
        self.closure.validate().map_err(NsfError::Sgs)?;
        if !(self.turbulent_prandtl.is_finite() && self.turbulent_prandtl > 0.0) {
            return Err(NsfError::Input(
                "SGS turbulent Prandtl number must be finite and positive",
            ));
        }
        Ok(())
    }
}

/// Modeled transport at one accepted fluid cell. Molecular transport remains
/// owned by the original gas law. Multiply `eddy_viscosity_lu` by `dx^2/dt`
/// for m^2/s; multiply `turbulent_conductivity_lu` by
/// `rho_ref R dx^2/dt` for W/(m K).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct NsfSgsTransport {
    pub eddy_viscosity_lu: f64,
    pub turbulent_conductivity_lu: f64,
}

impl<L: Lattice3<27>> NsfGrid3<L> {
    /// Attach an explicit exploratory closure and restart the current
    /// conservative flux history. Periodic faces, stationary solids and
    /// stationary/moving halfway walls, characteristic faces and their localized
    /// transverse reservoir aperture are supported. Complete prescribed-reservoir
    /// faces remain unsupported in either configuration order.
    /// The zero-coefficient control leaves the molecular trajectory unchanged.
    ///
    /// All existing gas/Mach/temperature bounds still apply; the *total*
    /// local kinematic viscosity must also fit the existing viscosity envelope.
    /// That admission bound is not a qualification of variable SGS transport.
    /// No coefficient, excess viscosity or heat flux is silently clamped.
    /// No wall function is supplied: near-wall flow must be resolved. The
    /// modeled heat flux vanishes at an impermeable wall; the source-owned
    /// molecular half-cell heat resistance remains the wall heat-transfer law.
    /// At an open face the adjacent algebraic conductivity is extrapolated
    /// normally, and its heat flux uses the original face/ghost temperature.
    /// This is an explicit zero-normal-gradient transport boundary rung, not a
    /// supplied turbulent inflow spectrum or a transparency qualification.
    pub fn with_sgs(mut self, config: NsfSgsConfig) -> Result<Self, NsfError> {
        config.validate()?;
        let coefficient = match config.closure.model {
            SgsModel::Smagorinsky { cs } => cs,
            SgsModel::Wale { cw } => cw,
            SgsModel::Vreman { cv } => cv,
        };
        if coefficient == 0.0 {
            if self.sgs.is_some() {
                return Err(NsfError::Input(
                    "attach the zero SGS control to a molecular grid",
                ));
            }
            return Ok(self);
        }
        if self.reservoir_faces.iter().flatten().any(Option::is_some) {
            return Err(NsfError::Input(
                "NSF SGS does not support complete prescribed-reservoir faces",
            ));
        }
        if self.characteristic_thermal_policy.contains(&NsfCharacteristicThermalPolicy::EnthalpyFlux) {
            return Err(NsfError::Input(
                "enthalpy-flux characteristic thermal contacts do not yet support active SGS",
            ));
        }
        self.sgs = Some(config);
        let (rho, mom) = self.macro_fields(&self.f);
        let velocity = self.velocity_field(&rho, &mom, &self.force);
        let transport = self.sgs_field(
            &rho,
            &velocity,
            &self.theta,
            &self.material,
            self.open_faces.as_ref(),
        )?;
        let stress = self.stress_field(
            &self.f,
            &self.force,
            &rho,
            &self.theta,
            &self.material,
            &transport,
        );
        let (div, wall, open) = self.flux_divergence(
            &self.e,
            &rho,
            &velocity,
            &self.theta,
            &stress,
            &transport,
            self.open_faces.as_ref(),
        )?;
        self.sgs_transport = transport;
        self.div = div;
        self.wall_out = wall;
        self.open_out = open;
        Ok(self)
    }

    /// Active controls; `None` also denotes an explicitly selected exact-zero
    /// coefficient, which uses the unchanged molecular arithmetic.
    #[must_use]
    pub const fn sgs_config(&self) -> Option<NsfSgsConfig> {
        self.sgs
    }

    /// Accepted-state modeled coefficients; `None` when SGS is disabled.
    /// Out-of-domain and solid cells refuse rather than reporting gas values.
    pub fn sgs_transport(&self, x: [usize; 3]) -> Result<Option<NsfSgsTransport>, NsfError> {
        if x.iter().zip(self.spec.dims).any(|(&i, n)| i >= n) {
            return Err(NsfError::Input("SGS transport cell outside the grid"));
        }
        let i = self.index(x);
        if self.solid[i] {
            return Err(NsfError::Input("solid cell has no gas SGS transport"));
        }
        Ok(self.sgs_transport.get(i).copied())
    }

    pub(super) fn sgs_rates(
        &self,
        rho: f64,
        theta: f64,
        material: Option<&NsfCellMaterial>,
        transport: Option<&NsfSgsTransport>,
    ) -> (f64, f64) {
        let molecular = self.rates(rho, theta, material);
        let Some(transport) = transport.filter(|t| t.eddy_viscosity_lu != 0.0) else {
            return molecular;
        };
        let mu = material.map_or_else(|| self.law.mu_lattice(theta), |m| m.mu_lattice);
        (
            1.0 / (0.5 + (mu / rho + transport.eddy_viscosity_lu) / THETA0),
            molecular.1,
        )
    }

    // One flat field per predictor; no per-cell allocation or gas-source
    // requery. Existing force-corrected velocity and local material caches are
    // reused. Fixed per-cell arithmetic is independent of worker count.
    pub(super) fn sgs_field(
        &self,
        rho: &[f64],
        velocity: &[[f64; 3]],
        theta: &[f64],
        material: &[Option<NsfCellMaterial>],
        faces: Option<&OpenFaces>,
    ) -> Result<Vec<NsfSgsTransport>, NsfError> {
        let Some(config) = self.sgs else {
            return Ok(Vec::new());
        };
        let mut out = vec![Ok(NsfSgsTransport::default()); self.n];
        par_fill(&mut out, self.spec.workers, &|i| {
            if self.solid[i] {
                return Ok(NsfSgsTransport::default());
            }
            let local = material.get(i).and_then(Option::as_ref);
            let mu = local.map_or_else(|| self.law.mu_lattice(theta[i]), |m| m.mu_lattice);
            let molecular = mu / rho[i];
            if !(molecular.is_finite() && molecular > 0.0 && rho[i] > 0.0) {
                return Err(NsfError::Sgs(SgsError::InvalidMolecularViscosity {
                    cell: i,
                }));
            }
            let eddy = eddy_viscosity(
                config.closure.model,
                config.closure.filter_width_lu,
                self.sgs_gradient(i, velocity, faces),
            )
            .map_err(NsfError::Sgs)?;
            let ratio = eddy / molecular;
            if ratio > config.closure.maximum_eddy_to_molecular {
                return Err(NsfError::Sgs(SgsError::ViscosityRatioExceeded {
                    cell: i,
                    ratio,
                    maximum: config.closure.maximum_eddy_to_molecular,
                }));
            }
            if molecular + eddy > MAX_VERIFIED_LATTICE_VISCOSITY {
                return Err(NsfError::Viscosity {
                    observed: molecular + eddy,
                    maximum: MAX_VERIFIED_LATTICE_VISCOSITY,
                });
            }
            let cp_over_r = local.map_or(self.law.gamma / (self.law.gamma - 1.0), |m| {
                m.state.specific_heat_cp / self.law.gas_constant
            });
            let conductivity = rho[i] * cp_over_r * eddy / config.turbulent_prandtl;
            if !(conductivity.is_finite() && conductivity >= 0.0) {
                return Err(NsfError::Input("nonfinite SGS turbulent heat conductivity"));
            }
            Ok(NsfSgsTransport {
                eddy_viscosity_lu: eddy,
                turbulent_conductivity_lu: conductivity,
            })
        });
        out.into_iter().collect()
    }

    // Second-order centred differences in fluid or across a full-spacing
    // prescribed ghost. At a halfway wall/characteristic face use its actual
    // velocity, this centre and the opposite fluid centre: the quadratic
    // derivative is exact. An open face is never treated as a no-slip wall.
    fn sgs_gradient(
        &self,
        i: usize,
        velocity: &[[f64; 3]],
        faces: Option<&OpenFaces>,
    ) -> [[f64; 3]; 3] {
        let x = self.at(i);
        let wall = |axis, side| {
            if let Some(face) = self.characteristic_face(faces, i, axis, side) {
                face.u
            } else if self.shifted(x, axis, side).is_none() {
                self.edge_velocity[axis][usize::from(side > 0)]
            } else {
                [0.0; 3]
            }
        };
        let mut gradient = [[0.0; 3]; 3];
        for axis in 0..3 {
            let low = self.fluid_at(x, axis, -1);
            let high = self.fluid_at(x, axis, 1);
            let wall_low = if low.is_none() {
                wall(axis, -1)
            } else {
                [0.0; 3]
            };
            let wall_high = if high.is_none() {
                wall(axis, 1)
            } else {
                [0.0; 3]
            };
            for component in 0..3 {
                gradient[component][axis] = match (low, high) {
                    (Some(low), Some(high)) => {
                        0.5 * (velocity[high][component] - velocity[low][component])
                    }
                    (None, Some(high)) => {
                        if let Some(ghost) = self.open_at(i, axis, -1) {
                            0.5 * (velocity[high][component] - ghost.u[component])
                        } else {
                            velocity[i][component] + velocity[high][component] / 3.0
                                - (4.0 / 3.0) * wall_low[component]
                        }
                    }
                    (Some(low), None) => {
                        if let Some(ghost) = self.open_at(i, axis, 1) {
                            0.5 * (ghost.u[component] - velocity[low][component])
                        } else {
                            (4.0 / 3.0) * wall_high[component]
                                - velocity[i][component]
                                - velocity[low][component] / 3.0
                        }
                    }
                    (None, None) => wall_high[component] - wall_low[component],
                };
            }
        }
        gradient
    }
}

#[cfg(test)]
mod tests;
#[cfg(test)]
mod vreman_tests;
