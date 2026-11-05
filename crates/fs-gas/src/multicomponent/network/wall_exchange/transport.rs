//! Material-dependent developed pipe-wall losses in the existing wall transaction.

use super::{NetworkWallOperator, WallContact, WallDrag, WallExchangeError, WallThermal, poll};
use crate::multicomponent::dissipation::TransportLaw;
use crate::multicomponent::dissipation::evolve::fractions;
use crate::multicomponent::{Cell, MAX_SPECIES, Mesh1d, MixtureGas, allocated, positive};
use fs_qty::{Length, Pressure, Temperature};

mod churchill;
mod turbulent;
pub use turbulent::{CircularChurchillWall, CircularGnielinskiWall, CircularTurbulentWall};

/// Explicit circular, developed-laminar approximation for one complete branch.
///
/// The diameter must agree with the original mesh's constant circular area and
/// integrated volumes. The full wetted area of a cell is therefore `4 V / D`;
/// callers cannot accidentally supply the entire pipe area to every cell.
///
/// The Nusselt number is an authored correlation input, not inferred from a gas
/// label. For example, a fully developed constant-wall-temperature approximation
/// may declare 3.66 with its own admitted regime. This type does not establish
/// entrance length, negligible axial conduction, fully developed flow, or a
/// transient boundary-layer model. Zero Reynolds is admitted only if explicitly
/// included in the caller's approximation interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CircularLaminarWall {
    /// Positive diameter of the actual circular branch [m].
    pub diameter: Length,
    /// Nonnegative declared `h D / k`; zero explicitly disables wall heat flow.
    pub nusselt_number: f64,
    /// Inclusive supported Reynolds interval, contained in `[0, 2300]`.
    pub reynolds_number: [f64; 2],
    /// Positive upper Mach bound, at most 0.3 for this low-Mach wall closure.
    pub maximum_mach_number: f64,
}

impl CircularLaminarWall {
    fn admit_mesh(self, mesh: &Mesh1d) -> Result<(), WallExchangeError> {
        let [minimum, maximum] = self.reynolds_number;
        if !positive(self.diameter.value())
            || !self.nusselt_number.is_finite()
            || self.nusselt_number < 0.0
            || !minimum.is_finite()
            || minimum < 0.0
            || !positive(maximum)
            || minimum > maximum
            || maximum > 2300.0
            || !positive(self.maximum_mach_number)
            || self.maximum_mach_number > 0.3
        {
            return Err(WallExchangeError::Invalid(
                "circular wall geometry or applicability",
            ));
        }
        admit_circular_mesh(self.diameter, mesh)
    }
}

/// An explicitly selected developed-flow model on one circular branch.
///
/// Laminar and turbulent branches can share one original transport owner.
/// Regimes are never selected automatically. Churchill explicitly selects a
/// continuous engineering relation; the other models retain their own domains.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CircularWallModel {
    /// Circular Poiseuille drag and an authored Nusselt number.
    Laminar(CircularLaminarWall),
    /// Colebrook-White drag and an authored Nusselt number.
    Turbulent(CircularTurbulentWall),
    /// Smooth-pipe drag and heat transfer from actual Reynolds and Prandtl numbers.
    Gnielinski(CircularGnielinskiWall),
    /// Explicit all-regime pipe friction and an independently authored Nusselt number.
    Churchill(CircularChurchillWall),
}

impl CircularWallModel {
    fn admit_mesh(self, mesh: &Mesh1d) -> Result<(), WallExchangeError> {
        match self {
            Self::Laminar(model) => model.admit_mesh(mesh),
            Self::Turbulent(model) => {
                model.validate()?;
                admit_circular_mesh(model.diameter, mesh)
            }
            Self::Churchill(model) => {
                model.validate()?;
                admit_circular_mesh(model.diameter, mesh)
            }
            Self::Gnielinski(model) => {
                model.validate()?;
                admit_circular_mesh(model.diameter, mesh)?;
                let positions = mesh.positions();
                let length = positions[positions.len() - 1] - positions[0];
                admit(
                    "pipe length / diameter",
                    length / model.diameter.value(),
                    10.0,
                    1e6,
                )
            }
        }
    }

    fn parameters(self) -> (Length, [f64; 2], f64) {
        match self {
            Self::Laminar(model) => (
                model.diameter,
                model.reynolds_number,
                model.maximum_mach_number,
            ),
            Self::Turbulent(model) => (
                model.diameter,
                model.reynolds_number,
                model.maximum_mach_number,
            ),
            Self::Gnielinski(model) => (
                model.diameter,
                model.reynolds_number,
                model.maximum_mach_number,
            ),
            Self::Churchill(model) => (
                model.diameter,
                model.reynolds_number,
                model.maximum_mach_number,
            ),
        }
    }
}

fn admit_circular_mesh(diameter: Length, mesh: &Mesh1d) -> Result<(), WallExchangeError> {
    let area = core::f64::consts::FRAC_PI_4 * diameter.value() * diameter.value();
    let Some(mesh_area) = mesh.uniform_cross_section() else {
        return Err(WallExchangeError::Invalid(
            "developed circular wall requires constant area",
        ));
    };
    if !positive(area) || !same_geometry(area, mesh_area) {
        return Err(WallExchangeError::Invalid(
            "wall diameter differs from original branch area",
        ));
    }
    for (faces, &volume) in mesh.positions().windows(2).zip(mesh.volumes()) {
        if !same_geometry(volume, mesh_area * (faces[1] - faces[0])) {
            return Err(WallExchangeError::Invalid(
                "circular wall volume differs from original geometry",
            ));
        }
    }
    Ok(())
}

pub(super) struct CircularWallTransport {
    law: Box<dyn TransportLaw + Send + Sync>,
    branches: Vec<Option<CircularWallModel>>,
    binary: Vec<f64>,
}

impl NetworkWallOperator<'_> {
    /// Bind original transport data to selected developed circular branches.
    ///
    /// `Some` replaces that branch's authored drag and conductance coefficients
    /// with `32 mu/(rho D²)` and `Nu k (4 V/D)/D`. Wall temperatures, finite heat
    /// capacities, and adiabatic choices retain their existing ownership. `None`
    /// leaves the authored branch unchanged; plenums retain their explicit
    /// contacts because they have no developed axial flow.
    ///
    /// The same original `TransportLaw` used by molecular flow supplies current
    /// viscosity and conductivity. In particular, `MaterialTransport` retains
    /// source cards and mixture rules without copying gas-property scalars here.
    /// Full original species/caloric identity must match, not only species names.
    ///
    /// Every source half resolves the actual conservative mixture's T, p and
    /// mole fractions. Coefficients are frozen over that half (first-order local
    /// constitutive lag), and both endpoints must satisfy the original law and
    /// this explicit applicability domain. The ordinary `advance_for`, reacting,
    /// dissipating and source-composed network steps all use this binding.
    ///
    /// Bind before the first accepted wall step. This consumes the law but no
    /// gas/wall state; it does not silently change an existing accepted report.
    pub fn with_circular_laminar_transport<L: TransportLaw + Send + Sync + 'static>(
        self,
        law: L,
        branches: Vec<Option<CircularLaminarWall>>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, WallExchangeError> {
        poll(&mut cancelled)?;
        let mut models = allocated(branches.len(), None)?;
        for (model, branch) in models.iter_mut().zip(branches) {
            poll(&mut cancelled)?;
            *model = branch.map(CircularWallModel::Laminar);
        }
        self.with_circular_transport(law, models, cancelled)
    }

    /// Resolve laminar or turbulent wall friction from each branch's actual gas.
    ///
    /// The original transport law supplies current viscosity and conductivity;
    /// the selected model supplies geometry, roughness or the smooth-pipe heat
    /// relation, and an explicit applicability domain. The same wall transaction retains heat,
    /// moving-wall work, finite-support recoil and whole-network rollback.
    /// Coefficients are frozen per source interval and checked at both endpoints.
    pub fn with_circular_transport<L: TransportLaw + Send + Sync + 'static>(
        mut self,
        law: L,
        branches: Vec<Option<CircularWallModel>>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, WallExchangeError> {
        poll(&mut cancelled)?;
        if self.transport.is_some() || self.last_report.is_some() || self.trial_active {
            return Err(WallExchangeError::Invalid(
                "wall transport must bind before evolution",
            ));
        }
        if law.source() != self.gas.source()
            || branches.len() != self.branches.len()
            || !branches.iter().any(Option::is_some)
        {
            return Err(WallExchangeError::Invalid(
                "wall transport original source or branch axis",
            ));
        }
        for (model, shape) in branches.iter().zip(&self.branches) {
            poll(&mut cancelled)?;
            if let Some(model) = model {
                model.admit_mesh(&shape.mesh)?;
            }
        }
        let n = self.gas.species_count();
        let binary = allocated(n * n, 0.0)?;
        poll(&mut cancelled)?;
        self.transport = Some(CircularWallTransport {
            law: Box::new(law),
            branches,
            binary,
        });
        Ok(self)
    }

    /// Original immutable constitutive owner, including its species source.
    pub fn transport_law(&self) -> Option<&dyn TransportLaw> {
        self.transport
            .as_ref()
            .map(|transport| transport.law.as_ref() as &dyn TransportLaw)
    }

    /// Original developed-wall models, geometry and applicability, in branch order.
    pub fn circular_walls(&self) -> Option<&[Option<CircularWallModel>]> {
        self.transport
            .as_ref()
            .map(|transport| transport.branches.as_slice())
    }
}

impl CircularWallTransport {
    pub(super) fn resolve(
        &mut self,
        gas: &MixtureGas,
        branch: usize,
        state: Cell<'_>,
        volume: f64,
        mut contact: WallContact,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<WallContact, WallExchangeError> {
        let Some(model) = self.branches[branch] else {
            return Ok(contact);
        };
        poll(cancelled)?;
        if self.law.source() != gas.source() {
            return Err(WallExchangeError::Invalid(
                "wall transport changed its original source",
            ));
        }
        let primitive = gas.primitive(state, &mut *cancelled)?;
        let mut mole_fractions = [0.0; MAX_SPECIES];
        let x = &mut mole_fractions[..gas.species_count()];
        fractions(gas, state.partial_densities, x)?;
        let coefficients = self.law.evaluate(
            Temperature::new(primitive.temperature),
            Pressure::new(primitive.pressure),
            x,
            &mut self.binary,
            cancelled,
        )?;
        let viscosity = coefficients.dynamic_viscosity.value();
        let conductivity = coefficients.thermal_conductivity.value();
        if !positive(viscosity) || !conductivity.is_finite() || conductivity < 0.0 {
            return Err(WallExchangeError::Invalid(
                "wall law viscosity or conductivity",
            ));
        }
        let (diameter, reynolds_number, maximum_mach_number) = model.parameters();
        let motion = contact.drag.moving_velocity();
        let slip = primitive.velocity - motion.unwrap_or(0.0);
        let reynolds = primitive.density * slip.abs() * diameter.value() / viscosity;
        let mach = slip.abs() / primitive.sound_speed;
        let [minimum, maximum] = reynolds_number;
        admit("Reynolds number", reynolds, minimum, maximum)?;
        admit("Mach number", mach, 0.0, maximum_mach_number)?;
        let (friction_factor, nusselt_number) = match model {
            CircularWallModel::Laminar(laminar) => (None, laminar.nusselt_number),
            CircularWallModel::Churchill(churchill) => (None, churchill.nusselt_number),
            CircularWallModel::Turbulent(turbulent) => (
                Some(turbulent.darcy_factor(reynolds)?),
                turbulent.nusselt_number,
            ),
            CircularWallModel::Gnielinski(thermal) => {
                // Cp, not Cv, defines Pr. Use the original NASA mixture at the
                // current composition and temperature, preserving formation-
                // reference independence and avoiding a constant-gamma proxy.
                let (_, cv_density) =
                    gas.caloric(state.partial_densities, primitive.temperature)?;
                let cp =
                    (cv_density + primitive.pressure / primitive.temperature) / primitive.density;
                let prandtl = viscosity * cp / conductivity;
                let (factor, nu) = thermal.coefficients(reynolds, prandtl)?;
                (Some(factor), nu)
            }
        };
        contact.drag = match model {
            CircularWallModel::Churchill(churchill) => {
                let WallDrag::Linear { rate_per_s } = WallDrag::circular_laminar(
                    coefficients.dynamic_viscosity,
                    primitive.density,
                    diameter,
                )?
                else {
                    unreachable!("Poiseuille resistance is linear")
                };
                let excess = churchill.excess(reynolds)?;
                // Retain the exact laminar term. The nonnegative remainder
                // satisfies a*u + b*u*|u| = f_D(Re)*u*|u|/(2D) at this state,
                // while remaining finite and exactly linear at rest.
                let quadratic_per_m = if reynolds == 0.0 {
                    0.0
                } else {
                    32.0 * (excess / reynolds) / diameter.value()
                };
                if !quadratic_per_m.is_finite() {
                    return Err(WallExchangeError::Invalid(
                        "unrepresentable Churchill resistance",
                    ));
                }
                WallDrag::LinearQuadratic {
                    linear_rate_per_s: rate_per_s,
                    quadratic_per_m,
                    wall_velocity_m_s: motion.unwrap_or(0.0),
                }
            }
            CircularWallModel::Laminar(_) => {
                let drag = WallDrag::circular_laminar(
                    coefficients.dynamic_viscosity,
                    primitive.density,
                    diameter,
                )?;
                match (drag, motion) {
                    (WallDrag::Linear { rate_per_s }, Some(wall_velocity_m_s)) => {
                        WallDrag::MovingLinear {
                            rate_per_s,
                            wall_velocity_m_s,
                        }
                    }
                    _ => drag,
                }
            }
            CircularWallModel::Turbulent(_) | CircularWallModel::Gnielinski(_) => {
                let friction_factor =
                    friction_factor.expect("turbulent branch evaluated its factor");
                match motion {
                    Some(wall_velocity_m_s) => WallDrag::MovingDarcy {
                        friction_factor,
                        hydraulic_diameter_m: diameter.value(),
                        wall_velocity_m_s,
                    },
                    None => WallDrag::Darcy {
                        friction_factor,
                        hydraulic_diameter_m: diameter.value(),
                    },
                }
            }
        };
        // The original source solver owns wall storage and energy exchange;
        // this binding supplies only its instantaneous constitutive coefficients.
        let conductance =
            nusselt_number * conductivity * (4.0 * volume / diameter.value()) / diameter.value();
        if !conductance.is_finite() || conductance < 0.0 {
            return Err(WallExchangeError::Invalid(
                "wall conductance is unrepresentable",
            ));
        }
        match &mut contact.thermal {
            WallThermal::Adiabatic => {}
            WallThermal::Reservoir {
                conductance_w_per_k,
                ..
            }
            | WallThermal::Lumped {
                conductance_w_per_k,
                ..
            }
            | WallThermal::Caloric {
                conductance_w_per_k,
                ..
            } => {
                *conductance_w_per_k = conductance;
            }
        }
        poll(cancelled)?;
        Ok(contact)
    }
}

fn same_geometry(a: f64, b: f64) -> bool {
    a.is_finite() && b.is_finite() && (a - b).abs() <= 128.0 * f64::EPSILON * a.abs().max(b.abs())
}

fn admit(
    quantity: &'static str,
    value: f64,
    minimum: f64,
    maximum: f64,
) -> Result<(), WallExchangeError> {
    if value.is_finite() && (minimum..=maximum).contains(&value) {
        Ok(())
    } else {
        Err(WallExchangeError::TransportDomain {
            quantity,
            value,
            minimum,
            maximum,
        })
    }
}
