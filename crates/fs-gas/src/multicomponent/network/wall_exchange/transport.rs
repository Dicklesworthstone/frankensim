//! Material-dependent developed pipe-wall losses in the existing wall transaction.

use super::{NetworkWallOperator, WallContact, WallDrag, WallExchangeError, WallThermal, poll};
use crate::multicomponent::dissipation::TransportLaw;
use crate::multicomponent::dissipation::evolve::fractions;
use crate::multicomponent::{Cell, MAX_SPECIES, Mesh1d, MixtureGas, allocated, positive};
use fs_qty::{Length, Pressure, Temperature};

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
        let area = core::f64::consts::FRAC_PI_4 * self.diameter.value() * self.diameter.value();
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
}

pub(super) struct CircularWallTransport {
    law: Box<dyn TransportLaw + Send + Sync>,
    branches: Vec<Option<CircularLaminarWall>>,
    binary: Vec<f64>,
}

impl NetworkWallOperator {
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
        mut self,
        law: L,
        branches: Vec<Option<CircularLaminarWall>>,
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

    /// Original developed-wall geometry and applicability, in branch order.
    pub fn circular_laminar_walls(&self) -> Option<&[Option<CircularLaminarWall>]> {
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
        let diameter = model.diameter.value();
        let reynolds = primitive.density * primitive.velocity.abs() * diameter / viscosity;
        let mach = primitive.velocity.abs() / primitive.sound_speed;
        let [minimum, maximum] = model.reynolds_number;
        admit("Reynolds number", reynolds, minimum, maximum)?;
        admit("Mach number", mach, 0.0, model.maximum_mach_number)?;
        contact.drag = WallDrag::circular_laminar(
            coefficients.dynamic_viscosity,
            primitive.density,
            model.diameter,
        )?;
        // The original source solver owns wall storage and energy exchange;
        // this binding supplies only its instantaneous constitutive coefficients.
        let conductance =
            model.nusselt_number * conductivity * (4.0 * volume / diameter) / diameter;
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
