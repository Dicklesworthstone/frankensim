//! Thermoelectric transport and piezoelectric motion in one native circuit.
//!
//! The original Kelvin owner supplies source-integrated Seebeck voltage and
//! reciprocal junction/interior heat. The original piezoelectric descriptor
//! evolves charge, rod motion and every other circuit storage together. Their
//! actual midpoint current and Joule loss enter the same finite enthalpy solve.

use super::*;
use crate::circuit_enthalpy::piezoelectric::PiezoelectricEnthalpyReport;
use crate::material_circuit::piezoelectric::{PiezoelectricCircuit, PiezoelectricCircuitCandidate};

impl ThermoelectricDriver for PiezoelectricCircuitCandidate {
    fn native(&self) -> &StepRecord {
        &self.step().native
    }
}

/// Original physical transfers and the external supply of the complete system.
pub struct ThermoelectricPiezoelectricReport {
    /// Original rod/circuit, finite enthalpy and source-integrated Kelvin accounts.
    /// Piezoelectric strain/field storage and resistor Joule heat appear once.
    pub coupled: PiezoelectricEnthalpyReport<ThermoelectricSample>,
    /// Authored electrical work, excluding internal Seebeck source ports [J].
    pub external_electrical_work_j: f64,
    /// Held external force times the actual rod extension change [J].
    pub external_mechanical_work_j: f64,
    /// Authored thermal source/reservoir heat, excluding Joule and Kelvin heat [J].
    pub external_heat_j: f64,
    /// Circuit + rod + thermal storage change minus all external supply [J].
    pub total_energy_residual_j: f64,
}

/// One admitted circuit, mechanical and thermal endpoint; inputs stay unchanged.
pub struct ThermoelectricPiezoelectricCandidate<'m> {
    thermal: DepositionSimulation<'m>,
    model: PiezoelectricCircuit,
    report: ThermoelectricPiezoelectricReport,
}
impl<'m> ThermoelectricPiezoelectricCandidate<'m> {
    /// Proposed finite junction, leg and passive-cell caloric histories.
    pub const fn thermal(&self) -> &DepositionSimulation<'m> {
        &self.thermal
    }
    /// Proposed original circuit charge/flux and reciprocal rod history.
    pub const fn model(&self) -> &PiezoelectricCircuit {
        &self.model
    }
    /// Original sources, transfers and independent external-work balance.
    pub const fn report(&self) -> &ThermoelectricPiezoelectricReport {
        &self.report
    }
    /// Publish all histories together without replaying any material law.
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'m>,
        PiezoelectricCircuit,
        ThermoelectricPiezoelectricReport,
    ) {
        (self.thermal, self.model, self.report)
    }
}

/// Advance a finite thermocouple/Peltier network and its piezoelectric load.
///
/// Every thermal trial resolves the original pinned electrical/conductive and
/// Seebeck sources, then restarts the original piezoelectric circuit from the
/// same accepted charge, flux, extension and momentum. The ACTUAL Kirchhoff
/// source current sets Kelvin junction and bulk heat; native resistor Joule
/// heat reaches each original finite cell once. Mechanical strain and voltage
/// exchange reversible energy inside the original coupled piezoelectric storage.
///
/// `ports` uses the original circuit order, with zero placeholders at internal
/// Seebeck ports. The force is positive toward rod extension and is appended
/// only by the original piezoelectric owner. The thermal clock defines the
/// common duration. One candidate passes native rod/circuit, cell thermal and
/// complete external-work balances before any accepted history can advance.
/// Refusal or cancellation publishes no partial endpoint.
///
/// The existing solid thermoelectric legs and uniform small-signal longitudinal
/// ceramic scopes apply. Ceramic source context, geometry and inertia stay
/// fixed; thermal cells evolve the stationary transport network. No ceramic
/// thermal strain, temperature-dependent piezoelectric coefficients, dielectric
/// loss, depoling, electrical contacts or magnetic transport is inferred.
/// Midpoint electrical/mechanical motion with endpoint thermal coefficients is
/// first order in general; balance closure does not certify time accuracy.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate<'m>(
    material: &MaterialNetwork<'_>,
    thermal: &DepositionSimulation<'m>,
    target_s: f64,
    model: &PiezoelectricCircuit,
    ports: &[f64],
    external_force_n: f64,
    legs: &[ThermoelectricLeg],
    electrical_bases: &[(usize, ElectricalPropertyBasis)],
    controls: CircuitEnthalpySolve,
    guess_j: &[f64],
    mut cancelled: impl FnMut() -> bool,
) -> Result<ThermoelectricPiezoelectricCandidate<'m>, Error> {
    let (thermal, model, coupled) = coupled_candidate(
        material,
        thermal,
        target_s,
        model.circuit(),
        ports,
        legs,
        electrical_bases,
        controls,
        guess_j,
        |sources, inputs, dt, cancelled| {
            model.step_candidate(
                sources,
                inputs,
                external_force_n,
                dt,
                controls.circuit,
                cancelled,
            )
        },
        &mut cancelled,
    )?;
    let field = &coupled.electrical;
    let external_electrical_work_j =
        finite(field.electrical_supply_j - coupled.source.reversible_electrical_work_j)?;
    let external_mechanical_work_j = finite(field.mechanical_supply_j)?;
    let external_heat_j = finite(coupled.non_electrical_heat_j - coupled.source.reversible_heat_j)?;
    let total_energy_residual_j = finite(
        field.stored_after_j - field.stored_before_j + coupled.thermal.stored_energy_change_j
            - external_electrical_work_j
            - external_mechanical_work_j
            - external_heat_j,
    )?;
    if total_energy_residual_j.abs() > controls.total_energy_tolerance_j {
        return Err(invalid(
            "thermoelectric, piezoelectric and thermal storage do not close external supply",
        ));
    }
    if cancelled() {
        return Err(CircuitEnthalpyError::Cancelled);
    }
    Ok(ThermoelectricPiezoelectricCandidate {
        thermal,
        model,
        report: ThermoelectricPiezoelectricReport {
            coupled,
            external_electrical_work_j,
            external_mechanical_work_j,
            external_heat_j,
            total_energy_residual_j,
        },
    })
}
