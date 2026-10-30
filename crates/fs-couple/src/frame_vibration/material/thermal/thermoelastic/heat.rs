//! Complete source heat transfer on the existing frame/circuit transaction.
use super::*;
use crate::circuit_enthalpy::{CircuitThermalCoefficients, contacts::ContactThermalCoefficients};
use crate::circuit_enthalpy::contacts::heat::{
    HeatNetworkError, HeatNetworkSample, PreparedHeatNetwork,
};
use fs_conduction::transient::deposition::{
    contact::ContactLinkSample, radiation::RadiationCellSources,
};
use fs_exec::Cx;
use fs_material::thermal::conduction::network::LayerLinkSample;
use std::cell::RefCell;

/// Borrow the original prepared material heat laws and execution authority.
/// Keep the same preparation between steps to retain source pins. No physical
/// state is stored here, and selecting another preparation is an explicit
/// caller intervention, not an automatic source fallback.
#[derive(Clone, Copy)]
pub struct FrameHeatInput<'p, 'n, 'a, 'scope> {
    pub network: &'p PreparedHeatNetwork<'n, 'a>,
    pub cx: &'p Cx<'scope>,
}

/// Original thermal-interface, nonlinear layer and enclosure results at the
/// jointly accepted endpoint. Bulk snapshots remain in the usual `source`
/// report; they are neither cloned nor moved into a second material owner.
#[derive(Debug)]
pub struct FrameHeatTransfer {
    pub thermal_contacts: Vec<ContactLinkSample>,
    pub integrated_conduction: Vec<LayerLinkSample>,
    pub radiation: Option<RadiationCellSources>,
}
impl FrameHeatTransfer {
    pub(super) fn cell_power_w(&self, cell: CellId) -> f64 {
        self.radiation.as_ref()
            .and_then(|radiation| radiation.sources_w().get(&cell)).copied().unwrap_or(0.0)
    }

    pub(super) fn exchange(
        &self, thermal: &DepositionSimulation<'_>, dt: f64, gate: &CancelGate,
    ) -> Result<(f64, f64)> {
        let mut entropy = 0.0;
        if let Some(radiation) = &self.radiation {
            for (&cell, &power) in radiation.sources_w() {
                check(gate)?;
                entropy = finite(entropy + dt * power / thermal.temperature_k(cell)?)?;
            }
        }
        let net = finite(dt * self.radiation.as_ref().map_or(0.0, |r| r.net_heat_w()))?;
        Ok((net, entropy))
    }
}

/// Preserve typed source refusals and all original mechanical/thermal failures.
/// Cancellation remains `Frame(CircuitEnthalpyError::Cancelled)`.
#[derive(Debug)]
pub enum FrameHeatError {
    Heat(HeatNetworkError),
    Frame(FrameThermalError),
}
impl From<FrameThermalError> for FrameHeatError {
    fn from(error: FrameThermalError) -> Self { Self::Frame(error) }
}
impl From<HeatNetworkError> for FrameHeatError {
    fn from(error: HeatNetworkError) -> Self {
        match error {
            HeatNetworkError::Cancelled => Self::Frame(CircuitEnthalpyError::Cancelled),
            other => Self::Heat(other),
        }
    }
}
impl core::fmt::Display for FrameHeatError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "frame and complete material heat transport: {self:?}")
    }
}
impl std::error::Error for FrameHeatError {}

impl FrameThermoelasticSystem<'_, '_> {
    /// Advance the SAME native assembly with complete source-defined thermal
    /// interfaces, source-integrated layer conduction and gray-diffuse radiation.
    /// An electrically bound frame takes `Some(input)`; an unpowered frame takes
    /// `None`. Material C(T)/L(T), electrical contacts and reversible elastic/field
    /// heat continue through their original owners and acceptance boundary.
    ///
    /// The prepared heat network must borrow the original material-cell and
    /// thermal-link slices. It replaces the trial transport evaluation rather
    /// than being added on top of it. Resistor sources, electrical interfaces
    /// and heat partitions must agree with the accepted circuit before motion.
    /// No electrical work may disappear into a thermal-only invocation.
    ///
    /// The report excludes the enclosure's numerical net from external supply
    /// and checks physical storage under the SAME energy tolerance. Individual
    /// cell powers are not corrected. Its signed radiation entropy transfer is
    /// separate from authored-source entropy, conduction and material losses.
    /// Conservation is not an accuracy or experimental-validity certificate.
    ///
    /// Fixed geometry, surface areas/view factors, solid frame laws and the
    /// existing first-order coefficient split remain. Mechanical displacement
    /// does not silently move the enclosure or thermal contacts. Failures from
    /// either execution authority return before the original joint publication;
    /// a Cx cancellation also requests the supplied native mechanical gate.
    #[allow(clippy::too_many_arguments)]
    pub fn step_heat_network(
        &mut self,
        thermal: &mut DepositionSimulation<'_>,
        external_modal_forces: &[f64],
        heat: FrameHeatInput<'_, '_, '_, '_>,
        electrical_input: Option<FrameElectricalInput<'_>>,
        controls: FrameThermalSolve,
        energy_guess_j: &[f64],
        gate: &CancelGate,
    ) -> core::result::Result<FrameThermoelasticReport, FrameHeatError> {
        let trial = Trial { input: heat, failure: RefCell::new(None) };
        trial.checkpoint(gate)?;
        let admission: Result<()> = crate::circuit_enthalpy::admit_driver_unknowns(
            thermal, controls.driver(), energy_guess_j);
        admission?;
        if self.electrical.is_some() != electrical_input.is_some() {
            return Err(invalid("electrically bound frame histories require explicit held electrical inputs").into());
        }
        let material = heat.network.field_storage_material_network();
        if !core::ptr::eq(self.source.cells, material.cells)
            || !core::ptr::eq(self.source.links, material.links) {
            return Err(invalid("frame and complete heat network require the original cell/link slices").into());
        }
        let source_error = |error| {
            if matches!(&error, HeatNetworkError::Cancelled) { gate.request(); }
            FrameHeatError::from(error)
        };
        let mut initial = heat.network.resolve(heat.cx, thermal).map_err(source_error)?;
        let destinations = heat.network.resistor_destinations();
        let routes = heat.network.heat_destinations(thermal).map_err(source_error)?;
        electrical::admit_heat_sources(self.electrical.as_ref(), &mut initial,
            &destinations, &routes, thermal, gate)?;
        drop(initial);
        trial.checkpoint(gate)?;
        // The callback has the original frame error type. Retain a typed source
        // error across that seam; never replace it with a fabricated fallback.
        let result = self.step_inner(thermal, external_modal_forces, electrical_input,
            Some(&trial), controls, energy_guess_j, gate);
        match result {
            Err(CircuitEnthalpyError::Invalid(SOURCE_REFUSAL)) => {
                match trial.failure.into_inner() {
                    Some(error) => Err(FrameHeatError::from(error)),
                    None => Err(invalid("complete heat source refusal lost its original cause").into()),
                }
            }
            other => other.map_err(FrameHeatError::from),
        }
    }
}

// A private marker cannot be mistaken for a source law. The original failure
// is retained until the public boundary, including through numerical stencils.
const SOURCE_REFUSAL: &str = "original complete heat-network source refused frame trial";
pub(super) struct Trial<'p, 'n, 'a, 'scope> {
    input: FrameHeatInput<'p, 'n, 'a, 'scope>,
    failure: RefCell<Option<HeatNetworkError>>,
}
impl Trial<'_, '_, '_, '_> {
    pub(super) fn checkpoint(&self, gate: &CancelGate) -> Result<()> {
        if self.input.cx.checkpoint().is_err() { gate.request(); }
        check(gate)
    }

    pub(super) fn resolve(
        &self, thermal: &DepositionSimulation<'_>, gate: &CancelGate,
    ) -> Result<(ContactThermalCoefficients<MaterialNetworkSample>, FrameHeatTransfer)> {
        self.checkpoint(gate)?;
        self.failure.borrow_mut().take();
        let resolved = match self.input.network.resolve(self.input.cx, thermal) {
            Ok(resolved) => resolved,
            Err(HeatNetworkError::Cancelled) => {
                gate.request();
                return Err(CircuitEnthalpyError::Cancelled);
            }
            Err(error) => {
                *self.failure.borrow_mut() = Some(error);
                return Err(invalid(SOURCE_REFUSAL));
            }
        };
        self.checkpoint(gate)?;
        let ContactThermalCoefficients { bulk, contacts } = resolved;
        let HeatNetworkSample { material, thermal_contacts, integrated_conduction, radiation } = bulk.source;
        Ok((ContactThermalCoefficients {
            contacts,
            bulk: CircuitThermalCoefficients {
                resistors: bulk.resistors, sources_w: bulk.sources_w,
                conductances_w_k: bulk.conductances_w_k, source: material,
            },
        }, FrameHeatTransfer { thermal_contacts, integrated_conduction, radiation }))
    }
}
