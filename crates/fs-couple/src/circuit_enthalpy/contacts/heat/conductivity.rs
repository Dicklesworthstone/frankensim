//! Native nonlinear electrical response with the original complete heat owner.
use super::*;
use super::super::material::ContactConductionResponse;
use crate::circuit_enthalpy::material::CellConductionResponse;
use crate::circuit_enthalpy::motion::support::{self, SupportThermalBinding,
    SupportThermalCandidate, SupportThermalError};
use crate::circuit_enthalpy::StepCoefficients;
use crate::material_circuit::conductivity::ConductivityIteration;
use crate::material_circuit::motion::{MaterialFieldMotion,
    conductivity::FieldConductivityStep};

impl<'a> MaterialHeatNetwork<'a> {
    /// Select the original bulk, electrical-interface, thermal-interface and
    /// emissivity claims for nonlinear current/field response. The bulk and
    /// contact descriptors declare operating coordinates, not replacement cards.
    /// Directional Fourier bindings retain their original principal sources.
    /// Auxiliary thermal contacts and emissivity must be electrically independent
    /// in both value and applicability; their source owners do not solve fields.
    pub fn prepare_conductivity<'n>(
        &'n self, cx: &Cx<'_>, thermal: &DepositionSimulation<'_>,
        responses: &[CellConductionResponse<'_>], contacts: &[ContactConductionResponse<'_>],
        directions: &[DirectionalConductionLink],
    ) -> Result<PreparedHeatNetwork<'n, 'a>, HeatNetworkError> {
        check(cx)?;
        let electrical = self.electrical.prepare_conductivity(thermal, responses, contacts,
            directions, || cx.checkpoint().is_err())?;
        let initial = electrical.resolve_conductivity_seeds(thermal, responses, contacts,
            || cx.checkpoint().is_err())?;
        let prepared = self.prepare_from_electrical(cx, thermal, electrical, initial,
            self.thermal_contacts, self.radiation.map_or(&[], |radiation| radiation.cells), None)?;
        prepared.admit_conductivity_heat(cx, thermal, responses, contacts)?;
        Ok(prepared)
    }
}

impl PreparedHeatNetwork<'_, '_> {
    fn admit_conductivity_heat(
        &self, cx: &Cx<'_>, thermal: &DepositionSimulation<'_>,
        responses: &[CellConductionResponse<'_>], contacts: &[ContactConductionResponse<'_>],
    ) -> Result<(), HeatNetworkError> {
        check(cx)?;
        let invalid = HeatNetworkError::Invalid;
        if self.integrated_conduction.is_some()
            || self.network.electrical.bulk.links.len() != thermal.links().len() {
            return Err(invalid("nonlinear electrical heat requires original Fourier links without integrated layers or additional link owners"));
        }
        let mut axes = Vec::new();
        for response in responses {
            check(cx)?;
            let resistor = self.network.electrical.bulk.resistors.iter()
                .find(|resistor| resistor.destination.branch == response.branch)
                .ok_or(invalid("electrical response has no original bulk heat destination"))?;
            axes.extend(response.axes.iter().map(|axis| (resistor.destination.cell, axis.name)));
        }
        for response in contacts {
            check(cx)?;
            let link = thermal.links().get(&response.id)
                .ok_or(invalid("electrical response has no original interface heat destination"))?;
            for cell in [link.a, link.b] {
                axes.extend(response.axes.iter().map(|axis| (cell, axis.name)));
            }
        }
        self.admit_field_independent_boundaries(cx, &axes)
    }

    /// Complete pinned heat coefficients at the current thermal endpoint and
    /// authored electrical seeds. Use the returned resistors/contacts to build
    /// the original circuit; these seed fields are not accepted operating points.
    pub fn resolve_conductivity_seeds(
        &self, cx: &Cx<'_>, thermal: &DepositionSimulation<'_>,
        responses: &[CellConductionResponse<'_>], contacts: &[ContactConductionResponse<'_>],
    ) -> Result<ContactThermalCoefficients<HeatNetworkSample>, HeatNetworkError> {
        self.admit_conductivity_heat(cx, thermal, responses, contacts)?;
        let electrical = self.electrical.resolve_conductivity_seeds(thermal, responses, contacts,
            || cx.checkpoint().is_err())?;
        self.complete_heat_coefficients(cx, thermal, electrical, &[], &[])
    }

    fn refresh_conductivity_heat(
        &self, cx: &Cx<'_>, thermal: &DepositionSimulation<'_>,
        electrical: &FieldConductivityStep, coefficients: StepCoefficients<HeatNetworkSample>,
    ) -> Result<StepCoefficients<HeatNetworkSample>, HeatNetworkError> {
        check(cx)?;
        let HeatNetworkSample { material, thermal_contacts: _, integrated_conduction, radiation } = coefficients.payload;
        let mut bulk = self.electrical.refresh_conductivity_transport(thermal, &electrical.evaluated,
            StepCoefficients { sources_w: coefficients.sources_w,
                conductances_w_k: coefficients.conductances_w_k, payload: material },
            || cx.checkpoint().is_err())?;
        // Refresh from bulk-only solved-field G. The interface series resistance
        // is applied exactly once; reusing the seed contact G would lose E/J -> k.
        let thermal_contacts = apply_material_contacts(thermal, &self.thermal_contacts,
            &mut bulk.conductances_w_k, || cx.checkpoint().is_err())?;
        check(cx)?;
        // Radiation was already evaluated at this same trial T. Its admitted
        // claims consume no electrical coordinate, so retain its result/power
        // once instead of adding the closed-enclosure exchange a second time.
        Ok(StepCoefficients { sources_w: bulk.sources_w, conductances_w_k: bulk.conductances_w_k,
            payload: HeatNetworkSample { material: bulk.payload, thermal_contacts,
                integrated_conduction, radiation } })
    }

    /// Close nonlinear bulk/contact current, thermal transport and live-support
    /// motion through their original shared residual. Thermal interfaces are
    /// composed with solved-field Fourier G; enclosure exchange remains internal.
    /// All original source pins, interface directions and A/B heat routes survive.
    /// Native electrical iterations consume the same circuit work budget.
    ///
    /// Geometry, inertia and field storage retain the support owner's restrictions.
    /// Field-dependent support/caloric/auxiliary heat claims, integrated layers
    /// and externally owned additional links require their own coupled owners
    /// and refuse. A failed candidate publishes no physical history.
    #[allow(clippy::too_many_arguments)]
    pub fn field_motion_support_conductivity_candidate<'m>(
        &self, cx: &Cx<'_>, thermal: &DepositionSimulation<'m>, target_s: f64,
        motion: &MaterialFieldMotion, ports: &[f64], external_force_n: f64,
        responses: &[CellConductionResponse<'_>], contacts: &[ContactConductionResponse<'_>],
        binding: &SupportThermalBinding<'_>, controls: CircuitEnthalpySolve,
        conductivity: impl Into<ConductivityIteration>, guess_j: &[f64],
    ) -> Result<SupportThermalCandidate<'m, HeatNetworkSample, FieldConductivityStep>,
        CircuitEnthalpyError<SupportThermalError>> {
        let source_error = |error| match error {
            HeatNetworkError::Cancelled => CircuitEnthalpyError::Cancelled,
            other => CircuitEnthalpyError::Source(SupportThermalError::Heat(other)),
        };
        self.admit_conductivity_heat(cx, thermal, responses, contacts).map_err(source_error)?;
        let (bulk_sources, contact_sources) = self.electrical.conductivity_sources(responses, contacts,
            || cx.checkpoint().is_err()).map_err(HeatNetworkError::from).map_err(source_error)?;
        let routes = self.heat_destinations(thermal).map_err(source_error)?;
        let result = support::step_conductivity_candidate(self.network.electrical.bulk,
            thermal, target_s, motion, ports, external_force_n, &bulk_sources, Some(&contact_sources),
            &routes, binding, controls, conductivity, guess_j,
            |state| self.resolve_conductivity_seeds(cx, state, responses, contacts),
            |state, electrical, transport| self.refresh_conductivity_heat(cx, state, electrical, transport),
            HeatNetworkSample::radiation_net_power_w, || cx.checkpoint().is_err());
        match result {
            Err(CircuitEnthalpyError::Source(SupportThermalError::Heat(HeatNetworkError::Cancelled))) =>
                Err(CircuitEnthalpyError::Cancelled),
            other => other,
        }
    }
}
