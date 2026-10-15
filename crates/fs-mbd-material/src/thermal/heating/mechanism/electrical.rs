//! The existing electrical/enthalpy driver publishes source-bound mechanics too.
use fs_exec::CancelGate;
use fs_lbm::conjugate::enthalpy::{JouleEnthalpyError, JouleEnthalpyInput};
use fs_lbm::conjugate::enthalpy::transport::{EnthalpyTransport, EnthalpyTransportSample,
    JouleTransportReport, advance_joule_enthalpy_with_material_transport_observed};
use fs_thermal::enthalpy::{ContactPatch, DepositedCell, StepConfig};
use fs_thermal::observers::ObserverCouplingError;
use crate::thermal::heating::{HeatedBody, HeatingError};
use super::HeatedArticulated;

impl HeatedArticulated<'_> {
    /// Evolve current, material heat transport and ALL link configurations together.
    ///
    /// This delegates to the original spatial current/enthalpy driver. Its
    /// geometry map, current/contact power deposition, actual thermal CFL,
    /// endpoint support checks and work limits remain authoritative. Mechanical
    /// material admission runs on the SAME thermal substeps. Nothing is returned
    /// until final electrical, thermal and mechanical acceptance all succeed.
    /// No copied electrical solver or intermediate thermal publication is used.
    ///
    /// `resolve` should consume the original per-cell material card and current
    /// state; `regions()` exposes those bindings. Its original concrete payload
    /// is retained in the report. `resolve_link` consumes the SAME adjacent
    /// samples, with explicit thermal geometry/contact assumptions. Unmapped
    /// thermal cells may exchange heat but receive no electrical source.
    /// The thermal audit ALREADY includes Joule heat; do not add it again.
    ///
    /// The electrical and thermal grids stay at their declared reference
    /// geometry while solid free-expansion mechanics is configured. This is
    /// explicit stationary electrothermal coupling, not moving-grid current,
    /// induction, constrained deformation or live mechanical energy transfer.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn advance_joule<F, I>(
        &self, target_s: f64, config: StepConfig, input: JouleEnthalpyInput<'_>,
        link_ids: &[u64], gate: &CancelGate,
        resolve: impl Fn(usize, &DepositedCell, f64) -> Result<EnthalpyTransport<F>, JouleEnthalpyError>,
        resolve_link: impl Fn(&ContactPatch, [&EnthalpyTransportSample<F>; 2]) -> Result<(f64, I), JouleEnthalpyError>,
    ) -> Result<(Self, JouleTransportReport<F, I>), ObserverCouplingError<JouleEnthalpyError, HeatingError>> {
        let (thermal, mechanics, report) = advance_joule_enthalpy_with_material_transport_observed(
            &self.thermal, target_s, config, input, link_ids, gate, &self.mechanics,
            resolve, resolve_link,
        )?;
        Ok((Self { thermal, mechanics }, report))
    }
}

impl HeatedBody<'_> {
    /// The same shared spatial electrical/thermal advance for one material body.
    ///
    /// Geometry matching, source callbacks, concrete material receipts, exact-once
    /// Joule deposition, endpoint refusal and fixed-grid limitations are those of
    /// `HeatedArticulated::advance_joule`. The body and its complete thermal
    /// history are published together only after the original driver accepts.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn advance_joule<F, I>(
        &self, target_s: f64, config: StepConfig, input: JouleEnthalpyInput<'_>,
        link_ids: &[u64], gate: &CancelGate,
        resolve: impl Fn(usize, &DepositedCell, f64) -> Result<EnthalpyTransport<F>, JouleEnthalpyError>,
        resolve_link: impl Fn(&ContactPatch, [&EnthalpyTransportSample<F>; 2]) -> Result<(f64, I), JouleEnthalpyError>,
    ) -> Result<(Self, JouleTransportReport<F, I>), ObserverCouplingError<JouleEnthalpyError, HeatingError>> {
        let (thermal, mechanics, report) = advance_joule_enthalpy_with_material_transport_observed(
            &self.thermal, target_s, config, input, link_ids, gate, &self.mechanics,
            resolve, resolve_link,
        )?;
        Ok((Self { thermal, mechanics }, report))
    }
}
