//! Material-body diffraction on the original DG terminal-flow interface.
//!
//! The cold BEM source retains its actual finite apertures, passive surfaces
//! and exterior gas. This adapter binds those source coordinates to the same
//! simultaneous terminal solve used by the other admitted MIMO loads. It does
//! not change the source fit, divide by aperture area again, or insert a sample
//! of feedback delay. Finite side branches use their outer outlet flow.

use super::{AcousticRealizeError, DgBoundary, DgBranch, DgDuct, DgPortLoad, PortId, invalid};
use crate::bernoulli_aperture::radiation::environment::multiport::bake::boundary::{
    BoundaryMimoLoad, BoundaryMimoRecord,
};
use fs_bem::observation_transfers::SourceId;
use fs_duct::profile::gas_field::GasFieldDuctError;
use fs_exec::CancelGate;

mod observation;
pub use observation::DgBoundaryRadiationObserver;

impl DgPortLoad {
    /// Connect a quiescent material-diffraction realization to every original
    /// non-inlet DG port. Each `(port, source)` pair must follow the ORIGINAL
    /// baked matrix order; a name mismatch refuses instead of permuting a
    /// matrix or interpreting a structural coordinate as volume flow.
    ///
    /// The physical step comes from the admitted source. Every axial aperture
    /// retains the actual outlet area; a side port is an explicitly authored
    /// flush aperture. Exterior and local interior gases may differ in
    /// composition and temperature, but must share static equilibrium pressure.
    /// The original same-interval work-conjugate coupling retains all mutual
    /// terms, native histories and transaction semantics. Source fit errors,
    /// `GridOnly` sampled passivity and unresolved BEM spatial error remain as
    /// reported; this is not a closed-loop accuracy certificate.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_boundary_radiation(
        duct: DgDuct,
        exponential_tolerance: f64,
        maximum_dofs: usize,
        bindings: Vec<(PortId, SourceId)>,
        source: BoundaryMimoLoad,
        gate: &CancelGate,
    ) -> Result<Self, AcousticRealizeError> {
        Self::new_with_boundary_radiation_branches(
            duct,
            exponential_tolerance,
            maximum_dofs,
            bindings,
            source,
            Vec::new(),
            gate,
        )
    }

    /// The same complete exterior with finite DG side ducts. A binding's
    /// original bore-side ID denotes its branch OUTLET after attachment;
    /// internal junction flow is never substituted for radiated volume flow.
    /// Sealed branches retain their interior storage but have no exterior row.
    /// Source support area is checked against each actual bore/branch outlet.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_boundary_radiation_branches(
        duct: DgDuct,
        exponential_tolerance: f64,
        maximum_dofs: usize,
        bindings: Vec<(PortId, SourceId)>,
        source: BoundaryMimoLoad,
        branches: Vec<DgBranch>,
        gate: &CancelGate,
    ) -> Result<Self, AcousticRealizeError> {
        poll(gate)?;
        let dt = source.load().time_step_s();
        super::branch::admit_sources(&duct, &branches, dt, maximum_dofs)?;
        admit_bindings(&duct, &branches, &bindings, source.record(), gate)?;
        let ports = bindings.into_iter().map(|(port, _)| port).collect();
        let (terminal, record) = source.into_parts();
        let mut load = Self::new_with_terminal_branches(
            duct,
            dt,
            exponential_tolerance,
            maximum_dofs,
            ports,
            terminal,
            branches,
        )?;
        load.boundary_radiation_record = Some(record);
        poll(gate)?;
        Ok(load)
    }

    /// Original unbaffled source geometry/material identity, exact aperture
    /// ordering, numerical fit diagnostics and exterior gas. The corresponding
    /// DG port at each index is `terminal_ports()[index]`. A generic authored
    /// terminal or a baffled-piston load does not manufacture this record.
    #[must_use]
    pub fn boundary_radiation_record(&self) -> Option<&BoundaryMimoRecord> {
        self.boundary_radiation_record.as_ref()
    }
}

fn poll(gate: &CancelGate) -> Result<(), AcousticRealizeError> {
    if gate.is_requested() {
        Err(AcousticRealizeError::Cancelled)
    } else {
        Ok(())
    }
}

fn admit_bindings(
    duct: &DgDuct,
    branches: &[DgBranch],
    bindings: &[(PortId, SourceId)],
    source: &BoundaryMimoRecord,
    gate: &CancelGate,
) -> Result<(), AcousticRealizeError> {
    let sealed = branches
        .iter()
        .filter(|branch| branch.duct.outlet() == DgBoundary::Closed)
        .count();
    let expected = duct.side_ports().len()
        + usize::from(matches!(duct.outlet(), DgBoundary::Port(_)))
        - sealed;
    if !matches!(duct.inlet(), DgBoundary::Port(_))
        || bindings.is_empty()
        || bindings.len() != expected
        || bindings.len() != source.sources().len()
        || bindings.len() != source.aperture_support_areas_m2().len()
    {
        return Err(invalid(
            "boundary radiation must bind every unsealed non-inlet DG port to the complete original aperture basis",
        ));
    }
    for (index, (port, id)) in bindings.iter().enumerate() {
        poll(gate)?;
        if id != &source.sources()[index]
            || !matches!(id, SourceId::ApertureVolumeFlow(_))
            || duct.inlet() == DgBoundary::Port(*port)
            || bindings[..index].iter().any(|(other, _)| other == port)
            || branches
                .iter()
                .any(|branch| branch.port == *port && branch.duct.outlet() == DgBoundary::Closed)
        {
            return Err(invalid(
                "boundary radiation requires distinct original ports in the exact baked aperture order",
            ));
        }
        let endpoint = branches
            .iter()
            .find(|branch| branch.port == *port)
            .map(|branch| &branch.duct)
            .or_else(|| (duct.outlet() == DgBoundary::Port(*port)).then_some(duct));
        let gas = if let Some(endpoint) = endpoint {
            let area = source.aperture_support_areas_m2()[index];
            let actual = endpoint.outlet_area_m2();
            if !area.is_finite()
                || area <= 0.0
                || (area - actual).abs() > 128.0 * f64::EPSILON * actual
            {
                return Err(invalid(
                    "boundary aperture support must retain the actual DG bore or branch outlet area",
                ));
            }
            endpoint.outlet_gas()
        } else {
            let local = duct.side_port_source(*port).map_err(|error| match error {
                GasFieldDuctError::Duct(source) => AcousticRealizeError::Duct(source),
                GasFieldDuctError::Material(source) => {
                    AcousticRealizeError::Ambient(source.to_string())
                }
            })?;
            if local.gas.pressure.to_bits() != source.gas().pressure.to_bits() {
                return Err(invalid(
                    "boundary radiation and local DG aperture gas require equal equilibrium pressures",
                ));
            }
            continue;
        };
        if gas.pressure.to_bits() != source.gas().pressure.to_bits() {
            return Err(invalid(
                "boundary radiation and local DG aperture gas require equal equilibrium pressures",
            ));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests;
