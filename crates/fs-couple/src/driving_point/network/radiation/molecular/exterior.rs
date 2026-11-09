//! The original relaxing exterior supplies both the load and receiver bank.

use super::*;
use fs_bem::helmholtz::laplace::complex::ComplexOptions;
use fs_bem::observation_transfers::{
    MolecularBemObservation, MolecularObservationFormulation, MolecularResponseLimits,
    PreparedMolecularResponse,
};
use fs_bem::observation_transfers::{ReceiverSpec, SourceId};
use fs_material::relaxation::RelaxationState;

/// Explicit exterior material and original numerical owners for one held band.
/// The full gas state must match the network's exterior record; the borrowed
/// reservoir law is never inferred from density, sound speed or attenuation.
///
/// This models homogeneous molecular compression relaxation with inviscid
/// momentum. Exterior shear/thermal diffusion and distributed wall radiation
/// are outside its scope. Boundary work includes exterior absorption and
/// outgoing energy; it is not a lossless far-field radiated-power measure.
#[derive(Clone, Copy, Debug)]
pub struct MolecularExteriorMaterial<'a> {
    pub source: &'a RelaxationState,
    pub band_rad_s: [f64; 2],
    pub formulation: MolecularObservationFormulation,
    pub assembly: ComplexOptions,
    pub response: MolecularResponseLimits,
}

impl MolecularExteriorMaterial<'_> {
    pub(in super::super) fn same_source(&self, other: &Self) -> bool {
        core::ptr::eq(self.source, other.source)
            && self.band_rad_s == other.band_rad_s
            && self.formulation == other.formulation
    }

    pub(in super::super) fn prepare(
        &self,
        model: &RadiatingNetwork,
        omega: f64,
        gate: &CancelGate,
    ) -> Result<(CoupledTermination, PreparedMolecularResponse), ExteriorNetworkError> {
        poll(gate)?;
        if self.source.gas() != &model.gas {
            return Err(bad("molecular exterior must retain the original complete gas record"));
        }
        let source_ids: Vec<_> = model.patches.iter()
            .map(|patch| SourceId::ApertureVolumeFlow(patch.port.clone()))
            .collect();
        let traces: Vec<Vec<f64>> = model.traces.iter()
            .map(|trace| trace.iter().map(|value| value.re).collect())
            .collect();
        let receivers: Vec<_> = model.receivers.iter().enumerate()
            .map(|(i, &position_m)| ReceiverSpec {
                id: format!("network-receiver-{i}"),
                position_m,
            })
            .collect();
        let owner = MolecularBemObservation::new_with_formulation(
            &model.surface,
            (*self.source).clone(),
            &source_ids,
            &traces,
            &receivers,
            model.clearance_m,
            self.band_rad_s,
            self.formulation,
            self.assembly,
        )?;
        // At most 32 openings fit the original 256-source batch. One native
        // complex-K assembly/LU supplies both reaction and physical p/u fields.
        let response = owner.prepare_response(omega, self.response, || gate.is_requested())?;
        let boundary = CoupledTermination {
            omega,
            labels: model.driving.coupled_labels().to_vec(),
            z: response.acoustic_impedance()?.to_vec(),
            source_id: format!("fs-bem:molecular-network-openings:{}", owner.operator_hash()),
        };
        poll(gate)?;
        Ok((boundary, response))
    }
}
