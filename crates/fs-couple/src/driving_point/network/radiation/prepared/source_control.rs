//! Original generator coordinates through loaded openings into BEM control.

use super::*;
use fs_bem::observation_transfers::{
    HarmonicControlLimits, HarmonicSourceMapLimits, HarmonicVolumeQuadrature,
    PreparedHarmonicEnergyControl,
};
use fs_duct::network::NetworkSourceQuantity;

/// Original source controls and their complete loaded physical exterior map.
///
/// Pressure controls remain open-circuit pressures [Pa], before the full
/// generator impedance. Prescribed flow controls remain graph inflows [m³/s].
/// Opening radiation and all mutual source/load impedances have already been
/// solved by the borrowed original owner. No independent outlet actuation,
/// new radiation solve, radiation-power debit or optimizer state is created.
///
/// The retained composed observation has canonical source labels, order and
/// SI units. Its unchanged BEM energy/adjoint/HVP owners therefore differentiate
/// actual upstream generator coefficients. Each source scale or magnitude
/// constraint belongs to that same coordinate, including mixed-unit banks.
pub struct PreparedRadiatingSourceControl<'p, 'a> {
    owner: &'p PreparedRadiatingNetwork<'a>,
    observation: PreparedHarmonicObservation,
}

impl<'a> PreparedRadiatingNetwork<'a> {
    /// Cold, bounded preparation of physical exterior control from the
    /// original complete source bank. Reuses the retained source-to-opening
    /// map, with Q_out = -U_graph, and the exact original exterior field basis.
    ///
    /// All source coordinates are retained, including disturbance-only or
    /// acoustically dark inputs. Silent columns do not bypass admission.
    /// The new observation costs four products per source/opening/receiver;
    /// subsequent value/gradient and HVP calls allocate nothing. The BEM
    /// limits include the temporary source identities and supplied map;
    /// the fixed 32-by-32 stack scratch is an additional bounded temporary.
    pub fn source_control(
        &self,
        limits: HarmonicSourceMapLimits,
        gate: &CancelGate,
    ) -> Result<PreparedRadiatingSourceControl<'_, 'a>, ExteriorNetworkError> {
        poll(gate)?;
        let coordinates = self.source_coordinates();
        let (n, openings) = (coordinates.len(), self.opening_indices.len());
        if n == 0 || n > MAX_NETWORK_PORTS || openings == 0 || openings > MAX_NETWORK_PORTS {
            return Err(bad(
                "loaded source control requires its complete finite source/opening bank",
            ));
        }
        if n > limits.maximum_sources || self.observer.sources().len() != openings {
            return Err(bad(
                "loaded source control basis exceeds or differs from its original owner",
            ));
        }
        // Check identity storage before copying original labels. The complete
        // composition allowance is subsequently checked by the BEM owner.
        let names = coordinates
            .iter()
            .try_fold(0usize, |bytes, coordinate| {
                bytes.checked_add(coordinate.label.len())
            })
            .and_then(|bytes| bytes.checked_mul(2))
            .and_then(|bytes| bytes.checked_add(2 * n * size_of::<SourceId>()))
            .ok_or_else(|| bad("loaded source control identity extent overflow"))?;
        if names
            > limits
                .maximum_storage_bytes
                .min(HarmonicSourceMapLimits::MAXIMUM_STORAGE_BYTES)
        {
            return Err(bad(
                "loaded source control identity storage exceeds its allowance",
            ));
        }
        let mut sources = Vec::new();
        sources
            .try_reserve_exact(n)
            .map_err(|_| bad("loaded source control identity allocation"))?;
        for coordinate in coordinates {
            poll(gate)?;
            let mut label = String::new();
            label
                .try_reserve_exact(coordinate.label.len())
                .map_err(|_| bad("loaded source control label allocation"))?;
            label.push_str(&coordinate.label);
            sources.push(match coordinate.quantity {
                NetworkSourceQuantity::OpenCircuitPressure => SourceId::OpenCircuitPressure(label),
                NetworkSourceQuantity::VolumeFlow => SourceId::PrescribedVolumeFlow(label),
            });
        }
        let mut map = [C64::ZERO; FACTOR_ENTRIES];
        for (opening, &port) in self.opening_indices.iter().enumerate() {
            poll(gate)?;
            if !matches!(&self.observer.sources()[opening],
                SourceId::ApertureVolumeFlow(label) if label == &self.ports()[port])
            {
                return Err(bad(
                    "loaded source control opening identity differs from its exterior owner",
                ));
            }
            for source in 0..n {
                map[source * openings + opening] = -self.network.flow_transfer()[port * n + source];
            }
        }
        let observation =
            self.observer
                .compose_sources(&sources, &map[..n * openings], limits, || {
                    gate.is_requested()
                })?;
        poll(gate)?;
        Ok(PreparedRadiatingSourceControl {
            owner: self,
            observation,
        })
    }
}

impl<'p, 'a> PreparedRadiatingSourceControl<'p, 'a> {
    /// Exact original loaded owner, including physical work and source laws.
    /// Evaluating controls never adds another debit to that owner's ledger.
    pub fn owner(&self) -> &'p PreparedRadiatingNetwork<'a> {
        self.owner
    }
    pub fn source_coordinates(&self) -> &[NetworkSourceCoordinate] {
        self.owner.source_coordinates()
    }
    pub fn default_excitations(&self) -> &[C64] {
        self.owner.default_excitations()
    }
    /// Original BEM field operations now in canonical upstream source units.
    /// `source_pullback` and all receiver/microphone observers use this same
    /// composed field; none implicitly treats generator pressure as flow.
    pub fn observation(&self) -> &PreparedHarmonicObservation {
        &self.observation
    }
    /// Construct the unchanged BEM energy metric, gradient and HVP owner.
    ///
    /// Quadrature uses the original physical receiver indices and m³ weights.
    /// Targets have this exact frequency, frame and phase reference. This is
    /// exterior residual acoustic energy [J], not generator dissipation or an
    /// actuator effort penalty. Its complete source-major Cartesian complex
    /// gradient is suitable for the original ASCENT source-control adapter.
    pub fn energy_control<'s>(
        &'s self,
        quadrature: &'s [HarmonicVolumeQuadrature],
        targets: Option<&'s [HarmonicField]>,
        limits: HarmonicControlLimits,
        gate: &CancelGate,
    ) -> Result<PreparedHarmonicEnergyControl<'s>, ExteriorNetworkError> {
        Ok(PreparedHarmonicEnergyControl::new(
            &self.observation,
            quadrature,
            targets,
            limits,
            || gate.is_requested(),
        )?)
    }
}
