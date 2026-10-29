//! Explicit receivers for original native connector dissipation.
use super::*;

/// A declared partition of ONE native spring/dashpot's physical loss.
/// Fractions for each connection must sum to one within floating-point
/// roundoff. They are authored interface data, not a predicted contact or
/// effusivity model. Conservative spring work is never deposited as heat.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrameConnectionHeatDestination {
    pub connection: usize,
    pub cell: CellId,
    pub fraction: f64,
}

/// Actual integrated physical heat delivered to an original finite cell.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrameConnectionHeat {
    pub connection: usize,
    pub cell: CellId,
    pub heat_j: f64,
}

impl<'n, 'a> FrameThermoelasticSystem<'n, 'a> {
    /// Route original non-frame Kelvin-Voigt connector losses into the SAME
    /// finite thermal solve as the bound frames. A receiving cell may be a
    /// frame layer or another material body in the shared network. Its actual
    /// temperature then drives the existing transport and bound-frame source
    /// laws. No heat capacity, source coefficient or extra time step is invented.
    ///
    /// Each mapped connector is fully partitioned. Unmapped connections,
    /// component modal damping, Maxwell memory and acoustic loss keep their
    /// previous external-sink interpretation. Source-frame material ports
    /// already own their layer heat and cannot be mapped a second time.
    /// Frequency-fitted constant-loss connectors require their own thermal
    /// model and are refused here.
    ///
    /// Configuration is permitted once, at a jointly accepted clock boundary.
    /// It changes only the destination of FUTURE losses: no past heat, state,
    /// clock or stored energy is altered. Existing native/thermal work and
    /// conservation gates remain active. `maximum_routes` bounds sorting and
    /// per-step receipt storage; duplicate connection/cell rows refuse.
    pub fn with_connection_heat(
        mut self,
        thermal: &DepositionSimulation<'_>,
        mut destinations: Vec<FrameConnectionHeatDestination>,
        maximum_routes: usize,
        gate: &CancelGate,
    ) -> Result<Self> {
        check(gate)?;
        let expected_time = finite(self.origin_s
            + self.mechanical.samples_rendered() as f64 * self.mechanical.sample_period_s())?;
        if destinations.is_empty() || !self.connection_heat_destinations.is_empty()
            || destinations.len() > maximum_routes
            || thermal.time_s().to_bits() != expected_time.to_bits()
        {
            return Err(invalid("connection heat needs one nonempty bounded map at the jointly accepted clock"));
        }
        // Retain the original caloric/material/transport correspondence, even
        // when the receiver is a passive body outside all bound frames.
        self.prepared.resolve(thermal, || gate.is_requested()).map_err(network)?;
        for (member, state) in self.members.iter().zip(&self.states) {
            check(gate)?;
            if temperatures(&member.region, thermal)?.iter().zip(&state.layers)
                .any(|(&t, layer)| t.to_bits() != layer.temperature_k.to_bits())
            {
                return Err(invalid("connection heat receivers differ from accepted frame temperatures"));
            }
        }
        destinations.sort_unstable_by_key(|r| (r.connection, r.cell));
        let mut first = 0;
        while first < destinations.len() {
            check(gate)?;
            let connection = destinations[first].connection;
            let law = self.mechanical.connection_coefficients(connection)
                .ok_or_else(|| invalid("connection heat names an absent native port"))?;
            let owned = self.owned_connections.partition_point(|range| range.end <= connection);
            if self.owned_connections.get(owned).is_some_and(|range| range.contains(&connection))
                || law.damping_n_s_m <= 0.0
                || self.mechanical.constant_loss_connections().iter()
                    .any(|fit| fit.connection == connection)
            {
                return Err(invalid("connection heat requires a non-frame native Kelvin-Voigt loss port"));
            }
            let mut end = first;
            let mut total = 0.0;
            let mut previous_cell = None;
            while end < destinations.len() && destinations[end].connection == connection {
                check(gate)?;
                let route = destinations[end];
                if !route.fraction.is_finite() || route.fraction <= 0.0 || route.fraction > 1.0
                    || !thermal.cells().contains_key(&route.cell)
                    || previous_cell == Some(route.cell)
                {
                    return Err(invalid("positive finite heat fractions and distinct active receivers required"));
                }
                total = finite(total + route.fraction)?;
                previous_cell = Some(route.cell);
                end += 1;
            }
            // Admission permits summation roundoff, never a missing heat sink.
            // Each step separately compares the actual partition with original
            // connector work in joules under its declared energy tolerance.
            if (total - 1.0).abs() > 8.0 * f64::EPSILON * (end - first) as f64 {
                return Err(invalid("each mapped connection must partition all of its physical heat"));
            }
            first = end;
        }
        check(gate)?;
        self.connection_heat_destinations = destinations;
        Ok(self)
    }

    /// Canonical connection/cell ordering of the admitted future heat map.
    pub fn connection_heat_destinations(&self) -> &[FrameConnectionHeatDestination] {
        &self.connection_heat_destinations
    }
}
