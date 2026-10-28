//! Source-retaining media on the original network and characteristic-line owner.

use super::*;
use fs_duct::network::{NetworkMedia, network_port_impedance_regional};

impl NetworkDrivingPoint {
    /// Solve a branched graph with a source record for each edge and an explicit
    /// exterior atmosphere. Gas-dependent scalar loads use that exterior.
    /// Geometry, material records, source domains and common static pressure are
    /// checked by the original network owner; no missing region is defaulted.
    pub fn response_regional(
        &self,
        media: &NetworkMedia,
        omega: f64,
        loss: LossModel,
    ) -> Result<NetworkDriveResponse, DuctError> {
        if !self.coupled_labels.is_empty() {
            return Err(DuctError::BadParameter {
                what: "network driving point with coupled ports needs their coupled termination",
            });
        }
        network_port_impedance_regional(&self.network, media, omega, loss)?
            .drive_bound(&self.active_port, &self.terminations)
    }

    /// The same regional propagation closed by one full, ordered multiport
    /// boundary. A supplied matrix still owns its own medium/geometry evidence;
    /// [`super::radiation::RadiatingNetwork::new_regional`] computes it from the
    /// retained exterior gas and uses that same BEM field at every microphone.
    pub fn response_coupled_regional(
        &self,
        media: &NetworkMedia,
        omega: f64,
        loss: LossModel,
        coupled: &CoupledTermination,
    ) -> Result<NetworkDriveResponse, DuctError> {
        self.admit_coupled_labels(coupled)?;
        network_port_impedance_regional(&self.network, media, omega, loss)?
            .drive_coupled_bound(&self.active_port, &self.terminations, coupled)
    }

    /// Realize a static, regionally filled graph through the original causal
    /// response ladder. Each edge contributes its actual L/c to the longest
    /// shortest travel path; composition/temperature therefore change both the
    /// sampled operator and its extracted propagation delay. `zc` is the caller's
    /// positive characteristic reference impedance, not an inferred new medium.
    /// This does not evolve mixing, heat transfer or a mean background flow.
    #[allow(clippy::too_many_arguments)]
    pub fn realize_line_regional(
        &self,
        media: &NetworkMedia,
        loss: LossModel,
        sample_rate_hz: u32,
        zc: f64,
        ladder: &LineLadder,
        mut coupled_at: Option<CoupledSource<'_>>,
    ) -> Result<(DelayedFilter, LineRealization), DrivingPointError> {
        if self.coupled_labels.is_empty() != coupled_at.is_none() {
            return Err(DrivingPointError::Invalid {
                what: "network line needs a coupled source exactly when coupled ports are declared",
            });
        }
        let geo_delay = self.regional_round_trip_seconds(media)? * f64::from(sample_rate_hz);
        let solve = |omega: f64| -> Result<C64, DuctError> {
            let response = match coupled_at.as_mut() {
                Some(source) => self.response_coupled_regional(media, omega, loss, &source(omega)?)?,
                None => self.response_regional(media, omega, loss)?,
            };
            Ok(reflectance(response.impedance, zc))
        };
        realize_response_line(geo_delay, sample_rate_hz, zc, ladder, solve)
    }

    /// Nominal bulk round trip [s] used by the regional causal realization:
    /// twice the largest shortest travel time from the active port to any node.
    /// Each edge contributes its geometric length divided by its own sound
    /// speed. This is a delay extraction scale, not a resonant decay time or a
    /// bound on all multiple-reflection paths in a graph with cycles.
    pub fn regional_round_trip_seconds(
        &self,
        media: &NetworkMedia,
    ) -> Result<f64, DrivingPointError> {
        media.validate_for(&self.network).map_err(DrivingPointError::Duct)?;
        Ok(2.0 * self.longest_path_by(|edge| {
            let length: f64 = edge.duct.segments.iter().map(Segment::axial_length).sum();
            length / media.edge_gas(&edge.id).expect("complete media binding admitted").sound_speed
        })?)
    }
}
