//! Native BDM1 volume transport for the original miscibility control volumes.
//!
//! The same sign-split normal trace used by momentum convection supplies both
//! upwind streams. An explicit depth converts the two-dimensional flow to SI
//! volume rates. Cell volumes, face areas and adjacency must match the original
//! mesh; authored diffusion distances are retained without an orthogonality claim.
//! The velocity is prescribed in the fixed volume frame. This adapter neither
//! advances momentum nor supplies capillary backreaction or a mass-frame change.

use super::advection::{BinaryMiscibilityAdvection, BinaryMiscibilityAdvectionFace};
use super::{BinaryMiscibilityError, BinaryMiscibilityNetwork, Sum, poll, positive};
use crate::ns::FluxSystem;
use crate::ns::transport::FluxTransportError;

/// An incompatible mesh/flow or a refusal by the original transport owner.
#[derive(Clone, Debug, PartialEq)]
pub enum BinaryMiscibilityFlowError {
    /// Mesh geometry or explicit depth does not describe the original graph.
    Geometry(&'static str),
    /// A closed miscibility graph cannot accept either exterior stream.
    OpenBoundary { edge: usize },
    /// Continuity of the original mean moments, before directional subtraction.
    VolumeDivergence {
        cell: usize,
        net_outward_m3_s: f64,
        maximum_stream_m3_s: f64,
    },
    /// The existing BDM1 trace could not supply representable physical rates.
    Transport(FluxTransportError),
    /// The finite-volume owner refused the graph/flow, or cancellation occurred.
    Miscibility(BinaryMiscibilityError),
}

impl core::fmt::Display for BinaryMiscibilityFlowError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "BDM1 miscibility transport: {self:?}")
    }
}
impl std::error::Error for BinaryMiscibilityFlowError {}
impl From<BinaryMiscibilityError> for BinaryMiscibilityFlowError {
    fn from(error: BinaryMiscibilityError) -> Self {
        Self::Miscibility(error)
    }
}
impl From<FluxTransportError> for BinaryMiscibilityFlowError {
    fn from(error: FluxTransportError) -> Self {
        Self::Transport(error)
    }
}

type Result<T> = core::result::Result<T, BinaryMiscibilityFlowError>;

impl BinaryMiscibilityNetwork {
    /// Bind the native normal velocity moments to this closed mixture graph.
    ///
    /// Cells follow the original triangle indices. Every interior edge must
    /// correspond to exactly one graph face, with the same area and adjacent
    /// cells; volumes and areas agree within 64 eps relative roundoff. No
    /// diffusion distance, material property or velocity is reconstructed.
    /// Both directions of a reversing face are retained, including zero-net
    /// counterflow. A boundary must have exactly zero transport in BOTH directions.
    ///
    /// Continuity is checked on the original mean moments as well as on the
    /// prepared directional rates, with the supplied allowance in [0,1e-12].
    /// This prevents small mean fluxes from disappearing through subtraction of
    /// large opposing streams. No continuity projection or rate repair occurs.
    /// The prepared result can be reused while this prescribed flow is frozen.
    ///
    /// Geometry/rate mapping is O(edges log faces), with O(cells + edges) scratch.
    /// Cancellation is checked around the existing native trace extraction and
    /// at each cell/edge of the mapping. A refusal publishes no prepared flow.
    pub fn prepare_advection_from_bdm(
        &self,
        system: &FluxSystem<'_>,
        velocity: &[f64],
        depth_m: f64,
        maximum_relative_volume_divergence: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BinaryMiscibilityAdvection> {
        poll(&mut cancelled)?;
        if !positive(depth_m)
            || !maximum_relative_volume_divergence.is_finite()
            || !(0.0..=1e-12).contains(&maximum_relative_volume_divergence)
        {
            return Err(BinaryMiscibilityFlowError::Geometry(
                "positive depth and bounded relative continuity allowance required",
            ));
        }
        let mesh = system.mesh;
        let n = self.owner.cells.len();
        if mesh.tris.len() != n
            || mesh.areas.len() != n
            || mesh.edges.len() > 3 * n
            || mesh.edges.len().checked_mul(2) != Some(system.n_u)
        {
            return Err(BinaryMiscibilityFlowError::Geometry(
                "one original finite volume per BDM1 triangle required",
            ));
        }
        for (cell, &area) in self.owner.cells.iter().zip(&mesh.areas) {
            poll(&mut cancelled)?;
            if !same_positive_measure(cell.volume_m3, area * depth_m) {
                return Err(BinaryMiscibilityFlowError::Geometry(
                    "finite-volume cell must retain original triangle area times depth",
                ));
            }
        }
        let transport = system.edge_volume_transport(velocity)?;
        poll(&mut cancelled)?;
        let mut rates = vec![
            BinaryMiscibilityAdvectionFace {
                forward_m3_s: 0.0,
                reverse_m3_s: 0.0,
            };
            self.owner.faces.len()
        ];
        let mut seen = vec![false; rates.len()];
        let mut net = vec![Sum::default(); n];
        let mut incoming = vec![Sum::default(); n];
        let mut outgoing = vec![Sum::default(); n];
        for (edge_index, (edge, flow)) in mesh.edges.iter().zip(&transport).enumerate() {
            poll(&mut cancelled)?;
            let (i, j) = edge.tris;
            if i >= n || (j != usize::MAX && (j >= n || i == j)) {
                return Err(BinaryMiscibilityFlowError::Geometry(
                    "valid original BDM1 edge adjacency required",
                ));
            }
            if j == usize::MAX {
                if flow.along_normal_m2_s() != 0.0
                    || flow.against_normal_m2_s() != 0.0
                    || flow.net_m2_s() != 0.0
                {
                    return Err(BinaryMiscibilityFlowError::OpenBoundary { edge: edge_index });
                }
                continue;
            }
            let pair = [i.min(j), i.max(j)];
            let index = self
                .owner
                .faces
                .binary_search_by_key(&pair, |face| face.physical.cells)
                .map_err(|_| {
                    BinaryMiscibilityFlowError::Geometry(
                        "every interior BDM1 edge needs its original miscibility face",
                    )
                })?;
            if seen[index]
                || !same_positive_measure(
                    self.owner.faces[index].physical.area_m2,
                    edge.len * depth_m,
                )
            {
                return Err(BinaryMiscibilityFlowError::Geometry(
                    "unique face with original edge length times depth required",
                ));
            }
            seen[index] = true;
            let along = dimensional_rate(flow.along_normal_m2_s(), depth_m)?;
            let against = dimensional_rate(flow.against_normal_m2_s(), depth_m)?;
            let mean = dimensional_rate(flow.net_m2_s(), depth_m)?;
            rates[index] = if i == pair[0] {
                BinaryMiscibilityAdvectionFace {
                    forward_m3_s: along,
                    reverse_m3_s: against,
                }
            } else {
                BinaryMiscibilityAdvectionFace {
                    forward_m3_s: against,
                    reverse_m3_s: along,
                }
            };
            net[i].add(mean)?;
            net[j].add(-mean)?;
            outgoing[i].add(along)?;
            incoming[j].add(along)?;
            outgoing[j].add(against)?;
            incoming[i].add(against)?;
        }
        if seen.iter().any(|present| !present) {
            return Err(BinaryMiscibilityFlowError::Geometry(
                "miscibility graph contains a face absent from the original BDM1 mesh",
            ));
        }
        for cell in 0..n {
            poll(&mut cancelled)?;
            let net_outward_m3_s = net[cell].value()?;
            let maximum_stream_m3_s = incoming[cell].value()?.max(outgoing[cell].value()?);
            if net_outward_m3_s.abs() > maximum_relative_volume_divergence * maximum_stream_m3_s {
                return Err(BinaryMiscibilityFlowError::VolumeDivergence {
                    cell,
                    net_outward_m3_s,
                    maximum_stream_m3_s,
                });
            }
        }
        Ok(self.prepare_advection(rates, maximum_relative_volume_divergence, cancelled)?)
    }
}

fn same_positive_measure(actual: f64, expected: f64) -> bool {
    positive(actual)
        && positive(expected)
        && (actual - expected).abs() / actual.max(expected) <= 64.0 * f64::EPSILON
}

fn dimensional_rate(rate_m2_s: f64, depth_m: f64) -> Result<f64> {
    let rate = rate_m2_s * depth_m;
    if !rate.is_finite() || (rate_m2_s != 0.0 && rate == 0.0) {
        return Err(BinaryMiscibilityError::Arithmetic("BDM1 physical volume rate").into());
    }
    Ok(rate)
}
