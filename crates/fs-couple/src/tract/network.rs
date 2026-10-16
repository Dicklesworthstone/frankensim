//! A tract chart supplies graph DATA to the generic frequency-domain owner.
//! A branch is not presumed to represent a human specimen: its source is explicit.

use super::TractChart;
use fs_duct::network::{DuctNetwork, JunctionCorrections, NetworkEdge, NetworkNode, NodeBoundary};
use fs_duct::{Duct, DuctError, Segment};
use fs_math::det;

#[cfg(test)]
mod tests;

/// A declared branch at a chart-section boundary.
#[derive(Debug, Clone)]
pub struct TractBranch {
    /// Stable graph identity.
    pub id: String,
    /// 0 is the glottis end, sections().len() the oral end.
    pub after_section: usize,
    /// Complete branch duct geometry [m].
    pub duct: Duct,
    /// Explicit chain-local correction policy.
    pub corrections: JunctionCorrections,
    /// Label of the branch's external port.
    pub port_label: String,
    /// Citation/specimen identity, or explicitly declared synthetic geometry.
    pub source_id: String,
}

impl TractChart {
    /// Produce a generic N-port graph with ports "glottis", "oral" and every
    /// supplied branch label. All external loads are supplied by the consumer.
    /// No tract-specific propagation or junction formula is introduced.
    pub fn network_with_branches(
        &self,
        branches: &[TractBranch],
    ) -> Result<DuctNetwork, DuctError> {
        let count = self.sections().len();
        if count.saturating_add(1).saturating_add(branches.len())
            > fs_duct::network::MAX_NETWORK_NODES
            || branches
                .iter()
                .fold(count, |n, b| n.saturating_add(b.duct.segments.len()))
                > fs_duct::network::MAX_NETWORK_SEGMENTS
        {
            return Err(DuctError::BadParameter {
                what: "tract graph exceeds network node cap before allocation",
            });
        }
        let mut nodes = Vec::with_capacity(count + 1 + branches.len());
        for i in 0..=count {
            nodes.push(NetworkNode {
                id: format!("tract-section-{i:04}"),
                boundary: if i == 0 {
                    NodeBoundary::Port {
                        label: "glottis".into(),
                    }
                } else if i == count {
                    NodeBoundary::Port {
                        label: "oral".into(),
                    }
                } else {
                    NodeBoundary::Junction
                },
            });
        }
        let mut edges = Vec::with_capacity(count + branches.len());
        for (i, s) in self.sections().iter().enumerate() {
            edges.push(NetworkEdge {
                id: format!("tract-edge-{i:04}"),
                from: i,
                to: i + 1,
                duct: Duct {
                    segments: vec![Segment::Cylinder {
                        radius: det::sqrt(s.area_m2 / core::f64::consts::PI),
                        length: s.length_m,
                    }],
                },
                corrections: JunctionCorrections::Disabled,
            });
        }
        for b in branches {
            if b.after_section > count || b.id.trim().is_empty() || b.source_id.trim().is_empty() {
                return Err(DuctError::BadParameter {
                    what: "tract branch needs an in-range section boundary, identity and declared source",
                });
            }
            let end = nodes.len();
            nodes.push(NetworkNode {
                id: format!("branch-{}", b.id),
                boundary: NodeBoundary::Port {
                    label: b.port_label.clone(),
                },
            });
            edges.push(NetworkEdge {
                id: format!("branch-edge-{}", b.id),
                from: b.after_section,
                to: end,
                duct: b.duct.clone(),
                corrections: b.corrections,
            });
        }
        let network = DuctNetwork { nodes, edges };
        network.validate()?;
        Ok(network)
    }
}
