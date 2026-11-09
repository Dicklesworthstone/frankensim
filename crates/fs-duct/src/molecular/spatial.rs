//! Retain the original homogeneous molecular chain for interior observations.
//!
//! Block cuts and constitutive waves belong to the cold geometry binding.
//! Splitting a cylinder or spherical cone therefore uses its original wave,
//! while a profile probe uses an original Magnus-panel boundary. No new
//! molecular state, radial spectrum or interpolated material is introduced.

use super::local::absorption::MolecularCylinder;
use super::*;

pub(crate) struct PreparedMolecularSpatial<'source> {
    pub(crate) chain: PreparedChain,
    pub(crate) segment_boundaries: Vec<usize>,
    // None denotes an indivisible Magnus panel; observations there use cuts.
    waves: Vec<Option<SegmentWave>>,
    // Only authored cylinders have a homogeneous finite-volume loss law.
    cylinders: Vec<Option<MolecularCylinder<'source>>>,
}

impl<'source> PreparedMolecularSpatial<'source> {
    pub(crate) fn wave(&self, block: usize) -> Result<SegmentWave, MolecularDuctError> {
        self.waves
            .get(block)
            .copied()
            .flatten()
            .ok_or_else(|| bad("molecular observation needs an original local-wave block").into())
    }

    pub(crate) fn cylinder(&self, segment: usize) -> Option<MolecularCylinder<'source>> {
        self.cylinders.get(segment).copied().flatten()
    }
}

impl MolecularDuct {
    pub(crate) fn prepare_spatial_controlled(
        &self,
        omega: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedMolecularSpatial<'_>, MolecularDuctError> {
        poll(cancelled, 0)?;
        let boundaries = self.segment_boundaries.as_ref().ok_or_else(|| {
            bad("molecular graph observations require an original homogeneous edge")
        })?;
        if self.compliant_wall.is_some()
            || self.pieces.iter().any(|p| matches!(p, Piece::ToneHole(_)))
        {
            return Err(bad(
                "molecular graph observations require rigid walls and explicit chimney branches",
            )
            .into());
        }
        let mut capacities = Vec::with_capacity(self.points.len());
        let (waves, min_shear) =
            self.waves_observed_controlled(omega, cancelled, &mut |_index, cp| {
                capacities.push(cp)
            })?;
        let (chain, _) = self.chain_from_waves_controlled(
            omega,
            &[],
            ComposeOptions::default(),
            &waves,
            min_shear,
            cancelled,
        )?;
        let block_waves = self
            .pieces
            .iter()
            .map(|piece| match piece {
                Piece::Wave { point, .. } => Some(waves[*point]),
                Piece::Magnus { .. } => None,
                Piece::ToneHole(_) => unreachable!("compact holes refused before wave preparation"),
            })
            .collect();
        let cylinders = boundaries
            .windows(2)
            .map(|cut| match &self.pieces[cut[0]] {
                Piece::Wave {
                    geometry: Segment::Cylinder { radius, .. },
                    point,
                } if cut[1] == cut[0] + 1 => {
                    let original = &self.points[*point];
                    Some(MolecularCylinder {
                        state: self.state_at_point(*point),
                        wall: original.wall.as_ref(),
                        radius_m: *radius,
                        x_m: original.x_m,
                        omega,
                        cp: capacities[*point],
                        wave: waves[*point],
                    })
                }
                _ => None,
            })
            .collect();
        poll(cancelled, chain.blocks.len())?;
        Ok(PreparedMolecularSpatial {
            chain,
            segment_boundaries: boundaries.clone(),
            waves: block_waves,
            cylinders,
        })
    }
}
