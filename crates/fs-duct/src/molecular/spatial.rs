//! Retain the original homogeneous molecular chain for interior observations.
//!
//! Block cuts and constitutive waves belong to the cold geometry binding.
//! Splitting a cylinder or spherical cone therefore uses its original wave,
//! while a profile probe uses an original Magnus-panel boundary. No new
//! molecular state, radial spectrum or interpolated material is introduced.

use super::*;

pub(crate) struct PreparedMolecularSpatial {
    pub(crate) chain: PreparedChain,
    pub(crate) segment_boundaries: Vec<usize>,
    // None denotes an indivisible Magnus panel; observations there use cuts.
    waves: Vec<Option<SegmentWave>>,
}

impl PreparedMolecularSpatial {
    pub(crate) fn wave(&self, block: usize) -> Result<SegmentWave, MolecularDuctError> {
        self.waves
            .get(block)
            .copied()
            .flatten()
            .ok_or_else(|| bad("molecular observation needs an original local-wave block").into())
    }
}

impl MolecularDuct {
    pub(crate) fn prepare_spatial_controlled(
        &self,
        omega: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedMolecularSpatial, MolecularDuctError> {
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
        let (waves, min_shear) = self.waves_controlled(omega, cancelled)?;
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
        poll(cancelled, chain.blocks.len())?;
        Ok(PreparedMolecularSpatial {
            chain,
            segment_boundaries: boundaries.clone(),
            waves: block_waves,
        })
    }
}
