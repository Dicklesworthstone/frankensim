//! Retain the original homogeneous molecular chain for interior observations.
//!
//! Block cuts and constitutive waves belong to the cold geometry binding.
//! Splitting a cylinder or spherical cone therefore uses its original wave,
//! while a profile probe uses an original Magnus-panel boundary. No new
//! molecular state, radial spectrum or interpolated material is introduced.

use super::local::absorption::{MolecularCone, MolecularCylinder};
use super::*;

pub(crate) struct PreparedMolecularSpatial<'source> {
    pub(crate) chain: PreparedChain,
    pub(crate) segment_boundaries: Vec<usize>,
    // None denotes an indivisible Magnus panel; observations there use cuts.
    waves: Vec<Option<SegmentWave>>,
    // Only authored cylinders have a homogeneous finite-volume loss law.
    cylinders: Vec<Option<MolecularCylinder<'source>>>,
    // Indexed by the ORIGINAL wave block, never by a new observation grid.
    cones: Vec<Option<MolecularCone<'source>>>,
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

    pub(crate) fn cone(&self, block: usize) -> Option<MolecularCone<'source>> {
        self.cones.get(block).copied().flatten()
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
        if self.pieces.iter().any(|p| matches!(p, Piece::ToneHole(_))) {
            return Err(
                bad("molecular graph observations require explicit chimney branches").into(),
            );
        }
        let mut capacities = Vec::with_capacity(self.points.len());
        let (waves, min_shear) =
            self.waves_observed_controlled(omega, cancelled, &mut |_index, cp| {
                capacities.push(cp)
            })?;
        let mut wall_waves = BTreeMap::new();
        let (chain, _) = self.chain_from_waves_cached_controlled(
            omega,
            &[],
            ComposeOptions::default(),
            &waves,
            min_shear,
            &mut wall_waves,
            cancelled,
        )?;
        let mut block_waves = Vec::with_capacity(self.pieces.len());
        let mut cones = Vec::with_capacity(self.pieces.len());
        for (index, piece) in self.pieces.iter().enumerate() {
            poll(cancelled, index)?;
            block_waves.push(match piece {
                Piece::Wave { geometry, point } => Some(self.wall_wave(
                    *point,
                    self.wave_slant(geometry),
                    omega,
                    &waves,
                    &mut wall_waves,
                )?),
                Piece::Magnus { .. } => None,
                Piece::ToneHole(_) => unreachable!("compact holes refused before wave preparation"),
            });
            cones.push(match piece {
                Piece::Wave {
                    geometry:
                        geometry @ Segment::Cone {
                            inlet_radius,
                            outlet_radius,
                            length,
                        },
                    point,
                } => {
                    let original = &self.points[*point];
                    let response = self.wall_response(
                        *point,
                        self.wave_slant(geometry),
                        omega,
                        &waves,
                        &mut wall_waves,
                    )?;
                    Some(MolecularCone {
                        inlet_radius_m: *inlet_radius,
                        outlet_radius_m: *outlet_radius,
                        length_m: *length,
                        original: MolecularCylinder {
                            state: self.state_at_point(*point),
                            wall: original.wall.as_ref(),
                            radius_m: original.radius,
                            x_m: original.x_m,
                            omega,
                            cp: capacities[*point],
                            wave: response.wave,
                            mechanical_admittance_per_metre: response.admittance_per_metre,
                        },
                    })
                }
                _ => None,
            });
        }
        let mut cylinders = Vec::with_capacity(boundaries.len().saturating_sub(1));
        for (index, cut) in boundaries.windows(2).enumerate() {
            poll(cancelled, index)?;
            cylinders.push(match &self.pieces[cut[0]] {
                Piece::Wave {
                    geometry: Segment::Cylinder { radius, .. },
                    point,
                } if cut[1] == cut[0] + 1 => {
                    let original = &self.points[*point];
                    let response =
                        self.wall_response(*point, 1.0, omega, &waves, &mut wall_waves)?;
                    Some(MolecularCylinder {
                        state: self.state_at_point(*point),
                        wall: original.wall.as_ref(),
                        radius_m: *radius,
                        x_m: original.x_m,
                        omega,
                        cp: capacities[*point],
                        wave: response.wave,
                        mechanical_admittance_per_metre: response.admittance_per_metre,
                    })
                }
                _ => None,
            });
        }
        poll(cancelled, chain.blocks.len())?;
        Ok(PreparedMolecularSpatial {
            chain,
            segment_boundaries: boundaries.clone(),
            waves: block_waves,
            cylinders,
            cones,
        })
    }
}
