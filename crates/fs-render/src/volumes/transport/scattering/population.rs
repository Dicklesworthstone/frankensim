//! Complete elastic walks for absorption by evolving emitter populations.
//!
//! Native monochromatic free flights sample only the original elastic
//! scattering. Dynamic molecular and passive absorption are subsequently
//! integrated along every ordered segment by the finite-population owner.
//! This is the scattering-process/absorption-survival decomposition: no
//! frozen absorption probability or spectral mixture weight enters the walk.

use super::prepared::{Monochromatic, PaletteOptics, PreparedOptics, sample_extinction_observed};
use super::{PhaseFunction, ScatteringGrid, VolumeTransportError, unit_variate};
use crate::volumes::transport::cells::validate_ray;
use crate::volumes::transport::deposition::BeamPopulationSegment;
use crate::volumes::transport::{AbsorbingGrid, checkpoint, positive_product};
use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_geom::{Point3, Vec3};

/// Prepare each distinct native optical model once, independently of the
/// number of emitted histories. The original cell phase palette is borrowed.
pub(crate) struct PopulationPathSampler<'m, 'g, const M: usize> {
    optics: PaletteOptics<'m, 'g, M>,
    volume_m3: f64,
}

impl<'m, 'g, const M: usize> PopulationPathSampler<'m, 'g, M> {
    pub(crate) fn try_new(
        cx: &Cx<'_>,
        medium: &'m ScatteringGrid<'g>,
        wavelengths_nm: [f64; M],
    ) -> Result<Self, VolumeTransportError> {
        checkpoint(cx)?;
        let (lo, hi) = medium.wavelength_bounds_nm();
        if M == 0 || wavelengths_nm.iter().any(|w| !w.is_finite() || !(lo..=hi).contains(w)) {
            return Err(VolumeTransportError::Sample);
        }
        let optics = PaletteOptics::try_new(medium, wavelengths_nm, &mut || checkpoint(cx))?;
        let cell = medium.absorption_grid().grid.cell;
        let volume_m3 = positive_product(cell[0], cell[1], cell[2]);
        if !volume_m3.is_finite() || volume_m3 <= 0.0 {
            return Err(VolumeTransportError::Unrepresentable);
        }
        Ok(Self { optics, volume_m3 })
    }

    pub(crate) fn passive_absorption(
        &self,
        cell: usize,
        lane: usize,
    ) -> Result<f64, VolumeTransportError> {
        self.optics.coefficient(cell, lane).map(|(a, _, _)| a)
    }

    /// Keep all original cell visits, including returns to a previous cell and
    /// zero-density cells that can contain independently authored molecules.
    /// Both native DDA visits and retained segments are bounded by `max_cells`;
    /// exhausting either budget refuses the entire walk, never an escape.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn sample_path(
        &self,
        cx: &Cx<'_>,
        lane: usize,
        mut origin: Point3,
        mut direction: Vec3,
        maximum_scattering_events: u32,
        max_cells: usize,
        random: &mut impl FnMut() -> f64,
    ) -> Result<Vec<BeamPopulationSegment>, VolumeTransportError> {
        checkpoint(cx)?;
        validate_ray(origin, direction, f64::INFINITY)?;
        if lane >= M {
            return Err(VolumeTransportError::Sample);
        }
        let mono = Monochromatic::<_, M> {
            medium: &self.optics,
            lane,
            wavelength: [self.optics.wavelengths()[lane]],
        };
        let proposal = ScatteringOnly(&mono);
        let mut path = Vec::new();
        let mut visits = 0usize;
        let mut events = 0u32;
        loop {
            checkpoint(cx)?;
            let distance = random();
            if !unit_variate(distance) {
                return Err(VolumeTransportError::Sample);
            }
            let allowance = u32::try_from(max_cells - visits).unwrap_or(u32::MAX);
            let (flight, used) = sample_extinction_observed(
                &proposal,
                &mut || checkpoint(cx),
                origin,
                direction,
                f64::INFINITY,
                0.0,
                distance,
                false,
                &[],
                allowance,
                false,
                &mut |cell_index, length_m, _, _, _, _| {
                    if length_m == 0.0 {
                        return Ok(());
                    }
                    if path.len() == max_cells {
                        return Err(VolumeTransportError::Budget);
                    }
                    if path.len() == path.capacity() {
                        path.try_reserve_exact(path.capacity().max(16).min(max_cells - path.len()))
                            .map_err(|_| VolumeTransportError::Budget)?;
                    }
                    path.push(BeamPopulationSegment {
                        cell_index,
                        length_m,
                        cell_volume_m3: self.volume_m3,
                    });
                    Ok(())
                },
            )?;
            visits = visits
                .checked_add(used as usize)
                .filter(|&count| count <= max_cells)
                .ok_or(VolumeTransportError::Budget)?;
            // One physical wavelength has no spectral proposal correction.
            // Refuse a violated invariant rather than silently dropping weight.
            if flight.weight != [1.0] {
                return Err(VolumeTransportError::Unrepresentable);
            }
            let Some(event) = flight.event else {
                checkpoint(cx)?;
                return Ok(path);
            };
            if events == maximum_scattering_events {
                return Err(VolumeTransportError::Budget);
            }
            events += 1;
            let phase = proposal.phase(event.cell_index).sample_packet(
                proposal.wavelengths(),
                direction,
                random(),
                random(),
            )?;
            if phase.weight != [1.0] {
                return Err(VolumeTransportError::Unrepresentable);
            }
            origin = event.position;
            direction = phase.direction;
        }
    }
}

/// Borrow original scattering without changing the physical snapshot. Native
/// DDA, wavelength-dependent phase laws and free-flight inversion are shared.
struct ScatteringOnly<'a, O>(&'a O);

impl<O: PreparedOptics<1>> PreparedOptics<1> for ScatteringOnly<'_, O> {
    fn geometry(&self) -> &AbsorbingGrid<'_> {
        self.0.geometry()
    }
    fn identity(&self) -> ContentHash {
        self.0.identity()
    }
    fn bounds(&self) -> (f64, f64) {
        self.0.bounds()
    }
    fn wavelengths(&self) -> &[f64; 1] {
        self.0.wavelengths()
    }
    fn coefficient(&self, cell: usize, lane: usize) -> Result<(f64, f64, f64), VolumeTransportError> {
        let (_, scattering, _) = self.0.coefficient(cell, lane)?;
        Ok((0.0, scattering, scattering))
    }
    fn phase(&self, cell: usize) -> &PhaseFunction {
        self.0.phase(cell)
    }
    fn has_scattering(&self) -> bool {
        self.0.has_scattering()
    }
    fn lane_has_scattering(&self, lane: usize) -> bool {
        self.0.lane_has_scattering(lane)
    }
    fn lane_has_extinction(&self, lane: usize) -> bool {
        self.0.lane_has_scattering(lane)
    }
    fn has_emission(&self) -> bool {
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::tests::{budget, close, with_cx};
    use super::super::MassScatteringSpectrum;
    use crate::volumes::VolumeGrid;
    use crate::volumes::transport::MassAbsorptionSpectrum;

    fn medium<'a>(cx: &Cx<'_>, density: &'a [f64], passive: f64) -> ScatteringGrid<'a> {
        let absorption = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([2, 1, 1], density, [0.0; 3], [1.0; 3]),
            MassAbsorptionSpectrum::constant(passive).unwrap(),
            None,
            budget(8),
        )
        .unwrap();
        ScatteringGrid::try_new(
            cx,
            absorption,
            MassScatteringSpectrum::table(&[
                (380.0, 1.0), (500.0, 1.0), (650.0, 2.0), (780.0, 2.0),
            ])
            .unwrap(),
            PhaseFunction::Isotropic,
        )
        .unwrap()
    }

    #[test]
    fn g1_elastic_path_clips_collision_and_retains_return_through_vacuum() {
        with_cx(|cx, _| {
            let density = [0.0, 1.0];
            for passive in [0.0, 7.0] {
                let medium = medium(cx, &density, passive);
                let sampler = PopulationPathSampler::try_new(cx, &medium, [500.0, 650.0]).unwrap();
                for lane in 0..2 {
                    // Cross the transparent source cell, collide at optical
                    // depth 1/4, reverse exactly, then escape without another
                    // collision. Passive absorption must not change this walk.
                    let mut random = [(-0.25_f64).exp(), 0.0, 0.0, 0.0].into_iter();
                    let path = sampler
                        .sample_path(
                            cx, lane, Point3::new(0.25, 0.5, 0.5),
                            Vec3::new(1.0, 0.0, 0.0), 1, 4,
                            &mut || random.next().expect("four native flight/phase draws"),
                        )
                        .unwrap();
                    assert!(random.next().is_none());
                    assert_eq!(path.iter().map(|p| p.cell_index).collect::<Vec<_>>(), [0, 1, 1, 0]);
                    let local = 0.25 / (lane + 1) as f64;
                    for (segment, length) in path.iter().zip([0.75, local, local, 1.0]) {
                        close(segment.length_m, length, 2e-14);
                        assert_eq!(segment.cell_volume_m3, 1.0);
                    }
                }
            }
        });
    }

    #[test]
    fn g0_unfinished_elastic_paths_refuse_collision_and_visit_caps() {
        with_cx(|cx, _| {
            let density = [0.0, 1.0];
            let medium = medium(cx, &density, 3.0);
            let sampler = PopulationPathSampler::try_new(cx, &medium, [500.0]).unwrap();
            for (events, cells) in [(0, 4), (1, 3)] {
                let mut random = [(-0.25_f64).exp(), 0.0, 0.0, 0.0].into_iter();
                assert!(matches!(
                    sampler.sample_path(
                        cx, 0, Point3::new(0.25, 0.5, 0.5), Vec3::new(1.0, 0.0, 0.0),
                        events, cells, &mut || random.next().unwrap(),
                    ),
                    Err(VolumeTransportError::Budget)
                ));
            }
            let path = sampler
                .sample_path(
                    cx, 0, Point3::new(0.25, 0.5, 0.5), Vec3::new(1.0, 0.0, 0.0),
                    0, 2, &mut || 0.0,
                )
                .unwrap();
            assert_eq!(path.len(), 2);
            close(path[0].length_m, 0.75, 2e-14);
            close(path[1].length_m, 1.0, 2e-14);
        });
    }
}
