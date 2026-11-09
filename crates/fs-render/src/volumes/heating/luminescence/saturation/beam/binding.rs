//! Original native geometry and immutable beam declarations.

use super::*;
use std::collections::{BTreeMap, BTreeSet};

/// One independent, mutually incoherent two-wavelength pump. Powers are watts
/// already integrated over the authored beam and spectral quadrature weights.
/// Repeated source identities do not merge rays. The optional half-open pulse
/// intersects the owner's global pump window.
#[derive(Clone, Debug)]
pub struct DepletingPumpBeam<'a> {
    pub source_identity: ContentHash,
    pub beam: IncidentBeam<2>,
    pub cross_sections: &'a [BeamEmitterCrossSections],
    pub pump_window_s: Option<(f64, f64)>,
}

impl<const M: usize> SaturableLuminescentHeating<M> {
    /// Bind one directed pump, preserving the original one-column numerical
    /// path. Use Self::with_depleting_beams for simultaneous illumination.
    #[allow(clippy::too_many_arguments)]
    pub fn with_depleting_beam(
        self,
        cx: &Cx<'_>,
        medium: &AbsorbingGrid<'_>,
        source_identity: ContentHash,
        beam: IncidentBeam<2>,
        cross_sections: &[BeamEmitterCrossSections],
        thermal: &DepositionSimulation<'_>,
        max_path_cells: usize,
        limits: PumpColumnLimits,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, OpticalHeatingError> {
        self.with_depleting_beams(
            cx, medium,
            &[DepletingPumpBeam { source_identity, beam, cross_sections, pump_window_s: None }],
            thermal, 1, max_path_cells,
            PumpBeamsLimits { column: limits, ..PumpBeamsLimits::default() },
            checkpoint,
        )
    }

    /// Bind simultaneous pumps to one shared finite inventory per original
    /// cell. Every beam follows its original DDA order, including opposing or
    /// crossing directions; a cell is never independently advanced per beam.
    /// Bind before advancing, with all prescribed local pump hazards zero.
    ///
    /// Each beam has its own wavelengths, source cross sections and optional
    /// pulse. Both pulse edges split native thermal intervals. Illumination is
    /// intensity-additive: optical phase, interference and coherent population
    /// dynamics are outside this rate model. The cell-uniform approximation
    /// requires spatial refinement when different beam footprints address
    /// materially different populations inside one cell.
    ///
    /// Host mass absorption is passive and must exclude these same active
    /// emitters. Each traversed emitter needs explicit molecular cross sections
    /// for that beam. A passive absorbing cell needs a native thermal receiver,
    /// which may be declared with zero populations. Transparent cells without
    /// inventories need no binding. The total path cap counts ALL traversed
    /// segments, including repeated visits by different beams.
    #[allow(clippy::too_many_arguments)]
    pub fn with_depleting_beams(
        mut self,
        cx: &Cx<'_>,
        medium: &AbsorbingGrid<'_>,
        beams: &[DepletingPumpBeam<'_>],
        thermal: &DepositionSimulation<'_>,
        max_beams: usize,
        max_total_path_cells: usize,
        limits: PumpBeamsLimits,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, OpticalHeatingError> {
        self.check_state(thermal, true, &mut checkpoint)?;
        if self.energy.duration_s != 0.0 || self.pending.is_some() || self.beam.is_some()
            || beams.is_empty() || beams.len() > max_beams
            || medium.content_identity() != self.medium_identity
            || self.emitters.iter().any(|e| e.binding.pump.rates_s != [0.0; 2])
        {
            return Err(invalid("invalid simultaneous pump binding, medium or prescribed hazards"));
        }
        let mut hash = DomainHasher::new("org.frankensim.render.depleting-luminescence-beams.v1");
        hash.update(self.source_identity.as_bytes());
        hash.update(&(beams.len() as u64).to_le_bytes());
        let mut paths = Vec::with_capacity(beams.len());
        let mut visited = 0usize;
        let mut emitters = BTreeSet::new();
        for declaration in beams {
            checkpoint()?;
            let (path, count) = compile_path(
                &self, cx, medium, declaration, max_total_path_cells - visited,
                &mut hash, &mut checkpoint,
            )?;
            visited += count;
            emitters.extend(path.cells.iter().map(|cell| cell.cell_index));
            paths.push(path);
        }
        let cells: Vec<_> = emitters.into_iter().collect();
        for path in &mut paths {
            for cell in &mut path.cells {
                checkpoint()?;
                cell.cell_index = cells.binary_search(&cell.cell_index)
                    .expect("collected original pump inventory");
            }
        }
        let accounts = BeamAccounts::new(&paths);
        let mut pump = DepletingBeam {
            model: Arc::new(BeamModel { cells, paths, limits }),
            prepared: None,
            accounts,
        };
        pump.prepared = Some(pump.prepare(&self, thermal, &mut checkpoint)?);
        checkpoint()?;
        self.source_identity = hash.finalize();
        self.beam = Some(pump);
        Ok(self)
    }

    /// Original one-beam account. Absent for prescribed local fields and for
    /// multiple beams, whose different wavelengths require separate accounts.
    pub fn depleting_pump_energy(&self) -> Option<&DepletingPumpEnergy> {
        self.beam.as_ref().and_then(|beam| beam.accounts.single())
    }

    /// Simultaneous accepted pump accounts in original declaration order,
    /// with whole-system work counted once. Also present for one beam.
    pub fn depleting_pump_report(&self) -> Option<&DepletingPumpReport> {
        self.beam.as_ref().map(|beam| beam.accounts.report())
    }
}

#[allow(clippy::too_many_arguments)]
fn compile_path<const M: usize>(
    owner: &SaturableLuminescentHeating<M>,
    cx: &Cx<'_>,
    medium: &AbsorbingGrid<'_>,
    declaration: &DepletingPumpBeam<'_>,
    remaining_segments: usize,
    hash: &mut DomainHasher,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<(BeamPath, usize), OpticalHeatingError> {
    let beam = &declaration.beam;
    if declaration.source_identity == ContentHash([0; 32])
        || !beam.wavelengths_nm.into_iter().all(|v| v.is_finite() && v > 0.0)
        || !beam.incident_power_w.into_iter().all(nonnegative)
        || declaration.cross_sections.len() > remaining_segments
        || declaration.pump_window_s.is_some_and(|(start, end)| {
            !start.is_finite() || !end.is_finite() || start < owner.time_s || end <= start
        })
    {
        return Err(invalid("invalid beam identity, spectrum, pulse or path budget"));
    }
    let segments = medium.population_beam_cell_path(
        cx, beam.origin, beam.unit_direction, beam.max_distance_m, remaining_segments,
    ).map_err(|error| match error {
        VolumeTransportError::Cancelled => OpticalHeatingError::Thermal(DepositionError::Cancelled),
        _ => invalid("depleting pump geometry or aggregate traversal budget"),
    })?;
    let count = segments.len();
    let mut bindings = BTreeMap::new();
    for binding in declaration.cross_sections {
        checkpoint()?;
        if bindings.insert(binding.grid_cell, binding.clone()).is_some() {
            return Err(OpticalHeatingError::DuplicateCell(binding.grid_cell));
        }
    }
    let mut photon_energies_j = [0.0; 2];
    let mut incident_photon_rates_s = [0.0; 2];
    hash.update(declaration.source_identity.as_bytes());
    hash.update(&[u8::from(declaration.pump_window_s.is_some())]);
    if let Some((start, end)) = declaration.pump_window_s {
        hash.update(&start.to_bits().to_le_bytes());
        hash.update(&end.to_bits().to_le_bytes());
    }
    for lane in 0..2 {
        photon_energies_j[lane] = finite(super::super::HC_J_NM / beam.wavelengths_nm[lane])?;
        if photon_energies_j[lane] <= 0.0 {
            return Err(invalid("unrepresentable pump photon energy"));
        }
        incident_photon_rates_s[lane] = finite(beam.incident_power_w[lane] / photon_energies_j[lane])?;
        if beam.incident_power_w[lane] > 0.0 && incident_photon_rates_s[lane] == 0.0 {
            return Err(invalid("positive pump photon rate is below representation"));
        }
        hash.update(&beam.wavelengths_nm[lane].to_bits().to_le_bytes());
        hash.update(&beam.incident_power_w[lane].to_bits().to_le_bytes());
    }
    hash.update(&(count as u64).to_le_bytes());
    let mut cells = Vec::with_capacity(count);
    for segment in segments {
        checkpoint()?;
        let coefficients = medium.coefficients_at_cell(segment.cell_index, &beam.wavelengths_nm)
            .map_err(|_| invalid("passive host does not support pump wavelength"))?;
        let density = medium.cell_density_kg_m3(segment.cell_index)
            .ok_or(OpticalHeatingError::UnmappedCell(segment.cell_index))?;
        let mut passive_optical_depth = [0.0; 2];
        for lane in 0..2 {
            if coefficients[lane] > 0.0 && density > 0.0 {
                passive_optical_depth[lane] = finite(positive_product(
                    coefficients[lane], density, segment.length_m,
                ))?;
                if passive_optical_depth[lane] == 0.0 {
                    return Err(invalid("positive passive pump optical depth is below representation"));
                }
            }
        }
        let cell_index = match owner.emitters.binary_search_by_key(
            &segment.cell_index, |e| e.binding.grid_cell,
        ) {
            Ok(index) => index,
            Err(_) if passive_optical_depth == [0.0; 2] => continue,
            Err(_) => return Err(OpticalHeatingError::UnmappedCell(segment.cell_index)),
        };
        let source = bindings.remove(&segment.cell_index)
            .ok_or(OpticalHeatingError::UnmappedCell(segment.cell_index))?;
        hash.update(&(segment.cell_index as u64).to_le_bytes());
        hash.update(&segment.length_m.to_bits().to_le_bytes());
        hash.update(&segment.cell_volume_m3.to_bits().to_le_bytes());
        hash.update(source.identity().as_bytes());
        for depth in passive_optical_depth { hash.update(&depth.to_bits().to_le_bytes()); }
        match &source.source {
            CrossSections::Authored { values, .. } => {
                for value in values.iter().flatten() { hash.update(&value.to_bits().to_le_bytes()); }
            }
            #[cfg(feature = "material-state")]
            CrossSections::Material { .. } => {}
        }
        cells.push(BeamCell {
            cell_index, length_m: segment.length_m, volume_m3: segment.cell_volume_m3,
            passive_optical_depth, cross_sections: source,
        });
    }
    if !bindings.is_empty() {
        return Err(invalid("pump cross-section binding is outside the traversed inventory"));
    }
    // Delimit the inventory-bearing path from skipped transparent segments.
    hash.update(&(cells.len() as u64).to_le_bytes());
    Ok((BeamPath {
        source_identity: declaration.source_identity,
        wavelengths_nm: beam.wavelengths_nm,
        photon_energies_j,
        incident_photon_rates_s,
        pump_window_s: declaration.pump_window_s,
        cells,
    }, count))
}
