//! Prepared gray and finite-band optical exchange on original thermal cells.

use std::collections::BTreeMap;

use fs_thermal::enthalpy::backward_euler::{StepCoefficients, TransportStepCoefficients};
use fs_thermal::enthalpy::{
    CellGeometry, CellId, DepositionError, DepositionObserver, DepositionSimulation, Event,
    HeatTransfer, MaterialId, ThermalLaw,
};

use super::super::STEFAN_BOLTZMANN_W_M2_K4 as SIGMA;
use super::super::spectral::planck::BlackbodyBand;
use super::{
    AbsorbingSlab, Arc, MAX_SLAB_LAYERS, SlabCellResponse, add, exterior_difference_product,
};

mod spectral;
pub use spectral::{
    MAX_SLAB_SPECTRAL_BANDS, SlabThermalBand, SlabThermalCellBand, SpectralSlabThermal,
};
mod irradiation;
use irradiation::add_beam_interval;
pub use irradiation::{
    MAX_SLAB_THERMAL_BEAMS, SlabBeamIllumination, SlabBeamInterval, SlabThermalBeam,
};

#[derive(Clone)]
struct BoundLayer<'a> {
    id: CellId,
    optical_layer: usize,
    material: MaterialId,
    segment: u64,
    body: u64,
    mass_kg: f64,
    geometry: CellGeometry,
    caloric: &'a dyn ThermalLaw,
    face_area_m2: [f64; 2],
    radiative: bool,
}

/// Full-spectrum gray absorption/emission on existing uniform-temperature cells.
///
/// Prepared Fresnel, trapped-ray and interlayer coefficients are shared with
/// the original slab. No new material, optical interface, thermal storage or
/// conduction link is introduced. External absorbed/escaped radiation is a
/// native source; each internal exchange is ONE debit/credit `HeatTransfer` and
/// never external supply or conductive surface coverage.
///
/// As a `DepositionObserver`, emission/exchange is evaluated at each substep's
/// start; attached prescribed beams use their exact accepted-interval means.
/// A global bound over the original optical temperature supports combines
/// with the native conduction timestep bound. This is an explicit/lagged source,
/// including when used with an implicit conduction observer driver; refine time.
/// `endpoint_transport` instead evaluates the SAME exchange at a candidate
/// endpoint for the existing backward-Euler residual. Neither path supplies a
/// radiation entropy certificate or a new nonlinear thermal solver.
#[derive(Clone)]
pub struct GraySlabThermal<'a> {
    slab: AbsorbingSlab,
    cell_response: Option<SlabCellResponse>,
    layers: Arc<[BoundLayer<'a>]>,
    pairs: Arc<[(usize, usize, f64)]>,
    derivative_bounds_w_k: Arc<BTreeMap<CellId, f64>>,
    area_m2: f64,
    irradiation_w_m2: [f64; 2],
    // Only the private spectral binder selects a finite band. The public gray
    // constructor continues to declare full-spectrum optical source coverage.
    emission_band: Option<BlackbodyBand>,
    beams: Option<SlabBeamIllumination>,
}

impl<'a> GraySlabThermal<'a> {
    /// Bind a declared FULL-SPECTRUM GRAY slab to one distinct active native
    /// cell and its ORIGINAL caloric law per optical layer, front to back.
    ///
    /// This is the caller's explicit declaration that each optical source covers
    /// the full thermal spectrum. A finite-band/spectral slab MUST NOT be passed
    /// here: its emission requires the corresponding Planck-window integration.
    /// Irradiance is fixed hemispherical W/m² in the slab's original sourced
    /// exterior medium. At equilibrium it equals `n_exterior² * σ * T⁴`;
    /// reciprocal thermal escape retains that SAME exterior mode measure.
    /// Each original cell volume must equal `area * layer thickness`; mass,
    /// located geometry, source identities and caloric bindings remain fixed.
    /// Positive coefficients/bounds that cannot be represented refuse.
    pub fn new_full_spectrum(
        slab: &AbsorbingSlab,
        area_m2: f64,
        irradiation_w_m2: [f64; 2],
        cells: &[(CellId, &'a dyn ThermalLaw)],
        thermal: &DepositionSimulation<'_>,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DepositionError> {
        Self::bind(
            slab,
            area_m2,
            irradiation_w_m2,
            cells,
            thermal,
            None,
            None,
            checkpoint,
        )
    }

    /// Full-spectrum gray exchange on independently heated depth cells inside
    /// the original layers. Each native volume is `area * interval thickness`;
    /// each caloric and optical source remains original. This constructor makes
    /// the SAME full-spectrum source declaration as `new_full_spectrum`.
    pub fn new_cells_full_spectrum(
        response: &SlabCellResponse,
        area_m2: f64,
        irradiation_w_m2: [f64; 2],
        cells: &[(CellId, &'a dyn ThermalLaw)],
        thermal: &DepositionSimulation<'_>,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DepositionError> {
        Self::bind(
            response.slab(),
            area_m2,
            irradiation_w_m2,
            cells,
            thermal,
            None,
            Some(response),
            checkpoint,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn bind(
        slab: &AbsorbingSlab,
        area_m2: f64,
        irradiation_w_m2: [f64; 2],
        cells: &[(CellId, &'a dyn ThermalLaw)],
        thermal: &DepositionSimulation<'_>,
        emission_band: Option<BlackbodyBand>,
        cell_response: Option<&SlabCellResponse>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DepositionError> {
        checkpoint()?;
        if !area_m2.is_finite()
            || area_m2 <= 0.0
            || cells.len() != cell_response.map_or(slab.layer_count(), SlabCellResponse::cell_count)
            || irradiation_w_m2.iter().any(|x| !x.is_finite() || *x < 0.0)
        {
            return Err(DepositionError::InvalidInput(
                "gray slab needs positive area, nonnegative irradiation and one original cell/law per layer",
            ));
        }
        let n = cells.len();
        let face_coupling =
            cell_response.map_or(slab.face_coupling(), SlabCellResponse::face_coupling);
        let pair_coupling =
            cell_response.map_or(slab.layer_coupling(), SlabCellResponse::cell_coupling);
        let exterior_index_squared = slab
            .exterior_refractive_indices()
            .map(|index| index * index);
        let mut layers = Vec::with_capacity(n);
        let mut derivative_bounds_w_k = BTreeMap::new();
        let mut pairs = Vec::new();
        for (i, &(id, caloric)) in cells.iter().enumerate() {
            checkpoint()?;
            let cell = thermal
                .cells()
                .get(&id)
                .ok_or(DepositionError::InvalidInput(
                    "gray slab needs active original cells",
                ))?;
            let optical_layer = cell_response.map_or(i, |response| response.cells()[i].layer);
            let material = &slab.layers()[optical_layer];
            let thickness = cell_response.map_or(material.thickness_m(), |response| {
                let [start, end] = response.cells()[i].interval_m;
                end - start
            });
            let volume = product(area_m2, thickness)?;
            if cells[..i].iter().any(|&(previous, _)| previous == id)
                || !thermal.uses_material_law(cell.material_id, caloric)
                || (cell.geometry.volume_m3 / volume - 1.0).abs() > 64.0 * f64::EPSILON
            {
                return Err(DepositionError::InvalidInput(
                    "gray slab requires distinct cells, their original calorics and matching layer volumes",
                ));
            }
            let face = face_coupling[i];
            let face_area_m2 = [product(area_m2, face[0])?, product(area_m2, face[1])?];
            layers.push(BoundLayer {
                id,
                optical_layer,
                material: cell.material_id,
                segment: cell.segment_id,
                body: cell.body_id,
                mass_kg: cell.mass_kg,
                geometry: cell.geometry,
                caloric,
                face_area_m2,
                radiative: false,
            });
            let mut coupling = product(exterior_index_squared[0], face_area_m2[0])?
                + product(exterior_index_squared[1], face_area_m2[1])?;
            let mut correction = 0.0;
            for j in 0..n {
                let coefficient = product(area_m2, pair_coupling[i * n + j])?;
                add(&mut coupling, &mut correction, coefficient);
                if j > i && coefficient > 0.0 {
                    if pairs.len() == thermal.budget().max_links {
                        return Err(DepositionError::Budget("gray slab internal transfers"));
                    }
                    pairs.push((i, j, coefficient));
                }
            }
            layers[i].radiative = coupling > 0.0;
            // A fixed-wavelength Planck power is convex in T: its second
            // derivative has the sign of x*coth(x/2)-2 >= 0. The SAME positive
            // integral's tangent at source Tmax therefore bounds the band;
            // 4*E_band/T would be incorrect at fixed wavelength endpoints.
            let high = material.temperature_support_k()[1];
            let derivative = if coupling == 0.0 {
                0.0
            } else if let Some(band) = emission_band {
                product(
                    coupling,
                    band.evaluate(high)
                        .map_err(|_| DepositionError::Arithmetic)?
                        .temperature_derivative_w_m2_k,
                )?
            } else {
                let mut derivative = product(4.0 * SIGMA, coupling)?;
                for _ in 0..3 {
                    derivative = product(derivative, high)?;
                }
                derivative
            };
            derivative_bounds_w_k.insert(id, derivative);
        }
        let result = Self {
            slab: slab.clone(),
            cell_response: cell_response.cloned(),
            layers: layers.into(),
            pairs: pairs.into(),
            derivative_bounds_w_k: Arc::new(derivative_bounds_w_k),
            area_m2,
            irradiation_w_m2,
            emission_band,
            beams: None,
        };
        result.temperatures(thermal, &mut checkpoint)?;
        Ok(result)
    }

    pub fn slab(&self) -> &AbsorbingSlab {
        &self.slab
    }
    pub fn cell_response(&self) -> Option<&SlabCellResponse> {
        self.cell_response.as_ref()
    }
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }
    pub const fn irradiation_w_m2(&self) -> [f64; 2] {
        self.irradiation_w_m2
    }

    /// Attach scheduled directional heating on these SAME gray layers/cells.
    /// Beam power is additional to the authored diffuse irradiation. All
    /// sources select band zero and retain the full-spectrum gray declaration.
    /// Fixed optical coefficients give prescribed beam heat no thermal
    /// derivative; the original emission/exchange timestep bound still applies.
    pub fn with_beams(
        mut self,
        beams: Vec<SlabThermalBeam>,
        maximum_schedule_knots: usize,
        thermal: &DepositionSimulation<'_>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DepositionError> {
        self.temperatures(thermal, &mut checkpoint)?;
        self.beams = Some(SlabBeamIllumination::prepare(
            beams,
            maximum_schedule_knots,
            std::slice::from_ref(&self),
            thermal.time_s(),
            &mut checkpoint,
        )?);
        Ok(self)
    }

    pub fn beam_illumination(&self) -> Option<&SlabBeamIllumination> {
        self.beams.as_ref()
    }

    /// Original native endpoint transport, ready for
    /// `BackwardEulerStep::evaluate_with_transport`. Physical sources and
    /// pairwise exchanges retain separate native energy-accounting owners.
    pub fn endpoint_transport(
        &self,
        thermal: &DepositionSimulation<'_>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<TransportStepCoefficients<()>, DepositionError> {
        if self.beams.is_some() {
            return Err(DepositionError::InvalidInput(
                "scheduled slab beams require endpoint_transport_over_interval",
            ));
        }
        self.thermal_endpoint_transport(thermal, &mut checkpoint)
    }

    /// Existing BE transport with emission at the candidate END temperature
    /// and beam power integrated over the ACTUAL native interval. Using the
    /// endpoint pulse amplitude would supply incorrect joules. The candidate
    /// clock must equal `end_s`; this query advances no optical history.
    pub fn endpoint_transport_over_interval(
        &self,
        thermal: &DepositionSimulation<'_>,
        start_s: f64,
        duration_s: f64,
        end_s: f64,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<TransportStepCoefficients<()>, DepositionError> {
        if thermal.time_s() != end_s {
            return Err(DepositionError::InvalidInput(
                "slab BE candidate must have the actual endpoint time",
            ));
        }
        let mut result = self.thermal_endpoint_transport(thermal, &mut checkpoint)?;
        add_beam_interval(
            self.beams.as_ref(),
            &mut result.heat.bulk.sources_w,
            start_s,
            duration_s,
            end_s,
            &mut checkpoint,
        )?;
        Ok(result)
    }

    fn thermal_endpoint_transport(
        &self,
        thermal: &DepositionSimulation<'_>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<TransportStepCoefficients<()>, DepositionError> {
        let temperatures = self.temperatures(thermal, &mut checkpoint)?;
        let sources_w = self.external_sources(&temperatures)?;
        let transfers = self.internal_transfers(&temperatures, &mut checkpoint)?;
        checkpoint()?;
        Ok(TransportStepCoefficients {
            heat: StepCoefficients {
                sources_w,
                conductances_w_k: BTreeMap::new(),
                payload: (),
            }
            .into(),
            transfers,
        })
    }

    fn temperatures(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<[f64; MAX_SLAB_LAYERS], DepositionError> {
        checkpoint()?;
        if self.pairs.len() > thermal.budget().max_links {
            return Err(DepositionError::Budget("gray slab internal transfers"));
        }
        let mut result = [0.0; MAX_SLAB_LAYERS];
        for (i, original) in self.layers.iter().enumerate() {
            let cell = thermal
                .cells()
                .get(&original.id)
                .ok_or(DepositionError::InvalidInput(
                    "gray slab original cell is no longer active",
                ))?;
            if cell.material_id != original.material
                || cell.segment_id != original.segment
                || cell.body_id != original.body
                || cell.mass_kg != original.mass_kg
                || cell.geometry != original.geometry
                || !thermal.uses_material_law(original.material, original.caloric)
            {
                return Err(DepositionError::InvalidInput(
                    "gray slab original thermal binding changed",
                ));
            }
            result[i] = thermal.temperature_k(original.id)?;
            let [lo, hi] = self.slab.layers()[original.optical_layer].temperature_support_k();
            if !result[i].is_finite() || !(lo..=hi).contains(&result[i]) {
                return Err(DepositionError::MaterialDomain);
            }
        }
        if let Some(beams) = &self.beams {
            beams.admit_time(thermal.time_s())?;
        }
        Ok(result)
    }

    fn external_sources(
        &self,
        temperatures: &[f64],
    ) -> Result<BTreeMap<CellId, f64>, DepositionError> {
        let mut sources = BTreeMap::new();
        for (layer, &temperature) in self.layers.iter().zip(temperatures) {
            if layer.face_area_m2 == [0.0; 2] {
                continue;
            }
            let power = self.exterior_power(layer, self.exitance(temperature)?)?;
            if power != 0.0 {
                sources.insert(layer.id, power);
            }
        }
        Ok(sources)
    }

    fn exitance(&self, temperature: f64) -> Result<f64, DepositionError> {
        if let Some(band) = self.emission_band {
            return band
                .evaluate(temperature)
                .map(|value| value.power_w_m2)
                .map_err(|_| DepositionError::Arithmetic);
        }
        let square = product(temperature, temperature)?;
        product(product(SIGMA, square)?, square)
    }

    fn exterior_power(
        &self,
        layer: &BoundLayer<'_>,
        emission: f64,
    ) -> Result<f64, DepositionError> {
        let exterior_index_squared = self
            .slab
            .exterior_refractive_indices()
            .map(|index| index * index);
        let mut power = 0.0;
        let mut correction = 0.0;
        for side in 0..2 {
            add(
                &mut power,
                &mut correction,
                exterior_difference_product(
                    layer.face_area_m2[side],
                    self.irradiation_w_m2[side],
                    emission,
                    exterior_index_squared[side],
                )
                .map_err(|_| DepositionError::Arithmetic)?,
            );
        }
        if power.is_finite() {
            Ok(power)
        } else {
            Err(DepositionError::Arithmetic)
        }
    }

    fn band_exitances(
        &self,
        temperatures: &[f64],
    ) -> Result<[f64; MAX_SLAB_LAYERS], DepositionError> {
        let mut values = [0.0; MAX_SLAB_LAYERS];
        if self.emission_band.is_some() {
            for (i, layer) in self.layers.iter().enumerate() {
                if layer.radiative {
                    values[i] = self.exitance(temperatures[i])?;
                }
            }
        }
        Ok(values)
    }

    fn pair_power(
        &self,
        coefficient: f64,
        a: f64,
        b: f64,
        emitted: [f64; 2],
    ) -> Result<f64, DepositionError> {
        if a == b {
            return Ok(0.0);
        }
        if let Some(band) = self.emission_band {
            let gap = emitted[0] - emitted[1];
            let scale = emitted[0].max(emitted[1]);
            // Well-separated cached powers subtract safely. Nearby powers use
            // the original positive Planck difference, preserving adjacent T.
            let difference = if scale > 0.0 && gap.abs() >= 0.5 * scale {
                gap
            } else {
                band.difference(a, b)
                    .map_err(|_| DepositionError::Arithmetic)?
            };
            return product(coefficient, difference);
        }
        // Preserve the original gray factorization and equilibrium arithmetic.
        let secant = product(product(product(coefficient, SIGMA)?, a + b)?, a * a + b * b)?;
        product(secant, a - b)
    }

    fn internal_transfers(
        &self,
        temperatures: &[f64],
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<HeatTransfer>, DepositionError> {
        let mut transfers = Vec::with_capacity(self.pairs.len());
        let emitted = self.band_exitances(temperatures)?;
        for (k, &(i, j, coefficient)) in self.pairs.iter().enumerate() {
            if k % 256 == 0 {
                checkpoint()?;
            }
            let (a, b) = (temperatures[i], temperatures[j]);
            if a == b {
                continue;
            }
            transfers.push(HeatTransfer {
                a: self.layers[i].id,
                b: self.layers[j].id,
                power_w: self.pair_power(coefficient, a, b, [emitted[i], emitted[j]])?,
            });
        }
        Ok(transfers)
    }
}

impl DepositionObserver for GraySlabThermal<'_> {
    type Error = DepositionError;
    fn requires_step_interval(&self) -> bool {
        self.beams.is_some()
    }
    fn next_transition_time_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Option<f64>, DepositionError> {
        checkpoint()?;
        if let Some(beams) = &self.beams {
            beams.admit_time(thermal.time_s())?;
            Ok(beams.next_transition(thermal.time_s()))
        } else {
            Ok(None)
        }
    }
    fn step_sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        start_s: f64,
        duration_s: f64,
        end_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, DepositionError> {
        if self.beams.is_none() {
            return Ok(BTreeMap::new());
        }
        self.temperatures(thermal, checkpoint)?;
        if thermal.time_s() != start_s {
            return Err(DepositionError::InvalidInput(
                "slab beam interval must start at native time",
            ));
        }
        let mut sources = BTreeMap::new();
        add_beam_interval(
            self.beams.as_ref(),
            &mut sources,
            start_s,
            duration_s,
            end_s,
            checkpoint,
        )?;
        Ok(sources)
    }
    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, DepositionError> {
        let temperatures = self.temperatures(thermal, checkpoint)?;
        self.external_sources(&temperatures)
    }
    fn heat_transfers(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<HeatTransfer>, DepositionError> {
        let temperatures = self.temperatures(thermal, checkpoint)?;
        self.internal_transfers(&temperatures, checkpoint)
    }
    fn limit_timestep_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        proposed_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, DepositionError> {
        self.temperatures(thermal, checkpoint)?;
        thermal.limit_timestep_for_additional_conductances(
            proposed_s,
            &self.derivative_bounds_w_k,
            checkpoint,
        )
    }
    fn event(
        &mut self,
        _: &Event,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<(), DepositionError> {
        self.temperatures(thermal, &mut || Ok(())).map(|_| ())
    }
    fn step(
        &mut self,
        _: f64,
        _: f64,
        _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DepositionError> {
        self.temperatures(thermal, checkpoint).map(|_| ())
    }
}

fn product(a: f64, b: f64) -> Result<f64, DepositionError> {
    let value = a * b;
    if !value.is_finite() || (value == 0.0 && a != 0.0 && b != 0.0) {
        Err(DepositionError::Arithmetic)
    } else {
        Ok(value)
    }
}
