//! Finite spectral windows sharing the native slab thermal bindings.

use super::*;

/// Maximum explicitly authored finite bands in one native thermal adapter.
pub const MAX_SLAB_SPECTRAL_BANDS: usize = 64;

/// One original optical slab, constant over an explicit vacuum-wavelength band.
///
/// The layer absorption/index sources and exterior-index source must support
/// THIS wavelength window and their declared temperature intervals. This is
/// supplied spectral data, never a reinterpretation of a total-gray material.
#[derive(Clone, Debug, PartialEq)]
pub struct SlabThermalBand {
    pub window: BlackbodyBand,
    pub optics: AbsorbingSlab,
    /// Physical incident hemispherical power integrated over this band [W/m²]
    /// in this band's original exterior medium, front then back.
    pub irradiation_w_m2: [f64; 2],
}

/// Native enthalpy heating from explicitly covered, nonoverlapping bands.
///
/// Each band retains its original coefficients, index, source and temperature
/// support. Every band binds the SAME original layer cells, thicknesses, area,
/// masses and caloric laws. The original Planck owner supplies finite-window
/// emission, exact fixed-window tangents and stable temperature differences.
/// Gaps and tails remain OMITTED physics, visible through `bands()`; no gray
/// continuation, interpolation or inferred transparency fills them.
///
/// Band powers are summed before thermal integration. Each cell pair supplies
/// ONE conservative native transfer; only exterior exchange is external heat.
/// All band stiffnesses combine before tightening the native timestep. Angular
/// responses stay prepared; current emission is cached once per layer and band.
/// The observer is an explicit/lagged source; `endpoint_transport` instead feeds
/// the existing backward-Euler residual at its actual candidate temperatures.
#[derive(Clone)]
pub struct SpectralSlabThermal<'a> {
    bands: Arc<[SlabThermalBand]>,
    drivers: Arc<[GraySlabThermal<'a>]>,
    pairs: Arc<[(usize, usize)]>,
    derivative_bounds_w_k: Arc<BTreeMap<CellId, f64>>,
}

impl<'a> SpectralSlabThermal<'a> {
    /// Bind 1..=64 ordered, disjoint finite bands to the original active cells.
    /// Adjacent bands may touch; explicit gaps are allowed but receive no heat.
    /// Bandwise absorption and dispersion may differ. Physical layer geometry
    /// and the native cells/laws may not. No optical or caloric extrapolation.
    pub fn new(
        bands: Vec<SlabThermalBand>,
        area_m2: f64,
        cells: &[(CellId, &'a dyn ThermalLaw)],
        thermal: &DepositionSimulation<'_>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DepositionError> {
        checkpoint()?;
        if !(1..=MAX_SLAB_SPECTRAL_BANDS).contains(&bands.len()) {
            return Err(DepositionError::Budget("spectral slab band count"));
        }
        let n = bands[0].optics.layer_count();
        let mut drivers = Vec::with_capacity(bands.len());
        let mut pair_present = [false; MAX_SLAB_LAYERS * MAX_SLAB_LAYERS];
        let mut bounds = [0.0; MAX_SLAB_LAYERS];
        let mut correction = [0.0; MAX_SLAB_LAYERS];
        for (b, band) in bands.iter().enumerate() {
            checkpoint()?;
            if (b > 0
                && band.window.lower_wavelength_m() < bands[b - 1].window.upper_wavelength_m())
                || band.optics.layer_count() != n
                || band
                    .optics
                    .layers()
                    .iter()
                    .zip(bands[0].optics.layers())
                    .any(|(a, b)| a.thickness_m() != b.thickness_m())
            {
                return Err(DepositionError::InvalidInput(
                    "spectral slab needs ordered nonoverlapping bands on identical original layer geometry",
                ));
            }
            let driver = GraySlabThermal::bind(
                &band.optics,
                area_m2,
                band.irradiation_w_m2,
                cells,
                thermal,
                Some(band.window),
                &mut checkpoint,
            )?;
            for (i, &(id, _)) in cells.iter().enumerate() {
                add(
                    &mut bounds[i],
                    &mut correction[i],
                    driver.derivative_bounds_w_k[&id],
                );
            }
            for &(i, j, _) in driver.pairs.iter() {
                pair_present[i * n + j] = true;
            }
            drivers.push(driver);
        }
        let mut pairs = Vec::new();
        for i in 0..n {
            for j in i + 1..n {
                if pair_present[i * n + j] {
                    if pairs.len() == thermal.budget().max_links {
                        return Err(DepositionError::Budget(
                            "spectral slab physical transfer pairs",
                        ));
                    }
                    pairs.push((i, j));
                }
            }
        }
        if bounds[..n].iter().any(|x| !x.is_finite()) {
            return Err(DepositionError::Arithmetic);
        }
        let derivative_bounds_w_k = cells
            .iter()
            .enumerate()
            .map(|(i, &(id, _))| (id, bounds[i]))
            .collect();
        checkpoint()?;
        Ok(Self {
            bands: bands.into(),
            drivers: drivers.into(),
            pairs: pairs.into(),
            derivative_bounds_w_k: Arc::new(derivative_bounds_w_k),
        })
    }

    /// Original finite windows, physical band irradiances and sourced optics.
    pub fn bands(&self) -> &[SlabThermalBand] {
        &self.bands
    }
    pub fn area_m2(&self) -> f64 {
        self.drivers[0].area_m2
    }

    pub fn endpoint_transport(
        &self,
        thermal: &DepositionSimulation<'_>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<TransportStepCoefficients<()>, DepositionError> {
        self.evaluate(thermal, true, true, &mut checkpoint)
    }

    fn evaluate(
        &self,
        thermal: &DepositionSimulation<'_>,
        external: bool,
        internal: bool,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<TransportStepCoefficients<()>, DepositionError> {
        checkpoint()?;
        if self.pairs.len() > thermal.budget().max_links {
            return Err(DepositionError::Budget(
                "spectral slab physical transfer pairs",
            ));
        }
        let layers = &self.drivers[0].layers;
        let n = layers.len();
        let mut source = [0.0; MAX_SLAB_LAYERS];
        let mut source_correction = [0.0; MAX_SLAB_LAYERS];
        let mut pair = [0.0; MAX_SLAB_LAYERS * MAX_SLAB_LAYERS];
        let mut pair_correction = [0.0; MAX_SLAB_LAYERS * MAX_SLAB_LAYERS];
        for driver in self.drivers.iter() {
            let temperatures = driver.temperatures(thermal, checkpoint)?;
            let emitted = driver.band_exitances(&temperatures)?;
            if external {
                for (i, layer) in driver.layers.iter().enumerate() {
                    add(
                        &mut source[i],
                        &mut source_correction[i],
                        driver.exterior_power(layer, emitted[i])?,
                    );
                }
            }
            if internal {
                for (k, &(i, j, coefficient)) in driver.pairs.iter().enumerate() {
                    if k % 256 == 0 {
                        checkpoint()?;
                    }
                    let power = driver.pair_power(
                        coefficient,
                        temperatures[i],
                        temperatures[j],
                        [emitted[i], emitted[j]],
                    )?;
                    add(&mut pair[i * n + j], &mut pair_correction[i * n + j], power);
                }
            }
        }
        let mut sources_w = BTreeMap::new();
        for (i, layer) in layers.iter().enumerate() {
            if !source[i].is_finite() {
                return Err(DepositionError::Arithmetic);
            }
            if source[i] != 0.0 {
                sources_w.insert(layer.id, source[i]);
            }
        }
        let mut transfers = Vec::new();
        if internal {
            for &(i, j) in self.pairs.iter() {
                let power_w = pair[i * n + j];
                if !power_w.is_finite() {
                    return Err(DepositionError::Arithmetic);
                }
                if power_w != 0.0 {
                    transfers.push(HeatTransfer {
                        a: layers[i].id,
                        b: layers[j].id,
                        power_w,
                    });
                }
            }
        }
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

    fn validate(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DepositionError> {
        for driver in self.drivers.iter() {
            driver.temperatures(thermal, checkpoint)?;
        }
        Ok(())
    }
}

impl DepositionObserver for SpectralSlabThermal<'_> {
    type Error = DepositionError;
    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, DepositionError> {
        Ok(self
            .evaluate(thermal, true, false, checkpoint)?
            .heat
            .bulk
            .sources_w)
    }
    fn heat_transfers(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<HeatTransfer>, DepositionError> {
        Ok(self.evaluate(thermal, false, true, checkpoint)?.transfers)
    }
    fn limit_timestep_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        proposed_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, DepositionError> {
        self.validate(thermal, checkpoint)?;
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
        self.validate(thermal, &mut || Ok(()))
    }
    fn step(
        &mut self,
        _: f64,
        _: f64,
        _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DepositionError> {
        self.validate(thermal, checkpoint)
    }
}
