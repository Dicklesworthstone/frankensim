//! Prepared gray optical exchange on the original native thermal cells.

use std::collections::BTreeMap;

use fs_thermal::enthalpy::backward_euler::{StepCoefficients, TransportStepCoefficients};
use fs_thermal::enthalpy::{
    CellGeometry, CellId, DepositionError, DepositionObserver, DepositionSimulation, Event,
    HeatTransfer, MaterialId, ThermalLaw,
};

use super::super::STEFAN_BOLTZMANN_W_M2_K4 as SIGMA;
use super::{AbsorbingSlab, Arc, MAX_SLAB_LAYERS, add, exterior_difference_product};

#[derive(Clone)]
struct BoundLayer<'a> {
    id: CellId,
    material: MaterialId,
    segment: u64,
    body: u64,
    mass_kg: f64,
    geometry: CellGeometry,
    caloric: &'a dyn ThermalLaw,
    face_area_m2: [f64; 2],
}

/// Full-spectrum gray absorption/emission on existing uniform-temperature cells.
///
/// Prepared Fresnel, trapped-ray and interlayer coefficients are shared with
/// the original slab. No new material, optical interface, thermal storage or
/// conduction link is introduced. External absorbed/escaped radiation is a
/// native source; each internal exchange is ONE debit/credit `HeatTransfer` and
/// never external supply or conductive surface coverage.
///
/// As a `DepositionObserver`, radiation is evaluated at each accepted substep's
/// start. A global bound over the original optical temperature supports combines
/// with the native conduction timestep bound. This is an explicit/lagged source,
/// including when used with an implicit conduction observer driver; refine time.
/// `endpoint_transport` instead evaluates the SAME exchange at a candidate
/// endpoint for the existing backward-Euler residual. Neither path supplies a
/// radiation entropy certificate or a new nonlinear thermal solver.
#[derive(Clone)]
pub struct GraySlabThermal<'a> {
    slab: AbsorbingSlab,
    layers: Arc<[BoundLayer<'a>]>,
    pairs: Arc<[(usize, usize, f64)]>,
    derivative_bounds_w_k: Arc<BTreeMap<CellId, f64>>,
    area_m2: f64,
    irradiation_w_m2: [f64; 2],
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
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DepositionError> {
        checkpoint()?;
        if !area_m2.is_finite()
            || area_m2 <= 0.0
            || cells.len() != slab.layer_count()
            || irradiation_w_m2.iter().any(|x| !x.is_finite() || *x < 0.0)
        {
            return Err(DepositionError::InvalidInput(
                "gray slab needs positive area, nonnegative irradiation and one original cell/law per layer",
            ));
        }
        let n = cells.len();
        let exterior_index_squared =
            slab.exterior_refractive_index() * slab.exterior_refractive_index();
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
            let volume = product(area_m2, slab.layers()[i].thickness_m())?;
            if cells[..i].iter().any(|&(previous, _)| previous == id)
                || !thermal.uses_material_law(cell.material_id, caloric)
                || (cell.geometry.volume_m3 / volume - 1.0).abs() > 64.0 * f64::EPSILON
            {
                return Err(DepositionError::InvalidInput(
                    "gray slab requires distinct cells, their original calorics and matching layer volumes",
                ));
            }
            let face = slab.face_coupling()[i];
            let face_area_m2 = [product(area_m2, face[0])?, product(area_m2, face[1])?];
            layers.push(BoundLayer {
                id,
                material: cell.material_id,
                segment: cell.segment_id,
                body: cell.body_id,
                mass_kg: cell.mass_kg,
                geometry: cell.geometry,
                caloric,
                face_area_m2,
            });
            let mut coupling = product(exterior_index_squared, face_area_m2[0] + face_area_m2[1])?;
            let mut correction = 0.0;
            for j in 0..n {
                let coefficient = product(area_m2, slab.layer_coupling()[i * n + j])?;
                add(&mut coupling, &mut correction, coefficient);
                if j > i && coefficient > 0.0 {
                    if pairs.len() == thermal.budget().max_links {
                        return Err(DepositionError::Budget("gray slab internal transfers"));
                    }
                    pairs.push((i, j, coefficient));
                }
            }
            // d(σ T^4)/dT is bounded over the SOURCE's full supported interval,
            // not merely at the current temperature or by a frozen secant.
            let high = slab.layers()[i].temperature_support_k()[1];
            let mut derivative = product(4.0 * SIGMA, coupling)?;
            for _ in 0..3 {
                derivative = product(derivative, high)?;
            }
            derivative_bounds_w_k.insert(id, derivative);
        }
        let result = Self {
            slab: slab.clone(),
            layers: layers.into(),
            pairs: pairs.into(),
            derivative_bounds_w_k: Arc::new(derivative_bounds_w_k),
            area_m2,
            irradiation_w_m2,
        };
        result.temperatures(thermal, &mut checkpoint)?;
        Ok(result)
    }

    pub fn slab(&self) -> &AbsorbingSlab {
        &self.slab
    }
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }
    pub const fn irradiation_w_m2(&self) -> [f64; 2] {
        self.irradiation_w_m2
    }

    /// Original native endpoint transport, ready for
    /// `BackwardEulerStep::evaluate_with_transport`. Physical sources and
    /// pairwise exchanges retain separate native energy-accounting owners.
    pub fn endpoint_transport(
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
            let [lo, hi] = self.slab.layers()[i].temperature_support_k();
            if !result[i].is_finite() || !(lo..=hi).contains(&result[i]) {
                return Err(DepositionError::MaterialDomain);
            }
        }
        Ok(result)
    }

    fn external_sources(
        &self,
        temperatures: &[f64],
    ) -> Result<BTreeMap<CellId, f64>, DepositionError> {
        let mut sources = BTreeMap::new();
        let exterior_index_squared =
            self.slab.exterior_refractive_index() * self.slab.exterior_refractive_index();
        for (layer, &temperature) in self.layers.iter().zip(temperatures) {
            if layer.face_area_m2 == [0.0; 2] {
                continue;
            }
            let square = product(temperature, temperature)?;
            let emission = product(product(SIGMA, square)?, square)?;
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
                        exterior_index_squared,
                    )
                    .map_err(|_| DepositionError::Arithmetic)?,
                );
            }
            if !power.is_finite() {
                return Err(DepositionError::Arithmetic);
            }
            if power != 0.0 {
                sources.insert(layer.id, power);
            }
        }
        Ok(sources)
    }

    fn internal_transfers(
        &self,
        temperatures: &[f64],
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<HeatTransfer>, DepositionError> {
        let mut transfers = Vec::with_capacity(self.pairs.len());
        for (k, &(i, j, coefficient)) in self.pairs.iter().enumerate() {
            if k % 256 == 0 {
                checkpoint()?;
            }
            let (a, b) = (temperatures[i], temperatures[j]);
            // Factored T^4 difference preserves nearby temperatures and gives
            // exactly zero at equilibrium before any large emission products.
            if a == b {
                continue;
            }
            let secant = product(product(product(coefficient, SIGMA)?, a + b)?, a * a + b * b)?;
            transfers.push(HeatTransfer {
                a: self.layers[i].id,
                b: self.layers[j].id,
                power_w: product(secant, a - b)?,
            });
        }
        Ok(transfers)
    }
}

impl DepositionObserver for GraySlabThermal<'_> {
    type Error = DepositionError;
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
