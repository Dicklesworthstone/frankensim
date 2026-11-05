//! Source-owned gas films in series with the original solid half-cell.
//!
//! The original mode holds conductivity at the actual cell phase. An explicit
//! pure-phase profile instead integrates the same pinned source k(T) through
//! the solid path. Both transfer the original gas film's actual signed heat;
//! the film owner retains every correlation, property convention and refusal.
//! Series heat balance: https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node118.html
use super::*;
use crate::air_path::input::convection::{
    GasConvectionBoundary, GasConvectionError, GasConvectionSample,
};
use fs_material::gas_species::input::{GasRecord, temperature::IsobaricGasLaw};
use fs_material::phase::SolidLiquidPhase;
use fs_material::state_point::{
    ScalarAdmissibility, ScalarPropertyRequirement, THERMAL_CONDUCTIVITY_DIMS,
    resolve_material_state_point,
};
use fs_material::thermal::conduction::{
    ConvectiveLayerError, LayerConductionConfig, ReservoirLayerGeometry,
};
use fs_qty::{
    QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_matdb::{ClaimSelection, PropertyLineSlice};

/// Conductivity through the existing solid half-layer. The integrated model
/// is homogeneous, source-free and quasi-steady; its original thermal cell
/// still owns all storage. It never infers a phase boundary inside the layer.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SolidGasConductivityProfile {
    /// Original constant-in-space coefficient at the current cell phase.
    /// A zero source conductivity explicitly represents insulation.
    CellPhase,
    /// Continuously supported positive k(T) on one explicitly selected pure
    /// phase of the SAME caloric curve. Only actual material temperatures must
    /// remain in this span; the gas bulk may lie outside it. Preparation pins
    /// the reference-selected source once and retains every original knot.
    IntegratedPurePhase {
        phase: SolidLiquidPhase,
        temperature_range_k: [f64; 2],
    },
}

/// Reproduce and pin the original bulk before varying only film temperature.
/// Shared by fixed-volume and moving-gas owners; a bare snapshot is insufficient.
pub(in crate::solid_enthalpy) fn endpoint_law(
    model: &fs_material::thermochemical_gas::ThermochemicalGasModel,
    point: &fs_matdb::QueryPoint,
    state: &fs_material::thermochemical_gas::ResolvedThermochemicalGasState,
) -> Result<IsobaricGasLaw, SolidGasError> {
    let law = IsobaricGasLaw::try_from_model(model.clone(), point.clone())
        .map_err(GasConvectionError::Material)?;
    if !matches!(law.reference(), GasRecord::Thermochemical(reference) if reference == state) {
        return Err(SolidGasError::Invalid(
            "film source law must reproduce the exact finite-gas endpoint",
        ));
    }
    Ok(law)
}

/// Prescribed exterior convection on one finite gas/solid patch.
///
/// Geometry/regime/speed/gravity stay prescribed; pressure, composition and
/// gas properties come from EACH actual finite-gas endpoint. This is a local
/// well-mixed bulk approximation, not a resolved boundary layer or flow solver.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SolidGasConvection {
    pub boundary: GasConvectionBoundary,
    pub solid_profile: SolidGasConductivityProfile,
    /// Bounded interior surface trials, 1..=128. Each polls cancellation.
    /// The integrated profile permits at most two additional film calls for
    /// its supported bracket; the held-k mode permits one limiting-endpoint call.
    pub maximum_surface_iterations: usize,
    /// Absolute |Q_solid - Q_film| allowance [W]. Zero requests exact balance.
    /// This is separate from the outer cell/gas joule residual and time error.
    pub surface_power_tolerance_w: f64,
}

impl SolidGasConvection {
    pub(super) fn validate(self, geometry: GasConductionGeometry) -> Result<(), SolidGasError> {
        self.boundary.validate()?;
        if self.boundary.exposed_area_m2() != geometry.area_m2
            || geometry.gas_distance_m != 0.0
            || !geometry.cell_distance_m.is_finite()
            || geometry.cell_distance_m < 0.0
            || !(1..=128).contains(&self.maximum_surface_iterations)
            || !self.surface_power_tolerance_w.is_finite()
            || self.surface_power_tolerance_w < 0.0
        {
            return Err(SolidGasError::Invalid(
                "convection needs the same patch area, no gas conduction layer, a nonnegative solid distance and bounded surface controls",
            ));
        }
        if let SolidGasConductivityProfile::IntegratedPurePhase { phase, temperature_range_k } =
            self.solid_profile
        {
            let [lo, hi] = temperature_range_k;
            if phase == SolidLiquidPhase::SolidLiquid
                || !lo.is_finite() || lo <= 0.0 || !hi.is_finite() || hi <= lo
                || geometry.cell_distance_m <= 0.0
            {
                return Err(SolidGasError::Invalid(
                    "integrated convection needs a positive solid path and an explicit pure-phase temperature span",
                ));
            }
        }
        Ok(())
    }
}

/// Exact accepted interface query; cloning shares immutable correlation evidence.
#[derive(Clone, Debug)]
pub struct SolidGasConvectionSample {
    /// Includes the ACTUAL surface temperature, endpoint bulk temperature,
    /// film/bulk/surface gas source queries, coefficient and correlation evidence.
    pub film: Arc<GasConvectionSample>,
    /// Fourier half-layer power, positive into the gas. The integrated profile
    /// retains the original inverse-integral flux coordinate, avoiding loss of
    /// thin-wall accuracy from reintegrating a rounded surface temperature.
    /// With zero distance the isothermal patch carries film power directly.
    pub solid_layer_power_w: f64,
    /// Actual solid-layer power minus actual film power [W]; independently gated.
    pub surface_balance_residual_w: f64,
    /// Instantaneous half-layer entropy production [W/K], separate from the
    /// film's own production. Neither is a time-integrated entropy audit.
    pub solid_layer_entropy_generation_w_k: f64,
    pub surface_iterations: usize,
    /// Original compiled source line, shared across endpoint evaluations and
    /// wall substeps. None denotes the original held cell-phase conductivity.
    pub solid_conductivity_profile: Option<Arc<PropertyLineSlice>>,
}

impl<'c, 'm> PreparedContact<'c, 'm> {
    /// Compile after reference source selection, before any nonlinear endpoint
    /// or explicit wall-substep trial. Trial fluxes borrow this original chart.
    pub(super) fn prepare_solid_profile(
        &self,
        source: &StructuralSource<'m>,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Option<Arc<ResolvedConductionLayer<'m>>>, SolidGasError> {
        let Some(film) = self.contact.convection else { return Ok(None) };
        let SolidGasConductivityProfile::IntegratedPurePhase { phase, temperature_range_k } =
            film.solid_profile else { return Ok(None) };
        film.validate(self.contact.geometry)?;
        let MaterialPropertySelection::PinnedByProperty(pins) = &self.selection else {
            return Err(SolidGasError::Invalid("integrated solid profile needs its reference claim pin"));
        };
        if pins.len() != 1 || pins[0].0 != self.contact.conductivity.name() {
            return Err(SolidGasError::Invalid("integrated solid profile has an incompatible conductivity pin"));
        }
        let layer = ResolvedConductionLayer::new(
            source.curve, source.point, &self.contact.conductivity,
            ClaimSelection::Pinned(pins[0].1), phase, temperature_range_k, checkpoint,
        )?;
        Ok(Some(Arc::new(layer)))
    }

    // Preserve the same source/quantity/phase/pin admission as the original
    // layer owner, without inventing a gas distance to access a material query.
    pub(super) fn convection_material(
        &self,
        source: &StructuralSource<'_>,
        phase: EquilibriumPhaseState,
        transport: Option<&MaterialNetworkSample>,
    ) -> Result<ResolvedMaterialStatePoint, SolidGasError> {
        let key = &self.contact.conductivity;
        if key.quantity() != QuantitySpec::dimensional(THERMAL_CONDUCTIVITY_DIMS)
            && key.quantity()
                != QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::ThermalConductivity,
                    ValueForm::Static,
                ))
        {
            return Err(SolidGasError::Invalid(
                "static scalar solid conductivity required",
            ));
        }
        let point = source.point.at(phase)?;
        let material = match transport::cell_material(transport, self.contact.cell)? {
            Some(material) => material.clone(),
            None => {
                let requirement =
                    ScalarPropertyRequirement::try_with_key(key, ScalarAdmissibility::NonNegative)
                        .map_err(ThermalMaterialError::Material)?;
                resolve_material_state_point(
                    source.point.material,
                    &point,
                    &[requirement],
                    self.selection.clone(),
                )
                .map_err(ThermalMaterialError::Material)?
            }
        };
        if material.card_identity() != source.point.material.content_hash()
            || material.axis_quantities() != point.axis_quantities()
            || material.query_point().len() != point.axes().len()
            || material
                .query_point()
                .iter()
                .zip(point.axes())
                .any(|((a, x), (b, y))| a != b || x.to_bits() != y.to_bits())
        {
            return Err(SolidGasError::Invalid(
                "convective solid snapshot must retain its exact cell phase and source conditions",
            ));
        }
        let property = material.property_by_key(key).ok_or(SolidGasError::Invalid(
            "convective solid conductivity is absent",
        ))?;
        if let MaterialPropertySelection::PinnedByProperty(pins) = &self.selection {
            if pins.len() != 1
                || pins[0].0 != key.name()
                || pins[0].1 != property.answer().receipt.selected
            {
                return Err(SolidGasError::Invalid(
                    "convective solid conductivity lost its reference claim pin",
                ));
            }
        }
        source.point.admit_static(property)?;
        if !property.value_si().is_finite() || property.value_si() < 0.0 {
            return Err(SolidGasError::Invalid(
                "solid conductivity must be finite and nonnegative",
            ));
        }
        Ok(material)
    }

    pub(super) fn resolve_convection(
        &self,
        source: &StructuralSource<'_>,
        law: &IsobaricGasLaw,
        phase: EquilibriumPhaseState,
        transport: Option<&MaterialNetworkSample>,
        model: SolidGasConvection,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<ContactSample, SolidGasError> {
        checkpoint()?;
        model.validate(self.contact.geometry)?;
        let material = self.convection_material(source, phase, transport)?;
        if matches!(model.solid_profile, SolidGasConductivityProfile::IntegratedPurePhase { .. }) {
            return self.resolve_integrated_convection(source, law, phase, material, model, checkpoint);
        }
        let conductivity = material
            .property_by_key(&self.contact.conductivity)
            .ok_or(SolidGasError::Invalid(
                "convective solid conductivity is absent",
            ))?
            .value_si();
        let cell_temperature = phase.temperature_k();
        let gas_temperature = law.reference().state().temperature;
        let mut evaluate = |temperature| -> Result<GasConvectionSample, SolidGasError> {
            checkpoint()?;
            let film = model.boundary.evaluate(law, temperature)?;
            checkpoint()?;
            Ok(film)
        };
        let geometry = self.contact.geometry;
        if geometry.cell_distance_m == 0.0 {
            // Explicit uniform/isothermal patch approximation: the cell is the
            // surface; there is no unreported half-layer temperature drop.
            let film = evaluate(cell_temperature)?;
            return finish(material, film, cell_temperature, None, 0);
        }
        if conductivity == 0.0 {
            // Insulation belongs to the actual solid source. Its gas-side face
            // equilibrates with the gas; a natural-law zero is not used as k=0.
            let film = evaluate(gas_temperature)?;
            return finish(material, film, cell_temperature, Some(0.0), 0);
        }
        let solid_g = finite((geometry.area_m2 / geometry.cell_distance_m) * conductivity)?;
        if solid_g <= 0.0 {
            return Err(SolidGasError::Invalid(
                "solid half-layer conductance is unrepresentable",
            ));
        }
        let first = evaluate(cell_temperature)?;
        if cell_temperature == gas_temperature {
            return finish(material, first, cell_temperature, Some(solid_g), 0);
        }
        let difference = cell_temperature - gas_temperature;
        // The held-coefficient series solution is a useful first trial, not a
        // frozen-film answer. Every accepted query still closes its actual flux.
        let film_g = first.conductance_w_k;
        let mut fraction = if solid_g >= film_g {
            1.0 / (1.0 + film_g / solid_g)
        } else {
            let ratio = solid_g / film_g;
            ratio / (1.0 + ratio)
        };
        if fraction <= 0.0 || fraction >= 1.0 {
            fraction = 0.5;
        }
        let (mut low, mut high) = (0.0, 1.0);
        // At the gas endpoint the film heat is exactly zero; this is a flux
        // value, not an assumed equilibrium coefficient or derivative. The
        // cell endpoint was already evaluated. An overflowing hypothetical
        // half-layer endpoint merely disables interpolation until narrowed;
        // every evaluated trial and accepted power remains finite.
        let mut low_residual = solid_g * difference;
        let mut high_residual = first.heat_into_wall_w;
        let mut interpolated_width = None;
        let mut residual = high_residual;
        for iteration in 1..=model.maximum_surface_iterations {
            let temperature =
                (1.0 - fraction).mul_add(gas_temperature, fraction * cell_temperature);
            let film = evaluate(temperature)?;
            let solid_power = finite(solid_g * (cell_temperature - temperature))?;
            residual = finite(solid_power + film.heat_into_wall_w)?;
            if residual.abs() <= model.surface_power_tolerance_w {
                return finish(material, film, cell_temperature, Some(solid_g), iteration);
            }
            if residual.is_sign_positive() == difference.is_sign_positive() {
                low = fraction;
                low_residual = residual;
            } else {
                high = fraction;
                high_residual = residual;
            }
            let width = high - low;
            let midpoint = low + 0.5 * width;
            // Safeguarded secant: reuse actual bracket residuals without an
            // extra source evaluation or numerical derivative. If the last
            // interpolation did not halve its bracket, take a midpoint; thus
            // the bracket contracts by at least half within two more trials.
            let next = if low_residual.is_finite()
                && high_residual.is_finite()
                && interpolated_width.is_none_or(|previous| width <= 0.5 * previous)
            {
                let a = low_residual.abs();
                let b = high_residual.abs();
                let weight = if a >= b {
                    1.0 / (1.0 + b / a)
                } else {
                    let ratio = a / b;
                    ratio / (1.0 + ratio)
                };
                let secant = low + width * weight;
                if secant > low && secant < high {
                    interpolated_width = Some(width);
                    secant
                } else {
                    interpolated_width = None;
                    midpoint
                }
            } else {
                interpolated_width = None;
                midpoint
            };
            if next <= low || next >= high || next == fraction {
                return Err(SolidGasError::SurfaceBalance {
                    residual_w: residual,
                    iterations: iteration,
                });
            }
            fraction = next;
        }
        Err(SolidGasError::SurfaceBalance {
            residual_w: residual,
            iterations: model.maximum_surface_iterations,
        })
    }

    fn resolve_integrated_convection(
        &self,
        source: &StructuralSource<'_>,
        law: &IsobaricGasLaw,
        phase: EquilibriumPhaseState,
        material: ResolvedMaterialStatePoint,
        model: SolidGasConvection,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<ContactSample, SolidGasError> {
        let layer = self.integrated_solid.as_ref().ok_or(SolidGasError::Invalid(
            "integrated convection needs its prepared source profile",
        ))?;
        if !std::ptr::eq(layer.caloric_curve(), source.curve) || layer.phase() != phase.phase() {
            return Err(SolidGasError::Invalid(
                "integrated solid profile must retain the original caloric object and pure phase",
            ));
        }
        // This also excludes a phase plateau or a cell outside the prepared
        // source support even when its temperature happens to be unchanged.
        layer.temperature_at_enthalpy(phase.specific_enthalpy_j_kg())?;
        let geometry = ReservoirLayerGeometry::planar(
            self.contact.geometry.area_m2, self.contact.geometry.cell_distance_m,
        )?;
        let result = geometry.evaluate_convective_material_to_bath(
            layer.conductivity_curve(),
            [phase.temperature_k(), law.reference().state().temperature],
            LayerConductionConfig {
                max_iterations: model.maximum_surface_iterations,
                absolute_heat_tolerance_w: model.surface_power_tolerance_w,
                relative_heat_tolerance: 0.0,
            },
            |surface| {
                let film = model.boundary.evaluate(law, surface)?;
                Ok::<_, GasConvectionError>((-film.heat_into_wall_w, film.conductance_w_k, film))
            },
            checkpoint,
        ).map_err(|error| match error {
            ConvectiveLayerError::Thermal(error) => SolidGasError::from(error),
            ConvectiveLayerError::Film(error) => SolidGasError::from(error),
            ConvectiveLayerError::SurfaceBalance { residual_w, iterations } =>
                SolidGasError::SurfaceBalance { residual_w, iterations },
        })?;
        finish_balance(
            material, result.film, phase.temperature_k(), result.solid_power_to_bath_w,
            result.conductance_w_k, result.iterations, Some(Arc::clone(layer.source_line())),
        )
    }
}

fn finish(
    source: ResolvedMaterialStatePoint,
    film: GasConvectionSample,
    cell_temperature: f64,
    solid_g: Option<f64>,
    iterations: usize,
) -> Result<ContactSample, SolidGasError> {
    let power_to_gas_w = -film.heat_into_wall_w;
    let solid_power = match solid_g {
        Some(g) => finite(g * (cell_temperature - film.wall_temperature_k))?,
        None => power_to_gas_w,
    };
    let film_g = film.conductance_w_k;
    if cell_temperature != film.reservoir_temperature_k && solid_g != Some(0.0) && film_g == 0.0 {
        return Err(SolidGasError::Invalid(
            "an equilibrium film sample cannot insulate a driven contact",
        ));
    }
    // Held-coefficient series value avoids dividing by a nearly zero bulk
    // temperature difference. An absent natural equilibrium coefficient still
    // declares only zero heat, never a temperature derivative.
    let conductance_w_k = match solid_g {
        None => film_g,
        Some(g) if g == 0.0 || film_g == 0.0 => 0.0,
        Some(g) if g >= film_g => finite(film_g / (1.0 + film_g / g))?,
        Some(g) => finite(g / (1.0 + g / film_g))?,
    };
    if solid_g != Some(0.0) && film_g > 0.0 && conductance_w_k <= 0.0 {
        return Err(SolidGasError::Invalid(
            "series convection conductance underflows",
        ));
    }
    finish_balance(source, film, cell_temperature, solid_power, conductance_w_k, iterations, None)
}

fn finish_balance(
    source: ResolvedMaterialStatePoint,
    film: GasConvectionSample,
    cell_temperature: f64,
    solid_power: f64,
    conductance_w_k: f64,
    iterations: usize,
    solid_conductivity_profile: Option<Arc<PropertyLineSlice>>,
) -> Result<ContactSample, SolidGasError> {
    let power_to_gas_w = -film.heat_into_wall_w;
    let residual = finite(solid_power - power_to_gas_w)?;
    let entropy = finite(
        (solid_power / film.wall_temperature_k)
            * ((cell_temperature - film.wall_temperature_k) / cell_temperature),
    )?;
    if entropy < 0.0
        || (solid_power != 0.0 && cell_temperature != film.wall_temperature_k && entropy == 0.0)
    {
        return Err(SolidGasError::Invalid(
            "solid half-layer entropy production is negative or underflows",
        ));
    }
    Ok(ContactSample {
        conductance_w_k,
        power_to_gas_w,
        source,
        convection: Some(SolidGasConvectionSample {
            film: Arc::new(film),
            solid_layer_power_w: solid_power,
            surface_balance_residual_w: residual,
            solid_layer_entropy_generation_w_k: entropy,
            surface_iterations: iterations,
            solid_conductivity_profile,
        }),
    })
}
