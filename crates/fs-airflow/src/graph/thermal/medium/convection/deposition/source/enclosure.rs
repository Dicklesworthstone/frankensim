//! Enclosed horizontal, inclined or vertical air layers on reciprocal links.
//!
//! Both endpoints are finite isothermal WALLS. The intervening gas layer has
//! negligible storage in this quasi-steady approximation; it is not a third
//! thermal cell or a prescribed-temperature reservoir. The original link owns
//! both energy transfers, area coverage and the outgoing-conductance bound.

/// Gas gaps in complete native heat, radiation and circuit assemblies.
pub mod heat;
/// Reevaluate the original material at each accepted explicit thermal substep.
pub mod observer;

use super::super::CellConvectionEvaluation;
use super::*;

/// One complete rectangular air gap between two parallel isothermal walls.
///
/// `lower_cell` declares the physically lower wall independently of the link's
/// endpoint ordering. At exactly 90 degrees it instead names the reference wall
/// for signed heat; either vertical wall may be hotter. `gap` is the normal wall
/// separation; `plate_spans` are the two WHOLE lateral spans, ordered as
/// [width, vertical height] at 90 degrees. They are not inferred from cell centers,
/// exposed-area remainders or contact patches. The native link must cover the
/// full rectangle on BOTH sides, and its authored conductance must be zero: the selected
/// Nusselt relation already includes the gas's Fourier conduction exactly once.
/// `inclination_degrees` is the actual wall-plane angle from horizontal, not
/// the angle between gravity and the wall plane. Zero retains the original
/// signed horizontal law; positive angles through 60 degrees select the full
/// inclined heated-from-below law, including its along-slope correction.
/// Exactly 90 degrees selects the vertical-cavity law with height/gap in
/// [5,110] and width/gap above 10. Intermediate angles above 60 degrees refuse.
/// The caller owns this physical orientation; mesh coordinates do not infer it.
///
/// Mean-temperature source properties are a Boussinesq approximation. The
/// caller explicitly bounds `|beta * (T_lower - T_upper)|` below one. In
/// addition, the original correlation admits only air, continuum transport,
/// its stated aspect-ratio and Rayleigh ranges and no
/// imposed throughflow or moving-wall shear. Matching Pr alone does not
/// establish that another fluid is supported by the air experiments.
#[derive(Clone, Copy, Debug)]
pub struct SourceEnclosedGasLayer {
    pub link: u64,
    pub lower_cell: CellId,
    pub gap: Length,
    pub plate_spans: [Length; 2],
    /// Nonnegative full gravity magnitude, not its wall-normal projection.
    pub gravity_m_s2: f64,
    /// Wall-plane angle from horizontal in [0, 60] or exactly 90 degrees.
    /// Positive angles through 60 refuse negative signed Rayleigh; exactly
    /// 90 uses its magnitude. Rebind after changing physical orientation.
    pub inclination_degrees: f64,
    /// Explicit relative-density-departure bound in [0, 1).
    pub maximum_relative_density_departure: f64,
}

/// Actual wall temperatures, retained original gas source and total gap heat.
#[derive(Debug)]
pub struct EnclosedGasLayerSample<F> {
    pub link: u64,
    pub lower_cell: CellId,
    pub upper_cell: CellId,
    pub lower_temperature_k: f64,
    pub upper_temperature_k: f64,
    pub property_temperature_k: f64,
    pub gap_m: f64,
    pub inclination_degrees: f64,
    pub area_m2: f64,
    pub conductance_w_k: f64,
    /// Signed lower-to-upper heat. The original link applies its exact opposite
    /// to the lower wall. This is internal exchange, never reservoir heat.
    pub heat_into_upper_w: f64,
    pub material: F,
    pub thermal: ThermalFluidTransportState,
    /// Includes Nu=1 at equilibrium and for stable horizontal stratification.
    /// A finite gap therefore retains its conductive coupling when delta T is zero.
    pub evaluation: CellConvectionEvaluation,
}

/// Step-local total conductances with original property and correlation owners.
#[derive(Debug)]
pub struct ResolvedEnclosedGasLayers<F> {
    conductances_w_k: BTreeMap<u64, f64>,
    samples: Vec<EnclosedGasLayerSample<F>>,
}
impl<F> ResolvedEnclosedGasLayers<F> {
    pub fn conductances_w_k(&self) -> &BTreeMap<u64, f64> {
        &self.conductances_w_k
    }
    pub fn samples(&self) -> &[EnclosedGasLayerSample<F>] {
        &self.samples
    }

    /// Install once in an original implicit trial's native conductance map.
    /// All ownership conflicts and cancellation are checked before mutation.
    pub fn apply_to(
        &self,
        destination: &mut BTreeMap<u64, f64>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), ConvectionSourceError> {
        check(&mut cancelled)?;
        for link in self.conductances_w_k.keys() {
            check(&mut cancelled)?;
            if destination.contains_key(link) {
                return Err(ConvectionSourceError::Invalid(
                    "enclosed gas layer already has a conductance owner",
                ));
            }
        }
        check(&mut cancelled)?;
        destination.extend(self.conductances_w_k.iter().map(|(&id, &g)| (id, g)));
        Ok(())
    }
}

/// Evaluate gas-layer transport between the ACTUAL two caloric wall states.
///
/// Use the actual proposed state inside `BackwardEulerStep` or a coupled heat
/// residual. Alternatively, the accompanying native observer resolves the
/// source at every explicit substep. The original thermal integrator alone
/// transfers `G * (T_lower - T_upper)` reciprocally and includes `G` in its CFL.
/// Do not also emit two source rates or add a second gas Fourier conductance.
/// Independent radiation can retain its own heat owner.
///
/// The resolver is queried once per layer in stable link-ID order, at the mean
/// of the original wall temperatures. It must preserve one original gas law,
/// composition and declared pressure model across all trials; the returned
/// owner is retained with the result. A prepared gas law may compute pressure
/// consistently from prescribed mass and fixed gap volume if its negligible
/// energy-storage approximation is justified. This adapter neither invents gas
/// mass/heat capacity nor advances pressure, composition or momentum.
///
/// Endpoint caloric temperatures must represent uniform isothermal wall-face
/// temperatures. Nonuniform wall conduction and gas thermal transients need
/// their own resolved owners; cell-center temperatures do not automatically
/// satisfy that approximation. Sidewall thermal bridges and radiation are
/// outside this gas relation. Every source/domain refusal or cancellation
/// returns no partial map and leaves the original thermal history untouched.
pub fn evaluate_source_enclosed_gas_layers<F: ThermalFluidTransport, E>(
    thermal: &DepositionSimulation<'_>,
    bindings: &[SourceEnclosedGasLayer],
    maximum_layers: usize,
    mut resolve: impl FnMut(
        &SourceEnclosedGasLayer,
        &DepositionSimulation<'_>,
        Temperature,
    ) -> Result<F, E>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ResolvedEnclosedGasLayers<F>, SourceCellConvectionError<E>> {
    check(&mut cancelled)?;
    if bindings.len() > maximum_layers || bindings.len() > thermal.links().len() {
        return Err(ConvectionSourceError::Invalid(
            "enclosed gas layers exceed their explicit work bound",
        )
        .into());
    }
    let mut ids = BTreeSet::new();
    let mut ordered = Vec::with_capacity(bindings.len());
    for binding in bindings {
        check(&mut cancelled)?;
        let link = thermal
            .links()
            .get(&binding.link)
            .ok_or(ConvectionSourceError::Invalid(
                "inactive enclosed gas-layer link",
            ))?;
        if !ids.insert(binding.link) {
            return Err(ConvectionSourceError::Invalid("duplicate enclosed gas layer").into());
        }
        if binding.lower_cell != link.a && binding.lower_cell != link.b {
            return Err(ConvectionSourceError::Invalid(
                "lower wall must be an original link endpoint",
            )
            .into());
        }
        if link.conductance_w_k != 0.0 {
            return Err(ConvectionSourceError::Invalid(
                "enclosed gas layers require zero authored conductance",
            )
            .into());
        }
        let gap = binding.gap.value();
        let [width, length] = binding.plate_spans.map(|span| span.value());
        if [gap, width, length]
            .into_iter()
            .any(|x| !x.is_finite() || x <= 0.0)
            || !binding.gravity_m_s2.is_finite()
            || binding.gravity_m_s2 < 0.0
            || !binding.inclination_degrees.is_finite()
            || (!(0.0..=60.0).contains(&binding.inclination_degrees)
                && binding.inclination_degrees != 90.0)
            || !binding.maximum_relative_density_departure.is_finite()
            || !(0.0..1.0).contains(&binding.maximum_relative_density_departure)
        {
            return Err(ConvectionSourceError::Invalid(
                "gas layers require positive dimensions, nonnegative gravity, inclination in [0,60] or exactly 90 degrees and a density bound in [0,1)",
            ).into());
        }
        let area = finite(width * length)?;
        if area <= 0.0
            || !link.area_m2.is_finite()
            || link.area_m2 <= 0.0
            || (link.area_m2 - area).abs() > 32.0 * f64::EPSILON * area.max(link.area_m2)
        {
            return Err(ConvectionSourceError::Invalid(
                "gas-layer native area must cover the complete declared plate rectangle",
            )
            .into());
        }
        let patch = thermal
            .contact_patch(binding.link)
            .ok_or(ConvectionSourceError::Invalid(
                "inactive enclosed gas-layer coverage",
            ))?;
        if patch.covered_area_m2.into_iter().any(|covered| {
            !covered.is_finite()
                || covered <= 0.0
                || (covered - area).abs() > 32.0 * f64::EPSILON * area.max(covered)
        }) {
            return Err(ConvectionSourceError::Invalid(
                "gas-layer coverage must occupy the complete rectangle on both walls",
            )
            .into());
        }
        if binding.inclination_degrees == 90.0 {
            let width_over_gap = finite(width / gap)?;
            let height_over_gap = finite(length / gap)?;
            if width_over_gap <= 10.0 || !(5.0..=110.0).contains(&height_over_gap) {
                return Err(ConvectionSourceError::Invalid(
                    "vertical gas layers require whole width/gap above 10 and whole height/gap in [5,110]",
                )
                .into());
            }
        } else {
            let aspect = finite(gap / width.min(length))?;
            if aspect <= 0.0 || aspect >= 0.1 {
                return Err(ConvectionSourceError::Invalid(
                    "gas layers require gap divided by each whole span below 0.1",
                )
                .into());
            }
        }
        ordered.push(binding);
    }
    ordered.sort_by_key(|binding| binding.link);
    let mut conductances_w_k = BTreeMap::new();
    let mut samples = Vec::with_capacity(ordered.len());
    for binding in ordered {
        check(&mut cancelled)?;
        let link = &thermal.links()[&binding.link];
        let lower_cell = binding.lower_cell;
        let upper_cell = if lower_cell == link.a { link.b } else { link.a };
        let lower_temperature_k = thermal.temperature_k(lower_cell)?;
        let upper_temperature_k = thermal.temperature_k(upper_cell)?;
        if [lower_temperature_k, upper_temperature_k]
            .into_iter()
            .any(|t| !t.is_finite() || t <= 0.0)
        {
            return Err(ConvectionSourceError::Invalid(
                "enclosed gas layers require positive finite absolute wall temperatures",
            )
            .into());
        }
        let property_temperature_k = 0.5 * lower_temperature_k + 0.5 * upper_temperature_k;
        let primary = material_at(lower_cell, property_temperature_k, || {
            resolve(binding, thermal, Temperature::new(property_temperature_k))
        })?;
        let expansion =
            primary
                .thermal
                .volumetric_expansion_per_k()
                .ok_or(ConvectionSourceError::Invalid(
                    "enclosed gas layers require same-source volumetric expansion",
                ))?;
        let difference = finite(lower_temperature_k - upper_temperature_k)?;
        let density_departure = finite(expansion * difference)?.abs();
        if density_departure > binding.maximum_relative_density_departure {
            return Err(ConvectionSourceError::Invalid(
                "enclosed gas layer exceeds its declared relative density departure",
            )
            .into());
        }
        let transport = primary.thermal.transport();
        let signed_rayleigh = transport
            .rayleigh(
                binding.gravity_m_s2,
                expansion,
                difference,
                binding.gap.value(),
            )
            .map_err(|error| ConvectionSourceError::Material(error.into()))?;
        let minimum_span = binding.plate_spans[0]
            .value()
            .min(binding.plate_spans[1].value());
        let aspect = binding.gap.value() / minimum_span;
        let (correlation, inputs) = if binding.inclination_degrees == 90.0 {
            (
                CorrelationId::WrightVerticalAirLayer,
                CorrelationInputs::vertical_air_layer(
                    signed_rayleigh.abs(),
                    transport.prandtl(),
                    binding.plate_spans[1].value() / binding.gap.value(),
                ),
            )
        } else if binding.inclination_degrees == 0.0 {
            (
                CorrelationId::HollandsHorizontalAirLayer,
                CorrelationInputs::horizontal_air_layer(
                    signed_rayleigh,
                    transport.prandtl(),
                    aspect,
                ),
            )
        } else {
            (
                CorrelationId::HollandsInclinedAirLayer,
                CorrelationInputs::inclined_air_layer(
                    signed_rayleigh,
                    transport.prandtl(),
                    aspect,
                    binding.inclination_degrees,
                ),
            )
        };
        let nusselt = fs_convection::evaluate(correlation, inputs)
            .map_err(|error| ConvectionSourceError::Material(error.into()))?;
        let coefficient = nusselt
            .heat_transfer_coefficient(
                ThermalConductivity::new(transport.thermal_conductivity_w_m_k()),
                binding.gap,
            )
            .map_err(|error| ConvectionSourceError::Material(error.into()))?;
        let conductance_w_k = finite(link.area_m2 * coefficient.value.value())?;
        let heat_into_upper_w = finite(conductance_w_k * difference)?;
        if conductance_w_k <= 0.0 || (difference != 0.0 && heat_into_upper_w == 0.0) {
            return Err(ConvectionSourceError::Invalid(
                "enclosed gas-layer conductance or nonzero heat underflowed to zero",
            )
            .into());
        }
        conductances_w_k.insert(binding.link, conductance_w_k);
        samples.push(EnclosedGasLayerSample {
            link: binding.link,
            lower_cell,
            upper_cell,
            lower_temperature_k,
            upper_temperature_k,
            property_temperature_k,
            gap_m: binding.gap.value(),
            inclination_degrees: binding.inclination_degrees,
            area_m2: link.area_m2,
            conductance_w_k,
            heat_into_upper_w,
            material: primary.material,
            thermal: primary.thermal,
            evaluation: CellConvectionEvaluation {
                nusselt,
                coefficient,
            },
        });
    }
    check(&mut cancelled)?;
    Ok(ResolvedEnclosedGasLayers {
        conductances_w_k,
        samples,
    })
}
