//! Explicit source-caloric ingress for the existing CHT energy owner.

use super::{FluidSource, J, Result, bad, fields};
use fs_lbm::conjugate::AffineSpecificEnthalpy;
use fs_lbm::material::conjugate::buoyant::enthalpy::{
    ConvectionCaloricLaw, MaterialEnthalpyConfig, MaterialEnthalpyState,
    MaterialEnthalpyStepReport, MaterialEnthalpyTransientConfig,
};
use fs_matdb::{ClaimSelection, PropertyKey, SelectionPolicy};
use fs_material::MaterialError;
use fs_material::gas_species::input::temperature::IsobaricGasLaw;
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::{
    DENSITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS, SPECIFIC_HEAT_CAPACITY_PROPERTY,
};
use fs_material::thermal::transport::isobaric::{
    IsobaricHeatCapacityLaw, IsobaricHeatCapacityQuery,
};
use fs_qty::{Density, Temperature};

/// The source span and approximation/acceptance tolerances are scene inputs.
#[derive(Clone, Debug)]
pub(in crate::cht_command) struct Input {
    minimum_temperature_k: f64,
    maximum_temperature_k: f64,
    maximum_interpolation_error_j_kg: Option<f64>,
    maximum_flux_defect_w: f64,
    transient: Option<TransientControls>,
    /// Opt-in exact original Cp integration; the existing chart path remains unchanged.
    basis: Basis,
}

#[derive(Clone, Debug)]
enum Basis {
    Existing,
    IsobaricCp { pressure_axis: Option<String> },
}

#[derive(Clone, Copy, Debug)]
struct TransientControls {
    maximum_passes_per_step: usize,
    maximum_energy_solves: usize,
    temperature_tolerance_k: f64,
    maximum_storage_defect_w: f64,
    maximum_energy_closure_j: f64,
    temperature_relaxation: f64,
}

pub(in crate::cht_command) fn parse(
    value: Option<&J>,
    source: Option<&FluidSource>,
    is_transient: bool,
) -> Result<Option<Input>> {
    let Some(value) = value else {
        return Ok(None);
    };
    let source = source
        .filter(|source| source.temperature_dependent())
        .ok_or_else(|| bad("fluid_enthalpy requires a temperature-dependent fluid source"))?;
    fields(
        value,
        &[
            "minimum_temperature_k",
            "maximum_temperature_k",
            "maximum_interpolation_error_j_kg",
            "maximum_flux_defect_w",
            "transient_controls",
            "basis",
            "pressure_axis",
        ],
    )?;
    let basis = match value.get("basis") {
        None => {
            if value.get("pressure_axis").is_some() {
                return Err(bad("fluid_enthalpy.pressure_axis requires basis isobaric-cp"));
            }
            Basis::Existing
        }
        Some(basis) if basis.as_str() == Some("isobaric-cp") => {
            let pressure_axis = match source {
                FluidSource::Liquid { .. } => Some(super::string(value, "pressure_axis")?.into()),
                FluidSource::Gas { .. } => {
                    if value.get("pressure_axis").is_some() {
                        return Err(bad("gas isobaric-cp retains its original pressure; pressure_axis is only for liquid cards"));
                    }
                    None
                }
            };
            Basis::IsobaricCp { pressure_axis }
        }
        Some(_) => return Err(bad("fluid_enthalpy.basis must be isobaric-cp when specified")),
    };
    let number = |name| super::super::number(value, name, "fluid_enthalpy");
    let minimum_temperature_k = number("minimum_temperature_k")?;
    let maximum_temperature_k = number("maximum_temperature_k")?;
    let maximum_flux_defect_w = number("maximum_flux_defect_w")?;
    let maximum_interpolation_error_j_kg = if matches!(basis, Basis::IsobaricCp { .. }) {
        if value.get("maximum_interpolation_error_j_kg").is_some() {
            return Err(bad("isobaric-cp uses the exact original Cp integral; no interpolation-error field applies"));
        }
        None
    } else { match source {
        FluidSource::Liquid { .. } => Some(number("maximum_interpolation_error_j_kg")?),
        FluidSource::Gas { .. } => {
            if value.get("maximum_interpolation_error_j_kg").is_some() {
                return Err(bad(
                    "gas fluid_enthalpy uses its original analytic caloric law; no interpolation-error field applies",
                ));
            }
            None
        }
    }};
    let reference = source
        .evolving_reference_temperature()
        .expect("evolving source admitted");
    if minimum_temperature_k <= 0.0
        || maximum_temperature_k <= minimum_temperature_k
        || reference < minimum_temperature_k
        || reference > maximum_temperature_k
        || maximum_flux_defect_w < 0.0
        || maximum_interpolation_error_j_kg.is_some_and(|v| v <= 0.0)
    {
        return Err(bad(
            "fluid_enthalpy requires a positive increasing temperature span containing the fluid reference, nonnegative flux tolerance and positive liquid interpolation tolerance",
        ));
    }
    let transient = match (is_transient, value.get("transient_controls")) {
        (true, Some(value)) => Some(TransientControls::parse(value)?),
        (true, None) => {
            return Err(bad(
                "transient fluid_enthalpy requires explicit transient_controls",
            ));
        }
        (false, Some(_)) => {
            return Err(bad(
                "fluid_enthalpy.transient_controls requires a transient scene",
            ));
        }
        (false, None) => None,
    };
    Ok(Some(Input {
        minimum_temperature_k,
        maximum_temperature_k,
        maximum_interpolation_error_j_kg,
        maximum_flux_defect_w,
        transient,
        basis,
    }))
}

impl TransientControls {
    fn parse(value: &J) -> Result<Self> {
        fields(
            value,
            &[
                "maximum_passes_per_step",
                "maximum_energy_solves",
                "temperature_tolerance_k",
                "maximum_storage_defect_w",
                "maximum_energy_closure_j",
                "temperature_relaxation",
            ],
        )?;
        let number = |name| super::super::number(value, name, "fluid_enthalpy.transient_controls");
        let count = |name, maximum| -> Result<usize> {
            let value = number(name)?;
            if value < 1.0 || value > maximum as f64 || value.fract() != 0.0 {
                return Err(bad(format!(
                    "fluid_enthalpy.transient_controls.{name} must be an integer in 1..={maximum}"
                )));
            }
            Ok(value as usize)
        };
        let result = Self {
            maximum_passes_per_step: count("maximum_passes_per_step", 1_000_000)?,
            maximum_energy_solves: count("maximum_energy_solves", 10_000_000)?,
            temperature_tolerance_k: number("temperature_tolerance_k")?,
            maximum_storage_defect_w: number("maximum_storage_defect_w")?,
            maximum_energy_closure_j: number("maximum_energy_closure_j")?,
            temperature_relaxation: number("temperature_relaxation")?,
        };
        if result.temperature_tolerance_k < 0.0
            || result.maximum_storage_defect_w < 0.0
            || result.maximum_energy_closure_j < 0.0
            || result.temperature_relaxation <= 0.0
            || result.temperature_relaxation > 1.0
        {
            return Err(bad(
                "fluid_enthalpy transient tolerances must be nonnegative and temperature_relaxation in (0,1]",
            ));
        }
        Ok(result)
    }
}

enum Source<'a> {
    Gas(&'a IsobaricGasLaw),
    Liquid(ResolvedHeatCapacityCurve),
    IsobaricCp(IsobaricHeatCapacityLaw),
}

pub(super) struct Bound<'a> {
    source: Source<'a>,
    input: Input,
}

impl Input {
    pub(super) fn controls(&self, reference_temperature_k: f64) -> MaterialEnthalpyConfig {
        MaterialEnthalpyConfig {
            reference_temperature_k,
            maximum_flux_defect_w: self.maximum_flux_defect_w,
        }
    }

    pub(super) fn transient_controls(
        &self,
        reference_temperature_k: f64,
    ) -> Result<MaterialEnthalpyTransientConfig> {
        let controls = self
            .transient
            .ok_or_else(|| bad("missing fluid_enthalpy transient_controls"))?;
        Ok(MaterialEnthalpyTransientConfig {
            caloric: self.controls(reference_temperature_k),
            maximum_passes_per_step: controls.maximum_passes_per_step,
            maximum_energy_solves: controls.maximum_energy_solves,
            temperature_tolerance_k: controls.temperature_tolerance_k,
            maximum_storage_defect_w: controls.maximum_storage_defect_w,
            maximum_energy_closure_j: controls.maximum_energy_closure_j,
            temperature_relaxation: controls.temperature_relaxation,
        })
    }

    pub(super) fn bind<'a>(&self, fluid: &'a FluidSource) -> Result<Bound<'a>> {
        let source = match fluid {
            FluidSource::Gas {
                temperature_law: Some(law),
                ..
            } => {
                // Original transport/phase support remains authoritative at
                // the declared endpoints as well as every later trial.
                law.at_temperature(Temperature::new(self.minimum_temperature_k))
                    .and_then(|_| law.at_temperature(Temperature::new(self.maximum_temperature_k)))
                    .map_err(|error| bad(error.to_string()))?;
                if matches!(self.basis, Basis::IsobaricCp { .. }) {
                    law.admit_caloric_interval(
                        Temperature::new(self.minimum_temperature_k),
                        Temperature::new(self.maximum_temperature_k),
                    ).map_err(|error| bad(error.to_string()))?;
                }
                Source::Gas(law)
            }
            FluidSource::Liquid {
                state,
                pack,
                temperature_axis: Some(axis),
                ..
            } => {
                let point = state.material().query();
                let lower = point
                    .clone()
                    .with(axis, self.minimum_temperature_k)
                    .map_err(|error| bad(error.to_string()))?;
                let upper = point
                    .with(axis, self.maximum_temperature_k)
                    .map_err(|error| bad(error.to_string()))?;
                if let Basis::IsobaricCp { pressure_axis: Some(pressure_axis) } = &self.basis {
                    let cp = state.material().property(SPECIFIC_HEAT_CAPACITY_PROPERTY)
                        .ok_or_else(|| bad("isobaric-cp requires the original resolved fluid Cp claim"))?;
                    let selected = cp.answer().receipt.selected;
                    let claim = pack.card().claims().claim(selected)
                        .ok_or_else(|| bad("original fluid Cp claim is missing from its retained card"))?;
                    Source::IsobaricCp(IsobaricHeatCapacityLaw::try_new(
                        pack.card(),
                        IsobaricHeatCapacityQuery {
                            specific_heat_cp: &claim.key,
                            temperature_axis: axis,
                            pressure_axis,
                            lower: &lower,
                            upper: &upper,
                            selection: ClaimSelection::Pinned(selected),
                        },
                        fluid.evolving_reference_temperature().expect("evolving fluid admitted"),
                    ).map_err(|error| bad(error.to_string()))?)
                } else {
                Source::Liquid(
                    ResolvedHeatCapacityCurve::try_new(
                        pack.card(),
                        HeatCapacitySourceQuery {
                            heat_capacity: &PropertyKey::new(
                                SPECIFIC_HEAT_CAPACITY_PROPERTY,
                                SPECIFIC_HEAT_CAPACITY_DIMS,
                            ),
                            density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
                            temperature_axis: axis,
                            lower: &lower,
                            upper: &upper,
                            selection: SelectionPolicy::SingleClaimOnly,
                        },
                        SolidLiquidPhase::Liquid,
                        0.0,
                        self.maximum_interpolation_error_j_kg
                            .expect("liquid chart tolerance admitted"),
                    )
                    .map_err(|error| bad(error.to_string()))?,
                )
                }
            }
            _ => {
                return Err(bad(
                    "fluid_enthalpy needs its original evolving material source",
                ));
            }
        };
        Ok(Bound {
            source,
            input: self.clone(),
        })
    }
}

impl ConvectionCaloricLaw for Bound<'_> {
    fn tangent_at_temperature(
        &self,
        reference_temperature_k: f64,
        temperature_k: f64,
    ) -> std::result::Result<AffineSpecificEnthalpy, MaterialError> {
        for t in [reference_temperature_k, temperature_k] {
            if !t.is_finite()
                || t < self.input.minimum_temperature_k
                || t > self.input.maximum_temperature_k
            {
                return Err(MaterialError::State {
                    what: "CHT temperature is outside the declared fluid_enthalpy span".into(),
                });
            }
        }
        match &self.source {
            Source::Gas(law) => law.tangent_at_temperature(reference_temperature_k, temperature_k),
            Source::Liquid(law) => {
                law.tangent_at_temperature(reference_temperature_k, temperature_k)
            }
            Source::IsobaricCp(law) => {
                law.tangent_at_temperature(reference_temperature_k, temperature_k)
            }
        }
    }
}

impl Bound<'_> {
    pub(super) fn transient_report(
        &self,
        reference_temperature_k: f64,
        state: &MaterialEnthalpyState,
        records: &[MaterialEnthalpyStepReport],
    ) -> Result<String> {
        let controls = self.input.transient_controls(reference_temperature_k)?;
        let mut report = self.report(reference_temperature_k, state)?;
        report.pop(); // Extend the existing caloric evidence object.
        let number = super::super::num;
        let worst = |value: fn(&MaterialEnthalpyStepReport) -> f64| {
            records.iter().map(value).fold(0.0f64, f64::max)
        };
        report.push_str(&format!(concat!(
            ",\"storage\":\"source-enthalpy-backward-euler\",\"energy_solves\":{},",
            "\"maximum_passes_in_step\":{},\"maximum_passes_per_step\":{},\"maximum_energy_solves\":{},",
            "\"temperature_tolerance_k\":{},\"storage_tolerance_w\":{},\"energy_closure_tolerance_j\":{},",
            "\"temperature_relaxation\":{},\"maximum_temperature_defect_k\":{},\"maximum_storage_defect_w\":{},",
            "\"maximum_step_face_flux_defect_w\":{},\"maximum_step_cell_flux_defect_w\":{}}}"),
            records.iter().map(|r| r.energy_solves).sum::<usize>(),
            records.iter().map(|r| r.passes).max().unwrap_or(0),
            controls.maximum_passes_per_step, controls.maximum_energy_solves,
            number(controls.temperature_tolerance_k)?, number(controls.maximum_storage_defect_w)?,
            number(controls.maximum_energy_closure_j)?, number(controls.temperature_relaxation)?,
            number(worst(|r| r.temperature_defect_k))?, number(worst(|r| r.maximum_storage_defect_w))?,
            number(worst(|r| r.advection.maximum_face_flux_defect_w))?,
            number(worst(|r| r.advection.maximum_cell_flux_defect_w))?,
        ));
        Ok(report)
    }

    pub(super) fn report(
        &self,
        reference_temperature_k: f64,
        state: &MaterialEnthalpyState,
    ) -> Result<String> {
        let number = super::super::num;
        let (identity, exact) = match &self.source {
            Source::Gas(law) => (
                "null".into(), // Original gas receipts remain in material_inputs.
                if matches!(self.input.basis, Basis::IsobaricCp { .. }) {
                    format!(",\"basis\":\"isobaric-cp\",\"pressure_pa\":{}",
                        number(law.reference().state().pressure)?)
                } else { String::new() },
            ),
            Source::Liquid(law) => (super::super::quote(&law.curve().identity().to_hex()), String::new()),
            Source::IsobaricCp(law) => (
                "null".into(),
                format!(",\"basis\":\"isobaric-cp\",\"caloric_identity\":{},\"pressure_pa\":{}",
                    super::super::quote(&law.identity().to_hex()), number(law.pressure_pa())?),
            ),
        };
        let interpolation = self
            .input
            .maximum_interpolation_error_j_kg
            .map(number)
            .transpose()?
            .unwrap_or_else(|| "null".into());
        Ok(format!(
            concat!(
                "{{\"advection\":\"source-enthalpy-first-order-upwind\",",
                "\"reference_temperature_k\":{},\"minimum_temperature_k\":{},\"maximum_temperature_k\":{},",
                "\"maximum_interpolation_error_j_kg\":{},\"chart_identity\":{},",
                "\"flux_tolerance_w\":{},\"maximum_face_flux_defect_w\":{},\"maximum_cell_flux_defect_w\":{}{}}}"
            ),
            number(reference_temperature_k)?,
            number(self.input.minimum_temperature_k)?,
            number(self.input.maximum_temperature_k)?,
            interpolation,
            identity,
            number(self.input.maximum_flux_defect_w)?,
            number(state.report.maximum_face_flux_defect_w)?,
            number(state.report.maximum_cell_flux_defect_w)?,
            exact,
        ))
    }
}
