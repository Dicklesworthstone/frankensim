//! Source-explicit soot radiation through two isothermal layers between plates.
//!
//! cargo run -p fs-airflow --features compartment-fire --example soot_layer_radiation --
//!   n k rho_particle max_fv T_source_min T_source_max
//!   c_lower c_upper T_gas_lower T_gas_upper T_spectral_lower T_spectral_upper
//!   area length_lower length_upper T_wall_0 T_wall_1 emissivity_0 emissivity_1 source
//!
//! Units: densities/concentrations kg/m³, temperatures K, area m², lengths m;
//! n, k, maximum volume fraction and emissivities are dimensionless. Quote a
//! source label containing spaces. All 20 arguments are required. `k` is the
//! nonnegative absorptive part in `m=n-i k`. The optical source must support
//! the supplied density, constant indices, spectral-temperature range/loading,
//! dilute Rayleigh particles and negligible scattering; this example never
//! supplies a universal soot optical constant or predicts a soot yield.
//!
//! Geometry is the infinite-parallel-plate view-factor limit. Its area is a
//! reference patch area, while ONE authored normal chord traverses both layers.
//! This conservative angular approximation is not a hemispherical slab solution.
//! Emissivities are explicitly supplied diffuse-gray engineering inputs with
//! unstated uncertainty, bound through the original material-card owner.
//!
//! Spectral temperatures weight the Planck mean independently of actual layer
//! temperatures. Equal values choose local emission means. Other incident
//! spectra and finite-path spectral filtering remain a gray-model approximation.
//! Particles and gas share each layer's actual temperature. Reported layer heat
//! belongs to that thermally equilibrated mixture; no soot moles or pressure,
//! combustion heat, kinetics, deposition, or molecular opacity are inferred.
//!
//! An optional suffix enables a wavelength-resolved comparison on the same
//! geometry and actual temperatures:
//! `spectral lambda_min_m lambda_max_m heat_tol_W radiosity_tol_W_m2 max_evaluations max_omitted_fraction`.
//! The supplied Planck weighting temperatures then affect only the retained
//! gray comparison. Molecular absorption is explicitly zero. Gray boundary
//! emissivities and the normal-chord angular approximation still apply.

use fs_airflow::compartment::radiation::soot::{RayleighSoot, SootOpticalDomain};
use fs_airflow::compartment::radiation::spectral::{
    SootPopulation, SpectralBand, SpectralOptics, SpectralSolve,
};
use fs_airflow::compartment::radiation::{
    AbsorptionCoefficient, LayerAbsorption, OpticalPath, OpticalSegment, ParticipatingEnclosure,
    RadiationSolve,
};
use fs_conduction::radiation::{
    EMISSIVITY_DIMS, SURFACE_EMISSIVITY_PROPERTY, SurfaceEmissivity, ViewFactorMatrix,
};
use fs_evidence::ValidityDomain;
use fs_gas::multicomponent::compartment::Layer;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, SelectionPolicy, UncertaintyModel,
};
use fs_qty::{Area, Density, Length, Temperature};
use std::error::Error;

const USAGE: &str = "required arguments:\n\
  n k rho_particle_kg_m3 max_volume_fraction T_source_min_K T_source_max_K\n\
  soot_lower_kg_m3 soot_upper_kg_m3 T_gas_lower_K T_gas_upper_K\n\
  T_spectral_lower_K T_spectral_upper_K area_m2 path_lower_m path_upper_m\n\
  T_wall0_K T_wall1_K emissivity0 emissivity1 optical_source\n\
All optical inputs are source declarations; this is an estimated gray normal-chord slab.\n\
Illustrative smoke case (authored values, not a validated soot-property source):\n\
  1.75 1.03 1800 0.001 200 2500 0.00018 0.0018 600 1000 600 1000 1 0.25 0.75 300 400 0.8 0.6 authored-demo\n\
Set both soot concentrations to zero for the transparent comparison.\n\
Optional wavelength-resolved comparison suffix (all seven values together):\n\
  spectral lambda_min_m lambda_max_m heat_tol_W radiosity_tol_W_m2 max_evaluations max_omitted_fraction\n\
Example suffix: spectral 1e-7 1e-3 1e-5 1e-5 8192 1e-4\n\
This keeps the gray result for comparison, then resolves Rayleigh soot at each wavelength.\n\
Each wall/layer emits at its own actual temperature; the two T_spectral inputs only weight the gray comparison.\n\
Molecular absorption is declared transparent. Emissivities remain explicitly gray.\n\
The finite wavelength interval is a caller-declared optical-source assumption.\n\
Reported omitted blackbody emission is separate from quadrature error and is not a bound on omitted net heat.";

fn surface(
    value: f64,
    temperature: f64,
    source: &str,
    index: usize,
) -> Result<SurfaceEmissivity, Box<dyn Error>> {
    if !temperature.is_finite() || temperature <= 0.0 {
        return Err("wall temperatures must be finite and positive".into());
    }
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
        value: PropertyValue::Scalar {
            value,
            dims: EMISSIVITY_DIMS,
        },
        // A single queried state, without extrapolating the supplied input.
        validity: ValidityDomain::unconstrained().with("T", temperature, temperature),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: Vec::new(),
        provenance: Provenance {
            source: format!("caller-supplied diffuse gray surface {index}; scenario {source}"),
            license: "caller-supplied-example-input".into(),
            artifact: None,
        },
    })?;
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "caller-supplied opaque boundary".into(),
            phase: "solid".into(),
            process: format!("diffuse gray surface {index}"),
            revision: 0,
        },
        claims,
        Vec::new(),
    )?;
    Ok(SurfaceEmissivity::from_card(
        &format!("wall-{index}"),
        &card,
        temperature,
        SelectionPolicy::SingleClaimOnly,
    )?)
}

fn main() -> Result<(), Box<dyn Error>> {
    let arguments: Vec<String> = std::env::args().skip(1).collect();
    if arguments.len() == 1 && matches!(arguments[0].as_str(), "-h" | "--help") {
        println!("{USAGE}");
        return Ok(());
    }
    if arguments.len() != 20 && arguments.len() != 27 {
        return Err(USAGE.into());
    }
    let spectral = if arguments.len() == 27 {
        if arguments[20] != "spectral" {
            return Err(USAGE.into());
        }
        Some((
            [
                Length::new(arguments[21].parse::<f64>()?),
                Length::new(arguments[22].parse::<f64>()?),
            ],
            SpectralSolve {
                heat_tolerance_w: arguments[23].parse()?,
                radiosity_tolerance_w_m2: arguments[24].parse()?,
                maximum_evaluations: arguments[25].parse()?,
                maximum_omitted_blackbody_fraction: arguments[26].parse()?,
                ..SpectralSolve::default()
            },
        ))
    } else {
        None
    };
    let numeric = arguments[..19]
        .iter()
        .map(|s| s.parse::<f64>())
        .collect::<Result<Vec<_>, _>>()?;
    if numeric.iter().any(|v| !v.is_finite()) {
        return Err("all numeric inputs must be finite".into());
    }
    let [
        n,
        k,
        rho,
        max_fv,
        t_min,
        t_max,
        c_lower,
        c_upper,
        t_lower,
        t_upper,
        s_lower,
        s_upper,
        area,
        l_lower,
        l_upper,
        w0,
        w1,
        e0,
        e1,
    ]: [f64; 19] = numeric.try_into().map_err(|_| USAGE)?;
    let source = &arguments[19];
    let optical = RayleighSoot::new(
        n,
        k,
        Density::new(rho),
        SootOpticalDomain {
            spectral_temperature: [Temperature::new(t_min), Temperature::new(t_max)],
            maximum_volume_fraction: max_fv,
        },
        source,
    )?;
    let soot = [
        optical.planck_from_mass_concentration(Density::new(c_lower), Temperature::new(s_lower))?,
        optical.planck_from_mass_concentration(Density::new(c_upper), Temperature::new(s_upper))?,
    ];
    let absorption = [soot[0].layer_absorption()?, soot[1].layer_absorption()?];
    let enclosure = ParticipatingEnclosure::new(
        ViewFactorMatrix::infinite_parallel_plates(area)?,
        vec![OpticalPath {
            surfaces: [0, 1],
            exchange_area: Area::new(area),
            segments: vec![
                OpticalSegment {
                    layer: Layer::Lower,
                    length: Length::new(l_lower),
                },
                OpticalSegment {
                    layer: Layer::Upper,
                    length: Length::new(l_upper),
                },
            ],
        }],
        "authored reciprocal normal chord between infinite parallel plates, ordered lower then upper",
        1e-12,
        || false,
    )?;
    let gas = [Temperature::new(t_lower), Temperature::new(t_upper)];
    let walls = [Temperature::new(w0), Temperature::new(w1)];
    let emissivities = [surface(e0, w0, source, 0)?, surface(e1, w1, source, 1)?];
    let report = enclosure.evaluate(
        &walls,
        &emissivities,
        gas,
        &absorption,
        RadiationSolve::default(),
        &mut enclosure.workspace()?,
        || false,
    )?;

    println!(
        "model=constant-index Rayleigh soot, Planck-weighted gray, reciprocal normal-chord slab"
    );
    println!("optical_source={source}");
    println!("source_n={n:.12e}; source_absorptive_k={k:.12e}; particle_density_kg_m3={rho:.12e}");
    println!(
        "source_spectral_temperature_K=[{t_min:.12e},{t_max:.12e}]; source_maximum_volume_fraction={max_fv:.12e}"
    );
    println!(
        "all_heat_powers_are_positive_into_the_named_owner; geometry=normal_chord_angular_approximation"
    );
    println!(
        "layer,gas_temperature_K,spectral_temperature_K,soot_kg_m3,soot_volume_fraction,kappa_m_inv,path_m,optical_depth,heat_W"
    );
    for (index, name) in ["lower", "upper"].iter().enumerate() {
        let length = [l_lower, l_upper][index];
        let optical_depth = soot[index].coefficient.value() * length;
        println!(
            "{name},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{length:.12e},{optical_depth:.12e},{:.12e}",
            gas[index].value(),
            soot[index].spectral_temperature.value(),
            soot[index].mass_concentration.value(),
            soot[index].volume_fraction,
            soot[index].coefficient.value(),
            report.gas_heat_power[index].value()
        );
    }
    println!("surface,temperature_K,emissivity,area_m2,heat_W,outgoing_radiosity_W_m2");
    for index in 0..2 {
        println!(
            "wall{index},{:.12e},{:.12e},{area:.12e},{:.12e},{:.12e}",
            walls[index].value(),
            [e0, e1][index],
            report.surface_heat_power[index].value(),
            report.surface_radiosity_w_m2[index]
        );
    }
    println!(
        "surface_equation_max_residual_W_m2={:.12e}",
        report.maximum_surface_residual_w_m2
    );
    println!(
        "enclosure_energy_residual_W={:.12e}",
        report.energy_residual_w
    );
    println!(
        "model_scope=estimated absorption; no soot yield, molecular absorption, scattering, or gas pressure contribution"
    );
    if let Some((wavelengths, controls)) = spectral {
        let molecular = [
            LayerAbsorption::new(
                AbsorptionCoefficient::new(0.0),
                "authored transparent lower-layer molecular background; soot-only comparison",
            )?,
            LayerAbsorption::new(
                AbsorptionCoefficient::new(0.0),
                "authored transparent upper-layer molecular background; soot-only comparison",
            )?,
        ];
        let optics = SpectralOptics::new(
            vec![SpectralBand {
                wavelengths,
                absorption: molecular,
            }],
            soot.map(|s| {
                Some(SootPopulation {
                    model: &optical,
                    volume_fraction: s.volume_fraction,
                })
            }),
            format!(
                "caller-supplied diffuse gray engineering wall inputs, held constant over the declared wavelength window; {source}"
            ),
        )?;
        let resolved = enclosure.evaluate_spectral(
            &walls,
            &emissivities,
            gas,
            &optics,
            controls,
            &mut enclosure.spectral_workspace()?,
            || false,
        )?;
        println!(
            "spectral_model=constant-index Rayleigh soot with individual Planck spectra; same normal-chord geometry and gray boundaries"
        );
        println!(
            "spectral_window_m=[{:.12e},{:.12e}]",
            wavelengths[0].value(),
            wavelengths[1].value()
        );
        println!(
            "spectral_heat_tolerance_W={:.12e}; radiosity_tolerance_W_m2={:.12e}; maximum_evaluations={}; maximum_omitted_blackbody_fraction={:.12e}",
            controls.heat_tolerance_w,
            controls.radiosity_tolerance_w_m2,
            controls.maximum_evaluations,
            controls.maximum_omitted_blackbody_fraction
        );
        println!(
            "spectral_wavelength_evaluations={}; factorizations={}; accepted_panels={}",
            resolved.wavelength_evaluations, resolved.factorizations, resolved.accepted_panels
        );
        println!("layer,spectral_heat_W,gray_comparison_heat_W,quadrature_error_estimate_W");
        for (i, name) in ["lower", "upper"].iter().enumerate() {
            println!(
                "{name},{:.12e},{:.12e},{:.12e}",
                resolved.gas_heat_power[i].value(),
                report.gas_heat_power[i].value(),
                resolved.gas_heat_error_estimate_w[i]
            );
        }
        println!(
            "surface,spectral_heat_W,gray_comparison_heat_W,quadrature_error_estimate_W,spectral_radiosity_W_m2,radiosity_error_estimate_W_m2"
        );
        for i in 0..2 {
            println!(
                "wall{i},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e}",
                resolved.surface_heat_power[i].value(),
                report.surface_heat_power[i].value(),
                resolved.surface_heat_error_estimate_w[i],
                resolved.surface_radiosity_w_m2[i],
                resolved.radiosity_error_estimate_w_m2[i]
            );
        }
        println!(
            "node,temperature_K,participating,omitted_shortwave_fraction,omitted_longwave_fraction,omitted_blackbody_emission_W_m2"
        );
        for (name, coverage) in ["wall0", "wall1", "lower", "upper"]
            .iter()
            .zip(&resolved.coverage)
        {
            println!(
                "{name},{:.12e},{},{:.12e},{:.12e},{:.12e}",
                coverage.temperature.value(),
                coverage.participating,
                coverage.shortwave_fraction,
                coverage.longwave_fraction,
                coverage.omitted_blackbody_w_m2
            );
        }
        println!(
            "spectral_surface_residual_integral_W_m2={:.12e}; enclosure_energy_residual_W={:.12e}",
            resolved.surface_residual_integral_w_m2, resolved.energy_residual_w
        );
        println!(
            "spectral_scope=estimated finite-window quadrature; supplied optical constants; no molecular spectroscopy, soot-yield prediction or scattering; omitted blackbody emission is not omitted net heat"
        );
    }
    Ok(())
}
