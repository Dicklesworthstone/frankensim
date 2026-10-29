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

use fs_airflow::compartment::radiation::soot::{RayleighSoot, SootOpticalDomain};
use fs_airflow::compartment::radiation::{
    OpticalPath, OpticalSegment, ParticipatingEnclosure, RadiationSolve,
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
Set both soot concentrations to zero for the transparent comparison.";

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
    if arguments.len() != 20 {
        return Err(USAGE.into());
    }
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
    Ok(())
}
