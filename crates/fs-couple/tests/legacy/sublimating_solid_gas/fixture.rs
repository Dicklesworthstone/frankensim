//! Analytic synthetic data through the ORIGINAL material/card/gas constructors.
use fs_couple::liquid_gas::evaporation::momentum::{
    MovingGasParcel, boundary::MovingGasBoundarySolve,
};
use fs_couple::solid_gas::kinetic::{
    KineticSublimationCandidate, KineticSublimationControls, KineticSublimationError,
    MovingSolidSurface, advance_solid_with_moving_gas,
};
use fs_material::sublimation::{
    ClausiusClapeyronSublimation,
    gas::{
        GasSurfaceBinding,
        chamber::GasChamberSurface,
        momentum::{MaxwellReflection, MolecularSurfaceFrame},
    },
};
use fs_material::thermochemical_gas::{
    ThermochemicalGasModel,
    energy::{
        GasEnergySolve,
        parcel::{ClosedGasParcel, transfer::species::boundary::GasSpeciesBoundarySolve},
    },
};
use fs_qty::Temperature;
use fs_thermal::sublimation::{SolidSurfaceState, SublimationError};
use fs_thermochem::{
    Composition, CompositionBasis, Nasa9RegionV1, Nasa9StandardStateModelV1, SpeciesId,
    UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R, mixture::FrozenIdealGasMixtureModelV1,
};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
#[path = "../../../../fs-material/examples/material_gas_spring/source.rs"]
mod source;

pub const LATENT: f64 = 200_000.0;
pub const MASS: f64 = 0.01;
pub const AREA: f64 = 1e-7;
pub const DT: f64 = 0.01;
pub const MOLAR_MASS: f64 = 0.018_015;

pub fn species() -> SpeciesId {
    SpeciesId::new("B").unwrap()
}
pub fn catalog(reference_j_kg: f64) -> ThermochemicalGasModel {
    let original = source::load(
        "A 0.99 0.02897 3.5 0 1.7e-5 0.025\nB 0.01 0.018015 4 0 2e-5 0.03",
        "synthetic finite sublimating recoil fixture",
    )
    .unwrap();
    if reference_j_kg == 0.0 {
        return original;
    }
    let species_models = original
        .thermodynamics()
        .components()
        .iter()
        .map(|component| {
            if component.species() != &species() {
                return component.clone();
            }
            let regions = component
                .regions()
                .iter()
                .map(|region| {
                    let mut card = region.card().clone();
                    card.parameters.get_mut("a7").unwrap().value += reference_j_kg * MOLAR_MASS / R;
                    Nasa9RegionV1::from_card(card).unwrap()
                })
                .collect();
            Nasa9StandardStateModelV1::new(
                component.species().clone(),
                component.molar_mass(),
                component.convention().clone(),
                regions,
            )
            .unwrap()
        })
        .collect();
    ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            species_models,
            original.thermodynamics().mole_composition().clone(),
        )
        .unwrap(),
        original.transport_materials().to_vec(),
        original.water_species(),
    )
    .unwrap()
}
pub fn binding(offset_j_kg: f64) -> GasSurfaceBinding {
    GasSurfaceBinding {
        species: species(),
        accommodation_coefficient: 0.1,
        solid_enthalpy_offset_j_kg: offset_j_kg,
    }
}
pub fn chart(catalog: &ThermochemicalGasModel, binding: &GasSurfaceBinding) -> (f64, f64) {
    let component = catalog
        .thermodynamics()
        .components()
        .iter()
        .find(|s| s.species() == &binding.species)
        .unwrap();
    let cp = component
        .evaluate(Temperature::new(400.0))
        .unwrap()
        .properties()
        .cp()
        .value()
        / component.molar_mass().value();
    let h = catalog
        .species_specific_enthalpy_j_kg(&binding.species, 400.0)
        .unwrap()
        - LATENT
        - binding.solid_enthalpy_offset_j_kg;
    (h, cp)
}
pub fn temperature(h: f64, chart: (f64, f64)) -> std::result::Result<f64, SublimationError> {
    let t = 400.0 + (h - chart.0) / chart.1;
    if t.is_finite() && (200.0..=900.0).contains(&t) {
        Ok(t)
    } else {
        Err(SublimationError::Invalid(
            "synthetic pure-solid caloric domain exceeded",
        ))
    }
}
pub fn enthalpy_range(chart: (f64, f64)) -> [f64; 2] {
    [chart.0 - 200.0 * chart.1, chart.0 + 500.0 * chart.1]
}
pub fn phase() -> ClausiusClapeyronSublimation {
    ClausiusClapeyronSublimation::try_new(400.0, 1000.0, LATENT, MOLAR_MASS, 917.0, [200.0, 900.0])
        .unwrap()
}
pub fn solid(t: f64, chart: (f64, f64), velocity: [f64; 3]) -> MovingSolidSurface {
    MovingSolidSurface::new(
        SolidSurfaceState {
            mass_kg: MASS,
            area_m2: AREA,
            specific_enthalpy_j_kg: chart.0 + chart.1 * (t - 400.0),
        },
        velocity,
    )
    .unwrap()
}
pub fn gas(
    catalog: &ThermochemicalGasModel,
    t: f64,
    vapor_pressure: f64,
    velocity: [f64; 3],
) -> MovingGasParcel {
    let fraction = vapor_pressure / 100_000.0;
    let model = catalog
        .with_composition(
            Composition::new(
                CompositionBasis::MoleFraction,
                vec![1.0 - fraction, fraction],
            )
            .unwrap(),
        )
        .unwrap();
    MovingGasParcel::new(
        ClosedGasParcel::new(model, source::point(t, 100_000.0).unwrap(), 0.01).unwrap(),
        velocity,
    )
    .unwrap()
}
pub fn frame() -> MolecularSurfaceFrame {
    MolecularSurfaceFrame {
        surface_id: 7,
        outward_normal: [1.0, 0.0, 0.0],
        centroid_m: [0.0; 3],
    }
}
pub fn controls() -> KineticSublimationControls {
    KineticSublimationControls {
        gas: MovingGasBoundarySolve {
            gas: GasSpeciesBoundarySolve {
                gas: GasEnergySolve {
                    minimum_temperature_k: 200.0,
                    maximum_temperature_k: 900.0,
                    energy_tolerance_j_kg: 1e-7,
                    maximum_iterations: 128,
                },
                mass_tolerance_kg: 1e-13,
                energy_tolerance_j: 5e-8,
            },
            momentum_absolute_n_s: 1e-12,
            momentum_relative_tolerance: 1e-10,
        },
        maximum_duration_s: 0.05,
        maximum_gross_inventory_fraction: 0.25,
        maximum_temperature_change_k: 5.0,
        latent_compatibility_tolerance_j_kg: 1e-5,
    }
}
#[allow(clippy::too_many_arguments)]
pub fn advance(
    solid: &MovingSolidSurface,
    gas: &MovingGasParcel,
    catalog: &ThermochemicalGasModel,
    binding: &GasSurfaceBinding,
    chart: (f64, f64),
    scattering: Option<&MaxwellReflection>,
    dt: f64,
    heat: f64,
    controls: KineticSublimationControls,
) -> std::result::Result<KineticSublimationCandidate, KineticSublimationError> {
    let law = phase();
    let caloric = |h| temperature(h, chart);
    let surface = GasChamberSurface {
        id: 7,
        law: &law,
        binding,
        solid: solid.solid(),
        temperature_from_enthalpy: &caloric,
    };
    advance_solid_with_moving_gas(
        &surface,
        solid,
        gas,
        catalog,
        frame(),
        enthalpy_range(chart),
        scattering,
        dt,
        heat,
        controls,
        || false,
    )
}
pub fn total_energy(solid: &MovingSolidSurface, gas: &MovingGasParcel, offset: f64) -> f64 {
    solid.solid().mass_kg * (solid.solid().specific_enthalpy_j_kg + offset)
        + solid.kinetic_energy_j()
        + gas.parcel().mass_kg() * gas.parcel().state().specific_internal_energy_j_kg()
        + gas.kinetic_energy_j()
}
pub fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e} (allowance {tolerance:.3e})"
    );
}
