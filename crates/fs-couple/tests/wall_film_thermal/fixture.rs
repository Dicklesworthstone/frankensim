//! Manufactured original caloric sources for finite wall-film coupling tests.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_conduction::lumped::{LumpedThermalEnvironment, LumpedThermalTransport};
use fs_couple::{cavity_phs::piston::thermal::HeatExchangeSolve, finite_gas_thermal::wall_film::*};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::QueryPoint;
use fs_material::{
    evaporation::gas::{GasEvaporationSolve, GasLiquidBinding, wall_film::GasWallFilmInventory},
    phase::{
        EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase, UniformEnthalpyStep,
        UniformEnthalpyStepInput,
    },
    thermochemical_gas::{
        ThermochemicalGasModel,
        energy::{
            GasEnergySolve, parcel::ClosedGasParcel,
            parcel::transfer::species::boundary::GasSpeciesBoundarySolve,
        },
    },
};
use fs_qty::{
    QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermal::evaporation::{EvaporationError, EvaporationLaw, EvaporationProperties};
use fs_thermochem::{
    Composition, CompositionBasis, SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K,
    mixture::FrozenIdealGasMixtureModelV1,
};

#[allow(dead_code)]
#[path = "../material_exchange/gas.rs"]
mod gas_data;

pub const M: f64 = 0.01801528;
pub const R: f64 = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / M;
pub const CP_LIQUID: f64 = 4200.0;
pub const CP_WALL: f64 = 800.0;
pub const LATENT: f64 = 2.44e6;
pub const AREA: f64 = 0.01;
pub const GAS_VOLUME: f64 = 1e-4;

pub struct Phase;
pub static PHASE: Phase = Phase;
impl EvaporationLaw for Phase {
    fn temperature_range_k(&self) -> [f64; 2] {
        [280.0, 340.0]
    }
    fn properties(&self, t: f64) -> Result<EvaporationProperties, EvaporationError> {
        let dc = 3.5 * R - CP_LIQUID;
        let latent = LATENT + dc * (t - 300.0);
        let log_pressure = (LATENT - dc * 300.0) / R * (1.0 / 300.0 - 1.0 / t)
            + dc / R * fs_math::det::ln(t / 300.0);
        let pressure = 3536.8 * fs_math::det::exp(log_pressure);
        Ok(EvaporationProperties {
            saturation_pressure_pa: pressure,
            saturation_pressure_derivative_pa_k: pressure * latent / (R * t * t),
            latent_heat_j_kg: latent,
            liquid_density_kg_m3: 1000.0,
            molar_mass_kg_mol: M,
        })
    }
}
pub fn liquid_h(t: f64) -> Result<f64, EvaporationError> {
    if (280.0..=340.0).contains(&t) {
        Ok(CP_LIQUID * (t - 300.0))
    } else {
        Err(EvaporationError::Invalid(
            "original fixture liquid h(T) domain",
        ))
    }
}
pub fn liquid_t(h: f64) -> Result<f64, EvaporationError> {
    let t = 300.0 + h / CP_LIQUID;
    if (280.0..=340.0).contains(&t) {
        Ok(t)
    } else {
        Err(EvaporationError::Invalid(
            "original fixture liquid T(h) domain",
        ))
    }
}
pub fn water() -> SpeciesId {
    SpeciesId::new("water").unwrap()
}
pub fn carrier() -> SpeciesId {
    SpeciesId::new("carrier").unwrap()
}

pub fn source() -> ThermochemicalGasModel {
    let (carrier, carrier_transport) = gas_data::component("carrier", 0.028, 0.0);
    let (water, water_transport) = gas_data::component("water", M, -13.4e6 * M);
    ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![carrier, water],
            Composition::new(CompositionBasis::MoleFraction, vec![0.95, 0.05]).unwrap(),
        )
        .unwrap(),
        vec![carrier_transport, water_transport],
        Some(&self::water()),
    )
    .unwrap()
}
fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
pub fn gas(source: &ThermochemicalGasModel, t: f64, vapor_pressure: f64) -> ClosedGasParcel {
    let x = vapor_pressure / 100000.0;
    ClosedGasParcel::new(
        source
            .with_composition(
                Composition::new(CompositionBasis::MoleFraction, vec![1.0 - x, x]).unwrap(),
            )
            .unwrap(),
        QueryPoint::new()
            .with_quantity(
                "temperature",
                quantity(QuantityKind::AbsoluteTemperature),
                t,
            )
            .unwrap()
            .with_quantity("pressure", quantity(QuantityKind::Pressure), 100000.0)
            .unwrap()
            .with("wall-film-lot", 13.0)
            .unwrap(),
        GAS_VOLUME,
    )
    .unwrap()
}
pub fn binding(source: &ThermochemicalGasModel) -> GasLiquidBinding {
    GasLiquidBinding {
        species: water(),
        liquid_enthalpy_offset_j_kg: source
            .species_specific_enthalpy_j_kg(&water(), 300.0)
            .unwrap()
            - LATENT,
    }
}
pub fn curve() -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_single_phase(
        ContentHash([77; 32]),
        SolidLiquidPhase::Solid,
        [280.0, 340.0]
            .map(|t| EnthalpyPhaseKnot {
                temperature_k: t,
                specific_enthalpy_j_kg: CP_WALL * (t - 300.0),
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 1000.0,
            })
            .to_vec(),
    )
    .unwrap()
}
pub fn body(htc: f64) -> LumpedThermalEnvironment {
    LumpedThermalEnvironment {
        temperature_k: 999.0, // Must be replaced by actual finite gas.
        radiation_temperature_k: 280.0,
        convection_w_per_m2_k: htc,
        transport: LumpedThermalTransport::try_declared(100.0, 0.0).unwrap(),
        enthalpy_tolerance_j_kg: 1e-8,
        maximum_thermal_residual_j: 1e-9,
    }
}
pub fn input(
    curve: &EquilibriumEnthalpyPhaseCurve,
    t: f64,
    dt: f64,
) -> UniformEnthalpyStepInput<'_> {
    UniformEnthalpyStepInput {
        curve,
        initial: curve
            .state_at_specific_enthalpy(CP_WALL * (t - 300.0))
            .unwrap(),
        mass_kg: 0.002,
        volume_m3: 2e-6,
        surface_area_m2: AREA,
        internal_heat_j: 0.0,
        duration_s: dt,
    }
}
pub fn controls() -> WallFilmThermalSolve {
    WallFilmThermalSolve {
        phase: GasEvaporationSolve {
            gas: GasSpeciesBoundarySolve {
                gas: GasEnergySolve {
                    minimum_temperature_k: 280.0,
                    maximum_temperature_k: 340.0,
                    energy_tolerance_j_kg: 1e-8,
                    maximum_iterations: 128,
                },
                mass_tolerance_kg: 1e-14,
                energy_tolerance_j: 1e-10,
            },
            maximum_inventory_fraction: 0.2,
            maximum_temperature_increment_k: 10.0,
            latent_compatibility_tolerance_j_kg: 1e-5,
        },
        sensible: HeatExchangeSolve {
            heat_bracket_j: [-1.0, 1.0],
            source_tolerance_j: 1e-9,
            energy_tolerance_j: 1e-9,
            maximum_iterations: 100,
        },
        wall: HeatExchangeSolve {
            heat_bracket_j: [-4.0, 4.0],
            source_tolerance_j: 1e-9,
            energy_tolerance_j: 1e-8,
            maximum_iterations: 100,
        },
        temperature_tolerance_k: 1e-7,
    }
}
pub fn boundary(binding: &GasLiquidBinding) -> WallFilmThermalBoundary<'_> {
    WallFilmThermalBoundary {
        law: &PHASE,
        binding,
        area_m2: AREA,
        mass_transfer_coefficient_m_s: 0.02,
        liquid_enthalpy_at_temperature: &liquid_h,
        liquid_temperature_from_enthalpy: &liquid_t,
    }
}
pub fn environment<'a>(
    gas: &'a ClosedGasParcel,
    catalog: &'a ThermochemicalGasModel,
    film: GasWallFilmInventory,
    binding: &'a GasLiquidBinding,
    body: &'a LumpedThermalEnvironment,
) -> WallFilmThermalEnvironment<'a> {
    WallFilmThermalEnvironment {
        gas,
        catalog,
        film,
        boundary: boundary(binding),
        body,
        gas_heat_j: 0.0,
        solve: controls(),
    }
}
pub fn scope<R>(run: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 71,
                kernel_id: 0x57414c4c,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        run(&cx, &gate)
    })
}
pub fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tolerance {tolerance:e}"
    );
}
pub fn film_energy(film: GasWallFilmInventory, binding: &GasLiquidBinding) -> f64 {
    match film {
        GasWallFilmInventory::Dry => 0.0,
        GasWallFilmInventory::Wet {
            mass_kg,
            specific_enthalpy_j_kg,
        } => mass_kg * (specific_enthalpy_j_kg + binding.liquid_enthalpy_offset_j_kg),
    }
}
pub fn complete_inventory(
    input: UniformEnthalpyStepInput<'_>,
    gas: &ClosedGasParcel,
    film: GasWallFilmInventory,
    binding: &GasLiquidBinding,
    gas_heat_j: f64,
    next: &UniformEnthalpyStep<WallFilmThermalReport>,
) {
    let report = &next.report;
    let film_change = report.next_film.mass_kg() - film.mass_kg();
    near(
        film_change + report.next_gas.mass_kg() - gas.mass_kg(),
        0.0,
        1e-14,
    );
    near(
        film_change + report.next_gas.species_mass_kg(&water()).unwrap()
            - gas.species_mass_kg(&water()).unwrap(),
        0.0,
        1e-14,
    );
    near(
        report.next_gas.species_mass_kg(&carrier()).unwrap(),
        gas.species_mass_kg(&carrier()).unwrap(),
        1e-18,
    );
    let wall_energy = input.mass_kg
        * (next.state.specific_enthalpy_j_kg() - input.initial.specific_enthalpy_j_kg());
    let gas_energy = report.next_gas.mass_kg()
        * report.next_gas.state().specific_internal_energy_j_kg()
        - gas.mass_kg() * gas.state().specific_internal_energy_j_kg();
    let film_energy = film_energy(report.next_film, binding) - film_energy(film, binding);
    near(
        wall_energy + gas_energy + film_energy,
        input.internal_heat_j + gas_heat_j + report.radiation_heat_j,
        controls().wall.energy_tolerance_j,
    );
    near(
        wall_energy,
        input.internal_heat_j + next.external_heat_j,
        controls().wall.energy_tolerance_j,
    );
    assert!(report.sum_absolute_energy_residual_j <= controls().wall.energy_tolerance_j);
    assert_eq!(report.next_gas.volume_m3(), gas.volume_m3());
    assert_eq!(
        report.next_gas.point().axes().get("wall-film-lot"),
        Some(&13.0)
    );
}
