use super::*;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, Cx, ExecMode, StreamKey};
use fs_lbm::d3q19::freesurface3::ContactModel3;
use fs_lbm::material::free_surface::enthalpy::FreeSurfaceEnthalpyControls;
use fs_lbm::material::free_surface::{
    MaterialFreeSurfaceDomain3, MaterialFreeSurfaceInitialCell3, plan_free_surface_from_material,
};
use fs_matdb::*;
use fs_material::liquid::surface::{LiquidCapillarity, resolve_liquid_surface_state};
use fs_material::liquid_temperature::{LiquidTemperatureLaw, LiquidViscosityBasis};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use fs_material::state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec, Temperature};
use fs_thermal::evaporation::EvaporationProperties;
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::{Composition, CompositionBasis, SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K};

#[allow(dead_code)]
#[path = "../material_exchange/gas.rs"]
mod gas_data;

pub const M: f64 = 0.01801528;
pub const R: f64 = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / M;
pub const CP: f64 = 4200.0;
pub const LATENT: f64 = 2.44e6;
pub const DT: f64 = 0.001;
pub const INTERFACE: [usize; 3] = [3, 3, 3];

/// Manufactured water-like phase law with the SAME original gas/liquid caloric
/// gap. It tests coupling; these authored coefficients are not empirical water.
pub struct Phase;
impl EvaporationLaw for Phase {
    fn temperature_range_k(&self) -> [f64; 2] {
        [280.0, 340.0]
    }
    fn properties(&self, t: f64) -> Result<EvaporationProperties, EvaporationError> {
        let dc = 3.5 * R - CP;
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

pub fn water() -> SpeciesId {
    SpeciesId::new("water").unwrap()
}
pub fn carrier() -> SpeciesId {
    SpeciesId::new("carrier").unwrap()
}
pub fn source(reference: f64) -> ThermochemicalGasModel {
    let (carrier, carrier_transport) = gas_data::component("carrier", 0.028, 0.0);
    let (water, water_transport) = gas_data::component("water", M, reference * M);
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
            .with("specimen-lot", 27.0)
            .unwrap(),
        1e-8,
    )
    .unwrap()
}
pub fn binding(source: &ThermochemicalGasModel, datum: f64) -> GasLiquidBinding {
    GasLiquidBinding {
        species: water(),
        liquid_enthalpy_offset_j_kg: source
            .species_specific_enthalpy_j_kg(&water(), 300.0)
            .unwrap()
            - datum
            - LATENT,
    }
}
pub fn surface<'a>(
    law: &'a dyn EvaporationLaw,
    binding: &'a GasLiquidBinding,
) -> LiquidGasEvaporationSurface<'a> {
    LiquidGasEvaporationSurface {
        liquid: INTERFACE,
        gas: [4, 3, 3],
        area_m2: 1e-6,
        mass_transfer_coefficient_m_s: 0.02,
        law,
        binding,
    }
}
pub fn controls() -> EvaporatingLiquidGasSolve {
    EvaporatingLiquidGasSolve {
        phase: GasEvaporationSolve {
            gas: GasSpeciesBoundarySolve {
                gas: GasEnergySolve {
                    minimum_temperature_k: 280.0,
                    maximum_temperature_k: 340.0,
                    energy_tolerance_j_kg: 1e-6,
                    maximum_iterations: 128,
                },
                mass_tolerance_kg: 1e-18,
                energy_tolerance_j: 1e-11,
            },
            maximum_inventory_fraction: 0.2,
            maximum_temperature_increment_k: 1.0,
            latent_compatibility_tolerance_j_kg: 1e-4,
        },
        native: FreeSurfaceMassExchangeControls {
            mass_tolerance_kg: 1e-18,
            energy_tolerance_j: 1e-11,
            momentum_tolerance_kg_m_s: 1e-18,
        },
        mass_tolerance_kg: 1e-15,
        energy_tolerance_j: 1e-9,
    }
}

fn liquid_card() -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, values) in [
        ("density", Dims([-3, 1, 0, 0, 0, 0]), [1000.0; 2]),
        (
            "specific_heat_capacity",
            Dims([2, 0, -2, -1, 0, 0]),
            [CP; 2],
        ),
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, [0.5; 2]),
        ("dynamic_viscosity", Dims([-1, 1, -1, 0, 0, 0]), [0.1, 0.02]),
    ] {
        let key = if name == "thermal_conductivity" {
            PropertyKey::with_quantity(name, quantity(QuantityKind::ThermalConductivity))
        } else {
            PropertyKey::new(name, dims)
        };
        claims
            .insert_claim(PropertyClaim {
                key,
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(280.0, values[0]), (340.0, values[1])],
                    dims,
                },
                validity: ValidityDomain::unconstrained().with("T", 280.0, 340.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured finite-phase native liquid".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured liquid".into(),
            phase: "liquid".into(),
            process: "authored".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

pub fn with_liquid<R>(
    datum: f64,
    fill: f64,
    run: impl FnOnce(&MaterialFreeSurfaceEnthalpy3<'_, '_>, &CancelGate) -> R,
) -> R {
    let card = liquid_card();
    let point = QueryPoint::new().with("T", 300.0).unwrap();
    let law = LiquidTemperatureLaw::try_new(
        &card,
        &point,
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let source = resolve_liquid_surface_state(
        &card,
        &point,
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let caloric = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        card.content_hash(),
        SolidLiquidPhase::Liquid,
        [280.0, 340.0]
            .map(|temperature_k| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: datum + CP * (temperature_k - 300.0),
                temperature_k,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: 1000.0,
            })
            .to_vec(),
    )
    .unwrap();
    let mut cells = Vec::new();
    for z in 1..=6 {
        for y in 1..=6 {
            for x in 0..=3 {
                let coordinates = [x, y, z];
                cells.push(MaterialFreeSurfaceInitialCell3 {
                    coordinates,
                    fill_fraction: if coordinates == INTERFACE { fill } else { 1.0 },
                    velocity_m_s: [0.0, 0.003 * (x as f64 - 1.5), 0.0],
                });
            }
        }
    }
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 71,
                kernel_id: 29,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let flow = plan_free_surface_from_material(&source, 0.001, DT, 0.1)
            .unwrap()
            .free_surface3_from_cells(
                MaterialFreeSurfaceDomain3 {
                    cells: [8; 3],
                    gravity_m_s2: [0.0; 3],
                    contact: ContactModel3::Neutral,
                    solid_tiles: &[],
                    max_cells: 512,
                },
                &cells,
                &cx,
            )
            .unwrap();
        let liquid = flow
            .with_enthalpy(
                &law,
                &caloric,
                FreeSurfaceEnthalpyControls::default(),
                |_| datum,
                &gate,
            )
            .unwrap();
        run(&liquid, &gate)
    })
}
