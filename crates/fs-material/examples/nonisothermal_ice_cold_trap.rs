//! Two ice films exchange mass AND molecular energy with one evolving vapor parcel.
//!
//! Run: cargo run -p fs-material --example nonisothermal_ice_cold_trap
//! Pass 0.1 as an optional argument to tighten local time-error targets tenfold.
//! Adaptive step doubling accepts the conservative two-half-step history; its
//! local estimate is not a bound on global or experimental error.
//!
//! Source/domain assumptions:
//! * The original WaterIceSublimation owner supplies IAPWS ice pressure and
//!   ideal-vapor Clausius-Clapeyron latent heat (solid volume is neglected).
//! * The original NASA-9 evaluator supplies water-vapor h, u and cp. Coefficients
//!   below are the 200--1000 K H2O region in NASA/TP-2002-211556, Appendix D,
//!   printed p.131; also data/matdb/seed-v1/water-vapor-nasa9-mcbride2002/nasa9.tsv:
//!   <https://ntrs.nasa.gov/api/citations/20020085330/downloads/20020085330.pdf>.
//! * This example DERIVES an ice caloric approximation h_s(T)=h_v(T)-L(T)-datum,
//!   so cp_s=cp_v-dL/dT. It is thermodynamically compatible with these two owners;
//!   it is NOT an independently measured ice heat-capacity correlation or the
//!   full IAPWS ice Gibbs equation. The bounded 250--270 K chart has positive cp.
//! * Pure-vapor transport uses the dilute factors of IAPWS R12-08 Eq.(11) and
//!   R15-11 Eq.(16), sampled every kelvin. This is an explicit low-density
//!   continuation below the triple point, limited here to 250--270 K. These
//!   transport properties satisfy the common gas-material contract; the molecular
//!   exchange below does not add a continuum conductivity/film heat flux.
//!
//! This is a fixed-area, well-mixed local Maxwellian model. Each film is a
//! uniform-temperature lump; no spatial gas transport, wall heat, radiation,
//! melting, nucleation, frost porosity or geometric growth is resolved. Density
//! 917 kg/m3 and accommodation 1 are explicit scenario choices, not fitted data.
//! The original gas carries both streams' h_v(T)-R_specific*T/2; there is no
//! separate latent source to double count. All accepted solid and gas states
//! are published together. This reduced example is not a chamber-design claim.

use std::collections::BTreeMap;
use std::error::Error;

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy, LawId, LawParameter,
    MaterialCard, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    QueryPoint, UncertaintyModel,
};
use fs_material::{
    state_point::{
        MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY,
    },
    sublimation::{
        WaterIceSublimation,
        gas::{
            GasSurfaceBinding, GasSurfaceSolve,
            chamber::{
                GasChamberSolve, GasChamberSurface,
                advance::{GasChamberAdvanceSolve, advance_chamber_adaptive},
            },
        },
    },
    thermochemical_gas::{
        SpeciesTransportMaterial, ThermochemicalGasModel,
        energy::{
            GasEnergySolve,
            parcel::{ClosedGasParcel, transfer::species::boundary::GasSpeciesBoundarySolve},
        },
    },
};
use fs_qty::{
    DynViscosity, MolarMass, Pressure, QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermal::sublimation::{SolidSurfaceState, SublimationError, SublimationLaw};
use fs_thermochem::{
    Composition, CompositionBasis, ElementalReferenceIdV1, NASA9_LAW_ID_V1, NASA9_LAW_VERSION_V1,
    NASA9_STATE_SCHEMA_VERSION_V1, Nasa9RegionV1, Nasa9StandardStateModelV1,
    ReferenceEquationOfStateV1, SpeciesId, StandardStateConventionV1, StandardStatePhaseV1,
    mixture::FrozenIdealGasMixtureModelV1, required_nasa9_parameter_dims,
};

const MIN_T: f64 = 250.0;
const MAX_T: f64 = 270.0;
const DATUM_T: f64 = 256.0;
const NASA_SOURCE: &str =
    "https://ntrs.nasa.gov/api/citations/20020085330/downloads/20020085330.pdf";
const VISCOSITY_SOURCE: &str = "https://iapws.org/technical-guidance/release/viscosity.download";
const CONDUCTIVITY_SOURCE: &str = "https://iapws.org/technical-guidance/release/ThCond.download";

fn typed(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

fn water_vapor(species: &SpeciesId) -> Result<Nasa9StandardStateModelV1, Box<dyn Error>> {
    // Source coefficients only: evaluate all thermodynamics in the NASA-9 owner.
    let coefficients = [
        -3.947960830e4,
        5.755731020e2,
        9.317826530e-1,
        7.222712860e-3,
        -7.342557370e-6,
        4.955043490e-9,
        -1.336933246e-12,
        -3.303974310e4,
        1.724205775e1,
    ];
    let mut parameters = BTreeMap::new();
    for ((name, dims), value) in required_nasa9_parameter_dims()
        .into_iter()
        .take(9)
        .zip(coefficients)
    {
        parameters.insert(name.to_owned(), LawParameter { value, dims });
    }
    parameters.insert(
        "reference_pressure".into(),
        LawParameter {
            value: 100000.0,
            dims: Pressure::DIMS,
        },
    );
    let card = ConstitutiveModelCard {
        law: LawId(NASA9_LAW_ID_V1.into()),
        law_version: NASA9_LAW_VERSION_V1,
        parameters,
        state_schema_version: NASA9_STATE_SCHEMA_VERSION_V1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 200.0, 1000.0),
        sources: vec![],
        provenance: Provenance {
            source: format!("{NASA_SOURCE}, Appendix D, H2O gas, 200--1000 K"),
            license: "Work-of-the-US-Government-Public-Use-Permitted".into(),
            artifact: None,
        },
    };
    Ok(Nasa9StandardStateModelV1::new(
        species.clone(),
        MolarMass::new(0.01801528),
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(100000.0),
            ElementalReferenceIdV1::new("NASA-TP-2002-211556-reference-elements-298.15K-1bar")?,
        )?,
        vec![Nasa9RegionV1::from_card(card)?],
    )?)
}

fn dilute_transport(temperature: f64) -> [f64; 2] {
    // Published zero-density factors; Horner polynomials in inverse T/Tcrit.
    let reduced = temperature / 647.096;
    let inverse = 1.0 / reduced;
    let viscosity_denominator =
        ((-0.241605 * inverse + 0.6366564) * inverse + 2.20462) * inverse + 1.67752;
    let conductivity_denominator =
        (((0.0004096266 * inverse - 0.003454586) * inverse + 0.006770357) * inverse + 0.01323095)
            * inverse
            + 0.002443221;
    [
        1e-6 * 100.0 * reduced.sqrt() / viscosity_denominator,
        1e-3 * reduced.sqrt() / conductivity_denominator,
    ]
}

fn transport(species: &SpeciesId) -> Result<SpeciesTransportMaterial, Box<dyn Error>> {
    let mut claims = ClaimSet::new();
    for (column, name, dims, source) in [
        (0, "dynamic_viscosity", DynViscosity::DIMS, VISCOSITY_SOURCE),
        (
            1,
            THERMAL_CONDUCTIVITY_PROPERTY,
            THERMAL_CONDUCTIVITY_DIMS,
            CONDUCTIVITY_SOURCE,
        ),
    ] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Curve {
                abscissa: "temperature".into(),
                abscissa_dims: Temperature::DIMS,
                knots: (0..=20)
                    .map(|i| {
                        let temperature = MIN_T + f64::from(i);
                        (temperature, dilute_transport(temperature)[column])
                    })
                    .collect(),
                dims,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "temperature",
                    typed(QuantityKind::AbsoluteTemperature),
                    MIN_T,
                    MAX_T,
                )
                // Author-selected dilute-gas application rectangle; no alteration
                // or expansion of a checked-in source-pack validity envelope.
                .with_quantity("pressure", typed(QuantityKind::Pressure), 1e-3, 1000.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: format!("{source}, zero-density term, 1 K linear samples, 250--270 K"),
                license: "IAPWS-reproduction-with-attribution".into(),
                artifact: None,
            },
        })?;
    }
    Ok(SpeciesTransportMaterial {
        species: species.clone(),
        material: MaterialCard::assemble(
            MaterialStateId {
                chemistry: "H2O".into(),
                phase: "gas".into(),
                process: "example-IAPWS-dilute-vapor".into(),
                revision: 0,
            },
            claims,
            vec![],
        )?,
        selection: MaterialPropertySelection::SingleClaimOnly,
    })
}

/// Example-local compatible caloric approximation, evaluated from original owners.
struct IceCaloric<'a> {
    law: &'a WaterIceSublimation,
    vapor: &'a ThermochemicalGasModel,
    species: &'a SpeciesId,
    datum: f64,
}
impl IceCaloric<'_> {
    fn enthalpy(&self, temperature: f64) -> Result<f64, SublimationError> {
        if !temperature.is_finite() || !(MIN_T..=MAX_T).contains(&temperature) {
            return Err(SublimationError::TemperatureOutOfRange {
                temperature_k: temperature,
                minimum_k: MIN_T,
                maximum_k: MAX_T,
            });
        }
        let vapor_h = self
            .vapor
            .species_specific_enthalpy_j_kg(self.species, temperature)
            .map_err(|_| SublimationError::Invalid("original vapor caloric query failed"))?;
        Ok(vapor_h - self.law.properties(temperature)?.latent_heat_j_kg - self.datum)
    }

    fn heat_capacity(&self, temperature: f64) -> Result<f64, SublimationError> {
        let source = &self.vapor.thermodynamics().components()[0];
        let vapor_cp = source
            .evaluate(Temperature::new(temperature))
            .map_err(|_| SublimationError::Invalid("original vapor cp query failed"))?
            .properties()
            .cp()
            .value()
            / source.molar_mass().value();
        let cp = vapor_cp
            - self
                .law
                .properties(temperature)?
                .latent_heat_derivative_j_kg_k;
        if !cp.is_finite() || cp <= 0.0 {
            return Err(SublimationError::Invalid(
                "derived ice caloric branch is not stable",
            ));
        }
        Ok(cp)
    }

    fn temperature(&self, target: f64) -> Result<f64, SublimationError> {
        let (mut lower, mut upper) = (MIN_T, MAX_T);
        let low_h = self.enthalpy(lower)?;
        let high_h = self.enthalpy(upper)?;
        if !target.is_finite() || target < low_h || target > high_h {
            return Err(SublimationError::Invalid(
                "ice enthalpy leaves the 250--270 K chart",
            ));
        }
        let mut temperature = lower + (upper - lower) * ((target - low_h) / (high_h - low_h));
        for _ in 0..64 {
            let residual = self.enthalpy(temperature)? - target;
            if residual.abs() <= 1e-7 {
                return Ok(temperature);
            }
            if residual < 0.0 {
                lower = temperature;
            } else {
                upper = temperature;
            }
            let newton = temperature - residual / self.heat_capacity(temperature)?;
            temperature = if newton > lower && newton < upper {
                newton
            } else {
                lower + 0.5 * (upper - lower)
            };
        }
        Err(SublimationError::NonConverged)
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    let tolerance_scale = std::env::args()
        .nth(1)
        .map(|value| value.parse::<f64>())
        .transpose()?
        .unwrap_or(1.0);
    if !tolerance_scale.is_finite() || tolerance_scale <= 0.0 {
        return Err("local error-target scale must be finite and positive".into());
    }
    let duration_s = 0.002;
    let water = SpeciesId::new("H2O")?;
    let catalog = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![water_vapor(&water)?],
            Composition::new(CompositionBasis::MoleFraction, vec![1.0])?,
        )?,
        vec![transport(&water)?],
        Some(&water),
    )?;
    let law = WaterIceSublimation::try_new(917.0)?;
    let datum = catalog.species_specific_enthalpy_j_kg(&water, DATUM_T)?
        - law.properties(DATUM_T)?.latent_heat_j_kg;
    let caloric = IceCaloric {
        law: &law,
        vapor: &catalog,
        species: &water,
        datum,
    };
    // Exercise positivity and inversion over the full admitted branch, not just
    // the initial points. This is a numerical check, not an empirical cp claim.
    for i in 0..=20 {
        let temperature = MIN_T + f64::from(i);
        assert!(caloric.heat_capacity(temperature)? > 0.0);
        assert!((caloric.temperature(caloric.enthalpy(temperature)?)? - temperature).abs() < 1e-8);
    }
    let inverse = |h| caloric.temperature(h);
    let binding = GasSurfaceBinding {
        species: water.clone(),
        accommodation_coefficient: 1.0,
        solid_enthalpy_offset_j_kg: datum,
    };
    let area_m2 = 1e-3;
    let initial_thickness_m = 1e-5;
    let initial_mass = 917.0 * area_m2 * initial_thickness_m;
    let initial = [
        SolidSurfaceState {
            mass_kg: initial_mass,
            specific_enthalpy_j_kg: caloric.enthalpy(260.0)?,
            area_m2,
        },
        SolidSurfaceState {
            mass_kg: initial_mass,
            specific_enthalpy_j_kg: caloric.enthalpy(252.0)?,
            area_m2,
        },
    ];
    let surfaces = [
        GasChamberSurface {
            id: 10,
            law: &law,
            binding: &binding,
            solid: initial[0],
            temperature_from_enthalpy: &inverse,
        },
        GasChamberSurface {
            id: 20,
            law: &law,
            binding: &binding,
            solid: initial[1],
            temperature_from_enthalpy: &inverse,
        },
    ];
    let point = QueryPoint::new()
        .with_quantity(
            "temperature",
            typed(QuantityKind::AbsoluteTemperature),
            DATUM_T,
        )?
        .with_quantity("pressure", typed(QuantityKind::Pressure), 130.0)?;
    let initial_gas = ClosedGasParcel::new(catalog.clone(), point, 100e-6)?;
    let initial_gas_mass = initial_gas.mass_kg();
    let initial_gas_u = initial_gas.state().specific_internal_energy_j_kg();
    let controls = GasChamberAdvanceSolve {
        step: GasChamberSolve {
            surface: GasSurfaceSolve {
                gas: GasSpeciesBoundarySolve {
                    gas: GasEnergySolve {
                        minimum_temperature_k: MIN_T,
                        maximum_temperature_k: MAX_T,
                        energy_tolerance_j_kg: 1e-7,
                        maximum_iterations: 128,
                    },
                    mass_tolerance_kg: 1e-14,
                    energy_tolerance_j: 1e-8,
                },
                maximum_gross_inventory_fraction: 0.1,
                latent_compatibility_tolerance_j_kg: 1e-5,
            },
            maximum_surfaces: 2,
        },
        absolute_mass_tolerance_kg: 2e-11 * tolerance_scale,
        absolute_temperature_tolerance_k: 0.002 * tolerance_scale,
        relative_tolerance: 1e-6 * tolerance_scale,
        minimum_dt_s: 1e-10,
        // Deliberately start with the full duration. The original gross-stream
        // inventory bound and physical-state error estimate choose smaller steps.
        maximum_dt_s: duration_s,
        maximum_trials: 10000,
    };
    let proposal = advance_chamber_adaptive(
        &surfaces,
        &initial_gas,
        &catalog,
        duration_s,
        controls,
        || false,
    )?;
    let rows = proposal.endpoint().surfaces();
    assert_eq!(rows.len(), 2);
    assert_eq!((rows[0].id, rows[1].id), (10, 20));
    let solids = [rows[0].solid, rows[1].solid];
    let gas = proposal.gas();
    let audit = proposal.audit();

    let warm_loss = initial[0].mass_kg - solids[0].mass_kg;
    let cold_gain = solids[1].mass_kg - initial[1].mass_kg;
    let mass_residual = (gas.mass_kg() - initial_gas_mass)
        + (solids[0].mass_kg - initial[0].mass_kg)
        + (solids[1].mass_kg - initial[1].mass_kg);
    // Independent whole-history storage differences preserve the original NASA
    // formation-energy datum. The endpoint().step() report covers only the last
    // half-step; audit() contains the complete accepted molecular streams.
    let solid_energy_change: f64 = solids
        .iter()
        .zip(initial)
        .map(|(new, old)| {
            new.mass_kg * (new.specific_enthalpy_j_kg - old.specific_enthalpy_j_kg)
                + (new.mass_kg - old.mass_kg) * (old.specific_enthalpy_j_kg + datum)
        })
        .sum();
    let gas_energy_change = gas.mass_kg()
        * (gas.state().specific_internal_energy_j_kg() - initial_gas_u)
        + (gas.mass_kg() - initial_gas_mass) * initial_gas_u;
    let energy_residual = solid_energy_change + gas_energy_change;
    let final_temperatures = [
        caloric.temperature(solids[0].specific_enthalpy_j_kg)?,
        caloric.temperature(solids[1].specific_enthalpy_j_kg)?,
    ];
    assert!(warm_loss > 0.0 && cold_gain > 0.0);
    assert!(final_temperatures[0] < 260.0 && final_temperatures[1] > 252.0);
    assert!((gas.state().state().temperature - DATUM_T).abs() > 0.01);
    assert!(
        mass_residual.abs() < 1e-14,
        "whole-run mass residual: {mass_residual:e}"
    );
    assert!(
        energy_residual.abs() < 1e-8,
        "whole-run energy residual: {energy_residual:e}"
    );
    assert!(audit.mass_residual_kg.abs() < 1e-14 && audit.energy_residual_j.abs() < 1e-8);
    assert!((audit.surfaces[0].transferred_mass_kg - warm_loss).abs() < 1e-14);
    assert!((audit.surfaces[1].transferred_mass_kg + cold_gain).abs() < 1e-14);
    assert!(audit.maximum_accepted_error_ratio <= 1.0);
    assert_eq!(audit.elapsed_s, duration_s);
    assert!(warm_loss / initial_mass < 0.05 && cold_gain / initial_mass < 0.05);

    println!("Closed nonisothermal ice cold trap: {duration_s:.6} s");
    println!(
        "Warm ice: 260 K -> {:.6} K; loss {:.6} micrograms",
        final_temperatures[0],
        warm_loss * 1e9
    );
    println!(
        "Cold ice: 252 K -> {:.6} K; frost gain {:.6} micrograms",
        final_temperatures[1],
        cold_gain * 1e9
    );
    println!(
        "Vapor: 256 K -> {:.6} K; 130 Pa -> {:.6} Pa",
        gas.state().state().temperature,
        gas.species_partial_pressure_pa(&water)?
    );
    println!(
        "Adaptive solve: {} accepted intervals (two half-steps each), {} rejected trials, {} kernel calls",
        audit.accepted_intervals, audit.rejected_trials, audit.kernel_evaluations
    );
    println!(
        "Largest accepted local error ratio: {:.6}",
        audit.maximum_accepted_error_ratio
    );
    for surface in &audit.surfaces {
        println!(
            "Surface {} molecular streams: outgoing {:.6}, incoming {:.6} micrograms",
            surface.id,
            surface.outgoing_mass_kg * 1e9,
            surface.incoming_mass_kg * 1e9
        );
    }
    println!("Whole-run mass residual: {mass_residual:.3e} kg");
    println!("Whole-run energy residual: {energy_residual:.3e} J");
    println!(
        "Derived cp_ice at 250/270 K: {:.3}/{:.3} J/(kg K)",
        caloric.heat_capacity(MIN_T)?,
        caloric.heat_capacity(MAX_T)?
    );
    Ok(())
}
