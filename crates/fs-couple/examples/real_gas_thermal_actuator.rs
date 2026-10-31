//! Two finite real gases, one source-bound moving solid, and reciprocal heat.
//!
//! Run: cargo run -p fs-couple --example real_gas_thermal_actuator -- 8000 325 200 0.000002
//! Arguments: solid density [kg/m3], initial solid temperature [K], steps, dt [s].
//! Compare densities 1000 and 8000 with the SAME remaining inputs. The original
//! 0.1 x 0.1 x 0.01 m body supplies BOTH mechanical inertia and thermal mass.
//! CSV includes each gas inventory/density, body temperature, wall motion,
//! actual convection groups, and an independent original-law energy account.
//! There is no external heat, force or exterior-pressure reservoir. Signed
//! body/gas heat and both gases' pressure work are internal to the assembly.
//!
//! The two synthetic PR mixtures retain variable-Cp calorics and source transport.
//! Forced convection uses the ORIGINAL flat-plate correlation with properties at
//! each actual bulk PR energy midpoint: this is a deliberate bulk-property
//! approximation for modest temperature contrasts, not a film-state closure.
//! Prescribed 0.5 m/s relative flow is parallel to each 0.1 m-long face; it is not
//! piston velocity. This thermal correlation models no fan power, flow momentum,
//! boundary-layer storage, or fluid shear heating. The remaining solid area is
//! insulated. Every midpoint checks the source conductivity's lumped Biot limit.
//!
//! These are homogeneous, frozen-composition, supported single-phase gases and
//! an all-solid, rigid, constant-density body; there is no phase flash, thermal
//! expansion, melting or acoustic propagation. Halve dt and double steps for a
//! fixed-duration refinement comparison. Original clock/source/correlation/work
//! admissions remain active; synthetic data do not establish measured accuracy.

use fs_convection::{CorrelationId, CorrelationInputs, NusseltEvaluation, ThermalConductivity};
use fs_couple::cavity_phs::piston::PistonSolve;
use fs_couple::cavity_phs::piston::double_acting::{DoublePistonInput, DoublePistonSolve};
use fs_couple::cavity_phs::piston::material::MaterialDoubleActingGasPiston;
use fs_couple::cavity_phs::piston::thermal::double_acting::DoubleHeatExchangeSolve;
use fs_couple::cavity_phs::piston::thermal::{ConductanceSample, HeatExchangeSolve};
use fs_evidence::ValidityDomain;
use fs_ga::Vec3;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::fluid_transport::FluidTransport;
use fs_material::phase::{
    EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase, UniformEnthalpyStepInput,
};
use fs_material::state_point::{
    MaterialPropertySelection, ResolvedMaterialStatePoint, SPECIFIC_HEAT_CAPACITY_DIMS,
    ScalarAdmissibility, ScalarPropertyRequirement, THERMAL_CONDUCTIVITY_DIMS,
    resolve_material_state_point,
};
use fs_material::thermochemical_gas::energy::GasEnergySolve;
use fs_material::thermochemical_gas::nonideal::{
    PengRobinsonSpeciesBinding, ResolvedThermochemicalRealGasState,
    ThermochemicalPengRobinsonModel, energy::parcel::IsochoricRealGasParcel,
};
use fs_mbd_material::{RegionBinding, VolumeMoments, resolve_body};
use fs_qty::{Density, Length, Pressure, Temperature};
use fs_thermochem::{SpeciesId, eos::peng_robinson::PengRobinsonPhaseV1};
use std::io::Write;

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
type InterfaceSample = ConductanceSample<(NusseltEvaluation, ResolvedMaterialStatePoint)>;

#[path = "../../fs-material/examples/material_gas_spring/source.rs"]
mod source;

const AREAS_M2: [f64; 2] = [0.01, 0.008];
const LENGTH_M: f64 = 0.1;
const SPEED_M_S: f64 = 0.5;
const SOLID_TEMPERATURE_RANGE_K: [f64; 2] = [200.0, 500.0];
const ENERGY_TOLERANCE_J: f64 = 2e-7;

fn gas(x_a: f64, temperature_k: f64, pressure_pa: f64) -> Result<IsochoricRealGasParcel> {
    let data = format!(
        "A {x_a} 0.028 3.5 0.001 2e-5 0.03\nB {} 0.004 2.5 0 2e-5 0.03\n",
        1.0 - x_a,
    );
    let model = ThermochemicalPengRobinsonModel::new(
        source::load(&data, "real_gas_thermal_actuator manufactured binary")?,
        vec![
            PengRobinsonSpeciesBinding {
                species: SpeciesId::new("A")?,
                critical_temperature: Temperature::new(190.564),
                critical_pressure: Pressure::new(4.5992e6),
                acentric_factor: 0.01142,
            },
            PengRobinsonSpeciesBinding {
                species: SpeciesId::new("B")?,
                critical_temperature: Temperature::new(150.687),
                critical_pressure: Pressure::new(4.863e6),
                acentric_factor: -0.00219,
            },
        ],
        vec![0.0, 0.03, 0.03, 0.0],
    )?;
    Ok(IsochoricRealGasParcel::new(
        model,
        source::point(temperature_k, pressure_pa)?,
        0.001,
        PengRobinsonPhaseV1::Vapor,
    )?)
}

fn solid_card(density: f64) -> Result<MaterialCard> {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        ("density", Density::DIMS, density),
        ("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS, 500.0),
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, 200.0),
    ] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: SOLID_TEMPERATURE_RANGE_K.map(|t| (t, value)).to_vec(),
                dims,
            },
            validity: ValidityDomain::unconstrained().with("T", 200.0, 500.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "SYNTHETIC actuator solid: constant density, Cp and conductivity".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })?;
    }
    Ok(MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic actuator solid".into(),
            phase: "solid".into(),
            process: "authored constant properties".into(),
            revision: 0,
        },
        claims,
        vec![],
    )?)
}

fn property(state: &ResolvedMaterialStatePoint, name: &str) -> Result<f64> {
    Ok(state
        .property(name)
        .ok_or_else(|| format!("missing resolved solid property {name}"))?
        .value_si())
}

fn convection(
    gas: &ResolvedThermochemicalRealGasState,
    solid: &ResolvedMaterialStatePoint,
    area_m2: f64,
    solid_characteristic_length_m: f64,
) -> Result<InterfaceSample> {
    // One actual PR snapshot supplies rho, Cp, viscosity and conductivity.
    let transport = gas.transport_state()?;
    let re = transport.reynolds(SPEED_M_S, LENGTH_M)?;
    let nu = fs_convection::evaluate(
        CorrelationId::FlatPlateLaminarAverage,
        CorrelationInputs::forced(re, transport.prandtl()),
    )?;
    let coefficient = nu.heat_transfer_coefficient(
        ThermalConductivity::new(transport.thermal_conductivity_w_m_k()),
        Length::new(LENGTH_M),
    )?;
    let h = coefficient.value.value();
    let biot = h * solid_characteristic_length_m / property(solid, "thermal_conductivity")?;
    if !biot.is_finite() || biot > 0.1 {
        return Err(format!(
            "uniform solid refuses midpoint Biot={biot}; spatial heat storage is required"
        )
        .into());
    }
    Ok(ConductanceSample {
        conductance_w_k: area_m2 * h,
        source: (nu, solid.clone()),
    })
}

fn internal_energy_change(
    initial: &IsochoricRealGasParcel,
    gas: &IsochoricRealGasParcel,
) -> Result<f64> {
    Ok(initial.mass_kg()
        * initial
            .model()
            .specific_internal_energy_change_between_densities(
                initial.density_kg_m3(),
                initial.state().temperature_k(),
                gas.density_kg_m3(),
                gas.state().temperature_k(),
            )?)
}

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() > 4 {
        return Err("usage: real_gas_thermal_actuator [solid_density_kg_m3] [body_temperature_k] [steps] [dt_s]".into());
    }
    let density: f64 = args.first().map_or(Ok(8000.0), |v| v.parse())?;
    let temperature: f64 = args.get(1).map_or(Ok(325.0), |v| v.parse())?;
    let steps: usize = args.get(2).map_or(Ok(200), |v| v.parse())?;
    let dt: f64 = args.get(3).map_or(Ok(2e-6), |v| v.parse())?;
    if !density.is_finite()
        || density <= 0.0
        || !temperature.is_finite()
        || !(200.0..=500.0).contains(&temperature)
        || !(1..=100_000).contains(&steps)
        || !dt.is_finite()
        || dt <= 0.0
        || !(dt * steps as f64).is_finite()
    {
        return Err("positive finite density/dt, solid temperature 200..=500 K, and 1..=100000 steps required".into());
    }

    let card = solid_card(density)?;
    let point = QueryPoint::new().with("T", temperature)?;
    let geometry = VolumeMoments::cuboid(Vec3::new(0.1, 0.1, 0.01))?;
    let body = resolve_body(
        &[RegionBinding {
            key: "moving-solid",
            material: &card,
            point: &point,
            geometry,
            selection: MaterialPropertySelection::SingleClaimOnly,
        }],
        1,
        || false,
    )?;
    let requirements = [
        ("density", Density::DIMS),
        ("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
    ]
    .map(|(name, dims)| {
        ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
    })
    .into_iter()
    .collect::<std::result::Result<Vec<_>, _>>()?;
    let solid = resolve_material_state_point(
        &card,
        &point,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )?;
    let cp = property(&solid, "specific_heat_capacity")?;
    let solid_density = property(&solid, "density")?;
    // Integrate this source's declared constant Cp on its all-solid interval.
    // The same card/density supplies the canonical inertia and the enthalpy chart.
    let curve = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        card.content_hash(),
        SolidLiquidPhase::Solid,
        SOLID_TEMPERATURE_RANGE_K
            .map(|t| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: cp * (t - SOLID_TEMPERATURE_RANGE_K[0]),
                temperature_k: t,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: solid_density,
            })
            .to_vec(),
    )?;
    let solid_mass = body.inertia().mass();
    let solid_volume = geometry.volume_m3();
    let exposed_area = AREAS_M2.iter().sum::<f64>();
    let solid_length = solid_volume / exposed_area;
    let initial_enthalpy = cp * (temperature - SOLID_TEMPERATURE_RANGE_K[0]);
    let mut thermal = UniformEnthalpyStepInput {
        curve: &curve,
        initial: curve.state_at_specific_enthalpy(initial_enthalpy)?,
        mass_kg: solid_mass,
        volume_m3: solid_volume,
        surface_area_m2: exposed_area,
        internal_heat_j: 0.0,
        duration_s: dt,
    };
    let mut actuator = MaterialDoubleActingGasPiston::new(
        gas(0.75, 300.0, 5e6)?,
        gas(0.35, 350.0, 3.8e6)?,
        body,
        AREAS_M2,
        0.0,
    )?;
    let initial = actuator.clone();
    let gas_controls = GasEnergySolve {
        minimum_temperature_k: 250.0,
        maximum_temperature_k: 450.0,
        energy_tolerance_j_kg: 1e-7,
        maximum_iterations: 128,
    };
    let mut cumulative_heat = [0.0; 2];
    let mut re = [0.0; 2];
    let mut h = [0.0; 2];
    let mut step_residual = 0.0;
    let stdout = std::io::stdout();
    let mut output = std::io::BufWriter::new(stdout.lock());
    writeln!(
        output,
        "step,time_s,solid_density_kg_m3,solid_mass_kg,solid_temperature_k,wall_velocity_m_s,primary_mass_kg,primary_density_kg_m3,primary_volume_m3,primary_temperature_k,primary_pressure_pa,opposing_mass_kg,opposing_density_kg_m3,opposing_volume_m3,opposing_temperature_k,opposing_pressure_pa,heat_to_primary_j,heat_to_opposing_j,primary_midpoint_re,opposing_midpoint_re,primary_midpoint_h_w_m2_k,opposing_midpoint_h_w_m2_k,primary_energy_change_j,opposing_energy_change_j,solid_enthalpy_change_j,kinetic_energy_change_j,total_energy_residual_j,step_energy_residual_j"
    )?;
    for step in 0..=steps {
        let current_solid = solid.at_query_point(
            &card,
            &QueryPoint::new().with("T", thermal.initial.temperature_k())?,
        )?;
        let current_gases = actuator.gases();
        let samples = [
            convection(
                current_gases[0].state(),
                &current_solid,
                AREAS_M2[0],
                solid_length,
            )?,
            convection(
                current_gases[1].state(),
                &current_solid,
                AREAS_M2[1],
                solid_length,
            )?,
        ];
        if step == 0 {
            for i in 0..2 {
                re[i] = samples[i].source.0.groups()["Re"];
                h[i] = samples[i].conductance_w_k / AREAS_M2[i];
            }
        } else {
            // Bound only a small supported displacement of either chamber.
            let dv = 0.02
                * current_gases[0]
                    .volume_m3()
                    .min(current_gases[1].volume_m3() * AREAS_M2[0] / AREAS_M2[1]);
            let volume = current_gases[0].volume_m3();
            let faces = std::array::from_fn(|i| {
                let contrast = (thermal.initial.temperature_k()
                    - current_gases[i].state().temperature_k())
                .abs();
                // A 1K margin admits the motion's small midpoint temperature
                // change near equilibrium; out-of-bracket physics still refuses.
                let bound = 4.0 * dt * samples[i].conductance_w_k * (contrast + 1.0);
                HeatExchangeSolve {
                    heat_bracket_j: [-bound, bound],
                    source_tolerance_j: 1e-10,
                    energy_tolerance_j: ENERGY_TOLERANCE_J,
                    maximum_iterations: 128,
                }
            });
            let candidate = actuator
                .thermal_candidate(
                    DoublePistonInput {
                        duration_s: dt,
                        heat_j: [0.0; 2],
                        external_force_n: 0.0,
                    },
                    thermal,
                    DoublePistonSolve {
                        piston: PistonSolve {
                            volume_bracket_m3: [volume - dv, volume + dv],
                            gas: gas_controls,
                            momentum_tolerance_n_s: 1e-9,
                            energy_tolerance_j: ENERGY_TOLERANCE_J,
                            maximum_iterations: 128,
                        },
                        opposing_gas: gas_controls,
                        volume_tolerance_m3: 1e-18,
                    },
                    DoubleHeatExchangeSolve {
                        faces,
                        maximum_joint_trials: 256,
                    },
                    |gases, midpoint| -> Result<[InterfaceSample; 2]> {
                        let resolved = solid.at_query_point(
                            &card,
                            &QueryPoint::new().with("T", midpoint.temperature_k())?,
                        )?;
                        Ok([
                            convection(&gases[0], &resolved, AREAS_M2[0], solid_length)?,
                            convection(&gases[1], &resolved, AREAS_M2[1], solid_length)?,
                        ])
                    },
                    || false,
                )
                .map_err(|error| error.to_string())?;
            let (next, accepted_thermal, report) = candidate.into_parts();
            for i in 0..2 {
                cumulative_heat[i] += report.exchanged_heat_j[i];
                re[i] = accepted_thermal.report[i].0.groups()["Re"];
                h[i] = report.conductance_w_k[i] / AREAS_M2[i];
            }
            step_residual = report.total_energy_residual_j;
            thermal.initial = accepted_thermal.state;
            actuator = next;
        }
        let gases = actuator.gases();
        let gas_energy = [
            internal_energy_change(initial.gases()[0], gases[0])?,
            internal_energy_change(initial.gases()[1], gases[1])?,
        ];
        let body_energy =
            solid_mass * (thermal.initial.specific_enthalpy_j_kg() - initial_enthalpy);
        let velocity = actuator.piston().primary().velocity_m_s();
        let kinetic = 0.5 * actuator.body().inertia().mass() * velocity * velocity;
        // Independent endpoint audit: pressure work and convection cancel
        // internally; there is no exterior-pressure supply or second wall mass.
        let total = gas_energy[0] + gas_energy[1] + body_energy + kinetic;
        if !total.is_finite() || total.abs() > (step.max(1) as f64) * ENERGY_TOLERANCE_J {
            return Err(format!("step {step}: original-law assembly energy residual {total} J exceeds its accumulated budget").into());
        }
        write!(output, "{step}")?;
        for value in [
            step as f64 * dt,
            solid_density,
            solid_mass,
            thermal.initial.temperature_k(),
            velocity,
            gases[0].mass_kg(),
            gases[0].density_kg_m3(),
            gases[0].volume_m3(),
            gases[0].state().temperature_k(),
            gases[0].state().pressure_pa(),
            gases[1].mass_kg(),
            gases[1].density_kg_m3(),
            gases[1].volume_m3(),
            gases[1].state().temperature_k(),
            gases[1].state().pressure_pa(),
            cumulative_heat[0],
            cumulative_heat[1],
            re[0],
            re[1],
            h[0],
            h[1],
            gas_energy[0],
            gas_energy[1],
            body_energy,
            kinetic,
            total,
            step_residual,
        ] {
            write!(output, ",{value:.12e}")?;
        }
        writeln!(output)?;
    }
    output.flush()?;
    Ok(())
}
