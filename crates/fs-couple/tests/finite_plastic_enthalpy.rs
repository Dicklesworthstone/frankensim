//! Native finite plasticity, reciprocal heat, retained inertia and source identity.
//!
//! All material values are synthetic. Endpoint loads and boundary heat are
//! manufactured from the independently tested native point law; the coupled
//! solve must recover an actual free displacement and temperature together.

#[path = "finite_plastic_enthalpy/combined.rs"]
mod combined;
#[path = "finite_plastic_enthalpy/conduction.rs"]
mod conduction;
#[path = "finite_plastic_enthalpy/contact.rs"]
#[cfg(feature = "experimental-solid-contact-heat")]
mod contact;
#[path = "finite_plastic_enthalpy/creep.rs"]
mod creep;
#[path = "finite_plastic_enthalpy/ductile.rs"]
mod ductile;
#[path = "finite_plastic_enthalpy/ductile_sheet.rs"]
mod ductile_sheet;
#[path = "finite_plastic_enthalpy/rate.rs"]
mod rate;
#[path = "finite_plastic_enthalpy/sheet.rs"]
mod sheet;
#[path = "finite_plastic_enthalpy/tabulated.rs"]
mod tabulated;

use std::collections::BTreeMap;

use fs_couple::iqn_ils::{IqnIlsConfig, IqnIlsError, driver::CouplingMethod};
use fs_couple::solid_enthalpy::finite_plastic::{
    CaloricReference, FiniteThermoplasticCell, FiniteThermoplasticConfig, FiniteThermoplasticError,
    FiniteThermoplasticMaterial, FiniteThermoplasticQuery, FiniteThermoplasticSolid,
};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy, LawId, LawParameter,
    MaterialCard, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::plastic::finite::thermal::FiniteJ2ThermalIncrement;
use fs_material::plastic::finite::thermal::source::{
    ResolvedAffineThermalFiniteJ2, THERMAL_FINITE_J2_LAW_ID,
};
use fs_material::state_point::DENSITY_PROPERTY;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, QuantitySpec};
use fs_solid::hyper3::HyperTetSettings;
use fs_solid::hyper3::plastic::FinitePlasticStepControls;
use fs_solid::hyper3::plastic::thermal::FiniteThermoplasticControls;
use fs_solid::linear3::TetAssemblyBudget;
use fs_thermal::enthalpy::*;

const L: f64 = 0.01;
const V: f64 = L * L * L / 6.0;
const RHO: f64 = 1000.0;
const CP: f64 = 500.0;
const CELL: CellId = CellId(100);
const NODES: [[f64; 3]; 4] = [[0.0; 3], [L, 0.0, 0.0], [0.0, L, 0.0], [0.0, 0.0, L]];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const IDENTITY: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
const PA: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const PA_K: Dims = Dims([-1, 1, -2, -1, 0, 0]);
const K: Dims = Dims([0, 0, 0, 1, 0, 0]);
const CP_DIMS: Dims = Dims([2, 0, -2, -1, 0, 0]);
const CP_NAME: &str = "specific_heat_capacity";

fn temperature_quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

fn source_context(temperature_k: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("temperature", temperature_quantity(), temperature_k)
        .unwrap()
}

fn provenance() -> Provenance {
    Provenance {
        source: "synthetic finite thermoplastic coupled regression; not measured calibration"
            .into(),
        license: "test fixture".into(),
        artifact: None,
    }
}

fn model() -> ConstitutiveModelCard {
    ConstitutiveModelCard {
        law: LawId(THERMAL_FINITE_J2_LAW_ID.into()),
        law_version: 1,
        parameters: [
            ("reference_temperature_k", 300.0, K),
            ("minimum_temperature_k", 280.0, K),
            ("maximum_temperature_k", 420.0, K),
            ("shear_modulus_pa", 80.0e6, PA),
            ("bulk_modulus_pa", 160.0e6, PA),
            ("yield_stress_pa", 2.0e6, PA),
            ("hardening_modulus_pa", 8.0e6, PA),
            ("shear_modulus_slope_pa_k", -40.0e3, PA_K),
            ("bulk_modulus_slope_pa_k", -60.0e3, PA_K),
            ("yield_stress_slope_pa_k", -2.0e3, PA_K),
            ("hardening_modulus_slope_pa_k", -5.0e3, PA_K),
            (
                "logarithmic_expansion_per_k",
                1.0e-4,
                Dims([0, 0, 0, -1, 0, 0]),
            ),
            (
                "logarithmic_expansion_slope_per_k2",
                1.0e-7,
                Dims([0, 0, 0, -2, 0, 0]),
            ),
            ("minimum_mechanical_jacobian", 0.25, Dims::NONE),
            ("maximum_elastic_log_strain", 0.5, Dims::NONE),
            ("maximum_plastic_log_strain", 2.0, Dims::NONE),
            ("maximum_accumulated_plastic_strain", 3.0, Dims::NONE),
        ]
        .into_iter()
        .map(|(name, value, dims)| (name.into(), LawParameter { value, dims }))
        .collect(),
        state_schema_version: 1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with_quantity(
            "temperature",
            temperature_quantity(),
            280.0,
            420.0,
        ),
        sources: vec![fs_blake3::hash_bytes(
            b"synthetic finite thermoplastic material model fixture",
        )],
        provenance: provenance(),
    }
}

fn scalar(name: &str, dims: Dims, value: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained().with_quantity(
            "temperature",
            temperature_quantity(),
            280.0,
            420.0,
        ),
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: provenance(),
    }
}

struct Fixture {
    card: MaterialCard,
    caloric: ResolvedHeatCapacityCurve,
}
impl Fixture {
    fn new(cp: f64, identity: u32) -> Self {
        Self::with_density(cp, identity, scalar(DENSITY_PROPERTY, Density::DIMS, RHO))
    }
    fn with_density(cp: f64, identity: u32, density: PropertyClaim) -> Self {
        Self::with_model(cp, identity, density, model())
    }
    fn with_model(
        cp: f64,
        identity: u32,
        density: PropertyClaim,
        constitutive: ConstitutiveModelCard,
    ) -> Self {
        let mut claims = ClaimSet::new();
        claims.insert_claim(scalar(CP_NAME, CP_DIMS, cp)).unwrap();
        claims.insert_claim(density).unwrap();
        let card = MaterialCard::assemble(
            MaterialStateId {
                chemistry: format!("synthetic original-volume finite thermal J2 {identity}"),
                phase: "solid".into(),
                process: "numerical regression".into(),
                revision: 0,
            },
            claims,
            vec![constitutive],
        )
        .unwrap();
        let caloric = ResolvedHeatCapacityCurve::try_new(
            &card,
            HeatCapacitySourceQuery {
                heat_capacity: &PropertyKey::new(CP_NAME, CP_DIMS),
                density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
                temperature_axis: "temperature",
                lower: &source_context(280.0),
                upper: &source_context(420.0),
                selection: SelectionPolicy::SingleClaimOnly,
            },
            SolidLiquidPhase::Solid,
            -20.0 * cp,
            1.0e-6,
        )
        .unwrap();
        Self { card, caloric }
    }
    fn query(&self) -> FiniteThermoplasticQuery<'_> {
        FiniteThermoplasticQuery {
            card: &self.card,
            model_index: 0,
            caloric: &self.caloric,
            context: source_context(300.0),
            temperature_axis: "temperature".into(),
            caloric_reference: CaloricReference::StressFreeZeroPressure,
        }
    }
}

fn context<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 3,
                kernel_id: 0x464a3254,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn settings() -> HyperTetSettings {
    HyperTetSettings {
        load_steps: 1,
        relative_tolerance: 1.0e-11,
        ..HyperTetSettings::default()
    }
}
fn controls() -> FiniteThermoplasticControls {
    FiniteThermoplasticControls {
        material: FinitePlasticStepControls {
            equilibrium_work_tolerance_j: 1.0e-8,
            maximum_endpoint_work_defect_j: 5.0,
            yield_tolerance_pa: 1.0e-4,
        },
        maximum_numerical_kinetic_loss_j: 0.01,
    }
}
fn thermal_config() -> StepConfig {
    StepConfig {
        max_dt_s: 0.25,
        energy_tolerance_j: 1.0e-11,
    }
}
fn config() -> FiniteThermoplasticConfig {
    FiniteThermoplasticConfig {
        maximum_iterations: 100,
        temperature_tolerance_k: 1.0e-8,
        energy_tolerance_j: 2.0e-9,
        enthalpy_relaxation: 0.7,
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.16e}, expected={expected:.16e}, tolerance={tolerance:e}"
    );
}

fn initial_velocity(inertia: bool) -> [[f64; 3]; 4] {
    let mut velocity = [[0.0; 3]; 4];
    if inertia {
        // A resolvable numerical kinetic loss must remain outside material heat.
        velocity[1][0] = 1.0;
    }
    velocity
}

fn mass_action(value: &[[f64; 3]; 4]) -> [[f64; 3]; 4] {
    core::array::from_fn(|a| {
        core::array::from_fn(|axis| {
            (0..4)
                .map(|b| RHO * V / 20.0 * if a == b { 2.0 } else { 1.0 } * value[b][axis])
                .sum()
        })
    })
}
fn dot(a: &[[f64; 3]; 4], b: &[[f64; 3]; 4]) -> f64 {
    a.iter()
        .flatten()
        .zip(b.iter().flatten())
        .map(|(x, y)| x * y)
        .sum()
}

struct Endpoint {
    time_s: f64,
    temperature_k: f64,
    old_temperature_k: f64,
    deformation: [f64; 9],
    displacement: [[f64; 3]; 4],
    old_displacement: [[f64; 3]; 4],
    velocity: [[f64; 3]; 4],
    old_velocity: [[f64; 3]; 4],
    inertial_force: [[f64; 3]; 4],
    free_force_n: f64,
    boundary_heat_j: f64,
    increment: FiniteJ2ThermalIncrement,
}
impl Endpoint {
    fn prescribed(&self) -> Vec<(usize, f64)> {
        self.displacement
            .iter()
            .flatten()
            .enumerate()
            .filter(|(dof, _)| *dof != 3)
            .map(|(dof, value)| (dof, *value))
            .collect()
    }
}

fn path(fixture: &Fixture, inertia: bool) -> Vec<Endpoint> {
    let source = ResolvedAffineThermalFiniteJ2::try_from_model_card_at(
        &fixture.card.models()[0],
        &source_context(300.0),
        "temperature",
    )
    .unwrap();
    let law = source.law();
    let mut old_state = law.initial_state();
    let mut old_f = IDENTITY;
    let mut old_temperature = 300.0;
    let mut old_u = [[0.0; 3]; 4];
    let mut old_v = initial_velocity(inertia);
    let mut path = Vec::new();
    for (index, (temperature, axial_log)) in [(310.0, 0.16_f64), (330.0, 0.16), (315.0, -0.10)]
        .into_iter()
        .enumerate()
    {
        // Unit TOTAL volume, independent of the thermal stretch. The second
        // endpoint holds exactly the same total F while the yield radius falls.
        let f = [
            axial_log.exp(),
            0.0,
            0.0,
            0.0,
            (-0.5 * axial_log).exp(),
            0.0,
            0.0,
            0.0,
            (-0.5 * axial_log).exp(),
        ];
        let u = core::array::from_fn(|node| {
            core::array::from_fn(|axis| (f[3 * axis + axis] - 1.0) * NODES[node][axis])
        });
        let velocity = core::array::from_fn(|node| {
            core::array::from_fn(|axis| u[node][axis] - old_u[node][axis])
        });
        let acceleration = core::array::from_fn(|node| {
            core::array::from_fn(|axis| velocity[node][axis] - old_v[node][axis])
        });
        let inertial_force = if inertia {
            mass_action(&acceleration)
        } else {
            [[0.0; 3]; 4]
        };
        let increment = law
            .checked_increment(&old_f, &old_state, old_temperature, &f, temperature, 1.0e-4)
            .unwrap();
        let free_force_n = V / L * increment.mechanical.first_piola_pa[0] + inertial_force[1][0];
        let boundary_heat_j =
            RHO * V * CP * (temperature - old_temperature) - V * increment.caloric_heat_j_m3();
        old_state = increment.mechanical.state.clone();
        path.push(Endpoint {
            time_s: (index + 1) as f64,
            temperature_k: temperature,
            old_temperature_k: old_temperature,
            deformation: f,
            displacement: u,
            old_displacement: old_u,
            velocity,
            old_velocity: old_v,
            inertial_force,
            free_force_n,
            boundary_heat_j,
            increment,
        });
        old_f = f;
        old_temperature = temperature;
        old_u = u;
        old_v = velocity;
    }
    path
}

fn boundary(power_w: f64) -> Boundary {
    Boundary {
        ambient_k: 300.0,
        htc_w_m2_k: 0.0,
        reservoir: None,
        source_w: power_w,
    }
}

fn initial<'a>(fixture: &'a Fixture, endpoints: &[Endpoint]) -> DepositionSimulation<'a> {
    let mut events = vec![TimedEvent {
        time_s: 0.0,
        event: Event::Birth(DepositedCell {
            id: CELL,
            material_id: MaterialId(0),
            segment_id: 1,
            body_id: 7,
            mass_kg: RHO * V,
            specific_h_j_kg: 0.0,
            geometry: CellGeometry {
                start_m: [0.0; 3],
                end_m: [L; 3],
                volume_m3: V,
                surface_area_m2: 3.0 * L * L,
            },
            boundary: boundary(endpoints[0].boundary_heat_j),
        }),
    }];
    for (i, endpoint) in endpoints.iter().enumerate().skip(1) {
        events.push(TimedEvent {
            time_s: i as f64,
            event: Event::Boundary(CELL, boundary(endpoint.boundary_heat_j)),
        });
    }
    DepositionSimulation::new(
        events,
        &[(MaterialId(0), fixture.caloric.curve() as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 8,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(0.0, thermal_config(), || Ok(()))
    .unwrap()
    .0
}

fn bind<'a>(
    fixture: &'a Fixture,
    thermal: &DepositionSimulation<'_>,
    inertia: bool,
    cx: &Cx<'_>,
) -> FiniteThermoplasticSolid<'a> {
    let solid = FiniteThermoplasticSolid::try_new(
        thermal,
        &NODES,
        &TETS,
        vec![FiniteThermoplasticCell {
            cell: CELL,
            query: fixture.query(),
        }],
        TetAssemblyBudget::standard(),
        cx,
    )
    .unwrap();
    if inertia {
        solid.with_inertia(&initial_velocity(true), cx).unwrap()
    } else {
        solid
    }
}

#[test]
fn source_bound_free_dof_recovers_heating_reyield_reversal_and_separates_numerical_loss() {
    context(|cx| {
        let fixture = Fixture::new(CP, 0);
        for inertia in [false, true] {
            let feedback = FiniteThermoplasticConfig {
                enthalpy_relaxation: if inertia { 0.7 } else { 1.0 },
                ..config()
            };
            let endpoints = path(&fixture, inertia);
            let mut thermal = initial(&fixture, &endpoints);
            let mut solid = bind(&fixture, &thermal, inertia, cx);
            let mut previous_alpha = 0.0;
            for (index, expected) in endpoints.iter().enumerate() {
                let prescribed = expected.prescribed();
                let candidate = solid
                    .advance(
                        &thermal,
                        &(),
                        expected.time_s,
                        &prescribed,
                        &[(3, expected.free_force_n)],
                        [0.0; 3],
                        settings(),
                        controls(),
                        thermal_config(),
                        feedback,
                        cx,
                    )
                    .unwrap();
                let report = candidate.report();
                let actual = candidate.solid().mechanics();
                close(
                    candidate.thermal().temperature_k(CELL).unwrap(),
                    expected.temperature_k,
                    2.0e-7,
                );
                close(actual.temperatures_k()[0], expected.temperature_k, 2.0e-7);
                for (u, want) in actual
                    .displacement_m()
                    .iter()
                    .flatten()
                    .zip(expected.displacement.iter().flatten())
                {
                    close(*u, *want, 2.0e-11);
                }
                for (f, want) in actual.element_deformation()[0]
                    .iter()
                    .zip(expected.deformation)
                {
                    close(*f, want, 2.0e-9);
                }
                let history = &actual.element_states()[0];
                close(
                    history.alpha,
                    expected.increment.mechanical.state.alpha,
                    3.0e-9,
                );
                for (metric, want) in history.inverse_plastic_cauchy_green.iter().zip(
                    expected
                        .increment
                        .mechanical
                        .state
                        .inverse_plastic_cauchy_green,
                ) {
                    close(*metric, want, 3.0e-9);
                }
                assert!(history.alpha > previous_alpha + 1.0e-6);
                previous_alpha = history.alpha;
                if index == 0 {
                    assert!(history.alpha > 0.14);
                }
                if index == 1 {
                    close(report.material_energy.internal_endpoint_work_j, 0.0, 2.0e-8);
                    assert!(report.material_energy.plastic_dissipation_j > 1.0e-6);
                }
                if index == 2 {
                    assert!(history.inverse_plastic_cauchy_green[0] > 1.0);
                    assert!(history.alpha > 0.3);
                }
                let material = &actual.families()[0];
                assert_eq!(
                    material.query().card.content_hash(),
                    fixture.card.content_hash()
                );
                assert_eq!(material.family().model_card(), &fixture.card.models()[0]);
                assert!(std::ptr::eq(material.caloric_source(), &fixture.caloric));
                close(material.reference_density_kg_m3(), RHO, 0.0);
                close(candidate.thermal().cells()[&CELL].mass_kg, RHO * V, 0.0);
                close(
                    report.material_energy.plastic_dissipation_j,
                    V * expected.increment.mechanical.plastic_dissipation_j_m3,
                    2.0e-9,
                );
                close(
                    report.material_energy.reversible_heat_j,
                    V * expected.increment.reversible_heat_j_m3,
                    2.0e-9,
                );
                close(
                    report.material_energy.caloric_heat_j,
                    V * expected.increment.caloric_heat_j_m3(),
                    2.0e-9,
                );
                close(
                    report.material_energy.internal_energy_change_j,
                    V * expected.increment.internal_energy_change_j_m3(),
                    2.0e-9,
                );
                close(
                    report.material_energy.endpoint_work_defect_j,
                    V * expected.increment.endpoint_work_defect_j_m3,
                    2.0e-9,
                );
                close(report.external_heat_j, expected.boundary_heat_j, 2.0e-9);
                close(
                    report.thermal_audit.stored_energy_change_j,
                    RHO * V * CP * (expected.temperature_k - expected.old_temperature_k),
                    2.0e-8,
                );
                close(
                    report.thermal_audit.stored_energy_change_j - report.external_heat_j,
                    report.material_energy.plastic_dissipation_j
                        + report.material_energy.reversible_heat_j,
                    2.0e-9,
                );
                close(report.thermal_residual_j, 0.0, config().energy_tolerance_j);
                close(report.combined_residual_j, 0.0, config().energy_tolerance_j);
                assert!(report.iterations > 1);
                let mut total_work = expected.free_force_n
                    * (expected.displacement[1][0] - expected.old_displacement[1][0]);
                for &(dof, reaction) in &report.equilibrium.reactions_n {
                    let node = dof / 3;
                    let axis = dof % 3;
                    let material_force = if node == 0 {
                        -(0..3)
                            .map(|column| {
                                expected.increment.mechanical.first_piola_pa[3 * axis + column]
                            })
                            .sum::<f64>()
                            * V
                            / L
                    } else {
                        expected.increment.mechanical.first_piola_pa[3 * axis + node - 1] * V / L
                    };
                    close(
                        reaction,
                        material_force + expected.inertial_force[node][axis],
                        2.0e-5,
                    );
                    total_work += reaction
                        * (expected.displacement[node][axis]
                            - expected.old_displacement[node][axis]);
                }
                let inertial_account = if inertia {
                    let account = report.inertia.as_ref().unwrap();
                    let delta_v = core::array::from_fn(|node| {
                        core::array::from_fn(|axis| {
                            expected.velocity[node][axis] - expected.old_velocity[node][axis]
                        })
                    });
                    close(account.reference_mass_kg, RHO * V, 1.0e-18);
                    close(
                        account.kinetic_before_j,
                        0.5 * dot(&expected.old_velocity, &mass_action(&expected.old_velocity)),
                        2.0e-12,
                    );
                    close(
                        account.kinetic_after_j,
                        0.5 * dot(&expected.velocity, &mass_action(&expected.velocity)),
                        2.0e-12,
                    );
                    close(
                        account.numerical_loss_j,
                        0.5 * dot(&delta_v, &mass_action(&delta_v)),
                        2.0e-12,
                    );
                    if index == 0 {
                        assert!(account.numerical_loss_j > 1.0e-6);
                    }
                    for (v, want) in actual
                        .velocity_mps()
                        .unwrap()
                        .iter()
                        .flatten()
                        .zip(expected.velocity.iter().flatten())
                    {
                        close(*v, *want, 2.0e-11);
                    }
                    account.kinetic_change_j + account.numerical_loss_j
                } else {
                    assert!(report.inertia.is_none());
                    assert!(actual.velocity_mps().is_none());
                    0.0
                };
                close(
                    total_work + report.external_heat_j,
                    report.material_energy.internal_energy_change_j
                        + report.thermal_audit.stored_energy_change_j
                        + report.material_energy.endpoint_work_defect_j
                        + inertial_account,
                    2.0e-8,
                );
                // A successful proposal also leaves its borrowed predecessor intact.
                close(solid.time_s(), index as f64, 0.0);
                close(thermal.time_s(), index as f64, 0.0);
                let (next_thermal, next_solid, (), _) = candidate.into_parts();
                thermal = next_thermal;
                solid = next_solid;
            }
        }
    });
}


#[test]
fn g1_iqn_and_picard_recover_native_thermoplastic_heat_motion_and_history() {
    context(|cx| {
        let fixture = Fixture::new(CP, 0);
        for inertia in [false, true] {
            let endpoints = path(&fixture, inertia);
            let expected = &endpoints[0];
            let thermal = initial(&fixture, &endpoints);
            let solid = bind(&fixture, &thermal, inertia, cx);
            let original_cells = thermal.cells().clone();
            let original_history = solid.mechanics().element_states().to_vec();
            let run = |method| {
                solid
                    .advance_with_coupling(
                        &thermal,
                        &(),
                        expected.time_s,
                        &expected.prescribed(),
                        &[(3, expected.free_force_n)],
                        [0.0; 3],
                        settings(),
                        controls(),
                        thermal_config(),
                        config(),
                        method,
                        cx,
                    )
                    .unwrap()
            };
            let picard = run(CouplingMethod::RelaxedPicard);
            let accelerated = run(CouplingMethod::IqnIls(IqnIlsConfig::default()));
            for candidate in [&picard, &accelerated] {
                let actual = candidate.solid().mechanics();
                let report = candidate.report();
                // The point-law oracle fixes load and outside heat independently
                // of the coupled solver, including the free displacement.
                close(candidate.thermal().temperature_k(CELL).unwrap(),
                    expected.temperature_k, 2e-7);
                close(actual.temperatures_k()[0], expected.temperature_k, 2e-7);
                for (u, expected_u) in actual.displacement_m().iter().flatten()
                    .zip(expected.displacement.iter().flatten())
                {
                    close(*u, *expected_u, 2e-11);
                }
                close(actual.element_states()[0].alpha,
                    expected.increment.mechanical.state.alpha, 3e-9);
                close(report.material_energy.caloric_heat_j,
                    V * expected.increment.caloric_heat_j_m3(), 2e-9);
                close(report.external_heat_j, expected.boundary_heat_j, 2e-9);
                assert!(report.temperature_defect_k <= config().temperature_tolerance_k);
                assert!(report.heat_defect_j <= config().energy_tolerance_j);
                close(report.thermal_residual_j, 0.0, config().energy_tolerance_j);
                close(report.combined_residual_j, 0.0, config().energy_tolerance_j);
                assert!(std::ptr::eq(actual.families()[0].caloric_source(), &fixture.caloric));
                assert_eq!(actual.families()[0].query().card.content_hash(),
                    fixture.card.content_hash());
                assert_eq!(candidate.thermal().cells()[&CELL].mass_kg, RHO * V);
                if inertia {
                    let account = report.inertia.as_ref().unwrap();
                    close(account.reference_mass_kg, RHO * V, 1e-18);
                    assert!(account.numerical_loss_j > 1e-6);
                    for (v, expected_v) in actual.velocity_mps().unwrap().iter().flatten()
                        .zip(expected.velocity.iter().flatten())
                    {
                        close(*v, *expected_v, 2e-11);
                    }
                } else {
                    assert!(report.inertia.is_none());
                }
                assert!(report.iterations > 1);
            }
            assert_eq!(thermal.cells(), &original_cells);
            assert_eq!(solid.mechanics().element_states(), original_history);
            assert_eq!(thermal.time_s(), 0.0);
            assert_eq!(solid.time_s(), 0.0);
        }
    });
}

#[test]
fn iqn_refusal_budget_and_retry_preserve_nonvirgin_thermoplastic_owners() {
    context(|cx| {
        let fixture = Fixture::new(CP, 0);
        let endpoints = path(&fixture, true);
        let thermal = initial(&fixture, &endpoints);
        let solid = bind(&fixture, &thermal, true, cx);
        let first = &endpoints[0];
        let accepted = solid.advance(
            &thermal, &(), first.time_s, &first.prescribed(), &[(3, first.free_force_n)],
            [0.0; 3], settings(), controls(), thermal_config(), config(), cx,
        ).unwrap();
        let (thermal, solid, (), _) = accepted.into_parts();
        let old_cells = thermal.cells().clone();
        let old_history = solid.mechanics().element_states().to_vec();
        let old_u = solid.mechanics().displacement_m().to_vec();
        let old_v = solid.mechanics().velocity_mps().unwrap().to_vec();
        let expected = &endpoints[1];
        let method = CouplingMethod::IqnIls(IqnIlsConfig::default());
        let propose = |feedback, method| {
            solid.advance_with_coupling(
                &thermal, &(), expected.time_s, &expected.prescribed(),
                &[(3, expected.free_force_n)], [0.0; 3], settings(), controls(),
                thermal_config(), feedback, method, cx,
            )
        };
        assert!(matches!(
            propose(config(), CouplingMethod::IqnIls(IqnIlsConfig {
                max_history: 0,
                ..IqnIlsConfig::default()
            })),
            Err(FiniteThermoplasticError::Accelerator(IqnIlsError::InvalidHistoryLimit(0)))
        ));
        let observer = LateRefusal { completed_substeps: 0 };
        let late = solid.advance_with_coupling(
            &thermal, &observer, expected.time_s, &expected.prescribed(),
            &[(3, expected.free_force_n)], [0.0; 3], settings(), controls(),
            thermal_config(), config(), method, cx,
        );
        assert!(matches!(
            late,
            Err(FiniteThermoplasticError::Observer(DepositionError::Cancelled))
        ));
        assert_eq!(observer.completed_substeps, 0);
        assert!(matches!(
            propose(FiniteThermoplasticConfig { maximum_iterations: 1, ..config() }, method),
            Err(FiniteThermoplasticError::NotConverged { iterations: 1, .. })
        ));
        assert_eq!(thermal.cells(), &old_cells);
        assert_eq!(solid.mechanics().element_states(), old_history);
        assert_eq!(solid.mechanics().displacement_m(), old_u);
        assert_eq!(solid.mechanics().velocity_mps().unwrap(), old_v);
        assert_eq!(thermal.time_s(), 1.0);
        assert_eq!(solid.time_s(), 1.0);

        let retried = propose(config(), method).unwrap();
        close(retried.thermal().temperature_k(CELL).unwrap(),
            expected.temperature_k, 2e-7);
        close(retried.solid().mechanics().element_states()[0].alpha,
            expected.increment.mechanical.state.alpha, 3e-9);
        assert!(retried.solid().mechanics().element_states()[0].alpha > old_history[0].alpha);
        close(retried.report().combined_residual_j, 0.0, config().energy_tolerance_j);
        // Even an accepted proposal has not adopted either predecessor's state.
        assert_eq!(thermal.cells(), &old_cells);
        assert_eq!(solid.mechanics().element_states(), old_history);
    });
}

#[derive(Clone, Default)]
struct HeatPulse {
    watts: f64,
    delivered_j: f64,
    intervals: Vec<(f64, f64)>,
}
impl DepositionObserver for HeatPulse {
    type Error = DepositionError;
    fn next_transition_time_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Option<f64>, Self::Error> {
        checkpoint()?;
        Ok((thermal.time_s() < 0.125).then_some(0.125))
    }
    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        Ok([(
            CELL,
            if thermal.time_s() < 0.125 {
                self.watts
            } else {
                0.0
            },
        )]
        .into())
    }
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        Ok(())
    }
    fn step(
        &mut self,
        start_s: f64,
        duration_s: f64,
        _: &BTreeMap<CellId, f64>,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.intervals.push((start_s, duration_s));
        if start_s < 0.125 {
            self.delivered_j += self.watts * duration_s;
        }
        Ok(())
    }
}

#[test]
fn source_transition_survives_finite_plastic_feedback_without_smearing_heat() {
    context(|cx| {
        let fixture = Fixture::new(CP, 0);
        let endpoints = path(&fixture, false);
        let expected = &endpoints[0];
        let original = initial(&fixture, &endpoints);
        let mut cell = original.cells()[&CELL].clone();
        cell.boundary.source_w = 0.0;
        let thermal = DepositionSimulation::new(
            vec![TimedEvent {
                time_s: 0.0,
                event: Event::Birth(cell),
            }],
            &[(MaterialId(0), fixture.caloric.curve() as &dyn ThermalLaw)],
            DepositionBudget {
                max_cells: 1,
                max_links: 0,
                max_events: 1,
                max_substeps: 16,
            },
        )
        .unwrap()
        .advance_to(0.0, thermal_config(), || Ok(()))
        .unwrap()
        .0;
        let solid = bind(&fixture, &thermal, false, cx);
        // Deliver the same independently manufactured endpoint heat as the
        // native plasticity oracle, but finish the pulse inside a normal step.
        let observer = HeatPulse {
            watts: expected.boundary_heat_j / 0.125,
            ..HeatPulse::default()
        };
        let next = solid
            .advance(
                &thermal,
                &observer,
                expected.time_s,
                &expected.prescribed(),
                &[(3, expected.free_force_n)],
                [0.0; 3],
                settings(),
                controls(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        assert_eq!(
            next.observer().intervals,
            [
                (0.0, 0.125),
                (0.125, 0.25),
                (0.375, 0.25),
                (0.625, 0.25),
                (0.875, 0.125)
            ]
        );
        close(next.observer().delivered_j, expected.boundary_heat_j, 1e-12);
        close(
            next.report().external_heat_j,
            expected.boundary_heat_j,
            1e-9,
        );
        close(
            next.thermal().temperature_k(CELL).unwrap(),
            expected.temperature_k,
            2e-7,
        );
        close(
            next.solid().mechanics().element_states()[0].alpha,
            expected.increment.mechanical.state.alpha,
            3e-9,
        );
        close(
            next.report().combined_residual_j,
            0.0,
            config().energy_tolerance_j,
        );
        assert_eq!(observer.delivered_j, 0.0);
        assert!(observer.intervals.is_empty());
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(solid.time_s(), 0.0);
    });
}

#[derive(Clone)]
struct LateRefusal {
    completed_substeps: usize,
}
impl DepositionObserver for LateRefusal {
    type Error = DepositionError;
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        Ok(())
    }
    fn step(
        &mut self,
        _: f64,
        _: f64,
        _: &BTreeMap<CellId, f64>,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.completed_substeps += 1;
        if self.completed_substeps == 2 {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    }
}

#[test]
fn late_heat_refusal_and_unconverged_feedback_preserve_nonvirgin_motion_and_all_owners() {
    context(|cx| {
        let fixture = Fixture::new(CP, 0);
        let endpoints = path(&fixture, true);
        let thermal = initial(&fixture, &endpoints);
        let solid = bind(&fixture, &thermal, true, cx);
        let first = &endpoints[0];
        let accepted = solid
            .advance(
                &thermal,
                &(),
                first.time_s,
                &first.prescribed(),
                &[(3, first.free_force_n)],
                [0.0; 3],
                settings(),
                controls(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        let (thermal, solid, (), _) = accepted.into_parts();
        let old_cells = thermal.cells().clone();
        let old_u = solid.mechanics().displacement_m().to_vec();
        let old_v = solid.mechanics().velocity_mps().unwrap().to_vec();
        let old_f = solid.mechanics().element_deformation().to_vec();
        let old_history = solid.mechanics().element_states().to_vec();
        let old_stress = solid.mechanics().element_first_piola_pa().to_vec();
        let old_t = solid.mechanics().temperatures_k().to_vec();
        let old_free = solid.mechanics().free_energy_j();
        let old_internal = solid.mechanics().internal_energy_j();
        let next = &endpoints[1];
        let observer = LateRefusal {
            completed_substeps: 0,
        };
        let late = solid.advance(
            &thermal,
            &observer,
            next.time_s,
            &next.prescribed(),
            &[(3, next.free_force_n)],
            [0.0; 3],
            settings(),
            controls(),
            thermal_config(),
            config(),
            cx,
        );
        assert!(matches!(
            late,
            Err(FiniteThermoplasticError::Observer(
                DepositionError::Cancelled
            ))
        ));
        assert_eq!(observer.completed_substeps, 0);
        let short = FiniteThermoplasticConfig {
            maximum_iterations: 1,
            ..config()
        };
        let unconverged = solid.advance(
            &thermal,
            &(),
            next.time_s,
            &next.prescribed(),
            &[(3, next.free_force_n)],
            [0.0; 3],
            settings(),
            controls(),
            thermal_config(),
            short,
            cx,
        );
        assert!(matches!(
            unconverged,
            Err(FiniteThermoplasticError::NotConverged { iterations: 1, .. })
        ));
        assert_eq!(thermal.cells(), &old_cells);
        assert_eq!(
            solid.cells(),
            &old_cells.values().cloned().collect::<Vec<_>>()
        );
        assert_eq!(solid.mechanics().displacement_m(), old_u);
        assert_eq!(solid.mechanics().velocity_mps().unwrap(), old_v);
        assert_eq!(solid.mechanics().element_deformation(), old_f);
        assert_eq!(solid.mechanics().element_states(), old_history);
        assert_eq!(solid.mechanics().element_first_piola_pa(), old_stress);
        assert_eq!(solid.mechanics().temperatures_k(), old_t);
        assert_eq!(solid.mechanics().free_energy_j(), old_free);
        assert_eq!(solid.mechanics().internal_energy_j(), old_internal);
        assert_eq!(thermal.time_s(), 1.0);
        assert_eq!(solid.time_s(), 1.0);
        assert_eq!(solid.mechanics().reference_mass_kg(), Some(RHO * V));
        assert!(old_history[0].alpha > 0.14);
        // The same accepted pair remains usable after both discarded proposals.
        let retried = solid
            .advance(
                &thermal,
                &(),
                next.time_s,
                &next.prescribed(),
                &[(3, next.free_force_n)],
                [0.0; 3],
                settings(),
                controls(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        close(
            retried.thermal().temperature_k(CELL).unwrap(),
            330.0,
            2.0e-7,
        );
        assert!(retried.solid().mechanics().element_states()[0].alpha > old_history[0].alpha);
    });
}

#[test]
fn foreign_caloric_identity_replaced_objects_and_nonpositive_total_capacity_are_refused() {
    context(|cx| {
        let fixture = Fixture::new(CP, 0);
        let admitted = FiniteThermoplasticMaterial::try_new(fixture.query(), 0.0).unwrap();
        assert!(admitted.minimum_total_heat_capacity_j_m3_k() > 0.0);
        let mut pressured = fixture.query();
        pressured.context = pressured
            .context
            .with_quantity(
                "environment_pressure",
                QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::Pressure,
                    ValueForm::Static,
                )),
                1.0e5,
            )
            .unwrap();
        let refused = FiniteThermoplasticMaterial::try_new(pressured, 0.0).unwrap_err();
        assert!(
            refused
                .to_string()
                .contains("zero-pressure caloric reference")
        );
        let foreign = Fixture::new(CP, 1);
        let mut foreign_query = fixture.query();
        foreign_query.caloric = &foreign.caloric;
        assert!(matches!(
            FiniteThermoplasticMaterial::try_new(foreign_query, 0.0),
            Err(FiniteThermoplasticError::Invalid(_))
        ));
        let inadequate = Fixture::new(1.0, 0);
        let refused = FiniteThermoplasticMaterial::try_new(inadequate.query(), 0.0).unwrap_err();
        assert!(refused.to_string().contains("heat capacity"));
        let endpoints = path(&fixture, false);
        let thermal = initial(&fixture, &endpoints);
        let duplicate_curve = fixture.caloric.clone();
        let mut replaced = fixture.query();
        replaced.caloric = &duplicate_curve;
        // Same source values and hash do not replace the exact caloric object
        // already carried by the original deposition simulation.
        assert!(FiniteThermoplasticMaterial::try_new(replaced.clone(), 0.0).is_ok());
        let refused = FiniteThermoplasticSolid::try_new(
            &thermal,
            &NODES,
            &TETS,
            vec![FiniteThermoplasticCell {
                cell: CELL,
                query: replaced,
            }],
            TetAssemblyBudget::standard(),
            cx,
        );
        assert!(matches!(refused, Err(FiniteThermoplasticError::Invalid(_))));
        assert_eq!(thermal.time_s(), 0.0);
        close(thermal.temperature_k(CELL).unwrap(), 300.0, 0.0);

        let mut density = scalar(DENSITY_PROPERTY, Density::DIMS, RHO);
        density.value = PropertyValue::Curve {
            abscissa: "temperature".into(),
            abscissa_dims: K,
            knots: vec![(280.0, 1010.0), (420.0, 940.0)],
            dims: Density::DIMS,
        };
        density.interpolation = InterpolationPolicy::LinearInside;
        let variable = Fixture::with_density(CP, 0, density);
        let material = FiniteThermoplasticMaterial::try_new(variable.query(), 0.0).unwrap();
        close(material.reference_density_kg_m3(), RHO, 1.0e-9);
        // A source-derived sensible caloric segment keeps the claim's affine
        // rho(T) law (e838056254), so the chart's reference density is the
        // claim's value at the reference temperature, as the material's is.
        close(
            variable
                .caloric
                .curve()
                .state_at_specific_enthalpy(0.0)
                .unwrap()
                .bulk_density_kg_m3(),
            RHO,
            1.0e-9,
        );
        let endpoints = path(&variable, true);
        let thermal = initial(&variable, &endpoints);
        let solid = bind(&variable, &thermal, true, cx);
        let next = &endpoints[0];
        let candidate = solid
            .advance(
                &thermal,
                &(),
                next.time_s,
                &next.prescribed(),
                &[(3, next.free_force_n)],
                [0.0; 3],
                settings(),
                controls(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        close(
            candidate.thermal().temperature_k(CELL).unwrap(),
            310.0,
            2.0e-7,
        );
        close(
            candidate.solid().mechanics().reference_mass_kg().unwrap(),
            RHO * V,
            1.0e-16,
        );
        close(candidate.thermal().cells()[&CELL].mass_kg, RHO * V, 0.0);
        let warm_h = candidate.thermal().cells()[&CELL].specific_h_j_kg;
        assert!(
            (variable
                .caloric
                .curve()
                .state_at_specific_enthalpy(warm_h)
                .unwrap()
                .bulk_density_kg_m3()
                - RHO)
                .abs()
                > 1.0
        );
    });
}

#[test]
fn heterogeneous_cells_conduct_heat_with_opposing_reversible_sources_and_local_histories() {
    context(|cx| {
        let first = Fixture::new(CP, 0);
        let mut other_model = model();
        other_model
            .parameters
            .get_mut("yield_stress_pa")
            .unwrap()
            .value = 4.0e6;
        other_model
            .parameters
            .get_mut("hardening_modulus_pa")
            .unwrap()
            .value = 12.0e6;
        let second = Fixture::with_model(
            250.0,
            1,
            scalar(DENSITY_PROPERTY, Density::DIMS, RHO),
            other_model,
        );
        let fixtures = [&first, &second];
        let temperatures = [310.0, 320.0];
        let capacities = [RHO * V * CP, RHO * V * 250.0];
        let deformation: [[f64; 9]; 2] = core::array::from_fn(|i| {
            let (v, a): (f64, f64) = [(0.04, 0.14), (-0.04, 0.20)][i];
            [
                (v + a).exp(),
                0.0,
                0.0,
                0.0,
                (v - 0.5 * a).exp(),
                0.0,
                0.0,
                0.0,
                (v - 0.5 * a).exp(),
            ]
        });
        let increments: Vec<_> = fixtures
            .iter()
            .enumerate()
            .map(|(i, fixture)| {
                let source = ResolvedAffineThermalFiniteJ2::try_from_model_card_at(
                    &fixture.card.models()[0],
                    &source_context(300.0),
                    "temperature",
                )
                .unwrap();
                source
                    .law()
                    .checked_increment(
                        &IDENTITY,
                        &source.law().initial_state(),
                        300.0,
                        &deformation[i],
                        temperatures[i],
                        1.0e-4,
                    )
                    .unwrap()
            })
            .collect();
        let point_heat = increments
            .iter()
            .map(|e| V * e.caloric_heat_j_m3())
            .collect::<Vec<_>>();
        assert!(point_heat[0] < 0.0 && point_heat[1] > 0.0);

        // Independent two-capacity explicit conduction recurrence over exactly
        // four 0.25 s substeps. Total heat and temperature difference determine
        // the constant net powers needed for the two distinct target endpoints.
        let conductance = 0.01;
        let decay = conductance * (capacities[0].recip() + capacities[1].recip());
        let factor = 1.0 - decay * thermal_config().max_dt_s;
        let power_sum =
            capacities[0] * (temperatures[0] - 300.0) + capacities[1] * (temperatures[1] - 300.0);
        let difference_rate = decay * (temperatures[0] - temperatures[1]) / (1.0 - factor.powi(4));
        let power_first = (difference_rate + power_sum / capacities[1])
            / (capacities[0].recip() + capacities[1].recip());
        let net_power = [power_first, power_sum - power_first];
        let boundary_power = [net_power[0] - point_heat[0], net_power[1] - point_heat[1]];
        let nodes: [[f64; 3]; 8] = core::array::from_fn(|i| {
            let mut node = NODES[i % 4];
            if i >= 4 {
                node[0] += 2.0 * L;
            }
            node
        });
        let tetrahedra = [[0, 1, 2, 3], [4, 5, 6, 7]];
        let prescribed: Vec<_> = (0..24)
            .map(|dof| {
                let node = dof / 3;
                let axis = dof % 3;
                (
                    dof,
                    (deformation[node / 4][3 * axis + axis] - 1.0) * NODES[node % 4][axis],
                )
            })
            .collect();
        let mut events = Vec::new();
        for i in 0..2 {
            events.push(TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id: CellId(100 + i as u64),
                    material_id: MaterialId(i as u64),
                    segment_id: i as u64,
                    body_id: 7 + i as u64,
                    mass_kg: RHO * V,
                    specific_h_j_kg: 0.0,
                    geometry: CellGeometry {
                        start_m: [2.0 * L * i as f64, 0.0, 0.0],
                        end_m: [2.0 * L * i as f64 + L, L, L],
                        volume_m3: V,
                        surface_area_m2: 3.0 * L * L,
                    },
                    boundary: boundary(boundary_power[i]),
                }),
            });
        }
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 7,
                a: CellId(100),
                b: CellId(101),
                area_m2: L * L,
                conductance_w_k: conductance,
            }),
        });
        let laws: Vec<_> = fixtures
            .iter()
            .enumerate()
            .map(|(i, fixture)| {
                (
                    MaterialId(i as u64),
                    fixture.caloric.curve() as &dyn ThermalLaw,
                )
            })
            .collect();
        let thermal = DepositionSimulation::new(
            events,
            &laws,
            DepositionBudget {
                max_cells: 2,
                max_links: 1,
                max_events: 4,
                max_substeps: 1000,
            },
        )
        .unwrap()
        .advance_to(0.0, thermal_config(), || Ok(()))
        .unwrap()
        .0;
        let solid = FiniteThermoplasticSolid::try_new(
            &thermal,
            &nodes,
            &tetrahedra,
            fixtures
                .iter()
                .enumerate()
                .map(|(i, fixture)| FiniteThermoplasticCell {
                    cell: CellId(100 + i as u64),
                    query: fixture.query(),
                })
                .collect(),
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let candidate = solid
            .advance(
                &thermal,
                &(),
                1.0,
                &prescribed,
                &[],
                [0.0; 3],
                settings(),
                controls(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        let mut local_link_heat = [0.0; 2];
        for i in 0..2 {
            let id = CellId(100 + i as u64);
            close(
                candidate.thermal().temperature_k(id).unwrap(),
                temperatures[i],
                2.0e-7,
            );
            close(
                candidate.solid().mechanics().temperatures_k()[i],
                temperatures[i],
                2.0e-7,
            );
            let actual = &candidate.solid().mechanics().element_states()[i];
            close(actual.alpha, increments[i].mechanical.state.alpha, 3.0e-9);
            for (metric, want) in actual
                .inverse_plastic_cauchy_green
                .iter()
                .zip(increments[i].mechanical.state.inverse_plastic_cauchy_green)
            {
                close(*metric, want, 3.0e-9);
            }
            assert_eq!(
                candidate.solid().mechanics().families()[i]
                    .family()
                    .model_card(),
                &fixtures[i].card.models()[0]
            );
            let stored = candidate.thermal().cells()[&id].specific_h_j_kg * RHO * V;
            local_link_heat[i] = stored - boundary_power[i] - point_heat[i];
            close(
                local_link_heat[i],
                capacities[i] * (temperatures[i] - 300.0) - net_power[i],
                2.0e-8,
            );
        }
        assert!(local_link_heat[0] > 0.01 && local_link_heat[1] < -0.01);
        close(local_link_heat.iter().sum(), 0.0, 2.0e-9);
        let report = candidate.report();
        close(
            report.material_energy.caloric_heat_j,
            point_heat.iter().sum(),
            2.0e-9,
        );
        close(report.external_heat_j, boundary_power.iter().sum(), 2.0e-9);
        close(
            report.thermal_audit.stored_energy_change_j,
            power_sum,
            2.0e-8,
        );
        assert!(report.heat_defect_j <= config().energy_tolerance_j);
        close(report.thermal_residual_j, 0.0, config().energy_tolerance_j);
        close(report.combined_residual_j, 0.0, config().energy_tolerance_j);
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(solid.time_s(), 0.0);
        assert!(
            solid
                .mechanics()
                .element_states()
                .iter()
                .all(|state| state.alpha == 0.0)
        );
    });
}

#[test]
fn hot_source_bound_impact_retains_caloric_timeline_and_continues_with_new_momentum() {
    use fs_contact::frictional_multicontact::elliptic::EllipticFriction;
    use fs_contact::frictional_multicontact::{FrictionalContactConfig, NormalVelocityLaw};
    use fs_contact::interface_binding::thermal_friction::elliptic::SourceFrictionConditions;
    use fs_couple::finite_solid::contact::FiniteSolidRigidContactError;
    use fs_couple::inelastic_impact::multiple::{
        FiniteSolidContactCounterpart, InelasticMixedContact, InelasticMixedContactControls,
        RigidImpactPoint,
    };
    use fs_couple::solid_enthalpy::finite_plastic::FiniteThermoplasticImpactError;
    use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
    use fs_solid::hyper3::plastic::thermal::{SolidImpulseControls, TetImpulsePoint};
    use fs_solver::CheckedCgConfig;

    context(|cx| {
        let fixture = Fixture::new(CP, 0);
        let endpoints = path(&fixture, true);
        let original_thermal = initial(&fixture, &endpoints);
        let original_solid = bind(&fixture, &original_thermal, true, cx);
        let first = &endpoints[0];
        let accepted = original_solid
            .advance(
                &original_thermal,
                &(),
                first.time_s,
                &first.prescribed(),
                &[(3, first.free_force_n)],
                [0.0; 3],
                settings(),
                controls(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        let (thermal, solid, (), _) = accepted.into_parts();
        let old = solid.mechanics();
        let old_velocity = old.velocity_mps().unwrap().to_vec();
        assert!(old.element_states()[0].alpha > 0.14);
        close(old.temperatures_k()[0], 310.0, 2e-7);
        let point = TetImpulsePoint {
            tetrahedron: 0,
            weights: [0.25; 4],
        };
        let trace = old.point_kinematics(point, cx).unwrap();
        let mass = RHO * V;
        let properties = [MassProperties::new(
            mass,
            Vec3::ZERO,
            Vec3::new(mass * L * L, mass * L * L, mass * L * L),
        )
        .unwrap()];
        let bodies = [RigidBodyState::new(
            Pose::new(
                Vec3::new(
                    trace.position_m[0],
                    trace.position_m[1],
                    trace.position_m[2],
                ),
                UnitQuaternion::IDENTITY,
            )
            .unwrap(),
            Vec3::new(
                mass * (trace.velocity_mps[0] - 1.0),
                mass * trace.velocity_mps[1],
                mass * trace.velocity_mps[2],
            ),
            Vec3::ZERO,
        )
        .unwrap()];
        let contacts = [InelasticMixedContact {
            id: 17,
            point: RigidImpactPoint {
                body: 0,
                arm_body: Vec3::ZERO,
            },
            counterpart: FiniteSolidContactCounterpart::Solid(point),
            normal_world: Vec3::new(1.0, 0.0, 0.0),
            tangent_world: Vec3::new(0.0, 1.0, 0.0),
            activation_tolerance_m: 1e-12,
            friction_coefficient: 0.0,
            normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.5 },
        }];
        let impact_controls = InelasticMixedContactControls {
            maximum_rigid_bodies: 1,
            geometry_tolerance_m: 1e-12,
            velocity_tolerance_m_s: 1e-10,
            linear_tolerance_kg_m_s: 1e-12,
            angular_tolerance_kg_m2_s: 1e-12,
            energy_tolerance_j: 1e-12,
            solid_impulse: SolidImpulseControls {
                solve: CheckedCgConfig {
                    tolerance: 1e-13,
                    max_iterations: 64,
                    max_corrections: 2,
                },
                maximum_workspace_components: 4096,
                linear_tolerance_kg_m_s: 1e-12,
                angular_tolerance_kg_m2_s: 1e-12,
                work_tolerance_j: 1e-12,
            },
            contact: FrictionalContactConfig {
                maximum_contacts: 1,
                velocity_tolerance: 1e-10,
                impulse_tolerance: 1e-12,
                complementarity_tolerance: 1e-12,
                energy_tolerance: 1e-12,
                ..Default::default()
            },
        };
        let old_cells = thermal.cells().clone();
        let event = solid
            .mixed_contact_candidate(
                &thermal,
                &[],
                &bodies,
                &properties,
                &contacts,
                [0.0; 3],
                impact_controls,
                cx,
            )
            .unwrap();
        let elliptic = solid
            .mixed_elliptic_contact_candidate(
                &thermal,
                &[],
                &bodies,
                &properties,
                &contacts.map(|contact| {
                    contact.with_friction(EllipticFriction {
                        coefficients_mu: [0.0; 2],
                        principal_axis: [1.0, 0.0],
                    })
                }),
                [0.0; 3],
                impact_controls,
                cx,
            )
            .unwrap();
        assert_eq!(
            elliptic.solid.mechanics().velocity_mps(),
            event.solid.mechanics().velocity_mps()
        );
        assert_eq!(elliptic.solid.cells(), event.solid.cells());
        assert_eq!(elliptic.solid.time_s(), event.solid.time_s());
        assert_eq!(
            elliptic.rigid_bodies[0].linear_momentum_world(),
            event.rigid_bodies[0].linear_momentum_world()
        );
        close(
            elliptic.report.contact.dissipated_energy,
            event.report.contact.dissipated_energy,
            1e-15,
        );
        // At the centroid, the full consistent mass gives mobility 1/mass and
        // a uniform jump. Equal rigid mass, approach 1 and e=1/2 give J=3m/4.
        let impulse = 0.75 * mass;
        close(event.report.contact_mobility_per_kg[0], 2.0 / mass, 1e-8);
        close(event.report.contact.dissipated_energy, 0.1875 * mass, 1e-12);
        close(event.report.kinetic_change_j, -0.1875 * mass, 1e-12);
        close(event.report.support_work_j, 0.0, 1e-15);
        close(
            event.rigid_bodies[0].linear_momentum_world().x,
            bodies[0].linear_momentum_world().x + impulse,
            1e-12,
        );
        let after = event.solid.mechanics();
        for node in 0..4 {
            for axis in 0..3 {
                close(
                    after.velocity_mps().unwrap()[node][axis],
                    old.velocity_mps().unwrap()[node][axis] - if axis == 0 { 0.75 } else { 0.0 },
                    1e-10,
                );
            }
        }
        assert_eq!(event.solid.time_s(), solid.time_s());
        assert_eq!(event.solid.cells(), solid.cells());
        assert_eq!(thermal.cells(), &old_cells);
        assert_eq!(after.element_states(), old.element_states());
        assert_eq!(after.element_first_piola_pa(), old.element_first_piola_pa());
        assert_eq!(after.displacement_m(), old.displacement_m());
        assert_eq!(after.temperatures_k(), old.temperatures_k());
        assert_eq!(after.free_energy_j(), old.free_energy_j());
        assert_eq!(after.internal_energy_j(), old.internal_energy_j());
        assert_eq!(after.reference_mass_kg(), Some(mass));
        assert!(std::ptr::eq(
            after.families()[0].caloric_source(),
            &fixture.caloric
        ));
        assert_eq!(
            after.families()[0].query().card.content_hash(),
            fixture.card.content_hash()
        );

        // Neither a stale clock nor heat advanced without the coupled plastic
        // source can be substituted for this accepted caloric snapshot.
        let uncoupled = original_thermal
            .advance_to(first.time_s, thermal_config(), || Ok(()))
            .unwrap()
            .0;
        assert_eq!(uncoupled.time_s(), thermal.time_s());
        assert_ne!(
            uncoupled.cells()[&CELL].specific_h_j_kg.to_bits(),
            old_cells[&CELL].specific_h_j_kg.to_bits()
        );
        for wrong in [&original_thermal, &uncoupled] {
            assert!(matches!(
                solid.mixed_contact_candidate(
                    wrong,
                    &[],
                    &bodies,
                    &properties,
                    &contacts,
                    [0.0; 3],
                    impact_controls,
                    cx,
                ),
                Err(FiniteThermoplasticImpactError::Source(
                    FiniteThermoplasticError::Invalid(_)
                ))
            ));
            assert!(matches!(
                solid.mixed_source_elliptic_contact_candidate(
                    wrong,
                    &[],
                    &bodies,
                    &properties,
                    &contacts.map(|contact| contact.with_friction(())),
                    &[],
                    &[],
                    [0.0; 3],
                    impact_controls,
                    |_, _| -> Result<SourceFrictionConditions, ()> {
                        panic!("stale caloric state must refuse before interface sampling")
                    },
                    cx,
                ),
                Err(FiniteThermoplasticImpactError::Source(
                    FiniteThermoplasticError::Invalid(_)
                ))
            ));
        }
        assert!(matches!(
            solid.mixed_source_elliptic_contact_candidate(
                &thermal,
                &[],
                &bodies,
                &properties,
                &contacts.map(|contact| contact.with_friction(())),
                &[],
                &[],
                [0.0; 3],
                impact_controls,
                |_, _| -> Result<SourceFrictionConditions, ()> {
                    panic!("missing original interface histories must refuse before sampling")
                },
                cx,
            ),
            Err(FiniteThermoplasticImpactError::Contact(
                FiniteSolidRigidContactError::Invalid(_)
            ))
        ));
        assert_eq!(thermal.cells(), &old_cells);
        assert_eq!(old.velocity_mps().unwrap(), old_velocity);

        // Retain the next manufactured deformation/temperature, but correct
        // its free force by the ORIGINAL consistent M times the impact jump.
        let jump: [[f64; 3]; 4] = core::array::from_fn(|node| {
            core::array::from_fn(|axis| {
                after.velocity_mps().unwrap()[node][axis] - old.velocity_mps().unwrap()[node][axis]
            })
        });
        let next = &endpoints[1];
        let duration = next.time_s - event.solid.time_s();
        let continued = event
            .solid
            .advance(
                &thermal,
                &(),
                next.time_s,
                &next.prescribed(),
                &[(3, next.free_force_n - mass_action(&jump)[1][0] / duration)],
                [0.0; 3],
                settings(),
                controls(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        close(continued.solid().time_s(), next.time_s, 0.0);
        close(
            continued.thermal().temperature_k(CELL).unwrap(),
            next.temperature_k,
            2e-7,
        );
        close(
            continued
                .report()
                .inertia
                .as_ref()
                .unwrap()
                .kinetic_before_j,
            event.report.solid_impulse.kinetic_after_j,
            2e-12,
        );
        close(
            continued.report().external_heat_j,
            next.boundary_heat_j,
            2e-9,
        );
        close(
            continued.report().material_energy.caloric_heat_j,
            V * next.increment.caloric_heat_j_m3(),
            2e-9,
        );
        close(
            continued.solid().mechanics().element_states()[0].alpha,
            next.increment.mechanical.state.alpha,
            3e-9,
        );
        close(
            continued.report().combined_residual_j,
            0.0,
            config().energy_tolerance_j,
        );
    });
}

