//! G0/G1/G3: original-source porous crush, heat feedback and native transport.
//!
//! The authored coefficients are not measured calibration. Hydrostatic return
//! and caloric closure are solved independently below, including thermal and
//! plastic volume, energetic Voce storage and signed mechanical work defects.

#[path = "finite_porous_enthalpy/pressure.rs"]
mod pressure;

use std::collections::BTreeMap;

use fs_couple::solid_enthalpy::finite_plastic::porous::{
    FiniteThermalPorousCell, FiniteThermalPorousMaterialBinding, FiniteThermalPorousSolid,
};
use fs_couple::solid_enthalpy::finite_plastic::{
    CaloricReference, FiniteThermoplasticConfig, FiniteThermoplasticError, FiniteThermoplasticQuery,
};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy, LawId, LawParameter,
    MaterialCard, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::plastic::finite::porous::thermal::source::THERMAL_FINITE_POROUS_VOCE_LAW_ID;
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
fn source_context(temperature: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("temperature", temperature_quantity(), temperature)
        .unwrap()
}
fn provenance() -> Provenance {
    Provenance {
        source: "authored finite porous thermal regression; not experimental calibration".into(),
        license: "test fixture".into(),
        artifact: None,
    }
}
fn validity() -> ValidityDomain {
    ValidityDomain::unconstrained().with_quantity(
        "temperature",
        temperature_quantity(),
        280.0,
        420.0,
    )
}
fn model() -> ConstitutiveModelCard {
    ConstitutiveModelCard {
        law: LawId(THERMAL_FINITE_POROUS_VOCE_LAW_ID.into()),
        law_version: 1,
        parameters: [
            ("reference_temperature_k", 300.0, K),
            ("minimum_temperature_k", 280.0, K),
            ("maximum_temperature_k", 420.0, K),
            ("shear_modulus_pa", 80e6, PA),
            ("bulk_modulus_pa", 160e6, PA),
            ("yield_stress_pa", 2e6, PA),
            ("hardening_modulus_pa", 8e6, PA),
            ("pressure_sensitivity", 1.2, Dims::NONE),
            ("shear_modulus_slope_pa_k", -20_000.0, PA_K),
            ("bulk_modulus_slope_pa_k", -30_000.0, PA_K),
            ("yield_stress_slope_pa_k", -4_000.0, PA_K),
            ("hardening_modulus_slope_pa_k", -10_000.0, PA_K),
            (
                "logarithmic_expansion_per_k",
                2e-4,
                Dims([0, 0, 0, -1, 0, 0]),
            ),
            (
                "logarithmic_expansion_slope_per_k2",
                1e-7,
                Dims([0, 0, 0, -2, 0, 0]),
            ),
            ("minimum_mechanical_jacobian", 0.1, Dims::NONE),
            ("maximum_elastic_log_strain", 0.8, Dims::NONE),
            ("maximum_plastic_log_strain", 3.0, Dims::NONE),
            ("maximum_accumulated_plastic_strain", 3.0, Dims::NONE),
            ("voce_saturation_pa", 12e6, PA),
            ("voce_rate_per_strain", 3.0, Dims::NONE),
            ("voce_saturation_slope_pa_k", -30_000.0, PA_K),
        ]
        .into_iter()
        .map(|(name, value, dims)| (name.into(), LawParameter { value, dims }))
        .collect(),
        state_schema_version: 1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: validity(),
        sources: vec![fs_blake3::hash_bytes(
            b"authored finite porous heat fixture",
        )],
        provenance: provenance(),
    }
}
fn scalar(name: &str, dims: Dims, value: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: validity(),
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
        let mut claims = ClaimSet::new();
        claims.insert_claim(scalar(CP_NAME, CP_DIMS, cp)).unwrap();
        claims
            .insert_claim(scalar(DENSITY_PROPERTY, Density::DIMS, RHO))
            .unwrap();
        let card = MaterialCard::assemble(
            MaterialStateId {
                chemistry: format!("authored porous thermomechanics {identity}"),
                phase: "solid".into(),
                process: "numerical regression".into(),
                revision: 0,
            },
            claims,
            vec![model()],
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
            1e-6,
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
                kernel_id: 0x46504f54,
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
        relative_tolerance: 1e-11,
        ..Default::default()
    }
}
fn controls() -> FiniteThermoplasticControls {
    FiniteThermoplasticControls {
        material: FinitePlasticStepControls {
            equilibrium_work_tolerance_j: 1e-8,
            maximum_endpoint_work_defect_j: 5.0,
            yield_tolerance_pa: 1e-4,
        },
        maximum_numerical_kinetic_loss_j: 0.01,
    }
}
fn thermal_config() -> StepConfig {
    StepConfig {
        max_dt_s: 0.25,
        energy_tolerance_j: 1e-11,
    }
}
fn config() -> FiniteThermoplasticConfig {
    FiniteThermoplasticConfig {
        maximum_iterations: 100,
        temperature_tolerance_k: 1e-8,
        energy_tolerance_j: 2e-9,
        enthalpy_relaxation: 0.7,
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.16e}, expected={expected:.16e}, tolerance={tolerance:e}"
    );
}
fn prescribed(stretch: f64) -> Vec<(usize, f64)> {
    (0..12)
        .map(|dof| (dof, (stretch - 1.0) * NODES[dof / 3][dof % 3]))
        .collect()
}
fn cell(index: usize) -> DepositedCell {
    DepositedCell {
        id: CellId(100 + index as u64),
        material_id: MaterialId(index as u64),
        segment_id: index as u64,
        body_id: 7 + index as u64,
        mass_kg: RHO * V,
        specific_h_j_kg: 0.0,
        geometry: CellGeometry {
            start_m: [2.0 * L * index as f64, 0.0, 0.0],
            end_m: [2.0 * L * index as f64 + L, L, L],
            volume_m3: V,
            surface_area_m2: 3.0 * L * L,
        },
        boundary: Boundary {
            ambient_k: 300.0,
            htc_w_m2_k: 0.0,
            reservoir: None,
            source_w: 0.0,
        },
    }
}
fn initial<'a>(fixtures: &[&'a Fixture], conductance: Option<f64>) -> DepositionSimulation<'a> {
    let mut events: Vec<_> = (0..fixtures.len())
        .map(|i| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell(i)),
        })
        .collect();
    if let Some(conductance_w_k) = conductance {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 7,
                a: CELL,
                b: CellId(101),
                area_m2: L * L,
                conductance_w_k,
            }),
        });
    }
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
    DepositionSimulation::new(
        events,
        &laws,
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
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
    cx: &Cx<'_>,
) -> FiniteThermalPorousSolid<'a> {
    FiniteThermalPorousSolid::try_new(
        thermal,
        &NODES,
        &TETS,
        vec![FiniteThermalPorousCell {
            cell: CELL,
            query: fixture.query(),
        }],
        TetAssemblyBudget::standard(),
        cx,
    )
    .unwrap()
}

#[derive(Clone, Copy)]
struct HydrostaticEndpoint {
    alpha: f64,
    plastic_jacobian: f64,
    first_piola_pa: f64,
    free_energy_j_m3: f64,
    internal_energy_j_m3: f64,
    plastic_heat_j_m3: f64,
    reversible_heat_j_m3: f64,
    endpoint_work_defect_j_m3: f64,
    work_j_m3: f64,
}
impl HydrostaticEndpoint {
    fn caloric_heat_j_m3(self) -> f64 {
        self.plastic_heat_j_m3 + self.reversible_heat_j_m3
    }
}

/// Scalar hydrostatic oracle from the analytic pressure/volume equation. No
/// material response, finite-element assembly or coupled driver is called.
fn hydrostatic(stretch: f64, temperature: f64) -> HydrostaticEndpoint {
    let delta = temperature - 300.0;
    let log_theta = 2e-4 * delta + 0.5e-7 * delta * delta;
    let theta = log_theta.exp();
    let expansion = 2e-4 + 1e-7 * delta;
    let bulk = 160e6 - 30_000.0 * delta;
    let yield_stress = 2e6 - 4_000.0 * delta;
    let hardening = 8e6 - 10_000.0 * delta;
    let saturation = 12e6 - 30_000.0 * delta;
    let beta = 1.2 / (1.0_f64 + 1.2 * 1.2 / 9.0).sqrt();
    let trial_g = 3.0 * (stretch.ln() - log_theta);
    let radius =
        |alpha: f64| yield_stress + hardening * alpha - saturation * (-3.0 * alpha).exp_m1();
    let (mut lo, mut hi) = (0.0, -trial_g / beta);
    assert!(bulk * beta * -trial_g > radius(0.0));
    for _ in 0..100 {
        let alpha = 0.5 * (lo + hi);
        if bulk * beta * (-trial_g - beta * alpha) > radius(alpha) {
            lo = alpha;
        } else {
            hi = alpha;
        }
    }
    let alpha = 0.5 * (lo + hi);
    let g = trial_g + beta * alpha;
    let piola = bulk * g / stretch;
    let shape = alpha + (-3.0 * alpha).exp_m1() / 3.0;
    let free_energy = 0.5 * bulk * g * g + 0.5 * hardening * alpha * alpha + saturation * shape;
    let entropy =
        3.0 * bulk * expansion * g + 15_000.0 * g * g + 5_000.0 * alpha * alpha + 30_000.0 * shape;
    let plastic_heat = yield_stress * alpha;
    let thermal_work = 3.0 * piola * log_theta.exp_m1();
    HydrostaticEndpoint {
        alpha,
        plastic_jacobian: (-beta * alpha).exp(),
        first_piola_pa: piola,
        free_energy_j_m3: free_energy,
        internal_energy_j_m3: free_energy + temperature * entropy,
        plastic_heat_j_m3: plastic_heat,
        reversible_heat_j_m3: thermal_work - temperature * entropy,
        endpoint_work_defect_j_m3: 3.0 * theta * piola * (stretch / theta - 1.0)
            - free_energy
            - plastic_heat,
        work_j_m3: 3.0 * piola * (stretch - 1.0),
    }
}
fn adiabatic_temperature(stretch: f64, cp: f64) -> f64 {
    let residual = |temperature| {
        RHO * cp * (temperature - 300.0) - hydrostatic(stretch, temperature).caloric_heat_j_m3()
    };
    let (mut lo, mut hi) = (300.0, 420.0);
    assert!(residual(lo) < 0.0 && residual(hi) > 0.0);
    for _ in 0..100 {
        let t = 0.5 * (lo + hi);
        if residual(t) > 0.0 {
            hi = t;
        } else {
            lo = t;
        }
    }
    0.5 * (lo + hi)
}

#[test]
fn closed_adiabatic_crush_solves_free_motion_temperature_and_plastic_volume_together() {
    context(|cx| {
        let fixture = Fixture::new(CP, 0);
        let thermal = initial(&[&fixture], None);
        let solid = bind(&fixture, &thermal, cx);
        let stretch = 0.8;
        let temperature = adiabatic_temperature(stretch, CP);
        let expected = hydrostatic(stretch, temperature);
        let supports: Vec<_> = prescribed(stretch)
            .into_iter()
            .filter(|(dof, _)| *dof != 3)
            .collect();
        // One actual unknown displacement. Independent hydrostatic traction
        // determines the load; feedback must recover both F and T, not just T.
        let force = expected.first_piola_pa * V / L;
        let candidate = solid
            .advance(
                &thermal,
                &(),
                1.0,
                &supports,
                &[(3, force)],
                [0.0; 3],
                settings(),
                controls(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        let actual = candidate.solid().mechanics();
        let report = candidate.report();
        close(
            candidate.thermal().temperature_k(CELL).unwrap(),
            temperature,
            2e-6,
        );
        close(actual.temperatures_k()[0], temperature, 2e-6);
        close(actual.displacement_m()[1][0], (stretch - 1.0) * L, 2e-10);
        for (f, want) in actual.element_deformation()[0]
            .iter()
            .zip(IDENTITY.map(|x| x * stretch))
        {
            close(*f, want, 2e-8);
        }
        let state = &actual.element_states()[0];
        close(state.alpha, expected.alpha, 2e-8);
        close(
            actual.element_plastic_jacobians()[0],
            expected.plastic_jacobian,
            2e-8,
        );
        close(actual.element_jacobians()[0], stretch.powi(3), 2e-8);
        for i in 0..6 {
            close(
                state.inverse_plastic_cauchy_green[i],
                if i < 3 {
                    expected.plastic_jacobian.powf(-2.0 / 3.0)
                } else {
                    0.0
                },
                2e-7,
            );
        }
        close(
            candidate.solid().current_density_kg_m3(cx).unwrap()[0],
            RHO / stretch.powi(3),
            1e-4,
        );
        let bound = &actual.families()[0];
        let family = bound.family().law();
        let isothermal = family
            .at_temperature(300.0)
            .unwrap()
            .response(&IDENTITY.map(|x| x * stretch), &family.initial_state())
            .unwrap();
        assert!((state.alpha - isothermal.state.alpha).abs() > 1e-4);
        assert!(temperature > 301.0 && report.iterations > 1);
        assert!(std::ptr::eq(bound.caloric_source(), &fixture.caloric));
        assert_eq!(bound.family().model_card(), &fixture.card.models()[0]);
        close(bound.reference_density_kg_m3(), RHO, 0.0);
        close(candidate.thermal().cells()[&CELL].mass_kg, RHO * V, 0.0);
        close(actual.free_energy_j(), V * expected.free_energy_j_m3, 2e-7);
        close(
            actual.internal_energy_j(),
            V * expected.internal_energy_j_m3,
            2e-7,
        );
        close(
            report.material_energy.plastic_dissipation_j,
            V * expected.plastic_heat_j_m3,
            2e-8,
        );
        close(
            report.material_energy.reversible_heat_j,
            V * expected.reversible_heat_j_m3,
            2e-7,
        );
        close(
            report.material_energy.endpoint_work_defect_j,
            V * expected.endpoint_work_defect_j_m3,
            2e-7,
        );
        close(
            report.material_energy.internal_endpoint_work_j,
            V * expected.work_j_m3,
            2e-7,
        );
        close(report.external_heat_j, 0.0, 2e-9);
        let stored = RHO * V * CP * (temperature - 300.0);
        close(report.thermal_audit.stored_energy_change_j, stored, 2e-7);
        close(stored, V * expected.caloric_heat_j_m3(), 2e-9);
        // Endpoint force quadrature error is substantial here and is never heat.
        assert!(report.material_energy.endpoint_work_defect_j.abs() > 0.1);
        close(
            V * expected.work_j_m3,
            actual.internal_energy_j() + stored + report.material_energy.endpoint_work_defect_j,
            3e-7,
        );
        close(report.thermal_residual_j, 0.0, config().energy_tolerance_j);
        close(report.combined_residual_j, 0.0, config().energy_tolerance_j);
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(solid.mechanics().element_states()[0].alpha, 0.0);
    });
}

#[test]
fn two_original_cells_conduct_only_physical_caloric_transfer_and_retain_local_compaction() {
    context(|cx| {
        let first = Fixture::new(CP, 0);
        let second = Fixture::new(700.0, 1);
        let fixtures = [&first, &second];
        let stretches = [0.82, 0.94];
        let capacities = [RHO * V * CP, RHO * V * 700.0];
        let conductance = 0.01;
        // Independent two-capacity explicit recurrence for the native four
        // substeps. Constant source powers use only physical endpoint heat.
        let mut temperatures = [300.0; 2];
        let mut last_change = f64::INFINITY;
        for _ in 0..100 {
            let powers: [f64; 2] = core::array::from_fn(|i| {
                V * hydrostatic(stretches[i], temperatures[i]).caloric_heat_j_m3()
            });
            let mut next = [300.0; 2];
            for _ in 0..4 {
                let flow = conductance * (next[1] - next[0]);
                next = [
                    next[0] + 0.25 * (powers[0] + flow) / capacities[0],
                    next[1] + 0.25 * (powers[1] - flow) / capacities[1],
                ];
            }
            last_change = (next[0] - temperatures[0])
                .abs()
                .max((next[1] - temperatures[1]).abs());
            temperatures = next;
        }
        assert!(last_change < 1e-10);
        let endpoints: [HydrostaticEndpoint; 2] =
            core::array::from_fn(|i| hydrostatic(stretches[i], temperatures[i]));
        let nodes: [[f64; 3]; 8] = core::array::from_fn(|i| {
            let mut node = NODES[i % 4];
            if i >= 4 {
                node[0] += 2.0 * L;
            }
            node
        });
        let tets = [[0, 1, 2, 3], [4, 5, 6, 7]];
        let supports: Vec<_> = (0..24)
            .map(|dof| {
                let node = dof / 3;
                (dof, (stretches[node / 4] - 1.0) * NODES[node % 4][dof % 3])
            })
            .collect();
        let thermal = initial(&fixtures, Some(conductance));
        let solid = FiniteThermalPorousSolid::try_new(
            &thermal,
            &nodes,
            &tets,
            fixtures
                .iter()
                .enumerate()
                .map(|(i, fixture)| FiniteThermalPorousCell {
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
                &supports,
                &[],
                [0.0; 3],
                settings(),
                controls(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        let mut link_heat = [0.0; 2];
        let mut total_stored = 0.0;
        for i in 0..2 {
            let id = CellId(100 + i as u64);
            let cell = &candidate.thermal().cells()[&id];
            close(
                candidate.thermal().temperature_k(id).unwrap(),
                temperatures[i],
                2e-6,
            );
            close(
                candidate.solid().mechanics().temperatures_k()[i],
                temperatures[i],
                2e-6,
            );
            close(
                candidate.solid().mechanics().element_states()[i].alpha,
                endpoints[i].alpha,
                2e-8,
            );
            close(
                candidate.solid().mechanics().element_plastic_jacobians()[i],
                endpoints[i].plastic_jacobian,
                2e-8,
            );
            close(cell.mass_kg, RHO * V, 0.0);
            assert!(std::ptr::eq(
                candidate.solid().mechanics().families()[i].caloric_source(),
                &fixtures[i].caloric
            ));
            let stored = cell.mass_kg * cell.specific_h_j_kg;
            total_stored += stored;
            link_heat[i] = stored - V * endpoints[i].caloric_heat_j_m3();
            close(
                link_heat[i],
                capacities[i] * (temperatures[i] - 300.0) - V * endpoints[i].caloric_heat_j_m3(),
                2e-7,
            );
        }
        assert!(link_heat[0] < -1e-4 && link_heat[1] > 1e-4);
        close(link_heat.iter().sum(), 0.0, 2e-8);
        let report = candidate.report();
        let physical_heat = V * endpoints.iter().map(|x| x.caloric_heat_j_m3()).sum::<f64>();
        close(total_stored, physical_heat, 2e-8);
        close(report.material_energy.caloric_heat_j, physical_heat, 2e-8);
        close(report.external_heat_j, 0.0, 2e-9);
        close(
            report.material_energy.endpoint_work_defect_j,
            V * endpoints
                .iter()
                .map(|x| x.endpoint_work_defect_j_m3)
                .sum::<f64>(),
            2e-7,
        );
        assert!(report.material_energy.endpoint_work_defect_j.abs() > 0.1);
        close(report.combined_residual_j, 0.0, config().energy_tolerance_j);
        close(report.thermal_residual_j, 0.0, config().energy_tolerance_j);
        assert_eq!(candidate.thermal().links(), thermal.links());
        assert_eq!(thermal.time_s(), 0.0);
        assert!(
            solid
                .mechanics()
                .element_states()
                .iter()
                .all(|state| state.alpha == 0.0)
        );
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
fn late_observer_and_iteration_refusals_preserve_hot_compacted_motion_and_retry() {
    context(|cx| {
        let fixture = Fixture::new(CP, 0);
        let thermal = initial(&[&fixture], None);
        let solid = bind(&fixture, &thermal, cx)
            .with_inertia(&[[0.01, -0.005, 0.003]; 4], cx)
            .unwrap();
        let first = solid
            .advance(
                &thermal,
                &(),
                1.0,
                &prescribed(0.9),
                &[],
                [0.0; 3],
                settings(),
                controls(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        let (thermal, solid, (), _) = first.into_parts();
        let old_cells = thermal.cells().clone();
        let old_history = solid.mechanics().element_states().to_vec();
        let old_u = solid.mechanics().displacement_m().to_vec();
        let old_v = solid.mechanics().velocity_mps().unwrap().to_vec();
        let old_f = solid.mechanics().element_deformation().to_vec();
        let old_p = solid.mechanics().element_first_piola_pa().to_vec();
        let old_jp = solid.mechanics().element_plastic_jacobians().to_vec();
        let old_temperature = solid.mechanics().temperatures_k().to_vec();
        let old_internal = solid.mechanics().internal_energy_j();
        let old_free = solid.mechanics().free_energy_j();
        assert!(old_history[0].alpha > 0.2 && old_jp[0] < 0.8 && old_temperature[0] > 301.0);
        let supports = prescribed(0.82);
        let observer = LateRefusal {
            completed_substeps: 0,
        };
        let late = solid.advance(
            &thermal,
            &observer,
            2.0,
            &supports,
            &[],
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
        let short = solid.advance(
            &thermal,
            &(),
            2.0,
            &supports,
            &[],
            [0.0; 3],
            settings(),
            controls(),
            thermal_config(),
            FiniteThermoplasticConfig {
                maximum_iterations: 1,
                ..config()
            },
            cx,
        );
        assert!(matches!(
            short,
            Err(FiniteThermoplasticError::NotConverged { iterations: 1, .. })
        ));
        assert_eq!(thermal.cells(), &old_cells);
        assert_eq!(solid.mechanics().element_states(), old_history);
        assert_eq!(solid.mechanics().displacement_m(), old_u);
        assert_eq!(solid.mechanics().velocity_mps().unwrap(), old_v);
        assert_eq!(solid.mechanics().element_deformation(), old_f);
        assert_eq!(solid.mechanics().element_first_piola_pa(), old_p);
        assert_eq!(solid.mechanics().element_plastic_jacobians(), old_jp);
        assert_eq!(solid.mechanics().temperatures_k(), old_temperature);
        assert_eq!(solid.mechanics().internal_energy_j(), old_internal);
        assert_eq!(solid.mechanics().free_energy_j(), old_free);
        assert_eq!(thermal.time_s(), 1.0);
        assert_eq!(solid.time_s(), 1.0);
        close(
            solid.mechanics().reference_mass_kg().unwrap(),
            RHO * V,
            1e-18,
        );
        let retried = solid
            .advance(
                &thermal,
                &(),
                2.0,
                &supports,
                &[],
                [0.0; 3],
                settings(),
                controls(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        assert!(retried.solid().mechanics().element_states()[0].alpha > old_history[0].alpha);
        assert!(retried.solid().mechanics().element_plastic_jacobians()[0] < old_jp[0]);
        let inertia = retried.report().inertia.as_ref().unwrap();
        close(inertia.reference_mass_kg, RHO * V, 1e-18);
        assert!(inertia.numerical_loss_j > 0.0);
        close(retried.thermal().cells()[&CELL].mass_kg, RHO * V, 0.0);
        close(
            retried.report().thermal_audit.stored_energy_change_j,
            retried.report().material_energy.caloric_heat_j,
            2e-9,
        );
        close(
            retried.report().combined_residual_j,
            0.0,
            config().energy_tolerance_j,
        );
        assert_eq!(thermal.cells(), &old_cells);
        assert_eq!(solid.mechanics().element_states(), old_history);
        assert_eq!(solid.mechanics().velocity_mps().unwrap(), old_v);
    });
}

#[test]
fn source_caloric_capacity_exact_object_and_reference_inventory_are_admitted_together() {
    context(|cx| {
        let fixture = Fixture::new(CP, 0);
        let admitted = FiniteThermalPorousMaterialBinding::try_new(fixture.query(), 0.0).unwrap();
        assert!(admitted.minimum_total_heat_capacity_j_m3_k() > 0.0);
        let inadequate = Fixture::new(1.0, 0);
        let refusal =
            FiniteThermalPorousMaterialBinding::try_new(inadequate.query(), 0.0).unwrap_err();
        assert!(refusal.to_string().contains("heat capacity"));
        let foreign = Fixture::new(CP, 1);
        let mut query = fixture.query();
        query.caloric = &foreign.caloric;
        assert!(FiniteThermalPorousMaterialBinding::try_new(query, 0.0).is_err());
        let mut pressured = fixture.query();
        pressured.context = pressured
            .context
            .with_quantity(
                "environment_pressure",
                QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::Pressure,
                    ValueForm::Static,
                )),
                1e5,
            )
            .unwrap();
        let refusal = FiniteThermalPorousMaterialBinding::try_new(pressured, 0.0).unwrap_err();
        assert!(
            refusal
                .to_string()
                .contains("zero-pressure caloric reference")
        );
        assert!(FiniteThermalPorousMaterialBinding::try_new(fixture.query(), CP).is_err());
        let thermal = initial(&[&fixture], None);
        let duplicate = fixture.caloric.clone();
        let mut replaced = fixture.query();
        replaced.caloric = &duplicate;
        assert!(FiniteThermalPorousMaterialBinding::try_new(replaced.clone(), 0.0).is_ok());
        assert!(
            FiniteThermalPorousSolid::try_new(
                &thermal,
                &NODES,
                &TETS,
                vec![FiniteThermalPorousCell {
                    cell: CELL,
                    query: replaced
                }],
                TetAssemblyBudget::standard(),
                cx,
            )
            .is_err()
        );
        let mut wrong_mass = cell(0);
        wrong_mass.mass_kg *= 0.99;
        let inconsistent = DepositionSimulation::new(
            vec![TimedEvent {
                time_s: 0.0,
                event: Event::Birth(wrong_mass),
            }],
            &[(MaterialId(0), fixture.caloric.curve() as &dyn ThermalLaw)],
            DepositionBudget {
                max_cells: 1,
                max_links: 0,
                max_events: 1,
                max_substeps: 100,
            },
        )
        .unwrap()
        .advance_to(0.0, thermal_config(), || Ok(()))
        .unwrap()
        .0;
        assert!(
            FiniteThermalPorousSolid::try_new(
                &inconsistent,
                &NODES,
                &TETS,
                vec![FiniteThermalPorousCell {
                    cell: CELL,
                    query: fixture.query()
                }],
                TetAssemblyBudget::standard(),
                cx,
            )
            .is_err()
        );
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(thermal.cells()[&CELL].specific_h_j_kg, 0.0);
    });
}
