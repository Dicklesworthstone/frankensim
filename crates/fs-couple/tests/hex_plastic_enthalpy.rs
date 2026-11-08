//! G1/G3 coupled Hex8 deformation, heat and retained history regressions.
//! Values are authored physical fixtures, not measured material calibration.

use std::collections::BTreeMap;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_couple::iqn_ils::{IqnIlsConfig, driver::CouplingMethod};
use fs_couple::solid_enthalpy::hex::{
    HexThermoplasticConfig, HexThermoplasticDiscretization, HexThermoplasticError,
    ThermoplasticHexSolid,
};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::hill::source::{
    HILL_COEFFICIENT_PROPERTIES, HILL_HARDENING_MODULUS_PROPERTY,
    HILL_OVERSTRESS_VISCOSITY_PROPERTY, HILL_REFERENCE_YIELD_STRESS_PROPERTY,
};
use fs_material::plastic::hill::thermal::{ThermalHillSource, ThermalHillSourceQuery};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, POISSON_RATIO_PROPERTY, YOUNG_MODULUS_PROPERTY,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature};
use fs_solid::Formulation;
use fs_solid::plastic3::IncrementSettings;
use fs_solid::plastic3::hex::HexAssemblyBudget;
use fs_solid::plastic3::hex::thermal::ThermalHexTolerances;
use fs_thermal::enthalpy::*;

const L: f64 = 0.01;
const V: f64 = L * L * L;
const E: f64 = 1e9;
const NU: f64 = 0.3;
const Y: f64 = 1e6;
const H: f64 = 5e7;
const SOFTENING: f64 = 1e5;
const ETA: f64 = 5e6;
const RHO: f64 = 1000.0;
const CP: f64 = 20.0;
const DT: f64 = 0.1;
const NODES: [[f64; 3]; 8] = [
    [0.0, 0.0, 0.0],
    [L, 0.0, 0.0],
    [L, L, 0.0],
    [0.0, L, 0.0],
    [0.0, 0.0, L],
    [L, 0.0, L],
    [L, L, L],
    [0.0, L, L],
];
const HEX: [[usize; 8]; 1] = [[0, 1, 2, 3, 4, 5, 6, 7]];
const CP_NAME: &str = "specific_heat_capacity";
const CP_DIMS: Dims = Dims([2, 0, -2, -1, 0, 0]);

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}; tolerance {tolerance:e}"
    );
}
fn context(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 17,
                kernel_id: 0x48455854,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}
fn temperature_spec() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", temperature_spec(), t)
        .unwrap()
}
fn property(name: &str, dims: Dims, value: f64, slope: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: if slope == 0.0 {
            PropertyValue::Scalar { value, dims }
        } else {
            PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(280.0, value - 20.0 * slope), (308.0, value + 8.0 * slope)],
                dims,
            }
        },
        validity: ValidityDomain::unconstrained().with_quantity(
            "T",
            temperature_spec(),
            280.0,
            308.0,
        ),
        interpolation: if slope == 0.0 {
            InterpolationPolicy::ConstantWithinValidity
        } else {
            InterpolationPolicy::LinearInside
        },
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "authored anisotropic Hex8 thermal feedback".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}
struct Fixture {
    card: MaterialCard,
    caloric: ResolvedHeatCapacityCurve,
}
impl Fixture {
    fn new() -> Self {
        let mut claims = ClaimSet::new();
        for (name, dims, value, slope) in [
            (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, E, 0.0),
            (POISSON_RATIO_PROPERTY, Dims::NONE, NU, 0.0),
            (
                HILL_REFERENCE_YIELD_STRESS_PROPERTY,
                Pressure::DIMS,
                Y,
                -SOFTENING,
            ),
            (HILL_HARDENING_MODULUS_PROPERTY, Pressure::DIMS, H, 0.0),
            (
                HILL_OVERSTRESS_VISCOSITY_PROPERTY,
                Dims([-1, 1, -1, 0, 0, 0]),
                ETA,
                0.0,
            ),
            (DENSITY_PROPERTY, Density::DIMS, RHO, 0.0),
            (CP_NAME, CP_DIMS, CP, 0.0),
        ] {
            claims
                .insert_claim(property(name, dims, value, slope))
                .unwrap();
        }
        for (name, value) in HILL_COEFFICIENT_PROPERTIES
            .into_iter()
            .zip([0.35, 0.8, 0.45, 1.5, 1.1, 2.0])
        {
            claims
                .insert_claim(property(name, Dims::NONE, value, 0.0))
                .unwrap();
        }
        let card = MaterialCard::assemble(
            MaterialStateId {
                chemistry: "authored rolled metal".into(),
                phase: "solid".into(),
                process: "coupled numerical fixture".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap();
        let caloric = ResolvedHeatCapacityCurve::try_new(
            &card,
            HeatCapacitySourceQuery {
                heat_capacity: &PropertyKey::new(CP_NAME, CP_DIMS),
                density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
                temperature_axis: "T",
                lower: &point(280.0),
                upper: &point(308.0),
                selection: SelectionPolicy::SingleClaimOnly,
            },
            SolidLiquidPhase::Solid,
            -20.0 * CP,
            1e-8,
        )
        .unwrap();
        Self { card, caloric }
    }
    fn source(&self, rate: bool) -> ThermalHillSource<'_> {
        let source = ThermalHillSource::try_new(
            ThermalHillSourceQuery {
                card: &self.card,
                caloric: &self.caloric,
                context: point(300.0),
                temperature_axis: "T".into(),
                selection: MaterialPropertySelection::SingleClaimOnly,
                thermal_expansion: None,
                strain_limit: 0.04,
                material_to_spatial: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            },
            0.0,
        )
        .unwrap();
        if rate {
            source
                .with_linear_overstress(MaterialPropertySelection::SingleClaimOnly)
                .unwrap()
        } else {
            source
        }
    }
}
fn cell(point: usize) -> CellId {
    CellId(100 + point as u64)
}
fn thermal_config() -> StepConfig {
    StepConfig {
        max_dt_s: DT / 2.0,
        energy_tolerance_j: 1e-10,
    }
}
fn config() -> HexThermoplasticConfig {
    HexThermoplasticConfig {
        maximum_iterations: 80,
        temperature_tolerance_k: 1e-8,
        energy_tolerance_j: 1e-8,
        enthalpy_relaxation: 0.8,
    }
}
fn tolerances() -> ThermalHexTolerances {
    ThermalHexTolerances {
        material: EnergyTolerances {
            energy_j_m3: 1e-4,
            yield_pa: 1e-4,
        },
        energy_j: 1e-8,
    }
}
fn settings() -> IncrementSettings {
    IncrementSettings {
        relative_tolerance: 1e-11,
        linear_tolerance: 1e-12,
        max_newton_iterations: 60,
        ..IncrementSettings::default()
    }
}
fn discretization(formulation: Formulation) -> HexThermoplasticDiscretization {
    HexThermoplasticDiscretization {
        formulation,
        strain_limit: 0.04,
        rotation_limit_rad: 0.04,
        budget: HexAssemblyBudget::default(),
    }
}
fn thermal<'a>(
    fixture: &'a Fixture,
    volumes: &[f64; 8],
    external_rise_k: f64,
    conductance: f64,
) -> DepositionSimulation<'a> {
    let mut events: Vec<_> = volumes
        .iter()
        .enumerate()
        .map(|(i, &volume)| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: cell(i),
                material_id: MaterialId(0),
                segment_id: i as u64,
                body_id: 1,
                mass_kg: RHO * volume,
                specific_h_j_kg: 0.0,
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [L; 3],
                    volume_m3: volume,
                    surface_area_m2: L * L,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: RHO * volume * CP * external_rise_k / DT,
                },
            }),
        })
        .collect();
    if conductance > 0.0 {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 1,
                a: cell(0),
                b: cell(7),
                area_m2: L * L / 4.0,
                conductance_w_k: conductance,
            }),
        });
    }
    DepositionSimulation::new(
        events,
        &[(MaterialId(0), fixture.caloric.curve() as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 8,
            max_links: 1,
            max_events: 9,
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
    nodes: &'a [[f64; 3]],
    formulation: Formulation,
    rate: bool,
    cx: &Cx<'_>,
) -> ThermoplasticHexSolid<'a, ThermalHillSource<'a>> {
    ThermoplasticHexSolid::try_new(
        thermal,
        nodes,
        &HEX,
        vec![fixture.source(rate)],
        &(0..8).map(cell).collect::<Vec<_>>(),
        discretization(formulation),
        cx,
    )
    .unwrap()
}
fn symmetry() -> Vec<(usize, f64)> {
    (0..24)
        .filter(|d| NODES[d / 3][d % 3] == 0.0)
        .map(|d| (d, 0.0))
        .collect()
}

// Independent uniaxial associated Hill flow, coupled to one caloric scalar.
// q=sqrt(G+H_hill)*sigma; ep=alpha*[r,-H_hill/r,-G/r].
fn uniaxial_endpoint(old_t: f64, old_alpha: f64, sigma: f64) -> (f64, f64) {
    let q = 1.25_f64.sqrt() * sigma;
    let residual = |t: f64| {
        let y = Y - SOFTENING * (t - 300.0);
        let alpha = old_alpha.max((q - y) / H);
        RHO * CP * (t - old_t - 0.5) - y * (alpha - old_alpha)
    };
    let (mut lo, mut hi) = (old_t + 0.5, 308.0);
    assert!(residual(lo) <= 0.0 && residual(hi) >= 0.0);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if residual(mid) > 0.0 {
            hi = mid
        } else {
            lo = mid
        }
    }
    let t = 0.5 * (lo + hi);
    (t, old_alpha.max((q - Y + SOFTENING * (t - 300.0)) / H))
}

#[test]
fn force_held_heating_softens_then_unloads_to_permanent_anisotropic_set() {
    context(|cx| {
        let fixture = Fixture::new();
        for formulation in [Formulation::Standard, Formulation::BBar] {
            for method in [
                CouplingMethod::RelaxedPicard,
                CouplingMethod::IqnIls(IqnIlsConfig::default()),
            ] {
                let mut heat = thermal(&fixture, &[V / 8.0; 8], 0.5, 0.0);
                let mut solid = bind(&fixture, &heat, &NODES, formulation, false, cx);
                let fixed = symmetry();
                let (mut old_t, mut old_alpha) = (300.0, 0.0);
                for (interval, sigma) in [1.8e6, 1.8e6, 0.0].into_iter().enumerate() {
                    let forces: Vec<_> = (0..8)
                        .filter(|n| NODES[*n][0] == L)
                        .map(|n| (3 * n, sigma * L * L / 4.0))
                        .collect();
                    let old_solution = solid.solution().clone();
                    let old_cells = heat.cells().clone();
                    let target = (interval + 1) as f64 * DT;
                    if interval == 0 {
                        let refused = solid.advance(
                            &heat,
                            &(),
                            target,
                            &fixed,
                            &forces,
                            settings(),
                            tolerances(),
                            thermal_config(),
                            HexThermoplasticConfig {
                                maximum_iterations: 1,
                                ..config()
                            },
                            method,
                            cx,
                        );
                        assert!(matches!(
                            refused,
                            Err(HexThermoplasticError::NotConverged { .. })
                        ));
                        assert_eq!(solid.solution(), &old_solution);
                        assert_eq!(heat.cells(), &old_cells);
                    }
                    let step = solid
                        .advance(
                            &heat,
                            &(),
                            target,
                            &fixed,
                            &forces,
                            settings(),
                            tolerances(),
                            thermal_config(),
                            config(),
                            method,
                            cx,
                        )
                        .unwrap();
                    let (t, alpha) = uniaxial_endpoint(old_t, old_alpha, sigma);
                    let r = 1.25_f64.sqrt();
                    let ep = [r * alpha, -0.45 / r * alpha, -0.8 / r * alpha];
                    for i in 0..8 {
                        near(step.thermal().temperature_k(cell(i)).unwrap(), t, 5e-7);
                        near(
                            step.solid().solution().quadrature_states[i].alpha,
                            alpha,
                            2e-9,
                        );
                        for j in 0..3 {
                            let elastic = if j == 0 { sigma / E } else { -NU * sigma / E };
                            near(
                                step.solid().solution().displacement_m[i][j],
                                (elastic + ep[j]) * NODES[i][j],
                                3e-11,
                            );
                        }
                    }
                    let energy = &step.report().material_energy.material;
                    near(
                        energy.physical_heat_j,
                        V * (Y - SOFTENING * (t - 300.0)) * (alpha - old_alpha),
                        3e-9,
                    );
                    near(energy.reversible_heat_j, 0.0, 1e-12);
                    near(step.report().external_heat_j, RHO * V * CP * 0.5, 3e-9);
                    near(
                        step.report().combined_residual_j,
                        0.0,
                        config().energy_tolerance_j,
                    );
                    near(
                        step.report().thermal_residual_j,
                        0.0,
                        config().energy_tolerance_j,
                    );
                    assert_eq!(solid.solution(), &old_solution);
                    assert_eq!(heat.cells(), &old_cells);
                    if interval == 1 {
                        assert!(alpha > old_alpha);
                    }
                    if interval == 2 {
                        near(alpha, old_alpha, 1e-14);
                        near(energy.physical_heat_j, 0.0, 1e-12);
                        assert!(step.solid().solution().displacement_m[1][0] > 0.0);
                    }
                    let (next_heat, next_solid, _, _, _) = step.into_parts();
                    heat = next_heat;
                    solid = next_solid;
                    (old_t, old_alpha) = (t, alpha);
                }
            }
        }
    });
}

// Independent pure-shear rate return with its actual dt and endpoint heat.
fn rate_endpoint(old_t: f64, old_alpha: f64, gamma: f64, dt: f64) -> (f64, f64, f64) {
    let mu = E / (2.0 * (1.0 + NU));
    let r = 2.0;
    let alpha = |t: f64| {
        old_alpha.max(
            (r * mu * gamma - Y + SOFTENING * (t - 300.0) + ETA / dt * old_alpha)
                / (r * r * mu + H + ETA / dt),
        )
    };
    let residual = |t: f64| {
        let da = alpha(t) - old_alpha;
        RHO * CP * (t - old_t) - (Y - SOFTENING * (t - 300.0)) * da - ETA / dt * da * da
    };
    let (mut lo, mut hi) = (old_t, 308.0);
    assert!(residual(lo) < 0.0 && residual(hi) > 0.0);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if residual(mid) > 0.0 {
            hi = mid
        } else {
            lo = mid
        }
    }
    let t = 0.5 * (lo + hi);
    (t, alpha(t), mu * (gamma - r * alpha(t)))
}

#[test]
fn self_heated_rate_load_and_holds_consume_the_original_history_and_clock_once() {
    context(|cx| {
        let fixture = Fixture::new();
        let gamma = 0.008;
        let fixed: Vec<_> = NODES
            .iter()
            .flat_map(|x| [gamma * x[1], 0.0, 0.0])
            .enumerate()
            .collect();
        let mut heat = thermal(&fixture, &[V / 8.0; 8], 0.0, 0.0);
        let mut solid = bind(&fixture, &heat, &NODES, Formulation::BBar, true, cx);
        let (mut old_t, mut old_alpha, mut time) = (300.0, 0.0, 0.0);
        for (interval, dt) in [0.05, 0.2, 0.1].into_iter().enumerate() {
            time += dt;
            let (t, alpha, tau) = rate_endpoint(old_t, old_alpha, gamma, dt);
            let step = solid
                .advance(
                    &heat,
                    &(),
                    time,
                    &fixed,
                    &[],
                    settings(),
                    tolerances(),
                    thermal_config(),
                    config(),
                    CouplingMethod::RelaxedPicard,
                    cx,
                )
                .unwrap();
            assert!(step.report().iterations > 1);
            for i in 0..8 {
                near(step.thermal().temperature_k(cell(i)).unwrap(), t, 5e-7);
                near(
                    step.solid().solution().quadrature_states[i].alpha,
                    alpha,
                    2e-10,
                );
                near(
                    step.solid().solution().quadrature_stress_pa[i][3],
                    tau,
                    0.03,
                );
            }
            let da = alpha - old_alpha;
            let expected_heat = V * ((Y - SOFTENING * (t - 300.0)) * da + ETA / dt * da * da);
            near(
                step.report().material_energy.material.physical_heat_j,
                expected_heat,
                2e-10,
            );
            near(
                step.report().thermal_audit.stored_energy_change_j,
                expected_heat,
                2e-10,
            );
            near(step.report().external_heat_j, 0.0, 1e-12);
            near(
                step.report().combined_residual_j,
                0.0,
                config().energy_tolerance_j,
            );
            if interval > 0 {
                near(
                    step.report().material_energy.applied_work_j
                        + step.report().material_energy.support_work_j,
                    0.0,
                    1e-12,
                );
                assert!(alpha > old_alpha && t > old_t);
            }
            let (next_heat, next_solid, _, _, _) = step.into_parts();
            heat = next_heat;
            solid = next_solid;
            (old_t, old_alpha) = (t, alpha);
        }
    });
}

// Independent det(dX/dxi) of X=L*s*(1+.2t), Y=L*t*(1+.1u),
// Z=L*u*(1+.15s), evaluated at the eight Gauss points.
fn warped_volumes() -> [f64; 8] {
    let g = 1.0 / 3.0_f64.sqrt();
    core::array::from_fn(|point| {
        let s = 0.5 * (1.0 + if point / 4 == 0 { -g } else { g });
        let t = 0.5 * (1.0 + if point / 2 % 2 == 0 { -g } else { g });
        let u = 0.5 * (1.0 + if point % 2 == 0 { -g } else { g });
        V / 8.0
            * ((1.0 + 0.2 * t) * (1.0 + 0.1 * u) * (1.0 + 0.15 * s) + 0.2 * 0.1 * 0.15 * s * t * u)
    })
}

#[test]
fn warped_bbar_points_keep_distinct_inventory_and_conduct_generated_heat() {
    context(|cx| {
        let fixture = Fixture::new();
        let nodes = NODES.map(|[x, y, z]| {
            [
                x * (1.0 + 0.2 * y / L),
                y * (1.0 + 0.1 * z / L),
                z * (1.0 + 0.15 * x / L),
            ]
        });
        let volumes = warped_volumes();
        let conductance = 0.002;
        let heat = thermal(&fixture, &volumes, 0.0, conductance);
        let solid = bind(&fixture, &heat, &nodes, Formulation::BBar, false, cx);
        let fixed: Vec<_> = nodes
            .iter()
            .flat_map(|x| [0.008 * x[1] + 0.003 * x[1] * x[2] / L, -0.001 * x[0], 0.0])
            .enumerate()
            .collect();
        let step = solid
            .advance(
                &heat,
                &(),
                DT,
                &fixed,
                &[],
                settings(),
                tolerances(),
                thermal_config(),
                config(),
                CouplingMethod::RelaxedPicard,
                cx,
            )
            .unwrap();
        let generated: Vec<_> = step
            .quadrature_energy()
            .iter()
            .map(|e| e.physical_heat_j + e.reversible_heat_j)
            .collect();
        let mut expected = [300.0; 8];
        // The native explicit timeline takes two .05 s substeps here. Each
        // converged mechanical heat is one constant source on that interval.
        for _ in 0..2 {
            let mut rate: Vec<_> = generated.iter().map(|q| q / DT).collect();
            let transfer = conductance * (expected[7] - expected[0]);
            rate[0] += transfer;
            rate[7] -= transfer;
            for i in 0..8 {
                expected[i] += DT / 2.0 * rate[i] / (RHO * volumes[i] * CP);
            }
        }
        let alpha0 = step.solid().solution().quadrature_states[0].alpha;
        assert!(
            step.solid()
                .solution()
                .quadrature_states
                .iter()
                .any(|state| (state.alpha - alpha0).abs() > 1e-5)
        );
        for i in 0..8 {
            near(
                step.solid().solution().quadrature_volume_m3[i],
                volumes[i],
                2e-21,
            );
            near(
                step.thermal().temperature_k(cell(i)).unwrap(),
                expected[i],
                5e-7,
            );
        }
        near(
            step.report().thermal_audit.stored_energy_change_j,
            generated.iter().sum(),
            3e-9,
        );
        near(
            step.report().combined_residual_j,
            0.0,
            config().energy_tolerance_j,
        );
        assert!((expected[0] - (300.0 + generated[0] / (RHO * volumes[0] * CP))).abs() > 1e-5);

        // Swapping unequal point inventories keeps total mass unchanged but
        // must not bind to the wrong native material-point measures.
        let mut wrong = volumes;
        wrong.swap(0, 7);
        let wrong_heat = thermal(&fixture, &wrong, 0.0, 0.0);
        assert!(
            ThermoplasticHexSolid::try_new(
                &wrong_heat,
                &nodes,
                &HEX,
                vec![fixture.source(false)],
                &(0..8).map(cell).collect::<Vec<_>>(),
                discretization(Formulation::BBar),
                cx,
            )
            .is_err()
        );
        let foreign = Fixture::new();
        assert!(
            ThermoplasticHexSolid::try_new(
                &heat,
                &nodes,
                &HEX,
                vec![foreign.source(false)],
                &(0..8).map(cell).collect::<Vec<_>>(),
                discretization(Formulation::BBar),
                cx,
            )
            .is_err()
        );
    });
}

#[derive(Clone, Debug)]
struct Probe {
    steps: usize,
    fail: bool,
}
impl DepositionObserver for Probe {
    type Error = &'static str;
    fn step(
        &mut self,
        start: f64,
        dt: f64,
        _old: &BTreeMap<CellId, f64>,
        _thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        checkpoint().map_err(|_| "cancelled")?;
        self.steps += 1;
        if self.fail && start + dt > DT / 2.0 {
            Err("late caloric observer refusal")
        } else {
            Ok(())
        }
    }
}

#[test]
fn consistent_inertial_stopping_stays_out_of_heat_and_late_failure_preserves_all_owners() {
    context(|cx| {
        let fixture = Fixture::new();
        let heat = thermal(&fixture, &[V / 8.0; 8], 0.0, 0.0);
        let velocity = [[0.2, -0.1, 0.05]; 8];
        let solid = bind(&fixture, &heat, &NODES, Formulation::Standard, true, cx)
            .with_initial_velocity(&velocity, cx)
            .unwrap();
        let fixed: Vec<_> = (0..24).map(|d| (d, 0.0)).collect();
        let observer = Probe {
            steps: 0,
            fail: true,
        };
        let old_solution = solid.solution().clone();
        let old_cells = heat.cells().clone();
        let refused = solid.advance(
            &heat,
            &observer,
            DT,
            &fixed,
            &[],
            settings(),
            tolerances(),
            thermal_config(),
            config(),
            CouplingMethod::RelaxedPicard,
            cx,
        );
        assert!(matches!(
            refused,
            Err(HexThermoplasticError::Observer(
                "late caloric observer refusal"
            ))
        ));
        assert_eq!(observer.steps, 0);
        assert_eq!(solid.solution(), &old_solution);
        assert_eq!(solid.velocity_mps().unwrap(), &velocity);
        assert_eq!(heat.cells(), &old_cells);
        near(solid.time_s(), 0.0, 0.0);
        let accepted_observer = Probe {
            fail: false,
            ..observer.clone()
        };
        let step = solid
            .advance(
                &heat,
                &accepted_observer,
                DT,
                &fixed,
                &[],
                settings(),
                tolerances(),
                thermal_config(),
                config(),
                CouplingMethod::RelaxedPicard,
                cx,
            )
            .unwrap();
        let kinetic = 0.5 * RHO * V * (0.2_f64.powi(2) + 0.1_f64.powi(2) + 0.05_f64.powi(2));
        let inertia = step.report().inertia.unwrap();
        near(inertia.reference_mass_kg, RHO * V, 2e-17);
        near(inertia.kinetic_change_j, -kinetic, 2e-18);
        near(inertia.numerical_loss_j, kinetic, 2e-18);
        near(
            step.report().material_energy.material.physical_heat_j,
            0.0,
            0.0,
        );
        near(step.report().thermal_audit.stored_energy_change_j, 0.0, 0.0);
        near(
            step.report().combined_residual_j,
            0.0,
            config().energy_tolerance_j,
        );
        assert!(
            step.solid()
                .velocity_mps()
                .unwrap()
                .iter()
                .flatten()
                .all(|v| *v == 0.0)
        );
        assert_eq!(step.observer().steps, 2);
        assert_eq!(accepted_observer.steps, 0);
        for i in 0..8 {
            near(step.thermal().temperature_k(cell(i)).unwrap(), 300.0, 1e-12);
        }
    });
}
