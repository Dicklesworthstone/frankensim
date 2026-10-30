//! G1 loaded-tet oracle, G3 spatial source/conduction feedback, G0 work and
//! history preservation, G4 late refusal. All material numbers are synthetic.
#[path = "solid_plastic_enthalpy/contact.rs"]
mod contact;
#[path = "solid_plastic_enthalpy/expansion_curve.rs"]
mod expansion_curve;
#[path = "solid_plastic_enthalpy/modulus_curve.rs"]
mod modulus_curve;
#[path = "solid_plastic_enthalpy/kinematic.rs"]
mod kinematic;
#[path = "solid_plastic_enthalpy/dynamics.rs"]
mod dynamics;
#[path = "solid_plastic_enthalpy/overstress.rs"]
mod overstress;

use std::collections::BTreeMap;

use fs_couple::solid_enthalpy::plastic::{
    ThermoplasticCell, ThermoplasticConfig, ThermoplasticSolid,
};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::plastic::source::HardeningPropertyBasis;
use fs_material::plastic::thermal::SourceQuery;
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, POISSON_RATIO_PROPERTY, YIELD_STRESS_PROPERTY,
    YOUNG_MODULUS_PROPERTY,
};
use fs_qty::{Density, Dims, Pressure, Temperature};
use fs_solid::linear3::TetAssemblyBudget;
use fs_solid::plastic3::IncrementSettings;
use fs_thermal::enthalpy::*;

const L: f64 = 0.01;
const V: f64 = L * L * L / 6.0;
const E: f64 = 1.2e9;
const H: f64 = 1e8;
const RHO: f64 = 1000.0;
const YIELD: f64 = 3e6;
const SOFTENING: f64 = 1e5;
const CP: &str = "specific_heat_capacity";
const CP_DIMS: Dims = Dims([2, 0, -2, -1, 0, 0]);
const BASIS: HardeningPropertyBasis = HardeningPropertyBasis::PlasticStrainModulus;
const NODES: [[f64; 3]; 4] = [[0.0; 3], [L, 0.0, 0.0], [0.0, L, 0.0], [0.0, 0.0, L]];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const ROLLERS: [(usize, f64); 9] = [
    (0, 0.0),
    (1, 0.0),
    (2, 0.0),
    (4, 0.0),
    (5, 0.0),
    (6, 0.0),
    (8, 0.0),
    (9, 0.0),
    (10, 0.0),
];

fn scalar(name: &str, dims: Dims, value: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained(),
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic spatial J2 fixture".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}
struct Fixture {
    card: MaterialCard,
    caloric: ResolvedHeatCapacityCurve,
    expansion: bool,
}
impl Fixture {
    fn new(cp: f64) -> Self {
        Self::with_expansion(cp, None)
    }
    fn with_expansion(cp: f64, alpha: Option<f64>) -> Self {
        Self::with_expansion_claim(
            cp,
            alpha.map(|alpha| {
                scalar(
                    fs_material::state_point::LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
                    fs_material::state_point::INVERSE_TEMPERATURE_DIMS,
                    alpha,
                )
            }),
        )
    }
    fn with_expansion_claim(cp: f64, expansion: Option<PropertyClaim>) -> Self {
        Self::with_constitutive_claims(
            cp,
            expansion,
            scalar(YOUNG_MODULUS_PROPERTY, Pressure::DIMS, E),
            scalar(BASIS.property_name(), Pressure::DIMS, H),
        )
    }
    fn with_constitutive_claims(
        cp: f64,
        expansion: Option<PropertyClaim>,
        youngs: PropertyClaim,
        hardening: PropertyClaim,
    ) -> Self {
        Self::with_density_and_constitutive_claims(cp, expansion, youngs, hardening, RHO)
    }
    fn with_density_and_constitutive_claims(
        cp: f64, expansion: Option<PropertyClaim>, youngs: PropertyClaim,
        hardening: PropertyClaim, density: f64,
    ) -> Self {
        let has_expansion = expansion.is_some();
        let mut claims = ClaimSet::new();
        for p in [
            youngs,
            scalar(POISSON_RATIO_PROPERTY, Dims([0; 6]), 0.2),
            hardening,
            scalar(DENSITY_PROPERTY, Density::DIMS, density),
            scalar(CP, CP_DIMS, cp),
        ] {
            claims.insert_claim(p).unwrap();
        }
        if let Some(expansion) = expansion {
            claims.insert_claim(expansion).unwrap();
        }
        let mut yield_claim = scalar(YIELD_STRESS_PROPERTY, Pressure::DIMS, YIELD);
        yield_claim.value = PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(280.0, 4e6), (300.0, YIELD), (320.0, 1e6)],
            dims: Pressure::DIMS,
        };
        yield_claim.interpolation = InterpolationPolicy::LinearInside;
        yield_claim.validity = ValidityDomain::unconstrained().with("T", 280.0, 320.0);
        claims.insert_claim(yield_claim).unwrap();
        let card = MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic J2".into(),
                phase: "solid".into(),
                process: "test".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap();
        let lower = QueryPoint::new().with("T", 280.0).unwrap();
        let upper = QueryPoint::new().with("T", 320.0).unwrap();
        let caloric = ResolvedHeatCapacityCurve::try_new(
            &card,
            HeatCapacitySourceQuery {
                heat_capacity: &PropertyKey::new(CP, CP_DIMS),
                density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
                temperature_axis: "T",
                lower: &lower,
                upper: &upper,
                selection: SelectionPolicy::SingleClaimOnly,
            },
            SolidLiquidPhase::Solid,
            -20.0 * cp,
            1e-5,
        )
        .unwrap();
        Self {
            card,
            caloric,
            expansion: has_expansion,
        }
    }
    fn query(&self) -> SourceQuery<'_> {
        SourceQuery {
            card: &self.card,
            caloric: &self.caloric,
            context: QueryPoint::new().with("T", 300.0).unwrap(),
            temperature_axis: "T".into(),
            selection: MaterialPropertySelection::SingleClaimOnly,
            hardening_basis: BASIS,
            thermal_expansion: self
                .expansion
                .then_some(MaterialPropertySelection::SingleClaimOnly),
            strain_limit: 0.04,
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
                seed: 1,
                kernel_id: 0x4a335448,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}
fn thermal_config() -> StepConfig {
    StepConfig {
        max_dt_s: 0.00025,
        energy_tolerance_j: 1e-10,
    }
}
fn config() -> ThermoplasticConfig {
    ThermoplasticConfig {
        maximum_iterations: 100,
        temperature_tolerance_k: 1e-8,
        energy_tolerance_j: 1e-8,
        yield_tolerance_pa: 1e-3,
        enthalpy_relaxation: 0.7,
    }
}
fn initial<'a>(
    fixtures: &[&'a Fixture],
    links: &[(usize, usize, f64)],
    powers: &[f64],
) -> DepositionSimulation<'a> {
    let mut events = Vec::new();
    for (i, fixture) in fixtures.iter().enumerate() {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(100 + i as u64),
                material_id: MaterialId(i as u64),
                segment_id: i as u64,
                body_id: 7,
                mass_kg: fixture.caloric.curve().state_at_specific_enthalpy(0.0)
                    .unwrap().bulk_density_kg_m3() * V,
                specific_h_j_kg: 0.0,
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [L; 3],
                    volume_m3: V,
                    surface_area_m2: 3.0 * L * L,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: powers[i],
                },
            }),
        });
    }
    for (i, &(a, b, g)) in links.iter().enumerate() {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: i as u64,
                a: CellId(100 + a as u64),
                b: CellId(100 + b as u64),
                area_m2: 0.5 * L * L,
                conductance_w_k: g,
            }),
        });
    }
    let laws: Vec<_> = fixtures
        .iter()
        .enumerate()
        .map(|(i, f)| (MaterialId(i as u64), f.caloric.curve() as &dyn ThermalLaw))
        .collect();
    DepositionSimulation::new(
        events,
        &laws,
        DepositionBudget {
            max_cells: 16,
            max_links: 16,
            max_events: 32,
            max_substeps: 10000,
        },
    )
    .unwrap()
    .advance_to(0.0, thermal_config(), || Ok(()))
    .unwrap()
    .0
}
fn bind<'a>(
    thermal: &DepositionSimulation<'_>,
    fixtures: &[&'a Fixture],
    nodes: &'a [[f64; 3]],
    tets: &'a [[usize; 4]],
    cx: &Cx<'_>,
) -> ThermoplasticSolid<'a> {
    ThermoplasticSolid::try_new(
        thermal,
        nodes,
        tets,
        fixtures
            .iter()
            .enumerate()
            .map(|(i, f)| ThermoplasticCell {
                cell: CellId(100 + i as u64),
                query: f.query(),
            })
            .collect(),
        TetAssemblyBudget::standard(),
        cx,
    )
    .unwrap()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

#[test]
fn g1_material_expansion_changes_free_shape_and_constrained_heating_on_native_mesh() {
    context(|cx| {
        let cp = 500.0;
        let duration = 0.001;
        let heat = RHO * V * cp * 5.0;
        let fixed: Vec<_> = (0..12).map(|i| (i, 0.0)).collect();
        for alpha in [1e-5, 2e-5] {
            let fixture = Fixture::with_expansion(cp, Some(alpha));
            let thermal = initial(&[&fixture], &[], &[heat / duration]);
            let solid = bind(&thermal, &[&fixture], &NODES, &TETS, cx);
            let free = solid
                .advance(
                    &thermal,
                    &(),
                    duration,
                    &ROLLERS,
                    &[],
                    IncrementSettings::default(),
                    thermal_config(),
                    config(),
                    cx,
                )
                .unwrap();
            close(
                free.thermal().temperature_k(CellId(100)).unwrap(),
                305.0,
                1e-8,
            );
            for (node, axis) in [(1, 0), (2, 1), (3, 2)] {
                close(
                    free.solid().solution().displacement_m[node][axis],
                    L * alpha * 5.0,
                    1e-12,
                );
            }
            for stress in free.solid().solution().element_stress[0] {
                close(stress, 0.0, 1e-3);
            }
            close(free.report().reversible_heat_j, 0.0, 1e-9);
            close(free.report().combined_residual_j, 0.0, 1e-8);
            assert_eq!(free.thermal().cells()[&CellId(100)].mass_kg, RHO * V);

            let clamped = solid
                .advance(
                    &thermal,
                    &(),
                    duration,
                    &fixed,
                    &[],
                    IncrementSettings::default(),
                    thermal_config(),
                    config(),
                    cx,
                )
                .unwrap();
            // Fully constrained heating: (rho Cp-9K alpha^2 T0) theta=Q/V.
            let nine_bulk = 3.0 * E / (1.0 - 0.4);
            let theta = heat / V / (RHO * cp - nine_bulk * alpha * alpha * 300.0);
            close(
                clamped.thermal().temperature_k(CellId(100)).unwrap(),
                300.0 + theta,
                2e-8,
            );
            for stress in &clamped.solid().solution().element_stress[0][..3] {
                close(*stress, -nine_bulk / 3.0 * alpha * theta, 2e-3);
            }
            assert!(clamped.report().reversible_heat_j > 0.0);
            assert!(clamped.report().thermoelastic_storage_change_j < 0.0);
            assert_eq!(clamped.report().physical_heat_j, 0.0);
            close(clamped.report().combined_residual_j, 0.0, 1e-8);
            let continued = clamped
                .solid()
                .advance(
                    clamped.thermal(),
                    &(),
                    2.0 * duration,
                    &fixed,
                    &[],
                    IncrementSettings::default(),
                    thermal_config(),
                    config(),
                    cx,
                )
                .unwrap();
            let old_t = 300.0 + theta;
            let next_theta = heat / V / (RHO * cp - nine_bulk * alpha * alpha * old_t);
            close(
                continued.thermal().temperature_k(CellId(100)).unwrap(),
                old_t + next_theta,
                3e-8,
            );
            close(continued.report().combined_residual_j, 0.0, 1e-8);
            assert_eq!(solid.time_s(), 0.0);
        }
    });
}

#[test]
fn g3_thermal_expansion_drives_native_constrained_plasticity_without_an_authored_force() {
    context(|cx| {
        let cp = 1000.0;
        let alpha = 5e-4;
        let duration = 0.001;
        let fixture = Fixture::with_expansion(cp, Some(alpha));
        let thermal = initial(&[&fixture], &[], &[RHO * V * cp * 5.0 / duration]);
        let solid = bind(&thermal, &[&fixture], &NODES, &TETS, cx);
        let mut axial_clamp = ROLLERS.to_vec();
        axial_clamp.push((3, 0.0));
        let candidate = solid
            .advance(
                &thermal,
                &(),
                duration,
                &axial_clamp,
                &[],
                IncrementSettings::default(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        let solution = candidate.solid().solution();
        let state = &solution.element_states[0];
        let t = candidate.thermal().temperature_k(CellId(100)).unwrap();
        assert!(state.alpha > 0.0);
        assert!(state.plastic_strain[0] < 0.0);
        assert_eq!(solution.displacement_m[1][0], 0.0);
        close(solution.element_stress[0][1], 0.0, 1e-3);
        close(solution.element_stress[0][2], 0.0, 1e-3);
        let yield_pa = YIELD - SOFTENING * (t - 300.0);
        close(
            -solution.element_stress[0][0],
            yield_pa + H * state.alpha,
            1e-3,
        );
        close(
            candidate.report().physical_heat_j,
            V * yield_pa * state.alpha,
            1e-9,
        );
        assert!(candidate.report().reversible_heat_j > 0.0);
        close(candidate.report().mechanical_work_j, 0.0, 1e-10);
        close(candidate.report().combined_residual_j, 0.0, 1e-8);
    });
}

#[test]
fn loaded_tet_temperature_shape_and_heat_follow_independent_uniaxial_oracle() {
    context(|cx| {
        let sigma = 3.4e6;
        let mut responses = Vec::new();
        for cp in [10.0, 100.0] {
            let fixture = Fixture::new(cp);
            let thermal = initial(&[&fixture], &[], &[0.0]);
            let solid = bind(&thermal, &[&fixture], &NODES, &TETS, cx);
            let candidate = solid
                .advance(
                    &thermal,
                    &(),
                    0.001,
                    &ROLLERS,
                    &[(3, sigma * V / L)],
                    IncrementSettings::default(),
                    thermal_config(),
                    config(),
                    cx,
                )
                .unwrap();
            // Uniform uniaxial stress: alpha=(sigma-y0+a*dT)/H, and
            // rho Cp dT=(y0-a*dT)*alpha. This quadratic is independent of FEM/J2.
            let b = RHO * cp * H - SOFTENING * (2.0 * YIELD - sigma);
            let c = YIELD * (sigma - YIELD);
            let theta = 2.0 * c / (b + (b * b + 4.0 * SOFTENING * SOFTENING * c).sqrt());
            let alpha = (sigma - YIELD + SOFTENING * theta) / H;
            let eps = sigma / E + alpha;
            close(
                candidate.thermal().temperature_k(CellId(100)).unwrap(),
                300.0 + theta,
                2e-7,
            );
            close(
                candidate.solid().solution().displacement_m[1][0],
                L * eps,
                2e-11,
            );
            close(
                candidate.solid().solution().element_states[0].alpha,
                alpha,
                2e-9,
            );
            close(
                candidate.solid().solution().element_stress[0][0],
                sigma,
                2e-3,
            );
            let report = candidate.report();
            close(report.mechanical_work_j, sigma * V * eps, 1e-9);
            close(report.physical_heat_j, RHO * V * cp * theta, 1e-9);
            close(report.external_heat_j, 0.0, 1e-12);
            close(report.combined_residual_j, 0.0, 1e-8);
            assert!(report.numerical_loss_j > 0.0);
            assert!(report.element_energy[0].energy.final_state.hardening_j_m3 > 0.0);
            close(candidate.thermal().energy_j(), report.physical_heat_j, 1e-9);
            assert_eq!(
                candidate.thermal().cells()[&CellId(100)].boundary.source_w,
                0.0
            );
            assert!(
                candidate
                    .thermal()
                    .uses_material_law(MaterialId(0), fixture.caloric.curve())
            );
            responses.push((
                theta,
                candidate.solid().solution().displacement_m[1][0],
                candidate.solid().materials()[0].material().identity(),
            ));
        }
        assert!(responses[0].0 > responses[1].0);
        assert!(responses[0].1 > responses[1].1);
        assert_ne!(responses[0].2, responses[1].2);
    });
}

#[test]
fn force_unloading_retains_plastic_history_and_does_not_turn_numerical_loss_into_heat() {
    context(|cx| {
        let fixture = Fixture::new(100.0);
        let thermal = initial(&[&fixture], &[], &[0.0]);
        let solid = bind(&thermal, &[&fixture], &NODES, &TETS, cx);
        let loaded = solid
            .advance(
                &thermal,
                &(),
                0.001,
                &ROLLERS,
                &[(3, 3.4e6 * V / L)],
                IncrementSettings::default(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        let alpha = loaded.solid().solution().element_states[0].alpha;
        let loaded_history = loaded.solid().solution().element_states.clone();
        let unloaded = loaded
            .solid()
            .advance(
                loaded.thermal(),
                &(),
                0.002,
                &ROLLERS,
                &[],
                IncrementSettings::default(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        assert_eq!(unloaded.solid().solution().element_states, loaded_history);
        close(
            unloaded.solid().solution().displacement_m[1][0],
            L * alpha,
            1e-11,
        );
        close(unloaded.solid().solution().element_stress[0][0], 0.0, 1e-3);
        close(
            unloaded.thermal().energy_j(),
            loaded.thermal().energy_j(),
            1e-11,
        );
        close(unloaded.report().physical_heat_j, 0.0, 1e-15);
        assert!(unloaded.report().numerical_loss_j > 0.0);
        close(unloaded.report().combined_residual_j, 0.0, 1e-8);
        // The original accepted loaded candidate remains available unchanged.
        assert_eq!(loaded.solid().solution().element_states, loaded_history);
        assert_eq!(loaded.solid().time_s(), 0.001);
    });
}

fn cube() -> (Vec<[f64; 3]>, Vec<[usize; 4]>) {
    let nodes: Vec<_> = (0..8)
        .map(|i| {
            [
                L * (i & 1) as f64,
                L * ((i >> 1) & 1) as f64,
                L * ((i >> 2) & 1) as f64,
            ]
        })
        .collect();
    // Positive winding; neighboring entries share complete triangular faces.
    (
        nodes,
        vec![
            [0, 1, 3, 7],
            [0, 3, 2, 7],
            [0, 2, 6, 7],
            [0, 6, 4, 7],
            [0, 4, 5, 7],
            [0, 5, 1, 7],
        ],
    )
}
fn grips(nodes: &[[f64; 3]], strain: f64) -> Vec<(usize, f64)> {
    let mut dofs = Vec::new();
    for (i, node) in nodes.iter().enumerate() {
        dofs.push((3 * i, strain * node[0]));
        if node[1] == 0.0 {
            dofs.push((3 * i + 1, 0.0));
        }
        if node[2] == 0.0 {
            dofs.push((3 * i + 2, 0.0));
        }
    }
    dofs
}
fn tension(solid: &ThermoplasticSolid<'_>, prescribed: &[(usize, f64)]) -> f64 {
    prescribed
        .iter()
        .zip(&solid.solution().increments[0].reactions_n)
        .filter(|((dof, value), _)| dof % 3 == 0 && *value > 0.0)
        .map(|(_, reaction)| reaction)
        .sum()
}

#[test]
fn spatial_heating_conducts_between_actual_loaded_elements_and_relaxes_grip_force() {
    context(|cx| {
        let fixture = Fixture::new(20.0);
        let (nodes, tets) = cube();
        let prescribed = grips(&nodes, 0.006);
        let mut spreads = Vec::new();
        for conductance in [0.0, 0.03] {
            let links: Vec<_> = (0..6).map(|i| (i, (i + 1) % 6, conductance)).collect();
            let fixtures = vec![&fixture; 6];
            let thermal = initial(&fixtures, &links, &[0.04, 0.0, 0.0, 0.0, 0.0, 0.0]);
            let solid = bind(&thermal, &fixtures, &nodes, &tets, cx);
            let first = solid
                .advance(
                    &thermal,
                    &(),
                    0.001,
                    &prescribed,
                    &[],
                    IncrementSettings::default(),
                    thermal_config(),
                    config(),
                    cx,
                )
                .unwrap();
            let first_force = tension(first.solid(), &prescribed);
            let (mut thermal, mut solid, (), first_report) = first.into_parts();
            let mut work = first_report.mechanical_work_j;
            let mut external = first_report.external_heat_j;
            let mut numerical = first_report.numerical_loss_j;
            let mut storage = first_report.stored_mechanical_change_j;
            for step in 2..=25 {
                let (next_thermal, next_solid, (), report) = solid
                    .advance(
                        &thermal,
                        &(),
                        step as f64 * 0.001,
                        &prescribed,
                        &[],
                        IncrementSettings::default(),
                        thermal_config(),
                        config(),
                        cx,
                    )
                    .unwrap()
                    .into_parts();
                close(report.mechanical_work_j, 0.0, 1e-12);
                work += report.mechanical_work_j;
                external += report.external_heat_j;
                numerical += report.numerical_loss_j;
                storage += report.stored_mechanical_change_j;
                thermal = next_thermal;
                solid = next_solid;
            }
            assert!(tension(&solid, &prescribed) < first_force);
            close(
                thermal.energy_j() + storage + numerical - work - external,
                0.0,
                2e-7,
            );
            let temperatures: Vec<_> = (0..6)
                .map(|i| thermal.temperature_k(CellId(100 + i)).unwrap())
                .collect();
            let maximum = temperatures
                .iter()
                .copied()
                .fold(f64::NEG_INFINITY, f64::max);
            let minimum = temperatures.iter().copied().fold(f64::INFINITY, f64::min);
            spreads.push(maximum - minimum);
            for (i, material) in solid.materials().iter().enumerate() {
                close(
                    material.law().yield_stress,
                    YIELD - SOFTENING * (temperatures[i] - 300.0),
                    1e-6,
                );
                assert!(solid.solution().element_states[i].alpha > 0.0);
            }
        }
        assert!(spreads[1] < spreads[0]);
    });
}

#[derive(Clone)]
struct LateCancellation {
    steps: usize,
}
impl DepositionObserver for LateCancellation {
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
        self.steps += 1;
        if self.steps == 2 {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    }
}
#[test]
fn late_cancellation_source_exit_and_stale_clock_preserve_all_accepted_owners() {
    context(|cx| {
        let fixture = Fixture::new(10.0);
        let thermal = initial(&[&fixture], &[], &[0.0]);
        let solid = bind(&thermal, &[&fixture], &NODES, &TETS, cx);
        let original = solid.solution().clone();
        let cells = thermal.cells().clone();
        let observer = LateCancellation { steps: 0 };
        assert!(
            solid
                .advance(
                    &thermal,
                    &observer,
                    0.001,
                    &ROLLERS,
                    &[(3, 3.4e6 * V / L)],
                    IncrementSettings::default(),
                    thermal_config(),
                    config(),
                    cx
                )
                .is_err()
        );
        assert_eq!(observer.steps, 0);
        // The local caloric support ends at 320 K; externally heating through it
        // cannot be hidden by the converged mechanical endpoint or later cooling.
        let hot = initial(&[&fixture], &[], &[100.0]);
        assert!(
            solid
                .advance(
                    &hot,
                    &(),
                    0.001,
                    &ROLLERS,
                    &[],
                    IncrementSettings::default(),
                    thermal_config(),
                    config(),
                    cx
                )
                .is_err()
        );
        let shifted = thermal
            .advance_to(0.0005, thermal_config(), || Ok(()))
            .unwrap()
            .0;
        assert!(
            solid
                .advance(
                    &shifted,
                    &(),
                    0.001,
                    &ROLLERS,
                    &[],
                    IncrementSettings::default(),
                    thermal_config(),
                    config(),
                    cx
                )
                .is_err()
        );
        assert_eq!(solid.solution(), &original);
        assert_eq!(thermal.cells(), &cells);
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(solid.time_s(), 0.0);
    });
}
