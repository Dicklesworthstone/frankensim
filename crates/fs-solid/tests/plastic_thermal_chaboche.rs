//! Actual tetrahedral equilibrium with transported, source-defined cyclic history.
//! Synthetic G0/G1 data; prescribed temperature, small strain, no heat-transport
//! or experimental material-calibration claim.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::SmallStrainLaw;
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::plastic::chaboche::ChabocheState;
use fs_material::plastic::chaboche::source::BackstressProperties;
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::rate::{
    OVERSTRESS_EXPONENT_PROPERTY, OVERSTRESS_STRESS_SCALE_PROPERTY, OVERSTRESS_VISCOSITY_DIMS,
    OVERSTRESS_VISCOSITY_PROPERTY,
};
use fs_material::plastic::source::{HardeningPropertyBasis, KINEMATIC_HARDENING_MODULUS_PROPERTY};
use fs_material::plastic::thermal::chaboche::{ThermalChaboche, ThermalChabocheSource};
use fs_material::plastic::thermal::{SourceQuery, ThermoplasticSource};
use fs_material::plastic::voce::source::{
    VOCE_SATURATION_INCREMENT_PROPERTY as Q, VOCE_SATURATION_RATE_PROPERTY as B,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, POISSON_RATIO_PROPERTY, YIELD_STRESS_PROPERTY,
    YOUNG_MODULUS_PROPERTY,
};
use fs_material::tensor::contract;
use fs_qty::{Density, Dims, Pressure, Temperature};
use fs_solid::linear3::TetAssemblyBudget;
use fs_solid::plastic3::{
    IncrementSettings, SmallStrainLawField, SmallStrainTetProblem, SmallStrainTetSolution,
};

const CP: f64 = 500.0;
const CP_DIMS: Dims = Dims([2, 0, -2, -1, 0, 0]);
const BASIS: HardeningPropertyBasis = HardeningPropertyBasis::MixedPlasticStrainModuli;
const SIZE: [f64; 3] = [0.03, 0.01, 0.01];
const TOL: EnergyTolerances = EnergyTolerances {
    energy_j_m3: 2e-5,
    yield_pa: 1e-3,
};
type Solution = SmallStrainTetSolution<ChabocheState<2>>;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 13,
                kernel_id: 31,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}
fn schema() -> [BackstressProperties; 2] {
    [
        BackstressProperties::new("fast_C", "fast_gamma"),
        BackstressProperties::new("slow_C", "slow_gamma"),
    ]
}
fn scalar(name: &str, dims: Dims, value: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained(),
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic thermal cyclic continuum fixture".into(),
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
    fn new(variable_branches: bool) -> Self {
        Self::new_with_rate(variable_branches, None)
    }

    fn new_with_rate(variable_branches: bool, exponent: Option<f64>) -> Self {
        let mut properties = vec![
            scalar(YOUNG_MODULUS_PROPERTY, Pressure::DIMS, 1.2e9),
            scalar(POISSON_RATIO_PROPERTY, Dims([0; 6]), 0.2),
            scalar(YIELD_STRESS_PROPERTY, Pressure::DIMS, 3e6),
            scalar(BASIS.property_name(), Pressure::DIMS, 1e8),
            scalar(KINEMATIC_HARDENING_MODULUS_PROPERTY, Pressure::DIMS, 1.6e8),
            scalar("specific_heat_capacity", CP_DIMS, CP),
            scalar(DENSITY_PROPERTY, Density::DIMS, 1000.0),
            scalar(Q, Pressure::DIMS, 6e6),
            scalar(B, Dims([0; 6]), 300.0),
        ];
        if let Some(exponent) = exponent {
            let saturation = properties.iter_mut().find(|p| p.key.name() == Q).unwrap();
            saturation.validity = ValidityDomain::unconstrained().with("T", 280.0, 320.0);
            saturation.interpolation = InterpolationPolicy::LinearInside;
            saturation.value = PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: [280.0, 300.0, 320.0]
                    .map(|t| (t, 6e6 - 1e5 * (t - 300.0)))
                    .to_vec(),
                dims: Pressure::DIMS,
            };
            properties.extend([
                scalar(
                    OVERSTRESS_VISCOSITY_PROPERTY,
                    OVERSTRESS_VISCOSITY_DIMS,
                    1e6,
                ),
                scalar(OVERSTRESS_STRESS_SCALE_PROPERTY, Pressure::DIMS, 1e6),
                scalar(OVERSTRESS_EXPONENT_PROPERTY, Dims([0; 6]), exponent),
            ]);
        }
        for (i, names) in schema().iter().enumerate() {
            let c = [4e8, 1e8][i];
            let mut modulus = scalar(&names.modulus, Pressure::DIMS, c);
            if variable_branches {
                modulus.validity = ValidityDomain::unconstrained().with("T", 280.0, 320.0);
                modulus.interpolation = InterpolationPolicy::LinearInside;
                modulus.value = PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: [280.0, 300.0, 320.0]
                        .map(|t| (t, c + [-8e6, 2e6][i] * (t - 300.0)))
                        .to_vec(),
                    dims: Pressure::DIMS,
                };
            }
            properties.push(modulus);
            properties.push(scalar(&names.recovery, Dims([0; 6]), [120.0, 15.0][i]));
        }
        let mut claims = ClaimSet::new();
        for property in properties {
            claims.insert_claim(property).unwrap();
        }
        let card = MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic cyclic solid".into(),
                phase: "solid".into(),
                process: "numerical fixture".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap();
        let caloric = ResolvedHeatCapacityCurve::try_new(
            &card,
            HeatCapacitySourceQuery {
                heat_capacity: &PropertyKey::new("specific_heat_capacity", CP_DIMS),
                density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
                temperature_axis: "T",
                lower: &QueryPoint::new().with("T", 280.0).unwrap(),
                upper: &QueryPoint::new().with("T", 320.0).unwrap(),
                selection: SelectionPolicy::SingleClaimOnly,
            },
            SolidLiquidPhase::Solid,
            -20.0 * CP,
            1e-6,
        )
        .unwrap();
        Self { card, caloric }
    }
    fn source(&self) -> ThermalChabocheSource<'_, 2> {
        self.source_with_rate(false)
    }

    fn source_with_rate(&self, kinetic: bool) -> ThermalChabocheSource<'_, 2> {
        let base = ThermoplasticSource::try_new(
            SourceQuery {
                card: &self.card,
                caloric: &self.caloric,
                context: QueryPoint::new().with("T", 300.0).unwrap(),
                temperature_axis: "T".into(),
                selection: MaterialPropertySelection::SingleClaimOnly,
                hardening_basis: BASIS,
                thermal_expansion: None,
                strain_limit: 0.04,
            },
            0.0,
        )
        .unwrap()
        .with_voce_hardening(MaterialPropertySelection::SingleClaimOnly)
        .unwrap();
        let base = if kinetic {
            base.with_power_law_overstress(MaterialPropertySelection::SingleClaimOnly)
                .unwrap()
        } else {
            base
        };
        base.with_chaboche_branches(MaterialPropertySelection::SingleClaimOnly, schema())
            .unwrap()
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, allowance {tolerance:e}"
    );
}
fn member_mesh() -> (Vec<[f64; 3]>, [[usize; 4]; 6]) {
    let nodes = (0..8)
        .map(|i| {
            [
                SIZE[0] * (i & 1) as f64,
                SIZE[1] * ((i >> 1) & 1) as f64,
                SIZE[2] * ((i >> 2) & 1) as f64,
            ]
        })
        .collect();
    (
        nodes,
        [
            [0, 1, 3, 7],
            [0, 1, 7, 5],
            [0, 2, 7, 3],
            [0, 2, 6, 7],
            [0, 4, 5, 7],
            [0, 4, 7, 6],
        ],
    )
}
fn member_step(
    source: &ThermalChabocheSource<'_, 2>,
    old_h: f64,
    new_h: f64,
    axial_strain: f64,
    old: Option<&Solution>,
    cx: &Cx<'_>,
) -> (Solution, f64) {
    let (_, law) = source.resolve_enthalpy(old_h, new_h).unwrap();
    member_step_with_law(&law, axial_strain, old, cx)
}

fn member_step_with_law(
    law: &ThermalChaboche<2>,
    axial_strain: f64,
    old: Option<&Solution>,
    cx: &Cx<'_>,
) -> (Solution, f64) {
    assert!(!law.admissibility().tangent_symmetric);
    let (nodes, tets) = member_mesh();
    let mut supports = Vec::new();
    for (node, x) in nodes.iter().enumerate() {
        supports.push((3 * node, x[0]));
        if x[1] == 0.0 {
            supports.push((3 * node + 1, 0.0));
        }
        if x[2] == 0.0 {
            supports.push((3 * node + 2, 0.0));
        }
    }
    // One resolved thermal transition per solve; all Newton proposals borrow
    // exactly the old state. Reusing a transition for subsequent increments
    // would apply temperature transport twice.
    let problem = SmallStrainTetProblem {
        nodes_m: &nodes,
        tetrahedra: &tets,
        laws: SmallStrainLawField::Uniform(law),
        prescribed_m: &supports,
        nodal_forces_n: &[],
        load_path: &[axial_strain],
        budget: TetAssemblyBudget::standard(),
        settings: IncrementSettings::default(),
    };
    let previous = old
        .map(|s| s.element_states.clone())
        .unwrap_or_else(|| vec![law.initial_state(); tets.len()]);
    let solution = match old {
        Some(accepted) => problem
            .solve_from(
                &accepted
                    .displacement_m
                    .iter()
                    .flatten()
                    .copied()
                    .collect::<Vec<_>>(),
                &accepted.element_states,
                cx,
            )
            .unwrap(),
        None => problem.solve(cx).unwrap(),
    };
    let flat: Vec<f64> = solution.displacement_m.iter().flatten().copied().collect();
    let force = problem.internal_force(&flat, &previous, cx).unwrap();
    for (dof, value) in force.iter().enumerate() {
        if !supports.iter().any(|&(fixed, _)| fixed == dof) {
            near(*value, 0.0, 1e-7);
        }
    }
    let record = &solution.increments[0];
    let right: f64 = supports
        .iter()
        .zip(&record.reactions_n)
        .filter(|((dof, value), _)| dof % 3 == 0 && *value > 0.0)
        .map(|(_, f)| f)
        .sum();
    let left: f64 = supports
        .iter()
        .zip(&record.reactions_n)
        .filter(|((dof, value), _)| dof % 3 == 0 && *value == 0.0)
        .map(|(_, f)| f)
        .sum();
    near(left + right, 0.0, 1e-7);
    for stress in &solution.element_stress {
        near(stress[0] * SIZE[1] * SIZE[2], right, 1e-6);
        for component in &stress[1..] {
            near(*component, 0.0, 0.02);
        }
    }
    (solution, right)
}

#[test]
fn g1_heated_held_member_changes_tension_shape_and_cyclic_history() {
    let fixture = Fixture::new(true);
    let fixed_fixture = Fixture::new(false);
    let source = fixture.source();
    let fixed = fixed_fixture.source();
    with_cx(|cx| {
        let (cold, cold_force) = member_step(&source, 0.0, 0.0, 0.012, None, cx);
        let old_states = cold.element_states.clone();
        assert!(
            old_states
                .iter()
                .all(|s| s.plastic.alpha > 0.001 && s.backstress_pa.iter().all(|x| x[0] > 0.0))
        );
        let hot_h = CP * 15.0;
        let (hot, hot_force) = member_step(&source, 0.0, hot_h, 0.012, Some(&cold), cx);
        let (unchanged, fixed_force) = member_step(&fixed, 0.0, hot_h, 0.012, Some(&cold), cx);
        assert_eq!(cold.element_states, old_states);
        near(fixed_force, cold_force, 1e-6);
        assert!(
            hot_force < fixed_force - 0.5,
            "thermal branch evolution must alter actual reaction"
        );
        assert!(hot.displacement_m[7][1] < unchanged.displacement_m[7][1] - 1e-8);
        let mut work_j = 0.0;
        let mut coefficient_j = 0.0;
        let volume = SIZE.iter().product::<f64>() / 6.0;
        for e in 0..6 {
            assert!(hot.element_states[e].plastic.alpha > cold.element_states[e].plastic.alpha);
            let energy = source
                .advance_with_energy(
                    0.0,
                    hot_h,
                    &cold.element_strain[e],
                    &cold.element_states[e],
                    &hot.element_strain[e],
                    TOL,
                )
                .unwrap();
            assert_eq!(energy.state, hot.element_states[e]);
            let account = &energy.thermal;
            assert!(energy.recovery_dissipation_j_m3 > 0.0);
            assert!(account.mechanical.energy.physical_dissipation_j_m3 > 0.0);
            near(
                account.mechanical_work_j_m3
                    - account.stored_change_j_m3()
                    - account.thermoelastic_storage_change_j_m3
                    - account.caloric_heat_j_m3()
                    - account.mechanical.energy.numerical_dissipation_j_m3,
                0.0,
                2e-5,
            );
            work_j += volume * account.mechanical_work_j_m3;
            coefficient_j += volume * account.coefficient_storage_change_j_m3;
        }
        near(work_j, 0.0, 2e-10); // Held supports and traction-free lateral faces do no work.
        assert!(coefficient_j < 0.0);
        let (cooled, cooled_force) = member_step(&source, hot_h, 0.0, 0.012, Some(&hot), cx);
        near(cooled_force, hot_force, 1e-6);
        for e in 0..6 {
            near(
                cooled.element_states[e].plastic.alpha,
                hot.element_states[e].plastic.alpha,
                1e-14,
            );
            for branch in 0..2 {
                let ratio = [4e8 / 2.8e8, 1e8 / 1.3e8][branch];
                for component in 0..6 {
                    near(
                        cooled.element_states[e].backstress_pa[branch][component],
                        hot.element_states[e].backstress_pa[branch][component] * ratio,
                        1e-7,
                    );
                }
            }
        }
        let (reversed, reversal_force) = member_step(&source, 0.0, 0.0, -0.009, Some(&cooled), cx);
        assert!(reversal_force < 0.0);
        assert!(
            reversed
                .element_states
                .iter()
                .zip(&cooled.element_states)
                .all(|(new, old)| new.plastic.alpha > old.plastic.alpha)
        );
        println!(
            "held member: cold={cold_force:.9} N hot={hot_force:.9} N, lateral cold={:.9e} m hot={:.9e} m, reversal={reversal_force:.9} N",
            cold.displacement_m[7][1], hot.displacement_m[7][1]
        );
    });
}

#[test]
fn g1_timed_heating_reversal_and_hold_retain_branches_and_separate_physical_losses() {
    let fixture = Fixture::new_with_rate(true, Some(3.0));
    let source = fixture.source_with_rate(true);
    assert!(source.resolve_enthalpy(0.0, 0.0).is_err());
    let volume = SIZE.iter().product::<f64>() / 6.0;
    let hot = 15.0 * CP;
    let stages = [
        (0.0, 0.012, 0.001),
        (hot, 0.012, 0.004),
        (hot, 0.012, 0.008),
        (hot, -0.009, 0.003),
        (hot, -0.009, 0.006),
    ];
    with_cx(|cx| {
        let mut accepted: Option<Solution> = None;
        let mut old_h = 0.0;
        let mut old_axial = 0.0;
        let mut old_force: f64 = 0.0;
        for (step, (new_h, axial, duration)) in stages.into_iter().enumerate() {
            let (_, law) = source
                .resolve_enthalpy_for_step(old_h, new_h, duration)
                .unwrap();
            let original = accepted.clone();
            let (next, force) = member_step_with_law(&law, axial, accepted.as_ref(), cx);
            let native = law.law();
            let previous = accepted
                .as_ref()
                .map(|s| s.element_states.clone())
                .unwrap_or_else(|| vec![law.initial_state(); 6]);
            let mut physical_j = 0.0;
            let mut viscous_j = 0.0;
            let mut recovery_j = 0.0;
            let mut coefficient_j = 0.0;
            let mut work_j = 0.0;
            for e in 0..6 {
                let old_strain = accepted.as_ref().map_or([0.0; 6], |s| s.element_strain[e]);
                let energy = source
                    .advance_with_energy_for_step(
                        old_h,
                        new_h,
                        &old_strain,
                        &previous[e],
                        &next.element_strain[e],
                        TOL,
                        duration,
                    )
                    .unwrap();
                assert_eq!(energy.state, next.element_states[e]);
                let da = energy.state.plastic.alpha - previous[e].plastic.alpha;
                assert!(da > 0.0);
                // Independently evaluate the selected Perzyna flow relation.
                let overstress = 1e6 * (da / duration).powf(1.0 / 3.0);
                near(
                    native.yield_function_with_state(&next.element_stress[e], &energy.state),
                    overstress,
                    TOL.yield_pa,
                );
                let recovery: f64 = native
                    .branches()
                    .iter()
                    .zip(&energy.state.backstress_pa)
                    .map(|(branch, x)| {
                        1.5 * branch.recovery_per_strain() * da * contract(x, x)
                            / branch.modulus_pa()
                    })
                    .sum();
                near(energy.recovery_dissipation_j_m3, recovery, TOL.energy_j_m3);
                let account = &energy.thermal;
                near(
                    account.mechanical.energy.physical_dissipation_j_m3,
                    (3e6 + overstress) * da + recovery,
                    TOL.energy_j_m3,
                );
                // Q and both C branches contribute to frozen-history energy;
                // neither coefficient work nor numerical loss is plastic heat.
                let old_temperature = 300.0 + old_h / CP;
                let new_temperature = 300.0 + new_h / CP;
                let old_alpha = previous[e].plastic.alpha;
                let voce_shape = old_alpha + (-300.0 * old_alpha).exp_m1() / 300.0;
                let mut coefficient = -1e5 * (new_temperature - old_temperature) * voce_shape;
                for branch in 0..2 {
                    let slope = [-8e6, 2e6][branch];
                    let c0 = [4e8, 1e8][branch] + slope * (old_temperature - 300.0);
                    let c1 = [4e8, 1e8][branch] + slope * (new_temperature - 300.0);
                    let shape =
                        0.75 * contract(
                            &previous[e].backstress_pa[branch],
                            &previous[e].backstress_pa[branch],
                        ) / (c0 * c0);
                    coefficient += (c1 - c0) * shape;
                }
                near(
                    account.coefficient_storage_change_j_m3,
                    coefficient,
                    TOL.energy_j_m3,
                );
                near(
                    account.mechanical_work_j_m3
                        - account.stored_change_j_m3()
                        - account.thermoelastic_storage_change_j_m3
                        - account.caloric_heat_j_m3()
                        - account.mechanical.energy.numerical_dissipation_j_m3,
                    0.0,
                    TOL.energy_j_m3,
                );
                assert!(account.mechanical.energy.numerical_dissipation_j_m3 >= 0.0);
                physical_j += volume * account.mechanical.energy.physical_dissipation_j_m3;
                viscous_j += volume * overstress * da;
                recovery_j += volume * recovery;
                coefficient_j += volume * coefficient;
                work_j += volume * account.mechanical_work_j_m3;
            }
            assert!(physical_j > 0.0 && viscous_j > 0.0 && recovery_j > 0.0);
            near(work_j, force * SIZE[0] * (axial - old_axial), 1e-9);
            if axial == old_axial {
                assert!(force.abs() < old_force.abs());
                near(work_j, 0.0, 1e-9);
            }
            if step == 1 {
                assert!(coefficient_j < 0.0);
                // Identical frozen-history Newton proposals are repeatable;
                // retrying does not spend another interval of viscous time.
                let (retry, retry_force) = member_step_with_law(&law, axial, accepted.as_ref(), cx);
                assert_eq!(retry, next);
                assert_eq!(retry_force, force);
                let (_, shorter) = source
                    .resolve_enthalpy_for_step(old_h, new_h, duration / 10.0)
                    .unwrap();
                let (short, short_force) =
                    member_step_with_law(&shorter, axial, accepted.as_ref(), cx);
                assert!(force < short_force - 0.1);
                assert!(
                    next.element_states[0].plastic.alpha > short.element_states[0].plastic.alpha
                );
            }
            if step == 3 {
                assert!(old_force > 0.0 && force < 0.0);
                assert!(
                    next.element_states[0]
                        .backstress_pa
                        .iter()
                        .all(|x| x[0] < 0.0)
                );
            }
            assert_eq!(accepted, original);
            println!(
                "timed cyclic member step={step} dt={duration:.6} s force={force:.9} N alpha={:.9e}, viscous={viscous_j:.9e} J recovery={recovery_j:.9e} J",
                next.element_states[0].plastic.alpha
            );
            accepted = Some(next);
            old_h = new_h;
            old_axial = axial;
            old_force = force;
        }
    });
}

fn tet_displacement(e: [f64; 6]) -> [f64; 12] {
    [
        0.0,
        0.0,
        0.0,
        e[0],
        0.0,
        0.0,
        2.0 * e[3],
        e[1],
        0.0,
        2.0 * e[5],
        2.0 * e[4],
        e[2],
    ]
}

#[test]
fn g0_heated_nonproportional_traction_preserves_the_full_continuum_jacobian() {
    let fixture = Fixture::new(true);
    let source = fixture.source();
    let (_, cold) = source.resolve_enthalpy(0.0, 0.0).unwrap();
    let old_strain = [0.011, -0.004, -0.003, 0.0015, 0.0, 0.0];
    let old = cold.update_state(&old_strain, &cold.initial_state());
    let (_, hot) = source.resolve_enthalpy(0.0, CP * 15.0).unwrap();
    let target = [0.004, 0.009, -0.003, -0.004, 0.001, 0.002];
    let expected = hot.response(&target, &old).unwrap();
    let s = expected.stress_pa;
    // Unit-coordinate simplex, V=1/6. These are independently integrated
    // uniform boundary tractions, retaining all six stress components.
    let loads = [
        (3, s[0] / 6.0),
        (6, s[3] / 6.0),
        (7, s[1] / 6.0),
        (9, s[5] / 6.0),
        (10, s[4] / 6.0),
        (11, s[2] / 6.0),
    ];
    let supports = [(0, 0.0), (1, 0.0), (2, 0.0), (4, 0.0), (5, 0.0), (8, 0.0)];
    let nodes = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let problem = SmallStrainTetProblem {
        nodes_m: &nodes,
        tetrahedra: &[[0, 1, 2, 3]],
        laws: SmallStrainLawField::Uniform(&hot),
        prescribed_m: &supports,
        nodal_forces_n: &loads,
        load_path: &[1.0],
        budget: TetAssemblyBudget::standard(),
        settings: IncrementSettings::default(),
    };
    with_cx(|cx| {
        let u = tet_displacement(target);
        let incoming = [old.clone()];
        let matrix = problem.tangent_matrix(&u, &incoming, cx).unwrap();
        let mut asymmetry = 0.0_f64;
        for &(i, _) in &loads {
            for &(j, _) in &loads {
                asymmetry = asymmetry.max((matrix.get(i, j) - matrix.get(j, i)).abs());
            }
        }
        assert!(
            asymmetry > 1e3,
            "fixture must exercise a nonsymmetric free tangent"
        );
        for &(column, _) in &loads {
            let (mut plus, mut minus) = (u, u);
            plus[column] += 1e-7;
            minus[column] -= 1e-7;
            let above = problem.internal_force(&plus, &incoming, cx).unwrap();
            let below = problem.internal_force(&minus, &incoming, cx).unwrap();
            for row in 0..12 {
                near(
                    matrix.get(row, column),
                    (above[row] - below[row]) / 2e-7,
                    20.0,
                );
            }
        }
        let result = problem
            .solve_from(&tet_displacement(old_strain), &incoming, cx)
            .unwrap();
        for (actual, expected) in result.displacement_m.iter().flatten().zip(u) {
            near(*actual, expected, 2e-10);
        }
        for (actual, expected) in result.element_stress[0].iter().zip(s) {
            near(*actual, expected, 0.02);
        }
        assert_eq!(incoming[0], old);
        assert!(result.increments[0].residual_history.len() > 1);
        let flat: Vec<f64> = result.displacement_m.iter().flatten().copied().collect();
        let force = problem.internal_force(&flat, &incoming, cx).unwrap();
        for (dof, load) in loads {
            near(force[dof], load, 0.002);
        }
    });
}
