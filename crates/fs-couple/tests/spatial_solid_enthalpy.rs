//! Independent P1 diffusion and affine thermoelastic oracles on the SAME mesh.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_conduction::ConductionProblem;
use fs_conduction::bc::{ThermalBoundary, ThermalBoundaryBuilder};
use fs_conduction::field::ScalarField;
use fs_conduction::material::{CONDUCTIVITY_DIMS, LiquidMassFractionConductivity};
use fs_conduction::mesh::ConductionMesh;
use fs_conduction::transient::enthalpy::heterogeneous::{
    EnthalpyRadiation, HeterogeneousEnthalpyBackwardEuler, ReferenceEnthalpyMaterial,
};
use fs_conduction::transient::enthalpy::{EnthalpyBudget, EnthalpyStepConfig};
use fs_couple::solid_enthalpy::spatial::{
    SpatialMaterialSource, SpatialSolidEnthalpy, SpatialSolidEnthalpyError,
};
use fs_couple::solid_enthalpy::{
    SolidEquilibriumConfig, ThermalExpansionModel, ThermalStiffnessModel,
};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::phase::{
    EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, HeatCapacitySourceQuery,
    ResolvedHeatCapacityCurve, SolidLiquidPhase,
};
use fs_material::state_point::{
    INVERSE_TEMPERATURE_DIMS, MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS,
};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature};
use fs_rep_mesh::TetComplex;
use fs_solid::linear3::{TetAssemblyBudget, TetStaticSolveConfig};
use fs_solid::loads3::VolumeLoad3;
use fs_solver::NewtonKrylovConfig;

const NODES: [[f64; 3]; 4] = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
const FREE: [usize; 6] = [0, 1, 2, 4, 5, 8];

#[path = "spatial_solid_enthalpy/conductivity.rs"]
mod conductivity;
#[path = "spatial_solid_enthalpy/motion.rs"]
mod motion;
#[path = "spatial_solid_enthalpy/radiation.rs"]
mod radiation;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 3,
                kernel_id: 17,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t)
        .unwrap()
}
fn mesh() -> ConductionMesh {
    ConductionMesh::new(TetComplex::from_tets(4, vec![[0, 1, 2, 3]]), NODES.to_vec()).unwrap()
}
fn boundary(mesh: &ConductionMesh) -> ThermalBoundary {
    ThermalBoundaryBuilder::new(mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap()
}
fn thermal_config() -> EnthalpyStepConfig {
    EnthalpyStepConfig {
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-11,
            relative_tolerance: 1e-12,
            linear_restart: 4,
            max_linear_cycles: 8,
            forcing_minimum: 1e-12,
            forcing_maximum: 1e-3,
            ..NewtonKrylovConfig::default()
        },
        max_newton_iterations: 32,
        energy_tolerance_j: 1e-8,
    }
}
fn solid_config() -> SolidEquilibriumConfig {
    SolidEquilibriumConfig {
        static_solve: TetStaticSolveConfig::default(),
        mesh_frame: ContentHash([51; 32]),
        maximum_nodal_displacement_m: 0.1,
    }
}
fn native<'m, 'c>(
    mesh: &'m ConductionMesh,
    curve: &'c EquilibriumEnthalpyPhaseCurve,
    density: f64,
    cx: &Cx<'_>,
) -> HeterogeneousEnthalpyBackwardEuler<'m, 'c> {
    HeterogeneousEnthalpyBackwardEuler::new(
        cx,
        mesh,
        &[ReferenceEnthalpyMaterial {
            curve,
            reference_density_kg_m3: density,
        }],
        &[0],
        EnthalpyBudget {
            max_vertices: 4,
            max_elements: 1,
        },
    )
    .unwrap()
}
fn source<'a>(
    card: &'a MaterialCard,
    curve: &'a EquilibriumEnthalpyPhaseCurve,
    context: &'a QueryPoint,
) -> SpatialMaterialSource<'a> {
    SpatialMaterialSource {
        curve,
        point: ThermalMaterialPoint {
            material: card,
            context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        reference_specific_enthalpy_j_kg: 0.,
        temperature_interval_k: [300., 500.],
        conductivity: PropertyKey::new("thermal_conductivity", CONDUCTIVITY_DIMS),
        conductivity_selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        mechanical_selection: MaterialPropertySelection::SingleClaimOnly,
        stiffness_model: ThermalStiffnessModel::Isotropic,
        expansion_selection: MaterialPropertySelection::SingleClaimOnly,
        expansion_model: ThermalExpansionModel::Isotropic,
    }
}
fn bind<'a, 'm: 'a, 'c: 'a>(
    native: &'a HeterogeneousEnthalpyBackwardEuler<'m, 'c>,
    source: SpatialMaterialSource<'a>,
    cx: &Cx<'_>,
) -> Result<SpatialSolidEnthalpy<'a, 'm, 'c>, SpatialSolidEnthalpyError> {
    SpatialSolidEnthalpy::try_new(
        native,
        vec![source],
        TetAssemblyBudget::standard(),
        solid_config().mesh_frame,
        0.1,
        cx,
    )
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a} != {b}, tolerance {tolerance}"
    );
}

// rho(T) intentionally falls, while the reference inventory remains 6 kg/m3.
// alpha(T) doubles each 100 K, so free strain requires its integral, not alpha(T)*dT.
fn fixture(
    k: f64,
    young: f64,
    alpha: f64,
    cp: f64,
    mechanical_upper: f64,
    interior_young: f64,
) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let heat = PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS);
    let density = PropertyKey::new("density", Density::DIMS);
    let mut claims = ClaimSet::new();
    let mut add = |key: PropertyKey, knots: Vec<(f64, f64)>, high: f64| {
        let dims = key.quantity().dims();
        claims
            .insert_claim(PropertyClaim {
                key,
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots,
                    dims,
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    300.,
                    high,
                ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic P1 identity fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    };
    add(heat.clone(), vec![(300., cp), (500., cp)], 500.);
    add(density.clone(), vec![(300., 6.), (500., 2.)], 500.);
    add(
        PropertyKey::new("thermal_conductivity", CONDUCTIVITY_DIMS),
        vec![(300., k), (500., k)],
        500.,
    );
    add(
        PropertyKey::new("young_modulus", Pressure::DIMS),
        vec![(300., young), (400., interior_young), (500., 0.5 * young)],
        mechanical_upper,
    );
    add(
        PropertyKey::new("poisson_ratio", Dims([0; 6])),
        vec![(300., 0.25), (500., 0.25)],
        500.,
    );
    add(
        PropertyKey::new(
            "linear_thermal_expansion_coefficient",
            INVERSE_TEMPERATURE_DIMS,
        ),
        vec![(300., alpha), (500., 3. * alpha)],
        500.,
    );
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "source substitution fixture".into(),
            phase: "solid".into(),
            process: "synthetic".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &heat,
            density: &density,
            temperature_axis: "T",
            lower: &point(300.),
            upper: &point(500.),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.,
        1e-8,
    )
    .unwrap();
    (card, caloric)
}

#[test]
fn g1_spatial_heat_drives_retained_inertial_motion_and_source_stiffness_work() {
    with_cx(|cx| {
        let (card, caloric) = fixture(2., 60., 1e-5, 10., 500., 45.);
        let mesh = mesh();
        let native = native(&mesh, caloric.curve(), 6., cx);
        let context = point(300.);
        let model = bind(&native, source(&card, caloric.curve(), &context), cx).unwrap();
        let geometry = model.geometry(cx).unwrap();
        let loads = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let fixed: Vec<_> = (0..12).filter(|&d| d != 3).collect();
        let boundary = boundary(&mesh);
        let input = ScalarField::uniform("resolved volumetric heating", 60_000.).unwrap();
        let old = [200., 400., 0., 200.];
        let old_bits = old.map(f64::to_bits);
        let dt = 0.01;
        let mut fields = vec![old.to_vec()];
        for _ in 0..2 {
            let heat = native
                .advance(
                    cx,
                    ConductionProblem {
                        mesh: &mesh,
                        boundary: &boundary,
                        material: model.conductivity().model_for(0).unwrap(),
                        element_materials: Some(model.conductivity()),
                        source: &input,
                    },
                    None,
                    fields.last().unwrap(),
                    dt,
                    thermal_config(),
                )
                .unwrap();
            near(heat.stored_energy_change_j, 100., 1e-8);
            near(heat.energy_residual_j, 0., 1e-8);
            fields.push(heat.specific_enthalpy_j_kg);
        }
        let states: Vec<_> = fields
            .iter()
            .enumerate()
            .map(|(i, h)| model.resolve_at(h, i as f64 * dt, cx).unwrap())
            .collect();
        let problems: Vec<_> = states.iter().map(|state| state.problem(&fixed)).collect();
        let mut displacement = [[0.; 3]; 4];
        let mut velocity = [[0.; 3]; 4];
        displacement[1][0] = 0.001;
        velocity[1][0] = 0.02;
        let mut run = loads
            .begin_transient(
                &problems[0],
                Some(states[0].thermal_strains()),
                0.,
                &displacement,
                &velocity,
                fs_solid::loads3::transient::TransientConfig3 {
                    mass_damping_per_s: 0.4,
                    ..Default::default()
                },
                cx,
            )
            .unwrap();

        // Independent one-DOF P1 oracle: V=1/6, rho=6, nu=1/4 gives
        // M=1/10, K=E/5, F_thermal=E*free_strain/3. The retained source
        // gives E(T)=60-0.15*(T-300), and integral alpha=1e-5*(r+.005*r^2).
        let mass = 0.1;
        let damping = 0.04;
        let (mut q, mut v) = (0.001, 0.02);
        let mut young = 57.;
        let strain = 1e-5 * (20. + 0.005 * 20_f64.powi(2));
        let mut a = (young * strain / 3. - young * q / 5. - damping * v) / mass;
        near(run.state().acceleration[0], a, 1e-12);
        for i in 1..states.len() {
            let previous = run.state().clone();
            let rise = 20. + 10. * i as f64;
            let next_young = 60. - 0.15 * rise;
            let next_strain = 1e-5 * (rise + 0.005 * rise * rise);
            let stiffness = next_young / 5.;
            let b = dt * dt / 4.;
            let predicted_q = q + dt * v + b * a;
            let predicted_v = v + 0.5 * dt * a;
            let next_a =
                (next_young * next_strain / 3. - stiffness * predicted_q - damping * predicted_v)
                    / (mass + 0.5 * dt * damping + b * stiffness);
            let next_q = predicted_q + b * next_a;
            let next_v = predicted_v + 0.5 * dt * next_a;
            let (candidate, report) = states[i]
                .advance_transient_candidate(&states[i - 1], &run, &problems[i], &loads, cx)
                .unwrap();
            near(candidate.state().displacement[0], next_q, 1e-13);
            near(candidate.state().velocity[0], next_v, 1e-12);
            near(candidate.state().acceleration[0], next_a, 1e-10);
            near(
                report.native.stiffness_parameter_work_j,
                0.5 * (next_young - young) / 5. * q * next_q,
                1e-14,
            );
            let elastic = next_young / 15.
                * ((next_q - next_strain).powi(2)
                    + 2. * next_strain.powi(2)
                    + 0.5 * (next_q - 3. * next_strain).powi(2));
            near(
                report.kinetic_plus_elastic_after_j,
                0.5 * mass * next_v.powi(2) + elastic,
                1e-13,
            );
            near(report.physical_balance_error_j, 0., 1e-13);
            assert!(report.stiffness_material_work_j.abs() > 1e-10);
            assert!(report.thermal_material_work_j.abs() > 1e-10);
            assert!(report.native.damping_loss_j > 0.);
            near(
                states[i].samples()[0].phase.temperature_k(),
                300. + rise,
                1e-9,
            );
            let attenuation = (1. + 4. * 2. * dt / 60_f64).powi(-(i as i32));
            near(fields[i][1] / 10., rise + 20. * attenuation, 2e-8);
            near(fields[i][2] / 10., rise - 20. * attenuation, 2e-8);
            assert_eq!(states[i].materials()[0].density_kg_m3, 6.);
            assert!(states[i].samples()[0].phase.bulk_density_kg_m3() < 6.);
            let response = candidate
                .recover_response(solid_config().mesh_frame, 0.1, cx)
                .unwrap();
            near(
                response.stress.elements()[0].stress_mandel_pa()[0],
                1.2 * next_young * next_q - 2. * next_young * next_strain,
                1e-11,
            );
            assert_eq!(candidate.state().step, i);
            assert_eq!(Some(candidate.state().time_s), states[i].time_s());
            assert_eq!(run.state(), &previous);
            (q, v, a, young) = (next_q, next_v, next_a, next_young);
            run = candidate;
        }
        assert_eq!(old.map(f64::to_bits), old_bits);
        assert!(
            states[2]
                .advance_transient_candidate(&states[0], &run, &problems[2], &loads, cx,)
                .is_err()
        );
    });
}

#[test]
fn g4_spatial_inertial_handoff_rejects_stale_clock_fields_and_cancellation_atomically() {
    with_cx(|cx| {
        let (card, caloric) = fixture(2., 60., 1e-5, 10., 500., 45.);
        let mesh = mesh();
        let native = native(&mesh, caloric.curve(), 6., cx);
        let context = point(300.);
        let model = bind(&native, source(&card, caloric.curve(), &context), cx).unwrap();
        let geometry = model.geometry(cx).unwrap();
        let loads = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let fixed: Vec<_> = (0..12).filter(|&d| d != 3).collect();
        let old = [200., 400., 0., 200.];
        let heated = old.map(|h| h + 100.);
        let start = model.resolve_at(&old, 2., cx).unwrap();
        let end = model.resolve_at(&heated, 2.01, cx).unwrap();
        let start_problem = start.problem(&fixed);
        let end_problem = end.problem(&fixed);
        let run = loads
            .begin_transient(
                &start_problem,
                Some(start.thermal_strains()),
                2.,
                &[[0.; 3]; 4],
                &[[0.; 3]; 4],
                Default::default(),
                cx,
            )
            .unwrap();
        let original = run.state().clone();
        let stale_material = model.resolve_at(&heated, 2., cx).unwrap();
        let stale_time = model.resolve_at(&old, 1., cx).unwrap();
        let untimed = model.resolve(&old, cx).unwrap();
        for previous in [&stale_material, &stale_time, &untimed] {
            assert!(matches!(
                end.advance_transient_candidate(previous, &run, &end_problem, &loads, cx,),
                Err(SpatialSolidEnthalpyError::Invalid(_))
            ));
        }
        assert!(matches!(
            end.advance_transient_candidate(&start, &run, &start_problem, &loads, cx,),
            Err(SpatialSolidEnthalpyError::Invalid(_))
        ));
        assert!(
            start
                .advance_transient_candidate(&start, &run, &start_problem, &loads, cx,)
                .is_err()
        );
        assert!(model.resolve_at(&heated, f64::NAN, cx).is_err());
        let changed_supports: Vec<_> = (0..12).filter(|&d| d != 6).collect();
        let changed_problem = end.problem(&changed_supports);
        assert!(
            end.advance_transient_candidate(&start, &run, &changed_problem, &loads, cx,)
                .is_err()
        );
        assert_eq!(run.state(), &original);

        let gate = CancelGate::new();
        gate.request();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            let cancelled = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 3,
                    kernel_id: 17,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            assert!(matches!(
                end.advance_transient_candidate(&start, &run, &end_problem, &loads, &cancelled,),
                Err(SpatialSolidEnthalpyError::Cancelled)
            ));
        });
        assert_eq!(run.state(), &original);
        let (accepted, _) = end
            .advance_transient_candidate(&start, &run, &end_problem, &loads, cx)
            .unwrap();
        let (replayed, _) = end
            .advance_transient_candidate(&start, &run, &end_problem, &loads, cx)
            .unwrap();
        assert_eq!(accepted.state(), replayed.state());
        assert_eq!(accepted.state().time_s, 2.01);
        assert_eq!(run.state(), &original);
        assert_eq!(old, [200., 400., 0., 200.]);
    });
}

#[test]
fn g1_joint_spatial_step_keeps_exact_heat_duration_and_resumes_native_motion() {
    with_cx(|cx| {
        let (card, caloric) = fixture(2., 60., 1e-5, 10., 500., 45.);
        let mesh = mesh();
        let native = native(&mesh, caloric.curve(), 6., cx);
        let context = point(300.);
        let model = bind(&native, source(&card, caloric.curve(), &context), cx).unwrap();
        let geometry = model.geometry(cx).unwrap();
        let loads = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let fixed: Vec<_> = (0..12).filter(|&d| d != 3).collect();
        let boundary = boundary(&mesh);
        let input = ScalarField::uniform("joint-step volumetric heating", 60_000.).unwrap();
        let old = [200., 400., 0., 200.];
        let dt = 0.01_f64;
        let start_time = 2.0_f64;
        assert_ne!(((start_time + dt) - start_time).to_bits(), dt.to_bits());
        let start = model.resolve_at(&old, start_time, cx).unwrap();
        let start_problem = start.problem(&fixed);
        let mut displacement = [[0.; 3]; 4];
        let mut velocity = [[0.; 3]; 4];
        displacement[1][0] = 0.001;
        velocity[1][0] = 0.02;
        let run = loads
            .begin_transient(
                &start_problem,
                Some(start.thermal_strains()),
                start_time,
                &displacement,
                &velocity,
                Default::default(),
                cx,
            )
            .unwrap();
        let original_motion = run.state().clone();
        let first = model
            .advance_transient(
                &old,
                dt,
                &boundary,
                None,
                &input,
                thermal_config(),
                EnthalpyRadiation::None,
                &run,
                &loads,
                cx,
            )
            .unwrap();
        assert_eq!(first.mechanical.native.step_size_s.to_bits(), dt.to_bits());
        assert_eq!(first.state.time_s(), Some(start_time + dt));
        assert_eq!(first.motion.state().time_s, start_time + dt);
        assert_eq!(first.motion.state().step, 1);
        near(
            first.thermal.conduction().stored_energy_change_j,
            100.,
            1e-8,
        );
        near(first.thermal.physical_energy_residual_j(), 0., 1e-8);
        near(first.state.samples()[0].phase.temperature_k(), 330., 1e-9);
        assert_eq!(run.state(), &original_motion);

        let first_problem = first.state.problem(&fixed);
        let restored = first.motion.bind(&first_problem, cx).unwrap();
        // Compare with the original native owner given the ACTUAL heat duration,
        // not a duration reconstructed by subtracting rounded endpoint clocks.
        let (native_motion, native_work) = run
            .advance_thermal_candidate(
                &first_problem,
                first.state.thermal_strains(),
                &loads,
                dt,
                solid_config().mesh_frame,
                0.1,
                cx,
            )
            .unwrap();
        assert_eq!(restored.state(), native_motion.state());
        assert_eq!(first.mechanical, native_work);

        let original_h = first.thermal.conduction().specific_enthalpy_j_kg.clone();
        let second = model
            .advance_transient(
                &original_h,
                dt,
                &boundary,
                None,
                &input,
                thermal_config(),
                EnthalpyRadiation::None,
                &restored,
                &loads,
                cx,
            )
            .unwrap();
        let checkpoint = first.motion.clone();
        let replay = checkpoint.bind(&first_problem, cx).unwrap();
        let repeated = model
            .advance_transient(
                &original_h,
                dt,
                &boundary,
                None,
                &input,
                thermal_config(),
                EnthalpyRadiation::None,
                &replay,
                &loads,
                cx,
            )
            .unwrap();
        assert_eq!(second.motion.state(), repeated.motion.state());
        assert_eq!(second.mechanical, repeated.mechanical);
        assert_eq!(
            second.thermal.conduction().specific_enthalpy_j_kg,
            repeated.thermal.conduction().specific_enthalpy_j_kg
        );
        assert_eq!(second.motion.state().step, 2);
        assert_eq!(second.motion.state().time_s, (start_time + dt) + dt);
        assert_eq!(second.state.time_s(), Some(second.motion.state().time_s));
        assert_eq!(second.mechanical.native.step_size_s.to_bits(), dt.to_bits());
        near(
            second.thermal.conduction().stored_energy_change_j,
            100.,
            1e-8,
        );
        near(second.state.samples()[0].phase.temperature_k(), 340., 1e-9);
        assert_eq!(second.state.materials()[0].density_kg_m3, 6.);
        assert_ne!(
            first.state.materials()[0].material_state_identity,
            second.state.materials()[0].material_state_identity
        );
        near(second.mechanical.physical_balance_error_j, 0., 1e-13);
        assert_eq!(
            first.thermal.conduction().specific_enthalpy_j_kg,
            original_h
        );
        assert_eq!(restored.state(), first.motion.state());
        assert_eq!(old, [200., 400., 0., 200.]);
    });
}

#[test]
fn g4_joint_spatial_heat_and_motion_refuse_without_publishing_either_history() {
    with_cx(|cx| {
        let (card, caloric) = fixture(2., 60., 1e-5, 10., 500., 45.);
        let mesh = mesh();
        let native = native(&mesh, caloric.curve(), 6., cx);
        let context = point(300.);
        let model = bind(&native, source(&card, caloric.curve(), &context), cx).unwrap();
        let geometry = model.geometry(cx).unwrap();
        let loads = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let fixed: Vec<_> = (0..12).filter(|&d| d != 3).collect();
        let old = [0.; 4];
        let start = model.resolve_at(&old, 2., cx).unwrap();
        let problem = start.problem(&fixed);
        let zero = [[0.; 3]; 4];
        let tiny = loads
            .begin_transient(
                &problem,
                Some(start.thermal_strains()),
                2.,
                &zero,
                &zero,
                fs_solid::loads3::transient::TransientConfig3 {
                    maximum_displacement_m: 1e-12,
                    ..Default::default()
                },
                cx,
            )
            .unwrap();
        let run = loads
            .begin_transient(
                &problem,
                Some(start.thermal_strains()),
                2.,
                &zero,
                &zero,
                Default::default(),
                cx,
            )
            .unwrap();
        let original = run.state().clone();
        let original_tiny = tiny.state().clone();
        let boundary = boundary(&mesh);
        let input = ScalarField::uniform("joint-step heating", 60_000.).unwrap();
        // Heat is solved successfully before the native displacement gate fails.
        let refused = model.advance_transient(
            &old,
            0.01,
            &boundary,
            None,
            &input,
            thermal_config(),
            EnthalpyRadiation::None,
            &tiny,
            &loads,
            cx,
        );
        assert!(matches!(
            refused,
            Err(SpatialSolidEnthalpyError::Solid(
                fs_couple::solid_enthalpy::SolidEnthalpyError::Transient(
                    fs_solid::loads3::transient::TransientError3::Invalid(
                        "candidate exceeds displacement allowance"
                    )
                )
            ))
        ));
        assert_eq!(tiny.state(), &original_tiny);
        assert!(
            model
                .advance_transient(
                    &[100.; 4],
                    0.01,
                    &boundary,
                    None,
                    &input,
                    thermal_config(),
                    EnthalpyRadiation::None,
                    &run,
                    &loads,
                    cx,
                )
                .is_err()
        );
        assert!(
            model
                .advance_transient(
                    &old,
                    f64::MIN_POSITIVE,
                    &boundary,
                    None,
                    &input,
                    thermal_config(),
                    EnthalpyRadiation::None,
                    &run,
                    &loads,
                    cx,
                )
                .is_err()
        );
        let gate = CancelGate::new();
        gate.request();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            let cancelled = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 3,
                    kernel_id: 17,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            assert!(matches!(
                model.advance_transient(
                    &old,
                    0.01,
                    &boundary,
                    None,
                    &input,
                    thermal_config(),
                    EnthalpyRadiation::None,
                    &run,
                    &loads,
                    &cancelled,
                ),
                Err(SpatialSolidEnthalpyError::Cancelled)
            ));
        });
        assert_eq!(run.state(), &original);
        assert_eq!(old, [0.; 4]);
        let retry = model
            .advance_transient(
                &old,
                0.01,
                &boundary,
                None,
                &input,
                thermal_config(),
                EnthalpyRadiation::None,
                &run,
                &loads,
                cx,
            )
            .unwrap();
        assert!(retry.motion.state().displacement[0] > 1e-12);
        assert_eq!(retry.motion.state().step, 1);
        near(
            retry.thermal.conduction().stored_energy_change_j,
            100.,
            1e-8,
        );
        near(retry.mechanical.physical_balance_error_j, 0., 1e-13);
        assert_eq!(run.state(), &original);
        assert_eq!(old, [0.; 4]);
    });
}

#[test]
fn source_only_substitution_changes_native_diffusion_affine_motion_and_constrained_stress() {
    with_cx(|cx| {
        for (k, young, alpha, cp) in [(2., 1000., 1e-5, 10.), (8., 2000., 2e-5, 20.)] {
            let (card, caloric) = fixture(k, young, alpha, cp, 500., 0.75 * young);
            let mesh = mesh();
            let native = native(&mesh, caloric.curve(), 6., cx);
            let context = point(300.);
            let model = bind(&native, source(&card, caloric.curve(), &context), cx).unwrap();
            let geometry = model.geometry(cx).unwrap();
            let loads = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
            let old = [100. * cp, 120. * cp, 80. * cp, 100. * cp];
            let original_bits = old.map(f64::to_bits);
            let boundary = boundary(&mesh);
            let input = ScalarField::uniform("heat source", 60.).unwrap();
            let result = model
                .advance(
                    &old,
                    1.,
                    &boundary,
                    None,
                    &input,
                    thermal_config(),
                    &loads,
                    &FREE,
                    solid_config(),
                    cx,
                )
                .unwrap();
            // Unit-tet P1: K=k/6 G'G, m_i=rho/24. (0,1,-1,0)
            // is an exact diffusion eigenvector; uniform source raises the mean.
            let rise = 100. + 60. / (6. * cp);
            let attenuation = 1. / (1. + 4. * k / (6. * cp));
            for (actual, perturbation) in result.thermal.temperature.iter().zip([0., 20., -20., 0.])
            {
                near(*actual, 300. + rise + perturbation * attenuation, 2e-8);
            }
            near(result.thermal.stored_energy_change_j, 10., 1e-8);
            near(result.thermal.energy_residual_j, 0., 1e-8);
            let free_strain = alpha * (rise + 0.005 * rise * rise);
            for (actual, reference) in result.equilibrium.deformed_mesh.nodes_m().iter().zip(NODES)
            {
                for c in 0..3 {
                    near(actual[c], reference[c] * (1. + free_strain), 2e-9);
                }
            }
            for stress in result.equilibrium.solution.stress.elements()[0].stress_mandel_pa() {
                near(*stress, 0., 1e-7);
            }
            near(result.state.materials()[0].density_kg_m3, 6., 0.);
            assert!(result.state.samples()[0].phase.bulk_density_kg_m3() < 6.);
            assert!(
                model
                    .conductivity()
                    .receipts()
                    .iter()
                    .all(|r| card.claims().claim(r.selected).is_some())
            );
            let all_fixed: Vec<_> = (0..12).collect();
            let clamped = model
                .advance(
                    &old,
                    1.,
                    &boundary,
                    None,
                    &input,
                    thermal_config(),
                    &loads,
                    &all_fixed,
                    solid_config(),
                    cx,
                )
                .unwrap();
            let current_young = young * (1. - 0.0025 * rise);
            let stress = -current_young * free_strain / (1. - 2. * 0.25);
            for &actual in
                &clamped.equilibrium.solution.stress.elements()[0].stress_mandel_pa()[..3]
            {
                near(actual, stress, 1e-7);
            }
            assert_eq!(clamped.equilibrium.deformed_mesh.nodes_m(), &NODES);
            assert_eq!(old.map(f64::to_bits), original_bits);
            assert_eq!(native.reference_nodal_masses_kg(), &[0.25; 4]);
            assert!(core::ptr::eq(model.thermal().mesh(), &mesh));
        }
    });
}

#[test]
fn enthalpy_centroid_is_explicit_and_is_not_mean_nodal_temperature() {
    with_cx(|cx| {
        let (card, _) = fixture(2., 1000., 1e-5, 10., 500., 750.);
        let curve = EquilibriumEnthalpyPhaseCurve::try_single_phase(
            card.content_hash(),
            SolidLiquidPhase::Solid,
            [(0., 300.), (1000., 400.), (3000., 500.)]
                .into_iter()
                .map(|(h, t)| EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: h,
                    temperature_k: t,
                    liquid_mass_fraction: 0.,
                    bulk_density_kg_m3: 6.,
                })
                .collect(),
        )
        .unwrap();
        let mesh = mesh();
        let native = native(&mesh, &curve, 6., cx);
        let context = point(300.);
        let model = bind(&native, source(&card, &curve, &context), cx).unwrap();
        let state = model.resolve(&[0., 0., 0., 3000.], cx).unwrap();
        near(state.samples()[0].phase.specific_enthalpy_j_kg(), 750., 0.);
        near(state.samples()[0].phase.temperature_k(), 375., 0.);
        assert_ne!(state.samples()[0].phase.temperature_k(), 350.); // mean nodal T
    });
}

#[test]
fn every_corner_refuses_melt_and_source_escape_even_when_centroid_remains_admissible() {
    with_cx(|cx| {
        let (card, _) = fixture(2., 1000., 1e-5, 10., 500., 750.);
        let curve = EquilibriumEnthalpyPhaseCurve::try_new(
            card.content_hash(),
            vec![
                EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: 0.,
                    temperature_k: 300.,
                    liquid_mass_fraction: 0.,
                    bulk_density_kg_m3: 6.,
                },
                EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: 1000.,
                    temperature_k: 400.,
                    liquid_mass_fraction: 0.,
                    bulk_density_kg_m3: 4.,
                },
                EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: 1200.,
                    temperature_k: 400.,
                    liquid_mass_fraction: 1.,
                    bulk_density_kg_m3: 3.,
                },
            ],
        )
        .unwrap();
        let mesh = mesh();
        let native = native(&mesh, &curve, 6., cx);
        let context = point(300.);
        let model = bind(&native, source(&card, &curve, &context), cx).unwrap();
        assert!(matches!(
            model.resolve(&[0., 0., 0., 1100.], cx),
            Err(SpatialSolidEnthalpyError::Vertex { vertex: 3, .. })
        ));
        let mut restricted = source(&card, &curve, &context);
        restricted.temperature_interval_k = [300., 350.];
        let restricted = bind(&native, restricted, cx).unwrap();
        assert!(matches!(
            restricted.resolve(&[0., 0., 0., 600.], cx),
            Err(SpatialSolidEnthalpyError::Vertex { vertex: 3, .. })
        ));
        let geometry = model.geometry(cx).unwrap();
        let loads = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let initial = [950.; 4];
        let bits = initial.map(f64::to_bits);
        // Uniform heat crosses into the native latent plateau. The native heat
        // owner can solve it; this solid candidate must refuse its partial melt.
        let hot = ScalarField::uniform("heat source", 600.).unwrap();
        assert!(matches!(
            model.advance(
                &initial,
                1.,
                &boundary(&mesh),
                None,
                &hot,
                thermal_config(),
                &loads,
                &FREE,
                solid_config(),
                cx
            ),
            Err(SpatialSolidEnthalpyError::Vertex { .. })
        ));
        assert_eq!(initial.map(f64::to_bits), bits);
    });
}

#[test]
fn native_chart_inventory_and_continuous_source_support_are_admission_requirements() {
    with_cx(|cx| {
        let (card, caloric) = fixture(2., 1000., 1e-5, 10., 500., 750.);
        let (other, foreign) = fixture(4., 1000., 1e-5, 10., 500., 750.);
        let mesh = mesh();
        let native_good = native(&mesh, caloric.curve(), 6., cx);
        let native_wrong_mass = native(&mesh, caloric.curve(), 3., cx);
        let context = point(300.);
        assert!(
            bind(
                &native_wrong_mass,
                source(&card, caloric.curve(), &context),
                cx
            )
            .is_err()
        );
        assert!(bind(&native_good, source(&other, foreign.curve(), &context), cx).is_err());
        assert!(bind(&native_good, source(&other, caloric.curve(), &context), cx).is_err());
        let native_scaled = native(&mesh, caloric.curve(), 6., cx)
            .with_phase_conductivity(
                cx,
                vec![LiquidMassFractionConductivity::declared(2., 2.).unwrap()],
            )
            .unwrap();
        assert!(bind(&native_scaled, source(&card, caloric.curve(), &context), cx).is_err());
        // Valid reference and two endpoint constants are insufficient when the
        // source ends early, or a source's interior knot has an invalid modulus.
        for (high, interior) in [(350., 750.), (500., -1.)] {
            let (card, caloric) = fixture(2., 1000., 1e-5, 10., high, interior);
            let native = native(&mesh, caloric.curve(), 6., cx);
            assert!(bind(&native, source(&card, caloric.curve(), &context), cx).is_err());
        }
    });
}

#[test]
fn a_refused_mechanical_endpoint_does_not_publish_or_modify_the_heat_history() {
    with_cx(|cx| {
        let (card, caloric) = fixture(2., 1000., 1e-5, 10., 500., 750.);
        let mesh = mesh();
        let native = native(&mesh, caloric.curve(), 6., cx);
        let context = point(300.);
        let model = bind(&native, source(&card, caloric.curve(), &context), cx).unwrap();
        let geometry = model.geometry(cx).unwrap();
        let loads = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let old = [500.; 4];
        let bits = old.map(f64::to_bits);
        let input = ScalarField::uniform("heat source", 60.).unwrap();
        let tiny = SolidEquilibriumConfig {
            maximum_nodal_displacement_m: 1e-9,
            ..solid_config()
        };
        assert!(
            model
                .advance(
                    &old,
                    1.,
                    &boundary(&mesh),
                    None,
                    &input,
                    thermal_config(),
                    &loads,
                    &FREE,
                    tiny,
                    cx
                )
                .is_err()
        );
        assert_eq!(old.map(f64::to_bits), bits);
        // The same original state remains usable after refusal.
        assert!(
            model
                .advance(
                    &old,
                    1.,
                    &boundary(&mesh),
                    None,
                    &input,
                    thermal_config(),
                    &loads,
                    &FREE,
                    solid_config(),
                    cx
                )
                .is_ok()
        );
        let gate = CancelGate::new();
        gate.request();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            let cancelled = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 3,
                    kernel_id: 17,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            assert!(matches!(
                model.advance(
                    &old,
                    1.,
                    &boundary(&mesh),
                    None,
                    &input,
                    thermal_config(),
                    &loads,
                    &FREE,
                    solid_config(),
                    &cancelled
                ),
                Err(SpatialSolidEnthalpyError::Cancelled)
            ));
            assert!(matches!(
                bind(
                    &native,
                    source(&card, caloric.curve(), &context),
                    &cancelled
                ),
                Err(SpatialSolidEnthalpyError::Cancelled)
            ));
        });
        assert_eq!(old.map(f64::to_bits), bits);
    });
}

#[test]
fn heterogeneous_assignment_order_drives_distinct_elements_without_a_second_heat_store() {
    with_cx(|cx| {
        let (first, first_curve) = fixture(2., 1000., 1e-5, 10., 500., 750.);
        let (second, second_curve) = fixture(8., 2000., 3e-5, 20., 500., 1500.);
        let mut nodes = NODES.to_vec();
        nodes.extend(NODES.map(|p| [p[0] + 3., p[1], p[2]]));
        let mesh = ConductionMesh::new(
            TetComplex::from_tets(8, vec![[0, 1, 2, 3], [4, 5, 6, 7]]),
            nodes,
        )
        .unwrap();
        let native = HeterogeneousEnthalpyBackwardEuler::new(
            cx,
            &mesh,
            &[
                ReferenceEnthalpyMaterial {
                    curve: first_curve.curve(),
                    reference_density_kg_m3: 6.,
                },
                ReferenceEnthalpyMaterial {
                    curve: second_curve.curve(),
                    reference_density_kg_m3: 6.,
                },
            ],
            &[1, 0],
            EnthalpyBudget {
                max_vertices: 8,
                max_elements: 2,
            },
        )
        .unwrap();
        let context = point(300.);
        let model = SpatialSolidEnthalpy::try_new(
            &native,
            vec![
                source(&first, first_curve.curve(), &context),
                source(&second, second_curve.curve(), &context),
            ],
            TetAssemblyBudget::standard(),
            solid_config().mesh_frame,
            0.1,
            cx,
        )
        .unwrap();
        let geometry = model.geometry(cx).unwrap();
        let loads = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let mut fixed = FREE.to_vec();
        fixed.extend(FREE.map(|dof| dof + 12));
        let input =
            ScalarField::nodal("local heating", 8, vec![60., 60., 60., 60., 0., 0., 0., 0.])
                .unwrap();
        let result = model
            .advance(
                &[0.; 8],
                1.,
                &boundary(&mesh),
                None,
                &input,
                thermal_config(),
                &loads,
                &fixed,
                solid_config(),
                cx,
            )
            .unwrap();
        assert_eq!(
            result
                .state
                .samples()
                .iter()
                .map(|s| s.material_index)
                .collect::<Vec<_>>(),
            vec![1, 0]
        );
        assert_eq!(
            result.state.samples()[0].elastic.card_identity(),
            second.content_hash()
        );
        assert_eq!(
            result.state.samples()[1].elastic.card_identity(),
            first.content_hash()
        );
        for &t in &result.thermal.temperature[..4] {
            near(t, 300.5, 1e-8);
        }
        for &t in &result.thermal.temperature[4..] {
            near(t, 300., 1e-8);
        }
        let strain = 3e-5 * (0.5 + 0.005 * 0.5 * 0.5);
        near(
            result.equilibrium.deformed_mesh.nodes_m()[1][0],
            1. + strain,
            1e-9,
        );
        for vertex in 4..8 {
            assert_eq!(
                result.equilibrium.deformed_mesh.nodes_m()[vertex],
                mesh.positions()[vertex]
            );
        }
        near(result.thermal.stored_energy_change_j, 10., 1e-8);
    });
}

#[test]
fn changing_only_conductivity_propagates_through_spatial_heat_into_actual_thermal_stress() {
    with_cx(|cx| {
        // Equal-volume mirrored tetrahedra. The odd apex mode (-1,0,0,0,1)
        // has K/m eigenvalue 12*k/rho; the three shared-face nodes stay at 400 K.
        let mut nodes = NODES.to_vec();
        nodes.push([2. / 3.; 3]);
        let mesh = ConductionMesh::new(
            TetComplex::from_tets(5, vec![[0, 1, 2, 3], [4, 1, 2, 3]]),
            nodes,
        )
        .unwrap();
        let context = point(300.);
        let mut responses = Vec::new();
        for k in [1., 8.] {
            let (card, caloric) = fixture(k, 1000., 1e-5, 10., 500., 750.);
            let native = HeterogeneousEnthalpyBackwardEuler::new(
                cx,
                &mesh,
                &[ReferenceEnthalpyMaterial {
                    curve: caloric.curve(),
                    reference_density_kg_m3: 6.,
                }],
                &[0, 0],
                EnthalpyBudget {
                    max_vertices: 5,
                    max_elements: 2,
                },
            )
            .unwrap();
            let model = bind(&native, source(&card, caloric.curve(), &context), cx).unwrap();
            let geometry = model.geometry(cx).unwrap();
            let loads = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
            let old = [800., 1000., 1000., 1000., 1200.];
            let result = model
                .advance(
                    &old,
                    1.,
                    &boundary(&mesh),
                    None,
                    &ScalarField::Uniform(0.),
                    thermal_config(),
                    &loads,
                    &FREE,
                    solid_config(),
                    cx,
                )
                .unwrap();
            let attenuation = 1. / (1. + 12. * k / 60.);
            near(
                result.thermal.temperature[0],
                400. - 20. * attenuation,
                2e-8,
            );
            near(
                result.thermal.temperature[4],
                400. + 20. * attenuation,
                2e-8,
            );
            for &t in &result.thermal.temperature[1..4] {
                near(t, 400., 2e-8);
            }
            near(
                result.state.samples()[0].phase.temperature_k(),
                400. - 5. * attenuation,
                2e-8,
            );
            near(
                result.state.samples()[1].phase.temperature_k(),
                400. + 5. * attenuation,
                2e-8,
            );
            near(result.thermal.stored_energy_change_j, 0., 1e-8);
            let stress: Vec<_> = result
                .equilibrium
                .solution
                .stress
                .elements()
                .iter()
                .flat_map(|element| *element.stress_mandel_pa())
                .collect();
            assert!(stress.iter().any(|s: &f64| s.abs() > 1e-5));
            responses.push((stress, result.equilibrium.deformed_mesh.nodes_m().to_vec()));
        }
        assert!(
            responses[0]
                .0
                .iter()
                .zip(&responses[1].0)
                .any(|(a, b)| (a - b).abs() > 1e-5)
        );
        assert!(
            responses[0]
                .1
                .iter()
                .flatten()
                .zip(responses[1].1.iter().flatten())
                .any(|(a, b)| (a - b).abs() > 1e-6)
        );
    });
}
