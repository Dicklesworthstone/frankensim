//! Finite source gas loads the real elastic modes and shares their volume work.
#[path = "gas_motion/reactive.rs"]
mod reactive;
#[path = "gas_motion/exchange.rs"]
mod exchange;

use super::*;
use fs_couple::cavity_phs::modal::{GasValveMotion, ModalGasError, ModalGasSolve};
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter};
use fs_material::thermochemical_gas::ThermochemicalGasModel;
use fs_material::thermochemical_gas::energy::{GasEnergySolve, parcel::ClosedGasParcel};
use fs_math::c64::C64;
use fs_modal::time::{ModalAcousticState, ModalAcousticTimeBudget};
use fs_qty::MolarMass;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_solid::reduce::thermal::motion::ThermalValveModalMotion;
use fs_solver::NewtonKrylovConfig;
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::{
    Composition, CompositionBasis, NASA9_LAW_ID_V1, NASA9_LAW_VERSION_V1,
    NASA9_STATE_SCHEMA_VERSION_V1, Nasa9RegionV1, Nasa9StandardStateModelV1,
    required_nasa9_parameter_dims,
};

// Reuse existing transport-source construction; only the caloric capacities
// differ here so composition really changes the acoustic bulk modulus.
#[allow(dead_code)]
#[path = "../../material_exchange/gas.rs"]
mod gas_fixture;

const RATE: u32 = 2_000;
const VOLUME: f64 = 0.001;
const EXTERIOR: f64 = 100_000.0;
const AREA: f64 = 0.001;

fn gas(diatomic_fraction: f64) -> ClosedGasParcel {
    let mut components = Vec::new();
    let mut transport = Vec::new();
    for (name, mass, cp_r) in [("A", 0.028, 3.5), ("B", 0.004, 2.5)] {
        let (original, source) = gas_fixture::component(name, mass, 0.0);
        let mut coefficients = [0.0; 9];
        coefficients[2] = cp_r;
        let mut parameters = required_nasa9_parameter_dims()
            .into_iter()
            .take(9)
            .zip(coefficients)
            .map(|((key, dims), value)| (key.to_owned(), LawParameter { value, dims }))
            .collect::<std::collections::BTreeMap<_, _>>();
        parameters.insert(
            "reference_pressure".into(),
            LawParameter {
                value: EXTERIOR,
                dims: Pressure::DIMS,
            },
        );
        let region = Nasa9RegionV1::from_card(ConstitutiveModelCard {
            law: LawId(NASA9_LAW_ID_V1.into()),
            law_version: NASA9_LAW_VERSION_V1,
            parameters,
            state_schema_version: NASA9_STATE_SCHEMA_VERSION_V1,
            initial_state: InitialStatePolicy::ZeroInternalState,
            validity: ValidityDomain::unconstrained().with("T", 200.0, 2000.0),
            sources: vec![],
            provenance: Provenance {
                source: "synthetic monoatomic/diatomic caloric fixture".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
        components.push(
            Nasa9StandardStateModelV1::new(
                original.species().clone(),
                MolarMass::new(mass),
                original.convention().clone(),
                vec![region],
            )
            .unwrap(),
        );
        transport.push(source);
    }
    let composition = Composition::new(
        CompositionBasis::MoleFraction,
        vec![diatomic_fraction, 1.0 - diatomic_fraction],
    )
    .unwrap();
    let model = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(components, composition).unwrap(),
        transport,
        None,
    )
    .unwrap();
    let quantity = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
    let point = QueryPoint::new()
        .with_quantity(
            "temperature",
            quantity(QuantityKind::AbsoluteTemperature),
            400.0,
        )
        .unwrap()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), EXTERIOR)
        .unwrap();
    ClosedGasParcel::new(model, point, VOLUME).unwrap()
}

fn motion(
    young: f64,
    velocity_m_s: f64,
    loss: f64,
    area: f64,
    rate: u32,
    cx: &Cx<'_>,
) -> ThermalValveModalMotion {
    let fixed: Vec<_> = (0..24).filter(|&d| d != 21).collect();
    let mut pressure = [[0.0; 3]; 8];
    pressure[7][0] = area;
    let mut input = reduction(&fixed, &pressure);
    input.equilibrium_nodal_forces_n = &[];
    input.loss_factor = loss;
    let context = point(300.0);
    let (card, caloric) = material(young, 1e-5, 10.0, 0.0, false);
    let cells: Vec<_> = (0..6).map(|id| (id, 0.0)).collect();
    let thermal = initial(caloric.curve(), &cells, 1.0);
    let source = SolidEnthalpy::try_new(
        &thermal,
        &CUBE,
        &CUBE_TETS,
        (0..6)
            .map(|id| super::super::source(id, &card, caloric.curve(), &context))
            .collect(),
        TetAssemblyBudget::standard(),
        0.1,
        config().mesh_frame,
        cx,
    )
    .unwrap();
    let valve = source
        .resolve(&thermal, cx)
        .unwrap()
        .reduce_pressure_valve(&input, cx)
        .unwrap();
    assert_eq!(valve.valve().source_modes().len(), 1);
    let shape = valve.valve().source_modes()[0].phi[0];
    ThermalValveModalMotion::new(
        valve,
        rate,
        &[ModalAcousticState {
            displacement_m_sqrt_kg: 0.0,
            velocity_m_sqrt_kg_per_s: velocity_m_s / shape,
        }],
        &[C64::new(2.0 * shape, 0.0)],
        ModalAcousticTimeBudget::audible_reference(),
        cx,
    )
    .unwrap()
}

fn controls() -> ModalGasSolve {
    ModalGasSolve {
        gas: GasEnergySolve {
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 1800.0,
            energy_tolerance_j_kg: 1e-6,
            maximum_iterations: 128,
        },
        volume_interval_m3: [0.0007, 0.0015],
        volume_tolerance_m3: 1e-14,
        energy_tolerance_j: 1e-8,
        maximum_trials: 64,
        maximum_newton_iterations: 16,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-12,
            relative_tolerance: 0.0,
            linear_restart: 1,
            max_linear_cycles: 2,
            ..Default::default()
        },
    }
}

fn energy(motion: &ThermalValveModalMotion) -> f64 {
    motion
        .model()
        .states()
        .iter()
        .zip(motion.model().modes())
        .map(|(state, mode)| {
            0.5 * (state.velocity_m_sqrt_kg_per_s.powi(2)
                + (mode.angular_frequency_rad_s * state.displacement_m_sqrt_kg).powi(2))
        })
        .sum()
}

#[test]
fn equilibrium_keeps_the_original_gas_and_solid_states_on_one_clock() {
    with_cx(|cx| {
        let original_gas = gas(0.9);
        let original_motion = motion(1000.0, 0.0, 0.0, AREA, RATE, cx);
        let identity = original_motion.valve().identity();
        let states = original_motion.model().states().to_vec();
        let mut run = GasValveMotion::new(original_gas.clone(), original_motion, cx).unwrap();
        for step in 1..=8 {
            let report = run.step(0.0, EXTERIOR, controls(), cx).unwrap();
            assert_eq!(run.gas(), &original_gas);
            assert_eq!(run.motion().model().states(), states);
            assert_eq!(run.motion().valve().identity(), identity);
            assert_eq!(run.motion().accepted_steps(), step);
            assert_eq!(report.motion.time_s, step as f64 / RATE as f64);
            assert_eq!(report.motion.swept_volume_change_m3, 0.0);
            assert_eq!(report.total_energy_residual_j, 0.0);
        }
    });
}

#[test]
fn gas_heat_drives_actual_wall_sweep_and_closes_storage_loss_and_exterior_work() {
    with_cx(|cx| {
        let parcel = gas(0.9);
        let initial_u = parcel.state().specific_internal_energy_j_kg();
        let initial_mass = parcel.mass_kg();
        let solid = motion(1000.0, 0.0, 0.1, AREA, RATE, cx);
        let initial_energy = energy(&solid);
        let identity = solid.valve().identity();
        let mut run = GasValveMotion::new(parcel, solid, cx).unwrap();
        let report = run.step(1.0, EXTERIOR, controls(), cx).unwrap();
        let physical = run.motion().nodal_motion(cx).unwrap();
        assert!(physical.displacement_m[7][0] > 0.0);
        assert!(physical.velocity_m_s[7][0] > 0.0);
        assert!(run.gas().state().state().temperature > 400.0);
        assert!(run.gas().volume_m3() > VOLUME);
        assert!(report.motion.modal.viscous_dissipation_j > 0.0);
        let sweep = AREA * physical.displacement_m[7][0];
        near(run.gas().volume_m3() - VOLUME, sweep, 1e-14);
        near(report.motion.swept_volume_change_m3, sweep, 1e-16);
        near(report.exterior_pressure_work_j, EXTERIOR * sweep, 1e-11);
        let measured = initial_mass
            * (run.gas().state().specific_internal_energy_j_kg() - initial_u)
            + energy(run.motion())
            - initial_energy
            + report.motion.modal.viscous_dissipation_j
            + EXTERIOR * sweep
            - 1.0;
        assert!(measured.abs() <= controls().energy_tolerance_j);
        assert!(report.pressure_work_residual_j.abs() <= controls().energy_tolerance_j);
        assert!(report.total_energy_residual_j.abs() <= controls().energy_tolerance_j);
        assert_eq!(run.gas().mass_kg(), initial_mass);
        assert_eq!(run.motion().valve().identity(), identity);
        // This is the explicitly supplied diagnostic transfer, not a claim of
        // automatically rebuilt exterior radiation when the cavity gas changes.
        near(
            report.motion.modal.observer_pressure_pa,
            2.0 * physical.velocity_m_s[7][0],
            1e-10,
        );
    });
}

#[test]
fn gas_composition_and_solid_elastic_data_change_the_same_continuing_motion() {
    let run = |fraction, young| {
        with_cx(|cx| {
            let mut run =
                GasValveMotion::new(gas(fraction), motion(young, 0.02, 0.0, AREA, RATE, cx), cx)
                    .unwrap();
            let omega = run.motion().model().modes()[0].angular_frequency_rad_s;
            let gamma = run.gas().state().state().gamma;
            let mut observation = 0.0;
            for _ in 0..32 {
                let report = run.step(0.0, EXTERIOR, controls(), cx).unwrap();
                observation = report.motion.modal.observer_pressure_pa;
                assert!(report.total_energy_residual_j.abs() <= controls().energy_tolerance_j);
            }
            (
                run.motion().nodal_motion(cx).unwrap().velocity_m_s[7][0],
                omega,
                gamma,
                observation,
            )
        })
    };
    let diatomic = run(0.9, 1000.0);
    let monoatomic = run(0.1, 1000.0);
    let stiff = run(0.9, 2000.0);
    assert!(monoatomic.2 > diatomic.2);
    assert_eq!(monoatomic.1, diatomic.1);
    assert!(stiff.1 > diatomic.1);
    // The higher source gamma increases the restoring gas bulk modulus at
    // the same initial P,T,V. A different molar mass alone would not prove this.
    assert!(monoatomic.0 < diatomic.0 - 1e-8);
    assert!(stiff.0 < diatomic.0 - 1e-8);
    assert!((monoatomic.3 - diatomic.3).abs() > 1e-8);
    assert!((stiff.3 - diatomic.3).abs() > 1e-8);
}

#[test]
fn failed_material_work_volume_and_cancellation_trials_publish_neither_owner() {
    with_cx(|cx| {
        let mut run =
            GasValveMotion::new(gas(0.9), motion(1000.0, 0.0, 0.0, AREA, RATE, cx), cx).unwrap();
        let original = run.clone();
        let unchanged = |run: &GasValveMotion| {
            assert_eq!(run.gas(), original.gas());
            assert_eq!(
                run.motion().model().states(),
                original.motion().model().states()
            );
            assert_eq!(
                run.motion().accepted_steps(),
                original.motion().accepted_steps()
            );
            assert_eq!(
                run.motion().valve().identity(),
                original.motion().valve().identity()
            );
        };
        let mut cold = controls();
        cold.gas.maximum_temperature_k = 400.01;
        assert!(matches!(
            run.step(1.0, EXTERIOR, cold, cx),
            Err(ModalGasError::Gas(_))
        ));
        unchanged(&run);
        let mut limited = controls();
        limited.maximum_trials = 1;
        assert!(matches!(
            run.step(1.0, EXTERIOR, limited, cx),
            Err(ModalGasError::TrialBudget)
        ));
        unchanged(&run);
        let mut narrow = controls();
        narrow.volume_interval_m3 = [VOLUME - 1e-16, VOLUME + 1e-16];
        assert!(run.step(1.0, EXTERIOR, narrow, cx).is_err());
        unchanged(&run);
        let gate = CancelGate::new_clock_free();
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
                run.step(1.0, EXTERIOR, controls(), &cancelled),
                Err(ModalGasError::Cancelled)
            ));
        });
        unchanged(&run);
        let mut clean = original;
        run.step(1.0, EXTERIOR, controls(), cx).unwrap();
        clean.step(1.0, EXTERIOR, controls(), cx).unwrap();
        assert_eq!(run.gas(), clean.gas());
        assert_eq!(
            run.motion().model().states(),
            clean.motion().model().states()
        );
        assert_eq!(run.motion().accepted_steps(), 1);
    });
}

#[test]
fn gas_loaded_sampling_gate_rejects_a_clock_that_admits_the_bare_wall() {
    with_cx(|cx| {
        let bare = motion(1000.0, 0.0, 0.0, 1.0, RATE, cx);
        let omega_dt =
            bare.model().modes()[0].angular_frequency_rad_s * bare.model().sample_period_s();
        assert!(omega_dt < fs_phs::linear_exact::MIDPOINT_ONE_CENT_ANGULAR_STEP);
        assert!(GasValveMotion::new(gas(0.9), bare, cx).is_err());
        assert!(
            GasValveMotion::new(gas(0.9), motion(1000.0, 0.0, 0.0, AREA, RATE, cx), cx).is_ok()
        );
    });
}
