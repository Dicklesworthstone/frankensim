//! Independent PR energy balances through the original wave/wall owners.

use super::*;
use fs_gas::dynamics::nonideal::{
    RealGasWallContact, RealGasWallOperator, SharedWall, SharedWallControl, SharedWallGas,
    SharedWallPort,
};
use fs_gas::multicomponent::wall_exchange::{
    MaterialWallCaloric, SharedWallStorage, WallContact, WallControl, WallDrag, WallThermal,
};
use fs_thermal::enthalpy::{
    Boundary as ThermalBoundary, CellGeometry, CellId, DepositedCell, DepositionError, MaterialId,
    ThermalLaw, ThermalStorageBasis,
};

const R: f64 = 8.314_462_618_153_24;

fn wall_controls() -> WallControl {
    WallControl {
        absolute_energy_tolerance_j: 1e-7,
        relative_energy_tolerance: 1e-10,
        maximum_iterations: 128,
    }
}

fn contact(thermal: WallThermal) -> RealGasWallContact<'static> {
    RealGasWallContact {
        contact: WallContact {
            drag: WallDrag::None,
            thermal,
        },
        heat_storage: None,
    }
}

fn uniform(gas: RealGas, temperature: f64, velocity: f64) -> RealGasEuler1d {
    let q = state(&gas, 20.0, velocity, temperature);
    uniform_mesh_field(gas, vec![q; 3])
}

// Independent closed PR expression for this equal-mole binary fixture.
// sqrt(a_mix) = A - B*sqrt(T); the original species have mean Cp/R=4.
// No production caloric evaluator or nonlinear heat solve supplies expectations.
fn internal_energy(density: f64, temperature: f64) -> f64 {
    let mut a = 0.0;
    let mut b_temperature = 0.0;
    let mut covolume = 0.0;
    for (tc, pc, omega) in [(190.564, 4.5992e6, 0.01142), (305.322, 4.8722e6, 0.0995)] {
        let kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega;
        let root_a = (0.45724_f64 * R * R * tc * tc / pc).sqrt();
        a += 0.5 * root_a * (1.0 + kappa);
        b_temperature += 0.5 * root_a * kappa / tc.sqrt();
        covolume += 0.5 * 0.07780 * R * tc / pc;
    }
    let v = common::MIXTURE_MASS / density;
    let sqrt_two = 2.0_f64.sqrt();
    let logarithm = ((v + (1.0 + sqrt_two) * covolume) / (v + (1.0 - sqrt_two) * covolume)).ln();
    let residual =
        (a * b_temperature * temperature.sqrt() - a * a) * logarithm / (2.0 * sqrt_two * covolume);
    (3.0 * R * temperature + residual) / common::MIXTURE_MASS
}

fn root(mut low: f64, mut high: f64, equation: impl Fn(f64) -> f64) -> f64 {
    assert!(equation(low) <= 0.0 && equation(high) >= 0.0);
    for _ in 0..100 {
        let mid = 0.5 * low + 0.5 * high;
        if equation(mid) > 0.0 {
            high = mid;
        } else {
            low = mid;
        }
    }
    0.5 * low + 0.5 * high
}

fn check_wall_balance(flow: &RealGasEuler1d, walls: &RealGasWallOperator<'_>, report: StepReport) {
    let thermal = walls.last_report().unwrap();
    assert_eq!(thermal.elapsed, report.elapsed);
    assert_eq!(report.wall_heat, thermal.heat_to_gas_j);
    assert_eq!(report.wall_work, thermal.mechanical_work_to_gas_j);
    assert_eq!(report.wall_momentum_impulse, -thermal.wall_impulse_ns);
    near(report.conservation_residual.mass, 0.0, 1e-12);
    near(report.conservation_residual.momentum, 0.0, 2e-10);
    near(report.conservation_residual.energy, 0.0, 2e-7);
    near(thermal.source_storage_residual_j, 0.0, 2e-7);
    near(thermal.gas_wall_storage_residual_j, 0.0, 3e-7);
    near(
        report.after.energy - report.before.energy + thermal.wall_storage_change_j,
        report.boundary_change.energy
            + thermal.external_heat_to_gas_j
            + thermal.mechanical_work_to_gas_j,
        3e-7,
    );
    assert_eq!(flow.time(), report.elapsed);
}

#[test]
fn dense_reservoir_heating_and_cooling_use_original_variable_cv_and_accepted_duration() {
    for (initial, reservoir) in [(360.0, 850.0), (850.0, 360.0)] {
        let gas = gas();
        let initial_source = gas
            .thermodynamics(state(&gas, 20.0, 0.0, initial), || false)
            .unwrap();
        let mut flow = uniform(gas, initial, 0.0);
        let conductance = 1e7;
        let mut walls = RealGasWallOperator::new(
            &flow,
            vec![
                contact(WallThermal::Reservoir {
                    temperature_k: reservoir,
                    conductance_w_per_k: conductance,
                });
                3
            ],
            || false,
        )
        .unwrap();
        let report = flow
            .step_with_walls(
                Time::new(10.0),
                controls(),
                &mut walls,
                wall_controls(),
                || false,
            )
            .unwrap();
        assert!(report.elapsed.value() < 0.001);
        let mass = 20.0 * flow.mesh().volumes()[0];
        let integrated = conductance * report.elapsed.value();
        let expected = root(initial.min(reservoir), initial.max(reservoir), |t| {
            mass * (internal_energy(20.0, t) - internal_energy(20.0, initial))
                - integrated * (reservoir - t)
        });
        let cv = initial_source.specific_heat_cv_j_kg_k();
        let frozen = (mass * cv * initial + integrated * reservoir) / (mass * cv + integrated);
        assert!(
            (frozen - expected).abs() > 0.01,
            "fixture must reject frozen Cv"
        );
        let cp_minus_r = initial_source.specific_heat_cp_j_kg_k() - R / common::MIXTURE_MASS;
        assert!((cp_minus_r - cv).abs() > 10.0, "dense Cv is not Cp-R");
        for (q, primitive) in flow.states().iter().zip(flow.primitives()) {
            near(primitive.temperature, expected, 3e-8);
            assert_eq!(q.density, 20.0);
            let source = flow.gas().thermodynamics(*q, || false).unwrap();
            assert!((source.specific_heat_cv_j_kg_k() - cv).abs() > 0.1);
            near(primitive.pressure, source.pressure_pa(), 1e-6);
            near(primitive.sound_speed, source.sound_speed_m_s(), 1e-9);
        }
        let expected_heat =
            3.0 * mass * (internal_energy(20.0, expected) - internal_energy(20.0, initial));
        near(report.wall_heat, expected_heat, 2e-5);
        assert_eq!(
            walls.last_report().unwrap().external_heat_to_gas_j,
            report.wall_heat
        );
        check_wall_balance(&flow, &walls, report);
    }
}

#[test]
fn caloric_wall_trials_cross_transport_limits_but_accepted_gas_must_remain_supported() {
    let mut wider = domain();
    wider.energy.minimum_temperature_k = 250.0;
    let gas = bind(common::model([0.0; 2], [0.0; 2], [200.0, 2000.0]), wider);
    // The original caloric branch admits the cold reservoir. Its dilute-gas
    // transport cards start at 300 K, but no gas endpoint need reach 250 K.
    assert!(gas.model().evaluate_density_temperature(&point(250.0), 20.0).is_err());
    assert!(gas.model().specific_internal_energy_change_and_cv_between_densities(
        20.0, 350.0, 20.0, 250.0,
    ).is_ok());
    let mut flow = uniform(gas, 350.0, 0.0);
    let conductance = 1e6;
    let make_wall = |conductance| contact(WallThermal::Reservoir {
        temperature_k: 250.0, conductance_w_per_k: conductance,
    });
    let mut walls = RealGasWallOperator::new(
        &flow, vec![make_wall(conductance); 3], || false,
    ).unwrap();
    let mut no_retry = controls();
    no_retry.maximum_halvings = 0;
    let report = flow.step_with_walls(
        Time::new(1e-6), no_retry, &mut walls, wall_controls(), || false,
    ).unwrap();
    assert_eq!(report.elapsed, Time::new(1e-6));
    let mass = 20.0 * flow.mesh().volumes()[0];
    let expected = root(250.0, 350.0, |t| {
        mass * (internal_energy(20.0, t) - internal_energy(20.0, 350.0))
            - conductance * report.elapsed.value() * (250.0 - t)
    });
    assert!(expected > 300.0 && expected < 350.0);
    for primitive in flow.primitives() {
        near(primitive.temperature, expected, 3e-8);
    }
    check_wall_balance(&flow, &walls, report);

    // Stronger heat exchange would finish below the transport domain. The
    // final full decode must still refuse, without publishing any owner.
    let mut outside = RealGasWallOperator::new(
        &flow, vec![make_wall(1e12); 3], || false,
    ).unwrap();
    let snapshot = (flow.states().to_vec(), flow.primitives().to_vec(), flow.time(),
        outside.contacts().to_vec());
    assert!(matches!(flow.step_with_walls(
        Time::new(1e-6), no_retry, &mut outside, wall_controls(), || false,
    ), Err(RealGasFlowError::StepRetries { .. })));
    assert_eq!((flow.states().to_vec(), flow.primitives().to_vec(), flow.time(),
        outside.contacts().to_vec()), snapshot);
    assert_eq!(outside.last_report(), None);
}

// Original enthalpy chart with Cp(T)=T. This is deliberately nonlinear even
// though its T(h) remains elementary; a cached capacity times delta T is wrong.
struct QuadraticEnthalpy;
impl ThermalLaw for QuadraticEnthalpy {
    fn storage_basis(&self) -> ThermalStorageBasis {
        ThermalStorageBasis::Enthalpy
    }
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !(45_000.0..=500_000.0).contains(&h) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok((2.0 * h).sqrt())
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / 300.0
    }
}

fn wall_cell(energy: f64) -> DepositedCell {
    DepositedCell {
        id: CellId(17),
        segment_id: 3,
        body_id: 8,
        material_id: MaterialId(29),
        mass_kg: 0.75,
        specific_h_j_kg: energy,
        geometry: CellGeometry {
            start_m: [0.0; 3],
            end_m: [1.0, 0.0, 0.0],
            volume_m3: 0.001,
            surface_area_m2: 1.0,
        },
        boundary: ThermalBoundary {
            ambient_k: 300.0,
            htc_w_m2_k: 0.0,
            reservoir: None,
            source_w: 0.0,
        },
    }
}

fn finite_wall(temperature: f64, conductance: f64) -> RealGasWallContact<'static> {
    let source = MaterialWallCaloric::from_thermal_law(
        QuadraticEnthalpy,
        wall_cell(0.5 * temperature * temperature),
        [45_000.0, 500_000.0],
    )
    .unwrap();
    RealGasWallContact {
        contact: source.initial_contact(conductance, WallDrag::None).unwrap(),
        heat_storage: Some(source),
    }
}

fn specific_wall_energy(walls: &RealGasWallOperator<'_>, cell: usize) -> f64 {
    let WallThermal::Caloric {
        specific_energy_j_kg,
        ..
    } = walls.contacts()[cell].thermal
    else {
        panic!("original finite enthalpy inventory required");
    };
    specific_energy_j_kg
}

fn shared_wall(input: RealGasWallContact<'static>, conductances: &[f64]) -> SharedWall<'static> {
    SharedWall {
        storage: SharedWallStorage {
            thermal: input.contact.thermal,
            heat_storage: input.heat_storage,
        },
        ports: conductances.iter().enumerate().map(|(cell, &conductance_w_per_k)| SharedWallPort {
            gas: SharedWallGas::Branch { branch: 0, cell }, conductance_w_per_k,
        }).collect(),
    }
}

#[test]
fn shared_nonlinear_wall_uses_each_current_real_gas_density_and_only_one_heat_inventory() {
    let gas = gas();
    let initial = [(12.0, 30.0, 350.0), (20.0, -10.0, 650.0), (30.0, 15.0, 800.0)]
        .into_iter().map(|(rho, velocity, t)| state(&gas, rho, velocity, t)).collect::<Vec<_>>();
    let mesh = Mesh1d::new(vec![0.0, 0.2, 0.5, 1.0], vec![0.07; 4],
        vec![0.014, 0.021, 0.035]).unwrap();
    let mut wave = RealGasEuler1d::new_on_mesh(gas.clone(), mesh.clone(), initial.clone(),
        RealGasBoundaries::Periodic, || false).unwrap();
    let mut flow = RealGasEuler1d::new_on_mesh(gas, mesh, initial,
        RealGasBoundaries::Periodic, || false).unwrap();
    let conductances = [1e8, 2e8, 7e7];
    let wall_initial = 500.0;
    let mut walls = RealGasWallOperator::new(&flow,
        vec![contact(WallThermal::Adiabatic); 3], || false).unwrap()
        .with_shared_walls(vec![shared_wall(finite_wall(wall_initial, 0.0), &conductances)],
            SharedWallControl::default(), || false).unwrap();
    let requested = Time::new(1e-6);
    let transport = wave.step(requested, controls(), || false).unwrap();
    let heat_controls = WallControl { relative_energy_tolerance: 0.0, ..wall_controls() };
    let report = flow.step_with_walls(requested, controls(), &mut walls,
        heat_controls, || false).unwrap();
    assert_eq!(report.elapsed, requested);
    assert_eq!(report.elapsed, transport.elapsed);
    assert_eq!(report.boundary_change, transport.boundary_change);
    // Independent nested scalar balance using the closed PR expression above.
    // The wave's actual candidate densities, momenta and original cell volumes
    // enter the heat solve; no reference-density or frozen-Cv approximation does.
    let gas_endpoint = |wall_t: f64, i: usize| {
        let rho = wave.states()[i].density;
        let t0 = wave.primitives()[i].temperature;
        let mass = rho * wave.mesh().volumes()[i];
        root(t0.min(wall_t), t0.max(wall_t), |t| {
            mass * (internal_energy(rho, t) - internal_energy(rho, t0))
                - report.elapsed.value() * conductances[i] * (wall_t - t)
        })
    };
    let gas_heat = |wall_t: f64| (0..3).map(|i| {
        report.elapsed.value() * conductances[i] * (wall_t - gas_endpoint(wall_t, i))
    }).sum::<f64>();
    let lower = wave.primitives().iter().fold(wall_initial, |t, p| t.min(p.temperature));
    let upper = wave.primitives().iter().fold(wall_initial, |t, p| t.max(p.temperature));
    let wall_change = |t: f64| 0.75 * 0.5 * (t - wall_initial) * (t + wall_initial);
    let expected_wall = root(lower, upper, |t| wall_change(t) + gas_heat(t));
    let duplicated_wall = root(lower, upper, |t| 3.0 * wall_change(t) + gas_heat(t));
    assert!((expected_wall - duplicated_wall).abs() > 1.0,
        "fixture must distinguish one finite wall from three copies");
    let shared = walls.shared_wall(0).unwrap();
    let source = shared.storage.heat_storage.as_ref().unwrap();
    near(source.temperature_k(shared.storage.thermal).unwrap(), expected_wall, 5e-8);
    let mut heating = false;
    let mut cooling = false;
    for i in 0..3 {
        assert_eq!(flow.states()[i].density, wave.states()[i].density);
        assert_eq!(flow.states()[i].momentum, wave.states()[i].momentum);
        near(flow.primitives()[i].temperature, gas_endpoint(expected_wall, i), 5e-8);
        heating |= flow.primitives()[i].temperature > wave.primitives()[i].temperature;
        cooling |= flow.primitives()[i].temperature < wave.primitives()[i].temperature;
    }
    assert!(heating && cooling);
    near(report.wall_heat, gas_heat(expected_wall), 5e-5);
    near(walls.last_report().unwrap().wall_storage_change_j, wall_change(expected_wall), 5e-5);
    let shared_report = walls.last_shared_wall_report().unwrap();
    assert!(shared_report.source_evaluations > 0);
    assert!(shared_report.contact_solves >= 3);
    assert_eq!(shared_report.elapsed, report.elapsed);
    near(shared_report.energy_residual_j, 0.0, 3e-5);
    check_wall_balance(&flow, &walls, report);
}

#[test]
fn finite_source_cp_retains_nonlinear_enthalpy_and_cancels_formation_datums() {
    for (initial, wall_initial) in [(360.0, 800.0), (800.0, 360.0)] {
        let conductance = 1e7;
        let mut accepted = None;
        for datums in [[0.0; 2], [1e22, -3e21]] {
            let gas = bind(
                common::model(datums, [2e18, -7e17], [200.0, 2000.0]),
                domain(),
            );
            let mut flow = uniform(gas, initial, 0.0);
            let mut walls = RealGasWallOperator::new(
                &flow,
                vec![finite_wall(wall_initial, conductance); 3],
                || false,
            )
            .unwrap();
            let report = flow
                .step_with_walls(
                    Time::new(10.0),
                    controls(),
                    &mut walls,
                    wall_controls(),
                    || false,
                )
                .unwrap();
            let mass = 20.0 * flow.mesh().volumes()[0];
            let integrated = conductance * report.elapsed.value();
            // Eliminate gas temperature in the independent wall-temperature
            // coordinate; its original quadratic enthalpy supplies every joule.
            let expected_wall = root(
                initial.min(wall_initial),
                initial.max(wall_initial),
                |wall_t| {
                    let heat = 0.75 * 0.5 * (wall_initial * wall_initial - wall_t * wall_t);
                    let gas_t = wall_t - heat / integrated;
                    if gas_t < initial.min(wall_initial) {
                        return -f64::MAX;
                    }
                    if gas_t > initial.max(wall_initial) {
                        return f64::MAX;
                    }
                    mass * (internal_energy(20.0, gas_t) - internal_energy(20.0, initial)) - heat
                },
            );
            let heat = 0.75 * 0.5 * (wall_initial * wall_initial - expected_wall * expected_wall);
            let expected_gas = expected_wall - heat / integrated;
            for (i, p) in flow.primitives().iter().enumerate() {
                near(p.temperature, expected_gas, 5e-8);
                near(
                    specific_wall_energy(&walls, i),
                    0.5 * expected_wall * expected_wall,
                    3e-5,
                );
                assert_eq!(
                    walls.heat_storage()[i].as_ref().unwrap().storage_basis(),
                    ThermalStorageBasis::Enthalpy
                );
            }
            near(report.wall_heat, 3.0 * heat, 3e-5);
            assert_eq!(walls.last_report().unwrap().external_heat_to_gas_j, 0.0);
            check_wall_balance(&flow, &walls, report);
            let snapshot = (flow.states().to_vec(), walls.contacts().to_vec(), report);
            if let Some(original) = &accepted {
                assert_eq!(&snapshot, original);
            }
            accepted = Some(snapshot);
        }
    }
}

struct LatentEnthalpy;
impl ThermalLaw for LatentEnthalpy {
    fn storage_basis(&self) -> ThermalStorageBasis {
        ThermalStorageBasis::Enthalpy
    }
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !(0.0..=200_000.0).contains(&h) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(if h < 50_000.0 {
            450.0 + h / 500.0
        } else if h <= 150_000.0 {
            550.0
        } else {
            550.0 + (h - 150_000.0) / 500.0
        })
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / 500.0
    }
}

#[test]
fn shared_latent_wall_keeps_one_phase_inventory_while_every_dense_gas_cell_heats() {
    let mut flow = uniform(gas(), 530.0, 0.0);
    let source = MaterialWallCaloric::from_thermal_law(
        LatentEnthalpy, wall_cell(100_000.0), [0.0, 200_000.0],
    ).unwrap();
    let input = RealGasWallContact {
        contact: source.initial_contact(0.0, WallDrag::None).unwrap(),
        heat_storage: Some(source),
    };
    let conductance = 1e7;
    let mut walls = RealGasWallOperator::new(&flow,
        vec![contact(WallThermal::Adiabatic); 3], || false).unwrap()
        .with_shared_walls(vec![shared_wall(input, &[conductance; 3])],
            SharedWallControl::default(), || false).unwrap();
    let heat_controls = WallControl { relative_energy_tolerance: 0.0, ..wall_controls() };
    let report = flow.step_with_walls(Time::new(1e-5), controls(), &mut walls,
        heat_controls, || false).unwrap();
    let mass = 20.0 * flow.mesh().volumes()[0];
    let expected = root(530.0, 550.0, |t| {
        mass * (internal_energy(20.0, t) - internal_energy(20.0, 530.0))
            - conductance * report.elapsed.value() * (550.0 - t)
    });
    let heat = 3.0 * mass * (internal_energy(20.0, expected) - internal_energy(20.0, 530.0));
    let shared = walls.shared_wall(0).unwrap();
    let source = shared.storage.heat_storage.as_ref().unwrap();
    assert_eq!(source.temperature_k(shared.storage.thermal).unwrap(), 550.0);
    let WallThermal::Caloric { specific_energy_j_kg, .. } = shared.storage.thermal else {
        panic!("one original latent inventory required");
    };
    near(specific_energy_j_kg, 100_000.0 - heat / 0.75, 2e-5);
    assert!(specific_energy_j_kg > 50_000.0 && specific_energy_j_kg < 100_000.0);
    for primitive in flow.primitives() {
        near(primitive.temperature, expected, 3e-8);
    }
    near(report.wall_heat, heat, 2e-5);
    check_wall_balance(&flow, &walls, report);
}

#[test]
fn shared_wall_binding_work_and_late_cancellation_preserve_every_accepted_owner() {
    let input = || shared_wall(contact(WallThermal::Lumped {
        temperature_k: 700.0, heat_capacity_j_per_k: 400.0, conductance_w_per_k: 0.0,
    }), &[1e7; 3]);
    let make = || {
        let flow = uniform(gas(), 400.0, 0.0);
        let walls = RealGasWallOperator::new(&flow,
            vec![contact(WallThermal::Adiabatic); 3], || false).unwrap()
            .with_shared_walls(vec![input()], SharedWallControl::default(), || false).unwrap();
        (flow, walls)
    };
    let mut no_retry = controls();
    no_retry.maximum_halvings = 0;
    let heat_controls = WallControl { relative_energy_tolerance: 0.0, ..wall_controls() };
    let (mut good, mut good_walls) = make();
    let (mut refused, mut refused_walls) = make();
    for (flow, walls) in [(&mut good, &mut good_walls), (&mut refused, &mut refused_walls)] {
        let report = flow.step_with_walls(Time::new(1e-6), no_retry, walls,
            heat_controls, || false).unwrap();
        let WallThermal::Lumped { temperature_k, .. } = walls.shared_wall(0).unwrap().storage.thermal else {
            panic!("one authored finite wall expected");
        };
        near(walls.last_report().unwrap().wall_storage_change_j, 400.0 * (temperature_k - 700.0), 1e-8);
        check_wall_balance(flow, walls, report);
    }
    let snapshot = (
        refused.states().to_vec(), refused.primitives().to_vec(), refused.time(),
        refused_walls.contacts().to_vec(), refused_walls.shared_wall(0).unwrap().clone(),
        refused_walls.last_report(), refused_walls.last_shared_wall_report(),
    );
    let mut polls = 0;
    good.step_with_walls(Time::new(1e-6), no_retry, &mut good_walls, heat_controls, || {
        polls += 1; false
    }).unwrap();
    let mut seen = 0;
    assert_eq!(refused.step_with_walls(Time::new(1e-6), no_retry, &mut refused_walls,
        heat_controls, || { seen += 1; seen == polls }), Err(RealGasFlowError::Cancelled));
    assert_eq!((
        refused.states().to_vec(), refused.primitives().to_vec(), refused.time(),
        refused_walls.contacts().to_vec(), refused_walls.shared_wall(0).unwrap().clone(),
        refused_walls.last_report(), refused_walls.last_shared_wall_report(),
    ), snapshot);

    let mut limited = RealGasWallOperator::new(&refused,
        vec![contact(WallThermal::Adiabatic); 3], || false).unwrap()
        .with_shared_walls(vec![input()], SharedWallControl { maximum_source_evaluations: 1 },
            || false).unwrap();
    let initial_wall = limited.shared_wall(0).unwrap().clone();
    assert!(matches!(refused.step_with_walls(Time::new(1e-6), no_retry, &mut limited,
        heat_controls, || false), Err(RealGasFlowError::StepRetries { .. })));
    assert_eq!(refused.states(), snapshot.0);
    assert_eq!(refused.primitives(), snapshot.1);
    assert_eq!(refused.time(), snapshot.2);
    assert_eq!(limited.shared_wall(0).unwrap(), &initial_wall);
    assert_eq!(limited.last_report(), None);
    assert_eq!(limited.last_shared_wall_report(), None);

    for location in [SharedWallGas::Plenum(0), SharedWallGas::Branch { branch: 1, cell: 0 },
        SharedWallGas::Branch { branch: 0, cell: 3 }, SharedWallGas::Branch { branch: 0, cell: 1 }]
    {
        let mut invalid = input();
        invalid.ports[0].gas = location;
        let owner = RealGasWallOperator::new(&refused,
            vec![contact(WallThermal::Adiabatic); 3], || false).unwrap();
        assert!(owner.with_shared_walls(vec![invalid], SharedWallControl::default(), || false).is_err());
    }
    let owner = RealGasWallOperator::new(&refused, vec![contact(WallThermal::Reservoir {
        temperature_k: 700.0, conductance_w_per_k: 1.0,
    }); 3], || false).unwrap();
    assert!(owner.with_shared_walls(vec![input()], SharedWallControl::default(), || false).is_err());
}

#[test]
fn latent_inventory_changes_while_wall_temperature_stays_on_its_original_plateau() {
    let mut flow = uniform(gas(), 530.0, 0.0);
    let source = MaterialWallCaloric::from_thermal_law(
        LatentEnthalpy,
        wall_cell(100_000.0),
        [0.0, 200_000.0],
    )
    .unwrap();
    let conductance = 1e7;
    let input = RealGasWallContact {
        contact: source.initial_contact(conductance, WallDrag::None).unwrap(),
        heat_storage: Some(source),
    };
    let mut walls = RealGasWallOperator::new(&flow, vec![input; 3], || false).unwrap();
    let report = flow
        .step_with_walls(
            Time::new(1.0),
            controls(),
            &mut walls,
            wall_controls(),
            || false,
        )
        .unwrap();
    let integrated = conductance * report.elapsed.value();
    let mass = 20.0 * flow.mesh().volumes()[0];
    let expected = root(530.0, 550.0, |t| {
        mass * (internal_energy(20.0, t) - internal_energy(20.0, 530.0)) - integrated * (550.0 - t)
    });
    let heat = mass * (internal_energy(20.0, expected) - internal_energy(20.0, 530.0));
    for (i, p) in flow.primitives().iter().enumerate() {
        near(p.temperature, expected, 3e-8);
        assert_eq!(
            walls.heat_storage()[i]
                .as_ref()
                .unwrap()
                .temperature_k(walls.contacts()[i].thermal)
                .unwrap(),
            550.0
        );
        near(
            specific_wall_energy(&walls, i),
            100_000.0 - heat / 0.75,
            2e-5,
        );
        assert!(specific_wall_energy(&walls, i) < 100_000.0);
    }
    check_wall_balance(&flow, &walls, report);
}

#[test]
fn sliding_drag_retains_original_real_gas_internal_heating_and_counts_drive_work_once() {
    let mut flow = uniform(gas(), 430.0, 90.0);
    let mut input = contact(WallThermal::Adiabatic);
    input.contact.drag = WallDrag::MovingLinear {
        rate_per_s: 2000.0,
        wall_velocity_m_s: 20.0,
    };
    let mut walls = RealGasWallOperator::new(&flow, vec![input; 3], || false).unwrap();
    let report = flow
        .step_with_walls(
            Time::new(1.0),
            controls(),
            &mut walls,
            wall_controls(),
            || false,
        )
        .unwrap();
    let velocity = 20.0 + 70.0 * (-2000.0 * report.elapsed.value()).exp();
    let specific_work = 20.0 * (velocity - 90.0);
    let internal_change = specific_work - 0.5 * (velocity * velocity - 90.0 * 90.0);
    let expected = root(430.0, 500.0, |t| {
        internal_energy(20.0, t) - internal_energy(20.0, 430.0) - internal_change
    });
    for p in flow.primitives() {
        near(p.velocity, velocity, 1e-12);
        near(p.temperature, expected, 3e-8);
    }
    near(report.wall_work, report.before.mass * specific_work, 2e-9);
    near(
        walls.last_report().unwrap().drag_heating_j,
        0.5 * report.before.mass * (70.0 * 70.0 - (velocity - 20.0).powi(2)),
        2e-9,
    );
    assert_eq!(report.wall_heat, 0.0);
    check_wall_balance(&flow, &walls, report);
}

#[test]
fn exact_thermal_equilibrium_retains_the_original_wave_step_even_at_stiff_conductance() {
    let mut original = uniform(gas(), 400.0, 0.0);
    let mut coupled = uniform(gas(), 400.0, 0.0);
    let mut walls = RealGasWallOperator::new(
        &coupled,
        vec![
            finite_wall(400.0, 1e20),
            contact(WallThermal::Adiabatic),
            contact(WallThermal::Reservoir {
                temperature_k: 400.0,
                conductance_w_per_k: 1e20,
            }),
        ],
        || false,
    )
    .unwrap();
    let contacts = walls.contacts().to_vec();
    let expected = original.step(Time::new(1.0), controls(), || false).unwrap();
    let report = coupled
        .step_with_walls(
            Time::new(1.0),
            controls(),
            &mut walls,
            wall_controls(),
            || false,
        )
        .unwrap();
    assert_eq!(coupled.states(), original.states());
    assert_eq!(coupled.primitives(), original.primitives());
    assert_eq!(walls.contacts(), contacts);
    assert_eq!(report, expected);
    assert_eq!(report.wall_heat, 0.0);
    check_wall_balance(&coupled, &walls, report);
}

#[test]
fn compression_then_heat_uses_each_candidate_density_and_keeps_boundary_transport_separate() {
    let gas = gas();
    let initial: Vec<_> = [
        (16.0, 30.0, 420.0),
        (20.0, -15.0, 460.0),
        (24.0, 10.0, 520.0),
    ]
    .into_iter()
    .map(|(rho, u, t)| state(&gas, rho, u, t))
    .collect();
    let mut transport = uniform_mesh_field(gas.clone(), initial.clone());
    let mut coupled = uniform_mesh_field(gas, initial.clone());
    let conductance = 1e7;
    let mut walls = RealGasWallOperator::new(
        &coupled,
        vec![
            contact(WallThermal::Reservoir {
                temperature_k: 700.0,
                conductance_w_per_k: conductance,
            }),
            finite_wall(360.0, conductance),
            contact(WallThermal::Adiabatic),
        ],
        || false,
    )
    .unwrap();
    let requested = Time::new(1e-5);
    let wave = transport.step(requested, controls(), || false).unwrap();
    let report = coupled
        .step_with_walls(requested, controls(), &mut walls, wall_controls(), || false)
        .unwrap();
    assert_eq!(report.elapsed, wave.elapsed);
    assert_eq!(report.boundary_change, wave.boundary_change);
    assert_eq!(report.left_boundary_change, wave.left_boundary_change);
    assert_eq!(report.right_boundary_change, wave.right_boundary_change);
    let mut heat = 0.0;
    for i in 0..3 {
        let before = transport.states()[i];
        let after = coupled.states()[i];
        assert_eq!(after.density, before.density);
        assert_eq!(after.momentum, before.momentum);
        assert!((before.density - initial[i].density).abs() > 1e-4);
        let volume = coupled.mesh().volumes()[i];
        let actual = (after.energy - before.energy) * volume;
        let expected = before.density
            * volume
            * (internal_energy(before.density, coupled.primitives()[i].temperature)
                - internal_energy(before.density, transport.primitives()[i].temperature));
        near(actual, expected, 3e-7);
        if i == 0 {
            near(
                actual,
                conductance
                    * report.elapsed.value()
                    * (700.0 - coupled.primitives()[i].temperature),
                3e-6,
            );
        } else if i == 1 {
            let wall_temperature = walls.heat_storage()[i]
                .as_ref()
                .unwrap()
                .temperature_k(walls.contacts()[i].thermal)
                .unwrap();
            near(
                actual,
                conductance
                    * report.elapsed.value()
                    * (wall_temperature - coupled.primitives()[i].temperature),
                3e-6,
            );
            near(
                actual + 0.75 * (specific_wall_energy(&walls, i) - 0.5 * 360.0 * 360.0),
                0.0,
                1e-7,
            );
        } else {
            assert_eq!(after, before);
        }
        heat += actual;
    }
    near(report.wall_heat, heat, 1e-7);
    check_wall_balance(&coupled, &walls, report);
}

#[test]
fn source_domain_binding_and_late_cancellation_leave_gas_walls_reports_and_clock_unchanged() {
    let make = || {
        let flow = uniform(gas(), 400.0, 0.0);
        let walls =
            RealGasWallOperator::new(&flow, vec![finite_wall(700.0, 1e7); 3], || false).unwrap();
        (flow, walls)
    };
    let mut no_retry = controls();
    no_retry.maximum_halvings = 0;
    let (mut good, mut good_walls) = make();
    let (mut refused, mut refused_walls) = make();
    for (flow, walls) in [
        (&mut good, &mut good_walls),
        (&mut refused, &mut refused_walls),
    ] {
        flow.step_with_walls(Time::new(1.0), no_retry, walls, wall_controls(), || false)
            .unwrap();
    }
    let snapshot = (
        refused.states().to_vec(),
        refused.primitives().to_vec(),
        refused.time(),
        refused_walls.contacts().to_vec(),
        refused_walls.last_report(),
    );
    let mut polls = 0;
    good.step_with_walls(
        Time::new(1.0),
        no_retry,
        &mut good_walls,
        wall_controls(),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let mut seen = 0;
    assert_eq!(
        refused.step_with_walls(
            Time::new(1.0),
            no_retry,
            &mut refused_walls,
            wall_controls(),
            || {
                seen += 1;
                seen == polls
            }
        ),
        Err(RealGasFlowError::Cancelled)
    );
    assert_eq!(
        (
            refused.states().to_vec(),
            refused.primitives().to_vec(),
            refused.time(),
            refused_walls.contacts().to_vec(),
            refused_walls.last_report()
        ),
        snapshot
    );

    let other = bind(
        common::model([1e20, -2e20], [0.0; 2], [200.0, 2000.0]),
        domain(),
    );
    let other_flow = uniform(other, 400.0, 0.0);
    let mut mismatched =
        RealGasWallOperator::new(&other_flow, vec![finite_wall(700.0, 1e7); 3], || false).unwrap();
    assert!(matches!(
        refused.step_with_walls(
            Time::new(1.0),
            no_retry,
            &mut mismatched,
            wall_controls(),
            || false
        ),
        Err(RealGasFlowError::Invalid(_))
    ));
    assert_eq!(refused.states(), snapshot.0);
    assert_eq!(refused.time(), snapshot.2);
    assert_eq!(mismatched.last_report(), None);

    let mut outside = RealGasWallOperator::new(
        &refused,
        vec![
            contact(WallThermal::Reservoir {
                temperature_k: 10_000.0,
                conductance_w_per_k: 1e12,
            });
            3
        ],
        || false,
    )
    .unwrap();
    let initial_contacts = outside.contacts().to_vec();
    assert!(matches!(
        refused.step_with_walls(
            Time::new(1.0),
            no_retry,
            &mut outside,
            wall_controls(),
            || false
        ),
        Err(RealGasFlowError::StepRetries { .. })
    ));
    assert_eq!(refused.states(), snapshot.0);
    assert_eq!(refused.primitives(), snapshot.1);
    assert_eq!(refused.time(), snapshot.2);
    assert_eq!(outside.contacts(), initial_contacts);
    assert_eq!(outside.last_report(), None);
    assert!(
        RealGasWallOperator::new(
            &refused,
            vec![
                contact(WallThermal::Caloric {
                    mass_kg: 1.0,
                    specific_energy_j_kg: 1e5,
                    conductance_w_per_k: 1.0,
                });
                3
            ],
            || false
        )
        .is_err()
    );
}
