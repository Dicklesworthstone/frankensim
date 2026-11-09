//! Native finite caloric feedback, original memory heat and global source phase.
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};
use fs_mbd::rotordynamics::RotorBearingModel;
use fs_mbd::rotordynamics::periodic::{
    RotorFloquetBudget, RotorHarmonicForce, RotorIntervalForcing, RotorMaxwellState,
};
use fs_mbd::rotordynamics::relaxation::MaxwellDampingBranch;
use fs_mbd_material::rotor_thermal::{
    RotorThermalConfig, RotorThermalError, RotorThermalInput, advance_candidate,
};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, ThermalLink,
    ThermalStorageBasis, TimedEvent,
};

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "{actual:e} != {expected:e}"
    );
}

// Explicit synthetic sensible capacity; no Cp/Cv or source accuracy is inferred.
struct Capacity(f64);
impl ThermalLaw for Capacity {
    fn storage_basis(&self) -> ThermalStorageBasis {
        ThermalStorageBasis::DeclaredSensibleStorage
    }
    fn temperature_k(&self, e: f64) -> Result<f64, DepositionError> {
        let temperature = 300.0 + e / self.0;
        if (250.0..=800.0).contains(&temperature) {
            Ok(temperature)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        self.0.recip()
    }
}

fn boundary(source_w: f64) -> Boundary {
    Boundary {
        ambient_k: 300.0,
        htc_w_m2_k: 0.0,
        reservoir: None,
        source_w,
    }
}

fn timeline<'a>(
    law: &'a dyn ThermalLaw,
    cells: &[(u64, f64, f64)],
    extra: Vec<TimedEvent>,
) -> DepositionSimulation<'a> {
    let mut events: Vec<_> = cells
        .iter()
        .map(|&(id, mass_kg, specific_h_j_kg)| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(id),
                segment_id: id,
                body_id: 1,
                material_id: MaterialId(1),
                mass_kg,
                specific_h_j_kg,
                geometry: CellGeometry {
                    start_m: [0.0, 0.0, id as f64],
                    end_m: [1.0, 0.0, id as f64],
                    volume_m3: 0.001,
                    surface_area_m2: 1.0,
                },
                boundary: boundary(0.0),
            }),
        })
        .collect();
    events.extend(extra);
    DepositionSimulation::new(
        events,
        &[(MaterialId(1), law)],
        DepositionBudget {
            max_cells: 4,
            max_links: 4,
            max_events: 16,
            max_substeps: 4096,
        },
    )
    .unwrap()
}

fn shift(validity: [f64; 2]) -> TemperatureShift {
    TemperatureShift::new(
        TemperatureShiftModel::WilliamsLandelFerry {
            c1: 4.0,
            c2_k: 100.0,
        },
        300.0,
        validity,
    )
    .unwrap()
}
fn branch(k: f64, tau: f64) -> MaxwellDampingBranch {
    MaxwellDampingBranch {
        strain: vec![0.0],
        stiffness: k,
        relaxation_time_s: tau,
    }
}
fn config(dt: f64) -> RotorThermalConfig {
    RotorThermalConfig {
        rotor: RotorFloquetBudget::default(),
        relative_generator_residual_tolerance: 1e-11,
        relative_map_residual_tolerance: 1e-9,
        thermal: StepConfig {
            max_dt_s: dt,
            energy_tolerance_j: 1e-8,
        },
        maximum_intervals: 2048,
        energy_tolerance_j: 1e-8,
    }
}

fn isolated_reference(time: f64) -> f64 {
    // Independent scalar ODE: sigma'=-sigma/tau(T), with the exact conserved
    // relation T=300+(100-sigma²)/4 for k=2 and total heat capacity 1 J/K.
    // RK4 integrates this nonlinear physical feedback, not the adapter's split.
    let rhs = |stress: f64| {
        let delta = (100.0 - stress * stress) / 4.0;
        let ratio = 10.0_f64.powf(-4.0 * delta / (100.0 + delta));
        -stress / ratio
    };
    let dt = time / 16384.0;
    let mut stress = 10.0;
    for _ in 0..16384 {
        let a = rhs(stress);
        let b = rhs(stress + 0.5 * dt * a);
        let c = rhs(stress + 0.5 * dt * b);
        let d = rhs(stress + dt * c);
        stress += dt * (a + 2.0 * b + 2.0 * c + d) / 6.0;
    }
    stress
}

#[test]
fn self_heating_feedback_converges_to_independent_nonlinear_memory_ode() {
    let law = Capacity(1.0);
    let thermal = timeline(&law, &[(1, 1.0, 0.0)], vec![]);
    let model = RotorBearingModel::new(1, &[1.0], &[0.2], &[1.0], &[0.0], || false).unwrap();
    let branches = [branch(2.0, 1.0)];
    let shifts = [shift([290.0, 400.0])];
    let cells = [CellId(1)];
    let state = RotorMaxwellState {
        displacement: vec![0.0],
        velocity: vec![0.0],
        branch_stress: vec![10.0],
    };
    let input = RotorThermalInput {
        time_s: 0.0,
        state: &state,
        model: &model,
        spin_rad_s: 0.0,
        branches: &branches,
        temperature_shifts: &shifts,
        thermal_cells: &cells,
        forcing: RotorIntervalForcing::Homogeneous,
    };
    let target = 0.4;
    let exact_stress = isolated_reference(target);
    let mut errors = Vec::new();
    for dt in [0.05, 0.0125, 0.003125] {
        let next = advance_candidate(input, &thermal, target, config(dt), || false).unwrap();
        let stress = next.state().branch_stress[0];
        let heat = (100.0 - stress * stress) / 4.0;
        near(
            next.thermal().temperature_k(CellId(1)).unwrap(),
            300.0 + heat,
            2e-11,
        );
        near(next.report().branch_heat_j[0], heat, 2e-10);
        near(next.report().thermal.stored_energy_change_j, heat, 2e-10);
        near(next.report().external_thermal_heat_j, 0.0, 1e-10);
        near(next.report().combined_energy_residual_j, 0.0, 1e-8);
        // Heating speeds the original clock; it cannot equal a frozen 300 K decay.
        assert!(stress < 10.0 * (-target).exp());
        assert_eq!(next.state().displacement, [0.0]);
        assert_eq!(next.state().velocity, [0.0]);
        errors.push((stress - exact_stress).abs());
    }
    assert!(errors[1] < 0.4 * errors[0], "{errors:?}");
    assert!(errors[2] < 0.4 * errors[1], "{errors:?}");
    assert!(errors[2] < 0.02 * exact_stress.abs(), "{errors:?}");
    assert!(thermal.cells().is_empty()); // Original t=0 births were not consumed.
    assert_eq!(state.branch_stress, [10.0]);
}

#[test]
fn original_branch_heat_routes_through_shared_cells_and_native_conduction() {
    let law = Capacity(2.0);
    let thermal = timeline(
        &law,
        &[(1, 2.0, 0.0), (2, 1.0, 60.0)],
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 3,
                a: CellId(1),
                b: CellId(2),
                area_m2: 0.25,
                conductance_w_k: 0.1,
            }),
        }],
    );
    let model = RotorBearingModel::new(1, &[1.0], &[0.2], &[1.0], &[0.0], || false).unwrap();
    let branches = [branch(2.0, 0.5), branch(4.0, 0.7), branch(3.0, 0.9)];
    let shifts = [shift([290.0, 400.0]); 3];
    let cells = [CellId(1), CellId(1), CellId(2)];
    let state = RotorMaxwellState {
        displacement: vec![0.0],
        velocity: vec![0.0],
        branch_stress: vec![3.0, -4.0, 2.0],
    };
    let input = RotorThermalInput {
        time_s: 0.0,
        state: &state,
        model: &model,
        spin_rad_s: 0.0,
        branches: &branches,
        temperature_shifts: &shifts,
        thermal_cells: &cells,
        forcing: RotorIntervalForcing::Homogeneous,
    };
    let dt = 0.2;
    let next = advance_candidate(input, &thermal, dt, config(dt), || false).unwrap();
    let ratios = [1.0, 1.0, 10.0_f64.powf(-4.0 * 30.0 / 130.0)];
    let mut heat = [0.0; 3];
    for j in 0..3 {
        let expected =
            state.branch_stress[j] * (-dt / (branches[j].relaxation_time_s * ratios[j])).exp();
        heat[j] =
            (state.branch_stress[j].powi(2) - expected.powi(2)) / (2.0 * branches[j].stiffness);
        near(next.state().branch_stress[j], expected, 2e-11);
        near(next.report().branch_heat_j[j], heat[j], 2e-11);
    }
    let conducted = dt * 0.1 * (330.0 - 300.0);
    near(
        next.thermal().temperature_k(CellId(1)).unwrap(),
        300.0 + (heat[0] + heat[1] + conducted) / 4.0,
        2e-11,
    );
    near(
        next.thermal().temperature_k(CellId(2)).unwrap(),
        330.0 + (heat[2] - conducted) / 2.0,
        2e-11,
    );
    near(next.report().thermal.incoming_energy_j, 60.0, 1e-12);
    near(
        next.report().thermal.external_heat_j,
        heat.iter().sum(),
        2e-11,
    );
    near(next.report().external_thermal_heat_j, 0.0, 2e-11);
    near(next.report().combined_energy_residual_j, 0.0, 1e-8);
    assert_eq!(next.report().thermal.substeps, 1);
}

#[test]
fn thermal_event_splits_and_continuation_preserve_physical_harmonic_phase() {
    let law = Capacity(2.0);
    let thermal = timeline(
        &law,
        &[(1, 1.0, 0.0)],
        vec![
            TimedEvent {
                time_s: 0.0,
                event: Event::Boundary(CellId(1), boundary(2.0)),
            },
            TimedEvent {
                time_s: 0.31,
                event: Event::Boundary(CellId(1), boundary(5.0)),
            },
        ],
    );
    let (m, c, k, omega) = (1.0, 0.2, 1.7, 1.3);
    let model = RotorBearingModel::new(1, &[m], &[c], &[k], &[0.0], || false).unwrap();
    let real = k - m * omega * omega;
    let imaginary = c * omega;
    let q_re = real / (real * real + imaginary * imaginary);
    let q_im = -imaginary / (real * real + imaginary * imaginary);
    let state = RotorMaxwellState {
        displacement: vec![q_re],
        velocity: vec![-omega * q_im],
        branch_stress: vec![],
    };
    let rows = [RotorHarmonicForce {
        cosine: &[1.0],
        sine: &[0.0],
    }];
    let input = RotorThermalInput {
        time_s: 0.0,
        state: &state,
        model: &model,
        spin_rad_s: 0.0,
        branches: &[],
        temperature_shifts: &[],
        thermal_cells: &[],
        forcing: RotorIntervalForcing::Harmonic {
            rows: &rows,
            angular_frequency_rad_s: omega,
        },
    };
    let end = 1.13;
    let whole = advance_candidate(input, &thermal, end, config(0.13), || false).unwrap();
    let first = advance_candidate(input, &thermal, 0.37, config(0.07), || false).unwrap();
    let continued = advance_candidate(
        RotorThermalInput {
            time_s: first.thermal().time_s(),
            state: first.state(),
            ..input
        },
        first.thermal(),
        end,
        config(0.07),
        || false,
    )
    .unwrap();
    let phase = omega * end;
    let expected_q = q_re * phase.cos() - q_im * phase.sin();
    let expected_v = -omega * (q_re * phase.sin() + q_im * phase.cos());
    for candidate in [&whole, &continued] {
        near(candidate.state().displacement[0], expected_q, 2e-10);
        near(candidate.state().velocity[0], expected_v, 2e-10);
        near(
            candidate.thermal().temperature_k(CellId(1)).unwrap(),
            300.0 + (2.0 * 0.31 + 5.0 * (end - 0.31)) / 2.0,
            2e-11,
        );
        assert!(candidate.report().velocity_force_work_j > 0.0);
        assert_eq!(candidate.report().deposited_branch_heat_j, 0.0);
        near(candidate.report().combined_energy_residual_j, 0.0, 1e-8);
    }
    near(
        whole.report().external_thermal_heat_j,
        2.0 * 0.31 + 5.0 * (end - 0.31),
        2e-11,
    );
}

#[test]
fn driven_maxwell_self_heating_resonance_converges_in_motion_temperature_and_work() {
    let (mass, damping, stiffness, branch_k, tau_reference, capacity, amplitude) =
        (1.0, 0.03, 1.0, 0.6, 1.0, 0.05, 1.2);
    // The real reference-temperature dynamic stiffness vanishes here; finite
    // material/viscous loss retains a resolved, strongly driven resonance.
    let omega = (0.3 + 1.09_f64.sqrt()).sqrt();
    let end = 3.7;
    let oracle = |frozen: bool| {
        // Original coordinates [q,v,sigma,T,force work,C work,branch heat].
        // This independently integrates the continuously coupled physical ODE,
        // including temperature-dependent dissipation, without held-T splitting.
        let rhs = |time: f64, x: &[f64; 7]| {
            let delta = if frozen { 0.0 } else { x[3] - 300.0 };
            let tau = tau_reference * 10.0_f64.powf(-4.0 * delta / (100.0 + delta));
            let force = amplitude * (omega * time).cos();
            let heat = x[2] * x[2] / (branch_k * tau);
            [
                x[1],
                (force - damping * x[1] - stiffness * x[0] - x[2]) / mass,
                branch_k * x[1] - x[2] / tau,
                heat / capacity,
                force * x[1],
                damping * x[1] * x[1],
                heat,
            ]
        };
        let dt = end / 16384.0;
        let mut x = [0.0, 0.0, 0.0, 300.0, 0.0, 0.0, 0.0];
        for i in 0..16384 {
            let time = i as f64 * dt;
            let a = rhs(time, &x);
            let b = rhs(
                time + dt / 2.0,
                &std::array::from_fn(|j| x[j] + dt * a[j] / 2.0),
            );
            let c = rhs(
                time + dt / 2.0,
                &std::array::from_fn(|j| x[j] + dt * b[j] / 2.0),
            );
            let d = rhs(time + dt, &std::array::from_fn(|j| x[j] + dt * c[j]));
            for j in 0..7 {
                x[j] += dt * (a[j] + 2.0 * b[j] + 2.0 * c[j] + d[j]) / 6.0;
            }
        }
        x
    };
    let reference = oracle(false);
    let frozen = oracle(true);
    let law = Capacity(capacity);
    let thermal = timeline(&law, &[(1, 1.0, 0.0)], vec![]);
    let model =
        RotorBearingModel::new(1, &[mass], &[damping], &[stiffness], &[0.0], || false).unwrap();
    let branches = [MaxwellDampingBranch {
        strain: vec![1.0],
        stiffness: branch_k,
        relaxation_time_s: tau_reference,
    }];
    let shifts = [shift([290.0, 400.0])];
    let state = RotorMaxwellState {
        displacement: vec![0.0],
        velocity: vec![0.0],
        branch_stress: vec![0.0],
    };
    let rows = [RotorHarmonicForce {
        cosine: &[amplitude],
        sine: &[0.0],
    }];
    let input = RotorThermalInput {
        time_s: 0.0,
        state: &state,
        model: &model,
        spin_rad_s: 0.0,
        branches: &branches,
        temperature_shifts: &shifts,
        thermal_cells: &[CellId(1)],
        forcing: RotorIntervalForcing::Harmonic {
            rows: &rows,
            angular_frequency_rad_s: omega,
        },
    };
    let mut errors: Vec<[f64; 7]> = Vec::new();
    let mut finest = [0.0; 7];
    for dt in [0.08, 0.04, 0.02] {
        let next = advance_candidate(input, &thermal, end, config(dt), || false).unwrap();
        let current = next.state();
        let report = next.report();
        let observed = [
            current.displacement[0],
            current.velocity[0],
            current.branch_stress[0],
            next.thermal().temperature_k(CellId(1)).unwrap(),
            report.applied_work_j,
            report.velocity_force_work_j,
            report.branch_heat_j[0],
        ];
        errors.push(std::array::from_fn(|j| (observed[j] - reference[j]).abs()));
        finest = observed;
        near(observed[3], 300.0 + observed[6] / capacity, 2e-10);
        let mechanical_storage = 0.5 * mass * observed[1].powi(2)
            + 0.5 * stiffness * observed[0].powi(2)
            + observed[2].powi(2) / (2.0 * branch_k);
        near(
            mechanical_storage + observed[5] + observed[6] - observed[4],
            0.0,
            1e-8,
        );
        near(report.combined_energy_residual_j, 0.0, 1e-8);
        near(report.external_thermal_heat_j, 0.0, 1e-10);
        assert!(report.branch_heat_j[0] > 0.0);
    }
    for j in 0..7 {
        assert!(
            errors[1][j] < 0.8 * errors[0][j],
            "coordinate {j}: {errors:?}"
        );
        assert!(
            errors[2][j] < 0.8 * errors[1][j],
            "coordinate {j}: {errors:?}"
        );
        // Scale temperature accuracy by its actual rise, not the 300 K datum.
        let scale = if j == 3 {
            reference[j] - 300.0
        } else {
            reference[j].abs()
        };
        assert!(
            errors[2][j] < 0.02 * (1.0 + scale),
            "coordinate {j}: {errors:?}"
        );
    }
    assert!((finest[1] - frozen[1]).abs() > 0.05);
    assert!((finest[2] - frozen[2]).abs() > 0.2);
    assert!(finest[3] > frozen[3] + 1.5);
}

#[test]
fn unsupported_endpoint_work_exhaustion_and_late_cancellation_are_atomic() {
    let law = Capacity(1.0);
    let thermal = timeline(&law, &[(1, 1.0, 0.0)], vec![])
        .advance_to(0.0, config(0.1).thermal, || Ok(()))
        .unwrap()
        .0;
    let original_cells = thermal.cells().clone();
    let model = RotorBearingModel::new(1, &[1.0], &[0.2], &[1.0], &[0.0], || false).unwrap();
    let branches = [branch(2.0, 1.0)];
    let shifts = [shift([290.0, 400.0])];
    let narrow = [shift([299.0, 300.01])];
    let state = RotorMaxwellState {
        displacement: vec![0.0],
        velocity: vec![0.0],
        branch_stress: vec![10.0],
    };
    let input = RotorThermalInput {
        time_s: 0.0,
        state: &state,
        model: &model,
        spin_rad_s: 0.0,
        branches: &branches,
        temperature_shifts: &shifts,
        thermal_cells: &[CellId(1)],
        forcing: RotorIntervalForcing::Homogeneous,
    };
    assert!(matches!(
        advance_candidate(
            RotorThermalInput {
                temperature_shifts: &narrow,
                ..input
            },
            &thermal,
            0.2,
            config(0.1),
            || false
        ),
        Err(RotorThermalError::Temperature(_))
    ));
    let mut limited = config(0.05);
    limited.maximum_intervals = 1;
    assert!(matches!(
        advance_candidate(input, &thermal, 0.2, limited, || false),
        Err(RotorThermalError::Thermal(DepositionError::Budget(_)))
    ));
    assert!(
        advance_candidate(
            RotorThermalInput {
                thermal_cells: &[CellId(9)],
                ..input
            },
            &thermal,
            0.2,
            config(0.1),
            || false
        )
        .is_err()
    );
    let rows = [RotorHarmonicForce {
        cosine: &[1.0],
        sine: &[0.0],
    }];
    for forcing in [
        RotorIntervalForcing::Held(&[]),
        RotorIntervalForcing::Harmonic {
            rows: &rows,
            angular_frequency_rad_s: f64::NAN,
        },
    ] {
        assert!(matches!(
            advance_candidate(
                RotorThermalInput { forcing, ..input },
                &thermal,
                0.0,
                config(0.1),
                || false
            ),
            Err(RotorThermalError::Rotor(_))
        ));
    }
    let mut polls = 0;
    let _ = advance_candidate(input, &thermal, 0.2, config(0.05), || {
        polls += 1;
        false
    })
    .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut seen = 0;
        assert!(
            advance_candidate(input, &thermal, 0.2, config(0.05), || {
                seen += 1;
                seen == stop
            })
            .is_err()
        );
    }
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.cells(), &original_cells);
    assert_eq!(state.branch_stress, [10.0]);
    assert_eq!(branches[0].relaxation_time_s, 1.0);
}
