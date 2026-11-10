//! Force-driven startup and physical heating on original tetrahedral inventory.
use super::*;
use fs_couple::tensor_material::prepared::PreparedTensorControls;
use fs_couple::tensor_material::solid::prepared::PreparedTensorSolidBody;

const FREQUENCIES: [f64; 2] = [17.0, 29.0];

fn tone_body<'s>(
    source: &'s ViscoElements,
    basis: &[Vec<f64>],
    duration: f64,
    cx: &Cx<'_>,
    gate: &CancelGate,
) -> PreparedTensorSolidBody<'s> {
    let projection = source.project_basis(
        basis,
        SolidProjectionBudget {
            max_coordinates: 8,
            max_basis_scalars: 4096,
            max_material_points: 64,
            max_scalar_terms: 1_000_000,
        },
        cx,
    ).unwrap();
    PreparedTensorSolidBody::new_tones(
        projection, duration, &FREQUENCIES, 1.0,
        PreparedTensorControls {
            maximum_coordinates: 8,
            maximum_material_points: 64,
            maximum_state_coordinates: 64,
            maximum_heat_groups: 4,
            maximum_preparation_bytes: 64 * 1024 * 1024,
            maximum_preparation_multiply_terms: 1_000_000_000_000,
            maximum_step_multiply_terms: 100_000_000,
            maximum_relative_operator_error: 1e-9,
            maximum_work_error_j: 1e-10,
            maximum_steps: 128,
        },
        gate,
    ).unwrap()
}

fn nodal_tones(a: [f64; 2], b: [f64; 2]) -> (Vec<Vec<[f64; 3]>>, Vec<Vec<[f64; 3]>>) {
    let field = |amplitudes: [f64; 2]| amplitudes.into_iter().map(|value| {
        let mut force = vec![[0.0; 3]; 4];
        force[3][2] = value;
        force
    }).collect();
    (field(a), field(b))
}

// Independent original-coordinate equations for the one free Tet4 z DOF:
// m=integral rho*N3^2=1/30 kg, k=integral Czzzz=20 N/m.
// z is the physical branch strain, not an exported prepared state coordinate.
fn oracle_step(
    mut state: [f64; 5],
    start: f64,
    duration: f64,
    a: [f64; 2],
    b: [f64; 2],
    viscosity: f64,
    branch_stiffness: f64,
) -> [f64; 5] {
    let derivative = |t: f64, y: [f64; 5]| -> [f64; 5] {
        let force: f64 = (0..2).map(|i| {
            a[i] * (FREQUENCIES[i] * t).cos()
                + b[i] * (FREQUENCIES[i] * t).sin()
        }).sum();
        let [q, v, z, _, _] = y;
        [
            v,
            30.0 * (force - 20.0 * q - viscosity * v - branch_stiffness * z),
            v - z / 0.08,
            viscosity * v * v + branch_stiffness * z * z / 0.08,
            force * v,
        ]
    };
    let count = 4096;
    let h = duration / count as f64;
    for i in 0..count {
        let t = start + i as f64 * h;
        let k1 = derivative(t, state);
        let k2 = derivative(t + h / 2.0, core::array::from_fn(|j| state[j] + h * k1[j] / 2.0));
        let k3 = derivative(t + h / 2.0, core::array::from_fn(|j| state[j] + h * k2[j] / 2.0));
        let k4 = derivative(t + h, core::array::from_fn(|j| state[j] + h * k3[j]));
        state = core::array::from_fn(|j| {
            state[j] + h * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]) / 6.0
        });
    }
    state
}

#[test]
fn continuous_force_startup_and_retained_decay_heat_the_original_cell() {
    scope(|cx, gate| {
        let mat = material(0.01, 0.5);
        let solid = source(
            SolidMeshRef::Tet4 { nodes_m: &NODES, elements: &TETS },
            TetViscoField::Uniform(&mat), &FIXED, cx,
        );
        let dt = 0.0625;
        let body = tone_body(&solid, &[vec![1.0]], dt, cx, gate);
        let initial = body.history(&[0.0], &[0.0], InitialMemory::Relaxed, gate).unwrap();
        let chart = curve(1.0);
        let initial_heat = thermal(&chart, &[30], vec![], 32);
        let bindings = cells(&chart, &[30]);
        let (a, b) = ([0.08, -0.03], [0.04, 0.06]);
        let (mut cosine, mut sine) = nodal_tones(a, b);
        let mut accepted = initial.clone();
        let mut caloric = initial_heat.clone();
        let mut reference = [0.0; 5];

        for step in 0..4 {
            let old_heat = reference[3];
            let old_work = reference[4];
            let (drive_a, drive_b) = if step < 2 { (a, b) } else { ([0.0; 2], [0.0; 2]) };
            reference = oracle_step(reference, step as f64 * dt, dt, drive_a, drive_b, 0.2, 10.0);
            let next = body.advance_tones_thermal_candidate(
                &accepted, &caloric, &(), &cosine, &sine, &bindings, controls(), gate,
            ).unwrap();
            let report = next.report();
            near(next.solid().material_history().displacement_m()[0], reference[0], 2e-10);
            near(next.solid().material_history().velocity_mps()[0], reference[1], 2e-9);
            near(report.material_heat_j, reference[3] - old_heat, 2e-10);
            near(report.mechanical.motion.material.energy.supplied_work_j, reference[4] - old_work, 2e-10);
            near(report.mechanical.motion.material.energy.memory_after_j, 5.0 * reference[2].powi(2), 2e-10);
            near(next.thermal().temperature_k(CellId(30)).unwrap(), 300.0 + reference[3] / MASS, 2e-9);
            near(report.combined_energy_residual_j, 0.0, controls().energy_tolerance_j);
            near(report.heat_handoff_residual_j, 0.0, controls().energy_tolerance_j);
            assert_eq!(next.solid().material_history().steps(), (step + 1) as u64);
            assert_eq!(report.thermal.substeps, 2);
            assert!(report.material_heat_j > 0.0);
            if step >= 2 {
                near(report.mechanical.motion.material.energy.supplied_work_j, 0.0, 1e-14);
            }
            let (history, heat, _, _) = next.into_parts();
            accepted = history;
            caloric = heat;
            if step == 1 {
                (cosine, sine) = nodal_tones([0.0; 2], [0.0; 2]);
            } else {
                (cosine, sine) = body.advance_force_phase(&cosine, &sine, gate).unwrap();
            }
        }
        // Both scalar source integrations start from rest; their sum omits the
        // actual physical cross-tone loss present in the simultaneous trajectory.
        let together = oracle_step([0.0; 5], 0.0, dt, a, b, 0.2, 10.0)[3];
        let first = oracle_step([0.0; 5], 0.0, dt, [a[0], 0.0], [b[0], 0.0], 0.2, 10.0)[3];
        let second = oracle_step([0.0; 5], 0.0, dt, [0.0, a[1]], [0.0, b[1]], 0.2, 10.0)[3];
        assert!((together - first - second).abs() > 1e-6);
        assert_eq!(initial.material_history().time_s(), 0.0);
        assert_eq!(initial.material_history().dissipated_heat_j(), 0.0);
        assert_eq!(initial_heat.cells()[&CellId(30)].specific_h_j_kg, 0.0);
    });
}

#[test]
fn element_tone_heat_preserves_native_thermal_events_and_transport_hooks() {
    scope(|cx, gate| {
        let mut nodes = NODES.to_vec();
        nodes.extend(NODES.map(|mut p| { p[0] += 2.0; p }));
        let tets = [[0, 1, 2, 3], [4, 5, 6, 7]];
        let fixed: Vec<_> = (0..24).filter(|&d| d != 11 && d != 23).collect();
        let materials = [material(0.01, 0.0), material(0.0, 0.5)];
        let solid = source(
            SolidMeshRef::Tet4 { nodes_m: &nodes, elements: &tets },
            TetViscoField::PerElement(&materials), &fixed, cx,
        );
        let dt = 0.0625;
        let body = tone_body(&solid, &[vec![1.0, 0.0], vec![0.0, 1.0]], dt, cx, gate);
        let initial = body.history(&[0.001, 0.002], &[0.02, -0.01], InitialMemory::Instantaneous, gate).unwrap();
        let chart = curve(1.0);
        let switch = 0.0234375;
        let ambient = Boundary { ambient_k: 310.0, htc_w_m2_k: 0.7, ..boundary(0.0) };
        let reservoir = Boundary {
            reservoir: Some(Reservoir { temperature_k: 320.0, area_m2: 0.2, conductance_w_k: 9.0 }),
            ..boundary(0.0)
        };
        let heat = thermal(&chart, &[30, 10], vec![
            TimedEvent { time_s: 0.0, event: Event::Contact(ThermalLink {
                id: 55, a: CellId(30), b: CellId(10), area_m2: 0.1, conductance_w_k: 0.0,
            }) },
            TimedEvent { time_s: 0.0, event: Event::Boundary(CellId(30), ambient) },
            TimedEvent { time_s: 0.0, event: Event::Boundary(CellId(10), reservoir) },
            TimedEvent { time_s: switch, event: Event::Boundary(CellId(30), Boundary { source_w: 0.004, ..ambient }) },
        ], 16);
        let a = [[0.08, -0.03], [-0.02, 0.06]];
        let b = [[0.04, 0.06], [0.05, -0.02]];
        let mut cosine = vec![vec![[0.0; 3]; 8]; 2];
        let mut sine = cosine.clone();
        for tone in 0..2 {
            for (element, node) in [3, 7].into_iter().enumerate() {
                cosine[tone][node][2] = a[element][tone];
                sine[tone][node][2] = b[element][tone];
            }
        }
        let reference = [
            oracle_step([0.001, 0.02, 0.001, 0.0, 0.0], 0.0, dt, a[0], b[0], 0.2, 0.0),
            oracle_step([0.002, -0.01, 0.002, 0.0, 0.0], 0.0, dt, a[1], b[1], 0.0, 10.0),
        ];
        let observer = OtherPhysics::default();
        let candidate = body.advance_tones_thermal_candidate(
            &initial, &heat, &observer, &cosine, &sine, &cells(&chart, &[30, 10]), controls(), gate,
        ).unwrap();
        let report = candidate.report();
        let obs = candidate.observer();
        assert_eq!(candidate.solid().material_history().steps(), 1);
        assert_eq!(report.thermal.substeps, 5);
        assert_eq!(obs.events, 1);
        assert_eq!((obs.sources, obs.conductances, obs.transfers, obs.limits, obs.steps), (5, 5, 5, 5, 5));
        let mut temperatures = [300.0; 2];
        let mut external = 0.0;
        for &(start, h) in &obs.intervals {
            let transfer = 0.5 * (temperatures[1] - temperatures[0]) + 0.0002;
            let outside_a = 0.2 * (310.0 - temperatures[0]) + 0.001
                + if start >= switch { 0.004 } else { 0.0 };
            let outside_b = 0.3 * (320.0 - temperatures[1]) + 0.002;
            temperatures[0] += h * (reference[0][3] / dt + outside_a + transfer) / MASS;
            temperatures[1] += h * (reference[1][3] / dt + outside_b - transfer) / MASS;
            external += h * (outside_a + outside_b);
        }
        for (i, id) in [30, 10].into_iter().enumerate() {
            near(report.element_heat_j[i], reference[i][3], 2e-10);
            near(candidate.solid().material_history().displacement_m()[i], reference[i][0], 2e-10);
            near(candidate.thermal().temperature_k(CellId(id)).unwrap(), temperatures[i], 2e-9);
        }
        near(report.other_external_heat_j, external, 2e-10);
        near(report.combined_energy_residual_j, 0.0, controls().energy_tolerance_j);
        assert!((report.element_heat_j[0] - report.element_heat_j[1]).abs() > 1e-8);
        assert_eq!(candidate.thermal().links()[&55].conductance_w_k, 0.0);
        assert_eq!(candidate.thermal().cells()[&CellId(30)].boundary.htc_w_m2_k, 0.7);
        assert_eq!(observer.steps, 0);
        assert_eq!(initial.material_history().steps(), 0);
    });
}

#[derive(Clone, Default)]
struct RefuseLate { steps: usize }
impl DepositionObserver for RefuseLate {
    type Error = DepositionError;
    fn step(
        &mut self, _: f64, _: f64, _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.steps += 1;
        if thermal.time_s() > 0.04 {
            Err(DepositionError::InvalidInput("synthetic late tone-heat refusal"))
        } else { Ok(()) }
    }
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        Ok(())
    }
}

#[test]
fn foreign_owner_late_refusal_and_retry_preserve_both_accepted_histories() {
    scope(|cx, gate| {
        let mat = material(0.01, 0.5);
        let solid = source(
            SolidMeshRef::Tet4 { nodes_m: &NODES, elements: &TETS },
            TetViscoField::Uniform(&mat), &FIXED, cx,
        );
        let body = tone_body(&solid, &[vec![1.0]], 0.125, cx, gate);
        let foreign = tone_body(&solid, &[vec![1.0]], 0.125, cx, gate);
        let initial = body.history(&[0.001], &[0.03], InitialMemory::Instantaneous, gate).unwrap();
        let saved = initial.material_history().state().to_vec();
        let chart = curve(1.0);
        let heat = thermal(&chart, &[30], vec![], 16);
        let bindings = cells(&chart, &[30]);
        let (cosine, sine) = nodal_tones([0.08, -0.03], [0.04, 0.06]);
        assert!(foreign.advance_tones_thermal_candidate(
            &initial, &heat, &(), &cosine, &sine, &bindings, controls(), gate,
        ).is_err());
        let chart_copy = chart.clone();
        assert!(body.advance_tones_thermal_candidate(
            &initial, &heat, &(), &cosine, &sine, &cells(&chart_copy, &[30]), controls(), gate,
        ).is_err());
        let mut narrow = bindings.clone();
        narrow[0].temperature_interval_k = [300.0, 300.000000001];
        assert!(body.advance_tones_thermal_candidate(
            &initial, &heat, &(), &cosine, &sine, &narrow, controls(), gate,
        ).is_err());
        let short = thermal(&chart, &[30], vec![], 1);
        assert!(matches!(body.advance_tones_thermal_candidate(
            &initial, &short, &(), &cosine, &sine, &bindings, controls(), gate,
        ), Err(SolidThermalError::Thermal(DepositionError::Budget(_)))));
        let observer = RefuseLate::default();
        assert!(matches!(body.advance_tones_thermal_candidate(
            &initial, &heat, &observer, &cosine, &sine, &bindings, controls(), gate,
        ), Err(SolidThermalError::Observer(_))));
        assert_eq!(observer.steps, 0);
        assert_eq!(initial.material_history().state(), saved.as_slice());
        assert_eq!(initial.material_history().time_s(), 0.0);
        assert_eq!(initial.material_history().dissipated_heat_j(), 0.0);
        assert_eq!(heat.cells()[&CellId(30)].specific_h_j_kg, 0.0);

        let retry = body.advance_tones_thermal_candidate(
            &initial, &heat, &(), &cosine, &sine, &bindings, controls(), gate,
        ).unwrap();
        let replay = body.advance_tones_thermal_candidate(
            &initial, &heat, &(), &cosine, &sine, &bindings, controls(), gate,
        ).unwrap();
        assert_eq!(retry.solid().material_history().state(), replay.solid().material_history().state());
        assert_eq!(retry.thermal().cells()[&CellId(30)].specific_h_j_kg,
            replay.thermal().cells()[&CellId(30)].specific_h_j_kg);
        assert!(retry.report().deposited_heat_j > 0.0);
        assert!(body.advance_tones_thermal_candidate(
            retry.solid(), &heat, &(), &cosine, &sine, &bindings, controls(), gate,
        ).is_err());
        gate.request();
        assert!(matches!(body.advance_tones_thermal_candidate(
            &initial, &heat, &(), &cosine, &sine, &bindings, controls(), gate,
        ), Err(SolidThermalError::Cancelled)));
        assert_eq!(initial.material_history().state(), saved.as_slice());
    });
}
