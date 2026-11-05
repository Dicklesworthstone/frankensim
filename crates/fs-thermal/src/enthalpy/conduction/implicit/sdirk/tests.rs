use super::super::adaptive::{AdaptiveImplicitConductionControls, AdaptiveImplicitConductionError};
use super::*;
use crate::enthalpy::advection::reconstruction::{AxialCellStencil, LimitedEnthalpyReconstruction};
use crate::enthalpy::advection::{Advected, MassBoundary, MassFlow};

#[path = "source_schedule_tests.rs"]
mod source_schedule;

struct Law {
    cp: f64,
    k: f64,
    datum: f64,
    nonlinear: bool,
    latent: bool,
}
impl Law {
    fn constant(cp: f64, k: f64) -> Self {
        Self {
            cp,
            k,
            datum: 0.0,
            nonlinear: false,
            latent: false,
        }
    }
    fn coordinate(&self, h: f64) -> f64 {
        (h - self.datum) / self.cp
    }
}
impl ThermalLaw for Law {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let x = self.coordinate(h);
        if !x.is_finite() || !(100.0..=1000.0).contains(&x) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(if !self.latent || x < 300.0 {
            x
        } else if x <= 500.0 {
            300.0
        } else {
            x - 200.0
        })
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / self.cp
    }
    fn specific_entropy_change_j_kg_k(
        &self,
        a: f64,
        b: f64,
    ) -> Result<Option<f64>, DepositionError> {
        let primitive = |h| -> Result<f64, DepositionError> {
            let t = self.temperature_k(h)?;
            let x = self.coordinate(h);
            Ok(self.cp
                * if !self.latent || x < 300.0 {
                    t.ln()
                } else if x <= 500.0 {
                    300.0_f64.ln() + (x - 300.0) / 300.0
                } else {
                    t.ln() + 200.0 / 300.0
                })
        };
        Ok(Some(primitive(b)? - primitive(a)?))
    }
}
impl ConductivityLaw for Law {
    fn conductivity_w_m_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = self.temperature_k(h)?;
        Ok(if self.nonlinear {
            self.k * t / 300.0
        } else {
            self.k
        })
    }
}
impl ImplicitConductionLaw for Law {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self
    }
    fn conductivity_law(&self) -> &dyn ConductivityLaw {
        self
    }
    fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        [self.datum + 100.0 * self.cp, self.datum + 1000.0 * self.cp]
    }
    fn state_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<ImplicitConductionState, DepositionError> {
        let slope = if self.latent && (300.0..=500.0).contains(&self.coordinate(h)) {
            0.0
        } else {
            1.0 / self.cp
        };
        Ok(ImplicitConductionState {
            temperature_k: self.temperature_k(h)?,
            temperature_slope_kg_k_j: slope,
            conductivity_w_m_k: self.conductivity_w_m_k(h)?,
            conductivity_slope_kg_w_m_k_j: if self.nonlinear {
                self.k * slope / 300.0
            } else {
                0.0
            },
        })
    }
}

fn id(i: usize) -> CellId {
    CellId(7 + 3 * i as u64)
}
fn boundary(bath: Option<f64>) -> Boundary {
    Boundary {
        ambient_k: 300.0,
        htc_w_m2_k: 0.0,
        source_w: 0.0,
        reservoir: bath.map(|temperature_k| Reservoir {
            area_m2: 1.0,
            conductance_w_k: 0.0,
            temperature_k,
        }),
    }
}
fn network<'a>(
    law: &'a Law,
    coordinates: &[f64],
    masses: &[f64],
    edges: &[(usize, usize, f64)],
    bath: Option<f64>,
    future: Vec<TimedEvent>,
) -> (DepositionSimulation<'a>, ConductionObserver<'a>) {
    let mut events: Vec<_> = coordinates
        .iter()
        .enumerate()
        .map(|(i, &x)| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: id(i),
                material_id: MaterialId(0),
                body_id: i as u64,
                segment_id: 0,
                mass_kg: masses[i],
                specific_h_j_kg: law.datum + law.cp * x,
                geometry: CellGeometry {
                    start_m: [i as f64, 0.0, 0.0],
                    end_m: [i as f64 + 1.0, 0.0, 0.0],
                    volume_m3: masses[i],
                    surface_area_m2: 4.0,
                },
                boundary: boundary(if i == 0 { bath } else { None }),
            }),
        })
        .collect();
    let mut links = Vec::new();
    for (i, &(a, b, distance)) in edges.iter().enumerate() {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: i as u64,
                a: id(a),
                b: id(b),
                area_m2: 1.0,
                conductance_w_k: 0.0,
            }),
        });
        links.push(ConductiveLink {
            link_id: i as u64,
            half_distances_m: [0.5 * distance; 2],
            contact_resistance_m2_k_w: 0.0,
        });
    }
    events.extend(future);
    let initial = DepositionSimulation::new(
        events,
        &[(MaterialId(0), law)],
        DepositionBudget {
            max_cells: 128,
            max_links: 256,
            max_events: 512,
            max_substeps: 100_000,
        },
    )
    .unwrap();
    let (thermal, _) = initial
        .advance_to(
            0.0,
            StepConfig {
                max_dt_s: 1.0,
                energy_tolerance_j: 0.0,
            },
            || Ok(()),
        )
        .unwrap();
    let reservoirs = bath.map(|_| ConductiveReservoir {
        cell_id: id(0),
        half_distance_m: 1.0,
        surface_resistance_m2_k_w: 0.0,
    });
    let observer = ConductionObserver::new(
        &thermal,
        &[(MaterialId(0), law)],
        &links,
        &reservoirs.into_iter().collect::<Vec<_>>(),
        || Ok(()),
    )
    .unwrap();
    (thermal, observer)
}
fn chain<'a>(
    thermal: &DepositionSimulation<'_>,
    observer: &ConductionObserver<'a>,
    law: &'a Law,
) -> ImplicitConductionChain<'a> {
    ImplicitConductionChain::new(thermal, observer, &[(MaterialId(0), law)], || Ok(())).unwrap()
}
fn controls(dt: f64) -> ImplicitConductionControls {
    ImplicitConductionControls {
        maximum_dt_s: dt,
        cell_tolerance_j: 1e-10,
        total_tolerance_j: 1e-8,
        ..Default::default()
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.16e} != {expected:.16e}"
    );
}
fn amplification(decay_dt: f64) -> f64 {
    (1.0 - (1.0 - 2.0 * GAMMA) * decay_dt) / (1.0 + GAMMA * decay_dt).powi(2)
}

fn throughflow(cells: usize, kg_s: f64, incoming_h: f64) -> Advected<()> {
    Advected::new(
        (),
        (0..cells - 1)
            .map(|i| MassFlow {
                id: i as u64,
                a: id(i),
                b: id(i + 1),
                kg_s,
            })
            .collect(),
        cells,
    )
    .unwrap()
    .with_boundaries(
        vec![
            MassBoundary::Inlet {
                id: cells as u64,
                cell: id(0),
                kg_s,
                material_id: MaterialId(0),
                specific_h_j_kg: incoming_h,
            },
            MassBoundary::Outlet {
                id: cells as u64 + 1,
                cell: id(cells - 1),
                kg_s,
            },
        ],
        cells + 1,
    )
    .unwrap()
}

#[test]
fn unequal_capacities_conserve_native_heat_and_match_two_stage_relaxation() {
    let law = Law::constant(3.0, 1.0);
    let masses = [2.0, 5.0];
    let (thermal, observer) = network(&law, &[600.0, 200.0], &masses, &[(0, 1, 1.0)], None, vec![]);
    let direct = chain(&thermal, &observer, &law);
    let graph = ImplicitConductionGraph::new(
        &thermal,
        &observer,
        &[(MaterialId(0), &law)],
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
    let dt = 4.0;
    let capacity = [masses[0] * law.cp, masses[1] * law.cp];
    let mean = (capacity[0] * 600.0 + capacity[1] * 200.0) / capacity.iter().sum::<f64>();
    let decay_dt = dt * (1.0 / capacity[0] + 1.0 / capacity[1]);
    let stage1 = [
        mean + (600.0 - mean) / (1.0 + GAMMA * decay_dt),
        mean + (200.0 - mean) / (1.0 + GAMMA * decay_dt),
    ];
    let expected = [
        mean + (600.0 - mean) * amplification(decay_dt),
        mean + (200.0 - mean) * amplification(decay_dt),
    ];
    let production = |t: [f64; 2]| (t[0] - t[1]).powi(2) / (t[0] * t[1]);
    let expected_generation =
        dt * ((1.0 - GAMMA) * production(stage1) + GAMMA * production(expected));
    for sparse in [false, true] {
        let (next, report) = if sparse {
            graph.advance_sdirk2_to(&thermal, dt, controls(dt), &mut workspace, || Ok(()))
        } else {
            direct.advance_sdirk2_to(&thermal, dt, controls(dt), &mut workspace, || Ok(()))
        }
        .unwrap();
        for (i, &t) in expected.iter().enumerate() {
            near(next.temperature_k(id(i)).unwrap(), t, 1e-10);
        }
        near(report.energy.external_heat_j, 0.0, 0.0);
        near(report.energy.stored_energy_change_j, 0.0, 2e-10);
        near(
            report.entropy.conduction_generation_j_k,
            expected_generation,
            1e-11,
        );
        near(
            report.entropy.thermal_storage_change_j_k,
            capacity[0] * (expected[0] / 600.0).ln() + capacity[1] * (expected[1] / 200.0).ln(),
            1e-12,
        );
        assert!(!report.entropy.has_unclassified_sources);
        assert_eq!(report.energy.substeps, 1);
        assert!(report.maximum_cell_residual_j <= 1e-10);
        assert_eq!(next.time_s(), dt);
    }
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.temperature_k(id(0)).unwrap(), 600.0);
}

#[test]
fn bath_heat_and_entropy_use_both_physical_stage_temperatures() {
    let law = Law::constant(3.0, 1.0);
    let (mut thermal, observer) = network(&law, &[600.0], &[2.0], &[], Some(350.0), vec![]);
    // Three exposed square metres give 0.5 W/K to a separate 250 K ambient;
    // the material-backed reservoir contributes 1 W/K at 350 K.
    thermal.cells.get_mut(&id(0)).unwrap().boundary.ambient_k = 250.0;
    thermal.cells.get_mut(&id(0)).unwrap().boundary.htc_w_m2_k = 1.0 / 6.0;
    let solver = chain(&thermal, &observer, &law);
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    let dt = 3.0;
    let equilibrium = (350.0 + 0.5 * 250.0) / 1.5;
    let t1 = equilibrium + (600.0 - equilibrium) / (1.0 + GAMMA * dt * 1.5 / 6.0);
    let t2 = equilibrium + (600.0 - equilibrium) * amplification(dt * 1.5 / 6.0);
    let power = |t: f64| (350.0 - t) + 0.5 * (250.0 - t);
    let bath_entropy = |t: f64| (350.0 - t) / 350.0 + 0.5 * (250.0 - t) / 250.0;
    let heat = dt * ((1.0 - GAMMA) * power(t1) + GAMMA * power(t2));
    let production =
        |t: f64| (t - 350.0).powi(2) / (350.0 * t) + 0.5 * (t - 250.0).powi(2) / (250.0 * t);
    let (next, report) = solver
        .advance_sdirk2_to(&thermal, dt, controls(dt), &mut workspace, || Ok(()))
        .unwrap();
    near(next.temperature_k(id(0)).unwrap(), t2, 1e-11);
    near(report.energy.external_heat_j, heat, 1e-10);
    near(
        report.energy.stored_energy_change_j,
        6.0 * (t2 - 600.0),
        1e-10,
    );
    near(
        report.entropy.reservoir_entropy_input_j_k,
        dt * ((1.0 - GAMMA) * bath_entropy(t1) + GAMMA * bath_entropy(t2)),
        1e-12,
    );
    near(
        report.entropy.conduction_generation_j_k,
        dt * ((1.0 - GAMMA) * production(t1) + GAMMA * production(t2)),
        1e-12,
    );
    near(
        report.entropy.thermal_storage_change_j_k,
        6.0 * (t2 / 600.0).ln(),
        1e-12,
    );
    near(
        report.entropy.temporal_defect_j_k,
        report.entropy.thermal_storage_change_j_k
            - report.entropy.reservoir_entropy_input_j_k
            - report.entropy.conduction_generation_j_k,
        1e-15,
    );
    assert!(!report.entropy.has_unclassified_sources);
    // The accepted endpoint retains its authored reservoir, not a stage weight.
    assert_eq!(
        next.cells()[&id(0)].boundary,
        thermal.cells()[&id(0)].boundary
    );
}

#[test]
fn nonlinear_material_diffusion_has_second_order_time_refinement() {
    let law = Law {
        nonlinear: true,
        ..Law::constant(1.0, 1.0)
    };
    let (thermal, observer) = network(&law, &[600.0], &[1.0], &[], Some(300.0), vec![]);
    let solver = chain(&thermal, &observer, &law);
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    // dT/dt = -T*(T-300)/300 has (T-300)/T = 0.5*exp(-t).
    let exact = 300.0 / (1.0 - 0.5 * (-1.0_f64).exp());
    let mut errors = Vec::new();
    for steps in [8, 16, 32] {
        let (next, report) = solver
            .advance_sdirk2_to(
                &thermal,
                1.0,
                controls(1.0 / steps as f64),
                &mut workspace,
                || Ok(()),
            )
            .unwrap();
        errors.push((next.temperature_k(id(0)).unwrap() - exact).abs());
        assert_eq!(report.energy.substeps, steps);
        assert!(report.energy.residual_j.abs() <= 1e-8);
        assert!(!report.entropy.has_unclassified_sources);
    }
    for pair in errors.windows(2) {
        assert!(
            (3.7..4.4).contains(&(pair[0] / pair[1])),
            "nonlinear errors: {errors:?}"
        );
    }
}

#[test]
fn cyclic_convection_diffusion_refines_at_second_order_and_preserves_mean() {
    let law = Law::constant(1.0, 0.1);
    let n = 12;
    let dx = 1.0 / n as f64;
    let shape: Vec<_> = (0..n)
        .map(|i| (2.0 * core::f64::consts::PI * (i as f64 + 0.5) / n as f64).cos())
        .collect();
    let temperature: Vec<_> = shape.iter().map(|&v| 400.0 + 50.0 * v).collect();
    let edges: Vec<_> = (0..n).map(|i| (i, (i + 1) % n, dx)).collect();
    let (thermal, observer) = network(&law, &temperature, &vec![dx; n], &edges, None, vec![]);
    let rate = 0.6;
    let flow = Advected::new(
        (),
        (0..n)
            .map(|i| MassFlow {
                id: i as u64,
                a: id(i),
                b: id((i + 1) % n),
                kg_s: rate,
            })
            .collect(),
        n,
    )
    .unwrap();
    let solver = ImplicitConductionGraph::new(
        &thermal,
        &observer,
        &[(MaterialId(0), &law)],
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap()
    .with_implicit_advection(&thermal, &flow, || Ok(()))
    .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(n).unwrap();
    let duration = 0.25;
    // Exact complex eigenvalue of the periodic upwind/conduction operator
    // isolates time error from spatial diffusion and phase error. The directed
    // part translates the Fourier mode; a symmetric-conductance surrogate fails.
    let theta = core::f64::consts::TAU / n as f64;
    let lambda =
        4.0 * law.k * (0.5 * theta).sin().powi(2) / dx.powi(2) + rate * (1.0 - theta.cos()) / dx;
    let frequency = rate * theta.sin() / dx;
    let mut errors = Vec::new();
    for steps in [4, 8, 16] {
        let (next, report) = solver
            .advance_sdirk2_to(
                &thermal,
                duration,
                controls(duration / steps as f64),
                &mut workspace,
                || Ok(()),
            )
            .unwrap();
        let mut error = 0.0_f64;
        let mut mean = 0.0;
        for i in 0..n {
            let actual = next.temperature_k(id(i)).unwrap();
            let exact = 400.0
                + 50.0
                    * (-lambda * duration).exp()
                    * (theta * (i as f64 + 0.5) - frequency * duration).cos();
            error = error.max((actual - exact).abs());
            mean += actual / n as f64;
        }
        errors.push(error);
        near(mean, 400.0, 1e-10);
        near(report.energy.external_heat_j, 0.0, 0.0);
        assert!(report.energy.residual_j.abs() <= 1e-8);
        assert!(report.linear_iterations > 0);
        assert!(report.entropy.has_unclassified_sources);
        assert_eq!(report.energy.substeps, steps);
    }
    for pair in errors.windows(2) {
        assert!(
            (3.7..4.5).contains(&(pair[0] / pair[1])),
            "Fourier errors: {errors:?}"
        );
    }
}

#[test]
fn unresolved_stiff_modes_decay_without_clipping_the_signed_amplification() {
    let law = Law::constant(1.0, 1.0);
    let (thermal, observer) = network(
        &law,
        &[400.0, 200.0],
        &[1.0, 1.0],
        &[(0, 1, 1.0)],
        None,
        vec![],
    );
    let solver = chain(&thermal, &observer, &law);
    let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
    let mut deviations = Vec::new();
    for dt in [100.0, 1000.0] {
        let (next, report) = solver
            .advance_sdirk2_to(&thermal, dt, controls(dt), &mut workspace, || Ok(()))
            .unwrap();
        let signed = next.temperature_k(id(0)).unwrap() - 300.0;
        near(signed, 100.0 * amplification(2.0 * dt), 1e-10);
        assert!(signed < 0.0); // L-stable does not imply monotone modal decay.
        deviations.push(signed.abs());
        assert!(report.entropy.thermal_storage_change_j_k > 0.0);
        assert_eq!(report.energy.substeps, 1);
        assert_eq!(report.rejected_intervals, 0);
    }
    assert!(deviations[1] < 0.11 * deviations[0]);
}

#[test]
fn latent_plateau_retains_original_enthalpy_storage_and_entropy() {
    let law = Law {
        latent: true,
        datum: 1_000_000.0,
        ..Law::constant(1.0, 1.0)
    };
    let (thermal, observer) = network(&law, &[450.0], &[1.0], &[], Some(350.0), vec![]);
    let solver = chain(&thermal, &observer, &law);
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    let (next, report) = solver
        .advance_sdirk2_to(&thermal, 0.5, controls(0.5), &mut workspace, || Ok(()))
        .unwrap();
    near(
        next.cells()[&id(0)].specific_h_j_kg,
        law.datum + 475.0,
        1e-9,
    );
    near(next.temperature_k(id(0)).unwrap(), 300.0, 0.0);
    near(report.energy.external_heat_j, 25.0, 1e-12);
    near(report.energy.stored_energy_change_j, 25.0, 1e-9);
    near(
        report.entropy.thermal_storage_change_j_k,
        25.0 / 300.0,
        1e-12,
    );
    near(report.entropy.temporal_defect_j_k, 0.0, 1e-12);
    assert!(!report.entropy.has_unclassified_sources);
}

#[test]
fn stage_two_domain_and_work_refusals_do_not_publish_first_stage() {
    let law = Law::constant(1.0, 1.0);
    let (thermal, observer) = network(&law, &[800.0], &[1.0], &[], Some(100.0), vec![]);
    let solver = chain(&thermal, &observer, &law);
    let before = thermal.cells().clone();
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    let result = solver.advance_sdirk2_to(
        &thermal,
        10.0,
        ImplicitConductionControls {
            maximum_halvings: 0,
            ..controls(10.0)
        },
        &mut workspace,
        || Ok(()),
    );
    // The first stage is in-domain; the signed SDIRK endpoint would be below
    // the caloric domain. No positivity clamp or first-stage partial commit.
    assert!(matches!(
        result,
        Err(ImplicitConductionError::NotConverged { .. })
    ));
    assert_eq!(thermal.cells(), &before);
    assert_eq!(thermal.time_s(), 0.0);

    let (recovered, report) = solver
        .advance_sdirk2_to(&thermal, 10.0, controls(10.0), &mut workspace, || Ok(()))
        .unwrap();
    assert!(report.rejected_intervals > 0);
    assert!(report.energy.substeps > 1);
    assert_eq!(recovered.time_s(), 10.0);
    assert!(recovered.temperature_k(id(0)).unwrap() >= 100.0);
    near(
        report.energy.external_heat_j,
        recovered.cells()[&id(0)].specific_h_j_kg - before[&id(0)].specific_h_j_kg,
        1e-8,
    );
    assert_eq!(thermal.cells(), &before);

    let (thermal, observer) = network(&law, &[600.0], &[1.0], &[], Some(300.0), vec![]);
    let solver = chain(&thermal, &observer, &law);
    let before = thermal.cells().clone();
    let result = solver.advance_sdirk2_to(
        &thermal,
        1.0,
        ImplicitConductionControls {
            maximum_residual_evaluations: 3,
            ..controls(1.0)
        },
        &mut workspace,
        || Ok(()),
    );
    assert!(matches!(
        result,
        Err(ImplicitConductionError::Thermal(DepositionError::Budget(
            "implicit residual evaluations"
        )))
    ));
    assert_eq!(thermal.cells(), &before);
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn constant_sources_preserve_two_cell_energy_and_native_entropy_quadrature() {
    let law = Law::constant(3.0, 1.0);
    let masses = [2.0, 5.0];
    let powers = [30.0, -5.0];
    let (mut thermal, observer) =
        network(&law, &[300.0, 450.0], &masses, &[(0, 1, 1.0)], None, vec![]);
    for (i, &power) in powers.iter().enumerate() {
        thermal.cells.get_mut(&id(i)).unwrap().boundary.source_w = power;
    }
    let direct = chain(&thermal, &observer, &law);
    let graph = ImplicitConductionGraph::new(
        &thermal,
        &observer,
        &[(MaterialId(0), &law)],
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
    let capacity = [law.cp * masses[0], law.cp * masses[1]];
    let total_capacity = capacity[0] + capacity[1];
    let total_power = powers[0] + powers[1];
    let initial_mean = (capacity[0] * 300.0 + capacity[1] * 450.0) / total_capacity;
    let decay = 1.0 / capacity[0] + 1.0 / capacity[1];
    let equilibrium_difference = (powers[0] / capacity[0] - powers[1] / capacity[1]) / decay;
    let dt = 4.0;
    let temperatures = |mean: f64, difference: f64| {
        [
            mean + capacity[1] / total_capacity * difference,
            mean - capacity[0] / total_capacity * difference,
        ]
    };
    let first = temperatures(
        initial_mean + GAMMA * dt * total_power / total_capacity,
        equilibrium_difference + (-150.0 - equilibrium_difference) / (1.0 + GAMMA * dt * decay),
    );
    let expected = temperatures(
        initial_mean + dt * total_power / total_capacity,
        equilibrium_difference + (-150.0 - equilibrium_difference) * amplification(dt * decay),
    );
    let source_entropy = |t: [f64; 2]| powers[0] / t[0] + powers[1] / t[1];
    let production = |t: [f64; 2]| (t[0] - t[1]).powi(2) / (t[0] * t[1]);
    for sparse in [false, true] {
        let (next, report) = if sparse {
            graph.advance_sdirk2_to(&thermal, dt, controls(dt), &mut workspace, || Ok(()))
        } else {
            direct.advance_sdirk2_to(&thermal, dt, controls(dt), &mut workspace, || Ok(()))
        }
        .unwrap();
        for (i, &temperature) in expected.iter().enumerate() {
            near(next.temperature_k(id(i)).unwrap(), temperature, 2e-10);
            near(next.cells()[&id(i)].boundary.source_w, powers[i], 0.0);
        }
        near(report.energy.external_heat_j, total_power * dt, 1e-12);
        near(
            report.energy.stored_energy_change_j,
            total_power * dt,
            2e-10,
        );
        near(
            report.entropy.unclassified_source_at_cell_j_k,
            dt * ((1.0 - GAMMA) * source_entropy(first) + GAMMA * source_entropy(expected)),
            1e-11,
        );
        near(
            report.entropy.conduction_generation_j_k,
            dt * ((1.0 - GAMMA) * production(first) + GAMMA * production(expected)),
            1e-11,
        );
        assert!(report.entropy.has_unclassified_sources);
        assert_eq!(report.energy.substeps, 1);
    }
}

#[test]
fn steady_heating_to_a_bath_has_second_order_and_adaptive_temperature_accuracy() {
    use super::super::adaptive::AdaptiveImplicitConductionControls;
    let law = Law::constant(1.0, 1.0);
    let (mut thermal, observer) = network(&law, &[400.0], &[1.0], &[], Some(300.0), vec![]);
    thermal.cells.get_mut(&id(0)).unwrap().boundary.source_w = 50.0;
    let solver = chain(&thermal, &observer, &law);
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    // C dT/dt = P - G(T - Tb), with C=G=1, P=50, Tb=300.
    let exact = 350.0 + 50.0 * (-1.0_f64).exp();
    let mut errors = Vec::new();
    for steps in [4, 8, 16] {
        let (next, report) = solver
            .advance_sdirk2_to(
                &thermal,
                1.0,
                controls(1.0 / steps as f64),
                &mut workspace,
                || Ok(()),
            )
            .unwrap();
        errors.push((next.temperature_k(id(0)).unwrap() - exact).abs());
        near(
            report.energy.external_heat_j,
            next.temperature_k(id(0)).unwrap() - 400.0,
            1e-9,
        );
        assert!(report.entropy.has_unclassified_sources);
        assert_eq!(report.energy.substeps, steps);
    }
    for pair in errors.windows(2) {
        assert!(
            (3.7..4.4).contains(&(pair[0] / pair[1])),
            "heated errors: {errors:?}"
        );
    }
    let (next, report) = solver
        .advance_sdirk2_adaptive_to(
            &thermal,
            1.0,
            AdaptiveImplicitConductionControls {
                step: controls(1.0),
                initial_dt_s: 1.0,
                enthalpy_absolute_tolerance_j_kg: 1e-5,
                temperature_absolute_tolerance_k: 1e-5,
                relative_tolerance: 0.0,
                ..Default::default()
            },
            &mut workspace,
            || Ok(()),
        )
        .unwrap();
    near(next.temperature_k(id(0)).unwrap(), exact, 4e-4);
    near(
        report.solve.energy.external_heat_j,
        next.temperature_k(id(0)).unwrap() - 400.0,
        1e-8,
    );
    assert!(report.accuracy_rejections > 0);
    assert!(report.solve.entropy.has_unclassified_sources);
    assert!(report.maximum_accepted_normalized_error <= 1.0);
}

#[test]
fn constant_heating_on_a_latent_plateau_stores_every_joule_once() {
    let law = Law {
        latent: true,
        datum: 1_000_000.0,
        ..Law::constant(1.0, 1.0)
    };
    let (mut thermal, observer) = network(&law, &[350.0], &[2.0], &[], None, vec![]);
    thermal.cells.get_mut(&id(0)).unwrap().boundary.source_w = 40.0;
    let solver = chain(&thermal, &observer, &law);
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    let (next, report) = solver
        .advance_sdirk2_to(&thermal, 2.0, controls(2.0), &mut workspace, || Ok(()))
        .unwrap();
    near(
        next.cells()[&id(0)].specific_h_j_kg,
        law.datum + 390.0,
        1e-9,
    );
    near(next.temperature_k(id(0)).unwrap(), 300.0, 0.0);
    near(report.energy.external_heat_j, 80.0, 1e-12);
    near(report.energy.stored_energy_change_j, 80.0, 1e-9);
    near(
        report.entropy.unclassified_source_at_cell_j_k,
        80.0 / 300.0,
        1e-12,
    );
    near(
        report.entropy.thermal_storage_change_j_k,
        80.0 / 300.0,
        1e-12,
    );
    near(report.entropy.temporal_defect_j_k, 0.0, 1e-12);
    assert!(report.entropy.has_unclassified_sources);
}

#[test]
fn scheduled_boundary_changes_refuse_explicitly() {
    let law = Law::constant(1.0, 1.0);
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    let (thermal, observer) = network(
        &law,
        &[400.0],
        &[1.0],
        &[],
        Some(300.0),
        vec![TimedEvent {
            time_s: 0.5,
            event: Event::Boundary(id(0), boundary(Some(350.0))),
        }],
    );
    let solver = chain(&thermal, &observer, &law);
    assert!(matches!(
        solver.advance_sdirk2_to(&thermal, 1.0, controls(1.0), &mut workspace, || Ok(())),
        Err(ImplicitConductionError::Invalid(_))
    ));
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(
        thermal.cells()[&id(0)]
            .boundary
            .reservoir
            .unwrap()
            .temperature_k,
        300.0
    );
}

#[test]
fn attached_convection_uses_both_stages_without_conductive_entropy_relabeling() {
    let law = Law::constant(1.0, 1.0);
    let (thermal, observer) = network(
        &law,
        &[400.0, 300.0],
        &[1.0, 1.0],
        &[(0, 1, 1.0)],
        None,
        vec![],
    );
    let flow = Advected::new(
        (),
        vec![
            MassFlow {
                id: 0,
                a: id(0),
                b: id(1),
                kg_s: 1.0,
            },
            MassFlow {
                id: 1,
                a: id(1),
                b: id(0),
                kg_s: 1.0,
            },
        ],
        2,
    )
    .unwrap();
    let direct = chain(&thermal, &observer, &law)
        .with_implicit_advection(&thermal, &flow, || Ok(()))
        .unwrap();
    let graph = ImplicitConductionGraph::new(
        &thermal,
        &observer,
        &[(MaterialId(0), &law)],
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap()
    .with_implicit_advection(&thermal, &flow, || Ok(()))
    .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
    // The balanced upwind loop doubles the contrast decay rate to 4/s, but its
    // signed Q/T remains transport, separate from Fourier entropy production.
    let first = [
        350.0 + 50.0 / (1.0 + 4.0 * GAMMA),
        350.0 - 50.0 / (1.0 + 4.0 * GAMMA),
    ];
    let expected = [
        350.0 + 50.0 * amplification(4.0),
        350.0 - 50.0 * amplification(4.0),
    ];
    let production = |t: [f64; 2]| (t[0] - t[1]).powi(2) / (t[0] * t[1]);
    let weighted = (1.0 - GAMMA) * production(first) + GAMMA * production(expected);
    for sparse in [false, true] {
        let (next, report) = if sparse {
            graph.advance_sdirk2_to(&thermal, 1.0, controls(1.0), &mut workspace, || Ok(()))
        } else {
            direct.advance_sdirk2_to(&thermal, 1.0, controls(1.0), &mut workspace, || Ok(()))
        }
        .unwrap();
        for (i, value) in expected.into_iter().enumerate() {
            near(next.temperature_k(id(i)).unwrap(), value, 1e-9);
        }
        near(report.energy.external_heat_j, 0.0, 0.0);
        near(report.energy.residual_j, 0.0, 1e-9);
        near(report.entropy.reservoir_entropy_input_j_k, 0.0, 0.0);
        near(report.entropy.conduction_generation_j_k, weighted, 1e-12);
        near(
            report.entropy.unclassified_source_at_cell_j_k,
            weighted,
            1e-12,
        );
        assert!(report.entropy.has_unclassified_sources);
        assert_eq!(report.energy.substeps, 1);
    }
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(flow.audit().substeps, 0);
}

#[test]
fn open_convection_weights_stream_source_and_bath_energy_in_original_coordinates() {
    for datum in [0.0, 1e9] {
        let law = Law {
            datum,
            ..Law::constant(3.0, 1.0)
        };
        let (mut thermal, observer) = network(&law, &[600.0], &[2.0], &[], Some(300.0), vec![]);
        thermal.cells.get_mut(&id(0)).unwrap().boundary.source_w = 75.0;
        let flow = throughflow(1, 0.7, datum + law.cp * 350.0);
        let solver = chain(&thermal, &observer, &law)
            .with_implicit_advection(&thermal, &flow, || Ok(()))
            .unwrap();
        let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
        let dt = 1.25;
        // 6*T' = 2.1*(350-T) + (300-T) + 75. Neither the absolute
        // enthalpy datum nor an internal reference subtraction is physical heat.
        let equilibrium = (2.1 * 350.0 + 300.0 + 75.0) / 3.1;
        let first = equilibrium + (600.0 - equilibrium) / (1.0 + GAMMA * dt * 3.1 / 6.0);
        let expected = equilibrium + (600.0 - equilibrium) * amplification(dt * 3.1 / 6.0);
        let stream_and_source = |t: f64| 2.1 * (350.0 - t) + 75.0;
        let weighted = |a: f64, b: f64| dt * ((1.0 - GAMMA) * a + GAMMA * b);
        let external = weighted(
            stream_and_source(first) + 300.0 - first,
            stream_and_source(expected) + 300.0 - expected,
        );
        let (next, report) = solver
            .advance_sdirk2_to(
                &thermal,
                dt,
                ImplicitConductionControls {
                    cell_tolerance_j: 2e-6,
                    total_tolerance_j: 5e-6,
                    ..controls(dt)
                },
                &mut workspace,
                || Ok(()),
            )
            .unwrap();
        near(next.temperature_k(id(0)).unwrap(), expected, 1e-6);
        near(report.energy.external_heat_j, external, 5e-6);
        near(
            report.energy.stored_energy_change_j,
            6.0 * (expected - 600.0),
            5e-6,
        );
        near(
            report.entropy.unclassified_source_at_cell_j_k,
            weighted(
                stream_and_source(first) / first,
                stream_and_source(expected) / expected,
            ),
            1e-8,
        );
        near(
            report.entropy.reservoir_entropy_input_j_k,
            weighted((300.0 - first) / 300.0, (300.0 - expected) / 300.0),
            1e-8,
        );
        near(
            report.entropy.conduction_generation_j_k,
            weighted(
                (first - 300.0).powi(2) / (300.0 * first),
                (expected - 300.0).powi(2) / (300.0 * expected),
            ),
            1e-8,
        );
        near(
            report.entropy.thermal_storage_change_j_k,
            6.0 * (expected / 600.0).ln(),
            1e-8,
        );
        assert!(report.entropy.has_unclassified_sources);
        assert_eq!(report.energy.substeps, 1);
        assert_eq!(
            next.cells()[&id(0)].boundary,
            thermal.cells()[&id(0)].boundary
        );
        assert_eq!(flow.audit().substeps, 0);
    }
}

#[test]
fn limited_latent_convection_reconstructs_each_stage_and_refines_time_error() {
    let law = Law {
        latent: true,
        ..Law::constant(1.0, 1.0)
    };
    let (thermal, observer) = network(&law, &[320.0, 360.0, 440.0], &[1.0; 3], &[], None, vec![]);
    let flow = throughflow(3, 1.0, 300.0);
    let reconstruction = LimitedEnthalpyReconstruction::try_new(
        &thermal,
        &[AxialCellStencil {
            previous: id(0),
            cell: id(1),
            next: id(2),
        }],
        1,
        || Ok(()),
    )
    .unwrap();
    let graph = ImplicitConductionGraph::new(
        &thermal,
        &observer,
        &[(MaterialId(0), &law)],
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap()
    .with_implicit_advection_reconstruction(&thermal, &flow, &reconstruction, || Ok(()))
    .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(3).unwrap();
    // In excess-enthalpy coordinates x=h-300, the left minmod slope is active
    // throughout: x' = [-x0, 1.5*(x0-x1), 1.5*x1-0.5*x0-x2]. These scalar
    // stage solves are independent of the native transport residual/Jacobian.
    let stage = |rhs: [f64; 3], weight: f64| {
        let a = rhs[0] / (1.0 + weight);
        let b = (rhs[1] + 1.5 * weight * a) / (1.0 + 1.5 * weight);
        let c = (rhs[2] + weight * (1.5 * b - 0.5 * a)) / (1.0 + weight);
        [a, b, c]
    };
    let initial = [20.0, 60.0, 140.0];
    let first = stage(initial, GAMMA);
    let rates = [
        -first[0],
        1.5 * (first[0] - first[1]),
        1.5 * first[1] - 0.5 * first[0] - first[2],
    ];
    let expected = stage(
        core::array::from_fn(|i| initial[i] + (1.0 - GAMMA) * rates[i]),
        GAMMA,
    );
    for x in [first, expected] {
        assert!(0.0 < x[1] - x[0] && x[1] - x[0] < x[2] - x[1]);
    }
    let (next, report) = graph
        .advance_sdirk2_to(&thermal, 1.0, controls(1.0), &mut workspace, || Ok(()))
        .unwrap();
    for (i, value) in expected.into_iter().enumerate() {
        near(next.cells()[&id(i)].specific_h_j_kg, 300.0 + value, 1e-9);
        near(next.temperature_k(id(i)).unwrap(), 300.0, 0.0);
    }
    let external = -((1.0 - GAMMA) * first[2] + GAMMA * expected[2]);
    near(report.energy.external_heat_j, external, 1e-9);
    near(report.energy.stored_energy_change_j, external, 1e-9);
    near(report.entropy.conduction_generation_j_k, 0.0, 0.0);
    near(
        report.entropy.unclassified_source_at_cell_j_k,
        external / 300.0,
        1e-12,
    );
    near(report.entropy.temporal_defect_j_k, 0.0, 1e-12);
    assert_eq!(report.energy.substeps, 1);

    // Exact semi-discrete solution at t=1: all cells remain on the latent
    // plateau, while h and its reconstructed face trace evolve continuously.
    let exact = [20.0, 60.0, 220.0].map(|x| 300.0 + x * (-1.0_f64).exp());
    let error = |state: &DepositionSimulation<'_>| {
        (0..3)
            .map(|i| (state.cells()[&id(i)].specific_h_j_kg - exact[i]).abs())
            .fold(0.0_f64, f64::max)
    };
    let mut errors = Vec::new();
    for steps in [4, 8, 16] {
        let (next, report) = graph
            .advance_sdirk2_to(
                &thermal,
                1.0,
                controls(1.0 / steps as f64),
                &mut workspace,
                || Ok(()),
            )
            .unwrap();
        errors.push(error(&next));
        near(report.energy.residual_j, 0.0, 1e-8);
        assert_eq!(report.energy.substeps, steps);
    }
    for pair in errors.windows(2) {
        assert!(
            (3.7..4.5).contains(&(pair[0] / pair[1])),
            "limited time errors: {errors:?}"
        );
    }
    let mut adaptive_errors = Vec::new();
    let mut comparisons = Vec::new();
    for tolerance in [1e-3, 1e-5] {
        let (next, report) = graph
            .advance_sdirk2_adaptive_to(
                &thermal,
                1.0,
                AdaptiveImplicitConductionControls {
                    step: controls(1.0),
                    initial_dt_s: 1.0,
                    enthalpy_absolute_tolerance_j_kg: tolerance,
                    temperature_absolute_tolerance_k: 1e-8,
                    relative_tolerance: 0.0,
                    ..Default::default()
                },
                &mut workspace,
                || Ok(()),
            )
            .unwrap();
        adaptive_errors.push(error(&next));
        comparisons.push(report.accepted_comparisons);
        assert!(report.accuracy_rejections > 0);
        assert!(report.maximum_accepted_normalized_error <= 1.0);
        near(
            report.last_error.unwrap().maximum_temperature_difference_k,
            0.0,
            0.0,
        );
        near(
            report.solve.energy.external_heat_j,
            (0..3)
                .map(|i| {
                    next.cells()[&id(i)].specific_h_j_kg - thermal.cells()[&id(i)].specific_h_j_kg
                })
                .sum(),
            1e-8,
        );
        assert_eq!(
            report.solve.energy.substeps,
            2 * report.accepted_comparisons
        );
    }
    assert!(
        adaptive_errors[1] < 0.1 * adaptive_errors[0],
        "{adaptive_errors:?}"
    );
    assert!(comparisons[1] > comparisons[0]);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(flow.audit().substeps, 0);
}

#[test]
fn convecting_source_events_keep_stage_timing_and_late_cancellation_atomic() {
    let law = Law::constant(1.0, 1.0);
    let source = |time_s, source_w| TimedEvent {
        time_s,
        event: Event::Boundary(
            id(0),
            Boundary {
                source_w,
                ..boundary(None)
            },
        ),
    };
    let (thermal, observer) = network(
        &law,
        &[400.0],
        &[1.0],
        &[],
        None,
        vec![source(0.35, 20.0), source(0.8, -10.0), source(1.2, 50.0)],
    );
    let before = thermal.cells().clone();
    let flow = throughflow(1, 1.0, 300.0);
    let solver = chain(&thermal, &observer, &law)
        .with_implicit_advection(&thermal, &flow, || Ok(()))
        .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    let mut expected = 400.0;
    let mut exact = 400.0;
    for (span, equilibrium) in [(0.35, 300.0), (0.45, 320.0), (0.4, 290.0)] {
        expected = equilibrium + (expected - equilibrium) * amplification(span);
        exact = equilibrium + (exact - equilibrium) * (-span).exp();
    }
    let mut checks = 0;
    let (next, report) = solver
        .advance_sdirk2_to(&thermal, 1.2, controls(1.0), &mut workspace, || {
            checks += 1;
            Ok(())
        })
        .unwrap();
    near(next.temperature_k(id(0)).unwrap(), expected, 1e-9);
    near(report.energy.external_heat_j, expected - 400.0, 1e-9);
    assert_eq!(report.energy.substeps, 3);
    assert_eq!(report.energy.events, 3);
    assert_eq!(next.cells()[&id(0)].boundary.source_w, 50.0);
    let mut cancel_checks = 0;
    assert!(matches!(
        solver.advance_sdirk2_to(&thermal, 1.2, controls(1.0), &mut workspace, || {
            cancel_checks += 1;
            if cancel_checks == checks - 1 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(ImplicitConductionError::Thermal(DepositionError::Cancelled))
    ));

    let (next, report) = solver
        .advance_sdirk2_adaptive_to(
            &thermal,
            1.2,
            AdaptiveImplicitConductionControls {
                step: controls(1.0),
                initial_dt_s: 1.0,
                enthalpy_absolute_tolerance_j_kg: 1e-5,
                temperature_absolute_tolerance_k: 1e-5,
                relative_tolerance: 0.0,
                ..Default::default()
            },
            &mut workspace,
            || Ok(()),
        )
        .unwrap();
    near(next.temperature_k(id(0)).unwrap(), exact, 1e-3);
    near(
        report.solve.energy.external_heat_j,
        next.temperature_k(id(0)).unwrap() - 400.0,
        1e-8,
    );
    assert_eq!(report.solve.energy.events, 3);
    assert_eq!(next.cells()[&id(0)].boundary.source_w, 50.0);
    assert_eq!(thermal.cells(), &before);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(flow.audit().substeps, 0);

    // A later bath change is still a whole-window refusal even with admitted
    // flow. It must win over the deliberately insufficient first-stage budget.
    let (thermal, observer) = network(
        &law,
        &[400.0],
        &[1.0],
        &[],
        None,
        vec![TimedEvent {
            time_s: 0.9,
            event: Event::Boundary(
                id(0),
                Boundary {
                    ambient_k: 250.0,
                    ..boundary(None)
                },
            ),
        }],
    );
    let solver = chain(&thermal, &observer, &law)
        .with_implicit_advection(&thermal, &flow, || Ok(()))
        .unwrap();
    let step = ImplicitConductionControls {
        maximum_residual_evaluations: 1,
        ..controls(0.1)
    };
    assert!(matches!(
        solver.advance_sdirk2_to(&thermal, 1.0, step, &mut workspace, || Ok(())),
        Err(ImplicitConductionError::Invalid(_))
    ));
    assert!(matches!(
        solver.advance_sdirk2_adaptive_to(
            &thermal,
            1.0,
            AdaptiveImplicitConductionControls {
                step,
                initial_dt_s: 0.1,
                ..Default::default()
            },
            &mut workspace,
            || Ok(())
        ),
        Err(AdaptiveImplicitConductionError::Implicit(
            ImplicitConductionError::Invalid(_)
        ))
    ));
}
