//! G1/G3: original thermal tables and actual-time flow under follower pressure.
//! Coaxial return, consistent inertia, caloric closure and volume work below are
//! independent scalar formulas; all coefficients are synthetic, not calibration.

use super::*;
use fs_material::plastic::finite::thermal::source::tabulated::THERMAL_PERZYNA_FINITE_J2_TABULATED_LAW_ID;
use fs_solid::hyper3::pressure::ClosedTetPressureSurface;

const ETA: f64 = 4000.0;
const ACTIVATION: f64 = 1800.0;

fn pressure_fixture() -> Fixture {
    let original = sourced_fixture(90);
    let mut constitutive = original.card.models()[0].clone();
    constitutive.law = LawId(THERMAL_PERZYNA_FINITE_J2_TABULATED_LAW_ID.into());
    for (name, value, dims) in [
        ("overstress_reference_viscosity_pa_s", ETA, Dims([-1, 1, -1, 0, 0, 0])),
        ("overstress_stress_scale_pa", 2.0e6, PA),
        ("overstress_exponent", 1.0, Dims::NONE),
        ("overstress_activation_temperature_k", ACTIVATION, K),
    ] {
        constitutive.parameters.insert(name.into(), LawParameter { value, dims });
    }
    let mut claims = ClaimSet::new();
    claims.insert_claim(scalar(CP_NAME, CP_DIMS, CP)).unwrap();
    claims.insert_claim(scalar(DENSITY_PROPERTY, Density::DIMS, RHO)).unwrap();
    claims.insert_claim(surface()).unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic heated pressure forming with tabulated resistance".into(),
            phase: "solid".into(),
            process: "numerical regression".into(),
            revision: 0,
        },
        claims,
        vec![constitutive],
    ).unwrap();
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
        -20.0 * CP,
        1e-6,
    ).unwrap();
    Fixture { card, caloric }
}

#[derive(Clone, Copy)]
struct PressurePoint {
    time: f64,
    temperature: f64,
    axial: f64,
    transverse: f64,
    plastic: f64,
    alpha: f64,
    displacement: [[f64; 3]; 4],
    velocity: [[f64; 3]; 4],
    piola: [f64; 3],
    pressure: f64,
    free: f64,
    entropy: f64,
    physical_heat: f64,
    viscous_heat: f64,
    reversible_heat: f64,
    boundary_heat_j: f64,
    defect: f64,
    kinetic_change_j: f64,
    kinetic_loss_j: f64,
}

impl PressurePoint {
    fn virgin() -> Self {
        Self {
            time: 0.0, temperature: 300.0, axial: 1.0, transverse: 1.0,
            plastic: 0.0, alpha: 0.0, displacement: [[0.0; 3]; 4],
            velocity: [[0.0; 3]; 4], piola: [0.0; 3], pressure: 0.0,
            free: 0.0, entropy: 0.0, physical_heat: 0.0, viscous_heat: 0.0,
            reversible_heat: 0.0, boundary_heat_j: 0.0, defect: 0.0,
            kinetic_change_j: 0.0, kinetic_loss_j: 0.0,
        }
    }

    fn at(old: Self, axial: f64, transverse: f64, temperature: f64, time: f64, inertia: bool) -> Self {
        let dt = time - old.time;
        let y = ordinates(temperature);
        let eta = ETA * (ACTIVATION * (1.0 / temperature - 1.0 / 300.0)).exp();
        let trial = 2.0 * MU * ((axial / transverse).ln() - 1.5 * old.plastic);
        let mut increment = 0.0;
        if trial.abs() > interpolate(old.alpha, y) {
            let (mut lower, mut upper) = (0.0, trial.abs() / (3.0 * MU));
            for _ in 0..80 {
                let middle = 0.5 * (lower + upper);
                let resistance = interpolate(old.alpha + middle, y)
                    + (3.0 * MU + eta / dt) * middle;
                if resistance < trial.abs() {
                    lower = middle;
                } else {
                    upper = middle;
                }
            }
            increment = 0.5 * (lower + upper);
        }
        let alpha = old.alpha + increment;
        let plastic = old.plastic + trial.signum() * increment;
        let q = trial - 3.0 * MU * trial.signum() * increment;
        let log_volume = (axial * transverse * transverse).ln();
        let bulk = 160e6;
        let piola = [
            (bulk * log_volume + 2.0 * q / 3.0) / axial,
            (bulk * log_volume - q / 3.0) / transverse,
            (bulk * log_volume - q / 3.0) / transverse,
        ];
        let deformation = [
            axial, 0.0, 0.0, 0.0, transverse, 0.0, 0.0, 0.0, transverse,
        ];
        let displacement = displacement(&deformation);
        let velocity: [[f64; 3]; 4] = core::array::from_fn(|node| {
            core::array::from_fn(|axis| {
                (displacement[node][axis] - old.displacement[node][axis]) / dt
            })
        });
        let delta_v: [[f64; 3]; 4] = core::array::from_fn(|node| {
            core::array::from_fn(|axis| velocity[node][axis] - old.velocity[node][axis])
        });
        let mass_delta_v = super::super::mass_action(&delta_v);
        let inertial_force = if inertia { mass_delta_v[1][0] / dt } else { 0.0 };
        let pressure = -(piola[0] + L * inertial_force / V) / (transverse * transverse);
        let free = q * q / (6.0 * MU)
            + 0.5 * bulk * log_volume * log_volume + primitive(alpha, y);
        let entropy = -primitive(alpha, SLOPES);
        let coefficient = (temperature - old.temperature) * primitive(old.alpha, SLOPES);
        let delta_ts = temperature * entropy - old.temperature * old.entropy;
        let viscous_heat = eta * increment * increment / dt;
        let physical_heat = y[0] * increment + viscous_heat;
        let reversible_heat = -coefficient - delta_ts;
        let work = piola[0] * (axial - old.axial)
            + 2.0 * piola[1] * (transverse - old.transverse);
        Self {
            time, temperature, axial, transverse, plastic, alpha, displacement, velocity,
            piola, pressure, free, entropy, physical_heat, viscous_heat, reversible_heat,
            boundary_heat_j: RHO * V * CP * (temperature - old.temperature)
                - V * (physical_heat + reversible_heat),
            defect: work - (free - old.free - coefficient) - physical_heat,
            kinetic_change_j: 0.5 * (
                super::super::dot(&velocity, &super::super::mass_action(&velocity))
                    - super::super::dot(&old.velocity, &super::super::mass_action(&old.velocity))
            ),
            kinetic_loss_j: 0.5 * super::super::dot(&delta_v, &mass_delta_v),
        }
    }

    fn under_pressure(
        old: Self, transverse: f64, temperature: f64, time: f64,
        pressure: f64, inertia: bool,
    ) -> Self {
        let (mut lower, mut upper) = (0.6, 1.3);
        assert!(Self::at(old, lower, transverse, temperature, time, inertia).pressure > pressure);
        assert!(Self::at(old, upper, transverse, temperature, time, inertia).pressure < pressure);
        for _ in 0..80 {
            let axial = 0.5 * (lower + upper);
            if Self::at(old, axial, transverse, temperature, time, inertia).pressure > pressure {
                lower = axial;
            } else {
                upper = axial;
            }
        }
        let mut point = Self::at(
            old, 0.5 * (lower + upper), transverse, temperature, time, inertia,
        );
        point.pressure = pressure;
        point
    }

    fn prescribed(self) -> Vec<(usize, f64)> {
        self.displacement.as_flattened().iter().enumerate()
            .filter(|(dof, _)| *dof != 3)
            .map(|(dof, value)| (dof, *value)).collect()
    }

    fn jacobian(self) -> f64 {
        self.axial * self.transverse * self.transverse
    }

    fn caloric_heat(self) -> f64 {
        self.physical_heat + self.reversible_heat
    }
}

fn pressure_path(inertia: bool) -> Vec<PressurePoint> {
    let mut old = PressurePoint::virgin();
    let held = PressurePoint::at(old, 0.85, 1.0, 310.0, 0.00011, inertia).pressure;
    let mut points = Vec::new();
    for (dt, transverse, temperature, pressure) in [
        (0.00011, 1.0, 310.0, held),
        (0.00017, 1.07, 400.0, held),
        (0.00019, 1.07, 320.0, 0.0),
    ] {
        let point = PressurePoint::under_pressure(
            old, transverse, temperature, old.time + dt, pressure, inertia,
        );
        assert!(point.alpha > old.alpha && point.viscous_heat > 0.0);
        points.push(point);
        old = point;
    }
    assert!(points[0].alpha > 0.04 && points[0].alpha < 0.12);
    assert!(points[1].alpha > 0.12 && points[1].alpha < 0.3);
    assert!(points[2].axial < 0.9);
    points
}

fn pressure_initial(fixture: &Fixture) -> DepositionSimulation<'_> {
    DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CELL, material_id: MaterialId(0), segment_id: 1, body_id: 7,
                mass_kg: RHO * V, specific_h_j_kg: 0.0,
                geometry: CellGeometry {
                    start_m: [0.0; 3], end_m: [L; 3],
                    volume_m3: V, surface_area_m2: 3.0 * L * L,
                },
                boundary: boundary(0.0),
            }),
        }],
        &[(MaterialId(0), fixture.caloric.curve() as &dyn ThermalLaw)],
        DepositionBudget { max_cells: 1, max_links: 0, max_events: 1, max_substeps: 1000 },
    ).unwrap().advance_to(0.0, thermal_config(), || Ok(())).unwrap().0
}

#[derive(Clone, Default)]
struct PressureHeat {
    watts: f64,
    delivered_j: f64,
    steps: usize,
    fail_after: Option<usize>,
}

impl DepositionObserver for PressureHeat {
    type Error = DepositionError;
    fn sources(
        &mut self, _: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        checkpoint()?;
        Ok([(CELL, self.watts)].into())
    }
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        Ok(())
    }
    fn step(
        &mut self, _: f64, dt: f64, _: &BTreeMap<CellId, f64>, _: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        checkpoint()?;
        self.delivered_j += dt * self.watts;
        self.steps += 1;
        if self.fail_after == Some(self.steps) {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    }
}

fn run_pressure_path(inertia: bool, method: CouplingMethod) {
    context(|cx| {
        let fixture = pressure_fixture();
        let mut thermal = pressure_initial(&fixture);
        let mut solid = bind(&fixture, &thermal, cx);
        if inertia {
            solid = solid.with_inertia(&[[0.0; 3]; 4], cx).unwrap();
        }
        let surface = ClosedTetPressureSurface::new(
            &NODES, &TETS, TetAssemblyBudget::standard(), 16, cx,
        ).unwrap();
        let mut observer = PressureHeat::default();
        let mut old = PressurePoint::virgin();
        for (index, expected) in pressure_path(inertia).into_iter().enumerate() {
            let dt = expected.time - old.time;
            let heat_config = StepConfig { max_dt_s: dt / 3.0, ..thermal_config() };
            observer.watts = expected.boundary_heat_j / dt;
            let supports = expected.prescribed();
            if inertia && index == 1 {
                let cells = thermal.cells().clone();
                let histories = solid.mechanics().element_states().to_vec();
                let displacements = solid.mechanics().displacement_m().to_vec();
                let velocities = solid.mechanics().velocity_mps().unwrap().to_vec();
                let temperatures = solid.mechanics().temperatures_k().to_vec();
                let mut failed_observer = observer.clone();
                failed_observer.fail_after = Some(observer.steps + 2);
                let failed = solid.advance_with_pressure(
                    &thermal, &failed_observer, expected.time, &supports, &[], [0.0; 3],
                    &surface, expected.pressure, settings(), controls(), heat_config, config(), cx,
                );
                assert!(matches!(failed, Err(FiniteThermoplasticError::Observer(DepositionError::Cancelled))));
                let failed = solid.advance_with_pressure(
                    &thermal, &observer, expected.time, &supports, &[], [0.0; 3],
                    &surface, expected.pressure, settings(), controls(), heat_config,
                    FiniteThermoplasticConfig { maximum_iterations: 1, ..config() }, cx,
                );
                assert!(matches!(failed, Err(FiniteThermoplasticError::NotConverged { iterations: 1, .. })));
                // Actual pressure correction opposes the material defect here.
                // Their small signed sum must not conceal two large errors.
                let mut strict = controls();
                strict.material.maximum_endpoint_work_defect_j = 0.08;
                let pressure_work = -expected.pressure * V
                    * (expected.jacobian() - old.jacobian());
                let pressure_endpoint = -expected.pressure * V
                    * (expected.transverse.powi(2) * (expected.axial - old.axial)
                        + 2.0 * expected.axial * expected.transverse
                            * (expected.transverse - old.transverse));
                let correction = pressure_work - pressure_endpoint;
                assert!((V * expected.defect).abs() < 0.08);
                assert!((V * expected.defect + correction).abs() < 0.08);
                assert!((V * expected.defect).abs() + correction.abs() > 0.08);
                let refused = solid.mechanics().advance_with_pressure(
                    &[expected.temperature], &supports, &[], [0.0; 3], dt,
                    &surface, expected.pressure, settings(), strict, cx,
                ).unwrap_err();
                assert!(refused.to_string().contains("pressure/constitutive endpoint"));

                assert_eq!(failed_observer.steps, observer.steps);
                assert_eq!(failed_observer.delivered_j, observer.delivered_j);
                assert_eq!(thermal.cells(), &cells);
                assert_eq!(solid.mechanics().element_states(), histories);
                assert_eq!(solid.mechanics().displacement_m(), displacements);
                assert_eq!(solid.mechanics().velocity_mps().unwrap(), velocities);
                assert_eq!(solid.mechanics().temperatures_k(), temperatures);
                assert_eq!(thermal.time_s(), old.time);
                assert_eq!(solid.time_s(), old.time);
            }
            let candidate = solid.advance_with_pressure_and_coupling(
                &thermal, &observer, expected.time, &supports, &[], [0.0; 3],
                &surface, expected.pressure, settings(), controls(), heat_config,
                config(), method, cx,
            ).unwrap();
            let mechanics = candidate.solid().mechanics();
            let report = candidate.report();
            let energy = &report.material_energy;
            let pressure_work = -expected.pressure * V * (expected.jacobian() - old.jacobian());
            let endpoint_work = -expected.pressure * V
                * (expected.transverse.powi(2) * (expected.axial - old.axial)
                    + 2.0 * expected.axial * expected.transverse * (expected.transverse - old.transverse));
            close(candidate.thermal().temperature_k(CELL).unwrap(), expected.temperature, 2e-6);
            assert_eq!(mechanics.temperatures_k()[0], candidate.thermal().temperature_k(CELL).unwrap());
            close(mechanics.displacement_m()[1][0], expected.displacement[1][0], 2e-10);
            close(mechanics.element_states()[0].alpha, expected.alpha, 3e-8);
            close(mechanics.free_energy_j(), V * expected.free, 3e-8);
            for axis in 0..3 {
                close(mechanics.element_first_piola_pa()[0][4 * axis], expected.piola[axis], 1.0);
                close(
                    mechanics.element_states()[0].inverse_plastic_cauchy_green[axis],
                    (if axis == 0 { -2.0 * expected.plastic } else { expected.plastic }).exp(),
                    4e-8,
                );
            }
            close(energy.material.plastic_dissipation_j, V * expected.physical_heat, 3e-8);
            close(energy.material.reversible_heat_j, V * expected.reversible_heat, 3e-8);
            close(energy.material.caloric_heat_j, V * expected.caloric_heat(), 3e-8);
            close(energy.material.endpoint_work_defect_j, V * expected.defect, 3e-8);
            close(energy.pressure.work_j, pressure_work, 3e-8);
            close(energy.pressure.endpoint_work_j, endpoint_work, 3e-8);
            close(energy.pressure.endpoint_work_correction_j, pressure_work - endpoint_work, 3e-8);
            close(
                energy.pressure.external_work_defect_j,
                V * expected.defect + pressure_work - endpoint_work,
                5e-8,
            );
            close(report.applied_heat_j, V * expected.caloric_heat(), 3e-8);
            close(report.external_heat_j, expected.boundary_heat_j, 3e-8);
            close(
                report.thermal_audit.stored_energy_change_j,
                RHO * V * CP * (expected.temperature - old.temperature),
                3e-8,
            );
            close(candidate.observer().delivered_j, observer.delivered_j + expected.boundary_heat_j, 3e-9);
            close(report.thermal_residual_j, 0.0, config().energy_tolerance_j);
            close(report.combined_residual_j, 0.0, config().energy_tolerance_j);
            if inertia {
                let motion = report.inertia.as_ref().unwrap();
                close(motion.reference_mass_kg, RHO * V, 1e-18);
                close(motion.kinetic_change_j, expected.kinetic_change_j, 3e-9);
                close(motion.numerical_loss_j, expected.kinetic_loss_j, 3e-9);
                assert!(motion.numerical_loss_j > 1e-5);
                close(mechanics.velocity_mps().unwrap()[1][0], expected.velocity[1][0], 2e-6);
            } else {
                assert!(report.inertia.is_none());
            }
            let material = &mechanics.families()[0];
            assert!(std::ptr::eq(material.caloric_source(), &fixture.caloric));
            assert_eq!(material.family().model_card(), &fixture.card.models()[0]);
            let source = material.family().tabulated_hardening_source().unwrap().source();
            assert_eq!(source.card_identity(), fixture.card.content_hash());
            assert_eq!(source.observations().len(), SAMPLE_LIMIT);
            if index == 1 {
                assert!(energy.pressure.endpoint_work_correction_j < -0.05);
            }
            let (next_thermal, next_solid, next_observer, _) = candidate.into_parts();
            thermal = next_thermal;
            solid = next_solid;
            observer = next_observer;
            old = expected;
        }
        assert!(solid.mechanics().displacement_m()[1][0] < -0.001);
    });
}

#[test]
fn source_table_pressure_heating_and_release_recover_both_feedback_methods() {
    run_pressure_path(false, CouplingMethod::RelaxedPicard);
    run_pressure_path(false, CouplingMethod::IqnIls(IqnIlsConfig::default()));
}

#[test]
fn pressure_forming_retains_actual_time_inertia_and_atomic_retry() {
    run_pressure_path(true, CouplingMethod::RelaxedPicard);
}
