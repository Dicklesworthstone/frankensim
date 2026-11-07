//! G1/G3/G4: heated forming, release and retained inertia under rigid tooling.
//! Diagonal logarithmic J2 return, original P1 mass and heat accounts below are
//! independent scalar formulas. These material values are synthetic.

use super::*;
use fs_solid::hyper3::plastic::thermal::contact::{
    PlaneContactControls, RigidPlaneTool, TetPlaneContact,
};

fn fixture() -> Fixture {
    let mut law = model();
    for name in [
        "shear_modulus_slope_pa_k",
        "bulk_modulus_slope_pa_k",
        "logarithmic_expansion_per_k",
        "logarithmic_expansion_slope_per_k2",
    ] {
        law.parameters.get_mut(name).unwrap().value = 0.0;
    }
    Fixture::with_model(CP, 301, scalar(DENSITY_PROPERTY, Density::DIMS, RHO), law)
}

#[derive(Clone, Copy)]
struct Point {
    time: f64,
    temperature: f64,
    axial: f64,
    transverse: f64,
    plastic: f64,
    alpha: f64,
    free: f64,
    entropy: f64,
    displacement: [[f64; 3]; 4],
    velocity: [[f64; 3]; 4],
    reaction: f64,
    physical_heat: f64,
    reversible_heat: f64,
    boundary_heat: f64,
    defect: f64,
    kinetic_loss: f64,
}

impl Point {
    fn virgin() -> Self {
        Self {
            time: 0.0,
            temperature: 300.0,
            axial: 1.0,
            transverse: 1.0,
            plastic: 0.0,
            alpha: 0.0,
            free: 0.0,
            entropy: 0.0,
            displacement: [[0.0; 3]; 4],
            velocity: [[0.0; 3]; 4],
            reaction: 0.0,
            physical_heat: 0.0,
            reversible_heat: 0.0,
            boundary_heat: 0.0,
            defect: 0.0,
            kinetic_loss: 0.0,
        }
    }

    fn at(old: Self, axial: f64, temperature: f64, time: f64, inertia: bool) -> Self {
        let transverse: f64 = 1.05;
        let mu = 80.0e6;
        let bulk = 160.0e6;
        let yield_stress = 2.0e6 - 2000.0 * (temperature - 300.0);
        let hardening = 8.0e6 - 5000.0 * (temperature - 300.0);
        let trial = 2.0 * mu * ((axial / transverse).ln() - 1.5 * old.plastic);
        let da = ((trial.abs() - yield_stress - hardening * old.alpha) / (3.0 * mu + hardening))
            .max(0.0);
        let alpha = old.alpha + da;
        let plastic = old.plastic + trial.signum() * da;
        let q = trial - 3.0 * mu * trial.signum() * da;
        let log_volume = (axial * transverse * transverse).ln();
        let piola = [
            (bulk * log_volume + 2.0 * q / 3.0) / axial,
            (bulk * log_volume - q / 3.0) / transverse,
        ];
        let free = q * q / (6.0 * mu)
            + 0.5 * bulk * log_volume * log_volume
            + 0.5 * hardening * alpha * alpha;
        let entropy = 2500.0 * alpha * alpha;
        let coefficient = -2500.0 * (temperature - old.temperature) * old.alpha * old.alpha;
        let physical_heat = yield_stress * da;
        let reversible_heat =
            -coefficient - (temperature * entropy - old.temperature * old.entropy);
        let stretch = [axial, transverse, transverse];
        let displacement = core::array::from_fn(|node| {
            core::array::from_fn(|axis| (stretch[axis] - 1.0) * NODES[node][axis])
        });
        let dt = time - old.time;
        let velocity = core::array::from_fn(|node| {
            core::array::from_fn(|axis| {
                (displacement[node][axis] - old.displacement[node][axis]) / dt
            })
        });
        let delta_v = core::array::from_fn(|node| {
            core::array::from_fn(|axis| velocity[node][axis] - old.velocity[node][axis])
        });
        let mass_delta = mass_action(&delta_v);
        let inertial_force = if inertia { mass_delta[1][0] / dt } else { 0.0 };
        let work = piola[0] * (axial - old.axial) + 2.0 * piola[1] * (transverse - old.transverse);
        Self {
            time,
            temperature,
            axial,
            transverse,
            plastic,
            alpha,
            free,
            entropy,
            displacement,
            velocity,
            reaction: -V / L * piola[0] - inertial_force,
            physical_heat,
            reversible_heat,
            boundary_heat: RHO * V * CP * (temperature - old.temperature)
                - V * (physical_heat + reversible_heat),
            defect: V * (work - (free - old.free - coefficient) - physical_heat),
            kinetic_loss: if inertia {
                0.5 * dot(&delta_v, &mass_delta)
            } else {
                0.0
            },
        }
    }

    fn supports(self) -> Vec<(usize, f64)> {
        self.displacement
            .as_flattened()
            .iter()
            .enumerate()
            .filter(|(dof, _)| *dof != 3)
            .map(|(dof, value)| (dof, *value))
            .collect()
    }
}

fn path(inertia: bool) -> Vec<Point> {
    let mut old = Point::virgin();
    let mut points = Vec::new();
    for (i, (temperature, dt)) in [(310.0, 0.00011), (380.0, 0.00017), (320.0, 0.00019)]
        .into_iter()
        .enumerate()
    {
        let time = old.time + dt;
        let axial = if i < 2 {
            0.85
        } else {
            let (mut lo, mut hi) = (0.6, 1.3);
            assert!(Point::at(old, lo, temperature, time, inertia).reaction > 0.0);
            assert!(Point::at(old, hi, temperature, time, inertia).reaction < 0.0);
            for _ in 0..80 {
                let mid = 0.5 * (lo + hi);
                if Point::at(old, mid, temperature, time, inertia).reaction > 0.0 {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            0.5 * (lo + hi)
        };
        let next = Point::at(old, axial, temperature, time, inertia);
        assert!(next.alpha > old.alpha);
        points.push(next);
        old = next;
    }
    points
}

fn thermal<'a>(fixture: &'a Fixture, points: &[Point]) -> DepositionSimulation<'a> {
    let mut events = vec![TimedEvent {
        time_s: 0.0,
        event: Event::Birth(DepositedCell {
            id: CELL,
            material_id: MaterialId(0),
            segment_id: 1,
            body_id: 7,
            mass_kg: RHO * V,
            specific_h_j_kg: 0.0,
            geometry: CellGeometry {
                start_m: [0.0; 3],
                end_m: [L; 3],
                volume_m3: V,
                surface_area_m2: 3.0 * L * L,
            },
            boundary: boundary(points[0].boundary_heat / points[0].time),
        }),
    }];
    for (old, next) in points.iter().zip(points.iter().skip(1)) {
        events.push(TimedEvent {
            time_s: old.time,
            event: Event::Boundary(CELL, boundary(next.boundary_heat / (next.time - old.time))),
        });
    }
    DepositionSimulation::new(
        events,
        &[(MaterialId(0), fixture.caloric.curve() as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 8,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(0.0, heat_config(), || Ok(()))
    .unwrap()
    .0
}

fn heat_config() -> StepConfig {
    StepConfig {
        max_dt_s: 1e-5,
        energy_tolerance_j: 1e-10,
    }
}

fn contact_controls() -> PlaneContactControls {
    PlaneContactControls {
        maximum_normal_work_defect_j: 1e-7,
        work_balance_tolerance_j: 1e-8,
        ..PlaneContactControls::default()
    }
}

fn tool(index: usize) -> RigidPlaneTool {
    let before = if index == 0 { L } else { 0.85 * L };
    let after = if index < 2 { 0.85 * L } else { 1.1 * L };
    RigidPlaneTool {
        normal: [-1.0, 0.0, 0.0],
        point_before_m: [before, 0.0, 0.0],
        translation_m: [after - before, 0.0, 0.0],
    }
}

fn contacts() -> [TetPlaneContact; 1] {
    [TetPlaneContact {
        node: 1,
        plane: 0,
        initial_force_n: 0.0,
        scale_n_per_m: 1e6,
    }]
}

#[test]
fn heated_tool_forming_and_release_retain_permanent_metric_and_close_both_work_ports() {
    context(|cx| {
        let fixture = fixture();
        let points = path(false);
        let mut thermal = thermal(&fixture, &points);
        let mut solid = bind(&fixture, &thermal, false, cx);
        for (i, expected) in points.iter().enumerate() {
            let next = solid
                .advance_with_plane_contact(
                    &thermal,
                    &(),
                    expected.time,
                    &expected.supports(),
                    &[],
                    [0.0; 3],
                    &[tool(i)],
                    &contacts(),
                    settings(),
                    controls(),
                    contact_controls(),
                    heat_config(),
                    config(),
                    cx,
                )
                .unwrap();
            let mechanics = next.solid().mechanics();
            let report = next.report();
            let energy = &report.material_energy;
            close(
                next.thermal().temperature_k(CELL).unwrap(),
                expected.temperature,
                2e-7,
            );
            assert_eq!(
                mechanics.temperatures_k()[0],
                next.thermal().temperature_k(CELL).unwrap()
            );
            close(mechanics.element_deformation()[0][0], expected.axial, 2e-9);
            let state = &mechanics.element_states()[0];
            close(state.alpha, expected.alpha, 3e-9);
            close(
                state.inverse_plastic_cauchy_green[0],
                (-2.0 * expected.plastic).exp(),
                5e-9,
            );
            close(
                state.inverse_plastic_cauchy_green[1],
                expected.plastic.exp(),
                5e-9,
            );
            close(
                energy.contact.solve.forces_n[0][0],
                expected.reaction.max(0.0),
                2e-5,
            );
            assert!(energy.contact.minimum_gap_m >= 0.0);
            close(energy.contact.dead_load_work_j, 0.0, 1e-10);
            close(
                energy.material.dead_load_endpoint_work_j,
                energy.contact.body_work_j,
                1e-10,
            );
            close(
                energy.contact.tool_work_j,
                -energy.contact.solve.forces_n[0][0] * tool(i).translation_m[0],
                1e-10,
            );
            close(
                energy.external_work_defect_j,
                energy.material.endpoint_work_defect_j + energy.contact.normal_work_defect_j,
                1e-12,
            );
            close(
                energy.material.endpoint_work_defect_j,
                expected.defect,
                2e-8,
            );
            close(
                energy.material.plastic_dissipation_j,
                V * expected.physical_heat,
                2e-9,
            );
            close(
                energy.material.caloric_heat_j,
                V * (expected.physical_heat + expected.reversible_heat),
                2e-9,
            );
            close(report.external_heat_j, expected.boundary_heat, 2e-9);
            close(report.combined_residual_j, 0.0, config().energy_tolerance_j);
            close(next.thermal().cells()[&CELL].mass_kg, RHO * V, 0.0);
            assert!(std::ptr::eq(
                mechanics.families()[0].caloric_source(),
                &fixture.caloric
            ));
            assert_eq!(
                mechanics.families()[0].query().card.content_hash(),
                fixture.card.content_hash()
            );
            if i == 1 {
                // Genuine hot re-yield at fixed total deformation, not a reset.
                assert!(state.alpha > points[0].alpha);
                close(energy.material.internal_endpoint_work_j, 0.0, 1e-8);
            }
            if i == 2 {
                assert!(mechanics.element_deformation()[0][0] < 0.9);
                close(energy.contact.solve.forces_n[0][0], 0.0, 1e-8);
                assert!(energy.contact.solve.paths_m[0][0] > 0.002);
            }
            (thermal, solid, (), _) = next.into_parts();
        }
    });
}

#[test]
fn inertial_contact_uses_actual_interval_and_late_heat_failure_preserves_accepted_forming() {
    context(|cx| {
        let fixture = fixture();
        let points = path(true);
        let mut thermal = thermal(&fixture, &points);
        let mut solid = bind(&fixture, &thermal, false, cx)
            .with_inertia(&[[0.0; 3]; 4], cx)
            .unwrap();
        for (i, expected) in points.iter().enumerate() {
            if i == 1 {
                let cells = thermal.cells().clone();
                let history = solid.mechanics().element_states().to_vec();
                let displacement = solid.mechanics().displacement_m().to_vec();
                let velocity = solid.mechanics().velocity_mps().unwrap().to_vec();
                let observer = LateRefusal {
                    completed_substeps: 0,
                };
                let refused = solid.advance_with_plane_contact(
                    &thermal,
                    &observer,
                    expected.time,
                    &expected.supports(),
                    &[],
                    [0.0; 3],
                    &[tool(i)],
                    &contacts(),
                    settings(),
                    controls(),
                    contact_controls(),
                    heat_config(),
                    config(),
                    cx,
                );
                assert!(matches!(
                    refused,
                    Err(FiniteThermoplasticError::Observer(
                        DepositionError::Cancelled
                    ))
                ));
                assert_eq!(observer.completed_substeps, 0);
                assert_eq!(thermal.cells(), &cells);
                assert_eq!(solid.mechanics().element_states(), history);
                assert_eq!(solid.mechanics().displacement_m(), displacement);
                assert_eq!(solid.mechanics().velocity_mps().unwrap(), velocity);
                assert_eq!(solid.time_s(), points[0].time);
            }
            let next = solid
                .advance_with_plane_contact_and_coupling(
                    &thermal,
                    &(),
                    expected.time,
                    &expected.supports(),
                    &[],
                    [0.0; 3],
                    &[tool(i)],
                    &contacts(),
                    settings(),
                    controls(),
                    contact_controls(),
                    heat_config(),
                    config(),
                    CouplingMethod::IqnIls(IqnIlsConfig::default()),
                    cx,
                )
                .unwrap();
            let report = next.report();
            close(
                next.solid().mechanics().element_deformation()[0][0],
                expected.axial,
                2e-9,
            );
            close(
                report.material_energy.contact.solve.forces_n[0][0],
                expected.reaction.max(0.0),
                2e-5,
            );
            close(
                report.inertia.as_ref().unwrap().numerical_loss_j,
                expected.kinetic_loss,
                2e-10,
            );
            assert!(report.inertia.as_ref().unwrap().numerical_loss_j > 1e-5);
            // Heat excludes both contact quadrature defect and numerical inertia.
            close(
                report.material_energy.material.caloric_heat_j,
                V * (expected.physical_heat + expected.reversible_heat),
                2e-9,
            );
            close(
                next.thermal().temperature_k(CELL).unwrap(),
                expected.temperature,
                2e-7,
            );
            close(report.combined_residual_j, 0.0, config().energy_tolerance_j);
            (thermal, solid, (), _) = next.into_parts();
        }
    });
}

#[test]
fn moving_tool_and_current_surface_conduction_share_one_actual_thermal_endpoint() {
    use fs_thermal::enthalpy::conduction::{
        ConductionObserver, ConductiveReservoir, ConductivityLaw, CurrentConductionGeometry,
        CurrentConductiveBoundary,
    };

    struct ConstantConductivity;
    impl ConductivityLaw for ConstantConductivity {
        fn conductivity_w_m_k(&self, h: f64) -> Result<f64, DepositionError> {
            if h.is_finite() && (-20.0 * CP..=120.0 * CP).contains(&h) {
                Ok(0.5)
            } else {
                Err(DepositionError::MaterialDomain)
            }
        }
    }

    context(|cx| {
        let fixture = fixture();
        let expected = Point::at(Point::virgin(), 0.85, 310.0, 1.0, false);
        let area = 0.5 * L * L;
        let distance = 0.25 * L;
        let bath = 340.0;
        let g = 0.5 * area * expected.transverse.powi(2) / (distance * expected.axial);
        let r4 = (1.0 - 0.25 * g / (RHO * V * CP)).powi(4);
        let total_power = g * (expected.temperature - bath - r4 * (300.0 - bath)) / (1.0 - r4);
        let mut cell = thermal(&fixture, &[expected]).cells()[&CELL].clone();
        cell.boundary = Boundary {
            ambient_k: 300.0,
            htc_w_m2_k: 0.0,
            reservoir: Some(Reservoir {
                area_m2: area,
                conductance_w_k: 0.0,
                temperature_k: bath,
            }),
            source_w: total_power - V * (expected.physical_heat + expected.reversible_heat),
        };
        let thermal = DepositionSimulation::new(
            vec![TimedEvent {
                time_s: 0.0,
                event: Event::Birth(cell),
            }],
            &[(MaterialId(0), fixture.caloric.curve() as &dyn ThermalLaw)],
            DepositionBudget {
                max_cells: 1,
                max_links: 0,
                max_events: 2,
                max_substeps: 100,
            },
        )
        .unwrap()
        .advance_to(0.0, thermal_config(), || Ok(()))
        .unwrap()
        .0;
        let solid = bind(&fixture, &thermal, false, cx);
        let conductivity = ConstantConductivity;
        let conduction = ConductionObserver::new(
            &thermal,
            &[(MaterialId(0), &conductivity)],
            &[],
            &[ConductiveReservoir {
                cell_id: CELL,
                half_distance_m: distance,
                surface_resistance_m2_k_w: 0.0,
            }],
            || Ok(()),
        )
        .unwrap();
        let (next, geometry_defect) = solid
            .advance_with_plane_contact_and_deformed_conduction(
                &thermal,
                &conduction,
                expected.time,
                &expected.supports(),
                &[],
                [0.0; 3],
                &[tool(0)],
                &contacts(),
                settings(),
                controls(),
                contact_controls(),
                thermal_config(),
                config(),
                CouplingMethod::RelaxedPicard,
                1e-10,
                |mechanics, cx| {
                    cx.checkpoint().map_err(|_| DepositionError::Cancelled)?;
                    let f = mechanics.element_deformation()[0];
                    Ok(CurrentConductionGeometry {
                        links: BTreeMap::new(),
                        reservoirs: [(
                            CELL,
                            CurrentConductiveBoundary {
                                area_m2: area * f[4] * f[8],
                                half_distance_m: distance * f[0],
                            },
                        )]
                        .into(),
                        ambient_faces: BTreeMap::new(),
                    })
                },
                cx,
            )
            .unwrap();
        let actual = next.solid().mechanics();
        let report = next.report();
        close(actual.element_deformation()[0][0], expected.axial, 2e-9);
        close(actual.element_states()[0].alpha, expected.alpha, 3e-9);
        close(
            next.thermal().temperature_k(CELL).unwrap(),
            expected.temperature,
            2e-7,
        );
        assert_eq!(
            actual.temperatures_k()[0],
            next.thermal().temperature_k(CELL).unwrap()
        );
        close(
            report.material_energy.contact.solve.forces_n[0][0],
            expected.reaction,
            2e-5,
        );
        close(
            report.material_energy.material.caloric_heat_j,
            V * (expected.physical_heat + expected.reversible_heat),
            2e-9,
        );
        close(report.external_heat_j, expected.boundary_heat, 2e-9);
        close(report.combined_residual_j, 0.0, config().energy_tolerance_j);
        assert!(geometry_defect <= 1e-10);
        assert_eq!(report.thermal_audit.substeps, 4);
        let metric = next.observer().current_geometry().unwrap().reservoirs[&CELL];
        close(0.5 * metric.area_m2 / metric.half_distance_m, g, 2e-11);
        assert!((g - 0.5 * area / distance).abs() > 0.002);
        assert_eq!(thermal.time_s(), 0.0);
        assert!(conduction.current_geometry().is_none());
    });
}
