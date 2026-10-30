//! G1/G3/G4: actual native traction, accelerated material residence and finite heat.
use fs_couple::friction_heat::moving::trajectory::QuadraticContactPath;
use fs_couple::friction_heat::{FrictionHeatTolerances, moving::*};
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::*;
use fs_tribo::tangential::elliptic::EllipticCoulomb;
use fs_tribo::tangential::elliptic::driven::*;
use fs_tribo::tangential::elliptic::relaxation::EllipticRelaxationControls;
use fs_tribo::{TangentialLaw, TangentialSolveConfig};
use std::cell::Cell;
use std::collections::BTreeMap;
use std::convert::Infallible;

struct Sensible;
impl ThermalLaw for Sensible {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if h.is_finite() && h > -300.0 {
            Ok(300.0 + h)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
    fn specific_entropy_change_j_kg_k(
        &self,
        a: f64,
        b: f64,
    ) -> Result<Option<f64>, DepositionError> {
        let t = self.temperature_k(a)?;
        self.temperature_k(b)?;
        Ok(Some(((b - a) / t).ln_1p()))
    }
}
static SENSIBLE: Sensible = Sensible;
fn close(a: f64, b: f64) {
    assert!((a - b).abs() <= 2e-11, "{a:.16e} != {b:.16e}");
}
fn position(x: f64) -> MovingContactPosition {
    MovingContactPosition {
        body_offsets_m: [x, 0.0],
        footprint_offset_m: x,
    }
}
fn fixture() -> (DepositionSimulation<'static>, MovingFrictionRoute) {
    let faces = [
        vec![StripSurfaceCell {
            cell: CellId(0),
            span_m: [0.0, 1.0],
        }],
        vec![
            StripSurfaceCell {
                cell: CellId(1),
                span_m: [0.0, 1.0],
            },
            StripSurfaceCell {
                cell: CellId(2),
                span_m: [1.0, 2.0],
            },
        ],
    ];
    let mut events = Vec::new();
    for (side, surface) in faces.iter().enumerate() {
        for face in surface {
            let centre = f64::midpoint(face.span_m[0], face.span_m[1]);
            events.push(TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id: face.cell,
                    segment_id: face.cell.0,
                    body_id: side as u64,
                    material_id: MaterialId(0),
                    mass_kg: 1.0,
                    specific_h_j_kg: 0.0,
                    geometry: CellGeometry {
                        start_m: [centre, 0.0, 0.0],
                        end_m: [centre, if side == 0 { -1.0 } else { 1.0 }, 0.0],
                        volume_m3: 1.0,
                        surface_area_m2: 6.0,
                    },
                    boundary: Boundary {
                        ambient_k: 300.0,
                        htc_w_m2_k: 0.0,
                        reservoir: None,
                        source_w: 0.0,
                    },
                }),
            });
        }
    }
    let state = DepositionSimulation::new(
        events,
        &[(MaterialId(0), &SENSIBLE as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 3,
            max_links: 0,
            max_events: 3,
            max_substeps: 100,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-12,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    let route = MovingFrictionRoute::new(
        &state,
        StripFrame {
            origin_m: [0.0; 3],
            tangent: [1.0, 0.0, 0.0],
            normal: [0.0, 1.0, 0.0],
            width_m: 1.0,
        },
        faces,
        [0.0, 1.0],
        MovingFrictionBudget {
            maximum_surface_cells: 3,
            maximum_overlap_tests: 4,
            maximum_patches: 4,
        },
        || Ok(()),
    )
    .unwrap();
    (state, route)
}
fn native_controls() -> EllipticDrivenControls {
    EllipticDrivenControls {
        maximum_step_s: 1.0,
        minimum_step_s: 1e-10,
        maximum_substeps: 1000,
        rest_velocity_tolerance_m_s: 1e-12,
        solve: TangentialSolveConfig {
            max_iterations: 128,
            maximum_force_n: 100.0,
            force_absolute_tolerance_n: 1e-12,
            force_relative_tolerance: 1e-12,
        },
        relaxation: EllipticRelaxationControls {
            time_tolerance_s: 1e-11,
            travel_tolerance_m: 1e-11,
            maximum_evaluations: 40_000,
            maximum_root_iterations: 96,
        },
    }
}
fn input(work: &[f64]) -> MovingFrictionInput<()> {
    MovingFrictionInput {
        dissipated_j: work.to_vec(),
        conductivity_w_m_k: (0..3).map(|i| (CellId(i), 1.0)).collect(),
        transport: StepCoefficients {
            sources_w: BTreeMap::new(),
            conductances_w_k: BTreeMap::new(),
            payload: (),
        },
    }
}
fn tolerance() -> FrictionHeatTolerances {
    FrictionHeatTolerances {
        cell_j: 1e-10,
        thermal_j: 1e-10,
        interface_j: 1e-10,
        total_j: 1e-10,
    }
}
fn path() -> QuadraticContactPath {
    QuadraticContactPath {
        start: position(0.0),
        midpoint: position(0.25),
        end: position(1.0),
    }
}
fn thermal_endpoint() -> ([f64; 3], [f64; 2]) {
    // Independent elimination of two massless interfaces:
    // 2x-(2/3)y-(1/3)z=1/2, -(2/3)x+(5/3)y=1/4,
    // -(1/3)x+(4/3)z=1/4. These are the three original one-joule/K stores.
    let energy = [53.0 / 132.0, 41.0 / 132.0, 38.0 / 132.0];
    let interfaces = [
        300.0 + 0.5 * (energy[0] + energy[1]) + 3.0 / 16.0,
        300.0 + 0.5 * (energy[0] + energy[2]) + 3.0 / 8.0,
    ];
    (energy, interfaces)
}

#[test]
fn native_accelerating_traction_heats_the_actual_material_residence() {
    let (thermal, route) = fixture();
    // Original one-kilogram translating slider against a prescribed guide:
    // external force=3 N, native traction=1 N, so v=2t and x=t².
    let mass = thermal.cells()[&CellId(0)].mass_kg;
    let law = EllipticCoulomb::new([1.0; 2], [1.0; 2], [1.0, 0.0]).unwrap();
    let mut segments = Vec::new();
    let native = law
        .advance_driven_segments(
            &law.initial_state(),
            &EllipticDrivenInput {
                velocity_m_s: [0.0; 2],
                mobility_per_kg: [[1.0 / mass, 0.0], [0.0, 1.0 / mass]],
                driving_acceleration_m_s2: [3.0 / mass, 0.0],
                mode: EllipticDrivenMode::Resting,
                normal_n: 1.0,
                dt_s: 1.0,
                maximum_force_n: 100.0,
                contact_temperature_k: None,
                nominal_pressure_pa: None,
            },
            native_controls(),
            &mut || false,
            &mut |segment| {
                segments.push(*segment);
                Ok::<(), Infallible>(())
            },
        )
        .unwrap();
    assert_eq!(segments.len(), 1);
    let segment = segments[0];
    let acceleration = 3.0 / mass - segment.traction_n[0] / mass;
    let displacement = |t: f64| segment.velocity_before_m_s[0] * t + 0.5 * acceleration * t * t;
    let actual_path = QuadraticContactPath {
        start: position(0.0),
        midpoint: position(displacement(0.5 * segment.dt_s)),
        end: position(displacement(segment.dt_s)),
    };
    close(actual_path.midpoint.footprint_offset_m, 0.25);
    close(actual_path.end.footprint_offset_m, 1.0);
    close(segment.dissipated_j, 1.0);
    let step =
        MovingFrictionHeatStep::new_quadratic(&thermal, segment.dt_s, &route, actual_path, || {
            Ok(())
        })
        .unwrap();
    let powers = [
        segment.traction_n[0] * segment.velocity_before_m_s[0],
        segment.traction_n[0] * segment.velocity_after_m_s[0],
    ];
    let work = step.affine_power_work_j(powers).unwrap();
    assert_eq!(step.patches().len(), 2);
    close(step.patches()[0].mean_area_m2, 2.0 / 3.0);
    close(step.patches()[1].mean_area_m2, 1.0 / 3.0);
    close(work[0], 0.5);
    close(work[1], 0.5);
    close(work.iter().sum(), native.dissipated_j);
    let uniform = step.uniform_work_j(native.dissipated_j).unwrap();
    assert!((uniform[0] - work[0]).abs() > 0.16);
    let (energy, interfaces) = thermal_endpoint();
    let trial = step
        .evaluate_with_entropy(
            &energy,
            &interfaces,
            |_, _, _| Ok::<_, DepositionError>(input(&work)),
            || Ok(()),
        )
        .unwrap();
    for residual in trial.residual_j() {
        close(*residual, 0.0);
    }
    let (accepted, (), report) = trial.accept(tolerance(), || Ok(())).unwrap();
    close(report.dissipated_j, native.dissipated_j);
    close(report.thermal.stored_energy_change_j, 1.0);
    close(report.total_energy_residual_j, 0.0);
    close(report.body_heat_j[0], energy[0]);
    close(report.body_heat_j[1], energy[1] + energy[2]);
    // Conductance/storage determine the body heat exchange; no equal split.
    assert!((report.body_heat_j[0] - 0.5).abs() > 0.09);
    close(
        native.kinetic_energy_change_j + report.thermal.stored_energy_change_j,
        native.drive_work_j,
    );
    for (i, expected) in energy.into_iter().enumerate() {
        close(
            accepted.temperature_k(CellId(i as u64)).unwrap(),
            300.0 + expected,
        );
        close(thermal.temperature_k(CellId(i as u64)).unwrap(), 300.0);
    }
    let entropy = report.entropy.unwrap();
    close(
        entropy.thermal_storage_change_j_k,
        energy.into_iter().map(|h| (h / 300.0).ln_1p()).sum(),
    );
    close(
        entropy.friction_generation_j_k,
        work[0] / interfaces[0] + work[1] / interfaces[1],
    );
    close(
        entropy.thermal_storage_change_j_k,
        entropy.reservoir_entropy_input_j_k
            + entropy.network_conduction_generation_j_k
            + entropy.interface_conduction_generation_j_k
            + entropy.friction_generation_j_k
            + entropy.unclassified_source_at_cell_j_k
            + entropy.interface_residual_entropy_j_k
            + entropy.temporal_defect_j_k,
    );
}

#[test]
fn quadratic_reversals_and_linear_reduction_use_complete_original_surface_coverage() {
    let (thermal, route) = fixture();
    let reversal = QuadraticContactPath {
        start: position(0.0),
        midpoint: position(1.0),
        end: position(0.0),
    };
    let step =
        MovingFrictionHeatStep::new_quadratic(&thermal, 1.0, &route, reversal, || Ok(())).unwrap();
    close(step.patches()[0].mean_area_m2, 1.0 / 3.0);
    close(step.patches()[1].mean_area_m2, 2.0 / 3.0);
    assert_eq!(step.patches()[0].cells, [CellId(0), CellId(1)]);
    assert_eq!(step.patches()[1].cells, [CellId(0), CellId(2)]);
    let static_endpoints = MovingFrictionHeatStep::new(
        &thermal,
        1.0,
        &route,
        MovingContactPath {
            start: reversal.start,
            end: reversal.end,
        },
        || Ok(()),
    )
    .unwrap();
    assert_eq!(static_endpoints.patches().len(), 1);
    // The interior excursion leaves body B; agreeing endpoints cannot admit it.
    assert!(
        MovingFrictionHeatStep::new_quadratic(
            &thermal,
            1.0,
            &route,
            QuadraticContactPath {
                midpoint: position(1.1),
                ..reversal
            },
            || Ok(())
        )
        .is_err()
    );
    // Vertex t=.625 leaves B even though start/midpoint/end all fit.
    assert!(
        MovingFrictionHeatStep::new_quadratic(
            &thermal,
            1.0,
            &route,
            QuadraticContactPath {
                start: position(0.0),
                midpoint: position(1.0),
                end: position(2.0 / 3.0)
            },
            || Ok(())
        )
        .is_err()
    );
    let linear_path = MovingContactPath {
        start: position(0.0),
        end: position(1.0),
    };
    let linear =
        MovingFrictionHeatStep::new(&thermal, 1.0, &route, linear_path, || Ok(())).unwrap();
    let quadratic = MovingFrictionHeatStep::new_quadratic(
        &thermal,
        1.0,
        &route,
        QuadraticContactPath {
            start: linear_path.start,
            midpoint: position(0.5),
            end: linear_path.end,
        },
        || Ok(()),
    )
    .unwrap();
    let power = quadratic.affine_power_work_j([1.0; 2]).unwrap();
    let uniform = linear.uniform_work_j(1.0).unwrap();
    for (a, b) in power.into_iter().zip(uniform) {
        close(a, b);
    }
    for (a, b) in linear.patches().iter().zip(quadratic.patches()) {
        assert_eq!(a.cells, b.cells);
        close(a.mean_area_m2, b.mean_area_m2);
    }
    assert!(quadratic.affine_power_work_j([-1.0, 2.0]).is_err());
}

#[test]
fn late_source_cancellation_and_residual_refusal_preserve_original_thermal_owner() {
    let (thermal, route) = fixture();
    let step =
        MovingFrictionHeatStep::new_quadratic(&thermal, 1.0, &route, path(), || Ok(())).unwrap();
    let work = step.affine_power_work_j([0.0, 2.0]).unwrap();
    let (energy, interfaces) = thermal_endpoint();
    let source_failure = step.evaluate(
        &energy,
        &interfaces,
        |_, _, _| Err::<MovingFrictionInput<()>, _>(DepositionError::MaterialDomain),
        || Ok(()),
    );
    assert!(matches!(
        source_failure,
        Err(DepositionError::MaterialDomain)
    ));
    let cancelled = Cell::new(false);
    let interrupted = step.evaluate(
        &energy,
        &interfaces,
        |_, _, _| {
            cancelled.set(true);
            Ok::<_, DepositionError>(input(&work))
        },
        || {
            if cancelled.get() {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        },
    );
    assert!(matches!(interrupted, Err(DepositionError::Cancelled)));
    let wrong = step
        .evaluate(
            &[0.0; 3],
            &interfaces,
            |_, _, _| Ok::<_, DepositionError>(input(&work)),
            || Ok(()),
        )
        .unwrap();
    assert!(wrong.accept(tolerance(), || Ok(())).is_err());
    assert_eq!(thermal.time_s(), 0.0);
    for cell in thermal.cells().values() {
        assert_eq!(cell.specific_h_j_kg, 0.0);
    }
    let mut expected = None;
    for _ in 0..2 {
        let trial = step
            .evaluate(
                &energy,
                &interfaces,
                |_, _, _| Ok::<_, DepositionError>(input(&work)),
                || Ok(()),
            )
            .unwrap();
        let (accepted, (), report) = trial.accept(tolerance(), || Ok(())).unwrap();
        let values = (
            accepted.time_s(),
            accepted
                .cells()
                .values()
                .map(|c| c.specific_h_j_kg)
                .collect::<Vec<_>>(),
            report.body_heat_j,
            report.dissipated_j,
        );
        if let Some(previous) = &expected {
            assert_eq!(&values, previous);
        } else {
            expected = Some(values);
        }
    }
}
