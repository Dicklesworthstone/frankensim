//! G0/G1/G3: conservative moving bands and an independent half-space integral.
//! Thermal coefficients are analytical fixtures, not measured material cards.
use fs_couple::friction_heat::{FrictionHeatTolerances, moving::*};
use fs_math::det;
use fs_solver::{LinearOp, krylov::CgState};
use fs_sparse::precond::Precond;
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::*;
use fs_tribo::{FrictionLaw, TangentialLaw, TangentialSolveConfig, TangentialStepInput};
use std::collections::BTreeMap;

struct Caloric;
impl ThermalLaw for Caloric {
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
}
static CALORIC: Caloric = Caloric;
const K: [f64; 2] = [1.0, 2.0]; // rho*c = k: equal diffusivity, unequal effusivity.
fn controls() -> FrictionHeatTolerances {
    FrictionHeatTolerances {
        cell_j: 1e-10,
        thermal_j: 1e-10,
        interface_j: 1e-10,
        total_j: 1e-10,
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.16e} != {b:.16e}, tolerance {tolerance:.3e}"
    );
}
fn position(x: f64) -> MovingContactPosition {
    MovingContactPosition {
        body_offsets_m: [0.0; 2],
        footprint_offset_m: x,
    }
}
fn frame() -> StripFrame {
    StripFrame {
        origin_m: [0.0; 3],
        tangent: [1.0, 0.0, 0.0],
        normal: [0.0, 1.0, 0.0],
        width_m: 1.0,
    }
}
fn budget(nx: usize) -> MovingFrictionBudget {
    MovingFrictionBudget {
        maximum_surface_cells: 2 * nx,
        maximum_overlap_tests: 4 * nx * nx,
        maximum_patches: nx * nx,
    }
}
fn grid(nx: usize, ny: usize) -> (DepositionSimulation<'static>, [Vec<StripSurfaceCell>; 2]) {
    let x: Vec<_> = (0..=nx)
        .map(|i| -2.0 + 5.0 * i as f64 / nx as f64)
        .collect();
    // Graded normal grid; all rungs keep the same finite-domain boundaries.
    let y: Vec<_> = (0..=ny)
        .map(|j| 1.5 * ((3.0 * j as f64 / ny as f64).exp() - 1.0) / (3.0_f64.exp() - 1.0))
        .collect();
    let mut events = Vec::new();
    let mut faces: [Vec<StripSurfaceCell>; 2] = std::array::from_fn(|_| Vec::new());
    let id = |side: usize, i: usize, j: usize| CellId((side * nx * ny + j * nx + i) as u64);
    for side in 0..2 {
        for j in 0..ny {
            for i in 0..nx {
                let dx = x[i + 1] - x[i];
                let dy = y[j + 1] - y[j];
                let xc = f64::midpoint(x[i], x[i + 1]);
                let (ya, yb) = if side == 0 {
                    (-y[j + 1], -y[j])
                } else {
                    (y[j], y[j + 1])
                };
                events.push(TimedEvent {
                    time_s: 0.0,
                    event: Event::Birth(DepositedCell {
                        id: id(side, i, j),
                        segment_id: id(side, i, j).0,
                        body_id: side as u64,
                        material_id: MaterialId(0),
                        mass_kg: K[side] * dx * dy,
                        specific_h_j_kg: 0.0,
                        geometry: CellGeometry {
                            start_m: [xc, ya, 0.0],
                            end_m: [xc, yb, 0.0],
                            volume_m3: dx * dy,
                            surface_area_m2: 2.0 * (dx * dy + dx + dy),
                        },
                        boundary: Boundary {
                            ambient_k: 300.0,
                            htc_w_m2_k: 0.0,
                            reservoir: None,
                            source_w: 0.0,
                        },
                    }),
                });
                if j == 0 {
                    faces[side].push(StripSurfaceCell {
                        cell: id(side, i, j),
                        span_m: [x[i], x[i + 1]],
                    });
                }
            }
        }
    }
    let mut link = 0;
    for side in 0..2 {
        for j in 0..ny {
            for i in 0..nx {
                let dx = x[i + 1] - x[i];
                let dy = y[j + 1] - y[j];
                if i + 1 < nx {
                    events.push(TimedEvent {
                        time_s: 0.0,
                        event: Event::Contact(ThermalLink {
                            id: link,
                            a: id(side, i, j),
                            b: id(side, i + 1, j),
                            area_m2: dy,
                            conductance_w_k: K[side] * dy / (0.5 * (x[i + 2] - x[i])),
                        }),
                    });
                    link += 1;
                }
                if j + 1 < ny {
                    events.push(TimedEvent {
                        time_s: 0.0,
                        event: Event::Contact(ThermalLink {
                            id: link,
                            a: id(side, i, j),
                            b: id(side, i, j + 1),
                            area_m2: dx,
                            conductance_w_k: K[side] * dx / (0.5 * (y[j + 2] - y[j])),
                        }),
                    });
                    link += 1;
                }
            }
        }
    }
    let count = events.len();
    let state = DepositionSimulation::new(
        events,
        &[(MaterialId(0), &CALORIC as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 2 * nx * ny,
            max_links: link as usize,
            max_events: count,
            max_substeps: 10000,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.1,
            energy_tolerance_j: 1e-12,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    (state, faces)
}
struct Matrix {
    storage: Vec<f64>,
    diagonal: Vec<f64>,
    edges: Vec<(usize, usize, f64)>,
}
impl LinearOp for Matrix {
    fn n(&self) -> usize {
        self.storage.len()
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        for ((yi, &a), &v) in y.iter_mut().zip(&self.storage).zip(x) {
            *yi = a * v;
        }
        for &(a, b, g) in &self.edges {
            let q = g * (x[a] - x[b]);
            y[a] += q;
            y[b] -= q;
        }
    }
}
impl Precond for Matrix {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        for ((zi, &ri), &d) in z.iter_mut().zip(r).zip(&self.diagonal) {
            *zi = ri / d;
        }
    }
}
// Independent linear algebra in temperature-rise coordinates. The production
// owner accepts only its original joule residuals, never this solver's estimate.
fn solve(
    thermal: &DepositionSimulation<'_>,
    step: &MovingFrictionHeatStep<'_, '_>,
    work: f64,
) -> (Vec<f64>, Vec<f64>, Vec<f64>, BTreeMap<CellId, f64>) {
    let dt = step.duration_s();
    let storage: Vec<_> = thermal.cells().values().map(|c| c.mass_kg / dt).collect();
    let mut matrix = Matrix {
        diagonal: storage.clone(),
        storage,
        edges: Vec::new(),
    };
    let mut rhs: Vec<_> = thermal
        .cells()
        .values()
        .zip(&matrix.storage)
        .map(|(c, &m)| m * c.specific_h_j_kg)
        .collect();
    let mut add = |a: usize, b: usize, g: f64| {
        matrix.edges.push((a, b, g));
        matrix.diagonal[a] += g;
        matrix.diagonal[b] += g;
    };
    for link in thermal.links().values() {
        add(link.a.0 as usize, link.b.0 as usize, link.conductance_w_k);
    }
    let heat = step.uniform_work_j(work).unwrap();
    let mut conductivities = BTreeMap::new();
    for (p, &h) in step.patches().iter().zip(&heat) {
        let ks = p.cells.map(|id| K[thermal.cells()[&id].body_id as usize]);
        let g = [
            p.mean_area_m2 * ks[0] / p.distances_m[0],
            p.mean_area_m2 * ks[1] / p.distances_m[1],
        ];
        add(
            p.cells[0].0 as usize,
            p.cells[1].0 as usize,
            g[0] * g[1] / (g[0] + g[1]),
        );
        for side in 0..2 {
            rhs[p.cells[side].0 as usize] += h / dt * g[side] / (g[0] + g[1]);
            conductivities.insert(p.cells[side], ks[side]);
        }
    }
    let mut cg = CgState::new(&matrix, &matrix, &rhs);
    let result = cg.run(&matrix, &matrix, 1e-13, 400);
    assert!(result.converged, "CG report {result:?}");
    let energy = thermal
        .cells()
        .values()
        .zip(&cg.x)
        .map(|(c, &rise)| c.mass_kg * (rise - c.specific_h_j_kg))
        .collect();
    let temperatures = step
        .patches()
        .iter()
        .zip(&heat)
        .map(|(p, &h)| {
            let g = p.cells.map(|id| {
                p.mean_area_m2 * conductivities[&id]
                    / p.distances_m[if id == p.cells[0] { 0 } else { 1 }]
            });
            300.0
                + (h / dt + g[0] * cg.x[p.cells[0].0 as usize] + g[1] * cg.x[p.cells[1].0 as usize])
                    / (g[0] + g[1])
        })
        .collect();
    (energy, temperatures, heat, conductivities)
}
fn input(work: Vec<f64>, k: BTreeMap<CellId, f64>) -> MovingFrictionInput<()> {
    MovingFrictionInput {
        dissipated_j: work,
        conductivity_w_m_k: k,
        transport: StepCoefficients {
            sources_w: BTreeMap::new(),
            conductances_w_k: BTreeMap::new(),
            payload: (),
        },
    }
}

#[test]
fn g0_fast_continuous_footprint_crossings_conserve_original_work_and_reverse_geometrically() {
    let (thermal, faces) = grid(10, 1);
    let route =
        MovingFrictionRoute::new(&thermal, frame(), faces, [-0.25, 0.25], budget(10), || {
            Ok(())
        })
        .unwrap();
    let path = MovingContactPath {
        start: position(-1.0),
        end: MovingContactPosition {
            body_offsets_m: [0.0, 0.2],
            footprint_offset_m: 2.0,
        },
    };
    let step = MovingFrictionHeatStep::new(&thermal, 0.1, &route, path, || Ok(())).unwrap();
    let reverse = MovingFrictionHeatStep::new(
        &thermal,
        0.1,
        &route,
        MovingContactPath {
            start: path.end,
            end: path.start,
        },
        || Ok(()),
    )
    .unwrap();
    assert!(
        step.patches().len() > 6,
        "whole swept trail must receive heat"
    );
    close(
        step.patches().iter().map(|p| p.mean_area_m2).sum(),
        0.5,
        1e-14,
    );
    assert_eq!(step.patches().len(), reverse.patches().len());
    for (a, b) in step.patches().iter().zip(reverse.patches()) {
        assert_eq!(a.cells, b.cells);
        close(a.mean_area_m2, b.mean_area_m2, 1e-14);
    }
    let law = FrictionLaw::Coulomb {
        static_mu: 0.5,
        kinetic_mu: 0.4,
    };
    let accepted = law.initial_state();
    let friction = law
        .advance(
            &accepted,
            &TangentialStepInput {
                normal_n: 2.0,
                slip_increment_m: [0.01, 0.0],
                dt_s: 0.1,
                contact_temperature_k: Some(300.0),
                nominal_pressure_pa: Some(4.0),
                compliance_m_per_n: 0.0,
                solve: TangentialSolveConfig {
                    max_iterations: 32,
                    maximum_force_n: 10.0,
                    force_absolute_tolerance_n: 1e-12,
                    force_relative_tolerance: 1e-12,
                },
            },
        )
        .unwrap();
    let (energy, t, heat, k) = solve(&thermal, &step, friction.dissipated_j);
    close(heat.iter().sum(), friction.dissipated_j, 1e-16);
    let trial = step
        .evaluate(
            &energy,
            &t,
            |_, _, _| Ok::<_, DepositionError>(input(heat, k)),
            || Ok(()),
        )
        .unwrap();
    let (next, (), report) = trial.accept(controls(), || Ok(())).unwrap();
    close(next.energy_j(), friction.dissipated_j, 1e-12);
    close(
        report.body_heat_j.iter().sum(),
        friction.dissipated_j,
        1e-12,
    );
    assert_eq!(report.position, path.end);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.energy_j(), 0.0);
    assert_eq!(accepted.slip_m(), [0.0; 2]);
}

#[test]
fn g0_missing_coverage_wrong_work_and_cancelled_trials_leave_original_heat_unchanged() {
    let (thermal, faces) = grid(4, 1);
    let route =
        MovingFrictionRoute::new(&thermal, frame(), faces, [-1.0, 1.0], budget(4), || Ok(()))
            .unwrap();
    assert!(
        MovingFrictionHeatStep::new(
            &thermal,
            0.1,
            &route,
            MovingContactPath {
                start: position(0.0),
                end: position(3.0)
            },
            || Ok(())
        )
        .is_err()
    );
    let step = MovingFrictionHeatStep::new(
        &thermal,
        0.1,
        &route,
        MovingContactPath {
            start: position(0.0),
            end: position(0.5),
        },
        || Ok(()),
    )
    .unwrap();
    let (energy, t, mut heat, k) = solve(&thermal, &step, 0.1);
    heat[0] += 0.01;
    let trial = step
        .evaluate(
            &energy,
            &t,
            |_, _, _| Ok::<_, DepositionError>(input(heat, k)),
            || Ok(()),
        )
        .unwrap();
    assert!(trial.accept(controls(), || Ok(())).is_err());
    let (energy, t, heat, k) = solve(&thermal, &step, 0.1);
    assert!(matches!(
        step.evaluate(
            &energy,
            &t,
            |_, _, _| Ok::<_, DepositionError>(input(heat, k)),
            || Err(DepositionError::Cancelled)
        ),
        Err(DepositionError::Cancelled)
    ));
    assert_eq!(thermal.energy_j(), 0.0);
    assert_eq!(thermal.time_s(), 0.0);
}

// Transient Jaeger band-source Green-function integral, evaluated by independent
// Simpson quadrature. Integrate the Neumann half-space heat kernel first over
// source and observation intervals, then use tau=z² to remove its 1/sqrt(tau)
// singularity. q_a/k_a=1 and kappa=1. This is not a fitted temperature closure.
// Jaeger, Moving sources of heat and the temperature at sliding contacts,
// Proc. Royal Soc. NSW 76, 203-224, doi:10.5962/p.360338.
fn jaeger_average(speed: f64, time: f64) -> f64 {
    jaeger_interval(speed, time, [-1.0, 1.0])
}
fn jaeger_interval(speed: f64, time: f64, span: [f64; 2]) -> f64 {
    let n = 4096;
    let dz = time.sqrt() / n as f64;
    let root_pi = std::f64::consts::PI.sqrt();
    let f = |z: f64| {
        if z == 0.0 {
            return 2.0 / root_pi;
        }
        let c = speed * z * z;
        let d = 2.0 * z;
        let primitive = |x: f64| x * det::erf(x / d) + d / root_pi * (-(x / d).powi(2)).exp();
        (primitive(span[1] + c + 1.0) - primitive(span[0] + c + 1.0) - primitive(span[1] + c - 1.0)
            + primitive(span[0] + c - 1.0))
            / ((span[1] - span[0]) * root_pi)
    };
    let mut sum = f(0.0) + f(time.sqrt());
    for i in 1..n {
        sum += if i % 2 == 0 { 2.0 } else { 4.0 } * f(i as f64 * dz);
    }
    dz * sum / 3.0
}
fn band(
    thermal: &DepositionSimulation<'static>,
    route: &MovingFrictionRoute,
    speed: f64,
    steps: usize,
) -> (f64, f64, f64) {
    let end = 0.02;
    let dt = end / steps as f64;
    let mut state = thermal.clone();
    let mut temperature = 0.0;
    let mut field_error = 0.0;
    for i in 0..steps {
        let path = MovingContactPath {
            start: position(speed * i as f64 * dt),
            end: position(speed * (i + 1) as f64 * dt),
        };
        let step = MovingFrictionHeatStep::new(&state, (i + 1) as f64 * dt, route, path, || Ok(()))
            .unwrap();
        let work = 6.0 * step.duration_s(); // flux q_a=1, q_b=2 over a two-metre band
        let (energy, t, heat, k) = solve(&state, &step, work);
        let (next, (), report) = step
            .evaluate(
                &energy,
                &t,
                |_, _, _| Ok::<_, DepositionError>(input(heat, k)),
                || Ok(()),
            )
            .unwrap()
            .accept(controls(), || Ok(()))
            .unwrap();
        close(report.body_heat_j[0] / work, 1.0 / 3.0, 1e-8);
        temperature = report
            .patches
            .iter()
            .map(|p| p.mean_area_m2 * (p.interface_temperature_k - 300.0))
            .sum::<f64>()
            / 2.0;
        if i + 1 == steps {
            let mut error = 0.0;
            let mut reference = 0.0;
            for p in &report.patches {
                let cell = &state.cells()[&p.cells[0]];
                let depth = (cell.geometry.end_m[1] - cell.geometry.start_m[1]).abs();
                let dx = cell.geometry.volume_m3 / depth;
                let center = cell.geometry.end_m[0] - speed * end;
                let span = [(center - 0.5 * dx).max(-1.0), (center + 0.5 * dx).min(1.0)];
                if span[1] <= span[0] {
                    continue;
                }
                let exact = jaeger_interval(speed, end, span);
                error += (span[1] - span[0]) * (p.interface_temperature_k - 300.0 - exact).powi(2);
                reference += (span[1] - span[0]) * exact * exact;
            }
            field_error = (error / reference).sqrt();
        }
        state = next;
    }
    close(state.energy_j(), 6.0 * end, 1e-9);
    let surface = |x: f64| {
        state
            .cells()
            .values()
            .filter(|c| c.body_id == 0 && c.geometry.end_m[1] == 0.0)
            .min_by(|a, b| {
                (a.geometry.end_m[0] - x)
                    .abs()
                    .total_cmp(&(b.geometry.end_m[0] - x).abs())
            })
            .unwrap()
            .specific_h_j_kg
    };
    let wake = surface(speed * end - 1.2) - surface(speed * end + 1.2);
    (temperature, wake, field_error)
}

#[test]
fn g1_moving_bands_match_jaeger_at_three_peclet_numbers_and_refine_in_space_and_time() {
    let pe = [0.1, 1.0, 10.0]; // Pe=v*a/(2*kappa), a=1, kappa=1.
    let references = pe.map(|p| jaeger_average(2.0 * p, 0.02));
    let mut errors = Vec::new();
    let mut finest = None;
    let mut finest_value = 0.0;
    for (nx, ny) in [(20, 8), (40, 16), (80, 32)] {
        let (thermal, faces) = grid(nx, ny);
        let route =
            MovingFrictionRoute::new(&thermal, frame(), faces, [-1.0, 1.0], budget(nx), || Ok(()))
                .unwrap();
        let mut rung = [0.0; 3];
        for j in 0..3 {
            let (actual, wake, field) = band(&thermal, &route, 2.0 * pe[j], 128);
            rung[j] = field;
            let mean_error = (actual / references[j] - 1.0).abs();
            eprintln!(
                "Pe={} grid={}x{} average={} reference={} mean_error={} field_L2={}",
                pe[j], nx, ny, actual, references[j], mean_error, field
            );
            if nx == 80 {
                assert!(mean_error < 0.02);
            }
            if j == 2 && nx == 80 {
                finest_value = actual;
                assert!(
                    wake > 0.005,
                    "the travelled material must retain a thermal wake"
                );
            }
        }
        errors.push(rung);
        if nx == 80 {
            finest = Some((thermal, route));
        }
    }
    for j in 0..3 {
        assert!(errors[2][j] < 0.02, "Pe={} errors {:?}", pe[j], errors);
        assert!(
            errors[1][j] < errors[0][j] && errors[2][j] < errors[1][j],
            "live spatial ladder {:?}",
            errors
        );
    }
    // Use the finest transition of the uniform bound over the declared Pe
    // range. A coarse cut face need not yet be in its asymptotic regime.
    let worst: Vec<f64> = errors
        .iter()
        .map(|r| r.iter().copied().fold(0.0, f64::max))
        .collect();
    let order = (worst[1] / worst[2]).log2();
    eprintln!(
        "spatial worst-field errors {:?}; fine-transition order {}",
        worst, order
    );
    assert!(order > 0.9);
    let (thermal, route) = finest.unwrap();
    // Refine each spatial axis independently as well: neither an axial wake
    // nor the near-surface penetration field may be frozen on the ladder.
    for rungs in [[(80, 8), (80, 16)], [(20, 32), (40, 32)]] {
        let mut values = Vec::new();
        for (nx, ny) in rungs {
            let (state, faces) = grid(nx, ny);
            let bound = MovingFrictionRoute::new(
                &state,
                frame(),
                faces,
                [-1.0, 1.0],
                budget(nx),
                || Ok(()),
            )
            .unwrap();
            values.push(band(&state, &bound, 20.0, 128).0);
        }
        values.push(finest_value);
        let ratio = ((values[0] - values[1]) / (values[1] - values[2])).abs();
        eprintln!(
            "axis {:?} temperature ladder {:?}; observed order {}",
            rungs,
            values,
            ratio.log2()
        );
        assert!(ratio > 1.0 && (values[1] - values[2]).abs() > 1e-7);
    }
    let temporal = [
        band(&thermal, &route, 20.0, 32).0,
        band(&thermal, &route, 20.0, 64).0,
        finest_value,
    ];
    let ratio = ((temporal[0] - temporal[1]) / (temporal[1] - temporal[2])).abs();
    eprintln!(
        "time ladder {:?}; observed order {}",
        temporal,
        ratio.log2()
    );
    assert!(ratio.log2() > 0.9 && (temporal[0] - temporal[1]).abs() > 1e-7);
    // The physical equal-diffusivity effusivity share is 1/(1+2); band()
    // checks it on every accepted step. Pe->0 recovers the stationary band.
    let stationary = jaeger_average(0.0, 0.02);
    close(jaeger_average(2e-5, 0.02), stationary, 1e-10);
    let low = band(&thermal, &route, 2e-5, 128).0;
    assert!((low / stationary - 1.0).abs() < 0.02);
}
