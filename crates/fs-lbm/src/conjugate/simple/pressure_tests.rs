//! G1/G3 pressure-driven flow through the existing MAC momentum and heat owners.

use super::*;
use crate::conjugate::{
    EnergyConfig, ThermalFace, ThermalSetup, UnsteadyConfig, simple_unsteady, solve_energy,
};

fn fluid() -> FluidProperties {
    FluidProperties {
        density_kg_m3: 1.0,
        specific_heat_j_kg_k: 2.0,
        conductivity_w_m_k: 0.01,
        kinematic_viscosity_m2_s: 0.25,
    }
}

fn channel(n: usize, axis: usize, inlet: f64, outlet: f64) -> (VoxelDomain, SimpleConfig) {
    let transverse = (axis + 1) % 3;
    let mut widths = [vec![1.0], vec![1.0], vec![1.0]];
    widths[axis] = vec![0.25; 8];
    widths[transverse] = vec![1.0 / n as f64; n];
    let domain = VoxelDomain::graded(widths).unwrap();
    let mut faces = [FvBoundary::Symmetry; 6];
    faces[2 * axis] = FvBoundary::Pressure { pressure_pa: inlet };
    faces[2 * axis + 1] = FvBoundary::Pressure {
        pressure_pa: outlet,
    };
    faces[2 * transverse] = FvBoundary::wall();
    faces[2 * transverse + 1] = FvBoundary::wall();
    let mut config = SimpleConfig::new(faces);
    config.tolerance = 1e-10;
    config.pressure_tolerance = 1e-12;
    (domain, config)
}

#[test]
fn pressure_opening_poiseuille_refines_and_preserves_the_face_pressure_datum() {
    let gate = CancelGate::new();
    let mut errors = Vec::new();
    for n in [4, 8] {
        let (domain, config) = channel(n, 0, 7.5, 7.0);
        let flow = simple_flow(&domain, &fluid(), &config, &gate).unwrap();
        // Fully developed planar flow, H=1, L=2, mu=1/4, depth=1:
        // u(y)=dp/(2 mu L)*y*(H-y), Q=dp*H^3/(12 mu L).
        let exact_q = 0.5 / (12.0 * 0.25 * 2.0);
        let q = flow.field.boundary_outflow(&domain, Face3::XMax);
        errors.push((q - exact_q).abs());
        assert!(q > 0.0);
        assert!((q / exact_q - (1.0 + 2.0 / (n * n) as f64)).abs() < 1e-7);
        assert!((flow.report.inflow_m3_s - q).abs() < 1e-12);
        assert!((flow.report.outflow_m3_s - q).abs() < 1e-12);
        assert!(flow.report.max_divergence_m3_s < 1e-12);
        for c in 0..domain.cell_count() {
            let at = domain.coords(c);
            let x = domain.center(at[0], at[1], at[2])[0];
            assert!((flow.pressure_pa[c] - (7.5 - 0.25 * x)).abs() < 2e-7);
        }
    }
    assert!((errors[0] / errors[1] - 4.0).abs() < 1e-5, "{errors:?}");
}

#[test]
fn pressure_opening_offsets_and_equal_pressure_do_not_create_motion() {
    let gate = CancelGate::new();
    let (domain, config) = channel(4, 0, 0.5, 0.0);
    let base = simple_flow(&domain, &fluid(), &config, &gate).unwrap();
    for offset in [-1024.0, 1_048_576.0] {
        let (_, shifted) = channel(4, 0, 0.5 + offset, offset);
        let flow = simple_flow(&domain, &fluid(), &shifted, &gate).unwrap();
        assert_eq!(flow.velocity_m_s, base.velocity_m_s);
        assert_eq!(flow.report, base.report);
        for (pressure, original) in flow.pressure_pa.iter().zip(&base.pressure_pa) {
            assert!((pressure - offset - original).abs() < 2e-10);
        }
        let (_, rest) = channel(4, 0, offset, offset);
        let rest = simple_flow(&domain, &fluid(), &rest, &gate).unwrap();
        assert!(rest.velocity_m_s.iter().all(|v| *v == [0.0; 3]));
        assert!(rest.pressure_pa.iter().all(|p| *p == offset));
        assert_eq!(rest.report.inflow_m3_s, 0.0);
        assert_eq!(rest.report.outflow_m3_s, 0.0);
    }
    let mut zero_outlet = config;
    zero_outlet.faces[1] = FvBoundary::Outlet;
    let flow = simple_flow(&domain, &fluid(), &zero_outlet, &gate).unwrap();
    assert_eq!(flow, base);
}

#[test]
fn pressure_opening_flow_reverses_on_every_axis_and_respects_material_viscosity() {
    let gate = CancelGate::new();
    let mut throughputs = Vec::new();
    for axis in 0..3 {
        for (left, right, sign) in [(0.5, 0.0, 1.0), (0.0, 0.5, -1.0)] {
            let (domain, config) = channel(4, axis, left, right);
            let mut rates = Vec::new();
            for viscosity in [0.25, 0.5] {
                let flow = simple_flow_with_viscosity(
                    &domain,
                    &fluid(),
                    &vec![viscosity; domain.cell_count()],
                    &config,
                    &gate,
                )
                .unwrap();
                let q = flow
                    .field
                    .boundary_outflow(&domain, Face3::ALL[2 * axis + 1]);
                assert!(q * sign > 0.0);
                assert!(flow.report.max_divergence_m3_s < 1e-12);
                assert!((flow.report.inflow_m3_s - q.abs()).abs() < 1e-12);
                rates.push(q.abs());
            }
            assert!(
                (rates[0] / rates[1] - 2.0).abs() < 2e-7,
                "axis {axis}: {rates:?}"
            );
            throughputs.push(rates[0]);
        }
    }
    assert!(
        throughputs
            .iter()
            .all(|q| (q / throughputs[0] - 1.0).abs() < 2e-7)
    );
}

#[test]
fn pressure_opening_backflow_uses_the_existing_thermal_reservoir() {
    let gate = CancelGate::new();
    for (left, right, temperature) in [(0.5, 0.0, 320.0), (0.0, 0.5, 280.0)] {
        let (domain, config) = channel(4, 0, left, right);
        let flow = simple_flow(&domain, &fluid(), &config, &gate).unwrap();
        let mut faces = [ThermalFace::Adiabatic; 6];
        faces[0] = ThermalFace::Outflow {
            backflow_temperature: 320.0,
        };
        faces[1] = ThermalFace::Outflow {
            backflow_temperature: 280.0,
        };
        let energy = solve_energy(
            &domain,
            &fluid(),
            &[],
            &flow.field,
            &ThermalSetup::new(faces),
            &EnergyConfig::default(),
            &gate,
        )
        .unwrap();
        assert!(
            energy
                .temperature
                .iter()
                .all(|t| (t - temperature).abs() < 1e-6)
        );
        assert!(energy.report.balance.residual_w.abs() < 1e-8);
    }
}

#[test]
fn pressure_opening_drives_the_original_transient_momentum_owner() {
    let gate = CancelGate::new();
    let (domain, config) = channel(4, 0, 0.5, 0.0);
    let steady = simple_flow(&domain, &fluid(), &config, &gate).unwrap();
    let mut time = UnsteadyConfig::new(0.1, 12);
    time.inner_tolerance = 1e-9;
    time.inner_iterations = 300;
    let started = simple_unsteady(&domain, &fluid(), &config, &time, |_| 1.0, &gate).unwrap();
    let energies: Vec<_> = started.records.iter().map(|r| r.kinetic_energy_j).collect();
    assert!(energies[0] > 0.0);
    assert!(energies.windows(2).all(|pair| pair[1] > pair[0]));
    let q = started.flow.field.boundary_outflow(&domain, Face3::XMax);
    assert!(q > 0.9 * steady.report.outflow_m3_s && q < steady.report.outflow_m3_s);
    assert!(started.flow.report.max_divergence_m3_s < 1e-12);
    let (_, shifted) = channel(4, 0, 1024.5, 1024.0);
    let shifted = simple_unsteady(&domain, &fluid(), &shifted, &time, |_| 1.0, &gate).unwrap();
    assert_eq!(started.records, shifted.records);
    assert_eq!(started.flow.velocity_m_s, shifted.flow.velocity_m_s);
}

#[test]
fn pressure_opening_refuses_nonfinite_and_unrepresentable_boundary_differences() {
    let gate = CancelGate::new();
    for pressure in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let (domain, config) = channel(4, 0, pressure, 0.0);
        assert!(matches!(
            simple_flow(&domain, &fluid(), &config, &gate),
            Err(ChtError::InvalidInput {
                field: "simple.pressure_pa",
                ..
            })
        ));
    }
    let (domain, config) = channel(4, 0, f64::MAX, -f64::MAX);
    assert!(matches!(
        simple_flow(&domain, &fluid(), &config, &gate),
        Err(ChtError::InvalidInput {
            field: "simple.pressure_difference_pa",
            ..
        })
    ));
    let (_, config) = channel(4, 0, 0.5, 0.0);
    gate.request();
    assert!(matches!(
        simple_flow(&domain, &fluid(), &config, &gate),
        Err(ChtError::Cancelled)
    ));
}

#[test]
fn pressure_projection_polishes_quiet_component_before_caloric_transport() {
    use crate::conjugate::{
        AdvectiveEnthalpyField, AffineSpecificEnthalpy, CellSink, SolidMaterial, Voxel,
    };

    // Five fluid cells: one drained cell and a separated 2x2 recirculating
    // pocket. A single global pressure solve can accept a large LOCAL defect
    // in that pocket because its flux is tiny relative to the other component.
    let mut domain = VoxelDomain::new(4, 2, 1, 1.0).unwrap();
    for [x, y] in [[1, 0], [1, 1], [0, 1]] {
        domain.set(x, y, 0, Voxel::Solid(0));
    }
    let mut faces = [FvBoundary::Symmetry; 6];
    faces[0] = FvBoundary::Outlet;
    faces[1] = FvBoundary::Outlet;
    let mut config = SimpleConfig::new(faces);
    config.pressure_solver = PressureSolver::IluBicgstab;
    let gate = CancelGate::new();
    let fluid = fluid();
    let epsilon = 1e-12;
    let seed = |solver: &mut Solver<'_>| {
        for axis in 0..3 {
            for &face in &solver.comps[axis].unknowns {
                solver.d[axis][face] = 1.0;
            }
        }
        for (axis, at, speed) in [
            (0, [0, 0, 0], 1.0),
            (0, [3, 0, 0], 2.0 * epsilon),
            (0, [3, 1, 0], -epsilon),
            (1, [2, 1, 0], -epsilon),
            (1, [3, 1, 0], epsilon),
        ] {
            let face = solver.comps[axis].index(at);
            solver.vel[axis][face] = speed;
        }
    };
    let mut single = Solver::new(&domain, &fluid, &config);
    seed(&mut single);
    single.correct(1e-8, &gate).unwrap();
    // The actual ILU(0) first correction is epsilon*[-5/8,1/8,-3/8,-1/16]
    // in the pocket. Its global true residual is ~1.13e-13, but the worst
    // cell retains imbalance 3*epsilon/16 and relative continuity defect 3/43.
    let field = single.field();
    let cell = domain.index(3, 0, 0);
    let magnitude: f64 = (0..6)
        .map(|face| field.outward(&domain, cell, face).abs())
        .sum();
    assert!((field.net_outflow(&domain, cell).abs() / magnitude - 3.0 / 43.0).abs() < 1e-12);

    let mut thermal = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    thermal.faces[0] = ThermalFace::Outflow {
        backflow_temperature: 300.0,
    };
    thermal.faces[1] = ThermalFace::Outflow {
        backflow_temperature: 300.0,
    };
    thermal.cell_sinks.push(CellSink {
        cell: 0,
        conductance_w_k: 1.0,
        temperature_k: 300.0,
    });
    let caloric: Vec<_> = (0..domain.cell_count())
        .map(|cell| {
            domain.is_fluid(cell).then_some(AffineSpecificEnthalpy {
                temperature_k: 300.0,
                enthalpy_j_kg: 0.0,
                heat_capacity_j_kg_k: 2.0,
            })
        })
        .collect();
    let heat = |field: &FlowField| {
        EnergyConfig {
            scheme: ConvectionScheme::Upwind,
            ..EnergyConfig::default()
        }
        .solve_with_enthalpy(
            &domain,
            &fluid,
            &[SolidMaterial::new("partition", 1.0)],
            field,
            &thermal,
            &vec![[1.0; 3]; domain.cell_count()],
            AdvectiveEnthalpyField {
                cells: &caloric,
                boundary_j_kg: [Some(0.0), Some(0.0), None, None, None, None],
            },
            &gate,
        )
    };
    assert!(matches!(
        heat(&field),
        Err(ChtError::InvalidInput {
            field: "thermal.advective_enthalpy",
            ..
        })
    ));

    for pressure_solver in [PressureSolver::IluBicgstab, PressureSolver::AmgCg] {
        config.pressure_solver = pressure_solver;
        let mut solver = Solver::new(&domain, &fluid, &config);
        seed(&mut solver);
        solver.project(&gate).unwrap();
        let velocities = solver.vel.clone();
        let pressure = solver.pressure.clone();
        let iterations = solver.pressure_krylov;
        solver.project(&gate).unwrap();
        assert_eq!(solver.vel, velocities);
        assert_eq!(solver.pressure, pressure);
        assert_eq!(solver.pressure_krylov, iterations);
        let field = solver.field();
        for cell in 0..domain.cell_count() {
            let magnitude: f64 = (0..6)
                .map(|face| field.outward(&domain, cell, face).abs())
                .sum();
            assert!(field.net_outflow(&domain, cell).abs() <= 128.0 * f64::EPSILON * magnitude);
        }
        // Pressure corrects the native faces; it does not clip the quiet
        // circulation to rest to manufacture continuity.
        assert!(field.flux_x(3, 0, 0).abs() > epsilon);
        let energy = heat(&field).unwrap();
        assert!(
            energy
                .temperature
                .iter()
                .all(|temperature| (*temperature - 300.0).abs() < 1e-8)
        );
        assert!(energy.report.balance.residual_w.abs() < 1e-8);
    }
}

#[test]
fn pressure_projection_refuses_incompatible_fixed_flux_without_clipping_it() {
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let mut faces = [FvBoundary::Symmetry; 6];
    faces[0] = FvBoundary::Inlet {
        velocity: [1.0, 0.0, 0.0],
    };
    let config = SimpleConfig::new(faces);
    let mut solver = Solver::new(&domain, &fluid(), &config);
    let original = solver.vel.clone();
    assert!(matches!(
        solver.project(&CancelGate::new()),
        Err(ChtError::SolverNotConverged {
            system: "pressure continuity",
            ..
        })
    ));
    assert_eq!(solver.vel, original);
    let pressure = solver.pressure.clone();
    let gate = CancelGate::new();
    gate.request();
    assert!(matches!(solver.project(&gate), Err(ChtError::Cancelled)));
    assert_eq!(solver.vel, original);
    assert_eq!(solver.pressure, pressure);
}

#[cfg(feature = "enthalpy-coupling")]
#[test]
fn pressure_opening_native_momentum_restart_retains_its_numerical_datum() {
    let (domain, config) = channel(4, 0, 1024.5, 1024.0);
    let fluid = fluid();
    let gate = CancelGate::new();
    let mut solver = Solver::new(&domain, &fluid, &config);
    for _ in 0..5 {
        solver.sweep(&gate).unwrap();
    }
    let state = solver.momentum_state();
    let mut restored = Solver::new(&domain, &fluid, &config);
    restored.restore_momentum(&state).unwrap();
    assert_eq!(restored.pressure, solver.pressure);
    assert_eq!(restored.pressure_reference_pa, solver.pressure_reference_pa);
    assert_eq!(restored.vel, solver.vel);
}
