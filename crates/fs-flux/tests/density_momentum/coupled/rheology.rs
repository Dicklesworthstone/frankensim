//! Native transported-density, material-rheology and dissipation closure.

use super::*;
use fs_flux::ns::boundary::StokesBoundaryError;
use fs_material::MaterialError;
use fs_material::fluid_transport::rheology::CarreauYasuda;
use std::cell::{Cell as Counter, RefCell};

#[path = "rheology/ports.rs"]
mod ports;

fn curve(cell: &Cell<1>) -> CarreauYasuda {
    CarreauYasuda::new(
        0.03 + 0.002 * cell.density_kg_m3 + 0.001 * cell.specific[0],
        0.006,
        4.0,
        0.45,
        2.0,
        [0.0, 100.0],
    )
    .unwrap()
}

#[test]
fn rheology_and_transport_share_the_accepted_endpoint_and_native_energy_balance() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let internal = mesh
        .edges
        .iter()
        .position(|e| e.tris.1 != usize::MAX)
        .unwrap();
    let mut previous = vec![0.0; system.n_u];
    previous[2 * internal + 1] = 0.05;
    let history = previous.clone();
    let original = [
        Cell {
            density_kg_m3: 2.0,
            specific: [3.0],
        },
        Cell {
            density_kg_m3: 8.0,
            specific: [19.0],
        },
    ];
    let inflow = vec![None; mesh.edges.len()];
    let zero = |_: [f64; 2]| [0.0; 2];
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|&e| mesh.is_boundary(e))
        .map(|edge| Boundary::Velocity { edge, value: &zero })
        .collect();
    let native: Vec<_> = (0..mesh.edges.len())
        .filter(|&e| mesh.is_boundary(e))
        .map(|edge| StokesBoundary::Velocity { edge, value: &zero })
        .collect();
    let dt = 0.06;
    let depth = 1.7;
    let mut first_velocity: Option<Vec<f64>> = None;
    for reference in [1.0, 11.0] {
        let calls = Counter::new(0);
        let samples = RefCell::new(vec![(Cell::default(), 0.0, 0.0); original.len()]);
        let source = |cell, point, state: &Cell<1>, rate| {
            assert_eq!(point, mesh.centroids[cell]);
            calls.set(calls.get() + 1);
            let mu = curve(state).viscosity_pa_s(rate)?;
            samples.borrow_mut()[cell] = (*state, rate, mu);
            Ok(mu)
        };
        let solved = system
            .density_momentum_coupled_bdf1_rheology_step(
                numerics(),
                &source,
                &|_, _, _| [0.0; 2],
                &boundaries,
                &previous,
                &original,
                &inflow,
                dt,
                controls(reference, depth),
                nonlinear(),
                || false,
            )
            .unwrap();
        let endpoint = &solved.step;
        let step = &endpoint.step;
        assert!(endpoint.iterations > 1);
        assert!(endpoint.nonlinear_residual <= nonlinear().relative_tolerance);
        assert!(calls.get() >= original.len() * (endpoint.iterations + 1));
        assert_eq!(calls.get() % original.len(), 0);
        assert!(step.cells[0].density_kg_m3 > original[0].density_kg_m3);
        assert!(step.cells[1].density_kg_m3 < original[1].density_kg_m3);
        let plan = system
            .prepare_density_transport(
                &step.flow.field.x[..system.n_u],
                depth,
                1e-12,
                controls(reference, depth).maximum_transport_storage_bytes,
                || false,
            )
            .unwrap();
        let mut workspace = plan.workspace(1 << 22).unwrap();
        let mut transported = vec![Cell::default(); original.len()];
        let report = plan
            .advance_into(
                &original,
                &inflow,
                dt,
                controls(reference, depth).maximum_courant,
                &mut transported,
                &mut workspace,
                || false,
            )
            .unwrap();
        assert_eq!(step.cells, transported);
        assert_eq!(step.transport, report);
        near(report.mass_balance_residual_kg, 0.0, 2e-12);
        near(report.carried_balance_residual[0], 0.0, 2e-11);
        let mut power = 0.0;
        for (cell, output) in solved.cells.iter().enumerate() {
            let rate = system
                .try_cell_kinematics(&step.flow.field.x, cell)
                .unwrap()
                .shear_rate();
            let recorded = samples.borrow()[cell];
            assert_eq!(recorded.0, step.cells[cell]);
            assert_eq!(recorded.1, output.strain_rate_per_s);
            assert_eq!(recorded.2, output.dynamic_viscosity_pa_s);
            near(output.strain_rate_per_s, rate, 1e-14);
            near(output.viscous_heating_w_m3, recorded.2 * rate * rate, 1e-14);
            power += output.viscous_heating_w_m3 * mesh.areas[cell] * depth;
        }
        near(solved.viscous_heating_w, power, 1e-14);
        assert!(power > 0.0);
        // Native SIP work includes the discrete stress coupling. Molecular heat
        // is reported separately and never silently put in the carried slot.
        let viscous = system
            .assemble_stokes_mixed_variable_viscosity(
                FluxParams {
                    sigma: numerics().penalty,
                    ..FluxParams::default()
                },
                &|cell, _| solved.cells[cell].dynamic_viscosity_pa_s / reference,
                &zero,
                &native,
            )
            .unwrap();
        let mut velocity = step.flow.field.x.clone();
        velocity[system.n_u..].fill(0.0);
        let mut action = vec![0.0; system.n];
        viscous.matrix().spmv(&velocity, &mut action);
        let viscous_work = velocity
            .iter()
            .zip(action)
            .map(|(u, au)| u * au)
            .sum::<f64>()
            * reference
            * depth
            * dt;
        let obs = step.observations;
        near(
            obs.kinetic_energy_after_j - obs.kinetic_energy_before_j
                + obs.temporal_dissipation_j
                + obs.interior_upwind_dissipation_j
                + viscous_work,
            0.0,
            2e-9,
        );
        let plateau = system
            .density_momentum_coupled_bdf1_step(
                numerics(),
                &|_, _, state| curve(state).viscosity_pa_s(0.0).unwrap(),
                &|_, _, _| [0.0; 2],
                &boundaries,
                &previous,
                &original,
                &inflow,
                dt,
                controls(reference, depth),
                nonlinear(),
                || false,
            )
            .unwrap();
        assert!(
            step.observations.kinetic_energy_after_j
                > plateau.step.observations.kinetic_energy_after_j
        );
        if let Some(expected) = &first_velocity {
            for (&a, &b) in step.flow.field.x[..system.n_u].iter().zip(expected) {
                near(a, b, 2e-10);
            }
        } else {
            first_velocity = Some(step.flow.field.x[..system.n_u].to_vec());
        }
    }
    assert_eq!(previous, history);
    assert_eq!(
        original[0],
        Cell {
            density_kg_m3: 2.0,
            specific: [3.0]
        }
    );
}

#[test]
fn constant_rheology_preserves_open_reservoir_density_pressure_and_reactions() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 1));
    let system = FluxSystem::new(&mesh);
    let uniform = |_: [f64; 2]| [0.4, -0.2];
    let traction = |_: [f64; 2], normal: [f64; 2]| normal.map(|n| -7.0 * n);
    let selected: Vec<_> = (0..mesh.edges.len())
        .filter(|&e| mesh.is_boundary(e))
        .collect();
    let boundaries: Vec<_> = selected
        .iter()
        .map(|&edge| Boundary::OpenReservoir {
            edge,
            traction: &traction,
            velocity: &uniform,
        })
        .collect();
    let original: Vec<_> = mesh
        .centroids
        .iter()
        .map(|p| Cell {
            density_kg_m3: if p[0] < 0.5 { 2.0 } else { 80.0 },
            specific: [123.0],
        })
        .collect();
    let inflow: Vec<_> = (0..mesh.edges.len())
        .map(|e| {
            mesh.is_boundary(e).then_some(Cell {
                density_kg_m3: 8.0,
                specific: [123.0],
            })
        })
        .collect();
    let previous = moments(&system, &uniform);
    let old = system
        .density_momentum_coupled_bdf1_step(
            numerics(),
            &|_, _, _| 0.2,
            &|_, _, _| [0.0; 2],
            &boundaries,
            &previous,
            &original,
            &inflow,
            0.04,
            controls(11.0, 1.7),
            nonlinear(),
            || false,
        )
        .unwrap();
    let solved = system
        .density_momentum_coupled_bdf1_rheology_step(
            numerics(),
            &|_, _, _, _| Ok(0.2),
            &|_, _, _| [0.0; 2],
            &boundaries,
            &previous,
            &original,
            &inflow,
            0.04,
            controls(11.0, 1.7),
            nonlinear(),
            || false,
        )
        .unwrap();
    let step = &solved.step.step;
    assert_eq!(step.cells, old.step.cells);
    assert_eq!(step.transport, old.step.transport);
    assert_eq!(step.flow.field.x, old.step.flow.field.x);
    assert_eq!(step.pressure_pa, old.step.pressure_pa);
    let actual = step
        .flow
        .reactions
        .resultant(&selected, [0.0; 2], &[None])
        .unwrap();
    let expected = old
        .step
        .flow
        .reactions
        .resultant(&selected, [0.0; 2], &[None])
        .unwrap();
    assert_eq!(actual, expected);
    assert!(solved.viscous_heating_w < 1e-22);
}

#[test]
fn rigid_rotation_uses_symmetric_strain_and_produces_no_molecular_heat() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let omega = 0.2;
    let rotation = |p: [f64; 2]| [-omega * p[1], omega * p[0]];
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|&e| mesh.is_boundary(e))
        .map(|edge| Boundary::Velocity {
            edge,
            value: &rotation,
        })
        .collect();
    let original = vec![
        Cell {
            density_kg_m3: 3.0,
            specific: [42.0]
        };
        mesh.tris.len()
    ];
    let inflow: Vec<_> = (0..mesh.edges.len())
        .map(|e| mesh.is_boundary(e).then_some(original[0]))
        .collect();
    let previous = moments(&system, &rotation);
    let law = CarreauYasuda::new(0.2, 0.01, 4.0, 0.4, 2.0, [0.0, 100.0]).unwrap();
    let solved = system
        .density_momentum_coupled_bdf1_rheology_step(
            numerics(),
            &|_, _, _, rate| law.viscosity_pa_s(rate),
            &|_, p, state| {
                [
                    -state.density_kg_m3 * omega * omega * p[0],
                    -state.density_kg_m3 * omega * omega * p[1],
                ]
            },
            &boundaries,
            &previous,
            &original,
            &inflow,
            0.03,
            controls(7.0, 1.7),
            nonlinear(),
            || false,
        )
        .unwrap();
    let (error, divergence) = system.velocity_error(&solved.step.step.flow.field.x, &rotation);
    assert!(
        error < 2e-10 && divergence < 2e-10,
        "{error:e}, {divergence:e}"
    );
    for output in &solved.cells {
        assert!(output.strain_rate_per_s < 1e-10);
        near(output.dynamic_viscosity_pa_s, 0.2, 1e-12);
    }
    assert!(solved.viscous_heating_w < 1e-18);
}

#[test]
fn original_material_refusal_on_a_later_trial_cancels_the_whole_candidate() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let shear = |p: [f64; 2]| [0.4 * p[1], 0.0];
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|&e| mesh.is_boundary(e))
        .map(|edge| Boundary::Velocity {
            edge,
            value: &shear,
        })
        .collect();
    let original = vec![
        Cell {
            density_kg_m3: 3.0,
            specific: [42.0]
        };
        mesh.tris.len()
    ];
    let inflow: Vec<_> = (0..mesh.edges.len())
        .map(|e| mesh.is_boundary(e).then_some(original[0]))
        .collect();
    let previous = vec![0.0; system.n_u];
    let saved = (previous.clone(), original.clone(), inflow.clone());
    let law = CarreauYasuda::new(0.2, 0.01, 4.0, 0.4, 2.0, [0.0, 0.1]).unwrap();
    let calls = Counter::new(0);
    let source = |_, _, _: &Cell<1>, rate| {
        calls.set(calls.get() + 1);
        law.viscosity_pa_s(rate)
    };
    let error = system
        .density_momentum_coupled_bdf1_rheology_step(
            numerics(),
            &source,
            &|_, _, _| [0.0; 2],
            &boundaries,
            &previous,
            &original,
            &inflow,
            0.03,
            controls(7.0, 1.7),
            nonlinear(),
            || false,
        )
        .err()
        .unwrap();
    assert!(calls.get() > original.len());
    assert!(
        matches!(error, Error::Material { error, .. } if error == law.viscosity_pa_s(0.4).unwrap_err())
    );
    assert_eq!((previous.clone(), original.clone(), inflow.clone()), saved);
    let calibration = MaterialError::Calibration {
        what: "original measured source unavailable".into(),
    };
    let error = system
        .density_momentum_coupled_bdf1_rheology_step(
            numerics(),
            &|_, _, _, _| Err(calibration.clone()),
            &|_, _, _| [0.0; 2],
            &boundaries,
            &previous,
            &original,
            &inflow,
            0.03,
            controls(7.0, 1.7),
            nonlinear(),
            || false,
        )
        .err()
        .unwrap();
    assert_eq!(
        error,
        Error::Material {
            cell: 0,
            error: calibration
        }
    );
    for mu in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        let error = system
            .density_momentum_coupled_bdf1_rheology_step(
                numerics(),
                &|_, _, _, _| Ok(mu),
                &|_, _, _| [0.0; 2],
                &boundaries,
                &previous,
                &original,
                &inflow,
                0.03,
                controls(7.0, 1.7),
                nonlinear(),
                || false,
            )
            .err()
            .unwrap();
        assert!(matches!(
            error,
            Error::Flow(StokesBoundaryError::InvalidViscosity { cell: 0, .. })
        ));
    }
    let queried = Counter::new(false);
    let error = system
        .density_momentum_coupled_bdf1_rheology_step(
            numerics(),
            &|_, _, _, _| {
                queried.set(true);
                Ok(0.2)
            },
            &|_, _, _| [0.0; 2],
            &boundaries,
            &previous,
            &original,
            &inflow,
            0.03,
            controls(7.0, 1.7),
            nonlinear(),
            || queried.get(),
        )
        .err()
        .unwrap();
    assert_eq!(error, Error::Cancelled);
}
