use super::*;
use crate::compressible::mixture as planar;
use crate::compressible::mixture::tests::{close, gas, source, standard};
use crate::compressible::{AxisBoundaries, Boundaries2d, CartesianGrid as Grid2d};
use fs_gas::multicomponent::{Boundaries, Boundary, MixtureEuler1d};
use fs_qty::Length;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

const PERIODIC: Boundaries3d = Boundaries3d {
    x: AxisBoundaries::Periodic,
    y: AxisBoundaries::Periodic,
    z: AxisBoundaries::Periodic,
};
const WALLS: AxisBoundaries = AxisBoundaries::Ends {
    lower: BoundaryCondition::Reflecting,
    upper: BoundaryCondition::Reflecting,
};
const OPEN: AxisBoundaries = AxisBoundaries::Ends {
    lower: BoundaryCondition::Transmissive,
    upper: BoundaryCondition::Transmissive,
};

fn field(
    gas: MixtureGas,
    grid: CartesianGrid,
    boundaries: Boundaries3d,
    initial: impl Fn([f64; 3]) -> ([f64; 3], [f64; 3], f64),
) -> MixtureEuler3d {
    let count: usize = grid.shape().iter().product();
    let mut states = vec![0.0; count * 7];
    for (index, row) in states.chunks_exact_mut(7).enumerate() {
        let (densities, velocity, pressure) = initial(grid.cell_center(index).unwrap());
        state_into(&gas, &densities, velocity, pressure, row).unwrap();
    }
    let bytes = MixtureEuler3d::required_workspace_bytes(&gas, &grid).unwrap();
    MixtureEuler3d::new(gas, grid, states, boundaries, bytes, || false).unwrap()
}

fn balance(report: StepReport) {
    assert!(report.maximum_cfl <= 0.40000000000001);
    for residual in report.conservation_residual.species_masses {
        close(residual, 0.0, 3e-12);
    }
    for residual in report.conservation_residual.momentum {
        close(residual, 0.0, 3e-12);
    }
    close(report.conservation_residual.energy, 0.0, 3e-12);
    close(report.absolute_energy_residual_j, 0.0, 3e-10);
}

fn run(flow: &mut MixtureEuler3d, end: f64) {
    let mut steps = 0;
    while flow.time().value() < end {
        balance(
            flow.step(
                Time::new(end - flow.time().value()),
                StepControl::default(),
                || false,
            )
            .unwrap(),
        );
        steps += 1;
        assert!(steps < 500);
    }
}

#[test]
fn g1_three_axis_species_contact_refines_with_original_local_calorics() {
    // Equal molar Cp/R gives common gamma, but distinct molecular masses make
    // the evolving temperature composition-dependent. The exact solution is
    // a smooth material contact, averaged independently over physical cells.
    let end = 0.1;
    let velocity = [0.4, 0.25, -0.15];
    let mut errors = Vec::new();
    for n in [12, 24] {
        let dx = 1.0 / n as f64;
        let z = std::f64::consts::PI * dx;
        let average = (z.sin() / z).powi(3);
        let exact = |x: [f64; 3], time: f64| {
            0.55 + 0.2
                * average
                * (2.0 * std::f64::consts::PI * (x[0] + x[1] + x[2] - 0.5 * time)).sin()
        };
        let grid = CartesianGrid::new(n, n, n, dx, dx, dx).unwrap();
        let mut flow = field(standard(0.0), grid, PERIODIC, |x| {
            let first = exact(x, 0.0);
            ([first, 1.0 - first, 0.0], velocity, 1.0)
        });
        let before = flow.totals();
        run(&mut flow, end);
        let mut error = 0.0;
        for index in 0..flow.cell_count() {
            let q = flow.cell(index).unwrap();
            let p = q.primitive(flow.gas(), || false).unwrap();
            let [a, b, absent] = <[f64; 3]>::try_from(q.partial_densities).unwrap();
            assert_eq!(absent, 0.0);
            close(p.density, 1.0, 3e-11);
            close(p.pressure, 1.0, 3e-10);
            close(p.temperature, 1.0 / (a + 0.5 * b), 3e-10);
            close(p.sound_speed.powi(2), 1.4, 3e-10);
            for axis in 0..3 {
                close(p.velocity[axis], velocity[axis], 3e-10);
            }
            error += (a - exact(grid.cell_center(index).unwrap(), end)).abs();
        }
        errors.push(error / flow.cell_count() as f64);
        let after = flow.totals();
        for species in 0..3 {
            close(
                after.species_masses[species],
                before.species_masses[species],
                3e-12,
            );
        }
        close(
            after.absolute_energy(flow.gas()).unwrap(),
            before.absolute_energy(flow.gas()).unwrap(),
            3e-12,
        );
    }
    let order = (errors[0] / errors[1]).ln() / 2.0_f64.ln();
    eprintln!("3D species contact L1 on 12^3/24^3: {errors:?}; order={order}");
    assert!(order > 1.5, "{errors:?}");
    assert!(errors[1] < 1.5e-3);
}

#[test]
fn g0_compatible_species_transport_preserves_uniform_inert_fraction() {
    // All three sources have gamma=1.4; the distinct molecular masses still
    // change the exact original temperature with composition. Constant p,u and
    // independently advected rho,Y are therefore an exact material-contact
    // solution. An inert constant mass fraction must remain constant even when
    // density and the other fractions have different spatial gradients.
    let gas = gas(vec![
        source("A", R, 3.5, 0.0, 30.0),
        source("B", 2.0 * R, 3.5, 0.0, -60.0),
        source("C", 0.5 * R, 3.5, 0.0, 7.0),
    ]);
    let n = 8;
    let dx = 1.0 / n as f64;
    let velocity = [0.4, -0.25, 0.15];
    let sinc = |wave: f64| {
        let x = std::f64::consts::PI * wave * dx;
        if x == 0.0 { 1.0 } else { x.sin() / x }
    };
    let mode = |x: [f64; 3], time: f64, wave: [f64; 3]| {
        let phase = 2.0
            * std::f64::consts::PI
            * (0..3)
                .map(|a| wave[a] * (x[a] - velocity[a] * time))
                .sum::<f64>();
        let average = wave.into_iter().map(sinc).product::<f64>();
        (average * phase.sin(), average * phase.cos())
    };
    let exact = |x: [f64; 3], time: f64| {
        let density = 1.0 + 0.3 * mode(x, time, [1.0, 1.0, 1.0]).0;
        // Average rho*Y_A, not the product of separate rho and Y averages.
        let a = 0.4 * density
            + 0.2 * mode(x, time, [1.0, -1.0, 1.0]).1
            + 0.03 * (mode(x, time, [2.0, 0.0, 2.0]).0 + mode(x, time, [0.0, 2.0, 0.0]).0);
        let c = 0.2 * density;
        [a, density - a - c, c]
    };
    let grid = CartesianGrid::new(n, n, n, dx, dx, dx).unwrap();
    let mut flow = field(gas, grid, PERIODIC, |x| (exact(x, 0.0), velocity, 1.0));
    let before = flow.totals();
    for _ in 0..5 {
        balance(
            flow.step(Time::new(0.003), StepControl::default(), || false)
                .unwrap(),
        );
    }
    let mut maximum_inert_error: f64 = 0.0;
    let mut species_error = 0.0;
    for i in 0..flow.cell_count() {
        let q = flow.cell(i).unwrap();
        let p = q.primitive(flow.gas(), || false).unwrap();
        maximum_inert_error =
            maximum_inert_error.max((q.partial_densities[2] / p.density - 0.2).abs());
        let expected = exact(grid.cell_center(i).unwrap(), flow.time().value());
        for k in 0..3 {
            species_error += (q.partial_densities[k] - expected[k]).abs();
            close(p.velocity[k], velocity[k], 3e-11);
        }
        close(p.pressure, 1.0, 3e-11);
        close(
            p.temperature,
            1.0 / (q.partial_densities[0]
                + 0.5 * q.partial_densities[1]
                + 2.0 * q.partial_densities[2]),
            3e-11,
        );
        close(p.sound_speed.powi(2), 1.4 / p.density, 3e-11);
    }
    let after = flow.totals();
    for k in 0..3 {
        close(after.species_masses[k], before.species_masses[k], 3e-12);
    }
    close(
        after.absolute_energy(flow.gas()).unwrap(),
        before.absolute_energy(flow.gas()).unwrap(),
        3e-12,
    );
    let species_error = species_error / flow.cell_count() as f64;
    eprintln!(
        "compatible 3-D species transport: maximum inert-Y error={maximum_inert_error:e}, total-species L1={species_error:e}"
    );
    assert!(species_error < 0.01);
    assert!(
        maximum_inert_error < 5e-13,
        "uniform inert mass fraction drifted by {maximum_inert_error:e}"
    );
}

#[test]
fn g2_variable_caloric_shocks_on_every_axis_match_original_one_dimensional_flow() {
    let n = 24;
    let dx = 1.0 / n as f64;
    let initial = |position: f64| {
        if position < 0.5 {
            ([0.9, 0.1, 0.0], 1.0)
        } else {
            ([0.025, 0.1, 0.0], 0.1)
        }
    };
    for axis in 0..3 {
        let mut shape = [2, 3, 2];
        let mut spacing = [0.2, 0.15, 0.3];
        let mut boundaries = [AxisBoundaries::Periodic; 3];
        shape[axis] = n;
        spacing[axis] = dx;
        boundaries[axis] = OPEN;
        let grid = CartesianGrid::new(
            shape[0], shape[1], shape[2], spacing[0], spacing[1], spacing[2],
        )
        .unwrap();
        let gas = standard(0.35);
        let mut flow = field(
            gas.clone(),
            grid,
            Boundaries3d {
                x: boundaries[0],
                y: boundaries[1],
                z: boundaries[2],
            },
            |x| {
                let (rho, p) = initial(x[axis]);
                (rho, [0.0; 3], p)
            },
        );
        let states = (0..n)
            .map(|i| {
                let (rho, p) = initial((i as f64 + 0.5) * dx);
                gas.conserved(&rho, 0.0, p).unwrap()
            })
            .collect();
        let area = (0..3)
            .filter(|&d| d != axis)
            .map(|d| shape[d] as f64 * spacing[d])
            .product();
        let mut axial = MixtureEuler1d::new(
            gas,
            Length::new(dx),
            area,
            states,
            Boundaries::Ends {
                left: Boundary::Transmissive,
                right: Boundary::Transmissive,
            },
            || false,
        )
        .unwrap();
        while flow.time().value() < 0.06 {
            let report = flow
                .step(
                    Time::new(0.06 - flow.time().value()),
                    StepControl::default(),
                    || false,
                )
                .unwrap();
            balance(report);
            let reference = axial
                .step(report.elapsed, StepControl::default(), || false)
                .unwrap();
            close(reference.elapsed.value(), report.elapsed.value(), 1e-14);
        }
        for index in 0..flow.cell_count() {
            let coordinates = [
                index % shape[0],
                (index / shape[0]) % shape[1],
                index / (shape[0] * shape[1]),
            ];
            let q = flow.cell(index).unwrap();
            let reference = axial.cell(coordinates[axis]).unwrap();
            for species in 0..3 {
                close(
                    q.partial_densities[species],
                    reference.partial_densities[species],
                    5e-9,
                );
            }
            for d in 0..3 {
                close(
                    q.momentum[d],
                    if d == axis { reference.momentum } else { 0.0 },
                    5e-9,
                );
            }
            close(q.energy, reference.energy, 5e-9);
        }
    }
}

#[test]
fn g3_nonuniform_planar_mixture_extrusion_preserves_the_existing_owner() {
    let gas = standard(0.35);
    let (nx, ny, nz) = (9, 7, 3);
    let (dx, dy, dz) = (1.0 / nx as f64, 1.0 / ny as f64, 0.2);
    let grid = CartesianGrid::new(nx, ny, nz, dx, dy, dz).unwrap();
    let initial = |x: [f64; 3]| {
        let phase = 2.0 * std::f64::consts::PI;
        let a = 0.5 + 0.1 * (phase * x[0]).sin() * (phase * x[1]).cos();
        (
            [a, 1.0 - a, 0.0],
            [
                0.2 + 0.03 * (phase * x[1]).sin(),
                -0.1 + 0.02 * (phase * x[0]).cos(),
                0.0,
            ],
            1.0 + 0.08 * (phase * (x[0] + x[1])).cos(),
        )
    };
    let mut volume = field(gas.clone(), grid, PERIODIC, initial);
    let plane = Grid2d::new(nx, ny, dx, dy, nz as f64 * dz).unwrap();
    let mut states = vec![0.0; nx * ny * 6];
    for (i, row) in states.chunks_exact_mut(6).enumerate() {
        let [x, y] = plane.cell_center(i).unwrap();
        let (rho, velocity, pressure) = initial([x, y, 0.0]);
        planar::state_into(&gas, &rho, [velocity[0], velocity[1]], pressure, row).unwrap();
    }
    let bytes = planar::MixtureEuler2d::required_workspace_bytes(&gas, &plane).unwrap();
    let mut flat = planar::MixtureEuler2d::new(
        gas,
        plane,
        states,
        Boundaries2d {
            x: AxisBoundaries::Periodic,
            y: AxisBoundaries::Periodic,
        },
        bytes,
        || false,
    )
    .unwrap();
    for _ in 0..4 {
        let report = volume
            .step(Time::new(0.003), StepControl::default(), || false)
            .unwrap();
        balance(report);
        let reference = flat
            .step(report.elapsed, StepControl::default(), || false)
            .unwrap();
        close(reference.elapsed.value(), report.elapsed.value(), 1e-14);
    }
    for index in 0..volume.cell_count() {
        let q = volume.cell(index).unwrap();
        let r = flat.cell(index % (nx * ny)).unwrap();
        for species in 0..3 {
            close(
                q.partial_densities[species],
                r.partial_densities[species],
                3e-10,
            );
        }
        close(q.momentum[0], r.momentum[0], 3e-10);
        close(q.momentum[1], r.momentum[1], 3e-10);
        close(q.momentum[2], 0.0, 3e-12);
        close(q.energy, r.energy, 3e-10);
    }
}

#[test]
fn g2_six_wall_mixture_impulses_close_species_full_momentum_and_chemical_energy() {
    let grid = CartesianGrid::new(8, 6, 5, 0.125, 0.13, 0.17).unwrap();
    let mut flow = field(
        standard(0.35),
        grid,
        Boundaries3d {
            x: WALLS,
            y: WALLS,
            z: WALLS,
        },
        |[x, y, z]| {
            let a = 0.3 + 0.2 * x + 0.1 * y + 0.1 * z;
            let hot = (x - 0.3).powi(2) + (y - 0.4).powi(2) + (z - 0.5).powi(2) < 0.2_f64.powi(2);
            (
                [a, 1.0 - a, 0.0],
                [1.2, -0.8, 0.6],
                if hot { 2.0 } else { 0.8 },
            )
        },
    );
    let initial = flow.totals();
    let mut outward = [0.0; 3];
    while flow.time().value() < 0.03 {
        let report = flow
            .step(
                Time::new(0.03 - flow.time().value()),
                StepControl::default(),
                || false,
            )
            .unwrap();
        balance(report);
        assert_eq!(report.boundary_change.energy, 0.0);
        assert!(
            report
                .boundary_change
                .species_masses
                .iter()
                .all(|&q| q == 0.0)
        );
        for side in BoundarySide::ALL {
            assert_eq!(flow.boundary_impulses(side).len(), grid.face_count(side));
            for face in flow.boundary_impulses(side) {
                assert_eq!(face.species_masses, &[0.0; 3]);
                assert_eq!(face.energy, 0.0);
                for axis in 0..3 {
                    if axis != side.axis() {
                        assert_eq!(face.momentum[axis], 0.0);
                    }
                    outward[axis] += face.momentum[axis];
                }
            }
        }
        for index in 0..flow.cell_count() {
            let q = flow.cell(index).unwrap();
            assert!(q.partial_densities.iter().all(|&rho| rho >= 0.0));
            assert_eq!(q.partial_densities[2], 0.0);
            assert!(q.primitive(flow.gas(), || false).unwrap().pressure > 0.0);
        }
    }
    let after = flow.totals();
    for axis in 0..3 {
        close(
            after.momentum[axis] + outward[axis],
            initial.momentum[axis],
            3e-12,
        );
    }
    for species in 0..3 {
        close(
            after.species_masses[species],
            initial.species_masses[species],
            3e-12,
        );
    }
    close(after.energy, initial.energy, 3e-12);
    close(
        after.absolute_energy(flow.gas()).unwrap(),
        initial.absolute_energy(flow.gas()).unwrap(),
        3e-12,
    );
    assert!(after.momentum[0] < initial.momentum[0]);
    assert!(after.momentum[1] > initial.momentum[1]);
    assert!(after.momentum[2] < initial.momentum[2]);
}

#[test]
fn g4_complete_workspace_budget_and_last_checkpoint_cancellation_are_atomic() {
    let grid = CartesianGrid::new(5, 4, 3, 0.2, 0.15, 0.25).unwrap();
    let make = || {
        field(standard(0.2), grid, PERIODIC, |[x, y, z]| {
            let a = 0.45 + 0.1 * x.sin() * y.cos() * z.sin();
            ([a, 1.0 - a, 0.0], [0.3, -0.2, 0.1], 1.0 + 0.02 * x.cos())
        })
    };
    let mut flow = make();
    flow.step(Time::new(0.001), StepControl::default(), || false)
        .unwrap();
    let bytes = flow.workspace_bytes();
    let before = flow.states().to_vec();
    assert!(matches!(
        MixtureEuler3d::new(
            flow.gas().clone(),
            grid,
            before.clone(),
            PERIODIC,
            bytes - 1,
            || false
        ),
        Err(Error::Invalid(_))
    ));
    let mut excess_capacity = Vec::with_capacity(before.len() + 64);
    excess_capacity.extend_from_slice(&before);
    assert!(matches!(
        MixtureEuler3d::new(
            flow.gas().clone(),
            grid,
            excess_capacity,
            PERIODIC,
            bytes,
            || false
        ),
        Err(Error::Invalid(_))
    ));
    let mut negative_species = before.clone();
    negative_species[2] = -1e-12;
    assert!(
        MixtureEuler3d::new(
            flow.gas().clone(),
            grid,
            negative_species,
            PERIODIC,
            bytes,
            || false
        )
        .is_err()
    );
    let mut probe = make();
    probe
        .step(Time::new(0.001), StepControl::default(), || false)
        .unwrap();
    let mut polls = 0;
    let expected = probe
        .step(Time::new(0.002), StepControl::default(), || {
            polls += 1;
            false
        })
        .unwrap();
    let time = flow.time();
    let receipt = *flow.last_step().unwrap();
    let loads = flow.boundary_impulses.clone();
    let mut seen = 0;
    assert!(matches!(
        flow.step(Time::new(0.002), StepControl::default(), || {
            seen += 1;
            seen == polls
        }),
        Err(Error::Cancelled)
    ));
    assert_eq!(seen, polls);
    assert_eq!(flow.states(), before);
    assert_eq!(flow.time(), time);
    assert_eq!(flow.last_step(), Some(&receipt));
    assert_eq!(flow.boundary_impulses, loads);
    let replay = flow
        .step(Time::new(0.002), StepControl::default(), || false)
        .unwrap();
    assert_eq!(replay, expected);
    assert_eq!(flow.states(), probe.states());
    assert_eq!(flow.boundary_impulses, probe.boundary_impulses);
}

#[test]
fn g4_both_transverse_shear_energies_participate_in_original_source_admission() {
    let grid = CartesianGrid::new(2, 2, 2, 0.5, 0.5, 0.5).unwrap();
    let make = |second: f64| {
        field(standard(0.0), grid, PERIODIC, |[x, _, _]| {
            if x < 0.5 {
                ([1.0, 0.0, 0.0], [-2.0, 12.0, second], 0.4)
            } else {
                ([0.0, 1.0, 0.0], [2.0, -12.0, -second], 0.4)
            }
        })
    };
    let controls = StepControl {
        reconstruction: Reconstruction::FirstOrder,
        maximum_halvings: 2,
        ..StepControl::default()
    };
    let mut one_shear = make(0.0);
    one_shear.step(Time::new(1e-4), controls, || false).unwrap();
    let mut both = make(12.0);
    let before = both.states().to_vec();
    let loads = both.boundary_impulses.clone();
    let error = both.step(Time::new(1e-4), controls, || false).unwrap_err();
    assert!(
        matches!(error, Error::StepRetries { attempts: 3, last } if matches!(*last, Error::Invalid("mixture HLL intermediate state outside source domain")))
    );
    assert_eq!(both.states(), before);
    assert_eq!(both.boundary_impulses, loads);
    assert_eq!(both.time().value(), 0.0);
    assert_eq!(both.last_step(), None);
}
