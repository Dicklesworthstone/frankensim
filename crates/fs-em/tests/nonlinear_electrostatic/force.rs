//! Independent held-control shape derivatives and nonlinear exterior stress.
use super::*;
use fs_em::{force::MaxwellForceLimits, transfer::PlanarFrame};
use fs_ga::{Vec3, Wrench};

fn frame() -> PlanarFrame {
    PlanarFrame::new([0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]).unwrap()
}
fn force_limits() -> MaxwellForceLimits {
    MaxwellForceLimits {
        max_cells: 2400,
        max_edges: 4000,
    }
}
fn field_controls() -> NonlinearElectrostaticControls {
    let mut settings = controls();
    settings.linear.max_iterations = 512;
    settings.max_operator_applications = 8192;
    settings
}
fn slab_problem(gap: f64, thickness: f64) -> PlanarNonlinearElectrostatics {
    let original = slab_mesh();
    let mesh = TriComplex2::from_indexed_triangles(
        "perturbed-nonlinear-capacitor",
        original
            .vertices()
            .iter()
            .map(|&[x, y]| [x * gap / GAP, y])
            .collect(),
        original.faces().to_vec(),
        Metric2::planar(thickness).unwrap(),
    )
    .unwrap();
    PlanarNonlinearElectrostatics::new(
        mesh,
        vec![DielectricMaterial::Reversible(nonlinear()); 4],
        limits(),
        &mut || false,
    )
    .unwrap()
}
fn slab_functionals(gap: f64, thickness: f64, charge_controlled: bool) -> [f64; 2] {
    let problem = slab_problem(gap, thickness);
    let electrodes = problem
        .prepare_electrodes(&[&[0, 3], &[2, 5]], &mut || false)
        .unwrap();
    let control = if charge_controlled {
        ElectrodeControl::Charge(2.5 * Q0)
    } else {
        ElectrodeControl::Voltage(2.0 * V0)
    };
    let field = electrodes
        .solve(
            &[0.0; 4],
            &[ElectrodeControl::Voltage(0.0), control],
            &[0.0, V0, 2.0 * V0, 0.0, V0, 2.0 * V0],
            controls(),
            &mut || false,
        )
        .unwrap();
    [field.energy_j(), field.coenergy_j()]
}

#[test]
fn nonlinear_shape_work_matches_distinct_held_voltage_and_charge_functionals() {
    let exact_force = -1.75 * AREA * E0 * D0;
    let gap_increment = GAP * 1e-5;
    let thickness_increment = THICKNESS * 1e-5;
    for charge_controlled in [false, true] {
        let problem = slab_problem(GAP, THICKNESS);
        let electrodes = problem
            .prepare_electrodes(&[&[0, 3], &[2, 5]], &mut || false)
            .unwrap();
        let control = if charge_controlled {
            ElectrodeControl::Charge(2.5 * Q0)
        } else {
            ElectrodeControl::Voltage(2.0 * V0)
        };
        let field = electrodes
            .solve(
                &[0.0; 4],
                &[ElectrodeControl::Voltage(0.0), control],
                &[0.0, V0, 2.0 * V0, 0.0, V0, 2.0 * V0],
                controls(),
                &mut || false,
            )
            .unwrap();
        let displacement: Vec<_> = problem
            .mesh()
            .vertices()
            .iter()
            .map(|&[x, _]| [gap_increment * x / GAP, 0.0])
            .collect();
        let work = field
            .virtual_work(&displacement, 0.0, &mut || false)
            .unwrap();
        near(work / gap_increment, exact_force, 1e-5);
        let forces = field.shape_forces(&mut || false).unwrap();
        let nodal_work: f64 = forces
            .nodal_force_n()
            .iter()
            .zip(&displacement)
            .map(|(f, u)| f[0].mul_add(u[0], f[1] * u[1]))
            .sum();
        near(nodal_work, work, U0 * 1e-5);
        near(forces.thickness_force_n(), 3.25 * U0 / THICKNESS, 1e-5);
        near(
            field
                .virtual_work(&[[0.0; 2]; 6], thickness_increment, &mut || false)
                .unwrap(),
            3.25 * U0 * thickness_increment / THICKNESS,
            U0 * 1e-5,
        );
        near(
            field
                .virtual_work(&[[3e-4, -2e-4]; 6], 0.0, &mut || false)
                .unwrap(),
            0.0,
            U0,
        );

        // Solve the perturbed geometry independently. At fixed Q mechanical
        // work is -dW; at fixed V it is +dW*, including the source's work.
        let functional = |value: [f64; 2]| {
            if charge_controlled {
                -value[0]
            } else {
                value[1]
            }
        };
        let gap_difference = (functional(slab_functionals(
            GAP + gap_increment,
            THICKNESS,
            charge_controlled,
        )) - functional(slab_functionals(
            GAP - gap_increment,
            THICKNESS,
            charge_controlled,
        ))) / (2.0 * gap_increment);
        assert!((gap_difference - exact_force).abs() < 2e-8 * exact_force.abs());
        let depth_difference = (functional(slab_functionals(
            GAP,
            THICKNESS + thickness_increment,
            charge_controlled,
        )) - functional(slab_functionals(
            GAP,
            THICKNESS - thickness_increment,
            charge_controlled,
        ))) / (2.0 * thickness_increment);
        assert!((depth_difference - forces.thickness_force_n()).abs() < 2e-8 * 1e-5);

        let mut polls = 0;
        let baseline = field
            .virtual_work(&displacement, thickness_increment, &mut || {
                polls += 1;
                false
            })
            .unwrap();
        for stop in [1, polls / 2, polls] {
            let mut seen = 0;
            assert!(
                field
                    .virtual_work(&displacement, thickness_increment, &mut || {
                        seen += 1;
                        seen >= stop
                    })
                    .is_err()
            );
            assert_eq!(seen, stop);
        }
        assert_eq!(
            field
                .virtual_work(&displacement, thickness_increment, &mut || false)
                .unwrap(),
            baseline
        );
        assert!(field.virtual_work(&[], 0.0, &mut || false).is_err());
    }
    // Using -dW while holding voltage fixed produces the wrong sign here.
    let wrong = -(slab_functionals(GAP + gap_increment, THICKNESS, false)[0]
        - slab_functionals(GAP - gap_increment, THICKNESS, false)[0])
        / (2.0 * gap_increment);
    assert!((wrong - 0.25 * AREA * E0 * D0).abs() < 2e-8 * 1e-5);
    assert!(wrong > 0.0 && exact_force < 0.0);
}

#[test]
fn shape_work_transports_original_cell_charge_through_the_held_control_functional() {
    let density = 1e-3;
    let increment = GAP * 1e-5;
    for charge_controlled in [false, true] {
        let control = if charge_controlled {
            ElectrodeControl::Charge(2.0 * Q0)
        } else {
            ElectrodeControl::Voltage(2.0 * V0)
        };
        let problem = slab_problem(GAP, THICKNESS);
        let electrodes = problem
            .prepare_electrodes(&[&[0, 3], &[2, 5]], &mut || false)
            .unwrap();
        let field = electrodes
            .solve(
                &[density; 4],
                &[ElectrodeControl::Voltage(0.0), control],
                &[0.0, V0, 2.0 * V0, 0.0, V0, 2.0 * V0],
                controls(),
                &mut || false,
            )
            .unwrap();
        let displacement: Vec<_> = problem
            .mesh()
            .vertices()
            .iter()
            .map(|&[x, _]| [increment * x / GAP, 0.0])
            .collect();
        let work = field
            .virtual_work(&displacement, 0.0, &mut || false)
            .unwrap();
        let functional = |gap| {
            let perturbed = slab_problem(gap, THICKNESS);
            let terminals = perturbed
                .prepare_electrodes(&[&[0, 3], &[2, 5]], &mut || false)
                .unwrap();
            // Each original triangle retains its total material charge.
            // Its density changes inversely with the affine volume change.
            let state = terminals
                .solve(
                    &[density * GAP / gap; 4],
                    &[ElectrodeControl::Voltage(0.0), control],
                    field.relative_voltage_v(),
                    controls(),
                    &mut || false,
                )
                .unwrap();
            near(state.volume_charge_c(), density * AREA * GAP, Q0);
            let voltage_source_work = if charge_controlled {
                0.0
            } else {
                2.0 * V0 * state.terminal_charge_c()[1]
            };
            -state.energy_j() + voltage_source_work
        };
        let derivative =
            (functional(GAP + increment) - functional(GAP - increment)) / (2.0 * increment);
        assert!(
            (work / increment - derivative).abs() < 2e-8 * 1e-5,
            "transported cell charge: analytic={}, independent={derivative:e}, charge control={charge_controlled}",
            work / increment,
        );
    }
}

struct Sector {
    mesh: TriComplex2,
    inner: Vec<usize>,
    outer: Vec<usize>,
    enclosed: Vec<usize>,
    charged: Vec<usize>,
    initial: Vec<f64>,
}
fn radial_potential(radius: f64) -> f64 {
    V0 * (3.0 * radius - 5.0 * radius.ln() - 3.0)
}
fn sector(n: usize, reverse: bool) -> Sector {
    let mut vertices = Vec::new();
    let mut initial = Vec::new();
    for radial in 0..=n {
        let radius = 0.9 + 0.3 * radial as f64 / n as f64;
        for angular in 0..=n {
            let angle = -0.3 + 0.6 * angular as f64 / n as f64;
            vertices.push([GAP * radius * angle.cos(), GAP * radius * angle.sin()]);
            initial.push(radial_potential(radius) - radial_potential(0.9));
        }
    }
    let mut faces = Vec::new();
    let mut enclosed = Vec::new();
    let mut charged = Vec::new();
    for radial in 0..n {
        for angular in 0..n {
            let a = (radial * (n + 1) + angular) as u32;
            let b = a + 1;
            let d = a + (n + 1) as u32;
            let c = d + 1;
            if radial >= n / 3 && radial < 2 * n / 3 && angular >= n / 3 && angular < 2 * n / 3 {
                enclosed.extend([faces.len(), faces.len() + 1]);
                if radial > n / 3
                    && radial + 1 < 2 * n / 3
                    && angular > n / 3
                    && angular + 1 < 2 * n / 3
                {
                    charged.extend([faces.len(), faces.len() + 1]);
                }
            }
            faces.extend(if reverse {
                [[a, c, d], [a, b, c]]
            } else {
                [[a, d, c], [a, c, b]]
            });
        }
    }
    Sector {
        mesh: TriComplex2::from_indexed_triangles(
            "radial-dielectric-sector",
            vertices,
            faces,
            Metric2::planar(THICKNESS).unwrap(),
        )
        .unwrap(),
        inner: (0..=n).collect(),
        outer: (n * (n + 1)..=(n + 1) * (n + 1) - 1).collect(),
        enclosed,
        charged,
        initial,
    }
}
fn sector_limits() -> WeightedAssemblyLimits {
    WeightedAssemblyLimits {
        max_cells: 2400,
        max_dofs: 4000,
        max_triplets: 21600,
    }
}

#[test]
fn source_free_nonlinear_radial_stress_converges_and_falsifies_half_field_pairing() {
    // On the middle source branch E/E0=5/r-3 and D/D0=2.5/r solve div D=0.
    // The internal r=[1,1.1], theta=[-0.1,0.1] contour encloses no charge.
    // Replacing coenergy pressure by E.D/2 would instead give this nonzero Fx.
    let wrong_half_force = -7.5 * E0 * D0 * THICKNESS * GAP * 0.1_f64.sin() * 1.1_f64.ln();
    let mut previous = f64::INFINITY;
    for n in [6, 12, 24] {
        let fixture = sector(n, false);
        let cell_count = fixture.mesh.faces().len();
        let problem = PlanarNonlinearElectrostatics::new(
            fixture.mesh,
            vec![DielectricMaterial::Reversible(nonlinear()); cell_count],
            sector_limits(),
            &mut || false,
        )
        .unwrap();
        let electrodes = problem
            .prepare_electrodes(&[&fixture.inner, &fixture.outer], &mut || false)
            .unwrap();
        let field = electrodes
            .solve(
                &vec![0.0; cell_count],
                &[
                    ElectrodeControl::Voltage(0.0),
                    ElectrodeControl::Voltage(radial_potential(1.2) - radial_potential(0.9)),
                ],
                &fixture.initial,
                field_controls(),
                &mut || false,
            )
            .unwrap();
        let load = field
            .maxwell_wrench(
                &fixture.enclosed,
                &frame(),
                [0.0, 0.0, THICKNESS / 2.0],
                force_limits(),
                &mut || false,
            )
            .unwrap();
        let error = load.force.x.hypot(load.force.y);
        assert!(
            error < 0.75 * previous,
            "n={n}, force error={error:e}, previous={previous:e}, half-pairing force={wrong_half_force:e}"
        );
        previous = error;
        let mut reordered = fixture.enclosed.clone();
        reordered.reverse();
        assert_eq!(
            load,
            field
                .maxwell_wrench(
                    &reordered,
                    &frame(),
                    [0.0, 0.0, THICKNESS / 2.0],
                    force_limits(),
                    &mut || false
                )
                .unwrap()
        );
    }
    assert!(
        previous < 0.15 * wrong_half_force.abs(),
        "correct stress error={previous:e}, wrong half-pairing force={wrong_half_force:e}"
    );
}

#[test]
fn nonlinear_field_owner_preserves_original_linear_force_frames_and_cancellation() {
    let fixture = sector(12, true);
    let source = linear(2.0 * EPSILON0);
    let cell_count = fixture.mesh.faces().len();
    let original = PlanarElectrostatics::new(
        fixture.mesh.clone(),
        vec![source.clone(); cell_count],
        sector_limits(),
        &mut || false,
    )
    .unwrap();
    let problem = PlanarNonlinearElectrostatics::new(
        fixture.mesh,
        vec![DielectricMaterial::isotropic(source); cell_count],
        sector_limits(),
        &mut || false,
    )
    .unwrap();
    let mut rho = vec![0.0; cell_count];
    for cell in fixture.charged {
        rho[cell] = 1e-3;
    }
    let electrodes = problem
        .prepare_electrodes(&[&fixture.inner, &fixture.outer], &mut || false)
        .unwrap();
    let potential = radial_potential(1.2) - radial_potential(0.9);
    let field = electrodes
        .solve(
            &rho,
            &[
                ElectrodeControl::Voltage(0.0),
                ElectrodeControl::Voltage(potential),
            ],
            &fixture.initial,
            field_controls(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(field.free_charge_density_c_m3(), rho);
    let boundary: Vec<_> = fixture
        .inner
        .iter()
        .map(|&v| (v, 0.0))
        .chain(fixture.outer.iter().map(|&v| (v, potential)))
        .collect();
    let old = original
        .solve(
            &rho,
            &boundary,
            ElectrostaticControls {
                linear: field_controls().linear,
                charge_tolerance_c: 1e-20,
                energy_tolerance_j: 1e-16,
            },
            &mut || false,
        )
        .unwrap();
    let origin = [0.0, 0.0, THICKNESS / 2.0];
    let old_load = old
        .maxwell_wrench(
            &fixture.enclosed,
            &frame(),
            origin,
            force_limits(),
            &mut || false,
        )
        .unwrap();
    let mut polls = 0;
    let load = field
        .maxwell_wrench(
            &fixture.enclosed,
            &frame(),
            origin,
            force_limits(),
            &mut || {
                polls += 1;
                false
            },
        )
        .unwrap();
    for (actual, expected) in load.to_array().into_iter().zip(old_load.to_array()) {
        near(actual, expected, 1e-8);
    }
    let placed = PlanarFrame::new([2.0, 3.0, 4.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]).unwrap();
    let global = field
        .maxwell_wrench(
            &fixture.enclosed,
            &placed,
            [2.0 + THICKNESS / 2.0 + 0.2, 3.3, 4.4],
            force_limits(),
            &mut || false,
        )
        .unwrap();
    let force = Vec3::new(load.force.z, load.force.x, load.force.y);
    let torque = Vec3::new(load.torque.z, load.torque.x, load.torque.y)
        + Vec3::new(-0.2, -0.3, -0.4).cross(force);
    for (actual, expected) in global
        .to_array()
        .into_iter()
        .zip(Wrench::new(torque, force).to_array())
    {
        near(actual, expected, 1e-8);
    }
    for stop in [1, polls / 2, polls] {
        let mut seen = 0;
        assert!(
            field
                .maxwell_wrench(
                    &fixture.enclosed,
                    &frame(),
                    origin,
                    force_limits(),
                    &mut || {
                        seen += 1;
                        seen >= stop
                    }
                )
                .is_err()
        );
        assert_eq!(seen, stop);
    }
    assert_eq!(
        load,
        field
            .maxwell_wrench(
                &fixture.enclosed,
                &frame(),
                origin,
                force_limits(),
                &mut || false
            )
            .unwrap()
    );
    assert!(
        field
            .maxwell_wrench(
                &fixture.enclosed,
                &frame(),
                origin,
                MaxwellForceLimits {
                    max_cells: 1,
                    max_edges: 1
                },
                &mut || false
            )
            .is_err()
    );
}
