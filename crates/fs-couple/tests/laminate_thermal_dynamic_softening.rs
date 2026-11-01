//! G1/G3 changing-source inertia on the original shell and finite heat owners.
//! Authored cards and independent numerical references, not measured validation.

mod laminate_thermal_common;
use laminate_thermal_common::*;

use fs_couple::laminate_thermal::reciprocal::{
    LaminateThermoelasticControls, LaminateThermoelasticError, ReciprocalLayeredThermalShell,
    dynamic::DynamicLaminateState,
};
use fs_couple::laminate_thermal::{LaminateShellDiscretization, LayeredThermalShell};
use fs_material::laminate::ReferenceSurface;
use fs_material::laminate::thermal::source::ResolvedLaminateThermal;
use fs_material::thermal::transport::stack::transient::implicit::{
    LayeredThermalSolveControls, LayeredThermalWorkspace,
};
use fs_plate::shell::ShellMesh;
use fs_plate::shell::mitc::strain::{PreparedElementStrain, ShellStrainPoint};
use fs_plate::shell::mitc::thermal::LaminateThermalField;
use fs_plate::shell::mitc::{ConstraintFrame, MitcBudget, NodalDirectors, ShellConstraint};

const T0: f64 = 300.0;
const ALPHA: f64 = 1e-5;
const ALPHA_SLOPE: f64 = 1e-7;
const E_SLOPE: f64 = -0.002 * E;

fn controls() -> LaminateThermoelasticControls {
    LaminateThermoelasticControls {
        maximum_iterations: 64,
        heat_relaxation: 1.0,
        cell_heat_tolerance_j: 1e-8,
        energy_tolerance_j: 1e-7,
    }
}

fn field_internal(state: &DynamicLaminateState) -> f64 {
    state
        .cell_potentials()
        .iter()
        .map(|p| p.internal_energy_j)
        .sum::<f64>()
        + state.transverse_shear_internal_energy_j()
        + state.drilling_gauge_energy_j()
        + state.kinetic_energy_j()
}

/// Reflection in y preserves equal x displacements at the two free nodes.
/// The mode is u_x=q*x, with exact k=E*V and m=rho*V/2. All other nodal
/// components are supported; a symmetric single layer has no bubble forcing.
fn axial_triangle() -> LaminateShellDiscretization {
    let mut d = discretization();
    d.mesh = ShellMesh::new(
        vec![[0.0, 0.0, 0.0], [1.0, -0.6, 0.0], [1.0, 0.6, 0.0]],
        vec![[0, 1, 2]],
    )
    .unwrap();
    d.directors = NodalDirectors::new(&d.mesh, vec![[0.0, 0.0, 1.0]; 3], 0.9).unwrap();
    d.constraints = (0..3)
        .map(|node| ShellConstraint {
            node,
            mask: [node == 0, true, true, true, true, true],
            frame: ConstraintFrame::Global,
        })
        .collect();
    d.section_angles_rad = vec![0.0];
    d.budget = MitcBudget {
        max_nodes: 3,
        max_triangles: 1,
    };
    d
}

fn axial_load(force: f64) -> [f64; 18] {
    let mut out = [0.0; 18];
    out[6] = 0.5 * force;
    out[12] = 0.5 * force;
    out
}

fn stiffness(t: f64) -> f64 {
    AREA * H * (E + E_SLOPE * (t - T0))
}
fn expansion(t: f64) -> f64 {
    let x = t - T0;
    ALPHA * x + 0.5 * ALPHA_SLOPE * x * x
}

/// Independently differentiate Psi=.5*k*((q-a)^2+a^2) at fixed total q.
/// The second term retains the constrained transverse normal strain at nu=0.
fn potential(q: f64, t: f64) -> (f64, f64, f64) {
    let k = stiffness(t);
    let a = expansion(t);
    let alpha = ALPHA + ALPHA_SLOPE * (t - T0);
    let squared = (q - a).powi(2) + a * a;
    let psi = 0.5 * k * squared;
    let entropy = -0.5 * AREA * H * E_SLOPE * squared + k * alpha * (q - 2.0 * a);
    (psi, entropy, psi + t * entropy)
}

fn scalar_root(mut lo: f64, mut hi: f64, residual: impl Fn(f64) -> f64) -> f64 {
    assert!(residual(lo) < 0.0 && residual(hi) > 0.0);
    for _ in 0..96 {
        let mid = lo + 0.5 * (hi - lo);
        if residual(mid) > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    lo + 0.5 * (hi - lo)
}

#[test]
fn g1_heated_softening_oscillator_matches_independent_momentum_and_caloric_equations() {
    let cards = [card_with_coefficients(
        "inertial softening and affine expansion",
        E_SLOPE,
        ALPHA,
        ALPHA_SLOPE,
        CP,
    )];
    let body = prepare_body(&cards, &[1], 1.0);
    let laws = body.thermal_laws(false).unwrap();
    let bound = body.bind(&laws, ids()).unwrap();
    let dt = 1e-3;
    let heat_per_step = 1500.0;
    let mut thermal = initial_with_sources(&bound, &[T0], &[heat_per_step / dt]);
    let shell = LayeredThermalShell::try_new(
        &bound,
        mechanical(&cards),
        ReferenceSurface::Midsurface,
        axial_triangle(),
        limits(),
        || false,
    )
    .unwrap();
    let mut previous_force = 1200.0;
    let equilibrium = ReciprocalLayeredThermalShell::try_new(
        &shell,
        &thermal,
        &axial_load(previous_force),
        &[[280.0, 360.0]],
        0.005,
        controls(),
        || false,
    )
    .unwrap();
    let mut velocity = [0.0; 18];
    velocity[6] = 0.025;
    velocity[12] = 0.025;
    let mut coupled = equilibrium
        .begin_transient(&velocity, &[[0.0; 2]], 1e-7, || false)
        .unwrap();
    let mass = RHO * AREA * H / 2.0;
    let capacity = RHO * AREA * H * CP;
    let (mut temperature, mut q, mut v) = (T0, previous_force / stiffness(T0), velocity[6]);
    let initial_energy = thermal.energy_j() + field_internal(coupled.state());
    let mut external_energy = 0.0;
    let mut workspace = LayeredThermalWorkspace::new(1).unwrap();
    close(
        coupled.state().kinetic_energy_j(),
        0.5 * mass * v * v,
        1e-10,
    );
    close(coupled.state().response().displacement_6n[6], q, 1e-11);

    for (index, force) in [1200.0, 1000.0, 1600.0, 1100.0].into_iter().enumerate() {
        let k0 = stiffness(temperature);
        let a0 = expansion(temperature);
        let u0 = potential(q, temperature).2;
        // Eliminate q1,v1 from endpoint-average momentum and midpoint travel.
        // Only the independent scalar finite-heat equation remains to solve.
        let next_q = |t: f64| {
            ((4.0 * mass / (dt * dt) - k0) * q
                + 4.0 * mass * v / dt
                + previous_force
                + force
                + k0 * a0
                + stiffness(t) * expansion(t))
                / (4.0 * mass / (dt * dt) + stiffness(t))
        };
        let stress_work = |t: f64| {
            0.5 * (k0 * (q - a0) + stiffness(t) * (next_q(t) - expansion(t))) * (next_q(t) - q)
        };
        let t1 = scalar_root(temperature, temperature + 3.0, |t| {
            capacity * (t - temperature) - heat_per_step - stress_work(t)
                + potential(next_q(t), t).2
                - u0
        });
        let q1 = next_q(t1);
        let v1 = 2.0 * (q1 - q) / dt - v;
        let k1 = stiffness(t1);
        let (psi1, entropy1, internal1) = potential(q1, t1);
        let mechanical_work = 0.5 * (previous_force + force) * (q1 - q);
        let parameter_work = 0.5 * q * (k1 - k0) * q1;
        let kinetic_change = 0.5 * mass * (v1 * v1 - v * v);
        let step = coupled
            .advance_to(
                &thermal,
                (index + 1) as f64 * dt,
                &axial_load(force),
                conduction(),
                LayeredThermalSolveControls::default(),
                controls(),
                &mut workspace,
                || false,
            )
            .unwrap();
        close(coupled.state().temperatures_k()[0], t1, 3e-8);
        for dof in [6, 12] {
            close(coupled.state().response().displacement_6n[dof], q1, 2e-11);
            close(coupled.state().nodal_velocity_6n()[dof], v1, 5e-8);
        }
        close(
            coupled.state().cell_potentials()[0].free_energy_j,
            psi1,
            1e-8,
        );
        close(
            coupled.state().cell_potentials()[0].entropy_j_k,
            entropy1,
            1e-9,
        );
        close(
            coupled.state().cell_potentials()[0].internal_energy_j,
            internal1,
            1e-7,
        );
        close(
            coupled.state().kinetic_energy_j(),
            0.5 * mass * v1 * v1,
            1e-8,
        );
        close(step.report.kinetic_energy_change_j, kinetic_change, 1e-8);
        close(
            step.report.native.stiffness_parameter_work_j,
            parameter_work,
            1e-9,
        );
        close(
            step.report.native.quadratic_energy_change_j,
            0.5 * (k1 * q1 * q1 - k0 * q * q),
            1e-9,
        );
        close(step.report.native.energy_residual_j, 0.0, 1e-7);
        close(
            step.report.thermoelastic.external_mechanical_work_j,
            mechanical_work,
            1e-8,
        );
        close(
            step.report.thermoelastic.external_heat_j,
            heat_per_step,
            1e-7,
        );
        close(
            step.report.thermoelastic.reversible_heat_j,
            stress_work(t1) - internal1 + u0,
            1e-7,
        );
        close(
            step.report.thermoelastic.combined_energy_residual_j,
            0.0,
            1e-7,
        );
        close(
            step.report.thermoelastic.drilling_gauge_parameter_work_j,
            0.0,
            1e-10,
        );
        assert!(
            parameter_work.abs() > 1e-5,
            "the oracle must exercise actual stiffness work"
        );
        close(step.thermal.mass_kg(), RHO * AREA * H, 1e-12);
        external_energy += heat_per_step + mechanical_work;
        close(
            step.thermal.energy_j() + field_internal(coupled.state()) - initial_energy,
            external_energy,
            4e-7,
        );
        thermal = step.thermal;
        (temperature, q, v, previous_force) = (t1, q1, v1, force);
    }
    assert!(temperature > T0 + 4.0);
    assert!(stiffness(temperature) < 0.995 * stiffness(T0));
}

fn observers(
    shell: &LayeredThermalShell<'_, '_>,
    material: &ResolvedLaminateThermal,
) -> Vec<PreparedElementStrain> {
    let d = shell.discretization();
    d.mesh
        .tris
        .iter()
        .enumerate()
        .map(|(e, _)| {
            let field =
                LaminateThermalField::try_new(material.prepared().clone(), d.section_angles_rad[e])
                    .unwrap();
            PreparedElementStrain::new(&d.mesh, e, &d.directors, field.section(), &mut || false)
                .unwrap()
        })
        .collect()
}

fn components(p: ShellStrainPoint) -> [f64; 8] {
    [
        p.membrane[0],
        p.membrane[1],
        p.membrane[2],
        p.gradient_per_m[0],
        p.gradient_per_m[1],
        p.gradient_per_m[2],
        p.transverse_shear[0],
        p.transverse_shear[1],
    ]
}

/// Independently add the flat triangle's cubic director bubble to the frozen
/// initial nodal observation. For b=27*lambda0*lambda1*lambda2, the curvature
/// follows grad(b) and the MITC assumed shear is [r_y/2,-r_x/2]. This detects
/// silently reinterpreting the retained r in a temperature-dependent T map.
fn check_reference_motion(
    shell: &LayeredThermalShell<'_, '_>,
    reference: &[PreparedElementStrain],
    state: &DynamicLaminateState,
) {
    let mesh = &shell.discretization().mesh;
    for (e, tri) in mesh.tris.iter().enumerate() {
        let nodal =
            core::array::from_fn(|i| state.response().displacement_6n[6 * tri[i / 6] + i % 6]);
        let frozen = reference[e].real_field(&nodal, &mut || false).unwrap();
        let [rx, ry] = state.motion().bubble_coordinates()[e];
        let geometry = mesh.facet(e).unwrap();
        for lambda in [[0.2, 0.3, 0.5], [0.6, 0.3, 0.1], [1.0 / 3.0; 3]] {
            let mut expected = frozen.at_barycentric(lambda).unwrap();
            let derivative: [f64; 2] = core::array::from_fn(|axis| {
                27.0 * (lambda[1] * lambda[2] * geometry.gradient[0][axis]
                    + lambda[0] * lambda[2] * geometry.gradient[1][axis]
                    + lambda[0] * lambda[1] * geometry.gradient[2][axis])
            });
            expected.gradient_per_m[0] += ry * derivative[0];
            expected.gradient_per_m[1] -= rx * derivative[1];
            expected.gradient_per_m[2] += ry * derivative[1] - rx * derivative[0];
            expected.transverse_shear[0] += 0.5 * ry;
            expected.transverse_shear[1] -= 0.5 * rx;
            let actual = state.strain_fields()[e].at_barycentric(lambda).unwrap();
            for (actual, expected) in components(actual).into_iter().zip(components(expected)) {
                close(actual, expected, 2e-10);
            }
        }
    }
}

#[test]
fn g1_unequal_heating_retains_reference_director_motion_and_atomic_energy_history() {
    let cards = [
        card_with_coefficients("moving lower modulus", -0.002 * E, ALPHA, 0.0, CP),
        card_with_coefficients("moving upper modulus", -0.005 * E, 2.0 * ALPHA, 0.0, CP),
    ];
    let body = prepare_body(&cards, &[1, 1], 1.0);
    let laws = body.thermal_laws(false).unwrap();
    let bound = body.bind(&laws, ids()).unwrap();
    let heater = [3e6, 1e6];
    let thermal = initial_with_sources(&bound, &[310.0, 340.0], &heater);
    let original_cells = thermal.cells().clone();
    let shell = LayeredThermalShell::try_new(
        &bound,
        mechanical(&cards),
        ReferenceSurface::Midsurface,
        discretization(),
        limits(),
        || false,
    )
    .unwrap();
    let seed = || {
        ReciprocalLayeredThermalShell::try_new(
            &shell,
            &thermal,
            &[0.0; 24],
            &[[300.0, 350.0]; 2],
            0.005,
            controls(),
            || false,
        )
        .unwrap()
        .begin_transient(&[0.0; 24], &[[2.0, -1.0], [-0.5, 1.5]], 1e-7, || false)
        .unwrap()
    };
    let mut coupled = seed();
    let reference = observers(&shell, coupled.state().material());
    let initial_energy = thermal.energy_j() + field_internal(coupled.state());
    let initial_q = coupled.state().motion().coordinates().to_vec();
    let initial_v = coupled.state().motion().velocity().to_vec();
    let initial_temperatures = coupled.state().temperatures_k().to_vec();
    // Pure residual-director rates have no nodal/cross-inertia contribution.
    // Integral(b^2)/A=81/280; source rho(T) must not change the reference mass.
    let i2 = RHO * H.powi(3) / 12.0;
    let kinetic = 0.5
        * (AREA / 2.0)
        * i2
        * (81.0 / 280.0)
        * (2.0_f64.powi(2) + 1.0 + 0.5_f64.powi(2) + 1.5_f64.powi(2));
    close(coupled.state().kinetic_energy_j(), kinetic, 1e-11);
    check_reference_motion(&shell, &reference, coupled.state());
    let target = 1e-4;
    let mut workspace = LayeredThermalWorkspace::new(2).unwrap();
    assert!(matches!(
        coupled.advance_to(
            &thermal,
            target,
            &[0.0; 24],
            conduction(),
            LayeredThermalSolveControls::default(),
            LaminateThermoelasticControls {
                maximum_iterations: 1,
                ..controls()
            },
            &mut workspace,
            || false,
        ),
        Err(LaminateThermoelasticError::NotConverged { .. })
    ));
    assert_eq!(coupled.state().motion().coordinates(), initial_q);
    assert_eq!(coupled.state().motion().velocity(), initial_v);
    assert_eq!(coupled.state().temperatures_k(), initial_temperatures);
    close(coupled.time_s(), 0.0, 0.0);

    let step = coupled
        .advance_to(
            &thermal,
            target,
            &[0.0; 24],
            conduction(),
            LayeredThermalSolveControls::default(),
            controls(),
            &mut workspace,
            || false,
        )
        .unwrap();
    check_reference_motion(&shell, &reference, coupled.state());
    assert!(coupled.state().temperatures_k()[0] > initial_temperatures[0] + 0.1);
    assert!(
        coupled
            .state()
            .transverse_shear_entropy_j_k()
            .iter()
            .any(|s| s.abs() > 1e-10)
    );
    assert!(step.report.bubble_material_work_j.abs() > 1e-9);
    assert!(
        step.report
            .thermoelastic
            .transverse_shear_reversible_heat_j
            .abs()
            > 1e-8
    );
    close(
        step.report.thermoelastic.external_heat_j,
        heater.iter().sum::<f64>() * target,
        1e-7,
    );
    close(
        step.report.thermoelastic.external_mechanical_work_j,
        0.0,
        1e-12,
    );
    close(
        step.report.thermoelastic.drilling_gauge_parameter_work_j,
        0.0,
        1e-7,
    );
    close(step.report.native.energy_residual_j, 0.0, 1e-7);
    close(
        step.report.thermoelastic.combined_energy_residual_j,
        0.0,
        1e-7,
    );
    close(
        step.thermal.energy_j() + field_internal(coupled.state()) - initial_energy,
        heater.iter().sum::<f64>() * target,
        2e-7,
    );
    close(step.thermal.mass_kg(), RHO * AREA * H, 1e-12);
    assert_eq!(thermal.cells(), &original_cells);

    // Confirm this fixture really changes elastic condensation: using the
    // current nodal-only map on an arbitrary fixed probe gives different
    // curvature/shear. Accepted physical motion above must still use T0.
    let current = observers(&shell, coupled.state().material());
    let mut probe = [0.0; 18];
    probe[6] = 3e-4;
    probe[8] = 2e-3;
    probe[15] = 1e-3;
    let mut changed = 0.0_f64;
    for (old, new) in reference.iter().zip(&current) {
        let a = old.real_field(&probe, &mut || false).unwrap();
        let b = new.real_field(&probe, &mut || false).unwrap();
        for (a, b) in a
            .bernstein_controls()
            .into_iter()
            .zip(b.bernstein_controls())
        {
            for (a, b) in components(a).into_iter().zip(components(b)) {
                changed = changed.max((b - a).abs());
            }
        }
    }
    assert!(
        changed > 1e-12,
        "the heterogeneous fixture must change the condensation map"
    );

    let accepted_q = coupled.state().motion().coordinates().to_vec();
    let accepted_v = coupled.state().motion().velocity().to_vec();
    assert!(matches!(
        coupled.advance_to(
            &thermal,
            2.0 * target,
            &[0.0; 24],
            conduction(),
            LayeredThermalSolveControls::default(),
            controls(),
            &mut workspace,
            || false,
        ),
        Err(LaminateThermoelasticError::Invalid(_))
    ));
    assert_eq!(coupled.state().motion().coordinates(), accepted_q);
    assert_eq!(coupled.state().motion().velocity(), accepted_v);
    close(coupled.time_s(), target, 0.0);

    let mut measured = seed();
    let mut measured_workspace = LayeredThermalWorkspace::new(2).unwrap();
    let mut checkpoints = 0;
    measured
        .advance_to(
            &thermal,
            target,
            &[0.0; 24],
            conduction(),
            LayeredThermalSolveControls::default(),
            controls(),
            &mut measured_workspace,
            || {
                checkpoints += 1;
                false
            },
        )
        .unwrap();
    let mut cancelled = seed();
    let mut cancelled_workspace = LayeredThermalWorkspace::new(2).unwrap();
    let mut checkpoint = 0;
    assert!(matches!(
        cancelled.advance_to(
            &thermal,
            target,
            &[0.0; 24],
            conduction(),
            LayeredThermalSolveControls::default(),
            controls(),
            &mut cancelled_workspace,
            || {
                checkpoint += 1;
                checkpoint == checkpoints
            },
        ),
        Err(LaminateThermoelasticError::Cancelled)
    ));
    assert_eq!(checkpoint, checkpoints);
    assert_eq!(cancelled.state().motion().coordinates(), initial_q);
    assert_eq!(cancelled.state().motion().velocity(), initial_v);
    assert_eq!(cancelled.state().temperatures_k(), initial_temperatures);
    close(cancelled.state().kinetic_energy_j(), kinetic, 1e-11);
    close(cancelled.time_s(), 0.0, 0.0);
    assert_eq!(thermal.cells(), &original_cells);
}
