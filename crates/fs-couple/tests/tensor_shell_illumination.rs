//! G0/G1/G3: incident sound drives the original six-coordinate material shell.
//! Synthetic gas and structural laws isolate coupling and coherent field work;
//! these checks do not certify continuum or measured material accuracy.

use fs_bem::helmholtz::Medium;
use fs_bem::helmholtz::bc::{IncidentField, Monopole, PlaneWave};
use fs_bem::helmholtz::molecular::MolecularExterior;
use fs_bem::near_field::Options;
use fs_bem::reciprocal_radiation::open_screen::{
    MolecularScreenIlluminationLimits, MolecularScreenPortLimits, OpenOptions, OpenPanels,
    ScreenComplexValueLimits, ScreenPortLimits,
};
use fs_couple::tensor_shell_radiation::open::molecular::OpenTensorShellIlluminationControls;
use fs_couple::tensor_shell_radiation::open::{
    OpenTensorShellRadiation, OpenTensorShellRadiationControls,
    OpenTensorShellRadiationFrequencyControls,
};
use fs_exec::{CancelGate, Cx};
use fs_material::visco::tensor::{Mandel6, PronyBranch, TensorProny};
use fs_material::{
    gas::{GasSpec, GasState},
    relaxation::{RelaxingMedium, input::read_relaxation},
};
use fs_math::c64::C64;
use fs_modal::harmonic::{HarmonicPortControls, HarmonicResponseControls};
use fs_plate::shell::viscosity::ShellViscosityBudget;
use fs_plate::shell::viscosity::harmonic::{TensorPronyShell, TensorPronyShellSection};
use fs_plate::shell::{CstDktMass, ShellMesh, ShellSupport};

const OMEGA: f64 = 650.0;

fn with_cx(run: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        run(
            &Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 47,
                    kernel_id: 59,
                    tile: 0,
                    iteration: 0,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            ),
            &gate,
        );
    });
}

fn source(height: f64, cx: &Cx<'_>) -> TensorPronyShell {
    // Explicit numerical polymer-like material, not a sourced material claim.
    let c: Mandel6 = core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            1e6 * if i < 3 && j < 3 {
                4.0 + if i == j { 6.0 } else { 0.0 }
            } else if i == j {
                6.0
            } else {
                0.0
            }
        })
    });
    let scale = |s: f64| c.map(|r| r.map(|x| x * s));
    let law = TensorProny::new_mandel(
        c,
        vec![
            PronyBranch::new(scale(0.4), 0.001).unwrap(),
            PronyBranch::new(scale(0.2), 0.02).unwrap(),
        ],
        0.02,
    )
    .unwrap()
    .with_kelvin_viscosity_mandel(scale(1e-5))
    .unwrap()
    .with_band((0.0, 1500.0))
    .unwrap();
    let section = TensorPronyShellSection::new(law, 0.01, 1200.0).unwrap();
    let mesh = ShellMesh::new(
        vec![
            [-0.2, -0.15, 0.0],
            [0.2, -0.15, 0.0],
            [0.2, 0.15, 0.0],
            [-0.2, 0.15, 0.0],
            [0.0, 0.0, height],
        ],
        vec![[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]],
    )
    .unwrap();
    TensorPronyShell::assemble(
        mesh,
        vec![section; 4],
        &[0, 1, 2, 3],
        ShellSupport::Clamped,
        CstDktMass::Consistent,
        0.0,
        ShellViscosityBudget {
            max_elements: 4,
            max_full_dofs: 30,
            max_dofs: 6,
            max_modes: 6,
            max_work: 20_000_000,
        },
        cx,
    )
    .unwrap()
}

fn medium(exterior: &MolecularExterior) -> Medium {
    let gas = exterior.source().gas();
    Medium {
        density: gas.density,
        sound_speed: gas.sound_speed,
    }
}
fn geometry_controls(level: u32) -> OpenTensorShellRadiationControls {
    OpenTensorShellRadiationControls {
        acoustic_refinements: level,
        maximum_panels: 64,
        maximum_coordinates: 6,
        maximum_setup_bytes: 1 << 20,
        maximum_setup_work: 10_000_000,
    }
}
fn boundary_options() -> OpenOptions {
    OpenOptions {
        edge_level: 0,
        maximum_edge_level: 0,
        triangle_order: 8,
        angular_refinement: 0,
        maximum_unknowns: 64,
        maximum_cache_bytes: 16 << 20,
        maximum_work: 2_000_000_000,
        workers: 1,
    }
}
fn frequency_controls() -> OpenTensorShellRadiationFrequencyControls {
    OpenTensorShellRadiationFrequencyControls {
        mechanical: HarmonicResponseControls {
            maximum_coordinates: 6,
            maximum_payload_bytes: 16 << 20,
            maximum_relative_residual: 2e-11,
        },
        ports: HarmonicPortControls {
            maximum_ports: 6,
            maximum_projection_terms: 1_000_000,
            maximum_payload_bytes: 4 << 20,
        },
        boundary: ScreenPortLimits {
            maximum_ports: 6,
            maximum_matrix_entries: 36,
            maximum_jump_entries: 6 * 64,
            maximum_work: 2_000_000_000,
        },
        maximum_boundary_payload_bytes: 16 << 20,
        maximum_radiated_power_work: 512_000_000,
        power_relative_tolerance: 1e-8,
    }
}
fn force() -> [C64; 6] {
    [
        C64::new(1e-4, -2e-5),
        C64::new(-3e-5, 1e-5),
        C64::new(2e-5, 4e-5),
        C64::new(1e-6, -2e-7),
        C64::new(3e-7, 1e-7),
        C64::new(-5e-7, 1e-7),
    ]
}
fn near(a: C64, b: C64) {
    assert!((a - b).abs() <= 3e-8 * b.abs().max(1e-23), "{a:?} != {b:?}");
}
fn near_power(a: f64, b: f64) {
    assert!((a - b).abs() <= 3e-8 * b.abs().max(1e-26), "{a} != {b}");
}

// Independent closed-form center shape on each of the four original facets.
// It remains this shape on the shallow tent: height changes facet normals.
fn center_shape(point: [f64; 3]) -> f64 {
    1.0 - (point[0].abs() / 0.2).max(point[1].abs() / 0.15)
}
fn area_normal(points: [[f64; 3]; 3]) -> (f64, [f64; 3]) {
    let a: [f64; 3] = core::array::from_fn(|d| points[1][d] - points[0][d]);
    let b: [f64; 3] = core::array::from_fn(|d| points[2][d] - points[0][d]);
    let n = [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ];
    let length = n.iter().map(|x| x * x).sum::<f64>().sqrt();
    (0.5 * length, n.map(|x| x / length))
}
// Degree-two quadrature integrates P1 pressure times the ORIGINAL P1 velocity
// exactly. This does not read adapter rows or the production triangle mass map.
fn projected_effort(surface: &OpenPanels, jump: &[C64]) -> [C64; 6] {
    let mut effort = [C64::ZERO; 6];
    for triangle in surface.triangles() {
        let points = triangle.map(|v| surface.vertices()[v]);
        let (area, normal) = area_normal(points);
        for q in [
            [2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0],
            [1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0],
            [1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0],
        ] {
            let position = core::array::from_fn(|d| (0..3).map(|a| q[a] * points[a][d]).sum());
            let pressure = (0..3).fold(C64::ZERO, |sum, a| sum + jump[triangle[a]].scale(q[a]));
            for d in 0..3 {
                effort[d] =
                    effort[d] + pressure.scale(area / 3.0 * center_shape(position) * normal[d]);
            }
        }
    }
    effort
}

fn exterior(relaxing: bool, pressure: f64) -> MolecularExterior {
    let data = if relaxing {
        "declared 1 synthetic 1 vibration 1 900 1 power 101.325 293.15 0 200 500"
    } else {
        "declared 1 synthetic 0"
    };
    let law = read_relaxation(&mut data.split_ascii_whitespace()).unwrap();
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, pressure).unwrap();
    let state = RelaxingMedium::new(gas, vec![(law.species()[0].species.clone(), 1.0)]).unwrap();
    MolecularExterior::new(law.evaluate(&state).unwrap(), [100.0, 2000.0]).unwrap()
}
fn controls() -> OpenTensorShellIlluminationControls {
    let native = frequency_controls();
    OpenTensorShellIlluminationControls {
        mechanical: native.mechanical,
        ports: native.ports,
        boundary: MolecularScreenIlluminationLimits {
            boundary: MolecularScreenPortLimits {
                ports: native.boundary,
                operator: ScreenComplexValueLimits {
                    maximum_matrix_entries: 4096,
                    maximum_workspace_bytes: 4 << 20,
                    maximum_work: 2_000_000_000,
                    maximum_growth_exponent: 8.0,
                },
                maximum_workspace_bytes: 4 << 20,
            },
            maximum_incident_cases: 4,
            maximum_incident_jump_entries: 256,
            maximum_source_point_products: 1_000_000,
            minimum_source_clearance_m: 0.001,
            maximum_live_bytes: 16 << 20,
        },
        maximum_response_bytes: 4 << 20,
        maximum_response_work: 10_000_000,
        power_relative_tolerance: 1e-8,
    }
}
fn plane(amplitude: C64) -> IncidentField<'static> {
    IncidentField::PlaneWave(PlaneWave {
        direction: [0.36, 0.48, -0.8],
        amplitude_pa: amplitude,
    })
}
fn work(force: &[C64], velocity: &[C64]) -> f64 {
    0.5 * force
        .iter()
        .zip(velocity)
        .map(|(&f, &v)| (f * v.conj()).re)
        .sum::<f64>()
}
fn close(a: C64, b: C64, relative: f64, floor: f64) {
    assert!(
        (a - b).abs() <= relative * a.abs().max(b.abs()) + floor,
        "{a:?} != {b:?}"
    );
}
fn identity() -> Vec<f64> {
    (0..36)
        .map(|i| if i / 6 == i % 6 { 1.0 } else { 0.0 })
        .collect()
}

#[test]
fn zero_illumination_recovers_full_coordinate_lossless_radiation() {
    with_cx(|cx, _| {
        let gas = exterior(false, 101325.0);
        let shell = source(0.0, cx);
        let geometry =
            OpenTensorShellRadiation::new(&shell, medium(&gas), geometry_controls(1), cx).unwrap();
        let prepared = geometry.prepare(1200.0, boundary_options(), cx).unwrap();
        let native = prepared
            .prepare_frequency(OMEGA, frequency_controls(), cx)
            .unwrap();
        let loaded = prepared
            .prepare_molecular_illumination(&gas, OMEGA, &[plane(C64::ZERO)], controls(), cx)
            .unwrap();
        assert!(core::ptr::eq(loaded.shell_preparation().source(), &shell));
        assert_eq!(loaded.radiation().motion().ports(), 6);
        assert_eq!(loaded.radiation().motion().pressure_jumps_pa().len(), 6);
        let actual = loaded.solve(0, &force(), cx).unwrap();
        let expected = native.solve(&force(), cx).unwrap();
        assert_eq!(actual.motion().displacement.len(), 6);
        assert!(
            actual.motion().displacement[3..]
                .iter()
                .all(|z| z.abs() > 0.0)
        );
        assert!(
            actual
                .blocked_acoustic_force()
                .iter()
                .all(|&z| z == C64::ZERO)
        );
        for (a, b) in actual
            .motion()
            .displacement
            .iter()
            .zip(&expected.motion().displacement)
        {
            close(*a, *b, 2e-10, 1e-23);
        }
        for (&a, &b) in loaded.impedance().iter().zip(native.impedance()) {
            close(a, b, 2e-10, 1e-24);
        }
        near_power(
            actual.power().intrinsic_material_heat_w,
            expected.motion().average_viscous_power_w,
        );
        near_power(
            actual.power().total_exterior_work_w,
            expected.radiated_power_w(),
        );
        near_power(
            actual.power().motion_exterior_work_w,
            expected.radiated_power_w(),
        );
        let jump = actual.total_pressure_jump(cx).unwrap();
        for (a, b) in jump.iter().zip(expected.pressure_jump(cx).unwrap()) {
            close(*a, b, 2e-10, 1e-22);
        }
        let points = [[0.04, -0.03, 0.36], [0.04, -0.03, -0.36]];
        let receivers = loaded
            .prepare_receivers(&points, 0.001, Options::default(), 2 << 20, cx)
            .unwrap();
        let observed = receivers.evaluate(&actual, cx).unwrap();
        let expected_pressure = expected
            .pressure_at(&points, 0.001, Options::default(), cx)
            .unwrap();
        for (&a, &b) in observed.scalar.pressure.iter().zip(&expected_pressure) {
            close(a, b, 2e-6, 1e-20);
        }
        // The entire zero-offset rotational source is retained without adding
        // fictitious rotational normal velocity or any second wet-face debit.
        for d in [0, 1, 3, 4, 5] {
            assert!(
                loaded.radiation().motion().pressure_jumps_pa()[d]
                    .iter()
                    .all(|&z| z == C64::ZERO)
            );
        }
    });
}

#[test]
fn sound_drives_original_relaxing_shell_with_independent_signed_traction_work() {
    with_cx(|cx, _| {
        let gas = exterior(true, 101325.0);
        let shell = source(0.02, cx);
        let geometry =
            OpenTensorShellRadiation::new(&shell, medium(&gas), geometry_controls(1), cx).unwrap();
        let prepared = geometry.prepare(1200.0, boundary_options(), cx).unwrap();
        let loaded = prepared
            .prepare_molecular_illumination(
                &gas,
                OMEGA,
                &[plane(C64::new(2.0, 0.4))],
                controls(),
                cx,
            )
            .unwrap();
        assert!(loaded.radiation().motion().wavenumber_per_m().im > 0.0);
        assert!(core::ptr::eq(loaded.shell_preparation().source(), &shell));
        let zero = [C64::ZERO; 6];
        for external in [&zero, &force()] {
            let response = loaded.solve(0, external, cx).unwrap();
            let velocity = &response.motion().velocity;
            assert!(velocity.iter().any(|z| z.abs() > 0.0));
            let blocked = projected_effort(
                geometry.surface(),
                &loaded.radiation().blocked_pressure_jumps_pa()[0],
            )
            .map(|z| -z);
            for (&a, &b) in response.blocked_acoustic_force().iter().zip(&blocked) {
                close(a, b, 3e-11, 1e-22);
            }
            assert_eq!(&blocked[3..], &[C64::ZERO; 3]);
            let total_jump = response.total_pressure_jump(cx).unwrap();
            let reaction = projected_effort(geometry.surface(), &total_jump).map(|z| -z);
            for (&a, &b) in response.total_fluid_reaction().iter().zip(&reaction) {
                close(a, b, 3e-11, 1e-22);
            }
            // Independently apply the exact pressure-jump force to the native
            // TensorPronyShell's complete physical displacement/rotation pencil.
            let effective: Vec<_> = external
                .iter()
                .zip(&blocked)
                .map(|(&a, &b)| a + b)
                .collect();
            let original = shell
                .prepare_harmonic_response_with_port_impedance(
                    OMEGA,
                    &identity(),
                    6,
                    loaded.impedance(),
                    controls().mechanical,
                    controls().ports,
                    cx,
                )
                .unwrap()
                .solve(&effective, cx)
                .unwrap();
            for (&a, &b) in response
                .motion()
                .displacement
                .iter()
                .zip(&original.displacement)
            {
                close(a, b, 4e-11, 1e-23);
            }
            let p = response.power();
            assert!(p.intrinsic_material_heat_w > 0.0);
            near_power(
                p.intrinsic_material_heat_w,
                original.average_viscous_power_w,
            );
            near_power(p.external_mechanical_input_w, work(external, velocity));
            near_power(p.blocked_acoustic_input_w, work(&blocked, velocity));
            near_power(p.total_exterior_work_w, -work(&reaction, velocity));
            near_power(
                p.motion_exterior_work_w,
                loaded
                    .radiation()
                    .motion()
                    .boundary_work_w(velocity, || false)
                    .unwrap(),
            );
            let scale = p.apparent_power_scale_w;
            assert!(
                (p.total_exterior_work_w - p.motion_exterior_work_w + p.blocked_acoustic_input_w)
                    .abs()
                    < 1e-8 * scale
            );
            assert!(
                (p.external_mechanical_input_w
                    - p.intrinsic_material_heat_w
                    - p.total_exterior_work_w)
                    .abs()
                    < 1e-8 * scale
            );
            assert!(
                p.relative_balance < 1e-8
                    && p.relative_motion_work_defect < 1e-8
                    && p.relative_traction_defect < 1e-8
            );
            if external == &zero {
                assert_eq!(p.external_mechanical_input_w, 0.0);
                assert!(p.total_exterior_work_w < 0.0);
                near_power(-p.total_exterior_work_w, p.intrinsic_material_heat_w);
            }
        }
        // A compensating mechanical force keeps the original structure still.
        // Its nonzero blocked field remains: zero motion is not zero scattering.
        let cancellation: Vec<_> = loaded.radiation().blocked_forces()[0]
            .iter()
            .map(|&z| -z)
            .collect();
        let held = loaded.solve(0, &cancellation, cx).unwrap();
        assert!(held.motion().velocity.iter().all(|&z| z == C64::ZERO));
        assert_eq!(
            held.total_pressure_jump(cx).unwrap(),
            loaded.radiation().blocked_pressure_jumps_pa()[0]
        );
        assert!(
            held.total_pressure_jump(cx)
                .unwrap()
                .iter()
                .any(|z| z.abs() > 0.0)
        );
        assert_eq!(held.power().intrinsic_material_heat_w, 0.0);
        assert_eq!(held.power().total_exterior_work_w, 0.0);
    });
}

// Independently generated Gauss rule and triangle parameterization. The fixed
// P1 jump is observed with the analytic double-layer kernel and its derivative;
// no production receiver rows or incident samplers enter this reference.
fn gauss(order: usize) -> Vec<(f64, f64)> {
    (0..order)
        .map(|i| {
            let mut x = (std::f64::consts::PI * (i as f64 + 0.75) / (order as f64 + 0.5)).cos();
            let mut derivative = 0.0;
            for _ in 0..32 {
                let (mut a, mut b) = (1.0, x);
                for n in 2..=order {
                    let c = ((2 * n - 1) as f64 * x * b - (n - 1) as f64 * a) / n as f64;
                    a = b;
                    b = c;
                }
                derivative = order as f64 * (x * b - a) / (x * x - 1.0);
                let update = b / derivative;
                x -= update;
                if update.abs() < 2e-16 {
                    break;
                }
            }
            (
                0.5 * (1.0 + x),
                1.0 / ((1.0 - x * x) * derivative * derivative),
            )
        })
        .collect()
}
fn reference_field(
    surface: &OpenPanels,
    jump: &[C64],
    k: C64,
    rho: f64,
    incident: PlaneWave,
    point: [f64; 3],
) -> (C64, [C64; 3]) {
    let coordinate: f64 = (0..3).map(|i| incident.direction[i] * point[i]).sum();
    let argument = k.re * coordinate;
    let direct = incident.amplitude_pa
        * C64::new(argument.cos(), argument.sin()).scale((-k.im * coordinate).exp());
    let mut pressure = direct;
    let mut velocity = incident
        .direction
        .map(|d| (k * direct).scale(d / (rho * OMEGA)));
    let rule = gauss(16);
    let ik = C64::new(-k.im, k.re);
    for ids in surface.triangles() {
        let triangle = ids.map(|id| surface.vertices()[id]);
        let (area, normal) = area_normal(triangle);
        for &(a, wa) in &rule {
            for &(b, wb) in &rule {
                let shape = [1.0 - a, a * (1.0 - b), a * b];
                let position: [f64; 3] =
                    core::array::from_fn(|d| (0..3).map(|i| shape[i] * triangle[i][d]).sum());
                let r_vec: [f64; 3] = core::array::from_fn(|d| point[d] - position[d]);
                let r = r_vec.iter().map(|x| x * x).sum::<f64>().sqrt();
                let projection: f64 = (0..3).map(|d| normal[d] * r_vec[d]).sum();
                let phase = k.re * r;
                let g = C64::new(phase.cos(), phase.sin())
                    .scale((-k.im * r).exp() / (4.0 * std::f64::consts::PI * r));
                let p = (0..3).fold(C64::ZERO, |sum, i| sum + jump[ids[i]].scale(shape[i]));
                let weighted = p.scale(2.0 * area * a * wa * wb);
                pressure =
                    pressure + weighted * g * (C64::new(1.0 / r, 0.0) - ik).scale(projection / r);
                let normal_coefficient = C64::new(1.0 / (r * r), 0.0) - ik.scale(1.0 / r);
                let radial_coefficient = (k * k).scale(1.0 / (r * r)) + ik.scale(3.0 / r.powi(3))
                    - C64::new(3.0 / r.powi(4), 0.0);
                for d in 0..3 {
                    let gradient = g
                        * (normal_coefficient.scale(normal[d])
                            + radial_coefficient.scale(projection * r_vec[d]));
                    velocity[d] =
                        velocity[d] + weighted * gradient * C64::new(0.0, -1.0 / (rho * OMEGA));
                }
            }
        }
    }
    (pressure, velocity)
}

#[test]
fn coherent_incident_and_mechanical_forces_reuse_total_pressure_velocity_rows() {
    with_cx(|cx, _| {
        let gas = exterior(true, 101325.0);
        let shell = source(0.02, cx);
        let geometry =
            OpenTensorShellRadiation::new(&shell, medium(&gas), geometry_controls(1), cx).unwrap();
        let prepared = geometry.prepare(1200.0, boundary_options(), cx).unwrap();
        let amplitude = C64::new(2.0, 0.4);
        let phase = C64::new(-0.4, 0.7);
        let beta = C64::new(0.3, -0.4);
        let incidents = [plane(amplitude), plane(amplitude * phase), plane(C64::ZERO)];
        let loaded = prepared
            .prepare_molecular_illumination(&gas, OMEGA, &incidents, controls(), cx)
            .unwrap();
        let zero = [C64::ZERO; 6];
        let blocked = &loaded.radiation().blocked_forces()[0];
        let force: Vec<_> = blocked.iter().map(|&v| v * beta).collect();
        let sound = loaded.solve(0, &zero, cx).unwrap();
        let phased = loaded.solve(1, &zero, cx).unwrap();
        let forced = loaded.solve(2, &force, cx).unwrap();
        let both = loaded.solve(0, &force, cx).unwrap();
        let cancel: Vec<_> = blocked.iter().map(|&v| -v).collect();
        let held = loaded.solve(0, &cancel, cx).unwrap();
        for i in 0..6 {
            close(
                phased.motion().velocity[i],
                sound.motion().velocity[i] * phase,
                3e-10,
                1e-22,
            );
            close(
                both.motion().velocity[i],
                sound.motion().velocity[i] + forced.motion().velocity[i],
                3e-10,
                1e-22,
            );
        }
        near_power(
            phased.power().intrinsic_material_heat_w,
            phase.norm_sq() * sound.power().intrinsic_material_heat_w,
        );
        near_power(
            both.power().intrinsic_material_heat_w,
            (C64::ONE + beta).norm_sq() * sound.power().intrinsic_material_heat_w,
        );
        let incoherent =
            sound.power().intrinsic_material_heat_w + forced.power().intrinsic_material_heat_w;
        assert!(
            (both.power().intrinsic_material_heat_w - incoherent).abs()
                > 0.3 * sound.power().intrinsic_material_heat_w
        );

        let points = [[0.04, -0.03, 0.36], [0.04, -0.03, -0.36]];
        let receivers = loaded
            .prepare_receivers(&points, 0.001, Options::default(), 2 << 20, cx)
            .unwrap();
        let a = receivers.evaluate(&sound, cx).unwrap();
        let b = receivers.evaluate(&forced, cx).unwrap();
        let combined = receivers.evaluate(&both, cx).unwrap();
        let shifted = receivers.evaluate(&phased, cx).unwrap();
        let stationary = receivers.evaluate(&held, cx).unwrap();
        let repeated = receivers.evaluate(&both, cx).unwrap();
        assert_eq!(combined.scalar.pressure, repeated.scalar.pressure);
        assert_eq!(
            combined.particle_velocity_m_s,
            repeated.particle_velocity_m_s
        );
        for point in 0..points.len() {
            near(
                combined.scalar.pressure[point],
                a.scalar.pressure[point] + b.scalar.pressure[point],
            );
            near(
                shifted.scalar.pressure[point],
                a.scalar.pressure[point] * phase,
            );
            for d in 0..3 {
                close(
                    combined.particle_velocity_m_s[point][d],
                    a.particle_velocity_m_s[point][d] + b.particle_velocity_m_s[point][d],
                    3e-10,
                    1e-20,
                );
                close(
                    shifted.particle_velocity_m_s[point][d],
                    a.particle_velocity_m_s[point][d] * phase,
                    3e-10,
                    1e-20,
                );
            }
        }
        let IncidentField::PlaneWave(wave) = incidents[0] else {
            unreachable!()
        };
        for (response, observed) in [(&both, &combined), (&held, &stationary)] {
            let jump = response.total_pressure_jump(cx).unwrap();
            for (index, &point) in points.iter().enumerate() {
                let (p, v) = reference_field(
                    geometry.surface(),
                    &jump,
                    loaded.radiation().motion().wavenumber_per_m(),
                    medium(&gas).density,
                    wave,
                    point,
                );
                close(observed.scalar.pressure[index], p, 2e-7, 1e-18);
                for (d, value) in v.into_iter().enumerate() {
                    close(observed.particle_velocity_m_s[index][d], value, 2e-7, 1e-20);
                }
            }
        }
    });
}

#[test]
fn original_source_domains_resource_caps_and_receiver_binding_remain_active() {
    with_cx(|cx, gate| {
        let gas = exterior(true, 101325.0);
        let shell = source(0.02, cx);
        let geometry =
            OpenTensorShellRadiation::new(&shell, medium(&gas), geometry_controls(1), cx).unwrap();
        let prepared = geometry.prepare(1200.0, boundary_options(), cx).unwrap();
        let incident = [plane(C64::new(2.0, 0.4))];
        let different_gas = exterior(true, 202650.0);
        assert!(
            prepared
                .prepare_molecular_illumination(&different_gas, OMEGA, &incident, controls(), cx)
                .is_err()
        );
        let narrow_gas = MolecularExterior::new(gas.source().clone(), [100.0, 600.0]).unwrap();
        assert!(
            prepared
                .prepare_molecular_illumination(&narrow_gas, OMEGA, &incident, controls(), cx)
                .is_err()
        );
        for omega in [0.0, f64::NAN, 1201.0] {
            assert!(
                prepared
                    .prepare_molecular_illumination(&gas, omega, &incident, controls(), cx)
                    .is_err()
            );
        }
        for capped in [
            OpenTensorShellIlluminationControls {
                maximum_response_bytes: 1,
                ..controls()
            },
            OpenTensorShellIlluminationControls {
                maximum_response_work: 1,
                ..controls()
            },
            OpenTensorShellIlluminationControls {
                mechanical: HarmonicResponseControls {
                    maximum_coordinates: 5,
                    ..controls().mechanical
                },
                ..controls()
            },
            OpenTensorShellIlluminationControls {
                ports: HarmonicPortControls {
                    maximum_ports: 5,
                    ..controls().ports
                },
                ..controls()
            },
            OpenTensorShellIlluminationControls {
                boundary: MolecularScreenIlluminationLimits {
                    maximum_source_point_products: 1,
                    ..controls().boundary
                },
                ..controls()
            },
            OpenTensorShellIlluminationControls {
                boundary: MolecularScreenIlluminationLimits {
                    maximum_incident_jump_entries: 1,
                    ..controls().boundary
                },
                ..controls()
            },
        ] {
            assert!(
                prepared
                    .prepare_molecular_illumination(&gas, OMEGA, &incident, capped, cx)
                    .is_err()
            );
        }
        assert!(
            prepared
                .prepare_molecular_illumination(&gas, OMEGA, &[], controls(), cx)
                .is_err()
        );
        let emitter = [Monopole {
            position_m: [0.0, 0.0, 0.0205],
            volume_velocity_m3_s: C64::new(1e-4, 0.0),
        }];
        assert!(
            prepared
                .prepare_molecular_illumination(
                    &gas,
                    OMEGA,
                    &[IncidentField::Monopoles(&emitter)],
                    controls(),
                    cx
                )
                .is_err()
        );

        let loaded = prepared
            .prepare_molecular_illumination(&gas, OMEGA, &incident, controls(), cx)
            .unwrap();
        assert!(loaded.solve(1, &force(), cx).is_err());
        assert!(loaded.solve(0, &force()[..5], cx).is_err());
        let mut malformed = force();
        malformed[4] = C64::new(f64::NAN, 0.0);
        assert!(loaded.solve(0, &malformed, cx).is_err());
        // The original tensor strain domain still checks the entire harmonic cycle.
        assert!(loaded.solve(0, &[C64::new(1e10, 0.0); 6], cx).is_err());
        let response = loaded.solve(0, &force(), cx).unwrap();
        assert!(
            loaded
                .prepare_receivers(&[[0.0, 0.0, 0.02]], 0.001, Options::default(), 2 << 20, cx)
                .is_err()
        );
        let points = [[0.04, -0.03, 0.36], [0.04, -0.03, -0.36]];
        assert!(
            loaded
                .prepare_receivers(&points, 0.001, Options::default(), 1, cx)
                .is_err()
        );
        let receivers = loaded
            .prepare_receivers(&points, 0.001, Options::default(), 2 << 20, cx)
            .unwrap();
        let other = prepared
            .prepare_molecular_illumination(&gas, OMEGA, &incident, controls(), cx)
            .unwrap();
        let foreign = other.solve(0, &force(), cx).unwrap();
        assert!(receivers.evaluate(&foreign, cx).is_err());
        assert!(receivers.evaluate(&response, cx).is_ok());
        gate.request();
        assert!(loaded.solve(0, &force(), cx).is_err());
        assert!(response.total_pressure_jump(cx).is_err());
        assert!(receivers.evaluate(&response, cx).is_err());
        assert!(
            prepared
                .prepare_molecular_illumination(&gas, OMEGA, &incident, controls(), cx)
                .is_err()
        );
    });
}
