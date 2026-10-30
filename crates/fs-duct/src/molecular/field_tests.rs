use super::*;
use crate::HoleState;
use crate::profile::gas_field::tests::field;
use fs_matdb::QueryPoint;
use fs_material::gas_field::PrimitiveGasState;
use fs_material::relaxation::{
    InternalReservoir, RelaxationProcess, RelaxationRateLaw, SpeciesRelaxation, VibrationalMode,
};
use fs_thermochem::SpeciesId;

fn law(active: bool) -> RelaxationAbsorption {
    RelaxationAbsorption::new(
        ["analytic-field-a", "analytic-field-b"]
            .into_iter()
            .enumerate()
            .map(|(i, name)| SpeciesRelaxation {
                species: SpeciesId::new(name).unwrap(),
                processes: if active && i == 0 {
                    vec![RelaxationProcess {
                        reservoir: InternalReservoir::Vibrational(vec![VibrationalMode {
                            theta_k: 900.0,
                            degeneracy: 1,
                        }]),
                        rate: RelaxationRateLaw::PartnerIndependent {
                            p_tau_ref_pa_s: 20.0,
                            t_ref_k: 300.0,
                            exponent: 1.2,
                            t_min_k: 200.0,
                            t_max_k: 500.0,
                        },
                    }]
                } else {
                    Vec::new()
                },
            })
            .collect(),
    )
    .unwrap()
}

fn cylinder(radius: f64, field: &AxialGasField) -> Duct {
    Duct {
        segments: vec![Segment::Cylinder {
            radius,
            length: field.length_m(),
        }],
    }
}

fn bind(
    duct: &Duct,
    field: &AxialGasField,
    law: &RelaxationAbsorption,
    poles: Option<usize>,
    steps: usize,
) -> MolecularDuct {
    MolecularDuct::from_gas_field(
        duct,
        field,
        law,
        poles.map(|poles_per_field| MolecularWallResolution {
            poles_per_field,
            band_rad_s: [100.0, 4000.0],
            maximum_relative_error: 0.5,
        }),
        RiccatiSteps(steps),
        ConeSlices::Fixed(8),
        4096,
    )
    .unwrap()
}

fn relative(a: C64, b: C64) -> f64 {
    (a - b).abs() / b.abs()
}

// Independent radial heat/momentum equation: 2 I1(z)/(z I0(z)), evaluated
// from the convergent power series in z^2. This consumes no diffusion poles.
fn radial_average(z_squared: C64) -> C64 {
    let u = z_squared.scale(0.25);
    let mut term = C64::ONE;
    let mut i0 = C64::ONE;
    let mut i1 = C64::ONE;
    for n in 1..100 {
        term = (term * u).scale(1.0 / (n * n) as f64);
        i0 = i0 + term;
        i1 = i1 + term.scale(1.0 / (n + 1) as f64);
        if term.abs() < 1e-17 * i0.abs() {
            break;
        }
    }
    i1 / i0
}

// Integrate the physical mass and momentum equations backward from a finite
// terminal load with a separate classical RK4 rule. The local constitutive
// law is assembled from reservoir capacities/times and the radial heat PDE;
// no production wave, Maxwell spectrum or Magnus step enters this oracle.
fn oracle(
    source: &AxialGasField,
    material: &RelaxationAbsorption,
    radius: f64,
    omega: f64,
    load: C64,
    wall: bool,
) -> (C64, C64, C64) {
    let gas = *source.state_at(0.0).unwrap().state();
    let area = core::f64::consts::PI * radius * radius;
    let reference = gas.characteristic_impedance / area;
    let s = C64::new(0.0, -omega);
    let rhs = |x: f64, y: [C64; 2]| {
        let state = evaluate(source, material, x.max(0.0).min(source.length_m())).unwrap();
        let gas = state.gas();
        let cv = state
            .processes()
            .iter()
            .fold(C64::from_re(state.frozen_cv_j_kg_k()), |sum, p| {
                sum + C64::from_re(p.heat_capacity_j_kg_k)
                    / (C64::ONE + s.scale(p.relaxation_time_s))
            });
        let cp = cv + C64::from_re(gas.specific_gas_constant);
        let (rho, compliance) = if wall {
            let heat = radial_average(
                (s * cp).scale(gas.density * radius * radius / gas.thermal_conductivity),
            );
            let shear =
                radial_average(s.scale(gas.density * radius * radius / gas.dynamic_viscosity));
            (
                C64::from_re(gas.density) / (C64::ONE - shear),
                (C64::ONE - (C64::ONE - heat) / cp.scale(1.0 / gas.specific_gas_constant))
                    .scale(1.0 / gas.pressure),
            )
        } else {
            (
                C64::from_re(gas.density),
                (cv / cp).scale(1.0 / gas.pressure),
            )
        };
        [
            (-s * rho).scale(1.0 / (area * reference)) * y[1],
            (-s * compliance).scale(area * reference) * y[0],
        ]
    };
    let n = 2048;
    let h = -source.length_m() / n as f64;
    let mut y = [load.scale(1.0 / reference), C64::ONE];
    for j in 0..n {
        let x = source.length_m() * (1.0 - j as f64 / n as f64);
        let a = rhs(x, y);
        let b = rhs(
            x + 0.5 * h,
            [y[0] + a[0].scale(0.5 * h), y[1] + a[1].scale(0.5 * h)],
        );
        let c = rhs(
            x + 0.5 * h,
            [y[0] + b[0].scale(0.5 * h), y[1] + b[1].scale(0.5 * h)],
        );
        let d = rhs(x + h, [y[0] + c[0].scale(h), y[1] + c[1].scale(h)]);
        for i in 0..2 {
            y[i] = y[i] + (a[i] + b[i].scale(2.0) + c[i].scale(2.0) + d[i]).scale(h / 6.0);
        }
    }
    ((y[0] / y[1]).scale(reference), y[1].recip(), load / y[1])
}

#[test]
fn uniform_source_field_preserves_molecular_chain_bitwise() {
    let source = field([310.0, 310.0], [0.8, 0.8], 101325.0);
    let material = law(true);
    let mut duct = cylinder(0.0007, &source);
    duct.segments.push(Segment::ToneHole {
        hole_radius: 0.0003,
        bore_radius: 0.0007,
        chimney_height: 0.004,
        state: HoleState::Open,
    });
    let resolution = MolecularWallResolution {
        poles_per_field: 8,
        band_rad_s: [100.0, 4000.0],
        maximum_relative_error: 0.5,
    };
    let direct = MolecularDuct::new(
        &duct,
        evaluate(&source, &material, 0.0).unwrap(),
        Some(resolution),
        ConeSlices::Fixed(8),
        4096,
    )
    .unwrap();
    let varying = bind(&duct, &source, &material, Some(8), 13);
    for omega in [200.0, 1000.0, 3000.0] {
        let a = direct
            .response_with_holes(
                omega,
                Some(C64::new(2e6, -3e5)),
                &[C64::new(4e6, -1e6)],
                ComposeOptions::default(),
            )
            .unwrap();
        let b = varying
            .response_with_holes(
                omega,
                Some(C64::new(2e6, -3e5)),
                &[C64::new(4e6, -1e6)],
                ComposeOptions::default(),
            )
            .unwrap();
        assert_eq!(a.duct.z_in, b.duct.z_in);
        assert_eq!(a.duct.u_out_over_u_in, b.duct.u_out_over_u_in);
        assert_eq!(
            a.holes[0].p_junction_over_u_in,
            b.holes[0].p_junction_over_u_in
        );
        assert_eq!(a.holes[0].u_mouth_over_u_in, b.holes[0].u_mouth_over_u_in);
        assert_eq!(
            a.holes[0].dissipated_power_per_squared_inlet_flow,
            b.holes[0].dissipated_power_per_squared_inlet_flow
        );
    }
    assert_eq!(varying.material_states(), 1);
    assert_eq!(varying.material_points(), direct.material_points());
    assert_eq!(varying.gas_field(), Some(&source));
    assert_eq!(varying.field_relaxation_law(), Some(&material));
}

#[test]
fn molecular_temperature_and_composition_field_has_fourth_order_physical_limit() {
    let source = field([290.0, 390.0], [0.99, 0.45], 101325.0);
    let material = law(true);
    let radius = 0.001;
    let duct = cylinder(radius, &source);
    let omega = 4700.0;
    let load = C64::new(2e7, -1e6);
    let exact = oracle(&source, &material, radius, omega, load, false);
    let errors = [8, 16, 32].map(|n| {
        let response = bind(&duct, &source, &material, None, n)
            .response(omega, Some(load), ComposeOptions::default())
            .unwrap();
        relative(response.z_in, exact.0)
    });
    for pair in errors.windows(2) {
        assert!(
            pair[0] / pair[1] > 13.0 && pair[0] / pair[1] < 19.0,
            "{errors:?}"
        );
    }
    assert!(errors[2] < 2e-5, "{errors:?}");
    let fine = bind(&duct, &source, &material, None, 128);
    let response = fine
        .response(omega, Some(load), ComposeOptions::default())
        .unwrap();
    assert!(relative(response.u_out_over_u_in, exact.1) < 1e-7);
    assert!(relative(response.p_out_over_u_in, exact.2) < 1e-7);
    let homogeneous = MolecularDuct::new(
        &duct,
        evaluate(&source, &material, 0.0).unwrap(),
        None,
        ConeSlices::Fixed(1),
        4096,
    )
    .unwrap();
    assert!(
        relative(
            homogeneous
                .response(omega, Some(load), ComposeOptions::default())
                .unwrap()
                .z_in,
            response.z_in
        ) > 0.01
    );
    assert!(response.z_in.re > 0.0);
    assert_eq!(
        fine.outlet_state(),
        &evaluate(&source, &material, source.length_m()).unwrap()
    );
}

#[test]
fn joint_wall_and_molecular_field_converges_to_radial_heat_equation() {
    let source = field([290.0, 370.0], [0.95, 0.5], 101325.0);
    let material = law(true);
    let radius = 0.0007;
    let duct = cylinder(radius, &source);
    let omega = 1900.0;
    let load = C64::new(4e7, -2e6);
    let exact = oracle(&source, &material, radius, omega, load, true);
    let axial = [4, 8, 16].map(|steps| {
        let response = bind(&duct, &source, &material, Some(8), steps)
            .response(omega, Some(load), ComposeOptions::default())
            .unwrap();
        relative(response.z_in, exact.0)
    });
    for pair in axial.windows(2) {
        assert!(pair[0] / pair[1] > 12.0, "axial={axial:?}");
    }
    let wall = [2, 4, 8].map(|poles| {
        let response = bind(&duct, &source, &material, Some(poles), 128)
            .response(omega, Some(load), ComposeOptions::default())
            .unwrap();
        relative(response.z_in, exact.0)
    });
    assert!(
        wall[1] < wall[0] * 0.1 && wall[2] < wall[1] * 0.1,
        "wall={wall:?}"
    );
    assert!(wall[2] < 1e-7, "wall={wall:?}");
    let equilibrium = bind(&duct, &source, &law(false), Some(8), 64)
        .response(omega, Some(load), ComposeOptions::default())
        .unwrap();
    let slip = bind(&duct, &source, &material, None, 64)
        .response(omega, Some(load), ComposeOptions::default())
        .unwrap();
    assert!(relative(equilibrium.z_in, exact.0) > 1e-3);
    assert!(relative(slip.z_in, exact.0) > 1e-2);
}

#[test]
fn molecular_field_chimney_uses_its_exact_local_gas_and_heat() {
    let source = field([290.0, 390.0], [0.99, 0.45], 101325.0);
    let material = law(true);
    let x = 0.4 * source.length_m();
    let duct = Duct {
        segments: vec![
            Segment::Cylinder {
                radius: 0.001,
                length: x,
            },
            Segment::ToneHole {
                hole_radius: 0.0004,
                bore_radius: 0.001,
                chimney_height: 0.012,
                state: HoleState::Vent(0.6),
            },
            Segment::Cylinder {
                radius: 0.001,
                length: source.length_m() - x,
            },
        ],
    };
    let model = bind(&duct, &source, &material, Some(8), 16);
    let hole = model
        .pieces
        .iter()
        .find_map(|piece| match piece {
            Piece::ToneHole(hole) => Some(hole),
            _ => None,
        })
        .unwrap();
    assert_eq!(
        model.state_at_point(hole.point),
        &evaluate(&source, &material, x).unwrap()
    );
    assert_ne!(model.state_at_point(hole.point), model.state());
    let response = model
        .response_with_holes(
            2200.0,
            Some(C64::new(2e7, -1e6)),
            &[C64::new(4e6, -1e6)],
            ComposeOptions::default(),
        )
        .unwrap();
    assert!(response.duct.z_in.re > 0.0);
    assert_eq!(response.holes.len(), 1);
    assert!(response.holes[0].dissipated_power_per_squared_inlet_flow > 0.0);
    assert!(response.holes[0].mouth_power_per_squared_inlet_flow > 0.0);
}

#[test]
fn field_binding_reuses_constant_spans_and_admits_actual_work() {
    let original = field([290.0, 390.0], [0.99, 0.45], 101325.0);
    let primitive: Vec<PrimitiveGasState> = vec![
        original.primitive_states()[0].clone(),
        original.primitive_states()[0].clone(),
        original.primitive_states()[1].clone(),
        original.primitive_states()[1].clone(),
    ];
    let length = original.length_m();
    let source = AxialGasField::new(
        vec![0.0, length / 3.0, 2.0 * length / 3.0, length],
        original.static_pressure_pa(),
        primitive,
        original.model().clone(),
        QueryPoint::new(),
        original.phase(),
    )
    .unwrap();
    let material = law(true);
    let duct = cylinder(0.001, &source);
    let model = bind(&duct, &source, &material, None, 8);
    assert_eq!(model.blocks(), 24);
    assert_eq!(model.material_states(), 18);
    assert_eq!(model.material_points(), 18);
    assert!(
        MolecularDuct::from_gas_field(
            &duct,
            &source,
            &material,
            None,
            RiccatiSteps(8),
            ConeSlices::Fixed(1),
            23
        )
        .is_err()
    );
    let mut polls = 0;
    assert!(matches!(
        MolecularDuct::from_gas_field_controlled(
            &duct,
            &source,
            &material,
            None,
            RiccatiSteps(8),
            ConeSlices::Fixed(1),
            4096,
            &mut || {
                polls += 1;
                polls == 15
            }
        ),
        Err(MolecularDuctError::Duct(DuctError::Cancelled { .. }))
    ));
    let short = Duct {
        segments: vec![Segment::Cylinder {
            radius: 0.001,
            length: length * 0.9,
        }],
    };
    assert!(
        MolecularDuct::from_gas_field(
            &short,
            &source,
            &material,
            None,
            RiccatiSteps(8),
            ConeSlices::Fixed(1),
            4096
        )
        .is_err()
    );
}

#[test]
fn uniform_molecular_field_wavefront_retains_direct_arithmetic() {
    let source = field([300.0, 300.0], [0.75, 0.75], 101325.0);
    let duct = Duct {
        segments: vec![
            Segment::Cylinder {
                radius: 0.0073,
                length: 0.21,
            },
            Segment::ToneHole {
                hole_radius: 0.002,
                bore_radius: 0.0073,
                chimney_height: 0.004,
                state: HoleState::Vent(0.4),
            },
            Segment::Cylinder {
                radius: 0.0073,
                length: source.length_m() - 0.21,
            },
        ],
    };
    let model = bind(&duct, &source, &law(true), None, 8);
    let times = model.field_travel_times_s().unwrap();
    let speed = model.state().frozen_sound_speed_m_s();
    let site = model.open_hole_sites()[0];
    assert_eq!(times.terminal_s, source.length_m() / speed);
    assert_eq!(
        times.open_mouths_s,
        vec![(0.21 + site.open_chimney_length_m) / speed]
    );
    let direct = MolecularDuct::new(
        &duct,
        model.state().clone(),
        None,
        ConeSlices::Fixed(8),
        4096,
    )
    .unwrap();
    assert!(direct.field_travel_times_s().is_none());
}

#[test]
fn molecular_field_wavefront_integrates_frozen_slowness_and_local_chimney() {
    let source = field([280.0, 410.0], [0.80, 0.98], 101325.0);
    let material = law(true);
    let duct = Duct {
        segments: vec![
            Segment::Cylinder {
                radius: 0.0073,
                length: 0.21,
            },
            Segment::ToneHole {
                hole_radius: 0.002,
                bore_radius: 0.0073,
                chimney_height: 0.004,
                state: HoleState::Open,
            },
            Segment::Cone {
                inlet_radius: 0.0073,
                outlet_radius: 0.009,
                length: source.length_m() - 0.21,
            },
        ],
    };
    // Independent composite Simpson integration of the gas-frozen ideal-gas
    // acoustic modulus. It uses fresh source evaluations, no prepared points.
    let slowness = |x| {
        let state = evaluate(&source, &material, x).unwrap();
        let gas = state.gas();
        let cv = state.frozen_cv_j_kg_k();
        let modulus = gas.pressure * (1.0 + gas.specific_gas_constant / cv);
        (gas.density / modulus).sqrt()
    };
    let integral = |length: f64| {
        let n = 4096;
        let h = length / n as f64;
        let mut sum = slowness(0.0) + slowness(length);
        for i in 1..n {
            sum += (if i % 2 == 0 { 2.0 } else { 4.0 }) * slowness(i as f64 * h);
        }
        sum * h / 3.0
    };
    let exact = integral(source.length_m());
    let exact_prefix = integral(0.21);
    let mut errors = Vec::new();
    for steps in [4, 8, 16] {
        let model = bind(&duct, &source, &material, None, steps);
        let times = model.field_travel_times_s().unwrap();
        let site = model.open_hole_sites()[0];
        let mouth = exact_prefix + site.open_chimney_length_m * slowness(0.21);
        errors.push((times.terminal_s - exact).abs());
        assert!((times.open_mouths_s[0] - mouth).abs() < 2e-10);
    }
    assert!(
        errors[0] > 1e-12,
        "wavefront refinement must be live: {errors:?}"
    );
    for pair in errors.windows(2) {
        let order = (pair[0] / pair[1]).log2();
        assert!((3.8..4.2).contains(&order), "{errors:?}, order={order}");
    }
    assert!(errors[2] / exact < 2e-8, "{errors:?}");
    let inlet_only = source.length_m() * slowness(0.0);
    let equilibrium = GasFieldDuct::new(&duct, &source, RiccatiSteps(64))
        .unwrap()
        .one_way_time_s();
    assert!(
        (inlet_only - exact).abs() / exact > 0.01,
        "inlet={inlet_only}, integrated={exact}"
    );
    assert!((equilibrium - exact).abs() / exact > 0.001);
}
