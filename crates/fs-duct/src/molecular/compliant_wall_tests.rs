use super::*;
use crate::network::CoupledTermination;
use crate::profile::gas_field::tests::field;
use crate::{HoleState, RadiusProfile};
use fs_material::relaxation::{
    InternalReservoir, RelaxationProcess, RelaxationRateLaw, RelaxingMedium, SpeciesRelaxation,
    VibrationalMode,
};
use fs_thermochem::SpeciesId;

// Analytic material fixture, not a fit to a measured specimen. The cylindrical
// hoop stiffness scales with r^-2, so a changed physical radius must be live.
#[derive(Clone, Copy)]
struct Surface;
impl WallLaw for Surface {
    fn specific_impedance(&self, radius: f64, omega: f64) -> Result<C64, DuctError> {
        Ok(C64::new(
            90.0,
            600.0 / (radius * radius * omega) - 0.12 * omega,
        ))
    }
}
fn law() -> RelaxationAbsorption {
    RelaxationAbsorption::new(
        ["analytic-field-a", "analytic-field-b"]
            .into_iter()
            .map(|name| SpeciesRelaxation {
                species: SpeciesId::new(name).unwrap(),
                processes: vec![RelaxationProcess {
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
                }],
            })
            .collect(),
    )
    .unwrap()
}
fn state(source: &AxialGasField, material: &RelaxationAbsorption, x: f64) -> RelaxationState {
    material
        .evaluate(&RelaxingMedium::from_thermochemical(&source.state_at(x).unwrap()).unwrap())
        .unwrap()
}
fn relative(a: C64, b: C64) -> f64 {
    (a - b).abs() / b.abs().max(1e-20)
}
#[track_caller]
fn near(a: C64, b: C64, tol: f64) {
    assert!(
        relative(a, b) < tol,
        "{a:?} != {b:?}; relative={}",
        relative(a, b)
    );
}

// Independent circular heat/momentum PDE solution: 2 I1(z)/(z I0(z)).
// No production wall diffusion spectrum, wave map or Magnus owner is called.
fn radial(z2: C64) -> C64 {
    let mut term = C64::ONE;
    let (mut i0, mut i1) = (C64::ONE, C64::ONE);
    for n in 1..120 {
        term = (term * z2).scale(0.25 / (n * n) as f64);
        i0 = i0 + term;
        i1 = i1 + term.scale(1.0 / (n + 1) as f64);
        if term.abs() < 1e-17 * i0.abs() {
            break;
        }
    }
    i1 / i0
}
fn coefficients(
    state: &RelaxationState,
    radius: f64,
    slant: f64,
    omega: f64,
    thermal: bool,
) -> (C64, C64) {
    let gas = state.gas();
    let s = C64::new(0.0, -omega);
    let cv = state
        .processes()
        .iter()
        .fold(C64::from_re(state.frozen_cv_j_kg_k()), |v, p| {
            v + C64::from_re(p.heat_capacity_j_kg_k) / (C64::ONE + s.scale(p.relaxation_time_s))
        });
    let cp = cv + C64::from_re(gas.specific_gas_constant);
    let (rho, beta) = if thermal {
        let shear = radial(s.scale(gas.density * radius * radius / gas.dynamic_viscosity));
        let heat = radial((s * cp).scale(gas.density * radius * radius / gas.thermal_conductivity));
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
    let area = core::f64::consts::PI * radius * radius;
    let wall = Surface
        .specific_impedance(radius, omega)
        .unwrap()
        .recip()
        .scale(core::f64::consts::TAU * radius * slant);
    ((s * rho).scale(1.0 / area), (s * beta).scale(area) + wall)
}
fn line(z: C64, y: C64, length: f64) -> [C64; 4] {
    let k = (-z * y).sqrt();
    let phase = k.scale(length);
    let (s, c) = det::sin_cos(phase.re);
    let e = det::exp(phase.im);
    let (ch, sh) = (0.5 * (e + 1.0 / e), 0.5 * (e - 1.0 / e));
    let (sin, cos) = (C64::new(s * ch, c * sh), C64::new(c * ch, -s * sh));
    [cos, z * sin / k, y * sin / k, cos]
}

#[test]
fn molecular_compliant_cylinder_matches_physical_equations_and_work() {
    let source = field([310.0; 2], [0.8; 2], 101325.0);
    let material = law();
    let state = state(&source, &material, 0.0);
    let (radius, length) = (0.0007, 0.057);
    let duct = Duct {
        segments: vec![Segment::Cylinder { radius, length }],
    };
    let load = C64::new(4e7, -2e7);
    for thermal in [false, true] {
        let model = MolecularDuct::new_with_compliant_wall(
            &duct,
            state.clone(),
            thermal.then_some(MolecularWallResolution {
                poles_per_field: 16,
                band_rad_s: [50.0, 4000.0],
                maximum_relative_error: 0.005,
            }),
            Arc::new(Surface),
            ConeSlices::Fixed(8),
            64,
        )
        .unwrap();
        for omega in [100.0, 700.0, 3000.0] {
            let (z, y) = coefficients(&state, radius, 1.0, omega, thermal);
            let [a, b, c, d] = line(z, y, length);
            let inlet = c * load + d;
            let response = model
                .response(omega, Some(load), ComposeOptions::default())
                .unwrap();
            near(response.z_in, (a * load + b) / inlet, 2e-8);
            near(response.u_out_over_u_in, inlet.recip(), 2e-8);
            let mut integrated = 0.0;
            let n = 512;
            for j in 0..=n {
                let x = length * j as f64 / n as f64;
                let [a, b, c, d] = line(z, y, length - x);
                let p = (a * load + b) / inlet;
                let u = (c * load + d) / inlet;
                let weight = if j == 0 || j == n {
                    1.0
                } else if j % 2 == 0 {
                    2.0
                } else {
                    4.0
                };
                integrated += weight * 0.5 * (z.re * u.abs().powi(2) + y.re * p.abs().powi(2));
            }
            integrated *= length / (3 * n) as f64;
            let work = 0.5
                * (response.z_in.re
                    - (response.p_out_over_u_in * response.u_out_over_u_in.conj()).re);
            assert!(work > 0.0 && (work - integrated).abs() < 3e-7 * integrated);
        }
    }
}

// A separate classical RK4 solve of d[p,U]/dx = -[0,Z';Y',0][p,U].
// Geometry breaks are integrated separately; the wall area is the physical
// surface area, and the primitive field re-evaluates the real material owner.
fn oracle(
    source: &AxialGasField,
    material: &RelaxationAbsorption,
    duct: &Duct,
    omega: f64,
    load: C64,
) -> C64 {
    let reference = 1e7;
    let mut q = [load.scale(1.0 / reference), C64::ONE];
    let mut end = source.length_m();
    for segment in duct.segments.iter().rev() {
        let length = segment.axial_length();
        let start = end - length;
        let rhs = |x: f64, q: [C64; 2]| {
            let local = x.clamp(0.0, length);
            let (radius, slope) = match segment {
                Segment::Cone {
                    inlet_radius,
                    outlet_radius,
                    ..
                } => (
                    inlet_radius + (outlet_radius - inlet_radius) * local / length,
                    (outlet_radius - inlet_radius) / length,
                ),
                Segment::Profile { profile, .. } => (
                    profile.radius(local).unwrap(),
                    profile.slope(local).unwrap(),
                ),
                Segment::Cylinder { radius, .. } => (*radius, 0.0),
                _ => unreachable!(),
            };
            let material = state(
                source,
                material,
                (start + local).clamp(0.0, source.length_m()),
            );
            let (z, y) = coefficients(
                &material,
                radius,
                det::sqrt(1.0 + slope * slope),
                omega,
                false,
            );
            [
                (-z).scale(1.0 / reference) * q[1],
                (-y).scale(reference) * q[0],
            ]
        };
        let n = 2048;
        let h = -length / n as f64;
        for j in 0..n {
            let x = length * (1.0 - j as f64 / n as f64);
            let a = rhs(x, q);
            let b = rhs(
                x + h / 2.0,
                [q[0] + a[0].scale(h / 2.0), q[1] + a[1].scale(h / 2.0)],
            );
            let c = rhs(
                x + h / 2.0,
                [q[0] + b[0].scale(h / 2.0), q[1] + b[1].scale(h / 2.0)],
            );
            let d = rhs(x + h, [q[0] + c[0].scale(h), q[1] + c[1].scale(h)]);
            for i in 0..2 {
                q[i] = q[i] + (a[i] + b[i].scale(2.0) + c[i].scale(2.0) + d[i]).scale(h / 6.0);
            }
        }
        end = start;
    }
    (q[0] / q[1]).scale(reference)
}

#[test]
fn compliant_cones_refine_without_thermal_diffusion_and_keep_distinct_slopes() {
    let source = field([310.0; 2], [0.8; 2], 101325.0);
    let material = law();
    let length = source.length_m();
    let duct = Duct {
        segments: vec![
            Segment::Cone {
                inlet_radius: 1.0 / 256.0,
                outlet_radius: 1.0 / 128.0,
                length: length * 0.8,
            },
            Segment::Cone {
                inlet_radius: 1.0 / 128.0,
                outlet_radius: 1.0 / 256.0,
                length: length * 0.2,
            },
        ],
    };
    let load = C64::new(2e6, -1e6);
    let omega = 1100.0;
    let expected = oracle(&source, &material, &duct, omega, load);
    let mut errors = Vec::new();
    for n in [8, 16, 32] {
        let model = MolecularDuct::new_with_compliant_wall(
            &duct,
            state(&source, &material, 0.0),
            None,
            Arc::new(Surface),
            ConeSlices::Fixed(n),
            64,
        )
        .unwrap();
        assert_eq!(model.blocks(), 2 * n);
        let actual = model
            .response(omega, Some(load), ComposeOptions::default())
            .unwrap()
            .z_in;
        errors.push(relative(actual, expected));
    }
    assert!(
        errors[0] > 1e-8 && errors[2] < errors[0] / 12.0 && errors[2] < 3e-4,
        "{errors:?}"
    );
    let rigid = MolecularDuct::new(
        &duct,
        state(&source, &material, 0.0),
        None,
        ConeSlices::Fixed(32),
        64,
    )
    .unwrap();
    assert_eq!(rigid.blocks(), 2);
    assert!(
        relative(
            rigid
                .response(omega, Some(load), ComposeOptions::default())
                .unwrap()
                .z_in,
            expected
        ) > 0.01
    );
}

#[test]
fn compliant_primitive_field_uses_original_profile_quadrature_and_uniform_bits() {
    let material = law();
    let source = field([285.0, 390.0], [0.99, 0.45], 101325.0);
    let profile = Arc::new(RadiusProfile::exponential(0.002, 8.0, source.length_m()).unwrap());
    let duct = Duct {
        segments: vec![Segment::Profile {
            profile,
            steps: RiccatiSteps(16),
        }],
    };
    let load = C64::new(2e6, -1e6);
    let omega = 1200.0;
    let expected = oracle(&source, &material, &duct, omega, load);
    let mut errors = Vec::new();
    for n in [16, 32, 64] {
        let model = MolecularDuct::from_gas_field_with_compliant_wall(
            &duct,
            &source,
            &material,
            None,
            Arc::new(Surface),
            RiccatiSteps(n),
            ConeSlices::Fixed(8),
            128,
        )
        .unwrap();
        errors.push(relative(
            model
                .response(omega, Some(load), ComposeOptions::default())
                .unwrap()
                .z_in,
            expected,
        ));
    }
    let order = (errors[1] / errors[2]).log2();
    assert!(
        errors[0] > 1e-7 && (3.8..4.2).contains(&order) && errors[2] < 2e-5,
        "{errors:?}, order={order}"
    );
    let uniform = field([310.0; 2], [0.8; 2], 101325.0);
    let direct = MolecularDuct::new_with_compliant_wall(
        &duct,
        state(&uniform, &material, 0.0),
        None,
        Arc::new(Surface),
        ConeSlices::Fixed(8),
        128,
    )
    .unwrap();
    let field = MolecularDuct::from_gas_field_with_compliant_wall(
        &duct,
        &uniform,
        &material,
        None,
        Arc::new(Surface),
        RiccatiSteps(3),
        ConeSlices::Fixed(8),
        128,
    )
    .unwrap();
    let a = direct
        .response(omega, Some(load), ComposeOptions::default())
        .unwrap();
    let b = field
        .response(omega, Some(load), ComposeOptions::default())
        .unwrap();
    assert_eq!(a.z_in, b.z_in);
    assert_eq!(a.p_out_over_u_in, b.p_out_over_u_in);
    assert_eq!(a.u_out_over_u_in, b.u_out_over_u_in);
}

#[test]
fn compliant_chimneys_and_coupled_exterior_retain_physical_work_and_flows() {
    let source = field([310.0; 2], [0.8; 2], 101325.0);
    let material = law();
    let state = state(&source, &material, 0.0);
    let mut segments = Vec::new();
    for opening in [HoleState::Open, HoleState::Closed, HoleState::Vent(0.3)] {
        segments.push(Segment::Cylinder {
            radius: 0.004,
            length: 0.15,
        });
        segments.push(Segment::ToneHole {
            hole_radius: 0.0016,
            bore_radius: 0.004,
            chimney_height: 0.005,
            state: opening,
        });
    }
    segments.push(Segment::Cylinder {
        radius: 0.004,
        length: 0.12,
    });
    let model = MolecularDuct::new_with_compliant_wall(
        &Duct { segments },
        state.clone(),
        None,
        Arc::new(Surface),
        ConeSlices::Fixed(8),
        64,
    )
    .unwrap();
    let omega = 1300.0;
    let mouths = [C64::new(2e6, -1e6), C64::new(3e6, -5e5)];
    let terminal = C64::new(4e5, -2e5);
    let scalar = model
        .response_with_holes(omega, Some(terminal), &mouths, ComposeOptions::default())
        .unwrap();
    let labels = vec!["open".into(), "vent".into(), "bell".into()];
    let prepared = model
        .prepare_exterior_ports(omega, &labels, ComposeOptions::default())
        .unwrap();
    let mut matrix = CoupledTermination {
        omega,
        labels,
        z: vec![C64::ZERO; 9],
        source_id: "analytic-wall-test".into(),
    };
    matrix.z[0] = mouths[0];
    matrix.z[4] = mouths[1].scale(1.0 / 0.3);
    matrix.z[8] = terminal;
    let coupled = prepared.response(&matrix).unwrap();
    near(coupled.z_in, scalar.duct.z_in, 1e-10);
    near(
        coupled.exterior_flows_over_u_in[2],
        scalar.duct.u_out_over_u_in,
        1e-10,
    );
    for (hole, load) in [(&scalar.holes[0], mouths[0]), (&scalar.holes[2], mouths[1])] {
        let (z, y) = coefficients(&state, hole.site.hole_radius_m, 1.0, omega, false);
        let [a, b, c, d] = line(z, y, hole.site.open_chimney_length_m);
        let pressure = hole.p_junction_over_u_in;
        let physical = pressure / (a * load + b).scale(1.0 / hole.site.open_fraction);
        near(hole.u_mouth_over_u_in, physical, 1e-10);
        near(hole.u_open_in_over_u_in, (c * load + d) * physical, 1e-10);
        assert!(hole.dissipated_power_per_squared_inlet_flow > 0.0);
    }
    assert_eq!(scalar.holes[1].u_mouth_over_u_in, C64::ZERO);
    assert!(coupled.dissipated_power_per_squared_inlet_flow > 0.0);
    matrix.z[1] = C64::new(2e5, -1e5);
    matrix.z[3] = matrix.z[1];
    let mutual = prepared.response(&matrix).unwrap();
    assert!(relative(mutual.z_in, coupled.z_in) > 1e-4);
    assert!(
        mutual.exterior_power_per_squared_inlet_flow > 0.0
            && mutual.dissipated_power_per_squared_inlet_flow > 0.0
    );
}

#[test]
fn compliant_wall_admission_cancellation_caps_and_cache_are_live() {
    use std::sync::atomic::{AtomicUsize, Ordering};
    struct Counted(AtomicUsize);
    impl WallLaw for Counted {
        fn specific_impedance(&self, radius: f64, omega: f64) -> Result<C64, DuctError> {
            self.0.fetch_add(1, Ordering::Relaxed);
            if omega > 2000.0 {
                return Err(bad("analytic wall source frequency domain"));
            }
            Surface.specific_impedance(radius, omega)
        }
    }
    struct Active;
    impl WallLaw for Active {
        fn specific_impedance(&self, _: f64, _: f64) -> Result<C64, DuctError> {
            Ok(C64::new(-1.0, 10.0))
        }
    }
    let source = field([310.0; 2], [0.8; 2], 101325.0);
    let state = state(&source, &law(), 0.0);
    let duct = Duct {
        segments: vec![
            Segment::Cylinder {
                radius: 0.004,
                length: 0.01
            };
            16
        ],
    };
    let wall = Arc::new(Counted(AtomicUsize::new(0)));
    let model = MolecularDuct::new_with_compliant_wall(
        &duct,
        state.clone(),
        None,
        wall.clone(),
        ConeSlices::Fixed(8),
        16,
    )
    .unwrap();
    assert_eq!(wall.0.load(Ordering::Relaxed), 0);
    model
        .response(1000.0, Some(C64::ZERO), ComposeOptions::default())
        .unwrap();
    assert_eq!(wall.0.load(Ordering::Relaxed), 1);
    assert!(
        model
            .response(3000.0, Some(C64::ZERO), ComposeOptions::default())
            .is_err()
    );
    assert!(
        model
            .response_controlled(
                1000.0,
                Some(C64::ZERO),
                ComposeOptions::default(),
                &mut || true
            )
            .is_err()
    );
    let active = MolecularDuct::new_with_compliant_wall(
        &duct,
        state.clone(),
        None,
        Arc::new(Active),
        ConeSlices::Fixed(8),
        16,
    )
    .unwrap();
    assert!(
        active
            .response(1000.0, Some(C64::ZERO), ComposeOptions::default())
            .is_err()
    );
    let cone = Duct {
        segments: vec![Segment::Cone {
            inlet_radius: 0.004,
            outlet_radius: 0.008,
            length: 0.2,
        }],
    };
    assert!(
        MolecularDuct::new_with_compliant_wall(
            &cone,
            state.clone(),
            None,
            Arc::new(Surface),
            ConeSlices::Fixed(8),
            7
        )
        .is_err()
    );
    assert!(
        MolecularDuct::new_with_compliant_wall_controlled(
            &cone,
            state,
            None,
            Arc::new(Surface),
            ConeSlices::Fixed(8),
            8,
            &mut || true
        )
        .is_err()
    );
}
