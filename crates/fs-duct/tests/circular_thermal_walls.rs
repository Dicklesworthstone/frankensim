//! G1/G3: a regular circular axis, native Bessel reference and material-wall work.
//! Material values below are explicitly synthetic, not empirical solid data.
use fs_duct::DuctError;
use fs_duct::fluid::annular::thermal_wall::AnnularThermalWall as Wall;
use fs_duct::fluid::annular::thermal_wall::conduction::{
    RadialThermalWallAdmission, RadialThermalWallError, RadialWallBath as Bath,
    RadialWallFace as Face,
};
use fs_duct::fluid::annular::{AnnularAdmission, AnnularSection};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::gas::{GasSpec, GasState};
use fs_material::relaxation::{
    InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw, RelaxationState,
    RelaxingMedium, SpeciesRelaxation, VibrationalMode,
};
use fs_material::thermal::transport::{PreparedThermalTransport, ThermalTransportSourceQuery};
use fs_math::c64::C64;
use fs_qty::{Dims, QuantitySpec};
use fs_thermochem::SpeciesId;
use std::f64::consts::PI;

fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}
fn admission() -> AnnularAdmission {
    AnnularAdmission {
        maximum_transverse_phase: 0.8,
        maximum_cell_diffusion_phase: 0.9,
        maximum_relative_residual: 1e-7,
        maximum_payload_bytes: 1_000_000,
    }
}
fn solid_admission() -> RadialThermalWallAdmission {
    RadialThermalWallAdmission {
        maximum_cell_diffusion_phase: 0.9,
        maximum_relative_residual: 1e-7,
        maximum_payload_bytes: 1_000_000,
    }
}
fn solid_source(k: f64, capacity: f64) -> PreparedThermalTransport {
    let keys = [
        PropertyKey::new("synthetic-density", Dims([-3, 1, 0, 0, 0, 0])),
        PropertyKey::new(
            "synthetic-constrained-heat-capacity",
            Dims([2, 0, -2, -1, 0, 0]),
        ),
        PropertyKey::new("synthetic-radial-conductivity", Dims([1, 1, -3, -1, 0, 0])),
    ];
    let mut claims = ClaimSet::new();
    for (key, value) in keys.iter().zip([2500.0, capacity / 2500.0, k]) {
        claims
            .insert_claim(PropertyClaim {
                key: key.clone(),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Dims([0, 0, 0, 1, 0, 0]),
                    knots: vec![(200.0, value), (500.0, value)],
                    dims: key.quantity().dims(),
                },
                validity: ValidityDomain::unconstrained().with("T", 200.0, 500.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic radial thermal wall".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic thermal shell".into(),
            phase: "solid".into(),
            process: "radial source fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let point = |t| {
        QueryPoint::new()
            .with_quantity("T", QuantitySpec::dimensional(Dims([0, 0, 0, 1, 0, 0])), t)
            .unwrap()
    };
    PreparedThermalTransport::try_new(
        &card,
        ThermalTransportSourceQuery {
            density: &keys[0],
            specific_heat: &keys[1],
            conductivity: &keys[2],
            temperature_axis: "T",
            lower: &point(200.0),
            upper: &point(500.0),
            selections: [ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly); 3],
        },
    )
    .unwrap()
}
fn near(actual: C64, expected: C64, tol: f64) {
    assert!(
        (actual - expected).abs() <= tol * expected.abs().max(1e-300),
        "{actual:?} != {expected:?}, relative={:e}",
        (actual - expected).abs() / expected.abs()
    );
}
fn relative(actual: f64, expected: f64, tol: f64) {
    assert!(
        (actual - expected).abs() <= tol * expected.abs().max(1e-300),
        "{actual:e} != {expected:e}"
    );
}
fn material(count: usize) -> RelaxationState {
    let species = SpeciesId::new("synthetic").unwrap();
    let law = RelaxationAbsorption::new(vec![SpeciesRelaxation {
        species: species.clone(),
        processes: (0..count)
            .map(|i| RelaxationProcess {
                reservoir: InternalReservoir::Vibrational(vec![VibrationalMode {
                    theta_k: 1300.0 + 600.0 * i as f64,
                    degeneracy: 1,
                }]),
                rate: RelaxationRateLaw::PartnerIndependent {
                    p_tau_ref_pa_s: 101325.0 * (0.002 + 0.01 * i as f64),
                    t_ref_k: 293.15,
                    exponent: 0.0,
                    t_min_k: 200.0,
                    t_max_k: 500.0,
                },
            })
            .collect(),
    }])
    .unwrap();
    law.evaluate(&RelaxingMedium::new(gas(), vec![(species, 1.0)]).unwrap())
        .unwrap()
}

use fs_duct::fluid::annular::thermal_wall::conduction::layered::{
    PreparedLayeredRadialThermalWall as LayeredWall, RadialThermalLayer as Layer,
};
use fs_duct::fluid::circular::{CircularSection, CircularThermalWall as CWall};
use fs_phs::{zk_one_minus_f, zwikker_kosten_f};

#[test]
fn regular_axis_recovers_poiseuille_profile_flow_and_viscous_work() {
    let gas = gas();
    let radius = 0.0015;
    let omega = 1e-5;
    let q = C64::new(2e-8, -1e-8);
    let area = PI * radius * radius;
    let expected_center = q.scale(2.0 / area);
    let expected_loss = 4.0 * gas.dynamic_viscosity / (PI * radius.powi(4)) * q.norm_sq();
    let mut previous = f64::INFINITY;
    for cells in [24, 48, 96] {
        let section = CircularSection::new(radius, cells).unwrap();
        let fields = section
            .prepare_thermal_wall_fields(&gas, omega, CWall::Boundary(Wall::Adiabatic), admission())
            .unwrap();
        let center = fields.field_in_cell(0, 0.0, C64::ZERO, q).unwrap();
        assert_eq!(center.radius_m, 0.0);
        assert!(center.axial_velocity_m_s.abs() > 1.9 * q.abs() / area);
        let error = (center.axial_velocity_m_s - expected_center).abs();
        assert!(error < previous * 0.31, "N={cells}, {error:e}/{previous:e}");
        previous = error;
        let mut recovered = C64::ZERO;
        for cell in 0..cells {
            let dr =
                section.node_radius_m(cell + 1).unwrap() - section.node_radius_m(cell).unwrap();
            for fraction in [0.5 - 0.5 / 3.0_f64.sqrt(), 0.5 + 0.5 / 3.0_f64.sqrt()] {
                let field = fields.field_in_cell(cell, fraction, C64::ZERO, q).unwrap();
                recovered = recovered + field.axial_velocity_m_s.scale(PI * field.radius_m * dr);
            }
        }
        near(recovered, q, 4e-13);
        assert_eq!(
            fields
                .field_in_cell(cells - 1, 1.0, C64::ZERO, q)
                .unwrap()
                .axial_velocity_m_s,
            C64::ZERO
        );
        if cells == 96 {
            near(center.axial_velocity_m_s, expected_center, 4e-4);
            relative(
                fields
                    .dissipation_per_axial_length(C64::ZERO, q)
                    .unwrap()
                    .gas_viscous_w_m,
                expected_loss,
                4e-4,
            );
        }
    }
    assert!(AnnularSection::new(0.0, radius, 96).is_err());
}

#[test]
fn isothermal_circular_wave_and_thermal_mean_converge_to_original_bessel_owner() {
    let gas = gas();
    let radius = 0.003;
    let area = PI * radius * radius;
    for omega in [30.0, 500.0, 6000.0] {
        let rv = radius * (gas.density * omega / gas.dynamic_viscosity).sqrt();
        let ft = zwikker_kosten_f(rv * gas.prandtl.sqrt()).unwrap();
        let z = C64::new(0.0, omega * gas.density / area) / zk_one_minus_f(rv).unwrap();
        let y = C64::new(
            0.0,
            omega * area / (gas.density * gas.sound_speed * gas.sound_speed),
        ) * (C64::ONE + ft.scale(gas.gamma - 1.0));
        let gamma = (z * y).sqrt();
        let expected_k = C64::new(gamma.im, gamma.re);
        let expected_z = (z / y).sqrt().conj();
        let expected_h = C64::ONE - ft;
        let mut previous = f64::INFINITY;
        for cells in [96, 192, 384] {
            let section = CircularSection::new(radius, cells).unwrap();
            let fields = section
                .prepare_thermal_wall_fields(
                    &gas,
                    omega,
                    CWall::Boundary(Wall::Isothermal),
                    admission(),
                )
                .unwrap();
            let error = (fields.wave().wavenumber - expected_k).abs() / expected_k.abs()
                + (fields.mean_temperature_response() - expected_h).abs();
            assert!(
                error < previous * 0.3,
                "omega={omega}, N={cells}: {error:e}/{previous:e}"
            );
            previous = error;
            relative(fields.wave().shear_number, rv, 3e-15);
            if cells == 384 {
                near(fields.wave().wavenumber, expected_k, 3e-4);
                near(fields.wave().characteristic_impedance, expected_z, 3e-4);
                near(fields.mean_temperature_response(), expected_h, 3e-4);
            }
        }
    }
}

// Independent regular cylindrical solution f=I0(sqrt(z)*r), from its
// defining series. The area mean is sum a_n/(n+1); no FE matrix is used.
fn regular_basis(z: C64, radius: f64) -> [C64; 3] {
    let q = z.scale(radius * radius / 4.0);
    let mut term = C64::ONE;
    let (mut value, mut derivative, mut mean) = (C64::ONE, C64::ZERO, C64::ONE);
    for n in 1..=100 {
        let f = n as f64;
        term = term * q / C64::from_re(f * f);
        value = value + term;
        derivative = derivative + term.scale(2.0 * f / radius);
        mean = mean + term.scale(1.0 / (f + 1.0));
    }
    [value, derivative, mean]
}
fn regular_thermal(gas: &GasState, cp: C64, omega: f64, radius: f64, wall: Wall) -> [C64; 3] {
    let basis = regular_basis(
        C64::new(cp.im, cp.re).scale(omega * gas.density / gas.thermal_conductivity),
        radius,
    );
    let a = match wall {
        Wall::Isothermal => basis[0].recip(),
        Wall::Adiabatic => C64::ZERO,
        Wall::Finite {
            conductance_w_m2_k: g,
            areal_heat_capacity_j_m2_k: c,
        } => {
            let y = C64::new(g, omega * c);
            y / (basis[1].scale(gas.thermal_conductivity) + y * basis[0])
        }
    };
    [
        C64::ONE - a * basis[2],
        C64::ONE - a,
        C64::ONE - a * basis[0],
    ]
}

#[test]
fn molecular_circular_thermal_boundary_matches_independent_regular_series_with_refinement() {
    let material = material(2);
    let radius = 0.002;
    let omega = 400.0;
    let cp = material.caloric_response(omega).unwrap().cp_j_kg_k;
    for wall in [
        Wall::Isothermal,
        Wall::Adiabatic,
        Wall::Finite {
            conductance_w_m2_k: 25.0,
            areal_heat_capacity_j_m2_k: 0.08,
        },
    ] {
        let expected = regular_thermal(material.gas(), cp, omega, radius, wall);
        let mut previous = f64::INFINITY;
        for cells in [32, 64, 128] {
            let section = CircularSection::new(radius, cells).unwrap();
            let fields = section
                .prepare_molecular_thermal_wall_fields(
                    &material,
                    omega,
                    CWall::Boundary(wall),
                    admission(),
                )
                .unwrap();
            let error = (fields.mean_temperature_response() - expected[0]).abs();
            if wall != Wall::Adiabatic {
                assert!(
                    error < previous * 0.3,
                    "{wall:?} N={cells}: {error:e}/{previous:e}"
                );
            } else {
                near(fields.mean_temperature_response(), C64::ONE, 2e-14);
            }
            previous = error;
            if cells == 128 {
                near(fields.mean_temperature_response(), expected[0], 5e-4);
                let scale = cp.recip().scale(1.0 / material.gas().density);
                let center = fields.field_in_cell(0, 0.0, C64::ONE, C64::ZERO).unwrap();
                near(center.gas_temperature_k, expected[1].conj() * scale, 7e-4);
                let actual = fields
                    .wall_observation(C64::ONE)
                    .unwrap()
                    .interface_temperature_k;
                if wall == Wall::Isothermal {
                    assert_eq!(actual, C64::ZERO);
                } else {
                    near(actual, expected[2].conj() * scale, 7e-4);
                }
            }
        }
    }
}

// Regular gas series followed by independent dimensional solid ODEs and
// finite contact jump. This reference never reads a prepared wall admittance.
fn circular_two_solid_reference(gas: &GasState, cp: C64, omega: f64) -> [C64; 5] {
    let (a, b, c) = (0.002, 0.0023, 0.0027);
    let (k1, c1, k2, c2, resistance, g) = (0.4, 100_000.0, 1.5, 200_000.0, 0.001, 30.0);
    let basis = regular_basis(
        C64::new(cp.im, cp.re).scale(omega * gas.density / gas.thermal_conductivity),
        a,
    );
    let integrate = |mut v: [C64; 2], lo: f64, hi: f64, k: f64, capacity: f64| {
        let steps = 4096;
        let h = (hi - lo) / steps as f64;
        let z = C64::new(0.0, omega * capacity / k);
        let rhs = |r: f64, t: [C64; 2]| [t[1], z * t[0] - t[1].scale(1.0 / r)];
        let add = |v: [C64; 2], p: [C64; 2], h: f64| [v[0] + p[0].scale(h), v[1] + p[1].scale(h)];
        for i in 0..steps {
            let r = lo + i as f64 * h;
            let p = rhs(r, v);
            let q = rhs(r + 0.5 * h, add(v, p, 0.5 * h));
            let s = rhs(r + 0.5 * h, add(v, q, 0.5 * h));
            let t = rhs(r + h, add(v, s, h));
            v = core::array::from_fn(|j| {
                v[j] + (p[j] + q[j].scale(2.0) + s[j].scale(2.0) + t[j]).scale(h / 6.0)
            });
        }
        v
    };
    let mut samples = [[C64::ZERO; 5]; 2];
    let mut residuals = [C64::ZERO; 2];
    for j in 0..2 {
        let (temperature, gradient, mean) = if j == 0 {
            (C64::ONE, C64::ZERO, C64::ONE)
        } else {
            (basis[0], basis[1], basis[2])
        };
        let first = integrate(
            [temperature, gradient.scale(gas.thermal_conductivity / k1)],
            a,
            b,
            k1,
            c1,
        );
        let jump_temperature = first[0] + first[1].scale(k1 * resistance);
        let second = integrate([jump_temperature, first[1].scale(k1 / k2)], b, c, k2, c2);
        samples[j] = [mean, temperature, first[0], jump_temperature, second[0]];
        residuals[j] = second[1].scale(k2) + second[0].scale(g);
    }
    let amplitude = -residuals[0] / residuals[1];
    core::array::from_fn(|i| samples[0][i] + amplitude * samples[1][i])
}

#[test]
fn circular_molecular_gas_and_two_real_source_layers_converge_to_independent_coupled_reference() {
    let material = material(2);
    let source_a = solid_source(0.4, 100_000.0);
    let source_b = solid_source(1.5, 200_000.0);
    let omega = 200.0;
    let cp = material.caloric_response(omega).unwrap().cp_j_kg_k;
    let expected = circular_two_solid_reference(material.gas(), cp, omega);
    let mut previous = f64::INFINITY;
    for cells in [24, 48, 96] {
        let section = CircularSection::new(0.002, cells).unwrap();
        let a = AnnularSection::new(0.002, 0.0023, cells).unwrap();
        let b = AnnularSection::new(0.0023, 0.0027, cells).unwrap();
        let wall = LayeredWall::prepare(
            &[
                Layer {
                    section: &a,
                    source: &source_a,
                },
                Layer {
                    section: &b,
                    source: &source_b,
                },
            ],
            &[0.001],
            293.15,
            omega,
            Face::Inner,
            Bath::Conductance {
                conductance_w_m2_k: 30.0,
            },
            solid_admission(),
        )
        .unwrap();
        let fields = section
            .prepare_molecular_thermal_wall_fields(
                &material,
                omega,
                CWall::Layered(&wall),
                admission(),
            )
            .unwrap();
        let error = (fields.mean_temperature_response() - expected[0]).abs();
        assert!(error < previous * 0.3, "N={cells}: {error:e}/{previous:e}");
        previous = error;
        if cells == 96 {
            near(fields.mean_temperature_response(), expected[0], 5e-4);
            let scale = cp.recip().scale(1.0 / material.gas().density);
            let boundary = fields.wall_observation(C64::ONE).unwrap();
            let contact = fields.thermal_contact_observation(0, C64::ONE).unwrap();
            for (actual, reference) in [
                (boundary.interface_temperature_k, expected[1]),
                (contact.inner_layer_temperature_k, expected[2]),
                (contact.outer_layer_temperature_k, expected[3]),
                (boundary.bath_face_temperature_k, expected[4]),
            ] {
                near(actual, reference.conj() * scale, 1e-3);
            }
        }
    }
}

#[test]
fn local_disk_volume_quadrature_recovers_flow_and_all_positive_gas_solid_contact_work() {
    let material = material(2);
    let source = solid_source(0.8, 120_000.0);
    let source_b = solid_source(1.5, 210_000.0);
    let section = CircularSection::new(0.002, 192).unwrap();
    let a = AnnularSection::new(0.002, 0.0023, 96).unwrap();
    let b = AnnularSection::new(0.0023, 0.0027, 96).unwrap();
    let (p, q) = (C64::new(3.0, -0.8), C64::new(2e-8, 3e-8));
    for omega in [20.0, 500.0, 3000.0] {
        let homogeneous = a
            .prepare_radial_thermal_wall(
                &source,
                293.15,
                omega,
                Face::Inner,
                Bath::Conductance {
                    conductance_w_m2_k: 40.0,
                },
                solid_admission(),
            )
            .unwrap();
        let layered = LayeredWall::prepare(
            &[
                Layer {
                    section: &a,
                    source: &source,
                },
                Layer {
                    section: &b,
                    source: &source_b,
                },
            ],
            &[0.001],
            293.15,
            omega,
            Face::Inner,
            Bath::Conductance {
                conductance_w_m2_k: 40.0,
            },
            solid_admission(),
        )
        .unwrap();
        for wall in [CWall::Solid(&homogeneous), CWall::Layered(&layered)] {
            let fields = section
                .prepare_molecular_thermal_wall_fields(&material, omega, wall, admission())
                .unwrap();
            let mut integral = [0.0; 3];
            let mut flow = C64::ZERO;
            for cell in 0..section.elements() {
                let dr =
                    section.node_radius_m(cell + 1).unwrap() - section.node_radius_m(cell).unwrap();
                for fraction in [0.5 - 0.5 / 3.0_f64.sqrt(), 0.5 + 0.5 / 3.0_f64.sqrt()] {
                    let value = fields.field_in_cell(cell, fraction, p, q).unwrap();
                    let weight = PI * value.radius_m * dr;
                    flow = flow + value.axial_velocity_m_s.scale(weight);
                    for (total, loss) in integral.iter_mut().zip([
                        value.viscous_dissipation_w_m3,
                        value.conduction_dissipation_w_m3,
                        value.molecular_dissipation_w_m3,
                    ]) {
                        *total += weight * loss;
                    }
                }
            }
            near(flow, q, 6e-13);
            let losses = fields.dissipation_per_axial_length(p, q).unwrap();
            for (actual, reference) in integral.into_iter().zip([
                losses.gas_viscous_w_m,
                losses.gas_conduction_w_m,
                losses.gas_molecular_w_m,
            ]) {
                relative(actual, reference, 5e-13);
            }
            let wave = fields.wave();
            let acoustic_work =
                0.5 * (wave.wavenumber * wave.characteristic_impedance).im * q.norm_sq()
                    + 0.5 * (wave.wavenumber / wave.characteristic_impedance).im * p.norm_sq();
            relative(losses.total_w_m, acoustic_work, 3e-10);
            let boundary = fields.wall_observation(p).unwrap();
            assert!(boundary.thermal_storage_j_m > 0.0);
            assert!(losses.solid_conduction_w_m > 0.0 && losses.bath_exchange_w_m > 0.0);
            let solid = match wall {
                CWall::Solid(_) => {
                    assert_eq!(losses.thermal_contact_w_m, 0.0);
                    fields.solid_field_in_cell(0, 0.0, p).unwrap()
                }
                CWall::Layered(_) => {
                    assert!(losses.thermal_contact_w_m > 0.0);
                    fields.layered_solid_field_in_cell(0, 0, 0.0, p).unwrap()
                }
                _ => unreachable!(),
            };
            near(solid.temperature_k, boundary.interface_temperature_k, 3e-15);
            let gas_wall = fields
                .field_in_cell(section.elements() - 1, 1.0, p, q)
                .unwrap();
            assert_eq!(gas_wall.gas_temperature_k, boundary.interface_temperature_k);
        }
    }
}

#[test]
fn axis_gradient_converges_to_regular_symmetry_without_a_pinned_center_velocity() {
    let gas = gas();
    let radius = 0.002;
    let q = C64::from_re(1e-8);
    let mut previous = f64::INFINITY;
    for cells in [16, 32, 64] {
        let section = CircularSection::new(radius, cells).unwrap();
        let fields = section
            .prepare_thermal_wall_fields(&gas, 1e-4, CWall::Boundary(Wall::Adiabatic), admission())
            .unwrap();
        let axis = fields.field_in_cell(0, 0.0, C64::ONE, q).unwrap();
        let gradient = axis.radial_velocity_gradient_per_s.abs();
        assert!(
            gradient < previous * 0.51,
            "N={cells}: {gradient:e}/{previous:e}"
        );
        previous = gradient;
        near(
            axis.gas_temperature_k,
            C64::from_re(1.0 / (gas.density * gas.specific_heat_cp)),
            3e-15,
        );
        assert_eq!(axis.radial_temperature_gradient_k_m, C64::ZERO);
        assert!(axis.axial_velocity_m_s.re > 1.9 * q.re / section.area_m2());
    }
}

#[test]
fn empty_molecular_fields_and_thermal_boundary_limits_preserve_the_original_material() {
    let material = material(0);
    let section = CircularSection::new(0.002, 192).unwrap();
    for boundary in [
        Wall::Isothermal,
        Wall::Adiabatic,
        Wall::Finite {
            conductance_w_m2_k: 25.0,
            areal_heat_capacity_j_m2_k: 0.1,
        },
    ] {
        let a = section
            .prepare_thermal_wall_fields(
                material.gas(),
                300.0,
                CWall::Boundary(boundary),
                admission(),
            )
            .unwrap();
        let b = section
            .prepare_molecular_thermal_wall_fields(
                &material,
                300.0,
                CWall::Boundary(boundary),
                admission(),
            )
            .unwrap();
        assert_eq!(a.wave().wavenumber, b.wave().wavenumber);
        assert_eq!(
            a.dissipation_per_axial_length(C64::ONE, C64::from_re(1e-8))
                .unwrap(),
            b.dissipation_per_axial_length(C64::ONE, C64::from_re(1e-8))
                .unwrap()
        );
        assert!(std::ptr::eq(b.source(), &material));
    }
    let insulated = section
        .prepare_thermal_wall_fields(
            material.gas(),
            300.0,
            CWall::Boundary(Wall::Adiabatic),
            admission(),
        )
        .unwrap();
    let finite_zero = section
        .prepare_thermal_wall_fields(
            material.gas(),
            300.0,
            CWall::Boundary(Wall::Finite {
                conductance_w_m2_k: 0.0,
                areal_heat_capacity_j_m2_k: 0.0,
            }),
            admission(),
        )
        .unwrap();
    assert_eq!(insulated.wave().wavenumber, finite_zero.wave().wavenumber);
    let loss = insulated
        .dissipation_per_axial_length(C64::ONE, C64::ZERO)
        .unwrap();
    assert_eq!(loss.total_w_m, 0.0);
    let cold = section
        .prepare_thermal_wall_fields(
            material.gas(),
            300.0,
            CWall::Boundary(Wall::Isothermal),
            admission(),
        )
        .unwrap();
    let almost_cold = section
        .prepare_thermal_wall_fields(
            material.gas(),
            300.0,
            CWall::Boundary(Wall::Finite {
                conductance_w_m2_k: 1e9,
                areal_heat_capacity_j_m2_k: 0.0,
            }),
            admission(),
        )
        .unwrap();
    near(almost_cold.wave().wavenumber, cold.wave().wavenumber, 2e-7);
}

#[test]
fn axis_geometry_material_wall_matching_resolution_and_cancellation_refuse_transactionally() {
    for (radius, cells) in [
        (0.0, 64),
        (-0.001, 64),
        (f64::NAN, 64),
        (0.002, 1),
        (0.002, 4097),
    ] {
        assert!(CircularSection::new(radius, cells).is_err());
    }
    assert!(matches!(
        CircularSection::new_controlled(0.002, 64, || true),
        Err(DuctError::Cancelled { .. })
    ));
    let gas = gas();
    let section = CircularSection::new(0.002, 64).unwrap();
    for config in [
        AnnularAdmission {
            maximum_payload_bytes: 0,
            ..admission()
        },
        AnnularAdmission {
            maximum_cell_diffusion_phase: 1e-5,
            ..admission()
        },
        AnnularAdmission {
            maximum_transverse_phase: 1e-6,
            ..admission()
        },
    ] {
        assert!(
            section
                .prepare_thermal_wall_fields(&gas, 300.0, CWall::Boundary(Wall::Isothermal), config)
                .is_err()
        );
    }
    let source = solid_source(0.8, 120_000.0);
    let shell = AnnularSection::new(0.002, 0.0024, 64).unwrap();
    for (t, omega, face) in [
        (294.0, 300.0, Face::Inner),
        (293.15, 301.0, Face::Inner),
        (293.15, 300.0, Face::Outer),
    ] {
        let wall = shell
            .prepare_radial_thermal_wall(
                &source,
                t,
                omega,
                face,
                Bath::Adiabatic,
                solid_admission(),
            )
            .unwrap();
        assert!(
            section
                .prepare_thermal_wall_fields(&gas, 300.0, CWall::Solid(&wall), admission())
                .is_err()
        );
    }
    let gap = AnnularSection::new(0.0021, 0.0024, 64).unwrap();
    let detached = gap
        .prepare_radial_thermal_wall(
            &source,
            293.15,
            300.0,
            Face::Inner,
            Bath::Adiabatic,
            solid_admission(),
        )
        .unwrap();
    assert!(
        section
            .prepare_thermal_wall_fields(&gas, 300.0, CWall::Solid(&detached), admission())
            .is_err()
    );
    let mut incoherent = gas;
    incoherent.density *= 2.0;
    assert!(
        section
            .prepare_thermal_wall_fields(
                &incoherent,
                300.0,
                CWall::Boundary(Wall::Isothermal),
                admission()
            )
            .is_err()
    );
    assert!(matches!(
        section.prepare_thermal_wall_fields_controlled(
            &gas,
            300.0,
            CWall::Boundary(Wall::Isothermal),
            admission(),
            || true
        ),
        Err(RadialThermalWallError::Duct(DuctError::Cancelled { .. }))
    ));
    let mut calls = 0;
    section
        .prepare_thermal_wall_fields_controlled(
            &gas,
            300.0,
            CWall::Boundary(Wall::Isothermal),
            admission(),
            || {
                calls += 1;
                false
            },
        )
        .unwrap();
    let mut retry = 0;
    assert!(matches!(
        section.prepare_thermal_wall_fields_controlled(
            &gas,
            300.0,
            CWall::Boundary(Wall::Isothermal),
            admission(),
            || {
                retry += 1;
                retry == calls
            }
        ),
        Err(RadialThermalWallError::Fluid(_))
    ));
}
