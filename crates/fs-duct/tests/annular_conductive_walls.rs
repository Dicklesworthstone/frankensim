//! G1/G3: resolved source-bound solid walls, dimensional references and work.
//! Material values below are explicitly synthetic, not empirical solid data.
use fs_duct::DuctError;
use fs_duct::fluid::annular::thermal_wall::AnnularThermalWall as Wall;
use fs_duct::fluid::annular::thermal_wall::conduction::{
    AnnularConjugateWall as CoupledWall, RadialThermalWallAdmission, RadialThermalWallError,
    RadialWallBath as Bath, RadialWallFace as Face,
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
use std::f64::consts::{PI, TAU};

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

#[test]
fn quasistatic_cylindrical_resistance_and_each_contact_orientation_converge() {
    let source = solid_source(2.0, 1000.0);
    let (ri, ro, g) = (0.002_f64, 0.004_f64, 60.0);
    for face in [Face::Inner, Face::Outer] {
        let (rc, rb) = if face == Face::Inner {
            (ri, ro)
        } else {
            (ro, ri)
        };
        for bath in [
            Bath::Isothermal,
            Bath::Conductance {
                conductance_w_m2_k: g,
            },
        ] {
            let resistance = (ro / ri).ln() / (TAU * 2.0)
                + if bath == Bath::Isothermal {
                    0.0
                } else {
                    1.0 / (TAU * rb * g)
                };
            let expected = 1.0 / (TAU * rc * resistance);
            let mut previous = f64::INFINITY;
            for n in [32, 64, 128] {
                let shell = AnnularSection::new(ri, ro, n).unwrap();
                let wall = shell
                    .prepare_radial_thermal_wall(
                        &source,
                        293.15,
                        1e-5,
                        face,
                        bath,
                        solid_admission(),
                    )
                    .unwrap();
                let y = wall.interface_admittance_w_m2_k();
                let error = (y.re - expected).abs();
                assert!(
                    error < 0.28 * previous,
                    "{face:?} {bath:?} n={n}, {error:e}/{previous:e}"
                );
                previous = error;
                assert!(y.im < 0.0);
                let response = wall.observation(C64::ONE).unwrap();
                near(
                    response.fluid_to_wall_heat_flux_w_m2.scale(rc),
                    response.wall_to_bath_heat_flux_w_m2.scale(rb),
                    2e-6,
                );
                if n == 128 {
                    relative(y.re, expected, 3e-6);
                }
            }
        }
    }
}

// Independent dimensional coupled ODE shooting. The gas contains forcing
// z*(h-1); the solid contains z*h. Temperature and k*dT/dr are continuous.
// No production stiffness, mass, response admittance or linear solve is reused.
fn coupled_ode(gas: &GasState, cp: C64, omega: f64, k_s: f64, c_s: f64, g: f64) -> [C64; 3] {
    let (a, b, c) = (0.002_f64, 0.003_f64, 0.0036_f64);
    let integrate = |mut v: [C64; 3], lo: f64, hi: f64, z: C64, force: bool, area: bool| {
        let steps = 4096;
        let h = (hi - lo) / steps as f64;
        let rhs = |r: f64, u: [C64; 3]| {
            [
                u[1],
                z * (u[0] - if force { C64::ONE } else { C64::ZERO }) - u[1].scale(1.0 / r),
                if area { u[0].scale(TAU * r) } else { C64::ZERO },
            ]
        };
        let step =
            |u: [C64; 3], k: [C64; 3], dt: f64| core::array::from_fn(|i| u[i] + k[i].scale(dt));
        for i in 0..steps {
            let r = lo + i as f64 * h;
            let k1 = rhs(r, v);
            let k2 = rhs(r + h * 0.5, step(v, k1, h * 0.5));
            let k3 = rhs(r + h * 0.5, step(v, k2, h * 0.5));
            let k4 = rhs(r + h, step(v, k3, h));
            v = core::array::from_fn(|j| {
                v[j] + (k1[j] + k2[j].scale(2.0) + k3[j].scale(2.0) + k4[j]).scale(h / 6.0)
            });
        }
        v
    };
    let mut end = [[C64::ZERO; 3]; 2];
    let mut interface = [C64::ZERO; 2];
    for j in 0..2 {
        let v = integrate(
            [
                if j == 0 { C64::ZERO } else { C64::ONE },
                C64::ZERO,
                C64::ZERO,
            ],
            a,
            b,
            C64::new(cp.im, cp.re).scale(omega * gas.density / gas.thermal_conductivity),
            j == 0,
            true,
        );
        interface[j] = v[0];
        end[j] = integrate(
            [v[0], v[1].scale(gas.thermal_conductivity / k_s), v[2]],
            b,
            c,
            C64::new(0.0, omega * c_s / k_s),
            false,
            false,
        );
    }
    let residual = |v: [C64; 3]| v[1].scale(k_s) + v[0].scale(g);
    let scale = -residual(end[0]) / residual(end[1]);
    [
        (end[0][2] + scale * end[1][2]).scale(1.0 / (PI * (b * b - a * a))),
        interface[0] + scale * interface[1],
        end[0][0] + scale * end[1][0],
    ]
}

#[test]
fn molecular_gas_and_resolved_solid_converge_to_independent_coupled_radial_equations() {
    let material = material(2);
    let source = solid_source(2.0, 100_000.0);
    let omega = 400.0;
    let cp = material.caloric_response(omega).unwrap().cp_j_kg_k;
    let expected = coupled_ode(material.gas(), cp, omega, 2.0, 100_000.0, 30.0);
    let mut previous = f64::INFINITY;
    for n in [16, 32, 64] {
        let gas_section = AnnularSection::new(0.002, 0.003, n).unwrap();
        let shell = AnnularSection::new(0.003, 0.0036, n).unwrap();
        let solid = shell
            .prepare_radial_thermal_wall(
                &source,
                293.15,
                omega,
                Face::Inner,
                Bath::Conductance {
                    conductance_w_m2_k: 30.0,
                },
                solid_admission(),
            )
            .unwrap();
        let fields = gas_section
            .prepare_molecular_conjugate_wall_fields(
                &material,
                omega,
                [
                    CoupledWall::Boundary(Wall::Adiabatic),
                    CoupledWall::Solid(&solid),
                ],
                admission(),
            )
            .unwrap();
        let error = (fields.mean_temperature_response() - expected[0]).abs();
        assert!(error < previous * 0.3, "n={n}: {error:e}/{previous:e}");
        previous = error;
        if n == 64 {
            near(fields.mean_temperature_response(), expected[0], 4e-4);
            let wall = fields.wall_observation(1, C64::ONE).unwrap();
            let scale = cp.recip().scale(1.0 / material.gas().density);
            near(
                wall.interface_temperature_k,
                expected[1].conj() * scale,
                5e-4,
            );
            near(
                wall.bath_face_temperature_k,
                expected[2].conj() * scale,
                8e-4,
            );
            assert!(std::ptr::eq(fields.source(), &material));
        }
    }
}

#[test]
fn physical_solid_quadrature_and_complex_heat_flux_close_for_both_orientations() {
    let source = solid_source(0.6, 200_000.0);
    let section = AnnularSection::new(0.002, 0.0027, 192).unwrap();
    let (omega, temperature) = (180.0, C64::new(0.04, -0.07));
    for face in [Face::Inner, Face::Outer] {
        for bath in [
            Bath::Adiabatic,
            Bath::Isothermal,
            Bath::Conductance {
                conductance_w_m2_k: 45.0,
            },
        ] {
            let wall = section
                .prepare_radial_thermal_wall(&source, 293.15, omega, face, bath, solid_admission())
                .unwrap();
            let mut integrals = [0.0; 2];
            let mut temperature_integral = C64::ZERO;
            for cell in 0..section.elements() {
                let dr = section.radius_m(cell + 1).unwrap() - section.radius_m(cell).unwrap();
                for fraction in [0.5 - 0.5 / 3.0_f64.sqrt(), 0.5 + 0.5 / 3.0_f64.sqrt()] {
                    let value = wall.field_in_cell(cell, fraction, temperature).unwrap();
                    let weight = PI * value.radius_m * dr;
                    integrals[0] += weight * value.conduction_dissipation_w_m3;
                    integrals[1] += weight * value.thermal_storage_j_m3;
                    temperature_integral = temperature_integral + value.temperature_k.scale(weight);
                    near(
                        value.radial_heat_flux_w_m2,
                        value.radial_temperature_gradient_k_m.scale(-0.6),
                        2e-15,
                    );
                }
            }
            let observed = wall.observation(temperature).unwrap();
            relative(integrals[0], observed.solid_conduction_w_m, 3e-13);
            relative(integrals[1], observed.thermal_storage_j_m, 3e-13);
            near(
                temperature_integral,
                observed.temperature_integral_k_m2,
                3e-13,
            );
            let input = observed
                .fluid_to_wall_heat_flux_w_m2
                .scale(TAU * observed.interface_radius_m);
            let output = observed
                .wall_to_bath_heat_flux_w_m2
                .scale(TAU * observed.bath_radius_m);
            near(
                input - output,
                C64::new(0.0, -omega * 200_000.0) * temperature_integral,
                2e-10,
            );
            relative(
                (temperature.conj() * input).re / (2.0 * 293.15),
                observed.solid_conduction_w_m + observed.bath_exchange_w_m,
                2e-13,
            );
            relative(
                -(temperature.conj() * input).im / (4.0 * 293.15 * omega),
                observed.thermal_storage_j_m,
                3e-13,
            );
            if bath == Bath::Adiabatic {
                assert_eq!(observed.wall_to_bath_heat_flux_w_m2, C64::ZERO);
                assert_eq!(observed.bath_exchange_w_m, 0.0);
            }
            assert!(std::ptr::eq(wall.source(), &source));
            assert!(std::ptr::eq(wall.section(), &section));
            assert_eq!(wall.sample().path_identity(), source.identity());
            assert_eq!(wall.sample().conductivity_w_m_k(), 0.6);
        }
    }
}

#[test]
fn two_solid_walls_split_gas_solid_and_bath_work_in_original_acoustic_balance() {
    let material = material(2);
    let source_a = solid_source(0.8, 120_000.0);
    let source_b = solid_source(2.0, 200_000.0);
    let gas_section = AnnularSection::new(0.002, 0.003, 192).unwrap();
    let shell_a = AnnularSection::new(0.0015, 0.002, 192).unwrap();
    let shell_b = AnnularSection::new(0.003, 0.0037, 192).unwrap();
    let (p, q) = (C64::new(3.0, -0.7), C64::new(2e-8, 1e-8));
    for omega in [20.0, 400.0, 3000.0] {
        let a = shell_a
            .prepare_radial_thermal_wall(
                &source_a,
                293.15,
                omega,
                Face::Outer,
                Bath::Conductance {
                    conductance_w_m2_k: 70.0,
                },
                solid_admission(),
            )
            .unwrap();
        let b = shell_b
            .prepare_radial_thermal_wall(
                &source_b,
                293.15,
                omega,
                Face::Inner,
                Bath::Adiabatic,
                solid_admission(),
            )
            .unwrap();
        let fields = gas_section
            .prepare_molecular_conjugate_wall_fields(
                &material,
                omega,
                [CoupledWall::Solid(&a), CoupledWall::Solid(&b)],
                admission(),
            )
            .unwrap();
        let losses = fields.dissipation_per_axial_length(p, q).unwrap();
        let wave = fields.wave();
        let work = 0.5 * (wave.wavenumber * wave.characteristic_impedance).im * q.norm_sq()
            + 0.5 * (wave.wavenumber / wave.characteristic_impedance).im * p.norm_sq();
        relative(losses.total_w_m, work, 3e-10);
        assert!(
            losses.gas_viscous_w_m > 0.0
                && losses.gas_conduction_w_m > 0.0
                && losses.gas_molecular_w_m > 0.0
        );
        assert!(losses.solid_conduction_w_m.into_iter().all(|x| x > 0.0));
        assert!(losses.bath_exchange_w_m[0] > 0.0);
        assert_eq!(losses.bath_exchange_w_m[1], 0.0);
        for (side, solid) in [(0, &a), (1, &b)] {
            let observed = fields.wall_observation(side, p).unwrap();
            let expected = solid.observation(observed.interface_temperature_k).unwrap();
            assert_eq!(observed.solid_conduction_w_m, expected.solid_conduction_w_m);
            assert_eq!(observed.bath_exchange_w_m, expected.bath_exchange_w_m);
            assert_eq!(observed.thermal_storage_j_m, expected.thermal_storage_j_m);
            let gas_field = fields
                .field_in_cell(if side == 0 { 0 } else { 191 }, side as f64, p, q)
                .unwrap();
            let solid_field = fields
                .solid_field_in_cell(side, if side == 0 { 191 } else { 0 }, (1 - side) as f64, p)
                .unwrap();
            near(
                gas_field.gas_temperature_k,
                solid_field.temperature_k,
                3e-15,
            );
        }
    }
}

#[test]
fn thermally_thin_wall_recovers_area_correct_lumped_capacity_and_external_conductance() {
    let source = solid_source(100.0, 10_000.0);
    let gas = gas();
    let gas_section = AnnularSection::new(0.002, 0.003, 96).unwrap();
    let shell = AnnularSection::new(0.003, 0.0031, 32).unwrap();
    let solid = shell
        .prepare_radial_thermal_wall(
            &source,
            293.15,
            20.0,
            Face::Inner,
            Bath::Conductance {
                conductance_w_m2_k: 25.0,
            },
            solid_admission(),
        )
        .unwrap();
    let g = 25.0 * 0.0031 / 0.003;
    let c = 10_000.0 * shell.area_m2() / (TAU * 0.003);
    near(
        solid.interface_admittance_w_m2_k(),
        C64::new(g, -20.0 * c),
        5e-5,
    );
    let old = gas_section
        .prepare_thermal_wall_fields(
            &gas,
            20.0,
            [
                Wall::Adiabatic,
                Wall::Finite {
                    conductance_w_m2_k: g,
                    areal_heat_capacity_j_m2_k: c,
                },
            ],
            admission(),
        )
        .unwrap();
    let resolved = gas_section
        .prepare_conjugate_wall_fields(
            &gas,
            20.0,
            [
                CoupledWall::Boundary(Wall::Adiabatic),
                CoupledWall::Solid(&solid),
            ],
            admission(),
        )
        .unwrap();
    near(resolved.wave().wavenumber, old.wave().wavenumber, 1e-5);
    let observation = solid.observation(C64::ONE).unwrap();
    relative(
        observation.thermal_storage_j_m,
        10_000.0 * shell.area_m2() / (4.0 * 293.15),
        6e-5,
    );
    assert!(observation.solid_conduction_w_m > 0.0);
}

#[test]
fn thick_large_radius_wall_reaches_the_independent_thermal_effusivity_limit() {
    let source = solid_source(0.2, 10_000.0);
    let shell = AnnularSection::new(0.1, 0.105, 512).unwrap();
    let wall = shell
        .prepare_radial_thermal_wall(
            &source,
            293.15,
            200.0,
            Face::Inner,
            Bath::Adiabatic,
            solid_admission(),
        )
        .unwrap();
    // Planar semi-infinite solution exp[-sqrt(-i*omega/alpha)*x]:
    // Y = sqrt(k*C_vol)*sqrt(-i*omega). Curvature is O(delta/r).
    let root = (0.2_f64 * 10_000.0 * 200.0 / 2.0).sqrt();
    near(
        wall.interface_admittance_w_m2_k(),
        C64::new(root, -root),
        0.002,
    );
    let observed = wall.observation(C64::ONE).unwrap();
    assert!(observed.bath_temperature_k.abs() < 4e-5);
    assert_eq!(observed.bath_exchange_w_m, 0.0);
}

#[test]
fn empty_molecular_owner_matches_equilibrium_and_boundary_only_calls_preserve_original_owner() {
    let material = material(0);
    let source = solid_source(0.5, 80_000.0);
    let shell = AnnularSection::new(0.003, 0.0035, 96).unwrap();
    let section = AnnularSection::new(0.002, 0.003, 96).unwrap();
    let solid = shell
        .prepare_radial_thermal_wall(
            &source,
            293.15,
            300.0,
            Face::Inner,
            Bath::Isothermal,
            solid_admission(),
        )
        .unwrap();
    let walls = [
        CoupledWall::Boundary(Wall::Adiabatic),
        CoupledWall::Solid(&solid),
    ];
    let a = section
        .prepare_conjugate_wall_fields(material.gas(), 300.0, walls, admission())
        .unwrap();
    let b = section
        .prepare_molecular_conjugate_wall_fields(&material, 300.0, walls, admission())
        .unwrap();
    assert_eq!(a.wave().wavenumber, b.wave().wavenumber);
    assert_eq!(
        a.dissipation_per_axial_length(C64::ONE, C64::ZERO).unwrap(),
        b.dissipation_per_axial_length(C64::ONE, C64::ZERO).unwrap()
    );
    let walls = [
        Wall::Finite {
            conductance_w_m2_k: 30.0,
            areal_heat_capacity_j_m2_k: 0.1,
        },
        Wall::Isothermal,
    ];
    let old = section
        .prepare_thermal_wall_fields(material.gas(), 300.0, walls, admission())
        .unwrap();
    let same = section
        .prepare_conjugate_wall_fields(
            material.gas(),
            300.0,
            walls.map(CoupledWall::Boundary),
            admission(),
        )
        .unwrap();
    assert_eq!(same.wave().wavenumber, old.wave().wavenumber);
    for side in 0..2 {
        let a = same.wall_observation(side, C64::ONE).unwrap();
        let b = old.wall_observation(side, C64::ONE).unwrap();
        assert_eq!(a.bath_exchange_w_m, b.bath_dissipation_w_m);
        assert_eq!(a.thermal_storage_j_m, b.thermal_storage_j_m);
        assert_eq!(a.solid_conduction_w_m, 0.0);
    }
}

#[test]
fn incompatible_material_geometry_frequency_and_unbounded_work_refuse() {
    let source = solid_source(0.5, 80_000.0);
    let shell = AnnularSection::new(0.003, 0.0035, 96).unwrap();
    let section = AnnularSection::new(0.002, 0.003, 96).unwrap();
    for (t, w) in [
        (500.1, 300.0),
        (f64::NAN, 300.0),
        (293.15, 0.0),
        (293.15, f64::INFINITY),
    ] {
        assert!(
            shell
                .prepare_radial_thermal_wall(
                    &source,
                    t,
                    w,
                    Face::Inner,
                    Bath::Adiabatic,
                    solid_admission()
                )
                .is_err()
        );
    }
    for config in [
        RadialThermalWallAdmission {
            maximum_payload_bytes: 0,
            ..solid_admission()
        },
        RadialThermalWallAdmission {
            maximum_cell_diffusion_phase: 1e-5,
            ..solid_admission()
        },
    ] {
        assert!(
            shell
                .prepare_radial_thermal_wall(
                    &source,
                    293.15,
                    300.0,
                    Face::Inner,
                    Bath::Adiabatic,
                    config
                )
                .is_err()
        );
    }
    assert!(
        shell
            .prepare_radial_thermal_wall(
                &source,
                293.15,
                300.0,
                Face::Inner,
                Bath::Conductance {
                    conductance_w_m2_k: -1.0
                },
                solid_admission()
            )
            .is_err()
    );
    assert!(matches!(
        shell.prepare_radial_thermal_wall_controlled(
            &source,
            293.15,
            300.0,
            Face::Inner,
            Bath::Adiabatic,
            solid_admission(),
            || true
        ),
        Err(RadialThermalWallError::Duct(DuctError::Cancelled { .. }))
    ));
    let mut calls = 0;
    shell
        .prepare_radial_thermal_wall_controlled(
            &source,
            293.15,
            300.0,
            Face::Inner,
            Bath::Adiabatic,
            solid_admission(),
            || {
                calls += 1;
                false
            },
        )
        .unwrap();
    let mut retry = 0;
    assert!(matches!(
        shell.prepare_radial_thermal_wall_controlled(
            &source,
            293.15,
            300.0,
            Face::Inner,
            Bath::Adiabatic,
            solid_admission(),
            || {
                retry += 1;
                retry == calls
            }
        ),
        Err(RadialThermalWallError::Duct(DuctError::Cancelled { .. }))
    ));
    let gas = gas();
    for (temperature, omega, face) in [
        (293.15, 301.0, Face::Inner),
        (294.0, 300.0, Face::Inner),
        (293.15, 300.0, Face::Outer),
    ] {
        let solid = shell
            .prepare_radial_thermal_wall(
                &source,
                temperature,
                omega,
                face,
                Bath::Adiabatic,
                solid_admission(),
            )
            .unwrap();
        assert!(
            section
                .prepare_conjugate_wall_fields(
                    &gas,
                    300.0,
                    [
                        CoupledWall::Boundary(Wall::Adiabatic),
                        CoupledWall::Solid(&solid)
                    ],
                    admission()
                )
                .is_err()
        );
    }
    let separated = AnnularSection::new(0.0031, 0.0035, 96).unwrap();
    let solid = separated
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
            .prepare_conjugate_wall_fields(
                &gas,
                300.0,
                [
                    CoupledWall::Boundary(Wall::Adiabatic),
                    CoupledWall::Solid(&solid)
                ],
                admission()
            )
            .is_err()
    );
    assert!(solid.field_in_cell(96, 0.5, C64::ONE).is_err());
    assert!(solid.field_in_cell(0, -0.1, C64::ONE).is_err());
    assert!(solid.observation(C64::new(f64::INFINITY, 0.0)).is_err());
}
