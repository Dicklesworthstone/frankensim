//! G1/G3: actual radial material layers, thermal contact jumps and acoustic work.
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

use fs_duct::fluid::annular::thermal_wall::conduction::layered::{
    PreparedLayeredRadialThermalWall as LayeredWall, RadialThermalLayer as Layer,
};

#[test]
fn one_layer_and_perfectly_joined_identical_material_recover_original_resolved_wall() {
    let source = solid_source(0.7, 140_000.0);
    let (a, b) = (0.002_f64, 0.004_f64);
    let mid = (a * b).sqrt();
    let all = AnnularSection::new(a, b, 128).unwrap();
    let first = AnnularSection::new(a, mid, 64).unwrap();
    let second = AnnularSection::new(mid, b, 64).unwrap();
    let one = [Layer {
        section: &all,
        source: &source,
    }];
    let split = [
        Layer {
            section: &first,
            source: &source,
        },
        Layer {
            section: &second,
            source: &source,
        },
    ];
    for face in [Face::Inner, Face::Outer] {
        for bath in [
            Bath::Adiabatic,
            Bath::Isothermal,
            Bath::Conductance {
                conductance_w_m2_k: 40.0,
            },
        ] {
            let original = all
                .prepare_radial_thermal_wall(&source, 293.15, 100.0, face, bath, solid_admission())
                .unwrap();
            for (layers, contacts) in [(&one[..], &[][..]), (&split[..], &[0.0][..])] {
                let wall = LayeredWall::prepare(
                    layers,
                    contacts,
                    293.15,
                    100.0,
                    face,
                    bath,
                    solid_admission(),
                )
                .unwrap();
                near(
                    wall.interface_admittance_w_m2_k(),
                    original.interface_admittance_w_m2_k(),
                    2e-11,
                );
                let t = C64::new(0.03, -0.04);
                let actual = wall.observation(t).unwrap();
                let expected = original.observation(t).unwrap();
                relative(
                    actual.solid_conduction_w_m,
                    expected.solid_conduction_w_m,
                    3e-11,
                );
                relative(
                    actual.thermal_storage_j_m,
                    expected.thermal_storage_j_m,
                    3e-11,
                );
                assert_eq!(actual.thermal_contact_w_m, 0.0);
                if contacts.len() == 1 {
                    let contact = wall.contact_observation(0, t).unwrap();
                    assert_eq!(
                        contact.inner_layer_temperature_k,
                        contact.outer_layer_temperature_k
                    );
                    assert_eq!(contact.contact_dissipation_w_m, 0.0);
                    assert!(std::ptr::eq(wall.layer_source(0).unwrap(), &source));
                    assert!(std::ptr::eq(wall.layer_section(1).unwrap(), &second));
                }
            }
        }
    }
}

#[test]
fn multilayer_cylindrical_resistance_contact_jump_and_orientation_match_static_reference() {
    let a = solid_source(0.4, 8_000.0);
    let b = solid_source(2.0, 20_000.0);
    let radii = [0.002_f64, 0.003_f64, 0.004_f64];
    let resistance = 0.003;
    let g = 20.0;
    for face in [Face::Inner, Face::Outer] {
        let (rc, rb, sign) = if face == Face::Inner {
            (radii[0], radii[2], 1.0)
        } else {
            (radii[2], radii[0], -1.0)
        };
        let total = (radii[1] / radii[0]).ln() / (TAU * 0.4)
            + (radii[2] / radii[1]).ln() / (TAU * 2.0)
            + resistance / (TAU * radii[1])
            + 1.0 / (TAU * rb * g);
        let expected = 1.0 / (TAU * rc * total);
        let mut previous = f64::INFINITY;
        for cells in [16, 32, 64] {
            let left = AnnularSection::new(radii[0], radii[1], cells).unwrap();
            let right = AnnularSection::new(radii[1], radii[2], cells).unwrap();
            let wall = LayeredWall::prepare(
                &[
                    Layer {
                        section: &left,
                        source: &a,
                    },
                    Layer {
                        section: &right,
                        source: &b,
                    },
                ],
                &[resistance],
                293.15,
                1e-5,
                face,
                Bath::Conductance {
                    conductance_w_m2_k: g,
                },
                solid_admission(),
            )
            .unwrap();
            let error = (wall.interface_admittance_w_m2_k().re - expected).abs();
            assert!(
                error < previous * 0.3,
                "{face:?}, n={cells}: {error:e}/{previous:e}"
            );
            previous = error;
            if cells == 64 {
                relative(wall.interface_admittance_w_m2_k().re, expected, 8e-6);
                let contact = wall.contact_observation(0, C64::ONE).unwrap();
                near(
                    contact.outward_heat_flux_w_m2,
                    C64::from_re(sign / (TAU * radii[1] * total)),
                    2e-5,
                );
                near(
                    contact.inner_layer_temperature_k - contact.outer_layer_temperature_k,
                    contact.outward_heat_flux_w_m2.scale(resistance),
                    2e-14,
                );
            }
        }
    }
}

// Independent dimensional ODE through gas, two solids and a finite contact.
// No production basis, radial matrix, thermal admittance or face solve is reused.
fn coupled_layered_ode(gas: &GasState, cp: C64, omega: f64, resistance: f64) -> [C64; 5] {
    let (a, b, c, d) = (0.002_f64, 0.003_f64, 0.0034_f64, 0.0038_f64);
    let (k1, c1, k2, c2, g) = (0.3, 80_000.0, 2.0, 180_000.0, 35.0);
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
        let add =
            |u: [C64; 3], k: [C64; 3], dt: f64| core::array::from_fn(|i| u[i] + k[i].scale(dt));
        for step in 0..steps {
            let r = lo + step as f64 * h;
            let p = rhs(r, v);
            let q = rhs(r + h * 0.5, add(v, p, h * 0.5));
            let s = rhs(r + h * 0.5, add(v, q, h * 0.5));
            let t = rhs(r + h, add(v, s, h));
            v = core::array::from_fn(|i| {
                v[i] + (p[i] + q[i].scale(2.0) + s[i].scale(2.0) + t[i]).scale(h / 6.0)
            });
        }
        v
    };
    let mut results = [[C64::ZERO; 5]; 2];
    let mut residuals = [C64::ZERO; 2];
    for j in 0..2 {
        let gas_end = integrate(
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
        let first = integrate(
            [
                gas_end[0],
                gas_end[1].scale(gas.thermal_conductivity / k1),
                gas_end[2],
            ],
            b,
            c,
            C64::new(0.0, omega * c1 / k1),
            false,
            false,
        );
        let outer_temperature = first[0] + first[1].scale(k1 * resistance);
        let second = integrate(
            [outer_temperature, first[1].scale(k1 / k2), first[2]],
            c,
            d,
            C64::new(0.0, omega * c2 / k2),
            false,
            false,
        );
        results[j] = [
            gas_end[2].scale(1.0 / (PI * (b * b - a * a))),
            gas_end[0],
            first[0],
            outer_temperature,
            second[0],
        ];
        residuals[j] = second[1].scale(k2) + second[0].scale(g);
    }
    let scale = -residuals[0] / residuals[1];
    core::array::from_fn(|i| results[0][i] + scale * results[1][i])
}

#[test]
fn molecular_gas_and_two_source_layers_with_a_contact_converge_to_independent_ode() {
    let material = material(2);
    let first_source = solid_source(0.3, 80_000.0);
    let second_source = solid_source(2.0, 180_000.0);
    let omega = 200.0;
    let resistance = 0.0007;
    let cp = material.caloric_response(omega).unwrap().cp_j_kg_k;
    let expected = coupled_layered_ode(material.gas(), cp, omega, resistance);
    let mut previous = f64::INFINITY;
    for n in [16, 32, 64] {
        let gas_section = AnnularSection::new(0.002, 0.003, n).unwrap();
        let first = AnnularSection::new(0.003, 0.0034, n).unwrap();
        let second = AnnularSection::new(0.0034, 0.0038, n).unwrap();
        let wall = LayeredWall::prepare(
            &[
                Layer {
                    section: &first,
                    source: &first_source,
                },
                Layer {
                    section: &second,
                    source: &second_source,
                },
            ],
            &[resistance],
            293.15,
            omega,
            Face::Inner,
            Bath::Conductance {
                conductance_w_m2_k: 35.0,
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
                    CoupledWall::Layered(&wall),
                ],
                admission(),
            )
            .unwrap();
        let error = (fields.mean_temperature_response() - expected[0]).abs();
        assert!(error < previous * 0.3, "n={n}: {error:e}/{previous:e}");
        previous = error;
        if n == 64 {
            near(fields.mean_temperature_response(), expected[0], 4e-4);
            let scale = cp.recip().scale(1.0 / material.gas().density);
            let result = fields.wall_observation(1, C64::ONE).unwrap();
            let contact = fields.thermal_contact_observation(1, 0, C64::ONE).unwrap();
            for (actual, reference) in [
                (result.interface_temperature_k, expected[1]),
                (contact.inner_layer_temperature_k, expected[2]),
                (contact.outer_layer_temperature_k, expected[3]),
                (result.bath_face_temperature_k, expected[4]),
            ] {
                near(actual, reference.conj() * scale, 9e-4);
            }
            assert!(result.thermal_contact_w_m > 0.0);
        }
    }
}

#[test]
fn coherent_layer_quadrature_contact_loss_storage_and_original_heat_flux_close() {
    let sources = [
        solid_source(0.5, 100_000.0),
        solid_source(1.2, 240_000.0),
        solid_source(0.2, 70_000.0),
    ];
    let sections = [
        AnnularSection::new(0.002, 0.0023, 128).unwrap(),
        AnnularSection::new(0.0023, 0.0027, 128).unwrap(),
        AnnularSection::new(0.0027, 0.003, 128).unwrap(),
    ];
    let layers: Vec<_> = sections
        .iter()
        .zip(&sources)
        .map(|(section, source)| Layer { section, source })
        .collect();
    let temperature = C64::new(0.06, -0.03);
    let omega = 240.0;
    for face in [Face::Inner, Face::Outer] {
        let wall = LayeredWall::prepare(
            &layers,
            &[0.0, 0.001],
            293.15,
            omega,
            face,
            Bath::Conductance {
                conductance_w_m2_k: 50.0,
            },
            solid_admission(),
        )
        .unwrap();
        let mut losses = 0.0;
        let mut storage = 0.0;
        let mut energy = C64::ZERO;
        for (layer, section) in sections.iter().enumerate() {
            let sample = wall.layer_sample(layer).unwrap();
            let mut integrals = [0.0; 2];
            let mut local_energy = C64::ZERO;
            for cell in 0..section.elements() {
                let dr = section.radius_m(cell + 1).unwrap() - section.radius_m(cell).unwrap();
                for x in [0.5 - 0.5 / 3.0_f64.sqrt(), 0.5 + 0.5 / 3.0_f64.sqrt()] {
                    let field = wall.field_in_cell(layer, cell, x, temperature).unwrap();
                    let weight = PI * field.radius_m * dr;
                    integrals[0] += weight * field.conduction_dissipation_w_m3;
                    integrals[1] += weight * field.thermal_storage_j_m3;
                    local_energy = local_energy
                        + field
                            .temperature_k
                            .scale(weight * sample.volumetric_heat_capacity_j_m3_k());
                }
            }
            let result = wall.layer_observation(layer, temperature).unwrap();
            relative(integrals[0], result.conduction_dissipation_w_m, 4e-13);
            relative(integrals[1], result.thermal_storage_j_m, 4e-13);
            near(local_energy, result.oscillatory_thermal_energy_j_m, 4e-13);
            let input = result.inward_heat_flux_w_m2[0].scale(TAU * result.face_radii_m[0])
                + result.inward_heat_flux_w_m2[1].scale(TAU * result.face_radii_m[1]);
            near(input, C64::new(0.0, -omega) * local_energy, 3e-10);
            losses += integrals[0];
            storage += integrals[1];
            energy = energy + local_energy;
        }
        let mut contact_loss = 0.0;
        for contact in 0..2 {
            let value = wall.contact_observation(contact, temperature).unwrap();
            let left = wall.layer_observation(contact, temperature).unwrap();
            let right = wall.layer_observation(contact + 1, temperature).unwrap();
            near(
                value.outward_heat_flux_w_m2,
                -left.inward_heat_flux_w_m2[1],
                3e-10,
            );
            near(
                value.outward_heat_flux_w_m2,
                right.inward_heat_flux_w_m2[0],
                3e-10,
            );
            if contact == 0 {
                assert_eq!(
                    value.inner_layer_temperature_k,
                    value.outer_layer_temperature_k
                );
                assert_eq!(value.contact_dissipation_w_m, 0.0);
            } else {
                let jump = value.inner_layer_temperature_k - value.outer_layer_temperature_k;
                near(
                    jump,
                    value.outward_heat_flux_w_m2.scale(value.resistance_m2_k_w),
                    3e-14,
                );
                relative(
                    value.contact_dissipation_w_m,
                    PI * value.radius_m / 293.15 / value.resistance_m2_k_w * jump.norm_sq(),
                    3e-14,
                );
            }
            contact_loss += value.contact_dissipation_w_m;
        }
        let result = wall.observation(temperature).unwrap();
        relative(losses, result.solid_conduction_w_m, 4e-13);
        relative(contact_loss, result.thermal_contact_w_m, 4e-13);
        relative(storage, result.thermal_storage_j_m, 4e-13);
        near(energy, result.oscillatory_thermal_energy_j_m, 4e-13);
        let input = result
            .fluid_to_wall_heat_flux_w_m2
            .scale(TAU * result.interface_radius_m);
        let output = result
            .wall_to_bath_heat_flux_w_m2
            .scale(TAU * result.bath_radius_m);
        near(input - output, C64::new(0.0, -omega) * energy, 3e-10);
        relative(
            (temperature.conj() * input).re / (2.0 * 293.15),
            result.solid_conduction_w_m + result.thermal_contact_w_m + result.bath_exchange_w_m,
            4e-13,
        );
        relative(
            -(temperature.conj() * input).im / (4.0 * 293.15 * omega),
            result.thermal_storage_j_m,
            4e-13,
        );
    }
}

#[test]
fn layered_inner_and_outer_walls_close_molecular_acoustic_work_with_contact_heat_separate() {
    let sources = [solid_source(0.8, 120_000.0), solid_source(1.8, 210_000.0)];
    let inner = [
        AnnularSection::new(0.0012, 0.0016, 128).unwrap(),
        AnnularSection::new(0.0016, 0.002, 128).unwrap(),
    ];
    let outer = [
        AnnularSection::new(0.003, 0.0034, 128).unwrap(),
        AnnularSection::new(0.0034, 0.0039, 128).unwrap(),
    ];
    let inner_layers: Vec<_> = inner
        .iter()
        .zip(&sources)
        .map(|(section, source)| Layer { section, source })
        .collect();
    let outer_layers: Vec<_> = outer
        .iter()
        .zip(&sources)
        .map(|(section, source)| Layer { section, source })
        .collect();
    let section = AnnularSection::new(0.002, 0.003, 192).unwrap();
    let material = material(2);
    let (p, q) = (C64::new(2.0, -0.7), C64::new(2e-8, 1e-8));
    for omega in [50.0, 500.0, 2500.0] {
        let left = LayeredWall::prepare(
            &inner_layers,
            &[0.0008],
            293.15,
            omega,
            Face::Outer,
            Bath::Conductance {
                conductance_w_m2_k: 40.0,
            },
            solid_admission(),
        )
        .unwrap();
        let right = LayeredWall::prepare(
            &outer_layers,
            &[0.0012],
            293.15,
            omega,
            Face::Inner,
            Bath::Adiabatic,
            solid_admission(),
        )
        .unwrap();
        let fields = section
            .prepare_molecular_conjugate_wall_fields(
                &material,
                omega,
                [CoupledWall::Layered(&left), CoupledWall::Layered(&right)],
                admission(),
            )
            .unwrap();
        let loss = fields.dissipation_per_axial_length(p, q).unwrap();
        let wave = fields.wave();
        let power = 0.5 * (wave.wavenumber * wave.characteristic_impedance).im * q.norm_sq()
            + 0.5 * (wave.wavenumber / wave.characteristic_impedance).im * p.norm_sq();
        relative(loss.total_w_m, power, 4e-10);
        assert!(
            loss.gas_molecular_w_m > 0.0 && loss.solid_conduction_w_m.into_iter().all(|v| v > 0.0)
        );
        assert!(loss.thermal_contact_w_m.into_iter().all(|v| v > 0.0));
        assert_eq!(loss.bath_exchange_w_m[1], 0.0);
        for (side, wall) in [(0, &left), (1, &right)] {
            let observed = fields.wall_observation(side, p).unwrap();
            let expected = wall.observation(observed.interface_temperature_k).unwrap();
            assert_eq!(observed.thermal_contact_w_m, expected.thermal_contact_w_m);
            let actual = fields.thermal_contact_observation(side, 0, p).unwrap();
            assert_eq!(
                actual.contact_dissipation_w_m,
                loss.thermal_contact_w_m[side]
            );
            let local = fields
                .layered_solid_field_in_cell(
                    side,
                    if side == 0 { 1 } else { 0 },
                    if side == 0 { 127 } else { 0 },
                    if side == 0 { 1.0 } else { 0.0 },
                    p,
                )
                .unwrap();
            near(local.temperature_k, observed.interface_temperature_k, 4e-15);
        }
    }
}

#[test]
fn thermally_opaque_layers_do_not_require_inverting_vanishing_through_admittance() {
    let source = solid_source(0.2, 1_000_000.0);
    let first = AnnularSection::new(0.1, 0.12, 2048).unwrap();
    let second = AnnularSection::new(0.12, 0.14, 2048).unwrap();
    let config = RadialThermalWallAdmission {
        maximum_payload_bytes: 3_000_000,
        ..solid_admission()
    };
    let wall = LayeredWall::prepare(
        &[
            Layer {
                section: &first,
                source: &source,
            },
            Layer {
                section: &second,
                source: &source,
            },
        ],
        &[0.001],
        293.15,
        1000.0,
        Face::Inner,
        Bath::Isothermal,
        config,
    )
    .unwrap();
    let original = first
        .prepare_radial_thermal_wall(
            &source,
            293.15,
            1000.0,
            Face::Inner,
            Bath::Adiabatic,
            config,
        )
        .unwrap();
    near(
        wall.interface_admittance_w_m2_k(),
        original.interface_admittance_w_m2_k(),
        2e-11,
    );
    let contact = wall.contact_observation(0, C64::ONE).unwrap();
    assert!(contact.inner_layer_temperature_k.abs() < 1e-100);
    assert!(contact.outer_layer_temperature_k.abs() < 1e-100);
    let cold = wall.layer_observation(1, C64::ONE).unwrap();
    assert!(cold.thermal_storage_j_m.is_finite());
    assert!(cold.conduction_dissipation_w_m.is_finite());
    assert!(wall.observation(C64::ONE).unwrap().solid_conduction_w_m > 0.0);
}

#[test]
fn finite_contacts_approach_perfect_contact_and_isolated_first_layer_without_fitted_limits() {
    let source_a = solid_source(0.5, 100_000.0);
    let source_b = solid_source(2.0, 200_000.0);
    let first = AnnularSection::new(0.003, 0.0033, 96).unwrap();
    let second = AnnularSection::new(0.0033, 0.0037, 96).unwrap();
    let layers = [
        Layer {
            section: &first,
            source: &source_a,
        },
        Layer {
            section: &second,
            source: &source_b,
        },
    ];
    let perfect = LayeredWall::prepare(
        &layers,
        &[0.0],
        293.15,
        400.0,
        Face::Inner,
        Bath::Isothermal,
        solid_admission(),
    )
    .unwrap();
    let mut previous = f64::INFINITY;
    for resistance in [1e-6, 1e-7, 1e-8] {
        let wall = LayeredWall::prepare(
            &layers,
            &[resistance],
            293.15,
            400.0,
            Face::Inner,
            Bath::Isothermal,
            solid_admission(),
        )
        .unwrap();
        let error =
            (wall.interface_admittance_w_m2_k() - perfect.interface_admittance_w_m2_k()).abs();
        assert!(
            error < previous * 0.12,
            "R={resistance}: {error:e}/{previous:e}"
        );
        previous = error;
        assert!(wall.observation(C64::ONE).unwrap().thermal_contact_w_m > 0.0);
    }
    let isolated = first
        .prepare_radial_thermal_wall(
            &source_a,
            293.15,
            400.0,
            Face::Inner,
            Bath::Adiabatic,
            solid_admission(),
        )
        .unwrap();
    let wall = LayeredWall::prepare(
        &layers,
        &[1e8],
        293.15,
        400.0,
        Face::Inner,
        Bath::Isothermal,
        solid_admission(),
    )
    .unwrap();
    near(
        wall.interface_admittance_w_m2_k(),
        isolated.interface_admittance_w_m2_k(),
        2e-10,
    );
}

#[test]
fn layered_empty_reservoir_parity_and_invalid_sources_interfaces_budgets_cancel_cleanly() {
    let source = solid_source(0.6, 100_000.0);
    let first = AnnularSection::new(0.003, 0.0033, 96).unwrap();
    let second = AnnularSection::new(0.0033, 0.0037, 96).unwrap();
    let layers = [
        Layer {
            section: &first,
            source: &source,
        },
        Layer {
            section: &second,
            source: &source,
        },
    ];
    let material = material(0);
    let gas = material.gas();
    let section = AnnularSection::new(0.002, 0.003, 96).unwrap();
    let wall = LayeredWall::prepare(
        &layers,
        &[0.001],
        293.15,
        300.0,
        Face::Inner,
        Bath::Adiabatic,
        solid_admission(),
    )
    .unwrap();
    let boundaries = [
        CoupledWall::Boundary(Wall::Adiabatic),
        CoupledWall::Layered(&wall),
    ];
    let a = section
        .prepare_conjugate_wall_fields(gas, 300.0, boundaries, admission())
        .unwrap();
    let b = section
        .prepare_molecular_conjugate_wall_fields(&material, 300.0, boundaries, admission())
        .unwrap();
    assert_eq!(a.wave().wavenumber, b.wave().wavenumber);
    assert_eq!(
        a.dissipation_per_axial_length(C64::ONE, C64::ZERO).unwrap(),
        b.dissipation_per_axial_length(C64::ONE, C64::ZERO).unwrap()
    );
    for contacts in [&[][..], &[-0.001][..], &[f64::INFINITY][..]] {
        assert!(
            LayeredWall::prepare(
                &layers,
                contacts,
                293.15,
                300.0,
                Face::Inner,
                Bath::Adiabatic,
                solid_admission()
            )
            .is_err()
        );
    }
    for (t, omega) in [(500.1, 300.0), (293.15, 0.0), (f64::NAN, 300.0)] {
        assert!(
            LayeredWall::prepare(
                &layers,
                &[0.001],
                t,
                omega,
                Face::Inner,
                Bath::Adiabatic,
                solid_admission()
            )
            .is_err()
        );
    }
    let gap = AnnularSection::new(0.0034, 0.0037, 96).unwrap();
    assert!(
        LayeredWall::prepare(
            &[
                layers[0],
                Layer {
                    section: &gap,
                    source: &source
                }
            ],
            &[0.001],
            293.15,
            300.0,
            Face::Inner,
            Bath::Adiabatic,
            solid_admission()
        )
        .is_err()
    );
    assert!(
        LayeredWall::prepare(
            &layers,
            &[0.001],
            293.15,
            300.0,
            Face::Inner,
            Bath::Adiabatic,
            RadialThermalWallAdmission {
                maximum_payload_bytes: 0,
                ..solid_admission()
            }
        )
        .is_err()
    );
    assert!(
        LayeredWall::prepare(
            &[layers[0]; 65],
            &[0.0; 64],
            293.15,
            300.0,
            Face::Inner,
            Bath::Adiabatic,
            solid_admission()
        )
        .is_err()
    );
    assert!(matches!(
        LayeredWall::prepare_controlled(
            &layers,
            &[0.001],
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
    LayeredWall::prepare_controlled(
        &layers,
        &[0.001],
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
        LayeredWall::prepare_controlled(
            &layers,
            &[0.001],
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
    assert!(
        section
            .prepare_conjugate_wall_fields(gas, 301.0, boundaries, admission())
            .is_err()
    );
    let wrong_face = LayeredWall::prepare(
        &layers,
        &[0.001],
        293.15,
        300.0,
        Face::Outer,
        Bath::Adiabatic,
        solid_admission(),
    )
    .unwrap();
    assert!(
        section
            .prepare_conjugate_wall_fields(
                gas,
                300.0,
                [
                    CoupledWall::Boundary(Wall::Adiabatic),
                    CoupledWall::Layered(&wrong_face)
                ],
                admission()
            )
            .is_err()
    );
    assert!(wall.field_in_cell(2, 0, 0.5, C64::ONE).is_err());
    assert!(wall.contact_observation(1, C64::ONE).is_err());
    assert!(
        a.layered_solid_field_in_cell(2, 0, 0, 0.5, C64::ONE)
            .is_err()
    );
    assert!(a.thermal_contact_observation(0, 0, C64::ONE).is_err());
}
