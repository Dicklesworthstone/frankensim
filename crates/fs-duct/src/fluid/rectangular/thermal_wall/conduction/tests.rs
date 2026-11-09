//! G1/G3 source-to-field integration. Material claims below are synthetic.
use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::gas::GasSpec;
use fs_material::relaxation::{
    InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw, RelaxingMedium,
    SpeciesRelaxation, VibrationalMode,
};
use fs_material::thermal::transport::{PreparedThermalTransport, ThermalTransportSourceQuery};
use fs_qty::{Dims, QuantitySpec};
use fs_thermochem::SpeciesId;

fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}
fn source(principal: [f64; 2], capacity: f64) -> PreparedOrthotropicThermalTransport<2> {
    let keys = [
        PropertyKey::new("synthetic-density", Dims([-3, 1, 0, 0, 0, 0])),
        PropertyKey::new("synthetic-constrained-capacity", Dims([2, 0, -2, -1, 0, 0])),
        PropertyKey::new("synthetic-principal-k1", Dims([1, 1, -3, -1, 0, 0])),
        PropertyKey::new("synthetic-principal-k2", Dims([1, 1, -3, -1, 0, 0])),
    ];
    let mut claims = ClaimSet::new();
    for (key, value) in keys
        .iter()
        .zip([2500.0, capacity / 2500.0, principal[0], principal[1]])
    {
        claims
            .insert_claim(PropertyClaim {
                key: key.clone(),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Dims([0, 0, 0, 1, 0, 0]),
                    knots: vec![(200.0, value), (500.0, value)],
                    dims: key.quantity().dims(),
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    QuantitySpec::dimensional(Dims([0, 0, 0, 1, 0, 0])),
                    200.0,
                    500.0,
                ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic 2D solid jacket".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic orthotropic solid".into(),
            phase: "solid".into(),
            process: "explicit in-plane material axes".into(),
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
    let paths = core::array::from_fn(|axis| {
        PreparedThermalTransport::try_new(
            &card,
            ThermalTransportSourceQuery {
                density: &keys[0],
                specific_heat: &keys[1],
                conductivity: &keys[2 + axis],
                temperature_axis: "T",
                lower: &point(200.0),
                upper: &point(500.0),
                selections: [ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly); 3],
            },
        )
        .unwrap()
    });
    let (s, c) = 0.37_f64.sin_cos();
    PreparedOrthotropicThermalTransport::try_new(paths, [[c, -s], [s, c]]).unwrap()
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
fn admission() -> RectangularAdmission {
    RectangularAdmission {
        maximum_relative_tail: 2e-5,
        maximum_transverse_phase: 0.8,
    }
}
fn limits() -> RectangularThermalLimits {
    RectangularThermalLimits {
        cells: [24, 32],
        maximum_cell_diffusion_phase: 0.9,
        maximum_relative_residual: 1e-7,
        maximum_payload_bytes: 512 * 1024 * 1024,
        sparse: fs_sparse::direct::DirectLimits {
            max_symbolic_entries: 8_000_000,
            max_analysis_work: 500_000_000,
            max_numeric_work: 1_000_000_000_000,
            max_front_dim: 4096,
        },
    }
}
fn layer(source: &PreparedOrthotropicThermalTransport<2>) -> RectangularSolidLayer<'_> {
    RectangularSolidLayer {
        source,
        thickness_m: [0.0001; 4],
        cells: [4; 4],
    }
}
fn relative(a: C64, b: C64) -> f64 {
    (a - b).abs() / b.abs().max(1e-300)
}

#[test]
fn original_gas_molecular_and_rotated_solid_sources_close_physical_port_work() {
    let material = material(2);
    let inner = source([1.8, 0.7], 1e6);
    let outer = source([0.8, 1.4], 7e5);
    let section = RectangularSection::new(0.002, 0.003, 256).unwrap();
    let layers = [
        layer(&inner),
        RectangularSolidLayer {
            source: &outer,
            thickness_m: [0.00005; 4],
            cells: [2; 4],
        },
    ];
    let baths = [
        RectangularJacketBath::Conductance {
            conductance_w_m2_k: 25.0,
        },
        RectangularJacketBath::Conductance {
            conductance_w_m2_k: 80.0,
        },
        RectangularJacketBath::Adiabatic,
        RectangularJacketBath::Isothermal,
    ];
    let fields = section
        .prepare_molecular_conjugate_jacket_fields(
            &material,
            40.0,
            &layers,
            baths,
            admission(),
            limits(),
        )
        .unwrap();
    assert!(core::ptr::eq(fields.source(), &material));
    assert!(core::ptr::eq(fields.layers()[0].source, &inner));
    assert_eq!(fields.solid_samples()[0].path_identity(), inner.identity());
    let p = C64::new(1.2, -0.4);
    let q = C64::new(2e-8, 1e-8);
    let work = fields.dissipation_per_axial_length(p, q).unwrap();
    let wave = fields.wave();
    let series = C64::new(0.0, 1.0) * (wave.wavenumber * wave.characteristic_impedance).conj();
    let shunt = C64::new(0.0, 1.0) * (wave.wavenumber / wave.characteristic_impedance).conj();
    let expected = 0.5 * (series.re * q.norm_sq() + shunt.re * p.norm_sq());
    assert!((work.total_w_m - expected).abs() < 1e-7 * expected);
    assert!(
        [
            work.viscous_w_m,
            work.gas_conduction_w_m,
            work.molecular_w_m,
            work.solid_conduction_w_m,
            work.baths_w_m[0],
            work.baths_w_m[1]
        ]
        .into_iter()
        .all(|v| v > 0.0)
    );
    assert_eq!(work.baths_w_m[2], 0.0);
    assert_eq!(work.baths_w_m[3], 0.0);
    let a = fields.layer_observation(0, p).unwrap();
    let b = fields.layer_observation(1, p).unwrap();
    assert!(
        (a.conduction_dissipation_w_m + b.conduction_dissipation_w_m - work.solid_conduction_w_m)
            .abs()
            < 1e-12 * work.solid_conduction_w_m
    );
    assert!(a.thermal_storage_j_m > 0.0 && b.thermal_storage_j_m > 0.0);
    let gas = fields.gas_field_at([0.0, 0.0013], p, q).unwrap();
    let solid = fields.solid_field_at(0, [0.0, 0.0013], p).unwrap();
    assert!(relative(gas.gas_temperature_k, solid.temperature_k) < 1e-12);
    let corner = fields.solid_field_at(0, [-0.00003, -0.00004], p).unwrap();
    assert!(
        corner
            .temperature_gradient_k_m
            .iter()
            .all(|v| v.abs() > 1e-3)
    );
    let k = fields.solid_samples()[0].conductivity_w_m_k();
    assert!(k[0][1].abs() > 0.1);
    assert!(
        relative(
            corner.heat_flux_w_m2[0],
            -corner.temperature_gradient_k_m[0].scale(k[0][0])
        ) > 1e-3
    );
    for side in 0..4 {
        for fraction in [0.0, 0.43, 1.0] {
            let bath = fields.outer_wall_observation(side, fraction, p).unwrap();
            if side < 2 {
                assert!(
                    relative(
                        bath.wall_to_bath_heat_flux_w_m2,
                        bath.temperature_k
                            .scale(if side == 0 { 25.0 } else { 80.0 })
                    ) < 1e-12
                        || bath.temperature_k == C64::ZERO
                );
            }
            if side == 2 {
                assert_eq!(bath.wall_to_bath_heat_flux_w_m2, C64::ZERO);
            }
            if side == 3 {
                assert!(bath.temperature_k.abs() < 1e-18);
            }
        }
    }
    assert!(fields.total_wall_to_bath_heat_flux_w_m(p).unwrap().abs() > 0.0);
    // Independent finite-length port power versus the positive field forms.
    let length = 0.1;
    let reference = 1e8_f64;
    let port = fields.two_port(length, reference).unwrap();
    assert_eq!(port.scattering[1], port.scattering[2]);
    let incident = [C64::new(0.0001, 0.00004), C64::new(-0.00003, 0.00002)];
    let outgoing = [
        port.scattering[0] * incident[0] + port.scattering[1] * incident[1],
        port.scattering[2] * incident[0] + port.scattering[3] * incident[1],
    ];
    let port_work = 0.5
        * (incident[0].norm_sq() + incident[1].norm_sq()
            - outgoing[0].norm_sq()
            - outgoing[1].norm_sq());
    assert!(port_work > 0.0);
    let p0 = (incident[0] + outgoing[0]).scale(reference.sqrt());
    let q0 = (incident[0] - outgoing[0]).scale(1.0 / reference.sqrt());
    let pl = (incident[1] + outgoing[1]).scale(reference.sqrt());
    let ql = (outgoing[1] - incident[1]).scale(1.0 / reference.sqrt());
    let forward = (p0 + wave.characteristic_impedance * q0).scale(0.5);
    let backward = (pl - wave.characteristic_impedance * ql).scale(0.5);
    let propagator = |distance: f64| {
        let phase = wave.wavenumber.scale(distance);
        C64::new(phase.re.cos(), phase.re.sin()).scale((-phase.im).exp())
    };
    let panels = 128;
    let mut integral = 0.0;
    for i in 0..=panels {
        let x = length * i as f64 / panels as f64;
        let f = forward * propagator(x);
        let b = backward * propagator(length - x);
        let loss = fields
            .dissipation_per_axial_length(f + b, (f - b) / wave.characteristic_impedance)
            .unwrap()
            .total_w_m;
        integral += loss
            * if i == 0 || i == panels {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
    }
    integral *= length / (3.0 * panels as f64);
    assert!((integral - port_work).abs() < 1e-7 * port_work);
}

#[test]
fn empty_reservoirs_and_original_axis_order_preserve_all_source_bound_fields() {
    let gas = gas();
    let empty = material(0);
    let solid = source([1.8, 0.7], 1e6);
    let a = RectangularSection::new(0.002, 0.003, 256).unwrap();
    let b = RectangularSection::new(0.003, 0.002, 256).unwrap();
    let layers = [layer(&solid)];
    let baths = [RectangularJacketBath::Adiabatic; 4];
    let x = a
        .prepare_conjugate_jacket_fields(&gas, 40.0, &layers, baths, admission(), limits())
        .unwrap();
    let y = b
        .prepare_molecular_conjugate_jacket_fields(
            &empty,
            40.0,
            &layers,
            baths,
            admission(),
            limits(),
        )
        .unwrap();
    assert_eq!(
        x.mean_gas_temperature_per_pressure_k_pa(),
        y.mean_gas_temperature_per_pressure_k_pa()
    );
    assert_eq!(x.wave().wavenumber, y.wave().wavenumber);
    let p = C64::new(0.8, 0.3);
    let q = C64::new(1e-8, -2e-8);
    assert_eq!(
        x.gas_field_at([0.0007, 0.0013], p, q).unwrap(),
        y.gas_field_at([0.0007, 0.0013], p, q).unwrap()
    );
    assert_eq!(
        x.solid_field_at(0, [-0.00004, 0.0013], p).unwrap(),
        y.solid_field_at(0, [-0.00004, 0.0013], p).unwrap()
    );
    assert_eq!(
        x.dissipation_per_axial_length(C64::ZERO, C64::ZERO)
            .unwrap(),
        RectangularConjugateDissipation::default()
    );
    assert_eq!(
        x.layer_observation(0, C64::ZERO)
            .unwrap()
            .thermal_storage_j_m,
        0.0
    );
    assert_eq!(x.total_wall_to_bath_heat_flux_w_m(p).unwrap(), C64::ZERO);
}

#[test]
fn original_material_support_and_preparation_refusals_preserve_sources() {
    let gas = gas();
    let solid = source([1.8, 0.7], 1e6);
    let identity = solid.identity();
    let section = RectangularSection::new(0.002, 0.003, 256).unwrap();
    let layers = [layer(&solid)];
    let baths = [RectangularJacketBath::Isothermal; 4];
    let hot = GasState::try_new(&GasSpec::dry_air_ussa1976(), 600.0, 101325.0).unwrap();
    assert!(
        section
            .prepare_conjugate_jacket_fields(&hot, 40.0, &layers, baths, admission(), limits())
            .is_err()
    );
    let mut stale = gas;
    stale.temperature += 1.0;
    assert!(
        section
            .prepare_conjugate_jacket_fields(&stale, 40.0, &layers, baths, admission(), limits())
            .is_err()
    );
    let mut small = limits();
    small.maximum_payload_bytes = 32;
    assert!(
        section
            .prepare_conjugate_jacket_fields(&gas, 40.0, &layers, baths, admission(), small)
            .is_err()
    );
    let mut count = 0;
    assert!(
        section
            .prepare_conjugate_jacket_fields_controlled(
                &gas,
                40.0,
                &layers,
                baths,
                admission(),
                limits(),
                || {
                    count += 1;
                    count > 20
                }
            )
            .is_err()
    );
    assert_eq!(solid.identity(), identity);
    let a = section
        .prepare_conjugate_jacket_fields(&gas, 40.0, &layers, baths, admission(), limits())
        .unwrap();
    let b = section
        .prepare_conjugate_jacket_fields(&gas, 40.0, &layers, baths, admission(), limits())
        .unwrap();
    assert_eq!(a.wave().wavenumber, b.wave().wavenumber);
    assert!(a.solid_field_at(0, [0.001, 0.001], C64::ONE).is_err());
    assert!(
        a.gas_field_at([-0.00001, 0.001], C64::ONE, C64::ZERO)
            .is_err()
    );
    assert!(a.outer_wall_observation(4, 0.5, C64::ONE).is_err());
}
