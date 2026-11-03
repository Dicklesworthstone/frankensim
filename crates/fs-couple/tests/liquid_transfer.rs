//! G1/G3/G4: actual finite liquid inventories exchange momentum and mixing heat.
//! Constitutive data are synthetic; conservation references are independent.

use fs_couple::liquid_transfer::{LiquidMassTransfer, transfer_interface_mass};
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_lbm::d3q19::freesurface3::ContactModel3;
use fs_lbm::material::free_surface::enthalpy::{
    FreeSurfaceEnthalpyControls, FreeSurfaceMassExchangeControls, MaterialFreeSurfaceEnthalpy3,
};
use fs_lbm::material::free_surface::{
    MaterialFreeSurfaceDomain3, MaterialFreeSurfaceInitialCell3, plan_free_surface_from_material,
};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::liquid::surface::{LiquidCapillarity, resolve_liquid_surface_state};
use fs_material::liquid_temperature::{LiquidTemperatureLaw, LiquidViscosityBasis};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec, Temperature};

const CELL: [usize; 3] = [3; 3];
const ORIGIN: [f64; 3] = [0.125, -0.25, 0.5];
const DX: f64 = 0.001;

fn controls() -> FreeSurfaceMassExchangeControls {
    FreeSurfaceMassExchangeControls {
        mass_tolerance_kg: 1e-18,
        momentum_tolerance_kg_m_s: 1e-18,
        energy_tolerance_j: 1e-12,
    }
}

fn transfer(mass_kg: f64) -> LiquidMassTransfer {
    LiquidMassTransfer {
        donor: CELL,
        receiver: CELL,
        donor_origin_m: ORIGIN,
        receiver_origin_m: ORIGIN,
        mass_kg,
    }
}

fn with_pair(
    datum: f64,
    boost: [f64; 3],
    receiver_fill: f64,
    receiver_speed_bound: f64,
    f: impl FnOnce(
        &mut MaterialFreeSurfaceEnthalpy3<'_, '_>,
        &mut MaterialFreeSurfaceEnthalpy3<'_, '_>,
        &CancelGate,
    ),
) {
    let mut claims = ClaimSet::new();
    for (name, dims, values) in [
        ("density", Dims([-3, 1, 0, 0, 0, 0]), [1000.0; 2]),
        (
            "specific_heat_capacity",
            Dims([2, 0, -2, -1, 0, 0]),
            [2000.0; 2],
        ),
        ("thermal_conductivity", Dims([1, 1, -3, -1, 0, 0]), [0.5; 2]),
        ("dynamic_viscosity", Dims([-1, 1, -1, 0, 0, 0]), [0.1, 0.02]),
    ] {
        let quantity = if name == "thermal_conductivity" {
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::ThermalConductivity,
                ValueForm::Static,
            ))
        } else {
            QuantitySpec::dimensional(dims)
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, quantity),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(300.0, values[0]), (500.0, values[1])],
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with("T", 300.0, 500.0)
                    .with("pressure", 100000.0, 100000.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic closed liquid transfer".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic liquid".into(),
            phase: "liquid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let point = QueryPoint::new()
        .with("T", 400.0)
        .unwrap()
        .with("pressure", 100000.0)
        .unwrap();
    let law = LiquidTemperatureLaw::try_new(
        &card,
        &point,
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let source = resolve_liquid_surface_state(
        &card,
        &point,
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let chart = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        card.content_hash(),
        SolidLiquidPhase::Liquid,
        vec![
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: datum,
                temperature_k: 300.0,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: 1000.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: datum + 400_000.0,
                temperature_k: 500.0,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: 1000.0,
            },
        ],
    )
    .unwrap();
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = fs_exec::Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 71,
                kernel_id: 29,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        let make = |fill, velocity: [f64; 3], h, maximum_speed| {
            let velocity_m_s = core::array::from_fn(|d| velocity[d] + boost[d]);
            let mut cells = Vec::new();
            for z in 1..=6 {
                for y in 1..=6 {
                    for x in 1..=3 {
                        cells.push(MaterialFreeSurfaceInitialCell3 {
                            coordinates: [x, y, z],
                            fill_fraction: if [x, y, z] == CELL { fill } else { 1.0 },
                            velocity_m_s,
                        });
                    }
                }
            }
            plan_free_surface_from_material(&source, DX, 0.001, maximum_speed)
                .unwrap()
                .free_surface3_from_cells(
                    MaterialFreeSurfaceDomain3 {
                        cells: [8; 3],
                        gravity_m_s2: [0.0; 3],
                        contact: ContactModel3::Neutral,
                        solid_tiles: &[],
                        max_cells: 512,
                    },
                    &cells,
                    &cx,
                )
                .unwrap()
                .with_enthalpy(
                    &law,
                    &chart,
                    FreeSurfaceEnthalpyControls::default(),
                    |_| datum + h,
                    &gate,
                )
                .unwrap()
        };
        let mut donor = make(0.6, [0.04, -0.02, 0.01], 200_000.0, 0.1);
        let mut receiver = make(
            receiver_fill,
            [-0.03, 0.02, -0.01],
            100_000.0,
            receiver_speed_bound,
        );
        f(&mut donor, &mut receiver, &gate);
    });
}

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        a.is_finite() && b.is_finite() && (a - b).abs() <= tolerance,
        "{a:e} != {b:e}"
    );
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

// Direct original state readback, independently of native transfer receipts.
fn totals(liquid: &MaterialFreeSurfaceEnthalpy3<'_, '_>) -> [f64; 8] {
    let mut result = [0.0; 8]; // mass, H+K, P[3], L[3]
    for z in 0..8 {
        for y in 0..8 {
            for x in 0..8 {
                if let Some(cell) = liquid.flow().sample_at_si(x, y, z).unwrap() {
                    let mass = cell.density_kg_m3 * cell.fill_fraction * DX.powi(3);
                    let h = liquid.specific_enthalpy_at([x, y, z]).unwrap();
                    let p = cell.velocity_m_s.map(|v| mass * v);
                    let r = core::array::from_fn(|d| ([x, y, z][d] as f64 + 0.5) * DX + ORIGIN[d]);
                    let angular = cross(r, p);
                    result[0] += mass;
                    result[1] +=
                        mass * (h + 0.5 * cell.velocity_m_s.iter().map(|v| v * v).sum::<f64>());
                    for d in 0..3 {
                        result[d + 2] += p[d];
                        result[d + 5] += angular[d];
                    }
                }
            }
        }
    }
    result
}

fn combined(
    a: &MaterialFreeSurfaceEnthalpy3<'_, '_>,
    b: &MaterialFreeSurfaceEnthalpy3<'_, '_>,
) -> [f64; 8] {
    let (a, b) = (totals(a), totals(b));
    core::array::from_fn(|d| a[d] + b[d])
}

fn unchanged(a: &MaterialFreeSurfaceEnthalpy3<'_, '_>, b: &MaterialFreeSurfaceEnthalpy3<'_, '_>) {
    assert_eq!(
        a.flow().lattice().grid().state_bits(),
        b.flow().lattice().grid().state_bits()
    );
    assert_eq!(
        a.flow().ledger_mass_kg().unwrap(),
        b.flow().ledger_mass_kg().unwrap()
    );
    assert_eq!(a.total_enthalpy_j().unwrap(), b.total_enthalpy_j().unwrap());
    assert_eq!(
        a.flow().sample_at_si(3, 3, 3).unwrap(),
        b.flow().sample_at_si(3, 3, 3).unwrap()
    );
    assert_eq!(a.material_at(CELL), b.material_at(CELL));
    assert_eq!(a.specific_enthalpy_at(CELL), b.specific_enthalpy_at(CELL));
    assert_eq!(a.physical_time_s(), b.physical_time_s());
    assert_eq!(a.accepted_steps(), b.accepted_steps());
}

#[test]
fn g1_unequal_finite_liquids_exchange_real_momentum_and_thermalize_mixing_once() {
    with_pair(0.0, [0.0; 3], 0.3, 0.1, |donor, receiver, gate| {
        let initial = combined(donor, receiver);
        let a = donor.flow().sample_at_si(3, 3, 3).unwrap().unwrap();
        let b = receiver.flow().sample_at_si(3, 3, 3).unwrap().unwrap();
        let m = b.density_kg_m3 * b.fill_fraction * DX.powi(3);
        let report =
            transfer_interface_mass(donor, receiver, transfer(1e-7), controls(), gate).unwrap();
        let dm = report.receiver.exchange.mass_change_kg;
        near(dm, -report.donor.exchange.mass_change_kg, 1e-21);
        let after = receiver.flow().sample_at_si(3, 3, 3).unwrap().unwrap();
        let heat = 0.5 * m * dm / (m + dm)
            * (0..3)
                .map(|d| (a.velocity_m_s[d] - b.velocity_m_s[d]).powi(2))
                .sum::<f64>();
        for d in 0..3 {
            near(
                after.velocity_m_s[d],
                (m * b.velocity_m_s[d] + dm * a.velocity_m_s[d]) / (m + dm),
                2e-14,
            );
        }
        near(
            receiver.specific_enthalpy_at(CELL).unwrap(),
            (m * 100_000.0 + dm * 200_000.0 + heat) / (m + dm),
            5e-10,
        );
        near(report.receiver.mixing_heat_j, heat, 1e-23);
        assert!(report.receiver.mixing_heat_j > 0.0);
        assert_eq!(report.donor.mechanical_work_j, 0.0);
        assert_eq!(report.donor.exchange.heat_j, 0.0);
        assert_eq!(report.receiver.exchange.heat_j, 0.0);
        assert_eq!(
            donor
                .flow()
                .sample_at_si(3, 3, 3)
                .unwrap()
                .unwrap()
                .velocity_m_s,
            a.velocity_m_s
        );
        assert_eq!(donor.physical_time_s(), 0.0);
        assert_eq!(receiver.physical_time_s(), 0.0);
        for (d, (a, b)) in combined(donor, receiver)
            .into_iter()
            .zip(initial)
            .enumerate()
        {
            near(a, b, if d == 1 { 5e-12 } else { 2e-18 });
        }
        assert!(report.momentum_defect_kg_m_s.iter().all(|x| x.is_finite()));
        assert!(report.enthalpy_kinetic_defect_j.is_finite());
    });
}

#[test]
fn g3_transfer_preserves_boosts_caloric_datums_and_chained_native_histories() {
    let mut reference_heat = None;
    for (datum, boost) in [(0.0, [0.0; 3]), (1e6, [0.02, -0.01, 0.01])] {
        with_pair(datum, boost, 0.3, 0.1, |donor, receiver, gate| {
            let before = combined(donor, receiver);
            let first =
                transfer_interface_mass(donor, receiver, transfer(1e-7), controls(), gate).unwrap();
            if let Some(heat) = reference_heat {
                near(first.receiver.mixing_heat_j, heat, 2e-23);
            } else {
                reference_heat = Some(first.receiver.mixing_heat_j);
            }
            let second =
                transfer_interface_mass(receiver, donor, transfer(5e-8), controls(), gate).unwrap();
            assert!(second.receiver.mixing_heat_j > 0.0);
            for (d, (a, b)) in combined(donor, receiver)
                .into_iter()
                .zip(before)
                .enumerate()
            {
                near(a, b, if d == 1 { 2e-10 } else { 2e-18 });
            }
            let energy = [
                donor.total_enthalpy_j().unwrap(),
                receiver.total_enthalpy_j().unwrap(),
            ];
            let a = donor.step(gate).unwrap();
            let b = receiver.step(gate).unwrap();
            near(a.enthalpy_before_j, energy[0], 1e-11);
            near(b.enthalpy_before_j, energy[1], 1e-11);
            assert_eq!(donor.accepted_steps(), 1);
            assert_eq!(receiver.accepted_steps(), 1);
            assert_eq!(donor.physical_time_s(), 0.001);
            assert_eq!(receiver.physical_time_s(), 0.001);
        });
    }
}

#[test]
fn g4_capacity_source_clock_geometry_and_cancellation_preserve_both_owners() {
    with_pair(0.0, [0.0; 3], 0.8, 0.1, |donor, receiver, gate| {
        let (a, b) = (donor.clone(), receiver.clone());
        assert!(
            transfer_interface_mass(donor, receiver, transfer(3e-7), controls(), gate).is_err()
        );
        let mut separated = transfer(1e-7);
        separated.receiver_origin_m[0] += DX;
        assert!(transfer_interface_mass(donor, receiver, separated, controls(), gate).is_err());
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(
            transfer_interface_mass(donor, receiver, transfer(1e-7), controls(), &cancelled)
                .is_err()
        );
        unchanged(donor, &a);
        unchanged(receiver, &b);
        let mut shifted = receiver.clone().with_time_origin(0.5).unwrap();
        let shift_before = shifted.clone();
        assert!(
            transfer_interface_mass(donor, &mut shifted, transfer(1e-7), controls(), gate).is_err()
        );
        unchanged(donor, &a);
        unchanged(&shifted, &shift_before);
        with_pair(0.0, [0.0; 3], 0.3, 0.1, |other, _, _| {
            let before = other.clone();
            assert!(!donor.shares_mass_transfer_source(other));
            assert!(
                transfer_interface_mass(donor, other, transfer(1e-7), controls(), gate).is_err()
            );
            unchanged(donor, &a);
            unchanged(other, &before);
        });
    });
    // The donor candidate succeeds, then the receiver rejects the outgoing
    // stream's speed. The original donor must still retain its complete mass.
    with_pair(0.0, [0.0; 3], 0.3, 0.04, |donor, receiver, gate| {
        let (a, b) = (donor.clone(), receiver.clone());
        assert!(
            transfer_interface_mass(donor, receiver, transfer(1e-7), controls(), gate).is_err()
        );
        unchanged(donor, &a);
        unchanged(receiver, &b);
    });
}
