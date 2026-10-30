//! G1/G3/G4: spatial heat survives the same geometric mass/momentum handoff.
//! Synthetic source data isolate the equations; no measured validation claim.

use super::*;
use super::super::enthalpy::{
    TetEnthalpyTransferError, TetEnthalpyTransferLimits, TetrahedralLiquidEnthalpyState,
};
use crate::material::free_surface::enthalpy::FreeSurfaceEnthalpyControls;
use fs_material::liquid_temperature::LiquidTemperatureLaw;
use fs_material::phase::{
    EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{QuantitySpec, Temperature};

mod components;
mod rheology;

fn thermal_card(rho: f64, hot_nu: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, values) in [
        ("density", Dims([-3, 1, 0, 0, 0, 0]), [rho; 2]),
        ("dynamic_viscosity", Dims([-1, 1, -1, 0, 0, 0]), [rho * 1e-4, rho * hot_nu]),
        ("specific_heat_capacity", Dims([2, 0, -2, -1, 0, 0]), [2000.0; 2]),
        ("thermal_conductivity", Dims([1, 1, -3, -1, 0, 0]), [0.5; 2]),
    ] {
        let quantity = if name == "thermal_conductivity" {
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::ThermalConductivity, ValueForm::Static,
            ))
        } else {
            QuantitySpec::dimensional(dims)
        };
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::with_quantity(name, quantity),
            value: PropertyValue::Curve {
                abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                knots: vec![(300.0, values[0]), (500.0, values[1])], dims,
            },
            validity: ValidityDomain::unconstrained().with("T", 300.0, 500.0)
                .with("pressure", 100000.0, 100000.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "authored constant-density P1 thermal-transfer fixture".into(),
                license: "test-only".into(), artifact: None,
            },
        }).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: format!("synthetic thermal liquid rho={rho}; hot nu={hot_nu}"),
            phase: "liquid".into(), process: "analytic fixture".into(), revision: 0,
        }, claims, vec![],
    ).unwrap()
}

fn point() -> QueryPoint {
    QueryPoint::new().with("T", 400.0).unwrap().with("pressure", 100000.0).unwrap()
}
fn chart(card: &MaterialCard, rho: f64) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_single_phase(
        card.content_hash(), SolidLiquidPhase::Liquid,
        vec![
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 0.0, temperature_k: 300.0,
                liquid_mass_fraction: 1.0, bulk_density_kg_m3: rho,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 400000.0, temperature_k: 500.0,
                liquid_mass_fraction: 1.0, bulk_density_kg_m3: rho,
            },
        ],
    ).unwrap()
}
fn thermal_limits() -> TetEnthalpyTransferLimits {
    TetEnthalpyTransferLimits {
        mechanical: limits(), max_source_queries: 64, enthalpy_tolerance_j: 1e-11,
    }
}
fn affine_h(p: [f64; 3]) -> f64 {
    200000.0 + 5e6 * (p[0] - 0.004) + 2.5e6 * (p[1] - 0.004)
}
fn source<'a>(
    p: &'a [[f64; 3]], t: &'a [[u32; 4]], v: &'a [[f64; 3]], h: &'a [f64],
) -> TetrahedralLiquidEnthalpyState<'a> {
    TetrahedralLiquidEnthalpyState {
        motion: TetrahedralLiquidState {
            positions_m: p, tetrahedra: t, velocities_m_s: v,
        },
        specific_enthalpy_j_kg: h,
    }
}

#[test]
fn g1_spinning_heated_cube_preserves_independent_integrals_and_advances_native_flow() {
    let card = thermal_card(1000.0, 2e-5);
    let law = LiquidTemperatureLaw::try_new(
        &card, &point(), "T", LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    ).unwrap();
    let surface = resolve_liquid_surface_state(
        &card, &point(), LiquidViscosityBasis::Dynamic, LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    ).unwrap();
    let caloric = chart(&card, 1000.0);
    let (p, t, v) = cube([0.0; 3]);
    let h: Vec<_> = p.iter().copied().map(affine_h).collect();
    let unchanged = (p.clone(), v.clone(), h.clone());
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let scaling = plan_free_surface_from_material(&surface, 0.001, 0.001, 0.1).unwrap();
        let (mut thermal, report) = scaling.free_surface3_from_tetrahedra_with_enthalpy(
            domain(), [0.0; 3], source(&p, &t, &v, &h), &law, &caloric,
            thermal_limits(), FreeSurfaceEnthalpyControls::default(), cx, &gate,
        ).unwrap();
        let mass = 1000.0 * 0.004_f64.powi(3);
        let expected_h = mass * affine_h([0.004; 3]);
        close(report.source_enthalpy_j, expected_h, 1e-13);
        close(report.target_enthalpy_j, expected_h, 1e-12);
        close(report.enthalpy_defect_j, 0.0, 1e-12);
        close(report.mechanical.source.mass_kg, mass, 2e-18);
        close(report.mechanical.target.mass_kg, mass, 2e-18);
        let translation = [0.01, -0.005, 0.002];
        let expected_p = scale(translation, mass);
        let inertia = mass * 0.004_f64.powi(2) / 6.0;
        let expected_l = add(cross([0.004; 3], expected_p), [0.0, 0.0, inertia * 3.0]);
        for axis in 0..3 {
            close(report.mechanical.target.momentum_kg_m_s[axis], expected_p[axis], 2e-20);
            close(report.mechanical.target.angular_momentum_kg_m2_s[axis], expected_l[axis], 2e-22);
        }
        // The required spin correction changes resolved K; no part of this
        // numerical representation change enters the thermal enthalpy.
        assert!(report.mechanical.kinetic_energy_change_j > 0.0);
        for z in 2..6 {
            for y in 2..6 {
                for x in 2..6 {
                    let coordinates = [x, y, z];
                    let expected = affine_h(coordinates.map(|n| (n as f64 + 0.5) * 0.001));
                    close(thermal.specific_enthalpy_at(coordinates).unwrap(), expected, 2e-8);
                    let temperature = 300.0 + expected / 2000.0;
                    assert_eq!(
                        thermal.material_at(coordinates).unwrap().material().card_identity(),
                        card.content_hash(),
                    );
                    close(
                        thermal.material_at(coordinates).unwrap().dynamic_viscosity_pa_s(),
                        0.1 + (temperature - 300.0) * (0.02 - 0.1) / 200.0,
                        2e-14,
                    );
                }
            }
        }
        let stepped = thermal.step(&gate).unwrap();
        assert!(stepped.mass_routes > 0);
        close(stepped.mass_after_kg, mass, 1e-14);
        close(stepped.enthalpy_before_j, expected_h, 1e-12);
        close(stepped.enthalpy_after_j, expected_h, 1e-9);
        close(stepped.external_heat_j, 0.0, 0.0);
        close(stepped.viscous_heat_j, 0.0, 0.0);
        close(stepped.balance_defect_j, 0.0, 1e-9);
    });
    assert_eq!((p, v, h), unchanged);
}

#[test]
fn g1_cut_tet_integrates_affine_enthalpy_over_actual_intersections() {
    let card = thermal_card(1000.0, 2e-5);
    let law = LiquidTemperatureLaw::try_new(
        &card, &point(), "T", LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    ).unwrap();
    let surface = resolve_liquid_surface_state(
        &card, &point(), LiquidViscosityBasis::Dynamic, LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    ).unwrap();
    let caloric = chart(&card, 1000.0);
    let p = [
        [0.0023, 0.0022, 0.0021], [0.0057, 0.0022, 0.0021],
        [0.0023, 0.0056, 0.0021], [0.0023, 0.0022, 0.0055],
    ];
    let t = [[0, 1, 2, 3]];
    let v = p.map(|x| [0.005 + x[0], 2.0 * x[1], -x[2]]);
    let h = p.map(affine_h);
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let scaling = plan_free_surface_from_material(&surface, 0.001, 0.001, 0.1).unwrap();
        let mut bounds = thermal_limits();
        bounds.mechanical.maximum_centroid_shift_m = 0.001;
        let (thermal, report) = scaling.free_surface3_from_tetrahedra_with_enthalpy(
            domain(), [0.0; 3], source(&p, &t, &v, &h), &law, &caloric,
            bounds, FreeSurfaceEnthalpyControls::default(), cx, &gate,
        ).unwrap();
        let mass = 1000.0 * 0.0034_f64.powi(3) / 6.0;
        let expected = mass * h.into_iter().sum::<f64>() / 4.0;
        close(report.source_enthalpy_j, expected, 1e-14);
        close(thermal.total_enthalpy_j().unwrap(), expected, 1e-12);
        assert!(report.mechanical.centroid_shift_m.iter().any(|value| value.abs() > 1e-8));
        assert!(report.mechanical.wet_cells > 1);
    });
}

#[test]
fn g3_transferred_temperature_and_source_viscosity_change_actual_native_motion() {
    let (p, t, _) = cube([0.0; 3]);
    let v: Vec<_> = p.iter().map(|x| [0.01 + 2.0 * (x[1] - 0.004), 0.0, 0.0]).collect();
    let h: Vec<_> = p.iter().copied().map(affine_h).collect();
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let mut outcomes = Vec::new();
        for hot_nu in [2e-5, 8e-5] {
            let card = thermal_card(1000.0, hot_nu);
            let law = LiquidTemperatureLaw::try_new(
                &card, &point(), "T", LiquidViscosityBasis::Dynamic,
                MaterialPropertySelection::SingleClaimOnly,
            ).unwrap();
            let surface = resolve_liquid_surface_state(
                &card, &point(), LiquidViscosityBasis::Dynamic, LiquidCapillarity::Neglected,
                MaterialPropertySelection::SingleClaimOnly,
            ).unwrap();
            let caloric = chart(&card, 1000.0);
            let scaling = plan_free_surface_from_material(&surface, 0.001, 0.001, 0.1).unwrap();
            let (mut thermal, report) = scaling.free_surface3_from_tetrahedra_with_enthalpy(
                domain(), [0.0; 3], source(&p, &t, &v, &h), &law, &caloric,
                thermal_limits(), FreeSurfaceEnthalpyControls::default(), cx, &gate,
            ).unwrap();
            for _ in 0..2 { thermal.step(&gate).unwrap(); }
            close(thermal.total_enthalpy_j().unwrap(), report.source_enthalpy_j, 1e-9);
            outcomes.push(thermal.flow().lattice().grid().state_bits());
        }
        assert_ne!(outcomes[0], outcomes[1]);
    });
}

#[test]
fn g4_hidden_phase_density_source_budget_and_late_refusals_publish_nothing() {
    let card = thermal_card(1000.0, 2e-5);
    let law = LiquidTemperatureLaw::try_new(
        &card, &point(), "T", LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    ).unwrap();
    let surface = resolve_liquid_surface_state(
        &card, &point(), LiquidViscosityBasis::Dynamic, LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    ).unwrap();
    let caloric = chart(&card, 1000.0);
    let (p, t, v) = cube([0.0; 3]);
    let h: Vec<_> = p.iter().map(|x| if x[0] < 0.004 { 100000.0 } else { 300000.0 }).collect();
    let before = (p.clone(), v.clone(), h.clone());
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let scaling = plan_free_surface_from_material(&surface, 0.001, 0.001, 0.1).unwrap();
        let call = |curve: &EquilibriumEnthalpyPhaseCurve, bounds, controls| {
            scaling.free_surface3_from_tetrahedra_with_enthalpy(
                domain(), [0.0; 3], source(&p, &t, &v, &h), &law, curve,
                bounds, controls, cx, &gate,
            ).map(|(_, report)| report)
        };
        // Both nodal endpoints look liquid at the correct density. An interior
        // caloric density change still requires a volume-change owner.
        let hidden_density = EquilibriumEnthalpyPhaseCurve::try_single_phase(
            card.content_hash(), SolidLiquidPhase::Liquid,
            vec![
                EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: 100000.0, temperature_k: 350.0,
                    liquid_mass_fraction: 1.0, bulk_density_kg_m3: 1000.0,
                },
                EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: 200000.0, temperature_k: 400.0,
                    liquid_mass_fraction: 1.0, bulk_density_kg_m3: 900.0,
                },
                EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: 300000.0, temperature_k: 450.0,
                    liquid_mass_fraction: 1.0, bulk_density_kg_m3: 1000.0,
                },
            ],
        ).unwrap();
        assert!(call(&hidden_density, thermal_limits(), Default::default()).is_err());
        let solid = EquilibriumEnthalpyPhaseCurve::try_single_phase(
            card.content_hash(), SolidLiquidPhase::Solid,
            caloric.knots().iter().copied().map(|mut knot| {
                knot.liquid_mass_fraction = 0.0; knot
            }).collect(),
        ).unwrap();
        assert!(call(&solid, thermal_limits(), Default::default()).is_err());
        let wrong_card = thermal_card(1000.0, 8e-5);
        assert!(call(&chart(&wrong_card, 1000.0), thermal_limits(), Default::default()).is_err());
        // A density peak between every nodal, caloric AND receiving-cell
        // temperature must be caught by the original source-line owner.
        let mut claims = ClaimSet::new();
        for (_, original) in card.claims().claims_ordered() {
            let mut claim = original.clone();
            if claim.key.name() == "density" {
                claim.value = PropertyValue::Curve {
                    abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                    knots: vec![
                        (300.0, 1000.0), (416.0, 1000.0), (420.0, 900.0),
                        (424.0, 1000.0), (500.0, 1000.0),
                    ],
                    dims: Dims([-3, 1, 0, 0, 0, 0]),
                };
            }
            claims.insert_claim(claim).unwrap();
        }
        let hidden_card = MaterialCard::assemble(card.id().clone(), claims, vec![]).unwrap();
        let hidden_law = LiquidTemperatureLaw::try_new(
            &hidden_card, &point(), "T", LiquidViscosityBasis::Dynamic,
            MaterialPropertySelection::SingleClaimOnly,
        ).unwrap();
        for temperature in [350.0, 450.0, 362.5, 387.5, 412.5, 437.5] {
            close(hidden_law.at_temperature(Temperature::new(temperature)).unwrap()
                .density_kg_m3(), 1000.0, 0.0);
        }
        close(hidden_law.at_temperature(Temperature::new(420.0)).unwrap()
            .density_kg_m3(), 900.0, 0.0);
        let hidden_surface = resolve_liquid_surface_state(
            &hidden_card, &point(), LiquidViscosityBasis::Dynamic, LiquidCapillarity::Neglected,
            MaterialPropertySelection::SingleClaimOnly,
        ).unwrap();
        let hidden_chart = chart(&hidden_card, 1000.0);
        let hidden_scaling =
            plan_free_surface_from_material(&hidden_surface, 0.001, 0.001, 0.1).unwrap();
        assert!(matches!(
            hidden_scaling.free_surface3_from_tetrahedra_with_enthalpy(
                domain(), [0.0; 3], source(&p, &t, &v, &h), &hidden_law, &hidden_chart,
                thermal_limits(), Default::default(), cx, &gate,
            ),
            Err(TetEnthalpyTransferError::Transfer(TetLiquidTransferError::Invalid(
                "source thermal field requires a density/volume-change owner"
            )))
        ));
        let mut short = thermal_limits(); short.max_source_queries = 2;
        assert!(call(&caloric, short, Default::default()).is_err());
        // Geometry/momenta have already succeeded when thermal attachment
        // refuses this undersized native storage budget.
        let mut late = FreeSurfaceEnthalpyControls::default(); late.max_cells = 1;
        assert!(matches!(
            call(&caloric, thermal_limits(), late),
            Err(TetEnthalpyTransferError::Enthalpy(_))
        ));
        let a = call(&caloric, thermal_limits(), Default::default()).unwrap();
        let b = call(&caloric, thermal_limits(), Default::default()).unwrap();
        assert_eq!(a, b);
        gate.request();
        assert!(matches!(
            call(&caloric, thermal_limits(), Default::default()),
            Err(TetEnthalpyTransferError::Transfer(TetLiquidTransferError::Cancelled))
        ));
    });
    assert_eq!((p, v, h), before);
}
