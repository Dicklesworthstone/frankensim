//! A CTE source knot drives the actual compatible tetrahedron and caloric heat.
//! All coefficients are manufactured; no measured-material claim is made.
use super::*;

const EXPANSION: [(f64, f64); 5] = [
    (280.0, 0.0002),
    (300.0, 0.0003),
    (302.0, 0.0005),
    (360.0, 0.0004),
    (500.0, 0.0004),
];
const CP: f64 = 10.0;

fn materials() -> (MaterialCard, EquilibriumEnthalpyPhaseCurve) {
    let mut claims = ClaimSet::new();
    for (name, dims, knots) in [
        (
            "young_modulus",
            Pressure::DIMS,
            vec![(280.0, 960.0), (500.0, 1400.0)],
        ),
        (
            "poisson_ratio",
            Dims::NONE,
            vec![(280.0, 0.0), (500.0, 0.0)],
        ),
        (
            "linear_thermal_expansion_coefficient",
            INVERSE_TEMPERATURE_DIMS,
            EXPANSION.to_vec(),
        ),
        ("density", Density::DIMS, vec![(280.0, 6.0), (500.0, 6.0)]),
        (
            "specific_heat_capacity",
            SPECIFIC_HEAT_CAPACITY_DIMS,
            vec![(280.0, CP), (500.0, CP)],
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots,
                    dims,
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    quantity(),
                    280.0,
                    500.0,
                ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "authored table-heated solid".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic table-heated solid".into(),
            phase: "solid".into(),
            process: "native thermoelastic test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let knot = |temperature_k, specific_enthalpy_j_kg| EnthalpyPhaseKnot {
        temperature_k,
        specific_enthalpy_j_kg,
        liquid_mass_fraction: 0.0,
        bulk_density_kg_m3: 6.0,
    };
    let curve = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        card.content_hash(),
        SolidLiquidPhase::Solid,
        vec![knot(280.0, 0.0), knot(500.0, 220.0 * CP)],
    )
    .unwrap();
    (card, curve)
}

fn original_expansion(t: f64) -> (f64, f64) {
    let (lo, hi, sign) = if t >= 300.0 {
        (300.0, t, 1.0)
    } else {
        (t, 300.0, -1.0)
    };
    let mut free = 0.0;
    let mut alpha = None;
    for pair in EXPANSION.windows(2) {
        let [(a, x), (b, y)] = [pair[0], pair[1]];
        let slope = (y - x) / (b - a);
        let left = lo.max(a);
        let right = hi.min(b);
        if right > left {
            let width = right - left;
            free += (x + slope * (left - a)) * width + 0.5 * slope * width * width;
        }
        if a <= t && t <= b {
            alpha = Some(x + slope * (t - a));
        }
    }
    (alpha.unwrap(), sign * free)
}

fn clamped_internal_energy(t: f64) -> f64 {
    let (alpha, free) = original_expansion(t);
    let young = 1000.0 + 2.0 * (t - 300.0);
    let psi = 1.5 * VOLUME * young * free * free;
    let entropy = VOLUME * (-3.0 * young * alpha * free - 3.0 * free * free);
    psi + t * entropy
}

#[test]
fn original_tabulated_expansion_crosses_knots_in_free_and_clamped_native_heat_feedback() {
    with_cx(|cx| {
        for clamped in [false, true] {
            let (card, curve) = materials();
            let thermal = thermal(&[&curve], &[CP], &[400.0], 0.0);
            let context = context(300.0); // Exactly at a genuine CTE slope knot.
            let source = SolidEnthalpy::try_new(
                &thermal,
                &NODES,
                &TETS,
                vec![source(30, &card, &curve, &context)],
                TetAssemblyBudget::standard(),
                0.05,
                mechanical().mesh_frame,
                cx,
            )
            .unwrap();
            let loads = source
                .geometry()
                .assemble(&[], VolumeLoad3::None, &[], cx)
                .unwrap();
            let fixed = if clamped {
                (0..12).collect::<Vec<_>>()
            } else {
                FREE.to_vec()
            };
            let body = ThermoelasticSolid::try_new(
                &source,
                &thermal,
                &loads,
                &fixed,
                &[[280.0, 360.0]],
                mechanical(),
                solve(),
                cx,
            )
            .unwrap();
            let next = body
                .advance(&thermal, &(), 0.1, &loads, thermal_controls(), solve(), cx)
                .unwrap();
            let expected_t = if clamped {
                root(280.0, 360.0, |t| {
                    CP * (t - 300.0) + clamped_internal_energy(t) - 40.0
                })
            } else {
                304.0
            };
            let actual_t = next.thermal.temperature_k(IDS[0]).unwrap();
            near(actual_t, expected_t, 3e-8);
            assert!(
                actual_t > 302.0,
                "the accepted physical interval must cross the source knot"
            );
            let free = original_expansion(expected_t).1;
            let u = next
                .solid
                .equilibrium()
                .solution
                .displacement
                .displacement_m();
            for (node, point) in NODES.iter().enumerate() {
                for k in 0..3 {
                    near(
                        u[node][k],
                        if clamped { 0.0 } else { free * point[k] },
                        3e-10,
                    );
                }
            }
            let stress = next.solid.equilibrium().solution.stress.elements()[0].stress_mandel_pa();
            for &value in &stress[..3] {
                near(
                    value,
                    if clamped {
                        -(1000.0 + 2.0 * (expected_t - 300.0)) * free
                    } else {
                        0.0
                    },
                    3e-8,
                );
            }
            near(next.report.external_heat_j, 40.0, 1e-9);
            near(next.report.mechanical_work_j, 0.0, 1e-12);
            near(next.report.combined_residual_j, 0.0, 1e-9);
            near(
                next.report.internal_energy_change_j,
                if clamped {
                    clamped_internal_energy(actual_t)
                } else {
                    0.0
                },
                2e-9,
            );
            near(
                next.thermal.energy_j() - thermal.energy_j() + next.report.internal_energy_change_j,
                40.0,
                1e-9,
            );
            assert_eq!(
                next.thermal.cells()[&IDS[0]].mass_kg,
                thermal.cells()[&IDS[0]].mass_kg
            );
            assert_eq!(
                next.solid.state().samples()[0].elastic.card_identity(),
                card.content_hash()
            );
            assert_eq!(thermal.time_s(), 0.0);
            assert_eq!(body.time_s(), 0.0);
        }
    });
}
