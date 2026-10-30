//! G1/G3: original silver crystal data reach waves, pressure and native forces.
//! Table 9 is an adiabatic crystal reference. The sole overlapping density
//! query is 293.15 K; no temperature-dependent density or wire texture is made up.

use super::*;
use fs_blake3::ContentHash;
use fs_matdb::{MaterialStateId, NormalizedMaterialCardPack};
use fs_material::elastic::acoustic::ElasticAcousticState;
use fs_material::elastic::acoustic::sourced::cubic::{
    CubicElasticRequirements, resolve_cubic_elastic_state,
};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
};

fn at(t: f64) -> QueryPoint {
    temperature(
        flags(&[
            ("source-pressure-known", 0.0),
            ("source-measurement-frequency-known", 0.0),
            ("engineering-reference-model", 1.0),
        ]),
        "T",
        t,
    )
}

fn requirements() -> CubicElasticRequirements {
    let positive = |name, dims| {
        ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
            .unwrap()
    };
    CubicElasticRequirements::try_new(
        positive("density", RHO),
        positive("cubic_adiabatic_c11", PA),
        ScalarPropertyRequirement::try_new("cubic_adiabatic_c12", PA, ScalarAdmissibility::Finite)
            .unwrap(),
        positive("cubic_adiabatic_c44", PA),
    )
    .unwrap()
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-11 * expected.abs().max(1e-12),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g1_g3_silver_table_reaches_the_same_oriented_wave_and_solid_law() {
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
    use fs_solid::{
        TetAssemblyBudget, TetElasticMaterial, TetLinearElasticProblem, TetMaterialField,
    };

    let pack = compile("silver-pure-smith-fickett-nist1995", 16);
    // Independently transcribed Table 9 facts: T, C11, C12, C44, B [K, GPa].
    for [t, c11, c12, c44, bulk] in [
        [10.0, 131.44, 97.30, 51.08, 108.68],
        [25.0, 131.41, 97.33, 50.98, 108.69],
        [50.0, 131.12, 97.26, 50.72, 108.55],
        [75.0, 130.54, 97.04, 50.28, 108.21],
        [100.0, 129.80, 96.66, 49.82, 107.71],
        [125.0, 129.08, 96.28, 49.36, 107.21],
        [150.0, 128.35, 95.91, 48.90, 106.72],
        [175.0, 127.62, 95.54, 48.44, 106.23],
        [200.0, 126.91, 95.17, 47.97, 105.75],
        [225.0, 126.19, 94.81, 47.50, 105.27],
        [250.0, 125.46, 94.44, 47.04, 104.78],
        [275.0, 124.73, 94.05, 46.58, 104.28],
        [300.0, 123.99, 93.67, 46.12, 103.78],
    ] {
        for (name, value) in [
            ("cubic_adiabatic_c11", c11),
            ("cubic_adiabatic_c12", c12),
            ("cubic_adiabatic_c44", c44),
            ("adiabatic_bulk_modulus", bulk),
        ] {
            sample(&pack, name, PA, &at(t), value * 1e9);
        }
        // The separately printed dependent B column is rounded to 0.01 GPa.
        assert!(((c11 + 2.0 * c12) / 3.0 - bulk).abs() <= 0.0051);
    }
    for bad in [
        at(9.99),
        at(300.01),
        at(293.15).with("source-pressure-known", 1.0).unwrap(),
        at(293.15).with("engineering-reference-model", 0.0).unwrap(),
        at(293.15)
            .with("source-measurement-frequency-known", 1.0)
            .unwrap(),
    ] {
        refuses(&pack, "cubic_adiabatic_c11", PA, &bad);
    }
    let normalized = NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: "Ag; NIST pure-silver reference".into(),
            phase: "cubic crystal".into(),
            process:
                "Table 1 density plus Table 9 adiabatic stiffness; unmatched specimens acknowledged"
                    .into(),
            revision: 0,
        },
        pack,
    )
    .unwrap();
    let state = resolve_cubic_elastic_state(
        normalized.card(),
        &at(293.15),
        &requirements(),
        ContentHash([0xc1; 32]),
        0.001,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert_eq!(state.material().properties().len(), 4);
    for property in state.material().properties() {
        normalized
            .card()
            .claims()
            .verify_receipt(&property.answer().receipt)
            .unwrap();
    }
    // Independent interpolation of the 275/300 K source rows.
    let (rho, c11, c12, c44) = (10492.0, 124.19276e9, 93.77412e9, 46.24604e9);
    let bulk = (c11 + 2.0 * c12) / 3.0;
    let principal = state.acoustic();
    for (direction, expected) in [
        ([1.0, 0.0, 0.0], [c44, c44, c11]),
        (
            [1.0, 1.0, 0.0],
            [(c11 - c12) / 2.0, c44, (c11 + c12 + 2.0 * c44) / 2.0],
        ),
        (
            [1.0, 1.0, 1.0],
            [
                (c11 - c12 + c44) / 3.0,
                (c11 - c12 + c44) / 3.0,
                (c11 + 2.0 * c12 + 4.0 * c44) / 3.0,
            ],
        ),
    ] {
        let actual = principal.directional_waves(direction).unwrap();
        for (mode, modulus) in actual.modes.iter().zip(expected) {
            near(mode.phase_speed_m_s, (modulus / rho).sqrt());
        }
    }
    let pressure = principal.hydrostatic_response(1e6).unwrap();
    near(pressure.volumetric_strain, -1e6 / bulk);
    near(pressure.stored_energy_density_j_m3, 1e12 / (2.0 * bulk));
    // The source data have a narrow common domain. Querying warmed density
    // cannot silently keep the 20 C density or replace the original card.
    assert!(state.at_query_point(normalized.card(), &at(294.0)).is_err());
    let refreshed = state
        .at_query_point(normalized.card(), &at(293.15))
        .unwrap();
    // Refresh pins the original claims, changing the selection-policy receipt
    // without changing the physical law, crystallographic frame or query.
    assert_eq!(
        refreshed.acoustic().density_kg_m3(),
        state.acoustic().density_kg_m3()
    );
    assert_eq!(
        refreshed.acoustic().stiffness_mandel_pa(),
        state.acoustic().stiffness_mandel_pa()
    );
    assert_eq!(refreshed.crystal_frame(), state.crystal_frame());
    assert_eq!(
        refreshed.material().card_identity(),
        state.material().card_identity()
    );
    assert_eq!(refreshed.material().query(), state.material().query());
    assert_eq!(
        refreshed.material().properties().len(),
        state.material().properties().len()
    );
    for (actual, original) in refreshed
        .material()
        .properties()
        .iter()
        .zip(state.material().properties())
    {
        assert_eq!(actual.requirement(), original.requirement());
        assert_eq!(
            actual.answer().receipt.selected,
            original.answer().receipt.selected
        );
    }
    // Both operations now use the same pinned claims and complete query.
    assert_eq!(
        refreshed
            .at_query_point(normalized.card(), &at(293.15))
            .unwrap(),
        refreshed,
    );

    let nodes = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let tets = [[0, 1, 2, 3]];
    let gradients = [
        [-1.0, -1.0, -1.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let epsilon = 1e-4;
    let displacement = [
        0.0, 0.0, 0.0, epsilon, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    ];
    let mut energies = Vec::new();
    for angle in [0.0_f64, core::f64::consts::FRAC_PI_4] {
        let (s, c) = angle.sin_cos();
        let q = [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]];
        let transformed =
            TetElasticMaterial::from_resolved_cubic_state(&state, ContentHash([0x51; 32]), q)
                .unwrap();
        assert_eq!(transformed.source_material_identity(), state.identity());
        let acoustic = transformed.material().acoustic_response(0.001).unwrap();
        let expected_waves = principal.directional_waves([c, -s, 0.0]).unwrap();
        let actual_waves = acoustic.directional_waves([1.0, 0.0, 0.0]).unwrap();
        for (actual, expected) in actual_waves.modes.iter().zip(expected_waves.modes) {
            near(actual.phase_speed_m_s, expected.phase_speed_m_s);
        }
        let problem = TetLinearElasticProblem {
            nodes_m: &nodes,
            tetrahedra: &tets,
            materials: TetMaterialField::Uniform(transformed.material()),
            fixed_dofs: &[],
            budget: TetAssemblyBudget::standard(),
        };
        let gate = CancelGate::new();
        let pool = ArenaPool::new(ArenaConfig::default());
        let assembled = pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 1,
                    kernel_id: 10,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            problem.assemble(&cx).unwrap()
        });
        near(assembled.total_mass_kg, rho / 6.0);
        // Rotate the actual imposed strain into crystal axes, apply the
        // three-constant stress law, rotate back, and integrate P1 gradients.
        let local: [[f64; 3]; 3] = core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                let strain = epsilon * q[0][i] * q[0][j];
                if i == j {
                    c12 * epsilon + (c11 - c12) * strain
                } else {
                    2.0 * c44 * strain
                }
            })
        });
        let stress: [[f64; 3]; 3] = core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                (0..3)
                    .flat_map(|a| (0..3).map(move |b| q[i][a] * local[a][b] * q[j][b]))
                    .sum()
            })
        });
        let k = assembled.stiffness.to_dense();
        let force: [f64; 12] =
            core::array::from_fn(|i| (0..12).map(|j| k[i * 12 + j] * displacement[j]).sum());
        for (node, gradient) in gradients.iter().enumerate() {
            for component in 0..3 {
                let expected = (0..3)
                    .map(|j| stress[component][j] * gradient[j])
                    .sum::<f64>()
                    / 6.0;
                assert!((force[3 * node + component] - expected).abs() < 3e-10 * c11 * epsilon);
            }
        }
        let energy = 0.5
            * force
                .iter()
                .zip(displacement)
                .map(|(f, u)| f * u)
                .sum::<f64>();
        near(energy, 0.5 * epsilon * stress[0][0] / 6.0);
        energies.push(energy);
    }
    assert!(
        (energies[1] / energies[0] - 1.0).abs() > 0.20,
        "crystal orientation must affect the actual native mechanical operator"
    );
}

#[test]
fn g0_g1_cubic_stability_and_isotropic_limit_use_the_shared_tensor_owner() {
    let identity = ContentHash([1; 32]);
    let cubic =
        |c11, c12, c44| ElasticAcousticState::from_cubic(2000.0, c11, c12, c44, identity, 0.001);
    // Cubic stability permits negative C12, but requires C11-C12>0,
    // C11+2C12>0 and C44>0. Do not replace these with C12 positivity.
    assert!(cubic(10e9, -2e9, 3e9).is_ok());
    for (a, b, c) in [
        (10e9, 10e9, 3e9),
        (10e9, -5e9, 3e9),
        (10e9, 2e9, 0.0),
        (f64::NAN, 2e9, 3e9),
        (10e9, f64::INFINITY, 3e9),
    ] {
        assert!(cubic(a, b, c).is_err());
    }
    let isotropic =
        ElasticAcousticState::from_bulk_and_shear(2000.0, 8e9, 3e9, identity, 0.001).unwrap();
    let limit = cubic(12e9, 6e9, 3e9).unwrap();
    for n in [[1.0, 0.0, 0.0], [1.0, -2.0, 3.0]] {
        let a = isotropic.directional_waves(n).unwrap();
        let b = limit.directional_waves(n).unwrap();
        for (a, b) in a.modes.iter().zip(b.modes) {
            near(a.phase_speed_m_s, b.phase_speed_m_s);
        }
    }
}
