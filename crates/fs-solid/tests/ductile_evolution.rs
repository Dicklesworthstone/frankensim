//! Actual tetrahedral reactions and history under progressive ductile softening.
//! Synthetic coefficients verify the numerical model, not a calibrated alloy.
#![deny(unsafe_code)]

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::IsotropicElastic;
use fs_material::ductile::evolution::{DuctileJ2, DuctileJ2Settings};
use fs_material::ductile::{DuctileDomain, TabulatedDuctileLocus};
use fs_material::plastic::J2Plasticity;
use fs_material::plastic::energy::EnergyTolerances;
use fs_solid::brittle::TetCrackBandGeometry;
use fs_solid::ductile::evolution::{
    DuctileEvolutionControls, DuctileTetEvolution, DuctileTetEvolutionError,
    DuctileTetEvolutionStep,
};
use fs_solid::ductile::{DuctileStepConditions, DuctileTemperatures};
use fs_solid::loads3::TetLoadGeometry3;
use fs_solid::{IncrementSettings, TetAssemblyBudget};

const NODES: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const GP: f64 = 0.004;
const INITIATION_STRAIN: f64 = 0.002;
const DAMAGE_LIMIT: f64 = 0.04;

fn with_cx<R>(f: impl FnOnce(&CancelGate, &Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 31,
                kernel_id: 89,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&gate, &cx)
    })
}

fn model(width: f64) -> DuctileJ2<TabulatedDuctileLocus> {
    let native =
        J2Plasticity::new(IsotropicElastic::new(1200.0, 0.2, 0.05).unwrap(), 3.0, 0.0).unwrap();
    let domain = DuctileDomain {
        triaxiality: [-2.0, 2.0],
        strain_rate_s_inv: [1e-10, 10.0],
        temperature_k: [300.0, 300.0],
    };
    let locus = TabulatedDuctileLocus::new(
        vec![(-2.0, INITIATION_STRAIN), (2.0, INITIATION_STRAIN)],
        domain,
    )
    .unwrap();
    DuctileJ2::new(
        native,
        locus,
        DuctileJ2Settings {
            characteristic_length_m: width,
            post_initiation_plastic_work_j_m2: GP,
            maximum_damage_increment: DAMAGE_LIMIT,
            energy_tolerances: EnergyTolerances {
                energy_j_m3: 1e-10,
                yield_pa: 1e-8,
            },
        },
    )
    .unwrap()
}

fn controls() -> DuctileEvolutionControls {
    DuctileEvolutionControls {
        initial_fraction: 1.0,
        minimum_fraction: 1.0 / 8192.0,
        maximum_fraction: 1.0,
        maximum_accepted_steps: 256,
        maximum_trials: 2048,
        maximum_numerical_loss_j: 0.1,
        energy_tolerance_j: 1e-9,
    }
}

fn conditions(duration_s: f64) -> DuctileStepConditions<'static> {
    DuctileStepConditions {
        duration_s,
        temperatures_k: DuctileTemperatures::Uniform(300.0),
    }
}

/// Affine u_x=2*load_factor*y, u_y=u_z=0. Therefore epsilon_xy is exactly
/// load_factor in the native tensor-Voigt convention. All reactions are real
/// native assembled forces; no spring or constitutive mock replaces the tet.
fn prescribed(nodes: &[[f64; 3]]) -> Vec<(usize, f64)> {
    nodes
        .iter()
        .enumerate()
        .flat_map(|(node, p)| {
            [
                (3 * node, 2.0 * p[1]),
                (3 * node + 1, 0.0),
                (3 * node + 2, 0.0),
            ]
        })
        .collect()
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 5e-11 * (1.0 + expected.abs()),
        "actual={actual:e}, expected={expected:e}"
    );
}

fn alpha_at(shear: f64) -> f64 {
    // Independent perfect-plastic radial return: mu=500, sigma_y=3.
    ((1000.0 * 3.0_f64.sqrt() * shear - 3.0) / 1500.0).max(0.0)
}

fn shear_stress(shear: f64, width: f64) -> f64 {
    let alpha = alpha_at(shear);
    let integrity = (-width * 3.0 * (alpha - INITIATION_STRAIN).max(0.0) / GP).exp();
    (1000.0 * shear).min(3.0_f64.sqrt()) * integrity
}

/// Independent endpoint traction work for this exact affine positive-shear
/// path. Gp*(1-g) is a diagnostic; it is deliberately excluded from the sum of
/// storage, physical plastic loss, damage release and numerical loss.
fn audit_affine_work(
    step: &DuctileTetEvolutionStep<'_, TabulatedDuctileLocus>,
    bands: &[TetCrackBandGeometry],
    initial_factor: f64,
) {
    let mut previous = initial_factor;
    let mut external = 0.0;
    for record in &step.records {
        let expected: f64 = bands
            .iter()
            .map(|band| {
                2.0 * band.reference_volume_m3()
                    * shear_stress(record.load_factor, band.band_width_m())
                    * (record.load_factor - previous)
            })
            .sum();
        near(record.energy.external_work_j, expected);
        near(
            expected,
            record.energy.stored_change_j
                + record.energy.plastic_dissipation_j
                + record.energy.damage_release_j
                + record.energy.numerical_loss_j,
        );
        assert!(record.energy.plastic_dissipation_j >= 0.0);
        assert!(record.energy.damage_release_j >= 0.0);
        assert!(record.energy.numerical_loss_j >= 0.0);
        near(record.energy.residual_j, 0.0);
        external += expected;
        previous = record.load_factor;
    }
    near(step.energy.external_work_j, external);
    near(
        external,
        step.energy.stored_change_j
            + step.energy.plastic_dissipation_j
            + step.energy.damage_release_j
            + step.energy.numerical_loss_j,
    );
    near(step.energy.residual_j, 0.0);
    near(
        step.records.iter().map(|r| r.duration_s).sum(),
        step.duration_s,
    );
}

#[test]
fn native_specimen_crosses_initiation_and_loses_force_with_bounded_damage_steps() {
    with_cx(|_, cx| {
        let geometry =
            TetLoadGeometry3::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let bands = [TetCrackBandGeometry::from_element(&geometry, 0, [0.0, 1.0, 0.0]).unwrap()];
        let material = model(bands[0].band_width_m());
        let materials = [&material];
        let supports = prescribed(&NODES);
        let specimen = DuctileTetEvolution::new(
            &geometry,
            &materials,
            &bands,
            &supports,
            &[],
            TetAssemblyBudget::standard(),
            IncrementSettings::default(),
            cx,
        )
        .unwrap();
        let virgin = [material.initial_state()];
        let intact = specimen
            .advance_from(
                &[0.0; 12],
                &virgin,
                [0.0, 0.003],
                conditions(1.0),
                controls(),
                cx,
            )
            .unwrap();
        assert!(!intact.mechanics.element_states[0].initiation().initiated());
        assert_eq!(intact.mechanics.element_states[0].damage(), 0.0);
        near(intact.mechanics.element_stress[0][3], 3.0_f64.sqrt());
        audit_affine_work(&intact, &bands, 0.0);

        let old_u = intact.mechanics.displacement_m.as_flattened().to_vec();
        let old_states = intact.mechanics.element_states.clone();
        let softened = specimen
            .advance_from(
                &old_u,
                &old_states,
                [0.003, 0.012],
                conditions(3.0),
                controls(),
                cx,
            )
            .unwrap();
        assert!(
            softened.rejected_trials > 0,
            "a full damage jump must trigger refinement"
        );
        assert!(softened.records.len() > 1);
        let mut previous_damage = 0.0;
        let mut previous_fraction = 0.0;
        for record in &softened.records {
            assert!(record.maximum_damage >= previous_damage);
            assert!(record.maximum_damage - previous_damage <= DAMAGE_LIMIT + 1e-13);
            assert!(record.fraction > previous_fraction);
            previous_damage = record.maximum_damage;
            previous_fraction = record.fraction;
        }
        assert_eq!(previous_fraction, 1.0);
        let state = &softened.mechanics.element_states[0];
        assert!(state.initiation().initiated());
        assert!(state.damage() > 0.85 && state.integrity() > 0.0);
        near(state.native_state().alpha, alpha_at(0.012));
        // Only node 2 has nonzero prescribed shear displacement. Its actual
        // native reaction is tau*V/height; height is one metre in this fixture.
        let first_force = intact.mechanics.increments.last().unwrap().reactions_n[6];
        let last_force = softened.mechanics.increments.last().unwrap().reactions_n[6];
        near(
            last_force,
            bands[0].reference_volume_m3() * shear_stress(0.012, bands[0].band_width_m()),
        );
        assert!(last_force > 0.0 && last_force < 0.2 * first_force);
        audit_affine_work(&softened, &bands, 0.003);
        let volume = bands[0].reference_volume_m3();
        near(
            intact.energy.external_work_j + softened.energy.external_work_j,
            volume
                * (state.stored_energy_j_m3()
                    + state.dissipated_plastic_j_m3()
                    + state.dissipated_damage_j_m3()
                    + state.numerical_dissipation_j_m3()),
        );
        assert!(softened.energy.damage_release_j > 0.0);
        assert!(softened.energy.degraded_equivalent_plastic_work_j > 0.0);
        assert_eq!(old_states, intact.mechanics.element_states);
        assert_eq!(
            old_u.as_slice(),
            intact.mechanics.displacement_m.as_flattened()
        );
        assert_eq!(virgin[0].damage(), 0.0);
    });
}

#[test]
fn whole_interval_refusal_and_cancellation_preserve_restart_bit_for_bit() {
    with_cx(|gate, cx| {
        let geometry =
            TetLoadGeometry3::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let bands = [TetCrackBandGeometry::from_element(&geometry, 0, [0.0, 1.0, 0.0]).unwrap()];
        let material = model(bands[0].band_width_m());
        let materials = [&material];
        let supports = prescribed(&NODES);
        let specimen = DuctileTetEvolution::new(
            &geometry,
            &materials,
            &bands,
            &supports,
            &[],
            TetAssemblyBudget::standard(),
            IncrementSettings::default(),
            cx,
        )
        .unwrap();
        let prefix = specimen
            .advance_from(
                &[0.0; 12],
                &[material.initial_state()],
                [0.0, 0.003],
                conditions(1.0),
                controls(),
                cx,
            )
            .unwrap();
        let old_u = prefix.mechanics.displacement_m.as_flattened().to_vec();
        let old_states = prefix.mechanics.element_states.clone();
        let saved_u = old_u.clone();
        let saved_states = old_states.clone();
        let bounded = DuctileEvolutionControls {
            maximum_accepted_steps: 1,
            ..controls()
        };
        assert!(matches!(
            specimen.advance_from(
                &old_u,
                &old_states,
                [0.003, 0.012],
                conditions(3.0),
                bounded,
                cx
            ),
            Err(DuctileTetEvolutionError::WorkLimit {
                accepted_steps: 1,
                ..
            })
        ));
        assert_eq!(old_u, saved_u);
        assert_eq!(old_states, saved_states);
        let accepted = specimen
            .advance_from(
                &old_u,
                &old_states,
                [0.003, 0.012],
                conditions(3.0),
                controls(),
                cx,
            )
            .unwrap();

        gate.request();
        assert!(matches!(
            specimen.advance_from(
                &old_u,
                &old_states,
                [0.003, 0.012],
                conditions(3.0),
                controls(),
                cx
            ),
            Err(DuctileTetEvolutionError::Cancelled)
        ));
        assert_eq!(old_u, saved_u);
        assert_eq!(old_states, saved_states);
        let restarted = with_cx(|_, fresh_cx| {
            specimen.advance_from(
                &old_u,
                &old_states,
                [0.003, 0.012],
                conditions(3.0),
                controls(),
                fresh_cx,
            )
        })
        .unwrap();
        assert_eq!(restarted.mechanics, accepted.mechanics);
        assert_eq!(restarted.records, accepted.records);
        assert_eq!(restarted.energy, accepted.energy);
        assert_eq!(restarted.trials, accepted.trials);
        assert_eq!(restarted.rejected_trials, accepted.rejected_trials);
        assert_eq!(old_states, saved_states);
    });
}

#[test]
fn wrong_material_width_and_foreign_geometry_with_the_same_width_refuse() {
    with_cx(|_, cx| {
        let geometry =
            TetLoadGeometry3::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        // An ordinary non-axis-aligned input must survive the geometry owner's
        // normalization when the adapter independently reconstructs this band.
        let normal = [3.0, 4.0, 0.0];
        let bands = [TetCrackBandGeometry::from_element(&geometry, 0, normal).unwrap()];
        let supports = prescribed(&NODES);
        let wrong = model(2.0 * bands[0].band_width_m());
        let wrong_materials = [&wrong];
        assert!(matches!(
            DuctileTetEvolution::new(
                &geometry,
                &wrong_materials,
                &bands,
                &supports,
                &[],
                TetAssemblyBudget::standard(),
                IncrementSettings::default(),
                cx
            ),
            Err(DuctileTetEvolutionError::Invalid(_))
        ));

        // Doubling z doubles volume and projected area for this xy-normal,
        // preserving V/A. Equality of the scalar width cannot admit this mesh.
        let changed_nodes = NODES.map(|p| [p[0], p[1], 2.0 * p[2]]);
        let foreign =
            TetLoadGeometry3::new(&changed_nodes, &TETS, TetAssemblyBudget::standard(), 4, cx)
                .unwrap();
        let foreign_bands = [TetCrackBandGeometry::from_element(&foreign, 0, normal).unwrap()];
        near(foreign_bands[0].band_width_m(), bands[0].band_width_m());
        assert_ne!(
            foreign_bands[0].reference_volume_m3(),
            bands[0].reference_volume_m3()
        );
        let correct = model(bands[0].band_width_m());
        let materials = [&correct];
        assert!(matches!(
            DuctileTetEvolution::new(
                &geometry,
                &materials,
                &foreign_bands,
                &supports,
                &[],
                TetAssemblyBudget::standard(),
                IncrementSettings::default(),
                cx
            ),
            Err(DuctileTetEvolutionError::Invalid(_))
        ));
        assert!(
            DuctileTetEvolution::new(
                &geometry,
                &materials,
                &bands,
                &supports,
                &[],
                TetAssemblyBudget::standard(),
                IncrementSettings::default(),
                cx
            )
            .is_ok()
        );
    });
}

#[test]
fn two_native_elements_keep_distinct_widths_histories_and_projected_area_work() {
    with_cx(|_, cx| {
        let second = NODES.map(|p| [p[0] + 2.0, 0.5 * p[1], p[2]]);
        let nodes: Vec<_> = NODES.into_iter().chain(second).collect();
        let tets = [[0, 1, 2, 3], [4, 5, 6, 7]];
        let geometry =
            TetLoadGeometry3::new(&nodes, &tets, TetAssemblyBudget::standard(), 8, cx).unwrap();
        let bands = [0, 1].map(|element| {
            TetCrackBandGeometry::from_element(&geometry, element, [0.0, 1.0, 0.0]).unwrap()
        });
        let models = bands.map(|band| model(band.band_width_m()));
        let materials = models.each_ref();
        let states = models.each_ref().map(|model| model.initial_state());
        let supports = prescribed(&nodes);
        let specimen = DuctileTetEvolution::new(
            &geometry,
            &materials,
            &bands,
            &supports,
            &[],
            TetAssemblyBudget::standard(),
            IncrementSettings::default(),
            cx,
        )
        .unwrap();
        let result = specimen
            .advance_from(
                &[0.0; 24],
                &states,
                [0.0, 0.012],
                conditions(4.0),
                controls(),
                cx,
            )
            .unwrap();
        assert_eq!(result.mechanics.element_states.len(), 2);
        assert!(result.rejected_trials > 0);
        audit_affine_work(&result, &bands, 0.0);
        for element in 0..2 {
            let state = &result.mechanics.element_states[element];
            assert!(state.initiation().initiated());
            near(state.native_state().alpha, alpha_at(0.012));
            near(
                result.mechanics.element_stress[element][3],
                shear_stress(0.012, bands[element].band_width_m()),
            );
            assert_eq!(states[element].damage(), 0.0);
        }
        let first = &result.mechanics.element_states[0];
        let second = &result.mechanics.element_states[1];
        assert!(first.damage() > second.damage());
        near(second.integrity().powi(2), first.integrity());
        near(bands[0].projected_area_m2(), bands[1].projected_area_m2());
        let expected_area_work: f64 = bands
            .iter()
            .zip(&result.mechanics.element_states)
            .map(|(band, state)| band.projected_area_m2() * GP * state.damage())
            .sum();
        near(
            result.energy.degraded_equivalent_plastic_work_j,
            expected_area_work,
        );
        let expected_storage: f64 = bands
            .iter()
            .zip(&result.mechanics.element_states)
            .map(|(band, state)| band.reference_volume_m3() * state.stored_energy_j_m3())
            .sum();
        near(result.energy.stored_change_j, expected_storage);
        let reactions = &result.mechanics.increments.last().unwrap().reactions_n;
        let force = reactions[6] + reactions[18];
        near(
            force,
            (shear_stress(0.012, bands[0].band_width_m())
                + shear_stress(0.012, bands[1].band_width_m()))
                / 6.0,
        );
    });
}
