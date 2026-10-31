use super::*;
use crate::material_state::gas_scattering::tests::{budget, with_cx};
use fs_geom::{Point3, Vec3};
use fs_material::gas::{GasSpec, GasState};
use fs_material::gas_scattering::{MolecularRayleigh, MolecularRayleighMixture};
use fs_qty::semantic::{Composition, CompositionBasis};

fn model(fractions: &[f64], sources: &[u8]) -> MolecularRayleighMixture {
    let optical: Vec<_> = sources
        .iter()
        .map(|&source| {
            MolecularRayleigh::try_new(0.0003, 2.5e25, 550.0, 0.03, ContentHash([source; 32]))
                .unwrap()
        })
        .collect();
    MolecularRayleighMixture::try_new(
        &vec![GasSpec::dry_air_ussa1976(); fractions.len()],
        Composition::new(CompositionBasis::MoleFraction, fractions.to_vec()).unwrap(),
        &optical,
        None,
    )
    .unwrap()
}

fn density(states: &[RayleighGasState; 2]) -> [f64; 2] {
    states.map(|state| state.gas().density)
}

fn grid(density_kg_m3: &[f64]) -> VolumeGrid<'_> {
    VolumeGrid::new([1, 1, 2], density_kg_m3, [0.0; 3], [1.0, 1.0, 1000.0])
}

fn near(got: f64, expected: f64) {
    assert!(
        (got - expected).abs() < 3.0e-13 * expected.abs().max(1.0e-10),
        "{got:.17e} != {expected:.17e}"
    );
}

#[test]
fn g1_two_layer_gas_uses_original_states_and_exact_spectral_optical_depth() {
    with_cx(|cx, _| {
        let model = model(&[1.0], &[1]);
        let original = [
            GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100_000.0).unwrap(),
            GasState::try_new(&GasSpec::dry_air_ussa1976(), 450.0, 200_000.0).unwrap(),
        ];
        let states = original.map(|state| model.bind_gas_state(&state).unwrap());
        let density = density(&states);
        let binding = MaterialGasScatteringBinding::try_new(states[0]).unwrap();
        let field = binding
            .bind_cell_grid(cx, grid(&density), &states, budget())
            .unwrap();
        assert_eq!(field.states().as_ptr(), states.as_ptr());
        assert_eq!(field.density_kg_m3().as_ptr(), density.as_ptr());
        for (accepted, original) in field.states().iter().zip(&original) {
            assert_eq!(accepted.gas(), original);
        }
        assert!(!field.medium().has_thermal_emission());

        let wavelengths = [380.0, 417.3, 550.0, 709.7, 780.0];
        // The finite segment traverses 750 m of cell 0 and 1000 m of cell 1.
        let forward = field
            .medium()
            .transmittance(
                cx,
                Point3::new(0.5, 0.5, 250.0),
                Vec3::new(0.0, 0.0, 1.0),
                1750.0,
                &wavelengths,
            )
            .unwrap();
        let reverse = field
            .medium()
            .transmittance(
                cx,
                Point3::new(0.5, 0.5, 2000.0),
                Vec3::new(0.0, 0.0, -1.0),
                1750.0,
                &wavelengths,
            )
            .unwrap();
        for ((got, back), wavelength) in forward.into_iter().zip(reverse).zip(wavelengths) {
            let tau = 750.0 * states[0].scattering_per_m(wavelength).unwrap()
                + 1000.0 * states[1].scattering_per_m(wavelength).unwrap();
            near(got, (-tau).exp());
            near(back, got);
        }
        // A distinct reference temperature is immaterial to this accepted grid;
        // the declared fixed model, not the reference snapshot, owns its optics.
        let other_reference =
            MaterialGasScatteringBinding::try_new(model.evaluate(700.0, 300_000.0).unwrap())
                .unwrap();
        let rebound = other_reference
            .bind_cell_grid(cx, grid(&density), &states, budget())
            .unwrap();
        assert_eq!(field.content_identity(), rebound.content_identity());
    });
}

#[test]
fn g3_heterogeneous_pressure_and_temperature_scaling_preserves_state_provenance() {
    with_cx(|cx, _| {
        let model = model(&[1.0], &[1]);
        let states = [
            model.evaluate(300.0, 100_000.0).unwrap(),
            model.evaluate(450.0, 200_000.0).unwrap(),
        ];
        let density = density(&states);
        let binding = MaterialGasScatteringBinding::try_new(states[0]).unwrap();
        let baseline = binding
            .bind_cell_grid(cx, grid(&density), &states, budget())
            .unwrap();
        let wavelengths = [380.0, 550.0, 780.0];
        let transmission = |field: &MaterialGasScatteringCellGrid<'_>| {
            field
                .medium()
                .transmittance(
                    cx,
                    Point3::new(0.5, 0.5, 0.0),
                    Vec3::new(0.0, 0.0, 1.0),
                    2000.0,
                    &wavelengths,
                )
                .unwrap()
        };
        let initial = transmission(&baseline);
        for (temperature_factor, pressure_factor, optical_depth_factor) in
            [(1.0, 2.0, 2.0), (2.0, 1.0, 0.5), (2.0, 2.0, 1.0)]
        {
            let changed = states.map(|state| {
                model
                    .evaluate(
                        temperature_factor * state.gas().temperature,
                        pressure_factor * state.gas().pressure,
                    )
                    .unwrap()
            });
            let changed_density = changed.map(|state| state.gas().density);
            let field = binding
                .bind_cell_grid(cx, grid(&changed_density), &changed, budget())
                .unwrap();
            for (got, original) in transmission(&field).into_iter().zip(initial) {
                near(got, original.powf(optical_depth_factor));
            }
            assert_ne!(baseline.content_identity(), field.content_identity());
            if optical_depth_factor == 1.0 {
                assert_eq!(density, changed_density);
                assert_eq!(
                    baseline.medium().content_identity(),
                    field.medium().content_identity()
                );
            }
        }
    });
}

#[test]
fn g0_cell_state_count_density_layout_geometry_and_budget_must_match() {
    with_cx(|cx, _| {
        let model = model(&[1.0], &[1]);
        let states = [
            model.evaluate(300.0, 100_000.0).unwrap(),
            model.evaluate(450.0, 200_000.0).unwrap(),
        ];
        let density = density(&states);
        let binding = MaterialGasScatteringBinding::try_new(states[0]).unwrap();
        assert!(matches!(
            binding.bind_cell_grid(cx, grid(&density), &states[..1], budget()),
            Err(GasGridBindingError::StateCount)
        ));
        let reversed = [density[1], density[0]];
        assert!(matches!(
            binding.bind_cell_grid(cx, grid(&reversed), &states, budget()),
            Err(GasGridBindingError::DensityMismatch { index: 0 })
        ));
        let wrong_shape = VolumeGrid {
            dims: [usize::MAX; 3],
            data: &density,
            origin: [0.0; 3],
            cell: [1.0; 3],
        };
        assert!(matches!(
            binding.bind_cell_grid(cx, wrong_shape, &states, budget()),
            Err(GasGridBindingError::Transport(VolumeTransportError::Grid))
        ));
        assert!(matches!(
            binding.bind_cell_grid(
                cx,
                grid(&density),
                &states,
                AbsorbingGridBudget {
                    max_grid_cells: 1,
                    ..budget()
                }
            ),
            Err(GasGridBindingError::Transport(VolumeTransportError::Budget))
        ));
    });
}

#[test]
fn g0_equal_coefficients_do_not_flatten_composition_or_optical_source_changes() {
    with_cx(|cx, _| {
        // Identical thermo and optical numbers deliberately hide composition
        // changes numerically. Source-bound model identity must still reject.
        let original = model(&[0.25, 0.75], &[1, 2])
            .evaluate(300.0, 100_000.0)
            .unwrap();
        let composition = model(&[0.75, 0.25], &[1, 2])
            .evaluate(300.0, 100_000.0)
            .unwrap();
        let source = model(&[0.25, 0.75], &[1, 3])
            .evaluate(300.0, 100_000.0)
            .unwrap();
        let binding = MaterialGasScatteringBinding::try_new(original).unwrap();
        for other in [composition, source] {
            assert_eq!(original.gas().density, other.gas().density);
            assert_eq!(
                original.mass_scattering_reference_m2_kg(),
                other.mass_scattering_reference_m2_kg()
            );
            assert_ne!(original.model_identity(), other.model_identity());
            let states = [original, other];
            let density = density(&states);
            assert!(matches!(
                binding.bind_cell_grid(cx, grid(&density), &states, budget()),
                Err(GasGridBindingError::ModelMismatch { index: 1 })
            ));
        }
    });
}

#[test]
fn g4_cancelled_cell_binding_publishes_no_grid_and_preserves_original_buffers() {
    let model = model(&[1.0], &[1]);
    let states = [
        model.evaluate(300.0, 100_000.0).unwrap(),
        model.evaluate(450.0, 200_000.0).unwrap(),
    ];
    let density = density(&states);
    let binding = MaterialGasScatteringBinding::try_new(states[0]).unwrap();
    with_cx(|cx, gate| {
        gate.request();
        assert!(matches!(
            binding.bind_cell_grid(cx, grid(&density), &states, budget()),
            Err(GasGridBindingError::Transport(
                VolumeTransportError::Cancelled
            ))
        ));
    });
    with_cx(|cx, _| {
        let field = binding
            .bind_cell_grid(cx, grid(&density), &states, budget())
            .unwrap();
        assert_eq!(field.states(), &states);
        assert_eq!(field.density_kg_m3(), &density);
    });
}

#[test]
fn g1_composition_grid_uses_each_cells_spectrum_phase_and_reuses_fixed_models() {
    with_cx(|cx, _| {
        let make = |fractions| {
            MolecularRayleighMixture::try_new(
                &[GasSpec::dry_air_ussa1976(); 2],
                Composition::new(CompositionBasis::MoleFraction, fractions).unwrap(),
                &[
                    MolecularRayleigh::try_new(0.00015, 2.5e25, 550.0, 0.0, ContentHash([1; 32]))
                        .unwrap(),
                    MolecularRayleigh::try_new(0.0006, 2.5e25, 550.0, 0.3, ContentHash([2; 32]))
                        .unwrap(),
                ],
                None,
            )
            .unwrap()
        };
        let a = make(vec![0.75, 0.25]);
        let b = make(vec![0.25, 0.75]);
        let states = [
            a.evaluate(300.0, 100_000.0).unwrap(),
            b.evaluate(300.0, 100_000.0).unwrap(),
            a.evaluate(600.0, 100_000.0).unwrap(),
        ];
        let density = states.map(|state| state.gas().density);
        let field = MaterialGasScatteringCellGrid::try_new(
            cx,
            VolumeGrid::new([1, 1, 3], &density, [0.0; 3], [1.0, 1.0, 1000.0]),
            &states,
            CellPaletteBudget {
                max_models: 2,
                max_cells: 3,
            },
            AbsorbingGridBudget {
                max_grid_cells: 3,
                ..budget()
            },
        )
        .unwrap();
        assert_eq!(field.states().as_ptr(), states.as_ptr());
        assert_eq!(field.density_kg_m3().as_ptr(), density.as_ptr());
        let palette = field.medium().scattering_palette().unwrap();
        assert_eq!(palette.values().len(), 2);
        assert_eq!(palette.cell_indices().as_ref(), &[0, 1, 0]);
        for (index, state) in states.iter().enumerate() {
            for cosine in [-1.0, 0.0, 0.6] {
                near(
                    field
                        .medium()
                        .phase_at_cell(index)
                        .unwrap()
                        .eval_cos(cosine)
                        .unwrap(),
                    state.phase_per_sr(cosine).unwrap(),
                );
            }
        }
        assert_ne!(
            states[0].depolarization_ratio(),
            states[1].depolarization_ratio()
        );
        let wavelengths = [380.0, 417.3, 550.0, 780.0];
        let transmission = field
            .medium()
            .transmittance(
                cx,
                Point3::new(0.5, 0.5, 0.0),
                Vec3::new(0.0, 0.0, 1.0),
                3000.0,
                &wavelengths,
            )
            .unwrap();
        for (got, wavelength) in transmission.into_iter().zip(wavelengths) {
            let tau: f64 = states
                .iter()
                .map(|state| 1000.0 * state.scattering_per_m(wavelength).unwrap())
                .sum();
            near(got, (-tau).exp());
        }
    });
}

#[test]
fn g0_composition_grid_enforces_source_palette_budget_density_and_cancellation() {
    with_cx(|cx, gate| {
        let states = [
            model(&[1.0], &[1]).evaluate(300.0, 100_000.0).unwrap(),
            model(&[1.0], &[2]).evaluate(300.0, 100_000.0).unwrap(),
        ];
        let density = density(&states);
        let palette_budget = CellPaletteBudget {
            max_models: 2,
            max_cells: 2,
        };
        // Numerically equal optical data from different declared sources are
        // distinct models, not deduplicated by coefficient coincidence.
        assert!(matches!(
            MaterialGasScatteringCellGrid::try_new(
                cx,
                grid(&density),
                &states,
                CellPaletteBudget {
                    max_models: 1,
                    ..palette_budget
                },
                budget(),
            ),
            Err(GasGridBindingError::Transport(VolumeTransportError::Budget))
        ));
        let field = MaterialGasScatteringCellGrid::try_new(
            cx,
            grid(&density),
            &states,
            palette_budget,
            budget(),
        )
        .unwrap();
        assert_eq!(
            field.medium().scattering_palette().unwrap().values().len(),
            2
        );
        let invalid = [density[0], density[1] * 2.0];
        assert!(matches!(
            MaterialGasScatteringCellGrid::try_new(
                cx,
                grid(&invalid),
                &states,
                palette_budget,
                budget(),
            ),
            Err(GasGridBindingError::DensityMismatch { index: 1 })
        ));
        gate.request();
        assert!(matches!(
            MaterialGasScatteringCellGrid::try_new(
                cx,
                grid(&density),
                &states,
                palette_budget,
                budget(),
            ),
            Err(GasGridBindingError::Transport(
                VolumeTransportError::Cancelled
            ))
        ));
    });
}
