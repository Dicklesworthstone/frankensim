use super::*;
use crate::material_state::particle_scattering::MaterialParticleScattering;
use crate::volumes::transport::scattering::mie::{MiePopulationBudget, ParticleSizeDistribution};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::gas::GasSpec;
use fs_material::gas_scattering::{MolecularRayleigh, MolecularRayleighMixture};
use fs_material::state_point::*;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, QuantitySpec};
use std::sync::Arc;

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 87,
                    kernel_id: 5,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
            &gate,
        )
    });
}

fn gas(source: u8) -> MolecularRayleighMixture {
    MolecularRayleighMixture::try_new(
        &[GasSpec::dry_air_ussa1976()],
        Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
        &[
            MolecularRayleigh::try_new(0.0003, 2.5e25, 550.0, 0.03, ContentHash([source; 32]))
                .unwrap(),
        ],
        None,
    )
    .unwrap()
}

// Actual source-card admission and Mie preparation, not fabricated spectra or
// replacement material physics. These are declared synthetic fixture sources.
fn particles(cx: &Cx<'_>, variant: usize, temperature: f64) -> MaterialParticleSpectrum {
    let temperature_quantity = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let point = QueryPoint::new()
        .with_quantity("thermal", temperature_quantity, temperature)
        .unwrap();
    let mut values = vec![(
        DENSITY_PROPERTY,
        Density::DIMS,
        1000.0 * (variant + 1) as f64,
    )];
    for (index, wavelength) in VISIBLE_COMPLEX_IOR_WAVELENGTHS_NM.into_iter().enumerate() {
        values.push((
            VISIBLE_COMPLEX_IOR_ETA_PROPERTIES[index],
            Dims::NONE,
            1.4 + 0.3 * variant as f64 - 0.0001 * (wavelength - 380.0),
        ));
        values.push((
            VISIBLE_COMPLEX_IOR_K_PROPERTIES[index],
            Dims::NONE,
            0.02 + 0.05 * variant as f64 + 0.00001 * (wavelength - 380.0),
        ));
    }
    let mut claims = ClaimSet::new();
    for (name, dims, value) in values {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "thermal",
                    temperature_quantity,
                    temperature,
                    temperature,
                ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: Vec::new(),
                provenance: Provenance {
                    source: format!("synthetic coexisting population {variant}; not measured"),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: format!("synthetic optical fixture {variant}"),
            phase: "solid".into(),
            process: "declared sphere population".into(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .unwrap();
    let accepted = resolve_material_state_point(
        &card,
        &point,
        &[ScalarPropertyRequirement::try_new(
            DENSITY_PROPERTY,
            Density::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap()],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let optical = resolve_visible_optical_state_point(
        &card,
        &point,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let binding = MaterialParticleScattering::try_new(
        &accepted,
        &optical,
        &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
        ParticleSizeDistribution::try_new(&[(40e-9 + 60e-9 * variant as f64, 1.0)]).unwrap(),
        1.0,
    )
    .unwrap();
    let wavelengths: &[f64] = if variant == 0 {
        &[380.0, 550.0, 780.0]
    } else {
        &[380.0, 500.0, 650.0, 780.0]
    };
    binding
        .visible_spectrum(cx, wavelengths, MiePopulationBudget::default(), 0.003, 2048)
        .unwrap()
}

fn composition(fractions: &[f64]) -> Composition {
    Composition::new(CompositionBasis::MassFraction, fractions.to_vec()).unwrap()
}

fn mixture_budget() -> OpticalMixtureBudget {
    OpticalMixtureBudget { max_components: 3 }
}

fn grid_budget() -> AbsorbingGridBudget {
    AbsorbingGridBudget {
        max_grid_cells: 2,
        max_steps_per_segment: 8,
    }
}

fn palette<'a>(
    cx: &Cx<'_>,
    mixtures: Vec<MaterialOpticalMixture<'a>>,
    indices: &[u32],
) -> CellPalette<MaterialOpticalMixture<'a>> {
    let budget = CellPaletteBudget {
        max_models: mixtures.len(),
        max_cells: indices.len(),
    };
    CellPalette::try_new(cx, mixtures, Arc::from(indices), budget).unwrap()
}

fn grid(density: &[f64]) -> VolumeGrid<'_> {
    VolumeGrid::new([1, 1, density.len()], density, [0.0; 3], [1.0, 1.0, 1000.0])
}

fn total_density(
    gas: &RayleighGasState,
    composition: &Composition,
    a: &MaterialParticleSpectrum,
    b: &MaterialParticleSpectrum,
) -> f64 {
    let w = composition.fractions();
    1.0 / (w[0] / gas.gas().density
        + w[1] / a.particle_density_kg_m3()
        + w[2] / b.particle_density_kg_m3())
}

fn close(got: f64, expected: f64) {
    assert!(
        (got - expected).abs() <= 4e-13 * expected.abs().max(1e-300),
        "{got:e} != {expected:e}"
    );
}

#[test]
fn g1_gas_and_two_real_particle_bindings_keep_exact_spectra_and_weighted_phases() {
    with_cx(|cx, _| {
        let state = gas(1).evaluate(300.0, 100_000.0).unwrap();
        let carrier = MaterialGasScatteringBinding::try_new(state).unwrap();
        let a = particles(cx, 0, 300.0);
        let b = particles(cx, 1, 300.0);
        let components = [
            MaterialOpticalComponent::Gas(&carrier),
            MaterialOpticalComponent::Particles(&a),
            MaterialOpticalComponent::Particles(&b),
        ];
        let fractions = composition(&[1.0 - 3e-8, 1e-8, 2e-8]);
        let recipe =
            MaterialOpticalMixture::try_new(cx, &components, &fractions, mixture_budget()).unwrap();
        assert_eq!(recipe.components().as_ptr(), components.as_ptr());
        assert!(core::ptr::eq(recipe.composition(), &fractions));
        let density = [total_density(&state, &fractions, &a, &b)];
        let states = [state];
        let field = MaterialOpticalMixtureGrid::try_new(
            cx,
            grid(&density),
            palette(cx, vec![recipe], &[0]),
            Some(&states),
            None,
            1e-6,
            grid_budget(),
        )
        .unwrap();
        assert!(!field.medium().has_thermal_emission());
        assert_eq!(field.density_kg_m3().as_ptr(), density.as_ptr());
        assert_eq!(field.gas_states().unwrap().as_ptr(), states.as_ptr());
        let wavelengths = [380.0, 417.3, 550.0, 779.0];
        let transmission = field
            .medium()
            .transmittance(
                cx,
                Point3::new(0.5, 0.5, -1.0),
                Vec3::new(0.0, 0.0, 1.0),
                1002.0,
                &wavelengths,
            )
            .unwrap();
        for (got, wavelength) in transmission.into_iter().zip(wavelengths) {
            let weights = fractions.fractions();
            let scatter = [
                weights[0] * state.mass_scattering_m2_kg(wavelength).unwrap(),
                weights[1] * a.spectrum().scattering.at(wavelength).unwrap(),
                weights[2] * b.spectrum().scattering.at(wavelength).unwrap(),
            ];
            let absorption = weights[1] * a.spectrum().absorption.at(wavelength).unwrap()
                + weights[2] * b.spectrum().absorption.at(wavelength).unwrap();
            let total: f64 = scatter.iter().sum();
            close(got, (-1000.0 * density[0] * (absorption + total)).exp());
            for cosine in [-1.0, -0.2, 0.0, 0.8, 1.0] {
                let phase = [
                    state.phase_per_sr(cosine).unwrap(),
                    a.spectrum().phase.eval_at(wavelength, cosine).unwrap(),
                    b.spectrum().phase.eval_at(wavelength, cosine).unwrap(),
                ];
                let expected = scatter.iter().zip(phase).map(|(s, p)| s * p).sum::<f64>() / total;
                close(
                    field
                        .medium()
                        .phase_at_cell(0)
                        .unwrap()
                        .eval_at(wavelength, cosine)
                        .unwrap(),
                    expected,
                );
            }
        }
        close(
            field.maximum_volume_fraction(),
            density[0] * (1e-8 / a.particle_density_kg_m3() + 2e-8 / b.particle_density_kg_m3()),
        );
    });
}

#[test]
fn g1_layered_mixture_recipes_share_maps_and_integrate_only_true_lte_absorption() {
    with_cx(|cx, _| {
        let model = gas(1);
        let states = [
            model.evaluate(300.0, 100_000.0).unwrap(),
            model.evaluate(300.0, 50_000.0).unwrap(),
        ];
        let carrier = MaterialGasScatteringBinding::try_new(states[0]).unwrap();
        let a = particles(cx, 0, 300.0);
        let b = particles(cx, 1, 300.0);
        let components = [
            MaterialOpticalComponent::Gas(&carrier),
            MaterialOpticalComponent::Particles(&a),
            MaterialOpticalComponent::Particles(&b),
        ];
        let fractions = [
            composition(&[1.0 - 3e-8, 1e-8, 2e-8]),
            composition(&[1.0 - 6e-8, 4e-8, 2e-8]),
        ];
        let recipes = fractions
            .iter()
            .map(|fractions| {
                MaterialOpticalMixture::try_new(cx, &components, fractions, mixture_budget())
                    .unwrap()
            })
            .collect();
        let original = palette(cx, recipes, &[0, 1]);
        let density = [0, 1].map(|i| total_density(&states[i], &fractions[i], &a, &b));
        let temperatures = [300.0; 2];
        let field = MaterialOpticalMixtureGrid::try_new(
            cx,
            grid(&density),
            original.clone(),
            Some(&states),
            Some(&temperatures),
            1e-6,
            grid_budget(),
        )
        .unwrap();
        assert!(field.medium().has_thermal_emission());
        for map in [
            field.mixtures().cell_indices(),
            field
                .medium()
                .absorption_grid()
                .absorption_palette()
                .unwrap()
                .cell_indices(),
            field.medium().scattering_palette().unwrap().cell_indices(),
        ] {
            assert!(Arc::ptr_eq(map, original.cell_indices()));
        }
        let wavelength = 760.0;
        let result = field
            .medium()
            .transport(
                cx,
                Point3::new(0.5, 0.5, -1.0),
                Vec3::new(0.0, 0.0, 1.0),
                2002.0,
                &[wavelength],
            )
            .unwrap();
        let (mut tau, mut emissivity) = ([0.0; 2], [0.0; 2]);
        for i in 0..2 {
            let f = fractions[i].fractions();
            let absorption = f[1] * a.spectrum().absorption.at(wavelength).unwrap()
                + f[2] * b.spectrum().absorption.at(wavelength).unwrap();
            let scattering = f[0] * states[i].mass_scattering_m2_kg(wavelength).unwrap()
                + f[1] * a.spectrum().scattering.at(wavelength).unwrap()
                + f[2] * b.spectrum().scattering.at(wavelength).unwrap();
            tau[i] = 1000.0 * density[i] * (absorption + scattering);
            emissivity[i] = absorption / (absorption + scattering);
        }
        let planck = crate::thermal_surface::planck_radiance_w_m2_sr_nm(300.0, wavelength).unwrap();
        close(result.transmittance[0], (-(tau[0] + tau[1])).exp());
        close(
            result.emission_w_m2_sr_nm[0],
            planck
                * (emissivity[0] * -(-tau[0]).exp_m1()
                    + (-tau[0]).exp() * emissivity[1] * -(-tau[1]).exp_m1()),
        );
    });
}

#[test]
fn g0_gas_source_and_partial_density_stay_bound_with_only_roundoff_admitted() {
    with_cx(|cx, _| {
        let state = gas(1).evaluate(300.0, 100_000.0).unwrap();
        let carrier = MaterialGasScatteringBinding::try_new(state).unwrap();
        let other =
            MaterialGasScatteringBinding::try_new(gas(2).evaluate(300.0, 100_000.0).unwrap())
                .unwrap();
        let components = [MaterialOpticalComponent::Gas(&carrier)];
        let other_components = [MaterialOpticalComponent::Gas(&other)];
        let fractions = composition(&[1.0]);
        let recipe =
            MaterialOpticalMixture::try_new(cx, &components, &fractions, mixture_budget()).unwrap();
        let other_recipe =
            MaterialOpticalMixture::try_new(cx, &other_components, &fractions, mixture_budget())
                .unwrap();
        assert_ne!(recipe.source_identity(), other_recipe.source_identity());
        let recipes = palette(cx, vec![recipe], &[0]);
        let states = [state];
        let roundoff = [fs_math::next_up(state.gas().density)];
        MaterialOpticalMixtureGrid::try_new(
            cx,
            grid(&roundoff),
            recipes.clone(),
            Some(&states),
            None,
            1e-6,
            grid_budget(),
        )
        .unwrap();
        let invalid = [state.gas().density * (1.0 + 1e-8)];
        assert!(matches!(
            MaterialOpticalMixtureGrid::try_new(
                cx,
                grid(&invalid),
                recipes.clone(),
                Some(&states),
                None,
                1e-6,
                grid_budget()
            ),
            Err(MaterialOpticalMixtureError::GasDensityMismatch { cell: 0 })
        ));
        assert!(matches!(
            MaterialOpticalMixtureGrid::try_new(
                cx,
                grid(&roundoff),
                recipes.clone(),
                None,
                None,
                1e-6,
                grid_budget()
            ),
            Err(MaterialOpticalMixtureError::MissingGasStates)
        ));
        let wrong_source = [*other.state()];
        assert!(matches!(
            MaterialOpticalMixtureGrid::try_new(
                cx,
                grid(&roundoff),
                recipes,
                Some(&wrong_source),
                None,
                1e-6,
                grid_budget()
            ),
            Err(MaterialOpticalMixtureError::GasModelMismatch { cell: 0 })
        ));
        let wrong_basis = Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap();
        assert!(matches!(
            MaterialOpticalMixture::try_new(cx, &components, &wrong_basis, mixture_budget()),
            Err(MaterialOpticalMixtureError::Composition)
        ));
    });
}

#[test]
fn g0_common_lte_and_sum_of_particle_volumes_are_required_before_publication() {
    with_cx(|cx, gate| {
        let state = gas(1).evaluate(300.0, 100_000.0).unwrap();
        let carrier = MaterialGasScatteringBinding::try_new(state).unwrap();
        let a = particles(cx, 0, 300.0);
        let b = particles(cx, 0, 400.0);
        let components = [
            MaterialOpticalComponent::Gas(&carrier),
            MaterialOpticalComponent::Particles(&a),
            MaterialOpticalComponent::Particles(&b),
        ];
        let fractions = composition(&[0.8, 0.1, 0.1]);
        assert!(matches!(
            MaterialOpticalMixture::try_new(
                cx,
                &components,
                &fractions,
                OpticalMixtureBudget { max_components: 2 }
            ),
            Err(MaterialOpticalMixtureError::Transport(
                VolumeTransportError::Budget
            ))
        ));
        let recipe =
            MaterialOpticalMixture::try_new(cx, &components, &fractions, mixture_budget()).unwrap();
        let recipes = palette(cx, vec![recipe], &[0]);
        let states = [state];
        let density = [total_density(&state, &fractions, &a, &b)];
        // Each population is individually below this ceiling; their SUM is not.
        let one_fraction = density[0] * 0.1 / a.particle_density_kg_m3();
        assert!(matches!(
            MaterialOpticalMixtureGrid::try_new(
                cx,
                grid(&density),
                recipes.clone(),
                Some(&states),
                None,
                1.5 * one_fraction,
                grid_budget()
            ),
            Err(MaterialOpticalMixtureError::VolumeFraction { cell: 0 })
        ));
        assert!(matches!(
            MaterialOpticalMixtureGrid::try_new(
                cx,
                grid(&density),
                recipes.clone(),
                Some(&states),
                Some(&[300.0]),
                3.0 * one_fraction,
                grid_budget()
            ),
            Err(MaterialOpticalMixtureError::TemperatureMismatch {
                cell: 0,
                component: 2
            })
        ));
        let nonemitting = MaterialOpticalMixtureGrid::try_new(
            cx,
            grid(&density),
            recipes.clone(),
            Some(&states),
            None,
            3.0 * one_fraction,
            grid_budget(),
        )
        .unwrap();
        assert!(!nonemitting.medium().has_thermal_emission());
        let phi = nonemitting.maximum_volume_fraction();
        close(phi, 2.0 * one_fraction);
        close(phi + density[0] * 0.8 / state.gas().density, 1.0);
        close(
            density[0],
            state.gas().density * (1.0 - phi)
                + one_fraction * (a.particle_density_kg_m3() + b.particle_density_kg_m3()),
        );
        // Filling particle-occupied volume with gas overcounts gas mass. Even
        // at this dilute loading, the discrepancy is far beyond roundoff.
        let displaced_gas_ignored = [state.gas().density / 0.8];
        assert!(matches!(
            MaterialOpticalMixtureGrid::try_new(
                cx,
                grid(&displaced_gas_ignored),
                recipes.clone(),
                Some(&states),
                None,
                3.0 * one_fraction,
                grid_budget(),
            ),
            Err(MaterialOpticalMixtureError::GasDensityMismatch { cell: 0 })
        ));
        let no_gas_volume = [5000.0];
        assert!(matches!(
            MaterialOpticalMixtureGrid::try_new(
                cx,
                grid(&no_gas_volume),
                recipes.clone(),
                Some(&states),
                None,
                1.0,
                grid_budget(),
            ),
            Err(MaterialOpticalMixtureError::GasFreeVolume { cell: 0 })
        ));
        gate.request();
        assert!(matches!(
            MaterialOpticalMixture::try_new(cx, &components, &fractions, mixture_budget()),
            Err(MaterialOpticalMixtureError::Transport(
                VolumeTransportError::Cancelled
            ))
        ));
        assert!(matches!(
            MaterialOpticalMixtureGrid::try_new(
                cx,
                grid(&density),
                recipes,
                Some(&states),
                None,
                3.0 * one_fraction,
                grid_budget()
            ),
            Err(MaterialOpticalMixtureError::Transport(
                VolumeTransportError::Cancelled
            ))
        ));
    });
}
