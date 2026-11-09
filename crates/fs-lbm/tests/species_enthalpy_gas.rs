//! Original NASA9 source fixtures drive the native constituent/enthalpy owner.
//! These are independent analytic coupling checks, not measured-gas validation.
#[path = "nsf3_thermochemical/source.rs"]
mod source;

use fs_exec::CancelGate;
use fs_lbm::conjugate::{
    ChtError, DensityProjectionControls, FlowFace, FlowField, GasSpeciesEnthalpySource,
    MaterialAdvectionControls, MaterialAdvectionInput, MaterialAdvectionStep, MaterialVolumeSource,
    SpeciesDiffusionControls, Voxel, VoxelDomain,
};
use fs_material::thermochemical_gas::{ThermochemicalGasModel, energy::GasEnergySolve};
use fs_qty::{
    chemistry::SpeciesId,
    semantic::{Composition, CompositionBasis},
};

const R: f64 = 8.314_462_618_153_24;
const MOLAR_MASS: [f64; 2] = [0.028, 0.018_015_28];
const ENERGY_BUDGET: f64 = 1e-5;

fn energy() -> GasEnergySolve {
    GasEnergySolve {
        minimum_temperature_k: 200.0,
        maximum_temperature_k: 700.0,
        energy_tolerance_j_kg: 1e-7,
        maximum_iterations: 128,
    }
}
fn controls() -> SpeciesDiffusionControls<2> {
    SpeciesDiffusionControls {
        substeps: 1,
        maximum_pair_updates: 1 << 20,
        maximum_storage_bytes: 1 << 22,
        maximum_mass_defect_kg: [1e-11; 2],
        maximum_partition_defect_kg: 1e-11,
    }
}
fn ids() -> [SpeciesId; 2] {
    ["manufactured", "water"].map(|name| SpeciesId::new(name).unwrap())
}
fn carrier(domain: &VoxelDomain, mass: &[f64]) -> MaterialAdvectionStep {
    let faces = [FlowFace::Wall; 6];
    FlowField::from_face_velocity(domain, faces, |_, _| 0.0)
        .unwrap()
        .project_and_advect_material(
            domain,
            faces,
            MaterialAdvectionInput {
                mass_kg: mass,
                specific_enthalpy_j_kg: &vec![0.0; mass.len()],
                volume_sources: &vec![MaterialVolumeSource::None; mass.len()],
                boundary_inflow: [None; 6],
            },
            MaterialAdvectionControls {
                projection: DensityProjectionControls {
                    duration_s: 1.0,
                    relative_tolerance: 1e-12,
                    maximum_volume_defect_m3_s: 1e-10,
                    maximum_energy_defect_j: 1e-9,
                },
                max_cells: mass.len(),
                max_connections: 6 * mass.len(),
                max_iterations: 100,
                relative_tolerance: 1e-12,
                maximum_mass_defect_kg: 1e-11,
                maximum_enthalpy_defect_j: 1e-9,
            },
            &CancelGate::new_clock_free(),
        )
        .unwrap()
}
fn enthalpy(bank: &ThermochemicalGasModel, mass: [f64; 2], t: f64) -> f64 {
    let total = mass.iter().sum::<f64>();
    bank.with_composition(
        Composition::new(
            CompositionBasis::MassFraction,
            mass.map(|m| m / total).to_vec(),
        )
        .unwrap(),
    )
    .unwrap()
    .evaluate(&source::point(t, 100_000.0))
    .unwrap()
    .specific_enthalpy_j_kg()
}
fn partial_h(t: f64, a7: f64) -> [f64; 2] {
    MOLAR_MASS.map(|m| R / m * (3.0 * t + 0.0005 * t * t + a7))
}
fn capacity(mass: [f64; 2]) -> f64 {
    (0..2).map(|k| mass[k] * R / MOLAR_MASS[k]).sum()
}
fn temperature(total_h: f64, mass: [f64; 2], a7: f64) -> f64 {
    let z = total_h / capacity(mass) - a7;
    // Stable positive root of the original cp/R = 3 + 0.001*T integral.
    2.0 * z / (3.0 + (9.0 + 0.002 * z).sqrt())
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

#[test]
fn g1_variable_temperature_pair_closes_original_enthalpy_after_composition_changes() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let mass = [1.0, 2.0];
    let accepted = carrier(&domain, &mass);
    let initial = [[0.75, 0.25], [0.25, 1.75]];
    let initial_t = [320.0, 480.0];
    let diffusivity = [0.1, 0.4];
    let names = ids();
    let context = source::point(400.0, 100_000.0);
    let gate = CancelGate::new_clock_free();
    let conductance = 1.0 / (0.5 / (mass[0] * diffusivity[0]) + 0.5 / (mass[1] * diffusivity[1]));
    let harmonic_mass = 1.0 / (1.0 / mass[0] + 1.0 / mass[1]);
    let q = harmonic_mass * (0.75 - 0.125) * (-(-conductance * 0.2 / harmonic_mass).exp_m1());
    let mut reference: Option<[f64; 2]> = None;
    for a7 in [0.0, -4000.0, 10000.0] {
        let bank = source::moist_model(0.5, a7);
        let calorics =
            GasSpeciesEnthalpySource::new(&bank, &names, &context, &diffusivity, energy(), &gate)
                .unwrap();
        let h: [f64; 2] = core::array::from_fn(|c| enthalpy(&bank, initial[c], initial_t[c]));
        let before = accepted.clone();
        let result = accepted
            .diffuse_thermochemical_species(
                &domain,
                &initial,
                &h,
                0.2,
                &calorics,
                controls(),
                ENERGY_BUDGET,
                &gate,
            )
            .unwrap();
        let face_h = partial_h(400.0, a7);
        let exchanged = q * (face_h[0] - face_h[1]);
        for (c, sign) in [(0, -1.0), (1, 1.0)] {
            let expected_mass = [initial[c][0] + sign * q, initial[c][1] - sign * q];
            for k in 0..2 {
                near(
                    result.step.species.species_mass_kg[c][k],
                    expected_mass[k],
                    5e-14,
                );
            }
            near(
                result.step.diffusive_enthalpy_change_j[c],
                sign * exchanged,
                1e-7,
            );
            let expected_t = temperature(mass[c] * h[c] + sign * exchanged, expected_mass, a7);
            near(result.state.temperatures_k()[c], expected_t, 1e-8);
            assert!((expected_t - initial_t[c]).abs() > 0.1);
            let cp = capacity(expected_mass) / mass[c] * (3.0 + 0.001 * expected_t);
            near(
                result.state.transport(c).unwrap().specific_heat_j_kg_k(),
                cp,
                2e-8,
            );
        }
        let temperatures = [
            result.state.temperatures_k()[0],
            result.state.temperatures_k()[1],
        ];
        if let Some(previous) = reference {
            for c in 0..2 {
                near(temperatures[c], previous[c], 1e-8);
            }
        } else {
            reference = Some(temperatures);
        }
        assert!(result.report.absolute_enthalpy_defect_j <= ENERGY_BUDGET);
        assert_eq!(result.state.species(), &names);
        assert_eq!(result.state.source().source_packs(), bank.source_packs());
        assert_eq!(accepted, before);
    }
}

#[test]
fn g0_graded_face_uses_original_nasa_partials_at_the_physical_face_temperature() {
    let domain = VoxelDomain::graded([vec![0.5, 1.5], vec![1.0], vec![1.0]]).unwrap();
    let accepted = carrier(&domain, &[1.0; 2]);
    let masses = [[0.75, 0.25], [0.25, 0.75]];
    let bank = source::moist_model(0.5, -4000.0);
    let h = [
        enthalpy(&bank, masses[0], 320.0),
        enthalpy(&bank, masses[1], 480.0),
    ];
    let names = ids();
    let context = source::point(400.0, 100_000.0);
    let diffusivity = [0.1; 2];
    let gate = CancelGate::new_clock_free();
    let calorics =
        GasSpeciesEnthalpySource::new(&bank, &names, &context, &diffusivity, energy(), &gate)
            .unwrap();
    let field = calorics
        .prepare(&domain, &accepted, &masses, &h, &gate)
        .unwrap();
    let partials = calorics
        .face_partials(&field, &domain, [0, 1], 0, &gate)
        .unwrap();
    // Distances from the face are 0.25 and 0.75 m; arithmetic averaging is wrong.
    let expected = partial_h((1.5 * 320.0 + 0.5 * 480.0) / 2.0, -4000.0);
    for k in 0..2 {
        near(partials[k], expected[k], 2e-7);
    }
    assert!(field.absolute_caloric_residual_j() <= 2.0 * energy().energy_tolerance_j_kg);
}

#[test]
fn g3_exact_absence_preserves_isothermal_source_identity_and_refuses_atomically() {
    let mut domain = VoxelDomain::new(3, 1, 1, 1.0).unwrap();
    domain.set(2, 0, 0, Voxel::Solid(0));
    let accepted = carrier(&domain, &[1.0, 1.0, 0.0]);
    let masses = [[1.0, 0.0], [0.0, 1.0], [0.0; 2]];
    let bank = source::moist_model(0.5, 0.0);
    let h = [
        enthalpy(&bank, masses[0], 400.0),
        enthalpy(&bank, masses[1], 400.0),
        0.0,
    ];
    let names = ids();
    let context = source::point(400.0, 100_000.0);
    let diffusivity = [0.2, 0.2, 0.0];
    let gate = CancelGate::new_clock_free();
    let calorics =
        GasSpeciesEnthalpySource::new(&bank, &names, &context, &diffusivity, energy(), &gate)
            .unwrap();
    let before = (accepted.clone(), masses, h);
    let result = accepted
        .diffuse_thermochemical_species(
            &domain,
            &masses,
            &h,
            0.1,
            &calorics,
            controls(),
            ENERGY_BUDGET,
            &gate,
        )
        .unwrap();
    for c in 0..2 {
        assert!(
            result.step.species.species_mass_kg[c]
                .iter()
                .all(|&m| m > 0.0)
        );
        near(result.state.temperatures_k()[c], 400.0, 1e-8);
        assert!(result.state.transport(c).is_some());
    }
    assert_eq!(result.step.species.species_mass_kg[2], [0.0; 2]);
    assert!(result.state.transport(2).is_none());
    assert_eq!(result.state.source().source_packs(), bank.source_packs());
    assert_eq!(result.state.species(), &names);
    let reversed = [names[1].clone(), names[0].clone()];
    assert!(
        GasSpeciesEnthalpySource::new(&bank, &reversed, &context, &diffusivity, energy(), &gate)
            .is_err()
    );
    let mut unsupported = h;
    unsupported[0] = 1e50;
    assert!(
        calorics
            .prepare(&domain, &accepted, &masses, &unsupported, &gate)
            .is_err()
    );
    let mut wrong_partition = masses;
    wrong_partition[0][0] += 0.01;
    assert!(
        calorics
            .prepare(&domain, &accepted, &wrong_partition, &h, &gate)
            .is_err()
    );
    let insufficient = SpeciesDiffusionControls {
        maximum_storage_bytes: 1,
        ..controls()
    };
    assert!(
        accepted
            .diffuse_thermochemical_species(
                &domain,
                &masses,
                &h,
                0.1,
                &calorics,
                insufficient,
                ENERGY_BUDGET,
                &gate
            )
            .is_err()
    );
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        accepted.diffuse_thermochemical_species(
            &domain,
            &masses,
            &h,
            0.1,
            &calorics,
            controls(),
            ENERGY_BUDGET,
            &cancelled
        ),
        Err(ChtError::Cancelled)
    ));
    assert_eq!(accepted, before.0);
    assert_eq!(masses, before.1);
    assert_eq!(h, before.2);
}


#[path = "species_enthalpy_gas/adaptive.rs"]
mod adaptive;


#[path = "species_enthalpy_gas/boundary.rs"]
mod boundary;
