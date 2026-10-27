//! Actual source gas fields and source wall memory share one DG time operator.
use super::*;
use fs_couple::acoustic_realize::AcousticRealizeError;
use fs_couple::bernoulli_aperture::wall::material::{
    LocalMaterialWall, RegionalMaterialWallDg, WallRelaxation,
};
use fs_duct::dg::{DgBoundary, DgDuct, DgMesh};
use fs_material::gas_field::{AxialGasField, PrimitiveGasState};
use fs_material::gas_species::input::ThermochemicalPhase;
use fs_material::thermochemical_gas::packs::from_packs;
use fs_thermochem::{Composition, CompositionBasis};

// Reuse the authored normalized species/NASA-9/transport packs, including
// their receipts and explicit synthetic-data scope, without copying a fit.
#[path = "../../src/render/schedule/reed/field_input/tests/support.rs"]
mod gas_support;
#[path = "spatial/live.rs"]
mod live;

fn composition(fraction: f64) -> Composition {
    Composition::new(
        CompositionBasis::MoleFraction,
        vec![fraction, 1.0 - fraction],
    )
    .unwrap()
}
fn gas_field(
    outlet_fraction: f64,
    outlet_temperature: f64,
) -> Result<AxialGasField, fs_material::MaterialError> {
    let a = gas_support::packs("A", 0.028, 3.5, 200.0);
    let b = gas_support::packs("B", 0.044, 4.5, 200.0);
    // Both catalog species are admitted; the primitive inlet below still has
    // exactly zero B. Actual field samples replace this required positive seed.
    let model = from_packs(&[a.binding(), b.binding()], composition(0.5), None)?;
    AxialGasField::new(
        vec![0.0, 0.3],
        101325.0,
        vec![
            PrimitiveGasState {
                temperature_k: 300.0,
                composition: composition(1.0),
            },
            PrimitiveGasState {
                temperature_k: outlet_temperature,
                composition: composition(outlet_fraction),
            },
        ],
        model,
        QueryPoint::new(),
        ThermochemicalPhase::FrozenGas,
    )
}
fn source_wall(density: f64, modulus: f64) -> LocalMaterialWall {
    source_wall_with_arm(density, modulus, 0.8, 0.0005)
}
fn source_wall_with_arm(density: f64, modulus: f64, arm_fraction: f64, tau: f64) -> LocalMaterialWall {
    let mut values = properties(density, modulus, 2.0).to_vec();
    values.push(("wall_arm_modulus", Pressure::DIMS, arm_fraction * modulus));
    values.push(("wall_arm_time", Dims([0, 0, 1, 0, 0, 0]), tau));
    let source = state(&values, None, None);
    let law = WallRelaxation::from_state_point(
        &source,
        &[("wall_arm_modulus".into(), "wall_arm_time".into())],
        BAND,
        1.0,
        0.8,
    )
    .unwrap();
    LocalMaterialWall::new(
        source,
        0.0002,
        AxialConstraint::StressFree,
        WallMaterialLoss::Maxwell(law),
        BAND,
    )
    .unwrap()
}
fn duct(field: &AxialGasField) -> DgDuct {
    DgDuct::new_with_gas_field(
        &fs_duct::Duct {
            segments: vec![fs_duct::Segment::Cylinder {
                radius: 0.006,
                length: 0.3,
            }],
        },
        field,
        DgMesh {
            elements: 2,
            order: 2,
        },
        DgBoundary::Closed,
        DgBoundary::Closed,
    )
    .unwrap()
}

#[test]
fn changing_only_source_composition_changes_real_pressure_and_material_motion() {
    let mut histories = Vec::new();
    let mut common_inlet = None;
    let mut common_wall = None;
    for (fraction, density, modulus) in
        [(0.85, 1000.0, 1e6), (0.2, 1000.0, 1e6), (0.85, 1800.0, 4e6)]
    {
        let field = gas_field(fraction, 300.0).unwrap();
        let inlet = *field.state_at(0.0).unwrap().state();
        if let Some(previous) = common_inlet {
            assert_eq!(previous, inlet);
        }
        common_inlet = Some(inlet);
        assert_eq!(field.model().source_packs().len(), 2);
        let wall = source_wall(density, modulus);
        let identity = wall.material().identity();
        let bound = wall.bind_dg(duct(&field)).unwrap();
        assert_eq!(bound.wall().material().identity(), identity);
        let dg = bound.duct();
        assert_eq!(dg.gas_field(), Some(&field));
        if density == 1000.0 {
            let pin = dg.compliant_wall().unwrap();
            if let Some(previous) = common_wall {
                assert_eq!(pin, previous);
            }
            common_wall = Some(pin);
        }
        let propagation = dg.propagator(2e-6, 1e-13).unwrap();
        let mut current = dg.project_initial_state(|_| 1.0, |_| 0.0);
        let initial = current.stored_energy_j();
        let mut previous = initial;
        let (mut pressure, mut displacement, mut maximum_memory) =
            (Vec::new(), Vec::new(), 0.0_f64);
        for _ in 0..100 {
            propagation.step(&mut current, &[]).unwrap();
            propagation.sync_nodal(&mut current).unwrap();
            let motion = dg.compliant_wall_state(&current).unwrap();
            assert_eq!(motion.relaxation_pressure_pa.len(), 1);
            assert!(current.stored_energy_j() <= previous + 2e-12 * initial);
            previous = current.stored_energy_j();
            close(motion.stored_energy_j, current.wall_stored_energy_j());
            close(
                motion.dissipated_power_w,
                dg.wall_dissipated_power_w(&current).unwrap(),
            );
            maximum_memory = maximum_memory.max(motion.relaxation_stored_energy_j);
            pressure.push(*current.pressure().last().unwrap());
            displacement.push(*motion.displacement_m.last().unwrap());
        }
        assert!(maximum_memory > 0.0);
        assert!(current.stored_energy_j() < initial);
        histories.push((pressure, displacement));
    }
    for changed in &histories[1..] {
        for (actual, reference) in [(&changed.0, &histories[0].0), (&changed.1, &histories[0].1)] {
            let signal: f64 = reference.iter().map(|x| x * x).sum();
            let difference: f64 = actual
                .iter()
                .zip(reference)
                .map(|(a, b)| (a - b) * (a - b))
                .sum();
            assert!(difference > 1e-8 * signal);
        }
    }
}

#[test]
fn source_gaps_and_cancellation_refuse_without_publishing_a_wall_operator() {
    assert!(gas_field(0.2, 1100.0).is_err()); // Original NASA-9/transport domain.
    let a = gas_support::packs("A", 0.028, 3.5, 200.0);
    let b = gas_support::packs("B", 0.044, 4.5, 200.0);
    let mut missing = b.binding();
    missing.transport = None;
    assert!(from_packs(&[a.binding(), missing], composition(0.2), None).is_err());
    let field = gas_field(0.2, 300.0).unwrap();
    let original = duct(&field);
    let identity = original.operator_identity();
    for stop in [1, 20, 40] {
        let mut polls = 0;
        assert!(matches!(
            source_wall(1000.0, 1e6).bind_dg_controlled(original.clone(), &mut || {
                polls += 1;
                polls == stop
            }),
            Err(AcousticRealizeError::Cancelled)
        ));
        assert_eq!(original.operator_identity(), identity);
        assert!(original.compliant_wall().is_none());
        assert_eq!(original.gas_field(), Some(&field));
    }
}

fn kelvin_wall(rho: f64, young: f64, viscosity: f64, thickness: f64) -> LocalMaterialWall {
    LocalMaterialWall::new(
        state(&properties(rho, young, viscosity), None, None),
        thickness,
        AxialConstraint::StrainFixed,
        WallMaterialLoss::KelvinVoigt,
        BAND,
    )
    .unwrap()
}

/// G3: only one material assignment changes; the original regional gas,
/// geometry, initial pressure and numerical solver are held fixed.
#[test]
fn one_regional_material_swap_changes_the_same_gas_and_wall_trajectory() {
    let field = gas_field(0.2, 350.0).unwrap();
    let left = kelvin_wall(1000.0, 1e6, 2.0, 0.0002);
    let mut histories = Vec::new();
    for right in [
        kelvin_wall(1800.0, 4e6, 5.0, 0.00015),
        kelvin_wall(3500.0, 9e6, 10.0, 0.0001),
    ] {
        let bound =
            RegionalMaterialWallDg::new(vec![left.clone(), right.clone()], duct(&field)).unwrap();
        assert_eq!(bound.walls(), &[left.clone(), right.clone()]);
        assert_eq!(bound.band_hz(), BAND);
        let dg = bound.duct();
        assert_eq!(dg.gas_field(), Some(&field));
        assert_eq!(
            dg.compliant_wall_at_element(0),
            Some(left.base_patch(0.006, 0.15).unwrap().wall)
        );
        assert_eq!(
            dg.compliant_wall_at_element(1),
            Some(right.base_patch(0.006, 0.15).unwrap().wall)
        );
        let propagation = dg.propagator(2e-6, 1e-13).unwrap();
        let mut current = dg.project_initial_state(|_| 1.0, |_| 0.0);
        let initial = current.stored_energy_j();
        let mut previous = initial;
        let mut history = Vec::new();
        for _ in 0..100 {
            propagation.step(&mut current, &[]).unwrap();
            propagation.sync_nodal(&mut current).unwrap();
            let wall = dg.compliant_wall_state(&current).unwrap();
            close(wall.stored_energy_j, current.wall_stored_energy_j());
            close(
                wall.dissipated_power_w,
                dg.wall_dissipated_power_w(&current).unwrap(),
            );
            assert!(current.stored_energy_j() <= previous + 2e-12 * initial);
            previous = current.stored_energy_j();
            history.push((
                current.pressure()[0],
                wall.displacement_m[0],
                wall.displacement_m[5],
            ));
        }
        assert!(current.stored_energy_j() < initial);
        assert!(dg.compliant_wall_state(&current).unwrap().stored_energy_j > 0.0);
        histories.push(history);
    }
    // The change reaches the unchanged opposite wall through the physical gas,
    // as well as local wall motion and pressure. It is not isolated metadata.
    for observable in 0..3 {
        let value = |row: &(f64, f64, f64)| match observable {
            0 => row.0,
            1 => row.1,
            _ => row.2,
        };
        let signal: f64 = histories[0].iter().map(|row| value(row).powi(2)).sum();
        let difference: f64 = histories[0]
            .iter()
            .zip(&histories[1])
            .map(|(a, b)| (value(a) - value(b)).powi(2))
            .sum();
        assert!(difference > 1e-8 * signal);
    }
}

#[test]
fn regional_source_admission_keeps_uniform_and_mixed_material_history() {
    let field = gas_field(0.2, 300.0).unwrap();
    let original = duct(&field);
    let memory = source_wall(1000.0, 1e6);
    let uniform = memory.clone().bind_dg(original.clone()).unwrap();
    let regional = RegionalMaterialWallDg::new(vec![memory.clone(); 2], original.clone()).unwrap();
    assert_eq!(
        uniform.duct().operator_identity(),
        regional.duct().operator_identity()
    );
    assert_eq!(regional.walls(), &[memory.clone(), memory.clone()]);
    let elastic = kelvin_wall(1800.0, 4e6, 5.0, 0.00015);
    let mixed = RegionalMaterialWallDg::new(vec![memory, elastic.clone()], original.clone()).unwrap();
    let dg = mixed.duct();
    let mut current = dg.project_initial_state(|_| 1.0, |_| 0.0);
    dg.propagator(2e-6, 1e-13).unwrap().step(&mut current, &[]).unwrap();
    let observed = dg.compliant_wall_state(&current).unwrap();
    let local = dg.mesh().order + 1;
    assert_eq!(observed.relaxation_pressure_pa.len(), 1);
    assert!(observed.relaxation_pressure_pa[0][..local].iter().any(|p| *p != 0.0));
    assert!(observed.relaxation_pressure_pa[0][local..].iter().all(|p| *p == 0.0));
    assert!(RegionalMaterialWallDg::new(vec![elastic.clone()], original.clone()).is_err());
    let other = kelvin_wall(1000.0, 1e6, 2.0, 0.0002);
    for stop in [1, 10, 20] {
        let mut polls = 0;
        assert!(matches!(
            RegionalMaterialWallDg::new_controlled(
                vec![elastic.clone(), other.clone()],
                original.clone(),
                &mut || {
                    polls += 1;
                    polls == stop
                },
            ),
            Err(AcousticRealizeError::Cancelled)
        ));
        assert!(original.compliant_wall_at_element(0).is_none());
    }
}

/// Source-only relaxation-time changes traverse the original material binding,
/// spatial gas, native propagation and pressure observation without a new fit.
#[test]
fn regional_source_spectra_reach_pressure_with_every_local_history_retained() {
    let field = gas_field(0.2, 350.0).unwrap();
    let left = source_wall(1000.0, 1e6);
    let mut histories = Vec::new();
    for tau in [0.0002, 0.002] {
        let right = source_wall_with_arm(1800.0, 4e6, 0.35, tau);
        let source_ids = [left.material().identity(), right.material().identity()];
        let bound = RegionalMaterialWallDg::new(
            vec![left.clone(), right], duct(&field),
        ).unwrap();
        assert_eq!(bound.walls()[0].material().identity(), source_ids[0]);
        assert_eq!(bound.walls()[1].material().identity(), source_ids[1]);
        assert_eq!(bound.band_hz(), BAND);
        let dg = bound.into_duct();
        assert_eq!(dg.gas_field(), Some(&field));
        let local = dg.mesh().order + 1;
        let propagator = dg.propagator(2e-6, 1e-13).unwrap();
        let mut current = dg.project_initial_state(|_| 1.0, |_| 0.0);
        let initial = current.stored_energy_j();
        let mut previous = initial;
        let mut history = Vec::new();
        for _ in 0..100 {
            propagator.step(&mut current, &[]).unwrap();
            propagator.sync_nodal(&mut current).unwrap();
            let wall = dg.compliant_wall_state(&current).unwrap();
            assert_eq!(wall.relaxation_pressure_pa.len(), 2);
            assert!(wall.relaxation_pressure_pa[0][local..].iter().all(|p| *p == 0.0));
            assert!(wall.relaxation_pressure_pa[1][..local].iter().all(|p| *p == 0.0));
            assert!(wall.relaxation_stored_energy_j > 0.0);
            close(wall.stored_energy_j, current.wall_stored_energy_j());
            assert!(current.stored_energy_j() <= previous + 2e-12 * initial);
            previous = current.stored_energy_j();
            history.push(current.pressure()[0]);
        }
        assert!(current.stored_energy_j() < initial);
        histories.push(history);
    }
    let signal: f64 = histories[0].iter().map(|p| p * p).sum();
    let difference: f64 = histories[0].iter().zip(&histories[1]).map(|(a,b)| (a-b).powi(2)).sum();
    assert!(difference > 1e-10 * signal);
}
