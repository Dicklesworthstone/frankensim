//! Actual source gas fields and source wall memory share one DG time operator.
use super::*;
use fs_couple::acoustic_realize::AcousticRealizeError;
use fs_couple::bernoulli_aperture::wall::material::{LocalMaterialWall, WallRelaxation};
use fs_duct::dg::{DgBoundary, DgDuct, DgMesh};
use fs_material::gas_field::{AxialGasField, PrimitiveGasState};
use fs_material::gas_species::input::ThermochemicalPhase;
use fs_material::thermochemical_gas::packs::from_packs;
use fs_thermochem::{Composition, CompositionBasis};

// Reuse the authored normalized species/NASA-9/transport packs, including
// their receipts and explicit synthetic-data scope, without copying a fit.
#[path = "../../src/render/schedule/reed/field_input/tests/support.rs"]
mod gas_support;

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
    let mut values = properties(density, modulus, 2.0).to_vec();
    values.push(("wall_arm_modulus", Pressure::DIMS, 0.8 * modulus));
    values.push(("wall_arm_time", Dims([0, 0, 1, 0, 0, 0]), 0.0005));
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
