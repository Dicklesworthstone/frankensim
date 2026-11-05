//! G1/G3/G5: original thermal/composition fields in all retained radial modes.
use fs_duct::modal::{
    ModalGasFieldProfile, ModalMouthImpedance, ModalSourceField, ModeCount, MouthLoad,
    PreparedModalProfile,
};
use fs_duct::network::PortExcitation;
use fs_duct::profile::gas_field::{GasFieldDuct, GasFieldDuctError};
use fs_duct::{Duct, DuctError, LossModel, RadiusProfile, RiccatiSteps, Segment, Termination};
use fs_matdb::QueryPoint;
use fs_material::gas_field::{AxialGasField, PrimitiveGasState};
use fs_material::gas_species::input::ThermochemicalPhase;
use fs_material::thermochemical_gas::packs::from_packs;
use fs_math::c64::C64;
use fs_thermochem::{Composition, CompositionBasis};
use std::sync::Arc;

// Original, explicitly synthetic NASA/transport packs, not measured-gas data.
#[path = "../../fs-couple/src/render/schedule/reed/field_input/tests/support.rs"]
mod gas_support;

fn field(
    stations: Vec<f64>,
    temperatures: &[f64],
    fractions: &[f64],
    cellwise: bool,
) -> AxialGasField {
    let a = gas_support::packs("A", 0.028, 3.5, 200.0);
    let b = gas_support::packs("B", 0.044, 4.5, 200.0);
    let composition =
        |f| Composition::new(CompositionBasis::MoleFraction, vec![f, 1.0 - f]).unwrap();
    let model = from_packs(&[a.binding(), b.binding()], composition(0.5), None).unwrap();
    let states = temperatures
        .iter()
        .zip(fractions)
        .map(|(&temperature_k, &fraction)| PrimitiveGasState {
            temperature_k,
            composition: composition(fraction),
        })
        .collect();
    let new = if cellwise {
        AxialGasField::new_piecewise_constant
    } else {
        AxialGasField::new
    };
    new(
        stations,
        101_325.0,
        states,
        model,
        QueryPoint::new(),
        ThermochemicalPhase::FrozenGas,
    )
    .unwrap()
}

fn flows(n: usize) -> Vec<PortExcitation> {
    [
        C64::new(1.0e-5, -2.0e-6),
        C64::new(-2.0e-6, 1.0e-6),
        C64::new(0.5e-6, -1.0e-6),
    ][..n]
        .iter()
        .copied()
        .map(PortExcitation::VolumeFlow)
        .collect()
}

fn exact_complex(a: &[C64], b: &[C64]) {
    assert_eq!(a.len(), b.len());
    for (a, b) in a.iter().zip(b) {
        assert_eq!(a.re.to_bits(), b.re.to_bits());
        assert_eq!(a.im.to_bits(), b.im.to_bits());
    }
}

fn same_fields(a: &ModalSourceField, b: &ModalSourceField) {
    assert_eq!(a.field.x, b.field.x);
    for (a, b) in a.field.pressure.iter().zip(&b.field.pressure) {
        exact_complex(a, b);
    }
    for (a, b) in a.field.flow.iter().zip(&b.field.flow) {
        exact_complex(a, b);
    }
    exact_complex(&a.field.axis_pressure, &b.field.axis_pressure);
    assert_eq!(
        a.inlet.inlet_power_w.to_bits(),
        b.inlet.inlet_power_w.to_bits()
    );
}

fn close(a: C64, b: C64, scale: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * b.abs().max(scale),
        "{a:?} != {b:?}"
    );
}

fn power(p: &[C64], u: &[C64]) -> f64 {
    p.iter().zip(u).map(|(p, u)| 0.5 * (*p * u.conj()).re).sum()
}

#[test]
fn uniform_primitive_field_preserves_original_modal_fields_bitwise() {
    let length = 0.08;
    let profile = RadiusProfile::exponential(0.025, 6.0, length).unwrap();
    let field = field(vec![0.0, 0.03, length], &[300.0; 3], &[0.8; 3], false);
    let binding = ModalGasFieldProfile::new(&profile, &field, RiccatiSteps(32)).unwrap();
    let gas = binding.inlet_gas();
    let omega = 11_000.0;
    assert_eq!(
        binding.steps_used(),
        32,
        "uniform source keeps the original integration axis"
    );
    assert_eq!(binding.profile(), &profile);
    assert_eq!(binding.field(), &field);
    for loss in [LossModel::Lossless, LossModel::WideTube, LossModel::Bessel] {
        for load in [
            MouthLoad::Matched,
            MouthLoad::Analytic(Termination::Closed),
            MouthLoad::Analytic(Termination::IdealOpen),
        ] {
            let expected = PreparedModalProfile::new(
                &profile,
                gas,
                omega,
                loss,
                load,
                ModeCount(3),
                RiccatiSteps(32),
                &mut || false,
            )
            .unwrap();
            let actual = binding
                .prepare(omega, loss, load, ModeCount(3), &mut || false)
                .unwrap();
            assert_eq!(actual.gas_field(), Some(&field));
            assert_eq!(actual.gas_state(), gas);
            assert_eq!(actual.outlet_gas_state(), gas);
            let sources = flows(3);
            same_fields(
                &actual.drive_sources(&sources, &mut || false).unwrap(),
                &expected.drive_sources(&sources, &mut || false).unwrap(),
            );
            let refs = [actual.reference_impedance_pa_s_m3(); 3];
            exact_complex(
                &actual.at_references(&refs).unwrap(),
                &expected.at_references(&refs).unwrap(),
            );
        }
    }
}

#[test]
fn nonuniform_plane_image_matches_independent_webster_owner_and_actual_stations() {
    // G3: two independent propagation implementations, including Bessel loss,
    // a material knot not on the old uniform grid, and a separate spline knot.
    let length = 0.12;
    let profile =
        RadiusProfile::spline(vec![(0.0, 0.023), (0.071, 0.031), (length, 0.028)], None).unwrap();
    let field = field(
        vec![0.0, 0.037, length],
        &[280.0, 330.0, 540.0],
        &[0.8, 0.6, 0.3],
        false,
    );
    let resolution = RiccatiSteps(24);
    let binding = ModalGasFieldProfile::new(&profile, &field, resolution).unwrap();
    let duct = Duct {
        segments: vec![Segment::Profile {
            profile: Arc::new(profile),
            steps: resolution,
        }],
    };
    let reference = GasFieldDuct::new(&duct, &field, resolution).unwrap();
    assert_eq!(binding.steps_used(), 72);
    for loss in [LossModel::Lossless, LossModel::Bessel] {
        for terminal in [
            Termination::Closed,
            Termination::IdealOpen,
            Termination::UnflangedOpen,
        ] {
            let omega = 8_000.0;
            let prepared = binding
                .prepare(
                    omega,
                    loss,
                    MouthLoad::Analytic(terminal),
                    ModeCount(1),
                    &mut || false,
                )
                .unwrap();
            let response = prepared.drive_sources(&flows(1), &mut || false).unwrap();
            let expected = reference.response(omega, loss, terminal, None).unwrap();
            let input_flow = response.field.flow[0][0];
            let actual = response.field.pressure[0][0] / input_flow;
            let zref = prepared.reference_impedance_pa_s_m3();
            close(actual, expected.impedance, zref, 3.0e-10);
            close(
                *response.field.flow.last().unwrap().first().unwrap() / input_flow,
                expected.u_mouth_over_u_in,
                1.0,
                3.0e-10,
            );
            close(
                *response.field.pressure.last().unwrap().first().unwrap() / input_flow,
                expected.p_mouth_over_u_in,
                zref,
                3.0e-10,
            );
            assert_eq!(response.field.x[0], 0.0);
            assert_eq!(*response.field.x.last().unwrap(), length);
            assert!(response.field.x.contains(&0.037));
            assert!(response.field.x.contains(&0.071));
            assert!(response.field.x.windows(2).all(|x| x[1] > x[0]));
        }
    }
}

#[test]
fn higher_modes_retain_full_radiation_matrix_and_lossless_power_in_a_gas_gradient() {
    // G3: a taper couples all three signed modes, with propagating/evanescent
    // changes along the source field. Verify local port work at every station.
    let length = 0.12;
    let profile = RadiusProfile::exponential(0.025, 7.5, length).unwrap();
    let field = field(
        vec![0.0, 0.043, length],
        &[280.0, 370.0, 550.0],
        &[0.8, 0.6, 0.3],
        false,
    );
    let binding = ModalGasFieldProfile::new(&profile, &field, RiccatiSteps(32)).unwrap();
    let omega = 60_000.0;
    let modes = ModeCount(3);
    let exterior = *field.state_at(0.02).unwrap().state();
    let matrix = fs_duct::flanged::flanged_modal_impedance(
        &exterior,
        profile.outlet_radius(),
        omega,
        modes.0,
    )
    .unwrap();
    let load = ModalMouthImpedance::new(
        omega,
        profile.outlet_radius(),
        exterior,
        modes,
        matrix.clone(),
    )
    .unwrap();
    let full = binding
        .prepare_with_flanged_mouth(&exterior, omega, LossModel::Lossless, modes, &mut || false)
        .unwrap();
    let declared = binding
        .prepare_with_mouth_matrix(omega, LossModel::Lossless, &load, modes, &mut || false)
        .unwrap();
    let sources = flows(3);
    let actual = full.drive_sources(&sources, &mut || false).unwrap();
    same_fields(
        &actual,
        &declared.drive_sources(&sources, &mut || false).unwrap(),
    );
    let outgoing = power(
        actual.field.pressure.last().unwrap(),
        actual.field.flow.last().unwrap(),
    );
    assert!(outgoing > 1.0e-7);
    for (p, u) in actual.field.pressure.iter().zip(&actual.field.flow) {
        assert!((power(p, u) - outgoing).abs() < 2.0e-10 * outgoing);
    }
    let p = actual.field.pressure.last().unwrap();
    let u = actual.field.flow.last().unwrap();
    for i in 0..3 {
        let expected: C64 = (0..3).fold(C64::ZERO, |sum, j| sum + matrix[3 * i + j] * u[j]);
        close(p[i], expected, 1.0, 3.0e-11);
    }
    // The local gradient changes physics, rather than just its source metadata.
    let frozen = PreparedModalProfile::with_flanged_mouth(
        &profile,
        binding.inlet_gas(),
        &exterior,
        omega,
        LossModel::Lossless,
        modes,
        RiccatiSteps(64),
        &mut || false,
    )
    .unwrap();
    let frozen = frozen.drive_sources(&sources, &mut || false).unwrap();
    assert!(
        (frozen.field.pressure[0][0] - actual.field.pressure[0][0]).abs()
            > 0.01 * actual.field.pressure[0][0].abs()
    );
}

#[test]
fn continuous_material_and_modal_coupling_converge_on_the_same_magnus_axis() {
    // G1: fourth-order Magnus convergence with material and geometric coupling,
    // including positive Bessel wall loss. The independent mode count stays fixed.
    let profile = RadiusProfile::exponential(0.025, 8.0, 0.12).unwrap();
    let field = field(vec![0.0, 0.12], &[280.0, 540.0], &[0.8, 0.3], false);
    let sources = flows(2);
    let run = |steps| {
        ModalGasFieldProfile::new(&profile, &field, RiccatiSteps(steps))
            .unwrap()
            .prepare(
                30_000.0,
                LossModel::Bessel,
                MouthLoad::Matched,
                ModeCount(2),
                &mut || false,
            )
            .unwrap()
            .drive_sources(&sources, &mut || false)
            .unwrap()
    };
    let reference = run(128);
    let norm = reference.field.pressure[0]
        .iter()
        .map(|p| p.norm_sq())
        .sum::<f64>()
        .sqrt();
    let errors: Vec<_> = [8, 16, 32]
        .into_iter()
        .map(|steps| {
            let actual = run(steps);
            let outgoing = power(
                actual.field.pressure.last().unwrap(),
                actual.field.flow.last().unwrap(),
            );
            assert!(actual.inlet.inlet_power_w > outgoing && outgoing > 0.0);
            actual.field.pressure[0]
                .iter()
                .zip(&reference.field.pressure[0])
                .map(|(&a, &b)| (a - b).norm_sq())
                .sum::<f64>()
                .sqrt()
                / norm
        })
        .collect();
    eprintln!("modal gas-field Magnus errors: {errors:?}");
    assert!(
        errors[0] > 1.0e-7,
        "coarse rung must visibly exercise the axis"
    );
    for pair in errors.windows(2) {
        let order = (pair[0] / pair[1]).log2();
        assert!(
            (order - 4.0).abs() < 0.2,
            "fourth-order source/geometry convergence: {order}"
        );
    }
}

#[test]
fn cellwise_interfaces_remain_explicit_and_refusals_do_not_publish_partial_bindings() {
    let length = 0.1;
    let profile = RadiusProfile::exponential(0.03, 0.0, length).unwrap();
    let field = field(
        vec![0.0, 0.037, length],
        &[300.0, 300.0, 550.0],
        &[0.8, 0.8, 0.3],
        true,
    );
    let binding = ModalGasFieldProfile::new(&profile, &field, RiccatiSteps(1)).unwrap();
    let fine = ModalGasFieldProfile::new(&profile, &field, RiccatiSteps(8)).unwrap();
    let source = flows(2);
    let run = |binding: &ModalGasFieldProfile| {
        binding
            .prepare(
                50_000.0,
                LossModel::Lossless,
                MouthLoad::Matched,
                ModeCount(2),
                &mut || false,
            )
            .unwrap()
            .drive_sources(&source, &mut || false)
            .unwrap()
    };
    let coarse = run(&binding);
    let resolved = run(&fine);
    assert_eq!(coarse.field.x, vec![0.0, 0.037, length]);
    for (&a, &b) in coarse.field.pressure[0]
        .iter()
        .zip(&resolved.field.pressure[0])
    {
        close(a, b, 1.0, 2.0e-10);
    }
    let wrong_length = RadiusProfile::exponential(0.03, 0.0, length * 1.01).unwrap();
    assert!(matches!(
        ModalGasFieldProfile::new(&wrong_length, &field, RiccatiSteps(8)),
        Err(GasFieldDuctError::Duct(DuctError::BadParameter { .. }))
    ));
    assert!(
        ModalGasFieldProfile::new(
            &profile,
            &field,
            RiccatiSteps(fs_duct::profile::MAX_RICCATI_STEPS)
        )
        .is_err()
    );
    let mut calls = 0;
    assert!(matches!(
        ModalGasFieldProfile::new_controlled(&profile, &field, RiccatiSteps(8), &mut || {
            calls += 1;
            calls == 4
        }),
        Err(GasFieldDuctError::Duct(DuctError::Cancelled { .. }))
    ));
    assert!(matches!(
        binding.prepare(
            50_000.0,
            LossModel::Lossless,
            MouthLoad::Matched,
            ModeCount(2),
            &mut || true
        ),
        Err(DuctError::Cancelled { .. })
    ));
    assert!(
        binding
            .prepare(
                0.0,
                LossModel::Lossless,
                MouthLoad::Matched,
                ModeCount(2),
                &mut || false
            )
            .is_err()
    );
    assert!(
        binding
            .prepare(
                50_000.0,
                LossModel::Lossless,
                MouthLoad::Matched,
                ModeCount(0),
                &mut || false
            )
            .is_err()
    );
}
