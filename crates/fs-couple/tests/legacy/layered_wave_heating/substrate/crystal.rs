//! G1/G3: original coherent crystal fields drive real thermal cells.
//! Authored constitutive controls establish wiring, not measured material data.
use super::*;
use fs_couple::layered_wave_heating::substrate::{
    SubstrateWaveHeatingLimits,
    crystal::{PreparedCrystalIlluminationHeating, PreparedCrystalWaveHeating},
};
use fs_em::plane_wave::{
    optical::UniaxialOpticalResponse,
    uniaxial::oriented::{
        OrientedCrystalFieldSide, OrientedUniaxialInterface,
        absorbing::{AbsorbingCrystalIllumination, PreparedAbsorbingCrystalField},
    },
};
use fs_material::state_point::{SpectralComplexIndexKnot, SpectralComplexIndexStatePoint};

const C0: f64 = 299_792_458.0;
const WAVELENGTH_M: f64 = 550e-9;

#[path = "crystal/slab.rs"]
mod slab;
#[path = "crystal/stack.rs"]
mod stack;

fn crystal_card(name: &str, loss: [f64; 2]) -> MaterialCard {
    let base = card(name, 2.0, 0.0);
    let mut claims = base.claims().clone();
    for (axis, n, k) in [("o", 1.4, loss[0]), ("e", 2.1, loss[1])] {
        for (kind, value) in [("n", n), ("k", k)] {
            for knot in 0..2 {
                claims
                    .insert_claim(PropertyClaim {
                        key: PropertyKey::new(format!("{axis}_{kind}_{knot}"), Dims::NONE),
                        value: PropertyValue::Scalar {
                            value,
                            dims: Dims::NONE,
                        },
                        validity: ValidityDomain::unconstrained()
                            .with_quantity(
                                "T",
                                quantity(QuantityKind::AbsoluteTemperature),
                                280.0,
                                360.0,
                            )
                            .with_quantity(
                                "pressure",
                                quantity(QuantityKind::Pressure),
                                90_000.0,
                                200_000.0,
                            ),
                        interpolation: InterpolationPolicy::ConstantWithinValidity,
                        uncertainty: UncertaintyModel::Unstated,
                        observations: vec![],
                        provenance: Provenance {
                            source: "authored coherent crystal heat control".into(),
                            license: "test-only".into(),
                            artifact: None,
                        },
                    })
                    .unwrap();
            }
        }
    }
    MaterialCard::supersede(&base, claims, vec![]).unwrap()
}

fn axes(card: &MaterialCard, pressure: f64) -> [SpectralComplexIndexStatePoint; 2] {
    let mut requirements = vec![
        ScalarPropertyRequirement::try_new(MU_NAME, MU, ScalarAdmissibility::StrictlyPositive)
            .unwrap(),
    ];
    for axis in ["o", "e"] {
        for kind in ["n", "k"] {
            for knot in 0..2 {
                requirements.push(
                    ScalarPropertyRequirement::try_new(
                        format!("{axis}_{kind}_{knot}"),
                        Dims::NONE,
                        ScalarAdmissibility::Finite,
                    )
                    .unwrap(),
                );
            }
        }
    }
    let resolved = resolve_material_state_point(
        card,
        &at(300.0, pressure),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    ["o", "e"].map(|axis| {
        let knots = [400.0, 800.0].map(|wavelength_nm| {
            let knot = usize::from(wavelength_nm == 800.0);
            SpectralComplexIndexKnot {
                wavelength_nm,
                real_index_property: PropertyKey::new(format!("{axis}_n_{knot}"), Dims::NONE),
                extinction_index_property: PropertyKey::new(format!("{axis}_k_{knot}"), Dims::NONE),
            }
        });
        SpectralComplexIndexStatePoint::try_from_resolved(resolved.clone(), &knots).unwrap()
    })
}

fn limits() -> SubstrateWaveHeatingLimits {
    SubstrateWaveHeatingLimits {
        maximum_cells: 128,
        maximum_mode_cell_pairs: 1024,
        maximum_workspace_bytes: 1 << 20,
    }
}

// Independent depth quadrature of the point field. The production path uses
// neither this sampling nor normal-power subtraction to determine cell heat.
fn integrated_density(
    field: &PreparedAbsorbingCrystalField<'_, '_>,
    side: OrientedCrystalFieldSide,
    a: f64,
    b: f64,
) -> f64 {
    let n = 1024;
    let step = (b - a) / n as f64;
    let mut sum = 0.0;
    for i in 0..=n {
        let z = a + step * i as f64;
        let weight = if i == 0 || i == n {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        sum += weight
            * field
                .field(side, [0.0, 0.0, z])
                .unwrap()
                .absorbed_power_density_w_m3;
    }
    AREA * step * sum / 3.0
}

#[test]
fn coherent_crystal_cell_heat_drives_native_temperature_once() {
    let incident = state(&card("crystal incidence", 1.0, 0.0), false, 100_000.0);
    let card = crystal_card("heated crystal", [0.12, 0.42]);
    let [o, e] = axes(&card, 100_000.0);
    let crystal = UniaxialOpticalResponse::try_new(&o, &e, 410e-9).unwrap();
    let interface = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&incident),
        crystal,
        [1.0, 2.0, -0.6],
    )
    .unwrap();
    let field = interface
        .prepare_absorbing_field(
            core::f64::consts::TAU * C0 / WAVELENGTH_M,
            0.6,
            [C64::new(1.4e5, 0.4e5), C64::new(-0.6e5, 1.0e5)],
        )
        .unwrap();
    let width = 2.0 * WAVELENGTH_M;
    let start = 0.1 * WAVELENGTH_M;
    let native = body(&[&card], &[width], &[4], AREA);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    let authored = [1.0, 2.0, 3.0, 4.0];
    let initial = initial(&bound, &authored);
    let original = initial.cells().clone();
    let heat =
        PreparedCrystalWaveHeating::try_new(&bound, &field, start, limits(), || false).unwrap();
    assert!(std::ptr::eq(heat.source(), &field));
    assert!(std::ptr::eq(heat.thermal_body(), &bound));
    let mut independent_modes = 0.0;
    let mut powers = Vec::new();
    for (i, cell) in heat.cells().iter().enumerate() {
        assert_eq!(cell.cell_id, bound.cell_ids()[i]);
        let [a, b] = cell.depth_m;
        let expected = integrated_density(&field, OrientedCrystalFieldSide::Transmitted, a, b);
        near(cell.absorbed_w, expected, 2e-10 * expected);
        independent_modes += integrated_density(&field, OrientedCrystalFieldSide::Ordinary, a, b)
            + integrated_density(&field, OrientedCrystalFieldSide::Extraordinary, a, b);
        powers.push(expected + authored[i]);
    }
    // These actual tilted lossy modes interfere in the original loss tensor.
    assert!((heat.power().absorbed_w - independent_modes).abs() > 1e-4 * heat.power().absorbed_w);
    assert!(heat.power().residual_w.abs() < 3e-11 * heat.power().entering_w);
    let dt = 2e-9;
    let mut workspace = LayeredThermalWorkspace::new(4).unwrap();
    let step = heat
        .advance_to(
            &initial,
            dt,
            conduction(),
            LayeredThermalSolveControls::default(),
            &mut workspace,
            || false,
        )
        .unwrap();
    let temperatures = thermal_reference(&powers, width, AREA, dt);
    for (i, cell) in heat.cells().iter().enumerate() {
        near(
            step.thermal.temperature_k(cell.cell_id).unwrap(),
            temperatures[i],
            3e-8,
        );
        near(step.cell_energy_j[i].1, cell.absorbed_w * dt, 1e-14);
    }
    near(
        step.heat_report.energy.external_heat_j,
        step.optical.absorbed_j + authored.iter().sum::<f64>() * dt,
        1e-9,
    );
    near(
        step.heat_report.energy.stored_energy_change_j,
        step.heat_report.energy.external_heat_j,
        2e-9,
    );
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.cells(), &original);
}

#[test]
fn crystal_refinement_preserves_total_heat_and_lossless_light_is_not_heat() {
    let incident = state(
        &card("crystal refinement incidence", 1.0, 0.0),
        false,
        100_000.0,
    );
    for loss in [[0.12, 0.42], [0.0, 0.0]] {
        let card = crystal_card("crystal refinement", loss);
        let [o, e] = axes(&card, 100_000.0);
        let crystal = UniaxialOpticalResponse::try_new(&o, &e, 410e-9).unwrap();
        let interface = OrientedUniaxialInterface::try_new(
            DielectricResponse::Lossless(&incident),
            crystal,
            [1.0, 2.0, -0.6],
        )
        .unwrap();
        let field = interface
            .prepare_absorbing_field(
                core::f64::consts::TAU * C0 / WAVELENGTH_M,
                0.6,
                [C64::ONE, C64::new(0.3, 0.7)],
            )
            .unwrap();
        let mut totals = Vec::new();
        for count in [1, 4, 16] {
            let native = body(&[&card], &[WAVELENGTH_M], &[count], AREA);
            let laws = native.thermal_laws(false).unwrap();
            let bound = native.bind(&laws, ids()).unwrap();
            let heat = PreparedCrystalWaveHeating::try_new(&bound, &field, 0.0, limits(), || false)
                .unwrap();
            assert_eq!(heat.cells().len(), count);
            assert_eq!(heat.preparation_mode_cell_pairs(), 4 * (count + 2));
            assert!(heat.cells().iter().all(|c| c.absorbed_w >= 0.0));
            totals.push(heat.power().absorbed_w);
            if loss == [0.0, 0.0] {
                assert_eq!(heat.power().absorbed_w, 0.0);
                assert!(heat.power().entering_w > 0.0);
                near(heat.power().leaving_w, heat.power().entering_w, 1e-15);
                let initial = initial(&bound, &vec![0.0; count]);
                let step = heat
                    .advance_to(
                        &initial,
                        1e-9,
                        conduction(),
                        LayeredThermalSolveControls::default(),
                        &mut LayeredThermalWorkspace::new(count).unwrap(),
                        || false,
                    )
                    .unwrap();
                assert_eq!(step.optical.absorbed_j, 0.0);
                for cell in heat.cells() {
                    assert_eq!(step.thermal.temperature_k(cell.cell_id).unwrap(), 300.0);
                }
            }
        }
        for &total in &totals[1..] {
            near(total, totals[0], 1e-12 * totals[0]);
        }
    }
}

#[test]
fn crystal_heat_refuses_changed_material_context_budget_and_cancellation() {
    let incident = state(
        &card("crystal refusal incidence", 1.0, 0.0),
        false,
        100_000.0,
    );
    let card = crystal_card("crystal refusal", [0.12, 0.42]);
    let [o, e] = axes(&card, 100_000.0);
    let law = UniaxialOpticalResponse::try_new(&o, &e, 410e-9).unwrap();
    let interface = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&incident),
        law,
        [1.0, 2.0, -0.6],
    )
    .unwrap();
    let field = interface
        .prepare_absorbing_field(
            core::f64::consts::TAU * C0 / WAVELENGTH_M,
            0.6,
            [C64::ONE, C64::ZERO],
        )
        .unwrap();
    let native = body(&[&card], &[WAVELENGTH_M], &[3], AREA);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    let heat =
        PreparedCrystalWaveHeating::try_new(&bound, &field, 0.0, limits(), || false).unwrap();
    let exact = SubstrateWaveHeatingLimits {
        maximum_cells: 3,
        maximum_mode_cell_pairs: heat.preparation_mode_cell_pairs(),
        maximum_workspace_bytes: heat.workspace_bytes(),
    };
    assert!(PreparedCrystalWaveHeating::try_new(&bound, &field, 0.0, exact, || false).is_ok());
    for which in 0..3 {
        let mut short = exact;
        match which {
            0 => short.maximum_cells -= 1,
            1 => short.maximum_mode_cell_pairs -= 1,
            _ => short.maximum_workspace_bytes -= 1,
        }
        assert!(PreparedCrystalWaveHeating::try_new(&bound, &field, 0.0, short, || false).is_err());
    }
    for start in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            PreparedCrystalWaveHeating::try_new(&bound, &field, start, limits(), || false).is_err()
        );
    }
    assert!(matches!(
        PreparedCrystalWaveHeating::try_new(&bound, &field, 0.0, limits(), || true),
        Err(LayeredWaveHeatingError::Cancelled)
    ));
    let other = crystal_card("different crystal", [0.12, 0.42]);
    let wrong_native = body(&[&other], &[WAVELENGTH_M], &[3], AREA);
    let wrong_laws = wrong_native.thermal_laws(false).unwrap();
    let wrong = wrong_native.bind(&wrong_laws, ids()).unwrap();
    assert!(PreparedCrystalWaveHeating::try_new(&wrong, &field, 0.0, limits(), || false).is_err());
    let [po, pe] = axes(&card, 150_000.0);
    let wrong_context = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&incident),
        UniaxialOpticalResponse::try_new(&po, &pe, 410e-9).unwrap(),
        [1.0, 2.0, -0.6],
    )
    .unwrap();
    let wrong_field = wrong_context
        .prepare_absorbing_field(
            core::f64::consts::TAU * C0 / WAVELENGTH_M,
            0.6,
            [C64::ONE, C64::ZERO],
        )
        .unwrap();
    assert!(
        PreparedCrystalWaveHeating::try_new(&bound, &wrong_field, 0.0, limits(), || false).is_err()
    );
    let initial = initial(&bound, &[0.0; 3]);
    let original = initial.cells().clone();
    let mut workspace = LayeredThermalWorkspace::new(3).unwrap();
    assert!(matches!(
        heat.advance_to(
            &initial,
            1e-9,
            conduction(),
            LayeredThermalSolveControls::default(),
            &mut workspace,
            || true
        ),
        Err(LayeredWaveHeatingError::Cancelled)
    ));
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.cells(), &original);
}

#[test]
fn crystal_illumination_spectral_heat_drives_one_native_thermal_solve() {
    let incident = state(
        &card("spectral crystal incidence", 1.0, 0.0),
        false,
        100_000.0,
    );
    let card = crystal_card("spectrally heated crystal", [0.12, 0.42]);
    let [o, e] = axes(&card, 100_000.0);
    let interface = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&incident),
        UniaxialOpticalResponse::try_new(&o, &e, 410e-9).unwrap(),
        [1.0, 2.0, -0.6],
    )
    .unwrap();
    let factors = [
        interface
            .prepare_absorbing_field(
                core::f64::consts::TAU * C0 / 450e-9,
                0.6,
                [C64::new(8e4, 2e4), C64::new(-3e4, 7e4)],
            )
            .unwrap(),
        interface
            .prepare_absorbing_field_in_plane(
                core::f64::consts::TAU * C0 / 700e-9,
                0.9,
                [0.6, 0.8],
                [C64::new(2e4, -4e4), C64::new(6e4, 3e4)],
            )
            .unwrap(),
    ];
    let illumination =
        AbsorbingCrystalIllumination::try_new(&interface, &factors, 2, || false).unwrap();
    let width = 1.4e-6;
    let start = 0.15e-6;
    let native = body(&[&card], &[width], &[3], AREA);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    let authored = [1.0, 2.0, 3.0];
    let initial = initial(&bound, &authored);
    let original = initial.cells().clone();
    let heat =
        PreparedCrystalIlluminationHeating::try_new(&bound, &illumination, start, limits(), || {
            false
        })
        .unwrap();
    assert!(std::ptr::eq(heat.source(), &illumination));
    assert!(std::ptr::eq(heat.thermal_body(), &bound));
    assert_eq!(heat.preparation_mode_cell_pairs(), 4 * 2 * (3 + 2));
    let mut powers = Vec::new();
    for (i, cell) in heat.cells().iter().enumerate() {
        assert_eq!(cell.cell_id, bound.cell_ids()[i]);
        let [a, b] = cell.depth_m;
        let expected: f64 = factors
            .iter()
            .map(|field| integrated_density(field, OrientedCrystalFieldSide::Transmitted, a, b))
            .sum();
        near(cell.absorbed_w, expected, 3e-10 * expected);
        powers.push(expected + authored[i]);
    }
    assert!(heat.power().residual_w.abs() < 4e-11 * heat.power().entering_w);
    let dt = 2e-9;
    let step = heat
        .advance_to(
            &initial,
            dt,
            conduction(),
            LayeredThermalSolveControls::default(),
            &mut LayeredThermalWorkspace::new(3).unwrap(),
            || false,
        )
        .unwrap();
    let temperatures = thermal_reference(&powers, width, AREA, dt);
    for (i, cell) in heat.cells().iter().enumerate() {
        near(
            step.thermal.temperature_k(cell.cell_id).unwrap(),
            temperatures[i],
            3e-8,
        );
        near(step.cell_energy_j[i].1, cell.absorbed_w * dt, 1e-14);
    }
    near(
        step.heat_report.energy.external_heat_j,
        step.optical.absorbed_j + authored.iter().sum::<f64>() * dt,
        1e-9,
    );
    near(
        step.heat_report.energy.stored_energy_change_j,
        step.heat_report.energy.external_heat_j,
        2e-9,
    );
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.cells(), &original);
}

#[test]
fn crystal_illumination_unpolarized_heat_is_basis_independent() {
    let incident = state(
        &card("unpolarized crystal incidence", 1.0, 0.0),
        false,
        100_000.0,
    );
    let card = crystal_card("unpolarized heated crystal", [0.12, 0.42]);
    let [o, e] = axes(&card, 100_000.0);
    let interface = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&incident),
        UniaxialOpticalResponse::try_new(&o, &e, 410e-9).unwrap(),
        [1.0, 2.0, -0.6],
    )
    .unwrap();
    let omega = core::f64::consts::TAU * C0 / WAVELENGTH_M;
    // At normal incidence these equal physical transverse amplitudes form
    // an unpolarized coherency matrix. Rotate its orthonormal factor basis.
    let amplitude = 8e4;
    let diagonal = amplitude * core::f64::consts::FRAC_1_SQRT_2;
    let aligned = [
        interface
            .prepare_absorbing_field(omega, 1.0, [C64::new(amplitude, 0.0), C64::ZERO])
            .unwrap(),
        interface
            .prepare_absorbing_field(omega, 1.0, [C64::ZERO, C64::new(amplitude, 0.0)])
            .unwrap(),
    ];
    let rotated = [
        interface
            .prepare_absorbing_field(
                omega,
                1.0,
                [C64::new(diagonal, 0.0), C64::new(diagonal, 0.0)],
            )
            .unwrap(),
        interface
            .prepare_absorbing_field(
                omega,
                1.0,
                [C64::new(-diagonal, 0.0), C64::new(diagonal, 0.0)],
            )
            .unwrap(),
    ];
    let a = AbsorbingCrystalIllumination::try_new(&interface, &aligned, 2, || false).unwrap();
    let b = AbsorbingCrystalIllumination::try_new(&interface, &rotated, 2, || false).unwrap();
    let native = body(&[&card], &[1e-6], &[4], AREA);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    let a_heat =
        PreparedCrystalIlluminationHeating::try_new(&bound, &a, 0.0, limits(), || false).unwrap();
    let b_heat =
        PreparedCrystalIlluminationHeating::try_new(&bound, &b, 0.0, limits(), || false).unwrap();
    for (a_cell, b_cell) in a_heat.cells().iter().zip(b_heat.cells()) {
        near(
            a_cell.absorbed_w,
            b_cell.absorbed_w,
            2e-12 * a_cell.absorbed_w,
        );
    }
    near(
        a_heat.power().absorbed_w,
        b_heat.power().absorbed_w,
        2e-12 * a_heat.power().absorbed_w,
    );
}

#[test]
fn crystal_illumination_darkness_budgets_and_typed_cancellation() {
    let incident = state(
        &card("illumination admission incidence", 1.0, 0.0),
        false,
        100_000.0,
    );
    let card = crystal_card("illumination admission crystal", [0.12, 0.42]);
    let [o, e] = axes(&card, 100_000.0);
    let interface = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&incident),
        UniaxialOpticalResponse::try_new(&o, &e, 410e-9).unwrap(),
        [1.0, 2.0, -0.6],
    )
    .unwrap();
    let factors = [interface
        .prepare_absorbing_field(
            core::f64::consts::TAU * C0 / WAVELENGTH_M,
            0.6,
            [C64::ONE, C64::ZERO],
        )
        .unwrap()];
    let illumination =
        AbsorbingCrystalIllumination::try_new(&interface, &factors, 1, || false).unwrap();
    let native = body(&[&card], &[WAVELENGTH_M], &[3], AREA);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    let mut polls = 0;
    let heat =
        PreparedCrystalIlluminationHeating::try_new(&bound, &illumination, 0.0, limits(), || {
            polls += 1;
            false
        })
        .unwrap();
    // Every possible interruption, including inner factor queries and late
    // boundary evaluation, remains the typed public cancellation outcome.
    for stop in 1..=polls {
        let mut calls = 0;
        assert!(matches!(
            PreparedCrystalIlluminationHeating::try_new(
                &bound,
                &illumination,
                0.0,
                limits(),
                || {
                    calls += 1;
                    calls == stop
                },
            ),
            Err(LayeredWaveHeatingError::Cancelled)
        ));
    }
    let exact = SubstrateWaveHeatingLimits {
        maximum_cells: 3,
        maximum_mode_cell_pairs: heat.preparation_mode_cell_pairs(),
        maximum_workspace_bytes: heat.workspace_bytes(),
    };
    assert!(
        PreparedCrystalIlluminationHeating::try_new(&bound, &illumination, 0.0, exact, || false,)
            .is_ok()
    );
    let short = SubstrateWaveHeatingLimits {
        maximum_mode_cell_pairs: exact.maximum_mode_cell_pairs - 1,
        ..exact
    };
    assert!(
        PreparedCrystalIlluminationHeating::try_new(&bound, &illumination, 0.0, short, || false,)
            .is_err()
    );
    let dark = AbsorbingCrystalIllumination::try_new(&interface, &[], 0, || false).unwrap();
    let dark_limits = SubstrateWaveHeatingLimits {
        maximum_mode_cell_pairs: 0,
        ..limits()
    };
    let dark_heat =
        PreparedCrystalIlluminationHeating::try_new(&bound, &dark, 0.0, dark_limits, || false)
            .unwrap();
    assert_eq!(dark_heat.preparation_mode_cell_pairs(), 0);
    assert_eq!(dark_heat.power().entering_w, 0.0);
    assert_eq!(dark_heat.power().leaving_w, 0.0);
    assert_eq!(dark_heat.power().absorbed_w, 0.0);
    let authored = [1.0, 2.0, 3.0];
    let initial = initial(&bound, &authored);
    let original = initial.cells().clone();
    let dt = 2e-9;
    let step = dark_heat
        .advance_to(
            &initial,
            dt,
            conduction(),
            LayeredThermalSolveControls::default(),
            &mut LayeredThermalWorkspace::new(3).unwrap(),
            || false,
        )
        .unwrap();
    assert_eq!(step.optical.absorbed_j, 0.0);
    near(
        step.heat_report.energy.external_heat_j,
        authored.iter().sum::<f64>() * dt,
        1e-18,
    );
    assert!(matches!(
        heat.advance_to(
            &initial,
            dt,
            conduction(),
            LayeredThermalSolveControls::default(),
            &mut LayeredThermalWorkspace::new(3).unwrap(),
            || true,
        ),
        Err(LayeredWaveHeatingError::Cancelled)
    ));
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.cells(), &original);
}

#[path = "crystal/coherent.rs"]
mod coherent;
