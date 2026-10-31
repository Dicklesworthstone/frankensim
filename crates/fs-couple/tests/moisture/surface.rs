//! Authored surface-emission data exercise the original moisture/gas/mechanics.
use super::*;
use fs_couple::moisture::atmosphere::{AtmosphereMoistureError, DryBasisSorption};
use fs_couple::moisture::surface::*;

fn emission(values: [f64; 2], upper: f64) -> DryBasisSurfaceTransfer {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(DRY_BASIS_SURFACE_EMISSION_PROPERTY, SURFACE_EMISSION_DIMS),
            value: PropertyValue::Curve {
                abscissa: "relative_humidity".into(),
                abscissa_dims: Dims::NONE,
                knots: vec![(0.0, values[0]), (1.0, values[1])],
                dims: SURFACE_EMISSION_DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "temperature",
                    qty(QuantityKind::AbsoluteTemperature),
                    293.15,
                    293.15,
                )
                .with_quantity(
                    "relative_humidity",
                    qty(QuantityKind::Moisture(MoistureBasis::RelativeHumidity)),
                    0.0,
                    upper,
                ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "authored dry-basis surface-emission verification law".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "authored moisture surface".into(),
            phase: "solid".into(),
            process: "verification".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let point = QueryPoint::new()
        .with_quantity(
            "temperature",
            qty(QuantityKind::AbsoluteTemperature),
            293.15,
        )
        .unwrap()
        .with_quantity(
            "relative_humidity",
            qty(QuantityKind::Moisture(MoistureBasis::RelativeHumidity)),
            0.5,
        )
        .unwrap();
    DryBasisSurfaceTransfer::new(card, point, MaterialPropertySelection::SingleClaimOnly).unwrap()
}

pub(super) fn law(h: [f64; 2]) -> DryBasisSorption {
    atmosphere::sorption(0.25)
        .with_surface_transfer(emission(h, 1.0))
        .unwrap()
}

fn source_faces(law: DryBasisSorption, rh: f64) -> [MoistureFace; 2] {
    let face = MoistureFace::SourceSorption {
        law,
        relative_humidity: rh,
    };
    [face.clone(), face]
}

#[test]
fn g1_surface_resistance_matches_flux_and_actual_surface_with_sealed_and_instant_limits() {
    with_cx(|cx| {
        let initial = slab(card([0.01; 2], 0.3), 1.0, vec![0.1], cx);
        let dt = 2.0;
        let h = 0.004;
        let diffusion_speed = 0.02; // D / (L/2), independent half-cell resistance.
        let effective = 1.0 / (1.0 / h + 1.0 / diffusion_speed);
        let alpha = 2.0 * dt * effective;
        for rh in [0.8, 0.2] {
            let equilibrium = 0.25 * rh;
            let expected = (0.1 + alpha * equilibrium) / (1.0 + alpha);
            let expected_surface =
                (diffusion_speed * expected + h * equilibrium) / (diffusion_speed + h);
            let (next, report) = initial
                .advance(dt, &source_faces(law([h; 2]), rh), config(), cx)
                .unwrap();
            near(next.moisture()[0], expected, 1e-14);
            for side in 0..2 {
                near(
                    report.surface_moisture[side].unwrap(),
                    expected_surface,
                    1e-14,
                );
                near(report.reservoir_moisture[side].unwrap(), equilibrium, 1e-14);
                let flux_integral =
                    initial.dry_mass_kg() * h * dt * (equilibrium - expected_surface);
                near(report.boundary_water_kg[side], flux_integral, 1e-14);
                assert!(report.surface_sorption[side].is_some());
                assert!(report.surface_transfer[side].is_some());
            }
            near(report.water_balance_residual_kg, 0.0, 1e-14);
            assert_eq!(report.boundary_water_kg[0].is_sign_positive(), rh > 0.4);
            assert!(
                (next.moisture()[0] - 0.1).abs()
                    < (initial
                        .advance(dt, &faces(equilibrium), config(), cx)
                        .unwrap()
                        .0
                        .moisture()[0]
                        - 0.1)
                        .abs()
            );
        }
        let nonuniform = slab(card([0.01; 2], 0.3), 1.0, vec![0.1, 0.15, 0.05], cx);
        let zero = nonuniform
            .advance(dt, &source_faces(law([0.0; 2]), 0.8), config(), cx)
            .unwrap();
        let sealed = nonuniform
            .advance(
                dt,
                &[MoistureFace::Sealed, MoistureFace::Sealed],
                config(),
                cx,
            )
            .unwrap();
        assert_eq!(zero.0.moisture(), sealed.0.moisture());
        assert_eq!(zero.1.boundary_water_kg, [0.0; 2]);
        assert_eq!(
            zero.1.surface_moisture,
            [Some(zero.0.moisture()[0]), Some(zero.0.moisture()[2])]
        );
        let instant = initial
            .advance(
                dt,
                &source_faces(atmosphere::sorption(0.25), 0.8),
                config(),
                cx,
            )
            .unwrap();
        let original = initial
            .advance(
                dt,
                &faces(instant.1.reservoir_moisture[0].unwrap()),
                config(),
                cx,
            )
            .unwrap();
        assert_eq!(instant.0.moisture(), original.0.moisture());
        assert_eq!(instant.1.boundary_water_kg, original.1.boundary_water_kg);
        assert!(instant.1.surface_transfer.iter().all(Option::is_none));
        let fast = initial
            .advance(dt, &source_faces(law([1e12; 2]), 0.8), config(), cx)
            .unwrap();
        near(fast.0.moisture()[0], original.0.moisture()[0], 1e-14);
    });
}

#[test]
fn g1_g3_endpoint_rh_updates_surface_source_and_finite_vapor_with_exact_one_cell_balance() {
    with_cx(|cx| {
        let boundary = [Some(law([0.002, 0.008])), Some(law([0.002, 0.008]))];
        let dt = 2.0;
        for initial_rh in [0.8, 0.2] {
            let mut gas = atmosphere::atmosphere(initial_rh, 1.0);
            let mut moisture = slab(card([0.01; 2], 0.3), 1.0, vec![0.1], cx);
            for _ in 0..2 {
                // Eliminating M from the one-cell backward-Euler and finite
                // vapor balances gives a quadratic in RH for h=h0+h1*RH.
                let r0 = gas.relative_humidity().unwrap();
                let m0 = moisture.moisture()[0];
                let saturation = fs_material::gas::saturation_pressure_water_pa(293.15).unwrap()
                    * 0.018_015_28
                    * gas.gas().volume_m3()
                    / (fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * 293.15);
                let dry = moisture.dry_mass_kg();
                let (h0, h1, a, kd) = (0.002, 0.006, 0.25, 0.02);
                let c = 2.0 * dt * kd;
                let qa = saturation * (1.0 + c) * h1 + dry * c * h1 * a;
                let qb = saturation * (kd + (1.0 + c) * h0 - (1.0 + c) * h1 * r0)
                    + dry * c * (h0 * a - h1 * m0);
                let qc = -saturation * r0 * (kd + (1.0 + c) * h0) - dry * c * h0 * m0;
                let expected_rh = -2.0 * qc / (qb + (qb * qb - 4.0 * qa * qc).sqrt());
                let h = h0 + h1 * expected_rh;
                let alpha = c * h / (kd + h);
                let expected_m = (m0 + alpha * a * expected_rh) / (1.0 + alpha);
                let candidate = gas
                    .advance(
                        &moisture,
                        dt,
                        &boundary,
                        config(),
                        atmosphere::controls(),
                        cx,
                    )
                    .unwrap();
                near(
                    candidate.atmosphere().relative_humidity().unwrap(),
                    expected_rh,
                    2e-11,
                );
                near(candidate.moisture().moisture()[0], expected_m, 2e-12);
                near(
                    candidate.atmosphere().gas().mass_kg() + candidate.moisture().wet_mass_kg(),
                    gas.gas().mass_kg() + moisture.wet_mass_kg(),
                    2e-13,
                );
                let selected = candidate.step().transport.surface_transfer[0]
                    .as_ref()
                    .unwrap();
                near(
                    selected
                        .property(DRY_BASIS_SURFACE_EMISSION_PROPERTY)
                        .unwrap()
                        .value_si(),
                    h,
                    2e-13,
                );
                near(
                    selected.query().axes()["relative_humidity"],
                    expected_rh,
                    2e-11,
                );
                assert_ne!(
                    selected.query().axes()["relative_humidity"].to_bits(),
                    r0.to_bits()
                );
                near(
                    candidate.step().maximum_carrier_mass_residual_kg,
                    0.0,
                    2e-13,
                );
                let replay = gas
                    .advance(
                        &moisture,
                        dt,
                        &boundary,
                        config(),
                        atmosphere::controls(),
                        cx,
                    )
                    .unwrap();
                assert_eq!(candidate.atmosphere(), replay.atmosphere());
                assert_eq!(
                    candidate.moisture().moisture(),
                    replay.moisture().moisture()
                );
                assert_eq!(
                    candidate.step().transport.surface_transfer,
                    replay.step().transport.surface_transfer
                );
                (gas, moisture, _) = candidate.into_parts();
            }
        }
    });
}

#[test]
fn g3_surface_data_changes_finite_uptake_and_original_native_plate_response() {
    with_cx(|cx| {
        let initial = slab(card([0.01; 2], 0.3), 1.0, vec![0.1; 4], cx);
        let binding = plate(&initial, PlateMaterialModel::Isotropic, cx);
        let gas = atmosphere::atmosphere(0.8, 2.0);
        let mut masses = Vec::new();
        let mut eigenvalues = Vec::new();
        for source in [law([0.001; 2]), law([0.02; 2]), atmosphere::sorption(0.25)] {
            let boundary = [Some(source.clone()), Some(source)];
            let (candidate, mechanics) = binding
                .advance_finite_atmosphere(
                    &gas,
                    &initial,
                    2.0,
                    &boundary,
                    config(),
                    atmosphere::controls(),
                    cx,
                )
                .unwrap();
            near(
                mechanics.wet_mass_kg(),
                candidate.moisture().wet_mass_kg(),
                1e-13,
            );
            assert_eq!(mechanics.time_s(), candidate.moisture().time_s());
            for (source, &m) in mechanics
                .elastic_states()
                .iter()
                .zip(candidate.moisture().moisture())
            {
                near(
                    source.property("young_modulus").unwrap().value_si(),
                    1e9 * (1.0 - m),
                    1e-6,
                );
            }
            let modal = modes(mechanics.model(), (15.0, 90.0), &SliceOptions::default()).unwrap();
            assert_eq!(modal.modes.len(), 1);
            masses.push(mechanics.wet_mass_kg());
            eigenvalues.push(modal.modes[0].lambda);
        }
        assert!(masses[0] < masses[1] && masses[1] < masses[2]);
        assert!(eigenvalues[0] > eigenvalues[1] && eigenvalues[1] > eigenvalues[2]);
    });
}

#[test]
fn g4_surface_source_late_mechanical_and_cancel_refusals_preserve_accepted_states() {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 3,
                kernel_id: 17,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let initial = slab(card([0.01; 2], 0.10001), 1.0, vec![0.1; 4], &cx);
        let binding = plate(&initial, PlateMaterialModel::Isotropic, &cx);
        let gas = atmosphere::atmosphere(0.5, 2.0);
        let saved = gas.clone();
        let unsupported = atmosphere::sorption(0.25)
            .with_surface_transfer(emission([0.004; 2], 0.6))
            .unwrap();
        // Initial RH is supported; a later bracket trial leaves the actual
        // transfer source. Neither gas nor specimen may be partly advanced.
        assert!(
            gas.advance(
                &initial,
                2.0,
                &[Some(unsupported), None],
                config(),
                atmosphere::controls(),
                &cx
            )
            .is_err()
        );
        let boundary = [Some(law([0.004; 2])), Some(law([0.004; 2]))];
        let transported = gas
            .advance(
                &initial,
                2.0,
                &boundary,
                config(),
                atmosphere::controls(),
                &cx,
            )
            .unwrap();
        assert!(transported.moisture().moisture()[0] > 0.10001);
        assert!(
            binding
                .advance_finite_atmosphere(
                    &gas,
                    &initial,
                    2.0,
                    &boundary,
                    config(),
                    atmosphere::controls(),
                    &cx
                )
                .is_err()
        );
        assert!(
            initial
                .advance(2.0, &source_faces(law([0.004; 2]), 1.01), config(), &cx)
                .is_err()
        );
        gate.request();
        assert!(matches!(
            gas.advance(
                &initial,
                2.0,
                &boundary,
                config(),
                atmosphere::controls(),
                &cx
            ),
            Err(AtmosphereMoistureError::Moisture(MoistureError::Cancelled))
        ));
        assert_eq!(gas, saved);
        assert_eq!(initial.time_s(), 0.0);
        assert_eq!(initial.moisture(), &[0.1; 4]);
    });
}
