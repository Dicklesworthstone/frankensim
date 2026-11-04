//! G1/G3: original dispersive incident laws survive every native heat rebuild.

use super::*;

fn incident_card(long_wave_loss: bool) -> MaterialCard {
    let base = card("authored finite-band incident glass", 1.0, 0.0);
    let mut claims = base.claims().clone();
    let template = base.claims().claims_ordered().next().unwrap().1;
    for i in 0..9 {
        for (name, quantity, value) in [
            (
                SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES[i],
                QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0])),
                (8.5 + 0.625 * i as f64) * 1e-6,
            ),
            (
                SPECTRAL_DIELECTRIC_IOR_PROPERTIES[i],
                QuantitySpec::dimensional(Dims::NONE),
                1.75 - 0.05 * i as f64,
            ),
            (
                SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES[i],
                SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY,
                if long_wave_loss && i >= 7 { 1.0 } else { 0.0 },
            ),
        ] {
            let mut claim = template.clone();
            claim.key = PropertyKey::with_quantity(name, quantity);
            claim.value = PropertyValue::Scalar {
                value,
                dims: quantity.dims(),
            };
            claim.provenance.source = "authored dispersive incident-port regression".into();
            claims.insert_claim(claim).unwrap();
        }
    }
    MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap()
}

fn incident_index(wavelength_m: f64) -> f64 {
    1.75 - 0.05 * (wavelength_m - 8.5e-6) / 0.625e-6
}

#[test]
fn optical_only_incidence_retains_phase_and_both_port_admittances_in_split_cells() {
    let card = incident_card(false);
    let source = infrared_source(&card, true, true);
    assert!(source.resolved().property(EPS_NAME).is_none());
    let response = DielectricResponse::SpectralOptical(
        SpectralOpticalResponse::try_new((&source).into(), MAX_GAP_M).unwrap(),
    );
    let width = 2.3e-6;
    let native = body(&[&card], &[width], &[2], AREA);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    let input = [C64::new(2.0, 1.0), C64::new(-0.5, 0.2)];
    for wavelength in [8.75e-6, 13.25e-6] {
        let omega = core::f64::consts::TAU * LIGHT_SPEED / wavelength;
        let n = incident_index(wavelength);
        let uniform = LayeredWaveHeating::try_new_with_incident_response(
            &bound,
            response,
            response,
            omega,
            &[response],
            || false,
        )
        .unwrap();
        let cellwise = LayeredWaveHeating::try_new_cellwise_with_incident_response(
            &bound,
            response,
            response,
            omega,
            &[response; 2],
            || false,
        )
        .unwrap();
        // The ordinary constant-law constructor cannot use this optical-only
        // snapshot: a second static epsilon source is intentionally absent.
        assert!(
            LayeredWaveHeating::try_new(
                &bound,
                source.resolved(),
                response,
                omega,
                &[response],
                || false,
            )
            .is_err()
        );
        for owner in [&uniform, &cellwise] {
            let stack = owner.optical_layers();
            assert!(std::ptr::eq(
                stack.incident_response().material(),
                source.resolved()
            ));
            relative(
                stack.incident_permittivity_f_m(),
                n * n / (M * LIGHT_SPEED * LIGHT_SPEED),
                2e-14,
            );
            assert_eq!(stack.incident_permeability_h_m(), M);
            let cosine = 0.8;
            let kz = omega * n * cosine / LIGHT_SPEED;
            for polarization in [Polarization::Te, Polarization::Tm] {
                let scattering = stack.response(cosine, polarization).unwrap();
                relative(scattering.incident_normal_wavenumber_rad_m, kz, 2e-14);
                near(
                    scattering.transmitted_tangential_e.re,
                    (kz * width).cos(),
                    2e-13,
                );
                near(
                    scattering.transmitted_tangential_e.im,
                    -(kz * width).sin(),
                    2e-13,
                );
                assert!(scattering.reflected_power_fraction < 1e-26);
                assert_eq!(scattering.absorbed_power_fraction, 0.0);
                let y = n / (M * LIGHT_SPEED)
                    * if polarization == Polarization::Te {
                        cosine
                    } else {
                        1.0 / cosine
                    };
                let mut heat = [0.0; 2];
                let power = owner
                    .prepare(cosine, polarization, || false)
                    .unwrap()
                    .cell_power_into(LayeredWaveDrive::Coherent(input), &mut heat, || false)
                    .unwrap();
                assert_eq!(heat, [0.0; 2]);
                for side in 0..2 {
                    relative(
                        power.incoming_w_m2[side],
                        0.5 * y * input[side].norm_sq(),
                        2e-13,
                    );
                    relative(
                        power.outgoing_w_m2[side],
                        0.5 * y * input[1 - side].norm_sq(),
                        2e-13,
                    );
                }
            }
        }
    }
}

#[test]
fn monochromatic_and_broadband_feedback_keep_the_original_incident_law_per_wavelength() {
    let incident_card = incident_card(false);
    let incident_source = infrared_source(&incident_card, true, true);
    let incident = DielectricResponse::SpectralOptical(
        SpectralOpticalResponse::try_new((&incident_source).into(), MAX_GAP_M).unwrap(),
    );
    let card = infrared_card(false, 360.0, false);
    let source =
        SpectralThermoOpticalSource::try_new(infrared_source(&card, false, true), MAX_GAP_M)
            .unwrap();
    let width = 2e-6;
    let native = body(&[&card], &[width], &[2], AREA);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    let initial = initial_at(&bound, &[291.0, 327.0]);
    let before = initial.cells().clone();
    let temperatures: Vec<_> = bound
        .cell_ids()
        .iter()
        .map(|&id| initial.temperature_k(id).unwrap())
        .collect();
    let refreshed: Vec<_> = temperatures
        .iter()
        .map(|&t| source.at_query_point(&card, &at(t, 100_000.0)).unwrap())
        .collect();
    let lines = [
        line(8.75e-6, 0.83, Polarization::Te),
        line(13.25e-6, 0.74, Polarization::Tm),
    ];
    let inputs = [
        [C64::new(8e5, 0.0), C64::new(2e5, 1e5)],
        [C64::new(7e5, 1e5), C64::new(3e5, 0.0)],
    ];
    let drives = [
        LayeredWaveDrive::Coherent(inputs[0]),
        LayeredWaveDrive::Independent(inputs[1]),
    ];
    let sources = [ThermoOpticalLayer {
        source: &source,
        card: &card,
    }];
    let owner = BroadbandThermoOpticalWaveHeating::try_new_with_incident_response(
        &bound,
        incident,
        incident,
        &lines,
        &sources,
        LIMITS,
        || false,
    )
    .unwrap();
    assert!(std::ptr::eq(
        owner.incident_response().material(),
        incident_source.resolved()
    ));
    let dt = 1e-10;
    let controls = LayeredThermalSolveControls {
        cell_tolerance_j: 1e-11,
        total_tolerance_j: 1e-11,
        ..LayeredThermalSolveControls::default()
    };
    let mut workspace = LayeredThermalWorkspace::new(2).unwrap();
    let mut expected_power = [0.0; 2];
    let mut expected_lines = Vec::new();
    for (i, line) in lines.iter().enumerate() {
        // Direct original Maxwell stack, bypassing the heat adapter. Source
        // data and all cell temperatures are the actual accepted snapshots.
        let films: Vec<_> = refreshed
            .iter()
            .map(|state| DielectricLayer {
                response: state.wave_response(line.angular_frequency_rad_s).unwrap(),
                thickness_m: width / 2.0,
            })
            .collect();
        let reference = PlanarDielectricLayers::try_new_with_incident_response(
            incident,
            incident,
            line.angular_frequency_rad_s,
            &films,
        )
        .unwrap();
        let prepared = reference
            .prepare_illumination(line.incidence_cosine, line.polarization)
            .unwrap();
        let mut heat = [0.0; 2];
        let power = if i == 0 {
            prepared.evaluate(inputs[i], &mut heat).unwrap().power
        } else {
            prepared
                .evaluate_independent_sources(inputs[i], &mut heat)
                .unwrap()
        };
        expected_lines.push(power);
        let wavelength = core::f64::consts::TAU * LIGHT_SPEED / line.angular_frequency_rad_s;
        let y = incident_index(wavelength) / (M * LIGHT_SPEED)
            * if line.polarization == Polarization::Te {
                line.incidence_cosine
            } else {
                1.0 / line.incidence_cosine
            };
        for side in 0..2 {
            relative(
                power.incoming_w_m2[side],
                0.5 * y * inputs[i][side].norm_sq(),
                2e-13,
            );
        }
        let mono = ThermoOpticalWaveHeating::try_new_with_incident_response(
            &bound,
            incident,
            incident,
            line.angular_frequency_rad_s,
            line.incidence_cosine,
            line.polarization,
            &sources,
            || false,
        )
        .unwrap();
        assert!(std::ptr::eq(
            mono.incident_response().material(),
            incident_source.resolved()
        ));
        let one = mono
            .advance_to(
                &initial,
                dt,
                drives[i],
                conduction(),
                controls,
                &mut workspace,
                || false,
            )
            .unwrap();
        assert_eq!(one.applied_optical_states, refreshed);
        assert_eq!(one.step.optical_power, power);
        for cell in 0..2 {
            relative(
                one.step.cell_energy_j[cell].1,
                heat[cell] * AREA * dt,
                2e-13,
            );
            expected_power[cell] += heat[cell] * AREA;
        }
    }
    let result = owner
        .advance_to(
            &initial,
            dt,
            &drives,
            conduction(),
            controls,
            &mut workspace,
            || false,
        )
        .unwrap();
    assert_eq!(result.line_power, expected_lines);
    assert_eq!(result.applied_optical_states, refreshed);
    let capacity = RHO * AREA * (width / 2.0) * CP;
    let g_dt = AREA * K / (width / 2.0) * dt;
    let mean = 0.5 * (temperatures[0] + temperatures[1])
        + 0.5 * dt * (expected_power[0] + expected_power[1]) / capacity;
    let contrast = (0.5 * capacity * (temperatures[0] - temperatures[1])
        + 0.5 * dt * (expected_power[0] - expected_power[1]))
        / (capacity + 2.0 * g_dt);
    for (cell, expected) in [mean + contrast, mean - contrast].into_iter().enumerate() {
        let id = bound.cell_ids()[cell];
        near(
            result.step.thermal.temperature_k(id).unwrap(),
            expected,
            3e-9,
        );
        relative(
            result.step.cell_energy_j[cell].1,
            expected_power[cell] * dt,
            2e-13,
        );
    }
    near(
        result.step.heat_report.energy.stored_energy_change_j,
        result.step.optical.absorbed_j,
        3e-11,
    );
    let next = owner
        .advance_to(
            &result.step.thermal,
            2.0 * dt,
            &drives,
            conduction(),
            controls,
            &mut workspace,
            || false,
        )
        .unwrap();
    assert_eq!(next.applied_optical_states, result.endpoint_optical_states);
    assert_eq!(
        owner.incident_response().material().identity(),
        incident_source.resolved().identity()
    );
    assert_eq!(
        owner.incident_response().material().state_coordinate("T"),
        Some(300.0)
    );
    assert_eq!(initial.cells(), &before);
    assert_eq!(initial.time_s(), 0.0);
}

#[test]
fn original_static_incident_constructors_remain_exact_lossless_wrappers() {
    let card = infrared_card(true, 360.0, false);
    let source =
        SpectralThermoOpticalSource::try_new(infrared_source(&card, true, true), MAX_GAP_M)
            .unwrap();
    let outside = exterior();
    let incident = DielectricResponse::Lossless(&outside);
    let native = body(&[&card], &[2e-6], &[1], AREA);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    let initial = initial_at(&bound, &[300.0]);
    let sources = [ThermoOpticalLayer {
        source: &source,
        card: &card,
    }];
    let lines = [line(11e-6, 0.9, Polarization::Tm)];
    let drive = LayeredWaveDrive::Coherent([C64::new(1e6, 0.0), C64::new(2e5, -1e5)]);
    let mono_old = ThermoOpticalWaveHeating::try_new(
        &bound,
        &outside,
        incident,
        IR_OMEGA,
        0.9,
        Polarization::Tm,
        &sources,
        || false,
    )
    .unwrap();
    let mono_new = ThermoOpticalWaveHeating::try_new_with_incident_response(
        &bound,
        incident,
        incident,
        IR_OMEGA,
        0.9,
        Polarization::Tm,
        &sources,
        || false,
    )
    .unwrap();
    let broad_old = BroadbandThermoOpticalWaveHeating::try_new(
        &bound,
        &outside,
        incident,
        &lines,
        &sources,
        LIMITS,
        || false,
    )
    .unwrap();
    let broad_new = BroadbandThermoOpticalWaveHeating::try_new_with_incident_response(
        &bound,
        incident,
        incident,
        &lines,
        &sources,
        LIMITS,
        || false,
    )
    .unwrap();
    let mut workspace = LayeredThermalWorkspace::new(1).unwrap();
    let controls = LayeredThermalSolveControls::default();
    let a = mono_old
        .advance_to(
            &initial,
            1e-10,
            drive,
            conduction(),
            controls,
            &mut workspace,
            || false,
        )
        .unwrap();
    let b = mono_new
        .advance_to(
            &initial,
            1e-10,
            drive,
            conduction(),
            controls,
            &mut workspace,
            || false,
        )
        .unwrap();
    let c = broad_old
        .advance_to(
            &initial,
            1e-10,
            &[drive],
            conduction(),
            controls,
            &mut workspace,
            || false,
        )
        .unwrap();
    let d = broad_new
        .advance_to(
            &initial,
            1e-10,
            &[drive],
            conduction(),
            controls,
            &mut workspace,
            || false,
        )
        .unwrap();
    for step in [&b.step, &c.step, &d.step] {
        assert_eq!(step.optical_power, a.step.optical_power);
        assert_eq!(step.optical, a.step.optical);
        assert_eq!(step.cell_energy_j, a.step.cell_energy_j);
        assert_eq!(step.thermal.cells(), a.step.thermal.cells());
    }
}

#[test]
fn any_absorbing_or_out_of_band_incident_line_refuses_the_whole_request() {
    let incident_card = incident_card(true);
    let incident_source = infrared_source(&incident_card, true, true);
    let incident = DielectricResponse::SpectralOptical(
        SpectralOpticalResponse::try_new((&incident_source).into(), MAX_GAP_M).unwrap(),
    );
    let card = infrared_card(true, 360.0, false);
    let source =
        SpectralThermoOpticalSource::try_new(infrared_source(&card, true, true), MAX_GAP_M)
            .unwrap();
    let outside = exterior();
    let exit = DielectricResponse::Lossless(&outside);
    let native = body(&[&card], &[2e-6], &[1], AREA);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    let sources = [ThermoOpticalLayer {
        source: &source,
        card: &card,
    }];
    let accepted = line(8.75e-6, 1.0, Polarization::Te);
    assert!(
        ThermoOpticalWaveHeating::try_new_with_incident_response(
            &bound,
            incident,
            exit,
            accepted.angular_frequency_rad_s,
            accepted.incidence_cosine,
            accepted.polarization,
            &sources,
            || false,
        )
        .is_ok()
    );
    for refused in [
        line(13.25e-6, 1.0, Polarization::Te),
        line(8.25e-6, 1.0, Polarization::Te),
    ] {
        // The film supports both frequencies; the refusal belongs only to
        // this incident source's loss or its narrower authored spectral band.
        assert!(
            source
                .wave_response(refused.angular_frequency_rad_s)
                .is_ok()
        );
        assert!(matches!(
            BroadbandThermoOpticalWaveHeating::try_new_with_incident_response(
                &bound,
                incident,
                exit,
                &[accepted, refused],
                &sources,
                LIMITS,
                || false,
            ),
            Err(LayeredWaveHeatingError::Electromagnetic(_))
        ));
    }
    assert!(matches!(
        BroadbandThermoOpticalWaveHeating::try_new_with_incident_response(
            &bound,
            incident,
            exit,
            &[accepted],
            &sources,
            LIMITS,
            || true,
        ),
        Err(LayeredWaveHeatingError::Cancelled)
    ));
    assert_eq!(
        incident_source.resolved().state_coordinate("T"),
        Some(300.0)
    );
}
