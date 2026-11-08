//! Physical open-film mixing, changing inventory and original caloric receipts.
use super::super::tests::{
    GAP, LENGTH, WIDTH, controls, film, independent_root, insulated, near, viscosity,
};
use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::phase::{HeatCapacitySourceQuery, SolidLiquidPhase};
use fs_qty::{Dims, Temperature};

const AREA: f64 = 1e-4;
const DENSITY: f64 = 1000.0;

/// An original, source-qualified Cp(T) chart, not a surrogate total capacity.
fn caloric(
    affine: bool,
    datum: f64,
    changing_density: bool,
    phase: SolidLiquidPhase,
) -> ResolvedHeatCapacityCurve {
    let cp = PropertyKey::new("specific_heat_capacity", Dims([2, 0, -2, -1, 0, 0]));
    let rho = PropertyKey::new("density", Dims([-3, 1, 0, 0, 0, 0]));
    let mut claims = ClaimSet::new();
    for (key, dims, values) in [
        (
            cp.clone(),
            Dims([2, 0, -2, -1, 0, 0]),
            if affine {
                [1600.0, 2200.0, 3200.0]
            } else {
                [2000.0; 3]
            },
        ),
        (
            rho.clone(),
            Dims([-3, 1, 0, 0, 0, 0]),
            if changing_density {
                [1000.0, 962.5, 900.0]
            } else {
                [DENSITY; 3]
            },
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key,
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(280.0, values[0]), (310.0, values[1]), (360.0, values[2])],
                    dims,
                },
                validity: ValidityDomain::unconstrained().with("T", 280.0, 360.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic open-film caloric fixture, not measurements".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic open-film liquid".into(),
            phase: "single phase".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    ResolvedHeatCapacityCurve::try_new_exact(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &cp,
            density: &rho,
            temperature_axis: "T",
            lower: &QueryPoint::new().with("T", 280.0).unwrap(),
            upper: &QueryPoint::new().with("T", 360.0).unwrap(),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        phase,
        datum,
    )
    .unwrap()
}

fn body(
    source: &ResolvedHeatCapacityCurve,
    nonlinear: bool,
    supply: f64,
) -> OpenThermalCarreauFilm<'_> {
    OpenThermalCarreauFilm::new(
        film(300.0, nonlinear, nonlinear, supply),
        source,
        &[0.0; 3],
        4.0,
    )
    .unwrap()
}

/// Independent analytic integral of the authored affine Cp, in either direction.
fn dh(from: f64, to: f64, affine: bool) -> f64 {
    (to - from)
        * (2000.0
            + if affine {
                10.0 * ((from - 300.0) + (to - 300.0))
            } else {
                0.0
            })
}

/// Newtonian Reynolds solution with the ORIGINAL Arrhenius viscosity and
/// midpoint geometry. No production pressure/enthalpy helpers enter this oracle.
fn newtonian_oracle(
    normal_velocity: f64,
    duration: f64,
    supply_pressure: f64,
    supplies: [f64; 2],
) -> (f64, [f64; 2], f64) {
    let gap = GAP - 0.5 * duration * normal_velocity;
    let final_mass = DENSITY * AREA * (GAP - duration * normal_velocity);
    let squeeze = AREA * normal_velocity;
    let fields = |endpoint: f64| {
        let midpoint = 300.0 + 0.5 * (endpoint - 300.0);
        let conductance =
            WIDTH * gap.powi(3) / (12.0 * viscosity(midpoint, 0.0, false, false) * LENGTH);
        let flows = [
            0.5 * squeeze - conductance * supply_pressure,
            0.5 * squeeze + conductance * supply_pressure,
        ];
        let dissipation =
            squeeze.powi(2) / (2.0 * conductance) + 2.0 * conductance * supply_pressure.powi(2);
        (midpoint, flows, dissipation)
    };
    let endpoint = independent_root(280.0, 360.0, |temperature| {
        let (midpoint, flows, dissipation) = fields(temperature);
        let mut residual = final_mass * dh(300.0, temperature, false) - duration * dissipation;
        for i in 0..2 {
            let outward_mass = duration * DENSITY * flows[i];
            residual += outward_mass
                * dh(
                    300.0,
                    if outward_mass < 0.0 {
                        supplies[i]
                    } else {
                        midpoint
                    },
                    false,
                );
        }
        residual
    });
    let (_, flows, dissipation) = fields(endpoint);
    (endpoint, flows, duration * dissipation)
}

/// Check the full donor/storage balance directly, without the kernel's reduced
/// mixing equation. This also checks the persistent reference across steps.
fn audit(
    old: &OpenThermalCarreauFilm<'_>,
    step: &OpenThermalFilmStep<'_>,
    channels: &[CarreauFilmChannelReport; 2],
    supplies: [f64; 2],
    affine: bool,
) {
    let r = &step.report;
    let t = &r.thermal;
    let initial = old.temperature_k();
    let endpoint = step.next.temperature_k();
    let midpoint = initial + 0.5 * (endpoint - initial);
    let heat = t.dissipated_heat_j + t.prescribed_heat_j + t.bath_heat_j;
    let mut step_energy = heat;
    let mut port_energy = 0.0;
    let mut net_mass_into = 0.0;
    let mut supplied = 0.0;
    let mut returned = 0.0;
    assert_eq!(r.reservoirs.len(), 2);
    for (i, channel) in channels.iter().enumerate() {
        let mass_into = -t.duration_s * DENSITY * channel.volume_flow_m3_s;
        let donor = if mass_into > 0.0 {
            supplies[i]
        } else {
            midpoint
        };
        step_energy += mass_into * dh(initial, donor, affine);
        port_energy += mass_into * dh(old.reference_temperature_k(), donor, affine);
        net_mass_into += mass_into;
        if mass_into > 0.0 {
            supplied += mass_into;
        } else {
            returned -= mass_into;
        }
        let receipt = r.reservoirs[i];
        assert_eq!(receipt.channel, i);
        near(receipt.mass_into_film_kg, mass_into, 1e-21, 1e-13);
        near(
            receipt.sensible_enthalpy_into_film_j,
            mass_into * dh(old.reference_temperature_k(), donor, affine),
            1e-14,
            1e-12,
        );
        assert_eq!(
            receipt.donor_temperature_k,
            if mass_into == 0.0 { None } else { Some(donor) }
        );
    }
    near(r.initial_liquid_mass_kg, old.liquid_mass_kg(), 1e-21, 1e-13);
    near(
        r.final_liquid_mass_kg,
        step.next.liquid_mass_kg(),
        1e-21,
        1e-13,
    );
    near(
        r.final_liquid_mass_kg,
        DENSITY * AREA * (GAP - step.next.configuration()[0]),
        1e-21,
        1e-13,
    );
    near(r.supplied_liquid_mass_kg, supplied, 1e-21, 1e-12);
    near(r.returned_liquid_mass_kg, returned, 1e-21, 1e-12);
    near(
        r.mass_balance_residual_kg,
        r.final_liquid_mass_kg - r.initial_liquid_mass_kg - net_mass_into,
        2e-20,
        1e-8,
    );
    assert!(r.mass_balance_residual_kg.abs() < 1e-14);
    near(
        r.final_liquid_mass_kg * dh(initial, endpoint, affine),
        step_energy,
        5e-10,
        1e-9,
    );
    let stored = r.final_liquid_mass_kg * dh(old.reference_temperature_k(), endpoint, affine)
        - r.initial_liquid_mass_kg * dh(old.reference_temperature_k(), initial, affine);
    near(t.stored_energy_change_j, stored, 1e-14, 1e-12);
    near(t.advected_enthalpy_j, port_energy, 1e-14, 1e-12);
    near(stored, heat + port_energy, 5e-10, 1e-9);
    near(
        t.thermal_residual_j,
        stored - heat - port_energy,
        2e-14,
        1e-6,
    );
    assert!(t.thermal_residual_j.abs() <= t.thermal_residual_allowance_j);
    near(
        t.advective_continuity_defect_j,
        r.mass_balance_residual_kg * dh(old.reference_temperature_k(), midpoint, affine),
        1e-20,
        1e-10,
    );
    assert_eq!(
        step.next.reference_temperature_k(),
        old.reference_temperature_k()
    );
    assert!(std::ptr::eq(step.next.caloric(), old.caloric()));
}

#[test]
fn hot_supply_mixes_with_retained_liquid_and_pressure_reversal_changes_actual_donor() {
    let source = caloric(false, 0.0, false, SolidLiquidPhase::Liquid);
    let supplies = [330.0, 290.0];
    for supply in [2000.0, -2000.0] {
        let old = body(&source, false, supply);
        let mut force = [0.0; 3];
        let mut pressure = [0.0];
        let mut channels = [CarreauFilmChannelReport::default(); 2];
        let step = old
            .advance_into(
                &[0.0; 3],
                1.0,
                &supplies.map(Some),
                insulated(),
                controls(),
                &mut force,
                &mut pressure,
                &mut channels,
            )
            .unwrap();
        let (expected, flows, heat) = newtonian_oracle(0.0, 1.0, supply, supplies);
        near(step.next.temperature_k(), expected, 2e-8, 1e-11);
        near(step.report.thermal.dissipated_heat_j, heat, 1e-13, 1e-8);
        for i in 0..2 {
            near(channels[i].volume_flow_m3_s, flows[i], 1e-20, 1e-9);
        }
        assert_eq!(step.next.liquid_mass_kg(), old.liquid_mass_kg());
        near(
            step.report.supplied_liquid_mass_kg,
            step.report.returned_liquid_mass_kg,
            1e-19,
            1e-11,
        );
        let incoming = usize::from(supply < 0.0);
        assert!(step.report.reservoirs[incoming].mass_into_film_kg > 0.0);
        assert_eq!(
            step.report.reservoirs[incoming].donor_temperature_k,
            Some(supplies[incoming])
        );
        audit(&old, &step, &channels, supplies, false);
        assert_eq!(old.temperature_k(), 300.0);
        assert_eq!(old.time_s(), 4.0);
        assert_eq!(step.next.time_s(), 5.0);
    }
}

#[test]
fn squeeze_expels_film_enthalpy_and_refill_adds_both_reservoirs_mass_and_heat() {
    let source = caloric(false, 0.0, false, SolidLiquidPhase::Liquid);
    let supplies = [320.0, 340.0];
    for normal in [2e-4, -2e-4] {
        let old = body(&source, false, 0.0);
        let mut force = [0.0; 3];
        let mut pressure = [0.0];
        let mut channels = [CarreauFilmChannelReport::default(); 2];
        let step = old
            .advance_into(
                &[normal, 0.0, 0.0],
                0.4,
                &supplies.map(Some),
                insulated(),
                controls(),
                &mut force,
                &mut pressure,
                &mut channels,
            )
            .unwrap();
        let (expected, flows, heat) = newtonian_oracle(normal, 0.4, 0.0, supplies);
        near(step.next.temperature_k(), expected, 2e-8, 1e-11);
        near(step.report.thermal.dissipated_heat_j, heat, 1e-13, 1e-8);
        for i in 0..2 {
            near(channels[i].volume_flow_m3_s, flows[i], 1e-20, 1e-10);
        }
        audit(&old, &step, &channels, supplies, false);
        if normal > 0.0 {
            assert_eq!(step.report.supplied_liquid_mass_kg, 0.0);
            assert!(step.report.returned_liquid_mass_kg > 0.0);
            // Both ports are outflow: their unused supply temperatures have no
            // influence, even though all exterior supply states must be valid.
            let alternate = old
                .advance_into(
                    &[normal, 0.0, 0.0],
                    0.4,
                    &[Some(280.0), Some(360.0)],
                    insulated(),
                    controls(),
                    &mut force,
                    &mut pressure,
                    &mut channels,
                )
                .unwrap();
            assert_eq!(
                alternate.next.temperature_k().to_bits(),
                step.next.temperature_k().to_bits()
            );
        } else {
            assert_eq!(step.report.returned_liquid_mass_kg, 0.0);
            assert!(step.report.supplied_liquid_mass_kg > 0.0);
            assert!(step.next.temperature_k() > 308.0);
        }
    }
}

#[test]
fn nonlinear_flow_uses_native_affine_cp_and_datum_independent_continuing_inventory() {
    let zero = caloric(true, 0.0, false, SolidLiquidPhase::Liquid);
    let shifted = caloric(true, 1e12, false, SolidLiquidPhase::Liquid);
    let mut a = body(&zero, true, 2000.0);
    let mut b = body(&shifted, true, 2000.0);
    let mut first_storage = 0.0;
    for (index, (velocity, duration, supplies)) in [
        ([0.0; 3], 1.0, [350.0, 280.0]),
        ([1e-4, 0.0, 0.0], 0.2, [285.0, 340.0]),
    ]
    .into_iter()
    .enumerate()
    {
        let mut fa = [0.0; 3];
        let mut pa = [0.0];
        let mut ca = [CarreauFilmChannelReport::default(); 2];
        let sa = a
            .advance_into(
                &velocity,
                duration,
                &supplies.map(Some),
                insulated(),
                controls(),
                &mut fa,
                &mut pa,
                &mut ca,
            )
            .unwrap();
        let mut fb = [0.0; 3];
        let mut pb = [0.0];
        let mut cb = [CarreauFilmChannelReport::default(); 2];
        let sb = b
            .advance_into(
                &velocity,
                duration,
                &supplies.map(Some),
                insulated(),
                controls(),
                &mut fb,
                &mut pb,
                &mut cb,
            )
            .unwrap();
        assert_eq!((fa, pa, ca), (fb, pb, cb));
        assert_eq!(
            sa.next.temperature_k().to_bits(),
            sb.next.temperature_k().to_bits()
        );
        assert_eq!(
            sa.report.thermal.stored_energy_change_j.to_bits(),
            sb.report.thermal.stored_energy_change_j.to_bits()
        );
        assert_eq!(
            sa.report.thermal.advected_enthalpy_j.to_bits(),
            sb.report.thermal.advected_enthalpy_j.to_bits()
        );
        assert_eq!(sa.report.reservoirs, sb.report.reservoirs);
        audit(&a, &sa, &ca, supplies, true);
        audit(&b, &sb, &cb, supplies, true);
        if index == 0 {
            assert!(
                sa.next.temperature_k() > 310.0,
                "storage must cross the original Cp knot"
            );
            first_storage = sa.report.thermal.stored_energy_change_j;
        } else {
            near(
                first_storage + sa.report.thermal.stored_energy_change_j,
                sa.next.liquid_mass_kg() * dh(300.0, sa.next.temperature_k(), true),
                1e-13,
                1e-12,
            );
            assert!(sa.next.liquid_mass_kg() < a.liquid_mass_kg());
            assert_eq!(sa.next.reference_temperature_k(), 300.0);
        }
        a = sa.next;
        b = sb.next;
    }
}

#[test]
fn incompatible_caloric_sources_and_invalid_step_proposals_refuse_without_publication() {
    for (varying_density, phase) in [
        (true, SolidLiquidPhase::Liquid),
        (false, SolidLiquidPhase::Solid),
    ] {
        let source = caloric(false, 0.0, varying_density, phase);
        assert!(
            OpenThermalCarreauFilm::new(film(300.0, false, false, 0.0), &source, &[0.0; 3], 0.0,)
                .is_err()
        );
    }
    let source = caloric(false, 0.0, false, SolidLiquidPhase::Liquid);
    let old = body(&source, false, 2000.0);
    let mut limited = controls();
    limited.maximum_evaluations = 2;
    for (velocity, supplies, numerical) in [
        ([0.0; 3], [None, Some(300.0)], controls()),
        ([0.0; 3], [Some(400.0), Some(300.0)], controls()),
        ([GAP, 0.0, 0.0], [Some(300.0); 2], controls()),
        ([0.0; 3], [Some(330.0); 2], limited),
    ] {
        let mut force = [19.0; 3];
        let mut pressure = [23.0];
        let before = [CarreauFilmChannelReport {
            volume_flow_m3_s: 29.0,
            ..Default::default()
        }; 2];
        let mut channels = before;
        assert!(
            old.advance_into(
                &velocity,
                1.0,
                &supplies,
                insulated(),
                numerical,
                &mut force,
                &mut pressure,
                &mut channels,
            )
            .is_err()
        );
        assert_eq!(force, [19.0; 3]);
        assert_eq!(pressure, [23.0]);
        assert_eq!(channels, before);
        assert_eq!(old.configuration(), &[0.0; 3]);
        assert_eq!(old.temperature_k(), 300.0);
        assert_eq!(old.liquid_mass_kg(), DENSITY * AREA * GAP);
        assert_eq!(old.time_s(), 4.0);
    }
}

#[test]
fn late_cancellation_preserves_retained_mass_ports_and_exact_retry() {
    let source = caloric(true, 1e12, false, SolidLiquidPhase::Liquid);
    let old = body(&source, true, 2000.0);
    let velocity = [-1e-4, 0.0, 0.0];
    let supplies = [Some(330.0), Some(310.0)];
    let mut force = [0.0; 3];
    let mut pressure = [0.0];
    let mut channels = [CarreauFilmChannelReport::default(); 2];
    let mut calls = 0;
    let clean = old
        .advance_with_check_into(
            &velocity,
            0.2,
            &supplies,
            insulated(),
            controls(),
            &mut force,
            &mut pressure,
            &mut channels,
            || {
                calls += 1;
                Ok(())
            },
        )
        .unwrap();
    let expected = (force, pressure, channels);
    force = [19.0; 3];
    pressure = [23.0];
    channels = [CarreauFilmChannelReport {
        volume_flow_m3_s: 29.0,
        ..Default::default()
    }; 2];
    let before = channels;
    let mut at = 0;
    let error = old
        .advance_with_check_into(
            &velocity,
            0.2,
            &supplies,
            insulated(),
            controls(),
            &mut force,
            &mut pressure,
            &mut channels,
            || {
                at += 1;
                if at == calls {
                    Err(FilmError("open thermal final cancellation"))
                } else {
                    Ok(())
                }
            },
        )
        .unwrap_err();
    assert_eq!(error, FilmError("open thermal final cancellation"));
    assert_eq!(force, [19.0; 3]);
    assert_eq!(pressure, [23.0]);
    assert_eq!(channels, before);
    assert_eq!(old.configuration(), &[0.0; 3]);
    assert_eq!(old.temperature_k(), 300.0);
    assert_eq!(old.liquid_mass_kg(), DENSITY * AREA * GAP);
    assert_eq!(old.time_s(), 4.0);
    let retry = old
        .advance_into(
            &velocity,
            0.2,
            &supplies,
            insulated(),
            controls(),
            &mut force,
            &mut pressure,
            &mut channels,
        )
        .unwrap();
    assert_eq!((force, pressure, channels), expected);
    assert_eq!(
        retry.next.temperature_k().to_bits(),
        clean.next.temperature_k().to_bits()
    );
    assert_eq!(
        retry.next.liquid_mass_kg().to_bits(),
        clean.next.liquid_mass_kg().to_bits()
    );
    assert_eq!(retry.report.reservoirs, clean.report.reservoirs);
}
