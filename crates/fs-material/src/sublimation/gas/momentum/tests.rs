//! Original gas/chamber/caloric owners with independently manufactured states.
use super::*;
use crate::sublimation::gas::{
    GasSurfaceBinding,
    chamber::{
        advance_chamber_with_accommodation, advance_chamber_with_gas,
        mixed::advance_partitioned_chamber,
    },
    tests as fixture,
};
use fs_thermal::sublimation::{
    SolidSurfaceState, SublimationError, SublimationLaw, SublimationProperties,
};

fn controls() -> GasChamberSolve {
    GasChamberSolve {
        surface: fixture::controls(),
        maximum_surfaces: 8,
    }
}
fn frame(id: u64) -> MolecularSurfaceFrame {
    MolecularSurfaceFrame {
        surface_id: id,
        outward_normal: [0.0, 0.0, 1.0],
        centroid_m: [2.0, 0.0, 0.0],
    }
}
fn surface<'a>(
    id: u64,
    phase: &'a dyn SublimationLaw,
    binding: &'a GasSurfaceBinding,
    solid: SolidSurfaceState,
    temperature: &'a dyn Fn(f64) -> Result<f64, SublimationError>,
) -> GasChamberSurface<'a> {
    GasChamberSurface {
        id,
        law: phase,
        binding,
        solid,
        temperature_from_enthalpy: temperature,
    }
}
fn near3(a: [f64; 3], b: [f64; 3], tolerance: f64) {
    for axis in 0..3 {
        fixture::near(a[axis], b[axis], tolerance);
    }
}

#[test]
fn equilibrium_and_sealed_surfaces_retain_pressure_without_net_mass_or_source_replay() {
    struct CountedPhase {
        phase: fixture::Phase,
        calls: std::cell::Cell<usize>,
    }
    impl SublimationLaw for CountedPhase {
        fn temperature_range_k(&self) -> [f64; 2] {
            self.phase.temperature_range_k()
        }
        fn properties(&self, t: f64) -> Result<SublimationProperties, SublimationError> {
            self.calls.set(self.calls.get() + 1);
            self.phase.properties(t)
        }
    }
    let catalog = fixture::catalog(0.0);
    let gas = fixture::gas(&catalog, 400.0, 1000.0);
    let phase = CountedPhase {
        phase: fixture::phase(),
        calls: std::cell::Cell::new(0),
    };
    for beta in [0.0, 0.1, 1.0] {
        let mut binding = fixture::binding(0.0);
        binding.accommodation_coefficient = beta;
        let chart = fixture::chart(&catalog, &binding);
        let temperature = |h| fixture::temperature(h, chart);
        let rows = [surface(
            7,
            &phase,
            &binding,
            fixture::solid(400.0, chart),
            &temperature,
        )];
        phase.calls.set(0);
        let plain =
            advance_chamber_with_gas(&rows, &gas, &catalog, 0.001, controls(), || false).unwrap();
        let calls = phase.calls.replace(0);
        let result = advance_chamber_with_momentum(
            &rows,
            &[frame(7)],
            [0.0; 3],
            None,
            &gas,
            &catalog,
            0.001,
            controls(),
            || false,
        )
        .unwrap();
        assert_eq!(result.exchange(), &plain);
        assert_eq!(
            phase.calls.get(),
            calls,
            "load composition must not replay phase sources"
        );
        fixture::near(result.exchange().step().transferred_mass_kg, 0.0, 1e-20);
        let row = result.momentum().surfaces[&7];
        fixture::near(row.total_normal_impulse_n_s, 0.01, 2e-17);
        near3(row.impulse_on_surface_n_s, [0.0, 0.0, -0.01], 2e-17);
        near3(
            row.angular_impulse_on_surface_n_m_s,
            [0.0, 0.02, 0.0],
            4e-17,
        );
        near3(row.impulse_into_gas_n_s(), [0.0, 0.0, 0.01], 2e-17);
        assert_eq!(row.fixed_support_impulse_n_s(), row.impulse_into_gas_n_s());
        assert!(result.reflection().is_empty());
        assert!(result.phases().is_empty());
    }
}

#[test]
fn absent_volatile_keeps_carrier_pressure_and_adds_only_actual_emission_recoil() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let temperature = |h| fixture::temperature(h, chart);
    let phase = fixture::phase();
    let gas = fixture::gas(&catalog, 400.0, 0.0);
    let rows = [surface(
        1,
        &phase,
        &binding,
        fixture::solid(400.0, chart),
        &temperature,
    )];
    let candidate = advance_chamber_with_momentum(
        &rows,
        &[frame(1)],
        [0.0; 3],
        None,
        &gas,
        &catalog,
        0.001,
        controls(),
        || false,
    )
    .unwrap();
    let load = candidate.momentum().surfaces[&1];
    assert_eq!(load.sticking_normal_impulse_n_s, 0.0);
    fixture::near(
        load.outgoing_normal_impulse_n_s,
        0.5 * 0.1 * 1000.0 * 1e-4 * 0.001,
        2e-20,
    );
    fixture::near(
        load.reflected_normal_impulse_n_s,
        100000.0 * 1e-4 * 0.001,
        2e-17,
    );
    let carrier = SpeciesId::new("A").unwrap();
    fixture::near(
        candidate
            .exchange()
            .gas()
            .species_mass_kg(&carrier)
            .unwrap(),
        gas.species_mass_kg(&carrier).unwrap(),
        1e-14,
    );
    assert!(candidate.exchange().surfaces()[0].transferred_mass_kg > 0.0);
}

#[test]
fn explicit_maxwell_scattering_shares_original_reflected_heat_and_calorics() {
    let catalog = fixture::catalog(-1e7);
    let binding = fixture::binding(-1e7);
    let chart = fixture::chart(&catalog, &binding);
    let temperature = |h| fixture::temperature(h, chart);
    let phase = fixture::phase();
    let gas = fixture::gas(&catalog, 400.0, 1000.0);
    let rows = [surface(
        1,
        &phase,
        &binding,
        fixture::solid(500.0, chart),
        &temperature,
    )];
    let scattering = MaxwellReflection::from_diffuse_fractions(
        &catalog,
        &[
            (SpeciesId::new("A").unwrap(), 0.75),
            (fixture::species(), 0.25),
        ],
        [200.0, 900.0],
    )
    .unwrap();
    let thermal = BTreeMap::from([(1, scattering.thermal_accommodation().clone())]);
    let maps = BTreeMap::from([(1, scattering)]);
    let plain = advance_chamber_with_accommodation(
        &rows,
        &thermal,
        &gas,
        &catalog,
        0.001,
        controls(),
        || false,
    )
    .unwrap();
    let result = advance_chamber_with_momentum(
        &rows,
        &[frame(1)],
        [0.0; 3],
        Some(&maps),
        &gas,
        &catalog,
        0.001,
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(result.exchange(), &plain.0);
    assert_eq!(result.reflection(), &plain.1);
    assert!(result.reflection()[&1].heat_into_gas_j > 0.0);
    let ratio = (500.0_f64 / 400.0).sqrt();
    let expected_pressure =
        0.5 * 0.1 * (phase.properties(500.0).unwrap().saturation_pressure_pa + 1000.0)
            + 0.9 * 1000.0 * (1.0 + 0.5 * 0.25 * (ratio - 1.0))
            + 99000.0 * (1.0 + 0.5 * 0.75 * (ratio - 1.0));
    fixture::near(
        result.momentum().surfaces[&1].total_normal_impulse_n_s,
        expected_pressure * 1e-4 * 0.001,
        3e-17,
    );
}

#[test]
fn wrenches_follow_patch_rotation_origin_shift_and_canonical_order() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let temperature = |h| fixture::temperature(h, chart);
    let phase = fixture::phase();
    let gas = fixture::gas(&catalog, 400.0, 1000.0);
    let rows = [
        surface(
            9,
            &phase,
            &binding,
            fixture::solid(500.0, chart),
            &temperature,
        ),
        surface(
            2,
            &phase,
            &binding,
            fixture::solid(300.0, chart),
            &temperature,
        ),
    ];
    let frames = [
        frame(9),
        MolecularSurfaceFrame {
            surface_id: 2,
            outward_normal: [0.0, 0.0, -1.0],
            centroid_m: [-1.0, 0.0, 0.0],
        },
    ];
    let original = advance_chamber_with_momentum(
        &rows,
        &frames,
        [0.0; 3],
        None,
        &gas,
        &catalog,
        0.001,
        controls(),
        || false,
    )
    .unwrap();
    let reversed = advance_chamber_with_momentum(
        &[rows[1], rows[0]],
        &[frames[1], frames[0]],
        [0.0; 3],
        None,
        &gas,
        &catalog,
        0.001,
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(original, reversed);
    let rotate = |v: [f64; 3]| [v[2], v[0], v[1]];
    let shift = [3.0, -2.0, 1.0];
    let rotated = frames.map(|f| MolecularSurfaceFrame {
        surface_id: f.surface_id,
        outward_normal: rotate(f.outward_normal),
        centroid_m: core::array::from_fn(|a| rotate(f.centroid_m)[a] + shift[a]),
    });
    let result = advance_chamber_with_momentum(
        &rows,
        &rotated,
        shift,
        None,
        &gas,
        &catalog,
        0.001,
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(original.exchange(), result.exchange());
    near3(
        result.momentum().impulse_on_surfaces_n_s,
        rotate(original.momentum().impulse_on_surfaces_n_s),
        2e-17,
    );
    near3(
        result.momentum().angular_impulse_on_surfaces_n_m_s,
        rotate(original.momentum().angular_impulse_on_surfaces_n_m_s),
        4e-17,
    );
    let shifted_origin = advance_chamber_with_momentum(
        &rows,
        &frames,
        [1.0, 0.0, 0.0],
        None,
        &gas,
        &catalog,
        0.001,
        controls(),
        || false,
    )
    .unwrap();
    fixture::near(
        shifted_origin.momentum().angular_impulse_on_surfaces_n_m_s[1],
        original.momentum().angular_impulse_on_surfaces_n_m_s[1]
            + original.momentum().impulse_on_surfaces_n_s[2],
        4e-17,
    );
}

#[test]
fn sharp_fusion_uses_the_same_molecular_load_and_original_phase_candidate() {
    use crate::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let (hs, cp) = fixture::chart(&catalog, &binding);
    let curve = EquilibriumEnthalpyPhaseCurve::try_new(
        fs_blake3::ContentHash([74; 32]),
        vec![
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: hs - 200.0 * cp,
                temperature_k: 200.0,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 917.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: hs,
                temperature_k: 400.0,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 917.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: hs + 80000.0,
                temperature_k: 400.0,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: 1000.0,
            },
        ],
    )
    .unwrap();
    let caloric = SharpMeltingCaloric::new(&curve).unwrap();
    let temperature = |h| {
        curve
            .temperature_at_specific_enthalpy(h)
            .map_err(|_| SublimationError::Invalid("fixture phase chart"))
    };
    let phase = fixture::phase();
    let gas = fixture::gas(&catalog, 400.0, 0.0);
    let rows = [surface(
        1,
        &phase,
        &binding,
        SolidSurfaceState {
            mass_kg: 0.01,
            specific_enthalpy_j_kg: hs + 40000.0,
            area_m2: 1e-4,
        },
        &temperature,
    )];
    let sources = BTreeMap::from([(1, caloric)]);
    let plain = advance_partitioned_chamber(
        &rows,
        &sources,
        None,
        &gas,
        &catalog,
        0.001,
        controls(),
        || false,
    )
    .unwrap();
    let result = advance_partitioned_chamber_with_momentum(
        &rows,
        &sources,
        &[frame(1)],
        [0.0; 3],
        None,
        &gas,
        &catalog,
        0.001,
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(result.exchange(), plain.exchange());
    assert_eq!(result.phases(), plain.phases());
    assert!(result.phases()[&1].melted_mass_kg < 0.0);
    fixture::near(
        result.momentum().surfaces[&1].total_normal_impulse_n_s,
        0.010005,
        3e-17,
    );
}

#[test]
fn zero_exposure_is_zero_and_late_failure_keeps_every_physical_input() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let temperature = |h| fixture::temperature(h, chart);
    let phase = fixture::phase();
    let gas = fixture::gas(&catalog, 400.0, 1000.0);
    let before = gas.clone();
    let mut rows = [surface(
        1,
        &phase,
        &binding,
        fixture::solid(500.0, chart),
        &temperature,
    )];
    for duration in [0.0, 0.001] {
        if duration > 0.0 {
            rows[0].solid.area_m2 = 0.0;
        }
        let result = advance_chamber_with_momentum(
            &rows,
            &[frame(1)],
            [0.0; 3],
            None,
            &gas,
            &catalog,
            duration,
            controls(),
            || false,
        )
        .unwrap();
        assert_eq!(result.momentum().surfaces[&1].total_normal_impulse_n_s, 0.0);
    }
    rows[0].solid.area_m2 = 1e-4;
    let mut checks = 0;
    advance_chamber_with_momentum(
        &rows,
        &[frame(1)],
        [0.0; 3],
        None,
        &gas,
        &catalog,
        0.001,
        controls(),
        || {
            checks += 1;
            false
        },
    )
    .unwrap();
    let mut calls = 0;
    let error = advance_chamber_with_momentum(
        &rows,
        &[frame(1)],
        [0.0; 3],
        None,
        &gas,
        &catalog,
        0.001,
        controls(),
        || {
            calls += 1;
            calls == checks
        },
    )
    .unwrap_err();
    assert_eq!(error, GasSurfaceError::Surface(SublimationError::Cancelled));
    assert_eq!(gas, before);
    let mut invalid_frame = frame(1);
    invalid_frame.outward_normal[2] = 2.0;
    assert!(
        advance_chamber_with_momentum(
            &rows,
            &[invalid_frame],
            [0.0; 3],
            None,
            &gas,
            &catalog,
            0.001,
            controls(),
            || false
        )
        .is_err()
    );
    assert!(
        advance_chamber_with_momentum(
            &rows,
            &[frame(2)],
            [0.0; 3],
            None,
            &gas,
            &catalog,
            0.001,
            controls(),
            || false
        )
        .is_err()
    );
    assert_eq!(gas, before);
}
