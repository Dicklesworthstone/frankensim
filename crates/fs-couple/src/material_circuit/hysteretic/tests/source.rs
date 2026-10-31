//! G1/G3/G4 material fits alter native circuit work and actual finite core heat.
use super::*;

mod heat;
use fs_material::field_storage::{
    MAGNETIC_PERMEABILITY_DIMS,
    hysteresis::source::{HysteresisBranchSources, HysteresisSources, ResolvedScalarHysteresis},
    permanent::{FIELD_STRENGTH_DIMS, FLUX_DENSITY_DIMS},
};
use std::{cell::Cell, sync::Arc};

fn roles() -> HysteresisSources {
    HysteresisSources {
        background_permeability: PropertyKey::new("play_background", MAGNETIC_PERMEABILITY_DIMS),
        field_limit: PropertyKey::new("play_field_limit", FIELD_STRENGTH_DIMS),
        branches: vec![HysteresisBranchSources {
            susceptibility: PropertyKey::new("play_slope", MAGNETIC_PERMEABILITY_DIMS),
            pinning_field: PropertyKey::new("play_pinning", FIELD_STRENGTH_DIMS),
            saturation_polarization: PropertyKey::new("play_saturation", FLUX_DENSITY_DIMS),
        }],
    }
}

pub(super) fn core_material(
    pinning: f64,
    upper_temperature: f64,
) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let (base, _) = material(1.0, 0.0);
    let mut claims = base.claims().clone();
    let keys = roles();
    for (key, value) in [
        (keys.background_permeability, 0.01),
        (keys.field_limit, 30.0),
        (keys.branches[0].susceptibility.clone(), 0.2),
        (keys.branches[0].pinning_field.clone(), pinning),
        (keys.branches[0].saturation_polarization.clone(), 1.0),
    ] {
        let dims = key.dims();
        claims
            .insert_claim(PropertyClaim {
                key,
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    290.0,
                    upper_temperature,
                ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic sourced play coefficients".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap();
    let curve = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &point(200.0),
            upper: &point(2000.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-6,
    )
    .unwrap();
    (card, curve)
}

pub(super) fn sourced_core(card: &MaterialCard) -> UniformHystereticCore {
    UniformHystereticCore::from_source(
        ResolvedScalarHysteresis::resolve(
            Arc::new(card.clone()),
            &point(300.0),
            roles(),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap(),
        1.0,
        1.0,
        1.0,
    )
    .unwrap()
}

fn winding() -> [CellResistor; 1] {
    [CellResistor {
        destination: ResistorCell {
            branch: 1,
            cell: CellId(11),
        },
        selection: MaterialPropertySelection::SingleClaimOnly,
        length_m: 1.0,
        area_m2: 1.0,
    }]
}

#[test]
fn changing_only_core_data_changes_voltage_driven_current_and_both_heat_histories() {
    let mut results = Vec::new();
    for pinning in [1.0, 2.0] {
        let (a, ca) = core_material(pinning, 310.0);
        let (b, cb) = material(1.0, 0.0);
        let cards = [a, b];
        let curves = [ca, cb];
        let context = point(300.0);
        let cells = cells(&cards, &curves, &context);
        let resistors = winding();
        let network = MaterialNetwork {
            cells: &cells,
            resistors: &resistors,
            links: &[],
        };
        let thermal = thermal_initial(&curves);
        let sources = network.resolve(&thermal, || false).unwrap().resistors;
        let core = sourced_core(&cards[0]);
        let start = initial(&core, 0.0, vec![5.0], vec![]);
        let circuit = HystereticCircuit::new(rl_graph(), 2, core, sources, start).unwrap();
        let coupled = ThermalHystereticCircuit::bind(
            circuit,
            &thermal,
            &network,
            &[(1, ElectricalPropertyBasis::Resistivity)],
            CellId(10),
            [290.0, 310.0],
            || false,
        )
        .unwrap();
        let before = coupled.circuit().clone();
        let cells_before = thermal.cells().clone();
        let next = coupled
            .step_candidate(
                &thermal,
                0.05,
                &[5.0],
                heat_controls(),
                controls(),
                &[0.0; 2],
                || false,
            )
            .unwrap();

        // Independent first-monotone-segment integral. B=mu*H+chi*(H-hc).
        // Integral H dB = (B²+2*chi*hc*B-mu*chi*hc²)/(2*(mu+chi)).
        // With N=A=l=R=1, B=dt*(V-mean_I) gives this quadratic exactly.
        let (mu, chi, dt, voltage): (f64, f64, f64, f64) = (0.01, 0.2, 0.05, 5.0);
        let sum = mu + chi;
        let a = 1.0 + dt / (2.0 * sum);
        let b = dt * chi * pinning / sum - dt * voltage;
        let c = -dt * mu * chi * pinning * pinning / (2.0 * sum);
        let flux = (-b + (b * b - 4.0 * a * c).sqrt()) / (2.0 * a);
        let field = (flux + chi * pinning) / sum;
        let polarization = chi * (field - pinning);
        let core_heat = pinning * polarization;
        let stored = 0.5 * mu * field * field + 0.5 * polarization * polarization / chi;
        let mean_current = (stored + core_heat) / flux;
        let coil_heat = dt * mean_current * mean_current;
        let step = &next.report().electrical;
        near(next.circuit().circuit().flux_linkage_wb(), flux, 3e-9);
        near(next.circuit().circuit().current_a(), field, 3e-9);
        near(step.core.core_heat_j, core_heat, 3e-9);
        near(step.core.stored_energy_j, stored, 3e-9);
        near(step.circuit.heat_j, coil_heat, 3e-9);
        near(
            next.thermal().temperature_k(CellId(10)).unwrap(),
            300.0 + core_heat,
            2e-8,
        );
        near(
            next.thermal().temperature_k(CellId(11)).unwrap(),
            300.0 + coil_heat,
            2e-8,
        );
        near(next.report().total_energy_residual_j, 0.0, 2e-8);
        let source = next.circuit().circuit().core().source().unwrap();
        assert_eq!(
            source.material().query(),
            point(next.thermal().temperature_k(CellId(10)).unwrap())
        );
        assert_eq!(source.material().card_identity(), cards[0].content_hash());
        assert_eq!(coupled.circuit(), &before);
        assert_eq!(thermal.cells(), &cells_before);
        results.push((field, core_heat, coil_heat));

        // The next accepted interval continues the original polarization and
        // starts from the actual warmed source query, rather than a virgin fit.
        let (t, c, _) = next.into_parts();
        let memory = c.circuit().magnetic_state().clone();
        let continued = c
            .step_candidate(
                &t,
                0.1,
                &[5.0],
                heat_controls(),
                controls(),
                &[0.0; 2],
                || false,
            )
            .unwrap();
        assert_eq!(c.circuit().magnetic_state(), &memory);
        assert!(
            continued
                .circuit()
                .circuit()
                .magnetic_state()
                .polarizations_t()[0]
                > memory.polarizations_t()[0]
        );
        assert!(
            continued
                .circuit()
                .circuit()
                .core()
                .source()
                .unwrap()
                .material()
                .identity()
                != source_identity(&c)
        );
    }
    assert!(results[1].0 > results[0].0);
    assert!(results[1].1 > results[0].1 && results[1].2 > results[0].2);
}

fn source_identity(circuit: &ThermalHystereticCircuit<'_, '_>) -> fs_blake3::ContentHash {
    circuit
        .circuit()
        .core()
        .source()
        .unwrap()
        .material()
        .identity()
}

#[test]
fn wrong_core_card_unsupported_heating_and_late_cancel_preserve_complete_state() {
    let (a, ca) = core_material(2.0, 300.1);
    let (b, cb) = material(1.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let context = point(300.0);
    let cells = cells(&cards, &curves, &context);
    let resistors = winding();
    let network = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &[],
    };
    let thermal = thermal_initial(&curves);
    let sources = network.resolve(&thermal, || false).unwrap().resistors;
    let core = sourced_core(&cards[0]);
    let start = initial(&core, 0.0, vec![5.0], vec![]);
    let circuit = HystereticCircuit::new(rl_graph(), 2, core, sources.clone(), start).unwrap();
    let bases = [(1, ElectricalPropertyBasis::Resistivity)];
    // Caller bounds cannot widen a source fit's measured/admitted domain.
    assert!(
        ThermalHystereticCircuit::bind(
            circuit.clone(),
            &thermal,
            &network,
            &bases,
            CellId(10),
            [290.0, 310.0],
            || false
        )
        .is_err()
    );
    let (foreign, _) = core_material(1.0, 300.1);
    let foreign_core = sourced_core(&foreign);
    let start = initial(&foreign_core, 0.0, vec![5.0], vec![]);
    let wrong = HystereticCircuit::new(rl_graph(), 2, foreign_core, sources, start).unwrap();
    assert!(
        ThermalHystereticCircuit::bind(
            wrong,
            &thermal,
            &network,
            &bases,
            CellId(10),
            [290.0, 300.1],
            || false
        )
        .is_err()
    );

    let coupled = ThermalHystereticCircuit::bind(
        circuit,
        &thermal,
        &network,
        &bases,
        CellId(10),
        [290.0, 300.1],
        || false,
    )
    .unwrap();
    let before = coupled.circuit().clone();
    let saved = thermal.cells().clone();
    assert!(
        coupled
            .step_candidate(
                &thermal,
                0.05,
                &[5.0],
                heat_controls(),
                controls(),
                &[0.0; 2],
                || false
            )
            .is_err()
    ); // Independent first step needs T=300.23657 K.
    assert_eq!(coupled.circuit(), &before);
    assert_eq!(thermal.cells(), &saved);

    let calls = Cell::new(0usize);
    let admitted = coupled
        .step_candidate(
            &thermal,
            0.005,
            &[5.0],
            heat_controls(),
            controls(),
            &[0.0; 2],
            || {
                calls.set(calls.get() + 1);
                false
            },
        )
        .unwrap();
    let total = calls.get();
    assert!(total > 10);
    calls.set(0);
    assert!(
        coupled
            .step_candidate(
                &thermal,
                0.005,
                &[5.0],
                heat_controls(),
                controls(),
                &[0.0; 2],
                || {
                    calls.set(calls.get() + 1);
                    calls.get() >= total - 1
                }
            )
            .is_err()
    );
    assert_eq!(coupled.circuit(), &before);
    assert_eq!(thermal.cells(), &saved);
    let retry = coupled
        .step_candidate(
            &thermal,
            0.005,
            &[5.0],
            heat_controls(),
            controls(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    assert_eq!(retry.circuit().circuit(), admitted.circuit().circuit());
    assert_eq!(retry.thermal().cells(), admitted.thermal().cells());
}
