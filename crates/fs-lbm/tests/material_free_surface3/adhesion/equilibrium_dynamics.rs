//! Original equilibrium adhesion reaches both isothermal native rate owners.

use super::*;
use fs_lbm::freesurface::{
    ADVANCING_CONTACT_ANGLE_PROPERTY, ContactLineRegime2, DYNAMIC_WETTING_SPEED_AXIS, FreeSurface,
    PlanarWall2, RECEDING_CONTACT_ANGLE_PROPERTY,
};
use fs_lbm::material::free_surface::MaterialFreeSurfaceWetting2;
use fs_lbm::material::free_surface::wetting::{
    CoxVoinovControls, CoxVoinovLaw, EquilibriumWettingRateLaw, FreeSurfaceWettingError,
    FreeSurfaceWettingSample, MaterialFreeSurfaceWetting3, MolecularKineticControls,
    MolecularKineticError, MolecularKineticLaw,
};

fn source(
    bulk: &ResolvedLiquidSurfaceState,
    origin: EquilibriumAngleOrigin,
    pinning: bool,
) -> InterfaceSystemCard {
    let base = interface_card(bulk, None, Dims::NONE, false, None);
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(match origin {
            EquilibriumAngleOrigin::YoungDupre => work(0.05 * (1.0 + 1.2_f64.cos())),
            EquilibriumAngleOrigin::SourcedAngle => {
                claim(EQUILIBRIUM_CONTACT_ANGLE_PROPERTY, 1.2, Dims::NONE)
            }
        })
        .unwrap();
    if pinning {
        for (name, angle) in [
            (RECEDING_CONTACT_ANGLE_PROPERTY, 1.1),
            (ADVANCING_CONTACT_ANGLE_PROPERTY, 1.3),
        ] {
            claims.insert_claim(claim(name, angle, Dims::NONE)).unwrap();
        }
    }
    InterfaceSystemCard::assemble(
        base.surface_a().clone(),
        base.surface_b().clone(),
        base.context().clone(),
        claims,
        vec![],
    )
    .unwrap()
}

fn molecular(maximum_speed_m_s: f64) -> EquilibriumWettingRateLaw {
    EquilibriumWettingRateLaw::MolecularKinetic(
        MolecularKineticLaw::try_new(MolecularKineticControls {
            // Explicit synthetic mobility, not inferred from reversible work.
            jump_length_m: 5e-10,
            equilibrium_jump_frequency_per_s: 1e9,
            minimum_temperature_k: 290.0,
            maximum_temperature_k: 310.0,
            maximum_speed_m_s,
            minimum_angle_radians: 0.01,
            maximum_angle_radians: 3.0,
        })
        .unwrap(),
    )
}

fn rate(molecular_rate: bool, gate: &CancelGate) -> EquilibriumWettingRateLaw {
    if molecular_rate {
        molecular(100.0)
    } else {
        EquilibriumWettingRateLaw::CoxVoinov(
            CoxVoinovLaw::try_new(
                CoxVoinovControls {
                    microscopic_cutoff_m: 1e-9,
                    observation_length_m: 0.001,
                    maximum_capillary_number: 0.01,
                    minimum_angle_radians: 0.01,
                    maximum_angle_radians: 3.0,
                },
                gate,
            )
            .unwrap(),
        )
    }
}

fn dense<'s, 'm>(
    bulk: &'s ResolvedLiquidSurfaceState,
    card: &'m InterfaceSystemCard,
    origin: EquilibriumAngleOrigin,
    gate: &CancelGate,
) -> MaterialFreeSurfaceWetting2<'s, 'm> {
    plan_free_surface_from_material(bulk, 0.001, 0.00002, 0.1)
        .unwrap()
        .free_surface2(domain2(), drop2)
        .unwrap()
        .with_equilibrium_planar_wetting(
            PlanarWall2::YMin,
            card,
            origin,
            MaterialPropertySelection::SingleClaimOnly,
            256,
            gate,
        )
        .unwrap()
}

fn sparse<'s, 'm>(
    bulk: &'s ResolvedLiquidSurfaceState,
    card: &'m InterfaceSystemCard,
    origin: EquilibriumAngleOrigin,
    gate: &CancelGate,
) -> MaterialFreeSurfaceWetting3<'s, 'm> {
    plan_free_surface_from_material(bulk, 0.001, 0.00002, 0.1)
        .unwrap()
        .free_surface3(domain(), drop_shape)
        .unwrap()
        .with_equilibrium_planar_wetting(
            PlanarWall3::ZMin,
            card,
            origin,
            MaterialPropertySelection::SingleClaimOnly,
            256,
            gate,
        )
        .unwrap()
}

fn same_dense(a: &FreeSurface, b: &FreeSurface) {
    assert_eq!(a.grid.f, b.grid.f);
    assert_eq!(a.grid.flags, b.grid.flags);
    assert_eq!(a.mass, b.mass);
    assert_eq!(a.carry.to_bits(), b.carry.to_bits());
    assert_eq!(a.conversion_redistribution(), b.conversion_redistribution());
    assert_eq!(a.conversions.to_fluid, b.conversions.to_fluid);
    assert_eq!(a.conversions.to_gas, b.conversions.to_gas);
    assert_eq!(
        a.conversions.gas_to_interface,
        b.conversions.gas_to_interface
    );
    assert_eq!(
        a.conversions.fluid_to_interface,
        b.conversions.fluid_to_interface
    );
}

fn same_sparse(a: &FreeSurface3, b: &FreeSurface3) {
    assert_eq!(a.grid().state_bits(), b.grid().state_bits());
    assert_eq!(a.carry().to_bits(), b.carry().to_bits());
    assert_eq!(a.conversions(), b.conversions());
    assert_eq!(
        a.conversion_momentum_report(),
        b.conversion_momentum_report()
    );
    for z in 0..12 {
        for y in 0..16 {
            for x in 0..16 {
                assert_eq!(
                    a.region_mass(|i, j, k| [i, j, k] == [x, y, z]).to_bits(),
                    b.region_mass(|i, j, k| [i, j, k] == [x, y, z]).to_bits(),
                );
            }
        }
    }
}

fn forward_rate(sample: &FreeSurfaceWettingSample, origin: EquilibriumAngleOrigin) {
    assert_eq!(sample.equilibrium, Some(origin));
    assert_eq!(sample.material.properties().len(), 1);
    assert_eq!(
        sample.material.state_coordinate(DYNAMIC_WETTING_SPEED_AXIS),
        None
    );
    assert_eq!(sample.material.state_coordinate("T"), Some(300.0));
    assert_eq!(
        sample.material.state_coordinate("pressure"),
        Some(100_000.0)
    );
    assert!(sample.equilibrium_hysteresis.is_none());
    let theta = match origin {
        EquilibriumAngleOrigin::YoungDupre => (sample
            .material
            .property(EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY)
            .unwrap()
            .value_si()
            / 0.05
            - 1.0)
            .acos(),
        EquilibriumAngleOrigin::SourcedAngle => 1.2,
    };
    if let Some(prediction) = sample.molecular_kinetic {
        // Independent forward Blake--Haynes equation, not the production inverse.
        // Farrokhi et al., JFM (2025), equation (3.2), doi:10.1017/jfm.2025.10407.
        let force = 0.05 * (theta.cos() - sample.applied_angle_radians.cos());
        let argument = force * 5e-10_f64.powi(2) / (2.0 * 1.380_649e-23 * 300.0);
        let expected_speed = 2.0 * 1e9 * 5e-10 * argument.sinh();
        assert!((expected_speed - sample.signed_speed_m_s).abs() < 2e-13);
        assert!((prediction.driving_force_per_length_n_m - force).abs() < 3e-15);
        assert_eq!(prediction.signed_speed_m_s, sample.signed_speed_m_s);
        assert_eq!(prediction.temperature_k, 300.0);
        assert_eq!(prediction.surface_tension_n_m, 0.05);
        assert!(prediction.dissipation_per_length_w_m >= 0.0);
        assert_eq!(
            prediction.dissipation_per_length_w_m,
            prediction.driving_force_per_length_n_m * sample.signed_speed_m_s
        );
        assert!(sample.hydrodynamic.is_none());
    } else {
        let prediction = sample.hydrodynamic.unwrap();
        // Independent Simpson integration of the full-angle, inviscid-outer
        // forward Cox relation; doi:10.1017/jfm.2020.499, equations (1.1)-(1.3).
        let h = (sample.applied_angle_radians - theta) / 256.0;
        let integrand = |t: f64| (t - t.sin() * t.cos()) / (2.0 * t.sin());
        let mut sum = integrand(theta) + integrand(sample.applied_angle_radians);
        for i in 1..256 {
            sum += (if i % 2 == 0 { 2.0 } else { 4.0 }) * integrand(theta + f64::from(i) * h);
        }
        let ca = 0.1 * sample.signed_speed_m_s / 0.05;
        assert!((h * sum / 3.0 - ca * 1e6_f64.ln()).abs() < 3e-11);
        assert_eq!(prediction.dynamic_viscosity_pa_s, 0.1);
        assert_eq!(prediction.surface_tension_n_m, 0.05);
        assert_eq!(
            prediction.model_signed_capillary_number,
            prediction.signed_capillary_number
        );
        assert_eq!(prediction.outer_to_liquid_viscosity_ratio, 0.0);
    }
}

#[test]
fn g1_g3_equilibrium_sources_drive_both_native_kinetic_owners() {
    let bulk = liquid(0, 300.0, Some(0.05));
    let gate = CancelGate::new_clock_free();
    for origin in [
        EquilibriumAngleOrigin::YoungDupre,
        EquilibriumAngleOrigin::SourcedAngle,
    ] {
        let card = source(&bulk, origin, false);
        for molecular_rate in [false, true] {
            let mut equilibrium2 = dense(&bulk, &card, origin, &gate);
            let mut dynamic2 = dense(&bulk, &card, origin, &gate)
                .with_dynamic_equilibrium_wetting(rate(molecular_rate, &gate), "T", &gate)
                .unwrap();
            let mut equilibrium3 = sparse(&bulk, &card, origin, &gate);
            let mut dynamic3 = sparse(&bulk, &card, origin, &gate)
                .with_dynamic_equilibrium_wetting(rate(molecular_rate, &gate), "T", &gate)
                .unwrap();
            let mass2 = dynamic2.flow().ledger_mass_kg_per_m().unwrap();
            let mass3 = dynamic3.flow().ledger_mass_kg().unwrap();
            equilibrium2.step_with_gate(&gate).unwrap();
            dynamic2.step_with_gate(&gate).unwrap();
            equilibrium3.step_with_gate(&gate).unwrap();
            dynamic3.step_with_gate(&gate).unwrap();
            same_dense(equilibrium2.flow().lattice(), dynamic2.flow().lattice());
            assert_eq!(
                equilibrium3.flow().lattice().grid().state_bits(),
                dynamic3.flow().lattice().grid().state_bits()
            );
            for sample in dynamic2
                .samples()
                .values()
                .chain(dynamic3.samples().values())
            {
                forward_rate(sample, origin);
                assert_eq!(&sample.material, dynamic2.source());
            }
            assert!(
                dynamic2
                    .samples()
                    .values()
                    .any(|sample| sample.signed_speed_m_s != 0.0)
            );
            assert!(
                dynamic3
                    .samples()
                    .values()
                    .any(|sample| sample.signed_speed_m_s != 0.0)
            );
            equilibrium2.step_with_gate(&gate).unwrap();
            dynamic2.step_with_gate(&gate).unwrap();
            equilibrium3.step_with_gate(&gate).unwrap();
            dynamic3.step_with_gate(&gate).unwrap();
            assert_ne!(
                equilibrium2.flow().lattice().grid.f,
                dynamic2.flow().lattice().grid.f
            );
            assert_ne!(
                equilibrium3.flow().lattice().grid().state_bits(),
                dynamic3.flow().lattice().grid().state_bits()
            );
            assert!(
                (dynamic2.flow().ledger_mass_kg_per_m().unwrap() - mass2).abs() < mass2 * 1e-10
            );
            assert!((dynamic3.flow().ledger_mass_kg().unwrap() - mass3).abs() < mass3 * 1e-10);
        }
    }
}

#[test]
fn g1_g3_independent_pinning_preserves_equilibrium_potential_and_rate_order() {
    let bulk = liquid(0, 300.0, Some(0.05));
    let origin = EquilibriumAngleOrigin::YoungDupre;
    let card = source(&bulk, origin, true);
    let gate = CancelGate::new_clock_free();
    let mut a = dense(&bulk, &card, origin, &gate)
        .with_equilibrium_wetting_hysteresis(MaterialPropertySelection::SingleClaimOnly, 0.0, &gate)
        .unwrap()
        .with_dynamic_equilibrium_wetting(molecular(100.0), "T", &gate)
        .unwrap();
    let mut b = dense(&bulk, &card, origin, &gate)
        .with_dynamic_equilibrium_wetting(molecular(100.0), "T", &gate)
        .unwrap()
        .with_equilibrium_wetting_hysteresis(MaterialPropertySelection::SingleClaimOnly, 0.0, &gate)
        .unwrap();
    let mut c = sparse(&bulk, &card, origin, &gate)
        .with_equilibrium_wetting_hysteresis(MaterialPropertySelection::SingleClaimOnly, 0.0, &gate)
        .unwrap()
        .with_dynamic_equilibrium_wetting(molecular(100.0), "T", &gate)
        .unwrap();
    let mut d = sparse(&bulk, &card, origin, &gate)
        .with_dynamic_equilibrium_wetting(molecular(100.0), "T", &gate)
        .unwrap()
        .with_equilibrium_wetting_hysteresis(MaterialPropertySelection::SingleClaimOnly, 0.0, &gate)
        .unwrap();
    assert_eq!(a.samples(), b.samples());
    assert_eq!(c.samples(), d.samples());
    a.step_with_gate(&gate).unwrap();
    b.step_with_gate(&gate).unwrap();
    c.step_with_gate(&gate).unwrap();
    d.step_with_gate(&gate).unwrap();
    same_dense(a.flow().lattice(), b.flow().lattice());
    same_sparse(c.flow().lattice(), d.flow().lattice());
    assert_eq!(a.samples(), b.samples());
    assert_eq!(c.samples(), d.samples());
    let mut dissipating = 0;
    for sample in a.samples().values().chain(c.samples().values()) {
        let pinning = sample.equilibrium_hysteresis.as_ref().unwrap();
        let kinetic = sample.molecular_kinetic.unwrap();
        assert_eq!(&sample.material, a.source());
        assert_eq!(sample.material.properties().len(), 1);
        assert_eq!(pinning.thresholds.properties().len(), 2);
        assert!(
            pinning
                .thresholds
                .property(EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY)
                .is_none()
        );
        assert_eq!(
            pinning
                .thresholds
                .state_coordinate(DYNAMIC_WETTING_SPEED_AXIS),
            None
        );
        let branch = match sample.regime {
            ContactLineRegime2::Advancing => 1.3,
            ContactLineRegime2::Receding => 1.1,
            ContactLineRegime2::Pinned => sample.geometric_angle_radians,
            ContactLineRegime2::Equilibrium => pinning.equilibrium_angle_radians,
        };
        let force = 0.05 * (pinning.equilibrium_angle_radians.cos() - branch.cos());
        assert!((pinning.pinning_force_per_length_n_m - force).abs() < 2e-16);
        assert_eq!(kinetic.signed_speed_m_s, pinning.model_speed_m_s);
        assert!(pinning.pinning_dissipation_per_length_w_m >= 0.0);
        assert!(kinetic.dissipation_per_length_w_m >= 0.0);
        if pinning.pinning_dissipation_per_length_w_m > 0.0 {
            dissipating += 1;
        }
    }
    assert!(dissipating > 0);
}

#[test]
fn g4_equilibrium_rate_exit_rolls_back_both_native_owners() {
    let bulk = liquid(0, 300.0, Some(0.05));
    let origin = EquilibriumAngleOrigin::YoungDupre;
    let card = source(&bulk, origin, false);
    let gate = CancelGate::new_clock_free();
    let mut a = dense(&bulk, &card, origin, &gate)
        .with_dynamic_equilibrium_wetting(molecular(1e-18), "T", &gate)
        .unwrap();
    let mut b = sparse(&bulk, &card, origin, &gate)
        .with_dynamic_equilibrium_wetting(molecular(1e-18), "T", &gate)
        .unwrap();
    let before2 = a.flow().lattice().clone();
    let before3 = b.flow().lattice().clone();
    let samples2 = a.samples().clone();
    let samples3 = b.samples().clone();
    for result in [a.step_with_gate(&gate), b.step_with_gate(&gate)] {
        assert!(matches!(
            result,
            Err(FreeSurfaceWettingError::MolecularKinetic(
                MolecularKineticError::OutsideDomain(_)
            ))
        ));
    }
    same_dense(a.flow().lattice(), &before2);
    same_sparse(b.flow().lattice(), &before3);
    assert_eq!(a.samples(), &samples2);
    assert_eq!(b.samples(), &samples3);
    assert!(
        dense(&bulk, &card, origin, &gate)
            .with_dynamic_equilibrium_wetting(molecular(100.0), "pressure", &gate)
            .is_err()
    );
    assert!(
        sparse(&bulk, &card, origin, &gate)
            .with_dynamic_equilibrium_wetting(molecular(100.0), "T", &gate)
            .unwrap()
            .with_dynamic_equilibrium_wetting(rate(false, &gate), "T", &gate)
            .is_err()
    );
    assert!(
        dense(&bulk, &card, origin, &gate)
            .with_equilibrium_wetting_hysteresis(
                MaterialPropertySelection::SingleClaimOnly,
                0.0,
                &gate
            )
            .is_err()
    );
}

#[test]
fn g1_pure_equilibrium_admits_pressure_driven_flat_fill_without_velocity_closure() {
    let bulk = liquid(0, 300.0, Some(0.05));
    let origin = EquilibriumAngleOrigin::YoungDupre;
    let card = source(&bulk, origin, false);
    let gate = CancelGate::new_clock_free();
    let make = || {
        plan(&bulk)
            .with_pressure_reference(100_000.0)
            .unwrap()
            .free_surface2(domain2(), |_, y| y == 1)
            .unwrap()
            .with_equilibrium_planar_wetting(
                PlanarWall2::YMin,
                &card,
                origin,
                MaterialPropertySelection::SingleClaimOnly,
                256,
                &gate,
            )
            .unwrap()
    };
    let mut equilibrium = make();
    let mut dynamic = make()
        .with_dynamic_equilibrium_wetting(molecular(100.0), "T", &gate)
        .unwrap();
    let center = equilibrium.flow().lattice().grid.idx(8, 1);
    let initial = equilibrium.flow().lattice().fill(center);
    equilibrium
        .step_with_gas_pressure_and_momentum(&gate, 0.05, |_| 100_000.125)
        .unwrap();
    let native = equilibrium.flow().lattice();
    assert_ne!(native.fill(center), initial);
    assert_eq!(native.fill(center), native.fill(native.grid.idx(7, 1)));
    assert_eq!(native.fill(center), native.fill(native.grid.idx(9, 1)));
    assert!(
        equilibrium
            .samples()
            .values()
            .all(|s| s.signed_speed_m_s == 0.0)
    );
    let before = dynamic.flow().lattice().clone();
    let samples = dynamic.samples().clone();
    assert!(matches!(
        dynamic.step_with_gas_pressure_and_momentum(&gate, 0.05, |_| 100_000.125),
        Err(
            fs_lbm::material::free_surface::MaterialFreeSurfaceMomentumError2::Wetting(
                FreeSurfaceWettingError::Invalid(
                    "moving contact geometry has no resolved tangential gradient"
                )
            )
        ),
    ));
    same_dense(dynamic.flow().lattice(), &before);
    assert_eq!(dynamic.samples(), &samples);
}
