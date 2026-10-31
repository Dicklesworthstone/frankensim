//! Source-bound line friction and hysteresis in the existing transport paths.
use super::*;
use fs_flux::capillary::ContactLineMotion;
use fs_flux::capillary::material::contact_line::ContactLineTubeLimits;
use fs_material::liquid::surface::adhesion::EquilibriumAngleOrigin;
use fs_material::liquid::surface::contact_line::{
    ADVANCING_CONTACT_ANGLE_PROPERTY, CONTACT_LINE_FRICTION_PROPERTY,
    RECEDING_CONTACT_ANGLE_PROPERTY, resolve_contact_line_wetting,
};

pub(super) fn limits() -> ContactLineTubeLimits {
    ContactLineTubeLimits {
        maximum_speed_m_s: 0.01,
        maximum_angle_departure_radians: 0.4,
    }
}

pub(super) fn sourced_wetting(
    liquid: &ResolvedLiquidSurfaceState,
    young: f64,
    friction: f64,
    receding: f64,
    advancing: f64,
) -> ResolvedLiquidWettingState {
    let equilibrium = wetting(liquid, &liquid.material().query(), young);
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (EQUILIBRIUM_CONTACT_ANGLE_PROPERTY, Dims::NONE, young),
        (
            CONTACT_LINE_FRICTION_PROPERTY,
            DYNAMIC_VISCOSITY_DIMS,
            friction,
        ),
        (RECEDING_CONTACT_ANGLE_PROPERTY, Dims::NONE, receding),
        (ADVANCING_CONTACT_ANGLE_PROPERTY, Dims::NONE, advancing),
    ] {
        claims.insert_claim(claim(name, dims, value)).unwrap();
    }
    let card = InterfaceSystemCard::assemble(
        equilibrium.surfaces()[0].clone(),
        equilibrium.surfaces()[1].clone(),
        equilibrium.context().clone(),
        claims,
        vec![],
    )
    .unwrap();
    resolve_contact_line_wetting(
        &card,
        liquid,
        EquilibriumAngleOrigin::SourcedAngle,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

#[test]
fn shared_sources_drive_two_fluid_motion_and_each_distinct_loss() {
    let (liquid, _) = fixture(1000.0, 0.05, 0.03, 1.0);
    let wetting = sourced_wetting(&liquid, 1.0, 2.0, 0.8, 1.2);
    let tube = MaterialCapillaryTube::new(&liquid, &wetting)
        .unwrap()
        .with_contact_line(limits())
        .unwrap();
    let displaced = DisplacedFixture {
        label: "same-viscosity displaced material",
        transport: FluidTransportState::try_new(1000.0, 1000.0, 0.5, 0.05).unwrap(),
        reads: Cell::new(0),
        refuse: false,
    };
    let input = MaterialMeniscusInput {
        wet_length_m: 0.02,
        applied_pressure_pa: 300.0,
        ..input()
    };
    let gravity = 2.0;
    let step = tube
        .step_two_fluid_inclined(input, &displaced, gravity)
        .unwrap();
    let line = step.contact_line().unwrap();
    // Equal densities and viscosities make both total head and resistance
    // constant, providing an independent linear-motion and loss reference.
    let head = 1000.0 * gravity * input.tube_length_m;
    let speed = (input.applied_pressure_pa + 2.0 * 0.03 * 1.2_f64.cos() / input.radius_m - head)
        / (8.0 * 0.05 * input.tube_length_m / input.radius_m.powi(2) + 2.0 * 2.0 / input.radius_m);
    let displacement = speed * input.duration_s;
    let mean_length = input.wet_length_m + 0.5 * displacement;
    let volume = PI * input.radius_m.powi(2) * displacement;
    close(
        step.discrete().meniscus.wet_length_m,
        input.wet_length_m + displacement,
    );
    close(line.initial_speed_m_s, speed);
    close(line.final_speed_m_s, speed);
    close(
        line.contact_line_viscous_dissipation_j,
        2.0 * PI * input.radius_m * 2.0 * speed.powi(2) * input.duration_s,
    );
    close(
        line.pinning_dissipation_j,
        2.0 * 0.03 * (1.0_f64.cos() - 1.2_f64.cos()) / input.radius_m * volume,
    );
    close(
        step.discrete().liquid_viscous_dissipation_j,
        8.0 * PI * 0.05 * mean_length * speed.powi(2) * input.duration_s,
    );
    close(
        step.discrete().displaced_viscous_dissipation_j,
        8.0 * PI * 0.05 * (input.tube_length_m - mean_length) * speed.powi(2) * input.duration_s,
    );
    close(
        step.discrete().meniscus.wall_surface_energy_change_j,
        -2.0 * 0.03 * 1.0_f64.cos() / input.radius_m * volume,
    );
    close(
        step.inclined_discrete().outlet_potential_energy_j,
        head * volume,
    );
    close(
        step.total_dissipation_j(),
        step.discrete().meniscus.pressure_work_j
            - step.discrete().meniscus.wall_surface_energy_change_j
            - step.inclined_discrete().outlet_potential_energy_j,
    );
    assert_eq!(step.inclined_discrete().gravitational_storage_change_j, 0.0);
    assert_eq!(line.motion, ContactLineMotion::Advancing);
    assert!(line.final_apparent_angle_radians > 1.2);
    assert!(std::ptr::eq(
        step.tube().contact_line_source().unwrap().wetting(),
        &wetting
    ));
    assert!(std::ptr::eq(step.tube().liquid(), &liquid));
    assert!(std::ptr::eq(step.displaced_fluid(), &displaced));
    assert_eq!(wetting.resolved().properties().len(), 4);
    assert_eq!(displaced.reads.get(), 1);
}

#[test]
fn nonzero_drive_pins_then_pressure_reversal_advances_and_recedes() {
    let (liquid, _) = fixture(1000.0, 0.05, 0.03, 1.0);
    let wetting = sourced_wetting(&liquid, 1.0, 2.0, 0.8, 1.2);
    let tube = MaterialCapillaryTube::new(&liquid, &wetting)
        .unwrap()
        .with_contact_line(limits())
        .unwrap();
    let initial = MaterialMeniscusInput {
        wet_length_m: 0.02,
        ..input()
    };
    let young_pressure = 2.0 * 0.03 * 1.0_f64.cos() / initial.radius_m;
    let pinned = tube
        .step_inclined(
            MaterialMeniscusInput {
                applied_pressure_pa: -young_pressure + 50.0 + 1000.0 * 9.81 * initial.wet_length_m,
                ..initial
            },
            9.81,
        )
        .unwrap();
    let line = pinned.contact_line().unwrap();
    assert_eq!(line.motion, ContactLineMotion::Pinned);
    assert_eq!(pinned.meniscus().wet_length_m, initial.wet_length_m);
    assert_eq!(pinned.transferred_mass_kg(), 0.0);
    assert_eq!(pinned.total_dissipation_j(), 0.0);
    assert!(line.final_apparent_angle_radians > 1.0 && line.final_apparent_angle_radians < 1.2);
    let mut length = initial.wet_length_m;
    for (drive, motion, threshold_angle) in [
        (200.0, ContactLineMotion::Advancing, 1.2_f64),
        (-200.0, ContactLineMotion::Receding, 0.8_f64),
    ] {
        let step = tube
            .step_horizontal(MaterialMeniscusInput {
                wet_length_m: length,
                applied_pressure_pa: drive - young_pressure,
                ..initial
            })
            .unwrap();
        let line = step.contact_line().unwrap();
        assert_eq!(line.motion, motion);
        assert!(step.meniscus().transferred_volume_m3 * drive > 0.0);
        assert!(line.contact_line_viscous_dissipation_j > 0.0);
        assert!(line.pinning_dissipation_j > 0.0);
        close(
            line.final_apparent_angle_radians.cos(),
            threshold_angle.cos() - 2.0 * line.final_speed_m_s / 0.03,
        );
        close(
            step.total_dissipation_j(),
            step.meniscus().pressure_work_j - step.meniscus().wall_surface_energy_change_j,
        );
        assert_eq!(step.tube().wetting().contact_angle_radians(), 1.0);
        length = step.meniscus().wet_length_m;
    }
}

#[test]
fn missing_sources_and_operating_domain_fail_without_changing_originals() {
    let (liquid, equilibrium_only) = fixture(1000.0, 0.05, 0.03, 1.0);
    assert!(
        MaterialCapillaryTube::new(&liquid, &equilibrium_only)
            .unwrap()
            .with_contact_line(limits())
            .is_err()
    );
    let wetting = sourced_wetting(&liquid, 1.0, 2.0, 0.8, 1.2);
    let original_liquid = liquid.clone();
    let original_wetting = wetting.clone();
    for domain in [
        ContactLineTubeLimits {
            maximum_speed_m_s: 1e-12,
            ..limits()
        },
        ContactLineTubeLimits {
            maximum_angle_departure_radians: 0.01,
            ..limits()
        },
    ] {
        let tube = MaterialCapillaryTube::new(&liquid, &wetting)
            .unwrap()
            .with_contact_line(domain)
            .unwrap();
        assert!(
            tube.step_horizontal(MaterialMeniscusInput {
                wet_length_m: 0.02,
                ..input()
            })
            .is_err()
        );
    }
    assert_eq!(liquid, original_liquid);
    assert_eq!(wetting, original_wetting);
}
