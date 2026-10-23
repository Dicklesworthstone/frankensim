use super::*;
use fs_material::gas::GasSpec;

fn ambient() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}

fn read(text: &str) -> Result<DgSourceInput, DgSourceError> {
    let mut lines = text.lines().enumerate();
    let (line, header) = lines.next().unwrap();
    let mut fields = header.split_ascii_whitespace();
    let kind = fields.next().unwrap();
    let source = DgSourceInput::read(
        kind, line + 1, &mut fields,
        &mut || {
            let (line, text) = lines.next().ok_or_else(|| bad(0, "missing source row"))?;
            Ok((line + 1, text.split_ascii_whitespace()))
        },
        &ambient(), None, &mut Vec::new(), &CancelGate::new_clock_free(),
    )?;
    assert!(lines.next().is_none());
    Ok(source)
}

const TUBE: &str = "dg_tube 0.05 0.001 open 2 2 8 50 2000 0.01 1e-11 128";
const BORE: &str = "dg_bore 2 open 2 2 8 50 2000 0.01 1e-11 160 4\ndg_segment cylinder 0.001 0.025\ndg_segment cone 0.001 0.002 0.025";

#[test]
fn g0_uniform_source_retains_original_geometry_gas_and_exact_mass_wall_owner() {
    let load = read(TUBE).unwrap().prepare(
        &ambient(), 1.0 / 96000.0, &CancelGate::new_clock_free(),
    ).unwrap();
    assert_eq!(load.duct().inlet_gas(), &ambient());
    assert_eq!(load.duct().outlet_gas(), &ambient());
    assert_eq!(load.duct().total_length(), 0.05);
    assert_eq!(load.duct().total_dofs(), 108);
    assert!(load.duct().viscothermal_model().is_some());
    assert!(load.duct().spatial_viscothermal_admission().is_none());
    assert_eq!(load.state().step_count(), 0);
}

#[test]
fn g3_authored_bore_uses_one_existing_spatial_wall_operator() {
    let load = read(BORE).unwrap().prepare(
        &ambient(), 1.0 / 96000.0, &CancelGate::new_clock_free(),
    ).unwrap();
    let duct = load.duct();
    let wall = duct.spatial_viscothermal_admission().unwrap();
    assert_eq!(wall.resolution.quadrature_points, 4);
    assert_eq!(wall.resolution.poles_per_field, 8);
    assert_eq!(wall.memory_dofs, 128);
    assert_eq!(duct.total_dofs(), 140);
    assert_eq!(duct.total_length(), 0.05);
    assert_eq!(duct.inlet_area_m2(), core::f64::consts::PI * 0.001 * 0.001);
    assert_eq!(duct.outlet_area_m2(), core::f64::consts::PI * 0.002 * 0.002);
    assert_eq!(duct.viscothermal_band_rad_s(), Some(wall.resolution.band_rad_s));
    assert!(duct.viscothermal_model().is_none());
    assert!(duct.has_viscothermal_memory());
    let uniform = read(TUBE).unwrap().prepare(
        &ambient(), 1.0 / 96000.0, &CancelGate::new_clock_free(),
    ).unwrap();
    assert_ne!(duct.operator_identity(), uniform.duct().operator_identity());
    assert_ne!(load.impedance_pa_s_m3(), uniform.impedance_pa_s_m3());
}

#[test]
fn g0_source_refuses_incomplete_geometry_unbounded_memory_and_unresolved_band() {
    for text in [
        BORE.replace("2 open 2", "2 open 1"),
        BORE.replace("160 4", "139 4"),
        BORE.replace("160 4", "160 2"),
        BORE.replace("160 4", "2048 65"),
        BORE.replace("dg_segment cone", "cone"),
        BORE.replace("0.001 0.002", "-0.001 0.002"),
        BORE.replace("50 2000", "2000 50"),
        BORE.replace("160 4", "160 4 gas_field 2 mole 1 missing"),
        TUBE.replace("128", "107"),
    ] {
        assert!(read(&text).is_err());
    }
    let source = read(&BORE.replace("50 2000", "50 50000")).unwrap();
    assert!(source.prepare(&ambient(), 1.0 / 96000.0, &CancelGate::new_clock_free()).is_err());
    let source = read(BORE).unwrap();
    let gate = CancelGate::new_clock_free();
    gate.request();
    assert!(source.prepare(&ambient(), 1.0 / 96000.0, &gate).is_err());
}
