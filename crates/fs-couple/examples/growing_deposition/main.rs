//! Growing clamped road: actual births -> heat -> material history -> geometry.
mod driver;
mod fixture;

use fs_solid::bridge::drawing::Drawing;

fn main() -> Result<(),Box<dyn std::error::Error>> {
    let law = fixture::material();
    let mut slack = fixture::build(&law)?;
    let mut drawn = fixture::build_with_drawing(&law, Drawing::TensionN(0.002))?;
    println!("time_s,cells,length_m,mass_kg,slack_tip_z_m,drawn_tip_z_m,min_temperature_k,max_temperature_k,birth_equilibrium_pending");
    for target in [0.25,0.5,0.75,1.0,1.5,2.0,3.0] {
        let (next_slack,report) = slack.advance_to(target,fixture::stepping(0.03125),|| Ok(()))?;
        let (next_drawn,drawn_report) = drawn.advance_to(target,fixture::stepping(0.03125),|| Ok(()))?;
        assert_eq!(next_slack.thermal().cells(),next_drawn.thermal().cells());
        assert_eq!(report.thermal,drawn_report.thermal);
        let slack_tip = next_slack.beam().deformed_position_m(next_slack.beam().length_m(),0.0)?;
        let drawn_tip = next_drawn.beam().deformed_position_m(next_drawn.beam().length_m(),0.0)?;
        let temperatures: Vec<_> = next_slack.thermal().cells().keys()
            .map(|&id| next_slack.thermal().temperature_k(id)).collect::<Result<_,_>>()?;
        let min = temperatures.iter().copied().fold(f64::INFINITY,f64::min);
        let max = temperatures.iter().copied().fold(f64::NEG_INFINITY,f64::max);
        assert!((next_slack.beam().mass_kg()-next_slack.thermal().mass_kg()).abs() < 1e-14);
        assert!((next_drawn.beam().mass_kg()-next_drawn.thermal().mass_kg()).abs() < 1e-14);
        assert!(report.thermal.residual_j.abs() < 1e-9);
        if let Some(mechanics) = &drawn_report.last_mechanics {
            assert_eq!(mechanics.drawing,Drawing::TensionN(0.002));
            assert!((mechanics.drawing_force_n-0.002).abs()<1e-12);
        }
        println!("{target:.4},{},{:.8},{:.10e},{:.10e},{:.10e},{min:.6},{max:.6},{}",
            next_slack.beam().cells().len(),next_slack.beam().length_m(),next_slack.beam().mass_kg(),
            slack_tip[1],drawn_tip[1],report.birth_equilibrium_pending);
        slack=next_slack;
        drawn=next_drawn;
    }
    assert_eq!(slack.beam().cells().len(),8);
    let slack_tip=slack.beam().deformed_position_m(slack.beam().length_m(),0.0)?[1];
    let drawn_tip=drawn.beam().deformed_position_m(drawn.beam().length_m(),0.0)?[1];
    assert!(slack_tip<0.0);
    assert!(drawn_tip>slack_tip);
    Ok(())
}
