//! A single SI player force acting through a geometry-owned displacement row.
use super::{Error, Program, StickDrive};

pub struct SpatialInput {
    pub program: Program,
    /// Physical point displacement / generalized displacement; signed, not normalized.
    pub weights: Vec<f64>,
}
impl StickDrive {
    /// Cold attachment before the first accepted tick. Rows may overlap because
    /// two physical hands can act on the same body; their forces are summed.
    pub fn add_spatial_input(&mut self, input: SpatialInput) -> Result<(), Error> {
        if self.accepted != 0 || self.inputs.len() + self.spatial_inputs.len() >= 4
            || input.weights.len() != self.force.len()
            || input.weights.iter().any(|x| !x.is_finite()) {
            return Err("spatial player needs an unstarted clock, at most four inputs and a finite full-basis row".into());
        }
        let maximum = input.weights.iter().map(|x| x.abs()).fold(0.0_f64, f64::max);
        // Admit the largest component but RETAIN every original signed weight.
        // The mechanical owner separately admits the final sum of all forces.
        input.program.admit(self.dt_s, self.steps, maximum)?;
        self.spatial_inputs.push(input);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::Input;
    fn drive() -> StickDrive {
        StickDrive::new_inputs(vec![Input { program: Program::parse("0,0\n1,0").unwrap(),
            coordinate: 1, tip_weight: 1.0 }],0.5,2,4).unwrap()
    }
    fn hand(weights: Vec<f64>) -> SpatialInput {
        SpatialInput { program: Program::parse("0,0\n0.25,4\n0.5,0\n0.75,-4\n1,0").unwrap(), weights }
    }
    #[test]
    fn signed_hand_row_preserves_work_and_one_retryable_clock() {
        let mut d=drive();let row=vec![2.0,0.0,-3.0,0.5];
        d.add_spatial_input(hand(row.clone())).unwrap();
        let external=[1.0,7.0,-2.0,3.0];let velocity=[0.2,-0.4,0.1,0.8];
        let staged=d.forces(&external).unwrap().to_vec();
        assert_eq!(staged,[5.0,7.0,-8.0,4.0]);
        let power=staged.iter().zip(external).zip(velocity).map(|((f,e),v)|(f-e)*v).sum::<f64>();
        let hand_speed=row.iter().zip(velocity).map(|(b,v)|b*v).sum::<f64>();
        assert!((power-2.0*hand_speed).abs()<1e-14);
        assert!(d.forces(&[0.0]).is_err());assert_eq!(d.accepted,0);
        assert_eq!(d.forces(&external).unwrap(),staged);
        d.accept();assert_eq!(d.forces(&external).unwrap(),[-3.0,7.0,4.0,2.0]);
    }
    #[test]
    fn shape_clock_and_force_refusals_do_not_attach_partial_inputs() {
        let mut d=drive();
        for row in [vec![],vec![0.0;4],vec![f64::NAN;4],vec![f64::MAX;4]] {
            assert!(d.add_spatial_input(hand(row)).is_err());assert!(d.spatial_inputs.is_empty());
        }
        let late=SpatialInput {program:Program::parse("0,0\n2,0").unwrap(),weights:vec![1.0;4]};
        assert!(d.add_spatial_input(late).is_err());assert!(d.spatial_inputs.is_empty());
        d.add_spatial_input(hand(vec![1.0,0.0,0.0,0.0])).unwrap();
        d.accept();assert!(d.add_spatial_input(hand(vec![1.0;4])).is_err());
        assert_eq!(d.spatial_inputs.len(),1);
    }
}
