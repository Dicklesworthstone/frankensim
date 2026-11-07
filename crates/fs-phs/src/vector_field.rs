//! One original compensated row for internal flow, drive and physical reaction.

use crate::interconnection::accumulation::Sum;

/// Retain the low word for a subsequent flow or power contraction. Empty rows
/// represent zero-port input. Do not round before adding another physical force.
pub(crate) fn port_input_row(g: &[f64], input: impl Fn(usize) -> f64) -> Option<Sum> {
    let mut sum = Sum::default();
    for (port, &coefficient) in g.iter().enumerate() {
        if coefficient != 0.0 {
            sum.product(coefficient, input(port))?;
        }
    }
    Some(sum)
}

/// J e - R e + G u - S e_boundary, with one final rounding after every source.
/// The enclosing pHS owner admits dimensions and finite original inputs.
pub(crate) fn row(
    j: &[f64],
    r: &[f64],
    effort: &[f64],
    g: &[f64],
    input: &[f64],
    reaction_columns: &[f64],
    reaction_effort: &[f64],
) -> Option<f64> {
    if j.len() != effort.len()
        || r.len() != effort.len()
        || g.len() != input.len()
        || reaction_columns.len() != reaction_effort.len()
    {
        return None;
    }
    let mut sum = Sum::default();
    for (col, &e) in effort.iter().enumerate() {
        if j[col] != 0.0 {
            sum.product(j[col], e)?;
        }
        if r[col] != 0.0 {
            sum.product(-r[col], e)?;
        }
    }
    sum.add_scaled(port_input_row(g, |p| input[p])?, 1.0)?;
    if !reaction_columns.is_empty() {
        sum.add_scaled(
            port_input_row(reaction_columns, |p| reaction_effort[p])?,
            -1.0,
        )?;
    }
    sum.value()
}

#[cfg(test)]
mod tests {
    use super::row;

    #[test]
    fn causal_history_reaction_preserves_weak_physical_drive() {
        // H=x²/2, J=0, R=1e20, G=S=1, x=1, f=1. A passive
        // Debye load A=-1,B=C=1,D=0 at z=-1e20 returns the strong force.
        let unloaded = row(&[0.0], &[1e20], &[1.0], &[1.0], &[1.0], &[], &[]).unwrap();
        assert_eq!(unloaded + 1e20, 0.0); // The former two-stage seam loses f.
        assert_eq!(
            row(&[0.0], &[1e20], &[1.0], &[1.0], &[1.0], &[1.0], &[-1e20]),
            Some(1.0)
        );
    }

    #[test]
    fn direct_load_and_original_host_tangent_round_together() {
        // Passive R=diag(1e-20,0), skew J01=1e20, Q=[[1e20,1],[1,1]],
        // and passive D=1 at S=[1,0]. In the first Hessian column the host
        // gyroscopic term cancels the load, leaving the original R term -1.
        let effort = [1e20, 1.0];
        let unloaded = row(&[0.0, 1e20], &[1e-20, 0.0], &effort, &[], &[], &[], &[]).unwrap();
        assert_eq!(unloaded - 1e20, 0.0);
        assert_eq!(
            row(
                &[0.0, 1e20],
                &[1e-20, 0.0],
                &effort,
                &[],
                &[],
                &[1.0],
                &[1e20]
            ),
            Some(-1.0)
        );
    }

    #[test]
    fn native_drive_and_reaction_map_remainders_are_not_dropped() {
        let big = 2.0_f64.powi(54);
        assert_eq!(
            row(
                &[0.0],
                &[0.0],
                &[0.0],
                &[big, 1.0],
                &[1.0, 1.0],
                &[big],
                &[1.0]
            ),
            Some(1.0)
        );
        assert_eq!(
            row(
                &[0.0],
                &[0.0],
                &[0.0],
                &[big],
                &[1.0],
                &[big, 1.0],
                &[1.0, 1.0]
            ),
            Some(-1.0)
        );
        assert!(row(&[0.0], &[], &[1.0], &[], &[], &[], &[]).is_none());
        assert!(row(&[0.0], &[0.0], &[1.0], &[], &[], &[1e300], &[1e300]).is_none());
    }
}
