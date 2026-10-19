//! Physical center/width coordinates for the existing scan optimization study.
use super::{ScanSensitivityError, finite, invalid, positive_coordinate};
use fs_conduction::bc::history::march::moving::sensitivity::geometry::MovingHeatGeometryControl as Key;
use fs_conduction::bc::history::march::moving::sensitivity::schedule::MovingHeatScheduleGradient;

/// Geometry controls appended after selected log-power coordinates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ScanGeometryControl {
    /// A signed center value uses c = reference_m + scale_m*x. A positive
    /// scale keeps optimizer coordinates dimensionless without forbidding
    /// zero/negative world positions. Times and interpolation remain fixed.
    CenterValue { source: usize, axis: usize, value: usize, reference_m: f64, scale_m: f64 },
    /// A half width uses a = reference_m*exp(x), preserving strict positivity.
    /// This is a shared footprint parameter, not an independent per-step width.
    HalfWidth { source: usize, axis: usize, reference_m: f64 },
}
impl ScanGeometryControl {
    pub(super) fn source(self) -> usize {
        match self { Self::CenterValue {source,..} | Self::HalfWidth {source,..} => source }
    }
    pub(super) fn key(self) -> Key {
        match self {
            Self::CenterValue {axis,value,..} => Key::CenterValue {axis,value},
            Self::HalfWidth {axis,..} => Key::HalfWidth {axis},
        }
    }
    pub(super) fn initial(self, value: f64) -> Result<f64, ScanSensitivityError> {
        finite(value)?;
        match self {
            Self::CenterValue {reference_m,scale_m,..} => {
                if !reference_m.is_finite() || !scale_m.is_finite() || scale_m <= 0.0 {
                    return Err(invalid("center controls require a finite reference and positive metre scale"));
                }
                finite((value-reference_m)/scale_m)
            }
            Self::HalfWidth {reference_m,..} => {
                if !reference_m.is_finite() || reference_m <= 0.0 || value <= 0.0 {
                    return Err(invalid("width controls require positive initial and reference metres"));
                }
                finite(fs_math::det::ln(value)-fs_math::det::ln(reference_m))
            }
        }
    }
    pub(super) fn decode(self, x: f64) -> Result<f64, ScanSensitivityError> {
        finite(x)?;
        match self {
            Self::CenterValue {reference_m,scale_m,..} => finite(scale_m.mul_add(x,reference_m)),
            Self::HalfWidth {reference_m,..} => positive_coordinate(x,reference_m),
        }
    }
    pub(super) fn chain(self, physical: f64, gradient: &MovingHeatScheduleGradient)
        -> Result<f64, ScanSensitivityError>
    {
        let (derivative, factor) = match self {
            Self::CenterValue {axis,value,scale_m,..} => (
                gradient.center_values_m.get(axis).and_then(|g|g.get(value)).copied(), scale_m),
            Self::HalfWidth {axis,..} => (gradient.half_width_m.get(axis).copied(), physical),
        };
        finite(derivative.ok_or_else(||invalid("native geometry gradient changed its authored-control layout"))? * factor)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn g3_signed_centers_and_positive_widths_apply_the_correct_coordinate_chain() {
        let center=ScanGeometryControl::CenterValue {source:0,axis:1,value:0,reference_m:-2.,scale_m:0.25};
        let width=ScanGeometryControl::HalfWidth {source:0,axis:1,reference_m:0.01};
        let gradient=MovingHeatScheduleGradient {center_values_m:[vec![0.],vec![7.],vec![0.]],
            half_width_m:[0.,-3.,0.],absorbed_power_values_w:vec![0.]};
        for x in [-4.,0.,4.] {
            let c=center.decode(x).unwrap(); let w=width.decode(x).unwrap();
            assert!((center.initial(c).unwrap()-x).abs()<1e-14);
            assert!((width.initial(w).unwrap()-x).abs()<1e-14);
            assert_eq!(center.chain(c,&gradient).unwrap(),1.75);
            assert_eq!(width.chain(w,&gradient).unwrap(),-3.*w);
            assert!(w>0.);
        }
        let bad=ScanGeometryControl::CenterValue {source:0,axis:0,value:0,reference_m:0.,scale_m:0.};
        assert!(bad.initial(0.).is_err());
        assert!(center.decode(f64::NAN).is_err());
        assert!(width.decode(-1000.).is_err());
    }
}
