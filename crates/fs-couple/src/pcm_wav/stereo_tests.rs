use super::*;

fn u16_at(bytes: &[u8], start: usize) -> u16 {
    u16::from_le_bytes(bytes[start..start+2].try_into().unwrap())
}
fn u32_at(bytes: &[u8], start: usize) -> u32 {
    u32::from_le_bytes(bytes[start..start+4].try_into().unwrap())
}

#[test]
fn stereo_wave_counts_frames_and_preserves_left_right_pressure_order() {
    let (wav, clips)=encode_pcm16_wav_interleaved(&[0.0,0.25,0.5,-0.5,2.0,-2.0],48_000,2,1.0).unwrap();
    assert_eq!(&wav[..4],b"RIFF");assert_eq!(&wav[8..12],b"WAVE");
    assert_eq!(&wav[12..16],b"fmt ");assert_eq!(&wav[36..40],b"data");
    assert_eq!(u32_at(&wav,4),48);assert_eq!(u32_at(&wav,16),16);
    assert_eq!(u16_at(&wav,20),1);assert_eq!(u16_at(&wav,22),2);
    assert_eq!(u32_at(&wav,24),48_000);assert_eq!(u32_at(&wav,28),192_000);
    assert_eq!(u16_at(&wav,32),4);assert_eq!(u16_at(&wav,34),16);
    assert_eq!(u32_at(&wav,40),12);assert_eq!(wav.len(),56);assert_eq!(clips,2);
    let pcm:Vec<_>=wav[44..].chunks_exact(2).map(|b|i16::from_le_bytes([b[0],b[1]])).collect();
    assert_eq!(pcm,[0,8192,16384,-16384,32767,-32768]);
}

#[test]
fn stereo_uses_the_exact_mono_quantizer_without_channel_normalization() {
    let left=[0.0,0.5,-0.5,1.0,2.0];let right=[0.125,-0.125,-2.0,0.0,0.0];
    let interleaved:Vec<_>=left.iter().zip(right).flat_map(|(&l,r)|[l,r]).collect();
    let (stereo,clips)=encode_pcm16_wav_interleaved(&interleaved,44_100,2,1.0).unwrap();
    let (a,ca)=encode_pcm16_wav(&left,44_100,1.0).unwrap();
    let (b,cb)=encode_pcm16_wav(&right,44_100,1.0).unwrap();
    assert_eq!(clips,ca+cb);
    for (i,frame) in stereo[44..].chunks_exact(4).enumerate() {
        assert_eq!(&frame[..2],&a[44+2*i..46+2*i]);
        assert_eq!(&frame[2..],&b[44+2*i..46+2*i]);
    }
    // Golden legacy mono payload and channel/rate layout, not a new encoder oracle.
    let (mono,clips)=encode_pcm16_wav(&[0.0,0.5,-0.5,1.0],8_000,1.0).unwrap();
    assert_eq!(&mono[44..],&[0,0,0,64,0,192,255,127]);assert_eq!(clips,1);
    assert_eq!(u16_at(&mono,22),1);assert_eq!(u16_at(&mono,32),2);
    assert_eq!(u32_at(&mono,28),16_000);assert_eq!(u32_at(&mono,40),8);
    assert_eq!(mono,encode_pcm16_wav_interleaved(&[0.0,0.5,-0.5,1.0],8_000,1,1.0).unwrap().0);
}

#[test]
fn invalid_stereo_frames_and_channel_byte_rates_refuse_without_padding() {
    for channels in [0,3,u16::MAX] {
        assert!(encode_pcm16_wav_interleaved(&[0.0,0.0],48_000,channels,1.0).is_err());
    }
    for samples in [&[][..],&[0.0][..],&[0.0,0.0,0.0][..],&[0.0,f64::NAN][..],&[f64::INFINITY,0.0][..]] {
        assert!(encode_pcm16_wav_interleaved(samples,48_000,2,1.0).is_err());
    }
    for rate in [0,u32::MAX,u32::MAX/4+1] {
        assert!(encode_pcm16_wav_interleaved(&[0.0;2],rate,2,1.0).is_err());
    }
    for scale in [0.0,-1.0,f64::NAN,f64::INFINITY] {
        assert!(encode_pcm16_wav_interleaved(&[0.0;2],48_000,2,scale).is_err());
    }
}
