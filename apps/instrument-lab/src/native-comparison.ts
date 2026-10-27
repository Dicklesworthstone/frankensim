import {el} from './util';
import type {NativePerformanceResult} from './native-wind-client.mjs';
import {inspectNativeWave, validatePressureReceipt, type NativePressureReceipt} from './native-wave.mjs';
import {NativeComparisonPlayback} from './native-comparison-playback.mjs';
import './native-comparison.css';

interface Capture {
  wave: Blob;
  bytes: Uint8Array;
  receipt: NativePressureReceipt | null;
  label: string;
  details: string;
  url: string;
}
const errorText = (error: unknown): string => error instanceof Error ? error.message : String(error);
const pa = (value: number): string => `${value.toPrecision(5)} Pa`;

/** Two immutable native renders on one listening clock and one playback gain. */
export function createNativeComparisonPanel(pausePreview: () => void) {
  const player = new NativeComparisonPlayback();
  let reference: Capture | null = null, current: Capture | null = null;
  let generation = 0, timer: number | undefined, starting = false, admitted = false;
  const name = el('input', {type:'text', maxlength:'120', 'aria-label':'Current result label', placeholder:'Label this material or atmosphere variant'});
  const keep = el('button', {type:'button', class:'btn'}, 'Keep current result as A');
  const clear = el('button', {type:'button', class:'btn'}, 'Clear reference');
  const aLabel = el('p', {class:'small'}), bLabel = el('p', {class:'small'});
  const aDownload = el('a', {class:'btn', download:'reference-a.wav'}, 'Download A');
  const bDownload = el('a', {class:'btn', download:'candidate-b.wav'}, 'Download B');
  const aSource = el('pre', {class:'nc-source'}), bSource = el('pre', {class:'nc-source'});
  const metrics = el('p', {class:'nc-metrics small', role:'status', 'aria-live':'polite'});
  const note = el('p', {class:'small muted'});
  const play = el('button', {type:'button', class:'btn primary'}, 'Play comparison');
  const stopButton = el('button', {type:'button', class:'btn'}, 'Pause');
  const aButton = el('button', {type:'button', class:'btn', 'aria-pressed':'true'}, 'Listen A');
  const bButton = el('button', {type:'button', class:'btn', 'aria-pressed':'false'}, 'Listen B');
  const seek = el('input', {type:'range', min:'0', max:'0', step:'0.001', value:'0', 'aria-label':'Shared simulation time in seconds'});
  const position = el('span', {class:'small muted'});
  const volume = el('input', {type:'range', min:'0', max:'100', step:'1', value:'25', 'aria-label':'Shared A and B playback volume'});
  const volumeText = el('span', {}, '25%');
  const controls = (): void => {
    play.disabled = !admitted || starting || player.playing;
    stopButton.disabled = !starting && !player.playing;
    aButton.disabled = bButton.disabled = !admitted;
    seek.disabled = !admitted || starting;
    keep.disabled = !current?.receipt || current.receipt.clippedSamples !== 0 || starting;
    clear.disabled = !reference || starting;
    aButton.setAttribute('aria-pressed', String(player.side === 'a'));
    bButton.setAttribute('aria-pressed', String(player.side === 'b'));
    position.textContent = `${Number(seek.value).toFixed(3)} / ${player.duration.toFixed(3)} s`;
  };
  const stop = (): void => {
    if (player.playing) seek.value = String(player.positionSeconds);
    player.stop();
    if (timer !== undefined) window.clearInterval(timer);
    timer = undefined; controls();
  };
  const describe = (capture: Capture | null, role: string, label: HTMLElement,
    download: HTMLAnchorElement, source: HTMLElement): void => {
    label.textContent = capture ? `${role}: ${capture.label}` : `${role}: no render retained`;
    download.hidden = !capture;
    if (capture) download.href = capture.url; else download.removeAttribute('href');
    source.textContent = capture ? `${capture.details}\n\nNative pressure calibration:\n${JSON.stringify(capture.receipt, null, 2)}` : '';
  };
  const refresh = (): void => {
    stop(); player.clear(); admitted = false; seek.value = '0'; seek.max = '0';
    describe(reference, 'Reference A', aLabel, aDownload, aSource);
    describe(current, 'Current B', bLabel, bDownload, bSource);
    metrics.classList.remove('nw-refused');
    note.textContent = 'Measurements reconstruct quantized PCM pressure, not speaker SPL. No automatic level, onset or pitch matching is applied.';
    if (!current) metrics.textContent = 'Render a baseline, keep it as A, then replace material data or edit the atmosphere and render B.';
    else if (!current.receipt) metrics.textContent = 'This executable supplied no pressure calibration. The WAV is playable above, but cannot be used for a calibrated comparison.';
    else if (current.receipt.clippedSamples) {
      metrics.textContent = `The current native render clipped ${current.receipt.clippedSamples} samples. Increase the declared full-scale pressure and render both variants again.`;
      metrics.classList.add('nw-refused');
    } else if (!reference?.receipt) metrics.textContent = `Current render: ${current.receipt.fullScalePa} Pa full scale, no native clips. Keep this result as A to compare the next render.`;
    else {
      try {
        const result = player.setPair(reference.bytes, reference.receipt, current.bytes, current.receipt);
        admitted = true; seek.max = String(result.durationSeconds);
        const level = result.levelDifferenceDb === null ? 'level ratio undefined for a silent result'
          : `B − A level: ${result.levelDifferenceDb.toFixed(3)} dB`;
        metrics.textContent = `RMS A: ${pa(result.rmsA)} · B: ${pa(result.rmsB)} · ${level}. Peak A: ${pa(result.peakA)} · B: ${pa(result.peakB)}. RMS difference: ${pa(result.rmsDifference)}. ${result.identical ? 'Identical PCM in the compared interval.' : 'Different PCM in the compared interval.'}`;
        note.textContent += ` Both start at native simulation time zero. ${result.truncated
          ? `Unequal durations: playback and measurements use only the shared first ${result.durationSeconds.toFixed(6)} seconds; complete WAVs remain available.`
          : `Compared window: ${result.durationSeconds.toFixed(6)} seconds.`}`;
      } catch (error) { metrics.textContent = errorText(error); metrics.classList.add('nw-refused'); }
    }
    controls();
  };

  keep.addEventListener('click', () => {
    if (!current?.receipt || current.receipt.clippedSamples || starting) return;
    const next = {...current, label:name.value.trim() || current.label, url:URL.createObjectURL(current.wave)};
    if (reference) URL.revokeObjectURL(reference.url);
    reference = next; refresh();
  });
  clear.addEventListener('click', () => {
    if (reference) URL.revokeObjectURL(reference.url);
    reference = null; refresh();
  });
  name.addEventListener('input', () => {
    if (current) { current.label = name.value.trim() || 'Untitled result'; bLabel.textContent = `Current B: ${current.label}`; }
  });
  for (const [button, side] of [[aButton, 'a'], [bButton, 'b']] as const) {
    button.addEventListener('click', () => { player.select(side); controls(); });
  }
  volume.addEventListener('input', () => {
    player.setVolume(Number(volume.value) / 100); volumeText.textContent = `${volume.value}%`;
  });
  seek.addEventListener('input', () => {
    // Preserve the requested seek position rather than overwriting it with
    // the currently playing position. Resume is explicit after scrubbing.
    const requested = seek.value; stop(); seek.value = requested; controls();
  });
  stopButton.addEventListener('click', stop);
  play.addEventListener('click', () => {
    if (!admitted || starting || player.playing) return;
    pausePreview(); starting = true; controls();
    const offset = Number(seek.value) >= player.duration ? 0 : Number(seek.value);
    void player.play(offset).then(started => {
      if (!started) return;
      timer = window.setInterval(() => {
        seek.value = String(player.playing ? player.positionSeconds : player.duration);
        controls();
        if (!player.playing && timer !== undefined) { window.clearInterval(timer); timer = undefined; }
      }, 100);
    }).catch(error => { metrics.textContent = `Playback unavailable: ${errorText(error)}`; })
      .finally(() => { starting = false; controls(); });
  });
  const element = el('section', {class:'nc-panel', 'aria-label':'Calibrated material and atmosphere comparison'},
    el('h4', {}, 'Compare material and atmosphere changes'),
    el('p', {class:'small muted'}, 'Keep one native render as A. Later successful renders become B; editing inputs or a refused render never changes A. Both sides use one playback volume.'),
    el('div', {class:'btn-row'}, name, keep, clear),
    el('div', {class:'nc-captures'},
      el('div', {}, aLabel, aDownload, el('details', {}, el('summary', {}, 'A source and calibration'), aSource)),
      el('div', {}, bLabel, bDownload, el('details', {}, el('summary', {}, 'B source and calibration'), bSource))),
    metrics, note,
    el('div', {class:'btn-row'}, play, stopButton, aButton, bButton),
    el('label', {class:'nc-control'}, 'Shared simulation time ', seek, position),
    el('label', {class:'nc-control'}, 'Shared listening volume ', volume, volumeText),
    el('p', {class:'small muted'}, 'A/B switching uses a 5 ms crossfade on one running clock. It does not restart either note, normalize loudness, or alter downloaded audio. Device playback may resample; the measurements use original 48 kHz PCM.'));
  refresh();
  return {element, stop,
    async setCurrent(result: NativePerformanceResult, label: string, details: string, signal?: AbortSignal): Promise<void> {
      const version = ++generation;
      const bytes = new Uint8Array(await result.wave.arrayBuffer());
      const info = inspectNativeWave(bytes);
      const receipt = result.receipt ? validatePressureReceipt(result.receipt, info) : null;
      if (version !== generation || signal?.aborted) return;
      const next = {wave:result.wave, bytes, receipt, label, details, url:URL.createObjectURL(result.wave)};
      if (current) URL.revokeObjectURL(current.url);
      current = next; name.value = label; refresh();
    },
    dispose(): void {
      generation++; stop(); player.clear();
      if (reference) URL.revokeObjectURL(reference.url);
      if (current) URL.revokeObjectURL(current.url);
      reference = null; current = null; name.value = ''; refresh();
      void player.close().catch(() => {});
    },
  };
}
