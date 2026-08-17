// Capability probe + honest banner (plan §4.3: tier availability is selected
// by MEASURED capability and the banner reports the ACTUAL disabled features;
// cross-origin isolation is an enhancement, never a silent requirement).
// Bead wf-root-guzez.1.2 (E0.2).

export interface Capabilities {
  readonly crossOriginIsolated: boolean;
  readonly sharedArrayBuffer: boolean;
  readonly hardwareConcurrency: number;
  readonly webgl2: boolean;
}

export function probeCapabilities(): Capabilities {
  let webgl2 = false;
  try {
    const canvas = document.createElement("canvas");
    webgl2 = canvas.getContext("webgl2") !== null;
  } catch {
    webgl2 = false;
  }
  return {
    crossOriginIsolated: globalThis.crossOriginIsolated === true,
    sharedArrayBuffer: typeof SharedArrayBuffer !== "undefined",
    hardwareConcurrency: navigator.hardwareConcurrency ?? 1,
    webgl2,
  };
}

export function describeCapabilities(caps: Capabilities): {
  text: string;
  degraded: boolean;
} {
  const parts: string[] = [];
  let degraded = false;
  if (caps.crossOriginIsolated && caps.sharedArrayBuffer) {
    parts.push("cross-origin isolated: shared-memory transport available");
  } else {
    degraded = true;
    parts.push(
      "compatibility mode: no cross-origin isolation — transferable-buffer transport, single-threaded wasm instances",
    );
  }
  parts.push(`${caps.hardwareConcurrency} logical cores`);
  if (!caps.webgl2) {
    degraded = true;
    parts.push("WebGL2 unavailable");
  }
  return { text: parts.join(" · "), degraded };
}
