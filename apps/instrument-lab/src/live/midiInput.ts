/** The channel messages used by the live instruments. Channels are zero-based. */
export type MidiChannelMessage =
  | { kind: "noteon" | "noteoff"; channel: number; key: number; velocity: number }
  | { kind: "cc"; channel: number; controller: number; value: number };

/** Web MIDI supplies complete messages, so running status and partial packets refuse. */
export function decodeMidiMessage(data: ArrayLike<number> | null): MidiChannelMessage | null {
  if (!data || data.length !== 3) return null;
  const status = data[0], a = data[1], b = data[2];
  if (!Number.isInteger(status) || status! < 0x80 || status! > 0xef
    || !Number.isInteger(a) || a! < 0 || a! > 127
    || !Number.isInteger(b) || b! < 0 || b! > 127) return null;
  const channel = status! & 0x0f;
  switch (status! & 0xf0) {
    case 0x80: return { kind: "noteoff", channel, key: a!, velocity: b! };
    case 0x90: return { kind: b === 0 ? "noteoff" : "noteon", channel, key: a!, velocity: b! };
    case 0xb0: return { kind: "cc", channel, controller: a!, value: b! };
    default: return null;
  }
}

export interface MidiInputState {
  status: "disconnected" | "connecting" | "connected" | "waiting" | "error";
  inputs: { id: string; name: string }[];
  selectedId: string | null;
  message: string;
}

interface MidiInputOptions {
  onMessage: (message: MidiChannelMessage, timeStamp: number, inputId: string) => void;
  onReset: () => void;
  onState: (state: MidiInputState) => void;
  requestAccess?: () => Promise<MIDIAccess>;
}

/**
 * One selected MIDI input per instrument. Permission is requested only by connect().
 * Listeners belong to this subscription; ports may also serve other instruments.
 */
export class LiveMidiInput {
  private access: MIDIAccess | null = null;
  private input: MIDIInput | null = null;
  private messageListener: ((event: MIDIMessageEvent) => void) | null = null;
  private selectedId: string | null = null;
  private status: MidiInputState["status"] = "disconnected";
  private enabled = false;
  private disposed = false;
  private accessGeneration = 0;
  private inputGeneration = 0;
  private pending: Promise<void> | null = null;
  private cancelPending: (() => void) | null = null;
  private opening: Promise<void> | null = null;
  private failure: { input: MIDIInput; message: string } | null = null;
  private readonly requestAccess: () => Promise<MIDIAccess>;
  private readonly onAccessState = (): void => { void this.reconcile(); };

  constructor(private readonly options: MidiInputOptions) {
    this.requestAccess = options.requestAccess ?? (() => {
      if (typeof navigator === "undefined" || typeof navigator.requestMIDIAccess !== "function") {
        return Promise.reject(new Error("MIDI input is unavailable in this browser."));
      }
      return navigator.requestMIDIAccess({ sysex: false });
    });
    this.report("disconnected", "MIDI is disconnected.");
  }

  connect(): Promise<void> {
    if (this.disposed) return Promise.resolve();
    if (this.pending) return this.pending;
    this.enabled = true;
    const generation = ++this.accessGeneration;
    let cancel!: () => void;
    const cancelled = new Promise<void>((resolve) => { cancel = resolve; });
    const job = Promise.resolve().then(async () => {
      if (!this.currentAccess(generation)) return;
      if (!this.access) {
        const access = await this.requestAccess();
        if (!this.currentAccess(generation)) return;
        this.access = access;
        access.addEventListener("statechange", this.onAccessState);
      }
      await this.reconcile(true);
    }).catch((error: unknown) => {
      if (!this.currentAccess(generation)) return;
      this.enabled = false;
      this.access?.removeEventListener("statechange", this.onAccessState);
      this.access = null;
      this.detachInput();
      this.report("error", error instanceof Error ? error.message : "MIDI access was refused.");
    });
    const pending: Promise<void> = Promise.race([job, cancelled]).finally(() => {
      if (this.pending === pending) {
        this.pending = null;
        this.cancelPending = null;
      }
    });
    this.pending = pending;
    this.cancelPending = cancel;
    if (!this.access) this.report("connecting", "Waiting for MIDI input permission…");
    return pending;
  }

  /** Retain this identity even while unplugged; never silently switch to another device. */
  select(id: string): void {
    if (this.disposed || !id) return;
    if (id !== this.selectedId && this.access) this.cancelConnectionWait();
    this.selectedId = id;
    if (this.enabled && this.access) void this.reconcile(true);
    else this.report(this.status, this.enabled ? "Waiting for MIDI input permission…" : "MIDI is disconnected.");
  }

  disconnect(): void {
    if (!this.disposed) this.stop();
  }

  dispose(): void {
    if (this.disposed) return;
    this.disposed = true;
    this.stop();
  }

  private currentAccess(generation: number): boolean {
    return !this.disposed && this.enabled && generation === this.accessGeneration;
  }

  private connectedInputs(): MIDIInput[] {
    return this.access ? [...this.access.inputs.values()].filter((input) => input.state === "connected") : [];
  }

  private report(status: MidiInputState["status"], message: string): void {
    this.status = status;
    if (this.disposed) return;
    this.options.onState({ status, selectedId: this.selectedId, message,
      inputs: this.connectedInputs().map((input) => ({ id: input.id, name: input.name || "MIDI input" })) });
  }

  private cancelConnectionWait(): void {
    this.cancelPending?.();
    this.cancelPending = null;
    this.pending = null;
  }

  private stop(): void {
    const active = this.enabled || this.pending !== null || this.status !== "disconnected";
    this.enabled = false;
    ++this.accessGeneration;
    this.cancelConnectionWait();
    this.access?.removeEventListener("statechange", this.onAccessState);
    this.access = null;
    this.failure = null;
    if (!this.detachInput() && active) this.options.onReset();
    if (active) this.report("disconnected", "MIDI is disconnected.");
  }

  private detachInput(): boolean {
    ++this.inputGeneration;
    const input = this.input, listener = this.messageListener;
    this.input = null;
    this.messageListener = null;
    this.opening = null;
    if (!input) return false;
    if (listener) input.removeEventListener("midimessage", listener);
    // Do not close the shared MIDIPort: another live instrument may own a listener.
    this.options.onReset();
    return true;
  }

  private reconcile(retry = false): Promise<void> {
    if (!this.enabled || !this.access || this.disposed) return Promise.resolve();
    const inputs = this.connectedInputs();
    if (this.selectedId === null) this.selectedId = inputs[0]?.id ?? null;
    const input = inputs.find((candidate) => candidate.id === this.selectedId);
    if (input && input === this.input) {
      if (this.opening) {
        this.report("connecting", `Opening ${input.name || "MIDI input"}…`);
        return this.opening;
      }
      if (input.connection === "open") {
        this.report("connected", `Connected to ${input.name || "MIDI input"}.`);
        return Promise.resolve();
      }
    }
    this.detachInput();
    if (!input) {
      this.failure = null;
      this.report("waiting", this.selectedId === null
        ? "No MIDI input is connected. Plug in a keyboard or controller."
        : "Waiting for the selected MIDI input to reconnect.");
      return Promise.resolve();
    }
    if (this.failure?.input === input && !retry) {
      this.report("error", this.failure.message);
      return Promise.resolve();
    }
    this.failure = null;
    this.input = input;
    const generation = ++this.inputGeneration;
    const current = (): boolean => !this.disposed && this.enabled && this.input === input
      && generation === this.inputGeneration;
    this.messageListener = (event): void => {
      if (!current() || input.state !== "connected" || input.connection !== "open") return;
      const message = decodeMidiMessage(event.data);
      if (message) this.options.onMessage(message, event.timeStamp, input.id);
    };
    const opening: Promise<void> = Promise.resolve().then(async () => {
      if (!current()) return;
      await input.open();
      if (!current()) return;
      if (input.state !== "connected") {
        this.detachInput();
        this.report("waiting", "Waiting for the selected MIDI input to reconnect.");
      } else if (input.connection !== "open") {
        throw new Error("The selected MIDI input did not open.");
      } else {
        this.report("connected", `Connected to ${input.name || "MIDI input"}.`);
      }
    }).catch((error: unknown) => {
      if (!current()) return;
      this.detachInput();
      const message = error instanceof Error ? error.message : "The selected MIDI input could not open.";
      this.failure = { input, message };
      this.report("error", message);
    }).finally(() => {
      if (this.opening === opening) this.opening = null;
    });
    this.opening = opening;
    input.addEventListener("midimessage", this.messageListener);
    this.report("connecting", `Opening ${input.name || "MIDI input"}…`);
    return opening;
  }
}
