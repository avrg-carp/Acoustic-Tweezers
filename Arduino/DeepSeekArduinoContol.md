# TinyLev Acoustic Levitator — Arduino Controller Code Analysis

A detailed, section-by-section explanation of the firmware that drives the **TinyLev**
acoustic levitator ("acoustic tweezer") built in the
[Instructables "Acoustic Levitator" project](https://www.instructables.com/Acoustic-Levitator/).
The original TinyLev is a multi-emitter, single-axis ultrasonic levitator introduced by
Asier Marzo and colleagues (2017) and released as an open design for levitating small
objects, droplets and even Small-Angle X-ray Scattering (SAXS) samples.

---

## 1. The Big Picture

The device consists of:

- **72 ultrasonic transducers** (40 kHz), split into a **top array of 36** and a
  **bottom array of 36**. The two arrays face each other across an air gap.
- An **Arduino Nano** (ATmega328P) as the controller.
- An **L298N dual full-bridge (H-bridge) driver board** that actually powers the two
  transducer arrays.

The physics is a **standing wave**. Both arrays emit 40 kHz ultrasound. The two
counter-propagating waves interfere and form a standing wave whose **pressure nodes**
(regions of minimum acoustic pressure) act as tiny traps. A lightweight object such as a
polystyrene bead or a water droplet sits at one of these nodes and floats.

The clever part — and the reason this code exists — is that the trap position can be moved
along the axis by changing the **phase difference between the top and bottom arrays**.
The Arduino therefore generates **two 40 kHz square waves whose relative phase it can
control in small steps**, and the L298N amplifies those signals to drive the transducers.

The code is not using a PWM peripheral to produce the transducer signals (PWM is only used
for a separate *synchronisation* clock). Instead it **bit-bangs** the two square waves on
four output pins, updating them 24 times per 40 kHz period. That gives 15° of phase
resolution and, because the updates are synchronised to a hardware 40 kHz clock, the
output frequency is locked to exactly 40 kHz.

---

## 2. Constants and the `animation` Table

```c
#define N_PORTS 1
#define N_DIVS  24
#define N_FRAMES 24
```

- **`N_DIVS = 24`** — one 40 kHz period (25 µs) is sliced into 24 "divisions". Each
  division is a time slot of roughly 25 µs / 24 ≈ **1.04 µs**. Writing a value to the
  output port once per division is what synthesises the waveforms. 24 divisions means a
  phase can be quantised in steps of 360°/24 = **15°**.
- **`N_FRAMES = 24`** — there are 24 precomputed "frames". Each frame is a full 24-division
  period that encodes the two square waves with a *specific* phase relationship. Stepping
  through frames steps the phase difference, which moves the trap.
- **`N_PORTS = 1`** — only one 8-bit output port (PORTC) is used, and only its lower 4 bits.

The heart of the program is the `animation` table:

```c
static byte animation[N_FRAMES][N_DIVS] = { ... };
```

It is a 24 × 24 array of bytes. `animation[frame][division]` is the 4-bit pattern written
to PORTC at that division of that frame. Every entry is one of four nibbles:
`0x5`, `0x9`, `0xA`, `0x6`.

| Value | Binary | bit3 (A3) | bit2 (A2) | bit1 (A1) | bit0 (A0) |
|-------|--------|-----------|-----------|-----------|-----------|
| `0x5` | `0101` | 0 | 1 | 0 | 1 |
| `0xA` | `1010` | 1 | 0 | 1 | 0 |
| `0x9` | `1001` | 1 | 0 | 0 | 1 |
| `0x6` | `0110` | 0 | 1 | 1 | 0 |

These four patterns are exactly what is needed to drive **two complementary H-bridge
channels**. Each L298N H-bridge has two logic inputs; driving them with a signal and its
inverse produces a push-pull (bipolar) square wave across the transducer array:

- **Top array** = H-bridge A, driven by `bit0` (A0 = IN1) and `bit1` (A1 = IN2), where
  `bit1` is always the complement of `bit0`.
- **Bottom array** = H-bridge B, driven by `bit2` (A2 = IN3) and `bit3` (A3 = IN4), where
  `bit3` is always the complement of `bit2`.

So the *meaning* of each nibble is:

| Value | Top array state (bit0) | Bottom array state (bit2) | Interpretation |
|-------|------------------------|---------------------------|----------------|
| `0x5` | 1 (driving "high") | 1 (driving "high") | both arrays in phase, positive half |
| `0xA` | 0 (driving "low") | 0 (driving "low") | both arrays in phase, negative half |
| `0x9` | 1 | 0 | top high, bottom low (out of phase) |
| `0x6` | 0 | 1 | top low, bottom high (out of phase) |

`0xA` is the bitwise complement of `0x5`, and `0x6` is the complement of `0x9` — so the
`(bit0, bit1)` pair and the `(bit2, bit3)` pair are always complementary, exactly as an
H-bridge requires.

### How one frame encodes the two waves

Consider **frame 0** (the first row):

```c
{0x5,0x5,...,0x5,  0xA,0xA,...,0xA}   // 12 copies of 0x5, then 12 copies of 0xA
```

Reading the 24 divisions in time order and extracting `bit0` (top) and `bit2` (bottom):

- `bit0` = `1` for divisions 0–11, then `0` for divisions 12–23 → a **50 %-duty square wave**.
- `bit2` = `1` for divisions 0–11, then `0` for divisions 12–23 → **identical** square wave.

So in frame 0 the top and bottom arrays are **exactly in phase**. This produces a centred
standing wave with a trap at the midpoint between the arrays.

Now consider **frame 6**:

```c
{0x9×6, 0x5×6, 0x6×6, 0xA×6}
```

Extracting `bit2` (bottom array):

- divisions 0–5: `0x9` → bit2 = 0
- divisions 6–11: `0x5` → bit2 = 1
- divisions 12–17: `0x6` → bit2 = 1
- divisions 18–23: `0xA` → bit2 = 0

So the bottom wave is now **high for divisions 6–17**, i.e. shifted by 6 divisions
(6 × 15° = **90°**) relative to frame 0. The top wave (`bit0`) is unchanged. The bottom
array is therefore 90° out of phase with the top, and the standing-wave node has moved.

Across the 24 frames the bottom array's phase is advanced by 15° per frame, completing a
full 360° rotation over the whole table and then wrapping back to frame 0. Stepping frames
therefore smoothly sweeps the trap along the axis.

---

## 3. The Delay Macros

```c
#define WAIT_LOT(a) __asm__ __volatile__ ("nop"); __asm__ __volatile__ ("nop"); ... // 14 nops
#define WAIT_MID(a) ... // 13 nops
#define WAIT_LIT(a) ... // 9 nops
```

These are **fine-grained timing delays** made of inline `nop` (no-operation) instructions.
Each `nop` takes one clock cycle (62.5 ns at 16 MHz). They differ by a single `nop` so the
author can balance the per-division timing against the varying amount of *other* work done
in each division of the unrolled loop (see §5).

The macros take a dummy argument `a` that is simply ignored; it exists only so the macro
"looks like" a function and can be written the same way everywhere.

The `volatile` keyword on each `nop` prevents the compiler from optimising the sequence
away, guaranteeing the exact instruction count survives compilation.

---

## 4. The Output Macro

```c
#define OUTPUT_WAVE(pointer, d)  PORTC = pointer[d*N_PORTS + 0]
```

This writes the byte at `pointer[d]` straight to the **PORTC register**. Because only the
lower four bits of PORTC are configured as outputs (§6), writing one of `0x5/0xA/0x9/0x6`
sets the four transducer-drive pins (A0–A3) to the desired pattern for that division.

`d*N_PORTS` is a leftover generalisation: with `N_PORTS = 1` it is just `d`. It was written
so the code could be extended to more output ports without changing this macro.

---

## 5. The Main Bit-Banging Loop

The program never returns from `setup()`; the real work is an infinite loop labelled
`LOOP:` with a `goto LOOP;` at the bottom. (`loop(){}` is left empty because it is never
reached.)

```c
byte* emittingPointer = &animation[frame][0];
...
LOOP:
  while(PINB & 0b00001000);   // wait for pin 11 (PB3) to go LOW
```

**Synchronisation.** `PINB & 0b00001000` tests bit 3 of port B, which is **pin 11
(PB3)**. Pin 10 generates a hardware 40 kHz square wave and is wired to pin 11 (§6). The
loop therefore blocks until the **falling edge** of the 40 kHz sync signal. This is the
crucial trick: the hardware PWM is the *master timebase*. The software only has to finish
emitting its 24 divisions *faster* than 25 µs; the `while` wait then pads the remainder so
that each period starts exactly on the falling edge. The result is a rock-steady 40 kHz
output even though the `nop`-based delays are only approximate.

Next, the 24 divisions are written out **fully unrolled**:

```c
OUTPUT_WAVE(emittingPointer, 0);  buttonsPort = PIND;              WAIT_LIT();
OUTPUT_WAVE(emittingPointer, 1);  anyButtonPressed = (...);         WAIT_MID();
OUTPUT_WAVE(emittingPointer, 2);  buttonPressed[0] = ...;           WAIT_MID();
...
OUTPUT_WAVE(emittingPointer, 7);  buttonPressed[5] = ...;           WAIT_MID();
OUTPUT_WAVE(emittingPointer, 8);                                   WAIT_LOT();
OUTPUT_WAVE(emittingPointer, 9);                                   WAIT_LOT();
...
OUTPUT_WAVE(emittingPointer, 23);
```

Two things are happening simultaneously:

1. **Waveform generation** — each division writes the appropriate nibble to PORTC.
2. **Button sampling** — the button port is read once (during division 0) and the six
   individual button bits are latched during divisions 1–7.

The reason the delays are *different* per division is that the button-handling code adds
extra instructions. Divisions 1–7 do extra work (masking and storing a button bit), so they
use the shorter `WAIT_MID` (13 nops); division 0 does one `PIND` read so it uses the even
shorter `WAIT_LIT` (9 nops); divisions 8–23 do no button work so they use the longer
`WAIT_LOT` (14 nops). The `nop` counts are tuned so that **each division takes the same
total time** — a classic bit-banging technique where instruction count is balanced against
delay length.

The last division (23) has no trailing `WAIT` because the loop immediately re-enters the
`while(PINB & 0b00001000)` wait, which absorbs the remaining time until the next falling
edge.

### Timing estimate

At 16 MHz, one cycle ≈ 62.5 ns. A store to PORTC plus address arithmetic costs roughly
7–10 cycles; each `nop` adds one more. The 24 divisions sum to ≈ 24 × ~1.04 µs ≈ **25 µs**,
i.e. one 40 kHz period, consistent with the sync clock. The exact figures are empirical;
the sync-wait makes the design tolerant of small errors.

---

## 6. `setup()` — Pins, Sync Clock, Power

```c
DDRC = 0b00001111;   // PC0–PC3 = outputs  → Arduino pins A0–A3
PORTC = 0b00000000;  // start all low
```

Only the lower four bits of PORTC are outputs. These are the four L298N logic inputs
(IN1–IN4): two per H-bridge, one bridge per transducer array.

```c
pinMode(10, OUTPUT);        // PB2 → 40 kHz sync OUT
pinMode(11, INPUT_PULLUP);  // PB3 → sync IN
// please connect pin 10 to pin 11
```

Pin 10 is wired back to pin 11 (a loopback). In a **single-board** build this lets the
board generate its own master clock and self-synchronise. In a **multi-board** build, one
board is the master and its pin 10 fans out to pin 11 of all boards, keeping several
levitators phase-locked.

```c
for (int i = 2; i < 8; ++i) pinMode(i, INPUT_PULLUP);  // D2–D7 = buttons
```

Six buttons on pins 2–7, all with internal pull-ups (so a pressed button reads LOW).

### Generating the 40 kHz sync with Timer1

```c
noInterrupts();
TCCR1A = bit(WGM10) | bit(WGM11) | bit(COM1B1); // Fast PWM, clear OC1B on match
TCCR1B = bit(WGM12) | bit(WGM13) | bit(CS10);   // Fast PWM, no prescaler
OCR1A  = (F_CPU / 40000L) - 1;
OCR1B  = (F_CPU / 40000L) / 2;
interrupts();
```

- `WGM13:WGM12:WGM11:WGM10 = 1110` selects **Fast PWM mode 14**, where the timer counts up
  to `OCR1A` and then resets.
- `CS10 = 1` sets **no prescaler**, so the timer ticks at the system clock (16 MHz).
- `COM1B1 = 1` selects "clear OC1B on compare match, set at bottom": pin 10 (PB2 = OC1B)
  is high from the start of each cycle until the counter reaches `OCR1B`, then low until the
  counter wraps.

With `OCR1A = F_CPU/40000 − 1 = 399`, the period is 400 cycles = **25 µs → 40 kHz**.
With `OCR1B = 200`, the output is a 50 %-duty square wave. `noInterrupts()`/`interrupts()`
guard the register writes so a timer interrupt cannot fire mid-configuration.

### Power and peripheral savings

```c
ADCSRA = 0;               // disable ADC
power_adc_disable();
power_spi_disable();
power_twi_disable();
power_timer0_disable();
// power_usart0_disable();   // left enabled because Serial.begin is called below
Serial.begin(115200);
```

Disabling the ADC, SPI, TWI and Timer0 removes unused peripherals. There are two benefits:
(1) lower power draw, and (2) — more important here — **less electrical noise and fewer
interrupt/timing perturbations**, which helps keep the bit-banged 40 kHz waveform clean.
Timer0 can be disabled safely because the sketch never uses `millis()`/`delay()`. USART0 is
kept enabled (and `Serial.begin` is called) even though nothing is ever printed, so the
port is available for debugging.

> Note: the commented-out block at the top of `setup()` is a self-contained snippet that
> generates a default standing-wave pattern directly into `animation[frame]` — filling the
> first half with `0b11111111` and the second half with `0`, then toggling bit 0 every
> division. It is a reference/utility fragment and is not active in this build.

---

## 7. Button Handling and Debouncing

```c
#define BUTTON_SENS 2500
```

After the 24 divisions are written, the loop checks whether any button was pressed:

```c
if (anyButtonPressed) {
  ++buttonCounter;
  if (buttonCounter > BUTTON_SENS) {
    buttonCounter = 0;
    if      (!buttonPressed[0]) { /* frame-- */ }
    else if (!buttonPressed[1]) { /* frame++ */ }
    else if (!buttonPressed[2]) { frame = 0; }
    emittingPointer = &animation[frame][0];
  }
} else {
  buttonCounter = 0;
}
```

`buttonCounter` increments once per 25 µs period *while any button is held*. The press is
only *acted upon* after `BUTTON_SENS` (= 2500) consecutive held periods:

- 2500 × 25 µs ≈ **62.5 ms** of continuous press before the first action.

This serves two purposes at once:

1. **Debounce** — a quick tap shorter than ~62 ms does not register, so contact bounce and
   accidental taps are ignored.
2. **Auto-repeat** — after acting, the counter resets to zero, so *continuing to hold* the
   button re-triggers the action every ~62 ms. This gives a smooth "hold to keep moving"
   behaviour.

### Button-to-action mapping

The buttons are read from port D (`PIND`). `buttonsPort & 0b11111100` selects bits 2–7
(pins D2–D7); bits 0–1 are the UART and are masked out.

| Button | Pin | Port bit | Action |
|--------|-----|----------|--------|
| button 0 | D2 | bit 2 | decrement `frame` (wrap around) |
| button 1 | D3 | bit 3 | increment `frame` (wrap around) |
| button 2 | D4 | bit 4 | reset `frame = 0` |
| button 3 | D5 | bit 5 | read, but no action in this build |
| button 4 | D6 | bit 6 | read, but no action in this build |
| button 5 | D7 | bit 7 | read, but no action in this build |

Because the pins use internal pull-ups, a pressed button reads **0**. Hence
`!buttonPressed[0]` is true only when D2 is actually pressed. `STEP_SIZE = 1` means each
registered press advances the frame by one step (15° of phase), so the trap moves a small,
predictable amount per press. The wrap-around logic makes the frame index behave like a
rotary control that cycles through all 24 positions.

---

## 8. The Physics: Why Phase Shift Moves the Trap

This is the conceptual payoff of the whole program.

The top and bottom arrays emit counter-propagating 40 kHz waves. At a point on the axis,
the pressure is the superposition of a wave travelling down and a wave travelling up:

```
p(z,t) = p₀·cos(ωt − kz) + p₀·cos(ωt + kz + φ)
```

where `φ` is the phase offset applied to the bottom array, `k = 2π/λ` is the wavenumber,
and `ω = 2πf`. Using the sum-to-product identity:

```
p(z,t) = 2·p₀·cos(kz + φ/2)·cos(ωt + φ/2)
```

The **nodes** (pressure minima, where the trap sits) occur where `cos(kz + φ/2) = 0`:

```
kz + φ/2 = (2n+1)·π/2   ⇒   z = [(2n+1)·π/2 − φ/2] / k
```

Differentiating with respect to `φ` gives how far a node moves per unit phase change:

```
Δz = −Δφ / (2k) = −Δφ · λ / (4π)
```

Consequences:

- A full 360° phase rotation (`Δφ = 2π`) moves the node by `λ/2`.
- A 180° shift moves it by `λ/4`.
- The sign shows the trap moves in the direction opposite to the phase change.

At 40 kHz in air (speed of sound c ≈ 343 m/s):

```
λ = c/f = 343 / 40000 ≈ 8.58 mm   ⇒   λ/2 ≈ 4.29 mm
```

So sweeping through all 24 frames (a full 360° of phase) translates the trap by **≈ 4.3 mm**.
Each 15° step — one button press — moves it by `λ/48 ≈ 0.18 mm`. That fine, repeatable
motion is exactly what makes TinyLev a usable *tweezer*: it can nudge a levitated droplet or
bead in sub-millimetre increments, or sweep it smoothly by holding a button.

---

## 9. End-to-End Summary

Stepping back, the whole sketch implements a small **two-channel, phase-adjustable 40 kHz
function generator**:

1. **Timer1 (hardware PWM)** produces a clean 40 kHz square wave on pin 10 and loops it
   back to pin 11 to act as the master clock.
2. **The main loop** waits for the falling edge of that clock, then unrolls 24 outputs to
   PORTC, one per ~1 µs division, reconstructing two complementary square-wave pairs.
3. **The `animation` table** holds 24 precomputed periods; each one encodes the top and
   bottom arrays at a specific phase difference (0°, 15°, 30°, … 345°).
4. **The L298N H-bridges** amplify each pair into a bipolar square wave that drives the top
   36 and bottom 36 transducers respectively.
5. **The buttons** change which frame is active, advancing or retarding the phase
   difference (with ~62 ms debounce/auto-repeat), which moves the acoustic trap up or down.

In one sentence: *the code generates two phase-locked 40 kHz square waves, and a lookup
table plus three buttons rotate the phase between them in 15° increments so the operator
can move a levitated object along the standing-wave axis.*

---

## 10. Reference: Key Numbers

| Quantity | Value |
|----------|-------|
| Ultrasonic frequency | 40 kHz |
| Period | 25 µs |
| Divisions per period (`N_DIVS`) | 24 |
| Time per division | ≈ 1.04 µs |
| Phase resolution | 15° (360°/24) |
| Frames (`N_FRAMES`) | 24 |
| Trap travel per full sweep | ≈ λ/2 ≈ 4.3 mm |
| Trap travel per step | ≈ 0.18 mm |
| Wavelength (air, c ≈ 343 m/s) | ≈ 8.58 mm |
| Debounce / auto-repeat threshold | 2500 periods ≈ 62.5 ms |
| Transducers | 72 total (2 × 36) |
| Driver | L298N dual H-bridge |
| Controller | Arduino Nano (ATmega328P, 16 MHz) |

---

## 11. Sources

- [Instructables — Acoustic Levitator (UpnaLab / Asier Marzo)](https://www.instructables.com/Acoustic-Levitator/) — original project and firmware.
- [TinyLev: A multi-emitter single-axis acoustic levitator (ResearchGate)](https://www.researchgate.net/publication/319046300_TinyLev_A_multi-emitter_single-axis_acoustic_levitator) — the underlying paper describing the device.
- [Elektor Magazine — Acoustic Wave Hovering](https://www.elektormagazine.com/review/acoustic-wave-hovering) — confirms the L298N dual driver supplies the 40 kHz signals to the two arrays.
