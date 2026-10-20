# Driven Two-Level System: Rabi Oscillations Beyond the RWA

MATLAB simulation of a two-level quantum system (atom or qubit) driven by a
classical monochromatic field. The time-dependent Schrödinger equation is
integrated numerically **without the rotating-wave approximation (RWA)** and
compared with the analytical RWA solution (Rabi formula), in order to see where
the approximation holds and where it fails.

## Model

In dimensionless units (`τ = ω₀·t`, energies in units of `ħω₀`), the state
`ψ = (c_e, c_g)` evolves under

    i · dψ/dτ = [ ½·σ_z − λ·cos(ω̄·τ + φ)·σ_x ] ψ

| Symbol | Meaning |
|---|---|
| `λ = d·E₀ / (ħ·ω₀)` | Dimensionless coupling strength (dipole `d`, field amplitude `E₀`) |
| `ω̄ = ω / ω₀` | Drive frequency relative to the transition frequency |
| `Ω = λ·ω₀` | Rabi frequency |
| `Δ = ω − ω₀` | Detuning |
| `φ` | Phase of the drive |

The system starts in the ground state (`c_g = 1`, `c_e = 0`).

**RWA solution** (Rabi formula), with `Ω_R = √(Ω² + Δ²)`:

    P_e(τ) = (Ω² / Ω_R²) · sin²( Ω_R·τ / 2 )

The exact solution uses `ode45` with tight tolerances (`RelTol = AbsTol = 1e-8`).

## What it does

The script runs five studies, each comparing the exact numerical solution with
the RWA:

1. **Ideal RWA regime:** weak coupling (`λ = 0.01`) at resonance, over 8 Rabi
   periods.
2. **Progressive increase of intensity:** `λ = 0.05, 0.1, 0.2, 0.3` at
   resonance, one Rabi period each.
3. **Bloch-Siegert shift:** for `λ = 0.1`, scans the drive frequency `ω̄` and
   finds the value that maximizes the transfer to the excited state.
4. **Off-resonance driving:** `λ = 0.05` with `ω̄ = 0.8, 1.2, 1.5`.
5. **Ultra-strong regime:** `λ = 0.6`, where the RWA breaks down.

Each figure shows the populations `|c_e|²`, `|c_g|²` (numerical and RWA) and
the population inversion `P_e − P_g`.

## Usage

1. Open MATLAB and move to the project folder.
2. Run:

```matlab
practica1_FMC
```

3. Edit the parameters at the top of each section (`lambda`, `wbar`, `phi`)
   to explore other cases. Set `guardar = 0` to avoid saving figures.

> **Note:** the exact integration uses tight tolerances over long time spans
> (the ideal regime covers thousands of time units), and the Bloch-Siegert scan
> solves the equations for 101 different frequencies, so it can take a while.

## Generated files

Figures are saved to the working folder (`.png`, 350 dpi):

| File | Contents |
|---|---|
| `regimenideal.png` | Ideal RWA regime |
| `intensidad1..4.png` | Progressive increase of intensity |
| `maximatransferencia.png` | Frequency of maximum transfer (Bloch-Siegert) |
| `fueraderesonancia1..3.png` | Off-resonance driving |
| `ultrafuerte.png` | Ultra-strong regime |

<img width="700" height="610" alt="regimenideal" src="https://github.com/user-attachments/assets/8838e71b-9984-41c2-9ca6-0d2cbdb55b4f" />

<img width="700" height="610" alt="maximatransferencia" src="https://github.com/user-attachments/assets/06ea1aa6-00e8-4872-a223-266618f6fc1b" />


## Requirements

- MATLAB R2021a or later (uses `Name=Value` syntax in `odeset`,
  `exportgraphics` and `sgtitle`)
- No additional toolboxes required
