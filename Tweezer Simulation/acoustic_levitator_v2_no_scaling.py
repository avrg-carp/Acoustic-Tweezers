#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
acoustic_levitator.py
=====================

Physics simulation of a 72-transducer ultrasonic acoustic levitator
(TinyLev-style, 40 kHz), with support for phased-array control to move,
rotate and shape trapped particles.

This module accompanies the README / notebook and implements:

  1. The transducer array geometry (36 top + 36 bottom, in either
     "whole-array", "flower-petal" (6 sectors x 6) or "concentric-ring"
     (3 rings) grouping).
  2. The complex acoustic pressure field produced by the array, using a
     baffled-piston (circular piston directivity) model per transducer.
  3. The Gor'kov acoustic radiation potential and the resulting radiation
     force on a small spherical particle.
  4. Phase-control recipes: anti-phase standing wave, phase-focusing,
     phase ramps (steering / translation), and acoustic vortices (rotation
     / orbital angular momentum).

Physical model
--------------
Each transducer is modelled as a circular piston of radius `a` in an
infinite baffle.  For a piston centred at r_j with outward normal n_j and
drive phase phi_j, the complex pressure at a field point r is

    p_j(r) = A_j * D(theta_j) * exp( i ( k |r - r_j| + phi_j ) ) / |r - r_j|

where k = 2*pi/lambda is the wavenumber, theta_j is the angle between n_j
and the vector (r - r_j), and

    D(theta) = 2 * J1( k a sin(theta) ) / ( k a sin(theta) )

is the on-axis-normalised directivity of a baffled circular piston
(J1 = first-order Bessel function).  The total field is the coherent sum

    p(r) = sum_j p_j(r).

The particle velocity is obtained from the pressure via the linearised
momentum equation  v = -grad(p) / (i * omega * rho0), so |v|^2 is computed
from a finite-difference gradient of the complex pressure.

For a small spherical particle (radius << wavelength) the time-averaged
acoustic radiation force is the negative gradient of the Gor'kov potential

    U(r) = Vp * [ f1 * <p^2> / (2 rho0 c0^2) - f2 * (3/4) rho0 * <v^2> ]

    <p^2> = |p|^2 / 2,   <v^2> = |v|^2 / 2

    f1 = 1 - (rho0 c0^2) / (rho_p c_p^2)      (monopole term)
    f2 = 2 (rho_p - rho0) / (2 rho_p + rho0)  (dipole term)

    Vp = (4/3) pi a_p^3   (particle volume)

    F(r) = - grad U(r)

Particles are trapped at *minima* of U.  In air, for solid particles
(polystyrene, glass, water droplets) both f1 and f2 are ~1, so U is
minimised where pressure amplitude is low (a pressure node) and velocity
is low.

Run the demo with:  python acoustic_levitator.py

Dependencies: numpy, scipy, matplotlib.
"""

from __future__ import annotations

import numpy as np
from scipy.special import j1

__all__ = [
    "AIR", "POLYSTYRENE", "GLASS", "WATER",
    "Transducer", "LevitatorArray", "AcousticField", "PhaseControl",
]


# ---------------------------------------------------------------------------
# Physical constants and materials
# ---------------------------------------------------------------------------

AIR = dict(
    c=343.0,      # speed of sound [m/s] at 20 degC
    rho=1.204,    # density [kg/m^3]
)

# Common levitated-particle materials (c = longitudinal sound speed).
POLYSTYRENE = dict(rho=1050.0, c=2350.0)
GLASS       = dict(rho=2500.0, c=5640.0)
WATER       = dict(rho=1000.0, c=1480.0)

# Operating frequency of TM-2425H13T/R and most 40 kHz levitators.
FREQ_HZ = 40_000.0


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

class Transducer:
    """A single ultrasonic transducer (piston source)."""

    def __init__(self, position, normal, amplitude=1.0, phase=0.0,
                 radius=5.0e-3, group=None):
        # position / normal: numpy arrays [x, y, z] in metres.
        self.position = np.asarray(position, dtype=float)
        self.normal = np.asarray(normal, dtype=float)
        self.normal /= np.linalg.norm(self.normal) + 1e-30
        self.amplitude = float(amplitude)   # relative drive amplitude
        self.phase = float(phase)           # drive phase [rad]
        self.radius = float(radius)         # piston radius [m]
        self.group = group                  # grouping label (for phase control)


def _ring(center, axis, n, r, z, ring_phase_offset=0.0):
    """Return `n` transducers evenly spaced on a circle.

    center : 3-vector, centre of the ring.
    axis   : 3-vector, axis the ring is normal to (outward unit vector).
    r      : ring radius [m].
    z      : axial offset of the ring plane from `center` [m].
    """
    axis = np.asarray(axis, float)
    axis /= np.linalg.norm(axis) + 1e-30
    # Build two orthonormal tangent vectors.
    helper = np.array([1.0, 0.0, 0.0]) if abs(axis[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    t1 = np.cross(axis, helper); t1 /= np.linalg.norm(t1) + 1e-30
    t2 = np.cross(axis, t1);      t2 /= np.linalg.norm(t2) + 1e-30

    transducers = []
    for i in range(n):
        ang = 2 * np.pi * i / n + ring_phase_offset
        pos = np.asarray(center, float) + axis * z + r * (np.cos(ang) * t1 + np.sin(ang) * t2)
        # normal points back toward the array's focal point (inward)
        transducers.append(Transducer(pos, axis))
    return transducers


class LevitatorArray:
    """The 72-transducer levitator: 36 top + 36 bottom, facing each other.

    The two arrays are placed on concave spherical caps of curvature radius
    `Rc`, whose focal point sits at the origin (the levitation region).
    Transducers point inward toward the origin.

    Layout of each 36-element shell uses 3 concentric rings (6 + 12 + 18),
    matching the classic TinyLev transducer count.
    """

    # ring polar angles [rad] and element counts for a 36-element shell
    _RING_ANGLES = np.radians([16.0, 32.0, 48.0])
    _RING_COUNTS = [6, 12, 18]

    def __init__(self, curvature_radius=0.090, piston_radius=5.0e-3,
                 gap=0.0, freq=FREQ_HZ):
        """
        curvature_radius : spherical cap radius [m] (concavity of each dish).
        piston_radius    : transducer piston radius [m] (default 5 mm).
        gap              : extra axial separation [m] added between the two
                           caps (0 = caps share the same focal sphere).
        """
        self.Rc = float(curvature_radius)
        self.piston_radius = float(piston_radius)
        self.gap = float(gap)
        self.freq = float(freq)
        self.k = 2 * np.pi * self.freq / AIR["c"]
        self.lam = AIR["c"] / self.freq
        self.top = self._build_shell(+1)
        self.bottom = self._build_shell(-1)
        self.transducers = self.top + self.bottom
        self.grouping = "whole"      # 'whole' | 'petals' | 'rings'

    # -- geometry ----------------------------------------------------------
    def _build_shell(self, sign):
        """Build one 36-element concave shell. sign=+1 top, -1 bottom."""
        axis = np.array([0.0, 0.0, sign])          # outward axis
        center = np.array([0.0, 0.0, sign * self.gap / 2.0])
        shell = []
        for ang, n in zip(self._RING_ANGLES, self._RING_COUNTS):
            # position on the sphere: radius Rc, polar angle ang.
            # `_ring` already places the ring plane at +z *along* `axis`,
            # so the scalar offset is always +Rc*cos(ang); the sign of the
            # shell is carried by `axis`. (Previously the offset also carried
            # a -sign, which put BOTH shells at z<0.)
            r = self.Rc * np.sin(ang)              # ring radius in-plane
            z = self.Rc * np.cos(ang)              # axial offset (positive)
            shell.extend(_ring(center, axis, n, r, z))
        for tr in shell:
            tr.radius = self.piston_radius
            tr.normal = -(tr.position - np.array([0, 0, 0]))
            tr.normal /= np.linalg.norm(tr.normal) + 1e-30
        return shell

    # -- grouping / phase control -----------------------------------------
    def assign_groups(self, mode):
        """Tag each transducer with a group label used by phase control.

        mode == 'whole'  : all 72 transducers form one group.
        mode == 'petals' : each shell split into 6 angular sectors (6 x 6).
        mode == 'rings'  : each shell split into its 3 concentric rings.
        """
        self.grouping = mode
        for tr in self.transducers:
            tr.group = 0
        if mode == "whole":
            return

        for shell in (self.top, self.bottom):
            if mode == "rings":
                ring_radii = self.Rc * np.sin(self._RING_ANGLES)
                for tr in shell:
                    x, y, _ = tr.position
                    r = np.hypot(x, y)
                    # assign ring index by radial distance (nearest ring)
                    tr.group = int(np.argmin(np.abs(r - ring_radii)))
            elif mode == "petals":
                # Split into 6 equal angular sectors of 60 deg (6 x 6 per
                # shell).  Elements sit at every 10 deg of azimuth, so sort
                # by azimuth and chunk into 6 consecutive groups of 6.  This
                # is robust against floating-point boundary flips (the old
                # angle-binning formula produced 7 uneven groups).
                n_sectors = 6
                ordered = sorted(
                    shell,
                    key=lambda tr: np.arctan2(tr.position[1], tr.position[0]) % (2 * np.pi),
                )
                per = len(ordered) // n_sectors
                for i, tr in enumerate(ordered):
                    tr.group = i // per
        return

    def group_ids(self, shell="both"):
        """Sorted unique group labels for top, bottom, or both shells."""
        if shell == "top":
            ts = self.top
        elif shell == "bottom":
            ts = self.bottom
        else:
            ts = self.transducers
        return sorted({tr.group for tr in ts})

    def set_phase(self, phases, shell="both"):
        """Assign drive phases [rad] to transducers.

        phases : scalar (applied to all), a dict {group_id: phase}, or a
                 per-transducer sequence in the same order as the selected
                 `shell` (e.g. the arrays returned by the PhaseControl
                 recipes).
        shell  : 'top' | 'bottom' | 'both'.
        """
        ts = {"top": self.top, "bottom": self.bottom, "both": self.transducers}[shell]
        if isinstance(phases, dict):
            for tr in ts:
                tr.phase = float(phases[tr.group])
            return
        arr = np.asarray(phases)
        if arr.ndim == 0:                       # scalar -> apply to all
            for tr in ts:
                tr.phase = float(arr)
            return
        arr = arr.ravel()
        if arr.size != len(ts):
            raise ValueError(
                f"per-transducer phase array has {arr.size} entries "
                f"but shell '{shell}' has {len(ts)} transducers")
        for tr, ph in zip(ts, arr):
            tr.phase = float(ph)
        return

    def set_amplitude(self, amplitudes, shell="both"):
        ts = {"top": self.top, "bottom": self.bottom, "both": self.transducers}[shell]
        if np.isscalar(amplitudes):
            for tr in ts:
                tr.amplitude = float(amplitudes)
        else:
            for tr in ts:
                tr.amplitude = float(amplitudes[tr.group])
        return


# ---------------------------------------------------------------------------
# Field and potential computation
# ---------------------------------------------------------------------------

def _coord_gradients(f, X, Y, Z):
    """Return {coord_name: partial derivative} for each varying axis.

    Works for both 2-D axial slices and full 3-D grids.  For each axis of
    the scalar field `f`, identify which of X/Y/Z varies along it and take
    the gradient with respect to that coordinate (so the X-Z slice, where
    axis 1 is Z rather than Y, is handled automatically).
    """
    nd = f.ndim
    result = {}
    for a in range(nd):
        sl = [0] * nd
        sl[a] = slice(None)
        sl = tuple(sl)
        name = None
        for coord, nm in ((X, "x"), (Y, "y"), (Z, "z")):
            line = coord[sl]
            if line.ndim == 1 and line.size > 1 and np.ptp(line) > 1e-15:
                name = nm
                break
        if name is None:
            continue
        result[name] = np.gradient(f, line, axis=a, edge_order=1)
    return result


class AcousticField:
    """Computes pressure, velocity, Gor'kov potential and radiation force."""

    def __init__(self, array):
        self.array = array

    # -- pressure ----------------------------------------------------------
    def pressure(self, X, Y, Z, transducers=None):
        """Complex acoustic pressure on a grid (X, Y, Z in metres).

        Returns a complex numpy array of the same shape as X/Y/Z.
        """
        X, Y, Z = np.asarray(X), np.asarray(Y), np.asarray(Z)
        pts = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)  # (N,3)
        p = np.zeros(pts.shape[0], dtype=complex)
        k = self.array.k
        for tr in (transducers if transducers is not None else self.array.transducers):
            d = pts - tr.position                       # (N,3)
            R = np.linalg.norm(d, axis=1) + 1e-12       # distance
            cosang = (d @ tr.normal) / R                # cos(theta)
            sinang = np.sqrt(np.clip(1 - cosang ** 2, 0, 1))
            ka = k * tr.radius
            arg = ka * sinang
            D = np.where(arg < 1e-6, 1.0, 2 * j1(arg) / arg)   # piston directivity
            p += tr.amplitude * D * np.exp(1j * (k * R + tr.phase)) / R
        return p.reshape(X.shape)

    def velocity_sq(self, X, Y, Z, transducers=None):
        """|v|^2 of the particle velocity via finite differences of p.

        Uses the linearised momentum equation v = -grad(p)/(i w rho0),
        so |v|^2 = |grad(p)|^2 / (w^2 rho0^2).
        """
        p = self.pressure(X, Y, Z, transducers)
        grads = _coord_gradients(p, X, Y, Z)
        omega = 2 * np.pi * self.array.freq
        rho0 = AIR["rho"]
        return sum(np.abs(g) ** 2 for g in grads.values()) / (omega ** 2 * rho0 ** 2)

    # -- Gor'kov potential --------------------------------------------------
    def gorkov(self, X, Y, Z, particle=POLYSTYRENE, particle_radius=1.0e-3,
               transducers=None):
        """Gor'kov potential U [J] and radiation force F = -grad U [N]."""
        p = self.pressure(X, Y, Z, transducers)
        p2 = np.abs(p) ** 2 / 2.0                       # <p^2>
        v2 = self.velocity_sq(X, Y, Z, transducers) / 2.0   # <v^2>

        rho0, c0 = AIR["rho"], AIR["c"]
        rho_p, c_p = particle["rho"], particle["c"]
        f1 = 1.0 - (rho0 * c0 ** 2) / (rho_p * c_p ** 2)
        f2 = 2.0 * (rho_p - rho0) / (2.0 * rho_p + rho0)
        Vp = (4.0 / 3.0) * np.pi * particle_radius ** 3

        U = Vp * (f1 * p2 / (rho0 * c0 ** 2) - f2 * (3.0 / 4.0) * rho0 * v2)

        grads = _coord_gradients(U, X, Y, Z)
        zero = np.zeros_like(U)
        F = -np.stack([grads.get("x", zero), grads.get("y", zero),
                       grads.get("z", zero)], axis=0)   # (3, ...)
        return U, F

    def gorkov_coeffs(self, particle=POLYSTYRENE):
        """Return (f1, f2) for a given particle material."""
        rho0, c0 = AIR["rho"], AIR["c"]
        rho_p, c_p = particle["rho"], particle["c"]
        f1 = 1.0 - (rho0 * c0 ** 2) / (rho_p * c_p ** 2)
        f2 = 2.0 * (rho_p - rho0) / (2.0 * rho_p + rho0)
        return f1, f2


# ---------------------------------------------------------------------------
# Phase-control recipes
# ---------------------------------------------------------------------------

class PhaseControl:
    """Recipes that map target behaviour onto per-transducer phases.

    Each method returns a list of phases (one per transducer) that can be
    handed to LevitatorArray.set_phase().
    """

    def __init__(self, array):
        self.array = array

    def standing_wave(self, delta=0.0):
        """Anti-phase drive of top vs bottom (classic TinyLev).

        delta : differential phase offset [rad], applied with OPPOSITE sign
                to the two shells. This shifts the standing-wave nodes along
                the axis by +delta/k_z, where k_z = k*cos(alpha) is the
                effective axial wavenumber set by the mean polar angle alpha
                of the array (the mechanism for vertical movement). A common
                phase added to both shells would leave the nodes fixed,
                because it only multiplies the total field by a global phase
                factor e^{i*delta}.
        """
        phases = []
        for tr in self.array.transducers:
            sign = +1 if tr.position[2] > 0 else -1
            phases.append(sign * (np.pi / 2 + delta))
        return np.array(phases)

    def focus(self, target, phase_only=True):
        """Phase-conjugate focusing onto `target` [x,y,z] metres.

        phi_j = -k |target - r_j| makes all waves arrive in phase at the
        target, creating a pressure maximum there.  Dense particles are
        *repelled* by a pressure maximum, so single-focus traps are usually
        combined with a pressure minimum (see twin-trap / vortex).
        """
        target = np.asarray(target, float)
        phases = []
        for tr in self.array.transducers:
            R = np.linalg.norm(target - tr.position)
            phases.append(-self.array.k * R if phase_only else self.array.k * R)
        return np.array(phases)

    def steer(self, direction, magnitude=1.0):
        """Linear phase ramp -> beam / trap steering.

        direction : unit-ish 3-vector giving the steer direction.
        magnitude : scales the ramp (larger -> larger displacement).
        phi_j = -magnitude * k * (direction . r_j)
        """
        d = np.asarray(direction, float)
        d = d / (np.linalg.norm(d) + 1e-30)
        phases = []
        for tr in self.array.transducers:
            phases.append(-magnitude * self.array.k * np.dot(d, tr.position))
        return np.array(phases)

    def vortex(self, m=1, shell="both"):
        """Azimuthal vortex phase -> orbital angular momentum / rotation.

        m       : topological charge (+1/-1 sets spin direction, |m|>=2 for
                  higher-order vortices).
        shell   : apply to 'top', 'bottom', or 'both'.
        phi_j = m * atan2(y_j, x_j)
        """
        phases = np.zeros(len(self.array.transducers))
        for i, tr in enumerate(self.array.transducers):
            ts = self.array.transducers
            if shell == "top" and tr not in self.array.top:
                continue
            if shell == "bottom" and tr not in self.array.bottom:
                continue
            phases[i] = m * np.arctan2(tr.position[1], tr.position[0])
        return phases

    def petal_tilt(self, direction, magnitude=0.5):
        """Per-petal phase ramp (flower-petal grouping) for lateral tilt.

        Applies a phase proportional to each transducer's projection onto
        `direction`, which tilts the standing-wave axis.
        """
        d = np.asarray(direction, float)
        d = d / (np.linalg.norm(d) + 1e-30)
        phases = []
        for tr in self.array.transducers:
            phases.append(magnitude * np.dot(d, tr.position[:2]))
        return np.array(phases)


# ---------------------------------------------------------------------------
# Demo / visualisation
# ---------------------------------------------------------------------------

def _mesh_axis(zmin, zmax, nz, xmin, xmax, nx, y=0.0):
    """Return (X, Y, Z) for an X-Z axial slice at a fixed y (metres)."""
    xs = np.linspace(xmin, xmax, nx)
    zs = np.linspace(zmin, zmax, nz)
    X, Z = np.meshgrid(xs, zs, indexing="ij")
    Y = np.full_like(X, y)
    return X, Y, Z


def demo():
    """Run the simulation and save example figures to ../outputs."""
    import os
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "outputs")
    outdir = os.path.abspath(outdir)
    os.makedirs(outdir, exist_ok=True)

    mm = 1e-3
    lev = LevitatorArray(curvature_radius=0.090, piston_radius=5 * mm)
    field = AcousticField(lev)
    ctrl = PhaseControl(lev)

    print(f"Frequency     : {lev.freq/1e3:.1f} kHz")
    print(f"Wavelength    : {lev.lam*1e3:.2f} mm")
    print(f"Transducers   : {len(lev.transducers)} (top {len(lev.top)}, bottom {len(lev.bottom)})")
    f1, f2 = field.gorkov_coeffs(POLYSTYRENE)
    print(f"Gor'kov f1,f2 : {f1:.4f}, {f2:.4f}  (polystyrene in air)")
    print(f"k             : {lev.k:.2f} rad/m")

    # ---- 1. Standing wave (classic TinyLev) -----------------------------
    lev.assign_groups("whole")
    lev.set_phase(ctrl.standing_wave(0.0))
    X, Y, Z = _mesh_axis(-25 * mm, 25 * mm, 400, -25 * mm, 25 * mm, 400, y=0.0)
    p = field.pressure(X, Y, Z)
    U, F = field.gorkov(X, Y, Z, particle=POLYSTYRENE, particle_radius=1 * mm)

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))
    axes[0].pcolormesh(X * 1e3, Z * 1e3, np.abs(p) / np.abs(p).max(), cmap="inferno", shading="auto")
    axes[0].set_title("|p| (standing wave)")
    axes[0].set_xlabel("x [mm]"); axes[0].set_ylabel("z [mm]")
    axes[1].pcolormesh(X * 1e3, Z * 1e3, U / np.abs(U).max(), cmap="RdBu_r", shading="auto")
    axes[1].set_title("Gor'kov potential U (blue = trap)")
    axes[1].set_xlabel("x [mm]"); axes[1].set_ylabel("z [mm]")
    step = 12
    axes[2].quiver(X[::step, ::step] * 1e3, Z[::step, ::step] * 1e3,
                   F[0, ::step, ::step], F[2, ::step, ::step], color="k")
    axes[2].set_title("Radiation force F = -grad U")
    axes[2].set_xlabel("x [mm]"); axes[2].set_ylabel("z [mm]")
    fig.suptitle("72-transducer levitator: standing-wave trap (polystyrene, 1 mm bead)")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "01_standing_wave.png"), dpi=120)
    plt.close(fig)

    # ---- 2. Phase shift moves the trap along z ---------------------------
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))
    for ax, delta in zip(axes, [0.0, np.pi / 4, np.pi / 2]):
        lev.set_phase(ctrl.standing_wave(delta))
        p = field.pressure(X, Y, Z)
        U, _ = field.gorkov(X, Y, Z, particle=POLYSTYRENE, particle_radius=1 * mm)
        ax.pcolormesh(X * 1e3, Z * 1e3, U / np.abs(U).max(), cmap="RdBu_r", shading="auto")
        ax.set_title(f"delta = {delta:.2f} rad")
        ax.set_xlabel("x [mm]"); ax.set_ylabel("z [mm]")
    fig.suptitle("Axial translation: relative phase between top and bottom arrays")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "02_translation_phase.png"), dpi=120)
    plt.close(fig)

    # ---- 3. Vortex (rotation / orbital angular momentum) -----------------
    lev.assign_groups("whole")
    lev.set_phase(ctrl.vortex(m=1, shell="both"))
    # transverse slice at the levitation plane
    xs = np.linspace(-20 * mm, 20 * mm, 300)
    ys = np.linspace(-20 * mm, 20 * mm, 300)
    Xx, Yy = np.meshgrid(xs, ys, indexing="ij")
    Zz = np.zeros_like(Xx)
    p_xy = field.pressure(Xx, Yy, Zz)
    U_xy, F_xy = field.gorkov(Xx, Yy, Zz, particle=POLYSTYRENE, particle_radius=1 * mm)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    axes[0].pcolormesh(Xx * 1e3, Yy * 1e3, np.abs(p_xy) / np.abs(p_xy).max(), cmap="inferno", shading="auto")
    axes[0].set_title("|p| in levitation plane (vortex m=1)")
    axes[0].set_xlabel("x [mm]"); axes[0].set_ylabel("y [mm]")
    axes[1].pcolormesh(Xx * 1e3, Yy * 1e3, U_xy / np.abs(U_xy).max(), cmap="RdBu_r", shading="auto")
    step = 8
    axes[1].quiver(Xx[::step, ::step] * 1e3, Yy[::step, ::step] * 1e3,
                   F_xy[0, ::step, ::step], F_xy[1, ::step, ::step], color="k")
    axes[1].set_title("Gor'kov U + force (vortex trap)")
    axes[1].set_xlabel("x [mm]"); axes[1].set_ylabel("y [mm]")
    fig.suptitle("Vortex phase -> ring trap with orbital angular momentum")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "03_vortex_rotation.png"), dpi=120)
    plt.close(fig)

    # ---- 4. Flower-petal grouping (6 sectors x 6) ------------------------
    lev.assign_groups("petals")
    # colour each petal with a distinct phase to show the pattern
    lev.set_phase(ctrl.standing_wave(0.0))             # both shells: anti-phase
    petal_phases = {g: (np.pi / 3) * g for g in lev.group_ids("top")}
    lev.set_phase(petal_phases, shell="top")           # override top shell

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    # top-view of the top shell, coloured by phase
    top_pos = np.array([t.position for t in lev.top])
    top_ph = np.array([t.phase for t in lev.top])
    sc = axes[0].scatter(top_pos[:, 0] * 1e3, top_pos[:, 1] * 1e3, c=top_ph,
                         cmap="twilight", s=80)
    axes[0].set_aspect("equal")
    axes[0].set_title("Top shell: flower-petal phases (6 sectors)")
    axes[0].set_xlabel("x [mm]"); axes[0].set_ylabel("y [mm]")
    plt.colorbar(sc, ax=axes[0], label="phase [rad]")

    p = field.pressure(X, Y, Z)
    U, _ = field.gorkov(X, Y, Z, particle=POLYSTYRENE, particle_radius=1 * mm)
    axes[1].pcolormesh(X * 1e3, Z * 1e3, U / np.abs(U).max(), cmap="RdBu_r", shading="auto")
    axes[1].set_title("Gor'kov potential (petal phase pattern)")
    axes[1].set_xlabel("x [mm]"); axes[1].set_ylabel("z [mm]")
    fig.suptitle("Flower-petal grouping: 6 sectors x 6 transducers per shell")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "04_flower_petals.png"), dpi=120)
    plt.close(fig)

    # ---- 5. Concentric-ring grouping (3 rings per shell) ------------------
    lev.assign_groups("rings")
    lev.set_phase(ctrl.standing_wave(0.0))             # both shells: anti-phase
    ring_phases = {g: (np.pi / 2) * g for g in lev.group_ids("top")}
    lev.set_phase(ring_phases, shell="top")            # override top shell

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    top_pos = np.array([t.position for t in lev.top])
    top_ph = np.array([t.phase for t in lev.top])
    sc = axes[0].scatter(top_pos[:, 0] * 1e3, top_pos[:, 1] * 1e3, c=top_ph,
                         cmap="twilight", s=80)
    axes[0].set_aspect("equal")
    axes[0].set_title("Top shell: concentric-ring phases (3 rings)")
    axes[0].set_xlabel("x [mm]"); axes[0].set_ylabel("y [mm]")
    plt.colorbar(sc, ax=axes[0], label="phase [rad]")

    p = field.pressure(X, Y, Z)
    U, _ = field.gorkov(X, Y, Z, particle=POLYSTYRENE, particle_radius=1 * mm)
    axes[1].pcolormesh(X * 1e3, Z * 1e3, U / np.abs(U).max(), cmap="RdBu_r", shading="auto")
    axes[1].set_title("Gor'kov potential (ring phase pattern)")
    axes[1].set_xlabel("x [mm]"); axes[1].set_ylabel("z [mm]")
    fig.suptitle("Concentric-ring grouping: 3 rings per shell (6/12/18)")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "05_concentric_rings.png"), dpi=120)
    plt.close(fig)

    print(f"\nSaved figures to {outdir}/")
    return lev, field, ctrl


if __name__ == "__main__":
    demo()